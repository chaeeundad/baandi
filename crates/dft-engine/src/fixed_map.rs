//! Diagnostic accessors for same-input SCF-map comparisons. No SCF dispatcher
//! calls these methods. Full density/orbitals are materialized for inspection;
//! this is not a memory-scaling or timing interface.
use super::*;
use dft_basis::HalfSphere;
use dft_mpi::Collective;
use dft_numerics::{RealMatrix, pool::WorkerPool};
use dft_physics::distributed_real_gamma::DistributedRealHamiltonian;
use dft_solvers::{SolverError, distributed_davidson::DistributedRealOperator, ref_davidson};
use std::sync::Arc;

pub struct ReferenceMap {
    pub density: Array3<f64>,
    pub orbitals: ComplexMatrix,
    pub eigenvalues: Vec<f64>,
    pub residuals: Vec<f64>,
    pub occupations: Vec<f64>,
    pub energy: energy::EnergyAudit,
    pub converged: Vec<bool>,
    pub record: ref_davidson::Record,
}

fn numerical(e: impl std::fmt::Display) -> ScfError {
    ScfError::Numerical(e.to_string())
}

struct Operator {
    h: DistributedRealHamiltonian,
    pool: Arc<WorkerPool>,
}
impl DistributedRealOperator for Operator {
    fn global_len(&self) -> usize {
        self.h.global_len()
    }
    fn local_rows(&self) -> &[usize] {
        self.h.local_rows()
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply_block<C: Collective>(
        &mut self,
        comm: &C,
        x: &RealMatrix,
        y: &mut RealMatrix,
    ) -> Result<(), SolverError> {
        self.h
            .apply_block(comm, x, y)
            .map_err(|e| SolverError::Application(e.to_string()))
    }
}

impl ScfProblem {
    /// Fixed integer band occupations used by the diagnostic map.
    pub fn fixed_map_occupations(&self) -> &[f64] {
        &self.occupations
    }

    /// Production input fields in Ha, without solving or mixing. Collective
    /// under MPI (XC may use the registered communicator). No clipping or
    /// renormalization of the supplied density.
    pub fn fixed_map_fields(
        &mut self,
        n: &Array3<f64>,
    ) -> Result<(hartree::HartreeGrid, xc::XcGrid, Array3<f64>), ScfError> {
        if n.shape() != self.grid.dims() || n.iter().any(|v| !v.is_finite()) {
            return Err(ScfError::InvalidInput("fixed-map density shape/nonfinite"));
        }
        self.input_fields(n)
    }

    /// One opt-in reference-Davidson solve using the existing H, density and
    /// energy kernels. Threshold is a Ritz change in Ha; no residual refinement,
    /// SCF schedule, mixing or verification. Callers must inspect every flag.
    /// A cold call uses BAANDI's atomic orbitals; subsequent calls may supply
    /// the preceding returned block. This is NOT a QE wavefunction replay.
    pub fn reference_fixed_map(
        &mut self,
        n: &Array3<f64>,
        initial: Option<&ComplexMatrix>,
        ethr_ha: f64,
        first_step: bool,
    ) -> Result<ReferenceMap, ScfError> {
        if let Some(comm) = crate::parallel::SharedComm::active() {
            self.reference_fixed_map_inner(&comm, n, initial, ethr_ha, first_step)
        } else {
            self.reference_fixed_map_inner(&refscf::metric::Serial, n, initial, ethr_ha, first_step)
        }
    }

    fn reference_fixed_map_inner<C: Collective>(
        &mut self,
        comm: &C,
        n: &Array3<f64>,
        initial: Option<&ComplexMatrix>,
        ethr_ha: f64,
        first_step: bool,
    ) -> Result<ReferenceMap, ScfError> {
        if self.states.len() != 1
            || self.states[0].point.fractional != [0.; 3]
            || self.symmetry.is_some()
            || self.options.gamma_real == crate::GammaReal::Off
        {
            return Err(ScfError::InvalidInput(
                "reference fixed map requires real Gamma without symmetry",
            ));
        }
        let (h, x, effective) = self.fixed_map_fields(n)?;
        let atomic;
        let initial = match initial {
            Some(v) => v,
            None => {
                atomic = self.atomic_start_vectors(n)?.ok_or(ScfError::InvalidInput(
                    "reference fixed map requires atomic starting orbitals or an explicit block",
                ))?;
                &atomic[0]
            }
        };
        let state = &mut self.states[0];
        state.hamiltonian.set_local_potential(effective)?;
        self.owned.complex_stale = false;
        let bands = self.occupations.len();
        let half = HalfSphere::new(&state.basis, &state.grid).map_err(numerical)?;
        let dim = half.full_len();
        if initial.rows() != dim || initial.cols() != bands {
            return Err(ScfError::InvalidInput("reference fixed-map start shape"));
        }
        let mut full = RealMatrix::zeros(dim, bands).map_err(numerical)?;
        for b in 0..bands {
            half.pack(
                initial.column(b),
                &mut full.as_mut_slice()[b * dim..(b + 1) * dim],
            )
            .map_err(numerical)?;
        }
        let mut op = Operator {
            h: DistributedRealHamiltonian::new(
                &state.hamiltonian,
                &state.basis,
                bands,
                comm.rank(),
                comm.size(),
                self.pool.clone(),
            )
            .map_err(numerical)?,
            pool: self.pool.clone(),
        };
        let start = op.h.extract_local(&full).map_err(numerical)?;
        // Full plane-wave diagonal, in Ha. Each Gamma pair shares a value;
        // ORC packing has already accounted for the pair normalization.
        let source = &state.hamiltonian;
        let v0 = source.local_potential_values().iter().sum::<f64>()
            / source.local_potential_values().len() as f64;
        let mut diagonal = source.kinetic_diagonal().to_vec();
        for (g, value) in diagonal.iter_mut().enumerate() {
            *value += v0;
            for atom in source.nonlocal().atoms() {
                let beta = atom.columns();
                let count = atom.channel_count();
                for j in 0..count {
                    for i in 0..count {
                        *value += atom.coupling_hartree()[i + count * j]
                            * (beta.column(i)[g] * beta.column(j)[g].conj()).re;
                    }
                }
            }
        }
        let diagonal: Vec<_> = op
            .local_rows()
            .iter()
            .map(|&r| diagonal[half.plus()[r.div_ceil(2)] as usize])
            .collect();
        let mut opts = ref_davidson::Options::new(
            bands,
            ethr_ha,
            self.occupations.iter().filter(|&&f| f > 0.).count(),
        );
        opts.first_step = first_step;
        opts.max_iterations = self.options.davidson_max_iterations;
        opts.return_hx = true;
        let solved = ref_davidson::solve(comm, &mut op, &start, &diagonal, opts)?;
        let rows = op.local_rows().len();
        let hx = solved.hx.as_ref().expect("requested H X");
        let local_norms: Vec<_> = (0..bands)
            .map(|b| {
                (0..rows)
                    .map(|r| {
                        (hx.as_slice()[b * rows + r]
                            - solved.eigenvalues[b] * solved.eigenvectors.as_slice()[b * rows + r])
                            .powi(2)
                    })
                    .sum()
            })
            .collect();
        let mut residuals = vec![0.; bands];
        comm.tree_sum_f64(&local_norms, &mut residuals)
            .map_err(numerical)?;
        residuals.iter_mut().for_each(|r| *r = r.sqrt());
        full.as_mut_slice().fill(0.);
        op.h.scatter_local(&solved.eigenvectors, &mut full)
            .map_err(numerical)?;
        let mut gathered = RealMatrix::zeros(dim, bands).map_err(numerical)?;
        comm.tree_sum_f64(full.as_slice(), gathered.as_mut_slice())
            .map_err(numerical)?;
        let mut orbitals = ComplexMatrix::zeros(dim, bands).map_err(numerical)?;
        for b in 0..bands {
            half.unpack(
                &gathered.as_slice()[b * dim..(b + 1) * dim],
                orbitals.column_mut(b),
            )
            .map_err(numerical)?;
        }
        let coefficients: Vec<_> = (0..bands).map(|b| orbitals.column(b).to_vec()).collect();
        state.band_density.compute_keep_for_map(
            &coefficients,
            &self.occupations,
            &state.grid,
            self.cell.volume_bohr3(),
            true,
            None,
            false,
        )?;
        let density = state.band_density.density().clone();
        let output_xc = self.options.evaluate_xc(
            &density,
            self.core_density.as_ref(),
            &mut self.xc_density,
            &self.cell,
            &mut self.fft,
            self.density_floor,
        )?;
        let energy = energy::evaluate_with_core_xc(
            AuditInput {
                cell: &self.cell,
                basis: &state.basis,
                grid: &state.grid,
                ionic_g: &self.ionic_g,
                nonlocal: state.hamiltonian.nonlocal(),
                orbitals: &coefficients,
                occupations: &self.occupations,
                eigenvalues: &solved.eigenvalues,
                eigen_residuals: &residuals,
                input_density: n,
                output_density: &density,
                input_hartree: &h,
                input_xc: &x,
                ewald: self.ewald.energy.0,
            },
            &mut self.fft,
            self.options.functional,
            Some(&output_xc),
            self.core_density.as_ref(),
        )?;
        Ok(ReferenceMap {
            density,
            orbitals,
            eigenvalues: solved.eigenvalues,
            residuals,
            occupations: self.occupations.clone(),
            energy,
            converged: solved.converged,
            record: solved.record,
        })
    }
}
