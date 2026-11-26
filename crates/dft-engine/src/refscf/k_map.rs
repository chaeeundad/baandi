//! Sequential-k maps sharing the real driver's density functional and mixing.
use super::{
    map::{self, Fields},
    numerical,
};
use crate::{
    Functional,
    energy::EnergyAudit,
    model::PhaseTimings,
    model::ScfError,
    refstart::{KOrbitals, KStartState},
};
use dft_geometry::Cell;
use dft_mpi::Collective;
use dft_solvers::{ref_davidson, ref_davidson_complex};
use std::time::Instant;
pub(super) struct PerK {
    pub eigenvalues: Vec<f64>,
    pub residuals: Vec<f64>,
    pub converged: Vec<bool>,
    pub calls: Vec<super::record::SolveRecord>,
}
pub(super) struct Solution {
    pub eigenvectors: KOrbitals,
    pub hx: Option<KOrbitals>,
    pub per_k: Vec<PerK>,
}
pub(super) struct Map {
    pub solve: Solution,
    pub warm_start: bool,
    pub density: Vec<f64>,
    pub residuals: Vec<f64>,
    pub energy: EnergyAudit,
    pub timings: PhaseTimings,
    input_fields: Fields,
}
fn sum<K: Collective>(comm: &K, v: &[f64]) -> Result<Vec<f64>, ScfError> {
    let mut out = vec![0.; v.len()];
    comm.tree_sum_f64(v, &mut out).map_err(numerical)?;
    Ok(out)
}
pub(super) fn cached_fields<K: Collective>(
    state: &mut KStartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    density: &[f64],
    cache: &mut Option<(Vec<f64>, Fields)>,
) -> Result<Fields, ScfError> {
    map::cached_fields(&mut state.common, comm, cell, functional, density, cache)
}
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate<K: Collective>(
    state: &mut KStartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    input: &[f64],
    occupations: &[f64],
    ethr_ry: f64,
    first_step: bool,
    diagnostics: bool,
    cache: &mut Option<(Vec<f64>, Fields)>,
    lean: bool,
    kept: &mut super::k_kept::KeptHx,
    symmetry: &mut Option<super::symmetry::Workspace>,
) -> Result<Map, ScfError> {
    let mut timings = PhaseTimings::default();
    let clock = Instant::now();
    let input_fields = cached_fields(state, comm, cell, functional, input, cache)?;
    let potential: Vec<_> = state
        .ionic
        .iter()
        .zip(&input_fields.hxc)
        .map(|(a, b)| a + b)
        .collect();
    let dims = state.grid.dims();
    let points = dims.iter().product::<usize>();
    let v0 = sum(comm, &[potential.iter().sum::<f64>()])?[0] / points as f64;
    timings.potential = clock.elapsed();
    let bands = occupations.len();
    let mut vectors = Vec::new();
    let mut products = Vec::new();
    let mut per_k = Vec::new();
    let mut warm_start = false;
    let mut density = vec![0.; input.len()];
    let mut kinetic = 0.;
    let mut nonlocal = 0.;
    let mut band_sum = 0.;
    let mut residual_bound = 0.;
    let mut residuals = vec![0_f64; bands];
    // GK4: ordinary (lean) maps consume only the band sum, the double counting
    // and the residual norms; kinetic, nonlocal, local and zero-mode terms are
    // filled by `complete_audit`. All per-k residual partials and the energy
    // integrals share one element-wise tree sum per map.
    let deferred = dft_numerics::switches::on(dft_numerics::switches::Id::K_LEAN_ENERGY);
    let band_terms = !(deferred && lean);
    let mut residual_local = Vec::with_capacity(if deferred {
        state.blocks.len() * bands
    } else {
        0
    });
    for k in 0..state.blocks.len() {
        let mut planes = crate::grid_dist::PlaneMap::new(
            dims,
            comm.rank(),
            comm.size(),
            state.fields.first_axis(),
            state.blocks[k].operator.first_axis(),
        );
        let mut op_potential = Vec::new();
        planes.a_to_b(comm, &potential, &mut op_potential)?;
        let op = &mut state.blocks[k].operator;
        let epoch = op
            .owned_epoch()
            .unwrap_or(0)
            .checked_add(1)
            .ok_or(ScfError::InvalidInput("k potential epoch overflow"))?;
        op.set_owned_potential(&op_potential, epoch)
            .map_err(numerical)?;
        let diagonal: Vec<_> = op.diagonal_static().iter().map(|v| v + v0).collect();
        let initial_hx = kept.prepare(k, comm, op, &state.orbitals.0[k], &op_potential)?;
        warm_start |= initial_hx.is_some();
        let mut options = ref_davidson::Options::new(
            bands,
            ethr_ry / 2.,
            occupations.iter().filter(|&&f| f > 0.).count(),
        );
        // GK10B research (`K_FIRST_EMPTY_RULE`, default off): this driver's
        // occupations are fixed and known before the first solve, so the
        // first map may keep the empty-band threshold for empty bands
        // instead of converging every band to the occupied threshold. It
        // changes the map-1 density when no redo follows.
        options.first_step = first_step
            && !dft_numerics::switches::on(dft_numerics::switches::Id::K_FIRST_EMPTY_RULE);
        options.return_hx = true;
        options.profile = diagnostics;
        let clock = Instant::now();
        let mut calls = Vec::new();
        let mut solve = ref_davidson_complex::solve_with_initial_hx(
            comm,
            op,
            &state.orbitals.0[k],
            &diagonal,
            options,
            initial_hx.as_ref(),
        )
        .map_err(numerical)?;
        for retry in 0..6 {
            calls.push(super::k_driver::solve_record(
                &mut solve,
                k,
                ethr_ry,
                first_step && retry == 0,
                diagnostics,
            ));
            if solve.converged.iter().filter(|&&v| !v).count() <= 5.max(bands / 4) {
                break;
            }
            if retry == 5 {
                return Err(ScfError::InvalidInput(
                    "complex reference Davidson exhausted six calls",
                ));
            }
            options.first_step = false;
            solve = ref_davidson_complex::solve(comm, op, &solve.eigenvectors, &diagonal, options)
                .map_err(numerical)?;
        }
        timings.diagonalization += clock.elapsed();
        let clock = Instant::now();
        let weight = state.blocks[k].point.weight;
        let weights: Vec<_> = occupations.iter().map(|f| f * weight).collect();
        let op = &mut state.blocks[k].operator;
        let hx = solve.hx.as_ref().expect("HX requested");
        let local: Vec<_> = (0..bands)
            .map(|j| {
                (0..op.local_rows().len())
                    .map(|i| {
                        (hx[(i, j)] - solve.eigenvalues[j] * solve.eigenvectors[(i, j)]).norm_sqr()
                    })
                    .sum()
            })
            .collect();
        let rs: Vec<_> = if deferred {
            residual_local.extend(local);
            Vec::new()
        } else {
            sum(comm, &local)?.into_iter().map(f64::sqrt).collect()
        };
        for j in 0..bands {
            if !deferred {
                residuals[j] = residuals[j].max(rs[j]);
            }
            band_sum += weights[j] * solve.eigenvalues[j];
            if !deferred {
                residual_bound += weights[j] * rs[j];
            }
            if band_terms {
                kinetic += band_kinetic(op, &solve.eigenvectors, j, weights[j]);
            }
        }
        if band_terms {
            nonlocal += band_nonlocal(op, comm, &solve.eigenvectors, &weights)?;
        }
        let mut op_density = vec![0.; op.plane_points()];
        op.band_density_sum(
            comm,
            &solve.eigenvectors,
            &weights,
            cell.volume_bohr3(),
            &mut op_density,
        )
        .map_err(numerical)?;
        let mut local_density = Vec::new();
        planes.b_to_a(comm, &op_density, &mut local_density)?;
        for (a, b) in density.iter_mut().zip(local_density) {
            *a += b;
        }
        kept.retain(k, state, &solve.eigenvectors, hx, &op_potential);
        timings.density += clock.elapsed();
        per_k.push(PerK {
            eigenvalues: solve.eigenvalues,
            residuals: rs,
            converged: solve.converged,
            calls,
        });
        vectors.push(solve.eigenvectors);
        products.push(solve.hx.unwrap());
    }
    if let Some(s) = symmetry {
        let clock = Instant::now();
        s.apply(&mut state.fields, comm, std::slice::from_mut(&mut density))?;
        timings.density += clock.elapsed();
    }
    let clock = Instant::now();
    let dv = cell.volume_bohr3() / points as f64;
    let integrals = if deferred {
        let mut packed = residual_local;
        let r = packed.len();
        if lean {
            packed.push(double_counting_local(dv, &density, &input_fields.hxc));
        } else {
            packed.extend(integral_locals(
                kinetic,
                dv,
                &density,
                input,
                &state.ionic,
                points,
                &input_fields,
            ));
        }
        let mut total = sum(comm, &packed)?;
        for (k, (info, block)) in per_k.iter_mut().zip(&state.blocks).enumerate() {
            info.residuals = total[k * bands..(k + 1) * bands]
                .iter()
                .map(|v| v.sqrt())
                .collect();
            let weight = block.point.weight;
            for j in 0..bands {
                residuals[j] = residuals[j].max(info.residuals[j]);
                if !lean {
                    residual_bound += occupations[j] * weight * info.residuals[j];
                }
            }
        }
        if lean {
            let ewald = state.ewald.energy.0;
            // NaN is an unavailable audit slot (JSON null), never a stop input;
            // `complete_audit` fills every slot before delivery.
            let energy = EnergyAudit {
                xc_floor: None,
                nlcc: None,
                kinetic: f64::NAN,
                local: f64::NAN,
                nonlocal: f64::NAN,
                local_zero_mode: f64::NAN,
                hartree: f64::NAN,
                xc: f64::NAN,
                ewald,
                band_sum,
                double_counting: total[r],
                band_form: f64::NAN,
                harris_foulkes: band_sum - input_fields.hartree + input_fields.xc
                    - input_fields.xc_integral
                    + ewald,
                functional_difference_bound: f64::NAN,
            };
            timings.energy = clock.elapsed();
            return Ok(Map {
                solve: Solution {
                    eigenvectors: KOrbitals(vectors),
                    hx: Some(KOrbitals(products)),
                    per_k,
                },
                warm_start,
                density,
                residuals,
                energy,
                timings,
                input_fields,
            });
        }
        total.drain(..r);
        total
    } else {
        sum(
            comm,
            &integral_locals(
                kinetic,
                dv,
                &density,
                input,
                &state.ionic,
                points,
                &input_fields,
            ),
        )?
    };
    let output_fields = if lean {
        None
    } else {
        Some(cached_fields(
            state, comm, cell, functional, &density, cache,
        )?)
    };
    let (hartree, xc, nlcc) =
        output_fields.map_or((f64::NAN, f64::NAN, None), |f| (f.hartree, f.xc, f.nlcc));
    let ewald = state.ewald.energy.0;
    let energy = EnergyAudit {
        xc_floor: None,
        nlcc,
        kinetic: integrals[0],
        local: integrals[1],
        nonlocal,
        local_zero_mode: state.electrons * integrals[3],
        hartree,
        xc,
        ewald,
        band_sum,
        double_counting: integrals[2],
        band_form: band_sum - integrals[2] + hartree + xc + ewald,
        harris_foulkes: band_sum - input_fields.hartree + input_fields.xc
            - input_fields.xc_integral
            + ewald,
        functional_difference_bound: (hartree - input_fields.hartree - integrals[5]).abs()
            + (xc - input_fields.xc - (integrals[4] - integrals[5])).abs()
            + residual_bound,
    };
    timings.energy = clock.elapsed();
    Ok(Map {
        solve: Solution {
            eigenvectors: KOrbitals(vectors),
            hx: Some(KOrbitals(products)),
            per_k,
        },
        warm_start,
        density,
        residuals,
        energy,
        timings,
        input_fields,
    })
}
#[allow(clippy::too_many_arguments)]
pub(super) fn complete_audit<K: Collective>(
    step: &mut Map,
    state: &mut KStartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    input: &[f64],
    occupations: &[f64],
    cache: &mut Option<(Vec<f64>, Fields)>,
) -> Result<(), ScfError> {
    if step.energy.kinetic.is_nan() {
        // GK4 deferred lean map: the per-k band terms in the map's k and band
        // order, then the integrals in the map's slot order (element-wise
        // tree sums, so each slot equals the eager map's value bit for bit).
        let mut kinetic = 0.;
        let mut nonlocal = 0.;
        for (block, x) in state.blocks.iter_mut().zip(&step.solve.eigenvectors.0) {
            let weights: Vec<_> = occupations.iter().map(|f| f * block.point.weight).collect();
            for (j, &w) in weights.iter().enumerate() {
                kinetic += band_kinetic(&block.operator, x, j, w);
            }
            nonlocal += band_nonlocal(&mut block.operator, comm, x, &weights)?;
        }
        let points = state.grid.dims().iter().product::<usize>();
        let dv = cell.volume_bohr3() / points as f64;
        let integrals = sum(
            comm,
            &integral_locals(
                kinetic,
                dv,
                &step.density,
                input,
                &state.ionic,
                points,
                &step.input_fields,
            ),
        )?;
        let fields = cached_fields(state, comm, cell, functional, &step.density, cache)?;
        let residual = residual_term(&step.solve.per_k, state, occupations);
        let e = &mut step.energy;
        e.kinetic = integrals[0];
        e.local = integrals[1];
        e.nonlocal = nonlocal;
        e.local_zero_mode = state.electrons * integrals[3];
        e.hartree = fields.hartree;
        e.xc = fields.xc;
        e.nlcc = fields.nlcc;
        e.band_form = e.band_sum - e.double_counting + e.hartree + e.xc + e.ewald;
        e.functional_difference_bound = (e.hartree - step.input_fields.hartree - integrals[5])
            .abs()
            + (e.xc - step.input_fields.xc - (integrals[4] - integrals[5])).abs()
            + residual;
    } else if step.energy.hartree.is_nan() {
        let fields = cached_fields(state, comm, cell, functional, &step.density, cache)?;
        let dv = cell.volume_bohr3() / state.grid.dims().iter().product::<usize>() as f64;
        let differences = sum(
            comm,
            &[
                dv * step
                    .density
                    .iter()
                    .zip(input)
                    .zip(&step.input_fields.hxc)
                    .map(|((a, b), v)| (a - b) * v)
                    .sum::<f64>(),
                dv * step
                    .density
                    .iter()
                    .zip(input)
                    .zip(&step.input_fields.vh)
                    .map(|((a, b), v)| (a - b) * v)
                    .sum::<f64>(),
            ],
        )?;
        let residual = residual_term(&step.solve.per_k, state, occupations);
        let e = &mut step.energy;
        e.hartree = fields.hartree;
        e.xc = fields.xc;
        e.nlcc = fields.nlcc;
        e.band_form = e.band_sum - e.double_counting + e.hartree + e.xc + e.ewald;
        e.functional_difference_bound = (e.hartree - step.input_fields.hartree - differences[1])
            .abs()
            + (e.xc - step.input_fields.xc - (differences[0] - differences[1])).abs()
            + residual;
    }
    if !step.energy.total_hartree().is_finite() {
        return Err(ScfError::NonFinite);
    }
    Ok(())
}

/// Residual term of the functional-difference bound as the completed audit
/// has always formed it (per-k weight times the occupied residual sum).
fn residual_term(per_k: &[PerK], state: &KStartState, occupations: &[f64]) -> f64 {
    per_k
        .iter()
        .zip(&state.blocks)
        .map(|(s, k)| {
            k.point.weight
                * s.residuals
                    .iter()
                    .zip(occupations)
                    .map(|(r, f)| r * f)
                    .sum::<f64>()
        })
        .sum::<f64>()
}
/// Weighted kinetic energy of band `j` on this rank's rows (summed later).
fn band_kinetic(
    op: &dft_physics::distributed_complex::DistributedComplexHamiltonian,
    x: &dft_physics::data::ComplexMatrix,
    j: usize,
    weight: f64,
) -> f64 {
    weight
        * op.kinetic_diagonal()
            .iter()
            .enumerate()
            .map(|(i, t)| t * x[(i, j)].norm_sqr())
            .sum::<f64>()
}
/// Weighted nonlocal band energy of one k block (collective).
fn band_nonlocal<K: Collective>(
    op: &mut dft_physics::distributed_complex::DistributedComplexHamiltonian,
    comm: &K,
    x: &dft_physics::data::ComplexMatrix,
    weights: &[f64],
) -> Result<f64, ScfError> {
    let mut nl = vec![0.; weights.len()];
    op.nonlocal_band_energies(comm, x, &mut nl)
        .map_err(numerical)?;
    Ok(nl.iter().zip(weights).map(|(e, f)| e * f).sum::<f64>())
}
fn double_counting_local(dv: f64, density: &[f64], hxc: &[f64]) -> f64 {
    dv * density.iter().zip(hxc).map(|(n, v)| n * v).sum::<f64>()
}
/// Rank-local partials of the map's energy integrals, in the fixed slot order
/// kinetic, local, double counting, ionic mean, hxc difference, vh difference.
fn integral_locals(
    kinetic: f64,
    dv: f64,
    density: &[f64],
    input: &[f64],
    ionic: &[f64],
    points: usize,
    fields: &Fields,
) -> [f64; 6] {
    [
        kinetic,
        dv * density.iter().zip(ionic).map(|(n, v)| n * v).sum::<f64>(),
        double_counting_local(dv, density, &fields.hxc),
        ionic.iter().sum::<f64>() / points as f64,
        dv * density
            .iter()
            .zip(input)
            .zip(&fields.hxc)
            .map(|((a, b), v)| (a - b) * v)
            .sum::<f64>(),
        dv * density
            .iter()
            .zip(input)
            .zip(&fields.vh)
            .map(|((a, b), v)| (a - b) * v)
            .sum::<f64>(),
    ]
}

super::checkpoint::record!(PerK {
    eigenvalues,
    residuals,
    converged,
    calls
});
super::checkpoint::record!(Solution {
    eigenvectors,
    hx,
    per_k
});
super::checkpoint::record!(Map {
    solve,
    warm_start,
    density,
    residuals,
    energy,
    timings,
    input_fields
});

impl Map {
    pub(super) fn validate_checkpoint(
        &self,
        state: &KStartState,
        planes: usize,
    ) -> Result<(), ScfError> {
        self.input_fields.validate_checkpoint(planes)?;
        if self.density.len() != planes
            || self
                .density
                .iter()
                .chain(&self.residuals)
                .any(|v| !v.is_finite())
            || self.solve.eigenvectors.0.len() != state.blocks.len()
            || self.solve.per_k.len() != state.blocks.len()
        {
            return Err(super::checkpoint::invalid("pending k map shape/values"));
        }
        for (k, (x, b)) in self
            .solve
            .eigenvectors
            .0
            .iter()
            .zip(&state.blocks)
            .enumerate()
        {
            let bands = state.orbitals.0[k].cols();
            let info = &self.solve.per_k[k];
            if x.rows() != b.operator.local_rows().len()
                || x.cols() != bands
                || info.eigenvalues.len() != bands
                || info.residuals.len() != bands
                || info.converged.len() != bands
                || x.as_slice()
                    .iter()
                    .any(|v| !v.re.is_finite() || !v.im.is_finite())
                || info
                    .eigenvalues
                    .iter()
                    .chain(&info.residuals)
                    .any(|v| !v.is_finite())
            {
                return Err(super::checkpoint::invalid("pending k orbital shape/values"));
            }
        }
        if let Some(hx) = &self.solve.hx
            && (hx.0.len() != state.blocks.len()
                || hx.0.iter().zip(&self.solve.eigenvectors.0).any(|(h, x)| {
                    h.rows() != x.rows()
                        || h.cols() != x.cols()
                        || h.as_slice()
                            .iter()
                            .any(|v| !v.re.is_finite() || !v.im.is_finite())
                }))
        {
            return Err(super::checkpoint::invalid("pending k HX shape/values"));
        }
        Ok(())
    }
}
