//! Explicit final-cutoff reference initialization. No session, cutoff ladder,
//! serial Hamiltonian, serial band-density workspace or global orbital block.
//!
//! Inputs must be identical on every communicator rank. All field values are
//! owned; replicated basis/grid *metadata* is retained. See the Start section
//! of docs/design/reference-scf.md for the payload-accounting boundary.
use crate::{Functional, ScfError, ScfOptions, SpeciesPotential};
use dft_basis::{FftGridSpec, HalfSphere, PlaneWaveBasis};
use dft_geometry::{Atom, Cell};
use dft_mpi::Collective;
use dft_numerics::{RealGemmWorkspace, RealMatrix, RealOp, pool::WorkerPool};
use dft_physics::{
    data::{Array3, Complex64 as C},
    distributed_real_gamma::{DistributedRealHamiltonian, owned_projector_indices},
    ewald::{self, EwaldOptions, Ion},
    hamiltonian::{NonlocalOperator, Species, structure_factor},
};
use dft_solvers::distributed_fft::DistributedPairFft;
use dft_upf::radial::{projector_transform_with_rule, real_spherical_harmonic, spherical_bessel};
use std::{
    collections::BTreeMap,
    f64::consts::{PI, SQRT_2},
    sync::Arc,
};

fn err(e: impl std::fmt::Display) -> ScfError {
    ScfError::Numerical(e.to_string())
}

/// Parsed, immutable Gamma input. This API deliberately does not accept an
/// already-constructed ScfProblem (which has paid the replicated setup cost).
/// The caller explicitly selects this API; no environment/profile is changed.
pub struct StartInput<'a> {
    pub cell: &'a Cell,
    pub atoms: &'a [Atom],
    pub potentials: &'a [SpeciesPotential<'a>],
    pub options: &'a ScfOptions,
}

/// Per-rank principal array payloads, in bytes, NOT allocator/RSS measurements.
/// The peak covers the named orbital/field buffers; excludes FFT/GEMM backend
/// workspaces, basis/layout metadata, radial caches, allocator and MPI storage.
#[derive(Debug, Clone, Default)]
pub struct StartMemory {
    pub wave_rows: usize,
    pub owned_rows: usize,
    pub grid_points: usize,
    pub owned_grid_coefficients: usize,
    pub owned_plane_points: usize,
    pub atomic_columns: usize,
    pub bands: usize,
    pub complex_projector_bytes: usize,
    pub real_projector_bytes: usize,
    pub retained_field_bytes: usize,
    pub start_block_bytes: usize,
    pub ewald_terms_bytes: usize,
    pub principal_peak_bound_bytes: usize,
    /// Same named buffers in a final-cutoff replicated legacy start; excludes
    /// the coarse problem and counts neither opaque workspaces nor RSS.
    pub legacy_principal_peak_bound_bytes: usize,
    pub legacy_idle_band_density_bytes: usize,
}

/// Owned values at the first Davidson boundary. `density_g` is forward
/// normalized on `fields.reciprocal_indices()` (full FFT sticks, zeros/high-G
/// ringing included). G=0 is exactly N_el/Omega. `density` and `potential` use
/// `fields.plane_offsets()` order; the operator receives its own plane layout.
pub struct StartState {
    pub basis: PlaneWaveBasis,
    pub grid: FftGridSpec,
    pub half: HalfSphere,
    pub operator: DistributedRealHamiltonian,
    pub fields: StartFields,
    pub density: Vec<f64>,
    pub density_g: Vec<C>,
    pub ionic: Vec<f64>,
    pub potential: Vec<f64>,
    /// Kinetic plus nonlocal plane-wave diagonal in Ha, in owned ORC row order.
    /// Add the current mean local potential at each solve. No pair weights.
    pub diagonal_static: Vec<f64>,
    pub orbitals: RealMatrix,
    pub ritz_values_hartree: Vec<f64>,
    pub electrons: f64,
    pub ewald: ewald::EwaldEnergy,
    pub memory: StartMemory,
}

/// Reusable owned full-grid transform. Only the index list used during plan
/// creation is full length; no full-length field values are created here.
pub struct StartFields {
    fft: DistributedPairFft,
    dims: [usize; 3],
}
impl StartFields {
    fn new<K: Collective>(dims: [usize; 3], comm: &K) -> Result<Self, ScfError> {
        let indices: Vec<_> = (0..dims[0])
            .flat_map(|i| (0..dims[1]).flat_map(move |j| (0..dims[2]).map(move |k| [i, j, k])))
            .collect();
        Ok(Self {
            fft: DistributedPairFft::new(dims, &indices, comm.rank(), comm.size()).map_err(err)?,
            dims,
        })
    }
    pub fn reciprocal_indices(&self) -> &[usize] {
        self.fft.local_positions()
    }
    pub fn first_axis(&self) -> usize {
        self.fft.first_axis()
    }
    pub fn plane_offsets(&self) -> Vec<usize> {
        (0..self.fft.plane_count())
            .flat_map(|i| {
                (0..self.fft.rows()).flat_map(move |y| {
                    let start = self.fft.grid_start(i, y);
                    (start..start + self.dims[2]).collect::<Vec<_>>()
                })
            })
            .collect()
    }
    pub fn inverse<K: Collective>(&mut self, comm: &K, coef: &[C]) -> Result<Vec<f64>, ScfError> {
        self.fft.begin_to_real(comm, coef).map_err(err)?;
        let mut out = Vec::with_capacity(self.fft.plane_points());
        let mut bad = 0.;
        for i in 0..self.fft.plane_count() {
            self.fft.to_real_plane(i).map_err(err)?;
            for z in self.fft.plane() {
                if !z.re.is_finite() || !z.im.is_finite() {
                    bad = 1.;
                }
                out.push(z.re);
            }
        }
        agree(comm, bad, "reference start: nonfinite field")?;
        Ok(out)
    }
    pub fn forward<K: Collective>(&mut self, comm: &K, planes: &[f64]) -> Result<Vec<C>, ScfError> {
        agree(
            comm,
            f64::from(planes.len() != self.fft.plane_points()),
            "reference start: plane shape",
        )?;
        let len = self.fft.rows() * self.dims[2];
        for i in 0..self.fft.plane_count() {
            for (z, &v) in self
                .fft
                .plane_mut()
                .iter_mut()
                .zip(&planes[i * len..(i + 1) * len])
            {
                *z = C::new(v, 0.);
            }
            self.fft.to_reciprocal_plane(i).map_err(err)?;
        }
        let mut out = vec![C::default(); self.fft.local_len()];
        self.fft.finish_to_reciprocal(comm, &mut out).map_err(err)?;
        Ok(out)
    }
    /// Explicit diagnostic/delivery boundary, never called by `start`.
    pub fn gather<K: Collective>(&self, comm: &K, planes: &[f64]) -> Result<Array3<f64>, ScfError> {
        let mut full = vec![0.; self.dims.iter().product()];
        self.fft
            .gather_planes(comm, planes, &mut Vec::new(), &mut full)
            .map_err(err)?;
        Array3::from_shape_vec(self.dims, full).map_err(err)
    }
}
fn agree<K: Collective>(comm: &K, bad: f64, message: &'static str) -> Result<(), ScfError> {
    let mut all = [0.];
    comm.tree_sum_f64(&[bad], &mut all).map_err(err)?;
    if all[0] != 0. {
        Err(ScfError::InvalidInput(message))
    } else {
        Ok(())
    }
}
fn sum<K: Collective>(comm: &K, v: &[f64]) -> Result<Vec<f64>, ScfError> {
    let mut out = vec![0.; v.len()];
    comm.tree_sum_f64(v, &mut out).map_err(err)?;
    Ok(out)
}

/// Build the final-cutoff cold start, collective on `comm` (including size 1).
/// Atomic charge/orbitals are required, never silently replaced with uniform
/// density or random waves. Currently LDA-PZ and library LDA-PW without NLCC;
/// PBE, core charge and frozen-strain inputs fail explicitly before setup.
/// Randomization rule: none. Atomic columns retain species/orbital/atom/m
/// order; missing bands are completed by deterministic kinetic-ordered ORC
/// cosine/sine waves, as in the legacy atomic start. Seed is not consumed.
pub fn start<K: Collective>(input: StartInput<'_>, comm: &K) -> Result<StartState, ScfError> {
    let StartInput {
        cell,
        atoms,
        potentials,
        options: o,
    } = input;
    if comm.size() == 0
        || comm.rank() >= comm.size()
        || atoms.is_empty()
        || potentials.is_empty()
        || o.strain_reference.is_some()
        || o.functional == Functional::Pbe
        || o.pw92_prefactor != crate::Pw92Prefactor::Modified
    {
        return Err(ScfError::InvalidInput(
            "reference start: unsupported problem",
        ));
    }
    let mut positions = vec![Vec::new(); potentials.len()];
    for (i, p) in potentials.iter().enumerate() {
        if potentials[..i].iter().any(|q| q.label == p.label)
            || p.potential.core_density().is_some()
            || (!o.allow_functional_override && !o.functional.matches_upf(p.potential.functional()))
            || p.potential
                .atomic_wavefunctions()
                .is_none_or(|w| w.is_empty() || w.iter().any(|c| c.l > 3))
            || p.potential.raw().rhoatom.values.len() != p.potential.mesh().r().len()
            || p.potential.raw().rhoatom.values.iter().all(|&v| v == 0.)
        {
            return Err(ScfError::InvalidInput(
                "reference start: requires compatible NC atomic data without core charge",
            ));
        }
    }
    let mut ions = Vec::with_capacity(atoms.len());
    for atom in atoms {
        let i = potentials
            .iter()
            .position(|p| p.label == atom.species)
            .ok_or(ScfError::InvalidInput("unknown species"))?;
        if atom.fractional_position.iter().any(|r| !r.is_finite()) {
            return Err(ScfError::NonFinite);
        }
        positions[i].push(
            cell.fractional_to_cartesian_bohr(atom.fractional_position.map(|r| r.rem_euclid(1.))),
        );
        ions.push(Ion {
            charge: potentials[i].potential.z_valence(),
            fractional_position: atom.fractional_position,
        });
    }
    let electrons: f64 = ions.iter().map(|i| i.charge).sum();
    if !electrons.is_finite()
        || electrons <= 0.
        || (electrons - electrons.round()).abs() > 1e-10
        || electrons > usize::MAX as f64
        || !(electrons.round() as usize).is_multiple_of(2)
    {
        return Err(ScfError::InvalidInput(
            "reference start: positive even electron count required",
        ));
    }
    let occupied = electrons.round() as usize / 2;
    let bands = o.bands.unwrap_or(
        occupied
            .checked_add(o.extra_bands)
            .ok_or(ScfError::InvalidInput("band overflow"))?,
    );
    let species: Vec<_> = potentials
        .iter()
        .zip(&positions)
        .map(|(p, r)| Species {
            potential: p.potential,
            positions_bohr: r,
        })
        .collect();
    let basis = o.basis(cell, [0.; 3], 1.)?;
    let db = o.basis(cell, [0.; 3], 4.)?;
    let grid = FftGridSpec::new(&basis, &db, o.grid_override)?;
    let dg = FftGridSpec::new(&db, &db, Some(grid.dims()))?;
    let half = HalfSphere::new(&basis, &grid).map_err(err)?;
    let j: usize = species
        .iter()
        .map(|s| {
            s.positions_bohr.len()
                * s.potential
                    .atomic_wavefunctions()
                    .unwrap()
                    .iter()
                    .map(|c| 2 * c.l + 1)
                    .sum::<usize>()
        })
        .sum();
    if bands < occupied || bands > basis.vectors().len() || j > basis.vectors().len() {
        return Err(ScfError::InvalidInput(
            "reference start: incompatible atomic rank/bands",
        ));
    }
    let (ewald, ewald_terms_bytes) = ewald_owned(cell, &ions, comm)?;
    drop(ions);
    let mut fields = StartFields::new(grid.dims(), comm)?;
    let volume = cell.volume_bohr3();
    let points = grid.dims().iter().product::<usize>();
    let indices = owned_projector_indices(&basis, &grid, comm.rank(), comm.size()).map_err(err)?;
    let pool = WorkerPool::new(o.threads).map_err(err)?;
    let nl = NonlocalOperator::new_owned_with_rule(
        &basis,
        volume,
        &species,
        o.radial_even_rule,
        &indices,
    );
    agree(
        comm,
        f64::from(nl.is_err()),
        "reference start: projector construction failed",
    )?;
    let nl = nl?;
    drop(indices);
    // H is applied in chunks of B during the initial Ritz solve; allocating
    // capacity J here would retain atomic-span-sized FFT batches for the run.
    let capacity = bands;
    let mut operator = DistributedRealHamiltonian::new_owned(
        &basis,
        &grid,
        &nl,
        capacity,
        comm.rank(),
        comm.size(),
        pool.clone(),
    )
    .map_err(err)?;
    let projector_count = nl.projector_count();
    let complex_projector_bytes = nl.table_bytes();
    let mut diagonal_static = Vec::with_capacity(operator.local_len());
    for &row in operator.local_rows() {
        let g = if row == 0 {
            0
        } else {
            half.plus()[row.div_ceil(2)] as usize
        };
        let mut value = basis.vectors()[g].norm_squared / 2.;
        let mut offset = 0;
        for atom in nl.atoms() {
            let count = atom.channel_count();
            for j in 0..count {
                for i in 0..count {
                    let bi = nl.owned_value(offset + i, g).expect("owned representative");
                    let bj = nl.owned_value(offset + j, g).expect("owned representative");
                    value += atom.coupling_hartree()[i + count * j] * (bi * bj.conj()).re;
                }
            }
            offset += count;
        }
        diagonal_static.push(value);
    }
    drop(nl);
    let coeff = atomic_fields(&fields, &db, &dg, volume, &species, o);
    agree(
        comm,
        f64::from(coeff.is_err()),
        "reference start: atomic fields failed",
    )?;
    let (rho_hat, ion_hat) = coeff?;
    drop(db);
    drop(dg);
    let ionic = fields.inverse(comm, &ion_hat)?;
    let mut density = fields.inverse(comm, &rho_hat)?;
    drop(rho_hat);
    drop(ion_hat);
    for x in &mut density {
        *x = x.max(0.);
    }
    let charge = sum(comm, &[density.iter().sum::<f64>()])?[0] * volume / points as f64;
    if !charge.is_finite() || charge <= 0. {
        return Err(ScfError::InvalidInput(
            "reference start: invalid atomic charge",
        ));
    }
    for x in &mut density {
        *x *= electrons / charge;
    }
    let mut density_g = fields.forward(comm, &density)?;
    for (z, &g) in density_g.iter_mut().zip(fields.reciprocal_indices()) {
        if g == 0 {
            *z = C::new(electrons / volume, 0.);
        }
    }
    let b = cell.reciprocal_columns_inverse_bohr();
    let dims = grid.dims();
    let vh_hat: Vec<_> = density_g
        .iter()
        .zip(fields.reciprocal_indices())
        .map(|(&z, &g)| {
            let ijk = [g / (dims[1] * dims[2]), g / dims[2] % dims[1], g % dims[2]];
            let m: [f64; 3] = std::array::from_fn(|a| {
                if ijk[a] <= (dims[a] - 1) / 2 {
                    ijk[a] as f64
                } else {
                    ijk[a] as f64 - dims[a] as f64
                }
            });
            let cart: [f64; 3] = std::array::from_fn(|a| (0..3).map(|k| b[a][k] * m[k]).sum());
            let g2: f64 = cart.iter().map(|v| v * v).sum();
            if g == 0 {
                C::default()
            } else {
                z * (4. * PI / g2)
            }
        })
        .collect();
    let mut potential = fields.inverse(comm, &vh_hat)?;
    drop(vh_hat);
    let xc = local_xc(o.functional, &density);
    agree(comm, f64::from(xc.is_err()), "reference start: XC failed")?;
    let xc = xc?;
    for ((v, &ion), &x) in potential.iter_mut().zip(&ionic).zip(&xc) {
        *v = ion + (*v + x);
    }
    drop(xc);
    let mut map = crate::grid_dist::PlaneMap::new(
        dims,
        comm.rank(),
        comm.size(),
        fields.first_axis(),
        operator.plane_layout().0,
    );
    let mut op_planes = Vec::new();
    map.a_to_b(comm, &potential, &mut op_planes)?;
    operator.set_owned_potential(&op_planes, 0).map_err(err)?;
    drop(op_planes);
    drop(map);
    let (orbitals, ritz_values_hartree) = atomic_block(
        cell,
        &basis,
        &half,
        &species,
        o,
        bands,
        &mut operator,
        comm,
        pool,
    )?;
    let l = operator.local_len();
    let f = fields.fft.plane_points();
    let g = fields.fft.local_len();
    let n = basis.vectors().len();
    let k = j.max(bands);
    let retained_field_bytes = 8 * (3 * f + points + l) + 16 * g;
    let start_block_bytes = 8 * l * bands;
    // Named array live upper envelope; Q + next-Q/HQ + raw backup + selected,
    // Gram/EVD/rotation matrices, owned fields and both projector formats.
    let principal_peak_bound_bytes = (retained_field_bytes
        + 8 * l * (3 * k + bands)
        + 64 * k * k
        + 8 * l * projector_count
        + complex_projector_bytes
        + 32 * g
        + 16 * f)
        .max(ewald_terms_bytes);
    let legacy_idle_band_density_bytes = 8 * points * (bands + 1);
    let legacy_principal_peak_bound_bytes = (legacy_idle_band_density_bytes
        + 16 * n * (2 * k + bands)
        + 8 * n * bands
        + 64 * k * k
        + 64 * points
        + complex_projector_bytes
        + 8 * l * projector_count)
        .max(ewald_terms_bytes);
    let memory = StartMemory {
        wave_rows: n,
        owned_rows: l,
        grid_points: points,
        owned_grid_coefficients: g,
        owned_plane_points: f,
        atomic_columns: j,
        bands,
        complex_projector_bytes,
        real_projector_bytes: 8 * l * projector_count,
        retained_field_bytes,
        start_block_bytes,
        ewald_terms_bytes,
        principal_peak_bound_bytes,
        legacy_principal_peak_bound_bytes,
        legacy_idle_band_density_bytes,
    };
    Ok(StartState {
        basis,
        grid,
        half,
        operator,
        fields,
        density,
        density_g,
        ionic,
        potential,
        diagonal_static,
        orbitals,
        ritz_values_hartree,
        electrons,
        ewald,
        memory,
    })
}

fn local_xc(functional: Functional, n: &[f64]) -> Result<Vec<f64>, ScfError> {
    if n.is_empty() {
        return Ok(Vec::new());
    }
    match functional {
        Functional::LdaPz => n
            .iter()
            .map(|&v| {
                dft_physics::xc::lda_pz(v)
                    .map(|p| p.v_xc())
                    .map_err(ScfError::from)
            })
            .collect(),
        Functional::LdaPw => Ok(dft_physics::spin::unpolarized_pw(
            &Array3::from_shape_vec([n.len(), 1, 1], n.to_vec()).map_err(err)?,
            1.,
        )?
        .potential
        .into_raw_vec_and_offset()
        .0),
        Functional::Pbe => Err(ScfError::InvalidInput(
            "reference start: PBE requires distributed gradients",
        )),
    }
}

fn atomic_fields(
    fields: &StartFields,
    db: &PlaneWaveBasis,
    dg: &FftGridSpec,
    volume: f64,
    species: &[Species<'_>],
    o: &ScfOptions,
) -> Result<(Vec<C>, Vec<C>), ScfError> {
    let dims = dg.dims();
    // Only owned lookup entries: scan the already-built density metadata once.
    // Do not allocate another full-density-sphere map on every rank.
    let owned: BTreeMap<_, _> = fields
        .reciprocal_indices()
        .iter()
        .enumerate()
        .map(|(slot, &flat)| (flat, slot))
        .collect();
    let represented: Vec<_> = dg
        .g_to_fft()
        .iter()
        .zip(db.vectors())
        .filter_map(|(&[i, j, k], g)| {
            owned
                .get(&((i * dims[1] + j) * dims[2] + k))
                .map(|&slot| (slot, g))
        })
        .collect();
    drop(owned);
    let mut rho = vec![C::default(); fields.fft.local_len()];
    let mut ion = rho.clone();
    for s in species {
        let mesh = s.potential.mesh();
        let samples = &s.potential.raw().rhoatom.values;
        let local = s
            .potential
            .local_transform_with_conventions(volume, o.local_mesh, o.radial_even_rule)
            .map_err(err)?;
        let mut cache = BTreeMap::new();
        let mut work = vec![0.; samples.len()];
        for &(i, g) in &represented {
            let (r, v) = if let Some(&v) = cache.get(&g.norm_squared.to_bits()) {
                v
            } else {
                let q = g.norm_squared.sqrt();
                for ((v, &r), &f) in work.iter_mut().zip(mesh.r()).zip(samples) {
                    *v = if q == 0. {
                        f
                    } else {
                        f * spherical_bessel(0, q * r).map_err(err)?
                    };
                }
                let r = mesh.integrate(&work).map_err(err)? / volume;
                let v = if q == 0. {
                    local.alpha_z()
                } else {
                    local.at_g(q)
                }
                .map_err(err)?;
                cache.insert(g.norm_squared.to_bits(), (r, v));
                (r, v)
            };
            // Preserve the charge builder's atom accumulation order.
            for &pos in s.positions_bohr {
                rho[i] += structure_factor(g.cartesian, &[pos])? * r;
            }
            ion[i] += v * structure_factor(g.cartesian, s.positions_bohr)?;
        }
    }
    Ok((rho, ion))
}

fn gram<K: Collective>(
    comm: &K,
    ws: &mut RealGemmWorkspace,
    a: &RealMatrix,
    b: &RealMatrix,
) -> Result<RealMatrix, ScfError> {
    let mut g = RealMatrix::zeros(a.cols(), b.cols()).map_err(err)?;
    if a.rows() > 0 {
        ws.multiply_into(a.as_ref(), RealOp::Transpose, b.as_ref(), g.as_mut())
            .map_err(err)?;
    }
    let out = sum(comm, g.as_slice())?;
    g.as_mut_slice().copy_from_slice(&out);
    Ok(g)
}
fn rotate(
    ws: &mut RealGemmWorkspace,
    a: &RealMatrix,
    t: &RealMatrix,
) -> Result<RealMatrix, ScfError> {
    let mut q = RealMatrix::zeros(a.rows(), t.cols()).map_err(err)?;
    if a.rows() > 0 {
        ws.multiply_into(a.as_ref(), RealOp::Identity, t.as_ref(), q.as_mut())
            .map_err(err)?;
    }
    Ok(q)
}

/// Owned form of the legacy guarded CholQR2: two small Gram reductions and
/// local GEMMs. Reuses its column scaling, pivot guards and factorization.
/// A rejected pivot takes two-pass distributed MGS with the legacy threshold.
fn orthogonalize<K: Collective>(
    comm: &K,
    ws: &mut RealGemmWorkspace,
    raw: &RealMatrix,
) -> Result<RealMatrix, ScfError> {
    let mut q = raw.clone();
    for first in [true, false] {
        let g = gram(comm, ws, &q, &q)?;
        let Some(t) = crate::AtomicStartReport::owned_qr_transform(&g, first)? else {
            return mgs(comm, raw);
        };
        q = rotate(ws, &q, &t)?;
    }
    Ok(q)
}

fn mgs<K: Collective>(comm: &K, raw: &RealMatrix) -> Result<RealMatrix, ScfError> {
    let mut q = raw.clone();
    let l = q.rows();
    for j in 0..q.cols() {
        let orig = sum(
            comm,
            &[q.as_slice()[j * l..(j + 1) * l].iter().map(|v| v * v).sum()],
        )?[0]
            .sqrt();
        if !orig.is_finite() || orig == 0. {
            return Err(ScfError::InvalidInput(
                "rank deficient atomic starting orbitals",
            ));
        }
        for i in 0..l {
            q[(i, j)] /= orig;
        }
        for _ in 0..2 {
            for k in 0..j {
                let dot = sum(comm, &[(0..l).map(|i| q[(i, k)] * q[(i, j)]).sum()])?[0];
                for i in 0..l {
                    q[(i, j)] -= dot * q[(i, k)];
                }
            }
        }
        let norm = sum(comm, &[(0..l).map(|i| q[(i, j)] * q[(i, j)]).sum()])?[0].sqrt();
        if !norm.is_finite() || norm <= dft_numerics::linalg::DEPENDENCE_TOLERANCE {
            return Err(ScfError::InvalidInput(
                "rank deficient atomic starting orbitals",
            ));
        }
        for i in 0..l {
            q[(i, j)] /= norm;
        }
    }
    Ok(q)
}

#[allow(clippy::too_many_arguments)]
fn atomic_block<K: Collective>(
    cell: &Cell,
    basis: &PlaneWaveBasis,
    half: &HalfSphere,
    species: &[Species<'_>],
    o: &ScfOptions,
    bands: usize,
    op: &mut DistributedRealHamiltonian,
    comm: &K,
    pool: Arc<WorkerPool>,
) -> Result<(RealMatrix, Vec<f64>), ScfError> {
    let rows = op.local_rows().to_vec();
    let l = rows.len();
    let j: usize = species
        .iter()
        .map(|s| {
            s.positions_bohr.len()
                * s.potential
                    .atomic_wavefunctions()
                    .unwrap()
                    .iter()
                    .map(|c| 2 * c.l + 1)
                    .sum::<usize>()
        })
        .sum();
    let mut raw = RealMatrix::zeros(l, j).map_err(err)?;
    let mut col = 0;
    let filled = (|| -> Result<(), ScfError> {
        for sp in species {
            for c in sp.potential.atomic_wavefunctions().unwrap() {
                let mut radial = BTreeMap::new();
                for &pos in sp.positions_bohr {
                    for m in -(c.l as isize)..=c.l as isize {
                        let il = C::new(0., -1.).powu(c.l as u32);
                        for (i, &row) in rows.iter().enumerate() {
                            let g = if row == 0 {
                                0
                            } else {
                                half.plus()[row.div_ceil(2)] as usize
                            };
                            let g = &basis.momenta()[g];
                            let r = if let Some(&v) = radial.get(&g.norm_squared.to_bits()) {
                                v
                            } else {
                                let v = projector_transform_with_rule(
                                    sp.potential.mesh(),
                                    c.l,
                                    &c.radial_times_r,
                                    g.norm_squared.sqrt(),
                                    o.radial_even_rule,
                                )
                                .map_err(err)?
                                    / cell.volume_bohr3().sqrt();
                                radial.insert(g.norm_squared.to_bits(), v);
                                v
                            };
                            let y = if g.norm_squared == 0. {
                                if c.l == 0 { 1. / (4. * PI).sqrt() } else { 0. }
                            } else {
                                real_spherical_harmonic(c.l, m, g.cartesian).map_err(err)?
                            };
                            let z = il * structure_factor(g.cartesian, &[pos])? * (r * y);
                            raw[(i, col)] = if row == 0 {
                                z.re
                            } else if row % 2 == 1 {
                                SQRT_2 * z.re
                            } else {
                                SQRT_2 * z.im
                            };
                        }
                        col += 1;
                    }
                }
            }
        }
        Ok(())
    })();
    agree(
        comm,
        f64::from(filled.is_err()),
        "reference start: atomic orbital transform failed",
    )?;
    filled?;
    // All ranks agree on a bad local fill before the Gram reductions.
    agree(
        comm,
        f64::from(raw.as_slice().iter().any(|v| !v.is_finite())),
        "nonfinite atomic span",
    )?;
    let mut ws = RealGemmWorkspace::with_pool(l, j.max(bands), pool).map_err(err)?;
    let mut q = orthogonalize(comm, &mut ws, &raw)?;
    drop(raw);
    if j < bands {
        let mut data = RealMatrix::zeros(l, bands).map_err(err)?;
        data.as_mut_slice()[..l * j].copy_from_slice(q.as_slice());
        let mut count = j;
        // Legacy completion skips the later full-sphere partner; its sign is
        // immaterial to the resulting subspace but retained here explicitly.
        let mut order: Vec<_> = (0..basis.vectors().len()).collect();
        order.sort_by(|&a, &b| {
            basis.momenta()[a]
                .norm_squared
                .total_cmp(&basis.momenta()[b].norm_squared)
        });
        for g in order {
            if g > half.mirror()[g] as usize {
                continue;
            }
            for component in 0..if g == 0 { 1 } else { 2 } {
                if count == bands {
                    break;
                }
                let h = if g == 0 {
                    0
                } else {
                    half.plus()
                        .iter()
                        .position(|&p| p as usize == g || half.mirror()[p as usize] as usize == g)
                        .unwrap()
                };
                let target = if g == 0 { 0 } else { 2 * h - 1 + component };
                for (i, &row) in rows.iter().enumerate() {
                    data[(i, count)] = if row == target {
                        if component == 1 && half.plus()[h] as usize != g {
                            -1.
                        } else {
                            1.
                        }
                    } else {
                        0.
                    };
                }
                for _ in 0..2 {
                    for k in 0..count {
                        let dot = sum(
                            comm,
                            &[(0..l).map(|i| data[(i, k)] * data[(i, count)]).sum()],
                        )?[0];
                        for i in 0..l {
                            data[(i, count)] -= dot * data[(i, k)];
                        }
                    }
                }
                let norm = sum(comm, &[(0..l).map(|i| data[(i, count)].powi(2)).sum()])?[0].sqrt();
                if norm > dft_numerics::linalg::DEPENDENCE_TOLERANCE {
                    for i in 0..l {
                        data[(i, count)] /= norm;
                    }
                    count += 1;
                }
            }
            if count == bands {
                break;
            }
        }
        if count != bands {
            return Err(ScfError::InvalidInput("insufficient atomic starting rank"));
        }
        q = data;
    }
    let k = q.cols();
    let mut small = RealMatrix::zeros(k, k).map_err(err)?;
    for first in (0..k).step_by(bands) {
        let width = bands.min(k - first);
        let x = RealMatrix::new(
            l,
            width,
            q.as_slice()[first * l..(first + width) * l].to_vec(),
        )
        .map_err(err)?;
        let mut hx = RealMatrix::zeros(l, width).map_err(err)?;
        op.apply_block(comm, &x, &mut hx).map_err(err)?;
        let part = gram(comm, &mut ws, &q, &hx)?;
        small.as_mut_slice()[first * k..(first + width) * k].copy_from_slice(part.as_slice());
    }
    for j in 0..k {
        for i in 0..=j {
            let v = 0.5 * (small[(i, j)] + small[(j, i)]);
            small[(i, j)] = v;
            small[(j, i)] = v;
        }
    }
    let (mut values, rot) = dft_numerics::real::symmetric_eigh(&small).map_err(err)?;
    let selected = RealMatrix::new(k, bands, rot.as_slice()[..k * bands].to_vec()).map_err(err)?;
    let out = rotate(&mut ws, &q, &selected)?;
    values.truncate(bands);
    Ok((out, values))
}

// Same candidate partition and chained summation as the adopted legacy Ewald
// split. Its local term arrays are dropped before fields/operator allocation.
fn ewald_owned<K: Collective>(
    cell: &Cell,
    ions: &[Ion],
    comm: &K,
) -> Result<(ewald::EwaldEnergy, usize), ScfError> {
    if comm.size() == 1 {
        return Ok((ewald::energy(cell, ions, EwaldOptions::default())?, 0));
    }
    let plan = ewald::EwaldPlan::new(cell, ions, EwaldOptions::default())?;
    let block = |n| {
        let counts = dft_mpi::block_counts(n, comm.size());
        let first = counts[..comm.rank()].iter().sum::<usize>();
        first..first + counts[comm.rank()]
    };
    let (mut real, mut reciprocal) = (Vec::new(), Vec::new());
    let failed = plan
        .real_terms(block(plan.real_candidates()), &mut real)
        .is_err();
    if !failed {
        plan.reciprocal_terms(block(plan.reciprocal_candidates()), &mut reciprocal);
    }
    let bytes = 8 * (real.len() + reciprocal.len());
    let mut state = [0.; 3];
    for root in 0..comm.size() {
        if root == comm.rank() {
            if failed {
                state[2] = 1.;
            } else if state[2] == 0. {
                state[0] = ewald::accumulate(state[0], &real);
                state[1] = ewald::accumulate(state[1], &reciprocal);
            }
        }
        comm.bcast_f64(&mut state, root).map_err(err)?;
    }
    if state[2] != 0. {
        return Err(ScfError::InvalidInput("reference start: Ewald failed"));
    }
    Ok((plan.finish_sums(state[0], state[1])?, bytes))
}
