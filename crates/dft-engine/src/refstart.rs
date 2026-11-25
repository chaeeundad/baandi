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
use dft_numerics::{GemmSite, RealGemmWorkspace, RealMatrix, RealOp, pool::WorkerPool};
use dft_physics::{
    data::{Array3, Complex64 as C},
    distributed_real_gamma::{DistributedRealHamiltonian, owned_projector_indices},
    ewald::{self, EwaldOptions, Ion},
    hamiltonian::{NonlocalOperator, Species, structure_factor},
};
use dft_solvers::distributed_fft::{DistributedPairFft, PlanCache};
use dft_upf::radial::{projector_transform_with_rule, real_spherical_harmonic, spherical_bessel};
use std::{
    collections::BTreeMap,
    f64::consts::{PI, SQRT_2},
    sync::Arc,
    time::Instant,
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
    pub setup: crate::refscf::record::SetupTimings,
    pub basis: PlaneWaveBasis,
    pub grid: FftGridSpec,
    pub half: HalfSphere,
    pub operator: DistributedRealHamiltonian,
    pub fields: StartFields,
    pub density: Vec<f64>,
    pub density_g: Vec<C>,
    pub ionic: Vec<f64>,
    /// Immutable NLCC density on the owned field planes, in electrons/Bohr^3.
    pub core_density: Option<Vec<f64>>,
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
impl<K: Collective> dft_physics::gga::distributed::Fields<K> for StartFields {
    type Error = ScfError;
    fn dims(&self) -> [usize; 3] {
        self.dims
    }
    fn reciprocal_indices(&self) -> &[usize] {
        self.reciprocal_indices()
    }
    fn plane_points(&self) -> usize {
        self.fft.plane_points()
    }
    fn forward(&mut self, comm: &K, n: &[f64]) -> Result<Vec<C>, ScfError> {
        self.forward(comm, n)
    }
    fn inverse(&mut self, comm: &K, n: &[C]) -> Result<Vec<f64>, ScfError> {
        self.inverse(comm, n)
    }
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
    pub(crate) fn shared_fft(&self) -> Option<DistributedPairFft> {
        self.fft.shared_workspace()
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
        self.inverse_checked(comm, coef, false)
    }
    fn inverse_checked<K: Collective>(
        &mut self,
        comm: &K,
        coef: &[C],
        core: bool,
    ) -> Result<Vec<f64>, ScfError> {
        self.fft.begin_to_real(comm, coef).map_err(err)?;
        let mut out = Vec::with_capacity(self.fft.plane_points());
        let mut bad = 0.;
        let mut re_max = 1_f64;
        let mut im_max = 0_f64;
        for i in 0..self.fft.plane_count() {
            self.fft.to_real_plane(i).map_err(err)?;
            for z in self.fft.plane() {
                if !z.re.is_finite() || !z.im.is_finite() {
                    bad = 1.;
                }
                if core {
                    re_max = re_max.max(z.re.abs());
                    im_max = im_max.max(z.im.abs());
                }
                out.push(z.re);
            }
        }
        agree(comm, bad, "reference start: nonfinite field")?;
        if core {
            // Match the legacy core's global reality guard without gathering a
            // field. Empty plane owners still participate in the scalar gather.
            let mut maxima = vec![0.; 2 * comm.size()];
            comm.allgatherv_f64(&[re_max, im_max], &vec![2; comm.size()], &mut maxima)
                .map_err(err)?;
            let scale = maxima
                .as_chunks::<2>()
                .0
                .iter()
                .map(|p| p[0])
                .fold(1., f64::max);
            if maxima
                .as_chunks::<2>()
                .0
                .iter()
                .any(|p| p[1] > 1e-12 * scale)
            {
                return Err(ScfError::InvalidInput(
                    "reference start: non-real core density",
                ));
            }
        }
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

/// Timing adapter used only at setup/delivery boundaries. It delegates exactly
/// one operation, without barriers or changes to the reduction association.
pub(crate) struct TimedCollective<'a, K> {
    inner: &'a K,
    seconds: std::cell::Cell<f64>,
}
impl<'a, K> TimedCollective<'a, K> {
    pub(crate) fn new(inner: &'a K) -> Self {
        Self {
            inner,
            seconds: std::cell::Cell::new(0.),
        }
    }
    pub(crate) fn seconds(&self) -> f64 {
        self.seconds.get()
    }
    fn timed<T>(&self, f: impl FnOnce() -> T) -> T {
        let clock = Instant::now();
        let out = f();
        self.seconds
            .set(self.seconds.get() + clock.elapsed().as_secs_f64());
        out
    }
}
impl<K: Collective> Collective for TimedCollective<'_, K> {
    fn rank(&self) -> usize {
        self.inner.rank()
    }
    fn size(&self) -> usize {
        self.inner.size()
    }
    fn barrier(&self) -> Result<(), dft_mpi::MpiError> {
        self.timed(|| self.inner.barrier())
    }
    fn tree_sum_f64(&self, local: &[f64], out: &mut [f64]) -> Result<(), dft_mpi::MpiError> {
        self.timed(|| self.inner.tree_sum_f64(local, out))
    }
    fn allgatherv_f64(
        &self,
        local: &[f64],
        counts: &[usize],
        out: &mut [f64],
    ) -> Result<(), dft_mpi::MpiError> {
        self.timed(|| self.inner.allgatherv_f64(local, counts, out))
    }
    fn bcast_f64(&self, data: &mut [f64], root: usize) -> Result<(), dft_mpi::MpiError> {
        self.timed(|| self.inner.bcast_f64(data, root))
    }
    fn alltoall_f64(
        &self,
        send: &[f64],
        block: usize,
        recv: &mut [f64],
    ) -> Result<(), dft_mpi::MpiError> {
        self.timed(|| self.inner.alltoall_f64(send, block, recv))
    }
}

/// Build the final-cutoff cold start, collective on `comm` (including size 1).
/// Atomic charge/orbitals are required, never silently replaced with uniform
/// density or random waves. LDA-PZ, library LDA-PW (both with optional NLCC)
/// and distributed PBE (without NLCC for now: PBE + core charge fails
/// explicitly until the core gradient is wired, ticket G2B-PBE-NLCC).
/// Randomization rule: none. Atomic columns retain species/orbital/atom/m
/// order; missing bands are completed by deterministic kinetic-ordered ORC
/// cosine/sine waves, as in the legacy atomic start. Seed is not consumed.
pub fn start<K: Collective>(input: StartInput<'_>, comm: &K) -> Result<StartState, ScfError> {
    start_prolonged(input, comm, None)
}
pub(crate) fn start_prolonged<K: Collective>(
    input: StartInput<'_>,
    comm: &K,
    seed: Option<crate::refscf::ladder::Seed>,
) -> Result<StartState, ScfError> {
    let clock = Instant::now();
    let share = dft_numerics::switches::on(dft_numerics::switches::Id::REFSTART_FFT_SHARE);
    let partitions =
        dft_numerics::switches::on(dft_numerics::switches::Id::REFSTART_PARTITION_CACHE);
    let cache = (share || partitions).then(|| PlanCache::with_partition_cache(share, partitions));
    let _plans = cache.as_ref().map(PlanCache::enter);
    let timed = TimedCollective::new(comm);
    let mut state = start_inner(input, &timed, seed, false)?;
    state.setup.total_s = clock.elapsed().as_secs_f64();
    state.setup.collective_s = timed.seconds();
    Ok(state)
}
fn start_inner<K: Collective>(
    input: StartInput<'_>,
    comm: &K,
    seed: Option<crate::refscf::ladder::Seed>,
    fields_only: bool,
) -> Result<StartState, ScfError> {
    let mut setup = crate::refscf::record::SetupTimings::default();
    let clock = Instant::now();
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
        || o.pw92_prefactor != crate::Pw92Prefactor::Modified
    {
        return Err(ScfError::InvalidInput(
            "reference start: unsupported problem",
        ));
    }
    let mut positions = vec![Vec::new(); potentials.len()];
    for (i, p) in potentials.iter().enumerate() {
        if potentials[..i].iter().any(|q| q.label == p.label)
            || (!o.allow_functional_override && !o.functional.matches_upf(p.potential.functional()))
            || p.potential
                .atomic_wavefunctions()
                .is_none_or(|w| w.is_empty() || w.iter().any(|c| c.l > 3))
            || p.potential.raw().rhoatom.values.len() != p.potential.mesh().r().len()
            || p.potential.raw().rhoatom.values.iter().all(|&v| v == 0.)
        {
            return Err(ScfError::InvalidInput(
                "reference start: requires compatible NC atomic charge and l<=3 orbitals",
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
    setup.record("validate_species", clock);
    let clock = Instant::now();
    let basis = o.basis(cell, [0.; 3], 1.)?;
    let db = o.basis(cell, [0.; 3], 4.)?;
    let grid = FftGridSpec::new(&basis, &db, o.grid_override)?;
    let compact_g = dft_numerics::switches::on(dft_numerics::switches::Id::REFSTART_COMPACT_G);
    let dg = if compact_g {
        None
    } else {
        Some(FftGridSpec::new(&db, &db, Some(grid.dims()))?)
    };
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
    if bands < occupied
        || (!fields_only && (bands > basis.vectors().len() || j > basis.vectors().len()))
    {
        return Err(ScfError::InvalidInput(
            "reference start: incompatible atomic rank/bands",
        ));
    }
    setup.record("basis_gvectors", clock);
    let clock = Instant::now();
    let (ewald, ewald_terms_bytes) = ewald_owned(cell, &ions, comm)?;
    drop(ions);
    setup.record("ewald", clock);
    let clock = Instant::now();
    let mut fields = StartFields::new(grid.dims(), comm)?;
    let volume = cell.volume_bohr3();
    let points = grid.dims().iter().product::<usize>();
    setup.record("field_fft_plan", clock);
    let clock = Instant::now();
    let indices = owned_projector_indices(&basis, &grid, comm.rank(), comm.size()).map_err(err)?;
    let pool = WorkerPool::new(o.threads).map_err(err)?;
    setup.record("distribution_owned_rows", clock);
    let clock = Instant::now();
    let nl = NonlocalOperator::new_owned_with_rule(
        &basis,
        volume,
        if fields_only { &[] } else { &species },
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
    setup.record("nonlocal_tables", clock);
    let clock = Instant::now();
    // H is applied in chunks of B during the initial Ritz solve; allocating
    // capacity J here would retain atomic-span-sized FFT batches for the run.
    let capacity = if fields_only { 1 } else { bands };
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
    setup.record("operator_fft_plans_and_tables", clock);
    let clock = Instant::now();
    let projector_count = nl.projector_count();
    let complex_projector_bytes = nl.table_bytes();
    let mut diagonal_static = Vec::with_capacity(operator.local_len());
    let reuse_diagonal =
        dft_numerics::switches::on(dft_numerics::switches::Id::REFSTART_DIAGONAL_REUSE);
    // One atom's owned projector values, not another G-by-projector table.
    let mut beta = if reuse_diagonal {
        vec![
            C::default();
            nl.atoms()
                .iter()
                .map(|a| a.channel_count())
                .max()
                .unwrap_or(0)
        ]
    } else {
        Vec::new()
    };
    let mut previous = None;
    for &row in operator.local_rows() {
        let g = if row == 0 {
            0
        } else {
            half.plus()[row.div_ceil(2)] as usize
        };
        // The real and imaginary ORC rows of one representative have the
        // same diagonal. Compare indices explicitly, including empty owners.
        if reuse_diagonal && previous == Some(g) {
            diagonal_static.push(*diagonal_static.last().unwrap());
            continue;
        }
        let mut value = basis.vectors()[g].norm_squared / 2.;
        let mut offset = 0;
        for atom in nl.atoms() {
            let count = atom.channel_count();
            if reuse_diagonal {
                for (i, b) in beta[..count].iter_mut().enumerate() {
                    *b = nl.owned_value(offset + i, g).expect("owned representative");
                }
            }
            for j in 0..count {
                for i in 0..count {
                    let (bi, bj) = if reuse_diagonal {
                        (beta[i], beta[j])
                    } else {
                        (
                            nl.owned_value(offset + i, g).expect("owned representative"),
                            nl.owned_value(offset + j, g).expect("owned representative"),
                        )
                    };
                    value += atom.coupling_hartree()[i + count * j] * (bi * bj.conj()).re;
                }
            }
            offset += count;
        }
        diagonal_static.push(value);
        previous = Some(g);
    }
    drop(beta);
    drop(nl);
    setup.record("nonlocal_diagonal", clock);
    let clock = Instant::now();
    let coeff = atomic_fields(&fields, &db, dg.as_ref(), volume, &species, o);
    agree(
        comm,
        f64::from(coeff.is_err()),
        "reference start: atomic fields failed",
    )?;
    let AtomicFields {
        rho: rho_hat,
        ion: ion_hat,
        core: core_hat,
    } = coeff?;
    drop(db);
    drop(dg);
    setup.record("atomic_local_tables_structure_factors", clock);
    let clock = Instant::now();
    let ionic = fields.inverse(comm, &ion_hat)?;
    drop(ion_hat);
    if o.functional == Functional::Pbe && core_hat.is_some() {
        return Err(ScfError::InvalidInput(
            "reference start: PBE with core charge (NLCC) is not supported yet",
        ));
    }
    let core_density = core_hat
        .map(|c| fields.inverse_checked(comm, &c, true))
        .transpose()?;
    let prolonged = seed
        .map(|s| s.inject(&basis, &grid, &half, operator.local_rows(), &fields, comm))
        .transpose()?;
    let (density, density_g, warm_orbitals) = if let Some((density_g, orbitals)) = prolonged {
        drop(rho_hat);
        let density = fields.inverse(comm, &density_g)?;
        (density, density_g, Some(orbitals))
    } else {
        let mut density = fields.inverse(comm, &rho_hat)?;
        drop(rho_hat);
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
        (density, density_g, None)
    };
    setup.record("atomic_density", clock);
    let clock = Instant::now();
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
    // Prolongation is not clipped or rescaled. Signed Fourier roundoff uses
    // the same absolute-density LDA evaluation as a reference SCF map.
    let xc = if o.functional == Functional::Pbe {
        dft_physics::gga::distributed::evaluate(&mut fields, comm, cell, &density, &density_g)
            .map(|x| x.potential)
    } else if let Some(core) = &core_density {
        local_xc(
            o.functional,
            &density
                .iter()
                .zip(core)
                .map(|(n, c)| (n + c).abs())
                .collect::<Vec<_>>(),
        )
    } else if warm_orbitals.is_some() {
        local_xc(
            o.functional,
            &density.iter().map(|v| v.abs()).collect::<Vec<_>>(),
        )
    } else {
        local_xc(o.functional, &density)
    };
    agree(comm, f64::from(xc.is_err()), "reference start: XC failed")?;
    let xc = xc?;
    for ((v, &ion), &x) in potential.iter_mut().zip(&ionic).zip(&xc) {
        *v = ion + (*v + x);
    }
    drop(xc);
    setup.record("initial_potential", clock);
    let clock = Instant::now();
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
    setup.record("potential_distribution", clock);
    let clock = Instant::now();
    let (orbitals, ritz_values_hartree) = if fields_only {
        (RealMatrix::zeros(0, bands).map_err(err)?, Vec::new())
    } else if let Some(orbitals) = warm_orbitals {
        (orbitals, Vec::new())
    } else {
        atomic_block(
            cell,
            &basis,
            &half,
            &species,
            o,
            bands,
            &mut operator,
            comm,
            pool,
            &mut setup,
        )?
    };
    setup.record("atomic_block_total", clock);
    let l = operator.local_len();
    let f = fields.fft.plane_points();
    let g = fields.fft.local_len();
    let n = basis.vectors().len();
    let k = j.max(bands);
    let core_bytes = core_density.as_ref().map_or(0, |c| 8 * c.len());
    let retained_field_bytes = 8 * (3 * f + points + l) + 16 * g + core_bytes;
    let start_block_bytes = 8 * l * bands;
    // Named array live upper envelope; Q + next-Q/HQ + raw backup + selected,
    // Gram/EVD/rotation matrices, owned fields and both projector formats.
    let phase_scratch_bytes =
        if dft_numerics::switches::on(dft_numerics::switches::Id::REFSCF_SETUP_PHASE_REUSE) {
            16 * l.max(atoms.len())
        } else {
            0
        };
    let principal_peak_bound_bytes = phase_scratch_bytes
        + (retained_field_bytes
            + 8 * l * (3 * k + bands)
            + 64 * k * k
            + 8 * l * projector_count
            + complex_projector_bytes
            + 32 * g
            + 16 * f
            + if core_density.is_some() {
                16 * g + 8 * f
            } else {
                0
            })
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
        setup,
        basis,
        grid,
        half,
        operator,
        fields,
        density,
        density_g,
        ionic,
        core_density,
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
            "reference start: PBE must use the distributed GGA evaluator",
        )),
    }
}

struct AtomicFields {
    rho: Vec<C>,
    ion: Vec<C>,
    core: Option<Vec<C>>,
}

fn core_radial(s: &Species<'_>, q: f64, volume: f64, o: &ScfOptions) -> Result<f64, ScfError> {
    s.potential.core_density().map_or(Ok(0.), |core| {
        dft_upf::radial::spherical_density_transform_with_rule(
            s.potential.mesh(),
            core,
            q,
            o.radial_even_rule,
        )
        .map(|v| v / volume)
        .map_err(err)
    })
}

fn atomic_fields(
    fields: &StartFields,
    db: &PlaneWaveBasis,
    dg: Option<&FftGridSpec>,
    volume: f64,
    species: &[Species<'_>],
    o: &ScfOptions,
) -> Result<AtomicFields, ScfError> {
    let mut core = species
        .iter()
        .any(|s| s.potential.core_density().is_some())
        .then(|| vec![C::default(); fields.fft.local_len()]);
    let dims = fields.dims;
    let represented: Vec<_> = if let Some(dg) = dg {
        let owned: BTreeMap<_, _> = fields
            .reciprocal_indices()
            .iter()
            .enumerate()
            .map(|(slot, &flat)| (flat, slot))
            .collect();
        dg.g_to_fft()
            .iter()
            .zip(db.vectors())
            .filter_map(|(&[i, j, k], g)| {
                owned
                    .get(&((i * dims[1] + j) * dims[2] + k))
                    .map(|&slot| (slot, g))
            })
            .collect()
    } else {
        // The wave grid already validated the complete density sphere. Reuse
        // its dimensions and the FFT's sorted local positions; no F-sized
        // reverse map or S_r-entry tree is needed to find the used G values.
        let mut represented = Vec::new();
        for g in db.vectors() {
            let [i, j, k] = FftGridSpec::map_miller(g.miller, dims)?;
            let flat = (i * dims[1] + j) * dims[2] + k;
            if let Ok(slot) = fields.reciprocal_indices().binary_search(&flat) {
                represented.push((slot, g));
            }
        }
        represented
    };
    if dft_numerics::switches::on(dft_numerics::switches::Id::REFSTART_LOCAL_SHELLS) {
        // Geometry is common to every species. Build owned shell IDs and sqrt
        // values once, keeping first-encounter ordering and exact norm bits.
        let mut shell_ids = BTreeMap::new();
        let mut q = Vec::new();
        let shells: Vec<_> = represented
            .iter()
            .map(|&(_, g)| {
                *shell_ids
                    .entry(g.norm_squared.to_bits())
                    .or_insert_with(|| {
                        q.push(g.norm_squared.sqrt());
                        q.len() - 1
                    })
            })
            .collect();
        let mut rho = vec![C::default(); fields.fft.local_len()];
        let mut ion = rho.clone();
        let mut angles = vec![0.; represented.len()];
        let mut factors = vec![C::default(); represented.len()];
        for s in species {
            let mesh = s.potential.mesh();
            let samples = &s.potential.raw().rhoatom.values;
            let local = s
                .potential
                .local_transform_with_conventions(volume, o.local_mesh, o.radial_even_rule)
                .map_err(err)?;
            let mut work = vec![0.; samples.len()];
            let radial: Vec<_> = q
                .iter()
                .map(|&q| {
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
                    Ok((r, v, core_radial(s, q, volume, o)?))
                })
                .collect::<Result<_, ScfError>>()?;
            factors.fill(C::default());
            for &pos in s.positions_bohr {
                // Separate contiguous angle and accumulation passes expose
                // independent G rows to vectorization, without changing the
                // dot-product order or using factored/interpolated phases.
                for (angle, &(_, g)) in angles.iter_mut().zip(&represented) {
                    *angle = g.cartesian.iter().zip(pos).map(|(g, r)| g * r).sum();
                }
                if angles.iter().any(|v| !v.is_finite()) {
                    return Err(err("non-finite atomic structure-factor angle"));
                }
                for (a, (&(i, _), &shell)) in represented.iter().zip(&shells).enumerate() {
                    let phase = C::from_polar(1., -angles[a]);
                    // The old single-atom structure_factor sums from zero.
                    rho[i] += (C::default() + phase) * radial[shell].0;
                    factors[a] += phase;
                }
            }
            for (a, (&(i, _), &shell)) in represented.iter().zip(&shells).enumerate() {
                ion[i] += radial[shell].1 * factors[a];
                if let Some(core) = &mut core {
                    core[i] += radial[shell].2 * factors[a];
                }
            }
        }
        return Ok(AtomicFields { rho, ion, core });
    }
    let reuse = dft_numerics::switches::on(dft_numerics::switches::Id::REFSCF_SETUP_PHASE_REUSE);
    let mut phases = if reuse {
        Vec::with_capacity(
            species
                .iter()
                .map(|s| s.positions_bohr.len())
                .max()
                .unwrap_or(0),
        )
    } else {
        Vec::new()
    };
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
            let (r, v, c) = if let Some(&v) = cache.get(&g.norm_squared.to_bits()) {
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
                let c = core_radial(s, q, volume, o)?;
                cache.insert(g.norm_squared.to_bits(), (r, v, c));
                (r, v, c)
            };
            // Preserve every charge multiply/add and the ionic sum order.
            // Only reuse direct phases: no axis factorization or interpolation.
            if reuse {
                phases.clear();
            }
            for &pos in s.positions_bohr {
                let phase = structure_factor(g.cartesian, &[pos])?;
                rho[i] += phase * r;
                if reuse {
                    phases.push(phase);
                }
            }
            let factor = if reuse {
                phases.iter().copied().sum()
            } else {
                structure_factor(g.cartesian, s.positions_bohr)?
            };
            ion[i] += v * factor;
            if let Some(core) = &mut core {
                core[i] += c * factor;
            }
        }
    }
    Ok(AtomicFields { rho, ion, core })
}

/// REFSTART_MKL: route the atomic-start Gram and rotation products through the
/// dispatching GEMM entries (MKL panels in mkl builds) instead of the direct
/// pure-Rust kernel (dev08 si128: 7 calls, 17 GFLOP per rank on pure_direct).
fn mkl_start() -> bool {
    dft_numerics::switches::on(dft_numerics::switches::Id::REFSTART_MKL)
}

fn gram<K: Collective>(
    comm: &K,
    ws: &mut RealGemmWorkspace,
    a: &RealMatrix,
    b: &RealMatrix,
) -> Result<RealMatrix, ScfError> {
    let mut g = RealMatrix::zeros(a.cols(), b.cols()).map_err(err)?;
    if a.rows() > 0 {
        if mkl_start() {
            let product =
                if dft_numerics::switches::on(dft_numerics::switches::Id::REFSTART_SINGLE_GEMM) {
                    RealGemmWorkspace::multiply_slices_single_into
                } else {
                    RealGemmWorkspace::multiply_slices_into
                };
            product(
                ws,
                GemmSite::Other,
                a.as_slice(),
                a.rows(),
                a.cols(),
                a.rows(),
                RealOp::Transpose,
                b.as_slice(),
                b.cols(),
                b.rows(),
                g.as_mut(),
            )
            .map_err(err)?;
        } else {
            ws.multiply_into(a.as_ref(), RealOp::Transpose, b.as_ref(), g.as_mut())
                .map_err(err)?;
        }
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
        if mkl_start() {
            ws.multiply_slices_wide_into(
                GemmSite::Rotate,
                a.as_slice(),
                a.rows(),
                a.cols(),
                a.rows(),
                RealOp::Identity,
                t.as_slice(),
                t.cols(),
                t.rows(),
                q.as_mut(),
            )
            .map_err(err)?;
        } else {
            ws.multiply_into(a.as_ref(), RealOp::Identity, t.as_ref(), q.as_mut())
                .map_err(err)?;
        }
    }
    Ok(q)
}

/// Replicated small-matrix congruence T^t A T. Do not call gram(): these
/// matrices are identical on every rank and must not be reduced again.
fn congruence(
    ws: &mut RealGemmWorkspace,
    a: &RealMatrix,
    t: &RealMatrix,
) -> Result<RealMatrix, ScfError> {
    let at = rotate(ws, a, t)?;
    let mut out = RealMatrix::zeros(t.cols(), t.cols()).map_err(err)?;
    if mkl_start() {
        ws.multiply_slices_into(
            GemmSite::Other,
            t.as_slice(),
            t.rows(),
            t.cols(),
            t.rows(),
            RealOp::Transpose,
            at.as_slice(),
            at.cols(),
            at.rows(),
            out.as_mut(),
        )
        .map_err(err)?;
    } else {
        ws.multiply_into(t.as_ref(), RealOp::Transpose, at.as_ref(), out.as_mut())
            .map_err(err)?;
    }
    for j in 0..out.cols() {
        for i in 0..j {
            let v = 0.5 * (out[(i, j)] + out[(j, i)]);
            out[(i, j)] = v;
            out[(j, i)] = v;
        }
    }
    Ok(out)
}

/// One distributed S and the existing scaled Cholesky inverse. Check its
/// metric on the small matrices before applying H. Ill-conditioned spans use
/// the unchanged two-pass QR/MGS path, not a truncated generalized problem.
fn overlap_transform<K: Collective>(
    comm: &K,
    ws: &mut RealGemmWorkspace,
    raw: &RealMatrix,
    setup: &mut crate::refscf::record::SetupTimings,
) -> Result<Option<RealMatrix>, ScfError> {
    let clock = Instant::now();
    let s = gram(comm, ws, raw, raw)?;
    setup.record("atomic_overlap_gram", clock);
    let clock = Instant::now();
    let t = crate::AtomicStartReport::owned_qr_transform(&s, true)?;
    setup.record("atomic_overlap_factor_inverse", clock);
    let Some(t) = t else { return Ok(None) };
    let clock = Instant::now();
    let metric = congruence(ws, &s, &t)?;
    let good = (0..metric.cols()).all(|j| {
        (0..metric.rows()).all(|i| {
            let defect = (metric[(i, j)] - f64::from(i == j)).abs();
            defect.is_finite() && defect <= dft_numerics::linalg::HERMITIAN_TOLERANCE
        })
    });
    setup.record("atomic_overlap_metric_check", clock);
    Ok(good.then_some(t))
}

/// Owned form of the legacy guarded CholQR2: two small Gram reductions and
/// local GEMMs. Reuses its column scaling, pivot guards and factorization.
/// A rejected pivot takes two-pass distributed MGS with the legacy threshold.
fn orthogonalize<K: Collective>(
    comm: &K,
    ws: &mut RealGemmWorkspace,
    raw: &RealMatrix,
    setup: &mut crate::refscf::record::SetupTimings,
) -> Result<RealMatrix, ScfError> {
    let mut q = raw.clone();
    for first in [true, false] {
        let clock = Instant::now();
        let g = gram(comm, ws, &q, &q)?;
        setup.record("atomic_qr_gram", clock);
        let clock = Instant::now();
        let transform = crate::AtomicStartReport::owned_qr_transform(&g, first)?;
        setup.record("atomic_qr_factor_inverse", clock);
        let Some(t) = transform else {
            let clock = Instant::now();
            let result = mgs(comm, raw);
            setup.record("atomic_qr_mgs_fallback", clock);
            return result;
        };
        let clock = Instant::now();
        q = rotate(ws, &q, &t)?;
        setup.record("atomic_qr_rotate", clock);
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

/// One owned-row amplitude table per species and one direct phase per atom/G.
/// Output ordering remains species, radial orbital, atom, m. This changes only
/// repeated scalar work, not the radial quadrature or complex expression order.
#[allow(clippy::too_many_arguments)]
fn fill_species_orbitals(
    cell: &Cell,
    basis: &PlaneWaveBasis,
    half: &HalfSphere,
    species: &[Species<'_>],
    o: &ScfOptions,
    rows: &[usize],
    raw: &mut RealMatrix,
) -> Result<(), ScfError> {
    let l = rows.len();
    let representatives: Vec<_> = rows
        .iter()
        .map(|&row| {
            if row == 0 {
                0
            } else {
                half.plus()[row.div_ceil(2)] as usize
            }
        })
        .collect();
    let mut phases = vec![C::default(); l];
    let mut angles = vec![0.; l];
    let mut species_offset = 0;
    for sp in species {
        if sp.positions_bohr.is_empty() {
            continue;
        }
        let channels = sp.potential.atomic_wavefunctions().unwrap();
        let width: usize = channels.iter().map(|c| 2 * c.l + 1).sum();
        let mut amplitudes = RealMatrix::zeros(l, width).map_err(err)?;
        let mut channel_offset = 0;
        for c in channels {
            let mut radial = BTreeMap::new();
            for m in -(c.l as isize)..=c.l as isize {
                let column = channel_offset + (m + c.l as isize) as usize;
                for (i, &g) in representatives.iter().enumerate() {
                    if i > 0 && representatives[i - 1] == g {
                        amplitudes[(i, column)] = amplitudes[(i - 1, column)];
                        continue;
                    }
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
                    amplitudes[(i, column)] = r * y;
                }
            }
            channel_offset += 2 * c.l + 1;
        }
        for (atom, &pos) in sp.positions_bohr.iter().enumerate() {
            // Contiguous independent angles; the direct structure_factor dot
            // product and its sum-from-zero convention are preserved exactly.
            for (angle, &g) in angles.iter_mut().zip(&representatives) {
                *angle = basis.momenta()[g]
                    .cartesian
                    .iter()
                    .zip(pos)
                    .map(|(g, r)| g * r)
                    .sum();
            }
            if angles.iter().any(|v| !v.is_finite()) {
                return Err(err("non-finite atomic structure-factor angle"));
            }
            for i in 0..l {
                phases[i] = if i > 0 && representatives[i - 1] == representatives[i] {
                    phases[i - 1]
                } else {
                    C::default() + C::from_polar(1., -angles[i])
                };
            }
            let mut channel_offset = 0;
            for c in channels {
                let count = 2 * c.l + 1;
                let il = C::new(0., -1.).powu(c.l as u32);
                for m in 0..count {
                    let column = species_offset
                        + channel_offset * sp.positions_bohr.len()
                        + atom * count
                        + m;
                    let amplitude = &amplitudes.as_slice()
                        [(channel_offset + m) * l..(channel_offset + m + 1) * l];
                    let out = &mut raw.as_mut_slice()[column * l..(column + 1) * l];
                    for (((out, &amplitude), &phase), &row) in
                        out.iter_mut().zip(amplitude).zip(&phases).zip(rows)
                    {
                        let z = il * phase * amplitude;
                        *out = if row == 0 {
                            z.re
                        } else if row % 2 == 1 {
                            SQRT_2 * z.re
                        } else {
                            SQRT_2 * z.im
                        };
                    }
                }
                channel_offset += count;
            }
        }
        species_offset += width * sp.positions_bohr.len();
    }
    Ok(())
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
    setup: &mut crate::refscf::record::SetupTimings,
) -> Result<(RealMatrix, Vec<f64>), ScfError> {
    let clock = Instant::now();
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
    let reuse = dft_numerics::switches::on(dft_numerics::switches::Id::REFSCF_SETUP_PHASE_REUSE);
    let reuse_orbitals =
        dft_numerics::switches::on(dft_numerics::switches::Id::REFSTART_ORBITAL_REUSE);
    let species_fill =
        dft_numerics::switches::on(dft_numerics::switches::Id::REFSTART_SPECIES_FILL);
    let mut phases = if reuse && !species_fill {
        vec![C::default(); l]
    } else {
        Vec::new()
    };
    let mut raw = RealMatrix::zeros(l, j).map_err(err)?;
    let mut col = 0;
    let filled = (|| -> Result<(), ScfError> {
        if species_fill {
            return fill_species_orbitals(cell, basis, half, species, o, &rows, &mut raw);
        }
        for sp in species {
            for c in sp.potential.atomic_wavefunctions().unwrap() {
                let mut radial = BTreeMap::new();
                // r(|G|) Y_lm(G) is independent of the atom position. Keep
                // exactly the first atom's arithmetic and column order. No
                // interpolation, symmetry tolerance or changed quadrature.
                let cache = reuse_orbitals && sp.positions_bohr.len() > 1;
                let mut amplitudes = if cache {
                    vec![0.; l * (2 * c.l + 1)]
                } else {
                    Vec::new()
                };
                for (atom, &pos) in sp.positions_bohr.iter().enumerate() {
                    if reuse {
                        for (i, &row) in rows.iter().enumerate() {
                            let g = if row == 0 {
                                0
                            } else {
                                half.plus()[row.div_ceil(2)] as usize
                            };
                            phases[i] = structure_factor(basis.momenta()[g].cartesian, &[pos])?;
                        }
                    }
                    for m in -(c.l as isize)..=c.l as isize {
                        let il = C::new(0., -1.).powu(c.l as u32);
                        for (i, &row) in rows.iter().enumerate() {
                            let g = if row == 0 {
                                0
                            } else {
                                half.plus()[row.div_ceil(2)] as usize
                            };
                            let g = &basis.momenta()[g];
                            let slot = (m + c.l as isize) as usize * l + i;
                            let amplitude = if cache && atom > 0 {
                                amplitudes[slot]
                            } else {
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
                                let amplitude = r * y;
                                if cache {
                                    amplitudes[slot] = amplitude;
                                }
                                amplitude
                            };
                            let phase = if reuse {
                                phases[i]
                            } else {
                                structure_factor(g.cartesian, &[pos])?
                            };
                            let z = il * phase * amplitude;
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
    drop(phases);
    // All ranks agree on a bad local fill before the Gram reductions.
    agree(
        comm,
        f64::from(raw.as_slice().iter().any(|v| !v.is_finite())),
        "nonfinite atomic span",
    )?;
    setup.record("atomic_orbital_fill", clock);
    let clock = Instant::now();
    let generalized =
        dft_numerics::switches::on(dft_numerics::switches::Id::REFSTART_GENERALIZED) && j >= bands;
    let mut ws =
        RealGemmWorkspace::with_pool(if generalized { l.max(j) } else { l }, j.max(bands), pool)
            .map_err(err)?;
    let cholesky = dft_numerics::switches::on(dft_numerics::switches::Id::REFSTART_CHOLESKY);
    if generalized && cholesky {
        let clock = Instant::now();
        let s = gram(comm, &mut ws, &raw, &raw)?;
        setup.record("atomic_overlap_gram", clock);
        let clock = Instant::now();
        // Same unit-column pivot guard as the legacy atomic QR transform.
        let factor = dft_numerics::linalg::RealGeneralizedFactor::new(&s, 1e-6).map_err(err)?;
        setup.record("atomic_overlap_factor_solve", clock);
        if let Some(factor) = factor {
            let clock = Instant::now();
            let defect = factor.metric_defect(&s).map_err(err)?;
            setup.record("atomic_overlap_metric_check", clock);
            if defect <= dft_numerics::linalg::HERMITIAN_TOLERANCE {
                let small = projected_h(comm, &mut ws, &raw, bands, op, setup)?;
                let clock = Instant::now();
                let (values, coefficients) = factor.solve_lowest(&small, bands).map_err(err)?;
                let out = rotate(&mut ws, &raw, &coefficients)?;
                setup.record("atomic_generalized_diagonalize_rotate", clock);
                return Ok((out, values));
            }
        }
    }
    if generalized
        && !cholesky
        && let Some(t) = overlap_transform(comm, &mut ws, &raw, setup)?
    {
        let small = projected_h(comm, &mut ws, &raw, bands, op, setup)?;
        let clock = Instant::now();
        let reduced = congruence(&mut ws, &small, &t)?;
        let (mut values, rot) = dft_numerics::real::symmetric_eigh(&reduced).map_err(err)?;
        let selected =
            RealMatrix::new(j, bands, rot.as_slice()[..j * bands].to_vec()).map_err(err)?;
        let coefficients = rotate(&mut ws, &t, &selected)?;
        let out = rotate(&mut ws, &raw, &coefficients)?;
        values.truncate(bands);
        setup.record("atomic_generalized_diagonalize_rotate", clock);
        return Ok((out, values));
    }
    // All ranks see the same reduced S and make the same fallback decision.
    // No rank deletion, looser pivot, or alternate completion is introduced.
    let mut q = orthogonalize(comm, &mut ws, &raw, setup)?;
    drop(raw);
    let completion_clock = Instant::now();
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
    setup.record("atomic_completion", completion_clock);
    setup.record("atomic_orthonormalize_complete", clock);
    let small = projected_h(comm, &mut ws, &q, bands, op, setup)?;
    let clock = Instant::now();
    let k = q.cols();
    let (mut values, rot) = dft_numerics::real::symmetric_eigh(&small).map_err(err)?;
    let selected = RealMatrix::new(k, bands, rot.as_slice()[..k * bands].to_vec()).map_err(err)?;
    let out = rotate(&mut ws, &q, &selected)?;
    values.truncate(bands);
    setup.record("atomic_diagonalize_rotate", clock);
    Ok((out, values))
}

/// H is still applied to every column in the same Ritz space, in B-wide chunks.
fn projected_h<K: Collective>(
    comm: &K,
    ws: &mut RealGemmWorkspace,
    q: &RealMatrix,
    bands: usize,
    op: &mut DistributedRealHamiltonian,
    setup: &mut crate::refscf::record::SetupTimings,
) -> Result<RealMatrix, ScfError> {
    let l = q.rows();
    let clock = Instant::now();
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
        let h_clock = Instant::now();
        op.apply_block(comm, &x, &mut hx).map_err(err)?;
        setup.record("atomic_h_apply", h_clock);
        let gram_clock = Instant::now();
        let part = gram(comm, ws, q, &hx)?;
        setup.record("atomic_h_gram", gram_clock);
        small.as_mut_slice()[first * k..(first + width) * k].copy_from_slice(part.as_slice());
    }
    for j in 0..k {
        for i in 0..=j {
            let v = 0.5 * (small[(i, j)] + small[(j, i)]);
            small[(i, j)] = v;
            small[(j, i)] = v;
        }
    }
    setup.record("atomic_h_apply_gram", clock);
    Ok(small)
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

/// Complex per-k owned blocks. The real Gamma entry above is unchanged.
#[derive(Clone)]
pub struct KOrbitals(pub Vec<dft_numerics::linalg::ComplexMatrix>);
impl KOrbitals {
    pub fn cols(&self) -> usize {
        self.0.first().map_or(0, |x| x.cols())
    }
}
pub struct KBlock {
    pub point: dft_geometry::kpoints::KPoint,
    pub basis: PlaneWaveBasis,
    pub operator: dft_physics::distributed_complex::DistributedComplexHamiltonian,
}
pub struct KStartState {
    pub initial_h_columns: usize,
    pub common: StartState,
    pub blocks: Vec<KBlock>,
    pub orbitals: KOrbitals,
}
impl std::ops::Deref for KStartState {
    type Target = StartState;
    fn deref(&self) -> &StartState {
        &self.common
    }
}
impl std::ops::DerefMut for KStartState {
    fn deref_mut(&mut self) -> &mut StartState {
        &mut self.common
    }
}

/// Initialize a validated (possibly time-reversal-reduced) list. All k points
/// share the density grid; their kinetic/projector supports use k+G. The common
/// atomic-field initializer skips Gamma orbital/projector construction. A small
/// unused Gamma FFT adapter remains in the shared field state; it imposes no
/// Gamma band-rank constraint on the actual k blocks.
pub fn start_kpoints<K: Collective>(
    input: StartInput<'_>,
    points: &[dft_geometry::kpoints::KPoint],
    comm: &K,
) -> Result<KStartState, ScfError> {
    use dft_numerics::linalg::{ComplexMatrix as Matrix, GemmWorkspace, MatrixOp, hermitian_eigh};
    dft_geometry::kpoints::validate(points).map_err(err)?;
    let o = input.options;
    let cell = input.cell;
    let db = o.basis(cell, [0.; 3], 4.)?;
    let bases: Vec<_> = points
        .iter()
        .map(|p| o.basis(cell, p.fractional, 1.))
        .collect::<Result<_, _>>()?;
    let gamma = o.basis(cell, [0.; 3], 1.)?;
    let mut dims = FftGridSpec::new(&gamma, &db, o.grid_override)?.dims();
    for b in &bases {
        let d = FftGridSpec::new(b, &db, o.grid_override)?.dims();
        for a in 0..3 {
            dims[a] = dims[a].max(d[a]);
        }
    }
    let mut common_options = o.clone();
    common_options.gamma_real = crate::GammaReal::On;
    common_options.grid_override = Some(dims);
    let common_clock = Instant::now();
    let timed = TimedCollective::new(comm);
    let mut common = start_inner(
        StartInput {
            cell,
            atoms: input.atoms,
            potentials: input.potentials,
            options: &common_options,
        },
        &timed,
        None,
        true,
    )?;
    common.setup.total_s = common_clock.elapsed().as_secs_f64();
    common.setup.collective_s = timed.seconds();
    let bands = common.orbitals.cols();
    let mut positions = vec![Vec::new(); input.potentials.len()];
    for a in input.atoms {
        let i = input
            .potentials
            .iter()
            .position(|p| p.label == a.species)
            .ok_or(ScfError::InvalidInput("unknown species"))?;
        positions[i].push(
            cell.fractional_to_cartesian_bohr(a.fractional_position.map(|r| r.rem_euclid(1.))),
        );
    }
    let species: Vec<_> = input
        .potentials
        .iter()
        .zip(&positions)
        .map(|(p, r)| Species {
            potential: p.potential,
            positions_bohr: r,
        })
        .collect();
    let pool = WorkerPool::new(o.threads).map_err(err)?;
    let mut blocks = Vec::new();
    let mut orbitals = Vec::new();
    let mut initial_h_columns = 0;
    let clock = Instant::now();
    for (point, basis) in points.iter().zip(bases) {
        if bands > basis.vectors().len() {
            return Err(ScfError::InvalidInput(
                "reference k start: too many bands for a k basis",
            ));
        }
        let grid = FftGridSpec::new(&basis, &db, Some(dims))?;
        let layout = DistributedPairFft::new(dims, grid.g_to_fft(), comm.rank(), comm.size())
            .map_err(err)?;
        let rows = layout.local_positions().to_vec();
        drop(layout);
        let nl = NonlocalOperator::new_owned_with_rule(
            &basis,
            cell.volume_bohr3(),
            &species,
            o.radial_even_rule,
            &rows,
        )?;
        let mut op = dft_physics::distributed_complex::DistributedComplexHamiltonian::new_owned(
            &basis,
            &grid,
            &nl,
            (2 * bands)
                .max(common.memory.atomic_columns)
                .min(basis.vectors().len()),
            comm.rank(),
            comm.size(),
            pool.clone(),
        )
        .map_err(err)?;
        let mut planes = crate::grid_dist::PlaneMap::new(
            dims,
            comm.rank(),
            comm.size(),
            common.fields.first_axis(),
            op.first_axis(),
        );
        let mut potential = Vec::new();
        planes.a_to_b(comm, &common.potential, &mut potential)?;
        op.set_owned_potential(&potential, 1).map_err(err)?;
        let mut columns = Vec::new();
        for sp in &species {
            for channel in sp.potential.atomic_wavefunctions().unwrap() {
                for pos in sp.positions_bohr {
                    for m in -(channel.l as isize)..=channel.l as isize {
                        let mut col = Vec::with_capacity(rows.len());
                        for &g in &rows {
                            let q = &basis.momenta()[g];
                            let radial = projector_transform_with_rule(
                                sp.potential.mesh(),
                                channel.l,
                                &channel.radial_times_r,
                                q.norm_squared.sqrt(),
                                o.radial_even_rule,
                            )
                            .map_err(err)?
                                / cell.volume_bohr3().sqrt();
                            let y = if q.norm_squared == 0. {
                                if channel.l == 0 {
                                    1. / (4. * PI).sqrt()
                                } else {
                                    0.
                                }
                            } else {
                                real_spherical_harmonic(channel.l, m, q.cartesian).map_err(err)?
                            };
                            let angle: f64 = q.cartesian.iter().zip(pos).map(|(q, r)| q * r).sum();
                            col.push(
                                C::new(0., -1.).powu(channel.l as u32)
                                    * C::from_polar(1., -angle)
                                    * (radial * y),
                            );
                        }
                        columns.push(col);
                    }
                }
            }
        }
        // Two-pass owned-row MGS, then lowest kinetic coordinate complements.
        // Complements are chosen from the actual k+G kinetic order, never G index.
        let mut order: Vec<_> = (0..basis.vectors().len()).collect();
        order.sort_by(|&a, &b| {
            basis.momenta()[a]
                .norm_squared
                .total_cmp(&basis.momenta()[b].norm_squared)
                .then(a.cmp(&b))
        });
        let target = columns.len().max(bands).min(basis.vectors().len());
        let mut q = Matrix::zeros(rows.len(), 0).map_err(err)?;
        q.reserve_columns(target).map_err(err)?;
        for mut col in columns.into_iter().chain(order.into_iter().map(|g| {
            rows.iter()
                .map(|&r| C::new(f64::from(r == g), 0.))
                .collect()
        })) {
            for _ in 0..2 {
                for j in 0..q.cols() {
                    let z: C = q
                        .column(j)
                        .iter()
                        .zip(&col)
                        .map(|(a, b)| a.conj() * b)
                        .sum();
                    let z = sum(comm, &[z.re, z.im])?;
                    let z = C::new(z[0], z[1]);
                    for (a, b) in col.iter_mut().zip(q.column(j)) {
                        *a -= z * b;
                    }
                }
            }
            let norm = sum(comm, &[col.iter().map(|z| z.norm_sqr()).sum()])?[0].sqrt();
            if norm <= dft_numerics::linalg::DEPENDENCE_TOLERANCE {
                continue;
            }
            let j = q.cols();
            q.resize_columns(j + 1).map_err(err)?;
            for (a, b) in q.column_mut(j).iter_mut().zip(col) {
                *a = b / norm;
            }
            if q.cols() == target {
                break;
            }
        }
        if q.cols() < bands {
            return Err(ScfError::InvalidInput(
                "reference complex atomic start has insufficient rank",
            ));
        }
        let mut hq = Matrix::zeros(rows.len(), q.cols()).map_err(err)?;
        op.apply_block(comm, &q, &mut hq).map_err(err)?;
        initial_h_columns += q.cols();
        let mut gemm = GemmWorkspace::with_pool(
            rows.len().max(q.cols()),
            q.cols(),
            o.threads.max(1),
            pool.clone(),
        )
        .map_err(err)?;
        let mut h = Matrix::zeros(q.cols(), q.cols()).map_err(err)?;
        gemm.multiply_into(&q, MatrixOp::Adjoint, &hq, MatrixOp::Identity, &mut h)
            .map_err(err)?;
        let local: Vec<_> = h.as_slice().iter().flat_map(|z| [z.re, z.im]).collect();
        let all = sum(comm, &local)?;
        for (z, p) in h.as_mut_slice().iter_mut().zip(all.as_chunks::<2>().0) {
            *z = C::new(p[0], p[1]);
        }
        for j in 0..h.cols() {
            h[(j, j)].im = 0.;
            for i in 0..j {
                let z = (h[(i, j)] + h[(j, i)].conj()) * 0.5;
                h[(i, j)] = z;
                h[(j, i)] = z.conj();
            }
        }
        let (_, mut c) = hermitian_eigh(&h).map_err(err)?;
        c.resize_columns(bands).map_err(err)?;
        let mut x = Matrix::zeros(rows.len(), bands).map_err(err)?;
        gemm.multiply_into(&q, MatrixOp::Identity, &c, MatrixOp::Identity, &mut x)
            .map_err(err)?;
        orbitals.push(x);
        blocks.push(KBlock {
            point: *point,
            basis,
            operator: op,
        });
    }
    common.setup.record("complex_k_atomic_starts", clock);
    common.orbitals = RealMatrix::zeros(0, bands).map_err(err)?;
    Ok(KStartState {
        initial_h_columns,
        common,
        blocks,
        orbitals: KOrbitals(orbitals),
    })
}
