//! Explicit final-cutoff reference initialization. No session, cutoff ladder,
//! serial Hamiltonian, serial band-density workspace or global orbital block.
//!
//! Inputs must be identical on every communicator rank. All field values are
//! owned; replicated basis/grid *metadata* is retained. See the Start section
//! of docs/design/reference-scf.md for the payload-accounting boundary.
use crate::{
    Functional,
    model::SpeciesPotential,
    model::{ScfError, ScfOptions},
};
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
use dft_upf::radial::{
    RadialTransformTable, projector_transform_with_rule, real_spherical_harmonic, spherical_bessel,
};
use std::{
    collections::BTreeMap,
    f64::consts::{PI, SQRT_2},
    sync::Arc,
    time::Instant,
};

fn err(e: impl std::fmt::Display) -> ScfError {
    ScfError::Numerical(e.to_string())
}

pub(crate) fn validate_strain_grid(options: &ScfOptions) -> Result<(), ScfError> {
    if options.strain_reference.is_some() && options.grid_override.is_none() {
        return Err(ScfError::InvalidInput(
            "strain reference requires a frozen FFT grid",
        ));
    }
    Ok(())
}

// Atomic data are optional initial guesses, independent of the NC operator.
fn atomic_orbitals(
    p: &dft_upf::Pseudopotential,
) -> impl Iterator<Item = &dft_upf::AtomicWavefunction> {
    p.atomic_wavefunctions()
        .unwrap_or(&[])
        .iter()
        .filter(|c| c.l <= 3)
}

fn has_atomic_charge(p: &dft_upf::Pseudopotential) -> bool {
    let rho = &p.raw().rhoatom.values;
    rho.len() == p.mesh().r().len() && rho.iter().any(|&v| v != 0.)
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

/// Origin of the density and orbital guesses, independent of SCF convergence.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct StartReport {
    pub density: DensityStartKind,
    pub orbitals: OrbitalStartKind,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DensityStartKind {
    Atomic,
    Uniform,
    Supplied,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OrbitalStartKind {
    Atomic,
    AtomicWithKineticCompletion,
    Kinetic,
    Supplied,
}
fn orbital_kind(columns: usize, bands: usize) -> OrbitalStartKind {
    if columns == 0 {
        OrbitalStartKind::Kinetic
    } else if columns < bands {
        OrbitalStartKind::AtomicWithKineticCompletion
    } else {
        OrbitalStartKind::Atomic
    }
}

/// Owned values at the first Davidson boundary. `density_g` is forward
/// normalized on `fields.reciprocal_indices()` (full FFT sticks, zeros/high-G
/// ringing included). G=0 is exactly N_el/Omega. `density` and `potential` use
/// `fields.plane_offsets()` order; the operator receives its own plane layout.
pub struct StartState {
    /// Run-owned unpolarized PBE convention, shared by startup and every map.
    pub pw92_prefactor: crate::Pw92Prefactor,
    pub report: StartReport,
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
    /// PBE-only core Fourier cache and reusable combined-density inputs.
    pub pbe_core: Option<PbeCore>,
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

/// Owned NLCC inputs for PBE. The core coefficients are prepared once from
/// the real core used by XC, so the spectral derivative matches the legacy path.
pub struct PbeCore {
    density_g: Vec<C>,
    total_density: Vec<f64>,
    total_density_g: Vec<C>,
}
impl PbeCore {
    /// Normalized core coefficients in `StartFields::reciprocal_indices` order.
    pub fn density_g(&self) -> &[C] {
        &self.density_g
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn evaluate<K: Collective>(
        &mut self,
        prefactor: crate::Pw92Prefactor,
        fields: &mut StartFields,
        comm: &K,
        cell: &Cell,
        density: &[f64],
        density_g: &[C],
        core: &[f64],
    ) -> Result<dft_physics::gga::distributed::XcSlab, ScfError> {
        agree(
            comm,
            f64::from(
                density.len() != self.total_density.len()
                    || core.len() != density.len()
                    || density_g.len() != self.density_g.len(),
            ),
            "reference PBE: core/valence field shape mismatch",
        )?;
        for ((total, n), c) in self.total_density.iter_mut().zip(density).zip(core) {
            *total = n + c;
        }
        for ((total, n), c) in self
            .total_density_g
            .iter_mut()
            .zip(density_g)
            .zip(&self.density_g)
        {
            *total = n + c;
        }
        dft_physics::gga::distributed::evaluate_pw92(
            prefactor,
            fields,
            comm,
            cell,
            &self.total_density,
            &self.total_density_g,
        )
    }
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
    fn real_indices(&self) -> Option<Vec<usize>> {
        Some(self.plane_offsets())
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
/// and distributed PBE with optional NLCC (including the core gradient).
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

    let cache = PlanCache::with_partition_cache(true, true);
    let _plans = cache.enter();
    let timed = TimedCollective::new(comm);
    let mut state = start_inner(input, &timed, seed, false, false)?;
    state.setup.total_s = clock.elapsed().as_secs_f64();
    state.setup.collective_s = timed.seconds();
    Ok(state)
}
fn start_inner<K: Collective>(
    input: StartInput<'_>,
    comm: &K,
    seed: Option<crate::refscf::ladder::Seed>,
    fields_only: bool,
    general_occupations: bool,
) -> Result<StartState, ScfError> {
    let mut setup = crate::refscf::record::SetupTimings::default();
    let clock = Instant::now();
    let StartInput {
        cell,
        atoms,
        potentials,
        options: o,
    } = input;
    validate_strain_grid(o)?;
    if comm.size() == 0
        || comm.rank() >= comm.size()
        || atoms.is_empty()
        || potentials.is_empty()
        || (o.pw92_prefactor == crate::Pw92Prefactor::Original && o.functional != Functional::Pbe)
    {
        return Err(ScfError::InvalidInput(
            "reference start: unsupported problem",
        ));
    }
    let mut positions = vec![Vec::new(); potentials.len()];
    for (i, p) in potentials.iter().enumerate() {
        if potentials[..i].iter().any(|q| q.label == p.label)
            || (!o.allow_functional_override && !o.functional.matches_upf(p.potential.functional()))
        {
            return Err(ScfError::InvalidInput(
                "reference start: duplicate species or incompatible functional",
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
        || (!general_occupations && !(electrons.round() as usize).is_multiple_of(2))
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

    let half = HalfSphere::new(&basis, &grid).map_err(err)?;
    let j: usize = species
        .iter()
        .map(|s| {
            s.positions_bohr.len()
                * atomic_orbitals(s.potential)
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

    // One atom's owned projector values, not another G-by-projector table.
    let mut beta = {
        vec![
            C::default();
            nl.atoms()
                .iter()
                .map(|a| a.channel_count())
                .max()
                .unwrap_or(0)
        ]
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
        if previous == Some(g) {
            diagonal_static.push(*diagonal_static.last().unwrap());
            continue;
        }
        let mut value = basis.vectors()[g].norm_squared / 2.;
        let mut offset = 0;
        for atom in nl.atoms() {
            let count = atom.channel_count();
            {
                for (i, b) in beta[..count].iter_mut().enumerate() {
                    *b = nl.owned_value(offset + i, g).expect("owned representative");
                }
            }
            for j in 0..count {
                for i in 0..count {
                    let (bi, bj) = { (beta[i], beta[j]) };
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
    let coeff = atomic_fields(&fields, &db, volume, &species, o);
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
    setup.record("atomic_local_tables_structure_factors", clock);
    let clock = Instant::now();
    let ionic = fields.inverse(comm, &ion_hat)?;
    drop(ion_hat);
    let core_density = core_hat
        .map(|c| fields.inverse_checked(comm, &c, true))
        .transpose()?;
    let mut pbe_core = if o.functional == Functional::Pbe {
        core_density
            .as_ref()
            .map(|core| {
                let density_g = fields.forward(comm, core)?;
                Ok::<_, ScfError>(PbeCore {
                    total_density_g: vec![C::default(); density_g.len()],
                    total_density: vec![0.; core.len()],
                    density_g,
                })
            })
            .transpose()?
    } else {
        None
    };
    let supplied = seed.is_some();
    let report = StartReport {
        density: if supplied {
            DensityStartKind::Supplied
        } else if species.iter().any(|s| !has_atomic_charge(s.potential)) {
            DensityStartKind::Uniform
        } else {
            DensityStartKind::Atomic
        },
        orbitals: if supplied {
            OrbitalStartKind::Supplied
        } else {
            orbital_kind(j, bands)
        },
    };
    let prolonged = seed
        .map(|s| s.inject(&basis, &grid, &half, operator.local_rows(), &fields, comm))
        .transpose()?;
    let (density, density_g, warm_orbitals) = if let Some((density_g, orbitals)) = prolonged {
        drop(rho_hat);
        let density = fields.inverse(comm, &density_g)?;
        (density, density_g, Some(orbitals))
    } else if report.density == DensityStartKind::Uniform {
        drop(rho_hat);
        // Match the legacy all-species fallback; keep the neutral zero mode exact.
        let density = vec![electrons / volume; fields.fft.plane_points()];
        let density_g = fields
            .reciprocal_indices()
            .iter()
            .map(|&g| C::new(if g == 0 { electrons / volume } else { 0. }, 0.))
            .collect();
        (density, density_g, None)
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
    let _xc_scope = dft_physics::xc::screening_scope("start", false, None);
    let screening = dft_physics::xc::screen_negative();
    let screened_total = screening.then(|| {
        density
            .iter()
            .enumerate()
            .map(|(i, n)| n + core_density.as_ref().map_or(0., |c| c[i]))
            .collect::<Vec<_>>()
    });
    if let Some(total) = &screened_total {
        dft_physics::gga::distributed::audit(
            &fields,
            comm,
            total,
            cell.volume_bohr3() / dims.iter().product::<usize>() as f64,
            "total",
            Some(electrons),
        )?;
    }
    let _local_scope = dft_physics::xc::screening_scope("start", screening, Some(electrons));
    let xc = if o.functional != Functional::Pbe && screening {
        local_xc(o.functional, screened_total.as_ref().unwrap())
    } else if o.functional == Functional::Pbe {
        if let Some(core) = &mut pbe_core {
            core.evaluate(
                o.pw92_prefactor,
                &mut fields,
                comm,
                cell,
                &density,
                &density_g,
                core_density.as_ref().unwrap(),
            )
        } else {
            dft_physics::gga::distributed::evaluate_pw92(
                o.pw92_prefactor,
                &mut fields,
                comm,
                cell,
                &density,
                &density_g,
            )
        }
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
    dft_physics::xc::agree_error(
        comm,
        xc.as_ref()
            .err()
            .map(|e| format!("reference start: XC failed: {e}")),
    )?;
    let xc = xc?;
    for ((v, &ion), &x) in potential.iter_mut().zip(&ionic).zip(&xc) {
        *v = ion + (*v + x);
    }
    drop(xc);
    drop(_local_scope);
    drop(_xc_scope);
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
    let core_bytes = core_density.as_ref().map_or(0, |c| 8 * c.len())
        + pbe_core.as_ref().map_or(0, |c| {
            16 * (c.density_g.len() + c.total_density_g.len()) + 8 * c.total_density.len()
        });
    let retained_field_bytes = 8 * (3 * f + points + l) + 16 * g + core_bytes;
    let start_block_bytes = 8 * l * bands;
    // Named array live upper envelope; Q + next-Q/HQ + raw backup + selected,
    // Gram/EVD/rotation matrices, owned fields and both projector formats.
    let phase_scratch_bytes = { 16 * l.max(atoms.len()) };
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
        report,
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
        pbe_core,
        pw92_prefactor: o.pw92_prefactor,
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
    volume: f64,
    species: &[Species<'_>],
    o: &ScfOptions,
) -> Result<AtomicFields, ScfError> {
    let mut core = species
        .iter()
        .any(|s| s.potential.core_density().is_some())
        .then(|| vec![C::default(); fields.fft.local_len()]);
    let dims = fields.dims;
    let represented: Vec<_> = {
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

    let mut phases = {
        Vec::with_capacity(
            species
                .iter()
                .map(|s| s.positions_bohr.len())
                .max()
                .unwrap_or(0),
        )
    };
    let mut rho = vec![C::default(); fields.fft.local_len()];
    let mut ion = rho.clone();
    for s in species {
        let mesh = s.potential.mesh();
        let samples = &s.potential.raw().rhoatom.values;
        let atomic_charge = has_atomic_charge(s.potential);
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
                let r = if atomic_charge {
                    mesh.integrate(&work).map_err(err)? / volume
                } else {
                    0.
                };
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
            {
                phases.clear();
            }
            for &pos in s.positions_bohr {
                let phase = structure_factor(g.cartesian, &[pos])?;
                rho[i] += phase * r;
                {
                    phases.push(phase);
                }
            }
            let factor: C = phases.iter().copied().sum();
            ion[i] += v * factor;
            if let Some(core) = &mut core {
                core[i] += c * factor;
            }
        }
    }
    Ok(AtomicFields { rho, ion, core })
}

fn gram<K: Collective>(
    comm: &K,
    ws: &mut RealGemmWorkspace,
    a: &RealMatrix,
    b: &RealMatrix,
) -> Result<RealMatrix, ScfError> {
    let mut g = RealMatrix::zeros(a.cols(), b.cols()).map_err(err)?;
    if a.rows() > 0 {
        {
            let product = { RealGemmWorkspace::multiply_slices_single_into };
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
        {
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
        }
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
    setup: &mut crate::refscf::record::SetupTimings,
) -> Result<RealMatrix, ScfError> {
    let mut q = raw.clone();
    for first in [true, false] {
        let clock = Instant::now();
        let g = gram(comm, ws, &q, &q)?;
        setup.record("atomic_qr_gram", clock);
        let clock = Instant::now();
        let transform = crate::model::AtomicStartReport::owned_qr_transform(&g, first)?;
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
        let channels: Vec<_> = atomic_orbitals(sp.potential).collect();
        let width: usize = channels.iter().map(|c| 2 * c.l + 1).sum();
        let mut amplitudes = RealMatrix::zeros(l, width).map_err(err)?;
        let mut channel_offset = 0;
        for c in &channels {
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
            for c in &channels {
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
                * atomic_orbitals(s.potential)
                    .map(|c| 2 * c.l + 1)
                    .sum::<usize>()
        })
        .sum();
    let mut raw = RealMatrix::zeros(l, j).map_err(err)?;
    let filled = fill_species_orbitals(cell, basis, half, species, o, &rows, &mut raw);
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
    setup.record("atomic_orbital_fill", clock);
    let clock = Instant::now();
    let generalized = j >= bands;
    let mut ws =
        RealGemmWorkspace::with_pool(if generalized { l.max(j) } else { l }, j.max(bands), pool)
            .map_err(err)?;

    if generalized {
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

    // All ranks see the same reduced S and make the same fallback decision.
    // No rank deletion, looser pivot, or alternate completion is introduced.
    let mut q = if j == 0 {
        raw.clone()
    } else {
        orthogonalize(comm, &mut ws, &raw, setup)?
    };
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
    /// GK10B (`K_START_KEPT_HX`): per k, the start Ritz orbitals X of the
    /// first channel, their H X (HQ C) and the operator potential it was
    /// applied with; the k driver seeds its kept HX with them so map 1 does
    /// not apply H to the start again. The kept cache hits only while the
    /// orbitals are still bit-identical to X (a seed or ladder replaces
    /// them). Empty otherwise.
    pub start_hx: Vec<(
        dft_numerics::linalg::ComplexMatrix,
        dft_numerics::linalg::ComplexMatrix,
        Vec<f64>,
    )>,
    pub common: StartState,
    /// The k blocks this rank holds: every k point (rows distributed over
    /// the ranks) without pools, the rank's own whole k points with GK9
    /// pools (`pool`, global indices `pool.local()`).
    pub blocks: Vec<KBlock>,
    pub orbitals: KOrbitals,
    /// GK9 k-point pools (`BAANDI_REFSCF_K_POOLS`); `None`: distributed layout.
    pub(crate) pool: Option<crate::refscf::k_pool::KPool>,
}
impl KStartState {
    /// Whether this start uses GK9 k-point pools.
    pub fn pooled(&self) -> bool {
        self.pool.is_some()
    }
    /// Plane-wave basis size of every global k point.
    pub(crate) fn basis_sizes(&self) -> Vec<usize> {
        self.pool.as_ref().map_or_else(
            || {
                self.blocks
                    .iter()
                    .map(|b| b.basis.vectors().len())
                    .collect()
            },
            |p| p.basis_len.clone(),
        )
    }
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
/// Reuse the density grid; their kinetic/projector supports use k+G. The common
/// atomic-field initializer skips Gamma orbital/projector construction. A small
/// unused Gamma FFT adapter remains in the shared field state; it imposes no
/// Gamma band-rank constraint on the actual k blocks.
pub fn start_kpoints<K: Collective>(
    input: StartInput<'_>,
    points: &[dft_geometry::kpoints::KPoint],
    comm: &K,
) -> Result<KStartState, ScfError> {
    Ok(start_kpoints_inner(input, points, comm, false, 1, true)?.remove(0))
}
/// Spin/FD entry: the caller has validated the sector and explicit band count.
/// The existing fixed unpolarized entry retains its even-electron restriction.
pub(crate) fn start_electronic<K: Collective>(
    input: StartInput<'_>,
    points: &[dft_geometry::kpoints::KPoint],
    comm: &K,
    pools: bool,
) -> Result<KStartState, ScfError> {
    Ok(start_kpoints_inner(input, points, comm, true, 1, pools)?.remove(0))
}
/// Collinear spin/FD entry for `channels` identical channel starts. With
/// `BAANDI_K_START_SPIN_SHARE=on` (default) the channel-independent work
/// (bases, grids, projector tables, atomic columns, orthonormalization, start
/// H and Ritz) runs once and each channel gets its own operator and a copy of
/// the orbitals; every channel equals `start_electronic` bit for bit. `off`
/// calls `start_electronic` once per channel. `pools`: the caller supports
/// GK9 k-point pools (used when the switch and the rank/k counts allow).
pub(crate) fn start_electronic_channels<K: Collective>(
    input: StartInput<'_>,
    points: &[dft_geometry::kpoints::KPoint],
    comm: &K,
    channels: usize,
    pools: bool,
) -> Result<Vec<KStartState>, ScfError> {
    if channels > 1 && dft_numerics::switches::on(dft_numerics::switches::Id::K_START_SPIN_SHARE) {
        start_kpoints_inner(input, points, comm, true, channels, pools)
    } else {
        (0..channels)
            .map(|_| {
                start_electronic(
                    StartInput {
                        cell: input.cell,
                        atoms: input.atoms,
                        potentials: input.potentials,
                        options: input.options,
                    },
                    points,
                    comm,
                    pools,
                )
            })
            .collect()
    }
}
/// Automatic FD band growth: append columns to every k block until it holds
/// `bands` orthonormal columns. The candidates are the start's deterministic
/// complements, plane-wave unit vectors in ascending |k+G|^2 order (ties by G
/// index), each orthogonalised twice against the current block (classical
/// Gram-Schmidt with one reduction per pass); dependent candidates are skipped.
/// The existing columns are kept bit for bit. Operator scratch grows to the
/// start's rule of two columns per band (bounded by the basis size).
pub(crate) fn extend_k_orbitals<K: Collective>(
    state: &mut KStartState,
    bands: usize,
    comm: &K,
) -> Result<(), ScfError> {
    use dft_numerics::linalg::ComplexMatrix as Matrix;
    // GK9: a pooled rank holds whole k blocks; their sums are rank-local.
    let comm = &crate::refscf::k_pool::KComm::new(comm, state.pool.is_some());
    for (block, x) in state.blocks.iter_mut().zip(&mut state.orbitals.0) {
        let have = x.cols();
        if bands <= have {
            continue;
        }
        let total = block.basis.vectors().len();
        if bands > total {
            return Err(ScfError::InvalidInput(
                "reference k start: too many bands for a k basis",
            ));
        }
        let rows = block.operator.local_rows().to_vec();
        let mut position = vec![usize::MAX; total];
        for (r, &g) in rows.iter().enumerate() {
            position[g] = r;
        }
        let mut order: Vec<_> = (0..total).collect();
        order.sort_by(|&a, &b| {
            block.basis.momenta()[a]
                .norm_squared
                .total_cmp(&block.basis.momenta()[b].norm_squared)
                .then(a.cmp(&b))
        });
        let mut q = Matrix::zeros(rows.len(), have).map_err(err)?;
        q.reserve_columns(bands).map_err(err)?;
        q.as_mut_slice().copy_from_slice(x.as_slice());
        let mut col = vec![C::default(); rows.len()];
        let mut local = Vec::with_capacity(2 * bands);
        for g in order {
            col.fill(C::default());
            if position[g] != usize::MAX {
                col[position[g]] = C::new(1., 0.);
            }
            for _ in 0..2 {
                local.clear();
                for j in 0..q.cols() {
                    let z: C = q
                        .column(j)
                        .iter()
                        .zip(&col)
                        .map(|(a, b)| a.conj() * b)
                        .sum();
                    local.extend([z.re, z.im]);
                }
                let z = sum(comm, &local)?;
                for (j, z) in z.as_chunks::<2>().0.iter().enumerate() {
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
            for (a, b) in q.column_mut(j).iter_mut().zip(&col) {
                *a = b / norm;
            }
            if q.cols() == bands {
                break;
            }
        }
        if q.cols() < bands {
            return Err(ScfError::InvalidInput(
                "reference band growth: insufficient complement rank",
            ));
        }
        block
            .operator
            .ensure_capacity((2 * bands).min(total))
            .map_err(err)?;
        *x = q;
    }
    Ok(())
}
fn start_kpoints_inner<K: Collective>(
    input: StartInput<'_>,
    points: &[dft_geometry::kpoints::KPoint],
    comm: &K,
    general_occupations: bool,
    channels: usize,
    pools: bool,
) -> Result<Vec<KStartState>, ScfError> {
    dft_geometry::kpoints::validate(points).map_err(err)?;
    let o = input.options;
    let cell = input.cell;
    // Check the caller's grid before common_options installs the k-union grid.
    validate_strain_grid(o)?;
    let bases_clock = Instant::now();
    let db = o.basis(cell, [0.; 3], 4.)?;
    let bases: Vec<_> = points
        .iter()
        .map(|p| o.basis(cell, p.fractional, 1.))
        .collect::<Result<_, _>>()?;
    if general_occupations
        && o.bands
            .is_none_or(|n| n == 0 || bases.iter().any(|b| n > b.vectors().len()))
    {
        return Err(ScfError::InvalidInput(
            "reference spin/FD start requires an explicit band count within every k basis",
        ));
    }
    let gamma = o.basis(cell, [0.; 3], 1.)?;
    let mut dims = FftGridSpec::new(&gamma, &db, o.grid_override)?.dims();
    for b in &bases {
        let d = FftGridSpec::new(b, &db, o.grid_override)?.dims();
        for a in 0..3 {
            dims[a] = dims[a].max(d[a]);
        }
    }
    let bases_s = bases_clock.elapsed().as_secs_f64();
    let mut common_options = o.clone();
    common_options.gamma_real = crate::model::GammaReal::On;
    common_options.grid_override = Some(dims);
    let mut commons = Vec::with_capacity(channels);
    for _ in 0..channels.max(1) {
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
            general_occupations,
        )?;
        common.setup.total_s = common_clock.elapsed().as_secs_f64();
        common.setup.collective_s = timed.seconds();
        *common
            .setup
            .phases
            .entry("complex_k_bases_grids")
            .or_default() += if commons.is_empty() { bases_s } else { 0. };
        commons.push(common);
    }
    let pool = crate::refscf::k_pool::KPool::decide(
        points,
        bases.iter().map(|b| b.vectors().len()).collect(),
        comm,
        pools,
    );
    finish_kpoints(
        input,
        points,
        comm,
        bases,
        &db,
        commons,
        !general_occupations,
        pool,
    )
}

/// Original k atomic columns (`BAANDI_K_START_MEMO=off`): one radial
/// quadrature per (owned row, atom, m). Kept as the bit reference of the memo.
fn atomic_k_columns_direct(
    species: &[Species<'_>],
    basis: &PlaneWaveBasis,
    rows: &[usize],
    volume: f64,
    rule: crate::RadialEvenRule,
    columns: &mut Vec<Vec<C>>,
) -> Result<(), ScfError> {
    for sp in species {
        for channel in atomic_orbitals(sp.potential) {
            for pos in sp.positions_bohr {
                for m in -(channel.l as isize)..=channel.l as isize {
                    let mut col = Vec::with_capacity(rows.len());
                    for &g in rows {
                        let q = &basis.momenta()[g];
                        let radial = projector_transform_with_rule(
                            sp.potential.mesh(),
                            channel.l,
                            &channel.radial_times_r,
                            q.norm_squared.sqrt(),
                            rule,
                        )
                        .map_err(err)?
                            / volume.sqrt();
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
    Ok(())
}

/// Reused buffers of the memoized k column build (GK1).
#[derive(Default)]
struct KColumnScratch {
    /// |k+G|^2 bits -> radial transform / sqrt(volume), per channel.
    radial_by_bits: std::collections::HashMap<u64, f64>,
    /// radial * Y_lm per (m, owned row), m-major.
    radial_y: Vec<f64>,
    /// (-i)^l exp(-i q.tau) per owned row for the current atom.
    phase: Vec<C>,
}

/// GK1: bit-identical to `atomic_k_columns_direct`. The radial transform
/// depends only on (species channel, |k+G|), so it is evaluated once per
/// distinct |k+G|^2 bit pattern (the argument is `norm_squared.sqrt()`, so equal
/// bits give equal arguments); `radial * y` is formed once per (m, row) and the
/// atomic phase once per (atom, row). Every element keeps the original
/// operands and product order `((-i)^l * phase) * (radial * y)`, and columns
/// come out in the original (species, channel, atom, m) order. Species without
/// atoms are skipped, as in the direct loop; on invalid input both return an
/// error (the message may come from a different row).
#[allow(clippy::too_many_arguments)]
fn atomic_k_columns_memo(
    species: &[Species<'_>],
    basis: &PlaneWaveBasis,
    rows: &[usize],
    volume: f64,
    rule: crate::RadialEvenRule,
    tables: Option<&[Vec<RadialTransformTable>]>,
    s: &mut KColumnScratch,
    columns: &mut Vec<Vec<C>>,
) -> Result<(), ScfError> {
    use std::collections::hash_map::Entry;
    let momenta = basis.momenta();
    let n = rows.len();
    for (si, sp) in species.iter().enumerate() {
        if sp.positions_bohr.is_empty() {
            continue;
        }
        for (ci, channel) in atomic_orbitals(sp.potential).enumerate() {
            let l = channel.l;
            let ms = 2 * l + 1;
            // GK14: interpolate from the run table when one is given.
            let table = tables.and_then(|t| t.get(si)).and_then(|t| t.get(ci));
            s.radial_by_bits.clear();
            s.radial_y.clear();
            s.radial_y.resize(ms * n, 0.);
            for (r, &g) in rows.iter().enumerate() {
                let q = &momenta[g];
                let tabulated = table.and_then(|t| t.eval(q.norm_squared.sqrt()));
                let radial = match tabulated {
                    Some(v) => v / volume.sqrt(),
                    None => match s.radial_by_bits.entry(q.norm_squared.to_bits()) {
                        Entry::Occupied(e) => *e.get(),
                        Entry::Vacant(e) => *e.insert(
                            projector_transform_with_rule(
                                sp.potential.mesh(),
                                l,
                                &channel.radial_times_r,
                                q.norm_squared.sqrt(),
                                rule,
                            )
                            .map_err(err)?
                                / volume.sqrt(),
                        ),
                    },
                };
                for (mi, m) in (-(l as isize)..=l as isize).enumerate() {
                    let y = if q.norm_squared == 0. {
                        if l == 0 { 1. / (4. * PI).sqrt() } else { 0. }
                    } else {
                        real_spherical_harmonic(l, m, q.cartesian).map_err(err)?
                    };
                    s.radial_y[mi * n + r] = radial * y;
                }
            }
            let prefactor = C::new(0., -1.).powu(l as u32);
            for pos in sp.positions_bohr {
                s.phase.clear();
                s.phase.extend(rows.iter().map(|&g| {
                    let q = &momenta[g];
                    let angle: f64 = q.cartesian.iter().zip(pos).map(|(q, r)| q * r).sum();
                    prefactor * C::from_polar(1., -angle)
                }));
                for mi in 0..ms {
                    let ry = &s.radial_y[mi * n..(mi + 1) * n];
                    columns.push(s.phase.iter().zip(ry).map(|(&p, &v)| p * v).collect());
                }
            }
        }
    }
    Ok(())
}

/// GK14 run tables: per species (in `species` order), one table per projector
/// (`nl`) and one per atomic-orbital channel (`start`); empty for species
/// without atoms (their rows keep the quadrature). Nodes are split round-robin
/// over the ranks and combined by one sum; each node has exactly one nonzero
/// contribution, so the node values equal the single-rank quadrature bits for
/// any rank count. A quadrature failure on any rank fails every rank.
type KRadialTables = (
    Vec<Vec<RadialTransformTable>>,
    Vec<Vec<RadialTransformTable>>,
);
fn k_radial_tables<K: Collective>(
    species: &[Species<'_>],
    q_max: f64,
    rule: crate::RadialEvenRule,
    comm: &K,
    nl: bool,
    start: bool,
) -> Result<KRadialTables, ScfError> {
    use dft_upf::radial::{RADIAL_TABLE_POINTS as POINTS, RADIAL_TABLE_STEP as STEP};
    let nodes = RadialTransformTable::nodes(q_max, STEP, POINTS).map_err(err)?;
    // (species, projector table?, l, samples)
    let mut channels: Vec<(usize, bool, usize, &[f64])> = Vec::new();
    for (s, sp) in species.iter().enumerate() {
        if sp.positions_bohr.is_empty() {
            continue;
        }
        if nl {
            for p in sp.potential.projectors() {
                channels.push((s, true, p.l(), p.radial_times_r()));
            }
        }
        if start {
            for c in atomic_orbitals(sp.potential) {
                channels.push((s, false, c.l, &c.radial_times_r));
            }
        }
    }
    let total = channels.len() * nodes;
    let mut values = vec![0.; total + 1];
    let (rank, size) = (comm.rank(), comm.size().max(1));
    let mut failed = false;
    for j in (rank..total).step_by(size) {
        let (c, i) = (j / nodes, j % nodes);
        let (s, _, l, f) = channels[c];
        match projector_transform_with_rule(
            species[s].potential.mesh(),
            l,
            f,
            i as f64 * STEP,
            rule,
        ) {
            Ok(v) => values[j] = v,
            Err(_) => failed = true,
        }
    }
    values[total] = f64::from(u8::from(failed));
    let values = sum(comm, &values)?;
    if values[total] > 0. {
        return Err(ScfError::Numerical(
            "reference k start: radial table quadrature failed".into(),
        ));
    }
    let mut out_nl: Vec<Vec<RadialTransformTable>> = species.iter().map(|_| Vec::new()).collect();
    let mut out_start: Vec<Vec<RadialTransformTable>> =
        species.iter().map(|_| Vec::new()).collect();
    for (c, &(s, projector, l, _)) in channels.iter().enumerate() {
        let table = RadialTransformTable::from_nodes(
            l,
            q_max,
            STEP,
            POINTS,
            values[c * nodes..(c + 1) * nodes].to_vec(),
        )
        .map_err(err)?;
        if projector {
            out_nl[s].push(table);
        } else {
            out_start[s].push(table);
        }
    }
    Ok((out_nl, out_start))
}

// Shared atomic-orbital Ritz construction; SCF and fixed-density starts use
// exactly the same columns, completion order and distributed operators.
#[allow(clippy::too_many_arguments)]
fn finish_kpoints<K: Collective>(
    input: StartInput<'_>,
    points: &[dft_geometry::kpoints::KPoint],
    comm: &K,
    bases: Vec<PlaneWaveBasis>,
    db: &PlaneWaveBasis,
    mut commons: Vec<StartState>,
    keep_hx: bool,
    pool: Option<crate::refscf::k_pool::KPool>,
) -> Result<Vec<KStartState>, ScfError> {
    use dft_numerics::linalg::{ComplexMatrix as Matrix, GemmWorkspace, MatrixOp, hermitian_eigh};
    let o = input.options;
    let cell = input.cell;
    let channels = commons.len();
    // Channels are built from identical inputs (GK1B share); the first one
    // carries the stage timers and the start H column count.
    let dims = commons[0].grid.dims();
    let bands = commons[0].orbitals.cols();
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
    let workers = WorkerPool::new(o.threads).map_err(err)?;
    let mut blocks: Vec<Vec<KBlock>> = (0..channels).map(|_| Vec::new()).collect();
    let mut orbitals: Vec<Vec<Matrix>> = (0..channels).map(|_| Vec::new()).collect();
    let mut initial_h_columns = 0;
    // Only the fixed-occupation k driver consumes it (`start_kpoints`).
    let keep_start_hx = keep_hx
        && channels == 1
        && dft_numerics::switches::on(dft_numerics::switches::Id::K_START_KEPT_HX);
    let mut start_hx = Vec::new();
    let memo = dft_numerics::switches::on(dft_numerics::switches::Id::K_START_MEMO);
    let nl_rows_only = dft_numerics::switches::on(dft_numerics::switches::Id::K_NL_TABLE_MEMO);
    let mut memo_scratch = KColumnScratch::default();
    let clock = Instant::now();
    // GK14: radial tables once per start (per cell), covering every k sphere.
    // NL tables act on the rows-only projector path, start tables on the memo
    // columns; the off arms of those switches keep the quadrature reference.
    let nl_table =
        nl_rows_only && dft_numerics::switches::on(dft_numerics::switches::Id::RADIAL_TABLE_NL);
    let start_table =
        memo && dft_numerics::switches::on(dft_numerics::switches::Id::RADIAL_TABLE_START);
    let (nl_tables, start_tables) = if nl_table || start_table {
        let stage = Instant::now();
        let q_max = bases
            .iter()
            .flat_map(|b| b.momenta())
            .map(|m| m.norm_squared)
            .fold(0., f64::max)
            .sqrt();
        let (n, s) = k_radial_tables(
            &species,
            q_max,
            o.radial_even_rule,
            comm,
            nl_table,
            start_table,
        )?;
        commons[0].setup.record("k_radial_tables", stage);
        (nl_table.then_some(n), start_table.then_some(s))
    } else {
        (None, None)
    };
    // GK9 pools: the potentials are gathered once per channel (collective,
    // every rank) and each rank builds only its own whole k points with the
    // serial collective; the outcome is agreed before the next collective.
    let owned = pool.as_ref().map(|p| p.local());
    let kc = &crate::refscf::k_pool::KComm::new(comm, pool.is_some());
    let mut full_potentials = pool
        .is_some()
        .then(|| {
            commons
                .iter()
                .map(|c| {
                    crate::refscf::k_pool::FullField::gather(
                        comm,
                        dims,
                        c.fields.first_axis(),
                        &c.potential,
                    )
                })
                .collect::<Result<Vec<_>, _>>()
        })
        .transpose()?;
    let built = (|| -> Result<(), ScfError> {
        for (point, basis) in points
            .iter()
            .zip(bases)
            .enumerate()
            .filter(|(k, _)| owned.as_ref().is_none_or(|r| r.contains(k)))
            .map(|(_, pair)| pair)
        {
            let common = &mut commons[0];
            if bands > basis.vectors().len() {
                return Err(ScfError::InvalidInput(
                    "reference k start: too many bands for a k basis",
                ));
            }
            // Per-k stage split of complex_k_atomic_starts (one clock per stage
            // per k; diagnostics only, no effect on arithmetic).
            let stage = Instant::now();
            let grid = FftGridSpec::new(&basis, db, Some(dims))?;
            let layout = DistributedPairFft::new(dims, grid.g_to_fft(), kc.rank(), kc.size())
                .map_err(err)?;
            let rows = layout.local_positions().to_vec();
            drop(layout);
            common.setup.record("k_start_grid_layout", stage);
            let stage = Instant::now();
            let nl = if nl_rows_only {
                // GK1B: owned-row factors only. Failures are now rank-local, so
                // the ranks agree on the outcome before any further collective.
                let local = match &nl_tables {
                    Some(tables) => NonlocalOperator::new_owned_rows_only_tabulated(
                        &basis,
                        cell.volume_bohr3(),
                        &species,
                        o.radial_even_rule,
                        &rows,
                        tables,
                    ),
                    None => NonlocalOperator::new_owned_rows_only_with_rule(
                        &basis,
                        cell.volume_bohr3(),
                        &species,
                        o.radial_even_rule,
                        &rows,
                    ),
                };
                let failed = sum(kc, &[f64::from(u8::from(local.is_err()))])?[0];
                let nl = local?;
                if failed > 0. {
                    return Err(ScfError::Numerical(
                        "reference k start: projector table failed on another rank".into(),
                    ));
                }
                nl
            } else {
                NonlocalOperator::new_owned_with_rule(
                    &basis,
                    cell.volume_bohr3(),
                    &species,
                    o.radial_even_rule,
                    &rows,
                )?
            };
            commons[0].setup.record("k_start_nonlocal_tables", stage);
            let stage = Instant::now();
            let mut ops = Vec::with_capacity(channels);
            let mut start_potential = Vec::new();
            for (c, channel) in commons.iter().enumerate() {
                let mut op =
                    dft_physics::distributed_complex::DistributedComplexHamiltonian::new_owned(
                        &basis,
                        &grid,
                        &nl,
                        (2 * bands)
                            .max(channel.memory.atomic_columns)
                            .min(basis.vectors().len()),
                        kc.rank(),
                        kc.size(),
                        workers.clone(),
                    )
                    .map_err(err)?;
                let mut potential = Vec::new();
                if let Some(full) = full_potentials.as_mut() {
                    // GK9: the channel potential was gathered once for all k.
                    potential.extend_from_slice(full[c].for_axis(op.first_axis())?);
                } else {
                    let mut planes = crate::grid_dist::PlaneMap::new(
                        dims,
                        comm.rank(),
                        comm.size(),
                        channel.fields.first_axis(),
                        op.first_axis(),
                    );
                    planes.a_to_b(comm, &channel.potential, &mut potential)?;
                }
                op.set_owned_potential(&potential, 1).map_err(err)?;
                ops.push(op);
                start_potential = potential;
            }
            drop(nl);
            let common = &mut commons[0];
            common.setup.record("k_start_operator_potential", stage);
            let stage = Instant::now();
            let mut columns = Vec::new();
            if memo {
                atomic_k_columns_memo(
                    &species,
                    &basis,
                    &rows,
                    cell.volume_bohr3(),
                    o.radial_even_rule,
                    start_tables.as_deref(),
                    &mut memo_scratch,
                    &mut columns,
                )?;
            } else {
                atomic_k_columns_direct(
                    &species,
                    &basis,
                    &rows,
                    cell.volume_bohr3(),
                    o.radial_even_rule,
                    &mut columns,
                )?;
            }
            common.setup.record("k_start_atomic_columns", stage);
            let stage = Instant::now();
            // Two-pass owned-row MGS, then lowest kinetic coordinate complements.
            // Complements are chosen from the actual k+G kinetic order, never G index.
            let mut order: Vec<_> = (0..basis.vectors().len()).collect();
            order.sort_by(|&a, &b| {
                basis.momenta()[a]
                    .norm_squared
                    .total_cmp(&basis.momenta()[b].norm_squared)
                    .then(a.cmp(&b))
            });
            let real_start = point.fractional == [0.; 3]
                && dft_numerics::switches::on(dft_numerics::switches::Id::REFSCF_REAL_START);
            let half = real_start
                .then(|| HalfSphere::new(&basis, &grid))
                .transpose()
                .map_err(err)?;
            let target = columns.len().max(bands).min(basis.vectors().len());
            let mut q = Matrix::zeros(rows.len(), 0).map_err(err)?;
            q.reserve_columns(target).map_err(err)?;
            let rows_ref = &rows;
            let complements = order.into_iter().flat_map(|g| {
                let mirror = half.as_ref().map(|h| h.mirror()[g] as usize);
                let count = match mirror {
                    Some(m) if g > m => 0,
                    Some(m) if g != m => 2,
                    _ => 1,
                };
                (0..count).map(move |imaginary| {
                    rows_ref
                        .iter()
                        .map(|&r| {
                            let Some(m) = mirror.filter(|&m| m != g) else {
                                return C::new(f64::from(r == g), 0.);
                            };
                            let v = if imaginary == 0 {
                                C::new(std::f64::consts::FRAC_1_SQRT_2, 0.)
                            } else {
                                C::new(0., std::f64::consts::FRAC_1_SQRT_2)
                            };
                            if r == g {
                                v
                            } else if r == m {
                                v.conj()
                            } else {
                                C::default()
                            }
                        })
                        .collect::<Vec<_>>()
                })
            });
            for mut col in columns.into_iter().chain(complements) {
                for _ in 0..2 {
                    for j in 0..q.cols() {
                        let z: C = q
                            .column(j)
                            .iter()
                            .zip(&col)
                            .map(|(a, b)| a.conj() * b)
                            .sum();
                        let z = sum(kc, &[z.re, z.im])?;
                        let z = C::new(z[0], z[1]);
                        for (a, b) in col.iter_mut().zip(q.column(j)) {
                            *a -= z * b;
                        }
                    }
                }
                let norm = sum(kc, &[col.iter().map(|z| z.norm_sqr()).sum()])?[0].sqrt();
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
            common.setup.record("k_start_orthonormalize", stage);
            let stage = Instant::now();
            let mut hq = Matrix::zeros(rows.len(), q.cols()).map_err(err)?;
            ops[0].apply_block(kc, &q, &mut hq).map_err(err)?;
            common.setup.record("k_start_h_apply", stage);
            let stage = Instant::now();
            initial_h_columns += q.cols();
            let mut gemm = GemmWorkspace::with_pool(
                rows.len().max(q.cols()),
                q.cols(),
                o.threads.max(1),
                workers.clone(),
            )
            .map_err(err)?;
            let mut h = Matrix::zeros(q.cols(), q.cols()).map_err(err)?;
            gemm.multiply_into(&q, MatrixOp::Adjoint, &hq, MatrixOp::Identity, &mut h)
                .map_err(err)?;
            let local: Vec<_> = h.as_slice().iter().flat_map(|z| [z.re, z.im]).collect();
            let all = sum(kc, &local)?;
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
            if real_start {
                for z in h.as_mut_slice() {
                    if z.im.abs() > 1e-12 {
                        return Err(err("Gamma start Hamiltonian is not real"));
                    }
                    z.im = 0.;
                }
            }
            let (values, mut c) = hermitian_eigh(&h).map_err(err)?;
            if dft_numerics::switches::on(dft_numerics::switches::Id::REFSCF_CANONICAL_CLUSTER) {
                let mut scratch =
                    vec![C::default(); c.rows() * c.cols() + c.cols() * c.cols() + c.cols()];
                dft_solvers::ref_davidson_complex::canonicalize_clusters(
                    &values,
                    &mut c,
                    &mut scratch,
                );
            }
            c.resize_columns(bands).map_err(err)?;
            let mut x = Matrix::zeros(rows.len(), bands).map_err(err)?;
            gemm.multiply_into(&q, MatrixOp::Identity, &c, MatrixOp::Identity, &mut x)
                .map_err(err)?;
            if keep_start_hx {
                // H X = (H Q) C: one GEMM instead of an H application in map 1.
                let mut hx = Matrix::zeros(rows.len(), bands).map_err(err)?;
                gemm.multiply_into(&hq, MatrixOp::Identity, &c, MatrixOp::Identity, &mut hx)
                    .map_err(err)?;
                start_hx.push((x.clone(), hx, std::mem::take(&mut start_potential)));
            }
            common.setup.record("k_start_ritz", stage);
            let mut ops = ops.into_iter();
            let first = ops.next().expect("at least one channel");
            for (c, op) in ops.enumerate() {
                orbitals[c + 1].push(x.clone());
                blocks[c + 1].push(KBlock {
                    point: *point,
                    basis: basis.clone(),
                    operator: op,
                });
            }
            orbitals[0].push(x);
            blocks[0].push(KBlock {
                point: *point,
                basis,
                operator: first,
            });
        }
        Ok(())
    })();
    if pool.is_some() {
        crate::refscf::k_pool::agree(comm, built)?;
    } else {
        built?;
    }
    commons[0].setup.record("complex_k_atomic_starts", clock);
    let mut out = Vec::with_capacity(channels);
    for (c, ((mut common, blocks), orbitals)) in
        commons.into_iter().zip(blocks).zip(orbitals).enumerate()
    {
        common.orbitals = RealMatrix::zeros(0, bands).map_err(err)?;
        out.push(KStartState {
            initial_h_columns: if c == 0 { initial_h_columns } else { 0 },
            start_hx: if c == 0 {
                std::mem::take(&mut start_hx)
            } else {
                Vec::new()
            },
            common,
            blocks,
            orbitals: KOrbitals(orbitals),
            pool: pool.clone(),
        });
    }
    Ok(out)
}

/// Fixed-density initialization, without an atomic density or an SCF start.
/// `density` is one total-density grid or two collinear (up/down) grids. All
/// inputs are replicated and identical on every rank. Returns the k start and
/// the immutable channel potentials on the common owned field planes.
/// Atomic orbitals and deterministic kinetic complements reuse the SCF helper.
pub fn start_kpoints_from_density<K: Collective>(
    input: StartInput<'_>,
    density: &[Array3<f64>],
    points: &[dft_geometry::kpoints::KPoint],
    comm: &K,
) -> Result<(KStartState, Vec<Vec<f64>>), ScfError> {
    let o = input.options;
    let cell = input.cell;
    dft_geometry::kpoints::validate(points).map_err(err)?;
    if comm.size() == 0
        || comm.rank() >= comm.size()
        || input.atoms.is_empty()
        || input.potentials.is_empty()
        || o.strain_reference.is_some()
        || (o.pw92_prefactor == crate::Pw92Prefactor::Original
            && (o.functional != Functional::Pbe || density.len() != 1))
        || !(1..=2).contains(&density.len())
        || (density.len() == 2 && o.functional == Functional::LdaPz)
    {
        return Err(ScfError::InvalidInput("reference NSCF: unsupported input"));
    }
    let dims = [
        density[0].shape()[0],
        density[0].shape()[1],
        density[0].shape()[2],
    ];
    if dims.contains(&0)
        || density
            .iter()
            .any(|n| n.shape() != dims || n.iter().any(|v| !v.is_finite() || *v < -1e-14))
    {
        return Err(ScfError::InvalidInput("reference NSCF: invalid density"));
    }
    let mut positions = vec![Vec::new(); input.potentials.len()];
    for (i, p) in input.potentials.iter().enumerate() {
        if input.potentials[..i].iter().any(|q| p.label == q.label)
            || (!o.allow_functional_override && !o.functional.matches_upf(p.potential.functional()))
        {
            return Err(ScfError::InvalidInput(
                "reference NSCF: duplicate species or incompatible functional",
            ));
        }
    }
    let mut electrons = 0.;
    for atom in input.atoms {
        let i = input
            .potentials
            .iter()
            .position(|p| p.label == atom.species)
            .ok_or(ScfError::InvalidInput("reference NSCF: unknown species"))?;
        if atom.fractional_position.iter().any(|r| !r.is_finite()) {
            return Err(ScfError::NonFinite);
        }
        positions[i].push(
            cell.fractional_to_cartesian_bohr(atom.fractional_position.map(|r| r.rem_euclid(1.))),
        );
        electrons += input.potentials[i].potential.z_valence();
    }
    if !electrons.is_finite()
        || electrons <= 0.
        || electrons > usize::MAX as f64
        || (electrons - electrons.round()).abs() > 1e-10
        || (density.len() == 1 && !(electrons.round() as usize).is_multiple_of(2))
    {
        return Err(ScfError::InvalidInput(
            "reference NSCF: invalid electron count",
        ));
    }
    let volume = cell.volume_bohr3();
    let charges: Vec<_> = density
        .iter()
        .map(|n| n.sum() * volume / n.len() as f64)
        .collect();
    if (charges.iter().sum::<f64>() - electrons).abs() > 1e-8 * electrons {
        return Err(ScfError::InvalidInput(
            "density charge does not match system",
        ));
    }
    // Spin densities need not integrate to integer sectors (e.g. a metal).
    // Round only charge-check-sized noise before taking the occupied ceiling.
    let occupied = if density.len() == 1 {
        electrons as usize / 2
    } else {
        charges
            .iter()
            .map(|q| (q - 1e-8 * electrons).max(0.).ceil() as usize)
            .max()
            .unwrap()
    };
    let minimum = occupied
        .checked_add(4)
        .ok_or(ScfError::InvalidInput("band overflow"))?;
    let bands = o.bands.ok_or(ScfError::InvalidInput(
        "reference NSCF requires explicit bands",
    ))?;
    if bands < minimum {
        return Err(ScfError::InvalidInput("at least four empty bands required"));
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
    let basis = o.basis(cell, [0.; 3], 1.)?;
    let db = o.basis(cell, [0.; 3], 4.)?;
    let grid = FftGridSpec::new(&basis, &db, Some(dims))?;
    let bases = points
        .iter()
        .map(|p| {
            let b = o.basis(cell, p.fractional, 1.)?;
            FftGridSpec::new(&b, &db, Some(dims))?;
            if bands > b.vectors().len() {
                return Err(ScfError::InvalidInput(
                    "reference NSCF: too many bands for a k basis",
                ));
            }
            Ok(b)
        })
        .collect::<Result<Vec<_>, ScfError>>()?;
    let half = HalfSphere::new(&basis, &grid).map_err(err)?;
    let mut fields = StartFields::new(dims, comm)?;
    let atomic = atomic_fields(&fields, &db, volume, &species, o);
    agree(
        comm,
        f64::from(atomic.is_err()),
        "reference NSCF: ionic/core fields failed",
    )?;
    let atomic = atomic?;
    let ionic = fields.inverse(comm, &atomic.ion)?;
    let core_density = atomic
        .core
        .map(|c| fields.inverse_checked(comm, &c, true))
        .transpose()?;
    let pbe_core = if o.functional == Functional::Pbe {
        core_density
            .as_ref()
            .map(|core| {
                let density_g = fields.forward(comm, core)?;
                Ok::<_, ScfError>(PbeCore {
                    total_density_g: vec![C::default(); density_g.len()],
                    total_density: vec![0.; core.len()],
                    density_g,
                })
            })
            .transpose()?
    } else {
        None
    };
    let offsets = fields.plane_offsets();
    let local: Vec<Vec<f64>> = density
        .iter()
        .map(|n| {
            let full: Vec<_> = n.iter().copied().collect();
            offsets.iter().map(|&i| full[i]).collect()
        })
        .collect();
    let total: Vec<_> = (0..offsets.len())
        .map(|i| local.iter().map(|n| n[i]).sum())
        .collect();
    let density_g = fields.forward(comm, &total)?;
    // The common real adapter is metadata for the existing field interface;
    // no Gamma eigenproblem is solved for fixed-density bands.
    let rows = owned_projector_indices(&basis, &grid, comm.rank(), comm.size()).map_err(err)?;
    let nl = NonlocalOperator::new_owned_with_rule(&basis, volume, &[], o.radial_even_rule, &rows)?;
    let operator = DistributedRealHamiltonian::new_owned(
        &basis,
        &grid,
        &nl,
        1,
        comm.rank(),
        comm.size(),
        WorkerPool::new(o.threads).map_err(err)?,
    )
    .map_err(err)?;
    let common = StartState {
        report: StartReport {
            density: DensityStartKind::Supplied,
            orbitals: orbital_kind(
                species
                    .iter()
                    .map(|s| {
                        s.positions_bohr.len()
                            * atomic_orbitals(s.potential)
                                .map(|c| 2 * c.l + 1)
                                .sum::<usize>()
                    })
                    .sum(),
                bands,
            ),
        },
        setup: Default::default(),
        basis,
        grid,
        half,
        operator,
        fields,
        density: total,
        density_g,
        ionic,
        core_density,
        pbe_core,
        pw92_prefactor: o.pw92_prefactor,
        potential: Vec::new(),
        diagonal_static: Vec::new(),
        orbitals: RealMatrix::zeros(0, bands).map_err(err)?,
        ritz_values_hartree: Vec::new(),
        electrons,
        ewald: ewald_owned(
            cell,
            &input
                .atoms
                .iter()
                .map(|a| Ion {
                    charge: input
                        .potentials
                        .iter()
                        .find(|p| p.label == a.species)
                        .unwrap()
                        .potential
                        .z_valence(),
                    fractional_position: a.fractional_position,
                })
                .collect::<Vec<_>>(),
            comm,
        )?
        .0,
        memory: StartMemory {
            atomic_columns: species
                .iter()
                .map(|s| {
                    s.positions_bohr.len()
                        * atomic_orbitals(s.potential)
                            .map(|c| 2 * c.l + 1)
                            .sum::<usize>()
                })
                .sum(),
            bands,
            ..Default::default()
        },
    };
    let mut state = KStartState {
        initial_h_columns: 0,
        start_hx: Vec::new(),
        common,
        blocks: Vec::new(),
        orbitals: KOrbitals(Vec::new()),
        pool: None,
    };
    let potentials = crate::refscf::nscf::potentials(&mut state, comm, cell, o.functional, &local)?;
    state.common.potential = potentials[0].clone();
    let state = finish_kpoints(
        input,
        points,
        comm,
        bases,
        &db,
        vec![state.common],
        false,
        None,
    )?
    .remove(0);
    Ok((state, potentials))
}
