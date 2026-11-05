//! Ionic and Kleinman–Bylander terms, and a reusable Bloch Hamiltonian.
use dft_basis::{FftGridSpec, GVector, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_numerics::pool::WorkerPool;
use dft_numerics::pool::analysis::{Piece, Span};
use dft_numerics::{
    Fft3d, FftPlan3d,
    fused::{FusedOperator, FusedPlan},
    linalg::ComplexMatrix,
};
use dft_solvers::{LinearOperator, SolverError};
use dft_upf::{Pseudopotential, radial::real_spherical_harmonic};
use ndarray::Array3;
use num_complex::Complex64;
use std::cell::RefCell;
use std::sync::Arc;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum HamiltonianError {
    #[error("incompatible basis, grid, or vector dimensions")]
    Shape,
    #[error("nonfinite input or result")]
    NonFinite,
    #[error("volume must be finite and positive")]
    Volume,
    #[error("D couples different angular momenta")]
    AngularCoupling,
    #[error("Hamiltonian scratch is already borrowed")]
    Busy,
    #[error(transparent)]
    Basis(#[from] dft_basis::BasisError),
    #[error(transparent)]
    Fft(#[from] dft_numerics::FftError),
    #[error(transparent)]
    Radial(#[from] dft_upf::radial::RadialError),
}
fn finite(z: &Complex64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn volume(v: f64) -> Result<(), HamiltonianError> {
    if v.is_finite() && v > 0.0 {
        Ok(())
    } else {
        Err(HamiltonianError::Volume)
    }
}
fn phase(g: [f64; 3], position: [f64; 3]) -> Result<Complex64, HamiltonianError> {
    let angle: f64 = g.iter().zip(position).map(|(g, r)| g * r).sum();
    if !angle.is_finite() {
        return Err(HamiltonianError::NonFinite);
    }
    Ok(Complex64::from_polar(1.0, -angle))
}
/// Positions are Cartesian, in Bohr. Multiple entries may use the same potential.
pub struct Species<'a> {
    pub potential: &'a Pseudopotential,
    pub positions_bohr: &'a [[f64; 3]],
}
/// S(G) = sum_I exp(-i G.tau_I).
pub fn structure_factor(
    g: [f64; 3],
    positions_bohr: &[[f64; 3]],
) -> Result<Complex64, HamiltonianError> {
    positions_bohr.iter().map(|&r| phase(g, r)).sum()
}
/// Coefficients in the supplied density-sphere order; zero outside that sphere.
/// The caller supplies the density basis belonging to its wavefunction basis/cell.
pub fn ionic_local(
    density: &PlaneWaveBasis,
    volume_bohr3: f64,
    species: &[Species<'_>],
) -> Result<Vec<Complex64>, HamiltonianError> {
    volume(volume_bohr3)?;
    let mut result = vec![Complex64::default(); density.vectors().len()];
    for s in species {
        let local = s.potential.local_transform(volume_bohr3, None)?;
        for (out, g) in result.iter_mut().zip(density.vectors()) {
            let v = if g.miller == [0; 3] {
                local.alpha_z()?
            } else {
                local.at_g(g.norm_squared.sqrt())?
            };
            *out += v * structure_factor(g.cartesian, s.positions_bohr)?;
        }
    }
    if !result.iter().all(finite) {
        return Err(HamiltonianError::NonFinite);
    }
    Ok(result)
}
/// A single atom's projector columns, ordered by radial index then m=-l..l.
/// D is expanded only over channels with equal l and m.
pub struct AtomProjectors {
    columns: Vec<Vec<Complex64>>,
    coupling: Vec<f64>,
}
impl AtomProjectors {
    pub fn columns(&self) -> &[Vec<Complex64>] {
        &self.columns
    }
    pub fn coupling_hartree(&self) -> &[f64] {
        &self.coupling
    }
}
/// O(N_G * N_projectors) storage, never an N_G by N_G matrix.
pub struct NonlocalOperator {
    basis: Vec<GVector>,
    atoms: Vec<AtomProjectors>,
    max_channels: usize,
}
impl NonlocalOperator {
    /// Build projectors at q=k+G from the basis momenta, including exp(-i q.tau).
    /// Gamma uses exactly the same radial/angular and structure-factor path.
    pub fn new(
        basis: &PlaneWaveBasis,
        volume_bohr3: f64,
        species: &[Species<'_>],
    ) -> Result<Self, HamiltonianError> {
        volume(volume_bohr3)?;
        let mut atoms = Vec::new();
        let mut max_channels = 0;
        for s in species {
            let pp = s.potential;
            let np = pp.projectors().len();
            for i in 0..np {
                for j in 0..np {
                    if pp.projectors()[i].l() != pp.projectors()[j].l()
                        && pp.dij_hartree()[i * np + j] != 0.0
                    {
                        return Err(HamiltonianError::AngularCoupling);
                    }
                }
            }
            let channels: Vec<_> = pp
                .projectors()
                .iter()
                .enumerate()
                .flat_map(|(i, p)| {
                    let l = p.l() as isize;
                    (-l..=l).map(move |m| (i, p.l(), m))
                })
                .collect();
            let nc = channels.len();
            max_channels = max_channels.max(nc);
            let mut table = vec![vec![Complex64::default(); basis.momenta().len()]; nc];
            // Radial transforms are evaluated once per radial projector and G,
            // independent of the number of atoms and magnetic channels.
            for (ig, g) in basis.momenta().iter().enumerate() {
                let radial: Vec<_> = (0..np)
                    .map(|i| pp.projector_transform(i, g.norm_squared.sqrt()))
                    .collect::<Result<_, _>>()?;
                // Direction at zero is arbitrary: all l>0 radial transforms vanish.
                let direction = if g.norm_squared == 0.0 {
                    [0.0, 0.0, 1.0]
                } else {
                    g.cartesian
                };
                for (c, &(i, l, m)) in channels.iter().enumerate() {
                    table[c][ig] = Complex64::new(0.0, -1.0).powu(l as u32)
                        * real_spherical_harmonic(l, m, direction)?
                        * radial[i]
                        / volume_bohr3.sqrt();
                }
            }
            let mut coupling = vec![0.0; nc * nc];
            for (a, &(i, l, m)) in channels.iter().enumerate() {
                for (b, &(j, ll, mm)) in channels.iter().enumerate() {
                    if l == ll && m == mm {
                        coupling[a * nc + b] = pp.dij_hartree()[i * np + j];
                    }
                }
            }
            for &r in s.positions_bohr {
                let phases = basis
                    .momenta()
                    .iter()
                    .map(|g| phase(g.cartesian, r))
                    .collect::<Result<Vec<_>, _>>()?;
                let columns: Vec<Vec<_>> = table
                    .iter()
                    .map(|col| {
                        col.iter()
                            .zip(&phases)
                            .map(|(p, phase)| p * phase)
                            .collect()
                    })
                    .collect();
                if !columns.iter().flatten().all(finite) {
                    return Err(HamiltonianError::NonFinite);
                }
                atoms.push(AtomProjectors {
                    columns,
                    coupling: coupling.clone(),
                });
            }
        }
        Ok(Self {
            basis: basis.momenta().to_vec(),
            atoms,
            max_channels,
        })
    }
    pub fn atoms(&self) -> &[AtomProjectors] {
        &self.atoms
    }
    fn add(&self, x: &[Complex64], y: &mut [Complex64], scratch: &mut [Complex64]) {
        let (projections, weights) = scratch.split_at_mut(self.max_channels);
        for atom in &self.atoms {
            let nc = atom.columns.len();
            if nc == 0 {
                continue;
            }
            for (p, col) in projections.iter_mut().zip(&atom.columns) {
                *p = col.iter().zip(x).map(|(b, c)| b.conj() * c).sum();
            }
            for (i, weight) in weights[..nc].iter_mut().enumerate() {
                *weight = (0..nc)
                    .map(|j| atom.coupling[i * nc + j] * projections[j])
                    .sum();
            }
            // Original loop structure (one pass over y per projector channel).
            // A 256-element blocked traversal (P10I) was measured 2.9x slower in
            // this function on dev08 (0.929 s versus 0.323 s, MKL T=1 si8).
            for (col, &weight) in atom.columns.iter().zip(&weights[..nc]) {
                for (out, b) in y.iter_mut().zip(col) {
                    *out += b * weight;
                }
            }
        }
    }
    // Column-major band slices. Each band's reduction is the scalar path's
    // left fold starting at complex +0; no reassociation or fused multiply-add.
    // Scratch is channel-major, with adjacent independent band accumulators.
    fn add_block<const G: usize>(
        &self,
        x: &[Complex64],
        y: &mut [Complex64],
        rows: usize,
        bands: usize,
        scratch: &mut [Complex64],
    ) {
        debug_assert!(bands <= G);
        // A one-band tail has no column reuse; retain the scalar traversal.
        if bands == 1 {
            self.add(x, y, scratch);
            return;
        }
        let (projections, weights) = scratch.split_at_mut(self.max_channels * G);
        for atom in &self.atoms {
            let nc = atom.columns.len();
            projections[..nc * G].fill(Complex64::default());
            for (c, col) in atom.columns.iter().enumerate() {
                let p = &mut projections[c * G..(c + 1) * G];
                for (i, b) in col.iter().enumerate() {
                    let conjugate = b.conj();
                    for band in 0..bands {
                        p[band] += conjugate * x[band * rows + i];
                    }
                }
            }
            for i in 0..nc {
                for band in 0..bands {
                    weights[i * G + band] = (0..nc)
                        .map(|j| atom.coupling[i * nc + j] * projections[j * G + band])
                        .sum();
                }
            }
            for (c, col) in atom.columns.iter().enumerate() {
                for (i, b) in col.iter().enumerate() {
                    for band in 0..bands {
                        y[band * rows + i] += b * weights[c * G + band];
                    }
                }
            }
        }
    }

    /// Projection storage required by `apply_into`.
    pub fn scratch_len(&self) -> usize {
        2 * self.max_channels
    }
    /// Caller-owned projection scratch; no allocation per orbital.
    pub fn apply_into(
        &self,
        x: &[Complex64],
        y: &mut [Complex64],
        scratch: &mut [Complex64],
    ) -> Result<(), HamiltonianError> {
        if x.len() != self.basis.len() || y.len() != x.len() || scratch.len() != self.scratch_len()
        {
            return Err(HamiltonianError::Shape);
        }
        if !x.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        y.fill(Complex64::default());
        self.add(x, y, scratch);
        if !y.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        Ok(())
    }
    pub fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), HamiltonianError> {
        if x.len() != self.basis.len() || y.len() != x.len() {
            return Err(HamiltonianError::Shape);
        }
        if !x.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        y.fill(Complex64::default());
        self.add(x, y, &mut vec![Complex64::default(); 2 * self.max_channels]);
        if !y.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        Ok(())
    }
}

/// Filter represented Fourier coefficients to the density sphere, then invert.
/// Already aliased samples cannot be recovered by filtering.
pub fn truncate_to_density_sphere(
    values: &Array3<f64>,
    density: &PlaneWaveBasis,
) -> Result<Array3<f64>, HamiltonianError> {
    let dims = [values.shape()[0], values.shape()[1], values.shape()[2]];
    let grid = FftGridSpec::new(density, density, Some(dims))?;
    if !values.iter().all(|x| x.is_finite()) {
        return Err(HamiltonianError::NonFinite);
    }
    let mut fft = FftPlan3d::new(dims)?;
    let mut real = Array3::from_shape_fn(dims, |idx| Complex64::new(values[idx], 0.0));
    let mut reciprocal = Array3::zeros(dims);
    fft.to_reciprocal(&real, &mut reciprocal)?;
    for (z, index) in reciprocal.iter_mut().zip(grid.fft_to_g()) {
        if index.is_none() {
            *z = Complex64::default();
        }
    }
    fft.to_real(&reciprocal, &mut real)?;
    Ok(real.mapv(|z| z.re))
}
// Fixed independently of the thread budget; tails retain the same sum order.
const NONLOCAL_BANDS: usize = 4;

struct Worker {
    scratch: std::sync::Mutex<Scratch>,
}
struct Scratch {
    fft: FusedPlan,
    projections: Vec<Complex64>,
}
/// Immutable physical data with interior scratch reuse; intentionally not Sync.
/// Local grid values are in Ha and must use the internally derived grid dimensions.
pub struct Hamiltonian {
    data: HamiltonianData,
    scratch: RefCell<Vec<Worker>>,
    pool: Arc<WorkerPool>,
}
struct HamiltonianData {
    grid: FftGridSpec,
    kinetic: Vec<f64>,
    potential: Array3<f64>,
    nonlocal: NonlocalOperator,
}
impl Hamiltonian {
    /// Bloch Hamiltonian for the supplied basis: |k+G|^2/2 + V_local + V_NL(k).
    /// FFTs act on the periodic part; local samples are on the common density grid.
    pub fn new(
        cell: &Cell,
        basis: &PlaneWaveBasis,
        potential: Array3<f64>,
        nonlocal: NonlocalOperator,
    ) -> Result<Self, HamiltonianError> {
        if basis.momenta() != nonlocal.basis {
            return Err(HamiltonianError::Shape);
        }
        let dims = [
            potential.shape()[0],
            potential.shape()[1],
            potential.shape()[2],
        ];
        let grid = FftGridSpec::new(basis, &basis.density_basis(cell)?, Some(dims))?;
        if !potential.iter().all(|x| x.is_finite()) {
            return Err(HamiltonianError::NonFinite);
        }
        let scratch = RefCell::new(vec![Worker {
            scratch: std::sync::Mutex::new(Scratch {
                fft: FusedPlan::new(dims, grid.g_to_fft())?,
                projections: vec![Complex64::default(); 2 * nonlocal.max_channels * NONLOCAL_BANDS],
            }),
        }]);
        Ok(Self {
            data: HamiltonianData {
                grid,
                kinetic: basis
                    .momenta()
                    .iter()
                    .map(|g| 0.5 * g.norm_squared)
                    .collect(),
                potential: if potential.is_standard_layout() {
                    potential
                } else {
                    potential.as_standard_layout().into_owned()
                },
                nonlocal,
            },
            scratch,
            pool: WorkerPool::new(1).map_err(|_| HamiltonianError::Busy)?,
        })
    }
    /// Set the intra-k budget once during construction. Zero means one.
    /// Numerical buffers and persistent workers are reused.
    pub fn with_band_threads(self, threads: usize) -> Result<Self, HamiltonianError> {
        let pool = WorkerPool::new(threads).map_err(|_| HamiltonianError::Busy)?;
        self.with_band_pool(threads, pool)
    }
    /// Share the SCF executor; the band budget limits tasks, not pool size.
    pub fn with_band_pool(
        mut self,
        threads: usize,
        executor: Arc<WorkerPool>,
    ) -> Result<Self, HamiltonianError> {
        self.pool = executor;
        let dims = self.data.grid.dims();
        let pool = self.scratch.get_mut();
        for _ in pool.len()..threads.max(1) {
            pool.push(Worker {
                scratch: std::sync::Mutex::new(Scratch {
                    fft: FusedPlan::new(dims, self.data.grid.g_to_fft())?,
                    projections: vec![
                        Complex64::default();
                        2 * self.data.nonlocal.max_channels * NONLOCAL_BANDS
                    ],
                }),
            });
        }
        pool.truncate(threads.max(1));
        Ok(self)
    }
    /// Replace only the effective local samples, retaining FFT and projector plans.
    /// Invalid input leaves the previous potential intact.
    pub fn set_local_potential(&mut self, potential: Array3<f64>) -> Result<(), HamiltonianError> {
        if potential.shape() != self.data.potential.shape() {
            return Err(HamiltonianError::Shape);
        }
        if !potential.iter().all(|v| v.is_finite()) {
            return Err(HamiltonianError::NonFinite);
        }
        self.data.potential = if potential.is_standard_layout() {
            potential
        } else {
            potential.as_standard_layout().into_owned()
        };
        Ok(())
    }
    /// Copy local samples into retained storage, without replacing its allocation.
    pub fn copy_local_potential(
        &mut self,
        potential: &Array3<f64>,
    ) -> Result<(), HamiltonianError> {
        if potential.shape() != self.data.potential.shape() {
            return Err(HamiltonianError::Shape);
        }
        if potential.iter().any(|v| !v.is_finite()) {
            return Err(HamiltonianError::NonFinite);
        }
        self.data.potential.assign(potential);
        Ok(())
    }
    /// Read-only access for term-by-term orbital energy audits.
    pub fn nonlocal(&self) -> &NonlocalOperator {
        &self.data.nonlocal
    }
    pub fn kinetic_diagonal(&self) -> &[f64] {
        &self.data.kinetic
    }
}
impl HamiltonianData {
    fn apply_group(
        &self,
        x: &ComplexMatrix,
        y: &mut [Complex64],
        first: usize,
        end: usize,
        scratch: &mut Scratch,
    ) -> Result<(), SolverError> {
        for j in first..end {
            let out = &mut y[(j - first) * x.rows()..(j - first + 1) * x.rows()];
            self.apply_local_with(x.column(j), out, scratch)
                .map_err(|e| SolverError::Application(format!("band index {j}: {e}")))?;
        }
        self.apply_nonlocal_group(x, y, first, end, scratch)
    }

    fn apply_nonlocal_group(
        &self,
        x: &ComplexMatrix,
        y: &mut [Complex64],
        first: usize,
        end: usize,
        scratch: &mut Scratch,
    ) -> Result<(), SolverError> {
        self.nonlocal.add_block::<NONLOCAL_BANDS>(
            &x.as_slice()[first * x.rows()..end * x.rows()],
            y,
            x.rows(),
            end - first,
            &mut scratch.projections,
        );
        if !y.iter().all(finite) {
            return Err(application(HamiltonianError::NonFinite));
        }
        Ok(())
    }

    fn apply_piece(
        &self,
        x: &ComplexMatrix,
        y: &mut [Complex64],
        first: usize,
        end: usize,
        scratch: &mut Scratch,
        local_only: bool,
    ) -> Result<(), SolverError> {
        if local_only {
            for j in first..end {
                self.apply_local_with(
                    x.column(j),
                    &mut y[(j - first) * x.rows()..(j - first + 1) * x.rows()],
                    scratch,
                )
                .map_err(|e| SolverError::Application(format!("band index {j}: {e}")))?;
            }
            Ok(())
        } else {
            self.apply_nonlocal_group(x, y, first, end, scratch)
        }
    }

    fn apply_with(
        &self,
        x: &[Complex64],
        y: &mut [Complex64],
        scratch: &mut Scratch,
    ) -> Result<(), HamiltonianError> {
        self.apply_local_with(x, y, scratch)?;
        self.nonlocal.add(x, y, &mut scratch.projections);
        if !y.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        Ok(())
    }
    fn apply_local_with(
        &self,
        x: &[Complex64],
        y: &mut [Complex64],
        scratch: &mut Scratch,
    ) -> Result<(), HamiltonianError> {
        if x.len() != self.kinetic.len() || y.len() != x.len() {
            return Err(HamiltonianError::Shape);
        }
        if !x.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        scratch.fft.apply_local(
            x,
            self.potential.as_slice().ok_or(HamiltonianError::Shape)?,
            y,
        )?;
        for (i, value) in y.iter_mut().enumerate() {
            *value += self.kinetic[i] * x[i];
        }
        Ok(())
    }
}
fn application(error: HamiltonianError) -> SolverError {
    SolverError::Application(error.to_string())
}
impl LinearOperator for Hamiltonian {
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn inner_threads(&self) -> usize {
        self.scratch.borrow().len()
    }
    fn dimension(&self) -> usize {
        self.data.kinetic.len()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        let mut scratch = self
            .scratch
            .try_borrow_mut()
            .map_err(|_| application(HamiltonianError::Busy))?;
        self.data
            .apply_with(
                x,
                y,
                scratch[0]
                    .scratch
                    .get_mut()
                    .map_err(|_| application(HamiltonianError::Busy))?,
            )
            .map_err(application)
    }
    fn apply_block(&self, x: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
        if x.rows() != self.dimension() || y.rows() != x.rows() || y.cols() != x.cols() {
            return Err(application(HamiltonianError::Shape));
        }
        let mut scratch = self
            .scratch
            .try_borrow_mut()
            .map_err(|_| application(HamiltonianError::Busy))?;
        let _region = Span::region("hamiltonian", scratch.len().min(x.cols()).max(1));
        if x.cols() == 0 {
            return Ok(());
        }
        if scratch.len() == 1 || x.cols() <= 1 {
            let scratch = scratch[0]
                .scratch
                .get_mut()
                .map_err(|_| application(HamiltonianError::Busy))?;
            for first in (0..x.cols()).step_by(NONLOCAL_BANDS) {
                let end = (first + NONLOCAL_BANDS).min(x.cols());
                let _piece = Piece::new();
                self.data.apply_group(
                    x,
                    &mut y.as_mut_slice()[first * x.rows()..end * x.rows()],
                    first,
                    end,
                    scratch,
                )?;
            }
            return Ok(());
        }
        // Preserve the FFT's individually stealable bands. Grouping expensive
        // local work with nonlocal work reduced parallelism in the Mac trial.
        // The join between phases completes local+kinetic before any NL update.
        for local_only in [true, false] {
            let width = if local_only { 1 } else { NONLOCAL_BANDS };
            let output = y.as_mut_slice();
            if !local_only && x.cols() <= NONLOCAL_BANDS {
                self.data.apply_nonlocal_group(
                    x,
                    output,
                    0,
                    x.cols(),
                    scratch[0]
                        .scratch
                        .get_mut()
                        .map_err(|_| application(HamiltonianError::Busy))?,
                )?;
            } else if scratch.len() < self.pool.workers_created() || x.cols() <= scratch.len() {
                // Outer-k parallelism uses disjoint static scratch, as before.
                let count = scratch.len().min(x.cols().div_ceil(width));
                let data = &self.data;
                let workers = &mut scratch[..count];
                self.pool.install(|| {
                    apply_static_workers(data, workers, x, output, 0, x.cols(), local_only)
                })?;
            } else {
                let workers = &scratch[..];
                let data = &self.data;
                self.pool
                    .install(|| apply_workers(data, workers, x, output, 0, x.cols(), local_only))?;
            }
        }
        Ok(())
    }
}

// Disjoint column-major output slices need neither a mutex nor a copy.
// Joins precede error propagation in ascending band order.
fn apply_workers(
    data: &HamiltonianData,
    workers: &[Worker],
    x: &ComplexMatrix,
    y: &mut [Complex64],
    first: usize,
    end: usize,
    local_only: bool,
) -> Result<(), SolverError> {
    let width = if local_only { 1 } else { NONLOCAL_BANDS };
    // One local band or fixed nonlocal band group per piece. Scratch belongs to a worker slot and is
    // allocated at Hamiltonian construction, never per band. No nested work
    // runs while holding a slot. Modulo bounds concurrency under outer-k work.
    if end - first <= width {
        let slot = WorkerPool::current_worker_index().unwrap_or(0) % workers.len();
        let mut scratch = workers[slot]
            .scratch
            .lock()
            .map_err(|_| application(HamiltonianError::Busy))?;
        let _piece = Piece::new();
        return data
            .apply_piece(x, y, first, end, &mut scratch, local_only)
            .map_err(|e| SolverError::Application(format!("band index {first}: {e}")));
    }
    let middle = first + (end - first).div_ceil(width) / 2 * width;
    let (left_output, right_output) = y.split_at_mut((middle - first) * x.rows());
    let (left, right) = WorkerPool::join(
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                apply_workers(data, workers, x, left_output, first, middle, local_only)
            }))
        },
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                apply_workers(data, workers, x, right_output, middle, end, local_only)
            }))
        },
    );
    left.map_err(|_| SolverError::Application(format!("band index {first}: worker panicked")))??;
    right.map_err(|_| SolverError::Application(format!("band index {middle}: worker panicked")))?
}

fn apply_static_workers(
    data: &HamiltonianData,
    workers: &mut [Worker],
    x: &ComplexMatrix,
    y: &mut [Complex64],
    first: usize,
    end: usize,
    local_only: bool,
) -> Result<(), SolverError> {
    let width = if local_only { 1 } else { NONLOCAL_BANDS };
    if workers.len() == 1 {
        let _piece = Piece::new();
        let scratch = workers[0]
            .scratch
            .get_mut()
            .map_err(|_| application(HamiltonianError::Busy))?;
        for begin in (first..end).step_by(width) {
            let stop = (begin + width).min(end);
            data.apply_piece(
                x,
                &mut y[(begin - first) * x.rows()..(stop - first) * x.rows()],
                begin,
                stop,
                scratch,
                local_only,
            )?;
        }
        return Ok(());
    }
    let left_count = workers.len() / 2;
    let middle = first + (end - first).div_ceil(width) * left_count / workers.len() * width;
    let (left, right) = workers.split_at_mut(left_count);
    let (left_output, right_output) = y.split_at_mut((middle - first) * x.rows());
    let (left, right) = WorkerPool::join(
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                apply_static_workers(data, left, x, left_output, first, middle, local_only)
            }))
        },
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                apply_static_workers(data, right, x, right_output, middle, end, local_only)
            }))
        },
    );
    left.map_err(|_| SolverError::Application(format!("band index {first}: worker panicked")))??;
    right.map_err(|_| SolverError::Application(format!("band index {middle}: worker panicked")))?
}

#[cfg(test)]
mod pool_tests {
    use super::*;
    #[test]
    fn scratch_addresses_survive_blocks_and_errors() {
        let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
        let basis = PlaneWaveBasis::new(&cell, dft_units::Hartree(2.5)).unwrap();
        let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
        let h = Hamiltonian::new(
            &cell,
            &basis,
            Array3::zeros(grid.dims()),
            NonlocalOperator::new(&basis, cell.volume_bohr3(), &[]).unwrap(),
        )
        .unwrap()
        .with_band_threads(3)
        .unwrap();
        let addresses = || {
            h.scratch
                .borrow()
                .iter()
                .map(|w| {
                    let scratch = w.scratch.lock().unwrap();
                    (
                        &scratch.fft as *const FusedPlan as usize,
                        scratch.projections.as_ptr() as usize,
                    )
                })
                .collect::<Vec<_>>()
        };
        let before = addresses();
        let mut x = ComplexMatrix::zeros(basis.vectors().len(), 5).unwrap();
        let mut y = x.clone();
        for _ in 0..3 {
            h.apply_block(&x, &mut y).unwrap();
            assert_eq!(before, addresses());
        }
        x.column_mut(1)[0].re = f64::NAN;
        assert!(h.apply_block(&x, &mut y).is_err());
        assert_eq!(before, addresses());
    }
}

#[cfg(test)]
mod fused_tests {
    use super::*;
    use dft_numerics::SpherePattern;
    use std::{hint::black_box, time::Instant};
    struct Legacy {
        fft: FftPlan3d,
        pattern: SpherePattern,
        reciprocal: Array3<Complex64>,
        real: Array3<Complex64>,
        projections: Vec<Complex64>,
    }
    impl Legacy {
        fn new(h: &Hamiltonian) -> Self {
            let dims = h.data.grid.dims();
            Self {
                fft: FftPlan3d::new(dims).unwrap(),
                pattern: SpherePattern::from_indices(dims, h.data.grid.g_to_fft()).unwrap(),
                reciprocal: Array3::zeros(dims),
                real: Array3::zeros(dims),
                projections: vec![Complex64::default(); 2 * h.data.nonlocal.max_channels],
            }
        }
        fn apply(&mut self, h: &Hamiltonian, x: &[Complex64], y: &mut [Complex64]) {
            assert!(x.iter().all(finite));
            self.reciprocal.fill(Complex64::default());
            for (&i, &z) in h.data.grid.g_to_fft().iter().zip(x) {
                self.reciprocal[i] = z;
            }
            self.fft
                .to_real_pruned(&self.reciprocal, &self.pattern, &mut self.real)
                .unwrap();
            for (z, &v) in self.real.iter_mut().zip(&h.data.potential) {
                *z *= v;
            }
            self.fft
                .to_reciprocal_pruned(&self.real, &self.pattern, &mut self.reciprocal)
                .unwrap();
            for (i, &idx) in h.data.grid.g_to_fft().iter().enumerate() {
                y[i] = self.reciprocal[idx] + h.data.kinetic[i] * x[i];
            }
            h.data.nonlocal.add(x, y, &mut self.projections);
            assert!(y.iter().all(finite));
        }
    }
    fn vector(n: usize, seed: usize) -> Vec<Complex64> {
        (0..n)
            .map(|i| {
                Complex64::new(
                    ((i * 17 + seed) as f64).sin(),
                    ((i * 31 + seed) as f64).cos(),
                )
            })
            .collect()
    }
    fn bits(v: &[Complex64]) -> Vec<(u64, u64)> {
        v.iter().map(|z| (z.re.to_bits(), z.im.to_bits())).collect()
    }
    fn synthetic(shape: [usize; 3], k: [f64; 3], channels: usize) -> Hamiltonian {
        let cell =
            Cell::from_user_rows_bohr([[5., 0., 0.], [0.7, 6., 0.], [0.2, 0.3, 7.]]).unwrap();
        let basis = PlaneWaveBasis::at_k(&cell, dft_units::Hartree(2.), k).unwrap();
        let nonlocal = NonlocalOperator {
            basis: basis.momenta().to_vec(),
            max_channels: channels,
            atoms: vec![AtomProjectors {
                columns: (0..channels)
                    .map(|i| vector(basis.vectors().len(), i))
                    .collect(),
                coupling: (0..channels * channels)
                    .map(|i| {
                        if i / channels == i % channels {
                            0.03
                        } else {
                            0.01
                        }
                    })
                    .collect(),
            }],
        };
        let potential = Array3::from_shape_vec(
            shape,
            vector(shape.iter().product(), 4)
                .iter()
                .map(|z| z.re)
                .collect(),
        )
        .unwrap();
        Hamiltonian::new(&cell, &basis, potential, nonlocal).unwrap()
    }
    // Independent pre-P10J traversal oracle, retained for bitwise comparison.
    fn legacy_nonlocal(
        op: &NonlocalOperator,
        x: &[Complex64],
        y: &mut [Complex64],
        p: &mut [Complex64],
    ) {
        for atom in &op.atoms {
            let nc = atom.columns.len();
            for (p, col) in p.iter_mut().zip(&atom.columns) {
                *p = col.iter().zip(x).map(|(b, c)| b.conj() * c).sum();
            }
            for (i, col) in atom.columns.iter().enumerate() {
                let weight: Complex64 = (0..nc).map(|j| atom.coupling[i * nc + j] * p[j]).sum();
                for (out, b) in y.iter_mut().zip(col) {
                    *out += b * weight;
                }
            }
        }
    }

    fn nonlocal_fixture(atoms: usize, rows: usize, channels: usize) -> NonlocalOperator {
        NonlocalOperator {
            basis: Vec::new(),
            max_channels: channels,
            atoms: (0..atoms)
                .map(|a| AtomProjectors {
                    columns: (0..channels)
                        .map(|c| vector(rows, a * channels + c))
                        .collect(),
                    coupling: (0..channels * channels)
                        .map(|i| if i % (channels + 1) == 0 { 0.7 } else { 0.03 })
                        .collect(),
                })
                .collect(),
        }
    }

    fn blocked_bits<const G: usize>(op: &NonlocalOperator, rows: usize) {
        for bands in 1..=G {
            let x = vector(rows * bands, 71);
            let mut expected = vector(rows * bands, 31);
            let mut actual = expected.clone();
            let mut scratch = vec![Complex64::default(); op.scratch_len() * G];
            for b in 0..bands {
                op.add(
                    &x[b * rows..(b + 1) * rows],
                    &mut expected[b * rows..(b + 1) * rows],
                    &mut scratch,
                );
            }
            op.add_block::<G>(&x, &mut actual, rows, bands, &mut scratch);
            assert_eq!(bits(&expected), bits(&actual), "G={G}, bands={bands}");
        }
    }

    #[test]
    fn band_block_preserves_every_bit_and_tail() {
        for (atoms, rows, channels) in [(2, 67, 4), (8, 129, 8), (20, 73, 12), (3, 71, 0)] {
            let op = nonlocal_fixture(atoms, rows, channels);
            blocked_bits::<1>(&op, rows);
            blocked_bits::<2>(&op, rows);
            blocked_bits::<4>(&op, rows);
            blocked_bits::<8>(&op, rows);
            blocked_bits::<16>(&op, rows);
        }
    }

    fn block_measure<const G: usize>(op: &NonlocalOperator, rows: usize) {
        let bands: usize = 32;
        let x = vector(rows * bands, 71);
        let mut y = vec![Complex64::default(); rows * bands];
        let mut scratch = vec![Complex64::default(); op.scratch_len() * G.max(1)];
        let mut samples = [0.0; 31];
        for sample in &mut samples {
            let start = Instant::now();
            for _ in 0..4 {
                y.fill(Complex64::default());
                for first in (0..bands).step_by(G.max(1)) {
                    let end = (first + G.max(1)).min(bands);
                    if G == 0 {
                        op.add(
                            &x[first * rows..end * rows],
                            &mut y[first * rows..end * rows],
                            &mut scratch,
                        );
                    } else {
                        op.add_block::<G>(
                            black_box(&x[first * rows..end * rows]),
                            &mut y[first * rows..end * rows],
                            rows,
                            end - first,
                            &mut scratch,
                        );
                    }
                }
                black_box(&y);
            }
            *sample = start.elapsed().as_secs_f64() * 1e6 / (4 * bands) as f64;
        }
        samples.sort_by(f64::total_cmp);
        let footprint: usize = op
            .atoms
            .iter()
            .flat_map(|a| &a.columns)
            .map(|c| c.len() * 16)
            .sum();
        println!(
            "{},{rows},{},{G},{:.6},{footprint},{}",
            op.atoms.len(),
            op.max_channels,
            samples[15],
            2 * footprint * bands.div_ceil(G.max(1)) / bands
        );
    }

    #[test]
    #[ignore] // cargo test --release -p dft-physics band_block_microbenchmark -- --ignored --nocapture
    fn band_block_microbenchmark() {
        println!(
            "atoms,rows,channels,g,us_per_vector,projector_bytes,projector_read_bytes_per_vector"
        );
        for (atoms, rows, channels) in [(2, 749, 4), (8, 2969, 4), (20, 8193, 12)] {
            let op = nonlocal_fixture(atoms, rows, channels);
            block_measure::<0>(&op, rows);
            block_measure::<1>(&op, rows);
            block_measure::<2>(&op, rows);
            block_measure::<4>(&op, rows);
            block_measure::<8>(&op, rows);
            block_measure::<16>(&op, rows);
        }
    }

    #[test]
    fn grouped_hamiltonian_matches_scalar_across_workers() {
        for channels in [0, 4, 8, 12] {
            let h = synthetic([12, 12, 12], [0.13, 0.07, 0.11], channels);
            let rows = h.dimension();
            let x = ComplexMatrix::new(rows, 17, vector(rows * 17, 91)).unwrap();
            let mut expected = x.clone();
            for b in 0..x.cols() {
                h.apply(x.column(b), expected.column_mut(b)).unwrap();
            }
            for threads in [1, 2, 3, 8] {
                let h = synthetic([12, 12, 12], [0.13, 0.07, 0.11], channels)
                    .with_band_threads(threads)
                    .unwrap();
                let mut actual = x.clone();
                h.apply_block(&x, &mut actual).unwrap();
                assert_eq!(bits(expected.as_slice()), bits(actual.as_slice()));
            }
            // Also exercise a bounded inner budget on a larger shared executor.
            let h = synthetic([12, 12, 12], [0.13, 0.07, 0.11], channels)
                .with_band_pool(3, WorkerPool::new(8).unwrap())
                .unwrap();
            let mut actual = x.clone();
            h.apply_block(&x, &mut actual).unwrap();
            assert_eq!(bits(expected.as_slice()), bits(actual.as_slice()));
        }
    }

    #[test]
    fn nonlocal_traversal_preserves_bits() {
        for (atoms, rows, channels) in [
            (4, 257, 1),
            (8, 513, 3),
            (8, 1537, 4),
            (20, 769, 9),
            (3, 259, 0),
        ] {
            let op = nonlocal_fixture(atoms, rows, channels);
            let x = vector(rows, 71);
            let mut a = vector(rows, 31);
            let mut b = a.clone();
            let mut scratch = vec![Complex64::default(); 2 * channels];
            legacy_nonlocal(&op, &x, &mut a, &mut scratch);
            op.add(&x, &mut b, &mut scratch);
            assert_eq!(bits(&a), bits(&b));
        }
    }

    #[test]
    #[ignore] // cargo test --release -p dft-physics nonlocal_traversal -- --ignored --nocapture
    fn nonlocal_traversal() {
        println!("atoms,rows,channels,old_us,new_us");
        for (atoms, rows, channels) in [
            (4, 257, 1),
            (2, 769, 4),
            (4, 1025, 4),
            (4, 4097, 4),
            (8, 2969, 4),
            (8, 4097, 3),
            (20, 16385, 9),
        ] {
            let op = nonlocal_fixture(atoms, rows, channels);
            let x = vector(rows, 71);
            let mut a = vector(rows, 31);
            let mut b = a.clone();
            let mut scratch = vec![Complex64::default(); 2 * channels];
            legacy_nonlocal(&op, &x, &mut a, &mut scratch);
            op.add(&x, &mut b, &mut scratch);
            assert_eq!(bits(&a), bits(&b));
            let mut times = [0.0; 2];
            for (mode, time) in times.iter_mut().enumerate() {
                let mut samples = Vec::new();
                for _ in 0..5 {
                    let start = Instant::now();
                    for _ in 0..50 {
                        b.fill(Complex64::default());
                        if mode == 0 {
                            legacy_nonlocal(black_box(&op), &x, &mut b, &mut scratch);
                        } else {
                            black_box(&op).add(&x, &mut b, &mut scratch);
                        }
                        black_box(&b);
                    }
                    samples.push(start.elapsed().as_secs_f64() * 1e6 / 50.0);
                }
                samples.sort_by(f64::total_cmp);
                *time = samples[2];
            }
            println!("{atoms},{rows},{channels},{:.3},{:.3}", times[0], times[1]);
        }
    }

    #[test]
    fn ordered_complete_hamiltonian_shapes_and_projectors() {
        for shape in [
            [24, 24, 24],
            [25, 25, 25],
            [27, 27, 27],
            [36, 36, 36],
            [15, 20, 24],
            [15, 17, 19],
        ] {
            for (k, channels) in [([0.; 3], 0), ([0.25, -0.25, 0.125], 3), ([0.125; 3], 5)] {
                let h = synthetic(shape, k, channels);
                let mut old = Legacy::new(&h);
                let x = vector(h.dimension(), 12);
                let mut expected = x.clone();
                let mut actual = x.clone();
                old.apply(&h, &x, &mut expected);
                h.apply(&x, &mut actual).unwrap();
                assert_eq!(bits(&actual), bits(&expected));
            }
        }
    }
    #[test]
    #[ignore = "P10C_SI_UPF=/external/Si.upf cargo test --release -p dft-physics complete_hamiltonian_timing -- --ignored --nocapture"]
    fn complete_hamiltonian_timing() {
        let pp =
            dft_upf::from_file(std::env::var("P10C_SI_UPF").expect("external Si UPF required"))
                .unwrap();
        let a = 10.261212856717933;
        let primitive = Cell::from_user_rows_bohr([
            [0., a / 2., a / 2.],
            [a / 2., 0., a / 2.],
            [a / 2., a / 2., 0.],
        ])
        .unwrap();
        let conventional =
            Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap();
        let mut cases = Vec::new();
        for (label, cell, sites, k) in [
            (
                "Si-Gamma",
                primitive.clone(),
                vec![[0.; 3], [0.25; 3]],
                [0.; 3],
            ),
            ("Si-shifted", primitive, vec![[0.; 3], [0.25; 3]], [0.25; 3]),
            (
                "Si8",
                conventional,
                vec![
                    [0.; 3],
                    [0.25; 3],
                    [0., 0.5, 0.5],
                    [0.5, 0., 0.5],
                    [0.5, 0.5, 0.],
                    [0.25, 0.75, 0.75],
                    [0.75, 0.25, 0.75],
                    [0.75, 0.75, 0.25],
                ],
                [0.; 3],
            ),
        ] {
            let basis = PlaneWaveBasis::at_k(&cell, dft_units::Hartree(15.), k).unwrap();
            let grid =
                FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
            let positions: Vec<_> = sites
                .iter()
                .map(|&r| cell.fractional_to_cartesian_bohr(r))
                .collect();
            let nonlocal = NonlocalOperator::new(
                &basis,
                cell.volume_bohr3(),
                &[Species {
                    potential: &pp,
                    positions_bohr: &positions,
                }],
            )
            .unwrap();
            let potential = Array3::from_shape_vec(
                grid.dims(),
                vector(grid.dims().iter().product(), 4)
                    .iter()
                    .map(|z| z.re)
                    .collect(),
            )
            .unwrap();
            cases.push((
                label,
                Hamiltonian::new(&cell, &basis, potential, nonlocal).unwrap(),
            ));
        }
        cases.push((
            "synthetic-skew-3-projectors",
            synthetic([15, 20, 24], [0.25, -0.25, 0.125], 3),
        ));
        for (label, h) in cases {
            let mut old = Legacy::new(&h);
            let x = vector(h.dimension(), 17);
            let mut out = x.clone();
            let mut expected = x.clone();
            old.apply(&h, &x, &mut expected);
            h.apply(&x, &mut out).unwrap();
            assert_eq!(bits(&expected), bits(&out));
            let mut before = Vec::new();
            let mut after = Vec::new();
            for _ in 0..5 {
                let t = Instant::now();
                for _ in 0..128 {
                    old.apply(&h, black_box(&x), &mut out);
                    black_box(&out);
                }
                before.push(t.elapsed().as_secs_f64() * 1e6 / 128.);
                let t = Instant::now();
                for _ in 0..128 {
                    h.apply(black_box(&x), &mut out).unwrap();
                    black_box(&out);
                }
                after.push(t.elapsed().as_secs_f64() * 1e6 / 128.);
            }
            println!(
                "{label} shape={:?} npw={} before_us={before:?} after_us={after:?}",
                h.data.grid.dims(),
                h.dimension()
            );
        }
    }
}
