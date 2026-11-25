//! Isolated real Gamma Hamiltonian. Derived from real-gamma-design sections 4–7;
//! no engine caller or change to the complex traversal (ledger A43).
use crate::distributed_real_gamma::lane_profile::{
    Profile as NonlocalProfile, Stage as NonlocalStage,
};
use crate::hamiltonian::{Hamiltonian, NonlocalOperator};
use dft_basis::{HalfSphere, PlaneWaveBasis};
use dft_numerics::fused_pair_plan::RealPairPlan;
use dft_numerics::work_timers::{self as timers, Slot};
use dft_numerics::{
    PairMap, RealGemmWorkspace, RealMatrix, RealOp, fused_pair::PairLocalWorkspace,
    pool::WorkerPool,
};
use num_complex::Complex64;
use std::sync::Arc;
use thiserror::Error;

/// Construction or application failure; callers may later dispatch to complex H.
#[derive(Debug, Error)]
pub enum RealGammaError {
    #[error("incompatible basis, block dimensions, or capacity")]
    Shape,
    #[error("nonfinite input or result")]
    NonFinite,
    #[error("projector G=0 imaginary component exceeds 1e-14 times column maximum")]
    ComplexZero,
    #[error(transparent)]
    Half(#[from] dft_basis::HalfSphereError),
    #[error(transparent)]
    Fft(#[from] dft_numerics::FftError),
    #[error(transparent)]
    Algebra(#[from] dft_numerics::linalg::LinalgError),
}

pub(crate) fn pack_projector(
    half: &HalfSphere,
    beta: &[Complex64],
    out: &mut [f64],
) -> Result<(), RealGammaError> {
    if beta.len() != half.full_len() || out.len() != half.full_len() {
        return Err(RealGammaError::Shape);
    }
    // Only representatives participate in conversion or the scale check.
    let mut scale = 0.0_f64;
    for &p in half.plus() {
        let z = beta[p as usize];
        if !z.re.is_finite() || !z.im.is_finite() || !z.norm().is_finite() {
            return Err(RealGammaError::NonFinite);
        }
        scale = scale.max(z.norm());
    }
    // Zero columns are valid. Equality is allowed for their exact zero origin.
    if beta[0].im.abs() > 1e-14 * scale {
        return Err(RealGammaError::ComplexZero);
    }
    out[0] = beta[0].re;
    for (i, &p) in half.plus().iter().enumerate().skip(1) {
        out[2 * i - 1] = std::f64::consts::SQRT_2 * beta[p as usize].re;
        out[2 * i] = std::f64::consts::SQRT_2 * beta[p as usize].im;
    }
    if out.iter().any(|x| !x.is_finite()) {
        return Err(RealGammaError::NonFinite);
    }
    Ok(())
}

/// ISSUE-20: the rows `rows` (ORC row indices) of `pack_projector(half, beta)`
/// written directly to `out` (one value per row, same coefficient mapping and
/// arithmetic: row 0 = Re beta(G=0), rows 2i-1 and 2i = sqrt(2) Re and Im of
/// beta at `plus[i]`), without the full-length vector. The validation is that
/// of `pack_projector` over ALL representatives, so the verdict (and the
/// error) does not depend on which rows a rank owns. The final finiteness
/// check is equivalent: sqrt(2) * x rounds monotonically in |x|, so some
/// scaled component overflows exactly when the largest one does.
pub(crate) fn pack_projector_rows(
    half: &HalfSphere,
    beta: &[Complex64],
    rows: &[usize],
    out: &mut [f64],
) -> Result<(), RealGammaError> {
    if beta.len() != half.full_len() || out.len() != rows.len() {
        return Err(RealGammaError::Shape);
    }
    let plus = half.plus();
    let mut scale = 0.0_f64;
    let mut largest = 0.0_f64;
    for (i, &p) in plus.iter().enumerate() {
        let z = beta[p as usize];
        if !z.re.is_finite() || !z.im.is_finite() || !z.norm().is_finite() {
            return Err(RealGammaError::NonFinite);
        }
        scale = scale.max(z.norm());
        if i > 0 {
            largest = largest.max(z.re.abs()).max(z.im.abs());
        }
    }
    if beta[0].im.abs() > 1e-14 * scale {
        return Err(RealGammaError::ComplexZero);
    }
    if !(std::f64::consts::SQRT_2 * largest).is_finite() {
        return Err(RealGammaError::NonFinite);
    }
    for (o, &r) in out.iter_mut().zip(rows) {
        *o = if r == 0 {
            beta[0].re
        } else {
            let z = beta[*plus.get(r.div_ceil(2)).ok_or(RealGammaError::Shape)? as usize];
            if r % 2 == 1 {
                std::f64::consts::SQRT_2 * z.re
            } else {
                std::f64::consts::SQRT_2 * z.im
            }
        };
    }
    Ok(())
}

struct Coupling {
    start: usize,
    count: usize,
    values: Vec<f64>,
}

/// Per column of `projections` (count rows) the energy sum over atoms, then i,
/// then j of d_ij p_i p_j (Ha), serial and in a fixed order.
fn coupled_energies(
    coupling: &[Coupling],
    projections: &[f64],
    count: usize,
    out: &mut [f64],
) -> Result<(), RealGammaError> {
    for (band, value) in out.iter_mut().enumerate() {
        let p = &projections[band * count..(band + 1) * count];
        let mut energy = 0.0;
        for atom in coupling {
            let nc = atom.count;
            for i in 0..nc {
                for j in 0..nc {
                    let d = atom.values[i * nc + j];
                    if d == 0.0 {
                        continue;
                    }
                    energy += d * (p[atom.start + i] * p[atom.start + j]);
                }
            }
        }
        *value = energy;
    }
    if out.iter().any(|v| !v.is_finite()) {
        return Err(RealGammaError::NonFinite);
    }
    Ok(())
}

/// P D P^T on ORC blocks, with construction-sized scratch and deterministic GEMM.
/// The source must describe the same Gamma basis as `half`; use
/// `RealHamiltonian::from_complex` to validate the complete basis identity.
pub struct RealNonlocal {
    projectors: RealMatrix,
    coupling: Vec<Coupling>,
    projections: RealMatrix,
    weights: RealMatrix,
    image: RealMatrix,
    gemm: RealGemmWorkspace,
    capacity: usize,
    /// Profiling counters for the frozen wide-panel projector products.
    lane: NonlocalProfile,
}
impl RealNonlocal {
    /// Copy the public complex columns and each atom's row-major D block.
    pub fn new(
        source: &NonlocalOperator,
        half: &HalfSphere,
        capacity: usize,
        pool: Arc<WorkerPool>,
    ) -> Result<Self, RealGammaError> {
        let n = half.full_len();
        let count = source.atoms().iter().map(|a| a.columns().len()).sum();
        let mut projectors = RealMatrix::zeros(n, count)?;
        let mut coupling = Vec::new();
        let mut start = 0;
        for atom in source.atoms() {
            let count = atom.columns().len();
            for (j, beta) in atom.columns().iter().enumerate() {
                pack_projector(
                    half,
                    beta,
                    &mut projectors.as_mut_slice()[(start + j) * n..(start + j + 1) * n],
                )?;
            }
            coupling.push(Coupling {
                start,
                count,
                values: atom.coupling_hartree().to_vec(),
            });
            start += count;
        }
        let lane = NonlocalProfile::new(coupling.iter().map(|c| c.count)).with_sparse_audit(|| {
            coupling
                .iter()
                .flat_map(|c| &c.values)
                .filter(|&&v| v != 0.)
                .count()
        });
        Ok(Self {
            projectors,
            coupling,
            projections: RealMatrix::zeros(count, capacity)?,
            weights: RealMatrix::zeros(count, capacity)?,
            image: RealMatrix::zeros(n, capacity)?,
            gemm: RealGemmWorkspace::with_pool(n, count.max(capacity), pool)?,
            capacity,
            lane,
        })
    }
    /// Add the nonlocal image to caller-owned output; no allocation on application.
    pub fn add_block(&mut self, x: &RealMatrix, y: &mut RealMatrix) -> Result<(), RealGammaError> {
        let n = self.projectors.rows();
        let m = x.cols();
        if x.rows() != n || y.rows() != n || y.cols() != m || m > self.capacity {
            return Err(RealGammaError::Shape);
        }
        if x.as_slice()
            .iter()
            .chain(y.as_slice())
            .any(|x| !x.is_finite())
        {
            return Err(RealGammaError::NonFinite);
        }
        self.add_block_core(x, y)
    }
    /// `add_block` without the input scans (APPLY-LEAN: the caller scanned X; a
    /// non-finite Y entry survives the addition and fails the final output scan
    /// in `add_from_projections` with the same error).
    fn add_block_core(&mut self, x: &RealMatrix, y: &mut RealMatrix) -> Result<(), RealGammaError> {
        let n = self.projectors.rows();
        let m = x.cols();
        if x.rows() != n || y.rows() != n || y.cols() != m || m > self.capacity {
            return Err(RealGammaError::Shape);
        }
        let count = self.projectors.cols();
        let clock = timers::start();
        let lane = self.lane.start(NonlocalStage::Projection, n, count, m);
        {
            self.gemm.multiply_slices_wide_into(
                dft_numerics::GemmSite::Other,
                self.projectors.as_slice(),
                n,
                count,
                n,
                RealOp::Transpose,
                x.as_slice(),
                m,
                n,
                self.projections.as_mut().subcols_mut(0, m),
            )?;
        }
        drop(lane);
        timers::stop(clock, Slot::NonlocalProjectionNs);
        let clock = timers::start();
        let result = self.add_from_projections(y, m);
        timers::stop(clock, Slot::NonlocalBackNs);
        result
    }

    /// Number of projector columns (rows of the projection block P^T X).
    pub fn projector_count(&self) -> usize {
        self.projectors.cols()
    }

    /// Projections P^T X of the first `m` columns of the most recent
    /// application, column-major (`projector_count()` rows per column). B6:
    /// a solver may cache these for its basis vectors.
    pub fn last_projections(&self, m: usize) -> Result<&[f64], RealGammaError> {
        if m > self.capacity {
            return Err(RealGammaError::Shape);
        }
        Ok(&self.projections.as_slice()[..self.projectors.cols() * m])
    }

    /// Cached-projector certificate: add the nonlocal image of a block whose
    /// projections P^T X the caller already holds (column-major, `projector_count()`
    /// rows, `y.cols()` columns), skipping the first projector product. For
    /// X = V C with cached P^T V the supplied block is (P^T V) C, which equals
    /// P^T (V C) to rounding (different summation order; not bit-identical).
    pub fn add_block_with_projections(
        &mut self,
        projections: &[f64],
        y: &mut RealMatrix,
    ) -> Result<(), RealGammaError> {
        let n = self.projectors.rows();
        let m = y.cols();
        let count = self.projectors.cols();
        if y.rows() != n || m > self.capacity || projections.len() != count * m {
            return Err(RealGammaError::Shape);
        }
        if projections
            .iter()
            .chain(y.as_slice())
            .any(|x| !x.is_finite())
        {
            return Err(RealGammaError::NonFinite);
        }
        let clock = timers::start();
        self.projections.as_mut_slice()[..count * m].copy_from_slice(projections);
        let result = self.add_from_projections(y, m);
        timers::stop(clock, Slot::NonlocalBackNs);
        result
    }

    /// Real nonlocal energy: per column the
    /// nonlocal expectation x^T P D P^T x in Ha from one projector product
    /// (no back product, no output block). Per column the sum runs over atoms,
    /// then i, then j of d_ij p_i p_j, as the complex `expectation_values`;
    /// equal to it to rounding (ORC packing is an isometry), independent of the
    /// thread count (fixed-order GEMM, serial per-column sums).
    pub fn expectation_values(
        &mut self,
        x: &RealMatrix,
        out: &mut [f64],
    ) -> Result<(), RealGammaError> {
        let n = self.projectors.rows();
        let m = x.cols();
        if x.rows() != n || out.len() != m || m > self.capacity {
            return Err(RealGammaError::Shape);
        }
        let count = self.projectors.cols();
        let lane = self.lane.start(NonlocalStage::Projection, n, count, m);
        {
            self.gemm.multiply_slices_wide_into(
                dft_numerics::GemmSite::Other,
                self.projectors.as_slice(),
                n,
                count,
                n,
                RealOp::Transpose,
                x.as_slice(),
                m,
                n,
                self.projections.as_mut().subcols_mut(0, m),
            )?;
        }
        drop(lane);
        let _lane = self.lane.start(NonlocalStage::Audit, n, count, m);
        coupled_energies(
            &self.coupling,
            &self.projections.as_slice()[..count * m],
            count,
            out,
        )
    }

    /// C3a (`BAANDI_LEAN_PROJ_CACHE=on` callers): the per-column nonlocal
    /// expectation values of `expectation_values` from projections P^T X the
    /// caller already holds (column-major, `projector_count()` rows, one column
    /// per entry of `out`), e.g. (P^T V) C kept by the eigensolver. The coupling
    /// sum is the same code; with the same projections the result is bit for
    /// bit that of `expectation_values`, with (P^T V) C equal to rounding.
    pub fn expectation_values_from_projections(
        &self,
        projections: &[f64],
        out: &mut [f64],
    ) -> Result<(), RealGammaError> {
        let count = self.projectors.cols();
        if projections.len() != count * out.len() {
            return Err(RealGammaError::Shape);
        }
        let _lane = self.lane.start(
            NonlocalStage::Audit,
            self.projectors.rows(),
            count,
            out.len(),
        );
        coupled_energies(&self.coupling, projections, count, out)
    }

    /// Coupling and back product from `self.projections` (first `m` columns).
    fn add_from_projections(&mut self, y: &mut RealMatrix, m: usize) -> Result<(), RealGammaError> {
        self.image_from_projections(m)?;
        self.add_image(y)
    }

    /// Y += image (first `Y.cols()` columns), then the output scan.
    fn add_image(&self, y: &mut RealMatrix) -> Result<(), RealGammaError> {
        let _lane = self.lane.start(
            NonlocalStage::Add,
            self.projectors.rows(),
            self.projectors.cols(),
            y.cols(),
        );
        for (y, v) in y.as_mut_slice().iter_mut().zip(self.image.as_slice()) {
            *y += v;
        }
        if y.as_slice().iter().any(|x| !x.is_finite()) {
            return Err(RealGammaError::NonFinite);
        }
        Ok(())
    }

    /// Coupling and back product into `self.image` (first `m` columns).
    fn image_from_projections(&mut self, m: usize) -> Result<(), RealGammaError> {
        let n = self.projectors.rows();
        let count = self.projectors.cols();
        let lane = self.lane.start(NonlocalStage::Coupling, n, count, m);
        for atom in &self.coupling {
            for band in 0..m {
                for j in 0..atom.count {
                    let mut value = 0.;
                    for k in 0..atom.count {
                        value += atom.values[j * atom.count + k]
                            * self.projections[(atom.start + k, band)];
                    }
                    self.weights[(atom.start + j, band)] = value;
                }
            }
        }
        drop(lane);
        let _lane = self.lane.start(NonlocalStage::Back, n, count, m);
        {
            self.gemm.multiply_slices_wide_into(
                dft_numerics::GemmSite::Other,
                self.projectors.as_slice(),
                n,
                count,
                n,
                RealOp::Identity,
                &self.weights.as_slice()[..count * m],
                m,
                count,
                self.image.as_mut().subcols_mut(0, m),
            )?;
        }
        Ok(())
    }
}

/// Snapshot of a complex Gamma H, acting on column-major ORC blocks.
/// Scratch capacity is fixed at construction; nonzero paired band norms must
/// differ by at most 1e3 (R2 debug assertion). Fixed local pair pieces and GEMM
/// share the supplied persistent pool.
pub struct RealHamiltonian {
    source_identity: Arc<()>,
    half: HalfSphere,
    kinetic: Vec<f64>,
    potential: Vec<f64>,
    local: PairLocalWorkspace<RealPairPlan>,
    nonlocal: RealNonlocal,
    capacity: usize,
    grid_shape: [usize; 3],
}
impl RealHamiltonian {
    /// Snapshot potential and public projectors; reject a mismatched or non-Gamma basis.
    pub fn from_complex(
        source: &Hamiltonian,
        basis: &PlaneWaveBasis,
        capacity: usize,
        pool: Arc<WorkerPool>,
    ) -> Result<Self, RealGammaError> {
        if source.basis_momenta() != basis.momenta() {
            return Err(RealGammaError::Shape);
        }
        let grid = source.fft_grid();
        let half = HalfSphere::new(basis, grid)?;
        let map = PairMap::new(grid.dims(), grid.g_to_fft(), half.plus(), half.minus())?;
        let local = PairLocalWorkspace::new(map, capacity, pool.clone(), || {
            RealPairPlan::new(grid.dims(), grid.g_to_fft())
        })?;
        let mut kinetic = vec![0.; half.full_len()];
        kinetic[0] = source.kinetic_diagonal()[0];
        for (i, &p) in half.plus().iter().enumerate().skip(1) {
            kinetic[2 * i - 1] = source.kinetic_diagonal()[p as usize];
            kinetic[2 * i] = kinetic[2 * i - 1];
        }
        let nonlocal = RealNonlocal::new(source.nonlocal(), &half, capacity, pool.clone())?;
        let h = Self {
            source_identity: source.operator_identity(),
            half,
            kinetic,
            potential: source.local_potential_values().to_vec(),
            local,
            nonlocal,
            capacity,
            grid_shape: grid.dims(),
        };
        Ok(h)
    }
    /// Refresh only the local samples when immutable operator data match.
    /// False asks the owner to rebuild; no scratch or cached data is changed.
    pub fn update_from_complex(&mut self, source: &Hamiltonian, basis: &PlaneWaveBasis) -> bool {
        if !Arc::ptr_eq(&self.source_identity, &source.operator_identity())
            || source.basis_momenta() != basis.momenta()
            || basis.k_fractional() != [0.; 3]
            || self.potential.len() != source.local_potential_values().len()
        {
            return false;
        }
        self.potential
            .copy_from_slice(source.local_potential_values());
        true
    }
    /// ORC mapping for conversion at the library boundary.
    pub fn half_sphere(&self) -> &HalfSphere {
        &self.half
    }
    /// Overwrite Y with H X, reusing all FFT, projection and GEMM scratch.
    pub fn apply_block(
        &mut self,
        x: &RealMatrix,
        y: &mut RealMatrix,
    ) -> Result<(), RealGammaError> {
        let clock = timers::start();
        let result = self.apply_block_untimed(x, y);
        record_apply(clock, x.cols());
        result
    }
    fn apply_block_untimed(
        &mut self,
        x: &RealMatrix,
        y: &mut RealMatrix,
    ) -> Result<(), RealGammaError> {
        if x.rows() != self.half.full_len()
            || y.rows() != x.rows()
            || y.cols() != x.cols()
            || x.cols() > self.capacity
        {
            return Err(RealGammaError::Shape);
        }
        if x.as_slice().iter().any(|x| !x.is_finite()) {
            return Err(RealGammaError::NonFinite);
        }
        {}
        self.apply_local(x, y)?;
        {
            self.add_kinetic(x, y);
            self.nonlocal.add_block_core(x, y)
        }
    }
    /// Local part into Y (overwrite): pairs on the leading even columns, and
    /// with HAPPLY-FOLD the odd tail through the folded plan; the two run as
    /// the branches of one pool join (disjoint outputs, serial arithmetic each).
    fn apply_local(&mut self, x: &RealMatrix, y: &mut RealMatrix) -> Result<(), RealGammaError> {
        let clock = timers::start();
        let result = self
            .local
            .apply_block(x.as_slice(), &self.potential, y.as_mut_slice())
            .map_err(RealGammaError::from);
        timers::stop(clock, Slot::LocalNs);
        result
    }
    /// APPLY-LEAN kinetic term: the same `y += t x` per element, column by column
    /// instead of an integer modulo per element (bit-identical).
    fn add_kinetic(&self, x: &RealMatrix, y: &mut RealMatrix) {
        let n = self.kinetic.len();
        for (xc, yc) in x.as_slice().chunks(n).zip(y.as_mut_slice().chunks_mut(n)) {
            for ((y, &x), &t) in yc.iter_mut().zip(xc).zip(&self.kinetic) {
                *y += t * x;
            }
        }
    }
    /// Local potential samples of the current snapshot (grid C order, Ha).
    pub fn potential(&self) -> &[f64] {
        &self.potential
    }
    /// FFT grid dimensions of the local part.
    pub fn grid_dims(&self) -> [usize; 3] {
        self.grid_shape
    }
    /// REUSE-PSI: overwrite Y with the local image `potential * X` of an ORC
    /// block for a caller-supplied potential (the potential change between two
    /// SCF maps), without kinetic or nonlocal parts. Pairs (0,1), (2,3), ...
    /// whose `kept` flag is set use `grids[p]`, the unscaled inverse transform
    /// psi_a + i psi_b of exactly those two columns (kept from the density
    /// build), and need one forward transform; the others run the full pair
    /// transform. Equal to the local part of `apply_block` with `potential` to
    /// rounding. The snapshot potential is not used or changed.
    pub fn apply_local_kept(
        &mut self,
        x: &RealMatrix,
        grids: &[Vec<Complex64>],
        kept: &[bool],
        potential: &[f64],
        y: &mut RealMatrix,
    ) -> Result<(), RealGammaError> {
        if x.rows() != self.half.full_len()
            || y.rows() != x.rows()
            || y.cols() != x.cols()
            || x.cols() > self.capacity
            || potential.len() != self.potential.len()
        {
            return Err(RealGammaError::Shape);
        }
        if potential.iter().any(|v| !v.is_finite()) {
            return Err(RealGammaError::NonFinite);
        }
        self.local
            .apply_block_kept(x.as_slice(), grids, kept, potential, y.as_mut_slice())?;
        Ok(())
    }
    /// Number of nonlocal projector columns.
    pub fn projector_count(&self) -> usize {
        self.nonlocal.projector_count()
    }
    /// NL-ENERGY-REAL: per-column nonlocal expectation values (Ha) of an ORC
    /// block; see `RealNonlocal::expectation_values`.
    pub fn nonlocal_expectation_values(
        &mut self,
        x: &RealMatrix,
        out: &mut [f64],
    ) -> Result<(), RealGammaError> {
        if x.rows() != self.half.full_len() {
            return Err(RealGammaError::Shape);
        }
        self.nonlocal.expectation_values(x, out)
    }
    /// C3a: per-column nonlocal expectation values (Ha) from projections P^T X
    /// the caller holds; see `RealNonlocal::expectation_values_from_projections`.
    pub fn nonlocal_expectation_values_from_projections(
        &self,
        projections: &[f64],
        out: &mut [f64],
    ) -> Result<(), RealGammaError> {
        self.nonlocal
            .expectation_values_from_projections(projections, out)
    }
    /// Projections P^T X of the most recent `apply_block` (first `m` columns,
    /// column-major); see `RealNonlocal::last_projections`.
    pub fn last_projections(&self, m: usize) -> Result<&[f64], RealGammaError> {
        self.nonlocal.last_projections(m)
    }
    /// B6: overwrite Y with H X where the caller supplies P^T X (for X = V C,
    /// the cached (P^T V) C). Local and kinetic parts are applied as in
    /// `apply_block`; only the first nonlocal product is skipped. Equal to
    /// `apply_block` to rounding, not bit for bit.
    pub fn apply_block_with_projections(
        &mut self,
        x: &RealMatrix,
        projections: &[f64],
        y: &mut RealMatrix,
    ) -> Result<(), RealGammaError> {
        let clock = timers::start();
        let result = self.apply_block_with_projections_untimed(x, projections, y);
        record_apply(clock, x.cols());
        result
    }
    fn apply_block_with_projections_untimed(
        &mut self,
        x: &RealMatrix,
        projections: &[f64],
        y: &mut RealMatrix,
    ) -> Result<(), RealGammaError> {
        if x.rows() != self.half.full_len()
            || y.rows() != x.rows()
            || y.cols() != x.cols()
            || x.cols() > self.capacity
        {
            return Err(RealGammaError::Shape);
        }
        if x.as_slice().iter().any(|x| !x.is_finite()) {
            return Err(RealGammaError::NonFinite);
        }
        {}
        self.apply_local(x, y)?;
        {
            self.add_kinetic(x, y);
        }
        self.nonlocal.add_block_with_projections(projections, y)
    }
}

/// P0.2 work timers: one whole block application of `apply_block` or
/// `apply_block_with_projections` (no-op when the timers are off).
fn record_apply(clock: Option<std::time::Instant>, columns: usize) {
    if clock.is_some() {
        timers::stop(clock, Slot::ApplyNs);
        timers::add(Slot::ApplyCalls, 1);
        timers::add(Slot::ApplyColumns, columns as u64);
    }
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn projector_origin_guard() {
        let cell =
            dft_geometry::Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]])
                .unwrap();
        let b = PlaneWaveBasis::new(&cell, dft_units::Hartree(1.)).unwrap();
        let g = dft_basis::FftGridSpec::new(&b, &b, None).unwrap();
        let h = HalfSphere::new(&b, &g).unwrap();
        let mut beta = vec![Complex64::new(1., 0.); h.full_len()];
        let mut out = vec![0.; h.full_len()];
        beta[0].im = 2e-14;
        assert!(matches!(
            pack_projector(&h, &beta, &mut out),
            Err(RealGammaError::ComplexZero)
        ));
        beta[0].im = 0.5e-14;
        pack_projector(&h, &beta, &mut out).unwrap();
        beta.fill(Complex64::default());
        pack_projector(&h, &beta, &mut out).unwrap();
        beta[0].re = f64::NAN;
        assert!(pack_projector(&h, &beta, &mut out).is_err());
    }

    /// ISSUE-20: owned-row packing equals the rows of the full packing bit
    /// for bit, and gives the same verdict for every row subset (rank
    /// independent), including invalid columns.
    #[test]
    fn projector_rows_match_full_packing() {
        let cell = dft_geometry::Cell::from_user_rows_bohr([
            [6.2, 0., 0.],
            [1.1, 5.7, 0.],
            [0.4, 0.7, 6.6],
        ])
        .unwrap();
        let b = PlaneWaveBasis::new(&cell, dft_units::Hartree(1.5)).unwrap();
        let g = dft_basis::FftGridSpec::new(&b, &b, None).unwrap();
        let h = HalfSphere::new(&b, &g).unwrap();
        let n = h.full_len();
        let mut beta: Vec<Complex64> = (0..n)
            .map(|k| Complex64::new(((k * 7 + 1) as f64).sin(), ((k * 3 + 2) as f64).cos()))
            .collect();
        beta[0].im = 0.;
        let subsets: Vec<Vec<usize>> = vec![
            (0..n).collect(),
            (0..n).filter(|r| r % 3 == 1).collect(),
            (1..n).step_by(5).collect(),
            vec![0],
            Vec::new(),
        ];
        let check = |beta: &[Complex64]| {
            let mut full = vec![0.; n];
            let whole = pack_projector(&h, beta, &mut full).map_err(|e| e.to_string());
            for rows in &subsets {
                let mut out = vec![0.; rows.len()];
                let part = pack_projector_rows(&h, beta, rows, &mut out).map_err(|e| e.to_string());
                assert_eq!(whole.is_ok(), part.is_ok());
                assert_eq!(whole, part);
                if whole.is_ok() {
                    for (o, &r) in out.iter().zip(rows) {
                        assert_eq!(o.to_bits(), full[r].to_bits(), "row {r}");
                    }
                }
            }
        };
        check(&beta);
        // Invalid columns: imaginary origin, NaN, overflow of the sqrt(2) scaling.
        let mut bad = beta.clone();
        bad[0].im = 1e-3;
        check(&bad);
        let mut bad = beta.clone();
        bad[h.plus()[2] as usize].im = f64::NAN;
        check(&bad);
        let mut bad = beta.clone();
        bad[h.plus()[3] as usize] = Complex64::new(f64::MAX * 0.9, 0.);
        check(&bad);
        let mut bad = beta.clone();
        bad[h.plus()[3] as usize] = Complex64::new(f64::MAX * 0.6, f64::MAX * 0.6);
        check(&bad);
    }
}
