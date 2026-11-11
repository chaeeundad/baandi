//! Isolated real Gamma Hamiltonian. Derived from real-gamma-design sections 4–7;
//! no engine caller or change to the complex traversal (ledger A43).
use crate::hamiltonian::{Hamiltonian, NonlocalOperator};
use dft_basis::{HalfSphere, PlaneWaveBasis};
use dft_numerics::fused_fold::FoldedRealPlan;
use dft_numerics::fused_pair_plan::RealPairPlan;
use dft_numerics::{
    FftPlan3d, PairMap, RealGemmWorkspace, RealMatrix, RealOp, fused_pair::PairLocalWorkspace,
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

struct Coupling {
    start: usize,
    count: usize,
    values: Vec<f64>,
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
    /// MKL builds, default ON since the user's adopt-what-improves rule (2026-11-10; dev08 NONLOCAL-DGEMM-AB: all 40 seed pairs faster, -6.4/-4.2 percent at T=1, -2.0/-1.5 percent at T=8); `BAANDI_MKL_NONLOCAL=off` restores the Rust kernel: both
    /// projector products go through dgemm panels with the panel grouping fixed.
    dgemm_products: bool,
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
        Ok(Self {
            projectors,
            coupling,
            projections: RealMatrix::zeros(count, capacity)?,
            weights: RealMatrix::zeros(count, capacity)?,
            image: RealMatrix::zeros(n, capacity)?,
            gemm: RealGemmWorkspace::with_pool(n, count.max(capacity), pool)?,
            capacity,
            dgemm_products: dft_numerics::switches::on(dft_numerics::switches::Id::MKL_NONLOCAL),
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
        if self.dgemm_products {
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
        } else {
            self.gemm.multiply_into(
                self.projectors.as_ref(),
                RealOp::Transpose,
                x.as_ref(),
                self.projections.as_mut().subcols_mut(0, m),
            )?;
        }
        self.add_from_projections(y, m)
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

    /// B6 (`BAANDI_B6=on` callers): add the nonlocal image of a block whose
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
        self.projections.as_mut_slice()[..count * m].copy_from_slice(projections);
        self.add_from_projections(y, m)
    }

    /// NL-ENERGY-REAL (engine switch `BAANDI_NL_ENERGY_REAL=on`): per column the
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
        if self.dgemm_products {
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
        } else {
            self.gemm.multiply_into(
                self.projectors.as_ref(),
                RealOp::Transpose,
                x.as_ref(),
                self.projections.as_mut().subcols_mut(0, m),
            )?;
        }
        for (band, value) in out.iter_mut().enumerate() {
            let p = &self.projections.as_slice()[band * count..(band + 1) * count];
            let mut energy = 0.0;
            for atom in &self.coupling {
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

    /// Coupling and back product from `self.projections` (first `m` columns).
    fn add_from_projections(&mut self, y: &mut RealMatrix, m: usize) -> Result<(), RealGammaError> {
        self.image_from_projections(m)?;
        self.add_image(y)
    }

    /// HAPPLY-OVERLAP: P^T X, coupling and back product into the private image
    /// (the arithmetic of `add_block_core` without touching Y).
    fn image_block(&mut self, x: &RealMatrix) -> Result<(), RealGammaError> {
        let n = self.projectors.rows();
        let m = x.cols();
        if x.rows() != n || m > self.capacity {
            return Err(RealGammaError::Shape);
        }
        let count = self.projectors.cols();
        if self.dgemm_products {
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
        } else {
            self.gemm.multiply_into(
                self.projectors.as_ref(),
                RealOp::Transpose,
                x.as_ref(),
                self.projections.as_mut().subcols_mut(0, m),
            )?;
        }
        self.image_from_projections(m)
    }

    /// HAPPLY-OVERLAP with caller projections (B6): the arithmetic of
    /// `add_block_with_projections` without touching Y.
    fn image_with_projections(
        &mut self,
        projections: &[f64],
        m: usize,
    ) -> Result<(), RealGammaError> {
        let count = self.projectors.cols();
        if m > self.capacity || projections.len() != count * m {
            return Err(RealGammaError::Shape);
        }
        if projections.iter().any(|x| !x.is_finite()) {
            return Err(RealGammaError::NonFinite);
        }
        self.projections.as_mut_slice()[..count * m].copy_from_slice(projections);
        self.image_from_projections(m)
    }

    /// Y += image (first `Y.cols()` columns), then the output scan.
    fn add_image(&self, y: &mut RealMatrix) -> Result<(), RealGammaError> {
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
        if self.dgemm_products {
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
        } else {
            self.gemm.multiply_into(
                self.projectors.as_ref(),
                RealOp::Identity,
                self.weights.as_ref().subcols(0, m),
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
    /// APPLY-LEAN (opt-in `BAANDI_APPLY_LEAN=on`, default off; bit-identical):
    /// kinetic term without a per-element modulo and no repeated input scans in
    /// the nonlocal add (X is scanned once at entry, Y once at the end).
    lean: bool,
    pool: Arc<WorkerPool>,
    /// HAPPLY-OVERLAP (opt-in `BAANDI_HAPPLY_OVERLAP=on`, default off;
    /// bit-identical): the local part and the nonlocal image run as the two
    /// branches of one pool join (the nonlocal GEMMs keep their fixed panel
    /// grouping), then Y += kinetic and Y += image in the serial order.
    overlap: bool,
    /// HAPPLY-FOLD (opt-in `BAANDI_HAPPLY_FOLD=on`, default off; rounding
    /// class for odd blocks, even blocks unchanged): the odd tail column goes
    /// through the half-length folded transform instead of a pair with a zero
    /// partner. `None` when off or when the last grid axis is odd.
    fold: Option<FoldedRealPlan<FftPlan3d>>,
    grid_shape: [usize; 3],
    grid_indices: Vec<[usize; 3]>,
}
impl RealHamiltonian {
    /// Test hook: force APPLY-LEAN on or off regardless of the environment.
    pub fn set_apply_lean(&mut self, on: bool) {
        self.lean = on;
    }
    /// Test hook: force HAPPLY-OVERLAP on or off regardless of the environment.
    pub fn set_apply_overlap(&mut self, on: bool) {
        self.overlap = on;
    }
    /// Test hook: force HAPPLY-FOLD on or off regardless of the environment.
    /// Returns whether the fold is active (false for an odd last grid axis).
    /// Building the plan allocates; call it outside numerical loops.
    pub fn set_apply_fold(&mut self, on: bool) -> Result<bool, RealGammaError> {
        self.fold = if on {
            match self.fold.take() {
                Some(plan) => Some(plan),
                None => FoldedRealPlan::new(
                    self.grid_shape,
                    &self.grid_indices,
                    self.half.plus(),
                    self.half.minus(),
                )?,
            }
        } else {
            None
        };
        Ok(self.fold.is_some())
    }
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
        let mut h = Self {
            source_identity: source.operator_identity(),
            half,
            kinetic,
            potential: source.local_potential_values().to_vec(),
            local,
            nonlocal,
            capacity,
            lean: dft_numerics::switches::on(dft_numerics::switches::Id::APPLY_LEAN),
            pool,
            overlap: dft_numerics::switches::on(dft_numerics::switches::Id::HAPPLY_OVERLAP),
            fold: None,
            grid_shape: grid.dims(),
            grid_indices: grid.g_to_fft().to_vec(),
        };
        h.set_apply_fold(dft_numerics::switches::on(
            dft_numerics::switches::Id::HAPPLY_FOLD,
        ))?;
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
        if self.overlap {
            return self.apply_overlapped(x, None, y);
        }
        self.apply_local(x, y)?;
        if self.lean {
            self.add_kinetic(x, y);
            return self.nonlocal.add_block_core(x, y);
        }
        for (j, (x, y)) in x.as_slice().iter().zip(y.as_mut_slice()).enumerate() {
            *y += self.kinetic[j % self.kinetic.len()] * x;
        }
        self.nonlocal.add_block(x, y)
    }
    /// Local part into Y (overwrite): pairs on the leading even columns, and
    /// with HAPPLY-FOLD the odd tail through the folded plan; the two run as
    /// the branches of one pool join (disjoint outputs, serial arithmetic each).
    fn apply_local(&mut self, x: &RealMatrix, y: &mut RealMatrix) -> Result<(), RealGammaError> {
        local_part(
            &mut self.local,
            self.fold.as_mut(),
            &self.pool,
            x,
            &self.potential,
            y,
        )
    }
    /// HAPPLY-OVERLAP: local part and nonlocal image in one pool join, then
    /// the kinetic and nonlocal additions in the serial order (bit-identical
    /// to the non-overlapped path for any thread budget).
    fn apply_overlapped(
        &mut self,
        x: &RealMatrix,
        projections: Option<&[f64]>,
        y: &mut RealMatrix,
    ) -> Result<(), RealGammaError> {
        let m = x.cols();
        let Self {
            local,
            fold,
            pool,
            potential,
            nonlocal,
            ..
        } = self;
        let (local_result, image_result) = pool.install(|| {
            WorkerPool::join(
                || local_part(local, fold.as_mut(), pool, x, potential, y),
                || match projections {
                    None => nonlocal.image_block(x),
                    Some(p) => nonlocal.image_with_projections(p, m),
                },
            )
        });
        local_result?;
        image_result?;
        self.add_kinetic(x, y);
        self.nonlocal.add_image(y)
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
        if self.overlap {
            return self.apply_overlapped(x, Some(projections), y);
        }
        self.apply_local(x, y)?;
        if self.lean {
            self.add_kinetic(x, y);
        } else {
            for (j, (x, y)) in x.as_slice().iter().zip(y.as_mut_slice()).enumerate() {
                *y += self.kinetic[j % self.kinetic.len()] * x;
            }
        }
        self.nonlocal.add_block_with_projections(projections, y)
    }
}

/// Local potential image of a block (overwrite). Without a fold plan this is
/// exactly `PairLocalWorkspace::apply_block`. With one, an odd block keeps the
/// pairs (0,1), (2,3), ... for its leading columns and sends the tail column
/// through the folded plan; the outputs are disjoint, so the join changes no bit.
fn local_part(
    local: &mut PairLocalWorkspace<RealPairPlan>,
    fold: Option<&mut FoldedRealPlan<FftPlan3d>>,
    pool: &WorkerPool,
    x: &RealMatrix,
    potential: &[f64],
    y: &mut RealMatrix,
) -> Result<(), RealGammaError> {
    let m = x.cols();
    let fold = match fold {
        Some(f) if m % 2 == 1 => f,
        _ => {
            local.apply_block(x.as_slice(), potential, y.as_mut_slice())?;
            return Ok(());
        }
    };
    let n = x.rows();
    let split = (m - 1) * n;
    let (xp, xt) = x.as_slice().split_at(split);
    let (yp, yt) = y.as_mut_slice().split_at_mut(split);
    if split == 0 {
        fold.apply_local_real(xt, potential, yt)?;
        return Ok(());
    }
    let (a, b) = pool.install(|| {
        WorkerPool::join(
            || local.apply_block(xp, potential, yp),
            || fold.apply_local_real(xt, potential, yt),
        )
    });
    a?;
    b?;
    Ok(())
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
}
