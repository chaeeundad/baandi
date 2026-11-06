//! Isolated real Gamma Hamiltonian. Derived from real-gamma-design sections 4–7;
//! no engine caller or change to the complex traversal (ledger A43).
use crate::hamiltonian::{Hamiltonian, NonlocalOperator};
use dft_basis::{HalfSphere, PlaneWaveBasis};
use dft_numerics::{
    PairMap, PairPlan, RealGemmWorkspace, RealMatrix, RealOp, fused::FusedPlan, pool::WorkerPool,
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

fn pack_projector(
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
        self.gemm.multiply_into(
            self.projectors.as_ref(),
            RealOp::Transpose,
            x.as_ref(),
            self.projections.as_mut().subcols_mut(0, m),
        )?;
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
        self.gemm.multiply_into(
            self.projectors.as_ref(),
            RealOp::Identity,
            self.weights.as_ref().subcols(0, m),
            self.image.as_mut().subcols_mut(0, m),
        )?;
        for (y, v) in y.as_mut_slice().iter_mut().zip(self.image.as_slice()) {
            *y += v;
        }
        if y.as_slice().iter().any(|x| !x.is_finite()) {
            return Err(RealGammaError::NonFinite);
        }
        Ok(())
    }
}

/// Snapshot of a complex Gamma H, acting on column-major ORC blocks.
/// Scratch capacity is fixed at construction; nonzero paired band norms must
/// differ by at most 1e3 (R2 debug assertion). Local pairs execute serially;
/// GEMM uses the supplied persistent pool. No SCF dispatcher is installed.
pub struct RealHamiltonian {
    half: HalfSphere,
    kinetic: Vec<f64>,
    potential: Vec<f64>,
    local: PairPlan<FusedPlan>,
    nonlocal: RealNonlocal,
    capacity: usize,
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
        let local = PairPlan::new(FusedPlan::new(grid.dims(), grid.g_to_fft())?, map)?;
        let mut kinetic = vec![0.; half.full_len()];
        kinetic[0] = source.kinetic_diagonal()[0];
        for (i, &p) in half.plus().iter().enumerate().skip(1) {
            kinetic[2 * i - 1] = source.kinetic_diagonal()[p as usize];
            kinetic[2 * i] = kinetic[2 * i - 1];
        }
        let nonlocal = RealNonlocal::new(source.nonlocal(), &half, capacity, pool)?;
        Ok(Self {
            half,
            kinetic,
            potential: source.local_potential_values().to_vec(),
            local,
            nonlocal,
            capacity,
        })
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
        self.local
            .apply_block(x.as_slice(), &self.potential, y.as_mut_slice())?;
        for (j, (x, y)) in x.as_slice().iter().zip(y.as_mut_slice()).enumerate() {
            *y += self.kinetic[j % self.kinetic.len()] * x;
        }
        self.nonlocal.add_block(x, y)
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
}
