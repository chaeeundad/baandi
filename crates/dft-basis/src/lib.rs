//! Full complex Bloch plane waves and alias-free density FFT grids.
pub mod half;
use dft_geometry::Cell;
use dft_units::Hartree;
pub use half::{HalfSphere, HalfSphereError, ORC_TOLERANCE};
use ndarray::Array3;
use num_complex::Complex64;
use std::f64::consts::TAU;
use thiserror::Error;

/// Relative inclusion allowance at the kinetic-energy cutoff.
pub const CUTOFF_TOLERANCE: f64 = 1e-12;

#[derive(Debug, Clone, PartialEq, Error)]
pub enum BasisError {
    #[error("cutoff must be finite, positive, and representable in Hartree")]
    InvalidCutoff,
    #[error("nonfinite fractional k point")]
    InvalidKPoint,
    #[error("basis bounds or FFT storage size cannot be represented")]
    SizeOverflow,
    #[error("FFT override {provided:?} is smaller than required {required:?}")]
    GridTooSmall {
        provided: [usize; 3],
        required: [usize; 3],
    },
    #[error("Miller index {miller:?} is outside the canonical frequencies for {dims:?}")]
    IndexOutOfRange { miller: [i32; 3], dims: [usize; 3] },
    #[error("coefficient count must be {expected}, got {actual}")]
    CoefficientCount { expected: usize, actual: usize },
    #[error("FFT array shape must be {expected:?}")]
    ArrayShape { expected: [usize; 3] },
}

#[derive(Debug, Clone, PartialEq)]
pub struct GVector {
    pub miller: [i32; 3],
    /// Cartesian components in Bohr^-1.
    pub cartesian: [f64; 3],
    /// Squared magnitude in Bohr^-2.
    pub norm_squared: f64,
}

#[derive(Debug, Clone)]
pub struct PlaneWaveBasis {
    cutoff: Hartree,
    vectors: Vec<GVector>,
    momenta: Vec<GVector>,
    k_fractional: [f64; 3],
    k_cartesian: [f64; 3],
    density_reference: Option<Box<Cell>>,
}

fn bounding_box(cell: &Cell, cutoff: Hartree) -> Result<[i32; 3], BasisError> {
    let radius = (2.0 * cutoff.0 * (1.0 + CUTOFF_TOLERANCE)).sqrt();
    if !cutoff.0.is_finite() || cutoff.0 <= 0.0 || !radius.is_finite() {
        return Err(BasisError::InvalidCutoff);
    }
    let a = cell.lattice_columns_bohr();
    let mut bounds = [0; 3];
    for i in 0..3 {
        // Duality gives m_i = a_i dot G / (2*pi). Cauchy-Schwarz bounds
        // |m_i| by |a_i|*Gmax/(2*pi), including for skewed cells. Ceil
        // and one extra integer protect enumeration against rounding.
        let bound = (a[0][i].hypot(a[1][i]).hypot(a[2][i]) * radius / TAU).ceil() + 1.0;
        if !bound.is_finite() || bound >= i32::MAX as f64 {
            return Err(BasisError::SizeOverflow);
        }
        bounds[i] = bound as i32;
    }
    Ok(bounds)
}

/// The original enumeration: every Miller index of the box `bounds`, in
/// (m1, m2, m3) ascending order, with the inclusion test of `at_k`. Kept as
/// the reference of [`for_each_in_sphere`] and as its fallback.
fn for_each_in_box(
    b: [[f64; 3]; 3],
    cutoff: Hartree,
    bounds: [i32; 3],
    k_cartesian: [f64; 3],
    mut f: impl FnMut(GVector),
) -> Result<(), BasisError> {
    for m1 in -bounds[0]..=bounds[0] {
        for m2 in -bounds[1]..=bounds[1] {
            for m3 in -bounds[2]..=bounds[2] {
                include(b, cutoff, [m1, m2, m3], k_cartesian, &mut f)?;
            }
        }
    }
    Ok(())
}

/// Inclusion test of one Miller index (the arithmetic of the original loop).
#[inline]
fn include(
    b: [[f64; 3]; 3],
    cutoff: Hartree,
    miller: [i32; 3],
    k_cartesian: [f64; 3],
    f: &mut impl FnMut(GVector),
) -> Result<(), BasisError> {
    let cartesian =
        std::array::from_fn(|i| (0..3).map(|j| b[i][j] * f64::from(miller[j])).sum::<f64>());
    let norm_squared = cartesian.iter().map(|x| x * x).sum::<f64>();
    if !norm_squared.is_finite() {
        return Err(BasisError::SizeOverflow);
    }
    let q2: f64 = (0..3)
        .map(|i| (cartesian[i] + k_cartesian[i]).powi(2))
        .sum();
    if 0.5 * q2 <= cutoff.0 * (1.0 + CUTOFF_TOLERANCE) {
        f(GVector {
            miller,
            cartesian,
            norm_squared,
        });
    }
    Ok(())
}

/// GK18: the admissible Miller indices of the box, in the box order, with the
/// same per-index arithmetic as [`for_each_in_box`], but for each (m1, m2)
/// only an m3 window is tested. With c = m1 b1 + m2 b2 + k, |c + t b3|^2 =
/// d^2 + |b3|^2 (t - t0)^2; the window is [t0 - w - 2, t0 + w + 2] with
/// |b3|^2 w^2 = max(0, 2 E (1 + tol) - d^2), so every index outside it has
/// |k + G|^2 >= 2 E (1 + tol) + 4 |b3|^2 exactly; the guard below keeps
/// 4 |b3|^2 far above the rounding of the test, so the skipped indices
/// would all have been rejected and the output (set and order) is unchanged.
/// The finiteness error of the box loop cannot occur when the largest
/// possible |G|^2 of the box is finite (checked first); otherwise, and for
/// degenerate metrics, the whole box is enumerated.
fn for_each_in_sphere(
    cell: &Cell,
    cutoff: Hartree,
    bounds: [i32; 3],
    k_cartesian: [f64; 3],
    mut f: impl FnMut(GVector),
) -> Result<(), BasisError> {
    let b = cell.reciprocal_columns_inverse_bohr();
    let limit = 2.0 * cutoff.0 * (1.0 + CUTOFF_TOLERANCE);
    let b3 = [b[0][2], b[1][2], b[2][2]];
    let b3sq: f64 = b3.iter().map(|x| x * x).sum();
    let largest: f64 = (0..3)
        .map(|i| {
            (0..3)
                .map(|j| b[i][j].abs() * f64::from(bounds[j]))
                .sum::<f64>()
                + k_cartesian[i].abs()
        })
        .map(|x| x * x)
        .sum();
    let prune = largest.is_finite()
        && limit.is_finite()
        && b3sq.is_finite()
        && 4.0 * b3sq > 1e-6 * limit
        && b.iter().flatten().all(|x| x.is_finite());
    if !prune {
        return for_each_in_box(b, cutoff, bounds, k_cartesian, f);
    }
    for m1 in -bounds[0]..=bounds[0] {
        for m2 in -bounds[1]..=bounds[1] {
            let c: [f64; 3] = std::array::from_fn(|i| {
                b[i][0] * f64::from(m1) + b[i][1] * f64::from(m2) + k_cartesian[i]
            });
            let cb: f64 = (0..3).map(|i| c[i] * b3[i]).sum();
            let cc: f64 = c.iter().map(|x| x * x).sum();
            let t0 = -cb / b3sq;
            let d2 = cc - cb * cb / b3sq;
            let w = ((limit - d2).max(0.0) / b3sq).sqrt();
            let lo = (t0 - w - 2.0).floor().max(-f64::from(bounds[2]));
            let hi = (t0 + w + 2.0).ceil().min(f64::from(bounds[2]));
            if lo > hi {
                continue;
            }
            for m3 in lo as i32..=hi as i32 {
                include(b, cutoff, [m1, m2, m3], k_cartesian, &mut f)?;
            }
        }
    }
    Ok(())
}

/// GK18: what callers need of a k basis without building it: the plane-wave
/// count, the Miller extent per axis and the largest |k+G|^2 (the value
/// `PlaneWaveBasis::momenta` holds, same arithmetic), from the enumeration of
/// `at_k` (no allocation, no sort).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KBasisSummary {
    pub count: usize,
    pub miller_min: [i32; 3],
    pub miller_max: [i32; 3],
    pub max_momentum_norm_squared: f64,
    pub k_fractional: [f64; 3],
}

impl KBasisSummary {
    /// Same input checks and canonical k as [`PlaneWaveBasis::at_k`].
    pub fn at_k(cell: &Cell, cutoff: Hartree, fractional: [f64; 3]) -> Result<Self, BasisError> {
        let k_fractional = dft_geometry::kpoints::wrap_fractional(fractional)
            .map_err(|_| BasisError::InvalidKPoint)?;
        let bounds = bounding_box(cell, cutoff)?.map(|b| b + 1);
        let b = cell.reciprocal_columns_inverse_bohr();
        let k_cartesian: [f64; 3] =
            std::array::from_fn(|i| (0..3).map(|j| b[i][j] * k_fractional[j]).sum());
        let mut s = Self {
            count: 0,
            miller_min: [i32::MAX; 3],
            miller_max: [i32::MIN; 3],
            max_momentum_norm_squared: 0.0,
            k_fractional,
        };
        for_each_in_sphere(cell, cutoff, bounds, k_cartesian, |g| {
            s.count += 1;
            for i in 0..3 {
                s.miller_min[i] = s.miller_min[i].min(g.miller[i]);
                s.miller_max[i] = s.miller_max[i].max(g.miller[i]);
            }
            let cartesian: [f64; 3] = std::array::from_fn(|i| g.cartesian[i] + k_cartesian[i]);
            let n: f64 = cartesian.iter().map(|x| x * x).sum();
            s.max_momentum_norm_squared = s.max_momentum_norm_squared.max(n);
        })?;
        Ok(s)
    }
}

impl PlaneWaveBasis {
    pub fn new(cell: &Cell, cutoff: Hartree) -> Result<Self, BasisError> {
        Self::at_k(cell, cutoff, [0.0; 3])
    }
    /// Bloch basis with |k+G|^2/2 <= cutoff (Ha). Fractional reciprocal k is
    /// canonicalized modulo integers to [-1/2,1/2); internally k is in Bohr^-1.
    /// Coefficients describe the periodic part of the Bloch wavefunction.
    pub fn at_k(cell: &Cell, cutoff: Hartree, fractional: [f64; 3]) -> Result<Self, BasisError> {
        let k_fractional = dft_geometry::kpoints::wrap_fractional(fractional)
            .map_err(|_| BasisError::InvalidKPoint)?;
        let bounds = bounding_box(cell, cutoff)?.map(|b| b + 1);
        let b = cell.reciprocal_columns_inverse_bohr();
        let k_cartesian: [f64; 3] =
            std::array::from_fn(|i| (0..3).map(|j| b[i][j] * k_fractional[j]).sum());
        let mut vectors = Vec::new();
        for_each_in_sphere(cell, cutoff, bounds, k_cartesian, |g| {
            vectors.push(g);
        })?;
        vectors.sort_by(|a, b| {
            a.norm_squared
                .total_cmp(&b.norm_squared)
                .then(a.miller.cmp(&b.miller))
        });
        let momenta = vectors
            .iter()
            .map(|g| {
                let cartesian = std::array::from_fn(|i| g.cartesian[i] + k_cartesian[i]);
                GVector {
                    miller: g.miller,
                    norm_squared: cartesian.iter().map(|x| x * x).sum(),
                    cartesian,
                }
            })
            .collect();
        Ok(Self {
            cutoff,
            vectors,
            momenta,
            k_fractional,
            k_cartesian,
            density_reference: None,
        })
    }
    /// Deform the reciprocal metric while retaining the exact ordered Miller
    /// set and fractional k. The recorded cutoff is the reference cutoff,
    /// not an assertion that the deformed vectors obey a cutoff sphere.
    pub fn in_cell(&self, cell: &Cell) -> Result<Self, BasisError> {
        let b = cell.reciprocal_columns_inverse_bohr();
        let cart = |m: [f64; 3]| std::array::from_fn(|i| (0..3).map(|j| b[i][j] * m[j]).sum());
        let k_cartesian = cart(self.k_fractional);
        let make = |g: &GVector, shift: [f64; 3]| {
            let cartesian: [f64; 3] =
                cart(std::array::from_fn(|j| f64::from(g.miller[j]) + shift[j]));
            GVector {
                miller: g.miller,
                norm_squared: cartesian.iter().map(|v| v * v).sum(),
                cartesian,
            }
        };
        let vectors: Vec<_> = self.vectors.iter().map(|g| make(g, [0.; 3])).collect();
        let momenta: Vec<_> = self
            .vectors
            .iter()
            .map(|g| make(g, self.k_fractional))
            .collect();
        if vectors
            .iter()
            .chain(&momenta)
            .any(|g| !g.norm_squared.is_finite())
        {
            return Err(BasisError::SizeOverflow);
        }
        Ok(Self {
            cutoff: self.cutoff,
            vectors,
            momenta,
            k_fractional: self.k_fractional,
            k_cartesian,
            density_reference: self.density_reference.clone(),
        })
    }
    /// Deform wavefunctions and retain the supplied reference density sphere.
    /// Hamiltonian setup must use the same frozen density support as the engine,
    /// including when deformation crosses a nominal-cutoff Miller shell.
    pub fn in_cell_with_frozen_density(
        &self,
        cell: &Cell,
        density_reference: &Cell,
    ) -> Result<Self, BasisError> {
        let mut deformed = self.in_cell(cell)?;
        deformed.density_reference = Some(Box::new(density_reference.clone()));
        Ok(deformed)
    }
    /// Cartesian k+G and squared magnitude, in the same ordering as `vectors`.
    pub fn momenta(&self) -> &[GVector] {
        &self.momenta
    }
    /// Canonical fractional reciprocal k.
    pub fn k_fractional(&self) -> [f64; 3] {
        self.k_fractional
    }
    /// Cartesian k in Bohr^-1.
    pub fn k_cartesian(&self) -> [f64; 3] {
        self.k_cartesian
    }
    pub fn cutoff(&self) -> Hartree {
        self.cutoff
    }
    /// Reciprocal lattice vectors G (not k+G), sorted by G norm then Miller index.
    pub fn vectors(&self) -> &[GVector] {
        &self.vectors
    }
    /// Unshifted density sphere at four times Ecut, independent of k.
    pub fn density_basis(&self, cell: &Cell) -> Result<Self, BasisError> {
        match &self.density_reference {
            Some(reference) => Self::new(reference, Hartree(4.0 * self.cutoff.0))?.in_cell(cell),
            None => Self::new(cell, Hartree(4.0 * self.cutoff.0)),
        }
    }
}

/// Immutable grid and forward/inverse mappings for one basis.
/// Unoccupied FFT entries map to None in `fft_to_g`.
#[derive(Debug, Clone)]
pub struct FftGridSpec {
    dims: [usize; 3],
    required: [usize; 3],
    overridden: bool,
    g_to_fft: Vec<[usize; 3]>,
    fft_to_g: Vec<Option<usize>>,
}

fn smooth(mut n: usize) -> bool {
    for p in [2, 3, 5] {
        while n.is_multiple_of(p) {
            n /= p;
        }
    }
    n == 1
}

impl FftGridSpec {
    /// `density` determines the minimum grid; `basis` determines coefficient ordering.
    pub fn new(
        basis: &PlaneWaveBasis,
        density: &PlaneWaveBasis,
        override_dims: Option<[usize; 3]>,
    ) -> Result<Self, BasisError> {
        let (dims, required, size) = Self::grid_dims(density, override_dims)?;
        let g_to_fft = basis
            .vectors()
            .iter()
            .map(|g| {
                if basis.k_fractional == [0.0; 3] {
                    Self::map_miller(g.miller, dims)
                } else {
                    // A shifted sphere can have an absolute G outside the
                    // canonical FFT interval in a very skewed cell. Only
                    // differences G-G' enter density/local matrix elements;
                    // those remain inside the unshifted density sphere.
                    Ok(std::array::from_fn(|i| {
                        i128::from(g.miller[i]).rem_euclid(dims[i] as i128) as usize
                    }))
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        let mut fft_to_g = vec![None; size];
        for (g, &[i, j, k]) in g_to_fft.iter().enumerate() {
            let index = (i * dims[1] + j) * dims[2] + k;
            if fft_to_g[index].is_some() {
                return Err(BasisError::IndexOutOfRange {
                    miller: basis.vectors()[g].miller,
                    dims,
                });
            }
            fft_to_g[index] = Some(g);
        }
        Ok(Self {
            dims,
            required,
            overridden: override_dims.is_some(),
            g_to_fft,
            fft_to_g,
        })
    }
    /// Required and selected dimensions and the grid size of `new` (its
    /// checks, in its order); independent of the wavefunction basis.
    fn grid_dims(
        density: &PlaneWaveBasis,
        override_dims: Option<[usize; 3]>,
    ) -> Result<([usize; 3], [usize; 3], usize), BasisError> {
        let mut required = [1; 3];
        for g in density.vectors() {
            for (i, r) in required.iter_mut().enumerate() {
                *r = (*r).max(2 * g.miller[i].unsigned_abs() as usize + 1);
            }
        }
        let dims = match override_dims {
            Some(dims) => {
                if (0..3).any(|i| dims[i] < required[i]) {
                    return Err(BasisError::GridTooSmall {
                        provided: dims,
                        required,
                    });
                }
                dims
            }
            None => {
                let mut dims = required;
                for d in &mut dims {
                    while !smooth(*d) {
                        *d = d.checked_add(1).ok_or(BasisError::SizeOverflow)?;
                    }
                }
                dims
            }
        };
        let size = dims
            .iter()
            .try_fold(1usize, |s, &n| s.checked_mul(n))
            .ok_or(BasisError::SizeOverflow)?;
        if size
            > isize::MAX as usize
                / std::mem::size_of::<Complex64>().max(std::mem::size_of::<Option<usize>>())
        {
            return Err(BasisError::SizeOverflow);
        }
        Ok((dims, required, size))
    }

    /// GK18: `FftGridSpec::new(basis, density, override_dims)?.dims()` for a k
    /// basis known by its [`KBasisSummary`]: the same dimension checks, and the
    /// coefficient map is proven valid from the Miller extents (canonical
    /// interval at Gamma, a window shorter than the grid otherwise, so no two
    /// coefficients share an FFT index). When the extents do not prove it, the
    /// basis is built (`build`) and `new` decides, with its exact error.
    pub fn dims_for_summary(
        summary: &KBasisSummary,
        density: &PlaneWaveBasis,
        override_dims: Option<[usize; 3]>,
        build: impl FnOnce() -> Result<PlaneWaveBasis, BasisError>,
    ) -> Result<[usize; 3], BasisError> {
        let (dims, _, _) = Self::grid_dims(density, override_dims)?;
        let fits = summary.count == 0
            || (0..3).all(|i| {
                let n = dims[i] as i64;
                let (lo, hi) = (
                    i64::from(summary.miller_min[i]),
                    i64::from(summary.miller_max[i]),
                );
                if summary.k_fractional == [0.0; 3] {
                    lo >= -(n / 2) && hi <= (n - 1) / 2
                } else {
                    hi - lo < n
                }
            });
        if fits {
            Ok(dims)
        } else {
            Ok(Self::new(&build()?, density, override_dims)?.dims())
        }
    }

    /// Canonical interval is [-n/2, (n-1)/2]; even grids use negative Nyquist.
    pub fn map_miller(miller: [i32; 3], dims: [usize; 3]) -> Result<[usize; 3], BasisError> {
        let mut index = [0; 3];
        for i in 0..3 {
            let n = dims[i] as i128;
            let m = i128::from(miller[i]);
            if n == 0 || m < -(n / 2) || m > (n - 1) / 2 {
                return Err(BasisError::IndexOutOfRange { miller, dims });
            }
            index[i] = if m >= 0 { m } else { n + m } as usize;
        }
        Ok(index)
    }
    pub fn dims(&self) -> [usize; 3] {
        self.dims
    }
    pub fn required_dims(&self) -> [usize; 3] {
        self.required
    }
    /// Construct an automatically selected grid with additional constraints
    /// (for example crystal symmetry), retaining automatic provenance.
    /// Dimensions undergo the same alias-free checks as an explicit grid.
    pub fn automatic_with_dims(
        basis: &PlaneWaveBasis,
        density: &PlaneWaveBasis,
        dims: [usize; 3],
    ) -> Result<Self, BasisError> {
        let mut grid = Self::new(basis, density, Some(dims))?;
        grid.overridden = false;
        Ok(grid)
    }
    /// True for user-specified dimensions, false for automatic selection.
    pub fn overridden(&self) -> bool {
        self.overridden
    }
    pub fn g_to_fft(&self) -> &[[usize; 3]] {
        &self.g_to_fft
    }
    /// C-order flat FFT index -> coefficient index, or None outside this basis.
    pub fn fft_to_g(&self) -> &[Option<usize>] {
        &self.fft_to_g
    }
    pub fn scatter(&self, coefficients: &[Complex64]) -> Result<Array3<Complex64>, BasisError> {
        if coefficients.len() != self.g_to_fft.len() {
            return Err(BasisError::CoefficientCount {
                expected: self.g_to_fft.len(),
                actual: coefficients.len(),
            });
        }
        let mut array = Array3::zeros(self.dims);
        for (&index, &c) in self.g_to_fft.iter().zip(coefficients) {
            array[index] = c;
        }
        Ok(array)
    }
    pub fn gather(&self, array: &Array3<Complex64>) -> Result<Vec<Complex64>, BasisError> {
        if array.shape() != self.dims {
            return Err(BasisError::ArrayShape {
                expected: self.dims,
            });
        }
        Ok(self.g_to_fft.iter().map(|&index| array[index]).collect())
    }
}

#[cfg(test)]
mod tests;
