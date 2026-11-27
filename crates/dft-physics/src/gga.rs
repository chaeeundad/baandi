//! Spin-unpolarized PBE with a discrete variational spectral gradient.
pub mod distributed;
use crate::{PhysicsError, xc::XcGrid};
use dft_geometry::Cell;
use dft_numerics::Fft3d;
use dft_units::Hartree;
use ndarray::{Array1, Array3};
use num_complex::Complex64;
use xcx::{Functional, FunctionalId, Spin, XcInput};

/// Bohr^-3. Above both xcx component screens; avoids singular vacuum arithmetic.
/// The existing LDA/QE-comparison candidate threshold is also 1e-10.
pub const DENSITY_THRESHOLD: f64 = 1e-10;

/// Cartesian spectral derivative. Zero each even-axis Nyquist Miller component
/// BEFORE multiplying by B, including on skew cells. This real operator is
/// skew-adjoint under the uniform-grid inner product.
fn derivative(
    n: &Array3<f64>,
    cell: &Cell,
    axis: usize,
    fft: &mut impl Fft3d,
) -> Result<Array3<f64>, PhysicsError> {
    let dims = fft.shape();
    let b = cell.reciprocal_columns_inverse_bohr();
    let mut reciprocal = Array3::zeros(n.raw_dim());
    fft.to_reciprocal(&n.mapv(|x| Complex64::new(x, 0.0)), &mut reciprocal)?;
    for ((i, j, k), z) in reciprocal.indexed_iter_mut() {
        let indices = [i, j, k];
        let g: f64 = (0..3)
            .map(|a| {
                let m = if dims[a].is_multiple_of(2) && indices[a] == dims[a] / 2 {
                    0.0
                } else if indices[a] <= (dims[a] - 1) / 2 {
                    indices[a] as f64
                } else {
                    indices[a] as f64 - dims[a] as f64
                };
                b[axis][a] * m
            })
            .sum();
        *z *= Complex64::new(0.0, g);
    }
    let mut real = Array3::zeros(n.raw_dim());
    fft.to_real(&reciprocal, &mut real)?;
    Ok(real.mapv(|x| x.re))
}

/// GD2: the three Cartesian gradients of [`derivative`] from ONE forward
/// transform of `n` (the forward transform of the same input is
/// deterministic, so every gradient keeps its bits).
fn gradients_once(
    n: &Array3<f64>,
    cell: &Cell,
    fft: &mut impl Fft3d,
) -> Result<Vec<Array3<f64>>, PhysicsError> {
    let dims = fft.shape();
    let b = cell.reciprocal_columns_inverse_bohr();
    let mut transformed = Array3::zeros(n.raw_dim());
    fft.to_reciprocal(&n.mapv(|x| Complex64::new(x, 0.0)), &mut transformed)?;
    let mut reciprocal = Array3::zeros(n.raw_dim());
    let mut real = Array3::zeros(n.raw_dim());
    let mut out = Vec::with_capacity(3);
    for row in &b {
        for (((i, j, k), z), &source) in reciprocal.indexed_iter_mut().zip(&transformed) {
            let indices = [i, j, k];
            let g: f64 = (0..3)
                .map(|a| {
                    let m = if dims[a].is_multiple_of(2) && indices[a] == dims[a] / 2 {
                        0.0
                    } else if indices[a] <= (dims[a] - 1) / 2 {
                        indices[a] as f64
                    } else {
                        indices[a] as f64 - dims[a] as f64
                    };
                    row[a] * m
                })
                .sum();
            *z = source;
            *z *= Complex64::new(0.0, g);
        }
        fft.to_real(&reciprocal, &mut real)?;
        out.push(real.mapv(|x| x.re));
    }
    Ok(out)
}

/// GD2: splits the pointwise functional evaluation of
/// [`evaluate_for_derivatives`] over participants (MPI ranks). Implementations
/// must be collective: every participant calls `gather` once per evaluation.
pub trait PointSplit {
    /// The contiguous block of the `len` grid points this participant evaluates.
    fn block(&self, len: usize) -> std::ops::Range<usize>;
    /// The `width` values per point of every participant's block (`local`,
    /// point-major), concatenated in point order (`len * width` values), on
    /// every participant. A participant with `failed` makes every participant
    /// return an error.
    fn gather(
        &self,
        local: &[f64],
        len: usize,
        width: usize,
        failed: bool,
    ) -> Result<Vec<f64>, PhysicsError>;
}

/// GD2: result of [`evaluate_for_derivatives`].
pub struct DerivativeXc {
    pub grid: XcGrid,
    /// The explicit gradient-metric stress of [`stress_metric`] (or
    /// [`stress_metric_pw92`]) of the same density, when requested and when
    /// the evaluation's gradients are those of the density itself (no point
    /// was screened); `None` otherwise (the caller computes it separately).
    pub metric: Option<[[f64; 3]; 3]>,
}

/// GD2: unpolarized PBE for the post-SCF derivatives, bit-identical to
/// [`evaluate`] (`prefactor` None: library correlation) or [`evaluate_pw92`]
/// (`Some`), with three cost changes: the gradients come from one forward
/// transform, the pointwise functional evaluation is split by `split` (each
/// participant evaluates its block; the per-point energy density, vrho and
/// vsigma sums are gathered; the energy is summed over the same values in
/// the same order), and with `metric` the GGA stress metric is formed from
/// the same gradients and vsigma instead of a second evaluation.
pub fn evaluate_for_derivatives(
    prefactor: Option<crate::pbe_pw92::Pw92Prefactor>,
    n: &Array3<f64>,
    cell: &Cell,
    fft: &mut impl Fft3d,
    split: Option<&dyn PointSplit>,
    metric: bool,
) -> Result<DerivativeXc, PhysicsError> {
    if n.is_empty() || n.shape() != fft.shape() {
        return Err(PhysicsError::Shape);
    }
    crate::xc::audit_density(n, cell.volume_bohr3(), "total", 0.)?;
    let screening = crate::xc::screen_negative();
    let projected = screening.then(|| n.mapv(crate::xc::screened_argument));
    let gradients = gradients_once(projected.as_ref().unwrap_or(n), cell, fft)?;
    // The metric uses the gradients of n itself: equal bits unless a point
    // was screened to zero.
    let same_gradients =
        projected.is_none() || n.iter().all(|&x| crate::xc::screened_argument(x) == x);
    let mut sigma = vec![0.0; n.len()];
    for g in &gradients {
        for (s, x) in sigma.iter_mut().zip(g) {
            *s += x * x;
        }
    }
    let rho: Vec<_> = n
        .iter()
        .map(|&x| if x <= DENSITY_THRESHOLD { 0.0 } else { x })
        .collect();
    for (s, r) in sigma.iter_mut().zip(&rho) {
        if *r == 0.0 {
            *s = 0.0;
        }
    }
    let len = n.len();
    let block = split.map_or(0..len, |s| s.block(len));
    // Per point of the block: n (x + c) exc, x + c vrho, x + c vsigma.
    let local = (|| -> Result<Vec<f64>, PhysicsError> {
        if block.is_empty() {
            return Ok(Vec::new());
        }
        let (r, g) = (&rho[block.clone()], &sigma[block.clone()]);
        let input = XcInput::gga(r, g);
        let eval = |id| {
            Functional::new(id, Spin::Unpolarized)
                .and_then(|f| f.eval(r.len(), &input))
                .map_err(|e| PhysicsError::Xc(e.to_string()))
        };
        let x = eval(FunctionalId::GgaXPbe)?;
        let c = match prefactor {
            None => {
                let c = eval(FunctionalId::GgaCPbe)?;
                CorrelationGrid {
                    exc: c.exc,
                    vrho: c.vrho,
                    vsigma: c.vsigma,
                }
            }
            Some(prefactor) => correlation_grid(r, g, prefactor),
        };
        let mut out = Vec::with_capacity(3 * r.len());
        for (i, n) in r.iter().enumerate() {
            out.push(n * (x.exc[i] + c.exc[i]));
            out.push(x.vrho[i] + c.vrho[i]);
            out.push(x.vsigma[i] + c.vsigma[i]);
        }
        Ok(out)
    })();
    let values = match split {
        Some(split) => {
            let failed = local.is_err();
            let gathered = split.gather(local.as_deref().unwrap_or(&[]), len, 3, failed);
            local?;
            gathered?
        }
        None => local?,
    };
    if values.len() != 3 * len {
        return Err(PhysicsError::Shape);
    }
    let weight = cell.volume_bohr3() / len as f64;
    // ndarray uses multiple accumulators, reducing long-grid summation error.
    let energy = weight * Array1::from_iter(values.iter().step_by(3).copied()).sum();
    let vsigma: Vec<f64> = values.iter().skip(2).step_by(3).copied().collect();
    let mut potential = Array3::zeros(n.raw_dim());
    for (v, &vr) in potential.iter_mut().zip(values.iter().skip(1).step_by(3)) {
        *v = vr;
    }
    let metric = (metric && same_gradients).then(|| {
        std::array::from_fn(|a| {
            std::array::from_fn(|b| {
                -2. / len as f64
                    * gradients[a]
                        .iter()
                        .zip(&gradients[b])
                        .enumerate()
                        .map(|(i, (a, b))| vsigma[i] * a * b)
                        .sum::<f64>()
            })
        })
    });
    for (axis, mut flux) in gradients.into_iter().enumerate() {
        for (f, s) in flux.iter_mut().zip(&vsigma) {
            *f *= 2.0 * s;
        }
        potential -= &derivative(&flux, cell, axis, fft)?;
    }
    if screening {
        for (v, &density) in potential.iter_mut().zip(n) {
            if density <= DENSITY_THRESHOLD {
                *v = 0.;
            }
        }
    }
    let integral = weight * (n * &potential).sum();
    if !energy.is_finite() || !integral.is_finite() || potential.iter().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(DerivativeXc {
        grid: XcGrid {
            potential,
            energy: Hartree(energy),
            potential_integral: Hartree(integral),
        },
        metric,
    })
}

/// Negative densities (including roundoff) and nonfinite inputs are rejected.
/// Energy density is zero at n<=threshold, with zero local derivatives. Gradients
/// use the original density; the divergence can be nonzero at screened points.
/// The derivative identity holds away from threshold crossings.
pub fn evaluate(
    n: &Array3<f64>,
    cell: &Cell,
    fft: &mut impl Fft3d,
) -> Result<XcGrid, PhysicsError> {
    if n.is_empty() || n.shape() != fft.shape() {
        return Err(PhysicsError::Shape);
    }
    crate::xc::audit_density(n, cell.volume_bohr3(), "total", 0.)?;
    let screening = crate::xc::screen_negative();
    let projected = screening.then(|| n.mapv(crate::xc::screened_argument));
    let gradients = (0..3)
        .map(|a| derivative(projected.as_ref().unwrap_or(n), cell, a, fft))
        .collect::<Result<Vec<_>, _>>()?;
    let mut sigma = vec![0.0; n.len()];
    for g in &gradients {
        for (s, x) in sigma.iter_mut().zip(g) {
            *s += x * x;
        }
    }
    let rho: Vec<_> = n
        .iter()
        .map(|&x| if x <= DENSITY_THRESHOLD { 0.0 } else { x })
        .collect();
    for (s, r) in sigma.iter_mut().zip(&rho) {
        if *r == 0.0 {
            *s = 0.0;
        }
    }
    let input = XcInput::gga(&rho, &sigma);
    let eval = |id| {
        Functional::new(id, Spin::Unpolarized)
            .and_then(|f| f.eval(n.len(), &input))
            .map_err(|e| PhysicsError::Xc(e.to_string()))
    };
    let x = eval(FunctionalId::GgaXPbe)?;
    let c = eval(FunctionalId::GgaCPbe)?;
    let weight = cell.volume_bohr3() / n.len() as f64;
    // ndarray uses multiple accumulators, reducing long-grid summation error.
    let energy = weight
        * Array1::from_iter(
            rho.iter()
                .zip(&x.exc)
                .zip(&c.exc)
                .map(|((n, x), c)| n * (x + c)),
        )
        .sum();
    let mut potential = Array3::zeros(n.raw_dim());
    for ((v, x), c) in potential.iter_mut().zip(&x.vrho).zip(&c.vrho) {
        *v = x + c;
    }
    for (axis, mut flux) in gradients.into_iter().enumerate() {
        for ((f, x), c) in flux.iter_mut().zip(&x.vsigma).zip(&c.vsigma) {
            *f *= 2.0 * (x + c);
        }
        potential -= &derivative(&flux, cell, axis, fft)?;
    }
    if screening {
        for (v, &density) in potential.iter_mut().zip(n) {
            if density <= DENSITY_THRESHOLD {
                *v = 0.;
            }
        }
    }
    let integral = weight * (n * &potential).sum();
    if !energy.is_finite() || !integral.is_finite() || potential.iter().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(XcGrid {
        potential,
        energy: Hartree(energy),
        potential_integral: Hartree(integral),
    })
}

/// Explicit gradient-metric contribution at fixed fractional density samples.
/// Density/volume response is accounted for separately through the XC potential.
pub fn stress_metric(
    n: &Array3<f64>,
    cell: &Cell,
    fft: &mut impl Fft3d,
) -> Result<[[f64; 3]; 3], PhysicsError> {
    let gradients = (0..3)
        .map(|a| derivative(n, cell, a, fft))
        .collect::<Result<Vec<_>, _>>()?;
    let rho: Vec<_> = n
        .iter()
        .map(|&v| if v <= DENSITY_THRESHOLD { 0. } else { v })
        .collect();
    let sigma: Vec<f64> = (0..n.len())
        .map(|i| {
            if rho[i] == 0. {
                0.
            } else {
                gradients
                    .iter()
                    .map(|g| g.as_slice().unwrap()[i].powi(2))
                    .sum()
            }
        })
        .collect();
    let input = XcInput::gga(&rho, &sigma);
    let eval = |id| {
        Functional::new(id, Spin::Unpolarized)
            .and_then(|f| f.eval(n.len(), &input))
            .map_err(|e| PhysicsError::Xc(e.to_string()))
    };
    let x = eval(FunctionalId::GgaXPbe)?;
    let c = eval(FunctionalId::GgaCPbe)?;
    Ok(std::array::from_fn(|a| {
        std::array::from_fn(|b| {
            -2. / n.len() as f64
                * gradients[a]
                    .iter()
                    .zip(&gradients[b])
                    .enumerate()
                    .map(|(i, (a, b))| (x.vsigma[i] + c.vsigma[i]) * a * b)
                    .sum::<f64>()
        })
    }))
}

// Separate opt-in assembly preserves every operation of the legacy path.
pub fn evaluate_pw92(
    prefactor: crate::pbe_pw92::Pw92Prefactor,
    n: &Array3<f64>,
    cell: &Cell,
    fft: &mut impl Fft3d,
) -> Result<XcGrid, PhysicsError> {
    if n.is_empty() || n.shape() != fft.shape() {
        return Err(PhysicsError::Shape);
    }
    crate::xc::audit_density(n, cell.volume_bohr3(), "total", 0.)?;
    let screening = crate::xc::screen_negative();
    let projected = screening.then(|| n.mapv(crate::xc::screened_argument));
    let gradients = (0..3)
        .map(|a| derivative(projected.as_ref().unwrap_or(n), cell, a, fft))
        .collect::<Result<Vec<_>, _>>()?;
    let mut sigma = vec![0.0; n.len()];
    for g in &gradients {
        for (s, x) in sigma.iter_mut().zip(g) {
            *s += x * x;
        }
    }
    let rho: Vec<_> = n
        .iter()
        .map(|&x| if x <= DENSITY_THRESHOLD { 0.0 } else { x })
        .collect();
    for (s, r) in sigma.iter_mut().zip(&rho) {
        if *r == 0.0 {
            *s = 0.0;
        }
    }
    let input = XcInput::gga(&rho, &sigma);
    let eval = |id| {
        Functional::new(id, Spin::Unpolarized)
            .and_then(|f| f.eval(n.len(), &input))
            .map_err(|e| PhysicsError::Xc(e.to_string()))
    };
    let x = eval(FunctionalId::GgaXPbe)?;
    let c = correlation_grid(&rho, &sigma, prefactor);
    let weight = cell.volume_bohr3() / n.len() as f64;
    // ndarray uses multiple accumulators, reducing long-grid summation error.
    let energy = weight
        * Array1::from_iter(
            rho.iter()
                .zip(&x.exc)
                .zip(&c.exc)
                .map(|((n, x), c)| n * (x + c)),
        )
        .sum();
    let mut potential = Array3::zeros(n.raw_dim());
    for ((v, x), c) in potential.iter_mut().zip(&x.vrho).zip(&c.vrho) {
        *v = x + c;
    }
    for (axis, mut flux) in gradients.into_iter().enumerate() {
        for ((f, x), c) in flux.iter_mut().zip(&x.vsigma).zip(&c.vsigma) {
            *f *= 2.0 * (x + c);
        }
        potential -= &derivative(&flux, cell, axis, fft)?;
    }
    if screening {
        for (v, &density) in potential.iter_mut().zip(n) {
            if density <= DENSITY_THRESHOLD {
                *v = 0.;
            }
        }
    }
    let integral = weight * (n * &potential).sum();
    if !energy.is_finite() || !integral.is_finite() || potential.iter().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(XcGrid {
        potential,
        energy: Hartree(energy),
        potential_integral: Hartree(integral),
    })
}

/// Explicit gradient-metric contribution at fixed fractional density samples.
/// Density/volume response is accounted for separately through the XC potential.
pub fn stress_metric_pw92(
    prefactor: crate::pbe_pw92::Pw92Prefactor,
    n: &Array3<f64>,
    cell: &Cell,
    fft: &mut impl Fft3d,
) -> Result<[[f64; 3]; 3], PhysicsError> {
    let gradients = (0..3)
        .map(|a| derivative(n, cell, a, fft))
        .collect::<Result<Vec<_>, _>>()?;
    let rho: Vec<_> = n
        .iter()
        .map(|&v| if v <= DENSITY_THRESHOLD { 0. } else { v })
        .collect();
    let sigma: Vec<f64> = (0..n.len())
        .map(|i| {
            if rho[i] == 0. {
                0.
            } else {
                gradients
                    .iter()
                    .map(|g| g.as_slice().unwrap()[i].powi(2))
                    .sum()
            }
        })
        .collect();
    let input = XcInput::gga(&rho, &sigma);
    let eval = |id| {
        Functional::new(id, Spin::Unpolarized)
            .and_then(|f| f.eval(n.len(), &input))
            .map_err(|e| PhysicsError::Xc(e.to_string()))
    };
    let x = eval(FunctionalId::GgaXPbe)?;
    let c = correlation_grid(&rho, &sigma, prefactor);
    Ok(std::array::from_fn(|a| {
        std::array::from_fn(|b| {
            -2. / n.len() as f64
                * gradients[a]
                    .iter()
                    .zip(&gradients[b])
                    .enumerate()
                    .map(|(i, (a, b))| (x.vsigma[i] + c.vsigma[i]) * a * b)
                    .sum::<f64>()
        })
    }))
}

struct CorrelationGrid {
    exc: Vec<f64>,
    vrho: Vec<f64>,
    vsigma: Vec<f64>,
}
fn correlation_grid(
    rho: &[f64],
    sigma: &[f64],
    prefactor: crate::pbe_pw92::Pw92Prefactor,
) -> CorrelationGrid {
    let mut c = CorrelationGrid {
        exc: Vec::with_capacity(rho.len()),
        vrho: Vec::with_capacity(rho.len()),
        vsigma: Vec::with_capacity(rho.len()),
    };
    for (&n, &s) in rho.iter().zip(sigma) {
        let point = if n == 0.0 {
            crate::pbe_pw92::Correlation {
                exc: 0.0,
                vrho: 0.0,
                vsigma: 0.0,
            }
        } else {
            crate::pbe_pw92::correlation(n, s, prefactor)
        };
        c.exc.push(point.exc);
        c.vrho.push(point.vrho);
        c.vsigma.push(point.vsigma);
    }
    c
}
