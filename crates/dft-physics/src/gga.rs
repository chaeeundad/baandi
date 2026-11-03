//! Spin-unpolarized PBE with a discrete variational spectral gradient.
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
    if let Some(&x) = n.iter().find(|x| !x.is_finite() || **x < 0.0) {
        return Err(PhysicsError::GgaDensity(x));
    }
    let gradients = (0..3)
        .map(|a| derivative(n, cell, a, fft))
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
