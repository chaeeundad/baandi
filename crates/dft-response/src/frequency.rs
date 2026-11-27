//! Imaginary-axis frequency grids (plan 2.4, decision D5): Gauss-Legendre on
//! t in (0, 1) mapped by `w = w0 t / (1 - t)`, weights for `int_0^inf dw`.
use crate::{ResponseError, input};

/// Gauss-Legendre nodes (ascending) and weights on (-1, 1), Newton iteration on
/// the three-term recurrence. Exact for polynomials of degree <= 2n - 1.
pub fn gauss_legendre(n: usize) -> Result<(Vec<f64>, Vec<f64>), ResponseError> {
    if n == 0 || n > 4096 {
        return Err(input("Gauss-Legendre order must be 1..=4096"));
    }
    let mut nodes = vec![0.; n];
    let mut weights = vec![0.; n];
    for i in 0..n.div_ceil(2) {
        let mut x = (std::f64::consts::PI * (i as f64 + 0.75) / (n as f64 + 0.5)).cos();
        let mut derivative = 1.;
        for _ in 0..100 {
            let (mut p0, mut p1) = (1., x);
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            if n == 1 {
                p1 = x;
                p0 = 1.;
            }
            derivative = n as f64 * (x * p1 - p0) / (x * x - 1.);
            let step = p1 / derivative;
            x -= step;
            if step.abs() <= 1e-16 * x.abs().max(1.) {
                break;
            }
        }
        let w = 2. / ((1. - x * x) * derivative * derivative);
        nodes[i] = -x;
        nodes[n - 1 - i] = x;
        weights[i] = w;
        weights[n - 1 - i] = w;
    }
    if n % 2 == 1 {
        nodes[n / 2] = 0.;
    }
    Ok((nodes, weights))
}

/// Integration grid on the positive imaginary axis (Ha), ascending.
#[derive(Clone, Debug, PartialEq)]
pub struct FrequencyGrid {
    /// Frequencies w (Ha), ascending, all positive.
    pub points: Vec<f64>,
    /// Weights of `int_0^inf f(w) dw ~ sum_i weights_i f(points_i)` (Ha).
    pub weights: Vec<f64>,
    /// Mapping scale w0 (Ha): half of the nodes lie below it.
    pub omega0: f64,
}
impl FrequencyGrid {
    /// `n` Gauss-Legendre nodes on t in (0, 1), `w = w0 t / (1 - t)`,
    /// weight `w_GL / 2 * w0 / (1 - t)^2`.
    pub fn mapped_gauss_legendre(n: usize, omega0: f64) -> Result<Self, ResponseError> {
        if !omega0.is_finite() || omega0 <= 0. {
            return Err(input("frequency scale omega0 must be positive"));
        }
        let (x, w) = gauss_legendre(n)?;
        let mut points = Vec::with_capacity(n);
        let mut weights = Vec::with_capacity(n);
        for (x, w) in x.iter().zip(&w) {
            let t = 0.5 * (x + 1.);
            points.push(omega0 * t / (1. - t));
            weights.push(0.5 * w * omega0 / ((1. - t) * (1. - t)));
        }
        Ok(Self {
            points,
            weights,
            omega0,
        })
    }
}

/// Rule of plan 2.4 (from the computed spectrum, never a material constant):
/// w0 = mean of the lowest direct transition and the largest transition
/// included.
pub fn auto_omega0(lowest_direct: f64, largest: f64) -> Result<f64, ResponseError> {
    if !(lowest_direct.is_finite() && largest.is_finite())
        || lowest_direct <= 0.
        || largest < lowest_direct
    {
        return Err(input("transition energies for omega0"));
    }
    Ok(0.5 * (lowest_direct + largest))
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn gauss_legendre_is_exact_for_polynomials() {
        for n in [1, 2, 5, 16, 33] {
            let (x, w) = gauss_legendre(n).unwrap();
            for p in 0..2 * n {
                let sum: f64 = x.iter().zip(&w).map(|(x, w)| w * x.powi(p as i32)).sum();
                let exact = if p % 2 == 1 { 0. } else { 2. / (p + 1) as f64 };
                assert!((sum - exact).abs() < 1e-13, "n={n} p={p} {sum} {exact}");
            }
            assert!(x.windows(2).all(|p| p[0] < p[1]));
        }
    }
    #[test]
    fn mapped_grid_integrates_lorentzians() {
        // int_0^inf a / (a^2 + w^2) dw = pi / 2 for every a > 0.
        let grid = FrequencyGrid::mapped_gauss_legendre(32, 0.5).unwrap();
        for a in [0.2, 0.5, 2.0] {
            let sum: f64 = grid
                .points
                .iter()
                .zip(&grid.weights)
                .map(|(w, g)| g * a / (a * a + w * w))
                .sum();
            assert!(
                (sum - std::f64::consts::FRAC_PI_2).abs() < 1e-8,
                "a={a} {sum}"
            );
        }
    }
}
