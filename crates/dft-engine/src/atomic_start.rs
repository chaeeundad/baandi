//! Shared atomic-start factorization contract, moved unchanged from the legacy start.
use crate::model::{AtomicStartReport, ScfError};
use dft_physics::data::Complex64;

pub(crate) fn numerical(e: impl std::fmt::Display) -> ScfError {
    ScfError::Numerical(e.to_string())
}

/// Smallest Cholesky pivot (relative remaining norm of a unit column) the
/// block path accepts. CholQR squares the condition number, so it is reliable
/// only well above sqrt(machine epsilon) ~ 1.5e-8; below this guard the exact
/// two-pass Gram-Schmidt (which also owns the dependence rule) runs instead.
/// A numerical-method validity bound, independent of the system.
pub(crate) const BLOCK_PIVOT_GUARD: f64 = 1e-6;

/// Upper Cholesky factor R (G = R^H R, column-major n x n) of a Hermitian
/// matrix from its upper triangle; None for a pivot that is nonfinite or
/// below `guard`. Replicated: every rank factors the same bits.
pub(crate) fn cholesky_upper(g: &[Complex64], n: usize, guard: f64) -> Option<Vec<Complex64>> {
    let mut r = vec![Complex64::default(); n * n];
    for j in 0..n {
        let mut d = g[j * n + j].re;
        for k in 0..j {
            d -= r[j * n + k].norm_sqr();
        }
        if !d.is_finite() || d <= guard * guard {
            return None;
        }
        let rjj = d.sqrt();
        r[j * n + j] = Complex64::new(rjj, 0.0);
        for i in j + 1..n {
            // G_ji = sum_k conj(R_kj) R_ki
            let mut v = g[i * n + j];
            for k in 0..j {
                v -= r[j * n + k].conj() * r[i * n + k];
            }
            r[i * n + j] = v / rjj;
        }
    }
    Some(r)
}

/// Inverse of an upper triangular column-major matrix (back substitution).
pub(crate) fn upper_inverse(r: &[Complex64], n: usize) -> Vec<Complex64> {
    let mut x = vec![Complex64::default(); n * n];
    for j in 0..n {
        x[j * n + j] = Complex64::new(1.0, 0.0) / r[j * n + j];
        for i in (0..j).rev() {
            let mut s = Complex64::default();
            for k in i + 1..=j {
                s += r[k * n + i] * x[j * n + k];
            }
            x[j * n + i] = -s / r[i * n + i];
        }
    }
    x
}

// R-SCF-START: the owned-row builder reuses the legacy factor/guard contract.
// This changes no call or arithmetic in the legacy start.
impl AtomicStartReport {
    pub(crate) fn owned_qr_transform(
        gram: &dft_numerics::RealMatrix,
        first: bool,
    ) -> Result<Option<dft_numerics::RealMatrix>, ScfError> {
        let n = gram.cols();
        let scale: Vec<f64> = (0..n)
            .map(|j| if first { gram[(j, j)].sqrt() } else { 1. })
            .collect();
        if scale.iter().any(|s| !s.is_finite() || *s == 0.) {
            return Ok(None);
        }
        let mut g: Vec<_> = gram
            .as_slice()
            .iter()
            .map(|&v| Complex64::new(v, 0.))
            .collect();
        for j in 0..n {
            for i in 0..n {
                g[j * n + i] /= scale[i] * scale[j];
            }
        }
        let Some(r) = cholesky_upper(&g, n, if first { BLOCK_PIVOT_GUARD } else { 0.5 }) else {
            return Ok(None);
        };
        let mut t = upper_inverse(&r, n);
        for j in 0..n {
            for i in 0..n {
                t[j * n + i] /= scale[i];
            }
        }
        dft_numerics::RealMatrix::new(n, n, t.into_iter().map(|v| v.re).collect())
            .map(Some)
            .map_err(numerical)
    }
}
