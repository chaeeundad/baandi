//! B5: bordered (secular-equation) update of the reduced eigenproblem of the
//! serial real Davidson (frozen default).
//!
//! When exactly one column is appended to a basis whose reduced matrix T has the
//! eigendecomposition T = C diag(theta) C^T, the new reduced matrix in the basis
//! [C, e_last] is the arrowhead matrix [[diag(theta), z], [z^T, d]] with
//! z = C^T b (b = the old part of the new column, d its diagonal entry). Its
//! eigenvalues are the roots of the secular equation
//!   F(lambda) = lambda - d + sum_k z_k^2 / (theta_k - lambda) = 0,
//! one per interval between consecutive poles and one beyond each end.
//! Implementation (general literature: Golub 1973; Bunch, Nielsen, Sorensen
//! 1978; Gu and Eisenstat 1995; standard divide-and-conquer deflation):
//! * deflation of negligible z_k and of close pole pairs by a Givens rotation
//!   (tolerance 8 eps max(|theta|, |d|, ||z||)),
//! * roots found relative to the nearer pole (origin + tau, so differences to
//!   the poles keep full relative accuracy) by Newton steps safeguarded by
//!   bisection inside a shrinking bracket,
//! * eigenvectors from the recomputed weights z_hat (Loewner formula), which
//!   makes them numerically orthogonal,
//! * new coefficients C' = [[C_rot, 0], [0, 1]] Y by the fixed-order GEMM, then
//!   the sign convention of `symmetric_eigh` (largest-magnitude entry positive,
//!   first index wins ties).
//!
//! Any conditioning doubt returns `false` and the caller runs the full dense
//! solve: a root without convergence, a non-finite value, broken interlacing,
//! or recomputed weights that differ from z by more than 1e-9 ||z||.
use dft_numerics::{
    linalg::LinalgError,
    real::{RealGemmWorkspace, RealMatrix, RealOp},
};

/// Bordered updates allowed in a row before a full solve refreshes the base
/// decomposition (bounds accumulated rounding drift).
pub(crate) const MAX_STREAK: usize = 8;
const MAX_ROOT_STEPS: usize = 200;
const WEIGHT_AGREEMENT: f64 = 1e-9;

/// Base decomposition and scratch, sized for the subspace capacity.
pub(crate) struct Bordered {
    /// Order of the stored base decomposition (0 = none).
    pub(crate) order: usize,
    /// Consecutive bordered updates since the last full solve.
    pub(crate) streak: usize,
    values: Vec<f64>,
    /// Base coefficients, column-major with leading dimension `order`.
    base: Vec<f64>,
    /// Base coefficients with the deflation rotations applied.
    rotated: Vec<f64>,
    y: Vec<f64>,
    theta: Vec<f64>,
    z: Vec<f64>,
    z_hat: Vec<f64>,
    deflated: Vec<bool>,
    poles: Vec<usize>,
    origin: Vec<f64>,
    tau: Vec<f64>,
    order_of: Vec<(f64, usize)>,
}
impl Bordered {
    pub(crate) fn new(capacity: usize) -> Self {
        Self {
            order: 0,
            streak: 0,
            values: vec![0.; capacity],
            base: vec![0.; capacity * capacity],
            rotated: vec![0.; capacity * capacity],
            y: vec![0.; capacity * capacity],
            theta: vec![0.; capacity],
            z: vec![0.; capacity],
            z_hat: vec![0.; capacity],
            deflated: vec![false; capacity],
            poles: Vec::with_capacity(capacity),
            origin: vec![0.; capacity + 1],
            tau: vec![0.; capacity + 1],
            order_of: Vec::with_capacity(capacity),
        }
    }
    pub(crate) fn invalidate(&mut self) {
        self.order = 0;
        self.streak = 0;
    }
    /// Store the decomposition of an m x m reduced matrix as the next base.
    pub(crate) fn store(&mut self, values: &[f64], coefficients: &RealMatrix) {
        let m = values.len();
        self.values[..m].copy_from_slice(values);
        self.base[..m * m].copy_from_slice(&coefficients.as_slice()[..m * m]);
        self.order = m;
    }

    /// Try the bordered update for the m x m reduced matrix `t` (its leading
    /// (m-1) x (m-1) block must be the matrix of the stored base). On success
    /// writes `values[..m]` and `coefficients` (m x m) and returns true.
    pub(crate) fn update(
        &mut self,
        gemm: &mut RealGemmWorkspace,
        t: &RealMatrix,
        values: &mut [f64],
        coefficients: &mut RealMatrix,
    ) -> Result<bool, LinalgError> {
        let m = t.rows();
        let n0 = m - 1;
        if self.order != n0 || n0 == 0 {
            return Ok(false);
        }
        let column = &t.as_slice()[n0 * m..m * m];
        let (b, d) = (&column[..n0], column[n0]);
        // z = C^T b, theta = base values, rotated = base.
        let base = &self.base[..n0 * n0];
        for k in 0..n0 {
            self.z[k] = dot(&base[k * n0..(k + 1) * n0], b);
        }
        self.theta[..n0].copy_from_slice(&self.values[..n0]);
        self.rotated[..n0 * n0].copy_from_slice(base);
        let z = &mut self.z[..n0];
        let theta = &mut self.theta[..n0];
        if !z.iter().chain(theta.iter()).all(|x| x.is_finite()) || !d.is_finite() {
            return Ok(false);
        }
        let z_norm = z.iter().fold(0_f64, |a, &x| a.hypot(x));
        let scale = theta
            .iter()
            .fold(d.abs().max(z_norm), |a, &x| a.max(x.abs()));
        let tol = 8. * f64::EPSILON * scale;
        // Deflation: negligible weights, then close pole pairs (Givens).
        let deflated = &mut self.deflated[..n0];
        for k in 0..n0 {
            deflated[k] = z[k].abs() <= tol;
        }
        let mut previous: Option<usize> = None;
        for j in 0..n0 {
            if deflated[j] {
                continue;
            }
            if let Some(i) = previous {
                let r = z[i].hypot(z[j]);
                let (c, s) = (z[j] / r, z[i] / r);
                if (c * s * (theta[j] - theta[i])).abs() <= tol {
                    // u_i = c e_i - s e_j, u_j = s e_i + c e_j: z_i -> 0, z_j -> r.
                    let (ti, tj) = (theta[i], theta[j]);
                    theta[i] = c * c * ti + s * s * tj;
                    theta[j] = s * s * ti + c * c * tj;
                    z[i] = 0.;
                    z[j] = r;
                    deflated[i] = true;
                    let (lo, hi) = self.rotated.split_at_mut(j * n0);
                    let ci = &mut lo[i * n0..(i + 1) * n0];
                    let cj = &mut hi[..n0];
                    for (a, b) in ci.iter_mut().zip(cj.iter_mut()) {
                        let (x, y) = (*a, *b);
                        *a = c * x - s * y;
                        *b = s * x + c * y;
                    }
                }
            }
            previous = Some(j);
        }
        self.poles.clear();
        self.poles.extend((0..n0).filter(|&k| !deflated[k]));
        let poles = &self.poles;
        let kk = poles.len();
        for w in poles.windows(2) {
            if theta[w[0]] >= theta[w[1]] {
                return Ok(false);
            }
        }
        // Secular roots r = 0..=kk; root r lies between poles r-1 and r.
        let lower = poles.first().map_or(d, |&p| theta[p].min(d)) - z_norm - tol;
        let upper = poles.last().map_or(d, |&p| theta[p].max(d)) + z_norm + tol;
        let roots = Roots {
            d,
            theta,
            z,
            poles,
            lower,
            upper,
        };
        let (origin, tau) = (&mut self.origin[..=kk], &mut self.tau[..=kk]);
        let solved = roots.solve_range(0, origin, tau);
        if !solved {
            return Ok(false);
        }
        // Root minus pole, accurate relative to the nearer pole.
        let gap = |r: usize, k: usize, o: &[f64], t: &[f64]| (o[r] - theta[k]) + t[r];
        // Interlacing check.
        for r in 0..=kk {
            if r > 0 && gap(r, poles[r - 1], &self.origin, &self.tau) <= 0. {
                return Ok(false);
            }
            if r < kk && gap(r, poles[r], &self.origin, &self.tau) >= 0. {
                return Ok(false);
            }
        }
        // Loewner weights: z_hat_k^2 = -(l_k - p_k)(l_K - p_k) prod_{j != k, j < K}
        // (l_j - p_k) / (p_j - p_k), every factor positive by interlacing.
        let mut worst = 0_f64;
        let weights = Weights {
            theta,
            z,
            poles,
            origin: &self.origin[..=kk],
            tau: &self.tau[..=kk],
        };
        // Signed weights per pole number into the scratch `y`, then placed.
        let scratch = &mut self.y[..kk];
        let w_ok = weights.range(0, scratch);
        if !w_ok {
            return Ok(false);
        }
        for (q, &k) in poles.iter().enumerate() {
            let zh = self.y[q];
            worst = worst.max((zh - z[k]).abs());
            self.z_hat[k] = zh;
        }
        if worst > WEIGHT_AGREEMENT * z_norm {
            return Ok(false);
        }
        // Eigenpairs in ascending order: deflated (theta_k, e_k), secular roots.
        self.order_of.clear();
        for k in 0..n0 {
            if deflated[k] {
                self.order_of.push((theta[k], k));
            }
        }
        for r in 0..=kk {
            self.order_of
                .push((self.origin[r] + self.tau[r], n0 + 1 + r));
        }
        self.order_of
            .sort_by(|a, b| a.0.total_cmp(&b.0).then(a.1.cmp(&b.1)));
        for (col, &(value, _)) in self.order_of.iter().enumerate() {
            values[col] = value;
        }
        let vectors = Vectors {
            theta,
            poles,
            z_hat: &self.z_hat[..n0],
            origin: &self.origin[..=kk],
            tau: &self.tau[..=kk],
            order_of: &self.order_of,
            n0,
        };
        let y = &mut self.y[..m * m];
        vectors.columns(0, y);
        if !values[..m].iter().all(|v| v.is_finite()) || !y.iter().all(|v| v.is_finite()) {
            return Ok(false);
        }
        // C' = [[C_rot, 0], [0, 1]] Y.
        coefficients.reshape_within_storage(m, m)?;
        {
            gemm.multiply_slices_into(
                dft_numerics::GemmSite::Other,
                &self.rotated[..n0 * n0],
                n0,
                n0,
                n0,
                RealOp::Identity,
                y,
                m,
                m,
                coefficients.as_mut().submatrix_mut(0, 0, n0, m),
            )?;
        }
        let out = coefficients.as_mut_slice();
        for col in 0..m {
            out[col * m + n0] = y[col * m + n0];
            // Sign convention of symmetric_eigh.
            let column = &mut out[col * m..(col + 1) * m];
            let mut big = 0;
            for (i, v) in column.iter().enumerate() {
                if v.abs() > column[big].abs() {
                    big = i;
                }
            }
            if column[big] < 0. {
                for v in column.iter_mut() {
                    *v = -*v;
                }
            }
        }
        if !out[..m * m].iter().all(|v| v.is_finite()) {
            return Ok(false);
        }
        Ok(true)
    }
}

/// Root minus pole, accurate relative to the nearer pole.
fn gap_of(origin: &[f64], tau: &[f64], theta: &[f64], r: usize, k: usize) -> f64 {
    (origin[r] - theta[k]) + tau[r]
}
/// Pieces of a parallel range (fixed, so the split never depends on the thread
/// count; each item's arithmetic is serial and unchanged).
/// Secular roots: root r between poles r-1 and r, as origin + tau.
struct Roots<'a> {
    d: f64,
    theta: &'a [f64],
    z: &'a [f64],
    poles: &'a [usize],
    lower: f64,
    upper: f64,
}
impl Roots<'_> {
    /// Roots `first..first + origin.len()`; false if any root failed.
    fn solve_range(&self, first: usize, origin: &mut [f64], tau: &mut [f64]) -> bool {
        for (i, (o, t)) in origin.iter_mut().zip(tau.iter_mut()).enumerate() {
            match self.root(first + i) {
                Some((origin, tau)) => {
                    *o = origin;
                    *t = tau;
                }
                None => return false,
            }
        }
        true
    }
    /// Newton steps safeguarded by bisection inside a shrinking bracket.
    fn root(&self, r: usize) -> Option<(f64, f64)> {
        let (d, theta, z, poles) = (self.d, self.theta, self.z, self.poles);
        let kk = poles.len();
        if kk == 0 {
            return Some((d, 0.));
        }
        let pole = |r: usize| theta[poles[r]];
        let (origin, mut lo, mut hi) = if r == 0 {
            (pole(0), self.lower - pole(0), 0.)
        } else if r == kk {
            (pole(kk - 1), 0., self.upper - pole(kk - 1))
        } else {
            let half = 0.5 * (pole(r) - pole(r - 1));
            let (f, _, _) = secular(pole(r - 1), half, d, theta, z, poles);
            if f >= 0. {
                (pole(r - 1), 0., half)
            } else {
                (pole(r), -half, 0.)
            }
        };
        let mut tau = 0.5 * (lo + hi);
        for _ in 0..MAX_ROOT_STEPS {
            let (f, df, size) = secular(origin, tau, d, theta, z, poles);
            if !f.is_finite() || !df.is_finite() {
                return None;
            }
            if f.abs() <= 8. * f64::EPSILON * size {
                return Some((origin, tau));
            }
            if f < 0. {
                lo = tau;
            } else {
                hi = tau;
            }
            let newton = tau - f / df;
            let next = if newton > lo && newton < hi {
                newton
            } else {
                0.5 * (lo + hi)
            };
            if (hi - lo) <= 2. * f64::EPSILON * lo.abs().max(hi.abs()) || next == tau {
                return Some((origin, next));
            }
            tau = next;
        }
        None
    }
}
/// Loewner weights z_hat_k^2 = -(l_k - p_k)(l_K - p_k) prod_{j != k, j < K}
/// (l_j - p_k) / (p_j - p_k), every factor positive by interlacing.
struct Weights<'a> {
    theta: &'a [f64],
    z: &'a [f64],
    poles: &'a [usize],
    origin: &'a [f64],
    tau: &'a [f64],
}
impl Weights<'_> {
    /// Signed weights of pole numbers `first..first + out.len()`; false on a
    /// non-finite or negative square.
    fn range(&self, first: usize, out: &mut [f64]) -> bool {
        let kk = self.poles.len();
        let gap = |r: usize, k: usize| gap_of(self.origin, self.tau, self.theta, r, k);
        for (i, dst) in out.iter_mut().enumerate() {
            let q = first + i;
            let k = self.poles[q];
            let mut w = -gap(q, k) * gap(kk, k);
            for (j, &p) in self.poles.iter().enumerate() {
                if j != q {
                    w *= gap(j, k) / (self.theta[p] - self.theta[k]);
                }
            }
            if !(w.is_finite() && w >= 0.) {
                return false;
            }
            *dst = w.sqrt().copysign(self.z[k]);
        }
        true
    }
}
/// Eigenvectors of the arrowhead matrix in the ascending order of `order_of`.
struct Vectors<'a> {
    theta: &'a [f64],
    poles: &'a [usize],
    z_hat: &'a [f64],
    origin: &'a [f64],
    tau: &'a [f64],
    order_of: &'a [(f64, usize)],
    n0: usize,
}
impl Vectors<'_> {
    /// Columns `first..` of Y (n0 + 1 rows each) into `y`.
    fn columns(&self, first: usize, y: &mut [f64]) {
        let m = self.n0 + 1;

        for (i, dst) in y.chunks_exact_mut(m).enumerate() {
            let id = self.order_of[first + i].1;
            dst.fill(0.);
            if id < self.n0 {
                dst[id] = 1.;
            } else {
                let r = id - self.n0 - 1;
                let mut sum = 1.;
                for &k in self.poles {
                    let v = self.z_hat[k] / gap_of(self.origin, self.tau, self.theta, r, k);
                    dst[k] = v;
                    sum += v * v;
                }
                dst[self.n0] = 1.;
                let inv = 1. / sum.sqrt();
                for v in dst.iter_mut() {
                    *v *= inv;
                }
            }
        }
    }
}

/// F(origin + tau), F', and the magnitude of its terms (for the stop test).
fn secular(
    origin: f64,
    tau: f64,
    d: f64,
    theta: &[f64],
    z: &[f64],
    poles: &[usize],
) -> (f64, f64, f64) {
    let mut f = (origin - d) + tau;
    let mut size = (origin - d).abs() + tau.abs();
    let mut df = 1.;
    for &k in poles {
        let delta = (theta[k] - origin) - tau;
        let q = z[k] / delta;
        f += z[k] * q;
        size += (z[k] * q).abs();
        df += q * q;
    }
    (f, df, size)
}

fn dot(x: &[f64], y: &[f64]) -> f64 {
    let mut sums = [0.; 4];
    for (a, b) in x.chunks(4).zip(y.chunks(4)) {
        for k in 0..a.len() {
            sums[k] += a[k] * b[k];
        }
    }
    (sums[0] + sums[1]) + (sums[2] + sums[3])
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_numerics::{pool::WorkerPool, real::symmetric_eigh};

    fn lcg(state: &mut u64) -> f64 {
        *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
        ((*state >> 11) as f64) / (1_u64 << 53) as f64 - 0.5
    }

    /// Bordered result against the full solve on matrices with clusters,
    /// exact degeneracies, tiny and zero couplings.
    #[test]
    fn bordered_update_matches_full_solve() {
        let mut state = 99;
        let mut used = 0;
        for case in 0..40 {
            let m = [2, 3, 7, 20, 41, 96][case % 6];
            let n0 = m - 1;
            let mut t = RealMatrix::zeros(m, m).unwrap();
            for j in 0..n0 {
                for i in 0..=j {
                    let v = if i == j {
                        // Clusters and exact degeneracies in the base spectrum.
                        match case % 4 {
                            0 => 0.1 * i as f64,
                            1 => 0.3 * (i / 3) as f64,
                            2 => 1e-9 * i as f64 + (i / 5) as f64,
                            _ => lcg(&mut state),
                        }
                    } else if case % 4 == 3 {
                        0.1 * lcg(&mut state)
                    } else {
                        0.
                    };
                    t[(i, j)] = v;
                    t[(j, i)] = v;
                }
            }
            for i in 0..m {
                let v = match case % 5 {
                    0 => lcg(&mut state),
                    1 => 1e-12 * lcg(&mut state),
                    2 if i % 2 == 0 => 0.,
                    _ => 1e-3 * lcg(&mut state),
                };
                t[(i, n0)] = v;
                t[(n0, i)] = v;
            }
            let mut base_t = RealMatrix::zeros(n0, n0).unwrap();
            for j in 0..n0 {
                for i in 0..n0 {
                    base_t[(i, j)] = t[(i, j)];
                }
            }
            let (bv, bc) = symmetric_eigh(&base_t).unwrap();
            let mut b = Bordered::new(m);
            b.store(&bv, &bc);
            let mut gemm = RealGemmWorkspace::with_pool(m, m, WorkerPool::new(1).unwrap()).unwrap();
            let mut values = vec![0.; m];
            let mut c = RealMatrix::zeros(m, m).unwrap();
            if !b.update(&mut gemm, &t, &mut values, &mut c).unwrap() {
                continue;
            }
            used += 1;
            let (ev, _) = symmetric_eigh(&t).unwrap();
            let norm = t.as_slice().iter().fold(1_f64, |a, &x| a.max(x.abs()));
            for (a, e) in values.iter().zip(&ev) {
                assert!(
                    (a - e).abs() <= 1e-13 * norm * m as f64,
                    "case {case}: {a} vs {e}"
                );
            }
            for j in 0..m {
                // Eigen-residual and orthonormality.
                let mut r = 0_f64;
                for i in 0..m {
                    let tc: f64 = (0..m).map(|k| t[(i, k)] * c[(k, j)]).sum();
                    r = r.max((tc - values[j] * c[(i, j)]).abs());
                }
                assert!(r <= 1e-13 * norm * m as f64, "case {case}: residual {r}");
                for k in 0..=j {
                    let o: f64 = (0..m).map(|i| c[(i, j)] * c[(i, k)]).sum();
                    let e = if j == k { 1. } else { 0. };
                    assert!(
                        (o - e).abs() <= 1e-13 * m as f64,
                        "case {case}: overlap {o}"
                    );
                }
            }
        }
        assert!(used >= 30, "bordered path used in {used} of 40 cases");
    }
}
