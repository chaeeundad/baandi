//! RPA correlation energy from the adiabatic-connection fluctuation-dissipation
//! theorem on the imaginary axis (R3 of `docs/roadmap/rpa-crpa-plan.md`, 2.1,
//! 2.3-2.6; ADR 0055):
//!
//! `E_c = sum_q w_q (1 / 2 pi) int_0^inf dw  Tr[ln(1 - X(q, i w)) + X(q, i w)]`,
//! `X = v^1/2 chi0 v^1/2`, `v^1/2(p) = sqrt(4 pi) / |p|` (Hartree per cell, q
//! weights summing to one).
//!
//! chi0(i w) is Hermitian negative semidefinite (accumulated as `-(Q Q^H)`),
//! so `1 - X` is Hermitian positive definite with eigenvalues >= 1:
//! `ln det(1 - X) = 2 sum ln diag(L)` from its Cholesky factor. Failure of
//! the factorization is a hard error (it would mean chi0 lost its sign).
//!
//! q -> 0: with head `chi0_00 = q^2 d.L.d` and wings `chi0_G0 = q d.W_G`
//! (unit direction d), `X(d)` is finite and depends on d only:
//! `X_00 = 4 pi d.L.d`, `X_G0 = sqrt(4 pi) v^1/2(G) d.W_G`. With the body
//! `B = 1 - X_BB` (G, G' != 0) the Schur complement gives
//! `ln det(1 - X(d)) = ln det B + ln(1 - 4 pi d.L.d - d.Q.d)`,
//! `Q_ab = 4 pi sum conj(U_Ga) [B^-1]_GG' U_G'b`, `U_Gb = v^1/2(G) W_Gb`,
//! i.e. `ln det B + ln eps_M^LF(d)`, and `Tr X(d) = Tr X_BB + 4 pi d.L.d`.
//! The q = 0 term of the mesh sum is the average of this direction function
//! over a fixed Lebedev set (the mini-BZ average of the q -> 0 limit; the
//! |q| dependence inside the mini zone is the O(N_q^-1) error that the
//! k-mesh series measures).
use crate::{ResponseError, chi0::Chi0, input, numerical};
use dft_numerics::linalg::{ComplexMatrix, hermitian_cholesky, hermitian_eigh};
use num_complex::Complex64 as C;
use std::f64::consts::PI;

/// `Tr[ln(1 - X) + X]` of one finite-q matrix and its diagnostics.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TraceTerm {
    /// `ln det(1 - X) + Tr X` (dimensionless, <= 0).
    pub value: f64,
    /// `ln det(1 - X)` (>= 0).
    pub log_det: f64,
    /// `Tr X` (<= 0).
    pub trace: f64,
    /// Smallest squared Cholesky pivot of `1 - X` (> 0).
    pub min_pivot: f64,
    /// `max |X_ij - conj(X_ji)| / max |X_ij|` before Hermitization.
    pub hermiticity: f64,
}

/// `X = v^1/2 chi0 v^1/2` restricted to `rows` and Hermitized; returns the
/// relative Hermiticity defect too.
fn coulomb_scaled(
    body: &ComplexMatrix,
    sv: &[f64],
    rows: &[usize],
) -> Result<(ComplexMatrix, f64), ResponseError> {
    let n = rows.len();
    let mut x = ComplexMatrix::zeros(n, n).map_err(numerical)?;
    let (mut defect, mut scale) = (0_f64, 0_f64);
    for (j, &gj) in rows.iter().enumerate() {
        for (i, &gi) in rows.iter().enumerate() {
            let a = body[(gi, gj)] * sv[gi] * sv[gj];
            let b = body[(gj, gi)].conj() * sv[gi] * sv[gj];
            defect = defect.max((a - b).norm());
            scale = scale.max(a.norm());
            x[(i, j)] = 0.5 * (a + b);
        }
    }
    Ok((x, if scale > 0. { defect / scale } else { 0. }))
}

fn one_minus(x: &ComplexMatrix) -> Result<ComplexMatrix, ResponseError> {
    let mut m = x.clone();
    for z in m.as_mut_slice() {
        *z = -*z;
    }
    for i in 0..m.rows() {
        m[(i, i)] += 1.;
    }
    Ok(m)
}

fn log_det_term(x: &ComplexMatrix, hermiticity: f64) -> Result<TraceTerm, ResponseError> {
    let trace: f64 = (0..x.rows()).map(|i| x[(i, i)].re).sum();
    let m = one_minus(x)?;
    let factor = hermitian_cholesky(&m).map_err(numerical)?.ok_or_else(|| {
        numerical("1 - v^1/2 chi0 v^1/2 is not positive definite (chi0 lost its sign)")
    })?;
    let log_det = factor.log_det();
    Ok(TraceTerm {
        value: log_det + trace,
        log_det,
        trace,
        min_pivot: factor.min_pivot(),
        hermiticity,
    })
}

/// Finite q (no q + G = 0): `qg_norms[g] = |q + G|` of the response sphere.
pub fn finite_q_term(body: &ComplexMatrix, qg_norms: &[f64]) -> Result<TraceTerm, ResponseError> {
    let n = qg_norms.len();
    if body.rows() != n || body.cols() != n || qg_norms.iter().any(|g| !g.is_finite() || *g <= 0.) {
        return Err(input("ACFDT finite q: shapes or norms"));
    }
    let sv: Vec<f64> = qg_norms.iter().map(|g| (4. * PI).sqrt() / g).collect();
    let rows: Vec<usize> = (0..n).collect();
    let (x, h) = coulomb_scaled(body, &sv, &rows)?;
    log_det_term(&x, h)
}

/// Independent path for checks: `sum_i ln(1 - lambda_i) + lambda_i` over
/// the eigenvalues of the Hermitized X (EVD instead of Cholesky).
pub fn finite_q_term_eigen(body: &ComplexMatrix, qg_norms: &[f64]) -> Result<f64, ResponseError> {
    let n = qg_norms.len();
    if body.rows() != n || qg_norms.iter().any(|g| !g.is_finite() || *g <= 0.) {
        return Err(input("ACFDT finite q: shapes or norms"));
    }
    let sv: Vec<f64> = qg_norms.iter().map(|g| (4. * PI).sqrt() / g).collect();
    let rows: Vec<usize> = (0..n).collect();
    let (x, _) = coulomb_scaled(body, &sv, &rows)?;
    let (lambda, _) = hermitian_eigh(&x).map_err(numerical)?;
    if lambda.iter().any(|l| *l >= 1.) {
        return Err(numerical("X has an eigenvalue >= 1"));
    }
    Ok(lambda.iter().map(|l| (-l).ln_1p() + l).sum())
}

/// The q -> 0 term as a function of the direction d at one frequency.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OpticalTerm {
    /// `ln det B`, `B = 1 - X_BB` (G != 0 body).
    pub body_log_det: f64,
    /// `Tr X_BB`.
    pub body_trace: f64,
    /// Hermitian part of `4 pi L` (head).
    pub head: [[f64; 3]; 3],
    /// Hermitian part of Q (local-field correction of the head).
    pub local_fields: [[f64; 3]; 3],
    /// Smallest squared Cholesky pivot of B.
    pub min_pivot: f64,
    /// Relative Hermiticity defect of X_BB and of the head tensor.
    pub hermiticity: f64,
}
impl OpticalTerm {
    /// `eps_M^LF(d) = 1 - 4 pi d.L.d - d.Q.d` (must be > 1 on the imaginary axis).
    pub fn macroscopic(&self, d: [f64; 3]) -> f64 {
        1. - quadratic(&self.head, d) - quadratic(&self.local_fields, d)
    }
    /// `Tr[ln(1 - X(d)) + X(d)]` along the unit direction d.
    pub fn value(&self, d: [f64; 3]) -> Result<f64, ResponseError> {
        let e = self.macroscopic(d);
        if !(e.is_finite() && e > 0.) {
            return Err(numerical(format!(
                "q -> 0: eps_M along {d:?} is {e:e} (1 - X not positive definite)"
            )));
        }
        Ok(self.body_log_det + e.ln() + self.body_trace + quadratic(&self.head, d))
    }
    /// The q -> 0 term without the G = 0 row and column (diagnostic).
    pub fn body_only(&self) -> f64 {
        self.body_log_det + self.body_trace
    }
    /// Weighted direction average (weights summing to one).
    pub fn average(&self, directions: &[([f64; 3], f64)]) -> Result<f64, ResponseError> {
        let mut sum = 0.;
        for (d, w) in directions {
            sum += w * self.value(*d)?;
        }
        Ok(sum)
    }
}

fn quadratic(t: &[[f64; 3]; 3], d: [f64; 3]) -> f64 {
    (0..3)
        .map(|a| (0..3).map(|b| d[a] * t[a][b] * d[b]).sum::<f64>())
        .sum()
}

/// q -> 0 term at frequency index `w` of an optical accumulation;
/// `g_norms[g] = |G|` of the response sphere at q = 0 (entry 0 is G = 0).
pub fn optical_term(chi0: &Chi0, w: usize, g_norms: &[f64]) -> Result<OpticalTerm, ResponseError> {
    let n = g_norms.len();
    if w >= chi0.body.len()
        || chi0.head.len() != chi0.body.len()
        || chi0.body[w].rows() != n
        || chi0.wings[w].cols() != 3
        || n == 0
        || g_norms[0] != 0.
        || g_norms[1..].iter().any(|g| !g.is_finite() || *g <= 0.)
    {
        return Err(input("ACFDT q -> 0: shapes or norms"));
    }
    let l = &chi0.head[w];
    let mut head = [[0.; 3]; 3];
    let mut head_defect = 0_f64;
    let mut head_scale = 0_f64;
    for a in 0..3 {
        for b in 0..3 {
            head[a][b] = 4. * PI * 0.5 * (l[a][b] + l[b][a].conj()).re;
            head_defect = head_defect.max((l[a][b] - l[b][a].conj()).norm());
            head_scale = head_scale.max(l[a][b].norm());
        }
    }
    let sv: Vec<f64> = g_norms
        .iter()
        .map(|g| if *g > 0. { (4. * PI).sqrt() / g } else { 0. })
        .collect();
    let rows: Vec<usize> = (1..n).collect();
    let mut local_fields = [[0.; 3]; 3];
    let (mut body_log_det, mut body_trace, mut min_pivot, mut hermiticity) =
        (0., 0., f64::INFINITY, 0.);
    if !rows.is_empty() {
        let (x, h) = coulomb_scaled(&chi0.body[w], &sv, &rows)?;
        let term = log_det_term(&x, h)?;
        body_log_det = term.log_det;
        body_trace = term.trace;
        min_pivot = term.min_pivot;
        hermiticity = h;
        let b = one_minus(&x)?;
        let factor = hermitian_cholesky(&b)
            .map_err(numerical)?
            .ok_or_else(|| numerical("q -> 0 body not positive definite"))?;
        // U_Gb = sqrt(4 pi) v^1/2(G) W_Gb; Q = U^H B^-1 U.
        let wings = &chi0.wings[w];
        let m = rows.len();
        let mut u = ComplexMatrix::zeros(m, 3).map_err(numerical)?;
        for (i, &g) in rows.iter().enumerate() {
            for c in 0..3 {
                u[(i, c)] = wings[(g, c)] * sv[g] * (4. * PI).sqrt();
            }
        }
        let mut solved = u.clone();
        factor.solve_in_place(&mut solved).map_err(numerical)?;
        let mut q = [[C::default(); 3]; 3];
        for (a, row) in q.iter_mut().enumerate() {
            for (c, z) in row.iter_mut().enumerate() {
                *z = u
                    .column(a)
                    .iter()
                    .zip(solved.column(c))
                    .map(|(x, y)| x.conj() * y)
                    .sum();
            }
        }
        for a in 0..3 {
            for c in 0..3 {
                local_fields[a][c] = 0.5 * (q[a][c] + q[c][a].conj()).re;
            }
        }
    }
    if head_scale > 0. {
        hermiticity = f64::max(hermiticity, head_defect / head_scale);
    }
    Ok(OpticalTerm {
        body_log_det,
        body_trace,
        head,
        local_fields,
        min_pivot,
        hermiticity,
    })
}

/// Lebedev-Laikov direction sets (orders 6, 14, 26, 50; exact for spherical
/// polynomials of degree 3, 5, 7, 11). Weights sum to one. The values are the
/// closed-form points and weights of the published rules (V. I. Lebedev,
/// 1975-1976); `lebedev_rules_are_exact` checks them against the exact
/// sphere moments.
pub fn lebedev(order: usize) -> Result<Vec<([f64; 3], f64)>, ResponseError> {
    let mut out = Vec::with_capacity(order);
    let a1 = |out: &mut Vec<([f64; 3], f64)>, w: f64| {
        for i in 0..3 {
            for s in [1., -1.] {
                let mut d = [0.; 3];
                d[i] = s;
                out.push((d, w));
            }
        }
    };
    let a2 = |out: &mut Vec<([f64; 3], f64)>, w: f64| {
        let r = std::f64::consts::FRAC_1_SQRT_2;
        for (i, j) in [(0, 1), (0, 2), (1, 2)] {
            for si in [1., -1.] {
                for sj in [1., -1.] {
                    let mut d = [0.; 3];
                    d[i] = si * r;
                    d[j] = sj * r;
                    out.push((d, w));
                }
            }
        }
    };
    let a3 = |out: &mut Vec<([f64; 3], f64)>, w: f64| {
        let r = 1. / 3_f64.sqrt();
        for sx in [1., -1.] {
            for sy in [1., -1.] {
                for sz in [1., -1.] {
                    out.push(([sx * r, sy * r, sz * r], w));
                }
            }
        }
    };
    match order {
        6 => a1(&mut out, 1. / 6.),
        14 => {
            a1(&mut out, 1. / 15.);
            a3(&mut out, 3. / 40.);
        }
        26 => {
            a1(&mut out, 1. / 21.);
            a2(&mut out, 4. / 105.);
            a3(&mut out, 9. / 280.);
        }
        50 => {
            a1(&mut out, 4. / 315.);
            a2(&mut out, 64. / 2835.);
            a3(&mut out, 27. / 1280.);
            // (l, l, m) with l = 1/sqrt(11), m = 3/sqrt(11) and permutations.
            let (l, m) = (1. / 11_f64.sqrt(), 3. / 11_f64.sqrt());
            for axis in 0..3 {
                for s0 in [1., -1.] {
                    for s1 in [1., -1.] {
                        for s2 in [1., -1.] {
                            let mut d = [l, l, l];
                            d[axis] = m;
                            out.push(([s0 * d[0], s1 * d[1], s2 * d[2]], 14641. / 725760.));
                        }
                    }
                }
            }
        }
        _ => {
            return Err(input(
                "q0_directions must be one of 6, 14, 26, 50 (Lebedev orders)",
            ));
        }
    }
    Ok(out)
}

/// Product rule (Gauss-Legendre in cos theta x uniform phi) with
/// `n x 2n` points: a high-order reference for the direction convergence.
pub fn product_directions(n: usize) -> Result<Vec<([f64; 3], f64)>, ResponseError> {
    let (x, w) = crate::frequency::gauss_legendre(n)?;
    let nphi = 2 * n;
    let mut out = Vec::with_capacity(n * nphi);
    for (c, wc) in x.iter().zip(&w) {
        let s = (1. - c * c).max(0.).sqrt();
        for j in 0..nphi {
            let phi = 2. * PI * (j as f64 + 0.5) / nphi as f64;
            out.push(([s * phi.cos(), s * phi.sin(), *c], wc / (2. * nphi as f64)));
        }
    }
    Ok(out)
}

/// `(1 / 2 pi) sum_w weight_w term_w` (Ha when the weights are in Ha).
pub fn frequency_integral(weights: &[f64], terms: &[f64]) -> Result<f64, ResponseError> {
    if weights.len() != terms.len() {
        return Err(input("frequency integral: lengths"));
    }
    Ok(weights.iter().zip(terms).map(|(w, t)| w * t).sum::<f64>() / (2. * PI))
}

/// Linear least-squares fit `E(c) = E_inf + A c^{-3/2}` over a cutoff series
/// (plan 2.5, Harl and Kresse).
#[derive(Clone, Debug, PartialEq)]
pub struct Extrapolation {
    /// Extrapolated value (same unit as the input values).
    pub value: f64,
    /// Slope A (value unit x Ha^{3/2}).
    pub slope: f64,
    /// `E_i - fit(c_i)` per input point.
    pub residuals: Vec<f64>,
    /// Largest |residual| among the three largest cutoffs (plan 2.5 linearity
    /// check); `None` with fewer than three points.
    pub max_residual_last3: Option<f64>,
}
pub fn extrapolate(cutoffs: &[f64], values: &[f64]) -> Result<Extrapolation, ResponseError> {
    let n = cutoffs.len();
    if n < 2
        || values.len() != n
        || cutoffs.iter().any(|c| !c.is_finite() || *c <= 0.)
        || values.iter().any(|v| !v.is_finite())
    {
        return Err(input(
            "extrapolation needs at least two positive cutoffs and finite values",
        ));
    }
    let x: Vec<f64> = cutoffs.iter().map(|c| c.powf(-1.5)).collect();
    let xm = x.iter().sum::<f64>() / n as f64;
    let ym = values.iter().sum::<f64>() / n as f64;
    let sxx: f64 = x.iter().map(|x| (x - xm).powi(2)).sum();
    if sxx <= 0. {
        return Err(input("extrapolation needs distinct cutoffs"));
    }
    let sxy: f64 = x.iter().zip(values).map(|(x, y)| (x - xm) * (y - ym)).sum();
    let slope = sxy / sxx;
    let value = ym - slope * xm;
    let residuals: Vec<f64> = x
        .iter()
        .zip(values)
        .map(|(x, y)| y - (value + slope * x))
        .collect();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|a, b| cutoffs[*b].total_cmp(&cutoffs[*a]));
    let max_residual_last3 = (n >= 3).then(|| {
        order[..3]
            .iter()
            .map(|&i| residuals[i].abs())
            .fold(0., f64::max)
    });
    Ok(Extrapolation {
        value,
        slope,
        residuals,
        max_residual_last3,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn double_factorial(n: i32) -> f64 {
        if n <= 0 {
            1.
        } else {
            n as f64 * double_factorial(n - 2)
        }
    }
    /// Exact sphere average of x^a y^b z^c.
    fn moment(p: [i32; 3]) -> f64 {
        if p.iter().any(|e| e % 2 == 1) {
            return 0.;
        }
        p.iter().map(|e| double_factorial(e - 1)).product::<f64>()
            / double_factorial(p.iter().sum::<i32>() + 1)
    }

    #[test]
    fn lebedev_rules_are_exact() {
        for (order, degree) in [(6, 3), (14, 5), (26, 7), (50, 11)] {
            let rule = lebedev(order).unwrap();
            assert_eq!(rule.len(), order);
            for (d, _) in &rule {
                assert!((d.iter().map(|x| x * x).sum::<f64>() - 1.).abs() < 1e-15);
            }
            for a in 0..=degree {
                for b in 0..=degree - a {
                    for c in 0..=degree - a - b {
                        let sum: f64 = rule
                            .iter()
                            .map(|(d, w)| w * d[0].powi(a) * d[1].powi(b) * d[2].powi(c))
                            .sum();
                        let exact = moment([a, b, c]);
                        assert!(
                            (sum - exact).abs() < 1e-14,
                            "order {order} x^{a} y^{b} z^{c}: {sum} vs {exact}"
                        );
                    }
                }
            }
        }
        assert!(lebedev(7).is_err());
        let product = product_directions(12).unwrap();
        let sum: f64 = product.iter().map(|(d, w)| w * d[0].powi(4)).sum();
        assert!((sum - 0.2).abs() < 1e-14);
    }

    #[test]
    fn extrapolation_recovers_exact_lines() {
        let c = [2., 3., 4.5, 6.];
        let y: Vec<f64> = c.iter().map(|c: &f64| -0.3 + 0.7 * c.powf(-1.5)).collect();
        let fit = extrapolate(&c, &y).unwrap();
        assert!((fit.value + 0.3).abs() < 1e-14 && (fit.slope - 0.7).abs() < 1e-13);
        assert!(fit.max_residual_last3.unwrap() < 1e-14);
        assert!(extrapolate(&[2.], &[1.]).is_err());
        assert!(extrapolate(&[2., 2.], &[1., 1.]).is_err());
    }

    fn random_nsd(n: usize, rank: usize, seed: u64) -> ComplexMatrix {
        let mut state = seed;
        let mut next = || {
            state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((state >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        };
        let mut p = ComplexMatrix::zeros(n, rank).unwrap();
        for z in p.as_mut_slice() {
            *z = C::new(next(), next()) * 0.4;
        }
        let mut m = ComplexMatrix::zeros(n, n).unwrap();
        for i in 0..n {
            for j in 0..n {
                m[(i, j)] = -(0..rank).map(|r| p[(i, r)] * p[(j, r)].conj()).sum::<C>();
            }
        }
        m
    }

    #[test]
    fn cholesky_log_det_equals_eigenvalue_path() {
        for (n, rank, seed) in [(1, 1, 1), (5, 3, 2), (17, 40, 3), (40, 12, 4)] {
            let body = random_nsd(n, rank, seed);
            let norms: Vec<f64> = (0..n).map(|i| 0.7 + 0.31 * i as f64).collect();
            let chol = finite_q_term(&body, &norms).unwrap();
            let eig = finite_q_term_eigen(&body, &norms).unwrap();
            println!("n={n}: cholesky {:.15e} eigen {eig:.15e}", chol.value);
            assert!(chol.value <= 0. && chol.log_det >= 0. && chol.trace <= 0.);
            assert!(chol.min_pivot >= 1. - 1e-12, "{}", chol.min_pivot);
            assert!(
                (chol.value - eig).abs() <= 1e-12 * chol.log_det.abs().max(1.),
                "n={n}: {} vs {eig}",
                chol.value
            );
            assert!(chol.hermiticity < 1e-15);
        }
        // A matrix with chi0 of the wrong sign rejects.
        let mut bad = random_nsd(4, 4, 9);
        for z in bad.as_mut_slice() {
            *z *= -40.;
        }
        assert!(finite_q_term(&bad, &[0.5; 4]).is_err());
    }

    #[test]
    fn optical_schur_complement_equals_full_determinant() {
        // Optical accumulation with random pairs and slopes: the Schur form
        // equals the direct log det of the full (G = 0 included) X(d).
        use crate::chi0::Accumulator;
        let n = 9;
        let t = 14;
        let body = random_nsd(n, t, 11);
        let mut acc = Accumulator::new(n, &[0.3], true, 1).unwrap();
        // Recover a factor P from a random NSD matrix is not needed: use random
        // pairs directly.
        let mut pairs = ComplexMatrix::zeros(n, t).unwrap();
        let mut slopes = Vec::new();
        let mut delta = Vec::new();
        for j in 0..t {
            for g in 0..n {
                pairs[(g, j)] = body[(g, j % n)] + C::new(0.01 * j as f64, -0.02 * g as f64);
            }
            slopes.push([
                C::new(0.3, 0.1 * j as f64),
                C::new(-0.2, 0.05),
                C::new(0.1 * (j % 3) as f64, -0.4),
            ]);
            delta.push(-0.2 - 0.07 * j as f64);
        }
        acc.add(&pairs, &delta, 0.05, Some(&slopes)).unwrap();
        let chi0 = acc.reduce(&dft_mpi::Comm::Single).unwrap();
        let norms: Vec<f64> = (0..n).map(|i| 0.9 * i as f64).collect();
        let term = optical_term(&chi0, 0, &norms).unwrap();
        for d in [
            [1., 0., 0.],
            [0., 0.6, 0.8],
            [1. / 3_f64.sqrt(), -1. / 3_f64.sqrt(), 1. / 3_f64.sqrt()],
        ] {
            // Full X(d) with the G = 0 row from head and wings.
            let l = &chi0.head[0];
            let wings = &chi0.wings[0];
            let mut x = ComplexMatrix::zeros(n, n).unwrap();
            let sv = |g: usize| (4. * PI).sqrt() / norms[g];
            for i in 1..n {
                for j in 1..n {
                    x[(i, j)] =
                        0.5 * (chi0.body[0][(i, j)] + chi0.body[0][(j, i)].conj()) * sv(i) * sv(j);
                }
                let wing: C = (0..3).map(|b| d[b] * wings[(i, b)]).sum();
                x[(i, 0)] = (4. * PI).sqrt() * sv(i) * wing;
                x[(0, i)] = x[(i, 0)].conj();
            }
            x[(0, 0)] = C::new(
                4. * PI
                    * (0..3)
                        .map(|a| (0..3).map(|b| d[a] * l[a][b].re * d[b]).sum::<f64>())
                        .sum::<f64>(),
                0.,
            );
            let (lambda, _) = hermitian_eigh(&x).unwrap();
            let direct: f64 = lambda.iter().map(|l| (-l).ln_1p() + l).sum();
            let schur = term.value(d).unwrap();
            println!("d={d:?}: schur {schur:.15e} direct {direct:.15e}");
            assert!((schur - direct).abs() < 1e-12 * direct.abs().max(1.));
            assert!(term.macroscopic(d) > 1.);
        }
    }
}
