//! Symmetrized dielectric matrix `eps~_GG' = delta_GG' - v^1/2(q+G) chi0_GG' v^1/2(q+G')`,
//! `v^1/2(p) = sqrt(4 pi) / |p|`, and the macroscopic dielectric function
//! `eps_M = 1 / [eps~^-1]_00` (plan 2.2, 2.3).
//!
//! Optical limit q -> 0 along q^: with `chi0_00 = q^2 q^.L.q^` and
//! `chi0_G0 = q q^.W_G`, the head is `1 - 4 pi q^.L.q^` and the wings are
//! `-4 pi q^.W_G / |G|`; the Schur complement over the body B (G != 0) gives
//! `eps_M(q^) = q^.E.q^` with the tensor
//! `E_ab = delta_ab - 4 pi L_ab - (4 pi)^2 sum_GG' conj(W_Ga) [eps~_BB^-1]_GG' W_G'b / (|G| |G'|)`
//! (with local fields) and `E = 1 - 4 pi L` without local fields.
use crate::{ResponseError, chi0::Chi0, input, numerical};
use dft_numerics::linalg::{ComplexMatrix, hermitian_eigh};
use num_complex::Complex64 as C;
use std::f64::consts::PI;

/// Macroscopic tensors at one frequency.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OpticalEpsilon {
    /// Real symmetric part of E without local fields (head only).
    pub without_local_fields: [[f64; 3]; 3],
    /// Real symmetric part of E with local fields (Schur complement).
    pub with_local_fields: [[f64; 3]; 3],
    /// Largest |Im E| and |E - E^T| entries of either tensor (time reversal
    /// makes both zero up to roundoff).
    pub max_asymmetry: f64,
    /// Smallest eigenvalue of the body eps~_BB (>= 1 for chi0 <= 0).
    pub min_body_eigenvalue: f64,
}

fn real_symmetric(e: &[[C; 3]; 3]) -> ([[f64; 3]; 3], f64) {
    let mut out = [[0.; 3]; 3];
    let mut worst = 0_f64;
    for a in 0..3 {
        for b in 0..3 {
            out[a][b] = 0.5 * (e[a][b].re + e[b][a].re);
            worst = worst
                .max(e[a][b].im.abs())
                .max((e[a][b].re - e[b][a].re).abs());
        }
    }
    (out, worst)
}

/// Hermitian part of `delta - sv chi0 sv` restricted to `rows` (indices into
/// the response sphere), `sv[g] = sqrt(4 pi) / |q+G|`.
fn symmetrized(
    body: &ComplexMatrix,
    sv: &[f64],
    rows: &[usize],
) -> Result<ComplexMatrix, ResponseError> {
    let n = rows.len();
    let mut e = ComplexMatrix::zeros(n, n).map_err(numerical)?;
    for (j, &gj) in rows.iter().enumerate() {
        for (i, &gi) in rows.iter().enumerate() {
            let a = body[(gi, gj)];
            let b = body[(gj, gi)].conj();
            e[(i, j)] = C::new(f64::from(i == j), 0.) - 0.5 * (a + b) * sv[gi] * sv[gj];
        }
    }
    Ok(e)
}

/// Optical tensors at frequency index `w`. `g_norms[g] = |G|` of the response
/// sphere at q = 0 (entry 0 is G = 0). `chi0` must come from an optical
/// accumulation.
pub fn optical(chi0: &Chi0, w: usize, g_norms: &[f64]) -> Result<OpticalEpsilon, ResponseError> {
    let n = g_norms.len();
    if w >= chi0.body.len()
        || chi0.head.len() != chi0.body.len()
        || chi0.body[w].rows() != n
        || chi0.wings[w].cols() != 3
        || g_norms[0] != 0.
        || g_norms[1..].iter().any(|g| !g.is_finite() || *g <= 0.)
    {
        return Err(input("optical epsilon: shapes or norms"));
    }
    let l = &chi0.head[w];
    let mut e = [[C::default(); 3]; 3];
    for a in 0..3 {
        for b in 0..3 {
            e[a][b] = C::new(f64::from(a == b), 0.) - 4. * PI * l[a][b];
        }
    }
    let (nlf, asym_nlf) = real_symmetric(&e);
    let sv: Vec<f64> = g_norms
        .iter()
        .map(|g| if *g > 0. { (4. * PI).sqrt() / g } else { 0. })
        .collect();
    let rows: Vec<usize> = (1..n).collect();
    let mut min_body = f64::INFINITY;
    if !rows.is_empty() {
        let body = symmetrized(&chi0.body[w], &sv, &rows)?;
        let (lambda, v) = hermitian_eigh(&body).map_err(numerical)?;
        min_body = lambda.first().copied().unwrap_or(f64::INFINITY);
        if lambda.iter().any(|x| !x.is_finite() || *x <= 0.) {
            return Err(numerical(format!(
                "eps~ body not positive definite (smallest eigenvalue {min_body:e})"
            )));
        }
        // Y = V^H U, U_b[G] = sqrt(4 pi) W_Gb / |G| (so the (4 pi)^2 factor is in U).
        let wings = &chi0.wings[w];
        let m = rows.len();
        let mut y = vec![[C::default(); 3]; m];
        for (k, yk) in y.iter_mut().enumerate() {
            for (i, &g) in rows.iter().enumerate() {
                let vc = v[(i, k)].conj();
                for b in 0..3 {
                    yk[b] += vc * wings[(g, b)] * sv[g];
                }
            }
        }
        for (k, yk) in y.iter().enumerate() {
            for a in 0..3 {
                for b in 0..3 {
                    e[a][b] -= 4. * PI * yk[a].conj() * yk[b] / lambda[k];
                }
            }
        }
    }
    let (lf, asym_lf) = real_symmetric(&e);
    Ok(OpticalEpsilon {
        without_local_fields: nlf,
        with_local_fields: lf,
        max_asymmetry: asym_nlf.max(asym_lf),
        min_body_eigenvalue: min_body,
    })
}

/// Finite q (no G with q+G = 0): `(eps_M without local fields = eps~_00,
/// eps_M with local fields = 1 / [eps~^-1]_00, smallest eigenvalue of eps~)`.
/// `qg_norms[g] = |q+G|`, entry 0 the G = 0 term.
pub fn finite_q(body: &ComplexMatrix, qg_norms: &[f64]) -> Result<(f64, f64, f64), ResponseError> {
    let n = qg_norms.len();
    if body.rows() != n || body.cols() != n || qg_norms.iter().any(|g| !g.is_finite() || *g <= 0.) {
        return Err(input("finite-q epsilon: shapes or norms"));
    }
    let sv: Vec<f64> = qg_norms.iter().map(|g| (4. * PI).sqrt() / g).collect();
    let rows: Vec<usize> = (0..n).collect();
    let e = symmetrized(body, &sv, &rows)?;
    let (lambda, v) = hermitian_eigh(&e).map_err(numerical)?;
    if lambda.iter().any(|x| !x.is_finite() || *x <= 0.) {
        return Err(numerical("eps~ not positive definite"));
    }
    let inverse_00: f64 = (0..n).map(|k| v[(0, k)].norm_sqr() / lambda[k]).sum();
    Ok((e[(0, 0)].re, 1. / inverse_00, lambda[0]))
}
