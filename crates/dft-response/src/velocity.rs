//! Velocity matrix elements `p_nm = <u_n| dH_k/dk |u_m>` (Ha Bohr) of the
//! periodic parts in the plane-wave basis of one k (plan 2.2):
//! `dH_k/dk_b = (k+G)_b delta_GG' + d V_NL(k) / dk_b`, where
//! `V_NL = sum_atom B D B^H` and `dV_NL/dk_b = sum_atom (G_b D B^H + B D G_b^H)`
//! with `G_b = dB/dk_b` from `dft_physics::stress::projector_k_gradient`. The
//! nonlocal term is mandatory for NC pseudopotentials (the commutator
//! `i[V_NL, r]`).
use crate::{ResponseError, input, numerical};
use dft_numerics::linalg::{ComplexMatrix, MatrixOp, gemm};
use dft_physics::stress::ProjectorKGradient;

/// Atom-block-diagonal `D X` (X: projections, channels x columns).
fn couple(
    projectors: &ProjectorKGradient,
    x: &ComplexMatrix,
) -> Result<ComplexMatrix, ResponseError> {
    let mut y = ComplexMatrix::zeros(x.rows(), x.cols()).map_err(numerical)?;
    for col in 0..x.cols() {
        let mut offset = 0;
        for (&count, d) in projectors.channels.iter().zip(&projectors.couplings) {
            for i in 0..count {
                y[(offset + i, col)] = (0..count)
                    .map(|j| d[i * count + j] * x[(offset + j, col)])
                    .sum();
            }
            offset += count;
        }
    }
    Ok(y)
}

/// `[p_x, p_y, p_z]`, each `left.cols() x right.cols()`, `p_b[(n, m)] =
/// <left_n| dH/dk_b |right_m>`. `momenta[g]` is k+G (Cartesian, Bohr^-1) of
/// row g; `projectors` (rows in the same order) may be `None` for a purely
/// local Hamiltonian.
pub fn velocity(
    momenta: &[[f64; 3]],
    projectors: Option<&ProjectorKGradient>,
    left: &ComplexMatrix,
    right: &ComplexMatrix,
) -> Result<[ComplexMatrix; 3], ResponseError> {
    let rows = momenta.len();
    if left.rows() != rows || right.rows() != rows {
        return Err(input("velocity: orbital rows differ from the basis"));
    }
    let mut out = Vec::with_capacity(3);
    for b in 0..3 {
        let mut scaled = right.clone();
        for col in 0..scaled.cols() {
            for (z, k) in scaled.column_mut(col).iter_mut().zip(momenta) {
                *z *= k[b];
            }
        }
        out.push(gemm(left, MatrixOp::Adjoint, &scaled, MatrixOp::Identity).map_err(numerical)?);
    }
    if let Some(p) = projectors {
        if p.values.rows() != rows || p.gradients.iter().any(|g| g.rows() != rows) {
            return Err(input("velocity: projector rows differ from the basis"));
        }
        if p.values.cols() > 0 {
            let project = |a: &ComplexMatrix, x: &ComplexMatrix| {
                gemm(a, MatrixOp::Adjoint, x, MatrixOp::Identity).map_err(numerical)
            };
            let b_left = project(&p.values, left)?;
            let db_right = couple(p, &project(&p.values, right)?)?;
            for (b, o) in out.iter_mut().enumerate() {
                let g_left = project(&p.gradients[b], left)?;
                let dg_right = couple(p, &project(&p.gradients[b], right)?)?;
                let first = gemm(&g_left, MatrixOp::Adjoint, &db_right, MatrixOp::Identity)
                    .map_err(numerical)?;
                let second = gemm(&b_left, MatrixOp::Adjoint, &dg_right, MatrixOp::Identity)
                    .map_err(numerical)?;
                for ((z, a), c) in o
                    .as_mut_slice()
                    .iter_mut()
                    .zip(first.as_slice())
                    .zip(second.as_slice())
                {
                    *z += a + c;
                }
            }
        }
    }
    let [x, y, z]: [ComplexMatrix; 3] = out
        .try_into()
        .map_err(|_| numerical("velocity components"))?;
    for m in [&x, &y, &z] {
        m.validate_finite().map_err(numerical)?;
    }
    Ok([x, y, z])
}

/// `<x| V_NL |x>` per column (Ha) from the same projector table: the H
/// consistency check of the table used for the velocity.
pub fn nonlocal_expectations(
    projectors: &ProjectorKGradient,
    x: &ComplexMatrix,
) -> Result<Vec<f64>, ResponseError> {
    let b =
        gemm(&projectors.values, MatrixOp::Adjoint, x, MatrixOp::Identity).map_err(numerical)?;
    let db = couple(projectors, &b)?;
    Ok((0..x.cols())
        .map(|col| {
            b.column(col)
                .iter()
                .zip(db.column(col))
                .map(|(a, c)| (a.conj() * c).re)
                .sum()
        })
        .collect())
}
