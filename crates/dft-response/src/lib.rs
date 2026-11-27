//! Linear density response of insulators on the imaginary frequency axis
//! (R2 of `docs/roadmap/rpa-crpa-plan.md`): pair densities, the independent-
//! particle response chi0(q, i omega) by the sum over states (decision D1),
//! head and wings at q -> 0 with the velocity operator d H_k / d k including
//! the nonlocal commutator, and the symmetrized dielectric matrix with the
//! macroscopic dielectric tensor with and without local fields.
//!
//! Numerical core only: no file I/O, serde, logging or transport. Hartree
//! atomic units throughout. Orchestration (NSCF bands, MPI distribution of
//! transitions) lives in `dft-engine::response`, writers in `dft-cli`.
//!
//! Conventions (plan 2.1, 2.2):
//! - pair density `rho_nm(q+G) = <psi_nk| e^{-i(q+G).r} |psi_m,k+q>`
//!   `= sum_G1 conj(c_nk(G1)) c_m,k+q(G1 + G)` (dimensionless, cell-normalized
//!   orbitals, `rho_nn(0) = 1`);
//! - `chi0_GG'(q, i w) = sum_k s_k sum_{n occ, m emp} rho_nm(q+G) conj(rho_nm(q+G'))
//!   2 D / (D^2 + w^2)`, `D = e_nk - e_m,k+q < 0`, `s_k = 2 w_k / Omega`
//!   (spin factor 2, k weights summing to one); valid for a time-reversal
//!   symmetric k set (the Adler-Wiser terms with n empty and m occupied are
//!   the time-reversed partners);
//! - q -> 0: `rho_nm(q) = q . a_nm + O(q^2)`, `a_nm = <u_n| dH_k/dk |u_m> / (e_m - e_n)`,
//!   head `chi0_00 = q^2 q^.L.q^`, wing `chi0_G0 = q q^.W_G`;
//! - symmetrized `eps~ = 1 - v^1/2 chi0 v^1/2`, `v(q+G) = 4 pi / |q+G|^2`.
pub mod chi0;
pub mod dielectric;
pub mod frequency;
pub mod pair;
pub mod transitions;
pub mod velocity;

pub use num_complex::Complex64;

#[derive(Debug, thiserror::Error)]
pub enum ResponseError {
    #[error("response input: {0}")]
    Input(String),
    #[error("response numerics: {0}")]
    Numerical(String),
}
pub(crate) fn numerical(e: impl std::fmt::Display) -> ResponseError {
    ResponseError::Numerical(e.to_string())
}
pub(crate) fn input(e: impl std::fmt::Display) -> ResponseError {
    ResponseError::Input(e.to_string())
}

/// Live column count of a scratch matrix, growing its capacity when needed
/// (reused across calls; no allocation once grown).
pub(crate) fn set_columns(
    m: &mut dft_numerics::linalg::ComplexMatrix,
    cols: usize,
) -> Result<(), ResponseError> {
    if cols > m.cols() {
        m.reserve_columns(cols).map_err(numerical)?;
    }
    m.resize_columns(cols).map_err(numerical)
}
