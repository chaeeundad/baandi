//! Track B Gate 4 Benchmark Plan — Suite C: Solver Acceleration R&D.
//!
//! Isolated standalone benchmark comparing two targeted upgrades against the
//! Gate 2 baselines (`DirectPcg` and `FisicaroPcg`) on sharp VASPsol-erfc water
//! cavities (`32^3, 48^3, 64^3`, optional `--96`, `L = 16 Bohr`, `eps_b = 78.4`,
//! `kappa in {0.0, 0.15}`):
//!
//! 1. **Experiment C1 — Consistent Fisicaro Polarization RHS**:
//!    - `(a) DirectPcg`: exact discrete SPD baseline (`A_disc phi_pol = b_pol_disc`).
//!    - `(b) FisicaroPcg (Gate 2 b_pol_disc)`: Gate 2 RHS mixing `A_Fis` on LHS
//!      with `A_disc` on RHS.
//!    - `(c) FisicaroPcg (Total 4pi*rho_sol)`: consistent total-potential solve
//!      `A_Fis phi_solv = 4*pi*rho_sol`, then `phi_pol = phi_solv - phi_vac`,
//!      eliminating `(A_Fis - A_disc) phi_vac` identically (and saving the 4 FFTs
//!      used to form `b_pol_disc`).
//!    - `(d) FisicaroPcg (Consistent b_pol_Fis)`: consistent polarization-potential
//!      solve `A_Fis phi_pol = 4*pi*rho_sol - A_Fis phi_vac`.
//!
//! 2. **Experiment C2 — Inexact Fisicaro Preconditioner inside Flexible Krylov**:
//!    Targets the **exact discrete SPD operator `A_disc`** (`|E - E_direct| < 1e-11 Ha`)
//!    using `ws.apply_operator(...)` in the outer loop and `ws.solve_pde(..., FisicaroPcg,
//!    max_iterations = inner_k)` (with unit-norm RHS scaling and diagonal fallback) as
//!    an inexact preconditioner:
//!    - `FPCG-PR` (Polak-Ribiere flexible PCG with exact `A_disc`-norm line search)
//!    - `FCG(4)` (exact `A_disc`-orthogonal truncated GCR / Flexible CG over the last 4
//!      directions `{p_j, A p_j}`, guaranteeing monotone energy-error decrease even for
//!      variable or coarse-grid non-SPD preconditioners)
//!      for `inner_k in {4, 8}`, across cold (`tol = 1e-10`) and `warm(1e-3)@1e-5` modes.
//!
//! Run:
//! `cargo run --release -p dft-physics --example solvation_suite_c`
//! `cargo run --release -p dft-physics --example solvation_suite_c -- --repeats 3`
//! `cargo run --release -p dft-physics --example solvation_suite_c -- --grids 32,48 --repeats 1`
//! `cargo run --release -p dft-physics --example solvation_suite_c -- --96`

use std::f64::consts::PI;
use std::time::Instant;

use dft_geometry::Cell;
use dft_mpi::Collective;
use dft_physics::solvation::{
    IonicPseudoCharge, SolvationError, SolvationParams, SolvationWorkspace, SolveReport,
    SolverConfig, SolverKind, build_ionic_pseudo_charge_g, evaluate_cavity_fields,
    skew_adjoint_g_cart,
};
use num_complex::Complex64;

const FOUR_PI: f64 = 4.0 * PI;
const BOX_LENGTH_BOHR: f64 = 16.0;
const EPS_BULK: f64 = 78.4;
const COLD_TOL: f64 = 1.0e-10;
const WARM_TOL: f64 = 1.0e-5;
const INNER_TOL: f64 = 1.0e-2;
const MAX_DIRECT_ITERS: usize = 300;
const MAX_FLEX_OUTER_ITERS: usize = 80;
const FCG_WINDOW: usize = 4;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum FlexUpdate {
    PolakRibiere,
    FcgWindow(usize),
}

impl FlexUpdate {
    fn label(self) -> String {
        match self {
            Self::PolakRibiere => "FPCG-PR".to_string(),
            Self::FcgWindow(m) => format!("FCG({m})"),
        }
    }
}

struct GridSpectralMeta {
    volume: f64,
    total_points: usize,
    g2: Vec<f64>,
    zero_g2_indices: Vec<usize>,
}

impl GridSpectralMeta {
    fn new(cell: &Cell, dims: [usize; 3], local_positions: &[usize]) -> Self {
        let recip_cols = cell.reciprocal_columns_inverse_bohr();
        let mut g2 = Vec::with_capacity(local_positions.len());
        let mut zero_g2_indices = Vec::new();
        for (k, &g_idx) in local_positions.iter().enumerate() {
            let gc = skew_adjoint_g_cart(dims, g_idx, recip_cols);
            let norm2 = gc[0].mul_add(gc[0], gc[1].mul_add(gc[1], gc[2] * gc[2]));
            g2.push(norm2);
            if norm2 == 0.0 {
                zero_g2_indices.push(k);
            }
        }
        Self {
            volume: cell.volume_bohr3(),
            total_points: dims[0] * dims[1] * dims[2],
            g2,
            zero_g2_indices,
        }
    }

    fn zero_nullspace(&self, has_kappa: bool, vec_g: &mut [Complex64]) {
        if !has_kappa {
            for &k in &self.zero_g2_indices {
                vec_g[k] = Complex64::default();
            }
        }
    }

    fn inner<K: Collective>(&self, comm: &K, a: &[Complex64], b: &[Complex64]) -> f64 {
        let mut loc = 0.0_f64;
        for (&x, &y) in a.iter().zip(b) {
            loc += x.re * y.re + x.im * y.im;
        }
        let mut glob = [0.0_f64];
        comm.tree_sum_f64(&[loc], &mut glob).unwrap();
        self.volume * glob[0]
    }

    fn norm<K: Collective>(&self, comm: &K, a: &[Complex64]) -> f64 {
        let loc: f64 = a.iter().map(Complex64::norm_sqr).sum();
        let mut glob = [0.0_f64];
        comm.tree_sum_f64(&[loc], &mut glob).unwrap();
        (self.volume * glob[0]).sqrt()
    }

    fn electrostatic_energy_ha<K: Collective>(
        &self,
        comm: &K,
        rho_sol_g: &[Complex64],
        phi_pol_g: &[Complex64],
    ) -> f64 {
        0.5 * self.inner(comm, rho_sol_g, phi_pol_g)
    }
}

struct PreparedProblem {
    eps: Vec<f64>,
    kappa2: Vec<f64>,
    sqrt_eps: Vec<f64>,
    has_kappa: bool,
    eps_bar: f64,
    eps_mean: f64,
    kappa2_bar: f64,
    rho_sol_g: Vec<Complex64>,
    rhs_total_g: Vec<Complex64>,
    phi_vac_g: Vec<Complex64>,
    b_pol_disc_g: Vec<Complex64>,
    b_pol_fis_g: Vec<Complex64>,
}

#[derive(Clone, Debug)]
struct FlexSolveStats {
    converged: bool,
    outer_iters: usize,
    inner_iters: usize,
    outer_ffts: usize,
    inner_ffts: usize,
    fallbacks: usize,
    rel_residual: f64,
}

impl FlexSolveStats {
    fn total_ffts(&self) -> usize {
        self.outer_ffts + self.inner_ffts
    }
}

struct WarmState {
    phi_g: Vec<Complex64>,
    rhs_g: Vec<Complex64>,
}

fn cubic_cell(length_bohr: f64) -> Cell {
    Cell::from_user_rows_bohr([
        [length_bohr, 0.0, 0.0],
        [0.0, length_bohr, 0.0],
        [0.0, 0.0, length_bohr],
    ])
    .unwrap()
}

fn water_ions() -> [IonicPseudoCharge; 2] {
    [
        IonicPseudoCharge {
            charge: 2.0,
            fractional_position: [0.43, 0.49, 0.50],
            width_bohr: 0.60,
        },
        IonicPseudoCharge {
            charge: 2.0,
            fractional_position: [0.57, 0.51, 0.50],
            width_bohr: 0.60,
        },
    ]
}

fn vaspsol_params(kappa: f64) -> SolvationParams {
    SolvationParams {
        epsilon_bulk: EPS_BULK,
        kappa_inv_bohr: kappa,
        tau_ha_bohr2: 2.0e-5,
        ..SolvationParams::vaspsol_water()
    }
}

fn grid_frac(global_idx: usize, dims: [usize; 3]) -> [f64; 3] {
    let iz = global_idx % dims[2];
    let iy = (global_idx / dims[2]) % dims[1];
    let ix = global_idx / (dims[1] * dims[2]);
    [
        ix as f64 / dims[0] as f64,
        iy as f64 / dims[1] as f64,
        iz as f64 / dims[2] as f64,
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn build_global_molecular_density(cell: &Cell, dims: [usize; 3], target_ne: f64) -> Vec<f64> {
    let n_grid = dims[0] * dims[1] * dims[2];
    let dv = cell.volume_bohr3() / n_grid as f64;
    let lobes = [
        ([0.43, 0.49, 0.50], 2.2, 1.15),
        ([0.57, 0.51, 0.50], 1.8, 1.15),
    ];
    let bg_floor = 1.0e-5;
    let mut raw = vec![0.0_f64; n_grid];

    for (idx, val) in raw.iter_mut().enumerate() {
        let u = grid_frac(idx, dims);
        let mut acc = 0.0;
        for &(center, weight, sigma) in &lobes {
            let inv_2s2 = 0.5 / (sigma * sigma);
            let du = [
                u[0] - center[0] - (u[0] - center[0]).round(),
                u[1] - center[1] - (u[1] - center[1]).round(),
                u[2] - center[2] - (u[2] - center[2]).round(),
            ];
            let dr = cell.fractional_to_cartesian_bohr(du);
            let r2 = dot3(dr, dr);
            acc += weight * (-r2 * inv_2s2).exp();
        }
        *val = acc;
    }

    let raw_int: f64 = raw.iter().sum::<f64>() * dv;
    let scale = (target_ne - bg_floor * cell.volume_bohr3()) / raw_int;
    for val in &mut raw {
        *val = bg_floor + *val * scale;
    }
    raw
}

fn perturb_density(base: &[f64], owned_indices: &[usize], dims: [usize; 3], amp: f64) -> Vec<f64> {
    base.iter()
        .zip(owned_indices)
        .map(|(&val, &idx)| {
            let [u, v, w] = grid_frac(idx, dims);
            let wave = 0.5 * (2.0 * PI * u).cos() - 0.3 * (2.0 * PI * v).sin()
                + 0.2 * (2.0 * PI * w).cos();
            val * (1.0 + amp * wave)
        })
        .collect()
}

fn median(mut values: Vec<f64>) -> f64 {
    values.sort_by(f64::total_cmp);
    let mid = values.len() / 2;
    if values.len().is_multiple_of(2) {
        0.5 * (values[mid - 1] + values[mid])
    } else {
        values[mid]
    }
}

/// Applies the exact discrete Fisicaro product-rule operator
/// `A_Fis phi = sqrt(eps) * (L_0 + q) (sqrt(eps) * phi)` matching `SolverKind::FisicaroPcg`.
fn apply_fisicaro_operator<K: Collective>(
    ws: &mut SolvationWorkspace,
    comm: &K,
    meta: &GridSpectralMeta,
    sqrt_eps: &[f64],
    q_fis: &[f64],
    phi_g: &[Complex64],
    out_g: &mut [Complex64],
) {
    let pp = sqrt_eps.len();
    let ll = phi_g.len();
    let mut phi_r = vec![0.0_f64; pp];
    let mut phi_tilde_r = vec![0.0_f64; pp];
    let mut phi_tilde_g = vec![Complex64::default(); ll];
    let mut l0_phi_tilde_g = vec![Complex64::default(); ll];
    let mut l0_phi_tilde_r = vec![0.0_f64; pp];

    ws.inverse_to_real(comm, phi_g, &mut phi_r).unwrap();
    for ((dst, &s), &p) in phi_tilde_r.iter_mut().zip(sqrt_eps).zip(&phi_r) {
        *dst = s * p;
    }
    ws.forward_real(comm, &phi_tilde_r, &mut phi_tilde_g)
        .unwrap();
    for ((dst, &g2), &pt) in l0_phi_tilde_g.iter_mut().zip(&meta.g2).zip(&phi_tilde_g) {
        *dst = g2 * pt;
    }
    ws.inverse_to_real(comm, &l0_phi_tilde_g, &mut l0_phi_tilde_r)
        .unwrap();
    for (((dst, &s), &l0), (&q, &pt)) in phi_r
        .iter_mut()
        .zip(sqrt_eps)
        .zip(&l0_phi_tilde_r)
        .zip(q_fis.iter().zip(&phi_tilde_r))
    {
        *dst = s * (l0 + q * pt);
    }
    ws.forward_real(comm, &phi_r, out_g).unwrap();
}

fn prepare_problem<K: Collective>(
    ws: &mut SolvationWorkspace,
    comm: &K,
    cell: &Cell,
    meta: &GridSpectralMeta,
    density_planes: &[f64],
    ions: &[IonicPseudoCharge],
    params: &SolvationParams,
) -> PreparedProblem {
    let pp = ws.plane_points();
    let ll = ws.local_positions().len();

    let mut s = vec![0.0; pp];
    let mut ds_dn = vec![0.0; pp];
    let mut eps = vec![0.0; pp];
    let mut deps_dn = vec![0.0; pp];
    let mut kappa2 = vec![0.0; pp];
    let mut dkappa2_dn = vec![0.0; pp];
    evaluate_cavity_fields(
        density_planes,
        params,
        &mut s,
        &mut ds_dn,
        &mut eps,
        &mut deps_dn,
        &mut kappa2,
        &mut dkappa2_dn,
    )
    .unwrap();

    let eps_min = eps.iter().copied().fold(f64::INFINITY, f64::min);
    let eps_max = eps.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    let eps_sum: f64 = eps.iter().copied().sum();
    let kap_sum: f64 = kappa2.iter().copied().sum();
    let mut sums = [0.0_f64; 2];
    comm.tree_sum_f64(&[eps_sum, kap_sum], &mut sums).unwrap();
    let has_kappa = sums[1] > 0.0;
    let eps_mean = sums[0] / meta.total_points as f64;
    let kappa2_bar = sums[1] / meta.total_points as f64;
    let eps_bar = (eps_min * eps_max).sqrt();

    let sqrt_eps: Vec<f64> = eps.iter().map(|&e| e.sqrt()).collect();
    let mut sqrt_eps_g = vec![Complex64::default(); ll];
    ws.forward_real(comm, &sqrt_eps, &mut sqrt_eps_g).unwrap();
    let mut neg_l0_s_g = vec![Complex64::default(); ll];
    for ((dst, &g2), &sg) in neg_l0_s_g.iter_mut().zip(&meta.g2).zip(&sqrt_eps_g) {
        *dst = -g2 * sg;
    }
    let mut neg_l0_s_r = vec![0.0_f64; pp];
    ws.inverse_to_real(comm, &neg_l0_s_g, &mut neg_l0_s_r)
        .unwrap();
    let q_fis: Vec<f64> = neg_l0_s_r
        .iter()
        .zip(&sqrt_eps)
        .zip(&kappa2)
        .map(|((&nl0, &se), &k2)| (nl0 + k2 / se) / se)
        .collect();

    let mut rho_ion_g = vec![Complex64::default(); ll];
    let _ =
        build_ionic_pseudo_charge_g(cell, ws.dims(), ws.local_positions(), ions, &mut rho_ion_g)
            .unwrap();
    let mut n_g = vec![Complex64::default(); ll];
    ws.forward_real(comm, density_planes, &mut n_g).unwrap();

    let rho_sol_g: Vec<Complex64> = rho_ion_g
        .iter()
        .zip(&n_g)
        .map(|(&ri, &ne)| ri - ne)
        .collect();
    let rhs_total_g: Vec<Complex64> = rho_sol_g.iter().map(|&r| FOUR_PI * r).collect();
    let phi_vac_g: Vec<Complex64> = meta
        .g2
        .iter()
        .zip(&rho_sol_g)
        .map(|(&g2, &rho)| {
            if g2 == 0.0 {
                Complex64::default()
            } else {
                (FOUR_PI / g2) * rho
            }
        })
        .collect();

    // Gate 2 discrete polarization RHS: b_pol_disc = 4*pi*rho_sol - A_disc * phi_vac
    let eps_minus_1: Vec<f64> = eps.iter().map(|&e| e - 1.0).collect();
    let mut b_pol_disc_g = vec![Complex64::default(); ll];
    ws.apply_operator(comm, &eps_minus_1, &kappa2, &phi_vac_g, &mut b_pol_disc_g)
        .unwrap();
    for v in &mut b_pol_disc_g {
        *v = -*v;
    }
    for &k in &meta.zero_g2_indices {
        b_pol_disc_g[k] += FOUR_PI * rho_sol_g[k];
    }

    // Consistent Fisicaro polarization RHS: b_pol_fis = 4*pi*rho_sol - A_Fis * phi_vac
    let mut a_fis_phi_vac_g = vec![Complex64::default(); ll];
    apply_fisicaro_operator(
        ws,
        comm,
        meta,
        &sqrt_eps,
        &q_fis,
        &phi_vac_g,
        &mut a_fis_phi_vac_g,
    );
    let b_pol_fis_g: Vec<Complex64> = rhs_total_g
        .iter()
        .zip(&a_fis_phi_vac_g)
        .map(|(&rhs, &afv)| rhs - afv)
        .collect();

    PreparedProblem {
        eps,
        kappa2,
        sqrt_eps,
        has_kappa,
        eps_bar,
        eps_mean,
        kappa2_bar,
        rho_sol_g,
        rhs_total_g,
        phi_vac_g,
        b_pol_disc_g,
        b_pol_fis_g,
    }
}

fn apply_operator_ffts(ws: &SolvationWorkspace, has_kappa: bool) -> usize {
    if ws.mirror_paired() {
        4
    } else if has_kappa {
        6
    } else {
        5
    }
}

/// Applies the inexact Fisicaro preconditioner `z = M_Fis^{-1}(r)` by scaling `r`
/// to unit L2 norm (avoiding `solve_pde`'s `max(1.0, ||b||_2)` floor), invoking
/// `ws.solve_pde(..., SolverKind::FisicaroPcg, max_iterations = inner_k)`, and
/// reconstructing `phi = phi_tilde / sqrt(eps)` when `NotConverged` stops early
/// at `inner_k` steps. Falls back to the diagonal G-space preconditioner on breakdown
/// or non-descent (`<r, z> <= 0`).
#[allow(clippy::too_many_arguments)]
fn apply_inexact_fisicaro_precond<K: Collective>(
    ws: &mut SolvationWorkspace,
    comm: &K,
    meta: &GridSpectralMeta,
    prob: &PreparedProblem,
    diag_precond: &[f64],
    r_g: &[Complex64],
    r_norm: f64,
    inner_k: usize,
    r_unit_scratch: &mut [Complex64],
    tmp_r_scratch: &mut [f64],
    z_out_g: &mut [Complex64],
) -> (usize, usize, bool) {
    if r_norm == 0.0 || !r_norm.is_finite() {
        z_out_g.fill(Complex64::default());
        return (0, 0, false);
    }

    let inv_norm = 1.0 / r_norm;
    for (dst, &r) in r_unit_scratch.iter_mut().zip(r_g) {
        *dst = r * inv_norm;
    }

    let inner_cfg = SolverConfig {
        kind: SolverKind::FisicaroPcg,
        tolerance: INNER_TOL,
        max_iterations: inner_k,
        diis_history: 8,
        warm_start: false,
    };

    let (iters, ffts, broke_down) = match ws.solve_pde(
        comm,
        &prob.eps,
        &prob.kappa2,
        r_unit_scratch,
        &inner_cfg,
        z_out_g,
    ) {
        Ok(rep) => (rep.iterations, rep.fft_count, false),
        Err(SolvationError::NotConverged { iterations, .. }) => {
            // `solve_fisicaro_pcg` leaves the `iterations`-step `phi_tilde_g` in `z_out_g`
            // before returning `NotConverged`. Reconstruct `phi = phi_tilde / sqrt(eps)`.
            ws.inverse_to_real(comm, z_out_g, tmp_r_scratch).unwrap();
            for (v, &se) in tmp_r_scratch.iter_mut().zip(&prob.sqrt_eps) {
                *v /= se;
            }
            ws.forward_real(comm, tmp_r_scratch, z_out_g).unwrap();
            meta.zero_nullspace(prob.has_kappa, z_out_g);
            // 3 setup + 2 initial residual + 2*iterations inner loop + 2 reconstruction.
            (iterations, 7 + 2 * iterations, false)
        }
        Err(SolvationError::Breakdown(_)) => (0, 5, true),
        Err(other) => panic!("unexpected inner Fisicaro preconditioner error: {other}"),
    };

    if !broke_down {
        for z in z_out_g.iter_mut() {
            *z *= r_norm;
        }
        meta.zero_nullspace(prob.has_kappa, z_out_g);
        let rz = meta.inner(comm, r_g, z_out_g);
        if rz.is_finite() && rz > 0.0 {
            return (iters, ffts, false);
        }
    }

    // Guaranteed-SPD diagonal fallback.
    for ((z, &pre), &r) in z_out_g.iter_mut().zip(diag_precond).zip(r_g) {
        *z = pre * r;
    }
    meta.zero_nullspace(prob.has_kappa, z_out_g);
    (iters, ffts, true)
}

/// Solves `A_disc * phi = b_pol_disc` via Flexible PCG (Polak-Ribiere) or
/// `A_disc`-orthogonal truncated Flexible CG / GCR(m) using inexact `FisicaroPcg`
/// preconditioning.
#[allow(clippy::too_many_arguments)]
fn solve_flexible_exact_disc<K: Collective>(
    ws: &mut SolvationWorkspace,
    comm: &K,
    meta: &GridSpectralMeta,
    prob: &PreparedProblem,
    update: FlexUpdate,
    inner_k: usize,
    tol: f64,
    warm: Option<&WarmState>,
    out_phi_g: &mut [Complex64],
) -> FlexSolveStats {
    let ll = prob.b_pol_disc_g.len();
    let pp = prob.eps.len();
    let op_ffts = apply_operator_ffts(ws, prob.has_kappa);

    let mut rhs_eff_g = prob.b_pol_disc_g.clone();
    meta.zero_nullspace(prob.has_kappa, &mut rhs_eff_g);
    let b_norm = meta.norm(comm, &rhs_eff_g);
    let b_scale = b_norm.max(1.0);
    let tol_abs = tol * b_scale;

    let diag_precond: Vec<f64> = meta
        .g2
        .iter()
        .map(|&g2| {
            if !prob.has_kappa && g2 == 0.0 {
                0.0
            } else {
                1.0 / (prob.eps_bar * g2 + prob.kappa2_bar)
            }
        })
        .collect();

    if let Some(prev) = warm {
        for (((phi, &prev_phi), &g2), (&rhs, &prev_rhs)) in out_phi_g
            .iter_mut()
            .zip(&prev.phi_g)
            .zip(&meta.g2)
            .zip(rhs_eff_g.iter().zip(&prev.rhs_g))
        {
            let bulk_pre = if !prob.has_kappa && g2 == 0.0 {
                0.0
            } else {
                1.0 / (prob.eps_mean * g2 + prob.kappa2_bar)
            };
            *phi = prev_phi + bulk_pre * (rhs - prev_rhs);
        }
    } else {
        for ((phi, &pre), &rhs) in out_phi_g.iter_mut().zip(&diag_precond).zip(&rhs_eff_g) {
            *phi = pre * rhs;
        }
    }
    meta.zero_nullspace(prob.has_kappa, out_phi_g);

    let mut w_g = vec![Complex64::default(); ll];
    ws.apply_operator(comm, &prob.eps, &prob.kappa2, out_phi_g, &mut w_g)
        .unwrap();
    meta.zero_nullspace(prob.has_kappa, &mut w_g);
    let mut outer_ffts = op_ffts;
    let mut inner_iters = 0usize;
    let mut inner_ffts = 0usize;
    let mut fallbacks = 0usize;

    let mut r_g: Vec<Complex64> = rhs_eff_g
        .iter()
        .zip(&w_g)
        .map(|(&rhs, &w)| rhs - w)
        .collect();
    let mut r_norm = meta.norm(comm, &r_g);
    if r_norm <= tol_abs {
        return FlexSolveStats {
            converged: true,
            outer_iters: 0,
            inner_iters,
            outer_ffts,
            inner_ffts,
            fallbacks,
            rel_residual: r_norm / b_scale,
        };
    }

    let mut r_unit_scratch = vec![Complex64::default(); ll];
    let mut tmp_r_scratch = vec![0.0_f64; pp];
    let mut z_g = vec![Complex64::default(); ll];
    let mut z_prev_g = vec![Complex64::default(); ll];
    let mut p_g = vec![Complex64::default(); ll];
    let mut rz_prev = 0.0_f64;

    // Circular buffer for exact A-orthogonal FCG(m) directions.
    let win_cap = match update {
        FlexUpdate::PolakRibiere => 0,
        FlexUpdate::FcgWindow(m) => m,
    };
    let mut hist_p: Vec<Vec<Complex64>> = Vec::with_capacity(win_cap);
    let mut hist_w: Vec<Vec<Complex64>> = Vec::with_capacity(win_cap);
    let mut hist_gamma: Vec<f64> = Vec::with_capacity(win_cap);
    let mut hist_cursor = 0usize;

    for outer in 1..=MAX_FLEX_OUTER_ITERS {
        let (it_k, fft_k, fb) = apply_inexact_fisicaro_precond(
            ws,
            comm,
            meta,
            prob,
            &diag_precond,
            &r_g,
            r_norm,
            inner_k,
            &mut r_unit_scratch,
            &mut tmp_r_scratch,
            &mut z_g,
        );
        inner_iters += it_k;
        inner_ffts += fft_k;
        if fb {
            fallbacks += 1;
        }

        match update {
            FlexUpdate::PolakRibiere => {
                let rz = meta.inner(comm, &r_g, &z_g);
                if outer == 1 || rz_prev <= 0.0 {
                    p_g.copy_from_slice(&z_g);
                } else {
                    // Polak-Ribiere: beta = max(0, <r_k, z_k - z_{k-1}> / <r_{k-1}, z_{k-1}>).
                    let rz_prev_dir = meta.inner(comm, &r_g, &z_prev_g);
                    let beta = ((rz - rz_prev_dir) / rz_prev).max(0.0);
                    for (p, &z) in p_g.iter_mut().zip(&z_g) {
                        *p = z + beta * *p;
                    }
                }
                z_prev_g.copy_from_slice(&z_g);
                rz_prev = rz;
            }
            FlexUpdate::FcgWindow(_) => {
                // Exact A-orthogonalization of z_k against the last <= m stored directions.
                p_g.copy_from_slice(&z_g);
                for ((pj_vec, wj_vec), &gamma_j) in hist_p.iter().zip(&hist_w).zip(&hist_gamma) {
                    let pw = meta.inner(comm, &p_g, wj_vec);
                    let coeff = pw / gamma_j;
                    for (p, &pj) in p_g.iter_mut().zip(pj_vec) {
                        *p -= coeff * pj;
                    }
                }
                // Guard against near-collinearity with history.
                let rp_check = meta.inner(comm, &r_g, &p_g);
                if !rp_check.is_finite() || rp_check <= 0.0 {
                    p_g.copy_from_slice(&z_g);
                    hist_p.clear();
                    hist_w.clear();
                    hist_gamma.clear();
                    hist_cursor = 0;
                }
            }
        }

        ws.apply_operator(comm, &prob.eps, &prob.kappa2, &p_g, &mut w_g)
            .unwrap();
        meta.zero_nullspace(prob.has_kappa, &mut w_g);
        outer_ffts += op_ffts;

        let pap = meta.inner(comm, &p_g, &w_g);
        let rp = meta.inner(comm, &r_g, &p_g);
        if !pap.is_finite() || pap <= 0.0 || !rp.is_finite() {
            return FlexSolveStats {
                converged: false,
                outer_iters: outer,
                inner_iters,
                outer_ffts,
                inner_ffts,
                fallbacks,
                rel_residual: r_norm / b_scale,
            };
        }

        let alpha = rp / pap;
        for ((phi, r), (&p, &w)) in out_phi_g.iter_mut().zip(&mut r_g).zip(p_g.iter().zip(&w_g)) {
            *phi += alpha * p;
            *r -= alpha * w;
        }
        r_norm = meta.norm(comm, &r_g);

        if r_norm <= tol_abs {
            meta.zero_nullspace(prob.has_kappa, out_phi_g);
            return FlexSolveStats {
                converged: true,
                outer_iters: outer,
                inner_iters,
                outer_ffts,
                inner_ffts,
                fallbacks,
                rel_residual: r_norm / b_scale,
            };
        }

        if let FlexUpdate::FcgWindow(m) = update
            && m > 0
        {
            if hist_p.len() < m {
                hist_p.push(p_g.clone());
                hist_w.push(w_g.clone());
                hist_gamma.push(pap);
            } else {
                hist_p[hist_cursor].copy_from_slice(&p_g);
                hist_w[hist_cursor].copy_from_slice(&w_g);
                hist_gamma[hist_cursor] = pap;
                hist_cursor = (hist_cursor + 1) % m;
            }
        }
    }

    meta.zero_nullspace(prob.has_kappa, out_phi_g);
    FlexSolveStats {
        converged: false,
        outer_iters: MAX_FLEX_OUTER_ITERS,
        inner_iters,
        outer_ffts,
        inner_ffts,
        fallbacks,
        rel_residual: r_norm / b_scale,
    }
}

#[allow(clippy::too_many_arguments)]
fn time_solve_pde<K: Collective>(
    ws: &mut SolvationWorkspace,
    comm: &K,
    eps: &[f64],
    kappa2: &[f64],
    rhs_g: &[Complex64],
    cfg: &SolverConfig,
    repeats: usize,
    out_phi_g: &mut [Complex64],
) -> (Result<SolveReport, SolvationError>, f64) {
    let mut times = Vec::with_capacity(repeats);
    let mut last_res = None;
    for _ in 0..repeats {
        ws.reset_warm_start();
        let t0 = Instant::now();
        let res = ws.solve_pde(comm, eps, kappa2, rhs_g, cfg, out_phi_g);
        times.push(t0.elapsed().as_secs_f64() * 1.0e3);
        last_res = Some(res);
    }
    (last_res.unwrap(), median(times))
}

fn run_experiment_c1(grids: &[usize], repeats: usize) {
    println!(
        "## Experiment C1 — Consistent Fisicaro Polarization RHS (`tol = 1e-10`, cold, median of {repeats} run(s))"
    );
    println!();
    println!(
        "Comparison on sharp VASPsol-erfc water (`L = {BOX_LENGTH_BOHR:.0} Bohr`, `eps_b = {EPS_BULK:.1}`, `tau = 2e-5`):"
    );
    println!("- `DirectPcg`: exact discrete SPD operator (`A_disc phi_pol = b_pol_disc`)");
    println!("- `Fisicaro (Gate 2 b_pol_disc)`: Gate 2 RHS mixing `A_Fis` LHS with `A_disc` RHS");
    println!(
        "- `Fisicaro (Total 4pi*rho_sol)`: consistent total-potential solve `A_Fis phi_solv = 4*pi*rho_sol`, then `phi_pol = phi_solv - phi_vac` (saves 4 `b_pol` setup FFTs)"
    );
    println!(
        "- `Fisicaro (Consistent b_pol_Fis)`: consistent polarization solve `A_Fis phi_pol = 4*pi*rho_sol - A_Fis phi_vac`"
    );
    println!();
    println!(
        "| Grid (`dx` Bohr) | `kappa` | Arm | Status | Iters | PDE FFTs | Wall (ms) | Speedup | `E_es` (Ha) | `|E - E_direct|` (Ha) | Rel Res |"
    );
    println!("| --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |");

    let cell = cubic_cell(BOX_LENGTH_BOHR);
    let ions = water_ions();

    for &n_axis in grids {
        let dims = [n_axis, n_axis, n_axis];
        let dx = BOX_LENGTH_BOHR / n_axis as f64;
        let global_n0 = build_global_molecular_density(&cell, dims, 4.0);

        for kappa in [0.0_f64, 0.15_f64] {
            let params = vaspsol_params(kappa);
            dft_mpi::simulate(1, |comm| {
                let mut ws = SolvationWorkspace::new(&cell, dims, comm).unwrap();
                let owned = ws.owned_real_grid_indices();
                let n0: Vec<f64> = owned.iter().map(|&idx| global_n0[idx]).collect();
                let meta = GridSpectralMeta::new(&cell, dims, ws.local_positions());
                let prob = prepare_problem(&mut ws, comm, &cell, &meta, &n0, &ions, &params);
                let ll = ws.local_positions().len();

                let cfg_direct = SolverConfig {
                    kind: SolverKind::DirectPcg,
                    tolerance: COLD_TOL,
                    max_iterations: MAX_DIRECT_ITERS,
                    diis_history: 8,
                    warm_start: false,
                };
                let cfg_fis = SolverConfig {
                    kind: SolverKind::FisicaroPcg,
                    tolerance: COLD_TOL,
                    max_iterations: MAX_DIRECT_ITERS,
                    diis_history: 8,
                    warm_start: false,
                };

                // Sanity cross-check that our `DirectPcg` `E_es` matches `ws.evaluate` identically.
                ws.reset_warm_start();
                let eval_direct = ws.evaluate(comm, &n0, &ions, &params, &cfg_direct).unwrap();

                // (a) DirectPcg
                let mut phi_direct_g = vec![Complex64::default(); ll];
                let (direct_res, direct_ms) = time_solve_pde(
                    &mut ws,
                    comm,
                    &prob.eps,
                    &prob.kappa2,
                    &prob.b_pol_disc_g,
                    &cfg_direct,
                    repeats,
                    &mut phi_direct_g,
                );
                let direct_rep = direct_res.unwrap();
                let e_direct = meta.electrostatic_energy_ha(comm, &prob.rho_sol_g, &phi_direct_g);
                assert!(
                    (e_direct - eval_direct.electrostatic_ha).abs() < 1.0e-13,
                    "DirectPcg E_es mismatch vs ws.evaluate: {e_direct} vs {}",
                    eval_direct.electrostatic_ha
                );
                println!(
                    "| `{n_axis}^3` (`{dx:.3}`) | {kappa:.2} | `DirectPcg` | OK | {} | {} | {direct_ms:.2} | 1.00x | {e_direct:.10} | `0` | `{:.2e}` |",
                    direct_rep.iterations, direct_rep.fft_count, direct_rep.relative_residual
                );

                // (b) Gate 2 FisicaroPcg (b_pol_disc)
                let mut phi_g = vec![Complex64::default(); ll];
                let (res_b, ms_b) = time_solve_pde(
                    &mut ws,
                    comm,
                    &prob.eps,
                    &prob.kappa2,
                    &prob.b_pol_disc_g,
                    &cfg_fis,
                    repeats,
                    &mut phi_g,
                );
                print_c1_row(
                    n_axis,
                    dx,
                    kappa,
                    "`Fisicaro (Gate 2 b_pol_disc)`",
                    &res_b,
                    ms_b,
                    direct_ms,
                    meta.electrostatic_energy_ha(comm, &prob.rho_sol_g, &phi_g),
                    e_direct,
                );

                // (c) Consistent-RHS FisicaroPcg (Total 4pi*rho_sol -> phi_solv - phi_vac)
                let mut phi_solv_g = vec![Complex64::default(); ll];
                let (res_c, ms_c) = time_solve_pde(
                    &mut ws,
                    comm,
                    &prob.eps,
                    &prob.kappa2,
                    &prob.rhs_total_g,
                    &cfg_fis,
                    repeats,
                    &mut phi_solv_g,
                );
                for ((dst, &ps), &pv) in phi_g.iter_mut().zip(&phi_solv_g).zip(&prob.phi_vac_g) {
                    *dst = ps - pv;
                }
                print_c1_row(
                    n_axis,
                    dx,
                    kappa,
                    "`Fisicaro (Total 4pi*rho_sol)`",
                    &res_c,
                    ms_c,
                    direct_ms,
                    meta.electrostatic_energy_ha(comm, &prob.rho_sol_g, &phi_g),
                    e_direct,
                );

                // (d) Consistent-RHS FisicaroPcg (Explicit b_pol_Fis = 4pi*rho_sol - A_Fis phi_vac)
                let (res_d, ms_d) = time_solve_pde(
                    &mut ws,
                    comm,
                    &prob.eps,
                    &prob.kappa2,
                    &prob.b_pol_fis_g,
                    &cfg_fis,
                    repeats,
                    &mut phi_g,
                );
                print_c1_row(
                    n_axis,
                    dx,
                    kappa,
                    "`Fisicaro (Consistent b_pol_Fis)`",
                    &res_d,
                    ms_d,
                    direct_ms,
                    meta.electrostatic_energy_ha(comm, &prob.rho_sol_g, &phi_g),
                    e_direct,
                );
            });
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn print_c1_row(
    n_axis: usize,
    dx: f64,
    kappa: f64,
    arm: &str,
    res: &Result<SolveReport, SolvationError>,
    ms: f64,
    direct_ms: f64,
    e_es: f64,
    e_direct: f64,
) {
    let speedup = direct_ms / ms.max(1.0e-9);
    match res {
        Ok(rep) => println!(
            "| `{n_axis}^3` (`{dx:.3}`) | {kappa:.2} | {arm} | OK | {} | {} | {ms:.2} | {speedup:.2}x | {e_es:.10} | `{:.2e}` | `{:.2e}` |",
            rep.iterations,
            rep.fft_count,
            (e_es - e_direct).abs(),
            rep.relative_residual
        ),
        Err(SolvationError::NotConverged {
            iterations,
            residual_norm,
        }) => println!(
            "| `{n_axis}^3` (`{dx:.3}`) | {kappa:.2} | {arm} | NOT_CONV | {iterations} | - | {ms:.2} | {speedup:.2}x | - | - | `{residual_norm:.2e}` |"
        ),
        Err(SolvationError::Breakdown(reason)) => println!(
            "| `{n_axis}^3` (`{dx:.3}`) | {kappa:.2} | {arm} | BREAKDOWN({reason}) | - | - | {ms:.2} | - | - | - | - |"
        ),
        Err(err) => println!(
            "| `{n_axis}^3` (`{dx:.3}`) | {kappa:.2} | {arm} | ERR({err}) | - | - | {ms:.2} | - | - | - | - |"
        ),
    }
}

fn run_experiment_c2(grids: &[usize], repeats: usize) {
    println!();
    println!(
        "## Experiment C2 — Inexact Fisicaro Preconditioner inside Flexible Krylov Targeting Exact `A_disc` (median of {repeats} run(s))"
    );
    println!();
    println!(
        "Outer loop evaluates exact residual `r_k = b_pol_disc - A_disc phi_k` via `ws.apply_operator`; inner preconditioner calls `ws.solve_pde(..., FisicaroPcg, max_iterations = k, tol = {INNER_TOL:.0e})` on unit-scaled residual with diagonal fallback. Energy difference `|E - E_direct(1e-10)|` is measured against exact `1e-10`-converged `DirectPcg` on the same density."
    );
    println!();
    println!(
        "| Grid (`dx` Bohr) | `kappa` | Mode | Solver | Status | Outer | Inner | Total FFTs | Wall (ms) | vs Direct | `|E - E_direct|` (Ha) | Rel Res | FB |"
    );
    println!(
        "| --- | ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |"
    );

    let cell = cubic_cell(BOX_LENGTH_BOHR);
    let ions = water_ions();
    let flex_arms = [
        (FlexUpdate::PolakRibiere, 4usize),
        (FlexUpdate::PolakRibiere, 8usize),
        (FlexUpdate::FcgWindow(FCG_WINDOW), 4usize),
        (FlexUpdate::FcgWindow(FCG_WINDOW), 8usize),
    ];

    for &n_axis in grids {
        let dims = [n_axis, n_axis, n_axis];
        let dx = BOX_LENGTH_BOHR / n_axis as f64;
        let global_n0 = build_global_molecular_density(&cell, dims, 4.0);

        for kappa in [0.0_f64, 0.15_f64] {
            let params = vaspsol_params(kappa);
            dft_mpi::simulate(1, |comm| {
                let mut ws = SolvationWorkspace::new(&cell, dims, comm).unwrap();
                let owned = ws.owned_real_grid_indices();
                let n0: Vec<f64> = owned.iter().map(|&idx| global_n0[idx]).collect();
                let n1 = perturb_density(&n0, &owned, dims, 1.0e-3);
                let meta = GridSpectralMeta::new(&cell, dims, ws.local_positions());
                let prob0 = prepare_problem(&mut ws, comm, &cell, &meta, &n0, &ions, &params);
                let prob1 = prepare_problem(&mut ws, comm, &cell, &meta, &n1, &ions, &params);
                let ll = ws.local_positions().len();

                let cfg_direct_cold = SolverConfig {
                    kind: SolverKind::DirectPcg,
                    tolerance: COLD_TOL,
                    max_iterations: MAX_DIRECT_ITERS,
                    diis_history: 8,
                    warm_start: false,
                };
                let cfg_direct_warm_seed = SolverConfig {
                    kind: SolverKind::DirectPcg,
                    tolerance: COLD_TOL,
                    max_iterations: MAX_DIRECT_ITERS,
                    diis_history: 8,
                    warm_start: true,
                };
                let cfg_direct_warm_run = SolverConfig {
                    tolerance: WARM_TOL,
                    ..cfg_direct_warm_seed
                };

                // Exact reference energies on n0 and n1 at tol = 1e-10.
                let mut phi0_direct_g = vec![Complex64::default(); ll];
                let (rep0_direct, ms0_direct) = time_solve_pde(
                    &mut ws,
                    comm,
                    &prob0.eps,
                    &prob0.kappa2,
                    &prob0.b_pol_disc_g,
                    &cfg_direct_cold,
                    repeats,
                    &mut phi0_direct_g,
                );
                let rep0_direct = rep0_direct.unwrap();
                let e0_direct =
                    meta.electrostatic_energy_ha(comm, &prob0.rho_sol_g, &phi0_direct_g);

                let mut phi1_ref_g = vec![Complex64::default(); ll];
                ws.reset_warm_start();
                ws.solve_pde(
                    comm,
                    &prob1.eps,
                    &prob1.kappa2,
                    &prob1.b_pol_disc_g,
                    &cfg_direct_cold,
                    &mut phi1_ref_g,
                )
                .unwrap();
                let e1_direct_ref =
                    meta.electrostatic_energy_ha(comm, &prob1.rho_sol_g, &phi1_ref_g);

                // Warm DirectPcg baseline on n1 @ 1e-5.
                let mut phi1_direct_warm_g = vec![Complex64::default(); ll];
                let mut direct_warm_times = Vec::with_capacity(repeats);
                let mut rep1_direct_warm = None;
                for _ in 0..repeats {
                    ws.reset_warm_start();
                    let mut seed_phi = vec![Complex64::default(); ll];
                    ws.solve_pde(
                        comm,
                        &prob0.eps,
                        &prob0.kappa2,
                        &prob0.b_pol_disc_g,
                        &cfg_direct_warm_seed,
                        &mut seed_phi,
                    )
                    .unwrap();
                    let tw = Instant::now();
                    let rep = ws
                        .solve_pde(
                            comm,
                            &prob1.eps,
                            &prob1.kappa2,
                            &prob1.b_pol_disc_g,
                            &cfg_direct_warm_run,
                            &mut phi1_direct_warm_g,
                        )
                        .unwrap();
                    direct_warm_times.push(tw.elapsed().as_secs_f64() * 1.0e3);
                    rep1_direct_warm = Some(rep);
                }
                let rep1_direct_warm = rep1_direct_warm.unwrap();
                let ms1_direct_warm = median(direct_warm_times);
                let e1_direct_warm =
                    meta.electrostatic_energy_ha(comm, &prob1.rho_sol_g, &phi1_direct_warm_g);

                let mut rhs0_eff_g = prob0.b_pol_disc_g.clone();
                meta.zero_nullspace(prob0.has_kappa, &mut rhs0_eff_g);
                let warm_state = WarmState {
                    phi_g: phi0_direct_g.clone(),
                    rhs_g: rhs0_eff_g,
                };

                // --- Cold @ 1e-10 ---
                println!(
                    "| `{n_axis}^3` (`{dx:.3}`) | {kappa:.2} | `cold@1e-10` | `DirectPcg` | OK | {} | 0 | {} | {ms0_direct:.2} | 1.00x | `0` | `{:.2e}` | 0 |",
                    rep0_direct.iterations, rep0_direct.fft_count, rep0_direct.relative_residual
                );
                for &(update, inner_k) in &flex_arms {
                    let mut phi_out = vec![Complex64::default(); ll];
                    let mut times = Vec::with_capacity(repeats);
                    let mut last_stats = None;
                    for _ in 0..repeats {
                        let t0 = Instant::now();
                        let stats = solve_flexible_exact_disc(
                            &mut ws,
                            comm,
                            &meta,
                            &prob0,
                            update,
                            inner_k,
                            COLD_TOL,
                            None,
                            &mut phi_out,
                        );
                        times.push(t0.elapsed().as_secs_f64() * 1.0e3);
                        last_stats = Some(stats);
                    }
                    let stats = last_stats.unwrap();
                    let ms = median(times);
                    let e_es = meta.electrostatic_energy_ha(comm, &prob0.rho_sol_g, &phi_out);
                    print_c2_row(
                        n_axis,
                        dx,
                        kappa,
                        "`cold@1e-10`",
                        &format!("`{} (Fis k={inner_k})`", update.label()),
                        &stats,
                        ms,
                        ms0_direct,
                        (e_es - e0_direct).abs(),
                    );
                }

                // --- Warm(1e-3) @ 1e-5 ---
                println!(
                    "| `{n_axis}^3` (`{dx:.3}`) | {kappa:.2} | `warm(1e-3)@1e-5` | `DirectPcg` | OK | {} | 0 | {} | {ms1_direct_warm:.2} | 1.00x | `{:.2e}` | `{:.2e}` | 0 |",
                    rep1_direct_warm.iterations,
                    rep1_direct_warm.fft_count,
                    (e1_direct_warm - e1_direct_ref).abs(),
                    rep1_direct_warm.relative_residual
                );
                for &(update, inner_k) in &flex_arms {
                    let mut phi_out = vec![Complex64::default(); ll];
                    let mut times = Vec::with_capacity(repeats);
                    let mut last_stats = None;
                    for _ in 0..repeats {
                        let t0 = Instant::now();
                        let stats = solve_flexible_exact_disc(
                            &mut ws,
                            comm,
                            &meta,
                            &prob1,
                            update,
                            inner_k,
                            WARM_TOL,
                            Some(&warm_state),
                            &mut phi_out,
                        );
                        times.push(t0.elapsed().as_secs_f64() * 1.0e3);
                        last_stats = Some(stats);
                    }
                    let stats = last_stats.unwrap();
                    let ms = median(times);
                    let e_es = meta.electrostatic_energy_ha(comm, &prob1.rho_sol_g, &phi_out);
                    print_c2_row(
                        n_axis,
                        dx,
                        kappa,
                        "`warm(1e-3)@1e-5`",
                        &format!("`{} (Fis k={inner_k})`", update.label()),
                        &stats,
                        ms,
                        ms1_direct_warm,
                        (e_es - e1_direct_ref).abs(),
                    );
                }
            });
        }
    }
}

#[allow(clippy::too_many_arguments)]
fn print_c2_row(
    n_axis: usize,
    dx: f64,
    kappa: f64,
    mode: &str,
    solver: &str,
    stats: &FlexSolveStats,
    ms: f64,
    direct_ms: f64,
    de_ha: f64,
) {
    let status = if stats.converged { "OK" } else { "NOT_CONV" };
    let speedup = direct_ms / ms.max(1.0e-9);
    println!(
        "| `{n_axis}^3` (`{dx:.3}`) | {kappa:.2} | {mode} | {solver} | {status} | {} | {} | {} | {ms:.2} | {speedup:.2}x | `{de_ha:.2e}` | `{:.2e}` | {} |",
        stats.outer_iters,
        stats.inner_iters,
        stats.total_ffts(),
        stats.rel_residual,
        stats.fallbacks
    );
}

fn parse_args() -> (Vec<usize>, usize) {
    let args: Vec<String> = std::env::args().collect();
    let mut grids = vec![32usize, 48, 64];
    let mut repeats = 3usize;
    let mut i = 1usize;
    while i < args.len() {
        match args[i].as_str() {
            "--96" => {
                if !grids.contains(&96) {
                    grids.push(96);
                }
            }
            "--repeats" => {
                i += 1;
                if let Some(val) = args.get(i) {
                    repeats = val
                        .parse::<usize>()
                        .expect("--repeats must be a positive integer")
                        .max(1);
                }
            }
            "--grids" => {
                i += 1;
                if let Some(val) = args.get(i) {
                    grids = val
                        .split(',')
                        .map(|s| {
                            s.trim()
                                .parse::<usize>()
                                .expect("--grids expects comma-separated integers")
                        })
                        .collect();
                }
            }
            other => {
                eprintln!("warning: ignoring unrecognized argument `{other}`");
            }
        }
        i += 1;
    }
    (grids, repeats)
}

fn main() {
    let (grids, repeats) = parse_args();
    println!("# Track B Gate 4 — Suite C: Implicit Solvation Solver Acceleration R&D");
    println!();
    run_experiment_c1(&grids, repeats);
    run_experiment_c2(&grids, repeats);
}
