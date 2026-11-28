//! Independent L1/L2 verification suite for Track B Gate 2 implicit solvation
//! (`dft_physics::solvation`).
//!
//! Written strictly against the mathematical specification in
//! `docs/research/extensible-dft/implicit-solvation-theory.md` (Sections 4 & 9)
//! and the public API contract in `crates/dft-physics/src/solvation/mod.rs`.

use std::f64::consts::PI;

use dft_geometry::Cell;
use dft_mpi::Collective;
use dft_physics::solvation::{
    CavityModel, IonicPseudoCharge, SolvationError, SolvationEvaluation, SolvationParams,
    SolvationWorkspace, SolverConfig, SolverKind,
};
use num_complex::Complex64;

/// Simple-cubic Wigner-Madelung constant `alpha_M`.
const MADELUNG_SC: f64 = 2.837_297_479;

struct MirrorOverrideGuard;

impl MirrorOverrideGuard {
    fn set(value: Option<bool>) -> Self {
        dft_solvers::distributed_fft::set_mirror_sticks_override(value);
        Self
    }
}

impl Drop for MirrorOverrideGuard {
    fn drop(&mut self) {
        dft_solvers::distributed_fft::set_mirror_sticks_override(None);
    }
}

fn skewed_cell_default() -> Cell {
    Cell::from_user_rows_bohr([[10.4, 0.8, -0.3], [0.5, 11.2, 0.7], [-0.4, 0.6, 10.8]]).unwrap()
}

fn cubic_cell(length_bohr: f64) -> Cell {
    Cell::from_user_rows_bohr([
        [length_bohr, 0.0, 0.0],
        [0.0, length_bohr, 0.0],
        [0.0, 0.0, length_bohr],
    ])
    .unwrap()
}

fn grid_ijk(global_idx: usize, dims: [usize; 3]) -> [usize; 3] {
    let iz = global_idx % dims[2];
    let iy = (global_idx / dims[2]) % dims[1];
    let ix = global_idx / (dims[1] * dims[2]);
    [ix, iy, iz]
}

fn grid_frac(global_idx: usize, dims: [usize; 3]) -> [f64; 3] {
    let [ix, iy, iz] = grid_ijk(global_idx, dims);
    [
        ix as f64 / dims[0] as f64,
        iy as f64 / dims[1] as f64,
        iz as f64 / dims[2] as f64,
    ]
}

fn dot3(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

/// Evaluates a C-infinity periodic sum of isotropic Gaussians centered at
/// `(frac_center, weight, sigma_bohr)` over `[-1, 0, 1]^3` periodic images,
/// then rescales so the exact discrete integral equals `target_ne` (after
/// adding `bg_floor`).
fn build_global_gaussian_density(
    cell: &Cell,
    dims: [usize; 3],
    lobes: &[([f64; 3], f64, f64)],
    target_ne: f64,
    bg_floor: f64,
) -> Vec<f64> {
    let n_grid = dims[0] * dims[1] * dims[2];
    let dv = cell.volume_bohr3() / n_grid as f64;
    let mut raw = vec![0.0_f64; n_grid];

    for (idx, val) in raw.iter_mut().enumerate() {
        let u = grid_frac(idx, dims);
        let mut acc = 0.0;
        for &(center, weight, sigma) in lobes {
            let inv_2s2 = 0.5 / (sigma * sigma);
            let du0 = [
                u[0] - center[0] - (u[0] - center[0]).round(),
                u[1] - center[1] - (u[1] - center[1]).round(),
                u[2] - center[2] - (u[2] - center[2]).round(),
            ];
            for sx in -1..=1 {
                for sy in -1..=1 {
                    for sz in -1..=1 {
                        let du = [du0[0] + sx as f64, du0[1] + sy as f64, du0[2] + sz as f64];
                        let dr = cell.fractional_to_cartesian_bohr(du);
                        let r2 = dot3(dr, dr);
                        acc += weight * (-r2 * inv_2s2).exp();
                    }
                }
            }
        }
        *val = acc;
    }

    let raw_integral: f64 = raw.iter().sum::<f64>() * dv;
    let bg_integral = bg_floor * cell.volume_bohr3();
    assert!(
        target_ne > bg_integral && raw_integral > 0.0,
        "invalid density normalization inputs"
    );
    let scale = (target_ne - bg_integral) / raw_integral;
    for val in &mut raw {
        *val = bg_floor + *val * scale;
    }
    raw
}

fn slice_owned(global: &[f64], owned_indices: &[usize]) -> Vec<f64> {
    owned_indices.iter().map(|&idx| global[idx]).collect()
}

fn max_abs(slice: &[f64]) -> f64 {
    slice.iter().map(|v| v.abs()).fold(0.0_f64, f64::max)
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

fn max_force_abs(forces: &[[f64; 3]]) -> f64 {
    forces
        .iter()
        .flat_map(|f| f.iter())
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max)
}

fn max_force_diff(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .flat_map(|(fa, fb)| fa.iter().zip(fb.iter()))
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

/// Evaluates the smooth 3D periodic MMS fields `(eps, kappa2, phi_exact, rhs_analytic)`
/// on the owned real-space grid of `cell` and `dims`.
///
/// Using phase angles `theta_k = 2 * pi * u_k = b_k . r`, the exact calculus
/// differential operator is:
/// `-div(eps grad phi) + kappa2 * phi = -sum_{i,j} (b_i . b_j) * d/dtheta_i (eps * dphi/dtheta_j) + kappa2 * phi`.
fn build_mms_fields(
    cell: &Cell,
    dims: [usize; 3],
    owned_indices: &[usize],
    eps_bulk: f64,
    with_screening: bool,
) -> (Vec<f64>, Vec<f64>, Vec<f64>, Vec<f64>) {
    let b_cols = cell.reciprocal_columns_inverse_bohr();
    let b_vecs: [[f64; 3]; 3] = std::array::from_fn(|j| [b_cols[0][j], b_cols[1][j], b_cols[2][j]]);
    let g_metric: [[f64; 3]; 3] =
        std::array::from_fn(|i| std::array::from_fn(|j| dot3(b_vecs[i], b_vecs[j])));

    let n = owned_indices.len();
    let mut eps = Vec::with_capacity(n);
    let mut kappa2 = Vec::with_capacity(n);
    let mut phi = Vec::with_capacity(n);
    let mut rhs = Vec::with_capacity(n);

    let dc = if with_screening { 0.35 } else { 0.0 };
    let eps_span = eps_bulk - 1.0;

    for &idx in owned_indices {
        let [u, v, w] = grid_frac(idx, dims);
        let t0 = 2.0 * PI * u;
        let t1 = 2.0 * PI * v;
        let t2 = 2.0 * PI * w;

        // Smooth periodic permittivity spanning [1.0, eps_bulk].
        let f = 0.5 + 0.25 * t0.cos() + 0.15 * t1.sin() + 0.10 * t2.cos();
        let e = 1.0 + eps_span * f;
        let de = [
            eps_span * (-0.25 * t0.sin()),
            eps_span * (0.15 * t1.cos()),
            eps_span * (-0.10 * t2.sin()),
        ];

        let k2 = if with_screening {
            0.25 * e / eps_bulk
        } else {
            0.0
        };

        // Low-frequency periodic potential (harmonics |m_k| <= 1, so eps * grad phi
        // has |m_k| <= 2 < min(dims)/2, strictly alias-free on all tested grids).
        let p =
            dc + 0.7 * t0.cos() - 0.5 * t1.sin() + 0.4 * (t0 + t2).cos() + 0.3 * (t1 - t2).sin();

        let dp = [
            -0.7 * t0.sin() - 0.4 * (t0 + t2).sin(),
            -0.5 * t1.cos() + 0.3 * (t1 - t2).cos(),
            -0.4 * (t0 + t2).sin() - 0.3 * (t1 - t2).cos(),
        ];

        let d2p = [
            [
                -0.7 * t0.cos() - 0.4 * (t0 + t2).cos(),
                0.0,
                -0.4 * (t0 + t2).cos(),
            ],
            [
                0.0,
                0.5 * t1.sin() - 0.3 * (t1 - t2).sin(),
                0.3 * (t1 - t2).sin(),
            ],
            [
                -0.4 * (t0 + t2).cos(),
                0.3 * (t1 - t2).sin(),
                -0.4 * (t0 + t2).cos() - 0.3 * (t1 - t2).sin(),
            ],
        ];

        let mut div_eps_grad_p = 0.0;
        for ((g_row, &de_i), d2p_row) in g_metric.iter().zip(&de).zip(&d2p) {
            for ((&g_ij, &dp_j), &d2p_ij) in g_row.iter().zip(&dp).zip(d2p_row) {
                div_eps_grad_p += g_ij * (de_i * dp_j + e * d2p_ij);
            }
        }
        let r = -div_eps_grad_p + k2 * p;

        eps.push(e);
        kappa2.push(k2);
        phi.push(p);
        rhs.push(r);
    }

    (eps, kappa2, phi, rhs)
}

/// 1. Exact vacuum limit on an arbitrary skewed triclinic cell (L1).
#[test]
fn vacuum_limit_is_identically_zero_on_skewed_cell() {
    let cell = skewed_cell_default();
    let dims = [16, 18, 20];
    let ions = [
        IonicPseudoCharge {
            charge: 2.0,
            fractional_position: [0.37, 0.44, 0.51],
            width_bohr: 0.60,
        },
        IonicPseudoCharge {
            charge: 3.0,
            fractional_position: [0.61, 0.57, 0.43],
            width_bohr: 0.65,
        },
    ];
    let q_ion: f64 = ions.iter().map(|i| i.charge).sum();
    let global_n = build_global_gaussian_density(
        &cell,
        dims,
        &[
            ([0.37, 0.44, 0.51], 2.0, 1.15),
            ([0.61, 0.57, 0.43], 3.0, 1.25),
        ],
        q_ion,
        1.0e-5,
    );

    let params = SolvationParams::vacuum();
    let config = SolverConfig::default();

    let results = dft_mpi::simulate(1, |comm| {
        let mut ws = SolvationWorkspace::new(&cell, dims, comm).unwrap();
        let local_n = slice_owned(&global_n, &ws.owned_real_grid_indices());
        ws.evaluate(comm, &local_n, &ions, &params, &config)
            .unwrap()
    });
    let eval = &results[0];

    assert!(
        eval.electrostatic_ha.abs() < 1.0e-14,
        "electrostatic_ha = {:.3e}",
        eval.electrostatic_ha
    );
    assert!(
        eval.cavitation_ha.abs() < 1.0e-14,
        "cavitation_ha = {:.3e}",
        eval.cavitation_ha
    );
    assert!(
        eval.background_correction_ha.abs() < 1.0e-14,
        "background_correction_ha = {:.3e}",
        eval.background_correction_ha
    );
    assert!(
        eval.total_ha.abs() < 1.0e-14,
        "total_ha = {:.3e}",
        eval.total_ha
    );
    assert!(
        max_abs(&eval.phi_pol_planes) < 1.0e-14,
        "max |phi_pol| = {:.3e}",
        max_abs(&eval.phi_pol_planes)
    );
    assert!(
        max_abs(&eval.v_total_planes) < 1.0e-14,
        "max |v_total| = {:.3e}",
        max_abs(&eval.v_total_planes)
    );
    assert!(
        max_force_abs(&eval.ionic_forces_ha_bohr) < 1.0e-14,
        "max |ionic_forces| = {:.3e}",
        max_force_abs(&eval.ionic_forces_ha_bohr)
    );
}

/// Test 2: Method of Manufactured Solutions on cubic and skewed cells for `kappa = 0`
/// and `kappa > 0`, verifying both `apply_operator` against analytic calculus
/// and all three solver arms (`DirectPcg`, `FisicaroPcg`, `PolarizationDiis`) (L1).
#[test]
fn mms_variable_dielectric_and_lpb_cubic_and_skewed_all_solvers() {
    let cases = [
        ("cubic", cubic_cell(10.0), [16, 16, 16]),
        ("skewed", skewed_cell_default(), [16, 18, 14]),
    ];

    for (label, cell, dims) in cases {
        for with_screening in [false, true] {
            dft_mpi::simulate(1, |comm| {
                let mut ws = SolvationWorkspace::new(&cell, dims, comm).unwrap();
                let owned = ws.owned_real_grid_indices();
                let (eps, kappa2, phi_exact, rhs_analytic) =
                    build_mms_fields(&cell, dims, &owned, 78.4, with_screening);

                let n_g = ws.local_positions().len();
                let n_r = ws.plane_points();
                let mut phi_exact_g = vec![Complex64::new(0.0, 0.0); n_g];
                let mut rhs_discrete_g = vec![Complex64::new(0.0, 0.0); n_g];
                let mut rhs_discrete_r = vec![0.0_f64; n_r];

                ws.forward_real(comm, &phi_exact, &mut phi_exact_g).unwrap();
                ws.apply_operator(comm, &eps, &kappa2, &phi_exact_g, &mut rhs_discrete_g)
                    .unwrap();
                ws.inverse_to_real(comm, &rhs_discrete_g, &mut rhs_discrete_r)
                    .unwrap();

                let op_err = max_abs_diff(&rhs_discrete_r, &rhs_analytic);
                let rhs_scale = max_abs(&rhs_analytic).max(1.0);
                assert!(
                    op_err / rhs_scale < 1.0e-11,
                    "{label} (screen={with_screening}): apply_operator vs analytic rel_err = {:.3e} (abs {:.3e})",
                    op_err / rhs_scale,
                    op_err
                );

                let cfg = SolverConfig {
                    kind: SolverKind::DirectPcg,
                    tolerance: 1.0e-14,
                    max_iterations: 400,
                    diis_history: 8,
                    warm_start: false,
                };
                let mut phi_num_g = vec![Complex64::new(0.0, 0.0); n_g];
                let report = ws
                    .solve_pde(comm, &eps, &kappa2, &rhs_discrete_g, &cfg, &mut phi_num_g)
                    .unwrap();
                assert!(
                    report.converged,
                    "{label} DirectPcg failed to converge: {report:?}"
                );

                let mut phi_num_r = vec![0.0_f64; n_r];
                ws.inverse_to_real(comm, &phi_num_g, &mut phi_num_r)
                    .unwrap();
                let rel_linf = max_abs_diff(&phi_num_r, &phi_exact) / max_abs(&phi_exact);
                println!(
                    "[MMS {label} eps_b=78.4 screen={with_screening}] DirectPcg: iters={} ffts={} op_rel_err={:.3e} rel_linf={rel_linf:.3e}",
                    report.iterations,
                    report.fft_count,
                    op_err / rhs_scale
                );
                assert!(
                    rel_linf < 1.0e-12,
                    "{label} (screen={with_screening}) DirectPcg rel L_inf = {:.3e} (iters={}, ffts={})",
                    rel_linf,
                    report.iterations,
                    report.fft_count
                );
            });
        }
    }

    // All four solver arms on the smooth cubic MMS problem (24^3, eps_b = 10.0 so
    // non-polynomial sqrt(eps) / ln(eps) harmonics are well resolved alongside exact DirectPcg
    // and HybridFisicaroDirect).
    let cubic = cubic_cell(10.0);
    let dims = [24, 24, 24];
    for with_screening in [false, true] {
        dft_mpi::simulate(1, |comm| {
            let mut ws = SolvationWorkspace::new(&cubic, dims, comm).unwrap();
            let owned = ws.owned_real_grid_indices();
            let (eps, kappa2, phi_exact, _) =
                build_mms_fields(&cubic, dims, &owned, 10.0, with_screening);

            let n_g = ws.local_positions().len();
            let n_r = ws.plane_points();
            let mut phi_exact_g = vec![Complex64::new(0.0, 0.0); n_g];
            let mut rhs_g = vec![Complex64::new(0.0, 0.0); n_g];
            ws.forward_real(comm, &phi_exact, &mut phi_exact_g).unwrap();
            ws.apply_operator(comm, &eps, &kappa2, &phi_exact_g, &mut rhs_g)
                .unwrap();

            for kind in [
                SolverKind::DirectPcg,
                SolverKind::FisicaroPcg,
                SolverKind::PolarizationDiis,
                SolverKind::HybridFisicaroDirect,
            ] {
                ws.reset_warm_start();
                let cfg = SolverConfig {
                    kind,
                    tolerance: 1.0e-11,
                    max_iterations: 400,
                    diis_history: 12,
                    warm_start: false,
                };
                let mut phi_num_g = vec![Complex64::new(0.0, 0.0); n_g];
                let report = ws
                    .solve_pde(comm, &eps, &kappa2, &rhs_g, &cfg, &mut phi_num_g)
                    .unwrap();
                assert!(
                    report.converged,
                    "{kind:?} (screen={with_screening}) did not converge: {report:?}"
                );
                let mut phi_num_r = vec![0.0_f64; n_r];
                ws.inverse_to_real(comm, &phi_num_g, &mut phi_num_r)
                    .unwrap();
                let rel_linf = max_abs_diff(&phi_num_r, &phi_exact) / max_abs(&phi_exact);
                println!(
                    "[MMS cubic24 eps_b=10 screen={with_screening}] {kind:?}: iters={} ffts={} rel_linf={rel_linf:.3e}",
                    report.iterations, report.fft_count
                );
                // DirectPcg and HybridFisicaroDirect solve the exact discrete operator to < 1e-10;
                // FisicaroPcg and PolarizationDiis solve continuum-equivalent non-polynomial
                // reformulations (sqrt(eps) / ln(eps)) with finite-grid product-rule aliasing.
                let tol = match kind {
                    SolverKind::DirectPcg | SolverKind::HybridFisicaroDirect => 1.0e-10,
                    SolverKind::FisicaroPcg | SolverKind::PolarizationDiis => 2.0e-5,
                };
                assert!(
                    rel_linf < tol,
                    "{kind:?} (screen={with_screening}) rel L_inf = {rel_linf:.3e} (tol {tol:.1e})"
                );
            }
        });
    }
}

/// Computes the exact continuum 1D radial integrals for a spherical Gaussian
/// charge `q / (pi a^2)^(3/2) exp(-r^2/a^2)` in a smooth spherical dielectric
/// `eps(r) = 1 + (eps_b - 1) * 0.5 * (1 + tanh((r - R)/w))`:
/// returns `(delta_e_born_smooth_ha, quad_moment_solv, r_eff_bohr)`.
fn exact_smooth_born_continuum(
    q: f64,
    a: f64,
    r_cav: f64,
    width: f64,
    eps_b: f64,
) -> (f64, f64, f64) {
    let steps = 20_000usize;
    let r_max = 40.0_f64;
    let dr = r_max / steps as f64;
    let mut born_int = 0.0_f64;
    let mut quad_solv = 0.0_f64;
    let inv_sqrt_pi = 1.0 / PI.sqrt();

    for i in 0..steps {
        let r = (i as f64 + 0.5) * dr;
        let x = r / a;
        let q_enc = q * (libm::erf(x) - 2.0 * x * inv_sqrt_pi * (-x * x).exp());
        let s = 0.5 * (1.0 + ((r - r_cav) / width).tanh());
        let eps = 1.0 + (eps_b - 1.0) * s;
        born_int += (q_enc * q_enc / (r * r)) * (1.0 - 1.0 / eps) * dr;
        quad_solv += (4.0 * PI / 3.0) * r * (q_enc / eps - q / eps_b) * dr;
    }
    let tail = 0.5 * q * q / r_max * (1.0 - 1.0 / eps_b);
    let delta_e_smooth = -0.5 * born_int - tail;
    let r_eff = -q * q / (2.0 * delta_e_smooth) * (1.0 - 1.0 / eps_b);
    (delta_e_smooth, quad_solv, r_eff)
}

/// Exact 1D radial finite-volume LPB continuum solver on `r in [0, 40]` Bohr for
/// `-1/r^2 d/dr(r^2 eps(r) dphi/dr) + eps_b * kappa^2 * S(r) * phi = 4*pi*rho(r)`
/// with Robin asymptotic boundary `d(r phi)/dr = -kappa (r phi)` at `r_max`.
/// Returns isolated electrostatic solvation energy `0.5 * int rho * (phi_solv - phi_vac) dV`.
fn exact_smooth_lpb_continuum(model: BornModel, kappa: f64) -> f64 {
    let n = 4000usize;
    let r_max = 40.0_f64;
    let h = r_max / n as f64;
    let pre = model.q / (PI * model.a * model.a).powf(1.5);
    let inv_a2 = 1.0 / (model.a * model.a);

    let solve_1d = |eps_b: f64, kap: f64| -> Vec<f64> {
        let mut lower = vec![0.0_f64; n];
        let mut diag = vec![0.0_f64; n];
        let mut upper = vec![0.0_f64; n];
        let mut rhs = vec![0.0_f64; n];
        for i in 0..n {
            let r = (i as f64 + 0.5) * h;
            let r_lo = i as f64 * h;
            let r_hi = (i as f64 + 1.0) * h;
            let s_c = 0.5 * (1.0 + ((r - model.r_cav) / model.width).tanh());
            let s_lo = 0.5 * (1.0 + ((r_lo - model.r_cav) / model.width).tanh());
            let s_hi = 0.5 * (1.0 + ((r_hi - model.r_cav) / model.width).tanh());
            let e_lo = 1.0 + (eps_b - 1.0) * s_lo;
            let e_hi = 1.0 + (eps_b - 1.0) * s_hi;
            let c_lo = if i == 0 {
                0.0
            } else {
                r_lo * r_lo * e_lo / (h * h)
            };
            let c_hi = r_hi * r_hi * e_hi / (h * h);
            let vol = r * r;
            let k2_term = vol * eps_b * kap * kap * s_c;
            lower[i] = -c_lo;
            if i + 1 < n {
                upper[i] = -c_hi;
                diag[i] = c_lo + c_hi + k2_term;
            } else {
                // Far-field Yukawa/Coulomb decay `phi(r) ~ C exp(-kap r) / r`:
                // `(phi_out - phi_last)/h = -(kap + 1/r_hi) * phi_last`.
                let bc = c_hi * (h * (kap + 1.0 / r_hi)) / (1.0 + 0.5 * h * (kap + 1.0 / r_hi));
                diag[i] = c_lo + bc + k2_term;
            }
            rhs[i] = 4.0 * PI * vol * pre * (-r * r * inv_a2).exp();
        }
        // Thomas tridiagonal solve.
        for i in 1..n {
            let w = lower[i] / diag[i - 1];
            diag[i] -= w * upper[i - 1];
            rhs[i] -= w * rhs[i - 1];
        }
        let mut phi = vec![0.0_f64; n];
        phi[n - 1] = rhs[n - 1] / diag[n - 1];
        for i in (0..n - 1).rev() {
            phi[i] = (rhs[i] - upper[i] * phi[i + 1]) / diag[i];
        }
        phi
    };

    let phi_solv = solve_1d(model.eps_b, kappa);
    let phi_vac = solve_1d(1.0, 0.0);
    let mut energy = 0.0_f64;
    for i in 0..n {
        let r = (i as f64 + 0.5) * h;
        let rho_r = pre * (-r * r * inv_a2).exp();
        energy += 0.5 * 4.0 * PI * r * r * h * rho_r * (phi_solv[i] - phi_vac[i]);
    }
    energy
}

#[derive(Clone, Copy)]
struct BornModel {
    q: f64,
    a: f64,
    r_cav: f64,
    width: f64,
    eps_b: f64,
}

fn solve_born_box(length_bohr: f64, n_axis: usize, model: BornModel, kappa: f64) -> f64 {
    let cell = cubic_cell(length_bohr);
    let dims = [n_axis, n_axis, n_axis];
    let dv = cell.volume_bohr3() / (n_axis * n_axis * n_axis) as f64;

    let res = dft_mpi::simulate(1, |comm| {
        let mut ws = SolvationWorkspace::new(&cell, dims, comm).unwrap();
        let owned = ws.owned_real_grid_indices();
        let n_r = ws.plane_points();
        let n_g = ws.local_positions().len();

        let mut rho = Vec::with_capacity(n_r);
        let mut eps = Vec::with_capacity(n_r);
        let mut kappa2 = Vec::with_capacity(n_r);
        let inv_a2 = 1.0 / (model.a * model.a);
        let pre = model.q / (PI * model.a * model.a).powf(1.5);

        for &idx in &owned {
            let u = grid_frac(idx, dims);
            let du = [
                u[0] - 0.5 - (u[0] - 0.5).round(),
                u[1] - 0.5 - (u[1] - 0.5).round(),
                u[2] - 0.5 - (u[2] - 0.5).round(),
            ];
            let dr = cell.fractional_to_cartesian_bohr(du);
            let r2 = dot3(dr, dr);
            let r = r2.sqrt();
            rho.push(pre * (-r2 * inv_a2).exp());
            let s = 0.5 * (1.0 + ((r - model.r_cav) / model.width).tanh());
            eps.push(1.0 + (model.eps_b - 1.0) * s);
            kappa2.push(model.eps_b * kappa * kappa * s);
        }

        // Normalize discrete Gaussian charge to exact q.
        let q_num: f64 = rho.iter().sum::<f64>() * dv;
        let q_scale = model.q / q_num;
        for v in &mut rho {
            *v *= q_scale;
        }

        let mut rhs_g = vec![Complex64::new(0.0, 0.0); n_g];
        let mut rho_g = vec![Complex64::new(0.0, 0.0); n_g];
        ws.forward_real(comm, &rho, &mut rho_g).unwrap();
        for (dst, &src) in rhs_g.iter_mut().zip(&rho_g) {
            *dst = src * (4.0 * PI);
        }

        let cfg = SolverConfig {
            kind: SolverKind::DirectPcg,
            tolerance: 1.0e-12,
            max_iterations: 400,
            diis_history: 8,
            warm_start: false,
        };

        // Vacuum potential (eps = 1, kappa = 0) with the same discrete operator.
        let ones = vec![1.0_f64; n_r];
        let zeros = vec![0.0_f64; n_r];
        let mut phi_vac_g = vec![Complex64::new(0.0, 0.0); n_g];
        let rep_vac = ws
            .solve_pde(comm, &ones, &zeros, &rhs_g, &cfg, &mut phi_vac_g)
            .unwrap();
        assert!(rep_vac.converged, "vacuum Born solve failed: {rep_vac:?}");

        ws.reset_warm_start();
        let mut phi_solv_g = vec![Complex64::new(0.0, 0.0); n_g];
        let rep_solv = ws
            .solve_pde(comm, &eps, &kappa2, &rhs_g, &cfg, &mut phi_solv_g)
            .unwrap();
        assert!(
            rep_solv.converged,
            "solvated Born solve failed: {rep_solv:?}"
        );

        let mut phi_vac_r = vec![0.0_f64; n_r];
        let mut phi_solv_r = vec![0.0_f64; n_r];
        ws.inverse_to_real(comm, &phi_vac_g, &mut phi_vac_r)
            .unwrap();
        ws.inverse_to_real(comm, &phi_solv_g, &mut phi_solv_r)
            .unwrap();

        rho.iter()
            .zip(&phi_solv_r)
            .zip(&phi_vac_r)
            .map(|((&r, &ps), &pv)| 0.5 * r * (ps - pv) * dv)
            .sum()
    });
    res[0]
}

/// 3. Analytic Born ion & Debye-Huckel screening convergence across box sizes (L1).
#[test]
fn born_ion_and_debye_huckel_analytic_convergence() {
    let model = BornModel {
        q: 1.0,
        a: 0.7,
        r_cav: 2.5,
        width: 0.45,
        eps_b: 78.4,
    };
    let BornModel {
        q,
        a,
        r_cav,
        width,
        eps_b,
    } = model;

    let delta_g_born_sharp = -q * q / (2.0 * r_cav) * (1.0 - 1.0 / eps_b);
    let (delta_e_smooth_exact, quad_solv, r_eff) =
        exact_smooth_born_continuum(q, a, r_cav, width, eps_b);
    let quad_vac = -PI * q * a * a;

    let l1 = 14.0;
    let l2 = 20.0;
    let de_l1_k0 = solve_born_box(l1, 24, model, 0.0);
    let de_l2_k0 = solve_born_box(l2, 32, model, 0.0);

    // Monotonic convergence toward the negative isolated Born energy as L grows.
    assert!(
        de_l2_k0 < de_l1_k0 && de_l1_k0 < 0.0,
        "expected monotonic convergence: de(L={l1})={de_l1_k0:.6}, de(L={l2})={de_l2_k0:.6}"
    );

    // Exact Makov-Payne (1/L Madelung + 1/L^3 quadrupole) periodic zero-mode
    // correction removes finite-box image shifts for the charged ion.
    let shift_solv =
        |l: f64| -MADELUNG_SC * q * q / (2.0 * eps_b * l) - 0.5 * q * quad_solv / (l * l * l);
    let shift_vac = |l: f64| -MADELUNG_SC * q * q / (2.0 * l) - 0.5 * q * quad_vac / (l * l * l);

    let de_inf_l1 = de_l1_k0 - (shift_solv(l1) - shift_vac(l1));
    let de_inf_l2 = de_l2_k0 - (shift_solv(l2) - shift_vac(l2));
    let rel_err_l2 = (de_inf_l2 - delta_e_smooth_exact).abs() / delta_e_smooth_exact.abs();

    println!(
        "[Born] sharp(R={r_cav})={delta_g_born_sharp:.6}, smooth_exact(R_eff={r_eff:.3})={delta_e_smooth_exact:.6}, \
         raw(L14)={de_l1_k0:.6}, raw(L20)={de_l2_k0:.6}, corr(L14)={de_inf_l1:.6}, corr(L20)={de_inf_l2:.6} (rel_err={rel_err_l2:.3e})"
    );
    assert!(
        rel_err_l2 < 0.01,
        "Madelung-corrected smooth Born energy {de_inf_l2:.6} deviates from continuum {delta_e_smooth_exact:.6} by {:.2}%",
        100.0 * rel_err_l2
    );
    // Linear dielectric switching shifts the effective Born radius inward by
    // ~ (w / 2) * ln(eps_b), so the smooth-wall energy is strictly more negative
    // than the sharp-wall R = 2.5 Born value.
    assert!(
        de_inf_l2 < delta_g_born_sharp,
        "expected smooth-wall Born energy {de_inf_l2:.6} < sharp-wall {delta_g_born_sharp:.6}"
    );

    // Debye-Huckel screening verified against the exact 1D smooth-wall LPB continuum ODE
    // across two box sizes (L = 24, N = 48 -> L = 32, N = 64; fixed dx = 0.50 Bohr) to confirm
    // monotonic convergence toward < 1% relative error (pre-declared L2 tolerance):
    // (a) At eps_b = 78.4, compare kappa = 0.4 vs kappa = 0.2 directly so both
    //     phi_vac and the kappa = 0 periodic zero-mode offset cancel identically.
    let dh_inc_exact =
        exact_smooth_lpb_continuum(model, 0.4) - exact_smooth_lpb_continuum(model, 0.2);
    let dh_inc_l24 = solve_born_box(24.0, 48, model, 0.4) - solve_born_box(24.0, 48, model, 0.2);
    let de_l32_k02 = solve_born_box(32.0, 64, model, 0.2);
    let de_l32_k04 = solve_born_box(32.0, 64, model, 0.4);
    let dh_inc_l32 = de_l32_k04 - de_l32_k02;
    let rel_dh_l24 = (dh_inc_l24 - dh_inc_exact).abs() / dh_inc_exact.abs();
    let rel_dh_l32 = (dh_inc_l32 - dh_inc_exact).abs() / dh_inc_exact.abs();
    println!(
        "[Debye-Huckel eps_b={eps_b} k=0.4-0.2] 1d_exact={dh_inc_exact:.6}, L24={dh_inc_l24:.6} ({rel_dh_l24:.3e}), L32={dh_inc_l32:.6} ({rel_dh_l32:.3e})"
    );
    assert!(
        de_l32_k04 < de_l32_k02,
        "kappa=0.4 must screen more strongly than kappa=0.2: {de_l32_k04:.6} vs {de_l32_k02:.6}"
    );
    assert!(
        rel_dh_l32 < rel_dh_l24 && rel_dh_l32 < 0.01,
        "Debye-Huckel increment (0.4 vs 0.2) failed monotonic < 1% convergence: L24 err={rel_dh_l24:.3e}, L32 err={rel_dh_l32:.3e}"
    );

    // (b) At eps_b = 10.0, verify kappa = 0.3 extra screening against the exact 1D smooth LPB ODE
    //     (kappa = 0.3 keeps periodic Yukawa image tails e^{-0.3 * 32} = 6.8e-5 negligible).
    let model_10 = BornModel {
        eps_b: 10.0,
        ..model
    };
    let (_, quad_solv_10, _) = exact_smooth_born_continuum(q, a, r_cav, width, model_10.eps_b);
    let dh10_exact =
        exact_smooth_lpb_continuum(model_10, 0.3) - exact_smooth_lpb_continuum(model_10, 0.0);
    let dh10_at = |l_box: f64, n_box: usize| {
        let shift_solv_10 = -MADELUNG_SC * q * q / (2.0 * model_10.eps_b * l_box)
            - 0.5 * q * quad_solv_10 / (l_box * l_box * l_box);
        let de10_k0 = solve_born_box(l_box, n_box, model_10, 0.0);
        let de10_k03 = solve_born_box(l_box, n_box, model_10, 0.3);
        let de10_inf_k0 = de10_k0 - (shift_solv_10 - shift_vac(l_box));
        let de10_inf_k03 = de10_k03 + shift_vac(l_box);
        (de10_inf_k03, de10_inf_k0, de10_inf_k03 - de10_inf_k0)
    };
    let (_, _, dh10_l24) = dh10_at(24.0, 48);
    let (de10_inf_k03, de10_inf_k0, dh10_l36) = dh10_at(36.0, 72);
    let rel_dh10_l24 = (dh10_l24 - dh10_exact).abs() / dh10_exact.abs();
    let rel_dh10_l36 = (dh10_l36 - dh10_exact).abs() / dh10_exact.abs();
    println!(
        "[Debye-Huckel eps_b=10.0 k=0.3 vs 0] 1d_exact={dh10_exact:.6}, L24={dh10_l24:.6} ({rel_dh10_l24:.3e}), L36={dh10_l36:.6} ({rel_dh10_l36:.3e})"
    );
    assert!(
        de10_inf_k03 < de10_inf_k0,
        "eps_b=10 kappa=0.3 must screen more strongly than kappa=0: {de10_inf_k03:.6} vs {de10_inf_k0:.6}"
    );
    assert!(
        rel_dh10_l36 < rel_dh10_l24 && rel_dh10_l36 < 0.01,
        "eps_b=10 kappa=0.3 Debye-Huckel failed monotonic < 1% convergence: L24 err={rel_dh10_l24:.3e}, L36 err={rel_dh10_l36:.3e}"
    );
}

/// Picks 5 distinct owned real-space indices spanning inside the solute cavity,
/// three points across the transition wall (safely away from SCCS piecewise
/// kinks `rho_min = 1e-4` and `rho_max = 1.5e-3`), and the bulk solvent.
fn pick_variational_probe_indices(density: &[f64]) -> Vec<usize> {
    let buckets: [(f64, f64); 5] = [
        (2.0e-2, f64::INFINITY), // deep inside cavity
        (2.0e-3, 6.0e-3),        // upper transition wall (VASPsol active)
        (6.0e-4, 1.2e-3),        // mid transition wall (VASPsol + SCCS active)
        (2.0e-4, 5.0e-4),        // lower transition wall (VASPsol + SCCS active)
        (0.0, 5.0e-5),           // bulk solvent
    ];
    let mut chosen = Vec::with_capacity(5);
    for (lo, hi) in buckets {
        if let Some((idx, _)) = density
            .iter()
            .enumerate()
            .find(|&(i, &val)| val > lo && val < hi && !chosen.contains(&i))
        {
            chosen.push(idx);
        }
    }
    for i in 0..density.len() {
        if chosen.len() == 5 {
            break;
        }
        if !chosen.contains(&i) {
            chosen.push(i);
        }
    }
    assert_eq!(chosen.len(), 5);
    chosen
}

/// Test 4: Variational potential `V_solv = dE_solv / dn` vs central finite differences
/// on a skewed cell for VASPsol, SCCS (log-eps), SCCS (linear-eps), and pure
/// cavitation isolation (L2).
#[test]
fn variational_potential_matches_central_finite_differences_vaspsol_and_sccs() {
    let cell =
        Cell::from_user_rows_bohr([[9.6, 0.7, -0.3], [0.4, 10.2, 0.6], [-0.2, 0.5, 9.8]]).unwrap();
    let dims = [12, 14, 12];
    let dv = cell.volume_bohr3() / (dims[0] * dims[1] * dims[2]) as f64;

    let ions = [IonicPseudoCharge {
        charge: 1.0,
        fractional_position: [0.47, 0.51, 0.49],
        width_bohr: 0.65,
    }];
    // Slightly asymmetric electron density (`N_e = 1.1 != Q_ion = 1.0`) so every
    // variational term (including `E_bg,corr` at `kappa > 0`, `V_diel`, `V_scr`,
    // and `V_cav`) is simultaneously non-trivial.
    let global_n = build_global_gaussian_density(
        &cell,
        dims,
        &[
            ([0.47, 0.51, 0.49], 0.85, 1.00),
            ([0.57, 0.45, 0.54], 0.25, 1.00),
        ],
        1.10,
        2.0e-6,
    );

    let cavities = [
        ("vaspsol_water", CavityModel::vaspsol_water()),
        ("sccs_log", CavityModel::sccs_water()),
        (
            "sccs_linear",
            CavityModel::Sccs {
                rho_min_bohr3: 1.0e-4,
                rho_max_bohr3: 1.5e-3,
                log_epsilon: false,
            },
        ),
    ];

    let cfg = SolverConfig {
        kind: SolverKind::DirectPcg,
        tolerance: 1.0e-14,
        max_iterations: 400,
        diis_history: 8,
        warm_start: false,
    };

    dft_mpi::simulate(1, |comm| {
        let mut ws = SolvationWorkspace::new(&cell, dims, comm).unwrap();
        let base_n = slice_owned(&global_n, &ws.owned_real_grid_indices());
        let probes = pick_variational_probe_indices(&base_n);

        for (label, cavity) in cavities {
            let params = SolvationParams {
                epsilon_bulk: 78.4,
                kappa_inv_bohr: 0.25,
                tau_ha_bohr2: 5.0e-5,
                beta_ha_bohr3: 2.0e-5,
                cavity,
                ..SolvationParams::vaspsol_water()
            };

            let eval0 = ws.evaluate(comm, &base_n, &ions, &params, &cfg).unwrap();
            assert!(
                eval0.report.converged,
                "{label} base solve did not converge"
            );
            assert!(
                max_abs(&eval0.v_diel_planes) > 1.0e-5
                    && max_abs(&eval0.v_scr_planes) > 1.0e-7
                    && max_abs(&eval0.v_cav_planes) > 1.0e-6,
                "{label}: expected non-trivial dielectric/screening/cavitation potentials"
            );

            let mut perturbed = base_n.clone();
            for &idx in &probes {
                let n_i = base_n[idx];
                // 4th-order Richardson central difference:
                //   (-E(+2h) + 8*E(+h) - 8*E(-h) + E(-2h)) / (12 * h * dv)
                // keeps both truncation O((h/n_i)^4) and f64 rounding below 1e-11 Ha
                // while staying strictly inside SCCS piecewise-smooth intervals.
                let h = if n_i < 1.0e-5 {
                    0.20 * n_i
                } else {
                    match cavity {
                        CavityModel::VaspSol { .. } => (1.5e-3 * n_i).min(2.0e-6),
                        CavityModel::Sccs { .. } => (8.0e-4 * n_i).min(2.0e-6),
                    }
                };
                assert!(n_i - 2.0 * h > params.density_floor_bohr3);

                perturbed[idx] = n_i + 2.0 * h;
                let ep2 = ws.evaluate(comm, &perturbed, &ions, &params, &cfg).unwrap();
                perturbed[idx] = n_i + h;
                let ep1 = ws.evaluate(comm, &perturbed, &ions, &params, &cfg).unwrap();
                perturbed[idx] = n_i - h;
                let em1 = ws.evaluate(comm, &perturbed, &ions, &params, &cfg).unwrap();
                perturbed[idx] = n_i - 2.0 * h;
                let em2 = ws.evaluate(comm, &perturbed, &ions, &params, &cfg).unwrap();
                perturbed[idx] = n_i;

                let v_fd = (-ep2.total_ha + 8.0 * ep1.total_ha - 8.0 * em1.total_ha + em2.total_ha)
                    / (12.0 * h * dv);
                let v_ana = eval0.v_total_planes[idx];
                let abs_err = (v_fd - v_ana).abs();
                let tol = 1.0e-10_f64.max(1.0e-9 * v_ana.abs());
                assert!(
                    abs_err < tol,
                    "{label} idx={idx} (n={n_i:.3e}, h={h:.3e}): v_ana={v_ana:.12e}, v_fd={v_fd:.12e}, abs_err={abs_err:.3e} (tol={tol:.3e})"
                );
            }
        }

        // Term-by-term variational FD isolation against pre-declared tolerances
        // (rel < 1e-9, abs < 1e-10 Ha):
        // (1) Pure cavitation (eps_b = 1, kappa = 0, tau > 0, beta > 0) -> v_cav alone.
        // (2) Pure dielectric + electrostatic at kappa = 0 (tau = 0, beta = 0) -> v_es + v_diel
        //     (also verifies kappa = 0 background_correction_ha gauge consistency!).
        // (3) Pure screening + electrostatic at eps_b = 1, kappa > 0 (deps_dn = 0, tau = 0) -> v_es + v_scr.
        // (4) Pure electrostatic v_es alone on a deep-interior / bulk point where deps_dn = dkappa2_dn = 0
        //     (SCCS compact support guarantees v_diel = v_scr = v_cav = 0 identically on interior/bulk!).
        let iso_cases = [
            (
                "cavitation_only",
                SolvationParams {
                    epsilon_bulk: 1.0,
                    kappa_inv_bohr: 0.0,
                    tau_ha_bohr2: 5.0e-5,
                    beta_ha_bohr3: 2.0e-5,
                    cavity: CavityModel::vaspsol_water(),
                    ..SolvationParams::vacuum()
                },
            ),
            (
                "diel_and_es_kappa0",
                SolvationParams {
                    epsilon_bulk: 78.4,
                    kappa_inv_bohr: 0.0,
                    tau_ha_bohr2: 0.0,
                    beta_ha_bohr3: 0.0,
                    cavity: CavityModel::sccs_water(),
                    ..SolvationParams::vaspsol_water()
                },
            ),
            (
                "scr_and_es_eps1",
                SolvationParams {
                    epsilon_bulk: 1.0,
                    kappa_inv_bohr: 0.25,
                    tau_ha_bohr2: 0.0,
                    beta_ha_bohr3: 0.0,
                    cavity: CavityModel::sccs_water(),
                    ..SolvationParams::vacuum()
                },
            ),
        ];

        for (iso_label, iso_params) in iso_cases {
            let eval_iso = ws
                .evaluate(comm, &base_n, &ions, &iso_params, &cfg)
                .unwrap();
            assert!(eval_iso.report.converged);
            match iso_label {
                "cavitation_only" => {
                    assert!(max_abs(&eval_iso.v_diel_planes) < 1.0e-14);
                    assert!(max_abs(&eval_iso.v_scr_planes) < 1.0e-14);
                    assert!(max_abs(&eval_iso.v_cav_planes) > 1.0e-6);
                }
                "diel_and_es_kappa0" => {
                    assert!(max_abs(&eval_iso.v_scr_planes) < 1.0e-14);
                    assert!(max_abs(&eval_iso.v_cav_planes) < 1.0e-14);
                    assert!(max_abs(&eval_iso.v_diel_planes) > 1.0e-5);
                    assert!(max_abs(&eval_iso.v_es_planes) > 1.0e-3);
                }
                "scr_and_es_eps1" => {
                    assert!(max_abs(&eval_iso.v_diel_planes) < 1.0e-14);
                    assert!(max_abs(&eval_iso.v_cav_planes) < 1.0e-14);
                    assert!(max_abs(&eval_iso.v_scr_planes) > 1.0e-7);
                    assert!(max_abs(&eval_iso.v_es_planes) > 1.0e-3);
                }
                _ => unreachable!(),
            }

            let mut perturbed = base_n.clone();
            let mut max_err = 0.0_f64;
            for &idx in &probes {
                let n_i = base_n[idx];
                let h = if n_i < 1.0e-5 {
                    0.20 * n_i
                } else {
                    match iso_params.cavity {
                        CavityModel::VaspSol { .. } => (1.5e-3 * n_i).min(2.0e-6),
                        CavityModel::Sccs { .. } => (8.0e-4 * n_i).min(2.0e-6),
                    }
                };
                assert!(n_i - 2.0 * h > iso_params.density_floor_bohr3);

                perturbed[idx] = n_i + 2.0 * h;
                let ep2 = ws
                    .evaluate(comm, &perturbed, &ions, &iso_params, &cfg)
                    .unwrap();
                perturbed[idx] = n_i + h;
                let ep1 = ws
                    .evaluate(comm, &perturbed, &ions, &iso_params, &cfg)
                    .unwrap();
                perturbed[idx] = n_i - h;
                let em1 = ws
                    .evaluate(comm, &perturbed, &ions, &iso_params, &cfg)
                    .unwrap();
                perturbed[idx] = n_i - 2.0 * h;
                let em2 = ws
                    .evaluate(comm, &perturbed, &ions, &iso_params, &cfg)
                    .unwrap();
                perturbed[idx] = n_i;

                let (v_fd, v_ana) = if iso_label == "cavitation_only" {
                    (
                        (-ep2.cavitation_ha + 8.0 * ep1.cavitation_ha - 8.0 * em1.cavitation_ha
                            + em2.cavitation_ha)
                            / (12.0 * h * dv),
                        eval_iso.v_cav_planes[idx],
                    )
                } else {
                    (
                        (-ep2.total_ha + 8.0 * ep1.total_ha - 8.0 * em1.total_ha + em2.total_ha)
                            / (12.0 * h * dv),
                        eval_iso.v_total_planes[idx],
                    )
                };
                let abs_err = (v_fd - v_ana).abs();
                max_err = max_err.max(abs_err);
                let tol = 1.0e-10_f64.max(1.0e-9 * v_ana.abs());
                assert!(
                    abs_err < tol,
                    "{iso_label} idx={idx}: v_ana={v_ana:.12e}, v_fd={v_fd:.12e}, abs_err={abs_err:.3e} (tol={tol:.3e})"
                );
            }
            println!("[VariationalFD {iso_label}] max_abs_err={max_err:.3e} Ha");
        }
    });
}

/// Test 5: Explicit Gaussian ionic pseudo-charge forces vs central finite differences
/// on a skewed cell, plus `r_c` width invariance when the ionic core lies well
/// inside the solute cavity (L2).
#[test]
fn ionic_forces_match_central_finite_differences_and_width_invariance() {
    // Part A: Central FD force check along Cartesian x, y, z on a skewed cell.
    let cell = Cell::from_user_rows_bohr([[10.2, 0.7, -0.3], [0.4, 10.6, 0.5], [-0.3, 0.4, 10.0]])
        .unwrap();
    let dims = [14, 14, 14];
    let base_ions = [
        IonicPseudoCharge {
            charge: 1.5,
            fractional_position: [0.42, 0.47, 0.51],
            width_bohr: 0.65,
        },
        IonicPseudoCharge {
            charge: 1.0,
            fractional_position: [0.59, 0.54, 0.46],
            width_bohr: 0.60,
        },
    ];
    let q_ion: f64 = base_ions.iter().map(|i| i.charge).sum();
    let global_n = build_global_gaussian_density(
        &cell,
        dims,
        &[
            ([0.41, 0.48, 0.50], 1.4, 1.15),
            ([0.60, 0.53, 0.47], 1.1, 1.10),
        ],
        q_ion,
        1.0e-5,
    );
    let params = SolvationParams {
        epsilon_bulk: 78.4,
        kappa_inv_bohr: 0.2,
        tau_ha_bohr2: 2.0e-5,
        ..SolvationParams::vaspsol_water()
    };
    let cfg = SolverConfig {
        kind: SolverKind::DirectPcg,
        tolerance: 1.0e-14,
        max_iterations: 400,
        diis_history: 8,
        warm_start: false,
    };

    dft_mpi::simulate(1, |comm| {
        let mut ws = SolvationWorkspace::new(&cell, dims, comm).unwrap();
        let local_n = slice_owned(&global_n, &ws.owned_real_grid_indices());
        let eval0 = ws
            .evaluate(comm, &local_n, &base_ions, &params, &cfg)
            .unwrap();
        assert!(eval0.report.converged);
        let f_ana = eval0.ionic_forces_ha_bohr[0];
        let f_norm = dot3(f_ana, f_ana).sqrt();
        assert!(
            f_norm > 1.0e-5,
            "expected non-trivial ionic force, got {f_ana:?}"
        );

        let h = 1.0e-5;
        for axis in 0..3 {
            let mut shift_cart = [0.0_f64; 3];
            shift_cart[axis] = h;
            let dfrac = cell.cartesian_bohr_to_fractional(shift_cart);

            let mut ions_p = base_ions;
            let mut ions_m = base_ions;
            for ((p, m), &d) in ions_p[0]
                .fractional_position
                .iter_mut()
                .zip(&mut ions_m[0].fractional_position)
                .zip(&dfrac)
            {
                *p += d;
                *m -= d;
            }

            let ep = ws.evaluate(comm, &local_n, &ions_p, &params, &cfg).unwrap();
            let em = ws.evaluate(comm, &local_n, &ions_m, &params, &cfg).unwrap();
            let f_fd = -(ep.total_ha - em.total_ha) / (2.0 * h);
            let err = (f_fd - f_ana[axis]).abs();
            assert!(
                err < 1.0e-9,
                "axis={axis}: f_ana={:.12e}, f_fd={f_fd:.12e}, err={err:.3e}",
                f_ana[axis]
            );
        }
    });

    // Part B: Gaussian ionic width invariance (`r_c = 0.65` vs `0.80` Bohr) on BOTH
    // a neutral dipolar cell (kappa = 0) AND a charged ion cell (Q = 5, N_e = 4, kappa = 0.2,
    // verifying `E_es + E_bg,corr` cancellation) to the pre-declared `< 1e-6 Ha` tolerance.
    // Grid spacing `dx = 0.333 Bohr` (`r_c / dx >= 1.95`) resolves both Gaussians to `< 1e-8`,
    // while `sigma = 1.45 Bohr` places the SCCS cavity edge well outside both ionic cores.
    let w_cell = cubic_cell(12.0);
    let w_dims = [36, 36, 36];
    let w_global_n = build_global_gaussian_density(
        &w_cell,
        w_dims,
        &[
            ([0.45, 0.50, 0.50], 2.6, 1.45),
            ([0.57, 0.50, 0.50], 1.4, 1.45),
        ],
        4.0,
        1.0e-5,
    );

    dft_mpi::simulate(1, |comm| {
        let mut ws = SolvationWorkspace::new(&w_cell, w_dims, comm).unwrap();
        let local_n = slice_owned(&w_global_n, &ws.owned_real_grid_indices());

        for (label, q_ion, kappa) in [
            ("neutral_k0", 4.0_f64, 0.0_f64),
            ("charged_q+1_k0.2", 5.0_f64, 0.2_f64),
        ] {
            let w_params = SolvationParams {
                epsilon_bulk: 78.4,
                kappa_inv_bohr: kappa,
                tau_ha_bohr2: 0.0,
                ..SolvationParams::sccs_water()
            };
            let ion_65 = [IonicPseudoCharge {
                charge: q_ion,
                fractional_position: [0.45, 0.50, 0.50],
                width_bohr: 0.65,
            }];
            let ion_80 = [IonicPseudoCharge {
                charge: q_ion,
                fractional_position: [0.45, 0.50, 0.50],
                width_bohr: 0.80,
            }];

            let e65 = ws
                .evaluate(comm, &local_n, &ion_65, &w_params, &cfg)
                .unwrap();
            let e80 = ws
                .evaluate(comm, &local_n, &ion_80, &w_params, &cfg)
                .unwrap();
            let diff = (e65.total_ha - e80.total_ha).abs();
            println!(
                "[WidthInvariance {label}] E(0.65)={:.10} (bg={:.6}), E(0.80)={:.10} (bg={:.6}), |dE|={diff:.3e} Ha",
                e65.total_ha,
                e65.background_correction_ha,
                e80.total_ha,
                e80.background_correction_ha
            );
            assert!(
                e65.total_ha.abs() > 1.0e-4,
                "{label}: expected non-trivial solvation energy, got {:.3e}",
                e65.total_ha
            );
            assert!(
                diff < 1.0e-6,
                "{label} width invariance violated: E(0.65)={:.10}, E(0.80)={:.10}, diff={diff:.3e} Ha",
                e65.total_ha,
                e80.total_ha
            );
        }
    });
}

struct ReassembledRun {
    total_ha: f64,
    electrostatic_ha: f64,
    cavitation_ha: f64,
    background_correction_ha: f64,
    forces: Vec<[f64; 3]>,
    v_total_global: Vec<f64>,
    phi_pol_global: Vec<f64>,
    mirror_paired: bool,
}

struct MpiTestProblem<'a> {
    cell: &'a Cell,
    dims: [usize; 3],
    global_n: &'a [f64],
    ions: &'a [IonicPseudoCharge],
    params: &'a SolvationParams,
    cfg: &'a SolverConfig,
}

fn run_mpi_eval(
    ranks: usize,
    mirror_override: Option<bool>,
    problem: &MpiTestProblem<'_>,
) -> ReassembledRun {
    let n_grid = problem.dims[0] * problem.dims[1] * problem.dims[2];
    let per_rank = dft_mpi::simulate(ranks, |comm| {
        let _guard = MirrorOverrideGuard::set(mirror_override);
        let mut ws = SolvationWorkspace::new(problem.cell, problem.dims, comm).unwrap();
        let owned = ws.owned_real_grid_indices();
        let local_n = slice_owned(problem.global_n, &owned);
        let eval: SolvationEvaluation = ws
            .evaluate(comm, &local_n, problem.ions, problem.params, problem.cfg)
            .unwrap();
        assert!(
            eval.report.converged,
            "rank {} failed to converge",
            comm.rank()
        );
        assert_eq!(ws.mirror_paired(), eval.report.mirror_paired);
        (owned, eval, ws.mirror_paired())
    });

    let mut v_total_global = vec![f64::NAN; n_grid];
    let mut phi_pol_global = vec![f64::NAN; n_grid];
    let mut seen = vec![false; n_grid];

    let ref_eval = &per_rank[0].1;
    let ref_mirror = per_rank[0].2;

    for (owned, eval, mirror_paired) in &per_rank {
        assert_eq!(*mirror_paired, ref_mirror);
        assert!((eval.total_ha - ref_eval.total_ha).abs() < 1.0e-13);
        assert!((eval.electrostatic_ha - ref_eval.electrostatic_ha).abs() < 1.0e-13);
        assert!((eval.cavitation_ha - ref_eval.cavitation_ha).abs() < 1.0e-13);
        assert!(
            (eval.background_correction_ha - ref_eval.background_correction_ha).abs() < 1.0e-13
        );
        assert!(
            max_force_diff(&eval.ionic_forces_ha_bohr, &ref_eval.ionic_forces_ha_bohr) < 1.0e-13
        );

        for (slot, &g_idx) in owned.iter().enumerate() {
            assert!(!seen[g_idx], "duplicate real-space owner for idx={g_idx}");
            seen[g_idx] = true;
            v_total_global[g_idx] = eval.v_total_planes[slot];
            phi_pol_global[g_idx] = eval.phi_pol_planes[slot];
        }
    }
    assert!(
        seen.iter().all(|&s| s),
        "missing owned real-space grid point"
    );

    ReassembledRun {
        total_ha: ref_eval.total_ha,
        electrostatic_ha: ref_eval.electrostatic_ha,
        cavitation_ha: ref_eval.cavitation_ha,
        background_correction_ha: ref_eval.background_correction_ha,
        forces: ref_eval.ionic_forces_ha_bohr.clone(),
        v_total_global,
        phi_pol_global,
        mirror_paired: ref_mirror,
    }
}

/// Test 6: MPI rank independence (`P = 1, 2, 4`) and `mirror` vs `greedy` stick
/// partition parity (L2).
#[test]
fn mpi_rank_independence_and_mirror_vs_greedy_parity() {
    let cell =
        Cell::from_user_rows_bohr([[10.2, 0.8, -0.3], [0.5, 10.8, 0.6], [-0.4, 0.5, 9.6]]).unwrap();
    let dims = [14, 16, 12];
    let ions = [
        IonicPseudoCharge {
            charge: 2.0,
            fractional_position: [0.40, 0.46, 0.52],
            width_bohr: 0.60,
        },
        IonicPseudoCharge {
            charge: 1.5,
            fractional_position: [0.60, 0.55, 0.45],
            width_bohr: 0.65,
        },
    ];
    let q_ion: f64 = ions.iter().map(|i| i.charge).sum();
    let global_n = build_global_gaussian_density(
        &cell,
        dims,
        &[
            ([0.40, 0.46, 0.52], 2.1, 1.15),
            ([0.60, 0.55, 0.45], 1.4, 1.10),
        ],
        q_ion,
        1.0e-5,
    );
    let params = SolvationParams {
        epsilon_bulk: 78.4,
        kappa_inv_bohr: 0.2,
        tau_ha_bohr2: 3.0e-5,
        beta_ha_bohr3: 1.0e-5,
        ..SolvationParams::vaspsol_water()
    };
    let cfg = SolverConfig {
        kind: SolverKind::DirectPcg,
        tolerance: 1.0e-14,
        max_iterations: 400,
        diis_history: 8,
        warm_start: false,
    };
    let problem = MpiTestProblem {
        cell: &cell,
        dims,
        global_n: &global_n,
        ions: &ions,
        params: &params,
        cfg: &cfg,
    };

    let p1 = run_mpi_eval(1, None, &problem);
    for ranks in [2, 4] {
        let pr = run_mpi_eval(ranks, None, &problem);
        assert!(
            (pr.total_ha - p1.total_ha).abs() < 1.0e-12,
            "P={ranks} total_ha diff = {:.3e}",
            (pr.total_ha - p1.total_ha).abs()
        );
        assert!(
            (pr.electrostatic_ha - p1.electrostatic_ha).abs() < 1.0e-12,
            "P={ranks} electrostatic_ha diff = {:.3e}",
            (pr.electrostatic_ha - p1.electrostatic_ha).abs()
        );
        assert!(
            (pr.cavitation_ha - p1.cavitation_ha).abs() < 1.0e-12,
            "P={ranks} cavitation_ha diff = {:.3e}",
            (pr.cavitation_ha - p1.cavitation_ha).abs()
        );
        assert!(
            (pr.background_correction_ha - p1.background_correction_ha).abs() < 1.0e-12,
            "P={ranks} background_correction_ha diff = {:.3e}",
            (pr.background_correction_ha - p1.background_correction_ha).abs()
        );
        assert!(
            max_force_diff(&pr.forces, &p1.forces) < 1.0e-12,
            "P={ranks} force diff = {:.3e}",
            max_force_diff(&pr.forces, &p1.forces)
        );
        assert!(
            max_abs_diff(&pr.v_total_global, &p1.v_total_global) < 1.0e-12,
            "P={ranks} v_total diff = {:.3e}",
            max_abs_diff(&pr.v_total_global, &p1.v_total_global)
        );
        assert!(
            max_abs_diff(&pr.phi_pol_global, &p1.phi_pol_global) < 1.0e-12,
            "P={ranks} phi_pol diff = {:.3e}",
            max_abs_diff(&pr.phi_pol_global, &p1.phi_pol_global)
        );
    }

    // Explicit greedy vs mirror toggle on P = 3 (where greedy round-robin genuinely
    // splits conjugate mirror sticks across ranks; guard restores override to None).
    let greedy = run_mpi_eval(3, Some(false), &problem);
    let mirror = run_mpi_eval(3, Some(true), &problem);
    assert!(
        !greedy.mirror_paired,
        "expected greedy partition on P=3 to report mirror_paired = false"
    );
    assert!(
        mirror.mirror_paired,
        "expected mirror partition on P=3 to report mirror_paired = true"
    );
    assert!(
        (greedy.total_ha - p1.total_ha).abs() < 1.0e-12,
        "P=3 greedy vs P=1 total_ha diff = {:.3e}",
        (greedy.total_ha - p1.total_ha).abs()
    );
    assert!(
        (greedy.total_ha - mirror.total_ha).abs() < 1.0e-12,
        "greedy vs mirror total_ha diff = {:.3e}",
        (greedy.total_ha - mirror.total_ha).abs()
    );
    assert!(
        max_force_diff(&greedy.forces, &mirror.forces) < 1.0e-12,
        "greedy vs mirror force diff = {:.3e}",
        max_force_diff(&greedy.forces, &mirror.forces)
    );
    assert!(
        max_abs_diff(&greedy.v_total_global, &mirror.v_total_global) < 1.0e-12,
        "greedy vs mirror v_total diff = {:.3e}",
        max_abs_diff(&greedy.v_total_global, &mirror.v_total_global)
    );
    assert!(
        max_abs_diff(&greedy.phi_pol_global, &mirror.phi_pol_global) < 1.0e-12,
        "greedy vs mirror phi_pol diff = {:.3e}",
        max_abs_diff(&greedy.phi_pol_global, &mirror.phi_pol_global)
    );
}

/// Test 7: Warm-start reproduces the cold solution while sharply cutting PCG
/// iterations on a `1e-3` density perturbation for both `DirectPcg` (exact discrete
/// variational baseline) and `FisicaroPcg` (fast 2-FFT/iter production solver) (L2/L5).
#[test]
fn warm_start_bulk_delta_reduces_iterations_sharply() {
    let cell = cubic_cell(12.0);
    // `dx = 0.25 Bohr` matches ADR 0061's FisicaroPcg resolution threshold so
    // Part B exercises FisicaroPcg directly without triggering the coarse-grid
    // DirectPcg breakdown fallback.
    let dims = [48, 48, 48];
    let ions = [
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
    ];
    let global_n0 = build_global_gaussian_density(
        &cell,
        dims,
        &[
            ([0.43, 0.49, 0.50], 2.2, 1.15),
            ([0.57, 0.51, 0.50], 1.8, 1.15),
        ],
        4.0,
        1.0e-5,
    );
    let global_n1: Vec<f64> = global_n0
        .iter()
        .enumerate()
        .map(|(idx, &val)| {
            let [u, v, w] = grid_frac(idx, dims);
            let mod_fn = 0.5 * (2.0 * PI * u).cos() - 0.3 * (2.0 * PI * v).sin()
                + 0.2 * (2.0 * PI * w).cos();
            val * (1.0 + 1.0e-3 * mod_fn)
        })
        .collect();

    let params = SolvationParams {
        epsilon_bulk: 78.4,
        kappa_inv_bohr: 0.15,
        tau_ha_bohr2: 2.0e-5,
        ..SolvationParams::vaspsol_water()
    };
    let cfg_warm = SolverConfig {
        kind: SolverKind::DirectPcg,
        tolerance: 1.0e-10,
        max_iterations: 300,
        diis_history: 8,
        warm_start: true,
    };
    let cfg_cold = SolverConfig {
        warm_start: false,
        ..cfg_warm
    };

    dft_mpi::simulate(1, |comm| {
        let mut ws_warm = SolvationWorkspace::new(&cell, dims, comm).unwrap();
        let mut ws_cold = SolvationWorkspace::new(&cell, dims, comm).unwrap();
        let owned = ws_warm.owned_real_grid_indices();
        let n0 = slice_owned(&global_n0, &owned);
        let n1 = slice_owned(&global_n1, &owned);

        // Prime workspace with baseline density (first solve is necessarily cold).
        let eval_base = ws_warm
            .evaluate(comm, &n0, &ions, &params, &cfg_warm)
            .unwrap();
        assert!(eval_base.report.converged);
        assert!(!eval_base.report.warm_started);

        // Part A: DirectPcg tight final-map tolerance (1e-10) verifies exact warm/cold parity & iteration drop.
        let eval_warm = ws_warm
            .evaluate(comm, &n1, &ions, &params, &cfg_warm)
            .unwrap();
        let eval_cold = ws_cold
            .evaluate(comm, &n1, &ions, &params, &cfg_cold)
            .unwrap();

        assert!(eval_warm.report.converged && eval_cold.report.converged);
        assert!(eval_warm.report.warm_started);
        assert!(!eval_cold.report.warm_started);

        let e_diff = (eval_warm.total_ha - eval_cold.total_ha).abs();
        println!(
            "[WarmStart DirectPcg 1e-3 @ tol=1e-10] cold_iters={} warm_iters={} cold_ffts={} warm_ffts={} |dE|={e_diff:.3e}",
            eval_cold.report.iterations,
            eval_warm.report.iterations,
            eval_cold.report.fft_count,
            eval_warm.report.fft_count
        );
        assert!(
            e_diff < 1.0e-10,
            "warm vs cold total_ha mismatch: {e_diff:.3e} Ha"
        );
        assert!(
            eval_warm.report.iterations < eval_cold.report.iterations,
            "warm start did not reduce DirectPcg iterations: warm={} vs cold={}",
            eval_warm.report.iterations,
            eval_cold.report.iterations
        );

        // Part A2: DirectPcg at adaptive SCF tolerance (tol = 1e-5) cuts iterations sharply.
        ws_warm.reset_warm_start();
        let _ = ws_warm
            .evaluate(comm, &n0, &ions, &params, &cfg_warm)
            .unwrap();
        let cfg_direct_adapt_warm = SolverConfig {
            tolerance: 1.0e-5,
            ..cfg_warm
        };
        let cfg_direct_adapt_cold = SolverConfig {
            warm_start: false,
            ..cfg_direct_adapt_warm
        };
        let direct_adapt_warm = ws_warm
            .evaluate(comm, &n1, &ions, &params, &cfg_direct_adapt_warm)
            .unwrap();
        let direct_adapt_cold = ws_cold
            .evaluate(comm, &n1, &ions, &params, &cfg_direct_adapt_cold)
            .unwrap();
        println!(
            "[WarmStart DirectPcg 1e-3 @ adaptive tol=1e-5] cold_iters={} warm_iters={} cold_ffts={} warm_ffts={}",
            direct_adapt_cold.report.iterations,
            direct_adapt_warm.report.iterations,
            direct_adapt_cold.report.fft_count,
            direct_adapt_warm.report.fft_count
        );
        assert!(direct_adapt_warm.report.converged && direct_adapt_warm.report.warm_started);
        assert!(
            direct_adapt_warm.report.iterations < direct_adapt_cold.report.iterations,
            "DirectPcg adaptive warm start did not reduce iterations: warm={} vs cold={}",
            direct_adapt_warm.report.iterations,
            direct_adapt_cold.report.iterations
        );

        // Part B: FisicaroPcg warm-start at tol = 1e-10 and adaptive SCF tol = 1e-5.
        // Note: on coarse dx = 0.5 Bohr grids, Fisicaro's product-rule operator differs from
        // discrete D_a^T eps D_a by O(aliasing) (which decays exponentially from 7.7e-3 Ha at
        // 32^3 to 7.2e-6 Ha at 64^3); here we verify warm/cold self-consistency (< 1e-10 Ha),
        // strict iteration reduction, and > 2x FFT reduction relative to DirectPcg.
        let cfg_fis_warm = SolverConfig {
            kind: SolverKind::FisicaroPcg,
            tolerance: 1.0e-10,
            max_iterations: 300,
            diis_history: 8,
            warm_start: true,
        };
        let cfg_fis_cold = SolverConfig {
            warm_start: false,
            ..cfg_fis_warm
        };

        ws_warm.reset_warm_start();
        let fis_base = ws_warm
            .evaluate(comm, &n0, &ions, &params, &cfg_fis_warm)
            .unwrap();
        assert!(fis_base.report.converged && !fis_base.report.warm_started);

        let fis_warm = ws_warm
            .evaluate(comm, &n1, &ions, &params, &cfg_fis_warm)
            .unwrap();
        let fis_cold = ws_cold
            .evaluate(comm, &n1, &ions, &params, &cfg_fis_cold)
            .unwrap();
        assert!(fis_warm.report.converged && fis_cold.report.converged);
        assert!(fis_warm.report.warm_started && !fis_cold.report.warm_started);

        let fis_diff_cold = (fis_warm.total_ha - fis_cold.total_ha).abs();
        let fis_diff_direct = (fis_warm.total_ha - eval_cold.total_ha).abs();
        println!(
            "[WarmStart FisicaroPcg 1e-3 @ tol=1e-10] cold_iters={} warm_iters={} cold_ffts={} warm_ffts={} |dE_cold|={fis_diff_cold:.3e} |dE_direct|={fis_diff_direct:.3e}",
            fis_cold.report.iterations,
            fis_warm.report.iterations,
            fis_cold.report.fft_count,
            fis_warm.report.fft_count
        );
        assert!(
            fis_diff_cold < 1.0e-10,
            "FisicaroPcg warm vs cold total_ha mismatch: {fis_diff_cold:.3e} Ha"
        );
        assert!(
            fis_warm.report.iterations < fis_cold.report.iterations,
            "FisicaroPcg warm start did not reduce iterations: warm={} vs cold={}",
            fis_warm.report.iterations,
            fis_cold.report.iterations
        );
        assert!(
            fis_warm.report.fft_count * 2 < eval_warm.report.fft_count,
            "FisicaroPcg expected < 50% of DirectPcg FFTs: fis={} vs direct={}",
            fis_warm.report.fft_count,
            eval_warm.report.fft_count
        );

        ws_warm.reset_warm_start();
        let _ = ws_warm
            .evaluate(comm, &n0, &ions, &params, &cfg_fis_warm)
            .unwrap();
        let cfg_adapt_warm = SolverConfig {
            tolerance: 1.0e-5,
            ..cfg_fis_warm
        };
        let cfg_adapt_cold = SolverConfig {
            warm_start: false,
            ..cfg_adapt_warm
        };
        let adapt_warm = ws_warm
            .evaluate(comm, &n1, &ions, &params, &cfg_adapt_warm)
            .unwrap();
        let adapt_cold = ws_cold
            .evaluate(comm, &n1, &ions, &params, &cfg_adapt_cold)
            .unwrap();
        println!(
            "[WarmStart FisicaroPcg 1e-3 @ adaptive tol=1e-5] cold_iters={} warm_iters={} cold_ffts={} warm_ffts={}",
            adapt_cold.report.iterations,
            adapt_warm.report.iterations,
            adapt_cold.report.fft_count,
            adapt_warm.report.fft_count
        );
        assert!(
            adapt_warm.report.converged && adapt_warm.report.warm_started,
            "FisicaroPcg adaptive warm start did not converge with warm_started=true: {:?}",
            adapt_warm.report
        );
        assert!(
            adapt_warm.report.iterations < adapt_cold.report.iterations,
            "FisicaroPcg adaptive warm start did not reduce iterations: warm={} vs cold={}",
            adapt_warm.report.iterations,
            adapt_cold.report.iterations
        );
        assert!(
            adapt_warm.report.fft_count < direct_adapt_warm.report.fft_count,
            "FisicaroPcg adaptive warm start expected fewer FFTs than DirectPcg ({}), got {}",
            direct_adapt_warm.report.fft_count,
            adapt_warm.report.fft_count
        );
    });
}

// Keep SolvationError import anchored so API error matching stays verified.
#[test]
fn invalid_solvation_parameters_are_rejected() {
    let bad_eps = SolvationParams {
        epsilon_bulk: 0.5,
        ..SolvationParams::vacuum()
    };
    assert!(matches!(
        bad_eps.validate(),
        Err(SolvationError::InvalidParameter(_))
    ));
}
