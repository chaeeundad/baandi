//! Track B Gate 2 release benchmark comparing `DirectPcg`, `FisicaroPcg`, and
//! `PolarizationDiis` across smooth MMS and sharp VASPsol water cavities,
//! cold vs warm-start (`1e-3` and `1e-5` perturbations), and `dft_mpi::simulate`
//! multi-rank parity & timing (`P = 1, 2, 4`).
//!
//! Run:
//! `cargo run --release -p dft-physics --example solvation_bench`
//! `cargo run --release -p dft-physics --example solvation_bench -- --96`

use std::f64::consts::PI;
use std::time::Instant;

use dft_geometry::Cell;
use dft_physics::solvation::{
    IonicPseudoCharge, SolvationError, SolvationParams, SolvationWorkspace, SolveReport,
    SolverConfig, SolverKind,
};
use num_complex::Complex64;

fn cubic_cell(length_bohr: f64) -> Cell {
    Cell::from_user_rows_bohr([
        [length_bohr, 0.0, 0.0],
        [0.0, length_bohr, 0.0],
        [0.0, 0.0, length_bohr],
    ])
    .unwrap()
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

fn build_smooth_mms(
    dims: [usize; 3],
    owned_indices: &[usize],
    eps_bulk: f64,
    kappa: f64,
) -> (Vec<f64>, Vec<f64>, Vec<f64>) {
    let n = owned_indices.len();
    let mut eps = Vec::with_capacity(n);
    let mut kappa2 = Vec::with_capacity(n);
    let mut phi = Vec::with_capacity(n);
    let dc = if kappa > 0.0 { 0.35 } else { 0.0 };

    for &idx in owned_indices {
        let [u, v, w] = grid_frac(idx, dims);
        let t0 = 2.0 * PI * u;
        let t1 = 2.0 * PI * v;
        let t2 = 2.0 * PI * w;
        let f = 0.5 + 0.25 * t0.cos() + 0.15 * t1.sin() + 0.10 * t2.cos();
        let e = 1.0 + (eps_bulk - 1.0) * f;
        eps.push(e);
        kappa2.push(eps_bulk * kappa * kappa * f);
        phi.push(
            dc + 0.7 * t0.cos() - 0.5 * t1.sin() + 0.4 * (t0 + t2).cos() + 0.3 * (t1 - t2).sin(),
        );
    }
    (eps, kappa2, phi)
}

fn solver_label(kind: SolverKind) -> &'static str {
    match kind {
        SolverKind::DirectPcg => "DirectPcg",
        SolverKind::FisicaroPcg => "FisicaroPcg",
        SolverKind::PolarizationDiis => "PolarizationDiis",
        SolverKind::HybridFisicaroDirect => "HybridFisicaroDirect",
    }
}

fn format_report_row(
    grid: usize,
    case: &str,
    kappa: f64,
    kind: SolverKind,
    mode: &str,
    elapsed_ms: f64,
    outcome: &Result<SolveReport, SolvationError>,
) -> String {
    match outcome {
        Ok(rep) => format!(
            "| {grid}^3 | {case} | {kappa:.2} | {} | {mode} | OK | {} | {} | {elapsed_ms:.2} | {:.2e} | {:.2e} | {} |",
            solver_label(kind),
            rep.iterations,
            rep.fft_count,
            rep.residual_norm,
            rep.relative_residual,
            rep.warm_started
        ),
        Err(SolvationError::NotConverged {
            iterations,
            residual_norm,
        }) => format!(
            "| {grid}^3 | {case} | {kappa:.2} | {} | {mode} | NOT_CONV | {iterations} | - | {elapsed_ms:.2} | {residual_norm:.2e} | - | - |",
            solver_label(kind)
        ),
        Err(SolvationError::Breakdown(reason)) => format!(
            "| {grid}^3 | {case} | {kappa:.2} | {} | {mode} | BREAKDOWN({reason}) | - | - | {elapsed_ms:.2} | - | - | - |",
            solver_label(kind)
        ),
        Err(other) => format!(
            "| {grid}^3 | {case} | {kappa:.2} | {} | {mode} | ERR({other}) | - | - | {elapsed_ms:.2} | - | - | - |",
            solver_label(kind)
        ),
    }
}

fn benchmark_mms_grid(n_axis: usize) {
    let cell = cubic_cell(16.0);
    let dims = [n_axis, n_axis, n_axis];
    let solvers = [
        SolverKind::DirectPcg,
        SolverKind::FisicaroPcg,
        SolverKind::PolarizationDiis,
        SolverKind::HybridFisicaroDirect,
    ];

    for kappa in [0.0_f64, 0.15_f64] {
        dft_mpi::simulate(1, |comm| {
            let mut ws = SolvationWorkspace::new(&cell, dims, comm).unwrap();
            let owned = ws.owned_real_grid_indices();
            let (eps, kappa2, phi_exact) = build_smooth_mms(dims, &owned, 78.4, kappa);
            let n_g = ws.local_positions().len();
            let mut phi_exact_g = vec![Complex64::new(0.0, 0.0); n_g];
            let mut rhs_g = vec![Complex64::new(0.0, 0.0); n_g];
            ws.forward_real(comm, &phi_exact, &mut phi_exact_g).unwrap();
            ws.apply_operator(comm, &eps, &kappa2, &phi_exact_g, &mut rhs_g)
                .unwrap();

            for kind in solvers {
                ws.reset_warm_start();
                let cfg = SolverConfig {
                    kind,
                    tolerance: 1.0e-10,
                    max_iterations: 300,
                    diis_history: 8,
                    warm_start: false,
                };
                let mut phi_out_g = vec![Complex64::new(0.0, 0.0); n_g];
                let t0 = Instant::now();
                let res = ws.solve_pde(comm, &eps, &kappa2, &rhs_g, &cfg, &mut phi_out_g);
                let ms = t0.elapsed().as_secs_f64() * 1.0e3;
                println!(
                    "{}",
                    format_report_row(n_axis, "MMS(eps=78.4)", kappa, kind, "cold", ms, &res)
                );
            }
        });
    }
}

fn benchmark_vaspsol_grid(n_axis: usize) {
    let cell = cubic_cell(16.0);
    let dims = [n_axis, n_axis, n_axis];
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
    let global_n0 = build_global_molecular_density(&cell, dims, 4.0);
    let solvers = [
        SolverKind::DirectPcg,
        SolverKind::FisicaroPcg,
        SolverKind::PolarizationDiis,
        SolverKind::HybridFisicaroDirect,
    ];

    for kappa in [0.0_f64, 0.15_f64] {
        let params = SolvationParams {
            epsilon_bulk: 78.4,
            kappa_inv_bohr: kappa,
            tau_ha_bohr2: 2.0e-5,
            ..SolvationParams::vaspsol_water()
        };

        dft_mpi::simulate(1, |comm| {
            let mut ws = SolvationWorkspace::new(&cell, dims, comm).unwrap();
            let owned = ws.owned_real_grid_indices();
            let n0: Vec<f64> = owned.iter().map(|&idx| global_n0[idx]).collect();
            let n_1e3 = perturb_density(&n0, &owned, dims, 1.0e-3);
            let n_1e5 = perturb_density(&n0, &owned, dims, 1.0e-5);

            let mut direct_cold_ha: Option<f64> = None;

            for kind in solvers {
                ws.reset_warm_start();
                let cfg_cold = SolverConfig {
                    kind,
                    tolerance: 1.0e-10,
                    max_iterations: 300,
                    diis_history: 8,
                    warm_start: false,
                };
                let t0 = Instant::now();
                let cold_eval = ws.evaluate(comm, &n0, &ions, &params, &cfg_cold);
                let cold_ms = t0.elapsed().as_secs_f64() * 1.0e3;
                if let Ok(ref ev) = cold_eval {
                    if kind == SolverKind::DirectPcg {
                        direct_cold_ha = Some(ev.total_ha);
                    } else if let Some(ref_ha) = direct_cold_ha {
                        let d_ha = (ev.total_ha - ref_ha).abs();
                        eprintln!(
                            "[OperatorParity {n_axis}^3 kappa={kappa:.2} {}] E={:.10} Ha, |E - E_direct|={d_ha:.3e} Ha",
                            solver_label(kind),
                            ev.total_ha
                        );
                    }
                }
                let cold_outcome = cold_eval.map(|e| e.report);
                println!(
                    "{}",
                    format_report_row(
                        n_axis,
                        "VASPsol-erfc(water)",
                        kappa,
                        kind,
                        "cold",
                        cold_ms,
                        &cold_outcome
                    )
                );

                if cold_outcome.is_ok() {
                    let cfg_warm = SolverConfig {
                        kind,
                        tolerance: 1.0e-10,
                        max_iterations: 300,
                        diis_history: 8,
                        warm_start: true,
                    };
                    for (mode, n_pert, tol) in [
                        ("warm(1e-3)", &n_1e3, 1.0e-10),
                        ("warm(1e-5)", &n_1e5, 1.0e-10),
                        ("warm(1e-3)@1e-5", &n_1e3, 1.0e-5),
                    ] {
                        if tol > 1.0e-9
                            && !matches!(kind, SolverKind::DirectPcg | SolverKind::FisicaroPcg)
                        {
                            continue;
                        }
                        ws.reset_warm_start();
                        if ws.evaluate(comm, &n0, &ions, &params, &cfg_warm).is_ok() {
                            let cfg_run = SolverConfig {
                                tolerance: tol,
                                ..cfg_warm
                            };
                            let tw = Instant::now();
                            let w_res = ws
                                .evaluate(comm, n_pert, &ions, &params, &cfg_run)
                                .map(|e| e.report);
                            let w_ms = tw.elapsed().as_secs_f64() * 1.0e3;
                            println!(
                                "{}",
                                format_report_row(
                                    n_axis,
                                    "VASPsol-erfc(water)",
                                    kappa,
                                    kind,
                                    mode,
                                    w_ms,
                                    &w_res
                                )
                            );
                        }
                    }
                }
            }
        });
    }
}

fn benchmark_mpi_scaling(grids: &[usize]) {
    println!();
    println!(
        "## Simulated MPI Scaling & Parity (`DirectPcg` & `FisicaroPcg`, VASPsol-erfc water, `kappa = 0.15`)"
    );
    println!();
    println!(
        "| Grid | Solver | Ranks | Mode | Mirror | Iters | FFTs | Wall (ms) | Total (Ha) | `|E(P) - E(1)|` |"
    );
    println!("| --- | --- | ---: | --- | --- | ---: | ---: | ---: | ---: | ---: |");

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
    let params = SolvationParams {
        epsilon_bulk: 78.4,
        kappa_inv_bohr: 0.15,
        tau_ha_bohr2: 2.0e-5,
        ..SolvationParams::vaspsol_water()
    };

    for &n_axis in grids {
        let cell = cubic_cell(16.0);
        let dims = [n_axis, n_axis, n_axis];
        let global_n0 = build_global_molecular_density(&cell, dims, 4.0);

        for kind in [SolverKind::DirectPcg, SolverKind::FisicaroPcg] {
            let cfg = SolverConfig {
                kind,
                tolerance: 1.0e-10,
                max_iterations: 300,
                diis_history: 8,
                warm_start: true,
            };
            let mut ref_cold_ha: Option<f64> = None;
            let mut ref_warm_ha: Option<f64> = None;

            for ranks in [1usize, 2, 4] {
                let runs = dft_mpi::simulate(ranks, |comm| {
                    let mut ws = SolvationWorkspace::new(&cell, dims, comm).unwrap();
                    let owned = ws.owned_real_grid_indices();
                    let n0: Vec<f64> = owned.iter().map(|&idx| global_n0[idx]).collect();
                    let n1 = perturb_density(&n0, &owned, dims, 1.0e-3);

                    let tc = Instant::now();
                    let cold = ws.evaluate(comm, &n0, &ions, &params, &cfg).unwrap();
                    let cold_ms = tc.elapsed().as_secs_f64() * 1.0e3;

                    let tw = Instant::now();
                    let warm = ws.evaluate(comm, &n1, &ions, &params, &cfg).unwrap();
                    let warm_ms = tw.elapsed().as_secs_f64() * 1.0e3;

                    (
                        cold.total_ha,
                        cold.report,
                        cold_ms,
                        warm.total_ha,
                        warm.report,
                        warm_ms,
                    )
                });

                let (cold_ha, cold_rep, _, warm_ha, warm_rep, _) = runs[0];
                for r in &runs[1..] {
                    assert!((r.0 - cold_ha).abs() < 1.0e-12);
                    assert!((r.3 - warm_ha).abs() < 1.0e-12);
                }
                let cold_ms = runs.iter().map(|r| r.2).fold(0.0_f64, f64::max);
                let warm_ms = runs.iter().map(|r| r.5).fold(0.0_f64, f64::max);

                let e1_cold = *ref_cold_ha.get_or_insert(cold_ha);
                let e1_warm = *ref_warm_ha.get_or_insert(warm_ha);

                println!(
                    "| {n_axis}^3 | {} | {ranks} | cold | {} | {} | {} | {cold_ms:.2} | {cold_ha:.10} | {:.2e} |",
                    solver_label(kind),
                    cold_rep.mirror_paired,
                    cold_rep.iterations,
                    cold_rep.fft_count,
                    (cold_ha - e1_cold).abs()
                );
                println!(
                    "| {n_axis}^3 | {} | {ranks} | warm(1e-3) | {} | {} | {} | {warm_ms:.2} | {warm_ha:.10} | {:.2e} |",
                    solver_label(kind),
                    warm_rep.mirror_paired,
                    warm_rep.iterations,
                    warm_rep.fft_count,
                    (warm_ha - e1_warm).abs()
                );
            }
        }
    }
}

fn main() {
    let include_96 = std::env::args().any(|arg| arg == "--96");
    let mut grids = vec![32usize, 48, 64];
    if include_96 {
        grids.push(96);
    }

    println!("# Track B Gate 2 Solvation Solver Benchmark");
    println!();
    println!("## Solver Arms: Smooth MMS vs Sharp VASPsol Water (`eps_b = 78.4`, `tol = 1e-10`)");
    println!();
    println!(
        "| Grid | Case | `kappa` | Solver | Mode | Status | Iters | FFTs | Wall (ms) | `||r||_2` | Rel Res | Warm |"
    );
    println!("| --- | --- | ---: | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | --- |");

    for &grid in &grids {
        benchmark_mms_grid(grid);
        benchmark_vaspsol_grid(grid);
    }

    benchmark_mpi_scaling(&[32, 48]);
}
