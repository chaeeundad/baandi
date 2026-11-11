//! B1 (selected-root products) and the Davidson stage timers: the switches are
//! off by default and bit-identical when off; with B1 on, decisions, iteration
//! counts and H work match the all-band path on fixed operators (only GEMM column
//! grouping differs), every returned strict solution carries a fresh certificate,
//! and results do not depend on the thread budget.
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use dft_solvers::real_davidson::{
    B1Mode, RealEigenSolution, ResearchOptions, SketchKind, b1_counters, stage_seconds,
};
use dft_solvers::{
    BandTolerances, DavidsonPreconditioner, RealDavidson, RealLinearOperator, SolverError,
};
use std::sync::Arc;

struct Operator {
    a: RealMatrix,
    pool: Arc<WorkerPool>,
}
impl RealLinearOperator for Operator {
    fn dimension(&self) -> usize {
        self.a.rows()
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        let n = self.a.rows();
        y.fill(0.);
        for (j, &xj) in x.iter().enumerate() {
            for (y, &a) in y.iter_mut().zip(&self.a.as_slice()[j * n..(j + 1) * n]) {
                *y += a * xj;
            }
        }
        Ok(())
    }
}

/// Kinetic-like diagonal with a decaying off-diagonal coupling (Ha).
fn kinetic_fixture(n: usize, threads: usize) -> (Operator, Vec<f64>) {
    let kinetic: Vec<f64> = (0..n).map(|i| 0.2 * (i / 2) as f64).collect();
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let v = if i == j {
                kinetic[i] + 0.3
            } else {
                0.025 * ((i * 31 + j * 17 + 1) as f64).sin() / ((j - i) as f64).sqrt()
            };
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
    (
        Operator {
            a,
            pool: WorkerPool::new(threads).unwrap(),
        },
        kinetic,
    )
}

/// Known spectrum in groups of `group` degenerate values, rotated by
/// deterministic Householder reflections (cached: debug builds are slow).
fn spectrum_fixture(n: usize, group: usize, threads: usize) -> (Operator, Vec<f64>) {
    static CACHE: std::sync::Mutex<Vec<(usize, usize, RealMatrix)>> =
        std::sync::Mutex::new(Vec::new());
    let mut cache = CACHE.lock().unwrap();
    if !cache.iter().any(|c| c.0 == n && c.1 == group) {
        // Columns of q are kept contiguous: q[j] is column j.
        let mut q: Vec<Vec<f64>> = (0..n)
            .map(|j| (0..n).map(|i| f64::from(u8::from(i == j))).collect())
            .collect();
        let mut state = 753_u64;
        for _ in 0..24 {
            let v: Vec<_> = (0..n)
                .map(|_| {
                    state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
                    ((state >> 11) as f64) / (1_u64 << 53) as f64 - 0.5
                })
                .collect();
            let scale = 2. / v.iter().map(|x| x * x).sum::<f64>();
            for col in &mut q {
                let d = v.iter().zip(col.iter()).map(|(a, b)| a * b).sum::<f64>() * scale;
                for (c, &vi) in col.iter_mut().zip(&v) {
                    *c -= vi * d;
                }
            }
        }
        let values: Vec<f64> = (0..n)
            .map(|i| 0.05 * ((i / group) as f64).powf(1.3) - 0.3)
            .collect();
        // a = sum_k values[k] q_k q_k^T, with q_k the k-th column of Q^T.
        let rows: Vec<Vec<f64>> = (0..n).map(|k| q.iter().map(|c| c[k]).collect()).collect();
        let mut a = RealMatrix::zeros(n, n).unwrap();
        for (k, r) in rows.iter().enumerate() {
            for j in 0..n {
                let s = values[k] * r[j];
                for (dst, &ri) in a.as_mut_slice()[j * n..(j + 1) * n].iter_mut().zip(r) {
                    *dst += s * ri;
                }
            }
        }
        for j in 0..n {
            for i in 0..j {
                let v = a[(i, j)];
                a[(j, i)] = v;
            }
        }
        cache.push((n, group, a));
    }
    let a = cache
        .iter()
        .find(|c| c.0 == n && c.1 == group)
        .unwrap()
        .2
        .clone();
    let diag: Vec<f64> = (0..n).map(|i| a[(i, i)]).collect();
    (
        Operator {
            a,
            pool: WorkerPool::new(threads).unwrap(),
        },
        diag,
    )
}

fn solver(bands: usize, tolerance: f64) -> RealDavidson {
    RealDavidson {
        bands,
        tolerance,
        max_iterations: 400,
        max_subspace: 6 * bands,
        seed: 42,
    }
}

fn run(
    s: &RealDavidson,
    op: &Operator,
    diag: &[f64],
    policy: BandTolerances,
    ritz: bool,
    options: ResearchOptions,
) -> RealEigenSolution {
    s.solve_research(
        op,
        None,
        diag,
        policy,
        ritz,
        !ritz,
        6,
        &mut None,
        DavidsonPreconditioner::Default,
        options,
    )
    .unwrap()
}

fn b1(mode: B1Mode) -> ResearchOptions {
    sketched(mode, SketchKind::Dense)
}
fn sketched(mode: B1Mode, sketch: SketchKind) -> ResearchOptions {
    ResearchOptions {
        b1: mode,
        sketch,
        ..ResearchOptions::default()
    }
}

/// Debug-build tests use the kinetic fixture only (the clustered fixtures need
/// many iterations and take minutes unoptimized); the ignored generality test
/// runs all three in release.
fn fixtures(threads: usize, full: bool) -> Vec<(&'static str, Operator, Vec<f64>)> {
    if !full {
        let (k, kd) = kinetic_fixture(180, threads);
        let (s1, d1) = spectrum_fixture(140, 1, threads);
        return vec![("kinetic", k, kd), ("spectrum-g1", s1, d1)];
    }
    let (k, kd) = kinetic_fixture(300, threads);
    let (s1, d1) = spectrum_fixture(240, 1, threads);
    let (s3, d3) = spectrum_fixture(240, 3, threads);
    vec![
        ("kinetic", k, kd),
        ("spectrum-g1", s1, d1),
        ("spectrum-g3", s3, d3),
    ]
}

/// Exactly degenerate clusters (spectrum-g3) make Ritz vectors inside a cluster
/// non-unique, so a rounding-class change may move the trajectory by an
/// iteration (as MKL-EVD did on 2 of 30 dev08 cells); there the counts may
/// differ by at most 10 percent and eigenvalues agree to 1e-10 Ha.
fn same_decisions(name: &str, base: &RealEigenSolution, b1: &RealEigenSolution, tol: f64) {
    if name == "spectrum-g3" {
        let near = |a: usize, b: usize| a.abs_diff(b) * 10 <= a.max(1);
        assert!(near(base.iterations, b1.iterations), "{name}: iterations");
        assert!(
            near(base.operator_applications, b1.operator_applications),
            "{name}: H columns"
        );
        for (a, b) in base.eigenvalues.iter().zip(&b1.eigenvalues) {
            assert!((a - b).abs() <= 1e-10, "{name}: eigenvalue {a} vs {b}");
        }
        return;
    }
    assert_eq!(base.iterations, b1.iterations, "{name}: iterations");
    assert_eq!(
        base.operator_applications, b1.operator_applications,
        "{name}: H columns"
    );
    assert_eq!(base.restarts, b1.restarts, "{name}: restarts");
    assert_eq!(base.converged, b1.converged, "{name}: flags");
    for (a, b) in base.eigenvalues.iter().zip(&b1.eigenvalues) {
        assert!((a - b).abs() <= tol, "{name}: eigenvalue {a} vs {b}");
    }
}

fn default_path(full: bool) {
    for (name, op, diag) in fixtures(1, full) {
        let s = solver(16, 1e-9);
        let policy = BandTolerances {
            strict_bands: 12,
            buffer_tolerance: 1e-6,
        };
        let a = s
            .solve_with_band_tolerances(&op, None, &diag, policy)
            .unwrap();
        let b = run(&s, &op, &diag, policy, false, ResearchOptions::from_env());
        assert_eq!(a, b, "{name}");
        let a = s.solve_ritz_changes(&op, None, &diag, policy).unwrap();
        let b = run(&s, &op, &diag, policy, true, ResearchOptions::from_env());
        assert_eq!(a, b, "{name}");
    }
}

fn standard_profile(full: bool) {
    let before = b1_counters();
    for (name, op, diag) in fixtures(1, full) {
        let s = solver(16, 1e-10);
        let policy = BandTolerances {
            strict_bands: 12,
            buffer_tolerance: 1e-8,
        };
        let base = run(&s, &op, &diag, policy, true, ResearchOptions::default());
        // `strict` does not act on Ritz-change solves.
        let strict_only = run(&s, &op, &diag, policy, true, b1(B1Mode::Strict));
        assert_eq!(base, strict_only, "{name}");
        let modes: &[B1Mode] = if full {
            &[B1Mode::Standard, B1Mode::On]
        } else {
            &[B1Mode::On]
        };
        for &mode in modes {
            let fast = run(&s, &op, &diag, policy, true, b1(mode));
            same_decisions(name, &base, &fast, 1e-12);
            // Returned residuals are formed on the all-band path; they differ
            // only by rounding (operator norm up to about 60 Ha here).
            // Inside exactly degenerate clusters (spectrum-g3) the Ritz vectors
            // are not unique and Ritz-change convergence leaves them rotating,
            // so rounding changes individual residual norms there; only the
            // nondegenerate fixtures compare them.
            if name != "spectrum-g3" {
                for (a, b) in base.residuals.iter().zip(&fast.residuals) {
                    assert!((a - b).abs() <= 1e-12, "{name}: {a} vs {b}");
                }
            }
            assert!(!fast.residuals_certified);
        }
    }
    let after = b1_counters();
    assert!(after.narrow_iterations > before.narrow_iterations);
    assert!(after.avoided_columns > before.avoided_columns);
}

fn strict_profile(full: bool) {
    let before = b1_counters();
    for (name, op, diag) in fixtures(1, full) {
        let s = solver(16, 1e-9);
        let policy = BandTolerances {
            strict_bands: 12,
            buffer_tolerance: 1e-6,
        };
        let base = run(&s, &op, &diag, policy, false, ResearchOptions::default());
        let standard_only = run(&s, &op, &diag, policy, false, b1(B1Mode::Standard));
        assert_eq!(base, standard_only, "{name}");
        let modes: &[B1Mode] = if full {
            &[B1Mode::Strict, B1Mode::On]
        } else {
            &[B1Mode::Strict]
        };
        for (&mode, sketch) in modes
            .iter()
            .flat_map(|m| [(m, SketchKind::Dense), (m, SketchKind::Count)])
        {
            let fast = run(&s, &op, &diag, policy, false, sketched(mode, sketch));
            same_decisions(name, &base, &fast, 1e-12);
            assert!(fast.residuals_certified, "{name}");
            assert_eq!(base.certification_attempts, fast.certification_attempts);
            for (j, &r) in fast.residuals.iter().enumerate() {
                let limit = if j < policy.strict_bands {
                    s.tolerance
                } else {
                    policy.buffer_tolerance
                };
                assert!(r <= limit, "{name}: root {j} residual {r}");
            }
        }
    }
    let after = b1_counters();
    assert!(after.narrow_iterations > before.narrow_iterations);
    assert!(after.sketch_estimates > before.sketch_estimates);
    // A root flagged by its estimate never failed the explicit test.
    assert_eq!(after.sketch_disagreements, before.sketch_disagreements);
}

fn thread_identity(full: bool) {
    let policy = BandTolerances {
        strict_bands: 12,
        buffer_tolerance: 1e-6,
    };
    let s = solver(16, 1e-9);
    let one = fixtures(1, full);
    let four = fixtures(4, full);
    for ((name, op1, d1), (_, op4, d4)) in one.iter().zip(&four) {
        for (ritz, sketch) in [
            (false, SketchKind::Dense),
            (false, SketchKind::Count),
            (true, SketchKind::Dense),
        ] {
            let a = run(&s, op1, d1, policy, ritz, sketched(B1Mode::On, sketch));
            let b = run(&s, op4, d4, policy, ritz, sketched(B1Mode::On, sketch));
            assert_eq!(a, b, "{name} ritz={ritz} {sketch:?}");
        }
    }
}

#[test]
fn default_options_are_the_default_path() {
    // The environment switches are unset in the test process; B6 is on by default
    // since the adoption rule (dev08 MICRO-ARM1), everything else is off.
    assert_eq!(
        ResearchOptions::from_env(),
        ResearchOptions {
            b6: true,
            b1: B1Mode::On,
            sketch: SketchKind::Count,
            b5: true,
            ..ResearchOptions::default()
        }
    );
    default_path(false);
}
#[test]
fn standard_profile_keeps_every_decision() {
    standard_profile(false);
}
#[test]
fn strict_profile_keeps_decisions_and_certifies_every_return() {
    strict_profile(false);
}
#[test]
fn b1_results_do_not_depend_on_threads() {
    thread_identity(false);
}
/// Generality set (kinetic, nondegenerate and degenerate clustered spectra):
/// `cargo test --release -p dft-solvers --test b1_selected_roots -- --ignored`
#[test]
#[ignore = "release only: clustered fixtures take minutes unoptimized"]
fn b1_generality_set() {
    default_path(true);
    standard_profile(true);
    strict_profile(true);
    thread_identity(true);
}

#[test]
fn stage_timers_count_every_iteration_without_changing_results() {
    let (op, diag) = kinetic_fixture(140, 1);
    let s = solver(12, 1e-9);
    let policy = BandTolerances {
        strict_bands: 10,
        buffer_tolerance: 1e-6,
    };
    let before = stage_seconds();
    let base = run(&s, &op, &diag, policy, false, ResearchOptions::default());
    let timed = run(
        &s,
        &op,
        &diag,
        policy,
        false,
        ResearchOptions {
            stage_timers: true,
            ..ResearchOptions::default()
        },
    );
    assert_eq!(base, timed);
    let after = stage_seconds();
    assert!(after.timed_solves > before.timed_solves);
    let counted: u64 =
        after.width_iterations.iter().sum::<u64>() - before.width_iterations.iter().sum::<u64>();
    // Other tests may time solves concurrently; this solve alone adds its own.
    assert!(counted >= timed.iterations as u64);
    assert!(after.h_apply > before.h_apply);
    assert!(after.eigen > before.eigen);
}
