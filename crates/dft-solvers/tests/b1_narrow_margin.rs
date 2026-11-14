//! Gap-analysis rank 2 (docs/performance/davidson-algebra-gap-analysis.md):
//! `BAANDI_B1_NARROW` (narrow path for any proper subset of the bands),
//! `BAANDI_B1_MARGIN` and `BAANDI_B1_ROWS` (strict sketch margin and rows).
//! Off / default values reproduce the previous path bit for bit; with the
//! switches on the answer stays inside the accuracy class, every strict return
//! carries a fresh all-band certificate, a root skipped by the sketch is caught
//! by the explicit all-band residual before any return, and results do not
//! depend on the thread budget.
//!
//! Own test binary: the B1 counters are process-wide and the aggressive-margin
//! test produces sketch disagreements on purpose, which `b1_selected_roots`
//! asserts to be zero. Tests here take a lock so counter deltas are their own.
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use dft_solvers::real_davidson::{
    B1Mode, RealEigenSolution, ResearchOptions, SketchKind, b1_counters,
};
use dft_solvers::{
    BandTolerances, DavidsonPreconditioner, RealDavidson, RealLinearOperator, SolverError,
};
use std::sync::{Arc, Mutex};

static SERIAL: Mutex<()> = Mutex::new(());

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

/// Kinetic-like diagonal (Ha) with a decaying off-diagonal coupling; `coupling`
/// varies the coupling so two shapes of spectrum are covered.
fn fixture(n: usize, coupling: f64, threads: usize) -> (Operator, Vec<f64>) {
    let kinetic: Vec<f64> = (0..n).map(|i| 0.2 * (i / 2) as f64).collect();
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let v = if i == j {
                kinetic[i] + 0.3
            } else {
                coupling * ((i * 31 + j * 17 + 1) as f64).sin() / ((j - i) as f64).sqrt()
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

fn fixtures(threads: usize) -> Vec<(&'static str, Operator, Vec<f64>)> {
    let (a, da) = fixture(180, 0.025, threads);
    let (b, db) = fixture(150, 0.08, threads);
    vec![("weak-coupling", a, da), ("strong-coupling", b, db)]
}

fn solver() -> RealDavidson {
    RealDavidson {
        bands: 16,
        tolerance: 1e-9,
        max_iterations: 400,
        max_subspace: 96,
        seed: 42,
    }
}
const POLICY: BandTolerances = BandTolerances {
    strict_bands: 12,
    buffer_tolerance: 1e-6,
};

fn run(op: &Operator, diag: &[f64], ritz: bool, options: ResearchOptions) -> RealEigenSolution {
    solver()
        .solve_research(
            op,
            None,
            diag,
            POLICY,
            ritz,
            !ritz,
            6,
            &mut None,
            DavidsonPreconditioner::Default,
            options,
        )
        .unwrap()
}

/// B1 as in production (count sketch) with the rank-2 knobs.
fn b1(narrow: bool, margin: f64, rows: Option<usize>) -> ResearchOptions {
    ResearchOptions {
        b1: B1Mode::On,
        sketch: SketchKind::Count,
        b1_narrow: narrow,
        b1_margin: margin,
        b1_rows: rows,
        ..ResearchOptions::default()
    }
}

fn limit(j: usize) -> f64 {
    if j < POLICY.strict_bands {
        solver().tolerance
    } else {
        POLICY.buffer_tolerance
    }
}

/// ||A x_j - lambda_j x_j|| recomputed from the returned vectors (Ha).
fn true_residuals(op: &Operator, s: &RealEigenSolution) -> Vec<f64> {
    let n = op.dimension();
    let x = s.eigenvectors.as_slice();
    let mut y = vec![0.; n];
    (0..s.eigenvalues.len())
        .map(|j| {
            op.apply(&x[j * n..(j + 1) * n], &mut y).unwrap();
            y.iter()
                .zip(&x[j * n..(j + 1) * n])
                .map(|(&h, &v)| (h - s.eigenvalues[j] * v).powi(2))
                .sum::<f64>()
                .sqrt()
        })
        .collect()
}

/// Strict return: certified, every residual (reported and recomputed) inside
/// its threshold; recomputation differs from the reported value by rounding.
fn assert_certified(name: &str, s: &RealEigenSolution, op: &Operator) {
    assert!(s.residuals_certified, "{name}: not certified");
    for (j, (&r, t)) in s.residuals.iter().zip(true_residuals(op, s)).enumerate() {
        assert!(r <= limit(j), "{name}: root {j} reported {r}");
        assert!(
            t <= limit(j) * (1. + 1e-6) + 1e-14,
            "{name}: root {j} true {t}"
        );
    }
}

#[test]
fn defaults_reproduce_the_previous_path() {
    let _g = SERIAL.lock().unwrap();
    // Environment unset in tests: the new switches are off / 0.5 / automatic.
    let env = ResearchOptions::from_env();
    assert!(!env.b1_narrow);
    assert_eq!(env.b1_margin, 0.5);
    assert_eq!(env.b1_rows, None);
    // The automatic row rule gives the former constant 128 at margin 0.5 and
    // keeps the 8-sigma miss bound at wider margins.
    assert_eq!(env.b1_sketch_rows(), 128);
    assert_eq!(b1(false, 0.7, None).b1_sketch_rows(), 356);
    assert_eq!(b1(false, 0.75, None).b1_sketch_rows(), 512);
    for (name, op, diag) in fixtures(1) {
        for ritz in [false, true] {
            let base = run(&op, &diag, ritz, env);
            // Explicit former values: bit-identical.
            let explicit = run(&op, &diag, ritz, b1(false, 0.5, Some(128)));
            assert_eq!(base, explicit, "{name} ritz={ritz}");
        }
    }
}

#[test]
fn narrow_path_agrees_with_full_path_within_the_class() {
    let _g = SERIAL.lock().unwrap();
    let before = b1_counters();
    let mut narrow_extra = 0_u64;
    for (name, op, diag) in fixtures(1) {
        for ritz in [false, true] {
            let off = b1_counters();
            let base = run(&op, &diag, ritz, b1(false, 0.5, None));
            let mid = b1_counters();
            let fast = run(&op, &diag, ritz, b1(true, 0.5, None));
            let after = b1_counters();
            narrow_extra += (after.narrow_iterations - mid.narrow_iterations)
                .saturating_sub(mid.narrow_iterations - off.narrow_iterations);
            // Decisions are the same (exact Ritz-change flags; sketch flags
            // with the unchanged margin); only GEMM column grouping differs.
            assert_eq!(base.iterations, fast.iterations, "{name} ritz={ritz}");
            assert_eq!(
                base.operator_applications, fast.operator_applications,
                "{name} ritz={ritz}"
            );
            assert_eq!(base.converged, fast.converged, "{name} ritz={ritz}");
            for (a, b) in base.eigenvalues.iter().zip(&fast.eigenvalues) {
                assert!((a - b).abs() <= 1e-12, "{name} ritz={ritz}: {a} vs {b}");
            }
            if ritz {
                assert!(!fast.residuals_certified);
            } else {
                assert_certified(name, &fast, &op);
                assert_eq!(base.certification_attempts, fast.certification_attempts);
            }
        }
    }
    // The switch acted: more narrow iterations than the half-band gate.
    assert!(narrow_extra > 0, "B1_NARROW never widened the narrow path");
    assert_eq!(
        b1_counters().sketch_disagreements,
        before.sketch_disagreements
    );
}

/// A margin near 1 with a tiny sketch makes misses likely: roots above their
/// threshold get estimates below margin x threshold and skip their explicit
/// residual. They must be caught by the explicit all-band residual of the next
/// full-path iteration (counted as sketch disagreements) and never returned.
#[test]
fn skipped_roots_are_rechecked_before_return() {
    let _g = SERIAL.lock().unwrap();
    let before = b1_counters();
    for (name, op, diag) in fixtures(1) {
        let base = run(&op, &diag, false, ResearchOptions::default());
        for (narrow, margin, rows) in [(true, 0.95, Some(8)), (false, 0.9, Some(8))] {
            let risky = run(&op, &diag, false, b1(narrow, margin, rows));
            assert_certified(name, &risky, &op);
            for (a, b) in base.eigenvalues.iter().zip(&risky.eigenvalues) {
                assert!((a - b).abs() <= 1e-10, "{name}: {a} vs {b}");
            }
        }
    }
    let after = b1_counters();
    assert!(
        after.sketch_disagreements > before.sketch_disagreements,
        "the aggressive sketch never skipped a root above threshold"
    );
}

#[test]
fn narrow_and_margin_do_not_depend_on_threads() {
    let _g = SERIAL.lock().unwrap();
    let one = fixtures(1);
    let four = fixtures(4);
    for ((name, op1, d1), (_, op4, d4)) in one.iter().zip(&four) {
        for (ritz, options) in [
            (false, b1(true, 0.5, None)),
            (true, b1(true, 0.5, None)),
            (false, b1(true, 0.7, None)),
            (false, b1(true, 0.95, Some(8))),
        ] {
            let a = run(op1, d1, ritz, options);
            let b = run(op4, d4, ritz, options);
            assert_eq!(a, b, "{name} ritz={ritz} {options:?}");
        }
    }
}
