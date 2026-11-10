//! Algebra-stream switches of the serial real Davidson
//! (docs/performance/requests/algebra-stream.md), on the fixed operators of
//! `tests/real_davidson.rs` (kinetic-coupled, Householder spectra with known
//! eigenvalues, degenerate triples) and the clustered spectra of
//! `tests/distributed_davidson.rs`:
//! * RITZ-FUSED (`ritz_fused`) returns the default solution bit for bit;
//! * NARROW-PROJ, NARROW-REDUCED and B5=par (rounding class) keep every root
//!   within 1e-11 Ha of the default path, certified strict returns, an
//!   orthonormal result (defect below 1e-12) and comparable work;
//! * every switch gives the same bits at pool sizes 1 and 4 (the large-basis
//!   dispatch case is the ignored release test at the end).
//!
//! Eigenvalues in Ha; tolerances dimensionless where stated.
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use dft_solvers::real_davidson::{
    B1Mode, RealEigenSolution, ResearchOptions, SketchKind, algebra_counters,
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

fn lcg(state: &mut u64) -> f64 {
    *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
    ((*state >> 11) as f64) / (1_u64 << 53) as f64 - 0.5
}
/// Kinetic-like diagonal (degenerate pairs, Ha) with decaying coupling.
fn kinetic(n: usize, coupling: f64) -> (RealMatrix, Vec<f64>) {
    let diag: Vec<f64> = (0..n).map(|i| 0.2 * (i / 2) as f64).collect();
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let v = if i == j {
                diag[i] + 0.3
            } else {
                coupling * ((i * 31 + j * 17 + 1) as f64).sin() / ((j - i) as f64).sqrt()
            };
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
    (a, diag)
}
/// Q diag(values) Q^T with Q from `reflections` Householder reflections; the
/// diagonal of the matrix is the preconditioner. Returns the sorted spectrum.
fn spectrum(n: usize, reflections: usize, values: &[f64]) -> (RealMatrix, Vec<f64>, Vec<f64>) {
    let mut q = RealMatrix::zeros(n, n).unwrap();
    for i in 0..n {
        q[(i, i)] = 1.;
    }
    let mut state = 753_u64;
    for _ in 0..reflections {
        let v: Vec<_> = (0..n).map(|_| lcg(&mut state)).collect();
        let scale = 2. / v.iter().map(|x| x * x).sum::<f64>();
        for j in 0..n {
            let d = (0..n).map(|i| v[i] * q[(i, j)]).sum::<f64>() * scale;
            for i in 0..n {
                q[(i, j)] -= v[i] * d;
            }
        }
    }
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let v: f64 = (0..n).map(|k| q[(i, k)] * values[k] * q[(j, k)]).sum();
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
    let diag = (0..n).map(|i| a[(i, i)]).collect();
    let mut sorted = values.to_vec();
    sorted.sort_by(f64::total_cmp);
    (a, diag, sorted)
}

struct Case {
    name: &'static str,
    a: RealMatrix,
    diag: Vec<f64>,
    oracle: Option<Vec<f64>>,
    bands: usize,
    strict: usize,
}
/// Quick set for the debug test run; the full set is the ignored release test.
fn quick_cases() -> Vec<Case> {
    let (a, diag) = kinetic(120, 0.08);
    let values: Vec<f64> = (0..96)
        .map(|i| 0.05 * ((i / 4) as f64).powf(1.3) - 0.3)
        .collect();
    let (b, bdiag, sorted) = spectrum(96, 12, &values);
    vec![
        Case {
            name: "kinetic-120",
            a,
            diag,
            oracle: None,
            bands: 10,
            strict: 8,
        },
        Case {
            name: "clustered-96",
            a: b,
            diag: bdiag,
            oracle: Some(sorted),
            bands: 12,
            strict: 10,
        },
    ]
}
fn cases() -> Vec<Case> {
    let mut out = Vec::new();
    for (name, coupling) in [("kinetic-0.025", 0.025), ("kinetic-0.08", 0.08)] {
        let (a, diag) = kinetic(220, coupling);
        out.push(Case {
            name,
            a,
            diag,
            oracle: None,
            bands: 16,
            strict: 12,
        });
    }
    // Known spectra: spread, degenerate triples, clustered groups of four
    // (the clustered strict-prefix fixture of the distributed tests).
    let n = 150;
    let spread: Vec<f64> = (0..n).map(|i| i as f64 * 0.13 - 0.7).collect();
    let triples: Vec<f64> = (0..n).map(|i| (i / 3) as f64 * 0.2).collect();
    let clustered: Vec<f64> = (0..240)
        .map(|i| 0.05 * ((i / 4) as f64).powf(1.3) - 0.3)
        .collect();
    for (name, values, reflections, bands, strict) in [
        ("householder-spread", spread, n, 12, 10),
        ("householder-triples", triples, n, 13, 9),
        ("clustered-groups", clustered, 24, 24, 20),
    ] {
        let (a, diag, sorted) = spectrum(values.len(), reflections, &values);
        out.push(Case {
            name,
            a,
            diag,
            oracle: Some(sorted),
            bands,
            strict,
        });
    }
    out
}

fn solve(case: &Case, threads: usize, ritz: bool, options: ResearchOptions) -> RealEigenSolution {
    let op = Operator {
        a: case.a.clone(),
        pool: WorkerPool::new(threads).unwrap(),
    };
    solve_op(&op, &case.diag, case.bands, case.strict, ritz, options)
}
fn solve_op(
    op: &dyn RealLinearOperator,
    diag: &[f64],
    bands: usize,
    strict: usize,
    ritz: bool,
    options: ResearchOptions,
) -> RealEigenSolution {
    let s = RealDavidson {
        bands,
        tolerance: if ritz { 1e-11 } else { 1e-10 },
        max_iterations: 600,
        max_subspace: 6 * bands,
        seed: 42,
    };
    let policy = BandTolerances {
        strict_bands: strict,
        buffer_tolerance: if ritz { 1e-9 } else { 1e-7 },
    };
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

/// Sparse operator for the large-basis dispatch test: diagonal plus nearest
/// and long-range couplings (Ha), O(n) per column.
struct Sparse {
    pool: Arc<WorkerPool>,
    n: usize,
}
impl RealLinearOperator for Sparse {
    fn dimension(&self) -> usize {
        self.n
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        let n = self.n;
        for i in 0..n {
            y[i] = (i as f64 * 0.004 + 0.5) * x[i];
            if i > 0 {
                y[i] += 0.02 * x[i - 1];
            }
            if i + 1 < n {
                y[i] += 0.02 * x[i + 1];
            }
            y[i] += 0.003 * x[(i + 513) % n];
            y[i] += 0.003 * x[(i + n - 513) % n];
        }
        Ok(())
    }
}

fn orthogonality_defect(x: &RealMatrix) -> f64 {
    let mut defect = 0_f64;
    for j in 0..x.cols() {
        for k in 0..=j {
            let o: f64 = (0..x.rows()).map(|i| x[(i, j)] * x[(i, k)]).sum();
            defect = defect.max((o - if j == k { 1. } else { 0. }).abs());
        }
    }
    defect
}

fn switches() -> Vec<(&'static str, ResearchOptions)> {
    let d = ResearchOptions::default();
    vec![
        (
            "narrow-proj",
            ResearchOptions {
                narrow_proj: true,
                ..d
            },
        ),
        (
            "narrow-reduced",
            ResearchOptions {
                narrow_reduced: true,
                ..d
            },
        ),
        (
            "b5-par",
            ResearchOptions {
                b5: true,
                b5_parallel: true,
                ..d
            },
        ),
        (
            "all",
            ResearchOptions {
                ritz_fused: true,
                narrow_proj: true,
                narrow_reduced: true,
                b5: true,
                b5_parallel: true,
                ..d
            },
        ),
    ]
}

/// RITZ-FUSED changes no bit of the solution, with and without B1 and B5.
#[test]
fn ritz_fused_is_bit_identical() {
    check_fused(&quick_cases());
}
/// Both checks on the full fixture set (kinetic 220, Householder spectra with
/// known roots, degenerate triples, clustered groups).
/// `cargo test --release -p dft-solvers --test algebra_stream full_fixture_set -- --ignored`
#[test]
#[ignore = "release-mode fixture set (about 100 s in debug)"]
fn full_fixture_set() {
    let all = cases();
    check_fused(&all);
    check_rounding(&all);
}
fn check_fused(cases: &[Case]) {
    let before = algebra_counters();
    for case in cases {
        for ritz in [false, true] {
            for (b1, b5) in [
                (B1Mode::Off, false),
                (B1Mode::On, false),
                (B1Mode::Off, true),
            ] {
                let base = ResearchOptions {
                    b1,
                    sketch: SketchKind::Count,
                    b5,
                    ..ResearchOptions::default()
                };
                let plain = solve(case, 1, ritz, base);
                let fused = ResearchOptions {
                    ritz_fused: true,
                    ..base
                };
                assert_eq!(
                    solve(case, 1, ritz, fused),
                    plain,
                    "{} ritz={ritz} b1={b1:?} b5={b5}",
                    case.name
                );
                assert_eq!(solve(case, 4, ritz, fused), plain, "{} T=4", case.name);
            }
        }
    }
    let after = algebra_counters();
    assert!(after.fused_iterations > before.fused_iterations);
    assert!(
        after.deferred_hx > before.deferred_hx,
        "a contraction formed H X late"
    );
}

/// Rounding-class switches: roots within 1e-11 Ha of the default path (and of
/// the known spectrum), certified strict returns, orthonormal vectors, work
/// within 10 percent, thread-count identity.
#[test]
fn rounding_class_switches_keep_roots_and_thread_identity() {
    check_rounding(&quick_cases());
}
fn check_rounding(cases: &[Case]) {
    let before = algebra_counters();
    for case in cases {
        for ritz in [false, true] {
            let base = solve(case, 1, ritz, ResearchOptions::default());
            for (label, options) in switches() {
                let fast = solve(case, 1, ritz, options);
                assert_eq!(
                    solve(case, 4, ritz, options),
                    fast,
                    "{} {label} ritz={ritz}: thread identity",
                    case.name
                );
                for (j, (a, b)) in base.eigenvalues.iter().zip(&fast.eigenvalues).enumerate() {
                    assert!(
                        (a - b).abs() <= 1e-11,
                        "{} {label} ritz={ritz} root {j}: {a} vs {b}",
                        case.name
                    );
                    if let Some(oracle) = &case.oracle
                        && j < case.strict
                    {
                        assert!((b - oracle[j]).abs() <= 1e-10, "{} oracle {j}", case.name);
                    }
                }
                assert!(fast.converged.iter().all(|&c| c), "{} {label}", case.name);
                if !ritz {
                    assert!(fast.residuals_certified);
                    assert!(fast.residuals[..case.strict].iter().all(|&r| r <= 1e-10));
                }
                let defect = orthogonality_defect(&fast.eigenvectors);
                assert!(defect <= 1e-12, "{} {label}: defect {defect:e}", case.name);
                let near = |a: usize, b: usize| a.abs_diff(b) * 10 <= a.max(10);
                assert!(
                    near(base.iterations, fast.iterations),
                    "{} {label} ritz={ritz}: iterations {} vs {}",
                    case.name,
                    base.iterations,
                    fast.iterations
                );
                assert!(near(base.operator_applications, fast.operator_applications));
            }
        }
    }
    let after = algebra_counters();
    assert!(after.narrow_projections > before.narrow_projections);
    assert!(after.narrow_reduced > before.narrow_reduced);
    assert!(after.b5_parallel > before.b5_parallel);
}

/// Large basis (parallel dispatch of every kernel): bits equal at pool sizes
/// 1, 2, 4 and 8 for each switch and for all of them together.
/// `cargo test --release -p dft-solvers --test algebra_stream large_basis -- --ignored --nocapture`
#[test]
#[ignore = "release-mode dispatch test (about 20 s in debug)"]
fn large_basis_thread_identity() {
    let n = 4097;
    let diag: Vec<f64> = (0..n).map(|i| i as f64 * 0.004).collect();
    let (bands, strict) = (48, 40);
    let run = |threads: usize, ritz: bool, options: ResearchOptions| {
        let op = Sparse {
            pool: WorkerPool::new(threads).unwrap(),
            n,
        };
        solve_op(&op, &diag, bands, strict, ritz, options)
    };
    let mut all = switches();
    all.push((
        "fused",
        ResearchOptions {
            ritz_fused: true,
            ..ResearchOptions::default()
        },
    ));
    let before = algebra_counters();
    for ritz in [false, true] {
        let base = run(1, ritz, ResearchOptions::default());
        for (label, options) in all.iter().copied() {
            let first = run(1, ritz, options);
            for threads in [2, 4, 8] {
                assert_eq!(run(threads, ritz, options), first, "{label} T={threads}");
            }
            let worst = base
                .eigenvalues
                .iter()
                .zip(&first.eigenvalues)
                .map(|(a, b)| (a - b).abs())
                .fold(0., f64::max);
            assert!(worst <= 1e-11, "{label}: {worst:e}");
            println!(
                "ritz={ritz} {label}: iterations {} vs {}, H {} vs {}, max |d lambda| {worst:e} Ha",
                first.iterations,
                base.iterations,
                first.operator_applications,
                base.operator_applications
            );
        }
    }
    println!("{:?}", algebra_counters_since(before));
}
fn algebra_counters_since(
    before: dft_solvers::real_davidson::AlgebraCounters,
) -> dft_solvers::real_davidson::AlgebraCounters {
    let a = algebra_counters();
    dft_solvers::real_davidson::AlgebraCounters {
        fused_iterations: a.fused_iterations - before.fused_iterations,
        deferred_hx: a.deferred_hx - before.deferred_hx,
        narrow_projections: a.narrow_projections - before.narrow_projections,
        narrow_reduced: a.narrow_reduced - before.narrow_reduced,
        b5_parallel: a.b5_parallel - before.b5_parallel,
    }
}
