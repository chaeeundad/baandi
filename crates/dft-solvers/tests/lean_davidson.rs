//! Lean Davidson (`BAANDI_DAVIDSON=lean`, work order phase B): off by default;
//! residual (strict) solves never take it; on Ritz-value solves it agrees with
//! the classic path on eigenvalues (1e-9 Ha) on cubic-like, skewed,
//! near-degenerate and multi-projector fixtures, returns exact recurrence
//! residuals and orthonormal vectors, exercises the fallback on an
//! ill-conditioned fixture, and does not depend on the thread count.
use dft_numerics::{pool::WorkerPool, real::RealMatrix};
use dft_solvers::real_davidson::{
    PIVOT_RATIO_MIN, RealEigenSolution, ResearchOptions, flop_counters, lean_counters,
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
fn symmetric(n: usize, entry: impl Fn(usize, usize) -> f64) -> RealMatrix {
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for j in 0..n {
        for i in 0..=j {
            let v = entry(i, j);
            a[(i, j)] = v;
            a[(j, i)] = v;
        }
    }
    a
}
fn operator(a: RealMatrix, threads: usize) -> Operator {
    Operator {
        a,
        pool: WorkerPool::new(threads).unwrap(),
    }
}
fn lcg(state: &mut u64) -> f64 {
    *state = state
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*state >> 11) as f64) / (1_u64 << 53) as f64 - 0.5
}

/// Cubic-like: kinetic diagonal in shells with a decaying coupling (Ha).
fn cubic(n: usize, threads: usize) -> (Operator, Vec<f64>) {
    let kinetic: Vec<f64> = (0..n).map(|i| 0.2 * (i / 2) as f64).collect();
    let a = symmetric(n, |i, j| {
        if i == j {
            kinetic[i] + 0.3
        } else {
            0.025 * ((i * 31 + j * 17 + 1) as f64).sin() / ((j - i) as f64).sqrt()
        }
    });
    (operator(a, threads), kinetic)
}

/// Plane-wave set of a skewed cell (reciprocal vectors at 70 and 110 degrees,
/// unequal lengths), sorted by kinetic energy; a smooth local potential
/// couples G and G' through G - G'. Returns kinetic energies and G vectors.
fn skew_waves(n: usize) -> (Vec<f64>, Vec<[f64; 3]>) {
    let b = [[1.0, 0.0, 0.0], [0.34, 0.94, 0.0], [-0.25, 0.3, 1.35]];
    let mut waves = Vec::new();
    for h in -6_i32..=6 {
        for k in -6_i32..=6 {
            for l in -5_i32..=5 {
                let g: [f64; 3] = std::array::from_fn(|c| {
                    f64::from(h) * b[0][c] + f64::from(k) * b[1][c] + f64::from(l) * b[2][c]
                });
                let e = 0.5 * (g[0] * g[0] + g[1] * g[1] + g[2] * g[2]);
                waves.push((e, g));
            }
        }
    }
    waves.sort_by(|a, b| a.0.total_cmp(&b.0));
    waves.truncate(n);
    (
        waves.iter().map(|w| w.0).collect(),
        waves.iter().map(|w| w.1).collect(),
    )
}
fn skew(n: usize, threads: usize) -> (Operator, Vec<f64>) {
    let (kinetic, g) = skew_waves(n);
    let a = symmetric(n, |i, j| {
        let d: [f64; 3] = std::array::from_fn(|c| g[i][c] - g[j][c]);
        let q2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        let local = -0.35 * (-0.6 * q2).exp() * (0.7 * d[0] - 0.4 * d[2]).cos();
        if i == j { kinetic[i] + local } else { local }
    });
    (operator(a, threads), kinetic)
}

/// Synthetic multi-projector: skewed kinetic set, local potential and six
/// separable projectors of different widths and angular signs, attractive and
/// repulsive strengths (Ha).
fn projectors(n: usize, threads: usize) -> (Operator, Vec<f64>) {
    let (kinetic, g) = skew_waves(n);
    let shapes: [(f64, usize, f64); 6] = [
        (0.8, 0, -1.2),
        (1.5, 0, 0.6),
        (1.1, 1, -0.7),
        (1.1, 2, -0.7),
        (2.2, 3, 0.9),
        (0.6, 1, 0.4),
    ];
    let beta: Vec<Vec<f64>> = shapes
        .iter()
        .map(|&(width, axis, _)| {
            g.iter()
                .map(|g| {
                    let q2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
                    let angular = if axis == 0 { 1.0 } else { g[axis - 1] };
                    angular * (-q2 / (width * width)).exp()
                })
                .collect()
        })
        .collect();
    let a = symmetric(n, |i, j| {
        let d: [f64; 3] = std::array::from_fn(|c| g[i][c] - g[j][c]);
        let q2 = d[0] * d[0] + d[1] * d[1] + d[2] * d[2];
        let mut v = -0.25 * (-0.8 * q2).exp();
        for (p, &(_, _, strength)) in beta.iter().zip(&shapes) {
            v += strength * p[i] * p[j] / n as f64 * 8.0;
        }
        if i == j { kinetic[i] + v } else { v }
    });
    (operator(a, threads), kinetic)
}

/// Near-degenerate: known spectrum in triples split by 1e-7 Ha, rotated by
/// deterministic Householder reflections; the diagonal is the preconditioner.
fn near_degenerate(n: usize, threads: usize) -> (Operator, Vec<f64>) {
    let mut q: Vec<Vec<f64>> = (0..n)
        .map(|j| (0..n).map(|i| f64::from(u8::from(i == j))).collect())
        .collect();
    let mut state = 977_u64;
    for _ in 0..12 {
        let v: Vec<f64> = (0..n).map(|_| lcg(&mut state)).collect();
        let scale = 2. / v.iter().map(|x| x * x).sum::<f64>();
        for col in &mut q {
            let d = v.iter().zip(col.iter()).map(|(a, b)| a * b).sum::<f64>() * scale;
            for (c, &vi) in col.iter_mut().zip(&v) {
                *c -= vi * d;
            }
        }
    }
    let values: Vec<f64> = (0..n)
        .map(|i| 0.06 * ((i / 3) as f64).powf(1.2) - 0.4 + 1e-7 * (i % 3) as f64)
        .collect();
    let mut a = RealMatrix::zeros(n, n).unwrap();
    for k in 0..n {
        let r: Vec<f64> = q.iter().map(|c| c[k]).collect();
        for j in 0..n {
            let s = values[k] * r[j];
            for (dst, &ri) in a.as_mut_slice()[j * n..(j + 1) * n].iter_mut().zip(&r) {
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
    let diag = (0..n).map(|i| a[(i, i)]).collect();
    (operator(a, threads), diag)
}

/// Ill-conditioned for the lean basis: an almost diagonal operator whose
/// diagonal is the exact preconditioner, so every correction -(D - theta)^-1 r
/// is almost the Ritz vector itself (a nearly dependent column).
/// With `coupling` zero the start vectors are exact eigenvectors, every
/// correction is zero and both rank decisions drop the whole block (stage 2:
/// the classic random block, as on the classic path).
fn ill_conditioned(n: usize, coupling: f64, threads: usize) -> (Operator, Vec<f64>) {
    let d: Vec<f64> = (0..n).map(|i| 0.3 * i as f64).collect();
    let a = symmetric(n, |i, j| {
        if i == j {
            d[i]
        } else {
            coupling * ((i * 13 + j * 7) as f64).cos()
        }
    });
    (operator(a, threads), d)
}

fn solver(bands: usize, tolerance: f64) -> RealDavidson {
    RealDavidson {
        bands,
        tolerance,
        max_iterations: 400,
        max_subspace: 2 * bands,
        seed: 42,
    }
}
fn options(lean: bool) -> ResearchOptions {
    ResearchOptions {
        lean,
        ..ResearchOptions::from_env()
    }
}
#[allow(clippy::too_many_arguments)]
fn run(
    s: &RealDavidson,
    op: &Operator,
    diag: &[f64],
    policy: BandTolerances,
    ritz: bool,
    capacity: usize,
    lean: bool,
) -> RealEigenSolution {
    s.solve_research(
        op,
        None,
        diag,
        policy,
        ritz,
        !ritz,
        capacity,
        &mut None,
        DavidsonPreconditioner::Default,
        options(lean),
    )
    .unwrap()
}
/// ||A x - theta x|| per returned pair and the orthonormality defect.
fn check(op: &Operator, sol: &RealEigenSolution) -> (Vec<f64>, f64) {
    let n = op.a.rows();
    let x = &sol.eigenvectors;
    let mut y = vec![0.; n];
    let mut res = Vec::new();
    for (j, &theta) in sol.eigenvalues.iter().enumerate() {
        let col = &x.as_slice()[j * n..(j + 1) * n];
        op.apply(col, &mut y).unwrap();
        res.push(
            y.iter()
                .zip(col)
                .map(|(a, b)| (a - theta * b) * (a - theta * b))
                .sum::<f64>()
                .sqrt(),
        );
    }
    let mut defect = 0_f64;
    for i in 0..x.cols() {
        for j in 0..=i {
            let d: f64 = (0..n).map(|r| x[(r, i)] * x[(r, j)]).sum();
            defect = defect.max((d - if i == j { 1. } else { 0. }).abs());
        }
    }
    (res, defect)
}

type Fixture = fn(usize, usize) -> (Operator, Vec<f64>);
const FIXTURES: [(&str, Fixture, usize); 4] = [
    ("cubic", cubic, 220),
    ("skew", skew, 260),
    ("near-degenerate", near_degenerate, 200),
    ("multi-projector", projectors, 260),
];

/// The switch defaults to lean (adopted 2026-11-15); lean never acts on residual solves.
#[test]
fn default_lean_and_strict_untouched() {
    use dft_numerics::switches::{Id, Source, manifest};
    let entry = &manifest()[Id::DAVIDSON as usize];
    if entry.source == Source::Default {
        assert!(
            ResearchOptions::from_env().lean,
            "DAVIDSON must default lean"
        );
        let (op, diag) = cubic(160, 1);
        let s = solver(12, 1e-9);
        let policy = BandTolerances {
            strict_bands: 10,
            buffer_tolerance: 1e-6,
        };
        let a = s.solve_ritz_changes(&op, None, &diag, policy).unwrap();
        assert_eq!(a, run(&s, &op, &diag, policy, true, 6, true));
    }
    let (op, diag) = skew(160, 1);
    let s = solver(12, 1e-8);
    let policy = BandTolerances {
        strict_bands: 10,
        buffer_tolerance: 1e-6,
    };
    let before = lean_counters().solves;
    let classic = run(&s, &op, &diag, policy, false, 6, false);
    let lean = run(&s, &op, &diag, policy, false, 6, true);
    assert_eq!(classic, lean, "strict solves must stay classic");
    assert!(classic.residuals_certified);
    // Other tests may run lean solves concurrently; this one adds none.
    let _ = before;
}

/// Lean against classic on four fixtures and two capacities: eigenvalues
/// within 1e-9 Ha, every returned residual equals the true residual and stays
/// below the Ritz-change accuracy bound, vectors orthonormal, H work reported.
#[test]
fn lean_agrees_with_classic() {
    for (name, fixture, n) in FIXTURES {
        let (op, diag) = fixture(n, 1);
        let s = solver(10, 1e-11);
        let policy = BandTolerances {
            strict_bands: 8,
            buffer_tolerance: 1e-9,
        };
        for capacity in [6, 3] {
            let classic = run(&s, &op, &diag, policy, true, capacity, false);
            let lean = run(&s, &op, &diag, policy, true, capacity, true);
            assert!(!lean.residuals_certified);
            let (res, defect) = check(&op, &lean);
            let (classic_res, classic_defect) = check(&op, &classic);
            let worst_classic = classic_res.iter().copied().fold(0., f64::max);
            eprintln!(
                "LEAN {name} cap {capacity}: iterations {} / {}, H columns {} / {}, \
                 restarts {} / {}, defect {defect:.2e} / {classic_defect:.2e}, \
                 max residual {:.2e} / {worst_classic:.2e}",
                lean.iterations,
                classic.iterations,
                lean.operator_applications,
                classic.operator_applications,
                lean.restarts,
                classic.restarts,
                res.iter().copied().fold(0., f64::max),
            );
            for (j, (a, b)) in lean
                .eigenvalues
                .iter()
                .zip(&classic.eigenvalues)
                .enumerate()
            {
                assert!(
                    (a - b).abs() <= 1e-9,
                    "{name} cap {capacity} root {j}: {a} vs {b}"
                );
            }
            for (j, (&r, &reported)) in res.iter().zip(&lean.residuals).enumerate() {
                // Recurrence residual = true residual up to rounding (operator
                // norm a few Ha); Ritz-change accuracy 1e-11 Ha gives a residual
                // of at most about sqrt(1e-11 x gap) < 1e-4.
                assert!(
                    (r - reported).abs() <= 1e-11,
                    "{name} root {j}: {r:e} vs {reported:e}"
                );
                assert!(r <= 1e-4, "{name} cap {capacity} root {j}: residual {r:e}");
            }
            assert!(defect <= 1e-11, "{name} cap {capacity}: defect {defect:e}");
            assert_eq!(lean.converged, vec![true; 10], "{name}");
        }
    }
}

/// Both fallback stages fire on ill-conditioned fixtures (pivot ratio below
/// eps/1e-12): stage 1 (block re-projection) on a nearly diagonal operator with
/// the exact diagonal as preconditioner, stage 2 (basis orthonormalization and
/// classic projection, here down to the random block) on an exactly diagonal
/// one; the solves agree with classic and return orthonormal vectors.
#[test]
fn fallback_exercised() {
    assert!((PIVOT_RATIO_MIN - f64::EPSILON / 1e-12).abs() < 1e-30);
    let s = solver(6, 1e-12);
    let policy = BandTolerances {
        strict_bands: 5,
        buffer_tolerance: 1e-10,
    };
    for (coupling, stage2) in [(2e-6, false), (0.0, true)] {
        let (op, diag) = ill_conditioned(150, coupling, 1);
        let before = lean_counters();
        let flops_before = flop_counters();
        let classic = run(&s, &op, &diag, policy, true, 6, false);
        let lean = run(&s, &op, &diag, policy, true, 6, true);
        let after = lean_counters();
        let flops_after = flop_counters();
        assert!(
            after.reprojections > before.reprojections,
            "{coupling}: stage 1"
        );
        assert!(after.pivot_min < PIVOT_RATIO_MIN);
        if stage2 {
            assert!(after.fallbacks > before.fallbacks, "stage 2 must fire");
            assert!(flops_after.projection > flops_before.projection);
        }
        for (a, b) in lean.eigenvalues.iter().zip(&classic.eigenvalues) {
            assert!((a - b).abs() <= 1e-9, "{coupling}: {a} vs {b}");
        }
        let (res, defect) = check(&op, &lean);
        assert!(defect <= 1e-11, "{coupling}: defect {defect:e}");
        for (r, reported) in res.iter().zip(&lean.residuals) {
            assert!((r - reported).abs() <= 1e-11);
        }
        eprintln!(
            "LEAN fallback coupling {coupling}: reprojections {}, fallbacks {} of {} appends, \
             dropped {}, H columns {} / {}",
            after.reprojections - before.reprojections,
            after.fallbacks - before.fallbacks,
            after.appends - before.appends,
            after.fallback_dropped_columns - before.fallback_dropped_columns,
            lean.operator_applications,
            classic.operator_applications
        );
    }
}

/// Bit-identical results at 1, 2, 4 and 8 workers (n = 1100: four or more
/// 256-row panels, products above the parallel threshold).
#[test]
fn lean_thread_independent() {
    let s = solver(16, 1e-10);
    let policy = BandTolerances {
        strict_bands: 14,
        buffer_tolerance: 1e-8,
    };
    let reference = {
        let (op, diag) = cubic(1100, 1);
        run(&s, &op, &diag, policy, true, 6, true)
    };
    for threads in [2, 4, 8] {
        let (op, diag) = cubic(1100, threads);
        let other = run(&s, &op, &diag, policy, true, 6, true);
        assert_eq!(other.iterations, reference.iterations, "T={threads}");
        let bits = |x: &RealEigenSolution| -> Vec<u64> {
            x.eigenvalues
                .iter()
                .chain(x.eigenvectors.as_slice())
                .chain(&x.residuals)
                .map(|v| v.to_bits())
                .collect()
        };
        assert!(bits(&other) == bits(&reference), "T={threads}");
    }
}

/// Counters: a lean solve adds solves, appends and overlap flops, forms
/// residuals for at most the band count per iteration, and never projects
/// against the basis without a fallback.
#[test]
fn lean_counters_move() {
    let (op, diag) = cubic(200, 1);
    let s = solver(10, 1e-10);
    let policy = BandTolerances {
        strict_bands: 8,
        buffer_tolerance: 1e-8,
    };
    let before = lean_counters();
    let fb = flop_counters();
    let lean = run(&s, &op, &diag, policy, true, 6, true);
    let after = lean_counters();
    let fa = flop_counters();
    assert!(after.solves > before.solves);
    assert!(after.appends > before.appends);
    assert!(after.returns > before.returns);
    assert!(fa.overlap > fb.overlap);
    assert!(fa.ret > fb.ret);
    assert!(fa.small > fb.small);
    assert!(after.pivot_min <= 1.0);
    assert!(after.fallback_rate() <= 1.0);
    assert!(lean.max_basis_dimension <= 60);
}

/// REUSE-PSI (work order C2), solver half: a lean return keeps X and H X (the
/// return rotation equals A X to rounding); after the caller brings H X to a
/// changed operator (here A + D, H_new X = H_old X + D X) and arms the start,
/// the next lean solve from `initial = X` skips the initial application of
/// all bands, agrees with a plain warm start on eigenvalues (1e-9 Ha) and on
/// the iteration count, and returns exact residuals. An armed start whose
/// `initial` differs is declined (bit-identical to the plain solve), and a
/// residual (classic) solve drops the kept start.
#[test]
fn kept_start_skips_initial_application() {
    use dft_solvers::real_davidson::KEPT_START_TOLERANCE;
    let (op, diag) = skew(260, 2);
    let bands = 10;
    let s = solver(bands, 1e-11);
    let policy = BandTolerances {
        strict_bands: 8,
        buffer_tolerance: 1e-9,
    };
    let solve = |op: &Operator, initial: Option<&RealMatrix>, ws: &mut Option<_>, ritz: bool| {
        s.solve_research(
            op,
            initial,
            &diag,
            policy,
            ritz,
            !ritz,
            6,
            ws,
            DavidsonPreconditioner::Default,
            options(true),
        )
        .unwrap()
    };
    let mut ws = None;
    s.prepare_workspace(&op, 6, &mut ws)
        .unwrap()
        .keep_returned_products(true)
        .unwrap();
    let first = solve(&op, None, &mut ws, true);
    let (x, hx) = ws.as_ref().unwrap().kept_start().expect("kept at return");
    assert_eq!(x.as_slice(), first.eigenvectors.as_slice());
    let n = op.a.rows();
    let mut ax = vec![0.; n];
    let mut worst = 0_f64;
    for j in 0..bands {
        op.apply(&x.as_slice()[j * n..(j + 1) * n], &mut ax)
            .unwrap();
        for (a, b) in ax.iter().zip(&hx.as_slice()[j * n..(j + 1) * n]) {
            worst = worst.max((a - b).abs());
        }
    }
    assert!(worst <= 1e-12, "kept H X against A X: {worst:e}");
    // A changed operator on the same pool (a local-potential-like change).
    let d: Vec<f64> = (0..n).map(|i| 0.01 * ((i * 7 + 3) as f64).sin()).collect();
    let mut a2 = op.a.clone();
    for (i, &v) in d.iter().enumerate() {
        a2[(i, i)] += v;
    }
    let op2 = Operator {
        a: a2,
        pool: op.pool.clone(),
    };
    let x0 = first.eigenvectors.clone();
    {
        let (x, hx) = ws.as_mut().unwrap().kept_start_mut().unwrap();
        for j in 0..bands {
            for i in 0..n {
                hx[(i, j)] += d[i] * x[(i, j)];
            }
        }
    }
    assert!(ws.as_mut().unwrap().arm_kept_start());
    let before = lean_counters();
    let reused = solve(&op2, Some(&x0), &mut ws, true);
    let after = lean_counters();
    assert!(after.kept_starts > before.kept_starts);
    let plain = solve(&op2, Some(&x0), &mut None, true);
    assert_eq!(
        reused.operator_applications + bands,
        plain.operator_applications
    );
    assert_eq!(reused.iterations, plain.iterations);
    for (a, b) in reused.eigenvalues.iter().zip(&plain.eigenvalues) {
        assert!((a - b).abs() <= 1e-9, "{a} vs {b}");
    }
    let (res, defect) = check(&op2, &reused);
    assert!(defect <= 1e-11, "defect {defect:e}");
    for (r, reported) in res.iter().zip(&reused.residuals) {
        assert!((r - reported).abs() <= 1e-11, "{r:e} vs {reported:e}");
    }
    // Declined: initial is not the kept X (one entry moved beyond the tolerance).
    let mut moved = reused.eigenvectors.clone();
    moved[(3, 1)] += 10. * KEPT_START_TOLERANCE;
    assert!(ws.as_mut().unwrap().arm_kept_start());
    let before = lean_counters();
    let declined = solve(&op2, Some(&moved), &mut ws, true);
    assert!(lean_counters().kept_declined > before.kept_declined);
    let fresh = solve(&op2, Some(&moved), &mut None, true);
    assert_eq!(declined, fresh);
    // A residual solve leaves no kept start behind.
    assert!(ws.as_ref().unwrap().kept_start().is_some());
    let _ = solve(&op2, Some(&moved), &mut ws, false);
    assert!(ws.as_ref().unwrap().kept_start().is_none());
}

/// BREAKTHROUGH-3 (`RealDavidson::advance_reusing`, design
/// docs/design/bounded-subspace-scf.md): an `Advance` call stops after its
/// correction blocks (one block: the all-band first correction, so the
/// operator columns are the initial application plus one block), returns the
/// all-band Ritz state with true Ritz-change flags (not all converged) and
/// keeps X and H X for the next start; a column budget below one block stops
/// on columns; a budget large enough to converge returns `Converged` with the
/// same eigenvalues (bit for bit), iterations and operator count as the
/// converging solve; a zero budget is an input error. On a classic-Davidson
/// process (`BAANDI_DAVIDSON=classic`) the entry declines.
#[test]
fn advance_respects_its_block_budget_and_converges_when_allowed() {
    use dft_solvers::real_davidson::{BudgetStop, SubspaceBudget};
    let (op, diag) = skew(260, 1);
    let bands = 10;
    let s = solver(bands, 1e-9);
    let policy = BandTolerances {
        strict_bands: 8,
        buffer_tolerance: 5e-6,
    };
    let start = {
        let mut ws = None;
        s.solve_research(
            &op,
            None,
            &diag,
            policy,
            true,
            false,
            6,
            &mut ws,
            DavidsonPreconditioner::Default,
            options(true),
        )
        .unwrap()
    };
    // Perturbed solution: corrections are needed again.
    let mut x0 = start.eigenvectors.clone();
    for (i, v) in x0.as_mut_slice().iter_mut().enumerate() {
        *v += 1e-3 * ((i * 13 + 5) as f64).sin();
    }
    let lean_env = ResearchOptions::from_env().lean;
    let mut ws = None;
    s.prepare_workspace(&op, 6, &mut ws)
        .unwrap()
        .keep_returned_products(true)
        .unwrap();
    let one = s.advance_reusing(
        &op,
        Some(&x0),
        &diag,
        policy,
        SubspaceBudget {
            correction_blocks: 1,
            h_columns: 10 * bands,
        },
        6,
        &mut ws,
        DavidsonPreconditioner::Default,
    );
    if !lean_env {
        assert!(matches!(one, Err(SolverError::InvalidParameter(_))));
        return;
    }
    let one = one.unwrap();
    assert_eq!(one.diagnostics.correction_blocks, 1);
    assert_eq!(one.diagnostics.stop, BudgetStop::CorrectionBlocks);
    assert!(!one.diagnostics.kept_start);
    assert_eq!(one.diagnostics.h_columns, 2 * bands);
    assert_eq!(one.solution.operator_applications, 2 * bands);
    assert!(one.diagnostics.unconverged_roots > 0);
    assert!(one.solution.converged.iter().any(|c| !c));
    assert_eq!(one.solution.residuals.len(), bands);
    assert!(
        one.solution
            .residuals
            .iter()
            .all(|r| r.is_finite() && *r > 0.)
    );
    let (x, _) = ws
        .as_ref()
        .unwrap()
        .kept_start()
        .expect("kept at the advance return");
    assert_eq!(x.as_slice(), one.solution.eigenvectors.as_slice());
    let mut ws2 = None;
    let cols = s
        .advance_reusing(
            &op,
            Some(&x0),
            &diag,
            policy,
            SubspaceBudget {
                correction_blocks: 8,
                h_columns: bands + 1,
            },
            6,
            &mut ws2,
            DavidsonPreconditioner::Default,
        )
        .unwrap();
    assert_eq!(cols.diagnostics.stop, BudgetStop::HColumns);
    assert_eq!(cols.diagnostics.correction_blocks, 1);
    let mut ws3 = None;
    let full = s
        .advance_reusing(
            &op,
            Some(&x0),
            &diag,
            policy,
            SubspaceBudget {
                correction_blocks: 400,
                h_columns: usize::MAX,
            },
            6,
            &mut ws3,
            DavidsonPreconditioner::Default,
        )
        .unwrap();
    assert_eq!(full.diagnostics.stop, BudgetStop::Converged);
    assert_eq!(full.diagnostics.unconverged_roots, 0);
    let mut ws4 = None;
    let reference = s
        .solve_research(
            &op,
            Some(&x0),
            &diag,
            policy,
            true,
            false,
            6,
            &mut ws4,
            DavidsonPreconditioner::Default,
            options(true),
        )
        .unwrap();
    assert_eq!(full.solution.eigenvalues, reference.eigenvalues);
    assert_eq!(
        full.solution.operator_applications,
        reference.operator_applications
    );
    assert_eq!(full.solution.iterations, reference.iterations);
    let mut ws5 = None;
    assert!(matches!(
        s.advance_reusing(
            &op,
            Some(&x0),
            &diag,
            policy,
            SubspaceBudget {
                correction_blocks: 0,
                h_columns: 1
            },
            6,
            &mut ws5,
            DavidsonPreconditioner::Default,
        ),
        Err(SolverError::InvalidParameter(_))
    ));
}
