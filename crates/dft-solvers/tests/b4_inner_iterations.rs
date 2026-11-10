//! B4 evaluation harness (proposal "inner iterations for narrow blocks"); no
//! engine code. One block Davidson written for this test runs twice on the same
//! fixtures: as the baseline (one preconditioned residual per selected root per
//! full Rayleigh-Ritz step) and as the sequential-directions variant (when at
//! most `NARROW` roots are selected, each builds up to `INNER` directions, each
//! followed by one H application and a small Rayleigh-Ritz over {x_j, w_1..w_i}
//! that refreshes that root's residual; then all directions are appended and one
//! full step follows). Counts are exact: H columns, full steps, single-column
//! steps. The serial `RealDavidson` on the same fixture is reported for scale.
//! Decision rule: docs/performance/requests/b1-b4-b5-davidson.md (B4).
//!
//! `cargo test --release -p dft-solvers --test b4_inner_iterations -- --ignored --nocapture`
use dft_numerics::real::{RealMatrix, symmetric_eigh};
use dft_solvers::{BandTolerances, RealDavidson, RealLinearOperator, SolverError};

const NARROW: usize = 2;

struct Dense {
    n: usize,
    a: Vec<f64>,
}
impl Dense {
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        y.fill(0.);
        for (j, &xj) in x.iter().enumerate() {
            for (y, &a) in y.iter_mut().zip(&self.a[j * self.n..(j + 1) * self.n]) {
                *y += a * xj;
            }
        }
    }
}
impl RealLinearOperator for Dense {
    fn dimension(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        Dense::apply(self, x, y);
        Ok(())
    }
}

fn kinetic_fixture(n: usize, coupling: f64) -> (Dense, Vec<f64>) {
    let kinetic: Vec<f64> = (0..n).map(|i| 0.2 * (i / 2) as f64).collect();
    let mut a = vec![0.; n * n];
    for j in 0..n {
        for i in 0..=j {
            let v = if i == j {
                kinetic[i] + 0.3
            } else {
                coupling * ((i * 31 + j * 17 + 1) as f64).sin() / ((j - i) as f64).sqrt()
            };
            a[j * n + i] = v;
            a[i * n + j] = v;
        }
    }
    (Dense { n, a }, kinetic)
}

/// Known spectrum 0.05 (i / group)^1.3 - 0.3, rotated by Householder
/// reflections (as `spectrum_fixture` in distributed_davidson.rs).
fn spectrum_fixture(n: usize, group: usize, seed: u64) -> (Dense, Vec<f64>) {
    let values: Vec<f64> = (0..n)
        .map(|i| 0.05 * ((i / group) as f64).powf(1.3) - 0.3)
        .collect();
    rotated(n, &values, seed)
}
/// Tail fixture: the group-1 spectrum with roots `edge..` lifted so that root
/// `edge - 1` (last strict root) is separated from root `edge` by only `gap` Ha;
/// that single root converges slowly (a single-column tail by construction).
fn gap_fixture(n: usize, edge: usize, gap: f64, seed: u64) -> (Dense, Vec<f64>) {
    let mut values: Vec<f64> = (0..n).map(|i| 0.05 * (i as f64).powf(1.3) - 0.3).collect();
    let shift = values[edge - 1] + gap - values[edge];
    for v in &mut values[edge..] {
        *v += shift;
    }
    rotated(n, &values, seed)
}
fn rotated(n: usize, values: &[f64], seed: u64) -> (Dense, Vec<f64>) {
    let mut q: Vec<Vec<f64>> = (0..n)
        .map(|j| (0..n).map(|i| f64::from(u8::from(i == j))).collect())
        .collect();
    let mut state = seed;
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
    let rows: Vec<Vec<f64>> = (0..n).map(|k| q.iter().map(|c| c[k]).collect()).collect();
    let mut a = vec![0.; n * n];
    for (k, r) in rows.iter().enumerate() {
        for j in 0..n {
            let s = values[k] * r[j];
            for (dst, &ri) in a[j * n..(j + 1) * n].iter_mut().zip(r) {
                *dst += s * ri;
            }
        }
    }
    for j in 0..n {
        for i in 0..j {
            a[i * n + j] = a[j * n + i];
        }
    }
    let diag = (0..n).map(|i| a[i * n + i]).collect();
    (Dense { n, a }, diag)
}

fn dot(x: &[f64], y: &[f64]) -> f64 {
    x.iter().zip(y).map(|(a, b)| a * b).sum()
}
fn axpy(alpha: f64, x: &[f64], y: &mut [f64]) {
    for (y, &x) in y.iter_mut().zip(x) {
        *y += alpha * x;
    }
}
/// Same kinetic damping as the solver (0.1 Ha floor, sign preserved).
fn precondition(r: &[f64], diag: &[f64], theta: f64) -> Vec<f64> {
    r.iter()
        .zip(diag)
        .map(|(&r, &d)| {
            let delta = d - theta;
            let den = if delta < 0. {
                -delta.abs().max(0.1)
            } else {
                delta.max(0.1)
            };
            -r / den
        })
        .collect()
}
/// Two-pass Gram-Schmidt against `against`; None when dependent.
fn orthonormalize(mut w: Vec<f64>, against: &[&Vec<f64>]) -> Option<Vec<f64>> {
    let start = dot(&w, &w).sqrt();
    if start == 0. {
        return None;
    }
    for _ in 0..2 {
        for q in against {
            let p = dot(&w, q);
            axpy(-p, q, &mut w);
        }
    }
    let len = dot(&w, &w).sqrt();
    if len <= 1e-10 * start {
        return None;
    }
    w.iter_mut().for_each(|x| *x /= len);
    Some(w)
}

#[derive(Debug, Default, Clone)]
struct Stats {
    steps: usize,
    h_columns: usize,
    single_steps: usize,
    inner_rr: usize,
    converged: bool,
    values: Vec<f64>,
}

#[allow(clippy::too_many_arguments)]
fn harness(
    op: &Dense,
    diag: &[f64],
    bands: usize,
    strict: usize,
    tol: f64,
    buffer_tol: f64,
    capacity: usize,
    inner: usize,
) -> Stats {
    let n = op.n;
    let thr = |j: usize| if j < strict { tol } else { buffer_tol };
    let mut stats = Stats::default();
    // Free-electron start: unit vectors at the lowest diagonal entries.
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&i, &j| diag[i].total_cmp(&diag[j]).then(i.cmp(&j)));
    let mut v: Vec<Vec<f64>> = order[..bands]
        .iter()
        .map(|&i| {
            let mut e = vec![0.; n];
            e[i] = 1.;
            e
        })
        .collect();
    let mut hv: Vec<Vec<f64>> = v
        .iter()
        .map(|x| {
            let mut y = vec![0.; n];
            op.apply(x, &mut y);
            y
        })
        .collect();
    stats.h_columns += bands;
    for _ in 0..2000 {
        stats.steps += 1;
        let m = v.len();
        let mut t = RealMatrix::zeros(m, m).unwrap();
        for j in 0..m {
            for i in 0..=j {
                let x = dot(&v[i], &hv[j]);
                t[(i, j)] = x;
                t[(j, i)] = x;
            }
        }
        let (theta, c) = symmetric_eigh(&t).unwrap();
        let mut x = vec![vec![0.; n]; bands];
        let mut hx = vec![vec![0.; n]; bands];
        for j in 0..bands {
            for k in 0..m {
                axpy(c[(k, j)], &v[k], &mut x[j]);
                axpy(c[(k, j)], &hv[k], &mut hx[j]);
            }
        }
        let r: Vec<Vec<f64>> = (0..bands)
            .map(|j| {
                hx[j]
                    .iter()
                    .zip(&x[j])
                    .map(|(h, x)| h - theta[j] * x)
                    .collect()
            })
            .collect();
        let norms: Vec<f64> = r.iter().map(|r| dot(r, r).sqrt()).collect();
        let selected: Vec<usize> = (0..bands).filter(|&j| norms[j] > thr(j)).collect();
        if selected.is_empty() {
            stats.converged = true;
            stats.values = theta[..bands].to_vec();
            return stats;
        }
        if selected.len() == 1 {
            stats.single_steps += 1;
        }
        let mut new_v: Vec<Vec<f64>> = Vec::new();
        let mut new_hv: Vec<Vec<f64>> = Vec::new();
        if inner > 1 && selected.len() <= NARROW {
            for &j in &selected {
                let (mut xj, mut hxj, mut th, mut rj) =
                    (x[j].clone(), hx[j].clone(), theta[j], r[j].clone());
                let first = new_v.len();
                for _ in 0..inner {
                    let w = precondition(&rj, diag, th);
                    let against: Vec<&Vec<f64>> = v.iter().chain(new_v.iter()).collect();
                    let Some(w) = orthonormalize(w, &against) else {
                        break;
                    };
                    let mut hw = vec![0.; n];
                    op.apply(&w, &mut hw);
                    stats.h_columns += 1;
                    new_v.push(w);
                    new_hv.push(hw);
                    // Small Rayleigh-Ritz over {x_j, w_first..}; x_j lies in
                    // span(V), the new directions are orthogonal to it.
                    let q = 1 + new_v.len() - first;
                    let mut g = RealMatrix::zeros(q, q).unwrap();
                    g[(0, 0)] = th;
                    for a in 1..q {
                        let e = dot(&hxj, &new_v[first + a - 1]);
                        g[(0, a)] = e;
                        g[(a, 0)] = e;
                        for b in 1..=a {
                            let e = dot(&new_v[first + b - 1], &new_hv[first + a - 1]);
                            g[(a, b)] = e;
                            g[(b, a)] = e;
                        }
                    }
                    let (mu, y) = symmetric_eigh(&g).unwrap();
                    stats.inner_rr += 1;
                    // Continue root j: the eigenvector with the largest x_j weight.
                    let k = (0..q)
                        .max_by(|&a, &b| y[(0, a)].abs().total_cmp(&y[(0, b)].abs()))
                        .unwrap();
                    let mut nx: Vec<f64> = xj.iter().map(|v| v * y[(0, k)]).collect();
                    let mut nhx: Vec<f64> = hxj.iter().map(|v| v * y[(0, k)]).collect();
                    for a in 1..q {
                        axpy(y[(a, k)], &new_v[first + a - 1], &mut nx);
                        axpy(y[(a, k)], &new_hv[first + a - 1], &mut nhx);
                    }
                    xj = nx;
                    hxj = nhx;
                    th = mu[k];
                    rj = hxj.iter().zip(&xj).map(|(h, x)| h - th * x).collect();
                    if dot(&rj, &rj).sqrt() <= thr(j) {
                        break;
                    }
                }
            }
        } else {
            for &j in &selected {
                let w = precondition(&r[j], diag, theta[j]);
                let against: Vec<&Vec<f64>> = v.iter().chain(new_v.iter()).collect();
                if let Some(w) = orthonormalize(w, &against) {
                    let mut hw = vec![0.; n];
                    op.apply(&w, &mut hw);
                    stats.h_columns += 1;
                    new_v.push(w);
                    new_hv.push(hw);
                }
            }
        }
        if new_v.is_empty() {
            return stats;
        }
        if v.len() + new_v.len() > capacity {
            // Contraction restart to the Ritz vectors (exact H images). The new
            // directions are orthogonal to span(old V), which contains them.
            v = x;
            hv = hx;
        }
        v.extend(new_v);
        hv.extend(new_hv);
    }
    stats
}

#[test]
#[ignore = "B4 evaluation harness: release, --nocapture; prints the decision table"]
fn b4_sequential_directions() {
    let bands = 24;
    let strict = 20;
    let tol = 1e-9;
    let buffer_tol = 1e-5;
    let capacity = 6 * bands;
    let mut cases: Vec<(String, Dense, Vec<f64>)> = Vec::new();
    for (name, coupling) in [("kinetic-0.025", 0.025), ("kinetic-0.05", 0.05)] {
        let (op, d) = kinetic_fixture(600, coupling);
        cases.push((name.into(), op, d));
    }
    for (group, seed) in [(1, 753), (1, 911), (2, 753), (3, 753)] {
        let (op, d) = spectrum_fixture(500, group, seed);
        cases.push((format!("spectrum-g{group}-s{seed}"), op, d));
    }
    // Added after the first run: no fixture above reached the rule's 30 percent
    // single-column filter, so tail fixtures are built on purpose.
    for (gap, seed) in [(2e-3, 753), (5e-4, 753), (5e-4, 911), (1e-4, 753)] {
        let (op, d) = gap_fixture(500, strict, gap, seed);
        cases.push((format!("gap-{gap:e}-s{seed}"), op, d));
    }
    println!(
        "case | RealDavidson steps/H | base steps H single% | inner=2 steps H | inner=3 steps H | inner=4 steps H"
    );
    for (name, op, diag) in &cases {
        let solver = RealDavidson {
            bands,
            tolerance: tol,
            max_iterations: 2000,
            max_subspace: capacity,
            seed: 42,
        };
        let policy = BandTolerances {
            strict_bands: strict,
            buffer_tolerance: buffer_tol,
        };
        let real = solver
            .solve_with_band_tolerances_and_certification(op, None, diag, policy, false)
            .unwrap();
        let base = harness(op, diag, bands, strict, tol, buffer_tol, capacity, 1);
        assert!(base.converged, "{name}: baseline harness converged");
        let mut line = format!(
            "{name} | {}/{} | {} {} {:.0}%",
            real.iterations,
            real.operator_applications,
            base.steps,
            base.h_columns,
            100. * base.single_steps as f64 / base.steps as f64
        );
        for inner in [2, 3, 4] {
            let v = harness(op, diag, bands, strict, tol, buffer_tol, capacity, inner);
            assert!(v.converged, "{name}: inner={inner} converged");
            for (a, b) in base.values.iter().zip(&v.values) {
                assert!((a - b).abs() < 1e-8, "{name}: same roots");
            }
            line += &format!(" | {} {} (rr {})", v.steps, v.h_columns, v.inner_rr);
        }
        println!("{line}");
    }
}
