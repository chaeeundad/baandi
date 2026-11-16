//! C3a (`ResearchOptions::lean_proj_cache`, `BAANDI_LEAN_PROJ_CACHE=on`): the
//! lean Davidson keeps P^T V of its basis columns and leaves P^T X = (P^T V) C
//! of the returned vectors in the workspace. Checks against a direct P^T X of
//! the returned vectors on an operator with several projectors per atom,
//! off-diagonal D blocks and an odd band count: projections and per-band
//! nonlocal energies agree to rounding (energies within 1e-12 Ha), the solve
//! itself is bit-identical with the cache on and off, the cache follows
//! contractions and the kept start (REUSE-PSI), and the cached projections do
//! not depend on the thread count.
use dft_numerics::pool::WorkerPool;
use dft_numerics::real::RealMatrix;
use dft_solvers::real_davidson::{
    RealDavidsonWorkspace, RealEigenSolution, ResearchOptions, lean_counters,
};
use dft_solvers::{
    BandTolerances, DavidsonPreconditioner, RealDavidson, RealLinearOperator, SolverError,
};
use std::cell::RefCell;
use std::sync::Arc;

/// Projector columns per atom (several per atom, unequal blocks).
const ATOMS: [usize; 4] = [3, 2, 4, 3];

/// H = diag(d) + L + sum_atoms P_a D_a P_a^T with dense projector columns,
/// a weak dense symmetric local coupling L and symmetric D_a blocks with
/// off-diagonal entries (Ha).
struct Operator {
    d: Vec<f64>,
    local: Vec<f64>,
    p: Vec<f64>,
    /// Row-major D block of each atom.
    blocks: Vec<Vec<f64>>,
    count: usize,
    /// P^T Y of the latest block application; `expose` turns the hooks off.
    last: RefCell<Vec<f64>>,
    expose: bool,
    pool: Arc<WorkerPool>,
}
impl Operator {
    fn new(n: usize, threads: usize, expose: bool, shift: f64) -> Self {
        let count: usize = ATOMS.iter().sum();
        let d = (0..n)
            .map(|i| 0.12 * (i / 2) as f64 + 0.05 + shift * ((i * 5 + 1) as f64).sin())
            .collect();
        let mut local = vec![0.; n * n];
        for j in 0..n {
            for i in 0..=j {
                let v = if i == j {
                    0.
                } else {
                    0.004 * ((i * 13 + j * 7) as f64).cos() / (1. + (j - i) as f64)
                };
                local[i + j * n] = v;
                local[j + i * n] = v;
            }
        }
        let p = (0..n * count)
            .map(|i| ((i * 7 + 3) as f64 * 0.37).sin() / (1. + (i % n) as f64).sqrt())
            .collect();
        let blocks = ATOMS
            .iter()
            .enumerate()
            .map(|(a, &nc)| {
                let mut b = vec![0.; nc * nc];
                for i in 0..nc {
                    for j in 0..nc {
                        b[i * nc + j] = if i == j {
                            -0.5 + 0.07 * (i + a) as f64
                        } else {
                            0.08 * (((i + j + a) as f64) * 0.9).cos()
                        };
                    }
                }
                b
            })
            .collect();
        Self {
            d,
            local,
            p,
            blocks,
            count,
            last: RefCell::new(Vec::new()),
            expose,
            pool: WorkerPool::new(threads).unwrap(),
        }
    }
    /// P^T X, column-major, `count` rows (naive loops in their own order).
    fn projections(&self, x: &[f64], cols: usize) -> Vec<f64> {
        let n = self.d.len();
        let c = self.count;
        let mut out = vec![0.; c * cols];
        for col in 0..cols {
            let xc = &x[col * n..(col + 1) * n];
            for j in 0..c {
                out[col * c + j] = self.p[j * n..(j + 1) * n]
                    .iter()
                    .zip(xc)
                    .map(|(a, b)| a * b)
                    .sum();
            }
        }
        out
    }
    /// Per column sum_a sum_ij D_a,ij p_i p_j (Ha).
    fn energies(&self, proj: &[f64], cols: usize) -> Vec<f64> {
        (0..cols)
            .map(|col| {
                let p = &proj[col * self.count..(col + 1) * self.count];
                let mut e = 0.;
                let mut start = 0;
                for (block, &nc) in self.blocks.iter().zip(&ATOMS) {
                    for i in 0..nc {
                        for j in 0..nc {
                            e += block[i * nc + j] * p[start + i] * p[start + j];
                        }
                    }
                    start += nc;
                }
                e
            })
            .collect()
    }
}
impl RealLinearOperator for Operator {
    fn dimension(&self) -> usize {
        self.d.len()
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, _x: &[f64], _y: &mut [f64]) -> Result<(), SolverError> {
        unreachable!("block applications only")
    }
    fn apply_block(&self, x: &RealMatrix, y: &mut RealMatrix) -> Result<(), SolverError> {
        let n = self.d.len();
        let cols = x.cols();
        let proj = self.projections(x.as_slice(), cols);
        for (col, (xc, yc)) in x
            .as_slice()
            .chunks(n)
            .zip(y.as_mut_slice().chunks_mut(n))
            .enumerate()
        {
            for i in 0..n {
                let mut v = self.d[i] * xc[i];
                for (k, &xk) in xc.iter().enumerate() {
                    v += self.local[i + k * n] * xk;
                }
                yc[i] = v;
            }
            let p = &proj[col * self.count..(col + 1) * self.count];
            let mut start = 0;
            for (block, &nc) in self.blocks.iter().zip(&ATOMS) {
                for i in 0..nc {
                    let w: f64 = (0..nc).map(|j| block[i * nc + j] * p[start + j]).sum();
                    let column = &self.p[(start + i) * n..(start + i + 1) * n];
                    for (y, &pv) in yc.iter_mut().zip(column) {
                        *y += pv * w;
                    }
                }
                start += nc;
            }
        }
        *self.last.borrow_mut() = proj;
        Ok(())
    }
    fn projector_count(&self) -> usize {
        if self.expose { self.count } else { 0 }
    }
    fn copy_last_projections(&self, m: usize, out: &mut [f64]) -> bool {
        let last = self.last.borrow();
        if !self.expose || out.len() != m * self.count || last.len() != out.len() {
            return false;
        }
        out.copy_from_slice(&last[..]);
        true
    }
}

fn solver(bands: usize) -> RealDavidson {
    RealDavidson {
        bands,
        tolerance: 1e-10,
        max_iterations: 400,
        max_subspace: 0,
        seed: 11,
    }
}
const POLICY: BandTolerances = BandTolerances {
    strict_bands: 5,
    buffer_tolerance: 1e-7,
};
fn options(cache: bool) -> ResearchOptions {
    ResearchOptions {
        lean: true,
        lean_proj_cache: cache,
        ..ResearchOptions::default()
    }
}
fn solve(
    op: &Operator,
    bands: usize,
    cache: bool,
    capacity_per_band: usize,
    initial: Option<&RealMatrix>,
    ws: &mut Option<RealDavidsonWorkspace>,
) -> RealEigenSolution {
    solver(bands)
        .solve_research(
            op,
            initial,
            &op.d,
            POLICY,
            true,
            false,
            capacity_per_band,
            ws,
            DavidsonPreconditioner::Default,
            options(cache),
        )
        .unwrap()
}

/// Largest |cached - direct| of the projections and of the band energies (Ha),
/// and whether every projection entry is bit-identical.
fn compare(op: &Operator, sol: &RealEigenSolution, cached: &[f64]) -> (f64, f64, bool) {
    let bands = sol.eigenvectors.cols();
    let direct = op.projections(sol.eigenvectors.as_slice(), bands);
    assert_eq!(cached.len(), direct.len());
    let dp = cached
        .iter()
        .zip(&direct)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max);
    let ec = op.energies(cached, bands);
    let ed = op.energies(&direct, bands);
    assert!(
        ed.iter().any(|e| e.abs() > 1e-4),
        "nonlocal energies {ed:?}"
    );
    let de = ec
        .iter()
        .zip(&ed)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max);
    (dp, de, cached == &direct[..])
}

#[test]
fn returned_projections_match_direct_product() {
    let n = 230;
    for bands in [7, 9] {
        for capacity_per_band in [2, 6] {
            let mut reference: Option<Vec<f64>> = None;
            for threads in [1, 3] {
                let op = Operator::new(n, threads, true, 0.);
                let before = lean_counters();
                let mut ws = None;
                let cached = solve(&op, bands, true, capacity_per_band, None, &mut ws);
                let after = lean_counters();
                let p = ws
                    .as_ref()
                    .unwrap()
                    .returned_projections()
                    .expect("cache complete at return")
                    .to_vec();
                assert!(after.proj_cache_returns > before.proj_cache_returns);
                let (dp, de, bitwise) = compare(&op, &cached, &p);
                println!(
                    "bands {bands} cap/band {capacity_per_band} T={threads}: max |dP| {dp:e}, \
                     max |dE| {de:e} Ha, bit-identical {bitwise}, contractions {}",
                    cached.contracted_restarts
                );
                assert!(dp <= 1e-13, "projections differ by {dp:e}");
                assert!(de <= 1e-12, "band energies differ by {de:e} Ha");
                if capacity_per_band == 2 {
                    assert!(cached.contracted_restarts > 0, "fixture must contract");
                }
                // The cache does not touch the solve.
                let plain = solve(&op, bands, false, capacity_per_band, None, &mut None);
                assert_eq!(cached.eigenvalues, plain.eigenvalues);
                assert_eq!(
                    cached.eigenvectors.as_slice(),
                    plain.eigenvectors.as_slice()
                );
                assert_eq!(cached.iterations, plain.iterations);
                assert_eq!(cached.operator_applications, plain.operator_applications);
                // Thread-count independence of the cached projections.
                match &reference {
                    None => reference = Some(p),
                    Some(r) => assert_eq!(r, &p),
                }
            }
        }
    }
}

#[test]
fn no_projections_without_switch_or_support() {
    let op = Operator::new(200, 2, true, 0.);
    let mut ws = None;
    solve(&op, 7, false, 6, None, &mut ws);
    assert!(ws.as_ref().unwrap().returned_projections().is_none());
    // A cached solve, then a solve with the cache off on the same workspace:
    // the stale projections are not reported.
    solve(&op, 7, true, 6, None, &mut ws);
    assert!(ws.as_ref().unwrap().returned_projections().is_some());
    solve(&op, 7, false, 6, None, &mut ws);
    assert!(ws.as_ref().unwrap().returned_projections().is_none());
    // A residual (strict) solve takes the classic path: nothing returned.
    solve(&op, 7, true, 6, None, &mut ws);
    solver(7)
        .solve_research(
            &op,
            None,
            &op.d,
            POLICY,
            false,
            true,
            6,
            &mut ws,
            DavidsonPreconditioner::Default,
            options(true),
        )
        .unwrap();
    assert!(ws.as_ref().unwrap().returned_projections().is_none());
    // An operator without projector products.
    let hidden = Operator::new(200, 2, false, 0.);
    let mut ws = None;
    solve(&hidden, 7, true, 6, None, &mut ws);
    assert!(ws.as_ref().unwrap().returned_projections().is_none());
}

/// REUSE-PSI kept start: the kept X carries its P^T X into the next solve
/// (operator change in the local part only, projectors unchanged).
#[test]
fn kept_start_carries_projections() {
    let n = 230;
    let bands = 7;
    let op = Operator::new(n, 2, true, 0.);
    let mut ws = None;
    solver(bands)
        .prepare_workspace(&op, 6, &mut ws)
        .unwrap()
        .keep_returned_products(true)
        .unwrap();
    let first = solve(&op, bands, true, 6, None, &mut ws);
    let x0 = first.eigenvectors.clone();
    // The next operator differs by a diagonal shift (same pool, so the
    // workspace is kept); bring H X to it.
    let op2 = Operator {
        pool: op.pool.clone(),
        ..Operator::new(n, 1, true, 0.002)
    };
    {
        let (x, hx) = ws.as_mut().unwrap().kept_start_mut().expect("kept");
        for j in 0..bands {
            for i in 0..n {
                hx[(i, j)] += (op2.d[i] - op.d[i]) * x[(i, j)];
            }
        }
    }
    assert!(ws.as_mut().unwrap().arm_kept_start());
    let before = lean_counters();
    let second = solve(&op2, bands, true, 6, Some(&x0), &mut ws);
    let after = lean_counters();
    assert_eq!(after.kept_starts, before.kept_starts + 1);
    assert_eq!(
        after.proj_cache_kept_loads,
        before.proj_cache_kept_loads + 1
    );
    let p = ws
        .as_ref()
        .unwrap()
        .returned_projections()
        .expect("cache complete after a kept start")
        .to_vec();
    let (dp, de, bitwise) = compare(&op2, &second, &p);
    println!("kept start: max |dP| {dp:e}, max |dE| {de:e} Ha, bit-identical {bitwise}");
    assert!(dp <= 1e-13 && de <= 1e-12, "{dp:e} {de:e}");
}
