//! B6 (`ResearchOptions::b6`, `BAANDI_B6=on`): the return certificate of the
//! real Davidson takes P^T X from the cached products of the basis columns,
//! (P^T V) C, instead of a fresh projector product. Checks: same eigenpairs bit
//! for bit (the certificate only recomputes residuals), residuals equal to
//! rounding, identical work counts, the cached path is taken (also across
//! contracted restarts), and thread-count independence.
use dft_numerics::pool::WorkerPool;
use dft_numerics::real::RealMatrix;
use dft_solvers::real_davidson::{RealEigenSolution, ResearchOptions};
use dft_solvers::{
    BandTolerances, DavidsonPreconditioner, RealDavidson, RealLinearOperator, SolverError,
};
use std::cell::{Cell, RefCell};
use std::sync::Arc;

/// H = diag(d) + P D P^T with dense projector columns P (n x count).
struct Operator {
    d: Vec<f64>,
    p: Vec<f64>,
    coupling: Vec<f64>,
    count: usize,
    last: RefCell<Vec<f64>>,
    cached_calls: Cell<usize>,
    pool: Arc<WorkerPool>,
}
impl Operator {
    fn new(n: usize, count: usize, threads: usize) -> Self {
        let d = (0..n).map(|i| 0.15 * (i / 2) as f64 + 0.1).collect();
        let p = (0..n * count)
            .map(|i| ((i * 7 + 3) as f64 * 0.37).sin() / (1. + (i % n) as f64).sqrt())
            .collect();
        let coupling = (0..count * count)
            .map(|i| {
                let (r, c) = (i / count, i % count);
                if r == c {
                    -0.4 + 0.05 * r as f64
                } else {
                    0.02 * ((r + c) as f64).cos()
                }
            })
            .collect();
        Self {
            d,
            p,
            coupling,
            count,
            last: RefCell::new(vec![0.; count * n]),
            cached_calls: Cell::new(0),
            pool: WorkerPool::new(threads).unwrap(),
        }
    }
    fn add_nonlocal(&self, proj: &[f64], y: &mut [f64]) {
        let n = self.d.len();
        let c = self.count;
        for (col, yc) in y.chunks_mut(n).enumerate() {
            let pr = &proj[col * c..(col + 1) * c];
            for j in 0..c {
                let w: f64 = (0..c).map(|k| self.coupling[j * c + k] * pr[k]).sum();
                for (y, &p) in yc.iter_mut().zip(&self.p[j * n..(j + 1) * n]) {
                    *y += p * w;
                }
            }
        }
    }
    fn local(&self, x: &RealMatrix, y: &mut RealMatrix) {
        let n = self.d.len();
        for (xc, yc) in x.as_slice().chunks(n).zip(y.as_mut_slice().chunks_mut(n)) {
            for ((y, &x), &d) in yc.iter_mut().zip(xc).zip(&self.d) {
                *y = d * x;
            }
        }
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
        let c = self.count;
        self.local(x, y);
        let mut last = self.last.borrow_mut();
        if last.len() < c * x.cols() {
            last.resize(c * x.cols(), 0.);
        }
        for (col, xc) in x.as_slice().chunks(n).enumerate() {
            for j in 0..c {
                last[col * c + j] = self.p[j * n..(j + 1) * n]
                    .iter()
                    .zip(xc)
                    .map(|(a, b)| a * b)
                    .sum();
            }
        }
        self.add_nonlocal(&last[..c * x.cols()], y.as_mut_slice());
        Ok(())
    }
    fn projector_count(&self) -> usize {
        self.count
    }
    fn copy_last_projections(&self, m: usize, out: &mut [f64]) -> bool {
        let last = self.last.borrow();
        if out.len() != m * self.count || last.len() < out.len() {
            return false;
        }
        out.copy_from_slice(&last[..out.len()]);
        true
    }
    fn apply_block_with_projections(
        &self,
        x: &RealMatrix,
        projections: &[f64],
        y: &mut RealMatrix,
    ) -> Result<(), SolverError> {
        self.cached_calls.set(self.cached_calls.get() + 1);
        self.local(x, y);
        self.add_nonlocal(projections, y.as_mut_slice());
        Ok(())
    }
}

fn solve(op: &Operator, b6: bool, capacity_per_band: usize) -> RealEigenSolution {
    let n = op.dimension();
    let s = RealDavidson {
        bands: 12,
        tolerance: 1e-9,
        max_iterations: 400,
        max_subspace: 0,
        seed: 7,
    };
    let policy = BandTolerances {
        strict_bands: 9,
        buffer_tolerance: 1e-6,
    };
    let diagonal: Vec<f64> = op.d.clone();
    assert_eq!(diagonal.len(), n);
    s.solve_research(
        op,
        None,
        &diagonal,
        policy,
        false,
        true,
        capacity_per_band,
        &mut None,
        DavidsonPreconditioner::Default,
        ResearchOptions {
            b6,
            ..ResearchOptions::default()
        },
    )
    .unwrap()
}

#[test]
fn b6_certificate_matches_plain_apply() {
    for capacity_per_band in [2, 6] {
        let mut reference: Option<RealEigenSolution> = None;
        for threads in [1, 3] {
            let plain_op = Operator::new(240, 6, threads);
            let plain = solve(&plain_op, false, capacity_per_band);
            assert_eq!(plain_op.cached_calls.get(), 0);
            let op = Operator::new(240, 6, threads);
            let cached = solve(&op, true, capacity_per_band);
            // Every certificate took the cached path (the cache survived appends,
            // contractions and rebuilds).
            assert!(op.cached_calls.get() >= 1, "cached certificate not used");
            assert_eq!(op.cached_calls.get(), cached.certification_attempts);
            assert!(cached.residuals_certified && plain.residuals_certified);
            assert_eq!(cached.eigenvalues, plain.eigenvalues);
            assert_eq!(cached.eigenvectors.as_slice(), plain.eigenvectors.as_slice());
            assert_eq!(cached.iterations, plain.iterations);
            assert_eq!(cached.operator_applications, plain.operator_applications);
            assert_eq!(cached.certification_attempts, plain.certification_attempts);
            assert_eq!(cached.converged, plain.converged);
            for (a, b) in cached.residuals.iter().zip(&plain.residuals) {
                assert!((a - b).abs() <= 1e-13, "{a:e} vs {b:e}");
            }
            if capacity_per_band == 2 {
                assert!(cached.restarts > 0, "fixture must restart");
                assert!(cached.contracted_restarts > 0, "fixture must contract");
            }
            // Thread-count independence of the cached result.
            match &reference {
                None => reference = Some(cached),
                Some(r) => {
                    assert_eq!(r.eigenvalues, cached.eigenvalues);
                    assert_eq!(r.residuals, cached.residuals);
                }
            }
        }
    }
}
