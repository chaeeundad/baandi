//! Complex, owned-row reference block Davidson. Gamma keeps its real solver.
//! Uses the reference Ritz-change/empty-band rules and Hermitian overlap metric.
//! All ranks process the same k point; there are no k-parallel pools here.
use crate::ref_davidson::accounting::{Counted, clock, elapsed};
use crate::{
    SolverError,
    ref_davidson::{self, Iteration, Options, Record, ReducedPlacement},
};
use dft_mpi::Collective;
use dft_numerics::{
    linalg::{
        ComplexMatrix as Matrix, ComplexOverlapFactor, DEPENDENCE_TOLERANCE, GemmWorkspace,
        MatrixOp,
    },
    pool::WorkerPool,
};
use num_complex::Complex64 as C;
use std::sync::Arc;

pub trait DistributedComplexOperator {
    fn global_len(&self) -> usize;
    fn local_rows(&self) -> &[usize];
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        None
    }
    /// Implementations agree local failures before entering later internal collectives.
    fn apply_block<K: Collective>(
        &mut self,
        comm: &K,
        x: &Matrix,
        out: &mut Matrix,
    ) -> Result<(), SolverError>;
}
#[derive(Debug)]
pub struct Solution {
    pub eigenvalues: Vec<f64>,
    pub eigenvectors: Matrix,
    pub hx: Option<Matrix>,
    pub converged: Vec<bool>,
    pub ritz_changes_ha: Vec<Vec<f64>>,
    pub ritz_values_ha: Vec<Vec<[f64; 2]>>,
    pub band_thresholds_ha: Vec<f64>,
    pub record: Record,
}
fn error(e: impl std::fmt::Display) -> SolverError {
    SolverError::Application(e.to_string())
}
fn agree<K: Collective>(comm: &K, bad: bool) -> Result<(), SolverError> {
    let mut out = [0.];
    comm.tree_sum_f64(&[f64::from(bad)], &mut out)
        .map_err(error)?;
    if out[0] != 0. {
        return Err(error("complex reference Davidson invalid input on a rank"));
    }
    Ok(())
}
fn reduce_into<K: Collective>(
    comm: &K,
    a: &mut Matrix,
    local: &mut [f64],
    all: &mut [f64],
) -> Result<(), SolverError> {
    let len = 2 * a.as_slice().len();
    for (p, z) in local[..len]
        .as_chunks_mut::<2>()
        .0
        .iter_mut()
        .zip(a.as_slice())
    {
        *p = [z.re, z.im];
    }
    comm.tree_sum_f64(&local[..len], &mut all[..len])
        .map_err(error)?;
    for (z, p) in a
        .as_mut_slice()
        .iter_mut()
        .zip(all[..len].as_chunks::<2>().0)
    {
        *z = C::new(p[0], p[1]);
    }
    Ok(())
}
#[cfg(test)]
fn reduce<K: Collective>(comm: &K, a: &mut Matrix) -> Result<(), SolverError> {
    let len = 2 * a.as_slice().len();
    reduce_into(comm, a, &mut vec![0.; len], &mut vec![0.; len])
}
struct Products {
    gemm: GemmWorkspace,
    local: Vec<f64>,
    all: Vec<f64>,
}
impl Products {
    fn gram<K: Collective>(
        &mut self,
        comm: &K,
        a: &Matrix,
        b: &Matrix,
        out: &mut Matrix,
        r: &mut Iteration,
    ) -> Result<(), SolverError> {
        let t = clock(r.profile);
        out.reshape_within_capacity(a.cols(), b.cols())?;
        let product = self
            .gemm
            .multiply_into(a, MatrixOp::Adjoint, b, MatrixOp::Identity, out);
        agree(comm, product.is_err())?;
        product?;
        r.gemm("complex_gram", a.cols(), b.cols(), a.rows());
        if let Some(g) = r.gemms.last_mut() {
            g.flops *= 4;
        }
        reduce_into(comm, out, &mut self.local, &mut self.all)?;
        r.seconds[1] += elapsed(t);
        Ok(())
    }
}
fn hermitize(a: &mut Matrix) {
    for j in 0..a.cols() {
        a[(j, j)].im = 0.;
        for i in 0..j {
            let z = (a[(i, j)] + a[(j, i)].conj()) * 0.5;
            a[(i, j)] = z;
            a[(j, i)] = z.conj();
        }
    }
}
fn rotate<K: Collective>(
    comm: &K,
    g: &mut GemmWorkspace,
    a: &Matrix,
    c: &Matrix,
    out: &mut Matrix,
    r: &mut Iteration,
) -> Result<(), SolverError> {
    out.resize_columns(c.cols())?;
    let product = g.multiply_into(a, MatrixOp::Identity, c, MatrixOp::Identity, out);
    agree(comm, product.is_err())?;
    product?;
    r.gemm("complex_rotate", a.rows(), c.cols(), a.cols());
    if let Some(g) = r.gemms.last_mut() {
        g.flops *= 4;
    }
    Ok(())
}
fn identity(a: &mut Matrix, n: usize) -> Result<(), SolverError> {
    a.reshape_within_capacity(n, n)?;
    a.as_mut_slice().fill(C::default());
    for i in 0..n {
        a[(i, i)] = C::new(1., 0.);
    }
    Ok(())
}
fn factor(
    f: &mut ComplexOverlapFactor,
    s: &Matrix,
    r: &mut Iteration,
) -> Result<bool, SolverError> {
    let t = clock(r.profile);
    let ok = f.factor(s, crate::real_davidson::PIVOT_RATIO_MIN)?;
    r.seconds[2] += elapsed(t);
    Ok(ok)
}
#[allow(clippy::too_many_arguments)]
fn eigen<K: Collective>(
    comm: &K,
    f: &mut ComplexOverlapFactor,
    h: &Matrix,
    o: Options,
    r: &mut Iteration,
    values: &mut [f64],
    c: &mut Matrix,
    packed: &mut [f64],
) -> Result<(), SolverError> {
    let t = clock(o.profile);
    let root = o.placement == ReducedPlacement::Root;
    c.reshape_within_capacity(h.rows(), o.bands)?;
    let result = if !root || comm.rank() == 0 {
        r.backend_solves += 1;
        f.solve_lowest_into(h, values, c)
    } else {
        Ok(())
    };
    agree(comm, result.is_err())?;
    result?;
    if root {
        comm.bcast_f64(values, 0).map_err(error)?;
        let len = 2 * c.as_slice().len();
        for (p, z) in packed[..len]
            .as_chunks_mut::<2>()
            .0
            .iter_mut()
            .zip(c.as_slice())
        {
            *p = [z.re, z.im];
        }
        comm.bcast_f64(&mut packed[..len], 0).map_err(error)?;
        for (z, p) in c
            .as_mut_slice()
            .iter_mut()
            .zip(packed[..len].as_chunks::<2>().0)
        {
            *z = C::new(p[0], p[1]);
        }
    }
    r.reduced_solves += 1;
    r.reduced_order = h.rows();
    r.seconds[3] += elapsed(t);
    Ok(())
}
fn matrix(rows: usize, cols: usize, capacity: usize) -> Result<Matrix, SolverError> {
    let mut a = Matrix::zeros(rows, cols)?;
    a.reserve_columns(capacity)?;
    Ok(a)
}

pub fn solve<K: Collective, O: DistributedComplexOperator>(
    comm: &K,
    operator: &mut O,
    start: &Matrix,
    diagonal: &[f64],
    options: Options,
) -> Result<Solution, SolverError> {
    solve_with_initial_hx(comm, operator, start, diagonal, options, None)
}
pub fn solve_with_initial_hx<K: Collective, O: DistributedComplexOperator>(
    comm: &K,
    operator: &mut O,
    start: &Matrix,
    diagonal: &[f64],
    o: Options,
    initial_hx: Option<&Matrix>,
) -> Result<Solution, SolverError> {
    let total = clock(o.profile);
    let comm = &Counted::new(comm, o.profile);
    let n = o.bands;
    let dim = operator.global_len();
    let rows = operator.local_rows().len();
    let cap = ref_davidson::workspace::capacity(dim, n, o.capacity);
    agree(
        comm,
        n == 0
            || n > dim
            || cap < n.saturating_mul(2).min(dim)
            || o.occupied > n
            || o.max_iterations == 0
            || !o.threshold.is_finite()
            || o.threshold <= 0.
            || !o.band_threshold(n.saturating_sub(1)).is_finite()
            || start.rows() != rows
            || start.cols() != n
            || diagonal.len() != rows
            || start.validate_finite().is_err()
            || diagonal.iter().any(|x| !x.is_finite())
            || initial_hx
                .is_some_and(|x| x.rows() != rows || x.cols() != n || x.validate_finite().is_err()),
    )?;
    let mut has = [0.];
    comm.tree_sum_f64(&[f64::from(initial_hx.is_some())], &mut has)
        .map_err(error)?;
    if has[0] != 0. && has[0] != comm.size() as f64 {
        return Err(error("complex initial HX differs across ranks"));
    }
    let pool = match operator.worker_pool() {
        Some(p) => p,
        None => WorkerPool::new(1).map_err(error)?,
    };
    let mut products = Products {
        gemm: GemmWorkspace::with_pool(rows.max(cap), cap, pool.workers_created().max(1), pool)?,
        local: vec![0.; 2 * cap * cap],
        all: vec![0.; 2 * cap * cap],
    };
    let mut f = ComplexOverlapFactor::with_capacity(cap)?;
    let mut trial = ComplexOverlapFactor::with_capacity(cap)?;
    let mut s = Matrix::zeros(cap, cap)?;
    let mut h = s.clone();
    let mut c = s.clone();
    let mut selected = s.clone();
    let mut transform = s.clone();
    let mut proj = s.clone();
    let mut q = matrix(rows, cap, cap)?;
    let mut hq = matrix(rows, cap, cap)?;
    let mut values = vec![0.; n];
    let mut old = vec![0.; n];
    let mut v = matrix(rows, n, cap)?;
    v.as_mut_slice().copy_from_slice(start.as_slice());
    let mut hv = matrix(rows, n, cap)?;
    let mut x = matrix(rows, n, cap)?;
    let mut hx = matrix(rows, n, cap)?;
    let mut e = matrix(rows, n, n)?;
    let mut he = matrix(rows, n, n)?;
    let mut record = Record::new(o.max_iterations, cap, o.profile, o.placement);
    let t = clock(o.profile);
    if let Some(h) = initial_hx {
        hv.as_mut_slice().copy_from_slice(h.as_slice());
    } else {
        let applied = operator.apply_block(comm, &v, &mut hv);
        agree(comm, applied.is_err() || hv.validate_finite().is_err())?;
        applied?;
        record.setup.h_calls += 1;
        record.setup.h_columns += n;
    }
    record.setup.seconds[0] += elapsed(t);
    products.gram(comm, &v, &v, &mut s, &mut record.setup)?;
    hermitize(&mut s);
    if !factor(&mut f, &s, &mut record.setup)? {
        return Err(error(
            "complex reference Davidson start is numerically dependent",
        ));
    }
    products.gram(comm, &v, &hv, &mut h, &mut record.setup)?;
    hermitize(&mut h);
    eigen(
        comm,
        &mut f,
        &h,
        o,
        &mut record.setup,
        &mut values,
        &mut c,
        &mut products.local,
    )?;
    record.setup.communication = comm.snapshot();
    let mut active: Vec<_> = (0..n).collect();
    let mut converged = vec![false; n];
    let mut changes = if o.profile {
        vec![vec![0.; n]; o.max_iterations]
    } else {
        Vec::new()
    };
    let mut history = if o.profile {
        vec![vec![[0.; 2]; n]; o.max_iterations]
    } else {
        Vec::new()
    };
    let bounds = diagonal
        .iter()
        .fold((f64::INFINITY, f64::NEG_INFINITY), |(lo, hi), &v| {
            (lo.min(v), hi.max(v))
        });
    let mut contracted = false;
    let mut contractions = 0usize;
    let mut reset_disabled = false;
    let mut norms = vec![0.; n];
    let mut sums = vec![0.; n];
    let mut iterations = 0;
    for (step, r) in record.iterations.iter_mut().enumerate() {
        let before = comm.snapshot();
        let m = v.cols();
        r.basis_before = m;
        r.active_bands.extend_from_slice(&active);
        let t = clock(o.profile);
        selected.reshape_within_capacity(m, active.len())?;
        for (j, &b) in active.iter().enumerate() {
            selected.column_mut(j).copy_from_slice(c.column(b));
        }
        if contracted {
            x.resize_columns(active.len())?;
            hx.resize_columns(active.len())?;
            for (j, &band) in active.iter().enumerate() {
                x.column_mut(j).copy_from_slice(v.column(band));
                hx.column_mut(j).copy_from_slice(hv.column(band));
            }
        } else {
            rotate(comm, &mut products.gemm, &v, &selected, &mut x, r)?;
            rotate(comm, &mut products.gemm, &hv, &selected, &mut hx, r)?;
        }
        contracted = false;
        e.resize_columns(active.len())?;
        norms[..active.len()].fill(0.);
        for (j, &b) in active.iter().enumerate() {
            for (i, &d) in diagonal.iter().enumerate() {
                let denominator = if ref_davidson::stream::sqrt_safe(bounds, values[b]) {
                    ref_davidson::stream::sqrt_denominator(d, values[b])
                } else {
                    ref_davidson::preconditioner_denominator(d, values[b])
                };
                let z = (hx[(i, j)] - values[b] * x[(i, j)]) / denominator;
                e[(i, j)] = z;
                norms[j] += z.norm_sqr();
            }
        }
        comm.tree_sum_f64(&norms[..active.len()], &mut sums[..active.len()])
            .map_err(error)?;
        let mut kept = 0;
        for (j, &norm2) in sums[..active.len()].iter().enumerate() {
            let norm = norm2.sqrt();
            if !norm.is_finite() {
                return Err(error("nonfinite complex correction"));
            }
            if norm <= DEPENDENCE_TOLERANCE {
                r.dropped_corrections += 1;
                continue;
            }
            for i in 0..rows {
                e[(i, kept)] = e[(i, j)] / norm;
            }
            kept += 1;
        }
        e.resize_columns(kept)?;
        r.seconds[4] += elapsed(t);
        // Try the complete normalized correction block, then repair only if
        // its metric is dependent or the full physical space has been reached.
        let mut accepted = false;
        if m + kept <= cap && kept > 0 {
            v.append_columns(&e)?;
            products.gram(comm, &v, &v, &mut s, r)?;
            hermitize(&mut s);
            accepted = factor(&mut trial, &s, r)?;
        }
        if !accepted {
            let t = clock(o.profile);
            let child_seconds = r.seconds[1] + r.seconds[2];
            r.repairs += 1;
            v.resize_columns(m)?;
            // Whiten the OLD basis and HV without reapplying H. Preserve the
            // physical Ritz block through a freshly measured coordinate map.
            rotate(comm, &mut products.gemm, &v, &c, &mut x, r)?;
            f.basis_transform_into(&mut transform)?;
            rotate(comm, &mut products.gemm, &v, &transform, &mut q, r)?;
            rotate(comm, &mut products.gemm, &hv, &transform, &mut hq, r)?;
            products.gram(comm, &q, &x, &mut c, r)?;
            std::mem::swap(&mut v, &mut q);
            std::mem::swap(&mut hv, &mut hq);
            for _ in 0..2 {
                products.gram(comm, &v, &e, &mut proj, r)?;
                rotate(comm, &mut products.gemm, &v, &proj, &mut x, r)?;
                for (a, b) in e.as_mut_slice().iter_mut().zip(x.as_slice()) {
                    *a -= b;
                }
            }
            let input = e.cols();
            kept = 0;
            for j in 0..input {
                for _ in 0..2 {
                    for a in 0..kept {
                        let z: C = (0..rows).map(|i| e[(i, a)].conj() * e[(i, j)]).sum();
                        let mut sum = [0.; 2];
                        comm.tree_sum_f64(&[z.re, z.im], &mut sum).map_err(error)?;
                        let z = C::new(sum[0], sum[1]);
                        for i in 0..rows {
                            let q = e[(i, a)];
                            e[(i, j)] -= z * q;
                        }
                    }
                }
                let z: f64 = e.column(j).iter().map(|z| z.norm_sqr()).sum();
                let mut sum = [0.];
                comm.tree_sum_f64(&[z], &mut sum).map_err(error)?;
                let norm = sum[0].sqrt();
                if norm <= DEPENDENCE_TOLERANCE || m + kept == cap {
                    r.dropped_corrections += 1;
                    continue;
                }
                for i in 0..rows {
                    e[(i, kept)] = e[(i, j)] / norm;
                }
                kept += 1;
            }
            e.resize_columns(kept)?;
            v.append_columns(&e)?;
            products.gram(comm, &v, &v, &mut s, r)?;
            hermitize(&mut s);
            accepted = factor(&mut trial, &s, r)?;
            if !accepted {
                return Err(error("complex overlap remains dependent after repair"));
            }
            r.seconds[5] += (elapsed(t) - (r.seconds[1] + r.seconds[2] - child_seconds)).max(0.);
        }
        std::mem::swap(&mut f, &mut trial);
        old.copy_from_slice(&values);
        if kept > 0 {
            let t = clock(o.profile);
            he.resize_columns(kept)?;
            let applied = operator.apply_block(comm, &e, &mut he);
            agree(comm, applied.is_err() || he.validate_finite().is_err())?;
            applied?;
            r.h_calls += 1;
            r.h_columns += kept;
            r.seconds[0] += elapsed(t);
            hv.append_columns(&he)?;
            products.gram(comm, &v, &hv, &mut h, r)?;
            hermitize(&mut h);
            eigen(
                comm,
                &mut f,
                &h,
                o,
                r,
                &mut values,
                &mut c,
                &mut products.local,
            )?;
        } else {
            r.reduced_order = m;
        }
        active.clear();
        for i in 0..n {
            let delta = (values[i] - old[i]).abs();
            converged[i] = delta < o.band_threshold(i);
            if !converged[i] {
                active.push(i);
            }
            if o.profile {
                changes[step][i] = delta;
                history[step][i] = [old[i], values[i]];
            }
        }
        r.unconverged = active.len();
        iterations = step + 1;
        if !active.is_empty() && iterations < o.max_iterations && v.cols() + active.len() > cap {
            let t = clock(o.profile);
            let child_seconds = r.seconds[1] + r.seconds[2];
            rotate(comm, &mut products.gemm, &v, &c, &mut x, r)?;
            rotate(comm, &mut products.gemm, &hv, &c, &mut hx, r)?;
            {
                std::mem::swap(&mut v, &mut x);
                std::mem::swap(&mut hv, &mut hx);
            }
            contractions += 1;
            contracted = true;
            identity(&mut c, n)?;
            let mut ideal = false;
            {
                r.resets.attempts += 1;
                if reset_disabled {
                    r.resets.disabled_refreshes += 1;
                } else if contractions.is_multiple_of(ref_davidson::ideal_reset::REFRESH_INTERVAL) {
                    r.resets.periodic_refreshes += 1;
                } else {
                    for (j, norm) in norms.iter_mut().enumerate() {
                        *norm = v.column(j).iter().map(|z| z.norm_sqr()).sum();
                    }
                    comm.tree_sum_f64(&norms, &mut sums).map_err(error)?;
                    ideal = sums.iter().all(|v| {
                        v.is_finite()
                            && (v - 1.).abs() <= ref_davidson::ideal_reset::ORTHOGONALITY_LIMIT
                    });
                    if ideal {
                        r.resets.ideal += 1;
                    } else {
                        reset_disabled = true;
                        r.resets.guard_fallbacks += 1;
                    }
                }
            }
            if ideal {
                identity(&mut s, n)?;
                identity(&mut h, n)?;
                for j in 0..n {
                    h[(j, j)] = C::new(values[j], 0.);
                }
            } else {
                products.gram(comm, &v, &v, &mut s, r)?;
                hermitize(&mut s);
                products.gram(comm, &v, &hv, &mut h, r)?;
                hermitize(&mut h);
                if !reset_disabled {
                    let mut defect = 0_f64;
                    for j in 0..n {
                        for i in 0..n {
                            defect = defect.max((s[(i, j)] - C::new(f64::from(i == j), 0.)).norm());
                        }
                    }
                    r.resets.max_measured_defect = defect;
                    if !defect.is_finite()
                        || defect > ref_davidson::ideal_reset::ORTHOGONALITY_LIMIT
                    {
                        reset_disabled = true;
                        r.resets.guard_fallbacks += 1;
                    }
                }
            }
            if !factor(&mut f, &s, r)? {
                return Err(error("complex restart metric is dependent"));
            }
            r.restarts += 1;
            r.seconds[6] += (elapsed(t) - (r.seconds[1] + r.seconds[2] - child_seconds)).max(0.);
        }
        r.communication = comm.snapshot().difference(before);
        if active.is_empty() || iterations == o.max_iterations {
            break;
        }
    }
    record.iterations.truncate(iterations);
    changes.truncate(iterations);
    history.truncate(iterations);
    let t = clock(o.profile);
    rotate(comm, &mut products.gemm, &v, &c, &mut x, &mut record.finish)?;
    if o.return_hx {
        rotate(
            comm,
            &mut products.gemm,
            &hv,
            &c,
            &mut hx,
            &mut record.finish,
        )?;
    }
    record.finish.seconds[7] += elapsed(t);
    record.communication = comm.snapshot();
    record.total_seconds = elapsed(total);
    record.summarize();
    Ok(Solution {
        eigenvalues: values,
        eigenvectors: x,
        hx: o.return_hx.then_some(hx),
        converged,
        ritz_changes_ha: changes,
        ritz_values_ha: history,
        band_thresholds_ha: (0..n).map(|i| o.band_threshold(i)).collect(),
        record,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    struct Dense {
        h: Matrix,
        rows: Vec<usize>,
    }
    impl DistributedComplexOperator for Dense {
        fn global_len(&self) -> usize {
            self.h.rows()
        }
        fn local_rows(&self) -> &[usize] {
            &self.rows
        }
        fn apply_block<K: Collective>(
            &mut self,
            comm: &K,
            x: &Matrix,
            out: &mut Matrix,
        ) -> Result<(), SolverError> {
            let mut full = Matrix::zeros(self.h.rows(), x.cols())?;
            for j in 0..x.cols() {
                for (i, &g) in self.rows.iter().enumerate() {
                    full[(g, j)] = x[(i, j)];
                }
            }
            reduce(comm, &mut full)?;
            for j in 0..x.cols() {
                for (i, &g) in self.rows.iter().enumerate() {
                    out[(i, j)] = (0..self.h.rows())
                        .map(|a| self.h[(g, a)] * full[(a, j)])
                        .sum();
                }
            }
            Ok(())
        }
    }
    #[test]
    fn complex_reference_matches_hermitian_oracle_and_supplied_hx() {
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let dim = 11;
                let n = 3;
                let mut h = Matrix::zeros(dim, dim).unwrap();
                for i in 0..dim {
                    h[(i, i)] = C::new(i as f64 * 0.3 - 0.7, 0.);
                    if i > 0 {
                        h[(i - 1, i)] = C::new(0.09, 0.07);
                        h[(i, i - 1)] = h[(i - 1, i)].conj();
                    }
                }
                let (exact, _) = dft_numerics::linalg::hermitian_eigh(&h).unwrap();
                let rows: Vec<_> = (0..dim).filter(|i| i % ranks == comm.rank()).collect();
                let mut start = Matrix::zeros(rows.len(), n).unwrap();
                for j in 0..n {
                    for (i, &g) in rows.iter().enumerate() {
                        start[(i, j)] = C::new(f64::from(g == j), 0.)
                            + C::new(
                                ((g + 1) * (j + 2)) as f64 * 0.001,
                                ((g + 2) * (j + 1)) as f64 * 0.002,
                            );
                    }
                }
                let diagonal: Vec<_> = rows.iter().map(|&i| h[(i, i)].re).collect();
                let mut op = Dense { h, rows };
                for (placement, flags) in [
                    (ReducedPlacement::Root, "off"),
                    (ReducedPlacement::Replicated, "on"),
                ] {
                    let config = dft_numerics::switches::RunConfig::from_pairs(&[]).unwrap();
                    let _scope = config.enter();
                    let mut o = Options::new(n, 1e-13, n);
                    o.max_iterations = 100;
                    o.return_hx = true;
                    o.profile = true;
                    o.placement = placement;
                    let a = solve(comm, &mut op, &start, &diagonal, o).unwrap();
                    assert!(a.converged.iter().all(|b| *b));
                    if flags == "on" {
                        assert!(a.record.resets.attempts > 0);
                    }
                    for (e, want) in a.eigenvalues.iter().zip(&exact) {
                        assert!((e - want).abs() < 1e-11, "{e} {want}");
                    }
                    let mut hs = Matrix::zeros(start.rows(), n).unwrap();
                    op.apply_block(comm, &start, &mut hs).unwrap();
                    let b = solve_with_initial_hx(comm, &mut op, &start, &diagonal, o, Some(&hs))
                        .unwrap();
                    assert_eq!(a.eigenvalues, b.eigenvalues);
                    assert_eq!(a.eigenvectors, b.eigenvectors);
                    assert_eq!(a.record.h_columns, b.record.h_columns + n);
                }
                let mut bad = start.clone();
                bad.column_mut(1).copy_from_slice(start.column(0));
                assert!(solve(comm, &mut op, &bad, &diagonal, Options::new(n, 1e-9, n)).is_err());
            });
        }
    }
    #[test]
    fn complex_reference_agrees_nonfinite_operator_output() {
        struct Fault {
            rows: Vec<usize>,
        }
        impl DistributedComplexOperator for Fault {
            fn global_len(&self) -> usize {
                1
            }
            fn local_rows(&self) -> &[usize] {
                &self.rows
            }
            fn apply_block<K: Collective>(
                &mut self,
                _comm: &K,
                _x: &Matrix,
                out: &mut Matrix,
            ) -> Result<(), SolverError> {
                out.as_mut_slice().fill(C::new(f64::NAN, 0.));
                Ok(())
            }
        }
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let rows = if comm.rank() == 0 {
                    vec![0]
                } else {
                    Vec::new()
                };
                let mut start = Matrix::zeros(rows.len(), 1).unwrap();
                start.as_mut_slice().fill(C::new(1., 0.));
                let diagonal = vec![1.; rows.len()];
                let mut op = Fault { rows };
                assert!(solve(comm, &mut op, &start, &diagonal, Options::new(1, 1e-8, 1)).is_err());
            });
        }
    }
}
