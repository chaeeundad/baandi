use super::accounting::{clock, elapsed};
use super::{Iteration, Options, ReducedPlacement, error, preconditioner_denominator};
use crate::{SolverError, distributed_davidson::DistributedRealOperator, real_davidson::lean};
use dft_mpi::Collective;
use dft_numerics::{
    GemmSite, RealGemmWorkspace, RealMatrix, RealOp, pool::WorkerPool, real::RealEighWorkspace,
};

pub(super) struct Workspace {
    rows: usize,
    n: usize,
    cap: usize,
    pub v: RealMatrix,
    pub hv: RealMatrix,
    pub x: RealMatrix,
    pub hx: RealMatrix,
    pub e: RealMatrix,
    he: RealMatrix,
    pub values: Vec<f64>,
    pub old_values: Vec<f64>,
    lean: lean::LeanScratch,
    hr: RealMatrix,
    h: RealMatrix,
    y: RealMatrix,
    coeff: RealMatrix,
    border: RealMatrix,
    reduced: Vec<f64>,
    norms: Vec<f64>,
    sums: Vec<f64>,
    gemm: RealGemmWorkspace,
    evd: RealEighWorkspace,
}
impl Workspace {
    pub fn new<O: DistributedRealOperator>(
        op: &O,
        n: usize,
        cap: usize,
    ) -> Result<Self, SolverError> {
        let rows = op.local_rows().len();
        let pool = op
            .worker_pool()
            .map(Ok)
            .unwrap_or_else(|| WorkerPool::new(op.inner_threads()))
            .map_err(|e| SolverError::Application(e.to_string()))?;
        Ok(Self {
            rows,
            n,
            cap,
            v: RealMatrix::zeros(rows, cap)?,
            hv: RealMatrix::zeros(rows, cap)?,
            x: RealMatrix::zeros(rows, cap)?,
            hx: RealMatrix::zeros(rows, cap)?,
            e: RealMatrix::zeros(rows, n)?,
            he: RealMatrix::zeros(rows, n)?,
            values: vec![0.0; cap],
            old_values: vec![0.0; n],
            lean: lean::LeanScratch::new(cap, n)?,
            hr: RealMatrix::zeros(cap, cap)?,
            h: RealMatrix::zeros(cap, cap)?,
            y: RealMatrix::zeros(cap, cap)?,
            coeff: RealMatrix::zeros(cap, cap)?,
            border: RealMatrix::zeros(cap, cap)?,
            reduced: vec![0.0; cap * cap],
            norms: vec![0.0; cap + 1],
            sums: vec![0.0; cap + 1],
            gemm: RealGemmWorkspace::with_pool(rows.max(cap), cap, pool)?,
            evd: RealEighWorkspace::new(cap)?,
        })
    }
    pub fn apply<C: Collective, O: DistributedRealOperator>(
        &mut self,
        c: &C,
        op: &mut O,
        x: &RealMatrix,
        _initial: bool,
        r: &mut Iteration,
    ) -> Result<(), SolverError> {
        r.h_calls += 1;
        r.h_columns += x.cols();
        let result = op.apply_block(c, x, &mut self.hv);
        agree_apply(c, result, &self.hv)
    }
    pub fn apply_corrections<C: Collective, O: DistributedRealOperator>(
        &mut self,
        c: &C,
        op: &mut O,
        m: usize,
        r: &mut Iteration,
    ) -> Result<(), SolverError> {
        let k = self.e.cols();
        self.he.reshape_within_storage(self.rows, k)?;
        r.h_calls += 1;
        r.h_columns += k;
        let result = op.apply_block(c, &self.e, &mut self.he);
        agree_apply(c, result, &self.he)?;
        self.hv.reshape_within_storage(self.rows, m + k)?;
        self.hv.as_mut_slice()[m * self.rows..(m + k) * self.rows]
            .copy_from_slice(self.he.as_slice());
        Ok(())
    }
    /// Only the new columns; mirror their upper triangles. The leading square
    /// is never recomputed. S is formed before H so failed columns cost no H.
    pub fn border<C: Collective>(
        &mut self,
        c: &C,
        from: usize,
        to: usize,
        is_h: bool,
        r: &mut Iteration,
        profile: bool,
    ) -> Result<(), SolverError> {
        let t = clock(profile);
        let k = to - from;
        self.border.reshape_within_storage(to, k)?;
        let b = if is_h { &self.hv } else { &self.v };
        r.gemm(if is_h { "border_h" } else { "border_s" }, to, k, self.rows);
        self.gemm.multiply_slices_into(
            GemmSite::Other,
            self.v.as_slice(),
            self.rows,
            to,
            self.rows,
            RealOp::Transpose,
            &b.as_slice()[from * self.rows..],
            k,
            self.rows,
            self.border.as_mut(),
        )?;
        c.tree_sum_f64(self.border.as_slice(), &mut self.reduced[..to * k])
            .map_err(error)?;
        let target = if is_h { &mut self.hr } else { &mut self.lean.s };
        for j in from..to {
            for i in 0..=j {
                let value = self.reduced[i + (j - from) * to];
                target[(i, j)] = value;
                target[(j, i)] = value;
            }
        }
        r.seconds[1] += elapsed(t);
        Ok(())
    }
    /// Reuse the engine's guarded bordered Cholesky adapter; no new numerical
    /// factorization is hand-written here. Record its actual panel GEMM shapes
    /// (64 is the adapter's existing FACTOR_PANEL, not a solver tuning knob).
    pub fn factor(
        &mut self,
        from: usize,
        to: usize,
        r: &mut Iteration,
        profile: bool,
    ) -> Result<bool, SolverError> {
        let t = clock(profile);
        // Invoke one panel at a time so a rejected pivot never charges GEMMs
        // from later panels. This is the same factor arithmetic/order.
        let mut p = from;
        let mut ok = true;
        while p < to {
            let width = 64.min(to - p);
            if p > 0 {
                r.gemm("factor_w", p, width, p);
                r.gemm("factor_gram", width, width, p);
            }
            let (accepted, ratio) =
                lean::extend_factor(&mut self.gemm, &mut self.lean, p, p + width)?;
            if p > 0 && (accepted || ratio >= lean::PIVOT_RATIO_MIN) {
                r.gemm("factor_u", width, p, p);
                r.gemm("factor_inverse", width, p, width);
            }
            if !accepted {
                ok = false;
                break;
            }
            p += width;
        }
        r.seconds[2] += elapsed(t);
        Ok(ok)
    }
    pub fn eigen<C: Collective>(
        &mut self,
        c: &C,
        from: usize,
        m: usize,
        r: &mut Iteration,
        o: &Options,
    ) -> Result<(), SolverError> {
        let t = clock(o.profile);
        self.h.reshape_within_storage(m, m)?;
        for j in 0..m {
            for i in 0..m {
                self.h[(i, j)] = self.hr[(i, j)];
            }
        }
        // The small transformed border is replicated; only the library EVD is
        // placed on root. This also leaves the repair transform on every rank.
        r.gemm("reduced_transform", m, m - from, m);
        r.gemm("reduced_transform", m, m - from, m);
        lean::reduced_border(&mut self.gemm, &self.h, &mut self.lean, from, m)?;
        self.y.reshape_within_storage(m, m)?;
        let root = o.placement == ReducedPlacement::Root;
        let mut success = [1.0];
        if !root || c.rank() == 0 {
            r.backend_solves += 1;
            if self
                .evd
                .solve_into(&self.lean.a, &mut self.values[..m], &mut self.y)
                .is_err()
            {
                success[0] = 0.0;
            }
        }
        if root {
            c.bcast_f64(&mut success, 0).map_err(error)?;
        } else {
            let local = success;
            c.tree_sum_f64(&local, &mut success).map_err(error)?;
            success[0] = f64::from(success[0] == c.size() as f64);
        }
        if success[0] == 0.0 {
            return Err(SolverError::Application(
                "reference reduced eigensolve failed".into(),
            ));
        }
        if root {
            c.bcast_f64(&mut self.values[..self.n], 0).map_err(error)?;
            c.bcast_f64(&mut self.y.as_mut_slice()[..m * self.n], 0)
                .map_err(error)?;
        }
        r.gemm("ritz_coefficients", m, self.n, m);
        lean::ritz_coefficients(&mut self.gemm, &mut self.lean, &self.y, m, self.n, None)?;
        r.reduced_order = m;
        r.reduced_solves += 1;
        r.seconds[3] += elapsed(t);
        Ok(())
    }
    pub fn corrections<C: Collective>(
        &mut self,
        c: &C,
        hamiltonian_diagonal: &[f64],
        active: &[usize],
        m: usize,
        r: &mut Iteration,
        profile: bool,
    ) -> Result<(), SolverError> {
        let t = clock(profile);
        let k = active.len();
        self.coeff.reshape_within_storage(m, k)?;
        for (j, &band) in active.iter().enumerate() {
            self.coeff.as_mut_slice()[j * m..(j + 1) * m]
                .copy_from_slice(&self.lean.c.as_slice()[band * m..(band + 1) * m]);
        }
        rotate(
            &mut self.gemm,
            &self.v,
            &self.coeff,
            &mut self.x,
            m,
            k,
            r,
            "residual_x",
        )?;
        rotate(
            &mut self.gemm,
            &self.hv,
            &self.coeff,
            &mut self.hx,
            m,
            k,
            r,
            "residual_hx",
        )?;
        self.e.reshape_within_storage(self.rows, k)?;
        self.norms[..k].fill(0.0);
        for (j, &band) in active.iter().enumerate() {
            for (i, &diagonal) in hamiltonian_diagonal.iter().enumerate() {
                let z = (self.hx[(i, j)] - self.values[band] * self.x[(i, j)])
                    / preconditioner_denominator(diagonal, self.values[band]);
                self.e[(i, j)] = z;
                self.norms[j] += z * z;
            }
        }
        c.tree_sum_f64(&self.norms[..k], &mut self.sums[..k])
            .map_err(error)?;
        let mut kept = 0;
        for j in 0..k {
            let norm = self.sums[j].sqrt();
            if !norm.is_finite() {
                return Err(SolverError::Application(
                    "nonfinite reference correction".into(),
                ));
            }
            // Only exact/roundoff-zero directions are omitted. Scale is the
            // dimensionless preconditioned residual, not the requested ethr.
            if norm <= 32.0 * f64::EPSILON {
                r.dropped_corrections += 1;
                continue;
            }
            for i in 0..self.rows {
                self.e[(i, kept)] = self.e[(i, j)] / norm;
            }
            kept += 1;
        }
        self.e.reshape_within_storage(self.rows, kept)?;
        r.seconds[4] += elapsed(t);
        Ok(())
    }
    pub fn append_e(&mut self, m: usize) -> Result<(), SolverError> {
        self.v
            .reshape_within_storage(self.rows, m + self.e.cols())?;
        self.v.as_mut_slice()[m * self.rows..].copy_from_slice(self.e.as_slice());
        Ok(())
    }
    /// Rank repair is exceptional: whiten the OLD full basis with its valid
    /// inverse factor, preserving its span and HV, then two-pass projection of
    /// the correction block. No H reapplication or blind truncation of roots.
    pub fn repair<C: Collective>(
        &mut self,
        c: &C,
        m: usize,
        r: &mut Iteration,
        profile: bool,
    ) -> Result<usize, SolverError> {
        let t = clock(profile);
        r.repairs += 1;
        let input = self.e.cols();
        // Coordinates in Q=V T^T are T S C. Do not assume the old S is I.
        self.lean.p.reshape_within_storage(m, self.n)?;
        self.y.reshape_within_storage(m, self.n)?;
        r.gemm("repair_coefficients", m, self.n, m);
        r.gemm("repair_coefficients", m, self.n, m);
        self.gemm.multiply_slices_into(
            GemmSite::Other,
            self.lean.s.as_slice(),
            m,
            m,
            self.cap,
            RealOp::Identity,
            self.lean.c.as_slice(),
            self.n,
            m,
            self.lean.p.as_mut(),
        )?;
        self.gemm.multiply_slices_into(
            GemmSite::Other,
            self.lean.t.as_slice(),
            m,
            m,
            self.cap,
            RealOp::Identity,
            self.lean.p.as_slice(),
            self.n,
            m,
            self.y.as_mut(),
        )?;
        self.coeff.reshape_within_storage(m, m)?;
        for j in 0..m {
            for i in 0..m {
                self.coeff[(i, j)] = self.lean.t[(j, i)];
            }
        }
        rotate(
            &mut self.gemm,
            &self.v,
            &self.coeff,
            &mut self.x,
            m,
            m,
            r,
            "repair_basis",
        )?;
        rotate(
            &mut self.gemm,
            &self.hv,
            &self.coeff,
            &mut self.hx,
            m,
            m,
            r,
            "repair_hbasis",
        )?;
        self.v.reshape_within_storage(self.rows, m)?;
        self.hv.reshape_within_storage(self.rows, m)?;
        self.v.as_mut_slice().copy_from_slice(self.x.as_slice());
        self.hv.as_mut_slice().copy_from_slice(self.hx.as_slice());
        // Preserve the physical Ritz block through the whitening.
        self.lean.c.reshape_within_storage(m, self.n)?;
        self.lean
            .c
            .as_mut_slice()
            .copy_from_slice(&self.y.as_slice()[..m * self.n]);
        // Measured fresh borders prevent cancellation errors in repeated
        // contractions/whitenings from becoming a stale diagonal H or unit S.
        r.seconds[5] += elapsed(t);
        self.refresh(c, m, r, profile)?;
        let t = clock(profile);
        // Full-block projection retains the wide product shape. Repeating it
        // costs four GEMMs, not a per-band chain against the old basis.
        self.coeff.reshape_within_storage(m, input)?;
        self.x.reshape_within_storage(self.rows, input)?;
        for _ in 0..2 {
            r.gemm("repair_projection", m, input, self.rows);
            self.gemm.multiply_slices_into(
                GemmSite::Other,
                self.v.as_slice(),
                self.rows,
                m,
                self.rows,
                RealOp::Transpose,
                self.e.as_slice(),
                input,
                self.rows,
                self.coeff.as_mut(),
            )?;
            c.tree_sum_f64(self.coeff.as_slice(), &mut self.reduced[..m * input])
                .map_err(error)?;
            self.coeff
                .as_mut_slice()
                .copy_from_slice(&self.reduced[..m * input]);
            rotate(
                &mut self.gemm,
                &self.v,
                &self.coeff,
                &mut self.x,
                m,
                input,
                r,
                "repair_subtract",
            )?;
            for (e, x) in self.e.as_mut_slice().iter_mut().zip(self.x.as_slice()) {
                *e -= x;
            }
        }
        let mut kept = 0;
        for j in 0..input {
            for _ in 0..2 {
                self.norms[..kept].fill(0.0);
                for a in 0..kept {
                    for i in 0..self.rows {
                        self.norms[a] += self.e[(i, a)] * self.e[(i, j)];
                    }
                }
                if kept > 0 {
                    c.tree_sum_f64(&self.norms[..kept], &mut self.sums[..kept])
                        .map_err(error)?;
                    for a in 0..kept {
                        for i in 0..self.rows {
                            self.e[(i, j)] -= self.sums[a] * self.e[(i, a)];
                        }
                    }
                }
            }
            self.norms[0] = (0..self.rows).map(|i| self.e[(i, j)].powi(2)).sum();
            c.tree_sum_f64(&self.norms[..1], &mut self.sums[..1])
                .map_err(error)?;
            let norm = self.sums[0].sqrt();
            if norm <= dft_numerics::linalg::DEPENDENCE_TOLERANCE || m + kept == self.cap {
                r.dropped_corrections += 1;
                continue;
            }
            for i in 0..self.rows {
                self.e[(i, kept)] = self.e[(i, j)] / norm;
            }
            kept += 1;
        }
        self.e.reshape_within_storage(self.rows, kept)?;
        self.append_e(m)?;
        r.seconds[5] += elapsed(t);
        if kept > 0 {
            self.border(c, m, m + kept, false, r, profile)?;
            if !self.factor(m, m + kept, r, profile)? {
                return Err(SolverError::Application(
                    "reference overlap remains dependent after repair".into(),
                ));
            }
        }
        Ok(kept)
    }
    pub fn contract<C: Collective>(
        &mut self,
        c: &C,
        m: usize,
        r: &mut Iteration,
        profile: bool,
    ) -> Result<(), SolverError> {
        let t = clock(profile);
        rotate(
            &mut self.gemm,
            &self.v,
            &self.lean.c,
            &mut self.x,
            m,
            self.n,
            r,
            "restart_v",
        )?;
        rotate(
            &mut self.gemm,
            &self.hv,
            &self.lean.c,
            &mut self.hx,
            m,
            self.n,
            r,
            "restart_hv",
        )?;
        self.v.reshape_within_storage(self.rows, self.n)?;
        self.hv.reshape_within_storage(self.rows, self.n)?;
        self.v.as_mut_slice().copy_from_slice(self.x.as_slice());
        self.hv.as_mut_slice().copy_from_slice(self.hx.as_slice());
        self.lean.a.reshape_within_storage(self.n, self.n)?;
        self.lean.c.reshape_within_storage(self.n, self.n)?;
        self.y.reshape_within_storage(self.n, self.n)?;
        for j in 0..self.n {
            for i in 0..self.n {
                let id = f64::from(i == j);
                let h = if i == j { self.values[i] } else { 0.0 };
                self.hr[(i, j)] = h;
                self.lean.a[(i, j)] = h;
                self.lean.s[(i, j)] = id;
                self.lean.t[(i, j)] = id;
                self.lean.c[(i, j)] = id;
                self.y[(i, j)] = id;
            }
        }
        r.restarts += 1;
        r.seconds[6] += elapsed(t);
        self.refresh(c, self.n, r, profile)?;
        Ok(())
    }
    fn refresh<C: Collective>(
        &mut self,
        c: &C,
        m: usize,
        r: &mut Iteration,
        profile: bool,
    ) -> Result<(), SolverError> {
        self.border(c, 0, m, false, r, profile)?;
        self.border(c, 0, m, true, r, profile)?;
        if !self.factor(0, m, r, profile)? {
            return Err(SolverError::Application(
                "reference contracted basis is dependent".into(),
            ));
        }
        self.h.reshape_within_storage(m, m)?;
        for j in 0..m {
            for i in 0..m {
                self.h[(i, j)] = self.hr[(i, j)];
            }
        }
        let t = clock(profile);
        r.gemm("refresh_transform", m, m, m);
        r.gemm("refresh_transform", m, m, m);
        lean::reduced_border(&mut self.gemm, &self.h, &mut self.lean, 0, m)?;
        r.seconds[2] += elapsed(t);
        Ok(())
    }
    pub fn finish(
        &mut self,
        m: usize,
        hx: bool,
        r: &mut Iteration,
        profile: bool,
    ) -> Result<(), SolverError> {
        let t = clock(profile);
        rotate(
            &mut self.gemm,
            &self.v,
            &self.lean.c,
            &mut self.x,
            m,
            self.n,
            r,
            "return_x",
        )?;
        if hx {
            rotate(
                &mut self.gemm,
                &self.hv,
                &self.lean.c,
                &mut self.hx,
                m,
                self.n,
                r,
                "return_hx",
            )?;
        }
        r.seconds[7] += elapsed(t);
        Ok(())
    }
}
#[allow(clippy::too_many_arguments)]
fn rotate(
    g: &mut RealGemmWorkspace,
    v: &RealMatrix,
    c: &RealMatrix,
    out: &mut RealMatrix,
    m: usize,
    k: usize,
    r: &mut Iteration,
    stage: &'static str,
) -> Result<(), SolverError> {
    out.reshape_within_storage(v.rows(), k)?;
    r.gemm(stage, v.rows(), k, m);
    g.multiply_slices_into(
        GemmSite::Rotate,
        v.as_slice(),
        v.rows(),
        m,
        v.rows(),
        RealOp::Identity,
        c.as_slice(),
        k,
        m,
        out.as_mut(),
    )?;
    Ok(())
}
fn agree_apply<C: Collective>(
    c: &C,
    result: Result<(), SolverError>,
    y: &RealMatrix,
) -> Result<(), SolverError> {
    let failed = result.is_err() || !y.as_slice().iter().all(|v| v.is_finite());
    let mut any = [0.0];
    c.tree_sum_f64(&[f64::from(failed)], &mut any)
        .map_err(error)?;
    if any[0] != 0.0 {
        return Err(SolverError::Application(
            "reference H application failed on a rank".into(),
        ));
    }
    Ok(())
}
