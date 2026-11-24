use super::accounting::{clock, elapsed};
use super::work::{self, Stage, span};
use super::{Iteration, Options, ReducedPlacement, error, preconditioner_denominator, stream};
use crate::{SolverError, distributed_davidson::DistributedRealOperator, real_davidson::lean};
use dft_mpi::Collective;
use dft_numerics::{
    GemmSite, RealGemmWorkspace, RealMatrix, RealOp, pool::WorkerPool, real::RealEighWorkspace,
};

// One bounded cache per calling rank/thread. Only initialized, non-returned
// row matrices are retained; numerical factors and operator state never are.
struct Rows {
    key: (usize, usize, usize),
    v: RealMatrix,
    hv: RealMatrix,
    e: RealMatrix,
    he: RealMatrix,
}
thread_local! {
    static ROWS: std::cell::RefCell<Option<Rows>> = const { std::cell::RefCell::new(None) };
}

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
    single_gemm: bool,
    swap_buffers: bool,
    lean_reset: bool,
    ideal_reset: bool,
    reset_disabled: bool,
    contractions: usize,
    pub sqrt_precondition: bool,
    pub diagonal_bounds: (f64, f64),
    fused_append: bool,
    check_append: bool,
    reuse_rows: bool,
    staged_append: Option<usize>,
    pub reused_elements: usize,
}
impl Workspace {
    pub fn new<O: DistributedRealOperator>(
        op: &O,
        n: usize,
        cap: usize,
        single_gemm: bool,
    ) -> Result<Self, SolverError> {
        let rows = op.local_rows().len();
        let pool = op
            .worker_pool()
            .map(Ok)
            .unwrap_or_else(|| WorkerPool::new(op.inner_threads()))
            .map_err(|e| SolverError::Application(e.to_string()))?;
        use dft_numerics::switches::{Id, on};
        let reuse_rows = on(Id::REF_DAV_REUSE_ROWS);
        let cached = if reuse_rows {
            ROWS.with(|slot| slot.borrow_mut().take())
                .filter(|r| r.key == (rows, n, cap))
        } else {
            None
        };
        let reused_elements = if cached.is_some() {
            2 * rows * cap + 2 * rows * n
        } else {
            0
        };
        let Rows {
            mut v,
            mut hv,
            mut e,
            mut he,
            ..
        } = match cached {
            Some(buffers) => buffers,
            None => Rows {
                key: (rows, n, cap),
                v: RealMatrix::zeros(rows, cap)?,
                hv: RealMatrix::zeros(rows, cap)?,
                e: RealMatrix::zeros(rows, n)?,
                he: RealMatrix::zeros(rows, n)?,
            },
        };
        v.reshape_within_storage(rows, cap)?;
        hv.reshape_within_storage(rows, cap)?;
        e.reshape_within_storage(rows, n)?;
        he.reshape_within_storage(rows, n)?;
        Ok(Self {
            sqrt_precondition: on(Id::REF_DAV_SQRT),
            diagonal_bounds: (f64::INFINITY, f64::NEG_INFINITY),
            fused_append: on(Id::REF_DAV_FUSED_APPEND),
            check_append: on(Id::REF_DAV_CHECK_APPEND),
            reuse_rows,
            reused_elements,
            staged_append: None,
            rows,
            n,
            cap,
            single_gemm,
            ideal_reset: on(Id::REF_DAV_IDEAL_RESET),
            reset_disabled: false,
            contractions: 0,
            lean_reset: dft_numerics::switches::on(dft_numerics::switches::Id::REF_DAV_LEAN_RESET),
            swap_buffers: dft_numerics::switches::on(
                dft_numerics::switches::Id::REF_DAV_SWAP_BUFFERS,
            ),
            v,
            hv,
            x: RealMatrix::zeros(rows, cap)?,
            hx: RealMatrix::zeros(rows, cap)?,
            e,
            he,
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
    pub fn recycle_rows(&mut self) {
        if self.reuse_rows {
            let take = |matrix: &mut RealMatrix| {
                std::mem::replace(matrix, RealMatrix::zeros(0, 0).expect("empty matrix"))
            };
            let rows = Rows {
                key: (self.rows, self.n, self.cap),
                v: take(&mut self.v),
                hv: take(&mut self.hv),
                e: take(&mut self.e),
                he: take(&mut self.he),
            };
            ROWS.with(|slot| *slot.borrow_mut() = Some(rows));
        }
    }
    pub fn normalize_restart<C: Collective>(
        &mut self,
        c: &C,
        norms: &[f64],
        sums: &mut [f64],
        r: &mut Iteration,
        profile: bool,
    ) -> Result<(), SolverError> {
        self.staged_append = stream::normalize(
            c,
            &mut self.e,
            &mut self.v,
            norms,
            sums,
            self.n,
            self.cap,
            self.fused_append,
            r,
            profile,
        )?;
        Ok(())
    }
    pub fn apply<C: Collective, O: DistributedRealOperator>(
        &mut self,
        c: &C,
        op: &mut O,
        x: &RealMatrix,
        _initial: bool,
        r: &mut Iteration,
    ) -> Result<(), SolverError> {
        let _stage = span(r.profile, Stage::ApplyH);
        r.h_calls += 1;
        r.h_columns += x.cols();
        let result = op.apply_block(c, x, &mut self.hv);
        agree_apply(c, result, &self.hv, r.profile)
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
        let _stage = span(r.profile, Stage::ApplyH);
        r.h_calls += 1;
        r.h_columns += k;
        let result = op.apply_block(c, &self.e, &mut self.he);
        if self.check_append {
            let _check = span(r.profile, Stage::ApplyCheck);
            self.hv.reshape_within_storage(self.rows, m + k)?;
            let dst = &mut self.hv.as_mut_slice()[m * self.rows..(m + k) * self.rows];
            // Copy and finite classification share a load. Bitwise OR forces
            // all rows to be visited, allowing vectorization without short circuit.
            let mut failed = result.is_err();
            for (dst, &value) in dst.iter_mut().zip(self.he.as_slice()) {
                *dst = value;
                failed |= !value.is_finite();
            }
            work::copy(r.profile, self.rows * k);
            agree_failure(c, failed)?;
        } else {
            agree_apply(c, result, &self.he, r.profile)?;
            let _append = span(r.profile, Stage::Append);
            work::copy(r.profile, self.rows * k);
            self.hv.reshape_within_storage(self.rows, m + k)?;
            self.hv.as_mut_slice()[m * self.rows..(m + k) * self.rows]
                .copy_from_slice(self.he.as_slice());
        }
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
        let _stage = span(profile, if is_h { Stage::BorderH } else { Stage::BorderS });
        let k = to - from;
        self.border.reshape_within_storage(to, k)?;
        let b = if is_h { &self.hv } else { &self.v };
        r.gemm(if is_h { "border_h" } else { "border_s" }, to, k, self.rows);
        let entry = if self.single_gemm {
            RealGemmWorkspace::multiply_slices_single_into
        } else {
            RealGemmWorkspace::multiply_slices_into
        };
        entry(
            &mut self.gemm,
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
        work::add(profile, |s| {
            let entries = (to * k - k * (k - 1) / 2) as u64;
            s.copied_bytes += 16 * entries;
            s.streamed_bytes += 24 * entries;
        });
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
        let _stage = span(profile, Stage::Factor);
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
            work::add(profile, |s| {
                // Nominal full-panel scalar model of the shared factor adapter.
                // Early rejected panels are explicitly flagged, not exact work.
                let b = width as u64;
                let p = p as u64;
                s.scalar_flops += (2 * b * b * b + 3 * b * b + b) / 3
                    + if p > 0 { b * b + p * b } else { 0 }
                    + 2 * p * b
                    + b * (b + 1)
                    + 2 * b;
                s.sqrt_calls += b + b * (2 * p + b + 1) / 2;
                s.copied_bytes += 16 * b * b;
                s.zeroed_bytes += 8 * (p * b + b * (b - 1) / 2);
                s.streamed_bytes += 32 * b * b + 8 * (p * b + b * (b - 1) / 2);
                s.rejected_factor_panels += u64::from(!accepted);
            });
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
        let _stage = span(o.profile, Stage::Eigen);
        let transform = span(o.profile, Stage::Transform);
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
        work::copy(
            o.profile,
            m * m + from * from + m * (m - from) + (m - from) * (m + from - 1) / 2,
        );
        drop(transform);
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
        let _coefficients = span(o.profile, Stage::Coefficients);
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
        let _stage = span(profile, Stage::Corrections);
        let k = active.len();
        self.coeff.reshape_within_storage(m, k)?;
        for (j, &band) in active.iter().enumerate() {
            self.coeff.as_mut_slice()[j * m..(j + 1) * m]
                .copy_from_slice(&self.lean.c.as_slice()[band * m..(band + 1) * m]);
        }
        work::copy(profile, m * k);
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
        let precondition = span(profile, Stage::Precondition);
        if self.sqrt_precondition {
            let mut sqrt_columns = 0;
            for (j, &band) in active.iter().enumerate() {
                let rows = j * self.rows..(j + 1) * self.rows;
                let x = &self.x.as_slice()[rows.clone()];
                let hx = &self.hx.as_slice()[rows.clone()];
                let out = &mut self.e.as_mut_slice()[rows];
                let safe = stream::sqrt_safe(self.diagonal_bounds, self.values[band]);
                self.norms[j] = if safe {
                    stream::column::<true>(hamiltonian_diagonal, x, hx, out, self.values[band])
                } else {
                    stream::column::<false>(hamiltonian_diagonal, x, hx, out, self.values[band])
                };
                sqrt_columns += usize::from(safe);
            }
            stream::count(profile, self.rows, sqrt_columns, k - sqrt_columns);
        } else {
            self.norms[..k].fill(0.0);
            for (j, &band) in active.iter().enumerate() {
                for (i, &diagonal) in hamiltonian_diagonal.iter().enumerate() {
                    let z = (self.hx[(i, j)] - self.values[band] * self.x[(i, j)])
                        / preconditioner_denominator(diagonal, self.values[band]);
                    self.e[(i, j)] = z;
                    self.norms[j] += z * z;
                }
            }
            work::add(profile, |s| {
                s.scalar_flops += 13 * (self.rows * k) as u64;
                s.hypot_calls += (self.rows * k) as u64;
                s.zeroed_bytes += 8 * k as u64;
                s.streamed_bytes += 32 * (self.rows * k) as u64 + 8 * k as u64;
            });
        }
        drop(precondition);
        self.staged_append = stream::normalize(
            c,
            &mut self.e,
            &mut self.v,
            &self.norms,
            &mut self.sums,
            m,
            self.cap,
            self.fused_append,
            r,
            profile,
        )?;
        r.seconds[4] += elapsed(t);
        Ok(())
    }
    pub fn append_e(&mut self, m: usize, profile: bool) -> Result<(), SolverError> {
        if self.staged_append.take() == Some(m) {
            return Ok(());
        }
        let _stage = span(profile, Stage::Append);
        work::copy(profile, self.rows * self.e.cols());
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
        let _stage = span(profile, Stage::Repair);
        self.staged_append = None;
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
        if self.swap_buffers {
            std::mem::swap(&mut self.v, &mut self.x);
            std::mem::swap(&mut self.hv, &mut self.hx);
        } else {
            work::copy(profile, 2 * self.v.rows() * self.v.cols());
            self.v.as_mut_slice().copy_from_slice(self.x.as_slice());
            self.hv.as_mut_slice().copy_from_slice(self.hx.as_slice());
        }
        work::copy(profile, m * m + m * self.n);
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
            let entry = if self.single_gemm {
                RealGemmWorkspace::multiply_slices_single_into
            } else {
                RealGemmWorkspace::multiply_slices_into
            };
            entry(
                &mut self.gemm,
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
            work::copy(profile, m * input);
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
            work::add(profile, |s| {
                s.scalar_flops += (self.rows * input) as u64;
                s.streamed_bytes += 24 * (self.rows * input) as u64;
            });
            for (e, x) in self.e.as_mut_slice().iter_mut().zip(self.x.as_slice()) {
                *e -= x;
            }
        }
        let panel = span(profile, Stage::RepairPanel);
        let mut kept = 0;
        for j in 0..input {
            work::add(profile, |s| {
                s.scalar_flops += (8 * self.rows * kept + 2 * self.rows) as u64;
                s.sqrt_calls += 1;
                s.zeroed_bytes += 16 * kept as u64;
                s.streamed_bytes += (80 * self.rows * kept + 8 * self.rows + 16 * kept) as u64;
            });
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
        work::add(profile, |s| {
            s.scalar_flops += (self.rows * kept) as u64;
            s.streamed_bytes += 16 * (self.rows * kept) as u64;
        });
        drop(panel);
        self.append_e(m, profile)?;
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
        let _stage = span(profile, Stage::Restart);
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
        if self.swap_buffers {
            std::mem::swap(&mut self.v, &mut self.x);
            std::mem::swap(&mut self.hv, &mut self.hx);
        } else {
            work::copy(profile, 2 * self.v.rows() * self.v.cols());
            self.v.as_mut_slice().copy_from_slice(self.x.as_slice());
            self.hv.as_mut_slice().copy_from_slice(self.hx.as_slice());
        }
        self.lean.a.reshape_within_storage(self.n, self.n)?;
        self.lean.c.reshape_within_storage(self.n, self.n)?;
        self.y.reshape_within_storage(self.n, self.n)?;
        self.contractions += 1;
        let mut ideal = false;
        if self.ideal_reset {
            r.resets.attempts += 1;
            if self.reset_disabled {
                r.resets.disabled_refreshes += 1;
            } else if self
                .contractions
                .is_multiple_of(super::ideal_reset::REFRESH_INTERVAL)
            {
                r.resets.periodic_refreshes += 1;
            } else if self.reset_norm_guard(c, profile)? {
                ideal = true;
                r.resets.ideal += 1;
            } else {
                self.reset_disabled = true;
                r.resets.guard_fallbacks += 1;
            }
        }
        #[cfg(test)]
        super::ideal_reset::audit_restart(c, &self.v);
        for j in 0..self.n {
            for i in 0..self.n {
                let id = f64::from(i == j);
                if self.lean_reset && !ideal {
                    self.lean.c[(i, j)] = id;
                } else {
                    let h = if i == j { self.values[i] } else { 0.0 };
                    self.hr[(i, j)] = h;
                    self.lean.a[(i, j)] = h;
                    self.lean.s[(i, j)] = id;
                    self.lean.t[(i, j)] = id;
                    self.lean.c[(i, j)] = id;
                    self.y[(i, j)] = id;
                }
            }
        }
        work::add(profile, |s| {
            let matrices = if self.lean_reset && !ideal { 1 } else { 6 };
            s.zeroed_bytes += matrices * 8 * (self.n * (self.n - 1)) as u64;
            s.streamed_bytes += matrices * 8 * (self.n * self.n) as u64;
        });
        r.restarts += 1;
        r.seconds[6] += elapsed(t);
        if !ideal {
            self.refresh(c, self.n, r, profile)?;
            // Periodic full S checks also catch unit-norm but correlated columns.
            // Once unsafe, retain measured arithmetic for the rest of this solve.
            if self.ideal_reset && !self.reset_disabled {
                let t = clock(profile);
                let mut defect: f64 = 0.0;
                for j in 0..self.n {
                    for i in 0..self.n {
                        let d = (self.lean.s[(i, j)] - f64::from(i == j)).abs();
                        if !d.is_finite() || d > super::ideal_reset::ORTHOGONALITY_LIMIT {
                            self.reset_disabled = true;
                        }
                        defect = defect.max(d);
                    }
                }
                work::add(profile, |s| {
                    s.scalar_flops += (self.n * self.n) as u64;
                    s.streamed_bytes += 8 * (self.n * self.n) as u64;
                });
                r.resets.max_measured_defect = r.resets.max_measured_defect.max(defect);
                if self.reset_disabled {
                    r.resets.guard_fallbacks += 1;
                }
                r.seconds[6] += elapsed(t);
            }
        }
        Ok(())
    }
    /// One O(rows * bands) pass and one bands-wide sum; no scratch allocation.
    fn reset_norm_guard<C: Collective>(
        &mut self,
        c: &C,
        profile: bool,
    ) -> Result<bool, SolverError> {
        for j in 0..self.n {
            self.norms[j] = self.v.as_slice()[j * self.rows..(j + 1) * self.rows]
                .iter()
                .map(|v| v * v)
                .sum();
        }
        c.tree_sum_f64(&self.norms[..self.n], &mut self.sums[..self.n])
            .map_err(error)?;
        work::add(profile, |s| {
            s.scalar_flops += (2 * self.rows * self.n + self.n) as u64;
            s.streamed_bytes += 8 * (self.rows * self.n + 2 * self.n) as u64;
        });
        Ok(self.sums[..self.n]
            .iter()
            .all(|v| v.is_finite() && (v - 1.0).abs() <= super::ideal_reset::ORTHOGONALITY_LIMIT))
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
        let _transform = span(profile, Stage::Transform);
        work::copy(profile, 2 * m * m + m * (m - 1) / 2);
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
        let _stage = span(profile, Stage::Return);
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
    // The plain entry falls back to the pure-Rust kernel whenever the whole
    // output (rows x k) exceeds the scratch square, which is every rotation of
    // a real cell (dev08 si128: 186 GFLOP per rank on that route).
    let entry = if dft_numerics::switches::on(dft_numerics::switches::Id::REF_DAV_WIDE_ROTATE) {
        RealGemmWorkspace::multiply_slices_wide_into
    } else {
        RealGemmWorkspace::multiply_slices_into
    };
    entry(
        g,
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
    profile: bool,
) -> Result<(), SolverError> {
    let _stage = span(profile, Stage::ApplyCheck);
    work::add(profile, |s| {
        s.streamed_bytes += 8 * y.as_slice().len() as u64
    });
    let failed = result.is_err() || !y.as_slice().iter().all(|v| v.is_finite());
    agree_failure(c, failed)
}
fn agree_failure<C: Collective>(c: &C, failed: bool) -> Result<(), SolverError> {
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

#[cfg(test)]
mod single_gemm_tests {
    use super::*;
    use crate::ref_davidson::{Solution, solve};
    use dft_numerics::switches::RunConfig;

    // Same diagonal-plus-rank-one, clustered and cyclic-row fixtures as the
    // reference restart-correction diagnostics; add a multi-panel local shape.
    struct Operator {
        dim: usize,
        rows: Vec<usize>,
        clustered: bool,
        threads: usize,
    }
    impl Operator {
        fn diagonal(&self, row: usize) -> f64 {
            if self.clustered && row < 2 {
                1.0
            } else {
                1.0 + row as f64 * 0.4
            }
        }
    }
    impl DistributedRealOperator for Operator {
        fn global_len(&self) -> usize {
            self.dim
        }
        fn local_rows(&self) -> &[usize] {
            &self.rows
        }
        fn inner_threads(&self) -> usize {
            self.threads
        }
        fn apply_block<C: Collective>(
            &mut self,
            c: &C,
            x: &RealMatrix,
            y: &mut RealMatrix,
        ) -> Result<(), SolverError> {
            for j in 0..x.cols() {
                let dot: f64 = self
                    .rows
                    .iter()
                    .enumerate()
                    .map(|(i, &row)| ((row + 1) as f64).sin() * x[(i, j)])
                    .sum();
                let mut sum = [0.0];
                c.tree_sum_f64(&[dot], &mut sum).map_err(error)?;
                for (i, &row) in self.rows.iter().enumerate() {
                    y[(i, j)] =
                        self.diagonal(row) * x[(i, j)] + 0.03 * ((row + 1) as f64).sin() * sum[0];
                }
            }
            Ok(())
        }
    }
    fn bits(x: &[f64]) -> Vec<u64> {
        x.iter().map(|v| v.to_bits()).collect()
    }
    fn decisions(s: &Solution) -> String {
        let r = &s.record;
        let stages: Vec<_> = std::iter::once(&r.setup)
            .chain(&r.iterations)
            .chain(std::iter::once(&r.finish))
            .map(|r| {
                let mut r = r.clone();
                r.seconds.fill(0.0);
                r
            })
            .collect();
        format!(
            "{:?} {:?} {stages:?}",
            s.converged,
            (
                r.communication,
                r.h_calls,
                r.h_columns,
                r.reduced_solves,
                r.backend_solves,
                r.restarts,
                r.gemm_flops,
                r.placement
            )
        )
    }
    fn repeat(a: &Solution, b: &Solution) {
        assert_eq!(decisions(a), decisions(b));
        assert_eq!(bits(&a.eigenvalues), bits(&b.eigenvalues));
        assert_eq!(
            bits(a.eigenvectors.as_slice()),
            bits(b.eigenvectors.as_slice())
        );
        assert_eq!(
            bits(a.hx.as_ref().unwrap().as_slice()),
            bits(b.hx.as_ref().unwrap().as_slice())
        );
        assert_eq!(a.ritz_changes_ha.len(), b.ritz_changes_ha.len());
        for (a, b) in a.ritz_changes_ha.iter().zip(&b.ritz_changes_ha) {
            assert_eq!(bits(a), bits(b));
        }
        let ritz_bits = |s: &Solution| {
            s.ritz_values_ha
                .iter()
                .map(|step| {
                    step.iter()
                        .map(|pair| pair.map(f64::to_bits))
                        .collect::<Vec<_>>()
                })
                .collect::<Vec<_>>()
        };
        assert_eq!(ritz_bits(a), ritz_bits(b));
        assert_eq!(bits(&a.band_thresholds_ha), bits(&b.band_thresholds_ha));
    }

    #[test]
    fn single_gemm_decisions_columns_repeats_and_threads() {
        let off = RunConfig::from_pairs(&[("ref_dav_single_gemm", "off")]).unwrap();
        let on = RunConfig::from_pairs(&[("ref_dav_single_gemm", "on")]).unwrap();
        for (dim, n, clustered) in [
            (12, 3, false),
            (17, 5, true),
            (7, 4, true),
            (1025, 3, false),
        ] {
            for ranks in [1, 3, 8] {
                dft_mpi::simulate(ranks, |c| {
                    let mut op = Operator {
                        dim,
                        rows: (c.rank()..dim).step_by(ranks).collect(),
                        clustered,
                        threads: 1,
                    };
                    let mut x = RealMatrix::zeros(op.rows.len(), n).unwrap();
                    for j in 0..n {
                        for (i, &row) in op.rows.iter().enumerate() {
                            x[(i, j)] = if row == j {
                                1.0
                            } else {
                                0.03 * ((row * 7 + j * 13 + 1) as f64).sin()
                            };
                        }
                    }
                    let diag: Vec<_> = op
                        .rows
                        .iter()
                        .map(|&row| op.diagonal(row) + 0.03 * ((row + 1) as f64).sin().powi(2))
                        .collect();
                    for placement in [ReducedPlacement::Root, ReducedPlacement::Replicated] {
                        let mut o = Options::new(n, 1e-10, n - 1);
                        o.max_iterations = if dim > 1000 { 2 } else { 40 };
                        o.profile = true;
                        o.return_hx = true;
                        o.placement = placement;
                        let a = {
                            let _scope = off.enter();
                            solve(c, &mut op, &x, &diag, o).unwrap()
                        };
                        let _scope = on.enter();
                        let b = solve(c, &mut op, &x, &diag, o).unwrap();
                        assert_eq!(decisions(&a), decisions(&b), "dim={dim}, ranks={ranks}");
                        for (a, b) in a.eigenvalues.iter().zip(&b.eigenvalues) {
                            assert!((a - b).abs() < 1e-10, "eigenvalue difference in Ha");
                        }
                        if !dft_numerics::MKL_BACKEND {
                            repeat(&a, &b);
                        }
                        repeat(&b, &solve(c, &mut op, &x, &diag, o).unwrap());
                        op.threads = 4;
                        repeat(&b, &solve(c, &mut op, &x, &diag, o).unwrap());
                        op.threads = 1;
                    }
                });
            }
        }
    }

    #[test]
    fn single_gemm_repair_projection_preserves_filter_and_h_work() {
        dft_mpi::simulate(1, |c| {
            let op = Operator {
                dim: 257,
                rows: (0..257).collect(),
                clustered: false,
                threads: 1,
            };
            let mut prior = None;
            for single in [false, true, true] {
                let mut w = Workspace::new(&op, 2, 6, single).unwrap();
                w.v.reshape_within_storage(257, 2).unwrap();
                w.hv.reshape_within_storage(257, 2).unwrap();
                w.lean.c.reshape_within_storage(2, 2).unwrap();
                for j in 0..2 {
                    w.v[(j, j)] = 1.0;
                    w.hv[(j, j)] = (j + 1) as f64;
                    w.lean.s[(j, j)] = 1.0;
                    w.lean.t[(j, j)] = 1.0;
                    w.lean.c[(j, j)] = 1.0;
                    w.e[(0, j)] = 0.5;
                    w.e[(2, j)] = 0.5;
                }
                let mut r = super::super::Record::new(1, 6, false, ReducedPlacement::Replicated);
                assert_eq!(w.repair(c, 2, &mut r.setup, false).unwrap(), 1);
                assert_eq!(r.setup.h_columns, 0);
                assert_eq!(r.setup.h_calls, 0);
                assert_eq!(r.setup.dropped_corrections, 1);
                assert_eq!(
                    r.setup
                        .gemms
                        .iter()
                        .filter(|g| g.stage == "repair_projection")
                        .count(),
                    2
                );
                let state = (bits(w.e.as_slice()), format!("{:?}", r.setup));
                if let Some(ref prior) = prior {
                    assert_eq!(&state, prior);
                }
                prior = Some(state);
            }
        });
    }
    #[test]
    fn ideal_reset_guard_fallback_and_periodic_refresh() {
        for ranks in [1, 3, 8] {
            for correlated in [false, true] {
                dft_mpi::simulate(ranks, |c| {
                    let _scope = RunConfig::from_pairs(&[("ref_dav_ideal_reset", "on")])
                        .unwrap()
                        .enter();
                    let op = Operator {
                        dim: 7,
                        rows: (c.rank()..7).step_by(ranks).collect(),
                        clustered: false,
                        threads: 1,
                    };
                    let mut w = Workspace::new(&op, 4, 7, false).unwrap();
                    w.v.reshape_within_storage(op.rows.len(), 4).unwrap();
                    w.hv.reshape_within_storage(op.rows.len(), 4).unwrap();
                    w.lean.c.reshape_within_storage(4, 4).unwrap();
                    for j in 0..4 {
                        w.values[j] = 1.0;
                        w.lean.c[(j, j)] = 1.0;
                        for (i, &row) in op.rows.iter().enumerate() {
                            let v = if correlated && j == 1 {
                                if row == 0 {
                                    1e-4
                                } else if row == 1 {
                                    (1.0_f64 - 1e-8).sqrt()
                                } else {
                                    0.0
                                }
                            } else if row == j {
                                if !correlated && j == 0 { 1.1 } else { 1.0 }
                            } else {
                                0.0
                            };
                            w.v[(i, j)] = v;
                            w.hv[(i, j)] = v;
                        }
                    }
                    let mut record =
                        crate::ref_davidson::Record::new(8, 7, false, ReducedPlacement::Root);
                    // Unit-norm correlated columns pass the diagonal guard, then
                    // the periodic full S catches the off-diagonal defect.
                    let trip = if correlated { 4 } else { 1 };
                    for step in 0..trip + 1 {
                        w.contract(c, 4, &mut record.iterations[step], false)
                            .unwrap();
                    }
                    record.summarize();
                    let stats = record.resets;
                    for r in &record.iterations[..trip + 1] {
                        let refresh_products = r
                            .gemms
                            .iter()
                            .filter(|g| {
                                matches!(g.stage, "border_s" | "border_h" | "refresh_transform")
                            })
                            .count();
                        assert_eq!(refresh_products, if r.resets.ideal > 0 { 0 } else { 4 });
                        assert_eq!(r.h_columns, 0);
                        assert_eq!(r.reduced_solves, 0);
                    }
                    assert_eq!(stats.guard_fallbacks, 1);
                    assert_eq!(stats.disabled_refreshes, 1);
                    assert_eq!(stats.periodic_refreshes, usize::from(correlated));
                    assert_eq!(stats.ideal, if correlated { 3 } else { 0 });
                    assert!(w.reset_disabled);
                    // Fallback actually measured and factored S; it is not an
                    // identity reset with only a counter claiming a refresh.
                    if correlated {
                        assert!((w.lean.s[(0, 1)] - 1e-4).abs() < 1e-14);
                    } else {
                        assert!((w.lean.s[(0, 0)] - 1.21).abs() < 1e-14);
                    }
                    let expected = if correlated { 1.0 } else { 1.0 / 1.1 };
                    assert!((w.lean.t[(0, 0)] - expected).abs() < 1e-14);
                    // A non-finite norm is rejected collectively even on ranks
                    // with no local rows. This test calls only the norm guard.
                    if let Some(v) = w.v.as_mut_slice().first_mut() {
                        *v = f64::NAN;
                    }
                    assert!(!w.reset_norm_guard(c, false).unwrap());
                });
            }
        }
    }
}
