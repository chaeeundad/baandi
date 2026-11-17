//! Lean Davidson of the serial real Gamma path (`BAANDI_DAVIDSON=lean`, work
//! order beat-QE 2026-11-14 phase B; default `classic`). Acts on Ritz-value
//! (standard-profile) solves only; residual (strict) solves and the
//! certificate keep the classic path. Derived from the study notes
//! `docs/research/qe-study/qe-solvers-and-subspace.md` (sections 1, 2.6, 4) and
//! `davidson-and-scf-thresholds.md` section 1; no QE source was read.
//!
//! Per iteration, with basis V (n x m, unit-norm columns, not orthonormal), its
//! H images HV, the reduced matrix H_r = V^T H V and the overlap S = V^T V:
//!
//! * B3a, generalized reduced problem: S is kept beside H_r. Its Cholesky factor
//!   S = L L^T is kept as the explicit inverse T = L^{-1} (lower triangular),
//!   extended block-wise when columns are appended (bordered Cholesky in panels
//!   of `FACTOR_PANEL` columns, all products through the fixed-panel GEMM
//!   workspace). The reduced problem H_r c = theta S c is solved as the standard
//!   problem A y = theta y with A = T H_r T^T (dsyevd / faer, or the B5 secular
//!   update, because a one-column append borders A exactly as it borders H_r),
//!   and c = T^T y. Appended corrections are only normalized to unit length:
//!   the one extra border GEMM V^T E replaces the projection GEMMs Q^T X and
//!   Q S of the classic path (and its in-block Gram-Schmidt).
//! * B1 at every width: residual vectors of the unconverged roots only
//!   (X_s = V C_s, H X_s = HV C_s); Ritz-change flags are exact, so no residual
//!   of a converged root is needed. Ritz vectors and H images of all bands are
//!   formed only at the return (with the residual norms the standard stop
//!   guard reads) and at a contraction restart (QE's update / last split).
//! * Fallback (fixed numerical rule, no material constant): every new column's
//!   squared Cholesky pivot over its overlap diagonal (the squared sine of its
//!   angle to the earlier columns) must be at least
//!   `PIVOT_RATIO_MIN = eps / ORTHOGONALITY_GATE` (2.2e-4). A Cholesky
//!   reduction loses about eps / ratio of S-orthonormality in the Ritz
//!   vectors; the bound keeps that at the 1e-12 orthogonality gate the classic
//!   path applies at a contraction, and it keeps 11 to 12 correct digits in the
//!   pivot itself (computed as 1 - |w|^2). When an append fails the test the
//!   basis is orthonormalized once by rotating it to its full Ritz basis
//!   (V C with all m Ritz coefficient columns; span unchanged, no H
//!   application), and the raw correction block is projected with the classic
//!   projection, whose per-column rank decision keeps the meaning of
//!   `DEPENDENCE_TOLERANCE`. A fresh basis (start, contraction) failing the test
//!   is rebuilt by projection and H application (counted; expected zero).
//! * Capacity per band (B2) is the caller's `davidson_capacity_per_band`
//!   (input `scf.davidson_capacity_per_band`, default 6, unchanged). B4: the
//!   dense solve has the order of the basis; at 2 to 3 bands per basis it is
//!   the Cholesky reduction plus one dsyevd of order at most 3 x bands.
//!
//! Stage timers: overlap border and fallback rotation under `orthonormalize`,
//! factor extension, reduced border of A and the coefficient transform under
//! `eigen`, return and contraction rotations under `rotate`.
//!
//! Ledger premises (docs/performance/attempts-ledger.md): P37A-1 rejected the
//! non-orthogonal basis by a rule that covered all profiles (strict-level
//! Gram eigenvalues near 1e-16); with one standard path (user decision
//! 2026-11-14) and this module acting on Ritz-value solves only, that rule no
//! longer applies, and P37A-3 measured 45 to 50 percent fewer multiply-adds per
//! standard step. A13 / P31A (capacity) assumed all roots always rotated and an
//! orthonormal basis; B1 at all widths and B3a remove both premises. A52
//! (CholeskyQR2 inside the MGS) is not repeated: the Cholesky acts on the small
//! overlap with a pivot bound and an MGS fallback, and the classic MGS is
//! unchanged. PRT1's stale-flag failure concerned residual flags; Ritz-change
//! flags are recomputed exactly every iteration.
use super::*;
use crate::davidson_stages::{LeanCount, lean_add, lean_pivot};
use dft_numerics::GemmSite;

/// The orthogonality gate of the classic contraction (`orthogonality_defect`
/// <= 1e-12), reused as the target for the S-orthonormality of Ritz vectors.
pub(crate) const ORTHOGONALITY_GATE: f64 = 1e-12;
/// Smallest accepted squared Cholesky pivot over its overlap diagonal.
pub const PIVOT_RATIO_MIN: f64 = f64::EPSILON / ORTHOGONALITY_GATE;
/// Columns per bordered Cholesky panel (serial work inside a panel is
/// O(panel^3); everything across panels is GEMM). A shape constant.
const FACTOR_PANEL: usize = 64;

pub(crate) struct LeanScratch {
    pub(crate) cap: usize,
    /// Basis overlap S, cap x cap, leading dimension cap.
    pub(crate) s: RealMatrix,
    /// Inverse Cholesky factor T = L^{-1}, lower triangular, cap x cap.
    pub(crate) t: RealMatrix,
    /// Reduced matrix T H_r T^T, compact m x m.
    pub(crate) a: RealMatrix,
    /// Ritz coefficients in basis coordinates, compact m x k.
    pub(crate) c: RealMatrix,
    /// Gathered eigenvector columns, compact m x k.
    y: RealMatrix,
    /// Reduced-border scratch: Z = T[new rows, :]^T and H_r Z, m x c.
    pub(crate) z: RealMatrix,
    pub(crate) p: RealMatrix,
    /// Factor scratch: W = T S_12 (p0 x w), U = W^T T (w x p0), W^T W and
    /// the panel inverse (w x w).
    w: RealMatrix,
    u: RealMatrix,
    ww: RealMatrix,
    tinv: RealMatrix,
    g: Vec<f64>,
    l: Vec<f64>,
    sel_values: Vec<f64>,
    sel_norms: Vec<f64>,
    /// LEAN_BORDER=fused: the block [E | H E] (n x 2k, sized on first use)
    /// and the product V^T [E | H E] ((m + k) x 2k, compact), with the column
    /// range of the reduced border it holds for the next iteration.
    pair: RealMatrix,
    fb: RealMatrix,
    pending: Option<(usize, usize)>,
    /// R1-b (`BAANDI_RESTART_KEEP`): history and scratch of the kept
    /// directions, allocated at the first solve with the switch on.
    keep: Option<Box<KeepScratch>>,
}

/// R1-b history-preserving restart (ticket BREAKTHROUGH-1-RESTART,
/// docs/design/restart-keep.md). The history holds the basis coefficients D of
/// the latest Ritz vectors Y = V D of the unconverged roots (the previous
/// iteration's selected roots, or the columns of X right after a contraction
/// or a Ritz-basis fallback rotation). At a contraction with Ritz coefficients
/// C (X = V C, C^T S C = I) the kept directions are
///   A = C^T S D,  G = D - C A,  P = V G,  H P = (H V) G,
/// formed by the same rotation as X (no H application; the X component is
/// removed twice in coefficient space). Candidates are taken greedily in
/// ascending root order; a candidate is kept when |P| > DEPENDENCE_TOLERANCE
/// |Y| (the classic rank floor on a unit input) and its squared pivot against
/// the directions kept before it over |P|^2 is at least `PIVOT_RATIO_MIN` (the
/// pivot rule of every lean append); the kept columns are S-orthonormalized in
/// coefficient space (G L^{-T}). The caller then factors the true overlap of
/// [X | P] with the same pivot rule and falls back to X alone on a failure. All of this is small replicated arithmetic in
/// a fixed order (identical on every rank, independent of the thread count).
pub(crate) struct KeepScratch {
    /// Basis columns the history coefficients refer to (a prefix of the basis).
    hist_m: usize,
    /// The history is expressed in the current basis.
    valid: bool,
    /// Roots of the history columns, ascending.
    roots: Vec<usize>,
    /// History coefficients, compact hist_m x roots.len().
    hist: Vec<f64>,
    /// Per-candidate columns D, S D, S G (m each) and C^T S D (bands).
    d: Vec<f64>,
    sd: Vec<f64>,
    a: Vec<f64>,
    /// S G of the kept columns, compact m x kept.
    sg: Vec<f64>,
    /// Cholesky factor of the kept columns' Gram G^T S G, kept x kept.
    l: Vec<f64>,
}
impl KeepScratch {
    fn new(cap: usize, bands: usize) -> Self {
        Self {
            hist_m: 0,
            valid: false,
            roots: Vec::with_capacity(bands),
            hist: vec![0.; cap * bands],
            d: vec![0.; cap],
            sd: vec![0.; cap],
            a: vec![0.; bands],
            sg: vec![0.; cap * bands],
            l: vec![0.; bands * bands],
        }
    }
}
impl LeanScratch {
    /// Make the R1-b scratch available (first use only).
    pub(crate) fn keep_prepare(&mut self, bands: usize) {
        if self.keep.as_ref().is_none_or(|k| k.a.len() != bands) {
            self.keep = Some(Box::new(KeepScratch::new(self.cap, bands)));
        }
        if let Some(k) = self.keep.as_mut() {
            k.valid = false;
        }
    }
    /// Forget the history (the basis changed in a way it does not follow).
    pub(crate) fn keep_invalidate(&mut self) {
        if let Some(k) = self.keep.as_mut() {
            k.valid = false;
        }
    }
    /// The history becomes the selected Ritz coefficients of this iteration
    /// (`self.c`, compact m x roots.len(), from `ritz_coefficients` with
    /// `Some(roots)`).
    pub(crate) fn keep_record(&mut self, m: usize, roots: &[usize]) {
        let Some(k) = self.keep.as_mut() else {
            return;
        };
        let len = m * roots.len();
        k.hist[..len].copy_from_slice(&self.c.as_slice()[..len]);
        k.roots.clear();
        k.roots.extend_from_slice(roots);
        k.hist_m = m;
        k.valid = true;
    }
    /// The history becomes basis columns: root j is column j (right after a
    /// contraction, whose basis starts with X, or a Ritz-basis rotation).
    pub(crate) fn keep_identity(&mut self, m: usize, roots: &[usize]) {
        let Some(k) = self.keep.as_mut() else {
            return;
        };
        k.roots.clear();
        k.roots.extend(roots.iter().copied().filter(|&j| j < m));
        let len = m * k.roots.len();
        k.hist[..len].fill(0.);
        for (slot, &j) in k.roots.iter().enumerate() {
            k.hist[slot * m + j] = 1.;
        }
        k.hist_m = m;
        k.valid = true;
    }
    /// At a contraction of the m-column basis (S[:m, :m] in `self.s`, C =
    /// `self.c` compact m x bands): append up to `limit` kept direction
    /// coefficient columns G to `self.c` (compact m x (bands + kept)) for the
    /// roots in `roots` (ascending) that the history holds. Returns (kept,
    /// candidates dropped by the pivot rule).
    pub(crate) fn keep_directions(
        &mut self,
        m: usize,
        bands: usize,
        roots: &[usize],
        limit: usize,
    ) -> Result<(usize, usize), SolverError> {
        let cap = self.cap;
        let Some(k) = self.keep.as_mut() else {
            return Ok((0, 0));
        };
        if !k.valid || limit == 0 || k.hist_m > m || bands + limit > cap {
            return Ok((0, 0));
        }
        let hm = k.hist_m;
        let (s, c) = (&self.s, &mut self.c);
        c.reshape_within_storage(m, bands + limit)?;
        let mut kept = 0;
        let mut dropped = 0;
        for &root in roots {
            if kept == limit {
                break;
            }
            let Ok(slot) = k.roots.binary_search(&root) else {
                continue;
            };
            // D (zero below the history prefix), S D and |Y|^2 = D^T S D.
            let d = &mut k.d[..m];
            d[..hm].copy_from_slice(&k.hist[slot * hm..(slot + 1) * hm]);
            d[hm..].fill(0.);
            let sd = &mut k.sd[..m];
            for (i, x) in sd.iter_mut().enumerate() {
                *x = (0..hm).map(|l| s[(i, l)] * d[l]).sum();
            }
            let length: f64 = d.iter().zip(sd.iter()).map(|(a, b)| a * b).sum();
            // G = D - C A with A = C^T S D, twice (the second pass removes
            // the rounding left in the X component; written as column
            // bands + kept of c), then S G.
            let col = bands + kept;
            let (head, tail) = c.as_mut_slice().split_at_mut(col * m);
            let g = &mut tail[..m];
            g.copy_from_slice(d);
            let sg = &mut k.sg[kept * m..(kept + 1) * m];
            for pass in 0..2 {
                let source: &[f64] = if pass == 0 { sd } else { sg };
                for (b, a) in k.a.iter_mut().enumerate() {
                    *a = head[b * m..(b + 1) * m]
                        .iter()
                        .zip(source)
                        .map(|(x, y)| x * y)
                        .sum();
                }
                for (b, &a) in k.a.iter().enumerate() {
                    for (x, &cc) in g.iter_mut().zip(&head[b * m..(b + 1) * m]) {
                        *x -= cc * a;
                    }
                }
                for (i, x) in sg.iter_mut().enumerate() {
                    *x = (0..m).map(|l| s[(i, l)] * g[l]).sum();
                }
            }
            // |P|^2 and the new row of the Cholesky factor of G^T S G.
            let sg = &k.sg[..(kept + 1) * m];
            let norm2 = g
                .iter()
                .zip(&sg[kept * m..])
                .map(|(a, b)| a * b)
                .sum::<f64>();
            let mut pivot = norm2;
            let lw = bands;
            for p in 0..kept {
                let gp = &head[(bands + p) * m..(bands + p + 1) * m];
                let mut x: f64 = gp.iter().zip(&sg[kept * m..]).map(|(a, b)| a * b).sum();
                for q in 0..p {
                    x -= k.l[kept + q * lw] * k.l[p + q * lw];
                }
                let x = x / k.l[p + p * lw];
                k.l[kept + p * lw] = x;
                pivot -= x * x;
            }
            // Rank floor of a unit input (|Y| = 1 in the S norm) and the pivot
            // rule against the directions kept before.
            let ratio = pivot / norm2;
            if !(norm2.is_finite()
                && norm2.sqrt() > DEPENDENCE_TOLERANCE * length.sqrt()
                && ratio >= PIVOT_RATIO_MIN)
            {
                dropped += 1;
                continue;
            }
            k.l[kept + kept * lw] = pivot.sqrt();
            kept += 1;
        }
        // S-orthonormalize the kept columns: G <- G L^{-T}, column by column.
        let lw = bands;
        let data = c.as_mut_slice();
        for q in 0..kept {
            let (done, rest) = data.split_at_mut((bands + q) * m);
            let gq = &mut rest[..m];
            for p in 0..q {
                let lqp = k.l[q + p * lw];
                for (x, &y) in gq
                    .iter_mut()
                    .zip(&done[(bands + p) * m..(bands + p + 1) * m])
                {
                    *x -= lqp * y;
                }
            }
            let inv = 1. / k.l[q + q * lw];
            for x in gq.iter_mut() {
                *x *= inv;
            }
        }
        c.reshape_within_storage(m, bands + kept)?;
        Ok((kept, dropped))
    }
}
impl LeanScratch {
    pub(crate) fn new(cap: usize, bands: usize) -> Result<Self, SolverError> {
        let panel = FACTOR_PANEL.min(cap.max(1));
        Ok(Self {
            cap,
            s: RealMatrix::zeros(cap, cap)?,
            t: RealMatrix::zeros(cap, cap)?,
            a: RealMatrix::zeros(cap, cap)?,
            c: RealMatrix::zeros(cap, cap)?,
            y: RealMatrix::zeros(cap, cap)?,
            z: RealMatrix::zeros(cap, cap)?,
            p: RealMatrix::zeros(cap, cap)?,
            w: RealMatrix::zeros(cap, panel)?,
            u: RealMatrix::zeros(panel, cap)?,
            ww: RealMatrix::zeros(panel, panel)?,
            tinv: RealMatrix::zeros(panel, panel)?,
            g: vec![0.; panel * panel],
            l: vec![0.; panel * panel],
            sel_values: vec![0.; bands],
            sel_norms: vec![0.; bands],
            pair: RealMatrix::zeros(0, 0)?,
            fb: RealMatrix::zeros(cap, 2 * bands)?,
            pending: None,
            keep: None,
        })
    }
}

/// C3a (`BAANDI_LEAN_PROJ_CACHE=on`; ledger C3a, review F10): projector
/// products P^T V of the basis columns, kept beside V and H V. The operator
/// leaves P^T Y of its latest block application available
/// (`RealLinearOperator::copy_last_projections`); every apply of basis columns
/// copies them in (start, appended corrections, a rebuild), a contraction
/// rotates them with the same coefficients as V (P^T (V C) = (P^T V) C) and the
/// return forms P^T X = (P^T V) C for the caller's nonlocal energy, a product
/// of order n_proj x m x bands instead of n x n_proj x bands. Any basis change
/// the cache does not follow (the stage-2 fallback rotation) invalidates it,
/// and an invalid cache is only made valid again by a fresh apply of the whole
/// basis; the return then reports no projections. Products go through the
/// fixed-panel GEMM workspace, so they do not depend on the thread count.
/// The solve's arithmetic is untouched (the cache is only read at the return).
pub(crate) struct ProjCache {
    count: usize,
    /// P^T V, count x capacity (column j belongs to basis column j).
    pv: RealMatrix,
    /// P^T X of the latest return (and contraction scratch), count x bands.
    pub(crate) pc: RealMatrix,
    cols: usize,
    valid: bool,
}
impl ProjCache {
    /// The cache of this solve (allocated on the first use or a shape change),
    /// emptied; `None` when the operator exposes no projector products.
    fn activate<'a>(
        slot: &'a mut Option<ProjCache>,
        op: &dyn RealLinearOperator,
        capacity: usize,
        bands: usize,
    ) -> Result<Option<&'a mut ProjCache>, SolverError> {
        let count = op.projector_count();
        if count == 0 {
            return Ok(None);
        }
        if slot
            .as_ref()
            .is_none_or(|c| c.count != count || c.pv.cols() != capacity || c.pc.cols() != bands)
        {
            *slot = Some(ProjCache {
                count,
                pv: RealMatrix::zeros(count, capacity)?,
                pc: RealMatrix::zeros(count, bands)?,
                cols: 0,
                valid: false,
            });
        }
        let cache = slot.as_mut().unwrap();
        cache.cols = 0;
        cache.valid = false;
        Ok(Some(cache))
    }
    /// Basis columns `first..first + k` were just applied by `op` (a block of
    /// exactly those columns). `first == 0` restarts the cache.
    fn load(&mut self, op: &dyn RealLinearOperator, first: usize, k: usize) {
        if first == 0 {
            self.cols = 0;
            self.valid = true;
        }
        let range = first * self.count..(first + k) * self.count;
        self.valid = self.valid
            && first == self.cols
            && range.end <= self.pv.as_slice().len()
            && op.copy_last_projections(k, &mut self.pv.as_mut_slice()[range]);
        self.cols = first + k;
        if self.valid {
            lean_add(LeanCount::ProjCacheColumns, k);
        }
    }
    fn covers(&self, m: usize) -> bool {
        self.valid && self.cols == m
    }
    /// pc = (P^T V)[:, :m] C with C m x bands (column-major, leading dimension m).
    fn rotate(
        &mut self,
        gemm: &mut RealGemmWorkspace,
        c: &[f64],
        m: usize,
        bands: usize,
    ) -> Result<(), SolverError> {
        let count = self.count;
        flops_add(FlopSite::Small, count, m, bands);
        gemm.multiply_slices_into(
            GemmSite::Other,
            &self.pv.as_slice()[..count * m],
            count,
            m,
            count,
            RealOp::Identity,
            c,
            bands,
            m,
            self.pc.as_mut().subcols_mut(0, bands),
        )?;
        Ok(())
    }
    /// The basis became V C (contraction to `bands` columns).
    fn contract(
        &mut self,
        gemm: &mut RealGemmWorkspace,
        c: &[f64],
        m: usize,
        bands: usize,
    ) -> Result<(), SolverError> {
        if !self.covers(m) {
            self.valid = false;
            return Ok(());
        }
        self.rotate(gemm, c, m, bands)?;
        let len = self.count * bands;
        self.pv.as_mut_slice()[..len].copy_from_slice(&self.pc.as_slice()[..len]);
        self.cols = bands;
        Ok(())
    }
}
/// `ProjCache::load` on an optional cache.
fn proj_load(
    cache: &mut Option<&mut ProjCache>,
    op: &dyn RealLinearOperator,
    first: usize,
    k: usize,
) {
    if let Some(c) = cache.as_mut() {
        c.load(op, first, k);
    }
}

/// S[:, from..to] = V^T V[:, from..to] (one border GEMM) and its mirror.
fn overlap_border(
    gemm: &mut RealGemmWorkspace,
    basis: &RealMatrix,
    from: usize,
    to: usize,
    s: &mut RealMatrix,
) -> Result<(), SolverError> {
    if to <= from {
        return Ok(());
    }
    let n = basis.rows();
    flops_add(FlopSite::Overlap, to, n, to - from);
    gemm.multiply_slices_into(
        GemmSite::Gram,
        basis.as_slice(),
        n,
        to,
        n,
        RealOp::Transpose,
        &basis.as_slice()[from * n..],
        to - from,
        n,
        s.as_mut().submatrix_mut(0, from, to, to - from),
    )?;
    // The computed upper entries define the symmetric matrix.
    for j in from..to {
        for i in 0..j {
            let v = s[(i, j)];
            s[(j, i)] = v;
        }
    }
    Ok(())
}

/// LEAN_BORDER=fused (ticket PERCORE-KERNELS): with the basis at `from + k`
/// columns (E appended) and H E in `he`, one product V^T [E | H E] gives the
/// overlap border S[:, from..from + k] (written with its mirror, as
/// `overlap_border`) and the reduced border V^T (H E), kept in `lean.fb`
/// columns k..2k for `reduced_from_fused`. The basis is streamed once instead
/// of twice. On the aarch64 pure kernel every entry is the same fused
/// multiply-add chain as in the two separate products (bit-identical); on MKL
/// and x86 faer the wider product is rounding class.
fn fused_border(
    gemm: &mut RealGemmWorkspace,
    basis: &RealMatrix,
    e: &RealMatrix,
    he: &RealMatrix,
    lean: &mut LeanScratch,
    from: usize,
) -> Result<(), SolverError> {
    let n = basis.rows();
    let k = e.cols();
    let to = from + k;
    if basis.cols() != to || he.cols() != k || e.rows() != n || he.rows() != n {
        return Err(LinalgError::Dimension.into());
    }
    if lean.pair.as_slice().len() < n * 2 * k {
        // Sized once per shape (2 x bands columns), outside the iteration.
        let cols = lean.sel_values.len().max(k);
        lean.pair = RealMatrix::zeros(n, 2 * cols)?;
    }
    lean.pair.reshape_within_storage(n, 2 * k)?;
    let (left, right) = lean.pair.as_mut_slice().split_at_mut(n * k);
    left.copy_from_slice(e.as_slice());
    right.copy_from_slice(he.as_slice());
    lean.fb.reshape_within_storage(to, 2 * k)?;
    flops_add(FlopSite::Overlap, to, n, k);
    flops_add(FlopSite::Reduced, to, n, k);
    gemm.multiply_slices_into(
        GemmSite::Gram,
        basis.as_slice(),
        n,
        to,
        n,
        RealOp::Transpose,
        lean.pair.as_slice(),
        2 * k,
        n,
        lean.fb.as_mut(),
    )?;
    for j in 0..k {
        for i in 0..to {
            lean.s[(i, from + j)] = lean.fb[(i, j)];
        }
    }
    for j in from..to {
        for i in 0..j {
            let v = lean.s[(i, j)];
            lean.s[(j, i)] = v;
        }
    }
    lean.pending = Some((from, to));
    Ok(())
}

/// The reduced-matrix step of `reduced_matrix_into` for the columns
/// `prefix..m` whose border V^T (H E) `fused_border` already formed: the same
/// layout expansion, the border copied from `fb` columns k..2k, the same
/// symmetrization (bit-identical to the product it replaces on the kernel
/// where `fused_border` is).
fn reduced_from_fused(
    out: &mut RealMatrix,
    fb: &RealMatrix,
    prefix: usize,
    m: usize,
) -> Result<(), SolverError> {
    let k = m - prefix;
    if fb.rows() != m || fb.cols() != 2 * k {
        return Err(LinalgError::Dimension.into());
    }
    out.reshape_within_storage(m, m)?;
    for j in (0..prefix).rev() {
        out.as_mut_slice()
            .copy_within(j * prefix..(j + 1) * prefix, j * m);
    }
    for j in 0..k {
        out.as_mut_slice()[(prefix + j) * m..(prefix + j + 1) * m]
            .copy_from_slice(&fb.as_slice()[(k + j) * m..(k + j + 1) * m]);
    }
    for j in prefix..m {
        let (earlier, current) = out.as_mut_slice().split_at_mut(j * m);
        for (col, &entry) in earlier.chunks_exact_mut(m).zip(&current[..j]) {
            col[j] = entry;
        }
    }
    Ok(())
}

/// Extend T = L^{-1} from the leading `from` columns to `to` (bordered
/// Cholesky in panels). Returns (accepted, smallest pivot ratio of the new
/// columns); on rejection T is valid for the leading `from` columns only.
pub(crate) fn extend_factor(
    gemm: &mut RealGemmWorkspace,
    lean: &mut LeanScratch,
    from: usize,
    to: usize,
) -> Result<(bool, f64), SolverError> {
    let cap = lean.cap;
    let mut min_ratio = f64::INFINITY;
    let mut p0 = from;
    while p0 < to {
        let w = FACTOR_PANEL.min(to - p0);
        let (g, l) = (&mut lean.g[..w * w], &mut lean.l[..w * w]);
        for j in 0..w {
            for i in 0..w {
                g[i + j * w] = lean.s[(p0 + i, p0 + j)];
            }
        }
        if p0 > 0 {
            // W = T[:p0, :p0] S[:p0, p0..p0 + w], then G -= W^T W.
            lean.w.reshape_within_storage(p0, w)?;
            flops_add(FlopSite::Small, p0, p0, w);
            gemm.multiply_slices_into(
                GemmSite::Other,
                lean.t.as_slice(),
                p0,
                p0,
                cap,
                RealOp::Identity,
                &lean.s.as_slice()[p0 * cap..],
                w,
                cap,
                lean.w.as_mut(),
            )?;
            lean.ww.reshape_within_storage(w, w)?;
            flops_add(FlopSite::Small, w, p0, w);
            gemm.multiply_slices_into(
                GemmSite::Other,
                lean.w.as_slice(),
                p0,
                w,
                p0,
                RealOp::Transpose,
                lean.w.as_slice(),
                w,
                p0,
                lean.ww.as_mut(),
            )?;
            for (g, &x) in g.iter_mut().zip(lean.ww.as_slice()) {
                *g -= x;
            }
        }
        // Cholesky of the panel with the pivot-ratio test.
        for j in 0..w {
            let mut d = g[j + j * w];
            for k in 0..j {
                d -= l[j + k * w] * l[j + k * w];
            }
            let ratio = d / lean.s[(p0 + j, p0 + j)];
            if ratio.is_nan() {
                return Ok((false, 0.));
            }
            min_ratio = min_ratio.min(ratio);
            if ratio < PIVOT_RATIO_MIN {
                return Ok((false, min_ratio));
            }
            let ljj = d.sqrt();
            l[j + j * w] = ljj;
            for i in j + 1..w {
                let mut x = g[i + j * w];
                for k in 0..j {
                    x -= l[i + k * w] * l[j + k * w];
                }
                l[i + j * w] = x / ljj;
            }
        }
        // Panel inverse T22 = L22^{-1} by forward substitution.
        lean.tinv.reshape_within_storage(w, w)?;
        let tinv = lean.tinv.as_mut_slice();
        for c in 0..w {
            for i in 0..c {
                tinv[i + c * w] = 0.;
            }
            for i in c..w {
                let mut x = if i == c { 1. } else { 0. };
                for k in c..i {
                    x -= l[i + k * w] * tinv[k + c * w];
                }
                tinv[i + c * w] = x / l[i + i * w];
            }
        }
        for c in 0..w {
            for i in 0..p0 {
                lean.t[(i, p0 + c)] = 0.;
            }
            for i in 0..w {
                lean.t[(p0 + i, p0 + c)] = lean.tinv[(i, c)];
            }
        }
        if p0 > 0 {
            // T[p0.., :p0] = -T22 (W^T T[:p0, :p0]).
            lean.u.reshape_within_storage(w, p0)?;
            flops_add(FlopSite::Small, w, p0, p0);
            gemm.multiply_slices_into(
                GemmSite::Other,
                lean.w.as_slice(),
                p0,
                w,
                p0,
                RealOp::Transpose,
                lean.t.as_slice(),
                p0,
                cap,
                lean.u.as_mut(),
            )?;
            flops_add(FlopSite::Small, w, w, p0);
            gemm.multiply_slices_into(
                GemmSite::Other,
                lean.tinv.as_slice(),
                w,
                w,
                w,
                RealOp::Identity,
                lean.u.as_slice(),
                p0,
                w,
                lean.t.as_mut().submatrix_mut(p0, 0, w, p0),
            )?;
            for c in 0..p0 {
                for i in 0..w {
                    let v = lean.t[(p0 + i, c)];
                    lean.t[(p0 + i, c)] = -v;
                }
            }
        }
        if !lean.t.as_slice()[p0 * cap..(p0 + w) * cap]
            .iter()
            .all(|x| x.is_finite())
        {
            return Ok((false, 0.));
        }
        // Componentwise Gram rounding is amplified by the weighted row
        // sums of T: eps * (sum_j |T_ij| sqrt(S_jj))^2. A pivot alone
        // cannot bound this (a bidiagonal chain has uniformly good pivots).
        // Reject loss of half the working digits; this is an early guard,
        // not an orthogonality certificate (the actual return is checked).
        for i in p0..p0 + w {
            let growth: f64 = (0..=i)
                .map(|j| lean.t[(i, j)].abs() * lean.s[(j, j)].sqrt())
                .sum();
            if !growth.is_finite() || f64::EPSILON * growth * growth > f64::EPSILON.sqrt() {
                return Ok((false, min_ratio));
            }
        }
        p0 += w;
    }
    Ok((true, min_ratio))
}

/// Border of A = T H_r T^T for the columns prefix..m (H_r compact m x m); the
/// leading prefix block is kept bit for bit (B5 relies on it).
pub(crate) fn reduced_border(
    gemm: &mut RealGemmWorkspace,
    h: &RealMatrix,
    lean: &mut LeanScratch,
    prefix: usize,
    m: usize,
) -> Result<(), SolverError> {
    let cap = lean.cap;
    let a = &mut lean.a;
    a.reshape_within_storage(m, m)?;
    if prefix != m {
        for j in (0..prefix).rev() {
            a.as_mut_slice()
                .copy_within(j * prefix..(j + 1) * prefix, j * m);
        }
    }
    let c = m - prefix;
    if c == 0 {
        return Ok(());
    }
    // Z = T[prefix..m, :m]^T (m x c), P = H_r Z, A[:, prefix..] = T[:m, :m] P.
    lean.z.reshape_within_storage(m, c)?;
    for jj in 0..c {
        for i in 0..m {
            lean.z[(i, jj)] = lean.t[(prefix + jj, i)];
        }
    }
    lean.p.reshape_within_storage(m, c)?;
    flops_add(FlopSite::Small, m, m, 2 * c);
    gemm.multiply_slices_into(
        GemmSite::Other,
        h.as_slice(),
        m,
        m,
        m,
        RealOp::Identity,
        lean.z.as_slice(),
        c,
        m,
        lean.p.as_mut(),
    )?;
    gemm.multiply_slices_into(
        GemmSite::Other,
        lean.t.as_slice(),
        m,
        m,
        cap,
        RealOp::Identity,
        lean.p.as_slice(),
        c,
        m,
        a.as_mut().subcols_mut(prefix, c),
    )?;
    for j in prefix..m {
        let (earlier, current) = a.as_mut_slice().split_at_mut(j * m);
        for (col, &entry) in earlier.chunks_exact_mut(m).zip(&current[..j]) {
            col[j] = entry;
        }
    }
    Ok(())
}

/// lean.c (m x k) = T^T Y[:, cols] for the eigenvector columns `cols` of the
/// compact m x m `y` (`None`: the leading k columns).
pub(crate) fn ritz_coefficients(
    gemm: &mut RealGemmWorkspace,
    lean: &mut LeanScratch,
    y: &RealMatrix,
    m: usize,
    k: usize,
    cols: Option<&[usize]>,
) -> Result<(), SolverError> {
    let cap = lean.cap;
    let b = match cols {
        Some(cols) => {
            lean.y.reshape_within_storage(m, k)?;
            for (slot, &j) in cols.iter().enumerate() {
                lean.y.as_mut_slice()[slot * m..(slot + 1) * m]
                    .copy_from_slice(&y.as_slice()[j * m..(j + 1) * m]);
            }
            lean.y.as_slice()
        }
        None => &y.as_slice()[..m * k],
    };
    lean.c.reshape_within_storage(m, k)?;
    flops_add(FlopSite::Small, m, m, k);
    gemm.multiply_slices_into(
        GemmSite::Other,
        lean.t.as_slice(),
        m,
        m,
        cap,
        RealOp::Transpose,
        b,
        k,
        m,
        lean.c.as_mut(),
    )?;
    Ok(())
}

/// Gram defect of the returned Ritz block from its current, measured overlap:
/// X = V C implies X^T X = C^T S C. Only replicated inputs participate.
/// Fixed-order compensated reduced dot products reuse `p` and `z`; no GEMM dispatch, N-sized work,
/// allocation or communication. S must belong to exactly the basis of C.
///
/// The second result requests the direct Gram when cancellation could hide
/// a threshold crossing. Compensated sums remove the first-order m factor
/// of naive summation; the screening allowance is 4 eps (a + g^2),
/// where a = max |C|^T |S| |C|
/// and g = max_j sum_i |C_ij| sqrt(S_ii). Using |S| rather than a
/// matrix of norm bounds avoids treating an orthogonal basis as parallel.
/// This is a sensitivity screen, not a certificate for an arbitrarily stale S:
/// contractions/rebuilds must refresh S before calling this routine.
///
/// The distributed lean solver forms the same quantities column range by
/// column range on its ranks (`return_gram_growth`, `return_gram_columns`,
/// `return_gram_decision`; RETURN-GRAM-DISTRIBUTED): every column's
/// arithmetic is unchanged, so the maxima and the decision are bit-identical.
pub(crate) fn return_gram_defect(
    lean: &mut LeanScratch,
    m: usize,
    bands: usize,
    pool: &WorkerPool,
) -> Result<(f64, bool), SolverError> {
    let growth = return_gram_growth(lean, m, bands);
    if !growth.is_finite() {
        return Ok((f64::INFINITY, true));
    }
    let (defect, amplification) = return_gram_columns(lean, m, bands, 0, bands, pool)?;
    Ok(return_gram_decision(defect, amplification, growth))
}

/// g = max_j sum_i |C_ij| sqrt(S_ii) of `return_gram_defect`, in the original
/// order over i; infinite as soon as one column sum is not finite.
pub(crate) fn return_gram_growth(lean: &LeanScratch, m: usize, bands: usize) -> f64 {
    let cap = lean.cap;
    let s = lean.s.as_slice();
    let c = &lean.c.as_slice()[..m * bands];
    let mut growth = 0.0_f64;
    for j in 0..bands {
        let mut g = 0.0;
        for i in 0..m {
            g += c[i + j * m].abs() * s[i + i * cap].sqrt();
        }
        if !g.is_finite() {
            return f64::INFINITY;
        }
        growth = growth.max(g);
    }
    growth
}

/// The maxima of |C^T S C - I| (infinite on any non-finite entry) and of
/// |C|^T |S| |C| over the upper-triangle entries (i <= j) of the columns
/// `from..to`. Every column is formed by the same operations whichever range,
/// worker or rank forms it, so the maxima over a partition of 0..bands equal
/// those of one call over 0..bands bit for bit.
pub(crate) fn return_gram_columns(
    lean: &mut LeanScratch,
    m: usize,
    bands: usize,
    from: usize,
    to: usize,
    pool: &WorkerPool,
) -> Result<(f64, f64), SolverError> {
    lean.p.reshape_within_storage(m, bands)?;
    lean.z.reshape_within_storage(m, bands)?;
    let cap = lean.cap;
    let s = lean.s.as_slice();
    let c = &lean.c.as_slice()[..m * bands];
    let p = &mut lean.p.as_mut_slice()[from * m..to * m];
    let z = &mut lean.z.as_mut_slice()[from * m..to * m];
    pool.install(|| gram_sc_columns(s, cap, c, m, from, p, z));
    let p = &lean.p.as_slice()[..m * bands];
    let z = &lean.z.as_slice()[..m * bands];
    Ok(pool.install(|| gram_ctp_columns(c, p, z, m, from, to)))
}

/// The decision of `return_gram_defect` from the maxima (finite growth).
pub(crate) fn return_gram_decision(defect: f64, amplification: f64, growth: f64) -> (f64, bool) {
    if !defect.is_finite() {
        return (f64::INFINITY, true);
    }
    let allowance = 4.0 * f64::EPSILON * (amplification + growth * growth);
    (defect, defect + allowance >= 1e-12)
}

/// Column range `from..to` of rank `rank` out of `ranks` for the distributed
/// return Gram: contiguous ranges of nearly equal work (column j costs m^2
/// for S C and |S| |C| plus m (j + 1) for its triangle of C^T P), a function
/// of replicated sizes only.
pub(crate) fn return_gram_range(
    m: usize,
    bands: usize,
    rank: usize,
    ranks: usize,
) -> (usize, usize) {
    let ranks = ranks.max(1) as u128;
    let cost = |j: usize| (m as u128) * (m as u128 + j as u128 + 1);
    let total: u128 = (0..bands).map(cost).sum();
    // Boundary b(r): the first column whose preceding work reaches r/ranks of the total.
    let boundary = |r: u128| {
        if r == 0 {
            return 0;
        }
        if r >= ranks {
            return bands;
        }
        let target = total * r / ranks;
        let mut acc = 0u128;
        for j in 0..bands {
            if acc >= target {
                return j;
            }
            acc += cost(j);
        }
        bands
    };
    let r = rank as u128;
    (boundary(r), boundary(r + 1))
}

/// Columns `first..` of P = S C and Z = |S| |C| (`p`, `z` hold whole columns
/// of m rows). Entry (i, j) is the compensated sum over k of S(k, i) C(k, j),
/// k ascending, exactly as a scalar loop forms it. S is stored bitwise
/// symmetric (every writer mirrors), so S(k, i) is read as S(i, k) from the
/// contiguous column k; `GRAM_LANES` rows advance in lockstep, each with its
/// own unchanged operation sequence (vectorizable, bit-identical).
fn gram_sc_columns(
    s: &[f64],
    cap: usize,
    c: &[f64],
    m: usize,
    first: usize,
    p: &mut [f64],
    z: &mut [f64],
) {
    let count = p.len() / m.max(1);
    if count > 2 {
        let split = count / 2;
        let (pl, pr) = p.split_at_mut(split * m);
        let (zl, zr) = z.split_at_mut(split * m);
        WorkerPool::join(
            || gram_sc_columns(s, cap, c, m, first, pl, zl),
            || gram_sc_columns(s, cap, c, m, first + split, pr, zr),
        );
        return;
    }
    for (jj, (pj, zj)) in p.chunks_exact_mut(m).zip(z.chunks_exact_mut(m)).enumerate() {
        let j = first + jj;
        let cj = &c[j * m..(j + 1) * m];
        let mut i0 = 0;
        while i0 < m {
            let w = (m - i0).min(GRAM_LANES);
            let (pb, zb) = (&mut pj[i0..i0 + w], &mut zj[i0..i0 + w]);
            if w == GRAM_LANES {
                gram_lanes_s::<GRAM_LANES>(s, cap, cj, i0, pb, zb);
            } else {
                gram_lanes_s_tail(s, cap, cj, i0, pb, zb);
            }
            i0 += w;
        }
    }
}

/// Upper triangle (i <= j) of C^T P for the columns `from..to`: the maximum
/// of |C^T P - I| (infinite on any non-finite entry) and of |C|^T Z. Each
/// entry is the compensated sum over k, k ascending; eight columns i of C
/// share one column j of P and Z (independent chains).
fn gram_ctp_columns(
    c: &[f64],
    p: &[f64],
    z: &[f64],
    m: usize,
    from: usize,
    to: usize,
) -> (f64, f64) {
    if to - from > 2 {
        // Split the triangle into halves of equal work (j + 1 entries per column).
        let total = (to * (to + 1) - from * (from + 1)) / 2;
        let mut mid = from + 1;
        while mid < to - 1 && (mid * (mid + 1) - from * (from + 1)) / 2 < total / 2 {
            mid += 1;
        }
        let (a, b) = WorkerPool::join(
            || gram_ctp_columns(c, p, z, m, from, mid),
            || gram_ctp_columns(c, p, z, m, mid, to),
        );
        let defect = if a.0.is_finite() && b.0.is_finite() {
            a.0.max(b.0)
        } else {
            f64::INFINITY
        };
        return (defect, a.1.max(b.1));
    }
    let mut defect = 0.0_f64;
    let mut amplification = 0.0_f64;
    for j in from..to {
        let pj = &p[j * m..(j + 1) * m];
        let zj = &z[j * m..(j + 1) * m];
        let mut i0 = 0;
        while i0 <= j {
            let w = (j + 1 - i0).min(8);
            let mut gram = [0.0_f64; 8];
            let mut correction = [0.0_f64; 8];
            let mut absolute = [0.0_f64; 8];
            if w == 8 {
                let cols: [&[f64]; 8] = std::array::from_fn(|l| &c[(i0 + l) * m..(i0 + l + 1) * m]);
                for k in 0..m {
                    let (pk, zk) = (pj[k], zj[k]);
                    for l in 0..8 {
                        let ck = cols[l][k];
                        let adjusted = ck * pk - correction[l];
                        let next = gram[l] + adjusted;
                        correction[l] = (next - gram[l]) - adjusted;
                        gram[l] = next;
                        absolute[l] += ck.abs() * zk;
                    }
                }
            } else {
                for l in 0..w {
                    let ci = &c[(i0 + l) * m..(i0 + l + 1) * m];
                    for ((&ck, &pk), &zk) in ci.iter().zip(pj).zip(zj) {
                        let adjusted = ck * pk - correction[l];
                        let next = gram[l] + adjusted;
                        correction[l] = (next - gram[l]) - adjusted;
                        gram[l] = next;
                        absolute[l] += ck.abs() * zk;
                    }
                }
            }
            for l in 0..w {
                let i = i0 + l;
                amplification = amplification.max(absolute[l]);
                let error = (gram[l] - if i == j { 1.0 } else { 0.0 }).abs();
                if !error.is_finite() {
                    return (f64::INFINITY, amplification);
                }
                defect = defect.max(error);
            }
            i0 += w;
        }
    }
    (defect, amplification)
}

/// Rows per lockstep block of `return_gram_defect`'s S C pass.
const GRAM_LANES: usize = 16;

/// `L` compensated sums sum_k S(i0 + l, k) C(k, j) in lockstep (see
/// `return_gram_defect`); each lane's operations are those of one scalar sum.
#[inline(always)]
fn gram_lanes_s<const L: usize>(
    s: &[f64],
    cap: usize,
    cj: &[f64],
    i0: usize,
    p_out: &mut [f64],
    z_out: &mut [f64],
) {
    let mut sc = [0.0_f64; L];
    let mut correction = [0.0_f64; L];
    let mut absolute = [0.0_f64; L];
    for (k, &ck) in cj.iter().enumerate() {
        let row: &[f64; L] = s[k * cap + i0..k * cap + i0 + L]
            .try_into()
            .expect("lane block inside S");
        for l in 0..L {
            let product = row[l] * ck;
            let adjusted = product - correction[l];
            let next = sc[l] + adjusted;
            correction[l] = (next - sc[l]) - adjusted;
            sc[l] = next;
            absolute[l] += product.abs();
        }
    }
    p_out.copy_from_slice(&sc);
    z_out.copy_from_slice(&absolute);
}

/// Tail of `gram_lanes_s` (fewer than `GRAM_LANES` rows), same operations.
fn gram_lanes_s_tail(
    s: &[f64],
    cap: usize,
    cj: &[f64],
    i0: usize,
    p_out: &mut [f64],
    z_out: &mut [f64],
) {
    let w = p_out.len();
    let mut sc = [0.0_f64; GRAM_LANES];
    let mut correction = [0.0_f64; GRAM_LANES];
    let mut absolute = [0.0_f64; GRAM_LANES];
    for (k, &ck) in cj.iter().enumerate() {
        let row = &s[k * cap + i0..k * cap + i0 + w];
        for l in 0..w {
            let product = row[l] * ck;
            let adjusted = product - correction[l];
            let next = sc[l] + adjusted;
            correction[l] = (next - sc[l]) - adjusted;
            sc[l] = next;
            absolute[l] += product.abs();
        }
    }
    p_out.copy_from_slice(&sc[..w]);
    z_out.copy_from_slice(&absolute[..w]);
}

/// Fallback stage 1 for an append whose pivot test failed (basis at m + k
/// columns, unit-norm block E in `extra`, overlap border computed): one
/// explicit projection E' = E - V T^T T (V^T E) against the basis through the
/// factor (V^T E is the border already formed), columns with
/// |E'_j| <= DEPENDENCE_TOLERANCE dropped (the classic rank floor on a unit
/// input), the rest through the classic in-block Gram-Schmidt (with its rank
/// floor), appended and tested again. Returns the appended column count, or
/// `None` with the basis back at m columns and the surviving columns of E' in
/// `corrections` (spanning the same space modulo V as the raw corrections).
#[allow(clippy::too_many_arguments)]
fn reproject_block(
    gemm: &mut RealGemmWorkspace,
    projection: &mut Projection,
    basis: &mut RealMatrix,
    extra: &mut RealMatrix,
    corrections: &mut RealMatrix,
    scratch: &mut RealMatrix,
    lean: &mut LeanScratch,
    m: usize,
    k: usize,
) -> Result<Option<usize>, SolverError> {
    let n = basis.rows();
    let cap = lean.cap;
    basis.reshape_within_storage(n, m)?;
    // G = T^T (T S[:m, m..m+k]), m x k.
    lean.z.reshape_within_storage(m, k)?;
    lean.p.reshape_within_storage(m, k)?;
    flops_add(FlopSite::Small, m, m, 2 * k);
    gemm.multiply_slices_into(
        GemmSite::Other,
        lean.t.as_slice(),
        m,
        m,
        cap,
        RealOp::Identity,
        &lean.s.as_slice()[m * cap..],
        k,
        cap,
        lean.z.as_mut(),
    )?;
    gemm.multiply_slices_into(
        GemmSite::Other,
        lean.t.as_slice(),
        m,
        m,
        cap,
        RealOp::Transpose,
        lean.z.as_slice(),
        k,
        m,
        lean.p.as_mut(),
    )?;
    // V G, n x k (per-panel dgemm at any width on MKL builds).
    scratch.reshape_within_storage(n, k)?;
    flops_add(FlopSite::Projection, n, m, k);
    gemm.multiply_slices_wide_into(
        GemmSite::ProjectionApply,
        basis.as_slice(),
        n,
        m,
        n,
        RealOp::Identity,
        lean.p.as_slice(),
        k,
        m,
        scratch.as_mut(),
    )?;
    corrections.reshape_within_storage(n, k)?;
    let mut kept = 0;
    for j in 0..k {
        let dst = &mut corrections.as_mut_slice()[kept * n..(kept + 1) * n];
        for ((d, &e), &p) in dst.iter_mut().zip(column(extra, j)).zip(column(scratch, j)) {
            *d = e - p;
        }
        let len = norm(dst);
        finite(&[len])?;
        if len > DEPENDENCE_TOLERANCE {
            kept += 1;
        }
    }
    corrections.reshape_within_storage(n, kept)?;
    if kept == 0 {
        return Ok(None);
    }
    projection.run_into(corrections, None, gemm, extra)?;
    let fresh = extra.cols();
    if fresh == 0 {
        return Ok(None);
    }
    append_into(basis, extra)?;
    overlap_border(gemm, basis, m, m + fresh, &mut lean.s)?;
    let (ok, ratio) = extend_factor(gemm, lean, m, m + fresh)?;
    lean_pivot(ratio);
    if ok {
        Ok(Some(fresh))
    } else {
        basis.reshape_within_storage(n, m)?;
        Ok(None)
    }
}

/// Overlap and factor of a fresh basis (start, contraction, fallback
/// rotation). If the factor is rejected the basis is orthonormalized by the
/// classic projection and H is applied again (counted as a rebuild).
#[allow(clippy::too_many_arguments)]
fn factor_fresh(
    op: &dyn RealLinearOperator,
    gemm: &mut RealGemmWorkspace,
    projection: &mut Projection,
    basis: &mut RealMatrix,
    products: &mut RealMatrix,
    scratch: &mut RealMatrix,
    lean: &mut LeanScratch,
    applications: &mut usize,
    work: &mut crate::work_accounting::Trace,
    timers: Timers,
    proj: &mut Option<&mut ProjCache>,
) -> Result<bool, SolverError> {
    let m = basis.cols();
    lean.pending = None;
    let t = timers.start();
    overlap_border(gemm, basis, 0, m, &mut lean.s)?;
    timers.stop(t, Stage::Orthonormalize);
    let t = timers.start();
    let (ok, _) = extend_factor(gemm, lean, 0, m)?;
    timers.stop(t, Stage::Eigen);
    if ok {
        return Ok(false);
    }
    lean_add(LeanCount::Rebuilds, 1);
    // R1-b: the rebuilt basis is not followed by the history.
    lean.keep_invalidate();
    scratch.reshape_within_storage(basis.rows(), m)?;
    scratch.as_mut_slice().copy_from_slice(basis.as_slice());
    let t = timers.start();
    projection.run_into(scratch, None, gemm, basis)?;
    timers.stop(t, Stage::Orthonormalize);
    apply_into(op, basis, applications, work, products, timers)?;
    proj_load(proj, op, 0, basis.cols());
    lean_add(LeanCount::RebuildColumns, basis.cols());
    let m = basis.cols();
    overlap_border(gemm, basis, 0, m, &mut lean.s)?;
    if !extend_factor(gemm, lean, 0, m)?.0 {
        return Err(SolverError::Application(
            "lean Davidson: orthonormalized basis failed the overlap pivot test".into(),
        ));
    }
    Ok(true)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn solve(
    solver: &RealDavidson,
    op: &dyn RealLinearOperator,
    initial: Option<&RealMatrix>,
    diagonal: &[f64],
    policy: BandTolerances,
    capacity: usize,
    workspace: &mut RealDavidsonWorkspace,
    preconditioner: crate::DavidsonPreconditioner,
    options: ResearchOptions,
) -> Result<RealEigenSolution, SolverError> {
    solve_with_intent(
        solver,
        op,
        initial,
        diagonal,
        policy,
        capacity,
        workspace,
        preconditioner,
        options,
        SolveIntent::Converge,
    )
    .map(|a| a.solution)
}

/// The lean loop under an explicit intent (BREAKTHROUGH-3). `Converge` is the
/// ordinary solve bit for bit; `Advance(budget)` returns after the budget's
/// correction blocks or operator columns with the all-band Ritz state and its
/// true Ritz-change flags (the kept start is stored exactly as at a return).
#[allow(clippy::too_many_arguments)]
pub(super) fn solve_with_intent(
    solver: &RealDavidson,
    op: &dyn RealLinearOperator,
    initial: Option<&RealMatrix>,
    diagonal: &[f64],
    policy: BandTolerances,
    capacity: usize,
    workspace: &mut RealDavidsonWorkspace,
    preconditioner: crate::DavidsonPreconditioner,
    options: ResearchOptions,
    intent: SolveIntent,
) -> Result<AdvancedSubspace, SolverError> {
    let n = op.dimension();
    let bands = solver.bands;
    if workspace
        .lean
        .as_ref()
        .is_none_or(|l| l.cap != capacity || l.sel_values.len() != bands)
    {
        workspace.lean = Some(LeanScratch::new(capacity, bands)?);
    }
    let RealDavidsonWorkspace {
        gemm,
        projection,
        basis,
        products,
        reduced,
        eigensolver,
        values,
        coefficients,
        vectors,
        rotated_products,
        residual_block,
        residuals,
        corrections,
        extra,
        extra_products,
        random,
        selected,
        converged,
        previous,
        column_of,
        b5: b5_scratch,
        lean,
        kept,
        pool,
        lean_proj,
        returned_projections,
        ..
    } = workspace;
    let lean = lean.as_mut().unwrap();
    // R1-b: previous directions kept at a contraction (0 = the plain restart).
    let keep_r = options.restart_keep.directions(bands);
    if keep_r > 0 {
        lean.keep_prepare(bands);
    }
    // C3a: projector products of the basis columns (switch and operator support).
    let mut proj = if options.lean_proj_cache {
        ProjCache::activate(lean_proj, op, capacity, bands)?
    } else {
        None
    };
    *returned_projections = false;
    let pool = pool.clone();
    projection.narrow = options.narrow_proj;
    let timers = Timers::new(options.stage_timers);
    lean_add(LeanCount::Solves, 1);
    if options.b5 && b5_scratch.is_none() {
        *b5_scratch = Some(crate::bordered::Bordered::new(capacity));
    }
    let mut b5 = if options.b5 {
        b5_scratch.as_mut()
    } else {
        None
    };
    if let Some(b) = b5.as_mut() {
        b.invalidate();
    }
    previous.fill(f64::NAN);
    selected.clear();
    // REUSE-PSI (work order C2): an armed kept start is taken when `initial`
    // is the kept X (to the ORC packing tolerance); either way it is consumed.
    let from_kept = match kept.as_mut() {
        Some(k) if k.armed => {
            let take = k.valid
                && initial.is_some_and(|x| {
                    x.rows() == n
                        && x.cols() == bands
                        && x.as_slice()
                            .iter()
                            .zip(k.vectors.as_slice())
                            .all(|(a, b)| (a - b).abs() <= KEPT_START_TOLERANCE)
                });
            lean_add(
                if take {
                    LeanCount::KeptStarts
                } else {
                    LeanCount::KeptDeclined
                },
                1,
            );
            k.armed = false;
            k.valid = false;
            take
        }
        Some(k) => {
            k.valid = false;
            false
        }
        None => false,
    };
    // Start: identical to the classic path (orthonormal start basis), or the
    // kept X itself (S = X^T X is factored like any fresh basis).
    let mut rng = Random(solver.seed);
    if from_kept {
        let k = kept.as_ref().unwrap();
        basis.reshape_within_storage(n, bands)?;
        basis.as_mut_slice().copy_from_slice(k.vectors.as_slice());
        products.reshape_within_storage(n, bands)?;
        products
            .as_mut_slice()
            .copy_from_slice(k.products.as_slice());
        lean_add(LeanCount::KeptColumns, bands);
        // C3a: the kept P^T X belongs to the kept X (same operator).
        if let Some(c) = proj.as_mut() {
            let len = c.count * bands;
            c.cols = bands;
            c.valid = k.projection_count == c.count && k.projections.len() == len;
            if c.valid {
                c.pv.as_mut_slice()[..len].copy_from_slice(&k.projections);
                lean_add(LeanCount::ProjCacheKeptLoads, 1);
            }
        }
    } else if let Some(x) = initial {
        if x.rows() != n || x.cols() > capacity {
            return Err(LinalgError::Dimension.into());
        }
        projection.run_into(x, None, gemm, basis)?;
    } else {
        let starting_basis = free_electron_initial_guess(diagonal, bands)?;
        basis.reshape_within_storage(n, starting_basis.cols())?;
        basis
            .as_mut_slice()
            .copy_from_slice(starting_basis.as_slice());
    }
    while basis.cols() < bands {
        let extra = projection.run(&rng.block(n, bands - basis.cols())?, Some(basis), gemm)?;
        if extra.cols() == 0 {
            return Err(SolverError::InvalidParameter(
                "unable to complete initial basis",
            ));
        }
        append_into(basis, &extra)?;
    }
    // Ritz-change thresholds: the last mandatory root uses the buffer threshold.
    let threshold = |j: usize| {
        if j < policy.strict_bands - 1 {
            solver.tolerance
        } else {
            policy.buffer_tolerance
        }
    };
    let mut history = Vec::with_capacity(solver.max_iterations);
    let mut applications = 0;
    let mut work = crate::work_accounting::Trace::new(
        bands,
        policy.strict_bands,
        solver.tolerance,
        policy.buffer_tolerance,
    );
    let mut restarts = 0;
    let mut contractions = 0;
    let buffer = (0..bands).map(|j| j >= policy.strict_bands).collect();
    let mut result_values = vec![0.; bands];
    let mut result_vectors = RealMatrix::zeros(n, bands)?;
    let mut result_residuals = vec![0.; bands];
    let mut result_converged = vec![false; bands];
    residuals.fill(0.);
    if !from_kept {
        apply_into(op, basis, &mut applications, &mut work, products, timers)?;
        proj_load(&mut proj, op, 0, basis.cols());
    }
    // Columns whose H_r / A border is formed, and columns of S and T.
    let mut prefix = 0;
    let mut factored = 0;
    let mut max_basis_dimension = 0;
    // LEAN_BORDER=fused until the first failed pivot test of this solve; the
    // remaining appends take the separate path (bounds the discarded H
    // columns to one block per solve on systems where the test fails often).
    let mut fuse_allowed = options.lean_border_fused && !options.narrow_reduced;
    // BREAKTHROUGH-3: correction blocks appended in this call (the budget of
    // an `Advance` intent); `applications` counts its operator columns.
    let mut blocks_appended = 0usize;
    for iteration in 1..=solver.max_iterations {
        lean_add(LeanCount::Iterations, 1);
        let iteration_clock = timers.start();
        if factored < basis.cols() {
            factor_fresh(
                op,
                gemm,
                projection,
                basis,
                products,
                vectors,
                lean,
                &mut applications,
                &mut work,
                timers,
                &mut proj,
            )?;
            factored = basis.cols();
            prefix = 0;
        }
        let m = basis.cols();
        max_basis_dimension = max_basis_dimension.max(m);
        let t = timers.start();
        match lean.pending.take() {
            Some((from, to)) if from == prefix && to == m && prefix < m => {
                reduced_from_fused(reduced, &lean.fb, prefix, m)?
            }
            _ => reduced_matrix_into(
                basis,
                products,
                prefix,
                gemm,
                reduced,
                options.narrow_reduced,
            )?,
        }
        timers.stop(t, Stage::Reduced);
        let t = timers.start();
        reduced_border(gemm, reduced, lean, prefix, m)?;
        let old_prefix = prefix;
        prefix = m;
        coefficients.reshape_within_storage(m, m)?;
        let bordered = match b5.as_mut() {
            Some(b) if old_prefix + 1 == m && b.order == old_prefix && b.streak < MAX_STREAK => {
                let done = b.update(
                    gemm,
                    &lean.a,
                    &mut values[..m],
                    coefficients,
                    options.b5_parallel.then_some(&*pool),
                )?;
                if !done {
                    b5_add(B5Count::Declined);
                }
                done
            }
            _ => false,
        };
        if !bordered {
            eigensolver.solve_into(&lean.a, &mut values[..m], coefficients)?;
            eigen_full_add(m);
        }
        if let Some(b) = b5.as_mut() {
            if bordered {
                b5_add(B5Count::Bordered);
                b.streak += 1;
            } else {
                b5_add(B5Count::Full);
                b.streak = 0;
            }
            b.store(&values[..m], coefficients);
        }
        timers.stop(t, Stage::Eigen);
        for (j, flag) in converged.iter_mut().enumerate() {
            *flag = iteration > 1 && (values[j] - previous[j]).abs() < threshold(j);
        }
        previous.copy_from_slice(&values[..bands]);
        let finished = converged.iter().all(|&b| b)
            || (iteration == solver.max_iterations
                && converged[..policy.strict_bands].iter().all(|&b| b));
        // An `Advance` budget ends the call after its blocks or columns; the
        // `Converge` intent never sets it (identical control flow).
        let budget_stop = match intent {
            SolveIntent::Advance(b) if !finished && blocks_appended >= b.correction_blocks => {
                Some(BudgetStop::CorrectionBlocks)
            }
            SolveIntent::Advance(b) if !finished && applications >= b.h_columns => {
                Some(BudgetStop::HColumns)
            }
            _ => None,
        };
        if finished || budget_stop.is_some() || iteration == solver.max_iterations {
            // All-band Ritz vectors, H images and residual norms (the standard
            // stop guard reads them), once.
            let t = timers.start();
            ritz_coefficients(gemm, lean, coefficients, m, bands, None)?;
            timers.stop(t, Stage::Eigen);
            vectors.reshape_within_storage(n, bands)?;
            rotated_products.reshape_within_storage(n, bands)?;
            residual_block.reshape_within_storage(n, bands)?;
            let hx_formed = ritz_and_residuals(
                FlopSite::Return,
                options.ritz_fused,
                gemm,
                basis,
                products,
                &lean.c.as_slice()[..m * bands],
                &values[..bands],
                vectors,
                rotated_products,
                residual_block,
                residuals,
                &pool,
                timers,
            )?;
            let (mut defect, direct) = return_gram_defect(lean, m, bands, &pool)?;
            if direct {
                defect = orthogonality_defect(vectors);
            }
            if defect > 1e-12 {
                lean_add(LeanCount::Rebuilds, 1);
                projection.run_into(vectors, None, gemm, basis)?;
                if basis.cols() < bands {
                    return Err(SolverError::Application(
                        "lean Davidson: return rebuild lost required bands".into(),
                    ));
                }
                apply_into(op, basis, &mut applications, &mut work, products, timers)?;
                proj_load(&mut proj, op, 0, bands);
                lean.keep_invalidate();
                (factored, prefix) = (0, 0);
                previous.fill(f64::NAN);
                if iteration == solver.max_iterations {
                    return Err(SolverError::Application(
                        "lean Davidson: return Gram defect exceeds 1e-12".into(),
                    ));
                }
                continue;
            }
            let maximum = residuals.iter().copied().fold(0., f64::max);
            history.push(maximum);
            work.residuals(iteration, residuals);
            timers.iteration(iteration_clock, 0);
            if !finished && budget_stop.is_none() {
                result_residuals.copy_from_slice(residuals);
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: maximum,
                    residuals: result_residuals,
                });
            }
            lean_add(LeanCount::Returns, 1);
            // C3a: P^T X = (P^T V) C for the caller's nonlocal energy.
            let proj_ready = match proj.as_mut() {
                Some(c) if c.covers(m) => {
                    let t = timers.start();
                    c.rotate(gemm, &lean.c.as_slice()[..m * bands], m, bands)?;
                    timers.stop(t, Stage::Rotate);
                    lean_add(LeanCount::ProjCacheReturns, 1);
                    true
                }
                Some(_) => {
                    lean_add(LeanCount::ProjCacheMisses, 1);
                    false
                }
                None => false,
            };
            *returned_projections = proj_ready;
            // REUSE-PSI: keep X and H X (the return rotation) for the next solve.
            if let Some(k) = kept.as_mut().filter(|k| k.store && hx_formed) {
                k.vectors.as_mut_slice().copy_from_slice(vectors.as_slice());
                k.products
                    .as_mut_slice()
                    .copy_from_slice(rotated_products.as_slice());
                k.valid = true;
                // C3a: and P^T X with them (allocated once per shape).
                k.projection_count = 0;
                if let Some(c) = proj.as_ref().filter(|_| proj_ready) {
                    let len = c.count * bands;
                    k.projections.resize(len, 0.);
                    k.projections.copy_from_slice(&c.pc.as_slice()[..len]);
                    k.projection_count = c.count;
                }
                lean_add(LeanCount::KeptStored, 1);
            }
            work.returned();
            result_values.copy_from_slice(&values[..bands]);
            result_vectors
                .as_mut_slice()
                .copy_from_slice(vectors.as_slice());
            result_residuals.copy_from_slice(residuals);
            result_converged.copy_from_slice(converged);
            let diagnostics = SubspaceDiagnostics {
                correction_blocks: blocks_appended,
                h_columns: applications,
                stop: budget_stop.unwrap_or(BudgetStop::Converged),
                kept_start: from_kept,
                unconverged_roots: converged.iter().filter(|&&c| !c).count(),
            };
            return Ok(AdvancedSubspace {
                solution: RealEigenSolution {
                    eigenvalues: result_values,
                    eigenvectors: result_vectors,
                    residuals: result_residuals,
                    converged: result_converged,
                    buffer,
                    residuals_certified: false,
                    iterations: iteration,
                    history,
                    operator_applications: applications,
                    restarts,
                    contracted_restarts: contractions,
                    certification_attempts: 0,
                    max_basis_dimension,
                },
                diagnostics,
            });
        }
        if iteration == 1 && m == n {
            // Complete space: the Ritz pairs are exact; the next iteration
            // confirms them without corrections (as the classic path).
            history.push(0.);
            timers.iteration(iteration_clock, 0);
            continue;
        }
        // Residuals of the unconverged roots only (QE's update).
        selected.clear();
        selected.extend(
            (0..bands)
                .filter(|&j| !converged[j])
                .take(bands.min(capacity - bands)),
        );
        let k = selected.len();
        if k == 0 {
            return Err(SolverError::Application(
                "lean Davidson: no root to correct below a full capacity".into(),
            ));
        }
        for (slot, &j) in selected.iter().enumerate() {
            column_of[j] = slot;
            lean.sel_values[slot] = values[j];
        }
        let t = timers.start();
        ritz_coefficients(gemm, lean, coefficients, m, k, Some(&selected[..]))?;
        timers.stop(t, Stage::Eigen);
        // Whether the basis becomes the Ritz basis of the bands in this iteration.
        let contracted = m + k > capacity;
        if keep_r > 0 && !contracted {
            // R1-b history: this iteration's selected Ritz coefficients.
            lean.keep_record(m, selected);
        }
        vectors.reshape_within_storage(n, k)?;
        rotated_products.reshape_within_storage(n, k)?;
        residual_block.reshape_within_storage(n, k)?;
        ritz_and_residuals(
            FlopSite::RitzLoop,
            options.ritz_fused || options.ritz_fused_loop,
            gemm,
            basis,
            products,
            &lean.c.as_slice()[..m * k],
            &lean.sel_values[..k],
            vectors,
            rotated_products,
            residual_block,
            &mut lean.sel_norms[..k],
            &pool,
            timers,
        )?;
        lean_add(LeanCount::ResidualColumns, k);
        for (slot, &j) in selected.iter().enumerate() {
            residuals[j] = lean.sel_norms[slot];
        }
        history.push(lean.sel_norms[..k].iter().copied().fold(0., f64::max));
        work.residuals(iteration, residuals);
        corrections.reshape_within_storage(n, k)?;
        let t = timers.start();
        precondition_block(
            &pool,
            preconditioner,
            corrections,
            selected,
            &column_of[..],
            residual_block,
            vectors,
            diagonal,
            &values[..bands],
        )?;
        timers.stop(t, Stage::Precondition);
        if contracted {
            // Contraction restart: X and H X of all bands (QE's last).
            restarts += 1;
            contractions += 1;
            work.restart();
            lean_add(LeanCount::Contractions, 1);
            let t = timers.start();
            ritz_coefficients(gemm, lean, coefficients, m, bands, None)?;
            // R1-b: kept direction coefficients appended to C (room is left
            // for this iteration's correction block).
            let (kept, dropped) = if keep_r > 0 {
                lean.keep_directions(m, bands, selected, keep_r.min(capacity - bands - k))?
            } else {
                (0, 0)
            };
            lean_add(LeanCount::KeepDropped, dropped);
            timers.stop(t, Stage::Eigen);
            let cols = bands + kept;
            vectors.reshape_within_storage(n, cols)?;
            rotated_products.reshape_within_storage(n, cols)?;
            let t = timers.start();
            flops_add(FlopSite::Contraction, n, m, 2 * cols);
            gemm.rotate_pair_slices(
                basis.as_slice(),
                products.as_slice(),
                m,
                &lean.c.as_slice()[..m * cols],
                cols,
                vectors.as_mut_slice(),
                rotated_products.as_mut_slice(),
            )?;
            if let Some(c) = proj.as_mut() {
                if kept == 0 {
                    c.contract(gemm, &lean.c.as_slice()[..m * bands], m, bands)?;
                } else {
                    // The cache holds bands columns of P^T X only.
                    c.valid = false;
                }
            }
            timers.stop(t, Stage::Rotate);
            std::mem::swap(basis, vectors);
            std::mem::swap(products, rotated_products);
            let mut plain = kept == 0;
            if !plain {
                // Factor X, then the kept directions (pivot rule); any failure
                // falls back to the plain restart from X alone.
                lean.pending = None;
                let t = timers.start();
                overlap_border(gemm, basis, 0, cols, &mut lean.s)?;
                timers.stop(t, Stage::Orthonormalize);
                let t = timers.start();
                let (x_ok, _) = extend_factor(gemm, lean, 0, bands)?;
                let p_ok = x_ok && {
                    let (ok, ratio) = extend_factor(gemm, lean, bands, cols)?;
                    lean_pivot(ratio);
                    ok
                };
                timers.stop(t, Stage::Eigen);
                if p_ok {
                    lean_add(LeanCount::KeepContractions, 1);
                    lean_add(LeanCount::KeptDirections, kept);
                } else {
                    lean_add(LeanCount::KeepFallbacks, 1);
                    basis.reshape_within_storage(n, bands)?;
                    products.reshape_within_storage(n, bands)?;
                    plain = !x_ok;
                }
            }
            let rebuilt = plain
                && factor_fresh(
                    op,
                    gemm,
                    projection,
                    basis,
                    products,
                    vectors,
                    lean,
                    &mut applications,
                    &mut work,
                    timers,
                    &mut proj,
                )?;
            if keep_r > 0 && !rebuilt {
                // The history of the next iteration: the columns of X (a
                // rebuild has invalidated it).
                lean.keep_identity(bands, selected);
            }
            prefix = 0;
        }
        let m = basis.cols();
        // First appended column (the fallback may rotate or rebuild the basis).
        let mut first = m;
        // Unit-norm corrections, appended without projection.
        extra.reshape_within_storage(n, k)?;
        {
            let t = timers.start();
            let source: &RealMatrix = &*corrections;
            pool.install(|| {
                column_work(extra.as_mut_slice(), n, 0, &|j, dst| {
                    let len = norm(column(source, j));
                    finite(&[len])?;
                    for (dst, &src) in dst.iter_mut().zip(column(source, j)) {
                        *dst = if len == 0. { 0. } else { src / len };
                    }
                    Ok(())
                })
            })?;
            timers.stop(t, Stage::Orthonormalize);
        }
        // LEAN_BORDER=fused: H on E first, then both borders from one product.
        // A zero correction (exactly dependent) fails the pivot test anyway and
        // must not reach the operator: such a block takes the separate path.
        // Shape rule (traffic, not tuning): the fused product saves one stream
        // of the basis (n m) and costs copies of E and H E (about 4 n k), so it
        // is taken for blocks with 4 k <= m; wider blocks (the first append of
        // every solve) keep the separate products.
        let fused = fuse_allowed
            && 4 * k <= m
            && extra
                .as_slice()
                .chunks_exact(n)
                .all(|c| c.iter().any(|&x| x != 0.));
        if fused {
            apply_into(
                op,
                extra,
                &mut applications,
                &mut work,
                extra_products,
                timers,
            )?;
            let t = timers.start();
            append_into(basis, extra)?;
            fused_border(gemm, basis, extra, extra_products, lean, m)?;
            timers.stop(t, Stage::Orthonormalize);
        } else {
            let t = timers.start();
            append_into(basis, extra)?;
            overlap_border(gemm, basis, m, m + k, &mut lean.s)?;
            timers.stop(t, Stage::Orthonormalize);
        }
        let t = timers.start();
        let (ok, ratio) = extend_factor(gemm, lean, m, m + k)?;
        timers.stop(t, Stage::Eigen);
        lean_add(LeanCount::Appends, 1);
        lean_pivot(ratio);
        if fused {
            if ok {
                proj_load(&mut proj, op, m, k);
                append_into(products, extra_products)?;
                factored = basis.cols();
                lean_add(LeanCount::FusedBorders, 1);
                lean_add(LeanCount::AppendedColumns, k);
                blocks_appended += 1;
                timers.iteration(iteration_clock, k);
                continue;
            }
            // The images of E are discarded; the fallback applies H to its
            // own columns below (counted, expected zero).
            lean.pending = None;
            lean_add(LeanCount::FusedWasted, k);
            fuse_allowed = false;
        }
        let mut appended = ok;
        if !ok {
            // Stage 1: re-project only this block against V and test again.
            lean_add(LeanCount::Reprojections, 1);
            let t = timers.start();
            let kept = reproject_block(
                gemm,
                projection,
                basis,
                extra,
                corrections,
                rotated_products,
                lean,
                m,
                k,
            )?;
            timers.stop(t, Stage::Orthonormalize);
            if let Some(kept) = kept {
                lean_add(LeanCount::FallbackDropped, k - kept);
                appended = true;
            }
        }
        if !appended {
            // Stage 2: orthonormalize the basis once by rotating it to its
            // full Ritz basis, then project the block classically.
            lean_add(LeanCount::Fallbacks, 1);
            basis.reshape_within_storage(n, m)?;
            // After a contraction in this iteration the basis already is the
            // (orthonormal) Ritz basis of the bands and was factored afresh;
            // the decomposition in `coefficients` belongs to the old basis.
            if !contracted {
                let t = timers.start();
                ritz_coefficients(gemm, lean, coefficients, m, m, None)?;
                vectors.reshape_within_storage(n, m)?;
                rotated_products.reshape_within_storage(n, m)?;
                flops_add(FlopSite::Fallback, n, m, 2 * m);
                gemm.rotate_pair_slices(
                    basis.as_slice(),
                    products.as_slice(),
                    m,
                    &lean.c.as_slice()[..m * m],
                    m,
                    vectors.as_mut_slice(),
                    rotated_products.as_mut_slice(),
                )?;
                timers.stop(t, Stage::Orthonormalize);
                std::mem::swap(basis, vectors);
                std::mem::swap(products, rotated_products);
                // C3a: this m-column rotation is not followed (rare; measured
                // fallback rate 0): the cache stays invalid until a full apply.
                if let Some(c) = proj.as_mut() {
                    c.valid = false;
                }
                let rebuilt = factor_fresh(
                    op,
                    gemm,
                    projection,
                    basis,
                    products,
                    vectors,
                    lean,
                    &mut applications,
                    &mut work,
                    timers,
                    &mut proj,
                )?;
                if keep_r > 0 && !rebuilt {
                    // R1-b: the basis is the full Ritz basis; root j is column j.
                    lean.keep_identity(m, selected);
                }
            }
            prefix = 0;
            let m = basis.cols();
            first = m;
            let t = timers.start();
            projection.run_into(corrections, Some(basis), gemm, extra)?;
            if extra.cols() == 0 {
                random.reshape_within_storage(n, bands.min(capacity - m))?;
                for value in random.as_mut_slice() {
                    *value = rng.next();
                }
                projection.run_into(random, Some(basis), gemm, extra)?;
            }
            timers.stop(t, Stage::Orthonormalize);
            lean_add(LeanCount::FallbackDropped, k.saturating_sub(extra.cols()));
            if extra.cols() == 0 {
                timers.iteration(iteration_clock, 0);
                result_residuals.copy_from_slice(residuals);
                return Err(SolverError::NotConverged {
                    iterations: iteration,
                    max_residual: history.last().copied().unwrap_or(f64::NAN),
                    residuals: result_residuals,
                });
            }
            let k2 = extra.cols();
            append_into(basis, extra)?;
            overlap_border(gemm, basis, m, m + k2, &mut lean.s)?;
            if !extend_factor(gemm, lean, m, m + k2)?.0 {
                return Err(SolverError::Application(
                    "lean Davidson: projected corrections failed the overlap pivot test".into(),
                ));
            }
        }
        // H on the appended columns.
        let k = basis.cols() - first;
        let fresh: &RealMatrix = &*extra;
        debug_assert_eq!(fresh.cols(), k);
        apply_into(
            op,
            fresh,
            &mut applications,
            &mut work,
            extra_products,
            timers,
        )?;
        proj_load(&mut proj, op, first, k);
        append_into(products, extra_products)?;
        factored = basis.cols();
        lean_add(LeanCount::AppendedColumns, k);
        blocks_appended += 1;
        timers.iteration(iteration_clock, k);
    }
    unreachable!("validated positive iteration limit")
}

#[cfg(test)]
mod tests {
    use super::*;
    /// The scalar form of `return_gram_defect` (one compensated sum at a
    /// time, commit 7320cd8), kept as the bitwise reference of the lockstep form.
    fn return_gram_defect_reference(lean: &mut LeanScratch, m: usize, bands: usize) -> (f64, bool) {
        let mut growth = 0.0_f64;
        let mut p = vec![0.0; m * bands];
        let mut z = vec![0.0; m * bands];
        for j in 0..bands {
            let mut g = 0.0;
            for i in 0..m {
                g += lean.c[(i, j)].abs() * lean.s[(i, i)].sqrt();
                let (mut sc, mut correction, mut absolute) = (0.0, 0.0, 0.0);
                for k in 0..m {
                    let product = lean.s[(k, i)] * lean.c[(k, j)];
                    let adjusted = product - correction;
                    let next = sc + adjusted;
                    correction = (next - sc) - adjusted;
                    sc = next;
                    absolute += product.abs();
                }
                p[i + j * m] = sc;
                z[i + j * m] = absolute;
            }
            if !g.is_finite() {
                return (f64::INFINITY, true);
            }
            growth = growth.max(g);
        }
        let (mut defect, mut amplification) = (0.0_f64, 0.0_f64);
        for j in 0..bands {
            for i in 0..=j {
                let (mut gram, mut correction, mut absolute) = (0.0, 0.0, 0.0);
                for k in 0..m {
                    let c = lean.c[(k, i)];
                    let adjusted = c * p[k + j * m] - correction;
                    let next = gram + adjusted;
                    correction = (next - gram) - adjusted;
                    gram = next;
                    absolute += c.abs() * z[k + j * m];
                }
                amplification = amplification.max(absolute);
                let error = (gram - if i == j { 1.0 } else { 0.0 }).abs();
                if !error.is_finite() {
                    return (f64::INFINITY, true);
                }
                defect = defect.max(error);
            }
        }
        let allowance = 4.0 * f64::EPSILON * (amplification + growth * growth);
        (defect, defect + allowance >= 1e-12)
    }

    #[test]
    fn return_gram_lockstep_matches_scalar_reference_bitwise() {
        let mut state = 0x9e37_79b9_7f4a_7c15_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        for &(cap, m, bands) in &[
            (40, 37, 9),
            (64, 64, 17),
            (50, 16, 16),
            (30, 5, 3),
            (70, 49, 1),
        ] {
            let mut lean = LeanScratch::new(cap, bands).unwrap();
            // A near-identity symmetric S (mirrored as the solver stores it)
            // and C close to S^{-1/2}-like scaling: a realistic small defect.
            for j in 0..m {
                for i in 0..=j {
                    let v = if i == j {
                        1.0 + 1e-3 * next()
                    } else {
                        1e-3 * next()
                    };
                    lean.s[(i, j)] = v;
                    lean.s[(j, i)] = v;
                }
            }
            lean.c.reshape_within_storage(m, bands).unwrap();
            for j in 0..bands {
                for i in 0..m {
                    lean.c[(i, j)] = if i == j { 1.0 } else { 0.0 } + 1e-2 * next();
                }
            }
            let expected = return_gram_defect_reference(&mut lean, m, bands);
            for threads in [1, 3, 4] {
                let pool = WorkerPool::new(threads).unwrap();
                let got = return_gram_defect(&mut lean, m, bands, &pool).unwrap();
                assert_eq!(
                    got.0.to_bits(),
                    expected.0.to_bits(),
                    "cap {cap} m {m} bands {bands} T {threads}"
                );
                assert_eq!(got.1, expected.1);
            }
        }
    }

    /// RETURN-GRAM-DISTRIBUTED: the rank column ranges partition 0..bands and
    /// the maxima combined as the distributed solver combines them (disjoint
    /// reduction slots, zeros elsewhere) give the replicated defect and
    /// decision bit for bit, for every rank count, including a non-finite case.
    #[test]
    fn return_gram_column_ranges_reproduce_replicated_screen_bitwise() {
        let mut state = 0x2545_f491_4f6c_dd1d_u64;
        let mut next = || {
            state ^= state << 13;
            state ^= state >> 7;
            state ^= state << 17;
            (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
        };
        let pool = WorkerPool::new(2).unwrap();
        for &(cap, m, bands, poison) in &[
            (40, 37, 9, false),
            (64, 64, 17, false),
            (50, 16, 16, false),
            (30, 5, 3, false),
            (70, 49, 1, false),
            (40, 37, 9, true),
        ] {
            let mut lean = LeanScratch::new(cap, bands).unwrap();
            for j in 0..m {
                for i in 0..=j {
                    let v = if i == j {
                        1.0 + 1e-3 * next()
                    } else {
                        1e-3 * next()
                    };
                    lean.s[(i, j)] = v;
                    lean.s[(j, i)] = v;
                }
            }
            lean.c.reshape_within_storage(m, bands).unwrap();
            for j in 0..bands {
                for i in 0..m {
                    // Small perturbation: a defect around 1e-13, near the screen.
                    lean.c[(i, j)] = if i == j { 1.0 } else { 0.0 } + 1e-13 * next();
                }
            }
            if poison {
                lean.s[(3, 2)] = f64::NAN;
                lean.s[(2, 3)] = f64::NAN;
            }
            let expected = return_gram_defect(&mut lean, m, bands, &pool).unwrap();
            let growth = return_gram_growth(&lean, m, bands);
            for ranks in 1..=7 {
                let mut slots = vec![0.0_f64; 2 * ranks];
                let mut covered = 0;
                for rank in 0..ranks {
                    let (from, to) = return_gram_range(m, bands, rank, ranks);
                    assert_eq!(from, covered, "contiguous ranges");
                    covered = to;
                    let (d, a) = return_gram_columns(&mut lean, m, bands, from, to, &pool).unwrap();
                    slots[2 * rank] += d;
                    slots[2 * rank + 1] += a;
                }
                assert_eq!(covered, bands);
                let (mut dmax, mut amax) = (0.0_f64, 0.0_f64);
                for slot in slots.as_chunks::<2>().0 {
                    dmax = if dmax.is_finite() && slot[0].is_finite() {
                        dmax.max(slot[0])
                    } else {
                        f64::INFINITY
                    };
                    amax = amax.max(slot[1]);
                }
                let got = return_gram_decision(dmax, amax, growth);
                assert_eq!(
                    got.0.to_bits(),
                    expected.0.to_bits(),
                    "m {m} B {bands} np {ranks}"
                );
                assert_eq!(got.1, expected.1);
            }
        }
    }

    #[test]
    fn return_gram_uses_measured_overlap_and_screens_cancellation() {
        let serial = WorkerPool::new(1).unwrap();
        let mut lean = LeanScratch::new(2, 2).unwrap();
        lean.c.reshape_within_storage(2, 2).unwrap();
        lean.s[(0, 0)] = 4.0;
        lean.s[(1, 1)] = 9.0;
        lean.c[(0, 0)] = 0.5;
        lean.c[(1, 1)] = 1.0 / 3.0;
        let (defect, direct) = return_gram_defect(&mut lean, 2, 2, &serial).unwrap();
        assert!(defect <= 1e-15 && !direct);
        // Borderline cases retain the original direct 1e-12 check, rather
        // than trusting the rounded reduced Gram at the threshold.
        lean.c[(0, 0)] = 0.5 * (1.0 + 0.4995e-12);
        assert!(return_gram_defect(&mut lean, 2, 2, &serial).unwrap().1);
        lean.c[(0, 0)] = 0.5 * (1.0 + 2e-12);
        assert!(return_gram_defect(&mut lean, 2, 2, &serial).unwrap().0 > 1e-12);
        // A normalized difference of almost parallel vectors: the reduced
        // Gram can be close to one despite amplification of errors in S.
        let delta = 1e-6_f64;
        lean.s[(0, 0)] = 1.0;
        lean.s[(1, 1)] = 1.0;
        lean.s[(0, 1)] = 1.0 - delta;
        lean.s[(1, 0)] = 1.0 - delta;
        lean.c.reshape_within_storage(2, 1).unwrap();
        lean.c[(0, 0)] = (2.0 * delta).sqrt().recip();
        lean.c[(1, 0)] = -lean.c[(0, 0)];
        assert!(return_gram_defect(&mut lean, 2, 1, &serial).unwrap().1);
    }

    #[test]
    fn inverse_growth_rejects_chain_with_good_pivots() {
        let n = 8;
        let mut v = RealMatrix::zeros(n, n).unwrap();
        v[(0, 0)] = 1.0;
        for j in 1..n {
            v[(j - 1, j)] = 0.99_f64.sqrt();
            v[(j, j)] = 0.1;
        }
        let mut gemm = RealGemmWorkspace::with_pool(n, n, WorkerPool::new(1).unwrap()).unwrap();
        let mut lean = LeanScratch::new(n, n).unwrap();
        overlap_border(&mut gemm, &v, 0, n, &mut lean.s).unwrap();
        assert!(
            !extend_factor(&mut gemm, &mut lean, 0, n).unwrap().0,
            "every pivot passes, but the inverse amplifies overlap rounding"
        );
    }

    /// T = L^{-1} from panel-wise extension satisfies T S T^T = I for several
    /// split points of the same basis.
    #[test]
    fn factor_extension_inverts_overlap() {
        let (n, m, cap) = (200, 90, 96);
        let pool = WorkerPool::new(1).unwrap();
        let mut gemm = RealGemmWorkspace::with_pool(n, cap, pool).unwrap();
        let mut v = Random(5).block(n, m).unwrap();
        // Unit columns, mildly dependent.
        for j in 1..m {
            for i in 0..n {
                let prev = v[(i, j - 1)];
                v[(i, j)] += 0.6 * prev;
            }
        }
        for j in 0..m {
            let len = norm(column(&v, j));
            for x in column_mut(&mut v, j) {
                *x /= len;
            }
        }
        let mut basis = RealMatrix::zeros(n, cap).unwrap();
        basis.reshape_within_storage(n, m).unwrap();
        basis.as_mut_slice().copy_from_slice(v.as_slice());
        for splits in [vec![0, m], vec![0, 7, 70, 71, m], vec![0, 64, 65, m]] {
            let mut lean = LeanScratch::new(cap, 8).unwrap();
            for w in splits.windows(2) {
                overlap_border(&mut gemm, &basis, w[0], w[1], &mut lean.s).unwrap();
                let (ok, ratio) = extend_factor(&mut gemm, &mut lean, w[0], w[1]).unwrap();
                assert!(ok, "{splits:?}: ratio {ratio:e}");
            }
            let mut worst = 0_f64;
            for i in 0..m {
                for j in 0..m {
                    let mut x = 0.;
                    for a in 0..m {
                        for b in 0..m {
                            x += lean.t[(i, a)] * lean.s[(a, b)] * lean.t[(j, b)];
                        }
                    }
                    worst = worst.max((x - if i == j { 1. } else { 0. }).abs());
                }
            }
            assert!(worst < 1e-10, "{splits:?}: |T S T^T - I| = {worst:e}");
        }
    }
}
