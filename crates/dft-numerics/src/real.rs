//! Real column-major algebra for the isolated Gamma library path.
use crate::{linalg::LinalgError, pool::WorkerPool};
use faer::{Accum, MatMut, MatRef, Par, linalg::matmul::matmul};
use std::{
    ops::{Index, IndexMut},
    sync::Arc,
};

/// Borrowed real matrix; faer views preserve strides for submatrices.
pub type RealMatrixRef<'a> = MatRef<'a, f64>;
/// Exclusive borrowed real matrix, including strided submatrices.
pub type RealMatrixMut<'a> = MatMut<'a, f64>;

/// Owned column-major storage; entry (i,j) is at i+j*rows.
#[derive(Debug)]
pub struct RealMatrix {
    rows: usize,
    cols: usize,
    data: Vec<f64>,
}
fn length(rows: usize, cols: usize) -> Result<usize, LinalgError> {
    rows.checked_mul(cols)
        .filter(|&n| n <= isize::MAX as usize / size_of::<f64>())
        .ok_or(LinalgError::Dimension)
}
impl RealMatrix {
    /// Adopt storage after checking dimensions and addressable size.
    pub fn new(rows: usize, cols: usize, data: Vec<f64>) -> Result<Self, LinalgError> {
        if length(rows, cols)? != data.len() {
            return Err(LinalgError::Dimension);
        }
        Ok(Self { rows, cols, data })
    }
    /// Allocate a zero matrix with checked dimensions.
    pub fn zeros(rows: usize, cols: usize) -> Result<Self, LinalgError> {
        Self::new(rows, cols, vec![0.; length(rows, cols)?])
    }
    /// Change the active contiguous shape within initialized storage. No allocation,
    /// clearing or copying occurs; the caller must overwrite newly exposed entries.
    pub fn reshape_within_storage(&mut self, rows: usize, cols: usize) -> Result<(), LinalgError> {
        if length(rows, cols)? > self.data.len() {
            return Err(LinalgError::Dimension);
        }
        self.rows = rows;
        self.cols = cols;
        Ok(())
    }
    /// Reserve initialized backing storage before a solve starts, preserving entries.
    pub fn reserve_storage(&mut self, entries: usize) {
        if entries > self.data.len() {
            self.data.resize(entries, 0.);
        }
    }
    /// Row count.
    pub fn rows(&self) -> usize {
        self.rows
    }
    /// Column count.
    pub fn cols(&self) -> usize {
        self.cols
    }
    /// Contiguous column-major storage.
    pub fn as_slice(&self) -> &[f64] {
        &self.data[..self.rows * self.cols]
    }
    /// Whole backing storage regardless of the active shape (scratch tiles).
    #[allow(dead_code)]
    pub(crate) fn as_mut_slice_full(&mut self) -> &mut [f64] {
        &mut self.data
    }
    /// Mutable contiguous column-major storage.
    pub fn as_mut_slice(&mut self) -> &mut [f64] {
        &mut self.data[..self.rows * self.cols]
    }
    /// Borrow without copying or allocating.
    pub fn as_ref(&self) -> RealMatrixRef<'_> {
        MatRef::from_column_major_slice(self.as_slice(), self.rows, self.cols)
    }
    /// Borrow exclusively without copying or allocating.
    pub fn as_mut(&mut self) -> RealMatrixMut<'_> {
        MatMut::from_column_major_slice_mut(
            &mut self.data[..self.rows * self.cols],
            self.rows,
            self.cols,
        )
    }
}
impl Clone for RealMatrix {
    fn clone(&self) -> Self {
        Self {
            rows: self.rows,
            cols: self.cols,
            data: self.as_slice().to_vec(),
        }
    }
}
impl PartialEq for RealMatrix {
    fn eq(&self, other: &Self) -> bool {
        self.rows == other.rows && self.cols == other.cols && self.as_slice() == other.as_slice()
    }
}
impl Index<(usize, usize)> for RealMatrix {
    type Output = f64;
    fn index(&self, (i, j): (usize, usize)) -> &f64 {
        assert!(i < self.rows && j < self.cols);
        &self.data[i + j * self.rows]
    }
}
impl IndexMut<(usize, usize)> for RealMatrix {
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut f64 {
        assert!(i < self.rows && j < self.cols);
        &mut self.data[i + j * self.rows]
    }
}
fn finite(a: RealMatrixRef<'_>) -> Result<(), LinalgError> {
    if let Some(a) = a.try_as_col_major() {
        return if a
            .col_iter()
            .all(|col| col.as_slice().iter().all(|x| x.is_finite()))
        {
            Ok(())
        } else {
            Err(LinalgError::NonFinite)
        };
    }
    for j in 0..a.ncols() {
        for i in 0..a.nrows() {
            if !a[(i, j)].is_finite() {
                return Err(LinalgError::NonFinite);
            }
        }
    }
    Ok(())
}
/// Left operand operation; the right operand is untransposed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RealOp {
    /// A B, partitioned along output rows.
    Identity,
    /// A^T B, partitioned along the reduction dimension.
    Transpose,
}
const PANEL: usize = 256;
// Match the complex scheduling guard's 2,097,152 real flops: a real
// multiply-add costs two flops, hence 1,048,576 real multiply-adds.
// This is a scheduling heuristic, not a measured real-path speed claim.
const PARALLEL_WORK: usize = 1_048_576;

/// Call site of a real Davidson product (labels the optional per-site dgemm timers).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum GemmSite {
    Rotate,
    Gram,
    ProjectionOverlap,
    ProjectionApply,
    Reduced,
    Other,
}
#[cfg(feature = "mkl")]
impl From<GemmSite> for dft_mkl::DgemmSite {
    fn from(site: GemmSite) -> Self {
        match site {
            GemmSite::Rotate => Self::Rotate,
            GemmSite::Gram => Self::Gram,
            GemmSite::ProjectionOverlap => Self::ProjectionOverlap,
            GemmSite::ProjectionApply => Self::ProjectionApply,
            GemmSite::Reduced => Self::Reduced,
            GemmSite::Other => Self::Other,
        }
    }
}
/// Kernel used for the rotation products Q C (research switch; default Legacy).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum RotateKernel {
    /// Bounded 4 x 4 faer nano-GEMM calls (original, bit-preserving).
    #[default]
    Legacy,
    /// 8 x 4 register-tile kernel with one ascending-p chain per element.
    Blocked,
}
/// Pure-Rust product computed inside each fixed panel (`BAANDI_PURE_GEMM`,
/// P48B). The panel decomposition (256 output rows for A B, 256 reduction terms
/// for A^T B, summed in ascending panel order) is the same for every kernel and
/// depends only on the shapes, so every kernel gives results that do not depend
/// on the thread count. MKL builds with dgemm enabled keep dgemm for the sites
/// it covers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PureGemm {
    /// Bounded faer nano-GEMM calls (16 x 16 x 16 for A^T B, 4 x 4 x 256 for
    /// A B) inside each panel: the original arithmetic (default).
    #[default]
    Nano,
    /// One packed product per panel with the architecture's kernel: `Library`
    /// on x86_64, `Simd` elsewhere (see those variants for the reason).
    Packed,
    /// One faer `matmul` (`Par::Seq`) per panel. On x86_64 with AVX2 or
    /// AVX-512 faer runs it through its private x86 GEMM, which packs into a
    /// thread-local buffer created on a thread's first call and reused after
    /// that (no steady-state allocation). On other targets faer's backend
    /// allocates packing storage on every call, so this variant is used there
    /// only by tests and benchmarks.
    Library,
    /// Own register-tiled kernel through `pulp` (safe SIMD dispatch): every
    /// output entry is one ascending fused multiply-add chain over the whole
    /// reduction starting from +0, so its value does not depend on the tile,
    /// the reduction block, the panel or the vector width; packing storage is
    /// a fixed stack array (no allocation).
    Simd,
}
impl PureGemm {
    /// The kernel actually run for this variant on this target.
    pub fn resolve(self) -> Self {
        match self {
            Self::Packed if cfg!(target_arch = "x86_64") => Self::Library,
            Self::Packed => Self::Simd,
            other => other,
        }
    }
}
/// Reusable fixed-panel Gram scratch and shared executor. Construction allocates;
/// multiply_into neither allocates numerical storage nor creates worker threads.
/// The faer backend stays sequential within each panel, including MKL builds.
pub struct RealGemmWorkspace {
    long: usize,
    short: usize,
    tiles: Vec<RealMatrix>,
    /// Right operand packed in column blocks of four for the rotation kernel.
    pack: Vec<f64>,
    kernel: RotateKernel,
    /// Product inside each panel (`BAANDI_PURE_GEMM`).
    pure: PureGemm,
    /// MKL builds: rotation through dgemm unless `BAANDI_MKL_DGEMM=off` (A/B
    /// switch; adopted as default after the dev08 distribution check).
    #[cfg(feature = "mkl")]
    dgemm_rotation: bool,
    pool: Arc<WorkerPool>,
}
impl RealGemmWorkspace {
    /// Reserve up to `long` full rows and `short` subspace columns.
    pub fn with_pool(
        long: usize,
        short: usize,
        pool: Arc<WorkerPool>,
    ) -> Result<Self, LinalgError> {
        length(short, short)?;
        length(
            long.div_ceil(PANEL),
            short.checked_mul(short).ok_or(LinalgError::Dimension)?,
        )?;
        let tiles = (0..long.div_ceil(PANEL))
            .map(|_| RealMatrix::zeros(short, short))
            .collect::<Result<_, _>>()?;
        Ok(Self {
            long,
            short,
            tiles,
            pack: vec![0.; short.div_ceil(4) * 4 * short],
            kernel: match crate::switches::choice(crate::switches::Id::ROTATE_KERNEL) {
                "blocked" => RotateKernel::Blocked,
                _ => RotateKernel::Legacy,
            },
            pure: match crate::switches::choice(crate::switches::Id::PURE_GEMM) {
                "packed" => PureGemm::Packed,
                _ => PureGemm::Nano,
            },
            #[cfg(feature = "mkl")]
            dgemm_rotation: crate::switches::on(crate::switches::Id::MKL_DGEMM),
            pool,
        })
    }
    /// Slice form of `rotate_pair`: `v`, `hv` are m x k, `c` is k x n and the
    /// outputs m x n, all contiguous column-major. With the `mkl` feature each
    /// 256-row panel is one single-thread dgemm (panel grouping fixed, so the
    /// result does not depend on the thread budget); otherwise this is the
    /// ordinary `rotate_pair`.
    #[allow(clippy::too_many_arguments)]
    pub fn rotate_pair_slices(
        &mut self,
        v: &[f64],
        hv: &[f64],
        k: usize,
        c: &[f64],
        n: usize,
        out_v: &mut [f64],
        out_hv: &mut [f64],
    ) -> Result<(), LinalgError> {
        if k == 0 || n == 0 || !v.len().is_multiple_of(k) || c.len() != k * n {
            return Err(LinalgError::Dimension);
        }
        let m = v.len() / k;
        if hv.len() != m * k || out_v.len() != m * n || out_hv.len() != m * n {
            return Err(LinalgError::Dimension);
        }
        #[cfg(feature = "mkl")]
        {
            if self.dgemm_rotation
                && m <= self.long
                && k <= self.short
                && n <= self.short
                && m > 0
                && PANEL.min(m) * n <= self.short * self.short
            {
                // Non-finite inputs are caught by the output scan below: a NaN or
                // infinity in A, B propagates to every dependent entry of C.
                let count = m.div_ceil(PANEL);
                let parallel = count >= 4 && m.saturating_mul(n).saturating_mul(k) >= PARALLEL_WORK;
                let tiles_v = &mut self.tiles[..count];
                let mut ov = MatMut::from_column_major_slice_mut(out_v, m, n);
                let mut work = || {
                    dgemm_panels(
                        GemmSite::Rotate,
                        ov.as_mut(),
                        v,
                        m,
                        c,
                        k,
                        k,
                        tiles_v,
                        parallel,
                    )
                };
                if parallel {
                    self.pool.install(&mut work);
                } else {
                    work();
                }
                let tiles_h = &mut self.tiles[..count];
                let mut oh = MatMut::from_column_major_slice_mut(out_hv, m, n);
                let mut work = || {
                    dgemm_panels(
                        GemmSite::Rotate,
                        oh.as_mut(),
                        hv,
                        m,
                        c,
                        k,
                        k,
                        tiles_h,
                        parallel,
                    )
                };
                if parallel {
                    self.pool.install(&mut work);
                } else {
                    work();
                }
                return if out_v.iter().chain(out_hv.iter()).all(|x| x.is_finite()) {
                    Ok(())
                } else {
                    Err(LinalgError::NonFinite)
                };
            }
        }
        self.rotate_pair(
            MatRef::from_column_major_slice(v, m, k),
            MatRef::from_column_major_slice(hv, m, k),
            MatRef::from_column_major_slice(c, k, n),
            MatMut::from_column_major_slice_mut(out_v, m, n),
            MatMut::from_column_major_slice_mut(out_hv, m, n),
        )
    }
    /// Slice-input form of `multiply_into` (A, B column-major with leading
    /// dimensions; output any view). With the `mkl` feature each 256-row panel
    /// (Identity) or 256-term reduction panel (Transpose) is one sequential
    /// dgemm with the panel grouping fixed; otherwise it forwards to
    /// `multiply_into`.
    #[allow(clippy::too_many_arguments)]
    pub fn multiply_slices_into(
        &mut self,
        site: GemmSite,
        a: &[f64],
        a_rows: usize,
        a_cols: usize,
        lda: usize,
        op: RealOp,
        b: &[f64],
        b_cols: usize,
        ldb: usize,
        out: RealMatrixMut<'_>,
    ) -> Result<(), LinalgError> {
        self.multiply_slices_impl(site, false, a, a_rows, a_cols, lda, op, b, b_cols, ldb, out)
    }
    /// As `multiply_slices_into`, but a non-transposed product with a long output
    /// (the nonlocal back-application P W, rows of the output are the long
    /// dimension) is also taken through dgemm: only the per-panel scratch bound
    /// applies, not the whole-output bound.
    #[allow(clippy::too_many_arguments)]
    pub fn multiply_slices_wide_into(
        &mut self,
        site: GemmSite,
        a: &[f64],
        a_rows: usize,
        a_cols: usize,
        lda: usize,
        op: RealOp,
        b: &[f64],
        b_cols: usize,
        ldb: usize,
        out: RealMatrixMut<'_>,
    ) -> Result<(), LinalgError> {
        self.multiply_slices_impl(site, true, a, a_rows, a_cols, lda, op, b, b_cols, ldb, out)
    }
    #[allow(clippy::too_many_arguments)]
    fn multiply_slices_impl(
        &mut self,
        site: GemmSite,
        wide: bool,
        a: &[f64],
        a_rows: usize,
        a_cols: usize,
        lda: usize,
        op: RealOp,
        b: &[f64],
        b_cols: usize,
        ldb: usize,
        #[cfg_attr(not(feature = "mkl"), allow(unused_mut))] mut out: RealMatrixMut<'_>,
    ) -> Result<(), LinalgError> {
        #[cfg(not(feature = "mkl"))]
        let _ = wide;
        #[cfg(not(feature = "mkl"))]
        let _ = site;
        let fallback = |ws: &mut Self, out: RealMatrixMut<'_>| {
            ws.multiply_into(
                MatRef::from_column_major_slice_with_stride(a, a_rows, a_cols, lda),
                op,
                MatRef::from_column_major_slice_with_stride(
                    b,
                    if op == RealOp::Transpose {
                        a_rows
                    } else {
                        a_cols
                    },
                    b_cols,
                    ldb,
                ),
                out,
            )
        };
        #[cfg(feature = "mkl")]
        {
            let (m, k) = if op == RealOp::Transpose {
                (a_cols, a_rows)
            } else {
                (a_rows, a_cols)
            };
            let n = b_cols;
            let count = if op == RealOp::Transpose { k } else { m }.div_ceil(PANEL);
            if self.dgemm_rotation
                && m > 0
                && n > 0
                && k > 0
                && a_rows <= self.long
                && a_cols <= self.short
                && n <= self.short
                && ((wide && op == RealOp::Identity) || m * n <= self.short * self.short)
                && PANEL.min(m) * n <= self.short * self.short
                && lda >= a_rows
                && ldb >= k
                && out.nrows() == m
                && out.ncols() == n
            {
                // Inputs are not scanned: a non-finite input makes the output
                // non-finite (NaN and infinity propagate through every product and
                // sum), and the output scan below rejects it.
                let parallel = count >= 4 && m.saturating_mul(n).saturating_mul(k) >= PARALLEL_WORK;
                if op == RealOp::Transpose {
                    let tiles = &mut self.tiles[..count];
                    let mut work = || gram_dgemm(site, tiles, 0, a, lda, b, ldb, k, m, n, parallel);
                    if parallel {
                        self.pool.install(work);
                    } else {
                        work();
                    }
                    let mut sum = || assemble(out.as_mut(), tiles, 0, parallel);
                    if parallel {
                        self.pool.install(sum);
                    } else {
                        sum();
                    }
                } else {
                    let tiles = &mut self.tiles[..count];
                    let mut work =
                        || dgemm_panels(site, out.as_mut(), a, lda, b, ldb, k, tiles, parallel);
                    if parallel {
                        self.pool.install(work);
                    } else {
                        work();
                    }
                }
                return finite(out.as_ref());
            }
        }
        fallback(self, out)
    }
    pub fn set_rotate_kernel(&mut self, kernel: RotateKernel) {
        self.kernel = kernel;
    }
    /// Select the per-panel pure-Rust product (tests and A/B; the default
    /// comes from `BAANDI_PURE_GEMM`).
    pub fn set_pure_gemm(&mut self, pure: PureGemm) {
        self.pure = pure;
    }
    /// The per-panel pure-Rust product in use (as selected, not resolved).
    pub fn pure_gemm(&self) -> PureGemm {
        self.pure
    }
    /// RITZ-FUSED (algebra-stream item 2): Ritz vectors X = V C into `out_v` and
    /// residual vectors R = HV C - X diag(values) into `out_r`, without storing
    /// H X: each 256-row panel of HV C lives only in that panel's tile and is
    /// combined at once with the panel of X (still in cache). Every panel product
    /// is the call `rotate_pair_slices` makes for that panel and the residual is
    /// the same expression `h - lambda v`, so X and R are bit-identical to
    /// `rotate_pair_slices` followed by the elementwise residual (faer path, unit
    /// test; MKL: the same dgemm calls into the same tiles). Returns `Ok(false)`
    /// without touching the outputs when the research `Blocked` rotation kernel
    /// is selected or a 256-row panel of H X does not fit one tile (very few
    /// columns against a narrow workspace); the caller then takes the unfused
    /// path.
    #[allow(clippy::too_many_arguments)]
    pub fn rotate_residual_slices(
        &mut self,
        v: &[f64],
        hv: &[f64],
        k: usize,
        c: &[f64],
        n: usize,
        values: &[f64],
        out_v: &mut [f64],
        out_r: &mut [f64],
    ) -> Result<bool, LinalgError> {
        if k == 0 || n == 0 || !v.len().is_multiple_of(k) || c.len() != k * n || values.len() != n {
            return Err(LinalgError::Dimension);
        }
        let m = v.len() / k;
        if hv.len() != m * k
            || out_v.len() != m * n
            || out_r.len() != m * n
            || m > self.long
            || k > self.short
            || n > self.short
        {
            return Err(LinalgError::Dimension);
        }
        if m == 0 {
            return Ok(true);
        }
        let count = m.div_ceil(PANEL);
        let parallel = count >= 4 && m.saturating_mul(n).saturating_mul(k) >= PARALLEL_WORK;
        // A panel of H X must fit one tile (the bound of the MKL rotation).
        if PANEL.min(m) * n > self.short * self.short {
            return Ok(false);
        }
        #[cfg(feature = "mkl")]
        let mkl = self.dgemm_rotation;
        #[cfg(not(feature = "mkl"))]
        let mkl = false;
        if !mkl && self.kernel == RotateKernel::Blocked {
            return Ok(false);
        }
        let tiles = &mut self.tiles[..count];
        let job = FusedPanels {
            v,
            hv,
            m,
            k,
            c,
            values,
            mkl,
            pure: self.pure,
        };
        let ov = MatMut::from_column_major_slice_mut(out_v, m, n);
        let or = MatMut::from_column_major_slice_mut(out_r, m, n);
        if parallel {
            self.pool
                .install(|| fused_panels(&job, ov, or, 0, tiles, true));
        } else {
            fused_panels(&job, ov, or, 0, tiles, false);
        }
        if out_v.iter().chain(out_r.iter()).all(|x| x.is_finite()) {
            Ok(true)
        } else {
            Err(LinalgError::NonFinite)
        }
    }
    /// One side of `rotate_pair_slices`: out = A C with the kernel that call uses
    /// for each side (bit-identical to the corresponding output of the pair).
    /// Used to form H X late, only when a contraction restart needs it.
    pub fn rotate_slices(
        &mut self,
        a: &[f64],
        k: usize,
        c: &[f64],
        n: usize,
        out: &mut [f64],
    ) -> Result<(), LinalgError> {
        if k == 0 || n == 0 || !a.len().is_multiple_of(k) || c.len() != k * n {
            return Err(LinalgError::Dimension);
        }
        let m = a.len() / k;
        if out.len() != m * n || m > self.long || k > self.short || n > self.short {
            return Err(LinalgError::Dimension);
        }
        if m == 0 {
            return Ok(());
        }
        let count = m.div_ceil(PANEL);
        let parallel = count >= 4 && m.saturating_mul(n).saturating_mul(k) >= PARALLEL_WORK;
        #[cfg(feature = "mkl")]
        {
            if self.dgemm_rotation && PANEL.min(m) * n <= self.short * self.short {
                let tiles = &mut self.tiles[..count];
                let mut o = MatMut::from_column_major_slice_mut(out, m, n);
                let mut work =
                    || dgemm_panels(GemmSite::Rotate, o.as_mut(), a, m, c, k, k, tiles, parallel);
                if parallel {
                    self.pool.install(&mut work);
                } else {
                    work();
                }
                return finite(MatRef::from_column_major_slice(out, m, n));
            }
        }
        let av = MatRef::from_column_major_slice(a, m, k);
        let cv = MatRef::from_column_major_slice(c, k, n);
        let mut o = MatMut::from_column_major_slice_mut(out, m, n);
        let pure = self.pure;
        if self.kernel == RotateKernel::Blocked {
            pack_columns(&mut self.pack, cv);
            let pack = &self.pack[..n.div_ceil(4) * 4 * k];
            if parallel {
                self.pool.install(|| rotate(o.as_mut(), av, pack, k, true));
            } else {
                rotate(o.as_mut(), av, pack, k, false);
            }
        } else if parallel {
            self.pool
                .install(|| rotate_legacy(o.as_mut(), av, cv, pure, true));
        } else {
            rotate_legacy(o.as_mut(), av, cv, pure, false);
        }
        finite(o.as_ref())
    }
    /// NARROW (algebra-stream item 1): out (mq x k) = Q^T X for a narrow block X
    /// (n x k, contiguous; Q n x mq contiguous). Every entry is one eight-lane
    /// dot product in a fixed order, streaming each column of Q once from memory
    /// (the block stays in cache); columns of Q are processed in fixed pieces on
    /// the pool, so the result does not depend on the thread count. Rounding
    /// differs from the panel GEMM (rounding class).
    #[allow(clippy::too_many_arguments)]
    pub fn narrow_transpose_into(
        &self,
        q: &[f64],
        n: usize,
        mq: usize,
        x: &[f64],
        k: usize,
        mut out: RealMatrixMut<'_>,
    ) -> Result<(), LinalgError> {
        if q.len() < n * mq || x.len() < n * k || out.nrows() != mq || out.ncols() != k {
            return Err(LinalgError::Dimension);
        }
        if mq == 0 || k == 0 {
            return Ok(());
        }
        let parallel = mq > NARROW_PIECE && n.saturating_mul(mq).saturating_mul(k) >= PARALLEL_WORK;
        if parallel {
            self.pool
                .install(|| narrow_dots(q, n, x, k, out.as_mut(), 0, true));
        } else {
            narrow_dots(q, n, x, k, out.as_mut(), 0, false);
        }
        finite(out.as_ref())
    }
    /// NARROW: X -= Q S for a narrow S (mq x k). Every entry subtracts
    /// s_jc q_ij in ascending j, one rounded multiply and one rounded subtraction
    /// each (no fused multiply-add), so the arithmetic of an entry does not depend
    /// on the row blocks; rows go in fixed blocks of `NARROW_ROWS` on the pool.
    /// Q is read once per call; the block of X stays in cache.
    pub fn narrow_subtract(
        &self,
        q: &[f64],
        n: usize,
        mq: usize,
        s: RealMatrixRef<'_>,
        mut x: RealMatrixMut<'_>,
    ) -> Result<(), LinalgError> {
        if q.len() < n * mq || s.nrows() != mq || x.nrows() != n || x.ncols() != s.ncols() {
            return Err(LinalgError::Dimension);
        }
        if n == 0 || mq == 0 || s.ncols() == 0 {
            return Ok(());
        }
        let parallel =
            n > NARROW_ROWS && n.saturating_mul(mq).saturating_mul(s.ncols()) >= PARALLEL_WORK;
        if parallel {
            self.pool
                .install(|| narrow_update(q, n, mq, s, x.as_mut(), 0, true));
        } else {
            narrow_update(q, n, mq, s, x.as_mut(), 0, false);
        }
        finite(x.as_ref())
    }
    /// B5 parallel form (`BAANDI_B5=par`): out (rows x n, contiguous) = A B with
    /// the output columns in fixed groups of `COLUMN_GROUP` (a multiple of the
    /// four-column micro tile) computed in parallel on the pool. Each group's
    /// arithmetic is fixed (faer: the row panels of `rotate_legacy`; MKL: one
    /// sequential dgemm per group), so the result does not depend on the thread
    /// count.
    #[allow(clippy::too_many_arguments)]
    pub fn multiply_column_groups(
        &mut self,
        a: &[f64],
        rows: usize,
        k: usize,
        lda: usize,
        b: &[f64],
        n: usize,
        ldb: usize,
        out: &mut [f64],
    ) -> Result<(), LinalgError> {
        if out.len() != rows * n
            || lda < rows
            || ldb < k
            || (k > 0 && a.len() < (k - 1) * lda + rows)
            || (n > 0 && b.len() < (n - 1) * ldb + k)
        {
            return Err(LinalgError::Dimension);
        }
        if rows == 0 || n == 0 {
            return Ok(());
        }
        if k == 0 {
            out.fill(0.);
            return Ok(());
        }
        #[cfg(feature = "mkl")]
        let mkl = self.dgemm_rotation;
        #[cfg(not(feature = "mkl"))]
        let mkl = false;
        let job = GroupJob {
            a,
            rows,
            k,
            lda,
            b,
            ldb,
            mkl,
            pure: self.pure,
        };
        let parallel =
            n > COLUMN_GROUP && rows.saturating_mul(n).saturating_mul(k) >= PARALLEL_WORK;
        if parallel {
            self.pool.install(|| column_groups(&job, out, 0, true));
        } else {
            column_groups(&job, out, 0, false);
        }
        finite(MatRef::from_column_major_slice(out, rows, n))
    }
    /// Rotate V and HV by the same Ritz coefficients in one parallel region.
    /// Each side retains the original 256-row pieces and bounded kernel calls
    /// (24 independent pieces for a 2969-row pair). The guard is unchanged.
    pub fn rotate_pair(
        &mut self,
        v: RealMatrixRef<'_>,
        hv: RealMatrixRef<'_>,
        c: RealMatrixRef<'_>,
        mut out_v: RealMatrixMut<'_>,
        mut out_hv: RealMatrixMut<'_>,
    ) -> Result<(), LinalgError> {
        let (m, k, n) = (v.nrows(), v.ncols(), c.ncols());
        if hv.nrows() != m
            || hv.ncols() != k
            || c.nrows() != k
            || out_v.nrows() != m
            || out_hv.nrows() != m
            || out_v.ncols() != n
            || out_hv.ncols() != n
            || m > self.long
            || k > self.short
            || n > self.short
        {
            return Err(LinalgError::Dimension);
        }
        // Inputs are not scanned; non-finite values propagate to the outputs, which are checked.
        if m == 0 || n == 0 {
            return Ok(());
        }
        if k == 0 {
            out_v.fill(0.);
            out_hv.fill(0.);
            return Ok(());
        }
        let parallel =
            m.div_ceil(PANEL) >= 4 && m.saturating_mul(n).saturating_mul(k) >= PARALLEL_WORK;
        let pure = self.pure;
        if self.kernel == RotateKernel::Blocked {
            pack_columns(&mut self.pack, c);
            let pack = &self.pack[..n.div_ceil(4) * 4 * k];
            if parallel {
                self.pool.install(|| {
                    WorkerPool::join(
                        || rotate(out_v.as_mut(), v, pack, k, true),
                        || rotate(out_hv.as_mut(), hv, pack, k, true),
                    );
                });
            } else {
                rotate(out_v.as_mut(), v, pack, k, false);
                rotate(out_hv.as_mut(), hv, pack, k, false);
            }
        } else if parallel {
            self.pool.install(|| {
                WorkerPool::join(
                    || rotate_legacy(out_v.as_mut(), v, c, pure, true),
                    || rotate_legacy(out_hv.as_mut(), hv, c, pure, true),
                );
            });
        } else {
            rotate_legacy(out_v.as_mut(), v, c, pure, false);
            rotate_legacy(out_hv.as_mut(), hv, c, pure, false);
        }
        finite(out_v.as_ref())?;
        finite(out_hv.as_ref())
    }
    /// Overwrite output with A B or A^T B. Dimensions, capacity, and finiteness
    /// are checked. Zero reduction length produces zeros. Inputs may be views.
    pub fn multiply_into(
        &mut self,
        a: RealMatrixRef<'_>,
        op: RealOp,
        b: RealMatrixRef<'_>,
        mut out: RealMatrixMut<'_>,
    ) -> Result<(), LinalgError> {
        let (m, k) = if op == RealOp::Transpose {
            (a.ncols(), a.nrows())
        } else {
            (a.nrows(), a.ncols())
        };
        let n = b.ncols();
        if k != b.nrows()
            || out.nrows() != m
            || out.ncols() != n
            || a.nrows() > self.long
            || a.ncols() > self.short
            || n > self.short
        {
            return Err(LinalgError::Dimension);
        }
        // Inputs are not scanned; non-finite values propagate to the checked output.
        if m == 0 || n == 0 {
            return Ok(());
        }
        if k == 0 {
            out.fill(0.);
            return Ok(());
        }
        let count = if op == RealOp::Transpose { k } else { m }.div_ceil(PANEL);
        let parallel = count >= 4 && m.saturating_mul(n).saturating_mul(k) >= PARALLEL_WORK;
        let pure = self.pure;
        if op == RealOp::Transpose {
            let tiles = &mut self.tiles[..count];
            for tile in tiles.iter_mut() {
                // Restore the full scratch shape (a dgemm Gram reshapes tiles).
                tile.reshape_within_storage(self.short, self.short)?;
            }
            let mut work = || gram(tiles, 0, a, b, m, n, pure, parallel);
            if parallel {
                self.pool.install(work);
            } else {
                work();
            }
            let mut sum = || assemble(out.as_mut(), tiles, 0, parallel);
            if parallel {
                self.pool.install(sum);
            } else {
                sum();
            }
        } else {
            let mut work = || rotate_legacy(out.as_mut(), a, b, pure, parallel);
            if parallel {
                self.pool.install(work);
            } else {
                work();
            }
        }
        finite(out.as_ref())
    }
}
// One panel product out = a b (a may be a transposed view).
// Nano (default): general faer GEMM allocates opaque packing buffers on some
// targets. These call shapes stay within its allocation-free 4096-MAC nano-GEMM
// branch. Rotation uses the full reduction panel instead of repeated 16-term
// updates; Gram retains the shape that won its separate panel benchmark.
// Neither grouping depends on workers.
// Packed (P48B): the whole panel product in one packed kernel call (see
// `PureGemm`); the panel itself is the same, so thread independence holds.
fn bounded_matmul(
    mut out: RealMatrixMut<'_>,
    a: RealMatrixRef<'_>,
    b: RealMatrixRef<'_>,
    transpose: bool,
    pure: PureGemm,
) {
    match pure.resolve() {
        PureGemm::Library => {
            matmul(out, Accum::Replace, a, b, 1., Par::Seq);
            return;
        }
        PureGemm::Simd if simd_gemm(out.as_mut(), a, b) => return,
        _ => {}
    }
    let (rows, cols, reduction) = if transpose { (16, 16, 16) } else { (4, 4, 256) };
    for j in (0..out.ncols()).step_by(cols) {
        let n = cols.min(out.ncols() - j);
        for i in (0..out.nrows()).step_by(rows) {
            let m = rows.min(out.nrows() - i);
            for k in (0..a.ncols()).step_by(reduction) {
                let len = reduction.min(a.ncols() - k);
                matmul(
                    out.as_mut().submatrix_mut(i, j, m, n),
                    if k == 0 { Accum::Replace } else { Accum::Add },
                    a.submatrix(i, k, m, len),
                    b.submatrix(k, j, len, n),
                    1.,
                    Par::Seq,
                );
            }
        }
    }
}

// ---- P48B packed register-tile kernel (`PureGemm::Simd`) ----
//
// Blocking rule (stated, not tuned on one CPU). The register tile is MR x NR
// with NR = 4 output columns and MR = 8 rows: 4 x 4 accumulators of two lanes
// (NEON) plus four operand vectors and four broadcasts stay within the 32
// vector registers. The left operand is packed per MR-row strip and per
// reduction block of KC = 256 terms into a stack array of MR x KC doubles
// (16 KiB), which stays in any current L1d (32 to 128 KiB) while the strip is
// reused against every four-column group of the right operand; the right
// operand is read in place (contiguous columns, four sequential streams). The
// panel (256 rows or 256 reduction terms) is the caller's fixed decomposition.
// Every output entry is one ascending chain acc = fma(a(i,p), b(p,j), acc) from
// +0 over the whole reduction (a later reduction block reloads the stored
// partial, which is exact), so the value does not depend on MR, NR, KC, the
// panel, the thread count or the instruction set: every target with IEEE fused
// multiply-add gives the same bits. aarch64 uses NEON vector values (safe
// `core::arch` arithmetic, loads through `bytemuck`); other targets use the
// portable form of the same chain (a reference there: x86_64 runs `Library`).
const KC: usize = 256;
const NR: usize = 4;
const MR: usize = 8;
static ZERO_COLUMN: [f64; KC] = [0.; KC];

/// Panel product through the register-tile kernel. Returns false (nothing
/// written) when a column of `b` is not contiguous; the caller then uses the
/// nano path.
fn simd_gemm(mut out: RealMatrixMut<'_>, a: RealMatrixRef<'_>, b: RealMatrixRef<'_>) -> bool {
    // Same stride test as `ColRef::try_as_col_major` (shape shortcuts such as
    // "one row is always contiguous" are wrong for faer views).
    if b.row_stride() != 1 {
        return false;
    }
    let (m, n, k) = (out.nrows(), out.ncols(), a.ncols());
    if m == 0 || n == 0 {
        return true;
    }
    if k == 0 {
        out.fill(0.);
        return true;
    }
    let mut strip = [0_f64; MR * KC];
    for pc in (0..k).step_by(KC) {
        let kc = KC.min(k - pc);
        let strip = &mut strip[..MR * kc];
        for i in (0..m).step_by(MR) {
            let rows = MR.min(m - i);
            pack_strip(strip, a, i, rows, pc, kc);
            for j in (0..n).step_by(NR) {
                let cols = NR.min(n - j);
                let column = |c: usize| -> &[f64] {
                    if c < cols {
                        &b.col(j + c)
                            .try_as_col_major()
                            .expect("checked contiguous column")
                            .as_slice()[pc..pc + kc]
                    } else {
                        &ZERO_COLUMN[..kc]
                    }
                };
                let columns = [column(0), column(1), column(2), column(3)];
                let mut tile = [[0_f64; MR]; NR];
                if pc > 0 {
                    for (c, col) in tile.iter_mut().enumerate().take(cols) {
                        for (r, value) in col.iter_mut().enumerate().take(rows) {
                            *value = out[(i + r, j + c)];
                        }
                    }
                }
                let tile = register_tile(strip, columns, tile);
                for (c, col) in tile.iter().enumerate().take(cols) {
                    match out.as_mut().col_mut(j + c).try_as_col_major_mut() {
                        Some(dst) => dst.as_slice_mut()[i..i + rows].copy_from_slice(&col[..rows]),
                        None => {
                            for (r, &value) in col.iter().enumerate().take(rows) {
                                out[(i + r, j + c)] = value;
                            }
                        }
                    }
                }
            }
        }
    }
    true
}
/// Pack rows i..i+rows, reduction terms pc..pc+kc of `a` term-major
/// (strip[p * MR + r]), zero padding rows rows..MR.
fn pack_strip(
    strip: &mut [f64],
    a: RealMatrixRef<'_>,
    i: usize,
    rows: usize,
    pc: usize,
    kc: usize,
) {
    // Branch on the strides themselves (the exact test faer's
    // `try_as_col_major` / `try_as_row_major` apply), not on the shape: a
    // one-row transposed view has row stride != 1 and must not take the
    // column branch (P48B robustness-suite panic, syn-s1 and friends).
    if a.row_stride() == 1 {
        for (p, dst) in strip.as_chunks_mut::<MR>().0.iter_mut().enumerate() {
            let col = a.col(pc + p).try_as_col_major().expect("contiguous column");
            dst[..rows].copy_from_slice(&col.as_slice()[i..i + rows]);
            dst[rows..].fill(0.);
        }
    } else if a.col_stride() == 1 {
        // Transposed view (A^T of a column-major A): rows are contiguous.
        for r in 0..MR {
            if r < rows {
                let row = a.row(i + r).try_as_row_major().expect("contiguous row");
                for (dst, &x) in strip
                    .as_chunks_mut::<MR>()
                    .0
                    .iter_mut()
                    .zip(&row.as_slice()[pc..pc + kc])
                {
                    dst[r] = x;
                }
            } else {
                for dst in strip.as_chunks_mut::<MR>().0.iter_mut() {
                    dst[r] = 0.;
                }
            }
        }
    } else {
        for (p, dst) in strip.as_chunks_mut::<MR>().0.iter_mut().enumerate() {
            for (r, d) in dst.iter_mut().enumerate() {
                *d = if r < rows { a[(i + r, pc + p)] } else { 0. };
            }
        }
    }
}
/// tile[c][r] = fma chain over p of strip[p * MR + r] * columns[c][p],
/// continuing from the given partial values.
#[inline(always)]
fn register_tile(strip: &[f64], columns: [&[f64]; NR], tile: [[f64; MR]; NR]) -> [[f64; MR]; NR] {
    #[cfg(target_arch = "aarch64")]
    if let Some(simd) = pulp::aarch64::Neon::try_new() {
        return pulp::Simd::vectorize(simd, || neon_tile(simd, strip, columns, tile));
    }
    portable_tile(strip, columns, tile)
}
/// NEON form: accumulators are native two-lane vectors (`pulp`'s safe
/// wrappers of the `core::arch` intrinsics; value casts by `bytemuck`).
#[cfg(target_arch = "aarch64")]
#[inline(always)]
fn neon_tile(
    simd: pulp::aarch64::Neon,
    strip: &[f64],
    columns: [&[f64]; NR],
    tile: [[f64; MR]; NR],
) -> [[f64; MR]; NR] {
    use core::arch::aarch64::float64x2_t;
    use pulp::bytemuck::cast;
    const V: usize = MR / 2;
    let neon = simd.neon;
    let mut acc: [[float64x2_t; V]; NR] = tile.map(cast::<[f64; MR], [float64x2_t; V]>);
    let kc = strip.len() / MR;
    let [b0, b1, b2, b3] = columns;
    let (b0, b1, b2, b3) = (&b0[..kc], &b1[..kc], &b2[..kc], &b3[..kc]);
    for ((((ap, &x0), &x1), &x2), &x3) in strip
        .as_chunks::<MR>()
        .0
        .iter()
        .zip(b0)
        .zip(b1)
        .zip(b2)
        .zip(b3)
    {
        let av = cast::<[f64; MR], [float64x2_t; V]>(*ap);
        for v in 0..V {
            acc[0][v] = neon.vfmaq_n_f64(acc[0][v], av[v], x0);
            acc[1][v] = neon.vfmaq_n_f64(acc[1][v], av[v], x1);
            acc[2][v] = neon.vfmaq_n_f64(acc[2][v], av[v], x2);
            acc[3][v] = neon.vfmaq_n_f64(acc[3][v], av[v], x3);
        }
    }
    acc.map(cast::<[float64x2_t; V], [f64; MR]>)
}
/// Portable form of the same chain (targets without the NEON form).
#[inline(always)]
fn portable_tile(
    strip: &[f64],
    columns: [&[f64]; NR],
    mut tile: [[f64; MR]; NR],
) -> [[f64; MR]; NR] {
    for (p, ap) in strip.as_chunks::<MR>().0.iter().enumerate() {
        for (c, col) in tile.iter_mut().enumerate() {
            let x = columns[c][p];
            for (acc, &a) in col.iter_mut().zip(ap) {
                *acc = a.mul_add(x, *acc);
            }
        }
    }
    tile
}

#[allow(clippy::too_many_arguments)]
fn gram(
    tiles: &mut [RealMatrix],
    start: usize,
    a: RealMatrixRef<'_>,
    b: RealMatrixRef<'_>,
    m: usize,
    n: usize,
    pure: PureGemm,
    parallel: bool,
) {
    if parallel && tiles.len() > 1 {
        let half = tiles.len() / 2;
        let (left, right) = tiles.split_at_mut(half);
        WorkerPool::join(
            || gram(left, start, a, b, m, n, pure, true),
            || gram(right, start + half, a, b, m, n, pure, true),
        );
    } else {
        for (t, tile) in tiles.iter_mut().enumerate() {
            let row = (start + t) * PANEL;
            let len = PANEL.min(a.nrows() - row);
            bounded_matmul(
                tile.as_mut().submatrix_mut(0, 0, m, n),
                a.subrows(row, len).transpose(),
                b.subrows(row, len),
                true,
                pure,
            );
        }
    }
}
// Output columns are independent. Four-column pieces and ascending panel sums
// are fixed for every thread budget. Retain the initial +0 addition (including
// its signed-zero behavior) while overwriting, rather than clearing, output.
fn assemble(mut out: RealMatrixMut<'_>, tiles: &[RealMatrix], first: usize, parallel: bool) {
    let count = out.ncols();
    if parallel && count > 4 {
        let split = (count.div_ceil(4) / 2) * 4;
        let (left, right) = out.split_at_col_mut(split);
        WorkerPool::join(
            || assemble(left, tiles, first, true),
            || assemble(right, tiles, first + split, true),
        );
    } else {
        let m = out.nrows();
        for j in 0..count {
            for (t, tile) in tiles.iter().enumerate() {
                let src =
                    &tile.as_slice()[(first + j) * tile.rows()..(first + j) * tile.rows() + m];
                if let Some(col) = out.as_mut().col_mut(j).try_as_col_major_mut() {
                    for (dst, &value) in col.as_slice_mut().iter_mut().zip(src) {
                        *dst = if t == 0 { 0. + value } else { *dst + value };
                    }
                } else {
                    for (i, &value) in src.iter().enumerate() {
                        out[(i, j)] = if t == 0 {
                            0. + value
                        } else {
                            out[(i, j)] + value
                        };
                    }
                }
            }
        }
    }
}
fn rotate_legacy(
    mut out: RealMatrixMut<'_>,
    a: RealMatrixRef<'_>,
    b: RealMatrixRef<'_>,
    pure: PureGemm,
    parallel: bool,
) {
    let count = a.nrows().div_ceil(PANEL);
    if parallel && count > 1 {
        let split = (count / 2) * PANEL;
        let (left, right) = out.split_at_row_mut(split);
        WorkerPool::join(
            || rotate_legacy(left, a.subrows(0, split), b, pure, true),
            || rotate_legacy(right, a.subrows(split, a.nrows() - split), b, pure, true),
        );
    } else {
        for row in (0..a.nrows()).step_by(PANEL) {
            let len = PANEL.min(a.nrows() - row);
            bounded_matmul(
                out.as_mut().subrows_mut(row, len),
                a.subrows(row, len),
                b,
                false,
                pure,
            );
        }
    }
}
/// Inputs of the fused Ritz-and-residual panels (slices with leading dimension m).
struct FusedPanels<'a> {
    v: &'a [f64],
    hv: &'a [f64],
    m: usize,
    k: usize,
    c: &'a [f64],
    values: &'a [f64],
    #[cfg_attr(not(feature = "mkl"), allow(dead_code))]
    mkl: bool,
    pure: PureGemm,
}
/// Panels `first..` of the fused rotation; the split is the panel split of
/// `rotate_legacy` / `dgemm_panels`, the per-panel work is fixed.
fn fused_panels(
    job: &FusedPanels<'_>,
    out_v: RealMatrixMut<'_>,
    out_r: RealMatrixMut<'_>,
    first: usize,
    tiles: &mut [RealMatrix],
    parallel: bool,
) {
    let count = out_v.nrows().div_ceil(PANEL);
    if parallel && count > 1 {
        let half = count / 2;
        let split = half * PANEL;
        let (lv, rv) = out_v.split_at_row_mut(split);
        let (lr, rr) = out_r.split_at_row_mut(split);
        let (tl, tr) = tiles.split_at_mut(half);
        WorkerPool::join(
            || fused_panels(job, lv, lr, first, tl, true),
            || fused_panels(job, rv, rr, first + half, tr, true),
        );
        return;
    }
    let (mut out_v, mut out_r) = (out_v, out_r);
    let n = out_v.ncols();
    let cv = MatRef::from_column_major_slice(job.c, job.k, n);
    let v = MatRef::from_column_major_slice(job.v, job.m, job.k);
    let hv = MatRef::from_column_major_slice(job.hv, job.m, job.k);
    for (t, tile) in tiles.iter_mut().enumerate().take(count) {
        let local = t * PANEL;
        let row = (first + t) * PANEL;
        let len = PANEL.min(out_v.nrows() - local);
        let buffer = &mut tile.as_mut_slice_full()[..len * n];
        #[cfg(feature = "mkl")]
        if job.mkl {
            // The two dgemm calls of `rotate_pair_slices` for this panel, the
            // same tile, the same shapes.
            dft_mkl::dgemm_at(
                GemmSite::Rotate.into(),
                false,
                len,
                n,
                job.k,
                &job.v[row..],
                job.m,
                job.c,
                job.k,
                buffer,
            )
            .expect("checked dgemm shapes");
            for j in 0..n {
                out_v
                    .as_mut()
                    .col_mut(j)
                    .subrows_mut(local, len)
                    .try_as_col_major_mut()
                    .expect("contiguous column")
                    .as_slice_mut()
                    .copy_from_slice(&buffer[j * len..(j + 1) * len]);
            }
            dft_mkl::dgemm_at(
                GemmSite::Rotate.into(),
                false,
                len,
                n,
                job.k,
                &job.hv[row..],
                job.m,
                job.c,
                job.k,
                buffer,
            )
            .expect("checked dgemm shapes");
            combine_residual(&mut out_r, &out_v, buffer, local, len, job.values);
            continue;
        }
        bounded_matmul(
            out_v.as_mut().subrows_mut(local, len),
            v.subrows(row, len),
            cv,
            false,
            job.pure,
        );
        bounded_matmul(
            MatMut::from_column_major_slice_mut(buffer, len, n),
            hv.subrows(row, len),
            cv,
            false,
            job.pure,
        );
        combine_residual(&mut out_r, &out_v, buffer, local, len, job.values);
    }
}
/// r = h - lambda v for one panel (the residual expression of the solver).
fn combine_residual(
    out_r: &mut RealMatrixMut<'_>,
    out_v: &RealMatrixMut<'_>,
    h: &[f64],
    local: usize,
    len: usize,
    values: &[f64],
) {
    for (j, &lambda) in values.iter().enumerate() {
        let hc = &h[j * len..(j + 1) * len];
        let vc = out_v
            .as_ref()
            .col(j)
            .subrows(local, len)
            .try_as_col_major()
            .expect("contiguous column")
            .as_slice();
        let rc = out_r
            .as_mut()
            .col_mut(j)
            .subrows_mut(local, len)
            .try_as_col_major_mut()
            .expect("contiguous column")
            .as_slice_mut();
        for ((r, &hval), &v) in rc.iter_mut().zip(hc).zip(vc) {
            *r = hval - lambda * v;
        }
    }
}
/// Fixed column pieces of the narrow transpose product.
const NARROW_PIECE: usize = 16;
/// Fixed row blocks of the narrow update (8 KiB of one column).
const NARROW_ROWS: usize = 1024;
/// Output column group of `multiply_column_groups` (multiple of four).
const COLUMN_GROUP: usize = 32;
/// Eight independent lanes (element i adds into lane i % 8, ascending), combined
/// pairwise: a fixed order, so the value does not depend on the caller.
fn dot8(x: &[f64], y: &[f64]) -> f64 {
    const L: usize = 8;
    let mut s = [0.; L];
    let ((cx, tx), (cy, ty)) = (x.as_chunks::<L>(), y.as_chunks::<L>());
    for (a, b) in cx.iter().zip(cy) {
        for k in 0..L {
            s[k] += a[k] * b[k];
        }
    }
    for ((s, &a), &b) in s.iter_mut().zip(tx).zip(ty) {
        *s += a * b;
    }
    ((s[0] + s[4]) + (s[2] + s[6])) + ((s[1] + s[5]) + (s[3] + s[7]))
}
fn narrow_dots(
    q: &[f64],
    n: usize,
    x: &[f64],
    k: usize,
    out: RealMatrixMut<'_>,
    first: usize,
    parallel: bool,
) {
    let rows = out.nrows();
    if parallel && rows > NARROW_PIECE {
        let split = (rows.div_ceil(NARROW_PIECE) / 2) * NARROW_PIECE;
        let (top, bottom) = out.split_at_row_mut(split);
        WorkerPool::join(
            || narrow_dots(q, n, x, k, top, first, true),
            || narrow_dots(q, n, x, k, bottom, first + split, true),
        );
        return;
    }
    let mut out = out;
    for j in 0..rows {
        let qj = &q[(first + j) * n..(first + j + 1) * n];
        for c in 0..k {
            out[(j, c)] = dot8(qj, &x[c * n..(c + 1) * n]);
        }
    }
}
fn narrow_update(
    q: &[f64],
    n: usize,
    mq: usize,
    s: RealMatrixRef<'_>,
    x: RealMatrixMut<'_>,
    first: usize,
    parallel: bool,
) {
    let rows = x.nrows();
    if parallel && rows > NARROW_ROWS {
        let split = (rows.div_ceil(NARROW_ROWS) / 2) * NARROW_ROWS;
        let (top, bottom) = x.split_at_row_mut(split);
        WorkerPool::join(
            || narrow_update(q, n, mq, s, top, first, true),
            || narrow_update(q, n, mq, s, bottom, first + split, true),
        );
        return;
    }
    let mut x = x;
    for start in (0..rows).step_by(NARROW_ROWS) {
        let len = NARROW_ROWS.min(rows - start);
        let lo = first + start;
        let col = |j: usize| &q[j * n + lo..j * n + lo + len];
        let mut j = 0;
        while j + 4 <= mq {
            let (q0, q1, q2, q3) = (col(j), col(j + 1), col(j + 2), col(j + 3));
            for c in 0..s.ncols() {
                let (s0, s1, s2, s3) = (s[(j, c)], s[(j + 1, c)], s[(j + 2, c)], s[(j + 3, c)]);
                let xc = x
                    .as_mut()
                    .col_mut(c)
                    .subrows_mut(start, len)
                    .try_as_col_major_mut()
                    .expect("contiguous column")
                    .as_slice_mut();
                for i in 0..len {
                    let mut t = xc[i];
                    t -= s0 * q0[i];
                    t -= s1 * q1[i];
                    t -= s2 * q2[i];
                    t -= s3 * q3[i];
                    xc[i] = t;
                }
            }
            j += 4;
        }
        while j < mq {
            let qj = col(j);
            for c in 0..s.ncols() {
                let sj = s[(j, c)];
                let xc = x
                    .as_mut()
                    .col_mut(c)
                    .subrows_mut(start, len)
                    .try_as_col_major_mut()
                    .expect("contiguous column")
                    .as_slice_mut();
                for (xv, &qv) in xc.iter_mut().zip(qj) {
                    *xv -= sj * qv;
                }
            }
            j += 1;
        }
    }
}
struct GroupJob<'a> {
    a: &'a [f64],
    rows: usize,
    k: usize,
    lda: usize,
    b: &'a [f64],
    ldb: usize,
    #[cfg_attr(not(feature = "mkl"), allow(dead_code))]
    mkl: bool,
    pure: PureGemm,
}
/// Column groups `first..` (in units of columns) of A B into `out` (contiguous).
fn column_groups(job: &GroupJob<'_>, out: &mut [f64], first: usize, parallel: bool) {
    let n = out.len() / job.rows;
    if parallel && n > COLUMN_GROUP {
        let split = (n.div_ceil(COLUMN_GROUP) / 2) * COLUMN_GROUP;
        let (left, right) = out.split_at_mut(split * job.rows);
        WorkerPool::join(
            || column_groups(job, left, first, true),
            || column_groups(job, right, first + split, true),
        );
        return;
    }
    let a = MatRef::from_column_major_slice_with_stride(job.a, job.rows, job.k, job.lda);
    for g in (0..n).step_by(COLUMN_GROUP) {
        let w = COLUMN_GROUP.min(n - g);
        let col = first + g;
        let dst = &mut out[g * job.rows..(g + w) * job.rows];
        let b = &job.b[col * job.ldb..];
        #[cfg(feature = "mkl")]
        if job.mkl {
            dft_mkl::dgemm_at(
                GemmSite::Other.into(),
                false,
                job.rows,
                w,
                job.k,
                job.a,
                job.lda,
                b,
                job.ldb,
                dst,
            )
            .expect("checked dgemm shapes");
            continue;
        }
        let bv = MatRef::from_column_major_slice_with_stride(b, job.k, w, job.ldb);
        rotate_legacy(
            MatMut::from_column_major_slice_mut(dst, job.rows, w),
            a,
            bv,
            job.pure,
            false,
        );
    }
}
/// Row panels of C = A B, each one sequential dgemm into a private tile that is
/// then copied into place; panels are independent, their grouping is fixed.
#[cfg(feature = "mkl")]
#[allow(clippy::too_many_arguments)]
fn dgemm_panels(
    site: GemmSite,
    mut out: RealMatrixMut<'_>,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    k: usize,
    tiles: &mut [RealMatrix],
    parallel: bool,
) {
    let count = out.nrows().div_ceil(PANEL);
    if parallel && count > 1 {
        let half = count / 2;
        let split = half * PANEL;
        let (left, right) = out.split_at_row_mut(split);
        let (tl, tr) = tiles.split_at_mut(half);
        WorkerPool::join(
            || dgemm_panels(site, left, a, lda, b, ldb, k, tl, true),
            || dgemm_panels(site, right, &a[split..], lda, b, ldb, k, tr, true),
        );
    } else {
        let n = out.ncols();
        for (t, tile) in tiles.iter_mut().enumerate().take(count) {
            let row = t * PANEL;
            let len = PANEL.min(out.nrows() - row);
            let buffer = &mut tile.as_mut_slice_full()[..len * n];
            dft_mkl::dgemm_at(
                site.into(),
                false,
                len,
                n,
                k,
                &a[row..],
                lda,
                b,
                ldb,
                buffer,
            )
            .expect("checked dgemm shapes");
            // Column copies (contiguous when the output has unit row stride);
            // the same values as an element loop, without per-element checks.
            for j in 0..n {
                let src = &buffer[j * len..(j + 1) * len];
                match out
                    .as_mut()
                    .col_mut(j)
                    .subrows_mut(row, len)
                    .try_as_col_major_mut()
                {
                    Some(dst) => dst.as_slice_mut().copy_from_slice(src),
                    None => {
                        for (i, &v) in src.iter().enumerate() {
                            out[(row + i, j)] = v;
                        }
                    }
                }
            }
        }
    }
}
/// Reduction-panel tiles of A^T B, each one sequential dgemm straight into a
/// tile shaped m x n (so `assemble` reads it with leading dimension m).
#[cfg(feature = "mkl")]
#[allow(clippy::too_many_arguments)]
fn gram_dgemm(
    site: GemmSite,
    tiles: &mut [RealMatrix],
    start: usize,
    a: &[f64],
    lda: usize,
    b: &[f64],
    ldb: usize,
    k_total: usize,
    m: usize,
    n: usize,
    parallel: bool,
) {
    if parallel && tiles.len() > 1 {
        let half = tiles.len() / 2;
        let (left, right) = tiles.split_at_mut(half);
        WorkerPool::join(
            || gram_dgemm(site, left, start, a, lda, b, ldb, k_total, m, n, true),
            || {
                gram_dgemm(
                    site,
                    right,
                    start + half,
                    a,
                    lda,
                    b,
                    ldb,
                    k_total,
                    m,
                    n,
                    true,
                )
            },
        );
    } else {
        for (t, tile) in tiles.iter_mut().enumerate() {
            let row = (start + t) * PANEL;
            let len = PANEL.min(k_total - row);
            tile.reshape_within_storage(m, n).expect("tile capacity");
            dft_mkl::dgemm_at(
                site.into(),
                true,
                m,
                n,
                len,
                &a[row..],
                lda,
                &b[row..],
                ldb,
                tile.as_mut_slice(),
            )
            .expect("checked dgemm shapes");
        }
    }
}
/// Pack the k x n right operand into column blocks of four: block jb holds, for
/// each reduction index p, the four values b(p, 4 jb ..= 4 jb + 3) (zero padded).
fn pack_columns(pack: &mut [f64], b: RealMatrixRef<'_>) {
    let (k, n) = (b.nrows(), b.ncols());
    for jb in 0..n.div_ceil(4) {
        let block = &mut pack[jb * k * 4..(jb + 1) * k * 4];
        for c in 0..4 {
            let j = jb * 4 + c;
            if j < n {
                for p in 0..k {
                    block[p * 4 + c] = b[(p, j)];
                }
            } else {
                for p in 0..k {
                    block[p * 4 + c] = 0.;
                }
            }
        }
    }
}
/// Fused multiply-add where the target has it (single instruction, so the
/// kernel stays fast); otherwise separate multiply and add. The choice is a
/// compile-time property of the build, never of the thread budget.
#[inline(always)]
fn fma(a: f64, b: f64, acc: f64) -> f64 {
    if cfg!(any(target_feature = "fma", target_arch = "aarch64")) {
        a.mul_add(b, acc)
    } else {
        acc + a * b
    }
}
// Rotation C = A B with a fixed arithmetic order: every output element is one
// ascending-p chain acc += a(i,p) b(p,j) starting from +0, so it does not depend
// on the row partition, the thread budget or the tile it falls into. The 8 x 4
// register tile keeps 32 accumulators; A columns are contiguous.
fn rotate(
    mut out: RealMatrixMut<'_>,
    a: RealMatrixRef<'_>,
    pack: &[f64],
    k: usize,
    parallel: bool,
) {
    let count = a.nrows().div_ceil(PANEL);
    if parallel && count > 1 {
        let split = (count / 2) * PANEL;
        let (left, right) = out.split_at_row_mut(split);
        WorkerPool::join(
            || rotate(left, a.subrows(0, split), pack, k, true),
            || rotate(right, a.subrows(split, a.nrows() - split), pack, k, true),
        );
    } else {
        let rows = a.nrows();
        let n = out.ncols();
        for ib in (0..rows).step_by(8) {
            let rr = 8.min(rows - ib);
            for jb in 0..n.div_ceil(4) {
                let nc = 4.min(n - jb * 4);
                let block = &pack[jb * k * 4..(jb + 1) * k * 4];
                let mut acc = [[0_f64; 8]; 4];
                if rr == 8 {
                    for p in 0..k {
                        let column = a.col(p).try_as_col_major().expect("contiguous column");
                        let av: &[f64; 8] = column.as_slice()[ib..ib + 8].try_into().unwrap();
                        let bv: &[f64; 4] = block[p * 4..p * 4 + 4].try_into().unwrap();
                        for c in 0..4 {
                            for r in 0..8 {
                                acc[c][r] = fma(av[r], bv[c], acc[c][r]);
                            }
                        }
                    }
                } else {
                    for p in 0..k {
                        let column = a.col(p).try_as_col_major().expect("contiguous column");
                        let av = &column.as_slice()[ib..ib + rr];
                        let bv = &block[p * 4..p * 4 + 4];
                        for c in 0..4 {
                            for (r, &x) in av.iter().enumerate() {
                                acc[c][r] = fma(x, bv[c], acc[c][r]);
                            }
                        }
                    }
                }
                for c in 0..nc {
                    for r in 0..rr {
                        out[(ib + r, jb * 4 + c)] = 0. + acc[c][r];
                    }
                }
            }
        }
    }
}

/// Real symmetric eigensolve, ascending eigenvalues, sequential faer backend.
/// Largest-magnitude component of each eigenvector is positive (first index
/// wins ties). Degenerate eigenspaces are not assigned a unique basis.
/// Symmetry uses the existing adapter's 1e-12 relative entrywise allowance,
/// with unit absolute floor; faer consumes the lower triangle after validation.
pub fn symmetric_eigh(a: &RealMatrix) -> Result<(Vec<f64>, RealMatrix), LinalgError> {
    let mut ws = RealEighWorkspace::new(a.rows())?;
    let mut values = vec![0.; a.rows()];
    let mut vectors = RealMatrix::zeros(a.rows(), a.rows())?;
    ws.solve_into(a, &mut values, &mut vectors)?;
    Ok((values, vectors))
}
/// Caller-owned EVD storage. Every active size is covered at construction;
/// backend algorithm, strides and sequential arithmetic match symmetric_eigh.
pub struct RealEighWorkspace {
    capacity: usize,
    /// MKL `dsyevd` instead of faer when built with `mkl`, unless `BAANDI_MKL_EVD=off`
    /// (default on after the dev08 distribution evidence, ledger MKL-EVD).
    #[cfg(feature = "mkl")]
    mkl: Option<dft_mkl::SymmetricEigen>,
    diagonal: Vec<f64>,
    buffer: faer::dyn_stack::MemBuffer,
}
impl RealEighWorkspace {
    /// Reserve scratch for every possible active subspace size.
    pub fn new(capacity: usize) -> Result<Self, LinalgError> {
        length(capacity, capacity)?;
        let requirements: Vec<_> = (1..=capacity)
            .map(|n| {
                faer::linalg::evd::self_adjoint_evd_scratch::<f64>(
                    n,
                    faer::linalg::evd::ComputeEigenvectors::Yes,
                    Par::Seq,
                    Default::default(),
                )
            })
            .collect();
        Ok(Self {
            capacity,
            #[cfg(feature = "mkl")]
            mkl: crate::switches::on(crate::switches::Id::MKL_EVD)
                .then(|| dft_mkl::SymmetricEigen::new(capacity).ok())
                .flatten(),
            diagonal: vec![0.; capacity],
            buffer: faer::dyn_stack::MemBuffer::new(faer::dyn_stack::StackReq::any_of(
                &requirements,
            )),
        })
    }
    /// Reuse output and scratch storage. The faer backend can still allocate
    /// internal GEMM packing/bookkeeping; this is not a whole-call no-allocation
    /// guarantee (see the P33A counting-allocator diagnostic).
    pub fn solve_into(
        &mut self,
        a: &RealMatrix,
        values: &mut [f64],
        vectors: &mut RealMatrix,
    ) -> Result<(), LinalgError> {
        if a.rows != a.cols {
            return Err(LinalgError::NotSquare);
        }
        finite(a.as_ref())?;
        let allowed = crate::linalg::HERMITIAN_TOLERANCE
            * a.as_slice().iter().map(|v| v.abs()).fold(1., f64::max);
        let mut defect: f64 = 0.;
        for j in 0..a.cols {
            for i in 0..j {
                defect = defect.max((a[(i, j)] - a[(j, i)]).abs());
            }
        }
        if defect > allowed {
            return Err(LinalgError::NotHermitian { defect, allowed });
        }
        let n = a.rows;
        if n > self.capacity || values.len() != n || vectors.rows() != n || vectors.cols() != n {
            return Err(LinalgError::Dimension);
        }
        if n == 0 {
            return Ok(());
        }
        #[cfg(feature = "mkl")]
        let used_mkl = if let Some(mkl) = self.mkl.as_mut() {
            vectors.as_mut_slice().copy_from_slice(a.as_slice());
            mkl.solve(n, vectors.as_mut_slice(), values)
                .map_err(|_| LinalgError::Decomposition)?;
            true
        } else {
            false
        };
        #[cfg(not(feature = "mkl"))]
        let used_mkl = false;
        if !used_mkl {
            faer::linalg::evd::self_adjoint_evd(
                a.as_ref(),
                faer::diag::DiagMut::from_slice_mut(&mut self.diagonal[..n]),
                Some(vectors.as_mut()),
                Par::Seq,
                faer::dyn_stack::MemStack::new(&mut self.buffer),
                Default::default(),
            )
            .map_err(|_| LinalgError::Decomposition)?;
            values.copy_from_slice(&self.diagonal[..n]);
        }
        if values.iter().any(|v| !v.is_finite()) {
            return Err(LinalgError::NonFinite);
        }
        finite(vectors.as_ref())?;
        for j in 0..n {
            let mut pivot = 0;
            for i in 1..n {
                if vectors[(i, j)].abs() > vectors[(pivot, j)].abs() {
                    pivot = i;
                }
            }
            if vectors[(pivot, j)].is_sign_negative() {
                for i in 0..n {
                    vectors[(i, j)] = -vectors[(i, j)];
                }
            }
        }
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::linalg::{ComplexMatrix, GemmWorkspace, MatrixOp};
    use num_complex::Complex64;
    fn fixture(r: usize, c: usize) -> RealMatrix {
        RealMatrix::new(
            r,
            c,
            (0..r * c)
                .map(|i| ((i * 13 + 7) as f64).sin() / 32.)
                .collect(),
        )
        .unwrap()
    }
    #[test]
    fn products_reference_and_pool_identity() {
        // Includes empty, partial, exact panels and parallel-dispatch shapes.
        for (long, short) in [(0, 3), (1, 1), (255, 5), (256, 4), (257, 7), (1025, 33)] {
            for op in [RealOp::Identity, RealOp::Transpose] {
                let a = fixture(long, short);
                let b = fixture(if op == RealOp::Identity { short } else { long }, short);
                let m = if op == RealOp::Identity { long } else { short };
                let k = b.rows;
                let mut expected = RealMatrix::zeros(m, short).unwrap();
                for j in 0..short {
                    for i in 0..m {
                        for p in 0..k {
                            expected[(i, j)] += if op == RealOp::Identity {
                                a[(i, p)]
                            } else {
                                a[(p, i)]
                            } * b[(p, j)];
                        }
                    }
                }
                let ca = ComplexMatrix::new(
                    a.rows,
                    a.cols,
                    a.data.iter().map(|&x| Complex64::new(x, 0.)).collect(),
                )
                .unwrap();
                let cb = ComplexMatrix::new(
                    b.rows,
                    b.cols,
                    b.data.iter().map(|&x| Complex64::new(x, 0.)).collect(),
                )
                .unwrap();
                let mut cw = GemmWorkspace::new(long.max(short), short, 1).unwrap();
                let complex = cw
                    .multiply(
                        &ca,
                        if op == RealOp::Identity {
                            MatrixOp::Identity
                        } else {
                            MatrixOp::Adjoint
                        },
                        &cb,
                        MatrixOp::Identity,
                    )
                    .unwrap();
                let mut first = None;
                for threads in [1, 2, 4, 8] {
                    let pool = WorkerPool::new(threads).unwrap();
                    let mut ws = RealGemmWorkspace::with_pool(long, short, pool.clone()).unwrap();
                    let mut out = RealMatrix::zeros(m, short).unwrap();
                    ws.multiply_into(a.as_ref(), op, b.as_ref(), out.as_mut())
                        .unwrap();
                    // Unit-scale products: 1e-14 allows backend FMA/reduction rounding.
                    let defect = out
                        .data
                        .iter()
                        .zip(&expected.data)
                        .map(|(a, b)| (a - b).abs())
                        .fold(0., f64::max);
                    assert!(defect < 1e-14, "{defect:e}");
                    assert!(
                        out.data
                            .iter()
                            .zip(complex.as_slice())
                            .all(|(a, b)| (a - b.re).abs() < 1e-14 && b.im == 0.)
                    );
                    let bits: Vec<_> = out.data.iter().map(|x| x.to_bits()).collect();
                    if let Some(ref first) = first {
                        assert_eq!(&bits, first);
                    } else {
                        first = Some(bits);
                    }
                    // Resident/nested dispatch must retain the same arithmetic too.
                    pool.install(|| ws.multiply_into(a.as_ref(), op, b.as_ref(), out.as_mut()))
                        .unwrap();
                    assert_eq!(
                        out.data.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                        *first.as_ref().unwrap()
                    );
                    if threads == 1 {
                        println!("GEMM {op:?} {long}x{short}: naive defect={defect:e}");
                    }
                }
            }
        }
    }
    #[test]
    fn eigen_residual_sign_and_errors() {
        for n in [0, 1, 2, 9, 24] {
            let mut a = fixture(n, n);
            for j in 0..n {
                for i in 0..j {
                    a[(i, j)] = a[(j, i)];
                }
                a[(j, j)] += (j / 3) as f64;
            }
            let (values, v) = symmetric_eigh(&a).unwrap();
            assert!(values.windows(2).all(|w| w[0] <= w[1]));
            let mut residual: f64 = 0.;
            let mut orth: f64 = 0.;
            for j in 0..n {
                let mut pivot = 0;
                for i in 1..n {
                    if v[(i, j)].abs() > v[(pivot, j)].abs() {
                        pivot = i;
                    }
                }
                assert!(v[(pivot, j)] >= 0.);
                for i in 0..n {
                    residual = residual.max(
                        ((0..n).map(|k| a[(i, k)] * v[(k, j)]).sum::<f64>()
                            - values[j] * v[(i, j)])
                            .abs(),
                    );
                    orth = orth.max(
                        ((0..n).map(|k| v[(k, i)] * v[(k, j)]).sum::<f64>()
                            - if i == j { 1. } else { 0. })
                        .abs(),
                    );
                }
            }
            println!("EVD n={n}: residual={residual:e}, orthogonality={orth:e}");
            // Small matrices with O(1) spectrum: tight backward error allowance.
            assert!(residual < 1e-13 && orth < 1e-13);
            assert_eq!(symmetric_eigh(&a).unwrap(), (values, v));
        }
        assert!(RealMatrix::zeros(usize::MAX, 2).is_err());
        assert!(RealMatrix::new(2, 2, vec![0.]).is_err());
        assert!(symmetric_eigh(&fixture(2, 3)).is_err());
        let mut a = fixture(2, 2);
        a[(0, 1)] = 9.;
        assert!(symmetric_eigh(&a).is_err());
        a[(0, 0)] = f64::NAN;
        assert!(symmetric_eigh(&a).is_err());
        let pool = WorkerPool::new(1).unwrap();
        let mut ws = RealGemmWorkspace::with_pool(2, 2, pool).unwrap();
        let b = fixture(2, 2);
        let mut out = fixture(2, 2);
        assert!(
            ws.multiply_into(a.as_ref(), RealOp::Identity, b.as_ref(), out.as_mut())
                .is_err()
        );
        assert!(
            ws.multiply_into(
                b.as_ref().subrows(0, 1),
                RealOp::Transpose,
                b.as_ref(),
                out.as_mut()
            )
            .is_err()
        );
    }
}
