//! F11 (review 2026-11-16): effective backend, shape and fallback reason of
//! every real GEMM wrapper call of `real::RealGemmWorkspace`. Count-only: per
//! (site, route) a call count, 2 m n k flops and the largest m, n, k seen, in
//! relaxed atomics (no allocation, no timing, no I/O; the dispatch decision and
//! the arithmetic are unchanged). Exported in results.json
//! `convergence.work.gemm_dispatch`. In the source-only GEMM accounting build
//! (`P14B_GEMM_ACCOUNTING=1`) every call also becomes an accounting record
//! named `real.<site>.<route>` with its exact shape and time.
//!
//! Shapes: the output is m x n and the reduction length k (A B: A is m x k;
//! A^T B: A is k x m), whatever the wrapper's argument order.
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

/// Call site label (`real::GemmSite` plus `Unlabelled` for the direct
/// `multiply_into` / `rotate_pair` entries and the B5 column groups).
pub const SITES: [&str; 7] = [
    "rotate",
    "gram",
    "projection_overlap",
    "projection_apply",
    "reduced",
    "other",
    "unlabelled",
];
pub const UNLABELLED: usize = 6;

/// Effective route of one call.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Route {
    /// MKL dgemm panels (fixed 256-row or 256-term panels).
    Mkl = 0,
    /// MKL dgemm panels that the wide entry made possible: the whole output
    /// m n exceeded short^2 and only the per-panel bound held.
    MklWide = 1,
    /// Pure-Rust kernel: build without the `mkl` feature.
    PureBuild = 2,
    /// Pure: MKL build with `BAANDI_MKL_DGEMM=off`.
    PureSwitch = 3,
    /// Pure fallback: the long dimension of A exceeds the workspace rows.
    PureLong = 4,
    /// Pure fallback: A's columns or the output columns exceed `short`.
    PureShort = 5,
    /// Pure fallback: whole-output bound m n > short^2 on a non-wide call (or a
    /// transposed product, which the wide relaxation does not cover).
    PureOutput = 6,
    /// Pure fallback: one 256-row panel of the output exceeds a tile.
    PurePanel = 7,
    /// Pure fallback: empty shape, leading dimension or output view mismatch.
    PureLayout = 8,
    /// Direct call of the pure kernel (`multiply_into`, `rotate_pair`).
    PureDirect = 9,
    /// Streaming dot / update kernels of the NARROW switches (not a GEMM).
    Narrow = 10,
    /// One sequential MKL call over the complete local reduction.
    MklSingle = 11,
}
pub const ROUTES: [&str; 12] = [
    "mkl",
    "mkl_wide",
    "pure_build",
    "pure_switch",
    "pure_long",
    "pure_short",
    "pure_output_bound",
    "pure_panel_bound",
    "pure_layout",
    "pure_direct",
    "narrow",
    "mkl_single",
];
const S: usize = SITES.len();
const R: usize = ROUTES.len();
/// Per (site, route): calls, flops, max m, max n, max k.
const FIELDS: usize = 5;
type Table = [[[AtomicU64; FIELDS]; R]; S];
static COUNTS: Table = [const { [const { [const { AtomicU64::new(0) }; FIELDS] }; R] }; S];
/// Pure builds only: the route an MKL build with dgemm enabled would take for
/// the same call on the same workspace (same rule), so a pure run shows the
/// MKL fallbacks of its shapes.
static AS_IF_MKL: Table = [const { [const { [const { AtomicU64::new(0) }; FIELDS] }; R] }; S];
/// Calls through `multiply_slices_wide_into` per site (requested, whether or
/// not the relaxation was needed).
static WIDE_REQUESTS: [AtomicU64; S] = [const { AtomicU64::new(0) }; S];

/// Record one call (output m x n, reduction k).
#[inline]
pub(crate) fn record(site: usize, route: Route, m: usize, n: usize, k: usize) {
    add(&COUNTS, site, route, m, n, k);
}
/// Pure builds: record the route an MKL build would take.
#[inline]
pub(crate) fn record_as_if_mkl(site: usize, route: Route, m: usize, n: usize, k: usize) {
    add(&AS_IF_MKL, site, route, m, n, k);
}
#[inline]
fn add(table: &Table, site: usize, route: Route, m: usize, n: usize, k: usize) {
    let c = &table[site.min(S - 1)][route as usize];
    c[0].fetch_add(1, Relaxed);
    c[1].fetch_add(2 * (m as u64) * (n as u64) * (k as u64), Relaxed);
    c[2].fetch_max(m as u64, Relaxed);
    c[3].fetch_max(n as u64, Relaxed);
    c[4].fetch_max(k as u64, Relaxed);
}
#[inline]
pub(crate) fn record_wide_request(site: usize) {
    WIDE_REQUESTS[site.min(S - 1)].fetch_add(1, Relaxed);
}

/// One non-empty (site, route) cell.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DispatchRow {
    pub site: &'static str,
    pub route: &'static str,
    pub calls: u64,
    pub flops: u64,
    pub max_m: u64,
    pub max_n: u64,
    pub max_k: u64,
}
impl DispatchRow {
    /// Whether this route ran the vendor library.
    pub fn is_mkl(&self) -> bool {
        matches!(self.route, "mkl" | "mkl_wide" | "mkl_single")
    }
}

/// Every (site, route) with at least one call, in table order.
pub fn rows() -> Vec<DispatchRow> {
    table_rows(&COUNTS)
}
/// Pure builds: the same calls under the route an MKL build with dgemm
/// enabled would take (empty on MKL builds, whose `rows` are the real routes).
pub fn rows_as_if_mkl() -> Vec<DispatchRow> {
    table_rows(&AS_IF_MKL)
}
fn table_rows(table: &Table) -> Vec<DispatchRow> {
    let mut out = Vec::new();
    for (s, site) in SITES.iter().enumerate() {
        for (r, route) in ROUTES.iter().enumerate() {
            let c = &table[s][r];
            let calls = c[0].load(Relaxed);
            if calls > 0 {
                out.push(DispatchRow {
                    site,
                    route,
                    calls,
                    flops: c[1].load(Relaxed),
                    max_m: c[2].load(Relaxed),
                    max_n: c[3].load(Relaxed),
                    max_k: c[4].load(Relaxed),
                });
            }
        }
    }
    out
}

/// Wide-entry calls per site (site, calls), non-zero only.
pub fn wide_requests() -> Vec<(&'static str, u64)> {
    SITES
        .iter()
        .zip(&WIDE_REQUESTS)
        .map(|(s, c)| (*s, c.load(Relaxed)))
        .filter(|(_, c)| *c > 0)
        .collect()
}

/// Accounting-build name of a (site, route) cell (`real.<site>.<route>`).
/// Built once, only in the source-only accounting build.
pub(crate) fn accounting_name(site: usize, route: Route) -> &'static str {
    use std::sync::OnceLock;
    static NAMES: OnceLock<Vec<&'static str>> = OnceLock::new();
    let names = NAMES.get_or_init(|| {
        SITES
            .iter()
            .flat_map(|s| {
                ROUTES
                    .iter()
                    .map(move |r| &*Box::leak(format!("real.{s}.{r}").into_boxed_str()))
            })
            .collect()
    });
    names[site.min(S - 1) * R + route as usize]
}

/// Accounting record of one call (`P14B_GEMM_ACCOUNTING` builds only; a no-op
/// span otherwise).
#[inline]
pub(crate) fn span(
    site: usize,
    route: Route,
    m: usize,
    n: usize,
    k: usize,
) -> crate::linalg::accounting::Span {
    let name = if crate::linalg::accounting::ENABLED {
        accounting_name(site, route)
    } else {
        ""
    };
    crate::linalg::accounting::Span::new(0, 0, name, [m, n, k], 1, 0)
}
