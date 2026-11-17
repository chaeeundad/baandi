//! First-correction-block residual diagnostic (ticket
//! BREAKTHROUGH-0-INSTRUMENT; review 2026-11-17 sections 6 and 14, experiment
//! E3; read-only, the solver is not changed).
//!
//! `BAANDI_FIRST_BLOCK=on` (default off). Before every warm real Gamma solve
//! (a start X of B vectors exists) the engine evaluates, at the CURRENT
//! Hamiltonian and on the solver's own rows:
//!
//! * a small Rayleigh-Ritz on the B kept vectors: A = X^T H X, theta and C
//!   (the overlap defect max |X^T X - I| is recorded; X is orthonormal to the
//!   restart rule in every kept start), the occupied / empty boundary from the
//!   occupation count, and the boundary gap theta[n_o] - theta[n_o - 1];
//! * the true residuals of the Ritz pairs, r_j = H x_j - theta_j x_j, for all B
//!   roots, and how many exceed the solver's per-root threshold (occupied
//!   roots: the solve tolerance; the last mandatory and the buffer roots: the
//!   buffer tolerance), against the number the Ritz-change rule corrects in its
//!   first iteration (all B roots, capped by the capacity rule);
//! * the singular values of the occupied residual block R_o = H X_o - X_o A_o
//!   (from the eigenvalues of its Gram matrix, one reduction), its effective
//!   rank at 1e-1, 1e-2, 1e-3 and 1e-4 of the largest singular value, and the
//!   number of singular values above the solve tolerance.
//!
//! When the kept H X of the start is armed (M2 kept start, REUSE-PSI) it is
//! used as is; otherwise H is applied to X once (B columns, visible in the
//! operator counters, not in the solver's `operator_applications`). The GEMMs
//! of the diagnostic appear in `gemm_dispatch` as unlabelled rows.
//!
//! Records are kept on rank 0 (the quantities are replicated) and printed as
//! one table per run (`table`), one line per solve with the context (stage,
//! kind, SCF iteration) of `crate::profile::context`.
use crate::profile;
use dft_mpi::Collective;
use dft_numerics::pool::WorkerPool;
use dft_numerics::real::{RealGemmWorkspace, RealMatrix, RealOp, symmetric_eigh};
use std::sync::{Arc, Mutex, OnceLock};

// 0: the switch; 1: forced on; 2: forced off (tests).
static FORCED: std::sync::atomic::AtomicU8 = std::sync::atomic::AtomicU8::new(0);
/// Tests: force the diagnostic on or off for the process (`None`: the switch).
pub fn force(mode: Option<bool>) {
    FORCED.store(
        match mode {
            None => 0,
            Some(true) => 1,
            Some(false) => 2,
        },
        std::sync::atomic::Ordering::Relaxed,
    );
}
/// `BAANDI_FIRST_BLOCK=on`, read once.
pub fn enabled() -> bool {
    match FORCED.load(std::sync::atomic::Ordering::Relaxed) {
        1 => return true,
        2 => return false,
        _ => {}
    }
    static ON: OnceLock<bool> = OnceLock::new();
    *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::FIRST_BLOCK))
}

/// One solve's diagnostic.
#[derive(Clone, Debug, PartialEq)]
pub struct Record {
    pub stage: &'static str,
    pub kind: &'static str,
    pub iteration: usize,
    pub bands: usize,
    pub occupied: usize,
    pub tolerance: f64,
    pub buffer_tolerance: f64,
    /// Whether the kept H X was used (else H was applied to X here).
    pub kept_hx: bool,
    pub overlap_defect: f64,
    pub boundary_gap: f64,
    /// Ritz values of the small problem (ascending).
    pub ritz_values: Vec<f64>,
    /// True residual norms of the B Ritz pairs.
    pub residual_norms: Vec<f64>,
    pub true_unconverged_all: usize,
    pub true_unconverged_occupied: usize,
    /// Roots the Ritz-change rule corrects in the solver's first iteration.
    pub ritz_rule_unconverged: usize,
    /// Singular values of R_o (descending).
    pub singular_values: Vec<f64>,
    pub rank_1e1: usize,
    pub rank_1e2: usize,
    pub rank_1e3: usize,
    pub rank_1e4: usize,
    /// Singular values above the solve tolerance.
    pub rank_above_tolerance: usize,
}

static RECORDS: Mutex<Vec<Record>> = Mutex::new(Vec::new());

/// Records collected since the last call.
pub fn take_records() -> Vec<Record> {
    std::mem::take(&mut *RECORDS.lock().expect("first-block records"))
}

/// Inputs of one evaluation.
pub(crate) struct Input<'a> {
    pub x: &'a RealMatrix,
    pub hx: Option<&'a RealMatrix>,
    pub bands: usize,
    pub occupied: usize,
    pub tolerance: f64,
    pub strict_bands: usize,
    pub buffer_tolerance: f64,
    pub capacity: usize,
    pub pool: Arc<WorkerPool>,
}

fn reduce<C: Collective>(comm: Option<&C>, local: &[f64], out: &mut [f64]) -> Result<(), String> {
    match comm {
        Some(c) => c.tree_sum_f64(local, out).map_err(|e| e.to_string()),
        None => {
            out.copy_from_slice(local);
            Ok(())
        }
    }
}

/// Evaluate and record (rank 0 keeps the record). `apply` forms H X on the
/// rows of `x` when no kept H X is given. Collective when `comm` is given.
pub(crate) fn evaluate<C: Collective>(
    comm: Option<&C>,
    input: Input<'_>,
    apply: &mut dyn FnMut(&RealMatrix, &mut RealMatrix) -> Result<(), String>,
) -> Result<(), String> {
    let Input {
        x,
        hx,
        bands,
        occupied,
        tolerance,
        strict_bands,
        buffer_tolerance,
        capacity,
        pool,
    } = input;
    let (nl, b) = (x.rows(), x.cols());
    if b != bands || bands == 0 || occupied == 0 || occupied > bands {
        return Ok(());
    }
    let e = |err: dft_numerics::linalg::LinalgError| err.to_string();
    let mut gemm = RealGemmWorkspace::with_pool(nl.max(1), b, pool).map_err(e)?;
    let owned;
    let kept_hx = hx.is_some();
    let hx = match hx {
        Some(h) if h.rows() == nl && h.cols() == b => h,
        _ => {
            let mut y = RealMatrix::zeros(nl, b).map_err(e)?;
            apply(x, &mut y)?;
            owned = y;
            &owned
        }
    };
    // S = X^T X and A = X^T H X, reduced together.
    let mut s_local = RealMatrix::zeros(b, b).map_err(e)?;
    let mut a_local = RealMatrix::zeros(b, b).map_err(e)?;
    if nl > 0 {
        gemm.multiply_into(x.as_ref(), RealOp::Transpose, x.as_ref(), s_local.as_mut())
            .map_err(e)?;
        gemm.multiply_into(x.as_ref(), RealOp::Transpose, hx.as_ref(), a_local.as_mut())
            .map_err(e)?;
    }
    let mut local = Vec::with_capacity(2 * b * b);
    local.extend_from_slice(s_local.as_slice());
    local.extend_from_slice(a_local.as_slice());
    let mut sum = vec![0.0; 2 * b * b];
    reduce(comm, &local, &mut sum)?;
    let mut overlap_defect = 0.0_f64;
    let mut a = RealMatrix::zeros(b, b).map_err(e)?;
    for j in 0..b {
        for i in 0..b {
            let s = sum[i + j * b];
            overlap_defect = overlap_defect.max((s - if i == j { 1.0 } else { 0.0 }).abs());
            // Symmetrize A (rounding of the two triangles).
            a[(i, j)] = 0.5 * (sum[b * b + i + j * b] + sum[b * b + j + i * b]);
        }
    }
    let (theta, c) = symmetric_eigh(&a).map_err(e)?;
    // Ritz vectors and their H images on the rows: X C and (H X) C.
    let mut xc = RealMatrix::zeros(nl, b).map_err(e)?;
    let mut hxc = RealMatrix::zeros(nl, b).map_err(e)?;
    if nl > 0 {
        gemm.multiply_into(x.as_ref(), RealOp::Identity, c.as_ref(), xc.as_mut())
            .map_err(e)?;
        gemm.multiply_into(hx.as_ref(), RealOp::Identity, c.as_ref(), hxc.as_mut())
            .map_err(e)?;
    }
    // Residual columns r_j = (H X C)_j - theta_j (X C)_j, written over hxc.
    let mut norms_local = vec![0.0; b];
    for j in 0..b {
        let col = &mut hxc.as_mut_slice()[j * nl..(j + 1) * nl];
        let xcol = &xc.as_slice()[j * nl..(j + 1) * nl];
        let mut n2 = 0.0;
        for (r, &xv) in col.iter_mut().zip(xcol) {
            *r -= theta[j] * xv;
            n2 += *r * *r;
        }
        norms_local[j] = n2;
    }
    // Gram of the occupied residual block, reduced with the norms.
    let no = occupied;
    let mut g_local = RealMatrix::zeros(no, no).map_err(e)?;
    if nl > 0 {
        let r_o = RealMatrix::new(nl, no, hxc.as_slice()[..nl * no].to_vec()).map_err(e)?;
        gemm.multiply_into(
            r_o.as_ref(),
            RealOp::Transpose,
            r_o.as_ref(),
            g_local.as_mut(),
        )
        .map_err(e)?;
    }
    let mut local = Vec::with_capacity(b + no * no);
    local.extend_from_slice(&norms_local);
    local.extend_from_slice(g_local.as_slice());
    let mut sum = vec![0.0; b + no * no];
    reduce(comm, &local, &mut sum)?;
    let residual_norms: Vec<f64> = sum[..b].iter().map(|v| v.max(0.0).sqrt()).collect();
    let mut g = RealMatrix::zeros(no, no).map_err(e)?;
    for j in 0..no {
        for i in 0..no {
            g[(i, j)] = 0.5 * (sum[b + i + j * no] + sum[b + j + i * no]);
        }
    }
    let (lambda, _) = symmetric_eigh(&g).map_err(e)?;
    let mut singular_values: Vec<f64> = lambda.iter().map(|v| v.max(0.0).sqrt()).collect();
    singular_values.sort_by(|a, b| b.total_cmp(a));
    let sigma_max = singular_values.first().copied().unwrap_or(0.0);
    let rank_at = |t: f64| {
        singular_values
            .iter()
            .filter(|&&s| s > t * sigma_max)
            .count()
    };
    let threshold = |j: usize| {
        if j + 1 < strict_bands {
            tolerance
        } else {
            buffer_tolerance
        }
    };
    let true_unconverged_all = (0..b).filter(|&j| residual_norms[j] > threshold(j)).count();
    let true_unconverged_occupied = (0..no)
        .filter(|&j| residual_norms[j] > threshold(j))
        .count();
    let record = Record {
        stage: profile::context().stage,
        kind: profile::context().kind,
        iteration: profile::context().iteration,
        bands: b,
        occupied: no,
        tolerance,
        buffer_tolerance,
        kept_hx,
        overlap_defect,
        boundary_gap: if no < b {
            theta[no] - theta[no - 1]
        } else {
            f64::NAN
        },
        ritz_values: theta,
        residual_norms,
        true_unconverged_all,
        true_unconverged_occupied,
        ritz_rule_unconverged: b.min(capacity.saturating_sub(b)),
        rank_1e1: rank_at(1e-1),
        rank_1e2: rank_at(1e-2),
        rank_1e3: rank_at(1e-3),
        rank_1e4: rank_at(1e-4),
        rank_above_tolerance: singular_values.iter().filter(|&&s| s > tolerance).count(),
        singular_values,
    };
    if comm.is_none_or(|c| c.rank() == 0) {
        RECORDS.lock().expect("first-block records").push(record);
    }
    Ok(())
}

/// Compact table of the records.
pub fn table(records: &[Record]) -> String {
    use std::fmt::Write;
    let mut out = String::new();
    let _ = writeln!(
        out,
        "first-block residual diagnostic ({} solves; B bands, n_o occupied; sigma = singular values of R_o = H X_o - X_o A_o; rank@t = count sigma > t sigma_max; unconverged = roots whose true residual exceeds their threshold, rule = roots the Ritz-change rule corrects first)",
        records.len()
    );
    let _ = writeln!(
        out,
        "  {:<4} {:<12} {:<9} {:>3} {:>4} {:>4} {:>9} {:>4} {:>9} {:>9} {:>9} {:>9} {:>5} {:>5} {:>5} {:>5} {:>6} {:>7} {:>7} {:>5} {:>9} {:>9}",
        "seq",
        "stage",
        "kind",
        "it",
        "B",
        "n_o",
        "tol",
        "kept",
        "gap",
        "sig_max",
        "sig_min",
        "r_max_occ",
        "r1e-1",
        "r1e-2",
        "r1e-3",
        "r1e-4",
        "r>tol",
        "unc_all",
        "unc_occ",
        "rule",
        "defect",
        "r_max_emp"
    );
    for (i, r) in records.iter().enumerate() {
        let sig_min = r.singular_values.last().copied().unwrap_or(0.0);
        let r_occ = r.residual_norms[..r.occupied]
            .iter()
            .copied()
            .fold(0.0, f64::max);
        let r_emp = r.residual_norms[r.occupied..]
            .iter()
            .copied()
            .fold(0.0, f64::max);
        let _ = writeln!(
            out,
            "  {:<4} {:<12} {:<9} {:>3} {:>4} {:>4} {:>9.2e} {:>4} {:>9.2e} {:>9.2e} {:>9.2e} {:>9.2e} {:>5} {:>5} {:>5} {:>5} {:>6} {:>7} {:>7} {:>5} {:>9.1e} {:>9.2e}",
            i,
            r.stage,
            r.kind,
            r.iteration,
            r.bands,
            r.occupied,
            r.tolerance,
            if r.kept_hx { "yes" } else { "no" },
            r.boundary_gap,
            r.singular_values.first().copied().unwrap_or(0.0),
            sig_min,
            r_occ,
            r.rank_1e1,
            r.rank_1e2,
            r.rank_1e3,
            r.rank_1e4,
            r.rank_above_tolerance,
            r.true_unconverged_all,
            r.true_unconverged_occupied,
            r.ritz_rule_unconverged,
            r.overlap_defect,
            r_emp
        );
    }
    out
}
