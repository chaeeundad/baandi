//! Switch manifest (work order A1): every `BAANDI_*` environment switch of the
//! workspace in one static table, read once per process.
//!
//! Each entry states its name, default, allowed values and meaning. The first
//! accessor call reads all switches from the environment, parses them with
//! exactly the rules the call sites used before this module existed (including
//! the fallback on an invalid value) and keeps the effective value and its
//! source (`default`, `env`, or `invalid -> <fallback>`). [`manifest`] returns
//! the whole table for the result file (`scf.switches` in the CLI output).
//!
//! No other source file reads a `BAANDI_*` variable; a workspace test scans
//! `crates/*/src` for that. Values are fixed for the process lifetime: setting a
//! variable after the first accessor call has no effect.
use std::ffi::{OsStr, OsString};
use std::sync::OnceLock;

/// How a switch value is parsed.
#[derive(Debug, Clone, Copy)]
pub enum Kind {
    /// One of `allowed`; any other value (or no value) gives `default`.
    Choice {
        allowed: &'static [&'static str],
        default: &'static str,
    },
    /// A number accepted by `valid`; anything else gives `default`.
    Number {
        default: Option<f64>,
        valid: fn(f64) -> bool,
    },
    /// `BAANDI_ORTHO`: unset gives the DGKS ratio 0.1, `dgks` the classical
    /// 1/sqrt(2), `dgks:<r>` with r in [1e-3, 1] that ratio; any other value
    /// (`off` included) gives `None`, two orthogonalization passes always.
    Ortho,
    /// Overrides an input-file value; unset keeps the file value. A value that
    /// fails `valid` is an input error at the call site (no fallback).
    Override { valid: fn(&str) -> bool },
    /// Acts when the variable is set to any value.
    Presence,
    /// An absolute path; a relative or empty one is an error at the call site.
    Path,
}

/// The effective value of a switch.
#[derive(Debug, Clone, PartialEq)]
pub enum Value {
    Choice(&'static str),
    Number(Option<f64>),
    /// Raw override or path text; `None` keeps the input file / built-in default.
    Raw(Option<String>),
    Present(bool),
}

/// Where the effective value came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Source {
    Default,
    Env,
    /// The variable was set but not accepted; the call site falls back to the
    /// default value.
    InvalidToDefault,
    /// `BAANDI_ORTHO` only: an unrecognized value selects two passes always.
    InvalidToOff,
    /// Override and path switches: the call site reports an input error.
    InvalidToError,
}
impl Source {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Default => "default",
            Self::Env => "env",
            Self::InvalidToDefault => "invalid -> default",
            Self::InvalidToOff => "invalid -> off",
            Self::InvalidToError => "invalid -> error",
        }
    }
}

/// One table entry.
#[derive(Debug, Clone, Copy)]
pub struct Switch {
    pub id: Id,
    pub name: &'static str,
    pub kind: Kind,
    /// Default as text (what the result file shows when the variable is unset).
    pub default: &'static str,
    /// Allowed values as text.
    pub allowed: &'static str,
    pub meaning: &'static str,
}

/// One manifest row: the table entry with its effective value and source.
#[derive(Debug, Clone)]
pub struct Entry {
    pub switch: &'static Switch,
    pub value: Value,
    pub source: Source,
}

const ON_OFF: &[&str] = &["on", "off"];

fn eigen_factor_valid(v: f64) -> bool {
    v.is_finite() && (1e-3..=1.0).contains(&v)
}
fn buffer_tol_valid(v: f64) -> bool {
    v.is_finite() && (1e-5..=1e-1).contains(&v)
}
fn growth_reset_valid(v: f64) -> bool {
    v.is_finite() && v > 1.0
}
fn b1_margin_valid(v: f64) -> bool {
    v.is_finite() && v > 0.0 && v < 1.0
}
fn b1_rows_valid(v: f64) -> bool {
    v.fract() == 0.0 && (8.0..=65536.0).contains(&v)
}
fn coarse_stop_valid(v: f64) -> bool {
    v.is_finite() && v > 0.0 && v < 1.0
}
fn advance_blocks_valid(v: f64) -> bool {
    v.fract() == 0.0 && (1.0..=8.0).contains(&v)
}
fn cutoff_ladder_valid(v: &str) -> bool {
    matches!(v, "off" | "on")
        || v.parse::<f64>()
            .is_ok_and(|f| f.is_finite() && f > 0.0 && f < 1.0)
}
fn relax_extrapolation_valid(v: &str) -> bool {
    matches!(v, "none" | "first_order" | "second_order")
}

macro_rules! table {
    ($($id:ident: $kind:expr, $default:expr, $allowed:expr, $meaning:expr;)*) => {
        /// Switch identifiers; the variable name is `BAANDI_<id>`.
        #[allow(non_camel_case_types)]
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        pub enum Id { $($id),* }
        /// The whole table, in `Id` order.
        pub static TABLE: &[Switch] = &[$(Switch {
            id: Id::$id,
            name: concat!("BAANDI_", stringify!($id)),
            kind: $kind,
            default: $default,
            allowed: $allowed,
            meaning: $meaning,
        }),*];
    };
}

table! {
    // Engine (dft-engine).
    CERT_SCREEN: Kind::Choice { allowed: &["on", "off", "trace"], default: "on" }, "on", "on, off, trace",
        "certificate screen (trace prints the estimates)";
    POTENTIAL_REUSE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "certificate solve reuses the screen's Hartree/XC fields";
    MIX_GROWTH_RESET: Kind::Number { default: None, valid: growth_reset_valid }, "none", "number > 1",
        "reset the mixing history when the residual grows by more than this factor";
    EIGEN_FACTOR: Kind::Number { default: Some(0.1), valid: eigen_factor_valid }, "0.1", "number in [1e-3, 1]",
        "eigen tolerance = factor x density residual";
    EMPTY_LOOSE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "first empty band is a buffer band in steps above the tolerance floor";
    BUFFER_TOL: Kind::Number { default: Some(1e-3), valid: buffer_tol_valid }, "1e-3", "number in [1e-5, 1e-1]",
        "residual floor (Ha) of empty buffer bands above the tolerance floor";
    NL_ENERGY: Kind::Choice { allowed: &["fast", "exact"], default: "fast" }, "fast", "fast, exact",
        "nonlocal energy from blocked projections (exact: band-by-band operator)";
    ENERGY_IDENTITY: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "nonlocal energy of SCF iterations from the band-sum identity";
    CUTOFF_LADDER_STANDARD: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "cutoff ladder also for the standard profile (S1)";
    LADDER_STANDARD_SEED: Kind::Choice { allowed: &["coarse", "cap"], default: "coarse" }, "coarse",
        "coarse, cap", "first fine eigen tolerance of the standard ladder: from the coarse ratchet or the cold-start cap";
    LADDER_SECANT: Kind::Choice {
        allowed: &["off", "on", "1", "2", "3", "4", "6", "8", "16", "32", "64"],
        default: "off",
    }, "off", "off, on (4), 1..64 difference pairs", "secant pairs carried across the ladder hand-off (S2)";
    // D1 ladder hand-off scan (work order 2026-11-14); defaults are the adopted ladder.
    LADDER_COARSE_STOP: Kind::Number { default: Some(1e-4), valid: coarse_stop_valid }, "1e-4",
        "number in (0, 1)", "cutoff ladder: coarse stage stops when its density residual (Bohr^-3/2) is at or below this level";
    LADDER_COARSE_RULE: Kind::Choice { allowed: &["fixed", "atomic"], default: "fixed" }, "fixed",
        "fixed, atomic", "cutoff ladder coarse stop level: fixed (LADDER_COARSE_STOP), or atomic (the larger of that level and the distance between the full-cutoff and prolonged coarse atomic start densities, same norm)";
    LADDER_FIRST_TOL: Kind::Choice { allowed: &["auto", "cap", "coarse", "next"], default: "auto" }, "auto",
        "auto, cap, coarse, next", "first fine eigen tolerance after the ladder hand-off: auto (standard per LADDER_STANDARD_SEED, strict cold 1e-2), cap (both cold), coarse (the last coarse step's requested tolerance), next (the coarse schedule's next tolerance: standard ratchet, strict factor x last coarse residual)";
    STANDARD_STOP: Kind::Choice { allowed: &["plain", "confirm"], default: "plain" }, "plain",
        "plain, confirm", "standard-class stop rule: plain threshold, or confirm (also stop when the estimator is within 2x the threshold, the previous step within 20x, and it contracted by 0.3 or better)";
    STANDARD_RESIDUAL_GUARD: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "standard-class stop also requires sum_k w_k sum_occupied f_i |r_i|^2 x 1 Ha^-1 (empirical second-order band-energy estimator, not a bound; review F03) below the threshold";
    DENSITY_PAIR: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "two bands per inverse transform in the density build";
    DENSITY_REDUCE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "block-parallel band sum of the density, zero-occupation bands skipped";
    NL_ENERGY_REAL: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "audit nonlocal energies from one real projector product on the real Gamma path";
    ATOMIC_REAL: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "atomic-start Ritz step through the real Gamma Hamiltonian";
    MPI_M2: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "MPI stick-row (M2) path; off selects the M1 row-block path";
    M2_BATCH: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "batched FFT transposes of the M2 path";
    PHASE_THREADS: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "P48F: thread the atomic start (span build, two-pass Gram-Schmidt right-looking) on the SCF pool; bit-identical";
    REUSE_PSI: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "C2: a lean real Gamma Davidson solve that starts from the previous solve's orbitals X takes H X from that solve's rotated H images plus (V_new - V_old) X, formed from the psi(r) the density build kept (one forward transform per pair; no inverse transform, nonlocal product or kinetic pass); falls back to the plain apply when X, the operator or the grid changed or the cache exceeds its memory gate; rounding class";
    SCF_ADVANCE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "BREAKTHROUGH-3 (docs/design/bounded-subspace-scf.md): intermediate standard-path maps of the serial real Gamma lean path run a bounded subspace update (`Advance`: at most SCF_ADVANCE_BLOCKS correction blocks from the kept start) instead of a solve to the ratchet tolerance; the map is completed by the exact solve at the same density when its residual estimate exceeds 0.1 x the Coulomb SCF residual, when the SCF residual is below the threshold (convergence candidate) or grew, or when the occupied/empty Ritz values cross; the stop rule, the verification and the delivered state use exact maps only; trajectory change, standard class unchanged";
    SCF_ADVANCE_BLOCKS: Kind::Number { default: Some(2.0), valid: advance_blocks_valid }, "2", "integer in [1, 8]",
        "starting number of correction blocks per advanced map (SCF_ADVANCE=on); the session adapts it by one block from the inner/SCF residual ratio of the previous advanced map";
    SCF_ADVANCE_GROWTH: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "SCF_ADVANCE research: complete an advanced map whose Coulomb SCF residual grew against the previous map (off: growth is left to the mixer, as on the exact path)";
    SCF_ADVANCE_ADAPT: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "SCF_ADVANCE research: adapt the block budget from the previous advanced map's inner/SCF ratio (off: the budget stays at SCF_ADVANCE_BLOCKS)";
    SCF_ADVANCE_SCOPE: Kind::Choice { allowed: &["all", "fine"], default: "all" }, "all", "all, fine",
        "SCF_ADVANCE research: advanced maps in every standard session (all) or only outside the cutoff-ladder coarse stage (fine)";
    VERIFY_REUSE: Kind::Choice { allowed: &["off", "warm", "feed", "feed2", "replace", "half"], default: "half" }, "half", "off, warm, feed, feed2, replace, half",
        "reuse of a failed fixed-point verification solve (F at the delivered density rho_out): warm = the next SCF step starts from its orbitals (so the REUSE_PSI kept start is taken); feed = warm and the pair (rho_out, F(rho_out)) also enters the mixer after the step's own pair, so the next input density is mixed from it (the verification map counts as an SCF map); off = discarded; rounding class (trajectory changes, acceptance still needs a passed verification)";
    // MPI entry cost (docs/performance/mpi-entry-cost-analysis.md). Adopted 2026-11-15 on dev08 A/B: M2_DGEMM, M2_ORTHO=serial and M2_LEAN default on; M2_FFT stays rust (neutral).
    M2_DGEMM: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "distributed Davidson rotations and projections and the M2 nonlocal products through the per-panel dgemm at any width (MKL builds: rounding class; pure builds bit-identical)";
    M2_ORTHO: Kind::Choice { allowed: &["classic", "serial"], default: "serial" }, "serial",
        "classic, serial",
        "distributed Davidson orthonormalization: classic (three basis projections, two whole-block sweeps) or serial (the serial solver's order: one projection, DGKS second pass, panel Gram-Schmidt; rounding class)";
    M2_FFT: Kind::Choice { allowed: &["rust", "mkl"], default: "rust" }, "rust", "rust, mkl",
        "1D kernels of the M2 distributed FFT (mkl: MKL DFTI batch lanes on MKL builds, rounding class; pure builds stay rustfft)";
    M2_LEAN: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "distributed Davidson residuals and corrections by slice loops without per-element index checks; bit-identical";
    M2_CAPACITY: Kind::Choice { allowed: &["auto", "off"], default: "auto" }, "auto", "auto, off",
        "MPI runs with four or more ranks use 3 Davidson basis columns per band when scf.davidson_capacity_per_band is not set (auto), else the default 6 (off)";
    STAGE_TIMERS: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "F21: MPI stage and collective counters measure wall time (on, the earlier behaviour) or run count-only (off: calls and collective payload bytes counted, no clock read, stage seconds 0)";
    // MPI reduced eigenproblem (docs/performance/mpi-reduced-eigen-analysis.md); default off.
    M2_EIGEN: Kind::Choice { allowed: &["replicated", "root"], default: "replicated" }, "replicated",
        "replicated, root",
        "distributed Davidson reduced eigenproblem: replicated (every rank solves) or root (rank 0 solves and broadcasts the lowest pairs; bit-identical)";
    M2_B5: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "distributed Davidson bordered (secular) reduced eigen update when one column was appended (the serial B5; rounding class)";
    // MPI solver workspace (review F06; ledger row M2-WORKSPACE-F06); bit-identical.
    M2_WORKSPACE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "distributed Davidson keeps its basis, products, reduced matrix (fixed leading dimension, border update), Ritz, residual and orthonormalization buffers, GEMM and EVD workspaces in a per-rank workspace reused across iterations and solves (no per-iteration large allocation); bit-identical to off";
    // MPI selected-root rotation (review F05; ledger row M2-SELECTED-ROOTS-F05); default off until dev08.
    M2_SELECTED_ROOTS: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "distributed Davidson Ritz-change (standard) solves: unconverged roots decided from the Ritz-value change first, Ritz vectors, H images and residuals of those roots only, all bands rotated only at return and restarts (the serial lean rule); strict solves unchanged; implies the workspace path; rounding class";
    // MPI in-panel Gram-Schmidt reductions (docs/performance/profile-2026-11-16-mpi.md); default classic until dev08.
    M2_PANEL_GS: Kind::Choice { allowed: &["classic", "pythag", "fused"], default: "fused" }, "fused",
        "classic, pythag, fused",
        "distributed panel Gram-Schmidt reductions per column: classic (two dot reductions and one norm), pythag (norm rides with the second dots, Pythagorean update, 2 reductions), fused (input norm rides with the first dots, second pass only when the per-column DGKS rule of BAANDI_ORTHO fails, about 1 reduction); rounding class";
    M2_DENSITY: Kind::Choice { allowed: &["bands", "fused"], default: "fused" }, "fused",
        "bands, fused",
        "M2 distributed density on the standard path: bands (per-band plane contributions kept, then summed) or fused (plane density accumulated in ascending band order, no bands x plane buffer; the per-band data only serves the strict certificate screen); bit-identical";
    // MPI kept start (ticket MPI-KEPT-START, docs/design/mpi-kept-start.md); default off until dev08.
    M2_REUSE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "MPI-KEPT-START: a standard (Ritz-change) M2 stick-row solve that starts from the previous solve's orbitals X takes X and H X from that solve's local rows (the return rotation of the H images) plus (V_new - V_old) X, formed from the plane slabs the M2 density build kept (one forward transform per occupied pair, no inverse transform, nonlocal product or kinetic pass), when initial equals X to 1e-14 and the Gram defect of X is at most 1e-12; the local rows of X come from the solver (no extraction from the gathered block); falls back to the plain start otherwise; rounding class";
    // MPI scaling gap (ticket MPI-SCALING-GAP, docs/performance/mpi-scaling-gap.md); default off.
    M2_STICKS: Kind::Choice { allowed: &["greedy", "mirror"], default: "mirror" }, "mirror",
        "greedy, mirror",
        "stick partition of the M2 distributed FFT: greedy (sticks one by one to the least loaded rank) or mirror (a stick and its Gamma mirror stick go to the same rank as one unit, so no ORC representative is split and the two mirror alltoalls of every apply vanish; rounding class: rows are summed in a different order)";
    M2_START: Kind::Choice { allowed: &["ortho", "gram"], default: "gram" }, "gram",
        "ortho, gram",
        "start of a distributed Davidson solve (workspace path): ortho (panel Gram-Schmidt of the start block) or gram (one Gram matrix and one reduction; a start block with Gram defect at most 1e-12, the restart contraction rule, is taken as the basis without Gram-Schmidt, otherwise ortho; rounding class)";
    // MPI overlap-matrix Davidson (ticket MPI-OVERLAP-DAVIDSON, docs/design/mpi-overlap-davidson.md); default off.
    M2_LEAN_ALGEBRA: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "MPI-OVERLAP-DAVIDSON: standard (Ritz-change) distributed Davidson solves run the serial lean algebra: unit-norm corrections appended without Gram-Schmidt, the overlap S = V^T V kept beside H_r with its inverse Cholesky factor (pivot rule eps/1e-12, reprojection and Ritz-basis fallbacks), reduced problem T H_r T^T; the S border, the H border and the selected residual norms share one tree sum per iteration; strict solves unchanged; implies the workspace path; rounding class";
    // MPI lazy orbitals (ticket MPI-LAZY-ORBITALS, docs/design/mpi-lazy-orbitals.md); default off.
    M2_LAZY: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "MPI-LAZY-ORBITALS: a standard (Ritz-change) M2 solve inside the SCF session returns only this rank's stick rows of the eigenvectors (no N x B allgatherv, no complex unpack); the next solve starts from those rows, the density, nonlocal and kinetic energy are formed from them (kinetic: one reduction of B values), and the full complex block is gathered once, on first demand (result, checkpoint, ladder hand-off, spectral audit, explicit nonlocal audit), and cached; gathered blocks are bit-identical, the energy is rounding class (distributed kinetic sum)";
    MPI_TREE: Kind::Choice { allowed: &["binomial", "rd", "rhd"], default: "rd" }, "rd",
        "binomial, rd, rhd",
        "message pattern of the deterministic MPI tree sum: binomial (reduce to rank 0 plus broadcast, 2 log2 P rounds), rd (recursive doubling with MPI_Sendrecv, log2 P rounds) or rhd (REL-TREE-RHD: payloads of at least 2048 doubles by a recursive-halving reduce-scatter plus a recursive-doubling allgather, 2 (P-1)/P of the payload moved instead of log2 P times it; smaller payloads rd); all bit-identical for power-of-two rank counts, other counts stay binomial; read by the CLI at MPI start";
    // MPI coarse stage and certificate (docs/performance/mpi-coarse-certificate-analysis.md); default off.
    M2_SCREEN: Kind::Choice { allowed: &["replicated", "dist"], default: "dist" }, "dist",
        "replicated, dist",
        "certificate screen clusters on the M2 path: replicated (every rank transforms the cluster bands on the full grid) or dist (band values on each rank's FFT planes, local sums plus one tree_sum per cluster; rounding class of the screen estimate only)";
    M2_ATOMIC: Kind::Choice { allowed: &["replicated", "split"], default: "split" }, "split",
        "replicated, split",
        "atomic start on the M2 path: replicated, or split (span columns, two-pass Gram-Schmidt first passes, H Q and the Gram columns split over ranks with one broadcast per column and one allgatherv; Gram-Schmidt bit-identical, Gram block GEMM rounding class)";
    // Replicated grid stages outside the solve (ticket MPI-REMAINING, docs/performance/mpi-remaining.md); default off.
    STOP_METRIC: Kind::Choice { allowed: &["fft", "kept"], default: "kept" }, "kept", "fft, kept",
        "standard stop-rule Coulomb estimate E_H[rho_out - rho_in]: fft (transform the residual forward and back on the full grid) or kept (Parseval sum over the map's own Fourier coefficients of rho_in, kept from the Hartree potential, and rho_out, kept from the energy audit; no FFT; falls back to fft when a density was clamped or the kept coefficients do not belong to the two densities; rounding class of the estimate only)";
    // Solvers (dft-solvers).
    MIX_INCREMENTAL: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "Pulay mixer keeps its difference vectors (and their metric images) and the Gram matrix entries across calls and forms only the newest difference and Gram row per step (O(m N) instead of O(m^2 N) grid work, no per-step difference allocation; ticket MPI-REMAINING); bit-identical";
    B1: Kind::Choice { allowed: &["on", "off", "standard", "strict"], default: "on" }, "on",
        "on, off, standard, strict", "B1 selected-root Davidson iterations";
    B1_SKETCH: Kind::Choice { allowed: &["count", "dense"], default: "count" }, "count", "count, dense",
        "B1 residual sketch";
    B1_NARROW: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "B1 narrow iteration whenever the explicit root set is smaller than the band count (off: fewer than half the bands)";
    B1_MARGIN: Kind::Number { default: Some(0.5), valid: b1_margin_valid }, "0.5", "number in (0, 1)",
        "B1 strict: a root whose sketch estimate is below margin x its threshold skips its explicit residual for one iteration";
    B1_ROWS: Kind::Number { default: None, valid: b1_rows_valid }, "auto", "integer in [8, 65536]",
        "B1 strict sketch rows (auto: ceil(32 / (1 - margin)^2), 128 at margin 0.5, same 8-sigma miss bound at any margin)";
    B5: Kind::Choice { allowed: &["on", "off", "par"], default: "on" }, "on", "on, off, par",
        "bordered reduced eigenproblem update (par: threaded)";
    B6: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "B6 projector-product reuse";
    RITZ_FUSED: Kind::Choice { allowed: &["on", "off", "loop"], default: "off" }, "off", "on, off, loop",
        "fused Ritz rotation (loop: only the selected-root residuals inside the lean Davidson iteration, where H X is not kept; the return keeps the pair rotation that REUSE-PSI needs; PERCORE-KERNELS)";
    NARROW_PROJ: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "narrow-iteration projection kernel";
    NARROW_REDUCED: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "narrow-iteration reduced-matrix kernel";
    DAVIDSON_STAGES: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "Davidson stage timers";
    SECOND_VECTOR: Kind::Choice { allowed: &["off", "v1", "v2", "v3", "v4"], default: "off" }, "off",
        "off, v1, v2, v3, v4", "second correction vector variant";
    ORTHO: Kind::Ortho, "dgks:0.1", "off, dgks, dgks:<ratio in [1e-3, 1]>",
        "DGKS second-pass rule of the Davidson orthogonalization (value = ratio, null = two passes)";
    MGS: Kind::Choice {
        allowed: &["blocked", "classic"],
        default: if crate::MKL_BACKEND { "blocked" } else { "classic" },
    }, if crate::MKL_BACKEND { "blocked" } else { "classic" }, "blocked, classic",
        "Gram-Schmidt variant (default blocked on MKL builds, classic otherwise)";
    ORTHO_BLOCK: Kind::Choice {
        allowed: ON_OFF,
        default: if crate::MKL_BACKEND { "off" } else { "on" },
    }, if crate::MKL_BACKEND { "off" } else { "on" }, "on, off",
        "panel Gram-Schmidt of the correction block, L2 shape-rule panel width (P48D; overrides MGS; default on in pure builds where it replaces the serial classic sweep, off on MKL builds where MGS-B already blocks)";
    WIDE_DOT: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "wide dot products in the orthogonalization";
    PROJECT_WIDE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "basis-projection apply X -= Q S and MGS-B cross-panel apply through the per-panel dgemm at any block width (MKL builds; pure builds unchanged)";
    DAVIDSON: Kind::Choice { allowed: &["classic", "lean"], default: "lean" }, "lean", "classic, lean",
        "serial real Gamma Davidson of Ritz-value (standard) solves: classic (orthonormal basis, Ritz vectors of all bands in wide iterations) or lean (work order phase B: residuals of unconverged roots only, all-band Ritz vectors at return and contraction, unit-norm corrections without projection and a generalized reduced problem by Cholesky reduction; one orthonormalization when a pivot ratio of the basis overlap falls below eps/1e-12; rounding class)";
    LEAN_PROJ_CACHE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "C3a: the lean real Gamma Davidson keeps P^T V (projector products) beside its basis (copied after each H apply, rotated at contraction and return) and the audit's nonlocal energy reads P^T X = (P^T V) C instead of a fresh n x n_proj x bands product; falls back to the direct product whenever the cache is incomplete; rounding class";
    // Per-core kernels (ticket PERCORE-KERNELS, docs/performance/percore-kernels-2026-11-17.md); default off.
    LEAN_BORDER: Kind::Choice { allowed: &["separate", "fused"], default: "separate" }, "separate",
        "separate, fused",
        "PERCORE-KERNELS: borders of the lean real Gamma Davidson: separate (overlap border V^T E before the H apply, reduced border V^T (H E) in the next iteration, two transposed products) or fused (H applied to the unit corrections first, then one transposed product V^T [E | H E] of twice the width forms both borders, so the basis is streamed once; the pivot test follows it and a failing append discards the images, takes the separate fallback and ends fusing for the rest of that solve; only blocks with 4 k <= m (basis width m) fuse, wider blocks such as the first append of every solve keep the separate products; bit-identical on the aarch64 pure kernel, rounding class on MKL and x86)";
    // History-preserving restart (ticket BREAKTHROUGH-1-RESTART, docs/design/restart-keep.md); default off.
    RESTART_KEEP: Kind::Choice {
        allowed: &["off", "1", "2", "4", "8", "b/8", "b/4", "b/2"],
        default: "off",
    }, "off", "off, 1, 2, 4, 8, b/8, b/4, b/2",
        "BREAKTHROUGH-1-RESTART R1-b: at a contraction restart of the lean Davidson (serial lean path and the distributed lean algebra) keep, besides the Ritz vectors X and H X, up to r previous Ritz vectors of the unconverged roots with their X component removed (P = Y - X A, H P = (H V)(D - C A): no new H application); r is a count or a fraction of the band count (b/4 = bands / 4, at least 1), limited by the room left for the next correction block; directions whose squared sine to span(X and the kept directions) falls below the pivot rule eps/1e-12 are dropped, and the restart falls back to X alone when the factor of the new basis fails the pivot rule; off = the plain restart (bit-identical)";
    // Physics (dft-physics).
    MKL_NONLOCAL: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "nonlocal products by dgemm on the real Gamma path";
    APPLY_LEAN: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "lean real Gamma H apply (no redundant scans)";
    HAPPLY_OVERLAP: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "overlap local and nonlocal parts of the H apply";
    HAPPLY_FOLD: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "fold the local potential into the pair FFT";
    NONLOCAL_BLOCK: Kind::Choice { allowed: &["auto", "on", "off"], default: "auto" }, "auto", "auto, on, off",
        "complex (k-point) nonlocal of a band block by two GEMMs (P48C; auto: when n_pw x n_projectors x 16 B >= 256 KiB)";
    // Numerics (dft-numerics).
    ROTATE_KERNEL: Kind::Choice { allowed: &["legacy", "blocked"], default: "legacy" }, "legacy",
        "legacy, blocked", "rotation kernel of the tiled workspace";
    PURE_GEMM: Kind::Choice {
        allowed: &["nano", "packed"],
        default: if crate::MKL_BACKEND { "nano" } else { "packed" },
    }, if crate::MKL_BACKEND { "nano" } else { "packed" }, "nano, packed",
        "pure-Rust real GEMM panels: nano tiles or one packed product per panel (P48B; default packed in pure builds, nano on MKL builds where dgemm replaces these sites)";
    PACK_FIXED: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "PERCORE-KERNELS: the aarch64 pure GEMM kernel packs full eight-row strips and stores full tile columns by constant-length copies instead of one memmove call per eight doubles (bit-identical; acts on pure builds on aarch64 only)";
    MKL_DGEMM: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "rotation panels by MKL dgemm (MKL builds)";
    MKL_EVD: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "reduced symmetric eigensolve by MKL (MKL builds)";
    PAIR_FFT: Kind::Choice {
        allowed: &["grid", "pruned"],
        default: if crate::MKL_BACKEND { "pruned" } else { "grid" },
    }, if crate::MKL_BACKEND { "pruned" } else { "grid" }, "grid, pruned",
        "pair FFT plan of MKL builds (pruned: stick lanes)";
    SPIN_SHARE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "share exactly degenerate fixed-sector standard spin solves (bit-identical effective potentials and occupations; default on since 2026-11-17: dev08 si8 -18.1 %, si16 -23.8 %; mpi, strict and unequal-sector fallbacks verified clean)";
    // CLI (dft-cli).
    CUTOFF_LADDER: Kind::Override { valid: cutoff_ladder_valid }, "config", "off, on, fraction in (0, 1)",
        "overrides scf.cutoff_ladder / fraction of the input file";
    RELAX_EXTRAPOLATION: Kind::Override { valid: relax_extrapolation_valid }, "config",
        "none, first_order, second_order", "overrides the relaxation density extrapolation of the input file";
    REFERENCE_STRICT: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "ADR 0043 (2026-11-16): force the strict reference SCF path with its fixed-point certificate (accuracy gates and tests); off = the single standard path with verification";
    DGEMM_TIMERS: Kind::Presence, "off", "any value (presence)",
        "print per-site dgemm timers at exit (dgemm-timers builds)";
    PSEUDO_DIR: Kind::Path, "platform cache", "absolute path",
        "pseudopotential cache directory";
    // Measurement infrastructure (ticket BREAKTHROUGH-0-INSTRUMENT); all default off, results unchanged.
    PHASE_PROFILE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "exclusive nested phase profile per rank (init / coarse SCF / fine ordinary SCF / verification pass and fail / materialization; map phases and solve stages inside; maximum over ranks and imbalance); clocks follow BAANDI_STAGE_TIMERS (off: counts only); printed by dft-cli (results.json phase_profile, parallel.phases), m2_entry_cost and davidson_bench";
    M2_LIGHT_STATE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "BREAKTHROUGH-4: on the M2 path the replicated serial real Gamma Hamiltonian of a solve state (rebuilt at every new problem and ladder stage) is built with the minimal block capacity 2 instead of the Davidson capacity: one pair FFT plan and scratch instead of capacity / 2 full-grid plans per rank; it only supplies the ORC half sphere there and is never applied (the full operator is built if the distributed operator is unavailable); bit-identical";
    M2_LAZY_SERIAL: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "ISSUE-20: an M2 real Gamma solve state owns only the ORC half-sphere descriptor and the distributed operator; the serial real Gamma Hamiltonian (pair plans, real projector table) is built lazily only when a fallback applies it (distributed operator unavailable, threads path); supersedes M2_LIGHT_STATE on the M2 path when on; bit-identical";
    M2_LAYOUT_CACHE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "BREAKTHROUGH-4: the distributed solver reuses the row-index layout (row counts and global row lists of every rank) of its previous solve on the same thread when every rank's local rows are unchanged (one 1-double reduction decides, replicated), instead of two allgathervs (one of N doubles) and the index rebuild per solve; bit-identical";
    EWALD_CACHE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "BREAKTHROUGH-4: problem construction reuses the Ewald energy of the previous construction on the same thread when the lattice, the ion charges and the fractional positions are bit-equal (the cutoff-ladder coarse problem after the full-cutoff one); bit-identical (the same pure function of the same inputs)";
    PHASE_WAIT: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "wait probe of the phase profile (ticket BREAKTHROUGH-4-OTHER-AND-INIT): a timed barrier before every collective of the distributed solver and at the entry of every M2 solve, so that the time ranks spend waiting for the last rank (load imbalance, arrival skew) is reported apart from the collective's own transfer (`solve.wait`, `solve.entry_wait`, `coll.barrier`); changes the timing (extra barriers), never the arithmetic; diagnostic only, never an acceptance arm";
    FIRST_BLOCK: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "first-correction-block residual diagnostic (review E3): before every warm real Gamma solve a small Rayleigh-Ritz on the kept B vectors, the occupied residual block R_o = H X_o - X_o A_o, its singular values and effective ranks, and the unconverged-root counts of the Ritz-change rule against the true residual; read-only (the solve is unchanged; B extra H columns per solve when no kept H X exists)";
    SNAPSHOT_DIR: Kind::Path, "unset", "absolute directory",
        "SCF snapshot capture: every real Gamma solve of the run writes the state that replays it (input density, this rank's orbital rows, kept H X and its armed flag, tolerance, policy, layout hash, context) to this directory outside the timed map phases; replay with davidson_bench --snapshot <manifest.json>";
}

struct Resolved {
    raw: Option<OsString>,
    text: Option<String>,
    value: Value,
    source: Source,
}

/// Parse one switch from its raw environment value (`None` = unset). Pure; the
/// accessors and the default test use it.
pub fn resolve(switch: &Switch, raw: Option<&OsStr>) -> (Value, Source) {
    let text = raw.map(|r| r.to_str());
    // Legacy call sites used `std::env::var`, which fails on a non-Unicode
    // value exactly as on an unset one.
    let invalid_default = if raw.is_some() {
        Source::InvalidToDefault
    } else {
        Source::Default
    };
    match switch.kind {
        Kind::Choice { allowed, default } => match text {
            Some(Some(v)) => match allowed.iter().find(|a| **a == v) {
                Some(a) => (Value::Choice(a), Source::Env),
                None => (Value::Choice(default), Source::InvalidToDefault),
            },
            _ => (Value::Choice(default), invalid_default),
        },
        Kind::Number { default, valid } => match text {
            Some(Some(v)) => match v.parse::<f64>().ok().filter(|x| valid(*x)) {
                Some(x) => (Value::Number(Some(x)), Source::Env),
                None => (Value::Number(default), Source::InvalidToDefault),
            },
            _ => (Value::Number(default), invalid_default),
        },
        Kind::Ortho => match text {
            Some(Some(v)) => {
                if v == "dgks" {
                    return (
                        Value::Number(Some(std::f64::consts::FRAC_1_SQRT_2)),
                        Source::Env,
                    );
                }
                let ratio = v
                    .strip_prefix("dgks:")
                    .and_then(|r| r.parse::<f64>().ok())
                    .filter(|r| (1e-3..=1.0).contains(r));
                match ratio {
                    Some(r) => (Value::Number(Some(r)), Source::Env),
                    None if v == "off" => (Value::Number(None), Source::Env),
                    None => (Value::Number(None), Source::InvalidToOff),
                }
            }
            _ => (Value::Number(Some(0.1)), invalid_default),
        },
        Kind::Override { valid } => match text {
            Some(Some(v)) => (
                Value::Raw(Some(v.to_string())),
                if valid(v) {
                    Source::Env
                } else {
                    Source::InvalidToError
                },
            ),
            _ => (Value::Raw(None), invalid_default),
        },
        Kind::Presence => (
            Value::Present(raw.is_some()),
            if raw.is_some() {
                Source::Env
            } else {
                Source::Default
            },
        ),
        Kind::Path => match raw {
            Some(r) => (
                Value::Raw(Some(r.to_string_lossy().into_owned())),
                if std::path::Path::new(r).is_absolute() {
                    Source::Env
                } else {
                    Source::InvalidToError
                },
            ),
            None => (Value::Raw(None), Source::Default),
        },
    }
}

fn state() -> &'static [Resolved] {
    static STATE: OnceLock<Vec<Resolved>> = OnceLock::new();
    STATE.get_or_init(|| {
        TABLE
            .iter()
            .map(|s| {
                let raw = std::env::var_os(s.name);
                let (value, source) = resolve(s, raw.as_deref());
                let text = raw.as_ref().and_then(|r| r.to_str()).map(str::to_string);
                Resolved {
                    raw,
                    text,
                    value,
                    source,
                }
            })
            .collect()
    })
}

fn get(id: Id) -> &'static Resolved {
    &state()[id as usize]
}

/// Table entry of a switch.
pub fn switch(id: Id) -> &'static Switch {
    &TABLE[id as usize]
}

/// Effective value of a choice switch.
pub fn choice(id: Id) -> &'static str {
    match get(id).value {
        Value::Choice(v) => v,
        _ => panic!("{} is not a choice switch", switch(id).name),
    }
}

/// `true` when the effective value of an on/off switch is `on`.
pub fn on(id: Id) -> bool {
    choice(id) == "on"
}

/// Effective value of a number switch (`BAANDI_ORTHO`: DGKS ratio, `None` =
/// two passes always).
pub fn number(id: Id) -> Option<f64> {
    match get(id).value {
        Value::Number(v) => v,
        _ => panic!("{} is not a number switch", switch(id).name),
    }
}

/// The variable's text as read at startup (`None` when unset or not Unicode).
/// For override switches and legacy reports of the raw value.
pub fn text(id: Id) -> Option<&'static str> {
    get(id).text.as_deref()
}

/// The variable's raw value as read at startup (`None` when unset).
pub fn os(id: Id) -> Option<&'static OsStr> {
    get(id).raw.as_deref()
}

/// Whether the variable was set (any value).
pub fn present(id: Id) -> bool {
    get(id).raw.is_some()
}

/// The whole table with effective values and sources, in table order.
pub fn manifest() -> Vec<Entry> {
    TABLE
        .iter()
        .zip(state())
        .map(|(s, r)| Entry {
            switch: s,
            value: r.value.clone(),
            source: r.source,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn eval(id: Id, raw: Option<&str>) -> (Value, Source) {
        resolve(switch(id), raw.map(OsStr::new))
    }

    #[test]
    fn table_order_matches_ids_and_names_are_unique() {
        for (i, s) in TABLE.iter().enumerate() {
            assert_eq!(s.id as usize, i, "{}", s.name);
        }
        let mut names: Vec<_> = TABLE.iter().map(|s| s.name).collect();
        names.sort();
        names.dedup();
        assert_eq!(names.len(), TABLE.len());
    }

    #[test]
    fn legacy_parsing_rules() {
        use Source::*;
        // Default-on flags: anything but "off" is on.
        assert_eq!(eval(Id::B6, Some("off")), (Value::Choice("off"), Env));
        assert_eq!(
            eval(Id::B6, Some("yes")),
            (Value::Choice("on"), InvalidToDefault)
        );
        // Default-off flags: only "on" is on.
        assert_eq!(
            eval(Id::RITZ_FUSED, Some("1")),
            (Value::Choice("off"), InvalidToDefault)
        );
        assert_eq!(eval(Id::B5, Some("par")), (Value::Choice("par"), Env));
        assert_eq!(eval(Id::NL_ENERGY, Some("exact")).0, Value::Choice("exact"));
        assert_eq!(
            eval(Id::BUFFER_TOL, Some("1e-6")),
            (Value::Number(Some(1e-3)), InvalidToDefault)
        );
        assert_eq!(
            eval(Id::BUFFER_TOL, Some("1e-5")),
            (Value::Number(Some(1e-5)), Env)
        );
        assert_eq!(
            eval(Id::MIX_GROWTH_RESET, Some("1.0")),
            (Value::Number(None), InvalidToDefault)
        );
        assert_eq!(
            eval(Id::ORTHO, Some("dgks")).0,
            Value::Number(Some(std::f64::consts::FRAC_1_SQRT_2))
        );
        assert_eq!(
            eval(Id::ORTHO, Some("dgks:0.5")),
            (Value::Number(Some(0.5)), Env)
        );
        assert_eq!(eval(Id::ORTHO, Some("off")), (Value::Number(None), Env));
        assert_eq!(
            eval(Id::ORTHO, Some("dgks:2")),
            (Value::Number(None), InvalidToOff)
        );
        assert_eq!(
            eval(Id::LADDER_COARSE_STOP, Some("3e-3")),
            (Value::Number(Some(3e-3)), Env)
        );
        assert_eq!(
            eval(Id::LADDER_COARSE_STOP, Some("0")),
            (Value::Number(Some(1e-4)), InvalidToDefault)
        );
        assert_eq!(
            eval(Id::LADDER_COARSE_RULE, Some("atomic")),
            (Value::Choice("atomic"), Env)
        );
        assert_eq!(
            eval(Id::LADDER_FIRST_TOL, Some("coarse")),
            (Value::Choice("coarse"), Env)
        );
        assert_eq!(
            eval(Id::LADDER_FIRST_TOL, Some("last")),
            (Value::Choice("auto"), InvalidToDefault)
        );
        assert_eq!(eval(Id::CUTOFF_LADDER, Some("0.5")).1, Env);
        assert_eq!(eval(Id::CUTOFF_LADDER, Some("x")).1, InvalidToError);
        assert_eq!(eval(Id::PSEUDO_DIR, Some("relative")).1, InvalidToError);
        assert_eq!(eval(Id::DGEMM_TIMERS, Some("")).0, Value::Present(true));
    }

    #[test]
    fn mpi_entry_cost_switches_defaults() {
        use Source::*;
        for id in [Id::M2_DGEMM, Id::M2_LEAN] {
            assert_eq!(eval(id, None), (Value::Choice("on"), Default));
            assert_eq!(eval(id, Some("off")), (Value::Choice("off"), Env));
            assert_eq!(eval(id, Some("1")), (Value::Choice("on"), InvalidToDefault));
        }
        assert_eq!(eval(Id::M2_ORTHO, None).0, Value::Choice("serial"));
        assert_eq!(
            eval(Id::M2_ORTHO, Some("serial")).0,
            Value::Choice("serial")
        );
        assert_eq!(eval(Id::M2_WORKSPACE, None), (Value::Choice("on"), Default));
        assert_eq!(
            eval(Id::M2_WORKSPACE, Some("off")),
            (Value::Choice("off"), Env)
        );
        assert_eq!(
            eval(Id::M2_SELECTED_ROOTS, None),
            (Value::Choice("on"), Default)
        );
        assert_eq!(
            eval(Id::M2_SELECTED_ROOTS, Some("off")),
            (Value::Choice("off"), Env)
        );
        assert_eq!(eval(Id::M2_B5, None), (Value::Choice("off"), Default));
        assert_eq!(eval(Id::M2_B5, Some("on")), (Value::Choice("on"), Env));
        // MPI-SCALING-GAP switches (defaults mirror and gram since 0d87d6e).
        assert_eq!(
            eval(Id::M2_STICKS, None),
            (Value::Choice("mirror"), Default)
        );
        assert_eq!(
            eval(Id::M2_STICKS, Some("greedy")),
            (Value::Choice("greedy"), Env)
        );
        assert_eq!(
            eval(Id::M2_STICKS, Some("on")),
            (Value::Choice("mirror"), InvalidToDefault)
        );
        assert_eq!(eval(Id::M2_START, None), (Value::Choice("gram"), Default));
        assert_eq!(
            eval(Id::M2_START, Some("ortho")),
            (Value::Choice("ortho"), Env)
        );
        assert_eq!(
            eval(Id::M2_EIGEN, None),
            (Value::Choice("replicated"), Default)
        );
        assert_eq!(
            eval(Id::M2_EIGEN, Some("root")),
            (Value::Choice("root"), Env)
        );
        assert_eq!(
            eval(Id::M2_EIGEN, Some("on")),
            (Value::Choice("replicated"), InvalidToDefault)
        );
        assert_eq!(eval(Id::M2_SCREEN, None), (Value::Choice("dist"), Default));
        assert_eq!(
            eval(Id::M2_SCREEN, Some("replicated")),
            (Value::Choice("replicated"), Env)
        );
        assert_eq!(
            eval(Id::M2_SCREEN, Some("on")),
            (Value::Choice("dist"), InvalidToDefault)
        );
        assert_eq!(eval(Id::M2_ATOMIC, None), (Value::Choice("split"), Default));
        assert_eq!(
            eval(Id::M2_ATOMIC, Some("replicated")),
            (Value::Choice("replicated"), Env)
        );
        assert_eq!(
            eval(Id::M2_ATOMIC, Some("dist")),
            (Value::Choice("split"), InvalidToDefault)
        );
        assert_eq!(eval(Id::M2_FFT, None).0, Value::Choice("rust"));
        assert_eq!(eval(Id::M2_FFT, Some("mkl")).0, Value::Choice("mkl"));
        assert_eq!(
            eval(Id::M2_FFT, Some("fftw")),
            (Value::Choice("rust"), InvalidToDefault)
        );
    }
}
