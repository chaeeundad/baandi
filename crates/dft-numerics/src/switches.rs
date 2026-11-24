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
//!
//! RUN-CONFIG (2026-11-18): the same table is also settable from the input
//! file (`[performance]` and `[profile]` tables, key = name without `BAANDI_`
//! in lower case) with precedence environment > input > default, resolved once
//! per run into an immutable [`RunConfig`] (see the section comment there).
use std::cell::Cell;
use std::ffi::{OsStr, OsString};
use std::sync::atomic::{AtomicBool, AtomicU64, Ordering};
use std::sync::{Mutex, OnceLock};

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
    /// The input file's `[performance]` / `[profile]` table (RUN-CONFIG).
    Input,
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
            Self::Input => "input",
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
    /// The requested text: the environment value if set, else the input
    /// file's value (`None`: neither, the default applies).
    pub requested: Option<String>,
    /// Whether the switch was read under the configuration (by the time the
    /// manifest was taken).
    pub consulted: bool,
    /// The path the switch actually took and the fallback reason, when the
    /// code that acted on it (or the result writer) recorded one.
    pub note: Option<PathNote>,
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
    INIT_OWNED_STATES: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "owned Gamma projector rows and fixed-chunk atomic-start projections (rounding class)";
    CERT_SCREEN: Kind::Choice { allowed: &["on", "off", "trace"], default: "on" }, "on", "on, off, trace",
        "certificate screen (trace prints the estimates)";
    POTENTIAL_REUSE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "certificate solve reuses the screen's Hartree/XC fields";
    MIX_GROWTH_RESET: Kind::Number { default: None, valid: growth_reset_valid }, "none", "number > 1",
        "reset the mixing history when the residual grows by more than this factor";
    MIX_CELL_KERKER: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "G-MIXER-DIVERGENCE: replace plain Pulay by Kerker-Pulay only when the density basis has wavelengths longer than the equal-volume cube: q0^2 = max(0, (2*pi/volume^(1/3))^2 - min_nonzero_G^2); beta/history unchanged, charge mode removed, acceptance unchanged; replicated mixer (including on MPI), not a physical dielectric model";
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
    LADDER_FRACTION_RULE: Kind::Choice { allowed: &["off", "grid"], default: "off" }, "off",
        "off, grid", "experimental cutoff choice: fill the density FFT grid plateau of the requested fraction without increasing any grid dimension; no material constants";
    LADDER_STOP_RULE: Kind::Choice { allowed: &["off", "modes"], default: "off" }, "off",
        "off, modes", "experimental coarse residual stop: preserve the half-cutoff RMS budget per nonzero density Fourier mode (Parseval); never loosen the requested residual";
    LADDER_WATCH: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "diagnose early fine error above the contraction forecast from maps 1 and 2; record a path note after completion, never redo or discard work";
    LADDER_BASIS_CHECK: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "plan-time coarse basis/band/subspace check: skip an incompatible cutoff ladder before constructing the coarse problem, preserving the ordinary full-cutoff result";
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
    ATOMIC_LOCAL_BLOCK: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "LANE-INIT: guarded CholQR2 on the atomic starting span when no multi-rank communicator is active; reuses adopted atomic block algebra, falls back to original MGS, leaves MPI and Davidson corrections unchanged; rounding class, default off pending dev08";
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
    // Adopted 2026-11-24 (dev08 MPI8 A/B at 0e8b08da; attempts ledger "P3 timed").
    REF_DAV_RESTART_CORRECTIONS: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "reference Davidson only: form the first post-contraction corrections directly from retained V/HV, omitting identity-coefficient GEMMs; measured H/S refresh and band decisions unchanged";
    // Adopted 2026-11-24 (dev08 MPI8 A/B at 66d45e6: si64 -1.8%, si128 -3.4%, others neutral).
    REFSTART_MKL: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "reference SCF atomic start only: Gram and rotation products through the dispatching GEMM entries (MKL panels) instead of the direct pure-Rust kernel; rounding class";
    // Adopted 2026-11-24 (dev08 MPI8 A/B at 0e8b08da; attempts ledger "P3 timed").
    REF_DAV_WIDE_ROTATE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "reference Davidson only: route the V C / HV C rotations through the wide GEMM entry (MKL panels) instead of the whole-output-bound pure-Rust fallback; rounding class";
    REFSCF_VERIFY_REUSE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "reference SCF only: consume a failed verification map as the next ordinary Broyden iterate; passing criterion unchanged";
    // Adopted 2026-11-24 (dev08 MPI8 A/B at 0ce77fdc: si64 -5.4%, si128 -6.5%, elongated/rskew neutral; attempts ledger).
    REFSCF_VERIFY_GATE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "reference SCF only: defer predicted failing verification attempts using observed error/estimator scaling; never deliver without a passing attempt";
    REFSCF_VERIFY_GATE_FIRST: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "reference SCF only, requires REFSCF_VERIFY_GATE: predict the first verification from density/residual secants; undefined predictions attempt and the passing criterion is unchanged";
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
    M2_ATOMIC_BLOCK: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "INIT-SCALE: with M2_ATOMIC=split, every rank fills its contiguous block of span rows (same elementwise arithmetic) and the span is orthonormalized by a row-distributed CholQR2 (local Gram GEMM, one tree_sum, replicated Cholesky and triangular inverse, local update; twice; one allgatherv of the rows) instead of the column-cyclic two-pass Gram-Schmidt whose second pass is a serial chain over the columns; falls back to the split Gram-Schmidt when a pivot is below the dependence guard or data is nonfinite; rounding class (span and trajectory change)";
    // Replicated grid stages outside the solve (ticket MPI-REMAINING, docs/performance/mpi-remaining.md); default off.
    STOP_METRIC: Kind::Choice { allowed: &["fft", "kept"], default: "kept" }, "kept", "fft, kept",
        "standard stop-rule Coulomb estimate E_H[rho_out - rho_in]: fft (transform the residual forward and back on the full grid) or kept (Parseval sum over the map's own Fourier coefficients of rho_in, kept from the Hartree potential, and rho_out, kept from the energy audit; no FFT; falls back to fft when a density was clamped or the kept coefficients do not belong to the two densities; rounding class of the estimate only)";
    M2_GRID: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "M2-GRID: on the M2 path (2+ ranks, LDA PZ/PW, no core density, density floor or symmetry) the full-grid stages outside the solve run on a slab-decomposed full-grid FFT (owned real-space planes, owned sticks, two transposes per transform): the input density's Hartree FFT pair, the LDA XC and the effective potential on the planes with one allgatherv of the effective potential and the partial sums; the energy audit's output-density Hartree FFT and grid sums on the planes with one allgatherv of partial sums; the stop-rule Coulomb estimate from the kept owned coefficients; every sum in a canonical plane/stick order, so results do not depend on the rank count; rounding class against the replicated path";
    M2_MIX: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "M2-GRID: on the M2 path (2+ ranks) plain Pulay mixing (no Hartree metric, no Kerker) keeps its window of input densities, residuals and differences on the rank's block of fixed 4096-point chunks only: the newest difference, its Gram row and the right-hand side on the block with one allgatherv of per-chunk partial sums added in chunk order (rank-count independent), the replicated small eigenproblem, the update and the two charge corrections on the block, one allgatherv of the new density; exact checkpoints reject the mixer (opaque state); rounding class against the replicated mixer";
    M2_LOCAL_FIELDS: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "M2-OWNED (default on since 2026-11-18, dev08 si32/si64 np=8 wall -2.5 to -6.0 %; declines automatically for active-boundary consumers): with M2_GRID on the M2 standard path, the effective potential stays on the producing ranks' slab planes (the planes the distributed H apply multiplies; redistributed by one alltoall only when the slab and operator plane axes differ) and is fed to the distributed operator and its kept-start potential change (delta V) without the full-grid allgatherv; the complex Hamiltonian's full potential is gathered only for a consumer that needs it (cold atomic start, complex or M1 fallback), tagged with an epoch so a stale potential is rejected; bit-identical to M2_GRID alone";
    M2_LOCAL_DENSITY: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "M2-OWNED (default on since 2026-11-18, dev08 si32/si64 np=8 wall -2.5 to -6.0 %; declines automatically for active-boundary consumers): with M2_LOCAL_FIELDS and M2_MIX on, the output density stays on the owned slab planes through the energy audit, the stop rule, the session metrics (residual, charge, minimum: canonical per-plane sums), the distributed Pulay mixer (window on the owned planes) and the next potential; the full density is gathered only for a consumer that needs it (delivered result, checkpoint, ladder hand-off, declined map); rounding class against M2_MIX alone (per-plane instead of per-chunk canonical sums)";
    M2_OWNED_HANDOFF: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "G-HANDOFF-OWNED: experimental direct coarse-to-fine owned ORC transfer by Miller index; unchanged full-matrix fallback for unsupported consumers, default off pending dev08 measurement";
    // Solvers (dft-solvers).
    MIX_COULOMB_GUARD: Kind::Choice { allowed: &["off", "trace", "on"], default: "off" }, "off", "off, trace, on",
        "SCF-TRAJECTORY: plain Pulay (L2 fit) measures the residual its least squares predicts, R - dR gamma, in the Coulomb (Hartree-energy) metric; on: when that predicted residual is larger than the current one in the Coulomb metric, this step's coefficients are refitted in the Coulomb metric (no tuned constant: the rule is 'the fit must not predict a larger long-wavelength error'); trace: compute and record the ratio, never act; replicated PulayMixer and the M2 distributed mixer (chunk and owned-plane modes, Coulomb forms on the owned sticks; MIX-COULOMB-GUARD-DIST); not kerker-pulay, not mixing_metric = hartree";
    VERIFY_SCHEDULE: Kind::Choice { allowed: &["off", "trace", "on", "first", "first-trace"], default: "on" }, "on", "off, trace, on, first, first-trace",
        "F-VERIFY-SCHEDULE: geometric retry prediction from existing Coulomb errors; off preserves the original path, trace records predictions, on defers predicted failures with a bounded retry horizon; certificate and stop threshold unchanged; first adds an experimental two-map-bounded first-attempt lag/contraction forecast, first-trace records that forecast without acting on it (both retain on retry scheduling); default remains on";
    STOP_PREDICTOR: Kind::Choice { allowed: &["off", "trace", "on"], default: "off" }, "off", "off, trace, on",
        "SCF-TAIL: second estimate of the ADR 0043 fixed-point verification from the two maps it already has, -X = -E_H-form <rho_out - rho_in, F(rho_out) - rho_out> (one extra Coulomb-metric evaluation, no extra solve); in the linear Hartree model -X / (true energy error of the delivered state) = 1 + A >= 1, while the old criterion E_H[F(rho_out) - rho_out] reads A (1 + A) (A = Coulomb-weighted static response); measured on the Mac set -X / T >= 0.654, so it is not a bound; trace: print both and the verification energy change, never act; on: an attempt whose old criterion fails also passes when 0 < -X < threshold / 4 (safety factor from the held-out calibration), the stop estimate is below the threshold and the residual guard passes; replicated densities only (the owned-plane loop keeps the old criterion)";
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
    DAVIDSON_PACKED_REDUCE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "Pack only consumed upper-border entries and residual norms for distributed lean reductions; GEMMs unchanged, default off.";
    DAVIDSON_ROTATE_COLUMNS: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "Serial real Davidson Ritz pairs through existing fixed column-group GEMMs; rounding class, default off.";
    DAVIDSON_UPPER_BORDER: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "Compute only retained upper-border tiles of real Davidson Gram products; rounding class, default off.";
    DAVIDSON_PROFILE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "Bounded per-call Davidson lane trace with shapes and rank/map labels; diagnostic only.";
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
    BINARY_HASH_CACHE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "LANE-LAUNCHIO: cache executable SHA-256 in a metadata-keyed sibling sidecar on Unix; preserves provenance, falls back to fresh hashing; trusted build directories only";
    CUTOFF_LADDER: Kind::Override { valid: cutoff_ladder_valid }, "config", "off, on, fraction in (0, 1)",
        "overrides scf.cutoff_ladder / fraction of the input file";
    RELAX_EXTRAPOLATION: Kind::Override { valid: relax_extrapolation_valid }, "config",
        "none, first_order, second_order", "overrides the relaxation density extrapolation of the input file";
    REFERENCE_SCF: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "issue #25 independent reciprocal-density SCF driver (opt-in; patch 1 verification pending)";
    REFERENCE_STRICT: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "ADR 0043 (2026-11-16): force the strict reference SCF path with its fixed-point certificate (accuracy gates and tests); off = the single standard path with verification";
    DGEMM_TIMERS: Kind::Presence, "off", "any value (presence)",
        "print per-site dgemm timers at exit (dgemm-timers builds)";
    PSEUDO_DIR: Kind::Path, "platform cache", "absolute path",
        "pseudopotential cache directory";
    FFT_KEPT_FUSE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "LANE-FFT: fuse M2 kept-plane copy with potential multiplication; bit-identical, no separate kept_copy pass";
    // Measurement infrastructure (ticket BREAKTHROUGH-0-INSTRUMENT); all default off, results unchanged.
    COLLECTIVE_TRACE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "bounded per-rank collective arrivals in maps; requires PHASE_PROFILE, clocks follow STAGE_TIMERS";
    FFT_PROFILE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "LANE-FFT Phase 0: detailed FFT stages and logical traffic under PHASE_PROFILE; clocks follow STAGE_TIMERS; no numerical change";
    SERIAL_DENSITY_FUSED: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "LANE-GLUE: serial one-core standard Gamma density sink; preserves band addition order, no MPI dispatch";
    GLUE_PROFILE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "LANE-GLUE Phase 0: rank-local grid and mixer spans/work; requires PHASE_PROFILE; clocks follow STAGE_TIMERS; no numerical change";
    PHASE_PROFILE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "exclusive nested phase profile per rank (init / coarse SCF / fine ordinary SCF / verification pass and fail / materialization; map phases and solve stages inside; maximum over ranks and imbalance); clocks follow BAANDI_STAGE_TIMERS (off: counts only); printed by dft-cli (results.json phase_profile, parallel.phases), m2_entry_cost and davidson_bench";
    M2_LIGHT_STATE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "BREAKTHROUGH-4: on the M2 path the replicated serial real Gamma Hamiltonian of a solve state (rebuilt at every new problem and ladder stage) is built with the minimal block capacity 2 instead of the Davidson capacity: one pair FFT plan and scratch instead of capacity / 2 full-grid plans per rank; it only supplies the ORC half sphere there and is never applied (the full operator is built if the distributed operator is unavailable); bit-identical";
    M2_OWNED_INPUT: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "G-OWNED-ROWS: pack a replicated complex warm start directly into owned ORC rows for the standard M2 workspace path; retain full-column validation and coefficient arithmetic; other paths retain global packing; bit-identical, default off pending dev08";
    M2_LAZY_SERIAL: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "ISSUE-20: an M2 real Gamma solve state owns only the ORC half-sphere descriptor and the distributed operator; the serial real Gamma Hamiltonian (pair plans, real projector table) is built lazily only when a fallback applies it (distributed operator unavailable, threads path); supersedes M2_LIGHT_STATE on the M2 path when on; bit-identical";
    M2_LAYOUT_CACHE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "BREAKTHROUGH-4: the distributed solver reuses the row-index layout (row counts and global row lists of every rank) of its previous solve on the same thread when every rank's local rows are unchanged (one 1-double reduction decides, replicated), instead of two allgathervs (one of N doubles) and the index rebuild per solve; bit-identical";
    EWALD_CACHE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "BREAKTHROUGH-4: problem construction reuses the Ewald energy of the previous construction on the same thread when the lattice, the ion charges and the fractional positions are bit-equal (the cutoff-ladder coarse problem after the full-cutoff one); bit-identical (the same pure function of the same inputs)";
    EWALD_FAST: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "INIT-EWALD-FAST: balanced Ewald scheme (dft-physics ewald::Scheme::Balanced): splitting parameter from the real/reciprocal cost balance for the cell volume and ion count, cutoffs from rigorous tail bounds (truncation <= min(tolerance (sum |Z|)^2 / Omega^(1/3), tolerance * N) Ha), ion pairs i <= j and the reciprocal half space with per-ion phase tables; energy, forces and stress share the parameters; rounding class (different truncation and summation, error vs a tight reference far below 1e-12 Ha/atom); off = reference scheme (eta sqrt(pi)/Omega^(1/3), bit-identical to the previous evaluator)";
    IONIC_FAST: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "INIT-IONIC-FAST: the structure factor of the local ionic potential sum_I exp(-i G.tau_I) from per-ion, per-axis phase tables exp(-i m_k b_k.tau_I) (two complex products per ion and G) instead of one sin/cos pair per ion and G; falls back to the direct phase for any G whose Cartesian vector is not the Miller combination of the unit-Miller vectors; rounding class (phase rounding, about 1e-13 relative)";
    INIT_PROFILE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "LANE-INIT phase 0: rank-local startup call/work counters and nested timers; STAGE_TIMERS=off retains counts without clocks; no numerical changes, default off";
    PHASE_VALIDATE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "G-PHASE-VALIDATE: owned projector finite-dot validation uses a sufficient numerical range proof (finite momentum and position components with maxima <= 2^500); otherwise retains the original pair loop and errors; coefficients and lazy materialization unchanged; off until dev08 A/B";
    IMPORT_CHECK: Kind::Choice { allowed: &["split", "serial"], default: "split" }, "split", "split, serial",
        "LADDER-IMPORT-SPLIT: orthonormality check of an imported warm orbital block (cutoff-ladder and direct-minimization hand-off, `initialize_electronic`): split computes the same per-entry overlaps (each the sequential row sum, bit-identical) eight columns in lockstep and, with a registered communicator, only for this rank's columns j (j mod P), with one 1-double tree sum of the failure flags; serial is the former replicated scalar loop on every rank; same verdict, results bit-identical";
    M2_LAZY_RESULT: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "LAZY-RESULT: when the caller declared that it reads no delivered orbitals (dft-cli run without forces, stress, ELF or partial-charge output, classic mode), an M2 result keeps its lazy distributed orbital blocks instead of gathering the full N x B block on every rank at delivery; ScfResult::materialize_orbitals gathers them on demand (collective); energy, density, eigenvalues and history unchanged; off = always gather";
    PHASE_WAIT: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "wait probe of the phase profile (ticket BREAKTHROUGH-4-OTHER-AND-INIT): a timed barrier before every collective of the distributed solver and at the entry of every M2 solve, so that the time ranks spend waiting for the last rank (load imbalance, arrival skew) is reported apart from the collective's own transfer (`solve.wait`, `solve.entry_wait`, `coll.barrier`); changes the timing (extra barriers), never the arithmetic; diagnostic only, never an acceptance arm";
    VERIFY_PROFILE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "LANE-VERIFY per-attempt/per-rank size and phase records; requires PHASE_PROFILE, clocks only with STAGE_TIMERS; no numerical changes";
    FIRST_BLOCK: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "first-correction-block residual diagnostic (review E3): before every warm real Gamma solve a small Rayleigh-Ritz on the kept B vectors, the occupied residual block R_o = H X_o - X_o A_o, its singular values and effective ranks, and the unconverged-root counts of the Ritz-change rule against the true residual; read-only (the solve is unchanged; B extra H columns per solve when no kept H X exists)";
    SNAPSHOT_DIR: Kind::Path, "unset", "absolute directory",
        "SCF snapshot capture: every real Gamma solve of the run writes the state that replays it (input density, this rank's orbital rows, kept H X and its armed flag, tolerance, policy, layout hash, context) to this directory outside the timed map phases; replay with davidson_bench --snapshot <manifest.json>";
    // Adopted 2026-11-24 (dev08 MPI8 A/B at 0e8b08da; attempts ledger "P3 timed").
    REFSCF_SETUP_PHASE_REUSE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "P3-SETUP: reference start reuses identical direct structure factors across atomic charge/local potential and magnetic orbital channels; same arithmetic, default off pending dev08 A/B";
    // Provisionally adopted 2026-11-25 (development protocol: screen + partial P3d si128 -1.2%); precise measurement pending.
    REF_DAV_SINGLE_GEMM: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "reference Davidson only: transposed borders and repair projections use one sequential MKL GEMM over all local rows; rounding class, pure builds retain panel arithmetic";
    FFT_SPARSE_CLEAR: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P3-FFT: clear only holes before coefficient/middle/slab scatter using cached contiguous ranges; same FFT order, default off pending dev08";
    FFT_DIRECT_PACK: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P3-FFT: pack distributed real H input ORC rows directly into FFT sticks, avoiding the intermediate coefficient write/read; default off pending dev08";
    // Default off since 2026-11-25: dev08 one-run screen, energy change ~0.85x the Coulomb error, no earlier pass (ADR 0044 status).
    REFSCF_VERIFY_ENERGY: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "reference SCF only: accept the independent verification map when the variational energy change is below 5e-9 Ha; Coulomb/residual diagnostics retained, predictor gates inactive; off restores the old fixed-point rule";
    FFT_FUSED_ROWS: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P3-FFT2: stream distributed local-potential last-axis inverse/multiply/forward through reusable rows with immutable zero-hole input; rustfft lanes only, rounding class, default off pending dev08";
}

struct Resolved {
    raw: Option<OsString>,
    text: Option<String>,
    /// The input file's text for this switch (`[performance]` / `[profile]`).
    input: Option<String>,
    value: Value,
    source: Source,
}

/// Input-file table a switch belongs to (RUN-CONFIG, 2026-11-18).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Table {
    /// `[performance]`: numerical and path switches.
    Performance,
    /// `[profile]`: measurement switches (clocks, counters, captures); they
    /// never change the arithmetic.
    Profile,
    /// Not settable from the input file (environment only).
    EnvOnly,
}
impl Table {
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Performance => "performance",
            Self::Profile => "profile",
            Self::EnvOnly => "environment only",
        }
    }
}

/// Table of a switch in the input file.
pub fn table_of(id: Id) -> Table {
    use Id::*;
    match id {
        COLLECTIVE_TRACE | GLUE_PROFILE | FFT_PROFILE | PHASE_PROFILE | INIT_PROFILE
        | STAGE_TIMERS | DAVIDSON_STAGES | DAVIDSON_PROFILE | SNAPSHOT_DIR | PHASE_WAIT
        | FIRST_BLOCK | DGEMM_TIMERS | VERIFY_PROFILE => Table::Profile,
        // Overrides of input-file values (`scf.cutoff_ladder`, the relaxation
        // extrapolation) and the pseudopotential cache location.
        CUTOFF_LADDER | RELAX_EXTRAPOLATION | PSEUDO_DIR => Table::EnvOnly,
        _ => Table::Performance,
    }
}

/// Input-file key of a switch: the registry name without `BAANDI_`, lower case
/// (`BAANDI_M2_GRID` -> `m2_grid`).
pub fn input_key(id: Id) -> String {
    switch(id)
        .name
        .trim_start_matches("BAANDI_")
        .to_ascii_lowercase()
}

/// Whether a switch takes part in the cross-rank agreement hash: every switch
/// except locations and the exit-time dgemm report (a rank may differ there
/// without changing what the ranks compute or which collectives they enter).
pub fn in_agreement(id: Id) -> bool {
    !matches!(id, Id::PSEUDO_DIR | Id::SNAPSHOT_DIR | Id::DGEMM_TIMERS)
}

/// Parse an input-file value. Strict: a value the environment rule would
/// replace by the default is an input error here.
fn resolve_input(switch: &Switch, text: &str) -> Result<Value, String> {
    let bad = || {
        format!(
            "[{}] {} = {text:?}: allowed values {}",
            table_of(switch.id).as_str(),
            input_key(switch.id),
            switch.allowed
        )
    };
    match switch.kind {
        Kind::Choice { allowed, .. } => allowed
            .iter()
            .find(|a| **a == text)
            .map(|a| Value::Choice(a))
            .ok_or_else(bad),
        Kind::Number { valid, .. } => text
            .parse::<f64>()
            .ok()
            .filter(|x| valid(*x))
            .map(|x| Value::Number(Some(x)))
            .ok_or_else(bad),
        Kind::Ortho => match resolve(switch, Some(OsStr::new(text))) {
            (v, Source::Env) => Ok(v),
            _ => Err(bad()),
        },
        Kind::Override { valid } => valid(text)
            .then(|| Value::Raw(Some(text.to_string())))
            .ok_or_else(bad),
        Kind::Presence => match text {
            "on" => Ok(Value::Present(true)),
            "off" => Ok(Value::Present(false)),
            _ => Err(bad()),
        },
        Kind::Path => std::path::Path::new(text)
            .is_absolute()
            .then(|| Value::Raw(Some(text.to_string())))
            .ok_or_else(bad),
    }
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

// ---------------------------------------------------------------------------
// Run configuration (RUN-CONFIG, 2026-11-18).
//
// A `RunConfig` is the immutable resolution of every switch for one run:
// environment > input file (`[performance]` / `[profile]` tables) > registry
// default. It is leaked once (`&'static`, Copy, no reference count), so values
// and texts can be handed out with the `'static` lifetime the legacy accessors
// return. Besides the values it carries the run's mutable bookkeeping, all
// cheap relaxed atomics or a lock taken outside hot loops: which switches were
// consulted, the path a switch actually took and its fallback reason, the run
// ID and the run-local counters.
//
// The legacy free functions (`choice`, `on`, ...) read the CURRENT
// configuration: the thread's scope (`RunConfig::enter`, propagated by
// `pool::WorkerPool` into its workers) if any, else the process configuration
// (`install_process`, or the environment alone at the first read). New code
// resolves its switches from an explicit `RunConfig` at session or plan
// creation and never looks one up inside a hot loop.
// ---------------------------------------------------------------------------

/// Run-local counters (cheap relaxed atomics in the run's configuration, so
/// two runs in one process count independently). Process-global counters of
/// the same events are kept as before for the existing result blocks.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RunCounter {
    /// Stop-rule Coulomb estimates from kept Fourier coefficients.
    StopKept = 0,
    /// Stop-rule Coulomb estimates that transformed the residual.
    StopFft = 1,
    /// Ewald energies reused from the previous construction.
    EwaldHits = 2,
    /// Ewald energies computed.
    EwaldComputed = 3,
    /// Problems constructed under this configuration.
    Problems = 4,
    /// SCF session steps.
    SessionSteps = 5,
    /// Fixed-point verification maps (passed or failed).
    Verifications = 6,
    /// Lazy complex projector table constructions.
    StatesFullBuilds = 7,
    /// Owned complex projector bytes constructed over this run.
    StatesTableBytes = 8,
    /// Equivalent replicated complex table bytes.
    StatesReplicatedBytes = 9,
    /// Atom-independent reconstruction factor bytes.
    StatesRecipeBytes = 10,
    /// Owned constructors reaching phase validation setup (per rank/run).
    PhaseValidateConstructors = 11,
    /// Species validation calls, excluding factor/norm checks.
    PhaseValidateCalls = 12,
    /// Actual finite-dot pair checks, including the failing pair.
    PhaseValidatePairs = 13,
    /// Species accepted by the sufficient range proof.
    PhaseValidateFast = 14,
    /// Species taking the original loop after an enabled proof declines.
    PhaseValidateFallbacks = 15,
    /// Exclusive validation nanoseconds, only with PHASE_PROFILE on.
    PhaseValidateNanos = 16,
    /// INIT lane: preflight.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitPreflightCalls = 17,
    /// INIT lane: preflight.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitPreflightNanos = 18,
    /// INIT lane: structure.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitStructureCalls = 19,
    /// INIT lane: structure.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitStructureNanos = 20,
    /// INIT lane: upf.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitUpfCalls = 21,
    /// INIT lane: upf.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitUpfNanos = 22,
    /// INIT lane: prepare.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitPrepareCalls = 23,
    /// INIT lane: prepare.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitPrepareNanos = 24,
    /// INIT lane: check_basis.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitCheckBasisCalls = 25,
    /// INIT lane: check_basis.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitCheckBasisNanos = 26,
    /// INIT lane: check_grid.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitCheckGridCalls = 27,
    /// INIT lane: check_grid.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitCheckGridNanos = 28,
    /// INIT lane: projectors.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitProjectorsCalls = 29,
    /// INIT lane: projectors.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitProjectorsNanos = 30,
    /// INIT lane: factors.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitFactorsCalls = 31,
    /// INIT lane: factors.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitFactorsNanos = 32,
    /// INIT lane: projector_assembly.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitProjectorAssemblyCalls = 33,
    /// INIT lane: projector_assembly.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitProjectorAssemblyNanos = 34,
    /// INIT lane: owned_rows_build.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitOwnedRowsBuildCalls = 35,
    /// INIT lane: owned_rows_build.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitOwnedRowsBuildNanos = 36,
    /// INIT lane: ionic.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitIonicCalls = 37,
    /// INIT lane: ionic.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitIonicNanos = 38,
    /// INIT lane: phase_tables.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitPhaseTablesCalls = 39,
    /// INIT lane: phase_tables.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitPhaseTablesNanos = 40,
    /// INIT lane: problem.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitProblemCalls = 41,
    /// INIT lane: problem.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitProblemNanos = 42,
    /// INIT lane: state_grid.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitStateGridCalls = 43,
    /// INIT lane: state_grid.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitStateGridNanos = 44,
    /// INIT lane: hamiltonian.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitHamiltonianCalls = 45,
    /// INIT lane: hamiltonian.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitHamiltonianNanos = 46,
    /// INIT lane: band_density.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitBandDensityCalls = 47,
    /// INIT lane: band_density.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitBandDensityNanos = 48,
    /// INIT lane: core_density.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitCoreDensityCalls = 49,
    /// INIT lane: core_density.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitCoreDensityNanos = 50,
    /// INIT lane: start_density.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitStartDensityCalls = 51,
    /// INIT lane: start_density.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitStartDensityNanos = 52,
    /// INIT lane: ionic_fft.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitIonicFftCalls = 53,
    /// INIT lane: ionic_fft.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitIonicFftNanos = 54,
    /// INIT lane: worker_ffts.calls; cumulative per rank/run, only with INIT_PROFILE.
    InitWorkerFftsCalls = 55,
    /// INIT lane: worker_ffts.ns; cumulative per rank/run, only with INIT_PROFILE.
    InitWorkerFftsNanos = 56,
    /// INIT lane: preflight.wave_rows; cumulative per rank/run, only with INIT_PROFILE.
    InitPreflightWaveRows = 57,
    /// INIT lane: preflight.density_rows; cumulative per rank/run, only with INIT_PROFILE.
    InitPreflightDensityRows = 58,
    /// INIT lane: preflight.grid_points; cumulative per rank/run, only with INIT_PROFILE.
    InitPreflightGridPoints = 59,
    /// INIT lane: factors.rows; cumulative per rank/run, only with INIT_PROFILE.
    InitFactorRows = 60,
    /// INIT lane: factors.shells; cumulative per rank/run, only with INIT_PROFILE.
    InitFactorShells = 61,
    /// INIT lane: factors.radial_transforms; cumulative per rank/run, only with INIT_PROFILE.
    InitRadialTransforms = 62,
    /// INIT lane: factors.radial_samples; cumulative per rank/run, only with INIT_PROFILE.
    InitRadialSamples = 63,
    /// INIT lane: factors.harmonics; cumulative per rank/run, only with INIT_PROFILE.
    InitHarmonics = 64,
    /// INIT lane: factors.bytes; cumulative per rank/run, only with INIT_PROFILE.
    InitFactorBytes = 65,
    /// INIT lane: projectors.phase_pairs; cumulative per rank/run, only with INIT_PROFILE.
    InitProjectorPhasePairs = 66,
    /// INIT lane: projectors.materialized_bytes; cumulative per rank/run, only with INIT_PROFILE.
    InitProjectorBytes = 67,
    /// INIT lane: projectors.owned_calls; cumulative per rank/run, only with INIT_PROFILE.
    InitOwnedConstructors = 68,
    /// INIT lane: ionic.rows; cumulative per rank/run, only with INIT_PROFILE.
    InitIonicRows = 69,
    /// INIT lane: ionic.shells; cumulative per rank/run, only with INIT_PROFILE.
    InitIonicShells = 70,
    /// INIT lane: ionic.radial_samples; cumulative per rank/run, only with INIT_PROFILE.
    InitIonicRadialSamples = 71,
    /// INIT lane: ionic.phase_pairs; cumulative per rank/run, only with INIT_PROFILE.
    InitIonicPhasePairs = 72,
    /// INIT lane: ionic.table_species; cumulative per rank/run, only with INIT_PROFILE.
    InitIonicTableSpecies = 73,
    /// INIT lane: problem.atoms; cumulative per rank/run, only with INIT_PROFILE.
    InitAtoms = 74,
    /// INIT lane: problem.bands; cumulative per rank/run, only with INIT_PROFILE.
    InitBands = 75,
    /// INIT lane: problem.kpoints; cumulative per rank/run, only with INIT_PROFILE.
    InitKpoints = 76,
    /// INIT lane: problem.grid_points; cumulative per rank/run, only with INIT_PROFILE.
    InitGridPoints = 77,
    /// INIT lane: problem.wave_rows; cumulative per rank/run, only with INIT_PROFILE.
    InitWaveRows = 78,
    /// INIT lane: problem.density_rows; cumulative per rank/run, only with INIT_PROFILE.
    InitDensityRows = 79,
    /// INIT Phase 2: atomic.span.calls, local to the rank/run.
    InitAtomicSpanCalls = 80,
    /// INIT Phase 2: atomic.span.ns, local to the rank/run.
    InitAtomicSpanNanos = 81,
    /// INIT Phase 2: atomic.mgs.calls, local to the rank/run.
    InitAtomicMgsCalls = 82,
    /// INIT Phase 2: atomic.mgs.ns, local to the rank/run.
    InitAtomicMgsNanos = 83,
    /// INIT Phase 2: atomic.block_fill.calls, local to the rank/run.
    InitAtomicBlockFillCalls = 84,
    /// INIT Phase 2: atomic.block_fill.ns, local to the rank/run.
    InitAtomicBlockFillNanos = 85,
    /// INIT Phase 2: atomic.block_gram.calls, local to the rank/run.
    InitAtomicBlockGramCalls = 86,
    /// INIT Phase 2: atomic.block_gram.ns, local to the rank/run.
    InitAtomicBlockGramNanos = 87,
    /// INIT Phase 2: atomic.block_factor.calls, local to the rank/run.
    InitAtomicBlockFactorCalls = 88,
    /// INIT Phase 2: atomic.block_factor.ns, local to the rank/run.
    InitAtomicBlockFactorNanos = 89,
    /// INIT Phase 2: atomic.block_update.calls, local to the rank/run.
    InitAtomicBlockUpdateCalls = 90,
    /// INIT Phase 2: atomic.block_update.ns, local to the rank/run.
    InitAtomicBlockUpdateNanos = 91,
    /// INIT Phase 2: atomic.block_gather.calls, local to the rank/run.
    InitAtomicBlockGatherCalls = 92,
    /// INIT Phase 2: atomic.block_gather.ns, local to the rank/run.
    InitAtomicBlockGatherNanos = 93,
    /// INIT Phase 2: atomic.rows, local to the rank/run.
    InitAtomicRows = 94,
    /// INIT Phase 2: atomic.columns, local to the rank/run.
    InitAtomicColumns = 95,
    /// INIT Phase 2: atomic.input_bytes, local to the rank/run.
    InitAtomicInputBytes = 96,
    /// INIT Phase 2: atomic.mgs_projection_macs_upper_bound, local to the rank/run.
    InitAtomicMgsMacs = 97,
    /// INIT Phase 2: atomic.local_attempts, local to the rank/run.
    InitAtomicLocalAttempts = 98,
    /// INIT Phase 2: atomic.local_accepted, local to the rank/run.
    InitAtomicLocalAccepted = 99,
    /// INIT Phase 2: atomic.local_fallbacks, local to the rank/run.
    InitAtomicLocalFallbacks = 100,
    /// INIT Phase 2: atomic.block_gram_flops, local to the rank/run.
    InitAtomicBlockGramFlops = 101,
    /// INIT Phase 2: atomic.block_update_flops, local to the rank/run.
    InitAtomicBlockUpdateFlops = 102,
    /// Sum of FFT grid points over atomic span calls.
    InitAtomicGridPoints = 103,
    /// Sum of first FFT dimensions over atomic span calls.
    InitAtomicGridX = 104,
    /// Sum of second FFT dimensions over atomic span calls.
    InitAtomicGridY = 105,
    /// Sum of third FFT dimensions over atomic span calls.
    InitAtomicGridZ = 106,
    /// Sum of requested bands over atomic span calls.
    InitAtomicBands = 107,
    /// GLUE diagnostic: grid.potential.calls; local to this run/rank.
    GlueGridPotentialCalls = 108,
    /// GLUE diagnostic: grid.potential.ns; local to this run/rank.
    GlueGridPotentialNanos = 109,
    /// GLUE diagnostic: grid.audit.calls; local to this run/rank.
    GlueGridAuditCalls = 110,
    /// GLUE diagnostic: grid.audit.ns; local to this run/rank.
    GlueGridAuditNanos = 111,
    /// GLUE diagnostic: grid.forward.calls; local to this run/rank.
    GlueGridForwardCalls = 112,
    /// GLUE diagnostic: grid.forward.ns; local to this run/rank.
    GlueGridForwardNanos = 113,
    /// GLUE diagnostic: grid.hartree.calls; local to this run/rank.
    GlueGridHartreeCalls = 114,
    /// GLUE diagnostic: grid.hartree.ns; local to this run/rank.
    GlueGridHartreeNanos = 115,
    /// GLUE diagnostic: grid.inverse.calls; local to this run/rank.
    GlueGridInverseCalls = 116,
    /// GLUE diagnostic: grid.inverse.ns; local to this run/rank.
    GlueGridInverseNanos = 117,
    /// GLUE diagnostic: grid.xc.calls; local to this run/rank.
    GlueGridXcCalls = 118,
    /// GLUE diagnostic: grid.xc.ns; local to this run/rank.
    GlueGridXcNanos = 119,
    /// GLUE diagnostic: grid.assemble.calls; local to this run/rank.
    GlueGridAssembleCalls = 120,
    /// GLUE diagnostic: grid.assemble.ns; local to this run/rank.
    GlueGridAssembleNanos = 121,
    /// GLUE diagnostic: grid.pack.calls; local to this run/rank.
    GlueGridPackCalls = 122,
    /// GLUE diagnostic: grid.pack.ns; local to this run/rank.
    GlueGridPackNanos = 123,
    /// GLUE diagnostic: grid.gather.calls; local to this run/rank.
    GlueGridGatherCalls = 124,
    /// GLUE diagnostic: grid.gather.ns; local to this run/rank.
    GlueGridGatherNanos = 125,
    /// GLUE diagnostic: grid.reduce.calls; local to this run/rank.
    GlueGridReduceCalls = 126,
    /// GLUE diagnostic: grid.reduce.ns; local to this run/rank.
    GlueGridReduceNanos = 127,
    /// GLUE diagnostic: grid.audit_local.calls; local to this run/rank.
    GlueGridAuditLocalCalls = 128,
    /// GLUE diagnostic: grid.audit_local.ns; local to this run/rank.
    GlueGridAuditLocalNanos = 129,
    /// GLUE diagnostic: dist.mix.calls; local to this run/rank.
    GlueDistMixCalls = 130,
    /// GLUE diagnostic: dist.mix.ns; local to this run/rank.
    GlueDistMixNanos = 131,
    /// GLUE diagnostic: dist.prepare.calls; local to this run/rank.
    GlueDistPrepareCalls = 132,
    /// GLUE diagnostic: dist.prepare.ns; local to this run/rank.
    GlueDistPrepareNanos = 133,
    /// GLUE diagnostic: dist.partial.calls; local to this run/rank.
    GlueDistLocalCalls = 134,
    /// GLUE diagnostic: dist.partial.ns; local to this run/rank.
    GlueDistLocalNanos = 135,
    /// GLUE diagnostic: dist.gather.calls; local to this run/rank.
    GlueDistGatherCalls = 136,
    /// GLUE diagnostic: dist.gather.ns; local to this run/rank.
    GlueDistGatherNanos = 137,
    /// GLUE diagnostic: dist.reduce.calls; local to this run/rank.
    GlueDistReduceCalls = 138,
    /// GLUE diagnostic: dist.reduce.ns; local to this run/rank.
    GlueDistReduceNanos = 139,
    /// GLUE diagnostic: dist.solve.calls; local to this run/rank.
    GlueDistSolveCalls = 140,
    /// GLUE diagnostic: dist.solve.ns; local to this run/rank.
    GlueDistSolveNanos = 141,
    /// GLUE diagnostic: dist.update.calls; local to this run/rank.
    GlueDistUpdateCalls = 142,
    /// GLUE diagnostic: dist.update.ns; local to this run/rank.
    GlueDistUpdateNanos = 143,
    /// GLUE diagnostic: dist.charge.calls; local to this run/rank.
    GlueDistChargeCalls = 144,
    /// GLUE diagnostic: dist.charge.ns; local to this run/rank.
    GlueDistChargeNanos = 145,
    /// GLUE diagnostic: serial.mix.calls; local to this run/rank.
    GlueSerialMixCalls = 146,
    /// GLUE diagnostic: serial.mix.ns; local to this run/rank.
    GlueSerialMixNanos = 147,
    /// GLUE diagnostic: serial.fit.calls; local to this run/rank.
    GlueSerialFitCalls = 148,
    /// GLUE diagnostic: serial.fit.ns; local to this run/rank.
    GlueSerialFitNanos = 149,
    /// GLUE diagnostic: serial.charge.calls; local to this run/rank.
    GlueSerialChargeCalls = 150,
    /// GLUE diagnostic: serial.charge.ns; local to this run/rank.
    GlueSerialChargeNanos = 151,
    /// GLUE diagnostic: serial.norms.calls; local to this run/rank.
    GlueSerialNormsCalls = 152,
    /// GLUE diagnostic: serial.norms.ns; local to this run/rank.
    GlueSerialNormsNanos = 153,
    /// GLUE diagnostic: grid.local_points; local to this run/rank.
    GlueGridPoints = 154,
    /// GLUE diagnostic: grid.local_coefficients; local to this run/rank.
    GlueGridCoefficients = 155,
    /// GLUE diagnostic: grid.xc_points; local to this run/rank.
    GlueGridXcPoints = 156,
    /// GLUE diagnostic: grid.send_values; local to this run/rank.
    GlueGridSendValues = 157,
    /// GLUE diagnostic: grid.recv_values; local to this run/rank.
    GlueGridRecvValues = 158,
    /// GLUE diagnostic: dist.local_points; local to this run/rank.
    GlueDistPoints = 159,
    /// GLUE diagnostic: dist.history; local to this run/rank.
    GlueDistHistory = 160,
    /// GLUE diagnostic: dist.dot_elements; local to this run/rank.
    GlueDistDotElements = 161,
    /// GLUE diagnostic: dist.send_values; local to this run/rank.
    GlueDistSendValues = 162,
    /// GLUE diagnostic: dist.recv_values; local to this run/rank.
    GlueDistRecvValues = 163,
    /// GLUE diagnostic: dist.update_elements; local to this run/rank.
    GlueDistUpdateElements = 164,
    /// GLUE diagnostic: serial.points; local to this run/rank.
    GlueSerialPoints = 165,
    /// GLUE diagnostic: serial.history; local to this run/rank.
    GlueSerialHistory = 166,
    /// GLUE diagnostic: serial.update_elements; local to this run/rank.
    GlueSerialUpdateElements = 167,
    GlueDensityEligibleCalls = 168,
    GlueDensityFusedCalls = 169,
    GlueDensityPairedPoints = 170,
    GlueDensityLogicalBytes = 171,
    GlueDensitySavedBytes = 172,
    OwnedInputCalls = 173,
    OwnedInputDirectCalls = 174,
    OwnedInputAllocatedBytes = 175,
    OwnedInputBaselineBytes = 176,
    OwnedInputSavedBytes = 177,
}
pub const RUN_COUNTERS: usize = 178;
pub const RUN_COUNTER_LABELS: [&str; RUN_COUNTERS] = [
    "stop_kept",
    "stop_fft",
    "ewald_hits",
    "ewald_computed",
    "problems",
    "session_steps",
    "verifications",
    "owned.states_full_builds",
    "owned.states_table_bytes",
    "owned.states_replicated_bytes",
    "owned.states_recipe_bytes",
    "phase_validate.constructors",
    "phase_validate.calls",
    "phase_validate.pair_checks",
    "phase_validate.fast_accepts",
    "phase_validate.fallbacks",
    "phase_validate.nanoseconds",
    "init_lane.preflight.calls",
    "init_lane.preflight.ns",
    "init_lane.structure.calls",
    "init_lane.structure.ns",
    "init_lane.upf.calls",
    "init_lane.upf.ns",
    "init_lane.prepare.calls",
    "init_lane.prepare.ns",
    "init_lane.check_basis.calls",
    "init_lane.check_basis.ns",
    "init_lane.check_grid.calls",
    "init_lane.check_grid.ns",
    "init_lane.projectors.calls",
    "init_lane.projectors.ns",
    "init_lane.factors.calls",
    "init_lane.factors.ns",
    "init_lane.projector_assembly.calls",
    "init_lane.projector_assembly.ns",
    "init_lane.owned_rows_build.calls",
    "init_lane.owned_rows_build.ns",
    "init_lane.ionic.calls",
    "init_lane.ionic.ns",
    "init_lane.phase_tables.calls",
    "init_lane.phase_tables.ns",
    "init_lane.problem.calls",
    "init_lane.problem.ns",
    "init_lane.state_grid.calls",
    "init_lane.state_grid.ns",
    "init_lane.hamiltonian.calls",
    "init_lane.hamiltonian.ns",
    "init_lane.band_density.calls",
    "init_lane.band_density.ns",
    "init_lane.core_density.calls",
    "init_lane.core_density.ns",
    "init_lane.start_density.calls",
    "init_lane.start_density.ns",
    "init_lane.ionic_fft.calls",
    "init_lane.ionic_fft.ns",
    "init_lane.worker_ffts.calls",
    "init_lane.worker_ffts.ns",
    "init_lane.preflight.wave_rows",
    "init_lane.preflight.density_rows",
    "init_lane.preflight.grid_points",
    "init_lane.factors.rows",
    "init_lane.factors.shells",
    "init_lane.factors.radial_transforms",
    "init_lane.factors.radial_samples",
    "init_lane.factors.harmonics",
    "init_lane.factors.bytes",
    "init_lane.projectors.phase_pairs",
    "init_lane.projectors.materialized_bytes",
    "init_lane.projectors.owned_calls",
    "init_lane.ionic.rows",
    "init_lane.ionic.shells",
    "init_lane.ionic.radial_samples",
    "init_lane.ionic.phase_pairs",
    "init_lane.ionic.table_species",
    "init_lane.problem.atoms",
    "init_lane.problem.bands",
    "init_lane.problem.kpoints",
    "init_lane.problem.grid_points",
    "init_lane.problem.wave_rows",
    "init_lane.problem.density_rows",
    "init_lane.atomic.span.calls",
    "init_lane.atomic.span.ns",
    "init_lane.atomic.mgs.calls",
    "init_lane.atomic.mgs.ns",
    "init_lane.atomic.block_fill.calls",
    "init_lane.atomic.block_fill.ns",
    "init_lane.atomic.block_gram.calls",
    "init_lane.atomic.block_gram.ns",
    "init_lane.atomic.block_factor.calls",
    "init_lane.atomic.block_factor.ns",
    "init_lane.atomic.block_update.calls",
    "init_lane.atomic.block_update.ns",
    "init_lane.atomic.block_gather.calls",
    "init_lane.atomic.block_gather.ns",
    "init_lane.atomic.rows",
    "init_lane.atomic.columns",
    "init_lane.atomic.input_bytes",
    "init_lane.atomic.mgs_projection_macs_upper_bound",
    "init_lane.atomic.local_attempts",
    "init_lane.atomic.local_accepted",
    "init_lane.atomic.local_fallbacks",
    "init_lane.atomic.block_gram_flops",
    "init_lane.atomic.block_update_flops",
    "init_lane.atomic.grid_points",
    "init_lane.atomic.grid_x",
    "init_lane.atomic.grid_y",
    "init_lane.atomic.grid_z",
    "init_lane.atomic.bands",
    "glue_lane.grid.potential.calls",
    "glue_lane.grid.potential.ns",
    "glue_lane.grid.audit.calls",
    "glue_lane.grid.audit.ns",
    "glue_lane.grid.forward.calls",
    "glue_lane.grid.forward.ns",
    "glue_lane.grid.hartree.calls",
    "glue_lane.grid.hartree.ns",
    "glue_lane.grid.inverse.calls",
    "glue_lane.grid.inverse.ns",
    "glue_lane.grid.xc.calls",
    "glue_lane.grid.xc.ns",
    "glue_lane.grid.assemble.calls",
    "glue_lane.grid.assemble.ns",
    "glue_lane.grid.pack.calls",
    "glue_lane.grid.pack.ns",
    "glue_lane.grid.gather.calls",
    "glue_lane.grid.gather.ns",
    "glue_lane.grid.reduce.calls",
    "glue_lane.grid.reduce.ns",
    "glue_lane.grid.audit_local.calls",
    "glue_lane.grid.audit_local.ns",
    "glue_lane.dist.mix.calls",
    "glue_lane.dist.mix.ns",
    "glue_lane.dist.prepare.calls",
    "glue_lane.dist.prepare.ns",
    "glue_lane.dist.partial.calls",
    "glue_lane.dist.partial.ns",
    "glue_lane.dist.gather.calls",
    "glue_lane.dist.gather.ns",
    "glue_lane.dist.reduce.calls",
    "glue_lane.dist.reduce.ns",
    "glue_lane.dist.solve.calls",
    "glue_lane.dist.solve.ns",
    "glue_lane.dist.update.calls",
    "glue_lane.dist.update.ns",
    "glue_lane.dist.charge.calls",
    "glue_lane.dist.charge.ns",
    "glue_lane.serial.mix.calls",
    "glue_lane.serial.mix.ns",
    "glue_lane.serial.fit.calls",
    "glue_lane.serial.fit.ns",
    "glue_lane.serial.charge.calls",
    "glue_lane.serial.charge.ns",
    "glue_lane.serial.norms.calls",
    "glue_lane.serial.norms.ns",
    "glue_lane.grid.local_points",
    "glue_lane.grid.local_coefficients",
    "glue_lane.grid.xc_points",
    "glue_lane.grid.send_values",
    "glue_lane.grid.recv_values",
    "glue_lane.dist.local_points",
    "glue_lane.dist.history",
    "glue_lane.dist.dot_elements",
    "glue_lane.dist.send_values",
    "glue_lane.dist.recv_values",
    "glue_lane.dist.update_elements",
    "glue_lane.serial.points",
    "glue_lane.serial.history",
    "glue_lane.serial.update_elements",
    "glue_lane.density.eligible_calls",
    "glue_lane.density.fused_calls",
    "glue_lane.density.paired_band_points",
    "glue_lane.density.logical_bytes",
    "glue_lane.density.saved_bytes",
    "owned_lane.input.calls",
    "owned_lane.input.direct_calls",
    "owned_lane.input.allocated_bytes",
    "owned_lane.input.baseline_bytes",
    "owned_lane.input.saved_bytes",
];

/// Lane GLUE diagnostics, cached at plan/mixer construction. Disabled hooks
/// perform only an Option check: no clocks, counters, allocation or collectives.
/// Nested spans are inclusive; their values are never read by numerical code.
#[derive(Clone, Copy, Debug)]
pub struct GlueProfile {
    run: RunConfig,
    timed: bool,
}
impl GlueProfile {
    pub fn current() -> Option<Self> {
        let run = RunConfig::current();
        (run.on(Id::GLUE_PROFILE) && run.on(Id::PHASE_PROFILE)).then(|| Self {
            run,
            timed: run.on(Id::STAGE_TIMERS),
        })
    }
    pub fn count(self, counter: RunCounter, value: usize) {
        self.run.count(counter, value as u64);
    }
    /// Entered calls include failures; RAII closes their clocks on every exit.
    pub fn span(self, calls: RunCounter, nanos: RunCounter) -> InitSpan {
        self.run.count(calls, 1);
        InitSpan {
            run: self.run,
            nanos,
            start: self.timed.then(std::time::Instant::now),
        }
    }
}

/// Initialization-only diagnostics. The disabled path performs no clocks,
/// counter updates, allocations or size traversals. Callers keep this token
/// outside numerical loops; nested spans are inclusive, not additive.
#[derive(Clone, Copy)]
pub struct InitProfile {
    run: RunConfig,
    timed: bool,
}
impl InitProfile {
    pub fn current() -> Option<Self> {
        let run = RunConfig::current();
        run.on(Id::INIT_PROFILE).then(|| Self {
            run,
            timed: run.on(Id::STAGE_TIMERS),
        })
    }
    pub fn count(self, counter: RunCounter, value: usize) {
        self.run.count(counter, value as u64);
    }
    /// Count an entered call, including calls that later return an error.
    pub fn span(self, calls: RunCounter, nanos: RunCounter) -> InitSpan {
        self.run.count(calls, 1);
        InitSpan {
            run: self.run,
            nanos,
            start: self.timed.then(std::time::Instant::now),
        }
    }
}
/// RAII ensures early errors do not leave an open timer. No I/O or collectives.
pub struct InitSpan {
    run: RunConfig,
    nanos: RunCounter,
    start: Option<std::time::Instant>,
}
impl Drop for InitSpan {
    fn drop(&mut self) {
        if let Some(start) = self.start {
            self.run.count(
                self.nanos,
                start.elapsed().as_nanos().min(u64::MAX as u128) as u64,
            );
        }
    }
}

/// The path a consulted switch actually took (recorded by the code that acted
/// on it, or by the result writer from its counters).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PathNote {
    pub path: String,
    pub fallback: Option<String>,
}

pub struct RunState {
    run_id: AtomicU64,
    entries: Box<[Resolved]>,
    consulted: Box<[AtomicBool]>,
    /// The consulted set at the end of the computation (`freeze_consulted`),
    /// before result writers read switches for their report blocks.
    frozen: OnceLock<Box<[bool]>>,
    notes: Mutex<Vec<Option<PathNote>>>,
    counters: [AtomicU64; RUN_COUNTERS],
}
impl std::fmt::Debug for RunState {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("RunState")
            .field("run_id", &self.run_id.load(Ordering::Relaxed))
            .finish_non_exhaustive()
    }
}

/// Immutable switch values of one run (see the section comment).
#[derive(Clone, Copy, Debug)]
pub struct RunConfig(&'static RunState);

impl PartialEq for RunConfig {
    /// Same configuration object (not merely equal values: see
    /// [`RunConfig::agreement_hash`] for value equality).
    fn eq(&self, other: &Self) -> bool {
        std::ptr::eq(self.0, other.0)
    }
}

fn fresh_run_id() -> u64 {
    static SEQ: AtomicU64 = AtomicU64::new(0);
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos() as u64);
    let seq = SEQ.fetch_add(1, Ordering::Relaxed);
    // splitmix64 of time, pid and a sequence number: distinct per process run
    // and per configuration within a process.
    let mut z = nanos ^ (u64::from(std::process::id()) << 32) ^ seq.wrapping_mul(0x9e37);
    z = z.wrapping_add(0x9e37_79b9_7f4a_7c15);
    z = (z ^ (z >> 30)).wrapping_mul(0xbf58_476d_1ce4_e5b9);
    z = (z ^ (z >> 27)).wrapping_mul(0x94d0_49bb_1331_11eb);
    z ^ (z >> 31)
}

impl RunConfig {
    /// Resolve every switch from `env` (a variable lookup; the process
    /// environment for production runs) and the input file's switch values
    /// (`(table, key, text)`), with precedence environment > input > default.
    /// Environment values keep the legacy rules exactly (an invalid value falls
    /// back as before); input values are strict: an unknown key, a key in the
    /// wrong table or a value outside the allowed set is an error.
    pub fn resolve_with(
        env: impl Fn(&str) -> Option<OsString>,
        input: &[(Table, String, String)],
    ) -> Result<Self, String> {
        let mut texts: Vec<Option<String>> = vec![None; TABLE.len()];
        for (table, key, text) in input {
            let Some(s) = TABLE.iter().find(|s| input_key(s.id) == key.as_str()) else {
                return Err(format!(
                    "[{}] {key}: unknown switch (keys are the registry names without BAANDI_, lower case)",
                    table.as_str()
                ));
            };
            let home = table_of(s.id);
            if home != *table {
                return Err(match home {
                    Table::EnvOnly => format!(
                        "[{}] {key}: {} is not settable from the input file (environment only)",
                        table.as_str(),
                        s.name
                    ),
                    _ => format!(
                        "[{}] {key}: belongs to the [{}] table",
                        table.as_str(),
                        home.as_str()
                    ),
                });
            }
            if texts[s.id as usize].replace(text.clone()).is_some() {
                return Err(format!("[{}] {key}: given twice", table.as_str()));
            }
        }
        let mut entries = Vec::with_capacity(TABLE.len());
        for (s, input) in TABLE.iter().zip(texts) {
            let raw = env(s.name);
            let text = raw.as_ref().and_then(|r| r.to_str()).map(str::to_string);
            let (value, source) = match (&raw, &input) {
                (None, Some(t)) => (resolve_input(s, t)?, Source::Input),
                _ => resolve(s, raw.as_deref()),
            };
            entries.push(Resolved {
                raw,
                text,
                input,
                value,
                source,
            });
        }
        Ok(Self::leak(entries))
    }

    /// The process environment and the input file's values.
    pub fn resolve_env_and_input(input: &[(Table, String, String)]) -> Result<Self, String> {
        Self::resolve_with(|name| std::env::var_os(name), input)
    }

    /// Registry defaults and the given input values only (no environment):
    /// for tests and in-process comparisons that must not see the caller's
    /// environment.
    pub fn from_input(input: &[(Table, String, String)]) -> Result<Self, String> {
        Self::resolve_with(|_| None, input)
    }

    /// `(key, text)` pairs of the `[performance]` table plus `[profile]` ones,
    /// convenience for tests: `RunConfig::from_pairs(&[("stop_metric", "fft")])`
    /// (the table is taken from the registry).
    pub fn from_pairs(pairs: &[(&str, &str)]) -> Result<Self, String> {
        let input: Vec<_> = pairs
            .iter()
            .map(|(k, v)| {
                let table = TABLE
                    .iter()
                    .find(|s| input_key(s.id) == *k)
                    .map_or(Table::Performance, |s| table_of(s.id));
                (table, (*k).to_string(), (*v).to_string())
            })
            .collect();
        Self::from_input(&input)
    }

    fn leak(entries: Vec<Resolved>) -> Self {
        let consulted = (0..entries.len()).map(|_| AtomicBool::new(false)).collect();
        let notes = Mutex::new(vec![None; entries.len()]);
        Self(Box::leak(Box::new(RunState {
            run_id: AtomicU64::new(fresh_run_id()),
            entries: entries.into_boxed_slice(),
            consulted,
            frozen: OnceLock::new(),
            notes,
            counters: [const { AtomicU64::new(0) }; RUN_COUNTERS],
        })))
    }

    /// The configuration the calling thread sees: its scope, else the process
    /// configuration.
    pub fn current() -> Self {
        SCOPE.with(Cell::get).unwrap_or_else(process)
    }

    /// The thread's scope only (`None` outside [`RunConfig::enter`]).
    pub fn scoped() -> Option<Self> {
        SCOPE.with(Cell::get)
    }

    /// Make this the configuration of the calling thread until the guard is
    /// dropped (nested scopes restore the previous one). Legacy accessors on
    /// this thread, and on `WorkerPool` workers running its parallel regions,
    /// read this configuration.
    #[must_use]
    pub fn enter(self) -> ScopeGuard {
        ScopeGuard(SCOPE.with(|s| s.replace(Some(self))))
    }

    /// Run `f` with `scope` as the thread's configuration (`None`: unchanged).
    /// The worker pool wraps its tasks with it.
    pub fn with_scope<R>(scope: Option<Self>, f: impl FnOnce() -> R) -> R {
        match scope {
            Some(c) => {
                let _g = c.enter();
                f()
            }
            None => f(),
        }
    }

    fn get(self, id: Id) -> &'static Resolved {
        let i = id as usize;
        let c = &self.0.consulted[i];
        if !c.load(Ordering::Relaxed) {
            c.store(true, Ordering::Relaxed);
        }
        &self.0.entries[i]
    }

    /// Effective value of a choice switch.
    pub fn choice(self, id: Id) -> &'static str {
        match self.get(id).value {
            Value::Choice(v) => v,
            _ => panic!("{} is not a choice switch", switch(id).name),
        }
    }
    /// `true` when the effective value of an on/off switch is `on`.
    pub fn on(self, id: Id) -> bool {
        self.choice(id) == "on"
    }
    /// Effective value of a number switch.
    pub fn number(self, id: Id) -> Option<f64> {
        match self.get(id).value {
            Value::Number(v) => v,
            _ => panic!("{} is not a number switch", switch(id).name),
        }
    }
    /// The environment text, else the input text (`None` when neither is set).
    pub fn text(self, id: Id) -> Option<&'static str> {
        let r = self.get(id);
        if r.raw.is_some() {
            r.text.as_deref()
        } else {
            r.input.as_deref()
        }
    }
    /// The variable's raw environment value (`None` when unset).
    pub fn os(self, id: Id) -> Option<&'static OsStr> {
        self.get(id).raw.as_deref()
    }
    /// Whether the switch was set (environment, any value; input file: a
    /// presence switch set to `on`, any other switch given).
    pub fn present(self, id: Id) -> bool {
        let r = self.get(id);
        r.raw.is_some()
            || match (&r.input, &r.value) {
                (Some(_), Value::Present(p)) => *p,
                (Some(_), _) => true,
                (None, _) => false,
            }
    }
    /// Where the effective value came from.
    pub fn source(self, id: Id) -> Source {
        self.0.entries[id as usize].source
    }
    /// Whether the switch was read under this configuration.
    pub fn consulted(self, id: Id) -> bool {
        match self.0.frozen.get() {
            Some(f) => f[id as usize],
            None => self.0.consulted[id as usize].load(Ordering::Relaxed),
        }
    }
    /// Freeze the consulted set (the first call wins): the computation is
    /// over, later reads (result writers reporting a switch) do not count.
    pub fn freeze_consulted(self) {
        let _ = self.0.frozen.get_or_init(|| {
            self.0
                .consulted
                .iter()
                .map(|c| c.load(Ordering::Relaxed))
                .collect()
        });
    }

    /// Record the path a switch actually took and, when it did not act as
    /// requested, why (the last note wins). Outside hot loops only.
    pub fn note_path(self, id: Id, path: &str, fallback: Option<&str>) {
        if let Ok(mut n) = self.0.notes.lock() {
            n[id as usize] = Some(PathNote {
                path: path.to_string(),
                fallback: fallback.map(str::to_string),
            });
        }
    }
    pub fn path_note(self, id: Id) -> Option<PathNote> {
        self.0
            .notes
            .lock()
            .ok()
            .and_then(|n| n[id as usize].clone())
    }

    /// Add to a run-local counter.
    pub fn count(self, c: RunCounter, n: u64) {
        self.0.counters[c as usize].fetch_add(n, Ordering::Relaxed);
    }
    /// Run-local counters in [`RUN_COUNTER_LABELS`] order.
    pub fn counters(self) -> [u64; RUN_COUNTERS] {
        std::array::from_fn(|i| self.0.counters[i].load(Ordering::Relaxed))
    }

    /// Run ID (random per configuration; MPI runs adopt rank 0's value).
    pub fn run_id(self) -> u64 {
        self.0.run_id.load(Ordering::Relaxed)
    }
    pub fn run_id_hex(self) -> String {
        format!("{:016x}", self.run_id())
    }
    /// Adopt another run ID (MPI: rank 0's, after the agreement exchange).
    pub fn adopt_run_id(self, id: u64) {
        self.0.run_id.store(id, Ordering::Relaxed);
    }

    /// FNV-1a hash of the effective values of every switch in the agreement
    /// set ([`in_agreement`]), in table order. Equal on every rank of a run
    /// whose ranks will take the same branches and enter the same collectives.
    pub fn agreement_hash(self) -> u64 {
        let mut h: u64 = 0xcbf2_9ce4_8422_2325;
        let mut eat = |bytes: &[u8]| {
            for b in bytes {
                h ^= u64::from(*b);
                h = h.wrapping_mul(0x0100_0000_01b3);
            }
        };
        for (s, r) in TABLE.iter().zip(self.0.entries.iter()) {
            if !in_agreement(s.id) {
                continue;
            }
            eat(s.name.as_bytes());
            eat(b"=");
            match &r.value {
                Value::Choice(c) => eat(c.as_bytes()),
                Value::Number(Some(x)) => eat(&x.to_bits().to_le_bytes()),
                Value::Number(None) => eat(b"none"),
                Value::Raw(Some(t)) => eat(t.as_bytes()),
                Value::Raw(None) => eat(b"unset"),
                Value::Present(p) => eat(if *p { b"1" } else { b"0" }),
            }
            eat(b";");
        }
        h
    }

    /// The whole table with effective values, sources and the run's notes.
    pub fn manifest(self) -> Vec<Entry> {
        let notes = self
            .0
            .notes
            .lock()
            .map(|n| n.clone())
            .unwrap_or_else(|_| vec![None; TABLE.len()]);
        TABLE
            .iter()
            .zip(self.0.entries.iter())
            .zip(notes)
            .enumerate()
            .map(|(i, ((s, r), note))| Entry {
                switch: s,
                value: r.value.clone(),
                source: r.source,
                requested: if r.raw.is_some() {
                    Some(
                        r.text
                            .clone()
                            .unwrap_or_else(|| "<non-unicode>".to_string()),
                    )
                } else {
                    r.input.clone()
                },
                consulted: self.consulted(TABLE[i].id),
                note,
            })
            .collect()
    }
}

/// Restores the previous thread scope on drop.
pub struct ScopeGuard(Option<RunConfig>);
impl Drop for ScopeGuard {
    fn drop(&mut self) {
        let previous = self.0.take();
        SCOPE.with(|s| s.set(previous));
    }
}

thread_local! {
    static SCOPE: Cell<Option<RunConfig>> = const { Cell::new(None) };
}
static PROCESS: OnceLock<RunConfig> = OnceLock::new();

fn process() -> RunConfig {
    *PROCESS.get_or_init(|| {
        RunConfig::resolve_with(|name| std::env::var_os(name), &[])
            .expect("the environment alone always resolves")
    })
}

/// Make `config` the process configuration (the CLI does this once, before
/// the first switch is read, with the input file's tables). If the process
/// configuration exists already (a switch was read first, or another run of
/// this process installed one), the call succeeds only when the effective
/// values agree ([`RunConfig::agreement_hash`]); otherwise it reports the
/// conflict and the caller must use a thread scope.
pub fn install_process(config: RunConfig) -> Result<(), String> {
    let mut installed = false;
    let current = *PROCESS.get_or_init(|| {
        installed = true;
        config
    });
    if installed || current.agreement_hash() == config.agreement_hash() {
        Ok(())
    } else {
        Err("switches were resolved for this process before the input file's [performance]/[profile] tables; the run uses a thread-scoped configuration".to_string())
    }
}

/// The process configuration (resolving it from the environment if nothing
/// was installed).
pub fn process_config() -> RunConfig {
    process()
}

/// Table entry of a switch.
pub fn switch(id: Id) -> &'static Switch {
    &TABLE[id as usize]
}

/// Effective value of a choice switch (current configuration).
pub fn choice(id: Id) -> &'static str {
    RunConfig::current().choice(id)
}

/// `true` when the effective value of an on/off switch is `on`.
pub fn on(id: Id) -> bool {
    choice(id) == "on"
}

/// Effective value of a number switch (`BAANDI_ORTHO`: DGKS ratio, `None` =
/// two passes always).
pub fn number(id: Id) -> Option<f64> {
    RunConfig::current().number(id)
}

/// The variable's text as read at startup, else the input file's text (`None`
/// when neither is set or the variable is not Unicode). For override switches
/// and legacy reports of the raw value.
pub fn text(id: Id) -> Option<&'static str> {
    RunConfig::current().text(id)
}

/// The variable's raw value as read at startup (`None` when unset).
pub fn os(id: Id) -> Option<&'static OsStr> {
    RunConfig::current().os(id)
}

/// Whether the switch was set (environment or input file).
pub fn present(id: Id) -> bool {
    RunConfig::current().present(id)
}

/// The whole table with effective values and sources, in table order
/// (current configuration).
pub fn manifest() -> Vec<Entry> {
    RunConfig::current().manifest()
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

    fn env_of(pairs: &'static [(&'static str, &'static str)]) -> impl Fn(&str) -> Option<OsString> {
        move |name| {
            pairs
                .iter()
                .find(|(k, _)| *k == name)
                .map(|(_, v)| OsString::from(v))
        }
    }
    fn perf(key: &str, value: &str) -> (Table, String, String) {
        (Table::Performance, key.to_string(), value.to_string())
    }

    #[test]
    fn run_config_precedence_env_over_input_over_default() {
        let input = [perf("m2_grid", "off"), perf("eigen_factor", "0.05")];
        let c = RunConfig::resolve_with(env_of(&[("BAANDI_M2_GRID", "on")]), &input).unwrap();
        assert_eq!(c.choice(Id::M2_GRID), "on");
        assert_eq!(c.source(Id::M2_GRID), Source::Env);
        assert_eq!(c.number(Id::EIGEN_FACTOR), Some(0.05));
        assert_eq!(c.source(Id::EIGEN_FACTOR), Source::Input);
        assert_eq!(c.text(Id::EIGEN_FACTOR), Some("0.05"));
        assert_eq!(c.choice(Id::B6), "on");
        assert_eq!(c.source(Id::B6), Source::Default);
        // Legacy environment semantics are unchanged even with an input value.
        let c = RunConfig::resolve_with(env_of(&[("BAANDI_M2_GRID", "yes")]), &input).unwrap();
        assert_eq!(c.choice(Id::M2_GRID), "on");
        assert_eq!(c.source(Id::M2_GRID), Source::InvalidToDefault);
        let m = c.manifest();
        let e = &m[Id::M2_GRID as usize];
        assert_eq!(e.requested.as_deref(), Some("yes"));
        assert!(e.consulted);
        assert!(!m[Id::B5 as usize].consulted);
    }

    #[test]
    fn run_config_input_is_strict() {
        let err = |input: &[(Table, String, String)]| RunConfig::from_input(input).unwrap_err();
        assert!(err(&[perf("m2_grid", "maybe")]).contains("allowed values on, off"));
        assert!(err(&[perf("no_such", "on")]).contains("unknown switch"));
        assert!(err(&[perf("phase_profile", "on")]).contains("belongs to the [profile] table"));
        assert!(err(&[perf("pseudo_dir", "/x")]).contains("environment only"));
        assert!(err(&[perf("b6", "on"), perf("b6", "off")]).contains("given twice"));
        assert!(err(&[perf("buffer_tol", "1e-6")]).contains("allowed values"));
        let c = RunConfig::from_input(&[
            (Table::Profile, "phase_profile".into(), "on".into()),
            (Table::Profile, "dgemm_timers".into(), "on".into()),
            perf("ortho", "dgks:0.5"),
        ])
        .unwrap();
        assert!(c.on(Id::PHASE_PROFILE));
        assert!(c.present(Id::DGEMM_TIMERS));
        assert_eq!(c.number(Id::ORTHO), Some(0.5));
        let off = RunConfig::from_pairs(&[("dgemm_timers", "off")]).unwrap();
        assert!(!off.present(Id::DGEMM_TIMERS));
    }

    #[test]
    fn run_config_scopes_are_per_thread_and_follow_the_pool() {
        let a = RunConfig::from_pairs(&[("stop_metric", "fft")]).unwrap();
        let b = RunConfig::from_pairs(&[("stop_metric", "kept")]).unwrap();
        let pool = crate::pool::WorkerPool::new(3).unwrap();
        let run = |c: RunConfig, want: &'static str| {
            let pool = pool.clone();
            std::thread::spawn(move || {
                let _g = c.enter();
                for _ in 0..200 {
                    assert_eq!(choice(Id::STOP_METRIC), want);
                    let (l, r) = pool.install(|| {
                        crate::pool::WorkerPool::join(
                            || choice(Id::STOP_METRIC),
                            || choice(Id::STOP_METRIC),
                        )
                    });
                    assert_eq!((l, r), (want, want));
                    c.count(RunCounter::SessionSteps, 1);
                }
            })
        };
        let (ta, tb) = (run(a, "fft"), run(b, "kept"));
        ta.join().unwrap();
        tb.join().unwrap();
        assert_eq!(a.counters()[RunCounter::SessionSteps as usize], 200);
        assert_eq!(b.counters()[RunCounter::SessionSteps as usize], 200);
        assert!(RunConfig::scoped().is_none());
        assert_ne!(a.run_id(), b.run_id());
        {
            let _outer = a.enter();
            {
                let _inner = b.enter();
                assert_eq!(RunConfig::current(), b);
            }
            assert_eq!(RunConfig::current(), a);
        }
        assert!(RunConfig::scoped().is_none());
    }

    #[test]
    fn agreement_hash_tracks_effective_values_only() {
        let d1 = RunConfig::from_input(&[]).unwrap();
        let d2 = RunConfig::from_pairs(&[("m2_grid", "on")]).unwrap();
        assert_eq!(
            d1.agreement_hash(),
            d2.agreement_hash(),
            "on is the default"
        );
        let x = RunConfig::from_pairs(&[("m2_grid", "off")]).unwrap();
        assert_ne!(d1.agreement_hash(), x.agreement_hash());
        let s = RunConfig::from_pairs(&[("snapshot_dir", "/tmp/x")]).unwrap();
        assert_eq!(
            d1.agreement_hash(),
            s.agreement_hash(),
            "locations excluded"
        );
        let p = RunConfig::from_pairs(&[("phase_wait", "on")]).unwrap();
        assert_ne!(
            d1.agreement_hash(),
            p.agreement_hash(),
            "collective-gating profile switches are included"
        );
    }

    #[test]
    fn path_notes_are_reported_in_the_manifest() {
        let c = RunConfig::from_pairs(&[("m2_grid", "on")]).unwrap();
        c.note_path(Id::M2_GRID, "replicated", Some("one rank"));
        let e = &c.manifest()[Id::M2_GRID as usize];
        assert_eq!(
            e.note,
            Some(PathNote {
                path: "replicated".into(),
                fallback: Some("one rank".into())
            })
        );
    }
}
