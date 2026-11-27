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

fn ladder_ratio_valid(v: f64) -> bool {
    v.is_finite() && (1.0..=1e6).contains(&v)
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



    NL_ENERGY: Kind::Choice { allowed: &["fast", "exact"], default: "fast" }, "fast", "fast, exact",
        "nonlocal energy from blocked projections (exact: band-by-band operator)";
    ENERGY_IDENTITY: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "nonlocal energy of SCF iterations from the band-sum identity";
    // D1 ladder hand-off scan (work order 2026-11-14); defaults are the adopted ladder.
    STANDARD_RESIDUAL_GUARD: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "standard-class stop also requires sum_k w_k sum_occupied f_i |r_i|^2 x 1 Ha^-1 (empirical second-order band-energy estimator, not a bound; review F03) below the threshold; read only by the reference spin/FD driver; default off since 2026-11-27 (dev08 gen59: binds only on Au of six metals, one map fewer, dF 3.9e-9 Ha, unchanged against QE), aligned with the scalar k path and QE";




    // MPI entry cost (docs/performance/mpi-entry-cost-analysis.md). Adopted 2026-11-15 on dev08 A/B: M2_DGEMM, M2_ORTHO=serial and M2_LEAN default on; M2_FFT stays rust (neutral).
    M2_DGEMM: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "distributed Davidson rotations and projections and the M2 nonlocal products through the per-panel dgemm at any width (MKL builds: rounding class; pure builds bit-identical)";

    M2_FFT: Kind::Choice { allowed: &["rust", "mkl"], default: "rust" }, "rust", "rust, mkl",
        "1D kernels of the M2 distributed FFT (mkl: MKL DFTI batch lanes on MKL builds, rounding class; pure builds stay rustfft)";


    STAGE_TIMERS: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "F21: MPI stage and collective counters measure wall time (on, the earlier behaviour) or run count-only (off: calls and collective payload bytes counted, no clock read, stage seconds 0)";
    // MPI reduced eigenproblem (docs/performance/mpi-reduced-eigen-analysis.md); default off.

    // MPI solver workspace (review F06; ledger row M2-WORKSPACE-F06); bit-identical.

    // MPI selected-root rotation (review F05; ledger row M2-SELECTED-ROOTS-F05); default off until dev08.

    // MPI in-panel Gram-Schmidt reductions (docs/performance/profile-2026-11-16-mpi.md); default classic until dev08.

    // MPI kept start (ticket MPI-KEPT-START, docs/design/mpi-kept-start.md); default off until dev08.
    // MPI scaling gap (ticket MPI-SCALING-GAP, docs/performance/mpi-scaling-gap.md); default off.
    M2_STICKS: Kind::Choice { allowed: &["greedy", "mirror"], default: "mirror" }, "mirror",
        "greedy, mirror",
        "stick partition of the M2 distributed FFT: greedy (sticks one by one to the least loaded rank) or mirror (a stick and its Gamma mirror stick go to the same rank as one unit, so no ORC representative is split and the two mirror alltoalls of every apply vanish; rounding class: rows are summed in a different order)";

    // MPI overlap-matrix Davidson (ticket MPI-OVERLAP-DAVIDSON, docs/design/mpi-overlap-davidson.md); default off.

    // MPI lazy orbitals (ticket MPI-LAZY-ORBITALS, docs/design/mpi-lazy-orbitals.md); default off.
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
    M2_GRID: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "M2-GRID: on the M2 path (2+ ranks, LDA PZ/PW, no core density, density floor or symmetry) the full-grid stages outside the solve run on a slab-decomposed full-grid FFT (owned real-space planes, owned sticks, two transposes per transform): the input density's Hartree FFT pair, the LDA XC and the effective potential on the planes with one allgatherv of the effective potential and the partial sums; the energy audit's output-density Hartree FFT and grid sums on the planes with one allgatherv of partial sums; the stop-rule Coulomb estimate from the kept owned coefficients; every sum in a canonical plane/stick order, so results do not depend on the rank count; rounding class against the replicated path";

    M2_LOCAL_FIELDS: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "M2-OWNED (default on since 2026-11-18, dev08 si32/si64 np=8 wall -2.5 to -6.0 %; declines automatically for active-boundary consumers): with M2_GRID on the M2 standard path, the effective potential stays on the producing ranks' slab planes (the planes the distributed H apply multiplies; redistributed by one alltoall only when the slab and operator plane axes differ) and is fed to the distributed operator and its kept-start potential change (delta V) without the full-grid allgatherv; the complex Hamiltonian's full potential is gathered only for a consumer that needs it (cold atomic start, complex or M1 fallback), tagged with an epoch so a stale potential is rejected; bit-identical to M2_GRID alone";
    M2_LOCAL_DENSITY: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "M2-OWNED (default on since 2026-11-18, dev08 si32/si64 np=8 wall -2.5 to -6.0 %; declines automatically for active-boundary consumers): with M2_LOCAL_FIELDS and M2_MIX on, the output density stays on the owned slab planes through the energy audit, the stop rule, the session metrics (residual, charge, minimum: canonical per-plane sums), the distributed Pulay mixer (window on the owned planes) and the next potential; the full density is gathered only for a consumer that needs it (delivered result, checkpoint, ladder hand-off, declined map); rounding class against M2_MIX alone (per-plane instead of per-chunk canonical sums)";
    // Solvers (dft-solvers).









    DAVIDSON_STAGES: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "Davidson stage timers";



    // Per-core kernels (ticket PERCORE-KERNELS, docs/performance/percore-kernels-2026-11-17.md); default off.

    DAVIDSON_PROFILE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "Bounded per-call Davidson lane trace with shapes and rank/map labels; diagnostic only.";

    // History-preserving restart (ticket BREAKTHROUGH-1-RESTART, docs/design/restart-keep.md); default off.

    // Physics (dft-physics).

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
    // CLI (dft-cli).
    BINARY_HASH_CACHE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "LANE-LAUNCHIO: cache executable SHA-256 in a metadata-keyed sibling sidecar on Unix; preserves provenance, falls back to fresh hashing; trusted build directories only";
    RELAX_EXTRAPOLATION: Kind::Override { valid: relax_extrapolation_valid }, "config",
        "none, first_order, second_order", "overrides the relaxation density extrapolation of the input file";
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
    GLUE_PROFILE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "LANE-GLUE Phase 0: rank-local grid and mixer spans/work; requires PHASE_PROFILE; clocks follow STAGE_TIMERS; no numerical change";
    PHASE_PROFILE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "exclusive nested phase profile per rank (init / coarse SCF / fine ordinary SCF / verification pass and fail / materialization; map phases and solve stages inside; maximum over ranks and imbalance); clocks follow BAANDI_STAGE_TIMERS (off: counts only); printed by dft-cli (results.json phase_profile, parallel.phases), m2_entry_cost and davidson_bench";


    EWALD_FAST: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "INIT-EWALD-FAST: balanced Ewald scheme (dft-physics ewald::Scheme::Balanced): splitting parameter from the real/reciprocal cost balance for the cell volume and ion count, cutoffs from rigorous tail bounds (truncation <= min(tolerance (sum |Z|)^2 / Omega^(1/3), tolerance * N) Ha), ion pairs i <= j and the reciprocal half space with per-ion phase tables; energy, forces and stress share the parameters; rounding class (different truncation and summation, error vs a tight reference far below 1e-12 Ha/atom); off = reference scheme (eta sqrt(pi)/Omega^(1/3), bit-identical to the previous evaluator)";
    IONIC_FAST: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "INIT-IONIC-FAST: the structure factor of the local ionic potential sum_I exp(-i G.tau_I) from per-ion, per-axis phase tables exp(-i m_k b_k.tau_I) (two complex products per ion and G) instead of one sin/cos pair per ion and G; falls back to the direct phase for any G whose Cartesian vector is not the Miller combination of the unit-Miller vectors; rounding class (phase rounding, about 1e-13 relative)";
    INIT_PROFILE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "LANE-INIT phase 0: rank-local startup call/work counters and nested timers; STAGE_TIMERS=off retains counts without clocks; no numerical changes, default off";
    PHASE_VALIDATE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "G-PHASE-VALIDATE: owned projector finite-dot validation uses a sufficient numerical range proof (finite momentum and position components with maxima <= 2^500); otherwise retains the original pair loop and errors; coefficients and lazy materialization unchanged; off until dev08 A/B";

    M2_LAZY_RESULT: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "LAZY-RESULT: when the caller declared that it reads no delivered orbitals (dft-cli run without forces, stress, ELF or partial-charge output, classic mode), an M2 result keeps its lazy distributed orbital blocks instead of gathering the full N x B block on every rank at delivery; ScfResult::materialize_orbitals gathers them on demand (collective); energy, density, eigenvalues and history unchanged; off = always gather";
    PHASE_WAIT: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "wait probe of the phase profile (ticket BREAKTHROUGH-4-OTHER-AND-INIT): a timed barrier before every collective of the distributed solver and at the entry of every M2 solve, so that the time ranks spend waiting for the last rank (load imbalance, arrival skew) is reported apart from the collective's own transfer (`solve.wait`, `solve.entry_wait`, `coll.barrier`); changes the timing (extra barriers), never the arithmetic; diagnostic only, never an acceptance arm";
    SNAPSHOT_DIR: Kind::Path, "unset", "absolute directory",
        "SCF snapshot capture: every real Gamma solve of the run writes the state that replays it (input density, this rank's orbital rows, kept H X and its armed flag, tolerance, policy, layout hash, context) to this directory outside the timed map phases; replay with davidson_bench --snapshot <manifest.json>";
    FFT_SPARSE_CLEAR: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P3-FFT: clear only holes before coefficient/middle/slab scatter using cached contiguous ranges; same FFT order, default off pending dev08";
    FFT_DIRECT_PACK: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P3-FFT: pack distributed real H input ORC rows directly into FFT sticks, avoiding the intermediate coefficient write/read; default off pending dev08";
    FFT_FUSED_ROWS: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P3-FFT2: stream distributed local-potential last-axis inverse/multiply/forward through reusable rows with immutable zero-hole input; rustfft lanes only, rounding class, default off pending dev08";
    FFT_A2A_BALANCE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P3-A2A: cap rank stick counts to reduce padded FFT alltoall payload, preserving mirror units; use only if the block shrinks, row-reduction rounding class, default off pending dev08";
    // Adopted 2026-11-25 (dev08 one-run P3l: si128 35.67 -> 34.60 s; setup -0.35 s, output -0.66 s).
    NONLOCAL_SINGLE_GEMM: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "P3-NL: one MKL projection over all local rows instead of 256-row reduction panels; requires M2_DGEMM, rounding class, default off pending dev08";
    // Adopted 2026-11-25 (dev08 one-run P3m: si128 border -0.97 s, back_and_add -0.40 s; counts identical).
    NONLOCAL_FUSE_COUPLING: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "P3-ADDV: store DQ in the dead local projection buffer after rank reduction, omitting separate W storage; bit-identical, default off pending dev08";
    // Adopted 2026-11-25 (dev08 one-run P3m: si128 border -0.97 s, back_and_add -0.40 s; counts identical).
    NONLOCAL_BACK_ACCUM: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "P3-ADDV: one sequential MKL P W directly into HX with beta=1; requires M2_DGEMM and MKL_DGEMM, rounding class, default off pending dev08";

    REFSCF_RMM: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P3-RMM: bounded RMM-DIIS sweeps between Davidson SCF iterations 1, 2 and multiples of five; verification always Davidson";
    REFSCF_RMM_STEPS: Kind::Number { default: Some(4.), valid: |v| (1. ..=16.).contains(&v) && v.fract() == 0. }, "4", "integer 1..16",
        "P3-RMM: maximum DIIS steps per sweep; default four, experimental overrides require multi-system evidence";
    REFSCF_RMM_BLOCK: Kind::Number { default: Some(16.), valid: |v| (1. ..=64.).contains(&v) && v.fract() == 0. }, "16", "integer 1..64",
        "P3-RMM: band-order Gram-Schmidt panel width; default sixteen";

    FFT_KERNEL_BY_LENGTH: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P3-PERMAP: calibrate rustfft/MKL per actual lane length, batch and direction once per run; pure builds retain rustfft";
    REFSCF_PERMAP: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P3-PERMAP: cache unchanged Broyden Gram entries and fuse descf with map diagnostic passes/reduction";
    REFSCF_LADDER_RATIO: Kind::Number { default: Some(100.), valid: ladder_ratio_valid }, "100", "number in [1, 1e6]",
        "reference cutoff ladder: coarse stopping threshold = fine dr2 threshold times this ratio (material independent)";
    REFSCF_LADDER: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P3-LADDER: reference coarse cutoff 0.5, exact owned prolongation, fresh fine Broyden history";

    REFSCF_KEPT_RECOMPUTE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P4-DENSKEEP: keep Ritz HX but recompute f64 inverse planes for the local-potential delta; no retained slabs; retained HX enabled";
    REFSCF_MEMORY_TRIM: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P4-MEM: release dead start/map/delivery arrays and return band-sized Davidson storage; identical arithmetic; pending dev08";
    K_LEAN_ENERGY: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK4: k-point reference maps defer per-k kinetic/nonlocal band energies and the local/zero-mode integrals of ordinary (lean) maps to the final audit and reduce all per-k residual norms with the energy integrals in one tree sum per map; delivered energies, density, eigenvalues and counts bit-identical; intermediate audit fields kinetic/local/local_zero_mode/nonlocal become null; off = per-k terms and one residual reduction per k every map";

    FORCE_BAND_FIRST: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "G5: one projector/derivative GEMM over all bands per k for nonlocal forces; unmeasured, default off";
    REFSCF_SPIN_FD_MOMENT: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "R6a: reference collinear FD moment policies enabled by default; G4B 48/48 and dev08 gen15 fixed-moment Fe VERIFIED in 75.8 s; free-moment PBE remains G4D";

    XC_SCREEN_NEGATIVE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "G4E: screen XC-only small/signed vacuum densities with a screened-charge budget; accepted ADR 0046";
    REFSCF_REAL_START: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "G4G research: Hermitian Gamma kinetic complements and real starting Ritz matrix; rejected (centred fixture stops converging), default off";
    REFSCF_CANONICAL_CLUSTER: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "G4G research: canonical coordinates within roundoff-degenerate Ritz clusters; default off";
    CDAV_REUSE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK3: complex reference Davidson keeps its workspace across calls, merges error-flag agreements into its numerical collectives and skips finite scans of known-finite GEMM operands; bit-identical; off restores per-call allocation and per-product agreements";
    REFSCF_SPIN_PULAY: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "G4G research: per-channel real-space Pulay replaces the Broyden proposal (serial fixed-sector only); default off pending dev08 O run";
    K_START_MEMO: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK1: k atomic start computes each orbital radial transform once per (species channel, exact |k+G|) and reuses it over atoms and m; bit-identical; off restores the per-(row, atom, m) quadrature";
    K_NL_TABLE_MEMO: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK1B: k start builds projector radial/angular factors for the rank's owned rows only (rows-only table, outcome agreed by one collective); bit-identical; off restores the full k-sphere factor table on every rank";
    K_START_SPIN_SHARE: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK1B: collinear spin k start builds bases, grids, projector tables, atomic columns and the start Ritz once and gives every channel its own operator and a copy of the orbitals; bit-identical; off builds the whole start per channel";
    CDAV_BORDER_GRAM: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK10: complex reference Davidson keeps the projected overlap and Hamiltonian of the previous basis and computes only the border of the appended columns (one Gram product and one reduction each); rounding class; off restores the full recompute after every append";
    RADIAL_TABLE_START: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK14: k atomic start orbitals read their radial transforms from one uniform table per run and species channel (step 0.005 Bohr^-1, 8-point Lagrange) instead of a quadrature per |k+G|; start only (not bit-identical to the quadrature, delivered results inside the strict class); off restores the GK1 quadrature memo";
    RADIAL_TABLE_NL: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK14: k-point projector factors read their radial transforms from one uniform table per run and species projector (step 0.005 Bohr^-1, 8-point Lagrange, measured relative error <= 7.4e-15) instead of a quadrature per |k+G|; changes the k Hamiltonian at roundoff level (delivered energies inside the strict 1e-10 Ha class); Gamma-only path unchanged; off restores the quadrature";
    CDAV_BORDER_FACTOR: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK10B: complex reference Davidson extends the overlap Cholesky factor over the appended columns only (bordered Cholesky, leading factor kept) when the overlap was border-updated; rounding class; off refactors the whole overlap after every append";
    CDAV_STREAM: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK10B: complex reference Davidson computes its preconditioned corrections in one column pass with the preconditioner branch hoisted per column (as the real path's stream kernel); bit-identical; off restores the per-element loop";
    K_FIRST_EMPTY_RULE: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "GK10B research: the first k map of the fixed-occupation reference k driver applies the empty-band Ritz threshold to empty bands (occupations are known before the first solve) instead of converging every band to the occupied threshold; changes the map-1 density when no redo follows (SCF trajectory; one synthetic fixture left the standard class), default off pending a dev08 A/B on redo cases";
    K_START_KEPT_HX: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK10B: the unpolarized reference k start keeps H X of its Ritz orbitals (H Q C, one GEMM) and the k driver seeds its kept HX with it, so map 1 does not apply H to the start orbitals again (potential change corrected exactly); rounding class; off applies H in map 1";
    K_REDO_KEPT_HX: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK10B: the first-iteration redo of the reference k driver starts from the coarse solve's kept X and HX (same input potential, zero change) instead of applying H to the start again; rounding class; off invalidates the kept HX before the redo";
    REFSCF_K_POOLS: Kind::Choice { allowed: &["off", "auto", "on"], default: "auto" }, "auto", "off, auto, on",
        "GK9 k-point pools (ADR 0041 option D) for the reference k and spin/FD drivers: the ranks form npool pools (contiguous rank blocks, MPI_Comm_split), each pool owns a contiguous block of whole k points and distributes each k's plane-wave rows over its ranks exactly as the distributed layout does over all ranks; per map one potential allgatherv, one record allgatherv from the pool leaders and one density alltoall summed in global k order; on = one rank per pool (stage 1; with fewer k than ranks the largest divisor of the rank count not above the k count); auto = the pool count of the GK9B cost rule (divisors of the rank count not above the k count, minimum of max over pools of sum_k (N_pw,k / p + 32 ceil(log2 p)) rows, ties to more pools); off = distributed layout; npool 1 is the distributed layout bit for bit; rounding class across rank/pool counts, fixed (ranks, npool) bit-reproducible; REFSCF_K_NPOOL forces a count";
    REFSCF_K_NPOOL: Kind::Number { default: None, valid: |v| (1. ..=65536.).contains(&v) && v.fract() == 0. }, "unset", "integer 1..65536",
        "GK9B: force this many k-point pools (clamped to the rank count and the k count; 1 = distributed layout), overriding REFSCF_K_POOLS=auto/on (not off); pool sizes and k blocks follow dft_mpi::block_counts; for A/B tests";
    RADIAL_TABLE_FIELDS: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK16: the k-path common start (and the fixed-density k start) reads the atomic-charge, Coulomb-subtracted local and core-density |G| transforms of the ionic/atomic fields from one uniform radial table per species and function (step 0.005 Bohr^-1, 8-point Lagrange, nodes split over ranks) instead of a quadrature per distinct |G|; G = 0 stays exact; changes the local potential at roundoff level (rounding class, delivered energies inside the strict 1e-10 Ha class); Gamma-only path unchanged; off restores the per-shell quadrature";
    FFT_ONE_RANK_DIRECT: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GK17: the one-rank stick FFT of a complex k operator (every pooled k operator, every serial k run; the real Gamma operator keeps the transpose path) runs its stick transforms in per-vector batch buffers and reads/writes the middle stage there directly, without packing into the transpose send buffer, the self-alltoall copy and the unpack into the stick buffer; rustfft lanes only (vendor lanes keep the transpose path); bit-identical; off restores the transpose path";
    K_KEPT_SLABS: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "GK17 (plan GK15; GF1: also the spin/FD map, budget shared by its channels): the fixed-occupation reference k map keeps the real-space planes psi_b(r) that its density pass computes for the leading bands of each k (bounded by K_KEPT_SLAB_MB per rank), and the next map's kept-HX update reuses them, so (V_new - V_old) psi_b needs only the forward transform; bit-identical; costs the slab memory; off = inverse and forward transform per band";
    K_KEPT_SLAB_MB: Kind::Number { default: Some(1024.), valid: |v| (0. ..=1048576.).contains(&v) && v.fract() == 0. }, "1024", "integer MiB 0..1048576",
        "GK17: per-rank budget of K_KEPT_SLABS in MiB, split evenly over the rank's k points in k order (bands per k = min(bands, budget / (k count x plane bytes))); 0 keeps no slab";
    SPIN_FD_START_KEPT_HX: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GF1 (port of GK10B K_START_KEPT_HX to the spin/FD driver): the reference spin/FD k start keeps H X of its Ritz orbitals (H Q C, one GEMM; the FD buffer columns appended after the start get one H application each) and the driver seeds every channel's kept HX with it, so map 1 does not apply H to the start orbitals again (potential change corrected exactly); rounding class; off applies H in map 1";
    SPIN_FD_REDO_KEPT_HX: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GF1 (port of GK10B K_REDO_KEPT_HX to the spin/FD driver): the first-iteration redo starts from the coarse solve's kept X and HX (same input potential, zero change, no transforms) instead of applying H again; rounding class; off invalidates the kept HX before the redo";
    STRESS_FORCE_PASS: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "GD1: on = the frozen-Miller stress re-runs the complete force pass as its precondition check (former behaviour); off = it runs only the shared precondition checks (identity, certificate/verification, tolerances, ledger) without computing forces; stress bit-identical; precondition_force_band_passes reports 0 when off";
    DERIV_NL_GEMM: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GD1: nonlocal forces and nonlocal stress per k from one projector/gradient GEMM per atom (A = [beta | d beta/dq], B = [psi | q psi]) instead of the NonlocalOperator build plus per-G scalar loops; rounding class (reduction order); off restores the former per-k operator and scalar loops";
    DERIV_RADIAL_TABLES: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GD1: with two or more k blocks the derivative pass reads the projector transforms and their q derivatives from uniform radial tables (step 0.005 Bohr^-1, 8-point Lagrange, nodes split over ranks, node values equal to the quadrature) instead of one quadrature per distinct |k+G|; rounding class; off restores the quadrature";
    DERIV_SHELLS: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GD1: the local and core |G|-shell transforms of the derivative pass (local/core forces and stress, the rebuilt core density) are evaluated once per pass with the G-independent integrand (including erf) hoisted, split over the ranks in contiguous blocks and shared by one allgatherv, instead of per shell on every rank (and per atom for the core density); bit-identical; off restores the per-shell evaluation";
    FD_EMPTY_THRESHOLD: Kind::Choice { allowed: &["off", "weighted", "percent"], default: "weighted" }, "weighted", "off, weighted, percent",
        "GF2 (plan GK13): ordinary spin/FD SCF maps after the first solve converge a retained band whose previous-map occupation fraction (occupation / channel capacity) is at or below a cut with the empty-band Ritz threshold max(5 thr, 5e-6 Ha) instead of the occupied threshold thr; weighted: cut = sqrt(thr / empty threshold), so the band's density error, about f sqrt(threshold / gap), stays at or below that of a fully occupied band converged at thr; percent: cut = 0.01; the first map, verification maps and fixed-sector policies are unchanged; SCF trajectory change (accuracy class); off = every retained band at thr";
    FD_EMPTY_VERIFY: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "GF2: with FD_EMPTY_THRESHOLD on, the fixed-point verification map of the spin/FD driver also uses the empty-band rule with the occupations of the map it verifies (the verification criterion, the Coulomb-metric error of the fresh output density, is unchanged); off = the verification map converges every retained band at the occupied threshold";
    FD_BUFFER: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GF2: the spin/FD reference Davidson solves max(4, ceil(n / 10)) buffer columns beyond the n retained bands (empty-band threshold, never required, not in the density); off solves exactly n columns (experiment only: on the GF2 Mac fixtures Au, Al, Sr, Ru the top retained bands then converge to wrong eigenpairs, missed states up to 0.36 Ha, and Sr converges to a free energy 2.7e-10 Ha away, outside the strict class); SCF trajectory change";
    COLLINEAR_TIME_REVERSAL: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GT1 (ADR 0050): collinear spin runs accept basis.ktime_reversal (k ~ -k within each spin channel: each channel Hamiltonian has real local and nonlocal potentials without spin-orbit coupling, so H_s(-k) = conj H_s(k)); a k grid therefore reduces by time reversal by default as in unpolarized runs; off restores the former rejection (collinear input then needs basis.ktime_reversal = false)";
    FD_GROWTH_AT_STOP: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GB2: automatic FD band growth reacts to an empty-state coverage miss only on a map that meets the stop rule, the last map of the budget or a verification map (grow, solve the same input density again, continue); earlier maps that miss are recorded as transient (reference_scf.band_count_growth.transient_coverage_misses); off restores growth on every map that misses (pre-GB2)";
    FD_VERIFY_PREDICTOR: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "GB2: the spin/FD driver uses the closed-shell verification predictor: after a verification attempt with Coulomb error e at stop residual r0, a later stop at residual r < r0 skips its attempt (an ordinary map instead) when e r / r0 >= the threshold, never twice in a row; the pass criterion is unchanged and a run with at most one failed attempt is unchanged; off verifies at every stop";
    SPIN_FD_ATOMIC_DENSITY: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "GB2: collinear Fermi-Dirac runs start each channel from the atomic start density scaled to the channel count (N +- m)/2 when that density is strictly positive on every grid point (otherwise, and for fixed sectors, uniform channel counts as before); off restores the uniform channel start; SCF trajectory change (accuracy class)";
    K_REDUCTION: Kind::Choice { allowed: &["pointwise", "mesh", "unreduced"], default: "pointwise" }, "pointwise", "pointwise, mesh, unreduced",
        "GT2 (ADR 0051): symmetry.mode=auto k reduction. pointwise = orbits under every crystal operation, k and S k merged whenever S k lies on the mesh (QE-like); density, forces, stress and relaxation projectors use the whole crystal group; on a mesh that every operation preserves this is the mesh value bit for bit; on other meshes the k quadrature becomes the crystal-group average of the mesh (exactly the full mesh with crystal projection). mesh = orbits under the mesh-preserving subgroup only (former behaviour; exactly the unsymmetrized full mesh). unreduced = full mesh with the crystal-group projection (oracle of pointwise; no k saving)";
    GRID_INCOMPATIBLE_OPS: Kind::Choice { allowed: &["drop", "reciprocal"], default: "reciprocal" }, "reciprocal", "drop, reciprocal",
        "GS1D-B (ADR 0049): with an explicit basis.fft_grid, crystal operations whose rotation does not map the grid onto itself. drop = keep only the grid-permuting subgroup (GS1C/GS1D); reciprocal = keep them (QE-like): the density is projected onto the whole group in G space as before and the local potential (ionic + Hartree + XC, each spin channel) is projected onto the group in G space after every evaluation, so the Hamiltonian stays invariant although the real-space XC quadrature is not; k reduction, forces, stress and relaxation then use the whole group. Operations whose fractional translation is not a multiple of the grid spacing are dropped in both modes";
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
        | DGEMM_TIMERS => Table::Profile,
        // Overrides of input-file values (`scf.cutoff_ladder`, the relaxation
        // extrapolation) and the pseudopotential cache location.
        RELAX_EXTRAPOLATION | PSEUDO_DIR => Table::EnvOnly,
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
    /// convenience for tests: `RunConfig::from_pairs(&[("m2_screen", "fft")])`
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

/// Effective value of a number switch.
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
        assert_eq!(eval(Id::M2_DGEMM, Some("off")), (Value::Choice("off"), Env));
        assert_eq!(
            eval(Id::M2_DGEMM, Some("yes")),
            (Value::Choice("on"), InvalidToDefault)
        );
        // Default-off flags: only "on" is on.
        assert_eq!(
            eval(Id::PACK_FIXED, Some("1")),
            (Value::Choice("off"), InvalidToDefault)
        );

        assert_eq!(eval(Id::NL_ENERGY, Some("exact")).0, Value::Choice("exact"));

        assert_eq!(eval(Id::PSEUDO_DIR, Some("relative")).1, InvalidToError);
        assert_eq!(eval(Id::DGEMM_TIMERS, Some("")).0, Value::Present(true));
    }

    #[test]
    fn mpi_entry_cost_switches_defaults() {
        use Source::*;
        {
            let id = Id::M2_DGEMM;
            assert_eq!(eval(id, None), (Value::Choice("on"), Default));
            assert_eq!(eval(id, Some("off")), (Value::Choice("off"), Env));
            assert_eq!(eval(id, Some("1")), (Value::Choice("on"), InvalidToDefault));
        }

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
        let input = [perf("m2_grid", "off"), perf("refscf_ladder_ratio", "4")];
        let c = RunConfig::resolve_with(env_of(&[("BAANDI_M2_GRID", "on")]), &input).unwrap();
        assert_eq!(c.choice(Id::M2_GRID), "on");
        assert_eq!(c.source(Id::M2_GRID), Source::Env);
        assert_eq!(c.number(Id::REFSCF_LADDER_RATIO), Some(4.));
        assert_eq!(c.source(Id::REFSCF_LADDER_RATIO), Source::Input);
        assert_eq!(c.text(Id::REFSCF_LADDER_RATIO), Some("4"));
        assert_eq!(c.choice(Id::M2_DGEMM), "on");
        assert_eq!(c.source(Id::M2_DGEMM), Source::Default);
        // Legacy environment semantics are unchanged even with an input value.
        let c = RunConfig::resolve_with(env_of(&[("BAANDI_M2_GRID", "yes")]), &input).unwrap();
        assert_eq!(c.choice(Id::M2_GRID), "on");
        assert_eq!(c.source(Id::M2_GRID), Source::InvalidToDefault);
        let m = c.manifest();
        let e = &m[Id::M2_GRID as usize];
        assert_eq!(e.requested.as_deref(), Some("yes"));
        assert!(e.consulted);
        assert!(!m[Id::PACK_FIXED as usize].consulted);
    }

    #[test]
    fn run_config_input_is_strict() {
        let err = |input: &[(Table, String, String)]| RunConfig::from_input(input).unwrap_err();
        assert!(err(&[perf("m2_grid", "maybe")]).contains("allowed values on, off"));
        assert!(err(&[perf("no_such", "on")]).contains("unknown switch"));
        assert!(err(&[perf("phase_profile", "on")]).contains("belongs to the [profile] table"));
        assert!(err(&[perf("pseudo_dir", "/x")]).contains("environment only"));
        assert!(err(&[perf("m2_dgemm", "on"), perf("m2_dgemm", "off")]).contains("given twice"));

        let c = RunConfig::from_input(&[
            (Table::Profile, "phase_profile".into(), "on".into()),
            (Table::Profile, "dgemm_timers".into(), "on".into()),
        ])
        .unwrap();
        assert!(c.on(Id::PHASE_PROFILE));
        assert!(c.present(Id::DGEMM_TIMERS));

        let off = RunConfig::from_pairs(&[("dgemm_timers", "off")]).unwrap();
        assert!(!off.present(Id::DGEMM_TIMERS));
    }

    #[test]
    fn run_config_scopes_are_per_thread_and_follow_the_pool() {
        let a = RunConfig::from_pairs(&[("m2_screen", "replicated")]).unwrap();
        let b = RunConfig::from_pairs(&[("m2_screen", "dist")]).unwrap();
        let pool = crate::pool::WorkerPool::new(3).unwrap();
        let run = |c: RunConfig, want: &'static str| {
            let pool = pool.clone();
            std::thread::spawn(move || {
                let _g = c.enter();
                for _ in 0..200 {
                    assert_eq!(choice(Id::M2_SCREEN), want);
                    let (l, r) = pool.install(|| {
                        crate::pool::WorkerPool::join(
                            || choice(Id::M2_SCREEN),
                            || choice(Id::M2_SCREEN),
                        )
                    });
                    assert_eq!((l, r), (want, want));
                    c.count(RunCounter::SessionSteps, 1);
                }
            })
        };
        let (ta, tb) = (run(a, "replicated"), run(b, "dist"));
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
