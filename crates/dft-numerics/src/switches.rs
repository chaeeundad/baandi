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
    STANDARD_STOP: Kind::Choice { allowed: &["plain", "confirm"], default: "plain" }, "plain",
        "plain, confirm", "standard-class stop rule: plain threshold, or confirm (also stop when the estimator is within 2x the threshold, the previous step within 20x, and it contracted by 0.3 or better)";
    STANDARD_RESIDUAL_GUARD: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "standard-class stop also requires sum_k w_k sum_occupied f_i |r_i|^2 (Ha, second-order band-energy bound) below the threshold";
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
    PHASE_THREADS: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "P48F: thread the atomic start (span build, two-pass Gram-Schmidt right-looking) on the SCF pool; bit-identical";
    // Solvers (dft-solvers).
    B1: Kind::Choice { allowed: &["on", "off", "standard", "strict"], default: "on" }, "on",
        "on, off, standard, strict", "B1 selected-root Davidson iterations";
    B1_SKETCH: Kind::Choice { allowed: &["count", "dense"], default: "count" }, "count", "count, dense",
        "B1 residual sketch";
    B5: Kind::Choice { allowed: &["on", "off", "par"], default: "on" }, "on", "on, off, par",
        "bordered reduced eigenproblem update (par: threaded)";
    B6: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "B6 projector-product reuse";
    RITZ_FUSED: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "fused Ritz rotation";
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
    ORTHO_BLOCK: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "panel Gram-Schmidt of the correction block, L2 shape-rule panel width (P48D; overrides MGS)";
    WIDE_DOT: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "wide dot products in the orthogonalization";
    // Physics (dft-physics).
    MKL_NONLOCAL: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "nonlocal products by dgemm on the real Gamma path";
    APPLY_LEAN: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "lean real Gamma H apply (no redundant scans)";
    HAPPLY_OVERLAP: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "overlap local and nonlocal parts of the H apply";
    HAPPLY_FOLD: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "fold the local potential into the pair FFT";
    NONLOCAL_BLOCK: Kind::Choice { allowed: ON_OFF, default: "off" }, "off", "on, off",
        "complex (k-point) nonlocal of a band block by two GEMMs (P48C)";
    // Numerics (dft-numerics).
    ROTATE_KERNEL: Kind::Choice { allowed: &["legacy", "blocked"], default: "legacy" }, "legacy",
        "legacy, blocked", "rotation kernel of the tiled workspace";
    MKL_DGEMM: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "rotation panels by MKL dgemm (MKL builds)";
    MKL_EVD: Kind::Choice { allowed: ON_OFF, default: "on" }, "on", "on, off",
        "reduced symmetric eigensolve by MKL (MKL builds)";
    PAIR_FFT: Kind::Choice { allowed: &["grid", "pruned"], default: "grid" }, "grid", "grid, pruned",
        "pair FFT plan of MKL builds (pruned: stick lanes)";
    // CLI (dft-cli).
    CUTOFF_LADDER: Kind::Override { valid: cutoff_ladder_valid }, "config", "off, on, fraction in (0, 1)",
        "overrides scf.cutoff_ladder / fraction of the input file";
    RELAX_EXTRAPOLATION: Kind::Override { valid: relax_extrapolation_valid }, "config",
        "none, first_order, second_order", "overrides the relaxation density extrapolation of the input file";
    DGEMM_TIMERS: Kind::Presence, "off", "any value (presence)",
        "print per-site dgemm timers at exit (dgemm-timers builds)";
    PSEUDO_DIR: Kind::Path, "platform cache", "absolute path",
        "pseudopotential cache directory";
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
        assert_eq!(eval(Id::CUTOFF_LADDER, Some("0.5")).1, Env);
        assert_eq!(eval(Id::CUTOFF_LADDER, Some("x")).1, InvalidToError);
        assert_eq!(eval(Id::PSEUDO_DIR, Some("relative")).1, InvalidToError);
        assert_eq!(eval(Id::DGEMM_TIMERS, Some("")).0, Value::Present(true));
    }
}
