//! Bounded subspace updates for intermediate SCF maps (ticket
//! BREAKTHROUGH-3-BOUNDED-SUBSPACE, design `docs/design/bounded-subspace-scf.md`).
//!
//! `BAANDI_SCF_ADVANCE=on` asks the standard-path session to run intermediate
//! maps as `SolveIntent::Advance(SubspaceBudget)` on the serial real Gamma lean
//! Davidson: a deterministic number of correction blocks from the kept start,
//! returning the approximate occupied subspace from which the density is
//! built. The stop rule, the residual guard, the verification and the delivered
//! state are evaluated on exact (`Converge`) maps only; a map whose approximate
//! state fails one of the completion rules is completed by the exact solve at
//! the same input density, warm-started from the advanced orbitals. With the
//! switch off nothing here runs (bit-identical to the plain session).
//!
//! Scope of the prototype: closed-shell fixed integer occupations, one exact
//! Gamma point, serial (threads-only) real path with the lean solver. Every
//! other configuration declines with a recorded reason and runs the exact
//! solve. Count-only diagnostics; no allocation or I/O in the solve.
use dft_solvers::real_davidson::{SubspaceBudget, SubspaceDiagnostics};
use std::sync::atomic::{AtomicU64, Ordering::Relaxed};

/// Inner/SCF residual ratio above which an advanced map is completed: the
/// band-energy estimate B = sum_i f_i |r_i|^2 of the advanced state against
/// the Coulomb SCF residual E_H[rho_out - rho_in] of the map. 0.1 is the
/// factor the standard ratchet already uses (tau = 0.1 E / N_e, so the exact
/// path accepts a summed band-energy error of order 0.1 E); not a new tuning
/// constant and not material dependent.
pub const ADVANCE_INNER_RATIO: f64 = 0.1;
/// Bounds of the adapted block budget (fixed integers of the prototype).
pub const ADVANCE_BLOCKS_MIN: usize = 1;
pub const ADVANCE_BLOCKS_MAX: usize = 8;

/// Why an advanced map was completed by an exact solve at the same density.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CompletionTrigger {
    /// E_H[rho_out - rho_in] below the class threshold: the stop rule would
    /// fire, and it may only fire on an exact map.
    Candidate,
    /// B > ADVANCE_INNER_RATIO x E: the inner error is comparable to the SCF
    /// residual (the Pulay secant would learn solver noise).
    InnerNoise,
    /// E grew against the previous map.
    Growth,
    /// Highest occupied Ritz value not below the lowest empty one.
    GapInstability,
}
impl CompletionTrigger {
    pub fn name(self) -> &'static str {
        match self {
            Self::Candidate => "candidate",
            Self::InnerNoise => "inner_noise",
            Self::Growth => "growth",
            Self::GapInstability => "gap_instability",
        }
    }
}

/// Per-iteration record of an advanced map (`IterationRecord.advance`; None
/// on exact maps). Energies in Ha.
#[derive(Clone, Debug, PartialEq)]
pub struct AdvanceRecord {
    pub budget_blocks: usize,
    pub blocks_used: usize,
    /// Operator columns applied by the advance call.
    pub h_columns: usize,
    /// Kept start taken (no initial H application).
    pub kept_start: bool,
    pub unconverged_roots: usize,
    /// B of the advanced state (Ha).
    pub inner_estimate_ha: f64,
    /// E_H[rho_out - rho_in] of the advanced map (Ha).
    pub scf_error_ha: f64,
    pub completion: Option<CompletionTrigger>,
    /// Operator columns of the completing exact solve (0 without completion).
    pub completion_h_columns: usize,
}

/// Request and outcome slot carried through one map (`GammaRealDiagnostics`
/// and `ScfStep`): the session fills `budget`, the solve fills `outcome` or
/// `declined`.
#[derive(Clone, Debug)]
pub struct AdvanceSlot {
    pub budget: SubspaceBudget,
    pub outcome: Option<SubspaceDiagnostics>,
    pub declined: Option<&'static str>,
}
impl AdvanceSlot {
    pub fn new(budget: SubspaceBudget) -> Self {
        Self {
            budget,
            outcome: None,
            declined: None,
        }
    }
}

thread_local! {
    static OVERRIDE: std::cell::Cell<Option<(bool, usize)>> = const { std::cell::Cell::new(None) };
}
/// Override `BAANDI_SCF_ADVANCE` and `BAANDI_SCF_ADVANCE_BLOCKS` for the
/// calling thread (tests); `None` restores the environment.
pub fn set_override(value: Option<(bool, usize)>) {
    OVERRIDE.with(|o| o.set(value));
}
/// `BAANDI_SCF_ADVANCE=on` (read once per process).
pub fn enabled() -> bool {
    OVERRIDE.with(|o| o.get()).map_or_else(
        || {
            static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
            *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::SCF_ADVANCE))
        },
        |(on, _)| on,
    )
}
/// `BAANDI_SCF_ADVANCE_BLOCKS` (integer 1..=8, default 2).
pub fn starting_blocks() -> usize {
    OVERRIDE.with(|o| o.get()).map_or_else(
        || {
            static BLOCKS: std::sync::OnceLock<usize> = std::sync::OnceLock::new();
            *BLOCKS.get_or_init(|| {
                dft_numerics::switches::number(dft_numerics::switches::Id::SCF_ADVANCE_BLOCKS)
                    .map_or(2, |v| v as usize)
                    .clamp(ADVANCE_BLOCKS_MIN, ADVANCE_BLOCKS_MAX)
            })
        },
        |(_, blocks)| blocks.clamp(ADVANCE_BLOCKS_MIN, ADVANCE_BLOCKS_MAX),
    )
}

/// Research sub-switches (read once per process): the growth trigger
/// (`BAANDI_SCF_ADVANCE_GROWTH`), the budget adaptation
/// (`BAANDI_SCF_ADVANCE_ADAPT`) and the scope (`BAANDI_SCF_ADVANCE_SCOPE`:
/// `all` sessions or the `fine` stage only).
pub fn growth_trigger() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::SCF_ADVANCE_GROWTH))
}
pub fn adapt() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| dft_numerics::switches::on(dft_numerics::switches::Id::SCF_ADVANCE_ADAPT))
}
pub fn fine_only() -> bool {
    static ON: std::sync::OnceLock<bool> = std::sync::OnceLock::new();
    *ON.get_or_init(|| {
        dft_numerics::switches::choice(dft_numerics::switches::Id::SCF_ADVANCE_SCOPE) == "fine"
    })
}

/// Budget adaptation from the previous advanced map (section 4 of the design):
/// a completion by inner noise grows the budget by one block; a ratio below a
/// tenth of the rule shrinks it by one; otherwise unchanged.
pub fn adapt_blocks(blocks: usize, ratio: f64, completion: Option<CompletionTrigger>) -> usize {
    if !adapt() {
        return blocks.clamp(ADVANCE_BLOCKS_MIN, ADVANCE_BLOCKS_MAX);
    }
    let next = if completion == Some(CompletionTrigger::InnerNoise) {
        blocks + 1
    } else if ratio <= 0.1 * ADVANCE_INNER_RATIO {
        blocks.saturating_sub(1)
    } else {
        blocks
    };
    next.clamp(ADVANCE_BLOCKS_MIN, ADVANCE_BLOCKS_MAX)
}

#[derive(Clone, Copy)]
#[repr(usize)]
enum Count {
    MapsAdvanced,
    MapsExact,
    Declined,
    CompletionCandidate,
    CompletionInnerNoise,
    CompletionGrowth,
    CompletionGap,
    AdvanceHColumns,
    CompletionHColumns,
    AdvanceBlocks,
    KeptStarts,
    ConvergedWithinBudget,
    RatioBelowRule,
    RatioAboveRule,
    Last,
}
static COUNTS: [AtomicU64; Count::Last as usize] =
    [const { AtomicU64::new(0) }; Count::Last as usize];
fn add(c: Count, v: usize) {
    COUNTS[c as usize].fetch_add(v as u64, Relaxed);
}

/// Process totals (count only), written to results.json
/// `convergence.work.scf_advance`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct AdvanceCounters {
    pub maps_advanced: u64,
    pub maps_exact: u64,
    pub declined: u64,
    pub completion_candidate: u64,
    pub completion_inner_noise: u64,
    pub completion_growth: u64,
    pub completion_gap: u64,
    pub advance_h_columns: u64,
    pub completion_h_columns: u64,
    pub advance_blocks: u64,
    pub kept_starts: u64,
    pub converged_within_budget: u64,
    pub ratio_below_rule: u64,
    pub ratio_above_rule: u64,
}
impl AdvanceCounters {
    pub fn completions(&self) -> u64 {
        self.completion_candidate
            + self.completion_inner_noise
            + self.completion_growth
            + self.completion_gap
    }
}
pub fn counters() -> AdvanceCounters {
    let g = |c: Count| COUNTS[c as usize].load(Relaxed);
    AdvanceCounters {
        maps_advanced: g(Count::MapsAdvanced),
        maps_exact: g(Count::MapsExact),
        declined: g(Count::Declined),
        completion_candidate: g(Count::CompletionCandidate),
        completion_inner_noise: g(Count::CompletionInnerNoise),
        completion_growth: g(Count::CompletionGrowth),
        completion_gap: g(Count::CompletionGap),
        advance_h_columns: g(Count::AdvanceHColumns),
        completion_h_columns: g(Count::CompletionHColumns),
        advance_blocks: g(Count::AdvanceBlocks),
        kept_starts: g(Count::KeptStarts),
        converged_within_budget: g(Count::ConvergedWithinBudget),
        ratio_below_rule: g(Count::RatioBelowRule),
        ratio_above_rule: g(Count::RatioAboveRule),
    }
}
pub(crate) fn count_exact_map() {
    add(Count::MapsExact, 1);
}
pub(crate) fn count_declined() {
    add(Count::Declined, 1);
}
pub(crate) fn count_record(r: &AdvanceRecord, diagnostics: &SubspaceDiagnostics) {
    add(Count::MapsAdvanced, 1);
    add(Count::AdvanceHColumns, r.h_columns);
    add(Count::AdvanceBlocks, r.blocks_used);
    if r.kept_start {
        add(Count::KeptStarts, 1);
    }
    if diagnostics.stop == dft_solvers::real_davidson::BudgetStop::Converged {
        add(Count::ConvergedWithinBudget, 1);
    }
    if r.inner_estimate_ha > ADVANCE_INNER_RATIO * r.scf_error_ha {
        add(Count::RatioAboveRule, 1);
    } else {
        add(Count::RatioBelowRule, 1);
    }
    if let Some(t) = r.completion {
        add(
            match t {
                CompletionTrigger::Candidate => Count::CompletionCandidate,
                CompletionTrigger::InnerNoise => Count::CompletionInnerNoise,
                CompletionTrigger::Growth => Count::CompletionGrowth,
                CompletionTrigger::GapInstability => Count::CompletionGap,
            },
            1,
        );
        add(Count::CompletionHColumns, r.completion_h_columns);
    }
}

/// The completion decision of one advanced map (section 5 of the design),
/// in rule order: candidate, inner noise, growth, gap. `previous_error` is
/// the previous map's E (None on the first map of a history).
pub fn completion_trigger(
    scf_error_ha: f64,
    inner_estimate_ha: f64,
    threshold_ha: f64,
    previous_error_ha: Option<f64>,
    gap_stable: bool,
) -> Option<CompletionTrigger> {
    if !(scf_error_ha.is_finite() && inner_estimate_ha.is_finite()) {
        return Some(CompletionTrigger::InnerNoise);
    }
    if scf_error_ha < threshold_ha {
        Some(CompletionTrigger::Candidate)
    } else if inner_estimate_ha > ADVANCE_INNER_RATIO * scf_error_ha {
        Some(CompletionTrigger::InnerNoise)
    } else if growth_trigger() && previous_error_ha.is_some_and(|p| scf_error_ha > p) {
        Some(CompletionTrigger::Growth)
    } else if !gap_stable {
        Some(CompletionTrigger::GapInstability)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn completion_rules_in_order() {
        // Candidate wins over everything.
        assert_eq!(
            completion_trigger(1e-10, 1.0, 5e-9, Some(1e-12), false),
            Some(CompletionTrigger::Candidate)
        );
        // Inner noise: B above 0.1 E.
        assert_eq!(
            completion_trigger(1e-5, 2e-6, 5e-9, None, true),
            Some(CompletionTrigger::InnerNoise)
        );
        // Growth against the previous map.
        assert_eq!(
            completion_trigger(2e-5, 1e-7, 5e-9, Some(1e-5), true),
            Some(CompletionTrigger::Growth)
        );
        // Gap crossing.
        assert_eq!(
            completion_trigger(1e-5, 1e-7, 5e-9, Some(2e-5), false),
            Some(CompletionTrigger::GapInstability)
        );
        // Accepted.
        assert_eq!(completion_trigger(1e-5, 1e-7, 5e-9, Some(2e-5), true), None);
        // Nonfinite estimates complete.
        assert_eq!(
            completion_trigger(f64::NAN, 1e-7, 5e-9, None, true),
            Some(CompletionTrigger::InnerNoise)
        );
    }

    #[test]
    fn budget_adaptation_is_bounded() {
        assert_eq!(adapt_blocks(2, 0.5, Some(CompletionTrigger::InnerNoise)), 3);
        assert_eq!(
            adapt_blocks(ADVANCE_BLOCKS_MAX, 0.5, Some(CompletionTrigger::InnerNoise)),
            ADVANCE_BLOCKS_MAX
        );
        assert_eq!(adapt_blocks(2, 1e-3, None), 1);
        assert_eq!(adapt_blocks(1, 1e-3, None), 1);
        assert_eq!(adapt_blocks(2, 0.05, None), 2);
        assert_eq!(adapt_blocks(2, 0.05, Some(CompletionTrigger::Growth)), 2);
    }
}
