//! Opt-in, deterministic boundary diagnostics; never convergence or control authority.
//!
//! Implements the proposed v1 catalogue in `docs/research/self-healing-scf.md`.
//! Six residuals above t stagnate when max/min <= 1.2 and last/first >= 0.95.
//! Six energy differences (seven records) oscillate when all magnitudes exceed
//! the energy tolerance, all five adjacent signs differ, and r_last/r_second >=
//! 0.9. The range denominator uses epsilon=max(t*1e-6,1e-30); a zero
//! endpoint denominator gives insufficient evidence. Divergence requires the
//! configured number of successive r_i > factor*max(prior_best,t).
//! These dimensionless thresholds are proposed advisory definitions, not tuned
//! engine defaults or universally calibrated classifiers. Six samples reject
//! short transients; the range and endpoint checks distinguish a plateau from
//! consistent progress. See the release integration test for measured scope.
//!
//! Supply one calculation's authoritative history, starting at iteration 1;
//! tag each sample with the active policy revision. Windows reset on revision
//! changes. Gaps, duplicate/reordered iterations, or nonfinite values suppress
//! advisory conclusions. A lossy event stream cannot replace this history.
//! Fourier probes are optional caller-supplied measurements, never inferred
//! from a scalar residual. No probe instrumentation is installed by this module.
//! Candidate intents are deliberately non-executable: P11B policy types are
//! not yet available. No beta, tolerance, solver, or mixer is modified here.
use crate::events::{EventKind, FailureCode, ScfEvent};
use crate::{IterationRecord, ScfOptions};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rule {
    NonFinite,
    Divergence,
    Stagnation,
    Oscillation,
    ChargeSloshingHint,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Status {
    Triggered,
    Clear,
    InsufficientEvidence,
    NotComputed,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Hypothesis {
    ArithmeticFailure,
    UnstableIteration,
    PoorMixing,
    AlternatingResponse,
    LongWavelengthMode,
}
/// Advisory intentions only, pending integration with P11B validation/approval.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Candidate {
    ReduceBetaAndResetHistory,
    ReviewNewCalculation,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Field {
    EnergyHartree,
    DensityResidual,
    ChargeErrorRelative,
    MinDensity,
    EigenResidualHartree,
    InputChargeErrorRelative,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NonFiniteClass {
    Nan,
    PositiveInfinity,
    NegativeInfinity,
}
#[derive(Clone, Debug, PartialEq)]
pub enum Evidence {
    Unavailable,
    /// The event has no field/value detail; never parse numerical error text.
    TerminalFailure {
        state_version: u64,
    },
    FiniteScalars {
        count: usize,
    },
    NonFinite {
        iteration: usize,
        field: Field,
        classification: NonFiniteClass,
    },
    Stagnation {
        residuals: [f64; 6],
        range_ratio: f64,
        endpoint_ratio: f64,
    },
    Oscillation {
        delta_energy_hartree: [f64; 6],
        sign_changes: usize,
        residual_ratio: f64,
    },
    Divergence {
        residual: f64,
        prior_best: f64,
        factor: f64,
        count: usize,
    },
    Sloshing {
        fractions: [f64; 3],
        q_min_inverse_bohr: f64,
        basis_id: String,
        iterations: [usize; 3],
    },
}
#[derive(Clone, Debug, PartialEq)]
pub struct Finding {
    pub rule: Rule,
    pub rule_version: &'static str,
    pub status: Status,
    pub iteration_first: usize,
    pub iteration_last: usize,
    pub policy_revision: u64,
    pub evidence: Evidence,
    /// An uncalibrated hypothesis, never an observed cause.
    pub hypothesis: Option<Hypothesis>,
    pub candidates: Vec<Candidate>,
}
/// Scalars only; optional values were not present in routine P11A events.
#[derive(Clone, Debug)]
pub struct Sample {
    pub iteration: usize,
    pub policy_revision: u64,
    pub energy_hartree: f64,
    pub density_residual: f64,
    pub charge_error_relative: f64,
    pub min_density: Option<f64>,
    pub eigen_residual_hartree: Option<f64>,
    pub input_charge_error_relative: Option<f64>,
}
impl Sample {
    pub fn from_record(r: &IterationRecord, policy_revision: u64) -> Self {
        Self {
            iteration: r.iteration,
            policy_revision,
            energy_hartree: r.energy.total_hartree(),
            density_residual: r.density_residual,
            charge_error_relative: r.charge_error,
            min_density: Some(r.min_density),
            eigen_residual_hartree: Some(r.max_eigen_residual),
            input_charge_error_relative: Some(r.input_charge_error),
        }
    }
    /// Caller must keep calculation identities separate and supply the revision
    /// from the authoritative policy log (P11A events do not carry it).
    pub fn from_event(event: &ScfEvent, policy_revision: u64) -> Option<Self> {
        if let EventKind::IterationFinished {
            energy_hartree,
            density_residual,
            charge_error_relative,
            ..
        } = event.kind
        {
            Some(Self {
                iteration: event.iteration,
                policy_revision,
                energy_hartree,
                density_residual,
                charge_error_relative,
                min_density: None,
                eigen_residual_hartree: None,
                input_charge_error_relative: None,
            })
        } else {
            None
        }
    }
    fn scalars(&self) -> [(Field, Option<f64>); 6] {
        [
            (Field::EnergyHartree, Some(self.energy_hartree)),
            (Field::DensityResidual, Some(self.density_residual)),
            (Field::ChargeErrorRelative, Some(self.charge_error_relative)),
            (Field::MinDensity, self.min_density),
            (Field::EigenResidualHartree, self.eigen_residual_hartree),
            (
                Field::InputChargeErrorRelative,
                self.input_charge_error_relative,
            ),
        ]
    }
}
/// Explicit Fourier residual power measurements: exclude G=0 from both sums;
/// low_q_power sums 0 < |G| <= 2*q_min; total_power sums all |G| > 0.
#[derive(Clone, Debug)]
pub struct FourierProbe {
    pub iteration: usize,
    pub policy_revision: u64,
    pub low_q_power: f64,
    pub total_power: f64,
    pub q_min_inverse_bohr: f64,
    pub basis_id: String,
}
#[derive(Clone, Copy, Debug)]
pub struct Thresholds {
    pub density_tolerance: f64,
    pub energy_tolerance_hartree: f64,
    pub divergence_factor: f64,
    pub divergence_patience: usize,
}
impl From<&ScfOptions> for Thresholds {
    fn from(o: &ScfOptions) -> Self {
        Self {
            density_tolerance: o.density_tolerance,
            energy_tolerance_hartree: o.energy_tolerance,
            divergence_factor: o.divergence_factor,
            divergence_patience: o.divergence_patience,
        }
    }
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct InvalidThresholds;

/// Evaluate the latest policy segment. Output order is fixed: errors before
/// advisory rules. NaN/Inf are represented by classification, never in evidence.
/// No timestamps, randomness, I/O, or engine state mutations are used.
pub fn diagnose(
    history: &[Sample],
    probes: &[FourierProbe],
    t: Thresholds,
) -> Result<Vec<Finding>, InvalidThresholds> {
    if !t.density_tolerance.is_finite()
        || t.density_tolerance <= 0.0
        || !t.energy_tolerance_hartree.is_finite()
        || t.energy_tolerance_hartree <= 0.0
        || !t.divergence_factor.is_finite()
        || t.divergence_factor <= 1.0
        || t.divergence_patience == 0
    {
        return Err(InvalidThresholds);
    }
    let last = history.last();
    let revision = last.map_or(0, |s| s.policy_revision);
    let start = history
        .iter()
        .rposition(|s| s.policy_revision != revision)
        .map_or(0, |i| i + 1);
    let h = &history[start..];
    let make = |rule| Finding {
        rule,
        rule_version: "1.0",
        status: Status::InsufficientEvidence,
        iteration_first: h.first().map_or(0, |s| s.iteration),
        iteration_last: last.map_or(0, |s| s.iteration),
        policy_revision: revision,
        evidence: Evidence::Unavailable,
        hypothesis: None,
        candidates: vec![],
    };
    let mut out: Vec<_> = [
        Rule::NonFinite,
        Rule::Divergence,
        Rule::Stagnation,
        Rule::Oscillation,
        Rule::ChargeSloshingHint,
    ]
    .into_iter()
    .map(make)
    .collect();
    out[4].status = Status::NotComputed;
    if h.is_empty() {
        return Ok(out);
    }
    let mut count = 0;
    for s in h {
        for (field, value) in s.scalars() {
            if let Some(v) = value {
                count += 1;
                if !v.is_finite() {
                    let classification = if v.is_nan() {
                        NonFiniteClass::Nan
                    } else if v.is_sign_positive() {
                        NonFiniteClass::PositiveInfinity
                    } else {
                        NonFiniteClass::NegativeInfinity
                    };
                    out[0].evidence = Evidence::NonFinite {
                        iteration: s.iteration,
                        field,
                        classification,
                    };
                    trigger(
                        &mut out[0],
                        Hypothesis::ArithmeticFailure,
                        Candidate::ReviewNewCalculation,
                    );
                    return Ok(out);
                }
            }
        }
    }
    out[0].status = Status::Clear;
    out[0].evidence = Evidence::FiniteScalars { count };
    if h[0].iteration == 0
        || h.windows(2)
            .any(|w| w[0].iteration.checked_add(1) != Some(w[1].iteration))
        || h.iter().any(|s| s.density_residual < 0.0)
    {
        return Ok(out);
    }
    let eps = (t.density_tolerance * 1e-6).max(1e-30);
    // A prefix omitted from the initial revision could hide a lower prior best.
    if h.len() >= 2
        && (start > 0 && history[start - 1].iteration.checked_add(1) == Some(h[0].iteration)
            || h[0].iteration == 1)
    {
        let mut best = h[0].density_residual;
        let mut count = 0;
        let mut prior_best = best;
        for s in &h[1..] {
            prior_best = best;
            count = if s.density_residual > t.divergence_factor * best.max(t.density_tolerance) {
                count + 1
            } else {
                0
            };
            best = best.min(s.density_residual);
        }
        out[1].status = Status::Clear;
        out[1].evidence = Evidence::Divergence {
            residual: h.last().unwrap().density_residual,
            prior_best,
            factor: t.divergence_factor,
            count,
        };
        if count >= t.divergence_patience {
            trigger(
                &mut out[1],
                Hypothesis::UnstableIteration,
                Candidate::ReviewNewCalculation,
            );
        }
    }
    if h.len() >= 6 {
        let w = &h[h.len() - 6..];
        let residuals = std::array::from_fn(|i| w[i].density_residual);
        let min = residuals.iter().copied().fold(f64::INFINITY, f64::min);
        let max = residuals.iter().copied().fold(0.0, f64::max);
        let range_ratio = max / min.max(eps);
        let endpoint_ratio = if residuals[0] > 0.0 {
            residuals[5] / residuals[0]
        } else {
            f64::NAN
        };
        if range_ratio.is_finite() && endpoint_ratio.is_finite() {
            out[2].iteration_first = w[0].iteration;
            out[2].status = Status::Clear;
            out[2].evidence = Evidence::Stagnation {
                residuals,
                range_ratio,
                endpoint_ratio,
            };
            if min > t.density_tolerance && range_ratio <= 1.2 && endpoint_ratio >= 0.95 {
                trigger(
                    &mut out[2],
                    Hypothesis::PoorMixing,
                    Candidate::ReduceBetaAndResetHistory,
                );
            }
        }
    }
    if h.len() >= 7 {
        let w = &h[h.len() - 7..];
        let delta_energy_hartree: [f64; 6] =
            std::array::from_fn(|i| w[i + 1].energy_hartree - w[i].energy_hartree);
        let residual_ratio = if w[1].density_residual > 0.0 {
            w[6].density_residual / w[1].density_residual
        } else {
            f64::NAN
        };
        if residual_ratio.is_finite() && delta_energy_hartree.iter().all(|d| d.is_finite()) {
            let sign_changes = delta_energy_hartree
                .windows(2)
                .filter(|d| {
                    d[0] != 0.0 && d[1] != 0.0 && d[0].is_sign_positive() != d[1].is_sign_positive()
                })
                .count();
            out[3].iteration_first = w[0].iteration;
            out[3].status = Status::Clear;
            out[3].evidence = Evidence::Oscillation {
                delta_energy_hartree,
                sign_changes,
                residual_ratio,
            };
            if sign_changes == 5
                && residual_ratio >= 0.9
                && delta_energy_hartree
                    .iter()
                    .all(|d| d.abs() > t.energy_tolerance_hartree)
            {
                trigger(
                    &mut out[3],
                    Hypothesis::AlternatingResponse,
                    Candidate::ReduceBetaAndResetHistory,
                );
            }
        }
    }
    if !probes.is_empty() {
        out[4].status = Status::InsufficientEvidence;
        if probes.len() >= 3 && h.len() >= 3 {
            let p = &probes[probes.len() - 3..];
            let w = &h[h.len() - 3..];
            let valid = p.iter().zip(w).all(|(p, s)| {
                p.iteration == s.iteration
                    && p.policy_revision == revision
                    && p.total_power.is_finite()
                    && p.total_power > 0.0
                    && p.low_q_power.is_finite()
                    && p.low_q_power >= 0.0
                    && p.low_q_power <= p.total_power
                    && p.q_min_inverse_bohr.is_finite()
                    && p.q_min_inverse_bohr > 0.0
                    && p.q_min_inverse_bohr == probes[probes.len() - 1].q_min_inverse_bohr
                    && !p.basis_id.is_empty()
                    && p.basis_id == probes[probes.len() - 1].basis_id
            });
            if valid && matches!(out[2].status, Status::Clear | Status::Triggered) {
                let fractions = std::array::from_fn(|i| p[i].low_q_power / p[i].total_power);
                out[4].iteration_first = out[2].iteration_first;
                out[4].status = Status::Clear;
                out[4].evidence = Evidence::Sloshing {
                    fractions,
                    q_min_inverse_bohr: p[0].q_min_inverse_bohr,
                    basis_id: p[0].basis_id.clone(),
                    iterations: std::array::from_fn(|i| p[i].iteration),
                };
                if out[2].status == Status::Triggered && fractions.iter().all(|f| *f >= 0.8) {
                    trigger(
                        &mut out[4],
                        Hypothesis::LongWavelengthMode,
                        Candidate::ReduceBetaAndResetHistory,
                    );
                }
            }
        }
    }
    Ok(out)
}
fn trigger(f: &mut Finding, hypothesis: Hypothesis, candidate: Candidate) {
    f.status = Status::Triggered;
    f.hypothesis = Some(hypothesis);
    f.candidates.push(candidate);
}

/// Preserve an authoritative terminal event even when no iteration completed.
/// Other failures (including untyped solver errors) are not inferred here.
pub fn diagnose_failure(event: &ScfEvent, policy_revision: u64) -> Option<Finding> {
    let (rule, hypothesis) = match event.kind {
        EventKind::Failed {
            code: FailureCode::NonFinite,
        } => (Rule::NonFinite, Hypothesis::ArithmeticFailure),
        EventKind::Diverged
        | EventKind::Failed {
            code: FailureCode::Diverged,
        } => (Rule::Divergence, Hypothesis::UnstableIteration),
        _ => return None,
    };
    Some(Finding {
        rule,
        rule_version: "1.0",
        status: Status::Triggered,
        iteration_first: event.iteration,
        iteration_last: event.iteration,
        policy_revision,
        evidence: Evidence::TerminalFailure {
            state_version: event.state_version,
        },
        hypothesis: Some(hypothesis),
        candidates: vec![Candidate::ReviewNewCalculation],
    })
}
