//! Observational stall records; never change the SCF trajectory or stop rule.
/// Three consecutive bad transitions suppress isolated condition spikes and
/// one/two-map plateaus, without waiting another history-length cycle.
pub const CONSECUTIVE_MAPS: usize = 3;
pub const CONDITION_LIMIT: f64 = 1e12;
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Diagnostic {
    pub iteration: usize,
    pub channel_gap_ha: [Option<f64>; 2],
    /// Squared Coulomb-norm fraction in an explicitly selected orbital rotation
    /// span. None means unmeasured; it must never be displayed as zero.
    /// A rotation span is not, by itself, proof of a physical zero mode.
    pub orbital_rotation_fraction: Option<f64>,
}
impl Diagnostic {
    pub const RECOMMENDATION: &str = "Ill-conditioned history and no decrease of the stop estimate for three consecutive maps. Consider Fermi-Dirac occupations with fixed moment or a symmetry-breaking start. An unmeasured rotation projection is unavailable, not zero; finite gaps do not exclude a rotational soft mode.";
}
/// Newest first, including the current ordinary map. Redos and verification
/// maps are excluded by callers. No tolerance is used to turn a decrease into
/// a plateau: this is precisely the requested non-decrease predicate.
pub(super) fn detected(mut newest: impl Iterator<Item = (f64, Option<f64>)>) -> bool {
    let Some(mut current) = newest.next() else {
        return false;
    };
    for _ in 0..CONSECUTIVE_MAPS {
        let Some(previous) = newest.next() else {
            return false;
        };
        if !current.0.is_finite()
            || !previous.0.is_finite()
            || current.0 < previous.0
            || !current.1.is_some_and(|c| c > CONDITION_LIMIT)
        {
            return false;
        }
        current = previous;
    }
    true
}
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn needs_three_bad_transitions_and_resets_on_a_decrease_or_good_condition() {
        let bad = Some(2e12);
        assert!(detected(
            [(4., bad), (3., bad), (2., bad), (1., None)].into_iter()
        ));
        assert!(detected([(1., bad); 4].into_iter()));
        assert!(!detected([(1., bad); 3].into_iter()));
        assert!(!detected(
            [(4., bad), (2., bad), (3., bad), (1., None)].into_iter()
        ));
        assert!(!detected(
            [(4., bad), (3., Some(1e12)), (2., bad), (1., None)].into_iter()
        ));
        assert!(!detected([(f64::NAN, bad); 4].into_iter()));
    }
}
