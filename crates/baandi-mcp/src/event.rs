use dft_engine::events::EventKind;
use serde_json::{Value, json};
pub(crate) fn payload(kind: EventKind) -> Value {
    match kind {
        EventKind::Initialized => json!({"type":"initialized"}),
        EventKind::IterationFinished {
            energy_hartree,
            density_residual,
            charge_error_relative,
            mixing_applied,
        } => {
            json!({"type":"iteration_finished","energy_hartree":energy_hartree,"density_residual_bohr_minus_3_over_2":density_residual,"charge_error_relative":charge_error_relative,"mixing_applied":mixing_applied})
        }
        EventKind::CertificateEvaluated {
            passed,
            density_residual,
            occupied_eigenvalue_change_hartree,
        } => {
            json!({"type":"certificate_evaluated","passed":passed,"density_residual_bohr_minus_3_over_2":density_residual,"occupied_eigenvalue_change_hartree":occupied_eigenvalue_change_hartree})
        }
        EventKind::MixerReset => json!({"type":"mixer_reset"}),
        EventKind::Diverged => json!({"type":"diverged"}),
        EventKind::Converged => json!({"type":"converged"}),
        EventKind::NotConverged => json!({"type":"not_converged"}),
        EventKind::Failed { code } => json!({"type":"failed","code":format!("{code:?}")}),
        EventKind::EventsDropped {
            count,
            first_missing_seq,
            last_missing_seq,
        } => {
            json!({"type":"events_dropped","count":count,"first_missing_seq":first_missing_seq,"last_missing_seq":last_missing_seq})
        }
    }
}
