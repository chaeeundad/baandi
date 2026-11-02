//! Validated UPF v2 NC data and radial transforms in Hartree atomic units.
mod parser;
pub mod radial;
pub use parser::{Projector, Pseudopotential, UpfError, from_file, parse_upf};

/// Whether a UPF label explicitly specifies Slater exchange and PZ correlation
/// with no gradient corrections. Matching is case-sensitive and ignores only
/// ASCII whitespace; the original label is never modified.
pub fn functional_is_lda_pz(functional: &str) -> bool {
    let mut tokens = functional.split_ascii_whitespace();
    if tokens.next() != Some("SLA") || tokens.next() != Some("PZ") {
        return false;
    }
    match tokens.next() {
        None => true,
        Some("NOGX") => tokens.next() == Some("NOGC") && tokens.next().is_none(),
        _ => false,
    }
}

#[cfg(test)]
mod tests {
    use super::functional_is_lda_pz;

    #[test]
    fn explicit_lda_pz_labels() {
        for label in [
            "SLA PZ",
            "SLA PZ NOGX NOGC",
            " SLA  PZ   NOGX NOGC",
            " \tSLA\tPZ\r\n",
            "\r\nSLA\tPZ\r\nNOGX NOGC ",
        ] {
            assert!(functional_is_lda_pz(label), "{label:?}");
        }
        for label in [
            "SLA PW",
            "SLA PZ PBX PBC",
            "PBE",
            "SLA PZ NOGX",
            "SLA PZ NOGC",
            "PZ SLA",
            "SLA PZ NOGC NOGX",
            "",
            " \t\r\n",
            "sla pz",
            "SLA PZ NOGX NOGC EXTRA",
            "SLA PZ EXTRA",
            "SLA\u{00a0}PZ",
        ] {
            assert!(!functional_is_lda_pz(label), "{label:?}");
        }
    }
}
