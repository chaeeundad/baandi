//! Closed XC selection shared by input checks and both energy evaluations.
use dft_geometry::Cell;
use dft_physics::{PhysicsError, data::*, gga, xc};

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Functional {
    LdaPz,
    LdaPw,
    Pbe,
}
impl Functional {
    pub fn name(self) -> &'static str {
        match self {
            Self::LdaPz => "lda-pz",
            Self::LdaPw => "lda-pw",
            Self::Pbe => "pbe",
        }
    }
    pub fn matches_upf(self, label: &str) -> bool {
        match self {
            Self::LdaPz => dft_upf::functional_is_lda_pz(label),
            Self::LdaPw => matches!(
                label
                    .split_whitespace()
                    .collect::<Vec<_>>()
                    .join(" ")
                    .to_ascii_uppercase()
                    .as_str(),
                "SLA PW" | "SLA PW NOGX NOGC"
            ),
            Self::Pbe => matches!(
                label
                    .split_whitespace()
                    .map(str::to_ascii_uppercase)
                    .collect::<Vec<_>>()
                    .join(" ")
                    .as_str(),
                "PBE" | "SLA PW PBX PBC"
            ),
        }
    }
    pub fn evaluate(
        self,
        n: &Array3<f64>,
        cell: &Cell,
        fft: &mut FftPlan3d,
    ) -> Result<xc::XcGrid, PhysicsError> {
        match self {
            Self::LdaPz => xc::evaluate(n, cell.volume_bohr3()),
            Self::LdaPw => dft_physics::spin::unpolarized_pw(n, cell.volume_bohr3()),
            Self::Pbe => gga::evaluate(n, cell, fft),
        }
    }
}
