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
    /// XC alone sees the immutable core. Scratch is retained by the caller;
    /// the existing allocating XC output/library boundary is unchanged.
    pub fn evaluate_with_core(
        self,
        n: &Array3<f64>,
        core: Option<&Array3<f64>>,
        scratch: &mut Array3<f64>,
        cell: &Cell,
        fft: &mut FftPlan3d,
    ) -> Result<xc::XcGrid, PhysicsError> {
        let Some(core) = core else {
            return self.evaluate(n, cell, fft);
        };
        if n.shape() != core.shape() || scratch.shape() != n.shape() {
            return Err(PhysicsError::Shape);
        }
        for ((total, &valence), &c) in scratch.iter_mut().zip(n).zip(core) {
            *total = valence + c;
        }
        // Each functional applies its existing screening/negative-density rule
        // to the combined density. Never renormalize or clip the model core.
        let mut result = self.evaluate(scratch, cell, fft)?;
        result.potential_integral = Hartree(
            cell.volume_bohr3() / n.len() as f64
                * n.iter()
                    .zip(&result.potential)
                    .map(|(n, v)| n * v)
                    .sum::<f64>(),
        );
        Ok(result)
    }
    /// Evaluate E_xc[max(n + core, floor)]. The chain rule masks the
    /// potential below the floor after the GGA divergence is assembled.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_with_floor(
        self,
        n: &Array3<f64>,
        core: Option<&Array3<f64>>,
        scratch: &mut Array3<f64>,
        cell: &Cell,
        fft: &mut FftPlan3d,
        floor: Option<f64>,
    ) -> Result<xc::XcGrid, PhysicsError> {
        let Some(floor) = floor else {
            return self.evaluate_with_core(n, core, scratch, cell, fft);
        };
        if !floor.is_finite()
            || floor <= 0.0
            || scratch.shape() != n.shape()
            || core.is_some_and(|c| c.shape() != n.shape())
        {
            return Err(PhysicsError::Shape);
        }
        for (i, (out, &v)) in scratch.iter_mut().zip(n).enumerate() {
            let total = v + core.map_or(0.0, |c| c.as_slice().unwrap()[i]);
            if !total.is_finite() {
                return Err(PhysicsError::Density(total));
            }
            *out = total.max(floor);
        }
        let mut result = self.evaluate(scratch, cell, fft)?;
        for (i, (v, &valence)) in result.potential.iter_mut().zip(n).enumerate() {
            if valence + core.map_or(0.0, |c| c.as_slice().unwrap()[i]) < floor {
                *v = 0.0;
            }
        }
        result.potential_integral = Hartree(
            cell.volume_bohr3() / n.len() as f64
                * n.iter()
                    .zip(&result.potential)
                    .map(|(n, v)| n * v)
                    .sum::<f64>(),
        );
        Ok(result)
    }
}

#[cfg(test)]
#[path = "floor_tests.rs"]
mod floor_tests;
