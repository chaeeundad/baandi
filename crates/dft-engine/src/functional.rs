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
            // UPF writers separate the four tokens with blanks or hyphens
            // ("SLA PW PBX PBC", "SLA-PW-PBX-PBC" in HGH files).
            Self::Pbe => matches!(
                label
                    .split(|c: char| c.is_whitespace() || c == '-')
                    .filter(|t| !t.is_empty())
                    .map(str::to_ascii_uppercase)
                    .collect::<Vec<_>>()
                    .join(" ")
                    .as_str(),
                "PBE" | "SLA PW PBX PBC" | "SLA PW PBE PBE"
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
        let _scope = xc::screen_negative().then(|| {
            xc::screening_scope(
                xc::screening_stage(),
                xc::screening_prechecked(),
                Some(cell.volume_bohr3() / n.len() as f64 * n.sum()),
            )
        });
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

pub use dft_physics::pbe_pw92::Pw92Prefactor;

impl crate::ScfOptions {
    /// Default delegates directly to the unchanged library-based implementation.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn evaluate_xc(
        &self,
        n: &Array3<f64>,
        core: Option<&Array3<f64>>,
        scratch: &mut Array3<f64>,
        cell: &Cell,
        fft: &mut FftPlan3d,
        floor: Option<f64>,
    ) -> Result<xc::XcGrid, PhysicsError> {
        if self.pw92_prefactor == Pw92Prefactor::Modified {
            return self
                .functional
                .evaluate_with_floor(n, core, scratch, cell, fft, floor);
        }
        if self.functional != Functional::Pbe {
            return Err(PhysicsError::Xc(
                "pw92_prefactor=original requires unpolarized PBE".into(),
            ));
        }
        if scratch.shape() != n.shape()
            || core.is_some_and(|c| c.shape() != n.shape())
            || floor.is_some_and(|f| !f.is_finite() || f <= 0.0)
        {
            return Err(PhysicsError::Shape);
        }
        for (i, (out, &v)) in scratch.iter_mut().zip(n).enumerate() {
            let total = v + core.map_or(0.0, |c| c.as_slice().unwrap()[i]);
            if !total.is_finite() {
                return Err(PhysicsError::Density(total));
            }
            *out = floor.map_or(total, |f| total.max(f));
        }
        let _scope = xc::screen_negative().then(|| {
            xc::screening_scope(
                xc::screening_stage(),
                xc::screening_prechecked(),
                Some(cell.volume_bohr3() / n.len() as f64 * n.sum()),
            )
        });
        let mut result = gga::evaluate_pw92(self.pw92_prefactor, scratch, cell, fft)?;
        if let Some(floor) = floor {
            for (i, (v, &valence)) in result.potential.iter_mut().zip(n).enumerate() {
                if valence + core.map_or(0.0, |c| c.as_slice().unwrap()[i]) < floor {
                    *v = 0.0;
                }
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

/// GD2 (`DERIV_XC_SPLIT`): splits the pointwise PBE evaluation of the
/// derivative pass over the ranks of the active communicator: contiguous
/// point blocks (`dft_mpi::block_counts`), one allgatherv with a status value
/// per rank. Every gathered value is computed by exactly one rank with the
/// serial operations, so the result does not depend on the rank count.
struct CommSplit(crate::parallel::SharedComm);
impl gga::PointSplit for CommSplit {
    fn block(&self, len: usize) -> std::ops::Range<usize> {
        use dft_mpi::Collective;
        let counts = dft_mpi::block_counts(len, self.0.size());
        let first = dft_mpi::block_displacements(&counts)[self.0.rank()];
        first..first + counts[self.0.rank()]
    }
    fn gather(
        &self,
        local: &[f64],
        len: usize,
        width: usize,
        failed: bool,
    ) -> Result<Vec<f64>, PhysicsError> {
        use dft_mpi::Collective;
        let counts = dft_mpi::block_counts(len, self.0.size());
        let mut send = Vec::with_capacity(1 + counts[self.0.rank()] * width);
        send.push(f64::from(u8::from(failed)));
        if failed {
            send.resize(1 + counts[self.0.rank()] * width, 0.);
        } else {
            if local.len() != counts[self.0.rank()] * width {
                return Err(PhysicsError::Shape);
            }
            send.extend_from_slice(local);
        }
        let counts: Vec<_> = counts.iter().map(|c| 1 + c * width).collect();
        let mut all = vec![0.; counts.iter().sum()];
        self.0
            .allgatherv_f64(&send, &counts, &mut all)
            .map_err(|e| PhysicsError::Xc(format!("derivative XC gather: {e}")))?;
        let mut out = Vec::with_capacity(len * width);
        let mut at = 0;
        for count in counts {
            if all[at] != 0. {
                return Err(PhysicsError::Xc(
                    "derivative XC evaluation failed on a rank".into(),
                ));
            }
            out.extend_from_slice(&all[at + 1..at + count]);
            at += count;
        }
        Ok(out)
    }
}

/// XC result and, when available, the GGA stress metric of the same density.
pub(crate) type DerivativeXc = (xc::XcGrid, Option<[[f64; 3]; 3]>);

impl crate::ScfOptions {
    /// GD2: `evaluate_xc` for the post-SCF derivatives (NLCC core forces, XC
    /// stress). With `DERIV_XC_SPLIT` on and unpolarized PBE it runs
    /// `gga::evaluate_for_derivatives` with the same surrounding operations
    /// as `evaluate_xc` (core addition, floor, screening scope, masking and
    /// valence integral), bit-identical to it, and with `metric` also returns
    /// the GGA stress metric of the same density when available. Otherwise it
    /// is `evaluate_xc` and the metric is `None`.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn evaluate_xc_for_derivatives(
        &self,
        n: &Array3<f64>,
        core: Option<&Array3<f64>>,
        scratch: &mut Array3<f64>,
        cell: &Cell,
        fft: &mut FftPlan3d,
        floor: Option<f64>,
        metric: bool,
    ) -> Result<DerivativeXc, PhysicsError> {
        if self.functional != Functional::Pbe
            || !dft_numerics::switches::on(dft_numerics::switches::Id::DERIV_XC_SPLIT)
        {
            return self
                .evaluate_xc(n, core, scratch, cell, fft, floor)
                .map(|r| (r, None));
        }
        let comm = crate::parallel::SharedComm::active().map(CommSplit);
        let split = comm.as_ref().map(|c| c as &dyn gga::PointSplit);
        let prefactor =
            (self.pw92_prefactor != Pw92Prefactor::Modified).then_some(self.pw92_prefactor);
        let run = |input: &Array3<f64>, fft: &mut FftPlan3d| {
            gga::evaluate_for_derivatives(prefactor, input, cell, fft, split, metric)
                .map(|r| (r.grid, r.metric))
        };
        let valence_integral = |result: &mut xc::XcGrid| {
            result.potential_integral = Hartree(
                cell.volume_bohr3() / n.len() as f64
                    * n.iter()
                        .zip(&result.potential)
                        .map(|(n, v)| n * v)
                        .sum::<f64>(),
            );
        };
        if prefactor.is_some() {
            // `evaluate_xc` with the original PW92 prefactor.
            if scratch.shape() != n.shape()
                || core.is_some_and(|c| c.shape() != n.shape())
                || floor.is_some_and(|f| !f.is_finite() || f <= 0.0)
            {
                return Err(PhysicsError::Shape);
            }
            for (i, (out, &v)) in scratch.iter_mut().zip(n).enumerate() {
                let total = v + core.map_or(0.0, |c| c.as_slice().unwrap()[i]);
                if !total.is_finite() {
                    return Err(PhysicsError::Density(total));
                }
                *out = floor.map_or(total, |f| total.max(f));
            }
            let _scope = xc::screen_negative().then(|| {
                xc::screening_scope(
                    xc::screening_stage(),
                    xc::screening_prechecked(),
                    Some(cell.volume_bohr3() / n.len() as f64 * n.sum()),
                )
            });
            let (mut result, metric) = run(scratch, fft)?;
            if let Some(floor) = floor {
                for (i, (v, &valence)) in result.potential.iter_mut().zip(n).enumerate() {
                    if valence + core.map_or(0.0, |c| c.as_slice().unwrap()[i]) < floor {
                        *v = 0.0;
                    }
                }
            }
            valence_integral(&mut result);
            return Ok((result, metric));
        }
        // `Functional::evaluate_with_floor` (modified prefactor).
        if let Some(floor) = floor {
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
            let (mut result, metric) = run(scratch, fft)?;
            for (i, (v, &valence)) in result.potential.iter_mut().zip(n).enumerate() {
                if valence + core.map_or(0.0, |c| c.as_slice().unwrap()[i]) < floor {
                    *v = 0.0;
                }
            }
            valence_integral(&mut result);
            return Ok((result, metric));
        }
        // `Functional::evaluate_with_core`.
        let Some(core) = core else {
            return run(n, fft);
        };
        if n.shape() != core.shape() || scratch.shape() != n.shape() {
            return Err(PhysicsError::Shape);
        }
        for ((total, &valence), &c) in scratch.iter_mut().zip(n).zip(core) {
            *total = valence + c;
        }
        let _scope = xc::screen_negative().then(|| {
            xc::screening_scope(
                xc::screening_stage(),
                xc::screening_prechecked(),
                Some(cell.volume_bohr3() / n.len() as f64 * n.sum()),
            )
        });
        let (mut result, metric) = run(scratch, fft)?;
        valence_integral(&mut result);
        Ok((result, metric))
    }
}
