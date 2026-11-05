//! Direct orbital functional and input-density Harris-Foulkes audit, in Ha.
use crate::{Functional, ScfError};
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_physics::{data::*, hamiltonian::NonlocalOperator, hartree, kinetic, xc};

#[derive(Clone, Debug)]
pub struct EnergyAudit {
    pub xc_floor: Option<XcFloorAudit>,
    /// Explicit audit-only core contribution; XC energy includes core, all band
    /// double counting excludes it. None preserves the legacy model/ledger.
    pub nlcc: Option<CoreXcLedger>,
    pub kinetic: f64,
    pub local: f64,
    pub local_zero_mode: f64,
    pub nonlocal: f64,
    pub hartree: f64,
    pub xc: f64,
    pub ewald: f64,
    pub band_sum: f64,
    /// Integral n_out (V_H[n_in]+V_xc[n_in]) d^3r.
    pub double_counting: f64,
    /// E_band - double_counting + E_H[n_out] + E_xc[n_out] + E_ewald.
    pub band_form: f64,
    pub harris_foulkes: f64,
    /// Bound using exact Hartree remainder, absolute XC remainder, and
    /// sum f_b ||H c_b-eps_b c_b|| (normalized orbitals).
    pub functional_difference_bound: f64,
}
#[derive(Clone, Debug)]
pub struct XcFloorAudit {
    pub density_floor: f64,
    pub input_negative_charge: f64,
    pub output_negative_charge: f64,
}
#[derive(Clone, Debug)]
pub struct CoreXcLedger {
    pub core_electrons: f64,
    pub core_potential_integral_ha: f64,
    pub valence_potential_integral_ha: f64,
}
pub(crate) fn core_ledger(
    core: Option<&Array3<f64>>,
    valence: &Array3<f64>,
    x: &xc::XcGrid,
    volume: f64,
) -> Option<CoreXcLedger> {
    core.map(|core| {
        let weight = volume / valence.len() as f64;
        CoreXcLedger {
            core_electrons: core.sum() * weight,
            core_potential_integral_ha: weight
                * core
                    .iter()
                    .zip(&x.potential)
                    .map(|(c, v)| c * v)
                    .sum::<f64>(),
            valence_potential_integral_ha: x.potential_integral.0,
        }
    })
}
impl EnergyAudit {
    pub fn total_hartree(&self) -> f64 {
        self.kinetic + self.local + self.nonlocal + self.hartree + self.xc + self.ewald
    }
    pub fn total_rydberg(&self) -> f64 {
        2.0 * self.total_hartree()
    }
}

/// Full-grid Parseval integral, including G=0, with the production FFT convention.
pub fn local_energy(
    n: &Array3<f64>,
    ionic_g: &Array3<Complex64>,
    volume: f64,
    fft: &mut FftPlan3d,
) -> Result<f64, ScfError> {
    if n.shape() != ionic_g.shape() || !volume.is_finite() || volume <= 0.0 {
        return Err(ScfError::InvalidInput("local energy shape/volume"));
    }
    let mut ng = Array3::zeros(n.raw_dim());
    fft.to_reciprocal(&n.mapv(|n| Complex64::new(n, 0.0)), &mut ng)?;
    let e = volume
        * ng.iter()
            .zip(ionic_g)
            .map(|(n, v)| (n.conj() * v).re)
            .sum::<f64>();
    if !e.is_finite() {
        return Err(ScfError::NonFinite);
    }
    Ok(e)
}

/// All orbitals are normalized in basis order. Input fields must be those used
/// for the supplied eigenpairs; output density must be built from these orbitals.
pub struct AuditInput<'a> {
    pub cell: &'a Cell,
    pub basis: &'a PlaneWaveBasis,
    pub grid: &'a FftGridSpec,
    pub ionic_g: &'a Array3<Complex64>,
    pub nonlocal: &'a NonlocalOperator,
    pub orbitals: &'a [Vec<Complex64>],
    pub occupations: &'a [f64],
    pub eigenvalues: &'a [f64],
    pub eigen_residuals: &'a [f64],
    pub input_density: &'a Array3<f64>,
    pub output_density: &'a Array3<f64>,
    pub input_hartree: &'a hartree::HartreeGrid,
    pub input_xc: &'a xc::XcGrid,
    pub ewald: f64,
}
pub fn evaluate(a: AuditInput<'_>, fft: &mut FftPlan3d) -> Result<EnergyAudit, ScfError> {
    evaluate_with_functional(a, fft, Functional::LdaPz)
}
pub fn evaluate_with_functional(
    a: AuditInput<'_>,
    fft: &mut FftPlan3d,
    functional: Functional,
) -> Result<EnergyAudit, ScfError> {
    evaluate_with_core_xc(a, fft, functional, None, None)
}
pub(crate) fn evaluate_with_core_xc(
    a: AuditInput<'_>,
    fft: &mut FftPlan3d,
    functional: Functional,
    output_xc: Option<&xc::XcGrid>,
    core: Option<&Array3<f64>>,
) -> Result<EnergyAudit, ScfError> {
    let n = a.output_density;
    if a.orbitals.len() != a.occupations.len()
        || a.eigenvalues.len() != a.occupations.len()
        || a.eigen_residuals.len() != a.occupations.len()
        || n.shape() != a.grid.dims()
        || a.input_density.shape() != n.shape()
        || a.input_hartree.potential.shape() != n.shape()
        || a.input_xc.potential.shape() != n.shape()
    {
        return Err(ScfError::InvalidInput("energy dimensions"));
    }
    if a.eigenvalues.iter().any(|v| !v.is_finite())
        || a.eigen_residuals.iter().any(|r| !r.is_finite() || *r < 0.0)
        || !a.ewald.is_finite()
        || !a.input_hartree.energy.0.is_finite()
        || !a.input_xc.energy.0.is_finite()
        || !a.input_xc.potential_integral.0.is_finite()
    {
        return Err(ScfError::InvalidInput("invalid energy/eigenpair data"));
    }
    let volume = a.cell.volume_bohr3();
    let weight = volume / n.len() as f64;
    let h = hartree::evaluate(n, a.cell, fft)?;
    let owned_xc;
    let x = if let Some(x) = output_xc {
        x
    } else {
        owned_xc = functional.evaluate(n, a.cell, fft)?;
        &owned_xc
    };
    let kinetic = kinetic::energy(a.basis, a.orbitals, a.occupations)?.0;
    let local = local_energy(n, a.ionic_g, volume, fft)?;
    let charge = n.sum() * weight;
    let mut nonlocal = 0.0;
    for (c, f) in a.orbitals.iter().zip(a.occupations) {
        let mut v = vec![Complex64::default(); c.len()];
        a.nonlocal.apply(c, &mut v)?;
        nonlocal += f * c.iter().zip(v).map(|(c, v)| (c.conj() * v).re).sum::<f64>();
    }
    let band_sum = a
        .eigenvalues
        .iter()
        .zip(a.occupations)
        .map(|(e, f)| e * f)
        .sum::<f64>();
    let double_counting = weight
        * n.iter()
            .zip(&a.input_hartree.potential)
            .zip(&a.input_xc.potential)
            .map(|((n, h), x)| n * (h.re + x))
            .sum::<f64>();
    let harris_foulkes = band_sum - a.input_hartree.energy.0 + a.input_xc.energy.0
        - a.input_xc.potential_integral.0
        + a.ewald;
    let band_form = band_sum - double_counting + h.energy.0 + x.energy.0 + a.ewald;
    // Exact discrete functional remainders apply to LDA and GGA alike.
    // No local-density Hessian or quadratic GGA bound is assumed.
    // Bregman remainders. Hartree is quadratic, so its remainder is >=0.
    let linear_h = weight
        * n.iter()
            .zip(a.input_density)
            .zip(&a.input_hartree.potential)
            .map(|((o, i), h)| (o - i) * h.re)
            .sum::<f64>();
    let linear_x = weight
        * n.iter()
            .zip(a.input_density)
            .zip(&a.input_xc.potential)
            .map(|((o, i), x)| (o - i) * x)
            .sum::<f64>();
    let rh = h.energy.0 - a.input_hartree.energy.0 - linear_h;
    let rx = x.energy.0 - a.input_xc.energy.0 - linear_x;
    let eig_bound = a
        .eigen_residuals
        .iter()
        .zip(a.occupations)
        .map(|(r, f)| r * f)
        .sum::<f64>();
    let audit = EnergyAudit {
        xc_floor: None,
        nlcc: core_ledger(core, n, x, volume),
        kinetic,
        local,
        local_zero_mode: charge * a.ionic_g[[0, 0, 0]].re,
        nonlocal,
        hartree: h.energy.0,
        xc: x.energy.0,
        ewald: a.ewald,
        band_sum,
        double_counting,
        band_form,
        harris_foulkes,
        functional_difference_bound: rh.abs() + rx.abs() + eig_bound,
    };
    if [
        audit.total_hartree(),
        audit.band_sum,
        audit.double_counting,
        audit.band_form,
        audit.harris_foulkes,
        audit.functional_difference_bound,
    ]
    .iter()
    .any(|v| !v.is_finite())
    {
        return Err(ScfError::NonFinite);
    }
    Ok(audit)
}

// Each audit uses the SAME total density and input potential. Thus averaging
// adds k weights only to orbital terms; common nonlinear terms occur once.
pub(crate) fn weighted_audits(audits: &[(f64, EnergyAudit)]) -> EnergyAudit {
    macro_rules! sum {
        ($field:ident) => {
            audits.iter().map(|(w, a)| w * a.$field).sum()
        };
    }
    EnergyAudit {
        xc_floor: None,
        nlcc: audits.first().and_then(|(_, a)| a.nlcc.clone()),
        kinetic: sum!(kinetic),
        local: sum!(local),
        local_zero_mode: sum!(local_zero_mode),
        nonlocal: sum!(nonlocal),
        hartree: sum!(hartree),
        xc: sum!(xc),
        ewald: sum!(ewald),
        band_sum: sum!(band_sum),
        double_counting: sum!(double_counting),
        band_form: sum!(band_form),
        harris_foulkes: sum!(harris_foulkes),
        functional_difference_bound: sum!(functional_difference_bound),
    }
}
