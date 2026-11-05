//! Collinear XC with optional spin-neutral NLCC and reusable spectral workspaces.
use crate::{PhysicsError, data::*, hartree, xc};
use dft_geometry::Cell;
use std::f64::consts::PI;
use xcx::{Functional, FunctionalId, Spin, XcInput};

/// The single authoritative density; total and magnetization are derived.
#[derive(Clone, Debug)]
pub enum SpinDensity {
    Unpolarized {
        grid_id: String,
        total: Array3<f64>,
    },
    Collinear {
        grid_id: String,
        channels: [Array3<f64>; 2],
    },
}
impl SpinDensity {
    pub fn validate(&self) -> Result<(), PhysicsError> {
        let (id, first, second) = match self {
            Self::Unpolarized { grid_id, total } => (grid_id, total, None),
            Self::Collinear { grid_id, channels } => (grid_id, &channels[0], Some(&channels[1])),
        };
        if id.is_empty()
            || id.len() > 128
            || !id.is_ascii()
            || first.is_empty()
            || !first.is_standard_layout()
            || second.is_some_and(|s| s.shape() != first.shape() || !s.is_standard_layout())
        {
            return Err(PhysicsError::Shape);
        }
        for &n in first
            .iter()
            .chain(second.into_iter().flat_map(|s| s.iter()))
        {
            xc::checked_density(n)?;
        }
        Ok(())
    }
    pub fn total_into(&self, out: &mut Array3<f64>) -> Result<(), PhysicsError> {
        self.validate()?;
        let first = match self {
            Self::Unpolarized { total, .. } => total,
            Self::Collinear { channels, .. } => &channels[0],
        };
        if out.shape() != first.shape() {
            return Err(PhysicsError::Shape);
        }
        out.assign(first);
        if let Self::Collinear { channels, .. } = self {
            *out += &channels[1];
        }
        Ok(())
    }
    /// Signed electron difference/Bohr^3 along the result's recorded axis.
    pub fn magnetization_into(&self, out: &mut Array3<f64>) -> Result<(), PhysicsError> {
        self.validate()?;
        match self {
            Self::Unpolarized { total, .. } => {
                if total.shape() != out.shape() {
                    return Err(PhysicsError::Shape);
                }
                out.fill(0.0);
            }
            Self::Collinear { channels, .. } => {
                if channels[0].shape() != out.shape() {
                    return Err(PhysicsError::Shape);
                }
                out.assign(&channels[0]);
                *out -= &channels[1];
            }
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpinFunctional {
    LdaPw,
    Pbe,
}

/// Ionic derivatives are unavailable even when valence XC and NLCC are evaluated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum IonicDerivativeCapability {
    Unavailable { reason: &'static str },
}

pub struct SpinXcResponse {
    pub potentials: [Array3<f64>; 2],
    /// First density functional derivatives are the supplied potentials.
    /// They do not establish ionic force or stress capabilities.
    pub force_capability: IonicDerivativeCapability,
    pub stress_capability: IonicDerivativeCapability,
    pub energy_ha: f64,
    pub valence_potential_integral_ha: f64,
    /// Audit only: integral (rho_core/2)*(V_up+V_down), excluded from band double counting.
    pub core_potential_integral_ha: f64,
}
impl SpinXcResponse {
    pub fn new(dims: [usize; 3]) -> Self {
        Self {
            potentials: std::array::from_fn(|_| Array3::zeros(dims)),
            force_capability: IonicDerivativeCapability::Unavailable {
                reason: "not_implemented",
            },
            stress_capability: IonicDerivativeCapability::Unavailable {
                reason: "not_implemented",
            },
            energy_ha: 0.0,
            valence_potential_integral_ha: 0.0,
            core_potential_integral_ha: 0.0,
        }
    }
}

pub struct SpinXc {
    kind: SpinFunctional,
    exchange: Functional,
    correlation: Functional,
    rho: Vec<f64>,
    sigma: Vec<f64>,
    gradients: [[Array3<f64>; 3]; 2],
    flux: Array3<f64>,
    total_channel: Array3<f64>,
    derivative: Array3<f64>,
    spectral: Spectral,
    weight: f64,
}
impl SpinXc {
    pub fn new(kind: SpinFunctional, cell: &Cell, dims: [usize; 3]) -> Result<Self, PhysicsError> {
        let (x, c) = match kind {
            SpinFunctional::LdaPw => (FunctionalId::LdaX, FunctionalId::LdaCPw),
            SpinFunctional::Pbe => (FunctionalId::GgaXPbe, FunctionalId::GgaCPbe),
        };
        let spectral = Spectral::new(cell, dims)?;
        let len = spectral.real.len();
        Ok(Self {
            kind,
            exchange: Functional::new(x, Spin::Polarized)
                .map_err(|e| PhysicsError::Xc(e.to_string()))?,
            correlation: Functional::new(c, Spin::Polarized)
                .map_err(|e| PhysicsError::Xc(e.to_string()))?,
            rho: vec![0.0; 2 * len],
            sigma: vec![0.0; 3 * len],
            gradients: std::array::from_fn(|_| std::array::from_fn(|_| Array3::zeros(dims))),
            flux: Array3::zeros(dims),
            total_channel: Array3::zeros(dims),
            derivative: Array3::zeros(dims),
            spectral,
            weight: cell.volume_bohr3() / len as f64,
        })
    }
    /// xcx owns its result vectors (explicit architect-approved allocation boundary).
    /// All adapter buffers are retained; no evaluation occurs per orbital/k point.
    pub fn evaluate(
        &mut self,
        density: &[Array3<f64>; 2],
        out: &mut SpinXcResponse,
    ) -> Result<(), PhysicsError> {
        self.evaluate_with_core(density, None, out)
    }
    /// Immutable, spin-neutral model core. Apply the existing strict collinear
    /// nonnegative-density rule to each combined channel, without clipping.
    pub fn evaluate_with_core(
        &mut self,
        density: &[Array3<f64>; 2],
        core: Option<&Array3<f64>>,
        out: &mut SpinXcResponse,
    ) -> Result<(), PhysicsError> {
        if core.is_some_and(|c| c.shape() != self.flux.shape() || c.iter().any(|v| !v.is_finite()))
        {
            return Err(PhysicsError::Shape);
        }
        let len = self.rho.len() / 2;
        for (s, channel) in density.iter().enumerate() {
            if channel.shape() != self.flux.shape()
                || out.potentials[s].shape() != self.flux.shape()
            {
                return Err(PhysicsError::Shape);
            }
            if let Some(core) = core {
                for ((total, &valence), &c) in self.total_channel.iter_mut().zip(channel).zip(core)
                {
                    *total = valence + 0.5 * c;
                }
            }
            let combined = if core.is_some() {
                &self.total_channel
            } else {
                channel
            };
            for (i, &n) in combined.iter().enumerate() {
                if !n.is_finite() || n < 0.0 {
                    return Err(if core.is_some() {
                        PhysicsError::Xc(format!(
                            "NLCC combined spin density {n:e} Bohr^-3 is not finite/nonnegative; resolve density-grid truncation or valence feasibility without clipping core"
                        ))
                    } else {
                        PhysicsError::Density(n)
                    });
                }
                self.rho[2 * i + s] = n;
            }
            if self.kind == SpinFunctional::Pbe {
                for a in 0..3 {
                    self.spectral
                        .derivative(combined, a, &mut self.gradients[s][a])?;
                }
            }
        }
        self.sigma.fill(0.0);
        if self.kind == SpinFunctional::Pbe {
            for a in 0..3 {
                for (i, (&u, &d)) in self.gradients[0][a]
                    .iter()
                    .zip(&self.gradients[1][a])
                    .enumerate()
                {
                    self.sigma[3 * i] += u * u;
                    self.sigma[3 * i + 1] += u * d;
                    self.sigma[3 * i + 2] += d * d;
                }
            }
        }
        // Screen only the total density, preserving fully polarized nonvacuum points.
        for i in 0..len {
            if self.rho[2 * i] + self.rho[2 * i + 1] <= xc::DENSITY_THRESHOLD {
                self.rho[2 * i] = 0.0;
                self.rho[2 * i + 1] = 0.0;
                self.sigma[3 * i..3 * i + 3].fill(0.0);
            }
        }
        let input = XcInput::gga(&self.rho, &self.sigma);
        let x = self
            .exchange
            .eval(len, &input)
            .map_err(|e| PhysicsError::Xc(e.to_string()))?;
        let c = self
            .correlation
            .eval(len, &input)
            .map_err(|e| PhysicsError::Xc(e.to_string()))?;
        out.energy_ha = self.weight
            * (0..len)
                .map(|i| (self.rho[2 * i] + self.rho[2 * i + 1]) * (x.exc[i] + c.exc[i]))
                .sum::<f64>();
        for s in 0..2 {
            for (i, v) in out.potentials[s].iter_mut().enumerate() {
                *v = x.vrho[2 * i + s] + c.vrho[2 * i + s];
            }
            if self.kind == SpinFunctional::Pbe {
                for a in 0..3 {
                    for (i, ((f, &same), &other)) in self
                        .flux
                        .iter_mut()
                        .zip(&self.gradients[s][a])
                        .zip(&self.gradients[1 - s][a])
                        .enumerate()
                    {
                        *f = 2.0 * (x.vsigma[3 * i + 2 * s] + c.vsigma[3 * i + 2 * s]) * same
                            + (x.vsigma[3 * i + 1] + c.vsigma[3 * i + 1]) * other;
                    }
                    self.spectral
                        .derivative(&self.flux, a, &mut self.derivative)?;
                    out.potentials[s] -= &self.derivative;
                }
            }
        }
        out.valence_potential_integral_ha = self.weight
            * (0..2)
                .map(|s| {
                    density[s]
                        .iter()
                        .zip(&out.potentials[s])
                        .map(|(n, v)| n * v)
                        .sum::<f64>()
                })
                .sum::<f64>();
        out.core_potential_integral_ha = core.map_or(0.0, |core| {
            self.weight
                * core
                    .iter()
                    .zip(&out.potentials[0])
                    .zip(&out.potentials[1])
                    .map(|((&c, &u), &d)| 0.5 * c * (u + d))
                    .sum::<f64>()
        });
        if !out.core_potential_integral_ha.is_finite()
            || !out.energy_ha.is_finite()
            || !out.valence_potential_integral_ha.is_finite()
            || out.potentials.iter().flatten().any(|v| !v.is_finite())
        {
            return Err(PhysicsError::Unrepresentable);
        }
        Ok(())
    }
}

/// Reused FFT storage; Cartesian derivatives zero even-axis Nyquist components
/// before applying B, matching ADR 0018 on skew cells.
struct Spectral {
    fft: FftPlan3d,
    real: Array3<Complex64>,
    reciprocal: Array3<Complex64>,
    g: [Array3<f64>; 3],
}
impl Spectral {
    fn new(cell: &Cell, dims: [usize; 3]) -> Result<Self, PhysicsError> {
        let fft = FftPlan3d::new(dims)?;
        let b = cell.reciprocal_columns_inverse_bohr();
        let g = std::array::from_fn(|a| {
            Array3::from_shape_fn(dims, |(i, j, k)| {
                [i, j, k]
                    .iter()
                    .enumerate()
                    .map(|(d, &index)| {
                        let m = if dims[d].is_multiple_of(2) && index == dims[d] / 2 {
                            0.0
                        } else if index <= (dims[d] - 1) / 2 {
                            index as f64
                        } else {
                            index as f64 - dims[d] as f64
                        };
                        b[a][d] * m
                    })
                    .sum()
            })
        });
        Ok(Self {
            fft,
            real: Array3::zeros(dims),
            reciprocal: Array3::zeros(dims),
            g,
        })
    }
    fn derivative(
        &mut self,
        n: &Array3<f64>,
        a: usize,
        out: &mut Array3<f64>,
    ) -> Result<(), PhysicsError> {
        for (z, &v) in self.real.iter_mut().zip(n) {
            *z = Complex64::new(v, 0.0);
        }
        self.fft.to_reciprocal(&self.real, &mut self.reciprocal)?;
        for (z, &g) in self.reciprocal.iter_mut().zip(&self.g[a]) {
            *z *= Complex64::new(0.0, g);
        }
        self.fft.to_real(&self.reciprocal, &mut self.real)?;
        for (v, z) in out.iter_mut().zip(&self.real) {
            *v = z.re;
        }
        Ok(())
    }
}

/// Allocation-free Poisson evaluation after construction; G=0 is zero.
pub struct SpinHartree {
    fft: FftPlan3d,
    real: Array3<Complex64>,
    reciprocal: Array3<Complex64>,
    g2: Array3<f64>,
    volume: f64,
    pub potential: Array3<f64>,
    pub energy_ha: f64,
}
impl SpinHartree {
    pub fn new(cell: &Cell, dims: [usize; 3]) -> Result<Self, PhysicsError> {
        Ok(Self {
            fft: FftPlan3d::new(dims)?,
            real: Array3::zeros(dims),
            reciprocal: Array3::zeros(dims),
            g2: hartree::squared_frequencies(cell, dims)?,
            volume: cell.volume_bohr3(),
            potential: Array3::zeros(dims),
            energy_ha: 0.0,
        })
    }
    pub fn evaluate(&mut self, n: &Array3<f64>) -> Result<(), PhysicsError> {
        if n.shape() != self.real.shape() {
            return Err(PhysicsError::Shape);
        }
        for (z, &v) in self.real.iter_mut().zip(n) {
            *z = Complex64::new(xc::checked_density(v)?, 0.0);
        }
        self.fft.to_reciprocal(&self.real, &mut self.reciprocal)?;
        let mut energy = 0.0;
        for ((index, z), &g) in self.reciprocal.indexed_iter_mut().zip(&self.g2) {
            if index == (0, 0, 0) {
                *z = Complex64::default();
            } else {
                let kernel = 4.0 * PI / g;
                energy += kernel * z.norm_sqr();
                *z *= kernel;
            }
        }
        self.energy_ha = energy * 0.5 * self.volume;
        self.fft.to_real(&self.reciprocal, &mut self.real)?;
        for (v, z) in self.potential.iter_mut().zip(&self.real) {
            *v = z.re;
        }
        if !self.energy_ha.is_finite() || self.potential.iter().any(|v| !v.is_finite()) {
            return Err(PhysicsError::Unrepresentable);
        }
        Ok(())
    }
}

/// Additive unpolarized PW92 path; legacy PZ/PBE arithmetic remains untouched.
pub fn unpolarized_pw(n: &Array3<f64>, volume: f64) -> Result<xc::XcGrid, PhysicsError> {
    crate::volume(volume)?;
    if n.is_empty() {
        return Err(PhysicsError::Shape);
    }
    let mut rho = Vec::with_capacity(n.len());
    for &v in n {
        let v = xc::checked_density(v)?;
        rho.push(if v <= xc::DENSITY_THRESHOLD { 0.0 } else { v });
    }
    let eval = |id| {
        Functional::new(id, Spin::Unpolarized)
            .and_then(|f| f.eval(rho.len(), &XcInput::lda(&rho)))
            .map_err(|e| PhysicsError::Xc(e.to_string()))
    };
    let x = eval(FunctionalId::LdaX)?;
    let c = eval(FunctionalId::LdaCPw)?;
    let mut potential = Array3::zeros(n.raw_dim());
    for (i, v) in potential.iter_mut().enumerate() {
        *v = x.vrho[i] + c.vrho[i];
    }
    let weight = volume / n.len() as f64;
    let energy = weight
        * rho
            .iter()
            .enumerate()
            .map(|(i, v)| v * (x.exc[i] + c.exc[i]))
            .sum::<f64>();
    let integral = weight * n.iter().zip(&potential).map(|(n, v)| n * v).sum::<f64>();
    if !energy.is_finite() || !integral.is_finite() || potential.iter().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(xc::XcGrid {
        potential,
        energy: Hartree(energy),
        potential_integral: Hartree(integral),
    })
}
