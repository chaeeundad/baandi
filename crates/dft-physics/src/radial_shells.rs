//! GD1: per-|G| shell values of the local and core radial transforms used by
//! the post-SCF derivatives (local/core forces, local/core stress and the
//! core density they rebuild).
//!
//! [`ShellEvaluator`] evaluates each value with the G-independent part of the
//! integrand (including the `erf` of the local integrand) formed once per
//! potential instead of once per shell; every value equals the per-shell
//! expression of the former code bit for bit (same operands, same order of
//! operations, same index-coordinate quadrature). [`RadialShellValues`] holds
//! values computed elsewhere (for example split over MPI ranks) for lookup by
//! the derivative kernels, which compute a missing value with their own
//! former expression.
use crate::hamiltonian::HamiltonianError;
use dft_upf::{
    Pseudopotential,
    radial::{RadialError, RadialMesh, spherical_bessel},
};
use std::collections::HashMap;
use std::f64::consts::PI;

/// Which radial function of a potential a shell value belongs to.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum ShellKind {
    /// `LocalTransform::at_g(q)` in Ha (q > 0).
    LocalValue,
    /// d/dq of the local transform, the frozen-Miller local stress factor
    /// (`4 pi / V` included; q > 0).
    LocalSlope,
    /// `spherical_density_transform(mesh, core, q)` in electrons (no volume).
    CoreValue,
    /// d/dq of the core transform divided by the volume (core stress factor).
    CoreSlope,
}

/// Identity of a potential inside one derivative evaluation.
pub fn potential_id(pp: &Pseudopotential) -> usize {
    pp as *const _ as *const () as usize
}

/// Shell values keyed by (potential, kind, bits of q).
#[derive(Debug, Default, Clone)]
pub struct RadialShellValues {
    map: HashMap<(usize, ShellKind, u64), f64>,
}
impl RadialShellValues {
    pub fn insert(&mut self, pp: &Pseudopotential, kind: ShellKind, q: f64, value: f64) {
        self.map
            .insert((potential_id(pp), kind, q.to_bits()), value);
    }
    #[inline]
    pub fn get(&self, pp: &Pseudopotential, kind: ShellKind, q: f64) -> Option<f64> {
        self.map
            .get(&(potential_id(pp), kind, q.to_bits()))
            .copied()
    }
    pub fn len(&self) -> usize {
        self.map.len()
    }
    pub fn is_empty(&self) -> bool {
        self.map.is_empty()
    }
}

/// Hoisted integrands of one potential at one cell volume.
pub struct ShellEvaluator<'a> {
    pp: &'a Pseudopotential,
    volume: f64,
    /// Local prefix mesh (the endpoint of `local_transform(volume, None)`).
    prefix: RadialMesh,
    /// r^2 V(r) + Z r erf(r) on the prefix.
    local_base: Vec<f64>,
    /// -(r^2 V(r) + Z r erf(r)) r on the prefix.
    local_slope_base: Vec<f64>,
    /// r^2 rho_core(r) on the full mesh.
    core_base: Option<Vec<f64>>,
    /// -r^3 rho_core(r) on the full mesh.
    core_slope_base: Option<Vec<f64>>,
    scratch: Vec<f64>,
}
impl<'a> ShellEvaluator<'a> {
    pub fn new(pp: &'a Pseudopotential, volume: f64) -> Result<Self, HamiltonianError> {
        let local = pp.local_transform(volume, None)?;
        let mesh = pp.mesh();
        let n = mesh.r().partition_point(|r| *r <= local.effective_r_max());
        let prefix = RadialMesh::new(mesh.r()[..n].to_vec(), mesh.rab()[..n].to_vec())?;
        let z = pp.z_valence();
        let local_base: Vec<f64> = (0..n)
            .map(|i| {
                let r = mesh.r()[i];
                r * r * pp.local_hartree()[i] + z * r * libm::erf(r)
            })
            .collect();
        let local_slope_base = (0..n).map(|i| -local_base[i] * mesh.r()[i]).collect();
        let core_base = pp
            .core_density()
            .map(|c| mesh.r().iter().zip(c).map(|(&r, &n)| r * r * n).collect());
        let core_slope_base = pp.core_density().map(|c| {
            mesh.r()
                .iter()
                .zip(c)
                .map(|(&r, &n)| -r * r * r * n)
                .collect()
        });
        Ok(Self {
            pp,
            volume,
            prefix,
            local_base,
            local_slope_base,
            core_base,
            core_slope_base,
            scratch: Vec::with_capacity(mesh.r().len()),
        })
    }
    /// The value of `kind` at `q` (> 0 for the local kinds); `None` for a core
    /// kind of a potential without a core density. Equal bit for bit to
    /// `finish(kind, q, integral(kind, q))`.
    pub fn eval(&mut self, kind: ShellKind, q: f64) -> Result<Option<f64>, HamiltonianError> {
        if kind == ShellKind::LocalValue && (!q.is_finite() || q <= 0.0) {
            return Err(RadialError::Argument("G must be finite and positive").into());
        }
        match self.integral(kind, q)? {
            Some(integral) => self.finish(kind, q, integral).map(Some),
            None => Ok(None),
        }
    }
    /// GD2: the q-dependent radial quadrature of `kind` alone (the short-range
    /// local integral without the analytic Gaussian-Coulomb term, the core
    /// integral without 4 pi and volume factors), defined for every finite
    /// q >= 0; even in q for the value kinds and odd for the slope kinds
    /// (spherical Bessel j0 and j1). `None` for a core kind without a core.
    pub fn integral(&mut self, kind: ShellKind, q: f64) -> Result<Option<f64>, HamiltonianError> {
        let mesh = self.pp.mesh();
        let fill = |scratch: &mut Vec<f64>, base: &[f64], l: usize, r: &[f64]| {
            scratch.clear();
            for (b, &r) in base.iter().zip(r) {
                scratch.push(b * spherical_bessel(l, q * r)?);
            }
            Ok::<(), RadialError>(())
        };
        Ok(Some(match kind {
            ShellKind::LocalValue => {
                fill(&mut self.scratch, &self.local_base, 0, self.prefix.r())?;
                self.prefix.integrate(&self.scratch)?
            }
            ShellKind::LocalSlope => {
                fill(
                    &mut self.scratch,
                    &self.local_slope_base,
                    1,
                    self.prefix.r(),
                )?;
                self.prefix.integrate(&self.scratch)?
            }
            ShellKind::CoreValue => {
                let Some(base) = &self.core_base else {
                    return Ok(None);
                };
                fill(&mut self.scratch, base, 0, mesh.r())?;
                mesh.integrate(&self.scratch)?
            }
            ShellKind::CoreSlope => {
                let Some(base) = &self.core_slope_base else {
                    return Ok(None);
                };
                fill(&mut self.scratch, base, 1, mesh.r())?;
                mesh.integrate(&self.scratch)?
            }
        }))
    }
    /// GD2: the shell value of `kind` at `q` from its radial `integral`
    /// (exact or tabulated): the same operations as the former per-shell
    /// expression after the quadrature. The local kinds need q > 0.
    pub fn finish(&self, kind: ShellKind, q: f64, integral: f64) -> Result<f64, HamiltonianError> {
        let z = self.pp.z_valence();
        Ok(match kind {
            ShellKind::LocalValue => {
                let v = 4.0 * PI / self.volume * (integral - z * (-q * q / 4.0).exp() / (q * q));
                if !v.is_finite() {
                    return Err(RadialError::Argument("unrepresentable result").into());
                }
                v
            }
            ShellKind::LocalSlope => {
                4. * PI / self.volume
                    * (integral + z * (-q * q / 4.).exp() * (0.5 / q + 2. / q.powi(3)))
            }
            ShellKind::CoreValue => {
                let v = 4.0 * PI * integral;
                if !v.is_finite() {
                    return Err(RadialError::Argument("unrepresentable result").into());
                }
                v
            }
            ShellKind::CoreSlope => 4. * PI / self.volume * integral,
        })
    }
}
