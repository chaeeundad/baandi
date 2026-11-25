//! PBE on owned full-box Fourier coefficients and real slabs. No gathers.
use super::DENSITY_THRESHOLD;
use crate::PhysicsError;
use dft_geometry::Cell;
use dft_mpi::Collective;
use ndarray::Array1;
use num_complex::Complex64 as C;
use xcx::{Functional, FunctionalId, Spin, XcInput};

/// Existing run-owned FFT plan, with normalized forward / unscaled inverse.
/// Implementations must participate collectively even with empty ownership.
pub trait Fields<K: Collective> {
    type Error: From<PhysicsError>;
    fn dims(&self) -> [usize; 3];
    fn reciprocal_indices(&self) -> &[usize];
    fn plane_points(&self) -> usize;
    /// Global C-order indices of the owned real samples, in complete z rows.
    /// Required by opt-in screening and precise invalid-density diagnostics.
    fn real_indices(&self) -> Option<Vec<usize>> {
        None
    }
    fn forward(&mut self, comm: &K, n: &[f64]) -> Result<Vec<C>, Self::Error>;
    fn inverse(&mut self, comm: &K, n: &[C]) -> Result<Vec<f64>, Self::Error>;
}

fn agree<K: Collective>(comm: &K, bad: bool, message: &'static str) -> Result<(), PhysicsError> {
    let mut result = [0.];
    comm.tree_sum_f64(&[f64::from(bad)], &mut result)
        .map_err(|e| PhysicsError::Xc(e.to_string()))?;
    if result[0] != 0. {
        return Err(PhysicsError::Xc(message.into()));
    }
    Ok(())
}

/// All ranks enter this boundary before any transform, including empty owners.
pub fn audit<K: Collective, F: Fields<K>>(
    fields: &F,
    comm: &K,
    density: &[f64],
    weight: f64,
    channel: &'static str,
    electrons: Option<f64>,
) -> Result<(), PhysicsError> {
    if crate::xc::screening_prechecked() {
        return Ok(());
    }
    let enabled = crate::xc::screen_negative();
    let mut bad = [0.];
    comm.tree_sum_f64(
        &[f64::from(density.iter().any(|v| !v.is_finite() || *v < 0.))],
        &mut bad,
    )
    .map_err(|e| PhysicsError::Xc(e.to_string()))?;
    if !enabled && bad[0] == 0. {
        return Ok(());
    }
    let indices = fields.real_indices();
    agree(
        comm,
        indices.as_ref().is_none_or(|v| v.len() != density.len()),
        "XC screening requires real-grid indices",
    )?;
    let indices = indices.unwrap();
    if !enabled {
        let mut local = [f64::INFINITY, 0.];
        for (&n, &i) in density.iter().zip(&indices) {
            if !n.is_finite() || n < local[0] {
                local = [n, i as f64];
            }
        }
        let mut minimum = [f64::INFINITY, 0.];
        for rank in 0..comm.size() {
            let mut values = [0.; 2];
            comm.tree_sum_f64(
                &if rank == comm.rank() { local } else { [0.; 2] },
                &mut values,
            )
            .map_err(|e| PhysicsError::Xc(e.to_string()))?;
            if values[0].is_nan()
                || values[0] < minimum[0]
                || (values[0] == minimum[0] && values[1] < minimum[1])
            {
                minimum = values;
            }
        }
        return Err(PhysicsError::Xc(format!(
            "{}: XC requires finite nonnegative density; channel {} minimum {:.17e} Bohr^-3 at grid point {} (screening off)",
            crate::xc::screening_stage(),
            channel,
            minimum[0],
            minimum[1] as usize
        )));
    }
    let rows = crate::xc::screening_rows(density, &indices, fields.dims());
    let mut global = vec![0.; rows.len()];
    comm.tree_sum_f64(&rows, &mut global)
        .map_err(|e| PhysicsError::Xc(e.to_string()))?;
    let r = crate::xc::finish_screening(&global, weight, channel, electrons)?;
    if !enabled && r.minimum_bohr3 < 0. {
        return Err(PhysicsError::Xc(format!(
            "{}: XC channel {} minimum {:.17e} Bohr^-3 at grid point {} (screening off)",
            r.stage, channel, r.minimum_bohr3, r.minimum_grid_index
        )));
    }
    Ok(())
}

/// Cartesian multiplier of the same real, skew-adjoint derivative as serial
/// GGA: remove each even-axis Nyquist Miller component BEFORE applying B.
fn multiplier(dims: [usize; 3], index: usize, b: [[f64; 3]; 3], axis: usize) -> C {
    let ijk = [
        index / (dims[1] * dims[2]),
        index / dims[2] % dims[1],
        index % dims[2],
    ];
    C::new(
        0.,
        (0..3)
            .map(|a| {
                let m = if dims[a].is_multiple_of(2) && ijk[a] == dims[a] / 2 {
                    0.
                } else if ijk[a] <= (dims[a] - 1) / 2 {
                    ijk[a] as f64
                } else {
                    ijk[a] as f64 - dims[a] as f64
                };
                b[axis][a] * m
            })
            .sum(),
    )
}

/// Three Cartesian gradients (Bohr^-4), reusing the caller's rho(G).
pub fn gradient<K: Collective, F: Fields<K>>(
    fields: &mut F,
    comm: &K,
    cell: &Cell,
    density_g: &[C],
) -> Result<[Vec<f64>; 3], F::Error> {
    agree(
        comm,
        density_g.len() != fields.reciprocal_indices().len()
            || density_g
                .iter()
                .any(|z| !z.re.is_finite() || !z.im.is_finite()),
        "invalid distributed PBE coefficients on a rank",
    )?;
    let b = cell.reciprocal_columns_inverse_bohr();
    let mut out = std::array::from_fn(|_| Vec::new());
    let mut work = vec![C::default(); density_g.len()];
    for (axis, values) in out.iter_mut().enumerate() {
        for ((v, &n), &g) in work
            .iter_mut()
            .zip(density_g)
            .zip(fields.reciprocal_indices())
        {
            *v = n * multiplier(fields.dims(), g, b, axis);
        }
        *values = fields.inverse(comm, &work)?;
    }
    Ok(out)
}

/// Owned potential (Ha) and LOCAL energy integral (Ha). The caller reduces the
/// energy together with its other field integrals. Density must be nonnegative,
/// exactly as in serial GGA; no absolute-value or clipping policy is added.
pub struct XcSlab {
    pub potential: Vec<f64>,
    pub energy: f64,
}

pub fn evaluate<K: Collective, F: Fields<K>>(
    fields: &mut F,
    comm: &K,
    cell: &Cell,
    density: &[f64],
    density_g: &[C],
) -> Result<XcSlab, F::Error> {
    agree(
        comm,
        density.len() != fields.plane_points(),
        "distributed PBE requires finite nonnegative density with matching slab shape",
    )?;
    audit(
        fields,
        comm,
        density,
        cell.volume_bohr3() / fields.dims().iter().product::<usize>() as f64,
        "total",
        None,
    )?;
    let screening = crate::xc::screen_negative();
    let mut screened_points = [0.];
    if screening {
        comm.tree_sum_f64(
            &[density.iter().filter(|&&n| n <= DENSITY_THRESHOLD).count() as f64],
            &mut screened_points,
        )
        .map_err(|e| PhysicsError::Xc(e.to_string()))?;
    }
    let projected_g = if screening && screened_points[0] > 0. {
        Some(
            fields.forward(
                comm,
                &density
                    .iter()
                    .map(|&n| crate::xc::screened_argument(n))
                    .collect::<Vec<_>>(),
            )?,
        )
    } else {
        None
    };
    let gradients = gradient(
        fields,
        comm,
        cell,
        projected_g.as_deref().unwrap_or(density_g),
    )?;
    let rho: Vec<_> = density
        .iter()
        .map(|&n| if n <= DENSITY_THRESHOLD { 0. } else { n })
        .collect();
    let mut sigma = vec![0.; density.len()];
    for grad in &gradients {
        for (s, g) in sigma.iter_mut().zip(grad) {
            *s += g * g;
        }
    }
    for (s, n) in sigma.iter_mut().zip(&rho) {
        if *n == 0. {
            *s = 0.;
        }
    }
    // Empty slab owners skip only the local library call, never transforms.
    let local = (|| -> Result<_, PhysicsError> {
        if rho.is_empty() {
            return Ok((Vec::new(), Vec::new(), 0.));
        }
        let input = XcInput::gga(&rho, &sigma);
        let eval = |id| {
            Functional::new(id, Spin::Unpolarized)
                .and_then(|f| f.eval(rho.len(), &input))
                .map_err(|e| PhysicsError::Xc(e.to_string()))
        };
        let x = eval(FunctionalId::GgaXPbe)?;
        let c = eval(FunctionalId::GgaCPbe)?;
        let energy = cell.volume_bohr3() / fields.dims().iter().product::<usize>() as f64
            * Array1::from_iter(
                rho.iter()
                    .zip(&x.exc)
                    .zip(&c.exc)
                    .map(|((n, x), c)| n * (x + c)),
            )
            .sum();
        let potential: Vec<_> = x.vrho.iter().zip(&c.vrho).map(|(x, c)| x + c).collect();
        let flux_factor: Vec<_> = x
            .vsigma
            .iter()
            .zip(&c.vsigma)
            .map(|(x, c)| 2. * (x + c))
            .collect();
        if !energy.is_finite() || potential.iter().chain(&flux_factor).any(|v| !v.is_finite()) {
            return Err(PhysicsError::Unrepresentable);
        }
        Ok((potential, flux_factor, energy))
    })();
    agree(
        comm,
        local.is_err(),
        "distributed PBE library evaluation failed on a rank",
    )?;
    let (mut potential, factor, energy) = local?;
    let b = cell.reciprocal_columns_inverse_bohr();
    let mut divergence_g = vec![C::default(); density_g.len()];
    for (axis, mut flux) in gradients.into_iter().enumerate() {
        for (f, a) in flux.iter_mut().zip(&factor) {
            *f *= a;
        }
        let coeff = fields.forward(comm, &flux)?;
        for ((sum, z), &g) in divergence_g
            .iter_mut()
            .zip(coeff)
            .zip(fields.reciprocal_indices())
        {
            *sum += z * multiplier(fields.dims(), g, b, axis);
        }
    }
    // Three flux forwards share one inverse of the accumulated divergence.
    let divergence = fields.inverse(comm, &divergence_g)?;
    for (v, d) in potential.iter_mut().zip(divergence) {
        *v -= d;
    }
    if screening {
        for (v, &n) in potential.iter_mut().zip(density) {
            if n <= DENSITY_THRESHOLD {
                *v = 0.;
            }
        }
    }
    agree(
        comm,
        potential.iter().any(|v| !v.is_finite()),
        "nonfinite distributed PBE potential on a rank",
    )?;
    Ok(XcSlab { potential, energy })
}

/// Owned polarized PBE response. Integrals are LOCAL (Ha); the caller reduces
/// them with its energy ledger. The core integral is diagnostic, not band DC.
pub struct SpinXcSlab {
    pub potentials: [Vec<f64>; 2],
    pub energy: f64,
    pub valence_potential_integral: f64,
    pub core_potential_integral: f64,
}

/// Run-owned polarized PBE adapter. All adapter scratch and response storage
/// are retained across maps. The existing Fields transform API and xcx library
/// return owned vectors; their boundary allocations remain, outside H/orbitals.
pub struct SpinPbe {
    exchange: Functional,
    correlation: Functional,
    dims: [usize; 3],
    indices: Vec<usize>,
    multipliers: [Vec<C>; 3],
    rho: Vec<f64>,
    sigma: Vec<f64>,
    channel: Vec<f64>,
    gradients: [[Vec<f64>; 3]; 2],
    flux: Vec<f64>,
    derivative_g: Vec<C>,
    weight: f64,
    response: SpinXcSlab,
}
impl SpinPbe {
    /// Bind the workspace to one cell and owned FFT layout for the entire run.
    pub fn new<K: Collective, F: Fields<K>>(fields: &F, cell: &Cell) -> Result<Self, PhysicsError> {
        let dims = fields.dims();
        if dims.contains(&0) {
            return Err(PhysicsError::Shape);
        }
        let len = fields.plane_points();
        let indices = fields.reciprocal_indices().to_vec();
        let b = cell.reciprocal_columns_inverse_bohr();
        Ok(Self {
            exchange: Functional::new(FunctionalId::GgaXPbe, Spin::Polarized)
                .map_err(|e| PhysicsError::Xc(e.to_string()))?,
            correlation: Functional::new(FunctionalId::GgaCPbe, Spin::Polarized)
                .map_err(|e| PhysicsError::Xc(e.to_string()))?,
            multipliers: std::array::from_fn(|axis| {
                indices
                    .iter()
                    .map(|&g| multiplier(dims, g, b, axis))
                    .collect()
            }),
            derivative_g: vec![C::default(); indices.len()],
            dims,
            indices,
            rho: vec![0.; 2 * len],
            sigma: vec![0.; 3 * len],
            channel: vec![0.; len],
            gradients: std::array::from_fn(|_| std::array::from_fn(|_| vec![0.; len])),
            flux: vec![0.; len],
            weight: cell.volume_bohr3() / dims.iter().product::<usize>() as f64,
            response: SpinXcSlab {
                potentials: std::array::from_fn(|_| vec![0.; len]),
                energy: 0.,
                valence_potential_integral: 0.,
                core_potential_integral: 0.,
            },
        })
    }

    /// Differentiate the unscreened combined channels n_s+c/2, as serial SpinXc
    /// does. A signed core is permitted, but each combined channel must be finite
    /// and nonnegative. Only total density <= the serial spin threshold screens
    /// local library inputs; no clipping or per-channel vacuum floor is applied.
    pub fn evaluate<K: Collective, F: Fields<K>>(
        &mut self,
        fields: &mut F,
        comm: &K,
        density: [&[f64]; 2],
        core: Option<&[f64]>,
    ) -> Result<&SpinXcSlab, F::Error> {
        let len = self.channel.len();
        agree(
            comm,
            fields.dims() != self.dims
                || fields.reciprocal_indices() != self.indices
                || fields.plane_points() != len
                || density.iter().any(|n| n.len() != len)
                || core.is_some_and(|c| c.len() != len || c.iter().any(|v| !v.is_finite())),
            "distributed spin PBE field/core shape mismatch",
        )?;
        let screening = crate::xc::screen_negative();
        let electrons = if screening && !crate::xc::screening_prechecked() {
            let indices = fields.real_indices();
            agree(
                comm,
                indices.as_ref().is_none_or(|v| v.len() != len),
                "XC screening requires real-grid indices",
            )?;
            for (i, n) in self.channel.iter_mut().enumerate() {
                *n = density[0][i] + density[1][i];
            }
            let rows = crate::xc::screening_rows(&self.channel, &indices.unwrap(), self.dims);
            let local: Vec<_> = rows.as_chunks::<5>().0.iter().map(|r| r[2]).collect();
            let mut global = vec![0.; local.len()];
            comm.tree_sum_f64(&local, &mut global)
                .map_err(|e| PhysicsError::Xc(e.to_string()))?;
            Some(self.weight * global.iter().sum::<f64>())
        } else {
            None
        };
        for (s, n) in density.iter().enumerate() {
            for (i, &n) in n.iter().enumerate() {
                self.channel[i] = n + core.map_or(0., |c| 0.5 * c[i]);
            }
            audit(
                fields,
                comm,
                &self.channel,
                self.weight,
                if s == 0 { "up" } else { "down" },
                electrons,
            )?;
            for (i, &n) in self.channel.iter().enumerate() {
                self.rho[2 * i + s] = if screening {
                    crate::xc::screened_argument(n)
                } else {
                    n
                };
            }
        }
        for s in 0..2 {
            for (i, v) in self.channel.iter_mut().enumerate() {
                *v = self.rho[2 * i + s];
            }
            // Forward the actual combined samples, matching serial rounding and
            // avoiding any cached-core/valence Fourier cancellation at vacuum.
            let coeff = fields.forward(comm, &self.channel)?;
            for axis in 0..3 {
                for ((d, n), g) in self
                    .derivative_g
                    .iter_mut()
                    .zip(&coeff)
                    .zip(&self.multipliers[axis])
                {
                    *d = n * g;
                }
                self.gradients[s][axis].copy_from_slice(&fields.inverse(comm, &self.derivative_g)?);
            }
        }
        self.sigma.fill(0.);
        for axis in 0..3 {
            for i in 0..len {
                let u = self.gradients[0][axis][i];
                let d = self.gradients[1][axis][i];
                self.sigma[3 * i] += u * u;
                self.sigma[3 * i + 1] += u * d;
                self.sigma[3 * i + 2] += d * d;
            }
        }
        for i in 0..len {
            if self.rho[2 * i] + self.rho[2 * i + 1] <= crate::xc::DENSITY_THRESHOLD {
                self.rho[2 * i..2 * i + 2].fill(0.);
                self.sigma[3 * i..3 * i + 3].fill(0.);
            }
        }
        let local = (|| -> Result<_, PhysicsError> {
            if len == 0 {
                return Ok(None);
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
            if x.exc
                .iter()
                .chain(&c.exc)
                .chain(&x.vrho)
                .chain(&c.vrho)
                .chain(&x.vsigma)
                .chain(&c.vsigma)
                .any(|v| !v.is_finite())
            {
                return Err(PhysicsError::Unrepresentable);
            }
            Ok(Some((x, c)))
        })();
        agree(
            comm,
            local.is_err(),
            "distributed spin PBE library evaluation failed on a rank",
        )?;
        let local = local?;
        self.response.energy = local.as_ref().map_or(0., |(x, c)| {
            self.weight
                * (0..len)
                    .map(|i| (self.rho[2 * i] + self.rho[2 * i + 1]) * (x.exc[i] + c.exc[i]))
                    .sum::<f64>()
        });
        for s in 0..2 {
            if let Some((x, c)) = &local {
                for i in 0..len {
                    self.response.potentials[s][i] = x.vrho[2 * i + s] + c.vrho[2 * i + s];
                }
            }
            for axis in 0..3 {
                if let Some((x, c)) = &local {
                    for i in 0..len {
                        self.flux[i] = 2.
                            * (x.vsigma[3 * i + 2 * s] + c.vsigma[3 * i + 2 * s])
                            * self.gradients[s][axis][i]
                            + (x.vsigma[3 * i + 1] + c.vsigma[3 * i + 1])
                                * self.gradients[1 - s][axis][i];
                    }
                }
                let coeff = fields.forward(comm, &self.flux)?;
                for ((d, n), g) in self
                    .derivative_g
                    .iter_mut()
                    .zip(&coeff)
                    .zip(&self.multipliers[axis])
                {
                    *d = n * g;
                }
                // Match serial's axis-by-axis subtraction, including empty owners.
                let derivative = fields.inverse(comm, &self.derivative_g)?;
                for (v, d) in self.response.potentials[s].iter_mut().zip(derivative) {
                    *v -= d;
                }
            }
        }
        if screening {
            for s in 0..2 {
                for (i, v) in self.response.potentials[s].iter_mut().enumerate() {
                    if self.rho[2 * i + s] == 0. {
                        *v = 0.;
                    }
                }
            }
        }
        self.response.valence_potential_integral = self.weight
            * (0..2)
                .map(|s| {
                    density[s]
                        .iter()
                        .zip(&self.response.potentials[s])
                        .map(|(n, v)| n * v)
                        .sum::<f64>()
                })
                .sum::<f64>();
        self.response.core_potential_integral = core.map_or(0., |core| {
            self.weight
                * core
                    .iter()
                    .zip(&self.response.potentials[0])
                    .zip(&self.response.potentials[1])
                    .map(|((&c, &u), &d)| 0.5 * c * (u + d))
                    .sum::<f64>()
        });
        agree(
            comm,
            !self.response.energy.is_finite()
                || !self.response.valence_potential_integral.is_finite()
                || !self.response.core_potential_integral.is_finite()
                || self
                    .response
                    .potentials
                    .iter()
                    .flatten()
                    .any(|v| !v.is_finite()),
            "nonfinite distributed spin PBE response on a rank",
        )?;
        Ok(&self.response)
    }
}
