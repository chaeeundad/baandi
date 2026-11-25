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
        density.len() != fields.plane_points() || density.iter().any(|n| !n.is_finite() || *n < 0.),
        "distributed PBE requires finite nonnegative density with matching slab shape",
    )?;
    let gradients = gradient(fields, comm, cell, density_g)?;
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
    agree(
        comm,
        potential.iter().any(|v| !v.is_finite()),
        "nonfinite distributed PBE potential on a rank",
    )?;
    Ok(XcSlab { potential, energy })
}
