//! Explicit frozen-Miller metric derivatives. Tensile-positive Ha/Bohr^3.
//! Integer phases at fixed fractional ions/k do not differentiate under strain.
use crate::{
    PhysicsError,
    data::*,
    hamiltonian::{HamiltonianError, NonlocalOperator, Species, structure_factor},
};
use dft_basis::PlaneWaveBasis;
use dft_geometry::Cell;
use dft_upf::radial::{RadialError, RadialMesh, spherical_bessel};
use std::f64::consts::PI;
pub type Tensor = [[f64; 3]; 3];

pub fn hartree(n: &Array3<f64>, cell: &Cell, fft: &mut FftPlan3d) -> Result<Tensor, PhysicsError> {
    let mut rho = Array3::zeros(n.raw_dim());
    fft.to_reciprocal(&n.mapv(|v| Complex64::new(v, 0.)), &mut rho)?;
    let dims = fft.shape();
    let b = cell.reciprocal_columns_inverse_bohr();
    let mut out = [[0.; 3]; 3];
    for ((i, j, k), z) in rho.indexed_iter() {
        if (i, j, k) == (0, 0, 0) {
            continue;
        }
        let idx = [i, j, k];
        let m: [f64; 3] = std::array::from_fn(|a| {
            if idx[a] <= (dims[a] - 1) / 2 {
                idx[a] as f64
            } else {
                idx[a] as f64 - dims[a] as f64
            }
        });
        let q: [f64; 3] = std::array::from_fn(|a| (0..3).map(|c| b[a][c] * m[c]).sum());
        let q2: f64 = q.iter().map(|v| v * v).sum();
        for (a, row) in out.iter_mut().enumerate() {
            for (b, v) in row.iter_mut().enumerate() {
                *v += 2. * PI * z.norm_sqr() / q2
                    * (2. * q[a] * q[b] / q2 - if a == b { 1. } else { 0. });
            }
        }
    }
    if out.iter().flatten().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(out)
}
fn integrate_prefix(mesh: &RadialMesh, values: Vec<f64>) -> Result<f64, RadialError> {
    if values.len() == mesh.r().len() {
        mesh.integrate(&values)
    } else {
        let n = values.len();
        RadialMesh::new(mesh.r()[..n].to_vec(), mesh.rab()[..n].to_vec())?.integrate(&values)
    }
}
fn bessel_prime(l: usize, x: f64) -> Result<f64, RadialError> {
    if x == 0. {
        return Ok(if l == 1 { 1. / 3. } else { 0. });
    }
    if l == 0 {
        Ok(-spherical_bessel(1, x)?)
    } else {
        Ok(spherical_bessel(l - 1, x)? - (l + 1) as f64 / x * spherical_bessel(l, x)?)
    }
}
/// Local full derivative, or NLCC derivative against the mean XC potential.
/// Fields have forward 1/N Fourier normalization. Includes all zero modes.
pub fn local_or_core(
    basis: &PlaneWaveBasis,
    volume: f64,
    species: &[Species<'_>],
    field: &[Complex64],
    core: bool,
) -> Result<Tensor, HamiltonianError> {
    if field.len() != basis.vectors().len() {
        return Err(HamiltonianError::Shape);
    }
    if !volume.is_finite() || volume <= 0. {
        return Err(HamiltonianError::Volume);
    }
    let mut out = [[0.; 3]; 3];
    // The radial value and its derivative depend only on (potential, |G|):
    // atoms of one species and G vectors on one shell share them. Keyed by the
    // exact bit pattern of |G|, so a hit returns the identical f64 results.
    let mut cache: std::collections::HashMap<(usize, u64), (f64, f64)> =
        std::collections::HashMap::new();
    for s in species {
        let pp = s.potential;
        let mesh = pp.mesh();
        let local = pp.local_transform(volume, None)?;
        let count = mesh.r().partition_point(|r| *r <= local.effective_r_max());
        let potential_id = pp as *const _ as *const () as usize;
        for (g, f) in basis.vectors().iter().zip(field) {
            let q = g.norm_squared.sqrt();
            let key = (potential_id, q.to_bits());
            let (v, prime) = if let Some(hit) = cache.get(&key) {
                *hit
            } else if core {
                let Some(n) = pp.core_density() else { continue };
                let v = dft_upf::radial::spherical_density_transform(mesh, n, q)? / volume;
                let values = mesh
                    .r()
                    .iter()
                    .zip(n)
                    .map(|(&r, &n)| Ok(-r * r * r * n * spherical_bessel(1, q * r)?))
                    .collect::<Result<Vec<_>, RadialError>>()?;
                let pair = (v, 4. * PI / volume * mesh.integrate(&values)?);
                cache.insert(key, pair);
                pair
            } else if q == 0. {
                (local.alpha_z()?, 0.)
            } else {
                let values = (0..count)
                    .map(|i| {
                        let r = mesh.r()[i];
                        Ok(
                            -(r * r * pp.local_hartree()[i] + pp.z_valence() * r * libm::erf(r))
                                * r
                                * spherical_bessel(1, q * r)?,
                        )
                    })
                    .collect::<Result<Vec<_>, RadialError>>()?;
                let prime = 4. * PI / volume
                    * (integrate_prefix(mesh, values)?
                        + pp.z_valence() * (-q * q / 4.).exp() * (0.5 / q + 2. / q.powi(3)));
                let pair = (local.at_g(q)?, prime);
                cache.insert(key, pair);
                pair
            };
            let product = (f.conj() * structure_factor(g.cartesian, s.positions_bohr)?).re;
            for (a, row) in out.iter_mut().enumerate() {
                for (b, x) in row.iter_mut().enumerate() {
                    *x += product
                        * (-if a == b { v } else { 0. }
                            - if q > 0. {
                                prime * g.cartesian[a] * g.cartesian[b] / q
                            } else {
                                0.
                            });
                }
            }
        }
    }
    if out.iter().flatten().any(|v| !v.is_finite()) {
        return Err(HamiltonianError::NonFinite);
    }
    Ok(out)
}
// Regular solid harmonics via polynomial recurrence. Derivatives are propagated
// alongside values, avoiding polar-axis singularities and numerical differencing.
#[derive(Clone, Copy)]
struct Poly {
    v: Complex64,
    d: [Complex64; 3],
}
impl Poly {
    fn constant(v: f64) -> Self {
        Self {
            v: Complex64::new(v, 0.),
            d: [Complex64::default(); 3],
        }
    }
    fn scale(self, v: f64) -> Self {
        Self {
            v: self.v * v,
            d: self.d.map(|x| x * v),
        }
    }
    fn mul(self, b: Self) -> Self {
        Self {
            v: self.v * b.v,
            d: std::array::from_fn(|i| self.d[i] * b.v + self.v * b.d[i]),
        }
    }
    fn sub(self, b: Self) -> Self {
        Self {
            v: self.v - b.v,
            d: std::array::from_fn(|i| self.d[i] - b.d[i]),
        }
    }
}
fn harmonic(l: usize, m: isize, q: [f64; 3]) -> (f64, [f64; 3]) {
    let a = m.unsigned_abs();
    let one = Complex64::new(1., 0.);
    let zero = Complex64::default();
    let xy = Poly {
        v: Complex64::new(q[0], q[1]),
        d: [one, Complex64::new(0., 1.), zero],
    };
    let z = Poly {
        v: Complex64::new(q[2], 0.),
        d: [zero, zero, one],
    };
    let r2 = q.iter().map(|v| v * v).sum::<f64>();
    let rr = Poly {
        v: Complex64::new(r2, 0.),
        d: q.map(|v| Complex64::new(2. * v, 0.)),
    };
    let mut p = Poly::constant(1.);
    for k in 1..=a {
        p = p.mul(xy).scale((2 * k - 1) as f64);
    }
    if l > a {
        let mut prev = p;
        p = p.mul(z).scale((2 * a + 1) as f64);
        for k in a + 2..=l {
            let next = p
                .mul(z)
                .scale((2 * k - 1) as f64)
                .sub(prev.mul(rr).scale((k + a - 1) as f64))
                .scale(1. / (k - a) as f64);
            prev = p;
            p = next;
        }
    }
    let norm = (((2 * l + 1) as f64 / (4. * PI))
        * ((l - a + 1)..=l + a).fold(1., |v, k| v / k as f64))
    .sqrt()
        * if a == 0 { 1. } else { 2_f64.sqrt() };
    let part = |z: Complex64| if m < 0 { z.im } else { z.re };
    let scale = norm / r2.powf(l as f64 / 2.);
    (
        part(p.v) * scale,
        std::array::from_fn(|i| scale * (part(p.d[i]) - l as f64 * part(p.v) * q[i] / r2)),
    )
}
/// Full KB strain derivative, all radial channels and Hermitian D couplings.
pub fn nonlocal(
    basis: &PlaneWaveBasis,
    volume: f64,
    species: &[Species<'_>],
    orbitals: &ComplexMatrix,
    occupations: &[f64],
) -> Result<Tensor, HamiltonianError> {
    if orbitals.rows() != basis.vectors().len() || orbitals.cols() != occupations.len() {
        return Err(HamiltonianError::Shape);
    }
    if occupations.iter().any(|f| !f.is_finite() || *f < 0.) {
        return Err(HamiltonianError::NonFinite);
    }
    let operator = NonlocalOperator::new(basis, volume, species)?;
    let mut out = [[0.; 3]; 3];
    let mut atom_index = 0;
    let mut metrics: std::collections::HashMap<usize, Vec<Vec<[[Complex64; 3]; 3]>>> =
        std::collections::HashMap::new();
    for s in species {
        let pp = s.potential;
        let channels: Vec<_> = pp
            .projectors()
            .iter()
            .enumerate()
            .flat_map(|(i, p)| (-(p.l() as isize)..=p.l() as isize).map(move |m| (i, p.l(), m)))
            .collect();
        let nc = channels.len();
        let potential_id = pp as *const _ as *const () as usize;
        // `?` inside the computation prevents the entry API.
        #[allow(clippy::map_entry)]
        if !metrics.contains_key(&potential_id) {
            let mut radials: std::collections::HashMap<(u64, usize), (f64, f64)> =
                std::collections::HashMap::new();
            let mut metric = vec![vec![[[Complex64::default(); 3]; 3]; basis.vectors().len()]; nc];
            for (ig, g) in basis.momenta().iter().enumerate() {
                let q = g.norm_squared.sqrt();
                if q == 0. {
                    continue;
                }
                for (c, &(i, l, m)) in channels.iter().enumerate() {
                    let p = &pp.projectors()[i];
                    // Radial value and derivative depend only on (|G|, projector):
                    // reuse them across G vectors with the identical |G| bit pattern.
                    let (radial, prime) = if let Some(hit) = radials.get(&(q.to_bits(), i)) {
                        *hit
                    } else {
                        let radial = pp.projector_transform(i, q)?;
                        let prime = 4.
                            * PI
                            * integrate_prefix(
                                pp.mesh(),
                                p.radial_times_r()
                                    .iter()
                                    .zip(pp.mesh().r())
                                    .map(|(&u, &r)| Ok(r * r * u * bessel_prime(l, q * r)?))
                                    .collect::<Result<Vec<_>, RadialError>>()?,
                            )?;
                        radials.insert((q.to_bits(), i), (radial, prime));
                        (radial, prime)
                    };
                    let (y, dy) = harmonic(l, m, g.cartesian);
                    let phase = Complex64::new(0., -1.).powu(l as u32) / volume.sqrt();
                    for (a, row) in metric[c][ig].iter_mut().enumerate() {
                        for (b, d) in dy.iter().enumerate() {
                            row[b] = -phase
                                * g.cartesian[a]
                                * (prime * y * g.cartesian[b] / q + radial * d);
                        }
                    }
                }
            }
            metrics.insert(potential_id, metric);
        }
        let metric = &metrics[&potential_id];
        for position in s.positions_bohr {
            let atom = &operator.atoms()[atom_index];
            atom_index += 1;
            let phases = basis
                .momenta()
                .iter()
                .map(|g| structure_factor(g.cartesian, std::slice::from_ref(position)))
                .collect::<Result<Vec<_>, _>>()?;
            let mut overlaps = vec![Complex64::default(); nc];
            for (band, &f) in occupations.iter().enumerate() {
                if f == 0. {
                    continue;
                }
                for (overlap, col) in overlaps.iter_mut().zip(atom.columns()) {
                    *overlap = col
                        .iter()
                        .enumerate()
                        .map(|(ig, b)| b.conj() * orbitals[(ig, band)])
                        .sum();
                }
                for (c, metric) in metric.iter().enumerate() {
                    let weighted: Complex64 = (0..nc)
                        .map(|d| atom.coupling_hartree()[c * nc + d] * overlaps[d])
                        .sum();
                    for (a, row) in out.iter_mut().enumerate() {
                        for (b, v) in row.iter_mut().enumerate() {
                            let deriv: Complex64 = basis.momenta().iter().enumerate().try_fold(
                                Complex64::default(),
                                |sum, (ig, _g)| {
                                    let phase = phases[ig];
                                    let d = metric[ig][a][b] * phase
                                        - if a == b {
                                            0.5 * atom.columns()[c][ig]
                                        } else {
                                            Complex64::default()
                                        };
                                    Ok::<_, HamiltonianError>(sum + d.conj() * orbitals[(ig, band)])
                                },
                            )?;
                            *v += 2. * f * (deriv.conj() * weighted).re / volume;
                        }
                    }
                }
            }
        }
    }
    if out.iter().flatten().any(|v| !v.is_finite()) {
        return Err(HamiltonianError::NonFinite);
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn solid_harmonic_gradients_include_poles_and_all_supported_channels() {
        for q in [[0.3, -0.4, 0.7], [0., 0., 1.], [0., 1., 0.]] {
            for l in 0..=3 {
                for m in -(l as isize)..=l as isize {
                    let (v, d) = harmonic(l, m, q);
                    let reference = dft_upf::radial::real_spherical_harmonic(l, m, q).unwrap();
                    assert!((v - reference).abs() < 1e-14);
                    for (a, derivative) in d.iter().enumerate() {
                        // Avoid the legacy polar acos cancellation in the oracle;
                        // off-axis derivatives still test every channel and axis.
                        if q[0] == 0. && q[1] == 0. {
                            continue;
                        }
                        let mut p = q;
                        let mut n = q;
                        p[a] += 1e-5;
                        n[a] -= 1e-5;
                        let fd = (dft_upf::radial::real_spherical_harmonic(l, m, p).unwrap()
                            - dft_upf::radial::real_spherical_harmonic(l, m, n).unwrap())
                            / 2e-5;
                        assert!((fd - derivative).abs() < 1e-8, "l={l} m={m} axis={a}");
                    }
                }
            }
        }
    }
}
