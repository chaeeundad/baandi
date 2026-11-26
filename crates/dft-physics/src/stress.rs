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
    local_or_core_cached(basis, volume, species, field, core, None)
}
/// GD1: [`local_or_core`] reading the per-shell radial values and slopes from
/// `shells` when present (identical expressions computed elsewhere); a
/// missing pair is computed here. Bit-identical to [`local_or_core`].
pub fn local_or_core_cached(
    basis: &PlaneWaveBasis,
    volume: f64,
    species: &[Species<'_>],
    field: &[Complex64],
    core: bool,
    shells: Option<&crate::radial_shells::RadialShellValues>,
) -> Result<Tensor, HamiltonianError> {
    use crate::radial_shells::ShellKind;
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
            let shared = if q == 0. {
                None
            } else if core {
                shells.and_then(|c| {
                    Some((
                        c.get(pp, ShellKind::CoreValue, q)? / volume,
                        c.get(pp, ShellKind::CoreSlope, q)?,
                    ))
                })
            } else {
                shells.and_then(|c| {
                    Some((
                        c.get(pp, ShellKind::LocalValue, q)?,
                        c.get(pp, ShellKind::LocalSlope, q)?,
                    ))
                })
            };
            let (v, prime) = if let Some(hit) = cache.get(&key) {
                *hit
            } else if let Some(pair) = shared {
                if core && pp.core_density().is_none() {
                    continue;
                }
                cache.insert(key, pair);
                pair
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
/// d/dq of `Pseudopotential::projector_transform(index, q)`:
/// `4 pi int r^2 u(r) j_l'(q r) dr` over the projector samples (u = r beta),
/// the same index-coordinate quadrature. Exactly the expression of the
/// frozen-Miller nonlocal stress (`nonlocal`), so a table built from it
/// reproduces that quadrature at its nodes bit for bit.
pub fn projector_transform_slope(
    pp: &dft_upf::Pseudopotential,
    index: usize,
    q: f64,
) -> Result<f64, HamiltonianError> {
    let p = pp
        .projectors()
        .get(index)
        .ok_or(RadialError::Argument("projector index"))?;
    let l = p.l();
    Ok(4.
        * PI
        * integrate_prefix(
            pp.mesh(),
            p.radial_times_r()
                .iter()
                .zip(pp.mesh().r())
                .map(|(&u, &r)| Ok(r * r * u * bessel_prime(l, q * r)?))
                .collect::<Result<Vec<_>, RadialError>>()?,
        )?)
}

/// GD1: uniform radial tables of one species' projector transforms
/// (`value[i]`, `Pseudopotential::projector_transform`) and their q
/// derivatives (`slope[i]`, [`projector_transform_slope`]), one per projector
/// `i`. An entry outside a table's range falls back to the quadrature.
#[derive(Debug, Clone)]
pub struct ProjectorRadialTables {
    pub value: Vec<dft_upf::radial::RadialTransformTable>,
    pub slope: Vec<dft_upf::radial::RadialTransformTable>,
}

/// GD1: nonlocal ionic forces (`with_stress == false`) or the frozen-Miller
/// nonlocal strain derivative (`with_stress == true`) of one k block.
///
/// Same physics as [`crate::forces::nonlocal`] and [`nonlocal`], reorganised
/// so that the band sum is done by two GEMMs per atom over the occupied bands:
/// `O = beta^H psi` (nc x nb), `W_cn = f_n sum_d D_cd O_dn`, and the
/// band-summed weight field `Y_c(G) = sum_n conj(psi_n(G)) W_cn` (one GEMM
/// `psi W^H`, conjugated). Then, with `beta_c(G) = b_c(G) e^{-i q.tau}` and
/// its momentum gradient `D^b_c(G) = d beta_c / d q_b`,
/// `F_a = -2 Im sum_G q_a sum_c beta_c Y_c` and
/// `sigma_ab = (2/V) Re[-sum_G q_a sum_c D^b_c Y_c - delta_ab/2 sum_cn conj(O_cn) W_cn]`.
/// The work is 2 ng nc nb multiply-adds per atom plus O(ng nc) for the field
/// sums, instead of (4 or 10) ng nc nb. Only the reduction order differs from
/// the per-G scalar loops (rounding class). `tables[s]` (when given) supplies
/// the radial values and slopes of `species[s]`; otherwise, and for |q|
/// outside a table, the quadrature. Occupations include the k weight and the
/// spin capacity. Returns the forces in species/position order and the stress
/// tensor; the quantity not requested is zero.
pub fn nonlocal_k_derivatives(
    basis: &PlaneWaveBasis,
    volume: f64,
    species: &[Species<'_>],
    tables: &[Option<&ProjectorRadialTables>],
    orbitals: &ComplexMatrix,
    occupations: &[f64],
    with_stress: bool,
) -> Result<(Vec<[f64; 3]>, Tensor), HamiltonianError> {
    use dft_numerics::linalg::{MatrixOp, complex_gemm_slices_into, gemm};
    let ng = basis.momenta().len();
    if orbitals.rows() != ng || orbitals.cols() != occupations.len() {
        return Err(HamiltonianError::Shape);
    }
    if !volume.is_finite() || volume <= 0. {
        return Err(HamiltonianError::Volume);
    }
    if occupations.iter().any(|f| !f.is_finite() || *f < 0.) {
        return Err(HamiltonianError::NonFinite);
    }
    let atoms_total: usize = species.iter().map(|s| s.positions_bohr.len()).sum();
    let mut forces = vec![[0.; 3]; atoms_total];
    let mut out = [[0.; 3]; 3];
    let occupied: Vec<usize> = (0..occupations.len())
        .filter(|&b| occupations[b] != 0.)
        .collect();
    let nb = occupied.len();
    if nb == 0 || ng == 0 {
        return Ok((forces, out));
    }
    let gemm_error = |_| HamiltonianError::NonFinite;
    // Occupied bands, column-major (ng x nb).
    let mut psi = ComplexMatrix::zeros(ng, nb).map_err(gemm_error)?;
    for (j, &band) in occupied.iter().enumerate() {
        psi.column_mut(j).copy_from_slice(orbitals.column(band));
    }
    psi.validate_finite().map_err(gemm_error)?;
    let volume_root = volume.sqrt();
    // Per-potential factors without the structure phase: columns beta_c (nc)
    // and, with stress, D^x_c, D^y_c, D^z_c (3 nc).
    let mut factors: Vec<(usize, Vec<Complex64>)> = Vec::new();
    let mut harmonics = vec![0.; 16];
    let mut gradients = [(0., [0.; 3]); 16];
    let mut atom_index = 0;
    let mut left = Vec::new();
    let mut overlaps = Vec::new();
    for (si, s) in species.iter().enumerate() {
        let pp = s.potential;
        let np = pp.projectors().len();
        let channels: Vec<_> = pp
            .projectors()
            .iter()
            .enumerate()
            .flat_map(|(i, p)| (-(p.l() as isize)..=p.l() as isize).map(move |m| (i, p.l(), m)))
            .collect();
        let nc = channels.len();
        if nc == 0 {
            atom_index += s.positions_bohr.len();
            continue;
        }
        for i in 0..np {
            for j in 0..np {
                if pp.projectors()[i].l() != pp.projectors()[j].l()
                    && pp.dij_hartree()[i * np + j] != 0.
                {
                    return Err(HamiltonianError::AngularCoupling);
                }
            }
        }
        let ncol = if with_stress { 4 * nc } else { nc };
        let lmax = channels.iter().map(|c| c.1).max().unwrap_or(0);
        if lmax > 3 {
            return Err(RadialError::Argument("projector l > 3").into());
        }
        let keys: Vec<usize> = channels
            .iter()
            .map(|&(_, l, m)| (l * l + l).wrapping_add_signed(m))
            .collect();
        let potential_id = pp as *const _ as *const () as usize;
        if !factors.iter().any(|f| f.0 == potential_id) {
            let table = tables
                .get(si)
                .copied()
                .flatten()
                .filter(|t| t.value.len() == np && (!with_stress || t.slope.len() == np));
            let prefactor: Vec<Complex64> = channels
                .iter()
                .map(|&(_, l, _)| Complex64::new(0., -1.).powu(l as u32) / volume_root)
                .collect();
            let mut distinct: Vec<(usize, usize, isize)> = Vec::new();
            for (&key, &(_, l, m)) in keys.iter().zip(&channels) {
                if !distinct.iter().any(|d| d.0 == key) {
                    distinct.push((key, l, m));
                }
            }
            let mut f = vec![Complex64::default(); ng * ncol];
            let mut radial = vec![0.; np];
            let mut slope = vec![0.; np];
            type Pair = (Vec<f64>, Vec<f64>);
            let mut fallback: std::collections::HashMap<u64, Pair> =
                std::collections::HashMap::new();
            for (ig, g) in basis.momenta().iter().enumerate() {
                let q = g.norm_squared.sqrt();
                let strained = with_stress && q > 0.;
                let mut missing = table.is_none();
                if let Some(t) = table {
                    for i in 0..np {
                        let value = t.value[i].eval(q);
                        let prime = if strained {
                            t.slope[i].eval(q)
                        } else {
                            Some(0.)
                        };
                        match (value, prime) {
                            (Some(v), Some(p)) => {
                                radial[i] = v;
                                slope[i] = p;
                            }
                            _ => missing = true,
                        }
                    }
                }
                if missing {
                    let entry = match fallback.entry(q.to_bits()) {
                        std::collections::hash_map::Entry::Occupied(e) => e.into_mut(),
                        std::collections::hash_map::Entry::Vacant(e) => {
                            let mut v = vec![0.; np];
                            let mut p = vec![0.; np];
                            for i in 0..np {
                                v[i] = pp.projector_transform(i, q)?;
                                if strained {
                                    p[i] = projector_transform_slope(pp, i, q)?;
                                }
                            }
                            e.insert((v, p))
                        }
                    };
                    radial.copy_from_slice(&entry.0);
                    slope.copy_from_slice(&entry.1);
                }
                // Direction at zero is arbitrary: all l>0 radial transforms vanish.
                let direction = if g.norm_squared == 0. {
                    [0., 0., 1.]
                } else {
                    g.cartesian
                };
                dft_upf::radial::real_spherical_harmonics_upto(lmax, direction, &mut harmonics)?;
                if strained {
                    for &(key, l, m) in &distinct {
                        gradients[key] = harmonic(l, m, g.cartesian);
                    }
                }
                for (c, &(i, _, _)) in channels.iter().enumerate() {
                    f[c * ng + ig] = prefactor[c] * harmonics[keys[c]] * radial[i];
                    if strained {
                        let (y, dy) = gradients[keys[c]];
                        for (b, d) in dy.iter().enumerate() {
                            f[(nc + b * nc + c) * ng + ig] =
                                prefactor[c] * (slope[i] * y * g.cartesian[b] / q + radial[i] * d);
                        }
                    }
                }
            }
            factors.push((potential_id, f));
        }
        let f = &factors
            .iter()
            .find(|f| f.0 == potential_id)
            .expect("cached factors")
            .1;
        let mut coupling = vec![0.; nc * nc];
        for (a, &(i, l, m)) in channels.iter().enumerate() {
            for (b, &(j, ll, mm)) in channels.iter().enumerate() {
                if l == ll && m == mm {
                    coupling[a * nc + b] = pp.dij_hartree()[i * np + j];
                }
            }
        }
        left.resize(ng * nc, Complex64::default());
        overlaps.resize(nc * nb, Complex64::default());
        let mut weights = ComplexMatrix::zeros(nc, nb).map_err(gemm_error)?;
        for position in s.positions_bohr {
            let phases = basis
                .momenta()
                .iter()
                .map(|g| structure_factor(g.cartesian, std::slice::from_ref(position)))
                .collect::<Result<Vec<_>, _>>()?;
            for c in 0..nc {
                for ((z, factor), phase) in left[c * ng..(c + 1) * ng]
                    .iter_mut()
                    .zip(&f[c * ng..(c + 1) * ng])
                    .zip(&phases)
                {
                    *z = factor * phase;
                }
            }
            // O_cn = <beta_c|psi_n> at (c + n nc).
            complex_gemm_slices_into(&left, ng, nc, true, psi.as_slice(), nb, &mut overlaps)
                .map_err(gemm_error)?;
            // W_cn = f_n sum_d D_cd O_dn; E = sum_cn conj(O_cn) W_cn.
            let mut energy = Complex64::default();
            for (j, &band) in occupied.iter().enumerate() {
                let occupation = occupations[band];
                for p in 0..nc {
                    let w: Complex64 = (0..nc)
                        .map(|q| coupling[p * nc + q] * overlaps[q + j * nc])
                        .sum::<Complex64>()
                        * occupation;
                    weights[(p, j)] = w;
                    energy += overlaps[p + j * nc].conj() * w;
                }
            }
            // conj(Y_c(G)) = sum_n psi_n(G) conj(W_cn), (ng x nc).
            let field =
                gemm(&psi, MatrixOp::Identity, &weights, MatrixOp::Adjoint).map_err(gemm_error)?;
            let force = &mut forces[atom_index];
            atom_index += 1;
            if with_stress {
                let mut sum = [[Complex64::default(); 3]; 3];
                for (ig, g) in basis.momenta().iter().enumerate() {
                    let mut s_b = [Complex64::default(); 3];
                    for c in 0..nc {
                        let y = field[(ig, c)].conj();
                        for (b, s) in s_b.iter_mut().enumerate() {
                            *s += f[(nc + b * nc + c) * ng + ig] * y;
                        }
                    }
                    for (a, row) in sum.iter_mut().enumerate() {
                        for (b, v) in row.iter_mut().enumerate() {
                            *v += g.cartesian[a] * phases[ig] * s_b[b];
                        }
                    }
                }
                for (a, row) in out.iter_mut().enumerate() {
                    for (b, v) in row.iter_mut().enumerate() {
                        let mut z = -sum[a][b];
                        if a == b {
                            z -= 0.5 * energy;
                        }
                        *v += 2. * z.re / volume;
                    }
                }
            } else {
                let mut sum = [Complex64::default(); 3];
                for (ig, g) in basis.momenta().iter().enumerate() {
                    let s: Complex64 = (0..nc)
                        .map(|c| left[c * ng + ig] * field[(ig, c)].conj())
                        .sum();
                    for (axis, v) in sum.iter_mut().enumerate() {
                        *v += g.cartesian[axis] * s;
                    }
                }
                for (axis, component) in force.iter_mut().enumerate() {
                    *component -= 2. * sum[axis].im;
                }
            }
        }
    }
    if forces
        .iter()
        .flatten()
        .chain(out.iter().flatten())
        .any(|v| !v.is_finite())
    {
        return Err(HamiltonianError::NonFinite);
    }
    Ok((forces, out))
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
