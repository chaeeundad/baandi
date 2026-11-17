//! Point-charge Ewald sum with a uniform neutralizing background.
use crate::PhysicsError;
use dft_geometry::Cell;
use dft_units::Hartree;
use nalgebra::{Matrix3, Vector3};
use num_complex::Complex64;
use std::f64::consts::{PI, TAU};

#[derive(Clone, Copy, Debug)]
pub struct Ion {
    pub charge: f64,
    pub fractional_position: [f64; 3],
}

#[derive(Clone, Copy, Debug)]
pub struct EwaldOptions {
    /// Inverse Bohr. None chooses sqrt(pi)/Omega^(1/3).
    pub eta: Option<f64>,
    /// Absolute truncation budget = tolerance * (sum |Z|)^2/Omega^(1/3).
    /// This is a scale-relative bound, not relative to a possibly zero energy.
    pub tolerance: f64,
}
impl Default for EwaldOptions {
    fn default() -> Self {
        Self {
            eta: None,
            tolerance: 1e-13,
        }
    }
}

#[derive(Clone, Debug)]
pub struct EwaldEnergy {
    pub energy: Hartree,
    pub real: Hartree,
    pub reciprocal: Hartree,
    pub self_energy: Hartree,
    pub background: Hartree,
    pub eta: f64,
    pub real_cutoff_bohr: f64,
    pub reciprocal_cutoff_inverse_bohr: f64,
    /// Conservative analytic truncation bound; excludes floating-point error.
    pub tail_bound: Hartree,
}

fn matrix(a: [[f64; 3]; 3]) -> Matrix3<f64> {
    Matrix3::from_fn(|i, j| a[i][j])
}

/// Tensile-positive ion/background strain derivative in Ha/Bohr^3, at fixed
/// fractional ion positions. Eta is held fixed in the derivative: its response
/// cancels between the four Ewald terms in the infinite sums. The energy
/// evaluator's truncation spheres are reused; its energy tail bound is not a
/// stress error bound. This is an isolated contribution, not a total stress.
pub fn stress(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
) -> Result<[[f64; 3]; 3], PhysicsError> {
    let parameters = energy(cell, ions, options)?;
    let eta = parameters.eta;
    let omega = cell.volume_bohr3();
    let a = matrix(cell.lattice_columns_bohr());
    let b = matrix(cell.reciprocal_columns_inverse_bohr());
    let rb = bounds(&b, parameters.real_cutoff_bohr, 1.0)?;
    let gb = bounds(&a, parameters.reciprocal_cutoff_inverse_bohr, 0.0)?;
    let positions: Vec<_> = ions
        .iter()
        .map(|ion| Vector3::from(ion.fractional_position.map(|x| x.rem_euclid(1.0))))
        .collect();
    let mut result = [[0.0; 3]; 3];
    for (i, ion) in ions.iter().enumerate() {
        for (j, other) in ions.iter().enumerate() {
            enumerate(
                |m| {
                    if i == j && m == [0; 3] {
                        return Ok(());
                    }
                    let r = a * (positions[i] - positions[j] + Vector3::from(m.map(f64::from)));
                    let distance = r.norm();
                    if distance <= parameters.real_cutoff_bohr {
                        let factor = -0.5 * ion.charge * other.charge / omega
                            * (libm::erfc(eta * distance) / distance.powi(3)
                                + 2.0 * eta / PI.sqrt() * (-(eta * distance).powi(2)).exp()
                                    / distance.powi(2));
                        for (a, row) in result.iter_mut().enumerate() {
                            for (b, value) in row.iter_mut().enumerate() {
                                *value += factor * r[a] * r[b];
                            }
                        }
                    }
                    Ok(())
                },
                rb,
            )?;
        }
    }
    enumerate(
        |m| {
            if m == [0; 3] {
                return Ok(());
            }
            let mv = Vector3::from(m.map(f64::from));
            let g = b * mv;
            let g2 = g.norm_squared();
            if g2 <= parameters.reciprocal_cutoff_inverse_bohr.powi(2) {
                let structure: Complex64 = ions
                    .iter()
                    .zip(&positions)
                    .map(|(ion, s)| Complex64::from_polar(ion.charge, TAU * mv.dot(s)))
                    .sum();
                let weight = 2.0 * PI / omega.powi(2) * (-g2 / (4.0 * eta * eta)).exp() / g2
                    * structure.norm_sqr();
                for (a, row) in result.iter_mut().enumerate() {
                    for (b, value) in row.iter_mut().enumerate() {
                        *value += weight
                            * ((2.0 / g2 + 0.5 / eta.powi(2)) * g[a] * g[b]
                                - if a == b { 1.0 } else { 0.0 });
                    }
                }
            }
            Ok(())
        },
        gb,
    )?;
    // The self term has no explicit fixed-eta strain derivative. Background
    // scales as 1/Omega and must not be mistaken for a zero-mode constant.
    for (a, row) in result.iter_mut().enumerate() {
        row[a] -= parameters.background.0 / omega;
    }
    if result.iter().flatten().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(result)
}
fn covering_radius(a: &Matrix3<f64>) -> f64 {
    (0..3).map(|i| a.column(i).norm()).sum::<f64>() / 2.0
}
fn bounds(dual: &Matrix3<f64>, radius: f64, shift: f64) -> Result<[i32; 3], PhysicsError> {
    let mut bounds = [0; 3];
    for (i, b) in bounds.iter_mut().enumerate() {
        let value = (dual.column(i).norm() * radius / TAU + shift).ceil() + 1.0;
        if !value.is_finite() || value >= i32::MAX as f64 {
            return Err(PhysicsError::Unrepresentable);
        }
        *b = value as i32;
    }
    Ok(bounds)
}
fn enumerate(
    mut visit: impl FnMut([i32; 3]) -> Result<(), PhysicsError>,
    b: [i32; 3],
) -> Result<(), PhysicsError> {
    for i in -b[0]..=b[0] {
        for j in -b[1]..=b[1] {
            for k in -b[2]..=b[2] {
                visit([i, j, k])?;
            }
        }
    }
    Ok(())
}

/// Parameters shared by the serial evaluation and the distributed plan.
struct Prepared {
    eta: f64,
    a: Matrix3<f64>,
    b: Matrix3<f64>,
    omega: f64,
    r_cut: f64,
    g_cut: f64,
    tail: f64,
    rb: [i32; 3],
    gb: [i32; 3],
    positions: Vec<Vector3<f64>>,
}

fn prepare(cell: &Cell, ions: &[Ion], options: EwaldOptions) -> Result<Prepared, PhysicsError> {
    let length = cell.volume_bohr3().cbrt();
    let eta = options.eta.unwrap_or(PI.sqrt() / length);
    if !eta.is_finite()
        || eta <= 0.0
        || !options.tolerance.is_finite()
        || options.tolerance <= 0.0
        || options.tolerance >= 1.0
        || ions
            .iter()
            .any(|i| !i.charge.is_finite() || i.fractional_position.iter().any(|x| !x.is_finite()))
    {
        return Err(PhysicsError::EwaldInput);
    }
    let a = matrix(cell.lattice_columns_bohr());
    let b = matrix(cell.reciprocal_columns_inverse_bohr());
    let rho_r = covering_radius(&a);
    let rho_g = covering_radius(&b);
    let omega = cell.volume_bohr3();
    let omega_g = TAU.powi(3) / omega;
    // Bounds are per (sum |Z|)^2, avoiding division by zero for empty/zero charges.
    let budget = options.tolerance / length;
    let mut t = (-options.tolerance.ln()).sqrt();
    let (r_cut, g_cut, tail) = loop {
        let u = t / eta;
        let v = 2.0 * eta * t;
        let real_bound =
            PI.sqrt() / (omega * eta.powi(3) * u) * (1.0 + rho_r / u).powi(2) * (-t * t).exp();
        let reciprocal_bound = (2.0 * PI / omega)
            * (8.0 * PI * eta * eta / (omega_g * v))
            * (1.0 + rho_g / v).powi(2)
            * (-t * t).exp();
        if !u.is_finite()
            || !v.is_finite()
            || !real_bound.is_finite()
            || !reciprocal_bound.is_finite()
            || budget == 0.0
        {
            return Err(PhysicsError::Unrepresentable);
        }
        if real_bound + reciprocal_bound <= budget {
            break (
                u + 2.0 * rho_r,
                v + 2.0 * rho_g,
                real_bound + reciprocal_bound,
            );
        }
        t += 1.0;
    };
    let rb = bounds(&b, r_cut, 1.0)?;
    let gb = bounds(&a, g_cut, 0.0)?;
    let positions: Vec<Vector3<f64>> = ions
        .iter()
        .map(|i| Vector3::from(i.fractional_position.map(|s| s.rem_euclid(1.0))))
        .collect();
    Ok(Prepared {
        eta,
        a,
        b,
        omega,
        r_cut,
        g_cut,
        tail,
        rb,
        gb,
        positions,
    })
}

/// Real-space term of ion pair (i, j) and lattice translation m; None when the
/// candidate is excluded (self term) or outside the cutoff sphere.
fn real_term(
    p: &Prepared,
    ions: &[Ion],
    i: usize,
    j: usize,
    m: [i32; 3],
) -> Result<Option<f64>, PhysicsError> {
    if i == j && m == [0; 3] {
        return Ok(None);
    }
    let delta = p.positions[i] - p.positions[j];
    let distance = (p.a * (delta + Vector3::from(m.map(f64::from)))).norm();
    if distance == 0.0 {
        return Err(PhysicsError::CoincidentIons);
    }
    Ok((distance <= p.r_cut)
        .then(|| 0.5 * ions[i].charge * ions[j].charge * libm::erfc(p.eta * distance) / distance))
}

/// Reciprocal term of reciprocal lattice point m; None for m = 0 or outside the
/// cutoff sphere.
fn reciprocal_term(p: &Prepared, ions: &[Ion], m: [i32; 3]) -> Option<f64> {
    if m == [0; 3] {
        return None;
    }
    let mv = Vector3::from(m.map(f64::from));
    let g2 = (p.b * mv).norm_squared();
    // The serial predicate exactly (a NaN g2 is skipped, as before).
    if !(g2 <= p.g_cut * p.g_cut) {
        return None;
    }
    let structure: Complex64 = ions
        .iter()
        .zip(&p.positions)
        .map(|(ion, s)| Complex64::from_polar(ion.charge, TAU * mv.dot(s)))
        .sum();
    Some((2.0 * PI / p.omega) * (-g2 / (4.0 * p.eta * p.eta)).exp() / g2 * structure.norm_sqr())
}

fn finish(
    p: &Prepared,
    ions: &[Ion],
    real: f64,
    reciprocal: f64,
) -> Result<EwaldEnergy, PhysicsError> {
    let eta = p.eta;
    let omega = p.omega;
    let total_charge: f64 = ions.iter().map(|i| i.charge).sum();
    let abs_charge: f64 = ions.iter().map(|i| i.charge.abs()).sum();
    let self_energy = -eta / PI.sqrt() * ions.iter().map(|i| i.charge * i.charge).sum::<f64>();
    let background = -PI * total_charge * total_charge / (2.0 * omega * eta * eta);
    let energy = real + reciprocal + self_energy + background;
    if !energy.is_finite() || !(p.tail * abs_charge * abs_charge).is_finite() {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(EwaldEnergy {
        energy: Hartree(energy),
        real: Hartree(real),
        reciprocal: Hartree(reciprocal),
        self_energy: Hartree(self_energy),
        background: Hartree(background),
        eta,
        real_cutoff_bohr: p.r_cut,
        reciprocal_cutoff_inverse_bohr: p.g_cut,
        tail_bound: Hartree(p.tail * abs_charge * abs_charge),
    })
}

pub fn energy(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
) -> Result<EwaldEnergy, PhysicsError> {
    let p = prepare(cell, ions, options)?;
    let mut real = 0.0;
    for i in 0..ions.len() {
        for j in 0..ions.len() {
            enumerate(
                |m| {
                    if let Some(term) = real_term(&p, ions, i, j, m)? {
                        real += term;
                    }
                    Ok(())
                },
                p.rb,
            )?;
        }
    }
    let mut reciprocal = 0.0;
    enumerate(
        |m| {
            if let Some(term) = reciprocal_term(&p, ions, m) {
                reciprocal += term;
            }
            Ok(())
        },
        p.gb,
    )?;
    finish(&p, ions, real, reciprocal)
}

/// Points of the enumeration box `b` (the order of `enumerate`).
fn box_points(b: [i32; 3]) -> Result<usize, PhysicsError> {
    b.iter()
        .map(|&e| usize::try_from(2 * i64::from(e) + 1).ok())
        .try_fold(1usize, |acc, n| n.and_then(|n| acc.checked_mul(n)))
        .ok_or(PhysicsError::Unrepresentable)
}

/// Point `k` of the enumeration box `b`, in the order of `enumerate`.
fn box_point(b: [i32; 3], k: usize) -> [i32; 3] {
    let n1 = (2 * b[1] + 1) as usize;
    let n2 = (2 * b[2] + 1) as usize;
    [
        (k / (n1 * n2)) as i32 - b[0],
        ((k / n2) % n1) as i32 - b[1],
        (k % n2) as i32 - b[2],
    ]
}

/// The Ewald sum split into independent terms (INIT-SCALE-3): the candidates of
/// the real-space sum (ion pair i, j in order, then the translations in the
/// order of the serial enumeration) and of the reciprocal sum are numbered;
/// `real_terms` / `reciprocal_terms` return the accepted terms of a candidate
/// range in that order; `finish` adds a complete term list sequentially from
/// +0.0. Skipped candidates add nothing in the serial loop, so the terms of
/// contiguous candidate ranges, concatenated in range order and summed by
/// `finish`, give [`energy`] bit for bit for any partition.
pub struct EwaldPlan<'a> {
    ions: &'a [Ion],
    p: Prepared,
    real_box: usize,
    real_candidates: usize,
    reciprocal_candidates: usize,
}
impl<'a> EwaldPlan<'a> {
    pub fn new(cell: &Cell, ions: &'a [Ion], options: EwaldOptions) -> Result<Self, PhysicsError> {
        let p = prepare(cell, ions, options)?;
        let real_box = box_points(p.rb)?;
        let real_candidates = ions
            .len()
            .checked_mul(ions.len())
            .and_then(|n| n.checked_mul(real_box))
            .ok_or(PhysicsError::Unrepresentable)?;
        let reciprocal_candidates = box_points(p.gb)?;
        Ok(Self {
            ions,
            p,
            real_box,
            real_candidates,
            reciprocal_candidates,
        })
    }
    pub fn real_candidates(&self) -> usize {
        self.real_candidates
    }
    pub fn reciprocal_candidates(&self) -> usize {
        self.reciprocal_candidates
    }
    /// Appends the accepted real-space terms of candidates `range` to `out`.
    pub fn real_terms(
        &self,
        range: std::ops::Range<usize>,
        out: &mut Vec<f64>,
    ) -> Result<(), PhysicsError> {
        let n = self.ions.len();
        for c in range.start..range.end.min(self.real_candidates) {
            let pair = c / self.real_box;
            let m = box_point(self.p.rb, c % self.real_box);
            if let Some(term) = real_term(&self.p, self.ions, pair / n, pair % n, m)? {
                out.push(term);
            }
        }
        Ok(())
    }
    /// Appends the accepted reciprocal terms of candidates `range` to `out`.
    pub fn reciprocal_terms(&self, range: std::ops::Range<usize>, out: &mut Vec<f64>) {
        for c in range.start..range.end.min(self.reciprocal_candidates) {
            if let Some(term) = reciprocal_term(&self.p, self.ions, box_point(self.p.gb, c)) {
                out.push(term);
            }
        }
    }
    /// The energy from complete term lists (all candidates, serial order).
    pub fn finish(&self, real: &[f64], reciprocal: &[f64]) -> Result<EwaldEnergy, PhysicsError> {
        self.finish_sums(accumulate(0.0, real), accumulate(0.0, reciprocal))
    }
    /// The energy from the two sums accumulated in the serial order (for
    /// example by `accumulate` chained over contiguous candidate blocks).
    pub fn finish_sums(&self, real: f64, reciprocal: f64) -> Result<EwaldEnergy, PhysicsError> {
        finish(&self.p, self.ions, real, reciprocal)
    }
}

/// `start` plus `terms` added one at a time in order: the serial loop's
/// association, so chaining over the blocks of a partition from +0.0 gives the
/// serial sum bit for bit.
pub fn accumulate(start: f64, terms: &[f64]) -> f64 {
    let mut sum = start;
    for t in terms {
        sum += t;
    }
    sum
}

/// Cartesian -dE/dR in Ha/Bohr, including real and reciprocal terms.
/// Self and uniform-background energies have no position derivative at fixed cell.
/// Uses exactly the energy evaluator's truncation spheres and eta.
pub fn forces(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
) -> Result<Vec<[f64; 3]>, PhysicsError> {
    let parameters = energy(cell, ions, options)?;
    let eta = parameters.eta;
    let a = matrix(cell.lattice_columns_bohr());
    let b = matrix(cell.reciprocal_columns_inverse_bohr());
    let rb = bounds(&b, parameters.real_cutoff_bohr, 1.0)?;
    let gb = bounds(&a, parameters.reciprocal_cutoff_inverse_bohr, 0.0)?;
    let positions: Vec<_> = ions
        .iter()
        .map(|ion| Vector3::from(ion.fractional_position.map(|x| x.rem_euclid(1.0))))
        .collect();
    let mut result = vec![[0.0; 3]; ions.len()];
    for (i, ion) in ions.iter().enumerate() {
        for (j, other) in ions.iter().enumerate() {
            enumerate(
                |m| {
                    if i == j && m == [0; 3] {
                        return Ok(());
                    }
                    let r = a * (positions[i] - positions[j] + Vector3::from(m.map(f64::from)));
                    let distance = r.norm();
                    if distance <= parameters.real_cutoff_bohr {
                        let factor = ion.charge
                            * other.charge
                            * (libm::erfc(eta * distance) / distance.powi(3)
                                + 2.0 * eta / PI.sqrt() * (-(eta * distance).powi(2)).exp()
                                    / distance.powi(2));
                        for axis in 0..3 {
                            result[i][axis] += factor * r[axis];
                        }
                    }
                    Ok(())
                },
                rb,
            )?;
        }
    }
    enumerate(
        |m| {
            if m == [0; 3] {
                return Ok(());
            }
            let mv = Vector3::from(m.map(f64::from));
            let g = b * mv;
            let g2 = g.norm_squared();
            if g2 <= parameters.reciprocal_cutoff_inverse_bohr.powi(2) {
                let structure: Complex64 = ions
                    .iter()
                    .zip(&positions)
                    .map(|(ion, s)| Complex64::from_polar(ion.charge, TAU * mv.dot(s)))
                    .sum();
                let prefactor =
                    4.0 * PI / cell.volume_bohr3() * (-g2 / (4.0 * eta * eta)).exp() / g2;
                for (i, ion) in ions.iter().enumerate() {
                    let phase = Complex64::from_polar(ion.charge, TAU * mv.dot(&positions[i]));
                    for axis in 0..3 {
                        result[i][axis] += prefactor * (structure.conj() * phase).im * g[axis];
                    }
                }
            }
            Ok(())
        },
        gb,
    )?;
    if result.iter().flatten().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    // Verbatim copy of the evaluator before the INIT-SCALE-3 refactor (commit
    // a7205c9), the bit reference of the serial `energy` and of the plan.
    fn old_energy(
        cell: &Cell,
        ions: &[Ion],
        options: EwaldOptions,
    ) -> Result<EwaldEnergy, PhysicsError> {
        let length = cell.volume_bohr3().cbrt();
        let eta = options.eta.unwrap_or(PI.sqrt() / length);
        if !eta.is_finite()
            || eta <= 0.0
            || !options.tolerance.is_finite()
            || options.tolerance <= 0.0
            || options.tolerance >= 1.0
            || ions.iter().any(|i| {
                !i.charge.is_finite() || i.fractional_position.iter().any(|x| !x.is_finite())
            })
        {
            return Err(PhysicsError::EwaldInput);
        }
        let a = matrix(cell.lattice_columns_bohr());
        let b = matrix(cell.reciprocal_columns_inverse_bohr());
        let rho_r = covering_radius(&a);
        let rho_g = covering_radius(&b);
        let omega = cell.volume_bohr3();
        let omega_g = TAU.powi(3) / omega;
        // Bounds are per (sum |Z|)^2, avoiding division by zero for empty/zero charges.
        let budget = options.tolerance / length;
        let mut t = (-options.tolerance.ln()).sqrt();
        let (r_cut, g_cut, tail) = loop {
            let u = t / eta;
            let v = 2.0 * eta * t;
            let real_bound =
                PI.sqrt() / (omega * eta.powi(3) * u) * (1.0 + rho_r / u).powi(2) * (-t * t).exp();
            let reciprocal_bound = (2.0 * PI / omega)
                * (8.0 * PI * eta * eta / (omega_g * v))
                * (1.0 + rho_g / v).powi(2)
                * (-t * t).exp();
            if !u.is_finite()
                || !v.is_finite()
                || !real_bound.is_finite()
                || !reciprocal_bound.is_finite()
                || budget == 0.0
            {
                return Err(PhysicsError::Unrepresentable);
            }
            if real_bound + reciprocal_bound <= budget {
                break (
                    u + 2.0 * rho_r,
                    v + 2.0 * rho_g,
                    real_bound + reciprocal_bound,
                );
            }
            t += 1.0;
        };
        let rb = bounds(&b, r_cut, 1.0)?;
        let gb = bounds(&a, g_cut, 0.0)?;
        let positions: Vec<Vector3<f64>> = ions
            .iter()
            .map(|i| Vector3::from(i.fractional_position.map(|s| s.rem_euclid(1.0))))
            .collect();
        let mut real = 0.0;
        for (i, ion_i) in ions.iter().enumerate() {
            for (j, ion_j) in ions.iter().enumerate() {
                let delta = positions[i] - positions[j];
                enumerate(
                    |m| {
                        if i == j && m == [0; 3] {
                            return Ok(());
                        }
                        let distance = (a * (delta + Vector3::from(m.map(f64::from)))).norm();
                        if distance == 0.0 {
                            return Err(PhysicsError::CoincidentIons);
                        }
                        if distance <= r_cut {
                            real += 0.5 * ion_i.charge * ion_j.charge * libm::erfc(eta * distance)
                                / distance;
                        }
                        Ok(())
                    },
                    rb,
                )?;
            }
        }
        let mut reciprocal = 0.0;
        enumerate(
            |m| {
                if m == [0; 3] {
                    return Ok(());
                }
                let mv = Vector3::from(m.map(f64::from));
                let g2 = (b * mv).norm_squared();
                if g2 <= g_cut * g_cut {
                    let structure: Complex64 = ions
                        .iter()
                        .zip(&positions)
                        .map(|(ion, s)| Complex64::from_polar(ion.charge, TAU * mv.dot(s)))
                        .sum();
                    reciprocal += (2.0 * PI / omega) * (-g2 / (4.0 * eta * eta)).exp() / g2
                        * structure.norm_sqr();
                }
                Ok(())
            },
            gb,
        )?;
        let total_charge: f64 = ions.iter().map(|i| i.charge).sum();
        let abs_charge: f64 = ions.iter().map(|i| i.charge.abs()).sum();
        let self_energy = -eta / PI.sqrt() * ions.iter().map(|i| i.charge * i.charge).sum::<f64>();
        let background = -PI * total_charge * total_charge / (2.0 * omega * eta * eta);
        let energy = real + reciprocal + self_energy + background;
        if !energy.is_finite() || !(tail * abs_charge * abs_charge).is_finite() {
            return Err(PhysicsError::Unrepresentable);
        }
        Ok(EwaldEnergy {
            energy: Hartree(energy),
            real: Hartree(real),
            reciprocal: Hartree(reciprocal),
            self_energy: Hartree(self_energy),
            background: Hartree(background),
            eta,
            real_cutoff_bohr: r_cut,
            reciprocal_cutoff_inverse_bohr: g_cut,
            tail_bound: Hartree(tail * abs_charge * abs_charge),
        })
    }

    fn cases() -> Vec<(Cell, Vec<Ion>)> {
        let shapes = [
            [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]],
            [[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]],
            [[4., 0., 0.], [0., 5., 0.], [0., 0., 14.]],
            [[7., 0., 0.], [0., 8., 0.], [0.3, 0., 18.]],
        ];
        let mut out = Vec::new();
        for rows in shapes {
            for atoms in [1usize, 2, 5, 9] {
                let ions = (0..atoms)
                    .map(|i| Ion {
                        charge: if i % 3 == 2 { 6.0 } else { 4.0 },
                        fractional_position: [
                            ((i * 7 + 1) % 23) as f64 / 23. - 0.5,
                            ((i * 11 + 2) % 29) as f64 / 29.,
                            ((i * 13 + 3) % 31) as f64 / 31. + 1.0,
                        ],
                    })
                    .collect();
                out.push((Cell::from_user_rows_bohr(rows).unwrap(), ions));
            }
        }
        out
    }

    fn bits(e: &EwaldEnergy) -> [u64; 9] {
        [
            e.energy.0.to_bits(),
            e.real.0.to_bits(),
            e.reciprocal.0.to_bits(),
            e.self_energy.0.to_bits(),
            e.background.0.to_bits(),
            e.eta.to_bits(),
            e.real_cutoff_bohr.to_bits(),
            e.reciprocal_cutoff_inverse_bohr.to_bits(),
            e.tail_bound.0.to_bits(),
        ]
    }

    /// The refactored serial evaluator and every contiguous partition of the
    /// plan's candidates into 1 to 7 parts reproduce the old evaluator bit for
    /// bit (cubic, skewed, elongated, slab; 1 to 9 ions; mixed charges;
    /// positions outside [0, 1)).
    #[test]
    fn serial_and_partitioned_ewald_match_the_previous_evaluator_bit_for_bit() {
        let options = EwaldOptions::default();
        for (cell, ions) in cases() {
            let want = bits(&old_energy(&cell, &ions, options).unwrap());
            assert_eq!(bits(&energy(&cell, &ions, options).unwrap()), want);
            let plan = EwaldPlan::new(&cell, &ions, options).unwrap();
            for parts in 1..=7usize {
                let split = |n: usize, k: usize| {
                    let first = (0..k).map(|q| n / parts + usize::from(q < n % parts)).sum();
                    first..first + n / parts + usize::from(k < n % parts)
                };
                let (mut real, mut reciprocal) = (Vec::new(), Vec::new());
                // Chained form (the engine's): each block folded into the
                // running sums of the blocks before it.
                let (mut r, mut g) = (0.0, 0.0);
                for k in 0..parts {
                    let (mut rk, mut gk) = (Vec::new(), Vec::new());
                    plan.real_terms(split(plan.real_candidates(), k), &mut rk)
                        .unwrap();
                    plan.reciprocal_terms(split(plan.reciprocal_candidates(), k), &mut gk);
                    r = accumulate(r, &rk);
                    g = accumulate(g, &gk);
                    real.extend(rk);
                    reciprocal.extend(gk);
                }
                assert_eq!(
                    bits(&plan.finish(&real, &reciprocal).unwrap()),
                    want,
                    "{parts}"
                );
                assert_eq!(bits(&plan.finish_sums(r, g).unwrap()), want, "{parts}");
            }
        }
    }

    /// Coincident ions fail in the serial evaluator and in the plan.
    #[test]
    fn coincident_ions_fail_in_both_evaluators() {
        let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
        let ion = Ion {
            charge: 4.0,
            fractional_position: [0.1, 0.2, 0.3],
        };
        let ions = [ion, ion];
        let options = EwaldOptions::default();
        assert!(matches!(
            energy(&cell, &ions, options),
            Err(PhysicsError::CoincidentIons)
        ));
        let plan = EwaldPlan::new(&cell, &ions, options).unwrap();
        let mut out = Vec::new();
        assert!(matches!(
            plan.real_terms(0..plan.real_candidates(), &mut out),
            Err(PhysicsError::CoincidentIons)
        ));
    }
}
