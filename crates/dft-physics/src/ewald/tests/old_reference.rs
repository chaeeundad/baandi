// Verbatim copy of the Ewald evaluator of commit eecd57f (before INIT-SCALE-3 and
// INIT-EWALD-FAST), renamed with an old_ prefix: the bit reference of the
// reference scheme and, at a tight tolerance, the accuracy reference of the
// balanced scheme. Test-only.
use super::super::*;

fn old_covering_radius(a: &Matrix3<f64>) -> f64 {
    (0..3).map(|i| a.column(i).norm()).sum::<f64>() / 2.0
}
fn old_bounds(dual: &Matrix3<f64>, radius: f64, shift: f64) -> Result<[i32; 3], PhysicsError> {
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
fn old_enumerate(
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

pub(super) fn old_energy(
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
        || ions
            .iter()
            .any(|i| !i.charge.is_finite() || i.fractional_position.iter().any(|x| !x.is_finite()))
    {
        return Err(PhysicsError::EwaldInput);
    }
    let a = matrix(cell.lattice_columns_bohr());
    let b = matrix(cell.reciprocal_columns_inverse_bohr());
    let rho_r = old_covering_radius(&a);
    let rho_g = old_covering_radius(&b);
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
    let rb = old_bounds(&b, r_cut, 1.0)?;
    let gb = old_bounds(&a, g_cut, 0.0)?;
    let positions: Vec<Vector3<f64>> = ions
        .iter()
        .map(|i| Vector3::from(i.fractional_position.map(|s| s.rem_euclid(1.0))))
        .collect();
    let mut real = 0.0;
    for (i, ion_i) in ions.iter().enumerate() {
        for (j, ion_j) in ions.iter().enumerate() {
            let delta = positions[i] - positions[j];
            old_enumerate(
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
    old_enumerate(
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

pub(super) fn old_forces(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
) -> Result<Vec<[f64; 3]>, PhysicsError> {
    let parameters = old_energy(cell, ions, options)?;
    let eta = parameters.eta;
    let a = matrix(cell.lattice_columns_bohr());
    let b = matrix(cell.reciprocal_columns_inverse_bohr());
    let rb = old_bounds(&b, parameters.real_cutoff_bohr, 1.0)?;
    let gb = old_bounds(&a, parameters.reciprocal_cutoff_inverse_bohr, 0.0)?;
    let positions: Vec<_> = ions
        .iter()
        .map(|ion| Vector3::from(ion.fractional_position.map(|x| x.rem_euclid(1.0))))
        .collect();
    let mut result = vec![[0.0; 3]; ions.len()];
    for (i, ion) in ions.iter().enumerate() {
        for (j, other) in ions.iter().enumerate() {
            old_enumerate(
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
    old_enumerate(
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

pub(super) fn old_stress(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
) -> Result<[[f64; 3]; 3], PhysicsError> {
    let parameters = old_energy(cell, ions, options)?;
    let eta = parameters.eta;
    let omega = cell.volume_bohr3();
    let a = matrix(cell.lattice_columns_bohr());
    let b = matrix(cell.reciprocal_columns_inverse_bohr());
    let rb = old_bounds(&b, parameters.real_cutoff_bohr, 1.0)?;
    let gb = old_bounds(&a, parameters.reciprocal_cutoff_inverse_bohr, 0.0)?;
    let positions: Vec<_> = ions
        .iter()
        .map(|ion| Vector3::from(ion.fractional_position.map(|x| x.rem_euclid(1.0))))
        .collect();
    let mut result = [[0.0; 3]; 3];
    for (i, ion) in ions.iter().enumerate() {
        for (j, other) in ions.iter().enumerate() {
            old_enumerate(
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
    old_enumerate(
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
