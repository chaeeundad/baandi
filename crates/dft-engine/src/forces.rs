//! Post-SCF derivatives of E, or F at fixed temperature/electron constraints.
//! No occupation-response term is added: stationarity of F cancels it together
//! with orbital/density response. The supplied model must be the state model.
use crate::{
    ScfError, ScfOptions, ScfResult, SpeciesPotential, collinear::CollinearResult,
    finite_temperature::FiniteTResult,
};
use dft_basis::FftGridSpec;
use dft_geometry::{Atom, Cell};
use dft_physics::{
    data::*,
    ewald, forces as terms,
    hamiltonian::{NonlocalOperator, Species},
};

#[derive(Clone, Debug)]
pub struct DerivativeEntry {
    pub active: bool,
    pub energy_ha: Option<f64>,
    pub included_in: &'static str,
    pub dependencies: Vec<&'static str>,
    pub term_id: &'static str,
    pub owner: &'static str,
    pub force_available: bool,
    pub stress_available: bool,
    pub reason: &'static str,
}
/// Every physical term in the supported fixed-cell NC functional. Audit-only
/// band sums and double-counting integrals are not independent energy terms.
pub fn derivative_ledger() -> Vec<DerivativeEntry> {
    [
        (
            "kinetic",
            "physics::kinetic",
            "stationary fixed plane-wave basis",
        ),
        (
            "local",
            "physics::forces",
            "ionic structure-factor derivative",
        ),
        (
            "nonlocal",
            "physics::forces",
            "Hermitian projector derivative",
        ),
        (
            "hartree",
            "physics::hartree",
            "no explicit ionic dependence",
        ),
        (
            "xc",
            "physics::xc",
            "valence response cancels by stationarity",
        ),
        (
            "core",
            "physics::forces",
            "moving core against XC potential",
        ),
        ("ewald", "physics::ewald", "real and reciprocal derivatives"),
        (
            "entropy",
            "engine::occupations",
            "stationary free energy at fixed theta,N,M",
        ),
    ]
    .into_iter()
    .map(|(term_id, owner, reason)| DerivativeEntry {
        active: true,
        energy_ha: None,
        included_in: "E",
        dependencies: Vec::new(),
        term_id,
        owner,
        reason,
        force_available: true,
        stress_available: false,
    })
    .collect()
}
/// Check the complete active ledger before calculating any derivative. Hosts
/// adding a physical term must append its capability here, never omit it.
pub fn require_complete(ledger: &[DerivativeEntry]) -> Result<(), String> {
    let missing: Vec<_> = ledger
        .iter()
        .filter(|e| e.active && !e.force_available)
        .map(|e| e.term_id)
        .collect();
    if missing.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "observables.forces: incomplete derivative coverage: {}; disable this request or select a validated combination",
            missing.join(", ")
        ))
    }
}

#[derive(Debug)]
pub struct ForceResult {
    pub ledger: Vec<DerivativeEntry>,
    pub units: &'static str,
    pub quantity: &'static str,
    pub values: Vec<[f64; 3]>,
    pub ewald: Vec<[f64; 3]>,
    pub local: Vec<[f64; 3]>,
    pub nonlocal: Vec<[f64; 3]>,
    pub core: Vec<[f64; 3]>,
    pub net_force: [f64; 3],
    pub maximum_absolute_component: f64,
    pub elapsed: std::time::Duration,
    pub band_passes: usize,
}
/// Explicit borrowed stationary electronic state/model view. This API deliberately
/// does not accept a density-only checkpoint as a stationary electronic state.
/// Callers must retain the exact cell, atoms, potentials and options used to
/// produce `electronic`; matching dimensions alone do not establish identity.
pub enum ElectronicState<'a> {
    Fixed(&'a ScfResult),
    Collinear(&'a CollinearResult),
    FiniteTemperature(&'a FiniteTResult),
}
pub struct ForceState<'a> {
    pub cell: &'a Cell,
    pub atoms: &'a [Atom],
    pub potentials: &'a [SpeciesPotential<'a>],
    pub options: &'a ScfOptions,
    pub electronic: ElectronicState<'a>,
    pub density_floor: Option<f64>,
    /// Exact spatial operations used by the SCF, empty when symmetry is off.
    pub symmetry: &'a [dft_geometry::symmetry::SymmetryOperation],
}

pub fn forces(state: &ForceState<'_>) -> Result<ForceResult, ScfError> {
    if state.options.local_mesh != crate::LocalMesh::Baandi
        || state.options.radial_even_rule != crate::RadialEvenRule::Simpson3_8
    {
        return Err(ScfError::InvalidInput(
            "research pseudo conventions do not support analytic derivatives",
        ));
    }

    let start = std::time::Instant::now();
    use crate::force_identity::{ModelDescription, OccupationIdentity};
    let (identity, policy, density, points) = match &state.electronic {
        ElectronicState::Fixed(r) => (
            &r.force_identity,
            OccupationIdentity::new(&[&r.occupations], None, None),
            &r.density,
            r.k_points.iter().map(|k| k.point).collect(),
        ),
        ElectronicState::FiniteTemperature(r) => (
            &r.force_identity,
            OccupationIdentity::new(&[], Some(&r.occupations), None),
            &r.state.density,
            r.state.k_points.iter().map(|k| k.point).collect(),
        ),
        ElectronicState::Collinear(r) => {
            let dft_physics::spin::SpinDensity::Collinear { channels, .. } = &r.density else {
                return Err(ScfError::InvalidInput(
                    "observables.forces.model: spin model mismatch; rerun SCF with the requested model",
                ));
            };
            if r.k_points[0].iter().map(|k| k.point).collect::<Vec<_>>()
                != r.k_points[1].iter().map(|k| k.point).collect::<Vec<_>>()
            {
                return Err(ScfError::InvalidInput(
                    "observables.forces.model: spin k-point mismatch; rerun SCF with the requested model",
                ));
            }
            (
                &r.force_identity,
                OccupationIdentity::new(
                    &[&r.occupations[0], &r.occupations[1]],
                    r.finite_occupations.as_ref(),
                    Some(&r.sector),
                ),
                &channels[0],
                r.k_points[0].iter().map(|k| k.point).collect(),
            )
        }
    };
    let dims: [usize; 3] = density
        .shape()
        .try_into()
        .map_err(|_| ScfError::InvalidInput("observables.forces.model: invalid grid; rerun SCF"))?;
    let mut model = ModelDescription::new(
        state.cell,
        state.atoms,
        state.potentials,
        state.options,
        dims,
        points,
    );
    model.floor = state.density_floor;
    model.symmetry = state.symmetry.to_vec();
    if !identity.matches(&model, &policy) || state.options.grid_override.is_some_and(|g| g != dims)
    {
        return Err(ScfError::InvalidInput(
            "observables.forces.model: model or occupation policy differs from the certified state; use its original model or rerun SCF with the requested model",
        ));
    }

    let (certified, standard, densities, blocks, quantity) = match &state.electronic {
        ElectronicState::Fixed(r) => (
            r.certified
                && r.history
                    .last()
                    .and_then(|h| h.fixed_point_certificate.as_ref())
                    .is_some(),
            r.standard.is_some(),
            vec![&r.density],
            r.k_points
                .iter()
                .map(|k| (k, r.occupations.as_slice()))
                .collect::<Vec<_>>(),
            "E",
        ),
        ElectronicState::FiniteTemperature(r) => (
            r.converged && r.certificate.as_ref().is_some_and(|c| c.passed),
            false,
            vec![&r.state.density],
            r.state
                .k_points
                .iter()
                .enumerate()
                .map(|(k, p)| (p, r.occupations.block(0, k)))
                .collect(),
            "F",
        ),
        ElectronicState::Collinear(r) => {
            let dft_physics::spin::SpinDensity::Collinear { channels, .. } = &r.density else {
                return Err(ScfError::InvalidInput(
                    "forces.spin: collinear density required",
                ));
            };
            let certified = r.converged
                && if r.finite_occupations.is_some() {
                    r.finite_certificate.as_ref().is_some_and(|c| c.passed)
                } else {
                    r.certificate.as_ref().is_some_and(|c| c.passed)
                };
            let blocks = r
                .k_points
                .iter()
                .enumerate()
                .flat_map(|(s, ks)| {
                    ks.iter().enumerate().map(move |(k, p)| {
                        (
                            p,
                            r.finite_occupations
                                .as_ref()
                                .map_or(r.occupations[s].as_slice(), |f| f.block(s, k)),
                        )
                    })
                })
                .collect();
            (
                certified,
                false,
                vec![&channels[0], &channels[1]],
                blocks,
                if r.finite_occupations.is_some() {
                    "F"
                } else {
                    "E"
                },
            )
        }
    };
    if state.options.standard.is_some() || standard {
        return Err(ScfError::InvalidInput(
            "observables.forces: strict convergence profile required; rerun with convergence_profile = strict",
        ));
    }
    if state.options.density_tolerance > 1e-7
        || state.options.energy_tolerance > 1e-8
        || state.options.eigen_tolerance > 1e-8
    {
        return Err(ScfError::InvalidInput(
            "observables.forces: certificate tolerances exceed strict caps; tighten density<=1e-7 and energy/eigen<=1e-8 Ha",
        ));
    }
    if !certified {
        return Err(ScfError::InvalidInput(
            "observables.forces: passing fixed-point certificate required; converge the state",
        ));
    }
    if state.atoms.is_empty()
        || blocks.is_empty()
        || densities[0].is_empty()
        || densities.iter().any(|d| d.shape() != densities[0].shape())
    {
        return Err(ScfError::InvalidInput(
            "observables.forces: incomplete state arrays; use the complete converged electronic state",
        ));
    }
    let mut density = densities[0].clone();
    if densities.len() == 2 {
        density += densities[1];
    }
    let (energy, entropy) = match &state.electronic {
        ElectronicState::Fixed(r) => (&r.energy, None),
        ElectronicState::FiniteTemperature(r) => {
            (&r.state.energy, Some(r.thermodynamics.entropy_term_ha))
        }
        ElectronicState::Collinear(r) => (&r.energy, r.thermodynamics.map(|t| t.entropy_term_ha)),
    };
    // Intentionally exhaustive: adding an energy/audit field must trigger a
    // force-coverage review instead of silently retaining a partial total.
    let crate::energy::EnergyAudit {
        xc_floor: _,
        nlcc,
        kinetic,
        local,
        local_zero_mode: _,
        nonlocal,
        hartree,
        xc,
        ewald,
        band_sum: _,
        double_counting: _,
        band_form: _,
        harris_foulkes: _,
        functional_difference_bound: _,
    } = energy;
    let mut ledger = derivative_ledger();
    for entry in &mut ledger {
        entry.energy_ha = Some(match entry.term_id {
            "kinetic" => *kinetic,
            "local" => *local,
            "nonlocal" => *nonlocal,
            "hartree" => *hartree,
            "xc" => *xc,
            "ewald" => *ewald,
            "entropy" => {
                entry.active = entropy.is_some();
                entry.included_in = "correction";
                -entropy.unwrap_or(0.0)
            }
            "core" => {
                entry.active = nlcc.is_some();
                entry.included_in = "audit_only";
                entry.dependencies.push("xc");
                nlcc.as_ref().map_or(0.0, |c| c.core_potential_integral_ha)
            }
            _ => unreachable!("closed NC ledger"),
        });
    }
    require_complete(&ledger).map_err(ScfError::Numerical)?;
    let volume = state.cell.volume_bohr3();
    let density_basis = state.options.basis(state.cell, [0.; 3], 4.0)?;
    let dims: [usize; 3] = density
        .shape()
        .try_into()
        .map_err(|_| ScfError::InvalidInput("forces.density: invalid shape"))?;
    let dg = FftGridSpec::new(&density_basis, &density_basis, Some(dims))?;
    // One Species entry per input atom prevents grouping from permuting output.
    let positions: Vec<_> = state
        .atoms
        .iter()
        .map(|a| {
            state
                .cell
                .fractional_to_cartesian_bohr(a.fractional_position.map(|x| x.rem_euclid(1.0)))
        })
        .collect();
    let species: Vec<_> = state
        .atoms
        .iter()
        .zip(&positions)
        .map(|(a, p)| {
            let potential = state
                .potentials
                .iter()
                .find(|s| s.label == a.species)
                .ok_or(ScfError::InvalidInput("forces.species: missing potential"))?
                .potential;
            Ok(Species {
                potential,
                positions_bohr: std::slice::from_ref(p),
            })
        })
        .collect::<Result<_, ScfError>>()?;
    let mut fft = FftPlan3d::new(dims)?;
    let mut reciprocal = Array3::zeros(dims);
    fft.to_reciprocal(&density.mapv(|x| Complex64::new(x, 0.0)), &mut reciprocal)?;
    let mut local = terms::local_or_core(
        &density_basis,
        volume,
        &species,
        &dg.gather(&reciprocal)?,
        false,
    )?;
    let core_density =
        dft_physics::density::core_density(&density_basis, &dg, volume, &species, &mut fft)?;
    let mut core = if let Some(core_density) = core_density {
        let mut scratch = Array3::zeros(dims);
        let potential = if densities.len() == 1 {
            state
                .options
                .evaluate_xc(
                    &density,
                    Some(&core_density),
                    &mut scratch,
                    state.cell,
                    &mut fft,
                    state.density_floor,
                )?
                .potential
        } else {
            use dft_physics::spin::{SpinFunctional, SpinXc, SpinXcResponse};
            if state.density_floor.is_some() {
                return Err(ScfError::InvalidInput(
                    "forces.spin: density floor unavailable",
                ));
            }
            let kind = match state.options.functional {
                crate::Functional::LdaPw => SpinFunctional::LdaPw,
                crate::Functional::Pbe => SpinFunctional::Pbe,
                _ => {
                    return Err(ScfError::InvalidInput(
                        "forces.spin: polarized PZ unavailable",
                    ));
                }
            };
            let mut xc = SpinXc::new(kind, state.cell, dims)?;
            let mut out = SpinXcResponse::new(dims);
            xc.evaluate_with_core(
                &[densities[0].clone(), densities[1].clone()],
                Some(&core_density),
                &mut out,
            )?;
            (&out.potentials[0] + &out.potentials[1]) * 0.5
        };
        fft.to_reciprocal(&potential.mapv(|x| Complex64::new(x, 0.0)), &mut reciprocal)?;
        terms::local_or_core(
            &density_basis,
            volume,
            &species,
            &dg.gather(&reciprocal)?,
            true,
        )?
    } else {
        vec![[0.0; 3]; state.atoms.len()]
    };
    let ions: Vec<_> = state
        .atoms
        .iter()
        .zip(&species)
        .map(|(a, s)| ewald::Ion {
            charge: s.potential.z_valence(),
            fractional_position: a.fractional_position,
        })
        .collect();
    let mut ewald = ewald::forces(state.cell, &ions, ewald::EwaldOptions::default())?;
    let mut nonlocal = vec![[0.0; 3]; state.atoms.len()];
    let mut band_passes = 0;
    for (k, occupations) in blocks {
        let basis = state.options.basis(state.cell, k.point.fractional, 1.0)?;
        let operator = NonlocalOperator::new(&basis, volume, &species)?;
        let occupations: Vec<_> = occupations.iter().map(|f| f * k.point.weight).collect();
        let part = terms::nonlocal(&operator, &basis, &k.solution.eigenvectors, &occupations)?;
        band_passes += occupations.len();
        for (sum, part) in nonlocal.iter_mut().zip(part) {
            for axis in 0..3 {
                sum[axis] += part[axis];
            }
        }
    }
    for term in [&mut ewald, &mut local, &mut nonlocal, &mut core] {
        symmetrize(state.cell, term, state.symmetry)?;
    }
    let values: Vec<[f64; 3]> = (0..state.atoms.len())
        .map(|i| std::array::from_fn(|a| ewald[i][a] + local[i][a] + nonlocal[i][a] + core[i][a]))
        .collect();
    if values.iter().flatten().any(|f| !f.is_finite()) {
        return Err(ScfError::NonFinite);
    }
    let net_force: [f64; 3] = std::array::from_fn(|a| values.iter().map(|f| f[a]).sum());
    if net_force.iter().any(|f| !f.is_finite()) {
        return Err(ScfError::NonFinite);
    }
    let maximum_absolute_component = values.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max);
    Ok(ForceResult {
        ledger,
        units: "Ha/Bohr",
        quantity,
        values,
        ewald,
        local,
        nonlocal,
        core,
        net_force,
        maximum_absolute_component,
        elapsed: start.elapsed(),
        band_passes,
    })
}

/// Polar-vector transformation A R A^-1 and the SCF atom permutation. No
/// optimizer constraint or subtraction of a residual net force is performed.
fn symmetrize(
    cell: &Cell,
    values: &mut Vec<[f64; 3]>,
    ops: &[dft_geometry::symmetry::SymmetryOperation],
) -> Result<(), ScfError> {
    if ops.is_empty() {
        return Ok(());
    }
    let mut result = vec![[0.0; 3]; values.len()];
    for op in ops {
        if op.atom_permutation.len() != values.len() {
            return Err(ScfError::InvalidInput(
                "forces.symmetry: invalid permutation",
            ));
        }
        let mut seen = vec![false; values.len()];
        for (i, v) in values.iter().enumerate() {
            let target = op.atom_permutation[i];
            if target >= values.len() || seen[target] {
                return Err(ScfError::InvalidInput(
                    "forces.symmetry: invalid permutation",
                ));
            }
            seen[target] = true;
            let f = cell.cartesian_bohr_to_fractional(*v);
            let rotated = cell.fractional_to_cartesian_bohr(std::array::from_fn(|a| {
                (0..3).map(|b| op.rotation[a][b] as f64 * f[b]).sum()
            }));
            for (a, component) in rotated.iter().enumerate() {
                result[target][a] += component / ops.len() as f64;
            }
        }
    }
    *values = result;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn symmetry_rotates_polar_vectors_and_permutes_atoms() {
        use dft_geometry::symmetry::{SpinAction, SymmetryOperation};
        let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
        let identity = SymmetryOperation {
            rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            translation: [0.; 3],
            atom_permutation: vec![0, 1],
            spin_action: SpinAction::Identity,
            antiunitary: false,
            model_invariance_proof: "test",
        };
        let inversion = SymmetryOperation {
            rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
            atom_permutation: vec![1, 0],
            ..identity.clone()
        };
        let mut values = vec![[1., 2., 3.], [2., 4., 6.]];
        symmetrize(&cell, &mut values, &[identity, inversion]).unwrap();
        assert_eq!(values, vec![[-0.5, -1., -1.5], [0.5, 1., 1.5]]);
    }
}
