//! Post-convergence frozen-basis stress and separately labelled cutoff estimate.
use crate::{
    Functional, ScfError,
    forces::{DerivativeEntry, ElectronicState, ForceState},
};
use dft_basis::FftGridSpec;
use dft_physics::{data::*, hamiltonian::Species, stress as terms};
pub type Tensor = terms::Tensor;
type Blocks<'a> = Vec<(&'a crate::KPointSolution, &'a [f64])>;
pub type StressState<'a> = ForceState<'a>;
#[derive(Debug)]
pub struct PulayAssessment {
    pub method: &'static str,
    /// Estimated fixed-cutoff minus frozen-set isotropic stress, Ha/Bohr^3.
    pub correction: f64,
    pub uncertainty: f64,
    pub cutoff_steps_ha: [f64; 2],
    pub warning_threshold: f64,
    pub warning: bool,
    pub shells_resolved: bool,
    pub elapsed: std::time::Duration,
}
#[derive(Debug)]
pub struct StressResult {
    pub units: &'static str,
    pub quantity: &'static str,
    pub basis_protocol: &'static str,
    pub cutoff_ha: f64,
    /// Empty: the production derivative is analytic, with no strain samples.
    pub strain_steps: &'static [f64],
    pub values: Tensor,
    pub symmetric: Tensor,
    pub antisymmetric: Tensor,
    pub antisymmetric_norm: f64,
    pub pressure: f64,
    pub terms: Vec<(&'static str, Tensor)>,
    pub ledger: Vec<DerivativeEntry>,
    pub pulay: Option<PulayAssessment>,
    pub elapsed: std::time::Duration,
    pub band_passes: usize,
    pub precondition_force_band_passes: usize,
}
/// Complete analytic derivative at fixed Miller set. Explicitly does not claim
/// cutoff convergence; `stress` additionally estimates the cutoff correction.
pub fn frozen_stress(state: &StressState<'_>) -> Result<StressResult, ScfError> {
    if state.options.local_mesh != crate::LocalMesh::Baandi
        || state.options.radial_even_rule != crate::RadialEvenRule::Simpson3_8
    {
        return Err(ScfError::InvalidInput(
            "research pseudo conventions do not support analytic derivatives",
        ));
    }

    let start = std::time::Instant::now();
    // Shared force path enforces the independent setup identity, active ledger,
    // strict tolerances and fresh certificate before accepting the state.
    let validated = crate::forces::forces(state)?;
    let (densities, blocks): (Vec<&Array3<f64>>, Blocks<'_>) = match &state.electronic {
        ElectronicState::Fixed(r) => (
            vec![&r.density],
            r.k_points
                .iter()
                .map(|k| (k, r.occupations.as_slice()))
                .collect(),
        ),
        ElectronicState::FiniteTemperature(r) => (
            vec![&r.state.density],
            r.state
                .k_points
                .iter()
                .enumerate()
                .map(|(k, p)| (p, r.occupations.block(0, k)))
                .collect(),
        ),
        ElectronicState::Collinear(r) => {
            let dft_physics::spin::SpinDensity::Collinear { channels, .. } = &r.density else {
                return Err(ScfError::InvalidInput("stress: collinear density required"));
            };
            (
                channels.iter().collect(),
                r.k_points
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
                    .collect(),
            )
        }
    };
    let volume = state.cell.volume_bohr3();
    let dims: [usize; 3] = densities[0]
        .shape()
        .try_into()
        .map_err(|_| ScfError::InvalidInput("stress: invalid density shape"))?;
    let density = if densities.len() == 2 {
        densities[0] + densities[1]
    } else {
        densities[0].clone()
    };
    let positions: Vec<_> = state
        .atoms
        .iter()
        .map(|a| {
            state
                .cell
                .fractional_to_cartesian_bohr(a.fractional_position.map(|v| v.rem_euclid(1.)))
        })
        .collect();
    let species: Vec<_> = state
        .atoms
        .iter()
        .zip(&positions)
        .map(|(a, p)| {
            Ok(Species {
                potential: state
                    .potentials
                    .iter()
                    .find(|s| s.label == a.species)
                    .ok_or(ScfError::InvalidInput("stress: missing species"))?
                    .potential,
                positions_bohr: std::slice::from_ref(p),
            })
        })
        .collect::<Result<_, ScfError>>()?;
    let db = state.options.basis(state.cell, [0.; 3], 4.)?;
    let dg = FftGridSpec::new(&db, &db, Some(dims))?;
    let mut fft = FftPlan3d::new(dims)?;
    let mut reciprocal = Array3::zeros(dims);
    fft.to_reciprocal(&density.mapv(|v| Complex64::new(v, 0.)), &mut reciprocal)?;
    let local = terms::local_or_core(&db, volume, &species, &dg.gather(&reciprocal)?, false)?;
    let hartree = terms::hartree(&density, state.cell, &mut fft)?;
    let core_density = dft_physics::density::core_density(&db, &dg, volume, &species, &mut fft)?;
    let (energy, integral, core_potential, mut xc) = if densities.len() == 1 {
        let mut scratch = Array3::zeros(dims);
        let result = state.options.evaluate_xc(
            &density,
            core_density.as_ref(),
            &mut scratch,
            state.cell,
            &mut fft,
            state.density_floor,
        )?;
        let mut combined = density.clone();
        if let Some(c) = &core_density {
            combined += c;
        }
        if let Some(floor) = state.density_floor {
            combined.mapv_inplace(|v| v.max(floor));
        }
        let metric = if state.options.functional == Functional::Pbe {
            if state.options.pw92_prefactor == crate::Pw92Prefactor::Original {
                dft_physics::gga::stress_metric_pw92(
                    state.options.pw92_prefactor,
                    &combined,
                    state.cell,
                    &mut fft,
                )?
            } else {
                dft_physics::gga::stress_metric(&combined, state.cell, &mut fft)?
            }
        } else {
            [[0.; 3]; 3]
        };
        (
            result.energy.0,
            result.potential_integral.0,
            result.potential,
            metric,
        )
    } else {
        use dft_physics::spin::{SpinFunctional, SpinXc, SpinXcResponse};
        let kind = match state.options.functional {
            Functional::LdaPw => SpinFunctional::LdaPw,
            Functional::Pbe => SpinFunctional::Pbe,
            _ => return Err(ScfError::InvalidInput("stress: polarized PZ unavailable")),
        };
        let mut adapter = SpinXc::new(kind, state.cell, dims)?;
        let mut result = SpinXcResponse::new(dims);
        adapter.evaluate_with_core(
            &[densities[0].clone(), densities[1].clone()],
            core_density.as_ref(),
            &mut result,
        )?;
        (
            result.energy_ha,
            result.valence_potential_integral_ha,
            (&result.potentials[0] + &result.potentials[1]) * 0.5,
            adapter.stress_metric()?,
        )
    };
    for (a, row) in xc.iter_mut().enumerate() {
        row[a] += (energy - integral) / volume;
    }
    let core = if core_density.is_some() {
        fft.to_reciprocal(
            &core_potential.mapv(|v| Complex64::new(v, 0.)),
            &mut reciprocal,
        )?;
        terms::local_or_core(&db, volume, &species, &dg.gather(&reciprocal)?, true)?
    } else {
        [[0.; 3]; 3]
    };
    let mut kinetic = [[0.; 3]; 3];
    let mut nonlocal = [[0.; 3]; 3];
    let mut band_passes = 0;
    for (k, occupations) in blocks {
        let basis = state.options.basis(state.cell, k.point.fractional, 1.)?;
        let f: Vec<_> = occupations.iter().map(|v| v * k.point.weight).collect();
        let c: Vec<_> = (0..k.solution.eigenvectors.cols())
            .map(|b| k.solution.eigenvectors.column(b).to_vec())
            .collect();
        add(
            &mut kinetic,
            &dft_physics::kinetic::stress(&basis, &c, &f, volume)?,
        );
        add(
            &mut nonlocal,
            &terms::nonlocal(&basis, volume, &species, &k.solution.eigenvectors, &f)?,
        );
        band_passes += f.len();
    }
    let ions: Vec<_> = state
        .atoms
        .iter()
        .zip(&species)
        .map(|(a, s)| dft_physics::ewald::Ion {
            charge: s.potential.z_valence(),
            fractional_position: a.fractional_position,
        })
        .collect();
    let ewald = dft_physics::ewald::stress(state.cell, &ions, Default::default())?;
    let mut parts = vec![
        ("kinetic", kinetic),
        ("local", local),
        ("nonlocal", nonlocal),
        ("hartree", hartree),
        ("xc", xc),
        ("core", core),
        ("ewald", ewald),
        ("entropy", [[0.; 3]; 3]),
    ];
    let mut raw = [[0.; 3]; 3];
    for (_, part) in &parts {
        add(&mut raw, part);
    }
    let antisymmetric =
        std::array::from_fn(|a| std::array::from_fn(|b| 0.5 * (raw[a][b] - raw[b][a])));
    let antisymmetric_norm = norm(&antisymmetric);
    if !antisymmetric_norm.is_finite() || antisymmetric_norm >= 1e-10 {
        return Err(ScfError::Numerical(format!(
            "stress: antisymmetric norm {antisymmetric_norm:e} Ha/Bohr^3 exceeds 1e-10"
        )));
    }
    let mut values = [[0.; 3]; 3];
    for (_, part) in &mut parts {
        *part = symmetrize(state.cell, part, state.symmetry);
        let sym = std::array::from_fn(|a| std::array::from_fn(|b| 0.5 * (part[a][b] + part[b][a])));
        *part = sym;
        add(&mut values, part);
    }
    if values.iter().flatten().any(|v| !v.is_finite()) {
        return Err(ScfError::NonFinite);
    }
    let mut ledger = validated.ledger;
    for entry in &mut ledger {
        entry.stress_available = true;
        entry.owner = match entry.term_id {
            "local" | "nonlocal" | "core" | "hartree" => "physics::stress",
            "xc" => "physics::xc/gga/spin",
            _ => entry.owner,
        };
        entry.reason = if entry.term_id == "entropy" {
            "stationary F at fixed theta,N,M: no explicit entropy strain term"
        } else {
            "complete NC frozen-Miller strain derivative"
        };
    }
    Ok(StressResult {
        units: "Ha/Bohr^3",
        quantity: validated.quantity,
        basis_protocol: "frozen_miller_set",
        cutoff_ha: state.options.cutoff_hartree,
        strain_steps: &[],
        values,
        symmetric: values,
        antisymmetric,
        antisymmetric_norm,
        pressure: -(0..3).map(|a| values[a][a]).sum::<f64>() / 3.,
        terms: parts,
        ledger,
        pulay: None,
        elapsed: start.elapsed(),
        band_passes,
        precondition_force_band_passes: validated.band_passes,
    })
}
fn add(a: &mut Tensor, b: &Tensor) {
    for (a, b) in a.iter_mut().flatten().zip(b.iter().flatten()) {
        *a += b;
    }
}
fn norm(a: &Tensor) -> f64 {
    a.iter().flatten().map(|v| v * v).sum::<f64>().sqrt()
}
fn symmetrize(
    cell: &dft_geometry::Cell,
    t: &Tensor,
    ops: &[dft_geometry::symmetry::SymmetryOperation],
) -> Tensor {
    if ops.is_empty() {
        return *t;
    }
    let mut result = [[0.; 3]; 3];
    for op in ops {
        let columns: [[f64; 3]; 3] = std::array::from_fn(|j| {
            let f = cell.cartesian_bohr_to_fractional(std::array::from_fn(|i| {
                if i == j { 1. } else { 0. }
            }));
            cell.fractional_to_cartesian_bohr(std::array::from_fn(|i| {
                (0..3).map(|k| op.rotation[i][k] as f64 * f[k]).sum()
            }))
        });
        for (a, row) in result.iter_mut().enumerate() {
            for (b, v) in row.iter_mut().enumerate() {
                *v += (0..3)
                    .flat_map(|i| (0..3).map(move |j| columns[i][a] * t[i][j] * columns[j][b]))
                    .sum::<f64>()
                    / ops.len() as f64;
            }
        }
    }
    result
}
/// Stress with a separately reported isotropic cutoff-derivative Pulay estimate.
/// No anisotropic convergence claim follows from this estimate.
pub fn stress(state: &StressState<'_>) -> Result<StressResult, ScfError> {
    let mut result = frozen_stress(state)?;
    let start = std::time::Instant::now();
    let c = state.options.cutoff_hartree;
    let steps = [0.05 * c, 0.10 * c];
    let mut energies = [0.; 3];
    let original_dims: [usize; 3] = match &state.electronic {
        ElectronicState::Fixed(r) => r.density.shape(),
        ElectronicState::FiniteTemperature(r) => r.state.density.shape(),
        ElectronicState::Collinear(r) => match &r.density {
            dft_physics::spin::SpinDensity::Collinear { channels, .. } => channels[0].shape(),
            _ => unreachable!(),
        },
    }
    .try_into()
    .map_err(|_| ScfError::InvalidInput("stress: grid shape"))?;
    let mut high = state.options.clone();
    high.strain_reference = None;
    high.cutoff_hartree = c + steps[1];
    let db = high.basis(state.cell, [0.; 3], 4.)?;
    let required = FftGridSpec::new(&db, &db, None)?.dims();
    let scan_dims = std::array::from_fn(|i| original_dims[i].max(required[i]));
    let points: Vec<_> = match &state.electronic {
        ElectronicState::Fixed(r) => r.k_points.iter().map(|k| k.point).collect(),
        ElectronicState::FiniteTemperature(r) => r.state.k_points.iter().map(|k| k.point).collect(),
        ElectronicState::Collinear(r) => r.k_points[0].iter().map(|k| k.point).collect(),
    };
    let mut shells_resolved = false;
    for point in points {
        shells_resolved |= high
            .basis(state.cell, point.fractional, 1.)?
            .vectors()
            .len()
            > state
                .options
                .basis(state.cell, point.fractional, 1.)?
                .vectors()
                .len();
    }
    // Recompute the center too: one grid adequate for the largest cutoff
    // isolates basis effects from FFT-grid changes. Production SCF is unchanged.
    for (i, delta) in [0., steps[0], steps[1]].into_iter().enumerate() {
        let mut options = state.options.clone();
        options.strain_reference = None;
        options.grid_override = Some(scan_dims);
        options.cutoff_hartree = c + delta;
        energies[i] = energy_at(state, state.cell, options)?;
    }
    let estimates = [0, 1].map(|i| {
        2. * c / (3. * state.cell.volume_bohr3()) * (energies[i + 1] - energies[0]) / steps[i]
    });
    let uncertainty = (estimates[0] - estimates[1]).abs()
        + 4. * state.options.energy_tolerance / (3. * state.cell.volume_bohr3() * 0.05);
    let threshold = 1e-5;
    result.pulay = Some(PulayAssessment {
        method: "isotropic_cutoff_derivative_5_10_percent",
        correction: estimates[0],
        uncertainty,
        cutoff_steps_ha: steps,
        warning_threshold: threshold,
        warning: !shells_resolved || estimates[0].abs() + uncertainty > threshold,
        shells_resolved,
        elapsed: start.elapsed(),
    });
    Ok(result)
}
/// Independent SCF energy for strain/cutoff scans. Preserves N, theta, spin
/// sector, fractional k and atoms; the caller specifies the basis protocol.
pub fn energy_at(
    state: &StressState<'_>,
    cell: &dft_geometry::Cell,
    options: crate::ScfOptions,
) -> Result<f64, ScfError> {
    use crate::{ScfProblem, collinear::CollinearProblem, finite_temperature::FiniteTProblem};
    match &state.electronic {
        ElectronicState::Fixed(r) => {
            let points: Vec<_> = r.k_points.iter().map(|k| k.point).collect();
            let mut p = if points.len() == 1 && points[0].fractional == [0.; 3] {
                ScfProblem::new(cell, state.atoms, state.potentials, options)?
            } else {
                ScfProblem::with_kpoints(
                    cell,
                    state.atoms,
                    state.potentials,
                    options,
                    &points,
                    crate::OccupationPolicy::Insulator {
                        minimum_gap_hartree: 1e-8,
                    },
                )?
            };
            if let Some(f) = state.density_floor {
                p.set_density_floor(f)?;
            }
            Ok(p.run()?.total_hartree())
        }
        ElectronicState::FiniteTemperature(r) => {
            let points: Vec<_> = r.state.k_points.iter().map(|k| k.point).collect();
            let r = FiniteTProblem::new(
                cell,
                state.atoms,
                state.potentials,
                options,
                &points,
                r.occupations.theta_ha(),
                r.occupations.bands,
            )?
            .run()?;
            if !r.converged {
                return Err(ScfError::InvalidInput(
                    "stress scan: finite-T certificate failed",
                ));
            }
            Ok(r.thermodynamics.free_energy_ha)
        }
        ElectronicState::Collinear(r) => {
            let points: Vec<_> = r.k_points[0].iter().map(|k| k.point).collect();
            let p = if let Some(f) = &r.finite_occupations {
                CollinearProblem::new_fermi_dirac(
                    cell,
                    state.atoms,
                    state.potentials,
                    options,
                    &points,
                    r.sector.clone(),
                    f.theta_ha(),
                    f.constraint(),
                    f.bands,
                )?
            } else {
                CollinearProblem::new(
                    cell,
                    state.atoms,
                    state.potentials,
                    options,
                    &points,
                    r.sector.clone(),
                )?
            };
            let r = p.run()?;
            if !r.converged {
                return Err(ScfError::InvalidInput(
                    "stress scan: spin certificate failed",
                ));
            }
            Ok(r.thermodynamics
                .map_or(r.energy.total_hartree(), |t| t.free_energy_ha))
        }
    }
}
