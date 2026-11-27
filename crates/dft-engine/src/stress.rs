//! Post-convergence frozen-basis stress and separately labelled cutoff estimate.
use crate::{
    Functional, ScfError,
    forces::{DerivativeEntry, ElectronicState, ForceState},
};
use dft_basis::FftGridSpec;
use dft_physics::{data::*, hamiltonian::Species, stress as terms};
pub type Tensor = terms::Tensor;
/// Standard-profile residual threshold of the cutoff-scan SCF runs, Ha. The
/// scan differences enter the estimate scaled by 2c/(3 V dc), about 1e-2 per Ha
/// for a 30 Ry silicon cell, so a 1e-9 Ha scan energy error moves the estimate
/// by about 1e-11 Ha/Bohr^3, orders below its reported uncertainty.
///
/// ONEPATH-2 (ADR 0043): for a standard caller the scan threshold is
/// `min(SCAN_THRESHOLD_HA, caller threshold)`, so a scan is never looser than
/// before and follows a caller who tightens `scf.standard_threshold_ha`.
pub const SCAN_THRESHOLD_HA: f64 = 1e-9;
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
    /// How the three scan SCF runs were converged (ONEPATH-2):
    /// `standard_then_strict` for a strict caller (standard profile at
    /// `SCAN_THRESHOLD_HA` without verification, a failed scan repeated with the
    /// caller's strict settings; unchanged reference behaviour) or
    /// `standard_verified` for a standard caller (every scan on the standard
    /// path with the fixed-point verification; no strict fallback).
    pub scan_profile: &'static str,
    /// Threshold each scan (centre, +5 %, +10 %) actually converged to, Ha;
    /// None for a scan that fell back to the caller's strict settings.
    pub scan_thresholds_ha: [Option<f64>; 3],
    /// SCF energy resolution entering `uncertainty`, Ha: the caller's
    /// `energy_tolerance`, and for a standard caller at least the loosest scan
    /// threshold actually used.
    pub scan_energy_resolution_ha: f64,
}
#[derive(Debug)]
pub struct StressResult {
    pub verification: &'static str,
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
    /// Band passes of the precondition force pass: 0 unless
    /// `STRESS_FORCE_PASS=on` (GD1; the stress itself never used them).
    pub precondition_force_band_passes: usize,
    /// GD1: wall seconds per stage of the frozen-Miller stress, in execution
    /// order (observation only; per rank).
    pub timings: Vec<(&'static str, f64)>,
}
/// Complete analytic derivative at fixed Miller set. Explicitly does not claim
/// cutoff convergence; `stress` additionally estimates the cutoff correction.
pub fn frozen_stress(state: &StressState<'_>) -> Result<StressResult, ScfError> {
    preconditions(state)?;

    // ADR 0043 (ONEPATH-2): a standard state is accepted exactly when forces
    // accept it, i.e. with a passed fixed-point verification; the shared force
    // path below rejects an unverified standard state.
    let start = std::time::Instant::now();
    // Shared force path enforces the independent setup identity, active ledger,
    // strict tolerances and fresh certificate (or passed verification) before
    // accepting the state.
    // GD1: the preconditions alone unless STRESS_FORCE_PASS asks for the
    // former complete force pass (its forces were never used here).
    let validated = if dft_numerics::switches::on(dft_numerics::switches::Id::STRESS_FORCE_PASS) {
        let f = crate::forces::forces(state)?;
        (
            crate::forces::Validated {
                verification: f.verification,
                ledger: f.ledger,
                quantity: f.quantity,
            },
            f.band_passes,
        )
    } else {
        (crate::forces::validate(state)?, 0)
    };
    let (validated, precondition_band_passes) = validated;
    frozen_pass(state, validated, precondition_band_passes, false, start).map(|r| r.0)
}

/// GD2 (`DERIV_JOINT`): forces and frozen-Miller stress of one state from one
/// pass (shared preconditions, density transform, shells, core density, XC
/// evaluation, core-potential transform, projector tables and k pass). Both
/// results are bit-identical to `forces::forces` and `frozen_stress`. `None`
/// when the joint pass does not apply (`DERIV_NL_GEMM` or `DERIV_SHELLS` off,
/// or `STRESS_FORCE_PASS` on); the caller then runs the separate passes. The
/// force record carries the joint stage clocks; the stress record reports a
/// single `joint_with_forces` stage of zero seconds.
pub fn forces_and_frozen_stress(
    state: &StressState<'_>,
) -> Result<Option<(crate::forces::ForceResult, StressResult)>, ScfError> {
    use dft_numerics::switches::{Id, on};
    if !on(Id::DERIV_JOINT)
        || !on(Id::DERIV_NL_GEMM)
        || !on(Id::DERIV_SHELLS)
        || on(Id::STRESS_FORCE_PASS)
    {
        return Ok(None);
    }
    let start = std::time::Instant::now();
    preconditions(state)?;
    let validated = crate::forces::validate(state)?;
    let (mut stress, forces) = frozen_pass(state, validated, 0, true, start)?;
    stress.timings = vec![("joint_with_forces", 0.)];
    stress.elapsed = std::time::Duration::ZERO;
    Ok(forces.map(|f| (f, stress)))
}

fn preconditions(state: &StressState<'_>) -> Result<(), ScfError> {
    if let ElectronicState::Fixed(r) = &state.electronic
        && !r.orbitals_delivered()
    {
        return Err(ScfError::InvalidInput(
            "stress needs the delivered orbitals: call ScfResult::materialize_orbitals",
        ));
    }
    if let ElectronicState::Reference(r) = &state.electronic
        && !r.orbitals_delivered()
    {
        return Err(ScfError::InvalidInput(
            "stress need the delivered orbitals: call electronic::Run::materialize_orbitals",
        ));
    }
    if state.options.local_mesh != crate::LocalMesh::Baandi
        || state.options.radial_even_rule != crate::RadialEvenRule::Simpson3_8
    {
        return Err(ScfError::InvalidInput(
            "research pseudo conventions do not support analytic derivatives",
        ));
    }
    Ok(())
}

/// The frozen-Miller stress after its preconditions; with `with_forces` also
/// the forces of `crate::forces::forces` from the same intermediate fields.
fn frozen_pass(
    state: &StressState<'_>,
    validated: crate::forces::Validated,
    precondition_band_passes: usize,
    with_forces: bool,
    start: std::time::Instant,
) -> Result<(StressResult, Option<crate::forces::ForceResult>), ScfError> {
    let mut timings = vec![("validate", start.elapsed().as_secs_f64())];
    let mut clock = std::time::Instant::now();
    let mut lap = |name: &'static str, timings: &mut Vec<(&'static str, f64)>| {
        timings.push((name, clock.elapsed().as_secs_f64()));
        clock = std::time::Instant::now();
    };
    let densities: Vec<&Array3<f64>> = state.electronic.densities()?;
    let blocks = state.electronic.blocks();
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
    lap("setup", &mut timings);
    let shells = crate::forces::shell_values(&species, &db, volume, true)?;
    lap("shells", &mut timings);
    let density_g = dg.gather(&reciprocal)?;
    let local =
        terms::local_or_core_cached(&db, volume, &species, &density_g, false, shells.as_ref())?;
    let mut local_forces = if with_forces {
        dft_physics::forces::local_or_core_cached(
            &db,
            volume,
            &species,
            &density_g,
            false,
            shells.as_ref(),
        )?
    } else {
        Vec::new()
    };
    lap("local", &mut timings);
    let hartree = terms::hartree(&density, state.cell, &mut fft)?;
    lap("hartree", &mut timings);
    let core_density = dft_physics::density::core_density_cached(
        &db,
        &dg,
        volume,
        &species,
        &mut fft,
        shells.as_ref(),
    )?;
    let (energy, integral, core_potential, mut xc) = if densities.len() == 1 {
        let mut scratch = Array3::zeros(dims);
        let (result, metric) = state.options.evaluate_xc_for_derivatives(
            &density,
            core_density.as_ref(),
            &mut scratch,
            state.cell,
            &mut fft,
            state.density_floor,
            true,
        )?;
        let metric = if let Some(metric) = metric {
            // GD2: the metric of the same PBE evaluation (bit-identical).
            metric
        } else if state.options.functional == Functional::Pbe {
            let mut combined = density.clone();
            if let Some(c) = &core_density {
                combined += c;
            }
            if let Some(floor) = state.density_floor {
                combined.mapv_inplace(|v| v.max(floor));
            }
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
        if with_forces && core_density.is_some() && state.density_floor.is_some() {
            // The force pass rejects this combination (`forces::forces`).
            return Err(ScfError::InvalidInput(
                "forces.spin: density floor unavailable",
            ));
        }
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
    lap("xc", &mut timings);
    let mut core_forces = vec![[0.0; 3]; if with_forces { state.atoms.len() } else { 0 }];
    let core = if core_density.is_some() {
        fft.to_reciprocal(
            &core_potential.mapv(|v| Complex64::new(v, 0.)),
            &mut reciprocal,
        )?;
        let potential_g = dg.gather(&reciprocal)?;
        if with_forces {
            core_forces = dft_physics::forces::local_or_core_cached(
                &db,
                volume,
                &species,
                &potential_g,
                true,
                shells.as_ref(),
            )?;
        }
        terms::local_or_core_cached(&db, volume, &species, &potential_g, true, shells.as_ref())?
    } else {
        [[0.; 3]; 3]
    };
    lap("core", &mut timings);
    let mut kinetic = [[0.; 3]; 3];
    let mut nonlocal = [[0.; 3]; 3];
    let mut band_passes = 0;
    let gemm_kernel = dft_numerics::switches::on(dft_numerics::switches::Id::DERIV_NL_GEMM);
    let tables = if gemm_kernel {
        crate::forces::projector_tables(state, &species, blocks.len(), true)?
    } else {
        Vec::new()
    };
    let per_atom = crate::forces::tables_per_atom(&species, &tables);
    lap("nonlocal_tables", &mut timings);
    // MPI-DERIV: per-k kinetic and nonlocal tensors distributed over ranks
    // (one gather), accumulated in k order as before; symmetrization follows.
    let natoms = state.atoms.len();
    let mut nonlocal_forces = vec![[0.0; 3]; if with_forces { natoms } else { 0 }];
    let stride = 18 + if with_forces { 3 * natoms } else { 0 };
    crate::forces::accumulate_blocks(
        &blocks,
        stride,
        |(point, orbitals, occupations), out| {
            let basis = state.options.basis(state.cell, point.fractional, 1.)?;
            let f: Vec<_> = occupations.iter().map(|v| v * point.weight).collect();
            let c: Vec<&[Complex64]> = (0..orbitals.cols()).map(|b| orbitals.column(b)).collect();
            let k = dft_physics::kinetic::stress(&basis, &c, &f, volume)?;
            let (nf, n) = if gemm_kernel {
                dft_physics::stress::nonlocal_k_derivatives(
                    &basis, volume, &species, &per_atom, orbitals, &f, true,
                )?
            } else {
                (
                    Vec::new(),
                    terms::nonlocal(&basis, volume, &species, orbitals, &f)?,
                )
            };
            if with_forces && nf.len() != natoms {
                return Err(ScfError::InvalidInput("nonlocal forces: atom count"));
            }
            for (slot, value) in out.iter_mut().zip(
                k.iter()
                    .flatten()
                    .chain(n.iter().flatten())
                    .chain(nf.iter().flatten().take(stride - 18)),
            ) {
                *slot = *value;
            }
            Ok(())
        },
        |part| {
            for (a, b) in kinetic.iter_mut().flatten().zip(&part[..9]) {
                *a += b;
            }
            for (a, b) in nonlocal.iter_mut().flatten().zip(&part[9..18]) {
                *a += b;
            }
            for (sum, part) in nonlocal_forces
                .iter_mut()
                .zip(part[18..].as_chunks::<3>().0)
            {
                for axis in 0..3 {
                    sum[axis] += part[axis];
                }
            }
        },
    )?;
    band_passes += blocks.iter().map(|b| b.2.len()).sum::<usize>();
    lap("kinetic_nonlocal", &mut timings);
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
    lap("ewald", &mut timings);
    let forces = if with_forces {
        let mut ewald_forces = dft_physics::ewald::forces(
            state.cell,
            &ions,
            dft_physics::ewald::EwaldOptions::default(),
        )?;
        for term in [
            &mut ewald_forces,
            &mut local_forces,
            &mut nonlocal_forces,
            &mut core_forces,
        ] {
            crate::forces::symmetrize(state.cell, term, state.symmetry)?;
        }
        lap("forces_ewald_symmetrize", &mut timings);
        Some((ewald_forces, local_forces, nonlocal_forces, core_forces))
    } else {
        None
    };
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
    let forces = match forces {
        Some((ewald, local, nonlocal, core)) => Some(crate::forces::assemble(
            crate::forces::Validated {
                verification: validated.verification,
                ledger: validated.ledger.clone(),
                quantity: validated.quantity,
            },
            ewald,
            local,
            nonlocal,
            core,
            band_passes,
            timings.clone(),
            start,
        )?),
        None => None,
    };
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
    Ok((
        StressResult {
            verification: validated.verification,
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
            precondition_force_band_passes: precondition_band_passes,
            timings,
        },
        forces,
    ))
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
/// Whether a stress evaluation also runs the isotropic cutoff-derivative Pulay
/// assessment (three extra SCF runs at 100 %, 105 % and 110 % cutoff). The
/// assessment is a diagnostic: the stress tensor is the frozen-Miller tensor
/// bit for bit either way. PULAY-OPTIN (user decision 2026-11-26, ADR 0033
/// amendment): off by default; adapters opt in explicitly.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum PulayAssessmentMode {
    /// Frozen-Miller stress only (`StressResult::pulay` is None, no extra SCF).
    #[default]
    Off,
    /// Frozen-Miller stress plus the three-SCF cutoff scan.
    CutoffScan,
}
/// Stress with the assessment selected explicitly: `Off` is `frozen_stress`,
/// `CutoffScan` is `stress` (or `stress_prepared` when `prepared` is given).
pub fn stress_with(
    state: &StressState<'_>,
    prepared: Option<&crate::refscf::symmetry::Prepared>,
    mode: PulayAssessmentMode,
) -> Result<StressResult, ScfError> {
    match (mode, prepared) {
        (PulayAssessmentMode::Off, _) => frozen_stress(state),
        (PulayAssessmentMode::CutoffScan, Some(p)) => stress_prepared(state, p),
        (PulayAssessmentMode::CutoffScan, None) => stress(state),
    }
}
/// GD2: `stress_with` given an already computed frozen-Miller stress of the
/// same state (`frozen`, from the joint force pass), which replaces the
/// frozen pass; `None` is `stress_with`. The Pulay assessment, when
/// selected, runs exactly as in `stress_with`.
pub(crate) fn stress_with_frozen(
    state: &StressState<'_>,
    prepared: Option<&crate::refscf::symmetry::Prepared>,
    mode: PulayAssessmentMode,
    frozen: Option<StressResult>,
) -> Result<StressResult, ScfError> {
    let Some(frozen) = frozen else {
        return stress_with(state, prepared, mode);
    };
    match (mode, prepared) {
        (PulayAssessmentMode::Off, _) => Ok(frozen),
        (PulayAssessmentMode::CutoffScan, Some(p)) => {
            if state.symmetry != p.operations(state.options.crystal_density_group) {
                return Err(ScfError::InvalidInput(
                    "reference stress: derivative and scan symmetry differ",
                ));
            }
            pulay_scan(state, Some(p), frozen)
        }
        (PulayAssessmentMode::CutoffScan, None) => pulay_scan(state, None, frozen),
    }
}
/// Stress with a separately reported isotropic cutoff-derivative Pulay estimate.
/// No anisotropic convergence claim follows from this estimate. Always runs
/// the assessment; see `stress_with` for the explicit switch.
pub fn stress(state: &StressState<'_>) -> Result<StressResult, ScfError> {
    stress_inner(state, None)
}
/// Retain the checked mesh and density projector in reference cutoff scans.
pub fn stress_prepared(
    state: &StressState<'_>,
    prepared: &crate::refscf::symmetry::Prepared,
) -> Result<StressResult, ScfError> {
    if state.symmetry != prepared.operations(state.options.crystal_density_group) {
        return Err(ScfError::InvalidInput(
            "reference stress: derivative and scan symmetry differ",
        ));
    }
    stress_inner(state, Some(prepared))
}
fn stress_inner(
    state: &StressState<'_>,
    prepared: Option<&crate::refscf::symmetry::Prepared>,
) -> Result<StressResult, ScfError> {
    let result = frozen_stress(state)?;
    pulay_scan(state, prepared, result)
}
fn pulay_scan(
    state: &StressState<'_>,
    prepared: Option<&crate::refscf::symmetry::Prepared>,
    mut result: StressResult,
) -> Result<StressResult, ScfError> {
    let start = std::time::Instant::now();
    let c = state.options.cutoff_hartree;
    let steps = [0.05 * c, 0.10 * c];
    let mut energies = [0.; 3];
    let original_dims: [usize; 3] = state.electronic.densities()?[0]
        .shape()
        .try_into()
        .map_err(|_| ScfError::InvalidInput("stress: grid shape"))?;
    let mut high = state.options.clone();
    high.strain_reference = None;
    high.cutoff_hartree = c + steps[1];
    let db = high.basis(state.cell, [0.; 3], 4.)?;
    let required = FftGridSpec::new(&db, &db, None)?.dims();
    let mut scan_dims = std::array::from_fn(|i| original_dims[i].max(required[i]));
    if let Some(prepared) = prepared {
        scan_dims = dft_geometry::symmetry::compatible_fft_grid(
            scan_dims,
            prepared.operations(state.options.crystal_density_group),
        )
        .map_err(|e| ScfError::Numerical(e.to_string()))?;
    }
    let points = state.electronic.points();
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
    // The three scan SCF runs are independent and their results do not depend on
    // the thread budget, so they run concurrently with a third of the budget each
    // (rounded up); the energies are unchanged bit for bit.
    let budget = state.options.threads.max(1);
    let per_scan = budget.div_ceil(3).max(1);
    let caller = state.options.standard;
    let scan_threshold =
        caller.map_or(SCAN_THRESHOLD_HA, |p| p.threshold_ha.min(SCAN_THRESHOLD_HA));
    // Each scan returns its energy and the standard threshold it converged to
    // (None: the strict fallback of a strict caller).
    let run_scan = |delta: f64| -> Result<(f64, Option<f64>), ScfError> {
        let mut options = state.options.clone();
        options.strain_reference = None;
        options.grid_override = Some(scan_dims);
        options.cutoff_hartree = c + delta;
        options.threads = per_scan;
        let mut fast = options.clone();
        fast.standard = Some(crate::StandardProfile {
            threshold_ha: scan_threshold,
        });
        match caller {
            // Strict caller (the reference path): standard profile first; a
            // scan that fails to converge there (stall, exhausted budget) is
            // repeated with the caller's strict settings. Unchanged.
            // ONEPATH-3: spin and Fermi-Dirac problems used to reject the
            // standard profile, so their strict-caller scans always fell back
            // to strict; they go there directly to keep the reference bit for bit.
            None if !matches!(state.electronic, ElectronicState::Fixed(_)) => {
                energy_at_prepared(state, state.cell, options, prepared).map(|e| (e, None))
            }
            None => match energy_at_prepared(state, state.cell, fast, prepared) {
                Ok(e) => Ok((e, Some(scan_threshold))),
                Err(_) => {
                    energy_at_prepared(state, state.cell, options, prepared).map(|e| (e, None))
                }
            },
            // Standard caller (ONEPATH-2): every scan runs the standard path
            // with the fixed-point verification (a scan energy is returned only
            // from a converged run whose verification passed). A scan that
            // fails at the scan threshold is repeated once at the caller's own
            // threshold when that is looser; there is no strict fallback.
            Some(p) => {
                fast.verification = crate::Verification::FixedPoint;
                match energy_at_prepared(state, state.cell, fast, prepared) {
                    Ok(e) => Ok((e, Some(scan_threshold))),
                    Err(first) if p.threshold_ha > scan_threshold => {
                        options.verification = crate::Verification::FixedPoint;
                        energy_at_prepared(state, state.cell, options, prepared)
                            .map(|e| (e, Some(p.threshold_ha)))
                            .map_err(|_| first)
                    }
                    Err(first) => Err(first),
                }
            }
        }
    };
    let mut used = [None; 3];
    if budget >= 3
        && !matches!(state.electronic, ElectronicState::Fixed(r) if r.force_identity.reference_delivery)
        && !matches!(state.electronic, ElectronicState::Reference(_))
    {
        let results: Vec<Result<(f64, Option<f64>), ScfError>> = std::thread::scope(|scope| {
            let handles: Vec<_> = [0., steps[0], steps[1]]
                .into_iter()
                .map(|delta| {
                    let run = &run_scan;
                    scope.spawn(move || run(delta))
                })
                .collect();
            handles
                .into_iter()
                .map(|h| h.join().expect("scan thread"))
                .collect()
        });
        for (i, r) in results.into_iter().enumerate() {
            (energies[i], used[i]) = r?;
        }
    } else {
        for (i, delta) in [0., steps[0], steps[1]].into_iter().enumerate() {
            (energies[i], used[i]) = run_scan(delta)?;
        }
    }
    let estimates = [0, 1].map(|i| {
        2. * c / (3. * state.cell.volume_bohr3()) * (energies[i + 1] - energies[0]) / steps[i]
    });
    // Strict caller: the configured energy tolerance (unchanged). Standard
    // caller: never below the loosest threshold a scan actually converged to.
    let resolution = if caller.is_some() {
        used.iter()
            .flatten()
            .fold(state.options.energy_tolerance, |a, t| a.max(*t))
    } else {
        state.options.energy_tolerance
    };
    let uncertainty = (estimates[0] - estimates[1]).abs()
        + 4. * resolution / (3. * state.cell.volume_bohr3() * 0.05);
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
        scan_profile: if caller.is_some() {
            "standard_verified"
        } else {
            "standard_then_strict"
        },
        scan_thresholds_ha: used,
        scan_energy_resolution_ha: resolution,
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
    energy_at_prepared(state, cell, options, None)
}
fn energy_at_prepared(
    state: &StressState<'_>,
    cell: &dft_geometry::Cell,
    options: crate::ScfOptions,
    prepared: Option<&crate::refscf::symmetry::Prepared>,
) -> Result<f64, ScfError> {
    #[cfg(feature = "strict-oracle")]
    use crate::{ScfProblem, collinear::CollinearProblem, finite_temperature::FiniteTProblem};
    match &state.electronic {
        ElectronicState::Fixed(r) if r.force_identity.reference_delivery => {
            let points: Vec<_> = r.k_points.iter().map(|k| k.point).collect();
            let mut profile = crate::refscf::profile::Profile::from_standard_ha(
                options
                    .standard
                    .ok_or(ScfError::InvalidInput(
                        "reference stress scan requires standard profile",
                    ))?
                    .threshold_ha,
            )?;
            // Preserve the caller's reference mixer for independent cutoff
            // scans as well as its convergence threshold.
            match options.mixing {
                crate::Mixing::Pulay { beta, history } => {
                    profile.beta = beta;
                    profile.history = history;
                }
                crate::Mixing::Linear { beta } => {
                    profile.beta = beta;
                    profile.history = 0;
                }
                _ => {}
            }
            let input = crate::refstart::StartInput {
                cell,
                atoms: state.atoms,
                potentials: state.potentials,
                options: &options,
            };
            let run = if let Some(prepared) = prepared {
                crate::refscf::k_driver::run_prepared(input, profile, prepared)
            } else {
                if !state.symmetry.is_empty() {
                    return Err(ScfError::InvalidInput(
                        "reference symmetric cutoff scan requires stress_prepared",
                    ));
                }
                crate::refscf::k_driver::run(input, profile, &points)
            }?;
            if !run.outcome.converged {
                return Err(ScfError::InvalidInput(
                    "reference stress scan did not converge",
                ));
            }
            Ok(run.outcome.result.total_hartree())
        }
        // Reference spin/FD scans: same occupation policy (theta, constraint,
        // band cap or fixed sector), reference mixer and threshold; the
        // variational energy is F (E for a fixed sector, where TS = 0).
        ElectronicState::Reference(r) => {
            let points = state.electronic.points();
            let mut profile = crate::refscf::profile::Profile::from_standard_ha(
                options
                    .standard
                    .ok_or(ScfError::InvalidInput(
                        "reference stress scan requires standard profile",
                    ))?
                    .threshold_ha,
            )?;
            match options.mixing {
                crate::Mixing::Pulay { beta, history } => {
                    profile.beta = beta;
                    profile.history = history;
                }
                crate::Mixing::Linear { beta } => {
                    profile.beta = beta;
                    profile.history = 0;
                }
                _ => {}
            }
            let input = crate::refstart::StartInput {
                cell,
                atoms: state.atoms,
                potentials: state.potentials,
                options: &options,
            };
            let run = if let Some(prepared) = prepared {
                crate::refscf::electronic::run_prepared(input, profile, prepared, r.policy.clone())
            } else {
                if !state.symmetry.is_empty() {
                    return Err(ScfError::InvalidInput(
                        "reference symmetric cutoff scan requires stress_prepared",
                    ));
                }
                crate::refscf::electronic::run(input, profile, &points, r.policy.clone())
            }?;
            if !run.converged
                || (options.verification == crate::Verification::FixedPoint
                    && !run.verification.as_ref().is_some_and(|v| v.passed()))
            {
                return Err(ScfError::InvalidInput(
                    "reference stress scan did not converge",
                ));
            }
            Ok(run.thermodynamics.free_energy_ha)
        }
        #[cfg(feature = "strict-oracle")]
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
        #[cfg(feature = "strict-oracle")]
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
        #[cfg(feature = "strict-oracle")]
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
        #[cfg(not(feature = "strict-oracle"))]
        _ => Err(ScfError::InvalidInput(
            "stress scan requires strict-oracle for legacy results",
        )),
    }
}
