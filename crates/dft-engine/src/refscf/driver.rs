//! Independent final-cutoff loop. The only reused SCF operation is one map.
//! TODO(issue #25 patch 2): extract the fixed-point verification contract from
//! session.rs without its scheduling, predictor and reuse switches. Until then
//! a requested verification fails closed at the reference stop: the output is
//! explicitly unverified and `ScfOutcome.converged` is false. `verify=off`
//! exposes the estimator stop for patch-1 trajectory comparisons.
use super::{
    broyden::{Broyden, Diagnostics},
    density::DensityPlan,
    invalid,
    metric::{Metric, Serial},
    numerical,
    profile::{Profile, next_ethr_ry},
    record::Record,
};
use crate::scf::*;
use dft_mpi::Collective;

pub struct Run {
    pub outcome: ScfOutcome,
    pub records: Vec<Record>,
    pub reference_stopped: bool,
    pub exit_reason: String,
}
pub fn run(problem: ScfProblem, profile: Profile) -> Result<Run, ScfError> {
    let profile = profile.validate()?;
    if let Some(comm) = crate::parallel::SharedComm::active() {
        drive(problem, profile, &comm)
    } else {
        drive(problem, profile, &Serial)
    }
}
fn drive<K: Collective>(
    mut problem: ScfProblem,
    profile: Profile,
    comm: &K,
) -> Result<Run, ScfError> {
    let start = Instant::now();
    if problem.states.len() != 1
        || problem.states[0].point.fractional != [0.0; 3]
        || problem.options.gamma_real != crate::GammaReal::On
        || problem.options.standard.is_none()
    {
        return Err(ScfError::InvalidInput(
            "reference SCF requires real Gamma and the standard accuracy class",
        ));
    }
    if problem.start_density.is_none() {
        return Err(ScfError::InvalidInput(
            "reference SCF requires an atomic valence superposition (PP_RHOATOM); no uniform fallback",
        ));
    }
    problem.options.standard = Some(StandardProfile {
        threshold_ha: profile.threshold_ry / 2.0,
    });
    let volume = problem.cell.volume_bohr3();
    let dims = problem.grid.dims();
    let mut plan = DensityPlan::new(
        comm,
        &problem.cell,
        dims,
        4.0 * problem.options.cutoff_hartree,
        4.0 * problem.options.cutoff_hartree,
    )?;
    let mut metric = Metric::new(volume, &plan.g2, &plan.smooth, comm.size(), 350)?;
    let mut mixer = Broyden::new(profile, &plan.g2, &plan.smooth, volume, problem.electrons)?;
    let mut input = problem.starting_density();
    let mut input_planes = vec![0.; plan.plane_len()];
    let mut output_planes = vec![0.; plan.plane_len()];
    let mut next_planes = vec![0.; plan.plane_len()];
    let mut ionic_planes = vec![0.; plan.plane_len()];
    plan.planes_from_full(&input, &mut input_planes)?;
    plan.planes_from_full(&problem.ionic_real, &mut ionic_planes)?;
    if comm.size() > 1 {
        if !problem.owned_density_eligible() {
            return Err(ScfError::InvalidInput(
                "reference SCF MPI requires the owned-density map kernel (LDA PZ/PW, no NLCC/floor/symmetry/snapshot; M2 grid/local fields/local density enabled); no replicated fallback",
            ));
        }
        if problem.grid_plan()?.1.first_axis() != plan.first_axis() {
            return Err(ScfError::InvalidInput(
                "reference/engine owned FFT layouts differ",
            ));
        }
        problem.owned.density_mode = true;
        problem.owned.input = Some(crate::grid_dist::Owned {
            epoch: problem.owned.next_epoch(),
            values: input_planes.clone(),
        });
        input = Array3::zeros((0, 0, 0));
    }
    let mut x = vec![Complex64::default(); plan.miller.len()];
    let mut y = x.clone();
    let mut f = x.clone();
    let mut next = x.clone();
    plan.forward(comm, &input_planes, &mut x)?;
    let mut correction_density = if comm.size() == 1 {
        Array3::zeros(dims)
    } else {
        Array3::zeros((0, 0, 0))
    };
    let mut history = Vec::with_capacity(problem.options.max_iterations);
    let mut records = Vec::with_capacity(problem.options.max_iterations);
    let mut previous: Option<Vec<KPointSolution>> = None;
    let mut ethr = 1e-2;
    let mut last_dr2 = f64::INFINITY;
    for iteration in 1..=problem.options.max_iterations {
        if iteration > 1 {
            ethr = next_ethr_ry(
                if iteration == 2 { 1e-2 } else { ethr },
                last_dr2,
                problem.electrons,
            );
        }
        let mut redo = false;
        let (mut dav, mut columns) = (0, 0);
        let (mut map_s, mut transform_s, mut metric_s) = (0., 0., 0.);
        let mut phases = PhaseTimings::default();
        let (mut step, dr2) = loop {
            let clock = Instant::now();
            let (warm, lazy) = split_warm(previous.as_deref());
            let step = problem.step_session(&input, warm.as_deref(), &lazy, ethr / 2.0, false)?;
            map_s += clock.elapsed().as_secs_f64();
            dav += step
                .k_points
                .iter()
                .map(|k| k.solution.iterations)
                .sum::<usize>();
            columns += step
                .k_points
                .iter()
                .map(|k| k.solution.operator_applications)
                .sum::<usize>();
            phases.potential += step.timings.potential;
            phases.diagonalization += step.timings.diagonalization;
            phases.density += step.timings.density;
            phases.energy += step.timings.energy;
            let clock = Instant::now();
            if let Some(o) = &step.owned_density {
                output_planes.copy_from_slice(&o.values);
            } else {
                plan.planes_from_full(&step.density, &mut output_planes)?;
            }
            plan.forward(comm, &output_planes, &mut y)?;
            for ((v, y), x) in f.iter_mut().zip(&y).zip(&x) {
                *v = y - x;
            }
            transform_s += clock.elapsed().as_secs_f64();
            let clock = Instant::now();
            let dr2 = metric.norm(comm, &f)?;
            metric_s += clock.elapsed().as_secs_f64();
            if iteration == 1 && !redo && dr2 < ethr * problem.electrons.max(1.0) {
                ethr = next_ethr_ry(ethr, dr2, problem.electrons);
                redo = true;
                previous = Some(step.k_points);
                continue;
            }
            break (step, dr2);
        };
        last_dr2 = dr2;
        let stopped = dr2 < profile.threshold_ry;
        let (residual, output_sum, input_sum, min_density) =
            problem.step_metrics(&input, &step.density, step.owned_density.as_ref())?;
        let mut diagnostics = Diagnostics {
            iteration,
            occupancy: mixer.occupancy(),
            ..Default::default()
        };
        let (mut mixing_s, mut correction_s, mut descf) = (0., 0., 0.);
        let mut variational = step.energy.band_form;
        let mut failure = None;
        if !stopped {
            let clock = Instant::now();
            match mixer.update(comm, &mut metric, &x, &f, &mut next) {
                Ok(d) => diagnostics = d,
                Err(e) => failure = Some(e.to_string()),
            }
            if failure.is_none() {
                plan.inverse(comm, &next, &mut next_planes)?;
            }
            mixing_s = clock.elapsed().as_secs_f64();
            phases.mixing += std::time::Duration::from_secs_f64(mixing_s);
            if failure.is_none() {
                let clock = Instant::now();
                let (h, xc) = if comm.size() == 1 {
                    plan.write_serial_full(&next_planes, &mut correction_density)?;
                    let next_full = &correction_density;
                    let (h, xc, _) = problem.input_fields(next_full)?;
                    descf = -volume / next_full.len() as f64
                        * next_full
                            .iter()
                            .zip(&step.density)
                            .zip(h.potential.iter().zip(&xc.potential))
                            .map(|((n, y), (h, x))| (n - y) * (h.re + x))
                            .sum::<f64>();
                    (h.energy.0, xc.energy.0)
                } else {
                    // Existing plane-local potential kernel, no full density gather.
                    let grid = problem.grid_dist.as_mut().expect("owned map built grid");
                    let (h, xc, _) = grid.potential_field(
                        comm,
                        crate::grid_dist::Field::Planes(&next_planes),
                        problem.ionic_real.as_slice().expect("standard"),
                        None,
                    )?;
                    let local = -volume / dims.iter().product::<usize>() as f64
                        * next_planes
                            .iter()
                            .zip(&output_planes)
                            .zip(grid.effective_planes().iter().zip(&ionic_planes))
                            .map(|((n, y), (v, ion))| (n - y) * (v - ion))
                            .sum::<f64>();
                    let mut total = [0.];
                    comm.tree_sum_f64(&[local], &mut total).map_err(numerical)?;
                    descf = total[0];
                    (h, xc)
                };
                variational = step.energy.band_sum - step.energy.double_counting
                    + h
                    + xc
                    + step.energy.ewald
                    + descf;
                correction_s = clock.elapsed().as_secs_f64();
            }
        }
        let energy_change = history
            .last()
            .map(|r: &IterationRecord| step.energy.total_hartree() - r.energy.total_hartree());
        history.push(IterationRecord {
            iteration,
            energy: step.energy.clone(),
            energy_change,
            standard_error_ha: Some(dr2 / 2.),
            standard_redo: redo,
            density_residual: residual,
            max_eigen_residual: step
                .k_points
                .iter()
                .flat_map(|k| &k.solution.residuals)
                .copied()
                .fold(0., f64::max),
            requested_eigen_tolerance: ethr / 2.,
            davidson_iterations: dav,
            operator_applications: columns,
            davidson_iterations_per_k: vec![dav],
            operator_applications_per_k: vec![columns],
            charge_error: output_sum * volume / dims.iter().product::<usize>() as f64
                - problem.electrons,
            input_charge_error: input_sum * volume / dims.iter().product::<usize>() as f64
                - problem.electrons,
            min_density,
            mixing_fraction: 1.0,
            mixer_reset: false,
            density_contraction: None,
            mixing_history_used: Some(diagnostics.occupancy),
            mixing_resets: 0,
            mixing_raw_recovery: false,
            fixed_point_certificate: None,
            timings: phases,
            advance: None,
            pulay: None,
        });
        records.push(Record {
            iteration,
            dr2_ry: dr2,
            ethr_ry: ethr,
            first_iteration_redo: redo,
            beta: profile.beta,
            mixing: diagnostics,
            davidson_iterations: dav,
            h_columns: columns,
            map_s,
            density_transform_s: transform_s,
            metric_s,
            mixing_s,
            correction_s,
            descf_ha: descf,
            variational_energy_ha: variational,
        });
        if stopped || failure.is_some() || iteration == problem.options.max_iterations {
            let verified_off = problem.options.verification == Verification::Off;
            let exit_reason = failure.unwrap_or_else(|| if stopped && !verified_off {
                "reference stop reached; fixed-point verification unavailable in patch 1 (TODO issue #25 patch 2); result is NOT verified".into()
            } else if stopped { "reference estimator stop reached; verification explicitly off".into() }
            else { "reference SCF maximum iterations reached".into() });
            if let Some(o) = &step.owned_density {
                step.density = problem.materialize_owned(o)?;
            }
            problem.materialize_orbitals(&mut step.k_points)?;
            let first = &step.k_points[0];
            let charge = step.density.sum() * volume / step.density.len() as f64;
            let result = ScfResult {
                atomic_start: problem.atomic_start.clone(),
                compat_mode: problem.options.compat_mode,
                crystal_density_group: problem.options.crystal_density_group,
                standard_stall_guard: false,
                standard_density_forcing: false,
                pw92_prefactor: problem.options.pw92_prefactor,
                local_mesh: problem.options.local_mesh,
                radial_even_rule: problem.options.radial_even_rule,
                gamma_real_used: problem.states[0].gamma_real.gamma_real_used,
                gamma_real_fallback_reason: problem.states[0].gamma_real.fallback_reason.clone(),
                force_identity: problem.force_model.with_policy(
                    crate::force_identity::OccupationIdentity::new(
                        &[&problem.occupations],
                        None,
                        None,
                    ),
                ),
                standard: problem.options.standard,
                certified: false,
                verification: (!verified_off).then(VerificationReport::default),
                spectral_completeness: None,
                eigenvalues: first.solution.eigenvalues.clone(),
                orbitals: first.solution.eigenvectors.clone(),
                occupations: problem.occupations.clone(),
                gap_hartree: problem.check_gap(&step.k_points)?,
                k_points: step.k_points,
                energy: step.energy,
                energy_identity: None,
                density: step.density,
                charge,
                history,
                ewald_breakdown: problem.ewald,
                setup_time: problem.setup_time,
                total_time: start.elapsed() + problem.setup_time,
            };
            return Ok(Run {
                outcome: ScfOutcome {
                    converged: stopped && verified_off,
                    result,
                },
                records,
                reference_stopped: stopped,
                exit_reason,
            });
        }
        std::mem::swap(&mut x, &mut next);
        input_planes.copy_from_slice(&next_planes);
        if comm.size() == 1 {
            plan.write_serial_full(&input_planes, &mut input)?;
        } else {
            let epoch = problem.owned.next_epoch();
            let owned = problem
                .owned
                .input
                .as_mut()
                .expect("owned input retained by map");
            owned.epoch = epoch;
            owned.values.copy_from_slice(&input_planes);
        }
        previous = Some(step.k_points);
    }
    Err(invalid("zero iteration budget").into())
}
