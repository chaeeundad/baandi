//! Ionic loop of the preconditioned L-BFGS (PRELAX-P2, ADR 0058).
//!
//! Same evaluation path as BFGS/FIRE (`Context::evaluate_with`: standard SCF
//! with verification, warm starts, optional C3 predictor, constraint
//! projection), same frames, restart and output files. Differences: a step may
//! evaluate several line-search trials; trials that are not accepted are kept
//! in `Restart::rejected` and republished at once; every attempted evaluation
//! carries a global `evaluation` index.
use super::{
    Context, Frame, Method, Restart, RunError, Settings, error, precon, predict, push_history,
    save, write_outputs,
};
use precon::{Judgement, LineSearch, Matrix, Model, State};
use serde_json::{Value, json};
use std::path::Path;

/// Preconditioner model of this input: lattice, per-atom element data and the
/// free-coordinate mask. Rejects elements without a covalent radius before any SCF.
pub(super) fn model(ctx: &Context, settings: &Settings) -> Result<Model, RunError> {
    debug_assert_eq!(settings.method, Method::PreconLbfgs);
    let kinds = ctx
        .config
        .system
        .atoms
        .iter()
        .map(|a| {
            let element = ctx
                .potentials
                .iter()
                .find(|(label, _)| *label == a.species)
                .map(|(_, p)| p.element().to_string())
                .ok_or_else(|| error(1, format!("no pseudopotential for {}", a.species)))?;
            precon::AtomKind::from_element(&element).map_err(|e| error(1, e))
        })
        .collect::<Result<Vec<_>, _>>()?;
    let mut free = vec![1.; 3 * kinds.len()];
    settings.project(&mut free);
    Ok(Model {
        lattice: ctx.cell.lattice_columns_bohr(),
        kinds,
        free: free.iter().map(|v| *v != 0.).collect(),
    })
}

fn build_record(matrix: &Matrix) -> Value {
    json!({"bonds": matrix.bonds, "min_bond_bohr": matrix.min_bond_bohr})
}

/// Optimizer memory at the initial frame (P built there, mu from the prior);
/// the initial frame gets its record (evaluation index and graph size).
pub(super) fn initial_state(
    model: &Model,
    positions: &[f64],
    evaluations: usize,
    frames: &mut [Frame],
) -> State {
    let matrix = Matrix::build(model, positions);
    let mut state = State::new(&matrix);
    state.evaluations = evaluations;
    frames[0].optimizer = Some(json!({
        "evaluation": evaluations, "mu_ha_bohr2": state.mu, "mu_source": state.mu_source,
        "line_search_trials": 0, "rejected_trials": 0, "rebuilt": true,
        "preconditioner": build_record(&matrix),
    }));
    state
}

/// results.json `relaxation.optimizer` block.
pub(super) fn summary(p: &State, state: &Restart) -> Value {
    json!({
        "method": "precon_lbfgs", "adr": "0058",
        "mu_ha_bohr2": p.mu, "mu_source": p.mu_source,
        "history_pairs": p.s.len(), "memory": precon::MEMORY,
        "rebuilds": p.rebuilds, "history_resets": p.history_resets,
        "attempted_evaluations": p.evaluations,
        "accepted_frames": state.frames.len(),
        "rejected_trials": p.rejected_trials,
        "failed_evaluations": state.failures.len(),
        "model": {
            "stiffness": "w = exp(-A (d/r0 - 1)), r0 = covalent radius sum (Cordero 2008)",
            "stiffness_decay_A": precon::STIFFNESS_DECAY,
            "stretch_fraction": precon::STRETCH_FRACTION,
            "stabiliser": precon::STABILISER,
            "bond_factor": precon::BOND_FACTOR,
            "mixed_bond_factor": precon::MIXED_BOND_FACTOR,
            "rebuild_fraction": precon::REBUILD_FRACTION,
            "prior_spring_ha_bohr2": precon::PRIOR_SPRING_HA_BOHR2,
            "mu_rule": "y.P1^-1 y / s.y of the latest stored pair",
            "armijo_c1": precon::ARMIJO_C1,
            "max_trials_per_direction": precon::MAX_TRIALS,
        },
    })
}

fn projected_gradient(settings: &Settings, frame: &Frame) -> Vec<f64> {
    let mut g: Vec<f64> = frame.forces_ha_bohr.iter().flatten().map(|f| -f).collect();
    settings.project(&mut g);
    // -0.0 from negating a zero force would make the step non-zero-signed;
    // normalise so that fixed coordinates stay exactly zero.
    for v in &mut g {
        if *v == 0. {
            *v = 0.;
        }
    }
    g
}

/// Run the preconditioned ionic loop from the accepted state (fresh or resumed).
pub(super) fn run(
    ctx: &mut Context,
    settings: &Settings,
    state: &mut Restart,
    model: &Model,
    out: &Path,
    restart_path: &Path,
    extrapolate: bool,
) -> Result<(), RunError> {
    let mut memory = state
        .precon
        .clone()
        .ok_or_else(|| error(1, "restart lacks the precon_lbfgs optimizer state"))?;
    let mut matrix = Matrix::build(model, &memory.build_positions);
    // Armijo noise tolerance: twice the SCF energy accuracy (standard
    // threshold, ADR 0057 default for relaxations; strict energy tolerance on
    // the oracle path).
    let noise = 2.
        * ctx
            .opts
            .standard
            .map_or(ctx.opts.energy_tolerance, |p| p.threshold_ha);
    let cap = settings
        .max_displacement_bohr
        .min(settings.trust_radius_bohr);
    while state.frames.len() - 1 < settings.max_steps && state.status == "running" {
        let started = std::time::Instant::now();
        let last = state.frames.last().unwrap().clone();
        ctx.position(&last.fractional);
        let source = if ctx.selection.reference() {
            None
        } else {
            Some(ctx.problem()?)
        };
        let g0 = projected_gradient(settings, &last);
        let mut search = LineSearch::begin(
            &mut memory,
            &matrix,
            &last.positions,
            last.energy_ha,
            &g0,
            cap,
        );
        let fallback_at_start = search.fallback;
        let (mut fresh, mut step_rejected, mut step_failed) = (false, 0usize, 0usize);
        let mut accepted = None;
        loop {
            let positions = search.trial();
            let fractional: Vec<_> = positions
                .chunks(3)
                .zip(last.positions.chunks(3))
                .zip(&last.fractional)
                .map(|((x, old), frac)| {
                    if x == old {
                        *frac
                    } else {
                        ctx.cell.cartesian_bohr_to_fractional([x[0], x[1], x[2]])
                    }
                })
                .collect();
            ctx.position(&fractional);
            let length = precon::max_atom_norm(&search.step());
            let warm = (!fresh && settings.warm_start)
                .then_some((source.as_ref(), state.density.as_slice()));
            let predicted = if extrapolate && warm.is_some() {
                predict(
                    &state.extrapolation,
                    &positions,
                    settings.start_extrapolation,
                )
            } else {
                None
            };
            memory.evaluations += 1;
            let evaluation = memory.evaluations;
            match ctx.evaluate_with(
                positions.clone(),
                length,
                warm,
                predicted.as_ref().map(|(d, kind)| (d.as_slice(), *kind)),
                extrapolate,
                settings,
            ) {
                Err(e) => {
                    step_failed += 1;
                    state
                        .failures
                        .push(json!({"step":state.frames.len(),"attempt":step_failed,
                        "evaluation":evaluation,"positions_bohr":positions,"fractional":fractional,
                        "reason":e.to_string()}));
                    if fresh {
                        // Second consecutive failure: as BFGS, keep the last good structure.
                        state.status = "failed".into();
                        break;
                    }
                    // Retry once from a fresh start at half the step.
                    fresh = true;
                    search.failed();
                }
                Ok((frame, density, n_at)) => {
                    fresh = false;
                    let g = projected_gradient(settings, &frame);
                    // A trial that already meets the force-only stop rule is a
                    // converged structure on real verified forces.
                    let stop = settings.energy_tol_ha.is_none() && settings.force_converged(&frame);
                    let judgement =
                        search.judge(&mut memory, &matrix, frame.energy_ha, &g, stop, noise);
                    if judgement == Judgement::Accept {
                        accepted = Some((frame, density, n_at, g, evaluation));
                        break;
                    }
                    step_rejected += 1;
                    let record = json!({"step":state.frames.len(),"trial":search.total_trials,
                        "evaluation":evaluation,"t":search.t,"positions_bohr":frame.positions,
                        "energy_ha":frame.energy_ha,"forces_ha_bohr":frame.forces_ha_bohr,
                        "max_force_ha_bohr":frame.max_force_ha_bohr,
                        "max_force_norm_ha_bohr":frame.max_force_norm_ha_bohr,
                        "energy_change_ha":frame.energy_ha - last.energy_ha,
                        "scf_iterations":frame.scf_iterations,"scf_h_vectors":frame.scf_h_vectors,
                        "provenance":frame.provenance,
                        "reason": if judgement == Judgement::Fail { "line_search_failed" } else { "armijo" }});
                    state.rejected.push(record.clone());
                    if judgement == Judgement::Fail {
                        state.status = "line_search_failed".into();
                        break;
                    }
                    ctx.position(&last.fractional);
                    write_outputs(ctx, settings, state, out, None, Some(&record))?;
                }
            }
        }
        if let Some((mut frame, density, n_at, g, evaluation)) = accepted {
            if extrapolate {
                push_history(&mut state.extrapolation, &frame.positions, &density, n_at);
            }
            let y: Vec<f64> = g.iter().zip(&g0).map(|(a, b)| a - b).collect();
            let stored = memory.store(&matrix, search.step(), y);
            let rebuilt = matrix.needs_rebuild(model, &frame.positions);
            let mut record = json!({
                "evaluation": evaluation, "mu_ha_bohr2": memory.mu, "mu_source": memory.mu_source,
                "t": search.t, "line_search_trials": search.total_trials,
                "rejected_trials": step_rejected, "failed_evaluations": step_failed,
                "fallback": search.fallback, "fallback_at_start": fallback_at_start,
                "pair_stored": stored, "history_pairs": memory.s.len(),
                "rebuilt": rebuilt, "armijo_noise_ha": noise,
            });
            if rebuilt {
                // L-BFGS pairs are kept: they are secant data of E, not of P (ADR 0058).
                matrix = Matrix::build(model, &frame.positions);
                memory.build_positions = frame.positions.clone();
                memory.rebuilds += 1;
                record["preconditioner"] = build_record(&matrix);
            }
            frame.optimizer = Some(record);
            let delta = frame.energy_ha - last.energy_ha;
            let (consecutive, converged) =
                settings.ionic_stop(&frame, Some(delta), state.consecutive);
            state.consecutive = consecutive;
            state.frames.push(frame);
            state.density = density;
            if converged {
                state.status = "converged".into();
            }
        } else {
            ctx.position(&last.fractional);
        }
        state.precon = Some(memory.clone());
        if state.status == "running" && state.frames.len() - 1 == settings.max_steps {
            state.status = "unconverged".into();
        }
        save(restart_path, state)?;
        write_outputs(ctx, settings, state, out, Some(started.elapsed()), None)?;
    }
    match state.status.as_str() {
        "converged" => Ok(()),
        "failed" => Err(error(
            4,
            "relax failed after one fresh, half-step retry; outputs retain last good structure",
        )),
        "line_search_failed" => Err(error(
            3,
            "relax line search found no acceptable point along the quasi-Newton and the preconditioned steepest-descent directions; outputs retain last good structure",
        )),
        _ => Err(error(3, "relax unconverged: max_steps reached")),
    }
}

#[cfg(test)]
mod tests {
    use super::super::bits;
    use super::*;
    /// Restart round trip (IEEE bits, as restart.json): the decoded optimizer
    /// memory is identical and gives a bit-identical next direction.
    #[test]
    fn restart_state_round_trip_is_bit_identical() {
        let l = [[9., 0.4, 0.], [0., 9.5, 0.3], [0.2, 0., 10.]];
        let si = precon::AtomKind::from_element("Si").unwrap();
        let model = Model {
            lattice: l,
            kinds: vec![si; 4],
            free: (0..12).map(|k| k != 4).collect(),
        };
        let x = vec![0.1, 0.2, 0.3, 2.3, 2.1, 2.4, 4.4, 0.3, 4.6, 0.2, 4.5, 4.3];
        let matrix = Matrix::build(&model, &x);
        let mut state = State::new(&matrix);
        // Two secant pairs from an anisotropic quadratic (Ha/Bohr^2).
        for k in 0..2 {
            let s: Vec<f64> = (0..12)
                .map(|i| {
                    if i == 4 {
                        0.
                    } else {
                        1e-2 * ((i + 3 * k) as f64).sin()
                    }
                })
                .collect();
            let y: Vec<f64> = s
                .iter()
                .enumerate()
                .map(|(i, v)| (0.05 + 0.01 * i as f64) * v)
                .collect();
            assert!(state.store(&matrix, s, y));
        }
        state.evaluations = 7;
        let mut value = serde_json::to_value(&state).unwrap();
        bits(&mut value, true).unwrap();
        let text = serde_json::to_string(&value).unwrap();
        let mut back: Value = serde_json::from_str(&text).unwrap();
        bits(&mut back, false).unwrap();
        let decoded: State = serde_json::from_value(back).unwrap();
        assert_eq!(decoded, state);
        let rebuilt = Matrix::build(&model, &decoded.build_positions);
        let g: Vec<f64> = (0..12)
            .map(|i| if i == 4 { 0. } else { 1e-3 * (i as f64 - 5.) })
            .collect();
        let a = state.inverse_hessian_apply(&matrix, &g);
        let b = decoded.inverse_hessian_apply(&rebuilt, &g);
        assert!(a.iter().zip(&b).all(|(a, b)| a.to_bits() == b.to_bits()));
        assert_eq!(a[4], 0.);
    }
}
