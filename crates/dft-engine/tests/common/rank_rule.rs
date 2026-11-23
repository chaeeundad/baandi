//! Test-only rank acceptance. Never used to choose a numerical solver branch.
use dft_engine::refscf::{driver::Run, record::SolveRecord};

#[derive(Debug)]
pub struct Branch {
    pub map: usize,
    pub verification: bool,
    pub call: usize,
    pub correction: usize,
    pub band: usize,
    pub changes: [f64; 2],
    pub thresholds: [f64; 2],
    pub margins: [f64; 2],
    pub rounding: f64,
}

fn require(ok: bool, reason: &str) -> Result<(), String> {
    if ok { Ok(()) } else { Err(reason.into()) }
}

fn close(a: f64, b: f64, atol: f64) -> bool {
    a.is_finite() && b.is_finite() && (a - b).abs() < atol + 1e-8 * b.abs()
}

/// Compare the chronological decisions, not only the final iteration totals.
/// Every predicate that differs at the first divergent correction must straddle
/// its own threshold within the same rounding budget. No material/rank special case.
pub fn calls(a: &[SolveRecord], b: &[SolveRecord]) -> Result<Vec<Branch>, String> {
    for (call, (a, b)) in a.iter().zip(b).enumerate() {
        require(a.first_step == b.first_step, "first-step branch differs")?;
        require(close(a.ethr_ry, b.ethr_ry, 1e-13), "call threshold differs")?;
        let bands = a.band_thresholds_ha.len();
        for c in [a, b] {
            require(
                c.band_thresholds_ha.len() == bands && bands > 0,
                "band count differs",
            )?;
            require(
                c.ritz_changes_ha.len() == c.iterations
                    && c.ritz_values_ha.len() == c.iterations
                    && c.corrections.len() == c.iterations,
                "incomplete decision trace",
            )?;
        }
        for iteration in 0..a.iterations.min(b.iterations) {
            let (aw, bw) = (&a.corrections[iteration], &b.corrections[iteration]);
            require(
                aw.active_bands == bw.active_bands
                    && aw.basis_before == bw.basis_before
                    && aw.h_columns == bw.h_columns
                    && aw.reduced_order == bw.reduced_order
                    && aw.repairs == bw.repairs
                    && aw.dropped_corrections == bw.dropped_corrections,
                "work diverged before Ritz decision",
            )?;
            let mut branches = Vec::new();
            for c in [a, b] {
                require(
                    c.ritz_changes_ha[iteration].len() == bands
                        && c.ritz_values_ha[iteration].len() == bands,
                    "incomplete band trace",
                )?;
            }
            for band in 0..bands {
                let changes = [
                    a.ritz_changes_ha[iteration][band],
                    b.ritz_changes_ha[iteration][band],
                ];
                let values = [
                    a.ritz_values_ha[iteration][band],
                    b.ritz_values_ha[iteration][band],
                ];
                let thresholds = [a.band_thresholds_ha[band], b.band_thresholds_ha[band]];
                require(
                    values
                        .iter()
                        .flatten()
                        .chain(&changes)
                        .chain(&thresholds)
                        .all(|v| v.is_finite())
                        && thresholds.iter().all(|t| *t > 0.),
                    "nonfinite/invalid Ritz evidence",
                )?;
                for j in 0..2 {
                    require(
                        changes[j].to_bits() == (values[j][1] - values[j][0]).abs().to_bits(),
                        "Ritz change inconsistent with eigenvalues",
                    )?;
                }
                // Cancellation subtracts eigenvalues, so scale by the operands,
                // not the tiny change. 32 eps is an audit envelope, not a solver tolerance.
                let rounding =
                    32. * f64::EPSILON * values.iter().flatten().fold(1_f64, |s, v| s.max(v.abs()));
                require(
                    close(changes[0], changes[1], 1e-11),
                    "Ritz trajectory differs before branch",
                )?;
                require(
                    (0..2).all(|j| close(values[0][j], values[1][j], 1e-11)),
                    "Ritz eigenvalues differ before branch",
                )?;
                let margins = [changes[0] - thresholds[0], changes[1] - thresholds[1]];
                if (margins[0] < 0.) != (margins[1] < 0.) {
                    require(
                        margins.iter().all(|m| m.abs() <= rounding)
                            && (changes[0] - changes[1]).abs() <= rounding
                            && (thresholds[0] - thresholds[1]).abs() <= rounding,
                        "unexplained Ritz branch away from threshold",
                    )?;
                    branches.push(Branch {
                        map: 0,
                        verification: false,
                        call: call + 1,
                        correction: iteration + 1,
                        band,
                        changes,
                        thresholds,
                        margins,
                        rounding,
                    });
                }
            }
            if !branches.is_empty() {
                return Ok(branches);
            }
            require(
                aw.restarts == bw.restarts,
                "restart differs without Ritz branch",
            )?;
        }
        require(
            a.iterations == b.iterations
                && a.active_bands == b.active_bands
                && a.h_columns == b.h_columns
                && a.reduced_orders == b.reduced_orders
                && a.restarts == b.restarts
                && a.repairs == b.repairs
                && a.dropped_corrections == b.dropped_corrections
                && a.unconverged == b.unconverged,
            "call work differs without Ritz branch",
        )?;
    }
    require(a.len() == b.len(), "call count differs without Ritz branch")?;
    Ok(Vec::new())
}

pub fn trajectory(a: &Run, b: &Run) -> Result<Vec<Branch>, String> {
    for (a_map, b_map) in a.records.iter().zip(&b.records) {
        require(
            a_map.iteration == b_map.iteration
                && close(a_map.ethr_requested_ry, b_map.ethr_requested_ry, 1e-13),
            "map input differs before branch",
        )?;
        for verification in [false, true] {
            let (ac, bc) = if verification {
                (
                    &a.verification_calls[a_map.verification_calls.clone()],
                    &b.verification_calls[b_map.verification_calls.clone()],
                )
            } else {
                (a_map.calls.as_slice(), b_map.calls.as_slice())
            };
            let mut branches = calls(ac, bc)
                .map_err(|e| format!("map={} verification={verification}: {e}", a_map.iteration))?;
            if !branches.is_empty() {
                for branch in &mut branches {
                    branch.map = a_map.iteration;
                    branch.verification = verification;
                }
                return Ok(branches);
            }
            if !verification {
                require(
                    close(a_map.dr2_ry, b_map.dr2_ry, 1e-11)
                        && close(a_map.ethr_ry, b_map.ethr_ry, 1e-13)
                        && a_map.first_iteration_redo == b_map.first_iteration_redo
                        && a_map.davidson_iterations == b_map.davidson_iterations
                        && a_map.h_columns == b_map.h_columns,
                    "map output differs without Ritz branch",
                )?;
            }
        }
    }
    require(
        a.records.len() == b.records.len(),
        "map count differs without Ritz branch",
    )?;
    Ok(Vec::new())
}

/// Gates remain mandatory even after an explained branch. All attempts and their
/// recovery remain in the driver's delivered wall; never time just the last pass.
pub fn delivered(run: &Run, strict_energy: f64, electrons: f64) -> bool {
    let r = &run.outcome.result;
    let Some(v) = &r.verification else {
        return false;
    };
    let Some(last) = v.attempts.last() else {
        return false;
    };
    run.outcome.converged
        && run.reference_stopped
        && v.passed()
        && last.iteration == run.records.len()
        && last.threshold_ha <= 5e-9
        && last.coulomb_error_ha.is_finite()
        && last.coulomb_error_ha <= last.threshold_ha
        && last
            .residual_bound_ha
            .is_some_and(|b| b.is_finite() && b <= last.threshold_ha)
        && (r.total_hartree() - strict_energy).abs() < 5e-9
        && (r.charge - electrons).abs() < 1e-10
        && r.density
            .iter()
            .chain(&r.eigenvalues)
            .all(|v| v.is_finite())
        && run
            .records
            .iter()
            .all(|r| r.energy.total_hartree().is_finite())
        && v.wall() <= r.total_time
        && v.attempts
            .iter()
            .map(|a| a.operator_applications)
            .sum::<usize>()
            == run
                .verification_calls
                .iter()
                .map(|c| c.h_columns)
                .sum::<usize>()
}

pub fn accept(a: &Run, b: &Run, strict_energy: f64, electrons: f64) -> Result<Vec<Branch>, String> {
    require(
        delivered(a, strict_energy, electrons) && delivered(b, strict_energy, electrons),
        "delivered-state checks failed",
    )?;
    require(
        a.outcome.result.density.shape() == b.outcome.result.density.shape(),
        "density shape differs",
    )?;
    let branches = trajectory(a, b)?;
    // A verification-only branch cannot explain a difference in the already
    // delivered density. Keep its cap unless an ordinary map diverged first.
    if branches.first().is_none_or(|b| b.verification) {
        require(
            density_delta(a, b) < 1e-9,
            "density differs without explained branch",
        )?;
    }
    Ok(branches)
}

pub fn density_delta(a: &Run, b: &Run) -> f64 {
    a.outcome
        .result
        .density
        .iter()
        .zip(&b.outcome.result.density)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max)
}

/// Debug's float formatting round-trips finite f64, including signed zero.
/// Strip only clocks; compare the full delivered arrays and scientific records.
pub fn fingerprint(run: &Run) -> String {
    let mut records = run.records.clone();
    for r in &mut records {
        r.map_s = 0.;
        r.solve_s = 0.;
        r.density_s = 0.;
        r.potential_s = 0.;
        r.energy_s = 0.;
        r.density_transform_s = 0.;
        r.metric_s = 0.;
        r.mixing_s = 0.;
        r.correction_s = 0.;
        r.verification_s = 0.;
    }
    let r = &run.outcome.result;
    let mut history = r.history.clone();
    for h in &mut history {
        h.timings = Default::default();
        if let Some(c) = &mut h.fixed_point_certificate {
            c.timings = Default::default();
        }
    }
    let mut verification = r.verification.clone();
    if let Some(v) = &mut verification {
        for a in &mut v.attempts {
            a.wall = Default::default();
        }
    }
    format!(
        "{:?}",
        (
            run.outcome.converged,
            run.reference_stopped,
            run.first_reference_stop,
            &run.exit_reason,
            records,
            &run.verification_calls,
            history,
            verification,
            (
                &r.energy,
                (
                    r.density.shape(),
                    r.density.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                ),
                &r.eigenvalues,
                &r.orbitals,
                &r.occupations,
                r.charge,
                &r.k_points
            )
        )
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_engine::refscf::record::CorrectionRecord;

    fn call(change: f64) -> SolveRecord {
        SolveRecord {
            ethr_ry: 1e-8,
            first_step: false,
            iterations: 1,
            active_bands: vec![1],
            h_columns: 2,
            reduced_orders: vec![1, 2],
            restarts: 0,
            unconverged: usize::from(change >= 5e-9),
            ritz_changes_ha: vec![vec![change]],
            inner_clocks: None,
            ritz_values_ha: vec![vec![[0., change]]],
            band_thresholds_ha: vec![5e-9],
            repairs: 0,
            dropped_corrections: 0,
            corrections: vec![CorrectionRecord {
                active_bands: vec![0],
                basis_before: 1,
                h_columns: 1,
                reduced_order: 2,
                repairs: 0,
                dropped_corrections: 0,
                restarts: 0,
            }],
        }
    }

    #[test]
    fn explained_rounding_branch_records_both_margins() {
        let a = call(5e-9 - 4e-16);
        let b = call(5e-9 + 5e-16);
        let branches = calls(&[a], &[b]).unwrap();
        assert_eq!(branches.len(), 1);
        assert!(branches[0].margins[0] < 0. && branches[0].margins[1] > 0.);
        assert_eq!(branches[0].correction, 1);
    }

    #[test]
    fn rejects_unexplained_divergence_away_from_threshold() {
        let error = calls(&[call(5e-9 - 1e-12)], &[call(5e-9 + 1e-12)]).unwrap_err();
        assert!(error.contains("away from threshold"), "{error}");
    }

    #[test]
    fn later_rounding_cannot_excuse_earlier_work_divergence() {
        let a = call(5e-9 - 4e-16);
        let mut b = call(5e-9 + 5e-16);
        b.corrections[0].h_columns += 1;
        assert!(
            calls(&[a], &[b])
                .unwrap_err()
                .contains("before Ritz decision")
        );
    }

    #[test]
    fn missing_or_nonfinite_evidence_fails_closed() {
        let a = call(1e-9);
        let mut b = a.clone();
        b.ritz_values_ha.clear();
        assert!(calls(std::slice::from_ref(&a), &[b]).is_err());
        assert!(calls(&[a], &[call(f64::NAN)]).is_err());
    }
}
