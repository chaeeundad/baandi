//! QE trajectory adapter. I/O and JSON stay outside numerical kernels.
use dft_cli::{
    config::Config,
    run::{RunError, reference_profile},
};
use dft_engine::refscf::{
    broyden::Broyden,
    density::{DensityPlan, g_squared},
    metric::{Metric, Serial},
    qe_reference::{QeDensity, printed_accuracy},
};
use dft_mpi::Collective;
use dft_physics::data::Complex64 as C;
use serde_json::{Value, json};
use std::{
    fs,
    path::{Path, PathBuf},
};

fn fail(s: impl std::fmt::Display) -> RunError {
    RunError {
        code: 4,
        cause: format!("refscf-check first failing term: {s}"),
    }
}
fn files(dir: &Path, name: &str, depth: usize, out: &mut Vec<PathBuf>) -> Result<(), RunError> {
    if dir.is_file() {
        if dir.file_name().is_some_and(|n| n == name) {
            out.push(dir.to_path_buf());
        }
        return Ok(());
    }
    if depth == 0 || !dir.is_dir() {
        return Ok(());
    }
    for entry in fs::read_dir(dir).map_err(fail)? {
        let entry = entry.map_err(fail)?;
        if entry.file_type().map_err(fail)?.is_symlink() {
            continue;
        }
        files(&entry.path(), name, depth - 1, out)?;
    }
    Ok(())
}
pub(crate) fn density(dir: &Path) -> Result<QeDensity, RunError> {
    let mut found = Vec::new();
    files(dir, "charge-density.dat", 4, &mut found)?;
    found.sort();
    if found.len() != 1 {
        return Err(fail(format!(
            "{}: expected exactly one charge-density.dat, found {}",
            dir.display(),
            found.len()
        )));
    }
    QeDensity::decode(&fs::read(&found[0]).map_err(fail)?).map_err(fail)
}
pub(crate) fn stdout(dir: &Path, root: &Path, k: usize) -> Result<Option<String>, RunError> {
    for name in [
        "stdout",
        "stdout.txt",
        "pw.out",
        "qe.out",
        "output.out",
        "run.out",
    ] {
        let path = dir.join(name);
        if path.is_file() {
            return fs::read_to_string(path).map(Some).map_err(fail);
        }
    }
    let snapshot = dir.file_name().and_then(|s| s.to_str());
    let fallback = format!("k{k}-mix");
    let snapshot = snapshot.unwrap_or(&fallback);
    for name in [format!("{snapshot}.out"), format!("{snapshot}.stdout")] {
        let path = root.join(name);
        if path.is_file() {
            return fs::read_to_string(path).map(Some).map_err(fail);
        }
    }
    Ok(None)
}
#[derive(Clone, Copy)]
pub struct Tolerance {
    pub absolute: f64,
    pub relative: f64,
}
pub fn execute(
    input: &Path,
    root: &Path,
    steps: usize,
    tolerance: Tolerance,
) -> Result<(), RunError> {
    if steps == 0
        || !tolerance.absolute.is_finite()
        || tolerance.absolute < 0.0
        || !tolerance.relative.is_finite()
        || tolerance.relative < 0.0
    {
        return Err(fail("invalid steps/tolerance"));
    }
    let config = Config::parse(&fs::read_to_string(input).map_err(fail)?).map_err(fail)?;
    if let Some(comm) = dft_engine::parallel::SharedComm::active() {
        replay(&config, root, steps, tolerance, &comm)
    } else {
        replay(&config, root, steps, tolerance, &Serial)
    }
}
fn replay<K: Collective>(
    config: &Config,
    root: &Path,
    steps: usize,
    tol: Tolerance,
    comm: &K,
) -> Result<(), RunError> {
    let profile = reference_profile(config)?;
    if profile.beta == 1.0 {
        return Err(fail(
            "rho_in(1) recovery requires beta != 1; beta=1 has no inverse startup",
        ));
    }
    let prepared = config.prepare().map_err(fail)?;
    let volume = prepared.cell.volume_bohr3();
    let first = density(&root.join("k1-out"))?;
    first
        .validate_cell(prepared.cell.reciprocal_columns_inverse_bohr())
        .map_err(fail)?;
    let order = first.order();
    let plan = DensityPlan::new(
        comm,
        &prepared.cell,
        prepared.grid.dims(),
        2.0 * config.basis.ecutwfc_ry,
        2.0 * config.basis.ecutwfc_ry,
    )
    .map_err(fail)?;
    // Validate complete membership first; then each rank retains only its grid
    // plan's G rows for the actual production mixer and history.
    let own: std::collections::HashSet<_> = plan.miller.iter().copied().collect();
    let local_indices: Vec<_> = order
        .iter()
        .enumerate()
        .filter_map(|(i, m)| own.contains(m).then_some(i))
        .collect();
    let mut counts = vec![0.; comm.size()];
    comm.allgatherv_f64(
        &[local_indices.len() as f64],
        &vec![1; comm.size()],
        &mut counts,
    )
    .map_err(fail)?;
    if counts.iter().sum::<f64>() != order.len() as f64 || local_indices.len() != plan.miller.len()
    {
        return Err(fail(
            "Miller membership: input density sphere and QE differ; no entries dropped",
        ));
    }
    let local_order: Vec<_> = local_indices.iter().map(|&i| order[i]).collect();
    let g2: Vec<_> = local_order
        .iter()
        .map(|&m| g_squared(first.reciprocal, m))
        .collect();
    let align = |d: QeDensity| -> Result<Vec<C>, RunError> {
        d.validate_cell(first.reciprocal).map_err(fail)?;
        let values = d.aligned(&order).map_err(fail)?;
        Ok(local_indices.iter().map(|&i| values[i]).collect())
    };
    let out1 = align(first.clone())?;
    let mix1 = align(density(&root.join("k1-mix"))?)?;
    let mut x: Vec<_> = mix1
        .iter()
        .zip(&out1)
        .map(|(m, y)| (*m - profile.beta * y) / (1.0 - profile.beta))
        .collect();
    let zero: Vec<_> = local_order
        .iter()
        .enumerate()
        .filter_map(|(i, m)| (*m == [0; 3]).then_some(i))
        .collect();
    let local_charge = zero.iter().map(|&i| x[i].re * volume).sum::<f64>();
    let mut charge = [0.];
    comm.tree_sum_f64(&[local_charge], &mut charge)
        .map_err(fail)?;
    if charge[0] <= 0.0 {
        return Err(fail("nonpositive recovered total charge"));
    }
    let mut metric =
        Metric::new(volume, &g2, &vec![true; g2.len()], comm.size(), 350).map_err(fail)?;
    let mut mixer =
        Broyden::new(profile, &g2, &vec![true; g2.len()], volume, charge[0]).map_err(fail)?;
    let mut f = vec![C::default(); g2.len()];
    let mut next = f.clone();
    let all: Vec<_> = (0..g2.len()).collect();
    let mut shell_g2 = Vec::<f64>::new();
    for &m in &order {
        let g = g_squared(first.reciprocal, m);
        if g > 0.0
            && shell_g2
                .last()
                .is_none_or(|p| (g - p).abs() > 1e-10 * g.max(1.0))
        {
            shell_g2.push(g);
            if shell_g2.len() == 10 {
                break;
            }
        }
    }
    let groups: Vec<_> = shell_g2
        .iter()
        .map(|&g| {
            g2.iter()
                .enumerate()
                .filter_map(|(i, &v)| ((v - g).abs() <= 1e-10 * g.max(1.0)).then_some(i))
                .collect::<Vec<_>>()
        })
        .collect();
    let mut beta_one_checked = false;
    for name in ["k1-mix-beta1.0", "k1-beta1", "k1-beta-1", "beta1"] {
        let path = root.join(name);
        if path.exists() {
            let v = align(density(&path)?)?;
            report_difference(comm, "startup beta=1 output", &out1, &v, &all, tol)?;
            beta_one_checked = true;
        }
    }
    let beta03 = root.join("k1-mix-beta0.3");
    if beta03.exists() {
        let mixed = align(density(&beta03)?)?;
        let recovered: Vec<_> = mixed
            .iter()
            .zip(&out1)
            .map(|(m, y)| (*m - 0.3 * y) / 0.7)
            .collect();
        report_difference(comm, "startup beta=0.3 recovery", &x, &recovered, &all, tol)?;
    }
    // An independently saved atomic/input density is a stronger recovery check.
    let initial = root.join("k1-in");
    if initial.exists() {
        let v = align(density(&initial)?)?;
        report_difference(comm, "recovered rho_in(1)", &x, &v, &all, tol)?;
    }
    if comm.rank() == 0 {
        println!(
            "{}",
            json!({"reader":"QE four-record Gamma/nspin=1","g_vectors":order.len(),"ranks":comm.size(),"charge_electrons":charge[0],"beta_one_cross_checked":beta_one_checked,"tolerance":{"absolute_e_bohr3":tol.absolute,"relative":tol.relative,"status":"provisional; calibrate with QE rank-count replay"}})
        );
    }
    for k in 1..=steps {
        let y = align(density(&root.join(format!("k{k}-out")))?)?;
        let expected = align(density(&root.join(format!("k{k}-mix")))?)?;
        for ((f, y), x) in f.iter_mut().zip(&y).zip(&x) {
            *f = y - x;
        }
        let dr2 = metric.norm(comm, &f).map_err(fail)?;
        let printed = stdout(&root.join(format!("k{k}-mix")), root, k)?
            .as_deref()
            .and_then(printed_accuracy);
        if let Some(p) = printed {
            let allowance = p.rounding_ry + tol.absolute + tol.relative * p.ry.abs();
            let bad = if p.upper_bound {
                dr2 > p.ry + allowance
            } else {
                (dr2 - p.ry).abs() > allowance
            };
            if bad {
                return Err(fail(format!(
                    "k={k} dr2: BAANDI {dr2:.16e} Ry, QE {} {:.16e} Ry (printed rounding {:.3e} Ry)",
                    if p.upper_bound { "<" } else { "=" },
                    p.ry,
                    p.rounding_ry
                )));
            }
        }
        let d = mixer
            .update(comm, &mut metric, &x, &f, &mut next)
            .map_err(|e| fail(format!("k={k} A/gamma: {e}")))?;
        let global = report_difference(comm, &format!("k={k} all G"), &next, &expected, &all, tol)?;
        let g0 = report_difference(comm, &format!("k={k} G=0"), &next, &expected, &zero, tol)?;
        let mut shells = Vec::new();
        for (s, indices) in groups.iter().enumerate() {
            let diff = report_difference(
                comm,
                &format!("k={k} shell {}", s + 1),
                &next,
                &expected,
                indices,
                tol,
            )?;
            shells.push(json!({"g2_bohr_minus2":shell_g2[s],"difference":diff}));
        }
        let local: Vec<_> = zero
            .iter()
            .map(|&i| [next[i].re * volume, expected[i].re * volume])
            .collect();
        let mut total = [0.; 2];
        comm.tree_sum_f64(local.first().unwrap_or(&[0.; 2]), &mut total)
            .map_err(fail)?;
        let charge_error = (total[0] - total[1]).abs();
        if charge_error > tol.absolute * volume + tol.relative * total[1].abs() {
            return Err(fail(format!(
                "k={k} total charge: {charge_error:e} electrons"
            )));
        }
        if comm.rank() == 0 {
            println!(
                "{}",
                json!({"iteration":k,"all_g":global,"g0":g0,"total_charge":{"actual_e":total[0],"expected_e":total[1],"difference_e":charge_error},"lowest_10_nonzero_shells":shells,"dr2_ry":dr2,"dr2_ha":dr2/2.0,"qe_printed_dr2_ry":printed.map(|p|p.ry),"qe_printed_is_upper_bound":printed.map(|p|p.upper_bound),"history_occupancy":d.occupancy,"gamma_l2":d.gamma_l2,"gamma_max":d.gamma_max,"a_condition_2":d.condition})
            );
        }
        // Teacher forcing: use QE's own x sequence, never our previous output.
        x = expected;
    }
    Ok(())
}
fn report_difference<K: Collective>(
    comm: &K,
    label: &str,
    a: &[C],
    b: &[C],
    indices: &[usize],
    tol: Tolerance,
) -> Result<Value, RunError> {
    let (mut max, mut sum, mut scale, mut norm) = (0.0_f64, 0., 0.0_f64, 0.);
    for &i in indices {
        let d = (a[i] - b[i]).norm();
        max = max.max(d);
        sum += d * d;
        scale = scale.max(b[i].norm());
        norm += b[i].norm_sqr();
    }
    let mut values = vec![0.; 5 * comm.size()];
    comm.allgatherv_f64(
        &[max, sum, scale, norm, indices.len() as f64],
        &vec![5; comm.size()],
        &mut values,
    )
    .map_err(fail)?;
    let (mut max, mut sum, mut scale, mut norm, mut n) = (0.0_f64, 0., 0.0_f64, 0., 0.);
    for v in values.as_chunks::<5>().0 {
        max = max.max(v[0]);
        sum += v[1];
        scale = scale.max(v[2]);
        norm += v[3];
        n += v[4];
    }
    let rms = (sum / n.max(1.)).sqrt();
    let relative = if scale > 0. { max / scale } else { max };
    if max > tol.absolute + tol.relative * scale {
        return Err(fail(format!(
            "{label}: max={max:.16e} e/bohr^3, RMS={rms:.16e} e/bohr^3, relative={relative:.6e}, limit={:.6e} e/bohr^3",
            tol.absolute + tol.relative * scale
        )));
    }
    Ok(
        json!({"max_e_bohr3":max,"rms_e_bohr3":rms,"relative_max":relative,"relative_rms":if norm>0. {(sum/norm).sqrt()} else {sum.sqrt()},"count":n as usize}),
    )
}

/// Offline reader: never constructs a problem or changes the verification verdict.
/// The file must be a completed results.json with the additive clock contract.
pub fn stop_to_verified(path: &Path) -> Result<(), RunError> {
    let doc: Value = serde_json::from_slice(&fs::read(path).map_err(fail)?).map_err(fail)?;
    let summary = stop_to_verified_summary(&doc)?;
    if dft_cli::parallel::rank() == 0 {
        println!("{summary}");
    }
    Ok(())
}

#[derive(serde::Deserialize)]
struct WallVerificationCost {
    gate_enabled: bool,
    skipped: bool,
    #[serde(default)]
    first_gate_enabled: bool,
    first_prediction_ha: Option<f64>,
    contraction_prediction_ha: Option<f64>,
    predicted_error_ha: Option<f64>,
    actual_error_ha: Option<f64>,
}
#[derive(serde::Deserialize)]
struct WallMap {
    iteration: usize,
    dr2_ry: f64,
    map_start_s: f64,
    estimator_s: f64,
    verification_start_s: Option<f64>,
    verification_end_s: Option<f64>,
    map_end_s: f64,
    verification_cost: Option<WallVerificationCost>,
}
#[derive(serde::Deserialize)]
struct WallAttempt {
    iteration: usize,
    passed: bool,
}
#[derive(serde::Deserialize)]
struct WallVerification {
    mode: String,
    passed: bool,
    attempts: Vec<WallAttempt>,
}
#[derive(serde::Deserialize)]
struct WallRun {
    stop_to_verified_clock: String,
    threshold_ry: f64,
    iterations: Vec<WallMap>,
    verification: Option<WallVerification>,
}

pub fn stop_to_verified_summary(doc: &Value) -> Result<Value, RunError> {
    let require = |ok, why| if ok { Ok(()) } else { Err(fail(why)) };
    require(doc["run"]["completed"] == true, "incomplete results record")?;
    let run: WallRun = serde_json::from_value(doc["reference_scf"].clone()).map_err(fail)?;
    require(
        run.stop_to_verified_clock == "reference-scf-wall-1",
        "unsupported or missing driver clock contract",
    )?;
    require(
        run.threshold_ry.is_finite() && run.threshold_ry > 0. && !run.iterations.is_empty(),
        "invalid threshold or empty map history",
    )?;
    let attempts = run
        .verification
        .as_ref()
        .map_or(&[][..], |v| v.attempts.as_slice());
    if let Some(v) = &run.verification {
        require(v.mode == "fixed_point", "verification is not fixed_point")?;
        require(
            v.passed == attempts.last().is_some_and(|a| a.passed),
            "inconsistent verification verdict",
        )?;
    }
    let mut first_stop = None;
    let mut passing_map = None;
    let mut interval_start = None;
    let mut interval_end = None;
    let mut previous_end = 0.;
    let mut attempted = 0;
    let mut skipped = 0;
    let mut gate_enabled = false;
    let mut failed = false;
    let mut verification_s = 0.;
    let mut recovery_s = 0.;
    for (i, map) in run.iterations.iter().enumerate() {
        require(
            map.iteration == i + 1,
            "missing, duplicate or unordered map",
        )?;
        require(passing_map.is_none(), "map after passing verification")?;
        require(
            [map.map_start_s, map.estimator_s, map.map_end_s, map.dr2_ry]
                .into_iter()
                .all(|x| x.is_finite() && x >= 0.),
            "nonfinite or negative map clock/estimator",
        )?;
        require(
            previous_end <= map.map_start_s
                && map.map_start_s <= map.estimator_s
                && map.estimator_s <= map.map_end_s,
            "overlapping maps or reversed map boundaries",
        )?;
        previous_end = map.map_end_s;
        let stopped = map.dr2_ry < run.threshold_ry;
        let skip = map.verification_cost.as_ref().is_some_and(|c| c.skipped);
        let gate = map
            .verification_cost
            .as_ref()
            .is_some_and(|c| c.gate_enabled);
        gate_enabled |= gate;
        require(
            !skip
                || (gate
                    && stopped
                    && run.verification.is_some()
                    && map.verification_start_s.is_none()
                    && map.verification_end_s.is_none()),
            "invalid skipped verification",
        )?;
        skipped += usize::from(skip);
        if stopped && first_stop.is_none() {
            first_stop = Some(map.iteration);
            interval_start = Some(map.estimator_s);
        }
        match (map.verification_start_s, map.verification_end_s) {
            (Some(start), Some(end)) => {
                require(
                    stopped && run.verification.is_some(),
                    "attempt without estimator stop or verification mode",
                )?;
                require(
                    start.is_finite()
                        && end.is_finite()
                        && map.estimator_s <= start
                        && start <= end
                        && end <= map.map_end_s,
                    "invalid verification boundaries",
                )?;
                let attempt = attempts
                    .get(attempted)
                    .ok_or_else(|| fail("missing verification verdict"))?;
                require(attempt.iteration == map.iteration, "attempt/map mismatch")?;
                if failed {
                    recovery_s += start - map.map_start_s;
                }
                verification_s += end - start;
                attempted += 1;
                if attempt.passed {
                    passing_map = Some(map.iteration);
                    interval_end = Some(end);
                } else {
                    failed = true;
                    recovery_s += map.map_end_s - end;
                }
            }
            (None, None) => {
                require(
                    !stopped || run.verification.is_none() || skip,
                    "missing verification timestamps",
                )?;
                if failed {
                    recovery_s += map.map_end_s - map.map_start_s;
                }
            }
            _ => return Err(fail("incomplete verification timestamp pair")),
        }
    }
    require(
        attempted == attempts.len(),
        "extra or unordered verification verdicts",
    )?;
    // Exhausted/off runs end at the last completed map, never a fabricated pass.
    let end = interval_start.map(|_| interval_end.unwrap_or(previous_end));
    let interval = interval_start.zip(end).map(|(a, b)| b - a);
    let unassigned = interval.map(|total| total - verification_s - recovery_s);
    let closure_error = interval
        .zip(unassigned)
        .map(|(total, other)| verification_s + recovery_s + other - total);
    if let Some(total) = interval {
        let tolerance = 64. * f64::EPSILON * end.unwrap().max(1.);
        require(
            unassigned.unwrap() >= -tolerance && closure_error.unwrap().abs() <= tolerance,
            "stop-to-verified parts do not close",
        )?;
        require(total >= 0., "negative stop-to-verified interval")?;
    }
    let mut summary = json!({
        "schema":"stop-to-verified-1",
        "status":if passing_map.is_some() {"VERIFIED"} else {"UNVERIFIED"},
        "first_stop_map":first_stop,"passing_map":passing_map,
        "attempt_count":attempted,
        "clock":"output rank monotonic seconds since driver start; no rank maxima",
        "interval_start_s":interval_start,"interval_end_s":end,"interval_wall_s":interval,
        "interval_end_kind":if passing_map.is_some() {"passing_verification_end"}
            else if first_stop.is_some() {"last_completed_map_unverified"} else {"no_estimator_stop"},
        "verification_attempts_s":interval.map(|_| verification_s),
        "ordinary_scf_after_failure_s":interval.map(|_| recovery_s),
        "unassigned_s":unassigned,"closure_error_s":closure_error,
        "parts_sum_checked":interval.is_some(),
        "diagnostic_on_off_delta":"separate measurement; may overlap this interval"
    });
    if gate_enabled {
        summary["skipped"] = json!(skipped);
    }
    if run.iterations.iter().any(|m| {
        m.verification_cost
            .as_ref()
            .is_some_and(|c| c.first_gate_enabled)
    }) {
        let mut predictions = Vec::new();
        for map in &run.iterations {
            if let Some(cost) = &map.verification_cost {
                require(
                    !cost.first_gate_enabled || cost.gate_enabled,
                    "first gate without verification gate",
                )?;
                for value in [
                    cost.first_prediction_ha,
                    cost.contraction_prediction_ha,
                    cost.predicted_error_ha,
                    cost.actual_error_ha,
                ]
                .into_iter()
                .flatten()
                {
                    require(
                        value.is_finite() && value >= 0.,
                        "invalid verification prediction/error",
                    )?;
                }
                if map.dr2_ry < run.threshold_ry {
                    predictions.push(json!({
                        "iteration":map.iteration, "skipped":cost.skipped,
                        "first_prediction_ha":cost.first_prediction_ha,
                        "contraction_prediction_ha":cost.contraction_prediction_ha,
                        "predicted_error_ha":cost.predicted_error_ha,
                        "actual_error_ha":cost.actual_error_ha
                    }));
                }
            }
        }
        summary["verification_predictions"] = json!(predictions);
    }
    Ok(summary)
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::f64::consts::{PI, TAU};
    fn write_density(path: &Path, reciprocal: [[f64; 3]; 3], miller: &[[i32; 3]], rho: &[C]) {
        fs::create_dir_all(path).unwrap();
        let mut bytes = Vec::new();
        let mut record = |p: Vec<u8>| {
            bytes.extend((p.len() as i32).to_le_bytes());
            bytes.extend(&p);
            bytes.extend((p.len() as i32).to_le_bytes());
        };
        record(
            [1, miller.len() as i32, 1]
                .into_iter()
                .flat_map(i32::to_le_bytes)
                .collect(),
        );
        record(
            (0..3)
                .flat_map(|j| (0..3).flat_map(move |i| reciprocal[i][j].to_le_bytes()))
                .collect(),
        );
        record(
            miller
                .iter()
                .flatten()
                .flat_map(|i| i.to_le_bytes())
                .collect(),
        );
        record(
            rho.iter()
                .flat_map(|c| [c.re, c.im])
                .flat_map(f64::to_le_bytes)
                .collect(),
        );
        fs::write(path.join("charge-density.dat"), bytes).unwrap();
    }
    #[test]
    fn synthetic_cli_replay_and_first_failing_term() {
        let dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
            .join("../../target/scratch")
            .join(format!("refscf-cli-test-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let text = format!(
            "[system]\nname='synthetic Fourier oracle'\nlattice_unit='bohr'\nlattice_vectors=[[{TAU},0,0],[0,{TAU},0],[0,0,{TAU}]]\n[[system.atoms]]\nspecies='Xx'\nfractional=[0,0,0]\n[pseudopotentials]\nXx='unused-synthetic.UPF'\n[basis]\necutwfc_ry=0.5\nkpoints=[[0,0,0]]\n[xc]\nfunctional='lda-pz'\n[scf]\ndriver='reference'\nmixing_beta=0.5\n"
        );
        let input = dir.join("input.toml");
        fs::write(&input, &text).unwrap();
        let config = Config::parse(&text).unwrap();
        let prep = config.prepare().unwrap();
        let plan = DensityPlan::new(&Serial, &prep.cell, prep.grid.dims(), 1., 1.).unwrap();
        let volume = prep.cell.volume_bohr3();
        let b = prep.cell.reciprocal_columns_inverse_bohr();
        let inputs = [[0., 0.], [1., 1.5], [1., 93. / 97.], [1., 1.]];
        for k in 1..=3 {
            let x = inputs[k - 1];
            let next = inputs[k];
            let f = [2. - 2. * x[0], 3. - 3. * x[1]];
            let values = |v: [f64; 2]| {
                plan.miller
                    .iter()
                    .map(|m| match m {
                        [0, 0, 0] => C::new(2. / volume, 0.),
                        [1, 0, 0] => C::new(1e-3 * v[0], 0.),
                        [0, 1, 0] => C::new(1e-3 * v[1], 0.),
                        _ => C::default(),
                    })
                    .collect::<Vec<_>>()
            };
            let y = values([x[0] + f[0], x[1] + f[1]]);
            let m = values(next);
            write_density(&dir.join(format!("k{k}-out")), b, &plan.miller, &y);
            let reverse_m: Vec<_> = plan.miller.iter().copied().rev().collect();
            let reverse_r: Vec<_> = m.iter().copied().rev().collect();
            let mix = dir.join(format!("k{k}-mix"));
            write_density(&mix, b, &reverse_m, &reverse_r);
            fs::write(
                mix.join("stdout"),
                format!(
                    "estimated scf accuracy = {:.16e} Ry\n",
                    8. * PI * volume * 1e-6 * (f[0] * f[0] + f[1] * f[1])
                ),
            )
            .unwrap();
        }
        let tol = Tolerance {
            absolute: 1e-14,
            relative: 1e-12,
        };
        execute(&input, &dir, 3, tol).unwrap();
        // Same production replay under simulated ranks; history remains owned.
        let parallel = dft_mpi::simulate(3, |comm| {
            let config = Config::parse(&text).unwrap();
            replay(&config, &dir, 3, tol, comm)
        });
        assert!(parallel.into_iter().all(|r| r.is_ok()));
        fs::write(
            dir.join("k2-mix/stdout"),
            "estimated scf accuracy = 99.000000000 Ry\n",
        )
        .unwrap();
        let error = execute(&input, &dir, 3, tol).unwrap_err();
        assert!(error.cause.contains("k=2 dr2"), "{}", error.cause);
        fs::remove_dir_all(&dir).unwrap();
    }
}
