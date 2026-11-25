//! Same-input density-map diagnostic. Teacher-forced QE densities never pass
//! through BAANDI's mixer. JSON and file access occur only between kernel calls.
use crate::refscf_check::{density, stdout};
use clap::{Args, ValueEnum};
use dft_cli::{
    config::Config,
    run::{self, RunError},
};
use dft_engine::{
    ScfProblem, SpeciesPotential,
    energy::EnergyAudit,
    input::{Atom, parse_upf},
    refscf::{
        density::DensityPlan,
        metric::{Metric, Serial},
        qe_reference::{QeDensity, difference},
    },
};
use dft_mpi::Collective;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::data::{Array3, Complex64 as C};
use serde_json::{Value, json};
use std::{
    fs,
    path::{Path, PathBuf},
};

// Four times the upper end of the user-supplied QE np8/np1 floor. This is an
// equality diagnostic, not a delivered-accuracy criterion for an SCF run.
const QE_FLOOR: f64 = 6.8e-17;
#[derive(Clone, Copy, Debug, ValueEnum)]
pub enum Solver {
    Legacy,
    Reference,
}
#[derive(Args)]
pub struct Arguments {
    pub input: PathBuf,
    #[arg(long)]
    pub qe: PathBuf,
    #[arg(long, default_value_t = 6)]
    pub steps: usize,
    /// Read each step's last printed ethr (Ry); also the default without --ethr.
    #[arg(long, conflicts_with = "ethr")]
    pub ethr_from_stdout: bool,
    /// Explicit occupied Ritz-change threshold in Ry, applied to every step.
    #[arg(long)]
    pub ethr: Option<f64>,
    #[arg(long, value_enum)]
    pub solver: Option<Solver>,
    #[arg(long, default_value_t = 4.0 * QE_FLOOR)]
    pub atol: f64,
    /// Restart atomic orbitals at every step instead of carrying BAANDI orbitals.
    #[arg(long)]
    pub cold_start: bool,
}
fn fail(s: impl std::fmt::Display) -> RunError {
    RunError {
        code: 4,
        cause: format!("refscf-map first failing term: {s}"),
    }
}

#[derive(Default, Debug)]
struct Printed {
    ethr: Vec<f64>,
    ethr_rounding: Vec<f64>,
    energies: Vec<(String, f64)>,
    bands_ev: Vec<f64>,
    occupations: Vec<f64>,
}
fn number(s: &str) -> Option<f64> {
    s.replace(['D', 'd'], "e")
        .parse::<f64>()
        .ok()
        .filter(|v| v.is_finite())
}
fn rounding(s: &str) -> f64 {
    let s = s.to_lowercase().replace('d', "e");
    let (m, e) = s
        .split_once('e')
        .map_or((s.as_str(), 0), |(m, e)| (m, e.parse::<i32>().unwrap_or(0)));
    let d = m.split_once('.').map_or(0, |(_, d)| d.len() as i32);
    0.5 * 10_f64.powi(e - d)
}
fn printed(text: &str, step: usize) -> Result<Printed, RunError> {
    let mut result = Printed::default();
    let mut current = 0;
    let mut reading = 0;
    for line in text.lines() {
        if let Some((_, tail)) = line.split_once("iteration #") {
            current = tail
                .split_whitespace()
                .next()
                .and_then(|s| s.parse().ok())
                .unwrap_or(0);
            reading = 0;
        }
        if current != step {
            continue;
        }
        if let Some((_, tail)) = line.split_once("ethr =") {
            let token = tail
                .split([',', ' '])
                .find(|s| !s.is_empty())
                .ok_or_else(|| fail("missing ethr"))?;
            let value = number(token)
                .filter(|v| *v > 0.)
                .ok_or_else(|| fail("invalid printed ethr"))?;
            result.ethr.push(value);
            result.ethr_rounding.push(rounding(token));
        }
        for (label, key) in [
            ("total energy", "total"),
            ("one-electron contribution", "one_electron"),
            ("hartree contribution", "hartree"),
            ("xc contribution", "xc"),
            ("ewald contribution", "ewald"),
            ("scf correction", "scf_correction"),
        ] {
            if line.contains(label)
                && let Some((_, tail)) = line.split_once('=')
            {
                let v = tail
                    .split_whitespace()
                    .next()
                    .and_then(number)
                    .ok_or_else(|| fail("invalid printed energy"))?;
                result.energies.retain(|(k, _)| k != key);
                result.energies.push((key.into(), v));
            }
        }
        if line.contains("bands (ev)") {
            reading = 1;
            result.bands_ev.clear();
            continue;
        }
        if line.contains("occupation numbers") {
            reading = 2;
            result.occupations.clear();
            continue;
        }
        if reading > 0 {
            let tokens: Vec<_> = line.split_whitespace().collect();
            if tokens.is_empty() {
                continue;
            }
            let values: Option<Vec<_>> = tokens.iter().map(|s| number(s)).collect();
            if let Some(values) = values {
                if reading == 1 {
                    result.bands_ev.extend(values);
                } else {
                    result.occupations.extend(values);
                }
            } else {
                reading = 0;
            }
        }
    }
    Ok(result)
}
fn diff(a: &[C], b: &[C], indices: &[usize], unit: &str) -> Value {
    let d = difference(a, b, indices);
    json!({"max":d.max,"rms":d.rms,"relative_max":d.relative_max,"relative_rms":d.relative_rms,
        "count":indices.len(),"unit":unit})
}
fn real_diff(a: &Array3<f64>, b: &Array3<f64>, unit: &str) -> Value {
    let max = a
        .iter()
        .zip(b)
        .map(|(a, b)| (a - b).abs())
        .fold(0., f64::max);
    let sum = a.iter().zip(b).map(|(a, b)| (a - b).powi(2)).sum::<f64>();
    json!({"max":max,"rms":(sum/a.len() as f64).sqrt(),"count":a.len(),"unit":unit})
}
fn inverse(
    plan: &mut DensityPlan,
    values: &[C],
    dims: [usize; 3],
) -> Result<Array3<f64>, RunError> {
    let mut planes = vec![0.; plan.plane_len()];
    plan.inverse(&Serial, values, &mut planes).map_err(fail)?;
    let mut full = Array3::zeros(dims);
    plan.write_serial_full(&planes, &mut full).map_err(fail)?;
    Ok(full)
}
fn forward(plan: &mut DensityPlan, full: &Array3<f64>) -> Result<Vec<C>, RunError> {
    let mut planes = vec![0.; plan.plane_len()];
    plan.planes_from_full(full, &mut planes).map_err(fail)?;
    let mut values = vec![C::default(); plan.miller.len()];
    plan.forward(&Serial, &planes, &mut values).map_err(fail)?;
    Ok(values)
}
struct Map {
    density: Array3<f64>,
    orbitals: ComplexMatrix,
    eigenvalues: Vec<f64>,
    residuals: Vec<f64>,
    occupations: Vec<f64>,
    energy: EnergyAudit,
    converged: Vec<bool>,
    iterations: usize,
    h_columns: usize,
}
fn solve(
    problem: &mut ScfProblem,
    n: &Array3<f64>,
    warm: Option<&ComplexMatrix>,
    ethr: f64,
    first: bool,
    selection: dft_cli::driver_select::Selection,
) -> Result<Map, RunError> {
    match selection.driver {
        dft_cli::driver_select::Driver::StrictOracle => {
            let starts = warm.map(std::slice::from_ref);
            let s = problem
                .step_kpoints_with_certification(n, starts, ethr / 2., selection.strict())
                .map_err(fail)?;
            Ok(Map {
                density: s.density,
                orbitals: s.solution.eigenvectors,
                eigenvalues: s.solution.eigenvalues,
                residuals: s.solution.residuals,
                occupations: problem.fixed_map_occupations().to_vec(),
                energy: s.energy,
                converged: s.solution.converged,
                iterations: s.solution.iterations,
                h_columns: s.solution.operator_applications,
            })
        }
        dft_cli::driver_select::Driver::Reference => {
            let s = problem
                .reference_fixed_map(n, warm, ethr / 2., first)
                .map_err(fail)?;
            Ok(Map {
                density: s.density,
                orbitals: s.orbitals,
                eigenvalues: s.eigenvalues,
                residuals: s.residuals,
                occupations: s.occupations,
                energy: s.energy,
                converged: s.converged,
                iterations: s.record.iterations.len(),
                h_columns: s.record.h_columns,
            })
        }
    }
}

pub fn execute(args: &Arguments) -> Result<(), RunError> {
    if args.steps == 0
        || !args.atol.is_finite()
        || args.atol < 0.
        || args.ethr.is_some_and(|e| !e.is_finite() || e <= 0.)
    {
        return Err(fail("invalid steps/threshold/tolerance"));
    }
    let text = fs::read_to_string(&args.input).map_err(fail)?;
    let config = Config::parse(&text).map_err(fail)?;
    let points = config.resolved_points().map_err(fail)?;
    if config.collinear()
        || config.finite_temperature().map_err(fail)?.is_some()
        || points.len() != 1
        || points[0].fractional != [0.; 3]
        || config.symmetry.mode != dft_cli::config::SymmetryMode::Off
    {
        return Err(fail(
            "requires unpolarized, fixed occupations, Gamma, symmetry off",
        ));
    }
    let profile = run::reference_profile(&config)?;
    if profile.beta == 1. {
        return Err(fail("rho_in(1) recovery requires mixing_beta != 1"));
    }
    let prepared = config.prepare().map_err(fail)?;
    let base = args.input.parent().unwrap_or(Path::new("."));
    let mut pps = Vec::new();
    let mut provenance = Vec::new();
    for (name, path) in &config.pseudopotentials {
        let data = fs::read(base.join(path)).map_err(fail)?;
        provenance.push(json!({"species":name,"path":base.join(path),"sha256":run::sha256(&data)}));
        pps.push((
            name.clone(),
            parse_upf(std::str::from_utf8(&data).map_err(fail)?).map_err(fail)?,
        ));
    }
    let refs: Vec<_> = pps
        .iter()
        .map(|(label, potential)| SpeciesPotential { label, potential })
        .collect();
    let atoms: Vec<_> = config
        .system
        .atoms
        .iter()
        .map(|a| Atom::new(&a.species, a.fractional))
        .collect();
    let selection = dft_cli::driver_select::select(
        &config,
        dft_cli::driver_select::Mode::Classic,
        dft_cli::driver_select::Command::RefscfMap {
            reference: args
                .solver
                .map(|solver| matches!(solver, Solver::Reference)),
        },
        &pps.iter().map(|(_, p)| p).collect::<Vec<_>>(),
    );
    selection.announce();
    let mut opts = run::options_for(&config, selection)?;
    // This diagnostic uses the oracle only for the shared map data.
    opts.standard = None;
    opts.verification = dft_engine::Verification::Off;
    opts.gamma_real = dft_engine::GammaReal::On;
    opts.start_wavefunctions = dft_engine::StartWavefunctions::Atomic;
    let mut problem = ScfProblem::new(&prepared.cell, &atoms, &refs, opts).map_err(fail)?;
    if let Some(floor) = problem.core_density().map(|_| 1e-10) {
        problem.set_density_floor(floor).map_err(fail)?;
    }
    let dims = problem.grid_dimensions();
    let volume = prepared.cell.volume_bohr3();
    // Replicated diagnostic transforms, even under MPI. Only the inner solve
    // uses distributed H/rows; no timing or distributed-memory claim is made.
    let mut plan = DensityPlan::new(
        &Serial,
        &prepared.cell,
        dims,
        2. * config.basis.ecutwfc_ry,
        2. * config.basis.ecutwfc_ry,
    )
    .map_err(fail)?;
    let align = |d: QeDensity| -> Result<Vec<C>, RunError> {
        d.validate_cell(prepared.cell.reciprocal_columns_inverse_bohr())
            .map_err(fail)?;
        d.aligned(&plan.miller).map_err(fail)
    };
    // Load/validate all snapshots before entering the maps. Never drop a G row.
    let mut expected = Vec::new();
    let mut mixed = Vec::new();
    let mut logs = Vec::new();
    for k in 1..=args.steps {
        expected.push(align(
            density(&args.qe.join(format!("k{k}-out"))).map_err(fail)?,
        )?);
        if k < args.steps || k == 1 {
            mixed.push(align(
                density(&args.qe.join(format!("k{k}-mix"))).map_err(fail)?,
            )?);
        }
        let log = stdout(&args.qe.join(format!("k{k}-out")), &args.qe, k).map_err(fail)?;
        let p = printed(log.as_deref().unwrap_or(""), k)?;
        if args.ethr.is_none() && p.ethr.is_empty() {
            return Err(fail(format!(
                "k={k} missing ethr in output stdout; supply --ethr Ry"
            )));
        }
        logs.push(p);
    }
    let initial: Vec<_> = mixed[0]
        .iter()
        .zip(&expected[0])
        .map(|(m, y)| (*m - profile.beta * y) / (1. - profile.beta))
        .collect();
    let all: Vec<_> = (0..plan.miller.len()).collect();
    let zero = plan
        .miller
        .iter()
        .position(|m| *m == [0; 3])
        .ok_or_else(|| fail("missing G=0"))?;
    let mut shell_g2 = plan
        .g2
        .iter()
        .copied()
        .filter(|g| *g > 0.)
        .collect::<Vec<_>>();
    shell_g2.sort_by(f64::total_cmp);
    shell_g2.dedup_by(|a, b| (*a - *b).abs() <= 1e-10 * a.max(*b).max(1.));
    shell_g2.truncate(10);
    let groups: Vec<Vec<_>> = shell_g2
        .iter()
        .map(|g| {
            plan.g2
                .iter()
                .enumerate()
                .filter_map(|(i, v)| ((v - g).abs() <= 1e-10 * g.max(1.)).then_some(i))
                .collect()
        })
        .collect();
    let mut metric = Metric::new(volume, &plan.g2, &plan.smooth, 1, 350).map_err(fail)?;
    let rank = dft_engine::parallel::SharedComm::active().map_or(0, |c| c.rank());
    let ranks = dft_engine::parallel::SharedComm::active().map_or(1, |c| c.size());
    if rank == 0 {
        println!(
            "{}",
            json!({"schema":"refscf-map-1","solver":if selection.reference() { "reference" } else { "legacy" },"scf":selection.json(),
        "input_sha256":run::sha256(text.as_bytes()),"pseudopotentials":provenance,"qe_directory":args.qe,
        "ranks":ranks,"grid":dims,"volume_bohr3":volume,"g_vectors":all.len(),"electrons":problem.electron_count(),
        "density_atol_e_bohr3":args.atol,"qe_rank_floor_max_e_bohr3":QE_FLOOR,
        "floor_provenance":"user/dev08 relay; QE np8 vs np1 1.7e-17..6.8e-17 e/bohr^3",
        "orbital_start":if args.cold_start {"BAANDI atomic at every step"} else {"BAANDI atomic at k1, preceding BAANDI orbitals thereafter"},
        "density_start":"QE input at EVERY step; never BAANDI mixed/output density",
        "limits":["QE wavefunctions are not supplied; loose solves can differ even at the same density and ethr",
            "stdout ethr is rounded; its printed precision is reported; --ethr accepts an exact value",
            "diagnostic density transforms and output orbitals are replicated; not a scaling benchmark"]})
        );
    }
    let mut warm = None;
    let mut first_failure = None;
    for k in 1..=args.steps {
        let x = if k == 1 { &initial } else { &mixed[k - 2] };
        let n = inverse(&mut plan, x, dims)?;
        let qout = inverse(&mut plan, &expected[k - 1], dims)?;
        let roundtrip = forward(&mut plan, &n)?;
        let roundtrip_diff = diff(&roundtrip, x, &all, "e/bohr^3");
        let p = &logs[k - 1];
        let ethr = args
            .ethr
            .unwrap_or_else(|| *p.ethr.last().expect("validated"));
        let (ih, _, _) = problem.fixed_map_fields(&n).map_err(fail)?;
        let ih_g = forward(&mut plan, &ih.potential.mapv(|v| v.re))?;
        let analytic_h: Vec<_> = x
            .iter()
            .zip(&plan.g2)
            .map(|(v, g)| {
                if *g == 0. {
                    C::default()
                } else {
                    4. * std::f64::consts::PI * v / g
                }
            })
            .collect();
        let hartree_check = diff(&ih_g, &analytic_h, &all, "Ha");
        // Replay every printed first-step redo on the identical input, preserving
        // the first returned BAANDI block; explicit --ethr performs one solve.
        let thresholds = if k == 1 && args.ethr.is_none() {
            p.ethr.clone()
        } else {
            vec![ethr]
        };
        if args.cold_start {
            warm = None;
        }
        let mut map = None;
        let (mut iterations, mut columns) = (0, 0);
        for t in thresholds.iter().copied() {
            let result = solve(&mut problem, &n, warm.as_ref(), t, k == 1, selection)?;
            iterations += result.iterations;
            columns += result.h_columns;
            warm = Some(result.orbitals.clone());
            map = Some(result);
        }
        let result = map.expect("nonempty thresholds");
        let out = forward(&mut plan, &result.density)?;
        let delta: Vec<_> = out
            .iter()
            .zip(&expected[k - 1])
            .map(|(a, b)| a - b)
            .collect();
        let map_delta: Vec<_> = out.iter().zip(x).map(|(a, b)| a - b).collect();
        let qe_delta: Vec<_> = expected[k - 1].iter().zip(x).map(|(a, b)| a - b).collect();
        let all_g = diff(&out, &expected[k - 1], &all, "e/bohr^3");
        let shells: Vec<_> = shell_g2.iter().zip(&groups).map(|(g,idx)|json!({"g2_bohr_minus2":g,"difference":diff(&out,&expected[k-1],idx,"e/bohr^3")})).collect();
        let charge = [
            x[zero].re * volume,
            out[zero].re * volume,
            expected[k - 1][zero].re * volume,
        ];
        let (oh, ox, _) = problem.fixed_map_fields(&result.density).map_err(fail)?;
        let (qh, qx, _) = problem.fixed_map_fields(&qout).map_err(fail)?;
        let energy = energy_comparison(&result.energy, p, qh.energy.0, qx.energy.0);
        let ev = result
            .eigenvalues
            .iter()
            .map(|e| dft_units::from_hartree(*e, dft_units::EnergyUnit::ElectronVolt))
            .collect::<Vec<_>>();
        let bands = if p.bands_ev.is_empty() {
            json!({"available":false,"reason":"no bands in this step stdout","baandi_ev":ev})
        } else {
            json!({"available":true,"baandi_ev":ev,"qe_ev":p.bands_ev,"difference_ev":ev.iter().zip(&p.bands_ev).map(|(a,b)|a-b).collect::<Vec<_>>(),"qe_print_rounding_ev":0.00005})
        };
        let qe_occ: Vec<_> = p.occupations.iter().map(|f| 2. * f).collect();
        let occupations_match = qe_occ.is_empty()
            || (qe_occ.len() == result.occupations.len()
                && qe_occ
                    .iter()
                    .zip(&result.occupations)
                    .all(|(a, b)| (a - b).abs() < 1e-12));
        let failure = if (charge[0] - problem.electron_count()).abs() > 1e-10 {
            Some("input total charge")
        } else if roundtrip_diff["max"].as_f64().unwrap() > args.atol {
            Some("density import roundtrip")
        } else if (charge[1] - charge[2]).abs() > 1e-10 {
            Some("output total charge")
        } else if hartree_check["max"].as_f64().unwrap() > 1e-12 {
            Some("input Hartree potential vs analytic")
        } else if !occupations_match {
            Some("occupations")
        } else if !p.bands_ev.is_empty() && p.bands_ev.len() != ev.len() {
            Some("band count")
        } else if result.converged.iter().any(|c| !*c) {
            Some("inner solve convergence")
        } else if all_g["max"].as_f64().unwrap() > args.atol {
            Some("rho_out all G")
        } else {
            None
        };
        if first_failure.is_none() {
            first_failure = failure.map(|f| format!("k={k} {f}"));
        }
        if rank == 0 {
            println!(
                "{}",
                json!({"iteration":k,"ethr_ry":ethr,"ethr_ha":ethr/2.,
            "ethr_source":if args.ethr.is_some(){"explicit"}else{"stdout (rounded)"},
            "qe_printed_ethr_ry":p.ethr,"qe_ethr_rounding_ry":p.ethr_rounding,"solve_thresholds_ry":thresholds,
            "first_failing_term":failure,"all_g":all_g,"g0":diff(&out,&expected[k-1],&[zero],"e/bohr^3"),
            "real_space":real_diff(&result.density,&qout,"e/bohr^3"),"input_roundtrip":roundtrip_diff,
            "total_charge":{"input_e":charge[0],"actual_e":charge[1],"expected_e":charge[2],"difference_e":charge[1]-charge[2]},
            "lowest_10_nonzero_shells":shells,"density_distance_ry":metric.norm(&Serial,&delta).map_err(fail)?,
            "density_distance_ha":metric.norm(&Serial,&delta).map_err(fail)?/2.,
            "baandi_dr2_ry":metric.norm(&Serial,&map_delta).map_err(fail)?,"qe_dr2_ry":metric.norm(&Serial,&qe_delta).map_err(fail)?,
            "potential":{"input_hartree_vs_analytic":hartree_check,"qe_input_xc_difference":null,"qe_input_local_difference":null,
                "reason":"QE potential grids absent; density-only Hartree formula is independently checkable",
                "output_density_response":{"meaning":"BAANDI V[rho_out_B] minus BAANDI V[rho_out_QE]; not cross-code potentials",
                    "hartree":real_diff(&oh.potential.mapv(|z|z.re),&qh.potential.mapv(|z|z.re),"Ha"),
                    "xc":real_diff(&ox.potential,&qx.potential,"Ha"),"local":{"max":0.,"rms":0.,"unit":"Ha","reason":"same density-independent ionic potential"}}},
            "energy":energy,"bands":bands,"occupations":{"baandi_e":result.occupations,"qe_e":qe_occ,"qe_spin_factor":2.,"available":!p.occupations.is_empty()},
            "inner_solve":{"iterations":iterations,"h_columns":columns,"converged":result.converged,"residuals_ha":result.residuals,"certified":false}})
            );
        }
    }
    if rank == 0 {
        println!(
            "{}",
            json!({"summary":"refscf-map","steps":args.steps,"first_failing_term":first_failure,"agreement":first_failure.is_none()})
        );
    }
    if let Some(f) = first_failure {
        Err(fail(f))
    } else {
        Ok(())
    }
}
fn energy_comparison(e: &EnergyAudit, p: &Printed, qe_density_h: f64, qe_density_x: f64) -> Value {
    let terms: Vec<_> = p
        .energies
        .iter()
        .map(|(key, q)| {
            let b = match key.as_str() {
                "total" => e.band_form,
                "one_electron" => e.band_sum - e.double_counting,
                "hartree" => e.hartree,
                "xc" => e.xc,
                "ewald" => e.ewald,
                "scf_correction" => 0.,
                _ => unreachable!(),
            } * 2.;
            json!({"term":key,"baandi_ry":b,"qe_ry":q,"difference_ry":b-q})
        })
        .collect();
    let qe_energy = |name: &str| p.energies.iter().find(|(k, _)| k == name).map(|(_, v)| *v);
    json!({"terms":terms,"available":!p.energies.is_empty(),"convention":"k-out forced-stop energy: band sum - integral rho_out*(VH+Vxc)[rho_in] + EH[rho_out] + Exc[rho_out] + Ewald; descf=0",
        "baandi_direct_orbital_total_ha":e.total_hartree(),"baandi_band_form_ha":e.band_form,
        "qe_output_density_in_baandi_kernels":{"hartree_ry":2.*qe_density_h,"xc_ry":2.*qe_density_x,
            "hartree_minus_qe_printed_ry":qe_energy("hartree").map(|v|2.*qe_density_h-v),
            "xc_minus_qe_printed_ry":qe_energy("xc").map(|v|2.*qe_density_x-v)},
        "status":"diagnostic differences; no energy tolerance inferred from a density floor"})
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn step_selection_redo_units_and_missing_terms() {
        let text = "iteration # 1\n ethr = 1.00D-02, avg # of iterations = 3.0\n ethr = 2.50d-04, avg # of iterations = 2.0\niteration # 2\n ethr = 1.37E-03, avg # of iterations = 1.0\n bands (ev):\n\n -1.25 2.5\n occupation numbers\n 1.0 0.0\n highest occupied, lowest unoccupied level (ev): 0 0\n! total energy = -3.5 Ry\n";
        let p = printed(text, 1).unwrap();
        assert_eq!(p.ethr, vec![1e-2, 2.5e-4]);
        assert!(p.energies.is_empty());
        let p = printed(text, 2).unwrap();
        assert_eq!(p.ethr, vec![1.37e-3]);
        assert!((p.ethr_rounding[0] - 5e-6).abs() < 1e-20);
        assert_eq!(p.bands_ev, vec![-1.25, 2.5]);
        assert_eq!(p.occupations, vec![1., 0.]);
        assert_eq!(p.energies, vec![("total".into(), -3.5)]);
        assert!(printed(text, 3).unwrap().ethr.is_empty());
        assert!(printed("iteration # 1\n ethr = NaN,", 1).is_err());
        assert!(printed("iteration # 1\n ethr = -0.1,", 1).is_err());
    }
    #[test]
    fn qe_stdout_is_read_from_the_output_snapshot() {
        let root =
            Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/reference/qe/i25-si8/si8-np1");
        for (k, ethr) in [1e-2, 1.37e-3, 6.11e-5].into_iter().enumerate() {
            let k = k + 1;
            let log = stdout(&root.join(format!("k{k}-out")), &root, k)
                .unwrap()
                .unwrap();
            let p = printed(&log, k).unwrap();
            assert_eq!(p.ethr, vec![ethr]);
            assert_eq!(p.bands_ev.len(), 20);
            assert_eq!(p.occupations.iter().sum::<f64>(), 16.);
            assert_eq!(p.energies.len(), 5);
        }
    }
}
