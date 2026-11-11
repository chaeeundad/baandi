#[path = "finite_run.rs"]
mod finite_run;
// File/console boundary for the sequential engine. No I/O occurs in SCF kernels.
use crate::{config::Config, grid_report};
#[path = "spin_run.rs"]
mod spin_run;
#[path = "volumetric_out.rs"]
pub mod volumetric_out;
use dft_engine::{
    IterationRecord, Mixing, PhaseTimings, ScfError, ScfOptions, ScfProblem, SpeciesPotential,
    energy::EnergyAudit,
    input::{Atom, parse_upf},
};
use dft_units::{EnergyUnit, from_hartree};
use serde_json::{Value, json};
use sha2::{Digest, Sha256};
use std::{
    fmt::Write as _,
    fs, io,
    path::{Path, PathBuf},
    time::{Instant, SystemTime, UNIX_EPOCH},
};

#[derive(Debug, thiserror::Error)]
#[error("{cause}")]
pub struct RunError {
    pub code: u8,
    pub cause: String,
}
fn error(code: u8, cause: impl std::fmt::Display) -> RunError {
    RunError {
        code,
        cause: cause.to_string(),
    }
}
impl From<io::Error> for RunError {
    fn from(e: io::Error) -> Self {
        error(5, e)
    }
}
pub(crate) fn engine_error(e: ScfError) -> RunError {
    let code = match e {
        ScfError::InvalidInput(_)
        | ScfError::OddElectronCount(_)
        | ScfError::InsufficientBands { .. } => 1,
        ScfError::NotConverged { .. } => 3,
        _ => 4,
    };
    error(code, e)
}
pub fn sha256(bytes: &[u8]) -> String {
    format!("{:x}", Sha256::digest(bytes))
}

/// Resolve only user-exposed controls, preserving all other engine defaults.
pub fn options(config: &Config) -> Result<ScfOptions, RunError> {
    let mut o = ScfOptions {
        compat_mode: if config.compat.reference() {
            dft_engine::CompatMode::Reference
        } else {
            dft_engine::CompatMode::None
        },
        crystal_density_group: config.symmetry.density_group
            == crate::config::DensityGroup::Crystal,
        functional: config.xc.functional.engine(),
        local_mesh: config
            .pseudo
            .local_mesh()
            .map_err(|e| error(1, e.to_string()))?,
        radial_even_rule: config
            .pseudo
            .radial_even_rule()
            .map_err(|e| error(1, e.to_string()))?,
        pw92_prefactor: config.xc.prefactor().map_err(|e| error(1, e.to_string()))?,
        allow_functional_override: config.xc.allow_functional_override,
        cutoff_hartree: config.basis.ecutwfc_ry / 2.,
        grid_override: if config.symmetry.mode == crate::config::SymmetryMode::Auto {
            Some(config.prepare().map_err(|e| error(1, e))?.grid.dims())
        } else {
            config.basis.fft_grid
        },
        ..ScfOptions::default()
    };
    if let Some(s) = &config.scf {
        o.standard = s.standard_profile().map_err(|e| error(1, e.to_string()))?;
        o.standard_stall_guard = s
            .standard_stall_guard()
            .map_err(|e| error(1, e.to_string()))?;
        o.standard_density_forcing = s
            .standard_density_forcing()
            .map_err(|e| error(1, e.to_string()))?;
        o.hartree_mixing = s.hartree_mixing().map_err(|e| error(1, e.to_string()))?;
        o.spectral_guard = s.spectral_guard().map_err(|e| error(1, e.to_string()))?;
        o.preconditioner = s.preconditioner().map_err(|e| error(1, e.to_string()))?;
        o.davidson_capacity_per_band = s
            .davidson_capacity_per_band()
            .map_err(|e| error(1, e.to_string()))?;
        o.start_wavefunctions = s
            .start_wavefunctions()
            .map_err(|e| error(1, e.to_string()))?;
        o.start_density = s.start_density().map_err(|e| error(1, e.to_string()))?;
        o.gamma_real = s.gamma_real().map_err(|e| error(1, e.to_string()))?;
        o.locking = s.locking().map_err(|e| error(1, e.to_string()))?;
        o.threads = s.threads.unwrap_or(1);
        if let Some(n) = s.max_iterations {
            o.max_iterations =
                usize::try_from(n).map_err(|_| error(1, "max_iterations must be positive"))?;
        }
        o.bands = s.nbands;
        if let Some(n) = s.extra_bands {
            o.extra_bands = n;
        }
        if let Some(x) = s.energy_tolerance_ha {
            o.energy_tolerance = x;
        }
        if let Some(x) = s.density_residual_tolerance {
            o.density_tolerance = x;
        }
        if let Some(x) = s.eigen_tolerance {
            o.eigen_tolerance = x;
        }
        if let Some(x) = s.seed {
            o.seed = x;
        }
        o.mixing = s.mixing_strategy().map_err(|e| error(1, e.to_string()))?;
    }
    if o.max_iterations == 0
        || [o.energy_tolerance, o.density_tolerance, o.eigen_tolerance]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.)
    {
        return Err(error(
            1,
            "iteration count and tolerances must be positive and finite",
        ));
    }
    Ok(o)
}
pub(crate) fn reference_metadata(o: &ScfOptions) -> Option<Value> {
    (o.compat_mode == dft_engine::CompatMode::Reference).then(|| {
        json!({
            "compat_mode": o.compat_mode.name(),
            "pw92_prefactor": o.pw92_prefactor.name(),
            "local_mesh": o.local_mesh.name(),
            "radial_even_rule": o.radial_even_rule.name(),
            "density_group": if o.crystal_density_group { "crystal" } else { "mesh" },
            "scope": crate::config::REFERENCE_SCOPE,
        })
    })
}
fn options_json(o: &ScfOptions) -> Value {
    let mixing = match o.mixing {
        Mixing::Linear { beta } => json!({"method":"linear","beta":beta}),
        Mixing::Pulay { beta, history } => json!({"method":"pulay","beta":beta,"history":history}),
        Mixing::KerkerPulay {
            beta,
            history,
            q0_inverse_bohr,
        } => {
            json!({"method":"kerker-pulay","beta":beta,"history":history,"q0_inverse_bohr":q0_inverse_bohr})
        }
    };
    let mut result = json!({"spectral_guard":o.spectral_guard.name(),"threads":o.threads,"cutoff_hartree":o.cutoff_hartree,"nbands":o.bands,"extra_bands":o.extra_bands,"grid_override":o.grid_override,"mixing":mixing,
        "density_residual_tolerance":o.density_tolerance,"energy_tolerance_ha":o.energy_tolerance,"eigen_tolerance":o.eigen_tolerance,
        "residual_eigen_schedule":o.residual_eigen_schedule,"free_electron_initial_guess":o.free_electron_initial_guess,"initial_eigen_tolerance":o.initial_eigen_tolerance,"eigen_tolerance_decay":o.eigen_tolerance_decay,"max_iterations":o.max_iterations,
        "preconditioner":o.preconditioner.name(),"davidson_max_iterations":o.davidson_max_iterations,"davidson_subspace":o.davidson_subspace,"davidson_subspace_min_per_band":o.davidson_capacity_per_band,"davidson_capacity_per_band":o.davidson_capacity_per_band,"seed":o.seed,
        "divergence_factor":o.divergence_factor,"divergence_patience":o.divergence_patience});
    if let Some(metadata) = reference_metadata(o) {
        result["compat_mode"] = json!(o.compat_mode.name());
        result["compat"] = metadata;
    }
    if o.local_mesh != dft_engine::LocalMesh::Baandi
        || o.radial_even_rule != dft_engine::RadialEvenRule::Simpson3_8
    {
        result["pseudo"] =
            json!({"local_mesh":o.local_mesh.name(),"radial_even_rule":o.radial_even_rule.name()});
    }
    if o.pw92_prefactor == dft_engine::Pw92Prefactor::Original {
        result["pw92_prefactor"] = json!(o.pw92_prefactor.name());
    }
    // Always recorded: the default is atomic (ADR 0042), so omitting the default
    // would hide the non-default free_electron.
    result["start_wavefunctions"] = json!(o.start_wavefunctions.name());
    if o.start_density != dft_engine::StartDensity::Uniform {
        result["start_density"] = json!(o.start_density.name());
    }
    // Always recorded: the default is Auto (ADR 0040), so omitting the default
    // would hide the non-default `off`.
    result["gamma_real"] = json!(o.gamma_real.name());
    if o.locking == dft_engine::Locking::Auto {
        result["locking"] = json!("auto");
    }
    if let Some(p) = o.standard {
        if o.standard_stall_guard {
            result["standard_stall_guard"] = json!("on");
        }
        if o.standard_density_forcing {
            result["standard_density_forcing"] = json!("on");
        }
        result["convergence_profile"] = json!("standard");
        result["standard_threshold_ha"] = json!(p.threshold_ha);
        result["eigensolver_stopping_rule"] = json!("ritz_value_change");
        result["standard_eigen_cap_ha"] = json!(dft_engine::StandardProfile::CAP_HA);
        result["standard_eigen_floor_ha"] = json!(dft_engine::StandardProfile::FLOOR_HA);
        result["empty_threshold_rule"] =
            json!("max(5 * occupied_threshold_ha, 5e-6 Ha); first empty required");
    }
    if o.hartree_mixing || o.standard.is_some() {
        result["mixing_metric"] = json!(if o.hartree_mixing { "hartree" } else { "l2" });
    }
    // Switch manifest (work order A1): every BAANDI_* switch with its effective
    // value and source, so a result shows which combination produced it.
    result["switches"] = switches_json();
    result
}
/// `scf.switches`: name -> {value, source, default, allowed, meaning}. Values are
/// strings, numbers, booleans or null (null: unset override / two-pass ORTHO /
/// no growth reset); source is `default`, `env`, or `invalid -> <fallback>`.
fn switches_json() -> Value {
    use dft_numerics::switches::{self, Value as V};
    let mut m = serde_json::Map::new();
    for e in switches::manifest() {
        let value = match e.value {
            V::Choice(s) => json!(s),
            V::Number(x) => json!(x),
            V::Raw(s) => json!(s),
            V::Present(b) => json!(b),
        };
        m.insert(
            e.switch.name.into(),
            json!({"value":value,"source":e.source.as_str(),"default":e.switch.default,
                "allowed":e.switch.allowed,"meaning":e.switch.meaning}),
        );
    }
    Value::Object(m)
}
fn energy_json(e: &EnergyAudit) -> Value {
    let mut m = serde_json::Map::new();
    for (name, x) in [
        ("kinetic", e.kinetic),
        ("local", e.local),
        ("local_zero_mode", e.local_zero_mode),
        ("nonlocal", e.nonlocal),
        ("hartree", e.hartree),
        ("xc", e.xc),
        ("ewald", e.ewald),
        ("band_sum", e.band_sum),
        ("double_counting", e.double_counting),
        ("band_form", e.band_form),
        ("harris_foulkes", e.harris_foulkes),
        ("functional_difference_bound", e.functional_difference_bound),
        ("total", e.total_hartree()),
    ] {
        m.insert(name.into(), json!({"ha":x,"ry":2.*x}));
    }
    if let Some(floor) = &e.xc_floor {
        m.insert(
            "xc_floor".into(),
            json!({"density_floor_electrons_bohr3":floor.density_floor,
            "input_negative_charge_electrons":floor.input_negative_charge,
            "output_negative_charge_electrons":floor.output_negative_charge}),
        );
    }
    if let Some(core) = &e.nlcc {
        m.insert("nlcc".into(), json!({"core_electrons":core.core_electrons,"core_potential_integral_ha":core.core_potential_integral_ha,"valence_potential_integral_ha":core.valence_potential_integral_ha,"energy_convention":"xc_includes_core; Hartree_and_band_double_counting_valence_only","force_stress":"unavailable"}));
    }
    Value::Object(m)
}
fn phase_json(t: &PhaseTimings) -> Value {
    json!({"potential_s":t.potential.as_secs_f64(),"diagonalization_s":t.diagonalization.as_secs_f64(),"density_s":t.density.as_secs_f64(),"total_energy_s":t.energy.as_secs_f64(),"mixing_s":t.mixing.as_secs_f64()})
}
fn phase_sum(t: &PhaseTimings) -> f64 {
    (t.potential + t.diagonalization + t.density + t.energy + t.mixing).as_secs_f64()
}
fn history_json(r: &IterationRecord) -> Value {
    let mut result = json!({"iteration":r.iteration,"energies":energy_json(&r.energy),"delta_energy_ha":r.energy_change,"density_residual":r.density_residual,
        "max_eigen_residual_ha":r.max_eigen_residual,"requested_eigen_tolerance_ha":r.requested_eigen_tolerance,"davidson_iterations":r.davidson_iterations,"operator_applications":r.operator_applications,"davidson_iterations_per_k":r.davidson_iterations_per_k,"operator_applications_per_k":r.operator_applications_per_k,
        "relative_charge_error":r.charge_error,"relative_input_charge_error":r.input_charge_error,"min_density_electrons_per_bohr3":r.min_density,
        "mixing_fraction":r.mixing_fraction,"mixer_reset":r.mixer_reset,"density_contraction":r.density_contraction,"mixing_history_used":r.mixing_history_used,"mixing_resets":r.mixing_resets,"mixing_raw_recovery":r.mixing_raw_recovery,"timings":phase_json(&r.timings),
        "fixed_point_certificate":r.fixed_point_certificate.as_ref().map(|c| json!({"strict_residuals_per_k_ha":c.strict_residuals_per_k,"buffer_residuals_per_k_ha":c.buffer_residuals_per_k,"density_residual":c.density_residual,"occupied_eigenvalue_change_ha":c.occupied_eigenvalue_change,"davidson_iterations":c.davidson_iterations,"operator_applications":c.operator_applications,"davidson_iterations_per_k":c.davidson_iterations_per_k,"operator_applications_per_k":c.operator_applications_per_k,"timings":phase_json(&c.timings)}))});
    if let Some(e) = r.standard_error_ha {
        result["standard_error_ha"] = json!(e);
        result["standard_first_iteration_redo"] = json!(r.standard_redo);
    }
    result
}
/// Gregorian UTC conversion from nonnegative Unix seconds, including leap years.
pub fn timestamp(seconds: u64) -> String {
    let days = (seconds / 86400) as i64;
    let z = days + 719468;
    let era = z / 146097;
    let doe = z - era * 146097;
    let yoe = (doe - doe / 1460 + doe / 36524 - doe / 146096) / 365;
    let mut year = yoe + era * 400;
    let doy = doe - (365 * yoe + yoe / 4 - yoe / 100);
    let mp = (5 * doy + 2) / 153;
    let day = doy - (153 * mp + 2) / 5 + 1;
    let month = mp + if mp < 10 { 3 } else { -9 };
    year += i64::from(month <= 2);
    let s = seconds % 86400;
    format!(
        "{year:04}-{month:02}-{day:02}T{:02}:{:02}:{:02}Z",
        s / 3600,
        s / 60 % 60,
        s % 60
    )
}

pub struct RunReport {
    pub console: String,
    pub converged: bool,
}

// Canonicalize existing destinations (including symlinks) before any writes.
// Unix inode checks additionally protect hard-linked input files.
fn protect_outputs(dir: &Path, inputs: &[PathBuf]) -> Result<(), RunError> {
    for name in [
        "results.json",
        "density.json",
        "density.bin",
        "density-up.cube",
        "density-down.cube",
        "density-total.cube",
        "magnetization.cube",
    ] {
        let dest = dir.join(name);
        if dest.exists() {
            let canonical = fs::canonicalize(&dest)?;
            for input in inputs {
                let mut same = canonical == *input;
                #[cfg(unix)]
                {
                    use std::os::unix::fs::MetadataExt;
                    let a = fs::metadata(&dest)?;
                    let b = fs::metadata(input)?;
                    same |= a.dev() == b.dev() && a.ino() == b.ino();
                }
                if same {
                    return Err(error(
                        5,
                        format!("refusing to overwrite input {}", dest.display()),
                    ));
                }
            }
        }
    }
    Ok(())
}
fn write_json(path: &Path, value: &Value) -> Result<(), RunError> {
    let data = serde_json::to_vec_pretty(value).map_err(|e| error(5, e))?;
    fs::write(path, data)?;
    Ok(())
}

pub fn execute(input: &Path, output_override: Option<&Path>) -> Result<RunReport, RunError> {
    execute_mode(
        input,
        output_override,
        &crate::telemetry::RunMode::default(),
    )
}

pub fn execute_mode(
    input: &Path,
    output_override: Option<&Path>,
    mode: &crate::telemetry::RunMode,
) -> Result<RunReport, RunError> {
    mode.validate()?;
    let start = Instant::now();
    let bytes = fs::read(input).map_err(|e| error(5, format!("{}: {e}", input.display())))?;
    let text = std::str::from_utf8(&bytes).map_err(|e| error(1, e))?;
    let config = Config::parse(text).map_err(|e| error(1, e))?;
    if (config.collinear()
        || config
            .finite_temperature()
            .map_err(|e| error(1, e))?
            .is_some())
        && (mode.observe || mode.agent)
    {
        return Err(error(
            1,
            "spin.execution: runtime control is not spin-aware; use Classic execution",
        ));
    }
    if config.symmetry.mode == crate::config::SymmetryMode::Auto && (mode.observe || mode.agent) {
        return Err(error(
            1,
            "symmetry.execution: checkpoint identity does not carry symmetry yet; use Classic execution",
        ));
    }
    let preflight = crate::config::preflight::check(&config, input).map_err(|e| error(1, e))?;
    if mode.events.as_deref() == Some("-") {
        eprintln!("{}", preflight.summary);
        for warning in &preflight.warnings {
            eprintln!("warning: {warning}");
        }
        for item in &preflight.overrides {
            eprintln!("override: {item}");
        }
    } else {
        preflight.print();
    }
    let prepared = config.prepare().map_err(|e| error(1, e))?;
    let opts = options(&config)?;
    let mut potentials = Vec::new();
    let mut species = Vec::new();
    let mut inputs = vec![fs::canonicalize(input)?];
    for (name, path) in &config.pseudopotentials {
        let path = input.parent().unwrap_or(Path::new(".")).join(path);
        let data = fs::read(&path).map_err(|e| error(5, format!("{}: {e}", path.display())))?;
        let p = parse_upf(std::str::from_utf8(&data).map_err(|e| error(2, e))?)
            .map_err(|e| error(2, format!("{}: {e}", path.display())))?;
        if !opts.functional.matches_upf(p.functional()) && !opts.allow_functional_override {
            if opts.functional == dft_engine::Functional::LdaPz {
                return Err(error(
                    2,
                    format!(
                        "unsupported UPF functional {:?}: requires SLA PZ or SLA PZ NOGX NOGC",
                        p.functional()
                    ),
                ));
            }
            return Err(error(
                2,
                format!(
                    "UPF functional {:?} mismatches requested XC; xc.allow_functional_override = true required",
                    p.functional()
                ),
            ));
        }
        let path = fs::canonicalize(path)?;
        species.push(json!({"name":name,"path":path,"sha256":sha256(&data),"z_valence":p.z_valence(),"functional":p.functional(),"n_projectors":p.projectors().len()}));
        inputs.push(path);
        potentials.push((name.clone(), p));
    }
    let xc_message = format!(
        "XC: {}; allow_functional_override={}; override_applied={}",
        opts.functional.name(),
        opts.allow_functional_override,
        potentials
            .iter()
            .any(|(_, p)| !opts.functional.matches_upf(p.functional()))
    );
    if mode.events.as_deref() == Some("-") {
        eprintln!("{xc_message}");
    } else {
        println!("{xc_message}");
    }
    let out = output_override.map(Path::to_owned).unwrap_or_else(|| {
        PathBuf::from(
            config
                .output
                .as_ref()
                .and_then(|o| o.directory.as_deref())
                .unwrap_or("results"),
        )
    });
    let out = if crate::parallel::rank() > 0 {
        out.join(format!("rank{}", crate::parallel::rank()))
    } else {
        out
    };
    fs::create_dir_all(&out)?;
    protect_outputs(&out, &inputs)?;
    volumetric_out::preflight(&config, &out).map_err(|e| error(5, e))?;
    let atoms: Vec<_> = config
        .system
        .atoms
        .iter()
        .map(|a| Atom::new(&a.species, a.fractional))
        .collect();
    let refs: Vec<_> = potentials
        .iter()
        .map(|(label, potential)| SpeciesPotential { label, potential })
        .collect();
    let points = config.resolved_points().map_err(|e| error(1, e))?;
    let electron_target: f64 = atoms
        .iter()
        .map(|a| {
            refs.iter()
                .find(|p| p.label == a.species)
                .unwrap()
                .potential
                .z_valence()
        })
        .sum();
    if let Some(theta) = config.finite_temperature().map_err(|e| error(1, e))? {
        return finite_run::execute(
            &config,
            &prepared,
            opts,
            &atoms,
            &refs,
            &points,
            &out,
            theta,
            electron_target,
            preflight,
            species,
            &bytes,
        );
    }
    if let Some(sector) = config
        .spin_sector(electron_target)
        .map_err(|e| error(1, e))?
    {
        return spin_run::execute(
            &config, &prepared, opts, &atoms, &refs, &points, sector, &out, preflight, species,
            &bytes,
        );
    }
    let mut problem = if config.basis.legacy_gamma() {
        ScfProblem::new(&prepared.cell, &atoms, &refs, opts.clone())
    } else {
        ScfProblem::with_kpoints(
            &prepared.cell,
            &atoms,
            &refs,
            opts.clone(),
            &points,
            dft_engine::OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-7,
            },
        )
    }
    .map_err(engine_error)?;
    if let Some(symmetry) = config.prepared_symmetry().map_err(|e| error(1, e))? {
        problem
            .enable_prepared_symmetry(&symmetry)
            .map_err(engine_error)?;
    }
    let density_floor = config
        .xc
        .density_floor
        .or_else(|| problem.core_density().map(|_| 1e-10));
    if let Some(floor) = density_floor {
        problem.set_density_floor(floor).map_err(engine_error)?;
    }
    let electrons = problem.electron_count();
    let setup = start.elapsed().as_secs_f64();
    let ladder = crate::config::Scf::cutoff_ladder(config.scf.as_ref())
        .map_err(|e| error(1, e.to_string()))?;
    let mut ladder_report = None;
    let outcome = if mode.observe || mode.agent {
        crate::telemetry::drive(problem, mode, &out)?
    } else {
        let run = match ladder {
            Some(fraction) => dft_engine::cutoff_ladder::run_ladder(
                problem,
                &refs,
                &dft_engine::cutoff_ladder::CutoffLadderOptions {
                    fraction,
                    ..Default::default()
                },
                true,
            )
            .map(|l| {
                ladder_report = Some(l.report);
                l.outcome
            }),
            None => problem.run_outcome(),
        };
        match run {
            Ok(outcome) => outcome,
            Err(e) => {
                if (density_floor.is_some() || opts.standard.is_some() || opts.hartree_mixing)
                    && let ScfError::Diverged { history, .. } | ScfError::NotConverged { history } =
                        &e
                {
                    let history: Vec<_> = history.iter().map(history_json).collect();
                    let mut failure = json!({
                        "convergence":{"converged":false,"iterations":history.len(),"history":history},
                        "failure":e.to_string(),"scf":options_json(&opts)});
                    if let Some(p) = opts.standard {
                        failure["convergence_profile"] = json!("standard");
                        failure["certified"] = json!(false);
                        failure["standard_threshold_ha"] = json!(p.threshold_ha);
                    }
                    write_json(&out.join("results.json"), &failure)?;
                    fs::write(
                        out.join("scf.log"),
                        format!(
                            "{}\n{}\n",
                            e,
                            serde_json::to_string_pretty(&history).map_err(|e| error(5, e))?
                        ),
                    )?;
                }
                return Err(engine_error(e));
            }
        }
    };
    let r = outcome.result;
    let ngrid = r.density.len();
    let bands = r.occupations.len();
    let npw = prepared.wavefunctions.vectors().len();
    let volume = prepared.cell.volume_bohr3();
    let cell = prepared.cell.lattice_columns_bohr();
    let density_bytes: Vec<u8> = r.density.iter().flat_map(|x| x.to_le_bytes()).collect();
    let integral = r.density.iter().sum::<f64>() * volume / ngrid as f64;
    let density_meta = json!({"schema_version":1,"fft_grid":prepared.grid.dims(),"lattice_columns_bohr":cell,"volume_bohr3":volume,
        "units":"electrons per Bohr^3","ordering":"C: (i*n2+j)*n3+k; k fastest","endianness":"little","dtype":"f64",
        "normalization":{"expected_electrons":electrons,"integral_electrons":integral},"sha256":sha256(&density_bytes)});
    let history: Vec<_> = r.history.iter().map(history_json).collect();
    let last = history.last().cloned();
    let mut phases = PhaseTimings::default();
    for h in &r.history {
        for t in std::iter::once(&h.timings)
            .chain(h.fixed_point_certificate.as_ref().map(|c| &c.timings))
        {
            phases.potential += t.potential;
            phases.diagonalization += t.diagonalization;
            phases.density += t.density;
            phases.energy += t.energy;
            phases.mixing += t.mixing;
        }
    }
    let total = start.elapsed().as_secs_f64();
    let energies = energy_json(&r.energy);
    let davidson_work: Vec<_> = dft_engine::work_accounting::take()
        .into_iter()
        .map(|w| {
            json!({
                "bands":w.bands,"strict_bands":w.strict_bands,"tolerance_ha":w.tolerance,
                "buffer_tolerance_ha":w.buffer_tolerance,"iterations":w.iterations,
                "block_calls":w.block_calls,"vector_applications":w.vector_applications,
                "buffer_origin_applications":w.buffer_origin_applications,
                "unattributed_applications":w.unattributed_applications,"restarts":w.restarts,
                "first_below":w.first_below,"last_above":w.last_above,"returned":w.returned,
                "locking":w.locking.map(|r| json!({"locked_roots":r.locked_roots,"unlocks":r.unlocks,"failures":r.failures}))
            })
        })
        .collect();
    let mut result = json!({"preflight":preflight,"schema_version":1,"tool":{"name":"dft-cli","version":env!("CARGO_PKG_VERSION"),"git_commit":env!("DFT_GIT_COMMIT"),"rustc_version":env!("DFT_RUSTC_VERSION")},
        "timestamp":timestamp(SystemTime::now().duration_since(UNIX_EPOCH).map_err(|e|error(5,e))?.as_secs()),"input_sha256":sha256(&bytes),"species":species,
        "cell":{"lattice_columns_bohr":cell,"volume_bohr3":volume,"reciprocal_columns_inverse_bohr":prepared.cell.reciprocal_columns_inverse_bohr()},
        "atoms":config.system.atoms.iter().map(|a|json!({"species":a.species,"fractional":a.fractional,"cartesian_bohr":prepared.cell.fractional_to_cartesian_bohr(a.fractional)})).collect::<Vec<_>>(),
        "cutoff":{"ecutwfc_ry":config.basis.ecutwfc_ry,"ecutwfc_ha":opts.cutoff_hartree},
        "basis":{"n_plane_waves":npw,"n_density_g":prepared.density.vectors().len(),"fft_grid":prepared.grid.dims(),"fft_grid_overridden":prepared.grid.overridden()},
        "electrons":electrons,"bands":bands,"occupations":r.occupations,"xc":opts.functional.name(),"scf":options_json(&opts),
        "convergence":{"converged":outcome.converged,"iterations":r.history.len(),"last_iteration":last,"history":history,
            "work":{
                "davidson_iterations_recorded":r.history.iter().map(|h|h.davidson_iterations).sum::<usize>(),
                "certificate_davidson_iterations":r.history.iter().filter_map(|h|h.fixed_point_certificate.as_ref()).map(|c|c.davidson_iterations).sum::<usize>(),
                "certificate_solves":r.history.iter().filter(|h|h.fixed_point_certificate.is_some()).count(),
                "second_vector_columns":dft_engine::second_vector_columns(),
                "davidson_stages":davidson_stages_json(),
                "b1":b1_json(),
                "b5":{"mode":dft_numerics::switches::choice(dft_numerics::switches::Id::B5),"bordered":dft_engine::b5_counters().bordered,"declined":dft_engine::b5_counters().declined,"full":dft_engine::b5_counters().full},
                "ortho":{"mode":dft_numerics::switches::text(dft_numerics::switches::Id::ORTHO).unwrap_or("dgks:0.1 (default)"),"second_pass_skipped":dft_engine::ortho_counters().0,"second_pass_repeated":dft_engine::ortho_counters().1},
                "algebra_stream":algebra_stream_json(),
                "h_vectors_recorded":r.history.iter().map(|h|h.operator_applications).sum::<usize>(),
                "certificate_h_vectors":r.history.iter().filter_map(|h|h.fixed_point_certificate.as_ref()).map(|c|c.operator_applications).sum::<usize>(),
                "note":"recorded = per-iteration solves in history; certificate = fixed-point certificate re-solves, reported only inside each history entry's fixed_point_certificate; total work is the sum"}},"energies":energies,
        "total_energy":{"ha":r.energy.total_hartree(),"ry":from_hartree(r.energy.total_hartree(), EnergyUnit::Rydberg),"ev":from_hartree(r.energy.total_hartree(), EnergyUnit::ElectronVolt)},
        "eigenvalues":{"ha":r.eigenvalues,"ev":r.eigenvalues.iter().map(|e|from_hartree(*e, EnergyUnit::ElectronVolt)).collect::<Vec<_>>()},
        "timings":{"phases":phase_json(&phases),"setup_s":setup,"engine_setup_s":r.setup_time.as_secs_f64(),"engine_total_s":r.total_time.as_secs_f64(),"total_s":total,"fft_hpsi":"included in enclosing phases; not separately measured"},
        "threading":{"threads":opts.threads.max(1),"threads_requested":opts.threads,"threads_used":opts.threads.max(1),"outer_threads":opts.threads.max(1).min(r.k_points.len()),"inner_threads":opts.threads.max(1),"backend":"persistent Rayon pool shared by k points and fixed inner pieces; sequential faer and serial FFT per worker"},
        "memory_estimate":{"basis_bytes":npw*std::mem::size_of_val(&prepared.wavefunctions.vectors()[0]),"density_basis_bytes":prepared.density.vectors().len()*std::mem::size_of_val(&prepared.density.vectors()[0]),"density_grid_bytes":ngrid*8,"one_complex_grid_bytes":ngrid*16,"orbital_bytes":npw*bands*16,"peak_rss_bytes":null,"note":"array payload estimates only; scratch, projectors, allocator overhead excluded; peak RSS unavailable in-process"},
        "units":"Ha, Ry, eV, Bohr, seconds and bytes as keyed; density_residual in Bohr^-3/2; electrons/occupations in electron counts; relative errors and fractions dimensionless; see docs/cli.md"});
    if let Some(l) = &ladder_report {
        result["convergence"]["work"]["cutoff_ladder"] = json!({
            "used":l.used,"fallback_reason":l.fallback_reason,"fraction":l.fraction,
            "coarse_cutoff_ha":l.coarse_cutoff_hartree,"coarse_fft_grid":l.coarse_grid,
            "fine_fft_grid":l.fine_grid,"coarse_plane_waves":l.coarse_plane_waves,
            "fine_plane_waves":l.fine_plane_waves,"coarse_iterations":l.coarse_iterations,
            "coarse_final_density_residual":l.coarse_final_density_residual,
            "coarse_davidson_iterations":l.coarse_davidson_iterations,
            "coarse_h_vectors":l.coarse_h_vectors,"coarse_start_h_vectors":l.coarse_start_h_vectors,
            "coarse_h_full_grid_equivalents":l.coarse_h_equivalents(),
            "coarse_setup_s":l.coarse_setup.as_secs_f64(),"coarse_stage_s":l.coarse_time.as_secs_f64(),
            "note":"coarse stage supplies the start only; convergence, energy audit and certificate are at the full cutoff (history and work above are full-cutoff only)"});
    }
    if let Some(c) = &r.energy_identity {
        // B2 (BAANDI_ENERGY_IDENTITY=on): history energies of non-final iterations
        // use the band-sum identity; the final energy is the explicit audit.
        result["convergence"]["work"]["energy_identity"] = json!({
            "history_energy_form":"band_sum_identity (non-final iterations)",
            "final_energy_form":"explicit_audit",
            "identity_nonlocal_ha":c.identity_nonlocal_ha,
            "audit_nonlocal_ha":c.audit_nonlocal_ha,
            "difference_ha":c.difference_ha});
    }
    if let Some(p) = opts.standard {
        if opts.standard_stall_guard {
            result["standard_stall_guard"] = json!("on");
        }
        if opts.standard_density_forcing {
            result["standard_density_forcing"] = json!("on");
        }
        result["convergence_profile"] = json!("standard");
        result["certified"] = json!(false);
        result["standard_threshold_ha"] = json!(p.threshold_ha);
        result["work"] = json!({"davidson_steps":r.history.iter().map(|h|h.davidson_iterations).sum::<usize>(),
            "h_vectors":r.history.iter().map(|h|h.operator_applications).sum::<usize>() + r.spectral_completeness.as_ref().map_or(0, |s| s.operator_applications),
            "first_iteration_redos":r.history.iter().filter(|h|h.standard_redo).count(),
            "certificate_solves":0});
        println!(
            "convergence_profile: standard; certified: false; standard_threshold_ha: {}",
            p.threshold_ha
        );
    }
    if dft_engine::work_accounting::ENABLED {
        result["davidson_work"] = json!(davidson_work);
    }
    if opts.start_wavefunctions == dft_engine::StartWavefunctions::Atomic {
        let a = &r.atomic_start;
        result["atomic_start"] = json!({
            "requested": opts.start_wavefunctions.name(), "fallback_reason": a.fallback_reason,
            "atomic_orbitals_used": a.atomic_orbitals,
            "initial_h_applications": a.initial_h_applications,
            "gram_rotation_complex_macs": a.gram_rotation_macs,
            "mgs_projection_complex_macs": a.mgs_projection_macs,
            "radial_samples": a.radial_samples, "evd_orders": a.evd_orders
        });
    }
    result["gamma_real_used"] = json!(r.gamma_real_used);
    if crate::parallel::size() > 1 {
        let m2_on = dft_engine::parallel::m2_enabled();
        let (m2_counts, m2_stages) = dft_engine::parallel::m2_stats();
        result["parallel"] = json!({
            "model": if m2_on { "mpi_stick_rows" } else { "mpi_row_block" },
            "m2": {
                "enabled": m2_on,
                "batched_transposes": dft_numerics::switches::on(dft_numerics::switches::Id::M2_BATCH),
                "counters": m2_counts.into_iter().map(|(n, v)| (n, json!(v))).collect::<serde_json::Map<_, _>>(),
                "stages": m2_stages
                    .into_iter()
                    .map(|(name, seconds)| json!({"stage": name, "seconds": seconds}))
                    .collect::<Vec<_>>(),
            },
            "ranks": crate::parallel::size(),
            "rank": crate::parallel::rank(),
            "distributed_solves": dft_engine::parallel::distributed_solves(),
            "stages": dft_engine::parallel::stage_seconds()
                .into_iter()
                .map(|(name, seconds)| json!({"stage": name, "seconds": seconds}))
                .collect::<Vec<_>>(),
            "collectives": dft_engine::parallel::collective_stats()
                .into_iter()
                .map(|(op, calls, seconds)| json!({"op": op, "calls": calls, "seconds": seconds}))
                .collect::<Vec<_>>(),
        });
    }
    result["gamma_real_fallback_reason"] = json!(r.gamma_real_fallback_reason);
    if let Some(report) = &r.spectral_completeness {
        result["spectral_completeness"] = spectral_json(report);
    }

    if let Some(symmetry) = config.symmetry_json().map_err(|e| error(1, e))? {
        result["symmetry"] = symmetry;
    }
    if opts.locking == dft_engine::Locking::Auto {
        result["locking"] = json!(r.k_points.iter().map(|k| k.solution.locking.as_ref().map(|r| json!({
            "locked_roots":r.locked_roots, "unlocked_roots":r.unlocks,
            "certificate_failures":r.failures.iter().map(|&(iteration,root,residual)| json!({"iteration":iteration,"root":root,"residual_ha":residual})).collect::<Vec<_>>()
        }))).collect::<Vec<_>>());
    }
    result["band_accuracy"] = json!(r.k_points[0].solution.residuals.iter()
        .zip(&r.k_points[0].solution.buffer).zip(&r.k_points[0].solution.converged)
        .map(|((&residual, &buffer), &converged)| json!({"residual_ha":residual,"buffer":buffer,"converged":converged}))
        .collect::<Vec<_>>());
    // Preserve the legacy LDA schema unless an override was explicitly enabled.
    if let Some(metadata) = reference_metadata(&opts) {
        result["compat_mode"] = json!(r.compat_mode.name());
        result["compat"] = metadata;
    }
    if opts.local_mesh != dft_engine::LocalMesh::Baandi
        || opts.radial_even_rule != dft_engine::RadialEvenRule::Simpson3_8
    {
        result["pseudo"] =
            json!({"local_mesh":r.local_mesh.name(),"radial_even_rule":r.radial_even_rule.name()});
    }
    if opts.pw92_prefactor == dft_engine::Pw92Prefactor::Original {
        result["pw92_prefactor"] = json!(r.pw92_prefactor.name());
    }
    if opts.functional != dft_engine::Functional::LdaPz || opts.allow_functional_override {
        result["allow_functional_override"] = json!(opts.allow_functional_override);
        result["functional_override_applied"] = json!(
            potentials
                .iter()
                .any(|(_, p)| !opts.functional.matches_upf(p.functional()))
        );
    }
    if !config.input_expressions.is_empty() {
        result["input_expressions"] = json!(config.input_expressions);
    }
    if !config.basis.legacy_gamma() {
        result["kpoints"] = json!({
            "grid": config.basis.kgrid, "shift": config.basis.kgrid.map(|_| config.basis.kshift.unwrap_or([0.0;3])),
            "explicit": config.basis.kgrid.is_none(), "count": r.k_points.len(),
            "time_reversal": config.basis.time_reversal(),
            "points": r.k_points.iter().map(|k| json!({"fractional": k.point.fractional,
                "weight": k.point.weight, "eigenvalues": {"ha": k.solution.eigenvalues,
                "ev": k.solution.eigenvalues.iter().map(|e| from_hartree(*e, EnergyUnit::ElectronVolt)).collect::<Vec<_>>()},
                "occupations": r.occupations,
                "band_accuracy": k.solution.residuals.iter().zip(&k.solution.buffer).zip(&k.solution.converged)
                    .map(|((&residual, &buffer), &converged)| json!({"residual_ha":residual,"buffer":buffer,"converged":converged})).collect::<Vec<_>>()
                })).collect::<Vec<_>>()
        });
    }
    let mut console = grid_report(&config, &prepared);
    writeln!(
        console,
        "Spectral completeness: {}",
        spectral_status(r.spectral_completeness.as_ref())
    )
    .unwrap();
    if let Some(floor) = density_floor {
        result["initial_density"] = json!({"method":"uniform", "expected_electrons":electrons,
            "atomic_table_used":false, "density_floor_electrons_bohr3":floor});
        for h in &r.history {
            if let Some(f) = &h.energy.xc_floor {
                writeln!(
                    console,
                    "XC floor iteration {}: clipped negative charge input={} output={} electrons",
                    h.iteration, f.input_negative_charge, f.output_negative_charge
                )
                .unwrap();
            }
        }
    }
    writeln!(
        console,
        "Electrons: {electrons}; bands: {bands}; mixing: {:?}",
        opts.mixing
    )
    .unwrap();
    writeln!(console,"iter E_total(Ha) E_total(Ry) delta_E(Ha) density_residual(Bohr^-3/2) max_eigen_residual(Ha) Davidson phase_wall(s)").unwrap();
    for h in &r.history {
        writeln!(
            console,
            "{} {:.15} {:.15} {} {:.6e} {:.6e} {} {:.6}",
            h.iteration,
            h.energy.total_hartree(),
            h.energy.total_rydberg(),
            h.energy_change
                .map(|x| format!("{x:.6e}"))
                .unwrap_or_else(|| "n/a".into()),
            h.density_residual,
            h.max_eigen_residual,
            h.davidson_iterations,
            phase_sum(&h.timings)
                + h.fixed_point_certificate
                    .as_ref()
                    .map(|c| phase_sum(&c.timings))
                    .unwrap_or(0.)
        )
        .unwrap();
    }
    for (term, v) in energies.as_object().unwrap() {
        if v.get("ha").is_none() {
            continue;
        }
        writeln!(
            console,
            "{term}: {:.15} Ha = {:.15} Ry",
            v["ha"].as_f64().unwrap(),
            v["ry"].as_f64().unwrap()
        )
        .unwrap();
    }
    let count = r.k_points.len();
    for (i, k) in r.k_points.iter().enumerate() {
        if count > 10 && i >= 3 && i < count - 3 {
            if i == 3 {
                writeln!(console, "... remaining k points: see results.json").unwrap();
            }
            continue;
        }
        writeln!(
            console,
            "k{} {:?}; weight {:.12}",
            i + 1,
            k.point.fractional,
            k.point.weight
        )
        .unwrap();
        for (b, (e, f)) in k
            .solution
            .eigenvalues
            .iter()
            .zip(&r.occupations)
            .enumerate()
        {
            writeln!(
                console,
                "band {}: {:.12} Ha = {:.12} eV; occupation {f}{}",
                b + 1,
                e,
                from_hartree(*e, EnergyUnit::ElectronVolt),
                if k.solution.buffer[b] {
                    " (buffer; less accurate)"
                } else {
                    ""
                }
            )
            .unwrap();
        }
    }
    let occupied = r.occupations.iter().filter(|f| **f > 0.).count();
    if bands > occupied {
        let mut homo = (0, f64::NEG_INFINITY);
        let mut lumo = (0, f64::INFINITY);
        let mut direct = (0, f64::INFINITY);
        for (i, k) in r.k_points.iter().enumerate() {
            let h = k.solution.eigenvalues[occupied - 1];
            let l = k.solution.eigenvalues[occupied];
            if h > homo.1 {
                homo = (i + 1, h);
            }
            if l < lumo.1 {
                lumo = (i + 1, l);
            }
            if l - h < direct.1 {
                direct = (i + 1, l - h);
            }
        }
        if density_floor.is_some() {
            result["band_gap"] = json!({"direct_ha":direct.1,"indirect_ha":lumo.1-homo.1});
        }
        writeln!(console, "HOMO/LUMO gap: direct {:.12} Ha = {:.12} eV at k{}; indirect {:.12} Ha = {:.12} eV (HOMO k{}, LUMO k{})",
            direct.1, from_hartree(direct.1, EnergyUnit::ElectronVolt), direct.0,
            lumo.1-homo.1, from_hartree(lumo.1-homo.1, EnergyUnit::ElectronVolt), homo.0, lumo.0).unwrap();
        if !outcome.converged && lumo.1 - homo.1 < 1e-3 {
            // Evidence: ledger ELONGATED-GAP / ELONGATED-SMEARING / SLAB-SMEARING.
            writeln!(console, "hint: the HOMO/LUMO gap is {:.3e} Ha (below 1e-3 Ha) at fixed occupations; an integer occupation cut inside a (near-)degenerate level cluster can stall the SCF. Use [occupations] mode = \"fermi_dirac\" with a smearing width above the cluster splitting (strict profile, band_cap above the occupied count).", lumo.1 - homo.1).unwrap();
        }
    }
    writeln!(
        console,
        "{}",
        if outcome.converged {
            "converged"
        } else {
            "NOT CONVERGED: iteration budget exhausted"
        }
    )
    .unwrap();
    writeln!(
        console,
        "Total wall: {total:.6} s; setup: {setup:.6} s; phases: {}",
        phase_json(&phases)
    )
    .unwrap();
    writeln!(
        console,
        "FFT/H psi included in enclosing phases; not separately measured."
    )
    .unwrap();
    volumetric_out::write_outputs(&config, &prepared, &refs, &r, &out).map_err(|e| error(5, e))?;
    fs::write(out.join("density.bin"), density_bytes)?;
    write_json(&out.join("density.json"), &density_meta)?;
    if let Some(force) = requested_forces(
        &config,
        &prepared.cell,
        &atoms,
        &refs,
        &opts,
        dft_engine::forces::ElectronicState::Fixed(&r),
        density_floor,
    )? {
        console.push_str(&force_table(&force));
        result["forces"] = force_json(&force);
    }
    if let Some(stress) = requested_stress(
        &config,
        &prepared.cell,
        &atoms,
        &refs,
        &opts,
        dft_engine::forces::ElectronicState::Fixed(&r),
        density_floor,
    )? {
        console.push_str(&stress_table(&stress));
        result["stress"] = stress_json(&stress);
    }
    write_json(&out.join("results.json"), &result)?;
    Ok(RunReport {
        console,
        converged: outcome.converged,
    })
}

/// Reporting boundary only; untested states never receive a passing flag.
pub(crate) fn spectral_status(report: Option<&dft_engine::SpectralReport>) -> &'static str {
    match report {
        None => "not verified (scf.spectral_guard=off)",
        Some(r) if r.passed => {
            "passed independent spectral check (see method and scope in results.json)"
        }
        Some(_) => "could not verify the reported state",
    }
}
pub(crate) fn spectral_json(report: &dft_engine::SpectralReport) -> Value {
    json!({"method":"final_reported_hamiltonian", "passed":report.passed,
        "retries":report.retries, "extra_operator_applications":report.operator_applications,
        "timings":{"wall_s":report.elapsed.as_secs_f64()},
        "verification_scope":"full-space dense for dimension <=128; independent seeded search otherwise, not an inertia proof",
        "status":if report.passed {"verified_by_method"} else {"could_not_verify_reported_state"},
        "events":report.events.iter().map(|e| json!({"iteration":e.iteration,"spin":e.spin,
            "k_point":e.k_point,"method":e.method,"passed":e.passed,"retries":e.retries,
            "extra_operator_applications":e.operator_applications,
            "maximum_discrepancy_ha":e.discrepancy_ha,"error":e.error})).collect::<Vec<_>>()})
}

#[allow(clippy::too_many_arguments)]
fn requested_forces(
    config: &Config,
    cell: &dft_geometry::Cell,
    atoms: &[Atom],
    potentials: &[SpeciesPotential<'_>],
    options: &ScfOptions,
    electronic: dft_engine::forces::ElectronicState<'_>,
    density_floor: Option<f64>,
) -> Result<Option<dft_engine::forces::ForceResult>, RunError> {
    if !config.output.as_ref().is_some_and(|o| o.forces) {
        return Ok(None);
    }
    let operations: Vec<_> =
        config
            .symmetry_map()
            .map_err(|e| error(1, e))?
            .map_or_else(Vec::new, |m| {
                m.operations
                    .into_iter()
                    .filter(|o| !o.antiunitary)
                    .collect()
            });
    dft_engine::forces::forces(&dft_engine::forces::ForceState {
        cell,
        atoms,
        potentials,
        options,
        electronic,
        density_floor,
        symmetry: &operations,
    })
    .map(Some)
    .map_err(engine_error)
}
fn force_json(f: &dft_engine::forces::ForceResult) -> Value {
    json!({"status":"available", "units":f.units,"frame":"cartesian_input_atom_order",
        "quantity":f.quantity,"values":f.values,
        "preconditions":{"convergence_profile":"strict","converged":true,"fixed_point_certificate_passed":true},
        "derivative_ledger":f.ledger.iter().map(|e|json!({"term_id":e.term_id,"owner":e.owner,"active":e.active,"energy_ha":e.energy_ha,"included_in":e.included_in,"dependencies":e.dependencies,"force_capability":if e.active {"available"} else {"not_applicable"},"stress_capability":"unavailable","reason":e.reason})).collect::<Vec<_>>(),"per_term":{"ewald":f.ewald,"local":f.local,"nonlocal":f.nonlocal,"core":f.core},
        "net_force":f.net_force,"maximum_absolute_component":f.maximum_absolute_component,
        "elapsed_seconds":f.elapsed.as_secs_f64(),"band_passes":f.band_passes,"unavailable_terms":[]})
}
fn force_table(f: &dft_engine::forces::ForceResult) -> String {
    let mut text = format!(
        "\nCartesian forces (-d{}/dR), Ha/Bohr\n atom               Fx               Fy               Fz\n",
        f.quantity
    );
    for (i, v) in f.values.iter().enumerate() {
        writeln!(
            text,
            "{:5} {:16.9e} {:16.9e} {:16.9e}",
            i + 1,
            v[0],
            v[1],
            v[2]
        )
        .unwrap();
    }
    text
}

#[allow(clippy::too_many_arguments)]
fn requested_stress(
    config: &Config,
    cell: &dft_geometry::Cell,
    atoms: &[Atom],
    potentials: &[SpeciesPotential<'_>],
    options: &ScfOptions,
    electronic: dft_engine::forces::ElectronicState<'_>,
    density_floor: Option<f64>,
) -> Result<Option<dft_engine::stress::StressResult>, RunError> {
    if !config.output.as_ref().is_some_and(|o| o.stress) {
        return Ok(None);
    }
    let operations: Vec<_> =
        config
            .symmetry_map()
            .map_err(|e| error(1, e))?
            .map_or_else(Vec::new, |m| {
                m.operations
                    .into_iter()
                    .filter(|o| !o.antiunitary)
                    .collect()
            });
    dft_engine::stress::stress(&dft_engine::forces::ForceState {
        cell,
        atoms,
        potentials,
        options,
        electronic,
        density_floor,
        symmetry: &operations,
    })
    .map(Some)
    .map_err(engine_error)
}

fn stress_json(s: &dft_engine::stress::StressResult) -> Value {
    let pulay=s.pulay.as_ref().map(|p|json!({"method":p.method,"isotropic_correction_ha_bohr3":p.correction,"uncertainty_ha_bohr3":p.uncertainty,"cutoff_steps_ha":p.cutoff_steps_ha,"warning_threshold_ha_bohr3":p.warning_threshold,"warning":p.warning,"shells_resolved":p.shells_resolved,"uncertainty_kind":"two_step_spread_not_bound","anisotropic_convergence":"not_established","elapsed_seconds":p.elapsed.as_secs_f64()}));
    json!({"status":"available","units":s.units,"frame":"cartesian_tensor","sign_convention":"tensile_positive","quantity":s.quantity,"values":s.values,"pressure_ha_bohr3":s.pressure,"symmetric":s.symmetric,"antisymmetric":s.antisymmetric,"antisymmetric_norm_ha_bohr3":s.antisymmetric_norm,"basis_protocol":s.basis_protocol,"cutoff_ha":s.cutoff_ha,"strain_steps":s.strain_steps,"cutoff_converged":false,"pulay_assessment":pulay,"per_term":s.terms.iter().map(|(name,t)|((*name).to_string(),json!(t))).collect::<serde_json::Map<_,_>>(),"derivative_ledger":s.ledger.iter().map(|e|json!({"term_id":e.term_id,"owner":e.owner,"energy_ha":e.energy_ha,"included_in":e.included_in,"dependencies":e.dependencies,"active":e.active,"stress_capability":if e.active {"available"}else{"not_applicable"},"reason":e.reason})).collect::<Vec<_>>(),"unavailable_terms":[],"preconditions":{"convergence_profile":"strict","converged":true,"fixed_point_certificate_passed":true},"elapsed_seconds":s.elapsed.as_secs_f64(),"band_passes":s.band_passes,"precondition_force_band_passes":s.precondition_force_band_passes,"cutoff_scf_runs":3})
}
fn stress_table(s: &dft_engine::stress::StressResult) -> String {
    let mut text = format!(
        "\nFrozen-Miller tensile-positive stress, Ha/Bohr^3: {:?}\nPressure (-trace/3): {:.12e} Ha/Bohr^3\n",
        s.values, s.pressure
    );
    if let Some(p) = &s.pulay {
        text.push_str(&format!("Isotropic Pulay estimate: {:.6e} +/- {:.6e} Ha/Bohr^3; cutoff convergence not established\n",p.correction,p.uncertainty));
        if p.warning {
            text.push_str("WARNING: Pulay estimate exceeds threshold or sampled shells are unresolved; increase cutoff and repeat.\n");
        }
    }
    text
}

/// Real-Davidson stage timers (`BAANDI_DAVIDSON_STAGES=on`; zeros otherwise),
/// cumulative over the process, seconds.
fn davidson_stages_json() -> serde_json::Value {
    let s = dft_engine::stage_seconds();
    let labels = dft_engine::DAVIDSON_WIDTH_LABELS;
    json!({
        "timed_solves":s.timed_solves,
        "iteration_s_by_correction_width":labels.iter().zip(s.width_seconds).map(|(l,v)|(l.to_string(),json!(v))).collect::<serde_json::Map<_,_>>(),
        "iterations_by_correction_width":labels.iter().zip(s.width_iterations).map(|(l,v)|(l.to_string(),json!(v))).collect::<serde_json::Map<_,_>>(),
        "stage_s":{"reduced":s.reduced,"eigen":s.eigen,"rotate":s.rotate,"residual":s.residual,"orthonormalize":s.orthonormalize,"h_apply":s.h_apply,"sketch":s.sketch,"precondition":s.precondition},
        "note":"BAANDI_DAVIDSON_STAGES=on; serial real Davidson only; width 0 = returning iteration; h_apply includes certificate and restart applications"
    })
}

/// B1 selected-root counters (`BAANDI_B1`; zeros when off), cumulative.
fn b1_json() -> serde_json::Value {
    let c = dft_engine::b1_counters();
    json!({
        "mode":dft_numerics::switches::choice(dft_numerics::switches::Id::B1),
        "sketch":dft_numerics::switches::choice(dft_numerics::switches::Id::B1_SKETCH),
        "solves":c.solves,"narrow_iterations":c.narrow_iterations,"full_iterations":c.full_iterations,
        "narrow_columns":c.narrow_columns,"avoided_columns":c.avoided_columns,"fallbacks":c.fallbacks,
        "sketch_estimates":c.sketch_estimates,"sketch_disagreements":c.sketch_disagreements
    })
}

/// Algebra-stream switches of the real Davidson and how often each acted
/// (docs/performance/requests/algebra-stream.md).
fn algebra_stream_json() -> serde_json::Value {
    let c = dft_engine::algebra_counters();
    use dft_numerics::switches::Id;
    let mode = |id: Id| dft_numerics::switches::text(id).unwrap_or("off");
    json!({
        "ritz_fused":mode(Id::RITZ_FUSED),"narrow_proj":mode(Id::NARROW_PROJ),
        "narrow_reduced":mode(Id::NARROW_REDUCED),
        "fused_iterations":c.fused_iterations,"deferred_hx":c.deferred_hx,
        "narrow_projections":c.narrow_projections,"narrow_reduced_borders":c.narrow_reduced,
        "b5_parallel_updates":c.b5_parallel
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn utc_dates() {
        for (s, t) in [
            (0, "1970-01-01T00:00:00Z"),
            (951782400, "2000-02-29T00:00:00Z"),
            (1709251199, "2024-02-29T23:59:59Z"),
            (4107542400, "2100-03-01T00:00:00Z"),
        ] {
            assert_eq!(timestamp(s), t);
        }
    }
    #[test]
    fn options_validation() {
        let example = include_str!("../../../examples/si_lda.toml");
        for fields in [
            "max_iterations=0",
            "mixing='bad'",
            "mixing_beta=nan",
            "energy_tolerance_ha=-1.0",
            "density_residual_tolerance=0.0",
            "eigen_tolerance=inf",
        ] {
            let c = Config::parse(&format!("{example}\n[scf]\n{fields}")).unwrap();
            assert_eq!(options(&c).unwrap_err().code, 1);
        }
        let c = Config::parse(example).unwrap();
        assert_eq!(options(&c).unwrap().seed, 2026);
    }
}
