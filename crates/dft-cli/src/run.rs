#[path = "finite_run.rs"]
mod finite_run;
// File/console boundary for the sequential engine. No I/O occurs in SCF kernels.
use crate::{config::Config, grid_report, run_config::launch};
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

// Bootstrap stamps are also retained when detailed profiling is disabled.
// Library callers have no main() boundary and report that interval as null.
thread_local! {
    static REFERENCE_LAUNCH: std::cell::Cell<Option<[Instant; 4]>> = const { std::cell::Cell::new(None) };
}
pub fn reference_launch_begin(start: Instant, mpi: Instant, args: Instant, config: Instant) {
    REFERENCE_LAUNCH.with(|v| v.set(Some([start, mpi, args, config])));
}
type PhaseSeconds = std::collections::BTreeMap<&'static str, f64>;
fn phase<T>(
    enabled: bool,
    seconds: &mut PhaseSeconds,
    name: &'static str,
    f: impl FnOnce() -> T,
) -> T {
    if !enabled {
        return f();
    }
    let clock = Instant::now();
    let out = f();
    *seconds.entry(name).or_default() += clock.elapsed().as_secs_f64();
    out
}
fn start_report_json(report: dft_engine::refstart::StartReport) -> Value {
    use dft_engine::refstart::{DensityStartKind, OrbitalStartKind};
    json!({
        "density": match report.density {
            DensityStartKind::Atomic => "atomic",
            DensityStartKind::Uniform => "uniform",
            DensityStartKind::Supplied => "supplied",
        },
        "orbitals": match report.orbitals {
            OrbitalStartKind::Atomic => "atomic",
            OrbitalStartKind::AtomicWithKineticCompletion => "atomic_with_kinetic_completion",
            OrbitalStartKind::Kinetic => "kinetic",
            OrbitalStartKind::Supplied => "supplied",
        }
    })
}

fn boundary_json(t: &dft_engine::refscf::record::SetupTimings) -> Value {
    json!({"phases":t.phases,"collectives_nested":t.collective_s,"total":t.total_s})
}

/// Write the large body once, then append the small timing trailer. The clocks
/// truthfully exclude that trailer and close; no second full serialization or
/// rewrite is performed to make a results file describe its own write.
fn write_reference_results(
    path: &Path,
    mut result: Value,
    mut output: Value,
) -> Result<(), RunError> {
    use std::io::Write;
    let clock = Instant::now();
    let reference = result
        .as_object_mut()
        .expect("results object")
        .remove("reference_scf")
        .expect("reference object");
    let body = serde_json::to_vec_pretty(&result).map_err(|e| error(5, e))?;
    let reference = serde_json::to_vec_pretty(&reference).map_err(|e| error(5, e))?;
    output["json_serialize_body"] = json!(clock.elapsed().as_secs_f64());
    let clock = Instant::now();
    let mut file = fs::File::create(path)?;
    file.write_all(&body[..body.len() - 1])?;
    file.write_all(b",\n\"reference_scf\":")?;
    file.write_all(&reference[..reference.len() - 1])?;
    output["json_write_body"] = json!(clock.elapsed().as_secs_f64());
    output["json_scope"] = json!(
        "body serialization and write including file create; excludes small output_s trailer, close, launch profile patch and MPI finalization"
    );
    file.write_all(b",\n\"output_s\":")?;
    let trailer = serde_json::to_vec(&output).map_err(|e| error(5, e))?;
    file.write_all(&trailer)?;
    // The existing launch-profile appender requires a final brace, no newline.
    file.write_all(b"}}")?;
    Ok(())
}

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
/// The SCF options of every command (`run` in every execution mode, `relax`,
/// `bands`, `dos`, `checkpoint`, `restart`, `explain`): the single standard
/// path with the verification, or the strict reference path under
/// `BAANDI_REFERENCE_STRICT=on` (ADR 0043, single path complete since
/// ONEPATH-4). The path fields are set explicitly from `Config::scf_path`.
pub fn options(config: &Config) -> Result<ScfOptions, RunError> {
    options_for(config, crate::driver_select::preliminary(config))
}
pub fn options_for(
    config: &Config,
    selection: crate::driver_select::Selection,
) -> Result<ScfOptions, RunError> {
    let path = config
        .scf_path_selected(selection)
        .map_err(|e| error(1, e))?;
    if selection.reference() {
        reference_options(config)
    } else {
        let mut options = oracle_options(config, path)?;
        // G4D accepts the reference mixer aliases/history but now uses the
        // strict oracle's beta rule.
        if selection.reference_mixing()
            && let Some(s) = &config.scf
            && matches!(s.mixing.as_deref(), None | Some("pulay"))
        {
            let profile = dft_engine::refscf::profile::Profile {
                beta: s.mixing_beta.unwrap_or(0.7),
                history: s.mixing_history.unwrap_or(8),
                ..Default::default()
            }
            .validate()
            .map_err(|e| error(1, e))?;
            // oracle_options already applies the class-dependent beta rule,
            // which also ran after the old reference mixing adapter.
            if let dft_engine::Mixing::Pulay { beta, .. } = options.mixing {
                options.mixing = dft_engine::Mixing::Pulay {
                    beta,
                    history: profile.history,
                };
            }
        }
        Ok(options)
    }
}
pub fn reference_requested(config: &Config) -> bool {
    crate::driver_select::preliminary(config).reference()
}
#[cfg(test)]
use crate::driver_select::G4D_FALLBACK;
#[cfg(test)]
fn reference_spin_fd_fallback(config: &Config) -> bool {
    reference_spin_fd_fallback_reason(config).is_some()
}
#[cfg(test)]
fn reference_spin_fd_fallback_reason(config: &Config) -> Option<&'static str> {
    crate::driver_select::preliminary(config)
        .reason
        .filter(|r| *r == G4D_FALLBACK)
}

/// Independent profile resolution: no legacy session tuning or ladder reads.
pub fn reference_profile(
    config: &Config,
) -> Result<dft_engine::refscf::profile::Profile, RunError> {
    use dft_engine::refscf::profile::{Preconditioning, Profile};
    let s = config.scf.as_ref();
    if s.and_then(|s| s.mixing.as_deref())
        .is_some_and(|m| !matches!(m, "plain" | "broyden" | "pulay"))
    {
        return Err(error(
            1,
            "reference SCF mixing must be plain/broyden (pulay is an alias for the difference-history method)",
        ));
    }
    Profile {
        beta: s.and_then(|s| s.mixing_beta).unwrap_or(0.7),
        history: s.and_then(|s| s.mixing_history).unwrap_or(8),
        threshold_ry: 2.0
            * s.and_then(|s| s.standard_threshold_ha).unwrap_or_else(|| {
                if config.wants_stress() {
                    dft_engine::StandardProfile::STRESS_THRESHOLD_HA
                } else {
                    dft_engine::StandardProfile::default().threshold_ha
                }
            }),
        preconditioning: if s.and_then(|s| s.reference_preconditioning.as_deref()) == Some("tf") {
            Preconditioning::ThomasFermi
        } else {
            Preconditioning::Plain
        },
    }
    .validate()
    .map_err(|e| error(1, e))
}

fn reference_symmetry(
    config: &Config,
) -> Result<Option<dft_engine::refscf::symmetry::Prepared>, RunError> {
    config
        .prepared_symmetry()
        .map_err(|e| error(1, e))?
        .map(|p| {
            let atoms: Vec<_> = config
                .system
                .atoms
                .iter()
                .map(|a| Atom::new(&a.species, a.fractional))
                .collect();
            dft_engine::refscf::symmetry::Prepared::from_prepared(
                &config.cell().map_err(|e| error(1, e))?,
                &atoms,
                &config.symmetry_full_mesh().map_err(|e| error(1, e))?,
                config.symmetry.tolerance_bohr,
                config.basis.time_reversal(),
                &p,
            )
            .map_err(engine_error)
        })
        .transpose()
}
fn reference_symmetry_json(report: &dft_engine::refscf::symmetry::Report) -> Value {
    json!({"group_order":report.group_order,"irreducible_k_count":report.irreducible_k_count,
        "channels":report.channels,"residual_units":"electrons/Bohr^(3/2)",
        "residual_definition":"pre-projection ||rho - P rho||; residual_norm aliases pre_projection_norm",
        "pre_projection_relative_units":"dimensionless",
        "pre_projection_relative_definition":"||rho - P rho|| / ||rho|| (zero for zero density); rho is the irreducible-k sum before projection",
        "post_projection_invariance_definition":"||P(P rho) - P rho||; diagnostic only, expected at rounding",
        "maps":report.maps.iter().map(|r| json!({"iteration":r.iteration,"verification":r.verification,
            "residual_norm":r.norm,"pre_projection_norm":r.norm,
            "pre_projection_relative_norm":r.pre_projection_relative_norm,
            "post_projection_invariance_norm":r.post_projection_invariance_norm})).collect::<Vec<_>>()})
}

fn reference_options(config: &Config) -> Result<ScfOptions, RunError> {
    let profile = reference_profile(config)?;
    let s = config.scf.as_ref();
    let finite = config.finite_temperature().map_err(|e| error(1, e))?;
    if config.collinear() && config.xc.functional.engine() == dft_engine::Functional::LdaPz {
        return Err(error(
            1,
            "reference spin XC requires lda-pw or pbe: polarized lda-pz is unavailable in legacy",
        ));
    }
    let points = config.resolved_points().map_err(|e| error(1, e))?;
    let mut o = ScfOptions {
        compat_mode: if config.compat.reference() {
            dft_engine::CompatMode::Reference
        } else {
            dft_engine::CompatMode::None
        },
        crystal_density_group: config.symmetry.density_group
            == crate::config::DensityGroup::Crystal,
        start_density: dft_engine::StartDensity::Atomic,
        start_wavefunctions: dft_engine::StartWavefunctions::Atomic,
        gamma_real: if !config.collinear()
            && finite.is_none()
            && points.len() == 1
            && points[0].fractional == [0.; 3]
            && s.and_then(|s| s.gamma_real.as_deref()) != Some("off")
        {
            dft_engine::GammaReal::On
        } else {
            dft_engine::GammaReal::Off
        },
        functional: config.xc.functional.engine(),
        cutoff_hartree: config.basis.ecutwfc_ry / 2.0,
        grid_override: if config.symmetry.mode == crate::config::SymmetryMode::Auto {
            Some(config.prepare().map_err(|e| error(1, e))?.grid.dims())
        } else {
            config.basis.fft_grid
        },
        pw92_prefactor: config.xc.prefactor().map_err(|e| error(1, e))?,
        local_mesh: config.pseudo.local_mesh().map_err(|e| error(1, e))?,
        radial_even_rule: config.pseudo.radial_even_rule().map_err(|e| error(1, e))?,
        allow_functional_override: config.xc.allow_functional_override,
        standard: Some(dft_engine::StandardProfile {
            threshold_ha: profile.threshold_ry / 2.0,
        }),
        mixing: dft_engine::Mixing::Pulay {
            beta: profile.beta,
            history: profile.history,
        },
        ..ScfOptions::default()
    };
    if let Some(s) = s {
        o.threads = s.threads.unwrap_or(1);
        o.bands = s.nbands;
        o.extra_bands = s.extra_bands.unwrap_or(o.extra_bands);
        o.seed = s.seed.unwrap_or(o.seed);
        o.max_iterations = usize::try_from(s.max_iterations.unwrap_or(o.max_iterations as i64))
            .map_err(|_| error(1, "max_iterations must be positive"))?;
        o.verification = if s.verify().map_err(|e| error(1, e))? == Some(false) {
            dft_engine::Verification::Off
        } else {
            dft_engine::Verification::FixedPoint
        };
    }
    if o.max_iterations == 0 {
        return Err(error(1, "max_iterations must be positive"));
    }
    Ok(o)
}
/// `options` with the reference switch given explicitly (`true` = the strict
/// reference path of `BAANDI_REFERENCE_STRICT=on`): the options a child
/// process with that environment resolves, for in-process comparisons.
pub fn options_with_reference(
    config: &Config,
    reference_strict: bool,
) -> Result<ScfOptions, RunError> {
    let path = config
        .scf_path_with(reference_strict)
        .map_err(|e| error(1, e.to_string()))?;
    if reference_strict {
        oracle_options(config, path)
    } else {
        reference_options(config)
    }
}
fn oracle_options(config: &Config, path: crate::config::ScfPath) -> Result<ScfOptions, RunError> {
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
        ..ScfOptions::strict_reference()
    };
    if let Some(s) = &config.scf {
        if crate::parallel::size() >= 4
            && dft_numerics::switches::choice(dft_numerics::switches::Id::M2_CAPACITY) == "auto"
        {
            o.davidson_capacity_per_band = 3;
        }
        o.start_wavefunctions = s
            .start_wavefunctions()
            .map_err(|e| error(1, e.to_string()))?;
        o.start_density = s.start_density().map_err(|e| error(1, e.to_string()))?;
        o.gamma_real = s.gamma_real().map_err(|e| error(1, e.to_string()))?;
        o.threads = s.threads.unwrap_or(1);
        if let Some(n) = s.max_iterations {
            o.max_iterations =
                usize::try_from(n).map_err(|_| error(1, "max_iterations must be positive"))?;
        }
        o.bands = s.nbands;
        if let Some(n) = s.extra_bands {
            o.extra_bands = n;
        }
        if let Some(x) = s.seed {
            o.seed = x;
        }
        o.mixing = s.mixing_strategy().map_err(|e| error(1, e.to_string()))?;
    }
    o.standard = path.standard;
    o.verification = path.verification;
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
    // Default mixing coefficient by class (dev08 mixing study, ledger MIXING-STUDY):
    // strict Gamma-only runs 0.3 (si32 -10.2 percent, si16 -4.3 percent on both
    // protocols, 10 of 10 geometries), everything else 0.5 (standard class: si8 +23 percent
    // at 0.3; k meshes and multi-species cells 4 to 8 percent better at 0.5).
    // An explicit mixing_beta, or a mixing method other than pulay, is never overridden.
    let explicit = config.scf.as_ref().is_some_and(|s| {
        s.mixing_beta.is_some() || s.mixing.as_deref().is_some_and(|m| m != "pulay")
    });
    if !explicit {
        let single_k = config.resolved_points().is_ok_and(|p| p.len() == 1);
        let beta = if o.standard.is_none() && single_k {
            0.3
        } else {
            0.5
        };
        if let dft_engine::Mixing::Pulay { history, .. } = o.mixing {
            o.mixing = dft_engine::Mixing::Pulay { beta, history };
        }
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
/// `single_path` block of every result (ADR 0043): which SCF path ran and why.
pub(crate) fn single_path_json(path: &crate::config::ScfPath, o: &ScfOptions) -> Value {
    json!({"status":path.status,"reason":path.reason,
        "effective_profile":if o.standard.is_some() { "standard" } else { "strict" },
        "standard_threshold_ha":o.standard.map(|p| p.threshold_ha),
        "threshold_source":path.threshold_source,
        "verification":if o.standard.is_some() { o.verification.name() } else { "strict_certificate" },
        "switch":"BAANDI_REFERENCE_STRICT","adr":"0043"})
}
/// Reason recorded for MPI runs with more than one k point (review F02).
pub(crate) const REPLICATED_MULTI_K: &str =
    "more than one k point: no k-point pools, every rank solves every k point";
/// Review F02: the parallel model that actually ran. `ranks > 1` without a
/// single distributed eigensolve means every rank repeated the whole
/// calculation (`mpi_replicated`, with the reason); such a run is not an MPI
/// measurement. One rank is `threads`.
pub(crate) fn parallel_model(
    ranks: usize,
    distributed_solves: usize,
    m2: bool,
    replicated_reason: Option<String>,
) -> (&'static str, Option<String>) {
    if ranks <= 1 {
        ("threads", None)
    } else if distributed_solves == 0 {
        (
            "mpi_replicated",
            Some(replicated_reason.unwrap_or_else(|| "no distributed eigensolve".into())),
        )
    } else if m2 {
        ("mpi_stick_rows", None)
    } else {
        ("mpi_row_block", None)
    }
}
/// What the solver actually did (filled from the SCF result).
pub(crate) struct Solver<'a> {
    pub gamma_real_used: bool,
    pub gamma_real_fallback_reason: Option<&'a str>,
    pub parallel_model: &'static str,
}
/// Review F01: requested and effective path, solver, verification, backend and
/// the switches whose value is not the default, so that a stored result (a
/// baseline) identifies the combination that produced it. The complete switch
/// table stays under `scf.switches`.
pub(crate) fn execution_json(
    config: &Config,
    path: &crate::config::ScfPath,
    o: &ScfOptions,
    s: Solver<'_>,
) -> Value {
    use dft_numerics::switches::{self, Id, Source, Value as V};
    let scf = config.scf.as_ref();
    let lean = dft_engine::lean_counters();
    let effective_solver = if dft_engine::parallel::distributed_solves() > 0 {
        if dft_engine::parallel::distributed_rows_solves() > 0 {
            "mpi_distributed_real_davidson_stick_rows"
        } else {
            "mpi_distributed_real_davidson_row_block"
        }
    } else if s.gamma_real_used {
        if lean.solves > 0 {
            "real_gamma_lean_davidson"
        } else {
            "real_gamma_classic_davidson"
        }
    } else {
        "complex_davidson"
    };
    let non_default: serde_json::Map<_, _> = switches::manifest()
        .into_iter()
        .filter(|e| e.source != Source::Default)
        .map(|e| {
            let value = match e.value {
                V::Choice(c) => json!(c),
                V::Number(x) => json!(x),
                V::Raw(r) => json!(r),
                V::Present(b) => json!(b),
            };
            (
                e.switch.name.to_owned(),
                json!({"value":value,"source":e.source.as_str()}),
            )
        })
        .collect();
    json!({
        "requested":{
            "standard_threshold_ha":scf.and_then(|s| s.standard_threshold_ha),
            "verify":scf.and_then(|s| s.verify.as_deref()),
            "gamma_real":scf.and_then(|s| s.gamma_real.as_deref()).unwrap_or("auto"),
            "threads":scf.and_then(|s| s.threads),
            "reference_strict":crate::driver_select::requested(config).strict(),
        },
        "effective_profile":if o.standard.is_some() { "standard" } else { "strict" },
        "path_status":path.status,
        "path_reason":path.reason,
        "standard_threshold_ha":o.standard.map(|p| p.threshold_ha),
        "verification":if o.standard.is_some() { o.verification.name() } else { "strict_certificate" },
        "effective_solver":effective_solver,
        "real_gamma":s.gamma_real_used,
        "real_gamma_fallback_reason":s.gamma_real_fallback_reason,
        "davidson_mode":switches::choice(Id::DAVIDSON),
        "lean_solves":lean.solves,
        "solver_note":"lean applies to Ritz-value (standard-path) real-Gamma solves; residual-criterion solves (strict path, certificate) use the classic real Davidson; process totals including the cutoff-ladder coarse stage",
        "backend":{
            "linear_algebra":if dft_numerics::MKL_BACKEND { "mkl" } else { "pure" },
            "mpi_feature":cfg!(feature = "mpi"),
            "ranks":crate::parallel::size(),
            "threads":o.threads.max(1),
        },
        "parallel_model":s.parallel_model,
        "non_default_switches":non_default,
    })
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
        let (path, fallback) = e.note.map_or((None, None), |n| (Some(n.path), n.fallback));
        m.insert(
            e.switch.name.into(),
            json!({"value":value,"source":e.source.as_str(),"default":e.switch.default,
                "allowed":e.switch.allowed,"meaning":e.switch.meaning,
                "requested":e.requested,"consulted":e.consulted,
                "path":path,"fallback_reason":fallback}),
        );
    }
    Value::Object(m)
}

/// RUN-CONFIG: the path a consulted switch actually took, derived from the
/// engine's counters after the run (process totals: one run per CLI process),
/// with the reason when it did not act as requested.
pub(crate) fn note_switch_paths(rc: dft_numerics::switches::RunConfig) {
    use dft_numerics::switches::{Id, RunCounter};
    let ranks = crate::parallel::size();
    let dist = dft_engine::parallel::distributed_solves();
    let rows = dft_engine::parallel::distributed_rows_solves();
    if rc.consulted(Id::MPI_M2) {
        let (path, why) = match (dist, rows) {
            (0, _) => (
                "none",
                Some(if ranks <= 1 {
                    "one rank (threads path)"
                } else {
                    "no distributed eigensolve"
                }),
            ),
            (_, 0) => (
                "m1_row_block",
                rc.on(Id::MPI_M2).then_some("M2 not eligible"),
            ),
            _ => ("m2_stick_rows", None),
        };
        rc.note_path(Id::MPI_M2, path, why);
    }
    if rc.consulted(Id::M2_GRID) {
        let [potential, _, _, _, declined] = dft_engine::grid_dist::counters();
        let (path, why) = if potential > 0 {
            (
                "distributed",
                (declined > 0).then_some("some maps not eligible (replicated)"),
            )
        } else if !rc.on(Id::M2_GRID) {
            ("replicated", None)
        } else if rows == 0 {
            (
                "replicated",
                Some("not the M2 path (one rank or row-block path)"),
            )
        } else {
            (
                "replicated",
                Some("not eligible (functional, core density, density floor or symmetry)"),
            )
        };
        rc.note_path(Id::M2_GRID, path, why);
    }
    if rc.consulted(Id::M2_MIX) {
        let [_, _, built] = dft_engine::mix_dist::counters();
        let (path, why) = if built > 0 {
            ("distributed", None)
        } else if !rc.on(Id::M2_MIX) {
            ("replicated", None)
        } else if rows == 0 {
            (
                "replicated",
                Some("not the M2 path (one rank or row-block path)"),
            )
        } else {
            (
                "replicated",
                Some("mixer not eligible (Hartree metric, Kerker or linear)"),
            )
        };
        rc.note_path(Id::M2_MIX, path, why);
    }
    // M2-OWNED: the owned modes' path and the reasons of every decline
    // (`owned.fields_declined.<reason>` / `owned.density_declined.<reason>`).
    let owned: std::collections::BTreeMap<String, usize> =
        dft_engine::grid_dist::owned_stats().0.into_iter().collect();
    for (id, used_key, prefix) in [
        (
            Id::M2_LOCAL_FIELDS,
            "owned.fields_owned_maps",
            "owned.fields_declined.",
        ),
        (
            Id::M2_LOCAL_DENSITY,
            "owned.density_owned_maps",
            "owned.density_declined.",
        ),
    ] {
        if !rc.consulted(id) {
            continue;
        }
        let used = owned.get(used_key).copied().unwrap_or(0);
        let reasons: Vec<&str> = owned
            .iter()
            .filter(|(k, v)| **v > 0 && k.starts_with(prefix))
            .map(|(k, _)| &k[prefix.len()..])
            .collect();
        let reason = (!reasons.is_empty()).then(|| reasons.join(","));
        let path = if used > 0 { "owned" } else { "gathered" };
        rc.note_path(id, path, reason.as_deref().filter(|_| rc.on(id)));
    }

    if rc.consulted(Id::REUSE_PSI) {
        let r = dft_engine::gamma_real::reuse_psi_counters();
        let (path, why) = if r.armed > 0 {
            ("kept_start", None)
        } else if rc.on(Id::REUSE_PSI) {
            (
                "plain_start",
                Some("no warm lean real Gamma solve took the kept start"),
            )
        } else {
            ("plain_start", None)
        };
        rc.note_path(Id::REUSE_PSI, path, why);
    }
    if rc.consulted(Id::EWALD_CACHE) {
        let hits = rc.counters()[RunCounter::EwaldHits as usize];
        rc.note_path(
            Id::EWALD_CACHE,
            if hits > 0 { "hit" } else { "computed" },
            None,
        );
    }
}

/// Hash with bounded scratch rather than allocating/copying the entire binary.
/// Read errors retain the previous caller's `None` identity behavior.
fn binary_sha256(path: &Path) -> io::Result<String> {
    use io::Read;
    let mut file = fs::File::open(path)?;
    let mut digest = Sha256::new();
    let mut buffer = [0_u8; 64 * 1024];
    loop {
        let size = match file.read(&mut buffer) {
            Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
            other => other?,
        };
        if size == 0 {
            break;
        }
        digest.update(&buffer[..size]);
    }
    Ok(format!("{:x}", digest.finalize()))
}

/// Optional cache of a content digest, never an alternative provenance value.
/// Unix metadata is a coherence key, not authentication. The executable and its
/// sibling sidecar must live in a trusted build directory; the measurement
/// master independently hashes the executable before and after each job.
#[cfg(unix)]
mod binary_hash_cache {
    use super::*;
    use io::{Read, Write};
    use std::os::unix::fs::{MetadataExt, OpenOptionsExt};
    use std::sync::atomic::{AtomicU64, Ordering};

    const SCHEMA: &str = "baandi-binary-sha256-2";
    const MAX_RECORD: u64 = 4096;
    // Two seconds covers whole-second ctime and coarse filesystem timestamps.
    // An entry recorded in the same tick remains suspect even after aging.
    const SAFETY_NS: i128 = 2_000_000_000;
    fn now_ns() -> io::Result<u64> {
        SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .map_err(io::Error::other)
            .and_then(|d| u64::try_from(d.as_nanos()).map_err(io::Error::other))
    }
    fn racy(stamp: &Stamp, recorded: u64, now: u64) -> bool {
        [
            (stamp.modified_s, stamp.modified_ns),
            (stamp.changed_s, stamp.changed_ns),
        ]
        .into_iter()
        .any(|(s, ns)| {
            let time = i128::from(s) * 1_000_000_000 + i128::from(ns);
            time >= i128::from(recorded) - SAFETY_NS || time >= i128::from(now) - SAFETY_NS
        }) || recorded > now
    }

    #[derive(Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
    #[serde(deny_unknown_fields)]
    struct Stamp {
        device: u64,
        inode: u64,
        bytes: u64,
        modified_s: i64,
        modified_ns: i64,
        changed_s: i64,
        changed_ns: i64,
        owner: u32,
        mode: u32,
    }
    impl Stamp {
        fn of(m: &fs::Metadata) -> Self {
            Self {
                device: m.dev(),
                inode: m.ino(),
                bytes: m.len(),
                modified_s: m.mtime(),
                modified_ns: m.mtime_nsec(),
                changed_s: m.ctime(),
                changed_ns: m.ctime_nsec(),
                owner: m.uid(),
                mode: m.mode(),
            }
        }
    }
    #[derive(serde::Serialize, serde::Deserialize)]
    #[serde(deny_unknown_fields)]
    struct Record {
        schema: String,
        stamp: Stamp,
        sha256: String,
        recorded_ns: u64,
        checksum: String,
    }
    fn checksum(stamp: &Stamp, digest: &str, recorded_ns: u64) -> String {
        // Detect accidental edits to even a syntactically valid digest; this is
        // deliberately not claimed to authenticate a malicious sidecar writer.
        sha256(&serde_json::to_vec(&(stamp, digest, recorded_ns)).expect("integer/string record"))
    }
    pub(super) fn sidecar(path: &Path) -> PathBuf {
        let mut name = path.as_os_str().to_owned();
        name.push(".baandi-sha256.json");
        PathBuf::from(name)
    }
    fn trusted_record(m: &fs::Metadata, owner: u32) -> bool {
        m.is_file() && m.uid() == owner && m.mode() & 0o022 == 0
    }
    fn read_record(path: &Path, stamp: &Stamp, now: u64) -> Option<String> {
        let link = fs::symlink_metadata(path).ok()?;
        if !trusted_record(&link, stamp.owner) || link.len() > MAX_RECORD {
            return None;
        }
        let file = fs::File::open(path).ok()?;
        let meta = file.metadata().ok()?;
        if Stamp::of(&link) != Stamp::of(&meta) {
            return None;
        }
        let mut bytes = Vec::new();
        file.take(MAX_RECORD + 1).read_to_end(&mut bytes).ok()?;
        if bytes.len() as u64 > MAX_RECORD {
            return None;
        }
        let record: Record = serde_json::from_slice(&bytes).ok()?;
        if record.schema != SCHEMA
            || racy(stamp, record.recorded_ns, now)
            || record.stamp != *stamp
            || record.sha256.len() != 64
            || !record
                .sha256
                .bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
            || record.checksum != checksum(&record.stamp, &record.sha256, record.recorded_ns)
        {
            return None;
        }
        Some(record.sha256)
    }
    fn publish(path: &Path, record: &Record) -> io::Result<()> {
        // Never follow or replace a symlink, directory or foreign/writable file.
        match fs::symlink_metadata(path) {
            Ok(m) if !trusted_record(&m, record.stamp.owner) => {
                return Err(io::Error::other("untrusted hash sidecar"));
            }
            Err(e) if e.kind() != io::ErrorKind::NotFound => return Err(e),
            _ => {}
        }
        static NEXT: AtomicU64 = AtomicU64::new(0);
        let mut name = path.as_os_str().to_owned();
        name.push(format!(
            ".tmp.{}.{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        let temporary = PathBuf::from(name);
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .mode(0o600)
            .open(&temporary)?;
        let result = (|| {
            let bytes = serde_json::to_vec(record)?;
            file.write_all(&bytes)?;
            drop(file);
            fs::rename(&temporary, path)
        })();
        if result.is_err() {
            // Only the temporary we successfully created is eligible for cleanup.
            let _ = fs::remove_file(&temporary);
        }
        result
    }
    pub(super) struct Outcome {
        pub digest: String,
        pub status: &'static str,
        pub hashed_bytes: u64,
        pub published: bool,
    }
    fn stable(file: &fs::File, path: &Path, stamp: &Stamp) -> io::Result<bool> {
        Ok(Stamp::of(&file.metadata()?) == *stamp && Stamp::of(&fs::metadata(path)?) == *stamp)
    }
    pub(super) fn digest(path: &Path) -> io::Result<Outcome> {
        digest_at(path, now_ns()?)
    }
    fn digest_at(path: &Path, now: u64) -> io::Result<Outcome> {
        let mut file = fs::File::open(path)?;
        let meta = file.metadata()?;
        if !meta.is_file() {
            return Err(io::Error::other("executable is not a regular file"));
        }
        let stamp = Stamp::of(&meta);
        let cache = sidecar(path);
        if let Some(digest) = read_record(&cache, &stamp, now)
            && stable(&file, path, &stamp)?
        {
            return Ok(Outcome {
                digest,
                status: "hit",
                hashed_bytes: 0,
                published: false,
            });
        }
        let mut hash = Sha256::new();
        let mut buffer = [0_u8; 64 * 1024];
        let mut hashed_bytes = 0;
        loop {
            let n = match file.read(&mut buffer) {
                Err(e) if e.kind() == io::ErrorKind::Interrupted => continue,
                other => other?,
            };
            if n == 0 {
                break;
            }
            hash.update(&buffer[..n]);
            hashed_bytes += n as u64;
        }
        if hashed_bytes != stamp.bytes || !stable(&file, path, &stamp)? {
            return Err(io::Error::other("executable changed while hashing"));
        }
        let digest = format!("{:x}", hash.finalize());
        let record = Record {
            schema: SCHEMA.into(),
            checksum: checksum(&stamp, &digest, now),
            recorded_ns: now,
            stamp,
            sha256: digest.clone(),
        };
        let published = publish(&cache, &record).is_ok();
        Ok(Outcome {
            digest,
            status: "miss",
            hashed_bytes,
            published,
        })
    }
    #[cfg(test)]
    mod tests {
        use super::*;
        use std::os::unix::fs::PermissionsExt;

        struct Fixture(PathBuf);
        impl Fixture {
            fn new(label: &str) -> Self {
                static NEXT: AtomicU64 = AtomicU64::new(0);
                let dir = std::env::temp_dir().join(format!(
                    "hash-cache-{label}-{}-{}",
                    std::process::id(),
                    NEXT.fetch_add(1, Ordering::Relaxed)
                ));
                fs::create_dir(&dir).unwrap();
                Self(dir)
            }
            fn binary(&self) -> PathBuf {
                self.0.join("binary")
            }
        }
        impl Drop for Fixture {
            fn drop(&mut self) {
                fs::remove_dir_all(&self.0).unwrap();
            }
        }
        #[test]
        fn cold_hit_and_same_length_rewrite_with_restored_mtime() {
            let f = Fixture::new("rewrite");
            let path = f.binary();
            fs::write(&path, b"original").unwrap();
            let modified = fs::metadata(&path).unwrap().modified().unwrap();
            let cold = digest(&path).unwrap();
            assert_eq!(cold.digest, sha256(b"original"));
            assert_eq!(cold.hashed_bytes, 8);
            assert_eq!(cold.status, "miss");
            assert!(cold.published);
            // Immediate re-read is deliberately unverified; a later rehash
            // records a safe timestamp and enables ordinary cache hits.
            assert_eq!(digest(&path).unwrap().status, "miss");
            let later = now_ns().unwrap() + 3_000_000_000;
            assert_eq!(digest_at(&path, later).unwrap().status, "miss");
            assert_eq!(digest_at(&path, later).unwrap().status, "hit");
            for i in 0..50 {
                let bytes: &[u8] = if i % 2 == 0 { b"replaced" } else { b"original" };
                fs::write(&path, bytes).unwrap();
                fs::File::options()
                    .write(true)
                    .open(&path)
                    .unwrap()
                    .set_times(fs::FileTimes::new().set_modified(modified))
                    .unwrap();
                // Force an exact metadata collision with a stale digest. This
                // does not depend on the host ctime resolution or scheduling.
                let stamp = Stamp::of(&fs::metadata(&path).unwrap());
                let recorded_ns = now_ns().unwrap();
                let stale = sha256(if i % 2 == 0 { b"original" } else { b"replaced" });
                let record = Record {
                    schema: SCHEMA.into(),
                    checksum: checksum(&stamp, &stale, recorded_ns),
                    stamp,
                    recorded_ns,
                    sha256: stale,
                };
                publish(&sidecar(&path), &record).unwrap();
                let changed = digest(&path).unwrap();
                assert_eq!(changed.status, "miss");
                assert_eq!(changed.digest, sha256(bytes));
                assert_eq!(changed.hashed_bytes, 8);
            }
        }
        #[test]
        fn safety_window_checks_both_clocks_and_future_timestamps() {
            let f = Fixture::new("clock");
            fs::write(f.binary(), b"x").unwrap();
            let mut stamp = Stamp::of(&fs::metadata(f.binary()).unwrap());
            stamp.modified_s = 10;
            stamp.modified_ns = 0;
            stamp.changed_s = 10;
            stamp.changed_ns = 0;
            assert!(!racy(&stamp, 13_000_000_000, 14_000_000_000));
            assert!(racy(&stamp, 12_000_000_000, 14_000_000_000));
            assert!(racy(&stamp, 13_000_000_000, 11_000_000_000));
            stamp.changed_s = 20;
            assert!(racy(&stamp, 13_000_000_000, 14_000_000_000));
        }
        #[test]
        fn replaced_inode_size_and_mode_invalidate() {
            let f = Fixture::new("identity");
            let path = f.binary();
            fs::write(&path, b"original").unwrap();
            digest(&path).unwrap();
            let replacement = f.0.join("replacement");
            fs::write(&replacement, b"replaced").unwrap();
            fs::rename(&replacement, &path).unwrap();
            assert_eq!(digest(&path).unwrap().status, "miss");
            fs::set_permissions(&path, fs::Permissions::from_mode(0o700)).unwrap();
            assert_eq!(digest(&path).unwrap().status, "miss");
            fs::write(&path, b"different length").unwrap();
            let out = digest(&path).unwrap();
            assert_eq!(out.status, "miss");
            assert_eq!(out.digest, sha256(b"different length"));
        }
        #[test]
        fn corrupt_oversized_and_writable_records_fall_back() {
            let f = Fixture::new("corrupt");
            let path = f.binary();
            fs::write(&path, b"payload").unwrap();
            digest(&path).unwrap();
            let cache = sidecar(&path);
            let original = fs::read(&cache).unwrap();
            let mut bad: Value = serde_json::from_slice(&original).unwrap();
            bad["sha256"] = json!("0".repeat(64)); // syntax valid, checksum invalid
            for bytes in [
                b"{".to_vec(),
                vec![b'x'; 4097],
                serde_json::to_vec(&bad).unwrap(),
            ] {
                fs::write(&cache, bytes).unwrap();
                let out = digest(&path).unwrap();
                assert_eq!(out.status, "miss");
                assert_eq!(out.digest, sha256(b"payload"));
                assert!(out.published);
            }
            fs::set_permissions(&cache, fs::Permissions::from_mode(0o666)).unwrap();
            let out = digest(&path).unwrap();
            assert_eq!(out.status, "miss");
            assert!(!out.published);
            assert_eq!(out.digest, sha256(b"payload"));
        }
        #[test]
        fn symlink_and_directory_sidecars_are_not_changed() {
            let f = Fixture::new("fallback");
            let path = f.binary();
            fs::write(&path, b"payload").unwrap();
            let target = f.0.join("unrelated");
            fs::write(&target, b"leave this alone").unwrap();
            let cache = sidecar(&path);
            std::os::unix::fs::symlink(&target, &cache).unwrap();
            let out = digest(&path).unwrap();
            assert_eq!(out.digest, sha256(b"payload"));
            assert!(!out.published);
            assert!(
                fs::symlink_metadata(&cache)
                    .unwrap()
                    .file_type()
                    .is_symlink()
            );
            assert_eq!(fs::read(&target).unwrap(), b"leave this alone");
            fs::remove_file(&cache).unwrap();
            fs::create_dir(&cache).unwrap();
            assert!(!digest(&path).unwrap().published);
            assert!(cache.is_dir());
            fs::remove_file(path).unwrap();
            assert!(digest(&f.binary()).is_err());
        }
        #[test]
        fn concurrent_cold_writers_publish_one_complete_record() {
            let f = Fixture::new("concurrent");
            let path = f.binary();
            let payload = vec![17; 131_081];
            fs::write(&path, &payload).unwrap();
            let barrier = std::sync::Arc::new(std::sync::Barrier::new(8));
            let threads: Vec<_> = (0..8)
                .map(|_| {
                    let (path, barrier) = (path.clone(), barrier.clone());
                    std::thread::spawn(move || {
                        barrier.wait();
                        digest(&path).unwrap().digest
                    })
                })
                .collect();
            for t in threads {
                assert_eq!(t.join().unwrap(), sha256(&payload));
            }
            // Concurrent publishers can all finish in the binary's ctime tick.
            // Their complete records remain unverified until a safe rehash.
            let later = now_ns().unwrap() + 3_000_000_000;
            assert_eq!(digest_at(&path, later).unwrap().status, "miss");
            assert_eq!(digest_at(&path, later).unwrap().status, "hit");
            assert_eq!(fs::read_dir(&f.0).unwrap().count(), 2);
        }
    }
}

/// SHA-256 of the running executable (fresh by default, cached only on request).
pub(crate) fn binary_identity() -> Value {
    // Do not use the process OnceLock for the opt-in path: library callers may
    // change scopes, and every cached use must revalidate the file metadata.
    if dft_numerics::switches::on(dft_numerics::switches::Id::BINARY_HASH_CACHE) {
        #[cfg(unix)]
        {
            let hash_clock = launch::clock();
            let path = std::env::current_exe().ok();
            let result = path
                .as_ref()
                .and_then(|p| binary_hash_cache::digest(p).ok());
            launch::elapsed("binary_sha256", hash_clock);
            return json!({"path":path.map(|p| p.display().to_string()),
                "sha256":result.as_ref().map(|r| &r.digest),
                "hash_cache":{"status":result.as_ref().map_or("unavailable", |r| r.status),
                    "hashed_bytes":result.as_ref().map(|r| r.hashed_bytes),
                    "published":result.as_ref().is_some_and(|r| r.published)}});
        }
        // Non-Unix platforms retain the original fresh digest, not an unchecked cache.
    }
    type Identity = (Option<String>, Option<String>);
    static ID: std::sync::OnceLock<Identity> = std::sync::OnceLock::new();
    let hash_clock = launch::clock();
    let (path, sha) = ID.get_or_init(|| {
        let path = std::env::current_exe().ok();
        let sha = path.as_ref().and_then(|p| binary_sha256(p).ok());
        (path.map(|p| p.display().to_string()), sha)
    });
    launch::elapsed("binary_sha256", hash_clock);
    json!({"path":path,"sha256":sha})
}

/// RUN-REPORT (`run` block of results.json): schema, run ID, exit reason, the
/// binary identity, the switch agreement hash, the run-local counters and the
/// completion marker. A results.json without `run.completed = true` is an
/// interrupted or foreign file (scripts/parse_results.py rejects it).
pub(crate) fn run_block(rc: dft_numerics::switches::RunConfig, exit_reason: &str) -> Value {
    let counters: serde_json::Map<_, _> = dft_numerics::switches::RUN_COUNTER_LABELS
        .iter()
        .zip(rc.counters())
        .map(|(k, v)| ((*k).to_string(), json!(v)))
        .collect();
    json!({"schema":"run-report-1","run_id":rc.run_id_hex(),"exit_reason":exit_reason,
        "source_commit":env!("DFT_GIT_COMMIT"),"binary":binary_identity(),
        "rank":crate::parallel::rank(),"ranks":crate::parallel::size(),
        "switch_agreement_hash":format!("{:016x}", rc.agreement_hash()),
        "counters":counters,"completed":true})
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
/// P48F: additive timing blocks that close the time budget of a run.
/// total = setup + phases + certificate outside phases + step overhead +
/// remainder (session start and result assembly).
fn closed_timings(
    phases: &PhaseTimings,
    history: &[IterationRecord],
    setup: f64,
    total: f64,
    clock: dft_engine::phase_clock::PhaseClock,
) -> Value {
    let fine = clock;
    let cert_phases: f64 = history
        .iter()
        .filter_map(|h| h.fixed_point_certificate.as_ref())
        .map(|c| phase_sum(&c.timings))
        .fold(0.0, |a, b| a + b);
    let cert_wall = fine.certificate.as_secs_f64();
    let cert_outside = (cert_wall - cert_phases).max(0.0);
    let phases_s = phase_sum(phases);
    let steps = fine.steps.as_secs_f64();
    let step_overhead = steps - phases_s - cert_outside;
    let remainder = total - setup - steps;
    json!({
        "certificate":{"blocks":fine.certificate_blocks,"wall_s":cert_wall,"screen_s":fine.screen.as_secs_f64(),
            "phases_s":cert_phases,"outside_phases_s":cert_outside,
            "note":"wall of the certificate blocks (screen, re-solve, checks); phases_s is already inside timings.phases"},
        "atomic_start_s":fine.atomic_start.as_secs_f64(),
        "steps_s":steps,
        "closure":{"setup_s":setup,"phases_s":phases_s,"certificate_outside_phases_s":cert_outside,
            "step_overhead_s":step_overhead,"remainder_s":remainder,"total_s":total,
            "note":"total = setup + phases + certificate_outside_phases + step_overhead + remainder; atomic_start_s is inside phases (potential of the first map)"}
    })
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
    if let Some(p) = &r.pulay {
        // SCF-TRAJECTORY: what the mixer did with this map's pair (diagnostic).
        result["pulay"] = json!({"differences":p.differences,"kept":p.kept,
            "gram_condition":p.gram_condition.is_finite().then_some(p.gram_condition),
            "coefficients":p.coefficients,"residual_l2_sum":p.residual_norm,
            "predicted_ratio":p.predicted_ratio,"step_l2_sum":p.step_norm,
            "coulomb_predicted_ratio":p.coulomb_predicted_ratio,"coulomb_refit":p.coulomb_refit});
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
    let report = execute_mode_inner(input, output_override, mode);
    launch::since("run_teardown", "outputs_complete");
    launch::mark("run_return");
    report
}

fn execute_mode_inner(
    input: &Path,
    output_override: Option<&Path>,
    mode: &crate::telemetry::RunMode,
) -> Result<RunReport, RunError> {
    mode.validate()?;
    let start = Instant::now();
    let mut setup_seconds = PhaseSeconds::new();
    let mut output_seconds = PhaseSeconds::new();
    let phase_clock = Instant::now();
    launch::mark("run_entry");
    let input_clock = launch::clock();
    let bytes = launch::timed("input_read", || fs::read(input))
        .map_err(|e| error(5, format!("{}: {e}", input.display())))?;
    let text = std::str::from_utf8(&bytes).map_err(|e| error(1, e))?;
    let config = launch::timed("input_parse", || Config::parse(text)).map_err(|e| error(1, e))?;
    launch::elapsed("input_read_parse", input_clock);
    setup_seconds.insert("input_read_parse", phase_clock.elapsed().as_secs_f64());
    let config_clock = launch::clock();
    // RUN-CONFIG: this run's switch configuration (the process object in the
    // CLI binary; a thread scope for in-process callers with other values).
    let (run_config, _run_scope) = crate::run_config::for_run(&config).map_err(|e| error(1, e))?;
    launch::elapsed("run_config", config_clock);
    let routing_mode = if mode.agent {
        crate::driver_select::Mode::Agent
    } else if mode.observe {
        crate::driver_select::Mode::Observe
    } else {
        crate::driver_select::Mode::Classic
    };
    let selection = crate::driver_select::select(
        &config,
        routing_mode,
        crate::driver_select::Command::Run,
        &[],
    );
    let reference = selection.reference();
    let input_digest =
        reference.then(|| phase(true, &mut setup_seconds, "input_hash", || sha256(&bytes)));
    // Phase 0 only: no storage or additional clock reads when profiling is off.
    // Keep rank-local file records; do not add collectives to measure redundant I/O.
    let mut upf_records = dft_engine::profile::enabled().then(Vec::new);
    dft_engine::profile::davidson_lane::reset();
    if !reference
        && (config.collinear()
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
    if !reference
        && config.symmetry.mode == crate::config::SymmetryMode::Auto
        && (mode.observe || mode.agent)
    {
        return Err(error(
            1,
            "symmetry.execution: checkpoint identity does not carry symmetry yet; use Classic execution",
        ));
    }
    let preflight_clock = launch::clock();
    let mut preflight = phase(
        reference,
        &mut setup_seconds,
        "preflight_cell_basis_upf",
        || crate::config::preflight::check(&config, input),
    )
    .map_err(|e| error(1, e))?;
    launch::elapsed("preflight_including_prepare_and_upf", preflight_clock);
    let prepare_clock = launch::clock();
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
    // INIT-CLI: reuse the preflight's `Config::prepare` result (built from the
    // same configuration) instead of building the bases and the grid again.
    launch::prepared(preflight.prepared.is_some());
    let prepared = match preflight.prepared.take() {
        Some(prepared) => prepared,
        None => config.prepare().map_err(|e| error(1, e))?,
    };
    launch::elapsed("prepare_reuse_and_preflight_print", prepare_clock);
    // ADR 0043 ONEPATH-4: the observe/agent runtime modes run the same single
    // standard path with the verification as a classic run (checkpoint v1
    // carries the verification option and its attempts).
    let opts = options_for(&config, selection)?;
    let mut potentials = Vec::new();
    let mut species = Vec::new();
    let mut inputs = vec![fs::canonicalize(input)?];
    let upf_clock = launch::clock();
    for (name, path) in &config.pseudopotentials {
        let path = input.parent().unwrap_or(Path::new(".")).join(path);
        let read_clock = upf_records.as_ref().map(|_| Instant::now());
        let data = phase(reference, &mut setup_seconds, "upf_read", || {
            fs::read(&path)
        })
        .map_err(|e| error(5, format!("{}: {e}", path.display())))?;
        let read_s = read_clock.map(|t| t.elapsed().as_secs_f64());
        let parse_clock = upf_records.as_ref().map(|_| Instant::now());
        let upf_text = std::str::from_utf8(&data).map_err(|e| error(2, e))?;
        let p = phase(reference, &mut setup_seconds, "upf_parse", || {
            parse_upf(upf_text)
        })
        .map_err(|e| error(2, format!("{}: {e}", path.display())))?;
        let parse_s = parse_clock.map(|t| t.elapsed().as_secs_f64());
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
        let hash_clock = upf_records.as_ref().map(|_| Instant::now());
        let digest = phase(reference, &mut setup_seconds, "upf_hash", || sha256(&data));
        let hash_s = hash_clock.map(|t| t.elapsed().as_secs_f64());
        if let Some(records) = &mut upf_records {
            records.push(
                json!({"species":name,"bytes":data.len(),"projectors":p.projectors().len(),
                "read_s":read_s,"parse_s":parse_s,"hash_s":hash_s}),
            );
        }
        species.push(json!({"name":name,"path":path,"sha256":digest,"z_valence":p.z_valence(),"functional":p.functional(),"n_projectors":p.projectors().len()}));
        inputs.push(path);
        potentials.push((name.clone(), p));
    }
    let upfs: Vec<_> = potentials.iter().map(|(_, p)| p).collect();
    let selection = crate::driver_select::select(
        &config,
        routing_mode,
        crate::driver_select::Command::Run,
        &upfs,
    );
    selection.announce();
    let scf_path = config
        .scf_path_selected(selection)
        .map_err(|e| error(1, e))?;
    let reference = selection.reference();
    let opts = options_for(&config, selection)?;
    if !reference
        && (mode.observe || mode.agent)
        && (config.collinear()
            || config
                .finite_temperature()
                .map_err(|e| error(1, e))?
                .is_some()
            || config.symmetry.mode == crate::config::SymmetryMode::Auto)
    {
        return Err(error(
            1,
            "runtime control for spin/FD or symmetry requires an eligible reference driver",
        ));
    }
    launch::elapsed("upf_read_parse_hash", upf_clock);
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
    if reference
        && (config.collinear()
            || config
                .finite_temperature()
                .map_err(|e| error(1, e))?
                .is_some())
    {
        return reference_electronic_run(
            &config,
            &prepared,
            &opts,
            &atoms,
            &refs,
            &points,
            &out,
            electron_target,
            preflight,
            species,
            &bytes,
            mode,
        );
    }
    // Collinear spin and Fermi-Dirac drivers: the standard path with the
    // verification since ONEPATH-3 (strict only behind BAANDI_REFERENCE_STRICT),
    // complex solver, never distributed (the real-Gamma path excludes both), so
    // MPI ranks replicate them (review F02).
    let driver_record = |reason: &'static str| -> Value {
        let (model, replicated) = parallel_model(
            crate::parallel::size(),
            0,
            false,
            Some(format!("{reason}: every rank runs the whole calculation")),
        );
        let mut v = json!({
            "single_path": single_path_json(&scf_path, &opts),
            "execution": execution_json(&config, &scf_path, &opts, Solver {
                gamma_real_used: false,
                gamma_real_fallback_reason: Some(reason),
                parallel_model: model,
            }),
        });
        if crate::parallel::size() > 1 {
            if crate::parallel::rank() == 0 {
                eprintln!(
                    "warning: MPI ranks replicate the whole calculation ({reason}); this run is not a distributed (MPI vs MPI) measurement"
                );
            }
            v["parallel"] = json!({"model":model,"replicated_reason":replicated,
                "ranks":crate::parallel::size(),"rank":crate::parallel::rank(),"distributed_solves":0});
        }
        v["scf"] = options_json(&opts);
        selection.record(&mut v);
        v
    };
    if let Some(theta) = config.finite_temperature().map_err(|e| error(1, e))? {
        let mut record = driver_record("Fermi-Dirac occupations");
        if let Some(reason) = selection.reason.filter(|_| !selection.strict()) {
            record["reference_driver_fallback_reason"] = json!(reason);
            record["execution"]["reference_driver_fallback_reason"] = json!(reason);
        }
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
            record,
        );
    }
    if let Some(sector) = config
        .spin_sector(electron_target)
        .map_err(|e| error(1, e))?
    {
        let record = driver_record("collinear spin");
        return spin_run::execute(
            &config, &prepared, opts, &atoms, &refs, &points, sector, &out, preflight, species,
            &bytes, record,
        );
    }
    // BREAKTHROUGH-4: the CLI's part of `init` (input, preflight, prepared
    // basis, pseudopotential reading) before the engine constructs the problem.
    if dft_engine::profile::enabled() {
        dft_engine::profile::init_part(dft_engine::profile::InitPart::Cli, Some(start));
    }
    launch::mark("engine_constructor_entry");
    let constructor_clock = launch::clock();
    // Reference initialization owns its final-cutoff state directly. Do not
    // construct the legacy problem (or its replicated band-density workspace).
    let mut problem = if reference {
        None
    } else {
        Some({
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
            let density_floor = problem.core_density().map(|_| 1e-10);
            if let Some(floor) = density_floor {
                problem.set_density_floor(floor).map_err(engine_error)?;
            }
            problem
        })
    };
    let density_floor = problem
        .as_ref()
        .and_then(|p| p.core_density().map(|_| 1e-10));
    let electrons = if let Some(p) = &problem {
        p.electron_count()
    } else {
        atoms
            .iter()
            .map(|a| {
                refs.iter()
                    .find(|p| p.label == a.species)
                    .expect("validated species")
                    .potential
                    .z_valence()
            })
            .sum()
    };
    launch::elapsed("engine_construction", constructor_clock);
    launch::mark("engine_ready");
    let setup = start.elapsed().as_secs_f64();
    // BREAKTHROUGH-0-INSTRUMENT: problem construction is the `init` bucket.
    dft_engine::profile::add_init(setup);
    let clock_start = dft_engine::phase_clock::snapshot();
    let reference_symmetry = if reference {
        reference_symmetry(&config)?
    } else {
        None
    };
    let mut reference_derivatives = None;
    let mut reference_symmetry_record = Value::Null;
    let mut reference_return = None;
    let mut reference_start_report = Value::Null;
    let mut reference_setup = Value::Null;
    let mut reference_output = Value::Null;
    let mut reference_records = Vec::new();
    let mut reference_status = None;
    let mut reference_summary = None;
    let mut reference_ladder = Value::Null;
    let mut reference_kept_hx = Value::Null;
    let mut reference_diagnostics = false;
    let mut reference_qe_clock = Vec::new();
    let mut reference_stages = None;
    if reference {
        dft_engine::refscf::driver::take_rmm_records();
        // Isolate the apply_h JSON from earlier in-process reference runs.
        dft_physics::distributed_real_gamma::take_apply_h_solves();
        if dft_numerics::switches::on(dft_numerics::switches::Id::DAVIDSON_STAGES) {
            dft_physics::distributed_real_gamma::ApplyH::take_davidson_work();
        }
    }
    let atomic_density_built = reference
        || problem
            .as_ref()
            .is_some_and(|p| p.atomic_superposition_density().is_some());
    // LAZY-RESULT: a classic run without forces, stress or orbital volumetric
    // output never reads the delivered orbitals, so an M2 result may keep them
    // distributed (no N x B gather on every rank). Replicated: every rank
    // parses the same input. Reset to the engine default after the SCF.
    let orbitals_needed = mode.observe
        || mode.agent
        || config.output.as_ref().is_some_and(|o| o.forces || o.stress)
        || volumetric_out::reads_orbitals(&config);
    struct RestoreDelivery;
    impl Drop for RestoreDelivery {
        fn drop(&mut self) {
            dft_engine::gamma_real::set_delivered_orbitals(true);
        }
    }
    let _restore_delivery = RestoreDelivery;
    dft_engine::gamma_real::set_delivered_orbitals(orbitals_needed);
    let scf_clock = launch::clock();
    launch::mark("scf_entry");
    let driver_entry = Instant::now();
    let driver_entry_s = driver_entry.duration_since(start).as_secs_f64();
    let mut outcome = if !reference && (mode.observe || mode.agent) {
        crate::telemetry::drive(problem.take().expect("legacy mode"), mode, &out)?
    } else {
        let run = if reference {
            let input = dft_engine::refstart::StartInput {
                cell: &prepared.cell,
                atoms: &atoms,
                potentials: &refs,
                options: &opts,
            };
            let profile = reference_profile(&config)?;
            let mut request = dft_engine::refscf::api::Request::new(profile, points.clone());
            request.symmetry = reference_symmetry;
            request.derivatives = dft_engine::refscf::api::Derivatives {
                forces: config.output.as_ref().is_some_and(|o| o.forces),
                stress: config.wants_stress(),
                deferred: true,
            };
            let delivered = if mode.observe || mode.agent {
                crate::telemetry::drive_reference(input, request, mode, &out)?
            } else {
                dft_engine::refscf::api::solve(input, request).map_err(engine_error)?
            };
            reference_derivatives = Some(delivered.derivatives);
            let result = delivered.run.into_scalar();
            result.map(|run| {
                reference_symmetry_record = run.symmetry.as_ref().map_or(Value::Null, reference_symmetry_json);
                reference_return = Some(Instant::now());
                reference_diagnostics = run.full_diagnostics;
                reference_qe_clock = run.qe_clock;
                reference_start_report = start_report_json(run.start_report);
                reference_setup = json!({
                    "schema":"reference-setup-1", "rank":crate::parallel::rank(),
                    "scope":"rank-local inclusive phases; atomic_block_total contains atomic_* children; collectives_nested overlaps phases; no barriers added",
                    "cli":setup_seconds,"cli_total":driver_entry_s,
                    "cli_other":driver_entry_s - setup_seconds.values().sum::<f64>(),
                    "preflight_scope":"includes cell/basis/grid and preflight UPF validation; upf_* is the subsequent engine parse; preflight internals are not separately timed",
                    "engine":boundary_json(&run.setup),
                    "run_entry_to_first_map":driver_entry_s + run.records.first().map_or(run.start_s, |r| r.map_start_s)
                });
                REFERENCE_LAUNCH.with(|v| {
                    let stamps = v.get();
                    reference_setup["main_to_run_entry"] = json!(stamps.map(|s| start.duration_since(s[0]).as_secs_f64()));
                    reference_setup["main_to_first_map"] = json!(stamps.map(|s| start.duration_since(s[0]).as_secs_f64() + driver_entry_s + run.records.first().map_or(run.start_s, |r| r.map_start_s)));
                    reference_setup["bootstrap"] = stamps.map_or(Value::Null, |s| json!({"mpi_init":(s[1]-s[0]).as_secs_f64(),"arguments":(s[2]-s[1]).as_secs_f64(),"switch_input_read_parse":(s[3]-s[2]).as_secs_f64(),"dispatch":start.duration_since(s[3]).as_secs_f64()}));
                });
                reference_output = boundary_json(&run.output);
                reference_output["driver_return_remainder"] = json!(reference_return.unwrap().duration_since(driver_entry).as_secs_f64()
                    - run.records.last().map_or(run.start_s, |r| r.map_end_s) - run.output_s);
                reference_output["driver_return_remainder_scope"] = json!("includes driver state teardown and dispatch/argument overhead outside the engine output clock");
                reference_ladder = run.ladder.as_ref().map_or(Value::Null, |l| json!({
                    "fraction":dft_engine::refscf::ladder::FRACTION,
                    "coarse_threshold_ry":l.coarse_threshold_ry,
                    "coarse_stopped":l.coarse_stopped,
                    "coarse_seconds":l.coarse_seconds,
                    "coarse_setup":boundary_json(&l.coarse_setup),
                    "coarse_maps":l.coarse_records.len(),
                    "fine_maps":run.records.len(),
                    "coarse_solve_calls":l.coarse_records.iter().map(|r|r.calls.len()).sum::<usize>(),
                    "fine_solve_calls":run.records.iter().map(|r|r.calls.len()).sum::<usize>(),
                    "iterations":l.coarse_records.iter().map(reference_record_json).collect::<Vec<_>>()
                }));
                let mut seconds = [0.; 8];
                for call in run.records.iter().chain(run.ladder.iter().flat_map(|l| &l.coarse_records)).flat_map(|r| &r.calls).chain(&run.verification_calls) {
                    if let Some((_, stages)) = call.inner_clocks {
                        for (sum, value) in seconds.iter_mut().zip(stages) {
                            *sum += value;
                        }
                    }
                }
                reference_stages = Some(json!({
                    "schema":"reference-davidson-stages-1",
                    "apply_h":reference_apply_h_json(),
                    "rmm":reference_rmm_json(),
                    "stage_s":dft_engine::refscf::record::INNER_STAGES.into_iter().zip(seconds).map(|(name, seconds)|(name.to_string(), json!(seconds))).collect::<serde_json::Map<_, _>>(),
                    "measured":run.full_diagnostics,
                    "note":"exclusive reference Davidson stages; includes every SCF/retry/verification call; rank max with min_over_ranks; total process wall includes output"
                }));
                reference_kept_hx = json!({
                    "cache_hits":run.kept_hx.cache_hits,
                    "full_h_starting_columns_replaced":run.kept_hx.full_h_starting_columns_replaced,
                    "delta_v_columns":run.kept_hx.delta_v_columns,
                    "retained_slab_bytes_per_rank":run.kept_hx.retained_slab_bytes_per_rank,
                    "scope":"whole run including coarse stage; slab bytes are per-rank high water, not RSS"
                });
                reference_summary = Some(json!({
                    "start_s":run.start_s,"materialize_s":run.output_s,"first_estimator_stop_iteration":run.first_reference_stop,
                    "memory_breakdown": {
                        "schema":"reference-memory-1",
                        "complete":false,
                        "units":"bytes",
                        "kind":"independent component maxima from allocated Vec capacities and initialized matrix backing; do not sum as peak RSS",
                        "scope":"fine stage; constructor and map boundaries; no allocator overhead or process high-water measurement",
                        "component_names":dft_engine::refscf::memory::NAMES,
                        "components_per_rank":run.memory_breakdown.components_per_rank,
                        "davidson_names":dft_engine::refscf::memory::DAVIDSON_NAMES,
                        "davidson_per_rank":run.memory_breakdown.davidson_per_rank,
                        "unmeasured":["atomic start temporaries", "Hamiltonian projectors and work buffers",
                            "GEMM tiles and packing", "FFT lanes and plans", "private LeanScratch and EVD backend",
                            "basis/layout metadata", "temporary XC and field transforms", "density gather scratch",
                            "verification candidate/check matrix overlap", "CLI serialization and output buffers",
                            "allocator, MPI, thread stacks and library internals"]
                    },
                    "peak_bytes_per_rank":run.peak_payload_bytes_per_rank,
                    "peak_kind":if run.ladder.is_some() { "fine-stage principal array payload estimate only; coarse stage and transfer workspace excluded; not a whole-run peak or RSS" } else { "conservative principal array payload envelope; not RSS; backend/FFT/allocator storage excluded" },
                    "peak_rss_bytes_per_rank":null,
                    "verification_calls":run.verification_calls.iter().enumerate().map(|(i,c)| {let mut v=reference_solve_json(c);if let Some(k)=run.verification_k_call_indices.get(i){v["k_index"]=json!(k);}v}).collect::<Vec<_>>()
                }));
                reference_records = run.records;
                reference_status = Some((run.reference_stopped, run.exit_reason));
                run.outcome
            })
        } else if config
            .scf
            .as_ref()
            .is_some_and(|s| s.cutoff_ladder.as_deref() == Some("off"))
        {
            // R10A retained this input pin but accidentally ignored it. The
            // strict oracle must keep the original full-cutoff cold start.
            problem.take().expect("legacy mode").run_outcome()
        } else {
            // Omitted pin: retain the historical fixed-default coarse start,
            // without restoring the retired ladder experiments or switches.
            dft_engine::cutoff_ladder::run_oracle(problem.take().expect("legacy mode"), &refs)
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
                    run_config.freeze_consulted();
                    note_switch_paths(run_config);
                    failure["scf"] = options_json(&opts);
                    selection.record(&mut failure);
                    failure["run"] = run_block(
                        run_config,
                        if matches!(e, ScfError::Diverged { .. }) {
                            "diverged"
                        } else {
                            "not_converged"
                        },
                    );
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
    launch::elapsed("scf_driver", scf_clock);
    launch::mark("scf_return");
    if let Some(derivatives) = &mut reference_derivatives {
        derivatives
            .prepare(
                dft_engine::refstart::StartInput {
                    cell: &prepared.cell,
                    atoms: &atoms,
                    potentials: &refs,
                    options: &opts,
                },
                &mut outcome,
            )
            .map_err(engine_error)?;
    }
    let r = outcome.result;
    // RUN-CONFIG: what the computation consulted (result writers below read
    // more switches for their report blocks; those reads do not count).
    run_config.freeze_consulted();
    // Everything from here to the results file is output materialization.
    let io_clock = Instant::now();
    if let Some(clock) = reference_return {
        output_seconds.insert(
            "driver_return_to_cli_output",
            io_clock.duration_since(clock).as_secs_f64(),
        );
    }
    let ngrid = r.density.len();
    let bands = r.occupations.len();
    let npw = prepared.wavefunctions.vectors().len();
    let volume = prepared.cell.volume_bohr3();
    let cell = prepared.cell.lattice_columns_bohr();
    let density_bytes: Vec<u8> = phase(reference, &mut output_seconds, "density_encode", || {
        launch::timed("density_encode", || {
            r.density.iter().flat_map(|x| x.to_le_bytes()).collect()
        })
    });
    let integral = r.density.iter().sum::<f64>() * volume / ngrid as f64;
    let density_meta = json!({"schema_version":1,"fft_grid":prepared.grid.dims(),"lattice_columns_bohr":cell,"volume_bohr3":volume,
        "units":"electrons per Bohr^3","ordering":"C: (i*n2+j)*n3+k; k fastest","endianness":"little","dtype":"f64",
        "normalization":{"expected_electrons":electrons,"integral_electrons":integral},"sha256":phase(reference, &mut output_seconds, "density_hash", || launch::timed("density_sha256", || sha256(&density_bytes)))});
    let history: Vec<_> = launch::timed("history_json", || {
        r.history
            .iter()
            .enumerate()
            .map(|(i, r)| {
                let mut row = history_json(r);
                if let Some(v) = reference_records.get(i) {
                    row["reference_scf"] = reference_record_json(v);
                }
                row
            })
            .collect()
    });
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
    let timings = closed_timings(
        &phases,
        &r.history,
        setup,
        total,
        dft_engine::phase_clock::snapshot().since(&clock_start),
    );
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
        "timestamp":timestamp(SystemTime::now().duration_since(UNIX_EPOCH).map_err(|e|error(5,e))?.as_secs()),"input_sha256":input_digest.unwrap_or_else(|| sha256(&bytes)),"species":species,
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
                "davidson_stages":phase(reference, &mut output_seconds, "report_collectives_with_pack", || over_ranks(davidson_stages_json()))?,
                "b1":b1_json(),
                "b5":{"mode":dft_numerics::switches::choice(dft_numerics::switches::Id::B5),"bordered":dft_engine::b5_counters().bordered,"declined":dft_engine::b5_counters().declined,"full":dft_engine::b5_counters().full},
                "ortho":{"mode":dft_numerics::switches::text(dft_numerics::switches::Id::ORTHO).unwrap_or("dgks:0.1 (default)"),"second_pass_skipped":dft_engine::ortho_counters().0,"second_pass_repeated":dft_engine::ortho_counters().1},
                "algebra_stream":algebra_stream_json(),
                "davidson_lean":davidson_lean_json(),
                "reuse_psi":reuse_psi_json(),
                "gemm_dispatch":gemm_dispatch_json(),
                "h_vectors_recorded":r.history.iter().map(|h|h.operator_applications).sum::<usize>(),
                "certificate_h_vectors":r.history.iter().filter_map(|h|h.fixed_point_certificate.as_ref()).map(|c|c.operator_applications).sum::<usize>(),
                "note":"recorded = per-iteration solves in history; certificate = fixed-point certificate re-solves, reported only inside each history entry's fixed_point_certificate; total work is the sum"}},"energies":energies,
        "total_energy":{"ha":r.energy.total_hartree(),"ry":from_hartree(r.energy.total_hartree(), EnergyUnit::Rydberg),"ev":from_hartree(r.energy.total_hartree(), EnergyUnit::ElectronVolt)},
        "eigenvalues":{"ha":r.eigenvalues,"ev":r.eigenvalues.iter().map(|e|from_hartree(*e, EnergyUnit::ElectronVolt)).collect::<Vec<_>>()},
        "timings":{"phases":phase_json(&phases),"setup_s":setup,"engine_setup_s":r.setup_time.as_secs_f64(),"engine_total_s":r.total_time.as_secs_f64(),"total_s":total,"fft_hpsi":"included in enclosing phases; not separately measured",
            "certificate":timings["certificate"],"atomic_start_s":timings["atomic_start_s"],"steps_s":timings["steps_s"],"closure":timings["closure"]},
        "threading":{"threads":opts.threads.max(1),"threads_requested":opts.threads,"threads_used":opts.threads.max(1),"outer_threads":opts.threads.max(1).min(r.k_points.len()),"inner_threads":opts.threads.max(1),"backend":"persistent Rayon pool shared by k points and fixed inner pieces; sequential faer and serial FFT per worker"},
        "memory_estimate":{"basis_bytes":npw*std::mem::size_of_val(&prepared.wavefunctions.vectors()[0]),"density_basis_bytes":prepared.density.vectors().len()*std::mem::size_of_val(&prepared.density.vectors()[0]),"density_grid_bytes":ngrid*8,"one_complex_grid_bytes":ngrid*16,"orbital_bytes":npw*bands*16,"peak_rss_bytes":null,"note":"array payload estimates only; scratch, projectors, allocator overhead excluded; peak RSS unavailable in-process"},
        "units":"Ha, Ry, eV, Bohr, seconds and bytes as keyed; density_residual in Bohr^-3/2; electrons/occupations in electron counts; relative errors and fractions dimensionless; see docs/cli.md"});
    selection.record(&mut result);
    if let Some((stopped, reason)) = &reference_status {
        result["convergence"]["work"]["davidson_stages"] = phase(
            reference,
            &mut output_seconds,
            "report_collectives_with_pack",
            || over_ranks(reference_stages.take().expect("reference stages")),
        )?;
        if reference_diagnostics {
            result["convergence"]["work"]["davidson_stages"]["work"] = phase(
                reference,
                &mut output_seconds,
                "report_collectives_with_pack",
                reference_davidson_work_json,
            )?;
        }
        let memory_breakdown = reference_summary
            .as_mut()
            .and_then(Value::as_object_mut)
            .and_then(|summary| summary.remove("memory_breakdown"));
        result["reference_scf"] = json!({"start":reference_start_report,"memory_breakdown":memory_breakdown,"setup_s":reference_setup,"reference_stopped":stopped,"exit_reason":reason,
            "stop_to_verified_clock":"reference-scf-wall-1",
            "cutoff_ladder":reference_ladder,
            "kept_hx":reference_kept_hx,
            "diagnostics":{"schema":"reference-scf-diagnostics-1",
                "level":if reference_diagnostics {"full"} else {"lightweight"},
                "setting":"profile.davidson_stages / BAANDI_DAVIDSON_STAGES",
                "ritz_changes":"per correction in Ha; empty when not collected",
                "inner_clocks":"null when not measured; per-call clocks are local to the output rank",
                "iteration_alias":"convergence.history[].reference_scf duplicates reference_scf.iterations[] for compatibility"},
            "mixing_input":config.scf.as_ref().and_then(|s|s.reference_mixing_input.as_ref().or(s.mixing.as_ref())),
            "mixing_implementation":"modified_broyden_difference_history",
            "verification":r.verification.as_ref().map(|v| reference_verification_json(v, &reference_records)),
            "iterations":reference_records.iter().map(reference_record_json).collect::<Vec<_>>(),
            "run":reference_summary,
            "threshold_ha":reference_profile(&config)?.threshold_ry/2.,
            "density_role":"last map output (also on exhaustion); energies refer to this output",
            "threshold_ry":reference_profile(&config)?.threshold_ry});
        if reference_records
            .iter()
            .any(|r| r.verification_cost.is_some())
        {
            let costs = || {
                reference_records
                    .iter()
                    .filter_map(|r| r.verification_cost.as_ref())
            };
            let reused = costs().filter(|c| c.reused_from.is_some()).count();
            let skipped = costs().filter(|c| c.skipped).count();
            let attempts = r.verification.as_ref().map_or(0, |v| v.attempts.len());
            let ordinary_requests = reference_records.len() - reused
                + reference_records
                    .iter()
                    .filter(|r| r.first_iteration_redo)
                    .count();
            let ordinary_maps = reference_records
                .iter()
                .map(|r| r.calls.len())
                .sum::<usize>();
            let verification_maps = result["reference_scf"]["run"]["verification_calls"]
                .as_array()
                .expect("reference verification calls")
                .len();
            result["reference_scf"]["verification_cost"] = json!({
                "gate_enabled":costs().any(|c| c.gate_enabled),
                "attempts":attempts,"skipped":skipped,"reused":reused,
                "logical_scf_iterations":reference_records.len(),
                "ordinary_map_requests":ordinary_requests,"map_requests":ordinary_requests + attempts,
                "evaluated_ordinary_maps":ordinary_maps,"verification_maps":verification_maps,
                "evaluated_maps":ordinary_maps + verification_maps,
                "solver_calls":ordinary_maps + verification_maps
            });
            if costs().any(|c| c.first_gate_enabled) {
                result["reference_scf"]["verification_cost"]["first_gate_enabled"] = json!(true);
            }
        }
    }
    result["effective_start"] = json!([
            {"stage":"single","density":if atomic_density_built {"atomic"} else {"uniform"},
             "wavefunctions":opts.start_wavefunctions.name(),
             "atomic_orbitals_used":(opts.start_wavefunctions == dft_engine::StartWavefunctions::Atomic).then_some(r.atomic_start.atomic_orbitals),
             "atomic_start_h_applications":(opts.start_wavefunctions == dft_engine::StartWavefunctions::Atomic).then_some(r.atomic_start.initial_h_applications),
             "atomic_fallback_reason":r.atomic_start.fallback_reason}]);

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
            "certificate_solves":r.verification.as_ref().map_or(0, |v| v.attempts.len())});
        println!(
            "convergence_profile: standard; certified: false; standard_threshold_ha: {}",
            p.threshold_ha
        );
    }
    result["single_path"] = single_path_json(&scf_path, &opts);
    if let Some(v) = &r.verification {
        result["verified"] = json!(v.passed());
        result["verification"] = if reference {
            reference_verification_json(v, &reference_records)
        } else {
            verification_json(v)
        };
        result["timings"]["verification_s"] = json!(v.wall().as_secs_f64());
        println!(
            "verification: fixed_point; passed: {}; attempts: {}; verification_s: {:.6}",
            v.passed(),
            v.attempts.len(),
            v.wall().as_secs_f64()
        );
    }
    if dft_engine::work_accounting::ENABLED {
        result["davidson_work"] = json!(davidson_work);
    }
    // F11: source-only GEMM accounting build (P14B_GEMM_ACCOUNTING=1): every
    // GEMM shape with its route (`real.<site>.<route>` for the real wrappers).
    if dft_numerics::linalg::accounting::ENABLED {
        result["gemm_accounting"] = json!(
            dft_numerics::linalg::accounting::take()
                .iter()
                .map(
                    |r| json!({"site":r.site,"m":r.m,"n":r.n,"k":r.k,"calls":r.calls,
                "gflop":r.flops as f64*1e-9,"seconds":r.nanoseconds as f64*1e-9})
                )
                .collect::<Vec<_>>()
        );
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
    if reference && !r.gamma_real_used {
        result["threading"]["outer_threads"] = json!(1);
        result["threading"]["backend"] = json!(
            "sequential k points on all ranks; G-distributed complex Hamiltonian and Davidson; shared inner worker pool"
        );
    }
    let replicated_reason = if points.len() > 1 && !reference {
        Some(REPLICATED_MULTI_K.to_owned())
    } else {
        r.gamma_real_fallback_reason.as_ref().map(|reason| {
            format!("real-Gamma path not taken ({reason}): every rank runs the complex solver")
        })
    };
    let (model, replicated) = parallel_model(
        crate::parallel::size(),
        if reference {
            reference_records.iter().map(|r| r.calls.len()).sum()
        } else {
            dft_engine::parallel::distributed_solves()
        },
        reference || dft_engine::parallel::m2_enabled(),
        replicated_reason,
    );
    if model == "mpi_replicated" && crate::parallel::rank() == 0 {
        eprintln!(
            "warning: MPI ranks replicate the whole calculation ({}); this run is not a distributed (MPI vs MPI) measurement",
            replicated.as_deref().unwrap_or("no distributed solve")
        );
    }
    result["execution"] = execution_json(
        &config,
        &scf_path,
        &opts,
        Solver {
            gamma_real_used: r.gamma_real_used,
            gamma_real_fallback_reason: r.gamma_real_fallback_reason.as_deref(),
            parallel_model: model,
        },
    );
    if reference {
        result["scf"]["mixing"] = json!({
            "method":"broyden", "beta":reference_profile(&config)?.beta,
            "history":reference_profile(&config)?.history
        });
        result["scf"]["mixing_metric"] = json!("hartree");
        result["execution"]["davidson_mode"] = json!("reference_generalized");
        result["execution"]["effective_solver"] = json!("reference_block_davidson");
        result["execution"]["solver_note"] = json!(
            "refstart and owned-row ref_davidson; no legacy session/ladder/solver switches; fresh-map fixed-point assessment from session.rs"
        );
        result["eigensolver_stopping_rule"] =
            json!("Ritz change: occupied ethr; empty max(5*ethr,1e-5 Ry); first call all occupied");
    }
    if crate::parallel::size() > 1 {
        let m2_on = reference || dft_engine::parallel::m2_enabled();
        let (m2_counts, m2_stages) = dft_engine::parallel::m2_stats();
        let s = dft_engine::parallel::reduced_solve_stats();
        let reduced_solves = json!({
            "solves": s.solves, "mean_order": s.mean_order, "p90_order": s.p90_order,
            "max_order": s.max_order, "sum_order_cubed": s.sum_cubed,
            "single_column": s.single_column, "bordered": s.bordered, "broadcast": s.broadcast,
        });
        let (reduction_solves, reduction_sites) = dft_engine::parallel::reduction_stats();
        let reductions = json!({
            "solves": reduction_solves,
            "buckets": dft_engine::parallel::REDUCTION_BUCKETS,
            "sites": reduction_sites
                .into_iter()
                .map(|r| json!({"site": r.site, "calls": r.calls, "payload": r.payload, "histogram": r.histogram}))
                .collect::<Vec<_>>(),
        });
        result["parallel"] = json!({
            "model": model,
            "replicated_reason": replicated,
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
            // F21: `off` = count-only (calls and payload bytes, no clock; seconds 0).
            "stage_timers": if dft_numerics::stage_clock::on() { "on" } else { "off" },
            "stages": dft_engine::parallel::stage_seconds()
                .into_iter()
                .zip(dft_engine::parallel::stage_calls())
                .map(|((name, seconds), (_, calls))| json!({"stage": name, "seconds": seconds, "calls": calls}))
                .collect::<Vec<_>>(),
            "reduced_solves": reduced_solves,
            "reductions": reductions,
            "collectives": dft_engine::parallel::collective_stats()
                .into_iter()
                .zip(dft_engine::parallel::collective_bytes())
                .map(|((op, calls, seconds), (_, bytes))| json!({"op": op, "calls": calls, "seconds": seconds, "payload_bytes": bytes}))
                .collect::<Vec<_>>(),
        });
    }
    // Only written when the calibration switch is on, so the default results.json
    // key set is unchanged.
    if dft_numerics::switches::on(dft_numerics::switches::Id::FFT_KERNEL_BY_LENGTH) {
        result["fft_kernel_choices"] = json!({
            "schema": "fft-kernel-choices-1",
            "rank": crate::parallel::rank(),
            "enabled": dft_numerics::switches::on(dft_numerics::switches::Id::FFT_KERNEL_BY_LENGTH),
            "scope": "shared and private distributed 1D lanes (start, H application, density); other ranks write rankN/results.json",
        "decision_policy": "per_run_rank_local_median",
        "rank_decisions_synchronized": false,
        "cross_run_bit_reproducibility_established": false,
            "choices": dft_numerics::fused_pair_plan::kernel_choices().iter().map(|c| json!({
                "length": c.length, "batch": c.batch,
                "direction": if c.inverse { "inverse" } else { "forward" },
                "kernel": if c.vendor { "mkl" } else { "rustfft" },
                "rust_seconds_per_batch": c.rust_seconds,
                "mkl_seconds_per_batch": c.mkl_seconds,
                "calibration_seconds": c.calibration_seconds,
            })).collect::<Vec<_>>(),
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
    if let Some((_, reason)) = &reference_status {
        writeln!(console, "Reference SCF: {reason}").unwrap();
        for record in &reference_records {
            writeln!(console, "reference iteration # {}", record.iteration).unwrap();
            for call in &record.calls {
                writeln!(
                    console,
                    "ethr = {:.6e} Ry, avg # of iterations = {:.1}",
                    call.ethr_ry, call.iterations as f64
                )
                .unwrap();
            }
            writeln!(
                console,
                "total energy = {:.12} Ry; estimated scf accuracy < {:.8e} Ry",
                2. * record.variational_energy_ha,
                record.dr2_ry
            )
            .unwrap();
        }
    }
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
    phase(
        reference,
        &mut output_seconds,
        "requested_volumetric_writes",
        || {
            launch::timed("volumetric_write", || {
                volumetric_out::write_outputs(&config, &prepared, &refs, &r, &out)
            })
        },
    )
    .map_err(|e| error(5, e))?;
    phase(reference, &mut output_seconds, "density_write", || {
        launch::timed("density_bin_write", || {
            fs::write(out.join("density.bin"), density_bytes)
        })
    })?;
    phase(
        reference,
        &mut output_seconds,
        "density_metadata_write",
        || {
            launch::timed("density_json_write", || {
                write_json(&out.join("density.json"), &density_meta)
            })
        },
    )?;
    let force = if let Some(derivatives) = &mut reference_derivatives {
        derivatives
            .evaluate_forces(
                dft_engine::refstart::StartInput {
                    cell: &prepared.cell,
                    atoms: &atoms,
                    potentials: &refs,
                    options: &opts,
                },
                &r,
            )
            .map_err(engine_error)?;
        derivatives.forces.take()
    } else {
        requested_forces(
            &config,
            &prepared.cell,
            &atoms,
            &refs,
            &opts,
            dft_engine::ElectronicState::Fixed(&r),
            density_floor,
        )?
    };
    if let Some(force) = force {
        console.push_str(&force_table(&force));
        result["forces"] = force_json(&force);
        if reference {
            result["forces"]["verification"] = json!(force.verification);
        }
        if let Some(p) = opts.standard {
            // ADR 0043: forces of a standard state rest on the passed verification
            // (the engine refuses forces otherwise). ONEPATH-1: the threshold and
            // the verification outcome are recorded next to the forces.
            let v = r.verification.as_ref();
            let last = v.and_then(|v| v.attempts.last());
            result["forces"]["preconditions"] = json!({"convergence_profile":"standard",
                "converged":outcome.converged,
                "fixed_point_verification_passed":v.is_some_and(|v| v.passed()),
                "standard_threshold_ha":p.threshold_ha,
                "verification_attempts":v.map_or(0, |v| v.attempts.len()),
                "verification_coulomb_error_ha":last.map(|a| a.coulomb_error_ha),
                "verification_residual_bound_ha":last.and_then(|a| a.residual_bound_ha),
                "state":"converged step's own state (ADR 0043 V4); forces use the same formulas as the strict path"});
        }
    }
    let stress = if let Some(derivatives) = &mut reference_derivatives {
        derivatives
            .evaluate_stress(
                dft_engine::refstart::StartInput {
                    cell: &prepared.cell,
                    atoms: &atoms,
                    potentials: &refs,
                    options: &opts,
                },
                &r,
            )
            .map_err(engine_error)?;
        derivatives.stress.take()
    } else {
        requested_stress(
            &config,
            &prepared.cell,
            &atoms,
            &refs,
            &opts,
            dft_engine::ElectronicState::Fixed(&r),
            density_floor,
        )?
    };
    if let Some(stress) = stress {
        console.push_str(&stress_table(&stress));
        result["stress"] = stress_json(&stress);
        if reference {
            result["stress"]["verification"] = json!(stress.verification);
        }
        if let Some(p) = opts.standard {
            // ONEPATH-2 (ADR 0043): stress of a standard state rests on its
            // passed verification (the engine refuses it otherwise); the cutoff
            // scans are verified standard runs (pulay_assessment.scan_*).
            let v = r.verification.as_ref();
            let last = v.and_then(|v| v.attempts.last());
            result["stress"]["preconditions"] = json!({"convergence_profile":"standard",
                "converged":outcome.converged,
                "fixed_point_verification_passed":v.is_some_and(|v| v.passed()),
                "standard_threshold_ha":p.threshold_ha,
                "verification_attempts":v.map_or(0, |v| v.attempts.len()),
                "verification_coulomb_error_ha":last.map(|a| a.coulomb_error_ha),
                "verification_residual_bound_ha":last.and_then(|a| a.residual_bound_ha),
                "state":"converged step's own state (ADR 0043 V4); stress uses the same formulas as the strict path"});
        }
    }
    // BREAKTHROUGH-0-INSTRUMENT: exclusive phase profile (BAANDI_PHASE_PROFILE=on)
    // instrumentation; absent otherwise.
    if launch::clock().is_some() {
        launch::elapsed("materialize_including_outputs", Some(io_clock));
    }
    dft_engine::profile::add_io(io_clock.elapsed().as_secs_f64());
    launch::mark("profile_finish");
    if let Some(mut mine) = dft_engine::profile::finish_rank(crate::parallel::rank()) {
        result["launchio"] = launchio_json(&mine, upf_records.take().unwrap_or_default());
        result["launchio"]["sizes"] = json!({"input_bytes":bytes.len(),
            "density_bytes":ngrid * 8,"grid_points":ngrid,"plane_waves":npw,
            "bands":bands,"history_entries":r.history.len(),"atoms":atoms.len()});
        // RUN-REPORT: this rank's own external wall closes its exclusive tree.
        mine.wall = start.elapsed().as_secs_f64();
        let gather_clock = launch::clock();
        let ranks = match dft_engine::parallel::SharedComm::active() {
            Some(comm) => dft_engine::profile::gather(&comm, &mine).map_err(|e| error(5, e))?,
            None => vec![mine],
        };
        launch::elapsed("profile_gather", gather_clock);
        let report_clock = launch::clock();
        let report = dft_engine::profile::Report::from_ranks(&ranks);
        if let Some(local) = report
            .per_rank
            .iter()
            .find(|r| r.rank == crate::parallel::rank())
        {
            launch::covered(local.covered);
        }
        let table = report.table(Some(total));
        console.push_str(&table);
        let phases = profile_json(&report, total);
        if crate::parallel::size() > 1 {
            result["parallel"]["phases"] = phases.clone();
        }
        result["phase_profile"] = phases;
        launch::elapsed("profile_report", report_clock);
    }
    if dft_engine::profile::davidson_lane::requested() {
        use dft_engine::profile::davidson_lane as lane;
        let mine = lane::take(crate::parallel::rank());
        let ranks = match dft_engine::parallel::SharedComm::active() {
            Some(comm) => lane::gather(&comm, mine).map_err(|e| error(5, e))?,
            None => vec![mine],
        };
        result["davidson_profile"] = json!({
            "schema": "davidson-profile-1", "capacity_per_rank": lane::CAPACITY,
            "fields": lane::FIELDS, "operations": lane::OPERATIONS,
            "stages": lane::STAGES, "kinds": lane::KINDS,
            "reduction_sites": lane::REDUCTION_SITES,
            "timing": "inclusive; nested rows must not be summed as exclusive wall time",
            "per_rank": ranks.iter().map(|r| json!({"rank":r.rank,
                "dropped":r.dropped,"records":r.records})).collect::<Vec<_>>()
        });
    }
    let manifest_clock = launch::clock();
    note_switch_paths(run_config);
    result["scf"]["switches"] = switches_json();
    launch::elapsed("switch_manifest", manifest_clock);
    let run_block_clock = launch::clock();
    let identity_clock = Instant::now();
    result["run"] = run_block(
        run_config,
        if outcome.converged {
            "converged"
        } else {
            "not_converged"
        },
    );
    output_seconds.insert(
        "binary_hash_and_run_manifest",
        identity_clock.elapsed().as_secs_f64(),
    );
    launch::elapsed("run_block", run_block_clock);
    if reference {
        let sum = |get: fn(&dft_engine::refscf::record::Record) -> f64| {
            reference_records.iter().map(get).sum::<f64>()
        };
        result["reference_scf"]["run"]["peak_rss_bytes_per_rank"] = json!(phase(
            reference,
            &mut output_seconds,
            "report_collectives_with_pack",
            reference_peak_rss
        )?);
        result["reference_scf"]["run"]["peak_rss_scope"] = json!(
            "Linux /proc/self/status VmHWM through output preparation, before final JSON serialization/write; null where unavailable; use external per-rank high-water measurements for the complete process"
        );
        result["reference_scf"]["wall_s"] = json!({
            "start":result["reference_scf"]["run"]["start_s"],
            "coarse_in_start":result["reference_scf"]["cutoff_ladder"]["coarse_seconds"],
            "start_scope":"inclusive run entry to first fine map; includes coarse stage and handoff when enabled",
            "solve":sum(|r|r.solve_s),"density":sum(|r|r.density_s),
            "potential":sum(|r|r.potential_s),"energy":sum(|r|r.energy_s),
            "mix":sum(|r|r.mixing_s),"density_transform":sum(|r|r.density_transform_s),
            "metric":sum(|r|r.metric_s),"energy_correction":sum(|r|r.correction_s),
            "verification":sum(|r|r.verification_s),
            "output":io_clock.elapsed().as_secs_f64()+result["reference_scf"]["run"]["materialize_s"].as_f64().unwrap_or(0.),
            "output_scope":"materialization and output preparation/writes before final results serialization; external wall includes the final write too",
            "elapsed_to_results_serialize":start.elapsed().as_secs_f64()
        });
    }
    if reference_diagnostics && crate::parallel::rank() == 0 {
        // Rank zero's raw values: deliberately never passed through over_ranks.
        // Output ends before the self-describing JSON write/console delivery.
        if let Some(row) = reference_qe_clock.iter_mut().find(|r| r.name == "output") {
            row.wall_s = result["reference_scf"]["wall_s"]["output"].as_f64();
        }
        writeln!(
            console,
            "\nQE-matched clocks (rank 0, wall seconds; CPU not collected)"
        )
        .unwrap();
        for row in &reference_qe_clock {
            if let (Some(wall), Some(calls)) = (row.wall_s, row.calls) {
                writeln!(
                    console,
                    " {:28}: {:12.6}s WALL ({:8} calls) [{}]",
                    row.name, wall, calls, row.boundary
                )
                .unwrap();
            } else {
                writeln!(console, " {:28}: unavailable (no exact boundary)", row.name).unwrap();
            }
        }
        result["reference_scf"]["qe_clock"] = json!(reference_qe_clock.iter().map(|r| json!({
            "name":r.name, "wall_s":r.wall_s, "calls":r.calls, "rank":r.rank, "boundary":r.boundary
        })).collect::<Vec<_>>());
    }
    add_xc_screening_json(&mut result);
    if reference && (mode.observe || mode.agent) {
        crate::telemetry::annotate_reference_result(&mut result, &out)?;
    }
    let result_path = out.join("results.json");
    let serialize_clock = launch::clock();
    if !reference_symmetry_record.is_null() {
        result["reference_scf"]["symmetry"] = reference_symmetry_record;
    }
    if reference {
        let total = io_clock.elapsed().as_secs_f64();
        let output = json!({"schema":"reference-output-1","rank":crate::parallel::rank(),
            "engine":reference_output,"cli":output_seconds,"cli_before_json":total + output_seconds.get("driver_return_to_cli_output").copied().unwrap_or(0.),
            "cli_other_before_json":total + output_seconds.get("driver_return_to_cli_output").copied().unwrap_or(0.) - output_seconds.values().sum::<f64>(),
            "scope":"rank-local; requested_volumetric_writes includes density/orbital formats selected by input; density.bin and metadata are mandatory outputs"});
        write_reference_results(&result_path, result, output)?;
        launch::elapsed("results_serialize_and_write", serialize_clock);
    } else {
        let data = serde_json::to_vec_pretty(&result).map_err(|e| error(5, e))?;
        launch::elapsed("results_serialize", serialize_clock);
        launch::timed("results_write", || fs::write(&result_path, data))?;
    }
    launch::output(&result_path);
    launch::mark("outputs_complete");
    Ok(RunReport {
        console,
        converged: outcome.converged,
    })
}

/// Process high water at the output boundary, never sampled inside a kernel.
fn reference_peak_rss() -> Result<Vec<Option<u64>>, RunError> {
    use dft_mpi::Collective;
    #[cfg(target_os = "linux")]
    let bytes = fs::read_to_string("/proc/self/status").ok().and_then(|s| {
        s.lines().find_map(|line| {
            let rest = line.strip_prefix("VmHWM:")?;
            rest.split_whitespace()
                .next()?
                .parse::<u64>()
                .ok()?
                .checked_mul(1024)
        })
    });
    #[cfg(target_os = "linux")]
    let local = [bytes.unwrap_or(0) as f64];
    #[cfg(not(target_os = "linux"))]
    let local = [0.];
    let values = if let Some(comm) = dft_engine::parallel::SharedComm::active() {
        let mut values = vec![0.; comm.size()];
        comm.allgatherv_f64(&local, &vec![1; comm.size()], &mut values)
            .map_err(|e| error(5, e))?;
        values
    } else {
        local.to_vec()
    };
    Ok(values
        .into_iter()
        .map(|n| (n > 0.).then_some(n as u64))
        .collect())
}

/// I/O for the additive reference spin/k/band result. Numerical work stays in
/// refscf::electronic; neither legacy spin nor finite-T drivers are constructed.
#[allow(clippy::too_many_arguments)]
fn reference_electronic_run(
    config: &Config,
    prepared: &crate::config::PreparedInput,
    opts: &ScfOptions,
    atoms: &[Atom],
    potentials: &[SpeciesPotential<'_>],
    points: &[dft_geometry::kpoints::KPoint],
    out: &Path,
    electrons: f64,
    preflight: crate::config::preflight::Report,
    species: Vec<Value>,
    input_bytes: &[u8],
    mode: &crate::telemetry::RunMode,
) -> Result<RunReport, RunError> {
    use dft_engine::refscf::electronic::Policy;
    let theta = config.finite_temperature().map_err(|e| error(1, e))?;
    let sector = config.spin_sector(electrons).map_err(|e| error(1, e))?;
    let cap = config
        .occupations
        .as_ref()
        .and_then(|o| o.band_cap)
        .unwrap_or(32);
    let policy = match (sector.clone(), theta) {
        (None, Some(theta_ha)) => Policy::FermiDirac {
            theta_ha,
            band_cap: cap,
        },
        (Some(sector), None) => Policy::FixedSpin(sector),
        (Some(sector), Some(theta_ha)) => Policy::SpinFermiDirac {
            sector,
            theta_ha,
            constraint: config.finite_constraint(),
            band_cap: cap,
        },
        _ => return Err(error(1, "reference electronic policy missing")),
    };
    if config
        .output
        .as_ref()
        .and_then(|o| o.volumetric.as_ref())
        .is_some()
    {
        return Err(error(
            1,
            "reference spin/FD writes total/channel density cubes; optional orbital/potential volumetric outputs are unavailable",
        ));
    }
    let profile = reference_profile(config)?;
    let symmetry = reference_symmetry(config)?;
    let input = dft_engine::refstart::StartInput {
        cell: &prepared.cell,
        atoms,
        potentials,
        options: opts,
    };
    let mut request = dft_engine::refscf::api::Request::new(profile, points.to_vec());
    request.symmetry = symmetry;
    request.policy = Some(policy);
    let delivery = if mode.observe || mode.agent {
        crate::telemetry::drive_reference(input, request, mode, out)?
    } else {
        dft_engine::refscf::api::solve(input, request).map_err(engine_error)?
    };
    let r = delivery.run.into_electronic().map_err(engine_error)?;
    let mut total = r.density[0].clone();
    if r.density.len() == 2 {
        total += &r.density[1];
    }
    let cube_atoms: Option<Vec<_>> = atoms
        .iter()
        .map(|a| {
            let p = potentials
                .iter()
                .find(|p| p.label == a.species)
                .unwrap()
                .potential;
            volumetric_out::atomic_number(p.element())
                .ok()
                .map(|atomic_number| volumetric_out::OutputAtom {
                    atomic_number,
                    valence_charge: p.z_valence(),
                    position_bohr: prepared
                        .cell
                        .fractional_to_cartesian_bohr(a.fractional_position),
                })
        })
        .collect();
    let cube = |name: &str, density: &dft_physics::data::Array3<f64>| -> Result<(), RunError> {
        volumetric_out::write_cube(
            std::io::BufWriter::new(fs::File::create(out.join(format!("{name}.cube")))?),
            &prepared.cell,
            cube_atoms.as_deref().unwrap_or(&[]),
            density,
            &format!("Reference {name}; electrons/Bohr^3"),
        )
        .map_err(|e| error(5, e))
    };
    cube("density-total", &total)?;
    if r.density.len() == 2 {
        cube("density-up", &r.density[0])?;
        cube("density-down", &r.density[1])?;
        cube("magnetization", &(&r.density[0] - &r.density[1]))?;
    }
    let bytes: Vec<_> = total.iter().flat_map(|v| v.to_le_bytes()).collect();
    fs::write(out.join("density.bin"), &bytes)?;
    write_json(
        &out.join("density.json"),
        &json!({"fft_grid":total.shape(),
        "lattice_columns_bohr":prepared.cell.lattice_columns_bohr(),"units":"electrons/Bohr^3",
        "ordering":"C_z_fastest","dtype":"f64","endianness":"little","sha256":sha256(&bytes)}),
    )?;
    let calls = |calls: &[(usize, usize, dft_engine::refscf::record::SolveRecord)]| -> Vec<Value> {
        calls
            .iter()
            .map(|(s, k, c)| {
                let mut v = reference_solve_json(c);
                v["spin_index"] = json!(s);
                v["k_index"] = json!(k);
                v
            })
            .collect()
    };
    let rows: Vec<_> = r.history.iter().map(|h| json!({"iteration":h.iteration,
        "dr2_ry":h.dr2_ry,"dr2_ha":h.dr2_ry/2.,"ethr_ry":h.ethr_ry,
        "first_iteration_redo":h.first_iteration_redo,"internal_energy_ha":h.energy_ha,
        "entropy_term_ha":h.entropy_term_ha,"free_energy_ha":h.free_energy_ha,
        "moment_mu_b":h.moment_mu_b,"residual_guard_ha":h.residual_guard_ha,
        "mixing":{"fraction":h.mixing_fraction,"history_used":h.mixing.occupancy,"condition":h.mixing.condition},"calls":calls(&h.calls)})).collect();
    let eigenvalues: Vec<_> = r
        .blocks
        .iter()
        .flatten()
        .flat_map(|b| b.eigenvalues_ha.iter().copied())
        .collect();
    let occupations: Vec<_> = r
        .blocks
        .iter()
        .flatten()
        .flat_map(|b| b.occupations.iter().copied())
        .collect();
    let residuals: Vec<_> = r
        .blocks
        .iter()
        .flatten()
        .flat_map(|b| b.residuals_ha.iter().copied())
        .collect();
    let chemical = r.occupations.map(|e| {
        e.roots.map(|r| {
            r.map(|r| json!({"mu_ha":r.mu_ha,"bracket_ha":r.bracket_ha,"iterations":r.iterations}))
        })
    });
    let shape = [
        r.density.len(),
        points.len(),
        r.blocks[0][0].occupations.len(),
    ];
    let constraint = if theta.is_some() {
        match config.finite_constraint() {
            dft_physics::occupations::Constraint::Total => json!({"kind":"none"}),
            dft_physics::occupations::Constraint::FixedMoment(m) => {
                json!({"kind":"fixed_total_moment","electron_difference":m})
            }
        }
    } else {
        json!({"kind":"fixed_sector","channel_electrons":sector.as_ref().map(|s| [s.n_up,s.n_down])})
    };
    let path = config.scf_path().map_err(|e| error(1, e))?;
    let (model, _) = parallel_model(crate::parallel::size(), 1, true, None);
    let mut result = json!({"schema_version":2,"output_profile":"reference_spin_fd_v1",
        "scf":options_json(opts),"preflight":preflight,"species":species,"xc":opts.functional.name(),
        "input_sha256":sha256(input_bytes),"gamma_real_used":false,"gamma_real_fallback_reason":"spin/FD uses complex Davidson",
        "tool":{"name":"dft-cli","version":env!("CARGO_PKG_VERSION"),"git_commit":env!("DFT_GIT_COMMIT")},
        "spin":{"mode":if sector.is_some(){"collinear"}else{"unpolarized"},"axis_cartesian":sector.as_ref().map(|s|s.axis_cartesian),
            "channel_order":if sector.is_some(){vec!["up","down"]}else{vec!["total"]},
            "channel_capacity_electrons":if sector.is_some(){vec![1,1]}else{vec![2]},
            "constraint":constraint,"moment_mu_b":r.moment_mu_b,"absolute_moment_mu_b":r.absolute_moment_mu_b,
            "channel_electrons":r.channel_electrons},
        "occupations":{"policy":if theta.is_some(){"fermi_dirac"}else{"fixed_insulator"},"shape":shape,
            "values":occupations,"weights_excluded":true,"theta_ha":theta,"chemical_potential":chemical,
            "chemical_potential_kind":if theta.is_none(){"not_unique_fixed_insulator"}else if matches!(config.finite_constraint(),dft_physics::occupations::Constraint::FixedMoment(_)){"per_spin"}else{"common"},
            "number_residual":r.occupations.map(|e|e.electron_residual),"moment_residual":r.occupations.and_then(|e|e.moment_residual),
            "edge_max_fraction":r.occupations.map(|e|e.edge_max_fraction),"edge_electrons":r.occupations.map(|e|e.edge_electrons)},
        "thermodynamics":{"internal_energy_ha":r.thermodynamics.internal_energy_ha,"entropy_term_ha":r.thermodynamics.entropy_term_ha,
            "entropy_over_kb":r.occupations.map_or(0.,|e|e.entropy_over_kb),"free_energy_ha":r.thermodynamics.free_energy_ha,
            "variational_quantity":if theta.is_some(){"free_energy"}else{"internal_energy"}},
        "energies":energy_json(&r.energy),"energy_audit":{"free_energy_band_ha":r.energy.band_form-r.thermodynamics.entropy_term_ha,
            "free_energy_direct_ha":r.thermodynamics.free_energy_ha},
        "band_data":{"shape":shape,"ordering":"spin_k_band","eigenvalues_ha":eigenvalues,"residuals_ha":residuals,
            "k_points":points.iter().map(|k|json!({"fractional":k.fractional,"weight":k.weight})).collect::<Vec<_>>()},
        "density":{"channel_values":r.density.iter().map(|n|n.iter().copied().collect::<Vec<_>>()).collect::<Vec<_>>()},
        "convergence":{"converged":r.converged,"iterations":r.history.len()},"verified":false,
        "single_path":single_path_json(&path,opts),
        "execution":execution_json(config,&path,opts,Solver { gamma_real_used:false,gamma_real_fallback_reason:Some("spin/FD"),parallel_model:model }),
        "parallel":{"model":model,"ranks":crate::parallel::size(),"rank":crate::parallel::rank(),"k_concurrency":1,"g_vectors":"distributed"},
        "timings":{"setup_s":r.setup_s,"engine_total_s":r.total_s},
        "reference_scf":{"start":start_report_json(r.start_report),"schema":"reference-spin-fd-1","iterations":rows,"verification_calls":calls(&r.verification_calls),
            "kept_hx":r.kept_hx.iter().enumerate().map(|(s,k)|json!({"spin_index":s,"cache_hits":k.cache_hits,"full_h_starting_columns_replaced":k.full_h_starting_columns_replaced})).collect::<Vec<_>>()},
        "limitations":["retained-basis occupation coverage is not an infinite-tail bound","polarized PZ unavailable; free-moment collinear FD with PBE routes to legacy pending G4D"]});
    if let Some(report) = &r.symmetry {
        result["reference_scf"]["symmetry"] = reference_symmetry_json(report);
    }
    if let Some(metadata) = config.symmetry_json().map_err(|e| error(1, e))? {
        result["symmetry"] = metadata;
    }
    crate::driver_select::select(
        config,
        crate::driver_select::Mode::Classic,
        crate::driver_select::Command::Run,
        &potentials.iter().map(|p| p.potential).collect::<Vec<_>>(),
    )
    .record(&mut result);
    result["execution"]["effective_solver"] = json!("reference_complex_davidson");
    let run_config = dft_numerics::switches::RunConfig::current();
    result["scf"]["switches"] = switches_json();
    result["run"] = run_block(
        run_config,
        if r.converged {
            "converged"
        } else {
            "not_converged"
        },
    );
    driver_standard_json(
        &mut result,
        Some(dft_engine::StandardProfile {
            threshold_ha: profile.threshold_ry / 2.,
        }),
        r.verification.as_ref(),
        &r.history
            .iter()
            .map(|h| Some(h.dr2_ry / 2.))
            .collect::<Vec<_>>(),
        r.converged,
    );
    add_xc_screening_json(&mut result);
    if mode.observe || mode.agent {
        crate::telemetry::annotate_reference_result(&mut result, out)?;
    }
    write_json(&out.join("results.json"), &result)?;
    Ok(RunReport {
        converged: r.converged,
        console: format!(
            "Reference spin/FD: E={:.16e} Ha; -TS={:.16e} Ha; F={:.16e} Ha; M={:.12e} mu_B; iterations={}; converged={}\n",
            r.thermodynamics.internal_energy_ha,
            -r.thermodynamics.entropy_term_ha,
            r.thermodynamics.free_energy_ha,
            r.moment_mu_b,
            r.history.len(),
            r.converged
        ),
    })
}

/// ONEPATH-3: standard-path records of the collinear and Fermi-Dirac drivers
/// (profile, threshold, verification attempts and wall) and the standard
/// preconditions of forces and stress, as on the closed-shell path.
pub(crate) fn driver_standard_json(
    result: &mut Value,
    standard: Option<dft_engine::StandardProfile>,
    verification: Option<&dft_engine::VerificationReport>,
    standard_errors: &[Option<f64>],
    converged: bool,
) {
    let Some(p) = standard else { return };
    result["convergence_profile"] = json!("standard");
    result["certified"] = json!(false);
    result["standard_threshold_ha"] = json!(p.threshold_ha);
    result["convergence"]["standard_error_ha"] = json!(standard_errors);
    result["convergence"]["estimator"] = json!(
        "E_H[d_n] + E_H[d_m] + E_0[d_m] (collinear; d_m = 0 without spin), Ha; docs/design/onepath3-spin-metals.md"
    );
    if let Some(v) = verification {
        result["verified"] = json!(v.passed());
        result["verification"] = verification_json(v);
        result["timings"]["verification_s"] = json!(v.wall().as_secs_f64());
    }
    let last = verification.and_then(|v| v.attempts.last());
    for key in ["forces", "stress"] {
        if result[key]["status"] == "available" {
            let pre = json!({"convergence_profile":"standard",
                "converged":converged,
                "fixed_point_verification_passed":verification.is_some_and(|v| v.passed()),
                "standard_threshold_ha":p.threshold_ha,
                "verification_attempts":verification.map_or(0, |v| v.attempts.len()),
                "verification_coulomb_error_ha":last.map(|a| a.coulomb_error_ha),
                "verification_residual_bound_ha":last.and_then(|a| a.residual_bound_ha),
                "state":"converged step's own state (ADR 0043 V4); same formulas as the strict path"});
            result[key]["preconditions"] = pre.clone();
            if result["state_contracts"][key]["status"] == "available" {
                result["state_contracts"][key]["preconditions"] = pre;
            }
        }
    }
}

/// Reference-only extension; the shared legacy report remains unchanged.
fn reference_verification_json(
    v: &dft_engine::VerificationReport,
    records: &[dft_engine::refscf::record::Record],
) -> Value {
    let mut result = verification_json(v);
    for (value, attempt) in result["attempts"]
        .as_array_mut()
        .unwrap()
        .iter_mut()
        .zip(&v.attempts)
    {
        value["warm_start"] = json!(
            records[attempt.iteration - 1]
                .verification_warm
                .expect("reference verification attempt")
        );
        value["h_columns"] = json!(attempt.operator_applications);
    }

    result
}

/// ADR 0043 verification record: every attempt, the pass criterion (Coulomb
/// fixed-point residual against the threshold) and information-only strict criteria.
fn verification_json(v: &dft_engine::VerificationReport) -> Value {
    let mut result = json!({"mode":"fixed_point","passed":v.passed(),"wall_s":v.wall().as_secs_f64(),
        "attempts":v.attempts.iter().map(|a| json!({"iteration":a.iteration,"passed":a.passed,
            "coulomb_error_ha":a.coulomb_error_ha,"threshold_ha":a.threshold_ha,
            "eigen_tolerance_ha":a.eigen_tolerance_ha,"residual_bound_ha":a.residual_bound_ha,
            "density_residual":a.density_residual,"occupied_eigenvalue_change_ha":a.occupied_eigenvalue_change_ha,
            "strict_criteria_met":a.strict_criteria_met,"davidson_iterations":a.davidson_iterations,
            "h_vectors":a.operator_applications,"wall_s":a.wall.as_secs_f64()})).zip(&v.attempts)
            .map(|(mut record, a)| {
                // ONEPATH-3 information keys (spin / Fermi-Dirac runs only).
                for (key, value) in [
                    ("spin_density_residual", a.spin_density_residual),
                    ("occupation_change_electrons", a.occupation_change_electrons),
                    ("chemical_potential_change_ha", a.chemical_potential_change_ha),
                    ("moment_change_mu_b", a.moment_change),
                ] {
                    if let Some(x) = value {
                        record[key] = json!(x);
                    }
                }
                record
            }).collect::<Vec<_>>(),
        "criterion":"E_H[F(rho_out) - rho_out] < threshold_ha (and the residual estimate residual_bound_ha, an empirical estimator, not a bound, when BAANDI_STANDARD_RESIDUAL_GUARD=on); strict_criteria_met is information only",
        "note":"delivered state is the converged step's own; a failed attempt continues SCF; work is also in the history entry's fixed_point_certificate"});
    if let Some(schedule) = &v.schedule {
        result["schedule"] = json!({
            "mode": schedule.mode,
            "prediction_evaluations": schedule.prediction_evaluations,
            "records": schedule.records.iter().map(|r| json!({
                "map_index": r.iteration,
                "predicted_additional_maps": r.predicted_additional_maps,
                "action": r.action,
                "would_wait": r.would_wait,
                "coulomb_error_ha": r.coulomb_error_ha,
                "ratio_to_previous_attempt": r.ratio_to_previous_attempt,
                "pass": r.passed,
            })).collect::<Vec<_>>()
        });
    }
    result
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
    let operations = if crate::driver_select::select(
        config,
        crate::driver_select::Mode::Classic,
        crate::driver_select::Command::Run,
        &potentials.iter().map(|p| p.potential).collect::<Vec<_>>(),
    )
    .reference()
    {
        config.density_operations().map_err(|e| error(1, e))?
    } else {
        config
            .symmetry_map()
            .map_err(|e| error(1, e))?
            .map_or_else(Vec::new, |m| {
                m.operations
                    .into_iter()
                    .filter(|o| !o.antiunitary)
                    .collect()
            })
    };
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
    let operations = if crate::driver_select::select(
        config,
        crate::driver_select::Mode::Classic,
        crate::driver_select::Command::Run,
        &potentials.iter().map(|p| p.potential).collect::<Vec<_>>(),
    )
    .reference()
    {
        config.density_operations().map_err(|e| error(1, e))?
    } else {
        config
            .symmetry_map()
            .map_err(|e| error(1, e))?
            .map_or_else(Vec::new, |m| {
                m.operations
                    .into_iter()
                    .filter(|o| !o.antiunitary)
                    .collect()
            })
    };
    let state = dft_engine::forces::ForceState {
        cell,
        atoms,
        potentials,
        options,
        electronic,
        density_floor,
        symmetry: &operations,
    };
    if crate::driver_select::select(
        config,
        crate::driver_select::Mode::Classic,
        crate::driver_select::Command::Run,
        &potentials.iter().map(|p| p.potential).collect::<Vec<_>>(),
    )
    .reference()
    {
        if let Some(prepared) = reference_symmetry(config)? {
            dft_engine::stress::stress_prepared(&state, &prepared)
        } else {
            dft_engine::stress::stress(&state)
        }
    } else {
        dft_engine::stress::stress(&state)
    }
    .map(Some)
    .map_err(engine_error)
}

fn stress_json(s: &dft_engine::stress::StressResult) -> Value {
    let pulay=s.pulay.as_ref().map(|p|json!({"method":p.method,"isotropic_correction_ha_bohr3":p.correction,"uncertainty_ha_bohr3":p.uncertainty,"cutoff_steps_ha":p.cutoff_steps_ha,"warning_threshold_ha_bohr3":p.warning_threshold,"warning":p.warning,"shells_resolved":p.shells_resolved,"uncertainty_kind":"two_step_spread_not_bound","anisotropic_convergence":"not_established","elapsed_seconds":p.elapsed.as_secs_f64(),"scan_profile":p.scan_profile,"scan_thresholds_ha":p.scan_thresholds_ha,"scan_energy_resolution_ha":p.scan_energy_resolution_ha}));
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

/// LANE-FFT Phase 0: same-rank records; never sum independent rank maxima.
fn fft_lane_json(report: &dft_engine::profile::Report) -> serde_json::Value {
    use dft_numerics::real::fft_profile::NAMES;
    json!({
        "schema": "lane-fft-1",
        "enabled": dft_numerics::switches::on(dft_numerics::switches::Id::FFT_PROFILE),
        "stage_timers": if dft_numerics::stage_clock::on() { "on" } else { "off" },
        "per_rank": report.fft_ranks.iter().map(|r| json!({
            "rank": r.rank,
            "buckets": r.buckets.iter().enumerate().map(|(b, bucket)| json!({
                "bucket": dft_engine::profile::BUCKET_LABELS[b],
                "maps": r.maps[b],
                "stages": NAMES.iter().zip(&bucket.0).map(|(name, v)| json!({
                    "stage": name, "seconds": v[0] as f64 * 1e-9,
                    "calls": v[1], "points": v[2], "logical_bytes": v[3],
                    "fft_flops_estimate": v[4],
                })).collect::<Vec<_>>(),
            })).collect::<Vec<_>>(),
        })).collect::<Vec<_>>(),
    })
}

/// Rank-local I/O payload sizes and existing between-map clocks (Phase 0).
fn launchio_json(mine: &dft_engine::profile::RankProfile, upf_records: Vec<Value>) -> Value {
    use dft_engine::profile::{BUCKET_LABELS, COUNT_LABELS, SEC_LABELS};
    // Export the existing local counters before the gather, rather than trying
    // to reconstruct a rank's closure from maxima belonging to different ranks.
    let buckets: Vec<_> = mine
        .buckets
        .iter()
        .zip(BUCKET_LABELS)
        .map(|(b, label)| {
            let seconds: serde_json::Map<_, _> = SEC_LABELS
                .iter()
                .zip(b.seconds)
                .filter(|(key, _)| **key == "wall" || key.starts_with("between.") || **key == "io")
                .map(|(key, value)| (key.to_string(), json!(value)))
                .collect();
            let counts: serde_json::Map<_, _> = COUNT_LABELS
                .iter()
                .zip(b.counts)
                .filter(|(key, _)| **key == "maps")
                .map(|(key, value)| (key.to_string(), json!(value)))
                .collect();
            json!({"bucket":label,"seconds":seconds,"counts":counts})
        })
        .collect();
    json!({"schema":"launchio-1","rank":mine.rank,"upf_records":upf_records,"buckets":buckets,
        "note":"Rank-local Phase 0 counters. UPF records cover run's second parse only; preflight's first parse remains inside launch.seconds.preflight_including_prepare_and_upf. Byte sizes are payloads, not allocator traffic or physical disk reads. Bucket walls are exclusive; io and between.* are subsets. Map counts normalize bucket totals, not individual map latency. Full lifecycle clocks including JSON write and MPI finalize are in launch. No new collectives."})
}

/// Exclusive phase profile (`BAANDI_PHASE_PROFILE=on`): per bucket the maximum
/// and minimum over ranks of every entry (seconds, then counts).
/// Same-host common epoch, expressed relative to the earliest rank anchor.
fn arrival_json(report: &dft_engine::profile::Report) -> serde_json::Value {
    use dft_numerics::stage_clock::arrival::CAPACITY;
    if report.arrival_ranks.is_empty() {
        return serde_json::Value::Null;
    }
    let epoch = report
        .arrival_ranks
        .iter()
        .map(|(_, t)| t.epoch_ns)
        .min()
        .unwrap_or(0);
    let timed = report.arrival_ranks.first().is_some_and(|(_, t)| t.timed);
    let solver_names = dft_engine::parallel::stage_calls();
    let reduction_names = dft_engine::parallel::reduction_stats().1;
    let stage_name = |id: u64| -> String {
        match id {
            0 => "unattributed".into(),
            1..=16 => format!("solver.{}", solver_names[id as usize - 1].0),
            101..=118 => format!(
                "fft.{}",
                dft_numerics::real::fft_profile::NAMES[id as usize - 101]
            ),
            201..=208 => format!(
                "m2.operator.{}",
                [
                    "mirror_exchange",
                    "pack_unpack_kinetic",
                    "fft",
                    "nonlocal",
                    "density",
                    "nonlocal_energy",
                    "gather_planes",
                    "band_values"
                ][id as usize - 201]
            ),
            301 => "nonlocal.projector_gemm".into(),
            401..=411 => format!("solver.reduce.{}", reduction_names[id as usize - 401].site),
            _ => format!("stage_slot.{}", id - 1),
        }
    };
    json!({
        "schema": "arrival-skew-2", "enabled": !report.arrival_ranks.is_empty(),
        "mode": if timed { "timed" } else { "counts" }, "rank_count": report.ranks,
        "capacity_per_rank": CAPACITY,
        "complete": report.arrival_ranks.len() == report.ranks && report.arrival_ranks.iter().all(|(_, t)| !t.overflow),
        "clock_domain": "same-host SystemTime anchored Instant",
        "clock_valid": timed && report.arrival_ranks.iter().all(|(_, t)| t.drift_ns <= 1_000_000),
        "clock_error_ns": report.arrival_ranks.iter().map(|(_, t)| t.error_ns + t.drift_ns).max().unwrap_or(0),
        "ranks": report.arrival_ranks.iter().map(|(rank, t)| {
            let offset = t.epoch_ns.saturating_sub(epoch);
            json!({"rank": rank, "clock_error_ns": t.error_ns, "clock_drift_ns": t.drift_ns,
                "events": t.events.iter().map(|e| json!({
                    "map_id": e[0], "bucket": dft_engine::profile::BUCKET_LABELS[e[1] as usize],
                    "seq": e[2], "op": (["allgatherv", "tree_sum", "bcast", "alltoall", "barrier"][e[3] as usize]),
                    "site": if e[17] == 301 { "nonlocal.reduce".to_string() } else { stage_name(e[17]) },
                    "payload_bytes": e[4], "enter_ns": timed.then_some(offset + e[5]),
                    "exit_ns": timed.then_some(offset + e[6]), "previous_stage": stage_name(e[7]),
                    "previous_stage_ns": timed.then_some(e[8]), "previous_stage_points": e[9],
                    "previous_stage_logical_bytes": e[10],
                    "previous_stage_end_ns": timed.then_some(offset + e[11]),
                    "previous_exit_ns": timed.then_some(offset + e[12]), "call_index": e[13],
                    "previous_stage_shape": [e[14], e[15], e[16]],
                    "stage_duration_is_exclusive": false
                })).collect::<Vec<_>>()})
        }).collect::<Vec<_>>()
    })
}

pub fn profile_json(report: &dft_engine::profile::Report, run_wall: f64) -> serde_json::Value {
    let covered = report.covered_wall();
    let closure = report.closure();
    json!({
        "schema": "phase-profile-2",
        "fft_lane": fft_lane_json(report),
        "arrival_trace": arrival_json(report),
        "per_rank": report.per_rank.iter().map(|r| json!({
            "rank": r.rank, "wall_s": r.wall, "covered_s": r.covered, "residual_s": r.residual,
        })).collect::<Vec<_>>(),
        "closure": closure.map(|c| json!({
            "covered_max_s": c.covered_max, "covered_min_s": c.covered_min,
            "critical_rank": c.critical_rank, "imbalance": c.imbalance,
            "wall_max_s": c.wall_max, "wall_min_s": c.wall_min,
            "residual_max_s": c.residual_max, "residual_min_s": c.residual_min,
            "note": "each rank's exclusive tree is closed on that rank (sum of its own bucket walls against its own external wall); max/min over ranks of the closed figures; the external wall stays the verdict",
        })),
        "ranks": report.ranks,
        "stage_timers": if dft_numerics::stage_clock::on() { "on" } else { "off" },
        "run_wall_s": run_wall,
        "covered_wall_s": covered,
        "covered_wall_note": "legacy: sum over buckets of the per-bucket maximum over ranks (maxima may belong to different ranks; not a closed tree, see closure)",
        "remainder_s": run_wall - covered,
        "buckets": report.entries().into_iter().map(|(label, entries)| {
            json!({
                "bucket": label,
                "max": entries.iter().map(|(k, mx, _)| (k.to_string(), json!(mx))).collect::<serde_json::Map<_, _>>(),
                "min": entries.iter().map(|(k, _, mn)| (k.to_string(), json!(mn))).collect::<serde_json::Map<_, _>>(),
            })
        }).collect::<Vec<_>>(),
        "note": "top-level buckets are mutually exclusive (init, coarse, fine, fine_redo, verify_pass, verify_fail, materialize); wall = potential + solve + density + energy (+ mixing, atomic start, io); solve.* are exclusive stages inside solve (solve.other = the remainder); coll.* and m2.* are inside the stages, not exclusive; seconds = max over ranks (critical path) and min over ranks (imbalance), never summed over ranks; counts are per rank (max); with BAANDI_STAGE_TIMERS=off only the counts and the map phase clocks are filled"
    })
}

/// Reference-only full diagnostics. Gather serialized rank records with variable
/// lengths: zero-row ranks and root-only EVD need not have identical shapes.
/// Reporting communication is outside the solve and never charged as solver work.
fn reference_davidson_work_json() -> Result<serde_json::Value, RunError> {
    use dft_physics::distributed_real_gamma::ApplyH;
    let names = ApplyH::DAVIDSON_WORK_STAGES;
    let solves = ApplyH::take_davidson_work();
    let mut shapes = std::collections::BTreeMap::new();
    let mut stages = serde_json::Map::new();
    for (i, name) in names.iter().enumerate() {
        let rows: Vec<_> = solves.iter().map(|s| &s.stages[i]).collect();
        let seconds: f64 = rows.iter().map(|s| s.seconds).sum();
        let flops: u64 = rows.iter().map(|s| s.gemm_flops + s.scalar_flops).sum();
        let streamed: u64 = rows.iter().map(|s| s.streamed_bytes).sum();
        let operands: u64 = rows.iter().map(|s| s.gemm_operand_bytes).sum();
        let kinds = ["barrier", "tree_sum", "broadcast", "allgatherv", "alltoall"];
        let communication: serde_json::Map<_, _> = kinds
            .iter()
            .enumerate()
            .map(|(k, kind)| {
                (
                    kind.to_string(),
                    json!({
                        "calls":rows.iter().map(|s| s.collective_calls[k]).sum::<u64>(),
                        "payload_bytes":rows.iter().map(|s| s.collective_bytes[k]).sum::<u64>(),
                    }),
                )
            })
            .collect();
        stages.insert(name.to_string(), json!({
            "seconds":seconds,
            "calls":rows.iter().map(|s| s.calls).sum::<u64>(),
            "nominal_flops":flops,
            "gemm_calls":rows.iter().map(|s| s.gemm_calls).sum::<u64>(),
            "gemm_flops":rows.iter().map(|s| s.gemm_flops).sum::<u64>(),
            "gemm_operand_bytes":operands,
            "scalar_flops":rows.iter().map(|s| s.scalar_flops).sum::<u64>(),
            "hypot_calls":rows.iter().map(|s| s.hypot_calls).sum::<u64>(),
            "sqrt_calls":rows.iter().map(|s| s.sqrt_calls).sum::<u64>(),
            "copied_bytes":rows.iter().map(|s| s.copied_bytes).sum::<u64>(),
            "zeroed_bytes":rows.iter().map(|s| s.zeroed_bytes).sum::<u64>(),
            "streamed_bytes":streamed,
            "collectives":communication,
            "collective_seconds":rows.iter().map(|s| s.collective_seconds).sum::<f64>(),
            "rejected_factor_panels":rows.iter().map(|s| s.rejected_factor_panels).sum::<u64>(),
            "nominal_gflop_s":if seconds > 0. { Some(flops as f64 / seconds / 1e9) } else { None },
            "logical_gb_s":if seconds > 0. { Some((operands + streamed) as f64 / seconds / 1e9) } else { None },
        }));
    }
    for solve in &solves {
        for g in &solve.gemms {
            *shapes
                .entry((g.stage, g.site, g.m, g.n, g.k))
                .or_insert(0_u64) += 1;
        }
    }
    let local = json!({
        "solves":solves.len(), "stages":stages,
        "gemm_shapes":shapes.iter().map(|(&(stage,site,m,n,k), &calls)| json!({
            "stage":names[stage], "site":site, "m":m, "n":n, "k":k, "calls":calls,
        })).collect::<Vec<_>>()
    });
    let comm = dft_engine::parallel::SharedComm::active();
    reference_davidson_work_ranks(comm.as_ref().map(|c| c as &dyn dft_mpi::Collective), local)
}

fn reference_davidson_work_ranks(
    comm: Option<&dyn dft_mpi::Collective>,
    local: Value,
) -> Result<Value, RunError> {
    let per_rank = if let Some(comm) = comm {
        let bytes = serde_json::to_vec(&local).map_err(|e| error(5, e))?;
        let mut sizes = vec![0.; comm.size()];
        comm.allgatherv_f64(&[bytes.len() as f64], &vec![1; comm.size()], &mut sizes)
            .map_err(|e| error(5, e))?;
        let counts: Vec<_> = sizes.iter().map(|&v| v as usize).collect();
        let mut all = vec![0.; counts.iter().sum()];
        let payload: Vec<_> = bytes.iter().map(|&v| f64::from(v)).collect();
        comm.allgatherv_f64(&payload, &counts, &mut all)
            .map_err(|e| error(5, e))?;
        let mut at = 0;
        let mut ranks = Vec::with_capacity(comm.size());
        for count in counts {
            let bytes: Vec<_> = all[at..at + count].iter().map(|&v| v as u8).collect();
            ranks.push(
                serde_json::from_slice::<serde_json::Value>(&bytes).map_err(|e| error(5, e))?,
            );
            at += count;
        }
        ranks
    } else {
        vec![local]
    };
    fn maxima(out: &mut serde_json::Value, next: &serde_json::Value) {
        if let (Some(a), Some(b)) = (out.as_object_mut(), next.as_object()) {
            for (key, value) in a {
                maxima(value, &b[key]);
            }
        } else if let (Some(a), Some(b)) = (out.as_u64(), next.as_u64()) {
            if b > a {
                *out = next.clone();
            }
        } else if next
            .as_f64()
            .is_some_and(|b| out.as_f64().is_none_or(|a| b > a))
        {
            *out = next.clone();
        }
    }
    let mut max = per_rank[0]["stages"].clone();
    for rank in &per_rank[1..] {
        maxima(&mut max, &rank["stages"]);
    }
    Ok(json!({
        "schema":"reference-davidson-work-1",
        "rank0":per_rank[0], "max":max, "per_rank":per_rank,
        "note":"Full diagnostics only; all completed reference solves including verification. Exclusive rank-local stage seconds; collective_seconds is a subset. Max is componentwise over rank-local totals/rates, not a summed critical path. GEMM shapes are logical adapter calls, not BLAS microtiles; H and EVD internals are opaque (their zero FLOPs/bytes are uncounted, not zero work). copied_bytes counts destination bytes; streamed_bytes includes read+write. Logical GB/s = (GEMM operand bytes + streamed bytes)/stage seconds, not hardware bandwidth. Scalar factor work is a nominal full-panel model; rejected panels flagged. Backend packing/scans/allocations and diagnostic storage are excluded."
    }))
}

#[cfg(test)]
mod davidson_work_json_tests {
    use super::*;
    use dft_mpi::Collective;

    #[test]
    fn variable_rank_shapes_and_integer_counts_survive_work_reporting() {
        let outputs = dft_mpi::simulate(3, |c| {
            let rank = c.rank();
            let local = json!({"solves":2,"stages":{"border_s":{
                "seconds":rank as f64 + 0.5,
                "nominal_flops":9_007_199_254_740_992_u64 + rank as u64,
                "logical_gb_s":if rank == 0 {None} else {Some(4./rank as f64)},
                "collectives":{"tree_sum":{"calls":1,"payload_bytes":8*(rank+1)}}}},
                "gemm_shapes":(0..rank).map(|i|json!({"m":i,"n":3,"k":rank})).collect::<Vec<_>>()});
            reference_davidson_work_ranks(Some(c), local).unwrap()
        });
        assert!(outputs.windows(2).all(|v| v[0] == v[1]));
        let w = &outputs[0];
        assert_eq!(w["rank0"]["stages"]["border_s"]["seconds"], 0.5);
        assert_eq!(w["max"]["border_s"]["seconds"], 2.5);
        assert_eq!(
            w["max"]["border_s"]["nominal_flops"],
            9_007_199_254_740_994_u64
        );
        assert_eq!(w["max"]["border_s"]["logical_gb_s"], 4.);
        assert_eq!(
            w["max"]["border_s"]["collectives"]["tree_sum"]["payload_bytes"],
            24
        );
        for rank in 0..3 {
            assert_eq!(
                w["per_rank"][rank]["gemm_shapes"].as_array().unwrap().len(),
                rank
            );
        }
        let serial = reference_davidson_work_ranks(None, w["rank0"].clone()).unwrap();
        assert_eq!(serial["rank0"], w["rank0"]);
        assert_eq!(serial["max"], w["rank0"]["stages"]);
    }
}

/// Real-Davidson stage timers (`BAANDI_DAVIDSON_STAGES=on`; zeros otherwise),
/// cumulative over the process, seconds.
fn davidson_stages_json() -> serde_json::Value {
    let s = dft_engine::stage_seconds();
    let mut out = stage_block_json(&s);
    out["note"] = json!(
        "BAANDI_DAVIDSON_STAGES=on; serial real Davidson only; width 0 = returning iteration; h_apply includes certificate and restart applications"
    );
    // P0.2 work timers (same switch): FFT by kind, h_apply parts, per run stage.
    let work = dft_numerics::work_timers::snapshot();
    let pooled = work
        .iter()
        .fold(dft_numerics::work_timers::WorkTimes::default(), |a, b| {
            a.plus(b)
        });
    out["fft"] = fft_json(&pooled);
    out["h_apply_parts"] = h_apply_parts_json(&pooled);
    out["collectives"] = collectives_json(&pooled);
    let stages = dft_engine::stage_seconds_by_run_stage();
    out["by_run_stage"] = dft_numerics::work_timers::RUN_STAGE_LABELS
        .iter()
        .zip(stages.iter().zip(&work))
        .map(|(label, (s, w))| {
            let mut block = stage_block_json(s);
            block["fft"] = fft_json(w);
            block["h_apply_parts"] = h_apply_parts_json(w);
            block["collectives"] = collectives_json(w);
            (label.to_string(), block)
        })
        .collect::<serde_json::Map<_, _>>()
        .into();
    out["work_timers_note"] = json!(
        "BAANDI_DAVIDSON_STAGES=on. fft and h_apply_parts: seconds summed over worker threads (thread-seconds; at T>1 compare with wall only after dividing by the busy threads), counts independent of the thread count. fft.pair_*: paired orbital transforms of the real Gamma h_apply (one inverse and one forward transform per pair, as QE's fftw count of wave transforms); pair_potential_s is the grid multiply between them, not FFT. fft.density_*: density transforms (one inverse each). fft.other_*: other wave transforms (complex-path local applies, folded odd tail), one inverse and one forward each. fft.kept_*: REUSE-PSI forward transforms of kept psi(r) times the potential change (one forward each, no inverse). h_apply_parts: real Gamma block applications (apply_s whole call, local_s wall of the local part on the calling thread including pack/unpack, nonlocal projection P^T X and back product); remainder_s = apply - local - nonlocal (kinetic, scans). by_run_stage splits every block by run stage (fine, ladder coarse, verification = certificate block); the pooled keys are the sums. MPI M2 path (distributed solver and operator): stage_s maps the distributed solver stages (reduced_gram -> reduced, reduced_eigen -> eigen, rotate, residual, orthonormalize, apply_total -> h_apply incl. the per-block finiteness scan outside it, precondition; sketch 0; iteration widths are not recorded); h_apply_parts: local_s = mirror exchange, pack, distributed pair FFTs with transposes and potential, mirror-h exchange, unpack and kinetic term; nonlocal_projection_s = local P^T X and its tree sum; nonlocal_back_s = coupling, P W and image update; fft.pair_unsplit_s = the distributed pair FFT call (transposes included), fft.density_* the band-density pair transforms (one inverse per pair), fft.kept_* the potential-change transforms (forward per pair; inverse of non-kept pairs inside the seconds); collectives.* = calls and seconds of every collective of the run stage (all) and of those inside H applications (in_h_apply). Seconds need BAANDI_STAGE_TIMERS=on (off: counts only). With two or more ranks every number is the MAXIMUM over ranks (critical path, entries derived per rank first), min_over_ranks the minimum (imbalance), never a sum; simulated ranks (threads of one process) share the process counters."
    );
    out
}
fn collectives_json(w: &dft_numerics::work_timers::WorkTimes) -> serde_json::Value {
    let block = |calls: &[u64], secs: &[f64]| {
        dft_numerics::work_timers::COLLECTIVE_LABELS
            .iter()
            .enumerate()
            .map(|(i, l)| (l.to_string(), json!({"calls":calls[i],"s":secs[i]})))
            .collect::<serde_json::Map<_, _>>()
    };
    json!({"all":block(&w.coll_calls, &w.coll),"in_h_apply":block(&w.apply_coll_calls, &w.apply_coll)})
}

/// With two or more ranks (registered communicator): every number of `value`
/// replaced by its maximum over ranks, plus `min_over_ranks` (the same tree
/// with the minima) and `ranks`. Collective: every rank must call with a tree
/// of the same shape (one allgatherv of the numbers in traversal order).
fn over_ranks(value: serde_json::Value) -> Result<serde_json::Value, RunError> {
    // Without the work timers every entry is 0: no collective (the switch is
    // replicated, so every rank takes the same branch).
    if !dft_numerics::work_timers::enabled() {
        return Ok(value);
    }
    match dft_engine::parallel::SharedComm::active() {
        Some(comm) => over_ranks_with(&comm, value),
        None => Ok(value),
    }
}
fn over_ranks_with<C: dft_mpi::Collective + ?Sized>(
    comm: &C,
    value: serde_json::Value,
) -> Result<serde_json::Value, RunError> {
    fn flatten(v: &serde_json::Value, out: &mut Vec<f64>) {
        match v {
            serde_json::Value::Number(n) => out.push(n.as_f64().unwrap_or(0.0)),
            serde_json::Value::Array(a) => a.iter().for_each(|x| flatten(x, out)),
            serde_json::Value::Object(m) => m.values().for_each(|x| flatten(x, out)),
            _ => {}
        }
    }
    fn refill(v: &mut serde_json::Value, vals: &[f64], at: &mut usize) {
        match v {
            serde_json::Value::Number(n) => {
                let x = vals[*at];
                *at += 1;
                *v = if n.is_u64() {
                    json!(x as u64)
                } else if n.is_i64() {
                    json!(x as i64)
                } else {
                    json!(x)
                };
            }
            serde_json::Value::Array(a) => a.iter_mut().for_each(|x| refill(x, vals, at)),
            serde_json::Value::Object(m) => m.values_mut().for_each(|x| refill(x, vals, at)),
            _ => {}
        }
    }
    let mut flat = Vec::new();
    flatten(&value, &mut flat);
    let (n, ranks) = (flat.len(), comm.size());
    let mut all = vec![0.0; n * ranks];
    comm.allgatherv_f64(&flat, &vec![n; ranks], &mut all)
        .map_err(|e| error(5, e))?;
    let reduce = |f: fn(f64, f64) -> f64| -> Vec<f64> {
        (0..n)
            .map(|i| (1..ranks).fold(all[i], |a, r| f(a, all[r * n + i])))
            .collect()
    };
    let (max, min) = (reduce(f64::max), reduce(f64::min));
    let mut hi = value.clone();
    refill(&mut hi, &max, &mut 0);
    let mut lo = value;
    refill(&mut lo, &min, &mut 0);
    hi["min_over_ranks"] = lo;
    hi["ranks"] = json!(ranks);
    Ok(hi)
}
fn stage_block_json(s: &dft_engine::StageSecondsReport) -> serde_json::Value {
    let labels = dft_engine::DAVIDSON_WIDTH_LABELS;
    json!({
        "timed_solves":s.timed_solves,
        "iteration_s_by_correction_width":labels.iter().zip(s.width_seconds).map(|(l,v)|(l.to_string(),json!(v))).collect::<serde_json::Map<_,_>>(),
        "iterations_by_correction_width":labels.iter().zip(s.width_iterations).map(|(l,v)|(l.to_string(),json!(v))).collect::<serde_json::Map<_,_>>(),
        "stage_s":{"reduced":s.reduced,"eigen":s.eigen,"rotate":s.rotate,"residual":s.residual,"orthonormalize":s.orthonormalize,"h_apply":s.h_apply,"sketch":s.sketch,"precondition":s.precondition},
    })
}
fn fft_json(w: &dft_numerics::work_timers::WorkTimes) -> serde_json::Value {
    let pair_fft = w.pair_inverse + w.pair_forward + w.pair_unsplit;
    json!({
        "pair_inverse_s":w.pair_inverse,"pair_forward_s":w.pair_forward,"pair_unsplit_s":w.pair_unsplit,
        "pair_potential_s":w.pair_potential,"pair_transforms":w.pair_transforms,"pair_columns":w.pair_columns,
        "density_s":w.density,"density_transforms":w.density_transforms,"density_columns":w.density_columns,
        "other_s":w.other,"other_transforms":w.other_transforms,"other_columns":w.other_columns,
        "kept_forward_s":w.kept_forward,"kept_transforms":w.kept_transforms,"kept_columns":w.kept_columns,
        "wave_fft_s":pair_fft + w.density + w.other + w.kept_forward,
        "wave_fft_calls":2 * w.pair_transforms + w.density_transforms + 2 * w.other_transforms + w.kept_transforms,
    })
}
fn h_apply_parts_json(w: &dft_numerics::work_timers::WorkTimes) -> serde_json::Value {
    json!({
        "apply_s":w.apply,"apply_calls":w.apply_calls,"apply_columns":w.apply_columns,"local_s":w.local,
        "nonlocal_projection_s":w.nonlocal_projection,"nonlocal_back_s":w.nonlocal_back,
        "remainder_s":w.apply - w.local - w.nonlocal_projection - w.nonlocal_back,
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

/// Lean Davidson (`BAANDI_DAVIDSON=lean`, work order phase B) counters and the
/// serial real Davidson GEMM flops by site (both paths; 2 m n k per product).
fn davidson_lean_json() -> serde_json::Value {
    let c = dft_engine::lean_counters();
    let f = dft_engine::flop_counters();
    json!({
        "mode":dft_numerics::switches::choice(dft_numerics::switches::Id::DAVIDSON),
        "solves":c.solves,"iterations":c.iterations,"residual_columns":c.residual_columns,
        "appends":c.appends,"appended_columns":c.appended_columns,
        "reprojections":c.reprojections,"reprojection_rate":c.reprojection_rate(),
        "fallbacks":c.fallbacks,"fallback_rate":c.fallback_rate(),"fallback_dropped_columns":c.fallback_dropped_columns,
        "rebuilds":c.rebuilds,"rebuild_columns":c.rebuild_columns,"contractions":c.contractions,
        "returns":c.returns,"pivot_min":c.pivot_min,"pivot_below_1e-1":c.pivot_below_1e_1,
        "pivot_below_1e-2":c.pivot_below_1e_2,
        "gemm_gflop":{"reduced":f.reduced as f64*1e-9,"overlap":f.overlap as f64*1e-9,
            "ritz_loop":f.ritz_loop as f64*1e-9,"return":f.ret as f64*1e-9,
            "contraction":f.contraction as f64*1e-9,"projection":f.projection as f64*1e-9,
            "gram_schmidt":f.gram_schmidt as f64*1e-9,"fallback":f.fallback as f64*1e-9,
            "small":f.small as f64*1e-9,"total":f.total() as f64*1e-9},
        "eigen_full_solves":f.eigen_full,"eigen_full_m3":f.eigen_full_m3,
        "proj_cache":{"mode":dft_numerics::switches::choice(dft_numerics::switches::Id::LEAN_PROJ_CACHE),
            "returns":c.proj_cache_returns,"misses":c.proj_cache_misses,
            "kept_loads":c.proj_cache_kept_loads,"columns_copied":c.proj_cache_columns,
            "audits_from_cache":dft_engine::nl_audit_counters().0,
            "audits_direct":dft_engine::nl_audit_counters().1},
        "restart_keep":{"mode":dft_numerics::switches::choice(dft_numerics::switches::Id::RESTART_KEEP),
            "keep_contractions":c.keep_contractions,"kept_directions":c.kept_directions,
            "dropped":c.keep_dropped,"fallbacks":c.keep_fallbacks},
        "note":"process totals of the serial real Davidson (all solves, both paths; MPI distributed solver not counted)"
    })
}

/// F11: effective route, shape and fallback reason of every real GEMM wrapper
/// call of this process (count-only; `dft_numerics::gemm_dispatch`).
fn gemm_dispatch_json() -> serde_json::Value {
    let rows = dft_numerics::gemm_dispatch::rows();
    let (mut mkl, mut total) = (0_u64, 0_u64);
    for r in &rows {
        total += r.flops;
        if r.is_mkl() {
            mkl += r.flops;
        }
    }
    json!({
        "backend":if dft_numerics::MKL_BACKEND {"mkl"} else {"pure"},
        "rows":rows.iter().map(|r|json!({"site":r.site,"route":r.route,"calls":r.calls,
            "gflop":r.flops as f64*1e-9,"max_m":r.max_m,"max_n":r.max_n,"max_k":r.max_k})).collect::<Vec<_>>(),
        "rows_as_if_mkl":dft_numerics::gemm_dispatch::rows_as_if_mkl().iter().map(|r|json!({"site":r.site,"route":r.route,"calls":r.calls,
            "gflop":r.flops as f64*1e-9,"max_m":r.max_m,"max_n":r.max_n,"max_k":r.max_k})).collect::<Vec<_>>(),
        "wide_requests":dft_numerics::gemm_dispatch::wide_requests().into_iter()
            .map(|(s,c)|json!({"site":s,"calls":c})).collect::<Vec<_>>(),
        "mkl_flop_fraction":if total > 0 {mkl as f64/total as f64} else {0.0},
        "note":"process totals of the real GEMM wrappers (serial and MPI real paths); output m x n, reduction k, gflop = 2 m n k; route mkl_wide = dgemm only because of the wide entry; pure_output_bound / pure_panel_bound / pure_long / pure_short / pure_layout = fallback reason on an MKL build; pure_direct = direct pure-kernel entry; narrow = streaming kernels; rows_as_if_mkl (pure builds only) = the same calls under the route an MKL build with dgemm on would take (same rule and workspace); complex GEMMs (atomic start, k points) are not included"
    })
}

/// REUSE-PSI (`BAANDI_REUSE_PSI`, work order C2) counters of this process.
fn reuse_psi_json() -> serde_json::Value {
    let c = dft_engine::reuse_psi_counters();
    let l = dft_engine::lean_counters();
    json!({
        "mode":dft_numerics::switches::choice(dft_numerics::switches::Id::REUSE_PSI),
        "armed_starts":c.armed,"kept_starts":l.kept_starts,"declined_starts":l.kept_declined,
        "kept_returns":l.kept_stored,"not_kept":c.not_kept,
        "unchanged_potential":c.unchanged_potential,
        "kept_pairs":c.kept_pairs,"full_pairs":c.full_pairs,
        "inverse_transforms_saved":c.transforms_saved(),
        "h_columns_not_applied":l.kept_columns,
        "nonlocal_columns_saved":c.columns_saved,
        "bytes_kept_max":c.bytes_kept,"memory_gate_bytes":dft_engine::PSI_KEEP_MAX_BYTES,"gated_states":c.gated,
        "note":"armed = solves whose kept H X was brought to the new potential; kept_starts = armed starts the lean solver took (its initial equals the kept X), declined = armed but initial differed; kept_pairs: one forward transform of kept psi(r) (inverse saved), full_pairs: pair transformed both ways; nonlocal_columns_saved: columns whose P^T X, back product and kinetic pass were skipped (counted when armed); h_columns_not_applied: columns of the initial H application the solver skipped"
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

/// Full-diagnostic RMM records, reduced across ranks by the enclosing writer.
fn reference_rmm_json() -> Value {
    use dft_engine::refscf::driver::{RMM_STAGES as STAGES, take_rmm_records as take_records};
    let records = take_records();
    let solves: Vec<_> = records.iter().map(|r| json!({
        "scf_iteration":r.scf_iteration, "solver":"rmm", "total_s":r.total_seconds,
        "stage_s":STAGES.iter().zip(r.seconds).map(|(k,v)| (k.to_string(),json!(v))).collect::<serde_json::Map<_,_>>(),
        "stage_gemm_flops":STAGES.iter().zip(r.gemm_flops).map(|(k,v)| (k.to_string(),json!(v))).collect::<serde_json::Map<_,_>>(),
        "rotation_h_columns":r.rotation_columns, "sweep_h_columns":r.step_columns,
        "h_columns":r.rotation_columns+r.step_columns.iter().sum::<usize>(),
        "gemm_flops":r.gemm_flops.iter().sum::<u64>(),
        "newest_fallbacks":r.newest_fallbacks, "clamped_low":r.clamped_low,
        "clamped_high":r.clamped_high, "gs_fallbacks":r.gs_fallbacks,
        "minimum_gs_norm":r.min_gs_norm, "unconverged_bands":r.unconverged,
        "communication":{"reductions":r.communication.reductions,"reduction_bytes":r.communication.reduction_bytes,
            "broadcasts":r.communication.broadcasts,"gathers":r.communication.gathers,"alltoalls":r.communication.alltoalls}
    })).collect();
    json!({"schema":"reference-rmm-stages-1", "solves":solves,
        "note":"One bounded sweep per listed SCF iteration; verification always Davidson. Rank-local exclusive seconds; GEMM FLOPs exclude scalar band work and EVD internals. Off/standard diagnostics retain no RMM records. Small DIIS solves are replicated. Legacy inner clocks do not describe RMM stages; use this table."})
}

/// Per-solve operator records are consumed once after the reference driver.
/// The surrounding over_ranks merge computes extrema of rank-local totals.
fn reference_apply_h_json() -> Value {
    use dft_physics::distributed_real_gamma::{APPLY_H_STAGES, ApplyH, take_apply_h_solves};
    let solves = take_apply_h_solves();
    let block = |s: &ApplyH| {
        let f = |indices: &[usize]| {
            indices
                .iter()
                .map(|&i| s.fft.0[i][0] as f64 * 1e-9)
                .sum::<f64>()
        };
        // The original solver counter records two dispatching entries per H.
        // P3-NL replaces only the projection entry; keep route labels honest.
        let single = if dft_numerics::switches::on(dft_numerics::switches::Id::NONLOCAL_SINGLE_GEMM)
        {
            s.wide_calls / 2
        } else {
            0
        };
        let mut row = json!({
            "stage_s": APPLY_H_STAGES.iter().zip(s.seconds).map(|(k,v)| (k.to_string(),json!(v))).collect::<serde_json::Map<_,_>>(),
            "fft": {
                "measured": s.fft.0.iter().any(|v| v[0] > 0),
                "inverse_s": f(&[0,1,2,3,4,5,6,7]),
                "forward_s": f(&[9,10,11,12,13,14,15,16]),
                "inverse_local_1d_s": f(&[1]), "inverse_local_2d_s": f(&[5,7]),
                "inverse_alltoall_s": f(&[3]), "inverse_transpose_pack_unpack_s": f(&[2,4]),
                "forward_local_1d_s": f(&[15]), "forward_local_2d_s": f(&[9,11]),
                "forward_alltoall_s": f(&[13]), "forward_transpose_pack_unpack_s": f(&[12,14]),
                "grid_pack_unpack_s": f(&[0,6,10,16]), "potential_multiply_s": f(&[8]),
                "stages": dft_numerics::real::fft_profile::NAMES.iter().zip(s.fft.0).map(|(k,v)| (k.to_string(),json!({"s":v[0] as f64*1e-9,"calls":v[1]}))).collect::<serde_json::Map<_,_>>()
            },
            "calls":s.columns_per_call.len(), "columns":s.columns_per_call.iter().sum::<usize>(),
            "columns_per_call":s.columns_per_call, "complex_pairs_per_direction":s.pairs,
            "two_band_pairs":s.pairs-s.odd_tails, "single_band_tails":s.odd_tails,
            "nonlocal_wide_entry_calls":s.wide_calls - single, "nonlocal_direct_pure_calls":s.direct_calls
        });
        if single > 0 {
            row["nonlocal_single_entry_calls"] = json!(single);
        }
        row
    };
    let mut total = ApplyH::default();
    for s in &solves {
        for (a, b) in total.seconds.iter_mut().zip(s.seconds) {
            *a += b;
        }
        total.fft.add(s.fft);
        total
            .columns_per_call
            .extend_from_slice(&s.columns_per_call);
        total.pairs += s.pairs;
        total.odd_tails += s.odd_tails;
        total.wide_calls += s.wide_calls;
        total.direct_calls += s.direct_calls;
    }
    json!({"totals":block(&total),"solves":solves.iter().map(block).collect::<Vec<_>>(),
        "note":"Full diagnostics only; DAVIDSON_STAGES enables FFT clocks independently of the other profile switches; FFT measured=false means no timed FFT work was recorded on this rank. Rank-local seconds, outer rank max/min after summing solves. SCF, retry and verification solves in execution order. total is inclusive; FFT child timers are exclusive. Pairing is adjacent bands with a zero partner for odd tails. Wide entry is a request, not proof of MKL dispatch; inspect gemm_dispatch for effective routes."})
}

fn reference_solve_json(c: &dft_engine::refscf::record::SolveRecord) -> Value {
    json!({"ethr_ry":c.ethr_ry,"ethr_ha":c.ethr_ry/2.,"first_step_all_occupied":c.first_step,
        "davidson_iterations":c.iterations,"avg_iterations":c.iterations as f64,
        "active_bands_per_iteration":c.active_bands,"active_bands_sum":c.active_bands.iter().sum::<usize>(),
        "h_columns":c.h_columns,"reduced_problem_orders":c.reduced_orders,
        "restarts":c.restarts,"unconverged_bands":c.unconverged,
        "ritz_changes_ha":c.ritz_changes_ha,
        "inner_clocks":c.inner_clocks.map(|(total, stages)| json!({"total_s":total,
            "stage_s":dft_engine::refscf::record::INNER_STAGES.into_iter().zip(stages)
                .map(|(name, seconds)|(name.to_string(), json!(seconds))).collect::<serde_json::Map<_, _>>()})),
        "repairs":c.repairs,"dropped_corrections":c.dropped_corrections})
}
fn reference_record_json(r: &dft_engine::refscf::record::Record) -> Value {
    let mut record = json!({"stage":r.stage,"iteration":r.iteration,"dr2_ry":r.dr2_ry,"dr2_ha":r.dr2_ry/2.,
        "ethr_requested_ry":r.ethr_requested_ry,"ethr_requested_ha":r.ethr_requested_ry/2.,
        "ethr_ry":r.ethr_ry,"ethr_ha":r.ethr_ry/2.,"ethr_used_ry":r.ethr_ry,"beta":r.beta,
        "first_iteration_redo":r.first_iteration_redo,"history_occupancy":r.mixing.occupancy,
        "history_inserted_slot":r.mixing.inserted_slot,"gamma_l2":r.mixing.gamma_l2,
        "gamma_max":r.mixing.gamma_max,"a_condition_2":r.mixing.condition,
        "davidson_iterations":r.davidson_iterations,"h_columns":r.h_columns,
        "active_bands_sum":r.calls.iter().flat_map(|c|&c.active_bands).sum::<usize>(),
        "reduced_problem_orders":r.calls.iter().flat_map(|c|&c.reduced_orders).collect::<Vec<_>>(),
        "restarts":r.calls.iter().map(|c|c.restarts).sum::<usize>(),
        "calls":r.calls.iter().enumerate().map(|(i,c)|{let mut v=reference_solve_json(c);if let Some(k)=r.k_call_indices.get(i){v["k_index"]=json!(k);}v}).collect::<Vec<_>>(),
        "descf_ha":r.descf_ha,"descf_ry":2.*r.descf_ha,
        "variational_energy_ha":r.variational_energy_ha,"variational_energy_ry":2.*r.variational_energy_ha,
        "energies_ha":{"total":r.energy.total_hartree(),"band":r.energy.band_sum,
            "hartree":r.energy.hartree,"xc":r.energy.xc,"double_counting":r.energy.double_counting,
            "mixed_hartree":r.variational_hartree_ha,"mixed_xc":r.variational_xc_ha},
        "wall_s":{"map":r.map_s,"solve":r.solve_s,"density":r.density_s,"potential":r.potential_s,
            "energy":r.energy_s,"density_transform":r.density_transform_s,"metric":r.metric_s,
            "mix":r.mixing_s,"energy_correction":r.correction_s,"verification":r.verification_s}});
    if dft_numerics::switches::on(dft_numerics::switches::Id::REFSCF_RMM) {
        let solver = if r.calls.is_empty() {
            "davidson_reused_verification"
        } else if dft_engine::refscf::driver::rmm_scheduled(r.iteration) {
            "rmm"
        } else {
            "davidson"
        };
        record["solver"] = json!(solver);
        if let Some(calls) = record["calls"].as_array_mut() {
            for call in calls {
                call["solver"] = json!(solver);
            }
        }
    }
    record["map_start_s"] = json!(r.map_start_s);
    record["estimator_s"] = json!(r.estimator_s);
    record["verification_start_s"] = json!(r.verification_start_s);
    record["verification_end_s"] = json!(r.verification_end_s);
    record["map_end_s"] = json!(r.map_end_s);
    if let Some(c) = &r.verification_cost {
        record["verification_cost"] = json!({
            "gate_enabled":c.gate_enabled,
            "skipped":c.skipped,
            "predicted_error_ha":c.predicted_error_ha,
            "actual_error_ha":c.actual_error_ha,
            "reused_from":c.reused_from,
        });
        if c.first_gate_enabled {
            record["verification_cost"]["first_gate_enabled"] = json!(true);
            record["verification_cost"]["first_prediction_ha"] = json!(c.first_prediction_ha);
            record["verification_cost"]["contraction_prediction_ha"] =
                json!(c.contraction_prediction_ha);
        }
    }
    record
}

#[cfg(test)]
mod tests {

    #[test]
    fn strict_cli_ladder_off_preserves_oracle() {
        let dir = std::env::temp_dir().join(format!("r10a-oracle-cli-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut xml = crate::test_fixtures::gaussian_with(2., 1., true, 161)
            .replace("element=\"Xx\"", "element=\"He\"");
        let pp = dft_upf::parse_upf(&xml).unwrap();
        let charge = pp
            .mesh()
            .r()
            .iter()
            .map(|&r| {
                format!(
                    "{:.17e}",
                    8. / std::f64::consts::PI.sqrt() * r * r * (-r * r).exp()
                )
            })
            .collect::<Vec<_>>()
            .join(" ");
        let start = xml.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
        let end = xml.find("</PP_RHOATOM>").unwrap();
        xml.replace_range(start..end, &charge);
        fs::write(dir.join("syn.xml"), &xml).unwrap();
        let pp = dft_upf::parse_upf(&xml).unwrap();
        // The dev08 reproducer selected STRICT despite driver=legacy and the
        // standard-only keys, and explicitly disabled the legacy ladder.
        // Smaller synthetic geometry keeps this full CLI execution debug-fast.
        let base = "[system]\nname='R10A oracle'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=5.0\nfft_grid=[12,12,12]\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nmixing='pulay'\nmixing_beta=0.7\nmixing_history=8\nstandard_threshold_ha=5e-9\nverify='fixed_point'\nthreads=1\nmax_iterations=200\nnbands=5\nstart_density='atomic'\n";
        let mut densities = Vec::new();
        for (name, pin) in [("off", "cutoff_ladder='off'\n"), ("default", "")] {
            let text = format!(
                "{base}{pin}[performance]\nreference_strict=true\n[profile]\ndavidson_stages=false\n"
            );
            let config = Config::parse(&text).unwrap();
            let (_, _scope) = crate::run_config::for_run(&config).unwrap();
            let opts = options(&config).unwrap();
            assert!(opts.standard.is_none());
            assert_eq!(opts.eigen_tolerance, 1e-8);
            let cell = config.cell().unwrap();
            let atoms = [dft_geometry::Atom::new("He", [0.13, 0.21, 0.32])];
            let potentials = [SpeciesPotential {
                label: "He",
                potential: &pp,
            }];
            let problem = ScfProblem::with_kpoints(
                &cell,
                &atoms,
                &potentials,
                opts,
                &config.resolved_points().unwrap(),
                dft_engine::OccupationPolicy::Insulator {
                    minimum_gap_hartree: 1e-7,
                },
            )
            .unwrap();
            // Independent of CLI dispatch: these are the two pre-R10A entry
            // points. Check both so respecting "off" cannot disable the default.
            let expected = if name == "off" {
                problem.run_outcome().unwrap()
            } else {
                dft_engine::cutoff_ladder::run_oracle(problem, &potentials).unwrap()
            };
            assert!(expected.converged);
            let expected = expected.result;
            let input = dir.join(format!("{name}.toml"));
            let out = dir.join(name);
            fs::write(&input, text).unwrap();
            execute(&input, Some(&out)).unwrap();
            let result: Value =
                serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
            assert_eq!(result["scf"]["driver"], "strict_oracle");
            assert_eq!(result["single_path"]["effective_profile"], "strict");
            assert_eq!(
                result["total_energy"]["ha"].as_f64().unwrap().to_bits(),
                expected.total_hartree().to_bits(),
                "{name}"
            );
            assert_eq!(result["convergence"]["iterations"], expected.history.len());
            let history = result["convergence"]["history"].as_array().unwrap();
            for (actual, expected) in history.iter().zip(&expected.history) {
                assert_eq!(actual["davidson_iterations"], expected.davidson_iterations);
                assert_eq!(
                    actual["requested_eigen_tolerance_ha"]
                        .as_f64()
                        .unwrap()
                        .to_bits(),
                    expected.requested_eigen_tolerance.to_bits()
                );
            }
            let density = fs::read(out.join("density.bin")).unwrap();
            let expected: Vec<_> = expected
                .density
                .iter()
                .flat_map(|v| v.to_le_bytes())
                .collect();
            assert_eq!(density, expected, "{name}");
            if name == "off"
                && (std::env::consts::OS, std::env::consts::ARCH) == ("macos", "aarch64")
                && !dft_numerics::MKL_BACKEND
            {
                // Captured with the pre-R10A a5ce9cf4 CLI, Rust 1.99.1.
                // The dev08 Si8 regression was c03f5261b8862aca -> ...2ad4
                // at MPI8; platform/backend bits are not interchangeable.
                assert_eq!(
                    result["total_energy"]["ha"].as_f64().unwrap().to_bits(),
                    0xbfe3_f626_c809_494f
                );
                assert_eq!(result["convergence"]["iterations"], 11);
                assert_eq!(
                    history
                        .iter()
                        .map(|h| h["davidson_iterations"].as_u64().unwrap())
                        .collect::<Vec<_>>(),
                    [3, 2, 3, 2, 2, 3, 1, 3, 3, 3, 2]
                );
                assert_eq!(
                    sha256(&density),
                    "310f2763cbc5f728e40f510460c39eb8ed19326a4f2484ee0b7351285fe84a1e"
                );
            }
            densities.push(density);
        }
        assert_ne!(
            densities[0], densities[1],
            "fixture must exercise the ladder"
        );
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn reference_k_mesh_records_and_unsupported_combinations() {
        let dir = std::env::temp_dir().join(format!("g3-k-cli-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut pp = crate::test_fixtures::gaussian_with(2., 1.8, true, 161)
            .replace("element=\"Xx\"", "element=\"He\"");
        let parsed = dft_upf::parse_upf(&pp).unwrap();
        let array = |f: &dyn Fn(f64) -> f64| {
            parsed
                .mesh()
                .r()
                .iter()
                .map(|&r| format!("{:.17e}", f(r)))
                .collect::<Vec<_>>()
                .join(" ")
        };
        let start = pp.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
        let end = pp.find("</PP_RHOATOM>").unwrap();
        pp.replace_range(
            start..end,
            &array(&|r| {
                8. / (std::f64::consts::PI.sqrt() * 1.8_f64.powi(3))
                    * r
                    * r
                    * (-(r / 1.8).powi(2)).exp()
            }),
        );
        pp = pp.replace("number_of_wfc=\"0\"", "number_of_wfc=\"1\"").replace("</UPF>",
            &format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1></PP_PSWFC></UPF>", array(&|r| r * (-0.6 * r * r).exp())));
        fs::write(dir.join("syn.xml"), pp).unwrap();

        let base = "[system]\nname='k mesh'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=2.0\nkgrid=[2,2,2]\nkshift=[0.5,0.5,0.5]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=3\nmax_iterations=2\ndriver='reference'\nverify='off'\n";
        let config = Config::parse(base).unwrap();
        assert_eq!(
            reference_options(&config).unwrap().gamma_real,
            dft_engine::GammaReal::Off
        );
        {
            let config =
                Config::parse_unvalidated(&format!("{base}[spin]\nmode='collinear'\n")).unwrap();
            let failure = reference_options(&config).unwrap_err();
            assert!(failure.cause.contains("spin XC"));
        }
        for suffix in [
            "[spin]\nmode='noncollinear'\n",
            "[occupations]\nmode='gaussian'\nsmearing_ha=0.01\n",
            "[occupations]\nmode='methfessel_paxton'\nsmearing_ha=0.01\n",
            "[occupations]\nmode='cold'\nsmearing_ha=0.01\n",
        ] {
            assert!(Config::parse(&format!("{base}{suffix}")).is_err());
        }
        let input = dir.join("input.toml");
        fs::write(&input, base).unwrap();
        let out = dir.join("output");
        execute(&input, Some(&out)).unwrap();
        let result: Value =
            serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
        assert_eq!(result["kpoints"]["count"], 4);
        assert_eq!(result["gamma_real_used"], false);
        let rows = result["reference_scf"]["iterations"]
            .as_array()
            .expect("reference iterations");
        for row in rows {
            for call in row["calls"].as_array().unwrap() {
                assert!(call["k_index"].as_u64().is_some_and(|k| k < 4));
            }
        }
        // G6: the CLI passes its verified preparation and exports per-map residuals.
        let symmetric = format!(
            "{}[symmetry]\nmode='auto'\n",
            base.replace("7.0", "5.0")
                .replace("fractional=[0.13,0.21,0.32]", "fractional=[0.25,0.25,0.25]")
        );
        fs::write(&input, symmetric).unwrap();
        execute(&input, Some(&out)).unwrap();
        let sym: Value =
            serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
        assert_eq!(sym["reference_scf"]["symmetry"]["group_order"], 2);
        assert_eq!(sym["reference_scf"]["symmetry"]["irreducible_k_count"], 4);
        assert!(
            sym["reference_scf"]["symmetry"]["maps"]
                .as_array()
                .unwrap()
                .iter()
                .all(|m| m["residual_norm"].as_f64().unwrap().is_finite())
        );
        // G5: actual complex-k derivative output, including the three reference
        // cutoff scans, and the established stress-specific default threshold.
        let derivatives = format!(
            "{}[output]\nforces=true\nstress=true\n",
            base.replace("max_iterations=2", "max_iterations=100")
                .replace("verify='off'", "verify='fixed_point'")
        );
        let config = Config::parse(&derivatives).unwrap();
        assert_eq!(
            reference_profile(&config).unwrap().threshold_ry / 2.,
            dft_engine::StandardProfile::STRESS_THRESHOLD_HA
        );
        fs::write(&input, derivatives).unwrap();
        execute(&input, Some(&out)).unwrap();
        let result: Value =
            serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
        for name in ["forces", "stress"] {
            assert_eq!(result[name]["verification"], "verified");
            assert_eq!(
                result[name]["preconditions"]["fixed_point_verification_passed"],
                true
            );
        }
        assert_eq!(
            result["stress"]["pulay_assessment"]["scan_profile"],
            "standard_verified"
        );
        let symmetric_derivatives = format!(
            "{}[symmetry]\nmode='auto'\n",
            base.replace("7.0", "5.0")
                .replace("fractional=[0.13,0.21,0.32]", "fractional=[0.25,0.25,0.25]")
                .replace("max_iterations=2", "max_iterations=100")
                .replace("verify='off'", "verify='fixed_point'")
        );
        fs::write(
            &input,
            format!("{symmetric_derivatives}[output]\nforces=true\nstress=true\n"),
        )
        .unwrap();
        execute(&input, Some(&out)).unwrap();
        let symmetric_result: Value =
            serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
        assert_eq!(symmetric_result["forces"]["verification"], "verified");
        assert_eq!(symmetric_result["stress"]["verification"], "verified");
        assert!(symmetric_result["stress"]["pulay_assessment"].is_object());
        assert!(symmetric_result["stress"]["pulay_unavailable_reason"].is_null());
        // G4 routes through the distributed reference map, with explicit axes.
        for (name, spin, finite) in [
            ("fd", false, true),
            ("fixed-spin", true, false),
            ("spin-fd", true, true),
            ("spin-fd-free", true, true),
            ("spin-fd-moment", true, true),
        ] {
            // Fallback is allowed only for an implicit reference request.
            let mut text = base.replace("driver='reference'\n", "");
            if spin {
                text = text.replace(
                    "kshift=[0.5,0.5,0.5]",
                    "kshift=[0.5,0.5,0.5]\nktime_reversal=false",
                );
                text = text.replace(
                    "functional='lda-pz'",
                    "functional='lda-pw'\nallow_functional_override=true",
                );
                text.push_str(match name {
                    "spin-fd-free" => "[spin]\nmode='collinear'\n",
                    "spin-fd-moment" => "[spin]\nmode='collinear'\ntotal_moment=0.3\n",
                    _ => "[spin]\nmode='collinear'\nn_up=1\nn_down=1\n",
                });
            }
            if finite {
                text.push_str("[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\nband_cap=12\n");
            }
            fs::write(&input, &text).unwrap();
            if spin && finite {
                // R6a keeps moment FD on reference for every A/B switch combination.
                let route_config = Config::parse(&text).unwrap();
                for mask in 0..2 {
                    let scope = dft_numerics::switches::RunConfig::from_pairs(&[(
                        "refscf_spin_fd_moment",
                        if mask & 1 != 0 { "on" } else { "off" },
                    )])
                    .unwrap();
                    let _guard = scope.enter();
                    assert!(!reference_spin_fd_fallback(&route_config));
                    let pbe =
                        Config::parse(&text.replace("functional='lda-pw'", "functional='pbe'"))
                            .unwrap();
                    assert_eq!(
                        reference_spin_fd_fallback_reason(&pbe),
                        if name == "spin-fd-free" {
                            Some(G4D_FALLBACK)
                        } else {
                            None
                        }
                    );
                }
                if name == "spin-fd-free" {
                    // L1/G4D remains the default fallback for free-moment PBE.
                    let pbe_text = text.replace("functional='lda-pw'", "functional='pbe'");
                    fs::write(&input, pbe_text).unwrap();
                    let dest = dir.join("pbe-free-legacy");
                    execute(&input, Some(&dest)).unwrap();
                    let result: Value =
                        serde_json::from_slice(&fs::read(dest.join("results.json")).unwrap())
                            .unwrap();
                    assert_eq!(result["reference_driver_fallback_reason"], G4D_FALLBACK);
                    assert_eq!(
                        result["execution"]["reference_driver_fallback_reason"],
                        G4D_FALLBACK
                    );
                    assert_eq!(result["scf"]["driver"], "strict_oracle");
                    assert!(result.get("reference_scf").is_none());
                }
                // Exercise adopted defaults, without an explicit switch override.
                fs::write(&input, text).unwrap();
            }
            let dest = dir.join(name);
            execute(&input, Some(&dest)).unwrap();
            let result: Value =
                serde_json::from_slice(&fs::read(dest.join("results.json")).unwrap()).unwrap();
            assert_eq!(result["scf"]["driver"], "reference");
            assert_eq!(result["reference_scf"]["schema"], "reference-spin-fd-1");
            assert_eq!(result["band_data"]["shape"][0], if spin { 2 } else { 1 });
            assert_eq!(result["band_data"]["shape"][1], if spin { 8 } else { 4 });
            assert_eq!(result["verified"], false);
            assert_eq!(result["run"]["completed"], true);
            let ns = if spin { 2 } else { 1 };
            let nk = if spin { 8 } else { 4 };
            let nb = result["band_data"]["shape"][2].as_u64().unwrap() as usize;
            let values = result["occupations"]["values"].as_array().unwrap();
            assert_eq!(values.len(), ns * nk * nb);
            assert_eq!(
                result["band_data"]["eigenvalues_ha"]
                    .as_array()
                    .unwrap()
                    .len(),
                values.len()
            );
            let charge: f64 = values
                .iter()
                .enumerate()
                .map(|(i, f)| {
                    f.as_f64().unwrap()
                        * result["band_data"]["k_points"][(i / nb) % nk]["weight"]
                            .as_f64()
                            .unwrap()
                })
                .sum();
            assert!((charge - 2.).abs() < 1e-10);
            if finite {
                assert!(
                    result["occupations"]["chemical_potential"][0]["mu_ha"]
                        .as_f64()
                        .unwrap()
                        .is_finite()
                );
                if spin && name != "spin-fd-free" {
                    assert!(
                        result["occupations"]["chemical_potential"][1]["mu_ha"]
                            .as_f64()
                            .unwrap()
                            .is_finite()
                    );
                    let target = if name == "spin-fd-moment" { 0.3 } else { 0. };
                    assert!(
                        (result["spin"]["moment_mu_b"].as_f64().unwrap() - target).abs() < 1e-10
                    );
                }
            }
            assert_eq!(
                result["execution"]["effective_solver"],
                "reference_complex_davidson"
            );
            assert!(dest.join("density-total.cube").exists());
            if spin {
                assert!(dest.join("magnetization.cube").exists());
            }
            let e = result["thermodynamics"]["internal_energy_ha"]
                .as_f64()
                .unwrap();
            let ts = result["thermodynamics"]["entropy_term_ha"]
                .as_f64()
                .unwrap();
            let f = result["thermodynamics"]["free_energy_ha"].as_f64().unwrap();
            assert_eq!((e - ts).to_bits(), f.to_bits());
            for row in result["reference_scf"]["iterations"].as_array().unwrap() {
                for call in row["calls"].as_array().unwrap() {
                    assert!(call["spin_index"].as_u64().is_some_and(|s| s < if spin {
                        2
                    } else {
                        1
                    }));
                    assert!(
                        call["k_index"]
                            .as_u64()
                            .is_some_and(|k| k < if spin { 8 } else { 4 })
                    );
                }
            }
        }
        fs::remove_dir_all(dir).unwrap();
    }
    #[test]
    fn reference_lean_audit_serializes_unavailable_terms_as_null() {
        let energy = dft_engine::energy::EnergyAudit {
            xc_floor: None,
            nlcc: None,
            kinetic: f64::NAN,
            local: f64::NAN,
            local_zero_mode: f64::NAN,
            nonlocal: f64::NAN,
            hartree: f64::NAN,
            xc: f64::NAN,
            ewald: 1.,
            band_sum: 2.,
            double_counting: 3.,
            band_form: f64::NAN,
            harris_foulkes: 4.,
            functional_difference_bound: f64::NAN,
        };
        let value = super::energy_json(&energy);
        for key in [
            "kinetic",
            "local",
            "local_zero_mode",
            "nonlocal",
            "hartree",
            "xc",
            "band_form",
            "functional_difference_bound",
            "total",
        ] {
            assert!(value[key]["ha"].is_null());
            assert!(value[key]["ry"].is_null());
        }
        assert_eq!(value["band_sum"]["ha"], 2.);
        assert_eq!(value["double_counting"]["ha"], 3.);
        assert_eq!(value["harris_foulkes"]["ha"], 4.);
    }

    use super::*;
    #[test]
    fn reference_lazy_output_preserves_json_science_and_density_bytes() {
        use dft_mpi::Collective;
        use std::sync::Arc;
        let dir = std::env::temp_dir().join(format!("p3-nl-cli-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut pp = crate::test_fixtures::gaussian_with(2., 1.8, true, 161)
            .replace("element=\"Xx\"", "element=\"He\"");
        let parsed = dft_upf::parse_upf(&pp).unwrap();
        let array = |f: &dyn Fn(f64) -> f64| {
            parsed
                .mesh()
                .r()
                .iter()
                .map(|&r| format!("{:.17e}", f(r)))
                .collect::<Vec<_>>()
                .join(" ")
        };
        let start = pp.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
        let end = pp.find("</PP_RHOATOM>").unwrap();
        pp.replace_range(
            start..end,
            &array(&|r| {
                8. / (std::f64::consts::PI.sqrt() * 1.8_f64.powi(3))
                    * r
                    * r
                    * (-(r / 1.8).powi(2)).exp()
            }),
        );
        pp = pp.replace("number_of_wfc=\"0\"", "number_of_wfc=\"1\"").replace("</UPF>",
            &format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1></PP_PSWFC></UPF>", array(&|r| r * (-0.6 * r * r).exp())));
        fs::write(dir.join("syn.xml"), pp).unwrap();
        let base = "[system]\nname='lazy output test'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=2.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=3\nmax_iterations=2\ndriver='reference'\nverify='off'\n[performance]\n";
        for arm in ["off", "on"] {
            fs::write(
                dir.join(format!("{arm}.toml")),
                format!("{base}m2_lazy_result='{arm}'\n"),
            )
            .unwrap();
        }
        dft_mpi::simulate(2, |comm| {
            dft_engine::parallel::set_communicator(Some(Arc::new(comm.clone())));
            let mut baseline = None;
            for arm in ["off", "on"] {
                let out = dir.join(format!("{arm}-{}", comm.rank()));
                execute(&dir.join(format!("{arm}.toml")), Some(&out)).unwrap();
                let result: Value =
                    serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
                // Full files contain clocks, paths, hashes and switch/counter
                // provenance. Compare every scalar scientific output and the
                // complete per-iteration scientific history, excluding clocks.
                fn without_clocks(value: &mut Value) {
                    match value {
                        Value::Object(object) => {
                            object.retain(|key, _| {
                                !key.ends_with("_s") && key != "timings" && key != "inner_clocks"
                            });
                            for value in object.values_mut() {
                                without_clocks(value);
                            }
                        }
                        Value::Array(values) => {
                            for value in values {
                                without_clocks(value);
                            }
                        }
                        _ => {}
                    }
                }
                let mut history = result["convergence"]["history"].clone();
                without_clocks(&mut history);
                let mut science = serde_json::Map::new();
                science.insert("history".into(), history);
                for key in [
                    "cell",
                    "atoms",
                    "cutoff",
                    "basis",
                    "electrons",
                    "bands",
                    "occupations",
                    "xc",
                    "energies",
                    "total_energy",
                    "eigenvalues",
                    "verification",
                ] {
                    science.insert(key.to_string(), result[key].clone());
                }
                // Check solver work explicitly as well as the history above.
                science.insert(
                    "iterations".into(),
                    result["convergence"]["iterations"].clone(),
                );
                science.insert(
                    "converged".into(),
                    result["convergence"]["converged"].clone(),
                );
                for key in ["davidson_iterations_recorded", "h_vectors_recorded"] {
                    science.insert(key.into(), result["convergence"]["work"][key].clone());
                }
                let density = fs::read(out.join("density.bin")).unwrap();
                let meta = fs::read(out.join("density.json")).unwrap();
                if let Some((old_science, old_density, old_meta)) = &baseline {
                    assert_eq!(&science, old_science);
                    assert_eq!(&density, old_density);
                    assert_eq!(&meta, old_meta);
                } else {
                    baseline = Some((science, density, meta));
                }
            }
            dft_engine::parallel::set_communicator(None);
        });
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn reference_ladder_records_and_fine_stop_clock() {
        let dir = std::env::temp_dir().join(format!("p3-ladder-cli-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let mut pp = crate::test_fixtures::gaussian_with(2., 1.8, false, 161)
            .replace("element=\"Xx\"", "element=\"He\"");
        let parsed = dft_upf::parse_upf(&pp).unwrap();
        let array = |f: &dyn Fn(f64) -> f64| {
            parsed
                .mesh()
                .r()
                .iter()
                .map(|&r| format!("{:.17e}", f(r)))
                .collect::<Vec<_>>()
                .join(" ")
        };
        let start = pp.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
        let end = pp.find("</PP_RHOATOM>").unwrap();
        pp.replace_range(
            start..end,
            &array(&|r| {
                8. / (std::f64::consts::PI.sqrt() * 1.8_f64.powi(3))
                    * r
                    * r
                    * (-(r / 1.8).powi(2)).exp()
            }),
        );
        pp = pp.replace("number_of_wfc=\"0\"", "number_of_wfc=\"1\"").replace("</UPF>",
            &format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1></PP_PSWFC></UPF>", array(&|r| r * (-0.6 * r * r).exp())));
        fs::write(dir.join("syn.xml"), pp).unwrap();

        let base = "[system]\nname='ladder'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=4.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=3\nmax_iterations=80\ndriver='reference'\n";
        for verify in ["fixed_point", "off"] {
            let input = dir.join(format!("{verify}.toml"));
            fs::write(
                &input,
                format!("{base}verify='{verify}'\n[performance]\nrefscf_ladder='on'\n"),
            )
            .unwrap();
            let out = dir.join(verify);
            execute(&input, Some(&out)).unwrap();
            let result: Value =
                serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
            assert_eq!(result["convergence"]["converged"], true);
            let r = &result["reference_scf"];
            let memory = &r["memory_breakdown"];
            assert_eq!(memory["schema"], "reference-memory-1");
            assert_eq!(memory["complete"], false);
            assert_eq!(memory["units"], "bytes");
            for (names, rows) in [
                ("component_names", "components_per_rank"),
                ("davidson_names", "davidson_per_rank"),
            ] {
                let width = memory[names].as_array().unwrap().len();
                let ranks = memory[rows].as_array().unwrap();
                assert_eq!(ranks.len(), 1);
                assert_eq!(ranks[0].as_array().unwrap().len(), width);
                assert!(
                    ranks[0]
                        .as_array()
                        .unwrap()
                        .iter()
                        .all(|n| n.as_u64().is_some())
                );
            }
            assert!(!memory["unmeasured"].as_array().unwrap().is_empty());
            let kept = &r["kept_hx"];
            let hits = kept["cache_hits"].as_u64().unwrap();
            assert!(hits > 0);
            assert_eq!(kept["full_h_starting_columns_replaced"], 3 * hits);
            assert_eq!(kept["delta_v_columns"], 3 * hits);
            let bytes = kept["retained_slab_bytes_per_rank"].as_array().unwrap();
            assert_eq!(bytes.len(), 1);
            assert!(bytes[0].as_u64().unwrap() > 0);
            assert!(bytes[0].as_u64().unwrap() <= 32 * 1024 * 1024);
            let coarse = r["cutoff_ladder"]["iterations"].as_array().unwrap();
            let fine = r["iterations"].as_array().unwrap();
            assert!(!coarse.is_empty() && !fine.is_empty());
            assert!(coarse.iter().all(|r| r["stage"] == "coarse"));
            assert!(fine.iter().all(|r| r["stage"] == "fine"));
            assert_eq!(
                fine.len(),
                result["convergence"]["history"].as_array().unwrap().len()
            );
            assert_eq!(r["cutoff_ladder"]["fine_maps"], fine.len());
            assert_eq!(r["cutoff_ladder"]["coarse_maps"], coarse.len());
            assert!(
                r["wall_s"]["start"].as_f64().unwrap()
                    >= r["cutoff_ladder"]["coarse_seconds"].as_f64().unwrap()
            );
        }
        fs::remove_dir_all(dir).unwrap();
    }

    #[test]
    fn reference_output_writer_records_body_cost_without_duplicate_json() {
        let path = std::env::temp_dir().join(format!("p3-setup-json-{}.json", std::process::id()));
        let input = json!({"science":[1.,2.,3.],"reference_scf":{"setup_s":{"total":0.25},"iterations":[1,2]}});
        write_reference_results(&path, input.clone(), json!({"engine":{"total":0.1}})).unwrap();
        let bytes = fs::read(&path).unwrap();
        assert_eq!(bytes.last(), Some(&b'}'), "launch appender contract");
        let mut read: Value = serde_json::from_slice(&bytes).unwrap();
        let output = read["reference_scf"]
            .as_object_mut()
            .unwrap()
            .remove("output_s")
            .unwrap();
        assert_eq!(read, input);
        assert!(output["json_serialize_body"].as_f64().unwrap() > 0.);
        assert!(output["json_write_body"].as_f64().unwrap() > 0.);
        assert_eq!(output["engine"]["total"], 0.1);
        fs::remove_file(path).unwrap();
    }
    #[test]
    fn streamed_binary_hash_matches_whole_file_at_chunk_boundaries() {
        let path = std::env::temp_dir().join(format!("binary-hash-{}", std::process::id()));
        for size in [0, 1, 65_535, 65_536, 65_537, 131_081] {
            let bytes: Vec<u8> = (0..size).map(|i| (i % 251) as u8).collect();
            fs::write(&path, &bytes).unwrap();
            assert_eq!(binary_sha256(&path).unwrap(), sha256(&bytes));
        }
        fs::remove_file(&path).unwrap();
        assert!(binary_sha256(&path).is_err());
    }
    /// M2-STAGE-TIMERS: the davidson_stages rank merge takes the maximum and
    /// the minimum of every number over ranks and keeps integers integral.
    #[test]
    fn davidson_stages_rank_merge_takes_max_and_min() {
        use dft_mpi::Collective;
        let out = dft_mpi::simulate(3, |comm: &dft_mpi::SimRank| {
            let r = comm.rank();
            let v = json!({"note":"x","a":{"s":0.5 + r as f64,"calls":10 - r as u64},
                "list":[r as f64 * 2.0, -(r as f64)],
                "apply_h":{"totals":{"stage_s":{"kinetic":0.25+r as f64},
                    "columns_per_call":[3,2], "complex_pairs_per_direction":3},
                    "solves":[{"fft":{"inverse_s":0.5+r as f64,"forward_s":2.5-r as f64}}]}});
            over_ranks_with(comm, v).unwrap()
        });
        for o in &out {
            assert_eq!(o["apply_h"]["totals"]["stage_s"]["kinetic"], 2.25);
            assert_eq!(
                o["min_over_ranks"]["apply_h"]["totals"]["stage_s"]["kinetic"],
                0.25
            );
            assert_eq!(o["apply_h"]["totals"]["columns_per_call"], json!([3, 2]));
            assert_eq!(o["apply_h"]["solves"][0]["fft"]["forward_s"], 2.5);
            assert_eq!(
                o["min_over_ranks"]["apply_h"]["solves"][0]["fft"]["inverse_s"],
                0.5
            );
            assert_eq!(o["ranks"], 3);
            assert_eq!(o["note"], "x");
            assert_eq!(o["a"]["s"], 2.5);
            assert_eq!(o["a"]["calls"], 10);
            assert!(o["a"]["calls"].is_u64());
            assert_eq!(o["list"], json!([4.0, 0.0]));
            assert_eq!(o["min_over_ranks"]["a"]["s"], 0.5);
            assert_eq!(o["min_over_ranks"]["a"]["calls"], 8);
            assert_eq!(o["min_over_ranks"]["list"], json!([0.0, -2.0]));
        }
    }
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
            if let Ok(c) = Config::parse(&format!("{example}\n[scf]\n{fields}")) {
                assert_eq!(options(&c).unwrap_err().code, 1);
            }
        }
        let c = Config::parse(example).unwrap();
        assert_eq!(options(&c).unwrap().seed, 2026);
    }
    /// Review F02: the reported parallel model follows the solves that ran,
    /// not the rank count or the M2 switch. MPI ranks without a distributed
    /// eigensolve (more than one k point, spin, Fermi-Dirac, or a real-Gamma
    /// fallback) are `mpi_replicated` with the reason.
    #[test]
    fn parallel_model_reports_replicated_mpi_runs() {
        assert_eq!(parallel_model(1, 0, true, None), ("threads", None));
        assert_eq!(
            parallel_model(1, 0, true, Some("x".into())),
            ("threads", None)
        );
        assert_eq!(
            parallel_model(4, 0, true, Some(REPLICATED_MULTI_K.into())),
            ("mpi_replicated", Some(REPLICATED_MULTI_K.into()))
        );
        assert_eq!(
            parallel_model(2, 0, false, None),
            ("mpi_replicated", Some("no distributed eigensolve".into()))
        );
        assert_eq!(parallel_model(8, 12, true, None), ("mpi_stick_rows", None));
        assert_eq!(parallel_model(8, 12, false, None), ("mpi_row_block", None));
    }
}

/// Serialize run-local XC field diagnostics only at the output boundary.
fn add_xc_screening_json(result: &mut serde_json::Value) {
    if dft_physics::xc::screen_negative() {
        result["xc_screening"] = json!({
            "schema": "xc-screening-1",
            "density_threshold_bohr_minus3": dft_physics::xc::DENSITY_THRESHOLD,
            "maximum_screened_charge_fraction_per_channel": dft_physics::xc::SCREENED_CHARGE_FRACTION,
            "scope": "XC field evaluations, including input/output audits; global C-order grid indices",
            "records": dft_physics::xc::screening_records().iter().enumerate().map(|(i, r)| json!({
                "evaluation_index": i, "stage": r.stage, "channel": r.channel,
                "screened_points": r.points, "screened_absolute_charge_electrons": r.charge_electrons,
                "gate_total_electrons": r.total_electrons, "minimum_density_bohr_minus3": r.minimum_bohr3,
                "minimum_grid_index": r.minimum_grid_index
            })).collect::<Vec<_>>()
        });
    }
}
