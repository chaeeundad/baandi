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
    options_with_reference(
        config,
        dft_numerics::switches::on(dft_numerics::switches::Id::REFERENCE_STRICT),
    )
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
        // Standard-only controls act on the standard path only (the reference
        // switch ignores them and records so).
        o.standard_stall_guard = path.standard.is_some()
            && s.standard_stall_guard()
                .map_err(|e| error(1, e.to_string()))?;
        o.standard_density_forcing = path.standard.is_some()
            && s.standard_density_forcing()
                .map_err(|e| error(1, e.to_string()))?;
        o.hartree_mixing = s.hartree_mixing().map_err(|e| error(1, e.to_string()))?;
        o.spectral_guard = s.spectral_guard().map_err(|e| error(1, e.to_string()))?;
        o.preconditioner = s.preconditioner().map_err(|e| error(1, e.to_string()))?;
        o.davidson_capacity_per_band = s
            .davidson_capacity_per_band()
            .map_err(|e| error(1, e.to_string()))?;
        // MPI capacity rule (dev08 capacity sweep 2026-11-15): with four or more
        // ranks the replicated dense reduced eigenproblem is a scaling ceiling,
        // and a basis of 3 columns per band (instead of 6) shrinks it at a
        // small cost in Davidson iterations (si64 strict np=4/8 -6.8/-13.9 %,
        // standard np=8 -19.8 %, MgO 132 bands np=8 -5.2 %). An explicit
        // scf.davidson_capacity_per_band always wins.
        if s.davidson_capacity_per_band.is_none()
            && crate::parallel::size() >= 4
            && dft_numerics::switches::choice(dft_numerics::switches::Id::M2_CAPACITY) == "auto"
        {
            o.davidson_capacity_per_band = 3;
        }
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
            "reference_strict":switches::on(Id::REFERENCE_STRICT),
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
    if rc.consulted(Id::STOP_METRIC) {
        let c = rc.counters();
        let (kept, fft) = (
            c[RunCounter::StopKept as usize],
            c[RunCounter::StopFft as usize],
        );
        let path = match (kept > 0, fft > 0) {
            (true, false) => "kept",
            (false, true) => "fft",
            (true, true) => "kept+fft",
            (false, false) => "unused",
        };
        let why = (rc.choice(Id::STOP_METRIC) == "kept" && fft > 0)
            .then_some("some estimates transformed: density clamped or coefficients not kept");
        rc.note_path(Id::STOP_METRIC, path, why);
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

/// SHA-256 of the running executable (computed once, after the timed run).
pub(crate) fn binary_identity() -> Value {
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
/// ladder coarse stage + remainder (fine-session start, result assembly).
/// `clock` is the process clock accumulated during the SCF (coarse stage
/// included); the coarse part is subtracted for the fine-session values.
fn closed_timings(
    phases: &PhaseTimings,
    history: &[IterationRecord],
    setup: f64,
    total: f64,
    clock: dft_engine::phase_clock::PhaseClock,
    ladder: Option<&dft_engine::cutoff_ladder::CutoffLadderReport>,
) -> Value {
    let coarse = ladder.map(|l| l.coarse_clock).unwrap_or_default();
    let fine = clock.since(&coarse);
    let cert_phases: f64 = history
        .iter()
        .filter_map(|h| h.fixed_point_certificate.as_ref())
        .map(|c| phase_sum(&c.timings))
        .fold(0.0, |a, b| a + b);
    let cert_wall = fine.certificate.as_secs_f64();
    let cert_outside = (cert_wall - cert_phases).max(0.0);
    let phases_s = phase_sum(phases);
    let steps = fine.steps.as_secs_f64();
    let ladder_s = ladder.map_or(0.0, |l| l.coarse_time.as_secs_f64());
    let step_overhead = steps - phases_s - cert_outside;
    let remainder = total - setup - ladder_s - steps;
    json!({
        "certificate":{"blocks":fine.certificate_blocks,"wall_s":cert_wall,"screen_s":fine.screen.as_secs_f64(),
            "phases_s":cert_phases,"outside_phases_s":cert_outside,
            "note":"wall of the certificate blocks (screen, re-solve, checks); phases_s is already inside timings.phases"},
        "atomic_start_s":fine.atomic_start.as_secs_f64(),
        "steps_s":steps,
        "ladder_coarse":ladder.filter(|l| l.coarse_time > std::time::Duration::ZERO).map(|l| json!({
            "setup_s":l.coarse_setup.as_secs_f64(),"session_s":l.coarse_session.as_secs_f64(),
            "prolong_s":l.coarse_prolong.as_secs_f64(),"total_s":l.coarse_time.as_secs_f64(),
            "phases":phase_json(&l.coarse_phases),"phases_sum_s":phase_sum(&l.coarse_phases),
            "atomic_start_s":l.coarse_clock.atomic_start.as_secs_f64(),"steps_s":l.coarse_clock.steps.as_secs_f64(),
            "certificate_s":l.coarse_clock.certificate.as_secs_f64()})),
        "closure":{"setup_s":setup,"phases_s":phases_s,"certificate_outside_phases_s":cert_outside,
            "step_overhead_s":step_overhead,"ladder_coarse_s":ladder_s,"remainder_s":remainder,"total_s":total,
            "note":"total = setup + phases + certificate_outside_phases + step_overhead + ladder_coarse + remainder; atomic_start_s is inside phases (potential of the first map)"}
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
    if let Some(a) = &r.advance {
        // BREAKTHROUGH-3: bounded subspace update of this map (research diagnostic).
        result["advance"] = json!({"budget_blocks":a.budget_blocks,"blocks_used":a.blocks_used,"h_columns":a.h_columns,
            "kept_start":a.kept_start,"unconverged_roots":a.unconverged_roots,"inner_estimate_ha":a.inner_estimate_ha,
            "scf_error_ha":a.scf_error_ha,"inner_scf_ratio":if a.scf_error_ha > 0.0 {a.inner_estimate_ha / a.scf_error_ha} else {f64::INFINITY},
            "completion":a.completion.map(|t| t.name()),"completion_h_columns":a.completion_h_columns});
    }
    result
}

/// BREAKTHROUGH-3 (`BAANDI_SCF_ADVANCE`) counters of this process (count only).
fn scf_advance_json() -> serde_json::Value {
    let c = dft_engine::advance_counters();
    json!({
        "mode":dft_numerics::switches::choice(dft_numerics::switches::Id::SCF_ADVANCE),
        "starting_blocks":dft_numerics::switches::number(dft_numerics::switches::Id::SCF_ADVANCE_BLOCKS),
        "inner_ratio_rule":dft_engine::ADVANCE_INNER_RATIO,
        "maps_advanced":c.maps_advanced,"maps_exact":c.maps_exact,"declined":c.declined,
        "completions":c.completions(),"completion_candidate":c.completion_candidate,"completion_inner_noise":c.completion_inner_noise,
        "completion_growth":c.completion_growth,"completion_gap":c.completion_gap,
        "advance_h_columns":c.advance_h_columns,"completion_h_columns":c.completion_h_columns,"advance_blocks":c.advance_blocks,
        "kept_starts":c.kept_starts,"converged_within_budget":c.converged_within_budget,
        "ratio_below_rule":c.ratio_below_rule,"ratio_above_rule":c.ratio_above_rule,
        "note":"maps_advanced = standard maps run as a bounded subspace update (Advance), maps_exact = eligible-path maps that ran the converging solve (first map, redo, candidates), declined = requests the solver could not honour; completions = advanced maps followed by the exact solve at the same density, by trigger; advance_h_columns / completion_h_columns = operator columns of the advance calls and of their completions (both inside h_vectors_recorded); ratio_* = advanced maps whose B / E was at or below / above the rule"
    })
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
    launch::mark("run_entry");
    let input_clock = launch::clock();
    let bytes = fs::read(input).map_err(|e| error(5, format!("{}: {e}", input.display())))?;
    let text = std::str::from_utf8(&bytes).map_err(|e| error(1, e))?;
    let config = Config::parse(text).map_err(|e| error(1, e))?;
    launch::elapsed("input_read_parse", input_clock);
    let config_clock = launch::clock();
    // RUN-CONFIG: this run's switch configuration (the process object in the
    // CLI binary; a thread scope for in-process callers with other values).
    let (run_config, _run_scope) = crate::run_config::for_run(&config).map_err(|e| error(1, e))?;
    launch::elapsed("run_config", config_clock);
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
    let preflight_clock = launch::clock();
    let mut preflight = crate::config::preflight::check(&config, input).map_err(|e| error(1, e))?;
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
    let opts = options(&config)?;
    let scf_path = config.scf_path().map_err(|e| error(1, e.to_string()))?;
    let mut potentials = Vec::new();
    let mut species = Vec::new();
    let mut inputs = vec![fs::canonicalize(input)?];
    let upf_clock = launch::clock();
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
        v
    };
    if let Some(theta) = config.finite_temperature().map_err(|e| error(1, e))? {
        let record = driver_record("Fermi-Dirac occupations");
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
    launch::elapsed("engine_construction", constructor_clock);
    launch::mark("engine_ready");
    let setup = start.elapsed().as_secs_f64();
    // BREAKTHROUGH-0-INSTRUMENT: problem construction is the `init` bucket.
    dft_engine::profile::add_init(setup);
    let clock_start = dft_engine::phase_clock::snapshot();
    let ladder = crate::config::Scf::cutoff_ladder(config.scf.as_ref())
        .map_err(|e| error(1, e.to_string()))?;
    let mut ladder_report = None;
    let atomic_density_built = problem.atomic_superposition_density().is_some();
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
                    run_config.freeze_consulted();
                    note_switch_paths(run_config);
                    failure["scf"] = options_json(&opts);
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
    let r = outcome.result;
    // RUN-CONFIG: what the computation consulted (result writers below read
    // more switches for their report blocks; those reads do not count).
    run_config.freeze_consulted();
    // Everything from here to the results file is output materialization.
    let io_clock = Instant::now();
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
    let timings = closed_timings(
        &phases,
        &r.history,
        setup,
        total,
        dft_engine::phase_clock::snapshot().since(&clock_start),
        ladder_report.as_ref(),
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
                "davidson_stages":over_ranks(davidson_stages_json())?,
                "b1":b1_json(),
                "b5":{"mode":dft_numerics::switches::choice(dft_numerics::switches::Id::B5),"bordered":dft_engine::b5_counters().bordered,"declined":dft_engine::b5_counters().declined,"full":dft_engine::b5_counters().full},
                "ortho":{"mode":dft_numerics::switches::text(dft_numerics::switches::Id::ORTHO).unwrap_or("dgks:0.1 (default)"),"second_pass_skipped":dft_engine::ortho_counters().0,"second_pass_repeated":dft_engine::ortho_counters().1},
                "algebra_stream":algebra_stream_json(),
                "davidson_lean":davidson_lean_json(),
                "reuse_psi":reuse_psi_json(),
                "scf_advance":scf_advance_json(),
                "gemm_dispatch":gemm_dispatch_json(),
                "h_vectors_recorded":r.history.iter().map(|h|h.operator_applications).sum::<usize>(),
                "certificate_h_vectors":r.history.iter().filter_map(|h|h.fixed_point_certificate.as_ref()).map(|c|c.operator_applications).sum::<usize>(),
                "note":"recorded = per-iteration solves in history; certificate = fixed-point certificate re-solves, reported only inside each history entry's fixed_point_certificate; total work is the sum"}},"energies":energies,
        "total_energy":{"ha":r.energy.total_hartree(),"ry":from_hartree(r.energy.total_hartree(), EnergyUnit::Rydberg),"ev":from_hartree(r.energy.total_hartree(), EnergyUnit::ElectronVolt)},
        "eigenvalues":{"ha":r.eigenvalues,"ev":r.eigenvalues.iter().map(|e|from_hartree(*e, EnergyUnit::ElectronVolt)).collect::<Vec<_>>()},
        "timings":{"phases":phase_json(&phases),"setup_s":setup,"engine_setup_s":r.setup_time.as_secs_f64(),"engine_total_s":r.total_time.as_secs_f64(),"total_s":total,"fft_hpsi":"included in enclosing phases; not separately measured",
            "certificate":timings["certificate"],"atomic_start_s":timings["atomic_start_s"],"steps_s":timings["steps_s"],"ladder_coarse":timings["ladder_coarse"],"closure":timings["closure"]},
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
            "standard":l.standard,"standard_seed_ha":l.standard_seed_ha,
            "secant_requested":l.secant_requested,"secant_seeded":l.secant_seeded,
            "secant_fallback":l.secant_fallback,
            "coarse_stop_rule":l.coarse_rule.map(|r| format!("{r:?}").to_lowercase()),
            "coarse_stop_level":l.coarse_stop_level,
            "atomic_truncation_distance":l.atomic_truncation_distance,
            "handoff_to_final_density_distance":l.handoff_to_final_density_distance,
            "first_tolerance":l.first_tolerance.map(|r| format!("{r:?}").to_lowercase()),
            "strict_seed_ha":l.strict_seed_ha,
            "fine_scf_steps":l.fine_iterations,"fine_h_vectors":l.fine_h_vectors,
            "fine_certificate_h_vectors":l.fine_certificate_h_vectors,
            "coarse_fine_step_equivalents":l.coarse_fine_step_equivalents(),
            "coarse_wall_fine_step_equivalents":(l.used && l.fine_iterations > 0).then(|| {
                let steps = timings["steps_s"].as_f64().unwrap_or(f64::NAN)
                    - timings["certificate"]["wall_s"].as_f64().unwrap_or(0.0);
                l.coarse_time.as_secs_f64() / (steps / l.fine_iterations as f64)
            }).filter(|x| x.is_finite()),
            "total_step_equivalents":l.coarse_fine_step_equivalents().map(|c| c + l.fine_iterations as f64),
            "total_h_full_grid_equivalents":l.coarse_h_equivalents() + (l.fine_h_vectors + l.fine_certificate_h_vectors) as f64,
            "coarse_history":l.coarse_history.iter().map(history_json).collect::<Vec<_>>(),
            "note":"coarse stage supplies the start only; convergence, energy audit and certificate are at the full cutoff (history and work above are full-cutoff only; coarse_history holds the coarse maps, energies at the coarse cutoff). D1: coarse_fine_step_equivalents = coarse_h_full_grid_equivalents / (fine_h_vectors / fine_scf_steps) (recorded H per fine step, certificate excluded); coarse_wall_fine_step_equivalents = coarse_stage_s / ((timings.steps_s - timings.certificate.wall_s) / fine_scf_steps); total_step_equivalents = fine_scf_steps + coarse_fine_step_equivalents; total_h_full_grid_equivalents = coarse equivalents + fine recorded + fine certificate H; distances in the density-residual norm (Bohr^-3/2)"});
    }
    // SCF-TRAJECTORY (ISSUE: start_density looked dead on dev08): the start
    // each stage actually used. With the ladder the coarse stage is the first
    // stage and takes the requested start; the full-cutoff stage always starts
    // from the prolonged coarse density and orbitals, so the top-level
    // `atomic_start` block (the full-cutoff session's) reads "supplied initial
    // vectors retained" by design.
    result["effective_start"] = match ladder_report.as_ref().filter(|l| l.used) {
        Some(l) => json!([
            {"stage":"coarse","density":l.coarse_start_density,
             "wavefunctions":opts.start_wavefunctions.name(),
             "atomic_orbitals_used":l.coarse_atomic_start.as_ref().map(|a| a.atomic_orbitals),
             "atomic_start_h_applications":l.coarse_atomic_start.as_ref().map(|a| a.initial_h_applications),
             "atomic_fallback_reason":l.coarse_atomic_start.as_ref().and_then(|a| a.fallback_reason.clone())},
            {"stage":"fine","density":"prolonged coarse density","wavefunctions":"prolonged coarse orbitals"}]),
        None => json!([
            {"stage":"single","density":if atomic_density_built {"atomic"} else {"uniform"},
             "wavefunctions":opts.start_wavefunctions.name(),
             "atomic_orbitals_used":(opts.start_wavefunctions == dft_engine::StartWavefunctions::Atomic).then_some(r.atomic_start.atomic_orbitals),
             "atomic_start_h_applications":(opts.start_wavefunctions == dft_engine::StartWavefunctions::Atomic).then_some(r.atomic_start.initial_h_applications),
             "atomic_fallback_reason":r.atomic_start.fallback_reason}]),
    };
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
        result["verification"] = verification_json(v);
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
    let replicated_reason = if points.len() > 1 {
        Some(REPLICATED_MULTI_K.to_owned())
    } else {
        r.gamma_real_fallback_reason.as_ref().map(|reason| {
            format!("real-Gamma path not taken ({reason}): every rank runs the complex solver")
        })
    };
    let (model, replicated) = parallel_model(
        crate::parallel::size(),
        dft_engine::parallel::distributed_solves(),
        dft_engine::parallel::m2_enabled(),
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
    if crate::parallel::size() > 1 {
        let m2_on = dft_engine::parallel::m2_enabled();
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
    launch::timed("volumetric_write", || {
        volumetric_out::write_outputs(&config, &prepared, &refs, &r, &out)
    })
    .map_err(|e| error(5, e))?;
    launch::timed("density_bin_write", || {
        fs::write(out.join("density.bin"), density_bytes)
    })?;
    launch::timed("density_json_write", || {
        write_json(&out.join("density.json"), &density_meta)
    })?;
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
    // and first-block diagnostic (BAANDI_FIRST_BLOCK=on); both absent otherwise.
    if launch::clock().is_some() {
        launch::elapsed("materialize_including_outputs", Some(io_clock));
    }
    dft_engine::profile::add_io(io_clock.elapsed().as_secs_f64());
    launch::mark("profile_finish");
    if let Some(mut mine) = dft_engine::profile::finish_rank(crate::parallel::rank()) {
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
    let first_block = dft_engine::first_block::take_records();
    if !first_block.is_empty() {
        console.push_str(&dft_engine::first_block::table(&first_block));
        result["convergence"]["work"]["first_block"] = first_block_json(&first_block);
    }
    let manifest_clock = launch::clock();
    note_switch_paths(run_config);
    result["scf"]["switches"] = switches_json();
    launch::elapsed("switch_manifest", manifest_clock);
    let run_block_clock = launch::clock();
    result["run"] = run_block(
        run_config,
        if outcome.converged {
            "converged"
        } else {
            "not_converged"
        },
    );
    launch::elapsed("run_block", run_block_clock);
    let result_path = out.join("results.json");
    let serialize_clock = launch::clock();
    let data = serde_json::to_vec_pretty(&result).map_err(|e| error(5, e))?;
    launch::elapsed("results_serialize", serialize_clock);
    launch::timed("results_write", || fs::write(&result_path, data))?;
    launch::output(&result_path);
    launch::mark("outputs_complete");
    Ok(RunReport {
        console,
        converged: outcome.converged,
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

/// Exclusive phase profile (`BAANDI_PHASE_PROFILE=on`): per bucket the maximum
/// and minimum over ranks of every entry (seconds, then counts).
pub fn profile_json(report: &dft_engine::profile::Report, run_wall: f64) -> serde_json::Value {
    let covered = report.covered_wall();
    let closure = report.closure();
    json!({
        "schema": "phase-profile-2",
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

/// First-block residual records (`BAANDI_FIRST_BLOCK=on`).
pub fn first_block_json(records: &[dft_engine::first_block::Record]) -> serde_json::Value {
    json!(records.iter().map(|r| json!({
        "stage": r.stage, "kind": r.kind, "iteration": r.iteration,
        "bands": r.bands, "occupied": r.occupied, "tolerance_ha": r.tolerance,
        "buffer_tolerance_ha": r.buffer_tolerance, "kept_hx": r.kept_hx,
        "overlap_defect": r.overlap_defect, "boundary_gap_ha": r.boundary_gap,
        "ritz_values_ha": r.ritz_values, "residual_norms": r.residual_norms,
        "electrons": r.electrons, "occupations": r.occupations,
        "mandatory_roots": r.mandatory_roots, "density_carrying_roots": r.occupied,
        "empty_buffer_roots": r.empty_buffer_roots,
        "convergence_criterion": r.criterion.name(), "root_thresholds_ha": r.root_thresholds,
        "ritz_changes_ha": r.ritz_changes,
        "diagnostic_operator_applications": r.diagnostic_operator_applications,
        "residual_unconverged_all": r.residual_unconverged_all,
        "residual_unconverged_occupied": r.residual_unconverged_occupied,
        "ritz_first_correction_capacity": r.ritz_first_correction_capacity,
        "singular_values": r.singular_values,
        "rank_1e1": r.rank_1e1, "rank_1e2": r.rank_1e2, "rank_1e3": r.rank_1e3, "rank_1e4": r.rank_1e4,
        "rank_above_residual_tolerance": r.rank_above_residual_tolerance,
    })).collect::<Vec<_>>())
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

#[cfg(test)]
mod tests {
    use super::*;
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
                "list":[r as f64 * 2.0, -(r as f64)]});
            over_ranks_with(comm, v).unwrap()
        });
        for o in &out {
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
            let c = Config::parse(&format!("{example}\n[scf]\n{fields}")).unwrap();
            assert_eq!(options(&c).unwrap_err().code, 1);
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
