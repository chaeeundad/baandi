//! File/console boundary for the sequential engine. No I/O occurs in SCF kernels.
use crate::{config::Config, grid_report};
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
fn engine_error(e: ScfError) -> RunError {
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
        functional: config.xc.functional.engine(),
        allow_functional_override: config.xc.allow_functional_override,
        cutoff_hartree: config.basis.ecutwfc_ry / 2.,
        grid_override: config.basis.fft_grid,
        ..ScfOptions::default()
    };
    if let Some(s) = &config.scf {
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
    json!({"threads":o.threads,"cutoff_hartree":o.cutoff_hartree,"nbands":o.bands,"extra_bands":o.extra_bands,"grid_override":o.grid_override,"mixing":mixing,
        "density_residual_tolerance":o.density_tolerance,"energy_tolerance_ha":o.energy_tolerance,"eigen_tolerance":o.eigen_tolerance,
        "residual_eigen_schedule":o.residual_eigen_schedule,"free_electron_initial_guess":o.free_electron_initial_guess,"initial_eigen_tolerance":o.initial_eigen_tolerance,"eigen_tolerance_decay":o.eigen_tolerance_decay,"max_iterations":o.max_iterations,
        "davidson_max_iterations":o.davidson_max_iterations,"davidson_subspace":o.davidson_subspace,"davidson_subspace_min_per_band":6,"seed":o.seed,
        "divergence_factor":o.divergence_factor,"divergence_patience":o.divergence_patience})
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
    Value::Object(m)
}
fn phase_json(t: &PhaseTimings) -> Value {
    json!({"potential_s":t.potential.as_secs_f64(),"diagonalization_s":t.diagonalization.as_secs_f64(),"density_s":t.density.as_secs_f64(),"total_energy_s":t.energy.as_secs_f64(),"mixing_s":t.mixing.as_secs_f64()})
}
fn phase_sum(t: &PhaseTimings) -> f64 {
    (t.potential + t.diagonalization + t.density + t.energy + t.mixing).as_secs_f64()
}
fn history_json(r: &IterationRecord) -> Value {
    json!({"iteration":r.iteration,"energies":energy_json(&r.energy),"delta_energy_ha":r.energy_change,"density_residual":r.density_residual,
        "max_eigen_residual_ha":r.max_eigen_residual,"requested_eigen_tolerance_ha":r.requested_eigen_tolerance,"davidson_iterations":r.davidson_iterations,"operator_applications":r.operator_applications,"davidson_iterations_per_k":r.davidson_iterations_per_k,"operator_applications_per_k":r.operator_applications_per_k,
        "relative_charge_error":r.charge_error,"relative_input_charge_error":r.input_charge_error,"min_density_electrons_per_bohr3":r.min_density,
        "mixing_fraction":r.mixing_fraction,"mixer_reset":r.mixer_reset,"density_contraction":r.density_contraction,"mixing_history_used":r.mixing_history_used,"mixing_resets":r.mixing_resets,"mixing_raw_recovery":r.mixing_raw_recovery,"timings":phase_json(&r.timings),
        "fixed_point_certificate":r.fixed_point_certificate.as_ref().map(|c| json!({"strict_residuals_per_k_ha":c.strict_residuals_per_k,"buffer_residuals_per_k_ha":c.buffer_residuals_per_k,"density_residual":c.density_residual,"occupied_eigenvalue_change_ha":c.occupied_eigenvalue_change,"davidson_iterations":c.davidson_iterations,"operator_applications":c.operator_applications,"davidson_iterations_per_k":c.davidson_iterations_per_k,"operator_applications_per_k":c.operator_applications_per_k,"timings":phase_json(&c.timings)}))})
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
    for name in ["results.json", "density.json", "density.bin"] {
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
    let start = Instant::now();
    let bytes = fs::read(input).map_err(|e| error(5, format!("{}: {e}", input.display())))?;
    let text = std::str::from_utf8(&bytes).map_err(|e| error(1, e))?;
    let config = Config::parse(text).map_err(|e| error(1, e))?;
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
    println!(
        "XC: {}; allow_functional_override={}; override_applied={}",
        opts.functional.name(),
        opts.allow_functional_override,
        potentials
            .iter()
            .any(|(_, p)| !opts.functional.matches_upf(p.functional()))
    );
    let out = output_override.map(Path::to_owned).unwrap_or_else(|| {
        PathBuf::from(
            config
                .output
                .as_ref()
                .and_then(|o| o.directory.as_deref())
                .unwrap_or("results"),
        )
    });
    fs::create_dir_all(&out)?;
    protect_outputs(&out, &inputs)?;
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
    let points = config.basis.resolve_kpoints().map_err(|e| error(1, e))?;
    let problem = if config.basis.legacy_gamma() {
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
    let electrons = problem.electron_count();
    let setup = start.elapsed().as_secs_f64();
    let outcome = problem.run_outcome().map_err(engine_error)?;
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
    let mut result = json!({"schema_version":1,"tool":{"name":"dft-cli","version":env!("CARGO_PKG_VERSION"),"git_commit":env!("DFT_GIT_COMMIT"),"rustc_version":env!("DFT_RUSTC_VERSION")},
        "timestamp":timestamp(SystemTime::now().duration_since(UNIX_EPOCH).map_err(|e|error(5,e))?.as_secs()),"input_sha256":sha256(&bytes),"species":species,
        "cell":{"lattice_columns_bohr":cell,"volume_bohr3":volume,"reciprocal_columns_inverse_bohr":prepared.cell.reciprocal_columns_inverse_bohr()},
        "atoms":config.system.atoms.iter().map(|a|json!({"species":a.species,"fractional":a.fractional,"cartesian_bohr":prepared.cell.fractional_to_cartesian_bohr(a.fractional)})).collect::<Vec<_>>(),
        "cutoff":{"ecutwfc_ry":config.basis.ecutwfc_ry,"ecutwfc_ha":opts.cutoff_hartree},
        "basis":{"n_plane_waves":npw,"n_density_g":prepared.density.vectors().len(),"fft_grid":prepared.grid.dims(),"fft_grid_overridden":prepared.grid.overridden()},
        "electrons":electrons,"bands":bands,"occupations":r.occupations,"xc":opts.functional.name(),"scf":options_json(&opts),
        "convergence":{"converged":outcome.converged,"iterations":r.history.len(),"last_iteration":last,"history":history},"energies":energies,
        "total_energy":{"ha":r.energy.total_hartree(),"ry":from_hartree(r.energy.total_hartree(), EnergyUnit::Rydberg),"ev":from_hartree(r.energy.total_hartree(), EnergyUnit::ElectronVolt)},
        "eigenvalues":{"ha":r.eigenvalues,"ev":r.eigenvalues.iter().map(|e|from_hartree(*e, EnergyUnit::ElectronVolt)).collect::<Vec<_>>()},
        "timings":{"phases":phase_json(&phases),"setup_s":setup,"engine_setup_s":r.setup_time.as_secs_f64(),"engine_total_s":r.total_time.as_secs_f64(),"total_s":total,"fft_hpsi":"included in enclosing phases; not separately measured"},
        "threading":{"threads":opts.threads.max(1).min(r.k_points.len()),"threads_requested":opts.threads,"threads_used":opts.threads.max(1).min(r.k_points.len()),"outer_threads":opts.threads.max(1).min(r.k_points.len()),"inner_threads":opts.threads.max(1) / opts.threads.max(1).min(r.k_points.len()),"backend":"std::thread::scope over k points and bands; sequential faer and serial FFT per worker"},
        "memory_estimate":{"basis_bytes":npw*std::mem::size_of_val(&prepared.wavefunctions.vectors()[0]),"density_basis_bytes":prepared.density.vectors().len()*std::mem::size_of_val(&prepared.density.vectors()[0]),"density_grid_bytes":ngrid*8,"one_complex_grid_bytes":ngrid*16,"orbital_bytes":npw*bands*16,"peak_rss_bytes":null,"note":"array payload estimates only; scratch, projectors, allocator overhead excluded; peak RSS unavailable in-process"},
        "units":"Ha, Ry, eV, Bohr, seconds and bytes as keyed; density_residual in Bohr^-3/2; electrons/occupations in electron counts; relative errors and fractions dimensionless; see docs/cli.md"});
    result["band_accuracy"] = json!(r.k_points[0].solution.residuals.iter()
        .zip(&r.k_points[0].solution.buffer).zip(&r.k_points[0].solution.converged)
        .map(|((&residual, &buffer), &converged)| json!({"residual_ha":residual,"buffer":buffer,"converged":converged}))
        .collect::<Vec<_>>());
    // Preserve the legacy LDA schema unless an override was explicitly enabled.
    if opts.functional == dft_engine::Functional::Pbe || opts.allow_functional_override {
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
        writeln!(console, "HOMO/LUMO gap: direct {:.12} Ha = {:.12} eV at k{}; indirect {:.12} Ha = {:.12} eV (HOMO k{}, LUMO k{})",
            direct.1, from_hartree(direct.1, EnergyUnit::ElectronVolt), direct.0,
            lumo.1-homo.1, from_hartree(lumo.1-homo.1, EnergyUnit::ElectronVolt), homo.0, lumo.0).unwrap();
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
    fs::write(out.join("density.bin"), density_bytes)?;
    write_json(&out.join("density.json"), &density_meta)?;
    write_json(&out.join("results.json"), &result)?;
    Ok(RunReport {
        console,
        converged: outcome.converged,
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
