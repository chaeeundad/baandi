//! Fixed-cell ionic driver. SCF kernels and derivative ownership stay in engine.
mod cell;
pub mod optimizer;
use crate::{
    config::Config,
    run::{RunError, options, sha256},
};
use dft_engine::{
    CompletionMode, ExecutionConfig, ScfOptions, ScfProblem, ScfSession, SpeciesPotential,
    StepStatus,
};
use dft_geometry::{Atom, Cell};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};
use std::{fs, io::Write, path::Path};

#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum Method {
    #[default]
    Bfgs,
    Fire,
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Fixed {
    /// Zero-based input atom index; omitted components freezes all components.
    pub atom: usize,
    #[serde(default = "all_components")]
    pub components: [bool; 3],
}
fn all_components() -> [bool; 3] {
    [true; 3]
}
#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Settings {
    pub method: Method,
    pub mode: cell::Mode,
    pub cell_dofs: cell::Dofs,
    pub fixed_axes: [bool; 3],
    pub pressure_ha_bohr3: f64,
    pub stress_tol_ha_bohr3: f64,
    pub pulay_warning_ha_bohr3: f64,
    pub rebuild_every: usize,
    pub force_tol_ha_bohr: f64,
    pub energy_tol_ha: f64,
    pub max_steps: usize,
    pub max_displacement_bohr: f64,
    pub trust_radius_bohr: f64,
    pub warm_start: bool,
    pub fixed: Vec<Fixed>,
}
impl Default for Settings {
    fn default() -> Self {
        Self {
            method: Method::Bfgs,
            mode: cell::Mode::Ions,
            cell_dofs: cell::Dofs::Full,
            fixed_axes: [false; 3],
            pressure_ha_bohr3: 0.,
            stress_tol_ha_bohr3: 1e-6,
            pulay_warning_ha_bohr3: 5e-6,
            rebuild_every: 0,
            force_tol_ha_bohr: 1e-3,
            energy_tol_ha: 1e-8,
            max_steps: 200,
            max_displacement_bohr: 0.2,
            trust_radius_bohr: 0.2,
            warm_start: true,
            fixed: vec![],
        }
    }
}
fn error(code: u8, e: impl std::fmt::Display) -> RunError {
    RunError {
        code,
        cause: e.to_string(),
    }
}
impl Settings {
    pub(crate) fn validate(&self, n: usize) -> Result<(), RunError> {
        if [
            self.force_tol_ha_bohr,
            self.stress_tol_ha_bohr3,
            self.pulay_warning_ha_bohr3,
            self.energy_tol_ha,
            self.max_displacement_bohr,
            self.trust_radius_bohr,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.)
            || !self.pressure_ha_bohr3.is_finite()
            || self.max_steps == 0
            || self.fixed.iter().any(|f| f.atom >= n)
        {
            return Err(error(
                1,
                "relax: positive finite tolerances/radii, positive max_steps and valid zero-based fixed atom indices required",
            ));
        }
        Ok(())
    }
    fn project(&self, v: &mut [f64]) {
        for fixed in &self.fixed {
            for axis in 0..3 {
                if fixed.components[axis] {
                    v[3 * fixed.atom + axis] = 0.;
                }
            }
        }
    }
}
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Frame {
    positions: Vec<f64>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    cell_columns_bohr: Option<[[f64; 3]; 3]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    stress_ha_bohr3: Option<[[f64; 3]; 3]>,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pulay_assessment: Option<Value>,
    fractional: Vec<[f64; 3]>,
    energy_ha: f64,
    forces_ha_bohr: Vec<[f64; 3]>,
    max_force_ha_bohr: f64,
    step_length_bohr: f64,
    scf_iterations: usize,
    provenance: String,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Restart {
    version: u32,
    identity: String,
    optimizer: optimizer::Optimizer,
    frames: Vec<Frame>,
    density: Vec<f64>,
    consecutive: usize,
    status: String,
    failures: Vec<Value>,
}
struct Context {
    config: Config,
    cell: Cell,
    opts: ScfOptions,
    potentials: Vec<(String, dft_upf::Pseudopotential)>,
    provenance: Value,
}
impl Context {
    fn position(&mut self, fractional: &[[f64; 3]]) {
        for (a, f) in self.config.system.atoms.iter_mut().zip(fractional) {
            a.fractional = *f;
        }
    }
    fn atoms(&self) -> Vec<Atom> {
        self.config
            .system
            .atoms
            .iter()
            .map(|a| Atom::new(&a.species, a.fractional))
            .collect()
    }
    fn refs(&self) -> Vec<SpeciesPotential<'_>> {
        self.potentials
            .iter()
            .map(|(label, potential)| SpeciesPotential { label, potential })
            .collect()
    }
    fn problem(&self) -> Result<ScfProblem, RunError> {
        let points = self.config.resolved_points().map_err(|e| error(1, e))?;
        let mut p = if self.config.basis.legacy_gamma() {
            ScfProblem::new(&self.cell, &self.atoms(), &self.refs(), self.opts.clone())
        } else {
            ScfProblem::with_kpoints(
                &self.cell,
                &self.atoms(),
                &self.refs(),
                self.opts.clone(),
                &points,
                dft_engine::OccupationPolicy::Insulator {
                    minimum_gap_hartree: 1e-7,
                },
            )
        }
        .map_err(|e| error(4, e))?;
        if self.config.symmetry.mode == crate::config::SymmetryMode::Auto {
            p.enable_symmetry(
                &self.config.symmetry_full_mesh().map_err(|e| error(1, e))?,
                self.config.symmetry.tolerance_bohr,
                self.config.basis.time_reversal(),
            )
            .map_err(|e| error(4, e))?;
        }
        if let Some(floor) = self.config.xc.density_floor {
            p.set_density_floor(floor).map_err(|e| error(4, e))?;
        }
        Ok(p)
    }
    fn evaluate(
        &self,
        positions: Vec<f64>,
        length: f64,
        warm: Option<(&ScfProblem, &[f64])>,
        settings: &Settings,
    ) -> Result<(Frame, Vec<f64>), RunError> {
        let problem = self.problem()?;
        let mut density = problem.uniform_density();
        let mut session =
            ScfSession::new(problem, CompletionMode::Strict, ExecutionConfig::default())
                .map_err(|e| error(4, e))?;
        if let Some((source, values)) = warm {
            if values.len() != density.len() {
                return Err(error(4, "warm_transfer: density size mismatch"));
            }
            for (d, v) in density.iter_mut().zip(values) {
                *d = *v;
            }
            if settings.mode == cell::Mode::Ions {
                session.initialize_density_transfer(source, &density)
            } else {
                session.initialize_cell_density_transfer(source, &density)
            }
            .map_err(|e| error(4, e))?;
        } else {
            session.initialize().map_err(|e| error(4, e))?;
        }
        while session.step().map_err(|e| error(4, e))?.status == StepStatus::Ready {}
        let result = session
            .into_result()
            .ok_or_else(|| error(4, "SCF unconverged; forces unavailable"))?;
        if !result.converged {
            return Err(error(4, "SCF unconverged; forces unavailable"));
        }
        let state = result.result;
        let ops = self
            .config
            .symmetry_map()
            .map_err(|e| error(1, e))?
            .map_or_else(Vec::new, |m| {
                m.operations
                    .into_iter()
                    .filter(|o| !o.antiunitary)
                    .collect()
            });
        let forces = dft_engine::forces(&dft_engine::ForceState {
            cell: &self.cell,
            atoms: &self.atoms(),
            potentials: &self.refs(),
            options: &self.opts,
            electronic: dft_engine::ElectronicState::Fixed(&state),
            density_floor: self.config.xc.density_floor,
            symmetry: &ops,
        })
        .map_err(|e| error(4, e))?;
        let mut projected: Vec<_> = forces.values.iter().flatten().copied().collect();
        settings.project(&mut projected);
        let mut pulay_assessment = None;
        let stress_ha_bohr3 = if settings.mode != cell::Mode::Ions {
            let atoms = self.atoms();
            let refs = self.refs();
            let view = dft_engine::ForceState {
                cell: &self.cell,
                atoms: &atoms,
                potentials: &refs,
                options: &self.opts,
                electronic: dft_engine::ElectronicState::Fixed(&state),
                density_floor: self.config.xc.density_floor,
                symmetry: &ops,
            };
            let stress = if self.opts.strain_reference.is_none() {
                dft_engine::stress(&view)
            } else {
                dft_engine::stress::frozen_stress(&view)
            }
            .map_err(|e| error(4, e))?;
            pulay_assessment = stress.pulay.map(|p| {
                json!({
                    "method":p.method, "correction_ha_bohr3":p.correction,
                    "uncertainty_ha_bohr3":p.uncertainty, "cutoff_steps_ha":p.cutoff_steps_ha,
                    "shells_resolved":p.shells_resolved, "warning":p.warning,
                })
            });
            Some(stress.symmetric)
        } else {
            None
        };
        let frame = Frame {
            cell_columns_bohr: (settings.mode != cell::Mode::Ions)
                .then(|| self.cell.lattice_columns_bohr()),
            stress_ha_bohr3,
            pulay_assessment,
            positions,
            fractional: self
                .config
                .system
                .atoms
                .iter()
                .map(|a| a.fractional)
                .collect(),
            energy_ha: state.total_hartree(),
            forces_ha_bohr: forces.values,
            max_force_ha_bohr: projected.iter().map(|v| v.abs()).fold(0., f64::max),
            step_length_bohr: length,
            scf_iterations: state.history.len(),
            provenance: if warm.is_some() {
                "warm_transfer"
            } else {
                "fresh"
            }
            .into(),
        };
        Ok((frame, state.density.iter().copied().collect()))
    }
}
// Decimal JSON parsing need not preserve binary64. Encode every floating state
// member as integer bits, independently of serde_json's float parsing features.
fn bits(value: &mut Value, encode: bool) -> Result<(), RunError> {
    match value {
        Value::Number(n) if encode && n.is_f64() => {
            *value = json!({"f64_bits":n.as_f64().unwrap().to_bits()});
        }
        Value::Object(m) if !encode && m.contains_key("f64_bits") => {
            if m.len() != 1 {
                return Err(error(1, "invalid restart float"));
            }
            let f = f64::from_bits(
                m["f64_bits"]
                    .as_u64()
                    .ok_or_else(|| error(1, "invalid float bits"))?,
            );
            if !f.is_finite() {
                return Err(error(1, "nonfinite restart float"));
            }
            *value = json!(f);
        }
        Value::Object(m) => {
            for v in m.values_mut() {
                bits(v, encode)?;
            }
        }
        Value::Array(a) => {
            for v in a {
                bits(v, encode)?;
            }
        }
        _ => {}
    }
    Ok(())
}
fn atomic(path: &Path, bytes: &[u8]) -> Result<(), RunError> {
    let tmp = path.with_extension(format!("tmp-{}", std::process::id()));
    let mut f = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&tmp)?;
    let write = (|| {
        f.write_all(bytes)?;
        f.sync_all()?;
        fs::rename(&tmp, path)?;
        fs::File::open(path.parent().unwrap())?.sync_all()
    })();
    if write.is_err() {
        let _ = fs::remove_file(&tmp);
    }
    write.map_err(Into::into)
}
fn save(path: &Path, restart: &Restart) -> Result<(), RunError> {
    let mut value = serde_json::to_value(restart).map_err(|e| error(5, e))?;
    bits(&mut value, true)?;
    let payload = serde_json::to_string(&value).map_err(|e| error(5, e))?;
    let envelope = json!({"sha256":sha256(payload.as_bytes()),"payload":payload});
    atomic(
        path,
        &serde_json::to_vec(&envelope).map_err(|e| error(5, e))?,
    )
}
fn load(path: &Path, identity: &str, n: usize) -> Result<Restart, RunError> {
    if fs::metadata(path)?.len() > 512 * 1024 * 1024 {
        return Err(error(1, "restart too large"));
    }
    let envelope: Value = serde_json::from_slice(&fs::read(path)?).map_err(|e| error(1, e))?;
    let payload = envelope["payload"]
        .as_str()
        .ok_or_else(|| error(1, "restart payload missing"))?;
    if envelope["sha256"].as_str() != Some(&sha256(payload.as_bytes())) {
        return Err(error(1, "restart checksum mismatch"));
    }
    let mut value: Value = serde_json::from_str(payload).map_err(|e| error(1, e))?;
    bits(&mut value, false)?;
    let state: Restart = serde_json::from_value(value).map_err(|e| error(1, e))?;
    if state.version != 1
        || state.identity != identity
        || state.frames.is_empty()
        || state.optimizer.h.len() != n * n
        || state.optimizer.velocity.len() != n
        || state.frames.iter().any(|f| {
            f.positions.len() != n || f.fractional.len() * 3 != n || f.forces_ha_bohr.len() * 3 != n
        })
        || (!state.optimizer.previous_x.is_empty()
            && (state.optimizer.previous_x.len() != n || state.optimizer.previous_force.len() != n))
    {
        return Err(error(
            1,
            "restart version, build/input identity or shape mismatch",
        ));
    }
    Ok(state)
}

/// Ionic iteration count excludes the initial evaluated structure. Exit 3 is
/// max-step exhaustion, exit 4 is an SCF/derivative failure, exit 0 converged.
pub fn execute(input: &Path, out: &Path, resume: bool) -> Result<(), RunError> {
    let bytes = fs::read(input)?;
    let config = Config::parse(std::str::from_utf8(&bytes).map_err(|e| error(1, e))?)
        .map_err(|e| error(1, e))?;
    if config.compat.reference() {
        return Err(error(1, crate::config::REFERENCE_SCOPE));
    }
    let settings = config.relax.clone().unwrap_or_default();
    settings.validate(config.system.atoms.len())?;
    let opts = options(&config)?;
    if opts.standard.is_some()
        || opts.energy_tolerance > 1e-8
        || opts.eigen_tolerance > 1e-8
        || opts.density_tolerance > 1e-7
    {
        return Err(error(
            1,
            "relax requires strict SCF force certificate tolerances; standard profile is unsupported",
        ));
    }
    if config
        .finite_temperature()
        .map_err(|e| error(1, e))?
        .is_some()
        || config.collinear()
    {
        return Err(error(
            1,
            "relax currently requires unpolarized fixed occupations",
        ));
    }
    let cell = config.cell().map_err(|e| error(1, e))?;
    // The Cartesian constraint projector must commute with every spatial
    // operation, including atom permutations. Otherwise auto symmetry cannot
    // preserve both the user's constraints and the selected space group.
    if let Some(map) = config.symmetry_map().map_err(|e| error(1, e))? {
        let mut free = vec![1.; config.system.atoms.len() * 3];
        settings.project(&mut free);
        for op in map.operations.iter().filter(|op| !op.antiunitary) {
            for (atom, target) in op.atom_permutation.iter().copied().enumerate() {
                for b in 0..3 {
                    let basis = cell.cartesian_bohr_to_fractional(std::array::from_fn(|i| {
                        if i == b { 1. } else { 0. }
                    }));
                    let rotated = cell.fractional_to_cartesian_bohr(std::array::from_fn(|a| {
                        (0..3).map(|i| op.rotation[a][i] as f64 * basis[i]).sum()
                    }));
                    for a in 0..3 {
                        if free[3 * target + a] != free[3 * atom + b] && rotated[a].abs() > 1e-12 {
                            return Err(error(
                                1,
                                "relax.fixed breaks the selected symmetry; use invariant constraints or symmetry off",
                            ));
                        }
                    }
                }
            }
        }
    }
    let mut identity_bytes = bytes.clone();
    let build_sha256 = sha256(&fs::read(std::env::current_exe()?)?);
    identity_bytes.extend(build_sha256.as_bytes());
    identity_bytes.extend(std::env::consts::ARCH.as_bytes());
    identity_bytes.extend(std::env::consts::OS.as_bytes());
    let mut potentials = Vec::new();
    let mut pp_hashes = serde_json::Map::new();
    let mut inputs = vec![fs::canonicalize(input)?];
    for (label, path) in &config.pseudopotentials {
        let path = input.parent().unwrap_or(Path::new(".")).join(path);
        let pp_bytes = fs::read(&path)?;
        let digest = sha256(&pp_bytes);
        identity_bytes.extend(digest.as_bytes());
        pp_hashes.insert(label.clone(), json!(digest));
        inputs.push(fs::canonicalize(path)?);
        let p = dft_upf::parse_upf(std::str::from_utf8(&pp_bytes).map_err(|e| error(2, e))?)
            .map_err(|e| error(2, e))?;
        potentials.push((label.clone(), p));
    }
    let identity = sha256(&identity_bytes);
    fs::create_dir_all(out)?;
    let names = [
        "restart.json",
        "trajectory.extxyz",
        "final.xyz",
        "final.cif",
        "POSCAR",
        "final.toml",
        "results.json",
        "relax.log",
    ];
    for name in names {
        let path = out.join(name);
        if let Ok(meta) = fs::symlink_metadata(&path) {
            if !meta.is_file()
                || meta.file_type().is_symlink()
                || inputs.contains(&fs::canonicalize(&path)?)
            {
                return Err(error(5, "refusing nonregular output or input overwrite"));
            }
            #[cfg(unix)]
            {
                use std::os::unix::fs::MetadataExt;
                if meta.nlink() > 1 {
                    return Err(error(5, "refusing hard-linked output"));
                }
            }
            if !resume {
                return Err(error(
                    5,
                    "relax outputs already exist; use a new directory or --resume",
                ));
            }
        }
    }
    let mut ctx = Context {
        config,
        cell,
        opts,
        potentials,
        provenance: json!({"input_sha256":sha256(&bytes),"executable_sha256":build_sha256,"pseudopotentials_sha256":pp_hashes,"platform":std::env::consts::OS,"architecture":std::env::consts::ARCH}),
    };
    if settings.mode != cell::Mode::Ions {
        return cell::execute(ctx, settings, identity, out, resume);
    }
    let n = ctx.config.system.atoms.len() * 3;
    let restart_path = out.join("restart.json");
    let mut state = if resume {
        load(&restart_path, &identity, n)?
    } else {
        let positions: Vec<_> = ctx
            .config
            .system
            .atoms
            .iter()
            .flat_map(|a| ctx.cell.fractional_to_cartesian_bohr(a.fractional))
            .collect();
        let mut failures = Vec::new();
        let initial = ctx
            .evaluate(positions.clone(), 0., None, &settings)
            .or_else(|e| {
                failures.push(
                    json!({"step":0,"attempt":1,"positions_bohr":positions,"reason":e.to_string()}),
                );
                ctx.evaluate(positions.clone(), 0., None, &settings)
            });
        let (frame, density) = match initial {
            Ok(v) => v,
            Err(e) => {
                failures.push(
                    json!({"step":0,"attempt":2,"positions_bohr":positions,"reason":e.to_string()}),
                );
                let report = json!({"schema_version":1,"calculation_type":"relaxation","provenance":ctx.provenance,"relaxation":{"status":"failed","steps":0,"last_good_structure":null,"final_max_force_ha_bohr":null,"constraints":settings.fixed,"failures":failures}});
                atomic(
                    &out.join("results.json"),
                    &serde_json::to_vec_pretty(&report).map_err(|e| error(5, e))?,
                )?;
                return Err(e);
            }
        };
        Restart {
            version: 1,
            identity,
            optimizer: optimizer::Optimizer::new(n, settings.trust_radius_bohr),
            frames: vec![frame],
            density,
            consecutive: 0,
            status: "running".into(),
            failures,
        }
    };
    ctx.position(&state.frames.last().unwrap().fractional);
    // Publish accepted state before derived outputs; resume regenerates outputs.
    save(&restart_path, &state)?;
    outputs(&ctx, &settings, &state, out)?;
    while state.frames.len() - 1 < settings.max_steps && state.status == "running" {
        let last = state.frames.last().unwrap().clone();
        ctx.position(&last.fractional);
        let source = ctx.problem()?;
        let mut f: Vec<_> = last.forces_ha_bohr.iter().flatten().copied().collect();
        settings.project(&mut f);
        let mut proposed = state.optimizer.clone();
        let mut step = proposed.step(
            &last.positions,
            &f,
            settings.method == Method::Fire,
            settings.max_displacement_bohr,
        );
        settings.project(&mut step);
        let mut accepted = None;
        for attempt in 0..2 {
            let scale = if attempt == 0 { 1. } else { 0.5 };
            let positions: Vec<_> = last
                .positions
                .iter()
                .zip(&step)
                .map(|(x, s)| if *s == 0. { *x } else { x + scale * s })
                .collect();
            // Preserve every unchanged fractional coordinate exactly for frozen atoms.
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
            let length = step
                .chunks(3)
                .map(|s| s.iter().map(|v| v * v).sum::<f64>().sqrt() * scale)
                .fold(0., f64::max);
            let warm = (attempt == 0 && settings.warm_start)
                .then_some((&source, state.density.as_slice()));
            match ctx.evaluate(positions.clone(), length, warm, &settings) {
                Ok((frame, density)) => {
                    accepted = Some((frame, density, attempt));
                    break;
                }
                Err(e) => state.failures.push(json!({"step":state.frames.len(),"attempt":attempt+1,"positions_bohr":positions,"fractional":fractional,"reason":e.to_string()})),
            }
        }
        if let Some((frame, density, attempt)) = accepted {
            let delta = frame.energy_ha - last.energy_ha;
            if frame.max_force_ha_bohr < settings.force_tol_ha_bohr
                && delta.abs() < settings.energy_tol_ha
            {
                state.consecutive += 1;
            } else {
                state.consecutive = 0;
            }
            if attempt > 0 || delta > settings.energy_tol_ha {
                proposed.radius *= 0.5;
            } else {
                proposed.radius = (proposed.radius * 1.1).min(settings.trust_radius_bohr);
            }
            if attempt > 0 {
                proposed.velocity.fill(0.);
            }
            state.optimizer = proposed;
            state.frames.push(frame);
            state.density = density;
            if state.consecutive >= 2 {
                state.status = "converged".into();
            }
        } else {
            state.status = "failed".into();
            ctx.position(&last.fractional);
        }
        if state.status == "running" && state.frames.len() - 1 == settings.max_steps {
            state.status = "unconverged".into();
        }
        save(&restart_path, &state)?;
        outputs(&ctx, &settings, &state, out)?;
    }
    match state.status.as_str() {
        "converged" => Ok(()),
        "failed" => Err(error(
            4,
            "relax failed after one fresh, half-step retry; outputs retain last good structure",
        )),
        _ => Err(error(3, "relax unconverged: max_steps reached")),
    }
}

#[allow(clippy::needless_range_loop)] // Lattice vectors are columns, file vectors are rows.
fn outputs(
    ctx: &Context,
    settings: &Settings,
    state: &Restart,
    out: &Path,
) -> Result<(), RunError> {
    use dft_units::{Angstrom, Bohr, EnergyUnit, from_hartree};
    use std::fmt::Write as _;
    let ang = Angstrom::from(Bohr(1.)).0;
    let ev = from_hartree(1., EnergyUnit::ElectronVolt);
    let labels: Vec<_> = ctx
        .config
        .system
        .atoms
        .iter()
        .map(|a| a.species.as_str())
        .collect();
    let elements: Vec<_> = labels
        .iter()
        .map(|label| {
            ctx.potentials
                .iter()
                .find(|(l, _)| l == label)
                .unwrap()
                .1
                .element()
        })
        .collect();
    let lattice = ctx.cell.lattice_columns_bohr();
    let mut xyz = String::new();
    let mut log = String::from(
        "step energy_Ha max_force_Ha_Bohr step_length_Bohr scf_iterations provenance\n",
    );
    for (i, frame) in state.frames.iter().enumerate() {
        let lattice = frame.cell_columns_bohr.unwrap_or(lattice);
        writeln!(xyz, "{}", labels.len()).unwrap();
        write!(xyz, "Lattice=\"").unwrap();
        for j in 0..3 {
            for row in &lattice {
                write!(xyz, "{:.17e} ", row[j] * ang).unwrap();
            }
        }
        writeln!(
            xyz,
            "\" Properties=species:S:1:pos:R:3:forces:R:3 energy={:.17e} pbc=\"T T T\"",
            frame.energy_ha * ev
        )
        .unwrap();
        for ((element, x), f) in elements
            .iter()
            .zip(frame.positions.chunks(3))
            .zip(&frame.forces_ha_bohr)
        {
            writeln!(
                xyz,
                "{element} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}",
                x[0] * ang,
                x[1] * ang,
                x[2] * ang,
                f[0] * ev / ang,
                f[1] * ev / ang,
                f[2] * ev / ang
            )
            .unwrap();
        }
        writeln!(
            log,
            "{i} {:.17e} {:.17e} {:.17e} {} {}",
            frame.energy_ha,
            frame.max_force_ha_bohr,
            frame.step_length_bohr,
            frame.scf_iterations,
            frame.provenance
        )
        .unwrap();
    }
    let last = state.frames.last().unwrap();
    let mut final_xyz = format!(
        "{}\npositions Angstrom; cell and forces omitted\n",
        labels.len()
    );
    for (el, x) in elements.iter().zip(last.positions.chunks(3)) {
        writeln!(
            final_xyz,
            "{el} {:.17e} {:.17e} {:.17e}",
            x[0] * ang,
            x[1] * ang,
            x[2] * ang
        )
        .unwrap();
    }
    let mut snippet = String::from("# Replace system.atoms in the original fixed-cell input.\n");
    for (label, f) in labels.iter().zip(&last.fractional) {
        writeln!(
            snippet,
            "[[system.atoms]]\nspecies = {label:?}\nfractional = [{:.17e}, {:.17e}, {:.17e}]",
            f[0], f[1], f[2]
        )
        .unwrap();
    }
    let mut poscar = String::from("BAANDI fixed-cell relaxation; Angstrom\n1.0\n");
    for j in 0..3 {
        writeln!(
            poscar,
            "{:.17e} {:.17e} {:.17e}",
            lattice[0][j] * ang,
            lattice[1][j] * ang,
            lattice[2][j] * ang
        )
        .unwrap();
    }
    // POSCAR groups by species; trajectory retains input order.
    let mut groups = Vec::new();
    for label in &labels {
        if !groups.contains(label) {
            groups.push(*label);
        }
    }
    writeln!(
        poscar,
        "{}",
        groups
            .iter()
            .map(|l| ctx
                .potentials
                .iter()
                .find(|(s, _)| s == l)
                .unwrap()
                .1
                .element())
            .collect::<Vec<_>>()
            .join(" ")
    )
    .unwrap();
    writeln!(
        poscar,
        "{}",
        groups
            .iter()
            .map(|l| labels.iter().filter(|a| *a == l).count().to_string())
            .collect::<Vec<_>>()
            .join(" ")
    )
    .unwrap();
    poscar.push_str("Cartesian\n");
    for group in groups {
        for (label, x) in labels.iter().zip(last.positions.chunks(3)) {
            if *label == group {
                writeln!(
                    poscar,
                    "{:.17e} {:.17e} {:.17e}",
                    x[0] * ang,
                    x[1] * ang,
                    x[2] * ang
                )
                .unwrap();
            }
        }
    }
    let result = json!({"schema_version":1,"calculation_type":"relaxation","provenance":ctx.provenance,"relaxation":{"status":state.status,"steps":state.frames.len()-1,"final_max_force_ha_bohr":last.max_force_ha_bohr,
        "energy_history_ha":state.frames.iter().map(|f|f.energy_ha).collect::<Vec<_>>(),"constraints":settings.fixed,"frames":state.frames,"failures":state.failures,
        "restart_semantics":"same_build_ionic_continuation; each transferred SCF is warm_transfer, not exact SCF replay","settings":settings}});
    for (name, text) in [
        ("trajectory.extxyz", xyz),
        ("final.xyz", final_xyz),
        ("POSCAR", poscar),
        ("final.toml", snippet),
        ("relax.log", log),
        (
            "results.json",
            serde_json::to_string_pretty(&result).map_err(|e| error(5, e))?,
        ),
    ] {
        atomic(&out.join(name), text.as_bytes())?;
    }
    println!(
        "step {} E={:.12e} Ha maxF={:.6e} Ha/Bohr displacement={:.6e} Bohr SCF={} {}",
        state.frames.len() - 1,
        last.energy_ha,
        last.max_force_ha_bohr,
        last.step_length_bohr,
        last.scf_iterations,
        state.status
    );
    Ok(())
}
