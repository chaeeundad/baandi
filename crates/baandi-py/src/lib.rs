//! Owned Python boundary objects; no Python callbacks enter numerical kernels.
//! PyO3 owns generated FFI and reference lifetimes. No handwritten unsafe code.
//! Enable `python-extension` through maturin to build the Python module.
//! Interpreter-dependent tests run against the installed wheel, not a Rust harness.
#![cfg(all(feature = "python-extension", not(test)))]

use dft_engine::{
    CompletionMode, ExecutionConfig, Mixing, ScfOptions, ScfProblem, ScfSession, SpeciesPotential,
    StepStatus,
    input::{Atom, Cell, kpoints, parse_upf},
    refscf::{api, profile::Profile, runtime, stepper::StepReport},
    refstart::StartInput,
};
use pyo3::{
    exceptions::{PyRuntimeError, PyValueError},
    prelude::*,
};
use std::{collections::BTreeMap, sync::mpsc};

fn value_error(e: impl std::fmt::Display) -> PyErr {
    PyValueError::new_err(e.to_string())
}
fn runtime_error(e: impl std::fmt::Display) -> PyErr {
    PyRuntimeError::new_err(e.to_string())
}

/// LDA settings. Reference runs use atomic starts, Broyden beta 0.7,
/// history 8 and fixed-point verification. Explicit legacy sessions retain
/// their former Pulay defaults and legacy step status strings.
#[pyclass(name = "Options")]
#[derive(Clone)]
struct Options {
    #[pyo3(get, set)]
    driver: String,
    /// Unshifted mesh, reduced by time reversal as in the CLI.
    #[pyo3(get, set)]
    kgrid: [usize; 3],
    #[pyo3(get, set)]
    reference_strict: bool,
    #[pyo3(get, set)]
    standard_threshold_ha: f64,
    #[pyo3(get, set)]
    verify: bool,
    #[pyo3(get, set)]
    cutoff_hartree: f64,
    #[pyo3(get, set)]
    bands: Option<usize>,
    #[pyo3(get, set)]
    extra_bands: usize,
    #[pyo3(get, set)]
    max_iterations: usize,
    #[pyo3(get, set)]
    density_tolerance: f64,
    #[pyo3(get, set)]
    energy_tolerance: f64,
    #[pyo3(get, set)]
    eigen_tolerance: f64,
    #[pyo3(get, set)]
    seed: u64,
    mixing_beta: Option<f64>,
    mixing_history: Option<usize>,
}
#[pymethods]
impl Options {
    #[new]
    fn new() -> Self {
        let o = ScfOptions::from_reference_switch();
        Self {
            driver: "reference".into(),
            kgrid: [1; 3],
            reference_strict: o.is_strict_reference(),
            standard_threshold_ha: dft_engine::StandardProfile::default().threshold_ha,
            verify: true,
            cutoff_hartree: o.cutoff_hartree,
            bands: o.bands,
            extra_bands: o.extra_bands,
            max_iterations: o.max_iterations,
            density_tolerance: o.density_tolerance,
            energy_tolerance: o.energy_tolerance,
            eigen_tolerance: o.eigen_tolerance,
            seed: o.seed,
            mixing_beta: None,
            mixing_history: None,
        }
    }
    #[getter(mixing_beta)]
    fn beta(&self) -> f64 {
        self.mixing_beta.unwrap_or(if self.driver == "reference" {
            0.7
        } else if self.reference_strict {
            0.3
        } else {
            0.5
        })
    }
    #[setter(mixing_beta)]
    fn set_beta(&mut self, beta: f64) {
        self.mixing_beta = Some(beta);
    }
    #[getter(mixing_history)]
    fn history(&self) -> usize {
        self.mixing_history
            .unwrap_or(if self.driver == "reference" { 8 } else { 6 })
    }
    #[setter(mixing_history)]
    fn set_history(&mut self, history: usize) {
        self.mixing_history = Some(history);
    }
}
impl Options {
    fn engine(&self) -> PyResult<ScfOptions> {
        if !matches!(self.driver.as_str(), "reference" | "legacy") {
            return Err(value_error("driver must be 'reference' or 'legacy'"));
        }
        if self.driver == "legacy" && self.kgrid != [1; 3] {
            return Err(value_error("non-Gamma kgrid requires driver='reference'"));
        }
        if self.driver == "reference" && self.reference_strict {
            return Err(value_error("reference_strict requires driver='legacy'"));
        }
        if self.driver == "reference" {
            self.profile()?;
        }
        if self.max_iterations == 0
            || [
                self.cutoff_hartree,
                self.density_tolerance,
                self.energy_tolerance,
                self.eigen_tolerance,
                self.standard_threshold_ha,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.)
        {
            return Err(value_error(
                "cutoff, tolerances, standard_threshold_ha and max_iterations must be positive and finite",
            ));
        }
        let base = if self.reference_strict {
            ScfOptions::strict_reference()
        } else {
            ScfOptions {
                standard: Some(dft_engine::StandardProfile {
                    threshold_ha: self.standard_threshold_ha,
                }),
                verification: if self.verify {
                    dft_engine::Verification::FixedPoint
                } else {
                    dft_engine::Verification::Off
                },
                ..ScfOptions::default()
            }
        };
        let mut options = ScfOptions {
            cutoff_hartree: self.cutoff_hartree,
            bands: self.bands,
            extra_bands: self.extra_bands,
            max_iterations: self.max_iterations,
            density_tolerance: self.density_tolerance,
            energy_tolerance: self.energy_tolerance,
            eigen_tolerance: self.eigen_tolerance,
            seed: self.seed,
            mixing: Mixing::Pulay {
                beta: self.beta(),
                history: self.history(),
            },
            ..base
        };
        if self.driver == "reference" {
            options.start_density = dft_engine::StartDensity::Atomic;
            options.start_wavefunctions = dft_engine::StartWavefunctions::Atomic;
            options.gamma_real = if self.kgrid == [1; 3] {
                dft_engine::GammaReal::On
            } else {
                dft_engine::GammaReal::Off
            };
        }
        Ok(options)
    }
    fn profile(&self) -> PyResult<Profile> {
        Profile {
            beta: self.beta(),
            history: self.history(),
            ..Profile::from_standard_ha(self.standard_threshold_ha).map_err(value_error)?
        }
        .validate()
        .map_err(value_error)
    }
}

/// Copied terminal output, independent of the session's lifetime.
#[pyclass(name = "Result", frozen, get_all)]
struct ResultView {
    converged: bool,
    energy_hartree: f64,
    eigenvalues_hartree: Vec<f64>,
    charge_electrons: f64,
    iterations: usize,
}

/// One reference map, including verification and retry maps.
#[pyclass(name = "StepReport", frozen, get_all)]
struct StepView {
    iteration: usize,
    kind: String,
    stage: String,
    estimator_ha: f64,
    dr2_ry: f64,
    energy_ha: Option<f64>,
    ethr_ry: f64,
    h_columns: usize,
    reused: bool,
    safe_point: bool,
}
impl From<StepReport> for StepView {
    fn from(r: StepReport) -> Self {
        Self {
            iteration: r.iteration,
            kind: format!("{:?}", r.kind),
            stage: r.stage.into(),
            estimator_ha: r.estimator_ha,
            dr2_ry: r.dr2_ry,
            energy_ha: r.energy_ha,
            ethr_ry: r.ethr_ry,
            h_columns: r.h_columns,
            reused: r.reused,
            safe_point: r.safe_point,
        }
    }
}

// The engine's mixer trait is not Send. Construct and retain the entire session
// on one worker; only owned boundary values cross channels. This also avoids
// unsafe Send implementations and releases the GIL while waiting for SCF work.
enum Worker<'a> {
    Legacy(Box<ScfSession>),
    Reference(Box<ReferenceSession<'a>>),
}

struct ReferenceSession<'a> {
    input: StartInput<'a>,
    runtime: Option<runtime::Runtime<'a, api::Delivered>>,
    snapshot: runtime::Snapshot,
    delivered: Option<api::Delivered>,
}
impl ReferenceSession<'_> {
    fn initialize(&mut self) -> PyResult<()> {
        let runtime = self
            .runtime
            .as_mut()
            .ok_or_else(|| runtime_error("already completed"))?;
        let result = runtime.initialize().map_err(runtime_error);
        self.snapshot = runtime.snapshot();
        result
    }
    fn step(&mut self) -> PyResult<StepView> {
        let runtime = self
            .runtime
            .as_mut()
            .ok_or_else(|| runtime_error("already completed"))?;
        let result = runtime.step().map_err(runtime_error);
        self.snapshot = runtime.snapshot();
        let report = result?;
        if runtime.result().is_some() {
            self.delivered = Some(
                self.runtime
                    .take()
                    .expect("live runtime")
                    .finish()
                    .map_err(runtime_error)?,
            );
        }
        Ok(report.into())
    }
    fn run(&mut self) -> PyResult<()> {
        use dft_engine::runtime::RuntimeState;
        if self.snapshot.state == RuntimeState::Created {
            self.initialize()?;
        }
        while self.snapshot.state == RuntimeState::Initialized {
            self.step()?;
        }
        if self.snapshot.state == RuntimeState::Failed {
            return Err(runtime_error("SCF failed"));
        }
        Ok(())
    }
    fn outcome(&self) -> PyResult<&dft_engine::ScfOutcome> {
        match &self
            .delivered
            .as_ref()
            .ok_or_else(|| runtime_error("terminal result not available; run first"))?
            .run
        {
            api::Run::Scalar(run) => Ok(&run.outcome),
            api::Run::Electronic(_) => Err(runtime_error("expected scalar reference result")),
        }
    }
    fn derivatives(&mut self, stress: bool) -> PyResult<&api::Observables> {
        let input = || StartInput {
            cell: self.input.cell,
            atoms: self.input.atoms,
            potentials: self.input.potentials,
            options: self.input.options,
        };
        let delivered = self
            .delivered
            .as_mut()
            .ok_or_else(|| runtime_error("derivatives require a converged result"))?;
        let api::Run::Scalar(run) = &mut delivered.run else {
            return Err(runtime_error("expected scalar reference result"));
        };
        let derivatives = &mut delivered.derivatives;
        derivatives
            .prepare(input(), &mut run.outcome)
            .map_err(runtime_error)?;
        derivatives
            .evaluate_forces(input(), &run.outcome.result)
            .map_err(runtime_error)?;
        if stress {
            derivatives
                .evaluate_stress(input(), &run.outcome.result)
                .map_err(runtime_error)?;
        }
        Ok(derivatives)
    }
}

type Job = Box<dyn FnOnce(&mut Worker<'_>) + Send>;
#[pyclass(name = "Session")]
struct Session {
    cell_rows_bohr: [[f64; 3]; 3],
    species: Vec<String>,
    positions_fractional: Vec<[f64; 3]>,
    pseudopotentials: BTreeMap<String, String>,
    force_options: ScfOptions,
    jobs: mpsc::Sender<Job>,
    volume: f64,
    reference: bool,
}
impl Session {
    fn worker_call<T: Send + 'static>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut Worker) -> PyResult<T> + Send + 'static,
    ) -> PyResult<T> {
        py.detach(|| {
            let (tx, rx) = mpsc::channel();
            self.jobs
                .send(Box::new(move |s| {
                    let _ = tx.send(f(s));
                }))
                .map_err(|_| runtime_error("session worker stopped"))?;
            rx.recv()
                .map_err(|_| runtime_error("session worker failed"))?
        })
    }
    fn call<T: Send + 'static>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut ScfSession) -> PyResult<T> + Send + 'static,
    ) -> PyResult<T> {
        self.worker_call(py, move |s| match s {
            Worker::Legacy(s) => f(s),
            Worker::Reference(_) => unreachable!("driver fixed at construction"),
        })
    }
    fn reference_call<T: Send + 'static>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut ReferenceSession) -> PyResult<T> + Send + 'static,
    ) -> PyResult<T> {
        self.worker_call(py, move |s| match s {
            Worker::Reference(s) => f(s),
            Worker::Legacy(_) => unreachable!("driver fixed at construction"),
        })
    }
}
#[pymethods]
impl Session {
    #[new]
    #[pyo3(signature = (cell_rows_bohr, species, positions_fractional, pseudopotentials, options=None))]
    fn new(
        py: Python<'_>,
        cell_rows_bohr: [[f64; 3]; 3],
        species: Vec<String>,
        positions_fractional: Vec<[f64; 3]>,
        pseudopotentials: BTreeMap<String, String>,
        options: Option<&Options>,
    ) -> PyResult<Self> {
        if species.len() != positions_fractional.len() || species.is_empty() {
            return Err(value_error(
                "species and positions must have the same nonzero length",
            ));
        }
        let settings = options.cloned().unwrap_or_else(Options::new);
        let reference = settings.driver == "reference";
        let options = settings.engine()?;
        let points = kpoints::reduce_time_reversal(
            &kpoints::monkhorst_pack(settings.kgrid, [0.; 3]).map_err(value_error)?,
        )
        .map_err(value_error)?;
        let profile = if reference {
            Some(settings.profile()?)
        } else {
            None
        };
        let retained = (
            species.clone(),
            positions_fractional.clone(),
            pseudopotentials.clone(),
            options.clone(),
        );
        py.detach(move || {
            let (jobs, requests) = mpsc::channel::<Job>();
            let (ready, response) = mpsc::channel();
            std::thread::Builder::new()
                .name("baandi-session".into())
                .spawn(move || {
                    let setup = || -> PyResult<_> {
                        let cell =
                            Cell::from_user_rows_bohr(cell_rows_bohr).map_err(value_error)?;
                        let volume = cell.volume_bohr3();
                        let atoms: Vec<_> = species
                            .into_iter()
                            .zip(positions_fractional)
                            .map(|(s, p)| Atom::new(s, p))
                            .collect();
                        let potentials: Vec<_> = pseudopotentials
                            .iter()
                            .map(|(s, text)| {
                                parse_upf(text).map(|p| (s.clone(), p)).map_err(value_error)
                            })
                            .collect::<PyResult<_>>()?;
                        Ok((cell, volume, atoms, potentials))
                    };
                    let (cell, volume, atoms, potentials) = match setup() {
                        Ok(inputs) => inputs,
                        Err(e) => {
                            let _ = ready.send(Err(e));
                            return;
                        }
                    };
                    let refs: Vec<_> = potentials
                        .iter()
                        .map(|(s, p)| SpeciesPotential {
                            label: s,
                            potential: p,
                        })
                        .collect();
                    let build = || -> PyResult<Worker<'_>> {
                        if let Some(profile) = profile {
                            let input = || StartInput {
                                cell: &cell,
                                atoms: &atoms,
                                potentials: &refs,
                                options: &options,
                            };
                            let mut request = api::Request::new(profile, points);
                            request.derivatives = api::Derivatives {
                                forces: true,
                                stress: true,
                                deferred: true,
                            };
                            let runtime = runtime::Runtime::new(api::stepper(input(), request));
                            let snapshot = runtime.snapshot();
                            return Ok(Worker::Reference(Box::new(ReferenceSession {
                                input: input(),
                                runtime: Some(runtime),
                                snapshot,
                                delivered: None,
                            })));
                        }
                        let problem = ScfProblem::new(&cell, &atoms, &refs, options.clone())
                            .map_err(value_error)?;
                        let session = ScfSession::new(
                            problem,
                            CompletionMode::RetainLast,
                            ExecutionConfig::default(),
                        )
                        .map_err(value_error)?;
                        Ok(Worker::Legacy(Box::new(session)))
                    };
                    match build() {
                        Ok(mut session) => {
                            if ready.send(Ok(volume)).is_ok() {
                                for job in requests {
                                    job(&mut session);
                                }
                            }
                        }
                        Err(e) => {
                            let _ = ready.send(Err(e));
                        }
                    }
                })
                .map_err(runtime_error)?;
            let volume = response
                .recv()
                .map_err(|_| runtime_error("session setup worker failed"))??;
            Ok(Self {
                jobs,
                volume,
                reference,
                cell_rows_bohr,
                species: retained.0,
                positions_fractional: retained.1,
                pseudopotentials: retained.2,
                force_options: retained.3,
            })
        })
    }
    fn initialize(&mut self, py: Python<'_>) -> PyResult<()> {
        if self.reference {
            return self.reference_call(py, |s| s.initialize());
        }
        self.call(py, |s| s.initialize().map_err(runtime_error))
    }
    /// Perform one reference map (ordinary, verification or retry); legacy
    /// sessions perform one legacy iteration. Never changes solver policy.
    fn step(&mut self, py: Python<'_>) -> PyResult<Py<PyAny>> {
        if self.reference {
            return Ok(Py::new(py, self.reference_call(py, |s| s.step())?)?.into_any());
        }
        let status = self.call(py, |s| {
            s.step()
                .map(|r| format!("{:?}", r.status))
                .map_err(runtime_error)
        })?;
        Ok(status.into_pyobject(py)?.unbind().into_any())
    }
    fn run(&mut self, py: Python<'_>) -> PyResult<()> {
        if self.reference {
            return self.reference_call(py, |s| s.run());
        }
        self.call(py, |s| {
            if s.snapshot().status == StepStatus::Created {
                s.initialize().map_err(runtime_error)?;
            }
            while s.snapshot().status == StepStatus::Ready {
                s.step().map_err(runtime_error)?;
            }
            if s.snapshot().status == StepStatus::Failed {
                return Err(runtime_error("SCF failed"));
            }
            Ok(())
        })
    }
    #[getter]
    fn status(&self, py: Python<'_>) -> PyResult<String> {
        if self.reference {
            return self.reference_call(py, |s| Ok(format!("{:?}", s.snapshot.state)));
        }
        self.call(py, |s| Ok(format!("{:?}", s.snapshot().status)))
    }
    #[getter]
    fn iteration(&self, py: Python<'_>) -> PyResult<usize> {
        if self.reference {
            return self
                .reference_call(py, |s| Ok(s.snapshot.last_step.map_or(0, |r| r.iteration)));
        }
        self.call(py, |s| Ok(s.snapshot().iteration))
    }
    fn result(&self, py: Python<'_>) -> PyResult<ResultView> {
        if self.reference {
            return self.reference_call(py, |s| {
                let o = s.outcome()?;
                Ok(ResultView {
                    converged: o.converged,
                    energy_hartree: o.result.total_hartree(),
                    eigenvalues_hartree: o.result.eigenvalues.clone(),
                    charge_electrons: o.result.charge,
                    iterations: o.result.history.len(),
                })
            });
        }
        self.call(py, |s| {
            let o = s.result().ok_or_else(|| {
                runtime_error("terminal result not available; initialize and step or run first")
            })?;
            Ok(ResultView {
                converged: o.converged,
                energy_hartree: o.result.total_hartree(),
                eigenvalues_hartree: o.result.eigenvalues.clone(),
                charge_electrons: o.result.charge,
                iterations: s.snapshot().iteration,
            })
        })
    }
    /// Cartesian forces in input atom order, Ha/Bohr. Reference delivery uses
    /// the facade verification contract, including explicit verify=False opt-out.
    /// Legacy delivery requires its strict certificate or passed verification.
    fn forces(&self, py: Python<'_>) -> PyResult<Vec<[f64; 3]>> {
        if self.reference {
            return self.reference_call(py, |s| {
                Ok(s.derivatives(false)?
                    .forces
                    .as_ref()
                    .expect("requested forces")
                    .values
                    .clone())
            });
        }
        let rows = self.cell_rows_bohr;
        let species = self.species.clone();
        let positions = self.positions_fractional.clone();
        let potentials = self.pseudopotentials.clone();
        let options = self.force_options.clone();
        self.call(py, move |s| {
            let outcome = s
                .result()
                .ok_or_else(|| runtime_error("forces require a converged result"))?;
            let cell = Cell::from_user_rows_bohr(rows).map_err(value_error)?;
            let atoms: Vec<_> = species
                .into_iter()
                .zip(positions)
                .map(|(s, p)| Atom::new(s, p))
                .collect();
            let parsed: Vec<_> = potentials
                .iter()
                .map(|(s, t)| parse_upf(t).map(|p| (s, p)).map_err(value_error))
                .collect::<PyResult<_>>()?;
            let refs: Vec<_> = parsed
                .iter()
                .map(|(s, p)| SpeciesPotential {
                    label: s,
                    potential: p,
                })
                .collect();
            dft_engine::forces::forces(&dft_engine::forces::ForceState {
                cell: &cell,
                atoms: &atoms,
                potentials: &refs,
                options: &options,
                electronic: dft_engine::forces::ElectronicState::Fixed(&outcome.result),
                density_floor: None,
                symmetry: &[],
            })
            .map(|r| r.values)
            .map_err(runtime_error)
        })
    }
    /// Tensile-positive Cartesian stress, Ha/Bohr^3; frozen Miller set. Uses
    /// the same delivery contract as forces, at a threshold of at most
    /// `StandardProfile::STRESS_THRESHOLD_HA` (5e-11 Ha, the stress class; the
    /// CLI default for stress inputs, ADR 0043 ONEPATH-2/4).
    fn stress(&self, py: Python<'_>) -> PyResult<[[f64; 3]; 3]> {
        if let Some(p) = self.force_options.standard
            && p.threshold_ha > dft_engine::StandardProfile::STRESS_THRESHOLD_HA
        {
            return Err(value_error(format!(
                "stress needs Options.standard_threshold_ha <= {:e} Ha (the stress class, ADR 0043 ONEPATH-2); this session ran at {:e} Ha; set it before creating the session",
                dft_engine::StandardProfile::STRESS_THRESHOLD_HA,
                p.threshold_ha
            )));
        }
        if self.reference {
            let (values, pulay) = self.reference_call(py, |s| {
                let r = s
                    .derivatives(true)?
                    .stress
                    .as_ref()
                    .expect("requested stress");
                Ok((
                    r.values,
                    r.pulay
                        .as_ref()
                        .map(|p| (p.correction, p.uncertainty, p.warning)),
                ))
            })?;
            warn_pulay(py, pulay)?;
            return Ok(values);
        }
        let rows = self.cell_rows_bohr;
        let species = self.species.clone();
        let positions = self.positions_fractional.clone();
        let potentials = self.pseudopotentials.clone();
        let options = self.force_options.clone();
        let (values, pulay) = self.call(py, move |s| {
            let outcome = s
                .result()
                .ok_or_else(|| runtime_error("stress requires a converged result"))?;
            let cell = Cell::from_user_rows_bohr(rows).map_err(value_error)?;
            let atoms: Vec<_> = species
                .into_iter()
                .zip(positions)
                .map(|(s, p)| Atom::new(s, p))
                .collect();
            let parsed: Vec<_> = potentials
                .iter()
                .map(|(s, t)| parse_upf(t).map(|p| (s, p)).map_err(value_error))
                .collect::<PyResult<_>>()?;
            let refs: Vec<_> = parsed
                .iter()
                .map(|(s, p)| SpeciesPotential {
                    label: s,
                    potential: p,
                })
                .collect();
            dft_engine::stress::stress(&dft_engine::forces::ForceState {
                cell: &cell,
                atoms: &atoms,
                potentials: &refs,
                options: &options,
                electronic: dft_engine::forces::ElectronicState::Fixed(&outcome.result),
                density_floor: None,
                symmetry: &[],
            })
            .map(|r| {
                (
                    r.values,
                    r.pulay.map(|p| (p.correction, p.uncertainty, p.warning)),
                )
            })
            .map_err(runtime_error)
        })?;
        warn_pulay(py, pulay)?;
        Ok(values)
    }
    /// Owned density values in ndarray iteration order, Bohr^-3, with grid shape.
    /// Copying is explicit and outside SCF maps. Reference requires delivery.
    fn density(&self, py: Python<'_>) -> PyResult<(Vec<usize>, Vec<f64>)> {
        if self.reference {
            return self.reference_call(py, |s| {
                let d = &s.outcome()?.result.density;
                Ok((d.shape().to_vec(), d.iter().copied().collect()))
            });
        }
        self.call(py, |s| {
            let d = s
                .density()
                .ok_or_else(|| runtime_error("density unavailable before initialization"))?;
            Ok((d.shape().to_vec(), d.iter().copied().collect()))
        })
    }
    /// (shape, minimum, maximum, mean [Bohr^-3], integrated electrons).
    fn density_summary(&self, py: Python<'_>) -> PyResult<(Vec<usize>, f64, f64, f64, f64)> {
        let volume = self.volume;
        if self.reference {
            return self.reference_call(py, move |s| {
                let d = &s.outcome()?.result.density;
                let mean = d.iter().sum::<f64>() / d.len() as f64;
                Ok((
                    d.shape().to_vec(),
                    d.iter().copied().fold(f64::INFINITY, f64::min),
                    d.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                    mean,
                    mean * volume,
                ))
            });
        }
        self.call(py, move |s| {
            let d = s
                .density()
                .ok_or_else(|| runtime_error("density unavailable before initialization"))?;
            let mean = d.iter().sum::<f64>() / d.len() as f64;
            Ok((
                d.shape().to_vec(),
                d.iter().copied().fold(f64::INFINITY, f64::min),
                d.iter().copied().fold(f64::NEG_INFINITY, f64::max),
                mean,
                mean * volume,
            ))
        })
    }
}

fn warn_pulay(py: Python<'_>, pulay: Option<(f64, f64, bool)>) -> PyResult<()> {
    if let Some((correction, uncertainty, true)) = pulay {
        py.import("warnings")?.call_method1("warn", (format!("Frozen-Miller stress is not cutoff converged: isotropic Pulay estimate {correction:e} +/- {uncertainty:e} Ha/Bohr^3 exceeds 1e-5 or sampled shells are unresolved; increase cutoff"),))?;
    }
    Ok(())
}

/// Use precisely the CLI's conversion operation, not ASE's CODATA constants.
#[pyfunction]
fn hartree_to_ev(value: f64) -> f64 {
    dft_units::from_hartree(value, dft_units::EnergyUnit::ElectronVolt)
}

#[pyfunction]
fn angstrom_to_bohr(value: f64) -> f64 {
    dft_units::Bohr::from(dft_units::Angstrom(value)).0
}

#[pymodule]
fn _baandi(m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add_function(wrap_pyfunction!(hartree_to_ev, m)?)?;
    m.add_function(wrap_pyfunction!(angstrom_to_bohr, m)?)?;
    m.add_class::<Options>()?;
    m.add_class::<Session>()?;
    m.add_class::<ResultView>()?;
    m.add_class::<StepView>()?;
    Ok(())
}
