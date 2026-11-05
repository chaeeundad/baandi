//! Owned Python boundary objects; no Python callbacks enter numerical kernels.
//! PyO3 owns generated FFI and reference lifetimes. No handwritten unsafe code.
//! Enable `python-extension` through maturin to build the Python module.
//! Interpreter-dependent tests run against the installed wheel, not a Rust harness.
#![cfg(all(feature = "python-extension", not(test)))]

use dft_engine::{
    CompletionMode, ExecutionConfig, Mixing, ScfOptions, ScfProblem, ScfSession, SpeciesPotential,
    StepStatus,
    input::{Atom, Cell, parse_upf},
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

/// Gamma-only LDA settings; omitted fields preserve engine defaults.
#[pyclass(name = "Options", get_all, set_all)]
#[derive(Clone)]
struct Options {
    cutoff_hartree: f64,
    bands: Option<usize>,
    extra_bands: usize,
    max_iterations: usize,
    density_tolerance: f64,
    energy_tolerance: f64,
    eigen_tolerance: f64,
    seed: u64,
    mixing_beta: f64,
    mixing_history: usize,
}
#[pymethods]
impl Options {
    #[new]
    fn new() -> Self {
        let o = ScfOptions::default();
        Self {
            cutoff_hartree: o.cutoff_hartree,
            bands: o.bands,
            extra_bands: o.extra_bands,
            max_iterations: o.max_iterations,
            density_tolerance: o.density_tolerance,
            energy_tolerance: o.energy_tolerance,
            eigen_tolerance: o.eigen_tolerance,
            seed: o.seed,
            mixing_beta: 0.3,
            mixing_history: 6,
        }
    }
}
impl Options {
    fn engine(&self) -> PyResult<ScfOptions> {
        if self.max_iterations == 0
            || [
                self.cutoff_hartree,
                self.density_tolerance,
                self.energy_tolerance,
                self.eigen_tolerance,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x <= 0.)
        {
            return Err(value_error(
                "cutoff, tolerances and max_iterations must be positive and finite",
            ));
        }
        Ok(ScfOptions {
            cutoff_hartree: self.cutoff_hartree,
            bands: self.bands,
            extra_bands: self.extra_bands,
            max_iterations: self.max_iterations,
            density_tolerance: self.density_tolerance,
            energy_tolerance: self.energy_tolerance,
            eigen_tolerance: self.eigen_tolerance,
            seed: self.seed,
            mixing: Mixing::Pulay {
                beta: self.mixing_beta,
                history: self.mixing_history,
            },
            ..ScfOptions::default()
        })
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

// The engine's mixer trait is not Send. Construct and retain the entire session
// on one worker; only owned boundary values cross channels. This also avoids
// unsafe Send implementations and releases the GIL while waiting for SCF work.
type Job = Box<dyn FnOnce(&mut ScfSession) + Send>;
#[pyclass(name = "Session")]
struct Session {
    cell_rows_bohr: [[f64; 3]; 3],
    species: Vec<String>,
    positions_fractional: Vec<[f64; 3]>,
    pseudopotentials: BTreeMap<String, String>,
    force_options: ScfOptions,
    jobs: mpsc::Sender<Job>,
    volume: f64,
}
impl Session {
    fn call<T: Send + 'static>(
        &self,
        py: Python<'_>,
        f: impl FnOnce(&mut ScfSession) -> PyResult<T> + Send + 'static,
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
        let options = options.cloned().unwrap_or_else(Options::new).engine()?;
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
                    let setup = || -> PyResult<(ScfSession, f64)> {
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
                            .map(|(s, text)| parse_upf(text).map(|p| (s, p)).map_err(value_error))
                            .collect::<PyResult<_>>()?;
                        let refs: Vec<_> = potentials
                            .iter()
                            .map(|(s, p)| SpeciesPotential {
                                label: s,
                                potential: p,
                            })
                            .collect();
                        let problem =
                            ScfProblem::new(&cell, &atoms, &refs, options).map_err(value_error)?;
                        let session = ScfSession::new(
                            problem,
                            CompletionMode::RetainLast,
                            ExecutionConfig::default(),
                        )
                        .map_err(value_error)?;
                        Ok((session, volume))
                    };
                    match setup() {
                        Ok((mut session, volume)) => {
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
                cell_rows_bohr,
                species: retained.0,
                positions_fractional: retained.1,
                pseudopotentials: retained.2,
                force_options: retained.3,
            })
        })
    }
    fn initialize(&mut self, py: Python<'_>) -> PyResult<()> {
        self.call(py, |s| s.initialize().map_err(runtime_error))
    }
    /// Perform exactly one complete engine iteration. Never changes solver policy.
    fn step(&mut self, py: Python<'_>) -> PyResult<String> {
        self.call(py, |s| {
            s.step()
                .map(|r| format!("{:?}", r.status))
                .map_err(runtime_error)
        })
    }
    fn run(&mut self, py: Python<'_>) -> PyResult<()> {
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
        self.call(py, |s| Ok(format!("{:?}", s.snapshot().status)))
    }
    #[getter]
    fn iteration(&self, py: Python<'_>) -> PyResult<usize> {
        self.call(py, |s| Ok(s.snapshot().iteration))
    }
    fn result(&self, py: Python<'_>) -> PyResult<ResultView> {
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
    /// Cartesian forces in input atom order, Ha/Bohr. Requires certification.
    fn forces(&self, py: Python<'_>) -> PyResult<Vec<[f64; 3]>> {
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
    /// Tensile-positive Cartesian stress, Ha/Bohr^3; frozen Miller set. Requires certification.
    fn stress(&self, py: Python<'_>) -> PyResult<[[f64; 3]; 3]> {
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
        if let Some((correction, uncertainty, true)) = pulay {
            py.import("warnings")?.call_method1("warn", (format!("Frozen-Miller stress is not cutoff converged: isotropic Pulay estimate {correction:e} +/- {uncertainty:e} Ha/Bohr^3 exceeds 1e-5 or sampled shells are unresolved; increase cutoff"),))?;
        }
        Ok(values)
    }
    /// (shape, minimum, maximum, mean [Bohr^-3], integrated electrons).
    fn density_summary(&self, py: Python<'_>) -> PyResult<(Vec<usize>, f64, f64, f64, f64)> {
        let volume = self.volume;
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
    Ok(())
}
