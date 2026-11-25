//! Reference SCF facade. Adapters own input parsing and output formatting; this
//! module owns dispatch, derivative delivery and fixed-cell continuation.
use super::{delivered, driver, electronic, profile::Profile, symmetry};
use crate::{ScfError, ScfOutcome, ScfResult, refstart::StartInput};
use dft_geometry::kpoints::KPoint;

pub use delivered::IonicSeed;

/// Derivatives are evaluated forces first, then stress. Deferred delivery lets
/// an adapter read SCF diagnostics before stress's cutoff scans consume/reset
/// counters or issue collectives. It must call prepare, forces, then stress at
/// its existing output boundaries. No background work or extra SCF is implied.
#[derive(Default, Clone, Copy)]
pub struct Derivatives {
    pub forces: bool,
    pub stress: bool,
    pub deferred: bool,
}

pub struct Request {
    pub profile: Profile,
    pub points: Vec<KPoint>,
    pub symmetry: Option<symmetry::Prepared>,
    pub policy: Option<electronic::Policy>,
    pub seed: Option<IonicSeed>,
    pub derivatives: Derivatives,
    /// Export the delivered density/orbitals after derivatives for ionic reuse.
    pub save_ionic_seed: bool,
}
impl Request {
    pub fn new(profile: Profile, points: Vec<KPoint>) -> Self {
        Self {
            profile,
            points,
            symmetry: None,
            policy: None,
            seed: None,
            derivatives: Derivatives::default(),
            save_ionic_seed: false,
        }
    }
}

/// Keep each driver's complete records and native result, without converting a
/// spin/FD state to a scalar state or copying numerical arrays.
pub enum Run {
    Scalar(Box<driver::Run>),
    Electronic(Box<electronic::Run>),
}
impl Run {
    pub fn into_scalar(self) -> Result<driver::Run, ScfError> {
        match self {
            Self::Scalar(run) => Ok(*run),
            Self::Electronic(_) => {
                Err(ScfError::InvalidInput("reference API: expected scalar run"))
            }
        }
    }
    pub fn into_electronic(self) -> Result<electronic::Run, ScfError> {
        match self {
            Self::Electronic(run) => Ok(*run),
            Self::Scalar(_) => Err(ScfError::InvalidInput(
                "reference API: expected electronic run",
            )),
        }
    }
}

pub struct Delivered {
    pub run: Run,
    pub derivatives: Observables,
    pub ionic_seed: Option<IonicSeed>,
}

/// Optional results plus the checked symmetry retained for deferred delivery.
/// Fields contain results immediately for ordinary solve requests. For deferred
/// requests each evaluator fills its field once; the caller may then take it.
pub struct Observables {
    pub forces: Option<crate::forces::ForceResult>,
    pub stress: Option<crate::stress::StressResult>,
    requested: Derivatives,
    symmetry: Option<symmetry::Prepared>,
    prepared: bool,
    forces_done: bool,
    stress_done: bool,
}
impl Observables {
    pub fn prepare(
        &mut self,
        input: StartInput<'_>,
        outcome: &mut ScfOutcome,
    ) -> Result<(), ScfError> {
        if !self.prepared && (self.requested.forces || self.requested.stress) {
            if let Some(symmetry) = &self.symmetry {
                delivered::force_state_prepared(input, outcome, symmetry)?;
            } else {
                delivered::force_state(input, outcome)?;
            }
            self.prepared = true;
        }
        Ok(())
    }

    fn state<'a>(
        &'a self,
        input: StartInput<'a>,
        result: &'a ScfResult,
    ) -> Result<crate::ForceState<'a>, ScfError> {
        if !self.prepared {
            return Err(ScfError::InvalidInput(
                "reference API: prepare derivatives first",
            ));
        }
        Ok(crate::ForceState {
            cell: input.cell,
            atoms: input.atoms,
            potentials: input.potentials,
            options: input.options,
            electronic: crate::ElectronicState::Fixed(result),
            density_floor: None,
            symmetry: self
                .symmetry
                .as_ref()
                .map_or(&[], |s| s.operations(input.options.crystal_density_group)),
        })
    }

    pub fn evaluate_forces(
        &mut self,
        input: StartInput<'_>,
        result: &ScfResult,
    ) -> Result<(), ScfError> {
        if self.requested.forces && !self.forces_done {
            self.forces = Some(crate::forces(&self.state(input, result)?)?);
            self.forces_done = true;
        }
        Ok(())
    }

    /// Compatibility boundary for ionic adapters whose derivative operations
    /// are selected independently of SCF preparation. The force kernel checks
    /// that these operations match the delivered model; no identity is changed.
    pub fn evaluate_forces_with_operations(
        &mut self,
        input: StartInput<'_>,
        result: &ScfResult,
        operations: &[dft_geometry::symmetry::SymmetryOperation],
    ) -> Result<(), ScfError> {
        if self.requested.forces && !self.forces_done {
            let mut state = self.state(input, result)?;
            state.symmetry = operations;
            self.forces = Some(crate::forces(&state)?);
            self.forces_done = true;
        }
        Ok(())
    }

    pub fn evaluate_stress(
        &mut self,
        input: StartInput<'_>,
        result: &ScfResult,
    ) -> Result<(), ScfError> {
        if self.requested.stress && !self.stress_done {
            if self.requested.forces && !self.forces_done {
                return Err(ScfError::InvalidInput(
                    "reference API: evaluate forces before stress",
                ));
            }
            let state = self.state(input, result)?;
            self.stress = Some(if let Some(symmetry) = &self.symmetry {
                crate::stress::stress_prepared(&state, symmetry)?
            } else {
                crate::stress(&state)?
            });
            self.stress_done = true;
        }
        Ok(())
    }
}

/// Same input and request on every rank. Existing drivers select the active
/// communicator; this facade neither changes switches nor adds collectives.
pub fn solve(input: StartInput<'_>, request: Request) -> Result<Delivered, ScfError> {
    if request.seed.is_some() && (request.symmetry.is_some() || request.policy.is_some()) {
        return Err(ScfError::InvalidInput(
            "reference API: ionic seeds with symmetry or electronic policy are unsupported",
        ));
    }
    if request.policy.is_some()
        && (request.derivatives.forces || request.derivatives.stress || request.save_ionic_seed)
    {
        return Err(ScfError::InvalidInput(
            "reference API: spin/FD derivatives and ionic seeds are unsupported",
        ));
    }
    if request.save_ionic_seed && (request.symmetry.is_some() || request.derivatives.deferred) {
        return Err(ScfError::InvalidInput(
            "reference API: ionic seed export requires immediate delivery without symmetry",
        ));
    }
    // StartInput contains immutable references; reconstruct it without cloning
    // the model or introducing a new numerical lifetime/ownership convention.
    let model = || StartInput {
        cell: input.cell,
        atoms: input.atoms,
        potentials: input.potentials,
        options: input.options,
    };
    let run = if let Some(policy) = request.policy {
        Run::Electronic(Box::new(if let Some(symmetry) = &request.symmetry {
            electronic::run_prepared(model(), request.profile, symmetry, policy)?
        } else {
            electronic::run(model(), request.profile, &request.points, policy)?
        }))
    } else {
        Run::Scalar(Box::new(if let Some(symmetry) = &request.symmetry {
            super::k_driver::run_prepared(model(), request.profile, symmetry)?
        } else {
            delivered::run(model(), request.profile, &request.points, request.seed)?
        }))
    };
    let mut output = Delivered {
        run,
        derivatives: Observables {
            forces: None,
            stress: None,
            requested: request.derivatives,
            symmetry: request.symmetry,
            prepared: false,
            forces_done: false,
            stress_done: false,
        },
        ionic_seed: None,
    };
    if let Run::Scalar(run) = &mut output.run {
        if !request.derivatives.deferred {
            output.derivatives.prepare(model(), &mut run.outcome)?;
            output
                .derivatives
                .evaluate_forces(model(), &run.outcome.result)?;
            output
                .derivatives
                .evaluate_stress(model(), &run.outcome.result)?;
        }
        if request.save_ionic_seed {
            let result = &mut run.outcome.result;
            result.materialize_orbitals()?;
            output.ionic_seed = Some(ionic_seed(model(), result)?);
        }
    }
    Ok(output)
}

/// Export fixed-cell continuation at an adapter's existing delivery boundary.
/// This reads replicated delivered arrays and retains the active communicator.
pub fn ionic_seed(input: StartInput<'_>, result: &ScfResult) -> Result<IonicSeed, ScfError> {
    let density: Vec<_> = result.density.iter().copied().collect();
    if let Some(comm) = crate::parallel::SharedComm::active() {
        delivered::ionic_seed(input, result, &density, &comm)
    } else {
        delivered::ionic_seed(input, result, &density, &super::metric::Serial)
    }
}
