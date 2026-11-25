//! Reference SCF facade. Adapters own input parsing and output formatting; this
//! module owns dispatch, derivative delivery and fixed-cell continuation.
use super::{delivered, driver, electronic, profile::Profile, symmetry};
use crate::{model::ScfError, model::ScfOutcome, model::ScfResult, refstart::StartInput};
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
    /// Export continuation at the adapter's delivery boundary, using the same
    /// checked symmetry as derivatives rather than discarding model identity.
    pub fn ionic_seed(
        &self,
        input: StartInput<'_>,
        result: &ScfResult,
    ) -> Result<IonicSeed, ScfError> {
        if let Some(prepared) = &self.symmetry {
            let density: Vec<_> = result.density.iter().copied().collect();
            delivered::ionic_seed_prepared(input, result, &density, prepared)
        } else {
            ionic_seed(input, result)
        }
    }

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
    stepper(input, request).run_to_end()
}

/// Incremental facade with the same dispatch and derivative delivery as solve.
pub fn stepper(input: StartInput<'_>, request: Request) -> super::stepper::Stepper<'_, Delivered> {
    let comm = crate::parallel::SharedComm::active();
    super::stepper::Stepper::new(move |boundary| async move {
        solve_continuation(input, request, comm, &boundary).await
    })
}
async fn solve_continuation(
    input: StartInput<'_>,
    mut request: Request,
    comm: Option<crate::parallel::SharedComm>,
    boundary: &super::stepper::Boundary,
) -> Result<Delivered, ScfError> {
    if request.seed.is_some() && request.policy.is_some() {
        return Err(ScfError::InvalidInput(
            "reference API: ionic seeds with electronic policy are unsupported",
        ));
    }
    if request.policy.is_some()
        && (request.derivatives.forces || request.derivatives.stress || request.save_ionic_seed)
    {
        return Err(ScfError::InvalidInput(
            "reference API: spin/FD derivatives and ionic seeds are unsupported",
        ));
    }
    if request.save_ionic_seed && request.derivatives.deferred {
        return Err(ScfError::InvalidInput(
            "reference API: ionic seed export requires immediate delivery",
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
    let run = if let Some(comm) = &comm {
        run_request(model(), &mut request, comm, Some(comm.clone()), boundary).await?
    } else {
        run_request(
            model(),
            &mut request,
            &super::metric::Serial,
            None,
            boundary,
        )
        .await?
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
            output.ionic_seed = Some(output.derivatives.ionic_seed(model(), result)?);
        }
    }
    Ok(output)
}

async fn run_request<K: dft_mpi::Collective>(
    input: StartInput<'_>,
    request: &mut Request,
    comm: &K,
    owned_comm: Option<crate::parallel::SharedComm>,
    boundary: &super::stepper::Boundary,
) -> Result<Run, ScfError> {
    let points = request
        .symmetry
        .as_ref()
        .map_or(request.points.as_slice(), |s| s.points());
    if let Some(policy) = request.policy.take() {
        return Ok(Run::Electronic(Box::new(
            electronic::run_continuation(
                input,
                request.profile,
                points,
                policy,
                comm,
                request.symmetry.as_ref(),
                None,
                boundary,
            )
            .await?,
        )));
    }
    let run = match request.seed.take() {
        Some(IonicSeed::Gamma(seed)) => {
            if request.symmetry.is_some()
                || points.len() != 1
                || points[0].fractional != [0.; 3]
                || input.options.gamma_real != crate::GammaReal::On
            {
                return Err(ScfError::InvalidInput(
                    "reference ionic seed: Gamma/k mismatch",
                ));
            }
            driver::run_seeded_continuation(input, request.profile, comm, *seed, boundary).await?
        }
        Some(IonicSeed::K(seed)) => {
            super::k_driver::run_continuation(
                input,
                request.profile,
                points,
                comm,
                request.symmetry.as_ref(),
                Some(seed),
                boundary,
            )
            .await?
        }
        None if request.symmetry.is_none()
            && points.len() == 1
            && points[0].fractional == [0.; 3]
            && input.options.gamma_real == crate::GammaReal::On =>
        {
            dft_geometry::kpoints::validate(points).map_err(super::numerical)?;
            driver::run_inner(input, request.profile, comm, owned_comm, boundary).await?
        }
        None => {
            super::k_driver::run_continuation(
                input,
                request.profile,
                points,
                comm,
                request.symmetry.as_ref(),
                None,
                boundary,
            )
            .await?
        }
    };
    Ok(Run::Scalar(Box::new(run)))
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

#[cfg(test)]
mod seeded_symmetry_tests {
    use super::*;
    use crate::{SpeciesPotential, model::ScfOptions};
    use dft_geometry::{Atom, Cell};

    #[test]
    fn prepared_seed_preserves_stepper_delivery() {
        let _scope = dft_numerics::switches::RunConfig::from_pairs(&[])
            .unwrap()
            .enter();
        let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
        let target =
            Cell::from_user_rows_bohr([[5.02, 0., 0.], [0., 5.02, 0.], [0., 0., 5.02]]).unwrap();
        let atoms = [Atom::new("Xx", [0.; 3])];
        let pp =
            crate::scf::pool_fixture::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
        let potentials = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let points = vec![KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }];
        let mut options = ScfOptions {
            cutoff_hartree: 1.,
            bands: Some(2),
            max_iterations: 100,
            strain_reference: Some(Box::new(cell.clone())),
            ..Default::default()
        };
        let prepared = symmetry::Prepared::new(&cell, &atoms, &points, 1e-7, false).unwrap();
        options.grid_override = Some(
            prepared
                .grid(&StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                })
                .unwrap(),
        );
        let input = |cell| StartInput {
            cell,
            atoms: &atoms,
            potentials: &potentials,
            options: &options,
        };
        let request = |cell, seed| {
            let mut request = Request::new(Profile::default(), points.clone());
            request.symmetry =
                Some(symmetry::Prepared::new(cell, &atoms, &points, 1e-7, false).unwrap());
            request.seed = seed;
            request.save_ionic_seed = true;
            request.derivatives.forces = true;
            request
        };
        let mut first = solve(input(&cell), request(&cell, None)).unwrap();
        let seed = first
            .ionic_seed
            .take()
            .unwrap()
            .in_cell(&cell, &target, &options)
            .unwrap();
        let a = solve(input(&target), request(&target, Some(seed.clone()))).unwrap();
        let mut stepping = stepper(input(&target), request(&target, Some(seed)));
        stepping.initialize().unwrap();
        let mut steps = 0;
        while !stepping.is_done() {
            stepping.step().unwrap();
            steps += 1;
        }
        assert!(steps > 0);
        let b = stepping.finish().unwrap();
        assert_eq!(
            a.derivatives.forces.as_ref().unwrap().values,
            b.derivatives.forces.as_ref().unwrap().values
        );
        assert!(a.ionic_seed.is_some() && b.ionic_seed.is_some());
        let a = a.run.into_scalar().unwrap();
        let b = b.run.into_scalar().unwrap();
        assert!(a.outcome.converged && b.outcome.converged);
        assert_eq!(
            a.outcome.result.total_hartree().to_bits(),
            b.outcome.result.total_hartree().to_bits()
        );
        assert_eq!(a.outcome.result.density, b.outcome.result.density);
        for (a, b) in a
            .outcome
            .result
            .k_points
            .iter()
            .zip(&b.outcome.result.k_points)
        {
            assert_eq!(
                a.solution.eigenvectors.as_slice(),
                b.solution.eigenvectors.as_slice()
            );
        }
        // The unsymmetrized export must not silently discard this identity.
        assert!(ionic_seed(input(&target), &a.outcome.result).is_err());
        let density: Vec<_> = a.outcome.result.density.iter().copied().collect();
        assert!(
            delivered::ionic_seed_prepared(input(&target), &a.outcome.result, &density, &prepared)
                .is_err()
        );
    }
}
