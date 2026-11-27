//! Fixed-cell ionic driver. SCF kernels and derivative ownership stay in engine.
pub(crate) mod cell;
pub mod optimizer;
use crate::{
    config::Config,
    run::{RunError, sha256},
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
/// Force measure of the stop rule (ADR 0056). `Norm`: the largest Cartesian
/// force norm |F_i| over atoms after constraint projection (ASE `fmax`, VASP
/// EDIFFG < 0). `Component`: the largest absolute free Cartesian component,
/// the rule before ADR 0056. Component <= norm <= sqrt(3) component, so at
/// equal tolerance the norm rule is the stricter one, by up to sqrt(3).
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "lowercase")]
pub enum ForceCriterion {
    #[default]
    Norm,
    Component,
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
/// C3 start-state predictor between ionic steps (ions mode only). The SCF
/// density is split as n = n_at(R) + d, with n_at the superposition of atomic
/// charges that moves rigidly with the atoms; only the bonding part d is
/// carried over (first order) or extrapolated along the ionic path (second
/// order). `none` keeps the plain previous-density warm transfer bit-identical.
#[derive(Clone, Copy, Debug, Default, Serialize, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "snake_case")]
pub enum Extrapolation {
    #[default]
    None,
    FirstOrder,
    SecondOrder,
}
impl Extrapolation {
    /// `BAANDI_RELAX_EXTRAPOLATION` (none, first_order, second_order) overrides
    /// the input file for A/B runs.
    pub(crate) fn effective(self) -> Result<Self, RunError> {
        match dft_numerics::switches::text(dft_numerics::switches::Id::RELAX_EXTRAPOLATION) {
            None => Ok(self),
            Some("none") => Ok(Self::None),
            Some("first_order") => Ok(Self::FirstOrder),
            Some("second_order") => Ok(Self::SecondOrder),
            Some(_) => Err(error(
                1,
                "BAANDI_RELAX_EXTRAPOLATION must be none, first_order or second_order",
            )),
        }
    }
}
/// One accepted ionic step of the predictor history: Cartesian positions and
/// the bonding part d = n_scf - n_at on the FFT grid.
#[derive(Clone, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct HistoryEntry {
    positions: Vec<f64>,
    bonding: Vec<f64>,
}
/// Predicted bonding density for the next step from up to three previous
/// steps (newest last). Second order fits the new displacement s = R_new - R_k
/// by a u + b v with u = R_k - R_{k-1}, v = R_{k-1} - R_{k-2} (least squares)
/// and applies the same combination to d. A Gram determinant below 1e-10 of
/// |u|^2 |v|^2 (nearly parallel steps) drops the v term. Returns the predicted
/// d and the start kind.
fn predict(
    history: &[HistoryEntry],
    new_positions: &[f64],
    mode: Extrapolation,
) -> Option<(Vec<f64>, &'static str)> {
    let k = history.last()?;
    if mode == Extrapolation::FirstOrder || history.len() == 1 {
        return Some((k.bonding.clone(), "extrapolated_first_order"));
    }
    let dot = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x * y).sum::<f64>();
    let diff = |a: &[f64], b: &[f64]| a.iter().zip(b).map(|(x, y)| x - y).collect::<Vec<_>>();
    let k1 = &history[history.len() - 2];
    let s = diff(new_positions, &k.positions);
    let u = diff(&k.positions, &k1.positions);
    let uu = dot(&u, &u);
    if uu == 0. {
        return Some((k.bonding.clone(), "extrapolated_first_order"));
    }
    let (a, b) = if history.len() >= 3 {
        let k2 = &history[history.len() - 3];
        let v = diff(&k1.positions, &k2.positions);
        let (vv, uv, su, sv) = (dot(&v, &v), dot(&u, &v), dot(&s, &u), dot(&s, &v));
        let det = uu * vv - uv * uv;
        if vv > 0. && det > 1e-10 * uu * vv {
            ((su * vv - sv * uv) / det, (sv * uu - su * uv) / det)
        } else {
            (su / uu, 0.)
        }
    } else {
        (dot(&s, &u) / uu, 0.)
    };
    let mut d = k.bonding.clone();
    for ((x, now), before) in d.iter_mut().zip(&k.bonding).zip(&k1.bonding) {
        *x += a * (now - before);
    }
    if b != 0. {
        let k2 = &history[history.len() - 3];
        for ((x, now), before) in d.iter_mut().zip(&k1.bonding).zip(&k2.bonding) {
            *x += b * (now - before);
        }
    }
    Some((d, "extrapolated_second_order"))
}
/// Append d = n_scf - n_at of an accepted step (at most three kept). Without
/// n_at (predictor off, or no atomic charge data) the history is cleared, so a
/// later prediction never mixes steps across a gap.
fn push_history(
    history: &mut Vec<HistoryEntry>,
    positions: &[f64],
    density: &[f64],
    n_at: Option<Vec<f64>>,
) {
    match n_at {
        Some(n_at) if n_at.len() == density.len() => {
            if history.len() == 3 {
                history.remove(0);
            }
            history.push(HistoryEntry {
                positions: positions.to_vec(),
                bonding: density.iter().zip(&n_at).map(|(n, a)| n - a).collect(),
            });
        }
        _ => history.clear(),
    }
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
    /// Force measure compared with `force_tol_ha_bohr` (default norm, ADR 0056).
    pub force_criterion: ForceCriterion,
    /// Optional energy-change requirement (ADR 0056). Unset (default) in ions
    /// mode: the force test alone stops. Set: |dE| < tol together with the force
    /// test on two consecutive accepted steps (the rule before ADR 0056). Cell
    /// modes keep that rule with 1e-8 Ha when unset ([`Settings::cell_energy_tol_ha`]).
    pub energy_tol_ha: Option<f64>,
    pub max_steps: usize,
    pub max_displacement_bohr: f64,
    pub trust_radius_bohr: f64,
    pub warm_start: bool,
    /// FIRE initial and maximum time step in ASE FIRE units, Angstrom/sqrt(eV)
    /// (unit masses, forces in eV/Angstrom): the same numbers move atoms as
    /// ase.optimize.FIRE does (ADR 0056). Defaults 0.1 and 1.0 (ASE defaults).
    pub fire_dt: f64,
    pub fire_dt_max: f64,
    /// C3 predictor; acts only with `warm_start` in ions mode. Default none.
    pub start_extrapolation: Extrapolation,
    pub fixed: Vec<Fixed>,
}
/// Default force tolerance (ADR 0056): 0.02 eV/Angstrom on the atomic norm,
/// 3.889e-4 Ha/Bohr. Without the energy test the force threshold alone sets
/// the accuracy; 0.02 eV/A is the PRELAX target and a common production value.
pub const DEFAULT_FORCE_TOL_EV_ANGSTROM: f64 = 0.02;
pub fn default_force_tol_ha_bohr() -> f64 {
    use dft_units::{Angstrom, Bohr, EnergyUnit, to_hartree};
    // eV/A -> Ha/Bohr: energy to Ha, then per Bohr = per A times (A per Bohr).
    to_hartree(DEFAULT_FORCE_TOL_EV_ANGSTROM, EnergyUnit::ElectronVolt) * Angstrom::from(Bohr(1.)).0
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
            force_tol_ha_bohr: default_force_tol_ha_bohr(),
            force_criterion: ForceCriterion::Norm,
            energy_tol_ha: None,
            max_steps: 200,
            max_displacement_bohr: 0.2,
            trust_radius_bohr: 0.2,
            warm_start: true,
            fire_dt: optimizer::FIRE_DT,
            fire_dt_max: optimizer::FIRE_DT_MAX,
            start_extrapolation: Extrapolation::None,
            fixed: vec![],
        }
    }
}
/// Per-step verification record (ONEPATH-1): threshold, outcome, attempts and
/// the last attempt's Coulomb fixed-point residual. None on the strict path.
/// Deterministic fields only (no wall time): frames are part of restart.json,
/// and a resumed relaxation must reproduce the uninterrupted files bit for bit.
fn verification_summary(state: &dft_engine::ScfResult, opts: &ScfOptions) -> Option<Value> {
    let p = opts.standard?;
    if opts.verification == dft_engine::Verification::Off {
        return Some(json!({
            "mode":"off", "status":"unverified", "passed":false,
            "threshold_ha":p.threshold_ha, "attempts":0,
            "coulomb_error_ha":null, "residual_bound_ha":null,
            "strict_criteria_met":null, "h_vectors":0,
        }));
    }
    let v = state.verification.as_ref()?;
    let last = v.attempts.last();
    Some(json!({
        "mode":"fixed_point","passed":v.passed(),"threshold_ha":p.threshold_ha,
        "attempts":v.attempts.len(),
        "coulomb_error_ha":last.map(|a| a.coulomb_error_ha),
        "residual_bound_ha":last.and_then(|a| a.residual_bound_ha),
        "strict_criteria_met":last.map(|a| a.strict_criteria_met),
        "h_vectors":v.attempts.iter().map(|a| a.operator_applications).sum::<usize>(),
    }))
}
/// Largest Euclidean norm of the per-atom triples of a flat Cartesian vector.
fn max_atom_norm(v: &[f64]) -> f64 {
    v.chunks(3)
        .map(|a| a.iter().map(|x| x * x).sum::<f64>().sqrt())
        .fold(0., f64::max)
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
            self.energy_tol_ha.unwrap_or(1.),
            self.max_displacement_bohr,
            self.trust_radius_bohr,
            self.fire_dt,
            self.fire_dt_max,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.)
            || self.fire_dt > self.fire_dt_max
            || !self.pressure_ha_bohr3.is_finite()
            || self.max_steps == 0
            || self.fixed.iter().any(|f| f.atom >= n)
        {
            return Err(error(
                1,
                "relax: positive finite tolerances/radii/FIRE time steps (fire_dt <= fire_dt_max), positive max_steps and valid zero-based fixed atom indices required",
            ));
        }
        Ok(())
    }
    /// Optimizer state with this input's FIRE time steps.
    fn optimizer(&self, n: usize) -> optimizer::Optimizer {
        optimizer::Optimizer::with_fire(n, self.trust_radius_bohr, self.fire_dt, self.fire_dt_max)
    }
    /// Force measure of `frame` used by the stop rule.
    fn force_measure(&self, frame: &Frame) -> f64 {
        match self.force_criterion {
            ForceCriterion::Norm => frame.max_force_norm_ha_bohr,
            ForceCriterion::Component => frame.max_force_ha_bohr,
        }
    }
    /// Force part of the stop rule (strict inequality, as ASE `fmax`).
    fn force_converged(&self, frame: &Frame) -> bool {
        self.force_measure(frame) < self.force_tol_ha_bohr
    }
    /// Energy tolerance of the cell modes: their stop rule and enthalpy
    /// backtracking keep the rule before ADR 0056 (default 1e-8 Ha).
    fn cell_energy_tol_ha(&self) -> f64 {
        self.energy_tol_ha.unwrap_or(1e-8)
    }
    /// Ions-mode stop rule after an accepted frame with energy change `delta`
    /// (None for the initial frame). Returns the new counter and whether the
    /// run converged. Without `energy_tol_ha` one frame that passes the force
    /// test stops (the initial one included); with it, the force test and
    /// |dE| < tol must hold on two consecutive accepted steps.
    fn ionic_stop(&self, frame: &Frame, delta: Option<f64>, consecutive: usize) -> (usize, bool) {
        let force = self.force_converged(frame);
        match self.energy_tol_ha {
            None => (usize::from(force), force),
            Some(tol) => {
                let next = if force && delta.is_some_and(|d| d.abs() < tol) {
                    consecutive + 1
                } else {
                    0
                };
                (next, next >= 2)
            }
        }
    }
    /// An energy rise above this halves the trust radius: the energy tolerance
    /// when given, otherwise twice the SCF energy accuracy (the standard
    /// threshold, or the strict energy tolerance), which equals the former
    /// default 1e-8 Ha at the default 5e-9 Ha threshold.
    fn energy_rise_tol_ha(&self, opts: &ScfOptions) -> f64 {
        self.energy_tol_ha.unwrap_or_else(|| {
            2. * opts
                .standard
                .map_or(opts.energy_tolerance, |p| p.threshold_ha)
        })
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
    /// Largest absolute free Cartesian component of the projected force.
    max_force_ha_bohr: f64,
    /// Largest atomic norm |F_i| of the projected force (ADR 0056).
    #[serde(default)]
    max_force_norm_ha_bohr: f64,
    step_length_bohr: f64,
    scf_iterations: usize,
    provenance: String,
    /// H vector applications of this step's SCF (solves plus certificates).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    scf_h_vectors: Option<usize>,
    /// Why a requested extrapolated start was replaced by the plain warm density.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    start_fallback: Option<String>,
    /// Charge fraction removed by clipping negative samples of an extrapolated
    /// start (then renormalised); absent when nothing was clipped.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    start_clipped_fraction: Option<f64>,
    /// ONEPATH-1: fixed-point verification of this step's SCF on the standard
    /// path (threshold, attempts, last Coulomb residual); absent on the strict
    /// path, whose forces rest on the strict certificate.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    verification: Option<Value>,
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
    /// C3 predictor history (at most three entries); empty when extrapolation is off.
    #[serde(default, skip_serializing_if = "Vec::is_empty")]
    extrapolation: Vec<HistoryEntry>,
}
/// Frame, converged density and (C3) the atomic superposition at its positions.
type Evaluated = (Frame, Vec<f64>, Option<Vec<f64>>);
struct Context {
    config: Config,
    cell: Cell,
    opts: ScfOptions,
    potentials: Vec<(String, dft_upf::Pseudopotential)>,
    provenance: Value,
    /// ADR 0043 path of every ionic step's SCF (recorded in results.json).
    scf_path: crate::config::ScfPath,
    selection: crate::driver_select::Selection,
    reference_seed: std::cell::RefCell<Option<(Cell, dft_engine::refscf::delivered::IonicSeed)>>,
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
        self.problem_with(false)
    }
    /// `atomic` also builds the atomic-charge superposition (C3 predictor). It
    /// sets only the cold start density, which a warm transfer replaces.
    fn problem_with(&self, atomic: bool) -> Result<ScfProblem, RunError> {
        let points = self.config.resolved_points().map_err(|e| error(1, e))?;
        let mut opts = self.opts.clone();
        if atomic {
            opts.start_density = dft_engine::StartDensity::Atomic;
        }
        if self.selection.reference() {
            // RELAX-CELL-DEFAULT: on the reference driver this problem is never
            // solved; it only supplies the FFT grid (frozen-Miller cell
            // segments) and the atomic superposition (C3 predictor). The strict
            // oracle constructor rejects the standard profile and verification
            // (R10BC), which do not affect either quantity, so build it with the
            // oracle's settings for those two fields.
            opts.standard = None;
            opts.verification = dft_engine::Verification::Off;
        }
        let mut p = if self.config.basis.legacy_gamma() {
            ScfProblem::new(&self.cell, &self.atoms(), &self.refs(), opts)
        } else {
            ScfProblem::with_kpoints(
                &self.cell,
                &self.atoms(),
                &self.refs(),
                opts,
                &points,
                dft_engine::OccupationPolicy::Insulator {
                    minimum_gap_hartree: 1e-7,
                },
            )
        }
        .map_err(|e| error(4, e))?;
        if self.config.symmetry.mode == crate::config::SymmetryMode::Auto {
            p.enable_symmetry_for_grid(
                &self.config.symmetry_full_mesh().map_err(|e| error(1, e))?,
                self.config.detection_tolerance(),
                self.config.basis.time_reversal(),
                self.config.basis.fft_grid,
            )
            .map_err(|e| error(4, e))?;
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
        self.evaluate_with(
            positions,
            length,
            warm.map(|(p, d)| (Some(p), d)),
            None,
            false,
            settings,
        )
        .map(|(frame, density, _)| (frame, density))
    }
    /// `predicted`: bonding density d and start kind for a C3 start
    /// n_at(R_new) + d (used only with a warm transfer); `atomic`: also return
    /// n_at at these positions for the predictor history.
    fn evaluate_with(
        &self,
        positions: Vec<f64>,
        length: f64,
        warm: Option<(Option<&ScfProblem>, &[f64])>,
        predicted: Option<(&[f64], &'static str)>,
        atomic: bool,
        settings: &Settings,
    ) -> Result<Evaluated, RunError> {
        let reference = self.selection.reference();
        let mut reference_derivatives = None;
        let (state, n_at, provenance, start_fallback, start_clipped_fraction) = if reference {
            let atoms = self.atoms();
            let refs = self.refs();
            let points = self.config.resolved_points().map_err(|e| error(1, e))?;
            // Reuse the existing atomic-superposition builder, without running
            // a legacy SCF. The predictor remains an opt-in ionic start policy.
            let n_at = if atomic {
                self.problem_with(true)?
                    .atomic_superposition_density()
                    .map(|a| a.iter().copied().collect::<Vec<_>>())
            } else {
                None
            };
            let mut seed = if warm.is_some() {
                self.reference_seed
                    .borrow()
                    .clone()
                    .map(|(source, seed)| {
                        if settings.mode == cell::Mode::Ions {
                            Ok(seed)
                        } else {
                            seed.in_cell(&source, &self.cell, &self.opts)
                        }
                    })
                    .transpose()
                    .map_err(|e| error(4, e))?
            } else {
                None
            };
            let mut provenance = if seed.is_some() {
                "warm_transfer"
            } else {
                "fresh"
            };
            let mut fallback = (warm.is_some() && seed.is_none())
                .then(|| "reference restart has no orbital seed; fresh atomic start".to_string());
            let mut clipped = None;
            if let Some((bonding, kind)) = predicted
                && let Some(current) = seed.take()
            {
                let mut replacement = None;
                match &n_at {
                    Some(n_at) if n_at.len() == bonding.len() => {
                        let mut start: Vec<_> =
                            n_at.iter().zip(bonding).map(|(a, d)| a + d).collect();
                        let total: f64 = start.iter().sum();
                        let negative: f64 = start.iter().filter(|v| **v < 0.).sum();
                        start.iter_mut().for_each(|v| *v = v.max(0.));
                        let kept: f64 = start.iter().sum();
                        if total.is_finite() && total > 0. && kept > 0. {
                            start.iter_mut().for_each(|v| *v *= total / kept);
                            replacement = Some(start);
                            provenance = kind;
                            if negative < 0. {
                                clipped = Some(-negative / total);
                            }
                        } else {
                            fallback = Some("extrapolated density not finite".into());
                        }
                    }
                    Some(_) => fallback = Some("predictor size mismatch".into()),
                    None => fallback = Some("no atomic charge density (PP_RHOATOM)".into()),
                }
                seed = Some(if let Some(density) = replacement {
                    let charge = warm.unwrap().1.iter().sum::<f64>() * self.cell.volume_bohr3()
                        / density.len() as f64;
                    current
                        .with_density(&density, self.cell.volume_bohr3(), charge)
                        .map_err(|e| error(4, e))?
                } else {
                    current
                });
            }
            let mut request = dft_engine::refscf::api::Request::new(
                crate::run::reference_profile(&self.config)?,
                points,
            );
            request.seed = seed;
            if self.config.symmetry.mode == crate::config::SymmetryMode::Auto {
                request.symmetry = Some(
                    dft_engine::refscf::symmetry::Prepared::new_for_grid(
                        &self.cell,
                        &atoms,
                        &self.config.symmetry_full_mesh().map_err(|e| error(1, e))?,
                        self.config.detection_tolerance(),
                        self.config.basis.time_reversal(),
                        self.config.basis.fft_grid,
                    )
                    .map_err(|e| error(4, e))?,
                );
            }
            request.derivatives.forces = true;
            request.derivatives.deferred = true;
            let mut delivery = crate::progress::solve(
                dft_engine::refstart::StartInput {
                    cell: &self.cell,
                    atoms: &atoms,
                    potentials: &refs,
                    options: &self.opts,
                },
                request,
            )
            .map_err(|e| error(4, e))?;
            let mut run = delivery.run.into_scalar().map_err(|e| error(4, e))?;
            delivery
                .derivatives
                .prepare(
                    dft_engine::refstart::StartInput {
                        cell: &self.cell,
                        atoms: &atoms,
                        potentials: &refs,
                        options: &self.opts,
                    },
                    &mut run.outcome,
                )
                .map_err(|e| error(4, e))?;
            reference_derivatives = Some(delivery.derivatives);
            (run.outcome.result, n_at, provenance, fallback, clipped)
        } else {
            let atomic_here = atomic && warm.is_some();
            let problem = self.problem_with(atomic_here)?;
            let n_at: Option<Vec<f64>> = if !atomic {
                None
            } else if atomic_here {
                problem
                    .atomic_superposition_density()
                    .map(|a| a.iter().copied().collect())
            } else {
                // A cold evaluation keeps its own start; n_at is built aside.
                self.problem_with(true)?
                    .atomic_superposition_density()
                    .map(|a| a.iter().copied().collect())
            };
            let mut density = problem.uniform_density();
            let mut session =
                ScfSession::new(problem, CompletionMode::Strict, ExecutionConfig::default())
                    .map_err(|e| error(4, e))?;
            let mut provenance = if warm.is_some() {
                "warm_transfer"
            } else {
                "fresh"
            };
            let mut start_fallback = None;
            let mut start_clipped_fraction = None;
            if let Some((source, values)) = warm {
                let source = source
                    .ok_or_else(|| error(4, "legacy warm transfer requires its source model"))?;
                if values.len() != density.len() {
                    return Err(error(4, "warm_transfer: density size mismatch"));
                }
                for (d, v) in density.iter_mut().zip(values) {
                    *d = *v;
                }
                if let Some((bonding, kind)) = predicted {
                    match &n_at {
                        Some(n_at)
                            if n_at.len() == bonding.len() && n_at.len() == density.len() =>
                        {
                            // Project onto the feasible set the same way the atomic
                            // start does: clip negative samples, restore the charge.
                            let mut start: Vec<f64> =
                                n_at.iter().zip(bonding).map(|(a, d)| a + d).collect();
                            let total: f64 = start.iter().sum();
                            let negative: f64 = start.iter().filter(|v| **v < 0.).sum();
                            start.iter_mut().for_each(|v| *v = v.max(0.));
                            let kept: f64 = start.iter().sum();
                            if total.is_finite() && total > 0. && kept > 0. {
                                for (d, v) in density.iter_mut().zip(&start) {
                                    *d = v * (total / kept);
                                }
                                provenance = kind;
                                if negative < 0. {
                                    start_clipped_fraction = Some(-negative / total);
                                }
                            } else {
                                start_fallback =
                                    Some("extrapolated density not finite".to_string());
                            }
                        }
                        Some(_) => start_fallback = Some("predictor size mismatch".into()),
                        None => {
                            start_fallback = Some("no atomic charge density (PP_RHOATOM)".into());
                        }
                    }
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
            (
                result.result,
                n_at,
                provenance,
                start_fallback,
                start_clipped_fraction,
            )
        };
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
        let forces = if let Some(derivatives) = &mut reference_derivatives {
            derivatives
                .evaluate_forces_with_operations(
                    dft_engine::refstart::StartInput {
                        cell: &self.cell,
                        atoms: &self.atoms(),
                        potentials: &self.refs(),
                        options: &self.opts,
                    },
                    &state,
                    &ops,
                )
                .map_err(|e| error(4, e))?;
            derivatives
                .forces
                .take()
                .expect("requested reference forces")
        } else {
            dft_engine::forces(&dft_engine::ForceState {
                cell: &self.cell,
                atoms: &self.atoms(),
                potentials: &self.refs(),
                options: &self.opts,
                electronic: dft_engine::ElectronicState::Fixed(&state),
                density_floor: None,
                symmetry: &ops,
            })
            .map_err(|e| error(4, e))?
        };
        // D3-DISPERSION (ADR 0048): the D3 term of this geometry joins the
        // physical energy, forces and stress before any projection.
        let dispersion =
            crate::run::dispersion_term(&self.config, &self.cell, &self.atoms(), &self.refs())?;
        let mut forces = forces;
        if let Some(d) = &dispersion {
            crate::run::add_dispersion_forces(&mut forces, d);
        }
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
                density_floor: None,
                symmetry: &ops,
            };
            // Nominal-cutoff frames (final rebuild, basis rebuilds) carry the
            // opt-in Pulay assessment (PULAY-OPTIN, default off).
            let mode = self.config.pulay_assessment();
            let stress = if self.opts.strain_reference.is_none() {
                let prepared = if mode == dft_engine::PulayAssessmentMode::CutoffScan
                    && reference
                    && self.config.symmetry.mode == crate::config::SymmetryMode::Auto
                {
                    Some(
                        dft_engine::refscf::symmetry::Prepared::new_for_grid(
                            &self.cell,
                            &atoms,
                            &self.config.symmetry_full_mesh().map_err(|e| error(1, e))?,
                            self.config.detection_tolerance(),
                            self.config.basis.time_reversal(),
                            self.config.basis.fft_grid,
                        )
                        .map_err(|e| error(4, e))?,
                    )
                } else {
                    None
                };
                pulay_assessment = (mode == dft_engine::PulayAssessmentMode::Off).then(
                    || json!({"status":"not_evaluated","reason":crate::run::PULAY_NOT_EVALUATED}),
                );
                dft_engine::stress::stress_with(&view, prepared.as_ref(), mode)
            } else {
                dft_engine::stress::frozen_stress(&view)
            }
            .map_err(|e| error(4, e))?;
            let mut stress = stress;
            if let Some(d) = &dispersion {
                crate::run::add_dispersion_stress(&mut stress, d);
            }
            if let Some(p) = stress.pulay {
                pulay_assessment = Some(json!({
                    "method":p.method, "correction_ha_bohr3":p.correction,
                    "uncertainty_ha_bohr3":p.uncertainty, "cutoff_steps_ha":p.cutoff_steps_ha,
                    "shells_resolved":p.shells_resolved, "warning":p.warning,
                }));
            }
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
            energy_ha: dispersion.as_ref().map_or(state.total_hartree(), |d| {
                state.total_hartree() + d.energy_ha
            }),
            forces_ha_bohr: forces.values,
            max_force_ha_bohr: projected.iter().map(|v| v.abs()).fold(0., f64::max),
            max_force_norm_ha_bohr: max_atom_norm(&projected),
            step_length_bohr: length,
            scf_iterations: state.history.len(),
            provenance: provenance.into(),
            scf_h_vectors: Some(
                state
                    .history
                    .iter()
                    .map(|h| {
                        h.operator_applications
                            + h.fixed_point_certificate
                                .as_ref()
                                .map_or(0, |c| c.operator_applications)
                    })
                    .sum(),
            ),
            start_fallback,
            start_clipped_fraction,
            verification: verification_summary(&state, &self.opts),
        };
        if reference && settings.warm_start {
            *self.reference_seed.borrow_mut() = Some((
                self.cell.clone(),
                reference_derivatives
                    .as_ref()
                    .expect("reference delivery")
                    .ionic_seed(
                        dft_engine::refstart::StartInput {
                            cell: &self.cell,
                            atoms: &self.atoms(),
                            potentials: &self.refs(),
                            options: &self.opts,
                        },
                        &state,
                    )
                    .map_err(|e| error(4, e))?,
            ));
        }
        Ok((frame, state.density.iter().copied().collect(), n_at))
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
    let config = Config::parse_at(std::str::from_utf8(&bytes).map_err(|e| error(1, e))?, input)
        .map_err(|e| error(1, e))?;
    if config.compat.reference() {
        return Err(error(1, crate::config::REFERENCE_SCOPE));
    }
    let mut settings = config.relax.clone().unwrap_or_default();
    // SI1: per-atom fixed flags of [structure] (or POSCAR selective
    // dynamics); parse rejects them together with relax.fixed.
    if settings.fixed.is_empty() {
        settings.fixed = config.structure_fixed.clone();
    }
    settings.start_extrapolation = settings.start_extrapolation.effective()?;
    settings.validate(config.system.atoms.len())?;
    // ONEPATH-1/ONEPATH-2 (ADR 0043 item 5): every relaxation mode runs each
    // ionic or cell step on the standard path with the fixed-point
    // verification; forces and stress come from the verified step.
    let scf_path = config.scf_path().map_err(|e| error(1, e))?;
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
    let upfs: Vec<_> = potentials.iter().map(|(_, p)| p).collect();
    let selection = crate::driver_select::select(
        &config,
        crate::driver_select::Mode::Classic,
        crate::driver_select::Command::Relax,
        &upfs,
    );
    selection.announce();
    let opts = crate::run::options_for(&config, selection)?;
    if opts.energy_tolerance > 1e-8 || opts.eigen_tolerance > 1e-8 || opts.density_tolerance > 1e-7
    {
        return Err(error(
            1,
            "relax requires SCF tolerances no looser than the force preconditions (energy 1e-8 Ha, eigen 1e-8 Ha, density 1e-7)",
        ));
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
        scf_path,
        selection,
        reference_seed: std::cell::RefCell::new(None),
        provenance: json!({"input_sha256":sha256(&bytes),"executable_sha256":build_sha256,"pseudopotentials_sha256":pp_hashes,"platform":std::env::consts::OS,"architecture":std::env::consts::ARCH}),
    };
    if settings.mode != cell::Mode::Ions {
        let result = cell::execute(ctx, settings, identity, out, resume);
        let path = out.join("results.json");
        if path.is_file() {
            let mut report: Value =
                serde_json::from_slice(&fs::read(&path)?).map_err(|e| error(5, e))?;
            selection.record(&mut report);
            atomic(
                &path,
                &serde_json::to_vec_pretty(&report).map_err(|e| error(5, e))?,
            )?;
        }
        return result;
    }
    let n = ctx.config.system.atoms.len() * 3;
    let extrapolate = settings.warm_start && settings.start_extrapolation != Extrapolation::None;
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
            .evaluate_with(positions.clone(), 0., None, None, extrapolate, &settings)
            .or_else(|e| {
                failures.push(
                    json!({"step":0,"attempt":1,"positions_bohr":positions,"reason":e.to_string()}),
                );
                ctx.evaluate_with(positions.clone(), 0., None, None, extrapolate, &settings)
            });
        let (frame, density, n_at) = match initial {
            Ok(v) => v,
            Err(e) => {
                failures.push(
                    json!({"step":0,"attempt":2,"positions_bohr":positions,"reason":e.to_string()}),
                );
                let mut report = json!({"schema_version":1,"calculation_type":"relaxation","provenance":ctx.provenance,"relaxation":{"status":"failed","steps":0,"last_good_structure":null,"final_max_force_ha_bohr":null,"constraints":settings.fixed,"failures":failures}});
                selection.record(&mut report);
                atomic(
                    &out.join("results.json"),
                    &serde_json::to_vec_pretty(&report).map_err(|e| error(5, e))?,
                )?;
                return Err(e);
            }
        };
        let mut extrapolation = Vec::new();
        push_history(&mut extrapolation, &frame.positions, &density, n_at);
        // ADR 0056: without the energy requirement an initial structure that
        // already meets the force test is converged after zero ionic steps.
        let (consecutive, converged) = settings.ionic_stop(&frame, None, 0);
        Restart {
            version: 1,
            identity,
            optimizer: settings.optimizer(n),
            frames: vec![frame],
            density,
            consecutive,
            status: if converged { "converged" } else { "running" }.into(),
            failures,
            extrapolation,
        }
    };
    ctx.position(&state.frames.last().unwrap().fractional);
    // Publish accepted state before derived outputs; resume regenerates outputs.
    save(&restart_path, &state)?;
    outputs(&ctx, &settings, &state, out)?;
    let rise_tol = settings.energy_rise_tol_ha(&ctx.opts);
    while state.frames.len() - 1 < settings.max_steps && state.status == "running" {
        // Wall of this ionic step (proposal, every attempt, outputs); printed on
        // the step line only: frames and relax.log stay deterministic for resume.
        let started = std::time::Instant::now();
        let last = state.frames.last().unwrap().clone();
        ctx.position(&last.fractional);
        let source = if ctx.selection.reference() {
            None
        } else {
            Some(ctx.problem()?)
        };
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
                .then_some((source.as_ref(), state.density.as_slice()));
            let predicted = if extrapolate && warm.is_some() {
                predict(
                    &state.extrapolation,
                    &positions,
                    settings.start_extrapolation,
                )
            } else {
                None
            };
            match ctx.evaluate_with(
                positions.clone(),
                length,
                warm,
                predicted.as_ref().map(|(d, kind)| (d.as_slice(), *kind)),
                extrapolate,
                &settings,
            ) {
                Ok((frame, density, n_at)) => {
                    accepted = Some((frame, density, n_at, attempt));
                    break;
                }
                Err(e) => state.failures.push(json!({"step":state.frames.len(),"attempt":attempt+1,"positions_bohr":positions,"fractional":fractional,"reason":e.to_string()})),
            }
        }
        if let Some((frame, density, n_at, attempt)) = accepted {
            if extrapolate {
                push_history(&mut state.extrapolation, &frame.positions, &density, n_at);
            }
            let delta = frame.energy_ha - last.energy_ha;
            let (consecutive, converged) =
                settings.ionic_stop(&frame, Some(delta), state.consecutive);
            state.consecutive = consecutive;
            if attempt > 0 || delta > rise_tol {
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
            if converged {
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
        outputs_with_wall(&ctx, &settings, &state, out, Some(started.elapsed()))?;
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

fn outputs(
    ctx: &Context,
    settings: &Settings,
    state: &Restart,
    out: &Path,
) -> Result<(), RunError> {
    outputs_with_wall(ctx, settings, state, out, None)
}
#[allow(clippy::needless_range_loop)] // Lattice vectors are columns, file vectors are rows.
fn outputs_with_wall(
    ctx: &Context,
    settings: &Settings,
    state: &Restart,
    out: &Path,
    wall: Option<std::time::Duration>,
) -> Result<(), RunError> {
    // UX1 progress stream (ADR 0053): one record per new frame or failure.
    if let Some(f) = state.frames.last() {
        crate::progress::relax_frame(
            state.frames.len(),
            f.energy_ha,
            f.max_force_ha_bohr,
            f.step_length_bohr,
            f.scf_iterations,
            &f.provenance,
            f.stress_ha_bohr3,
            &state.status,
            &state.failures,
        );
    }
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
    // Columns 1-6 unchanged since P21A; max_force_Ha_Bohr is the largest free
    // Cartesian component, max_force_norm_Ha_Bohr the largest atomic norm
    // (ADR 0056; the stop rule uses relaxation.stop_rule.force_criterion).
    let mut log = String::from(
        "step energy_Ha max_force_Ha_Bohr step_length_Bohr scf_iterations provenance max_force_norm_Ha_Bohr\n",
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
            "{i} {:.17e} {:.17e} {:.17e} {} {} {:.17e}",
            frame.energy_ha,
            frame.max_force_ha_bohr,
            frame.step_length_bohr,
            frame.scf_iterations,
            frame.provenance,
            frame.max_force_norm_ha_bohr
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
    let mut result = json!({"schema_version":1,"calculation_type":"relaxation","provenance":ctx.provenance,
        "single_path":crate::run::single_path_json(&ctx.scf_path, &ctx.opts),
        "scf":{"effective_profile":if ctx.opts.standard.is_some() { "standard" } else { "strict" },
            "standard_threshold_ha":ctx.opts.standard.map(|p| p.threshold_ha),
            "threshold_source":ctx.scf_path.threshold_source,
            "verification":if ctx.opts.standard.is_some() { ctx.opts.verification.name() } else { "strict_certificate" },
            "note":"every ionic or cell step runs this SCF; forces and stress of a standard step rest on its passed fixed-point verification (frames[].verification); tighten scf.standard_threshold_ha for force- or stress-sensitive runs"},
        "relaxation":{"status":state.status,"steps":state.frames.len()-1,"final_max_force_ha_bohr":last.max_force_ha_bohr,
        "final_max_force_norm_ha_bohr":last.max_force_norm_ha_bohr,"stop_rule":stop_rule_json(settings, false),
        "energy_history_ha":state.frames.iter().map(|f|f.energy_ha).collect::<Vec<_>>(),"constraints":settings.fixed,"frames":state.frames,"failures":state.failures,
        "restart_semantics":"same_build_ionic_continuation; each transferred SCF is warm_transfer, not exact SCF replay","settings":settings}});
    ctx.selection.record(&mut result);
    if let Some(d) = &ctx.config.dispersion {
        result["dispersion"] = json!({"method":d.method,
            "included_in":"every frame's energy_ha, forces_ha_bohr and stress_ha_bohr3 (ADR 0048)"});
    }
    if let Some(provenance) = ctx.config.system.provenance_json() {
        result["structure_provenance"] = provenance;
    }
    ctx.config.annotate_structure(&mut result);
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
        "step {} E={:.12e} Ha maxF={:.6e} Ha/Bohr max|F|={:.6e} Ha/Bohr displacement={:.6e} Bohr SCF={} {}{}",
        state.frames.len() - 1,
        last.energy_ha,
        last.max_force_ha_bohr,
        last.max_force_norm_ha_bohr,
        last.step_length_bohr,
        last.scf_iterations,
        state.status,
        wall.map_or(String::new(), |w| format!(" wall={:.3}s", w.as_secs_f64()))
    );
    Ok(())
}
/// The stop rule that applied, for results.json (ADR 0056).
fn stop_rule_json(settings: &Settings, cell: bool) -> Value {
    let energy = if cell {
        Some(settings.cell_energy_tol_ha())
    } else {
        settings.energy_tol_ha
    };
    json!({
        "force_criterion": settings.force_criterion,
        "force_measure": match settings.force_criterion {
            ForceCriterion::Norm => "max over atoms of |F_i| (projected), Ha/Bohr",
            ForceCriterion::Component => "max free Cartesian |F_ia| (projected), Ha/Bohr",
        },
        "force_tol_ha_bohr": settings.force_tol_ha_bohr,
        "energy_tol_ha": energy,
        "consecutive_steps_required": if energy.is_some() { 2 } else { 1 },
        "inequality": "strict",
        "adr": "0056",
    })
}

#[cfg(test)]
mod extrapolation_tests {
    use super::*;
    /// A bonding density linear in the positions, d = d0 + M R (four grid samples).
    fn bonding(r: &[f64]) -> Vec<f64> {
        let m = [[0.3, -0.2], [1.1, 0.4], [-0.7, 0.9], [0.05, -0.6]];
        let d0 = [0.01, -0.02, 0.03, 0.0];
        (0..4)
            .map(|g| d0[g] + m[g][0] * r[0] + m[g][1] * r[1])
            .collect()
    }
    fn entry(r: [f64; 2]) -> HistoryEntry {
        HistoryEntry {
            positions: r.to_vec(),
            bonding: bonding(&r),
        }
    }
    #[test]
    fn second_order_is_exact_when_the_step_lies_in_the_span_of_the_last_two() {
        // A curved ionic path R(t) = (t, t^2/2): s = R(3) - R(2) lies in the span
        // of u = R(2) - R(1) and v = R(1) - R(0), so the fit is exact and the
        // prediction equals d(R(3)) up to rounding.
        let path = |t: f64| [t, 0.5 * t * t];
        let history: Vec<_> = (0..3).map(|t| entry(path(t as f64))).collect();
        let (d, kind) = predict(&history, &path(3.), Extrapolation::SecondOrder).unwrap();
        assert_eq!(kind, "extrapolated_second_order");
        for (p, e) in d.iter().zip(bonding(&path(3.))) {
            assert!((p - e).abs() < 1e-12, "{p} vs {e}");
        }
        // First order carries d_k unchanged.
        let (d1, kind1) = predict(&history, &path(3.), Extrapolation::FirstOrder).unwrap();
        assert_eq!(kind1, "extrapolated_first_order");
        assert_eq!(d1, history[2].bonding);
    }
    #[test]
    fn short_or_degenerate_histories_reduce_the_order() {
        assert!(predict(&[], &[0., 0.], Extrapolation::SecondOrder).is_none());
        let one = [entry([0., 0.])];
        let (d, kind) = predict(&one, &[1., 0.], Extrapolation::SecondOrder).unwrap();
        assert_eq!(
            (d, kind),
            (one[0].bonding.clone(), "extrapolated_first_order")
        );
        // Two steps on a straight line: one-term fit, exact for linear d.
        let two = [entry([0., 0.]), entry([1., 1.])];
        let (d, _) = predict(&two, &[2.5, 2.5], Extrapolation::SecondOrder).unwrap();
        for (p, e) in d.iter().zip(bonding(&[2.5, 2.5])) {
            assert!((p - e).abs() < 1e-12);
        }
        // Three collinear steps: the Gram determinant vanishes, the v term is
        // dropped and the one-term fit is still exact.
        let three = [entry([0., 0.]), entry([1., 0.]), entry([3., 0.])];
        let (d, _) = predict(&three, &[4., 0.], Extrapolation::SecondOrder).unwrap();
        for (p, e) in d.iter().zip(bonding(&[4., 0.])) {
            assert!((p - e).abs() < 1e-12);
        }
        // A repeated position falls back to first order.
        let same = [entry([1., 0.]), entry([1., 0.])];
        let (_, kind) = predict(&same, &[2., 0.], Extrapolation::SecondOrder).unwrap();
        assert_eq!(kind, "extrapolated_first_order");
    }
    #[test]
    fn history_keeps_three_and_clears_without_atomic_density() {
        let mut h = Vec::new();
        for i in 0..5 {
            push_history(&mut h, &[i as f64], &[1., 2.], Some(vec![0.5, 0.5]));
        }
        assert_eq!(h.len(), 3);
        assert_eq!(h[0].positions, vec![2.]);
        assert_eq!(h[2].bonding, vec![0.5, 1.5]);
        push_history(&mut h, &[9.], &[1., 2.], None);
        assert!(h.is_empty());
        let settings: Settings = toml::from_str("start_extrapolation = \"second_order\"").unwrap();
        assert_eq!(settings.start_extrapolation, Extrapolation::SecondOrder);
        assert_eq!(Settings::default().start_extrapolation, Extrapolation::None);
        assert!(toml::from_str::<Settings>("start_extrapolation = \"cubic\"").is_err());
    }
}

#[cfg(test)]
mod stop_rule_tests {
    //! ADR 0056 stop rule: units Ha/Bohr for forces, Ha for energies; all
    //! comparisons below are exact (no tolerance) because the measures are
    //! computed from the given numbers only.
    use super::*;
    /// A frame whose projected force is `forces` under `settings` constraints,
    /// built the same way `evaluate_with` builds the recorded maxima.
    fn frame(settings: &Settings, forces: &[[f64; 3]], energy: f64) -> Frame {
        let mut projected: Vec<f64> = forces.iter().flatten().copied().collect();
        settings.project(&mut projected);
        Frame {
            positions: vec![0.; 3 * forces.len()],
            cell_columns_bohr: None,
            stress_ha_bohr3: None,
            pulay_assessment: None,
            fractional: vec![[0.; 3]; forces.len()],
            energy_ha: energy,
            forces_ha_bohr: forces.to_vec(),
            max_force_ha_bohr: projected.iter().map(|v| v.abs()).fold(0., f64::max),
            max_force_norm_ha_bohr: max_atom_norm(&projected),
            step_length_bohr: 0.,
            scf_iterations: 0,
            provenance: "test".into(),
            scf_h_vectors: None,
            start_fallback: None,
            start_clipped_fraction: None,
            verification: None,
        }
    }
    #[test]
    fn defaults_parse_and_validate() {
        let d = Settings::default();
        assert_eq!(d.force_criterion, ForceCriterion::Norm);
        assert_eq!(d.energy_tol_ha, None);
        assert_eq!((d.fire_dt, d.fire_dt_max), (0.1, 1.0));
        // 0.02 eV/A = 0.02 / (27.211386 eV / 0.52917721 A) Ha/Bohr.
        let expected = 0.02 * 0.529_177_210_9 / 27.211_386_245_988;
        assert!(
            (d.force_tol_ha_bohr - expected).abs() < 1e-8 * expected,
            "{}",
            d.force_tol_ha_bohr
        );
        assert!((d.force_tol_ha_bohr - 3.889e-4).abs() < 1e-7);
        // Old inputs keep parsing: energy_tol_ha now opts into the energy test.
        let old: Settings =
            toml::from_str("force_tol_ha_bohr = 1e-3\nenergy_tol_ha = 1e-8").unwrap();
        assert_eq!(old.energy_tol_ha, Some(1e-8));
        let s: Settings =
            toml::from_str("force_criterion = \"component\"\nfire_dt = 0.2\nfire_dt_max = 0.5")
                .unwrap();
        assert_eq!(s.force_criterion, ForceCriterion::Component);
        assert!(s.validate(1).is_ok());
        assert!(toml::from_str::<Settings>("force_criterion = \"rms\"").is_err());
        for bad in [
            "energy_tol_ha = 0.0",
            "energy_tol_ha = -1e-8",
            "fire_dt = 0.0",
            "fire_dt = 2.0",
            "fire_dt_max = -1.0",
        ] {
            let s: Settings = toml::from_str(bad).unwrap();
            assert!(s.validate(1).is_err(), "{bad}");
        }
    }
    #[test]
    fn norm_is_stricter_than_component_by_up_to_sqrt3() {
        // One atom with equal components c: component = c, norm = sqrt(3) c.
        let c = 1e-3 / 1.5;
        let norm = Settings {
            force_tol_ha_bohr: 1e-3,
            ..Settings::default()
        };
        let component = Settings {
            force_criterion: ForceCriterion::Component,
            ..norm.clone()
        };
        let f = frame(&norm, &[[c, -c, c], [1e-5, 0., 0.]], 0.);
        assert_eq!(f.max_force_ha_bohr, c);
        assert_eq!(f.max_force_norm_ha_bohr, (3. * c * c).sqrt());
        assert!(component.force_converged(&f));
        assert!(!norm.force_converged(&f));
        // Along one axis the two measures coincide.
        let axis = frame(&norm, &[[0., 9e-4, 0.]], 0.);
        assert!(norm.force_converged(&axis) && component.force_converged(&axis));
        // Strict inequality at the tolerance.
        let edge = frame(&norm, &[[1e-3, 0., 0.]], 0.);
        assert!(!norm.force_converged(&edge));
    }
    #[test]
    fn constraints_are_projected_before_the_norm() {
        let settings = Settings {
            fixed: vec![
                Fixed {
                    atom: 0,
                    components: [true; 3],
                },
                Fixed {
                    atom: 1,
                    components: [false, true, true],
                },
            ],
            force_tol_ha_bohr: 1e-3,
            ..Settings::default()
        };
        // Large forces on frozen components do not count; atom 1 keeps only x.
        let f = frame(
            &settings,
            &[[1., 1., 1.], [5e-4, 1., 1.], [3e-4, 4e-4, 0.]],
            0.,
        );
        assert_eq!(f.max_force_norm_ha_bohr, 5e-4);
        assert!(settings.force_converged(&f));
        // Raw physical forces stay in the frame.
        assert_eq!(f.forces_ha_bohr[0], [1., 1., 1.]);
    }
    #[test]
    fn energy_requirement_is_optional() {
        let off = Settings::default();
        let small = frame(&off, &[[1e-4, 0., 0.]], -1.);
        let large = frame(&off, &[[1e-2, 0., 0.]], -1.);
        // Off (default): the first frame that passes the force test stops,
        // the initial frame (no energy change) included.
        assert_eq!(off.ionic_stop(&small, None, 0), (1, true));
        assert_eq!(off.ionic_stop(&small, Some(-0.5), 0), (1, true));
        assert_eq!(off.ionic_stop(&large, Some(0.), 1), (0, false));
        // On: force and |dE| < tol on two consecutive accepted steps; the
        // initial frame never counts (the rule before ADR 0056).
        let on = Settings {
            energy_tol_ha: Some(1e-8),
            ..Settings::default()
        };
        assert_eq!(on.ionic_stop(&small, None, 0), (0, false));
        assert_eq!(on.ionic_stop(&small, Some(5e-9), 0), (1, false));
        assert_eq!(on.ionic_stop(&small, Some(-5e-9), 1), (2, true));
        assert_eq!(on.ionic_stop(&small, Some(2e-8), 1), (0, false));
        assert_eq!(on.ionic_stop(&large, Some(0.), 1), (0, false));
        // Trust-radius energy-rise tolerance: explicit value, else 2x the SCF
        // energy accuracy.
        let mut strict = ScfOptions::strict_reference();
        strict.energy_tolerance = 1e-9;
        assert_eq!(on.energy_rise_tol_ha(&strict), 1e-8);
        assert_eq!(off.energy_rise_tol_ha(&strict), 2e-9);
        let standard = ScfOptions::default();
        let threshold = standard.standard.unwrap().threshold_ha;
        assert_eq!(off.energy_rise_tol_ha(&standard), 2. * threshold);
        // Cell modes keep their energy rule (default 1e-8 Ha).
        assert_eq!(off.cell_energy_tol_ha(), 1e-8);
        assert_eq!(on.cell_energy_tol_ha(), 1e-8);
        let rule = stop_rule_json(&off, false);
        assert_eq!(rule["force_criterion"], "norm");
        assert!(rule["energy_tol_ha"].is_null());
        assert_eq!(rule["consecutive_steps_required"], 1);
        assert_eq!(stop_rule_json(&off, true)["consecutive_steps_required"], 2);
    }
}

#[cfg(test)]
mod reference_relax_tests {
    use super::*;

    #[test]
    #[ignore = "external Si2, 3 ionic steps; DFT_SI_UPF=/path/Si.pz-vbc.UPF cargo test --release -p dft-cli reference_relax_tests -- --ignored --nocapture"]
    fn reference_three_step_si2_matches_strict_oracle() {
        let pp = std::env::var("DFT_SI_UPF").expect("external Si UPF required");
        let root = std::env::temp_dir().join(format!("g5-si2-relax-{}", std::process::id()));
        fs::create_dir_all(&root).unwrap();
        let mut trajectories = Vec::new();
        for driver in ["strict_oracle", "reference"] {
            let switches = dft_numerics::switches::RunConfig::from_pairs(&[(
                "reference_strict",
                if driver == "strict_oracle" {
                    "on"
                } else {
                    "off"
                },
            )])
            .unwrap();
            let _scope = switches.enter();
            let input = format!(
                r#"
[system]
name='G5 Si2'
lattice_unit='bohr'
lattice_vectors=[[0.0,5.13,5.13],[5.13,0.0,5.13],[5.13,5.13,0.0]]
[[system.atoms]]
species='Si'
fractional=[0.0,0.0,0.0]
[[system.atoms]]
species='Si'
fractional=[0.253,0.25,0.25]
[pseudopotentials]
Si={pp:?}
[basis]
ecutwfc_ry=8.0
kpoints=[[0.0,0.0,0.0]]
[xc]
functional='lda-pz'
[scf]
driver='reference'
nbands=5
max_iterations=200
standard_threshold_ha=5e-13
start_density='atomic'
start_wavefunctions='atomic'
[relax]
max_steps=3
force_tol_ha_bohr=1e-12
energy_tol_ha=1e-12
"#
            );
            let path = root.join(format!("{driver}.toml"));
            fs::write(&path, input).unwrap();
            let out = root.join(driver);
            let result = execute(&path, &out, false);
            assert!(
                result.is_ok() || result.as_ref().is_err_and(|e| e.code == 3),
                "{driver}: {result:?}"
            );
            let restart: Value =
                serde_json::from_slice(&fs::read(out.join("restart.json")).unwrap()).unwrap();
            let mut payload: Value =
                serde_json::from_str(restart["payload"].as_str().unwrap()).unwrap();
            bits(&mut payload, false).unwrap();
            let state: Restart = serde_json::from_value(payload).unwrap();
            assert_eq!(state.frames.len(), 4);
            if driver == "reference" {
                for frame in &state.frames[1..] {
                    assert_eq!(frame.provenance, "warm_transfer");
                }
            }
            trajectories.push(state.frames);
        }
        for (step, (old, new)) in trajectories[0].iter().zip(&trajectories[1]).enumerate() {
            let de = (old.energy_ha - new.energy_ha).abs();
            let df = old
                .forces_ha_bohr
                .iter()
                .flatten()
                .zip(new.forces_ha_bohr.iter().flatten())
                .map(|(a, b)| (a - b).abs())
                .fold(0., f64::max);
            let dx = old
                .positions
                .iter()
                .zip(&new.positions)
                .map(|(a, b)| (a - b).abs())
                .fold(0., f64::max);
            println!("G5 step={step} dE={de:e} Ha dF={df:e} Ha/Bohr dR={dx:e} Bohr");
            assert!(de < 5e-9 && df < 1e-6 && dx < 1e-6);
        }
    }
}
