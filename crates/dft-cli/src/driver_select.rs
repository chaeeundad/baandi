//! Cold-path routing policy. Resolve again with parsed UPFs before constructing
//! an engine. An engine error is never a reason to retry with another driver.
use crate::config::Config;
use dft_numerics::switches::{self, Id, RunConfig};
use dft_upf::Pseudopotential;
use serde::Serialize;
use serde_json::{Value, json};

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
pub enum Driver {
    #[serde(rename = "reference")]
    Reference,
    #[serde(rename = "legacy")]
    LegacyStandard,
    #[serde(rename = "strict_oracle")]
    StrictOracle,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Source {
    Default,
    Env,
    Input,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Mode {
    Classic,
    Observe,
    Agent,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Command {
    Run,
    Relax,
    Bands,
    Dos,
    Checkpoint,
    Restart,
    Explain,
    RefscfMap { reference: Option<bool> },
    Session,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Selection {
    pub driver: Driver,
    pub reason: Option<&'static str>,
    pub source: Source,
}
impl Selection {
    pub fn reference(self) -> bool {
        self.driver == Driver::Reference
    }
    pub fn strict(self) -> bool {
        self.driver == Driver::StrictOracle
    }
    /// The G4B/G4D routes predate R1 and already accepted reference mixer
    /// aliases/defaults before falling back. Keep that input contract intact.
    pub(crate) fn reference_mixing(self) -> bool {
        self.reference() || matches!(self.reason, Some(G4B_FALLBACK | G4D_FALLBACK))
    }
    pub fn json(self) -> Value {
        json!({"driver":self.driver,"driver_reason":self.reason,"driver_source":self.source})
    }
    pub fn record(self, result: &mut Value) {
        if !result["scf"].is_object() {
            result["scf"] = json!({});
        }
        result["scf"]
            .as_object_mut()
            .unwrap()
            .extend(self.json().as_object().unwrap().clone());
        if self.driver == Driver::LegacyStandard && self.reason.is_some() {
            result["reference_driver_fallback_reason"] = json!(self.reason);
        }
    }
    pub fn announce(self) {
        eprintln!("SCF driver: {}", self.json());
    }
}

pub const STRICT: &str = "BAANDI_REFERENCE_STRICT=on (strict reference path)";
pub const G4D_FALLBACK: &str =
    "reference driver: free-moment collinear FD with PBE not yet qualified, G4D";
pub const G4B_FALLBACK: &str = "reference driver: collinear Fermi-Dirac not yet qualified, G4B";
pub const VARIABLE_CELL: &str =
    "L4: variable-cell relaxation requires legacy frozen-strain support";
pub const EXTRAPOLATION: &str = "L5: relaxation density extrapolation requires the legacy driver";
pub const CHECKPOINT: &str = "L6: checkpoint/restart v1 requires the legacy session";
pub const RUNTIME: &str = "L7: observe/agent runtime requires the legacy session";
pub const SESSION: &str = "L8: step/status/iteration APIs require the legacy session";
pub const LEGACY_INPUT: &str = "L9: explicit legacy-only SCF controls";
pub const COMPAT: &str = "L10: original PW92 prefactor requires the legacy driver";
pub const ATOMIC_DATA: &str =
    "L11: reference start requires nonzero PP_RHOATOM and PP_PSWFC with l<=3";

fn source(id: Id) -> Source {
    match RunConfig::current().source(id) {
        switches::Source::Input => Source::Input,
        switches::Source::Default => Source::Default,
        _ => Source::Env,
    }
}
/// Highest-priority request, without feature eligibility (also used by validation).
pub fn requested(config: &Config) -> Selection {
    if switches::on(Id::REFERENCE_STRICT) {
        return Selection {
            driver: Driver::StrictOracle,
            reason: Some(STRICT),
            source: source(Id::REFERENCE_STRICT),
        };
    }
    if let Some(driver) = config.scf.as_ref().and_then(|s| s.driver.as_deref()) {
        return Selection {
            driver: if driver == "reference" {
                Driver::Reference
            } else {
                Driver::LegacyStandard
            },
            reason: None,
            source: Source::Input,
        };
    }
    Selection {
        driver: if switches::on(Id::REFERENCE_SCF) {
            Driver::Reference
        } else {
            Driver::LegacyStandard
        },
        reason: None,
        source: source(Id::REFERENCE_SCF),
    }
}
fn spin_fd(c: &Config) -> bool {
    c.collinear() && matches!(c.finite_temperature(), Ok(Some(_)))
}
fn l1(c: &Config, _: Mode, _: Command, _: &[&Pseudopotential]) -> bool {
    spin_fd(c)
        && c.xc.functional.engine() == dft_engine::Functional::Pbe
        && matches!(
            c.finite_constraint(),
            dft_physics::occupations::Constraint::Total
        )
}
fn l2(c: &Config, _: Mode, _: Command, _: &[&Pseudopotential]) -> bool {
    spin_fd(c)
        && !(switches::on(Id::REFSCF_SPIN_FD_MOMENT)
            && switches::on(Id::REFSCF_SPIN_STIFFNESS)
            && switches::on(Id::REFSCF_HISTORY_CONDITIONING))
}
fn l4(c: &Config, _: Mode, _: Command, _: &[&Pseudopotential]) -> bool {
    c.relax
        .as_ref()
        .is_some_and(|r| r.mode != crate::relax::cell::Mode::Ions)
}
fn l5(c: &Config, _: Mode, command: Command, _: &[&Pseudopotential]) -> bool {
    (command == Command::Relax || c.relax.is_some())
        && match switches::text(Id::RELAX_EXTRAPOLATION) {
            Some(v) => v != "none",
            None => c
                .relax
                .as_ref()
                .is_some_and(|r| r.start_extrapolation != crate::relax::Extrapolation::None),
        }
}
fn l6(_: &Config, _: Mode, command: Command, _: &[&Pseudopotential]) -> bool {
    matches!(command, Command::Checkpoint | Command::Restart)
}
fn l7(_: &Config, mode: Mode, _: Command, _: &[&Pseudopotential]) -> bool {
    mode != Mode::Classic
}
fn l8(_: &Config, _: Mode, command: Command, _: &[&Pseudopotential]) -> bool {
    command == Command::Session
}
fn l9(c: &Config, _: Mode, _: Command, _: &[&Pseudopotential]) -> bool {
    c.xc.density_floor.is_some()
        || c.scf.as_ref().is_some_and(|s| {
            matches!(s.mixing.as_deref(), Some("linear" | "kerker-pulay"))
                || s.mixing_q0_inverse_bohr.is_some()
                || s.mixing_metric.is_some()
                || s.spectral_guard.is_some()
                || s.locking.is_some()
                || s.preconditioner.is_some()
                || s.davidson_capacity_per_band.is_some()
                || s.start_density.as_deref() == Some("uniform")
                || s.start_wavefunctions
                    .as_deref()
                    .is_some_and(|v| v != "atomic")
                || s.energy_tolerance_ha.is_some()
                || s.density_residual_tolerance.is_some()
                || s.eigen_tolerance.is_some()
                || s.standard_stall_guard.is_some()
                || s.standard_density_forcing.is_some()
                // "off" is the reference behaviour; only an explicit legacy ladder is legacy-only.
                || s.cutoff_ladder.as_deref() == Some("on")
                || s.cutoff_ladder_fraction.is_some()
        })
}
fn l10(c: &Config, _: Mode, _: Command, _: &[&Pseudopotential]) -> bool {
    c.compat.reference() || c.xc.pw92_prefactor.as_deref() == Some("original")
}
fn l11(_: &Config, _: Mode, _: Command, upfs: &[&Pseudopotential]) -> bool {
    upfs.iter().any(|p| {
        p.atomic_wavefunctions()
            .is_none_or(|w| w.is_empty() || w.iter().any(|c| c.l > 3))
            || p.raw().rhoatom.values.len() != p.mesh().r().len()
            || p.raw().rhoatom.values.iter().all(|&v| v == 0.)
    })
}
type Predicate = fn(&Config, Mode, Command, &[&Pseudopotential]) -> bool;
const LEGACY_ONLY: &[(Predicate, &str)] = &[
    (l1, G4D_FALLBACK),
    (l2, G4B_FALLBACK),
    (l4, VARIABLE_CELL),
    (l5, EXTRAPOLATION),
    (l6, CHECKPOINT),
    (l7, RUNTIME),
    (l8, SESSION),
    (l9, LEGACY_INPUT),
    (l10, COMPAT),
    (l11, ATOMIC_DATA),
];
/// An empty UPF slice is a preliminary input-only decision. Commands must pass
/// every parsed potential for the final decision before engine construction.
pub fn select(
    config: &Config,
    mode: Mode,
    command: Command,
    upfs: &[&Pseudopotential],
) -> Selection {
    let mut selection = requested(config);
    // The diagnostic's explicit --solver choice is an input request. STRICT
    // still wins; subsequent eligibility checks are shared with every command.
    if let Command::RefscfMap {
        reference: Some(reference),
    } = command
        && !selection.strict()
    {
        selection = Selection {
            driver: if reference {
                Driver::Reference
            } else {
                Driver::LegacyStandard
            },
            reason: None,
            source: Source::Input,
        };
    }
    if selection.reference() {
        if let Some((_, reason)) = LEGACY_ONLY
            .iter()
            .find(|(predicate, _)| predicate(config, mode, command, upfs))
        {
            selection.driver = Driver::LegacyStandard;
            selection.reason = Some(reason);
        }
        if spin_fd(config) {
            RunConfig::current().note_path(
                Id::REFSCF_SPIN_FD_MOMENT,
                if selection.reference() {
                    "reference"
                } else {
                    "legacy"
                },
                selection.reason,
            );
        }
    }
    selection
}
pub fn preliminary(config: &Config) -> Selection {
    select(config, Mode::Classic, Command::Run, &[])
}
