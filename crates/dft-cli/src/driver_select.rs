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
    /// The G4D route predates R1 and already accepted reference mixer
    /// aliases/defaults before falling back. Keep that input contract intact.
    pub(crate) fn reference_mixing(self) -> bool {
        self.reference() || matches!(self.reason, Some(G4D_FALLBACK))
    }
    pub fn validate(self, _config: &Config) -> Result<(), crate::config::ConfigError> {
        Ok(())
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
        if self.reason == Some(G4D_FALLBACK) {
            result["reference_driver_fallback_reason"] = json!(self.reason);
            result["execution"]["reference_driver_fallback_reason"] = json!(self.reason);
        }
    }
    pub fn announce(self) {
        eprintln!("SCF driver: {}", self.json());
    }
}

pub const STRICT: &str = "BAANDI_REFERENCE_STRICT=on (strict reference path)";
pub const G4D_FALLBACK: &str =
    "reference driver: free-moment collinear FD with PBE not yet qualified, G4D";
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
    Selection {
        driver: Driver::Reference,
        reason: None,
        source: if config
            .scf
            .as_ref()
            .and_then(|s| s.driver.as_deref())
            .is_some()
        {
            Source::Input
        } else {
            Source::Default
        },
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
                Driver::StrictOracle
            },
            reason: None,
            source: Source::Input,
        };
    }
    if selection.reference() {
        if l1(config, mode, command, upfs) {
            selection.driver = Driver::StrictOracle;
            selection.reason = Some(G4D_FALLBACK);
        }
        if spin_fd(config) {
            RunConfig::current().note_path(
                Id::REFSCF_SPIN_FD_MOMENT,
                if selection.reference() {
                    "reference"
                } else {
                    "strict_oracle"
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
