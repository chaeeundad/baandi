//! Host-owned inputs and one thread per continuation: reference solver scopes
//! must never be interleaved on a thread. Only scalar records cross the channel.
use dft_engine::{
    ScfOptions, SpeciesPotential,
    input::{Atom, Cell, Pseudopotential},
    refscf::{
        api,
        runtime::{Control, PolicyPatch, Runtime, Snapshot},
    },
    refstart::StartInput,
};
use std::sync::mpsc;

/// Parsed model registered by a trusted host. The engine validates the reference
/// start at creation; tool callers cannot supply paths or model changes.
pub struct ReferenceInput {
    pub cell: Cell,
    pub atoms: Vec<Atom>,
    pub potentials: Vec<(String, Pseudopotential)>,
    pub options: ScfOptions,
    pub request: api::Request,
}
#[derive(Clone, Copy, Debug)]
pub(crate) enum Action {
    Step,
    Control(Control),
    Policy(PolicyPatch),
}
pub(crate) struct Worker {
    commands: mpsc::Sender<Action>,
    replies: mpsc::Receiver<(Snapshot, Result<(), String>)>,
    pub snapshot: Snapshot,
}
impl Worker {
    pub fn new(input: ReferenceInput) -> Result<Self, String> {
        let (commands, receive) = mpsc::channel();
        let (send, replies) = mpsc::channel();
        std::thread::Builder::new()
            .name("baandi-mcp-reference".into())
            .spawn(move || {
                let refs: Vec<_> = input
                    .potentials
                    .iter()
                    .map(|(label, potential)| SpeciesPotential { label, potential })
                    .collect();
                let mut runtime = Runtime::new(api::stepper(
                    StartInput {
                        cell: &input.cell,
                        atoms: &input.atoms,
                        potentials: &refs,
                        options: &input.options,
                    },
                    input.request,
                ));
                let result = runtime.initialize().map_err(|e| e.to_string());
                if send.send((runtime.snapshot(), result)).is_err() {
                    return;
                }
                for action in receive {
                    let result = match action {
                        Action::Step => runtime.step().map(|_| ()),
                        Action::Control(control) => runtime.control(control),
                        Action::Policy(patch) => runtime.apply_policy(patch).map(|_| ()),
                    }
                    .map_err(|e| e.to_string());
                    if send.send((runtime.snapshot(), result)).is_err() {
                        break;
                    }
                }
            })
            .map_err(|e| e.to_string())?;
        let (snapshot, result) = replies.recv().map_err(|e| e.to_string())?;
        result?;
        Ok(Self {
            commands,
            replies,
            snapshot,
        })
    }
    pub fn apply(&mut self, action: Action) -> Result<(), String> {
        let reply = self
            .commands
            .send(action)
            .map_err(|_| "reference worker stopped")
            .and_then(|()| self.replies.recv().map_err(|_| "reference worker failed"));
        match reply {
            Ok((snapshot, result)) => {
                self.snapshot = snapshot;
                result
            }
            Err(message) => {
                self.snapshot.state = dft_engine::runtime::RuntimeState::Failed;
                self.snapshot.state_version += 1;
                self.snapshot.safe_point = false;
                Err(message.into())
            }
        }
    }
}
