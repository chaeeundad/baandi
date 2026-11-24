//! Opt-in map-boundary experiment. No transport, serialization or environment I/O.
use super::record::Record;
use crate::ScfError;

#[derive(Clone, Debug)]
pub struct Action {
    pub beta: f64,
    pub reset_history: bool,
    pub controller_status: String,
    pub controller_wall_s: f64,
    pub exploration_seed: Option<u64>,
    pub exploration_draw: Option<u64>,
}
impl Action {
    pub fn fixed(beta: f64) -> Self {
        Self {
            beta,
            reset_history: false,
            controller_status: "off".into(),
            controller_wall_s: 0.,
            exploration_seed: None,
            exploration_draw: None,
        }
    }
}

#[derive(Clone, Debug)]
pub struct MapExperiment {
    pub action: Action,
    pub action_applied: bool,
    pub estimator_stop: bool,
    pub verification_status: &'static str,
    pub verification_error_ha: Option<f64>,
    pub verification_h_columns: usize,
    pub verification_davidson_iterations: usize,
    pub reused_from: Option<usize>,
    pub mixing_error: Option<String>,
}

/// Called outside the numerical kernels. `choose` sees the unmixed map and
/// verification verdict; `complete` sees the action and mixed-energy correction.
/// Terminal maps are observed, but their action is never applied.
pub trait MapHook {
    fn choose(&mut self, observation: &Record) -> Action;
    fn complete(&mut self, record: &Record) -> Result<(), ScfError>;
}

/// SplitMix64, fixed wrapping arithmetic and top-53-bit inclusive conversion.
/// This is reproducible exploration, not a cryptographic generator.
pub struct Exploration {
    seed: u64,
    state: u64,
}
impl Exploration {
    pub fn new(seed: u64) -> Self {
        Self { seed, state: seed }
    }
    pub fn draw(&mut self) -> Action {
        self.state = self.state.wrapping_add(0x9e3779b97f4a7c15);
        let mut z = self.state;
        z = (z ^ (z >> 30)).wrapping_mul(0xbf58476d1ce4e5b9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94d049bb133111eb);
        z ^= z >> 31;
        let mut action = Action::fixed(0.2 + 0.8 * ((z >> 11) as f64 / ((1_u64 << 53) - 1) as f64));
        action.controller_status = "exploration".into();
        action.exploration_seed = Some(self.seed);
        action.exploration_draw = Some(z);
        action
    }
}
