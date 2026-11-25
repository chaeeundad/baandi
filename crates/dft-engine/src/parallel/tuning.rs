//! Shared engine-side distributed initialization switches.
/// Engine-side M2 switches of the coarse stage and the certificate
/// (docs/performance/mpi-coarse-certificate-analysis.md; both default off).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct EngineTuning {
    /// `BAANDI_M2_SCREEN=dist`: certificate-screen clusters from the distributed
    /// band values (plane-local sums, one tree_sum per cluster).
    pub screen_dist: bool,
    /// `BAANDI_M2_ATOMIC=split`: atomic-start span, Gram-Schmidt first passes,
    /// H Q and Gram columns split over ranks.
    pub atomic_split: bool,
    /// `BAANDI_M2_ATOMIC_BLOCK=on` (default off, rounding class): with
    /// `atomic_split`, the span is filled by row blocks and orthonormalized by
    /// a row-distributed CholQR2 instead of the column-cyclic two-pass
    /// Gram-Schmidt (whose second pass is a serial chain over the columns).
    pub atomic_block: bool,
}
impl EngineTuning {
    pub fn from_env() -> Self {
        use dft_numerics::switches::{Id, choice};
        Self {
            screen_dist: choice(Id::M2_SCREEN) == "dist",
            atomic_split: choice(Id::M2_ATOMIC) == "split",
            atomic_block: choice(Id::M2_ATOMIC_BLOCK) == "on",
        }
    }
}
thread_local! {
    static TUNING_OVERRIDE: std::cell::Cell<Option<EngineTuning>> = const { std::cell::Cell::new(None) };
}
/// Override the engine M2 switches for the calling thread (simulated ranks in
/// tests); `None` restores the environment.
pub fn set_tuning_override(t: Option<EngineTuning>) {
    TUNING_OVERRIDE.with(|o| o.set(t));
}
/// The engine M2 switches in force on this thread.
pub fn tuning() -> EngineTuning {
    TUNING_OVERRIDE.with(|o| o.get()).unwrap_or_else(|| {
        static ENV: std::sync::OnceLock<EngineTuning> = std::sync::OnceLock::new();
        *ENV.get_or_init(EngineTuning::from_env)
    })
}
