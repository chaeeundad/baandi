//! Shared ESLE/checkpoint data, independent of serialization and storage.
//!
//! ModelIdentity and ExecutionIdentity are supplied by the trusted host which
//! constructed the problem. In particular, the engine cannot hash the original
//! UPF files or attest the executable. The CLI computes those hashes itself.
use crate::energy::EnergyAudit;
use crate::events::ExecutionConfig;
use crate::runtime::{AuditRecord, RuntimeConfig, RuntimeState};
use crate::{CompletionMode, IterationRecord, ScfOptions};
pub use dft_solvers::kerker::{MixerState, MixingDiagnostics};
pub use dft_solvers::mixing::PulayState;
pub use dft_solvers::{HybridRmm, RmmDiis, SolverChoice};

pub const SCHEMA_VERSION: u32 = 1;
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CheckpointKind {
    ExactContinuation,
    WarmStart,
    ResultExport,
}
#[derive(Clone, Debug, PartialEq)]
pub struct ModelIdentity {
    pub cell_bohr: [[f64; 3]; 3],
    pub atoms_fractional: Vec<[f64; 3]>,
    pub species: Vec<String>,
    pub pseudopotential_sha256: Vec<String>,
    pub cutoff_hartree: f64,
    pub kpoints: Vec<[f64; 4]>,
    pub electrons: f64,
    pub xc: String,
    pub spin: String,
    pub environment: String,
}
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct ExecutionIdentity {
    /// SHA-256 of the actual running executable (includes linked numerical code).
    pub build_sha256: String,
    /// Host-supplied CPU/OS/backend identity; exact import requires equality.
    pub platform: String,
    pub threads: usize,
}
#[derive(Clone, Debug)]
pub struct Orbitals {
    pub rows: usize,
    pub cols: usize,
    /// Column major, interleaved real/imaginary f64 coefficients.
    pub coefficients: Vec<f64>,
}
#[derive(Clone, Debug)]
pub struct ElectronicState {
    pub grid: [usize; 3],
    /// C order, last coordinate contiguous, electrons/Bohr^3.
    pub density: Vec<f64>,
    pub orbitals: Vec<Orbitals>,
    pub occupations: Vec<f64>,
}
#[derive(Clone, Debug)]
pub struct ContinuationState {
    pub atomic_start: Option<crate::AtomicStartReport>,
    pub options: ScfOptions,
    pub completion: CompletionMode,
    pub config: ExecutionConfig,
    pub iteration: usize,
    pub state_version: u64,
    pub sequence: u64,
    pub mixer: MixerState,
    pub history: Vec<IterationRecord>,
    pub best: Option<f64>,
    pub growing: usize,
    pub mixing_resets: usize,
    pub minimum_gap: Option<f64>,
}
#[derive(Clone, Debug)]
pub struct RuntimeCheckpoint {
    pub state: RuntimeState,
    pub revision: u64,
    pub config: RuntimeConfig,
    pub audit: Vec<AuditRecord>,
    pub autonomous_total: usize,
    pub autonomy_enabled: bool,
    pub last_policy_iteration: Option<usize>,
    pub rule_counts: Vec<usize>,
    pub rule_last_iteration: Vec<Option<usize>>,
}
#[derive(Clone, Debug)]
pub struct ExportResult {
    pub converged: bool,
    pub energy: EnergyAudit,
    pub charge: f64,
    pub gap_hartree: Option<f64>,
    pub eigenvalues: Vec<Vec<f64>>,
    pub residuals: Vec<Vec<f64>>,
    pub history: Vec<IterationRecord>,
}
#[derive(Clone, Debug)]
pub struct Checkpoint {
    pub standard_stall_guard: bool,
    pub standard_density_forcing: bool,
    pub standard: Option<crate::StandardProfile>,
    pub hartree_mixing: bool,
    pub schema_version: u32,
    pub kind: CheckpointKind,
    pub calc_id: String,
    pub sequence: u64,
    pub timestamp: std::time::SystemTime,
    pub model: ModelIdentity,
    pub execution: ExecutionIdentity,
    pub electronic_state: ElectronicState,
    pub continuation: Option<ContinuationState>,
    pub runtime: Option<RuntimeCheckpoint>,
    pub result: Option<ExportResult>,
}
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Rejection {
    SchemaVersion,
    Kind,
    Cell,
    Atoms,
    Species,
    Pseudopotential,
    Cutoff,
    KMesh,
    ElectronCount,
    Xc,
    Spin,
    Environment,
    Build,
    Threads,
    Platform,
    Options,
    Grid,
    Orbitals,
    Occupations,
    Density,
    History,
    Mixer,
    Boundary,
    UndrainedEvents,
    PendingControls,
    ExternalSink,
    Runtime,
    Identity,
}
impl std::fmt::Display for Rejection {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "checkpoint rejected: {self:?}")
    }
}
impl std::error::Error for Rejection {}
impl ModelIdentity {
    pub fn validate(&self) -> Result<(), Rejection> {
        if self.atoms_fractional.is_empty()
            || self.atoms_fractional.len() != self.species.len()
            || self.species.len() != self.pseudopotential_sha256.len()
            || self
                .cell_bohr
                .iter()
                .flatten()
                .chain(self.atoms_fractional.iter().flatten())
                .chain(self.kpoints.iter().flatten())
                .any(|v| !v.is_finite())
            || self.kpoints.is_empty()
            || self.kpoints.iter().any(|k| k[3] <= 0.0)
            || !self.cutoff_hartree.is_finite()
            || self.cutoff_hartree <= 0.0
            || !self.electrons.is_finite()
            || self.electrons <= 0.0
            || self.species.iter().any(String::is_empty)
            || self.pseudopotential_sha256.iter().any(|s| !is_digest(s))
        {
            return Err(Rejection::Identity);
        }
        Ok(())
    }
    pub fn compatible(&self, target: &Self) -> Result<(), Rejection> {
        self.validate()?;
        target.validate()?;
        macro_rules! same {
            ($field:ident, $why:ident) => {
                if self.$field != target.$field {
                    return Err(Rejection::$why);
                }
            };
        }
        if self
            .cell_bohr
            .iter()
            .flatten()
            .map(|v| v.to_bits())
            .ne(target.cell_bohr.iter().flatten().map(|v| v.to_bits()))
        {
            return Err(Rejection::Cell);
        }
        if self
            .atoms_fractional
            .iter()
            .flatten()
            .map(|v| v.to_bits())
            .ne(target
                .atoms_fractional
                .iter()
                .flatten()
                .map(|v| v.to_bits()))
        {
            return Err(Rejection::Atoms);
        }
        same!(species, Species);
        same!(pseudopotential_sha256, Pseudopotential);
        if self.cutoff_hartree.to_bits() != target.cutoff_hartree.to_bits() {
            return Err(Rejection::Cutoff);
        }
        if self.kpoints.iter().flatten().map(|v| v.to_bits()).ne(target
            .kpoints
            .iter()
            .flatten()
            .map(|v| v.to_bits()))
        {
            return Err(Rejection::KMesh);
        }
        if self.electrons.to_bits() != target.electrons.to_bits() {
            return Err(Rejection::ElectronCount);
        }
        same!(xc, Xc);
        same!(spin, Spin);
        same!(environment, Environment);
        Ok(())
    }
}
pub fn is_digest(s: &str) -> bool {
    s.len() == 64
        && s.bytes()
            .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
}
impl Checkpoint {
    pub fn compatible(
        &self,
        kind: CheckpointKind,
        model: &ModelIdentity,
        execution: &ExecutionIdentity,
    ) -> Result<(), Rejection> {
        if self.schema_version != SCHEMA_VERSION {
            return Err(Rejection::SchemaVersion);
        }
        if self.kind != kind {
            return Err(Rejection::Kind);
        }
        self.model.compatible(model)?;
        if !is_digest(&self.execution.build_sha256) || !is_digest(&execution.build_sha256) {
            return Err(Rejection::Identity);
        }
        if kind == CheckpointKind::ExactContinuation {
            if self.execution.build_sha256 != execution.build_sha256 {
                return Err(Rejection::Build);
            }
            if self.execution.platform != execution.platform {
                return Err(Rejection::Platform);
            }
            if self.execution.threads != execution.threads {
                return Err(Rejection::Threads);
            }
        }
        Ok(())
    }
}

/// Explicit bit comparison of every option, including signed floating zero.
/// No formatting or serialization participates in numerical compatibility.
pub(crate) fn options_identical(a: &ScfOptions, b: &ScfOptions) -> bool {
    if a.strain_reference
        .as_ref()
        .map(|c| c.lattice_columns_bohr())
        != b.strain_reference
            .as_ref()
            .map(|c| c.lattice_columns_bohr())
    {
        return false;
    }
    use crate::Mixing;
    let mixing = match (&a.mixing, &b.mixing) {
        (Mixing::Linear { beta: a }, Mixing::Linear { beta: b }) => a.to_bits() == b.to_bits(),
        (
            Mixing::Pulay {
                beta: a,
                history: x,
            },
            Mixing::Pulay {
                beta: b,
                history: y,
            },
        ) => a.to_bits() == b.to_bits() && x == y,
        (
            Mixing::KerkerPulay {
                beta: a,
                history: x,
                q0_inverse_bohr: p,
            },
            Mixing::KerkerPulay {
                beta: b,
                history: y,
                q0_inverse_bohr: q,
            },
        ) => a.to_bits() == b.to_bits() && x == y && p.to_bits() == q.to_bits(),
        _ => false,
    };
    let solver = match (&a.eigensolver, &b.eigensolver) {
        (SolverChoice::Davidson, SolverChoice::Davidson) => true,
        (SolverChoice::Hybrid(a), SolverChoice::Hybrid(b)) => {
            a.warmup_tolerance.to_bits() == b.warmup_tolerance.to_bits()
                && a.refinement.max_sweeps == b.refinement.max_sweeps
                && a.refinement.history_size == b.refinement.history_size
                && a.refinement.stagnation_sweeps == b.refinement.stagnation_sweeps
        }
        _ => false,
    };
    macro_rules! eq { ($($f:ident),*) => { true $(&& a.$f == b.$f)* }; }
    macro_rules! bits { ($($f:ident),*) => { true $(&& a.$f.to_bits() == b.$f.to_bits())* }; }
    mixing
        && a.crystal_density_group == b.crystal_density_group
        && a.standard == b.standard
        && a.standard_stall_guard == b.standard_stall_guard
        && a.standard_density_forcing == b.standard_density_forcing
        && a.hartree_mixing == b.hartree_mixing
        && solver
        && eq!(
            compat_mode,
            gamma_real,
            spectral_guard,
            locking,
            threads,
            functional,
            pw92_prefactor,
            local_mesh,
            radial_even_rule,
            allow_functional_override,
            bands,
            extra_bands,
            grid_override,
            residual_eigen_schedule,
            free_electron_initial_guess,
            start_wavefunctions,
            max_iterations,
            davidson_max_iterations,
            davidson_subspace,
            davidson_capacity_per_band,
            seed,
            divergence_patience
        )
        && bits!(
            cutoff_hartree,
            density_tolerance,
            energy_tolerance,
            eigen_tolerance,
            initial_eigen_tolerance,
            eigen_tolerance_decay,
            divergence_factor
        )
}

#[cfg(test)]
mod density_group_tests {
    use super::*;
    #[test]
    fn exact_restart_compares_density_group() {
        let mesh = ScfOptions::default();
        let crystal = ScfOptions {
            crystal_density_group: true,
            ..mesh.clone()
        };
        assert!(options_identical(&mesh, &mesh));
        assert!(!options_identical(&mesh, &crystal));
        assert!(!options_identical(&crystal, &mesh));
    }
}
