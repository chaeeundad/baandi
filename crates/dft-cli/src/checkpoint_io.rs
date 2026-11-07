//! Checkpoint v1 container: 8-byte magic, u64 LE header length, u64 LE binary
//! length, UTF-8 JSON header, contiguous f64 LE section, SHA-256 of all preceding
//! bytes. Arrays have byte offsets/counts; complex orbitals are column-major
//! interleaved real/imaginary, densities are C-order. Scalar floats use IEEE bits
//! as JSON u64 so decoding never depends on decimal floating-point round trips.
//! The same envelope carries ESLE scientific state and exact continuation.
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    fs,
    io::{Read, Write},
    path::Path,
};

pub type Error = Box<dyn std::error::Error>;
const MAGIC: &[u8; 8] = b"BAANDICP";
const MAX_FILE: usize = 512 * 1024 * 1024;
const MAX_HEADER: usize = 32 * 1024 * 1024;
fn invalid(s: &str) -> Error {
    std::io::Error::new(std::io::ErrorKind::InvalidData, s).into()
}
trait Wire: Sized {
    type Repr: Serialize + serde::de::DeserializeOwned;
    fn encode(&self, bytes: &mut Vec<u8>) -> Result<Self::Repr, Error>;
    fn decode(w: Self::Repr, bytes: &[u8], cursor: &mut usize) -> Result<Self, Error>;
}
macro_rules! scalar {
    ($($t:ty),*) => { $(impl Wire for $t {
        type Repr = Self;
        fn encode(&self, _: &mut Vec<u8>) -> Result<Self, Error> { Ok(self.clone()) }
        fn decode(w: Self, _: &[u8], _: &mut usize) -> Result<Self, Error> { Ok(w) }
    })* };
}
scalar!(u64, u32, usize, bool, String);
impl Wire for f64 {
    type Repr = u64;
    fn encode(&self, _: &mut Vec<u8>) -> Result<u64, Error> {
        if !self.is_finite() {
            return Err(invalid("nonfinite scalar"));
        }
        Ok(self.to_bits())
    }
    fn decode(w: u64, _: &[u8], _: &mut usize) -> Result<Self, Error> {
        let v = Self::from_bits(w);
        if !v.is_finite() {
            return Err(invalid("nonfinite scalar"));
        }
        Ok(v)
    }
}
impl<T: Wire> Wire for Option<T> {
    type Repr = Option<T::Repr>;
    fn encode(&self, b: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        self.as_ref().map(|v| v.encode(b)).transpose()
    }
    fn decode(w: Self::Repr, b: &[u8], c: &mut usize) -> Result<Self, Error> {
        w.map(|v| T::decode(v, b, c)).transpose()
    }
}
impl<T: Wire, const N: usize> Wire for [T; N] {
    type Repr = Vec<T::Repr>;
    fn encode(&self, b: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        self.iter().map(|v| v.encode(b)).collect()
    }
    fn decode(w: Self::Repr, b: &[u8], c: &mut usize) -> Result<Self, Error> {
        if w.len() != N {
            return Err(invalid("fixed array length"));
        }
        w.into_iter()
            .map(|v| T::decode(v, b, c))
            .collect::<Result<Vec<_>, _>>()?
            .try_into()
            .map_err(|_| invalid("fixed array length"))
    }
}
macro_rules! wire_vec {
    ($t:ty) => {
        impl Wire for Vec<$t> {
            type Repr = Vec<<$t as Wire>::Repr>;
            fn encode(&self, b: &mut Vec<u8>) -> Result<Self::Repr, Error> {
                self.iter().map(|v| v.encode(b)).collect()
            }
            fn decode(w: Self::Repr, b: &[u8], c: &mut usize) -> Result<Self, Error> {
                w.into_iter().map(|v| <$t>::decode(v, b, c)).collect()
            }
        }
    };
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct ArrayRef {
    offset: usize,
    count: usize,
}
impl Wire for Vec<f64> {
    type Repr = ArrayRef;
    fn encode(&self, bytes: &mut Vec<u8>) -> Result<ArrayRef, Error> {
        let offset = bytes.len();
        let length = self
            .len()
            .checked_mul(8)
            .and_then(|v| v.checked_add(offset))
            .ok_or_else(|| invalid("array overflow"))?;
        if length > MAX_FILE {
            return Err(invalid("array quota"));
        }
        for v in self {
            if !v.is_finite() {
                return Err(invalid("nonfinite array"));
            }
            bytes.extend_from_slice(&v.to_le_bytes());
        }
        Ok(ArrayRef {
            offset,
            count: self.len(),
        })
    }
    fn decode(w: ArrayRef, bytes: &[u8], cursor: &mut usize) -> Result<Self, Error> {
        let end = w
            .count
            .checked_mul(8)
            .and_then(|v| v.checked_add(w.offset))
            .ok_or_else(|| invalid("array overflow"))?;
        if w.offset != *cursor || end > bytes.len() {
            return Err(invalid("array offset/length"));
        }
        let result = bytes[w.offset..end]
            .as_chunks::<8>()
            .0
            .iter()
            .map(|v| {
                let x = f64::from_le_bytes(*v);
                if x.is_finite() {
                    Ok(x)
                } else {
                    Err(invalid("nonfinite array"))
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        *cursor = end;
        Ok(result)
    }
}
impl Wire for Duration {
    type Repr = (u64, u32);
    fn encode(&self, _: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok((self.as_secs(), self.subsec_nanos()))
    }
    fn decode(w: Self::Repr, _: &[u8], _: &mut usize) -> Result<Self, Error> {
        if w.1 >= 1_000_000_000 {
            return Err(invalid("duration nanoseconds"));
        }
        Ok(Duration::new(w.0, w.1))
    }
}
impl Wire for SystemTime {
    type Repr = (bool, <Duration as Wire>::Repr);
    fn encode(&self, b: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self.duration_since(UNIX_EPOCH) {
            Ok(d) => (true, d.encode(b)?),
            Err(e) => (false, e.duration().encode(b)?),
        })
    }
    fn decode(w: Self::Repr, b: &[u8], c: &mut usize) -> Result<Self, Error> {
        let d = Duration::decode(w.1, b, c)?;
        (if w.0 {
            UNIX_EPOCH.checked_add(d)
        } else {
            UNIX_EPOCH.checked_sub(d)
        })
        .ok_or_else(|| invalid("timestamp overflow"))
    }
}
macro_rules! wire_struct {
    ($t:ident, $w:ident { $($(#[$attr:meta])* $f:ident: $ty:ty,)* } $(defaults { $($df:ident: $dv:expr,)* })?) => {
        #[derive(Serialize, Deserialize)]
        #[serde(deny_unknown_fields)]
        struct $w { $($(#[$attr])* $f: <$ty as Wire>::Repr,)* }
        impl Wire for $t {
            type Repr = $w;
            fn encode(&self, bytes: &mut Vec<u8>) -> Result<$w, Error> { Ok($w { $($f: self.$f.encode(bytes)?,)* }) }
            fn decode(w: $w, bytes: &[u8], cursor: &mut usize) -> Result<Self, Error> {
                Ok(Self { $($f: <$ty>::decode(w.$f, bytes, cursor)?,)* $($($df: $dv,)*)? })
            }
        }
    };
}

/// Encoding never opens files and can be used by an authenticated host adapter.
pub fn encode(checkpoint: &Checkpoint) -> Result<Vec<u8>, Error> {
    if checkpoint
        .continuation
        .as_ref()
        .is_some_and(|c| c.options.strain_reference.is_some())
    {
        return Err(invalid(
            "checkpoint v1 cannot preserve frozen-Miller strain reference",
        ));
    }
    if checkpoint.continuation.as_ref().is_some_and(|c| {
        c.options.spectral_guard == dft_engine::SpectralGuard::Final
            || c.options.locking != dft_engine::Locking::Off
    }) {
        return Err(invalid(
            "checkpoint v1 cannot preserve scf.spectral_guard=final audit history or scf.locking=auto policy",
        ));
    }
    let mut binary = Vec::new();
    let wire = checkpoint.encode(&mut binary)?;
    let header = serde_json::to_vec(&wire)?;
    if header.len() > MAX_HEADER || header.len() + binary.len() + 56 > MAX_FILE {
        return Err(invalid("checkpoint quota"));
    }
    let mut out = MAGIC.to_vec();
    out.extend_from_slice(&(header.len() as u64).to_le_bytes());
    out.extend_from_slice(&(binary.len() as u64).to_le_bytes());
    out.extend_from_slice(&header);
    out.extend_from_slice(&binary);
    let digest = Sha256::digest(&out);
    out.extend_from_slice(&digest);
    Ok(out)
}
pub fn decode(bytes: &[u8]) -> Result<Checkpoint, Error> {
    if bytes.len() < 56 || bytes.len() > MAX_FILE || &bytes[..8] != MAGIC {
        return Err(invalid("checkpoint magic/size"));
    }
    let head = usize::try_from(u64::from_le_bytes(bytes[8..16].try_into().unwrap()))?;
    let body = usize::try_from(u64::from_le_bytes(bytes[16..24].try_into().unwrap()))?;
    if head > MAX_HEADER
        || head.checked_add(body).and_then(|v| v.checked_add(56)) != Some(bytes.len())
    {
        return Err(invalid("checkpoint lengths"));
    }
    let payload_end = bytes.len() - 32;
    if Sha256::digest(&bytes[..payload_end])[..] != bytes[payload_end..] {
        return Err(invalid("checkpoint checksum"));
    }
    // Derived struct visitors reject unknown and duplicate fields.
    let wire: WCheckpoint = serde_json::from_slice(&bytes[24..24 + head])?;
    let mut cursor = 0;
    let cp = Checkpoint::decode(wire, &bytes[24 + head..payload_end], &mut cursor)?;
    if cursor != body {
        return Err(invalid("unreferenced binary data"));
    }
    if cp.schema_version != dft_engine::checkpoint::SCHEMA_VERSION {
        return Err(invalid("SchemaVersion"));
    }
    cp.model.validate()?;
    Ok(cp)
}
pub fn read(path: &Path) -> Result<Checkpoint, Error> {
    let file = fs::File::open(path)?;
    if !file.metadata()?.is_file() || file.metadata()?.len() > MAX_FILE as u64 {
        return Err(invalid("checkpoint file quota/type"));
    }
    let mut bytes = Vec::new();
    file.take(MAX_FILE as u64 + 1).read_to_end(&mut bytes)?;
    decode(&bytes)
}
/// Publish a complete immutable file without replacing an existing destination.
/// A same-directory hard link is the atomic no-clobber publication operation.
pub fn write(path: &Path, cp: &Checkpoint) -> Result<(), Error> {
    let bytes = encode(cp)?;
    let parent = path
        .parent()
        .filter(|p| !p.as_os_str().is_empty())
        .unwrap_or(Path::new("."));
    let temp = parent.join(format!(
        ".checkpoint-{}-{}.tmp",
        std::process::id(),
        SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos()
    ));
    let mut file = fs::OpenOptions::new()
        .write(true)
        .create_new(true)
        .open(&temp)?;
    let result = (|| -> Result<(), Error> {
        file.write_all(&bytes)?;
        file.sync_all()?;
        fs::hard_link(&temp, path)?;
        fs::File::open(parent)?.sync_all()?;
        Ok(())
    })();
    let _ = fs::remove_file(&temp);
    result
}

// Adapter-owned schema v1 mirrors. No serde in the numerical core.
use dft_engine::CompletionMode;
use dft_engine::FixedPointCertificate;
use dft_engine::Functional;
use dft_engine::IterationRecord;
use dft_engine::Mixing;
use dft_engine::PhaseTimings;
use dft_engine::{ScfOptions, StandardProfile};
fn is_false(v: &bool) -> bool {
    !v
}
use dft_engine::checkpoint::Checkpoint;
use dft_engine::checkpoint::CheckpointKind;
use dft_engine::checkpoint::ContinuationState;
use dft_engine::checkpoint::ElectronicState;
use dft_engine::checkpoint::ExecutionIdentity;
use dft_engine::checkpoint::ExportResult;
use dft_engine::checkpoint::HybridRmm;
use dft_engine::checkpoint::MixerState;
use dft_engine::checkpoint::MixingDiagnostics;
use dft_engine::checkpoint::ModelIdentity;
use dft_engine::checkpoint::Orbitals;
use dft_engine::checkpoint::PulayState;
use dft_engine::checkpoint::RmmDiis;
use dft_engine::checkpoint::RuntimeCheckpoint;
use dft_engine::checkpoint::SolverChoice;
use dft_engine::energy::EnergyAudit;
use dft_engine::events::ExecutionConfig;
use dft_engine::events::ExecutionMode;
use dft_engine::policy::PolicyChange;
use dft_engine::policy::PolicyError;
use dft_engine::policy::RuntimeEigensolver;
use dft_engine::policy::ToleranceFloors;
use dft_engine::runtime::AllowedAction;
use dft_engine::runtime::Approval;
use dft_engine::runtime::AuditOutcome;
use dft_engine::runtime::AuditRecord;
use dft_engine::runtime::Authority;
use dft_engine::runtime::AutonomyPolicy;
use dft_engine::runtime::AutonomyRule;
use dft_engine::runtime::BoundarySnapshot;
use dft_engine::runtime::Comparison;
use dft_engine::runtime::ControlAction;
use dft_engine::runtime::ControlError;
use dft_engine::runtime::ControlRequest;
use dft_engine::runtime::Expiry;
use dft_engine::runtime::InvalidTransition;
use dft_engine::runtime::Precondition;
use dft_engine::runtime::RuntimeConfig;
use dft_engine::runtime::RuntimeState;
use dft_engine::runtime::Scalar;
use dft_engine::runtime::StalePolicy;
use dft_engine::runtime::Trigger;
use std::time::{Duration, SystemTime, UNIX_EPOCH};
wire_vec!(AuditRecord);
wire_vec!(AutonomyRule);
wire_vec!(IterationRecord);
wire_vec!(Option<usize>);
wire_vec!(Orbitals);
wire_vec!(PolicyChange);
wire_vec!(Precondition);
wire_vec!(String);
wire_vec!(Vec<f64>);
wire_vec!([f64; 3]);
wire_vec!([f64; 4]);
wire_vec!(usize);
wire_struct!(Checkpoint, WCheckpoint {
    #[serde(default, skip_serializing_if = "is_false")]
    standard_stall_guard: bool,
    #[serde(default, skip_serializing_if = "is_false")]
    standard_density_forcing: bool,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    standard: Option<StandardProfile>,
    #[serde(default, skip_serializing_if = "is_false")]
    hartree_mixing: bool,
    schema_version: u32,
    kind: CheckpointKind,
    calc_id: String,
    sequence: u64,
    timestamp: SystemTime,
    model: ModelIdentity,
    execution: ExecutionIdentity,
    electronic_state: ElectronicState,
    continuation: Option<ContinuationState>,
    runtime: Option<RuntimeCheckpoint>,
    result: Option<ExportResult>,
});
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WCheckpointKind {
    ExactContinuation,
    WarmStart,
    ResultExport,
}
impl Wire for CheckpointKind {
    type Repr = WCheckpointKind;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::ExactContinuation => WCheckpointKind::ExactContinuation,
            Self::WarmStart => WCheckpointKind::WarmStart,
            Self::ResultExport => WCheckpointKind::ResultExport,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WCheckpointKind::ExactContinuation => Self::ExactContinuation,
            WCheckpointKind::WarmStart => Self::WarmStart,
            WCheckpointKind::ResultExport => Self::ResultExport,
        })
    }
}
wire_struct!(ModelIdentity, WModelIdentity {
    cell_bohr: [[f64; 3]; 3],
    atoms_fractional: Vec<[f64; 3]>,
    species: Vec<String>,
    pseudopotential_sha256: Vec<String>,
    cutoff_hartree: f64,
    kpoints: Vec<[f64; 4]>,
    electrons: f64,
    xc: String,
    spin: String,
    environment: String,
});
wire_struct!(
    ExecutionIdentity,
    WExecutionIdentity {
        build_sha256: String,
        platform: String,
        threads: usize,
    }
);
wire_struct!(ElectronicState, WElectronicState {
    grid: [usize; 3],
    density: Vec<f64>,
    orbitals: Vec<Orbitals>,
    occupations: Vec<f64>,
});
wire_struct!(Orbitals, WOrbitals {
    rows: usize,
    cols: usize,
    coefficients: Vec<f64>,
});
wire_struct!(ContinuationState, WContinuationState {
    options: ScfOptions,
    completion: CompletionMode,
    config: ExecutionConfig,
    iteration: usize,
    state_version: u64,
    sequence: u64,
    mixer: MixerState,
    history: Vec<IterationRecord>,
    best: Option<f64>,
    growing: usize,
    mixing_resets: usize,
    minimum_gap: Option<f64>,
});
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct WStandardProfile {
    convergence_profile: String,
    certified: bool,
    threshold_ha: <f64 as Wire>::Repr,
}
impl Wire for StandardProfile {
    type Repr = WStandardProfile;
    fn encode(&self, bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(WStandardProfile {
            convergence_profile: "standard".into(),
            certified: false,
            threshold_ha: self.threshold_ha.encode(bytes)?,
        })
    }
    fn decode(w: Self::Repr, bytes: &[u8], cursor: &mut usize) -> Result<Self, Error> {
        let threshold_ha = f64::decode(w.threshold_ha, bytes, cursor)?;
        if w.convergence_profile != "standard"
            || w.certified
            || !threshold_ha.is_finite()
            || threshold_ha <= 0.0
        {
            return Err(invalid("invalid standard profile provenance"));
        }
        Ok(Self { threshold_ha })
    }
}
impl Wire for dft_engine::GammaReal {
    type Repr = String;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<String, Error> {
        Ok(self.name().into())
    }
    fn decode(w: String, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        match w.as_str() {
            "off" => Ok(Self::Off),
            "auto" => Ok(Self::Auto),
            "on" => Ok(Self::On),
            _ => Err(invalid("invalid gamma_real option")),
        }
    }
}
fn gamma_real_off() -> String {
    "off".into()
}
fn is_gamma_real_off(s: &str) -> bool {
    s == "off"
}
impl Wire for dft_engine::Pw92Prefactor {
    type Repr = String;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<String, Error> {
        Ok(self.name().into())
    }
    fn decode(w: String, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        match w.as_str() {
            "modified" => Ok(Self::Modified),
            "original" => Ok(Self::Original),
            _ => Err(invalid("invalid pw92_prefactor option")),
        }
    }
}
fn pw92_modified() -> String {
    "modified".into()
}
fn is_pw92_modified(s: &str) -> bool {
    s == "modified"
}
impl Wire for dft_engine::LocalMesh {
    type Repr = String;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<String, Error> {
        Ok(self.name().into())
    }
    fn decode(w: String, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        match w.as_str() {
            "baandi" => Ok(Self::Baandi),
            "qe_like" => Ok(Self::QeLike),
            _ => Err(invalid("invalid local_mesh option")),
        }
    }
}
fn default_local_mesh() -> String {
    "baandi".into()
}
fn is_default_local_mesh(s: &str) -> bool {
    s == "baandi"
}
impl Wire for dft_engine::RadialEvenRule {
    type Repr = String;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<String, Error> {
        Ok(self.name().into())
    }
    fn decode(w: String, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        match w.as_str() {
            "simpson_3_8" => Ok(Self::Simpson3_8),
            "simpson_quadratic_end" => Ok(Self::SimpsonQuadraticEnd),
            _ => Err(invalid("invalid radial_even_rule option")),
        }
    }
}
fn default_radial_even_rule() -> String {
    "simpson_3_8".into()
}
fn is_default_radial_even_rule(s: &str) -> bool {
    s == "simpson_3_8"
}

fn default_davidson_capacity_per_band() -> usize {
    6
}
fn is_default_davidson_capacity_per_band(value: &usize) -> bool {
    *value == 6
}
wire_struct!(ScfOptions, WScfOptions {
    #[serde(default = "default_davidson_capacity_per_band", skip_serializing_if = "is_default_davidson_capacity_per_band")]
    davidson_capacity_per_band: usize,
    #[serde(default = "default_radial_even_rule", skip_serializing_if = "is_default_radial_even_rule")]
    radial_even_rule: dft_engine::RadialEvenRule,
    #[serde(default = "default_local_mesh", skip_serializing_if = "is_default_local_mesh")]
    local_mesh: dft_engine::LocalMesh,
    #[serde(default, skip_serializing_if = "is_false")]
    standard_stall_guard: bool,
    #[serde(default, skip_serializing_if = "is_false")]
    standard_density_forcing: bool,
    #[serde(default = "pw92_modified", skip_serializing_if = "is_pw92_modified")]
    pw92_prefactor: dft_engine::Pw92Prefactor,
    #[serde(default = "gamma_real_off", skip_serializing_if = "is_gamma_real_off")]
    gamma_real: dft_engine::GammaReal,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    standard: Option<StandardProfile>,
    #[serde(default, skip_serializing_if = "is_false")]
    hartree_mixing: bool,
    eigensolver: SolverChoice,
    threads: usize,
    functional: Functional,
    allow_functional_override: bool,
    cutoff_hartree: f64,
    bands: Option<usize>,
    extra_bands: usize,
    grid_override: Option<[usize; 3]>,
    mixing: Mixing,
    density_tolerance: f64,
    energy_tolerance: f64,
    eigen_tolerance: f64,
    initial_eigen_tolerance: f64,
    eigen_tolerance_decay: f64,
    residual_eigen_schedule: bool,
    free_electron_initial_guess: bool,
    max_iterations: usize,
    davidson_max_iterations: usize,
    davidson_subspace: usize,
    seed: u64,
    divergence_factor: f64,
    divergence_patience: usize,
} defaults { locking: dft_engine::Locking::Off, spectral_guard: dft_engine::SpectralGuard::Off, strain_reference: None, });
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WSolverChoice {
    Davidson,
    Hybrid(<HybridRmm as Wire>::Repr),
}
impl Wire for SolverChoice {
    type Repr = WSolverChoice;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Davidson => WSolverChoice::Davidson,
            Self::Hybrid(v) => WSolverChoice::Hybrid(v.encode(_bytes)?),
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WSolverChoice::Davidson => Self::Davidson,
            WSolverChoice::Hybrid(v) => Self::Hybrid(Wire::decode(v, _bytes, _cursor)?),
        })
    }
}
wire_struct!(
    HybridRmm,
    WHybridRmm {
        refinement: RmmDiis,
        warmup_tolerance: f64,
    }
);
wire_struct!(
    RmmDiis,
    WRmmDiis {
        max_sweeps: usize,
        history_size: usize,
        stagnation_sweeps: usize,
    }
);
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WFunctional {
    LdaPz,
    LdaPw,
    Pbe,
}
impl Wire for Functional {
    type Repr = WFunctional;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::LdaPz => WFunctional::LdaPz,
            Self::LdaPw => WFunctional::LdaPw,
            Self::Pbe => WFunctional::Pbe,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WFunctional::LdaPz => Self::LdaPz,
            WFunctional::LdaPw => Self::LdaPw,
            WFunctional::Pbe => Self::Pbe,
        })
    }
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WMixing {
    Linear {
        beta: <f64 as Wire>::Repr,
    },
    Pulay {
        beta: <f64 as Wire>::Repr,
        history: <usize as Wire>::Repr,
    },
    KerkerPulay {
        beta: <f64 as Wire>::Repr,
        history: <usize as Wire>::Repr,
        q0_inverse_bohr: <f64 as Wire>::Repr,
    },
}
impl Wire for Mixing {
    type Repr = WMixing;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Linear { beta } => WMixing::Linear {
                beta: beta.encode(_bytes)?,
            },
            Self::Pulay { beta, history } => WMixing::Pulay {
                beta: beta.encode(_bytes)?,
                history: history.encode(_bytes)?,
            },
            Self::KerkerPulay {
                beta,
                history,
                q0_inverse_bohr,
            } => WMixing::KerkerPulay {
                beta: beta.encode(_bytes)?,
                history: history.encode(_bytes)?,
                q0_inverse_bohr: q0_inverse_bohr.encode(_bytes)?,
            },
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WMixing::Linear { beta } => Self::Linear {
                beta: Wire::decode(beta, _bytes, _cursor)?,
            },
            WMixing::Pulay { beta, history } => Self::Pulay {
                beta: Wire::decode(beta, _bytes, _cursor)?,
                history: Wire::decode(history, _bytes, _cursor)?,
            },
            WMixing::KerkerPulay {
                beta,
                history,
                q0_inverse_bohr,
            } => Self::KerkerPulay {
                beta: Wire::decode(beta, _bytes, _cursor)?,
                history: Wire::decode(history, _bytes, _cursor)?,
                q0_inverse_bohr: Wire::decode(q0_inverse_bohr, _bytes, _cursor)?,
            },
        })
    }
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WCompletionMode {
    Strict,
    RetainLast,
}
impl Wire for CompletionMode {
    type Repr = WCompletionMode;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Strict => WCompletionMode::Strict,
            Self::RetainLast => WCompletionMode::RetainLast,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WCompletionMode::Strict => Self::Strict,
            WCompletionMode::RetainLast => Self::RetainLast,
        })
    }
}
wire_struct!(
    ExecutionConfig,
    WExecutionConfig {
        mode: ExecutionMode,
        calc_id: String,
        routine_capacity: usize,
        critical_capacity: usize,
    }
);
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WExecutionMode {
    Classic,
    Observed,
    AgentControlled,
}
impl Wire for ExecutionMode {
    type Repr = WExecutionMode;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Classic => WExecutionMode::Classic,
            Self::Observed => WExecutionMode::Observed,
            Self::AgentControlled => WExecutionMode::AgentControlled,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WExecutionMode::Classic => Self::Classic,
            WExecutionMode::Observed => Self::Observed,
            WExecutionMode::AgentControlled => Self::AgentControlled,
        })
    }
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WMixerState {
    Linear {
        beta: <f64 as Wire>::Repr,
    },
    Pulay(<PulayState as Wire>::Repr),
    Kerker {
        pulay: <PulayState as Wire>::Repr,
        mean: <Option<f64> as Wire>::Repr,
        previous: <Option<f64> as Wire>::Repr,
        diagnostics: <MixingDiagnostics as Wire>::Repr,
    },
}
impl Wire for MixerState {
    type Repr = WMixerState;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Linear { beta } => WMixerState::Linear {
                beta: beta.encode(_bytes)?,
            },
            Self::Pulay(v) => WMixerState::Pulay(v.encode(_bytes)?),
            Self::Kerker {
                pulay,
                mean,
                previous,
                diagnostics,
            } => WMixerState::Kerker {
                pulay: pulay.encode(_bytes)?,
                mean: mean.encode(_bytes)?,
                previous: previous.encode(_bytes)?,
                diagnostics: diagnostics.encode(_bytes)?,
            },
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WMixerState::Linear { beta } => Self::Linear {
                beta: Wire::decode(beta, _bytes, _cursor)?,
            },
            WMixerState::Pulay(v) => Self::Pulay(Wire::decode(v, _bytes, _cursor)?),
            WMixerState::Kerker {
                pulay,
                mean,
                previous,
                diagnostics,
            } => Self::Kerker {
                pulay: Wire::decode(pulay, _bytes, _cursor)?,
                mean: Wire::decode(mean, _bytes, _cursor)?,
                previous: Wire::decode(previous, _bytes, _cursor)?,
                diagnostics: Wire::decode(diagnostics, _bytes, _cursor)?,
            },
        })
    }
}
wire_struct!(PulayState, WPulayState {
    beta: f64,
    capacity: usize,
    inputs: Vec<Vec<f64>>,
    residuals: Vec<Vec<f64>>,
    mean: Option<f64>,
});
wire_struct!(MixingDiagnostics, WMixingDiagnostics {
    residual_norm: f64,
    contraction: Option<f64>,
    history_used: usize,
    resets: usize,
    raw_recoveries: usize,
    used_raw_recovery: bool,
});
wire_struct!(IterationRecord, WIterationRecord {
    #[serde(default, skip_serializing_if = "Option::is_none")]
    standard_error_ha: Option<f64>,
    #[serde(default, skip_serializing_if = "is_false")]
    standard_redo: bool,
    iteration: usize,
    energy: EnergyAudit,
    energy_change: Option<f64>,
    density_residual: f64,
    max_eigen_residual: f64,
    requested_eigen_tolerance: f64,
    davidson_iterations: usize,
    operator_applications: usize,
    davidson_iterations_per_k: Vec<usize>,
    operator_applications_per_k: Vec<usize>,
    charge_error: f64,
    input_charge_error: f64,
    min_density: f64,
    mixing_fraction: f64,
    mixer_reset: bool,
    density_contraction: Option<f64>,
    mixing_history_used: Option<usize>,
    mixing_resets: usize,
    mixing_raw_recovery: bool,
    fixed_point_certificate: Option<FixedPointCertificate>,
    timings: PhaseTimings,
});
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct WEnergyAudit {
    kinetic: <f64 as Wire>::Repr,
    local: <f64 as Wire>::Repr,
    local_zero_mode: <f64 as Wire>::Repr,
    nonlocal: <f64 as Wire>::Repr,
    hartree: <f64 as Wire>::Repr,
    xc: <f64 as Wire>::Repr,
    ewald: <f64 as Wire>::Repr,
    band_sum: <f64 as Wire>::Repr,
    double_counting: <f64 as Wire>::Repr,
    band_form: <f64 as Wire>::Repr,
    harris_foulkes: <f64 as Wire>::Repr,
    functional_difference_bound: <f64 as Wire>::Repr,
}
impl Wire for EnergyAudit {
    type Repr = WEnergyAudit;
    fn encode(&self, bytes: &mut Vec<u8>) -> Result<WEnergyAudit, Error> {
        if self.nlcc.is_some() || self.xc_floor.is_some() {
            return Err(invalid("checkpoint v1 does not carry NLCC ledger"));
        }
        Ok(WEnergyAudit {
            kinetic: self.kinetic.encode(bytes)?,
            local: self.local.encode(bytes)?,
            local_zero_mode: self.local_zero_mode.encode(bytes)?,
            nonlocal: self.nonlocal.encode(bytes)?,
            hartree: self.hartree.encode(bytes)?,
            xc: self.xc.encode(bytes)?,
            ewald: self.ewald.encode(bytes)?,
            band_sum: self.band_sum.encode(bytes)?,
            double_counting: self.double_counting.encode(bytes)?,
            band_form: self.band_form.encode(bytes)?,
            harris_foulkes: self.harris_foulkes.encode(bytes)?,
            functional_difference_bound: self.functional_difference_bound.encode(bytes)?,
        })
    }
    fn decode(w: WEnergyAudit, bytes: &[u8], cursor: &mut usize) -> Result<Self, Error> {
        Ok(Self {
            nlcc: None,
            xc_floor: None,
            kinetic: f64::decode(w.kinetic, bytes, cursor)?,
            local: f64::decode(w.local, bytes, cursor)?,
            local_zero_mode: f64::decode(w.local_zero_mode, bytes, cursor)?,
            nonlocal: f64::decode(w.nonlocal, bytes, cursor)?,
            hartree: f64::decode(w.hartree, bytes, cursor)?,
            xc: f64::decode(w.xc, bytes, cursor)?,
            ewald: f64::decode(w.ewald, bytes, cursor)?,
            band_sum: f64::decode(w.band_sum, bytes, cursor)?,
            double_counting: f64::decode(w.double_counting, bytes, cursor)?,
            band_form: f64::decode(w.band_form, bytes, cursor)?,
            harris_foulkes: f64::decode(w.harris_foulkes, bytes, cursor)?,
            functional_difference_bound: f64::decode(w.functional_difference_bound, bytes, cursor)?,
        })
    }
}
wire_struct!(FixedPointCertificate, WFixedPointCertificate {
    strict_residuals_per_k: Vec<Vec<f64>>,
    buffer_residuals_per_k: Vec<Vec<f64>>,
    density_residual: f64,
    occupied_eigenvalue_change: f64,
    davidson_iterations: usize,
    operator_applications: usize,
    davidson_iterations_per_k: Vec<usize>,
    operator_applications_per_k: Vec<usize>,
    timings: PhaseTimings,
});
wire_struct!(
    PhaseTimings,
    WPhaseTimings {
        potential: Duration,
        diagonalization: Duration,
        density: Duration,
        energy: Duration,
        mixing: Duration,
    }
);
wire_struct!(RuntimeCheckpoint, WRuntimeCheckpoint {
    state: RuntimeState,
    revision: u64,
    config: RuntimeConfig,
    audit: Vec<AuditRecord>,
    autonomous_total: usize,
    autonomy_enabled: bool,
    last_policy_iteration: Option<usize>,
    rule_counts: Vec<usize>,
    rule_last_iteration: Vec<Option<usize>>,
});
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WRuntimeState {
    Created,
    Initialized,
    Running,
    PauseRequested,
    Paused,
    Restarting,
    Converged,
    Exhausted,
    Failed,
    Cancelled,
}
impl Wire for RuntimeState {
    type Repr = WRuntimeState;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Created => WRuntimeState::Created,
            Self::Initialized => WRuntimeState::Initialized,
            Self::Running => WRuntimeState::Running,
            Self::PauseRequested => WRuntimeState::PauseRequested,
            Self::Paused => WRuntimeState::Paused,
            Self::Restarting => WRuntimeState::Restarting,
            Self::Converged => WRuntimeState::Converged,
            Self::Exhausted => WRuntimeState::Exhausted,
            Self::Failed => WRuntimeState::Failed,
            Self::Cancelled => WRuntimeState::Cancelled,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WRuntimeState::Created => Self::Created,
            WRuntimeState::Initialized => Self::Initialized,
            WRuntimeState::Running => Self::Running,
            WRuntimeState::PauseRequested => Self::PauseRequested,
            WRuntimeState::Paused => Self::Paused,
            WRuntimeState::Restarting => Self::Restarting,
            WRuntimeState::Converged => Self::Converged,
            WRuntimeState::Exhausted => Self::Exhausted,
            WRuntimeState::Failed => Self::Failed,
            WRuntimeState::Cancelled => Self::Cancelled,
        })
    }
}
wire_struct!(RuntimeConfig, WRuntimeConfig {
    tolerance_floors: Option<ToleranceFloors>,
    autonomy: Option<AutonomyPolicy>,
    max_audit_records: usize,
    max_pending: usize,
});
wire_struct!(
    ToleranceFloors,
    WToleranceFloors {
        density_residual: f64,
        energy_hartree: f64,
        eigen_hartree: f64,
    }
);
wire_struct!(AutonomyPolicy, WAutonomyPolicy {
    max_applied_total: usize,
    rules: Vec<AutonomyRule>,
});
wire_struct!(
    AutonomyRule,
    WAutonomyRule {
        action: AllowedAction,
        max_count: usize,
        cooldown_iterations: usize,
        max_version_lag: u64,
        allow_revalidation: bool,
    }
);
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WAllowedAction {
    Beta {
        minimum: <f64 as Wire>::Repr,
        maximum: <f64 as Wire>::Repr,
        only_decrease: <bool as Wire>::Repr,
    },
    ResetHistory,
    Eigensolver(<RuntimeEigensolver as Wire>::Repr),
    TightenTolerances,
    ExtendBudget {
        maximum: <usize as Wire>::Repr,
    },
    Pause,
    Resume,
    Cancel,
}
impl Wire for AllowedAction {
    type Repr = WAllowedAction;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Beta {
                minimum,
                maximum,
                only_decrease,
            } => WAllowedAction::Beta {
                minimum: minimum.encode(_bytes)?,
                maximum: maximum.encode(_bytes)?,
                only_decrease: only_decrease.encode(_bytes)?,
            },
            Self::ResetHistory => WAllowedAction::ResetHistory,
            Self::Eigensolver(v) => WAllowedAction::Eigensolver(v.encode(_bytes)?),
            Self::TightenTolerances => WAllowedAction::TightenTolerances,
            Self::ExtendBudget { maximum } => WAllowedAction::ExtendBudget {
                maximum: maximum.encode(_bytes)?,
            },
            Self::Pause => WAllowedAction::Pause,
            Self::Resume => WAllowedAction::Resume,
            Self::Cancel => WAllowedAction::Cancel,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WAllowedAction::Beta {
                minimum,
                maximum,
                only_decrease,
            } => Self::Beta {
                minimum: Wire::decode(minimum, _bytes, _cursor)?,
                maximum: Wire::decode(maximum, _bytes, _cursor)?,
                only_decrease: Wire::decode(only_decrease, _bytes, _cursor)?,
            },
            WAllowedAction::ResetHistory => Self::ResetHistory,
            WAllowedAction::Eigensolver(v) => Self::Eigensolver(Wire::decode(v, _bytes, _cursor)?),
            WAllowedAction::TightenTolerances => Self::TightenTolerances,
            WAllowedAction::ExtendBudget { maximum } => Self::ExtendBudget {
                maximum: Wire::decode(maximum, _bytes, _cursor)?,
            },
            WAllowedAction::Pause => Self::Pause,
            WAllowedAction::Resume => Self::Resume,
            WAllowedAction::Cancel => Self::Cancel,
        })
    }
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WRuntimeEigensolver {
    Davidson,
    Hybrid,
}
impl Wire for RuntimeEigensolver {
    type Repr = WRuntimeEigensolver;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Davidson => WRuntimeEigensolver::Davidson,
            Self::Hybrid => WRuntimeEigensolver::Hybrid,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WRuntimeEigensolver::Davidson => Self::Davidson,
            WRuntimeEigensolver::Hybrid => Self::Hybrid,
        })
    }
}
wire_struct!(AuditRecord, WAuditRecord {
    schema_version: u32,
    calc_id: String,
    sequence: u64,
    timestamp: SystemTime,
    request: ControlRequest,
    requester: String,
    approval: Option<Approval>,
    authority: Authority,
    replayed: bool,
    autonomy_rule_index: Option<usize>,
    before: BoundarySnapshot,
    after: BoundarySnapshot,
    outcome: AuditOutcome,
    history_reset: bool,
});
wire_struct!(ControlRequest, WControlRequest {
    schema_version: u32,
    calc_id: String,
    id: String,
    observed_state_version: u64,
    observed_policy_revision: u64,
    action: ControlAction,
    preconditions: Vec<Precondition>,
    expiry: Expiry,
    stale_policy: StalePolicy,
    reason: String,
});
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WControlAction {
    Policy(<Vec<PolicyChange> as Wire>::Repr),
    Pause,
    Resume,
    Cancel,
    SwitchMode(<ExecutionMode as Wire>::Repr),
    Checkpoint,
}
impl Wire for ControlAction {
    type Repr = WControlAction;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Policy(v) => WControlAction::Policy(v.encode(_bytes)?),
            Self::Pause => WControlAction::Pause,
            Self::Resume => WControlAction::Resume,
            Self::Cancel => WControlAction::Cancel,
            Self::SwitchMode(v) => WControlAction::SwitchMode(v.encode(_bytes)?),
            Self::Checkpoint => WControlAction::Checkpoint,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WControlAction::Policy(v) => Self::Policy(Wire::decode(v, _bytes, _cursor)?),
            WControlAction::Pause => Self::Pause,
            WControlAction::Resume => Self::Resume,
            WControlAction::Cancel => Self::Cancel,
            WControlAction::SwitchMode(v) => Self::SwitchMode(Wire::decode(v, _bytes, _cursor)?),
            WControlAction::Checkpoint => Self::Checkpoint,
        })
    }
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WPolicyChange {
    Mixing(<Mixing as Wire>::Repr),
    MixingBeta(<f64 as Wire>::Repr),
    MixingHistory(<usize as Wire>::Repr),
    ResetHistory,
    Eigensolver(<RuntimeEigensolver as Wire>::Repr),
    DensityTolerance(<f64 as Wire>::Repr),
    EnergyToleranceHartree(<f64 as Wire>::Repr),
    EigenToleranceHartree(<f64 as Wire>::Repr),
    MaxIterations(<usize as Wire>::Repr),
    ScheduleFlags {
        residual_eigen_schedule: <bool as Wire>::Repr,
    },
    DivergenceParameters {
        factor: <f64 as Wire>::Repr,
        patience: <usize as Wire>::Repr,
    },
}
impl Wire for PolicyChange {
    type Repr = WPolicyChange;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Mixing(v) => WPolicyChange::Mixing(v.encode(_bytes)?),
            Self::MixingBeta(v) => WPolicyChange::MixingBeta(v.encode(_bytes)?),
            Self::MixingHistory(v) => WPolicyChange::MixingHistory(v.encode(_bytes)?),
            Self::ResetHistory => WPolicyChange::ResetHistory,
            Self::Eigensolver(v) => WPolicyChange::Eigensolver(v.encode(_bytes)?),
            Self::DensityTolerance(v) => WPolicyChange::DensityTolerance(v.encode(_bytes)?),
            Self::EnergyToleranceHartree(v) => {
                WPolicyChange::EnergyToleranceHartree(v.encode(_bytes)?)
            }
            Self::EigenToleranceHartree(v) => {
                WPolicyChange::EigenToleranceHartree(v.encode(_bytes)?)
            }
            Self::MaxIterations(v) => WPolicyChange::MaxIterations(v.encode(_bytes)?),
            Self::ScheduleFlags {
                residual_eigen_schedule,
            } => WPolicyChange::ScheduleFlags {
                residual_eigen_schedule: residual_eigen_schedule.encode(_bytes)?,
            },
            Self::DivergenceParameters { factor, patience } => {
                WPolicyChange::DivergenceParameters {
                    factor: factor.encode(_bytes)?,
                    patience: patience.encode(_bytes)?,
                }
            }
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WPolicyChange::Mixing(v) => Self::Mixing(Wire::decode(v, _bytes, _cursor)?),
            WPolicyChange::MixingBeta(v) => Self::MixingBeta(Wire::decode(v, _bytes, _cursor)?),
            WPolicyChange::MixingHistory(v) => {
                Self::MixingHistory(Wire::decode(v, _bytes, _cursor)?)
            }
            WPolicyChange::ResetHistory => Self::ResetHistory,
            WPolicyChange::Eigensolver(v) => Self::Eigensolver(Wire::decode(v, _bytes, _cursor)?),
            WPolicyChange::DensityTolerance(v) => {
                Self::DensityTolerance(Wire::decode(v, _bytes, _cursor)?)
            }
            WPolicyChange::EnergyToleranceHartree(v) => {
                Self::EnergyToleranceHartree(Wire::decode(v, _bytes, _cursor)?)
            }
            WPolicyChange::EigenToleranceHartree(v) => {
                Self::EigenToleranceHartree(Wire::decode(v, _bytes, _cursor)?)
            }
            WPolicyChange::MaxIterations(v) => {
                Self::MaxIterations(Wire::decode(v, _bytes, _cursor)?)
            }
            WPolicyChange::ScheduleFlags {
                residual_eigen_schedule,
            } => Self::ScheduleFlags {
                residual_eigen_schedule: Wire::decode(residual_eigen_schedule, _bytes, _cursor)?,
            },
            WPolicyChange::DivergenceParameters { factor, patience } => {
                Self::DivergenceParameters {
                    factor: Wire::decode(factor, _bytes, _cursor)?,
                    patience: Wire::decode(patience, _bytes, _cursor)?,
                }
            }
        })
    }
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WPrecondition {
    State(<RuntimeState as Wire>::Repr),
    PolicyRevision(<u64 as Wire>::Repr),
    Mode(<ExecutionMode as Wire>::Repr),
    Scalar {
        field: <Scalar as Wire>::Repr,
        comparison: <Comparison as Wire>::Repr,
        value: <f64 as Wire>::Repr,
    },
    DiagnosticActive,
}
impl Wire for Precondition {
    type Repr = WPrecondition;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::State(v) => WPrecondition::State(v.encode(_bytes)?),
            Self::PolicyRevision(v) => WPrecondition::PolicyRevision(v.encode(_bytes)?),
            Self::Mode(v) => WPrecondition::Mode(v.encode(_bytes)?),
            Self::Scalar {
                field,
                comparison,
                value,
            } => WPrecondition::Scalar {
                field: field.encode(_bytes)?,
                comparison: comparison.encode(_bytes)?,
                value: value.encode(_bytes)?,
            },
            Self::DiagnosticActive => WPrecondition::DiagnosticActive,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WPrecondition::State(v) => Self::State(Wire::decode(v, _bytes, _cursor)?),
            WPrecondition::PolicyRevision(v) => {
                Self::PolicyRevision(Wire::decode(v, _bytes, _cursor)?)
            }
            WPrecondition::Mode(v) => Self::Mode(Wire::decode(v, _bytes, _cursor)?),
            WPrecondition::Scalar {
                field,
                comparison,
                value,
            } => Self::Scalar {
                field: Wire::decode(field, _bytes, _cursor)?,
                comparison: Wire::decode(comparison, _bytes, _cursor)?,
                value: Wire::decode(value, _bytes, _cursor)?,
            },
            WPrecondition::DiagnosticActive => Self::DiagnosticActive,
        })
    }
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WScalar {
    DensityResidual,
    ChargeErrorRelative,
    EnergyChangeAbsHartree,
}
impl Wire for Scalar {
    type Repr = WScalar;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::DensityResidual => WScalar::DensityResidual,
            Self::ChargeErrorRelative => WScalar::ChargeErrorRelative,
            Self::EnergyChangeAbsHartree => WScalar::EnergyChangeAbsHartree,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WScalar::DensityResidual => Self::DensityResidual,
            WScalar::ChargeErrorRelative => Self::ChargeErrorRelative,
            WScalar::EnergyChangeAbsHartree => Self::EnergyChangeAbsHartree,
        })
    }
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WComparison {
    Less,
    LessOrEqual,
    Greater,
    GreaterOrEqual,
}
impl Wire for Comparison {
    type Repr = WComparison;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Less => WComparison::Less,
            Self::LessOrEqual => WComparison::LessOrEqual,
            Self::Greater => WComparison::Greater,
            Self::GreaterOrEqual => WComparison::GreaterOrEqual,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WComparison::Less => Self::Less,
            WComparison::LessOrEqual => Self::LessOrEqual,
            WComparison::Greater => Self::Greater,
            WComparison::GreaterOrEqual => Self::GreaterOrEqual,
        })
    }
}
wire_struct!(Expiry, WExpiry {
    max_iteration: Option<usize>,
    deadline: Option<SystemTime>,
});
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WStalePolicy {
    Reject,
    Revalidate,
}
impl Wire for StalePolicy {
    type Repr = WStalePolicy;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Reject => WStalePolicy::Reject,
            Self::Revalidate => WStalePolicy::Revalidate,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WStalePolicy::Reject => Self::Reject,
            WStalePolicy::Revalidate => Self::Revalidate,
        })
    }
}
wire_struct!(
    Approval,
    WApproval {
        id: String,
        approver: String,
        expires_at: SystemTime,
        allow_revalidation: bool,
    }
);
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WAuthority {
    PendingApproval,
    Approved(<Approval as Wire>::Repr),
    Autonomous,
}
impl Wire for Authority {
    type Repr = WAuthority;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::PendingApproval => WAuthority::PendingApproval,
            Self::Approved(v) => WAuthority::Approved(v.encode(_bytes)?),
            Self::Autonomous => WAuthority::Autonomous,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WAuthority::PendingApproval => Self::PendingApproval,
            WAuthority::Approved(v) => Self::Approved(Wire::decode(v, _bytes, _cursor)?),
            WAuthority::Autonomous => Self::Autonomous,
        })
    }
}
wire_struct!(
    BoundarySnapshot,
    WBoundarySnapshot {
        state: RuntimeState,
        mode: ExecutionMode,
        iteration: usize,
        state_version: u64,
        policy_revision: u64,
        options: ScfOptions,
    }
);
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WAuditOutcome {
    Accepted,
    Applied,
    Rejected(<ControlError as Wire>::Repr),
}
impl Wire for AuditOutcome {
    type Repr = WAuditOutcome;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Accepted => WAuditOutcome::Accepted,
            Self::Applied => WAuditOutcome::Applied,
            Self::Rejected(v) => WAuditOutcome::Rejected(v.encode(_bytes)?),
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WAuditOutcome::Accepted => Self::Accepted,
            WAuditOutcome::Applied => Self::Applied,
            WAuditOutcome::Rejected(v) => Self::Rejected(Wire::decode(v, _bytes, _cursor)?),
        })
    }
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WControlError {
    InvalidSchema,
    WrongCalculation,
    StaleState,
    FutureState,
    StalePolicy,
    Expired,
    PreconditionFailed,
    ApprovalRequired,
    ApprovalExpired,
    ApprovalDenied,
    AutonomousDenied,
    Cooldown,
    Quota,
    Capacity,
    Conflict,
    ModeDisallowsControl,
    ModeUnavailableAfterStart,
    Terminal,
    NoChange,
    UnsupportedCapability,
    PrepareFailed,
    JournalFailed,
    ReplayMismatch,
    Superseded,
    Policy(<PolicyError as Wire>::Repr),
    InvalidTransition(<InvalidTransition as Wire>::Repr),
}
impl Wire for ControlError {
    type Repr = WControlError;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::InvalidSchema => WControlError::InvalidSchema,
            Self::WrongCalculation => WControlError::WrongCalculation,
            Self::StaleState => WControlError::StaleState,
            Self::FutureState => WControlError::FutureState,
            Self::StalePolicy => WControlError::StalePolicy,
            Self::Expired => WControlError::Expired,
            Self::PreconditionFailed => WControlError::PreconditionFailed,
            Self::ApprovalRequired => WControlError::ApprovalRequired,
            Self::ApprovalExpired => WControlError::ApprovalExpired,
            Self::ApprovalDenied => WControlError::ApprovalDenied,
            Self::AutonomousDenied => WControlError::AutonomousDenied,
            Self::Cooldown => WControlError::Cooldown,
            Self::Quota => WControlError::Quota,
            Self::Capacity => WControlError::Capacity,
            Self::Conflict => WControlError::Conflict,
            Self::ModeDisallowsControl => WControlError::ModeDisallowsControl,
            Self::ModeUnavailableAfterStart => WControlError::ModeUnavailableAfterStart,
            Self::Terminal => WControlError::Terminal,
            Self::NoChange => WControlError::NoChange,
            Self::UnsupportedCapability => WControlError::UnsupportedCapability,
            Self::PrepareFailed => WControlError::PrepareFailed,
            Self::JournalFailed => WControlError::JournalFailed,
            Self::ReplayMismatch => WControlError::ReplayMismatch,
            Self::Superseded => WControlError::Superseded,
            Self::Policy(v) => WControlError::Policy(v.encode(_bytes)?),
            Self::InvalidTransition(v) => WControlError::InvalidTransition(v.encode(_bytes)?),
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WControlError::InvalidSchema => Self::InvalidSchema,
            WControlError::WrongCalculation => Self::WrongCalculation,
            WControlError::StaleState => Self::StaleState,
            WControlError::FutureState => Self::FutureState,
            WControlError::StalePolicy => Self::StalePolicy,
            WControlError::Expired => Self::Expired,
            WControlError::PreconditionFailed => Self::PreconditionFailed,
            WControlError::ApprovalRequired => Self::ApprovalRequired,
            WControlError::ApprovalExpired => Self::ApprovalExpired,
            WControlError::ApprovalDenied => Self::ApprovalDenied,
            WControlError::AutonomousDenied => Self::AutonomousDenied,
            WControlError::Cooldown => Self::Cooldown,
            WControlError::Quota => Self::Quota,
            WControlError::Capacity => Self::Capacity,
            WControlError::Conflict => Self::Conflict,
            WControlError::ModeDisallowsControl => Self::ModeDisallowsControl,
            WControlError::ModeUnavailableAfterStart => Self::ModeUnavailableAfterStart,
            WControlError::Terminal => Self::Terminal,
            WControlError::NoChange => Self::NoChange,
            WControlError::UnsupportedCapability => Self::UnsupportedCapability,
            WControlError::PrepareFailed => Self::PrepareFailed,
            WControlError::JournalFailed => Self::JournalFailed,
            WControlError::ReplayMismatch => Self::ReplayMismatch,
            WControlError::Superseded => Self::Superseded,
            WControlError::Policy(v) => Self::Policy(Wire::decode(v, _bytes, _cursor)?),
            WControlError::InvalidTransition(v) => {
                Self::InvalidTransition(Wire::decode(v, _bytes, _cursor)?)
            }
        })
    }
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WPolicyError {
    Empty,
    NoChange,
    TooManyChanges,
    DuplicateField,
    NonFiniteOrNonPositive,
    InvalidBeta,
    HistoryQuota,
    HistoryRequiresPulay,
    IterationQuota,
    BudgetMustIncrease,
    MissingToleranceFloors,
    InvalidToleranceFloors,
    BelowToleranceFloor,
    ToleranceRelaxation,
    CertificateThreshold,
    RestrictedSchedule,
    RestrictedDivergence,
    LockingRequiresDavidson,
}
impl Wire for PolicyError {
    type Repr = WPolicyError;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Empty => WPolicyError::Empty,
            Self::NoChange => WPolicyError::NoChange,
            Self::TooManyChanges => WPolicyError::TooManyChanges,
            Self::DuplicateField => WPolicyError::DuplicateField,
            Self::NonFiniteOrNonPositive => WPolicyError::NonFiniteOrNonPositive,
            Self::InvalidBeta => WPolicyError::InvalidBeta,
            Self::HistoryQuota => WPolicyError::HistoryQuota,
            Self::HistoryRequiresPulay => WPolicyError::HistoryRequiresPulay,
            Self::IterationQuota => WPolicyError::IterationQuota,
            Self::BudgetMustIncrease => WPolicyError::BudgetMustIncrease,
            Self::MissingToleranceFloors => WPolicyError::MissingToleranceFloors,
            Self::InvalidToleranceFloors => WPolicyError::InvalidToleranceFloors,
            Self::BelowToleranceFloor => WPolicyError::BelowToleranceFloor,
            Self::ToleranceRelaxation => WPolicyError::ToleranceRelaxation,
            Self::CertificateThreshold => WPolicyError::CertificateThreshold,
            Self::RestrictedSchedule => WPolicyError::RestrictedSchedule,
            Self::RestrictedDivergence => WPolicyError::RestrictedDivergence,
            Self::LockingRequiresDavidson => WPolicyError::LockingRequiresDavidson,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WPolicyError::Empty => Self::Empty,
            WPolicyError::NoChange => Self::NoChange,
            WPolicyError::TooManyChanges => Self::TooManyChanges,
            WPolicyError::DuplicateField => Self::DuplicateField,
            WPolicyError::NonFiniteOrNonPositive => Self::NonFiniteOrNonPositive,
            WPolicyError::InvalidBeta => Self::InvalidBeta,
            WPolicyError::HistoryQuota => Self::HistoryQuota,
            WPolicyError::HistoryRequiresPulay => Self::HistoryRequiresPulay,
            WPolicyError::IterationQuota => Self::IterationQuota,
            WPolicyError::BudgetMustIncrease => Self::BudgetMustIncrease,
            WPolicyError::MissingToleranceFloors => Self::MissingToleranceFloors,
            WPolicyError::InvalidToleranceFloors => Self::InvalidToleranceFloors,
            WPolicyError::BelowToleranceFloor => Self::BelowToleranceFloor,
            WPolicyError::ToleranceRelaxation => Self::ToleranceRelaxation,
            WPolicyError::CertificateThreshold => Self::CertificateThreshold,
            WPolicyError::RestrictedSchedule => Self::RestrictedSchedule,
            WPolicyError::RestrictedDivergence => Self::RestrictedDivergence,
            WPolicyError::LockingRequiresDavidson => Self::LockingRequiresDavidson,
        })
    }
}
wire_struct!(
    InvalidTransition,
    WInvalidTransition {
        from: RuntimeState,
        trigger: Trigger,
    }
);
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
enum WTrigger {
    Initialize,
    Start,
    Complete,
    CertificatePassed,
    BudgetExhausted,
    Fail,
    Pause,
    Resume,
    ResumeComplete,
    Cancel,
}
impl Wire for Trigger {
    type Repr = WTrigger;
    fn encode(&self, _bytes: &mut Vec<u8>) -> Result<Self::Repr, Error> {
        Ok(match self {
            Self::Initialize => WTrigger::Initialize,
            Self::Start => WTrigger::Start,
            Self::Complete => WTrigger::Complete,
            Self::CertificatePassed => WTrigger::CertificatePassed,
            Self::BudgetExhausted => WTrigger::BudgetExhausted,
            Self::Fail => WTrigger::Fail,
            Self::Pause => WTrigger::Pause,
            Self::Resume => WTrigger::Resume,
            Self::ResumeComplete => WTrigger::ResumeComplete,
            Self::Cancel => WTrigger::Cancel,
        })
    }
    fn decode(w: Self::Repr, _bytes: &[u8], _cursor: &mut usize) -> Result<Self, Error> {
        Ok(match w {
            WTrigger::Initialize => Self::Initialize,
            WTrigger::Start => Self::Start,
            WTrigger::Complete => Self::Complete,
            WTrigger::CertificatePassed => Self::CertificatePassed,
            WTrigger::BudgetExhausted => Self::BudgetExhausted,
            WTrigger::Fail => Self::Fail,
            WTrigger::Pause => Self::Pause,
            WTrigger::Resume => Self::Resume,
            WTrigger::ResumeComplete => Self::ResumeComplete,
            WTrigger::Cancel => Self::Cancel,
        })
    }
}
wire_struct!(ExportResult, WExportResult {
    converged: bool,
    energy: EnergyAudit,
    charge: f64,
    gap_hartree: Option<f64>,
    eigenvalues: Vec<Vec<f64>>,
    residuals: Vec<Vec<f64>>,
    history: Vec<IterationRecord>,
});

#[derive(Clone, Copy, Debug, clap::ValueEnum)]
pub enum Kind {
    Exact,
    Warm,
    Result,
}
fn context(
    input: &Path,
) -> Result<(dft_engine::ScfProblem, ModelIdentity, ExecutionIdentity), Error> {
    use crate::{
        config::Config,
        run::{options, sha256},
    };
    use dft_engine::{OccupationPolicy, ScfProblem, SpeciesPotential};
    let config = Config::parse(&fs::read_to_string(input)?)?;
    if config.finite_temperature()?.is_some() {
        return Err(invalid(
            "checkpoint.occupation_policy: v1 fixed-occupation checkpoint cannot be resumed under fermi_dirac; start a new run",
        ));
    }
    if config.xc.density_floor.is_some() {
        return Err(invalid(
            "xc.density_floor: checkpoint v1 cannot encode the XC floor policy; use run without checkpoint/restart",
        ));
    }
    let cell = config.cell()?;
    let opts = options(&config)?;
    let atoms: Vec<_> = config
        .system
        .atoms
        .iter()
        .map(|a| dft_geometry::Atom::new(&a.species, a.fractional))
        .collect();
    let mut potentials = Vec::new();
    let mut hashes = std::collections::BTreeMap::new();
    for (label, path) in &config.pseudopotentials {
        let bytes = fs::read(input.parent().unwrap_or(Path::new(".")).join(path))?;
        hashes.insert(label.clone(), sha256(&bytes));
        potentials.push((
            label.clone(),
            dft_upf::parse_upf(std::str::from_utf8(&bytes)?)?,
        ));
    }
    let refs: Vec<_> = potentials
        .iter()
        .map(|(label, potential)| SpeciesPotential { label, potential })
        .collect();
    let points = config.basis.resolve_kpoints()?;
    let problem = if config.basis.legacy_gamma() {
        ScfProblem::new(&cell, &atoms, &refs, opts.clone())?
    } else {
        ScfProblem::with_kpoints(
            &cell,
            &atoms,
            &refs,
            opts.clone(),
            &points,
            OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-7,
            },
        )?
    };
    let model = ModelIdentity {
        cell_bohr: cell.lattice_columns_bohr(),
        atoms_fractional: atoms.iter().map(|a| a.fractional_position).collect(),
        species: atoms.iter().map(|a| a.species.clone()).collect(),
        pseudopotential_sha256: atoms
            .iter()
            .map(|a| {
                hashes
                    .get(&a.species)
                    .cloned()
                    .ok_or_else(|| invalid("species hash missing"))
            })
            .collect::<Result<_, _>>()?,
        cutoff_hartree: opts.cutoff_hartree,
        kpoints: points
            .iter()
            .map(|k| [k.fractional[0], k.fractional[1], k.fractional[2], k.weight])
            .collect(),
        electrons: problem.electron_count(),
        xc: opts.functional.name().into(),
        spin: "unpolarized".into(),
        environment: "none".into(),
    };
    // Hash the actual executable, including uncommitted code and linked backend.
    let executable = fs::read(std::env::current_exe()?)?;
    let host = std::process::Command::new("uname").arg("-a").output()?;
    if !host.status.success() {
        return Err(invalid("execution platform unavailable"));
    }
    let execution = ExecutionIdentity {
        build_sha256: sha256(&executable),
        platform: format!(
            "{}:{}:{}",
            std::env::consts::ARCH,
            std::env::consts::OS,
            String::from_utf8(host.stdout)?.trim()
        ),
        threads: opts.threads,
    };
    Ok((problem, model, execution))
}
fn finish(session: &mut dft_engine::ScfSession) -> Result<(), Error> {
    while session.snapshot().status == dft_engine::StepStatus::Ready {
        session.step()?;
    }
    if session.result().is_none() {
        return Err(invalid("no evaluated result"));
    }
    Ok(())
}
pub fn create(input: &Path, output: &Path, iteration: usize, kind: Kind) -> Result<(), Error> {
    let (problem, model, execution) = context(input)?;
    let mut session = dft_engine::ScfSession::new(
        problem,
        CompletionMode::RetainLast,
        ExecutionConfig::default(),
    )?;
    session.initialize()?;
    let kind = match kind {
        Kind::Exact => CheckpointKind::ExactContinuation,
        Kind::Warm => CheckpointKind::WarmStart,
        Kind::Result => CheckpointKind::ResultExport,
    };
    if kind == CheckpointKind::ResultExport {
        finish(&mut session)?;
    } else {
        for _ in 0..iteration {
            session.step()?;
        }
    }
    let cp = session.checkpoint(kind, model, execution)?;
    write(output, &cp)?;
    if let Some(p) = cp.standard {
        println!(
            "convergence_profile: standard; certified: false; standard_threshold_ha: {}",
            p.threshold_ha
        );
    }
    println!(
        "Wrote {:?}, iteration {}: {}",
        kind,
        session.snapshot().iteration,
        output.display()
    );
    Ok(())
}
pub fn restart(input: &Path, source: &Path, output: &Path, warm: bool) -> Result<(), Error> {
    let cp = read(source)?;
    if cp.runtime.is_some() {
        return Err(invalid(
            "runtime checkpoint requires an authenticated runtime host",
        ));
    }
    let (problem, model, execution) = context(input)?;
    let mut session = if warm {
        dft_engine::ScfSession::warm_start(
            problem,
            &cp,
            &model,
            &execution,
            CompletionMode::RetainLast,
            ExecutionConfig {
                calc_id: format!("warm-{}", std::process::id()),
                ..Default::default()
            },
        )?
    } else {
        dft_engine::ScfSession::restore_exact(problem, &cp, &model, &execution)?
    };
    finish(&mut session)?;
    let result = session.checkpoint(CheckpointKind::ResultExport, model, execution)?;
    write(output, &result)?;
    if let Some(p) = result.standard {
        println!(
            "convergence_profile: standard; certified: false; standard_threshold_ha: {}",
            p.threshold_ha
        );
    }
    let r = result.result.as_ref().unwrap();
    println!(
        "Energy {:.17} Ha, converged={}, iterations={}",
        r.energy.total_hartree(),
        r.converged,
        r.history.len()
    );
    if !r.converged {
        return Err(invalid(
            "NOT CONVERGED: result exported without convergence",
        ));
    }
    Ok(())
}

#[cfg(test)]
mod pw92_wire_tests {
    use super::*;
    #[test]
    fn capacity_wire_roundtrip_and_legacy_omission() {
        for multiplier in [2, 3, 4, 6, 12] {
            let options = ScfOptions {
                davidson_capacity_per_band: multiplier,
                ..Default::default()
            };
            let mut bytes = Vec::new();
            let json = serde_json::to_value(options.encode(&mut bytes).unwrap()).unwrap();
            assert_eq!(
                json.get("davidson_capacity_per_band").is_none(),
                multiplier == 6
            );
            let decoded =
                ScfOptions::decode(serde_json::from_value(json).unwrap(), &bytes, &mut 0).unwrap();
            assert_eq!(decoded.davidson_capacity_per_band, multiplier);
        }
    }
    #[test]
    fn standard_stall_guard_wire_roundtrip_and_legacy_omission() {
        for enabled in [false, true] {
            let options = ScfOptions {
                standard: Some(StandardProfile::default()),
                standard_stall_guard: enabled,
                ..Default::default()
            };
            let mut bytes = Vec::new();
            let wire = options.encode(&mut bytes).unwrap();
            let json = serde_json::to_value(wire).unwrap();
            assert_eq!(json.get("standard_stall_guard").is_some(), enabled);
            let decoded =
                ScfOptions::decode(serde_json::from_value(json).unwrap(), &bytes, &mut 0).unwrap();
            assert_eq!(decoded.standard_stall_guard, enabled);
        }
    }
    #[test]
    fn standard_density_forcing_wire_roundtrip_and_legacy_omission() {
        for enabled in [false, true] {
            let options = ScfOptions {
                standard: Some(StandardProfile::default()),
                standard_density_forcing: enabled,
                ..Default::default()
            };
            let mut bytes = Vec::new();
            let wire = options.encode(&mut bytes).unwrap();
            let json = serde_json::to_value(wire).unwrap();
            assert_eq!(json.get("standard_density_forcing").is_some(), enabled);
            let decoded =
                ScfOptions::decode(serde_json::from_value(json).unwrap(), &bytes, &mut 0).unwrap();
            assert_eq!(decoded.standard_density_forcing, enabled);
        }
    }
    #[test]
    fn convention_roundtrip_legacy_default_and_rejection() {
        for prefactor in [
            dft_engine::Pw92Prefactor::Modified,
            dft_engine::Pw92Prefactor::Original,
        ] {
            let options = ScfOptions {
                pw92_prefactor: prefactor,
                ..Default::default()
            };
            let mut bytes = Vec::new();
            let wire = options.encode(&mut bytes).unwrap();
            let mut json = serde_json::to_value(wire).unwrap();
            if prefactor == dft_engine::Pw92Prefactor::Modified {
                assert!(json.get("pw92_prefactor").is_none());
            } else {
                assert_eq!(json["pw92_prefactor"], "original");
            }
            let decoded = ScfOptions::decode(
                serde_json::from_value(json.clone()).unwrap(),
                &bytes,
                &mut 0,
            )
            .unwrap();
            assert_eq!(decoded.pw92_prefactor, prefactor);
            json["pw92_prefactor"] = serde_json::json!("unknown");
            assert!(
                ScfOptions::decode(serde_json::from_value(json).unwrap(), &bytes, &mut 0).is_err()
            );
        }
    }
}

#[cfg(test)]
mod pseudo_wire_tests {
    use super::*;
    #[test]
    fn defaults_legacy_and_research_roundtrip() {
        for local_mesh in [dft_engine::LocalMesh::Baandi, dft_engine::LocalMesh::QeLike] {
            for radial_even_rule in [
                dft_engine::RadialEvenRule::Simpson3_8,
                dft_engine::RadialEvenRule::SimpsonQuadraticEnd,
            ] {
                let options = ScfOptions {
                    local_mesh,
                    radial_even_rule,
                    ..Default::default()
                };
                let mut bytes = Vec::new();
                let mut json = serde_json::to_value(options.encode(&mut bytes).unwrap()).unwrap();
                assert_eq!(
                    json.get("local_mesh").is_none(),
                    local_mesh == dft_engine::LocalMesh::Baandi
                );
                assert_eq!(
                    json.get("radial_even_rule").is_none(),
                    radial_even_rule == dft_engine::RadialEvenRule::Simpson3_8
                );
                let decoded = ScfOptions::decode(
                    serde_json::from_value(json.clone()).unwrap(),
                    &bytes,
                    &mut 0,
                )
                .unwrap();
                assert_eq!(decoded.local_mesh, local_mesh);
                assert_eq!(decoded.radial_even_rule, radial_even_rule);
                for key in ["local_mesh", "radial_even_rule"] {
                    let old = json.clone();
                    json[key] = serde_json::json!("unknown");
                    assert!(
                        ScfOptions::decode(
                            serde_json::from_value(json.clone()).unwrap(),
                            &bytes,
                            &mut 0
                        )
                        .is_err()
                    );
                    json = old;
                }
            }
        }
    }
}
