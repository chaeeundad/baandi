//! Reference checkpoint v2. Pure bounded state encoding; persistence and SHA-256
//! belong to the host. A checkpoint contains one rank's owned rows and planes.
//! No compiler-generated continuation bytes, pointers, or scratch plans are saved.
use crate::model::ScfError;
use dft_numerics::{RealMatrix, linalg::ComplexMatrix};
use dft_physics::data::Complex64;

pub const SCHEMA_VERSION: u32 = 2;
pub const MAX_BYTES: usize = 512 * 1024 * 1024;
pub(crate) fn invalid(message: &str) -> ScfError {
    ScfError::Numerical(format!("reference checkpoint v2: {message}"))
}
/// Trusted host identity, computed from actual input/UPF/executable bytes.
/// Rank count, thread count, backend and effective switches must be identical.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Identity {
    pub input_sha256: String,
    pub upf_sha256: Vec<String>,
    pub binary_sha256: String,
    pub platform: String,
    pub threads: usize,
    pub ranks: usize,
}
impl Identity {
    pub fn validate(&self) -> Result<(), ScfError> {
        let digest = |s: &str| s.len() == 64 && s.bytes().all(|b| b.is_ascii_hexdigit());
        if !digest(&self.input_sha256)
            || !digest(&self.binary_sha256)
            || self.upf_sha256.is_empty()
            || !self.upf_sha256.iter().all(|s| digest(s))
            || self.platform.is_empty()
            || self.threads == 0
            || self.ranks == 0
        {
            return Err(invalid("invalid input/UPF/binary/execution identity"));
        }
        Ok(())
    }
}
#[derive(Clone, Debug)]
pub struct Checkpoint {
    pub schema_version: u32,
    pub identity: Identity,
    pub rank: usize,
    pub iteration: usize,
    /// Opaque versioned numerical payload, not opaque Rust memory.
    pub state: Vec<u8>,
}
impl Checkpoint {
    pub fn validate(&self, expected: &Identity, rank: usize) -> Result<(), ScfError> {
        self.identity.validate()?;
        expected.validate()?;
        if self.schema_version != SCHEMA_VERSION {
            return Err(invalid("unsupported version"));
        }
        if &self.identity != expected {
            return Err(invalid("input/UPF/binary/execution identity mismatch"));
        }
        if self.rank != rank || rank >= expected.ranks {
            return Err(invalid("rank/layout mismatch"));
        }
        if self.state.is_empty() || self.state.len() > MAX_BYTES {
            return Err(invalid("state size quota"));
        }
        Ok(())
    }
}

pub(super) struct Reader<'a> {
    data: &'a [u8],
    cursor: usize,
}
impl<'a> Reader<'a> {
    /// A reader over an in-memory encoding (GK9 k-pool record exchange).
    pub(super) fn new(data: &'a [u8]) -> Self {
        Self { data, cursor: 0 }
    }
    pub(super) fn remaining(&self) -> usize {
        self.data.len() - self.cursor
    }
    fn take(&mut self, n: usize) -> Result<&'a [u8], ScfError> {
        let end = self
            .cursor
            .checked_add(n)
            .filter(|&v| v <= self.data.len())
            .ok_or_else(|| invalid("truncated state or invalid length"))?;
        let out = &self.data[self.cursor..end];
        self.cursor = end;
        Ok(out)
    }
}
pub(super) trait Codec: Sized {
    fn put(&self, out: &mut Vec<u8>) -> Result<(), ScfError>;
    fn get(input: &mut Reader<'_>) -> Result<Self, ScfError>;
}
macro_rules! record {
    ($t:ty { $($f:ident),* $(,)? }) => {
        impl $crate::refscf::checkpoint::Codec for $t {
            fn put(&self, out: &mut Vec<u8>) -> Result<(), crate::model::ScfError> {
                $( $crate::refscf::checkpoint::Codec::put(&self.$f, out)?; )* Ok(())
            }
            fn get(input: &mut $crate::refscf::checkpoint::Reader<'_>) -> Result<Self, crate::model::ScfError> {
                Ok(Self { $( $f: $crate::refscf::checkpoint::Codec::get(input)?, )* })
            }
        }
    };
}
pub(super) use record;
impl Codec for u64 {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        if o.len() > MAX_BYTES - 8 {
            return Err(invalid("state quota"));
        }
        o.extend_from_slice(&self.to_le_bytes());
        Ok(())
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        Ok(Self::from_le_bytes(r.take(8)?.try_into().unwrap()))
    }
}
impl Codec for usize {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        (*self as u64).put(o)
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        Self::try_from(u64::get(r)?).map_err(|_| invalid("integer overflow"))
    }
}
impl Codec for f64 {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        self.to_bits().put(o)
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        Ok(Self::from_bits(u64::get(r)?))
    }
}
impl Codec for bool {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        u64::from(*self).put(o)
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        match u64::get(r)? {
            0 => Ok(false),
            1 => Ok(true),
            _ => Err(invalid("boolean tag")),
        }
    }
}
impl<T: Codec> Codec for Vec<T> {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        self.len().put(o)?;
        for v in self {
            v.put(o)?;
            if o.len() > MAX_BYTES {
                return Err(invalid("state quota"));
            }
        }
        Ok(())
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        let n = usize::get(r)?;
        if n > (r.data.len() - r.cursor) / 8 {
            return Err(invalid("array length quota"));
        }
        // Grow only after decoding each element; nested empty vectors cannot
        // turn a tiny hostile input into an unbounded allocation.
        let mut out = Vec::new();
        for _ in 0..n {
            out.push(T::get(r)?);
        }
        Ok(out)
    }
}
impl<T: Codec> Codec for Option<T> {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        self.is_some().put(o)?;
        if let Some(v) = self {
            v.put(o)?;
        }
        Ok(())
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        if bool::get(r)? {
            Ok(Some(T::get(r)?))
        } else {
            Ok(None)
        }
    }
}
impl<T: Codec, const N: usize> Codec for [T; N] {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        for v in self {
            v.put(o)?;
        }
        Ok(())
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        let mut v = Vec::new();
        for _ in 0..N {
            v.push(T::get(r)?);
        }
        v.try_into().map_err(|_| invalid("array length"))
    }
}
impl<A: Codec, B: Codec> Codec for (A, B) {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        self.0.put(o)?;
        self.1.put(o)
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        Ok((A::get(r)?, B::get(r)?))
    }
}
impl<A: Codec, B: Codec, C: Codec> Codec for (A, B, C) {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        self.0.put(o)?;
        self.1.put(o)?;
        self.2.put(o)
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        Ok((A::get(r)?, B::get(r)?, C::get(r)?))
    }
}
impl Codec for String {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        self.len().put(o)?;
        if self.len() > MAX_BYTES.saturating_sub(o.len()) {
            return Err(invalid("string quota"));
        }
        o.extend_from_slice(self.as_bytes());
        Ok(())
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        let n = usize::get(r)?;
        String::from_utf8(r.take(n)?.to_vec()).map_err(|_| invalid("UTF-8"))
    }
}
impl Codec for &'static str {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        self.to_string().put(o)
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        match String::get(r)?.as_str() {
            "fine" => Ok("fine"),
            "coarse" => Ok("coarse"),
            _ => Err(invalid("stage tag")),
        }
    }
}
impl Codec for std::time::Duration {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        self.as_secs().put(o)?;
        (self.subsec_nanos() as u64).put(o)
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        let s = u64::get(r)?;
        let n = u64::get(r)?;
        if n >= 1_000_000_000 {
            return Err(invalid("duration"));
        }
        Ok(Self::new(s, n as u32))
    }
}
record!(Complex64 { re, im });
record!(std::ops::Range<usize> { start, end });
macro_rules! matrix {
    ($t:ty, $element:ty) => {
        impl Codec for $t {
            fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
                self.rows().put(o)?;
                self.cols().put(o)?;
                for v in self.as_slice() {
                    v.put(o)?;
                }
                Ok(())
            }
            fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
                let rows = usize::get(r)?;
                let cols = usize::get(r)?;
                let n = rows
                    .checked_mul(cols)
                    .filter(|&n| n <= (r.data.len() - r.cursor) / std::mem::size_of::<$element>())
                    .ok_or_else(|| invalid("matrix shape quota"))?;
                let mut out = Self::zeros(rows, cols).map_err(super::numerical)?;
                for v in out.as_mut_slice().iter_mut().take(n) {
                    *v = <$element>::get(r)?;
                }
                Ok(out)
            }
        }
    };
}
matrix!(RealMatrix, f64);
matrix!(ComplexMatrix, Complex64);
impl Codec for crate::refstart::KOrbitals {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        self.0.put(o)
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        Ok(Self(Vec::get(r)?))
    }
}
/// An exchange is constructed only on explicit capture/restore. Drivers visit
/// it from a suspended boundary; normal map polling never constructs one.
pub(super) enum Exchange<'a> {
    Save(Vec<u8>),
    Load(Reader<'a>),
}
impl<'a> Exchange<'a> {
    pub fn save() -> Self {
        Self::Save(Vec::new())
    }
    pub fn load(data: &'a [u8]) -> Self {
        Self::Load(Reader { data, cursor: 0 })
    }
    pub fn loading(&self) -> bool {
        matches!(self, Self::Load(_))
    }
    pub fn value<T: Codec>(&mut self, v: &mut T) -> Result<(), ScfError> {
        match self {
            Self::Save(o) => v.put(o),
            Self::Load(r) => {
                *v = T::get(r)?;
                Ok(())
            }
        }
    }
    pub fn same<T: Codec>(&mut self, v: &T) -> Result<(), ScfError> {
        match self {
            Self::Save(o) => v.put(o),
            Self::Load(r) => {
                let saved = T::get(r)?;
                let mut a = Vec::new();
                let mut b = Vec::new();
                v.put(&mut a)?;
                saved.put(&mut b)?;
                if a != b {
                    return Err(invalid("model/options/layout mismatch"));
                }
                Ok(())
            }
        }
    }
    pub fn finish(self) -> Result<Vec<u8>, ScfError> {
        match self {
            Self::Save(o) => {
                if o.len() > MAX_BYTES {
                    Err(invalid("state quota"))
                } else {
                    Ok(o)
                }
            }
            Self::Load(r) => {
                if r.cursor != r.data.len() {
                    Err(invalid("trailing state bytes"))
                } else {
                    Ok(Vec::new())
                }
            }
        }
    }
    pub fn finite(&self, v: &[f64]) -> Result<(), ScfError> {
        if v.iter().any(|v| !v.is_finite()) {
            Err(invalid("nonfinite numerical state"))
        } else {
            Ok(())
        }
    }
    pub fn real(&mut self, v: &mut RealMatrix) -> Result<(), ScfError> {
        let shape = [v.rows(), v.cols()];
        self.value(v)?;
        if shape != [v.rows(), v.cols()] {
            return Err(invalid("orbital shape"));
        }
        self.finite(v.as_slice())
    }
    pub fn complex(&mut self, v: &mut ComplexMatrix) -> Result<(), ScfError> {
        let shape = [v.rows(), v.cols()];
        self.value(v)?;
        if shape != [v.rows(), v.cols()]
            || v.as_slice()
                .iter()
                .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err(invalid("complex orbital shape/values"));
        }
        Ok(())
    }
    pub fn planes(&mut self, v: &mut Vec<f64>) -> Result<(), ScfError> {
        let n = v.len();
        self.value(v)?;
        if n != v.len() {
            return Err(invalid("density plane shape"));
        }
        self.finite(v)
    }
    pub fn reciprocal(&mut self, v: &mut Vec<Complex64>) -> Result<(), ScfError> {
        let n = v.len();
        self.value(v)?;
        if n != v.len() || v.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) {
            return Err(invalid("reciprocal shape/values"));
        }
        Ok(())
    }
}
record!(crate::energy::EnergyAudit {
    xc_floor,
    nlcc,
    solvation,
    kinetic,
    local,
    local_zero_mode,
    nonlocal,
    hartree,
    xc,
    ewald,
    band_sum,
    double_counting,
    band_form,
    harris_foulkes,
    functional_difference_bound
});
record!(crate::energy::SolvationEnergyLedger {
    electrostatic_ha,
    background_correction_ha,
    cavitation_ha,
    total_ha,
    potential_integral_ha,
    iterations,
    relative_residual,
    fft_count
});
record!(crate::energy::XcFloorAudit {
    density_floor,
    input_negative_charge,
    output_negative_charge
});
record!(crate::energy::CoreXcLedger {
    core_electrons,
    core_potential_integral_ha,
    valence_potential_integral_ha
});
record!(crate::model::PhaseTimings {
    potential,
    diagonalization,
    density,
    energy,
    mixing
});
record!(crate::model::FixedPointCertificate {
    strict_residuals_per_k,
    buffer_residuals_per_k,
    density_residual,
    occupied_eigenvalue_change,
    davidson_iterations,
    operator_applications,
    davidson_iterations_per_k,
    operator_applications_per_k,
    timings
});
record!(crate::model::VerificationAttempt {
    iteration,
    passed,
    coulomb_error_ha,
    threshold_ha,
    eigen_tolerance_ha,
    residual_bound_ha,
    density_residual,
    occupied_eigenvalue_change_ha,
    strict_criteria_met,
    davidson_iterations,
    operator_applications,
    wall,
    spin_density_residual,
    occupation_change_electrons,
    chemical_potential_change_ha,
    moment_change
});
record!(super::record::CorrectionRecord {
    active_bands,
    basis_before,
    h_columns,
    reduced_order,
    repairs,
    dropped_corrections,
    restarts
});
record!(super::record::SolveRecord {
    ethr_ry,
    first_step,
    iterations,
    active_bands,
    h_columns,
    reduced_orders,
    restarts,
    unconverged,
    ritz_changes_ha,
    inner_clocks,
    ritz_values_ha,
    band_thresholds_ha,
    corrections,
    repairs,
    dropped_corrections
});
record!(super::record::VerificationCost {
    gate_enabled,
    first_gate_enabled,
    first_prediction_ha,
    contraction_prediction_ha,
    skipped,
    predicted_error_ha,
    actual_error_ha,
    reused_from
});
record!(super::record::VerificationEnergy {
    delivered_ha,
    checked_ha
});
record!(super::record::Record {
    k_call_indices,
    warm_start,
    stage,
    iteration,
    verification_cost,
    verification_energy,
    verification_warm,
    map_start_s,
    estimator_s,
    verification_start_s,
    verification_end_s,
    map_end_s,
    dr2_ry,
    ethr_requested_ry,
    ethr_ry,
    first_iteration_redo,
    beta,
    mixing,
    calls,
    verification_calls,
    davidson_iterations,
    h_columns,
    map_s,
    solve_s,
    density_s,
    potential_s,
    energy_s,
    density_transform_s,
    metric_s,
    mixing_s,
    correction_s,
    verification_s,
    descf_ha,
    variational_energy_ha,
    energy,
    variational_hartree_ha,
    variational_xc_ha
});
record!(super::broyden::Diagnostics {
    iteration,
    occupancy,
    inserted_slot,
    gamma_l2,
    gamma_max,
    condition
});
record!(super::broyden::ConditioningDiagnostics {
    dropped_pairs,
    dropped_conditions
});
record!(super::electronic::Iteration {
    iteration,
    dr2_ry,
    ethr_ry,
    first_iteration_redo,
    energy_ha,
    entropy_term_ha,
    free_energy_ha,
    moment_mu_b,
    residual_guard_ha,
    mixing,
    history_conditioning,
    mixing_fraction,
    calls
});
impl Codec for crate::model::IterationRecord {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        if self.pulay.is_some() {
            return Err(invalid("legacy iteration fields"));
        }
        self.standard_error_ha.put(o)?;
        self.standard_redo.put(o)?;
        self.iteration.put(o)?;
        self.energy.put(o)?;
        self.energy_change.put(o)?;
        self.density_residual.put(o)?;
        self.max_eigen_residual.put(o)?;
        self.requested_eigen_tolerance.put(o)?;
        self.davidson_iterations.put(o)?;
        self.operator_applications.put(o)?;
        self.davidson_iterations_per_k.put(o)?;
        self.operator_applications_per_k.put(o)?;
        self.charge_error.put(o)?;
        self.input_charge_error.put(o)?;
        self.min_density.put(o)?;
        self.mixing_fraction.put(o)?;
        self.mixer_reset.put(o)?;
        self.density_contraction.put(o)?;
        self.mixing_history_used.put(o)?;
        self.mixing_resets.put(o)?;
        self.mixing_raw_recovery.put(o)?;
        self.fixed_point_certificate.put(o)?;
        self.timings.put(o)?;
        Ok(())
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        Ok(Self {
            standard_error_ha: Codec::get(r)?,
            standard_redo: Codec::get(r)?,
            iteration: Codec::get(r)?,
            energy: Codec::get(r)?,
            energy_change: Codec::get(r)?,
            density_residual: Codec::get(r)?,
            max_eigen_residual: Codec::get(r)?,
            requested_eigen_tolerance: Codec::get(r)?,
            davidson_iterations: Codec::get(r)?,
            operator_applications: Codec::get(r)?,
            davidson_iterations_per_k: Codec::get(r)?,
            operator_applications_per_k: Codec::get(r)?,
            charge_error: Codec::get(r)?,
            input_charge_error: Codec::get(r)?,
            min_density: Codec::get(r)?,
            mixing_fraction: Codec::get(r)?,
            mixer_reset: Codec::get(r)?,
            density_contraction: Codec::get(r)?,
            mixing_history_used: Codec::get(r)?,
            mixing_resets: Codec::get(r)?,
            mixing_raw_recovery: Codec::get(r)?,
            fixed_point_certificate: Codec::get(r)?,
            timings: Codec::get(r)?,
            pulay: None,
        })
    }
}
impl Codec for crate::model::VerificationReport {
    fn put(&self, o: &mut Vec<u8>) -> Result<(), ScfError> {
        if self.schedule.is_some() {
            return Err(invalid("legacy verification schedule"));
        }
        self.attempts.put(o)
    }
    fn get(r: &mut Reader<'_>) -> Result<Self, ScfError> {
        Ok(Self {
            attempts: Codec::get(r)?,
            schedule: None,
        })
    }
}

pub(super) fn header<K: dft_mpi::Collective>(
    io: &mut Exchange<'_>,
    input: &crate::refstart::StartInput<'_>,
    profile: super::profile::Profile,
    comm: &K,
    kind: &str,
    dims: [usize; 3],
) -> Result<(), ScfError> {
    io.same(&SCHEMA_VERSION.to_string())?;
    io.same(&kind.to_string())?;
    io.same(&comm.rank())?;
    io.same(&comm.size())?;
    io.same(&dims)?;
    io.same(&input.cell.lattice_columns_bohr())?;
    io.same(&format!("{:?}", input.atoms))?;
    io.same(&format!("{:?}", input.options))?;
    io.same(&format!("{profile:?}"))?;
    io.same(&dft_numerics::switches::RunConfig::current().agreement_hash())
}
pub(super) fn gamma_operator(
    io: &mut Exchange<'_>,
    state: &mut crate::refstart::StartState,
) -> Result<(), ScfError> {
    io.same(&state.operator.local_rows().to_vec())?;
    io.same(&state.ionic)?;
    io.same(&state.diagonal_static)?;
    io.same(&state.electrons)?;
    let mut epoch = state.operator.owned_epoch();
    io.value(&mut epoch)?;
    if io.loading()
        && let Some(epoch) = epoch
    {
        // The next map installs its entire potential before H is used. Restore
        // the epoch now so the retained-HX eligibility check is unchanged.
        state
            .operator
            .set_owned_potential(&vec![0.; state.operator.plane_points()], epoch)
            .map_err(super::numerical)?;
    }
    io.real(&mut state.orbitals)
}
pub(super) fn k_state(
    io: &mut Exchange<'_>,
    state: &mut crate::refstart::KStartState,
) -> Result<(), ScfError> {
    io.same(&state.blocks.len())?;
    // GK9 k pools: a rank's checkpoint holds its own whole k points; the
    // switch is in the agreement hash, the owned range is checked here.
    if let Some(pool) = &state.pool {
        let (j, i, p) = pool.position();
        io.same(&[
            pool.local().start,
            pool.local().end,
            pool.nk,
            pool.npool(),
            j,
            i,
            p,
        ])?;
    }
    io.same(&state.ionic)?;
    io.same(&state.electrons)?;
    for (block, x) in state.blocks.iter_mut().zip(&mut state.orbitals.0) {
        io.same(&block.point.fractional)?;
        io.same(&block.point.weight)?;
        io.same(&block.operator.local_rows().to_vec())?;
        io.same(&block.operator.diagonal_static().to_vec())?;
        let mut epoch = block.operator.owned_epoch();
        io.value(&mut epoch)?;
        if io.loading()
            && let Some(epoch) = epoch
        {
            if block.operator.owned_epoch().is_none_or(|old| epoch > old) {
                block
                    .operator
                    .set_owned_potential(&vec![0.; block.operator.plane_points()], epoch)
                    .map_err(super::numerical)?;
            } else if block.operator.owned_epoch() != Some(epoch) {
                return Err(invalid("operator epoch"));
            }
        }
        io.complex(x)?;
    }
    Ok(())
}
pub(super) fn counters(
    iteration: usize,
    budget: usize,
    history: usize,
    ethr: f64,
    previous: f64,
) -> Result<(), ScfError> {
    if iteration >= budget
        || history != iteration
        || !ethr.is_finite()
        || ethr <= 0.
        || ethr > 1e-2
        || (iteration > 0 && (!previous.is_finite() || previous < 0.))
    {
        return Err(invalid("iteration/history/ethr state"));
    }
    Ok(())
}
record!(super::symmetry::Residual {
    iteration,
    verification,
    norm,
    pre_projection_relative_norm,
    post_projection_invariance_norm
});
record!(super::symmetry::Report {
    group_order,
    irreducible_k_count,
    channels,
    maps
});
pub(super) fn symmetry(
    io: &mut Exchange<'_>,
    s: &mut Option<super::symmetry::Workspace>,
) -> Result<(), ScfError> {
    io.same(&s.is_some())?;
    if let Some(s) = s {
        io.same(&s.report.group_order)?;
        io.same(&s.report.irreducible_k_count)?;
        io.same(&s.report.channels)?;
        io.value(&mut s.report)?;
        io.value(&mut s.iteration)?;
        io.value(&mut s.verification)?;
    }
    Ok(())
}
pub(super) fn cached_fields(
    cache: &Option<(Vec<f64>, super::map::Fields)>,
    planes: usize,
) -> Result<(), ScfError> {
    if let Some((n, fields)) = cache {
        if n.len() != planes || n.iter().any(|v| !v.is_finite()) {
            return Err(invalid("cached density shape/values"));
        }
        fields.validate_checkpoint(planes)?;
    }
    Ok(())
}
