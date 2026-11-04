//! Optional reciprocal-space residual filtering; the legacy mixers are unchanged.
use crate::{
    SolverError,
    mixing::{DensityMixer, PulayMixer},
};
use num_complex::Complex64;

/// In-place transform supplied by the FFT adapter. `inverse=false` means
/// real -> reciprocal (1/N normalization); `true` means reciprocal -> real.
/// Storage is C order, with the last coordinate contiguous. Implementations
/// must reuse scratch and overwrite every element, or return an error.
pub trait DensityTransform {
    fn transform(&mut self, values: &mut [Complex64], inverse: bool) -> Result<(), SolverError>;
}
impl<F> DensityTransform for F
where
    F: FnMut(&mut [Complex64], bool) -> Result<(), SolverError>,
{
    fn transform(&mut self, values: &mut [Complex64], inverse: bool) -> Result<(), SolverError> {
        self(values, inverse)
    }
}

/// A switchable residual operation. The input is the *raw* output-input residual;
/// the destination has identical length. Filtering must not decide convergence.
pub trait ResidualPreconditioner {
    /// True only when reconstruction from immutable settings captures all
    /// arithmetic state. Arbitrary external preconditioners fail closed.
    fn checkpoint_stateless(&self) -> bool {
        false
    }

    fn apply(&mut self, residual: &[f64], filtered: &mut [f64]) -> Result<(), SolverError>;
}

/// Kerker multiplier |G|^2/(|G|^2+q0^2), q0 in Bohr^-1, G=B*m.
/// B is indexed [Cartesian component][reciprocal column], including 2*pi.
/// The zero Fourier coefficient is set to exactly zero. On even grids, filters
/// are averaged over all signs of Nyquist coordinates to preserve conjugacy
/// even in skew cells (those aliased samples have no unique signed G).
pub struct Kerker<T> {
    transform: T,
    weights: Vec<f64>,
    scratch: Vec<Complex64>,
}
impl<T: DensityTransform> Kerker<T> {
    pub fn new(
        shape: [usize; 3],
        reciprocal: [[f64; 3]; 3],
        q0: f64,
        transform: T,
    ) -> Result<Self, SolverError> {
        let invalid = || {
            SolverError::InvalidParameter(
                "nonzero grid, nonsingular finite reciprocal cell and positive finite q0 required",
            )
        };
        let size = shape
            .iter()
            .try_fold(1usize, |a, &b| a.checked_mul(b))
            .ok_or_else(invalid)?;
        if size == 0
            || size > isize::MAX as usize / std::mem::size_of::<Complex64>()
            || !q0.is_finite()
            || q0 <= 0.0
            || !(q0 * q0).is_finite()
            || q0 * q0 == 0.0
            || reciprocal.iter().flatten().any(|x| !x.is_finite())
        {
            return Err(invalid());
        }
        let b = reciprocal;
        let determinant = b[0][0] * (b[1][1] * b[2][2] - b[1][2] * b[2][1])
            - b[0][1] * (b[1][0] * b[2][2] - b[1][2] * b[2][0])
            + b[0][2] * (b[1][0] * b[2][1] - b[1][1] * b[2][0]);
        if !determinant.is_finite() || determinant == 0.0 {
            return Err(invalid());
        }
        let mut weights = Vec::with_capacity(size);
        for i in 0..shape[0] {
            for j in 0..shape[1] {
                for k in 0..shape[2] {
                    let indices = [i, j, k];
                    let mut m = [0.0; 3];
                    let mut nyquist = Vec::new();
                    for a in 0..3 {
                        m[a] = if indices[a] <= shape[a] / 2 {
                            indices[a] as f64
                        } else {
                            indices[a] as f64 - shape[a] as f64
                        };
                        if shape[a].is_multiple_of(2) && indices[a] == shape[a] / 2 {
                            nyquist.push(a);
                        }
                    }
                    let count = 1usize << nyquist.len();
                    let mut weight = 0.0;
                    for signs in 0..count {
                        let mut signed = m;
                        for (bit, &axis) in nyquist.iter().enumerate() {
                            if signs & (1 << bit) != 0 {
                                signed[axis] = -signed[axis];
                            }
                        }
                        let g2: f64 = b
                            .iter()
                            .map(|row| {
                                row.iter()
                                    .zip(signed)
                                    .map(|(v, n)| v * n)
                                    .sum::<f64>()
                                    .powi(2)
                            })
                            .sum();
                        if !g2.is_finite() || !(g2 + q0 * q0).is_finite() {
                            return Err(invalid());
                        }
                        weight += g2 / (g2 + q0 * q0) / count as f64;
                    }
                    weights.push(weight);
                }
            }
        }
        weights[0] = 0.0;
        Ok(Self {
            transform,
            weights,
            scratch: vec![Complex64::default(); size],
        })
    }
    pub fn weights(&self) -> &[f64] {
        &self.weights
    }
}
impl<T: DensityTransform> ResidualPreconditioner for Kerker<T> {
    fn checkpoint_stateless(&self) -> bool {
        true
    }

    fn apply(&mut self, residual: &[f64], filtered: &mut [f64]) -> Result<(), SolverError> {
        if residual.len() != self.weights.len()
            || filtered.len() != residual.len()
            || residual.iter().any(|x| !x.is_finite())
        {
            return Err(SolverError::InvalidParameter(
                "finite residual on the configured grid required",
            ));
        }
        for (v, &r) in self.scratch.iter_mut().zip(residual) {
            *v = r.into();
        }
        self.transform.transform(&mut self.scratch, false)?;
        for (v, &w) in self.scratch.iter_mut().zip(&self.weights) {
            *v *= w;
        }
        self.scratch[0] = Complex64::default();
        self.transform.transform(&mut self.scratch, true)?;
        let scale = self.scratch.iter().map(|x| x.re.abs()).fold(1.0, f64::max);
        if self
            .scratch
            .iter()
            .any(|x| !x.re.is_finite() || !x.im.is_finite() || x.im.abs() > 1e-12 * scale)
        {
            return Err(SolverError::Application(
                "nonfinite or non-real preconditioned residual".into(),
            ));
        }
        for (out, v) in filtered.iter_mut().zip(&self.scratch) {
            *out = v.re;
        }
        Ok(())
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct MixingDiagnostics {
    /// Euclidean raw density residual, in electrons/Bohr^3 (no quadrature).
    pub residual_norm: f64,
    /// Raw residual norm divided by its previous value; None after a reset
    /// or when the previous residual is zero. Never a convergence certificate.
    pub contraction: Option<f64>,
    /// Number of differences presented to Pulay (not retained numerical rank).
    pub history_used: usize,
    pub resets: usize,
    /// Count of positivity recovery steps, retained across explicit resets.
    pub raw_recoveries: usize,
    pub used_raw_recovery: bool,
}

/// Type-II Anderson applied to the filtered fixed-point map x+P(F(x)-x).
/// Least squares therefore minimizes ||P R - Delta(P R) gamma||_2.
/// An infeasible history-free step uses raw convex recovery; other feasibility
/// damping, divergence and physical convergence remain the controller's responsibility. The existing Pulay implementation owns history and its
/// spectral regularization; this wrapper does not duplicate that algorithm.
pub struct PreconditionedPulay<P> {
    pulay: PulayMixer,
    recovery: crate::mixing::LinearMixer,
    preconditioner: P,
    residual: Vec<f64>,
    mapped: Vec<f64>,
    mean: Option<f64>,
    previous: Option<f64>,
    diagnostics: MixingDiagnostics,
}
impl<P: ResidualPreconditioner> PreconditionedPulay<P> {
    pub fn new(
        beta: f64,
        history: usize,
        size: usize,
        preconditioner: P,
    ) -> Result<Self, SolverError> {
        if size == 0 || size > isize::MAX as usize / std::mem::size_of::<f64>() {
            return Err(SolverError::InvalidParameter(
                "representable nonempty density grid required",
            ));
        }
        Ok(Self {
            pulay: PulayMixer::new(beta, history)?,
            recovery: crate::mixing::LinearMixer {
                beta: beta.min(1.0),
            },
            preconditioner,
            residual: vec![0.0; size],
            mapped: vec![0.0; size],
            mean: None,
            previous: None,
            diagnostics: MixingDiagnostics::default(),
        })
    }
    pub fn diagnostics(&self) -> MixingDiagnostics {
        self.diagnostics
    }
}
impl<P: ResidualPreconditioner> DensityMixer for PreconditionedPulay<P> {
    fn reset(&mut self) {
        self.pulay.reset();
        self.mean = None;
        self.previous = None;
        self.diagnostics = MixingDiagnostics {
            resets: self.diagnostics.resets + 1,
            raw_recoveries: self.diagnostics.raw_recoveries,
            ..Default::default()
        };
    }
    fn mix(&mut self, input: &[f64], output: &[f64], mean: f64) -> Result<Vec<f64>, SolverError> {
        if input.len() != self.residual.len()
            || output.len() != input.len()
            || !mean.is_finite()
            || mean <= 0.0
            || input.iter().chain(output).any(|x| !x.is_finite())
        {
            return Err(SolverError::InvalidParameter(
                "finite densities on configured grid and positive mean required",
            ));
        }
        if self.mean.is_some_and(|m| m != mean) {
            self.reset();
        }
        for ((r, &x), &y) in self.residual.iter_mut().zip(input).zip(output) {
            *r = y - x;
        }
        let norm = self.residual.iter().fold(0.0_f64, |a, &b| a.hypot(b));
        if !norm.is_finite() {
            return Err(dft_numerics::linalg::LinalgError::NonFinite.into());
        }
        self.preconditioner
            .apply(&self.residual, &mut self.mapped)?;
        for (y, &x) in self.mapped.iter_mut().zip(input) {
            *y += x;
        }
        let mut result = self.pulay.mix(input, &self.mapped, mean)?;
        // A reciprocal filter is not positivity preserving. With no secant
        // history, repeated feasibility damping can trap the iterate on the
        // boundary. A convex raw-density step leaves that boundary without
        // clipping samples or changing the requested charge. The sampled pair
        // remains valid history even though a different step was accepted.
        self.diagnostics.used_raw_recovery =
            self.pulay.history_len() == 1 && result.iter().any(|v| *v < 0.0);
        if self.diagnostics.used_raw_recovery {
            result = self.recovery.mix(input, output, mean)?;
            self.diagnostics.raw_recoveries += 1;
        }
        self.mean = Some(mean);
        self.diagnostics.residual_norm = norm;
        self.diagnostics.contraction = self.previous.filter(|v| *v > 0.0).map(|v| norm / v);
        self.diagnostics.history_used = self.pulay.history_len().saturating_sub(1);
        self.previous = Some(norm);
        Ok(result)
    }
}

/// Read-only history introspection without changing the legacy mixer interface.
/// A returned length counts stored input/residual pairs, not fitted rank.
/// Version 1 continuation memory. FFT scratch is overwritten on every apply;
/// the host reconstructs the identical Kerker filter from the model and options.
#[derive(Clone, Debug)]
pub enum MixerState {
    Linear {
        beta: f64,
    },
    Pulay(crate::mixing::PulayState),
    Kerker {
        pulay: crate::mixing::PulayState,
        mean: Option<f64>,
        previous: Option<f64>,
        diagnostics: MixingDiagnostics,
    },
}
pub trait MixerHistory: DensityMixer {
    /// Custom stateful preconditioners must explicitly opt in; never silently
    /// discard opaque memory when an exact continuation is requested.
    fn checkpoint(&self) -> Result<MixerState, SolverError> {
        Err(SolverError::InvalidParameter("opaque mixer state"))
    }
    fn restore(&mut self, _state: &MixerState, _size: usize) -> Result<(), SolverError> {
        Err(SolverError::InvalidParameter("opaque mixer state"))
    }

    fn history_len(&self) -> usize;
    fn used_raw_recovery(&self) -> bool {
        false
    }
}
impl MixerHistory for crate::mixing::LinearMixer {
    fn checkpoint(&self) -> Result<MixerState, SolverError> {
        Ok(MixerState::Linear { beta: self.beta })
    }
    fn restore(&mut self, state: &MixerState, _: usize) -> Result<(), SolverError> {
        match state {
            MixerState::Linear { beta } if beta.to_bits() == self.beta.to_bits() => Ok(()),
            _ => Err(SolverError::InvalidParameter("mixer checkpoint mismatch")),
        }
    }

    fn history_len(&self) -> usize {
        0
    }
}
impl MixerHistory for PulayMixer {
    fn checkpoint(&self) -> Result<MixerState, SolverError> {
        Ok(MixerState::Pulay(self.checkpoint()))
    }
    fn restore(&mut self, state: &MixerState, size: usize) -> Result<(), SolverError> {
        match state {
            MixerState::Pulay(p) => self.restore(p, size),
            _ => Err(SolverError::InvalidParameter("mixer checkpoint mismatch")),
        }
    }

    fn history_len(&self) -> usize {
        self.history_len()
    }
}
impl<P: ResidualPreconditioner> MixerHistory for PreconditionedPulay<P> {
    fn checkpoint(&self) -> Result<MixerState, SolverError> {
        if !self.preconditioner.checkpoint_stateless() {
            return Err(SolverError::InvalidParameter("opaque preconditioner state"));
        }
        Ok(MixerState::Kerker {
            pulay: self.pulay.checkpoint(),
            mean: self.mean,
            previous: self.previous,
            diagnostics: self.diagnostics,
        })
    }
    fn restore(&mut self, state: &MixerState, size: usize) -> Result<(), SolverError> {
        let MixerState::Kerker {
            pulay,
            mean,
            previous,
            diagnostics,
        } = state
        else {
            return Err(SolverError::InvalidParameter("mixer checkpoint mismatch"));
        };
        if !self.preconditioner.checkpoint_stateless()
            || size != self.residual.len()
            || mean.is_some_and(|v| !v.is_finite() || v <= 0.0)
            || previous.is_some_and(|v| !v.is_finite() || v < 0.0)
            || !diagnostics.residual_norm.is_finite()
            || diagnostics.residual_norm < 0.0
            || diagnostics
                .contraction
                .is_some_and(|v| !v.is_finite() || v < 0.0)
            || mean != &pulay.mean
            || mean.is_some() != previous.is_some()
            || diagnostics.history_used > pulay.capacity
            || diagnostics.resets == usize::MAX
            || diagnostics.raw_recoveries == usize::MAX
        {
            return Err(SolverError::InvalidParameter("invalid Kerker checkpoint"));
        }
        self.pulay.restore(pulay, size)?;
        self.mean = *mean;
        self.previous = *previous;
        self.diagnostics = *diagnostics;
        Ok(())
    }

    fn used_raw_recovery(&self) -> bool {
        self.diagnostics.used_raw_recovery
    }
    fn history_len(&self) -> usize {
        self.pulay.history_len()
    }
}
