//! Bounded, plan-local lane observations; no clocks, allocation or I/O.

/// All counters describe the adapter boundary, not opaque MKL internals.
/// `thread_limit_*` is queried from the committed descriptor on every profiled
/// call. A limit of one rules out multiple FFT workers; it is not a sampler of
/// OS threads. MXCSR is observed without changing the floating-point mode.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct LaneDiagnostics {
    /// Forward, backward DFTI calls, including failed computations.
    pub calls: [u64; 2],
    pub points: u64,
    /// Source bytes explicitly copied by this adapter, excluding FFT traffic.
    /// The direct caller-buffer path has no copy sites, so this stays zero.
    pub explicit_copy_bytes: u64,
    pub failed_calls: u64,
    pub thread_limit_min: u64,
    pub thread_limit_max: u64,
    pub mxcsr_samples: u64,
    pub mxcsr_control_changes: u64,
    pub ftz_off_calls: u64,
    pub daz_off_calls: u64,
    /// Calls that newly set the sticky denormal-operand flag. Not an operand
    /// count: an already-set flag cannot identify a subsequent denormal.
    pub new_denormal_flag_calls: u64,
    pub mxcsr_before_or: u32,
    pub mxcsr_after_or: u32,
}

impl LaneDiagnostics {
    pub(crate) fn record(
        &mut self,
        inverse: bool,
        points: usize,
        limit: u64,
        before: Option<u32>,
        after: Option<u32>,
        failed: bool,
    ) {
        if self.calls == [0; 2] {
            self.thread_limit_min = limit;
        }
        self.calls[usize::from(inverse)] += 1;
        self.points += points as u64;
        self.failed_calls += u64::from(failed);
        self.thread_limit_min = self.thread_limit_min.min(limit);
        self.thread_limit_max = self.thread_limit_max.max(limit);
        if let (Some(before), Some(after)) = (before, after) {
            self.mxcsr_samples += 1;
            // Ignore sticky exception flags (bits 0..5); include DAZ, exception
            // masks, rounding mode and FTZ. Do not clear the caller's flags.
            self.mxcsr_control_changes += u64::from((before ^ after) & 0xffc0 != 0);
            self.ftz_off_calls += u64::from(before & (1 << 15) == 0);
            self.daz_off_calls += u64::from(before & (1 << 6) == 0);
            self.new_denormal_flag_calls += u64::from(before & 2 == 0 && after & 2 != 0);
            self.mxcsr_before_or |= before;
            self.mxcsr_after_or |= after;
        }
    }
}

#[cfg(target_arch = "x86_64")]
pub(crate) fn mxcsr() -> Option<u32> {
    let mut value = 0u32;
    // SAFETY: SSE is mandatory on x86-64. STMXCSR writes exactly four bytes to
    // this live stack slot and changes no FP control bits or status flags.
    unsafe {
        std::arch::asm!("stmxcsr [{p}]", p = in(reg) &mut value,
            options(nostack, preserves_flags));
    }
    Some(value)
}

#[cfg(not(target_arch = "x86_64"))]
pub(crate) fn mxcsr() -> Option<u32> {
    None // Apple Silicon type-checks the MKL adapter but cannot execute it.
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn counts_control_changes_separately_from_sticky_flags() {
        let mut d = LaneDiagnostics::default();
        d.record(false, 60, 1, Some(0x1f80), Some(0x1fa0), false);
        d.record(true, 60, 1, Some(0x1f80), Some(0x9fc2), true);
        assert_eq!(d.calls, [1, 1]);
        assert_eq!(d.points, 120);
        assert_eq!(d.explicit_copy_bytes, 0);
        assert_eq!(d.failed_calls, 1);
        assert_eq!((d.thread_limit_min, d.thread_limit_max), (1, 1));
        assert_eq!(d.mxcsr_samples, 2);
        assert_eq!(d.mxcsr_control_changes, 1);
        assert_eq!((d.ftz_off_calls, d.daz_off_calls), (2, 2));
        assert_eq!(d.new_denormal_flag_calls, 1);
        d.record(false, 1, 2, None, None, false);
        assert_eq!(d.thread_limit_max, 2);
        assert_eq!(d.mxcsr_samples, 2);
    }

    #[test]
    fn observing_mxcsr_preserves_it() {
        assert_eq!(mxcsr(), mxcsr());
    }
}
