# Phase-0 integration record (completed)

The orchestrator explicitly approved instrumentation-only changes to
`stage_clock.rs`, `fused_pair_plan.rs`, `dft-engine/src/profile.rs`, and
`dft-cli/src/run.rs` on 2026-11-20.

- Serial and distributed plans share the default-off probe in `real.rs`.
- `stage_clock::local::Snapshot::fft` carries thread-local deltas.
- `profile::Accum::fft` and the separate `FftRankProfile` retain bucket/rank
  identity; existing seconds/count slot lists and wall-closure math are unchanged.
- The existing final gather carries extra words only with FFT_PROFILE enabled.
  Two 32-bit limbs encode each u64 exactly; a 3-rank test exercises values above
  2^53 and verifies rank-specific map counts.
- `fft_lane_json` is one localized serializer section in the CLI; the parser
  consumes actual CLI output, not a manually assembled surrogate.

Scope is one thread per rank, paired Gamma local-potential kernels and kept-grid
potential changes. Density/Hartree transforms are excluded. No numerical
algorithm, default, tolerance, FFT layout or communication in the calculation
was changed. See profile-request.md for exact field names and coverage limits.
