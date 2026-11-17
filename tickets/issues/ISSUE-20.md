# GitHub issue #20: Remove remaining replicated M2 construction after setup and fused-density default adoption

Reviewed against `1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4`, with the subsequent adoption commit `a430f8b7e86159c68ada328278b94666d509e82a` checked before finalizing this review (2026-11-17). Static review and committed measurements only; no new BAANDI/QE benchmark was run here. Related tracking issue: #13.

**Suggested priority:** P2 — remove remaining replicated construction after the newly adopted defaults.

### Update: the pending acceptance decision has now been made

Commit [a430f8b](https://github.com/chaeeundad/baandi/commit/a430f8b7e86159c68ada328278b94666d509e82a) landed during this review. It adopts:

- `M2_LIGHT_STATE=on`
- `EWALD_CACHE=on`
- `M2_DENSITY=fused`

`M2_LAYOUT_CACHE` remains off because the real-MPI experiment did not improve time. **Do not repeat the previously requested adoption screen for these completed decisions.** This issue now covers the remaining architectural work.

The new ledger row reports a fixed **10-step**, MPI np=8 screen on build `0fb09c0`, not a new full-convergence QE comparison:

| Arm | Committed observation |
|---|---|
| LIGHT_STATE alone | Si32 -4.47%; Si64 -5.19% |
| EWALD_CACHE alone | Si32 -1.98% |
| LAYOUT_CACHE alone | Si32 +0.23%, no time gain |
| All three setup switches | Si32 -6.40%; Si64 -7.83% / -7.96% |
| Fused density, separate experiments | -3.04% / -3.42% / -3.21%; density phase 0.2727→0.1690 s, with matching total-time reduction |

“All three setup switches” includes the layout switch that was not adopted. It is **not** a timing of the final default combination of light state + Ewald cache + fused density. Do not add or multiply these separate gains into a claimed new QE ratio. Treat the current hardened head and its full-convergence cost as the baseline for further work.

### Remaining work: remove redundant ownership

The [other/init profile](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/docs/performance/profile-2026-11-17-other-init.md) identifies work that remains even after light-state adoption.

1. Make the M2 state own a lightweight reciprocal-layout descriptor and the distributed operator. Construct a full serial real operator lazily only if fallback actually needs it. [gamma_real.rs:826](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-engine/src/gamma_real.rs#L826) still builds a serial `RealHamiltonian`; the adopted switch shrinks plans/scratch but does not remove its replicated geometry/projector conversion.
2. Construct owned kinetic/projector rows directly. [distributed_real_gamma.rs:262](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-physics/src/distributed_real_gamma.rs#L262) builds a full ORC kinetic vector and packs each projector into a full-length temporary before selecting local rows. Preserve the same coefficient mapping and arithmetic.
3. Separate geometry-only descriptors from density-dependent potential data and FFT-size-specific plans. Validate cell/position/charge, coarse/fine basis, layout and rank-count invalidation.
4. Treat distributed ionic-local construction and block atomic-start orthogonalization as subsequent independent experiments: the former can preserve per-G arithmetic, while the latter changes rounding and trajectory. Do not combine them into the first refactor.

### Expected scope

Incremental gains and lower per-rank memory/initialization cost. After the large plan/scratch waste has been removed, the remaining payoff may be small; measure it before expanding the refactor. This is not evidence that construction explains the historical 2–3x QE gap.

### Validation and adoption

Use the current two-lane fast-loop protocol and the **new defaults**. Start with one short case and an identity/correctness smoke. Include full SCF and verification for adoption, and Stage 1 for any trajectory change, with the existing accuracy/generality and same-model/same-rank comparisons. Retain measurable small gains under the existing adoption rule. Keep M2 fallback and invalidation tests; report negative results as well as wins.
