# GitHub issue #17: Fix first-block diagnostics: occupied count, electron metadata, and Ritz/residual semantics

Reviewed at `1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4` (2026-11-17). This is a static review of the implementation and committed experiment records; no new BAANDI/QE timing run was performed in this review. Related tracking issue: #13.

**Suggested priority:** P1 — diagnostic validity.

### Problem and evidence

Two diagnostic assumptions currently undermine the first-block optimization decision.

1. [gamma_real.rs:1051](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-engine/src/gamma_real.rs#L1051) sets `policy_occupied = policy.strict_bands.min(bands)`. [scf.rs:2916](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-engine/src/scf.rs#L2916) sets `strict_bands = occupied_count + 1` on the standard path, capped by the number of bands. The mandatory prefix includes a first empty state; it is not the occupation count. That value is passed to the residual spectrum diagnostic and snapshot, and [gamma_real.rs:1107](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-engine/src/gamma_real.rs#L1107) writes `electrons = 2 * policy_occupied` (also in the distributed branches).
2. [first_block.rs:251](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-engine/src/first_block.rs#L251) compares true residual norms and singular values directly with the solve tolerance and labels the resulting counts `true_unconverged_*`. On the standard path that tolerance controls **Ritz-value change**, not residual norm. The diagnostic also hard-codes the standard prefix rule without carrying the convergence-criterion kind.

For neutral four-valence-electron Si32, the physical values are 128 electrons and 64 doubly occupied orbitals. The [committed diagnostic report](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/docs/performance/first-block-residual-2026-11-17.md) instead uses 65 “occupied” states out of 68. Consequently its residual SVD includes the first empty state, and its `theta[n_o] - theta[n_o-1]` gap is between the first two empty states, not HOMO and LUMO. The same offset appears for Si16 and Si8.

This is a diagnostic/snapshot-metadata defect. The production density path uses the actual occupations; this finding is **not evidence that SCF ran with two extra electrons**.

### Why the interpretation matters

The report concludes that a late occupied residual is dominated by one direction and that every direction remains unconverged under the solver's own criterion. Neither conclusion follows from the current diagnostic. The dominant direction might be the included empty state; the real occupied spectrum must be recomputed. A Ritz-change tolerance is not a residual-norm tolerance even though both can be expressed in Ha.

An independently evaluated 2x2 algebra example illustrates the scale mismatch: for `H=diag(0,1)` and `x=(cos(1e-4),sin(1e-4))`, Rayleigh-quotient error is approximately `1e-8`, while residual norm is approximately `1e-4`. A numerical threshold `1e-7` classifies them differently. This does not turn Ritz change into a rigorous error bound; it disproves treating the quantities as interchangeable.

### Proposed change

- Pass actual occupation count and electron count, preferably the occupation vector, from `ScfProblem` into the Gamma diagnostic/capture context. Keep mandatory solver roots, density-carrying roots and empty buffers as distinct fields.
- Carry an explicit convergence-criterion enum into the diagnostic. Report residuals and Ritz changes separately. Remove or rename counts that imply a standard-path convergence verdict from a residual comparison.
- For fixed integer occupations, form the occupied residual spectrum over exactly the occupied block and locate the physical occupied/empty boundary. Do not manufacture integer-occupation metadata for fractional/spin cases.
- Recompute the affected Si tables and qualify/replace Findings 1–2. Preserve the old record as superseded evidence rather than silently changing its meaning.
- Only then decide whether first-block selection/compression merits a new experiment. The correction does not by itself establish a speedup.

### Acceptance checks

- A fixture with `n_occ=2`, `strict_bands=3`, `bands=5` reports 4 electrons and the gap between indices 1 and 2.
- Make the first empty residual dominate the spectrum; it must not affect the occupied-only spectrum.
- Separate strict-residual and standard-Ritz diagnostic expectations, including the last mandatory root.
- Diagnostics on/off leave energy/density bits and solver work unchanged, apart from explicitly separated diagnostic work.
- Updated real-SCF evidence states the criterion and actual occupations.
