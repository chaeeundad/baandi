# GitHub issue #23: Research: measure current-state discretization error before changing cutoff-ladder handoff

Reviewed at `1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4` (2026-11-17). This is a static review of the implementation and committed experiment records; no new BAANDI/QE timing run was performed in this review. Related tracking issue: #13.

**Suggested priority:** P3 — calibration before changing the adopted ladder.

### Preserve what the new evidence established

The ledger's current-stack real-MPI ladder experiment reports 3.480 s with the ladder versus 4.083 s without it on the screened Si32 case. Keep the adopted ladder. Do not repeat the earlier “remove the coarse stage” proposal.

[cutoff_ladder.rs](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-engine/src/cutoff_ladder.rs#L34) currently uses a fixed coarse residual threshold and offers an atomic-density truncation alternative. The research question is whether the **current orbitals'** unresolved fine-basis error can predict a better handoff or a better initial fine state across materials. This differs from a fixed fraction/threshold sweep or the existing atomic-start proxy.

### Research basis

Bordignon et al., **SIAM J. Sci. Comput. 47, A2881–A2905 (2025)**, [DOI](https://doi.org/10.1137/25M1735676), [full author manuscript](https://arxiv.org/html/2409.11769v2), separate discretization and SCF contributions to energy error for periodic Kohn–Sham models. Their rigorous guarantee relies on convex functionals and spectral conditions. For practical nonconvex LDA/PBE, they demonstrate useful estimates, not the same general guarantee. Approximate inverse-operator evaluations also introduce additional qualifications.

Use this as motivation for an error-budget experiment; do not label a BAANDI LDA/PBE heuristic “certified” by citing the paper.

### Proposed experiment, before a new controller

1. Reuse validated snapshots at candidate coarse stopping points. Prolong X to the **same final cutoff/grid** and evaluate the residual component outside the coarse orbital basis. Measure where the handoff defect lies: high-G discretization error, low-G density-response error, or poorly resolved empty states.
2. Compute a cheap kinetic-inverse-weighted residual indicator and calibrate it against an accurately completed fine solve. Compare energy-like quantities in a common norm/units; a coarse density-residual norm cannot simply be compared numerically with an energy error indicator.
3. Estimate the cost of the fine-grid probe itself. Start with offline data and at most a sparse trigger near a candidate handoff; evaluating a full fine H block at every coarse iteration may cost more than it saves.
4. Only if predictive, test a policy that balances remaining coarse-SCF work against current-state discretization error and observed fine-recovery cost. Keep the final cutoff, FFT grid, PP, final tolerance and verification unchanged.
5. Optional second experiment: approximate a high-G response correction to the prolonged orbitals. Compare against spending the **same H/FFT budget** on the existing first fine Davidson correction. A diagonal high-G correction that merely duplicates that step is not a new gain.

### Counterevidence and safeguards

BT3-COARSE already shows that cheaper coarse work can be repaid by more fine work. LADDER-SECANT had material-dependent failures. A candidate must beat the adopted ladder on **coarse + handoff + all fine + verification** time, not just coarse iterations. Include a harder/non-Si PP case from the supported validation set; do not learn a Si-only handoff constant.

Success yields a physically informed initialization strategy. Failure yields a validated indicator and a closed experiment, not another default switch.

### Validation and adoption

Use the repository's current two-lane fast-loop protocol. Start with one short case and an identity/correctness smoke; use full-SCF Stage 1 for any trajectory change. Final acceptance uses the designated MPI benchmark host, the current default stack, fixed-point verification included, the existing accuracy classes and same-model/same-rank QE comparisons. Measurable small improvements remain eligible under the adoption rule; do not impose a new minimum percentage. Report failed or negative experiments as well as wins.
