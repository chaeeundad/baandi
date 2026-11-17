# GitHub issue #18: Validate snapshot operator and row-layout identity before reusing captured H X

Reviewed at `1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4` (2026-11-17). This is a static review of the implementation and committed experiment records; no new BAANDI/QE timing run was performed in this review. Related tracking issue: #13.

**Suggested priority:** P1 — benchmark validity.

### Problem and trigger

The new replay harness can accept an incompatible snapshot and still arm its cached `H X`. This invalidates the invariant that cached products belong to the **current** Hamiltonian and current local-row interpretation.

- [Manifest in davidson_bench.rs:88](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-cli/examples/davidson_bench.rs#L88) does not parse the captured `potential_fnv`, `density_fnv`, convergence-profile flag, strict prefix or buffer tolerance.
- [Replay loading:749–787](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-cli/examples/davidson_bench.rs#L749) checks the FFT grid and band count. It reads the rank file's row hash and layout, but prints them without matching them to the newly constructed operator's rows.
- [ScfProblem::replay_solve](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-engine/src/scf.rs#L2247) injects `X,HX`; [gamma_real.rs](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-engine/src/gamma_real.rs#L1037) loads the kept products.
- [Reproduction check:928](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-cli/examples/davidson_bench.rs#L928) prints whether work/eigenvalues reproduced; a mismatch is not an enforceable replay identity gate.

For example, changing coordinates, pseudopotential, XC, or a same-sized stick assignment while retaining ranks/grid/bands can combine old products with new operator applications. A different strict/standard profile can also silently change the replay contract. Orthogonality of X does not detect this.

Scope: the measurement-only `--snapshot` path, not ordinary SCF restart correctness. It nevertheless matters because optimization decisions now depend on replay evidence.

### Proposed contract

Separate **exact replay** from an explicitly labeled **algorithm A/B on the same physical problem**.

- Record a versioned operator fingerprint: cell and positions, PP content hashes, XC/spin/occupation settings, cutoff and ordered reciprocal basis, density and effective-potential fingerprints, relevant normalization/layout convention.
- Before injecting products, compare each rank's ordered row mapping with the captured mapping. Validate all dimensions and metadata, not only column count.
- Exact replay restores or verifies tolerance policy, criterion, start state and kept-product status. Its identity check exits unsuccessfully on a mismatch.
- A/B mode may change solver/kernel settings, but must preserve the physical operator. If a permissible implementation change alters potential/product rounding, recompute current `HX` or explicitly validate its consistency; do not silently claim exact replay.
- A rank-local mismatch must become an agreed collective error **before** any rank enters a divergent solve/collective sequence. Preserve the new operator-failure agreement guarantees.
- Store enough metadata to diagnose a mismatch without relying on filenames or a whole-input-file hash alone.

### Acceptance checks

Reject snapshots with (a) moved atom, (b) different same-valence PP, (c) changed XC, (d) changed local row order/ownership, (e) modified density bytes, or (f) incompatible stopping profile, while ranks/grid/bands are unchanged. The original snapshot must reproduce its captured contract. Add an explicit invalidation/recompute test for a permitted rounding-changing kernel A/B.

A failed eigenvalue reproduction check is useful output, but detecting stale products must happen before the solver uses them.
