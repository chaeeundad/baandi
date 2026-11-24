# BAANDI work order — reduce full-cost SCF work

**Target** The reference SCF driver on dev08, real Γ, LDA, fixed occupations, MPI8. Screen two independent changes against the same verified baseline. Keep both changes opt-in during this work.

**Test budget** One focused local test per candidate. One shared release build on dev08. At most three complete Si64 SCFs, with a 120-second limit per run (six minutes of execution in total, excluding build and queue time). A failed local test removes that arm from the remote screen.

Reviewed against [`main@09f6628b`](https://github.com/chaeeundad/baandi/commit/09f6628be2bfb1d788166b4dd3540088fdb5d79f). If implementation has advanced, use the existing switch or patch and run only the missing check.

## H · Reuse the starting Hamiltonian products

**Why first.** The reference solver returns both `X` and `HX`, but the next ordinary map retains only `X` and applies the full Hamiltonian to every starting column. On the older MPI path, reuse removed 24.8–28.3% of full-H columns and reduced delivered wall by 1.2–5.1%. [Reference map](https://github.com/chaeeundad/baandi/blob/09f6628be2bfb1d788166b4dd3540088fdb5d79f/crates/dft-engine/src/refscf/map.rs#L150-L218) · [Davidson entry](https://github.com/chaeeundad/baandi/blob/09f6628be2bfb1d788166b4dd3540088fdb5d79f/crates/dft-solvers/src/ref_davidson.rs#L185-L205) · [Earlier dev08 result](https://github.com/chaeeundad/baandi/blob/09f6628be2bfb1d788166b4dd3540088fdb5d79f/docs/performance/attempts-ledger.md#L1046)

1. Add an opt-in reference-driver switch (proposed name `BAANDI_REFSCF_KEPT_HX`). Retain the prior Ritz block, its H products and the previous owned potential. Use the existing `band_density_sum_keeping` and `apply_local_change` with bounded slab storage. On a valid warm start, supply `H_old X + (V_new − V_old) X` to Davidson in place of its initial full-H application. Keep the overlap, rank and convergence checks.
2. Tie reuse to the exact retained orbital block, fixed operator, basis and row layout. The current map passes a constant potential epoch of zero. Track the actual old/new potentials for the update. Use a fresh H on the first map, cutoff/layout changes, retries and independent verification. All ranks must agree before entering the reuse branch.
3. **One local test.** In a tiny two-rank synthetic nonlocal fixture, change the local potential, then compare the reused starting product with a direct application. Require relative Frobenius error ≤ `1e-13` on a nonzero reference block and exactly `B` fewer full-H starting columns for one `B`-band hit. In the same test, invalidate the retained block and assert that the next entry takes the full-H path. [Existing delta-potential kernel](https://github.com/chaeeundad/baandi/blob/09f6628be2bfb1d788166b4dd3540088fdb5d79f/crates/dft-physics/src/distributed_real_gamma.rs#L1367-L1445)

## L · Use the existing cutoff-ladder ticket on the reference driver

**Why second.** The strongest prior result, 28–36% lower wall, was for displaced Si32 on the older strict single-core path. The new reference driver already reaches the Si64 estimator stop in 11 maps. This screen tests whether the saved fine-grid work pays for the coarse solve and handoff. Continue the existing ticket. [P3-LADDER](https://github.com/chaeeundad/baandi/blob/09f6628be2bfb1d788166b4dd3540088fdb5d79f/tickets/P3-LADDER.md) · [Prior measurements](https://github.com/chaeeundad/baandi/blob/09f6628be2bfb1d788166b4dd3540088fdb5d79f/docs/performance/attempts-ledger.md#L875)

1. Use `BAANDI_REFSCF_LADDER`, coarse cutoff `E₁ = 0.5 E₂`, a smaller compatible grid and one documented coarse stopping rule. Reuse the existing coarse rule if using the same density-residual norm and units. Its threshold of `1e-4 Bohr^(-3/2)` cannot be copied numerically into the reference driver’s `dr2` in Ry.
2. Inject orbitals by Miller index and Fourier-pad the density. Start the fine solve with empty Broyden history. Use its normal tolerances, audit and independent verification. Record coarse, handoff and fine time separately. Keep the fraction and stopping rule fixed throughout this screen.
3. **One local test.** On one tiny skewed Γ cell with two simulated ranks, check the copied coefficients and zero-filled new modes bitwise before any fine-stage reorthogonalization. Check Γ conjugate pairs, charge error ≤ `1e-12 e` and maximum change of an orbital Gram entry ≤ `1e-12`. Continue the same fixture to verified fine-grid delivery and compare its energy with the direct fine solve (≤ `1e-10 Ha` when both use tight settings). [Existing prolongation and coarse rule](https://github.com/chaeeundad/baandi/blob/09f6628be2bfb1d788166b4dd3540088fdb5d79f/crates/dft-engine/src/cutoff_ladder.rs)

## Run once and decide

Put the two focused tests in the existing `dft-engine` integration target `refscf`. Give them the prefix `fast_order_` and run `cargo test -p dft-engine --test refscf fast_order_`. Check the reported test count against the number implemented. Reuse the existing build cache.

If both local tests fail, skip dev08. Otherwise use the archived production Si64 MPI8 reference-driver input and matching strict energy from the verification campaign. Use one binary for all arms, eight bound physical cores on socket 0, one thread per rank, and identical input, UPF, seed, lightweight diagnostics, output and fixed-point verification. Capture external wall from MPI launch through exit. Read verification status offline with `dft-cli refscf-check --stop-to-verified <results.json>`. [Manifest and launch recipe](https://github.com/chaeeundad/baandi/blob/09f6628be2bfb1d788166b4dd3540088fdb5d79f/docs/performance/lanes/RSCF/stop-to-verified-dev08-request.md#L40-L76)

| Run | Switches | Evidence to retain |
|---|---|---|
| D · baseline | H off, L off | Verified delivery, wall, energy, fine maps, fine full-H columns, verification attempts and peak rank RSS |
| H · warm start | H on, L off | Same fields, plus cache hits, full-H starting columns replaced, delta-V columns and retained slab bytes/rank |
| L · ladder | H off, L on | Same baseline fields, plus coarse maps, coarse wall and handoff wall |

**Decision.** Stop the screen if D fails. Stop an individual arm on timeout, failed verification or standard-class energy error above `5e-9 Ha` against the matching strict result. If that strict result is unavailable, report the D-relative difference and mark strict-reference agreement *unverified*. For H, require a cache hit and fewer fine full-H columns. For L, require fewer fine maps and check that fine-stage time saved covers coarse and handoff cost. A verified arm with its expected work reduction and lower delivered wall is a **follow-up candidate**. If wall is within known noise, report **work reduced, speed unresolved**. A zero hit count or a missing coarse stage is a wiring failure.

Return one row per executed arm, including failures, with the commit, binary and input hashes, and one decision per arm. This development screen ends there. A combined arm, QE rerun, large-cell matrix and statistical timing campaign belong to subsequent validation of a surviving candidate.
