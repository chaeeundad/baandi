# Work order: BAANDI at least 10 percent faster than QE

Audience: orchestrator, `perf-optimizer`, implementers. Conventions: AGENTS.md (adoption rule, comparison rules of 2026-11-14) and `docs/performance/AGENT-BRIEF.md`.

## 0. Provenance and limits

- Written from a source read of BAANDI `40b6b84` and QE 7.6 (`KS_Solvers/Davidson/regterg.f90`, `PW/src/c_bands.f90`, `electrons.f90`, `h_psi.f90`, `vloc_psi_acc.f90`, `sum_band.f90`, `mix_rho.f90`, `LAXlib/rdiaghg.f90`, `FFTXlib`). **Nothing was built or run.**
- Labels: **M** dev08 number from the ledger or `qe-stage-decomposition-2026-11-14.md`; **D** arithmetic on M; **E** estimate; **H** hypothesis.
- Earlier estimates of this author were wrong more than once (Track A, B1, B5). Every item therefore starts with a count or a timer, and has a kill rule.

## 1. Goal, decided direction, and what "faster" means

**User decisions of 2026-11-14:**
1. One calculation path. The standard path (QE's stopping rule and eigenvalue criterion) becomes the only SCF path; accuracy is a threshold, as in QE. The fixed-point certificate becomes an optional verification step after convergence. The present strict path stays as a test reference tool.
2. Target: BAANDI wall at most 0.90 times QE wall.

**User decisions of 2026-11-14 (second set):**
3. The goal covers all three parallel models, each compared only with the same model: one core against one core, BAANDI threads against QE OpenMP at the same thread count, BAANDI MPI against QE MPI at the same rank count. No analysis, table, ratio or argument may put one model against another (no "threads against MPI", not even as information).
4. Every QE against BAANDI comparison uses identical conditions: same cell and coordinates, cutoff, FFT grid, pseudopotential file, band count, start, threshold, cores and binding.
5. The fixed-point verification (certificate) is ON by default in production runs. The timed BAANDI runs therefore include it; its cost is inside the target.
6. "QE" is the fastest QE configuration found in P0.3 for that cell and model.
7. GPU stays out of scope for now.

**Acceptance:**
- Cells si16, si32, si64 Gamma plus one k-point cell and one non-silicon cell; displacement f1e-3 (the headline amplitude) and unperturbed.
- Thresholds: QE `conv_thr` 1e-8 and 1e-10, BAANDI the same estimator at the matching threshold, verification on.
- Three separate pass/fail results, one per model (decision 3): BAANDI wall at most 0.90 times QE wall in each.
- Report the verification time as its own line in every table, so the price of the default is visible; do not subtract it from the headline.

## 2. The budget (si64, f1e-3, one core, conv_thr 1e-8)

QE total 26.55 s (M). Its split at 1e-8 is not measured; the 1e-10 split scaled to 26.55 s gives (E): `h_psi` 14.8 s (local FFT 11.3, nonlocal 3.4), Davidson algebra 7.3, `sum_band` 2.0, init 0.9, potential and mixing 0.5, rest 1.0.

BAANDI today 45.98 s (M): fine diagonalization 29.43, ladder coarse 10.43, density 2.44, energy 0.54, potential 0.42, setup 0.88, mixing 0.18, other 1.1.

Target 23.9 s. What must be true (D/E):

| Component | QE (E) | BAANDI target | How |
|---|---|---|---|
| H application | 14.8 | 11.5 to 12.5 | same column count as QE; cheaper column (C1, C2) |
| Davidson algebra | 7.3 | 6.0 to 7.0 | QE's operation count (B1 to B3), cheaper reduced solve (B4) |
| Density | 2.0 | 1.6 | pruned transforms (C1), shared with C2 |
| Setup, potential, mixing, energy | 1.4 | 1.3 | energy without operator work |
| SCF work in step equivalents | 11 | 9.5 to 10.5 | ladder with an efficient hand-off (D) |
| **Total** | **26.55** | **21.5 to 23.9** | |

Reading: there is no single lever. Parity needs B (algebra) and D (step count); the 10 percent margin has to come from C (cost per H column below QE's) and from the ladder, which QE does not have. If C1 and C2 both fail, the goal is not reachable at one core without a new idea; say so early (section 8).

## 3. Why BAANDI is 1.74 times QE today (standard path, summary of the source comparison)

- **Algebra per Davidson iteration.** QE (`regterg.f90`): residuals for the unconverged roots only (two GEMM, N x nbase x notcnv), reduced border (two GEMM), no orthogonalization (generalized problem, `DSYGVX` for the lowest nvec), basis at most 2 nbnd, Ritz vectors formed once at exit (`regterg:last`). BAANDI (`real_davidson.rs`): in every wide iteration Ritz vectors and H images of **all** bands (two GEMM, N x m x bands), projection of the correction block on the whole basis with one or two passes, basis up to 6 bands, full `dsyevd` of order m. Measured stage ratios (si64): rotate 12x QE's `last`, orthonormalize 7.5x QE's `overlap`, eigen 2.9x.
- **SCF work.** Fine stage 10 steps plus a coarse stage worth 3.5 to 4 fine steps, against QE's 11 steps (D).
- **Cost per H column** is about equal to QE's (within 15 percent, D; BAANDI has no FFT timer).

---

## Phase 0. Definitions and instruments (first; no speed change expected)

### P0.1 One path
- `ScfOptions.standard` becomes the default path for every run; the threshold is the accuracy control (`StandardProfile.threshold_ha`, `convergence.rs`; floor `FLOOR_HA` 5e-14).
- Certificate: `verify = "fixed_point"` (name free) runs the existing certificate once on the converged state and reports it; it no longer drives per-step tolerances. **Default: on** (user decision 5); `verify = "off"` exists for users who do not want it.
- Verification cost is now part of the target. Today one certificate solve costs about 0.18 s (si32) and 1.24 s (si64) at one core (M, P48A-OTHER), about 3 to 5 percent of QE's wall. Reduce it: (a) keep CERT_SCREEN and the reused potential; (b) the verification solve starts from converged orbitals, so with C2 its initial apply costs half a transform per pair; (c) define what the verification must establish on the single path (the fixed-point residual of the density at the delivered threshold, not eigen-residuals at 1e-9) and write that definition into the ADR before measuring. If a failed verification is found, the run continues SCF steps and verifies again, as today.
- The strict session stays callable for tests and references; all accuracy gates keep using it.
- Check the standard path at the 1e-10 threshold (never measured): convergence, delivered error against the strict reference, no stall, on the generality set.
- Check forces: maximum force difference standard (1e-8 and 1e-10) against strict on si16/si32 displaced and mgo8. Report; the user decides which threshold relaxations and MD use.

### P0.2 Timers that close
- FFT timer split by kind (pair inverse, pair forward, density) and a nonlocal timer inside `h_apply`, always available under `BAANDI_DAVIDSON_STAGES=on`; `qe-stage-decomposition` says this is missing.
- Stage times reported per stage (ladder coarse, fine, verification), not pooled.
- QE side: store all repetitions; extract `h_psi` columns from the `fftw` call count (wave transforms only: `dffts%wave_clock_label = 'fftw'`, `PW/src/data_structure.f90:79`) minus `sum_band` transforms, per run.

### P0.3 QE reference that cannot be questioned
- Add QE runs **without** `diago_full_acc` (QE's default is `.false.`; `make_qe_input.py:146` sets it true, which makes QE converge the empty bands at the full threshold) and report the faster of the two as "QE".
- Add one arm each with `diagonalization = 'rmm-davidson'` and `'ppcg'` on si32/si64. If one is faster than `david`, the target is that number. State this in the table.
- QE OpenMP at 1, 2, 4, 8 threads and QE MPI at 1, 2, 4, 8 ranks for every cell, same binding as the BAANDI runs. The fastest QE configuration per cell and model is "QE" (decision 6).
- Remove or relabel existing mixed-model columns: `qe-vs-baandi-2026-11-13.md` still carries "T=8 against np=8" ratios; they must not appear in any new table or ledger conclusion (decision 3).

---

## Phase B. Davidson algebra at QE's operation count

All items act on `crates/dft-solvers/src/real_davidson.rs` (and the complex counterpart afterwards), behind one switch `BAANDI_DAVIDSON=lean` (default off until dev08), Ritz-value criterion only (`ritz = true`).

### B1. Never rotate all bands inside the loop
- Today: `narrow` requires `2 * set.len() < bands`; otherwise `ritz_and_residuals` runs for all bands.
- Change: in lean mode always take the selected-root branch (residuals for the unconverged set only, any width). Ritz vectors of all bands are formed once, at return (X only; H X only when the caller needs it, see C2) and at a contraction restart.
- This is exactly QE's `update` / `last` split. Work count to report: GEMM 2mnk by site before and after.
- Ledger check: B1-PART-B (adopted for narrow iterations) is extended to all widths in the Ritz criterion, where no residual of an unselected root is needed; PRT1's stale-flag failure concerned the residual criterion.

### B2. Basis size and restart
- With B1 the cost of a large basis is the GEMM length only. Sweep capacity 2, 3, 4, 6 per band **with B1 on** (`davidson_capacity_per_band` exists), counts first (H columns, iterations, GEMM flops), then dev08.
- Ledger check: A13 and P31A rejected smaller capacities under "all roots always rotated and orthonormalised" and name "a non-orthonormal basis or locking" / "independent dev08 evidence" as the reopening condition. B1 and B3 change that premise.

### B3. Orthogonalization: two candidates, measured against each other
- **B3a, generalized problem (QE's way).** Append unit-norm corrections without projection; keep S = V^T V beside H_r (one more border GEMM, which replaces the projection GEMMs); solve the generalized reduced problem. Fallback: if the Cholesky factor of S fails or its smallest pivot ratio is below a fixed bound, orthonormalize the basis once and continue.
  - Ledger check: P37A rejected this by a rule that covered **all** profiles; the failures were at strict residual levels (Gram eigenvalues near 1e-16), and P37A-3 measured 45 to 50 percent lower multiply-add counts per step for the standard profile. With one standard path the all-profile rule no longer applies; re-register on the standard path at thresholds 1e-8 and 1e-10 with the same safety metric (fallback rate below 5 percent, no lost root).
- **B3b, orthonormal basis, one pass.** Keep the orthonormal basis but project once (DGKS already skips the second pass in most cases) and replace the in-block MGS by Cholesky QR of the small block Gram with a fallback (ledger A52 rejected CholeskyQR2 on a skew cell: keep the MGS fallback on any conditioning doubt).
- Decide by dev08 wall at equal H counts. Kill rule: drop B3a if its fallback rate exceeds 5 percent on any generality cell.

### B4. Reduced eigenproblem
- QE spends 2.0 of 30.4 s (si64) in `DSYGVX` (bisection and inverse iteration). With a basis of at most 2 to 3 bands, the full divide-and-conquer solve is small; for B3a use Cholesky reduction plus `dsyevd`. Keep B5 (bordered update) only if the stage timer still shows more than 5 percent.
- Ledger check: MKL-EVD (range solve slower than `dsyevd`) is respected.

### B5. Expected and gate
- Expected (E): BAANDI algebra from about 2.4 times QE's to 0.9 to 1.0 times; fine diagonalization si64 29.4 to about 21 s.
- Gates: same SCF count and energy within the class against the current path on the generality set, perturbed-start distributions, density identical across thread counts.
- Kill rule: if total H columns rise by more than 5 percent median on si32 and two non-silicon cells, the lean mode does not replace the default.

---

## Phase C. H column cheaper than QE's

### C1. Pruned transforms on MKL (measure first; the code exists)
- `BAANDI_PAIR_FFT=pruned` (`PairFusedPlan` on MKL 1D batch descriptors) was never measured on dev08. QE prunes (`cft_1z` on wave sticks, `cft_2xy` with `pl2ix`); the BAANDI MKL default is the full-grid adapter (`GridFusedPlan`, `fused.rs:499`). PFFT-5: per transform equivalent at 36^3, full MKL adapter 195 us (M), QE 135 us (M), MKL stick pipeline 97 us (E).
- Run the dev08 arm at T=1 and T=8 with the P0.2 FFT timer. If it wins, make the density transforms use the same plan.
- Kill rule: if the pruned plan is not at least 10 percent faster per transform than the full-grid plan at si32 and si64 grids, record and stop FFT layout work.

### C2. Do not transform the same orbitals twice (new)
- Observation: every SCF step transforms the converged orbitals to real space for the density, and the next Davidson solve starts by applying H to the same orbitals, which transforms them again. QE does the same (`regterg` starts with `h_psi` on all nvec; 15 x 132 of about 6100 columns, 32 percent of its H columns at si64, D).
- Change, two parts:
  1. At return, form H X by rotation of the stored H images (`rotated_products`, one GEMM; `hx_formed` already tracks it). Then H_new X = H_old X + (V_new - V_old) X: kinetic and nonlocal parts are unchanged, only the local potential difference acts.
  2. Keep psi(r) of the returned orbitals from the density build (pair layout). The initial "apply" of the next solve is then: multiply the cached psi(r) by (V_new - V_old), one forward transform per pair, add to the rotated H X. No inverse transform, no nonlocal products, no kinetic pass.
- Saving per SCF step (D on QE's si64 per-column split: local 2.15 ms, nonlocal 0.64 ms): 132 columns x (half the local part + nonlocal) is about 0.23 s of 0.37 s, less one GEMM of about 0.05 s; about 2 s over 11 steps, 7 percent of QE's wall.
- Requirements: memory gate (pair-complex grid per two bands; si64 about 0.8 GB, E), automatic fallback to the plain apply above the gate; exactness class: rounding (H X is a sum of two parts); the accumulated error does not grow because every solve rebuilds its basis from X.
- Also valid for the ladder coarse stage and for the first step after a relaxation extrapolation (then the nonlocal part changes: fall back).
- Stage 0: count the transforms and products saved on si32/si64; unit test H_new X against a plain apply (relative 1e-13).
- Ledger check: no row (text search "real-space cache", "psi(r)", "potential difference": none).

### C3. Small exact items
- Energy without operator work in ordinary steps: ENERGY-IDENTITY-DEV08 was neutral to slightly slower; do not repeat. With C2 the nonlocal expectation values come from the cached projections at no cost; re-evaluate only then.
- Remaining finite-value scans and copies inside the apply path: check with the P0.2 timers; act only above 2 percent.

### C4. Build and memory (cheap, measure once)
- Profile-guided optimization of the release binary (pure and MKL builds); `Cargo.toml` has no `[profile.release]` section (FIN1 measured `codegen-units = 1`): fix the profile in the repository.
- Transparent huge pages for the basis, H images and FFT buffers (`madvise` or a huge-page allocation for the few large arrays): streaming GEMM and stick scatter are TLB-heavy. Measure with the system setting first, then decide on code.
- Ledger check: no rows for PGO or huge pages.

---

## Phase D. SCF work below QE's step count

### D0. Measure the baseline that is missing
- BAANDI fine SCF steps **without** the ladder at f1e-3 and f1e-2 against QE's (QE 11 and 13 at si32, M). Per step: estimator, true error against the strict reference.
- Equal-mixer experiment: an arm that reproduces QE's mixer (G space, Hartree metric as `rho_ddot`: weight 1/G^2 on the smooth-grid vectors, beta 0.7, history 8, modified Broyden of `mix_rho.f90`; BAANDI has `hartree_mixing` and Pulay, which spans the same subspace) at the same thresholds. If the step count then equals QE's, the difference is the mixer; if not, it is the eigensolver accuracy per step.
- Ledger check: P40A-1 (mixer survey: case-dependent) was run before the atomic start, at 1e-10 relative density noise, with both profiles under one rule. This arm is a diagnostic at the headline amplitude on the single path, not a new default by itself.

### D1. Ladder hand-off
- Fact: unperturbed, the fine stage needs 5 steps (QE 7); at f1e-3 it needs 10 (QE 11) (M). The ladder's advantage disappears exactly where the headline is measured.
- Scan (counts on the Mac, then dev08): coarse stop level 1e-2, 1e-3, 1e-4 (today fixed 1e-4, `cutoff_ladder.rs`), coarse fraction 0.5 and 0.7, at f1e-3. Report coarse cost in fine-step equivalents and fine steps. Rule candidate: stop the coarse stage when its estimator falls below the measured coarse-to-fine density difference of the previous run of that cell class (a general rule, not a per-system constant, has to be shown on the generality set).
- First fine step: the standard ratchet starts from the cap; test starting it from the coarse stage's last tolerance (the orbitals are already converged to it).
- Ledger check: LADDER-S2-SECANT (rejected) is not repeated.

### D2. Decision
- If D0 shows the mixer is the cause, adopt the better mixer for the single path by the adoption rule. If D1 finds a hand-off with fine steps at or below QE's minus two at f1e-3, the ladder supplies the margin of section 2.
- Kill rule: if neither D0 nor D1 lowers total step equivalents below QE's count, the ladder is kept only where it is a measured net gain and the margin must come from Phase C alone.

---

## Phase E. Parallel models

Each model is its own target and its own table (decision 3). Nothing in this phase compares one model with another.

### E1. Threads: BAANDI threads against QE OpenMP, same thread count
- Known (M, si32, unperturbed, 1e-10): QE one thread 3.484 s, eight OpenMP threads 1.957 s. BAANDI numbers at the same thread counts on the single path do not exist yet: measure after P0.1.
- Work: per-stage time at 1 and 8 threads after Phase B (the wide algebra shrinks and the reduced solve is smaller); open arms `B5=par`, `HAPPLY_OVERLAP`, `ORTHO_BLOCK` (pure builds), `PHASE_THREADS` (on).
- Gate: 0.90 at 2, 4 and 8 threads on the acceptance cells.

### E2. MPI: BAANDI MPI against QE MPI, same rank count
- Known (M): QE at 8 ranks, f1e-3: si32 1.737 s (1e-8) and 2.138 s (1e-10); si64 6.449 s and 7.497 s. BAANDI MPI (M2d, `--features mpi`) has no measurement on the single path at these cells; the wiki freezes MPI at M2d "until multi-node".
- Step 1 (measurement): BAANDI MPI at 1, 2, 4, 8 ranks on si32 and si64, f1e-3, both thresholds, with the stage timers; QE MPI at the same ranks with its own timers. Produce the stage table rank for rank.
- Step 2 (after Phase B): the lean Davidson changes what has to be distributed. QE distributes plane-wave rows (and FFT sticks and planes) over ranks; all its Davidson GEMMs are local products followed by one reduction (`mp_sum` on `hr`, `sr`), and the reduced problem is solved on one rank. Bring the BAANDI MPI path to the same structure for the lean solver: row-distributed basis and H images, local GEMM plus one `Allreduce` per reduced border, serial reduced solve, distributed stick transforms. M2d already distributes rows (`solve_distributed_rows`); list what in the lean solver still runs replicated.
- Step 3: C2 in the distributed layout (cached psi(r) is per-rank slab data).
- Gate: 0.90 at 2, 4 and 8 ranks on the acceptance cells. This is expected to be the hardest of the three targets; it gets its own detailed work order once step 1 exists.

### E3. GPU
Out of scope for now (decision 7).

---

## 7. Order

| Order | Item | Kind |
|---|---|---|
| 1 | P0.1 to P0.3 | definitions, timers, honest QE reference |
| 2 | C1 dev08 arm, D0 measurements | measurement, no code |
| 3 | B1, then B2 sweep | implementation, counts first |
| 4 | B3a against B3b, B4 | implementation |
| 5 | C2 | implementation |
| 6 | D1 scan, D2 | measurement then adoption |
| 7 | C4 | measurement |
| 8 | E1 threads, E2 step 1 (MPI stage table) | measurement |
| 9 | E2 steps 2 and 3 | implementation, own work order |

## 8. Honest outlook

- Phase B plus D0/D1 should bring the one-core ratio from 1.74 to about 1.0 to 1.15 (E). That is parity, not the goal.
- The last 10 to 20 percent depends on C1 (pruned MKL transforms) and C2 (no double transforms), together about 12 to 18 percent of QE's wall if both work (E), plus the ladder.
- If P0.3 shows a faster QE configuration, the target moves by that amount.
- The default-on verification adds about 3 to 5 percent of QE's wall today (M); it must be reduced (P0.1) or the other items must cover it.
- The MPI target has no BAANDI baseline on the single path; its feasibility is unknown until E2 step 1.
- Checkpoints: after step 4, if the ratio is above 1.25 on si32 and si64, stop and re-plan; after step 6, if above 1.0, the goal needs a decision on scope (GPU, or a different eigensolver class).

## 9. Decisions

All four open decisions were answered on 2026-11-14 and are recorded in section 1 (decisions 3 to 7). Record them in AGENTS.md and `comparison-rules.md` before any new table is written. Open point for the ADR only: the exact statement the verification step must establish on the single path (P0.1 item c).

## 10. Not proposed (ledger respected)

Track A direct minimization (TRACK-A0), secant seeding (LADDER-S2-SECANT), hybrid RMM as default (A22), fp32 (P47B), band-batched or lane-staged FFT (A28, A11, A34), locking (P24A), preconditioner variants (A14, P47A), eigen schedule constants (A21, EIGEN-FACTOR-AB), Chebyshev filtering (needs about three times the H columns per step; H is two thirds of the wall).
