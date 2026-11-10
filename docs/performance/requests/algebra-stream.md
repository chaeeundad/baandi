# dev08 request: algebra-stream switches of the serial real Davidson (2026-11-11)

Status: implemented, opt-in, NOT yet measured on dev08. Adoption rule (AGENTS.md, user 2026-11-10): every switch that improves the dev08 median wall without breaking a correctness clause is adopted and stacked; there is no percentage threshold. Mac numbers below are observations on a shared Apple Silicon machine (load 4 to 11 from other sessions, pure Rust faer build, no MKL), never acceptance evidence. Ledger rows: ALGEBRA-STREAM-* (attempts-ledger.md, section "Algebra stream").

Trigger: dev08 stage profile (MKL, si32 strict, before the DGKS default): of the summed Davidson iteration time h_apply 38 percent, orthonormalize 33, rotate 21, eigen 6, reduced 4, residual 2; at T=8 eigen 16, rotate 14, orthonormalize 35. 58 percent of the iterations append one correction column and cost about 28 ms, a third of it H; the rest is streaming the basis (up to 408 x 11877, 39 MB). H apply belongs to another agent and is not touched here.

## Switches (all default off; each is bit-identical to the default path when off)

| Switch | Item | What it does | Exactness when ON |
|---|---|---|---|
| `BAANDI_RITZ_FUSED=on` | 2 and 4 | Ritz vectors X = V C and residuals R = HV C - X diag(theta) in one pass over V and HV: each 256-row panel of HV C lives only in that panel's GEMM tile and is combined at once with the panel of X; H X is no longer stored, the residual stage no longer re-reads X and H X. H X is formed late (one single-side rotation, the same kernel) only when a contraction restart needs it. Also used for the B1 narrow subset. | **bit-identical** (same per-panel product calls into the same tiles, same residual expression; unit test against the pair rotation at 1/2/4/8 threads, Davidson test `ritz_fused_is_bit_identical` with B1 and B5 on and off, Mac si16/si32 density SHA equal). MKL: same dgemm calls and tiles, expected bit-identical, dev08 confirms. |
| `BAANDI_NARROW_PROJ=on` | 1 | Projection of a correction block of at most `NARROW_WIDTH` = 4 columns against the basis: S = Q^T X as one fixed-order eight-lane dot per entry (each basis column streamed once, the block stays in cache, no tiles, no tile sum), then X -= Q S in one pass over fixed 1024-row blocks (ascending column order per entry, no FMA). Replaces the two panel GEMMs, the projected-block buffer and the separate subtraction. With the DGKS default most appends run one such projection. | rounding class; thread-count independent (fixed pieces) |
| `BAANDI_NARROW_REDUCED=on` | 4 | Border of the reduced matrix V^T H V for at most 4 new columns by the same streaming dot kernel instead of the reduction-panel GEMM and tile sum. | rounding class; thread-count independent |
| `BAANDI_B5=par` | 3 | B5 bordered update (as `BAANDI_B5=on`) with the secular roots, Loewner weights and eigenvector columns computed per root on the pool, and the m x m coefficient product in fixed 32-column groups on the pool (`multiply_column_groups`). Serial at T=1. | identical roots and vectors to `BAANDI_B5=on` bit for bit; product bit-identical on faer (test), a fixed different dgemm grouping on MKL (rounding class against `B5=on`); B5 itself is rounding class against off |

Unconditional, bit-identical: the MKL panel copy of `dgemm_panels` (every MKL rotation and Identity product, default path) now copies whole columns with `copy_from_slice` instead of an indexed element loop with bounds checks (same values; strided outputs keep the element loop). dev08 clause 0 below checks the default path.

Not changed: B1 (`BAANDI_B1`, existing, default off). Its own dev08 A/B (b1-b4-b5-davidson.md part B) is still pending; the Mac comparison below does not make it the default (see item 2).

Counters: `results.json` -> `convergence.work.algebra_stream` (`fused_iterations`, `deferred_hx`, `narrow_projections`, `narrow_reduced_borders`, `b5_parallel_updates`, and the three env modes). Also corrected: `convergence.work.ortho.mode` now reports `dgks:0.1 (default)` when `BAANDI_ORTHO` is unset (it said `two-pass` since the DGKS default).

## Analysis per item

1. Single-column projection. Classical Gram-Schmidt needs the whole overlap Q^T x before the update, so Q is read twice per projection pass; one read is impossible without a different orthogonalization. What the panel path wastes for narrow blocks is everything around the two reads: 47 reduction-panel tiles summed afterwards, a projected n x k buffer written and read back, a separate subtraction pass, and (pure build) nano-GEMM calls of shape 16 x 1 x 16. The narrow kernels remove these. Fusing the next iteration's reduced border into the projection is not possible: the border V^T H v_new needs H v_new, which exists only after the projection and the H application (chain reduced -> eigen -> rotate -> residual -> precondition -> project -> H -> reduced). The symmetric rewrite (HV)^T v_new = ((HV)^T x - A s)/l moves the stream from V to HV instead of removing it, and adds a 1/l cancellation (l down to 0.1 under DGKS), so it was not implemented.
2. Rotation. B1 removes rotation work (fewer columns); RITZ-FUSED removes memory traffic of the columns that remain (H X write, X and H X re-read). They stack (B1 narrow iterations use the fused pass too). B1 stays a switch: the Mac shows the count sketch positive on si32 strict but the dense sketch near neutral and the decision rule is the dev08 part B.
3. Eigen stage. dsyevr with a range was already measured slower than dsyevd (ledger MKL-EVD); a threaded dsyevd would make results depend on the thread count. B5 replaces about half of the eigensolves by an O(m^2) secular solve plus one m x m product; `par` puts both on the idle pool at T=8. Mac microbenchmark (`bordered_timing`, faer, order 408): full solve 14.3 ms, bordered serial 10.3 ms, par T=4 4.0 ms, par T=8 3.2 ms (order 306: 7.0 / 4.8 / 2.0 / 1.7).
4. Reduced and residual stages. The residual stage disappears into the fused rotation except the column norms (one read of R). The reduced border for one column is a dgemv shape: NARROW-REDUCED.

## Mac observations (pure Rust release, shared machine; observations only)

Kernel microbenchmark at the si32 shape (n = 11877, m = 408, 68 bands; `cargo test --release -p dft-numerics --test algebra_stream stream_timing -- --ignored --nocapture`, medians of 15 or 7):

| kernel | T=1 panel / new (ms) | T=8 panel / new (ms) |
|---|---|---|
| projection k=1 (overlap + update) | 3.64 / 1.60 | 0.99 / 0.52 |
| projection k=2 | 5.33 / 2.19 | 1.10 / 0.72 |
| projection k=4 | 5.15 / 4.01 | 1.37 / 1.11 |
| projection k=8 (not used: NARROW_WIDTH = 4) | 7.85 / 8.03 | 2.05 / 1.90 |
| reduced border k=1 | 1.11 / 0.56 | 0.40 / 0.24 |
| rotation + residual, 68 bands (pair vs fused) | 47.7 / 47.0 | 11.9 / 11.5 |

On faer the rotation is compute bound (about 28 GFLOP/s), so the fused pass saves only 1 to 3 percent there; on MKL the same products run about 3x faster and the saved 3 n b streams (about 20 MB per si32 iteration) are a larger share. dev08 decides.

In-solver observations, si32-gamma seed 601 (Mac, pure faer, release). Stage seconds divided by the `h_apply` stage of the same run (H apply is untouched by every switch; this removes most of the machine-load noise), medians of 3 runs (strict T=1) or 2 runs (others); "non-H" = summed iteration time minus h_apply, over h_apply. Davidson iterations and H vectors were identical in every arm (strict 362 / 4907, standard 45 / 2058), every run converged, energies within 6e-13 Ha of OFF, density SHA equal at T=1 and T=8 within every arm (si32 strict, si16 strict and standard), RITZ-FUSED SHA equal to OFF everywhere.

| arm | strict T=1: rotate+res / orth / reduced / eigen / non-H | strict T=8: same | standard T=1 non-H | standard T=8 non-H |
|---|---|---|---|---|
| OFF | 1.421 / 0.574 / 0.254 / 0.183 / 2.305 | 1.136 / 0.977 / 0.193 / 0.654 / 2.907 | 0.668 | 1.387 |
| RITZ_FUSED | 1.467 / 0.574 / 0.255 / 0.182 / 2.351 | 1.204 / 0.963 / 0.196 / 0.638 / 2.945 | 0.673 | 1.379 |
| NARROW_PROJ | 1.448 / 0.536 / 0.256 / 0.177 / 2.286 | 1.158 / 0.925 / 0.201 / 0.628 / 2.857 | 0.667 | 1.384 |
| NARROW_REDUCED | 1.425 / 0.586 / 0.241 / 0.180 / 2.300 | 1.150 / 0.965 / 0.191 / 0.629 / 2.882 | 0.657 | 1.372 |
| B5=on | 1.424 / 0.594 / 0.252 / 0.157 / 2.317 | 1.135 / 0.982 / 0.195 / 0.589 / 2.846 | 0.657 | 1.357 |
| B5=par | 1.439 / 0.579 / 0.254 / 0.163 / 2.320 | 1.143 / 0.978 / 0.193 / 0.492 / 2.748 | 0.670 | 1.335 |
| B1=on, count sketch | 1.002 / 0.576 / 0.253 / 0.177 / 1.904 | 0.799 / 0.990 / 0.196 / 0.660 / 2.680 | 0.468 | 1.187 |
| B1 + RITZ_FUSED | 1.016 / 0.566 / 0.253 / 0.177 / 1.910 | - | - | - |
| STACK (fused, narrow x2, B5=par) | 1.493 / 0.541 / 0.243 / 0.162 / 2.316 | 1.209 / 0.947 / 0.186 / 0.499 / 2.789 | 0.661 | 1.336 |
| STACK + B1 count | - | 0.857 / 0.954 / 0.185 / 0.503 / 2.535 | 0.463 | 1.149 |

Reading (Mac, faer): B1 with the count sketch is the only large effect (-17 percent non-H at strict T=1, -8 at T=8, -30 / -14 percent standard; 245 of 362 strict iterations narrow, 15 of 45 standard). B5=par cuts the eigen stage by 25 percent at T=8 (B5=on alone 10 percent). NARROW_PROJ cuts the orthonormalize stage by 5 to 7 percent and NARROW_REDUCED the reduced stage by 1 to 5 percent: in the solver the narrow kernels act only on the projection GEMMs, while the stage also holds normalisation and the scalar in-block MGS. RITZ_FUSED is neutral to slightly negative here (+2 to +6 percent on rotate+residual, inside the run-to-run spread of this machine): the faer rotation is compute bound, so removing three n x bands streams buys nothing measurable; it is aimed at the MKL rotation, which is 3x faster and therefore relatively more memory bound. Every number above is an observation; the dev08 A/B decides.

B1 default: not switched on. The Mac evidence is consistent with the earlier B1 row (count sketch positive, dense sketch near neutral), but the user's rule asks for a dev08 median-wall gain and B1's own pre-registered part B is still pending; it stays `BAANDI_B1=on BAANDI_B1_SKETCH=count` (arm 5 below).

## dev08 protocol

One MKL binary (`MKLROOT=/opt/intel/oneapi/mkl/2026.1`, `MKL_NUM_THREADS=1`) and one pure build at the branch head; environment-only A/B; node idle (load below 0.30). Cells `benchmarks/cases/si16-gamma` and `si32-gamma`, profiles strict and standard, perturbed seeds 601-610 (+-1e-7 Bohr as MPI-M1-GATE3), T=1 (taskset core 0) and T=8 (cores 0-7). Interleave OFF/ON per seed and compare paired per seed (sign test and median of the per-seed differences).

Arms (OFF = all switches below unset, every current default unchanged):
1. `BAANDI_RITZ_FUSED=on` alone (MKL and pure).
2. `BAANDI_NARROW_PROJ=on` alone.
3. `BAANDI_NARROW_REDUCED=on` alone.
4. `BAANDI_B5=on` and `BAANDI_B5=par`, each against OFF (this also answers b1-b4-b5-davidson.md part C if not yet run); one extra run per cell with `BAANDI_DAVIDSON_STAGES=on` for the `eigen` stage.
5. `BAANDI_B1=on BAANDI_B1_SKETCH=count` alone and together with `BAANDI_RITZ_FUSED=on` (part B of b1-b4-b5-davidson.md, if not yet run).
6. STACK = `BAANDI_RITZ_FUSED=on BAANDI_NARROW_PROJ=on BAANDI_NARROW_REDUCED=on BAANDI_B5=par` on all cells and seeds; STACK+B1 (count) as a second stack arm.
7. If an arm alone is slower in a cell, STACK without it.

Stage profile: one `BAANDI_DAVIDSON_STAGES=on` run per arm and cell (seed 601) for the mechanism (`rotate`, `residual`, `orthonormalize`, `reduced`, `eigen`).

## Correctness clauses (decisive; unchanged by the adoption rule)

0. Default path (all switches unset) of the branch binary: energy bits, density SHA and counts equal to a build of the parent commit bc33efe on every cell (MKL and pure); this covers the MKL panel-copy rewrite.
1. RITZ-FUSED (bit-identical class): energy bits, density SHA, SCF/Davidson/H counts and certificate verdicts equal to OFF in every run, T=1 and T=8, MKL and pure. Any difference is a failure of the switch (report the first differing run).
2. B5=par against B5=on: Davidson and H counts and certificate verdicts equal; pure build energy bits and density SHA equal; MKL rounding class.
3. Rounding-class arms (NARROW_PROJ, NARROW_REDUCED, B5, B5=par, STACK): every run converges and certifies; |E_on - E_off| within the class (strict 1e-10 Ha, standard 5e-9 Ha); density SHA reproducible run to run and equal between T=1 and T=8 within each arm (thread-count independence); no contraction may fail its orthogonality gate more often than OFF (report `restarts` and `contracted_restarts` if available).
4. Report the Davidson and H counts per seed (a rounding change may move a trajectory; not a failure by itself).

## What to send back

Per run: wall, phase timers, `convergence.work` (Davidson and H vectors recorded plus certificate, `algebra_stream`, `b1`, `b5`, `ortho`, `davidson_stages` where on), energy (17 digits), density SHA, exit code; binary sha256, commit, build command, uptime before and after, interleave order.
