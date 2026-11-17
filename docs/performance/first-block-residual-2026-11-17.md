# First-correction-block residual spectrum of real SCF solves (BREAKTHROUGH-0-INSTRUMENT, review E3)

Review 2026-11-17 section 6 (R2) argued that the first Davidson iteration of every warm solve corrects
all B roots because the Ritz-change rule has no previous value, and asked for a measurement before any
implementation: at real SCF states, what is the singular spectrum of the occupied residual block
R_o = H X_o - X_o A_o of the kept start, does its effective rank stay near B or fall in late SCF steps,
and how many roots are really unconverged at the current Hamiltonian. `BAANDI_FIRST_BLOCK=on`
(default off, read-only, the solver is untouched; module `dft-engine/src/first_block.rs`) records that
at every warm solve. Labels: **O** Mac observation (values are exact arithmetic of the run, the machine
does not matter), **C** count. The rows below are copies of the run tables (si32 np=4 simulated,
si16 np=4, si8 np=4, and the threads si32 T=1 run); the full tables are in results.json
(`convergence.work.first_block`) of a CLI run and in the harness console.

## Columns
`it` SCF iteration of the map; `tol` the solve tolerance (Ha) of that map; `kept` whether the kept H X was
used (else H was applied to X here, B columns); `gap` theta[n_o] - theta[n_o - 1] of the small
Rayleigh-Ritz on the kept block; `sig_max`, `sig_min` the largest and smallest singular value of R_o
(n_o = 65 occupied of B = 68 for si32); `r_max_occ` / `r_max_emp` the largest true residual norm of an
occupied / empty Ritz pair; `r1e-k` the effective rank = number of singular values above 10^-k sig_max;
`r>tol` singular values above the solve tolerance; `unc_all` / `unc_occ` roots whose true residual exceeds
their threshold (all B / occupied); `rule` the roots the Ritz-change rule corrects in the first iteration
(all B, capped by the capacity rule: 68 at capacity 3B); `defect` max |X^T X - I|.

## si32 np=4 production + ladder (O, seed 601 f1e-3; identical to the threads run's structure)
| seq | stage | kind | it | tol | kept | gap | sig_max | sig_min | r1e-1 | r1e-2 | r1e-3 | r1e-4 | r>tol | unc_all | unc_occ | rule |
|---|---|---|---:|---:|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 0 | coarse | cold (atomic start) | 1 | 5.0e-3 | no | 3.9e-2 | 2.05e-1 | 1.44e-1 | 65 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 1 | coarse | ordinary | 2 | 5.0e-3 | yes | 2.1e-3 | 9.1e-2 | 6.6e-2 | 65 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 2 | coarse | ordinary | 3 | 1.5e-3 | yes | 4.7e-3 | 7.6e-2 | 3.3e-2 | 65 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 3 | coarse | ordinary | 4 | 1.1e-4 | yes | 4.8e-3 | 6.3e-2 | 3.2e-3 | 57 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 4 | coarse | ordinary | 5 | 1.1e-4 | yes | 3.0e-3 | 5.2e-2 | 2.9e-3 | 49 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 5 | coarse | ordinary | 6 | 1.1e-4 | yes | 1.6e-3 | 2.0e-2 | 1.5e-3 | 58 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 6 | coarse | ordinary | 7 | 4.5e-5 | yes | 1.3e-4 | 1.4e-2 | 7.2e-4 | 48 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 7 | coarse | ordinary | 8 | 5.7e-6 | yes | 2.7e-4 | 9.3e-3 | 2.5e-4 | 37 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 8 | coarse | ordinary | 9 | 2.6e-6 | yes | 3.5e-4 | 3.7e-3 | 1.8e-4 | 50 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 9 | fine | ordinary (ladder hand-off) | 1 | 7.6e-7 | no | 4.5e-4 | 1.28e-1 | 6.5e-2 | 65 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 10 | fine | ordinary | 2 | 7.6e-7 | yes | 3.7e-4 | 2.4e-3 | 8.2e-4 | 65 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 11 | fine | ordinary | 3 | 7.6e-7 | yes | 3.7e-4 | 6.6e-3 | 9.5e-4 | 65 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 12 | fine | ordinary | 4 | 7.6e-7 | yes | 2.7e-4 | 6.9e-3 | 4.2e-4 | 48 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 13 | fine | ordinary | 5 | 1.3e-8 | yes | 3.9e-4 | 4.1e-3 | 7.1e-5 | 14 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 14 | fine | ordinary | 6 | 4.9e-9 | yes | 4.2e-4 | 1.4e-3 | 1.7e-5 | 4 | 65 | 65 | 65 | 65 | 68 | 65 | 68 |
| 15 | fine | ordinary | 7 | 4.2e-10 | yes | 4.1e-4 | 1.2e-3 | 7.6e-6 | 1 | 56 | 65 | 65 | 65 | 68 | 65 | 68 |
| 16 | fine | ordinary | 8 | 3.5e-10 | yes | 4.1e-4 | 8.4e-4 | 5.0e-6 | 1 | 46 | 65 | 65 | 65 | 68 | 65 | 68 |
| 17 | fine | ordinary | 9 | 3.5e-10 | yes | 4.1e-4 | 5.9e-4 | 2.9e-6 | 1 | 42 | 65 | 65 | 65 | 68 | 65 | 68 |
| 18 | fine | ordinary | 10 | 9.5e-12 | yes | 4.1e-4 | 3.9e-4 | 1.2e-6 | 1 | 35 | 65 | 65 | 65 | 68 | 65 | 68 |
| 19 | verification | verify (failed) | 10 | 2.0e-12 | yes | 4.1e-4 | 3.1e-4 | 8.7e-7 | 1 | 32 | 65 | 65 | 65 | 68 | 65 | 68 |
| 20 | fine | ordinary | 11 | 2.0e-12 | yes | 4.1e-4 | 2.3e-4 | 6.1e-7 | 1 | 34 | 65 | 65 | 65 | 68 | 65 | 68 |
| 21 | verification | verify (passed) | 11 | 5.0e-14 | yes | 4.1e-4 | 2.4e-4 | 1.8e-7 | 1 | 1 | 60 | 65 | 65 | 68 | 65 | 68 |

Overlap defect of every kept block 1e-14 to 3e-14 (C). The threads si32 T=1 run (CLI, REUSE-PSI kept
start) gives the same picture with its own trajectory: fine it 5 to 10 r1e-1 = 39, 26, 2, 1, 1, 1 and
r1e-2 = 65, 65, 65, 65, 55, 37; verification r1e-2 = 39; r1e-3 = 65 everywhere; sig_min 2e-5 down to
3e-7 against tolerances 6e-9 down to 8e-14. si16 np=4 (B = 36, n_o = 33): r1e-1 falls to 1 from fine
it 6, r1e-2 to 2 at it 7 and 13 to 14 in the last maps, r1e-3 stays 33 (24 only in the last
verification); si8 np=4 (B = 20, n_o = 17): r1e-1 = 1 from fine it 6, r1e-2 between 1 and 17 (1 at
it 7, 9, 10 and both verifications), r1e-3 = 17 except 16 and 13 at two solves.

## Findings (what the numbers say; no implementation decision is taken here)
1. **Late SCF steps: the occupied residual block is dominated by ONE direction, but its tail is far above
   the tolerance.** From fine iteration 6 or 7 on (si8, si16, si32 alike) one singular value carries the
   block at the 10 % level (r1e-1 = 1); at the 1 % level 32 to 56 of 65 directions remain (si32), and at
   0.1 % all 65. The smallest singular value (3e-7 to 1e-5 Ha) is 4 to 6 orders of magnitude above the
   solve tolerance (1e-12 to 1e-9 Ha): every occupied direction is still "unconverged" by the
   solver's own threshold (`r>tol` = 65 = `unc_occ` at every solve). So a rank-r compression of the
   first block at a relative threshold would drop directions that the solve must still correct later; the
   saving is in the FIRST correction block only, and it must be checked at the level of the whole solve
   (H columns, iterations) rather than assumed from the spectrum.
2. **The Ritz-change rule and the true residual agree on the count.** `unc_all` = 68 = `rule` at every
   warm solve: at the start of every map all B roots really exceed their thresholds, because the
   potential change (V_new - V_old) X moves every Ritz pair by far more than the tolerance (sig_max
   2e-4 to 7e-3 Ha in fine maps against tolerances below 1e-6). The review's concern that the rule
   over-corrects at the first block is not confirmed on this trajectory: it corrects what the current-H
   residual also flags. The place where the rule could be loose is the LAST iterations of a solve, not
   the first block; that is a different measurement (per-iteration, inside the solver) and is not made
   here.
3. **The ladder hand-off is the worst start of the run.** The first fine map (prolonged coarse orbitals,
   `kept = no`) has sig_max 1.28e-1 and sig_min 6.5e-2 with full rank at every threshold: the prolonged
   block is 2 orders of magnitude worse (in residual) than the kept start of the next map (2.4e-3), and
   its residual block has no low-rank structure at all. Any R2-style compression has nothing to gain on
   it; a better hand-off (R5) would, and the number to beat is sig_max 1.3e-1 against the 2e-3 of a
   kept start one map later.
4. **Early coarse and early fine maps have full rank at every threshold** (seq 0 to 2, 9 to 11): the
   compression window opens only once the SCF residual has fallen by about three orders (fine it 5 for
   si32), which is also where the solves become short (4 to 7 Davidson iterations, 2 to 3 correction
   blocks). The share of the run wall in those late maps (profile-2026-11-17-phases.md: fine maps 5 to
   11 plus the two verification maps) bounds what first-block compression can give.
5. **Verification maps start as far from convergence as ordinary maps** (sig_max 2.4e-4 to 3.1e-4
   against tolerances 2e-12 and 5e-14): the verification solve is a full warm solve, not a check; the
   failed attempt (seq 19) and the passed one (seq 21) have the same spectrum shape.
6. **The boundary gap is stable** (3.4e-4 to 4.5e-4 Ha for si32 from the hand-off on; 3.9e-2 only at the
   atomic start) and the occupied / empty boundary never moves: the occupied block is well defined in
   every map, which is the precondition the review states for working with X_o separately.
7. The same-H redo (R1-a) did not occur in any of these runs (no `fine_redo` map); the ledger row
   BREAKTHROUGH-1-RESTART records the same observation.

Summary for the R2 decision: the effective rank of the first correction block does fall in late SCF
steps, but only at a 10 % relative threshold (to rank 1) and partly at 1 % (about half of n_o at si32);
at 0.1 % it stays at n_o, and in absolute terms every direction is 4 to 6 orders above the tolerance.
A first-block compression therefore cannot be judged from the spectrum alone; the test is a replayed
late-step solve (BAANDI_SNAPSHOT_DIR + `davidson_bench --snapshot`) with and without the compression,
counting H columns and iterations to the same converged result.

Commands (O): `env M2_PROFILE_MODE=production M2_PROFILE_LADDER=on M2_PROFILE_SERIAL=off
M2_PROFILE_RANKS=4 M2_PROFILE_CASE=si32 BAANDI_FIRST_BLOCK=on target/release/deps/m2_entry_cost-<hash>
--ignored --nocapture m2_entry_cost_profile` (np=4 simulated) and `BAANDI_FIRST_BLOCK=on
target/release/dft-cli run si32.toml` (threads). The diagnostic adds B H columns at the two `kept = no`
solves of each run (visible in the operator counters, never in the solver's counts) and about five
nl x B x B GEMMs per solve; results are bit-identical with it on and off (profile-2026-11-17-phases.md,
identity section).
