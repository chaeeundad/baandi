# Exclusive phase profile of the production MPI path (BREAKTHROUGH-0-INSTRUMENT, 2026-11-17)

Review 2026-11-17 section 13.1 asked for a profile whose top-level buckets are mutually exclusive, whose
inner rows nest inside them, and whose numbers are the maximum over ranks (critical path) with the
imbalance beside them, never a sum of CPU seconds. `BAANDI_PHASE_PROFILE=on` (default off) now gives
that; this note is its first use. Labels: **O** Mac observation (pure build, simulated ranks = threads of
one process on a shared 12-CPU Mac at load 7 to 12, so a tree sum is two thread barriers and an alltoall
a memcpy; absolute seconds are not dev08 numbers, the shares and the counts are what to read), **C**
exact count. dev08 (label D) numbers come through `docs/performance/requests/breakthrough-instrument.md`.
MPI against MPI only: no threads number is compared with a rank number here (the threads table at the
end is a separate observation of the same switch on the threads path).

## What the profile is
- Top level (exclusive): `init` (problem construction + the atomic start, which the first map's potential
  phase would otherwise carry), `coarse` (ladder coarse maps and their mixing), `fine` (ordinary
  full-cutoff maps and mixing), `fine_redo` (discarded first-step solves), `verify_pass` / `verify_fail`
  (verification maps by outcome), `materialize` (lazy gathers, output writing when the caller reports it).
  The run wall minus the sum is the remainder (prolongation, stop-rule estimate, harness overhead).
- Inside a map bucket: `potential`, `solve`, `density`, `energy`, `mixing` (the SCF phase clocks).
- Inside `solve`, exclusive stages: `start`, `apply` (full H columns and delta-V columns, with the pair
  transform counts), `reduced_gram` (overlap / reduced border), `reduced_eigen` (factor / EVD), `rotate`,
  `residual`, `orthonormalize`, `restart` (new stage: Gram defect and contraction copy; a
  re-orthonormalization at a restart stays under `orthonormalize`), `return`, `other` (the solve wall not
  covered: the correction loop, copies, waits outside stages). Collectives are inside the stages and are
  listed beside them (tree sums mostly inside `orthonormalize` and `reduced_gram`, alltoalls inside
  `apply`).
- Every rank thread keeps its own accumulators (per-thread mirrors of the stage, collective and operator
  counters, `dft_numerics::stage_clock::local`); the report takes max and min over ranks. With
  `BAANDI_STAGE_TIMERS=off` the counts and the map clocks remain and the stages read 0 (count-only mode).
- Outputs: dft-cli `results.json` (`phase_profile`; `parallel.phases` on MPI runs, gathered with one
  allgatherv), `tests/m2_entry_cost.rs` (console table and `M2SCREEN_JSON` field `profile`),
  `davidson_bench` (console table, `DAVBENCH_JSON` field `phase_profile`).

## Identity (C): the switches never change the result
`m2_entry_cost` production + ladder, np=4 simulated, main d5d04b1 harness (sha256 d261ce52...) against this
branch (04cdfc9c...): si8 energy bits 0xc03f52607f2b73b5 / density FNV 0xea947c41239431a0, si32
0xc05fad02a3080b23 / 0x39b5916d4d315b41, identical for main, branch with every switch off, and branch
with `BAANDI_PHASE_PROFILE=on BAANDI_FIRST_BLOCK=on BAANDI_SNAPSHOT_DIR=...` all on. dft-cli si32 f1e-3
seed 601 T=1: density.bin SHA-256 9491950d470503bf... with and without the three switches (the value of
mpi-remaining.md section 7). Engine test `tests/instrument.rs`: energy and density bits equal with the
switches on and off on the threads path and on 2 simulated ranks.

## si32 f1e-3 seed 601, production (standard 5e-9 Ha, verification, C1 ladder, capacity 3), np=4 (O)
Run 1 (wall 3.100 s; run 2: 2.890 s, same counts; count-only mode run: 2.725 s). Seconds = max over
ranks, `[min]` = min over ranks; counts per rank (identical on every rank).

| bucket | wall | [min] | maps | share of wall | potential | solve | density | energy | mixing |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| init (construct 0.077 + atomic start 0.103) | 0.180 | 0.176 | 0 | 5.8 % | | | | | |
| coarse | 0.566 | 0.563 | 9 | 18.3 % | 0.026 | 0.453 | 0.049 | 0.028 | 0.010 |
| fine | 1.901 | 1.899 | 11 | 61.3 % | 0.076 | 1.567 | 0.139 | 0.087 | 0.036 |
| verify_pass | 0.177 | 0.177 | 1 | 5.7 % | 0.006 | 0.150 | 0.013 | 0.008 | |
| verify_fail | 0.177 | 0.177 | 1 | 5.7 % | 0.007 | 0.151 | 0.012 | 0.008 | |
| materialize (lazy gather) | 0.007 | 0.006 | 0 | 0.2 % | | | | | |
| covered / remainder | 3.008 / 0.093 | | | 97.0 / 3.0 % | | | | | |

Solve stages (exclusive; seconds max [min]):

| bucket | start | apply | reduced_gram | reduced_eigen | rotate | residual | orthonormalize | restart | return | other |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| coarse (0.453) | 0.006 | 0.179 [0.178] | 0.030 [0.029] | 0.065 [0.062] | 0.039 [0.033] | 0.017 [0.010] | 0.063 [0.063] | 0.001 | 0.000 | 0.053 [0.069] |
| fine (1.567) | 0.019 | 0.674 [0.671] | 0.135 [0.132] | 0.096 [0.095] | 0.189 [0.179] | 0.045 [0.034] | 0.241 [0.241] | 0.011 | 0.002 | 0.156 [0.184] |
| verify_pass (0.150) | 0.002 | 0.068 [0.066] | 0.014 [0.013] | 0.008 [0.007] | 0.019 [0.013] | 0.009 [0.002] | 0.025 [0.025] | 0.002 | 0.000 | 0.003 [0.021] |
| verify_fail (0.151) | 0.002 | 0.068 [0.066] | 0.015 [0.013] | 0.006 [0.005] | 0.020 [0.014] | 0.008 [0.002] | 0.027 [0.027] | 0.002 | 0.000 | 0.003 [0.021] |

Counts per rank (C): coarse 48 Davidson iterations, H 1723 columns (1179 full + 544 delta-V), restarts 2
(2 contractions), kept starts 8, FFT pairs 868 forward / 900 inverse, tree sums 1495 (0.058 s),
alltoalls 105 (0.035 s). Fine 67 iterations, H 2568 (1888 full + 680 delta-V), restarts 7, kept 10, FFT
pairs 1295 / 1327, tree sums 2464 (0.143 s), alltoalls 145 (0.098 s). Verification pass 7 iterations,
257 H (189 + 68); fail 5 iterations, 267 H (199 + 68); one tree sum chain of 260 / 287.

Reading (O): the fine solve is 50 % of the wall and the coarse solve 15 %; inside the fine solve the H
application is 43 %, orthonormalization 15 %, rotation 12 %, reduced border 9 %, EVD 6 %, residual 3 %,
the uncovered loop body 10 %. The two verification maps (one failed, one passed) are 11 % of the wall,
exactly the review's "failed verification attempts" line. The rank imbalance of the exclusive stages is
small (max against min within 2 to 5 percent) except `residual` and `other`, which absorb the waits at
the first synchronisation after the rotation (the slowest rank has the smallest `other` and the largest
`residual`: 0.045 against 0.034 s). The earlier mpi-remaining.md table (same case) listed "outside the
recorded fine phases 1.28 s"; that number is now resolved into coarse 0.566, verification 0.354, init
0.180, mixing 0.046 and a 0.09 s remainder.

## si16 np=4 (O; wall 0.759 s)
init 0.047 (atomic start 0.018), coarse 0.161 (10 maps, solve 0.113), fine 0.407 (9 maps, solve 0.306),
verify_pass 0.041 (1), verify_fail 0.083 (2 maps), materialize 0.001, remainder 0.019 (2.5 %). Fine
solve: apply 0.143, orthonormalize 0.051, rotate 0.024, reduced_gram 0.019, reduced_eigen 0.018,
residual 0.007, start 0.003, restart 0.002, other 0.039; 52 iterations, H 1195 (907 + 288), restarts 8,
tree sums 1456 (0.037 s), alltoalls 113 (0.025 s). Two failed verifications: 11 % of the wall in failed
attempts, 5 % in the passed one.

## si8 np=4 (O; wall 0.194 s)
init 0.016, coarse 0.029 (5 maps), fine 0.119 (10 maps, solve 0.073), verify_pass 0.012, verify_fail
0.012, remainder 0.006 (3.1 %). Fine solve: apply 0.033, orthonormalize 0.012, reduced_eigen 0.005,
reduced_gram 0.003, rotate 0.003, residual 0.002, other 0.014; 50 iterations, H 584 (404 + 180),
restarts 3, tree sums 780, alltoalls 110.

## Threads path, si32 T=1, dft-cli (O; separate observation of the same switch, not a rank comparison)
Wall 6.327 s: init 0.272 (atomic start 0.181), coarse 1.089 (7 maps, solve 0.924), fine 4.328 (10 maps,
solve 3.776), verify_pass 0.509 (1 map, solve 0.460), materialize 0.005 (output writing), remainder
0.124 (2.0 %). Coarse 34 Davidson iterations / 1018 H columns, fine 48 / 1717, verification 5 / 201; no
delta-V columns are counted on this path (REUSE-PSI's kept H X update is a forward transform of kept
psi(r), counted in `fft.kept_transforms`, not in the operator's change-column counter). The serial solve
stages fill only with `BAANDI_DAVIDSON_STAGES=on` (not set here: `solve.other` = the whole solve).

## What can now be seen that could not before
1. One table closes the wall: covered 97 to 98 percent, the rest named. The review's warning that
   "25 % coarse + 12 % orthogonalization + 6 % EVD" were overlapping shares is answered by exclusive
   numbers: at si32 np=4 (O) the coarse bucket is 18 %, and the orthonormalization inside the FINE solve
   alone is 7.8 % of the wall (0.241 s), inside the coarse solve 2.0 %.
2. The H application is split into full columns and delta-V columns per bucket (fine: 1888 + 680), with
   the forward / inverse pair counts beside them, so a "27 % fewer H columns" claim can be weighed by the
   transforms it removes (a delta-V column costs one forward transform, not a pair).
3. Verification cost by outcome: failed attempts are a separate line (si32: 5.7 %, si16: 11 %).
4. Max against min over ranks per stage shows where the waits sit (residual / other).

Raw logs: `/private/tmp/claude-501/.../scratchpad/b0-instrument/B-profile-*.log` (this session's
scratchpad; the tables above are complete copies of the numbers used). Commands:
```
export DFT_SI_UPF=/path/Si.pz-vbc.UPF
env M2_PROFILE_MODE=production M2_PROFILE_LADDER=on M2_PROFILE_SERIAL=off M2_PROFILE_RANKS=4 \
    M2_PROFILE_CASE=si32 BAANDI_PHASE_PROFILE=on \
    target/release/deps/m2_entry_cost-<hash> --ignored --nocapture m2_entry_cost_profile
BAANDI_PHASE_PROFILE=on target/release/dft-cli run si32.toml --output-dir out      # threads path
```
