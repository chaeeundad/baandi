# What was in `solve.other`, `init` and the remainder (BREAKTHROUGH-4-OTHER-AND-INIT, 2026-11-17)

Starting point (relayed, dev08, real MPI, si32 np=8 production, a220a85, wall 3.3345 s): `solve.other` of
the fine solve 0.2601 s (24.7 % of the fine solve, the second largest term), the same pattern in the coarse
(0.0907 of 0.3235) and verification solves, `init` 0.5606 s (atomic start 0.2048 s), 0.2552 s (7.7 %) of the
wall outside every bucket. This note classifies all three with finer labels of the exclusive phase profile
(`BAANDI_PHASE_PROFILE=on`, profile.rs) and records the cheap fixes the labels exposed. Labels: **O** Mac
observation (pure build, simulated ranks = threads of one process, shared 12-CPU Mac at load 4 to 10:
seconds are not dev08 numbers, the shares and the identity are what to read), **C** exact count. dev08
request: `docs/performance/requests/breakthrough-other-init.md`.

## New labels (all exclusive unless stated; count-only mode keeps the calls; results bit-identical)
Solver (`dft_solvers::distributed_davidson::STAGES` 11 to 15):
- `solve.setup`: the per-solve row-index exchange of `solve_rows_impl` (two allgathervs, one of N doubles,
  and the index rebuild) and the workspace preparation before `solve_start`.
- `solve.precondition`: convergence flags, selected-root bookkeeping and coefficient copies, the
  preconditioned correction block.
- `solve.append`: copies of the new columns into the basis and product storage, the finiteness scan of
  every applied block.
- `solve.solver_total`: the whole solver call (container). `solve.other` is now split into
  `other_solver` (the solver call minus all its stages) and `other_engine` (the rest of the solve window).
- `solve.wait` (NOT exclusive, `BAANDI_PHASE_WAIT=on` only): a timed barrier before every collective of
  the solver; the time a rank waits for the last one, reported apart from the collective's transfer.
Engine (m2 stages `kept_update`, `solve_entry`, `solve_finish`, `entry_wait`, `real_state`):
- `solve.real_state`: `ensure_state` of a real Gamma solve, i.e. the replicated SERIAL real Hamiltonian
  (reused with a potential copy, or rebuilt after a new problem / ladder stage).
- `solve.kept_update`: the MPI-KEPT-START update H_new X = H_old X + (V_new - V_old) X (the delta-V
  columns: forward transforms of the kept slabs, both ways for the pairs not kept).
- `solve.entry` / `solve.finish`: the interface between the preparation and the solver call (local
  diagonal, start rows, kept-start bookkeeping) and after it (lazy handle copy of the local rows, old
  potential copy). `solve.prepare` (`m2.solve_prepare`) and `solve.unpack` were already timed but were
  not subtracted from `other`; they are now.
- `solve.entry_wait` (`BAANDI_PHASE_WAIT=on` only): a timed barrier at the entry of every M2 solve
  (arrival skew of the ranks after the potential phase).
Init (`init` bucket): construction parts `cli` (input, preflight, prepared basis, UPF reading of the
CLI), `basis`, `ionic_local`, `densities` (FFT plan, core and start density, real ionic potential),
`states` (pool, k states: basis, nonlocal projectors, Hamiltonian, band-density workspace),
`atomic_prepare`, `ewald`, `rest`; atomic start `span` (= `radial` + `fill` + `orthonormalize` +
`finish`), `hq` (threads path, ATOMIC-REAL) and `ritz` (= `gram` + `eigh` + `rotate`). The cutoff-ladder
COARSE problem's construction now belongs to `init` too (it was in the remainder).
New bucket `between`: `ladder_handoff` (orbital and density prolongation), `ladder_truncation` (the
report's atomic truncation distance), `session_setup` (session construction, mixer validation, imported
state and seeds), `stop_rule` (stop-rule estimate, density residual and iteration record outside the
map), `verify_glue` (verification block outside its map).
The per-rank `solve.other` is now derived on each rank before the max/min over ranks (a sum of per-stage
maxima over different ranks exceeded the solve and made `other` read low; the old table's `[min] > max`
in `other` was this artefact).

## What was hiding in `solve.other` (O, si32 np=8 simulated, production + ladder, seed 601)
Same binary (branch, all switches off), profile on. Fine solve 1.333 s (11 maps), coarse 0.515 s (9 maps):

| label | fine | coarse | verify (pass + fail) | what it is |
|---|---:|---:|---:|---|
| `real_state` | 0.1095 | 0.0548 | 0.0004 | rebuild of the replicated serial real Hamiltonian at the first map of each ladder stage: capacity / 2 = 102 full-grid pair FFT plans and scratch, real projector table, per rank, never applied on the M2 path |
| `kept_update` | 0.0720 | 0.0264 | 0.0143 | delta-V columns of the kept start (680 fine, 544 coarse; real FFT work, already counted as `h_delta_columns`) |
| `setup` | 0.0238 [0.0098] | 0.0093 [0.0033] | 0.0013 | row-index exchange + workspace preparation; max/min 2.4: mostly arrival skew absorbed by the first collective |
| `prepare` | 0.0228 | 0.0084 | 0.0017 | operator update, `DistributedRealHamiltonian::new` at the first map of a stage, diagonal, lazy-start check |
| `finish` + `entry` + `append` + `precondition` | 0.0114 | 0.0024 | 0.0023 | lazy handle copy, local diagonal, column appends, correction block |
| remaining `other` | 0.0005 | 0.0003 | 0.0001 | (was 0.112 fine / 0.066 coarse before `real_state` was labelled) |

Before the labels (same case, earlier binary of the branch): fine `other` 0.188 s, coarse 0.095 s. On
dev08 the fine `other` was 0.260 s: the share is larger there, as expected from the replicated serial
operator build (memory allocation and FFT plan creation per rank at np=8 on one node). With
`BAANDI_PHASE_WAIT=on` (diagnostic run, extra barriers): `entry_wait` fine 0.022 s, coarse 0.016 s; the
solver's waits before its collectives 0.249 s fine (2547 probe barriers) and 0.144 s coarse (1562), inside
the stages that issue the collectives (simulated ranks on a loaded Mac: the waiting share is not a dev08
number). The wait split is the question Job 1 of the
request answers on real MPI.

After merging main b1c95a6 (`M2_LEAN_ALGEBRA=on` by default; the lean path got the same
`precondition`/`append` labels) the picture is the same (O, si32 np=8): fine `other` 0.0006 s,
`real_state` 0.118 s, `kept_update` 0.075 s, `prepare` 0.025 s, `setup` 0.016 s.

## Init and remainder (O, same run)
`init` 0.362 s = fine + coarse construction (basis 0.028, ionic_local 0.044, densities 0.004, states 0.113,
ewald 0.067, rest 0.013; max over ranks, noisy at 8 simulated ranks: `states` reads 0.03 to 0.15 s between
runs, allocation-sensitive) + atomic start 0.108 s (span 0.062 of which the distributed Gram-Schmidt chain
0.058; Ritz 0.047 = gram 0.030 + eigh 0.004 + rotate 0.013). `between` 0.078 s: ladder hand-off 0.028,
session setup 0.036, stop rule 0.019, verification glue 0.002. Remainder of the run wall 0.3 % (was 3.0 to
5.6 % on the Mac and 7.7 % on dev08). Threads path (CLI si32 T=1, O, wall 5.45 s): init 0.341 s = CLI
0.031 + engine construction (basis 0.007, ionic_local 0.026, states 0.021, ewald 0.041 for two
constructions) + atomic start 0.206 s (orthonormalize 0.069, H span 0.090, Ritz 0.033); between 0.031 s;
remainder 0.2 %.

## Cheap fixes (each behind a switch, default off, bit-identical; ledger rows BT4-*)
1. `BAANDI_M2_LIGHT_STATE=on`: on the M2 path the serial real Hamiltonian only supplies its ORC half
   sphere (packing) and is never applied, so it is built with block capacity 2 (one pair plan) instead of
   the Davidson capacity; if the distributed operator cannot be built (the M1 fallback applies the serial
   operator) the full one is built. `real_state` fine 0.1095 -> 0.0107 s, coarse 0.0548 -> 0.0026 s (O).
2. `BAANDI_EWALD_CACHE=on`: a construction with bit-equal lattice, charges and fractional positions on the
   same thread reuses the Ewald energy (the ladder's coarse problem after the full-cutoff one); `init.ewald`
   0.067 -> 0.043 s (O; one of the two evaluations).
3. `BAANDI_M2_LAYOUT_CACHE=on`: the solver reuses the row-index layout of its previous solve when every
   rank's local rows are unchanged (one 1-double reduction decides); allgatherv calls 117 -> 77 per rank and
   1.4 MB less payload per run, tree sums +22 (C). Wall neutral on simulated ranks (an allgatherv there is a
   memcpy); a latency item for dev08.

Wall, alternating pairs, profile off (O; same binary, environment A/B; energy bits and density FNV
identical in every run):

| case | arm B | pairs | A median s | B median s | B/A median | B/A per pair |
|---|---|---:|---:|---:|---:|---|
| si32 np=8 | stack (1+2+3) | 7 | 3.097 | 2.934 | 0.936 | 0.885 0.952 0.968 0.932 0.911 0.961 0.936 |
| si32 np=8 | light state | 5 | 3.154 | 2.895 | 0.918 | 0.866 1.073 0.915 0.961 0.918 |
| si32 np=8 | Ewald cache | 5 | 3.055 | 3.064 | 0.985 | 1.018 0.971 0.985 1.020 0.908 |
| si32 np=8 | layout cache | 5 | 3.098 | 3.063 | 1.000 | 1.015 1.000 0.980 1.015 0.925 |
| si32 np=8, main b1c95a6 defaults (M2_LEAN_ALGEBRA on) | stack | 3 | 2.993 | 2.846 | 0.951 | 0.951 0.939 1.017 |
| si16 np=4 | stack | 3 | 0.842 | 0.807 | 0.959 | 0.964 0.959 0.942 |
| si8 np=2 | stack | 3 | 0.298 | 0.275 | 0.930 | 0.930 0.922 0.957 |

Peak RSS of the si32 np=8 process (all 8 simulated ranks): 2.2 to 2.3 GB default, 2.4 to 2.6 GB light state
(three runs each, O): no memory claim; the untouched scratch of the full serial build is not resident.

## Not done here (backlog, with the numbers that justify them)
- `kept_update` (fine 0.072 s O) is real transform work (one forward pair transform per kept pair, both
  ways for the 2 pairs that are not kept: the partly occupied pair and the empty pair force the mirror
  exchange and a transpose to r). Keeping the slab of the partly occupied pair (it is transformed by the
  density build anyway) would leave one pair; the exchange stays.
- Atomic start at np=8: the distributed two-pass Gram-Schmidt of the 128 span columns is a serial chain of
  128 broadcasts (0.058 s O against 0.069 s serial on the threads path: the split barely helps); a block
  orthonormalization (Gram + Cholesky with the same dependence rule) would change the span's rounding and
  therefore the whole trajectory (rounding class, needs stage 1 to convergence). The Ritz rotation
  Q x C (npw x 128 x 68, replicated, 0.013 s O) can be split by band columns with one allgatherv.
- Replicated construction work per rank (O, two constructions): `ionic_local` 0.044 s (independent per
  G: a G-range split plus one allgatherv is bit-identical), the remaining Ewald evaluation (reciprocal sum
  over G: a split changes the summation order of the energy only, rounding class of E, the density is
  untouched), `DistributedRealHamiltonian::new` packs every projector over all N before keeping the local
  rows.
- `session_setup` 0.036 s and `stop_rule` 0.019 s (O): replicated full-grid passes (mixer validation,
  imported-state checks, density residual); small.

Raw logs: this session's scratchpad `bt4/` (b7-A.log, b7-B.log, b2-W.log, stackLEC.txt, alone*.txt,
si16np4.txt, si8np2.txt, cli-T1.log); the tables above are complete copies of the numbers used.
Commands:
```
export DFT_SI_UPF=/path/Si.pz-vbc.UPF
env M2_PROFILE_MODE=production M2_PROFILE_LADDER=on M2_PROFILE_SERIAL=off M2_PROFILE_RANKS=8 \
    M2_PROFILE_CASE=si32 BAANDI_PHASE_PROFILE=on [BAANDI_PHASE_WAIT=on] \
    target/release/deps/m2_entry_cost-<hash> --ignored --nocapture m2_entry_cost_profile
BAANDI_PHASE_PROFILE=on target/release/dft-cli run si32.toml --output-dir out      # threads path
```
