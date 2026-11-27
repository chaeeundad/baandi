# Reference-path switch retirement (R11a)

R12 status: [ADR 0047](../decisions/0047-reference-driver-single-path.md) records
the accepted D3 policy. The table below is the **R11a source-commit inventory**,
not the current switch manifest. R10c subsequently removed `REFERENCE_SCF` and
legacy standard controls. At R12 base `27d91ece`, R11b is ticketed but not merged;
its oracle/shared entries marked “Keep for R11b” must not be reported as removed.
Use `crates/dft-numerics/src/switches.rs` for actual accepted names. Recreate
retired experimental branches from git; D3 requires no rollback waiting period.

Source/recreation commit: `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1`. User decision D3 supersedes the two-generation sunset suggestion and the standing rollback policy. This is code retirement, not a new performance claim.

The inventory below covers the entire pre-change manifest, including legacy-only entries, so shared-kernel and indirect oracle readers cannot disappear from the audit. Reader locations are **pre-change** source locations; repeated reads in one file are grouped. All readers were searched in the checkout; `Id::NAME` alone misses bare imports, and string pins in tests were searched separately. Historical ledger rows remain unchanged.

## Decisions

- Remove 29 reference-only default-on rollback switches and six unadopted/rejected reference experiments. Keep the adopted arithmetic and its runtime safety fallbacks (rank/pivot guards, overflow-safe hypot, missing-cache cold solve, and backend dispatch).
- Keep all switches in shared numerical/physics kernels callable by the strict oracle, including FFT and nonlocal optimizations; R11b owns their retirement after R9c. No legacy driver or legacy-only solver file changes.
- Keep RMM, reference ladder, Kerker/TF and backend selection, diagnostics, infrastructure, REFERENCE_STRICT, REFERENCE_SCF and XC_SCREEN_NEGATIVE.
- Keep REFSCF_KEPT_RECOMPUTE: the ledger explicitly retains it as a memory-saving option, rather than rejecting it. REFSCF_PERMAP and REFSCF_MEMORY_TRIM have no rejection/supersession row; absence of adoption is insufficient to retire them.
- Keep REFSCF_SPIN_FD_MOMENT: a routing reader is in `dft-cli/src/driver_select.rs`, outside this ticket.
- Historical `VerificationEnergy` payload/accessors and the always-None pending-map checkpoint slot remain for checkpoint compatibility; no production energy-acceptance or failed-map-reuse branch remains. Restore rejects a nonempty retired pending-map slot.

## Complete inventory

| Switch (former/current default) | Class | Readers at source commit | Ledger evidence | Decision |
| --- | --- | --- | --- | --- |
| `INIT_OWNED_STATES` (on) | oracle | `crates/dft-engine/src/scf.rs:991,1017` | Mac observation and retained design (O, not adoption) (line 1611) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `CERT_SCREEN` (on) | oracle | `crates/dft-engine/src/session.rs:2666` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 875) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `POTENTIAL_REUSE` (on) | oracle | `crates/dft-engine/src/session.rs:2671` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 967) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `MIX_GROWTH_RESET` (none) | oracle | `crates/dft-engine/src/session.rs:2672` | P47B — fp32 H feasibility (2026-11-08; measurement only, Mac observations) (line 790) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `MIX_CELL_KERKER` (off) | algorithm module | `crates/dft-engine/src/scf.rs:905,916,922` | G-MIXER-DIVERGENCE (2026-11-23; Mac counts/correctness only) (line 3121) | Keep: Generality rule/backend selection. |
| `EIGEN_FACTOR` (0.1) | oracle | `crates/dft-cli/src/run_config.rs:160`; `crates/dft-engine/src/session.rs:2673` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 873) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `EMPTY_LOOSE` (on) | oracle | `crates/dft-engine/src/scf.rs:3698` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 901) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `BUFFER_TOL` (1e-3) | oracle | `crates/dft-engine/src/scf.rs:3708` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 901) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `NL_ENERGY` (fast) | oracle | `crates/dft-engine/src/energy.rs:177` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 939) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `ENERGY_IDENTITY` (off) | oracle | `crates/dft-engine/src/energy.rs:258` | LANE-GLUE serial density fusion — 2026-11-20, default-off, dev08 request (line 2898) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `STANDARD_STOP` (plain) | oracle | `crates/dft-engine/src/session.rs:2655` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 887) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `STANDARD_RESIDUAL_GUARD` (on) | oracle | `crates/dft-engine/src/session.rs:2656`; `crates/dft-engine/src/finite_temperature.rs:186`; `crates/dft-engine/src/collinear.rs:893,1309`; `crates/dft-engine/src/refscf/electronic.rs:509` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 885) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `DENSITY_PAIR` (on) | oracle | `crates/dft-engine/src/band_density.rs:19` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 938) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `DENSITY_REDUCE` (on) | oracle | `crates/dft-engine/src/band_density.rs:28` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 868) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `NL_ENERGY_REAL` (on) | oracle | `crates/dft-engine/src/gamma_real.rs:429` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 939) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `ATOMIC_LOCAL_BLOCK` (on) | oracle | `crates/dft-engine/src/atomic_start.rs:161` | LANE-INIT-ATOMIC-LOCAL-BLOCK-DEFAULT-ON (2026-11-19/20, dev08 via master; ledger by orchestrator) (line 2586) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `ATOMIC_REAL` (on) | oracle | `crates/dft-engine/src/gamma_real.rs:667` | Full complex table consumers at this revision (line 1371) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `MPI_M2` (on) | oracle | `crates/dft-cli/src/run.rs:683,695,699`; `crates/dft-engine/src/parallel.rs:90` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 871) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_BATCH` (on) | oracle | `crates/dft-cli/src/run.rs:2407`; `crates/dft-engine/src/m2.rs:60` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 880) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `PHASE_THREADS` (on) | oracle | `crates/dft-engine/src/scf.rs:2513`; `crates/dft-engine/src/collinear.rs:489` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 915) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `REUSE_PSI` (on) | oracle | `crates/dft-cli/src/run.rs:787,791,799,4082`; `crates/dft-engine/src/gamma_real.rs:146,1201,1210`; `crates/dft-engine/src/snapshot.rs:351,391` | Fast-loop screens (2026-11-16; threads/one-core owner, Mac observations O, `scripts/fast_screen.py`, si32-gamma f1e-3 seed 601 T=1 unless stated, main 62c821f pure build; profile docs/performance/profile-2026-11-16-threads.md; dev08 request docs/performance/requests/verify-reuse.md) (line 983) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `REF_DAV_RESTART_CORRECTIONS` (on) | rollback default-on | `crates/dft-solvers/src/ref_davidson.rs:225`; `crates/dft-solvers/src/ref_davidson_complex.rs:237` | P3 timed (line 3539) | Remove: D3; preserve default arithmetic. |
| `REFSTART_MKL` (on) | rollback default-on | `crates/dft-engine/src/refstart.rs:1136` | REFSTART_MKL (2026-11-24): ADOPTED (line 3561) | Remove: D3; preserve default arithmetic. |
| `REF_DAV_WIDE_ROTATE` (on) | rollback default-on | `crates/dft-solvers/src/ref_davidson_complex.rs:246`; `crates/dft-solvers/src/ref_davidson/workspace.rs:1119` | P3 timed (line 3539) | Remove: D3; preserve default arithmetic. |
| `REFSCF_VERIFY_REUSE` (off) | experiment default-off | `crates/dft-engine/src/refscf/driver.rs:266`; `crates/dft-engine/src/refscf/electronic.rs:332`; `crates/dft-engine/src/refscf/k_driver.rs:146` | I25-VERIFY-COST (line 3503) | Remove: D3; preserve default arithmetic. |
| `REFSCF_VERIFY_GATE` (on) | rollback default-on | `crates/dft-engine/src/refscf/driver.rs:268,270`; `crates/dft-engine/src/refscf/electronic.rs:339`; `crates/dft-engine/src/refscf/k_driver.rs:148,150` | I25-VERIFY-COST (line 3503) | Remove: D3; preserve default arithmetic. |
| `REFSCF_VERIFY_GATE_FIRST` (off) | experiment default-off | `crates/dft-engine/src/refscf/driver.rs:270,278`; `crates/dft-engine/src/refscf/electronic.rs:340`; `crates/dft-engine/src/refscf/k_driver.rs:150,158` | I25-FIRST-ATTEMPT (line 3525) | Remove: D3; preserve default arithmetic. |
| `VERIFY_REUSE` (half) | oracle | `crates/dft-engine/src/session.rs:2658` | Run-owned configuration, test harness, validated report (RUN-CONFIG / TEST-HARNESS / RUN-REPORT, 2026-11-18; Mac observations O, pure release build, si32-gamma seed 601 amp 1e-3, T=1, load 3 to 16; bits C; base main 1ddcc07; docs/performance/run-config.md, docs/performance/test-targets.md) (line 1236) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_DGEMM` (on) | oracle | `crates/dft-physics/src/distributed_real_gamma.rs:443`; `crates/dft-solvers/src/distributed_davidson.rs:812` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 930) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_ORTHO` (serial) | oracle | `crates/dft-solvers/src/distributed_davidson.rs:813` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 930) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_FFT` (rust) | algorithm module | `crates/dft-solvers/src/distributed_fft.rs:568` | Work-order screen (2026-11-26, user work order "reduce full-cost SCF work"), dev08 Si64 MPI8, one run per arm at bb416c73 (line 3680) | Keep: Generality rule/backend selection. |
| `M2_LEAN` (on) | oracle | `crates/dft-solvers/src/distributed_davidson.rs:816` | Measurement infrastructure (BREAKTHROUGH-0-INSTRUMENT, 2026-11-17; review 2026-11-17 sections 13.1 to 13.3 and E3; Mac observations O on simulated ranks np=4 and T=1, counts exact C; notes docs/performance/profile-2026-11-17-phases.md, docs/performance/first-block-residual-2026-11-17.md; dev08 request docs/performance/requests/breakthrough-instrument.md). No speed claim; three default-off switches whose on/off state leaves energy bits and density bits unchanged (main d5d04b1 harness against this branch, si8/si32 np=4 production + ladder; dft-cli si32 T=1 density SHA 9491950d...; test `dft-engine/tests/instrument.rs`). (line 1101) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_CAPACITY` (auto) | oracle | `crates/dft-cli/src/run.rs:393` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 958) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `STAGE_TIMERS` (on) | diagnostic | `crates/dft-numerics/src/real.rs:2456`; `crates/dft-numerics/src/stage_clock.rs:24`; `crates/dft-numerics/src/stage_clock/arrival.rs:53` | LANE-GLUE Phase 0 (2026-11-20): PROFILE_REQUEST (line 2699) | Keep: R11(e). |
| `M2_EIGEN` (replicated) | oracle | `crates/dft-solvers/src/distributed_davidson.rs:817` | INIT-SCALE: making `init` scale with MPI ranks (2026-11-18; Mac observations O on 8 simulated ranks, load 3 to 5, pure build; counters exact C; branch base main 0c4d652, commits 03fd272, 5262fcf, a7205c9, 529e485) (line 1164) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_B5` (off) | oracle | `crates/dft-solvers/src/distributed_davidson.rs:818` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 930) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_WORKSPACE` (on) | oracle | `crates/dft-solvers/src/distributed_davidson.rs:819` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 964) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_SELECTED_ROOTS` (on) | oracle | `crates/dft-solvers/src/distributed_davidson.rs:820` | Stage-0 Davidson screens (FASTBENCH-DAVIDSON, 2026-11-16; Mac observations O, counters exact; `scripts/fast_screen_dav.py`, one solve at a fixed Hamiltonian, si32-gamma f1e-3 seed 601 T=1, base density = converged `dft-cli run` density.bin, cold = atomic start at 5e-3 Ha, warm = kept start at the ratchet tolerance 3.6e-11 Ha, 3 alternating pairs x 5 repeats, pure build, load 3 to 5; protocol docs/performance/fast-loop-protocol.md stage 0) (line 1008) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_PANEL_GS` (fused) | oracle | `crates/dft-solvers/src/distributed_davidson.rs:821` | MPI lazy orbitals (MPI-LAZY-ORBITALS, 2026-11-17; Mac observations O on simulated ranks, counters exact C, relayed dev08 D, estimates E; design docs/design/mpi-lazy-orbitals.md, dev08 request docs/performance/requests/mpi-lazy-orbitals.md; review F07 second half) (line 1033) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_DENSITY` (fused) | oracle | `crates/dft-engine/src/scf.rs:2066` | ONEPATH-5-SPIN-SHARE (2026-11-17) (line 1127) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_REUSE` (on) | oracle | `crates/dft-engine/src/gamma_real.rs:251`; `crates/dft-engine/src/snapshot.rs:351,388` | 2026-11-20 — LANE-VERIFY Phase 2: closed without a qualifying lever (line 2851) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_STICKS` (mirror) | oracle | `crates/dft-solvers/src/distributed_fft.rs:261` | MPI remaining terms (MPI-REMAINING, 2026-11-17; Mac observations O on simulated ranks, counters exact C, relayed dev08 D, estimates E; analysis docs/performance/mpi-remaining.md, dev08 request docs/performance/requests/mpi-remaining.md) (line 1046) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_START` (gram) | oracle | `crates/dft-solvers/src/distributed_davidson.rs:873` | MPI overlap-matrix Davidson (MPI-OVERLAP-DAVIDSON, 2026-11-17; Mac observations O on simulated ranks, counters exact C, estimates E; design docs/design/mpi-overlap-davidson.md, dev08 request docs/performance/requests/mpi-overlap-davidson.md; review F05 step 2) (line 1063) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_LEAN_ALGEBRA` (on) | oracle | `crates/dft-engine/src/snapshot.rs:351,389`; `crates/dft-solvers/src/distributed_davidson.rs:826` | Measurement infrastructure (BREAKTHROUGH-0-INSTRUMENT, 2026-11-17; review 2026-11-17 sections 13.1 to 13.3 and E3; Mac observations O on simulated ranks np=4 and T=1, counts exact C; notes docs/performance/profile-2026-11-17-phases.md, docs/performance/first-block-residual-2026-11-17.md; dev08 request docs/performance/requests/breakthrough-instrument.md). No speed claim; three default-off switches whose on/off state leaves energy bits and density bits unchanged (main d5d04b1 harness against this branch, si8/si32 np=4 production + ladder; dft-cli si32 T=1 density SHA 9491950d...; test `dft-engine/tests/instrument.rs`). (line 1101) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_LAZY` (on) | oracle | `crates/dft-engine/src/gamma_real.rs:271`; `crates/dft-engine/src/snapshot.rs:351,390` | Between, verification glue and materialize (BETWEEN-MATERIALIZE, 2026-11-18; Mac observations O on 8 simulated ranks, pure build, loaded Mac: buckets only, walls not comparable; counters exact C; base main e902853, commits 5fe2796, b3f8509, 3f7d82b) (line 1188) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `MPI_TREE` (rd) | infra | `crates/dft-cli/src/main.rs:355` | INIT-SCALE: making `init` scale with MPI ranks (2026-11-18; Mac observations O on 8 simulated ranks, load 3 to 5, pure build; counters exact C; branch base main 0c4d652, commits 03fd272, 5262fcf, a7205c9, 529e485) (line 1165) | Keep: routing/infrastructure. |
| `M2_SCREEN` (dist) | oracle | `crates/dft-engine/src/parallel/tuning.rs:22` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 935) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_ATOMIC` (split) | oracle | `crates/dft-engine/src/parallel/tuning.rs:23` | INIT-SCALE: making `init` scale with MPI ranks (2026-11-18; Mac observations O on 8 simulated ranks, load 3 to 5, pure build; counters exact C; branch base main 0c4d652, commits 03fd272, 5262fcf, a7205c9, 529e485) (line 1159) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_ATOMIC_BLOCK` (on) | oracle | `crates/dft-engine/src/parallel/tuning.rs:24` | INIT-SCALE: making `init` scale with MPI ranks (2026-11-18; Mac observations O on 8 simulated ranks, load 3 to 5, pure build; counters exact C; branch base main 0c4d652, commits 03fd272, 5262fcf, a7205c9, 529e485) (line 1159) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `STOP_METRIC` (kept) | oracle | `crates/dft-cli/src/run.rs:771,783,785`; `crates/dft-engine/src/scf.rs:235` | Run-owned configuration, test harness, validated report (RUN-CONFIG / TEST-HARNESS / RUN-REPORT, 2026-11-18; Mac observations O, pure release build, si32-gamma seed 601 amp 1e-3, T=1, load 3 to 16; bits C; base main 1ddcc07; docs/performance/run-config.md, docs/performance/test-targets.md) (line 1231) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_GRID` (on) | oracle | `crates/dft-cli/src/run.rs:701,708,721`; `crates/dft-cli/src/run_config.rs:159`; `crates/dft-engine/src/grid_dist.rs:63` | Owned potential and density of the M2 path (M2-OWNED, 2026-11-18; base 6be8f09; Mac counters exact C, wall O; dev08 request docs/performance/requests/m2-owned.md) (line 1218) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_MIX` (on) | oracle | `crates/dft-cli/src/run.rs:723,727,740`; `crates/dft-engine/src/mix_dist.rs:51` | Owned potential and density of the M2 path (M2-OWNED, 2026-11-18; base 6be8f09; Mac counters exact C, wall O; dev08 request docs/performance/requests/m2-owned.md) (line 1226) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_LOCAL_FIELDS` (on) | oracle | `crates/dft-cli/src/run.rs:748`; `crates/dft-engine/src/grid_dist.rs:101` | G-FILTER-REPLAY — resident block Chebyshev, stage 1 (2026-11-23) (line 3097) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_LOCAL_DENSITY` (on) | oracle | `crates/dft-cli/src/run.rs:753`; `crates/dft-engine/src/grid_dist.rs:108` | Owned potential and density of the M2 path (M2-OWNED, 2026-11-18; base 6be8f09; Mac counters exact C, wall O; dev08 request docs/performance/requests/m2-owned.md) (line 1225) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `MIX_COULOMB_GUARD` (off) | oracle | `crates/dft-engine/src/session.rs:319,379` | G-MIXER-DIVERGENCE (2026-11-23; Mac counts/correctness only) (line 3122) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `VERIFY_SCHEDULE` (on) | oracle | `crates/dft-engine/src/session.rs:2674` | 2026-11-20 — LANE-VERIFY reopened: first-attempt contraction scheduling (line 2935) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `STOP_PREDICTOR` (off) | oracle | `crates/dft-engine/src/session.rs:2675` | G-MIXER-DIVERGENCE (2026-11-23; Mac counts/correctness only) (line 3122) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `MIX_INCREMENTAL` (on) | oracle | `crates/dft-solvers/src/mixing.rs:301` | MPI overlap-matrix Davidson (MPI-OVERLAP-DAVIDSON, 2026-11-17; Mac observations O on simulated ranks, counters exact C, estimates E; design docs/design/mpi-overlap-davidson.md, dev08 request docs/performance/requests/mpi-overlap-davidson.md; review F05 step 2) (line 1067) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `B1` (on) | oracle | `crates/dft-cli/src/run.rs:4013`; `crates/dft-solvers/src/real_davidson.rs:1066,1149`; `crates/dft-solvers/src/davidson_stages.rs:52` | G-PHASE-VALIDATE (2026-11-19): exact finite-dot range proof (line 2303) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `B1_SKETCH` (count) | oracle | `crates/dft-cli/src/run.rs:4014`; `crates/dft-solvers/src/davidson_stages.rs:252` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 869) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `B1_NARROW` (off) | oracle | `crates/dft-solvers/src/davidson_stages.rs:249` | No named row; manifest semantics / reader scope | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `B1_MARGIN` (0.5) | oracle | `crates/dft-solvers/src/davidson_stages.rs:250` | No named row; manifest semantics / reader scope | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `B1_ROWS` (auto) | oracle | `crates/dft-solvers/src/davidson_stages.rs:251` | No named row; manifest semantics / reader scope | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `B5` (on) | oracle | `crates/dft-cli/src/run.rs:2158`; `crates/dft-solvers/src/davidson_stages.rs:260,261` | LANE-DAVIDSON Phase 0 — requested 2026-11-20 (line 2422) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `B6` (on) | oracle | `crates/dft-solvers/src/davidson_stages.rs:268` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 945) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `RITZ_FUSED` (off) | oracle | `crates/dft-cli/src/run.rs:4102`; `crates/dft-solvers/src/davidson_stages.rs:262,263` | Per-core kernels (PERCORE-KERNELS, 2026-11-17; Mac observations O, pure build, counters exact C, relayed dev08 R; analysis docs/performance/percore-kernels-2026-11-17.md, dev08 request docs/performance/requests/percore-kernels.md; stage 0 = `scripts/fast_screen_dav.py`, si32/si64-gamma f1e-3 seed 601 T=1, converged base density, 3 alternating pairs x 5 repeats; load 2 to 21 during the session, spreads as printed) (line 1055) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `NARROW_PROJ` (off) | oracle | `crates/dft-cli/src/run.rs:4102`; `crates/dft-solvers/src/davidson_stages.rs:265` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 868) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `NARROW_REDUCED` (off) | oracle | `crates/dft-cli/src/run.rs:4103`; `crates/dft-solvers/src/davidson_stages.rs:266` | Per-core kernels (PERCORE-KERNELS, 2026-11-17; Mac observations O, pure build, counters exact C, relayed dev08 R; analysis docs/performance/percore-kernels-2026-11-17.md, dev08 request docs/performance/requests/percore-kernels.md; stage 0 = `scripts/fast_screen_dav.py`, si32/si64-gamma f1e-3 seed 601 T=1, converged base density, 3 alternating pairs x 5 repeats; load 2 to 21 during the session, spreads as printed) (line 1056) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `DAVIDSON_STAGES` (off) | diagnostic | `crates/dft-numerics/src/real.rs:2452`; `crates/dft-numerics/src/work_timers.rs:84`; `crates/dft-numerics/src/fused_pair_plan.rs:70`; `crates/dft-physics/src/distributed_real_gamma.rs:452`; `crates/dft-cli/src/run.rs:1859,4189`; `crates/dft-engine/src/refscf/driver.rs:255`; `crates/dft-engine/src/refscf/electronic.rs:510`; `crates/dft-engine/src/refscf/k_driver.rs:135`; `crates/dft-solvers/src/davidson_stages.rs:267` | I26-2-DIAGLEVEL dev08 A/B (2026-11-24): diagnostic level, no measurable wall effect (line 3488) | Keep: R11(e). |
| `SECOND_VECTOR` (off) | oracle | `crates/dft-solvers/src/real_davidson.rs:1326` | P47B — fp32 H feasibility (2026-11-08; measurement only, Mac observations) (line 819) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `ORTHO` (see manifest) | oracle | `crates/dft-cli/src/run.rs:2159`; `crates/dft-solvers/src/real_davidson.rs:2103,2968` | Stage-0 Davidson screens (FASTBENCH-DAVIDSON, 2026-11-16; Mac observations O, counters exact; `scripts/fast_screen_dav.py`, one solve at a fixed Hamiltonian, si32-gamma f1e-3 seed 601 T=1, base density = converged `dft-cli run` density.bin, cold = atomic start at 5e-3 Ha, warm = kept start at the ratchet tolerance 3.6e-11 Ha, 3 alternating pairs x 5 repeats, pure build, load 3 to 5; protocol docs/performance/fast-loop-protocol.md stage 0) (line 1003) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `MGS` (blocked, classic) | oracle | `crates/dft-solvers/src/real_davidson.rs:2142`; `crates/dft-solvers/src/distributed_davidson.rs:806` | LANE-INIT Phase 2 — local atomic block (2026-11-20) (line 2548) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `ORTHO_BLOCK` (on, off) | oracle | `crates/dft-solvers/src/real_davidson.rs:2149`; `crates/dft-solvers/src/distributed_davidson.rs:803` | Stage-0 Davidson screens (FASTBENCH-DAVIDSON, 2026-11-16; Mac observations O, counters exact; `scripts/fast_screen_dav.py`, one solve at a fixed Hamiltonian, si32-gamma f1e-3 seed 601 T=1, base density = converged `dft-cli run` density.bin, cold = atomic start at 5e-3 Ha, warm = kept start at the ratchet tolerance 3.6e-11 Ha, 3 alternating pairs x 5 repeats, pure build, load 3 to 5; protocol docs/performance/fast-loop-protocol.md stage 0) (line 1003) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `WIDE_DOT` (on) | oracle | `crates/dft-solvers/src/real_davidson.rs:2147` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 938) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `PROJECT_WIDE` (off) | oracle | `crates/dft-solvers/src/real_davidson.rs:2151,3126,3129` | Stage-0 Davidson screens (FASTBENCH-DAVIDSON, 2026-11-16; Mac observations O, counters exact; `scripts/fast_screen_dav.py`, one solve at a fixed Hamiltonian, si32-gamma f1e-3 seed 601 T=1, base density = converged `dft-cli run` density.bin, cold = atomic start at 5e-3 Ha, warm = kept start at the ratchet tolerance 3.6e-11 Ha, 3 alternating pairs x 5 repeats, pure build, load 3 to 5; protocol docs/performance/fast-loop-protocol.md stage 0) (line 1005) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `DAVIDSON` (lean) | oracle | `crates/dft-cli/src/run.rs:637,4027`; `crates/dft-engine/src/snapshot.rs:351,392`; `crates/dft-solvers/src/davidson_stages.rs:269` | I26-2-DIAGLEVEL dev08 A/B (2026-11-24): diagnostic level, no measurable wall effect (line 3488) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `LEAN_PROJ_CACHE` (off) | oracle | `crates/dft-cli/src/run.rs:4041`; `crates/dft-solvers/src/davidson_stages.rs:271` | Stage-0 Davidson screens (FASTBENCH-DAVIDSON, 2026-11-16; Mac observations O, counters exact; `scripts/fast_screen_dav.py`, one solve at a fixed Hamiltonian, si32-gamma f1e-3 seed 601 T=1, base density = converged `dft-cli run` density.bin, cold = atomic start at 5e-3 Ha, warm = kept start at the ratchet tolerance 3.6e-11 Ha, 3 alternating pairs x 5 repeats, pure build, load 3 to 5; protocol docs/performance/fast-loop-protocol.md stage 0) (line 1001) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `DAVIDSON_PACKED_REDUCE` (off) | oracle | `crates/dft-solvers/src/distributed_davidson/workspace/lean.rs:350` | DAVIDSON-PACKED-REDUCE — alternative pending dev08 (2026-11-20) (line 2687) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `DAVIDSON_ROTATE_COLUMNS` (off) | oracle | `crates/dft-solvers/src/real_davidson.rs:1894` | DAVIDSON-UPPER-BORDER-REJECTED — orchestrator decision (2026-11-20) (line 2636) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `DAVIDSON_UPPER_BORDER` (off) | oracle | `crates/dft-solvers/src/real_davidson.rs:1692` | LANE-DAVIDSON Phase 2 — upper-border GEMM tiles (2026-11-20) (line 2557) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `DAVIDSON_PROFILE` (off) | diagnostic | `crates/dft-solvers/src/real_davidson.rs:3355` | LANE-DAVIDSON Phase 0 — requested 2026-11-20 (line 2419) | Keep: R11(e). |
| `LEAN_BORDER` (separate) | oracle | `crates/dft-solvers/src/davidson_stages.rs:274` | MPI overlap-matrix Davidson (MPI-OVERLAP-DAVIDSON, 2026-11-17; Mac observations O on simulated ranks, counters exact C, estimates E; design docs/design/mpi-overlap-davidson.md, dev08 request docs/performance/requests/mpi-overlap-davidson.md; review F05 step 2) (line 1069) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `RESTART_KEEP` (off) | oracle | `crates/dft-cli/src/run.rs:4046`; `crates/dft-solvers/src/davidson_stages.rs:182` | Breakthrough restart (BREAKTHROUGH-1-RESTART, 2026-11-17; review 2026-11-17 R1 / E1 / E2, backlog B05; Mac observations O on simulated ranks and T=1, counters exact C; design docs/design/restart-keep.md, dev08 request docs/performance/requests/breakthrough-restart.md) (line 1077) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `MKL_NONLOCAL` (on) | algorithm module | `crates/dft-physics/src/real_gamma.rs:216` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 875) | Keep: Generality rule/backend selection. |
| `APPLY_LEAN` (on) | oracle | `crates/dft-physics/src/real_gamma.rs:614` | Algebra stream (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/algebra-stream.md) (line 967) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `HAPPLY_OVERLAP` (off) | oracle | `crates/dft-physics/src/real_gamma.rs:616` | H application (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/h-apply.md) (line 852) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `HAPPLY_FOLD` (off) | oracle | `crates/dft-physics/src/real_gamma.rs:622` | H application (2026-11-11; implementer, Mac observations; dev08 request docs/performance/requests/h-apply.md) (line 853) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `NONLOCAL_BLOCK` (auto) | oracle | `crates/dft-physics/src/hamiltonian.rs:704,2655` | Stage-0 Davidson screens (FASTBENCH-DAVIDSON, 2026-11-16; Mac observations O, counters exact; `scripts/fast_screen_dav.py`, one solve at a fixed Hamiltonian, si32-gamma f1e-3 seed 601 T=1, base density = converged `dft-cli run` density.bin, cold = atomic start at 5e-3 Ha, warm = kept start at the ratchet tolerance 3.6e-11 Ha, 3 alternating pairs x 5 repeats, pure build, load 3 to 5; protocol docs/performance/fast-loop-protocol.md stage 0) (line 1006) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `ROTATE_KERNEL` (legacy) | oracle | `crates/dft-numerics/src/real.rs:265` | P47B — fp32 H feasibility (2026-11-08; measurement only, Mac observations) (line 759) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `PURE_GEMM` (nano, packed) | oracle | `crates/dft-numerics/src/real.rs:269` | Stage-0 Davidson screens (FASTBENCH-DAVIDSON, 2026-11-16; Mac observations O, counters exact; `scripts/fast_screen_dav.py`, one solve at a fixed Hamiltonian, si32-gamma f1e-3 seed 601 T=1, base density = converged `dft-cli run` density.bin, cold = atomic start at 5e-3 Ha, warm = kept start at the ratchet tolerance 3.6e-11 Ha, 3 alternating pairs x 5 repeats, pure build, load 3 to 5; protocol docs/performance/fast-loop-protocol.md stage 0) (line 1004) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `PACK_FIXED` (off) | oracle | `crates/dft-numerics/src/real.rs:1352` | Per-core kernels (PERCORE-KERNELS, 2026-11-17; Mac observations O, pure build, counters exact C, relayed dev08 R; analysis docs/performance/percore-kernels-2026-11-17.md, dev08 request docs/performance/requests/percore-kernels.md; stage 0 = `scripts/fast_screen_dav.py`, si32/si64-gamma f1e-3 seed 601 T=1, converged base density, 3 alternating pairs x 5 repeats; load 2 to 21 during the session, spreads as printed) (line 1057) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `MKL_DGEMM` (on) | algorithm module | `crates/dft-numerics/src/real.rs:274` | P47B — fp32 H feasibility (2026-11-08; measurement only, Mac observations) (line 765) | Keep: Generality rule/backend selection. |
| `MKL_EVD` (on) | algorithm module | `crates/dft-numerics/src/real.rs:2273` | P47B — fp32 H feasibility (2026-11-08; measurement only, Mac observations) (line 782) | Keep: Generality rule/backend selection. |
| `PAIR_FFT` (grid, pruned) | algorithm module | `crates/dft-numerics/src/fused_pair_plan.rs:349` | F-LAUNCH-IO (2026-11-19; base ff0c83caaccdfc26b83598cc7b073d980f055d17; Mac observations, dev08 pending) (line 1279) | Keep: Generality rule/backend selection. |
| `SPIN_SHARE` (on) | oracle | `crates/dft-engine/src/collinear.rs:255` | ONEPATH-5-SPIN-SHARE (2026-11-17) (line 1133) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `BINARY_HASH_CACHE` (on) | infra | `crates/dft-cli/src/run.rs:1231` | LAUNCHIO-BINARY-HASH-CACHE-DEFAULT-ON (2026-11-20, dev08 via master; ledger by orchestrator) (line 3027) | Keep: routing/infrastructure. |
| `RELAX_EXTRAPOLATION` (config) | infra | `crates/dft-cli/src/relax.rs:52` | P47B — fp32 H feasibility (2026-11-08; measurement only, Mac observations) (line 815) | Keep: routing/infrastructure. |
| `REFERENCE_SCF` (on) | infra | `crates/dft-cli/src/driver_select.rs:156,162` | No named row; manifest semantics / reader scope | Keep: routing/infrastructure. |
| `REFERENCE_STRICT` (off) | infra | `crates/dft-cli/src/driver_select.rs:137,141`; `crates/dft-engine/src/model/options.rs:203` | R10a final verification (line 3892) | Keep: routing/infrastructure. |
| `DGEMM_TIMERS` (see manifest) | diagnostic | `crates/dft-cli/src/main.rs:374` | No named row; manifest semantics / reader scope | Keep: R11(e). |
| `PSEUDO_DIR` (see manifest) | infra | `crates/dft-cli/src/pseudo.rs:175` | No named row; manifest semantics / reader scope | Keep: routing/infrastructure. |
| `FFT_KEPT_FUSE` (off) | oracle | `crates/dft-solvers/src/distributed_fft.rs:684` | FFT-ARRIVAL-OWNERSHIP-CLOSED (2026-11-20) (line 3024) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `COLLECTIVE_TRACE` (off) | diagnostic | `crates/dft-numerics/src/stage_clock/arrival.rs:44` | LANE-FFT arrival trace implementation (2026-11-20) (line 2825) | Keep: R11(e). |
| `FFT_PROFILE` (off) | diagnostic | `crates/dft-numerics/src/real.rs:2453`; `crates/dft-numerics/src/stage_clock.rs:152`; `crates/dft-numerics/src/fused_pair_plan.rs:71`; `crates/dft-cli/src/run.rs:3555`; `crates/dft-engine/src/profile.rs:966` | LANE-FFT Phase 0 (2026-11-20; revised user instruction) (line 2438) | Keep: R11(e). |
| `SERIAL_DENSITY_FUSED` (on) | oracle | `crates/dft-engine/src/band_density.rs:102` | GLUE-SERIAL-DENSITY-FUSED-DEFAULT-ON (2026-11-20, dev08 via master; ledger by orchestrator) (line 3030) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `GLUE_PROFILE` (off) | diagnostic | No direct production reader (manifest/test-only or legacy removal) | LANE-GLUE serial density fusion — 2026-11-20, default-off, dev08 request (line 2905) | Keep: R11(e). |
| `PHASE_PROFILE` (off) | diagnostic | `crates/dft-numerics/src/real.rs:2453`; `crates/dft-numerics/src/fused_pair_plan.rs:72`; `crates/dft-physics/src/distributed_real_gamma.rs:1714`; `crates/dft-physics/src/hamiltonian.rs:366`; `crates/dft-cli/src/run_config.rs:162,163`; `crates/dft-engine/src/profile.rs:616` | LANE-GLUE serial density fusion — 2026-11-20, default-off, dev08 request (line 2905) | Keep: R11(e). |
| `M2_LIGHT_STATE` (on) | oracle | `crates/dft-engine/src/gamma_real.rs:293` | ONEPATH-5-SPIN-SHARE (2026-11-17) (line 1127) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_OWNED_INPUT` (off) | oracle | `crates/dft-engine/src/gamma_real.rs:1104` | R-SCF-START — final-cutoff owned initialization (2026-11-23) (line 3176) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_LAZY_SERIAL` (on) | oracle | `crates/dft-engine/src/gamma_real.rs:315` | INIT-SCALE: making `init` scale with MPI ranks (2026-11-18; Mac observations O on 8 simulated ranks, load 3 to 5, pure build; counters exact C; branch base main 0c4d652, commits 03fd272, 5262fcf, a7205c9, 529e485) (line 1163) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_LAYOUT_CACHE` (off) | oracle | `crates/dft-solvers/src/distributed_davidson.rs:120` | ONEPATH-5-SPIN-SHARE (2026-11-17) (line 1127) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `EWALD_CACHE` (on) | oracle | `crates/dft-cli/src/run.rs:801,804`; `crates/dft-engine/src/scf.rs:419` | INIT-SCALE: making `init` scale with MPI ranks (2026-11-18; Mac observations O on 8 simulated ranks, load 3 to 5, pure build; counters exact C; branch base main 0c4d652, commits 03fd272, 5262fcf, a7205c9, 529e485) (line 1161) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `EWALD_FAST` (on) | oracle | `crates/dft-physics/src/ewald.rs:81` | Ewald and ionic local construction (INIT-EWALD-FAST, 2026-11-18; Mac observations O, release, one thread, examples `ewald_bench` / `ionic_bench`; counts exact C; dev08 request docs/performance/requests/init-ewald-fast.md; base main 9d942a0) (line 1210) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `IONIC_FAST` (on) | oracle | `crates/dft-physics/src/hamiltonian.rs:79` | Ewald and ionic local construction (INIT-EWALD-FAST, 2026-11-18; Mac observations O, release, one thread, examples `ewald_bench` / `ionic_bench`; counts exact C; dev08 request docs/performance/requests/init-ewald-fast.md; base main 9d942a0) (line 1209) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `INIT_PROFILE` (off) | diagnostic | No direct production reader (manifest/test-only or legacy removal) | LANE-INIT Phase 2 — local atomic block (2026-11-20) (line 2536) | Keep: R11(e). |
| `PHASE_VALIDATE` (off) | diagnostic | `crates/dft-physics/src/hamiltonian.rs:367,3471` | G-PHASE-VALIDATE (2026-11-19): exact finite-dot range proof (line 2303) | Keep: R11(e). |
| `IMPORT_CHECK` (split) | oracle | `crates/dft-engine/src/session.rs:2657` | Between, verification glue and materialize (BETWEEN-MATERIALIZE, 2026-11-18; Mac observations O on 8 simulated ranks, pure build, loaded Mac: buckets only, walls not comparable; counters exact C; base main e902853, commits 5fe2796, b3f8509, 3f7d82b) (line 1186) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `M2_LAZY_RESULT` (on) | oracle | `crates/dft-engine/src/refscf/lazy.rs:27` | Between, verification glue and materialize (BETWEEN-MATERIALIZE, 2026-11-18; Mac observations O on 8 simulated ranks, pure build, loaded Mac: buckets only, walls not comparable; counters exact C; base main e902853, commits 5fe2796, b3f8509, 3f7d82b) (line 1188) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `PHASE_WAIT` (off) | diagnostic | `crates/dft-solvers/src/distributed_davidson.rs:144` | Other and init (BREAKTHROUGH-4-OTHER-AND-INIT, 2026-11-17; Mac observations O on simulated ranks np=8/4/2 and T=1, counts exact C; note docs/performance/profile-2026-11-17-other-init.md, dev08 request docs/performance/requests/breakthrough-other-init.md). Starting point (relayed dev08, a220a85, si32 np=8): fine `solve.other` 0.2601 s = 24.7 % of the fine solve, remainder 0.2552 s, init 0.5606 s. (line 1108) | Keep: R11(e). |
| `SNAPSHOT_DIR` (see manifest) | diagnostic | `crates/dft-engine/src/snapshot.rs:62` | Measurement infrastructure (BREAKTHROUGH-0-INSTRUMENT, 2026-11-17; review 2026-11-17 sections 13.1 to 13.3 and E3; Mac observations O on simulated ranks np=4 and T=1, counts exact C; notes docs/performance/profile-2026-11-17-phases.md, docs/performance/first-block-residual-2026-11-17.md; dev08 request docs/performance/requests/breakthrough-instrument.md). No speed claim; three default-off switches whose on/off state leaves energy bits and density bits unchanged (main d5d04b1 harness against this branch, si8/si32 np=4 production + ladder; dft-cli si32 T=1 density SHA 9491950d...; test `dft-engine/tests/instrument.rs`). (line 1094) | Keep: R11(e). |
| `REFSCF_SETUP_PHASE_REUSE` (on) | rollback default-on | `crates/dft-engine/src/refstart.rs:830,1060,1491` | P3 timed (line 3539) | Remove: D3; preserve default arithmetic. |
| `REF_DAV_SINGLE_GEMM` (on) | rollback default-on | `crates/dft-solvers/src/ref_rmm.rs:314`; `crates/dft-solvers/src/ref_davidson.rs:134`; `crates/dft-solvers/src/ref_davidson_complex.rs:245` | REF_DAV_SINGLE_GEMM provisionally ADOPTED (2026-11-25, development protocol) (line 3578) | Remove: D3; preserve default arithmetic. |
| `FFT_SPARSE_CLEAR` (off) | oracle | `crates/dft-solvers/src/distributed_fft.rs:634` | P3c timed (2026-11-25): REF_DAV_SINGLE_GEMM and FFT_DIRECT_PACK/FFT_SPARSE_CLEAR not adopted (inconclusive) (line 3570) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `FFT_DIRECT_PACK` (off) | oracle | `crates/dft-physics/src/distributed_real_gamma.rs:424` | P3c timed (2026-11-25): REF_DAV_SINGLE_GEMM and FFT_DIRECT_PACK/FFT_SPARSE_CLEAR not adopted (inconclusive) (line 3570) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `REFSCF_VERIFY_ENERGY` (off) | experiment default-off | `crates/dft-cli/src/run.rs:3307,3313`; `crates/dft-engine/src/refscf/driver.rs:267,274`; `crates/dft-engine/src/refscf/electronic.rs:331`; `crates/dft-engine/src/refscf/k_driver.rs:147,154` | P3f/P3g screens (2026-11-25, development protocol, one run each, dev08 MPI8) (line 3588) | Remove: D3; preserve default arithmetic. |
| `FFT_FUSED_ROWS` (off) | oracle | `crates/dft-solvers/src/distributed_fft.rs:645,648,654,1651` | P3f/P3g screens (2026-11-25, development protocol, one run each, dev08 MPI8) (line 3591) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `REF_DAV_SWAP_BUFFERS` (on) | rollback default-on | `crates/dft-solvers/src/ref_davidson_complex.rs:238`; `crates/dft-solvers/src/ref_davidson/workspace.rs:177` | P3i (line 3598) | Remove: D3; preserve default arithmetic. |
| `REF_DAV_LEAN_RESET` (on) | rollback default-on | `crates/dft-solvers/src/ref_davidson_complex.rs:251`; `crates/dft-solvers/src/ref_davidson/workspace.rs:175` | P3i (line 3598) | Remove: D3; preserve default arithmetic. |
| `REFSTART_DIAGONAL_REUSE` (on) | rollback default-on | `crates/dft-engine/src/refstart.rs:563` | P3i (line 3598) | Remove: D3; preserve default arithmetic. |
| `REFSTART_ORBITAL_REUSE` (on) | rollback default-on | `crates/dft-engine/src/refstart.rs:1493` | P3i (line 3598) | Remove: D3; preserve default arithmetic. |
| `REF_DAV_SQRT` (on) | rollback default-on | `crates/dft-solvers/src/ref_davidson_complex.rs:236`; `crates/dft-solvers/src/ref_davidson/workspace.rs:159` | P3j (line 3611) | Remove: D3; preserve default arithmetic. |
| `REF_DAV_FUSED_APPEND` (on) | rollback default-on | `crates/dft-solvers/src/ref_davidson_complex.rs:248`; `crates/dft-solvers/src/ref_davidson/workspace.rs:161` | P3j (line 3611) | Remove: D3; preserve default arithmetic. |
| `REF_DAV_CHECK_APPEND` (on) | rollback default-on | `crates/dft-solvers/src/ref_davidson_complex.rs:249`; `crates/dft-solvers/src/ref_davidson/workspace.rs:162` | P3j (line 3611) | Remove: D3; preserve default arithmetic. |
| `REF_DAV_REUSE_ROWS` (on) | rollback default-on | `crates/dft-solvers/src/ref_davidson_complex.rs:252`; `crates/dft-solvers/src/ref_davidson/workspace.rs:94` | P3j (line 3611) | Remove: D3; preserve default arithmetic. |
| `FFT_A2A_BALANCE` (off) | oracle | `crates/dft-solvers/src/distributed_fft.rs:276,359,1418` | No named row; manifest semantics / reader scope | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `REFSTART_GENERALIZED` (on) | rollback default-on | `crates/dft-engine/src/refstart.rs:1606` | P3l (line 3620) | Remove: D3; preserve default arithmetic. |
| `REFSTART_COMPACT_G` (on) | rollback default-on | `crates/dft-engine/src/refstart.rs:487` | P3l (line 3620) | Remove: D3; preserve default arithmetic. |
| `NONLOCAL_SINGLE_GEMM` (on) | oracle | `crates/dft-physics/src/distributed_real_gamma.rs:445`; `crates/dft-cli/src/run.rs:4145` | P3l (2026-11-25, development protocol): REFSTART_GENERALIZED, REFSTART_COMPACT_G, NONLOCAL_SINGLE_GEMM, REFSCF_LAZY_OUTPUT ADOPTED (line 3620) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `REFSCF_LAZY_OUTPUT` (on) | rollback default-on | `crates/dft-engine/src/refscf/driver.rs:1091` | P3l (line 3620) | Remove: D3; preserve default arithmetic. |
| `REF_DAV_IDEAL_RESET` (on) | rollback default-on | `crates/dft-solvers/src/ref_davidson_complex.rs:239`; `crates/dft-solvers/src/ref_davidson/workspace.rs:172`; `crates/dft-solvers/src/ref_davidson/ideal_reset.rs:46,59,184,246` | P3m (line 3628) | Remove: D3; preserve default arithmetic. |
| `NONLOCAL_FUSE_COUPLING` (on) | oracle | `crates/dft-physics/src/distributed_real_gamma.rs:434,448` | P3m (2026-11-25, development protocol): REF_DAV_IDEAL_RESET, NONLOCAL_FUSE_COUPLING, NONLOCAL_BACK_ACCUM ADOPTED (line 3628) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `NONLOCAL_BACK_ACCUM` (on) | oracle | `crates/dft-physics/src/distributed_real_gamma.rs:450` | P3m (2026-11-25, development protocol): REF_DAV_IDEAL_RESET, NONLOCAL_FUSE_COUPLING, NONLOCAL_BACK_ACCUM ADOPTED (line 3628) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `REFSCF_LEAN_ENERGY` (on) | rollback default-on | `crates/dft-engine/src/refscf/driver.rs:423`; `crates/dft-engine/src/refscf/electronic.rs:341`; `crates/dft-engine/src/refscf/k_driver.rs:297` | P3q (line 3638) | Remove: D3; preserve default arithmetic. |
| `REFSCF_RMM` (off) | algorithm module | `crates/dft-cli/src/run.rs:4224`; `crates/dft-engine/src/refscf/driver.rs:257`; `crates/dft-engine/src/refscf/electronic.rs:329`; `crates/dft-engine/src/refscf/k_driver.rs:138` | P3r (2026-11-26, development protocol): REFSCF_RMM (QE-style rmm-davidson) not adopted, stays off (line 3647) | Keep: Generality rule/backend selection. |
| `REFSCF_RMM_STEPS` (4) | algorithm module | `crates/dft-engine/src/refscf/driver.rs:259` | No named row; manifest semantics / reader scope | Keep: Generality rule/backend selection. |
| `REFSCF_RMM_BLOCK` (16) | algorithm module | `crates/dft-engine/src/refscf/driver.rs:260` | No named row; manifest semantics / reader scope | Keep: Generality rule/backend selection. |
| `REFSTART_FFT_SHARE` (on) | rollback default-on | `crates/dft-engine/src/refstart.rs:388` | P3t (line 3659) | Remove: D3; preserve default arithmetic. |
| `REFSTART_SINGLE_GEMM` (on) | rollback default-on | `crates/dft-engine/src/refstart.rs:1149` | P3t (line 3659) | Remove: D3; preserve default arithmetic. |
| `REFSTART_LOCAL_SHELLS` (off) | experiment default-off | `crates/dft-engine/src/refstart.rs:985` | P3t (line 3659) | Remove: D3; preserve default arithmetic. |
| `FFT_KERNEL_BY_LENGTH` (off) | oracle | `crates/dft-numerics/src/fused_pair_plan.rs:235`; `crates/dft-cli/src/run.rs:2435,2439`; `crates/dft-solvers/src/distributed_fft.rs:172,1589` | Work-order screen (2026-11-26, user work order "reduce full-cost SCF work"), dev08 Si64 MPI8, one run per arm at bb416c73 (line 3678) | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `REFSCF_PERMAP` (off) | experiment default-off | `crates/dft-engine/src/refscf/driver.rs:424`; `crates/dft-engine/src/refscf/k_driver.rs:298`; `crates/dft-engine/src/refscf/broyden.rs:148` | No named row; manifest semantics / reader scope | Keep: retained option or no rejection/supersession evidence. |
| `REFSCF_LADDER_RATIO` (100) | algorithm module | `crates/dft-cli/src/run_config.rs:161`; `crates/dft-engine/src/refscf/driver.rs:193` | Ladder retry with looser coarse stop (2026-11-26, dev08 MPI8, one run per arm at ac6aee9c): rejected (line 3703) | Keep: Generality rule/backend selection. |
| `REFSCF_LADDER` (off) | algorithm module | `crates/dft-engine/src/refscf/driver.rs:184,437,468`; `crates/dft-engine/src/refscf/electronic.rs:330`; `crates/dft-engine/src/refscf/k_driver.rs:138` | Ladder retry with looser coarse stop (2026-11-26, dev08 MPI8, one run per arm at ac6aee9c): rejected (line 3703) | Keep: Generality rule/backend selection. |
| `REFSCF_KEPT_HX` (on) | rollback default-on | `crates/dft-engine/src/refscf/driver.rs:393,427`; `crates/dft-engine/src/refscf/electronic.rs:485`; `crates/dft-engine/src/refscf/k_driver.rs:301` | KEPTHX2 (2026-11-26, development protocol): REFSCF_KEPT_HX ADOPTED (line 3682) | Remove: D3; preserve default arithmetic. |
| `REFSCF_VERIFY_WARM` (on) | rollback default-on | `crates/dft-engine/src/refscf/driver.rs:426`; `crates/dft-engine/src/refscf/electronic.rs:684`; `crates/dft-engine/src/refscf/k_driver.rs:300` | P3-VERWARM (2026-11-26, development protocol): REFSCF_VERIFY_WARM ADOPTED (line 3711) | Remove: D3; preserve default arithmetic. |
| `REFSTART_CHOLESKY` (on) | rollback default-on | `crates/dft-engine/src/refstart.rs:1610` | P4 screen (line 3724) | Remove: D3; preserve default arithmetic. |
| `REFSTART_SPECIES_FILL` (on) | rollback default-on | `crates/dft-engine/src/refstart.rs:1495` | P4 screen (line 3724) | Remove: D3; preserve default arithmetic. |
| `REFSTART_PARTITION_CACHE` (on) | rollback default-on | `crates/dft-engine/src/refstart.rs:390` | P4 screen (line 3724) | Remove: D3; preserve default arithmetic. |
| `REFSCF_KEPT_RECOMPUTE` (off) | experiment default-off | `crates/dft-engine/src/refscf/kept_hx.rs:73` | P4 screen (2026-11-26, dev08 MPI8 verify off, one run per arm at db29dd7) (line 3728) | Keep: retained option or no rejection/supersession evidence. |
| `REF_DAV_TRIANGLE` (on) | rollback default-on | `crates/dft-solvers/src/ref_davidson_complex.rs:247`; `crates/dft-solvers/src/ref_davidson/workspace.rs:156` | GEN-2432530 (line 3732) | Remove: D3; preserve default arithmetic. |
| `REF_DAV_STACK_ROTATE` (off) | experiment default-off | `crates/dft-solvers/src/ref_davidson_complex.rs:250`; `crates/dft-solvers/src/ref_davidson/workspace.rs:146` | GEN-2432530 (line 3732) | Remove: D3; preserve default arithmetic. |
| `REF_DAV_SETUP_FUSED` (off) | experiment default-off | `crates/dft-solvers/src/ref_davidson.rs:133`; `crates/dft-solvers/src/ref_davidson_complex.rs:244` | GEN-2432530 (line 3732) | Remove: D3; preserve default arithmetic. |
| `REFSCF_MEMORY_TRIM` (off) | experiment default-off | `crates/dft-engine/src/refscf/driver.rs:256,1373`; `crates/dft-engine/src/refscf/k_driver.rs:136`; `crates/dft-engine/src/refscf/map.rs:291,321`; `crates/dft-solvers/src/ref_davidson/workspace.rs:95` | No named row; manifest semantics / reader scope | Keep: retained option or no rejection/supersession evidence. |
| `FORCE_BAND_FIRST` (off) | oracle | `crates/dft-engine/src/forces.rs:515` | No named row; manifest semantics / reader scope | Keep for R11b: legacy or oracle-shared kernel/configuration. |
| `REFSCF_SPIN_FD_MOMENT` (on) | infra | `crates/dft-cli/src/driver_select.rs:303` | No named row; manifest semantics / reader scope | Keep: CLI routing reader outside ownership. |
| `REFSCF_SPIN_STIFFNESS` (on) | rollback default-on | `crates/dft-engine/src/refscf/electronic.rs:480` | R6a / G4B 48/48, dev08 gen15; ADR 0045 (manifest adoption record; no separate ledger row) | Remove: D3; preserve default arithmetic. |
| `REFSCF_HISTORY_CONDITIONING` (on) | rollback default-on | `crates/dft-engine/src/refscf/electronic.rs:481,525` | R6a / G4B 48/48, dev08 gen15; ADR 0045 (manifest adoption record; no separate ledger row) | Remove: D3; preserve default arithmetic. |
| `XC_SCREEN_NEGATIVE` (on) | oracle | `crates/dft-physics/src/xc.rs:120` | No named row; manifest semantics / reader scope | Keep: ADR 0046 comparison off-switch. |

## Test disposition

All remaining identity tests keep their numerical tolerances. Retired input/environment pins are removed; cold/warm repetitions, thread and simulated-rank comparisons, independent dense/scalar oracles, capacity fallbacks and cache invalidation remain. Work-count assertions for deleted off arms now check the adopted sqrt/append/cache path. The restart-correction identity test is retained, including its assertion that post-contraction identity GEMMs are absent.

The following 11 test functions exercised retired experiments and are deleted (helpers are not counted):

| File | Deleted test | Reason |
| --- | --- | --- |
| `crates/dft-cli/tests/checkpoint_v2.rs` | `verification_reuse_identity` | Retired failed-verification reuse / predictor A/B matrix |
| `crates/dft-cli/tests/checkpoint_v2.rs` | `pending_verification_file_roundtrip` | Retired failed-verification reuse / predictor A/B matrix |
| `crates/dft-engine/tests/refscf.rs` | `first_prediction_linear_oracle_and_read_only_history` | Retired first-attempt predictor |
| `crates/dft-engine/tests/refscf.rs` | `verification_cost_small_fixture` | Retired failed-verification reuse / predictor A/B matrix |
| `crates/dft-engine/tests/refscf.rs` | `first_prediction_eviction_and_preconditioner_independence` | Retired first-attempt predictor |
| `crates/dft-engine/tests/refscf.rs` | `verification_cost_failed_attempt_fixture` | Retired failed-verification reuse / predictor A/B matrix |
| `crates/dft-engine/tests/refscf.rs` | `verification_energy_attempts_and_failure_recovery` | Retired energy acceptance |
| `crates/dft-engine/tests/refscf.rs` | `verification_energy_threshold_is_strict_and_finite` | Retired energy acceptance |
| `crates/dft-engine/tests/refscf_checkpoint.rs` | `verification_reuse_identity` | Retired failed-verification reuse / predictor A/B matrix |
| `crates/dft-engine/tests/refscf_checkpoint.rs` | `pending_verification_map_is_restored_without_replay` | Retired failed-verification reuse / predictor A/B matrix |
| `crates/dft-solvers/src/ref_davidson/remainder_tests.rs` | `p4_stacked_rotation_scalar_oracle_traffic_and_x_only_return` | Rejected stacked rotation |

The generic verification-cost contract and ordinary checkpoint identity tests remain. Historical checkpoint payload layouts remain readable; nonempty failed-map reuse state now fails closed.

## Final validation (2026-11-25)

All counts compare source commit `6537cb6fa063a2b5ee382e1cbf72a92a2949feb1` with this uncommitted worktree. Rust lines and test annotations cover every tracked `crates/**/*.rs` file (including examples); runtime test totals are reported separately.

| Measure | Before | After | Change |
| --- | ---: | ---: | ---: |
| Manifest switches | 168 | 133 | -35 (29 default-on, 6 default-off) |
| Rust source lines | 210370 | 208157 | -2213 |
| `#[test]` annotations | 1258 | 1247 | -11 |

Final commands below all exited **0**. Logs are in `/private/tmp/r11a/`; these are local validation logs, not portable reference artifacts.

```sh
cargo check --workspace --tests
cargo check -p dft-solvers -p dft-numerics -p dft-mkl --tests --features dft-numerics/mkl
CHECK_SCCACHE=0 CHECK_NICE=0 BAANDI_TEST_TMP=/private/tmp/r11a/tests TMPDIR_990=/private/tmp/r11a/check scripts/check.sh
cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
cargo test -p dft-cli --test reference_identity -- --nocapture
cargo test -p dft-solvers --lib ref_davidson -- --nocapture
cargo build --release -p dft-cli
python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /Volumes/990-2TB/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --threads 1 --repeats 1 --baseline docs/performance/baseline-reference.json --require-identical --record /private/tmp/r11a/reference-final.json
BAANDI_REFERENCE_STRICT=on python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /Volumes/990-2TB/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --threads 1 --repeats 1 --baseline docs/performance/baseline-oracle.json --require-identical --record /private/tmp/r11a/oracle-final.json
python3 -I scripts/make_readme_figures.py
python3 -I scripts/tests/test_readme_figures.py
git diff --check
```

Full gate: **1010 passed, 0 failed, 228 ignored, 112 binaries**; build, test, fmt and clippy PASS (`full-final.log`). Reference Davidson: **22 passed**, including real/complex solves, supplied HX, restart corrections, guarded resets, independent dense/scalar checks, fused append and counted row traffic (`solver-final.log`). The full gate also runs kept-HX delta/invalidation, delivery fingerprint, complex-k kept-state and warm-verification recovery tests. `reference_identity` passed its unchanged golden checks. README regeneration has no diff; all nine freshness tests pass.

The first full gate failed on obsolete test pins/counts and off-arm work assertions, plus one unused test binding. Those were corrected without changing scientific tolerances or golden values, and the full gate was rerun successfully. The optional load probe cannot execute `ps` in this sandbox; its warning does not suppress any of the four successful gates. MKL was compile-checked only; real MPI/MKL runtime and dev08 performance measurements are outside this local validation.

### Stored baseline identity

Both files contain exactly the three cases below. The final release binary matches **energy bits, density SHA-256, iteration counts and recorded solver work** in the final recorded JSON. A separate Python assertion compared all seven fields (`total_energy_bits_hex`, `density_sha256`, `iterations`, `operator_applications`, `davidson_iterations`, `certificate_applications`, `certificate_iterations`) for all six path/case pairs and exited 0. The benchmark itself checks density identity and energy tolerance (`reference-final.log`, `oracle-final.log`). Delta energy is **0 Ha** in every case. These are BAANDI regression baselines, not new QE evidence or performance results. UPF SHA-256: `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.

| Path | Case | Energy (Ha) | Iterations | Density |
| --- | --- | ---: | ---: | --- |
| reference | gamma | -7.295530466374794 | 7 | bit-identical |
| reference | k222s | -7.919423412321121 | 7 | bit-identical |
| reference | k222 | -7.831230055901806 | 7 | bit-identical |
| strict oracle | gamma | -7.295530468210924 | 8 | bit-identical |
| strict oracle | k222s | -7.91942341246909 | 10 | bit-identical |
| strict oracle | k222 | -7.831230057121148 | 10 | bit-identical |

Baseline JSON, identity goldens, legacy drivers and README figures remain unchanged. No timing claim is inferred from concurrent local builds or baseline checks. The complete changes are left uncommitted for review.

## R11b: freeze oracle-only switches (D3)

Source/recreation commit: `1535810f4d6abffbe36483dd82c174ece6664b61`. Reader locations below are pre-change.
The manifest shrinks from **121 to 75 switches** (46 retired). Each candidate was
searched with `git grep Id::NAME`, followed by bare-name/string searches and a
caller audit. The serial real Gamma operator, legacy real/distributed Davidson
solvers and Pulay mixer are used by the oracle; reference Davidson/RMM share
traits and selected factorization/panel helpers, which remain. No reference
arithmetic, baseline, golden, or tolerance is changed. This is retirement,
not a performance adoption or a QE comparison.

| Retired switch | Frozen default | Readers at source commit | Decision |
| --- | --- | --- | --- |
| `INIT_OWNED_STATES` | on | `crates/dft-engine/src/oracle/scf.rs:775,801` | Freeze default; delete switch and alternative selection. |
| `CERT_SCREEN` | on | `crates/dft-engine/src/oracle/session.rs:1565` | Freeze default; delete switch and alternative selection. |
| `POTENTIAL_REUSE` | on | `crates/dft-engine/src/oracle/session.rs:1570` | Freeze default; delete switch and alternative selection. |
| `MIX_GROWTH_RESET` | None | `crates/dft-engine/src/oracle/session.rs:1571` | Freeze default; delete switch and alternative selection. |
| `MIX_CELL_KERKER` | off | `crates/dft-engine/src/oracle/scf.rs:689,700,706` | Freeze default; delete switch and alternative selection. |
| `EMPTY_LOOSE` | on | `crates/dft-engine/src/oracle/scf.rs:2835` | Freeze default; delete switch and alternative selection. |
| `BUFFER_TOL` | 1e-3 | `crates/dft-engine/src/oracle/scf.rs:2845` | Freeze default; delete switch and alternative selection. |
| `DENSITY_PAIR` | on | `crates/dft-engine/src/oracle/band_density.rs:17` | Freeze default; delete switch and alternative selection. |
| `DENSITY_REDUCE` | on | `crates/dft-engine/src/oracle/band_density.rs:26` | Freeze default; delete switch and alternative selection. |
| `NL_ENERGY_REAL` | on | `crates/dft-engine/src/oracle/gamma_real.rs:144` | Freeze default; delete switch and alternative selection. |
| `ATOMIC_LOCAL_BLOCK` | on | `crates/dft-engine/src/oracle/atomic_start.rs:158` | Freeze default; delete switch and alternative selection. |
| `ATOMIC_REAL` | on | `crates/dft-engine/src/oracle/gamma_real.rs:373` | Freeze default; delete switch and alternative selection. |
| `PHASE_THREADS` | on | `crates/dft-engine/src/oracle/collinear.rs:285`; `crates/dft-engine/src/oracle/scf.rs:1815` | Freeze default; delete switch and alternative selection. |
| `MIX_COULOMB_GUARD` | off | `crates/dft-engine/src/oracle/session.rs:310,364` | Freeze default; delete switch and alternative selection. |
| `M2_LIGHT_STATE` | on | `crates/dft-engine/src/oracle/gamma_real.rs:90` | Freeze default; delete switch and alternative selection. |
| `M2_LAZY_SERIAL` | on | `crates/dft-engine/src/oracle/gamma_real.rs:112` | Freeze default; delete switch and alternative selection. |
| `IMPORT_CHECK` | split | `crates/dft-engine/src/oracle/session.rs:1564` | Freeze default; delete switch and alternative selection. |
| `M2_DENSITY` | fused | No production reader (orphan after R10). | Freeze default; delete switch and alternative selection. |
| `M2_ORTHO` | serial | `crates/dft-solvers/src/distributed_davidson.rs:813` | Freeze default; delete switch and alternative selection. |
| `M2_LEAN` | on | `crates/dft-solvers/src/distributed_davidson.rs:816` | Freeze default; delete switch and alternative selection. |
| `M2_EIGEN` | replicated | `crates/dft-solvers/src/distributed_davidson.rs:817` | Freeze default; delete switch and alternative selection. |
| `M2_B5` | off | `crates/dft-solvers/src/distributed_davidson.rs:818` | Freeze default; delete switch and alternative selection. |
| `M2_WORKSPACE` | on | `crates/dft-solvers/src/distributed_davidson.rs:819` | Freeze default; delete switch and alternative selection. |
| `M2_SELECTED_ROOTS` | on | `crates/dft-solvers/src/distributed_davidson.rs:820` | Freeze default; delete switch and alternative selection. |
| `M2_PANEL_GS` | fused | `crates/dft-solvers/src/distributed_davidson.rs:821` | Freeze default; delete switch and alternative selection. |
| `M2_START` | gram | `crates/dft-solvers/src/distributed_davidson.rs:873` | Freeze default; delete switch and alternative selection. |
| `M2_LEAN_ALGEBRA` | on | `crates/dft-engine/src/oracle/snapshot.rs:389`; `crates/dft-solvers/src/distributed_davidson.rs:826` | Freeze default; delete switch and alternative selection. |
| `M2_LAYOUT_CACHE` | off | `crates/dft-solvers/src/distributed_davidson.rs:120` | Freeze default; delete switch and alternative selection. |
| `B1_NARROW` | off | `crates/dft-solvers/src/davidson_stages.rs:249` | Freeze default; delete switch and alternative selection. |
| `B1_MARGIN` | 0.5 | `crates/dft-solvers/src/davidson_stages.rs:250` | Freeze default; delete switch and alternative selection. |
| `B1_ROWS` | None | `crates/dft-solvers/src/davidson_stages.rs:251` | Freeze default; delete switch and alternative selection. |
| `B6` | on | `crates/dft-solvers/src/davidson_stages.rs:268` | Freeze default; delete switch and alternative selection. |
| `SECOND_VECTOR` | off | `crates/dft-solvers/src/real_davidson.rs:1326` | Freeze default; delete switch and alternative selection. |
| `MGS` | blocked / classic (MKL / pure) | `crates/dft-solvers/src/distributed_davidson.rs:806`; `crates/dft-solvers/src/real_davidson.rs:2142` | Freeze default; delete switch and alternative selection. |
| `ORTHO_BLOCK` | off / on (MKL / pure) | `crates/dft-solvers/src/distributed_davidson.rs:803`; `crates/dft-solvers/src/real_davidson.rs:2149` | Freeze default; delete switch and alternative selection. |
| `WIDE_DOT` | on | `crates/dft-solvers/src/real_davidson.rs:2147` | Freeze default; delete switch and alternative selection. |
| `PROJECT_WIDE` | off | `crates/dft-solvers/src/real_davidson.rs:2151,3126` | Freeze default; delete switch and alternative selection. |
| `DAVIDSON_ROTATE_COLUMNS` | off | `crates/dft-solvers/src/real_davidson.rs:1894` | Freeze default; delete switch and alternative selection. |
| `DAVIDSON_UPPER_BORDER` | off | `crates/dft-solvers/src/real_davidson.rs:1692` | Freeze default; delete switch and alternative selection. |
| `LEAN_BORDER` | separate | `crates/dft-solvers/src/davidson_stages.rs:274` | Freeze default; delete switch and alternative selection. |
| `DAVIDSON_PACKED_REDUCE` | off | `crates/dft-solvers/src/distributed_davidson/workspace/lean.rs:350` | Freeze default; delete switch and alternative selection. |
| `MIX_INCREMENTAL` | on | `crates/dft-solvers/src/mixing.rs:301` | Freeze default; delete switch and alternative selection. |
| `MKL_NONLOCAL` | on | `crates/dft-physics/src/real_gamma.rs:216` | Freeze default; delete switch and alternative selection. |
| `APPLY_LEAN` | on | `crates/dft-physics/src/real_gamma.rs:614` | Freeze default; delete switch and alternative selection. |
| `HAPPLY_OVERLAP` | off | `crates/dft-physics/src/real_gamma.rs:616` | Freeze default; delete switch and alternative selection. |
| `HAPPLY_FOLD` | off | `crates/dft-physics/src/real_gamma.rs:622` | Freeze default; delete switch and alternative selection. |

### Retained boundaries

- CLI readers (including reporting/validation) retain `EIGEN_FACTOR`, `MPI_M2`,
  `M2_BATCH`, `M2_CAPACITY`, `M2_MIX`, `B1`, `B1_SKETCH`, `B5`, `RITZ_FUSED`,
  `NARROW_PROJ`, `NARROW_REDUCED`, `ORTHO`, `DAVIDSON`, `LEAN_PROJ_CACHE`,
  `RESTART_KEEP`, and `EWALD_CACHE`. Their call sites are outside this ticket;
  deleting metadata readers is not inferred to be authorized.
- `M2_DGEMM`, `M2_FFT`, `M2_STICKS`, FFT/nonlocal shared-kernel switches,
  `NL_ENERGY`, `ENERGY_IDENTITY`, `STANDARD_RESIDUAL_GUARD`, shared parallel/grid
  settings and reference-driver switches retain non-oracle consumers.
- Profiling/diagnostics remain under R11(e), including `DAVIDSON_PROFILE`,
  `DAVIDSON_STAGES`, `PHASE_WAIT`, `SNAPSHOT_DIR`, `INIT_PROFILE`, `GLUE_PROFILE`,
  and the shared timing/infrastructure controls. `REFERENCE_STRICT` remains the
  oracle routing selector.
- `M2_DENSITY` had no reader and is removed as an orphan. Backend-dependent MGS
  and ORTHO_BLOCK defaults stay backend-dependent; no MKL/pure arithmetic is
  interchanged. Numerical eligibility, rank/Gram, missing-operator and
  certification fallbacks remain.

### Deleted code and test migration

The allocating distributed Davidson solver and its orthogonalization helpers,
second-vector corrections, alternative distributed root/bordered solves, packed
border reduction, fused lean border, serial H-apply overlap/fold/rollback,
non-incremental Pulay update, geometric mixer experiment and serial/distributed
Coulomb guards are removed. Shared solver traits, factorization/panel kernels and the
legacy default lean/workspace paths remain. Historical reporting counters remain
where CLI consumers require them.

Tests of deleted alternatives are retired; default-path shape, residual,
Hermiticity, thread/rank independence, dense-reference and identity coverage is
retained. B6's plain certificate comparison uses an operator without projector
cache capability, exercising the real default fallback. Serial orthonormality
checking is retained only as an independent test oracle for the split checker.
`b1_narrow_margin.rs` is a documented empty harness because its explicit Cargo
target declaration is outside the owned files.

### R11b verification

Local checks of the 46-switch patch pass. The broader **untunable oracle** goal
is not claimed complete: the 16 CLI-read switches listed above remain, including
legacy-only arithmetic controls (`B1`, `B5`, `DAVIDSON`, `ORTHO`, etc.). Their
non-oracle readers are largely reporting, not reference arithmetic. Removing
them requires extending ownership to CLI source (`run.rs`, plus the
`EIGEN_FACTOR` pin in `run_config.rs`); `MPI_M2` also has its selector in
`dft-engine/src/parallel.rs`. These files are not in this ticket's `OWNED_FILES`.
The remaining boundary needs the orchestrator's scope decision; it is not a
claim that reporting readers make an oracle tunable by design.

| Measure (tracked Rust files throughout the repository) | Before | After | Change |
| --- | ---: | ---: | ---: |
| Manifest switches | 121 | 75 | -46 |
| Rust source lines (423 tracked `.rs` files) | 197237 | 192409 | -4828 |
| Literal `#[test]` annotations | 1182 | 1158 | -24 |

The following commands exited **0**. Logs and full recorded benchmark JSON are
under `/tmp/r11b-evidence/` (local scratch, not portable reference evidence).

```sh
cargo check --workspace --tests
cargo check -p dft-engine --no-default-features
cargo check -p dft-solvers -p dft-numerics -p dft-mkl --tests --features dft-numerics/mkl
CHECK_SCCACHE=0 BAANDI_TEST_TMP=/tmp/r11b-evidence scripts/check.sh
cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --test oracle_identity --test reference_identity
cargo build --release -p dft-cli
BAANDI_REFERENCE_STRICT=on python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /Volumes/990-2TB/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --repeats 1 --baseline docs/performance/baseline-oracle.json --require-identical --record /tmp/r11b-evidence/oracle-final.json
python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /Volumes/990-2TB/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --repeats 1 --baseline docs/performance/baseline-reference.json --require-identical --record /tmp/r11b-evidence/reference-final.json
python3 -I scripts/make_readme_figures.py
python3 -I scripts/tests/test_readme_figures.py
git diff --check
```

Full gate: **930 passed, 0 failed, 217 ignored, 103 binaries**, with build, test,
fmt, clippy, production-check and production-clippy PASS. `oracle_identity` and
`reference_identity` each pass their one unchanged golden test. README generation
has no diff; all nine freshness tests pass. Early checks exposed obsolete
alternative-arm tests/pins and compiler leftovers; these were retired or migrated
to default-path/independent-reference coverage before the successful full gate.
No tolerance or golden was relaxed. The sandbox denies `ps` and `nice` load/priority
probes; the gate still runs all tests and checks. MKL is compile-checked only.

The benchmark command checks density identity and energy tolerance. An additional
Python assertion checked exact equality of **all seven** recorded fields for
all six path/case pairs: energy bits, density SHA-256, SCF iterations,
operator applications, Davidson iterations, certificate applications and
certificate iterations. All agree. UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.

| Path | Case | Energy (Ha) | SCF iterations | Delta energy (Ha) | Density |
| --- | --- | ---: | ---: | ---: | --- |
| oracle | gamma | -7.295530468210924 | 8 | 0 | bit-identical |
| oracle | k222s | -7.91942341246909 | 10 | 0 | bit-identical |
| oracle | k222 | -7.831230057121148 | 10 | 0 | bit-identical |
| reference | gamma | -7.295530466374794 | 7 | 0 | bit-identical |
| reference | k222s | -7.919423412321121 | 7 | 0 | bit-identical |
| reference | k222 | -7.831230055901806 | 7 | 0 | bit-identical |

These are BAANDI regression baselines, not QE reference runs or performance
results. Local synthetic fixtures include simulated ranks; real MPI/MKL execution
and broader dev08 RG-ORC evidence are not claimed. Baseline JSON, goldens, README
figures and non-owned artifact/dashboard files remain unchanged. Changes are
uncommitted for review.

## R11b part 2: remaining oracle arithmetic controls (2026-11-25)

Base: `6885f8c75c5426aa341a36c51b0a6bfac9eeccda`. The user extended ownership
to CLI run/configuration and engine parallel dispatch, resolving the Part 1
blocker above. All sixteen remaining controls are retired. CLI JSON fields retain
the former default as a constant where present; retired controls are absent from
the switch manifest and input configuration. Only `REFERENCE_STRICT` selects the
oracle. Shared kernels and reference controls remain registered.

The following is the pre-change `git grep Id::NAME` reader audit (exact ID
boundaries; declarations excluded). Test readers are included. Paths are relative
to the repository root.

| Removed switch | Frozen default | Former readers |
| --- | --- | --- |
| `EIGEN_FACTOR` | 0.1 | `crates/dft-cli/src/run_config.rs`, `crates/dft-engine/src/oracle/session.rs`, `crates/dft-numerics/src/switches.rs`, `crates/dft-numerics/tests/all/switches.rs` |
| `MPI_M2` | on | `crates/dft-cli/src/run.rs`, `crates/dft-engine/src/parallel.rs`, `crates/dft-numerics/tests/all/switches.rs` |
| `M2_BATCH` | on | `crates/dft-cli/src/run.rs`, `crates/dft-engine/src/oracle/m2.rs`, `crates/dft-numerics/tests/all/switches.rs` |
| `M2_CAPACITY` | auto (3 bands at >=4 ranks; otherwise 6) | `crates/dft-cli/src/run.rs`, `crates/dft-engine/tests/m2_entry_cost.rs`, `crates/dft-numerics/tests/all/switches.rs` |
| `M2_MIX` | on | `crates/dft-cli/src/run.rs`, `crates/dft-engine/src/oracle/mix_dist.rs`, `crates/dft-numerics/tests/all/switches.rs` |
| `B1` | on | `crates/dft-cli/src/run.rs`, `crates/dft-numerics/tests/all/switches.rs`, `crates/dft-solvers/src/davidson_stages.rs` |
| `B1_SKETCH` | count | `crates/dft-cli/src/run.rs`, `crates/dft-numerics/tests/all/switches.rs`, `crates/dft-solvers/src/davidson_stages.rs` |
| `B5` | on (serial) | `crates/dft-cli/src/run.rs`, `crates/dft-numerics/src/switches.rs`, `crates/dft-numerics/tests/all/switches.rs`, `crates/dft-solvers/src/davidson_stages.rs` |
| `RITZ_FUSED` | off | `crates/dft-cli/src/run.rs`, `crates/dft-numerics/src/switches.rs`, `crates/dft-numerics/tests/all/switches.rs`, `crates/dft-solvers/src/davidson_stages.rs` |
| `NARROW_PROJ` | off | `crates/dft-cli/src/run.rs`, `crates/dft-numerics/tests/all/switches.rs`, `crates/dft-solvers/src/davidson_stages.rs` |
| `NARROW_REDUCED` | off | `crates/dft-cli/src/run.rs`, `crates/dft-numerics/tests/all/switches.rs`, `crates/dft-solvers/src/davidson_stages.rs` |
| `ORTHO` | DGKS 0.1 | `crates/dft-cli/src/run.rs`, `crates/dft-numerics/src/switches.rs`, `crates/dft-numerics/tests/all/switches.rs`, `crates/dft-solvers/src/real_davidson.rs` |
| `DAVIDSON` | lean for Ritz-change; orthonormal for strict | `crates/dft-cli/src/run.rs`, `crates/dft-engine/src/oracle/snapshot.rs`, `crates/dft-numerics/tests/all/switches.rs`, `crates/dft-solvers/src/davidson_stages.rs`, `crates/dft-solvers/tests/lean_davidson.rs` |
| `LEAN_PROJ_CACHE` | off | `crates/dft-cli/src/run.rs`, `crates/dft-numerics/tests/all/switches.rs`, `crates/dft-solvers/src/davidson_stages.rs` |
| `RESTART_KEEP` | off | `crates/dft-cli/src/run.rs`, `crates/dft-numerics/tests/all/switches.rs`, `crates/dft-solvers/src/davidson_stages.rs` |
| `EWALD_CACHE` | on | `crates/dft-cli/src/run.rs`, `crates/dft-engine/src/oracle/scf.rs`, `crates/dft-engine/tests/all/run_config.rs`, `crates/dft-numerics/tests/all/switches.rs` |

Deleted alternative implementations include dense B1 sketches, parallel B5
bordered updates, fused/deferred Ritz residual paths, narrow projection/reduced
paths, classic Ritz-change Davidson selection, lean projector caching and restart
history. Strict residual certification, serial bordered updates and their dense
fallback, CountSketch, lean Ritz-change algebra, and numerical fallback test hooks
remain. Generic numerics kernels and shared distributed tuning helpers remain.
`ResearchOptions` now contains diagnostics only. Retired branch-only test targets
retain empty harness files because their explicit Cargo manifest entries are
outside ownership. Surviving tests use independent dense spectra, residuals, Gram
checks and thread/simulated-rank identity without changing tolerances or goldens.

### Part 2 validation

Both release baseline comparisons pass all three cases (`gamma,k222s,k222`), with
zero energy difference in Ha and bit-identical density. An additional exact
comparison passes all seven recorded fields for all six case/path pairs: energy
bits, density SHA-256, SCF iterations, operator applications, Davidson iterations,
certificate applications and certificate iterations. Energies and iteration counts
are exactly those in the Part 1 table above. Evidence is in `/tmp/r11b-part2/`
(`oracle-final.json`, `reference-final.json` and command logs). These local checks
are regression evidence, not performance measurements or new QE results.

Commands use the same UPF, cases and baselines as Part 1, with evidence output
paths under `/tmp/r11b-part2/`. Full verification results and source counts follow.

Part 2 source counts (all tracked `.rs` files, physical lines and literal
`#[test]` attributes, compared with the base commit): **192,409 → 189,411
lines** (−2,998); **1,158 → 1,143 test attributes** (−15). Registry entries:
**75 → 59** (−16); cumulative R11b retirement: 121 → 59 (−62).

The full check passes **920 tests, 0 failed, 212 ignored, 103 binaries**;
all six gates pass (build, test, fmt, clippy, production check/clippy).
`oracle_identity` and `reference_identity` each pass their unchanged golden test.
The workspace test compilation, engine no-default-features compilation, and
solver/numerics/MKL feature compilation pass. MKL is type-checked only because
`MKLROOT` is unset. README generation is unchanged; all nine freshness tests pass.
The sandbox denies the script's optional `ps` load probe; verification continues.
The initial full check exposed obsolete configuration assertions and unsuitable
new dense comparisons for loosely converged buffer bands. Those branch-comparison
tests were retired; existing independent residual and orthogonality checks remain,
and no tolerance or golden was relaxed.

Reproduction (all final commands exit 0):

```sh
cargo check --workspace --tests
cargo check -p dft-engine --no-default-features
cargo check -p dft-solvers -p dft-numerics -p dft-mkl --tests --features dft-numerics/mkl
CHECK_SCCACHE=0 CHECK_NICE=0 BAANDI_TEST_TMP=/tmp/r11b-part2 scripts/check.sh
cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings
cargo test --workspace --test oracle_identity --test reference_identity
cargo build --release -p dft-cli
BAANDI_REFERENCE_STRICT=on python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /Volumes/990-2TB/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --repeats 1 --baseline docs/performance/baseline-oracle.json --require-identical --record /tmp/r11b-part2/oracle-final.json
python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /Volumes/990-2TB/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222 --repeats 1 --baseline docs/performance/baseline-reference.json --require-identical --record /tmp/r11b-part2/reference-final.json
python3 -I scripts/make_readme_figures.py
python3 -I scripts/tests/test_readme_figures.py
git diff --check
```

Shared distributed Hamiltonian batching kernels remain because `refstart` and
`refscf/kept_hx` also use that Hamiltonian. The oracle always requests the default
batching. Diagnostic-only switches remain; they do not select arithmetic.
No real MPI/MKL execution, broader dev08 RG-ORC evidence or speed improvement is
claimed. Baselines, goldens, generated README figures and non-owned dashboard
artifacts remain unchanged. All changes are uncommitted for review.

The final caller audit also removes the two physics expectation-from-cached-
projection methods and their branch-only test; neither had a reference/default
caller after removal of the lean projector cache.

## CL1: speed-sprint switch consolidation (2026-11-27)

Base: main `d63fd1ac` (table 103, after GV1-VERIFY-COST). Source/recreation
commit for every removed switch: `d63fd1ac` (`git show d63fd1ac:<path>`
restores the switch entry, its reader and the off arm; the readers are listed
in the "Removed code" paragraph below). User decision 2026-11-27 (speed
sprint done, clean up before RPA) under refactor R decision D3: remove
rollback switches whose new default was measured and adopted on dev08 and is
bit-identical or accepted, and rejected experiments; keep (a) user/policy
options, (b) research switches under evaluation, (c) switches the strict
oracle or tests need as A/B oracles, (d) anything on the verification/stop
code of GV1. Pre-sprint MPI Gamma switches (M2_*) and shared backend
selections were not reviewed against dev08 evidence here; they stay and are
listed as follow-up candidates.

Result: **26 removed, 77 kept**; manifest 103 -> 77 (both switch-count
tests updated: `dft-numerics` `TABLE.len()` and the CLI `scf.switches` count).

| Switch | Decision | Reason | Evidence |
| --- | --- | --- | --- |
| `NL_ENERGY` | keep: (c) oracle option | exact band-by-band nonlocal energy is the independent oracle of the blocked path; shared energy kernel | R11a/R11b inventory |
| `ENERGY_IDENTITY` | keep: (b) research | default-off band-sum identity, no rejection row | ledger LANE-GLUE |
| `STANDARD_RESIDUAL_GUARD` | keep: (d) stop rule | stop-rule guard of the spin/FD driver (GV1 area) | ledger 2026-11-27 gen59 arm D |
| `M2_DGEMM` | keep: backend | dgemm vs pure panels (MKL/pure builds); generality/backend selection | R11a |
| `M2_FFT` | keep: backend | rustfft vs MKL lanes; backend selection | R11a; ledger work-order screen 2026-11-26 |
| `STAGE_TIMERS` | keep: diagnostic | measurement only | R11(e) |
| `M2_STICKS` | keep: outside CL1 | pre-sprint MPI Gamma kernel shared with the strict oracle; not reviewed in CL1 (follow-up candidate) | R11b boundary |
| `MPI_TREE` | keep: infrastructure | collective pattern selection | R11a |
| `M2_SCREEN` | keep: outside CL1 | pre-sprint MPI Gamma path; follow-up candidate | R11b boundary |
| `M2_ATOMIC` | keep: outside CL1 | pre-sprint MPI Gamma start; follow-up candidate | R11b boundary |
| `M2_ATOMIC_BLOCK` | keep: outside CL1 | pre-sprint MPI Gamma start; follow-up candidate | R11b boundary |
| `M2_GRID` | keep: outside CL1 | pre-sprint MPI grid stages, CLI routing reader; follow-up candidate | R11b boundary |
| `M2_LOCAL_FIELDS` | keep: outside CL1 | pre-sprint M2-OWNED; follow-up candidate | R11b boundary |
| `M2_LOCAL_DENSITY` | keep: outside CL1 | pre-sprint M2-OWNED; follow-up candidate | R11b boundary |
| `DAVIDSON_STAGES` | keep: diagnostic | stage timers | R11(e) |
| `DAVIDSON_PROFILE` | keep: diagnostic | lane trace | R11(e) |
| `NONLOCAL_BLOCK` | keep: backend | auto/on/off block nonlocal of the complex Hamiltonian, shared with the oracle | R11b boundary |
| `ROTATE_KERNEL` | keep: backend | shared numerics kernel | R11b boundary |
| `PURE_GEMM` | keep: backend | pure vs MKL build default | R11b boundary |
| `PACK_FIXED` | keep: backend | aarch64 pure kernel option | R11b boundary |
| `MKL_DGEMM` | keep: backend | MKL builds | R11a |
| `MKL_EVD` | keep: backend | MKL builds | R11a |
| `PAIR_FFT` | keep: backend | MKL builds | R11a |
| `BINARY_HASH_CACHE` | keep: infrastructure | launch I/O cache | R11a |
| `RELAX_EXTRAPOLATION` | keep: (a) user option | overrides an input key | R11a |
| `REFERENCE_STRICT` | keep: (c) oracle routing | selects the strict oracle | ADR 0043, R11a |
| `DGEMM_TIMERS` | keep: diagnostic | exit-time timers | R11(e) |
| `PSEUDO_DIR` | keep: infrastructure | cache location | R11a |
| `FFT_KEPT_FUSE` | keep: (b) closed neutral | bit-identical, closed neutral, default off; no rejection for correctness; follow-up candidate | ledger FFT-KEPT-FUSE-NEUTRAL |
| `COLLECTIVE_TRACE` | keep: diagnostic | - | R11(e) |
| `FFT_PROFILE` | keep: diagnostic | - | R11(e) |
| `GLUE_PROFILE` | keep: diagnostic | - | R11(e) |
| `PHASE_PROFILE` | keep: diagnostic | - | R11(e) |
| `EWALD_FAST` | keep: (c) oracle | off = reference Ewald scheme, used as the tight reference in tests | ledger INIT-EWALD-FAST |
| `IONIC_FAST` | keep: (c) oracle | off = direct phase reference (setup_memo_identity) | ledger INIT-EWALD-FAST |
| `INIT_PROFILE` | keep: diagnostic | - | R11(e) |
| `PHASE_VALIDATE` | keep: (b) research | default off pending dev08 | ledger G-PHASE-VALIDATE |
| `M2_LAZY_RESULT` | keep: (c) test oracle | off (eager delivery) is the A/B oracle of the CLI lazy-delivery tests (Gamma and pooled k) | ledger GK18/GK18B |
| `PHASE_WAIT` | keep: diagnostic | - | R11(e) |
| `SNAPSHOT_DIR` | keep: diagnostic | replay capture | R11(e) |
| `FFT_SPARSE_CLEAR` | keep: (b) research | P3c inconclusive, not rejected | ledger P3c |
| `FFT_DIRECT_PACK` | keep: (b) research | P3c inconclusive, not rejected | ledger P3c |
| `FFT_FUSED_ROWS` | **remove** (former default off) | P3-FFT2 fused z rows: rejected (slower on si64/elongated/rskew) | ledger P3f/P3g screens (2026-11-25) |
| `FFT_A2A_BALANCE` | keep: (b) research | default off pending dev08 | docs/performance/lanes/RSCF/a2a-dev08-request.md |
| `NONLOCAL_SINGLE_GEMM` | **remove** (former default on) | P3-NL single projection GEMM; rounding class | ledger P3l (2026-11-25, dev08 MPI8) |
| `NONLOCAL_FUSE_COUPLING` | **remove** (former default on) | P3-ADDV fused coupling storage; bit-identical | ledger P3m (2026-11-25, dev08 MPI8) |
| `NONLOCAL_BACK_ACCUM` | **remove** (former default on) | P3-ADDV back-projection accumulated into HX; rounding class | ledger P3m (2026-11-25, dev08 MPI8) |
| `REFSCF_RMM` | keep: algorithm module | Generality rule (P3r) | ledger P3r |
| `REFSCF_RMM_STEPS` | keep: algorithm module | RMM parameter | ledger P3r |
| `REFSCF_RMM_BLOCK` | keep: algorithm module | RMM parameter | ledger P3r |
| `FFT_KERNEL_BY_LENGTH` | keep: (b) research | default off, no rejection row | ledger work-order screen 2026-11-26 |
| `REFSCF_PERMAP` | keep: (b) research | no rejection row (R11a) | R11a |
| `REFSCF_LADDER_RATIO` | keep: algorithm module | reference ladder (D3 keep) | refactor-r-plan D3 |
| `REFSCF_LADDER` | keep: algorithm module | reference ladder (D3 keep) | refactor-r-plan D3 |
| `REFSCF_KEPT_RECOMPUTE` | keep: (a) memory option | memory-saving option | ledger P4 screen |
| `REFSCF_MEMORY_TRIM` | keep: (b) research | pending dev08, no rejection row | R11a |
| `K_LEAN_ENERGY` | **remove** (former default on) | GK4 deferred per-k energy terms; bit-identical | ledger GK4; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `FORCE_BAND_FIRST` | **remove** (former default off) | G5 band-first force GEMM: never measured, unreachable since GD1's GEMM kernel became the default (superseded experiment) | docs/design/reference-scf.md (G5); ledger GD1-DERIV-SPEED |
| `REFSCF_SPIN_FD_MOMENT` | keep: routing | CLI driver routing reader | R11a |
| `XC_SCREEN_NEGATIVE` | keep: physics policy / (c) | ADR 0046 policy; off is the historical fixture contract of derivative tests | ADR 0046 |
| `REFSCF_REAL_START` | **remove** (former default off) | G4G real Gamma start: rejected (centred fixture stalls 200 maps) | docs/design/reference-scf.md G4G decision (2026-11-26) |
| `REFSCF_CANONICAL_CLUSTER` | keep: (b) research | G4G research switch kept by decision | docs/design/reference-scf.md G4G |
| `CDAV_REUSE` | **remove** (former default on) | GK3 kept complex Davidson workspace and merged agreements; bit-identical | ledger GSPEED-PLAN/GK3 (solver test GK3); default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `REFSCF_SPIN_PULAY` | keep: (b) research | G4G candidate pending a dev08 O run | docs/design/reference-scf.md G4G |
| `K_START_MEMO` | **remove** (former default on) | GK1 memoized k atomic columns; bit-identical | ledger GK1-START-MEMO; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `K_NL_TABLE_MEMO` | **remove** (former default on) | GK1B rows-only k projector tables; bit-identical | ledger GK1B-NL-TABLES; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `K_START_SPIN_SHARE` | **remove** (former default on) | GK1B shared collinear spin start; bit-identical | ledger GK1B-NL-TABLES; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `CDAV_BORDER_GRAM` | **remove** (former default on) | GK10 border-only projected S/H; rounding class, strict-class agreement | ledger GK10-BORDER-GRAM; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `RADIAL_TABLE_START` | **remove** (former default on) | GK14 radial tables for the k start; strict class | ledger GK14-RADIAL-TABLES; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `RADIAL_TABLE_NL` | **remove** (former default on) | GK14 radial tables for k projectors; strict class | ledger GK14-RADIAL-TABLES; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `CDAV_BORDER_FACTOR` | **remove** (former default on) | GK10B bordered Cholesky; rounding class | ledger GK10B-CDAV-REALPATH; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `CDAV_STREAM` | **remove** (former default on) | GK10B streamed correction pass; bit-identical | ledger GK10B-CDAV-REALPATH; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `K_FIRST_EMPTY_RULE` | keep: (b) research | named by the user decision | ledger GK10B |
| `K_START_KEPT_HX` | **remove** (former default on) | GK10B kept start H X for k map 1; rounding class | ledger GK10B-CDAV-REALPATH; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `K_REDO_KEPT_HX` | **remove** (former default on) | GK10B kept HX on the first-map redo; rounding class | ledger GK10B-CDAV-REALPATH; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `REFSCF_K_POOLS` | keep: (a) user option | documented k-pool policy | docs/cli.md, ledger GK9 |
| `REFSCF_K_NPOOL` | keep: (a) user option | documented pool count | docs/cli.md, ledger GK9B |
| `RADIAL_TABLE_FIELDS` | **remove** (former default on) | GK16 radial tables for k-path atomic fields; strict class | ledger GK16-SETUP-TAIL; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `FFT_ONE_RANK_DIRECT` | **remove** (former default on) | GK17 one-rank direct stick FFT for k operators; bit-identical | ledger GK17-SCF-GAP; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `K_KEPT_SLABS` | keep: (b) research | named by the user decision | ledger GK17 |
| `K_KEPT_SLAB_MB` | keep: (b) research | budget of K_KEPT_SLABS | ledger GK17 |
| `SPIN_FD_START_KEPT_HX` | **remove** (former default on) | GF1 port of the kept start HX to the spin/FD driver; rounding class | ledger GF1-FD-PATH; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `SPIN_FD_REDO_KEPT_HX` | **remove** (former default on) | GF1 port of the kept redo HX to the spin/FD driver; rounding class | ledger GF1-FD-PATH; default since its merge; ran in the dev08 campaigns up to gen66 (final MP campaign at fb419f49, README gen66 summary) with QE-matched results |
| `STRESS_FORCE_PASS` | **remove** (former default off) | GD1: stress no longer re-runs the force pass; stress bit-identical | ledger GD1-DERIV-SPEED; dev08 gen60/gen66 derivative stages |
| `DERIV_NL_GEMM` | **remove** (former default on) | GD1 GEMM nonlocal derivative kernel; rounding class | ledger GD1-DERIV-SPEED; dev08 gen60/gen66 |
| `DERIV_RADIAL_TABLES` | **remove** (former default on) | GD1 tabulated projector radials in the derivative pass; rounding class | ledger GD1-DERIV-SPEED; dev08 gen60/gen66 |
| `DERIV_SHELLS` | **remove** (former default on) | GD1 shared local/core shell values; bit-identical | ledger GD1-DERIV-SPEED; dev08 gen60/gen66 |
| `DERIV_SHELL_TABLES` | keep: not yet measured | GD2 merged after gen66; dev08 pending | ledger GD2-DERIV-XC |
| `DERIV_XC_SPLIT` | keep: not yet measured | GD2 merged after gen66; dev08 pending | ledger GD2-DERIV-XC |
| `DERIV_JOINT` | keep: not yet measured | GD2 merged after gen66; dev08 pending | ledger GD2-DERIV-XC |
| `FD_EMPTY_THRESHOLD` | keep: (d) verification-coupled | FD_EMPTY_VERIFY depends on it; three-way policy | ledger GF2 |
| `FD_EMPTY_VERIFY` | keep: (d) verification | verification map (GV1 area) | ledger GF2, gen63 |
| `FD_BUFFER` | keep: (b) documented experiment | named by the user decision | ledger GF2 |
| `COLLINEAR_TIME_REVERSAL` | keep: (c) oracle reader | read by the strict oracle (oracle/collinear.rs) and CLI config | ADR 0050 |
| `FD_GROWTH_AT_STOP` | keep: (d) stop/verification | growth judged at stop and verification maps (GV1 area) | ledger GB2 |
| `FD_VERIFY_PREDICTOR` | keep: (d) verification | verification predictor (GV1 area) | ledger GB2 |
| `SPIN_FD_ATOMIC_DENSITY` | keep: (b) research | named by the user decision | ledger GB2 |
| `K_REDUCTION` | keep: (a) user option | ADR 0051; unreduced is the oracle of pointwise | ADR 0051 |
| `GRID_INCOMPATIBLE_OPS` | keep: (a) user option | ADR 0049 | ADR 0049 |
| `SCF_ETHR_GUARD` | keep: (b) research | named by the user decision | ledger GM1 |
| `VERIFY_FEED` | keep: (d) GV1 research | dev08 evaluation | ledger GV1-VERIFY-COST |
| `VERIFY_SKIPS` | keep: (d) GV1 research | dev08 evaluation | ledger GV1-VERIFY-COST |
| `VERIFY_RESPONSE_PREDICTOR` | keep: (d) GV1 research | dev08 evaluation | ledger GV1-VERIFY-COST |
| `RESPONSE_VELOCITY` | keep: diagnostic | R2-VEL-CHECK: q -> 0 dipoles with full / local-only / l <= 1 nonlocal velocity; default full; results of the default unchanged | tests/reference/response/r2/README.md |
| `RESPONSE_DIPOLE_DUMP` | keep: diagnostic | R2-VEL-CHECK: dipoles.csv of baandi dielectric for external dipole comparison; default off; results unchanged | tests/reference/response/r2/README.md |

### Removed code

- `dft-solvers/src/ref_davidson_complex.rs`: per-call workspace allocation,
  separate agreements and operand rescans (CDAV_REUSE off), full S/H
  recompute after an append (CDAV_BORDER_GRAM off), full refactor
  (CDAV_BORDER_FACTOR off), per-element correction loop (CDAV_STREAM off).
- `dft-engine/src/refstart.rs`: direct per-(row, atom, m) quadrature columns
  (K_START_MEMO), full-sphere projector tables (K_NL_TABLE_MEMO), per-channel
  spin start (K_START_SPIN_SHARE), quadrature arms of the k start, projector
  and field tables (RADIAL_TABLE_START/NL/FIELDS, also the fixed-density
  start), the switch selection of the kept start HX (K_START_KEPT_HX,
  SPIN_FD_START_KEPT_HX: the SCF starts always keep it), the real Gamma start
  with mirror complements (REFSCF_REAL_START).
- `dft-engine/src/refscf/k_map.rs`: per-k residual reduction and eager
  per-map energy terms (K_LEAN_ENERGY off).
- `dft-engine/src/refscf/k_driver.rs`, `refscf/electronic.rs`: the kept-HX
  invalidation before the first-map redo (K_REDO_KEPT_HX,
  SPIN_FD_REDO_KEPT_HX off); no verification or stop logic touched.
- `dft-solvers/src/distributed_fft.rs`: the transpose path selection for
  one-rank k plans (FFT_ONE_RANK_DIRECT; plain `new` plans keep the transpose
  path) and the whole fused-row workspace (FFT_FUSED_ROWS);
  `dft-physics/src/distributed_complex.rs::kept_planes_exact` is now always
  true.
- `dft-physics/src/distributed_real_gamma.rs`: panel projection route, the
  separate coupling buffer W and the separate back-projection add
  (NONLOCAL_SINGLE_GEMM/FUSE_COUPLING/BACK_ACCUM off; MKL-only arithmetic,
  pure builds unchanged); `dft-cli/src/run.rs` reports the single-GEMM route
  unconditionally.
- `dft-engine/src/forces.rs`, `stress.rs`: per-k operator build with scalar
  loops (DERIV_NL_GEMM off), quadrature projector radials
  (DERIV_RADIAL_TABLES off), per-shell evaluation (DERIV_SHELLS off; the
  kernels' `Option` interface stays), the complete force pass inside the
  stress (STRESS_FORCE_PASS on; `precondition_force_band_passes` stays in the
  record as 0), `forces::nonlocal_band_first` (FORCE_BAND_FIRST).

### Test migration

Tests that compared a retired arm with the default keep their evidence in the
ledger rows named in the table and now check the default path only: rank-count
independence within the unchanged classes, bit reproducibility, worker-count
independence, convergence, kept-HX hits. Deleted: `dft-engine/tests/k_radial_table.rs`
(only table-vs-quadrature arms), the two GK16 strict-class A/B tests and the
real IrI3 A/B in `k_setup_tail.rs` (the table accuracy scan stays),
`band_first_identity_complex_k_fractional_occupations`, the solver unit tests
`border_gram_matches_full_recompute` and `border_factor_and_stream_agree`
(replaced by `border_path_is_rank_independent`), the plan unit test
`fused_rows_respects_backend_selection` and
`fused_rows_preserve_transforms_and_remove_slab_passes`. The direct-stick
plan test now compares `new_direct` with `new` (same bits). The G4G stall
record test (`dft-cli/tests/reference_stalls.rs`) pinned three retired arms to
reach a rounding-dependent stall; with the adopted paths the same fixture
reaches it at mixing beta 0.9 / history 20 (input change only, no
tolerance changed). No tolerance or golden was relaxed.

### CL1 validation (Mac, pure Rust, 2026-11-27)

| Measure | Before (`d63fd1ac`) | After | Change |
| --- | ---: | ---: | ---: |
| Manifest switches | 103 | 77 | -26 |
| Rust lines under `crates/` (diff) | | | +424 / -2223 |
| Literal `#[test]` annotations | 1289 | 1277 | -12 |

All commands below exited **0** (logs in the worktree's `target/cl1-evidence/`,
local scratch, not portable evidence):

```sh
cargo clippy --workspace --all-targets -- -D warnings
cargo check -p dft-solvers -p dft-numerics -p dft-mkl -p dft-physics --tests --features dft-numerics/mkl
cargo check -p dft-engine --no-default-features
cargo check -p dft-cli --features mpi
CHECK_SCCACHE=0 CHECK_NICE=0 scripts/check.sh
cargo build --release -p dft-cli
BAANDI_REFERENCE_STRICT=on python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf .../Si.pz-vbc.UPF --cases gamma,k222s,k222 --repeats 1 --baseline docs/performance/baseline-oracle.json --require-identical
cargo test --release -p dft-engine --test k_setup_fingerprint -- --ignored --nocapture   # base and candidate, diff
cargo test --release -p dft-cli --test all relax_extrapolation -- --ignored
DFT_SI_UPF=.../Si.pz-vbc.UPF cargo test --release -p dft-cli --test reference_gaps reference_k_relax -- --ignored
cargo test --release -p dft-engine --test refscf_k_border_gram border_gram_scf_converges_on_k_generality_set -- --ignored
cargo test --release -p dft-engine --test refscf_fd_kept ports_within_class_at_eight_ranks -- --ignored
python3 -I scripts/tests/test_readme_figures.py
git diff --check
```

Full gate: **1028 passed, 0 failed, 239 ignored, 122 binaries**; build, test,
fmt, clippy, production-check and production-clippy PASS. (Two earlier runs
with `BAANDI_TEST_TMP` set to a long path failed only in environment-bound
tests: the pseudo cache tests and a Unix socket probe, "path must be shorter
than SUN_LEN"; both pass with the default temp directory.) MKL and MPI are
type-checked only.

Identity of the default behaviour:

- `k_setup_fingerprint` (release, 60 runs: unpolarized k, collinear spin and
  Fermi-Dirac; 1 to 3 ranks; pools on and off; energy bits, map count and a
  hash of every delivered bit): the base build at `d63fd1ac` and the candidate
  print identical lines (diff empty).
- Reference driver `bench_scf` (gamma, k222s, k222): the stored
  `baseline-reference.json` is already stale on main for the k cases (the base
  build `d63fd1ac` fails it in the same way: density differs, dE 0 and
  8.9e-16 Ha, from the sprint's rounding-class adoptions). Base and candidate
  agree exactly in all seven recorded fields (energy bits, density SHA-256,
  SCF iterations, operator applications, Davidson iterations, certificate
  applications and iterations) for all three cases.
- Strict oracle `bench_scf`: identical to `baseline-oracle.json` (dE 0,
  density identical) on all three cases.
- `reference_identity`, `oracle_identity` and the other identity tests pass in
  the full gate with unchanged goldens.

No speed claim; dev08 not run (the NONLOCAL_* arms are MKL-only arithmetic and
only type-checked here; their retirement keeps the dev08-adopted default).
