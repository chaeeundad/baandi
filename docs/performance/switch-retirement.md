# Reference-path switch retirement (R11a)

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
