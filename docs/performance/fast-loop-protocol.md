# Fast optimization loop: screen on one case, confirm on the set (user request, 2026-11-16)

Running the whole generality set for every idea is expensive (hours of dev08 time, most of it spent on ideas that do not pay). The loop is split into four stages. The Generality rule and the Adoption rule (AGENTS.md) are unchanged: a default is adopted only after stage 3.

User instruction 2026-11-16: split the verification at the DIAGONALIZATION level. A kernel idea is first benchmarked on one eigensolve at a fixed Hamiltonian (seconds, stage 0); the SCF convergence process (SCF steps, verification attempts, mixing) is tested only after that (stage 1). Whole-SCF screens are slower (5 to 60 s per run) and noisy, because the number of SCF steps and verification attempts varies with the seed.

| stage | level | cost | decides |
|---|---|---|---|
| 0 | one Davidson solve at a fixed Hamiltonian (`scripts/fast_screen_dav.py`) | seconds per process, a si32 screen under a minute | kernel ideas: counters first, then solve time |
| 1 | one SCF, one case (`scripts/fast_screen.py`, `scripts/fast_screen_mpi.py`) | minutes | trajectory ideas and the stage-0 survivors |
| 2 | five or six cases | tens of minutes | widening |
| 3 | the full generality set and accuracy classes | hours (dev08) | adoption |

## Stage 0: diagonalization level (seconds, one solve, one binary)
Harness: `crates/dft-cli/examples/davidson_bench.rs` (the benchmark binary) and `scripts/fast_screen_dav.py` (the A/B screen, python3 -I, standard library only). The solve-level mode of the MPI screen is `scripts/fast_screen_mpi.py --level solve` (same harness on simulated ranks).

What one process does (in-process, deterministic):
- Builds the case once from a CLI input (`fast_screen_dav.py` writes it from `benchmarks/cases/<case>.toml` with the stage-1 displacement stream: seed 601, amplitude 1e-3, T=1 by default) exactly as `dft-cli run` builds it: `dft_cli::run::options` (single standard path, Ritz-change solves, lean Davidson by default, REUSE-PSI, the MPI capacity rule), real Gamma path. The cutoff ladder is not used: every solve runs at the full cutoff.
- Base density rho0: the `density.bin` of a converged `dft-cli run` of the same input (`--cli <dft-cli binary>` runs it once and caches it in `--work`; or `--density FILE`), else the problem's start density (the superposition of atomic densities, or the uniform density when the UPF carries no atomic charge, as Si.pz-vbc). The JSON field `base_density` records which. Use the converged density for threads screens (late-SCF warm solve); the start density gives a step-2-like warm solve.
- Setup (untimed): the atomic-start vectors at rho0 (computed once by the SCF's routine, `ScfProblem::atomic_start_vectors`, measurement hook, not used by the SCF), one map at rho0 from them at `--prev-tol` (1e-8 Ha) giving X_prev and rho_out; warm density rho_w = rho0 + 0.5 (rho_out - rho0); warm tolerance = the standard ratchet after that map, ratchet(1e-8, E_H[rho_out - rho0], N_e).
- **cold** solve: the map at rho0 from the atomic start at the standard first-step tolerance 5e-3 Ha (`StandardProfile::CAP_HA`). The vectors are passed to the map, so the solve is identical to a map that recomputes them (`--cold-start internal`, checked: same eigenvalue FNV).
- **warm** solve: an untimed re-arm map at rho0 from X_prev (2 Davidson iterations; it leaves X, H X and psi(r) kept for REUSE-PSI as an SCF step does), then the timed map at rho_w from the re-arm orbitals at the ratchet tolerance: the typical mid-SCF solve with the kept start taken.
- One discarded warm-up and N=5 timed repeats per kind. The solve time is the diagonalization phase of the map (the solver call only); the map time (potential, solve, density, audit) is reported beside it. Counters are deltas over each timed map: H columns and Davidson iterations of the solve; GEMM calls and FLOP per site/route (`gemm_dispatch`); the serial Davidson GEMM FLOP and full reduced eigensolves (sum m^3); lean Davidson residual and appended columns, contractions, kept starts, projection-cache columns; REUSE-PSI and audit counters; on ranks: collective calls and bytes per operation, distributed stage calls (rotate, residual, reduced_gram, orthonormalize, reduced_eigen, apply), solver reductions (calls, doubles), reduced solves and sum m^3, rotation FLOP. Map-level counters include the density build and the energy audit of the map (for example the audit's projector GEMM, which `BAANDI_LEAN_PROJ_CACHE` removes); the solver's own counters are solve-only.
- FFT counts (`fft.pair_transforms`), `h_apply.*` and stage seconds are filled only with `BAANDI_DAVIDSON_STAGES=on` (timer cost; set it on BOTH arms and never use those runs for timing; without it these counters are reported as unusable). Allocation bytes are not available (no allocation counter in the engine, and the harness may not install a counting allocator, which needs `unsafe`); the screen records the peak RSS of each process instead (reported by the OS, not used by the verdict).
- Identity: eigenvalue FNV-1a per repeat (must be equal across repeats and processes of one arm: determinism), eigenvalues and residual norms, converged roots.

Screen (`fast_screen_dav.py`): three alternating pairs (A B, B A, A B) of the same binary with an environment A/B (switches are read once per process, so each arm is its own process; process start and setup are outside the timed solves). Per kind: per-process median solve time, paired ratios B/A, spread = max(1 %, half range of the ratios). Verdict REJECT (a failed run; eigenvalues not reproducible within an arm; an occupied eigenvalue moved by more than max(10 tol, 1e-10) Ha; fewer converged roots; an `--expect` counter moved the wrong way; median ratio above 1 + spread), CANDIDATE (median ratio below 1 - spread, or H columns / Davidson iterations / GEMM FLOP or an `--expect` counter down by at least 5 % with the time inside the spread), UNUSABLE (every `--expect` counter is zero on both arms: **counter zero on both arms = unusable**, never "no difference"), NEUTRAL otherwise. Fewer than three pairs: the time is not judged. Unit tests: `python3 -I scripts/tests/test_fast_screen_dav.py`.

Cost on the Mac (O, pure build, shared host, 2026-11-16): si8 one process 0.4 to 0.6 s, a three-pair screen 2.5 to 4 s; si32 one process 7 to 10 s (setup about 2 s, cold solve 0.26 s, warm solve 0.36 s at load 3 to 5), a three-pair screen 44 to 48 s at load 3 to 5 and 59 s at load 10 to 12; si32 simulated np=4 three pairs 27 to 29 s; si8 simulated np=4 about 4 s. A/A on si32: cold ratio 0.992 (spread 1 %), warm 0.989 (spread 4.6 %), bit-identical. First results: ledger section "Stage-0 Davidson screens".

Testable at stage 0 (kernel ideas, the solve or map changes at a fixed Hamiltonian): Davidson variant (`BAANDI_DAVIDSON`), Gram-Schmidt and GEMM kernels (`ORTHO_BLOCK`, `PURE_GEMM`, `MGS`, `PROJECT_WIDE` and `MKL_*` on MKL builds), projection cache (`LEAN_PROJ_CACHE`), REUSE-PSI, H apply (`APPLY_LEAN`, `HAPPLY_*`), MPI solver kernels on ranks (`M2_PANEL_GS`, `M2_SELECTED_ROOTS`, `M2_WORKSPACE`, `M2_EIGEN`, `M2_B5`, the capacity). **Not testable here** (they depend on the SCF trajectory, the stop rule or the verification): `VERIFY_REUSE`, mixing (beta, Kerker, metric, history), the cutoff ladder, the stop rule and its estimators, the verification tolerance, start-density ideas, SCF-step counts in general: these go to stage 1 directly. Latency-only MPI ideas (`MPI_TREE`) give counters only on simulated ranks (a simulated tree sum is two thread barriers); time them with real MPI on dev08 (below). Gamma-only: k-point code paths (`NONLOCAL_BLOCK`, the complex nonlocal) are not reached by this benchmark.

A stage-0 result is never adoption evidence. A CANDIDATE goes to stage 1 (does the SCF-level wall move?) and then stage 2/3; a REJECT row is written to the ledger (section "Stage-0 Davidson screens") and not repeated without new evidence.

Commands (Mac, pure build):
```
python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case si32 --cli target/release/dft-cli \
    --work <scratch dir> --b BAANDI_DAVIDSON=classic                       # threads, T=1
python3 -I scripts/fast_screen_dav.py --upf <UPF> --case si8 --b BAANDI_LEAN_PROJ_CACHE=on \
    --expect lean.proj_cache_columns:up --expect nl_audit.direct:down
python3 -I scripts/fast_screen_mpi.py --level solve --upf <UPF> --case si32 --ranks 4 \
    --b BAANDI_M2_PANEL_GS=fused --expect collectives.tree_sum.calls:down  # simulated ranks
```
Options of `fast_screen_dav.py`: `--case si8|si16|si32|si64|<path>`, `--seed`, `--amp`, `--threads`, `--ranks` (simulated), `--pairs`, `--repeats`, `--kinds cold,warm`, `--a/--b NAME=value[,..]`, `--expect counter:down|up` (summed over kinds), `--min-drop` (0.05), `--density`, `--cli`, `--features mkl`, `--binary`, `--bench-args` (passes `--prev-tol`, `--mix`, `--cold-tol`, `--warm-tol`, `--cold-start` to the binary), `--emit-input PATH` (write the CLI input and exit), `--json-out`, `--all-counters`.

### dev08 usage (master session; exact commands)
Build once per commit (pure, MKL, MPI; separate target directories as usual), record the sha256 of each binary:
```
cargo build --release -p dft-cli --example davidson_bench                                    # pure
MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo build --release -p dft-cli --features mkl --example davidson_bench
cargo build --release -p dft-cli --features mpi,mkl --example davidson_bench                   # MPI (+MKL)
cargo build --release -p dft-cli [--features mkl]                                             # dft-cli for the base density
sha256sum target/release/examples/davidson_bench
```
Threads (T=1 and T=8; `taskset -c 0-(T-1)`, `MKL_NUM_THREADS=1`, load below 0.30 as for every dev08 row):
```
python3 -I scripts/fast_screen_dav.py --upf <Si.pz-vbc.UPF> --case si32 --threads 1 \
    --binary target/release/examples/davidson_bench --cli target/release/dft-cli --work ~/fastbench \
    --b BAANDI_X=on --json-out ~/fastbench/X-si32-T1.json
taskset -c 0-7 python3 -I scripts/fast_screen_dav.py ... --threads 8 ...
```
Real MPI (np 2/4/8; one solve per process, the binary prints rank 0's JSON with the maximum over ranks of every repeat):
```
python3 -I scripts/fast_screen_dav.py --upf <UPF> --case si32 --emit-input ~/fastbench/si32.toml
target/release/dft-cli run ~/fastbench/si32.toml --output-dir ~/fastbench/si32-cli     # base density, once
for P in 2 4 8; do for ARM in A B; do
  env $( [ $ARM = B ] && echo BAANDI_X=on ) mpirun -np $P --bind-to core \
    target/release/examples/davidson_bench --input ~/fastbench/si32.toml \
    --density ~/fastbench/si32-cli/density.bin | grep DAVBENCH_JSON > ~/fastbench/X-np$P-$ARM.json
done; done
```
(alternate the arm order per pair as the screen does; three pairs per np). Each process takes seconds (si32) and prints `solve_ms` per repeat, the per-rank counters (`coll.*`, `dist.*`) and the eigenvalue FNV; `BAANDI_MPI_TREE=rd` is read by the binary at MPI start like the CLI. Report: binary sha256, the JSON lines, load before and after.

## Stage 1: screen (minutes, one case, one binary)
- One representative case per parallel model: threads = si32 standard f1e-3 seed 601 T=1 (and T=8 when the idea touches threading); MPI = si32 standard np=4 seed 601; add si64 only for ideas about large-N terms (EVD, memory).
- Same binary, environment A/B (switch off/on), three alternating pairs, load gate as before. Report the median per-arm wall, the spread and the deterministic counters.
- Prune first on the deterministic counters, which have no noise: Davidson/H columns, GEMM FLOPs, FFT counts, collective counts and bytes, allocations. An idea that does not move a counter in the intended direction is dropped without any wall measurement.
- Pass to stage 2 when the median wall improves beyond the spread, or the counter drop is large and the wall effect is below the spread of one case.

### Counter availability (checked 2026-11-16, source 006baaf; re-check after changes)
- Always filled: `convergence.work.h_vectors_recorded` + `certificate_h_vectors`; `parallel.collectives[].{calls,payload_bytes}`; `convergence.work.gemm_dispatch.rows[].{calls,gflop,route,site,max_m,max_n,max_k}` plus `mkl_flop_fraction`, `wide_requests` (count-only, no switch; present in builds from 69e13f1 on, absent before); `parallel.reductions` (MPI, count-only).
- Only with a switch: `davidson_lean.gemm_gflop.*` needs `BAANDI_DGEMM_TIMERS` (timer cost; use it on both arms and never for acceptance timing); `davidson_stages.fft.*` is all zero on the MPI M2 / lean paths seen so far (si64 strict MPI at 69e13f1): treat the FFT counters as unusable there.
- Not available: heap allocation counts (only the workspace allocation tests of the MPI solver). Use measured peak RSS from outside the process; add an engine counter before an allocation-based screen is relied on.
- Harness self-check (mandatory): a counter that is zero on both arms is reported as "unusable (all zero)" and excluded from the screen. A zero never means "no difference".

## Stage 2: widen (tens of minutes, five to six cases)
si16, si32, si64 (standard), k444s, mgo-40, one skewed or slab case; paired env A/B; both backends when the code path differs (MKL and pure). Check energy, density SHA and counts against the arm without the change.

## Stage 3: confirm and adopt (the full set, once per candidate)
Only candidates that survived stage 2 run the full generality set and the accuracy classes (strict 1e-10 Ha, standard 5e-9 Ha, certificate, thread independence). Stacked candidates are confirmed together in one run where the code paths are independent.

## Rules
- Implementers measure stage 1 locally where possible (Mac observation, label O) before asking dev08 for stage 2; dev08 requests list the stage explicitly.
- A stage-1 negative result is recorded in the ledger as a short row (idea, counter, wall, spread) so that it is not repeated.
- Screening results are never adoption evidence. Never skip stage 3 for a default change.

## Harness (MPI)
`scripts/fast_screen_mpi.py` (python3 -I, no network; added 2026-11-16) is the stage-1 MPI screen on this Mac. It builds the release test binary of `crates/dft-engine/tests/m2_entry_cost.rs` once (sha256 printed), then runs the SAME binary on N alternating pairs (A B, B A, ...) with an environment A/B. Default case: si32, np=4 simulated ranks (threads of one process), production mode (`ScfOptions::default()`: standard 5e-9 Ha, fixed-point verification, Pulay 0.5) with the CLI capacity rule (3 per band at four or more ranks) and the C1 ladder, i.e. the CLI default path. Options: `--case si8|si16|si32|si64`, `--ranks`, `--mode production|standard|strict`, `--ladder on|off`, `--pairs`, `--a/--b NAME=value[,..]`, `--expect counter:down|up` (repeatable), `--min-drop` (default 0.10), `--json-out`, `--binary` (skip the build), `--all-counters`.
- The harness prints one `M2SCREEN_JSON` line per run (per-rank values): wall, energy bits, density FNV-1a, SCF / Davidson / H columns, rotation FLOPs, reduced solves and sum of m^3, `parallel.reductions` per site (calls, doubles), `parallel.collectives` per operation (calls, payload bytes, seconds), Davidson and M2 stage seconds, M2 counters.
- Self-check (mandatory): a counter that is zero on both arms is listed as "unusable (all zero)"; an `--expect` counter that is unusable makes the verdict UNUSABLE, never "no difference". Counters must be identical between repeats of one arm (a note is printed otherwise); the density hash must be identical between repeats of one arm (determinism contract) or the verdict is REJECT.
- Verdict: REJECT (|dE| above the class limit of the mode: strict 1e-10, standard 5e-9 Ha; non-reproducible density; an expected counter that does not move in its direction; median wall worse than A by more than the spread), CANDIDATE (median wall better by more than the spread, or an expected counter moved by at least `--min-drop` with the wall inside the spread), NEUTRAL otherwise. Identity is reported as "bit-identical" (energy bits and density hash equal) or "rounding class".
- Calibration (A/A, 2026-11-16, load 4.2 to 4.6): si32 np=4 production 4.0 to 4.4 s per run, spread about 0.4 s (10 percent): simulated-rank wall decides only large effects; the counters decide the screen. Simulated ranks do not model MPI latency (a tree sum is two thread barriers), so latency-only ideas (message pattern) are not screenable here and go to a dev08 microbenchmark first.
- Example: `python3 -I scripts/fast_screen_mpi.py --upf <Si.pz-vbc.UPF> --b BAANDI_M2_PANEL_GS=fused --expect collectives.tree_sum.calls:down` (3 pairs, about 30 s). Unit tests of the verdict logic: `python3 -I scripts/tests/test_fast_screen_mpi.py`.
- Sizing knob of the harness only: `M2_PROFILE_VERIFY=off` (production without the verification, to measure its cost; never an adoption arm).

## Harness: `scripts/fast_screen.py` (stage 1, added 2026-11-16)
Standard library only (`python3 -I`), no network. One case, two arms, N alternating pairs (AB, BA, AB, ...; default 3).
```
python3 -I scripts/fast_screen.py --binary target/release/dft-cli --upf /path/Si.pz-vbc.UPF \
    --env-b BAANDI_X=on                              # environment A/B, same binary
python3 -I scripts/fast_screen.py --binary old/dft-cli --binary-b new/dft-cli --upf ...   # two binaries
python3 -I scripts/fast_screen.py ... --scf-b "davidson_capacity_per_band = 3"            # input A/B
```
- Defaults: `benchmarks/cases/si32-gamma.toml`, seed 601, fractional amplitude 1e-3 (every atom, every component, `random.Random(seed).uniform(-amp, amp)` in atom order, the documented MPI-M1-GATE3 stream; byte equality with dev08's p48lib geometry is not verified, so Mac screen walls are not comparable with dev08 tables), T=1, single-path defaults (no `[scf]` keys beyond seed, threads, max_iterations; `--scf`, `--scf-a`, `--scf-b` add lines).
- Per run: external wall (process), exit code, 1-min load, `converged`, `verified`, total energy, density.bin SHA-256, `timings`, and the deterministic counters: H vectors (recorded + certificate, and each part), coarse-stage H vectors, Davidson iterations (fine + verification, and coarse), SCF steps (fine, coarse), verification attempts, `gemm_dispatch.rows` calls and GFLOP.
- Output: per-arm median, min, max, range/median; paired ratios B/A with median and spread (half range, at least 1 %); energy difference; density SHA identical or not; counters with B/A (the mandatory self-check: a counter zero in both arms prints `unusable (all zero)` and takes no part in the verdict); whether counters are identical within each arm.
- Verdict: REJECT (a failed run, not converged, failed verification, |dE| above `--energy-gate` 5e-9 Ha, or median ratio above 1 + spread); CANDIDATE (median ratio below 1 - spread, or H vectors / Davidson iterations / GEMM GFLOP down by at least 5 % with the wall inside the spread); NEUTRAL otherwise. `--json` keeps every record.
- One si32 T=1 run takes 5 to 6 s on the Mac (2026-11-16, pure build), so a three-pair screen is under a minute. Pick a seed whose verification fails when the idea concerns the verification (si32 f1e-3 seed 603: 4 attempts, 13 fine steps; seed 601: 1 attempt).
