# Fast optimization loop: screen on one case, confirm on the set (user request, 2026-11-16)

Running the whole generality set for every idea is expensive (hours of dev08 time, most of it spent on ideas that do not pay). The loop is split into three stages. The Generality rule and the Adoption rule (AGENTS.md) are unchanged: a default is adopted only after stage 3.

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
