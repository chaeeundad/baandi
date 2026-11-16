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
