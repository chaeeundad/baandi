# Performance work: agent brief (read before any hot-path change)

Applies to Codex implementers and Claude agents. Rules from AGENTS.md (Performance work, Generality rule, QE source study policy) stay in force; this is the operational summary.

## Read order
1. AGENTS.md. 2. This file. 3. [attempts-ledger.md](attempts-ledger.md) (what was tried, section C pitfalls, section D bottleneck map). 4. [backlog.md](backlog.md). 5. [README.md](README.md) rules and tools. 6. Your ticket. 7. [optimization-log.md](optimization-log.md) sections of the rows you touch. 8. tests/reference/perf/dev08_final_relayed.md (parts 4 and 5 are the current dev08 state). QE: study notes in docs/research/qe-study/ only, never QE source.

## Process (in this order)
1. **Check the ledger.** Search it for your idea (ID, keyword, mechanism). If a row says rejected, superseded or gated, stop and write the new evidence that changes its premise into your ticket report before any code. "Different machine" or "different system" counts only if the row's mechanism depends on it.
2. **Measure first.** Use a profile of the current main on the machine that matters (dev08 for threading and MKL). The Codex sandbox cannot run `sample` (exit 255): ask the orchestrator for a profile, or use counters (`worker-analysis`, `P14A_WORK_ACCOUNTING=1`, `P14B_GEMM_ACCOUNTING=1`, source-only cfg builds, never used for acceptance timing). A profile ratio from fewer than 30 samples is directional only.
3. **One change at a time.** One mechanism per branch. Ablate (validation, resident branch, pairing...) when a change has several parts; keep each part only if its own number is positive.
4. **Work, not time, first.** Count operator applications, Davidson steps, GEMM 2mnk, copied bytes, region counts before and after; these are exact on any machine. Time is the second, noisier evidence.
5. **Three-point comparison.** Build and measure (a) the baseline commit, (b) current main, (c) your branch, each in its own target directory, binary sha256 recorded. The effect is (c) against (b); (b) against (a) shows what else moved.
6. **Mac is an observation.** Alternating pairs: at least 30 per condition for a 1-3% claim (7 pairs only to settle one outlier); load recorded; no builds overlapping. Mac scaling is not evidence for dev08 (see ledger C.1). Never claim acceptance from the Mac.
7. **dev08 acceptance goes through the orchestrator**, who asks the remote measurement session. You do not log into dev08. Give them the request in the ticket report (below).
8. **Update ledger and backlog** with new rows (numbers, mechanism, lesson, source), including negative results.

## Gates (every change)
- `scripts/check.sh` exit 0 (build, test, fmt, clippy -D warnings); tests longer than 10 s in debug are `#[ignore]` with the release command documented.
- `scripts/bench_scf.py` cases gamma, k222s, k222, k444s, si8, T=1,2,4,8,12, `--require-identical` against a same-HEAD build. Overhead-only change: energy and density bit-identical. Algorithm or rounding change: energy within 1e-12 Ha, SCF count and fixed-point certificate unchanged for realistic fixtures, density bit-identical across thread counts within each backend (pure and MKL separately), and the change must be stated and approved as a rounding change (P14B rule: converged toy fixtures keep the energy tolerance, iteration changes are reported, never a new convergence failure).
- Generality set (`scripts/bench_suite.py`, ticket P10G): cubic, skewed, elongated and slab cells, 4/8/20 atoms, Gamma and k meshes, synthetic projectors; no default accepted from Si alone. Toy iteration counts are unreliable (fixed-point multiplicity, ledger C.11); judge them by energy tolerance and convergence.
- No tolerance relaxed, no certificate weakened, no allocation, logging or I/O in numerical hot loops, `unsafe` only in FFI adapters. `scripts/make_readme_figures.py` and `scripts/tests/test_readme_figures.py` pass (README results block).

## Do-not-repeat list (short form; the ledger has mechanisms)
Capacity 2x/3x/4x/32/64 and block 4 (A08, A13); fuller preconditioner (A14); Ritz-value-change test with strict refinement (A15); multiplied SCF tolerance factors (A21); 1 or 2 extra bands as default (A19); FFT lane batching (A11) and packed-lane FFT threads (A34); scoped OS threads per GEMM or H block (A33, A35); leaf caps (A46); finer GEMM subtiles (A49); CholeskyQR2 (A52); blocked Nonlocal traversal (A43); band-grouped or many-band 3D FFT (A28); hybrid RMM as default (A22). Ideas not yet tried: backlog.md.

## Optimization-log entry (docs/performance/optimization-log.md)
Newest first at the top. Include: ticket and commit hashes (baseline, main, branch); profile or count that motivated it (with sample counts); the change in a short paragraph; table of before and after numbers per case and thread count (machine, load, pairs, medians, raw samples or where stored); equivalence evidence (energy deltas in Ha, density sha256, SCF counts, certificates, QE comparison script results with unchanged tolerances); generality results; rejected ideas and failed attempts with numbers and exit codes (never hide a failing run); exact commands; limits (what was not measured, relayed versus reproduced). Then add rows to the ledger. Derived changes from a QE study note are marked "derived from note X".

## Requesting dev08 measurements (ticket report, section "dev08 request")
Ask for: (1) the three commits and the exact branch HEAD; build commands for pure and `--features mkl` with `MKLROOT=/opt/intel/oneapi/mkl/2026.1`, binary sha256 for each of the three builds, separate target directories, build exit codes; (2) protocol: 5 cases (gamma, k222s, k222, k444s, si8) x pure, MKL x T=1,4,8,16, `taskset -c 0-(T-1)` on one socket, `MKL_NUM_THREADS=1`, node idle (load below 0.30 before each run), uptime before and after, 3-run medians (5 for a close call), spread reported, any batch above 5% repeated; (3) profiles: `perf`/eu-stack with at least 30 samples on the function you changed, caller-resolved, for si8 T=1 and T=8 (state samples per run); (4) work counters: SCF iterations, Davidson steps, operator applications per case; (5) accuracy: energy pure vs MKL (gate 1e-12 Ha), density hash identical across T per backend, energy of branch vs main; (6) for QE ratios: same node, pw.x and input hashes as part 4, conv_thr 1e-8 and 1e-12, QE internal time, entries below 0.3 s excluded (0.01 s resolution); (7) the original result files stored on dev08 with sha256, so the relay can be verified later.

## Rules about relayed numbers
- dev08 numbers reach the repository as text. Mark them "relayed" everywhere, with date, commit and binary hash; never write them as measured here or into a README claim beyond what dev08_final_relayed.md part 4 supports.
- Never quote a profile-derived ratio below about 30 samples; state sample counts.
- Do not mix denominators: QE steps are per block, ours per vector; QE ethr is a Ritz change in Ry, ours a residual in Ha; QE timer resolution is 0.01 s; a conv_thr must accompany every QE ratio.
- A relayed number with a transcription problem or without the source file is corrected or retracted in the same document when the original arrives (as with the NonlocalOperator::add reference, 0.323 s corrected to 0.474 s).
- Build hash embeds the git commit: identical code can have different hashes; record source hashes.

## Reporting
Final report to the orchestrator follows docs/schemas/report.schema.json: changed files, commands with real exit codes, numbers with units, limitations, what the orchestrator must ask dev08. Leave changes uncommitted (the sandbox cannot write git metadata); the orchestrator reviews and commits.
