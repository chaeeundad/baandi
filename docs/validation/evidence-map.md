# Claim-to-evidence map

Audit date: 2026-11-05. Checkout `b5d515a`; documentation-only P12A.
Paths in code spans are repository-relative. “Local” describes the original
BAANDI experiment or this comparison replay; it never means QE ran on this Mac.
Remote full output can be parsed locally; a relayed summary is an attestation,
not a verified remote file. Dates below are record dates unless a bundle supplies
an exact UTC timestamp. No new numerical or performance results are invented.

## README accuracy and performance inventory

| Claim (README unless noted) | Stored evidence | How produced; date | Scope and limits |
| --- | --- | --- | --- |
| LDA Gamma: 11 QE iterations, 1.28e-6 Ha/atom energy, 1.27e-6 Ha eigenvalue difference | `tests/reference/reports/si30_ours.json`, `si30_auto.json`; `tests/reference/si30_{ours,auto}/{qe.in,qe.out,provenance.json}` and `si30_ours/ours/results.json` | QE 7.6 dev08, full output; local verifier; 2026-11-02 | One Si cell/UPF, 30 Ry; no overlap matrix comparison |
| LDA 2x2x2: 8 k, 9 iterations, 1.28e-6 Ha/atom, 1.99e-6 Ha | `tests/reference/kpoints/k222_c/{qe.in,qe.in.meta.json,qe.out,ours_results.json}` | QE 7.6 dev08 full output; 2026-11-03 record; replay below | k weights normalized; no density/orbital comparison |
| LDA shifted 2x2x2: 4 k, 9 iterations, 1.28e-6 Ha/atom, 1.77e-6 Ha | `tests/reference/kpoints/k222s_c/` same files | Same provenance class/date | Same limits |
| LDA shifted 4x4x4: 32 k, 8 iterations, 1.28e-6 Ha/atom, 2.10e-6 Ha | `tests/reference/kpoints/k444s_c/{qe.in,qe.in.meta.json,qe_summary.json,ours_results.json}` | Dev08 text relay; 2026-11-03 record | Full qe.out absent; sender digest is not independently verified |
| PBE Gamma: 1 k, 12 iterations, 1.70e-6 Ha/atom, 1.27e-6 Ha | `tests/reference/pbe/gamma/{qe.in,qe.in.meta.json,qe_summary.json,ours_results.json}` | Dev08 text relay; 2026-11-03 record | LDA UPF with explicit PBE override; consistency only, not physically validated PBE |
| PBE shifted 2x2x2: 4 k, 9 iterations, 1.69e-6 Ha/atom, 1.82e-6 Ha | `tests/reference/pbe/k222s/` same files | Same provenance class/date | No density comparison, no full QE output |
| All six rows pass 1e-4 Ha/atom and 1e-3 Ha gates | Above files and `scripts/compare_kpoints_qe.py`, `scripts/verify_qe.py` | Local stored-data replay, 2026-11-05 | Passing comparison does not upgrade relay provenance |
| Seven cutoffs 20–80 Ry; 2.54–2.57 micro-Ha/cell offset and cutoff figure | `tests/reference/reports/scan_e{20,25,40,50,60,80}.json`, Gamma report; corresponding bundle directories | Dev08 full outputs/local BAANDI; R4 2026-11-03 record | Constant within printed precision; one Si cell |
| Offset explained by G=0 integral convention | `docs/validation.md` R4 narrative | Historical local diagnostic, 2026-11-03 | Diagnostic program/results not committed: explanation not independently reproducible from stored files; TODO |
| Checkpoint: 162 memory + 18 binary comparisons; warm max 8.22e-15 Ha; paused audit and cross-process restart | `tests/reference/checkpoint/{results.json,engine.txt,cli.txt,check.txt}` | Local synthetic tests, P11C 2026-11-04 | Commands/source hashes in results; same build, timings excluded; 1e-10 Ha warm tolerance |
| Self-healing: four arms twice; reference/rules/agent 30 iterations, -4.684561184277879 Ha, zero error; baseline fails at 1; applied/rejected 0/0,0/0,1/0,1/1 | `tests/reference/runtime/generated/results.json` and arm subdirectories (`arrays.json`, `audit.jsonl`, summaries) | Local `scripts/demo_self_healing.sh`, 2026-11-04 | Budget repair, not general instability recovery; 1e-10 Ha gate; invalid beta=1.1 rejection |
| PBE closed forms “2e-16 Ha”; LDA unchanged bits | ADR 0018 implementation results; `crates/dft-physics/tests/pbe.rs`; historical snapshot hashes in ADR | Local oracle/snapshot experiment, 2026-11-03 | Maximum is 2.220446049250313e-16 Ha; README rounds. Original snapshot harness/raw output not retained |
| Gamma k-point engine reproduces prior density bits | R6 narrative; `tests/reference/si30_ours/ours/density.bin`; engine k-point tests | Local historical regression, 2026-11-03 | Paired historical density/build records not stored together; test is not evidence of that exact historical pair |
| uom/arithmetic input leave real-Si density bits unchanged | README narrative and historical progress log | Historical local claim | No paired before/after result/hash bundle identified in tests/reference: TODO retain it; do not infer from current tests |
| Fast gates and Python regression tests pass | Current commands below; historical test logs under runtime/checkpoint | Local; current 2026-11-05 checks are scoped explicitly | Current README test has 3 tests, not the entire Python suite |
| Deterministic results “any thread count” | `tests/reference/perf/dev08_threads.json`; final relay Part 1; engine threading tests | Dev08 relay plus local synthetic tests, 2026-11-03–05 | Stored real-Si checks cover T=1,4,8,16; not an empirical proof for every count/platform/build |
| Current speed host: two Xeon Gold 6230, Rocky 9.7, one socket; MKL 2026.1, QE 7.6, fdce016 | `tests/reference/perf/dev08_final_relayed.md`, Parts 1/2/4 | Dev08 relayed, Part 4 2026-11-05 | Pure Rust and MKL BAANDI both shown; no independent reproduction here |
| Shifted 4x4x4 8–16 cores: 2.0–2.2x at 1e-8, 1.4–1.6x at 1e-12; 1 core 4.2–5.1x (1e-8); QE 5 vs 8 iterations | Final relay Part 4 table | Same-session BAANDI median/QE mean, five runs, dev08 | Ratio BAANDI/QE; varies by backend and convergence criterion |
| si8 Gamma 1 core 6.4–8.8x; 16 cores 15–17x at 1e-12 | Final relay Part 4 | Dev08 relayed, 2026-11-05 | One-core range spans both thresholds; README parenthesis is ambiguous. 16-core row is exactly 0.3000 s QE; exclude faster rows |
| si8 threads only 1.4x at eight threads; QE timer 0.01 s, exclude <0.3 s | Final relay Part 4 and protocol | Dev08 relay, 5.801/4.151 s pure Rust | Known weak scaling; timer exclusion is a reporting rule, not an accuracy gate |
| “Kernel measurements promising”, “end-to-end numbers follow” | `docs/performance/fft-backend-comparison.md`, final relay | Historical kernel and later remote SCF results | README wording stale: end-to-end results now stored; MKL is not uniformly faster |
| Buffer-band fix lets 8-atom Si converge | Optimization log P7PC final policy; final relay si8 rows | Local historical diagnosis and remote campaign | One demonstrated failure/fix, not universal robust convergence |

## Validation document inventory beyond README

Each row covers all numerical tables and derived claims in the named section;
protocol thresholds are requirements, not measurements.

| Claim/section | Evidence path | How produced; date | Limits / missing evidence |
| --- | --- | --- | --- |
| R1 energy decomposition, spectrum, charge, degeneracy, gap | `tests/reference/si30_ours/ours/results.json` | Local R1 build; 2026-11-02 | R1 wall time/raw command log not retained in bundle; build/input hashes in result/record |
| R2 cutoff/grid/tolerance tables, timings and relative convergence | Cutoff subset in `tests/reference/scan_e*/ours/results.json`; full table in `docs/validation.md` | Local external `runs/si_scan`, 2026-11-02/03 | Full grid/tolerance scan manifest and raw timing outputs absent from tests/reference; TODO |
| R3 energy/components/all eigenvalues/gap/charge/residual | `tests/reference/reports/si30_{ours,auto}.json` and bundles | Dev08 full outputs/local verifier | Orthogonality not assessed; unexplained-offset wording was historical, later R4 hypothesis |
| R4 automatic-grid equality and diagnostic convention experiment | R4 narrative; explicit-grid bundles only | Dev08 summary/local diagnostic | Automatic outputs and diagnostic code absent; “parallel to 2e-9” corrected to 2.1e-8 Ha from stored energies |
| R5 density rho(0), sum squares, sampled coefficients | `tests/reference/reports/si30_density.json`; `tests/reference/si30_ours/density/qe_rhoG_summary.json` and local density files | Dev08 summary/local comparison, 2026-11-03 | Not all 5985 coefficients; no real full-file density comparison |
| R6/R8 sampled gaps, terms, weights and iteration differences | kpoints/pbe files above | Dev08 full/summary + local results, 2026-11-03 | PBE XC offset unexplained; no converged gap/material claim; no unshifted 4x4x4 reference |
| R7 baseline, O1/O3/O4 updates, cross-platform digits, threshold sweep | `tests/reference/perf/dev08_baseline.json` | Dev08 relay, 2026-11-03 | Historical builds, unequal timing/threshold conventions; remote originals absent |
| R7 k-point threads/P7PA timing, density hashes, ratios | `tests/reference/perf/dev08_threads.json` | Dev08 relay, 2026-11-03 | Superseded current ratios; bit identity only tested cases |
| runtime-equivalence protocol targets / safety threat tables | `docs/validation/{runtime-equivalence,intervention-safety}.md`; engine runtime tests | Design requirements, then executed slices | Proposed security quotas/fault tests are not blanket delivered security claims |
| P11E recovery/safety executed slice | `tests/reference/runtime/generated/` | Local demo, 2026-11-04 | Recorded offline agent; no live LLM/network or general recovery claim |
| P11J all local timing, events, control, memory, checkpoint and historical-baseline tables | `tests/reference/runtime/{benchmark.json,p11j-checks.json,si-baseline-comparison.json,si-pre-refactor.json,si-provenance.json}` and logs | Local paired benchmark, 2026-11-04 | Old Si baseline gate failed; independently built pre-refactor comparison isolates extraction; not a new QE run |
| P11J remote overhead verdict | Final relay Part 3 | Dev08 relay, 2026-11-04/05 | 300-pair verdict provisional; originals missing, does not erase local failures |
| Stagnation table, extended skew/slab, diagnostics, small-gap negatives | `tests/reference/stagnation/{measured,extended-skewed,extended-slab,extended-slab-offset}/`; `inputs/sweep.json`, `verification.json`, `checks.json` | Local synthetic sweep, 2026-11-04 | Complete per-case JSON/input/CSV; negative results retained; no universal recovery defaults |
| Fixed-point multiplicity tables and backend exclusions | `docs/validation/scf-fixed-point-multiplicity.md` only | Local/remote historical narrative, 2026-11-04 | Explicitly no committed raw evidence; syn-s3/z14 unconverged, not proven non-unique; TODO |

## Commands, versions and integrity

R3/R4 bundles retain the exact remote launch in `provenance.json.command`, QE
version/host/date/MPI/OMP and input/UPF digests; use those records rather than
reconstructing commands. R6/R8 do **not** retain a full per-run remote command,
binary digest or UTC manifest. Their reported executable is
`~/fq/builds/orig-mpi-dfti/bin/pw.x`, QE 7.6, MPI x4 on dev08; exact shell command
and launcher flags must be recovered, not invented. Generator arguments are in
`qe.in.meta.json.generator_options`. The common UPF SHA-256 is
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`;
license/provenance remains unverified and no UPF is committed.

The R6 stored `tool.git_commit` is `c31a3d42ef7395f1b8c80f0fc4aaa309f593c8f6`,
not the approximate “dc28615-era” label. R8 stores
`47b2e899061a48dcdb307022e3f39a731c2f7fc1`. Both report dft-cli 0.1.0,
rustc 1.99.1 (58cf8f9c3 2026-10-13); PBE uses xcx 0.3.1 (ADR 0018).
A commit field alone does not establish a clean build. Summary SHA-256 below
checks the local transcription, never the unavailable original qe.out.

Local replay commands (2026-11-05, each exit 0):

```sh
python3 -I scripts/compare_kpoints_qe.py --ours tests/reference/kpoints/k222_c/ours_results.json --qe-out tests/reference/kpoints/k222_c/qe.out
python3 -I scripts/compare_kpoints_qe.py --ours tests/reference/kpoints/k222s_c/ours_results.json --qe-out tests/reference/kpoints/k222s_c/qe.out
python3 -I scripts/compare_kpoints_qe.py --ours tests/reference/kpoints/k444s_c/ours_results.json --qe-summary tests/reference/kpoints/k444s_c/qe_summary.json
python3 -I scripts/compare_kpoints_qe.py --ours tests/reference/pbe/gamma/ours_results.json --qe-summary tests/reference/pbe/gamma/qe_summary.json
python3 -I scripts/compare_kpoints_qe.py --ours tests/reference/pbe/k222s/ours_results.json --qe-summary tests/reference/pbe/k222s/qe_summary.json
scripts/check.sh
python3 -I scripts/tests/test_readme_figures.py
```

The workspace gate passed build/test/fmt/clippy. README regression: 3 tests passed,
with existing unclosed-file ResourceWarnings. No native MKL or QE process ran.
Checkpoint commands/build/source hashes are in `checkpoint/results.json.commands`
and `source_sha256`; runtime benchmark commands and versions are in its README
and JSON; stagnation commands/hashes are in `verification.json`. PBE oracle:
`cargo test --offline -p dft-physics --test pbe -- --nocapture` (historical command,
not separately executed in this audit). Historical density verifier invocation is
in `tests/reference/reports/si30_density.md`/R5. No absent command is marked passing.

### R6/R8 local file integrity (SHA-256, computed 2026-11-05)

All five input and all five result digests match their sidecar fields
`input_sha256` and `source_results_sha256`. Hashes are not producer signatures.

| File | SHA-256 |
| --- | --- |
| `tests/reference/kpoints/k222_c/qe.in` | `40609d7a9ae57fc79e9d720566dbcbf296f4e89a58919a2c637e9cdc3bb33bb8` |
| `tests/reference/kpoints/k222_c/qe.in.meta.json` | `660362f86a04a84e22803bfd580f54ac27a4927b1b6693551b871a7ed9aa1e4d` |
| `tests/reference/kpoints/k222_c/ours_results.json` | `477d25bf7e8c37611316f6fd9fb05a205ac2eefece357df3ec9b8da75feaf417` |
| `tests/reference/kpoints/k222_c/qe.out` | `f93516f6a488f870b5880a9d2f1c1da6ff7e4cde7dd391afdd338948169afaf2` |
| `tests/reference/kpoints/k222s_c/qe.in` | `5a3766085dddada0639c5e2ae5804a217cbc1e8d77e301cb9a8717ad588ed887` |
| `tests/reference/kpoints/k222s_c/qe.in.meta.json` | `13acfea6202893c7a94014d2352e620f0b8a22efff370be9b10c6f8b914b4135` |
| `tests/reference/kpoints/k222s_c/ours_results.json` | `2c84aa7296e9d3b5b137391b4fd95f92ec433ee42d196d3247bbacf6888015b2` |
| `tests/reference/kpoints/k222s_c/qe.out` | `162a89801e5ea6edace18e5e26635a9f9a10e0b84518f3a5316c0a3a21e45ea6` |
| `tests/reference/kpoints/k444s_c/qe.in` | `a45ed650b613c768bdd3ba999ade0941fed5a22325d9a625c69138823e296600` |
| `tests/reference/kpoints/k444s_c/qe.in.meta.json` | `56885b7614965ba0b3f2d4e7788d804f97f58db5972093d8def9b1a14006c7cf` |
| `tests/reference/kpoints/k444s_c/ours_results.json` | `f75c1a3091e40f20041aedf16f036cbf1d16010702a8f5d914e8d5669bdf59ae` |
| `tests/reference/kpoints/k444s_c/qe_summary.json` | `a1b671aee1c20904e4659bdd371cfb82a5f0d103fe76ad2f1e17df70e2e491bd` |
| `tests/reference/pbe/gamma/qe.in` | `7e249ba65e862a0915c91df846c52f94ed99d48f8b11c795f36dd31a099dc034` |
| `tests/reference/pbe/gamma/qe.in.meta.json` | `2cc3ed44f9046842eb3a3adba217153c04c231ed941034938a28ddcf395cdd59` |
| `tests/reference/pbe/gamma/ours_results.json` | `59eb94550327ec971ade840c94c9184373e856e4aa19dcf12213b1952d96fc68` |
| `tests/reference/pbe/gamma/qe_summary.json` | `7a387e456c6f707397d5f3ba49ee72560b0585472256f6c9668699adf8b65170` |
| `tests/reference/pbe/k222s/qe.in` | `d00bee3f9842a220a89a00bb49573acc151f830804a68f6fff8a57acc4f00da5` |
| `tests/reference/pbe/k222s/qe.in.meta.json` | `49507bfd2f7ec56399b7b5c49faa1cd48dc5affd7de00339c02f21f978319d48` |
| `tests/reference/pbe/k222s/ours_results.json` | `faec867ce7e268d98f172e6560cd704d0ea4305d8ced37ad143b342f5fb21b55` |
| `tests/reference/pbe/k222s/qe_summary.json` | `4cb22b2d0cd269ce38856c3210a532a398a2b7ab5cebee2f39521a7aa6b10fc8` |

Remote original hashes: k444s qe.out is recorded in its summary `_note`; PBE
original output and block digests are in `qe_summary.json`. Those cannot be
verified here. Final performance relay Part 4 records pure executable
`f5803c7a1004c7aa465ab8f6f9849d808d415d57df7cd09f00bf7f7c7d550e19`
and MKL executable
`41da71e9dbbd9df470f4e1458d4a5358739794211511d4262901bec621d7db6e`.
Remote protocol commands are in `docs/performance/final-measurement-protocol.md`;
a protocol is not an execution manifest. Part 4's exact shell commands/raw timing
files still need collection. The top-level relay digest refers to the earlier
remote dev08_final.json, not automatically to every later Part 3/4 measurement.

## TODOs requiring owners or missing evidence

- Obtain exact remote launch manifests, full qe.out for k444s/PBE, full R5
  charge-density file, and perf7/perf8/later fdce016 raw timing records. Verify
  sender hashes before promoting any relay to independently checked evidence.
- Retain R2 full grid/tolerance manifest, R4 diagnostic source/output, historical
  Gamma/PBE/uom before/after density records, and multiplicity raw inputs/results.
  Claims above without these records remain explicitly limited, not validated.
- README owner: remove “none of layers B to E exists”, “no runtime code yet”,
  runtime “design in progress”, MKL “end-to-end numbers follow”, and FFT
  “in progress”; update generated `docs/feature-status.json` through its owner.
  PBE oracle bound is rounded, and si8 one-core ratio range spans both thresholds.
  README is intentionally untouched; complete repository-wide prose consistency
  is therefore blocked on this out-of-scope handoff, despite passing generators.
- Performance document owner: `docs/performance/fft-backend-comparison.md` still
  presents stage-1 remote acceptance as pending. Its historical matrix needs a
  current-status annotation; only the performance README/log are owned here.
- Artifact owner: reconcile stale `I-P10B-REMOTE`, session notes and task registry
  with merged commits and provisional remote evidence; generate/upload dashboard
  datasets. P12A cannot mutate those files under its OWNED_FILES. No dashboard
  upload or session registration was performed here; no recurring timer is active.
- Retain unmet worker utilization/scaling goals and rejected CholeskyQR2 as open
  limits; do not convert “merged” into proof every performance target passed.

## Additional stored-data verification — 2026-11-05

R3/R4 full-output replay and input/output-hash verification: every corrected
command below exited 0. These commands write only temporary reports and execute
no QE or SCF binary. An initial si30_auto invocation incorrectly used a nonexistent
`si30_auto/ours/results.json` and exited 1; the corrected invocation uses the shared
si30_ours result, as the stored report does.

```sh
python3 -I scripts/verify_qe.py --ours tests/reference/si30_ours/ours/results.json --bundle tests/reference/si30_ours --json /tmp/p12a-si30_ours.json --markdown /tmp/p12a-si30_ours.md
python3 -I scripts/verify_qe.py --ours tests/reference/si30_ours/ours/results.json --bundle tests/reference/si30_auto --json /tmp/p12a-si30_auto.json --markdown /tmp/p12a-si30_auto.md
python3 -I scripts/verify_qe.py --ours tests/reference/scan_e20/ours/results.json --bundle tests/reference/scan_e20 --json /tmp/p12a-scan_e20.json --markdown /tmp/p12a-scan_e20.md
python3 -I scripts/verify_qe.py --ours tests/reference/scan_e25/ours/results.json --bundle tests/reference/scan_e25 --json /tmp/p12a-scan_e25.json --markdown /tmp/p12a-scan_e25.md
python3 -I scripts/verify_qe.py --ours tests/reference/scan_e40/ours/results.json --bundle tests/reference/scan_e40 --json /tmp/p12a-scan_e40.json --markdown /tmp/p12a-scan_e40.md
python3 -I scripts/verify_qe.py --ours tests/reference/scan_e50/ours/results.json --bundle tests/reference/scan_e50 --json /tmp/p12a-scan_e50.json --markdown /tmp/p12a-scan_e50.md
python3 -I scripts/verify_qe.py --ours tests/reference/scan_e60/ours/results.json --bundle tests/reference/scan_e60 --json /tmp/p12a-scan_e60.json --markdown /tmp/p12a-scan_e60.md
python3 -I scripts/verify_qe.py --ours tests/reference/scan_e80/ours/results.json --bundle tests/reference/scan_e80 --json /tmp/p12a-scan_e80.json --markdown /tmp/p12a-scan_e80.md
python3 -I scripts/compare_rhog.py --ours-dir tests/reference/si30_ours/ours --qe-summary tests/reference/si30_ours/density/qe_rhoG_summary.json --json-report /tmp/p12a-density.json --markdown-report /tmp/p12a-density.md
```

Density summary replay also exits 0: largest sampled coefficient difference
4.008255174783675e-9 electrons/Bohr^3; full-set verification remains false.
Subtracting the R4 `total energy per atom` checks, multiplying by two atoms and
subtracting the 80 Ry offset gives maximum 2.0414727686102196e-8 Ha, supporting
the corrected rounded upper bound of 2.1e-8 Ha. This is arithmetic from stored
reports, not a new reference run.
