# AGENTS.md: rust-dft (read this first)

Rust-first plane-wave DFT engine. Goal of the first release: a **verifiable** Si Gamma-only LDA SCF, not speed. This repo is AI-first: future changes are made by AI agents, so every contract below is written to be read and checked by machines and humans.

## Read order for any task
1. This file. 2. `docs/INDEX.md` (map of all docs). 3. `ownership.toml` (who may touch what, crate dependencies). 4. `crates/<name>/AGENTS.md` for each crate you edit. 5. Your ticket (`tickets/<TASK_ID>.md`). 6. Relevant ADRs in `docs/decisions/`.

## Roles
Claude Code = architect, orchestrator, reviewer. Codex = implementer. The user = final approver. A completion claim is never trusted: run `scripts/check.sh` and report real exit codes and numbers. Never report an untested command as passing.

## One-command verification
`scripts/check.sh [crate]` runs build, test, fmt, clippy (`-D warnings`) and prints PASS/FAIL per gate. Toolchain pinned in `rust-toolchain.toml`; add `/opt/homebrew/opt/rustup/bin` to PATH if cargo is missing.

## Source of truth for state
`artifacts/*.json` (milestones, issues, tests, agents, sessions), `docs/progress.md`, `docs/decisions/` (ADRs), tickets, and git commits. Update the JSON files when a gate is verified. Every Codex run is recorded in `artifacts/sessions.json` (PID, role, session id, status); use `scripts/delegate.sh`, and stop runs only with `scripts/delegate.sh stop <TASK_ID>`, never by process name.

## Physics and numerics conventions (details: `docs/physics-conventions.md`)
- Hartree atomic units internally (Bohr, Hartree). 1 Ha = 2 Ry. Convert only at I/O boundaries.
- Lattice vectors are the columns of A; B = 2*pi*inverse(A)^T; psi(r) = Omega^-1/2 * sum_G c(G) exp(iG.r), sum |c|^2 = 1.
- Gamma-only complex calculation, spin-unpolarized, norm-conserving, LDA. Do not relax tolerances to make a test pass; trace mismatches in this order: units, FFT normalization, G=0, UPF projectors, double counting, XC parameterization.

## Library policy (ADR 0003)
Adopted so far (ADR 0004): faer, ndrustfft/rustfft, ndarray, num-complex. Prefer maintained general-purpose crates over hand-written numerics or FFI (linear algebra, FFT, arrays, serde/toml, XC cross-checks). Hand-write only physics glue. New crates need license, maintenance and Apple Silicon checks and an entry in `docs/dependencies.md`; the Codex sandbox has no network, so Claude fetches crates and adds them to manifests.

## Hard rules
- Keep `scripts/check.sh` fast: a test that takes more than about 10 s in debug must be `#[ignore]` and documented with its release-mode command (ADR 0004).
- MVP approved 2026-11-03; P7 in progress (see `docs/p7-plan.md`). Still out of scope unless a P7 ticket says otherwise: CUDA, Metal, MPI, PAW, USPP, forces, stress, smearing. In scope only through the named P7 tickets: multiple k-points (P7A), PBE through a library (P7B), GPU backend design note (P7C, no GPU code).
- Scope decision U1 (user, 2026-11-05): option B, staged extension beyond the silicon baseline: collinear spin, metallic occupations (smearing), DFT+U, symmetry, forces, stress, relaxation, in the order of `docs/roadmap/master-plan-proposal.md` (MF02 to MF18), each stage through its own ticket, its own approved design contracts and independent validation oracles. U2 to U8 default to the proposal's recommendations until the user says otherwise. PAW/USPP, GPU, hybrids and SOC stay out of scope until a ticket and a measured need say otherwise. Scope decision U9 (user, 2026-11-09): MPI is approved as the next parallel layer, staged through `docs/decisions/0041-mpi-parallel-layer.md` (design contracts, determinism contract and pre-registered scaling gates); it lives behind a build feature, the threads-only default path stays bit-identical, and every stage is validated on dev08 through the master session. Scope decision U10 (user, 2026-11-10): MPI stage M2 (plane-wave/G-vector distributed 3D FFT with transposes and distribution of the replicated stages) is approved, to be planned in `docs/roadmap/mpi-m2-plan.md` with pre-registered gates before code; in parallel a separate sub-session works on single-node speed (threads-only path) in its own worktree. The opt-in `BAANDI_MKL_NONLOCAL` stays opt-in (user approved the recommendation, 2026-11-10).
- Never copy QE (GPL) source. Never commit UPF files. Never record numbers as "QE reference" unless they were produced by a real run with command, version, input and UPF SHA-256.
- QE reference runs are executed remotely (dev08 via the dft5090 session). Do not install pw.x on this Mac.
- No JSON, logging or file I/O in numerical hot loops. `unsafe` only in FFI adapters, with preconditions documented and tested.
- Stay inside `OWNED_FILES` of your ticket; `ownership.toml` lists frozen files.
- The sandbox cannot write git metadata for task worktrees: leave changes uncommitted and say so; Claude commits after review.
- No remote server changes, bulk installs or destructive git operations.

## Language
All agent-facing artifacts (docs, tickets, reports, code comments, commit messages) are in English.

## Reporting rules to the user (user rule, 2026-11-03, apply for the whole session)
The user reads Korean; these rules govern the final messages of the orchestrating agent. Silence from the user means "carry on": a task is considered in progress until reported otherwise, and nobody should ask whether it is still running.
0. (added 2026-11-03) Every report shows the tasks in progress at the very top, each with its ID and a plain Korean description of what it does (descriptions live in `artifacts/task_titles.json`; add an entry whenever a task is created), and every report is accompanied by an update of the dashboard artifact (`scripts/gen_agents_json.py`, then upload the four `artifacts/*.json` files and update the datasets).
1. A decision needed from the user is written as a `[승인필요]` block at the very BOTTOM of the message (user rule, 2026-11-13; it was at the top before), after the `이 세션의 예정사항` line, with the options and ONE recommendation. Nothing else is asked of the user. If no decision is needed the block is omitted.
2. Every other report is per milestone and has three lines: done / in progress / next.
3. A task expected to take more than 30 minutes gets an estimated duration when it starts, and its progress is updated on the dashboard.
4. Every report ends with `이 세션의 예정사항: 없음` or `이 세션의 예정사항: 있음 (what)` (user rule, 2026-11-13: renamed from "지금 내가 할 일" because "나" was unclear; it means what this orchestrating session itself will do next, never a task for the user, which goes in the `[승인필요]` line).
5. (added 2026-11-04, periodic reporting) While any task is running the orchestrating agent reports on a fixed rhythm without being asked: (a) at least every 30 minutes, (b) immediately when a task finishes, is blocked, or fails, (c) at the start of every session or when the user returns. Before each report it reads `artifacts/sessions.json` and the task's `.delegate/<TASK>/report.json` and the worktree status; a task is called "in progress" only if its registry status is `running`. A `blocked` or `failed` run is reported at once with the reason and either the fix already applied or a `[승인필요]`. The orchestrator keeps a recurring session timer for this and says so if the timer is not active. Each periodic report is the normal short report (rule 0 to 4) and updates the dashboard.
Keep reports short; put detail in the repository documents and link them.

## README results (user rule, 2026-11-03)
After every verified implementation, update the "Validation results" block of README.md: add the key test results (real numbers from stored evidence, with the reference code, scope and limits stated) and a figure where it helps. The block and the figures in `docs/img/` are generated by `python3 -I scripts/make_readme_figures.py` from `tests/reference`; extend that script for a new result type instead of typing numbers into README.md. `scripts/tests/test_readme_figures.py` fails when the README block or a figure is stale, so run the script and commit its output with the change that produced the evidence. Do not edit the rest of README.md (the user maintains it) except for a factual correction that you report.

## QE source study policy (user decision, 2026-11-03)
The original rule "never copy QE (GPL) source" stays. On 2026-11-03 the user also allowed reading the Quantum ESPRESSO source to find performance tuning points. Safeguards, because this project is MIT licensed:
- The QE 7.6 source is cloned OUTSIDE this repository (`/Users/mino/Develop/rust-dft-wt/qe-src`, shallow clone of tag qe-7.6); it is never vendored, committed, copied into a worktree or pasted into a ticket.
- Only Claude Code reads it. It writes analysis notes (`docs/research/qe-study/*.md`) in its own words: algorithms, parameters and their reasons as prose and mathematics. A note contains no code, no source excerpts, no comments, no identifiers beyond file and routine names used as references, no data tables copied from QE, and no structural mirroring that would amount to a transcription. Codex and other implementers receive only the notes and the existing rust-dft code; they must not be given QE source and must not clone it.
- Every change derived from a note is marked in `docs/performance/optimization-log.md` ("derived from note X") and is implemented against the engine's own design, validated by its own tests and by the existing QE output comparisons. Numerical choices (thresholds, subspace sizes) are justified by measurement on our code, not by "QE does it".
- Evidence independence: the QE comparisons recorded before this date (R1 to R8, README) were produced by an engine written without reading QE source and stay valid as independent evidence; statements in the docs that the QE source was not read are historical and are annotated, not deleted. Results of code changed after a note are still compared against QE outputs, and the notes list what was learned so that the influence is visible.
- If a note would require reproducing non-obvious expression from QE (for example a specific table or formula layout) it is not written: describe the principle instead.

## Generality rule (user, 2026-11-03)
BAANDI is a general-purpose engine, not a silicon code. Performance changes, default parameters (tolerances, schedules, subspace sizes, thresholds, mixing parameters) and algorithm choices must not be tuned or accepted on silicon alone. Evidence for adopting a default or an optimization must come from a set of different systems: cell shapes (cubic, skewed, elongated, slab-like), atom counts, k meshes and Gamma, band counts (small and large), pseudopotentials with different projector structures (synthetic fixtures when real files are missing), insulating gaps from small to large. A result that holds on silicon only is recorded as such and the default stays unchanged. Algorithms that matter for larger or different systems (for example RMM-DIIS, preconditioned mixing) are developed as isolated, switchable modules even when they do not pay off on the small benchmark cases. The benchmark set lives in ticket P10G and `docs/performance/`.

## Disk hygiene (user request, 2026-11-03)
Each task worktree builds its own `target/` (1 to 6 GB; the worktree folder had grown to 79 GB and filled the disk). After a task is merged or rejected run `python3 -I scripts/wt_cleanup.py` (it removes `../rust-dft-wt/<TASK>/target` only for tasks whose latest run is not running; `--dry-run` lists first). Never delete anything else in `rust-dft-wt` (vendor, pseudo, qe-src, runs) without the user's say.

## Output formats policy (user rule, 2026-11-03)
Scientific outputs are written in standard, widely readable formats first, and in the established codes' formats where a de-facto convention exists, so that existing tools and workflows work without conversion; the project's own JSON/CSV is the complete record, not the only format. Priorities:
- Volumetric data (density, ELF, potentials, partial charge): Gaussian cube and XSF first; VASP CHGCAR/LOCPOT style and VTK image data (`.vti`, ParaView) as additional writers; raw binary only as an internal exchange.
- Band structures and DOS: CSV with a documented header and plain `.dat` columns readable by gnuplot, xmgrace `.agr`, a QE-style `bands.dat.gnu` / `dos.dat` layout, and a VASP-style `EIGENVAL` / `DOSCAR` layout; JSON metadata always; an SVG figure for people who do not code.
- Structures: POSCAR, CIF, XYZ/extxyz, XSF for input and output when structure files are added.
- Every writer is validated by an independent reader or a documented format specification test (not only by the writer's own reader), states units and conventions (energy reference, origin, ordering, atomic numbers), and never invents data a format cannot carry (the document says what is lost).
- Writers live in the CLI/adapter layer; the numerical core does not know file formats (ADR 0017).
Formats of other codes are written from public format documentation only; no source of those codes is copied.

## Performance owner (user request, 2026-11-05)
Performance work has a standing owner role, the `perf-optimizer` agent (`.claude/agents/perf-optimizer.md`, assigned by the orchestrator). Every agent that touches a performance question (Codex implementers and Claude agents) reads `docs/performance/AGENT-BRIEF.md`, `docs/performance/attempts-ledger.md` (every optimization tried, adopted or rejected, with numbers and the reason) and `docs/performance/backlog.md` (ideas not yet tried) before starting, and never repeats a rejected idea without stating the new evidence that changes its premise. After finishing, add the attempt with its numbers and lesson to the ledger. Acceptance numbers come from dev08 through the orchestrator; Mac timings are observations.
Performance target (user, 2026-11-05): BAANDI must match and then exceed Quantum ESPRESSO 7.6 at equal core count and equal accuracy on the benchmark set (`docs/performance/benchmark-suite.md`, dev08), compared against QE's BEST configuration for each case (for Gamma-only runs QE's `K_POINTS gamma` real-wavefunction path, not only the complex Gamma input used so far). The user re-confirmed on 2026-11-05 that the QE source may be studied for everything useful (policy in the next section): Claude Code agents (never Codex) read the QE source and write analysis notes in their own words under `docs/research/qe-study/`; implementers receive only the notes.

## Goal and path decisions of 2026-11-14 (user; details: `docs/performance/requests/work-order-beat-qe-2026-11-14.md`)
1. One calculation path: the standard path (QE-style stopping rule and eigenvalue criterion; accuracy is a threshold) becomes the only SCF path. The fixed-point certificate becomes a verification step after convergence, ON by default in production runs (its cost is inside every timed BAANDI run and reported as its own table line); `verify = off` exists. The present strict path stays as a test and reference tool; accuracy gates keep using it as the reference.
2. Target: BAANDI wall at most 0.90 times the fastest QE configuration (P0.3) for each cell and each parallel model.
3. The target covers three models, each compared only with the same model: one core, BAANDI threads against QE OpenMP at the same thread count, BAANDI MPI against QE MPI at the same rank count. No table, ratio, argument or conclusion may put one model against another, not even as information.
4. Every QE against BAANDI comparison uses identical conditions: cell and coordinates, cutoff, FFT grid, pseudopotential file, band count, start, threshold, cores and binding.
5. GPU stays out of scope. MPI development is reopened for the target (work order phase E2) and the earlier freeze at M2d is lifted for that purpose.

## Speed comparison rule (user, 2026-11-14)
Speed comparisons use the same parallel model on both sides: multithreading against multithreading (BAANDI threads against QE OpenMP, one rank), MPI against MPI (BAANDI `--features mpi` against QE MPI at the same rank count). Never compare BAANDI threads with QE MPI ranks in a headline, ratio, conclusion or informational column (rule tightened 2026-11-14: a mixed column may not appear at all). One core may always be compared. Details and the list of affected tables: `docs/performance/comparison-rules.md`. Every report, table, ledger row and dev08 request involving QE follows it.

## Performance work (user rule, 2026-11-03)
Read `docs/performance/README.md` before changing any hot path (Hamiltonian apply, eigensolver, FFT, density). Measure and profile first, prove equivalence with `scripts/bench_scf.py` against `docs/performance/baseline.json` (energy within 1e-12 Ha; bit-identical density for overhead-only changes), keep allocation and logging out of hot loops, and record every accepted optimization (profile, change, before and after numbers, rejected ideas) in `docs/performance/optimization-log.md`. New code on the hot path is written to those rules from the start.

## Report format
Final message follows `docs/schemas/report.schema.json`: status, changed files, commands with exit codes, numerical outputs with units, limitations, open questions, commit.

## Adoption rule (user, 2026-11-10; replaces the per-change percentage thresholds of earlier pre-registrations)
Be aggressive: every change that measurably improves speed on dev08 is adopted, even by 1 to 3 percent, and improvements are collected and stacked. Nothing is rejected for a small gain. A change is rejected only for (a) physical hard-coding (a constant, threshold or branch tuned to one system or material instead of a general rule), or (b) a correctness failure (accuracy classes strict 1e-10 Ha and standard 5e-9 Ha, certificate, convergence, thread-count independence of results). dev08 measurement is the priority: implement first, measure on dev08 through the master session at once, adopt as the default when the median wall improves and nothing regresses beyond run-to-run noise. Defaults stay switchable by an environment variable for A/B checks and rollback. Per-change percentage thresholds in the plans and requests under docs/performance/requests and docs/roadmap are superseded by this rule; their correctness and identity clauses stand.

## Decisions of 2026-11-11 (user, answers to the work-order questions)
1. Strict gate: the strict ACCURACY CLASS is 1e-10 Ha (the number that decides adoption); cross-path agreement on unique-solution fixtures is checked at 1e-12 Ha as an internal regression check and reported, not used as an adoption gate.
2. Headline protocol of the QE comparison: after the displacement-amplitude scan (work order M2) the headline column uses a realistic displacement amplitude; the +-1e-7 Bohr column stays as a stress test. The ten-geometry protocol (seeds 801-810 at +-1e-12 Bohr, never repeated runs of one geometry) is the unperturbed-equivalent protocol.
3. Standard class estimator: match QE. If the equal-delivered-accuracy measurement (work order M6) shows that the standard estimator stops later than QE at equal delivered true error, the standard threshold/estimator is recalibrated so that BAANDI's delivered accuracy equals QE's (conv_thr 1e-8 equivalent), documented in an ADR with the measurement, and the QE comparison table states the definition.

