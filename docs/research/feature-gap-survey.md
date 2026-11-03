# Practical feature and workflow gap survey

| Priority | Researcher outcome | Verified current position | Recommendation |
|---|---|---|---|
| P0 / T1 | Trust a run beyond the smallest demonstration | Si 8-atom default run fails here; the fast gates do not exercise that real input | Complete P7PC with a maintained failure/regression matrix |
| P0 / T2 | Know what material and settings were actually calculated | Species are aliases; Fe can silently select the Si potential; `validate` accepts invalid SCF settings | Add full preflight, explicit element identity, actionable errors |
| P1 / T3 | Open results in a spreadsheet or viewer | JSON + raw density exist; bands/DOS and cube/XSF are tickets, not commands in this checkout | Deliver P8A/P8B, plus small CSV exports and recipes |
| P1 / T4 | Run oxides, cathodes and alloys with auditable inputs | One real Si potential in validation; UPF v2 reader has a three-projector ceiling | Curated metadata registry and generator compatibility work |
| P1 / T5 | Reuse results safely in automated workflows | Versioned JSON exists, but no formal result schema or transactional output bundle | Stabilize contracts and completion manifests |
| P1 / T6 | Prepare structures without writing coordinates by hand | Only TOML lattice/fractional atoms | External conversion adapter; ASE interoperability outside core |
| P1 / T7 | Choose converged cutoff and k sampling | Cutoff/FFT/tolerance scan script exists; no k scan or spacing input | Extend scans and expose accessible commands |
| P1 / T8 | Resume interrupted, costly supercells | No CLI restart/checkpoint reader | Align density restart with ESLE phase B |
| P1 / T9 | Catch regressions introduced by AI changes | Strong local unit/numerical tooling; no tracked hosted-CI workflow/fuzz harness | Tiered CI, mutation corpus and numerical stress cases |
| P2 / T10 | Install and obtain a first plot without a developer environment | Source build uses a sibling vendor tree; no tracked binary/container release recipe | Reproducible release package, onboarding and citation metadata |

Surveyed checkout: `1a9e20c61e6d5c0c0ed61aa0da13ed8902d013ac`, 2026-11-03. This is a proposal, not implemented functionality. Only this document is changed. Build products, generated inputs and experiment outputs are under `/tmp`. No Quantum ESPRESSO source was accessed; no real UPF was copied or modified.

## Scope, priorities and evidence method

The target applications are transition-metal oxides, battery cathodes, defects, alloys and magnetic materials, as stated in [P8D](../../tickets/P8D-XC-ROADMAP.md). Practical features help those users but do not remove the present spin/occupation/pseudopotential physics restrictions. Si is the current real-material test vehicle, not evidence that these target classes work. Physics expansion belongs to P8D; bands and volumetric implementation belong to [P8A](../../tickets/P8A-BANDS-DOS.md) and [P8B](../../tickets/P8B-VOLUMETRIC.md). Avoid creating competing implementations.

P0 means correctness or misleading-success risk before wider use; P1 means the next usable research workflow; P2 means subsequent scale or specialist work. Estimates below are engineering judgments in **agent-days**, including implementation, tests and documentation, excluding human review, dependency procurement and remote reference queues. They are not measured schedules. High risk means numerical, compatibility or scientific interpretation risk, not a demonstrated exploit. Values are specific to the targets above. Overlapping tickets must not have their estimates added as independent work.

Repository claims use the concrete implementation and tests listed below. Documentation and tickets establish intentions only. Absence claims mean no exposed implementation was found in this checkout, not that a remote branch or unpublished artifact cannot exist. Review commands included:

```sh
git rev-parse HEAD
git status --short
rg --files --hidden -g '!target/**' -g '!.git'
git ls-files '*CITATION*' '*Docker*' '*Container*' '.github/*' '.gitlab*' '*fuzz*' '*proptest*' '*pyproject*'
rg -n 'pub struct|enum |fn |deny_unknown|threads|memory|schema_version|println|write|element' crates/dft-cli/src/{config,run}.rs
rg -n 'restart|checkpoint|initial_density|initial_orbit|pub fn|struct ScfResult' crates/dft-engine/src/scf.rs
```

The `git ls-files` query returned no paths (exit 0). Full file enumeration, the complete CLI command enum and input structs were inspected, rather than treating a failed keyword search as proof. Files named below are repository-relative; symbols are given where more useful than transient line numbers.

## Verified inventory

| Area | What exists, with implementation evidence | Practical gap / boundary |
|---|---|---|
| Structure formats | `crates/dft-cli/src/main.rs::Command` has `eval`, `convert`, `run`, `inspect-upf`, `validate`, `grid`. `config.rs::System/Atom` reads TOML lattice vectors and fractional positions. | No CIF, POSCAR, XYZ or extxyz import command/reader. No Cartesian atom input field. Periodicity and cell requirements need explicit handling in an external importer; XYZ alone must not invent a cell. |
| Units and expressions | `config.rs::normalize`, `crates/dft-units/src/{lib,expr}.rs`, `crates/dft-units/tests/{units,expressions}.rs`: legacy numeric fields, unit strings and arithmetic expressions; normalization at the boundary. `run.rs` records expressions when present. | Good foundation, not a missing units system. Legacy cutoff is Ry, tolerances Ha; input rows are lattice vectors but result matrix vectors are columns. Tutorials/export headers must make this distinction visible. |
| Pseudopotentials | `crates/dft-upf/src/parser.rs::parse_upf/from_file`: validated UPF v2 NC subset, radial arrays, functional labels and typed errors. `run.rs::execute` hashes original bytes, resolves paths relative to input and records valence/functional. `inspect-upf` displays element/projectors. | No validated multi-element library, license registry or cutoff recommendation service. Candidate Si provenance/license remain unresolved in `docs/candidate-si-upf.md`. No packaged UPFs. Hash identity is not quality or licensing validation. |
| UPF v1 / ONCV | `parser.rs` version checks accept 2.0/2.0.0/2.0.1; reject >3 projectors. `docs/upf-support.md` specifies nominal projector units, no NLCC/PAW/USPP/SOC and required cutoff metadata. | UPF v1 unsupported. ONCV is a generator family, not a single extra format switch: some UPF2 files may meet this subset, others exceed projector limits or require rejected features. **No real ONCV file was tested here; blanket compatibility is unestablished.** |
| k points | `config.rs::Basis/resolve_kpoints`, `crates/dft-geometry/src/kpoints.rs`: explicit positive normalized weights, integer meshes, shifts 0/0.5, time reversal. | No reciprocal spacing or density selector. No space-group discovery/reduction or density symmetrization; time reversal already exists and must not be described as no symmetry at all. |
| Convergence studies | `scripts/convergence_scan.py::run_scan` emits `scan.csv`, `scan.json`, `scan.md`, manifest and optional QE inputs for `ecut`, `grid`, `tolerance`, `all`; tests in `scripts/tests/test_scan.py`. | `grid` means **FFT grid**, not k mesh. No built-in CLI scan command or k convergence scan. Script edits only supported single-line TOML values. Its gap calculation uses the first-k eigenvalue array, not a general sampled indirect gap. |
| SCF and numerical contract | `crates/dft-engine/src/scf.rs`: LDA/PBE functional selection, Pulay/linear mixing, scheduled eigen tolerance, energy/charge/density/eigen checks and fixed-point certificate. `run_outcome` retains budget-exhausted state. | Spin-unpolarized fixed integer occupations; odd electron counts rejected. Gap certification applies to new k input, with legacy Gamma compatibility explicitly retained in ADR 0013. Default Davidson behavior still fails the Si8 experiment below. |
| Restart / checkpoint | `ScfResult` exposes density and orbitals in memory; `step/step_kpoints` exist. CLI writes density plus metadata. `docs/roadmap/esle.md` phase B proposes `run_from`, save/load and compatibility checks. | Existing density files are outputs, not restart support. `Command`, `Config`, and engine `run_outcome` have no persisted restart path, periodic checkpoint or resume contract. Do not promise exact continuation from density alone. |
| Progress / logging | `run.rs::execute` prints XC selection before solving; `main.rs` writes returned console after execution. Iteration history/timings are assembled into JSON and final console after SCF. | No live iteration progress, structured event stream or cancellation checkpoint. A long solve can be quiet after its initial XC line. Numerical failure returns without new result history. |
| Machine-readable output | `run.rs` writes `results.json` and `density.json` with `schema_version=1`, density binary hash, energy audit, eigenvalues, occupations, geometry and convergence. `docs/cli.md` describes units/layout. ADR 0013 adds k data without breaking legacy keys. | No result JSON Schema file: `docs/schemas/` contains report schema only. Optional keys vary with legacy/new k input and functional overrides. Existing top-level eigenvalues and orbital memory are first-k views. A consumer must not interpret them as all-k results. |
| Filesystem robustness | `run.rs::protect_outputs` checks input aliases, including Unix hard links; tests in `crates/dft-cli/tests/run.rs`. Writes use `fs::write` separately. | No atomic multi-file commit; interrupted writes can mix generations. Failure before writing leaves old outputs. Readers need exit status and identity checks; directory existence is not success. |
| Errors / validation | `Config` tables use `deny_unknown_fields`; geometry checks finiteness and scaled determinant; `run.rs::engine_error` maps typed failures to exit classes. Missing species entries and UPF paths are diagnosed. | `validate` checks geometry/basis but does not open UPFs or validate SCF numerical options. Species alias is not checked against UPF element, and output species omit that element. Low-cutoff band incompatibility is wrapped as numerical exit 4. See experiments. |
| Threads / scaling | `ScfOptions::threads`, CLI `scf.threads`, `crates/dft-engine/src/{scf,band_density}.rs` support scoped CPU workers over k points/bands; `tests/threading.rs` compares numerical bits across layouts. `run.rs` records requested/outer/inner threads. | No MPI/GPU controls. `crates/dft-mkl` and ADR 0019 exist, but `crates/dft-numerics/Cargo.toml` has no MKL feature/dependency in this snapshot: this CLI build is not an enabled MKL backend. No current evidence of large distributed scaling. |
| Memory / resource planning | `grid` reports dimensions/basis; `run.rs` reports array payload bytes and null peak RSS; `scripts/benchmark.sh` uses macOS `/usr/bin/time -l`. Basis/grid code checks representability; k generation uses checked multiplication and fallible reservation. | Memory metadata is produced after work and excludes scratch/projectors/allocator overhead and all-k aggregation. `validate` itself enumerates basis vectors. No preallocation resource budget, reliable dry-run peak estimate or user-facing timeout setting. |
| Reproducibility | `run.rs`: input/UPF hashes, seed/options, tool version, timestamp, threading. `build.rs` captures HEAD/rustc. `tests/run.rs` checks deterministic output and density; stored references carry provenance. | No source-dirty fingerprint, full build feature/dependency fingerprint or complete restart identity. HEAD alone does not distinguish modified builds. Output paths can disclose local paths when results are shared. |
| Bands / DOS | Current `ScfResult`/JSON contains SCF k eigenvalues and occupations; `main.rs` lacks `bands`/`dos`. P8A specifies NSCF paths, CSV and DOS. | Current eigenvalues are not a band-path workflow or a converged DOS. Deliver P8A and mark energy reference, band accuracy and sampling explicitly. |
| Volumetric | `run.rs` emits little-endian f64 C-order `density.bin` and metadata, in electrons/Bohr^3. `scripts/compare_rhog.py` validates density against stored Fourier summaries. P8B specifies cube/XSF, ELF, potentials and partial charge. | No cube/XSF writer or end-user viewer recipe implemented in this checkout. Raw density requires custom reading. P8B is planned, not completed. |
| Charge analysis | Core has valence density/orbitals; UPF reader retains optional atomic wavefunctions through `raw()` (`docs/upf-support.md`). No analysis command in `main.rs` or corresponding analysis module in source inventory. | No Bader or Lowdin workflow. Bader needs a grid-convergence and valence-density interpretation contract; Lowdin needs an explicit localized projector basis, overlaps and normalization. Neither is a synonym for integrating the total cell charge. |
| Work function / dielectric | Potentials/Poisson exist in `crates/dft-physics/src`; `main.rs` and `ScfResult` expose no work-function or dielectric calculation. Coupling fields remain design in ADR 0017. | Work function needs a justified vacuum reference, potential definition, slab convergence and energy reference; total local KS potential is not automatically electrostatic voltage. Dielectric response needs a response/field formalism, and ionic contribution needs derivatives. These are not CSV-only tasks. |
| Bindings / interoperability | Rust crate API, CLI and Python subprocess tooling exist. No Python extension manifest or ASE calculator in tracked inventory. ADR 0017 keeps serialization/frameworks outside core; it is marked proposed/design-only. | Prefer an external file/process adapter first; advertise energy-only capabilities until forces/stress actually exist. ABI, in-process binding and array ownership can wait. |
| Documentation / tutorials | `docs/INDEX.md`, `getting-started.md`, `cli.md`, `validation.md`, examples, ADRs and generated per-crate guides exist. | Drift: CLI help says Gamma-only LDA; getting-started lists PBE out of scope despite `config.rs::Functional::Pbe`; early CLI documentation says one execution thread despite current worker controls. No completed noncoder plot/view walkthrough from P8A/P8B yet. |
| Benchmarks / tests / CI | `scripts/check.sh` builds/tests/formats/clippies; `bench_scf.py` includes Si8 and compares baseline energy/density; Python reference suite and ignored release tests exist. `tests/reference` has real Si LDA/k/PBE evidence and performance records. | No tracked hosted CI workflow, fuzz harness or coverage configuration found. Default gate omits ignored tests and Python suite. A named benchmark case is not a scheduled regression gate; Si8 can fail while default gates pass. |
| Distribution / license / citation | Root MIT `LICENSE`; workspace `publish=false`; pinned toolchain and `Cargo.lock`; `.cargo/config.toml` points to `../rust-dft-vendor`; source-build guide. | No tracked release binaries recipe, container recipe or `CITATION.cff`. Remote release assets were not checked. A self-contained source/release workflow must resolve the sibling vendor dependency without shipping restricted UPFs. |
| Parser security | UPF supplemental scanner validates tag nesting, duplicate attributes, unsupported declarations and dimensions; negative tests in `crates/dft-upf/tests/parser.rs`. Strict TOML and expression tests exist. | Whole-file reads and multiple parsed/normalized representations have no application-level byte/depth/resource budget. No fuzz campaign was found. This is an untested robustness boundary, **not evidence of arbitrary code execution or an XXE exploit**. |

## Ranked gap backlog and estimates

Each row is a bounded proposal. H/M/L are high/medium/low user value or risk. Dependencies are exit conditions, not permission to broaden this survey's ownership.

| ID | Gap and priority | User value for targets | Agent-days | Risk | Dependencies |
|---|---|---|---:|---|---|
| T1 | P0: Davidson robustness and systematic SCF stress regressions | H: defects/alloys require supercells; accidental success is insufficient | 3–6 | H, eigen/certificate semantics | P7PC diagnosis and existing performance/equivalence gates |
| T2 | P0: full preflight, element identity, coherent error classes | H: avoids wrong-material and wasted-job outcomes for all users | 2–4 | M, aliases/backward compatibility | Explicit species-versus-element contract; reuse engine checks |
| T3 | P1: plot/view outputs and noncoder walkthroughs | H: directly satisfies requested easy plotting | 6–10 total P8A/P8B envelope | H, NSCF/ELF/reference conventions | P8A/P8B, strict-band accuracy; stable metadata T5 |
| T4 | P1: PP registry plus expanded UPF2/ONCV compatibility | H: prerequisite for multi-element oxides/cathodes/alloys | 5–10 initial compatible subset | H, scientific validity/license | Exact licensed files, generator audit; NLCC and other physics remain P8D |
| T5 | P1: result schemas, run identity and transactional bundles | H: AI workflows must distinguish stale/partial/failed runs | 2–4 | M, compatibility and crash handling | ADR 0010/0013; existing readers and golden outputs |
| T6 | P1: external structure conversion and ASE energy adapter | H: noncoders, defect/supercell imports | 3–5 | M, cells/species/disorder | T2/T5; dependency/license review; ADR 0017 boundary |
| T7 | P1: k-spacing input and cutoff/k convergence commands | H: comparable target-material energies need converged sampling | 3–5 | M, skew cells/definitions | Existing scan tool; T5; smearing required for future metals (P8D) |
| T8 | P1: density checkpoint/restart | H: expensive supercell and multistart workflows | 4–7 | H, state compatibility | ESLE phase B/ADR 0016; T5; unchanged fresh-run gate |
| T9 | P1: tiered CI, adversarial inputs and broader numerical matrix | H: prevents repeated regressions during AI development | 3–6 initial setup | M, runtime/corpus availability | T1 cases; redistributable synthetic corpus; remote real PP access |
| T10 | P2: installable release, tutorial and citation package | H for noncoders; less immediate than correctness | 2–4 | M, platform/vendor reproducibility | T3/T5; platform builds and approved dependency provenance |
| G11 | P2: live SCF progress and structured failure history | H for long supercells and batch jobs | 2–3 | M, overhead/cancellation | Engine observer outside kernels; T5/T8; performance rules |
| G12 | P2: preallocation memory estimate/resource budgets | H for defects and high cutoffs | 2–4 | M, conservative bounds | Basis enumeration bound, all-k bands/scratch accounting |
| G13 | P2: space-group detection and reciprocal reduction | M/H for ordered crystals; lower for disordered defects | 4–7 | H, rotations, tolerances and weights | Library review, full-grid oracle, spin-aware deferral |
| G14 | P2: Bader interchange and charge-analysis guide | M/H for cathodes and bonding trends | 2–4 adapter/validation | H, valence versus total charge interpretation | P8B, external tool/license review, grid convergence |
| G15 | P2: Lowdin populations / projected DOS | H for transition metals, but physics interpretation matters | 5–9 | H, projector dependence | Validated atomic-wavefunction adapter, overlap machinery, P8A |
| G16 | P2: work-function workflow | M/H for surfaces/electrodes | 4–8 after potential export | H, vacuum/electrostatics/reference | P8B, slab and boundary treatment, occupation/reference policy |
| G17 | Deferred: electronic/ionic dielectric response | M/H for insulating oxides; not the next usability task | 10–20+ preliminary | H, new response physics | P8D/field design, polarization/response validation; derivatives for ions |
| G18 | P2: UPF v1 converter/import | M, useful for legacy collections; below a modern validated set | 3–6 | H, old metadata/units | Published format fixtures and source-independent parser; T4 |
| G19 | Deferred: in-process Python/other bindings | M once process adapter proven | 4–8 | H, ABI/ownership/threading | Stable engine state API, T6, ADR 0017 |
| G20 | Deferred: MPI/GPU scaling and portable performance controls | H eventually for large targets | 15–30+ architecture-dependent | H | Measured CPU/memory limits, backend contracts, separate approval/ticket |

T4's estimate does **not** include new NLCC/PAW/USPP physics or validating a whole periodic table. T8 first delivers a density warm start, then orbital checkpoints/exact continuation as a separately estimated extension (roughly 4–7 more agent-days). G17/G20 are order-of-magnitude placeholders, not commitments. Retain CPU thread controls already present.

### Quick wins: independently reviewable work under one agent-day each

| Small deliverable | Value | Agent-days | Risk | Dependency / proof |
|---|---|---:|---|---|
| Q1 Correct CLI help and contradictory introductory docs | H, all users | 0.25–0.5 | L | Current command/config implementation; help snapshot and docs review |
| Q2 Export existing SCF history/energy/eigenvalues as CSV from JSON | H, immediate spreadsheet plotting | 0.5–0.75 | L/M | Unit headers, convergence flag, k index and reference; test legacy and multi-k files; not new bands/DOS |
| Q3 Add a full-preflight recipe naming `validate` limitations and remedies for exits 1–5 | H | 0.25–0.5 | L | Reproduce experiments below; distinguish preparation from convergence |
| Q4 Include UPF header element in species output and console | H | 0.5–0.75 | M | Additive field test; warning/identity policy remains T2, don't forbid useful aliases blindly |
| Q5 Add synthetic odd-count/zero-budget/low-cutoff cases to an explicit smoke matrix | H, AI regression signal | 0.5–0.75 | L | Existing tests; actionable failure classes; no real UPF redistributed |
| Q6 Add citation metadata template and a source-build checklist | M | 0.25–0.5 | L | Maintainer-supplied authors/release information; never invent DOI/authors |

These are useful standalone steps; no claim that all of T2, T3, T9 or T10 fits in one day. Suggested first batch is Q1–Q3 while P7PC diagnosis proceeds; then T2/T5 and existing P8A/P8B. Wider target-material enablement depends on T4 and P8D together.

## Concrete robustness experiments

### Environment and reproduction

Built this checkout with the pinned toolchain, default pure Rust backend and one requested thread. The approved UPF was only read at `/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF`; SHA-256 was independently checked as `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`. Its unresolved provenance/license remains a limitation. Experiment timings are not benchmarks (builds/tests shared the machine).

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export RUSTUP_HOME=/private/tmp/p0c-rustup
export CARGO_NET_OFFLINE=true
CARGO_TARGET_DIR=/tmp/p8e-target cargo build --locked --offline --release -p dft-cli
mkdir -p /tmp/p8e-survey
# Save the complete Python listing below as /tmp/p8e-survey/experiments.py.
python3 -I /tmp/p8e-survey/experiments.py > /tmp/p8e-survey/experiment-log.jsonl
```

Both build and Python driver exited 0; driver exit 0 means experiments were recorded, not that child runs all succeeded. Each actual child argv and status is recorded in `records.json`; the table reports those child statuses. Reproduction should use an empty temporary directory or remove only prior experiment outputs first: the runner reports existing result files if reused after failure.

### Observed outcomes and classification

| Input / command variant | Actual exit | Observation | Classification |
|---|---:|---|---|
| Baseline two-atom Si, 30 Ry, legacy Gamma, `run` | 0 | 8 electrons, 12 SCF steps, total **-7.29553046821089 Ha** | Successful control, not a new QE comparison |
| Missing UPF, `validate` then `run` | 0 / 5 | Validator explicitly says file existence not checked; run names absent path | Acceptable documented limitation; full-preflight usability gap |
| Atom label Fe without matching map entry | 1 | `species "Fe" has no pseudopotential entry` | Acceptable rejection |
| Rename both atom and map label to Fe, retain Si UPF | 0 | Same 8 electrons, 12 steps, exactly same energy as control; no stderr warning | **Missing validation/identity visibility**, not evidence of an Fe calculation. Labels are aliases by current API; require separate chemical identity and expose header element |
| Cutoff 0.000001 Ry, `validate` then `run` | 0 / 4 | `Numerical("k index 0: InvalidInput(\"bands/subspace incompatible with k basis\")")` | **Bug in error classification** plus incomplete preflight; incompatible input should identify available basis/bands and exit 1, not imply solver divergence |
| kgrid [0,1,1] | 1 | Positive grid/weight/shift error | Acceptable rejection |
| kgrid [1,1,1] | 0 | Same energy and steps as legacy Gamma | Acceptable sampling input; convergence with respect to k mesh is not established by this success |
| Near-singular lattice ([1,0,0], [1,1e-13,0], [0,0,1]) Angstrom | 1 | Scaled determinant <=1e-12 | Acceptable geometry rejection |
| Skew lattice ([3,0,0], [29,0.7,0], [41,2.3,0.2]) Angstrom, `grid` only | 0 | Basis/grid inspection completes (details below) | Successful geometry preparation only; SCF accuracy on skew cells untested |
| Cube side 1e12 Angstrom | 1 | `basis bounds or FFT storage size cannot be represented` | Acceptable representability rejection; does not establish graceful handling of realistic memory exhaustion |
| `scf.max_iterations=0`, `validate` then `run` | 0 / 1 | Run says iteration count/tolerances must be positive and finite | **Missing full validation** in validator; deliberate shallow API today, confusing for users |
| `max_iterations=1` | 3 | Output saved with `converged=false`, 1 step, **-7.172342783970571 Ha** | Acceptable budget exhaustion, not a usable converged energy |
| Legacy Gamma + `extra_bands=0` | 0 | **-7.2955304682108695 Ha**, 12 steps | Acceptable compatibility limitation; no sampled-gap certificate |
| [1,1,1] grid + `extra_bands=0` | 1 | `gap certification requires an empty band at every k` | Acceptable documented new-API restriction (ADR 0013); tutorial should explain difference |
| Conventional 8-atom Si, 30 Ry, defaults | 4 | After 400 Davidson iterations, residual **3.953844866590969e-5 Ha**, no new result | **Numerical robustness bug**, reproduced independently of fast-gate success; P7PC owns diagnosis/fix |
| Odd electron synthetic engine test | 0 (test process) | `invalid_electrons_bands_and_options` passes; asserts `OddElectronCount(3.0)` for one synthetic Z=3 atom | Acceptable explicit unsupported occupation; not a real odd-electron material or CLI test |
| Synthetic DTD / truncated UPF through `inspect-upf` | 1 / 1 | `unsupported XML declaration` / `unterminated tag` | Acceptable malformed-input rejection; not a full security assessment |

The Si8 residual differs from the historical residual in P7PC. This is the current measured result for the exact atom order below; the ticket's numbers are not relabeled as reproduced measurements. A failure is reproduced, but its degeneracy mechanism is not proven by this survey. Do not weaken tolerances as a workaround.

Skew-cell preparation returned volume **2.834300487732 Bohr^3**, 7 wavefunctions, 65 density-sphere vectors and FFT dimensions **20 × 180 × 250 = 900,000 samples**. This shows why plane-wave count alone is a poor memory predictor for skew cells; it is not a measured peak RSS.

For the same Si8 input, changing only `extra_bands` to 1 gave exit 0, **-31.32180360103159 Ha**, 14 SCF steps; setting it to 8 gave exit 0, **-31.321803601031554 Ha**, 11 steps. Their difference is about **3.6e-14 Ha**. This supports a band-count-sensitive failure, not a proof of the specific stalled multiplet.

### Complete input generator and child commands

```python
from pathlib import Path
import subprocess, json, hashlib, time
root=Path('/tmp/p8e-survey')
binary='/tmp/p8e-target/release/dft-cli'
upf=Path('/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF')
assert hashlib.sha256(upf.read_bytes()).hexdigest()=='e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217'
base=Path('examples/si_lda.toml').read_text().replace('/absolute/path/to/Si.UPF',str(upf))
lat='[[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]'
cases={
'baseline':base,
'missing':base.replace(str(upf),'/tmp/p8e-survey/absent.UPF'),
'unknown_species':base.replace('species = "Si"','species = "Fe"'),
'wrong_element':base.replace('species = "Si"','species = "Fe"').replace('Si =','Fe ='),
'low_cutoff':base.replace('30.0','0.000001'),
'zero_kgrid':base.replace('kpoints = [[0.0, 0.0, 0.0]]','kgrid = [0, 1, 1]'),
'one_kgrid':base.replace('kpoints = [[0.0, 0.0, 0.0]]','kgrid = [1, 1, 1]'),
'singular':base.replace(lat,'[[1.0,0.0,0.0],[1.0,1e-13,0.0],[0.0,0.0,1.0]]'),
'skew':base.replace(lat,'[[3.0,0.0,0.0],[29.0,0.7,0.0],[41.0,2.3,0.2]]'),
'huge':base.replace(lat,'[[1e12,0.0,0.0],[0.0,1e12,0.0],[0.0,0.0,1e12]]'),
'scf_zero':base+'\n[scf]\nmax_iterations = 0\n',
'budget':base+'\n[scf]\nmax_iterations = 1\n',
'gamma_no_empty':base+'\n[scf]\nextra_bands = 0\n',
'grid_no_empty':base.replace('kpoints = [[0.0, 0.0, 0.0]]','kgrid = [1, 1, 1]')+'\n[scf]\nextra_bands = 0\n',
}
si8=base.replace(lat,'[[5.43,0.0,0.0],[0.0,5.43,0.0],[0.0,0.0,5.43]]')
extra=[[0,.5,.5],[.5,0,.5],[.5,.5,0],[.25,.75,.75],[.75,.25,.75],[.75,.75,.25]]
si8=si8.replace('[pseudopotentials]', ''.join('[[system.atoms]]\nspecies = "Si"\nfractional = '+str(p)+'\n\n' for p in extra)+'[pseudopotentials]')
cases['si8']=si8
records=[]
for name,content in cases.items():
    path=root/(name+'.toml'); path.write_text(content)
    commands=['validate','run'] if name in ['missing','low_cutoff','scf_zero'] else ['grid'] if name=='skew' else ['run']
    for command in commands:
        args=[binary,command,str(path)]
        out=root/(name+'-'+command)
        if command=='run': args+=['--output-dir',str(out)]
        start=time.monotonic()
        try:
            p=subprocess.run(args,capture_output=True,text=True,timeout=45)
            (root/(name+'-'+command+'.stdout')).write_text(p.stdout)
            (root/(name+'-'+command+'.stderr')).write_text(p.stderr)
            record=dict(name=name,command=args,exit_code=p.returncode,seconds=round(time.monotonic()-start,3),stderr=p.stderr.strip())
            if (out/'results.json').exists():
                r=json.loads((out/'results.json').read_text())
                record.update(energy_ha=r['energies']['total']['ha'],electrons=r['electrons'],iterations=r['convergence']['iterations'],converged=r['convergence']['converged'])
        except subprocess.TimeoutExpired:
            record=dict(name=name,command=args,timeout_seconds=45)
        records.append(record); print(json.dumps(record),flush=True)
(root/'records.json').write_text(json.dumps(records,indent=2)+'\n')
```

Additional experiments were run with these exact commands (Python drivers exited 0):

```sh
python3 -I - <<'PY'
from pathlib import Path
import subprocess,json
root=Path('/tmp/p8e-survey'); records=[]
for extra in [1,8]:
    name=f'si8-extra{extra}'
    path=root/(name+'.toml'); path.write_text((root/'si8.toml').read_text()+f'\n[scf]\nextra_bands = {extra}\n')
    args=['/tmp/p8e-target/release/dft-cli','run',str(path),'--output-dir',str(root/name)]
    p=subprocess.run(args,capture_output=True,text=True,timeout=45)
    (root/(name+'.stdout')).write_text(p.stdout)
    r=dict(command=args,exit_code=p.returncode,stderr=p.stderr.strip())
    if (root/name/'results.json').exists():
        j=json.loads((root/name/'results.json').read_text()); r.update(energy_ha=j['energies']['total']['ha'],iterations=j['convergence']['iterations'])
    records.append(r); print(json.dumps(r),flush=True)
(root/'si8-variants.json').write_text(json.dumps(records,indent=2))
PY

python3 -I - <<'PY'
from pathlib import Path
import subprocess,json
root=Path('/tmp/p8e-survey')
for name,text in [('dtd','<!DOCTYPE UPF [<!ENTITY x "x">]><UPF version="2.0.1"></UPF>'),('truncated','<UPF version="2.0.1"><PP_HEADER')]:
    path=root/(name+'.xml'); path.write_text(text)
    args=['/tmp/p8e-target/release/dft-cli','inspect-upf',str(path)]
    p=subprocess.run(args,capture_output=True,text=True)
    r=dict(command=args,exit_code=p.returncode,stderr=p.stderr.strip())
    (root/(name+'.json')).write_text(json.dumps(r,indent=2)); print(json.dumps(r))
PY

CARGO_TARGET_DIR=/tmp/p8e-check-target cargo test --locked --offline -p dft-engine --test scf invalid_electrons_bands_and_options -- --exact --nocapture
```

The odd-count command exited 0: 1 passed, 0 failed, 7 filtered out. The test explicitly constructs a synthetic Z=3 potential; the approved Si potential cannot produce an odd count by changing the number of neutral Si atoms. No real potential was edited to manufacture this test.

## Proposed top-ten ticket scopes and gates

These are proposed future ownership sets, not permission to edit them in this task. Frozen schema/manifests require an explicitly authorized ticket; no dependency is fetched by this survey. Each implementation ticket also runs `scripts/check.sh`, applicable Python tests, and records evidence. Verified new results must flow through the README figure generator under the repository rule. Coordinate shared CLI/config files serially or with explicit ownership; P8A/P8B already claim overlapping areas.

| Ticket / title | Proposed owned files | Acceptance gates beyond the common gate |
|---|---|---|
| T1 — Finish P7PC: robust buffer bands and supercell regression | Existing P7PC ownership: `crates/dft-solvers/**`, `crates/dft-engine/**`, focused `crates/dft-cli/**`, `scripts/bench_scf.py`, specified ADR appendices and `docs/performance/optimization-log.md` | Diagnose per-band residuals; reproduce default failure before fix; Si8 seeds 1/7 and threads 1/4 converge under unchanged strict-band/certificate tolerances; retain historical caller behavior; real benchmark/QE comparison gates exactly as P7PC; long cases ignored in debug with release command |
| T2 — Full preflight and species/element identity | `crates/dft-cli/src/{main,config,run}.rs`, new `preflight.rs`, `crates/dft-cli/tests/preflight.rs`, `docs/cli.md`; narrowly scoped engine error propagation if required | Full check opens/hashes UPFs, validates SCF and every k basis before solving; report alias plus header element; explicit expected-element mismatch fails, alias remains allowed; missing file exit 5, invalid settings/band capacity exit 1; all experiments have stable actionable diagnostics |
| T3 — Plot-ready bands/DOS/density delivery | Preserve exact P8A/P8B ownership rather than creating duplicate modules; Q2 additionally proposes `scripts/export_results_csv.py`, `scripts/tests/test_export_results_csv.py`, `docs/howto-results.md` | Independent CSV parse; units and energy reference explicit; nonconverged status survives export; k-point/band labels retained; P8A integration/normalization gates and P8B skew/left-handed cell, density-integral and round-trip gates; tutorial produces a first plot/view with exact commands |
| T4 — Validated PP metadata catalog and UPF2 compatibility | New `docs/pseudopotentials.md`, `tests/reference/pseudopotential-manifest.json`, `scripts/check_pseudopotentials.py` and tests; `crates/dft-upf/src/parser.rs`, parser tests, `docs/upf-support.md`; dependency change only through separately approved manifest ownership | No UPF committed; manifest has source URL/version, SHA-256, element, valence, XC, relativistic/core flags, generator, license evidence, supported status and measured cutoff-study references; unsupported features fail closed; synthetic >3-projector coverage plus licensed external real-file tests; inspect support never implies validated accuracy |
| T5 — Stable result contract and atomic completion manifest | New `docs/schemas/results.schema.json`, `density.schema.json`, `run-manifest.schema.json` (explicit frozen-path authorization); `crates/dft-cli/src/run.rs`, contract tests, `docs/cli.md`, affected Python readers | Old legacy Gamma/k/PBE outputs validate; additive evolution policy tested; consumers reject unsupported versions and failed/incomplete identity; inject failures between writes and detect mixed/stale bundles; hashes cover payloads; publish manifest last; distinguish source/build fingerprints from input hashes |
| T6 — External structure/ASE bridge, energy only | New external `adapters/ase/` or `scripts/import_structure.py` and adapter tests, `docs/howto-import.md`; dependency record if selected | Round-trip CIF/POSCAR/extxyz fixtures preserve cell/positions/species; XYZ without explicit periodic cell rejected; fractional/cartesian unit and handedness checks; disorder/partial occupancy rejected explicitly unless a defined conversion exists; subprocess failures propagate; forces/stress requests report unsupported, never zeros; no framework dependency in numerical crates |
| T7 — Reproducible cutoff/k scans and spacing input | `scripts/convergence_scan.py`, `scripts/tests/test_scan.py`, `crates/dft-cli/src/{main,config}.rs`, geometry k-grid helper/tests if shared, `docs/howto-convergence.md` | Define spacing in inverse Bohr or inverse Angstrom, rounding and shift; test anisotropic/skew cells; resolved mesh stored; independent cutoff and k series with per-atom energy and sampled gap; failed points retained, not silently skipped; scan CSV supports multiline inputs through a reliable editing contract |
| T8 — ESLE B density restart and checkpoint envelope | `crates/dft-engine/src/scf.rs`, new state/compatibility modules in existing engine crate, engine restart tests; new CLI persistence module, wiring/tests; ESLE architecture/CLI docs | Fresh path unchanged; interrupted write cannot load; wrong hash/cell/grid/k/XC/PP/electron identity rejected; same-model restart conserves charge and reaches certified state within measured noise; unknown schema fails closed; no claim of exact continuation without mixer/orbital/schedule state |
| T9 — Tiered CI and adversarial numerical/parser regressions | New `.github/workflows/ci.yml` (or approved CI equivalent), `scripts/check_release.sh`, parser/input mutation harness and corpus under tests, existing numerical test files, `docs/testing.md`; fuzz dependencies only with approval | Fast offline gate + Python suite on PR; release ignored tests scheduled; external-Si tests opt-in with hash check and no redistributed payload; Si8 matrix is a real gate; bounded byte/depth/time parser corpus asserts no panic; malformed k/cell/resource cases; selected platforms produce real results before portability claims; compiler/lint failures remain visible |
| T10 — Reproducible install and first-result tutorial | New release build/packaging scripts, optional `packaging/Dockerfile`, CI release workflow, `CITATION.cff`, `docs/install.md`, examples and tutorials; `.cargo/config.toml` change only with explicit ownership | Clean-machine install without developer sibling paths; checksum/version and dependency license notices; supported-platform smoke run with synthetic PP and separate user-provided real-PP guide; no restricted PP in package/container; documentation commands exercised; authors/citation metadata supplied by maintainers |

## Design guidance and ecosystem evidence

The ecosystem references below support candidate directions; they are not evidence that an adapter already works here. No third-party source code was copied, and no dependency adoption is approved by this survey.

- **Structure conversion outside the core:** ASE documents CIF, VASP, XYZ/extxyz and other file I/O; use an optional external adapter to convert to the engine's small input contract. Explicitly handle missing cells, nonperiodic structures and ambiguous occupancies. Review distribution/license terms before shipping an adapter; using an external process does not settle all packaging obligations. [ASE file I/O documentation](https://docs.ase-lib.org/ase/io/io.html).
- **Symmetry as a separately validated optimization:** spglib documents crystal symmetry and irreducible reciprocal meshes. A candidate adapter should retain the user's cell, store mapping/weights/tolerances, and compare with full-grid energies/densities before enabling reduction by default. Do not silently standardize a defective/alloy structure. [spglib C API](https://spglib.readthedocs.io/en/stable/api.html).
- **Pseudopotentials require quality evidence:** the PseudoDojo authors describe ONCV tables and cutoff-dependent testing. That motivates storing validated cutoff-study provenance, not trusting a nonzero header recommendation or assuming every published table works with this reader. This survey did not verify licenses of particular PseudoDojo/SG15 files and does not authorize download/redistribution. [PseudoDojo paper](https://arxiv.org/abs/1710.10138). The attempted website fetch timed out; use the paper only for this limited claim.

For T7, a reasonable proposed convention is `n_i = max(1, ceil(|b_i| / spacing))` with reciprocal vectors including 2π; document that edge spacing is not an isotropic nearest-neighbor guarantee in skew cells. Store the resolved mesh and allow overrides. Cutoff suggestions should name observable, tolerance, units and PP hash, with no suggestion when evidence is absent. For slabs, anisotropic/vacuum sampling must be user-controllable.

For G14–G17, require independent observable-specific validation before advertising scientific conclusions. Valence pseudo-density charge partitions are not all-electron oxidation states; projection populations depend on the chosen localized basis. A potential export needs its zero-mode and physical role, not only the label “potential.” A dielectric number must identify electronic versus ionic contribution and tensor/frame. These are proposed acceptance constraints, not claims that current files encode those observables.

## Coverage findings and what to test next

The existing tests are substantial: CLI determinism/exit/path safety (`crates/dft-cli/tests/run.rs`), parser negative and radial analytic cases (`crates/dft-upf/tests/{parser,radial}.rs`), solver true-residual checks (`crates/dft-solvers/src/product_tests.rs`), SCF invariances and fixed-point tests (`crates/dft-engine/tests/scf.rs`), k-point/threading/PBE suites and stored independent-reference scripts. The gap is breadth and execution policy, not an absence of tests.

The fast shell gate does not invoke `scripts/tests`, and Cargo does not run `#[ignore]` tests by default (`scripts/check.sh`; ignored annotations in engine/CLI tests). `scripts/bench_scf.py` includes `si8`, but its default cases are `gamma,k222`; it is not invoked by the fast gate. Add nightly/explicit release gates rather than forcing expensive physical runs into every debug test.

Prioritize a reproducible Cartesian test matrix: occupied/buffer count near degeneracies; seed and thread layouts; atom-order/permutation equivalence; primitive versus conventional supercells; cutoff where basis capacity is insufficient; left-handed/anisotropic/skew cells; zero/tiny/overflowing k grids; metallic/odd occupation rejection; missing/extra/mismatched species; malformed and oversized files. Existing tests cover parts of this matrix, so extend them rather than claim each category has zero coverage. For multi-element/transition-metal tests, obtain independently validated compatible potentials first. No line/branch coverage percentage was measured.

Static resource concerns need bounded tests, not destructive OOM experiments: preflight currently enumerates vectors before a practical allocation budget; UPF input is read whole and normalized; time-reversal reduction searches accumulated representatives. These paths merit size/depth/work limits and profiling. No hostile large-file run or memory exhaustion of the engine was attempted. Output transaction fault injection, dirty-build fingerprint tests and schema-consumer compatibility are especially valuable for AI-authored changes that otherwise pass numerical tests.

## Verification record and evidence limits

The following are real process results from this survey. The first full check failed for environment storage, so it is retained rather than rewritten as a pass. A subsequent retry completed all four gates without changing code.

| Command (repository root) | Exit | Result |
|---|---:|---|
| `env PATH=/opt/homebrew/opt/rustup/bin:$PATH RUSTUP_HOME=/private/tmp/p0c-rustup CARGO_TARGET_DIR=/tmp/p8e-target CARGO_NET_OFFLINE=true cargo build --locked --offline --release -p dft-cli` | 0 | Release CLI built from surveyed checkout |
| `python3 -I /tmp/p8e-survey/experiments.py > /tmp/p8e-survey/experiment-log.jsonl` | 0 | 18 child invocations recorded; no timeout; individual statuses above |
| `env PATH=/opt/homebrew/opt/rustup/bin:$PATH RUSTUP_HOME=/private/tmp/p0c-rustup CARGO_TARGET_DIR=/tmp/p8e-check-target CARGO_NET_OFFLINE=true cargo test --locked --offline -p dft-engine --test scf invalid_electrons_bands_and_options -- --exact --nocapture > /tmp/p8e-odd.log 2>&1` | 0 | 1 passed, including synthetic odd-count rejection |
| `env PYTHONDONTWRITEBYTECODE=1 python3 -I -B -m unittest discover -s scripts/tests > /tmp/p8e-python.log 2>&1` | 0 | 113 tests, 1 skipped; stored-reference checks included |
| `env PATH=/opt/homebrew/opt/rustup/bin:$PATH RUSTUP_HOME=/private/tmp/p0c-rustup CARGO_TARGET_DIR=/tmp/p8e-check-target CARGO_NET_OFFLINE=true scripts/check.sh > /tmp/p8e-check.log 2>&1` | 1 | Build/test/fmt PASS; clippy interrupted by `No space left on device` writing temporary dependency metadata |
| Same full-check command, redirecting to `/tmp/p8e-check-retry.log` | 0 | Build/test/fmt/clippy all PASS |

The Python suite rechecked existing QE artifacts, not a new QE execution. It printed stored k222/k222s maximum eigenvalue differences of **1.9874921852982487e-6 / 1.77285386153842e-6 Ha**, and total-energy differences of **1.2814394616533775e-6 / 1.281265464836423e-6 Ha/atom**. Those numbers belong to the stored reference fixtures processed by the tests; they are not fresh comparison results for the baseline or Si8 experiments above. No claim of general material accuracy follows.

Evidence limits:

- Current numerical experiments cover one approved real Si potential, two Si cell sizes, a geometry-only skew example, synthetic odd-electron rejection and two tiny malformed XML inputs. No real oxide, cathode, alloy, magnetic, charged-defect or surface was solved. No real UPF v1/ONCV corpus was obtained or checked.
- No new remote QE run, external viewer run, ASE/spglib integration, Bader/Lowdin calculation, work function or dielectric response was performed. External sources only informed proposals. Candidate dependencies still require license, maintenance and Apple Silicon review under ADR 0003.
- Si8 buffer sensitivity was measured once per setting, at default seed 2026/thread request 1. Seeds 1/7, thread sweeps, per-band residual history and a degeneracy diagnosis remain P7PC work. No new performance optimization or timing baseline is claimed.
- Full check success establishes the default fast build/test/lint gates only. Ignored release numerical suites were not all rerun. Python suite has one skip. No coverage percentage, fuzz campaign, parser resource-bound guarantee, sanitizers, OOM resilience or concurrency fault injection was established.
- The huge-cell rejection tests integer representability, not realistic large-supercell memory behavior. The skew case exercises preparation only; 900,000 FFT samples with 7 wavefunctions is a resource-planning warning, not an SCF correctness result.
- The transient build-storage failure was observed in temporary compiler output, not caused by an engine run. Disk availability later permitted a clean retry. This does not diagnose system storage policy or claim an engine bug.
- No public release service or remote CI settings were queried. “No release/CI recipe found” is a tracked-checkout finding. Packaging platforms beyond this local build were not validated.
- Effort/priorities are proposals for the stated user/material goals. Physics readiness may dominate workflow readiness; T4 cannot make unsupported magnetic/metallic models valid. Scheme aliases are legitimate, so the wrong-element experiment motivates explicit identity rather than banning all non-element labels.
- Only `docs/research/feature-gap-survey.md` is owned. Session registry, dashboard/artifact JSON, README, figures, tickets and code are unchanged, as explicitly required by this ticket. No implementation was made, so no new README validation result was generated. Changes are left **uncommitted** for review.

Recommendation: finish the reproduced P7PC failure and clarify preflight/identity first; deliver existing P8A/P8B alongside stable result contracts; then expand validated pseudopotentials, external structure import, convergence scans and density restart. Maintain an explicit separation between tested capabilities, accepted limitations and future physics.
