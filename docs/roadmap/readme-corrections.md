# README factual corrections applied (D6)

Prepared for R12-DOCS against README.md at `27d91ece` (2026-11-25).
**Application: applied (D6 approved 2026-11-26), option A, by R12B-README.**
The proposals below are retained as the approval record. Line numbers refer to the original file;
“current text” is copied exactly (line 84 uses a marked excerpt), so the reviewer
can check each replacement without repeating forbidden mixed-model speed ratios.

The generator owns **two** blocks: `status:begin/end` (from
`docs/feature-status.json`) and `validation:begin/end` (stored evidence).
The first table covers prose outside both blocks. The second includes the
status statements explicitly requested by R12 (lines 94/101/104/112) and other
stale status entries. R12B-README updated their source JSON and regenerated;
generated text was not hand-edited. Freshness does not establish semantic currency.
The full README was reviewed; unchanged research ambitions, license, build
statements and explicit material-validation limits need no correction.

## User-maintained prose

### README.md line 35

Status: applied (D6 approved 2026-11-26).

Current text:

```text
Interoperability is an architectural requirement. BAANDI aims to exchange structures, densities, potentials, energies and electronic-state metadata (and forces and stresses once they exist) with external tools through physically consistent interfaces: QM/MM coupling, explicit and implicit solvents, machine-learned interatomic potentials, multiscale simulations. External methods connect through adapters and drivers outside the numerical core. BAANDI does not aim to replace specialized molecular dynamics, machine learning or continuum engines.
```

Proposed text: Replace “(and forces and stresses once they exist)” with “including the supported force and stress outputs”; keep external coupling as a future integration direction.

Reason: Forces/stress are implemented; this does not claim QM/MM/solvent integration. See docs/validation/forces.md and stress.md.

### README.md line 49

Status: applied (D6 approved 2026-11-26).

Current text:

```text
Rules that hold for every layer: the stepwise path must give bit-identical results to the conventional path; a change to a running calculation is validated, applied only between SCF iterations, rejected without side effects if unsafe, and recorded; convergence is declared only by the existing fixed-point certificate, never by the control layer; diagnostics separate measured evidence from hypothesis. The first end-to-end demonstration is a **self-healing SCF**: a failing calculation is diagnosed, repaired through a validated policy change, continued and verified, with the intervention recorded. Design and plan: [docs/research/agent-native-directive-review.md](docs/research/agent-native-directive-review.md) and the P11 tickets; none of layers B to E exists yet.
```

Proposed text: Rules: stepwise and batch execution preserve the same numerical trajectory. Controls apply at explicit safe map boundaries and are audited. The engine alone decides convergence; independent reference verification is on by default and can be disabled for supported outputs. Typed state, reference runtime, checkpoint v2, Python/ASE and host-embedded MCP exist. The historical self-healing demo demonstrates legacy budget recovery; the reference runtime rejects budget extension. MLIP/external solver integrations remain designs.

Reason: R8a–R8d implement B/C/E and part of D. Candidate/verification boundaries can be unsafe even between map returns. The frozen oracle certificate differs from production verification (ADR 0047).

### README.md line 56

Status: applied (D6 approved 2026-11-26).

Current text:

```text
- Periodic boundary conditions and Monkhorst–Pack k-point sampling (time-reversal symmetry only)
```

Proposed text: - Periodic boundary conditions and Monkhorst–Pack k-point sampling, time reversal and opt-in compatible spatial symmetry reduction

Reason: P19A/G6 implement spatial symmetry; it is no longer time reversal only.

### README.md line 60

Status: applied (D6 approved 2026-11-26).

Current text:

```text
- Block-Davidson Kohn–Sham eigensolver, Pulay density mixing, a fixed-point convergence certificate and an energy audit
```

Proposed text: - Reference Davidson eigensolvers, Broyden/plain density mixing, default-on fixed-point verification and an energy audit; a feature-gated strict oracle for validation

Reason: Production mixer/driver changed; Pulay/strict certification describe historical/oracle behavior (ADR 0047).

### README.md line 64

Status: applied (D6 approved 2026-11-26).

Current text:

```text
In design or development: richer SCF algorithms and an optional Intel MKL backend selected at compile time, band-structure/DOS and volumetric outputs in standard formats, electronic-state exploration, forces and stresses, external-field coupling, spin polarization, GPU backends, and the agent-native runtime above.
```

Proposed text: Implemented with documented limits: optional Intel MKL, bands/DOS and standard volumetric outputs, spin/FD, forces/stress and reference runtime. Research/design directions include richer algorithms, electronic-state exploration, external coupling and GPU backends.

Reason: The original list puts implemented capabilities in development. No new physical validation is implied.

### README.md line 66

Status: applied (D6 approved 2026-11-26).

Current text:

```text
These are development directions, not claims of currently available functionality; the status table below is the authoritative list.
```

Proposed text: Current supported workflows and limits are indexed in docs/INDEX.md and ADR 0047; dated evidence retains the implementation that produced it.

Reason: The existing status table contains stale implementation statements; it must not override the current contract.

### README.md line 70

Status: applied (D6 approved 2026-11-26).

Current text:

```text
**BAANDI is in an early experimental stage.** A verified silicon-class workflow exists (insulators, fixed occupations, LDA and PBE, Gamma and k-point meshes), but it is not a production DFT code. It should not be used for scientific conclusions without independent validation. Features such as metals/smearing, forces, stress, spin and structure relaxation are not available yet.
```

Proposed text: BAANDI remains experimental and needs independent validation for scientific conclusions. In addition to the silicon insulating workflow, it implements limited-validation collinear spin, Fermi–Dirac occupations, forces/stress and unpolarized fixed-occupation relaxation. Free-moment PBE remains unqualified and uses the strict-oracle fallback.

Reason: Remove only the false absence claim, retaining the experimental/material-validation limits.

## Generated status source corrections

These rows were applied through `docs/feature-status.json` and regenerated.

### README.md line 84

Status: applied (D6 approved 2026-11-26).

Current text (exact excerpt; the full line is in README.md):

```text
it has no real-wavefunction Gamma path yet
```

The same paragraph presents 2026-11-05 thread-versus-MPI ratios as current speed;
those ratios are deliberately not reproduced here.

Proposed text: Replace the entire speed paragraph with: “Current timing claims use the dated same-model tables and matched-condition protocol in docs/performance/comparison-rules.md and docs/performance/qe-vs-baandi-current.md. The 2026-11-05 relay describes earlier implementations and is not current reference-driver performance. Real-Gamma and MPI paths are implemented. No new speed measurement is supplied by this documentation update.”

Reason: The old paragraph describes no real-Gamma path and mixes BAANDI threads with QE MPI; such ratios are forbidden even as informational columns. Do not restate those ratios as current evidence.

### README.md line 90

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| Physics | Kohn-Sham DFT with plane waves | Periodic cells, norm-conserving pseudopotentials (UPF v2 subset), spin-unpolarized. |
```

Proposed text: Periodic NC plane-wave Kohn–Sham DFT; unpolarized and supported collinear sectors, with documented occupation/XC limits.

Reason: The overview says spin-unpolarized although collinear support is implemented.

### README.md line 94

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| SCF | Pulay and linear density mixing, block Davidson eigensolver | Default: one standard SCF path (threshold 5e-9 Ha) with a fixed-point verification after convergence (ADR 0043). The strict path with its certificate and energy audit is a reference tool (BAANDI_REFERENCE_STRICT=on). Verified standard results are not strict-certified. |
```

Proposed text: Broyden/plain density mixing and reference Davidson; beta 0.7/history 8, threshold 5e-9 Ha and verification on by default. StrictOracle is an explicit validation route and the temporary G4D fallback. Verified standard is not strict-certified.

Reason: R2/R10 routing and mixer defaults supersede Pulay/linear production controls. `pulay` is a Broyden alias, not the old Pulay algorithm.

### README.md line 97

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| Speed | Threads over k points and over bands, bit-identical results for any thread count | `[scf] threads = N`; default is 1. |
```

Proposed text: Thread budget via [scf] threads = N (default 1); same-build/backend determinism is checked separately from MPI rank-count accuracy. Reference k/spin/FD work is distributed under MPI.

Reason: Scope the old unconditional “any thread count” statement to the documented determinism contract; rank-count agreement is not bit identity.

### README.md line 101

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| Robustness | Convergence fix for larger cells (8-atom silicon cell stalled at the default settings) | Buffer bands: the highest requested bands are loosely converged; strict bands keep the full tolerance. The 8-atom silicon cell now converges. |
```

Proposed text: Reference convergence uses its Davidson Ritz-change schedule, Broyden mixing and independent verification. Historical buffer-band recovery is oracle evidence; it is not a production tuning recipe.

Reason: The old EMPTY_LOOSE/BUFFER_TOL controls and legacy standard path must not be presented as the current convergence fix.

### README.md line 104

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| Physics | Fermi-Dirac occupations and chemical potential | Unpolarized and collinear finite-temperature SCF; analytic synthetic oracles and an Al smoke run, no QE metal comparison. Check empty-band convergence; finite-T restart unavailable. [Evidence](tests/reference/occupations/summary.json). |
```

Proposed text: Fermi–Dirac occupations and chemical potential for supported unpolarized/collinear models, with limited validation and empty-band checks. Reference fixed-cutoff electronic steppers support exact checkpoint v2; warm imports and legacy v1 FD restart remain unavailable. Free-moment PBE uses the G4D oracle fallback.

Reason: R8b added FD exact continuation. The original plan’s line-104 “still true” note predates that implementation; see tests/reference/checkpoint/reference-v2-results.json.

### README.md line 109

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| SCF | Standard convergence profile with fixed-point verification | Default for every command since 2026-11-16 (ADR 0043); strict is a reference tool only. The stored P17A evidence measured the unverified profile: less solver work in successful cases, skew/slab/mesh failures recorded. [Evidence](tests/reference/perf/p17a/final-report.json). |
```

Proposed text: Reference driver is the single production standard path with independent verification (ADR 0047). StrictOracle remains for explicit validation and the G4D exception. P17A numbers remain historical unverified legacy-profile evidence.

Reason: “Default for every command” and “strict is a reference tool only” omit the current G4D routing exception; historical P17A data cannot label reference performance.

### README.md line 111

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| Physics | Fixed-cell ionic relaxation | BFGS/FIRE, Cartesian constraints, warm starts and byte-identical restart. CLI spin-polarized and finite-temperature relaxation unsupported; Gamma-only MgO is not bulk qualification; comparison with QE pending. [Evidence](tests/reference/relaxation/manifest.json). |
```

Proposed text: BFGS/FIRE and Cartesian constraints with reference SCF and live warm seeds. Optimizer resume restores optimizer state but currently lacks the reference orbital seed and records a fresh atomic-start fallback; historical byte-identical resume evidence belongs to the legacy path. CLI spin/FD relaxation remains unavailable.

Reason: Current relax.rs does not persist reference_seed in restart.json; do not advertise bit-identical reference optimizer continuation.

### README.md line 112

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| Physics | Cell relaxation | Works for the stated synthetic cases and isotropic PseudoDojo Si/MgO shifted-mesh runs, with limits: rebuild plus fresh-basis polish verifies start agreement and final pressure. Full-tensor symmetry preservation on shifted meshes and cutoff convergence are not established; comparison with QE pending. [Evidence](tests/reference/relaxation/cell/validation.json). |
```

Proposed text: Variable-cell relaxation uses the reference driver with frozen Miller/FFT support and a final nominal-cutoff rebuild. Historical synthetic/isotropic Si/MgO rebuild-plus-polish evidence retains its original driver; gen19 reports reference/legacy lattice difference 6.3e-8 Bohr and pressure difference 2.3e-9 Ha/Bohr^3. Full-tensor shifted-mesh symmetry and cutoff convergence remain unestablished.

Reason: R7 implemented the port. The old table is not fresh reference-driver evidence and does not prove cutoff convergence.

### README.md line 118

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| Speed | Optional Intel MKL backend for FFT and matrix products, chosen at compile time | Merged behind a compile-time feature; the default build stays pure Rust. Kernel-level measurements on the test server look promising; end-to-end numbers follow. |
```

Proposed text: Optional Intel MKL FFT/matrix backend is implemented behind a compile-time feature; the default build is pure Rust. See dated same-model dev08 tables for end-to-end measurements and their limits.

Reason: End-to-end measurements already exist; this is no longer merely a promising kernel-level development item.

### README.md line 119

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| Speed | RMM-DIIS eigensolver for large band counts and Kerker-preconditioned mixing for hard cases (merged as switchable modules); cheaper FFT call structure (in progress) | Defaults are unchanged; a default changes only with evidence on several systems. |
```

Proposed text: Reference RMM and reference ladder remain isolated opt-in algorithms; reference preconditioning is plain or TF. Adopted reference rollback branches and rejected experiments were retired under D3; retained choices are listed in the current switch manifest.

Reason: Avoid describing old Kerker/Pulay controls as the current reference mixer, or implying every historical switch remains live.

### README.md line 122

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| Runtime | Observable, controllable calculations for people and AI agents: step-by-step SCF, live events, safe-point policy changes, checkpoint/restart, self-healing SCF demo, MCP server, Python/ASE | Design in progress (docs/research/agent-native-directive-review.md); no runtime code yet. The numerical core stays free of LLM and network dependencies. |
```

Proposed text: Reference stepper/runtime, observed and supervised CLI, checkpoint v2, Python stepping and host-embedded MCP are implemented with documented limits. Legacy budget-recovery evidence remains historical; reference budget-extension recovery is unsupported. Numerical core has no LLM/network dependency.

Reason: R8a–R8d make “no runtime code yet” false; do not claim general self-healing from the historical demo.

### README.md line 132

Status: applied (D6 approved 2026-11-26).

Current text:

```text
| Speed | GPU backends and MPI | Design study in progress; the CPU path stays the reference. |
```

Proposed text: MPI is implemented behind a build feature; its distributed reference path is validated under documented gates. GPU remains a design topic with no implementation.

Reason: ADR 0041/M2/Stage R implement MPI; GPU and MPI no longer share an unavailable status.

## Generated validation block: provenance follow-up, not replacement data

Status: applied (D6 approved 2026-11-26).

Do not replace stored numerical values with the gen14–gen28 summaries. The
existing QE table, v1 checkpoint counts and P11 self-healing numbers retain their
original inputs, drivers and limitations. An approved generator change should
label those drivers and add separately stored reference v2/runtime evidence.
In particular, the v1 checkpoint paragraph and budget-recovery demo do not
establish v2 runtime-audit restoration or reference budget recovery. The NLCC
checkpoint-unavailable statement in that block is historical v1 scope, not a
blanket statement about v2. The original R12-DOCS ticket did not own the generator or new evidence.
R12B-README owns the generator follow-up: historical drivers are now labelled,
and separate checkpoint v2/runtime sections read the existing evidence files.
Stored values are preserved, including the v2 report’s BLOCKED status and
explain failure; no new measurements are claimed.

## Additional review boundaries

Status: applied (D6 approved 2026-11-26). All four boundaries retained.

- Line 61’s deterministic-threading statement is a same-build/backend contract,
  not cross-platform or MPI rank-count bit identity; the proposed line-97 text
  supplies the scope for the status table.
- Reference `restart` writes an export container, not normal result files.
  Add this distinction if a new checkpoint feature row is approved.
- Leave magnetic-ground-state, dataset transferability, cross-code force/stress,
  spin/FD-relaxation and cutoff-convergence limits intact. Successful Fe/O runs
  do not establish an independent magnetic-material accuracy comparison.
- A future current-performance paragraph must use same-model comparisons only;
  this list supplies no new speed ratio or QE reference calculation.
