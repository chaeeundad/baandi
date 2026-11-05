# State contracts for staged materials support

Status: **proposed**, P13A-STATE-CONTRACTS / MF02, 2026-11-05. Documentation and
schema proposal only. Scope B is authorized in [AGENTS.md](../../AGENTS.md);
these contracts still require user approval before MF03. No feature, external
reference run or derivative implementation is delivered here.

## Authority and decision coverage

The [master plan](../roadmap/master-plan-proposal.md) supplies D1–D12 and
MF02–20; [feature matrix E1–E8](../research/dft-feature-matrix.md) and
[gap survey T4/T5/T8](../research/feature-gap-survey.md) motivate the boundaries.
Those surveys are historical snapshots: [evidence map](../validation/evidence-map.md)
and [checkpoint-format](../checkpoint-format.md) govern current capability claims.
[Physics conventions](../physics-conventions.md) remain unchanged. Every rule
below is a proposed choice under the cited decision, not an assertion of support.

| Decision | Concrete disposition in this proposal |
| --- | --- |
| D1 | [ADR 0031](../decisions/0031-spin-electronic-state.md): tagged spin, explicit axes, XC boundary |
| D2 | [ADR 0032](../decisions/0032-occupations-scf-certificate.md): fixed/FD policy, E/TS/F, fresh finite-T map |
| D3 | MF11 record below; full mesh is the oracle, magnetic reduction disabled initially |
| D4 | NC-first U3; dataset and Hubbard records below reserve capabilities without implementing augmentation |
| D5, D6, D8 | [ADR 0033](../decisions/0033-grids-derivative-ledger.md): named grids, ownership, units, complete derivatives |
| D7 | Assembly rejects unsupported differential/nonlocal XC; do not overload a scalar XC potential. Existing solver operator boundary is retained. |
| D9 | [ADR 0034](../decisions/0034-state-envelope-extension.md): extend concrete P11C, density-first U5 |
| D10, D11 | Per-worker scratch and one scheduler budget; optional MKL remains compile-time ([ADR 0019](../decisions/0019-optional-mkl-backend.md)); MPI/GPU require separate measured need |
| D12 | Immutable typed environment identity and separate inner/outer quality; [ADR 0017](../decisions/0017-coupling-commitments.md) remains proposed, no coupled total or derivative without all owners |

## Shared records

All logical numeric arrays contain finite f64 (complex128 for orbitals), with
checked dimensions and explicit units. IDs are bounded nonempty ASCII strings
up to 128 characters; digests are lowercase SHA-256 hex. Canonical fingerprints
use versioned typed field order, length-prefixed UTF-8/arrays, little-endian u64
counts and finite IEEE-754 f64 bits; normalize negative zero to positive zero
for identity, not stored scientific payload. No decimal rounding or atom/k
sorting before hashing. Store the canonicalization version. This is a proposed
identity rule (D1/D9), distinct from P11C's byte-integrity checksum.

| Record | Complete proposed contents |
| --- | --- |
| `KPointDescriptor` | `global_id:u64`, `fractional_reciprocal:[3]` canonical [-1/2,1/2), `weight:f64`, `basis_order_digest`; model order stable |
| `ModelIdentity` | canonicalization version; cell columns Bohr; ordered atom species and fractional positions; ordered per-species dataset descriptors; total electron count; XC ID/version and explicit mismatch override; spin mode/axis/constraints; full occupation policy including theta/width/order; Hubbard model; symmetry model; environment ID and boundary/G=0 convention |
| `DatasetDescriptor` | SHA-256 of actual bytes, element, valence electrons, XC label, relativistic tag, NC kind, NLCC flag, semicore label, ordered projector channel labels `(l,radial_index)`, source/version/license reference and qualification ID. Provenance/license strings are retained but physical digest uses bytes and physics tags. Unsupported overlap/augmentation/SOC rejects. |
| `BasisIdentity` | cell/grid descriptors, wavefunction/density cutoffs Ha, ordered k descriptors and weights, per-k ordered integer Miller vectors or their verified digest, band layout and orbital normalization; transfer identity when applicable |
| `ElectronicState` | schema version; state/calculation IDs; model/basis/environment digests plus their records; tagged spin density; occupation policy; optional band data and optional complete orbital payload; optional thermodynamic ledger; quality; provenance. Evaluated states require band data and ledger. Density-only unevaluated guesses have neither, explicitly. |
| `StateQuality` | `kind:initial_guess|evaluated`; `status:unevaluated|converged|unconverged|failed`; optional typed certificate with thresholds and measured residuals; inner SCF status; outer status `not_applicable|unconverged|converged`; unavailable-observable reasons. Only evaluated, freshly certified compatible states are rankable. |
| `StateProvenance` | parent state IDs, start recipe ID and parameters/seed, transfer ID or absent, attempt ID, execution/build digest, evidence IDs, migration ID/source checksum or absent. Timestamp and wall durations are audit-only. |
| `AttemptRecord` | attempt ID, model/basis digests, start recipe, execution identity, status, error `{key,code,message,remedy}` or absent, optional resulting state ID; retain failures with no state |

`model_digest` fingerprints ModelIdentity physical fields **and** basis digest
for strict energy-ranking compatibility, following [ADR 0016](../decisions/0016-esle-layering.md).
Separate basis digest permits diagnostic compatibility explanations; a different
cutoff/mesh is a convergence comparison, not automatic state ranking. Policy
parameters and fixed moment are model data; actual occupations and spontaneous
moment are solution data. State digest covers model/basis/environment digests,
authoritative density, evaluated band data including actual occupations, chemical
potentials and entropy, quality and optional orbital payload hash. Distinguish
payload identity from physical equivalence: orbital gauge or retention changes
can change state digest without proving a new branch. SAME/DISTINCT/UNRESOLVED
classification uses measured density/moment noise floors, never hashes/energy
alone (MF09, ADR 0016). Seeds and thread counts do not change model identity;
execution identity controls exact replay.

For MF10, propose `HubbardModel=none|dudarev{sites}` where each site specifies
input atom index, l, ordered m=-l..l, U_eff in Ha, projector source hash,
radial channel ID, normalization and orthogonalization convention, and
`double_counting:fully_localized_limit`. `HubbardState` contains per-site/spin
Hermitian complex occupation matrices `[2l+1,2l+1]`, basis digest and trace.
Proposed E_U=(U_eff/2) sum_s Tr[n_s-n_s^2], with V_U=U_eff(1/2 I-n_s),
requires MF10 approval and independent variational/gauge tests; this formula is
a proposal under MF10, not a claim of matching a code's projector convention.
No U_force capability follows from energy support. Reject nonorthogonal projectors
unless the accepted projector metric and dual definition are supplied.

For MF11, `SymmetryOperation` specifies integer fractional rotation `[3,3]`,
fractional translation `[3]`, atom permutation, spin action `identity|swap`,
antiunitary flag, and model-invariance proof ID. `KOrbitMap` records every full
k ID, representative ID, operation ID, reciprocal integer shift, phase convention
`exp(-i(k+G).tau)`, and multiplicity; weights sum to one. Transform density and
projectors with the same operation and origin. Active spin constraints, fields,
Hubbard projectors and species must be invariant; geometry alone is insufficient.
Initial collinear runs use full meshes. Any later antiunitary/time-reversal use
needs an accepted operation-specific proof and full-mesh tests, including the
spin swap. SOC/noncollinear remains unsupported. No silent cell standardization.

## Results extension and schema boundaries

[results-state.proposal.json](../schemas/results-state.proposal.json) validates
an **extension object**, proposed to live at `results.json.state_contracts` under
an explicitly selected future output profile. It does not pretend to validate
the entire historical results.json, for which the gap survey T5 called for a
separate full contract. Existing schema_version=1 output and keys remain unchanged
until a separately authorized adapter ticket. The extension carries
`proposal_status:proposed`, `extension_version:1`, spin/occupation/thermodynamic
blocks and unavailable force/stress blocks. Schema references are local files;
no network resolution is needed. JSON cannot express runtime array-product,
weighted sums, positivity, digest authentication or energy equations: the importer
must check the invariants in ADRs 0031–34 after schema validation. Reject duplicate
keys, NaN and Infinity before schema validation. No executable engine code or
placeholder outputs are added.

`spin.mode=noncollinear_reserved` and named future smearings are documented in
ADRs but deliberately fail this first executable-profile proposal schema.
Force/stress `status=available` and numerical arrays also fail until a future
version explicitly defines a validated availability branch. Unavailable is not
zero. Fixed-insulator mu is null with `fixed_insulator` reason; FD mu is required.
Three small schema files separate reusable spin, occupations and results envelope.

## Migration and invalid combinations

Apply [ADR 0022](../decisions/0022-input-preflight.md)'s key-path/reason/remedy
style in both validate and run, before expensive SCF allocation; engine entry
points still enforce the same typed rules. Example messages below are proposed
strings, not existing emitted diagnostics. Invalid input retains CLI exit 1;
numerical failure/certificate failure remains distinct per [CLI contract](../cli.md).

| Input or migration | Proposed outcome / message |
| --- | --- |
| Existing Gamma or fixed-k constructor/input without new fields | Same defaults, rejection behavior, arithmetic and first-k result view; additive new entry point only |
| Old result without extension | Valid historical output; no spin/force capability inferred |
| v1 checkpoint direct exact restart | Original strict v1 contract; preserve build/options checks |
| Explicit v1 -> v2 warm migration | Unpolarized fixed-policy mapping, recorded migration/checksum, new calc ID/reset; no automatic exact conversion |
| Collinear -> unpolarized or reverse direct restart | `state.spin: incompatible spin mode; request an explicit separately validated warm transfer` (not available initially) |
| noncollinear/SOC, PAW/USPP, local moment constraints | `model.capabilities: unsupported combination; use supported NC collinear settings or await the named feature ticket` |
| Unpolarized fixed odd/fractional N | `occupations.electrons: fixed unpolarized policy requires even integer N; select a supported explicit spin/occupation policy` |
| MF03 fixed collinear without integer channel counts | `spin.constraint: channel counts required for fixed occupations; specify a fixed-moment sector` |
| FD theta<=0, nonfinite width, unsupported Gaussian/MP/cold | `occupations.policy: unsupported parameter or scheme; use fixed_insulator or implemented positive-temperature Fermi-Dirac` |
| Fixed M with abs(M)>N; fixed channels exceed bands | `spin.constraint: electron counts outside capacity; correct moment or increase bands` |
| Wrong shape, spin capacity, weights, duplicate k ID, NaN | `state.band_data: invalid layout or values; supply explicit matching spin/k/band axes and normalized weights` |
| Legacy scalar view requested for collinear result | `result.spin: scalar legacy view is ambiguous; request spin-resolved data` |
| Spin XC or NLCC not validated for selected functional | `xc.spin: required polarized/core capability unavailable; select a validated functional/dataset combination` |
| Magnetic/field symmetry without model-invariance proof | `symmetry.operations: reduction is not valid for this model; use the full mesh` |
| Any active term lacks force/stress | Query returns unavailable listing terms; a required derivative request fails `observables.forces: incomplete derivative coverage; disable this request or select a validated combination` |
| Unconverged E or FD E supplied to optimizer as variational F | `relaxation.objective: requires certified matching variational quantity; converge the state and request its complete derivatives` |
| Changed cell/grid/band count exact restart | `checkpoint.basis: incompatible exact continuation; create a new calculation using an explicitly supported warm-transfer plan` |
| Wrong checksum/version/layout or incomplete publication | `checkpoint.integrity: invalid or unsupported artifact; retain original and obtain a complete supported checkpoint` |
| Density-only payload requested as exact or result as restart | `checkpoint.kind: requested capability absent; use an explicit warm artifact or an exact checkpoint` |
| Non-neutral/reservoir or density-dependent environment | `environment.model: unsupported thermodynamic functional; fixed-N energy field contracts do not define this calculation` |
| Compare different model/policy fingerprints as same branch | `comparison.model: incompatible fingerprints; report a labeled convergence/model study instead of ranking` |

## Questions for user review

These are proposed choices, not blockers to writing this ticket. Scope B, U3
NC-first and U5 density-first are already defaults; no repeat approval is sought.
Acceptance of this package is the MF03 gate; later tickets must own their files
and produce evidence. Each question has one recommendation and its consequence.

| Question | Recommendation | Consequence / alternative |
| --- | --- | --- |
| Approve ADRs 0031–34 as the MF03 foundation? | Approve after review of these explicit invariants | Allows scoped MF03 implementation; deferring keeps all new physics unavailable |
| First collinear fixed occupations: explicit sector or inferred allocation? | Explicit integer N_up/N_down sector | Deterministic pre-MF04 behavior; unconstrained moment equilibration waits for shared-mu FD |
| First finite-T policy and units? | FD with theta=k_B T in Ha internally, optional boundary Kelvin conversion | Physical S and F; other smearings wait for named generalized corrections |
| First Hubbard convention? | Dudarev U_eff with explicit orthonormalized projectors/FLL, reviewed in MF10 | Small auditable first scope; literature/QE runs must match that projector definition, not just U |
| Stress sign and reporting? | Tensile-positive full Cartesian tensor, pressure=-trace/3 | Avoids Voigt/shear ambiguity; external adapters must convert conventions explicitly |
| Persistence evolution? | Version 2 on P11C with explicit v1 warm migration | Old readers fail closed; exact replay remains version/build specific |
| First chemistry suite under U2? | Qualified NC O/Ni pair, Al metal, and neutral Si vacancy; MgO as nonmagnetic control if qualified | License/NLCC/semicore gates can delay real comparisons; synthetic tests cannot replace qualification |
| How to freeze new material tolerances? | Adopt the proposed caps in stage-oracles, require uncertainty <= one third of each cap before comparison | If not resolved, report inconclusive and refine inputs; never enlarge a cap to fit a result |

## Verification scope

Required commands: `scripts/check.sh` and
`python3 -I scripts/tests/test_readme_figures.py`; actual execution is reported
in the task response. Proposal validation also checks local links and schema
fixtures. No new numerical evidence means no generated README result is added.
Artifacts/session/dashboard files are outside this ticket's OWNED_FILES and
remain the invoking orchestrator's responsibility, as in the master-plan scope
note. All changes remain uncommitted for review.

Verification executed in this worktree on 2026-11-05:

| Command | Exit | Observed scope |
| --- | --- | --- |
| `scripts/check.sh` | 0 | PASS build, test, fmt, clippy; default workspace gates, not all ignored release suites |
| `python3 -I scripts/tests/test_readme_figures.py` | 0 | 3 tests passed; pre-existing unclosed-file ResourceWarnings; generated README/figures unchanged |
| `python3 -I /tmp/p13a_validate.py` | 0 | Temporary offline audit: 3 schema keyword/reference audits, 4 accepted and 20 rejected structural fixtures, 52 local link targets present |
| `git diff --check` | 0 | No whitespace errors |

The temporary audit is limited to the JSON Schema keywords actually used; a
full Draft 2020-12 metaschema validator is unavailable in this offline environment.
It is not a numerical importer or proof of cross-field sums/identities. Its
fixtures cover fixed unpolarized/collinear and FD common/per-spin mu; rejections
include reserved modes/schemes, wrong dimensions/capacities, missing mu, nonpositive
temperature, negative entropy, mismatched objective and fabricated derivatives.
No physical measurement, new QE comparison or performance claim was produced.
Session registry has no P13A entry and no local `.delegate/P13A-STATE-CONTRACTS/report.json`
was present at inspection; wrapper/orchestrator registration is not fabricated.
