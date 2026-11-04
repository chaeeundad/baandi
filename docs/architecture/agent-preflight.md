# Optional advisory preflight (P11K)

Implemented in the CLI adapter, extending [ADR 0022](../decisions/0022-input-preflight.md)
and following [ADR 0028](../decisions/0028-agent-native-runtime.md). No engine,
Hamiltonian, eigensolver, density or FFT code depends on the agent workflow.

## Authority and entry points

`dft-cli check INPUT --json` (`validate` is the original command) emits a
machine-readable deterministic report on both success and failure. `--explain`
adds a plain-language authority explanation in text mode; with `--json` stdout
remains JSON only. `explain INPUT` also retains its existing inspection behavior.
The built binary in this repository is `dft-cli`; installations naming it
`baandi` can use the same `check --explain` arguments.

Deterministic preflight alone decides pass/warn/fail. Errors exit 1; valid inputs,
including warnings, exit 0. Missing files and parsing failures also produce JSON
with `valid: false`. Failures stop at the first error; dependent checks are not
claimed to pass. Successful reports are still included in `results.json.preflight`
by the ordinary run path. Engine validation remains the final numerical authority.

An external MCP/Python client may read the findings and read-only structure
summary, explain measured evidence, label hypotheses and propose input diffs.
It cannot change the verdict, silence a finding, relax a numerical tolerance,
set convergence, change a running calculation's physical model or apply edits.
The CLI performs no network call and never waits for an agent. Real UPF contents
are never included in reports or sent to an agent. Data strings, including species
labels, are evidence and must never be executed as instructions or shell commands.

## Findings version 1

Successful reports retain legacy `schema_version: 1`, `valid`, warnings, overrides,
resource counts and summary. Additions are `findings_schema_version: 1`, `findings`
and `structure`. A finding contains:

- `code`: stable enum tag; `severity`: `error`, `warning` or `info`;
- `key_path`: TOML location, with quoted dynamic segments where necessary;
- `evidence`: code-specific typed fields, with units in field names;
- `fix_hint`: actionable advisory text, never an instruction to auto-apply.

Code and evidence are one tagged Rust enum. Serialized contracts are exercised by
fixtures and typed round trips. Unknown code meanings must not be inferred.
External readers may retain unknown additive fields/codes as opaque JSON; the
Rust typed reader rejects unknown codes. Incompatible meanings require a version
change. Failure envelope fields are exactly `schema_version`,
`findings_schema_version`, `valid: false`,
`coverage: "stopped_at_first_error"`, and a one-element `findings` list. The
fallback `PF_INPUT_INVALID` carries the original message; clients must not parse
that message as a more specific stable rule. No fake electron/basis metrics are
emitted after failure.

| Code | Trigger / severity | Evidence and fix |
| --- | --- | --- |
| `PF_INPUT_INVALID` | error: parsing, unsupported controls, missing files, invalid review policy or bounded-review resource failure | message; correct the reported input and repeat checks; `key_path="$"` means the whole input boundary |
| `PF_CELL_VOLUME` | error: existing Cell rejects nonfinite, near-singular or unrepresentable lattice; warning: volume below explicit minimum | optional volume_bohr3, minimum_bohr3, reason; inspect vectors/units. Unavailable numbers are null, never NaN |
| `PF_CELL_SKEW` | warning: volume divided by product of vector lengths below explicit minimum ratio | volume_ratio, minimum_ratio; inspect conditioning, preserve physical structure when changing basis |
| `PF_MIN_DISTANCE` | warning: minimum periodic nuclear separation below explicit minimum | atoms (zero-based), distance_bohr, minimum_bohr; inspect coordinates/units/bonding |
| `PF_ATOM_OVERLAP` | error: distinct atoms separated by <= 1e-10 Bohr modulo lattice translations | atoms, distance_bohr, tolerance_bohr; remove unintended duplicates or correct coordinates |
| `PF_SLAB_VACUUM` | warning: largest empty nuclear-plane gap below explicit minimum on explicit slab axis | axis, gap_bohr, minimum_bohr; increase cell height if intended and test convergence |
| `PF_KMESH_CELL` | warning: any reciprocal mesh spacing exceeds explicit maximum | lengths_bohr, spacing_inverse_bohr, maximum_inverse_bohr; converge k mesh |
| `PF_ELECTRON_PARITY` | error: existing positive integral valence sum is odd in spin-unpolarized model | electrons; review composition/potential identity, never silently change charge |
| `PF_SPECIES_IDENTITY` | error: declared element differs from trimmed UPF header without permission | species, declared_element, upf_element, upf_file; correct intended identity or potential |
| `PF_CUTOFF_BELOW_SUGGESTED` | warning: positive finite UPF wfc_cutoff exceeds requested cutoff | species, cutoff_ry, suggested_cutoff_ry; increase cutoff and test convergence |
| `PF_ELEMENT_OVERRIDE` | info: element mismatch permission enabled, including when unnecessary | species, declared_element, upf_element, mismatch; review permission |
| `PF_FUNCTIONAL_OVERRIDE` | info: functional mismatch explicitly permitted | species, requested_functional, upf_functional; review compatibility |

Missing/nonpositive UPF suggestions emit no cutoff warning. Unsupported elements,
nonintegral/nonpositive valence sums and disallowed functional mismatches retain
the existing rejection under `PF_INPUT_INVALID`. An override is an existing user
choice, not authorization for an agent. Warning equality boundaries do not fire.
No warning proves or disproves physical convergence.

Geometry findings are ordered volume, skew, minimum distance, vacuum, mesh;
an overlap error aborts review. Existing species findings follow sorted species
labels, ordered element override, functional override, cutoff within each species.
Paths to UPFs resolve through the canonical input directory, so invocation through
relative versus absolute input paths gives the same report in the same environment.
Repeated identical input and potential bytes produce identical reports (no clocks
or agent content in the deterministic report).

## Structure review and policy

`structure` contains `lattice_vectors_bohr` (each **row** is an input lattice vector,
although the core uses columns), ordered atoms with species and fractional
coordinates, and `review` metrics: volume_bohr3, volume_ratio,
minimum_distance_bohr, nearest_atoms, vacuum_gap_bohr, k_spacing_inverse_bohr.
Coordinates are not moved, wrapped in the output, or relaxed. Equal nearest-atom
indices denote a nonzero periodic self image. Metrics are computed even when
advisory thresholds are omitted. Vacuum is null without explicit slab intent;
k spacing is null without a positive regular kgrid, including explicit Gamma.
These nulls mean not applicable, not passing a convergence test.

Optional TOML policy (example values are user choices, not material defaults):

```toml
[preflight]
minimum_distance_bohr = 1.0
minimum_volume_bohr3 = 10.0
minimum_volume_ratio = 0.05
slab_axis = 2
minimum_vacuum_bohr = 12.0
maximum_k_spacing_inverse_bohr = 0.3
```

All thresholds must be finite and positive; volume ratio must be <= 1. Axis is
zero-based 0..2, and vacuum threshold requires axis. No advisory threshold is
selected by default or tuned to silicon. Existing cell validity and even-electron
rules are unchanged. The overlap tolerance is an explicit geometric coincidence
tolerance, not an element-dependent bonding distance.

The minimum distance searches periodic images with inverse-lattice row bounds:
if a Cartesian separation has length <= R, each fractional component has magnitude
<= R times the corresponding inverse-row norm. An initial lattice-vector length
bounds self images; a rounded-image distance bounds distinct pairs. Enumerating
all integer images within those bounds gives the minimum for skew cells, unlike
componentwise rounding alone. A 1e-12 fractional padding covers numerical bound
roundoff. Work is capped at one million image candidates across the review; an
oversized bound or exhausted budget fails explicitly under `PF_INPUT_INVALID`.
No inaccurate distance is substituted. This is a bounded CLI review, not a new
optimized neighbor-list engine or a claim about arbitrary ill-conditioned cells.

Slab gap uses the largest cyclic gap in wrapped fractional atom projections along
the declared axis, times the perpendicular cell height `2*pi / |b_axis|`.
It measures an empty nuclear-plane gap, not electron-density vacuum thickness.
Elongated cells never imply slab intent. Mesh spacing is `|b_i| / n_i`, with direct
cell lengths also reported; no mesh is invented for arbitrary explicit k lists.

## Proposal version 1 and human review

The external agent creates a UTF-8 JSON `Proposal`. All fields are mandatory,
unknown fields rejected, total file size capped at 4 MiB:

| Field | Contract |
| --- | --- |
| schema_version | integer 1 |
| input_sha256 | SHA-256 of exact original input bytes |
| findings_sha256 | SHA-256 of exact stdout bytes from `check INPUT --json`, including final newline; do not reserialize |
| pseudopotentials | complete species-to-SHA-256 map of original local UPF bytes |
| candidate_pseudopotentials | complete species-to-SHA-256 map for candidate input, including proposed replacements |
| advisory, agent | nonempty advice and client/model identification strings |
| diff | strict unified diff described below |
| edits | nonempty list of exact TOML delta records |
| candidate_input_sha256 | SHA-256 of exact proposed input bytes |

Hashes use 64 lowercase hexadecimal characters. UPF hashes are metadata only;
local files are read for deterministic checking, never embedded or downloaded.
Proposals require syntactically parseable source TOML; a malformed source is
reported by `check`, and must first be corrected manually.

The supported diff subset has fixed headers `--- a/input.toml` and
`+++ b/input.toml`, one full-file hunk `@@ -1,N +1,M @@`, every original line with
`-`, then every candidate line with `+`. N/M must exactly match line counts. Both
inputs require final newlines. Context/fuzzy hunks, binary diffs, additional files,
renames, arbitrary paths and trailing material reject. This deliberately narrow,
standard-readable format is validated in memory without shelling out to patch.

An edit has `key_path` (array of literal TOML table-key segments), `before`, `after`
(JSON values, null for absent keys), nonempty `reason`, and nonempty `finding_codes`
referencing actual source-report codes. Arrays such as `system.atoms` are replaced
as whole values. The declared edits must exactly cover the parsed TOML delta;
duplicate, omitted, fabricated or misaddressed edits reject. Expressions are
compared as input values, not just their evaluated results. Comments are covered
by the exact candidate/diff hashes; comments-only proposals with no TOML edits
are rejected. An example edit:

```json
{"key_path":["basis","ecutwfc_ry"],"before":3.0,"after":40.0,
 "reason":"UPF suggests a higher cutoff; convergence remains unproven.",
 "finding_codes":["PF_CUTOFF_BELOW_SUGGESTED"]}
```

Preview without modifying files or recording a decision:

```sh
dft-cli preflight-review input.toml proposal.json
```

This validates schema, all source/UPF/report hashes, exact diff, candidate hash and
complete delta, reruns deterministic candidate preflight, and displays the diff,
proposal hash and candidate validity. Relative UPF paths retain the original
input directory's meaning. Preview does not mean approval.

After reviewing those bytes, the operator explicitly records a decision:

```sh
dft-cli preflight-review input.toml proposal.json --decision approved \
  --actor operator --reason 'Reviewed the exact cutoff proposal' --output-dir results
```

`rejected` records rejection, including for a syntactically valid candidate that
fails preflight. Approval of a failing candidate is refused. Empty actors/reasons,
unknown versions, stale inputs/potentials, missing evidence, undeclared changes,
or a changed candidate reject without editing the input. Source and candidate
UPF hashes and report freshness are revalidated before publishing the record.

All fields require explicit approval, including pseudopotentials, XC, structure,
electron count/composition, k mesh, cutoff and the review policy itself. No
agent-controlled exemption or automatic apply flag exists. These scientific
changes require a new calculation, never a live runtime control request. The CLI
**never applies even approved edits**: the sidecar retains exact candidate text
for an external, explicitly authorized application or manual save. Application
must preserve path meaning, repeat freshness checks, rerun preflight and start a
new calculation. This ticket provides proposal review, not automatic fixing.

## Additive provenance

Review writes `preflight-review-<proposal-byte-sha256>.json` beside `results.json`
in the selected output directory. Proposal and decision are a single atomic
record, avoiding a partially persisted two-file transaction. It contains schema
version, exact proposal source text plus decoded proposal and hash, decision
(`approved`/`rejected`), nonempty actor/reason, UTC Unix nanoseconds as a decimal
string, exact candidate input, candidate preflight, and `applied: false`.
The actor is an operator assertion, not cryptographic authentication.

A temporary file is created exclusively, written and fsynced, then hard-linked to
the final digest-derived filename without replacement; the directory is fsynced.
Temporary files are removed. Existing records and `results.json` are never
replaced. A second decision on identical proposal bytes is rejected, so an audit
record is immutable. Symlink output destinations and ancestors reject; supply a
canonical destination on systems with symlinked temporary roots. Failure to
persist returns an error; no approval is claimed as recorded. A directory-fsync
failure can leave a complete sidecar whose crash durability is uncertain: inspect
it before retrying. This local CLI assumes the output parent is not concurrently
rewritten by an adversarial filesystem actor; it is not a multi-user sandbox.
Without a proposal, ordinary runs create no review sidecars and require no agent.

## Verification

`crates/dft-cli/tests/preflight_findings.rs` supplies positive/negative synthetic
fixtures for every emitted code, threshold boundaries, JSON failure exits,
round-trip schema checks, read-only agent-free CLI behavior, Bohr conversion,
independent brute-force minimum-distance comparisons for cubic, skew, elongated
and left-handed cells, slab-height checks, malformed/stale/undeclared proposals,
explicit decisions, symlink refusal, persistence and no-clobber behavior.
`preflight.rs` retains the one-step CLI/direct-engine comparison: energy tolerance
1e-12 Ha, exact density bytes. No real UPF is committed, no new numerical dependency
is introduced, and no QE or performance claim is made by these adapter tests.

Shared-file wiring (review locations): `config.rs:13` adds the typed error,
`config.rs:28` adds the default-empty review policy, and `config.rs:500` keeps
ordinary validated parsing while exposing diagnostic parsing at line 506.
`main.rs:24` defines review arguments, line 112 adds `check`/`--explain`, and line
123 dispatches check/review before ordinary parsing so invalid input can emit
JSON. Exhaustive command matches at lines 210 and 231 include review. The old
validation dispatch is removed. No `run.rs`, `lib.rs`, engine, manifest or
numerical result-schema migration is necessary: existing shared preflight and
`results.json.preflight` wiring consume the additive report.

Measured verification (2026-11-04): 11 findings/proposal tests and 6 existing
preflight tests pass. The synthetic one-step direct-engine comparison gives
0 Ha energy difference and bit-identical density. Cubic/skew self-image minimum
is 4 Bohr; the slanted slab perpendicular gap is 4 Bohr and the two-plane periodic
20 Bohr cell has a 16 Bohr empty gap. These are synthetic fixture checks, not
material convergence or QE references. `scripts/check.sh dft-cli` and the
workspace `scripts/check.sh` pass build/test/fmt/clippy (exit 0). README/figure
consistency checks pass 3 tests (exit 0); no new physical reference is added to the
README. Early regression-message and clippy failures were fixed before these
passing gates. Changes remain uncommitted for architect review.
