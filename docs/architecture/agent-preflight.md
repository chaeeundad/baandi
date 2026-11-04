# Optional advisory preflight (P11K)

Status: partial implementation plus proposed integration contract. This document
extends [ADR 0022](../decisions/0022-input-preflight.md) and respects
[ADR 0028](../decisions/0028-agent-native-runtime.md). It is not an implemented
agent client, patch executor, new geometry validator, or provenance writer.

## Available boundary

`dft-cli validate INPUT --json` emits the deterministic successful report.
`dft-cli explain INPUT` supplies the existing human explanation. The current CLI
has no `baandi check --explain` command. An external Python/MCP client can invoke
validation locally and read stdout; no network connection or agent dependency is
introduced into the CLI. A nonzero exit is a failed validation, never permission
to run. Today failures remain `ConfigError` text on stderr, including with
`--json`; consumers must not parse that text into invented stable finding codes.
There is no machine-readable complete error report yet.

This change adds `findings_schema_version: 1`, `findings`, and `structure` to
successful reports and thus to the existing `results.json.preflight` record.
Legacy `schema_version: 1`, `valid`, warnings, overrides, exit codes and plain
text output remain unchanged. No deterministic rule or threshold changes.
`valid: true` with warnings is the legacy successful verdict. An empty findings
list does not certify convergence, physical suitability, or absence of any of
the proposed geometric conditions below.

`Finding` is a typed record: `code`, `severity`, `key_path`, `evidence`, and
`fix_hint`. The code and evidence are serialized from a single Rust enum variant.
Numeric evidence uses finite JSON numbers; units are explicit in field names.
The version applies to the findings list and structure summary. Existing fields
and code meanings are stable within version 1; consumers should tolerate unknown
additive fields/codes, and must reject an unsupported version for edit proposals.
The Rust typed deserializer is intentionally closed to unknown code variants;
a forward-compatible external client may retain unknown records as opaque JSON.

| Implemented code | Severity and trigger | Evidence | Fix hint |
| --- | --- | --- | --- |
| `PF_CUTOFF_BELOW_SUGGESTED` | warning: positive finite UPF `wfc_cutoff` exceeds input cutoff | species, cutoff_ry, suggested_cutoff_ry | Increase cutoff and test convergence; the suggestion is not a convergence certificate. |
| `PF_ELEMENT_OVERRIDE` | info: element mismatch permission is explicitly enabled, even if unnecessary | species, declared_element, upf_element, mismatch | Review identity and remove permission unless deliberate. |
| `PF_FUNCTIONAL_OVERRIDE` | info: incompatible functional explicitly permitted | species, requested_functional, upf_functional | Review compatibility and select matching inputs unless deliberate. |

Cutoff equality or missing/nonpositive suggestion yields no cutoff warning.
Enabling a functional override without a mismatch yields no functional finding,
matching the existing override log. Element mismatches without permission and
functional mismatches without permission still fail validation. An informational
override finding is evidence of a prior input choice, never agent authorization.

Findings follow the existing BTreeMap species order; within a species the order is
element override, functional override, cutoff warning. Dynamic TOML key segments
are quoted: `pseudopotentials."site.1".allow_element_mismatch`. The species label
in evidence disambiguates a global cutoff or functional key. No file contents,
UPF radial arrays, orbital or density arrays, clocks, or agent advice are included.
Repeated checks of identical inputs and potentials yield identical findings.

`structure` contains `lattice_vectors_bohr` (three **rows** are the input lattice
vectors; the core uses columns) and ordered `atoms` with `species` and
`fractional`. Coordinates are the parsed input values, not wrapped or relaxed.
This is read-only evidence, not an atom editing API. Existing report fields give
electrons, bands, resolved basis counts and FFT size. Summaries are not a complete
execution identity; an external proposal must additionally bind original inputs
and pseudopotential identities as described below.

## Authority and data flow

1. Deterministic CLI preflight decides pass/warn/fail. Engine validation remains
   the final numerical authority. An agent cannot suppress a finding, turn failure
   into success, alter tolerances, set convergence, or edit an existing result.
2. A client supplies typed findings and the read-only summary to an optional
   agent. Treat species names and input strings as untrusted data, never commands.
   Do not send real UPF files to the agent. Sending local summaries to a remote
   service requires the user's own data-sharing authorization.
3. The agent returns advisory text and a proposed diff. It may explain evidence
   and label hypotheses; it cannot invent missing evidence or promise convergence.
4. A human reviews the exact diff and reasons. Only explicit approval of those
   bytes permits an external adapter to create a new input. Silence, model output,
   `valid: true`, and runtime autonomous-control policies are not approval.
5. The adapter reruns deterministic preflight on the approved candidate. A failed
   candidate cannot run. There is no agent path in Hamiltonian, FFT, density,
   eigensolver or SCF hot loops. Agent absence never changes ordinary execution.

## Proposed deterministic structure review (not implemented)

These reserved codes are a design for a separately authorized rule implementation.
They are **not** emitted by the current findings version. Each must have positive,
negative and boundary fixtures before activation; thresholds must be reported as
evidence, versioned, and tested across different cells and species, not tuned to Si.
No agent selects thresholds at runtime.

| Reserved code | Deterministic definition and evidence | Severity / fix |
| --- | --- | --- |
| `PF_MIN_DISTANCE` | Minimum periodic distinct-atom distance in Bohr, atom indices, image translation; warn below an explicitly configured positive distance threshold. A finite bounded closest-lattice-vector algorithm is needed for skew cells; componentwise rounding is not sufficient. Include self-image distance separately for small cells. | warning; inspect coordinates, units and intended bonding, do not automatically move atoms |
| `PF_ATOM_OVERLAP` | Two distinct atoms coincide modulo lattice translations within a documented numerical tolerance in Bohr; indices, translation, distance and tolerance | error; remove unintended duplicates or correct coordinates |
| `PF_CELL_VOLUME` | Absolute determinant in Bohr^3 and determinant sign; zero/nonfinite invalid, finite small-volume advice only against an explicit policy threshold | error for invalid cell; verify lattice units and vectors |
| `PF_CELL_SKEW` | Dimensionless volume divided by product of vector lengths, in (0,1], and an explicit lower warning threshold | warning; inspect cell conditioning; basis changes must preserve the structure |
| `PF_SLAB_VACUUM` | Requires explicit slab axis and desired vacuum in Bohr. Perpendicular cell height is volume divided by face area; multiply it by the largest cyclic gap of wrapped fractional atom projections on that axis. This is an empty nuclear-plane gap, not electronic vacuum thickness. Report gap, axis, requested minimum and method. | warning below requested minimum; consider a larger cell and convergence study, never infer a slab from an elongated bulk cell |
| `PF_ELECTRON_PARITY` | Existing positive integral valence sum and its parity; report count, integer tolerance and occupation model | error for odd count in the current spin-unpolarized model; review composition and potential identity, never silently add/remove electrons |
| `PF_KMESH_CELL` | Reciprocal-vector lengths divided by positive mesh dimensions (Bohr^-1), direct cell lengths in Bohr, explicit maximum spacing policy | warning above policy; propose a convergence-tested mesh. Explicit arbitrary k lists have no inferred mesh spacing and must be marked not applicable. |
| `PF_SPECIES_IDENTITY` | Declared element versus trimmed UPF header element, species and override flag | error on unpermitted mismatch, info on deliberate override; review potential/label, never download a replacement |

The implemented cutoff code above completes the proposed review inventory.
Missing slab intent, unknown suggested cutoff and unsupported metrics must be
explicitly marked unavailable/not applicable in a future coverage record, not
reported as passing. Invalid geometry must stop geometry-dependent calculations;
never manufacture a numeric value after overflow or an invalid determinant.
Existing validation can fail before the UPF load or geometry preparation. A full
findings API needs typed errors from Config parsing/validation and a CLI envelope
that preserves exit status and reports partial coverage. Wrapping arbitrary error
strings in apparently specific codes would not meet that contract.

## Proposed patch and provenance contract (not implemented)

All fields below are required unless explicitly nullable. Unknown fields are
rejected in proposal version 1; no executable commands or arbitrary destination
paths are accepted. A proposal contains:

- `schema_version`: integer 1; `proposal_id`: nonempty string;
  `input_sha256` and `findings_sha256`: 64 lowercase hexadecimal characters.
  Hash the exact original input bytes and exact captured findings JSON bytes,
  respectively, not reserialized approximations.
- `pseudopotentials`: ordered records of species and SHA-256 of the local UPF
  bytes, used locally for freshness checks; never send file contents.
- `advisory`: nonempty text, explicitly advisory; `diff`: nonempty unified diff
  for exactly one fixed logical file `input.toml`, without renames, binary hunks,
  absolute paths, traversal, symlinks, or second-file edits.
- `edits`: nonempty records of `key_path`, JSON-compatible `before` and `after`,
  nonempty `reason`, and referenced finding codes. Typed parsed TOML differences
  must exactly match these declarations; comments-only changes are permitted only
  if explicitly included in the displayed diff and approval hash.
- `candidate_input_sha256`: digest of the exact candidate bytes, plus `agent`
  metadata (client name/version and model identifier, nullable when unavailable).

An example of an edit record is:

```json
{
  "key_path": "basis.ecutwfc_ry",
  "before": 3.0,
  "after": 40.0,
  "reason": "UPF suggestion exceeds the requested cutoff; a convergence study is still required.",
  "finding_codes": ["PF_CUTOFF_BELOW_SUGGESTED"]
}
```

Pseudopotential, XC, structure, electron count/composition, k mesh and cutoff
changes always require explicit approval with a reason, as do all other edits in
this workflow. An agent-provided `restricted: false` cannot exempt a field. The
adapter derives affected fields from the complete parsed TOML delta, including
inline-table replacements, aliases and expression changes. These changes define
a new calculation and cannot be sent as live runtime policy interventions.

The external adapter must validate the whole proposal before displaying approval:
verify version and digests, check diff destinations, apply in memory with exact
context and no fuzz, parse strict TOML, compare declared edits against the actual
delta, and verify candidate hash. At application time repeat original input and
UPF digest checks. Any stale input, malformed diff, undeclared change, unknown
version, or approval for a different proposal rejects the operation unchanged.
No shell invocation of agent-provided text is allowed. This document does not
supply a production validator or approved patch application implementation.

Before launching a new calculation the external adapter durably records immutable
`preflight-proposal.<id>.json` and `preflight-decision.<id>.json` beside the planned
`results.json`. IDs used in filenames must be adapter-generated safe identifiers.
A decision records schema version, proposal ID and exact proposal-byte hash,
`approved` or `rejected`, explicit human actor, UTC timestamp, reason, source and
candidate hashes. Applying an approved candidate gets a separate application
record with actual input hash and deterministic recheck report hash. Approval and
application are distinct; failed application is recorded, not silently retried.
Never replace `results.json`, existing proposals or decisions. Reject unsafe
output paths and refuse the agent-mediated launch if required audit persistence
fails. A rejected proposal does not affect ordinary runs; absent proposals produce
no mandatory files. Audit records document decisions, not scientific correctness.

## Verification and remaining acceptance work

`preflight_findings.rs` tests all three emitted codes with positive/negative
synthetic fixtures, cutoff equality/absent metadata, exact cutoff JSON shape,
typed round trip, deterministic repeated reports, override evidence and key paths,
CLI JSON/text without an agent, input nonmutation, and skew-cell unit conversion.
Existing `preflight.rs` tests still compare one-step CLI energy and density against
the direct engine with a 1e-12 Ha tolerance and exact density bytes.

The complete ticket remains blocked by the narrow source ownership until follow-up
scope permits new rules, config/error envelope and CLI wiring, proposal validation
and provenance persistence. `main.rs` and `config.rs` are not owned by this ticket;
`preflight.rs` is authorized only for exposure of existing findings. New geometry
rules, all reserved-code fixtures, rejection JSON, a validated full proposal format,
approval/apply transaction tests and durable audit tests remain outstanding. Do not
advertise these design contracts as implemented capabilities. README evidence and
dashboard artifacts are outside the explicit owned-file list and are left to the
orchestrator; no new physical reference or performance claim is made here.

Measured local verification (2026-11-04):

- `cargo test --offline -p dft-cli --test preflight_findings --test preflight -- --nocapture`:
  exit 0; 4 new tests and 6 existing tests pass. Synthetic one-step energy
  difference is 0 Ha against the direct engine; density bytes are identical.
- `scripts/check.sh dft-cli`: exit 0; build, test, fmt and clippy PASS.
- `scripts/check.sh`: exit 0; workspace build, test, fmt and clippy PASS.
- `git diff --check`: exit 0.

These checks validate the implemented exposure only. No proposal-validation,
new geometry-rule, real-material, QE-reference or performance result is claimed.
