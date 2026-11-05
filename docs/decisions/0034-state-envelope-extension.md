# ADR 0034: Extend the P11C envelope for spin and occupations

- Status: **accepted by the user 2026-11-05** (contract approved; implementation capability still absent until its ticket delivers it); P13A-STATE-CONTRACTS / MF02.
- Basis: [master plan D9, MF07–08, U5](../roadmap/master-plan-proposal.md), [scientific-state API](../architecture/scientific-state-api.md), [implemented checkpoint v1](../checkpoint-format.md).

## One container, explicit version

Reuse P11C's `BAANDICP` container, bounded header/binary lengths, contiguous
array traversal, scalar f64 bit encoding, checksum and no-clobber publication.
Propose integer checkpoint `schema_version:2` for new spin/occupation semantics.
Do not confuse it with results.json's version or the broader draft API's string
v1.0. The implemented single-file P11C format supersedes that draft's external
array paths/manifest-last recipe: payload and header are written and fsynced in
an exclusive temporary file; atomic no-clobber publication is the commit point.
No second state store or path-referenced payload format is introduced.

v1 readers reject v2. Future v2 readers retain a separate strict v1 parser; no
unknown-field ignoring. Explicit migration of v1 to a v2 **warm** state records
source checksum and migration ID, maps scalar total density to unpolarized,
fixed occupations, zero entropy, and preserves cell/k/band/orbital ordering.
v1 exact continuation stays on its v1 restore path with matching build/options;
converting it to v2 is not a same-trajectory guarantee. No automatic spin split,
FD occupation invention or exact-to-warm downgrade. Existing destinations and
APIs are unchanged. No migration is implemented by this ADR.

## Proposed header extension

The pure shared records are specified in [state-contracts](../architecture/state-contracts.md#shared-records)
and ADRs [0031](0031-spin-electronic-state.md)/[0032](0032-occupations-scf-certificate.md).
For v2 the electronic payload contains `spin_mode`, `band_layout`,
`occupation_policy`, `spin_density`, `band_data` (optional for an unevaluated
initial guess), `orbital_blocks` (optional), `thermodynamics` (evaluated states
only), `quality`, `model_digest`, `basis_digest`, `environment_digest`,
`state_digest`, and `provenance`. All array descriptors retain v1 offset/count
references, and floats retain the v1 finite bit encoding. Shapes/units/axis
labels are determined by typed fields, not inferred from payload byte counts.

Traversal: spin-density authoritative arrays first, then band-data eigenvalues,
occupations and residuals in [s,k,b] C order, then orbital blocks in spin/k order
(each block retains v1 column-major complex interleaving). Exact continuation
then carries spin-aware mixer inputs/residuals in historical order, plus all
existing history/schedule/runtime fields. Density-first warm archival may omit
band data and orbitals if marked unevaluated; an evaluated export must retain
band data used for its certificate. Large orbital payload retention remains
optional for warm/result states, mandatory for exact continuation if they affect
future arithmetic. Do not infer exactness from a density hash.

| Kind | New required information and promise |
| --- | --- |
| Exact continuation | Next mixed spin density, previous spin/k orbitals, current occupations and any reused mu/root state, complete spin mixer state, spin/FD certificates and histories, edge-band policy and all numerical hints that affect arithmetic; identical model/basis/execution/policy identities. Capture only Ready boundary, with runtime audit/telemetry restrictions unchanged. Unknown opaque memory rejects. |
| Warm start | Validated authoritative spin density, model/basis metadata, occupation policy/constraints; optional complete orbitals. New calculation ID, reset roots/mixer/history/certificate. Direct reuse requires matching model and basis; changed band counts or geometry need separately accepted explicit transfer, never implicit coercion. |
| Result export | Evaluated state, E/TS/F, actual quality including failure, band data and optional orbitals. Inspection only. Certificate is historical evidence, never a restored live convergence flag. |

Import checks finite values, bounded sizes/checked products, supported modes,
channel capacity, exact axis counts/order, positive k weights and sum, per-channel
density feasibility under XC, charge and fixed-M sums, occupation sums, optional
orbital orthonormality, energy identity F=E-TS, fingerprints and checksum before
starting numerical work. Existing relative charge <=1e-8 and orthonormality
<=1e-8 import gates remain upper bounds; stricter policy electron/moment gates
also apply when importing evaluated FD data. A density-only guess does not
fabricate eigenvalues/entropy. Reject unknown versions, duplicate keys, invalid
k/spin layouts and corruption, including valid-checksum semantic inconsistency.

Model digest includes policy definition and constraints; state digest includes
actual occupations, mu and entropy. Execution identity retains P11C's trusted-host
build/platform/thread contract; hashes are integrity checks, not authentication.
FFT/backend handles are rebuilt, never serialized. Density-only warm restarts
promise no trajectory or magnetic-branch identity; measure their noise before
state classification. The [MF07–09 tests](../validation/stage-oracles.md) extend
P11C's independent byte-layout and fault tests instead of replacing them.
