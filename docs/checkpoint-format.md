# Checkpoint format v1 (P11C)

The CLI owns serialization; `dft-engine::checkpoint` owns pure state types shared
with ESLE. No serde, file access, checksum, or checkpoint allocation is added to
SCF iteration code. This implements the three distinct capabilities in ADR 0016
and ADR 0028, rather than inferring continuation from an exported density.

## Commands

```sh
# Stop after three completed iterations (including mixing).
dft-cli checkpoint input.toml boundary.cp --iteration 3 --kind exact
# Rebuild the same immutable problem, restore, finish, and export the result.
dft-cli restart input.toml boundary.cp result.cp
# An explicit initial-guess artifact permits different mixing/solver settings.
dft-cli checkpoint input.toml guess.cp --iteration 3 --kind warm
dft-cli restart different-mixing.toml guess.cp warm-result.cp --warm
# Evaluate a calculation and export its result, including unconverged outcomes.
dft-cli checkpoint input.toml evaluated.cp --kind result
```

Existing destinations are never overwritten. Exact and warm input kinds are not
interchangeable; there is no silent downgrade. An exact checkpoint must be Ready
and cannot capture a terminal or partially failed iteration. Result export needs
an evaluated outcome and preserves its actual convergence flag. It cannot be
used for restart. Restart writes an unconverged result on budget exhaustion, but
returns a failure exit status. All numerical work uses the existing fixed-point
certificate; deserialization never sets a live session to Converged.

## Byte layout

All offsets are byte offsets. A file contains exactly:

| Offset | Size | Content |
| --- | --- | --- |
| 0 | 8 | ASCII `BAANDICP` |
| 8 | 8 | Unsigned little-endian JSON header length H |
| 16 | 8 | Unsigned little-endian binary section length B |
| 24 | H | UTF-8 JSON object, no trailing non-whitespace content |
| 24 + H | B | Contiguous little-endian IEEE-754 binary64 arrays |
| 24 + H + B | 32 | Raw SHA-256 digest of all preceding bytes |

Readers require file length `56 + H + B` with checked arithmetic. The initial
adapter limits the complete file to 512 MiB and its header to 32 MiB. There are no
external payload paths. Truncation, trailing bytes, wrong magic, unsupported
schema, checksum errors, nonfinite floats, array offset/length errors and unknown
or duplicate JSON fields reject. No migrations are implemented.

Publication writes a same-directory exclusive temporary file, fsyncs it, then
creates an atomic **no-clobber hard link** at the destination and fsyncs the
directory. Finally it removes the temporary name. A filesystem without hard-link
support fails explicitly. A crash can leave a temporary file but cannot publish
a partially written destination. A directory-sync error may leave a complete
published file while returning an error; callers should inspect it before retry.
A checksum detects corruption, not malicious alteration or source authenticity.

## JSON and array conventions

The header mirrors the public `Checkpoint` record in
`crates/dft-engine/src/checkpoint.rs`. Its schema is the adapter-owned `WCheckpoint`
and the `wire_struct!`/enum definitions in `crates/dft-cli/src/checkpoint_io.rs`.
`schema_version` is integer 1. Enum unit variants are strings; payload variants
are single-key objects (Serde's external tagging). Options use JSON null. Integers
are unsigned; platform-sized lengths must fit the importing platform.

Scalar floats are their **IEEE-754 bits encoded as JSON u64**, not decimal floats.
This preserves every bit without requiring a decimal round-trip feature in serde.
They must decode to finite values. Fixed small arrays (cell, fractional atoms and
k points) use nested JSON arrays of these scalars. Durations are `[seconds,nanos]`
with nanos below 1e9. System timestamps are `[after_epoch,[seconds,nanos]]`, a
signed offset from the UTC Unix epoch; wall times do not participate in numerical
identity. This concrete v1 representation supersedes the draft's RFC3339 wire
sketch for this container only, not the event API.

Every `Vec<f64>` is an `{ "offset": N, "count": M }` reference into the binary
section. Offset is relative to that section and count is the number of f64
values. Arrays appear in declaration/element order of the typed header traversal,
contiguously, without overlap or padding. Each reference must start at the end of
the previous one, including empty arrays; no unreferenced bytes are allowed.
Reordering JSON object keys does not change the typed traversal order.

- Density: grid `[nx,ny,nz]`, C order, `z` contiguous, electrons/Bohr^3.
- Orbitals: per k point in model order, rows in the existing plane-wave basis
  order, columns in band order. Column-major coefficients interleave real then
  imaginary f64; count is `2 * rows * cols`. Coefficients are dimensionless and
  each orbital has norm one. The reader never reorders or renormalizes them.
- Pulay inputs and residuals: oldest to newest pairs, same grid ordering and
  density units. Each pair is independently referenced in the binary section.
- Energies, eigenvalues and residuals: Hartree; occupations: electron counts;
  cell: Bohr, lattice vectors in columns; k points: three fractional reciprocal
  coordinates followed by normalized weight. History uses the existing engine
  residual and timing conventions without unit conversion.

## Three kinds and restoration

`ExactContinuation` contains next-input density, previous per-k orbitals,
occupations, all completed iteration records and certificates, best/growth/reset
counters, options (including seed, schedules and solver choice), completion mode,
calculation identity, state version, event cursor, and version-1 mixer memory.
Linear stores beta. Pulay stores beta, capacity, mean and ordered input/residual
pairs. KerkerPulay additionally stores previous residual norm, mean and feasibility
diagnostics/counters. FFT/filter scratch is overwritten before use and is rebuilt
from the identical problem; no in-flight eigensolver state is captured. Opaque
external mixer/preconditioner memory rejects rather than disappearing.

The importing host constructs a fresh validated `ScfProblem` with the saved
current options, then calls `ScfSession::restore_exact`. Every option must match,
as must the insulating-gap policy, grid and occupation contract. Import validates
finite density, the existing negative-roundoff allowance (1e-14 electrons/Bohr^3),
relative charge error <= 1e-8, orbital dimensions and orthonormality <= 1e-8,
history lengths/indices and mixer dimensions/settings. No renormalization occurs.
The original SCF code resumes at the next complete iteration. Wall-clock timings
are intentionally excluded from bitwise continuation claims.

`WarmStart` contains only the same shared `ElectronicState` and identities. A new
calculation ID is required. Initialization resets mixer history, numerical history,
counters, policy and certification. Mixing, tolerances, solver, budget and thread
count may differ. Direct reuse still requires identical physical model, grid,
basis dimensions and occupations. Cell/grid transfer and changed band-count
orbital projection are **not implemented**. A density-only warm state may omit
orbitals. Exact artifacts are never implicitly treated as warm artifacts.

`ResultExport` contains the shared density/orbitals/occupations plus evaluated
energy audit, charge, sampled gap, per-k eigenvalues/residuals, convergence flag,
and completed history. It is an inspection payload, not a reconstructible live
`ScfResult`: setup/total wall time, Ewald component breakdown and solver-internal
result bookkeeping are not exported. Unconverged results remain unconverged.

## Identity and runtime boundary

Compatibility has named rejections for schema, kind, cell, atom positions,
species, per-atom pseudopotential SHA-256, cutoff, k mesh, electron count, XC,
spin, environment, executable build, platform, threads, options and array state.
The **trusted host** supplies model/ execution identities from the inputs used to
construct the problem. The engine cross-checks its retained cell, k points,
cutoff, electron count and XC; it cannot independently authenticate original UPF
files or recover discarded species/input metadata. Supplying false host identity
is outside this API's trust contract. The CLI hashes the actual UPF bytes and
running executable (including dirty code and linked backend), and records host
`uname -a`, architecture/OS and requested threads. This is conservative same-host
identity, not a portable bitwise restart guarantee.

Runtime callers use `Runtime::checkpoint` / `Runtime::restore_exact`, not the
read-only inner session's checkpoint API. The envelope also preserves the complete
bounded policy audit, revision, launch authority configuration, autonomy counts,
cooldowns, last numerical-policy boundary, revocation flag and paused/initialized
state. The session-only importer rejects a runtime envelope to prevent losing
policy history. Runtime import is a trusted-host operation; saved approvals are
historical provenance, not new approval tokens.

Exact capture requires all telemetry drained and durably published by the host.
Undrained queues/loss counters, pending controls, in-flight steps, failed journals
and external audit sinks reject explicitly. Arbitrary sinks and live pending
approval deadlines are not serializable capabilities. The CLI's local commands
use Classic sessions; runtime envelopes require the authenticated runtime host
API. The reserved agent `ControlAction::Checkpoint` still returns
UnsupportedCapability: no agent-to-filesystem authority is introduced here.

## Verification

```sh
scripts/check.sh
cargo test --release -p dft-engine --test checkpoint -- --include-ignored --nocapture --test-threads=1
cargo test --release -p dft-cli --test checkpoint_io -- --include-ignored --nocapture --test-threads=1
```

The tests use synthetic nonlocal cubic Gamma, four k points, and skew two-species
nonlocal-plus-local fixtures with Linear/Pulay/KerkerPulay and requested threads
1/4. Memory tests cover all three execution modes and boundaries 0/3/8 (including
history eviction); file tests serialize at boundary 3 and rebuild the problem.
Comparisons include complete density/orbital bits, scientific histories and
certificates, excluding timings. Separate cases exercise different-mixer warm
starts, paused runtime policy/audit restoration, consumed autonomy, named
compatibility failures, corrupted files and an independent byte-layout check.
These are engine self-equivalence tests, not new QE results or performance
measurements. Stored evidence is in `tests/reference/checkpoint/`.
