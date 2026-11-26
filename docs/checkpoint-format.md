# Checkpoint and restart: v2 reference, v1 legacy/oracle

The production reference driver uses **checkpoint v2**. The CLI owns serialization;
`dft-engine::refscf::checkpoint` owns explicit in-memory state. No serde, file I/O
or checkpoint allocation is added to ordinary SCF maps. v1 is retained for the
strict oracle and historical artifacts; it cannot resume a reference trajectory.
Routing and scope: [ADR 0047](decisions/0047-reference-driver-single-path.md).

## Commands and outputs

```sh
# Capture a nonterminal safe point after iteration 3 (including recovery/mixing).
dft-cli checkpoint input.toml boundary.cp --iteration 3 --kind exact
# Resume and finish; result.cp is a ResultExport container, NOT a results directory.
dft-cli restart input.toml boundary.cp result.cp
# Export a fresh calculation's evaluated result without a resumable boundary.
dft-cli checkpoint input.toml evaluated.cp --kind result
# Inspect either version without executing SCF.
dft-cli status boundary.cp
dft-cli explain result.cp
```

Each write also creates `<output>.results.json` with routing metadata, not a full
`run` result. The container's export summary has energy bits (Ha), convergence,
iteration/verification counts and density/orbital hashes. Neither the companion
nor `ResultExport` supplies normal density cubes, band arrays or a restartable
state. `run --output-dir DIR` is the normal scientific-output command.

Existing container and companion destinations are not overwritten. An exact
capture must be nonterminal and safe; result export preserves a nonconverged
outcome and restart then returns failure. v2 rejects `--kind warm` and `--warm`.
Use the same input bytes, UPFs, executable, backend, platform, threads, MPI ranks
and partition for exact restart; even an input comment change rejects.

## Reference checkpoint v2 (R8b)

`checkpoint` selects v2 when routing selects the reference driver. `restart`
selects its decoder from the checked container version: v1 always uses the
legacy session, while v2 requires the reference driver. There is no conversion
of a v1 continuation into a reference trajectory. There is no legacy standard
route; strict precedence and the G4D oracle exception are described in [ADR 0047](decisions/0047-reference-driver-single-path.md).

```sh
dft-cli checkpoint input.toml state.cp --iteration 3 --kind exact
dft-cli restart input.toml state.cp result.cp
dft-cli checkpoint input.toml uninterrupted.cp --kind result
```

The `BAANDICP` envelope described below, little-endian header/payload lengths, whole-file
SHA-256, 512 MiB total limit, 32 MiB header limit, and atomic no-clobber publication
apply. A **single container holds rank-indexed sections**. The header has
`schema_version: 2`, `kind`, `identity`, `iteration`, `states`, and `result`.
Each `states` entry is `{offset, count}`, both in **bytes**, relative to the
binary payload. Its array index is the rank. Sections are contiguous, nonempty,
and cover the entire payload. Unknown/duplicate header fields, unknown versions,
overflows, gaps, overlaps, trailing bytes, and bad checksums reject. The v1
reader remains separate and rejects v2.

The identity contains SHA-256 hashes of the actual input bytes, ordered UPF
bytes, and running executable, plus platform, thread count, and rank count.
The engine additionally matches the effective run configuration, options,
profile, cell/atoms, grid, local row layout, and immutable ionic data. Exact
restart requires the same rank count and partition, not merely the same global
density. Identity is supplied by a trusted host to the in-memory engine API;
the CLI computes it itself. Even an input comment change changes input identity.
A checksum detects corruption, not malicious rewriting or authenticity.

### Safe-point exchange and payload

`Stepper::checkpoint(identity, rank)` polls the suspended continuation solely
to visit its explicit state. It executes **zero SCF maps**. Capture is legal
after initialization (iteration 0) or a nonterminal report with `safe_point`
true. Candidate boundaries awaiting verification, first-map retry boundaries,
failed runs, and delivered terminal results reject. The CLI advances to the
safe point of the requested iteration, including its verification/recovery if
needed; it never publishes a terminal result as a resumable checkpoint.

`Stepper::restore` initializes a fresh stepper, rebuilding immutable basis,
FFT, operator, and solver workspaces, then installs the saved state before the
next map. Completed maps are never replayed. The suspended future is not copied
or serialized. Ordinary stepping creates no checkpoint arrays, serialization,
file I/O, or additional heap allocations. Encoding/cloning/gathering happens
only on an explicit checkpoint request. All ranks participate in lockstep;
rank 0 alone publishes the gathered file, and publication failure is agreed
before the command returns. Packing six bytes into an exactly representable
f64 integer permits byte transport through the existing collective interface;
there is no arithmetic reduction of payload values.

The binary schema is the typed traversal in `refscf/checkpoint.rs`, with the
field visits at the driver safe points and serialization helpers in
`broyden.rs`, `kept_hx.rs`, `k_kept.rs`, `map.rs`, and `k_map.rs`:

- Integers and boolean/option tags occupy unsigned little-endian 64-bit words;
  booleans/tags accept only 0 or 1. f64 values retain their raw IEEE-754 bits.
  Vectors have a u64 element count followed by their typed elements. Strings
  have a u64 byte length and UTF-8 bytes. Fixed arrays omit a length. Matrices
  carry rows and columns followed by column-major elements; complex values
  store real then imaginary components. Durations store seconds and nanoseconds.
- Both next-input density planes (electrons/Bohr^3) and reciprocal coefficients
  (electrons/Bohr^3) are stored in their owned order. Neither is reconstructed
  from the other by an FFT. Broyden stores previous input/residual, every
  difference-history slot in slot order, iteration/retention counters, and its
  incremental Gram cache. The cache is rebound to the rebuilt metric identity.
- Gamma ORC orbitals and complex spin/k orbitals, kept X/HX, owned potentials,
  potential epochs, and retained real-space pair slabs are saved exactly.
  Cache identities are rebound to the rebuilt fixed operators. The next map
  installs its complete potential before applying H; restoration establishes
  the old epoch without recomputing HX or changing the cache-hit decision.
- The ethr ratchet (Ry), previous dr2 (Ry), completed iteration, ordinary and
  verification histories, predictor state, first stopping iteration and cached
  fields are retained. Failed-verification reuse was removed in R11a; its reserved
  pending-map slot must be empty, and a nonempty retired slot rejects. Historical
  diagnostic NaN placeholders are preserved; live density, orbitals, mixer,
  and cache arrays are shape/finite checked. No normalization or clipping is
  performed during import. Start options/configuration and immutable model
  checks retain the deterministic atomic-start contract.

Wall-clock fields and allocator/cache high-water diagnostics are not part of
the numerical identity guarantee. Stored historical timings retain their
original process origin; new wall measurements start in the restored process.
Numerical comparison is against an uninterrupted run with the **same rank
count**, not between different rank counts.

### Result inspection and scope

A v2 `ResultExport` has no continuation sections. It contains evaluated energy
bits (Ha), convergence, iteration and verification-attempt counts, and SHA-256
hashes of delivered density and orbitals. This is an inspection summary, not a
warm start or a replacement for normal scientific output files. Hash order is
spin/channel, grid order for density, and spin/k/column-major real-imaginary
order for orbitals, using little-endian binary64 bytes. A nonconverged export
remains nonconverged, and `restart` returns a failure status after publishing it.

Exact v2 supports fixed-cutoff Gamma, complex k, and electronic spin/FD
steppers; prepared-symmetry report state is retained. A reference FD run with
automatic band count (`scf.nbands` omitted) additionally records the growth
mode, the current band count and the growth events after the policy record, and
a restart resizes the orbital blocks to that count before reading them; the
fixed-count state is unchanged. Capture after a runtime
policy mutation rejects because v2 does not persist its runtime audit. Gamma coarse-ladder or
external seeded continuations, external electronic trace sinks, warm imports,
and terminal exact captures currently reject explicitly. No cross-build,
cross-backend, changed-thread, changed-rank, or changed-model promise is made.
The existing v1 warm/result behavior remains available only through the
legacy/oracle format and its compatibility checks.

R8b evidence is recorded separately in
`tests/reference/checkpoint/reference-v2-results.json`; the historical v1
`results.json` is unchanged. Full identity matrices are intentionally ignored
in the fast debug gate; reproduce them with:

```sh
cargo test --release -p dft-engine --test refscf_checkpoint -- --include-ignored --nocapture
cargo test --release -p dft-cli --test checkpoint_v2 -- --include-ignored --nocapture
```

The rank 1/3/8 matrix uses the deterministic collective simulator. It establishes
same-rank continuation identity, not a dev08 MPI performance measurement or a
new Quantum ESPRESSO reference.

## Byte layout

All offsets are byte offsets. A file contains exactly:

| Offset | Size | Content |
| --- | --- | --- |
| 0 | 8 | ASCII `BAANDICP` |
| 8 | 8 | Unsigned little-endian JSON header length H |
| 16 | 8 | Unsigned little-endian binary section length B |
| 24 | H | UTF-8 JSON object, no trailing non-whitespace content |
| 24 + H | B | Version-specific binary payload: v1 binary64 arrays; v2 rank-indexed typed sections |
| 24 + H + B | 32 | Raw SHA-256 digest of all preceding bytes |

Readers require file length `56 + H + B` with checked arithmetic. The initial
adapter limits the complete file to 512 MiB and its header to 32 MiB. There are no
external payload paths. Truncation, trailing bytes, wrong magic, unsupported
schema, checksum errors, invalid live numerical arrays, array offset/length errors and unknown
or duplicate JSON fields reject. No migrations are implemented.

Publication writes a same-directory exclusive temporary file, fsyncs it, then
creates an atomic **no-clobber hard link** at the destination and fsyncs the
directory. Finally it removes the temporary name. A filesystem without hard-link
support fails explicitly. A crash can leave a temporary file but cannot publish
a partially written destination. A directory-sync error may leave a complete
published file while returning an error; callers should inspect it before retry.
A checksum detects corruption, not malicious alteration or source authenticity.

## Legacy/oracle v1 reference

The following is the historical v1 contract, not a production reference recipe.
Use `BAANDI_REFERENCE_STRICT=on` for oracle checkpoint creation/restart; v1 decoding
never converts state into v2. Old standard-profile artifacts may fail current
option/identity checks; readability does not imply restart compatibility.
The historical warm-start and runtime-audit capabilities below are not v2 features.

## Commands

```sh
export BAANDI_REFERENCE_STRICT=on
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

### Fixed-point verification (ADR 0043 amendment ONEPATH-4)

Version 1 carries the single standard path's fixed-point verification through
additive optional keys; the schema version stays 1. A key is omitted when it has
its default, so a file without verification is byte-identical to one written
before ONEPATH-4, and such files decode unchanged.

- `verification` (top level) and `continuation.options.verification`: `"off"`
  (default, omitted) or `"fixed_point"`. Exact restore requires both to match
  each other and the importing problem (`Options` rejection).
- `continuation.verification_attempts`: the attempts made before the boundary,
  oldest first (default empty, omitted). They are always failed attempts: a
  passed attempt ends the run, so an exact boundary never follows one. Each
  record has `iteration`, `passed`, `coulomb_error_ha`, `threshold_ha`,
  `eigen_tolerance_ha`, `residual_bound_ha` (the residual estimate, or null),
  `density_residual`, `occupied_eigenvalue_change_ha`, `strict_criteria_met`,
  `davidson_iterations`, `operator_applications`, `wall` (a duration, excluded
  from bitwise claims) and the optional spin/occupation fields of the engine's
  `VerificationAttempt`. Exact restore requires one attempt per history entry
  that holds a verification record (the standard history's
  `fixed_point_certificate`), in iteration order, not passed, with the session's
  threshold and finite nonnegative numbers (`History` rejection otherwise).
- `result.verification_attempts`: every attempt of an exported run (default
  empty, omitted); the last one passed when the export is converged on the
  verified standard path.

The restored session continues with the saved attempts, so its final
verification report equals the live session's continued after the checkpoint
(wall times excluded). As before, taking an exact checkpoint makes the live
session forget its kept REUSE_PSI state, so the restart reproduces the
checkpointed session; against a run that never checkpointed it is bit-identical
with `BAANDI_REUSE_PSI=off` and in the REUSE_PSI rounding class otherwise
(`crates/dft-cli/tests/all/onepath_last_fallbacks.rs`, bound 1e-10 Ha).

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
cargo test --release -p dft-engine --test all checkpoint:: -- --include-ignored --nocapture --test-threads=1
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
