# dft-engine

**Purpose:** SCF orchestration, occupations, energy audit and post-SCF derivatives.
No file I/O, serde, transport, network or logging in numerical loops.

**Before editing:** read `../../AGENTS.md`, `../../ownership.toml`, the ticket,
[ADR 0047](../../docs/decisions/0047-reference-driver-single-path.md), and relevant
physics/performance contracts. Stay inside ticket ownership. Shared contracts
and numerical conventions cannot change incidentally during a refactor.

**Verify:** `scripts/check.sh dft-engine`, then `scripts/check.sh`; report actual
exit codes. Tests over about 10 seconds in debug stay ignored with a documented
release command. Physics oracles remain independent and are never deleted merely
to remove an implementation dependency. Do not relax numerical gates.

## Production ownership and entry points

`refscf::api::{solve, stepper}` is the production facade, with Gamma, complex-k
and electronic spin/FD dispatch. `refstart` initializes the physical problem;
`model/` owns shared options, results, verification and spin contracts. These
modules are top-level, not children of the legacy `scf` facade. Use the reference
`Request`/`Profile` for Broyden beta 0.7/history 8 and threshold 1e-8 Ry = 5e-9 Ha.
Shared `ScfOptions` alone does not select a driver. Verification defaults on;
strict certification is a different oracle contract. `refscf::nscf` solves the
fixed-density band problem; it does not reuse legacy SCF orchestration.

The CLI has Reference or StrictOracle routing, strict precedence and the sole
G4D free-moment PBE fallback. Do not add hidden fallback after an engine error.
Fixed spin, fixed-moment FD, free-moment LDA-PW, NLCC/PBE, k meshes and spatial
symmetry use their qualified reference paths. Qualification limits are in
[reference design](../../docs/design/reference-scf.md) and ADR 0047.

Keep Hartree units in physical contracts; reference dr2/ethr and policy
`threshold_ry` are explicitly Ry. Preserve charge/moment conventions, owned-row
ordering, same-rank determinism and rank-count acceptance. No material-specific
thresholds or implicit renormalization. Every energy term is independently
checked and reported; do not double count core or verification energies.

## Frozen legacy oracle (R9c, R10, R11)

`src/oracle/` is gated by `strict-oracle`, enabled by default. The CLI explicitly
enables it. Root `scf`, `ScfProblem`, `ScfSession`, runtime/checkpoint compatibility
exports are oracle APIs, not the production facade. Legacy standard SCF branches
were removed in R10/R11, see git history; strict spin/FD and the MPI bridge remain.
The oracle's narrow fixed-default ladder and `cutoff_ladder="off"` compatibility
pin preserve its original start identity; they are not the removed tunable ladder.

Freeze default arithmetic. Any change needs an oracle-specific ticket, local
L-ORC and dev08 RG-ORC identity; shared changes also require reference identity.
R11b removes oracle-only switches/default-off arms while retaining the default
operations. **Integration boundary at R12 base `27d91ece`: R11b is a ticket,
not a merged implementation; oracle-only switches still exist.** Do not infer
retirement from the policy or retune these switches for an acceptance run.
Shared-kernel switches need reader analysis. R10d MPI trimming remains a separate
decision with a 1e-12 Ha oracle-versus-oracle prerequisite.

Reference/model/physics compile without the feature. The gate independently
checks and lints `dft-engine --no-default-features`, without workspace feature
unification. Oracle-dependent test harnesses declare
`required-features = ["strict-oracle"]`. Reference stepper/checkpoint and module
independence tests also compile without it. CLI `--no-default-features` does
not disable an explicitly enabled dependency feature.

## Stepper, runtime and checkpoint

`refscf::stepper::Stepper` owns the thread-affine synchronous continuation:
initialize once, step one logical map, finish only a completed delivery.
Ordinary/retry/verification reports can repeat an iteration number. `safe_point`
is the mutation/capture boundary, not every map return. Numerical failure is
terminal; exhausted delivery stays nonconverged. Never replay maps to implement
stepping or serialize compiler-generated future bytes. Keep the run configuration
scope alive on the owning thread. MPI ranks drive the same sequence collectively.

`refscf::runtime::Runtime` owns lifecycle and scalar audit. Policy allow-list:
beta `(0,1]`, history `1..=25`, tightening positive `threshold_ry`, verification.
Mixer changes clear history; density/orbitals and the ethr ratchet remain.
Legacy solver/schedule/divergence/budget controls reject. Pause/cancel waits for
a safe return; terminal delivery wins. The historical budget-recovery demo is
not supported reference recovery. Transport approval and persistence belong to
CLI/MCP, never to SCF. See [runtime contract](../../docs/design/reference-scf.md#runtime-control).

`refscf::checkpoint` captures explicit v2 numerical state only on request:
density planes/coefficients, Broyden cache/history, kept orbitals/HX, potential
epochs, ratchet and verification/history. Rebuild immutable workspace and install
saved state before the next map; capture runs zero maps. Same build/input/model,
backend/platform, threads, ranks and partition are required. v2 rejects warm,
terminal exact, coarse-ladder/external-seeded continuations and capture after
runtime policy mutation (audit not persisted). Retired failed-map slots must be
empty. v1 under `oracle/` is historical/oracle only; no automatic conversion.
Wire format and file publication belong to [CLI checkpoint I/O](../../docs/checkpoint-format.md).

Useful release matrices (synthetic, not new QE/MPI speed evidence):

```sh
cargo test --release -p dft-engine --test refscf_stepper -- --include-ignored --nocapture
cargo test --release -p dft-engine --test refscf_checkpoint -- --include-ignored --nocapture
```

## Derivative and ionic-state contracts

`delivered`/`api` expose deferred force/stress state and `IonicSeed`. Prepare
observables once, then compute derivatives from the delivered verified state;
verification never substitutes its map output for that state. No forces from
an uncertified/unverified or failed candidate. Preserve per-term ledgers,
ordered atoms, and unavailable values as unavailable (never zero).
Stress is tensile-positive Ha/Bohr^3; pressure is minus trace/3. Frozen-Miller
stress and its separate Pulay assessment do not prove cutoff convergence.
Reference stress scans remain available without the oracle feature.

Live ionic seeds carry compatible density/orbitals, invalidate changed-potential
caches and reset SCF history. Cell seeds retain frozen support and transform
with `IonicSeed::in_cell`; do not confuse a warm seed with exact checkpoint
continuation. CLI optimizer restart currently lacks the reference orbital seed.
Read [relaxation](../../docs/user/relaxation.md), [cell relaxation](../../docs/user/cell-relaxation.md),
[forces](../../docs/validation/forces.md), and [stress](../../docs/validation/stress.md).

## Retired implementation instructions

The old descriptions of tunable cutoff ladders, direct minimization, bounded
SCF advance, fixed/spectral maps, first-block/replay experiments, standard stall
and density-forcing branches, verification scheduling and reference rollback
arms were **removed in R10/R11, see git history**. They are not instructions for
new work. Retained reference RMM/ladder algorithms are separate modules.
Detailed retirement/recreation provenance is in
[switch-retirement.md](../../docs/performance/switch-retirement.md); historical
measurements remain in the ledger and stored evidence. Root AGENTS.md is unchanged.

## Measurement and evidence boundaries

Use `baseline-reference.json` for production identity and `baseline-oracle.json`
with `BAANDI_REFERENCE_STRICT=on` for oracle identity. Never use historical
`baseline-legacy.json` as a production baseline. Mac timing is observational;
dev08 acceptance belongs to the orchestrator and must obey same-model comparison.
No new performance claim follows from a move or deletion. Spin/FD results lack
scalar map-boundary clocks; report missing measurements honestly. Historical
P11 v1 runtime/checkpoint evidence does not establish current v2 capabilities.
