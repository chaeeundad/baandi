# ADR 0047: Reference driver as the single production SCF path

- Status: **Accepted**, user decisions D1–D6, 2026-11-25.
- Scope: Stage R; [decision record](../roadmap/refactor-r-plan.md#7-user-decisions-2026-11-25).
- Supersedes: the one-step legacy removal rule in the reference design; ADR 0043's
  legacy implementation/routing details; ADR 0044's optional energy-verification
  experiment (retired in R11a). ADRs 0045 and 0046 retain their physics contracts; their historical rollback
  switch availability is subject to the D3 retirement inventory.

## Production routing

The CLI decision is `Reference | StrictOracle`. The default is `Reference`;
`scf.driver = "reference"` is an optional no-op. `BAANDI_REFERENCE_STRICT=on`
takes precedence and selects the frozen strict oracle, including its strict
certificate. The only automatic fallback is free-moment collinear Fermi–Dirac
PBE (G4D): it selects **StrictOracle**, not a legacy standard driver. The recorded
reason is `reference driver: free-moment collinear FD with PBE not yet qualified, G4D`
(the code constant is `G4D_FALLBACK`). This fallback does not establish convergence
or a qualified magnetic ground state. Remove it only through G4D qualification
with independent accuracy and convergence evidence, not on a timeout or engine error.
There is no retry with another driver after a numerical failure.

`scf.driver`, `driver_reason` and `driver_source` identify the actual selection
in results and inspection output. Source (`default`, `env`, `input`) identifies
the winning request, including when eligibility changes its route. Every oracle
gate must assert `strict_oracle`; a zero difference without route evidence is
not an oracle comparison. `refscf-map --solver legacy` remains an explicit
oracle diagnostic. `compat.mode = "reference"` is a physical-convention option
(ADR 0038), not a driver selector.

Reference means the reference Broyden/plain mixer and reference Davidson maps:
beta 0.7, history 8, standard threshold 5e-9 Ha and fixed-point verification on.
`pulay` is an input alias for Broyden on this route. The old standard mixer
beta 0.5/history 6 and strict beta 0.3 are not production defaults. Explicit
beta has the reference Broyden meaning. Stress/variable-cell inputs default to
5e-11 Ha unless a threshold is supplied. `verify = "off"` skips independent
verification; derivative requests require verification. Stop/verification
estimators are not mathematical bounds on total-energy error. The adopted
verification predictor can defer an attempt; successful verified delivery still
requires an actual independent map, and all attempts/recovery count toward wall time. Accuracy classes
remain strict 1e-10 Ha and standard 5e-9 Ha; the 1e-12 Ha cross-path check on
unique-solution fixtures is an internal regression check.

## Oracle isolation and freeze (D2)

The full oracle lives under `crates/dft-engine/src/oracle/`, behind the
`strict-oracle` Cargo feature, enabled by default in the engine and explicitly
by the CLI dependency. `model/` holds shared options/results/verification/spin
contracts; `refscf` and `refstart` are top-level production modules. Compatibility
exports do not make `ScfProblem` or `ScfSession` production entry points.

The reference engine builds with `cargo check -p dft-engine --no-default-features`;
`check.sh` also checks Clippy in that configuration. This is an engine isolation
gate, not a claim that `dft-cli --no-default-features` excludes the oracle.
Oracle-dependent test harnesses declare `required-features = ["strict-oracle"]`.
MPI oracle support is retained. Optional R10d trimming is a separate decision,
requiring oracle-versus-oracle agreement within 1e-12 Ha before removing it.

Freeze the oracle's default arithmetic. Changes require an oracle-specific
ticket, local identity and dev08 RG-ORC gates, plus reference identity for shared
code. Oracle-only tuning branches are to be removed by R11b; do not retune the
oracle to match a production result. **Checkout boundary:** R12 was authored at
`27d91ece`, which contains the R11b ticket but not its implementation. Oracle-only
switches still exist here; the freeze is the accepted rule, not evidence that
R11b removal has landed. The [retirement inventory](../performance/switch-retirement.md)
distinguishes completed R11a removals from that follow-up.

## Input and switch retirement (D3, D4, D5)

D4 removes legacy-only keys outright: ordinary unknown-key errors, including
when strict mode is requested. The removed `[scf]` keys are `mixing_q0_inverse_bohr`,
`mixing_metric`, `hartree_mixing`, `spectral_guard`, `locking`, `preconditioner`,
`davidson_capacity_per_band`, `energy_tolerance_ha`, `density_residual_tolerance`,
`eigen_tolerance`, `standard_stall_guard`, `standard_density_forcing`, and
`cutoff_ladder_fraction`, along with removed tolerance unit aliases.
`xc.density_floor` is removed. `mixing = "linear" | "kerker-pulay"`,
`driver = "legacy"`, `start_density = "uniform"` and non-atomic
`start_wavefunctions` are unknown **values** of retained keys. Atomic start is
the production default; missing atomic charge/orbital data use neutral uniform
density/deterministic kinetic completion, not a legacy route. `convergence_profile`
retains its earlier dedicated removal error from ADR 0043.

The retained `scf.cutoff_ladder = "off"` is an **oracle compatibility pin**:
omission preserves its historical fixed-default ladder; other values reject.
It does not enable the reference ladder, which is a separate opt-in algorithm.
This exception preserves the oracle identity recovered after gen20.

D3 retires adopted reference rollback switches and their pre-adoption branches,
plus rejected default-off experiments. There is no two-generation wait; recreate
an experiment from git when needed. This supersedes Stage R's standing requirement
to keep every adopted default switchable. R11a removed 29 rollback switches and
six experiments. `REFERENCE_SCF` was removed in R10c; do not export it in new
commands. Registry entries, not historical recipes, define accepted overrides.
No general retired-variable warning/migration table is promised.

Keep diagnostics/profiling, infrastructure, `REFERENCE_STRICT`, protected
isolated algorithms (RMM, reference ladder, preconditioning), backend selection,
and explicitly retained options documented in the inventory. Shared/oracle
switches need reader analysis and oracle identity before removal; D3 is not
permission to change shared arithmetic. Energy-only verification and failed-map
reuse experiments are removed; compatibility payload slots do not reactivate them.
D5 keeps the ledger-backed relaxation extrapolation as an opt-in ionic seed
operation; `RELAX_EXTRAPOLATION` and `relax.start_extrapolation` remain supported.

## Runtime and persistence (D1)

Adapters call `refscf::api::{solve, stepper}`. `Stepper` owns a synchronous,
thread-affine continuation. Initialize once, step by logical map, then retain or
finish the terminal delivery. Ordinary, retry and verification maps have distinct
reports; iteration numbers may repeat. Interventions wait for `safe_point`;
stepping never replays completed maps or serializes the suspended future.
`refscf::runtime::Runtime` adds pause/resume/cancel, scalar snapshots and audit.

Reference policies are mixing beta `(0,1]`, history `1..=25`, tightening
`threshold_ry`, and verification. Mixer changes clear history; density/orbitals
and the eigenvalue-tolerance ratchet remain. Legacy solver/tolerance/schedule or
budget-extension controls reject. Observe/agent CLI modes, Python stepping and
the host-embedded MCP adapter use this runtime; transport approval remains in
adapters, never in numerical loops. The historical budget-recovery demo is not
reference recovery evidence: the reference demo records an unsupported-policy
rejection. Checkpoint after policy mutation is currently unsupported because v2
does not preserve the runtime audit. See [runtime contract](../design/reference-scf.md#runtime-control).

[Checkpoint v2](../checkpoint-format.md) is the primary reference format:
rank-indexed state in a checksum-protected BAANDICP container, including density
planes and coefficients, Broyden history/cache, kept orbitals/HX, potential epochs,
ratchet, verification and ordinary histories. Same-build/model/input/UPF/backend,
threads, ranks and partition are required for exact continuation. Capture executes
zero maps at a nonterminal safe point; restore rebuilds immutable workspace and
installs explicit state before the next map. No cross-rank bitwise promise follows.

v1 remains legacy/oracle only; there is no conversion into a reference trajectory.
v2 supports fixed-cutoff Gamma, complex-k and electronic spin/FD steppers; warm
imports, coarse-ladder/external-seeded continuation and terminal exact capture
reject. `restart input.toml state.cp result.cp` writes a **ResultExport container**
and routing companion, not a normal results directory or another continuation.

## Evidence: gen14–gen28

The following are relayed dev08 results, not runs performed by R12. Paths name
remote evidence on dev08; the RESULT files are not bundled in this checkout.
No new QE reference values are created by this ADR. Cross-code offsets cite
existing QE runs; speed comparisons must use the same parallel model and matched
conditions. Single wall samples are not a new speed-adoption campaign.

| Generation | Recorded result and limits | RESULT path |
| --- | --- | --- |
| gen14 (`d81ef8af`) | Eight RG-ACC cases VERIFIED within 5e-9 Ha of strict; Si8 force difference 8.2e-6 Ha/Bohr, stress 1.1e-9 Ha/Bohr^3. Five of six README cases unchanged within 4e-10 Ha/atom; PBE Gamma negative-density failure subsequently fixed by screening. | `/data/sig-dft/gen14/RESULT.md` |
| gen15 | Fixed-moment bcc Fe PBE+NLCC, MPI8: VERIFIED in 75.8 s. Legacy and strict timed out at 2700 s; no oracle agreement claim. | `/data/sig-dft/gen15/RESULT.md` |
| gen16–gen18 | No summaries or RESULT files supplied in the available section 8; no pass inferred. | `/data/sig-dft/gen16/RESULT.md`, `/data/sig-dft/gen17/RESULT.md`, `/data/sig-dft/gen18/RESULT.md` (not inspected) |
| gen19 (`22531ab1`) | R8a Si128 MPI8 energy/density bits, maps, iterations and H columns identical to `5d787b70`, verify off/on; single wall ratios 1.021 / 1.005. Bands/DOS both stages reference. PBE Gamma VERIFIED, archived QE offset +1.6983e-6 Ha/atom. O triplet VERIFIED in 9 iterations (oracle not converged in 200). Cell-relax lattice difference 6.3e-8 Bohr, pressure 2.3e-9 Ha/Bohr^3. Si32 map-5 checkpoint/restart bit-identical. | `/data/sig-dft/gen19/RESULT.md` |
| gen20–gen21 | R10A oracle regression bisected at `/data/sig-dft/bisect20`, fixed in `ab77ec10`: five RG-ORC cases again bit-identical to gen19 (MPI8/np1, forces/stress). Observe/plain identical. Fe moment-FD default identical to gen15; O fixed-moment FD PBE+NLCC VERIFIED after G4F, oracle growth-limit failure retained. | `/data/sig-dft/gen20/RESULT.md`, `/data/sig-dft/gen21/RESULT.md` |
| gen22–gen26 | Not present in the supplied plan's evidence log; no validation result inferred. | `/data/sig-dft/gen22/RESULT.md` through `/data/sig-dft/gen26/RESULT.md` (not inspected) |
| gen27 | R12 ticket reports preflight/result start-fallback reason mismatch; no numerical conclusion supplied. | `/data/sig-dft/gen27/RESULT.md` (not inspected) |
| gen28 | RESULT contents unavailable in this checkout; no pass inferred. | `/data/sig-dft/gen28/RESULT.md` (not inspected) |

Local reference-specific evidence remains separate from the legacy README data:
[checkpoint v2](../../tests/reference/checkpoint/reference-v2-results.json),
[runtime](../../tests/reference/runtime/generated/reference-results.json),
[Python](../../crates/baandi-py/tests/r5-reference-evidence.json), and
[MCP](../../crates/baandi-mcp/validation.json). These do not replace remote
performance or physical accuracy gates.

## Open items and documentation (D6)

- G4D free-moment PBE remains unqualified; the oracle fallback can itself fail.
- Restart currently exports a container, not normal scientific result files.
- Spin/FD results lack the scalar driver's map-boundary clock fields; missing
  clocks are unavailable measurements, not zero costs.
- dev08 gen27 reports different start-fallback explanations in preflight and
  results. Do not treat the preliminary reason as proof of the executed start.
- R11b implementation and gen22–gen28 RESULT contents must be reconciled by the
  integrating session; this checkout cannot establish their completion.

D6 requires reporting README corrections before editing user-maintained prose.
The [correction list](../roadmap/readme-corrections.md) is the approval artifact;
README.md is unchanged. Its generated block remains a historical evidence view,
with its original drivers and limitations; freshness is not current-path coverage.
