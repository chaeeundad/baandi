# EXP-TABPFN-SCF: external mixing controller experiment

Experiment branch only. No merge or default adoption without a separate decision.
Mac, one rank, small real-Gamma reference-driver cells. This implements the Rust
boundary and a standard-library Python **test fixture**, not TabPFN or training.
The orchestrator owns the learned Python executable and subsequent experiments.
No QE source or new numerical dependencies were used.

## Controls and numerical contract

All three environment switches are registered in the run manifest, environment
only, and disabled by default:

| Switch | Value | Effect |
| --- | --- | --- |
| `EXP_SCF_TRACE` | `on` / `off` (default) | Write `<output-dir>/scf_trace.jsonl` |
| `EXP_SCF_CONTROLLER` | Executable path; unset/empty disables | Start one child before SCF; implicitly enable trace |
| `EXP_SCF_BETA_EXPLORE` | Decimal unsigned 64-bit seed; unset disables | Deterministic beta draws; implicitly enable trace |

Controller and exploration together are an input error. Experiment runs require
`scf.driver='reference'` and `scf.verify='fixed_point'` (the default). The serial
engine hook also refuses verification-off and multiple ranks. No MPI protocol,
legacy SCF controller, checkpoint continuation, or relaxed stopping test is added.
Use fixed-rule inputs `mixing='broyden'`, `mixing_beta=0.7`, `mixing_history=8`.
Fallback restores the input profile's fixed beta on each subsequent update and
keeps its history. Ordinary non-experimental profile validation is unchanged.

The only numerical intervention is the beta of the current Broyden update and
an optional logical history clear, including the previous input/residual pair.
Valid action beta is finite and in `(0, 1.5]`; default profile input beta still
has its original range. Invalid actions are rejected before either mutation.
The ethr ratchet, first-step redo, estimator, verification residual checks,
`REFSCF_VERIFY_GATE`, first gate and verification reuse remain unchanged.
Verification is evaluated before asking for an action. A terminal map is still
observed, but its returned action is recorded as unapplied. Failure to converge
or verify retains the existing unsuccessful result/exit behavior. Trace data
are never a convergence certificate.

With verification reuse enabled, a failed check may supply the next map and its
input density supersedes the mixing proposal, as in the original driver.
`reused_from_iteration` identifies the reused map; work is charged once at the
original verification, not again at the reused iteration. Disable reuse in all
arms together if measuring a controller whose proposed updates must always be
used. Do not silently change the gate or reuse setting in just one arm.

## Protocol and causal ordering

An executable path is passed directly to `Command`, without shell expansion or
argument parsing. To run Python, make a wrapper executable containing
`exec python3 -I /absolute/path/controller.py` (quote paths as appropriate).
The child receives one newline-terminated UTF-8 JSON observation per outer SCF
iteration on stdin and must flush exactly one reply on stdout:

```json
{"beta": 0.7, "reset_history": false}
```

Both fields are required, with JSON number and boolean types respectively;
unknown/duplicate fields, nonfinite values, out-of-range beta, malformed UTF-8,
EOF, unterminated replies and replies larger than 64 KiB are failures. Stdout
is exclusively the protocol; use stderr for logging. There is a 30-second
exchange deadline, including a blocked stdin write. The first failure kills and
reaps the child, disables further requests, and records the reason on every
subsequent fixed-rule step. Spawn failure also falls back. The child is killed
and reaped at the end of the run, including error unwinding. The executable must
manage any descendants it creates; this adapter manages the direct child only.

There is an unavoidable causal distinction: the action **for this update**,
its Broyden diagnostics, mixed-density `descf`, and controller elapsed time do
not exist when requesting that action. Therefore request and completed trace
use the same `exp-scf-trace-1` schema with `phase='observation'` or `'complete'`.
Unavailable fields are null, never invented future values. The completed file
line embeds the exact observation JSON object under `observation`; it then adds
the returned/applied action and measured outcome. It is one file line per outer
map. A byte-identical *completed* line cannot be sent before choosing the action
that creates it. Training must use `observation` as features to avoid leaking
action-dependent future values. Prior completed outcomes can be joined offline
by iteration; the child can retain its own earlier observations and replies.

Flow: compute map and estimator -> unchanged verification/gate -> CLI observation
and reply -> Broyden update -> mixed-energy correction -> CLI completed trace and
flush -> next iteration or delivery. The engine owns only typed in-memory data
and callbacks; serialization, child process handling and file I/O live in
`dft-cli/src/run.rs`, outside eigensolver, FFT and mixing kernels.

## Trace schema and units

One outer map means one SCF iteration; first-step redo and bounded Davidson
reentries are aggregated in its work counts, and verification work is listed
separately. The trace is streamed at every diagnostic level. To obtain the last
Ritz correction, experimental solves opt into existing solver detail collection
even when normal result diagnostics are lightweight. This adds measurement
cost but uses the same arithmetic, tested against the non-experimental run.
This is an instrumented experiment, not a production speed claim.

| Keys | Meaning |
| --- | --- |
| `schema`, `phase`, `iteration` | Version, causal stage, one-based outer SCF map |
| `dr2_ry`, `ethr_ry`, `ethr_requested_ry` | Hartree-metric squared density residual in Ry; used/requested Ritz threshold in Ry |
| `first_iteration_redo` | Whether the first outer map needed its tighter redo |
| `beta_used`, `beta_selected`, `reset_history`, `action_applied` | Dimensionless beta used by mixing (null on terminal maps), reply/draw, history clear, whether update was attempted |
| `history_occupancy` | Request: history before this update; complete: secants used in this update |
| `broyden_gamma_l2`, `broyden_gamma_max`, `broyden_condition` | Dimensionless norms/condition from existing Diagnostics; condition null without a history solve; gamma fields null in observation |
| `davidson_iterations`, `h_columns` | Actual ordinary work, including redo/reentries, excluding separately charged verification |
| `ritz_last_max_change_ha`, `ritz_last_mean_change_ha` | Absolute Ritz changes across all bands at the final correction of the last ordinary solve, Ha; null if no new ordinary solve (reuse) |
| `band_energy_ha`, `hartree_energy_ha`, `xc_energy_ha`, `total_energy_ha` | Output-density audit terms; total is `EnergyAudit::total_hartree()`, Ha |
| `delta_band_energy_ha`, `delta_hartree_energy_ha`, `delta_xc_energy_ha`, `delta_total_energy_ha` | Current minus previous outer map's corresponding audit term, Ha; null on first map |
| `double_counting_ha`, `ewald_energy_ha` | Audit double counting and ion-ion term, Ha |
| `variational_energy_ha`, `delta_variational_energy_ha`, `variational_hartree_ha`, `variational_xc_ha`, `descf_ha` | Reference mixed-density energy/correction and its change, Ha; null until update complete; terminal map has zero descf |
| `estimator_stop` | Unchanged `dr2_ry < threshold_ry` decision, not verification |
| `verification_status`, `verification_error_ha` | `not_attempted`, `skipped`, `failed`, `passed`; actual fixed-point discrepancy in Ha when measured |
| `verification_h_columns`, `verification_davidson_iterations`, `verification_wall_s` | All verification work and inclusive elapsed seconds charged here |
| `reused_from_iteration` | Earlier iteration whose failed verification provided this map, else null |
| `map_wall_s` | Boundary-to-boundary seconds: ordinary map, verification, controller exchange, mixing and energy correction; excludes final trace serialization/flush and final state delivery |
| `controller_wall_s`, `controller_startup_s`, `controller_status` | Exchange/validation elapsed seconds; process/thread startup on first row only; `off`, `exploration`, `accepted`, or `fallback: <reason>` |
| `exploration_seed`, `exploration_draw_u64` | Exact decimal strings (avoid JSON consumer integer precision loss), null outside exploration |
| `mixing_error` | Existing Broyden failure text, null on success |
| `observation` | Exact object supplied to `choose`; present only in completed file lines |

Null means unavailable, not zero. Gamma norms on a terminal completed map are
zero because no update was performed; `action_applied=false` distinguishes it.
A fatal solver/FFT/verification numerical exception returns the ordinary engine
error; already completed rows remain flushed, but no completed row is fabricated
for a map that did not reach the callback. A trace-write failure fails the run.
Keep the CLI exit code and stderr with every trajectory, including failed ones.

## Exploration and planned evaluation

The generator is SplitMix64 with wrapping unsigned arithmetic. Each observed
outer map, including the terminal observation, consumes one draw `z`:
`beta = 0.2 + 0.8 * ((z >> 11) / (2^53 - 1))`, with the division in binary64.
This is uniform on a finite 53-bit grid covering `[0.2,1.0]`. The seed and full
64-bit draw are recorded. Seed 42's first draw is 13679457532755275413. A reset
does not rewind the generator; no global RNG or dependency is introduced.

Example commands after a debug build, using a reference-driver input:

```sh
EXP_SCF_TRACE=on target/debug/dft-cli run small.toml --output-dir runs/fixed
EXP_SCF_BETA_EXPLORE=42 target/debug/dft-cli run small.toml --output-dir runs/explore-42
EXP_SCF_CONTROLLER=/absolute/path/controller target/debug/dft-cli run small.toml --output-dir runs/controlled
```

For each of Si8 displaced, Si16, rskew, elongated 3x1x1, and the synthetic
nonlocal fixture, retain an explicit input and pseudopotential SHA-256, cutoff,
FFT grid, band count, start, threshold and all switches. Do not commit UPFs.
Generate training trajectories with seeds 0 through 31; reserve seeds 1000
through 1007 and disjoint geometries for evaluation. Include failures, timeouts,
resets and verification rejections; do not train only on successful trajectories.
Keep validation/test geometries out of fitting and controller selection.

Run a fixed beta 0.7/history 8 baseline, exploration, and the frozen learned
controller on identical input conditions. Also obtain a separate strict reference
at that geometry for `E - E_strict` (Ha); never confuse the ordinary reference
driver name with the strict accuracy class. Use an external monotonic process
clock to charge controller startup, all maps and failed verification attempts,
serialization, trace flush, teardown and final delivery. Report both total wall
and total wall minus `sum(controller_wall_s) + sum(controller_startup_s)`;
the latter still contains instrumentation and must be labelled accordingly.

Metrics: outer maps (plus redo count), sum ordinary and verification H columns,
sum ordinary and verification Davidson iterations, full delivered wall, controller
cost, `E - E_strict` in Ha, maximum absolute density difference in electrons/Bohr^3
and `sqrt(dV * sum((rho-rho_default)^2))` in electrons/Bohr^(3/2) against the
fixed-rule run's **delivered** density, verification attempts/skips/failures and
final pass status. Compare matched one-rank runs; no QE speed ratio is claimed.
An improvement in map count alone is not an improvement in delivered wall or
physics. No generality or performance adoption follows from the synthetic tests.

## Local validation

The ticket's debug commands are the verification scope; no full `check.sh`,
release matrix, remote execution or ML installation is authorized here. Tests
exercise the real CLI, Python isolated mode, generated temporary pseudopotentials,
bit-identical density bytes and map/verification scientific records, both normal
diagnostic levels, trace-only mode, garbage/EOF/range/missing-field/spawn fallback,
seed repeatability, history reset, beta above one, and verification-off rejection.
Observed debug synthetic result (Mac; no material or QE claim): fixed and echo
both used 5 outer maps, 37 ordinary H columns, 7 verification H columns and
10 ordinary Davidson corrections. Delivered energy was
-1.1050370150385769 Ha; the verification error was
7.287844397843669e-12 Ha. Energy difference was exactly 0 Ha and delivered
`density.bin` bytes matched at both normal diagnostic levels. Exploration seeds
42 and 43 produced different first draws, and repeated seed 42 reproduced the
full scientific trajectory and density bytes.

Executed gates: `cargo check --workspace --tests`; `cargo test -p dft-engine
--test refscf` (19 passed, one pre-existing ignored external-data test);
`cargo test -p dft-cli --test exp_scf_controller -- --nocapture` (9 passed);
`cargo test -p dft-numerics switches` (14 passed across two harnesses);
`cargo fmt --check` and scoped all-target clippy with `-D warnings`. All exited 0.
The controller suite took 1.17 s total in the observed debug test run; this is
only a fast-test check, not an SCF performance result. Budget exhaustion also
kept the partial trace and reported `verified=false`.

Shared dashboard/session registry, README generator and performance ledgers are
outside this ticket's `OWNED_FILES`; they are left to the orchestrator. The
experiment remains uncommitted and has no automatic main-branch adoption.
