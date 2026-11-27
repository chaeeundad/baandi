# dft-cli

**Purpose:** CLI, TOML input, output files. No numerical logic.

**Invariants to preserve and test:** run/validate/inspect-upf subcommands; nonzero exit with cause on failure; outputs carry reproducibility metadata.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-cli`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.

## Current routing (Stage R)

[ADR 0047](../../docs/decisions/0047-reference-driver-single-path.md) is the current
contract. `driver_select` resolves Reference or StrictOracle before construction;
STRICT wins, and only G4D free-moment PBE falls back to the oracle. Assert the
recorded driver in gates; never retry an engine error on another path.
Reference uses `refscf::api`, beta 0.7/history 8, standard verification on.
Removed D4 keys are ordinary unknown-key errors even in oracle mode. Retained
atomic-only start keys and the oracle `cutoff_ladder="off"` pin are documented
in [CLI usage](../../docs/cli.md). `REFERENCE_SCF`/legacy driver are removed.

## Collinear input/results (P13B / MF03)

`SpinConfig` and shared preflight validate explicit integer sectors. Reference
`run::reference_electronic_run` serializes version 2 `reference_spin_fd_v1`;
`spin_run` is the strict-oracle serializer. Preserve the extension,
full spin/k/band axes, fixed-sector certificate and conventional density cubes.
Fixed-occupation unpolarized results retain their scalar schema; reference FD
uses the electronic version 2 profile even without spin. Never route a collinear result through a
scalar legacy orbital view or an unpolarized-only postprocessor. Reference
spin/FD runtime and checkpoint v2 are supported; v1 is not spin restart.
Use `tests/spin.rs` and its independent standard-library schema/cube checker.
See `docs/user/collinear-spin.md` for supported combinations and format losses.

P13C: `inspect-upf --json --manifest PATH` matches provenance by SHA-256 and
reports header facts even when SCF support fails. The independent stdlib reader
is `scripts/check_upf_headers.py`; external files are selected by `DFT_NC_DIR`.
`tests/all/inspect_upf.rs::real_nc_headers` is ignored and checks inspection only,
not SCF or material qualification. Engine tests cover NLCC and real atom smoke;
see docs/cli.md. NLCC local-potential export includes the core in XC only.
NLCC ledger data cannot be silently serialized into checkpoint v1.

## MF04 finite-temperature boundary

`finite_run.rs` writes the strict-oracle FD profile; production uses the
reference electronic serializer and its explicit E/TS/F and spin/k/band axes. Unit conversion is
only in config. Unknown schemes, competing widths, finite-T v1 checkpoints and
unsupported output/runtime options reject. Model identity excludes spontaneous
moment and seeds; state identity includes the actual occupations, mu and entropy.
See `docs/validation/occupations.md` for caps and reproduction commands.

## MF11 crystal symmetry (P19A)

`[symmetry]` defaults off. Auto selects a mesh subgroup and compatible FFT dimensions, records multiplicities/operations/orbits and model identity, and rejects unsupported model terms. Preserve legacy output and arithmetic when off.

See `docs/validation/symmetry.md` and `tests/reference/symmetry/` for commands,
observed failures and validation limits. No QE symmetry-on oracle is available.

## Force output (P20A)

`output.forces` defaults false. Since ADR 0043 (2026-11-16) the default standard
path supplies forces from a passed fixed-point verification (`verify = "off"`
with forces rejects at preflight); the strict reference path
(`BAANDI_REFERENCE_STRICT=on`) keeps the certificate preconditions.
All paths use the engine's post-SCF derivative implementation. Force-on output
has a per-term block and text table. The results proposal accepts extension versions 1 (historical unavailable) and 2
(available forces with strict certificate preconditions); keep outer schema_version
unchanged. `tests/validate_force_schema.py` checks fixtures and actual CLI output. Run `cargo test -p dft-cli --test forces`.

## Fixed-cell ionic driver (P21A)

`relax` is the ticket-authorized driver-layer optimizer exception to the general
I/O-only rule. BFGS/FIRE and optimizer restart/output ownership stay here, never
in the electronic numerical core. Keep physical forces separate from constraint
projection. Since ONEPATH-1 (ions) and ONEPATH-2 (cell, ions_cell) every mode
runs the standard path with the fixed-point verification on every step
(`verify = "off"` rejects); spin/FD modes reject. Retry a failed
candidate exactly once fresh at half displacement, then retain the last good
structure; never query uncertified forces. The restart uses float bits and a
checksum, includes the next warm density and optimizer memory, and is not an
exact electronic checkpoint. Commands, limits and measured evidence are in
`docs/user/relaxation.md` and `docs/validation/relaxation.md`. Long CLI trajectory
oracles are ignored in debug and must run with `--release`.

## Variable-cell driver (P23A)

`relax/cell.rs` owns the fixed reference fractional/strain coordinate chart,
constraint/symmetry projection, enthalpy backtracking, segment rebuilds and cell
optimizer persistence. Preserve the default ionic branch. Every frame retains
its own cell and raw derivatives; final nominal-cutoff rebuilds are separately
labelled. No electronic exact-continuation claim follows from optimizer resume.
Read docs/user/cell-relaxation.md for conjugate gradients and rejected constraints.
Run the ignored `cell_relaxation` CLI tests in release; the external-file oracle
requires DFT_NC_DIR and reports Si/MgO cutoff dependence without reference claims.

The P23A rebuild follow-up separates reporting from acceptance:
`real_files_cutoff_report` retains frozen-only start/cutoff differences without
asserting equal minima for different initial Miller sets. The ignored
`real_files_rebuild_start_agreement` invokes the stdlib `rebuild_protocol.py`
runner for N=3/5 and 2–3 fresh-reference polish cycles. Its acceptance uses final
rebuilt pressure and the unchanged previously reported lattice tolerance bounds.
Polish is a new optimizer calculation, never a renamed exact restart. Preserve
all stage inputs, times and failures. `rebuilt_pressure_check` warns independently
of the Pulay tensor-difference warning; frozen-only convergence must not imply
nominal-cutoff equilibrium, and the default rebuild cadence remains zero.

## Reference runtime and checkpoint

Observe/agent runs drive `refscf::runtime::Runtime` and the production stepper.
Policies are beta/history, tightening threshold in Ry and verification; legacy
budget/solver/tolerance controls reject. Keep JSON, sockets and audit outside
maps. Approval intent precedes mutation, disposition follows; do not promise
rollback after mutation if journal publication fails. MPI rank zero owns intake;
all ranks apply controls and agree transport errors at safe points.

Checkpoint defaults to v2; v1 is legacy/oracle only. Restart dispatches by checked
version and validates identity. v2 restart writes an export container and a
routing-only `.results.json` companion, not a normal results directory. No warm
v2 import or conversion from v1. Safe-point state and kept HX must survive exactly;
never rerun maps to fake continuation. Wire/publication rules and reproduction:
[checkpoint format](../../docs/checkpoint-format.md).

## Retired measurement and standard-path instructions

Legacy standard routes, tunable ladder JSON, FIRST_BLOCK/replay experiments and
reference rollback arms were **removed in R10/R11, see git history**. Old
`ScfSession` runtime/standard-profile descriptions are oracle history, not
production instructions. `PHASE_PROFILE` diagnostics remain. R11b oracle-only
switch retirement is pending integration at R12 base `27d91ece`; consult the
[retirement inventory](../../docs/performance/switch-retirement.md) before
claiming a switch disappeared. Do not edit root AGENTS.md.

## Run configuration and run report (RUN-CONFIG / RUN-REPORT, 2026-11-18)

`[performance]` and `[profile]` input tables set switches by registry name
(`run_config.rs`). `main` installs the input's configuration as the process
configuration before the first switch read (`MPI_TREE`) and, with two or more
ranks, checks `agreement_hash` across ranks (one allgatherv; mismatch = the
same error and exit 1 on every rank; rank 0's run ID adopted). `execute_mode`
uses the process object or a thread scope (`run_config::for_run`). results.json:
`scf.switches.<name>` gains `requested`, `consulted` (frozen after the SCF, so
report blocks do not count), `path`, `fallback_reason` (from counters,
`note_switch_paths`) and source `input`; the `run` block holds schema
`run-report-1`, run ID, exit reason, binary path and SHA-256, agreement hash,
run-local counters and the completion marker `completed: true` (written last).
`phase_profile` is `phase-profile-2` with `per_rank` and `closure`.
`scripts/parse_results.py` validates a results file (rejects truncated, stale
binary, missing/zero clocks of an active path, impossible rows); unit tests
`scripts/tests/test_parse_results.py`, real-output checks in `tests/all/run_config.rs`.

## Structure input (SI1, 2026-11-27)

`config/structure_input.rs` rewrites `[structure]` (compact atoms block or
`file` = POSCAR/extended XYZ) into the explicit `[system]` TOML value before
deserialization, maps nonzero per-atom `magmom` to `spin.initial_moment_seed`
(sum) and keeps per-atom `fixed` flags in `Config::structure_fixed` (merged
into relax settings only when `relax.fixed` is empty; both given = error).
`[pseudopotentials] dir` fills missing species entries by documented rules.
Parse inputs with a path through `Config::parse_at`/`parse_unvalidated_at` so
relative files resolve against the input directory. `[system]` inputs keep
their arithmetic; results gain `structure` (and `pseudopotential_lookup`).
Readers are checked against ASE 3.26 readings in
`tests/fixtures/structure_input/` (`make_fixtures.py`, test-only);
`cargo test -p dft-cli --test all structure_input`.
