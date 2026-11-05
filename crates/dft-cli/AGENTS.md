# dft-cli

**Purpose:** CLI, TOML input, output files. No numerical logic.

**Invariants to preserve and test:** run/validate/inspect-upf subcommands; nonzero exit with cause on failure; outputs carry reproducibility metadata.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies) and the ADRs in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-cli`, then `scripts/check.sh` for the whole workspace.
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.

## Collinear input/results (P13B / MF03)

`SpinConfig` and shared preflight validate explicit integer sectors; `spin_run`
serializes the additive version 2 result profile, proposal-schema extension,
full spin/k/band axes, fixed-sector certificate and conventional density cubes.
Unpolarized results remain version 1. Never route a collinear result through a
scalar legacy orbital view or unpolarized runtime/checkpoint/postprocessing.
Use `tests/spin.rs` and its independent standard-library schema/cube checker.
See `docs/user/collinear-spin.md` for supported combinations and format losses.

P13C: `inspect-upf --json --manifest PATH` matches provenance by SHA-256 and
reports header facts even when SCF support fails. The independent stdlib reader
is `scripts/check_upf_headers.py`; external files are selected by `DFT_NC_DIR`.
`tests/inspect_upf.rs::real_nc_headers` is ignored and checks inspection only,
not SCF or material qualification. Engine tests cover NLCC and real atom smoke;
see docs/cli.md. NLCC local-potential export includes the core in XC only.
NLCC ledger data cannot be silently serialized into checkpoint v1.

## MF04 finite-temperature boundary

`finite_run.rs` writes the explicit FD profile and proposed state-contracts
extension; default/fixed outputs retain their existing path. Unit conversion is
only in config. Unknown schemes, competing widths, finite-T v1 checkpoints and
unsupported output/runtime options reject. Model identity excludes spontaneous
moment and seeds; state identity includes the actual occupations, mu and entropy.
See `docs/validation/occupations.md` for caps and reproduction commands.

## MF11 crystal symmetry (P19A)

`[symmetry]` defaults off. Auto selects a mesh subgroup and compatible FFT dimensions, records multiplicities/operations/orbits and model identity, and rejects unsupported model terms. Preserve legacy output and arithmetic when off.

See `docs/validation/symmetry.md` and `tests/reference/symmetry/` for commands,
observed failures and validation limits. No QE symmetry-on oracle is available.

## Force output (P20A)

`output.forces` defaults false. Standard profile rejects at shared preflight;
all paths use the engine's post-SCF derivative implementation. Force-on output
has a per-term block and text table. The results proposal accepts extension versions 1 (historical unavailable) and 2
(available forces with strict certificate preconditions); keep outer schema_version
unchanged. `tests/validate_force_schema.py` checks fixtures and actual CLI output. Run `cargo test -p dft-cli --test forces`.

## Fixed-cell ionic driver (P21A)

`relax` is the ticket-authorized driver-layer optimizer exception to the general
I/O-only rule. BFGS/FIRE and optimizer restart/output ownership stay here, never
in the electronic numerical core. Keep physical forces separate from constraint
projection. Standard profile and unsupported spin/FD modes reject. Retry a failed
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
