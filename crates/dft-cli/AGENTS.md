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
