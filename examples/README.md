# Runnable workflow examples

These tutorials use the single production reference driver ([ADR 0047](../docs/decisions/0047-reference-driver-single-path.md)); example 13 explicitly selects the strict oracle. They demonstrate interfaces, not cutoff-converged material predictions. Check process exit status and verification before interpreting output.

## Setup

Run all commands from the repository root. Use the pinned Rust toolchain and the
build prerequisites in [getting started](../docs/getting-started.md).

```sh
cargo build --release -p dft-cli
export PSEUDO_DIR=/absolute/path/outside/repository/pseudo
python3 examples/prepare.py 01-si-gamma-lda
input=target/tmp/examples/01-si-gamma-lda/input.toml
target/release/baandi validate "$input"
target/release/baandi run "$input" --output-dir results/01-si-gamma-lda
```

Obtain the named UPF in each tutorial first; **no UPF is committed**. Download and
extract files outside this repository, preserve their licenses/attribution, and
inspect each with `target/release/baandi inspect-upf "$PSEUDO_DIR/FILE"`.
Use NC scalar/nonrelativistic UPF v2 with matching XC, not PAW/USPP or SOC.
See [current UPF restrictions](../docs/upf-support.md) and the
[pseudopotential guide](../docs/howto-pseudopotentials.md).

`PSEUDO_DIR` is a convenience variable for `prepare.py` (Python 3.9+), **not a CLI
search path**. The helper substitutes absolute UPF paths and writes only
`target/tmp/examples/WORKFLOW/input.toml`; it overwrites that prepared copy on
rerun. The CLI does not expand shell variables in TOML. You may instead copy the
input yourself and edit its paths. Relative UPF paths resolve against the TOML;
output directories resolve against the working directory. The registry cache's
`BAANDI_PSEUDO_DIR` is a separate setting and does not resolve these input paths.

All directories below contain source inputs/documentation only. Keep generated
results in `results/` or another scratch directory. Use fresh output directories,
especially for relaxation, volumetric exports and checkpoints. Exit 3 means budget
exhaustion, not success; other nonzero statuses and their causes are in [CLI](../docs/cli.md).
Laptop times in the tutorials are unmeasured release-build planning estimates,
excluding compilation. Real runtimes and convergence depend on UPF, hardware,
cutoff, mesh and electronic state. Heavy workflows are **validate-only** in debug
CI: `cargo test -p dft-cli --test examples_validate` substitutes synthetic
potentials and never executes SCF, NSCF or relaxation. Passing preflight does not
establish numerical convergence or physical accuracy.

## Workflows

| Example | Workflow |
| --- | --- |
| [01-si-gamma-lda](01-si-gamma-lda/) | Si Gamma LDA energy and eigenvalues |
| [02-si-kmesh](02-si-kmesh/) | Si Monkhorst–Pack mesh with symmetry |
| [03-si-pbe](03-si-pbe/) | Si PBE and NLCC |
| [04-bands-dos](04-bands-dos/) | Band path and density of states |
| [05-al-metal-fd](05-al-metal-fd/) | Al metal with Fermi–Dirac occupations |
| [06-fe-spin](06-fe-spin/) | Fe collinear spin with a fixed moment |
| [07-forces-relax](07-forces-relax/) | Displaced Si8 forces, stress and ion relaxation |
| [08-cell-relax](08-cell-relax/) | Si variable-cell relaxation |
| [09-checkpoint-restart](09-checkpoint-restart/) | Map-boundary checkpoint and restart |
| [10-observe](10-observe/) | Observe the SCF event stream |
| [11-python](11-python/) | Python batch, forces and stepwise API |
| [12-mpi](12-mpi/) | MPI with the same Si input |
| [13-strict-oracle](13-strict-oracle/) | Strict oracle for validation |

## Earlier inputs

The existing [si_lda.toml](si_lda.toml), [si_lda_expr.toml](si_lda_expr.toml) and
[si_lda_kgrid.toml](si_lda_kgrid.toml) remain at their original paths. They still
need their placeholder UPF paths replaced; `validate` now reads UPF contents.
The expression example demonstrates `5.43/2`, `1/4` and energy units.
Each lattice-vector row is a vector; positions are fractional and energies are
converted at the input boundary (30 Ry = 15 Ha).

Read more: [getting started](../docs/getting-started.md), [CLI](../docs/cli.md),
[user guides](../docs/user/), [Python](../python/baandi/) and
[host-embedded MCP](../crates/baandi-mcp/README.md).
