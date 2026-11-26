# Getting started

Moved here from the project README (rewritten as a project overview on 2026-11-03) and refreshed to the current code. Project overview: [README](../README.md).

Current routing and limits: [ADR 0047](decisions/0047-reference-driver-single-path.md)
and [CLI](cli.md). The production reference driver supports NC LDA-PZ/LDA-PW/PBE,
Gamma and complex k points, NLCC, eligible collinear/FD models, symmetry,
derivatives, unpolarized fixed-occupation relaxation, and runtime/checkpoint v2.
Free-moment PBE remains unqualified and routes to StrictOracle. Optional MPI is
implemented; GPU, PAW/USPP, hybrids and SOC remain outside the current path.
The stored Si/QE comparisons are scoped evidence, not general material accuracy.

## Build and check

Use Rust 1.99.1 (pinned in [rust-toolchain.toml](../rust-toolchain.toml)), rustfmt,
clippy and a working native linker. The verified host is Apple Silicon macOS
with Xcode Command Line Tools; other platforms are not established by this
record. See [environment](environment.md) for the offline toolchain setup.
The local Cargo configuration uses a shared, untracked vendor directory; it must
be available for offline builds (see [dependencies](dependencies.md)).
Python 3.9+ is needed for the standard-library comparison tools.

From the repository root:

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export CARGO_NET_OFFLINE=true
scripts/check.sh
# The four gates, individually:
cargo build --workspace
cargo test --workspace
cargo fmt --check
cargo clippy --workspace --all-targets -- -D warnings
```

Long numerical tests are ignored in the fast debug gate; their release commands
are linked from [architecture](architecture.md).

## Run a workflow example

Start with the [13-workflow index](../examples/README.md), covering Gamma LDA,
k meshes/symmetry, PBE, bands/DOS, metals, spin, forces/relaxation, variable cell,
checkpoint/restart, observation, Python, MPI and the strict validation oracle.
Every tutorial has an input, UPF source, commands and output interpretation.
All debug example checks are validate-only; laptop runtime estimates are unmeasured.

```sh
cargo build --release -p dft-cli
# Obtain the named external UPF first; this is a helper variable, not a CLI setting.
export PSEUDO_DIR=/absolute/path/outside/repository/pseudo
python3 examples/prepare.py 01-si-gamma-lda
input=target/tmp/examples/01-si-gamma-lda/input.toml
target/release/baandi inspect-upf "$PSEUDO_DIR/Si.pz-vbc.UPF"
target/release/baandi validate "$input"
target/release/baandi grid "$input"
target/release/baandi run "$input" --output-dir results/si-gamma
```

The helper (Python 3.9+) copies the template into `target/tmp` with absolute UPF
paths. Repeating it overwrites the prepared copy. Download files outside the
repository; no UPF is distributed. File names and sources are in each example.
Check licenses and [supported UPF features](upf-support.md); the historical
Si validation file's provenance/license remain unresolved ([candidate notes](candidate-si-upf.md)).

The old [si_lda.toml](../examples/si_lda.toml), expression and k-grid inputs are
retained at their original paths. Replace their UPF placeholders manually.
`validate` opens UPFs and checks full preflight; `grid` inspects geometry/grid.
Relative UPF paths resolve against the TOML directory, output paths against the
working directory. TOML does not expand `$PSEUDO_DIR`.

An ordinary scalar run writes `results.json`, `density.bin` (little-endian f64)
and `density.json` (metadata/checksum). Existing ordinary outputs are overwritten;
use fresh directories for each tutorial. Exit 0 means convergence and successful
writes, with verification on by default; exit 3 means budget exhaustion.
Checkpoint v2 restart instead writes an export container and routing companion.
See [CLI](cli.md) for units and exit codes, and [user guides](user/) for workflow limits.

## Reproduce the comparisons

QE runs are executed remotely on dev08 through the master session, using the
[reference protocol](reference-environment.md); do not install or run QE
locally. The [reference bundles and tooling guide](../tests/reference/README.md)
define the required inputs, hashes, provenance and comparison gates.

```sh
python3 -I scripts/make_qe_input.py --results results/results.json \
  --pseudo-dir /remote/pseudo --upf-name Si.pz-vbc.UPF --out /tmp/qe.in
# Recheck the stored real bundle; this does not run QE:
python3 -I scripts/verify_qe.py --ours tests/reference/si30_ours/ours \
  --bundle tests/reference/si30_ours --json /tmp/si-comparison.json \
  --markdown /tmp/si-comparison.md
python3 -I scripts/compare_rhog.py --ours-dir tests/reference/si30_ours/ours \
  --qe-summary tests/reference/si30_ours/density/qe_rhoG_summary.json \
  --json-report /tmp/si-density.json --markdown-report /tmp/si-density.md
```

Use fresh density-report paths. [convergence_scan.py](../scripts/convergence_scan.py)
runs isolated cutoff/grid/SCF studies and can emit QE inputs; its invocation and
measured results are in [validation](validation.md). The example needs a
real UPF path before a scan. Regression tests use stored bundles without QE:
`python3 -I -m unittest discover -s scripts/tests -v`.

## Working in this repository

This is an AI-first repository. Follow [AGENTS.md](../AGENTS.md): read it first,
then [docs/INDEX.md](INDEX.md), [ownership.toml](../ownership.toml), relevant
crate guides, the [ticket](../tickets/) and relevant ADRs. Crates separate units,
geometry, numerics, basis, UPF, physics, solvers, engine and CLI; see
[architecture](architecture.md). Delegated runs use
[scripts/delegate.sh](../scripts/delegate.sh) and the session registry
[artifacts/sessions.json](../artifacts/sessions.json).

[ADR 0003](decisions/0003-prefer-general-purpose-libraries.md) prefers
maintained general-purpose libraries over hand-written numerics or FFI;
project code supplies physics glue. Versions, consumers and licenses are in
[dependencies](dependencies.md).

The project is licensed under the MIT License (see [LICENSE](../LICENSE)), selected by the user on 2026-11-03. The Si pseudopotential used for validation is not part of the repository and is not covered by this license.

