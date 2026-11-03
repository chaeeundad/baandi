# Getting started

Moved here from the project README (rewritten as a project overview on 2026-11-03) and refreshed to the current code. Project overview: [README](../README.md).

Current capabilities: CPU plane-wave Kohn-Sham DFT, spin-unpolarized, norm-conserving, insulating LDA-PZ, Bloch k points (Gamma, Monkhorst-Pack grids, explicit weighted lists; time-reversal reduction only), units and arithmetic expressions in the input (`"30 Ry"`, `"1/4"`). Si has been compared with Quantum ESPRESSO 7.6 at Gamma and on three k grids ([validation](validation.md) R3 to R6). Out of scope for now: CUDA, Metal, MPI, PAW, ultrasoft pseudopotentials, PBE, forces, stress, smearing, spin (see [AGENTS.md](../AGENTS.md) and [p7-plan](p7-plan.md)). Evidence covers one cell and one pseudopotential, not general material accuracy.

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

## Run the Si example

```sh
cargo run --offline -p dft-cli -- validate examples/si_lda.toml
cargo run --offline -p dft-cli -- grid examples/si_lda.toml
# k-point grid and expression inputs: examples/si_lda_kgrid.toml, examples/si_lda_expr.toml
cargo run --offline -p dft-cli -- convert 30 Ry Ha     # unit conversion
cargo run --offline -p dft-cli -- eval "atan(0.5*pi)"   # expression check
# Supply your UPF and replace the placeholder path in examples/si_lda.toml first:
cargo run --offline -p dft-cli -- inspect-upf /absolute/path/to/Si.UPF
cargo run --offline --release -p dft-cli -- run examples/si_lda.toml --output-dir results
```

[si_lda.toml](../examples/si_lda.toml) is the minimal input. `validate` and `grid`
check configuration and basis/grid settings without opening the UPF.
No UPF file is distributed: supply a supported pseudopotential yourself and
check its license. The Si file used for validation has unverified provenance
and an unknown license; its user-approved use was for internal validation only,
with no redistribution. See [candidate notes](candidate-si-upf.md) and
[UPF support](upf-support.md).

Relative UPF paths resolve against the input directory. Output defaults to
`results/` relative to the working directory; `--output-dir` overrides it.
The CLI writes `results.json` (energies, convergence history and provenance),
`density.bin` (little-endian f64 density) and `density.json` (metadata/checksum).
Existing output files are overwritten. Only certified convergence with successful
writes returns 0; budget exhaustion writes the final state and returns 3.
See [CLI documentation](cli.md) for controls, units and all exit codes.

## Reproduce the comparisons

QE runs are executed remotely on dev08 through the dft5090 session, using the
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

