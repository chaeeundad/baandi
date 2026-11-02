# rust-dft

Rust-first CPU plane-wave DFT engine for Gamma-only, spin-unpolarized,
norm-conserving, insulating LDA-PZ self-consistent calculations. The current
implementation uses complex wavefunctions, Davidson diagonalization, and linear
or Pulay density mixing, with an explicit total-energy audit.

Si diamond has been compared with Quantum ESPRESSO 7.6. The
[validation record](docs/validation.md) reports passing energy comparisons at
1e-4 Ha/atom, eigenvalues and gap at 1e-3 Ha, and electron count at 1e-8 relative.
It includes a cutoff scan and a Fourier-density comparison of sampled
coefficients and aggregates; full-set density agreement is not certified.
Evidence covers one cell and one pseudopotential, not general material accuracy.
Final MVP approval remains pending.

Out of scope until MVP approval: CUDA, Metal, MPI, PAW, ultrasoft
pseudopotentials, PBE, forces, stress, smearing and multiple k-points.

## Build and check

Use Rust 1.99.1 (pinned in [rust-toolchain.toml](rust-toolchain.toml)), rustfmt,
clippy and a working native linker. The verified host is Apple Silicon macOS
with Xcode Command Line Tools; other platforms are not established by this
record. See [environment](docs/environment.md) for the offline toolchain setup.
The local Cargo configuration uses a shared, untracked vendor directory; it must
be available for offline builds (see [dependencies](docs/dependencies.md)).
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
are linked from [architecture](docs/architecture.md).

## Run the Si example

```sh
cargo run --offline -p dft-cli -- validate examples/si_lda.toml
cargo run --offline -p dft-cli -- grid examples/si_lda.toml
# Supply your UPF and replace the placeholder path in examples/si_lda.toml first:
cargo run --offline -p dft-cli -- inspect-upf /absolute/path/to/Si.UPF
cargo run --offline --release -p dft-cli -- run examples/si_lda.toml --output-dir results
```

[si_lda.toml](examples/si_lda.toml) is the minimal input. `validate` and `grid`
check configuration and basis/grid settings without opening the UPF.
No UPF file is distributed: supply a supported pseudopotential yourself and
check its license. The Si file used for validation has unverified provenance
and an unknown license; its user-approved use was for internal validation only,
with no redistribution. See [candidate notes](docs/candidate-si-upf.md) and
[UPF support](docs/upf-support.md).

Relative UPF paths resolve against the input directory. Output defaults to
`results/` relative to the working directory; `--output-dir` overrides it.
The CLI writes `results.json` (energies, convergence history and provenance),
`density.bin` (little-endian f64 density) and `density.json` (metadata/checksum).
Existing output files are overwritten. Only certified convergence with successful
writes returns 0; budget exhaustion writes the final state and returns 3.
See [CLI documentation](docs/cli.md) for controls, units and all exit codes.

## Reproduce the comparisons

QE runs are executed remotely on dev08 through the dft5090 session, using the
[reference protocol](docs/reference-environment.md); do not install or run QE
locally. The [reference bundles and tooling guide](tests/reference/README.md)
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

Use fresh density-report paths. [convergence_scan.py](scripts/convergence_scan.py)
runs isolated cutoff/grid/SCF studies and can emit QE inputs; its invocation and
measured results are in [validation](docs/validation.md). The example needs a
real UPF path before a scan. Regression tests use stored bundles without QE:
`python3 -I -m unittest discover -s scripts/tests -v`.

## Working in this repository

This is an AI-first repository. Follow [AGENTS.md](AGENTS.md): read it first,
then [docs/INDEX.md](docs/INDEX.md), [ownership.toml](ownership.toml), relevant
crate guides, the [ticket](tickets/) and relevant ADRs. Crates separate units,
geometry, numerics, basis, UPF, physics, solvers, engine and CLI; see
[architecture](docs/architecture.md). Delegated runs use
[scripts/delegate.sh](scripts/delegate.sh) and the session registry
[artifacts/sessions.json](artifacts/sessions.json).

[ADR 0003](docs/decisions/0003-prefer-general-purpose-libraries.md) prefers
maintained general-purpose libraries over hand-written numerics or FFI;
project code supplies physics glue. Versions, consumers and licenses are in
[dependencies](docs/dependencies.md).

No project license has been selected yet; this is pending a user decision.
