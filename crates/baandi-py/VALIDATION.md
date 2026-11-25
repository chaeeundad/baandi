# P11G validation and handoff

Implementation complete after the approved `baandi-py -> dft-units` dependency
was merged. The unit conversion patch has been applied to `src/lib.rs` and the
patch file removed. ASE uses the CLI's exact energy conversion operation and
Angstrom-to-Bohr wrapper, without copying constants or using ASE's CODATA values.

## Measured acceptance

- ARM64 CPython 3.9 wheel builds with PyO3 0.30.1 and maturin 1.16.0 using the
  offline vendor directory. No Cargo wrapper or network access is required.
- Five Python test methods pass with no skips (1.066 s). They exercise native
  stepwise/batch equality, owned results, density summaries, invalid inputs,
  exhausted results, unsupported ASE properties and native/ASE/CLI equality.
- ASE and CLI energy values are bit-identical in Ha and eV on cubic two-ion
  Z=1 and skewed one-ion Z=2 systems with different projector sets. Eigenvalues
  also match. `ase-evidence.json` records the eV values and hexadecimal equality.
- `native-evidence.json` independently records native/CLI comparisons using
  cells specified in Bohr. The ASE cases specify different cells in Angstrom,
  so their total energies differ from those native-only cases.
- Fixtures are synthetic nonphysical NC-format data, not real-element
  pseudopotentials or QE references. No defaults or tolerances were relaxed.
- Without an importable extension, all five tests skip cleanly (exit 0).

The wheel is `target/wheels/baandi-0.1.0-cp39-cp39-macosx_11_0_arm64.whl`.
It was installed into `/private/tmp/p11g-wheel-test`, not the shared virtualenv.
ASE remains an optional extra. Unsupported forces, stress and free energy raise
NotImplementedError. Exhausted calculations raise ConvergenceError without
publishing an ASE energy; uncertified native diagnostics remain accessible.

## Commands

Run from the repository root unless a different directory is indicated:

```sh
export PYO3_PYTHON=/Users/mino/Develop/rust-dft-wt/pyenv/bin/python
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
cargo build -p dft-cli --offline --locked
# Exit 0.

# From python/baandi; do not override maturin's manifest/project discovery:
/Users/mino/Develop/rust-dft-wt/pyenv/bin/maturin build --offline --locked
# Exit 0, macosx_11_0_arm64 wheel.

# Back at repository root:
"$PYO3_PYTHON" -m pip install --no-deps --no-index --upgrade \
  --target /private/tmp/p11g-wheel-test \
  target/wheels/baandi-0.1.0-cp39-cp39-macosx_11_0_arm64.whl
# Exit 0.

PYTHONPATH=/private/tmp/p11g-wheel-test \
BAANDI_CLI="$PWD/target/debug/dft-cli" \
BAANDI_NATIVE_EVIDENCE="$PWD/crates/baandi-py/native-evidence.json" \
BAANDI_EVIDENCE="$PWD/crates/baandi-py/ase-evidence.json" \
"$PYO3_PYTHON" -m unittest discover -s python/baandi/tests -v
# Exit 0: five passed, zero skipped.
```

The tests use the repository's synthetic generator and must run from a checkout.
No UPFs are committed; fixtures exist only in temporary directories.

## Python-independent workspace gates

Both plain commands were run with no environment overrides and exited 0:

```sh
scripts/check.sh
scripts/check.sh baandi-py
```

Build, test, fmt and clippy all PASS in both runs. No `PYO3_PYTHON`,
`LIBRARY_PATH`, `RUSTFLAGS`, Cargo wrapper or Python development library is
required for the default workspace build. The earlier Xcode libpython linking
failure is resolved by making the Python ABI opt-in, not by changing the host.

`pyo3` is optional and the `python-extension` Cargo feature has no default
activation. Maturin enables it via `python/baandi/pyproject.toml`; it enables
PyO3's `extension-module` feature, suppressing libpython linkage. The crate's
build script emits macOS `-undefined dynamic_lookup` for the cdylib only,
equivalent to `pyo3_build_config::add_extension_module_link_args` on macOS.
The script needs no Python-probing build dependency. Python symbols are resolved
from the loading interpreter. Default Cargo artifacts contain no Python ABI;
use maturin for a usable Python module.

Rust unit-test builds omit the Python binding module, so no interpreter is
embedded. The adapter has zero Rust unit tests; the installed-wheel Python
suite exercises its actual FFI implementation. The full workspace still runs
all other Rust crate tests. To lint the enabled binding explicitly:

```sh
cargo clippy -p baandi-py --features python-extension --all-targets --offline --locked -- -D warnings
```

The feature-enabled clippy command also exited 0.

Additional inspection of the installed wheel with `otool -L` confirms that its
only system dylib dependencies are libc++, libiconv and libSystem: no libpython
or Python framework. `nm` shows unresolved Python symbols, and the successful
installed-wheel import and five tests verify runtime resolution. The default
`cargo tree -p baandi-py --edges normal` contains no PyO3 dependency. No manifest
lockfile change, system-library removal or interpreter configuration edit was
needed.

## Scope and handoff

Gamma-only LDA is supported here. Audited policy interventions, full workflow
schemas, forces/stress and other wheel platforms remain outside this binding.
The core is unchanged. All changes are intentionally uncommitted.
README generation and dashboard registry/dataset updates remain with the
orchestrator outside this ticket's OWNED_FILES. The local registry now contains
a P11G running entry; no local delegate report or recurring timer was created.

## ONEPATH-4 re-run (2026-11-16, Mac, debug wheel)

After the default flip (ADR 0043 amendment ONEPATH-4) the wheel was rebuilt with
`maturin build --offline --locked` and installed into a scratch target; with
`BAANDI_CLI` set to the same-tree debug `dft-cli`, `python -m unittest discover
-s python/baandi/tests -v` ran 9 tests, 9 passed, 0 skipped. Native, CLI and ASE
energies are bit-identical on the default standard path (CLI with
`cutoff_ladder = "off"`); `native-evidence.json` and `ase-evidence.json` were
regenerated by that run (standard path: fewer SCF iterations; the four
energies moved by 1.1e-12, 1.2e-12, 9.8e-12 and 7.6e-11 Ha from the former
strict values, inside the 5e-9 Ha standard class; loosely converged
empty bands beyond the first differ, as the standard band policy allows). The
ASE stress oracle runs at `standard_threshold_ha = 5e-11` (max FD error 8.8e-8
eV/Angstrom^3, bound 1e-5); a default-threshold session refuses `stress()`.


## R5 reference batch migration (2026-11-25)

The approved `dft-engine::input::kpoints` re-export lets the binding use the
same mesh and time-reversal helpers as the CLI without a new dependency.
`Options.kgrid` defaults to `[1, 1, 1]` (Gamma). Other positive integer meshes
use unshifted Monkhorst-Pack sampling and time-reversal reduction, matching CLI
`basis.kgrid` defaults. Gamma selects the real driver; other meshes select the
complex k driver. Legacy sessions remain Gamma-only and reject other meshes.

Native acceptance: all **6 tests passed, zero skips**, using the CPython 3.9
debug extension and same-build CLI. New measured evidence is in
[`tests/r5-reference-evidence.json`](tests/r5-reference-evidence.json); it records
the input TOML, CLI command/tool identity, external UPF SHA-256, and both sides'
hexadecimal results. Every comparison uses exact float bits, with no tolerance:

| Si2 mesh | Energy (Ha) | SCF iterations | Max energy difference (Ha) | Max force difference (Ha/Bohr) | Max stress difference (Ha/Bohr^3) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gamma | -7.233591259346283 | 8 | 0 | 0 | 0 |
| 2×2×2 | -7.776922790383691 | 7 | 0 | 0 | 0 |

These are displaced Si2 adapter-parity measurements at a low cutoff, not
converged physical predictions or new QE references. The native suite also
checks defaults/overrides, invalid options, explicit stepwise gating, rejection
of unsupported atomic-start data without fallback, legacy mesh rejection and
legacy batch/step/owned-result equivalence.

The measurements above and `native-evidence.json` / `ase-evidence.json` are
historical **legacy-driver** evidence. They are preserved, not re-labelled as
reference results. R5 changes the Python default driver to `"reference"`.

`Options()` now selects atomic density and wavefunction starts, reference
Broyden mixing beta **0.7**, history **8**, standard threshold **5e-9 Ha**, and
fixed-point verification **on**, matching the CLI reference facade. Explicit
mixing overrides are passed into the reference profile (threshold in Ry = twice
the Ha option). `driver = "legacy"` retains the prior session implementation,
Pulay beta **0.5**, history **6**, and stepwise control. Strict legacy defaults
use beta **0.3**. Mixing defaults resolve when read/used, so changing `driver`
restores the appropriate defaults without overwriting explicit overrides.

`reference_strict = True` (including its existing environment-derived default)
requires `driver = "legacy"`; it never silently redirects an explicit reference
session. Reference inputs need the atomic UPF data accepted by `refstart`;
unsupported inputs fail, with no legacy retry. The former synthetic tests have
zero atomic density and no atomic wavefunctions; they remain legacy fixtures.

`run()` calls `refscf::api::solve` once and retains its delivery. `result()` and
`density_summary()` read that delivery, while `forces()` and `stress()` use the
facade's deferred observable preparation/evaluation. Stress called before forces
still follows prepare → forces → stress. Repeated calls reuse the delivered
state. The GIL is released while the worker calculates. Reference `verify=False`
uses the facade's explicit unverified opt-out, as on the CLI.

Until the reference stepper (R8), `initialize()`, `step()`, `status` and
`iteration` raise exactly:

```text
stepwise control requires driver='legacy' until the reference stepper (R8)
```

Select `standard_threshold_ha = 5e-11` or tighter **before** creating a session
whose stress will be requested. Python cannot infer a later stress request at
construction, whereas CLI `[output] stress=true` selects this threshold by
default. The returned stress is tensile-positive, frozen-Miller Ha/Bohr^3, with
the existing separate Pulay warning; it is not a cutoff-convergence claim.

New native tests are in `tests/test_reference.py`, discoverable by pytest or
standard-library unittest. They skip when the extension is unavailable. Si2
parity additionally requires `BAANDI_CLI` and `BAANDI_SI_UPF` pointing to a
same-build CLI and an external Si.pz-vbc.UPF. The two cases use the same
slightly displaced Si2 cell with Gamma and an unshifted 2×2×2 mesh, 4 Ha cutoff,
8 bands, 5e-11 Ha threshold, and verification on. They assert the CLI driver is
reference and require hexadecimal identity for energy, all six force components
and all nine stress components. This is adapter parity, not independent physics
validation or a QE comparison. Stress adds three cutoff SCF scans; run these
explicitly outside the fast Cargo gate.

```sh
# Build CLI and extension with the same debug profile/toolchain.
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
cargo build -p dft-cli -p baandi-py --features baandi-py/python-extension --offline --locked
# Load a wheel built with maturin, or stage the cdylib in a scratch baandi package.
PYTHONPATH=/path/to/scratch/package BAANDI_CLI="$PWD/target/debug/dft-cli" \
BAANDI_SI_UPF=/path/to/Si.pz-vbc.UPF \
BAANDI_R5_EVIDENCE="$PWD/crates/baandi-py/tests/r5-reference-evidence.json" \
python -m unittest discover -s crates/baandi-py/tests -v
```

Default `cargo test -p baandi-py` still has zero binding tests because PyO3 is
omitted from Rust test builds. It cannot establish native acceptance. Always run
the feature-enabled check/clippy and the native suite as well. README generation,
existing top-level `python/baandi` compatibility tests/ASE option forwarding, and
registry/dashboard updates are outside this ticket's owned files and remain
orchestrator follow-up.

R5 local checks (2026-11-25):

- `cargo check --workspace --tests`: exit 0.
- `cargo test -p baandi-py`: exit 0, zero Rust tests and zero doctests.
- `cargo fmt --check && cargo clippy -p baandi-py --all-targets -- -D warnings`: exit 0.
- Feature-enabled Clippy (`--features python-extension --all-targets`): exit 0.
- Native build (`cargo build -p baandi-py -p dft-cli --features
  baandi-py/python-extension --offline --locked`): exit 0. CPython 3.9 loaded
  the cdylib from a scratch package under `target/r5-python`; no wheel-install
  acceptance is claimed.
- Native unittest command above, with the external Si UPF and built CLI:
  exit 0, 6 passed, zero skipped (1.304 s).

Logs are local, ignored files under `target/r5-*.log`. Native assertions and the
new stored evidence establish the adapter comparisons; zero-test Cargo passes
alone do not establish native acceptance.

The full workspace gate was also run with
`CHECK_SCCACHE=0 BAANDI_TEST_TMP="$PWD/target/tmp" scripts/check.sh`:
exit **1**, build/fmt/clippy PASS; 107 test binaries, **1012 passed, 3 failed,
221 ignored**, 335 s total (42 s test wall). All three failures are existing
CLI cache tests whose successful-cache fixtures reject the sandbox's checkout-
local scratch directory with `pseudopotential cache must be outside repositories`:

- `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`
- `pseudo::cache_resolution_and_no_network_side_effects`
- `pseudo::cached_fetch_checks_hash_and_provenance`

A direct rerun of the first test under the same scratch root exited 101 and
confirmed the same path rejection at `crates/dft-cli/src/pseudo.rs:611`. No cache
policy or unrelated tests were changed. The sandbox also denies the gate's
optional `ps` load probe and `nice` priority adjustment. These limitations are
separate from the six passing native binding tests. The no-extension run skips
all six methods cleanly (exit 0), which is not numerical acceptance.

The final scoped gate,
`CHECK_SCCACHE=0 BAANDI_TEST_TMP="$PWD/target/tmp" scripts/check.sh baandi-py`,
exited **0**: build/test/fmt/clippy all PASS (2 s). It remains a zero-Rust-test
crate gate; the separate native suite supplies the binding acceptance.


## R8d reference stepping (2026-11-25)

`tests/test_reference.py` now exercises the reference runtime through native
`initialize`, `step`, `status`, `iteration` and `run`. `step` returns all fields
of an owned `StepReport`. Reference status names are Created, Initialized,
Converged, Exhausted and Failed; legacy status names and step strings remain.
The explicit `density()` accessor returns `(shape, values)` in ndarray iteration
order, in Bohr^-3; it copies only on request and requires reference delivery.

The same-build native check passed six tests; the existing external Si/CLI test
was skipped (its external input variables were not provided). Gamma and k211
synthetic one-atom cases compare stepping, partially stepped `run()` and fresh
`run()` bit-for-bit: total energy, all eigenvalues and all density values, plus
forces for stepped versus fresh delivery. They take 8 maps / 7 iterations and
7 maps / 6 iterations respectively, including verification. New evidence is
`tests/r8d-reference-evidence.json`; these are not QE or real-material claims.
The extension-absent invocation skips all seven tests cleanly.

Reproduction from the worktree root (no wheel installation is required):

```sh
cargo build -p baandi-py --features python-extension
mkdir -p target/tmp/r8d-native/baandi
cp python/baandi/__init__.py target/tmp/r8d-native/baandi/__init__.py
cp target/debug/lib_baandi.dylib target/tmp/r8d-native/baandi/_baandi.so
PYTHONPATH="$PWD/target/tmp/r8d-native" \
BAANDI_R8D_EVIDENCE="$PWD/crates/baandi-py/tests/r8d-reference-evidence.json" \
python3 -m unittest discover -s crates/baandi-py/tests -p test_reference.py -v
```

The dylib copy is the macOS native-development route; this run does not claim
installed-wheel, cross-platform or external Si acceptance. The module must be
built for the Python interpreter used to load it.

R8d verification exits were all 0: `cargo check --workspace --tests`,
`cargo test -p baandi-py` (no native tests in the Rust harness),
`cargo test -p baandi-mcp` (13 server tests),
`cargo test -p dft-cli --test driver_routing` (8 tests), `cargo fmt --check`,
`cargo clippy -p baandi-py -p baandi-mcp --all-targets -- -D warnings`, and
`cargo clippy -p baandi-py --features python-extension --all-targets -- -D warnings`.
`BAANDI_TEST_TMP="$PWD/target/tmp" CHECK_SCCACHE=0 CHECK_TEST_JOBS=2
scripts/check.sh fast baandi-py baandi-mcp` exited 0 with build, test, fmt,
clippy, production-check and production-clippy PASS (110 s total). Its optional
`ps` load probe was denied by the sandbox; verification gates still completed.
Root README/dashboard integration is left to the orchestrator because those
files are outside this ticket's OWNED_FILES. Existing legacy evidence is retained.
