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
