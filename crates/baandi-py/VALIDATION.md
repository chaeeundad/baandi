# P11G validation and handoff

Status: BLOCKED on the additional `baandi-py -> dft-units` dependency decision.
The ticket owns the crate but `ownership.toml` currently declares only the
engine dependency. The user was asked asynchronously; no approval was received.
`units-integration.patch` is a concrete, unapplied proposal. Applying it also
requires the architect to synchronize ownership and Cargo.lock. No new external
library is needed. Do not replace it with copied conversion constants: ASE's
CODATA version need not match the CLI's uom conversion.

## Implemented and measured

- Native Options, Session and owned Result; initialized density summary and
  terminal energy/eigenvalues. Long calls release the GIL. Each worker owns its
  non-Send engine session from construction to destruction. No unsafe changes.
- Native run versus stepwise equality on cubic two-ion Z=1 and skewed one-ion
  Z=2 synthetic systems, with different projector sets. These are nonphysical
  parser/operator fixtures, not real pseudopotentials or QE references.
- CLI versus native energy and eigenvalues are bit-identical on both systems.
  `native-evidence.json` contains the measured Hartree outputs, hex encodings,
  charges and iteration counts. Neither tolerances nor defaults were relaxed.
- ASE unsupported-property requests raise NotImplementedError. The pending
  energy path raises an explicit error until the units bridge is available.
- macOS arm64, CPython 3.9, PyO3 0.30.1, maturin 1.16.0: mixed-package wheel
  `target/wheels/baandi-0.1.0-cp39-cp39-macosx_11_0_arm64.whl` builds and imports.
  ASE is an optional extra. No installation into the shared virtualenv occurred.

## Reproduction and environment workarounds

The local `.cargo/config.toml` replaces crates.io with a vendor directory that
lacks PyO3. The registry cache already has PyO3 and all required crates. Use a
small temporary Cargo wrapper which runs the real Cargo from `/private/tmp`
with the absolute manifest passed by maturin. This avoids changing shared
configuration or downloading anything:

```sh
# /private/tmp/p11g-tools/cargo (executable)
#!/bin/sh
cd /private/tmp || exit 1
exec /opt/homebrew/opt/rustup/bin/cargo "$@"
```

From `python/baandi`, run:

```sh
PATH=/private/tmp/p11g-tools:$PATH \
PYO3_PYTHON=/Users/mino/Develop/rust-dft-wt/pyenv/bin/python \
/Users/mino/Develop/rust-dft-wt/pyenv/bin/maturin build --offline --locked
```

Do not pass `--manifest-path` to maturin: it would select a different Python
project root. The pyproject points to the Rust crate. The tested wheel was
installed with `pip --no-deps --no-index --target /private/tmp/p11g-wheel-test`.
Run tests from the repository checkout (the synthetic generator is a repository
test utility, not a wheel runtime dependency):

```sh
PYTHONPATH=/private/tmp/p11g-wheel-test \
BAANDI_CLI="$PWD/target/debug/dft-cli" \
BAANDI_NATIVE_EVIDENCE="$PWD/crates/baandi-py/native-evidence.json" \
/Users/mino/Develop/rust-dft-wt/pyenv/bin/python \
-m unittest discover -s python/baandi/tests -v
```

With the current unapplied bridge, four test methods pass and the ASE/CLI test
errors on the missing unit conversion function; this is an unmet acceptance
criterion, not a skipped success. Without the extension, all tests skip cleanly.
The dedicated native/CLI test passes independently. After applying the bridge,
run the full suite with `BAANDI_EVIDENCE` pointing at a separate ASE evidence
JSON. No ASE energy equality has been claimed yet.

The ordinary `scripts/check.sh baandi-py` exited 1 (vendor lacks PyO3). An offline
cache wrapper ran the same script successfully (exit 0; build/test/fmt/clippy
all PASS). The script invocation is saved at `/private/tmp/p11g-check.sh` for
this session. Its exported Bash `cargo` function changes cwd to `/private/tmp`
and inserts the absolute workspace manifest before forwarding arguments;
`fmt` additionally passes `--all`. It sets:

- `PYO3_PYTHON=/Users/mino/Develop/rust-dft-wt/pyenv/bin/python`
- `CARGO_NET_OFFLINE=true`
- `LIBRARY_PATH=/Library/Developer/CommandLineTools/Library/Frameworks/Python3.framework/Versions/3.9/lib`
- `RUSTFLAGS=-C link-arg=-Wl,-rpath,/Library/Developer/CommandLineTools/Library/Frameworks`

The last two fix this Python installation's stale Xcode library path and allow
the Rust test executable to load the real framework. Initial cache attempts
failed linking and loading Python before these workarounds. Cargo test runs
zero Rust tests here; functional validation is the Python suite, not that gate.

## Outstanding integration

The architect must approve/synchronize the units dependency, apply the patch,
rebuild and run the ASE equality test. README result generation and dashboard
registry/dataset updates are outside this ticket's OWNED_FILES and were not
changed. The local session registry has no P11G run and no `.delegate` report;
no recurring orchestrator timer was started. The task title already exists.
No remote changes, UPF commits, git metadata writes or core edits were made.
The working tree is intentionally uncommitted.
