# P15C spectral completeness evidence

See [scope, findings and limitations](../../../docs/validation/spectral-completeness.md).
Current runtime implementation and measurements live in [`final/`](final/).
`[scf] spectral_guard = "final"` checks every reported k/spin state after the
fixed-point certificate passes, records cumulative repairs and fails closed.
Default is `"off"` because final-only work exceeds the 3% target. Off emits no
new top-level block and preserves schema_version. See the report for the seeded
large-space check's limits. Root-level `verification.json` and candidate records
are historical first-stage evidence, not the current implementation status.

Reproduce the original oracle scope (source-only diagnostic build):

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
cargo build --offline --release -p dft-cli
cp target/release/dft-cli /tmp/p15c-baseline-cli
P15C_ORACLE=1 cargo build --offline --release -p dft-cli
cp target/release/dft-cli /tmp/p15c-oracle-cli
python3 -I tests/reference/spectral/scope.py --binary /tmp/p15c-oracle-cli
python3 -I tests/reference/spectral/build_spin_fixture.py
python3 -I tests/reference/spectral/spin_scope.py /tmp/p15c-oracle-cli
python3 -I tests/reference/spectral/scope.py --binary /tmp/p15c-oracle-cli --cases si-30,si-40,mgo-50,mgo-60
python3 -I tests/reference/spectral/summarize.py
python3 -I tests/reference/spectral/test_scope.py
```

Scope inputs use the existing benchmark generators. DFT_NC_DIR defaults to the
local `../pseudo/pseudodojo/ex` directory, whose fixture hashes are checked.
The runner samples initial/mid/final ordinary solves and the final certificate,
at k indices 0, floor(nk/2), nk-1. All dimensions <=128 are checked every time.
A source-only `P15C_DENSE_IDS` environment variable holds exact sequential solve
IDs, stored per case. Stderr is parsed only after the run, outside numeric loops.

`P15C_ORACLE=1 P15C_CANDIDATES=1` measures warm-plus-random and independent-cold
starts on each selected frozen H. It returns only the unchanged native solution.
Use `--output tests/reference/spectral/candidates` to retain separate evidence.

Current final-policy checks use an ordinary release binary with diagnostic flags
unset:

```sh
cargo build --release --offline -p dft-cli
python3 -I tests/reference/spectral/final_policy.py --binary target/release/dft-cli --threads 1 --output /tmp/p15c-final-recheck
cargo test --release -p dft-solvers spectral_ -- --include-ignored --nocapture
P15C_ORACLE=1 cargo test --release -p dft-cli --test materials final_state_dense_spectra -- --ignored --nocapture
python3 -I tests/reference/spectral/final_summary.py
```

The ignored materials test requires readable `DFT_NC_DIR` with the P13D fixture
hashes and checks all six Si/MgO cutoffs. `P15C_DENSE_SCF=1 P15C_ORACLE=1` selects
the dense-every-step diagnostic used only for the small skew trajectory.
`P15C_CANDIDATES`, `P15C_FULL_GUARD`, and `P15C_RICH_START` remain source-only
experiments. The old `P15C_GUARD`/`P15C_UNGUARDED` flags no longer control code.

Stored files:

- `<case>.json`: input/hashes/commands, baseline identity, per-solve comparisons.
- `<case>-traces.json.gz`: full dense spectrum plus every native Ritz/residual
  row; candidate traces have variants 1 (warm+random) and 2 (cold random).
- `summary.json`, `table.md`: generated final ordinary comparisons. Per-k band
  sums are unweighted sample values, not an inferred full-mesh energy error.
- `spin.json`: exact P13B Gaussian fixture generated using the existing Rust
  test helper; `spin-generator-build.json` records the standalone link command.
  The earlier Python-erf approximation is retained under `candidates/`.
- `baseline/guard/screen/rich-t*.json`: original and candidate work/observational
  wall counts. `default-t*.json` and `determinism-commands.json` check the unchanged
  default at four thread counts. No speedup is inferred from overlapping runs.
- `guard-synthetic.json`, `rich-synthetic.json`: complete candidate SCF results,
  including nonconvergence. Exit 3 means exhausted, not passed ground state.
- `*-robustness.json`: existing robustness suite reports; suite exit 0 does not
  mean every synthetic case converged. Individual statuses are retained.
- `verification.json`, `logs/`: exact checks and failed attempts, source hashes
  and executable hashes. Binaries and UPF contents are never committed here.

The default checks must be run with all P15C diagnostic flags unset.
Current CLI tests cover additive output, unchanged schema_version, disabled
block absence, both spins, and invalid options. A failed or exhausted final-mode
run has `passed=false`; guard-off convergence is not a completeness claim.
Current final guard sessions reject checkpoint v1 rather than lose audit history.
