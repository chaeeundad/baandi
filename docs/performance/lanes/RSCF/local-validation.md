# R-SCF-1 local validation

Machine: Mac, pure Rust backend. Source base:
`c534359c7ca6a7b2ee18694f905c453839bcb764`; changes uncommitted by instruction.
Numeric evidence, input hashes, commands and raw synthetic replay rows:
[local-validation.json](local-validation.json). The four adjacent Si TOMLs
are reproduction inputs with a UPF placeholder; replace that one path with the
hash-matched external Si.pz-vbc.UPF. No UPF or fabricated QE fixture is committed.

## Reproduce

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export TMPDIR_990="$PWD/target/scratch"
mkdir -p "$TMPDIR_990"
cargo test -p dft-engine --test refscf -- --nocapture
cargo test -p dft-cli --bin dft-cli refscf_check::tests -- --nocapture
cargo build --release -p dft-cli
python3 -I docs/performance/lanes/RSCF/synthetic-replay.py \
  --binary target/release/dft-cli --scratch "$TMPDIR_990/synthetic-replay"
CHECK_SCCACHE=0 scripts/check.sh full
python3 -I scripts/make_readme_figures.py
python3 -I -m unittest discover -s scripts/tests -p test_readme_figures.py
```

The engine target has ten tests; one reports an explicit runtime skip because
real QE Si8 snapshots are absent. The other nine exercise analytic density,
metric, Broyden, singularity, ring/TF/high-G rules, binary parsing and the
serial/simulated-MPI driver. The CLI unit test covers full analytic trajectory
replay and injected first-failure detection. The subprocess script additionally
checks actual CLI dispatch, stdout records, exit 0 on the synthetic oracle and
exit 1 naming `k=2 dr2` on the deliberately corrupted stdout. Synthetic maximum
error is 2.168404344971009e-19 e/bohr^3; no scientific QE validation is implied.

For a real new-driver smoke, add `driver="reference"` and `verify="off"` to
an input's existing `[scf]` table, then `dft-cli run INPUT --output-dir OUTPUT`.
With the external Si PP, the displaced Si8/Si16 cases stop in 10/9 maps and exit
0. Omit `verify="off"` and both return exit 3 at exactly the same estimator
stop, with `reference_stopped=true`, an unpassed verification report and the
explicit patch-2 TODO. Environment-only selection also passes; a k222 input
with the switch on exits 1 with the real-Gamma scope error.

## Gate limitations requiring the owner

The full gate's **build, fmt and clippy pass**. The test phase reports
**813 passed, 3 failed, 181 ignored**. All three failures are pre-existing
pseudopotential cache tests:

- `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`
- `pseudo::cache_resolution_and_no_network_side_effects`
- `pseudo::cached_fetch_checks_hash_and_provenance`

Their caches derive from TMPDIR but the cache code rejects paths in a repository.
The ticket mandates worktree-local scratch; no external scratch, test exclusions,
or change to the out-of-scope cache policy was used to manufacture a full PASS.
The load-report subprocess also cannot execute `ps` in this sandbox; the script
continues to the four real gates.

README generation exits 0 and changes no generated file. The nine Python README
tests have eight passes and one pre-existing scope-text failure: the test still
requires `standard profile is unsupported`, while current generated text already
states that standard verified states are supported after ADR 0043. Generator,
test and input reference evidence are unchanged by this patch.

The task worktree's copied `artifacts/sessions.json` has no R-SCF-1 entry and
`.delegate/R-SCF-1/report.json` is absent. Registry/dashboard files and the
README generator's result types are outside OWNED_FILES. The orchestrator must
record this report and update/upload its dashboard from the master checkout;
this implementer did not start another delegated run or modify git metadata.

## Delivered QE regression (2026-11-23 follow-up)

Supersedes the earlier absent-fixture status. Copied all 46 files (2,806,532
bytes) byte-for-byte from the main worktree; no UPF is included. Provenance,
QE version, commands and hashes remain in the unmodified fixture README,
MANIFEST and SHA256SUMS. This is density-update correctness evidence, not
Mac performance evidence or independent full-SCF trajectory agreement.

`cargo test -p dft-engine --test refscf -- --nocapture` passes all 10 tests,
including the real production Broyden replay. Initial densities recovered with
beta 0.7 and 0.3 agree to 3.46945e-18 e/bohr^3 (relative max 1.17140e-16);
beta 1 output agrees to 3.83323e-20 e/bohr^3. Each iteration uses QE's own x/F;
history occupancy is 0, 1, 2. All 11,936 Miller indices are matched in BAANDI
norm/Miller ordering. The ten lowest nonzero shells are checked separately.

| k | max density error (e/bohr^3) | RMS density error (e/bohr^3) | G0 error (e/bohr^3) | max error over lowest 10 shells (e/bohr^3) | dr2 (Ry) | dr2 minus printed (Ry) |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 3.46944695e-18 | 3.18984922e-20 | 3.46944695e-18 | 1.53329342e-19 | 0.438555331969 | 1.96919914e-09 |
| 2 | 3.46944695e-18 | 3.35435183e-20 | 3.46944695e-18 | 8.67361738e-19 | 0.019564039712 | -2.88015854e-10 |
| 3 | 8.67361738e-19 | 1.10392361e-20 | 0.00000000e+00 | 8.67361738e-19 | 0.000853988046877 | -1.95312339e-09 |

Charge differences versus QE are -3.74849e-15, -3.74849e-15, 0 electrons;
absolute total charge agrees with 32 electrons within 1.43e-14 electrons.
Maximum shell relative max/RMS errors are 6.4000e-16 / 9.6936e-16.
Regression limits: max and RMS density errors <= 1e-19 + 1e-13 times the
reference group's maximum amplitude, and charge error < 1e-12 electrons.
These give rounding headroom over the measured errors without accepting
physics-scale discrepancies. dr2 must match the printed rounded value within
half its last decimal place plus 1e-14 Ry, even though QE prints '<'.
CLI general diagnostic defaults remain unchanged; this fixture is tested at
the tighter measured limits. Native MPI/MKL Si128 remains a dev08 check.

CLI command (replace the replay TOML's UPF placeholder with the local file):

```sh
target/release/dft-cli refscf-check target/scratch/si8-qe-replay.toml \
  --qe tests/reference/qe/i25-si8/si8-np1 --steps 3 --atol 1e-19 --rtol 1e-13
```

Exit 0; no failing term. Machine-readable per-shell errors and diagnostics:
[qe-si8-validation.json](qe-si8-validation.json), [qe-si8-cli.jsonl](qe-si8-cli.jsonl).
Final full gate: exit 1, 813 passed / 3 known cache failures / 181 ignored;
build, fmt and clippy pass. The orchestrator explicitly accepts those TMPDIR
failures as an environment limitation. Fixed-point verification remains the
accepted patch-2 TODO. No commit was created.
