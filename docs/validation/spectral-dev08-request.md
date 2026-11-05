# P15C dev08 measurement request

No remote changes or runs were performed by Codex. The orchestrator should run
this after reviewing/committing the uncommitted branch, and record that commit.
Purpose: measure the **final-only explicit option**, not the earlier source-only
candidates. Default remains off because the measured H-work target is missed.

## Builds and provenance

Use the same reviewed source and existing dependencies/UPFs; no installs.
Unset all P15C diagnostic build variables (`P15C_ORACLE`, `P15C_DENSE_SCF`,
`P15C_CANDIDATES`, `P15C_RICH_START`, `P15C_FULL_GUARD`) and profiling flags.
Use separate target directories, pure Rust, identical toolchain and release
profile. Build once with ordinary flags, once with `-C target-cpu=native`:

```sh
cargo build --offline --release -p dft-cli --target-dir /tmp/p15c-default-target
RUSTFLAGS='-C target-cpu=native' cargo build --offline --release -p dft-cli --target-dir /tmp/p15c-native-target
```

Record commit, source hash manifest, rustc version, build commands/exit codes,
CPU model, affinity, RUSTFLAGS, binary SHA-256, hostname, kernel and load.
Existing QE timing is not needed: this request compares guard on/off at the
same source. Do not call these wall times a QE comparison.

## Paired cases and core budgets

Cases: `si8,k444s,syn20-skew`. Use the existing readable Si UPF via DFT_SI_UPF;
the synthetic file is generated temporarily by the existing fixture generator.
Record the UPF hashes. Budgets T=1,2,4,8 on one socket, equal affinity for both
policies; no oversubscription. Use an idle node, 5 warmup pairs, then at least
30 alternating off/final pairs per case, build and thread budget. The runner
alternates policies within each repetition:

```sh
taskset -c 0-7 python3 -I tests/reference/spectral/final_policy.py \
  --binary /tmp/p15c-default-target/release/dft-cli --upf "$DFT_SI_UPF" \
  --cases si8,k444s,syn20-skew --threads 8 --repeats 30 \
  --output /tmp/p15c-default-t8
```

Repeat with the native executable and each budget/affinity. Use separate output
paths for warmups and measured batches. The runner stores exact commands,
inputs, hashes, complete result JSON, exit codes, wall samples and audit events.
Preserve files on dev08 and supply SHA-256 with the relay.

## Required returned numbers

For every case/build/thread/policy:

- Total ordinary H columns, fixed-point-certificate H columns, additional guard
  H columns, and their sum. Count **vector columns**, not apply-block calls.
- Guard attempts per k/spin, iterations of attempts, cumulative retries, method,
  passed/failure status and any probe error. On a complete run expect one event
  per block at the final iteration, no ordinary-iteration guard work.
- Whole-run and guard-only wall samples, median paired overhead, spread and
  outliers; repeat any batch with >5% spread. Report both absolute seconds and
  relative overhead. Do not use Mac speedup columns as an acceptance baseline.
- SCF iterations, convergence, energy (Ha), density SHA-256, spectrum and residual
  comparisons; guard on/off must preserve complete real-case results within
  1e-12 Ha and identical density bytes. Check T determinism within each build;
  report cross-build differences separately.

Expected portable-build work here: Si8 off 2265 columns, guard +90; k444s off
11181, guard +2618. Native rounding may change counts: report it, do not force
agreement. Skew off is a known **wrong** 6-step certified state; guard on finds a
missing set at step 6, adds 130 guard columns, and exhausts 100 SCF steps here.
Its off/on wall ratio is **not an equal-accuracy speed comparison**. Report
recovery work and trajectory separately. No tuning to restore its false success.

Adoption decision remains off unless a subsequent, independently validated
policy meets the <3% total-work target on the well-behaved systems. These timing
runs alone cannot turn the current seeded search into a rigorous inertia proof.
