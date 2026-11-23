# I25-STOP-TO-VERIFIED: Si128 and Si64 pure MPI8

Status: request prepared; no dev08 execution or speed claim by this implementer.
CLI dispatch, timestamp serialization and rank-rule clock normalization are
integrated under the orchestrator's 2026-11-24 ownership extension.
Operator: master session via `baandi-dev08`. Use the exact reviewed, pushed hash;
leave implementation changes uncommitted until review. Follow
[gate-on-dev08](../../requests/gate-on-dev08.md), including its 2026-11-23 update:
artifacts under `/scr`, gate on socket 1 via `/scr/mino/gate/gate.sh`, measurements
on socket 0, no measurement beside a gate. All jobs have a hard 600 s limit.

## Separate jobs

- `final-a0-gate-<hash>`: `scripts/check.sh full` through the gate script, with
  pinned toolchain and persistent warm target. Return every gate and exit code.
- `final-a0-stop-to-verified-si8-<hash>`: release correctness test below, external
  UPF, no timings used as performance evidence. One fixture, both diagnostics.
- `final-a0-stop-to-verified-si64-mpi8-<hash>`: one production-level Si64 case.
- `final-a0-stop-to-verified-si128-mpi8-<hash>`: one production-level Si128 case.

Build once before timing: `cargo build --release -p dft-cli --features mpi`
with the approved native build flags, **no MKL**. Record flags, dependency lock,
compiler, revision and binary SHA-256. Use that exact executable (same hash) for
both Si64 and Si128; do not rebuild between cases. Compilation/preflight and
queue delay are outside measured run wall. Split a cold build into its own job;
never hide timeout by bundling cases or silently truncating SCF iterations.

Release Si8 test (the full test is intentionally ignored on Mac):

```sh
DFT_SI_UPF=/absolute/external/Si.pz-vbc.UPF cargo test --release \
  -p dft-engine --test refscf stop_to_verified_si8_displaced \
  -- --ignored --nocapture --test-threads=1
```

It checks first stop 10, failed attempts 10/11, pass 12, three attempts,
nonzero recovery and identical delivered arrays/energies across diagnostic
levels. Return complete output and exit code, including any unexpected branch.
Do not relax these assertions to obtain a passing report.

## Inputs and production run

Use the planned Si128 pure-MPI8 input and its Si64 control from the existing
RSCF run package. Archive their exact TOML, cell/coordinates, band counts,
30 Ry cutoff and specified FFT grid, atomic start, LDA-PZ, Gamma, UPF SHA-256,
reference Broyden settings and unmodified standard 5e-9 Ha threshold. Confirm
`scf.driver="reference"`, `scf.threads=1`, `scf.verify="fixed_point"` (or omit it to
retain the default fixed-point verification) and `[profile] davidson_stages=false`.
The environment must not override verification or enable full diagnostics.
Do not substitute a different geometry when reusing the planned inputs.

Use 8 ranks on 8 bound physical cores on socket 0, one thread/rank;
`OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1`. Record actual MPI
launcher version, complete command, CPU binding and environment. With the
queue's established socket-0 binding, the command has this shape:

```sh
/usr/bin/time -p -o "$OUT/external-wall.txt" \
  timeout 590s mpirun -np 8 <approved socket-0 binding> \
  "$BINARY" run "$INPUT" --output-dir "$OUT" \
  >"$OUT/stdout" 2>"$OUT/stderr"
rc=$?
"$BINARY" refscf-check --stop-to-verified "$OUT/results.json" \
  >"$OUT/stop-to-verified.json" 2>"$OUT/reader.stderr"
reader_rc=$?
```

Create a fresh `$OUT` before timing. The external timer encloses launch through
MPI finalize/process exit, including all normal JSON, density, diagnostic and
stdout/stderr output; discard no output or failed attempt. The offline reader
runs afterward and is not part of delivered SCF wall. A timeout, nonzero run
exit, missing completed results or failed reader is reported explicitly; never
turn missing output into a zero interval. A completed exhausted run can produce
an `UNVERIFIED` reader report even if the SCF command exits nonzero. Retain it.

## Return record and gates

For **each case separately**, return binary/input/UPF hashes, job ID, command,
run and reader exit codes, complete artifacts, external wall seconds, diagnostic
level, verification mode, first stop map, passing map or `UNVERIFIED`, attempt
count, interval endpoints, interval wall, all three parts and closure error.
The reader's clock is rank 0's local driver clock. Never subtract clocks from
different ranks, sum rank times, or combine independently maximized stages.
Also retain all fixed-point errors, thresholds and verdicts, final energy,
charge and convergence state. Unverified output is not a delivered success.

Require parts to close and first-pass recovery to be exactly zero. Verification
remains the unchanged independent fixed-point density check. No energy-error
proxy, verification-off timing, threshold change or subtraction of failed work
is permitted. Full/default scientific identity and the complete full gate are
required before acceptance. Pre-change default and reference-run fingerprints
must exclude only clocks and explicitly optional diagnostic payloads; record
any mismatch rather than weakening the comparison.

The I26-2 diagnostic on/off delta is reported in its own table, using its
paired same-binary protocol and separate <=600 s one-case jobs if requested.
It may overlap both verification and ordinary recovery and is **not** a fourth
interval part or a correction to any part. This request needs no QE run and
makes no cross-model or cross-rank speed comparison.

## Local implementation checks (2026-11-24, integration complete)

- `cargo check --workspace --tests`: exit 0, no warnings.
- `cargo test -p dft-engine --test refscf`: exit 0; 11 passed, one ignored
  (the real displaced Si8 release run above). The live synthetic first-pass
  control has one attempt, first/passing map 1, exactly 0 s recovery, and exact
  scientific equality between full and lightweight diagnostics. Its intentionally
  large stop threshold is a control-flow fixture, not production evidence.
- `cargo test -p dft-cli stop_to_verified`: exit 0; six reader/CLI tests passed.
  Independent synthetic clock fixture: maps 10/11 fail, map 12 passes,
  interval 23 s = 6 s attempts + 12 s recovery + 5 s unassigned, closure 0 s.
  First-pass clock fixture: 3 s = 2 s attempts + 0 s recovery + 1 s unassigned.
  These invented clocks test arithmetic; they are not Si8 timing measurements.
- The actual CLI integration test runs a small synthetic SCF at default
  lightweight and explicit full diagnostic levels. Both outputs contain every
  timestamp in the canonical iteration records and history aliases, pass the
  offline reader, and produce the expected JSON via `--stop-to-verified`.
  The original replay mode still requires `--qe`. Six tests finish in 0.53 s
  on this Mac; this is a test-runtime observation, not a speed comparison.
- `cargo fmt --check && cargo clippy -p dft-engine -p dft-cli --all-targets -- -D warnings`:
  exit 0; no warning suppression added to production code.

Earlier engine compilation and unused-reader clippy errors were corrected
before these final checks. Integration initially exceeded the existing JSON
macro recursion limit; the new timestamp fields are now inserted after the
existing record object, without changing crate-wide compiler settings.
No full Mac gate, release matrix, dev08/QE run, baseline-default replay,
README regeneration or dashboard upload was performed. README/dashboard files
remain outside the ticket's owned paths. The rank-rule fingerprint now strips
all new clocks while preserving the presence/absence of verification timestamps.
The complete rank matrix remains the orchestrator's release validation.
