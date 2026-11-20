# PACKED_REDUCE Mac correctness evidence (2026-11-20)

No Mac performance measurements or conclusions. Default-off experimental
communication change; dev08 real MPI/MKL acceptance remains pending.

* Unit release tests:2 passed. Independent index-selection oracle covers
  orders0/1/2/17/132/264, absent/present S and H borders, empty/full norm tails.
  Reduced S/H and norm values compare by float bits OFF/ON at ranks1/2/3/4/8
  and WorkerPool threads1/3. Existing cross-rank upper-border oracle retains
  its1e-12 tolerance; no stronger cross-rank floating-bit guarantee is claimed.
* Full standard SCF: synthetic nonlocal PP, cubic and skewed cells,6 bands,
  simulated MPI2/3,1 thread/rank. Four OFF/ON pairs, all verified;
  energy bits, density bits, scientific history and verification match.
  Cubic work(SCF records,Davidson iterations,H applications,verification
  attempts)=(17,62,269,1); skewed=(11,43,197,1), both rank counts.
  Fixed-rank energy bits also agree across ranks. abs(dE)=0 Ha in every pair.
* Counter transport only: cubic3 ranks,116 border calls/rank in both arms;
  sent65520->43840 B/rank, prediction43840, saving21680 B/rank, all3 ranks.
  `packed-counter-evidence.json` stores these integer counters. The lane
  parser accepts both exported Mac trace envelopes with no dropped records.
  These small synthetic shapes are correctness evidence, not Si/dev08 speed.
* Fast gate including CLI: build/fmt/clippy PASS;362 tests passed,3 failed,
  87 ignored. Failures are the known pseudopotential-cache tests rejecting
  the required in-repository TMPDIR. No new scientific failure.
* Python tools:4 tests passed; packed-job.sh syntax passes; git diff check
  clean. The final core-only gate is recorded in packed-report.json.

Commands (scratch on990 volume):

```sh
TMPDIR="$PWD/target/scratch" cargo test --release -p dft-solvers --lib packed_border -- --nocapture
# Reproduce the full SCF oracle without leaving an engine-owned test edit:
test ! -e crates/dft-engine/tests/lane_packed_oracle.rs
cp docs/performance/lanes/DAVIDSON/packed-scf-oracle.rs crates/dft-engine/tests/lane_packed_oracle.rs
TMPDIR="$PWD/target/scratch" PACKED_ORACLE_OUT="$PWD/target/scratch" cargo test --release -p dft-engine --test lane_packed_oracle -- --nocapture
rm crates/dft-engine/tests/lane_packed_oracle.rs
python3 -I docs/performance/lanes/DAVIDSON/packed-counter.py target/scratch/packed-0-np3-false.json target/scratch/packed-0-np3-true.json
python3 -I scripts/lane_davidson_parse.py target/scratch/packed-0-np3-false.json target/scratch/packed-0-np3-true.json --json-out target/scratch/packed-tables.json
TMPDIR="$PWD/target/scratch" python3 -I docs/performance/lanes/DAVIDSON/test_tools.py
CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh fast dft-solvers dft-numerics dft-cli
CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh fast dft-solvers dft-numerics
```

The initial temporary SCF assertion compared Debug records including wall
metadata, and failed on those fields; exclude only PhaseTimings and wall
fields when comparing scientific records. The corrected oracle passed. A
trace-export compile error (shadowed communicator) was fixed before the final
passing run. Neither was a numerical mismatch. Permanent solver edits are
in distributed_davidson's existing workspace/lean.rs implementation, not a
new engine path. Temporary engine test file was removed after validation.

Limits: real MPI association/transport thresholds and MKL are checked on
dev08. Mac full SCF uses small synthetic cubic/skewed systems, not a new
real-PP/slab qualification. No changed tolerance, physics, GEMM shape, EVD
size or restart policy. Serial production solver does not consult this switch.
Both rejected switches remain off. No default adoption; no commit.
