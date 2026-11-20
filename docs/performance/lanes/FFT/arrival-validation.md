# Arrival trace: Mac correctness validation

No Mac profiling supports a performance conclusion. This is instrumentation
and identity validation only; dev08 has not been run from this session.

## Actual scientific outputs

Same binary/input, OFF/count-only/timed (COLLECTIVE_TRACE off/on/on,
STAGE_TIMERS off/off/on), PHASE_PROFILE=on, FFT_PROFILE=on,
DAVIDSON_STAGES=on, PHASE_WAIT=off, FFT_KEPT_FUSE=off.

- Serial CLI Si8: -31.321776563202498 Ha, verified in all three arms,
  energy/eigenvalue/occupation/scientific history bits equal; density bytes
  equal, SHA-256 43f2a33d8b0d0a5a90129076353d8e4ad26b2305695c3d834c13b7598fb7342b.
  Serial has no SharedComm events; this is the unchanged default-path control.
- Real-UPF Si8 simulated MPI SCF, ranks 2/3/8: all three arms have identical
  energy bits, density FNV and SCF/Davidson/H counts. Energies respectively
  -31.321784923608313 / -31.32178492360831 / -31.321784923656523 Ha; within each
  rank count delta=0 Ha. Counts respectively 11/52/437, 11/52/437, 10/50/404.
  The existing harness exports density FNV, not complete density/history files;
  do not claim those MPI bytes were directly compared.
- Real-UPF Si8 fixed-H warm solve, simulated MPI3: eigenvalue FNV
  0xfdff2330905ef17e, all eigenvalues/residuals/work counters equal in all arms.
  Actual CLI profile serializer produces 408 events/rank, including 27
  allgatherv, 130 alltoall, 251 tree_sum. All count/timed signatures, payloads,
  stage labels/work/shapes match exactly. The new trace adds no SCF collective.
  Diagnostic gather transport is independently exercised in the wrapper test.

The parser succeeds on both actual count and timed outputs (408 fine rows).
mac-arrival-trace.json is the exact trace object extracted from the timed
DAVBENCH_JSON record, not a fabricated fixture; its diagnostic times must not
be used for speed claims. Reported clock uncertainty in that retained run is
1001 ns; no production rank-wait conclusion is made from simulated ranks.
Some semantic sites/preceding stages are unattributed; the trace does not
invent their identities. arrival-equivalence.json records numerical evidence;
arrival-provenance.json identifies delivery sources and binaries.

## Commands executed

Build CLI and fixed-H harness, exit 0:

~~~sh
cargo build --offline --release -p dft-cli --bin dft-cli --example davidson_bench
~~~

Using target/scratch/fft-si8.toml, the retained Si8 Phase-0 input with external
Si UPF, for each ARM=off,counts,timed and flags above (all six executions exit 0):

~~~sh
target/release/dft-cli run target/scratch/fft-si8.toml \
  --output-dir "target/scratch/arrival-cli-$ARM"
target/release/examples/davidson_bench --input target/scratch/fft-si8.toml \
  --sim-ranks 3 --kinds warm --repeats 1
~~~

SCF harness build, exit 0:

~~~sh
cargo test --offline --release -p dft-engine --test m2_entry_cost \
  --no-run --message-format=json
~~~

Executable selected from compiler-artifact JSON; each of the three flag arms
exits 0 with DFT_SI_UPF pointing at the external file:

~~~sh
M2_PROFILE_CASE=si8 M2_PROFILE_RANKS=2,3,8 M2_PROFILE_SERIAL=off \
M2_PROFILE_MODE=production M2_PROFILE_LADDER=on \
"$FFT_SCF_TEST_BINARY" --ignored --nocapture --exact m2_entry_cost_profile
~~~

Verification, all exit 0:

~~~sh
python3 -I docs/performance/lanes/FFT/verify_arrival.py target/scratch
python3 -I docs/performance/lanes/FFT/wait_parse.py docs/performance/lanes/FFT/mac-arrival-trace.json
python3 -I docs/performance/lanes/FFT/test_wait_parse.py
~~~

Five parser tests: legacy schema contract/ties/invalid data, actual production
trace and invalid clock/count-only cases, production-format known last rank.
Wrapper unit test exercises all five collective kinds OFF/count/timed on three
simulated ranks, unchanged values, post-SCF gather and exact u64 limb transport
above 2^53. Numerics unit test checks disabled empty state, fixed-buffer pointer,
count-only zero clocks, and overflow without growth. Both pass in the final gate.

## Required gate

~~~sh
CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh fast \
  dft-numerics dft-solvers dft-physics dft-engine dft-cli
~~~

Final exit 1: **build/fmt/clippy PASS**, 705 tests passed, 178 ignored, three
known in-repository TMPDIR failures:

- pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites
- pseudo::cache_resolution_and_no_network_side_effects
- pseudo::cached_fetch_checks_hash_and_provenance

Earlier iterations found clippy requirements for constant chunks and item order;
those were fixed before the final gate. Sandbox blocks the optional ps load
report, but the script runs all gates. Final log:
target/scratch/arrival-delivery-check.log. No native MPI/MKL validation or full
push gate is claimed. Changes remain uncommitted.
