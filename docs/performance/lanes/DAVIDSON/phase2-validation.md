# Phase 2 local validation

Mac is used only for compilation, correctness, independent-oracle checks and
trace/parser schema verification. No Mac timing was inspected for a performance
conclusion. Stored evidence: phase2-correctness.json. Raw logs/results stay in
`target/scratch/phase2`; no UPF committed. Defaults remain unchanged.

Implementation: `upper_border_product` in real_davidson.rs calls the original
rectangle unchanged when DAVIDSON_UPPER_BORDER is off. When on, fixed 64-column
GEMMs omit rows below the last column of each tile; existing callers mirror
only the defined upper entries. It preserves the PW reduction axis, allocation
workspaces, collective grouping/payload, H work, tolerances and cluster handling.
Distributed lean border wiring necessarily edits its child module
`distributed_davidson/workspace/lean.rs`; a change only in the top-level
original lane ownership file would not reach this production path. Scripts
and lane docs provide the requested independent instrumentation/measurement.
No adapter, dependency, other lane or physics algorithm changed.

## Commands and actual outcomes

Commands used `PATH=/opt/homebrew/opt/rustup/bin:$PATH`, scratch
`TMPDIR=$PWD/target/scratch`, and `CHECK_SCCACHE=0` for the gate.

* `cargo build --release -p dft-cli`: exit 0 (build-final.log).
* `cargo test -p dft-solvers --lib upper_border_tiles`: exit 0, 1 new test.
* `BAANDI_DAVIDSON_UPPER_BORDER=on cargo test --release -p dft-solvers
  --test lean_davidson --test distributed_lean_algebra`: exit 0, 13 tests
  (solvers-on-final.log). Dense oracles, synthetic multi-projector/skew and
  clustered/degenerate operators, fallbacks and pool/rank independence.
* `CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh fast
  dft-solvers dft-numerics dft-engine dft-cli`: exit 1; build/fmt/clippy PASS,
  615 passed, 3 failed, 160 ignored (check-final.log). New direct-border test
  covers m=137, n=257, prefixes 0/17/64/136/137 and workers 1/3. The new
  distributed-border test covers simulated ranks 1/2/3/4/8, retained entries,
  residual norm tail and actual multi-tile borders. Both passed.
* The three failures are the known CLI cache cases rejecting local TMPDIR:
  `pseudo::cache_resolution_and_no_network_side_effects`,
  `pseudo::cached_fetch_checks_hash_and_provenance`, and
  `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`.
  No skips or looser assertions were introduced. Optional load-report `ps`
  is sandbox-denied; gates still execute. No permission bypass.
* `python3 -I scripts/tests/test_lane_davidson_parse.py`: exit 0, 2 tests.
* `python3 -I docs/performance/lanes/DAVIDSON/test_tools.py`: exit 0, 2
  synthetic counter/statistics tests. Their fabricated seconds are not data.
* `bash -n docs/performance/lanes/DAVIDSON/phase2-job.sh` and
  `git diff --check`: exit 0. dev08 commands have not been executed here.

## Complete SCF equivalence

Inputs generated with `python3 -I scripts/fast_screen_dav.py --case CASE
--upf /Volumes/990-2TB/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --emit-input INPUT`;
seed 601, displacement 0.001, T=1. CLI command in every arm:
`env BAANDI_DAVIDSON_UPPER_BORDER=off|on BAANDI_DAVIDSON_PROFILE=off
BAANDI_PHASE_PROFILE=off BAANDI_STAGE_TIMERS=off BAANDI_DAVIDSON_STAGES=off
BINARY run INPUT --output-dir OUTPUT`.
Cases/paths/energies/hashes and per-process exit codes are in the JSON.

Si8/Si32 original binary, final OFF, final ON: bit-identical energy,
eigenvalues, occupations, scientific history, verification and density bytes.
The serial skew case also converges and is bit-identical OFF/ON. Si32 ON at
workers 1 and 3 is bit-identical; no threading performance work was done.
Si8 strict (`BAANDI_REFERENCE_STRICT=on`) also converges with bit-identical
history including fixed-point certificates and density bytes, abs(dE)=0 Ha.
The original Phase 0 verifier requires the standard-only `verified` key and
therefore initially returned exit 1 on strict output; an explicit comparison
of the strict schema's history/certificates and other science then passed.
This was a comparison-tool schema limitation, not a solver discrepancy.

The slab fixture returns exit 3 (100-step non-convergence) in BOTH arms;
its scientific history and final density are identical. This is a coverage
limitation, not a passed convergence gate. Tolerances were not changed.

Additional Si32 OFF/ON runs enable only DAVIDSON_PROFILE to check the extended
schema. The parser and border-counter both exit 0; these runs remain
bit-identical to the non-profiled output. Shape arithmetic validates the
counter: 174 parent products in both arms, 258 tiles ON; logical
22,772,337,784 FLOP, predicted/executed ON 22,390,885,752 FLOP. These counts
are Mac fixture/schema evidence, NOT dev08 performance estimates or floors.
No timestamp from these files is used in the diagnosis.

## Remaining prerequisites

Original dev08 n/m/k/call/payload traces were not supplied with the relay.
The numeric roofline diagnosis is incomplete and explicitly marked in
`diagnosis.md`; no historical Mac trajectory is substituted. Real-MPI/MKL
rounding, performance, larger band SCF rank independence and converged slab
generality remain dev08/master gates. A/B request is ready, default stays OFF.
No README reference data changed: no new accepted scientific/performance
claim. Dashboard and final adoption belong to the orchestrator.
