# Serial density fusion: implementation and Mac correctness

2026-11-20, LANE-GLUE. Ready for dev08 A/B, default OFF. No adoption or closure.
Source before edits: clean `311887c`; changes are uncommitted. The explicit
follow-up authorizes serial band-density/FFT plumbing beyond the original
Phase-0 grid/mixer files. `parallel.rs` exposes communicator presence (including
rank one); `dft-mpi::MPI_ENABLED` exposes the backend feature without changing
MPI operations. This also excludes MPI rank-one CLI runs that do not register
a communicator. No manifest or unsafe code change.

## Arithmetic and traffic

The new FFT sink computes each split band's `+0 + weight * square / volume`,
then adds band a and band b separately to rho. It does not form a pair sum
before addition/division. The inverse FFT and REUSE_PSI writes are identical.
Unpaired occupied tails use existing contribution scratch and are folded at
their original band position. Zero occupations are skipped as before. Strict
maps restore all per-band arrays; a validity assertion prevents a future
certificate caller from silently reading discarded contributions.

The counter model deliberately conservatively charges old split operations
32 bytes/band-point versus new rho updates 16; it does not credit old output
reduction traffic. It is not a hardware memory counter. Counters are gated by
GLUE_PROFILE and PHASE_PROFILE, cached outside the loop, with clocks disabled
for these runs. Successful calls use actual grid/occupation sizes, including
coarse and verification calls. No Si-specific branch or constant is in code.

| fixture, standard one core | eligible/fused calls B | A logical bytes | B logical bytes | saved bytes |
|---|---:|---:|---:|---:|
| Si8 cubic, displaced | 18 | 350542336 | 175271168 | 175271168 |
| Si skew | 22 | 30211456 | 15105728 | 15105728 |
| synthetic S2 | 36 | 288000 | 144000 | 144000 |
| Si32 | 18 | 5249359872 | 2624679936 | 2624679936 |

Every row satisfies A−B = 16×paired_band_points exactly. In the independent
Si64 15-map model: 716636160 paired band-points and 11466178560 bytes saved,
priced at 0.884054 s with the relayed 12.97 GB/s STREAM. Actual Si64 counters
and wall gain remain dev08 work; no Mac performance inference is made.

## Scientific results

[Final evidence](fusion-mac-equivalence.json): pre-edit release binary versus
candidate OFF versus ON, all 12 runs exit 0, converge and verify. Final energy,
all energy terms, density.bin SHA-256, eigenvalues, occupations, SCF history,
band accuracy, work and verification data agree bit for bit (timings excluded).
Maximum |dE|=**0 Ha**. Inputs and PP hashes are embedded; UPF contents are not.

- Si8: −31.32177656320248 Ha.
- Skew Si: −7.2506046203814964 Ha.
- Synthetic S2: 0.6354522013776593 Ha.
- Si32: −126.70343204998849 Ha.

[Strict skew reference](fusion-mac-strict.json): 3 arms exit 0, converge,
−7.250604620384125 Ha, |dE|=0 Ha and identical certificate/band-accuracy data;
all fusion counters zero. The strict output schema has no `verified` field;
no standard verification flag is fabricated. No trajectory or acceptance
rule changes occurred, so no new accuracy calibration is claimed.

[Thread fallback](fusion-mac-t3.json): Si8 with three threads OFF/ON has
identical scientific bits, verified, and zero fusion counters. This checks
that the switch does not change the multi-worker lane; it is not a speed test.

[Shape stress](fusion-mac-shapes.json): elongated and slab runs exhaust as
before (exit 3 in all six arms), with identical density, energy and history.
Their energy values are −14.48485038813729 and −15.222453803896748 Ha. These
are equivalence checks on known failing fixtures, not convergence evidence.
Saved logical bytes: 147294720 and 194918400 respectively, exact model matches.

Candidate binary SHA-256:
`6498893ab3734569b42403ac1286bbffa827501e8f26d49d9144730fffdd928f`.
Pre-edit binary SHA-256:
`6ca4aeacffc2efe750dd5c59e1d20fc8061559e9a69d515e7c0e423b9e7c5221`.

## Gates actually run

- `cargo build --release -p dft-cli`: exit 0.
- Four `fusion_verify.py` invocations recorded in the evidence: exit 0 each;
  standard fixtures (12 arms), strict (3), thread fallback (2), shape stress (6).
  Stage clocks disabled, no timings analyzed or retained as performance data.
- `TMPDIR_990="$PWD/target/scratch" CHECK_SCCACHE=0 scripts/check.sh fast
  dft-mpi dft-numerics dft-engine dft-cli dft-physics dft-solvers`: exit 1;
  build/fmt/clippy PASS, 728 passed / 3 failed / 179 ignored in 74 binaries.
  Failures are the previously documented in-repo TMPDIR CLI pseudo-cache
  tests. The sandbox also rejects the script's `ps` load probe; gates continue.
  Initial clippy issues were fixed before this final run.
- `BAANDI_SERIAL_DENSITY_FUSED=on cargo test -p dft-engine --test m2_owned`:
  exit 0, 3 passed. Existing simulated-MPI owned/gathered density and standard
  accuracy assertions still pass with the new switch requested.
- New unit coverage in the fast gate: pure ordered and generic full-grid FFT
  sinks retain bits/kept grids, varied shapes and weights, validation failures;
  engine odd bands, changing occupations, zeros, strict restoration; dispatch
  rejects registered simulated ranks 1/2/3/8, strict and multi-worker maps.
  Real-MPI/MKL builds are not available on this Mac; dev08 must validate them.
- `python3 -I docs/performance/lanes/GLUE/test_fusion_tools.py`: exit 0,
  3 synthetic statistics/model tests. `test_parser.py`: exit 0, 3 tests.
- `python3 -I scripts/tests/test_readme_figures.py`: exit 1, 8 passed / 1
  pre-existing structural wording failure (`standard profile is unsupported`).
  Generated scientific references were not changed; no unmeasured speed result
  was added to README. The unrelated structural wording is left to its owner.
- Stored-output audit: scientific convergence work and coarse histories match
  bitwise across all 23 outputs; wall diagnostics excluded. The final verifier
  was also exercised again on synthetic S2 (three arms), exit 0.
- `git diff --check`: exit 0.

Detailed gate logs remain under `target/scratch/fusion-*.log`. This is Mac
correctness against the engine's baseline, not a new QE run. No remote change,
new PP, default adoption or commit. Next action belongs to the orchestrator:
[dev08-request.md](dev08-request.md), A/A then crossed A/B on Si64/Si32 one core.
