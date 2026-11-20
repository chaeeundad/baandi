# First-attempt scheduler: final Mac correctness/work evidence

Final model: geometric-amplitude tail, `VERIFY_SCHEDULE=first`. Default remains `on`. **Not approved for adoption**: correctness passes, but the supplemental Si64 case does not remove the failed attempt. No Mac timing was analyzed.

Base commit: `9bfda96c31b7172a3f9be98f73fd1994bace340d`. Raw output: `target/scratch/first-schedule/matrix/*/{parent,default,trace-tail,first-tail,strict,strict-first}`; full per-row counters/features, output/input hashes and binary hashes: [correctness JSON](first-schedule-correctness.json). The lag-only precursor is retained separately and is not the final candidate.

## Final 8-system × 3-seed × 2-threshold matrix

All 48 candidate runs converge and pass the unchanged certificate, including the Coulomb and enabled eigen-residual checks. The precursor has 48 direct default/parent and 24 strict-path switch identities. Subsequent predictor changes are confined to the opt-in branch. The final tail binary has 48 first-trace/parent scientific/density identities, plus four fresh direct default/strict identities on Si8 and synthetic spd. Cached direct-default and strict-switch rows are not mislabelled as fresh final-binary runs. The strict references also satisfy the strict occupied residual/eigenvalue-change criteria.

| System (six runs each) | Baseline first failures | Predicted failure skips | Wrong skips | Maps A→B | Attempts A→B | Total H A→B | Failed H A→B |
|---|---:|---:|---:|---|---|---|---|
| si8 | 6 | 6 | 0 | 67→72 | 12→7 | 4812→4646 | 299→53 |
| si16 | 6 | 6 | 0 | 65→66 | 12→6 | 10682→9919 | 578→0 |
| si32 | 4 | 4 | 2 | 71→75 | 10→7 | 19882→19776 | 670→192 |
| si-skew | 6 | 6 | 0 | 51→51 | 15→6 | 1655→1557 | 109→0 |
| mgo8 | 1 | 1 | 5 | 51→56 | 7→6 | 9285→9387 | 74→0 |
| syn-gaussian-fcc | 6 | 5 | 0 | 57→58 | 13→7 | 2207→2141 | 84→10 |
| syn-p3 | 0 | 0 | 6 | 33→39 | 6→6 | 81→87 | 0→0 |
| syn-spd | 0 | 0 | 3 | 24→27 | 6→6 | 72→75 | 0→0 |

Totals: maps **419→444**, attempts **81→51**, total H **48,676→47,588**, failed-attempt H **1,814→255**, failed-attempt Davidson iterations **151→16**. The 28 skipped baseline failures are counterfactual labels; 16 baseline passes would also be skipped, and one baseline failure is missed. These are not weighted flop or wall-time savings. Several controls acquire extra maps.

Energy against same-geometry strict reference: maximum **1.4630785472036223e-10 Ha** at the 5e-9 Ha standard class and **3.5704772471945034e-12 Ha** at the 1e-10 Ha tighter class. Maximum candidate/control difference **5.624656296276953e-11 Ha**. This is a trajectory change, not a rounding-only change; the requested accuracy-class gates, not a bit-equality claim for different stopping states, apply. These are internal references, not QE/absolute-accuracy proofs.

## Supplemental Si64, seed 601, standard class

| Policy | Maps | Attempt iterations | Failed H | Total H |
|---|---:|---|---:|---:|
| control | 13 | 11,13 | 395 | 6644 |
| lag-only | 13 | 12,13 | 394 | 6622 |
| geometric tail (final) | 14 | 13,14 | 394 | 6661 |

The final model skips maps 11 and 12, but its map-13 attempt still fails (7.621069462727505e-9 Ha); map 14 passes (1.575336923724383e-11 Ha). Final energy error against the strict reference is **2.3021584638627246e-12 Ha**. No certificate rule changed. Both models fail the intended Si64 attempt-count reduction locally. Omitting HALF changes the subsequent trajectory: baseline map 13 is not the same candidate state. This directly limits the proposed 2.52 s/0.54 s savings argument. No net speed gain is asserted.

## Commands and checks

```sh
python3 -I docs/performance/lanes/VERIFY/first_schedule_study.py \
  --records ../F-VERIFY-SCHEDULE/target/scratch/followup-ab \
  --out target/scratch/first-schedule/matrix \
  --binary target/scratch/first-schedule/tail-cli \
  --baseline target/scratch/first-schedule/baseline-cli --candidate-tag=-tail
cargo test --offline --release -p dft-engine --lib verification_schedule -- --include-ignored --test-threads=1
CHECK_SCCACHE=0 TMPDIR_990="$PWD/target/scratch" scripts/check.sh fast dft-engine dft-numerics dft-cli
git diff --check
```

- Release CLI builds and final matrix: exit **0**. Rebuilding the two binaries requires the base commit and current patch, respectively; the stored hashes identify the actual runs. The final trace/candidate runs use the frozen `tail-cli`, never a binary being rebuilt.
- Release scheduler tests: exit **0**, **6 passed**, including 2/3 simulated-rank certificate and bit-consistency checks; no MPI speed claim.
- Final fast gate: exit **1**, build/fmt/clippy pass; **473 passed, 3 failed, 137 ignored**. The three known CLI pseudo-cache tests reject the repository-local TMPDIR. The informational sandbox `ps` query also fails; gates continue. Logs: `target/scratch/first-schedule/check-tail.log`.
- Initial issues fixed: switch test supplied str instead of OsStr; clippy rejected an unnecessary unwrap; the study initially assumed a strict result had the standard-only `certified` JSON key. The first-cross horizon was hardened against rebound/re-crossing. Final checks above include these fixes.
- The verification-function source region is byte-identical to the base (hash in the JSON). No new thresholds, certificate predicates or solver tolerance changes.
- Pending: actual dev08 MPI/one-core A/A and A/B, real MPI rank evidence, hardware-specific net cost. The [request](dev08-request.md) is conditional, not a validated speedup claim.
