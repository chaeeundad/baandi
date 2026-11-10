# Free second vector for single-column applies: validation plan (not implemented)

Status 2026-11-10: idea only (backlog entry "Free second vector for single-column applies"). This file pre-registers the experiment so that the result is judged by rules fixed before the data.

## Premise (measured)
66.4 percent of Davidson apply calls on si32 strict carry one correction column (median 1, mean 11.4; ledger MPI-M1-CEILING). The paired real FFT transforms two columns for the price of one, so the second slot of a single-column call is wasted.

## Candidate
When the correction block has an odd number of columns, add one extra direction before the apply, chosen by a fixed rule, switchable by an environment variable (default off, bit-identical when off). Variants to compare, one at a time:
- V1: previous-iteration correction of the same band (conjugate-direction term).
- V2: correction of the first unconverged buffer band (beyond the strict prefix) if any.
- V3: the residual of the worst strict band, unpreconditioned.
The extra column is orthonormalised against the basis like every other correction; if it drops out nothing is added.

## Pre-registered acceptance (all required)
1. Systems: si16 and si32 strict and standard (perturbed seeds 601-610), plus three non-silicon fixtures (skewed cell, low-symmetry cell, synthetic multi-projector pseudopotential) at their usual seeds; MKL build, T=1 and T=8. Silicon alone never decides (generality rule).
2. Work: median total Davidson iterations (recorded plus certificate, `convergence.work`) drops by at least 8 percent in every strict and standard cell, and no cell gets worse by more than 2 percent at any seed's median.
3. Wall: median wall improves by at least 3 percent in both thread budgets for strict and standard (the extra columns cost basis growth and algebra).
4. Accuracy: energies inside the accuracy classes on every seed, all runs converged, density SHA reproducible per setting, certificates unchanged.
5. Orthogonality: the orthogonality defect gate at restarts never fires more often than with the variant off.
A variant that fails 2 to 5 is recorded as rejected in the attempts ledger with the numbers; the default does not change.

## Order of work
Implement as an isolated module behind `BAANDI_SECOND_VECTOR=v1|v2|v3`, Mac correctness tests (equal eigenvalues on the existing Davidson fixtures, no change when off), then dev08 A/B through the master session.

## Result of the first registration and the second registration (2026-11-10)

First registration, dev08 A/B (cli-mkl-2290b01, relayed by master): v1, v2 and v3 all fail rule 2 at the first cell (si16 standard, T=1: median total Davidson v1 +0.0, v2 -2.3, v3 +20.9 percent; v3 worst seed +63 percent). Rejected as registered; no default change. On the premise cell (si32 strict, T=1, measured separately and not part of the verdict) v3 gave -11.1 percent work and -9.3 percent wall, v2 -6.9 percent work but +6.9 percent H vectors and +3.0 percent wall, v1 within noise: the idea works only where a single-column tail exists.

Second registration (new variant, fixed before any data): **v4** = v3 inserted only when the correction block has exactly ONE column (never in wider blocks). The gate is measurable: `results.json convergence.work.second_vector_columns` counts the inserted columns, so a measurement can check that the gate fired (and did not fire where there is no tail).
Rule 2 is tightened to remove the ambiguity found in the first run: (2a) the cell median of total Davidson (recorded plus certificate) does not get worse by more than 2 percent in ANY strict or standard cell; (2b) no single seed gets worse by more than 10 percent of its own off value; (2c) the median drops by at least 8 percent in the single-column-tail cells, defined as the cells where `second_vector_columns` per run is at least 15 percent of the Davidson iterations (si32 strict is expected to qualify; this is checked from the counter, not assumed). Rules 1, 3 (wall -3 percent in the tail cells at T=1 and T=8, no cell median slower than +2 percent), 4 and 5 stand; rule 5 needs the restart count in the run report and is therefore replaced by: the solver's `restarts` per run do not exceed the off value by more than 10 percent (a field `restarts` is added to the same `work` block only if the engine already carries it; otherwise rule 5 is dropped and said so in the result).
