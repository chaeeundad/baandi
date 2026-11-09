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
