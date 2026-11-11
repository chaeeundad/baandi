# QE 7.6 against BAANDI after the 2026-11 adoption round (provisional)

Status 2026-11-11. PROVISIONAL: the QE numbers are the dev08 relays of `qe-vs-baandi-current.md` (unperturbed single start, QE best start `startingwfc = atomic`, strict = QE conv_thr 1e-10 with BAANDI's strict certified class, standard = conv_thr 1e-8). The BAANDI numbers are the dev08 campaign medians (ledger CAMPAIGN-1, MKL build, all adopted defaults plus the eight winners, perturbed-start protocol of +-1e-7 Bohr, seeds 601-605, load gated). The two sides were NOT measured in one session: a like-for-like run (same day, same node, same protocol for both codes) is requested from the master session. QE's SCF count does not depend on the perturbation, while BAANDI's does (the perturbed start is harder than the unperturbed one, e.g. unperturbed si32 strict needed about half the Davidson work), so the table is conservative for BAANDI.

## One core (the only like-for-like column: one process, one physical core)
| case | class | QE (s) | BAANDI now (s) | BAANDI / QE | before the round (doc `qe-vs-baandi-current.md`) |
|---|---|---|---|---|---|
| si16 | standard | 1.413 | 1.934 | 1.4x | 2.1x (pure build) |
| si16 | strict | 1.831 | 3.795 | 2.1x | 4.1x |
| si32 | standard | 5.028 | 8.438 | 1.7x | 2.8x |
| si32 | strict | 6.700 | 22.768 | 3.4x | 6.1x |
| si8 | standard, strict | 0.469 / 0.552 | not re-measured | | 2.1x / 3.3x |

## Eight cores (mixed models: QE = 8 MPI ranks, BAANDI = 8 threads; not like-for-like)
| case | class | QE (s) | BAANDI now (s) | BAANDI / QE | before |
|---|---|---|---|---|---|
| si16 | standard | 0.565 | 1.356 | 2.4x | 3.0x |
| si16 | strict | 0.688 | 2.412 | 3.5x | 6.0x |
| si32 | standard | 1.700 | 4.589 | 2.7x | 3.5x |
| si32 | strict | 1.984 | 10.662 | 5.4x | 8.1x |

## Reading
1. The round roughly halved the one-core gap: si32 strict 6.1x to 3.4x, si32 standard 2.8x to 1.7x, si16 standard 2.1x to 1.4x. QE is still faster in every cell.
2. The strict class includes BAANDI's fixed-point certificate, which QE does not do (about a third of the Davidson work on si32 before the certificate screen).
3. The remaining one-core gap is dominated by the Davidson iteration cost (H application 38 percent and orthonormalization about 20 percent of Davidson time on si32 strict, ledger MICRO-ORTHO and the stage profile), and by the per-SCF-step Davidson count (QE about 3 iterations per SCF step on these cells).
4. The 8-core gap is wider because QE scales better with MPI at this size (4.3x) than BAANDI scales with threads (2.1x); the MPI work (M1 to M2d) targets this and is still being measured.

## What to run next (request to the master session)
Same day, same node, load gated: QE best start and BAANDI (current defaults, MKL) on si8, si16, si32 in both classes, at 1 core and at 8 cores (QE 8 ranks; BAANDI 8 threads and, once measured, MPI M2), under two start protocols for BAANDI (unperturbed and perturbed +-1e-7 Bohr seeds 601-605) and QE on its single start plus the perturbed starts for completeness. Replace this table with those numbers.
