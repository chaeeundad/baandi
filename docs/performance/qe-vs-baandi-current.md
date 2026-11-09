# QE 7.6 vs BAANDI, current speed table (2026-11-09)

All wall seconds, dev08 (Xeon Gold 6230), real Gamma path, one measured start per cell (unperturbed input), same grids and UPF. QE: `K_POINTS gamma` real path, 1 core, conv_thr 1e-8 (standard class) and 1e-10 (strict class), stage-A runs. BAANDI: standard profile versus QE 1e-8, strict (certified) profile versus QE 1e-10. BAANDI numbers are the best recorded build at the stated commit; the newest changes (AUD1 and the pending MKL eigensolver) are noted separately. Sources: `tests/reference/perf/dev08_stage_decomposition_relayed.md`, `dev08_materials_qe_relayed.md`, the ledger rows FIN1, MGS-B, AUD1.

## 1 core
| case | class | QE (s) | BAANDI pure (s) | ratio | BAANDI MKL (s) | ratio | measured at |
|---|---|---|---|---|---|---|---|
| si8 | standard | 0.469 | 1.057 | 2.3x | 0.993 | 2.1x | c000d11 (FIN1) |
| si8 | strict | 0.552 | 2.097 | 3.8x | 1.993 | 3.6x | c000d11 |
| si16 | standard | 1.413 | 2.955 | 2.1x | 2.942 | 2.1x | pure c000d11, MKL e6d9b9f blocked (= 43bed40 default) |
| si16 | strict | 1.83 | 8.968 | 4.9x | 7.625 | 4.2x | same |
| si32 | standard | 5.04 | 16.194 | 3.2x | 14.054 | 2.8x | pure c000d11, MKL e6d9b9f blocked |
| si32 | strict | 6.70 | 46.846 | 7.0x | 51.043 | 7.6x | same (the MKL single start takes the long trajectory 14/464/6512; its distribution median is better, see below) |
| k444s (symmetry on, 10 k) | standard | 0.516 | 1.612 | 3.1x | - | - | 88cb469 (not re-measured) |

## 8 cores (BAANDI only; QE was measured on 1 core in this round, so no valid ratio)
| case | profile | BAANDI MKL (s) | note |
|---|---|---|---|
| si16 | strict | 4.491 | e6d9b9f blocked, before AUD1 (AUD1 -6.2 percent) |
| si16 | standard | 2.012 | before AUD1 (-11.4 percent) |
| si32 | strict | 20.637 | 43bed40 default plus AUD1 (reported post-AUD1 wall) |
| si32 | standard | about 6.0 | 8.284 before AUD1, AUD1 -27.6 percent (derived) |

## Progress since the 88cb469 re-measurement (1 core, ratio to QE)
| case | class | 88cb469 | now (best build) |
|---|---|---|---|
| si8 | standard / strict | 2.3x / 4.2x | 2.1x / 3.6x |
| si16 | standard / strict | 2.1x / 5.3x | 2.1x / 4.2x |
| si32 | standard / strict | 3.2x / 7.3x | 2.8x / 7.0x |

## Reading and limits
1. Each cell is one unperturbed start. Starts are draws of a chaotic trajectory distribution (known 1e-10 sensitivity): perturbed-start medians are the reliable measure and show larger MKL gains than the unperturbed single starts (for example si32 strict at 8 cores -21.6 percent, MGS-B R2 and the post-adoption check), but QE was not measured on the same perturbed inputs, so they are not used for ratios.
2. MKL versus pure: MKL builds win on si16 and si32 standard and si16 strict; pure is ahead on si32 strict in this single start only.
3. The remaining gap in the strict class comes from the Davidson algebra (stage decomposition results): about 64 to 86 percent of the gap; FFT is already at or below QE per call.
4. Not in this table: the stress pass (si8 forces plus stress 6.6 s against QE 0.49 s at 886a0fe-era numbers 2.8x the SCF time) and the MKL eigensolver (pending dev08 distribution evidence).
5. Missing for a full 8-core comparison: QE timings at 8 cores for the same cases.
