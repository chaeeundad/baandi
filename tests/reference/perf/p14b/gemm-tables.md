# P14B logical GEMM accounting

Flops use 2mnk; complex real-operation equivalents are four times these counts. Times are summed logical-call wall times in seconds (includes panel scheduling and assembly), not backend-only CPU time. Shared Mac observations are not acceptance timings.

The minimum column reports a conservative bound for the **same unlocked residual-based step**. Reduced matrices need new rows/columns only. Other categories preserve current certification and two-pass projection; their bounds do not assume unvalidated locking or adaptive QR.

## si8

| Phase | Site | Calls | 2mnk | Minimum 2mnk | Seconds | GEMM time % |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| certificate | projection_overlap_1 | 67 | 56921668 | 56921668 | 0.007914 | 0.45 |
| certificate | projection_overlap_2 | 67 | 56921668 | 56921668 | 0.007919 | 0.45 |
| certificate | projection_subtract_1 | 67 | 56921668 | 56921668 | 0.005675 | 0.32 |
| certificate | projection_subtract_2 | 67 | 56921668 | 56921668 | 0.005704 | 0.32 |
| certificate | reduced_hamiltonian | 68 | 1851509966 | 127405728 | 0.150047 | 8.44 |
| certificate | ritz_and_residual | 136 | 1000196720 | 1000196720 | 0.076524 | 4.30 |
| ordinary | projection_overlap_1 | 325 | 652235858 | 652235858 | 0.072740 | 4.09 |
| ordinary | projection_overlap_2 | 325 | 652235858 | 652235858 | 0.073004 | 4.11 |
| ordinary | projection_subtract_1 | 325 | 652235858 | 652235858 | 0.055810 | 3.14 |
| ordinary | projection_subtract_2 | 325 | 652235858 | 652235858 | 0.055575 | 3.13 |
| ordinary | reduced_hamiltonian | 340 | 10644803204 | 1511411016 | 0.859277 | 48.34 |
| ordinary | ritz_and_residual | 680 | 5303821600 | 5303821600 | 0.407512 | 22.92 |

Restart contraction: 0 extra GEMMs (saved Ritz H products). Internal block MGS: scalar dot/AXPY, 0 GEMMs. There is no separate generalized-overlap reduced matrix for this orthonormal-basis solver. Projection overlaps are listed explicitly.

## k444s

| Phase | Site | Calls | 2mnk | Minimum 2mnk | Seconds | GEMM time % |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| certificate | projection_overlap_1 | 109 | 6665072 | 6665072 | 0.001326 | 0.44 |
| certificate | projection_overlap_2 | 109 | 6665072 | 6665072 | 0.001195 | 0.39 |
| certificate | projection_subtract_1 | 109 | 6665072 | 6665072 | 0.000699 | 0.23 |
| certificate | projection_subtract_2 | 109 | 6665072 | 6665072 | 0.000676 | 0.22 |
| certificate | reduced_hamiltonian | 141 | 52941942 | 18909182 | 0.005432 | 1.79 |
| certificate | ritz_and_residual | 282 | 51184992 | 51184992 | 0.004252 | 1.40 |
| ordinary | projection_overlap_1 | 1366 | 193341010 | 193341010 | 0.025190 | 8.29 |
| ordinary | projection_overlap_2 | 1366 | 193341010 | 193341010 | 0.024397 | 8.03 |
| ordinary | projection_subtract_1 | 1366 | 193341010 | 193341010 | 0.017074 | 5.62 |
| ordinary | projection_subtract_2 | 1366 | 193341010 | 193341010 | 0.016948 | 5.58 |
| ordinary | reduced_hamiltonian | 1686 | 1500655894 | 492291096 | 0.135418 | 44.56 |
| ordinary | ritz_and_residual | 3372 | 895597536 | 895597536 | 0.071283 | 23.46 |

Restart contraction: 0 extra GEMMs (saved Ritz H products). Internal block MGS: scalar dot/AXPY, 0 GEMMs. There is no separate generalized-overlap reduced matrix for this orthonormal-basis solver. Projection overlaps are listed explicitly.

## syn20-skew

| Phase | Site | Calls | 2mnk | Minimum 2mnk | Seconds | GEMM time % |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| certificate | projection_overlap_1 | 1 | 3936 | 3936 | 0.000001 | 0.33 |
| certificate | projection_overlap_2 | 1 | 3936 | 3936 | 0.000001 | 0.27 |
| certificate | projection_subtract_1 | 1 | 3936 | 3936 | 0.000001 | 0.20 |
| certificate | projection_subtract_2 | 1 | 3936 | 3936 | 0.000000 | 0.16 |
| certificate | reduced_hamiltonian | 2 | 102664 | 55432 | 0.000011 | 3.69 |
| certificate | ritz_and_residual | 4 | 196800 | 196800 | 0.000016 | 5.55 |
| ordinary | projection_overlap_1 | 7 | 74292 | 74292 | 0.000011 | 3.74 |
| ordinary | projection_overlap_2 | 7 | 74292 | 74292 | 0.000010 | 3.25 |
| ordinary | projection_subtract_1 | 7 | 74292 | 74292 | 0.000010 | 3.24 |
| ordinary | projection_subtract_2 | 7 | 74292 | 74292 | 0.000007 | 2.54 |
| ordinary | reduced_hamiltonian | 13 | 824920 | 452230 | 0.000103 | 35.17 |
| ordinary | ritz_and_residual | 26 | 1409088 | 1409088 | 0.000123 | 41.85 |

Restart contraction: 0 extra GEMMs (saved Ritz H products). Internal block MGS: scalar dot/AXPY, 0 GEMMs. There is no separate generalized-overlap reduced matrix for this orthonormal-basis solver. Projection overlaps are listed explicitly.

## syn4-elongated

| Phase | Site | Calls | 2mnk | Minimum 2mnk | Seconds | GEMM time % |
| --- | --- | ---: | ---: | ---: | ---: | ---: |
| certificate | projection_overlap_1 | 2 | 1932 | 1932 | 0.000000 | 0.42 |
| certificate | projection_overlap_2 | 2 | 1932 | 1932 | 0.000000 | 0.39 |
| certificate | projection_subtract_1 | 2 | 1932 | 1932 | 0.000001 | 1.31 |
| certificate | projection_subtract_2 | 2 | 1932 | 1932 | 0.000000 | 0.27 |
| certificate | reduced_hamiltonian | 3 | 14868 | 7098 | 0.000002 | 1.73 |
| certificate | ritz_and_residual | 6 | 21504 | 21504 | 0.000002 | 2.08 |
| ordinary | projection_overlap_1 | 19 | 31416 | 31416 | 0.000006 | 5.78 |
| ordinary | projection_overlap_2 | 19 | 30744 | 30744 | 0.000005 | 4.93 |
| ordinary | projection_subtract_1 | 19 | 31416 | 31416 | 0.000008 | 7.63 |
| ordinary | projection_subtract_2 | 19 | 30744 | 30744 | 0.000004 | 3.58 |
| ordinary | reduced_hamiltonian | 28 | 206514 | 98952 | 0.000052 | 47.74 |
| ordinary | ritz_and_residual | 56 | 237216 | 237216 | 0.000026 | 24.13 |

Restart contraction: 0 extra GEMMs (saved Ritz H products). Internal block MGS: scalar dot/AXPY, 0 GEMMs. There is no separate generalized-overlap reduced matrix for this orthonormal-basis solver. Projection overlaps are listed explicitly.

## Initial incremental-Gram experiment (subsequently adopted)

| Case | Before 2mnk | Trial 2mnk | SCF before / trial | Energy change (Ha) | T1/T4 density identical |
| --- | ---: | ---: | --- | ---: | --- |
| si8 | 21636961594 | 9489226838 | 15 / 15 | -1.421085e-14 | True |
| k444s | 3300404692 | 2258007134 | 10 / 10 | 0.000000e+00 | True |
| syn20-skew | 2846384 | 2426462 | 6 / 6 | -2.842171e-13 | True |
| syn4-elongated | 612150 | 496818 | 9 / 9 | 1.776357e-15 | True |

The four requested cases pass the narrow iteration/energy test, with these toy iteration changes (subsequently allowed by the architect):

| Case | Before SCF | Trial SCF | Energy difference (Ha) |
| --- | ---: | ---: | ---: |
| syn-fft-primes | 14 | 16 | 0.000000e+00 |
| syn-p3 | 8 | 10 | -7.105427e-15 |
| syn-s1 | 13 | 14 | 5.329071e-15 |
| syn-s2 | 17 | 16 | -8.881784e-16 |

This table records the initial pre-adoption snapshots. The architect subsequently adopted the trial as default; see adopted/acceptance.md for the current rule and results. Exact shape histograms, including Ritz column widths, remain in each raw JSON stderr. Backend microkernel packing is not separately timed. No new QE run or cross-code performance claim is made.
