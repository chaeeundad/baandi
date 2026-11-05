# P14A work accounting

Vector applications count bands passed to H. Block calls count calls of H on an entire block; a block can shrink during Davidson. Steps count Rayleigh–Ritz iterations. None of these units is interchangeable.

Timings are shared-host observations, not accepted speedups. The source, binary hashes, inputs and units are recorded in the JSON files.

## Loosened schedule experiment

Only the density-residual term is multiplied. Initial cap 1e-2 Ha and floor 1e-9 Ha remain unchanged. Final fixed-point certificate is unchanged. Steps exclude the separate fixed-point certificate; H totals include it.

| Case | Factor | Converged | SCF count | Davidson steps | Total H vectors | Energy difference from factor 1 (Ha) | Same SCF count |
| --- | ---: | --- | ---: | ---: | ---: | ---: | --- |
| si8 | 1 | True | 15 | 340 | 2416 | 0.000000e+00 | True |
| k444s | 1 | True | 10 | 1686 | 11181 | 0.000000e+00 | True |
| syn20-skew | 1 | True | 6 | 13 | 255 | 0.000000e+00 | True |
| syn4-elongated | 1 | True | 9 | 28 | 171 | 0.000000e+00 | True |
| si8 | 3 | True | 16 | 242 | 2072 | 3.552714e-14 | False |
| k444s | 3 | True | 12 | 1627 | 10582 | 2.131628e-14 | False |
| syn20-skew | 3 | diverged (exit 4) | 6 completed | 13 | 174 | — | no |
| syn4-elongated | 3 | True | 17 | 49 | 270 | -3.552714e-15 | False |
| si8 | 10 | True | 16 | 264 | 1794 | -4.760636e-13 | False |
| k444s | 10 | True | 15 | 1745 | 10896 | 1.598721e-14 | False |
| syn20-skew | 10 | diverged (exit 4) | 6 completed | 13 | 170 | — | no |
| syn4-elongated | 10 | True | 19 | 46 | 261 | 0.000000e+00 | False |
| si8 | 30 | True | 21 | 278 | 1722 | -4.938272e-13 | False |
| k444s | 30 | True | 35 | 2655 | 16781 | 1.509903e-14 | False |
| syn20-skew | 30 | True | 8 | 14 | 296 | -1.136868e-13 | False |
| syn4-elongated | 30 | True | 27 | 49 | 309 | -1.776357e-15 | False |
| si8 | 100 | True | 57 | 288 | 2280 | -5.009326e-13 | False |
| k444s | 100 | True | 58 | 3042 | 22046 | 3.552714e-15 | False |
| syn20-skew | 100 | True | 9 | 16 | 319 | -5.684342e-14 | False |
| syn4-elongated | 100 | True | 32 | 73 | 371 | 0.000000e+00 | False |

No tested factor above 1 preserves the SCF count on any of the four cases. This is a measured bound within this capped, fixed-factor family, not a bound on every possible adaptive schedule. Equal iteration count alone would not establish an identical trajectory.

## H7: band-count experiments

| Case | Extra bands | Converged | SCF count | Total H vectors | Energy difference (Ha) | Engine seconds |
| --- | ---: | --- | ---: | ---: | ---: | ---: |
| si8 | 4 | True | 15 | 2416 | 0.000000e+00 | 9.130853 |
| k444s | 4 | True | 10 | 11181 | 0.000000e+00 | 8.594726 |
| syn20-skew | 4 | True | 6 | 255 | 0.000000e+00 | 0.031066 |
| syn4-elongated | 4 | True | 9 | 171 | 0.000000e+00 | 0.005322 |
| si8 | 1 | True | 15 | 2444 | -1.421085e-13 | 7.607048 |
| k444s | 1 | True | 10 | 8310 | 6.572520e-14 | 4.598597 |
| syn20-skew | 1 | False | 100 | 3260 | -2.793308e+01 | 0.100437 |
| syn4-elongated | 1 | True | 9 | 164 | 3.552714e-15 | 0.007258 |
| si8 | 2 | True | 17 | 3249 | 1.421085e-14 | 7.181081 |
| k444s | 2 | True | 10 | 9068 | 2.664535e-14 | 2.793621 |
| syn20-skew | 2 | True | 9 | 295 | 2.938398e+01 | 0.040231 |
| syn4-elongated | 2 | True | 9 | 149 | -3.552714e-15 | 0.011309 |

## Per-SCF baseline: si8

Early/late count band roots over all k points. Early means permanently below its threshold before the last Davidson iteration; late means first permanently below at the last iteration. Buffer fraction attributes corrections to their originating root before orthogonalization, not a physical decomposition of mixed vectors. Unknown lineage is reported separately.

| SCF | Tolerance (Ha residual) | Steps | H blocks | H vectors | H / band / k | Restarts | Early / late | Buffer-origin % | Unknown-origin H | Certificate H | Potential / diagonalization / density / energy / mixing (%) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |
| 1 | 1.0000e-02 | 9 | 9 | 104 | 5.2000 | 0 | 19 / 1 | 14.42 | 0 | 0 | 1.0 / 93.2 / 2.8 / 2.9 / 0.1 |
| 2 | 2.6159e-03 | 10 | 10 | 113 | 5.6500 | 0 | 18 / 2 | 18.58 | 0 | 0 | 0.8 / 94.2 / 2.2 / 2.6 / 0.1 |
| 3 | 1.3661e-03 | 10 | 10 | 123 | 6.1500 | 1 | 17 / 3 | 23.58 | 0 | 0 | 0.8 / 94.2 / 2.2 / 2.6 / 0.2 |
| 4 | 3.9618e-04 | 24 | 24 | 150 | 7.5000 | 1 | 19 / 1 | 20.00 | 0 | 0 | 0.6 / 95.2 / 1.8 / 2.1 / 0.3 |
| 5 | 2.8052e-04 | 11 | 11 | 125 | 6.2500 | 1 | 19 / 1 | 24.80 | 0 | 0 | 0.8 / 93.0 / 2.4 / 2.9 / 0.9 |
| 6 | 9.1345e-05 | 13 | 13 | 129 | 6.4500 | 1 | 19 / 1 | 24.81 | 0 | 0 | 0.7 / 94.1 / 1.9 / 2.6 / 0.7 |
| 7 | 2.9999e-05 | 15 | 15 | 146 | 7.3000 | 1 | 17 / 3 | 30.82 | 0 | 0 | 0.9 / 94.3 / 2.0 / 2.0 / 0.7 |
| 8 | 1.8234e-05 | 16 | 16 | 146 | 7.3000 | 1 | 19 / 1 | 31.51 | 0 | 0 | 0.6 / 94.3 / 1.9 / 2.4 / 0.8 |
| 9 | 6.0750e-06 | 33 | 33 | 212 | 10.6000 | 1 | 18 / 2 | 44.81 | 0 | 0 | 0.4 / 97.2 / 0.9 / 1.0 / 0.5 |
| 10 | 8.2866e-07 | 34 | 34 | 193 | 9.6500 | 1 | 19 / 1 | 38.34 | 0 | 0 | 0.5 / 96.2 / 1.3 / 1.5 / 0.5 |
| 11 | 4.8212e-07 | 18 | 18 | 135 | 6.7500 | 1 | 19 / 1 | 24.44 | 0 | 0 | 0.8 / 93.2 / 1.9 / 2.8 / 1.3 |
| 12 | 8.9188e-08 | 52 | 52 | 192 | 9.6000 | 1 | 19 / 1 | 24.48 | 0 | 0 | 0.6 / 97.0 / 1.0 / 1.1 / 0.4 |
| 13 | 6.8656e-08 | 45 | 45 | 165 | 8.2500 | 1 | 19 / 1 | 21.21 | 0 | 0 | 0.5 / 97.0 / 1.0 / 1.1 / 0.5 |
| 14 | 2.0070e-08 | 23 | 23 | 131 | 6.5500 | 1 | 19 / 1 | 18.32 | 0 | 0 | 0.6 / 95.5 / 1.6 / 1.7 / 0.6 |
| 15 | 6.0010e-09 | 27 | 28 | 148 | 7.4000 | 1 | 19 / 1 | 16.89 | 0 | 204 | 0.6 / 97.2 / 1.1 / 1.2 / 0.0 |

## Per-SCF baseline: k444s

Early/late count band roots over all k points. Early means permanently below its threshold before the last Davidson iteration; late means first permanently below at the last iteration. Buffer fraction attributes corrections to their originating root before orthogonalization, not a physical decomposition of mixed vectors. Unknown lineage is reported separately.

| SCF | Tolerance (Ha residual) | Steps | H blocks | H vectors | H / band / k | Restarts | Early / late | Buffer-origin % | Unknown-origin H | Certificate H | Potential / diagonalization / density / energy / mixing (%) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |
| 1 | 1.0000e-02 | 225 | 225 | 1410 | 5.5078 | 0 | 206 / 50 | 42.41 | 0 | 0 | 0.2 / 89.5 / 3.7 / 6.6 / 0.0 |
| 2 | 5.6798e-03 | 178 | 178 | 1129 | 4.4102 | 0 | 189 / 67 | 42.34 | 0 | 0 | 0.2 / 86.8 / 4.3 / 8.6 / 0.0 |
| 3 | 1.9872e-03 | 173 | 173 | 1159 | 4.5273 | 0 | 170 / 86 | 42.36 | 0 | 0 | 0.2 / 87.2 / 4.0 / 8.4 / 0.0 |
| 4 | 3.4344e-04 | 172 | 172 | 1003 | 3.9180 | 0 | 197 / 59 | 45.36 | 0 | 0 | 0.3 / 85.3 / 4.6 / 9.7 / 0.1 |
| 5 | 1.2501e-04 | 170 | 170 | 1144 | 4.4688 | 0 | 167 / 89 | 42.40 | 0 | 0 | 0.2 / 86.9 / 4.1 / 8.6 / 0.1 |
| 6 | 8.2253e-06 | 166 | 166 | 1017 | 3.9727 | 0 | 197 / 59 | 43.07 | 0 | 0 | 0.3 / 85.4 / 4.6 / 9.6 / 0.2 |
| 7 | 2.4059e-06 | 147 | 147 | 943 | 3.6836 | 0 | 159 / 97 | 33.83 | 0 | 0 | 0.3 / 84.2 / 5.0 / 10.3 / 0.2 |
| 8 | 1.6293e-07 | 148 | 148 | 749 | 2.9258 | 0 | 179 / 77 | 13.22 | 0 | 0 | 0.4 / 80.2 / 6.2 / 12.9 / 0.3 |
| 9 | 5.7187e-08 | 160 | 160 | 752 | 2.9375 | 0 | 224 / 32 | 12.77 | 0 | 0 | 0.4 / 80.5 / 6.1 / 12.8 / 0.3 |
| 10 | 4.5450e-09 | 147 | 179 | 984 | 3.8438 | 0 | 178 / 78 | 19.72 | 0 | 891 | 0.3 / 83.4 / 5.3 / 11.1 / 0.0 |

## Per-SCF baseline: syn20-skew

Early/late count band roots over all k points. Early means permanently below its threshold before the last Davidson iteration; late means first permanently below at the last iteration. Buffer fraction attributes corrections to their originating root before orthogonalization, not a physical decomposition of mixed vectors. Unknown lineage is reported separately.

| SCF | Tolerance (Ha residual) | Steps | H blocks | H vectors | H / band / k | Restarts | Early / late | Buffer-origin % | Unknown-origin H | Certificate H | Potential / diagonalization / density / energy / mixing (%) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |
| 1 | 1.0000e-02 | 3 | 3 | 35 | 1.4583 | 0 | 22 / 2 | 8.57 | 0 | 0 | 2.7 / 72.7 / 9.1 / 15.2 / 0.3 |
| 2 | 8.2126e-04 | 2 | 2 | 26 | 1.0833 | 0 | 22 / 2 | 11.54 | 0 | 0 | 3.4 / 60.0 / 14.1 / 21.3 / 1.2 |
| 3 | 4.1023e-04 | 2 | 2 | 29 | 1.2083 | 0 | 19 / 5 | 13.79 | 0 | 0 | 3.2 / 60.8 / 13.2 / 21.1 / 1.7 |
| 4 | 4.7553e-07 | 2 | 2 | 30 | 1.2500 | 0 | 18 / 6 | 10.00 | 0 | 0 | 3.0 / 62.5 / 12.8 / 19.2 / 2.5 |
| 5 | 5.1168e-08 | 2 | 2 | 28 | 1.1667 | 0 | 20 / 4 | 10.71 | 0 | 0 | 3.3 / 60.6 / 13.2 / 19.9 / 3.1 |
| 6 | 2.8679e-09 | 2 | 3 | 57 | 2.3750 | 0 | 15 / 9 | 10.53 | 0 | 50 | 2.2 / 74.9 / 9.1 / 13.9 / 0.0 |

## Per-SCF baseline: syn4-elongated

Early/late count band roots over all k points. Early means permanently below its threshold before the last Davidson iteration; late means first permanently below at the last iteration. Buffer fraction attributes corrections to their originating root before orthogonalization, not a physical decomposition of mixed vectors. Unknown lineage is reported separately.

| SCF | Tolerance (Ha residual) | Steps | H blocks | H vectors | H / band / k | Restarts | Early / late | Buffer-origin % | Unknown-origin H | Certificate H | Potential / diagonalization / density / energy / mixing (%) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- | ---: | ---: | ---: | --- |
| 1 | 1.0000e-02 | 2 | 2 | 12 | 1.5000 | 0 | 2 / 6 | 25.00 | 4 | 0 | 8.5 / 70.1 / 6.0 / 14.3 / 1.0 |
| 2 | 3.8550e-03 | 2 | 2 | 12 | 1.5000 | 0 | 4 / 4 | 50.00 | 0 | 0 | 10.7 / 57.5 / 8.1 / 19.4 / 4.2 |
| 3 | 1.8603e-03 | 3 | 3 | 16 | 2.0000 | 0 | 6 / 2 | 50.00 | 0 | 0 | 6.6 / 73.1 / 5.2 / 11.7 / 3.4 |
| 4 | 2.5139e-05 | 3 | 3 | 21 | 2.6250 | 0 | 3 / 5 | 33.33 | 0 | 0 | 6.2 / 73.9 / 4.8 / 10.8 / 4.4 |
| 5 | 5.5172e-06 | 3 | 3 | 15 | 1.8750 | 0 | 6 / 2 | 26.67 | 0 | 0 | 7.5 / 66.6 / 5.8 / 13.2 / 6.9 |
| 6 | 4.2065e-07 | 4 | 4 | 19 | 2.3750 | 0 | 6 / 2 | 15.79 | 0 | 0 | 5.7 / 74.3 / 4.4 / 9.9 / 5.7 |
| 7 | 1.1732e-07 | 4 | 4 | 18 | 2.2500 | 0 | 6 / 2 | 16.67 | 0 | 0 | 6.4 / 68.6 / 4.9 / 11.2 / 8.8 |
| 8 | 1.2438e-08 | 4 | 4 | 19 | 2.3750 | 0 | 5 / 3 | 15.79 | 0 | 0 | 6.3 / 69.0 / 5.3 / 11.4 / 8.0 |
| 9 | 1.0000e-09 | 3 | 4 | 18 | 2.2500 | 0 | 7 / 1 | 33.33 | 0 | 21 | 8.7 / 69.6 / 6.7 / 15.1 / 0.0 |

## si8 alignment with relayed QE data

QE data are relayed, not rerun here. QE uses eigenvalue-change thresholds in Ry; our residual norm thresholds in Ha are different tests. Local checkout differs from relayed fdce016. The fifteenth local iteration must not be dropped to claim equal SCF convergence.

| SCF | QE ethr (Ry) | QE steps | Relayed BAANDI steps | Local BAANDI steps | Local residual tolerance (Ha) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.000e-02 | 3 | 9 | 9 | 1.0000e-02 |
| 2 | 1.360e-03 | 4 | 10 | 10 | 2.6159e-03 |
| 3 | 6.030e-05 | 4 | 9 | 10 | 1.3661e-03 |
| 4 | 2.550e-06 | 12 | 15 | 24 | 3.9618e-04 |
| 5 | 2.550e-06 | 2 | 10 | 11 | 2.8052e-04 |
| 6 | 4.350e-08 | 7 | 14 | 13 | 9.1345e-05 |
| 7 | 4.350e-08 | 2 | 14 | 15 | 2.9999e-05 |
| 8 | 2.450e-08 | 2 | 16 | 16 | 1.8234e-05 |
| 9 | 5.820e-09 | 11 | 29 | 33 | 6.0750e-06 |
| 10 | 1.370e-09 | 1 | 19 | 34 | 8.2866e-07 |
| 11 | 3.360e-10 | 11 | 16 | 18 | 4.8212e-07 |
| 12 | 3.360e-10 | 2 | 22 | 52 | 8.9188e-08 |
| 13 | 9.940e-11 | 20 | 54 | 45 | 6.8656e-08 |
| 14 | 5.940e-11 | 2 | 21 | 23 | 2.0070e-08 |
| 15 | — | — | 32 | 27 | 6.0010e-09 |

QE totals: 83 steps and 98 h_psi block calls. QE block widths are not present in the relay, so vector-application totals cannot be reconstructed. The k444s relay is internally ambiguous about five versus six SCF iterations and omits one threshold/error pairing; it is not expanded into invented iteration rows.

## Stored k222 / k222s QE outputs

| Case | SCF | ethr (Ry) | Average iterations |
| --- | ---: | ---: | ---: |
| k222_c | 1 | 1.000e-02 | 3.9 |
| k222_c | 2 | 1.310e-03 | 2.0 |
| k222_c | 3 | 5.710e-05 | 3.0 |
| k222_c | 4 | 1.010e-06 | 2.8 |
| k222_c | 5 | 1.060e-08 | 6.5 |
| k222_c | 6 | 1.610e-10 | 4.6 |
| k222_c | 7 | 8.820e-12 | 5.0 |
| k222_c | 8 | 2.920e-13 | 4.5 |
| k222_c | 9 | 1.000e-13 | 1.5 |
| k222s_c | 1 | 1.000e-02 | 4.5 |
| k222s_c | 1 | 8.200e-04 | 3.8 |
| k222s_c | 2 | 8.810e-04 | 1.2 |
| k222s_c | 3 | 3.180e-05 | 2.5 |
| k222s_c | 4 | 8.560e-07 | 3.2 |
| k222s_c | 5 | 5.660e-09 | 6.2 |
| k222s_c | 6 | 4.870e-10 | 2.0 |
| k222s_c | 7 | 9.540e-12 | 4.2 |
| k222s_c | 8 | 3.690e-13 | 5.0 |
| k222s_c | 9 | 1.000e-13 | 2.0 |
