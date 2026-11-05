# QE 7.6 iteration and timing data on dev08 (relayed 2026-11-05)

Provenance: extracted by the remote measurement session from existing pw.x outputs (`verbosity='high'`) in `~/rustdft_p6/qe2/` on dev08: si8 `fin_1e-8_np1_r1.log`, `fin_1e-12_np1_r1.log`, k444s `fin_1e-8_np1_r1.log`; BAANDI side from `~/rustdft_p6/hist/out/results.json` (main fdce016, pure build, si8, T=1). Relayed as text, not independently verified here. QE runs: `K_POINTS automatic 1 1 1` (complex Gamma, not the real-wavefunction gamma trick), nbnd 20, same FFT grid as BAANDI, one MPI process.

Unit caution: QE `h_psi` is called once per Davidson step for the whole band block (nbnd 20); BAANDI `operator_applications` may count per band. Comparable unit: Davidson steps (QE 83 versus BAANDI 290).

## si8, N=1, conv_thr 1e-8, per SCF iteration

```
it   QE ethr      QE avg_it   QE est_scf_acc(Ry)   BAANDI dav   BAANDI ops   BAANDI density_residual
 1   1.00E-02      3.0         0.43619171             9           104        0.02618218329066498
 2   1.36E-03      4.0         0.01928503            10           116        0.013708363221295374
 3   6.03E-05      4.0         0.00081500             9           118        0.004179047899005566
 4   2.55E-06     12.0         0.00152229            15           136        0.003028842521435974
 5   2.55E-06      2.0         0.00001393            10           119        0.0009501788147661855
 6   4.35E-08      7.0         0.00002386            14           136        0.00032558660590850094
 7   4.35E-08      2.0         0.00000783            14           134        0.00019274637064643088
 8   2.45E-08      2.0         0.00000186            16           143        6.467064979225458e-05
 9   5.82E-09     11.0         0.00000044            29           183        1.0438126996325737e-05
10   1.37E-09      1.0         0.00000011            19           148        3.916624419500134e-06
11   3.36E-10     11.0         0.00000018            16           127        8.007704768342913e-07
12   3.36E-10      2.0         0.00000003            22           132        2.4702780148496874e-07
13   9.94E-11     20.0         0.00000002            54           200        2.0957092902446895e-07
14   5.94E-11      2.0         5.4E-09               21           128        3.738561712796787e-08
15      -           -              -                 32           146        1.1338673421536119e-08
```
QE: 14 SCF iterations, 83 Davidson steps in total (mean 5.93 per SCF iteration), total energy -62.64362757 Ry. BAANDI: 15 SCF iterations, 290 Davidson iterations (19.33 per SCF iteration), 2070 operator applications (7.14 per step), total energy -31.32180360103162 Ha = -62.64360720 Ry (2.04e-5 Ry from QE: different convergence thresholds and basis handling, so the work comparison is not at identical accuracy). Ratios BAANDI/QE: SCF iterations 1.07, Davidson steps 3.49, Davidson per SCF iteration 3.26.

## si8, N=1, conv_thr 1e-12 (QE)
First 14 iterations identical to the 1e-8 run; continues: it 15 ethr 1.70E-11 avg 12.0; 16: 1.69E-12, 2.0; 17: 9.62E-13, 15.0; 18: 1.00E-13, 1.0; 19: 1.00E-13, 2.0; 20: 1.00E-13, 1.0. 20 SCF iterations, 116 Davidson steps.

## QE time breakdown (wall, calls)
si8 conv 1e-8, N=1: PWSCF 0.68 s; init_run 0.05; electrons 0.61; c_bands 0.47 (14 calls); sum_band 0.07; v_of_rho 0.03; mix_rho 0.04; cegterg 0.45 (14); cdiaghg 0.03 (97); cegterg:over 0.03 (83); cegterg:upda 0.02 (83); cegterg:last 0.02 (39); h_psi 0.36 (98 calls); g_psi 0.01 (83); vloc_psi 0.32 (98); h_psi:calbec 0.02; add_vuspsi 0.02; fftw 0.29 s (2138 calls).
si8 conv 1e-12: electrons 0.84 s; c_bands 0.64 (20); cegterg 0.61; cdiaghg 0.04 (136); h_psi 0.48 (137 calls); vloc_psi 0.43.
k444s conv 1e-8, N=1: PWSCF 1.19 s; electrons 1.09; c_bands 0.95 (6 SCF iterations x 32 k points); cegterg 0.92 (192 calls); cdiaghg 0.07 (946); h_psi 0.82 (1010 calls, 31.6 per k point); vloc_psi 0.79; fftw 0.76 s (12262 calls). Per iteration: ethr 1.00E-02, 7.02E-04, 7.59E-04, 2.34E-05, 7.34E-07 (last 9.88E-10), avg iterations 4.7 3.8 1.2 3.5 4.5 7.0, est_scf_acc 0.06073066, 0.00187486, 0.00005869, 0.00000008, 6.0E-09; converged in 5 to 6 SCF iterations.

Observation (relaying session): QE starts with a loose ethr (1e-2) and tightens it with the density residual; its Davidson step count per SCF iteration is irregular (3, 4, 4, 12, 2, 7, 2, 2, 11, 1, 11, 2, 20, 2): large right after ethr tightens, 1 to 2 otherwise. BAANDI runs 9 to 54 Davidson iterations in every SCF iteration.
