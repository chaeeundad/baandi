# dev08 BAANDI versus QE 7.6 on Si and MgO bulk (relayed summary, 2026-11-07)

Relayed by the dft5090/master session as text; the raw run directories are on dev08 under `/scr/mino/rustdft/{baandi,baandi_def,qe}/`. This file is a summary, not the raw evidence: the comparison tool `scripts/compare_materials_qe.py` exited 3 on all six cases because QE 7.6 prints `Exchange-correlation= PBE` twice when `input_dft` is set (fixed in commit b820a3c: identical repeats accepted, conflicting values still rejected); the numbers below were aggregated by hand from the outputs by the relaying session, and the density comparison in G space has not been run yet. Nothing here is a validation claim or an acceptance decision.

Setup: BAANDI main `32770e2`, rustc 1.99.1; QE 7.6 `/home/mino/fq/builds/orig-mpi-dfti/bin/pw.x` (MKL DFTI), Open MPI 5.0.8, 1 rank, OMP/MKL threads 1; QE inputs from `make_qe_input.py --disk-io nowf --conv-thr 1e-12`, not edited by hand; PseudoDojo NC-SR v0.4 PBE standard UPF, SHA-256 Si `39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282`, Mg `15cdd9380687a47a1f4bc4ab5fe1f0c6d9c50891350fd8534c8f292b3ce0be14`, O `2a543889f6522e2a1f94e2c2ff1477e61cb89f25c67b27c54940e8263ec6fb11`. BAANDI pure binaries: default `26bc188c0e599040e8ab0f663ec12a9c864e18c2e8f7069f93c4cd84f08db1cc`, native `8dc93345efcbdc438afa6d29909fc027a4cf96c20af0f4c1d10c8fdb5229ea8f`. UTC 2026-11-07T02:45:29Z, host dev08.

## Total energy

| Case | BAANDI (Ha) | BAANDI (Ry) | QE (Ry) | BAANDI minus QE (Ry) |
| --- | ---: | ---: | ---: | ---: |
| si-20 | -8.462317836906 | -16.92463567 | -16.92463785 | +2.176e-06 |
| si-30 | -8.462329006518 | -16.92465801 | -16.92466016 | +2.147e-06 |
| si-40 | -8.462336822427 | -16.92467364 | -16.92467579 | +2.145e-06 |
| mgo-40 | -76.050880851711 | -152.10176170 | -152.10176985 | +8.147e-06 |
| mgo-50 | -76.050969240430 | -152.10193848 | -152.10194662 | +8.139e-06 |
| mgo-60 | -76.050977170773 | -152.10195434 | -152.10196248 | +8.138e-06 |

BAANDI is higher than QE in every case, and the difference is almost independent of the cutoff. This is an unexplained systematic difference under investigation (see docs/research/qe-study/qe-pseudo-treatment-differences.md when written); it is not a tolerance question and no cap is set here.

## Iterations, electrons, bands, grid

| Case | SCF BAANDI / QE | Electrons B / Q | Bands B / Q | FFT grid |
| --- | ---: | ---: | ---: | --- |
| si-20 | 15 / 8 | 8 / 8 | 8 / 8 | 30x30x30 |
| si-30 | 15 / 8 | 8 / 8 | 8 / 8 | 36x36x36 |
| si-40 | 15 / 8 | 8 / 8 | 8 / 8 | 45x45x45 |
| mgo-40 | 18 / 10 | 16 / 16 | 12 / 12 | 36x36x36 |
| mgo-50 | 19 / 10 | 16 / 16 | 12 / 12 | 36x36x36 |
| mgo-60 | 19 / 10 | 16 / 16 | 12 / 12 | 40x40x40 |

BAANDI (strict, certified, conv 1e-8 class) uses about twice QE's SCF iterations; the definitions of convergence differ (see ADR 0035).

## Eigenvalues and gaps

Largest |BAANDI - QE| eigenvalue difference, eV (QE prints four decimals, so these are at or below its rounding resolution, not real differences): si-20 4.482e-05 (band 5), si-30 3.973e-05 (1), si-40 3.609e-05 (4), mgo-40 5.719e-05 (3), mgo-50 4.551e-05 (9), mgo-60 4.273e-05 (4). Rounded to four decimals, BAANDI's si-20 eigenvalues equal QE's printed ones digit for digit.
Indirect gap BAANDI versus QE HOMO-LUMO (eV): si-20 1.179793 / 1.1798, si-30 1.179770 / 1.1797, si-40 1.179794 / 1.1798, mgo-40 6.383667 / 6.3837, mgo-50 6.383602 / 6.3836, mgo-60 6.383583 / 6.3836 (differences at most 7.0e-05, again at QE print resolution).

## QE cost (single run each, N=1, `/usr/bin/time -v`)

| Case | wall (s) | max RSS (kB) | QE iterations |
| --- | ---: | ---: | ---: |
| si-20 | 4.11 | 61636 | 8 |
| si-30 | 5.36 | 76756 | 8 |
| si-40 | 13.58 | 93548 | 8 |
| mgo-40 | 2.79 | 67916 | 10 |
| mgo-50 | 2.97 | 76932 | 10 |
| mgo-60 | 4.37 | 84668 | 10 |

(QE ran with conv_thr 1e-12 Ry; times are single runs, not a speed benchmark.)

## Native versus default BAANDI build

Total energies are bit-identical (difference 0.000e+00) in all six cases, including the NLCC and two-species cases.

## Inputs

qe.in SHA-256: si-20 `efebc34241d7abb1613dc449fd8b5292afb3f8f1ccaf4c30a0dbeab251c1cddc`, si-30 `3ce812a1b73d8a2332eebc349dbecda23ae771a91b554d662e9c51a3c1cc20f8`, si-40 `20926338327149d42f2bb8acbe1500cf23a816ff4acaa0177877c4ef1fac2a52`, mgo-40 `e35bd14949aea3f73e45b8942504b59a4337bff770fbc10fbf2f9cf7ddf37ff6`, mgo-50 `0a30f0a1114036d4a194cfe72fab29f9bb8def592313effd1557e039b4e1956a`, mgo-60 `f447fcc7805970fa41aa1d81b9a6dc22152471e5c4f69fcdb9110042484d3c0a`.

## Addendum 2026-11-07 (relayed): energy-term breakdown, density comparison, full-precision eigenvalues

Energy terms in Ry, BAANDI minus QE (QE terms are the literal `qe.out` lines; BAANDI one-electron = kinetic + local + nonlocal):

| Term | si-20 QE | si-20 BAANDI minus QE | mgo-40 QE | mgo-40 BAANDI minus QE |
| --- | ---: | ---: | ---: | ---: |
| one-electron | 4.96587879 | +1.90e-07 | -74.17982927 | +2.10e-07 |
| Hartree | 1.09536042 | -2.60e-07 | 41.45948608 | +1.00e-06 |
| xc | -6.18494746 | +2.220e-06 | -24.49976776 | +6.940e-06 |
| Ewald | -16.80092960 | +3.0e-08 | -94.88165890 | 0 |
| total | -16.92463785 | +2.180e-06 | -152.10176985 | +8.150e-06 |

The xc term alone carries 102 percent (Si) and 85 percent (MgO) of the total energy difference; kinetic + local + nonlocal agree to 2.1e-07 Ry and the Ewald sum to 3e-08 Ry (Si) and exactly (MgO). Lattice sums and operator application are therefore not the cause. Candidates in the xc path: NLCC core density handling, density floor (`density_floor_electrons_bohr3 = 1e-10`), GGA thresholds. Not yet separated by numbers.

Density comparison in G space (fixed script, all six cases exit 0, status REPORTED, acceptance not assessed): the G sets are identical in every case (only_ours 0, only_qe 0); largest coefficient difference 1.5e-09 to 6.4e-09 (all at Miller [1,1,1]); rho(G=0) difference at most 8.3e-17; relative L2 difference 1.0e-07 to 1.6e-07; integrated electrons within 1.1e-14 of the integer.

Eigenvalues from QE's `data-file-schema.xml` (full precision) versus BAANDI, largest |difference| in Ha: si-20 2.0026e-07, si-30 1.8929e-07, si-40 2.0445e-07, mgo-40 5.5721e-07, mgo-50 5.4659e-07, mgo-60 5.6127e-07 (k points 32 for Si, 14 for MgO, bands 8 and 12). The stdout-based differences of about 2e-06 Ha were dominated by QE's four-decimal printing. `scripts/compare_materials_qe.py` still reads eigenvalues from stdout; reading the XML would raise its precision by about a factor of ten (open item).

## Addendum 2 (relayed, 2026-11-07): effect of `[xc] pw92_prefactor = "original"` on the six cases

BAANDI main `946515e`, pure native build sha256 `7daad734656942b0e5d3f64abdedf0f60238aa1b92a32aec38c73ba0a5fb7896`; QE values as in the table above. The runner `run_materials_smoke.py` has no extra-key argument, so the relaying session copied the generated working inputs, inserted the single line `pw92_prefactor = "original"` below `[xc]`, and called `dft-cli run` directly; the "off" side was run the same way and its energies are bit-identical to the runner results reported earlier. 12 runs, all exit 0, SCF counts unchanged (Si 15/15/15, MgO 18/19/19). Raw output: `/scr/mino/rustdft/pw92/<case>_{off,original}/`.

| Case | BAANDI minus QE, option off (Ry) | option original (Ry) | reduction factor | remaining per atom (Ha) |
| --- | ---: | ---: | ---: | ---: |
| si-20 | 2.1762e-06 | 7.5937e-08 | 28.7 | 3.7968e-08 |
| si-30 | 2.1470e-06 | 4.6713e-08 | 46.0 | 2.3357e-08 |
| si-40 | 2.1451e-06 | 4.4898e-08 | 47.8 | 2.2449e-08 |
| mgo-40 | 8.1466e-06 | 1.2760e-06 | 6.4 | 6.3800e-07 |
| mgo-50 | 8.1391e-06 | 1.2685e-06 | 6.4 | 6.3426e-07 |
| mgo-60 | 8.1385e-06 | 1.2679e-06 | 6.4 | 6.3394e-07 |

These agree with the P27A fixed-density predictions (Si +7.97e-08, MgO +1.279e-06 Ry) to about 0.2 percent for MgO. Si's remainder falls with the cutoff (basis-convergence tail); MgO's remainder of 1.27e-06 Ry does not depend on the cutoff, so a second cutoff-independent term exists. P27C showed that BAANDI's local-potential alpha-Z term is very sensitive to the radial endpoint (about 8e-5 Ry for 0.5 Bohr), so a one-mesh-point endpoint convention difference would be enough; this is a hypothesis under study, not a finding.

# dev08 standard profile versus strict, and versus QE at equal accuracy (relayed, 2026-11-07)

Relayed text from the dft5090/master session; raw output on dev08 in `/scr/mino/rustdft/p4/` (BAANDI) and `/scr/mino/rustdft/p4qe/` (QE with `data-file-schema.xml`). BAANDI main `849e0ae` (without the later setup/audit speedups P28A/P28B), rustc 1.99.1; binaries pure default `0593835297aa9fc5ac8ed123620432a345bbc103b13fc38556823da694cd717a`, pure native `ed6700a0f52858de09ea461b7437e497c8e37a5612a34c5ee18a315de4b6c133`, MKL default `f66b219d30cb891610e08540cb7a8fc0874f68b5dad519ef2bf05e0e54a9eb39`, MKL native `84a8681a1d0dcf291aa0f38a2881a58eaef8e8c20a0e489b3a2ba4d28d1efc48`. Load below 0.30 before each combination, 3-run medians (two outliers re-measured with 9 runs). gamma_real was not used in this measurement (batch-7 B pending). Standard rows are labelled `standard`, never certified.

Wall time in seconds (strict / standard / speed-up):

| Case | Backend | Build | T | strict | standard | speed-up |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| si8 | pure | default | 1 | 5.0108 | 1.6567 | 3.03 |
| si8 | pure | default | 8 | 3.6242 | 1.2024 | 3.01 |
| si8 | pure | native | 1 | 3.8209 | 1.3342 (9 runs: 1.3195) | 2.86 |
| si8 | pure | native | 8 | 2.4963 | 0.9209 | 2.71 |
| si8 | MKL | default | 1 | 4.3793 | 2.3648 | 1.85 |
| si8 | MKL | default | 8 | 3.2754 | 1.7068 | 1.92 |
| si8 | MKL | native | 1 | 4.1156 | 2.2799 | 1.81 |
| si8 | MKL | native | 8 | 3.0818 | 1.6635 | 1.85 |
| k444s | pure | default | 1 | 5.2623 | 2.9678 | 1.77 |
| k444s | pure | default | 8 | 0.9829 | 0.6295 (9 runs: 0.5963) | 1.56 |
| k444s | pure | native | 1 | 4.1358 | 2.4020 | 1.72 |
| k444s | pure | native | 8 | 0.8802 | 0.5418 | 1.63 |
| k444s | MKL | default | 1 | 3.9841 | 2.3574 | 1.69 |
| k444s | MKL | default | 8 | 0.8375 | 0.5362 | 1.56 |
| k444s | MKL | native | 1 | 3.7765 | 2.2464 | 1.68 |
| k444s | MKL | native | 8 | 0.7911 | 0.5098 | 1.55 |

Spread below 5 percent except the two re-measured entries. Best combination: si8 pure native standard T=8 at 0.9209 s. Default and native builds give bit-identical energies and density hashes in all eight case/backend/profile combinations; density hashes are identical across thread counts. Note: this run had no k-point symmetry and the Mac observations of the standard profile (3.2x on si8) are confirmed on dev08 for the pure build (3.0x).

Work counts at T=1 (SCF iterations / Davidson steps / H vector applications):

| Case | Backend | strict | standard | ops reduction |
| --- | --- | --- | --- | ---: |
| si8 | pure | 15 / 279 / 2044 | 15 / 77 / 772 | 2.65 |
| si8 | MKL | 15 / 302 / 2139 | **27** / 110 / 1297 | 1.65 |
| k444s | pure (= MKL) | 10 / 1686 / 10290 | 9 / 866 / 5864 | 1.76 |

**Anomaly:** with the MKL build the standard profile on si8 needs 27 SCF iterations instead of 15 (identical for T = 1 and 8, default and native builds; energies equal in all of them), halving the profile's gain (1.85x versus 3.0x for pure). k444s shows no such effect (10 to 9 on both backends). Open: per-iteration histories (density residual, standard error estimate, eigen thresholds) for si8 standard, MKL versus pure, requested from dev08.

Total energies (Ha) strict / standard / difference: si8 pure -31.321803601031636 / -31.321803599478528 / 1.553e-09; si8 MKL -31.32180360103142 / -31.321803598673775 / 2.358e-09; k444s pure -7.92597612664809 / -7.92597612651481 / 1.333e-10; k444s MKL -7.9259761266480915 / -7.925976126514809 / 1.333e-10.

## QE conv_thr ladder on si8 (N=1, nowf, 3 runs; total energy from `data-file-schema.xml` in Ha, error against the conv_thr 1e-14 result of the same path)

| Path | conv_thr (Ry) | total energy (Ha) | abs error (Ha) | internal time (s) | iterations |
| --- | ---: | ---: | ---: | ---: | ---: |
| complex k=0 | 1e-6 | -31.32181368954285 | 9.703e-08 | 0.4600 | 9 |
| complex k=0 | 1e-8 | -31.32181378461762 | 1.955e-09 | 0.6800 | 14 |
| complex k=0 | 1e-10 | -31.32181378656577 | 7.141e-12 | 0.8400 | 17 |
| complex k=0 | 1e-12 | -31.32181378657216 | 7.496e-13 | 0.9133 | 20 |
| complex k=0 | 1e-14 | -31.32181378657291 | 0 | 1.6167 | 43 |
| real gamma | 1e-6 | -31.32181341221784 | 3.744e-07 | 0.2500 | 8 |
| real gamma | 1e-8 | -31.32181378453889 | 2.040e-09 | 0.3667 | 12 |
| real gamma | 1e-10 | -31.32181378656565 | 1.312e-11 | 0.4500 | 15 |
| real gamma | 1e-12 | -31.32181378657835 | 4.192e-13 | 0.5033 | 18 |
| real gamma | 1e-14 | -31.32181378657877 | 0 | 0.9100 | 42 |

The two QE paths' conv_thr 1e-14 energies differ by 5.858e-12 Ha. (QE's 8-digit stdout prints the same energy for every conv_thr from 1e-8 up, hence the XML.)

## Equal-accuracy comparison (si8, N=1 unless noted)

| Code / setting | energy error (Ha) | time (s) |
| --- | ---: | ---: |
| QE real gamma, conv_thr 1e-8 | 2.040e-09 | 0.3667 |
| QE complex k=0, conv_thr 1e-8 | 1.955e-09 | 0.6800 |
| BAANDI standard, pure native, T=1 | 1.553e-09 | 1.3342 (9 runs 1.3195) |
| BAANDI standard, MKL native, T=1 | 2.358e-09 | 2.2799 |
| BAANDI standard, pure native, T=8 | 1.553e-09 | 0.9209 |

At the same delivered-accuracy class (about 2e-9 Ha) QE real gamma is 3.6 times faster than BAANDI standard pure native on one core, and 2.5 times faster than BAANDI on eight cores; BAANDI is 1.9 times behind QE complex k=0 on one core. (BAANDI's gamma_real path and the setup/audit speedups were not in this measurement.)

# dev08 si8 real Gamma path on/off, T=1 (relayed, 2026-11-07, batch-7 B, first rows)

Relayed by the dft5090/master session. BAANDI main `b6a5416`, pure native, `taskset -c 0`, 30 runs per row, Si.pz-vbc.UPF SHA-256 `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`. si8, one thread.

| gamma_real | profile | n | wall median (s) | min | max | SCF | Davidson steps | H vectors | energy (Ha) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| off | strict | 30 | 3.764 | 3.744 | 3.823 | 15 | 279 | 2044 | -31.321803601032 |
| off | standard | 30 | 1.308 | 1.300 | 1.358 | 15 | 77 | 772 | -31.321803599479 |
| on | strict | 30 | 3.525 | 3.511 | 3.586 | 15 | 319 | 2125 | -31.321803601032 |
| on | standard | 30 | 1.300 | 1.292 | 1.355 | 17 | 85 | 864 | -31.321803599706 |

Energy difference on versus off: strict 0.0 Ha; standard 2.27e-10 Ha (the strict 1e-10 gate applies to strict only; standard uses the delivered-accuracy gate defined for P26G).

Observation: on dev08 the real Gamma path gives only -6.4 percent wall time with strict (and -0.6 percent with standard), whereas the Mac observation was 1.53x (strict) and a clear gain with standard. Davidson steps rise by 14 percent (strict) and 10 percent (standard), H vectors by 4 and 12 percent. The Mac gain does NOT transfer to dev08 for this option; the standard-profile gain (3.0x) does. Open: where the real path loses on dev08 (real GEMM kernels, tile size, FFT pair on x86 AVX-512, step count); the D0 profile of docs/performance/diagonalization-plan.md must cover real on versus off.

Note: `gamma_real_used` is only written to results.json when the option is not off (intentional, keeps off-path output unchanged); a missing key means off.

# dev08 batch-7 B complete: si8, si16, si32, real Gamma path on/off, T=1 and T=8 (relayed, 2026-11-07)

Relayed by the dft5090/master session. si8 and si16 from main `b6a5416` (pure native binary `5b827976e23902299a9cae9c86a8386a4870cb782e2812bb25a6bc2de67f3584`, 30 runs per row), si32 from main `8894b0b` (pure native `16cbfa76d5bb22252648fc3e1cab9f499b2297fc7aad13292b2d53aef50bed54`, 5 runs per row because one run is 24 to 120 s); rustc 1.99.1; `taskset -c 0` (T=1) or `0-7` (T=8), no other jobs. Si.pz-vbc.UPF SHA-256 `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`. si16: 2x1x1 conventional cell, 16 atoms, 30 Ry, Gamma; si32: 2x2x1, 32 atoms, 68 bands, 11877 plane waves, grid 72x72x36. lda-pz, Pulay, beta 0.5, etol 1e-8, dens 1e-7. SCF / Davidson steps / H vectors are identical across repeats inside a row and between T=1 and T=8; energies are identical across repeats and thread counts.

Wall time (s, median) with SCF / Davidson steps / H vectors:

| Case | T | off strict | on strict | off standard | on standard |
| --- | ---: | --- | --- | --- | --- |
| si8 | 1 | 3.764 (15/279/2044) | 3.525 (15/319/2125) | 1.308 (15/77/772) | 1.300 (17/85/864) |
| si8 | 8 | 2.466 | 3.417 | 0.924 | 1.304 |
| si16 | 1 | 15.124 (12/182/2700) | 13.319 (13/216/3053) | 5.866 (14/76/1350) | 5.205 (15/81/1424) |
| si16 | 8 | 7.456 | 10.823 | 3.143 | 4.470 |
| si32 | 1 | 119.223 (13/352/5946) | 97.475 (13/359/5874) | 38.357 (15/67/2806) | 24.002 (13/60/2380) |
| si32 | 8 | 50.735 | 68.437 | 16.845 | 19.869 |

Real path on versus off, wall ratio (negative = on faster): T=1: si8 strict -6.4 %, standard -0.6 %; si16 strict -11.9 %, standard -11.3 %; si32 strict -18.2 %, standard -37.4 % (SCF 15 versus 13: different trajectories, read as total cost to the same criterion). T=8: si8 +38.6 % / +41.1 %, si16 +45.2 % / +42.2 %, si32 +34.9 % / +18.0 % (strict / standard; the si32 standard value is partly masked by two fewer SCF iterations).
Thread scaling T=1 to T=8 (off strict / on strict / off standard / on standard): si8 1.53 / 1.03 / 1.42 / 1.00; si16 2.03 / 1.23 / 1.87 / 1.16; si32 2.35 / 1.42 / 2.28 / 1.21.
Energy difference on versus off: si8 strict 1.8e-14, si16 strict 1.2e-11, si32 strict 4.5e-11 Ha (all below 1e-10); standard 2.3e-10, 2.9e-10, 9.3e-10 Ha.

Reading: the real path gain at one thread grows with size (strict -6 %, -12 %, -18 %), but with eight threads the real path is 18 to 45 percent slower than the complex path and scales poorly in every combination. Mac check (this session, si8 strict T=8, 3 runs, median): off 1.106 s, on 1.310 s: the inversion reproduces on the Mac.

First QE batch-7 C result: si8 real gamma conv_thr 1e-8, np=1 wall 0.4659 to 0.4739 s (5 runs); np=8 0.2516 to 0.2933 s, below the 0.3 s timer-resolution rule and therefore marked invalid. QE pw.x sha256 `07c7cd7935f42b15157500ab719a5c3948825d8ce8a4487d992d7501c659658f`.

# dev08 batch-7 C: QE 7.6 real-wavefunction Gamma path, timers and conv_thr ladder (relayed, 2026-11-07)

Relayed by the dft5090/master session. QE 7.6 (`pw.x` sha256 `07c7cd7935f42b15157500ab719a5c3948825d8ce8a4487d992d7501c659658f`), `K_POINTS gamma` (the regterg/rdiaghg real path is used), `verbosity='high'`, `disk_io='nowf'`, `nosym`/`noinv`, FFT grid equal to BAANDI's (`--fft-grid ours`), OMP/MKL threads 1, np=1, 5 runs, medians; times below 0.3 s are invalid under the 0.01 s timer-resolution rule (all np=8 si8 values except conv_thr 1e-12 are therefore not used).

## Wall time and accuracy ladder (np=1)

| Case | conv_thr (Ry) | wall median (s) | min | max | SCF | estimated accuracy (Ry) | etot from XML (Ha) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| si8 | 1e-8 | 0.4669 | 0.4656 | 0.4739 | 12 | 6.6e-09 | -31.32181378565385 |
| si8 | 1e-9 | 0.5131 | 0.5120 | 0.5165 | 14 | 2.7e-10 | -31.32181378641102 |
| si8 | 1e-10 | 0.5509 | 0.5496 | 0.5520 | 15 | 5.1e-11 | -31.32181378655076 |
| si8 | 1e-12 | 0.6055 | 0.6024 | 0.6080 | 18 | 2.7e-13 | -31.32181378657689 |
| si16 | 1e-8 | 1.4132 | 1.4101 | 1.4975 | 10 | 6.3e-09 | -63.12926478535374 |
| si16 | 1e-9 | 1.6851 | 1.6838 | 1.6914 | 12 | 1.1e-10 | -63.12926478535374 |
| si16 | 1e-10 | 1.8261 | 1.8224 | 1.8308 | 13 | 2.6e-11 | -63.12926478561818 |
| si16 | 1e-12 | 1.9741 | 1.9726 | 1.9759 | 15 | 1.8e-13 | -63.12926478562362 |
| si32 | 1e-8 | 5.03 (single run) | | | | | |

Valid np=8 values: si8 conv_thr 1e-12 0.3383 s; si16 0.5559 / 0.5555 / 0.6481 / 0.6977 s for 1e-8 / 1e-9 / 1e-10 / 1e-12.

## QE internal timers (wall seconds, call counts)

si8 (nbnd 20, grid 36^3, SCF 12): electrons 0.31; c_bands 0.24 (12); regterg 0.23 (12); h_psi 0.18 (94); vloc_psi 0.17 (94); h_psi:calbec 0.01; add_vuspsi 0.00; rdiaghg 0.02 (93); regterg:over 0.01 (81); regterg:upda 0.00 (81); regterg:last 0.01 (34); g_psi 0.00 (81); sum_band 0.03; v_of_rho 0.03; mix_rho 0.02; fftw 0.13 (964); init_run 0.04; PWSCF total 0.37. (0.01 s quantization: low precision.)
si16 (nbnd 36, grid 72x36x36, SCF 11): electrons 1.20; c_bands 1.02 (11); regterg 0.99 (11); h_psi 0.76 (148); vloc_psi 0.66 (148); h_psi:calbec 0.06; add_vuspsi 0.02; rdiaghg 0.11 (146); regterg:over 0.06 (136); regterg:upda 0.04 (136); regterg:last 0.03 (32); g_psi 0.01 (136); sum_band 0.11; v_of_rho 0.04; mix_rho 0.03; fftw 0.55 (1612 calls, 341 us each); init_run 0.09; PWSCF total 1.30.
Shares of `electrons` (si16): local FFT (fftw) 45.8 %, nonlocal (calbec + add_vuspsi) 6.7 %, subspace algebra (rdiaghg + over + upda + last) 20.0 %, preconditioner 0.8 %, other (sum_band, v_of_rho, mix_rho, remainder) 17.5 %; h_psi alone is 63.3 % (the fftw timer also counts calls inside sum_band, so the first and last shares overlap).

## BAANDI tight reference (strict tightened, gamma_real off, T=1)

si8: 21 SCF, 427 Davidson steps, 2674 H vectors, 4.72 s, E = -31.321803601031771 Ha; si16: 18 / 276 / 3919, 19.74 s, E = -63.129244421617670 Ha; si32: 22 / 567 / 8995, 157.68 s, E = -126.703998096108023 Ha; all converged.

## Wall-to-wall gap at one thread (QE conv_thr 1e-8 versus BAANDI standard, both about 1e-9 Ha delivered; orchestrator arithmetic from the relayed tables, wall including process start for both)

| Case | QE real gamma (s) | BAANDI standard, best of real on/off (s) | BAANDI / QE |
| --- | ---: | ---: | ---: |
| si8 | 0.467 | 1.300 (on) | 2.8 |
| si16 | 1.413 | 5.205 (on) | 3.7 |
| si32 | 5.03 | 24.0 (on, 13 SCF) | 4.8 |

The gap grows with system size (2.8, 3.7, 4.8). Earlier statements that compared BAANDI wall time with QE's internal 0.366 s timer overstated the si8 gap (3.6x); the wall-to-wall figure is 2.8x. Equal-accuracy caveat: BAANDI standard results are uncertified and SCF counts differ between the codes. The BAANDI standard rows are single trajectories (see docs/performance/standard-robust-gain.md: the median Si8 H-work gain of the standard profile is 1.78x, so the wall gap against a median trajectory would be larger).

# dev08 batch-7 C complete: equal-accuracy comparison BAANDI versus QE real gamma, one core (relayed, 2026-11-07)

Relayed by the dft5090/master session; BAANDI T=1 pure native `8894b0b`, `taskset -c 0`, 30 runs (si32 5 runs); QE 7.6 `K_POINTS gamma`, np=1, OMP/MKL 1, 5-run medians. Each code is compared to its OWN tight reference because the absolute energies differ by a convention offset (BAANDI tight minus QE tight, default PW92 prefactor: si8 +1.019e-05, si16 +2.036e-05, si32 +4.073e-05 Ha, about 1.27e-06 Ha per atom, see the PW92 and radial-convention studies). Tight references: BAANDI strict tightened (energy 1e-12 Ha, density residual 1e-10, gamma_real off, T=1): si8 -31.321803601031771 (21 SCF, 427 Davidson steps, 2674 H vectors, 4.72 s), si16 -63.129244421617670 (18/276/3919, 19.74 s), si32 -126.703998096108023 (22/567/8995, 157.68 s); QE conv_thr 1e-12 Ry: si8 -31.32181378657689 (18 SCF), si16 -63.12926478562362 (15), si32 -126.7040388229817 (17); energies from `data-file-schema.xml`. Standard rows are `standard`, uncertified.

## BAANDI delivered accuracy (T=1)

| Case | Row | wall (s) | abs dE vs own tight (Ha) | density rel L2 | max abs drho | SCF | Davidson | H vectors |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| si8 | off strict | 3.764 | 1.350e-13 | 8.742e-08 | 2.028e-08 | 15 | 279 | 2044 |
| si8 | on strict | 3.525 | 1.172e-13 | 1.066e-07 | 2.180e-08 | 15 | 319 | 2125 |
| si8 | off standard | 1.308 | 1.553e-09 | 1.058e-05 | 1.707e-06 | 15 | 77 | 772 |
| si8 | on standard | 1.300 | 1.326e-09 | 1.204e-05 | 2.452e-06 | 17 | 85 | 864 |
| si16 | off strict | 15.124 | 1.707e-11 | 5.296e-07 | 9.461e-08 | 12 | 182 | 2700 |
| si16 | on strict | 13.319 | 5.379e-12 | 2.811e-07 | 6.101e-08 | 13 | 216 | 3053 |
| si16 | off standard | 5.866 | 1.268e-09 | 5.254e-06 | 9.420e-07 | 14 | 76 | 1350 |
| si16 | on standard | 5.205 | 9.787e-10 | 5.377e-06 | 9.568e-07 | 15 | 81 | 1424 |
| si32 | off strict | 119.223 | 2.942e-11 | 3.721e-07 | 6.149e-08 | 13 | 352 | 5946 |
| si32 | on strict | 97.475 | 7.392e-11 | 5.534e-07 | 8.887e-08 | 13 | 359 | 5874 |
| si32 | off standard | 38.357 | 7.067e-10 | 2.234e-06 | 3.469e-07 | 15 | 67 | 2806 |
| si32 | on standard | 24.002 | 1.640e-09 | 4.218e-06 | 6.634e-07 | 13 | 60 | 2380 |

Density differences are real-space `density.bin` against the same code's tight reference on the same FFT grid (electron counts 32/64/128 agree).

## QE ladder (np=1, 5-run medians; abs dE against QE's own conv_thr 1e-12 result)

| Case | conv_thr (Ry) | wall (s) | abs dE (Ha) | SCF | estimated accuracy (Ry) |
| --- | ---: | ---: | ---: | ---: | ---: |
| si8 | 1e-8 / 1e-9 / 1e-10 / 1e-12 | 0.467 / 0.513 / 0.551 / 0.606 | 9.230e-10 / 1.659e-10 / 2.613e-11 / 0 | 12 / 14 / 15 / 18 | 6.6e-9 / 2.7e-10 / 5.1e-11 / 2.7e-13 |
| si16 | 1e-8 / 1e-9 / 1e-10 / 1e-12 | 1.413 / 1.685 / 1.826 / 1.974 | 2.699e-10 / 2.699e-10 / 5.436e-12 / 0 | 10 / 12 / 13 / 15 | 6.3e-9 / 1.1e-10 / 2.6e-11 / 1.8e-13 |
| si32 | 1e-8 / 1e-9 / 1e-10 / 1e-12 | 5.031 / 6.676 / 6.689 / 7.867 | 2.588e-10 / 2.588e-10 / 1.820e-11 / 0 | 10 / 13 / 13 / 17 | 4.1e-9 / 6.9e-11 / 6.9e-11 / 7.4e-13 |

QE np=8: si8 conv_thr 1e-8/1e-9/1e-10 0.252/0.258/0.273 s are below 0.3 s and INVALID; si8 1e-12 0.338 s, si16 0.556 to 0.698 s and si32 1.693 to 2.263 s are valid.

## Equal-accuracy pairs, one core each

Delivered error of about 1e-9 Ha (BAANDI off standard versus QE conv_thr 1e-8): si8 1.308 s (1.55e-09) vs 0.467 s (9.23e-10): QE 2.80x faster; si16 5.866 s (1.27e-09) vs 1.413 s (2.70e-10): 4.15x; si32 38.357 s (7.07e-10) vs 5.031 s (2.59e-10): 7.62x. With gamma_real on: 2.78x, 3.68x, 4.77x (but the si32 on-standard error 1.64e-09 is 6.3x worse than QE's).
Delivered error of about 1e-11 Ha (BAANDI off strict versus QE conv_thr 1e-10): si8 3.764 s (1.35e-13) vs 0.551 s (2.61e-11): 6.83x; si16 15.124 s (1.71e-11) vs 1.826 s (5.44e-12): 8.28x; si32 119.223 s (2.94e-11) vs 6.689 s (1.82e-11): 17.82x.
Cleanest pair (accuracy practically equal): si16 on strict 13.319 s / 5.379e-12 Ha versus QE conv_thr 1e-10 1.826 s / 5.436e-12 Ha: QE 7.29x faster (BAANDI 3053 H vectors; QE 148 h_psi block calls and 1612 FFT calls).
Gap versus size (off strict): 6.83x, 8.28x, 17.82x for si8, si16, si32: the gap grows with the system size.

Tool limitation: `scripts/compare_materials_qe.py --qe-density` rejects the `K_POINTS gamma` charge-density file (`unsupported gamma_only: expected 0`, half-sphere G set), so no code-to-code density comparison exists for the real path; density differences above are against BAANDI's own tight reference. A complex-path (`K_POINTS automatic 1 1 1`) comparison remains possible.

# dev08 re-measurement after P30A (pool-distributed pair FFT) and P32A (real algebra kernels), main 72574f9 (relayed, 2026-11-07)

Relayed by the dft5090/master session. main `72574f9`, pure native `d9b64928cbe9ffbb84324677bf812823b29b94318dcbe6680ecced40fb03035b`, rustc 1.99.1, n=10 per row, `taskset -c 0` (T=1) or `0-7` (T=8), start load average 1.44 (earlier measurements 0.16 to 0.44: a confounder). Compared with `8894b0b` (results identical to `b6a5416`).

**Correction of an earlier statement:** the real path (gamma_real on) is NOT bit-identical to the previous implementation after P32A, because the GEMM blocking and slice loops change the summation order; Davidson then converges in a different number of steps. Only the complex path (gamma_real off) is bit-identical (8 rows, energy difference 0, same SCF/Davidson/H counts). Real-path changes: si8 strict energy difference 8.53e-14 Ha (Davidson steps 319 to 258, H vectors 2125 to 2014); si16 strict 4.59e-12 Ha (216 to 185, 3053 to 2895); si8 standard 7.85e-11 Ha (SCF 17 to 15, 85 to 78, 864 to 795); si16 standard 4.41e-10 Ha (81 to 70, 1424 to 1442). All within the strict 1e-10 Ha gate or the standard delivered-accuracy class; future regression checks must compare the real path against the new numbers, not the old ones.

| Case | T | row | wall median (s) | min | max | SCF | Davidson | H vectors | energy (Ha) |
| --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| si8 | 1 | off strict | 3.861 | 3.831 | 3.927 | 15 | 279 | 2044 | -31.321803601031636 |
| si8 | 1 | on strict | 3.134 | 3.126 | 3.178 | 15 | 258 | 2014 | -31.321803601031739 |
| si8 | 1 | off standard | 1.337 | 1.331 | 1.380 | 15 | 77 | 772 | -31.321803599478528 |
| si8 | 1 | on standard | 1.352 | 1.348 | 1.392 | 15 | 78 | 795 | -31.321803599784658 |
| si8 | 8 | off strict | 2.513 | 2.468 | 2.562 | 15 | 279 | 2044 | -31.321803601031636 |
| si8 | 8 | on strict | 2.525 | 2.508 | 2.613 | 15 | 258 | 2014 | -31.321803601031739 |
| si8 | 8 | off standard | 0.943 | 0.910 | 0.967 | 15 | 77 | 772 | -31.321803599478528 |
| si8 | 8 | on standard | 1.176 | 1.158 | 1.210 | 15 | 78 | 795 | -31.321803599784658 |
| si16 | 1 | off strict | 15.392 | 15.286 | 15.517 | 12 | 182 | 2700 | -63.129244421600596 |
| si16 | 1 | on strict | 11.694 | 11.669 | 11.740 | 13 | 185 | 2895 | -63.129244421616882 |
| si16 | 1 | off standard | 5.949 | 5.931 | 6.009 | 14 | 76 | 1350 | -63.129244420349295 |
| si16 | 1 | on standard | 5.423 | 5.412 | 5.464 | 15 | 70 | 1442 | -63.129244421079889 |
| si16 | 8 | off strict | 7.552 | 7.494 | 7.582 | 12 | 182 | 2700 | -63.129244421600596 |
| si16 | 8 | on strict | 7.576 | 7.524 | 7.594 | 13 | 185 | 2895 | -63.129244421616882 |
| si16 | 8 | off standard | 3.186 | 3.165 | 3.227 | 14 | 76 | 1350 | -63.129244420349295 |
| si16 | 8 | on standard | 4.028 | 4.016 | 4.064 | 15 | 70 | 1442 | -63.129244421079889 |

Same trajectory (SCF equal) between on and off: si8 both profiles (15 = 15; Davidson and H counts differ); si16 strict 12 versus 13 and standard 14 versus 15 differ.
Real path on versus off, wall ratio, previous (before P30A/P32A) to now: T=1: si8 strict -6.4 % to -18.8 %, si8 standard -0.6 % to +1.1 %, si16 strict -11.9 % to -24.0 %, si16 standard -11.3 % to -8.8 %; T=8: si8 strict +38.6 % to +0.5 %, si8 standard +41.1 % to +24.7 %, si16 strict +45.2 % to +0.3 %, si16 standard +42.2 % to +26.4 %.
Thread scaling T=1 to T=8, previous to now (off / on): si8 strict 1.53 to 1.54 / 1.03 to 1.24; si8 standard 1.42 to 1.42 / 1.00 to 1.15; si16 strict 2.03 to 2.04 / 1.23 to 1.54; si16 standard 1.87 to 1.87 / 1.16 to 1.35.
The T=8 inversion has practically disappeared for strict (+0.5 %, +0.3 %: a tie) and shrunk but remains for standard (+24.7 %, +26.4 %). The Mac observation (off 1.19 s, on 0.82 s at T=8) does not carry over: on dev08 strict T=8 is a tie.
Open: the complex path (off) is 1.3 to 2.6 percent slower in all eight rows with identical work counts; the load difference is a confounder and an alternating-run control (8894b0b versus 72574f9, 10 each, same session) is running. si32 n=3 is running.

## Addendum (relayed 2026-11-07): complex-path slowdown bisect and si32 re-measurement on main 72574f9

Complex path (gamma_real off) slowdown is real and reproducible on dev08: alternating runs (same session, 10 each, si8/si16 strict T=1): si8 new/old +2.55 % (3.811 versus 3.716 s median, no interval overlap), si16 +1.27 %. Bisect (si8 off strict T=1, five commits round-robin 10 rounds): 8894b0b 3.726 s, 74f4855 +0.07 %, 59102b6 -0.20 %, 25e0ff4 (merge of P30A) +1.98 %, 72574f9 (P32A) +2.78 %; paired steps 59102b6 to 25e0ff4 +2.06 % (10 of 10 rounds slower) and 25e0ff4 to 72574f9 +1.02 % (9 of 10). Orchestrator check: `git diff 59102b6 25e0ff4 -- crates` touches only real/pair files (gamma_real.rs, fused_pair.rs, real_gamma.rs and tests); no complex-path source changed. On the Mac (Apple Silicon) alternating runs of the 8894b0b and current binaries give new/old = 0.999 (10 each, si8 off strict T=1: 2.416 s versus 2.413 s). The effect is therefore x86-specific code generation or layout, not a logic change; results stay bit-identical. Status: tracked, about 2.8 percent on the complex path; a codegen-units=1 control on dev08 is optional.

si32 (main 72574f9, n=3, wall median s, SCF/Davidson/H vectors): T=1 off strict 120.288 (13/352/5946), on strict 103.603 (14/426/6481), off standard 38.771 (15/67/2806), on standard 22.562 (13/63/2430); T=8 off strict 51.149, on strict 57.958, off standard 16.918, on standard 15.202. Energies: off strict -126.703998096079, on strict -126.703998096098, off standard -126.703998095401, on standard -126.703998094163. No si32 row has equal SCF counts for on and off. Real path on versus off, previous to now: strict T=1 -18.2 % to -13.9 % (more work: Davidson 359 to 426), strict T=8 +34.9 % to +13.3 %, standard T=1 -37.4 % to -41.8 %, standard T=8 +18.0 % to -10.1 %. Thread scaling T=1 to T=8 for on: strict 1.42 to 1.79, standard 1.21 to 1.48 (off unchanged 2.35 and 2.29).

# dev08 P29A perturbation harness: standard-profile SCF distribution (relayed, 2026-11-07)

Relayed by the dft5090/master session. main `72574f9`, `P29A_SI_UPF=<Si.pz-vbc.UPF> P29A_CASES=si8-gamma,si16-gamma`, `RUSTFLAGS="-C target-cpu=native"`, `cargo test --release -p dft-engine --test standard_stall_trajectory sweep -- --ignored --nocapture --test-threads=1`, exit 0, 149 s; log `/scr/mino/rustdft/p29a/sweep.log` on dev08. Runs 0 and 1 are strict references (tolerance multipliers 1.0 and 0.1), runs 2 to 11 are standard with ten 1e-10 relative start-density perturbations; all 24 runs converged, no failures; standard results are uncertified.

si8-gamma standard SCF iterations (runs 2 to 11): 25, 23, 25, 25, 24, 20, 21, 25, 34, 21 (min 20, median 24.5, max 34; Mac: 20 / 25.5 / 32). Davidson steps 82 to 120, H work 959 to 1539 (median 1138); strict run 0: 18 SCF, 416 steps, 2162 H; run 1: 21, 557, 2634.
si16-gamma standard SCF: 20, 19, 18, 24, 19, 23, 22, 20, 20, 17 (min 17, median 20, max 24); H work 1491 to 2067 (median 1743.5); strict run 0: 18, 262, 3730; run 1: 21, 308, 4321.
H-work gain, strict run 0 over standard median: si8 1.90x, si16 2.14x (against the tighter run 1: 2.31x and 2.48x); consistent with the Mac study (Si8 median 1.78x).
Standard delivered error against strict run 1: si8 worst 3.06e-09, best 2.6e-10 Ha; si16 worst 1.74e-09, best 4.2e-10 Ha.

Resolution of the earlier question: the dev08 "always 15 SCF iterations" of the batch-7 B measurements is the single unperturbed trajectory (deterministic); with the 1e-10 start perturbation the dev08 distribution (20 / 24.5 / 34) matches the Mac (20 / 25.5 / 32). The earlier suspicion of a platform or build difference was wrong. Consequence: standard-profile wall times from the unperturbed batch-7 tables (si8 1.31 s, si16 5.9 s at T=1) belong to a short trajectory; a median trajectory costs about 1.5 times more H work for si8 (1138 versus 772 H vectors).
