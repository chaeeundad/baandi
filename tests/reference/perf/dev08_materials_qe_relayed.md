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

## Addendum (relayed 2026-11-07): codegen-units=1 control for the complex-path slowdown

si8 gamma_real off strict T=1, `taskset -c 0`, round-robin 10 rounds, both commits built with `RUSTFLAGS="-C target-cpu=native -C codegen-units=1"` (binaries cgu1-59102b6 `85a9d51578a8d0e6...`, cgu1-25e0ff4 `c5b653d3dc2631a4...`): 59102b6 median 3.758 s (min 3.726, max 3.832), 25e0ff4 3.761 s (3.731, 3.816); paired per-round difference median -0.21 % (min -2.65 %, max +2.41 %), slower in 4 of 10 rounds. The same pair with the default codegen-units (16): paired median +2.06 %, slower in 10 of 10 rounds. Conclusion: the +2.06 % at the P30A merge is a code-placement/codegen-unit effect, not a logic or path regression (results are bit-identical). Side observation (different sessions, weak): codegen-units=1 builds were not faster than the default build (59102b6 3.719 s default versus 3.758 s cgu1; 25e0ff4 3.800 versus 3.761 s), so "codegen-units=1 is always faster" does not hold for this code base on this CPU. A same-session default-versus-cgu1 alternation (about 3 minutes) is available if wanted.

# dev08 (5) k-point symmetry: BAANDI auto versus off, QE sym versus nosym (relayed, 2026-11-07)

Relayed by the dft5090/master session. main `72574f9`, pure native `d9b64928cbe9ffbb84324677bf812823b29b94318dcbe6680ecced40fb03035b`; inputs from `bench_scf.case_input` plus one `[symmetry] mode` line; n=5; `taskset -c 0 / 0-7 / 0-15`. QE 7.6 `pw.x` sha256 `07c7cd7935f42b15157500ab719a5c3948825d8ce8a4487d992d7501c659658f`, `disk_io=nowf`, conv_thr 1e-10 Ry, OMP/MKL 1.

BAANDI k222 (8 full to 3 irreducible, 48 mesh-preserving operations): T=1 off 0.998 s, auto 0.628 s (-37.1 %); T=8 off 0.299 s*, auto 0.390 s (+30.4 %); T=16 off 0.287 s*, auto 0.390 s (+35.9 %) (* below 0.3 s: invalid under the timer rule for QE only; shown for BAANDI as relayed). SCF 10 both; Davidson 408 to 149, H vectors 2543 to 934; E off -7.831230057121136, auto -7.831230076724104, difference 1.960e-08 Ha.
BAANDI k444s (64 full to 10 irreducible, 12 mesh-preserving operations): T=1 off 4.136 s, auto 1.662 s (-59.8 %); T=8 0.868 / 0.634 s (-27.0 %); T=16 0.631 / 0.505 s (-20.0 %). SCF 10 both; Davidson 1686 to 524, H vectors 10290 to 3201; E off -7.925976126648090, auto -7.925976134331592, difference **7.683e-09 Ha**; with tightened criteria (energy 1e-12, density 1e-10, T=8) the same difference (off -7.925976126648077, auto -7.925976134331563, 13 SCF each), i.e. not a convergence-floor artifact.
QE k444s shifted 4x4x4 (1 1 1): nosym 32 k points, 6 SCF, -7.925978690105397 Ha; sym 10 k points, 48 symmetry operations with inversion (24 fractional translations), 7 SCF, -7.925994229219341 Ha. Wall (5-run median): nosym np1 1.4576 s, np8 0.7779 s; sym np1 0.6225 s, np8 0.3575 s; QE symmetry gain -57.3 % (np1), -54.1 % (np8). QE sym minus nosym = **1.554e-05 Ha**.
Control (unshifted 4x4x4, np=8): nosym 36 points -7.919062122458864 Ha, sym 8 points -7.919062122460780 Ha, difference 1.9e-12 Ha: the QE difference appears only with the shifted mesh.
Cross-code differences: BAANDI off (64 points) versus QE nosym (32 points) 2.563e-06 Ha (1.28e-06 per atom, equal to the known systematic offset 1.27e-06 per atom); BAANDI auto versus QE sym 1.809e-05 Ha (9.05e-06 per atom).
Observation: with symmetry on, thread scaling worsens because the outer parallel width is limited by the number of irreducible k points (k222: 3 points, slower than off at T=8 and T=16).

Orchestrator check (Mac, same bench case `k444s`, default settings, T=1): off -7.9259761266480915, auto -7.925976134331593 Ha: the same 7.68e-09 difference as on dev08. The P19A validation recorded 1.2434e-14 Ha for the same geometry with its own tightened settings (docs/validation/symmetry.md); the two settings differ: the bench case `off` run uses a 25x25x25 FFT grid while `auto` selects 32x32x32 (a grid that can represent diamond's quarter translation as a permutation; P19A forced 32^3 on both sides, see docs/validation/symmetry.md), so the 7.68e-09 Ha is a grid-discretization difference between the two runs (checked: results.json `fft_grid` 25^3 versus 32^3), not a symmetry inconsistency. Equal-grid off/auto comparisons are the correct test of the symmetry reduction.

## Correction and additions (relayed 2026-11-07): k444s symmetry with equal FFT grid, QE k-point reduction check

BAANDI k444s with `[basis] fft_grid = [32, 32, 32]` forced on both sides (T=1, n=5, `taskset -c 0`): off 4.893 s (4.880 to 4.901), 10 SCF, 1686 Davidson steps, 10290 H vectors, E -7.925976134331574 Ha; auto 1.672 s (1.669 to 1.762), 10 SCF, 524, 3201, E -7.925976134331592 Ha. Energy difference **1.78e-14 Ha** (P19A: 1.2e-14), same trajectory; auto/off wall -65.8 %. (The earlier 7.683e-09 Ha was the 25^3 versus 32^3 grid difference: off [25,25,25] -7.925976126648090, auto [32,32,32] -7.925976134331592; the tight-criteria runs kept their own grids.) k222 with equal grid not yet re-measured.
QE k444s shifted mesh k points: sym 10 points, weights 0.0625 (x2), 0.1875 (x6), 0.375 (x2) summing to 2.0; nosym 32 points all 0.0625. Checks by the relaying session: weights times 16 give 1,3,3,3,3,6,6,3,1,3 (sum 32 = nosym point count); all 10 sym points occur in the nosym list (nosym indices 1,2,3,4,6,7,8,11,22,23): the sym set is a complete orbit decomposition of the nosym set, so QE does not sample a different set. QE FFT dimensions (25,25,25), 5985 G-vectors for shifted and unshifted, sym and nosym; 48 symmetry operations with inversion (24 fractional translations) in both sym runs. The converged states differ: k(1) eigenvalues (eV) sym -5.6806 4.4485 5.7545 5.7545 8.5169 9.0270 9.0270 10.6627, nosym -5.6805 4.4516 5.7551 5.7551 8.5037 9.0293 9.0293 10.6676 (maximum difference 0.0132 eV, band 5); both met their convergence criteria (sym 2.3e-12 Ry, nosym 6.6e-11 Ry). Unshifted control (same grid, same 48 operations): sym minus nosym 1.9e-12 Ha. QE shifted sym minus nosym 1.554e-05 Ha therefore remains UNRESOLVED; excluded: different k sampling, FFT grid difference between runs, convergence. Orchestrator hypothesis (untested): the 25^3 grid cannot represent diamond's quarter translation as a real-space permutation, which matters for the shifted run; control: QE with nr1=nr2=nr3=32 for both sym and nosym.
Fair-accuracy reading (orchestrator arithmetic): BAANDI auto equals BAANDI off to 1.8e-14 Ha, while QE sym differs from QE nosym by 1.55e-05 Ha; at full-mesh accuracy the relevant QE time is nosym (np1 1.4576 s) against BAANDI auto 1.672 s (1.15x) rather than QE sym 0.6225 s (BAANDI 2.7x slower).

# dev08 (6) P15B band-blocked nonlocal: isolation comparison (relayed, 2026-11-07)

Relayed by the dft5090/master session. Isolation base `ce64b57` (parent of P15B; commit hash embedded in the binary), candidate `963f94c` built fresh (the 5 November `main-*` binaries were main at that time, not 963f94c, and were not used). Variants pure/MKL x default/native; pure native sha256 prefixes: p15b-pure-def `30b7e43361fcb744`, p15b-pure-nat `49b84bd4cc8c2874`, p15b-mkl-def `a00776a9726c693d`, p15b-mkl-nat `786da1c893857756` (MKL link needs `/opt/intel/oneapi/mkl/2026.1/env/vars.sh` sourced, otherwise `undefined symbol: DftiErrorMessage`). Raw data on dev08 `/scr/mino/rustdft/p15b/`.
Accuracy (T=1, 4 variants x 2 cases): all 8 rows bit-identical (dE 0; si8 grid 36^3: pure 15 SCF / 279 Davidson / 2044 H, MKL 15 / 302 / 2139; k444s grid 25^3: 10 / 1686 / 10290).
Time (pure native, round-robin alternation, 10 rounds): si8 T=1 base 4.126 s versus candidate 4.376 s (+6.05 %); T=4 2.804 versus 3.041 (+8.46 %); T=8 2.878 versus 3.085 (+7.20 %); T=16 3.169 versus 3.409 (+7.55 %); k444s T=1 5.008 versus 5.246 (+4.74 %); T=8 1.068 versus 1.098 (+2.77 %). The candidate was faster in 0 of 10 paired rounds in every one of the six configurations; intervals do not overlap.
Profile (si8 T=1 pure native, eu-stack 20 ms, 15 runs): nonlocal bucket 12.4 % (238 samples) for the base versus 17.0 % (349) for the candidate; top symbol `hamiltonian::add` 10.69 % (207) versus `apply_nonlocal_group` 12.87 % (266) plus `add` 2.42 % (50); other buckets moved by less than one percentage point (local FFT 29.0 to 27.9, algebra 19.5 to 18.5, unattributable GEMM 18.2 to 18.3). Nonlocal time estimate 0.512 s versus 0.744 s: the nonlocal increase (0.232 s) explains 93 % of the total increase (0.250 s).
Reading: the band-blocked traversal reads projectors 3.2x less often but is slower on the Xeon Gold 6230 in every configuration; the Mac observation of the same commit (ledger A58: nonlocal -5.7 % per vector at T=1, but whole H +14.7 % at T=8) was already not clearly positive. Decision (orchestrator): the grouped path must not be the default.

## Resolution (orchestrator, 2026-11-07): reference-code shifted-mesh symmetry shift explained by ticket P34A

BAANDI's research option `[symmetry] density_group = "crystal"` (k-point reduction unchanged, density averaged over all 48 crystal operations) gives, on the shifted Si 4x4x4 mesh with tightened convergence, E(crystal) minus E(off) = -1.5538989465824216e-05 Ha; the reference code's symmetry-on minus symmetry-off energy is -1.554e-05 Ha (-1.55390e-05 at grid 32): agreement to 1.0e-09 Ha (0.0065 %), same sign. The unshifted mesh gives exactly zero difference for crystal versus mesh, and the unsymmetrized mesh density violates the 36 non-mesh operations by about 2.1e-3 electrons/Bohr^3. Interpretation: the reference code averages the density over the full crystal group although only 12 operations map the shifted mesh onto itself, which amounts to a different effective k sampling (the union of rotated meshes); it is not an error of either code, and neither sampling is claimed to be more correct. Details: docs/validation/symmetry-full-group.md. For equal-accuracy speed comparisons on shifted meshes, BAANDI with `density_group = "crystal"` is the like-for-like counterpart of the reference code's symmetry-on run.

# dev08 P33A comparison (relayed, 2026-11-07): real Davidson reusable workspace is -9.6 to -10.7 percent on x86

Relayed by the dft5090/master session. Base `73f59e3`, candidate `56c4127` (branch `task/P33A-REAL-WORKSPACE-PARALLEL`), both `RUSTFLAGS="-C target-cpu=native"` (plain sha256 prefixes base `0982b8b22f6341ac`, candidate `762d897c2acc6787`; frame-pointer variants `dd4ecbf2e289cf67`, `40563e1646c66630`). gamma_real on, strict, pure native, round-robin alternation, 10 rounds: si8 T=1 3.132 s [3.115, 3.184] versus 2.831 s [2.818, 2.837] (-9.61 %, paired -9.56 %, candidate faster in 10 of 10 rounds); si8 T=8 2.552 versus 2.278 s (-10.72 %, paired -10.43 %, 10/10); si16 T=1 11.778 versus 10.649 s (-9.59 %, paired -9.61 %, 10/10). Intervals separate completely. Accuracy: bit-identical in all three (si8 T=1/T=8 -31.321803601031739 Ha, 15 SCF, 258 Davidson, 2014 H; si16 -63.129244421616882, 13/185/2895).
Profile (si8 T=1 on strict, frame-pointer build, 15 runs): bucket shares of working samples base versus candidate: local FFT 32.5 vs 34.2 %, nonlocal 2.6* vs 5.0*, subspace algebra 46.4 vs 47.0 %, residual 0.8* vs 0.5*, other 17.6 vs 13.4 % (working samples 1177 vs 1101; `*` under 30 samples). memcpy self share 6.29 % (75 samples) versus 2.33 % (26 samples: below 30, direction only); memset 5.45 % (65) versus 6.98 % (78); sum 11.74 % versus 9.31 %.
Reading: a change within noise on the Mac (1.1845 versus 1.1826 s) is a clear gain on the Xeon, in line with the 12.6 % memset plus memcpy share seen in the dev08 profile. Orchestrator decision: merged to main (gates: bit-identical, 30 of 30 rounds consistent, cause symbol matches the intent).

# dev08 k444s: BAANDI auto+crystal versus QE sym, equal grid (relayed, 2026-11-07)

Relayed by the dft5090/master session. main `b022b61` (P33A included), pure native sha256 `c4fba2eca3cece046fd9dfea4f69525bc8d5cfe7e1b210262a4613d56b2879e1`, rustc 1.99.1; k444s (shifted 4x4x4, 2-atom Si cell, 30 Ry), `[basis] fft_grid = [32,32,32]` forced, T=1, `taskset -c 0`, n=5. tight = energy_tolerance_ha 1e-12, density_residual_tolerance 1e-10.

| BAANDI mode | tolerance | wall median (s) | min | max | SCF | Davidson | H vectors | energy (Ha) |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | --- |
| off | default | 4.790 | 4.765 | 4.802 | 10 | 1686 | 10290 | -7.925976134331574 |
| auto (mesh) | default | 1.688 | 1.683 | 1.763 | 10 | 524 | 3201 | -7.925976134331592 |
| auto + density_group crystal | default | 1.698 | 1.672 | 1.707 | 10 | 532 | 3256 | -7.925991673321009 |
| off | tight | 5.348 | 5.343 | 5.435 | 13 | 1889 | 12190 | -7.925976134331593 |
| auto (mesh) | tight | 1.870 | 1.838 | 1.882 | 13 | 586 | 3794 | -7.925976134331563 |
| auto + crystal | tight | 1.836 | 1.830 | 1.927 | 12 | 582 | 3686 | -7.925991673321009 |

E(crystal) minus E(off) = -1.5538989435e-05 Ha (independent dev08 reproduction of the P34A value -1.5538989e-05; crystal energies are identical at both tolerances to the last digit).
QE sym (same 32^3 grid, 10 irreducible points, np=1): -7.925994236779095 Ha, 7 SCF, accuracy 2.3e-12 Ry, wall median 0.805 s (n=3, min 0.801, max 0.875). BAANDI auto+crystal minus QE sym = 2.563458e-06 Ha = 1.2817e-06 Ha per atom, equal to the known convention offset of 1.27e-06 Ha per atom (batch-7 C). The earlier mismatched pair (auto mesh versus QE sym) was 1.810e-05 Ha, seven times the offset.
Speed at equal grid, one core: QE sym 0.805 s; BAANDI auto+crystal 1.698 s (default tolerance, 2.11x slower) or 1.836 s (tight, 2.28x); BAANDI auto mesh 1.688 s (2.10x). Targets differ (QE conv_thr 1e-10 Ry versus BAANDI default energy tolerance 1e-8 Ha; tight is closer to QE). QE sym on its 25^3 grid takes 0.6225 s (BAANDI 3.1x), but the equal-grid comparison is the fair one. The crystal option costs +0.6 % here (Davidson 524 to 532). This case (2 atoms, k points reduced 64 to 10) is relatively favourable to BAANDI compared with the Gamma cases (2.8 to 17.8x).

# dev08 batch-7 D: latest main b022b61 (P33A included), real path on/off, equal-accuracy gap to QE (relayed, 2026-11-07)

Relayed by the dft5090/master session. main `b022b61`, pure native `c4fba2eca3cece046fd9dfea4f69525bc8d5cfe7e1b210262a4613d56b2879e1`, rustc 1.99.1; Si.pz-vbc.UPF `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`; QE 7.6 `pw.x` `07c7cd7935f42b15157500ab719a5c3948825d8ce8a4487d992d7501c659658f`; si8, si16, si32 as defined earlier (Gamma, 30 Ry, lda-pz, Pulay beta 0.5); n=10 (si32 n=3); `taskset -c 0` or `0-7`; QE OMP/MKL 1, `disk_io=nowf`. Rows: wall median (s) with SCF / Davidson / H vectors (identical inside a row, and between T=1 and T=8).

| Case | T | off strict | on strict | off standard | on standard |
| --- | ---: | --- | --- | --- | --- |
| si8 | 1 | 3.797 (15/279/2044) | 2.961 (15/258/2014) | 1.307 (15/77/772) | 1.323 (15/78/795) |
| si8 | 8 | 2.465 | 2.464 | 0.920 | 1.200 |
| si16 | 1 | 15.198 (12/182/2700) | 11.115 (13/185/2895) | 5.850 (14/76/1350) | 5.294 (15/70/1442) |
| si16 | 8 | 7.450 | 6.975 | 3.129 | 3.802 |
| si32 | 1 | 119.525 (13/352/5946) | 97.046 (14/426/6481) | 38.326 (15/67/2806) | 22.140 (13/63/2430) |
| si32 | 8 | 50.552 | 50.284 | 16.857 | 14.828 |

Same trajectory (equal SCF) between on and off: si8 all rows; not si16 (strict 12 vs 13, standard 14 vs 15) and not si32 (strict 13 vs 14, standard 15 vs 13).
Real path on versus off wall ratio, b6a5416 / 72574f9 / b022b61: si8 strict T=1 -6.4 / -18.8 / -22.0 %, T=8 +38.6 / +0.5 / -0.04 %; si8 standard T=1 -0.6 / +1.1 / +1.2 %, T=8 +41.1 / +24.7 / +30.4 %; si16 strict T=1 -11.9 / -24.0 / -26.9 %, T=8 +45.2 / +0.3 / -6.4 %; si16 standard T=1 -11.3 / -8.8 / -9.5 %, T=8 +42.2 / +26.4 / +21.5 %; si32 strict T=1 -18.2 / -13.9 / -18.8 %, T=8 +34.9 / +13.3 / -0.5 %; si32 standard T=1 -37.4 / -41.8 / -42.2 %, T=8 +18.0 / -10.1 / -12.0 %. The strict T=8 inversion has disappeared at all three sizes; the standard T=8 inversion remains for si8 and si16 (trajectories differ, so the standard rows are not same-trajectory comparisons).
Thread scaling T=1 to T=8 (off / on), b6a5416 to b022b61: si8 strict 1.53 to 1.54 / 1.03 to 1.20; si8 standard 1.42 / 1.00 to 1.10; si16 strict 2.03 to 2.04 / 1.23 to 1.59; si16 standard 1.87 / 1.16 to 1.39; si32 strict 2.35 to 2.36 / 1.42 to 1.93; si32 standard 2.28 to 2.27 / 1.21 to 1.49.

## Equal-accuracy gap to QE real gamma (BAANDI T=1 versus QE np=1; each code against its own tight reference; BAANDI tight references unchanged because the complex path is bit-identical)

| Case | BAANDI best setting | wall (s) | abs dE (Ha) | QE conv_thr (Ry) | QE wall (s) | QE abs dE (Ha) | QE faster by |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| si8 | T=1 off standard | 1.307 | 1.553e-09 | 1e-8 | 0.467 | 9.230e-10 | 2.80x |
| si8 | T=1 on strict | 2.961 | 3.197e-14 | 1e-10 | 0.551 | 2.613e-11 | 5.37x |
| si16 | T=1 on standard | 5.294 | 5.378e-10 | 1e-8 | 1.413 | 2.699e-10 | 3.75x |
| si16 | T=1 on strict | 11.115 | 7.887e-13 | 1e-10 | 1.826 | 5.436e-12 | 6.09x |
| si32 | T=1 on standard | 22.140 | 1.945e-09 | 1e-8 | 5.031 | 2.588e-10 | 4.40x |
| si32 | T=1 on strict | 97.046 | 1.026e-11 | 1e-10 | 6.689 | 1.820e-11 | 14.51x |

Gap change against the previous measurement (strict class / standard class): si8 6.83x to 5.37x / 2.80x to 2.80x; si16 8.28x to 6.09x / 4.15x to 3.75x; si32 17.82x to 14.51x / 7.62x to 4.40x. The growth of the gap with system size persists. Standard rows are single unperturbed trajectories (the median trajectory costs about 1.5 times more H work, see the P29A/P29B studies).
Open items listed by the relaying session: the T=8 standard inversion (si8 +30.4 %, si16 +21.5 %); only about two of 17 pool threads busy at T=8 (unchanged by P30A; band-direction width at a single k point to be confirmed); si32 strict real path needs more work than off (SCF 13 to 14, Davidson 352 to 426) for unknown reasons.

# dev08 Q1 (a): QE parameter sweep (relayed, 2026-11-07)

Relayed by the dft5090/master session. QE 7.6 `pw.x` sha256 `07c7cd7935f42b15157500ab719a5c3948825d8ce8a4487d992d7501c659658f`, OMP/MKL threads 1, `disk_io=nowf`, `verbosity='high'`, real gamma runs with `nosym/noinv=.true.` on BAANDI's FFT grids (`--fft-grid ours`: si8 nbnd 20 nr 36, si16 36 / 72x36x36, si32 68 / 72x72x36), load below 0.30 before each run, 5-run medians; k444s = shifted 4x4x4, symmetry on, 10 irreducible points, conv_thr 1e-10 Ry, nr 25. Errors against the QE tight reference (conv_thr 1e-12) of each case. Only a4_k444s_nk2 (0.293 s) falls below 0.3 s and is marked invalid.
a1 nbnd = occupied versus current (wall med s, SCF, h_psi calls, fftw calls, |dE| Ha): si8 conv 1e-8: nb16 0.468, 16, 86, 938, 2.419e-09 versus nb20 0.468, 12, 94, 964, 2.038e-09; si8 1e-10: nb16 0.634 (25 SCF) versus nb20 0.552 (15); si16 1e-8: nb32 1.133 (12 SCF, 55 h_psi, 1432 fftw, 2.890e-09) versus nb36 1.418 (10, 148, 1612, 2.081e-09); si16 1e-10: 1.358 (15 SCF, 1.549e-11) versus 1.841 (13, 7.482e-12); si32 1e-8: nb64 3.749 (10, 45, 2406, 1.115e-09) versus nb68 4.988 (10, 123, 2768, 1.701e-09); si32 1e-10: 5.332 (15, 2.660e-11) versus 6.650 (13, 2.490e-11). Effect of occupied-only bands: si8 0 % / +14.9 %, si16 -20.1 % / -26.2 %, si32 -24.8 % / -19.8 % (1e-8 / 1e-10); h_psi block calls fall strongly but fftw calls only by about 13 % (si16 1612 to 1432, si32 2768 to 2406).
a2 diago_david_ndim (default 2 in this build; identical to the nd2 run): si8 nd2 0.468 s (12 SCF, 94 h_psi), nd3 0.445, nd4 0.456; si16 nd2 1.411, nd3 1.347, nd4 1.522; si32 nd2 5.037, nd3 5.261, nd4 5.633 (worse); k444s default 0.625, nd2 0.623, nd3 0.590, nd4 0.592 (h_psi 397, 340, 336).
a3 mixing_beta 0.5 versus 0.7: si8 0.480 s (13 SCF) vs 0.467 (12); si16 1.279 (9) vs 1.413 (10); si32 4.735 (10) vs 5.029 (10); k444s 0.592 (7) vs 0.621 (7).
a4 k444s np=8: plane-wave split nk1 0.358 s (7 SCF, h_psi 393, fftw 4912); nk2 0.293 (invalid, below 0.3 s); nk4 0.306 (124, 1484); nk8 0.335 (84, 992). Best valid: nk4, -14.5 % against nk1.
a5 forces and stress cost: si8 0.467 s to 0.494 s (+5.8 %; forces timer 0.01 s, stress 0.02 s), k444s 0.621 to 0.642 (+3.4 %); fft calls 38 to 43 and 24 to 29.
a6 verbosity high counts (default setting): si8 h_psi 94, fftw 964, fft 38; si16 148, 1612, 33; si32 123, 2768, 33; k444s (complex path, cegterg) 397, 4932, 24. fftw time per call: si8 135 us, si16 341 us, si32 about 795 us.
Raw logs on dev08: `/scr/mino/rustdft/q1a/`, summary `/scr/mino/rustdft/q1a_run.log`.
Orchestrator note: the plan document `docs/performance/qe-parity-plan.md` was committed but not yet pushed when the relaying session looked for it; pushed afterwards.

# dev08 re-measurement on main 42240e5 (P36A, P36B, P37A included) and Q1 (b) (relayed, 2026-11-07)

Relayed by the dft5090/master session. main `42240e5`, pure native plain sha256 `f089848b30a22be7ccf141357511ca7268094ca2c20b03c20f370531ee3c34de`, frame-pointer variant `d4d8f75dca8d3a1e3d9cd80d265973a39babddf063ed88da4261ffbb3cc5071b`; comparison base `b022b61` sha256 `c4fba2eca3cece046fd9dfea4f69525bc8d5cfe7e1b210262a4613d56b2879e1`. Round-robin alternation, 10 rounds; all rows bit-identical (energy, SCF, Davidson, H vectors).

Real path (gamma_real on), old b022b61 median versus new 42240e5 median (paired change, candidate faster in all 10 of 10 rounds in every row): si8 strict T=1 2.960 to 2.609 s (-11.8 %), T=8 2.452 to 1.941 (-20.9 %); si16 strict T=1 11.085 to 9.920 (-10.5 %), T=8 6.958 to 5.489 (-21.1 %); si8 standard T=1 1.324 to 1.060 (-19.9 %), T=8 1.244 to 0.840 (-32.7 %); si16 standard T=1 5.301 to 4.254 (-19.8 %), T=8 3.867 to 2.677 (-30.9 %). Density hashes identical on both sides (si8 strict `129b581a1760431b`, si16 strict `00c5a60103b76407`, si8 standard `7e403690df3662c5`, si16 standard `b182bcf8bbc729db`). The expected 5 to 7 % at T=8 was clearly exceeded. Thread scaling T=1 to T=8 (new versus b022b61): si8 strict 1.34x (1.21x), si16 strict 1.81x (1.59x), si8 standard 1.26x, si16 standard 1.59x. The earlier standard T=8 inversion is gone against the previous off timings (si8 on 0.840 s versus off 0.920 s from b022b61; si16 on 2.677 versus off 3.129), but off was not re-measured on 42240e5.
k-point cases (auto symmetry, k444s with fft_grid 32^3), old versus new: k444s T=4 0.736 to 0.706 s (-4.0 %), T=8 0.629 to 0.593 (-5.7 %), T=12 0.497 to 0.503 (+1.2 %, 1/10 faster), T=16 0.506 to 0.485 (-4.1 %); k222 T=4 0.397 to 0.390 (-1.9 %), T=8 0.386 to 0.377 (-2.1 %), T=12 0.377 to 0.379 (+0.4 %), T=16 0.387 to 0.383 (-1.1 %); SCF/Davidson/H unchanged (k444s 10/524/3201, k222 10/149/934). Density hashes constant across T = 4, 8, 12, 16 (k444s `5b3de638e70a7ad0`, k222 `b3b8b21ad8022d32`). k444s saturates near T=12 (10 irreducible points).
D0 profile real on + standard (20 runs, working-sample shares): si8 T=1 / T=8, si16 T=1 / T=8: local FFT 36.8 / 34.7 / 38.1 / 38.5 %, nonlocal 8.3 / 8.7 / 14.8 / 9.8 %, subspace algebra 31.2 / 30.6 / 33.8 / 34.5 %, residual 0.5* / 0.5* / 0.4* / 0.4* % (under 30 samples), other 22.9 / 25.5 / 12.9 / 16.7 %; working samples 590 / 369 / 2261 / 1379; waiting share 3.3 / 82.6 / 0.1 / 79.1 %. With the standard profile local FFT is the largest bucket (strict, earlier: algebra 48.3 %, FFT 31.9 %); nonlocal is larger than under strict (8 to 15 % versus 3 %). si16 T=1 top symbols: `hamiltonian::add` 11.75 %, `FusedOperator::apply_local` 11.36 %, `real_davidson::run_into` 11.14 %, `matmul_1_4_dyn` 7.47 %, `perform_fft_f64` 6.54 %, `fused` 5.30 %, memset 4.68 %.
Pair FFT cost per call (repository test `pair_microbenchmark`, dev08, native; the printed label said "Mac observation" but the numbers are dev08; cubic grids only): 25^3 pair 100.366 us/band versus single 200.144 (ratio 0.501); 36^3 212.856 versus 434.153 (0.490); 48^3 688.607 versus 1257.161 (0.548). MKL full (unpruned) 3D FFT, same grids (separate crate using `dft_mkl::Plan3d`): 36^3 195.01 us per transform (194.26 to 195.55), 72x36x36 488.53 (487.45 to 489.75), 25^3 75.03. QE `fftw` timer per call: si8 36^3 135 us, si16 72x36x36 341 us, si32 795 us. At 36^3: QE 135 us per call, MKL full 195 us per transform, BAANDI pair 212.9 us per band (which includes forward transform, potential multiply and inverse transform for one band of a pair).
Per-vector H time (upper bound, from FFT bucket share x non-waiting share x wall / vectors): si8 T=1 474 us per H vector (949 per pair), si16 T=1 1123 (2246 per pair).
Orchestrator correction in the repository: the pair microbenchmark now includes the non-cubic 72x36x36 grid and the label says `observation`.

# dev08 pair FFT with a non-cubic grid (relayed, 2026-11-07)

Relayed by the dft5090/master session. main `367e884`, `RUSTFLAGS="-C target-cpu=native"`, dev08, repository test `pair_microbenchmark` (label `observation`): 25x25x25 n=1021 pair 100.349812 us/band single 196.471150 ratio 0.510761; 36x36x36 n=3071 213.621162 / 434.794188 / 0.491316; 48x48x48 n=7153 688.117375 / 1263.263975 / 0.544714; 72x36x36 n=6097 486.944363 / 1003.323413 / 0.485331 (the first three reproduce the earlier cubic run within 1 %).
Same-grid comparison (us per transform; one pair call handles two bands with one forward and one inverse transform, so us per band equals us per transform): BAANDI pruned pair path 100.35 / 213.62 / 486.94 versus MKL full unpruned 3D transform 75.03 / 195.01 / 488.53 (25^3 / 36^3 / 72x36x36) and QE `fftw` per call -- / 135 / 341. BAANDI over MKL 1.34x / 1.10x / 1.00x; QE over MKL 0.69x / 0.70x. At 72x36x36 only 6097 of 93312 grid points (6.5 %) hold plane waves, yet the pruned pair path takes the same time as the unpruned MKL transform; QE is 1.43x faster (1.58x at 36^3). Whether one QE fftw call handles one pair (as assumed here: per band 341 us) or two bands each is not established by the relayed data. The pruned pair FFT implementation itself is therefore a candidate (analysis: docs/performance/pair-fft-efficiency.md when written).

# dev08 equal-accuracy gap to QE on main 367e884 (relayed, 2026-11-08)

Relayed by the dft5090/master session. main `367e884`, pure native sha256 `62e3447704ab67583aaa7e2b4a4bbe0b1a5df7991a863a74050efbc1f0ef80dc`, rustc 1.99.1; Si.pz-vbc.UPF `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`; QE `pw.x` `07c7cd7935f42b15157500ab719a5c3948825d8ce8a4487d992d7501c659658f`; `taskset -c 0` / `0-7`; si8 and si16 n=10, si32 n=3. SCF / Davidson / H vectors identical inside each row and between T=1 and T=8. BAANDI tight references unchanged (si8 -31.321803601031771, si16 -63.129244421617670, si32 -126.703998096108023).

Wall median (s) with |dE| against the code's own tight reference and SCF/dav/Hvec:
si8: T=1 off strict 3.846 (1.350e-13; 15/279/2044), on strict 2.613 (3.197e-14; 15/258/2014), off standard 1.340 (1.553e-09; 15/77/772), on standard 1.062 (1.247e-09; 15/78/795); T=8 off strict 2.514, on strict 1.948, off standard 0.925, on standard 0.835.
si16: T=1 off strict 15.344 (1.707e-11; 12/182/2700), on strict 9.993 (7.887e-13; 13/185/2895), off standard 5.932 (1.268e-09; 14/76/1350), on standard 4.257 (5.378e-10; 15/70/1442); T=8 7.573, 5.537, 3.180, 2.648.
si32: T=1 off strict 121.043 (2.942e-11; 13/352/5946), on strict 91.119 (1.026e-11; 14/426/6481), off standard 38.705 (7.067e-10; 15/67/2806), on standard 19.087 (1.945e-09; 13/63/2430); T=8 50.484, 42.667, 16.809, 11.039.
Real (on) versus complex (off) wall ratio, T=1 / T=8: strict si8 -32.1 / -22.5 %, si16 -34.9 / -26.9 %, si32 -24.7 / -15.5 %; standard si8 -20.7 / -9.7 %, si16 -28.2 / -16.7 %, si32 -50.7 / -34.3 %. All T=8 inversions are gone (the earlier standard inversions si8 +30.4 %, si16 +21.5 % were re-checked with off re-measured on the same commit: now -9.7 % and -16.7 %). Thread scaling T=1 to T=8 (off strict / on strict / off standard / on standard): si8 1.53 / 1.34 / 1.45 / 1.27; si16 2.03 / 1.80 / 1.87 / 1.61; si32 2.40 / 2.14 / 2.30 / 1.73.

Equal-accuracy gap to QE real gamma (BAANDI T=1 best setting versus QE np=1; QE faster by; previous value from the b022b61 table):
si8 on standard 1.062 s (1.247e-09) versus QE conv_thr 1e-8 0.467 s (9.230e-10): 2.27x (2.80x); si8 on strict 2.613 s (3.197e-14) versus conv_thr 1e-10 0.551 s (2.613e-11): 4.74x (5.37x); si16 on standard 4.257 (5.378e-10) versus 1.413 (2.699e-10): 3.01x (3.75x); si16 on strict 9.993 (7.887e-13) versus 1.826 (5.436e-12): 5.47x (6.09x); si32 on standard 19.087 (1.945e-09) versus 5.031 (2.588e-10): 3.79x (4.40x); si32 on strict 91.119 (1.026e-11) versus 6.689 (1.820e-11): 13.62x (14.51x).
Trend of the strict-class gap: si8 6.83 to 5.37 to 4.74x, si16 8.28 to 6.09 to 5.47x, si32 17.82 to 14.51 to 13.62x. The growth with system size persists; si32 strict stays furthest and is the only case where the real path does more work (SCF 13 to 14, Davidson 352 to 426). BAANDI delivers better accuracy than QE in the strict class (si8 817 times, si16 6.9 times; si32 comparable). Standard-class BAANDI rows are single unperturbed trajectories (see the P29B study for the median trajectory). Raw data on dev08: `/scr/mino/rustdft/gap/`.

# dev08 raw FFT rates: MKL pure call versus wrapper, 1D batched, rustfft, FFTW3 (relayed, 2026-11-08; corrects the earlier MKL figures)

Relayed by the dft5090/master session. Code outside the repository at `/scr/mino/rustdft/fftbench` (declares DFTI by FFI directly; `RUSTFLAGS="-C target-cpu=native" cargo build --release --features rfft,fftw`; run `taskset -c 0 ./target/release/fftbench`); idle machine (start load 0.29), 32 samples (34 rounds, first 2 dropped), each round averages 20 to 40 inner calls.
**Correction:** the earlier "MKL full 3D 195.0 / 488.5 us" came from `dft_mkl::Plan3d::transform`, which copies the input into an internal buffer and the result back to the output buffer. Pure in-place `DftiComputeBackward` (`DFTI_PLACEMENT=INPLACE`, one buffer, no re-initialization between calls), median / min / max us: 25^3 76.29 / 76.21 / 87.36 (wrapper 74.61: copy share -2.2 %); 36^3 147.25 / 147.15 / 155.22 (wrapper 191.13 / 190.38 / 261.99: copy share 23.0 %); 72x36x36 323.29 / 322.38 / 327.53 (wrapper 486.30 / 484.49 / 490.94: copy share 33.5 %).
MKL 1D complex batched inverse, in place (length, batch: us contiguous / us transposed (stride = batch, distance 1); Gflop/s contiguous): 25 x 1021: 70.86 / 78.89 (8.36); 36 x 1700: 134.30 / 149.25 (11.78); 36 x 1296: 101.33 / 113.32 (11.90); 48 x 1296: 115.99 / 132.07 (14.98); 72 x 1296: 183.28 / 248.33 (15.71); 25 x 625: 43.42 / 46.69; 72 x 2592: 359.38 / 562.77 (16.02). The transposed layout costs +11 % (25), +12 % (36), +14 % (48), +36 % (72; +57 % at batch 2592).
rustfft 6.4 1D complex inverse, same lengths and batches, single buffer (us, Gflop/s, versus MKL contiguous): 25 x 1021 56.46 (10.50, 1.26x faster); 36 x 1700 55.52 (28.50, 2.42x); 36 x 1296 41.62 (28.98, 2.43x); 48 x 1296 96.87 (17.93, 1.20x); 72 x 1296 172.58 (16.68, 1.06x). rustfft resolved offline (already in the registry from the repository build).
FFTW3 3.3.10 (static `~/fq/deps/lib/libfftw3.a`), 3D complex inverse in place, FFTW_MEASURE: 25^3 60.73 us (1.26x faster than MKL pure), 36^3 130.01 (1.13x), 72x36x36 288.08 (1.12x).
Per-call comparison (us): 36^3: QE fftw 135, MKL pure 147.3, FFTW3 130.0, BAANDI pruned pair 213.6; 72x36x36: QE 341, MKL pure 323.3, FFTW3 288.1, BAANDI pruned pair 486.9. QE's fftw timer value is at MKL-pure level (8 percent faster at 36^3, 5 percent slower at 72x36x36), taking one QE call as one band; if a call handles two bands the QE value per band is half. BAANDI's pruned pair path is 1.45x (36^3) and 1.51x (72x36x36) slower than the pure unpruned MKL call, an upper bound because the pair microbenchmark per band also contains half of the potential multiply and the pack/unpack. Replacing the pruned pair transform by pure MKL would cut 31 percent (36^3) and 34 percent (72x36x36) of the per-band figure; FFTW3 39 and 41 percent.
Orchestrator notes: (1) FFTW is GPL-licensed; BAANDI is MIT, so FFTW is a measurement reference only and not a candidate for linking into the project. (2) The `Plan3d` copies also sit in the complex path of MKL builds (GridFusedPlan through the adapter); removing them is bit-identical by construction.

# dev08 P39A (mirror x-slab pair FFT) on main 010b4e7 (relayed, 2026-11-08)

Relayed by the dft5090/master session. main `010b4e7`, pure native sha256 `af1841071dbbdbc7a23be8fe27cb5b24918ad84a0439034d221e006bceb313bf`; comparison `367e884` `62e3447704ab67583aaa7e2b4a4bbe0b1a5df7991a863a74050efbc1f0ef80dc`; round-robin alternation, 10 rounds, `taskset -c 0` / `0-7`.
Pair microbenchmark (same run, library plan versus mirror plan, us per band, label `mirror_observation`, 30 samples): 25^3 106.800 versus 74.370 (mirror/library 0.696; single 215.561); 36^3 224.734 versus 119.015 (0.530; single 464.233); 48^3 675.308 versus 401.946 (0.595; 1290.019); 72x36x36 538.942 versus 270.673 (0.502; 1010.415). The library plan's absolute values rose 5 to 11 percent against the previous run (three plans now share one run: cache interference), so only the ratios within one run are reliable; comparisons with the earlier MKL (pure 323.29 / 147.25 / 76.29 us) and QE (341 / 135 us) figures cross runs: mirror 270.67 versus MKL 323.29 at 72x36x36 (1.19x faster), 119.02 versus 147.25 at 36^3 (1.24x), 74.37 versus 76.29 at 25^3 (1.03x), versus QE 1.26x (72x36x36) and 1.13x (36^3) faster. The earlier statement that pruned beats nothing holds only for the previous plan.
SCF rows (real path on), 367e884 versus 010b4e7: counts and energies differ (rounding change, density hashes differ in all rows; expected); ms per H vector falls in all rows (-2.1 to -19.2 %). si8 strict T=1: SCF 15 / 15, Davidson 258 / 292, H 2014 / 2137, |dE| vs tight 3.197e-14 / 3.695e-13, ms/H 1.3004 / 1.0852 (-16.6 %); T=8 0.9687 / 0.8521 (-12.0 %); si8 standard T=1: SCF 15 / 19, Davidson 78 / 100, H 795 / 950, |dE| 1.247e-09 / 3.005e-09, ms/H 1.3358 / 1.1600 (-13.2 %); T=8 1.0566 / 0.9863 (-6.7 %); si16 strict T=1: 13 / 13, 185 / 242, 2895 / 3053, 7.887e-13 / 7.390e-13, 3.4518 / 3.1602 (-8.4 %); T=8 1.9136 / 1.8742 (-2.1 %); si16 standard T=1: 15 / 12, 70 / 65, 1442 / 1258, 5.378e-10 / 6.707e-10, 2.9730 / 2.4030 (-19.2 %); T=8 1.8391 / 1.6137 (-12.3 %).
Wall median (s) old to new (candidate faster in x of 10 rounds): si8 strict T=1 2.619 to 2.319 (-11.45 %, 10/10), T=8 1.951 to 1.821 (-6.7 %, 10/10); si8 standard T=1 1.062 to 1.102 (+3.75 %, 1/10), T=8 0.840 to 0.937 (+11.5 %, 0/10); si16 strict T=1 9.993 to 9.648 (-3.5 %, 10/10), T=8 5.540 to 5.722 (+3.3 %, 0/10); si16 standard T=1 4.287 to 3.023 (-29.5 %, 10/10), T=8 2.652 to 2.030 (-23.5 %, 10/10). Wall differences are dominated by trajectory changes (si8 standard SCF 15 to 19, si16 standard 15 to 12); the single-trajectory rows are not a robust comparison: the SCF-count distribution under perturbation must decide.

# dev08 real-path SCF distribution, P39A before/after (relayed, 2026-11-08)

Relayed by the dft5090/master session. Harness file of main `316aa50` (`P29A_GAMMA_REAL=on`) copied onto each checkout (`367e884` before P39A, `010b4e7` after), built and run in turn (`RUSTFLAGS="-C target-cpu=native"`, same command as before with `P29A_CASES=si8-gamma,si16-gamma`), load below 0.30 before each; the real path really ran (ROW/END/GAP diff 1137 lines, whereas the complex-path run gave 0 lines). Raw logs `/scr/mino/rustdft/p29r/real_367e884.log`, `real_010b4e7.log`.
si8-gamma standard SCF over ten 1e-10 perturbations: 367e884 [20, 22, 22, 22, 25, 26, 28, 29, 30, 31] (min 20, median 25.5, max 31); 010b4e7 [20, 20, 20, 25, 25, 26, 27, 29, 31, 32] (min 20, median 25.5, max 32). H-work median 1206.5 to 1190.5 (-1.3 %); strict run 0 SCF 18/18, H 2179 to 2236 (+2.6 %), run 1 SCF 20/20, 2543 to 2600 (+2.2 %); standard |dE| against strict run 1: min 3.566e-10 to 2.660e-10, median 1.420e-09 to 8.786e-10, max 3.963e-09 to 2.727e-09.
si16-gamma: 367e884 [15, 15, 15, 15, 17, 19, 20, 22, 24, 25] (min 15, median 18.0, max 25); 010b4e7 [15, 15, 15, 18, 18, 20, 20, 22, 25, 25] (median 19.0); H-work median 1614.0 to 1687.0 (+4.5 %); strict run 0 SCF 17 to 18, H 3353 to 3540 (+5.6 %), run 1 20/20, 3898 to 3913 (+0.4 %); |dE| min 2.786e-10 to 3.588e-10, median 7.544e-10 to 1.353e-09, max 1.794e-09 to 2.380e-09.
Reading: si8 shows no shift of the distribution (min and median equal); the unperturbed 15 to 19 iterations of the bench case is one trajectory. si16 shows a small worsening (median +1 iteration, H work +4.5 %, error median 1.8x), smaller than the per-H-vector gain of P39A (-8.4 to -19.2 %). Per-run scatter (same run index, si8 run 2: 25 to 20, run 3: 22 to 25; si16 run 11: 15 to 20) is wider than the single-trajectory difference seen earlier. Real versus complex path SCF medians (367e884 real versus complex distributions): si8 25.5 versus 24.5, si16 18.0 versus 20.0.

# dev08 P39B (MKL copy removal) tests and microbenchmark on main 5ee9b83 (relayed, 2026-11-08)

Relayed by the dft5090/master session. main `5ee9b83`, MKL 2026.1 (`MKL_NUM_THREADS=OMP_NUM_THREADS=1`), `RUSTFLAGS="-C target-cpu=native"`. Tests: `dft-mkl in_place_matches_wrapper_bits` ok (1 passed); `dft-numerics --test in_place` 4 passed (default in-place preserves bits and validates dimensions; rust grid preserves legacy local and density bits; mkl in-place preserves bits and validates dimensions; mkl grid preserves legacy pruned-pattern bits); `dft-numerics --features mkl` whole crate 39 passed, 0 failed, 6 ignored (all other test targets 0 failed).
In-place microbenchmark (`taskset -c 0`, 32 samples, load 0.14), wrapper median (min to max) versus in place, us: 25^3 inverse 74.596 (74.174 to 83.816) versus 60.266 (60.151 to 62.569) -19.2 %, forward 73.737 versus 59.437 -19.4 %; 36^3 inverse 204.494 (189.075 to 208.802) versus 133.037 (129.324 to 138.323) -34.9 %, forward 208.384 versus 137.218 -34.1 %; 72x36x36 inverse 487.529 (482.020 to 567.810) versus 292.588 (287.595 to 381.991) -40.0 %, forward 497.343 versus 302.477 -39.2 %. This agrees in direction and size with the independent raw measurement (wrapper 486.3 versus pure call 323.3 us at 72x36x36); the repository in-place path is 10 % faster than that pure call, probably because the independent run set the backward scale per call and initialized the buffer differently.
Bit identity of MKL SCF results: baseline MKL binary built from the parent commit `316aa50` (sha256 `001c37f3249b5b6e1a10c589d08ffc79cfe4f109a170c43f26c40e2ef8c3aa1f`) versus candidate MKL binary `efbaf8066d226563e0f2b49594ef0ad05516112049fa93fa543e9711f228ce21`: `bench_scf.py --cases gamma,k222s,k222,k444s,si8 --threads {1,2,4,8,12} --repeats 1`, energy hex, density sha256, SCF, operator applications and Davidson steps all equal: 25/25 bit-identical with unchanged work. SCF timing (base MKL versus candidate MKL versus pure native, alternating, 10 rounds) is running. The dev08 home partition has 6.1 GB left.
