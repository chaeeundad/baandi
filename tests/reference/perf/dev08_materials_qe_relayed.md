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
