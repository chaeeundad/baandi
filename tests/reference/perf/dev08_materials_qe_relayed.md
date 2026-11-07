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
