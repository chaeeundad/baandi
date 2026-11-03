# Validation record

Only measured results are recorded here, each with its command, inputs and checksums. Nothing is copied from memory or from a reference that was not actually run. QE reference numbers will be added only from real pw.x runs returned with their provenance (see docs/reference-environment.md).

## Pseudopotential provenance (applies to every run below)
`Si.pz-vbc.UPF`, sha256 `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`, 74552 bytes. Author, date, generator input and license are **blank or unknown** in the file; provenance is **unverified**. Use approved by the user on 2026-11-02 for internal validation only, no redistribution. The file is not in the repository. Facts: docs/candidate-si-upf.md.

## R1: Si diamond primitive cell, Gamma only, LDA-PZ, 30 Ry (first real-file run)
- Engine commit: `4fc299d3215cbfc7248105d0be5a886591f69009` (release build, Rust 1.99.1, single thread, faer sequential).
- Command: `./target/release/dft-cli run ../rust-dft-wt/runs/si_30ry/input.toml`
- Input sha256 `30d9dc99f2439f3997d9583bf67762d32b1397e6137300c0edb3a3abc2dcf2bd`; a = 5.43 Angstrom, atoms at (0,0,0) and (1/4,1/4,1/4); ecutwfc 30 Ry; Pulay mixing beta 0.5; tolerances: density 1e-7, energy 1e-8 Ha, eigen 1e-8 Ha.
- Basis: 749 plane waves, 5985 density G, FFT 25x25x25 (our rule, not necessarily QE's); 8 electrons, 8 bands (4 occupied, 4 empty).
- Result: converged in 11 iterations, wall time 1.2 s. Total energy -7.295530468210927 Ha = -14.591060936421854 Ry.
- Terms (Ha): kinetic 4.143257239948187; local -2.950156276146177 (including the G=0 alpha-Z part -0.145624326070663); nonlocal 1.598682517429664; Hartree 0.836046055268684; XC -2.523888138201652; Ewald -8.399471866509632. Harris-Foulkes functional -7.295530468210838 Ha (difference 9e-14 Ha); double-counting integral -1.620206350107395 Ha.
- Eigenvalues at Gamma (eV): -5.1733, 7.0378 (x3, degenerate to 1e-12 Ha), 9.1742 (x3), 10.2981; HOMO/LUMO gap 2.1364 eV.
- Observations (not claims of agreement): the triple degeneracies are preserved to 1e-12, valence bandwidth at Gamma (Gamma_1 to Gamma_25') is 12.21 eV, Gamma_2' to Gamma_25' is 3.26 eV. Gamma-only sampling of a two-atom cell is a coarse Brillouin-zone sample; the absolute energy is not comparable to k-converged literature values.
- Status when recorded: not yet compared with QE. Superseded: see R3 to R5 below for the QE comparisons.

## Real-world UPF deviations found by this run (fixed)
1. `relativistic="no"` (legacy value) instead of the specification's "nonrelativistic" (P2C).
2. Numbers separated by newlines only; the `upf` crate splits on spaces (P2C).
3. Functional label `SLA  PZ   NOGX NOGC` instead of exactly `SLA PZ` (P5C).
These were invisible to synthetic fixtures; the lesson is recorded in docs/progress.md.

## QE comparison protocol

P6 tooling produces no QE reference values. Generate inputs with
`python3 -I scripts/make_qe_input.py`, and validate returned bundles with
`python3 -I scripts/verify_qe.py`. The authoritative
[bundle format and comparison rules](../tests/reference/README.md) specify exact
CLI options, parser patterns, required provenance, and report semantics. The
remote execution protocol remains [reference-environment.md](reference-environment.md).

A bundle contains the executed `qe.in`, verbatim `qe.out`, and
`provenance.json` with command, printed version, host, MPI/OMP counts, UTC date,
input/UPF SHA-256, measured wall time, runner, notes, and fixture marker. Input
hash, UPF hash, and printed version must match before any comparison.
FORMAT FIXTURES always refuse comparison. No real reference is fabricated or
inferred from an internally consistent synthetic test.

Conditions must match: Bohr cell and positions (1e-12 Bohr), cutoffs in Ry,
electrons, bands, plane waves, fixed occupations, full complex Gamma sampling,
and reported LDA-PZ functional. Explicit FFT dimensions must match both our
results and QE output. Automatic QE grids are recorded beside ours with a
warning and require convergence-study interpretation. QE's density cutoff is
four times its wavefunction cutoff. Geometry uses 17 significant digits with
our matrix columns transposed to QE vector rows. The UPF's own functional is
used; `input_dft` is not set.

The total-energy tolerance is 1e-4 Ha/atom; eigenvalues and gap use 1e-3 Ha;
electron count uses 1e-8 relative. Decomposition terms use 1e-4 Ha and map
one-electron to kinetic + local + nonlocal, plus Hartree, XC, and Ewald directly.
An explicit documented convention difference may make the decomposition
informational, but cannot exempt total energy. The project also requires
orthogonality 1e-8 and eigen residual 1e-7 Ha. The verifier checks the exported
occupied eigen residual; it explicitly marks orthogonality unassessed because
these files lack an overlap matrix. It does not claim density comparison.
All used tolerances and nondefault overrides appear in reports. Tolerances
are never loosened to pass; trace units, FFT normalization, G=0, projectors,
double counting, and XC parameterization before changing physics.

Exit codes: 0 comparison passes; 1 usage/our-results/I/O error; 2 provenance
incomplete or mismatched; 3 unsupported/missing parser record or conditions
mismatch; 4 format fixture; 5 completed comparison exceeds a tolerance.
JSON and Markdown are written after all prerequisite checks, with provenance,
units, absolute differences, maxima and attaining items. For the candidate
Si digest, every report states "Pseudopotential provenance is unverified".
Only real returned, provenance-complete runs may later populate reference data.

## Convergence studies

```sh
python3 -I scripts/convergence_scan.py --input examples/si_lda.toml \
  --dft-cli target/release/dft-cli --out-dir runs/scan --scan all \
  --emit-qe-inputs --pseudo-dir /remote/pseudo --upf-name Si.pz-vbc.UPF
```

Supply a valid input with an accessible pseudopotential path; the example path
is a placeholder. The driver does not alter the source. Each point receives a
derived TOML, isolated output directory, stdout/stderr files, and a timed
`dft-cli run INPUT --output-dir DIR` command. Existing scan directories are
refused to prevent reuse of stale results. A nonzero exit, timeout (default
14400 seconds, `--timeout` override), or uncertified result aborts the scan and
records a failed point in `manifest.json`. Partial runs do not produce a success
table. Python 3.9 is supported without a TOML dependency: changed tables require
documented bare keys and single-line values; pseudopotential paths require
single-line quoted strings. Unrelated text is preserved, relative pseudo paths
are made absolute, and the real CLI validates derived TOML.

The independent series are:

- `ecut`: 20,25,30,40,50,60,80 Ry by default (`--ecuts` overrides). Remove any
  grid override so each cutoff gets its own automatic alias-free grid. Other
  settings remain those of the source input.
- `grid`: at the source cutoff, first run the automatic grid. Enumerate the
  density sphere using the first result's reciprocal matrix and cutoff (density
  radius squared = 4*ecutwfc_ry, relative boundary allowance 1e-12, matching the
  engine). Its maximum absolute Miller indices give the just-sufficient odd
  dimensions `2*max_abs_m+1`. Also run two larger grids, rounding each component
  of 1.25 and 1.5 times the automatic dimensions upward to a 2/3/5-smooth integer.
  Explicit minimum may equal automatic; retain both to exercise the override.
  `--grids 24x24x24,30x30x30` substitutes explicit candidates after automatic.
- `tolerance`: density residual 1e-5,1e-6,1e-7,1e-8,1e-9 Bohr^(-3/2).
  Set both energy tolerance in Ha and eigen tolerance in Ha to 0.1 times each
  density value. Cutoff and optional grid remain the source settings. These
  solver thresholds are recorded scan parameters, not relaxed acceptance gates.

`scan.json`, `scan.csv`, and `scan.md` record point ID, parameter, total Ha,
Ha/atom, signed total-energy difference to the finest point of that series,
HOMO-LUMO gap eV (null without an empty band), iterations, and CLI wall seconds.
Finest means largest cutoff, largest grid volume, or smallest density tolerance;
JSON/CSV retain its point ID and ordering rank. For arbitrary nonnested grids,
volume is only an ordering convention, not a proof of quadrature improvement.

`manifest.json` records source input SHA-256, invocation, each exact subprocess
argument list, derived input SHA-256, result SHA-256, exit status, elapsed wall
seconds, and the binary's recorded `tool.git_commit` (including `unknown` if
that is what the binary reports). A commit identifies the build's HEAD and does
not certify a clean worktree. With `--emit-qe-inputs`, every point emits both
explicit-our-grid and automatic-QE-grid inputs plus their sidecars. Manifest
entries bind each input hash to its matching results path/hash and expected UPF
digest. No QE execution or comparison occurs during a scan.

## R2: rust-dft convergence scans, real Si (ours only, no QE yet)
- Engine commit: `bf65e26f15a060f3c180f7fbbaeaf852a34d2368`; release build; single thread.
- Command: `python3 -I scripts/convergence_scan.py --input ../rust-dft-wt/runs/si_30ry/input.toml --dft-cli ./target/release/dft-cli --out-dir ../rust-dft-wt/runs/si_scan --scan all --emit-qe-inputs --pseudo-dir /home/mino/fq/pseudo --upf-name Si.pz-vbc.UPF` (32 s).
- Base input sha256 `30d9dc99f2439f3997d9583bf67762d32b1397e6137300c0edb3a3abc2dcf2bd`; UPF sha256 `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217` (provenance unverified, see top of this file).
- Results (energies in Ha, two atoms; the "difference to finest" column is relative to the last/finest point of each scan):

| point_id | parameter | total_energy_ha | energy_per_atom_ha | difference_to_finest_ha | gap_ev | iterations | wall_time_s |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 000-ecut-20Ry | 20Ry | -7.291089364151892 | -3.645544682075946 | 0.004951545907649901 | 2.132876256889901 | 11 | 0.475500333 |
| 001-ecut-25Ry | 25Ry | -7.294268118369173 | -3.6471340591845864 | 0.0017727916903691465 | 2.136868446135041 | 11 | 0.510810167 |
| 002-ecut-30Ry | 30Ry | -7.295530468210927 | -3.6477652341054636 | 0.0005104418486148177 | 2.136426740250426 | 11 | 1.172997791 |
| 003-ecut-40Ry | 40Ry | -7.295929509540116 | -3.647964754770058 | 0.00011140051942604146 | 2.136181653414427 | 11 | 1.6263766670000002 |
| 004-ecut-50Ry | 50Ry | -7.296019054644814 | -3.648009527322407 | 2.1855414727767197e-05 | 2.13599958636852 | 11 | 3.229204833 |
| 005-ecut-60Ry | 60Ry | -7.296038607887422 | -3.648019303943711 | 2.3021721196059275e-06 | 2.1359951394486782 | 11 | 3.177454084 |
| 006-ecut-80Ry | 80Ry | -7.296040910059542 | -3.648020455029771 | 0.0 | 2.135991484381121 | 11 | 6.002052917 |
| 007-grid-automatic | automatic | -7.295530468210927 | -3.6477652341054636 | 4.132113762267409e-08 | 2.136426740250426 | 11 | 1.189520917 |
| 008-grid-0-25x25x25 | 0-25x25x25 | -7.295530468210927 | -3.6477652341054636 | 4.132113762267409e-08 | 2.136426740250426 | 11 | 1.176850291 |
| 009-grid-1-32x32x32 | 1-32x32x32 | -7.295530546544072 | -3.647765273272036 | -3.7012006792735974e-08 | 2.1364271673422257 | 11 | 1.6912443750000001 |
| 010-grid-2-40x40x40 | 2-40x40x40 | -7.295530509532065 | -3.6477652547660324 | 0.0 | 2.136427387739426 | 11 | 3.505630833 |
| 011-tolerance-1e-05 | 1e-05 | -7.295530468209554 | -3.647765234104777 | 1.3873346915715956e-12 | 2.136426261911646 | 8 | 1.064756625 |
| 012-tolerance-1e-06 | 1e-06 | -7.295530468209568 | -3.647765234104784 | 1.3731238368563936e-12 | 2.136426261911551 | 8 | 1.100988959 |
| 013-tolerance-1e-07 | 1e-07 | -7.295530468210927 | -3.6477652341054636 | 1.4210854715202004e-14 | 2.136426740250426 | 11 | 1.186914625 |
| 014-tolerance-1e-08 | 1e-08 | -7.295530468210915 | -3.6477652341054574 | 2.6645352591003757e-14 | 2.13642680073045 | 15 | 1.245200917 |
| 015-tolerance-1e-09 | 1e-09 | -7.295530468210941 | -3.6477652341054707 | 0.0 | 2.1364268001950073 | 18 | 1.333065417 |

- Reading: (1) cutoff: the energy decreases monotonically; at 30 Ry it is 5.1e-4 Ha (2.6e-4 Ha/atom) above the 80 Ry value, at 50 Ry 2.2e-5 Ha, at 60 Ry 2.3e-6 Ha. A comparison with QE must use the same cutoff, not an absolute converged value. (2) FFT grid at 30 Ry: 25^3, 32^3, 40^3 differ by at most 7.8e-8 Ha (XC quadrature); the gap changes by 1e-6 eV. (3) SCF tolerance: density tolerance 1e-5 already reproduces the 1e-9 energy to 1.4e-12 Ha (the energy is stationary, error quadratic in the density residual); eigenvalue/gap digits keep changing at the 1e-7 eV level until the tolerance is 1e-7.
- Status: **ours only.** Matching QE inputs were generated for every scan point (22 files, manifest in the scan directory) but no QE run has been compared yet.

## R3: first QE comparison, Si 30 Ry Gamma-only (real pw.x output)
- Reference: QE PWSCF 7.6 (build orig-mpi-dfti, OpenMPI 5.0.8, MKL 2026.1) on dev08, 4 MPI processes, 1 thread each, run by the dft5090 session. Bundles: `tests/reference/si30_ours` (FFT 25x25x25 requested) and `tests/reference/si30_auto` (QE automatic grid; it chose 25x25x25 as well). Inputs generated by `scripts/make_qe_input.py` from our R1 results; input sha256 `017e91ea6dc1c1d0e382f3a8158a52248f155379ae360ceb0ef6f3d3f2554065` (ours) and `61f0499c2056a5668c1b3a5ad6004aa6c0e8f44abb3cbe93185b34b51d209887` (auto); stdout sha256 `9ee1a3a2feb9f5ddc26e5271ca1084a3b732fde363cbeaeb2ebefe308add2d41` and `d677d4e6724f801bf9d95d8fe4b38c985261062b4dfec698b5bc10c57944725e`; UPF sha256 `e8d93375...5217` (provenance unverified). Full provenance: the bundles' provenance.json.
- Our side: R1 results (engine commit recorded there; comparison tool at commit `9e26346022c2e2de56ea9b907eb364c92ccd2f11`).
- Command: `python3 -I scripts/verify_qe.py --ours tests/reference/si30_ours/ours/results.json --bundle tests/reference/si30_ours --json tests/reference/reports/si30_ours.json --markdown tests/reference/reports/si30_ours.md` (same for si30_auto). Exit code 0 for both: all tolerances passed.
- Conditions that matched (checked by the tool before comparing): cutoffs 30 Ry and 120 Ry, 8 electrons, 8 bands, 749 plane waves, 5985 dense G-vectors, FFT 25x25x25 (also QE's own automatic choice), lattice and positions within 1e-12 Bohr, complex Gamma sampling, fixed occupations, functional string "SLA  PZ   NOGX NOGC", 11 SCF iterations on both sides.
- Result table (Ha unless stated; tolerances are the project spec defaults and were not changed):

| Item | Ours | QE | Absolute difference | Unit | Tolerance | Result |
| --- | ---: | ---: | ---: | --- | ---: | --- |
| total energy per atom | -3.6477652341054636 | -3.6477665149999998 | 1.28089454e-06 | Ha/atom | 0.0001 | PASS |
| one_electron | 2.7917834812316733 | 2.7917809249999999 | 2.55623167e-06 | Ha | 0.0001 | PASS |
| hartree | 0.83604605526868414 | 0.836046065 | 9.73131586e-09 | Ha | 0.0001 | PASS |
| xc | -2.5238881382016523 | -2.523888135 | 3.20165228e-09 | Ha | 0.0001 | PASS |
| ewald | -8.3994718665096322 | -8.3994718850000005 | 1.84903683e-08 | Ha | 0.0001 | PASS |
| eigenvalue band 1 | -0.19011399644747087 | -0.19011526841131599 | 1.27196385e-06 | Ha | 0.001 | PASS |
| eigenvalue band 2 | 0.25863418733601878 | 0.2586343796078247 | 1.92271806e-07 | Ha | 0.001 | PASS |
| eigenvalue band 3 | 0.25863418733659288 | 0.2586343796078247 | 1.92271232e-07 | Ha | 0.001 | PASS |
| eigenvalue band 4 | 0.2586341873370463 | 0.2586343796078247 | 1.92270778e-07 | Ha | 0.001 | PASS |
| eigenvalue band 5 | 0.33714642191919358 | 0.33714563150389404 | 7.904153e-07 | Ha | 0.001 | PASS |
| eigenvalue band 6 | 0.3371464219196072 | 0.33714563150389404 | 7.90415713e-07 | Ha | 0.001 | PASS |
| eigenvalue band 7 | 0.33714642192003197 | 0.33714563150389404 | 7.90416138e-07 | Ha | 0.001 | PASS |
| eigenvalue band 8 | 0.37844791825995638 | 0.37844819469711266 | 2.76437156e-07 | Ha | 0.001 | PASS |
| HOMO-LUMO gap | 0.078512234582147289 | 0.078511251896069334 | 9.82686078e-07 | Ha | 0.001 | PASS |
| electron count | 8 | 8 | 0 | electrons | 8e-08 | PASS |
| our integrated relative charge error | 0 | 6.6613381477509392e-16 | 6.66133815e-16 | relative | 1e-08 | PASS |
| our occupied eigen residual | 0 | 4.3152503285188984e-09 | 4.31525033e-09 | Ha | 1e-07 | PASS |

- Reading: the total energy differs by 2.56e-6 Ha for the cell (1.28e-6 Ha/atom, tolerance 1e-4), 78 times below the tolerance. Hartree, XC and Ewald agree to 1e-8 Ha or better; the whole difference sits in the one-electron term (kinetic + local + nonlocal). The cause of the residual 2.6e-6 Ha is **not yet identified** (candidates, untested: how the radial integrals and tabulation of the pseudopotential are done in each code). Eigenvalues agree to 1.3e-6 Ha or better (QE prints them in eV with four decimals, which limits their resolution to about 1.8e-6 Ha). Orthogonality is not exported by either side and was not assessed in this comparison (it is tested inside our engine).
- Scope: one cell, one cutoff, Gamma only. This is not yet the cutoff-convergence comparison (R4).

## R4: QE cutoff scan and explanation of the constant offset
- QE runs: PWSCF 7.6 on dev08 through the dft5090 session, same build and settings as R3. Six explicit-grid runs (20, 25, 40, 50, 60, 80 Ry; FFT 24, 24, 30, 36, 36, 45 cubed, equal to ours) were returned in full and stored as bundles `tests/reference/scan_e20 ... scan_e80` (stdout verified by piece, stream, gzip, tar and per-file SHA-256; inputs by SHA-256 on the runner). Reports: `tests/reference/reports/scan_e*.md`, all exit code 0 from `scripts/verify_qe.py`. The same six points with QE's automatic FFT grid were also run; the runner returned only the grid, plane-wave, energy and iteration lines (full outputs kept on dev08, informational): **QE's automatic grid equals our grid at every cutoff** (24, 24, 30, 36, 36, 45), energies identical to the explicit-grid runs, 11 iterations everywhere.
- Our side: scan R2 (tooling commit `0228b282892c27b24418cb98ac65a5eb69281f00`).
- Total energies (QE values are the printed 8-decimal Ry values divided by 2):

| cutoff (Ry) | ours (Ha) | QE (Ha) | ours - QE (Ha) | ours rel. to 80 Ry (Ha) | QE rel. to 80 Ry (Ha) |
| --- | ---: | ---: | ---: | ---: | ---: |
| 20 | -7.291089364 | -7.291091920 | +2.556e-06 | +4.951546e-03 | +4.951535e-03 |
| 25 | -7.294268118 | -7.294270675 | +2.557e-06 | +1.772792e-03 | +1.772780e-03 |
| 30 | -7.295530468 | -7.295533030 | +2.562e-06 | +5.104418e-04 | +5.104250e-04 |
| 40 | -7.295929510 | -7.295932070 | +2.560e-06 | +1.114005e-04 | +1.113850e-04 |
| 50 | -7.296019055 | -7.296021620 | +2.565e-06 | +2.185541e-05 | +2.183500e-05 |
| 60 | -7.296038608 | -7.296041155 | +2.547e-06 | +2.302172e-06 | +2.300000e-06 |
| 80 | -7.296040910 | -7.296043455 | +2.545e-06 | +0.000000e+00 | +0.000000e+00 |

- Result: the difference ours - QE is **constant in cutoff**: mean 2.556e-6 Ha, range 2.545e-6 to 2.565e-6 Ha (QE prints 8 decimals in Ry, rounding +-2.5e-9 Ha) across 20 to 80 Ry, i.e. 1.28e-6 Ha/atom, 78 times below the 1e-4 Ha/atom target. The convergence curves of the two codes are parallel to better than 2e-9 Ha relative to the 80 Ry value (error against 80 Ry: 4.95e-3 Ha at 20 Ry, 5.10e-4 Ha at 30 Ry, 2.2e-5 Ha at 50 Ry in both codes).
- Explanation of the constant offset (measured on the real pseudopotential with a diagnostic program that is not part of the repository): the local potential's G=0 constant (alpha-Z) can be computed with an erf-subtracted integrand (our method, exact analytic tail pi*Z/Omega) or directly as 4*pi/Omega * integral r*(r*V + Z) dr on the mesh. On this file they give -9.101520379e-3 Ha and -9.101679477e-3 Ha per atom, which differ by -1.591e-7 Ha. The direct integrand behaves like Z*r near the origin and the radial mesh starts at r0 = 1.30826e-3 Bohr, so a mesh integral starting at the first point omits Z*r0^2/2 of the integral; times 4*pi/Omega this is 1.5925e-7 Ha per atom. With 8 electrons and 2 atoms the total shift is 8*2*1.5925e-7 = **2.548e-6 Ha**, which reproduces the observed 2.545e-6 to 2.565e-6 Ha offset to about 0.3 percent. Replacing our alpha by the direct-integral value moves our R1 energy from -7.295530468 Ha to -7.295533016 Ha, against QE's printed -7.29553303 Ha (difference 1.4e-8 Ha; QE's rounding alone is 2.5e-9 Ha). This explains the offset quantitatively, but **QE's internal procedure was not read**, so it is stated as consistent with, not proven for, QE's implementation. Our erf-subtracted integrand goes like r^2 at the origin, so the omitted segment is negligible and our value is the more accurate one.
- Status: the P6 comparison gates for energy, eigenvalues, electron count, iteration count and cutoff convergence are met on one cell. The density comparison (grid by grid) has not been done.

## R5: density comparison in Fourier space (QE charge-density.dat, no pp.x)
- QE source: the native `charge-density.dat` of run `qe_si30_ours` (QE 7.6, dev08, sha256 `7fdc2ce8d509ea587cc4fc55e0cc41666cc7b64d3e011902f03e51de14000f27`, 167696 bytes). pp.x does not exist in any of the runner's builds; the file's four Fortran records (markers 12, 72, 71820, 95760, no trailing bytes) hold gamma_only=0, ngm_g=5985, nspin=1, the Miller indices and rho(G) in e/Bohr^3. The runner parsed it read-only with a standard-library script and returned the summary stored in `tests/reference/si30_ours/density/qe_rhoG_summary.json` (rho(0), sum of |rho|^2, maximum coefficient, 12 sample coefficients, symmetry check: no missing -G, no conjugation violation beyond 1e-12, largest deviation 1.2e-18).
- Ours: `tests/reference/si30_ours/ours/density.bin` (R1), coefficients computed independently in pure Python by a separable DFT (`scripts/qe_rhog.py`).
- Command: `python3 -I scripts/compare_rhog.py --ours-dir tests/reference/si30_ours/ours --qe-summary tests/reference/si30_ours/density/qe_rhoG_summary.json --json-report tests/reference/reports/si30_density.json --markdown-report tests/reference/reports/si30_density.md` -> exit 0.
- Result (defaults, none changed): rho(G=0) difference 0 (both 0.029617874491946364 e/Bohr^3, equal to N/Omega); relative difference of sum |rho|^2 1.1e-8 (tolerance 1e-6); largest coefficient difference over the 12 samples 4.0e-9 e/Bohr^3 at G=(0,1,0) (tolerance 1e-7). Symmetry spread of the three cubic-equivalent coefficients (1,0,0), (0,1,0), (0,0,1): ours 2.1e-11, QE 6.8e-9 e/Bohr^3, i.e. the observed difference is at the size of QE's own symmetry-breaking noise from `nosym` and its SCF convergence.
- Limits: summary mode compares the sampled coefficients and aggregates, not all 5985 coefficients; full-set equality is not certified (the report says so). Full-file mode (`--qe-file`) is implemented and tested on synthetic files but has not been run on the real file because the 167 KB file was not transferred. One cell, one cutoff.

## R6: k-point comparison with QE 7.6 (Si, 30 Ry, three Monkhorst-Pack grids)
- Code: P7A (weighted Bloch SCF, ADR 0011) and P7A2 (CLI k-point input). Inputs: R1 input with `kgrid`/`kshift` instead of `kpoints = [[0,0,0]]`; QE inputs generated by `scripts/make_qe_input.py` (`K_POINTS automatic`, `nosym = .true.`; QE therefore reduces only by time reversal, as we do), FFT grid 25^3 as ours, `conv_thr` 1e-12, same UPF (sha256 `e8d93375...5d68d5217`, provenance unverified, see `docs/mvp-review.md`). The code under test is commit `dc28615`-era main (post P7A/P7A2); the real-Si Gamma run on the same code gives a bit-identical density file to R1.
- QE runs: pw.x 7.6 `~/fq/builds/orig-mpi-dfti/bin/pw.x`, MPI x4, executed on dev08 by the dft5090 session (relayed; the input files were checked by sha256 before running: `40609d7a...`, `5a376608...`, `a45ed650...`). Evidence: `tests/reference/kpoints/<run>/` (qe.in, qe.in.meta.json, ours_results.json; qe.out for the two small runs, transferred whole and verified by sha256 of the gz and of the decoded file: `f93516f6...afaf2` for k222_c, `162a8980...e45ea6` for k222s_c). For k444s_c the 27 kB output was not transferred whole; `qe_summary.json` holds the relayed plain-text k list and eigenvalue blocks (the full file stays on dev08, sha256 `c93d458b...b035`, 655 lines), so that run is a relayed attestation, not a re-verifiable file.
- Command: `python3 -I scripts/compare_kpoints_qe.py --ours tests/reference/kpoints/<run>/ours_results.json --qe-out tests/reference/kpoints/<run>/qe.out` (`--qe-summary .../qe_summary.json` for k444s_c). Defaults, none changed: total energy 1e-4 Ha/atom, eigenvalues 1e-3 Ha, normalized k weights 1e-6. QE prints weights summing to 2 (spin), ours sum to 1; both are normalized before comparing.

| run | grid and shift | k points ours / QE | iterations QE | total energy ours - QE (Ha/atom) | max eigenvalue diff (Ha) | k-weight diff |
|---|---|---|---|---|---|---|
| k222_c | 2x2x2, 0 | 8 / 8 | 9 | +1.2814e-6 | 1.99e-6 | 0 |
| k222s_c | 2x2x2, half step | 4 / 4 | 9 | +1.2813e-6 | 1.77e-6 | 0 |
| k444s_c | 4x4x4, half step | 32 / 32 | 8 | +1.2817e-6 | 2.10e-6 | 0 |

- All three exit 0 (tolerances met by factors of about 80 for the energy and 500 for the eigenvalues). The total-energy offset equals the Gamma-only offset of R3/R4 (1.28e-6 Ha/atom) to 3 digits, so the k-point implementation adds no measurable error beyond the already explained G=0 integration convention. The k point sets agree one by one up to k -> -k (matched on crystal coordinates), and every compared eigenvalue (all bands printed by QE, 8 bands) agrees within the 4-decimal eV printout (1e-4 eV = 3.7e-6 Ha).
- Energy components (QE, Ry): Ewald -16.79894377 in all three (same structure). Our components per run are in `ours_results.json`.
- Sampled gaps (ours; QE's HOMO/LUMO lines agree): k444s_c HOMO 5.7551 eV and LUMO 6.8875 eV in QE, ours 5.75509 and 6.88746 eV; the 2x2x2 grids give different indirect gaps (0.468 eV unshifted, 2.176 eV shifted) because the gap of a sampled set depends on which k points are included. These are sampling-dependent diagnostics, not a converged band gap.
- Limits: one cell, one pseudopotential, 30 Ry only; density in real space or Fourier space was not compared for the k-point runs; orbitals were not compared (only eigenvalues and total energy); the QE outputs are relayed attestations (two with whole-file hash checks, one without); no symmetry reduction is exercised; no performance comparison was made.

## R7: performance baseline against QE 7.6 and cross-platform reproducibility (dev08)
- Evidence: `tests/reference/perf/dev08_baseline.json` (relayed attestation from the dft5090 session; outputs stay on dev08 under `~/rustdft_p6/perf`, they were not transferred, so this record cannot be re-verified from the repository). Code: github commit `3274955` built with `cargo build --release -p dft-cli`, rustc 1.99.1 (58.13 s wall), on the node used for the QE references: 2 x Xeon Gold 6230 (20 cores each), Rocky Linux 9.7; node idle (load 0.4 to 1.8 before and after the batch), no concurrent jobs, three runs per case, single-core runs pinned with `taskset -c 6`. rust-dft runs on one thread (explicit sequential faer and serial FFT, recorded in `results.json`); QE 7.6 uses MKL 2026.1 with `OMP_NUM_THREADS=1`, once with one MPI process and once with four (as in R3 to R6). QE inputs are unedited copies of the verified ones (the Gamma 30 Ry input and the 4x4x4 shifted input `a45ed650...6600`), the same UPF (`e8d93375...5d68d5217`, provenance unverified).
- Wall times, mean of three (s): Gamma, 1 k point: rust-dft 2.28, QE 1 process 0.217, QE 4 processes 0.173 (rust-dft 10.5 times slower than one QE process, 13.2 than four). 4x4x4 shifted, 32 k points: rust-dft 94.3, QE 1 process 1.82, QE 4 processes 0.88 (51.9 and 107.1 times). Repetition spread under 1 percent. From 1 to 32 k points rust-dft grows 41.3 times and QE (one process) 8.4 times.
- Where rust-dft spends its time (its own phase timers): diagonalization is 91 percent of the engine time at Gamma (2.08 of 2.28 s) and 97 percent for 32 k points (91.8 of 94.5 s); density, potential, energy and mixing together are below 3 percent. Iterations: rust-dft 11 and 11, QE 11 and 8.
- Interpretation limits: this is the state before any optimization (the project rule was physics first); QE converges to `conv_thr` 1e-12 while rust-dft uses its default tolerances (density residual 1e-7, energy 1e-8 Ha) and certifies a fixed point with an extra map, so the times are not a like-for-like accuracy comparison; QE's own WALL lines were 0.12 s and 1.71 s (one process), the larger wall times include startup. No profile was taken, so the reasons for the gap in the diagonalization are not established; the per-k-point cost difference (QE reuses more work per k point) is a hypothesis to test by profiling, not a finding.
- Cross-platform: total energy on Linux x86 against macOS Apple Silicon: Gamma -7.295530468210918 against -7.295530468210927 Ha (difference +8.9e-15 Ha, a few ulp); 4x4x4 shifted -7.925976126648072 on both, bit-identical to all 16 printed digits. This supports same-source reproducibility across these two platforms to about 1e-14 Ha, but it is two cases only and says nothing about other compilers or CPUs.
- Toolchain note: `rust-toolchain.toml` pins 1.99.1 while the first install was `stable`; a `rustup show` call installed the pinned toolchain automatically (same compiler build, 617 MB more). The unused `stable` toolchain (1.5 GB) was left in place; everything is under `~/rustdft_p6/toolchain` (removal is `rm -rf ~/rustdft_p6/toolchain`).

## R8: PBE (GGA) comparison with QE 7.6 (code-to-code consistency, same LDA pseudopotential)
- Code: ticket P7B, ADR 0018: PBE exchange-correlation values and derivatives from `xcx` 0.3.1 (pure Rust, MPL-2.0 part for the libxc-derived functionals), plug-in glue in `dft-physics::gga` (spectral gradient with Nyquist zeroing, `V_xc = vrho - div(2 vsigma grad n)`, density threshold 1e-10, negative density rejected). The validation Si file was generated for LDA, so `[xc] functional = "pbe"` needs the explicit override `allow_functional_override = true`, which `results.json` and the console record (`override_applied=true`). QE received `input_dft = 'pbe'`, which it printed as "XC functional enforced from input: Exchange-correlation= PBE; any further DFT definition will be discarded". This tests that our PBE implementation agrees with QE's on identical input; the PBE result itself is not physically meaningful for an LDA pseudopotential.
- Library check before the SCF work (ADR 0018): xcx PBE against closed forms typed independently from the PBE96/PW92 papers at seven (n, sigma) points: energy per particle within 2.2e-16 Ha, derivatives against finite differences within 8e-11 relative.
- QE runs: pw.x 7.6, MPI x4, dev08 via the dft5090 session (relayed text blocks with sender-side sha256; the qe.out files stay on dev08; `tests/reference/pbe/<run>/qe_summary.json` holds the relayed values), inputs checked by sha256 before running (`7e249ba6...c034` Gamma, `d00bee3f...0da5` 2x2x2 shifted), same UPF (`e8d93375...5d68d5217`, provenance unverified), FFT grid 25^3, `conv_thr` 1e-12.
- Command: `python3 -I scripts/compare_kpoints_qe.py --ours tests/reference/pbe/<run>/ours_results.json --qe-summary tests/reference/pbe/<run>/qe_summary.json` (defaults: 1e-4 Ha/atom, 1e-3 Ha eigenvalues, 1e-6 weights), both exit 0.

| run | k points | iterations QE / ours | total energy ours - QE (Ha/atom) | largest eigenvalue difference (Ha) |
|---|---|---|---|---|
| PBE, Gamma, 30 Ry | 1 | 12 / 11 | +1.698e-6 | 1.27e-6 |
| PBE, 2x2x2 shifted, 30 Ry | 4 | 9 / 11 | +1.686e-6 | 1.82e-6 |

- Term by term (Gamma, ours minus QE, Ha): Hartree -1.9e-8, Ewald +1.8e-8, one-electron +2.55e-6, XC +8.4e-7. The one-electron offset is the known G=0 integral convention of R4 (2.55e-6 Ha for the LDA runs too). The XC difference is new: for the LDA runs XC agreed to 1.8e-8 Ha, for PBE it is 8.4e-7 Ha, which is the origin of the PBE offset being 4.2e-7 Ha/atom larger than the LDA offset (1.70e-6 against 1.28e-6 Ha/atom). It is 120 times below the energy tolerance and the eigenvalues still agree to 1.3e-6 Ha, but it is not explained: candidates are differences in how the density gradient is formed (our full FFT box with Nyquist zeroing against QE's treatment) and in density thresholds (ours 1e-10), neither of which was verified because QE's source was not read. It is recorded as a measured, unexplained, small difference.
- Limits: one cell, one pseudopotential used for a functional it was not generated for, Gamma and one k grid only, no 4x4x4, relayed text evidence (no whole-file hash verification), no density comparison, no forces, spin-unpolarized only; the iteration count differs (QE 12 against ours 11 at Gamma) with unequal convergence criteria.
- R7 update (2026-11-03, step O1 + O3, commit `f6e3679`, dev08 quiet node, `scripts/bench_scf.py`, median of 3, core 6, relayed attestation): wall times before and after, Gamma 2.278 s to 1.219 s (1.87x), 4 k points 10.97 s to 5.60 s (1.96x), 8 k points 21.24 s to 10.74 s (1.98x), 32 k points 98.98 s to 46.98 s (2.11x); total energy difference 0 and density files bit-identical for all four cases; 99 percent of the time saved was in the diagonalization. Against the QE 4-process times of R7 the factor changed from 13.2 to 7.0 (Gamma) and from 107 (by /usr/bin/time) to 53 (32 k points); QE was not re-measured in this batch. The benchmark-script medians are up to 5 percent higher than the earlier `/usr/bin/time` values because of the method. Later steps (Davidson product reuse, SCF tolerance schedule) are appended to `tests/reference/perf/dev08_baseline.json` and the README after they are pushed and measured on dev08; on the shared Mac the Davidson step took the 8-k-point case from 12.08 s to 2.06 s (indicative only).
- R7 update (2026-11-03, steps O4a/b at commit `3b94d53` and O4c at `ec77e7e`, dev08 quiet node, same protocol): wall times (s, median of 3) before optimization, after O1+O3, after the Davidson step, after the SCF schedule: Gamma 2.278, 1.219, 0.710, 0.427; 4 k points 10.97, 5.60, 2.81, 1.46; 8 k points 21.24, 10.74, 5.34, 2.67; 32 k points 98.98, 46.98, 22.41, 10.61. Cumulative speedup 5.3 (Gamma) to 9.3 (32 k points). Energy differences to the pre-optimization state at most 5.3e-15 Ha after the Davidson step and 1.5e-13 Ha after the SCF schedule (the gate is 1e-12 Ha; the headroom is now a factor 6.6, so the next arithmetic change must be checked against it); SCF iterations 11 became 10 with the schedule. QE 7.6 was re-measured in the same window (mean of 3, /usr/bin/time): Gamma 0.2133 s (1 process) and 0.1700 s (4 processes), 32 k points 1.8167 s and 0.8700 s, within 1.7 percent of the earlier values. Same-window ratios rust-dft over QE: Gamma 2.0 (1 process) and 2.5 (4 processes), 32 k points 5.8 and 12.2. The diagonalization share of the rust-dft time fell to 57 percent (Gamma) and 80 percent (32 k points). Relayed attestation; binaries sha256 in `tests/reference/perf/dev08_baseline.json`.
- R7 correction and equal-accuracy comparison (2026-11-03): the statement in earlier summaries that the QE reference runs were conducted under conditions favorable to QE was wrong: QE ran with `conv_thr` 1e-12 Ry, tighter than needed, which makes the times above favorable to rust-dft. QE `conv_thr` sweep on dev08 (same inputs, only that line changed, mean of 3, 1 process): Gamma 1e-6, 1e-8, 1e-10, 1e-12 gave 0.183, 0.190, 0.200, 0.210 s with 6, 7, 9, 11 iterations; 32 k points 1.003, 1.307, 1.463, 1.813 s with 4, 5, 6, 8 iterations. Total energy at 1e-8 equals the 1e-12 value to QE's printed 8 decimals (at 1e-6 it differs by 3e-8 Ry for Gamma and 4.0e-7 Ry for 32 k points). Against QE at 1e-8 the latest rust-dft step (O4c, 0.427 s and 10.61 s) is 2.2 times (Gamma) and 8.1 times (32 k points) slower; at 1e-6 2.3 and 10.6 times; four-process QE at 1e-8 gives 2.7 and 16. On the rust-dft side loosening the tolerances does not help: for the 4-k-point case, density, energy and eigen tolerances from 1e-8 down to 1e-5 changed the time from 0.55 s to 0.49 to 0.54 s with an energy change at most 4e-13 Ha (Mac, indicative), because the fixed-point certificate and the eigen-residual criterion still apply; the default already reaches an energy error of about 3e-14 Ha against a tight run. QE's `conv_thr` is a bound on an estimated energy error in Ry, rust-dft's gates are a density residual, an energy change and an eigen residual: the metrics are not identical, which is why the comparison above is made on achieved energy error and wall time.
- R7 update, deterministic k-point threads (2026-11-03, commit `1decf9c`, dev08 quiet node, one socket, `taskset -c 0-(T-1)`, `--threads T`, median of 3, relayed attestation, evidence `tests/reference/perf/dev08_threads.json`): wall times (s) for T = 1, 4, 8, 16: Gamma 0.403, 0.396, 0.396, 0.397; 4 k points 1.306, 0.702, 0.711, 0.713; 8 k points 2.370, 1.191, 1.016, 1.047; 32 k points 9.346, 4.270, 3.450, 3.057. The density file hashes are identical for T = 1, 4, 8, 16 in all four cases and the energy difference to the pre-optimization state is the same for every T (1.5e-13, 1.8e-15, 7.0e-14, 8.9e-15 Ha): the results are bit-identical for every thread count. The 32-k-point diagonalization went from 7.24 s to 0.85 s (8.5 times) but the engine time only from 9.39 s to 3.04 s: about 2.2 s (density accumulation, energy audit, per-k setup, certificate) is still sequential; on 4 k points the speedup stops at 4 threads because there are only 4 k points. QE 7.6 on the same node, 32 k points, `conv_thr` 1e-12, MPI processes 1, 4, 8, 16: 1.813, 0.873, 0.847, 0.697 s, so rust-dft is 5.2, 4.9, 4.1, 4.4 times slower with the same number of cores. The QE `conv_thr` 1e-8 runs at 8 and 16 processes had 31 and 19 percent spread (the node was not yet idle) and are not used; at 1 and 4 processes they gave 1.300 and 0.643 s.
- R7 update, parallel phases (2026-11-03, commit `20cb416`, ticket P7PA, same dev08 protocol, relayed attestation, evidence `tests/reference/perf/dev08_threads.json`): wall times (s) for T = 1, 4, 8, 16: Gamma 0.398, 0.392, 0.392, 0.395; 4 k points 1.295, 0.516, 0.539, 0.574; 8 k points 2.348, 0.819, 0.587, 0.650; 32 k points 9.270, 2.859, 1.751, 1.226. Density hashes identical across T and identical to the batch before P7PA (the k-ordered reductions keep the association); the energy differences did not change. For 32 k points and 16 threads the engine time fell from 3.043 s to 1.211 s (diagonalization 0.839 s, 69 percent of the run), a speedup of 80.7 against the pre-optimization 98.98 s. Against QE 7.6 with the same number of cores (`conv_thr` 1e-12): 5.1, 3.3, 2.1, 1.8 times slower for 1, 4, 8, 16 cores. With 4 and 8 k points the best thread count is 8 and 16 threads are slightly slower (few k points per thread; not investigated). The Gamma case is not parallelized yet (band-level threads: ticket P7PB in progress).
