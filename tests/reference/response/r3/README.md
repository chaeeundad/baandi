# R3 oracle protocol: RPA correlation energy (ACFDT) against ABINIT (dev08)

Status: **finite-q part validated against ABINIT 10.8.3 on dev08,
2026-11-28; q = 0 term reported separately and shown to be a convention
that vanishes as 1/N_q** (inputs prepared 2026-11-27, RPA-R3-ACFDT; q image
amendment R3-QIMAGE). Record: `dev08-2026-11-28/summary.json` (BAANDI tree
gen75 `be57c04c` and binary SHA-256, ABINIT SHA-256, input, UPF and
`case.abi` SHA-256 per case, commands; raw records on dev08 under
`/scr/mino/rustdft/r75/r3`). No pseudopotential file is committed. Plan:
`docs/roadmap/rpa-crpa-plan.md` ticket R3; conventions: ADR 0055 (on top of
ADR 0054).

Matched conventions: BAANDI `[rpa] response_sphere = "G"` (one |G| sphere
for every q, v = 4 pi / |q+G|^2) with each mesh q evaluated at its
first-Brillouin-zone image (ADR 0055 item 2), which is ABINIT's
convention; ABINIT `gwcalctyp 1`, `nfreqim 32`, `ixc -101130` (libxc PBE),
one dataset per (cutoff, IBZ q) for the per-q comparison. Si uses the
PseudoDojo PBE file (`si2-dojo-444`), because ABINIT 10.8.3 cannot read
Si.pz-vbc.UPF. The q = 0 term is not part of the gate: BAANDI averages the
k.p head and wings over directions, ABINIT uses its auxiliary-function
treatment (`gw_icutcoul 6`, default); ABINIT's own alternatives 3/7/16
change its Si q = 0 term by at most 0.5 %, so the difference is not an
ABINIT setting choice.

Gate: finite-q sum (all q != 0) within 1 meV/atom at each matched cutoff.

| Case | E_cut^chi (Ha) | finite-q sum diff (meV/atom) | largest single-q diff | q = 0 term diff, weight 1/N_q (meV/atom) |
| --- | --- | --- | --- | --- |
| Si2 dojo 4x4x4 | 6 / 8 / 10 / 12 | -0.004 / -0.006 / -0.010 / -0.014 | <= 0.004 | -104.69 at every cutoff |
| LiF 4x4x4 | 4 / 5 / 6 / 8 | -0.065 / -0.062 / -0.063 / -0.063 | <= 0.025 | -9.15 to -9.53 |
| C 4x4x4 | 4 / 5 / 6 / 8 | +0.273 / -0.005 / +0.037 / +0.010 | <= 0.112 | -62.82 to -63.61 |

All twelve points pass. q mesh series, Si dojo 6 Ha (SCF 4x4x4 in both
codes): total difference -104.7 / -21.9 / -8.8 meV/atom at 4x4x4 / 6x6x6 /
8x8x8; difference x N_q = -0.49 / -0.35 / -0.33 Ha per cell, tending to
about -0.33 Ha: the remaining difference falls as 1/N_q, a pure q = 0
convention (ABINIT q = 0 terms at 6x6x6 and 8x8x8 inferred as total minus
BAANDI's finite-q part, which matched to 0.014 meV/atom at 4x4x4).
Frequency grid: ABINIT nfreqim 16/32/64 spread at 6 Ha below 0.001
meV/atom. BAANDI-only MPI rank series (tree `03aab866`, Si.pz-vbc, 4 Ha,
4x4x4): np 2/4/8 against np 1 within 1.0e-13 Ha (gate 1e-9 Ha). Not in
this record: extrapolated-cutoff comparison (2 meV/atom gate), a BAANDI
frequency series on these inputs (Mac series: `mac-2026-11-27/`), Ne/Ar,
N2 in a box, skewed Si (later R3 items).

Source note: the ABINIT conventions (G sphere, first-zone q image, q = 0
treatment and `gw_icutcoul` choices) were established from numbers and
partly from the dev08 session's reading of the ABINIT source, for oracle
matching only; recorded here in prose, no oracle code in this repository.
The README results block is generated from the summary by
`scripts/make_readme_figures.py`.

## Required oracle software (decision D3: built on dev08 in user space through the master session; nothing is assumed to exist)

| Code | Role | Needed pieces | Notes |
| --- | --- | --- | --- |
| ABINIT >= 9.10 (10.x preferred) | primary R3 oracle | `abinit` with the GW/screening driver (`optdriver 3`) and `gwrpacorr` | GPL: run as an external program only; reads NC UPF (UPF1/UPF2) directly, or psp8 (PseudoDojo) |
| QE 7.6 ACFDT | optional second oracle | an ACFDT/RPA tool in the dev08 QE build | availability not verified; not assumed |

ABINIT build (user space, no system change), to be recorded here with the
exact `configure` line and `abinit --build` output:

- compilers: GNU Fortran/C (or Intel oneAPI) with MPI (Open MPI as used by
  BAANDI on dev08), OpenMP optional;
- mandatory libraries for ABINIT 9/10: LibXC (5.x/6.x), HDF5, NetCDF-C and
  NetCDF-Fortran (built in user space if absent);
- linear algebra: MKL (as on dev08) or OpenBLAS + LAPACK; FFT: MKL DFTI or
  FFTW3;
- GW is part of the default build; confirm `optdriver 3` and `gwrpacorr`
  are accepted by running the tiny Si case below with `nband` small first.

Before every run: start/done/failed/idle signal files on the shared NFS
directory (feedback rule), one request at a time through the master session.

## Identical conditions

Same cell and coordinates (`lattice_columns_bohr`, `structure` in
`rpa.json`), same pseudopotential file (SHA-256 in `rpa.json`
`pseudopotentials`; a psp8 conversion records both hashes), wavefunction
cutoff `ecut` = `cutoff_hartree`, FFT grid `ngfft` = `settings.fft_grid`,
Gamma-centred k mesh `ngkpt` = `settings.k_mesh` with `shiftk 0 0 0`,
fixed occupations, band count `nband` = `settings.bands_per_cutoff` per
cutoff (all bands enter chi0), response cutoff `ecuteps` =
`ecut_chi_series_hartree` per point, nonlocal commutator in the q -> 0
dipoles (`inclvkb 2`), no xc kernel (RPA), q = differences of the k mesh.
Frequencies: both codes use their own imaginary grids; each side is
converged separately in the number of points (BAANDI: `points` 16 -> 24 ->
32 must change E_c by <= 0.05 meV/atom; ABINIT: `nfreqim` likewise), so the
grids are not a source of difference above 0.05 meV/atom.

Known convention differences to check in the ABINIT log, not to tune away:

1. Response sphere at finite q. BAANDI default: `|q+G|^2/2 <= E_cut^chi`.
   If the ABINIT log prints the same `npweps` for every q, ABINIT uses one
   `|G|` sphere for all q: rerun BAANDI with `[rpa] response_sphere = "G"`
   (one Miller set, v still `4 pi/|q+G|^2`) for the matched comparison and
   record which convention was used. In G mode E_c(q) depends on the q
   vector: BAANDI evaluates every mesh q at its first-Brillouin-zone image
   (mean over tied images on zone faces; ADR 0055 item 2, R3-QIMAGE) and
   `rpa_q.csv` prints that vector. A per-q comparison must use the same
   vector on the ABINIT side (ABINIT's printed q may be a different G0
   image; only the shortest-image value is comparable).
2. q -> 0. BAANDI averages the k.p head/wings limit over a Lebedev set;
   ABINIT evaluates the head at a small q along `gw_qlwl`. For the cubic
   cases (Si, C, LiF) eps_M is isotropic and the two agree in the limit;
   the BAANDI record shows the direction spread (`q0_direction_check`).
   Anisotropic cases are compared only after this is resolved.
3. Ground state: ABINIT runs its own SCF from the same pseudopotential;
   record the total-energy and eigenvalue offsets (R1 method) next to the
   E_c difference.

## Gates (plan R3)

- `|E_c(BAANDI) - E_c(ABINIT)| <= 1 meV/atom` at each matched cutoff
  (same bands, k mesh, response sphere convention);
- extrapolated values (`(E_cut^chi)^{-3/2}` fit of the same cutoffs on both
  sides) within 2 meV/atom;
- frequency-grid error <= 0.05 meV/atom on both sides;
- BAANDI rank-count independence within 1e-9 Ha (dev08, `--features mpi`,
  np = 1/2/4/8 of `si2-444.toml`; np = 1 bit-identical to threads).

## Cases

| Case | Input | Pseudopotential | k | E_cut^chi series (Ha) | bands |
| --- | --- | --- | --- | --- | --- |
| Si2 | `si2-444.toml` | Si.pz-vbc.UPF (LDA PZ) | 4x4x4 | 6, 8, 10, 12 (32 points; Mac 4-8 failed the linearity check; 16 Ha would exceed the smallest k basis of 30 Ry) | auto = N_G |
| C diamond | `c2-444.toml` | PseudoDojo nc-sr-04 PBE standard C | 4x4x4 | 4, 5, 6, 8 | auto |
| LiF | `lif-444.toml` | PseudoDojo PBE Li, F | 4x4x4 | 4, 5, 6, 8 | auto |

Later in R3 (plan): Ne or Ar fcc, N2 in a box (atomization needs the N atom
with spin: collinear `[spin]`, lda-pw or pbe), a skewed Si cell (the R1/R2
skewed geometry is metallic at LDA on 4x4x4 and rejects; a milder strain is
needed), and k-mesh series n, n+2, n+4.

## Steps (dev08, through the master session)

```sh
B=target/release/baandi             # built from the R3 commit
R=tests/reference/response/r3
mkdir -p run && cp $R/si2-444.toml run/ && cp /path/Si.pz-vbc.UPF run/ && cd run
$B rpa si2-444.toml --output-dir baandi --checkpoint-dir baandi-q
python3 -I ../$R/make_abinit_r3.py baandi/rpa.json --nfreqim 16 --out si2_444.abi
abinit si2_444.abi > si2_444.log 2>&1      # ABINIT 9.10+/10.x, record `abinit --version`
# Write abinit-si2_444.json from the log (E_c per screening dataset, quoted lines),
# then:
python3 -I ../$R/compare_r3.py baandi/rpa.json abinit-si2_444.json --out si2_444-compare.json
```

Series runs use `points = 32` (Mac: with the series-wide auto omega0, 16
points were 0.2 meV/atom off). Frequency convergence on both sides: BAANDI `points` 16, 24, 32 (same
input, `[response.frequency] points`), ABINIT `nfreqim` 16, 24, 32.

Store per case: BAANDI `rpa.json`, `rpa_q.csv`, `rpa_extrapolation.csv`,
the ABINIT input, log excerpt with the E_c lines, `abinit --version` and
build line, pseudopotential hashes, host, exact commands and the
`*-compare.json`. Do not store WFK/SCR files.

## BAANDI observations on the Mac (2026-11-27, not oracle evidence)

See `mac-2026-11-27/README.md` (Si.pz-vbc, small cutoff series, frequency,
direction, band, q-symmetry and k-mesh checks).
