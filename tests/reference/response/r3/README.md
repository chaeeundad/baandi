# R3 oracle protocol: RPA correlation energy (ACFDT) against ABINIT (dev08)

Status: inputs and scripts prepared 2026-11-27 (RPA-R3-ACFDT); no
external-code run recorded yet. A number counts as "ABINIT reference" only
after a real dev08 run with the command, code version and build options,
input file, pseudopotential SHA-256 (UPF and, if converted, psp8) and machine
stored here. No pseudopotential file is committed. Plan:
`docs/roadmap/rpa-crpa-plan.md` ticket R3; conventions: ADR 0055 (on top of
ADR 0054).

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
