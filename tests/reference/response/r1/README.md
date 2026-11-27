# R1 oracle protocol: many-band NSCF against QE 7.6 `pw.x` nscf (dev08)

Status: inputs prepared 2026-11-27 (RPA-R1-MANY-BANDS); no QE run recorded
yet. Numbers count as "QE reference" only after a real dev08 run with the
command, QE version, input, UPF SHA-256 and machine stored here. No UPF is
committed. Plan: `docs/roadmap/rpa-crpa-plan.md` ticket R1.

Gate (plan R1): BAANDI vs QE eigenvalues <= 2e-5 Ha for every band up to
N_b at every k, after removing one constant offset (highest occupied
eigenvalue at the first k; the ground-state/convention offset is recorded
separately as `offset_ha`). Plane-wave counts per k must agree exactly.

## Identical conditions

Same cell and coordinates, UPF file (SHA-256 recorded by both inputs),
`ecutwfc` (ecutrho = 4 ecutwfc), FFT grid (`nr1..3` = BAANDI's grid from
`baandi grid`), explicit k list (BAANDI's Gamma-centered full mesh in crystal
coordinates, `nosym`, `noinv`), band count `nbnd` = `[nscf] bands`, fixed
occupations, `diago_full_acc = .true.`, `conv_thr = 1e-12` (QE nscf eigenvalue
threshold then about 1e-14 Ha/electron). Thread/rank counts are free: this is
an accuracy comparison, no timing is compared (decision D10).

## Cases

| Case | Input | UPF | k (full mesh) | N_b | Notes |
| --- | --- | --- | --- | --- | --- |
| Si2 4x4x4 | `si2-444.toml` | Si.pz-vbc.UPF (LDA PZ) | 64 | 400 | Mac smoke run: N_PW 749 max, grid 25^3, 64 k dense in about 1 min |
| Si2 8x8x8 | `si2-888.toml` | Si.pz-vbc.UPF | 512 | 400 | plan cost: under 5 min dense on 8 cores |
| LiF, AlN, skewed Si | build from `si2-444.toml` | PseudoDojo PBE Li/F, Al/N (nc-sr-04 standard) | 2x2x2 / Gamma+A | 40-80 | dev08 item: cutoffs per the dataset hints; BAANDI internal gates already pass on the Mac at low cutoff |

## Steps (dev08, through the master session, with start/done/failed/idle signal files)

```sh
B=target/release/baandi            # built from the R1 commit
R=tests/reference/response/r1
mkdir -p run && cp $R/si2-444.toml run/ && cp /path/Si.pz-vbc.UPF run/ && cd run
$B grid si2-444.toml | grep 'FFT dims'           # e.g. 25 x 25 x 25
$B nscf si2-444.toml --output-dir baandi --checkpoint-dir cp
python3 -I ../$R/make_qe_r1.py baandi/nscf.json --nr 25 25 25 --scf-kgrid 4 4 4 --prefix si2_444 --dir .
mpirun -np 8 pw.x -in si2_444-scf.in  > si2_444-scf.out
mpirun -np 8 pw.x -in si2_444-nscf.in > si2_444-nscf.out
python3 -I ../$R/compare_r1.py baandi/nscf.json qe-work/si2_444.save/data-file-schema.xml --json si2_444-compare.json
```

Repeat with `si2-888.toml` (`--scf-kgrid 8 8 8`, prefix `si2_888`). Store
here per case: the BAANDI `nscf.json` (and `nscf_eigenvalues.csv`), the two
QE inputs and outputs (`*.out`), `*-compare.json`, `pw.x` version line, UPF
SHA-256, host, and the exact commands. Do not store `qe-work/` or `cp/`.

Also record on dev08 (BAANDI only): `band_solver = "chebyshev"` on
`si2-444.toml` against the dense `nscf.json` (gate 1e-8 Ha, all certified
bands) and an MPI rank series 1/2/4/8 of the dense and Chebyshev runs
(`--features mpi`, gate 1e-10 Ha between rank counts; dense is expected
bit-identical).

## Tools

- `make_qe_r1.py`: writes `PREFIX-scf.in` / `PREFIX-nscf.in` from BAANDI's
  `nscf.json` (reads no UPF contents). Python 3.9+, standard library only.
- `compare_r1.py`: parses QE's `data-file-schema.xml` (`ks_energies`,
  Hartree) and compares per k and band; exit 0 only when npw matches and the
  aligned maximum is within `--gate` (default 2e-5 Ha). Parser path checked
  on the Mac against a synthetic XML built from a BAANDI run (offset 0.01 Ha
  recovered, aligned max 8.9e-16 Ha); not yet run on a real QE output.
