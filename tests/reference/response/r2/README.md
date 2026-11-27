# R2 oracle protocol: dielectric tensor (eps_M, chi0 head) against external codes (dev08)

Status: inputs prepared 2026-11-27 (RPA-R2-CHI0); no external-code run
recorded yet. Numbers count as "Yambo / BerkeleyGW / QE reference" only after
a real dev08 run with the command, code version, input, UPF SHA-256 and
machine stored here. No UPF is committed. Plan:
`docs/roadmap/rpa-crpa-plan.md` ticket R2; conventions: ADR 0054.

## Required oracle software (decision D3: built on dev08 in user space through the master session; none is assumed to exist)

| Code | Role | Needed pieces | Notes |
| --- | --- | --- | --- |
| QE 7.6 `pw.x` | ground state and NSCF for the oracles | `pw.x` (exists on dev08) | same UPF, cutoff, FFT grid, k list as BAANDI |
| Yambo 5.x | primary RPA oracle, eps_M with and without local fields | `p2y`, `yambo` (optics / `em1s` static screening, `Chimod = "Hartree"`) | GPL: run as an external program only |
| BerkeleyGW 4.x | second RPA oracle | `pw2bgw.x` (QE), `epsilon.x` | needs registration; only if the user accepts it |
| QE `ph.x` | DFPT eps_infinity (includes the LDA/PBE xc kernel) | `ph.x` (`epsil = .true.`, `trans = .false.`) | matched only to a future BAANDI ALDA check, never to plain RPA; recorded, not an R2 gate |
| ABINIT 9.x (fallback) | RPA screening (`optdriver 3`) eps_M | `abinit` with the same UPF (or psp8 conversion, both hashes recorded) | only if Yambo cannot be built |

Before any oracle run: the master session builds the codes in user space
(no system change), records their versions and build flags here, and
writes start/done/failed/idle signal files for every dev08 request.

## Identical conditions

Same cell and coordinates, UPF files (SHA-256 recorded), `ecutwfc`
(`ecutrho = 4 ecutwfc`), FFT grid (`nr1..3` from `baandi grid`), SCF on the
same Gamma-centred mesh, NSCF on the FULL Gamma-centred mesh (explicit k
list in crystal coordinates, `nosym`, `noinv`), band count N_b = BAANDI
`settings.bands` (all bands enter the sum, no band window), response cutoff
`E_cut^chi` = `settings.ecut_chi_hartree` (Yambo `NGsBlkXd` = 2 E_cut^chi in
Ry; BerkeleyGW `epsilon_cutoff` = 2 E_cut^chi in Ry), static limit omega = 0,
q -> 0 with the nonlocal commutator in the dipoles (Yambo default
`DipApproach = "G-space v"` with the nonlocal term on), RPA (no xc kernel),
fixed occupations. Threads and ranks are free: accuracy comparison only, no
timing (decision D10).

## Gates (plan R2)

- eps_M vs Yambo and BerkeleyGW at matched settings: 0.5 percent relative,
  without and with local fields separately (diagonal of the tensor; for
  AlN both the in-plane xx = yy and the axial zz component).
- chi0 head tensor (`head_tensor_real` at omega = 0, Bohr^-1 Ha^-1, i.e.
  1 - 4 pi L = eps without local fields): 0.5 percent.
- DFPT-matched ALDA value vs ph.x: 1 percent (later, needs a BAANDI ALDA
  kernel; R2 records the ph.x number only).
- Experiment (Si 11.9, C 5.7, LiF about 1.9) is reported with the RPA/LDA
  overestimate stated, never a gate.
- Ground-state offset (BAANDI vs QE eigenvalues on the same mesh, R1
  comparison) is recorded separately for each case.

## Cases

| Case | Input | UPF | k (full mesh) | N_b | E_cut^chi |
| --- | --- | --- | --- | --- | --- |
| Si2 | `si2-444.toml`, `si2-888.toml` | Si.pz-vbc.UPF (LDA PZ) | 64 / 512 | 100 | 3 Ha |
| C diamond | `c2-444.toml` | PseudoDojo nc-sr-04 PBE standard C | 64 | 100 | 4 Ha |
| LiF | `lif-444.toml` | PseudoDojo PBE Li, F | 64 | 60 | 3 Ha |
| AlN wurtzite (tensor) | `aln-444.toml` | PseudoDojo PBE Al, N | 48 | 120 | 3 Ha |
| skewed Si2 | `si2-skewed-444.toml` | Si.pz-vbc.UPF | 64 | 100 | 3 Ha |

The R1 skewed Si2 geometry is (semi)metallic at LDA on 4x4x4 (BAANDI
indirect gap -0.2 mHa; the response driver rejects it, decision D6): a
milder strain must be chosen before it can serve as an R2 case.

Meshes and cutoffs are oracle-comparison settings, not converged values.

## Steps (dev08, through the master session)

```sh
B=target/release/baandi            # built from the R2 commit
R=tests/reference/response/r2
mkdir -p run && cp $R/si2-444.toml run/ && cp /path/Si.pz-vbc.UPF run/ && cd run
$B grid si2-444.toml | grep 'FFT dims'
$B dielectric si2-444.toml --output-dir baandi
# QE ground state + NSCF on the same full mesh (reuse ../r1/make_qe_r1.py for
# the k list and grid: it reads only BAANDI's JSON), then per oracle:
#   Yambo:  p2y; yambo -i -V RL; yambo -o c (or em1s), NGsBlkXd, BndsRnXd = 1..N_b,
#           Chimod "Hartree" (LF) and "IP" (no LF), q -> 0 along x, y, z
#   BGW:    pw2bgw.x; epsilon.x with epsilon_cutoff, number_bands = N_b, q0 small
#   ph.x:   epsil = .true. (record eps_infinity only)
```

Store per case: BAANDI `epsilon.json` and `epsilon_M.csv`, the oracle
inputs and outputs, code version lines, UPF SHA-256, host, exact commands
and a `*-compare.json` (eps_M diagonal per direction, with and without local
fields, relative differences, pass/fail). Do not store wavefunction
databases or `*.save` directories.

Also record on dev08 (BAANDI only): an MPI rank series 1/2/4/8 of
`si2-444.toml` (`--features mpi`; gate: eps tensors within 1e-12 relative
between rank counts, np = 1 bit-identical to the threads-only run) and
`band_solver = "chebyshev"` against dense (gate 1e-8 relative on eps).

## BAANDI observations on the Mac (2026-11-27, not oracle evidence)

Uncommitted worktree RPA-R2-CHI0 on main `ee8779dc`, `cargo build --release
-p dft-cli`, Apple M6, threads default, one rank. Si.pz-vbc.UPF SHA-256
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`;
PseudoDojo files: SHA-256 in each `epsilon.json` (`pseudopotentials`).
Command: `baandi dielectric <input> --output-dir out`. Records:
`mac-2026-11-27/`.

| Case | N_G | bands | eps_inf with LF | eps_inf without LF | max NSCF residual | wall |
| --- | --- | --- | --- | --- | --- | --- |
| Si2 4x4x4 | 65 | 100 | 22.6814 | 24.8945 | 7.0e-14 Ha | 40 s |
| Si2 8x8x8 | 65 | 100 | 13.4268 | 14.8668 | 6.0e-14 Ha | 250 s |
| C 4x4x4 (PBE, 80 Ry, E_cut^chi 4 Ha) | 27 | 100 | 7.1525 | 7.4393 | 1.3e-13 Ha | 29 s cpu |
| LiF 4x4x4 (PBE, 60 Ry) | 27 | 60 | 1.9988 | 2.0893 | 9.1e-14 Ha | 24 s cpu |
| AlN wurtzite 4x4x3 (PBE, 60 Ry), xx = yy / zz | 67 | 120 | 4.5744 / 4.5431 | 4.9347 / 4.9010 | 9.0e-14 Ha | 279 s cpu |

Band convergence, Si2 4x4x4 (eps_inf with / without LF): 30 bands
22.6097 / 24.8906, 50 bands 22.6582 / 24.8933, 100 bands 22.6814 / 24.8945,
200 bands 22.6848 / 24.8947 (monotone, as every added transition lowers
chi0). The 4x4x4 Gamma-centred mesh overweights the small direct gaps
near Gamma; 8x8x8 brings eps_inf to 13.4 (LF) / 14.9 (no LF), the order
of published RPA-LDA values for Si (no source number is quoted here; the
matched Yambo/BerkeleyGW runs above are the comparison, and a k-mesh series
on dev08 is needed before any convergence statement). Experiment 11.9.
