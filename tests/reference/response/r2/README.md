# R2 oracle protocol: dielectric tensor (eps_M, chi0 head) against external codes (dev08)

Status: **validated against Yambo 5.3.0 and ABINIT 10.8.3 on dev08,
2026-11-28** (inputs prepared 2026-11-27, RPA-R2-CHI0). Records:
`dev08-2026-11-28/summary.json` (matched run, BAANDI tree gen75 `be57c04c`,
binary, Yambo/p2y/pw.x/ABINIT SHA-256, input and UPF SHA-256 per case,
ABINIT input hashes) and `dev08-2026-11-27/summary.json` (first run, BAANDI
`full` against Yambo, plus the BAANDI-only rank and Chebyshev checks). No
UPF is committed. Plan: `docs/roadmap/rpa-crpa-plan.md` ticket R2;
conventions: ADR 0054. BerkeleyGW was not run (Yambo plus ABINIT give two
independent oracles); ph.x (ALDA) is not an R2 gate.

Gate 0.5 % relative on eps_inf with and without local fields, every
diagonal component. Matched results (max relative difference over x, y, z
and LF / no LF):

| Case | BAANDI `full` eps_inf LF / no LF | `full` vs Yambo (11-27) | `no-g0-all` vs Yambo | `no-g0` vs ABINIT | gate |
| --- | --- | --- | --- | --- | --- |
| Si2 4x4x4 | 22.6814 / 24.8945 | 1.1e-5 | 8.0e-6 | (ABINIT cannot read Si.pz-vbc) | pass |
| Si2 8x8x8 | 13.4268 / 14.8668 | 1.8e-5 | not run | not run | pass (`full`) |
| C 4x4x4 | 7.1525 / 7.4393 | 4.3e-5 | 6.9e-6 | not run | pass |
| LiF 4x4x4 | 1.9988 / 2.0893 | 3.2e-3 | 1.9e-5 | 5.0e-5 | pass |
| AlN 4x4x3, xx = yy / zz | 4.5744 / 4.9347, 4.5431 / 4.9010 | 5.2e-3 | 1.3e-5 | 1.4e-5 | pass |

Convention (Gamma, k + G = 0 nonlocal velocity): Yambo drops the whole
k + G = 0 plane wave from [V_NL, r], reproduced by BAANDI's diagnostic
`BAANDI_RESPONSE_VELOCITY=no-g0-all`; ABINIT drops the l = 1 projector
gradient there, reproduced by `no-g0`. BAANDI's default `full` keeps the
term, which is part of the exact d H / d k (finite-difference proof at
k = 0, section "Velocity at Gamma" below, R2-GAMMA-VEL), and is the
physical result; the `full` LiF/AlN differences of 0.3 to 0.5 % are the
oracles' Gamma treatment on these coarse meshes. ABINIT observable: its
printed dielectric constant with / without local fields at q_red
(1, 2, 3)e-5 (4 to 5 digits), against BAANDI's tensor projected on that
direction. Diagnostics: BAANDI `full` dipoles equal Yambo `ndb.dipoles` to
2e-5 at every k except Gamma; per-band eigenvalues against QE: LiF
1.2e-5 Ha, AlN 2.4e-6 Ha. BAANDI-only (11-27): MPI np 1 bit-identical to
threads, np 2/4/8 within 6.2e-14 relative, Chebyshev against dense
7.0e-10 relative. The skewed Si2 4x4x4 case is rejected as noninsulating
(gap -1.97e-4 Ha, decision D6), as expected.

Source note: the ABINIT and Yambo conventions above were established from
the numbers (sections below) and partly from the dev08 session's reading of
the ABINIT source, for oracle matching only; this repository records them
in prose, contains no oracle code, and nothing in BAANDI's default path was
changed toward the oracles.

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
`E_cut^chi` = `settings.ecut_chi_hartree` (Yambo `NGsBlkXd` = BAANDI
`response_size` in RL units, not Ry: Yambo rounds an energy value up to the
next G shell, dev08 2026-11-27 Si2: 6 Ry became 89 G against BAANDI's 65; QE
nscf with `nbnd = N_b + 8` per the R1 lesson and Yambo `BndsRnXd` =
`DipBands` = 1..N_b; prefer N_b at a degeneracy closure, Yambo warns when N_b
splits a degenerate level; BerkeleyGW `epsilon_cutoff` = 2 E_cut^chi in Ry), static limit omega = 0,
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

## Velocity check R2-VEL-CHECK (Mac, 2026-11-27, BAANDI-internal, not oracle evidence)

`crates/dft-response/tests/velocity_fd.rs`: `<psi_i| dH_k/dk_b |psi_j>`
(i == j and i != j, x/y/z, several k) against central differences of H_k
at fixed Miller indices, H from the production `NonlocalOperator` rows.
Nonlocal part, Richardson (h, h/2), best over h: synthetic s/p/d/f with two
coupled projectors per l <= 2 on a skewed cell 8.9e-11; PseudoDojo LiF 30 Ha
6.2e-11; PseudoDojo AlN 30 Ha 9.0e-10 (plain central difference shows
clean h^2 convergence). The velocity is d H_k / d k of BAANDI's Hamiltonian.

Size of the nonlocal commutator (switch `BAANDI_RESPONSE_VELOCITY`, same
inputs and UPFs as above; eps_inf at omega = 0):

| Case | full (LF / no LF) | local only | l = 0, 1 nonlocal only (`sp`) |
| --- | --- | --- | --- |
| LiF (avg) | 1.99876 / 2.08929 | 2.40301 / 2.50069 | 2.00081 / 2.09183 |
| AlN xx | 4.57445 / 4.93470 | 5.64332 / 6.04854 | 4.60179 / 4.96448 |
| AlN zz | 4.54312 / 4.90099 | 5.58710 / 5.99460 | 4.56912 / 4.92919 |

LiF band series (full, LF / no LF avg): 60 1.99876 / 2.08929, 61 1.99879 /
2.08940, 64 1.99893 / 2.08972, 68 1.99910 / 2.09011, 100 1.99967 / 2.09146.

## Velocity at Gamma, k + G = 0 (R2-GAMMA-VEL, Mac 2026-11-28)

dev08 found BAANDI's full velocity equal to Yambo 5.3.0 `ndb.dipoles` to
2e-5 at every k except k = 0 (BAANDI Gamma sum of |p|^2 +1.7 % LiF, +1.9 %
AlN; LiF VBM[3-5] -> CBM[6] r 3.84 vs 2.97 Bohr^2). The only plane wave
special to Gamma is G = 0 with k + G = 0.

Analysis. The KB factor of a projector of angular momentum l is
`beta_l(|q|) Y_lm(q^)`, q = k + G, with `beta_l(q) = q^l (a_0 + a_2 q^2 + ...)`:
a solid harmonic (homogeneous polynomial of degree l in q) times an even
series, hence smooth at q = 0. Its q gradient at q = 0 vanishes for l = 0
and l >= 2 and equals `beta_1'(0) Y_1m(e_b)` for l = 1 (finite, direction
dependent). `projector_k_gradient` uses exactly that limit (value 0 for
l >= 1, `beta_0(0)` for l = 0). So at Gamma the commutator has a G = 0 piece
`B(0) D G_b(0)^H + G_b(0) D B(0)^H` that connects an s-like state with a
G = 0 coefficient to the p projectors of the atoms.

Checks (BAANDI-internal, all at k = 0 exactly and the full matrix, i != j
included):

| Check | max rel err vs central FD of H (best h) | with the k+G = 0 gradient dropped |
| --- | --- | --- |
| synthetic s/p/d/f (`velocity_fd`, debug) | 4.6e-9 | 3.5e-2 |
| PseudoDojo LiF 30 Ha (`velocity_fd`, release) | 1.2e-10 | 5.7e-3 |
| PseudoDojo AlN 30 Ha (`velocity_fd`, release) | 4.5e-10 | 1.7e-2 |
| LiF, production NSCF H, complete EVD, 24 bands (`response_velocity_hf gamma`), Richardson | 9.6e-11 | 4.2e-2 |
| AlN, same, 24 bands, Richardson | 4.1e-11 | 1.2e-1 |

Continuity: p(k = 1.8e-6 Bohr^-1) - p(0) with the same vectors is 3e-6
relative (synthetic), 3.0e-6 (LiF), 1.6e-6 (AlN): the velocity is smooth
through Gamma. BAANDI's Gamma velocity is d H / d k of its Hamiltonian;
nothing is changed in the default path.

What the oracles do (inferred from numbers only; no oracle source was read
for this Mac analysis. Annotation 2026-11-28: the dev08 session later read
parts of the ABINIT source for oracle matching only; see Status).
`BAANDI_RESPONSE_VELOCITY=no-g0` drops the k+G = 0 projector gradient,
`no-g0-all` drops the whole k+G = 0 plane wave from the commutator (also the
l = 0 projector value there). Same inputs, UPFs and settings as above;
record `mac-2026-11-28-gamma-vel/summary.json`.

| Case | eps_M LF / no LF: BAANDI full | no-g0 | no-g0-all | Yambo 5.3.0 | ABINIT 10.8.3 (q dir) |
| --- | --- | --- | --- | --- | --- |
| LiF (xx) | 1.998757 / 2.089293 | 1.995196 / 2.085782 | 1.992497 / 2.083124 | 1.992459 / 2.083084 | 1.9951 / 2.0857 |
| AlN xx | 4.574447 / 4.934704 | 4.555187 / 4.915699 | 4.554345 / 4.914863 | 4.554305 / 4.914821 | |
| AlN zz | 4.543120 / 4.900993 | 4.520322 / 4.878276 | 4.519763 / 4.877702 | 4.51971 / 4.877639 | |
| AlN along ABINIT q | 4.56588 / 4.92549 | 4.54566 / 4.90547 | 4.54489 / 4.90471 | | 4.5456 / 4.9054 |
| Si2 | 22.68135 / 24.89449 | 22.68176 / 24.89426 | 22.68129 / 24.89462 | 22.68111 / 24.89444 | |
| C | 7.152542 / 7.439297 | 7.152410 / 7.439281 | 7.152280 / 7.439388 | 7.152231 / 7.43934 | |

Gamma dipoles, full / no-g0-all: sum |p|^2 ratio 1.0177 (LiF), 1.0194 (AlN)
(dev08: +1.7 %, +1.9 %); LiF [3-5] -> [6] r 3.8411 / 2.9716 Bohr^2 (Yambo
2.97); AlN [6-7] -> [9] 17.631 / 15.900 and [8] -> [9] 9.756 / 8.789
(+10.9 %, +11.0 %; dev08 +11 %). Every other k is unchanged (ratio 1.00000).

Conclusion: with `no-g0-all` BAANDI reproduces Yambo to <= 2e-5 relative
on all four cases (LiF 1.9e-5, AlN 1.3e-5, Si 8e-6, C 7e-6), and with
`no-g0` it reproduces the ABINIT LiF and AlN values to their printed digits
(1e-4). Both treatments omit a term that is part of d H / d k at Gamma (the
FD checks above fail without it), so BAANDI's default (`full`) stays and the
R2 LiF/AlN differences to Yambo and ABINIT (0.3 to 0.5 %) are the oracles'
k + G = 0 treatment at Gamma on these coarse 4x4x4 / 4x4x3 meshes (Gamma
weight 1/64, 1/48), not a BAANDI error. The effect shrinks with the k mesh
(one k point of N_k); Si and C move by <= 4e-5 because their dominant Gamma
transitions connect p-like triplets, whose G = 0 coefficients vanish by
symmetry. Nothing was tuned
toward the oracles; the two diagnostic values are never production settings.

## Dipole dump for external comparison (R2-VEL-CHECK diagnostic)

`BAANDI_RESPONSE_DIPOLE_DUMP=gamma|all baandi dielectric <input>` (default
off; eps and every other output unchanged) also writes `dipoles.csv`:
`gamma` = mesh k index 0 (Gamma) and 1 (the first other k of the
Gamma-centred mesh, k-major order i, j, l with the last index fastest),
`all` = every k. One row per (V, C) pair of degenerate groups at one k: V
an occupied group, C an empty group, a group = consecutive bands whose
eigenvalue spacing is below 1e-5 Ha (chained). Gauge-invariant sums only.

| Column | Meaning |
| --- | --- |
| `k_index` | mesh index (0 = Gamma) |
| `k1,k2,k3` | k in the reciprocal-lattice basis (fractional, wrapped to [-1/2, 1/2)) |
| `kx,ky,kz` | k Cartesian (Bohr^-1), input lattice frame |
| `v_first,v_last,c_first,c_last` | band ranges, 1-based, inclusive, ascending energy |
| `e_v_ha,e_c_ha` | group mean eigenvalues (Ha) |
| `c_at_band_limit` | 1 if C ends at the last computed band (group possibly incomplete) |
| `s_x,s_y,s_z` | `sum_{v in V, c in C} abs(<u_v| dH_k/dk_a |u_c>)^2` (a.u.: (Ha Bohr)^2) |
| `r_x,r_y,r_z` | `sum abs(<u_v|dH_k/dk_a|u_c>)^2 / (e_c - e_v)^2 = sum abs(<v| r_a |c>)^2` (Bohr^2), per-band energies |
| `s_local_*`, `r_local_*` | the same with the kinetic velocity k + G only (no [V_NL, r]) |

Conventions: Hartree atomic units (hbar = m_e = 1, so dH/dk = velocity =
momentum for the local part); Cartesian a = x, y, z of the input frame;
orbitals cell-normalized (sum_G |c|^2 = 1); spin-unpolarized, one spin
channel, no factor 2 for spin and no k weight; `dH_k/dk = (k+G) + dV_NL/dk`
is the matrix of `-i[r, H]`, so `<v|r|c> = <v|dH/dk|c> / (i (e_c - e_v))`
up to a phase. The full velocity is what enters the BAANDI head
(`chi0_00 = q^2 sum 2 |q^.<v|r|c>|^2 ... / (e_v - e_c)` at omega = 0).

Mapping to Yambo `ndb.dipoles` (check the variable names of the installed
version with `ncdump -h`): Yambo stores `DIP_iR(c, v, k, spin) = i <c|r|v>`
(from the velocity divided by the energy difference, including [V_NL, r]
with `DipApproach = "G-space v"`), and depending on the version also
`DIP_P` (<c|p|v>, kinetic) and `DIP_v` (velocity incl. nonlocal). Compare
`r_a` with `sum_{v,c} |DIP_iR_a|^2` over the same groups (Yambo energies:
the QE eigenvalues, so the group partition may differ at the 1e-5 Ha level
near accidental degeneracies), `s_a` with `sum |DIP_v_a|^2`, and
`s_local_a` with `sum |DIP_P_a|^2`. Yambo's Cartesian axes are QE's
(the same input frame when the QE input uses the same cell vectors); its k
list is in its own order (match by k coordinates; Yambo internally uses
units of 2 pi / alat). Mac LiF example (Gamma and k index 1 = (0, 0, 1/4)):
`mac-2026-11-27/lif-444-dipoles.csv`.
