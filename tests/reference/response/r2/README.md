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
