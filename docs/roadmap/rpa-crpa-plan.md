# RPA and cRPA: design and validation plan (plan only, not implemented)

Status 2026-11-26: proposal under scope decision U11 (AGENTS.md). No code before
2026-12-01, and no stage starts before its ticket and the open decisions of
section 7 that it depends on are answered. Implementers are Codex; Claude
architects, reviews and runs gates; dev08 measurements go through the master
session. No QE source was read for this plan; it relies on public papers and
public documentation of the codes it names as oracles.

Anchors in the current code (read 2026-11-26):

- `crates/dft-engine/src/refscf/nscf.rs`: `nscf`/`nscf_with` solve a fixed
  supplied density at an ordered k list, every rank on every k in sequence,
  every band strict, a fresh `||H psi - eps psi||` certificate, at most 64
  restarts, at least N/2 + 4 bands. Gamma uses complex orbitals. Eigenvectors
  come back as owned plane-wave rows (`Bands::rows`).
- `crates/dft-solvers/src/ref_davidson_complex.rs`: owned-row block Davidson
  with fixed-tree reductions; Hermitian EVD and Cholesky through
  `dft-numerics::linalg::ComplexOverlapFactor` (faer).
- `crates/dft-physics/src/distributed_complex.rs`: `DistributedComplexHamiltonian`
  (owned G rows, `DistributedPairFft`, per-k immutable kinetic/projectors,
  `band_density_sum`, scratch bounded by a band capacity).
- `crates/dft-physics/src/stress.rs`: radial derivatives of projector form
  factors and real-harmonic gradients (`harmonic`), reusable for the
  k-derivative of the nonlocal operator.
- `crates/dft-geometry/src/{kpoints,symmetry}.rs`: Monkhorst-Pack meshes,
  time-reversal reduction, `reduce_mesh` -> `KOrbitMap` (representative,
  operation, reciprocal shift, multiplicity).
- ADR 0041 / `mpi-m2-plan.md`: `dft-mpi` FFI crate, deterministic
  `tree_sum_f64`, G-row (stick) distribution, determinism contract.
- ADR 0047: the reference driver is the only production SCF path; ground
  states for this plan come from it (verify on), never from the oracle.
- DFT+U is not implemented (`"hubbard":"none"` in result records); it is
  MF10 in `master-plan-proposal.md`, after MF04 (occupations) and MF06 (NC
  coverage).

## 1. Goals and non-goals

Goals, in order:

1. Dielectric response of insulators: chi0(q, i omega), epsilon_GG'(q, i omega)
   with local fields, the macroscopic constant epsilon_M (epsilon_infinity in
   the RPA sense), as a validated building block.
2. RPA correlation energy from the adiabatic-connection fluctuation-dissipation
   theorem (ACFDT) on the imaginary axis, for insulators first:

   E_c^RPA = (1/N_q) sum_q  integral_0^inf (d omega / 2 pi)
             Tr[ ln(1 - chi0(q, i omega) v(q)) + chi0(q, i omega) v(q) ]

   with v_GG'(q) = 4 pi delta_GG' / |q + G|^2 and N_q the q mesh size.
3. Optional later stage: EXX + RPA total energy, E = E_EXX[psi_KS] + E_c^RPA,
   with E_EXX the Hartree-Fock energy expression evaluated on KS orbitals.
4. cRPA: partially screened interaction W_r and its matrix elements U, U', J
   (and the full tensor U_ijkl) for a correlated subspace, to feed DFT+U
   (MF10) and to export for model solvers.

Non-goals now: GW, BSE, TDDFT kernels beyond an optional ALDA check, RPA
forces and stress, PAW/USPP, SOC/noncollinear, real-frequency U(omega) on a
dense real axis, 2D Coulomb truncation for slabs, GPU. Metals are deferred
(section 2.7). No speed target versus QE applies to these stages (QE has no
matching production RPA path); cost is measured and ledgered, correctness
first.

## 2. Physics and numerics contracts

Conventions: Hartree atomic units, Omega the cell volume, N_k k points of the
full mesh, spin-unpolarized factor 2 written explicitly, orbitals normalized
per `docs/physics-conventions.md`. All matrices are indexed by G in a
response basis {G : |q + G|^2 / 2 <= E_cut^chi}, N_G its size.

### 2.1 chi0: sum over states (Adler-Wiser) versus Sternheimer

Adler-Wiser on the imaginary axis (insulator, occupations f in {0,1}):

  chi0_GG'(q, i w) = (2 / (N_k Omega)) sum_k sum_{n occ} sum_{m emp}
      rho_nm^k(q + G) rho_nm^k(q + G')^*  *  2 (e_nk - e_m,k+q) / ((e_nk - e_m,k+q)^2 + w^2)

with pair densities rho_nm^k(q + G) = <psi_nk | e^{-i(q+G).r} | psi_m,k+q>,
computed by FFT of psi_nk^*(r) psi_m,k+q(r). The frequency factor is real and
negative, so chi0(i w) is Hermitian negative semidefinite, and
X = v^{1/2} chi0 v^{1/2} has eigenvalues in (-inf, 0]. Then
Tr ln(1 - chi0 v) = ln det(1 - X) = 2 sum ln diag(L), L the Cholesky factor of
the positive definite 1 - X. Cost per (q, all w): N_k N_occ N_emp pair FFTs plus
a GEMM of N_G x (N_occ N_emp N_k) times its adjoint, frequency weights applied
as column scaling, i.e. O(N_G^2 N_occ N_emp N_k N_w).

Sternheimer / DFPT alternative (Nguyen and de Gironcoli 2009; Wilson, Gygi,
Galli 2008, PDEP): apply chi0(i w) to a trial potential through linear solves
(H - e_n +/- i w) dpsi_n = -P_c dV psi_n for occupied n only, and obtain the
leading eigenpairs of v^{1/2} chi0 v^{1/2} iteratively; E_c converges with a
few eigenpairs per electron. No empty states; cost O(N_eig N_occ N_PW log N_PW
x iterations) per frequency.

Comparison for BAANDI:

| Aspect | Sum over states | Sternheimer / PDEP |
| --- | --- | --- |
| New numerics | pair-density FFT + GEMM; many-empty-band solve (R1) | complex shifted linear solver (COCG/MINRES with P_c), iterative eigensolver on chi0 v |
| Reuse | NSCF, Davidson or dense EVD, DistributedPairFft, GEMM dispatch | H apply, Davidson framework, P_c from occupied block |
| Empty-band convergence | explicit, a known slow error (2.5) | none |
| cRPA | natural: constrained chi0 removes transitions (2.9) | awkward: the exclusion needs explicit subspace states anyway |
| Head/wings | velocity matrix elements, standard | needs separate d/dk Sternheimer |
| Oracles (Yambo, BerkeleyGW, ABINIT, RESPACK) | all sum over states: same truncation conventions, comparable per stage | only QE-type ACFDT; hard to match truncation |
| Scaling for large cells | quartic | cubic-quartic with smaller prefactor |

Recommendation: sum over states is the primary, validated path for R2 to R5,
because cRPA needs the empty states anyway, every oracle we can actually run
uses it, and its building blocks already exist. Sternheimer/PDEP is a later,
isolated and switchable module (generality rule: modules that matter for large
systems are built even when they do not pay off on small cases), validated
against the sum-over-states path on the R3 set.

### 2.2 Head and wings at q -> 0

For q -> 0 the pair density of an interband transition is linear in q:
rho_nm(q) ~ -i q . <n| v_op |m> / (e_m - e_n), with the velocity operator
v_op = -i nabla + i [V_NL, r], which in a k-dependent plane-wave basis is
d H_k / d k. The nonlocal term is mandatory for NC pseudopotentials (its
omission changes epsilon_M by several percent). It is computed from
d/dk of beta_lm(k + G), i.e. radial form-factor derivatives and real-harmonic
gradients, the same functions as the stress nonlocal term (reuse after review;
ownership of `dft-physics` applies). Head and wings:

  chi0_00(q -> 0) = q^2 * (q_hat . L . q_hat), chi0_G0 ~ q (q_hat . w_G)

stored as the 3x3 tensor L and vectors w_G; epsilon_M(q_hat) = 1 / [epsilon^{-1}]_00.
Primary: k.p as above. Cross-check (internal oracle, not default): finite
q = 1e-3 to 1e-2 of a reciprocal vector with an NSCF on the shifted mesh and
quadratic extrapolation; the two must agree to 0.2 percent on Si and LiF.
Anisotropic cells use the full tensor; the RPA q = 0 term averages over
directions (2.3).

### 2.3 Coulomb singularity

- Response (epsilon, E_c): the q = 0, G = 0 element enters only through the
  finite product v chi0 via the head tensor. For E_c the q = 0 term of the mesh
  sum is replaced by its average over the mini Brillouin-zone cell around
  Gamma: the head of v chi0 is a direction-only function, averaged over a
  fixed Lebedev set (26 or 50 directions, chosen by a convergence check, not by
  material). The q = 0 term carries an O(N_q^{-1}) error that must vanish with
  the mesh; the k-mesh convergence series (2.6) records it.
- Bare and screened Coulomb integrals of localized orbitals (cRPA, EXX): the
  integrable 1/q^2 divergence of v(q + G) at G = 0 is treated by an auxiliary
  function (Gygi-Baldereschi type) integrated analytically over the BZ and
  numerically for the remainder, or by the Spencer-Alavi spherical truncation
  for EXX. Both are implemented and must agree to 1 meV on the bare U of an
  isolated-band subspace (SrVO3 t2g) as the mesh grows. Choice recorded in
  every output.

### 2.4 Frequency grids

Default: Gauss-Legendre on t in (0, 1), w = w0 t / (1 - t), weight scaled by
the Jacobian, N_w points; w0 chosen per system by the rule w0 = mean of the
lowest direct transition and the largest transition included (a rule from the
computed spectrum, not a material constant). Convergence protocol: doubling
N_w from 8 until |Delta E_c| < 0.05 meV/atom on every set member, recorded.
Minimax grids (Kaltak, Klimes, Kresse 2014) are a later optimization; we would
compute our own grids by a Remez fit for the spectrum ratio, never copy
tabulated grids. cRPA static U uses w = 0 exactly; U(i w) on the same grid is
an output (2.10).

### 2.5 Cutoffs and convergence protocol

Parameters: wavefunction cutoff E_cut (from the ground state), response cutoff
E_cut^chi <= 4 E_cut (density limit), empty bands N_b. Following Harl and Kresse
(PRB 2008, 2010), the correlation energy converges as

  E_c(E_cut^chi) = E_c^inf + A (E_cut^chi)^{-3/2}

equivalently linear in 1/N_G (N_G grows as (E_cut^chi)^{3/2}), when the band
count follows the basis: N_b = N_G(E_cut^chi) (all bands below the response
cutoff). Protocol: at least four E_cut^chi values spanning a factor >= 2,
N_b tied to N_G, linear fit in (E_cut^chi)^{-3/2}, report E_c^inf, the fit
residual and the slope; the fit must be checked for linearity (the residual of
the last three points below 0.5 meV/atom), else more points. Energy
differences (lattice constants, atomization) are extrapolated per point, not
taken from a single cutoff. Note: the user brief wrote "1/N_G^{3/2}"; the
published law is 1/N_G, i.e. E_cut^{-3/2}. epsilon_M convergence: in N_b and
E_cut^chi separately, until relative change < 1e-3.

### 2.6 q-point sampling and symmetry

Gamma-centered MP k mesh only (so k + q stays on the mesh; shifted meshes
reject). q runs over the same mesh; q is reduced by the crystal point group
plus time reversal using `reduce_mesh` weights, because E_c(q) and the trace
are invariant. Stage 1: chi0(q) for an irreducible q sums over the FULL k mesh,
so the NSCF is run on the full mesh (no wavefunction rotation). Stage 2
(optimization, own ticket): NSCF on the IBZ, orbitals at other k obtained by
the `KOrbitMap` operation and reciprocal shift with the documented phase
convention, validated bit-for-bit-equivalent (1e-12 relative) against stage 1.
k convergence: meshes n, n+2, n+4 recorded with E_c per atom.

### 2.7 Metals and spin

Metals add the intraband (Drude) head at q -> 0 and fractional occupation
factors (f_n - f_m) in chi0; Fermi-Dirac exists in the reference SCF. Deferred
to a separate stage after R5 (decision D6); until then a metallic ground state
rejects in the response driver with an explicit message. Collinear spin:
chi0 = chi0_up + chi0_down (no factor 2), v unchanged, E_c formula unchanged;
enabled in R3 for fixed-moment insulators (O2 triplet in a box is an
available reference-path case). cRPA for magnetic systems (NiO AFM) uses
spin-resolved constrained chi0 and spin-resolved U (R5b).

### 2.8 cRPA correlated subspace

Options:

- (a) Atomic-orbital projection: the UPF PSWFC chi_lm (already used by the
  reference atomic start) projected on a band window and Loewdin
  orthonormalized ("projected Wannier", no gauge optimization). Cheap,
  deterministic, and the same object as the DFT+U ortho-atomic projector, so U
  and its consumer share one convention.
- (b) Maximally localized Wannier functions: projection as the initial gauge,
  then spread minimization. Implemented externally through the Wannier90 file
  interface (BAANDI writes `.win`, `.amn`, `.mmn`, `.eig`; Wannier90 3.x runs as
  an external program; BAANDI reads `_u.mat` / `_u_dis.mat`). No linking to
  Wannier90 (GPL) from the MIT code base; an own MLWF minimizer is a possible
  later module.
- (c) Disentanglement for entangled bands (Souza-Marzari-Vanderbilt): outer and
  frozen windows; through (b), same file interface.

Recommendation, staged: R4a = (a) for isolated and entangled sets (it is
always defined); R4b = (b)+(c) through Wannier90 files, which also gives an
interoperable export users expect. cRPA (R5) supports both; the DFT+U handoff
(R6) uses (a) unless the user selects otherwise, so the projector of U and of
the DFT+U functional is identical.

### 2.9 Constrained chi0

chi0 = chi0_d + chi0_r; chi0_d contains transitions inside the correlated
subspace, W_r = (1 - v chi0_r)^{-1} v.

- Band-based exclusion (Aryasetiawan et al. 2004): remove pairs (n, m) both in
  the correlated band window. Exact only when the subspace equals an isolated
  band set (SrVO3 t2g).
- Projector-based weighting (Sasioglu, Friedrich, Bluegel 2011; Kaltak thesis
  2015): with p_nk = sum_i |<w_ik | psi_nk>|^2, weight each transition of chi0 by
  (1 - p_nk p_m,k+q) so that only the subspace-subspace part is removed;
  reduces to band exclusion for isolated sets.

Rule: both are implemented; the default is the projector weighting (defined for
every case and identical to band exclusion on isolated sets, which a test
asserts); band exclusion stays selectable for comparison with published
numbers. Recorded in every output.

### 2.10 Static versus frequency-dependent U, bare versus screened

Static U = W_r(w = 0) is the DFT+U input. U(i w) on the R2 frequency grid is
written as a diagnostic (it is smooth and real on the imaginary axis); real
axis U(w) needs analytic continuation or a real-frequency chi0 and is out of
scope. Every output gives both the bare V_ijkl (v) and the screened U_ijkl
(W_r), with matrix elements

  U_ijkl(R) = (1 / (N_q Omega)) sum_q sum_GG' rho_ij(q + G)^* W_r,GG'(q) rho_kl(q + G') e^{iq.R}

rho_ij the Fourier transform of w_i^*(r) w_j(r) over the supercell; on-site
R = 0 first, nearest-neighbour V later. Derived quantities: Hubbard-Kanamori
U = average U_iiii, U' = average_{i!=j} U_iijj, J = average_{i!=j} U_ijji; for full
l shells Slater F0, F2, F4 by least squares with the residual reported
(J = (F2 + F4)/14 for d). Dudarev U_eff = U - J is a derived field, never the
only output.

### 2.11 Outputs (output-format policy)

Writers live in `dft-cli` (ADR 0017); the core returns typed arrays.

- epsilon: `epsilon.json` (metadata, units, settings, head tensor), `epsilon_M.csv`
  and `.dat` (direction, with/without local fields), `chi0_head.dat` vs i w.
- RPA: `rpa.json` (E_c per q and per w, weights, extrapolation fit),
  `rpa_extrapolation.csv`, SVG convergence figure.
- Wannier: Wannier90 `.win/.amn/.mmn/.eig` (public format docs), and
  `seedname_hr.dat` style hoppings from projection.
- cRPA: VASP-documented text layout `UIJKL`, `VIJKL`, `UIJJI`, `UIIJJ`
  (i j k l Re Im, eV) from the public VASP wiki description; RESPACK-style
  `Uij`/`Jij` matrices per R from its public manual; `crpa.json` (Ha and eV,
  subspace definition, exclusion scheme, Coulomb treatment, frequency, orbital
  ordering and real-harmonic convention); `crpa_U_iw.csv`. Each writer gets an
  independent reader test (Python reader per public spec), per policy.

## 3. Software design

### 3.1 Crates and modules

- New crate `dft-response` (numerical core, no I/O, no serde): `pair` (pair
  densities on the owned-row layout, batched FFT via `DistributedPairFft`),
  `chi0` (accumulation GEMM, frequency weights, head/wings), `velocity`
  (d H_k / dk, incl. nonlocal), `coulomb` (v, mini-BZ average, auxiliary
  function, truncation), `dielectric` (epsilon, inversion, epsilon_M), `acfdt`
  (log det, frequency/q sums, extrapolation fit), `subspace` (projection,
  Loewdin, Wannier gauge application), `crpa` (constrained chi0, W_r, matrix
  elements, Slater fit). Depends on dft-numerics, dft-physics, dft-basis,
  dft-geometry, dft-mpi, dft-solvers. Needs an `ownership.toml` entry and a
  crate AGENTS.md (orchestrator).
- `dft-engine::response` (orchestration): ground state via
  `refscf::api::solve` (verify on), NSCF bands via an extended `refscf::nscf`,
  then the response pipeline; checkpoint capture.
- `dft-solvers`: a many-band eigensolver option (3.2).
- `dft-cli`: subcommands, TOML parsing, writers, Wannier90 file I/O.
- `dft-mkl`: optional ScaLAPACK FFI (3.4), behind features.

### 3.2 Reuse of NSCF, reference Davidson, distributed FFT, symmetry

- NSCF: reuse `start_kpoints_from_density` and the fresh-residual certificate.
  Changes needed (R1): (i) band counts in the hundreds to thousands: the
  64-restart all-band loop and the strict-all-bands mode are not designed for
  that; add a buffer of extra bands not certified and discarded, and locking
  of converged bands; (ii) a dense path: for N_PW up to about 6000 build H_k
  explicitly (kinetic diagonal, local V(G - G') from the potential FFT, nonlocal
  B D B^dagger from the existing projector matrix) and diagonalize with faer
  Hermitian EVD; all eigenpairs, exact up to EVD roundoff. This is both the
  production path for small cells and the R1 oracle; (iii) Chebyshev-filtered
  subspace iteration as the large-N_b iterative path, using `apply_block`,
  validated against (ii). The existing Davidson stays for the occupied block.
- Distributed FFT: pair densities psi_nk^* psi_m,k+q on each rank's planes;
  no new transposes beyond `DistributedPairFft`'s batched API.
- Symmetry: `reduce_mesh` for q weights (2.6); wavefunction rotation later.
- GEMM: `dft-numerics` dispatch (faer default, MKL feature) for the chi0
  accumulation; this is the hot loop and follows `docs/performance/README.md`
  (no allocation, logging or I/O inside).

### 3.3 Memory

chi0 is N_G^2 complex: 16 MB at N_G = 1000, 1.6 GB at 10^4. Pair densities for
one k: N_occ N_b N_G x 16 B (Si 4 x 400 x 400: 10 MB; SrVO3 with semicore
about 20 x 600 x 1500: 290 MB, streamed in band blocks). All frequencies are
accumulated from one pass over the pair densities (N_w matrices in memory or
frequency blocks). E_c needs only one chi0(q, w) at a time after accumulation.

### 3.4 Distribution under MPI (behind `--features mpi`)

- Level 1 (R2, R3): transitions distributed over ranks (k, n-block, m-block in
  a fixed assignment), each rank accumulates a partial replicated chi0 for all
  w, combined by `tree_sum_f64` (fixed association, ADR 0041 determinism
  contract: bit-reproducible at fixed rank count, np=1 bit-identical to
  threads). Ground-state and NSCF use the existing G-row layout; pair
  densities are formed on G rows, so the GEMM is G-row-local and the chi0 row
  block is owned; the replicated matrix is assembled by allgatherv.
- Level 2 (only on measured need: N_G^2 N_w exceeds about a quarter of rank
  memory): 2D block-cyclic chi0 with ScaLAPACK (pzpotrf for log det, pzgesv/
  pzpotri for epsilon^{-1}) through MKL ScaLAPACK + BLACS for Open MPI on
  dev08, as an FFI module in `dft-mkl` with documented preconditions.
  Requires ADR 0003 checks and a `docs/dependencies.md` entry. ELPA and an own
  distributed Cholesky are alternatives (decision D7).
- (q, w) parallelism is not chosen: all w share one pass over pair densities,
  and q counts are small after symmetry.

### 3.5 Checkpointing

Expensive units: NSCF bands per k (minutes to hours), chi0(q) for all w,
constrained chi0 and W_r. Checkpoint per completed unit in the BAANDICP
container (rank-indexed, checksummed, manifest-last), new payload kinds
`nscf_bands_v1` (eigenvalues, residuals, owned rows, k, input/UPF/backend
fingerprint) and `chi0_q_v1`. A response run resumes by skipping completed
units whose fingerprints match; mismatch rejects. Pair densities are
recomputed, never stored. Payload design needs an ADR extension of
`docs/checkpoint-format.md` (R1 ticket).

### 3.6 CLI and input

Subcommands: `baandi dielectric input.toml`, `baandi rpa input.toml`,
`baandi wannier-export input.toml`, `baandi crpa input.toml`, each running or
reading the ground state like `baandi bands`, `--output-dir`, `--resume DIR`.
Input (key names follow the explicit-unit convention of P7U):

```toml
[response]
ecut_chi_ha = 5.0            # response basis cutoff
bands = "auto"               # auto = N_G(ecut_chi); or an integer
band_solver = "dense"        # dense | chebyshev
k_mesh = [8, 8, 8]           # Gamma-centered, full-mesh NSCF
head = "kp"                  # kp | finite-q (check only)
local_fields = true
[response.frequency]
grid = "gauss-legendre"      # gauss-legendre | minimax (later)
points = 16
omega0_ha = "auto"
[rpa]
ecut_chi_ha_series = [3.0, 4.0, 5.0, 6.0]
q0_directions = 26
[subspace]
method = "atomic-projection" # atomic-projection | wannier90
species = "V"
l = 2
orbitals = ["dxy", "dyz", "dxz"]
window_ha = [-0.2, 0.15]     # relative to the Fermi level / VBM
[crpa]
exclusion = "projector"      # projector | band
frequencies = "static"       # static | grid
coulomb_q0 = "auxiliary"     # auxiliary | spherical-truncation
```

## 4. Staged tickets

Common rules for every stage: oracle numbers count only from real runs with
command, code version, input file, pseudopotential SHA-256 and machine stored
under `tests/reference/response/<stage>/`; no UPF committed; QE-family runs on
dev08 through the master session with start/done/failed/idle signal files;
tests over about 10 s debug are `#[ignore]` with their release command; Mac
timings are observations only. Generality set (rule of 2026-11-03), used from
R2 on: Si and C (diamond, cubic), LiF (ionic, large gap), wurtzite AlN
(hexagonal), a skewed/strained Si cell, and a molecule in a box (N2 or H2O,
elongated cell) for R3. Defaults (frequency count, buffer bands, direction
sets) are accepted only if they hold on the whole set. Oracle software is not
installed on dev08 without decision D3.

### R1: many-empty-band solve

Scope: NSCF extension (buffer bands, locking, dense path, Chebyshev path,
`nscf_bands_v1` checkpoint), result with all bands to N_b.
Oracles: (i) internal dense EVD of the explicit H_k in the complete basis
(N_b = N_PW); completeness: sum_n |psi_n><psi_n| = 1 to 1e-10; (ii) QE 7.6
`pw.x` nscf with the same UPF, cutoff, density-converged SCF, `nbnd` large,
`diago_full_acc = .true.`.
Cases: Si 2-atom (4x4x4 and 8x8x8), LiF, AlN, skewed Si.
Gates: iterative vs dense eigenvalues <= 1e-8 Ha for all certified bands,
residual certificate <= 1e-6 Ha; vs QE eigenvalues <= 2e-5 Ha for all bands
up to N_b (ground-state offset recorded separately); rank counts 1/2/4/8 within
1e-10 Ha; np=1 bit-identical to threads.
dev08 cost: Si 8x8x8, 512 k, N_PW about 750, 400 bands, dense: under 5 min on
8 cores; QE side similar.

### R2: chi0 and epsilon for semiconductors

Scope: pair densities, chi0 with frequency grid, head/wings with nonlocal
velocity, epsilon with and without local fields, epsilon_M tensor, MPI level 1.
Internal oracles: (a) Lindhard function of the free-electron gas on the
imaginary axis, analytic, from synthetic plane-wave states (no UPF); (b) static
chi0 columns by finite differences of the NSCF occupied density response to
a small periodic potential at q != 0 (independent code path), against the
sum over ALL states of the dense path: 1e-6 relative; (c) Hermiticity and
negative semidefiniteness on the imaginary axis; (d) the high-frequency limit
chi0_GG' -> -(q+G).(q+G') n(G - G') / w^2 as an empty-band completeness check;
(e) k.p head vs finite-q extrapolation within 0.2 percent.
External oracles (same UPF, cutoff, k mesh, N_b, E_cut^chi): Yambo (RPA, with
and without local fields, from a QE `pw.x` + `p2y` ground state), BerkeleyGW
`epsilon.x` (from `pw2bgw`), and QE `ph.x` (DFPT epsilon_infinity, which
includes the LDA xc kernel and is matched only to an optional BAANDI ALDA
kernel check, never to plain RPA).
Cases: Si, C, LiF, AlN (tensor), skewed Si.
Gates: epsilon_M vs Yambo and BerkeleyGW at matched settings: 0.5 percent
relative (no local fields and with local fields separately); chi0 head tensor
0.5 percent; DFPT-matched ALDA value vs ph.x: 1 percent. Experiment (Si 11.9,
C 5.7, LiF about 1.9) is reported with the RPA/LDA overestimate stated, never
a gate.
dev08 cost: Si 8x8x8, N_G about 400, 400 bands, one q (q -> 0): about 2e11
complex MAC, seconds to minutes; the full set with convergence series about
2 to 4 h on 8 ranks including oracle runs.

### R3: RPA correlation energy for small systems

Scope: ACFDT log det, q0 mini-BZ average, q symmetry weights, extrapolation in
(E_cut^chi)^{-3/2}, checkpoint per q, spin channels for fixed-moment
insulators. Optional R3b: EXX on KS orbitals and the EXX + RPA total energy.
Oracles: ABINIT (`gwrpacorr`, NC pseudopotential in UPF or converted psp8 with
both hashes recorded) at matched settings; QE ACFDT if a working
implementation is available in the dev08 QE 7.6 build (availability not
verified; not assumed); published values (Harl and Kresse 2010; Nguyen and de
Gironcoli 2009) only for energy differences, with their PAW/NC settings stated,
as sanity bands. Internal: frequency and direction convergence; E_c(N_q)
trend; Sternheimer path later.
Cases: Si, C, LiF, Ne or Ar fcc, N2 in a box (atomization needs N atom: spin),
skewed Si; lattice constant of Si from E_EXX+RPA in R3b.
Gates: |E_c(BAANDI) - E_c(ABINIT)| <= 1 meV/atom at matched cutoff, bands,
k, frequency grid; extrapolated values 2 meV/atom; frequency-grid error
<= 0.05 meV/atom; rank-count independence within 1e-9 Ha.
dev08 cost: Si 8x8x8, 29 irreducible q, N_G about 400, 16 w: about 6e13
complex MAC, 15 to 30 min on 8 cores per cutoff; a four-cutoff series and
two k meshes about 4 to 8 h per material; ABINIT comparable.

### R4: projection and Wannier subspace

R4a: atomic projection with Loewdin, band window, hoppings in real space.
R4b: Wannier90 file export and import of gauge matrices.
Oracles: Wannier90 3.x run on BAANDI exports versus Wannier90 on QE
`pw2wannier90` exports with the same UPF and windows; QE `projwfc.x` projected
populations for R4a.
Cases: Si sp3 (entangled, disentanglement), SrVO3 t2g (isolated), diamond,
AlN.
Gates: orthonormality 1e-12; `.amn/.mmn` overlaps vs QE export after a gauge
check: |M| elementwise within 1e-5; final spreads within 1 percent; hr.dat
hoppings within 1 meV; Wannier-interpolated bands within 1 meV in the frozen
window; projection populations within 1e-3.
dev08 cost: minutes per case.

### R5: cRPA

Scope: constrained chi0 (projector default, band option), W_r, bare and
screened U_ijkl on site (R = 0) and nearest neighbour, Kanamori and Slater
outputs, U(i w), Coulomb q0 treatment, writers. R5b: collinear spin (NiO AFM).
Oracles: RESPACK (cRPA from QE ground state, MLWF, same UPF) for SrVO3 t2g;
the same subspace in BAANDI through the Wannier90 gauge of R4b for a
like-for-like comparison, plus projection-only for the DFT+U convention.
Internal: bare V of two Gaussian charge densities against the analytic erf
expression; band exclusion equals projector weighting on isolated sets to
1e-10; U -> V when screening is switched off.
Cases: SrVO3 t2g and full d (entangled), a 3d metal-free test (Si sp3 cRPA as
generality fixture), NiO (R5b), published SrVO3 t2g values (Aryasetiawan et al.
2006; Vaugier, Jiang, Biermann 2012) as sanity only.
Gates vs RESPACK at matched Wannier functions (spread within 1 percent): bare V
within 1 percent, static U within 3 percent, J within 5 percent; the two q0
treatments within 1 meV on bare U at the largest mesh.
dev08 cost (estimate, to be refined at ticket time): SrVO3 with V semicore,
k 4x4x4, N_G about 600, 300 bands, static only: about 1e12 complex MAC per
irreducible q, under 1 h on 8 ranks; a k-convergence series to 6x6x6 several
hours; RESPACK runs similar.

### R6: DFT+U consumption

Depends on MF10 (DFT+U functional) being accepted; R6 does not implement DFT+U.
Scope: `crpa.json` -> DFT+U input with a projector fingerprint (UPF hash,
orbital set, orthonormalization, window) that MF10 checks and rejects on
mismatch; Dudarev and Liechtenstein forms from the same record; documented
round trip. Optional: self-consistent loop (cRPA on the DFT+U ground state).
Oracle: QE 7.6 `pw.x` DFT+U (ortho-atomic projector, same U and J values) on
SrVO3 or NiO, as required by MF10.
Gates: imported U equals exported U bit-for-bit; DFT+U energies and occupation
matrices vs QE per MF10 gates; mismatched projector rejects.
dev08 cost: as MF10.

## 5. Risks

- Cost and memory: chi0 scales as N_G^2 N_occ N_b N_k N_w, about N_atoms^4; a
  20-atom cell is near the limit of the level-1 distribution. Mitigation:
  pair-density streaming, frequency blocking, IBZ wavefunction rotation, level-2
  ScaLAPACK, later Sternheimer/PDEP or space-time methods, each as its own
  ticket with ledger rows.
- Empty-band and cutoff convergence: E_c converges slowly; a wrong protocol
  gives errors of tens of meV/atom that look converged at fixed N_b. Mitigation:
  N_b tied to N_G and mandatory extrapolation with linearity check; never a
  single-cutoff headline.
- Many-band solve: high bands of iterative solvers converge poorly; NC
  projectors can produce ghost states high in the spectrum. Mitigation: dense
  oracle, buffer bands, ghost scan with the UPF qualification tools.
- NC pseudopotentials: absolute E_c lacks core-valence correlation and depends
  on the dataset; only energy differences are compared with literature.
  Nonlocal commutator must be in the velocity (several percent in epsilon_M).
  3d systems need semicore and high cutoffs (cost).
- Head/wing and q0 errors dominate small k meshes; mitigation: k.p vs finite
  q, direction averaging, mesh series.
- Metals: intraband term, Fermi-surface sampling and smearing interplay;
  deferred by design.
- Validation gaps: no installed RPA-energy oracle on dev08 today; QE ACFDT
  availability unknown; BerkeleyGW needs registration; Yambo/BerkeleyGW/RESPACK
  start from a QE ground state, so cross-code differences include the
  ground-state offset (measured separately, about 1e-6 Ha/atom so far for LDA);
  cRPA comparisons depend on matching Wannier functions, which is itself a
  source of percent-level differences.
- Licensing: Wannier90, RESPACK, Yambo, ABINIT are GPL family: run as external
  programs through files only, never linked or copied. Published grids/tables
  (minimax) are not copied.
- Determinism: chi0 accumulation order must be fixed for reproducibility;
  GEMM backends with threading must use the existing deterministic dispatch.

## 6. Order and dependencies

R1 -> R2 -> R3 (R3b optional) ; R1 -> R4a -> R4b ; R2 + R4 -> R5 -> R5b ;
R5 + MF10 -> R6. R4a may run in parallel with R2/R3. Each stage: ticket with
OWNED_FILES, design ADR where a contract is new (response basis and head
conventions in R2; checkpoint payloads in R1; cRPA conventions in R5), Codex
implementation, `scripts/check.sh`, dev08 gates, README validation block
through `scripts/make_readme_figures.py` extension, wiki and ledger updates.

## 7. Open decisions for the user

D1 chi0 method. Options: (A) sum over states primary, Sternheimer/PDEP later
as an isolated module; (B) Sternheimer first; (C) both at once.
Recommendation: A.

D2 Wannier functions. Options: (A) atomic projection in BAANDI plus MLWF and
disentanglement through Wannier90 files (external program); (B) own MLWF
minimizer in Rust; (C) link the Wannier90 library. Recommendation: A (C is
excluded by the MIT/GPL boundary; B only on measured need).

D3 Oracle software on dev08. Options: (A) build Yambo, ABINIT, Wannier90 and
RESPACK in user space (no system change) through the master session, and check
`ph.x`, `pw2wannier90.x`, `projwfc.x` and any ACFDT tool in the QE 7.6 build;
BerkeleyGW only if registration is acceptable; (B) published values only.
Recommendation: A (B cannot meet the real-run oracle rule).

D4 cRPA exclusion default. Options: (A) projector weighting by default, band
exclusion selectable; (B) band exclusion by default. Recommendation: A.

D5 Frequency grid. Options: (A) mapped Gauss-Legendre now, own minimax later;
(B) minimax from the start. Recommendation: A.

D6 Metals. Options: (A) defer to a stage after R5; (B) include intraband terms
in R2/R3. Recommendation: A.

D7 Distributed dense algebra. Options: (A) replicated chi0 with
transition-parallel tree reduction, MKL ScaLAPACK only on measured need;
(B) ScaLAPACK from the start; (C) ELPA or own distributed Cholesky.
Recommendation: A.

D8 EXX + RPA total energy. Options: (A) optional R3b after R3; (B) skip for now.
Recommendation: A.

D9 cRPA ground state and DFT+U handoff. Options: (A) cRPA from the LDA/PBE
ground state, U exported with the ortho-atomic projector fingerprint, R6 after
MF10; (B) self-consistent cRPA with DFT+U from the start. Recommendation: A.

D10 Performance. Options: (A) no speed target for R1 to R6, cost measured and
ledgered, the 0.90x QE target unaffected; (B) a speed target against an oracle
code. Recommendation: A (no same-model, same-algorithm QE counterpart exists).

## 8. User decisions (2026-11-26)

D1-D10 approved as recommended (option A for each): sum over states first; Wannier90 only as an external program through files; oracle codes (Yambo, ABINIT, Wannier90, RESPACK) built in user space on dev08 through the master session; projector-weighted exclusion by default; mapped Gauss-Legendre grids; metals deferred after R5; replicated chi0 with tree reduction, ScaLAPACK on measured need; EXX+RPA as optional R3b; cRPA from the LDA/PBE ground state with R6 after MF10; no speed target for R1-R6.
