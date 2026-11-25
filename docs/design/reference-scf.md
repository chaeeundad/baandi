# Reference SCF path (issue #25): design

Status: design draft by the orchestrator, 2026-11-23. User decisions of the same day: all lane/A-B work stopped; a NEW SCF driver replaces the legacy one; numerical kernels are reused only where they are good, and bold rewrites are allowed wherever they promise speed (the old code stays in git history). Implementation goes through Codex tickets; QE conventions reach implementers only through the Claude note `docs/research/qe-study/qe-density-mixing-spec.md` and the published literature.

## Why a new driver
The dev08 decomposition (docs/performance/requests/si64-gap-attribution-2026-11-20.md, architecture brief addendum) shows the ratio to QE as per-map cost x map count x non-SCF overhead: si32 1.51 x 0.83 x 1.66, si64 1.56 x 1.00 x 1.49, si128 1.35 x 1.93 x 1.31 (A/D 3.586 measured), si216 1.21 x 1.50 x 1.47. The per-map cost is closing with size; the map count (si128: 29 maps vs QE 15, maps 3-13 diverge with the standard error growing 272x while no control reacts) and the overhead outside the SCF loop (coarse ladder 18-24 % of wall at every size, verification 9-11 %) are the levers. The legacy driver (`ScfSession::advance`, about 1100 lines, plus the cutoff ladder and verification scheduling, more than 130 switches) is where both live. A new driver that follows a published, QE-equivalent fixed-point scheme lets us (1) compare every iteration with QE term by term and (2) add BAANDI's own accelerations back one at a time, each with a measured gain.

## Scope of the new path
Real-Gamma, integer occupations, spin-unpolarized, norm-conserving, LDA/PBE: the cases of the benchmark suite (Si8...Si216, skew, slab, elongated, MgO Gamma, synthetic PPs). Everything else (k points, spin, smearing, DFT+U, forces/stress outputs) keeps using the legacy driver until the new one covers it; the CLI chooses by profile.

## Architecture
A new module `crates/dft-engine/src/refscf/` (no session/ladder control flow or legacy orchestration switches; the pure verification assessment is shared with `session.rs`):

1. `profile`: default-off execution profile `BAANDI_REFERENCE_SCF=on` (input `[scf] driver = "reference"`), resolved once, recorded in results.json (`scf.driver`). The CLI branch sits where `run.rs` chooses between `run_ladder` and `run_outcome`; output schema stays the same (`ScfOutcome`), plus a `reference_scf` block per iteration.
2. `density`: the mixed quantity in reciprocal space exactly as defined in the note (half-sphere storage for real Gamma, G=0 handling, which G shells are mixed and what happens to the rest, normalization, units), distributed by G vectors under MPI (owned G slabs from the existing grid plan; no full-grid replication).
3. `metric`: the Hartree-like inner product used for the residual norm (the SCF accuracy estimate) and for the history overlap matrix, with its exact weight in |G|^2, the G=0 rule and the half-sphere weights; one reduction per inner product, batched.
4. `broyden`: the modified Broyden/Anderson update (Johnson, PRB 38, 12807 (1988)) with the startup step, beta, history length, insertion/eviction order and the linear solve, as specified in the note; deterministic reductions (fixed order) so results are rank-count independent within rounding.
5. `driver`: `run(StartInput, Profile)` / `run_with(StartInput, Profile, Collective)` consume parsed input directly, build `refstart` once, and run the owned-row map with `ref_davidson`, the residual metric and patch-1 Broyden. No `ScfProblem`, session or ladder is constructed. After an estimator stop, a fresh map at the delivered output is assessed by the shared fixed-point function extracted from `session.rs`; all attempts are charged. See Driver below.
6. `record`: per iteration: dr2, threshold, Davidson iterations and H columns, beta and history occupancy, step norm, energy terms, wall per phase; plus run totals. This is the comparison surface with QE (estimated scf accuracy, ethr, avg # of iterations).

## Patches (issue #25)
- Patch 1 (runnable density update): modules 1-4 and the driver skeleton; a reproducible test that calls the production update routine on (a) synthetic arrays with closed-form answers and (b) QE reference arrays (input density, output density, history) generated on dev08 from a recorded QE 7.6 Si128 input; compare the next density, total charge and low-G components in common units against a documented tolerance; report the first failing term if any. Deliver: commit, CLI command, comparison table.
- Patch 2 (common start convergence): same atomic start at the final cutoff for both codes; compare output densities first (a mismatch points at potential, occupations or the inner solve), then next inputs (a mismatch points at mixing); show whether Si128's early excursion disappears and the map count at equal accuracy (QE's 15 is a comparator, never a stopping rule). Record the reference-style stop and BAANDI's certificate separately, both charged to the delivered wall.
- Patch 3 (speed): dev08 pure MPI8 matched runs, Si128 with Si64 and Si216 controls, A/A floor first, crossed A/B, effect >= 3x floor and |mean|/SE >= 3, uninstrumented external wall including startup, certificate and output; skewed/slab and non-Si PP coverage; strict 1e-10 / standard 5e-9 Ha before any default change. Then add BAANDI accelerations back one at a time (kept orbitals first, then the ladder only if it pays at equal accuracy), and rewrite kernels where the per-iteration profile shows waste (per-iteration non-column cost was 99 ms vs QE 69 ms at si64; Gram/rotate/update 2.53 s vs QE 0.78 s of non-eigen algebra).

## Rules
No QE source to implementers; equations and conventions from the note and the literature only, the note's INFERRED items verified against QE reference data. Same parallel model on both sides. No constant tuned to silicon. The legacy driver stays the default until the new path passes the gates on the generality set; then the old driver and its switches are removed in one step.

## Diagnostics level

I26-2-DIAGLEVEL: the reference driver resolves the existing run-owned
`[profile] davidson_stages = true|false` setting once at entry and passes it
explicitly to every inner solve, including bounded re-entry, first-step redo
and every independent verification attempt. `BAANDI_DAVIDSON_STAGES=on|off`
overrides the input, with the existing registry default `off`. In-process
callers use `RunConfig::enter`; no process-cached flag selects this level.
This changes diagnostic collection only. The reference driver remains opt-in;
all numerical defaults and the strict 1e-10 Ha / standard 5e-9 Ha classes stay
unchanged. Verification is independent of diagnostics and remains on by default.

- **Lightweight (default):** no per-band Ritz-history allocation and no inner
  Davidson stage-clock reads. `ritz_changes_ha` is an empty array and
  `inner_clocks` is null (unmeasured, not zero physical cost). Per-call iteration,
  active-width, H-column, reduced-order, restart, repair and unconverged-root
  counters remain available. Solver work/communication accounting remains on.
- **Full:** every correction has one Ritz change in Ha per band, including
  empty bands and unconverged/exhausted calls. Each call also reports inclusive
  solve time and eight exclusive stage times (seconds, local output rank).
  Ownership of the Ritz matrix moves into the driver record without a clone.
  `convergence.work.davidson_stages` aggregates all SCF and verification calls;
  MPI uses the existing rank maximum/minimum convention.

Coarse phase clocks, all scientific history, density, eigenvalues, energies,
verification attempts (including failures and recovery), run provenance and
`run.completed` remain present in both modes. The engine total includes all
attempts and recovery; it is not delivered wall. External timing must extend
through required files, console output, final JSON write, and process/MPI exit.

Reader contract: `reference_scf.diagnostics.schema` is
`reference-scf-diagnostics-1`, with explicit `level`. All existing iteration and
call keys, array types and locations are preserved. In full mode both historic
paths (`reference_scf.iterations[]` and
`convergence.history[].reference_scf`) still contain identical complete rows;
`convergence.last_iteration` retains its usual copy. The duplication is retained
for compatibility and is included in byte counts and external timings. An empty
Ritz array in lightweight mode means **not collected**, never zero iterations
or successful convergence. Consumers requiring per-band evidence must request
full mode; absence of the new metadata identifies the historical always-full
writer. `inner_clocks` and the diagnostics descriptor are additive fields.
No existing trace shape is removed or relocated. `refscf-check` and
`refscf-map` read QE density/stdout records, not this JSON; their four-record
density and `refscf-map-1` contracts are unaffected and covered by their tests.

Reproducer, byte definitions, Mac evidence and the pending dev08 A/B protocol:
[diaglevel-dev08-request.md](../performance/lanes/RSCF/diaglevel-dev08-request.md).

## Start

R-SCF-START (2026-11-23), explicit library entry
`dft_engine::refstart::start(StartInput, &impl Collective)`. Derived from the
initialization principles in `research/qe-study/qe-scf-cycle.md` (own-words
note); no QE source was read. No switch/default, legacy driver, CLI, or
`src/refscf/` change. The caller supplies parsed cell/atoms/NC potentials and
options, identically on every rank, and the communicator explicitly. Passing
an already constructed `ScfProblem` would pay the costs this entry removes.

### Initialization ownership audit

Reconciles `performance/lanes/INIT/diagnosis.md`,
`performance/lanes/LAUNCHIO/diagnosis.md`, and
`performance/lanes/OWNED/audit.md` against base `d97d864`.
Historical timing/RSS figures in those documents and the ticket are relayed
dev08 evidence, not measurements of this implementation. In particular,
INIT's atomic span/Ritz happens after problem construction; its nested times
must not be added to the inclusive init total. The LAUNCHIO startup and output
buckets are outside this library entry.

Notation: G full wavefunction rows, D density-sphere rows, F grid points,
L_r owned ORC rows, H_r owned representatives, F_r owned real plane points,
S_r owned full-grid reciprocal coefficients, Q projector columns, J atomic
columns, B bands, K=max(J,B), C operator capacity. Bytes use f64=8,
complex=16, usize=8 on supported 64-bit targets. Rank maxima use actual
partitions, not ceil(G/P); uneven and empty owners are tested.

| Stage | Legacy payload per rank and frequency | Reference requirement / omission |
| --- | --- | --- |
| Process, MPI, binary, input/UPF parsing | Once/process, O(input and radial mesh bytes); launcher/binary overhead is outside numerical init | Retain host responsibilities; no second parser or process launcher. SHA cache already exists; no new credit. |
| Problem/basis/grid preflight | Per problem (coarse and final), repeated G/D lists and grid maps in serial state, force identity and metric | One final G and D construction; G grid retained, D metadata dropped after fields. No coarse state, strain/force model, Coulomb stop metric or trial mixer. Replicated O(G+D+F) layout metadata remains. |
| Projectors | Already sparse complex 16 H_r Q, then real 8 L_r Q; coupling 8 sum(q_atom^2), per k/problem | Reuse `new_owned_with_rule`; construct real operator directly, then drop sparse complex source. No full G*Q table or serial Hamiltonian/FFT. Previously adopted sparse assembly is not a new saving. |
| Ionic phases/radials | Density coefficients 16D gathered, scattered 16F, real ionic 8F; per problem | Evaluate only owned density-sphere coefficients, shell cache per species; public local radial transform and structure-factor kernels. One distributed inverse; retain 8F_r ionic values. No full-valued ionic arrays or allgather. |
| Ewald | Candidate blocks, local accepted terms, P chained three-double broadcasts; once/problem, cache can avoid repeat | Same plan, scheme, candidate partition and exact ordered accumulation; run before operator/field allocations. Local term payload 8*(accepted real+reciprocal); no new scheme or reduction rounding. Rank one uses streaming serial evaluator. |
| Atomic density | Owned radial evaluation followed by full coefficient gather, 16F scatter/inverse, 8F clipped normalized density; per problem; starting_density clones 8F | Own coefficients, inverse to F_r planes, same nonnegative clipping, scalar charge reduction and rescaling. Forward once after clipping; retain 16S_r coefficients. No full field gather. |
| First potential | Serial Hartree FFT pair / XC and full potential replicas; first cold map | Reuse the density forward transform; 4pi/G^2 on owned full-grid coefficients, inverse to planes, pointwise/library LDA and ionic+(Hartree+XC). G=0 Hartree exactly zero. |
| Serial band density | 8F(B+1) idle contribution/total grids plus backend FFT workspace per problem | Entirely absent. At the OWNED audit's Si128 F=72*72*144, B=260 this alone is 1,558,683,648 logical bytes/rank; this is a shape calculation, not RSS. |
| Atomic Q, HQ, Ritz | Complex Q/HQ 16GK each, rotated 16GB; full gather remains after split CholQR2; once/cold start | Direct real ORC rows, Q 8L_r K, H applied in B-column chunks, Gram reductions K*width and local rotation. No full Q/HQ, no complex/ORC full input bridge. |
| M2_OWNED_INPUT | Merged, default off: removes 8GB full ORC temporary at eligible solve entries, still consumes full complex input | New start already returns 8L_r B; does not consult or alter this switch. Avoided complex source and atomic gather are additional, distinct savings. |
| Distributed H/plans | Local kinetic 8L_r, projectors 8L_r Q; image 8L_r C; hbuf 16ceil(C/2)*local_fft_rows; mirror send/recv; Gram scratch | Reuse existing kernel and deterministic layout; capacity B, atomic applications chunked at B. Consumer must also respect that capacity. Construct once. |
| Distributed H potential | 8F full-indexed buffer, unowned entries poisoned; owned plane writes | Still 8F: current FFT/apply/kept-potential API uses global indexing. Physical compacting requires a separately owned kernel change. All *meaningful values* are owned; do not claim all storage is distributed. |
| Diagnostics/delivery | Force model, mixer probes, snapshots, full output density/orbitals; once/problem or delivery | No diagnostics/materialization in `start`; optional `StartFields::gather` is an explicit test/delivery boundary only. No CLI output or outside-process speed claim. |

### Numerical and API contract

`StartState` returns the basis/grid/half-sphere, distributed operator, reusable
owned field transform, ionic/density/first-potential planes, owned forward
normalized density coefficients, owned orthonormal Ritz block, Ritz energies,
electron count, Ewald terms and per-rank payload counters. Plane offsets and
reciprocal indices are explicit; the density full-grid stick layout and the
wavefunction operator plane axis may differ. Existing `PlaneMap` moves the
potential between them without a full-grid gather. Operator epoch is zero.

Density samples use the existing PP_RHOATOM integral (including its existing
quadrature rule), negative ringing clipped at zero, then scaled to N_el.
The scalar charge sum is a fixed-tree reduction, reproducible for a fixed
rank count, rounding class across rank counts and against the serial sum.
`density_g[G=0]` is set to exactly N_el/Omega after the final forward transform;
other coefficients include the high-frequency components introduced by clipping.
Real density and its forward transform agree up to this zero-mode rounding.
Hartree uses the existing full-box signed-frequency convention (negative even
Nyquist), not a density-sphere truncation. As in the legacy Hartree evaluator,
the real part is used on skew even grids; a Nyquist imaginary component is not
silently used as a physical potential.

Atomic ordering is species, radial orbital, atom, m. Each owned coefficient
uses the existing radial transform, real spherical harmonic and ionic phase,
then the existing sqrt(2)*Re/Im ORC convention. Guarded CholQR2 reuses the legacy
column normalization, Cholesky/inverse and pivot guards through an additive
owned-form helper; no legacy calls change. Bad conditioning takes distributed
two-pass MGS with the unchanged 1e-12 relative dependence rule. Missing bands
use the legacy kinetic-ordered Hermitian cosine/sine completion. There is **no
random perturbation** and no PRNG consumption; the seed cannot affect this
start. Rank-deficient atomic data fail explicitly, not via a random repair.
The initial H/Ritz rotation is included before returning (and before the first
Davidson). Eigenvector phases and rotations inside degenerate eigenspaces are
not unique; compare orthogonality, Ritz energies and subspaces, not raw bits.

Current implemented envelope: exact Gamma, positive even integral valence,
fixed unpolarized occupations, NC atomic charge and l<=3 orbitals, LDA-PZ or
library LDA-PW with modified PW92. PBE/distributed gradients, NLCC,
frozen-strain and missing atomic-data inputs reject explicitly; their legacy
fallback belongs to the caller. These are remaining extensions of the broader
reference-driver design, not silently changed physics. Parsed input agreement
and communicator lifetime are caller contracts; mid-run remapping is not an API.
Local radial/projector/XC failures are agreed before later collectives.

### Memory-counter interpretation

`StartMemory` records exact logical lengths for output orbitals, sparse/real
projectors, retained meaningful fields plus the full-indexed operator potential,
and transient Ewald terms. `principal_peak_bound_bytes` and the corresponding
legacy field are **conservative named-array envelopes**, not sampled allocation
high-water marks or RSS. The reference envelope is the max of Ewald terms and

`8*(3F_r+F)+16S_r + 8L_r*(3K+B) + 64K^2 + 8L_r Q + 16H_r Q + 32S_r+16F_r`.

It covers overlapping raw/Q/update, chunked H and selected orbitals, reduced
matrices/factors, and field temporaries. The legacy comparator counts
`8F(B+1) + 16G*(2K+B) + 8GB + 64K^2 + 64F + 16H_r Q + 8L_r Q`, at the same
**final** cutoff (also maxed with the identical local Ewald term payload),
without claiming a coarse/final simultaneous RSS. The literal
buffer counts/phase lifetimes are auditable in the constructors; these envelopes
are not a byte-for-byte simulator of every legacy switch combination.

Excluded from both envelopes: basis/grid/FFT ownership metadata, radial caches,
FFT and GEMM workspaces, operator mirror/FFT scratch, library XC/EVD storage,
allocator capacity/alignment, thread stacks and MPI. In particular real GEMM
panel storage scales as 8*ceil(L_r/256)*K^2 plus its packed small matrix; it is
not zero. The request below requires all-rank RSS to establish the actual peak.
No allocator/RSS saving or startup speed is asserted from these counters alone.

### Validation and integration

`crates/dft-engine/tests/refstart.rs` is the independent legacy-comparison
surface. Fast tests cover local/nonlocal synthetic PPs, skew cells, kinetic
completion, rank-deficient rejection, and empty field/orbital owners at
1/2/3/8 simulated ranks. Release `silicon_matrix` compares displaced Si8/16/32
and skew Si8 at the final 15 Ha cutoff; all arms use the same deterministic
1e-4 fractional displacement to remove accidental Ritz degeneracy. It checks
charge, exact zero mode, field agreement, orthogonality, selected subspace and
Ritz energies against the legacy start. `legacy_driver_from_owned_start`
imports a gathered diagnostic copy through the unchanged warm-checkpoint API
and runs the real legacy standard driver, including verification. That gather
is test-only and is not part of `start`.

Commands and actual numerical results are recorded in the R-SCF-START ledger
entry and `crates/dft-engine/tests/refstart-evidence.json` (126 rank records).
The final scoped gate passes; the full gate retains three existing CLI cache
failures caused by requiring scratch inside the repository. No failing test
was weakened. Dev08 protocol: [start-dev08-request.md](../performance/lanes/RSCF/start-dev08-request.md).
R-SCF-1 must consume this state directly; rebuilding `ScfProblem` would defeat
the setup/memory purpose. No CLI profile integration is made by this ticket.

## Implemented conventions (R-SCF-1, patch 1)

The independent implementation is `dft_engine::refscf`. It is a child module
of `scf` solely to reuse private map/setup data without widening the engine API;
it neither constructs a session nor calls the cutoff ladder. The CLI selects it
with `BAANDI_REFERENCE_SCF=on` or `[scf] driver = "reference"`. Both are opt-in;
`driver = "legacy"`/an absent key with the switch off preserves the existing
numerical path. Environment `on` also selects reference with `driver="legacy"`;
explicit `driver="reference"` remains a request even with the switch off.
`[scf] mixing_beta` and `mixing_history` configure beta and m; their reference
defaults are **0.7** and **8** (m in 1..25). `reference_preconditioning="tf"`
enables Thomas-Fermi screening; `plain` is the default. The optional mixing
names `plain`, `broyden`, and `pulay` denote this one difference-history method.
There is no legacy mixer selection or damping/positivity recovery in this arm.
The config parser validates reference-only names and retains `mixing_input`,
then resolves `broyden`/`plain` to the historical `pulay` alias solely for the
shared preflight vocabulary. The reference driver always constructs its own
Broyden implementation; legacy input parsing is unchanged when the selector is off.

Only the own-words [density specification](../research/qe-study/qe-density-mixing-spec.md)
was used for QE conventions; no QE source was consulted. The method references
are D. D. Johnson, *Phys. Rev. B* **38**, 12807 (1988),
doi:10.1103/PhysRevB.38.12807; P. Pulay, *Chem. Phys. Lett.* **73**, 393–398
(1980), doi:10.1016/0009-2614(80)80396-4; G. Kresse and J. Furthmueller,
*Phys. Rev. B* **54**, 11169 (1996), doi:10.1103/PhysRevB.54.11169. The
implementation uses the specified difference-space least-squares equations;
it does not claim to reproduce every variant in these papers.

### Density and metric

For valence density only, in electrons/bohr^3,

\[
 \rho_G=N_r^{-1}\sum_r\rho(r)e^{-iG\cdot r},\qquad
 \rho(r)=\sum_G\rho_G e^{iG\cdot r},\qquad \rho_0=N_{el}/\Omega.
\]

Reciprocal vectors are columns of B, including 2 pi, with G=B m and integer
Miller triple m. Keep m1>0, or m1=0,m2>0, or m1=m2=0,m3>=0. This is a
*density* half sphere at four times the wavefunction energy cutoff, not the
smaller orbital half sphere. BAANDI basis order is retained locally. The
reference reader uses exact Miller keys and validates complete membership;
no positional correspondence or inferred QE shell-tie order is assumed.
A sphere reaching an FFT Nyquist frequency is rejected rather than aliased.

The smooth mask selects ngms history entries. All G=0 linear combinations are
carried, with zero metric weight. For ngms<ngm, the remaining dense entries
receive x+beta F and never history subtraction or TF screening. The current
CLI has only the 4x density cutoff, hence ngms=ngm; the general mask is tested.
NLCC core charge never enters the mixed vector. No mixing-step normalization,
clipping, or feasibility damping is performed. The existing atomic-start
kernel normalizes its superposition once (and clips initial Fourier ringing).

For the half sphere, excluding zero and high-G entries,

\[
 \langle a,b\rangle_{Ry}=8\pi\Omega
  \sum_{G\ne0}^{half,smooth}\frac{\Re(a_G^*b_G)}{|G|^2},
 \qquad \langle a,b\rangle_{Ha}=\tfrac12\langle a,b\rangle_{Ry}.
\]

Here Omega is bohr^3 and |G|^2 is bohr^-2. The factor two from the implicit
negative half is included. dr2 is this form on F=y-x, before any mixing or
history mutation. Each batch uses compensated local sums and **one** allgather
of scalar partials, summed in ascending rank order. Fixed rank counts receive
identical scalars; rank-count changes are permitted to change rounding. No
full density is gathered by the metric or Broyden module.

### Update and solve

At accepted mixing call k>1, slot (k-2) mod m receives the **unnormalized**
differences dF=F_previous-F and dx=x_previous-x; occupancy is min(k-1,m).
The previous input and residual are saved before the current update. k=1
has no differences and performs x_next=x+beta F. A redo or stop never calls
the mixer, hence never inserts/evicts history.

\[
 A_{ij}=\langle dF_i,dF_j\rangle,\quad c_i=\langle dF_i,F\rangle,\quad
 A\gamma=c,\qquad x_{next}=x+\beta F-\sum_i\gamma_i(dx_i+\beta dF_i).
\]

A's triangle and c share one reduction batch. The existing symmetric EVD
adapter (faer; MKL dsyevd on the optional MKL backend) supplies
A=Q Lambda Q^T and gamma=Q Lambda^-1 Q^T c. This is a symmetric spectral
factorization, different from the note's indefinite-factorization route but
solving the same system. cond2=lambda_max/lambda_min is recorded. A numerical
rank deficiency lambda_min<=p*epsilon*lambda_max, nonfinite result, or backend
failure stops with a named error; **no regularization, normalized secants,
pseudoinverse, automatic restart, or silently reduced history** is used.
Near-rank-deficiency rejection is deliberately stricter than a zero-pivot-only
check. Storage and numerical work arrays are reserved at construction.

TF acts only on the history-corrected residual, using G^2/(G^2+a), with
r_s=(3 Omega/(4 pi N_el))^(1/3), a=(12/pi)^(2/3)/r_s in bohr^-2. It does not
act on dx. The plain profile carries G=0 without special projection; TF's
zero factor applies to its residual, as in the stated formula.

### Driver (R-SCF-2, patch 2 part A)

`driver::run(StartInput, Profile)` selects the registered communicator or serial
collective; `run_with` accepts an explicit communicator for simulated-rank tests.
The CLI branches **before** legacy problem construction. `refstart::start` owns
the final-cutoff atomic charge/orbital initialization, local projector rows and
FFT planes. Neither the ordinary loop nor its verification constructs a
`ScfProblem`, `ScfSession`, cutoff ladder, legacy mixer or legacy eigensolver.
`map.rs` calls the distributed real Hamiltonian, plane-density and LDA kernels
explicitly, also at one rank. Previous output orbitals become the next start;
H products are not retained between calls. Full density/orbitals are gathered
only at delivery. Immutable basis/layout metadata remains replicated.

Only LDA-PZ / modified LDA-PW, real Gamma, fixed integral occupations and NC
atomic data are implemented. Missing atomic orbitals/charge, NLCC, PBE,
symmetry, density floors, spin, smearing, forces/stress/relaxation reject.
Hamiltonian/FFT/backend kernel implementations remain shared; legacy SCF
scheduling, reuse, mixing and solver switches are not consulted by this loop.
The default driver is unchanged. `mixing="broyden"` is accepted with either
reference selector, including general validation; `plain` and `pulay` remain
aliases. Legacy validation still rejects `broyden` when reference is off.

All schedule quantities here are **Ry** (solver boundary divides by two):

- First call: ethr=1e-2, all bands treated as occupied.
- First-step redo once if dr2 < ethr*max(1,N_el), with the **same density and
  potential**, empty Broyden history and the just-computed orbitals. The redo
  uses max(0.1*dr2/max(1,N_el),1e-13) and the actual occupied/empty flags.
- At step 2 reset the cap to 1e-2; from then on use
  max(min(previous,0.1*previous_dr2/max(1,N_el)),1e-13).
- Occupied bands use ethr; empty bands use max(5*ethr,1e-5 Ry). `ref_davidson`
  uses the full owned plane-wave diagonal: kinetic + nonlocal projector-pair
  diagonal + current mean local potential, all in Ha, without Gamma weights.
- Each call permits 20 corrections. As described in `qe-scf-cycle.md`, retry
  from returned orbitals if more than max(5,bands/4) flags remain, at most six
  calls; exhaustion fails explicitly. Every call and remaining flag is recorded.
- Estimator stop: dr2 < 2*`standard_threshold_ha`. Default **5e-9 Ha maps
  explicitly to conv_thr=1e-8 Ry**, matching the relayed dev08 Si128 input.
  This is an estimator threshold, not an absolute-energy proof.

The new map retains signed mixed density. LDA evaluates its absolute value,
as specified in the own-words SCF-cycle note; neither Broyden output nor its
charge is clipped, rescaled or damped. The direct energy audit uses the
nonnegative orbital output. Input double counting uses the original signed
input potential and output density. A low-cutoff synthetic negative-ringing
probe is therefore handled without changing the mixer or tolerances.

The QE-style energy per iteration is

`E_band - integral y(V_H[x]+V_xc[x]) + E_H[x_next]+E_xc[x_next]+E_ewald+descf`,

where `descf = -Omega/N_r sum (x_next-y)(V_H[x_next]+V_xc[x_next])` in Ha.
On a successful stop, x_next=y and descf=0. The ordinary result energy remains
the independently accumulated orbital functional of y (kinetic, ionic local,
nonlocal, output Hartree/XC, Ewald). The direct and variational energies are
both retained; intermediates must not be compared under the same label.
On exhaustion the delivered density is the last **output**, consistent with
its energy and orbitals, not QE's unconverged restart mixed density.

#### Independent verification and continuation

`session::assess_fixed_point` is the extracted, callable ADR 0043 assessment;
the legacy verifier calls it with its original map and retains its original
predictor/clipping/test overrides. No legacy ordering or numerical path changes.
The reference driver evaluates a fresh map at y, warm-started from the stop
orbitals, at the next ratchet threshold. It measures E_H[F(y)-y], the occupied
residual-square estimate, density L2 change and occupied eigenvalue change.
The same standard threshold and residual guard decide pass; strict certificate
criteria remain information, so `certified=false` and `verified=true` have
their established distinct meanings. The output state is never replaced by
the verification output. The reference solver's explicit empty-band policy
also applies during verification; no legacy solver switches are inherited.

ADR 0043 V5 applies: after a failure, discard the verification output, keep
Broyden history, mix the original input/output pair and continue, within the
configured iteration budget. Every attempt and the first estimator-stop index
are recorded. This is essential: the displaced Si8 validation first stops at
step 10 but needs step 12 to pass. `verify="off"` delivers the first estimator
stop and creates no passing certificate. No predictor/reuse/scheduling heuristic
is used. All verification maps and output materialization are in engine wall.

#### Records and memory boundary

`results.json.reference_scf.iterations` (also each ordinary history row's
`reference_scf`) records dr2 in Ry/Ha, initial requested ethr and final used
ethr, redo, every call's ethr/iterations/active-band widths and sums, H columns,
reduced orders, restarts, repair/drop counts and unconverged flags. Per-correction
Ritz changes are retained to diagnose threshold straddles across ranks.
Mixing records beta, history occupancy/slot, coefficients' norms and condition.
Energy keys distinguish output functional terms and the mixed-density
Hartree/XC terms used in the QE-style variational energy. Console reference
lines use `ethr`, `avg # of iterations`, `total energy` (Ry) and `estimated scf
accuracy`; Gamma has one k point, so a call's average is its iteration count.

Phase seconds separate start, solve, density, potential, audit, transform,
metric, mixing, energy correction, verification and output. `map_s` is inclusive;
do not add it to its phases. Output covers delivery plus CLI work before the
final JSON serialization; an external timer is required for the final write
and process lifetime. `peak_bytes_per_rank` is explicitly a conservative named
array payload envelope, **not RSS or allocator high water**. It includes owned
Davidson blocks/small matrices, reference start, fields, mixing history and
the current replicated output-materialization buffers; opaque backend/FFT/allocator
storage is excluded. On Linux, `peak_rss_bytes_per_rank` records `/proc/self/status`
VmHWM through output preparation, gathered once at the I/O boundary. On macOS it
is null (unavailable). External dev08 rank wrappers must measure the complete
process high water, including final serialization/write and MPI finalization.

Validation and remaining rank-tail exception: [patch-2 validation](../performance/lanes/RSCF/patch2-validation.md).
No speed claim or default promotion is made.

### Reference files and tolerance

`refscf-check INPUT --qe DIRECTORY --steps 6` reads `k{k}-out` and `k{k}-mix`,
recursing at most four directory levels to find exactly one charge-density.dat
per snapshot. It checks four-byte Fortran markers, both endiannesses, lengths,
Gamma/nspin, finite values, reciprocal basis, unique Miller entries and exact
membership. The dev08 probe supplied to this task reports Si128 ngm=191712,
5,368,052 bytes (=116+28*ngm), and rho0=512/Omega; those are **relayed probe
facts**, not newly measured here.

The initial input is (mix1-beta*out1)/(1-beta); beta=1 recovery rejects.
Optional `k1-in` cross-checks the recovered atomic input. Optional `k1-beta1`,
`k1-beta-1`, or `beta1` cross-checks out1 (a beta=1 output alone cannot uniquely
recover the initial density). Every later input is QE's preceding mixed file,
not BAANDI's replay result. The reader recognizes stdout, stdout.txt, pw.out,
qe.out, output.out, run.out in each k-mix directory, or k{k}-mix.out/.stdout
at the root. Missing stdout estimates are null, not invented. Printed '<'
is an upper bound; '=' uses half the last printed decimal's resolution.

Checks are ordered: file/cell/membership, optional initial checks, per-k dr2
when available, A/gamma solve, all G, zero mode, lowest ten distinct nonzero
shells, total charge. It exits at the first failed term and names it. Output
contains max/RMS errors in e/bohr^3, relative errors, charge in electrons,
shell G^2 in bohr^-2 and dr2 in Ry/Ha. Absolute/relative defaults (1e-12
e/bohr^3, 1e-9) are **provisional diagnostic tolerances**, adjustable by `--atol`
and `--rtol`, not a measured QE agreement claim. Real QE tolerances await the
rank-count/repeated-run evidence requested in [the dev08 protocol](../performance/lanes/RSCF/dev08-request.md).

Synthetic three-step replay (including reversed file order and simulated MPI3)
measured max density difference 2.168404344971009e-19 e/bohr^3, maximum relative
shell error 2.168404344971009e-16, and zero charge error. This is a closed-form
linear fixed-point oracle, **not QE data**. The real Si8 fixture test emits
`SKIP QE si8 reference` while the delivery branch/files are absent. Detailed
Mac evidence and default-path fingerprints are in
[local-validation.json](../performance/lanes/RSCF/local-validation.json).

## Delivered QE regression (2026-11-23 follow-up)

Supersedes the earlier absent-fixture status. Copied all 46 files (2,806,532
bytes) byte-for-byte from the main worktree; no UPF is included. Provenance,
QE version, commands and hashes remain in the unmodified fixture README,
MANIFEST and SHA256SUMS. This is density-update correctness evidence, not
Mac performance evidence or independent full-SCF trajectory agreement.

`cargo test -p dft-engine --test refscf -- --nocapture` passes all 10 tests,
including the real production Broyden replay. Initial densities recovered with
beta 0.7 and 0.3 agree to 3.46945e-18 e/bohr^3 (relative max 1.17140e-16);
beta 1 output agrees to 3.83323e-20 e/bohr^3. Each iteration uses QE's own x/F;
history occupancy is 0, 1, 2. All 11,936 Miller indices are matched in BAANDI
norm/Miller ordering. The ten lowest nonzero shells are checked separately.

| k | max density error (e/bohr^3) | RMS density error (e/bohr^3) | G0 error (e/bohr^3) | max error over lowest 10 shells (e/bohr^3) | dr2 (Ry) | dr2 minus printed (Ry) |
|---|---:|---:|---:|---:|---:|---:|
| 1 | 3.46944695e-18 | 3.18984922e-20 | 3.46944695e-18 | 1.53329342e-19 | 0.438555331969 | 1.96919914e-09 |
| 2 | 3.46944695e-18 | 3.35435183e-20 | 3.46944695e-18 | 8.67361738e-19 | 0.019564039712 | -2.88015854e-10 |
| 3 | 8.67361738e-19 | 1.10392361e-20 | 0.00000000e+00 | 8.67361738e-19 | 0.000853988046877 | -1.95312339e-09 |

Charge differences versus QE are -3.74849e-15, -3.74849e-15, 0 electrons;
absolute total charge agrees with 32 electrons within 1.43e-14 electrons.
Maximum shell relative max/RMS errors are 6.4000e-16 / 9.6936e-16.
The original Si8-only regression used 1e-19 + 1e-13 times amplitude.
R-SCF-2 replaces it, as explicitly authorized by the task, with
**1e-17 + 1e-12 times the reference group's maximum amplitude** for same-rank
next-density replay (max and RMS), charge error < 1e-12 electrons and unchanged
G=0. Evidence is the relayed patch-1 dev08 floor at main 29beda5: Si128/Si8,
k=1..6, max 4.9e-18 e/Bohr^3, RMS around 1e-20; the old failure was only a
Si128 k=5 lowest-shell coefficient (1.6e-18 absolute, 3.8e-13 relative,
condition 6689). Patch 1 is passed. For **cross-rank QE density comparisons**
use 1e-16 + 1e-12 times amplitude: the measured QE Si8 np8/np1 floor is
1.7e-17..6.8e-17 in output and 7.3e-18..1.7e-17 in mixed density. These are
user-relayed measurements, not runs performed on this Mac. This replay floor
is not an SCF energy tolerance or permission to conceal a late Ritz-branch
trajectory difference. dr2 must match the printed rounded value within
half its last decimal place plus 1e-14 Ry, even though QE prints '<'.
CLI general diagnostic defaults remain unchanged; this fixture is tested at
the tighter measured limits. Native MPI/MKL Si128 remains a dev08 check.

CLI command (replace the replay TOML's UPF placeholder with the local file):

```sh
target/release/dft-cli refscf-check target/scratch/si8-qe-replay.toml \
  --qe tests/reference/qe/i25-si8/si8-np1 --steps 3 --atol 1e-17 --rtol 1e-12
```

Exit 0; no failing term. Machine-readable per-shell errors and diagnostics:
[qe-si8-validation.json](../performance/lanes/RSCF/qe-si8-validation.json), [qe-si8-cli.jsonl](../performance/lanes/RSCF/qe-si8-cli.jsonl).
Final full gate: exit 1, 813 passed / 3 known cache failures / 181 ignored;
build, fmt and clippy pass. The orchestrator explicitly accepts those TMPDIR
failures as an environment limitation. That patch-1 verification TODO is superseded by the Driver section above. No commit was created.

## Inner solver

R-SCF-DAV adds `dft_solvers::ref_davidson::solve(comm, operator, start,
hamiltonian_diagonal, options)`. It is a standalone call, wired directly into the reference driver by R-SCF-2. No existing solver, dispatch default, or SCF stopping rule
is changed. The algorithm organization is derived only from
`research/qe-study/davidson-iteration-structure-gamma.md`, not QE source.

`operator` implements `DistributedRealOperator`; one rank is the serial case.
`start` is local rows by requested bands, in `operator.local_rows()` order.
`hamiltonian_diagonal` is the full local-row plane-wave Hamiltonian diagonal in **Ha**, passed separately because
that trait intentionally has no diagonal accessor. Both use the existing
orthonormal real coordinates (ORC): the Gamma pair weights and single G=0
coordinate are already encoded. Never add another half-sphere/G=0 weight.
The caller guarantees a complete disjoint row partition and identical options
on every rank. Local invalid-input/H failures are agreed before proceeding.

`Options::new(n, ethr_ha, occupied)` selects capacity `min(2*n, global_len)`,
20 correction iterations, replicated reduced EVD, no returned H X, and clocks
off. `capacity`, `first_step`, `max_iterations`, `placement`, `return_hx` and
`profile` are explicit call-owned fields; no environment-selected algorithm.
A nondefault capacity must accommodate the first full correction block, i.e.
at least `min(2*n, global_len)`. `Solution` returns the lowest n values (Ha),
local eigenvector rows, optional local H X, per-band convergence flags and the
record. There is no output gather. A limit returns the approximate vectors
with unconverged flags; the caller must inspect them and decide on a new call.

Every call applies H to **all n start columns** and solves the n-dimensional
initial generalized problem, including first-step calls. This API has no
"already rotated start" bypass. The first correction attempts all n bands;
later corrections use only bands whose last Ritz-value change failed the
strict `<` test. Occupied thresholds are `ethr`; empty thresholds are
`max(5*ethr, 5e-6 Ha)` (the note's `1e-5 Ry`); `first_step` treats every band
as occupied. Flags are recomputed rather than locked. Residuals are constructed
with two selected-column GEMMs from V and HV, preconditioned, normalized, and
normally appended without projection. H/S borders alone are added. H is
applied only to accepted new columns. Exact/roundoff-zero directions cannot
be appended to S and are counted as dropped; the idealized 2n minimum H count
does not apply to such exact-start cases.

The implemented smooth positive denominator follows the 2026-11-23 addendum
"the diagonal preconditioner" in the study note on main:
`d_Ha = (1+x+sqrt(1+(x-1)^2))/4`, `x=2*(h_Ha-e_Ha)`.
The caller supplies `h(G) = |G|^2/2 + v0 + sum_atom,ij D_ij beta_i(G) conj(beta_j(G))`
(real part), with D and v0 in Ha. `v0` is the G=0 local-potential coefficient,
computed in snapshot replay as the average of the assembled total local potential;
that interpretation remains marked INFERRED in the study note. NC overlap is 1.
The complex plane-wave diagonal is shared by both real Gamma coordinates of a
pair; no extra pair weight is inserted. The legacy replay arm retains its
kinetic-only diagonal. Assembly is outside the timed solve and numerical loop.
The denominator uses a cancellation-free algebraic equivalent for negative x;
it tends to h-e at large positive x and 0.5 Ha at large negative x. Division
by d_Ha differs from division by d_Ry only by a uniform factor two, removed
by correction normalization. There is no fitted parameter. The confirmed
formula resolves the earlier provenance question; QE iteration identity still
requires matched dev08 measurements.

The reduced generalized solve reuses the existing guarded bordered Cholesky
helpers and `RealEighWorkspace` (faer, or the existing MKL backend). It
transforms H by the inverse overlap factor and returns the lowest n roots;
the available adapter computes a full symmetric EVD internally. One logical
reduced solve is performed at setup and after each nonempty correction block.
`ReducedPlacement::Root` broadcasts status, n eigenvalues and m-by-n coefficient
vectors; overlap factor and small transformed border remain replicated.
`Replicated` is the initial choice, **not a measured speed winner**. No sliced
EVD or rejected GEMM-shape/reduction experiment is reopened.

Contraction uses two GEMMs, V C and (HV) C, and no new H application. A necessary
numerical refinement to the note's ideal identity/diagonal reset is retained:
measure the fresh H/S borders from these stored blocks, refactor S and rebuild
the small transformed matrix, **without an additional EVD**. An early version
that assumed exact identity/diagonal matrices accumulated 1.0154e-11 Ha error
on a synthetic degenerate fixture; the measured refresh removes that drift
without relaxing the test. Those extra GEMMs/reductions are charged explicitly.
A failed overlap pivot invokes full old-basis whitening with the existing
factor and two-pass correction projection/rank filtering. This preserves the
old span and HV and is counted as a repair, not hidden as normal work. No H is
spent on rejected directions. Invalid/dependent starts fail explicitly.

`Record` owns setup, each correction iteration, and finish records: incoming
active band indices, H blocks/columns, reduced order and logical/backend EVD
calls, restarts, repairs, dropped columns, collective calls/input payload
bytes per rank, and GEMM `(m,n,k,2mnk)` shapes. Operator collectives are included
in the enclosing iteration; opaque FFT/nonlocal internal GEMMs are outside the
solver GEMM total. Payload bytes are not a network-wire traffic model. Eight
exclusive stage clocks are populated only with `profile=true`; total time
includes construction/glue. Traces and workspaces are reserved before the
iteration loop, with no I/O or serialization there. Backend-internal faer
packing can still allocate, as in the existing adapter.

Ritz changes are a stopping definition, **not an absolute eigenvalue bound**,
a residual certificate, or a spectral-completeness proof for an arbitrary
invariant start missing a lower state. Validation checks sorted low roots and
occupied subspaces independently, including a complete degenerate cluster at
the occupied edge. The SCF driver still owns independent final verification.
Mac evidence, the replay interface, exact counts and the outstanding dev08
comparison are in [the request](../performance/lanes/RSCF/dav-dev08-request.md)
and [stored validation](../../crates/dft-solvers/tests/ref_davidson_evidence.json).

## Inner solver cost accounting

P3-DAV-DENSE, 2026-11-24. This accounts for the reference solver at
`2aee083`, specifically `ref_davidson/workspace.rs`, not the older
`distributed_davidson/workspace.rs`. The reference solver shares the operator
trait and lean factor adapters; it does not call that older workspace.
Derived from the own-words notes
[davidson-iteration-structure-gamma.md](../research/qe-study/davidson-iteration-structure-gamma.md)
(including the preconditioner addendum),
[davidson-gemm-accounting.md](../research/qe-study/davidson-gemm-accounting.md)
and [qe-solvers-and-subspace.md](../research/qe-study/qe-solvers-and-subspace.md).
No QE source was read. Historical complex-solver counts in the latter two
notes are not measurements of this real reference implementation.

Let n be requested bands, m the current basis width, r this rank's ORC rows,
u incoming unconverged roots, k accepted corrections (k <= u <= n), and
q=m+k. Default capacity is min(global rows, 2n); explicitly larger capacities
are supported and the same formulas hold. All products are real, one
multiply/add counts as two FLOPs. A GEMM (a x b)(b x c) has 2abc nominal
FLOPs and a logical operand/output payload of
`8*(ab+bc+ac)` bytes (each operand once, aliases counted separately, output
write once). This is not measured DRAM traffic: adapter packing, tile rereads, zero fills, cache reuse
and backend scratch change physical traffic. A copy of t doubles moves
16t logical read+write bytes. Collective payload is input bytes per rank,
not wire traffic; all ranks replicate small matrices.

| Stage, baseline | Dense work per occurrence, per rank | Logical data movement / collectives | Algorithmic minimum and retained constraints |
| --- | --- | --- | --- |
| Corrections | Two GEMMs `(r,m)(m,u)`: 4rmu FLOPs. Pack selected C columns, precondition and normalize u residuals; no converged-band rotations. | GEMM payload `16*(rm+mu+ru)` bytes; C copy 16mu; residual pass reads X/HX/diagonal and writes E (32ru bytes, scalar caches excluded), normalization/compaction 16rk bytes. One sum of 8u bytes. | General C needs both products. Immediately after restart C=I, so stored V/HV columns already are X/HX; neither GEMM nor packing is necessary. Denominator, norm, zero-direction cutoff and accepted-root order must remain unchanged. |
| New borders | S: `V_q^T E`, H: `V_q^T HE`, each `(q,r)(r,k)`, total 4rqk FLOPs. | GEMM payload `16*(rq+rk+qk)` bytes; two sums of 8qk bytes each; upper triangle mirrored to both halves (48d read/store bytes for both matrices, d=qk-k(k-1)/2, counting the two diagonal stores). Appending E/HE adds 32rk copy bytes (charged to glue/apply, not the border clock). | Already incremental and one rectangle per matrix. k(k-1)/2 unused lower-triangle entries per product are computed, but splitting those products repeats the rejected upper-border experiment; not reopened here. S and its pivot decision precede H application, so their normal-step collectives cannot simply be fused without changing failure work. |
| Restart | Rotate V and HV by C: 4rmn FLOPs. Then fresh S/H: 4rn^2 FLOPs; factor S; two small refresh-transform GEMMs: 4n^3 FLOPs. No H call or EVD. | Rotation payload `16*(rm+mn+rn)` bytes; copies into V/HV: 32rn bytes. Fresh Gram payload `16*(2rn+n^2)` bytes and two sums totaling 16n^2 bytes, plus 48n^2 reset-store bytes, 16n^2 HR-to-H copy bytes and adapter-owned small-matrix/factor work. | Required physical rotations: two, not three. The note's ideal diagonal-H/unit-S reset eliminates the measured refresh. Historically rejected after 1.0154e-11 Ha degenerate-fixture drift under the old tighter regression premise; P3-RESET explicitly reopens guarded contraction under the strict 1e-10 Ha / standard 5e-9 Ha classes (below). Measured refresh remains the default and every repair/whitening still uses it. |
| Return | X=VC: 2rmn FLOPs; optional HX adds another 2rmn. | Per product `8*(rm+mn+rn)` bytes; no collective or gather. Owned local output buffers are returned. | Already minimal for a general final C. No unconditional final restart or reapplication of H. Optional HX is caller-required work, not removable duplication. |

Elementwise correction work, outside GEMM totals, is 2ru multiply/subtract
operations for HX-eX, ru residual divisions, 2ru norm multiply/adds and rk
normalization divisions, plus u square roots after reduction. Per entry the
unchanged denominator evaluates x=2(h-e), z=abs(x-1), one hypot(z,1), a
branch-dependent multiply/add, a reciprocal and final scaling. Hypot is a
library operation, not assigned a fictitious fixed FLOP count; comparisons,
abs, indexing and copies are not FLOPs. Mirroring borders, restart matrix
resets and return ownership transfer add no arithmetic FLOPs.

The `restart` timer charges the two large rotations, copies and resets only.
Fresh contracted Gram work belongs to `border`; refresh factor and transformed
H belong to `factor`. Thus the ticket's restart/border seconds cannot be read
as independent algorithmic causes or compared one-for-one with QE labels.
The ticket-relayed Si128 MPI8 numbers (8.500 / 5.495 / 5.953 / 3.874 s for
restart/border/corrections/return) motivate this work, not a reproduced timing
claim. Rank maxima from different stages need not add to critical-rank wall.

Setup costs n H columns, two n-square Gram GEMMs (4rn^2 FLOPs), factor and
initial reduced solve. Each ordinary step applies H to k accepted columns,
not u attempted ones. The default cap forces a restart after the first full
n-column expansion if another correction is needed. Exact-zero drops and
rank repair invalidate an unconditional 2n minimum-H claim.

For completeness, at a factor panel starting at p with width b<=64, the
adapter records W=2p^2b, Gram=2pb^2, U=2p^2b and inverse=2pb^2 nominal
GEMM FLOPs when reached. The first panel's small Cholesky/inverse uses the
existing library/factor adapter. Those scalar operations are not included
in `gemm_flops`. The reduced border transforms cost 4q^2k and Ritz
coefficients cost 2q^2n; full symmetric EVD backend work is not a GEMM count.
Repair separately whitens V/HV (4rm^2), transforms coefficients (4m^2n),
remeasures the m-square borders (4rm^2), refactors/transforms (4m^3 plus
factor), and projects the correction block twice (8rmu), followed by
in-panel two-pass dots/subtractions/norms. Its reductions and dropped
columns stay visible; repair must not be disguised as restart savings.

### Implemented switch and boundary

`BAANDI_REF_DAV_RESTART_CORRECTIONS=on|off` (input key
`ref_dav_restart_corrections`) is resolved once at reference solve entry and
**defaults off**. It affects only the first correction immediately after a
successful contraction: compute the unchanged preconditioned residual from
the selected retained V/HV columns directly. A solve-local flag is cleared
before any rank repair or reduced solve; fresh starts and subsequent ordinary
steps always use the general-coefficient path. Fresh S/H, guarded overlap
factor, active-band selection, thresholds, capacity, iteration limit, all H
applications and final rotations are retained. Two norm buffers (16n bytes)
are reserved once only when enabled; no new per-iteration allocation,
collective, logging or I/O is added. The existing backend may still allocate.

Per affected correction this removes **4rnu FLOPs**, two rotation calls,
16nu bytes of coefficient copying and `16*(rn+nu+ru)` logical GEMM payload
bytes. The residual still reads the same number of selected row entries,
now directly from V/HV, with the same row summation order and one u-norm sum.
These are arithmetic/data-movement counts, not predicted speedups. Existing
GEMM records omit the unexecuted `residual_x` and `residual_hx` calls, so the
saving is observable without inventing a new timing category.

Skipping identity multiplication can change signed zero behavior and depends
on the floating-point backend; no cross-backend bit-identity claim is made.
The switch stays off pending dev08. Fixed-layout tests compare every value's
bits, every incoming active set, Ritz operands, drop/repair/restart decisions,
iteration/H-column counts and communication. They also sum actual removed
GEMM records against 4rnu. Any dev08 divergence requires the first-branch
operands/margins and the existing rank-rule envelope; agreement of final
energies alone is insufficient. The old distributed and strict solver paths
never consult this switch.

The workspace's private buffers and methods lie outside the original ticket
ownership. No workspace copy/fusion experiment is implemented by this change;
removing its measured Gram refresh is deliberately excluded by prior evidence.
See the [dev08 request](../performance/lanes/RSCF/dav-dense-dev08-request.md)
for independent arms, full-run verification and remaining measurement gates.

### P3-DAVREM: full-diagnostic remainder attribution and copy removal

The new `convergence.work.davidson_stages.work` object is emitted only for
full reference diagnostics (`BAANDI_DAVIDSON_STAGES=on`). It contains `rank0`,
`per_rank` (rank-index order), and componentwise `max` stage totals. All
completed solves, including verification/retries, contribute. `rank0.solves`
is the denominator for per-solve rates; it is not the SCF map count.
There is no extra solver collective. At CLI reporting, two variable-length
allgathers collect rank records; their cost belongs to output, outside Davidson.
Shapes are retained separately on each rank, including zero-row ranks; maxima
never fabricate a GEMM by combining dimensions from different ranks.

Exclusive stages are setup, apply_h, apply_check, append, border_s, border_h,
factor, transform, reduced_eigen, coefficients, corrections, precondition,
normalize, repair, repair_panel, restart, return, decisions and glue. Nested
stages pause their parent. Thus the restart includes rotations and reset only;
its refresh is in border/factor/transform. The old eight clocks and QE-matched
clock boundaries are unchanged. The new reduced_eigen stage isolates the
backend plus status/eigenpair communication; transform and coefficients have
their own stages. Do not subtract this new eigen time from the old QE clock
and call the result the old remainder: the boundaries differ.

Every stage records seconds, invocations, logical GEMM calls and nominal
`2*m*n*k` FLOPs, scalar nominal FLOPs, sqrt/hypot counts, copy destination
bytes, zero-fill bytes, logical streamed bytes, and calls/input payload bytes
for barrier, fixed-tree sum, broadcast, allgatherv and alltoall. Collective
seconds are a subset of stage seconds, not an additive stage. GEMM shape
histograms identify stage, site, m, n, k and calls. These count adapter
invocations, not MKL internal blocks or pure-kernel microtiles.

`nominal_gflop_s = (gemm_flops + scalar_flops)/seconds/1e9` and
`logical_gb_s = (gemm_operand_bytes + streamed_bytes)/seconds/1e9` are computed
per rank before taking maxima. Copy bytes count destinations; their read and
write sides contribute twice that amount to streamed bytes. GEMM operand
bytes count each logical operand and output once. Neither metric measures
DRAM bandwidth. Backend packing, finite scans, hidden allocations, FFT work,
EVD internals and diagnostic storage are not modeled. The H and EVD clocks
remain useful, but their zero recorded FLOPs do not mean no arithmetic.

The preconditioner model charges 13 nominal scalar FLOPs/entry (upper branch;
the lower branch has one fewer multiplication), plus one hypot, and separates
normalization and its sqrt. The factor adapter uses a full-panel scalar/copy
model with width b and prefix p; rejected panels increment
`rejected_factor_panels`, so their modeled totals are upper bounds, not exact
executed scalar work. Successful-panel scalar FLOPs count Cholesky, inverse,
Schur subtraction, inverse negation and growth guards. Its scalar access
traffic is not a cache/DRAM model. Explicit repair dots/subtractions/norms
are counted from executed panel widths. No tolerance or decision reads these
counters. Off adds cached-flag branches only: no new clocks, TLS counter
access, work storage allocation, publication, or JSON work object. Existing
always-collected legacy records remain unchanged.

Non-GEMM inventory and decisions (r rows, n bands, m basis, u active roots):

| Work | Scaling and treatment |
| --- | --- |
| Start and workspace zero fill | Four r*capacity buffers, two r*n buffers, small matrices/vectors are allocated once per solve. Explicit workspace zeros are counted; shared adapter scratch is excluded. Cross-solve pooling is deferred because of ownership/lifetime changes. |
| Selected coefficient packing | m*u doubles once per general correction; already omitted by adopted restart-correction reuse when C=I. Retained, as noncontiguous selected roots need compact coefficients. |
| Preconditioning and norms | One r*u pass forms residuals, divides by the unchanged denominator and sums squares in row order. There is already no separate norm pass here. Normalization/compaction is a second r*k pass after the reduction and cannot precede the global norm. |
| Appending E and HE | Two r*k copies remain: operator and basis own separate contiguous matrices. Eliminating these requires an operator-view contract, outside this ticket. |
| Restart and repair rotation copies | Two r*n (restart) or r*m (repair) copies removed by `REF_DAV_SWAP_BUFFERS=on`: swap V/X and HV/HX, which have equal r*capacity backing storage. Numerical matrices, shape, column order and backend calls are identical. |
| Restart resets | `REF_DAV_LEAN_RESET=on` initializes only C=I. It removes five n*n stores (40*n*n bytes) overwritten by measured S/H refresh, complete factor/transform, or the next eigen/repair output. Fresh S/H and guarded factor are retained; this is not the rejected diagonal-H/unit-S approximation. |
| Small replicated work | H packing, S/T factor panels, transformed border, Ritz coefficients, reset C all execute on every rank; EVD placement remains explicit Root/Replicated. Separate clocks now expose this floor. Root-only factorization would require new broadcasts and fail-consensus handling; deferred pending its measured cost. |
| Repair | Two full projection passes plus in-panel dots/subtractions and a separate norm pass. Required when rank repair occurs; no reduction in pass count or threshold. Large post-whitening copies are covered by SWAP_BUFFERS. |
| Collectives | Norm reduction, S border, H-error consensus, H border, EVD status/eigenpairs. S/pivot must precede H, so fusing S/H changes failure work. Restart S/H could be packed together, but volume-only packed reduction was neutral in the previous lane; no speculative fusion is repeated. |
| Allocation/packing in iterations | Solver numerical buffers and diagnostic shape capacity are reserved before the loop. Existing EVD/backend packing may allocate; opaque backend work needs the requested dev08 stack profile. No new numerical per-iteration allocation. |

Both new switches default **off**, use the run registry/input precedence and
MPI agreement contract, and are resolved once per solve. No GEMM shape,
reduction order, H column count, active-root decision or threshold changes.
Mac correctness probes compare both independent arms and their combination,
full/off diagnostics, Root/Replicated EVD, cyclic partitions with 1/3/8 ranks
(including zero local rows), scattered/clustered spectra, capacity 2n and 3n,
and forced rank repair. See the request for commands and recorded counts.
The ticket's 9.25 s / 18 versus 5.29 s / 16 motivation is relayed one-run
MPI8 evidence; these Mac counters cannot assign that dev08 time to a cause.
No speedup or default adoption is claimed before the requested profiles.


### P3-STREAM: contiguous kernels and fused stores

Four independent `BAANDI_REF_DAV_*` switches default off: `SQRT`,
`FUSED_APPEND`, `CHECK_APPEND`, `REUSE_ROWS`. They affect only the reference
Davidson and use the run-owned switch registry. Derived from the own-words
preconditioner addendum in `research/qe-study/davidson-iteration-structure-gamma.md`;
no QE source was read. The off path retains the baseline arithmetic.

`SQRT` uses contiguous column slices and fixed 32-row stack chunks. Residual,
denominator division and squares are formed in one independent-row loop;
squares are then added in ascending row order from that small stack block.
This exposes vector sqrt/divide without reassociating the norm sum, changing
MPI reductions, or allocating in an iteration. The tail and zero-row ranks
use the same order. The pre-existing residual/division/norm fusion is retained;
there was no separate full residual or norm pass to remove.

For x=2(h-e), z=|x-1|, the Ha denominator remains
`0.25 * (2*max(x,1) + 1/(sqrt(1+z*z)+z))`. This is the baseline's positive,
cancellation-free rationalization, not the literal subtracting expression
for negative x. The changed operation is `hypot(z,1)` to `sqrt(1+z*z)`.
For finite computed x/z and z <= 2^500, round-to-nearest binary64 with unit
roundoff u=2^-53 gives sqrt relative error <=2u+O(u^2); the positive add,
reciprocal and final add give denominator error <=5u+O(u^2), relative to the
same computed x/z with exact remaining arithmetic. Scaling by 2 and 1/4 is
exact in this range. If the platform hypot has relative error <=2u (one
ulp), the two computed denominators differ by <=10u+O(u^2); the local
comparison gate is the conservative 8*EPSILON=16u relative bound. This is
an explicit floating-point model, not a portable libm accuracy guarantee.
Underflow of z*z is harmless: its contribution to 1 rounds away. The
common h-e/x/z rounding is unchanged and excluded from this bound.

A once-per-solve diagonal min/max and per-column eigenvalue check select the
sqrt kernel only when **every** row satisfies the bound. Extreme/nonfinite
columns use the original hypot kernel for the entire column; the ordinary
row loop has no data-dependent branch. Fallback is visible in precondition
hypot counts; sqrt and hypot counts include only rows actually executed.
Neither path changes band thresholds or silently accepts nonfinite norms.
The nominal sqrt model charges 15 scalar FLOPs/row, versus baseline 13 plus
hypot, with sqrt counted separately. The 32-row scratch stores are register/
stack temporaries, excluded from the logical matrix-traffic model. The sqrt
path removes the norm-vector zero fill; it writes each completed norm once.

`FUSED_APPEND` writes each normalized value to E and to the future V column
in the same pass, after the unchanged global norm reduction. It eliminates
the later E reread, saving 8*r*k logical bytes for k retained columns.
The second store is still counted as copy-destination bytes, attributed to
normalize. The append call consumes the staged-column marker. Capacity
exhaustion leaves the original basis intact and goes through repair; repair
invalidates any staged append and recomputes its own normalized columns.
It retains the existing exceptional normalization/copy passes. Exact-zero
corrections, compaction and failure decisions use unchanged tests.

`CHECK_APPEND` classifies finite HE entries while copying them to HV,
saving 8*r*k logical bytes per correction block. It retains the same
collective agreement, including operator errors and failures on another
rank. The copy is attributed to apply_check; no separate append scan runs.
Initial H validation stays separate because its producer belongs to the
operator layer. No finite check is removed or delegated without a contract.

`REUSE_ROWS` retains initialized V/HV/E/HE storage on the calling thread
between successful solves with identical (local rows, bands, capacity).
One cache entry is bounded to the most recent shape; a mismatch replaces it.
All active data are overwritten by start copy, H application, correction or
append before use. Factors, reduced matrices, operator state and returned
X/HX are never cached. No unsafe/uninitialized memory is introduced. Failed
solves drop their acquired buffers and publish no cache entry. On a warm
same-shape solve, setup avoids 16*r*(capacity+bands) zero-fill bytes; the
cold solve retains all fills. The two returned row matrices and shared
adapter scratch still allocate/initialize. Memory retention is a deliberate
tradeoff (four buffers per calling thread, released on thread exit or shape
replacement), recorded with the request rather than hidden as free storage.

No default is promoted. See the [one-run dev08 request and local attempt
record](../performance/lanes/RSCF/stream-dev08-request.md) for independent
arms, fixtures, rounding branches and validation limits. Timings await dev08;
Mac correctness and work-count evidence is not performance acceptance.

### P3-RESET: guarded ideal contraction (default off)

`BAANDI_REF_DAV_IDEAL_RESET=on` (input `ref_dav_ideal_reset = "on"`)
replaces only the contraction refresh with S=T=C=Y=I, H_r=A=diag(theta).
The rotated physical V/HV are retained. It is derived from the own-words
iteration note cited above; no QE source was read. The new premise is the
explicit P3-RESET work order and accuracy-class adoption rule, not evidence
that the historical unguarded error disappeared. Repair/whitening still
measures both Grams and factors/transforms them; no convergence threshold,
certificate or default changes.

A candidate reset checks global squared column norms in O(r*n) work with
one n-double tree sum; nonfinite norms or |norm²-1| > 1e-12 cause the measured
refresh and disable ideal resets for the rest of this solve. Every fourth
contraction instead performs the complete measured refresh and tests
max|S-I|, including off-diagonal entries, against the same bound. A trip
also disables subsequent ideal resets. The bound follows the existing
restart orthogonality contract and is dimensionless, not an energy class.
Four is a conservative experiment refresh interval, not a material-tuned
choice or a proof of accumulated-error control. Unit-norm correlations can
escape the norm test until the periodic check; independent verification is
still required. Guard state resets at each new solve.

`Record.resets` and each iteration's `resets` expose attempts, ideal,
periodic_refreshes, guard_fallbacks, disabled_refreshes and the maximum
periodic measured defect. A periodic check can also trip the guard, so these
counts overlap; actual measured restarts are attempts minus ideal. Existing
run-manifest path notes publish cumulative successful-solve counts per run
ID and calling rank/thread after each solve, outside numerical loops. Failed
solves publish no partial count. Fine-stage accounting charges norm scans and
their sum to restart. Each successful ideal reset saves the two full Gram
products (4*r*n² FLOPs), factor, two refresh transforms (4*n³ FLOPs), and two
n² sums, at the cost of a 2*r*n+n scalar norm check, one n-wide sum and five
additional n² reset stores versus adopted LEAN_RESET. The small-matrix stores
cannot be removed with LEAN_RESET on because they are no longer overwritten.

Mac fixtures: 24 A/B pairs (four fixtures, ranks 1/3/8, root/replicated EVD),
including the original 32-row degenerate DCT spectrum at capacities 16 and
24 and a seven-row/eight-rank empty-owner case. All converge, all active-band,
order/restart/repair/drop decisions and H calls/columns match. Maximum root
change is 2.886579864025407e-15 Ha; maximum doubly occupied band-sum change is
9.325873406851315e-15 Ha; maximum candidate/full dense-EVD root difference is
7.838174553853605e-14 Ha. These are fixed-H quantities, not SCF total energies.
A test-only full Gram measures every restart without changing its numerical
path; max|S-I| is 1.3880785409980945e-11 in both arms, before the measured
fallback. The original 2n degenerate case trips the norm guard immediately
on all tested ranks, so it does not demonstrate unguarded ideal-reset
accuracy; the same degenerate spectrum at 3n performs five ideal resets and
one periodic refresh on all ranks. No tolerance was loosened to hide that
limitation. Injected scaled and unit-norm correlated bases prove diagonal
and periodic off-diagonal fallback, sticky disablement and nonfinite rejection.

Full per-fixture counts, reproducible restart defect traces, gate results and
the exactly-one-profiled-run-per-arm MPI8 request are in
[the P3-RESET attempt record](../performance/lanes/RSCF/reset-dev08-request.md).
A separate temporary SCF probe also verifies cubic/skew/elongated two-atom
synthetic cells at ranks 1/3/8: all 18 off/on runs converge and pass independent
fixed-point verification, maximum total-energy change 2.7533531010703882e-14 Ha,
with unchanged SCF/Davidson/H counts. Source, reproduction and results are in
the attempt record. Real MPI/MKL accuracy and dev08 performance remain unmeasured. The conservative guard's fallback rate must
be reported; passing energy tests alone is not a speed or adoption claim.


### P4-DAVREM: triangular borders, stacked rotations and setup fusion

Three independent registry switches: `REF_DAV_TRIANGLE` (default **on** since
2026-11-25, dev08 one-run screen at 2432530: si216 MPI8 border stages 9.95 to
8.35 s, wall -1.8 s; si128 -0.26 s; map, iteration and H counts identical;
energies within 5e-13 Ha), `REF_DAV_STACK_ROTATE` and `REF_DAV_SETUP_FUSED` (default **off**) (environment prefix `BAANDI_`,
input keys lower case). No threshold, capacity, rank rule, collective ordering
or H-column policy changes. Derived from the own-words
[davidson-gemm-accounting.md](../research/qe-study/davidson-gemm-accounting.md)
and this section's reference-path accounting, not QE source.

**Triangle.** MKL computes a top m-by-k rectangle and a k-by-k upper triangle
of the appended border, using one DGEMM and one DGEMMT. For an initial/full
refresh (m=0), only DGEMMT is needed. The S pivot still precedes the H call;
failed corrections never incur speculative H work. The complete upper H
triangle uses V and HV as distinct operands, not a SYRK substitution.
The lower triangle is never read from GEMMT; it is zeroed before the existing
full-border reduction, then the upper result is mirrored as before. Pure Rust
and `MKL_DGEMM=off` retain the original rectangular product and report zero
savings. API/storage contract:
[Intel GEMMT reference](https://www.intel.com/content/www/us/en/docs/onemkl/developer-reference-c/2026-0/cblas-gemmt.html).

Per S/H pair, work falls from `4*r*(m+k)*k` to
`4*r*m*k + 2*r*k*(k+1)`: **2*r*k*(k-1) FLOPs saved**. Border calls increase
from two to four when m>0 and remain two when m=0. The top-output
scatter moves 32*m*k bytes per pair versus 32*(m+k)*k in the adopted single-GEMM
baseline; both adapter scatters remain excluded from the existing stage model.
Lower-zero stores add 8*k*(k-1) bytes per pair and are counted. GEMMT records retain square
m/n shapes, use site suffix `_gemmt`, and carry triangular FLOPs in the
iteration and stage totals. A histogram consumer must use `m*(m+1)*k`
for that suffix, not dense `2*m*n*k`. GEMMT is absent from the older dense
`gemm_dispatch` histogram; MKL's optional adapter timer includes its actual
triangular FLOPs under `other`. No hardware-internal FLOP claim is made.

The earlier UPPER_BORDER experiment saved 16.38% FLOPs but lost 0.57% wall,
with 12.7% more calls and narrower DGEMMs losing throughput. This explicitly
requested new candidate changes the square kernel to GEMMT (no recursive
upper tiles), not merely the old tile width. It still adds a rectangular call
and might lose for the same reason. The dev08 request measures that risk;
FLOP savings alone do not establish speed. Combining S/H into V^T[E|HE]
before the pivot would violate failed-append H-work identity, so is not used.

**Stacked rotations.** Pack the separate V and HV allocations into a persistent
2r-by-m matrix; multiply by C once, then unpack X and HX. Restart always uses
this route; return uses it only when HX is requested. An X-only return retains
its single original product. MKL calls one sequential DGEMM directly into
contiguous output, with no 256-row call subdivision. The pure backend uses
its existing ordered kernel. FLOPs remain **4*r*m*n**, logical calls fall
**2 to 1**. Explicit packing/unpacking adds **32*r*(m+n) bytes** per rotation
pair. Extra initialized storage is **16*r*(capacity+n) bytes** per solve,
counted in setup zero fills. No scratch allocation occurs in the iteration.
The logical GEMM payload saves one C read (8*m*n bytes); that does not erase
the packing cost or imply a net bandwidth saving. Repairs retain their
original whitening route and buffer swaps.

**Setup fusion.** After local shape validation, copy/finite-check the start
and optional supplied HX in one traversal, and finite-check/compute diagonal
bounds in one traversal. The unchanged two-value collective still agrees
invalid inputs and optional-HX presence before any operator call. Invalid
local shapes allocate no workspace. Setup saves **8*r*n + 8*r bytes** per
solve with the adopted sqrt preconditioner; supplied HX saves another
**8*r*n bytes**. With sqrt off the bounds pass was absent, so only the start/HX
read savings apply. Accounting now includes the previously omitted baseline
supplied-HX finite scan. Numerical results and decisions are bit-identical.

Remaining streamed traffic, with adopted P3 fusions, is explicit:

| Stage | Remaining logical bytes; no physical-DRAM claim |
| --- | --- |
| Setup row storage | Warm row cache leaves 16*r*capacity bytes of X/HX initialization; cold cache adds 16*r*(capacity+n). Add small arrays, adapter-excluded work, and 16*r*n start copy + 8*r diagonal read (+16*r*n if supplied HX). Stacked scratch adds the initialization above. |
| Normalize / append E | 24*r*k per ordinary accepted block: read raw E, write compact E and V. Global norms must precede normalization; removing another store requires changing the operator's owned-matrix contract. Repair retains its exceptional passes. |
| Apply check / append HE | 16*r*k per block: one HE read with finite classification and one HV store, then unchanged collective error agreement. A separate append scan is already absent. |
| Initial apply check | 8*r*n finite scan for newly produced H(start). Its operator owns production; this ticket cannot fuse that producer. Supplied HX uses the setup copy/check instead. |

Divide accumulated `streamed_bytes` by 1e9 for GB. The ticket's si216 setup
5.5 GB and 2.11 s are historical relayed observations, not enough to infer
new total GB without executed shapes/solve counts. The request requires
all-rank stage totals, same-rank GFLOP/s and GB/s, extra RSS and decisions.
See [P4 dev08 request](../performance/lanes/RSCF/p4-davrem-dev08-request.md)
for validation evidence and exact one-run arms. No default is promoted here.

## Rank-count acceptance for the reference driver (orchestrator, 2026-11-23)
R-SCF-2 validation: 24 system/rank combinations converge and pass the independent fixed-point verification, maximum energy difference to the strict legacy reference 1.69e-9 Ha (inside the standard class 5e-9 Ha); legacy default bit-identical. One combination (si16, 3 simulated ranks) took a different discrete branch in the last Davidson call (2 versus 3 iterations: a Ritz-value change sits at the threshold, and the rank-dependent rounding of the reductions decides which side), ending with a final density difference of 8.2e-8 e/bohr^3 while still converging and certifying. QE has the same property (its own np8 versus np1 floor flips discrete decisions near thresholds). Acceptance rule from now on: across rank counts, the delivered result must stay inside the accuracy class and pass the certificate; trajectories must agree within rounding EXCEPT at documented branch points where a discrete test sits within rounding of its threshold (recorded with the margin). Bit-identical trajectories across rank counts are not required for the reference driver.

Issue #26 / I26-1-RANKRULE makes this rule executable in
`crates/dft-engine/tests/reference_driver.rs` and its test-only
`common/rank_rule.rs`. Compare every layout to one rank, including verification
maps in chronological order. Before the first divergent predicate, require the
existing dr2 (`1e-11 Ry + 1e-8 relative`) and ethr (`1e-13 Ry + 1e-8 relative`)
agreement, Ritz-change agreement (`1e-11 Ha + 1e-8 relative`), identical incoming
active-band indices, correction work, repairs, drops, reduced orders, redo and
restart decisions. A missing trace or unexplained count difference fails.

At the first divergent Davidson correction, **every** differing band predicate
must straddle the strict `<` test. For each, retain both previous/current Ritz
values, both absolute changes, both actual band thresholds (including the
empty-band/first-step rules), and signed margins `change - threshold`, in Ha.
The test-only rounding envelope is
`32 * f64::EPSILON * max(1 Ha, abs(all four Ritz values))`.
Both margins, the cross-rank change difference and the threshold difference
must fit this envelope. It scales with the subtraction operands rather than the
small cancelled difference, is independent of material and rank count, and is
not an eigensolver tolerance or an absolute error bound. The test reports map,
ordinary/verification phase, call and correction (one-based), and band
(zero-based). A near-threshold event later in a trajectory cannot excuse an
earlier work or numerical mismatch.

Only after this evidence passes may subsequent trajectories/counts and the
cross-rank density maximum differ. Without an explained ordinary-map branch the existing
`1e-9 e/Bohr^3` density check remains: a verification-only branch cannot
explain a change in the density already delivered. In both cases **each** delivered run must
converge, pass its independent fixed-point verification on the final map, remain
within the standard `5e-9 Ha` energy class against the tight legacy oracle,
conserve charge within `1e-10 e`, and have finite delivered density/eigenvalues
and map energies. The strict `1e-10 Ha` class is unchanged. Failed verification
attempts, retries and recovery stay in the delivered wall and work; the matrix
prints total wall, summed verification wall and attempt/work counts.

Each of the six cases is run **twice at each of 1, 2, 3 and 8 simulated ranks**.
The fixed-layout assertion compares full density, orbitals, eigenvalues,
occupations, energy, verification and scientific histories/decision traces
exactly, stripping clocks only. This is repeatability at a fixed build/layout;
it does not promise cross-rank or cross-backend bit identity. Fast negative
tests reject a threshold crossing outside the rounding envelope, earlier work
divergence, and missing/nonfinite evidence. No default or numerical decision is
changed; the added production fields only record existing operands and work.

Release command (the external UPF is never committed):

```sh
DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine \
  --test reference_driver -- --include-ignored --nocapture --test-threads=1
```

`RSCF_MATRIX` emits one delivered-state row per case/layout and `RSCF_BRANCH`
emits the first divergent correction's per-band margins. `RSCF_RANKS` and
`RSCF_MATRIX_CASE` may narrow diagnostic runs; a full acceptance run leaves both
unset (one rank is always included as the baseline).


### I26-1-RANKRULE local release evidence

Release matrix **exit 0**, six tests passed (four fast acceptance tests and two
full-SCF tests), 24 case/layout combinations, each repeated twice. macOS arm64,
pure Rust, simulated collectives, rustc 1.99.1; base `c2a143d5` plus this patch.
Exact commands, exit code, UPF SHA-256, full precision values and delivered/verification wall
seconds are in [the evidence record](../../crates/dft-engine/tests/common/rank_rule_evidence.json).
These are local correctness observations, not native MPI/dev08 measurements or
QE reference results. No default is promoted or numerical threshold changed.

The following differences are against the tight legacy energy; verification is
the independent fixed-point Coulomb error of the delivered state. `PASS` means
energy class, convergence, final independent verification, finite state and charge
checks all passed; all rows also passed exact fixed-layout repetition. Margins
are `one rank / compared rank`. A dash means no discrete divergence, so no
rounding exception was used. Map/call/correction are one-based; band is zero-based.

| Case | Ranks | First map/call/correction (band) | Margins (Ha) | Energy difference (Ha) | Verification (Ha) | Charge error (e) | Delivered / repeat |
| --- | ---: | --- | --- | ---: | ---: | ---: | --- |
| si8 | 1 | none | — | +1.437996389e-10 | 1.776502657e-09 | +2.132e-14 | PASS / bits equal |
| si8 | 2 | none | — | +1.438138497e-10 | 1.776502656e-09 | +3.553e-14 | PASS / bits equal |
| si8 | 3 | none | — | +1.438422714e-10 | 1.776502658e-09 | +2.132e-14 | PASS / bits equal |
| si8 | 8 | none | — | +1.438849040e-10 | 1.776502657e-09 | +3.553e-14 | PASS / bits equal |
| si16 | 1 | none | — | +1.089119905e-10 | 3.087484819e-09 | -1.137e-13 | PASS / bits equal |
| si16 | 2 | none | — | +1.091464696e-10 | 3.087484818e-09 | -7.816e-14 | PASS / bits equal |
| si16 | 3 | 12/1/2 (31, SCF) | +8.648192699e-16 / -5.111474669e-17 | +1.002860017e-10 | 2.585725270e-09 | +0.000e+00 | PASS / bits equal |
| si16 | 8 | none | — | +1.094662139e-10 | 3.087484819e-09 | -1.421e-13 | PASS / bits equal |
| si32 | 1 | none | — | +3.232969448e-11 | 1.029302385e-09 | +3.126e-13 | PASS / bits equal |
| si32 | 2 | 12/1/4 (55, verify) | +4.873920448e-16 / -3.996389187e-17 | +3.754507816e-11 | 1.028119408e-09 | +3.126e-13 | PASS / bits equal |
| si32 | 3 | 12/1/3 (50, verify) | +2.098362887e-16 / -6.771946748e-17 | +3.761613243e-11 | 1.027209407e-09 | +3.126e-13 | PASS / bits equal |
| si32 | 8 | none | — | +3.809930149e-11 | 1.029302385e-09 | +2.274e-13 | PASS / bits equal |
| skew | 1 | none | — | +7.590550410e-11 | 9.689380252e-11 | +3.553e-15 | PASS / bits equal |
| skew | 2 | none | — | +7.594724849e-11 | 9.689372736e-11 | -2.665e-15 | PASS / bits equal |
| skew | 3 | none | — | +7.593570217e-11 | 9.689353486e-11 | -4.441e-15 | PASS / bits equal |
| skew | 8 | none | — | +7.593037310e-11 | 9.689345296e-11 | +1.776e-15 | PASS / bits equal |
| elongated | 1 | none | — | +5.967137895e-11 | 4.657857559e-09 | +8.527e-14 | PASS / bits equal |
| elongated | 2 | none | — | +5.965716809e-11 | 4.657857558e-09 | +9.948e-14 | PASS / bits equal |
| elongated | 3 | none | — | +6.069456049e-11 | 4.657857564e-09 | +8.527e-14 | PASS / bits equal |
| elongated | 8 | none | — | +6.046718681e-11 | 4.657857558e-09 | +9.948e-14 | PASS / bits equal |
| synthetic | 1 | none | — | +1.692470719e-09 | 3.366433472e-10 | +8.882e-16 | PASS / bits equal |
| synthetic | 2 | none | — | +1.692470386e-09 | 3.366433473e-10 | +0.000e+00 | PASS / bits equal |
| synthetic | 3 | none | — | +1.692470386e-09 | 3.366433473e-10 | +0.000e+00 | PASS / bits equal |
| synthetic | 8 | none | — | +1.692470608e-09 | 3.366433473e-10 | -2.220e-16 | PASS / bits equal |

si16/3 first diverged: the two Ritz changes were `1.98435712306377354e-12 Ha` / `1.98344118906845779e-12 Ha`, at thresholds `1.98349230379391778e-12 Ha` / `1.98349230381514659e-12 Ha`. Both signed margins and their difference fit the `7.10542735760100186e-15 Ha` rounding envelope.

si32/2 first diverged: the two Ritz changes were `5.04873920448289937e-14 Ha` / `4.99600361081320443e-14 Ha`, at thresholds `5.00000000000000015e-14 Ha` / `5.00000000000000015e-14 Ha`. Both signed margins and their difference fit the `7.10542735760100186e-15 Ha` rounding envelope.

si32/3 first diverged: the two Ritz changes were `5.02098362886727045e-14 Ha` / `4.99322805325164154e-14 Ha`, at thresholds `5.00000000000000015e-14 Ha` / `5.00000000000000015e-14 Ha`. Both signed margins and their difference fit the `7.10542735760100186e-15 Ha` rounding envelope.

Si16/3 has cross-rank density maximum `8.174729508592105e-8 e/Bohr^3`; this is
accepted only with the recorded threshold evidence and delivered-state checks.
All other layouts retain the density cap, including the two Si32
verification-only branches. Verification was actually retried in this matrix (the evidence includes every attempt count); failed attempts and
recovery remain charged to each run's total wall and H-column count. Each
layout also rejects a copy whose final verification is deliberately marked
failed, including the Si16 rounding-exception case.


Final full gate: `TMPDIR="$PWD/target/tmp" TMPDIR_990="$PWD/target/tmp"
CHECK_SCCACHE=0 CARGO_NET_OFFLINE=true scripts/check.sh full` exits **1**:
build, fmt and clippy pass; 832 tests pass, 3 fail, 191 are ignored.
The three existing CLI cache tests reject a cache inside a repository, while
this ticket requires scratch inside the worktree. A targeted rerun of
`pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`
exits 101 with the same explicit cache-location error. The other two are
`pseudo::cache_resolution_and_no_network_side_effects` and
`pseudo::cached_fetch_checks_hash_and_provenance`. No cache policy or test was
changed. The gate's informational load probe also cannot run `ps` in this
sandbox; the four gates still execute. README/dashboard/ledger files are outside
this ticket's owned paths and remain unchanged.

Default-path preservation was also checked against a pristine `git archive` of
`c2a143d5` built inside `target/tmp/pristine` with the same debug toolchain and
vendor dependencies. Archived Rust files were touched before building to avoid
Cargo reusing newer artifacts by mtime. Four identical CLI inputs (Si2, skewed
Si2, elongated Si2, Si8; 30 Ry; default SCF settings) have byte-identical complete
`density.bin` files and bit-identical parsed floating-point energy audits,
eigenvalues, occupations, convergence histories and independent-verification
records. Only clocks and the clock-derived coarse/fine wall equivalent were
excluded. Three cases converge and verify; elongated Si2 exhausts the default
budget identically in both builds (exit 3). This preserves that failure; it is
not a new convergence result. Default CLI runs do not emit orbital arrays, so
that comparison makes no orbital-file claim; the reference matrix separately
compares its full delivered orbitals. Exact inputs and hashes are in the same
evidence JSON. No solver/default arithmetic or threshold changed.

## Stop-to-verified interval

Integration includes CLI option wiring, timestamp JSON serialization and
clock exclusion in the existing rank-rule fingerprint, authorized by the
orchestrator on 2026-11-24. Verification and numerical defaults are unchanged.

I25 adds the additive `reference-scf-wall-1` contract, named by
`reference_scf.stop_to_verified_clock`. Every iteration (both the canonical
`reference_scf.iterations` array and its historical history alias) contains
`map_start_s`, `estimator_s`, `verification_start_s`, `verification_end_s`, and
`map_end_s`. These are seconds from **one output rank's monotonic driver-start
clock**, retained with both lightweight and full diagnostics. No MPI reduction,
barrier, rank maximum, logging or numerical-kernel timer is added. They do not
share an origin across ranks or with an external launcher clock. The new
contract is additive; old phase durations and scientific schemas are unchanged.

Map start precedes the threshold schedule and ordinary solve, including any
redo/retry. Estimator time is sampled immediately after the final ordinary-map
`dr2 < threshold_ry` decision. Verification start precedes cloning the candidate
orbitals and setting up the independent solve; verification end follows the
unchanged fixed-point density assessment and recording its verdict. Both are
null when no attempt ran. Map end includes failed-attempt mixing, energy
correction, bookkeeping and state handoff, but excludes final public output
materialization. Existing `wall_s.map` is a different, narrower solve-map cost;
it is not the difference between the new map boundaries.

The offline reader is `dft-cli refscf-check --stop-to-verified results.json`.
It prints one `stop-to-verified-1` JSON row with first stop map, passing map,
attempt count and seconds. It uses the stored independent verdict, never a new
energy-based passing rule. The interval starts at the first stopping
estimator's `estimator_s` and ends at the passing attempt's
`verification_end_s`. Its disjoint parts are:

1. **Verification attempts:** sum of every attempt's `[start,end]`, including
   failures and the passing attempt.
2. **Ordinary SCF after failure:** each failed map's `[verification_end,map_end]`,
   followed by each later map's `[map_start,verification_start]` when it attempts
   verification or `[map_start,map_end]` when it does not. This includes mixing
   after failure on the first stopping map. It is exactly zero on a first pass.
3. **Unassigned:** interval minus the two parts above, including the first
   estimator-to-verification gap and gaps between maps. No overlapping legacy
   phase duration is added. The reader checks ordering, complete attempt/verdict
   pairing and closure within `64 * EPSILON * max(1 second, interval_end_s)`.

If no attempt passes, status is `UNVERIFIED`, passing map is null, and the
interval ends at the last completed map (explicitly labelled, not time to a
pass). When no estimator stop occurred, interval fields are null and no closure
is claimed. Verification-off runs cannot be labelled verified. Incomplete
results, old records without the clock contract, missing boundaries, overlapping
maps, unordered attempts and maps after a pass are errors, not zero durations.

All driver work, failed attempts, recovery and all output remain in externally
measured delivered wall. The I26-2 full/lightweight diagnostic delta is a
**separate** paired-run measurement; it can overlap this interval and must never
be added to or subtracted from these parts. Local tests use independent exact
wall fixtures and a live first-pass synthetic control at both diagnostic levels.
The Si8 displaced full SCF is release-only on dev08 under the Mac test policy;
its test asserts attempts fail at maps 10 and 11, pass at map 12 and positive
recovery. This is a pending measurement, not new Si8 timing evidence. See the
[dev08 request](../performance/lanes/RSCF/stop-to-verified-dev08-request.md).

## Verification cost variants

I25-VERIFY-COST (2026-11-24) adds two independent, run-owned switches,
`BAANDI_REFSCF_VERIFY_REUSE=on|off` and `BAANDI_REFSCF_VERIFY_GATE=on|off`
(input `[performance] refscf_verify_reuse/refscf_verify_gate`). Reuse defaults
to **off**; the empirical gate defaults to **on** following the 2026-11-24
dev08 adoption (attempts ledger, I25-VERIFY-COST). Both affect only this
reference driver and can be combined for correctness checks. The legacy `VERIFY_REUSE`/`VERIFY_SCHEDULE` switches do
not select these experiments. The new first-attempt variant below remains experimental and default off.

The shared `assess_fixed_point` call is unchanged: a candidate's own output
state is deliverable only after a fresh map at that density, at
`next_ethr_ry(ethr, dr2, electrons)`, passes the Coulomb error test
`E_H[F(rho_out)-rho_out] < threshold_ry/2` and the existing occupied residual
estimate guard. The default threshold remains **5e-9 Ha**. Energy, density,
eigenvalues and orbitals all belong to that candidate, never to its passing
verification output. Prediction cannot mark convergence. Exhaustion, mixer
failure, and `verify=off` retain their explicit existing meanings.

### V1: consume a failed verification as the next iterate

For an ordinary pair `(x,y=F(x))` whose verification `z=F(y)` fails:

1. Run the ordinary Broyden update for `(x,y-x)` exactly once. This saves the
   pair and inserts its secant against the preceding pair in the existing
   circular history; do not clear, replace or regularize any history.
2. Choose `x_next=y`, overriding the computed mixed proposal. Use the actual
   real-space `y` and its already transformed reciprocal coefficients, with
   the usual fixed G=0 charge. Retain the verification map `z` and its
   eigenpairs/audit in memory. No new solve is done at this boundary.
3. The next logical SCF iteration consumes that map as `(y,z-y)` instead of
   evaluating a new map. Its ethr is the tolerance actually used in the
   failed verification, with **no second ratchet** or iteration-two reset.
   Its `dr2` is evaluated by the ordinary transform/metric path. From then on
   the ordinary next-ratchet rule applies to this iterate's ethr and dr2.
4. If this iterate continues, append `(y,z-y)` with the ordinary Broyden
   update, including normal history eviction and conditioning/error checks.
   If its estimator stops, independently verify **z** before delivery. A
   failed verification can itself be consumed next. As for all final
   converged iterates, no unused final mixing/history insertion is performed.

Energy accounting follows the input actually selected: for the source
iterate the existing correction `-integral (x_next-y) V_Hxc[x_next]` is zero
because `x_next=y`; the Hartree/XC terms are evaluated there. For the consumed
iterate, use its own map's band sum and double counting, then the ordinary
correction at its selected next input. Its map energy remains the energy of
its own output. Never reuse the source iterate's energy/audit for the new
iterate. The max-iteration budget counts logical SCF iterates, including
consumed maps; a failed final-budget check is not queued for nonexistent work.

The expensive map and its solver retries are counted/timed **only** in the
source verification attempt. The consumed iteration has no new solver calls,
zero map/solve/density/potential/energy work clocks, and `reused_from` points
to the source iteration. Its transform, metric, mixing, correction and any
new verification are charged normally. This preserves wall-interval closure
without counting the already computed map twice. Logical iteration count is
therefore distinct from evaluated-map count and solver-call count.

### V2: empirical gate with mandatory refresh

Let the last actual attempt have SCF estimator `d_a` (Ry) and measured
verification Coulomb error `e_a` (Ha). At a new eligible stopping candidate
with estimator `d`, predict `p=e_a*(d/d_a)` Ha. This estimates a locally
constant verification-error/SCF-estimator ratio. It uses no material, size,
rank or fitted constant. It is undefined until an actual attempt exists, or
if any input is nonfinite, `d_a<=0`, `e_a<0`, `d<=0`, or `d>=d_a`
(noncontracting relative to the calibration). Undefined means **attempt**.
No claim of an upper bound is made; the residual guard may fail independently.

If defined and `p>=threshold_ry/2`, skip this attempt and do ordinary SCF
mixing. After one skipped eligible candidate, the next eligible candidate
always attempts, even if predicted to fail, refreshing calibration and
preventing indefinite deferral by a stale model. Non-stopping iterations do
not consume that refresh. Every actual attempt updates `(d_a,e_a)`; passing
still depends solely on the unchanged verifier. This rule alone cannot save the first failed attempt without empirical calibration.

### V3: first-attempt secant amplification (I25-FIRST-ATTEMPT)

`BAANDI_REFSCF_VERIFY_GATE_FIRST=on|off` (input
`[performance] refscf_verify_gate_first`) defaults to **off** and is effective
only with `REFSCF_VERIFY_GATE=on`. It changes only the pre-calibration branch
of V2. The adopted error/estimator predictor still controls every candidate
after the first actual attempt. Both branches share the no-consecutive-skips
state: a skipped first candidate forces an actual attempt at the next eligible
candidate, regardless of its prediction. Undefined means attempt. `verify=off`
performs no prediction or verification.

Write `f(x)=F(x)-x`, `y=x+f`, and `J=df/dx`. Verification needs `f(y)`, locally
`f + J f`, not `J^-1 f`. The candidate predictors were evaluated as follows:

1. **Inverse Broyden:** its inverse-Jacobian model is
   `H=-beta I+(S+beta D)(D^T W D)^-1 D^T W`, with stored density/residual
   differences `S,D` and Hartree metric `W`. Direct `H f` is a Newton
   displacement, not the required residual. To estimate the latter, the
   evaluated prototype used Woodbury:
   `(D^T W S) alpha=D^T W f`,
   `r_hat=(1-1/beta)f+(D+S/beta)alpha`. It used scaled normal equations and
   the existing symmetric EVD adapter, with the existing numerical-rank rule
   and no regularization. The failed-first synthetic fixture's first prediction
   was **1.504757613e-10 Ha** versus actual **1.783487664e-8 Ha**; the first-pass
   control was undefined by the rank guard. No first skip or map saving.
2. **Forward directional secant fit:**
   `alpha=(S^T W S)^-1 S^T W f` with the same expression for `r_hat`, i.e.
   a least-change forward Jacobian with prior `-I/beta`. First predictions
   were **6.375661735e-10 Ha** (failed-first) and **5.666906701e-11 Ha**
   (first-pass, actual **1.649390463e-10 Ha**). Again no first skip or saving.
   Both directional prototypes were exact to 1.12e-16 Ha on a spanning,
   nonnormal linear two-dimensional oracle, but that did not predict the
   incomplete-history SCF behavior. These prototypes are removed, not hidden
   alternative production paths. Neither inverse application nor inversion
   is promoted on that evidence.
3. **Scalar dr2 contraction:** record `p_contract=(d/d_previous)*d/2` Ha
   for positive finite decreasing consecutive ordinary-map dr2 values. These
   belong to *mixed* iterates, not the unmixed verification. Once
   `d<threshold_ry`, this predictor is necessarily below `threshold_ry/2`.
   It cannot detect a first failure. Min with it suppresses useful skips;
   max adds no failing prediction. It is retained for audit only.
4. **Selected: observed unmixed-map amplification envelope.** Each Broyden
   secant already supplies `s=x_previous-x`, `q=f_previous-f`, hence
   `s+q=F(x_previous)-F(x)`. Unlike a mixed dr2 ratio, their norm ratio measures
   the response of F itself. Use the largest squared amplification in the
   retained history, including the current virtual secant:

   ```
   a_i = <s_i,s_i>_W
   b_i = <s_i,q_i>_W
   c_i = <q_i,q_i>_W
   L2 = max_i (a_i + 2 b_i + c_i) / a_i
   p_first = L2 * <f,f>_W / 2                       [Ha]
   ```

The exact identity `f(y)=F(y)-F(x)` follows from `y=F(x)`, while the
input separation is `y-x=f(x)`. Thus the needed error is the response of F
to a known separation, motivating scaling its squared norm by an observed
squared gain. For a scalar linear map this is exact. For multiple modes it
uses the largest *observed* gain rather than assuming the current residual lies in a reliably
reconstructed Jacobian subspace. It is neither a global Lipschitz bound nor a
certificate; unobserved directions, old secants, eigensolve noise and
nonlinearity can invalidate it. No material, size, rank, safety multiplier or
fitted constant enters. History length is the existing user's mixing history.
Both plain and Thomas-Fermi mixing supply raw secants, so the formula is
independent of the preconditioner. It replaces the two directional models
based on the synthetic evidence: first prediction **5.504054419e-9 Ha** skips
the measured failing candidate, while **3.951894732e-10 Ha** permits the
first-pass control. The passing refresh is badly overpredicted
(**1.755352352e-8** predicted versus **4.378277166e-11 Ha** actual), which
explicitly demonstrates why the mandatory attempt must remain.

Virtually evict exactly the slot the next ordinary mixer update would replace.
The predictor borrows the mixer immutably: no insertion, gamma/factor change,
iteration advance, density/orbital mutation or extra mixing. Two complex
vectors and 76 scalars are allocated/reserved once only when V3 is effective.
One batched scalar collective forms the three products per secant plus the
current norm (at most 76 products); no extra map, FFT, factorization, JSON or
I/O is introduced. All ranks inspect the same gathered products. With V3 off
there is no added predictor allocation or collective and default arithmetic
and work decisions are unchanged.

No history, any nonpositive density-secant norm, negative reconstructed output
norm (including cancellation), nonfinite product/ratio/prediction yields null
and an actual attempt. Collective/shape failures remain errors. A singular
multisecant matrix is irrelevant: no matrix is inverted in the selected rule.
The same model is evaluated on every eligible V3 candidate for audit, including
later actual attempts whose scheduling remains V2. There is no fitted waiting
count. This is different from the legacy LANE-VERIFY first-tail experiment
(2026-11-20 ledger: 16 wrong skips, 419 -> 444 maps, negative Si64 evidence):
that model used scalar mixed-iterate contraction; this one measures unmixed
output changes from reference-driver Broyden vector history.

Local counters and the staged protocol are in
[first-attempt-dev08-request.md](../performance/lanes/RSCF/first-attempt-dev08-request.md).
Two synthetic fixtures are insufficient for promotion. V3 remains off; no
strict-reference accuracy or dev08 speed claim follows from these observations.

### Records and validation boundary

When reuse or the gate is on, each iteration has a `verification_cost`
record: `gate_enabled`, `skipped`, `predicted_error_ha` (null when undefined),
`actual_error_ha` (null without an attempt), and `reused_from` (null unless
consumed). Effective V3 additionally exports `first_gate_enabled=true`,
`first_prediction_ha` and `contraction_prediction_ha` on every record. Predictions
are computed before the attempt on every eligible candidate; undefined values
stay null. `predicted_error_ha` is the value actually used for the decision:
secant amplification before any actual attempt, V2 scaling afterwards. A forced
refresh may run despite a predicted failure. Skips have no actual error; a
baseline-arm counterfactual is labelled separately, never inserted as a measured
attempt in the skipped arm. The CLI interval reader exports
`verification_predictions` with the three estimates, actual error, iteration
and skip flag. With V3 off these additive JSON fields are absent, preserving
the previous reference-history cost fields and numerical/work decisions. The existing
verification report retains only actual attempts;
skipped candidates have no verification call range or attempt timestamps.
The CLI exports these fields in both reference-history copies and an aggregate
`reference_scf.verification_cost` object with attempts/skipped/reused, logical
SCF iterations, map requests, evaluated maps and solver calls. Map requests
exclude bounded retries; evaluated maps include each retry (one solver call
per evaluation). Reuse is charged only at its source attempt.
With both switches off the optional record is absent and numerical arithmetic
and work decisions follow the previous path.

The stop-to-verified interval still starts at the first estimator stop,
including a skipped one, and ends at the passing verification's end. A
skipped candidate after a failure is ordinary recovery; pre-first-attempt
skipped work belongs to the existing unassigned interval part. The reader
accepts missing attempt clocks only with an explicit valid skipped record and
adds `skipped` when gate diagnostics are present. No rank maxima are mixed
into this output-rank interval. All output remains inside external wall.

Small-fixture and optional displaced-Si8 checks, exact counters, limitations,
and the pending dev08 protocol are recorded in
[verify-cost-dev08-request.md](../performance/lanes/RSCF/verify-cost-dev08-request.md).

## Setup and output cost

P3-SETUP adds coarse `Instant` clocks at **every** diagnostic level, including
`DAVIDSON_STAGES=off`, `PHASE_PROFILE=off`, and `STAGE_TIMERS=off`. No barrier,
rank reduction, JSON, logging or file I/O is added to a numerical loop. These
are rank-local elapsed seconds, including arrival waits in existing collectives;
they are not CPU time or synchronized inter-rank timestamps.

`reference_scf.setup_s` records CLI input read/parse/hash, the inclusive
preflight (cell/basis/grid plus UPF validation), the subsequent UPF
read/parse/hash, CLI remainder, and the driver setup. Preflight internals are
outside this ticket's ownership and remain explicitly grouped. `bootstrap`
records MPI initialization, arguments, initial switch-input parsing and dispatch.
`main_to_first_map` starts at the first statement in Rust `main`, not at exec
or mpirun: dynamic loading/process launch is outside that interval. Library
callers without `main` stamps receive null bootstrap/main intervals.

The engine phases distinguish species validation, basis/G-vector metadata,
Ewald, field FFT planning, owned rows, nonlocal tables, operator FFT plans and
real tables, the nonlocal diagonal, local/atomic coefficient tables and structure
factors, atomic density transforms/normalization, initial potential, potential
redistribution, and the atomic block. Its children separate raw orbital fill,
orthonormalization/completion, H application/Gram, and diagonalization/rotation.
Density-plan construction and mixer/metric buffers plus the initial transform
close the interval before the first SCF map. `atomic_block_total` contains its
four children; `collectives_nested` overlaps the phases. Never sum parents and
children or add nested communication to the wall total. Plan construction
necessarily includes its internal layout work. Tiny validation/clock overhead
can remain between coarse spans; retain the inclusive total as denominator.

`reference_scf.output_s` separates engine density gather, orbital gather/unpack,
result assembly, nested collectives, and CLI delivery/encoding/hash, requested
volumetric writes, mandatory density.bin/metadata writes, binary identity/run
manifest and reporting collectives (including their packing). CLI remainder
includes ordinary JSON-value/console construction, optional postprocessing and
profiling. Existing map `energy_s` remains separate; energy evaluation is not
part of setup or delivery and no energy formula changes here.

The reference writer serializes and writes the large result body once, then
appends `output_s`. `json_serialize_body` and `json_write_body` are actual measured
spans, including file creation in the latter. The small timing trailer, file
close, optional launch-profile patch, console delivery, MPI finalization and
process exit are excluded explicitly. This avoids a second whole-file rewrite
merely to time the first one. External process wall is still the acceptance
denominator. Every rank's existing results file carries its own setup/output
record; collect all rank files for min/max/critical-rank analysis.

Validation and remote campaign: [P3-SETUP request](../performance/lanes/RSCF/setup-dev08-request.md).

The engine output object also reports `driver_return_remainder`: host-measured
driver-call wall minus its last map boundary and engine delivery span. This
exposes state destruction after constructing the return value, plus the small
dispatch/argument overhead; it is not attributed to a numerical kernel.

`BAANDI_REFSCF_SETUP_PHASE_REUSE=on` is an explicit **default-off** experiment.
It retains the direct structure factor for one atom/G pair and reuses it for
the atomic-charge and ionic sums, and across the magnetic channels of each
atomic orbital. Multiplication, atom/column order, radial quadrature and
orthonormalization are unchanged; there is no phase factorization or radial
interpolation. Scratch is bounded by `16 * max(owned_rows, atom_count)` bytes
and included conservatively in the named-array peak envelope. The atomic
coefficient builder still works on owned G rows only. Real MPI/MKL timing and
identity remain dev08 gates even though the pure-Rust local tests are exact.


### P3-SETUP2: remaining work and exact reuse experiments

At base `0e5cfa6`, `REFSTART_MKL` and `REFSCF_SETUP_PHASE_REUSE` are already
**on** (the earlier default-off paragraph above describes their original
experiments). This task adds two independent **default-off** switches:

- `REFSTART_DIAGONAL_REUSE`: load each atom's owned complex projector values
  once per representative, then preserve the original coupling traversal,
  complex product and scalar accumulation. Copy the result for the adjacent
  ORC partner row only when its representative index is identical. No zero
  coupling is skipped. Scratch: 16 max(c_a) bytes, dropped with the sparse
  complex projectors; no full-G table or new collective.
- `REFSTART_ORBITAL_REUSE`: retain the first atom's exact `r(|G|)*Y_lm(G)`
  values and use them at subsequent atoms of that species. Scratch per radial
  channel: 8 L (2l+1) bytes, freed before QR; bypassed for a single atom. The
  phase multiplication, column order, radial quadrature and exact-bit shell
  keys are unchanged. No interpolation, harmonic approximation or shell
  tolerance. This scratch, like existing radial caches, is outside
  `StartMemory.principal_peak_bound_bytes`; it is reported analytically here.

Both preserve starting values and verified fixture trajectories bit for bit
on the tested pure Mac build. They do not change any legacy start. Exact
state equality on dev08 MKL/MPI remains a required measurement, not a claim.

Notation extends Start above: c_a projector channels on atom a,
Q=sum(c_a), C2=sum(c_a^2), A_s atoms of species s, M_c=2l_c+1,
U_c distinct *bit patterns* of squared momentum on this rank, R_c radial
mesh points, W=ceil(K/B), K=max(J,B). FLOPs use 2mnk per real GEMM.
Byte counts below are logical array payload/read-write traffic, not cache
misses, allocator peaks or MPI wire bytes. A tree reduction's stated bytes
are its local input payload (the result has the same size).

| Setup item | Arithmetic / scaling | Principal bytes and collectives; remaining excess |
| --- | --- | --- |
| Species validation | O(atoms + radial metadata); no dense algebra | Small replicated metadata; scalar agreement. |
| Basis/G-vectors and ownership | Enumerate full wave/density integer boxes, norms and sort/layout; O(G+D+F) metadata on every rank | Replicated basis/grid/index arrays, no field values; no collective in enumeration. Layouts are rebuilt, not cached between runs. |
| Ewald | Real/reciprocal candidate enumeration and accepted pair sums; exp/erfc are not assigned fictional FLOPs | 8 times accepted local terms for P>1, P broadcasts of 24 bytes in fixed rank order; serial streams. No algorithm change. |
| Field/operator/density FFT plans | Three distinct role-specific layouts; no transform FLOPs during planning | Layout-dependent O(F+G) indices, local sticks/planes and padded send/receive buffers. No plan-construction collectives. RustLanes creates a new FftPlanner per plan, shares equal axis lengths only within that planner; no cross-run cache. MKL descriptors likewise belong to each plan. Safe reusable immutable plans require adapter ownership, outside this ticket. |
| Nonlocal tables | Exact shell transforms per species/projector plus angular/phase evaluation on owned representatives (including mirrors needed by H) | Sparse complex 16 H_r Q then real 8 L Q; complex source released after diagonal. Scalar failure agreement. Already distributed; not G*Q on every rank. |
| Nonlocal diagonal | Baseline L C2 coupling terms, **2 L C2 binary searches** of the owned-index list. Real part of each complex term needs 5 useful real FLOPs including coupling multiply/add (full complex expression may execute more). | Baseline logical projector loads 32 L C2 bytes, coupling loads 8 L C2, output 8 L. D arm: H_r' C2 terms and H_r' Q searches/loads (16 H_r' Q bytes), where H_r' is this rank's representative count; ~half the terms and far fewer searches. No reductions in either arm. Not an all-G replicated diagonal. |
| Local/atomic field tables and phases | Exact-bit shell cache per species; O(U_s R_s) radial work, O(D_r atoms) direct phases and coefficient accumulation | 32 S_r output bytes; radial scratch O(R_s+U_s), phase scratch 16 max(A_s). One failure agreement. Existing phase reuse retained. Skew cells have many more U_s; interpolating them would change the start and is deferred. |
| Atomic density | Two inverse field FFTs, clipping/charge normalization, one forward FFT | Owned 8 F_r real and 16 S_r complex arrays; three transpose all-to-alls, each with padded send payload 16 P max_sticks max_planes bytes. Three scalar transform agreements plus one charge sum. G=0 fixed exactly. |
| Initial potential and redistribution | O(S_r) Hartree multiplier/metric geometry, one inverse FFT, O(F_r) LDA and pointwise additions | One padded field transpose, inverse-finiteness and XC scalar agreements; existing PlaneMap redistribution. No extra density FFT. |
| Atomic orbital fill | Baseline U_c radial transforms per channel (already memoized, **not recomputed at every atom**); sum_c A_s M_c L shell lookups/harmonics. O arm keeps radial transform count, reduces lookups/harmonics to sum_c M_c L for repeated species. | Raw 8 L J bytes; cache writes/reads 8 L M_c and 8 (A_s-1) L M_c per channel. Phase work unchanged from the existing phase-reuse choice. Two scalar failure/finite agreements. |
| CholQR2 Gram / rotations | Two Grams + two rotations: **8 L J^2 FLOPs**. Each Gram and rotation costs 2 L J^2. | Two reductions of 8 J^2 bytes each. Per Gram logical GEMM payload 8(2LJ+J^2); per rotation 8(2LJ+J^2), plus Q clone and reduced-result copies. Do not remove the second pass: conditioning/orthogonality contract still applies. |
| CholQR2 factor/inverse | Replicated O(J^3), currently the legacy complex scalar Cholesky/back-substitution helper even for real Gram data. Two passes each perform binomial(J,3) Cholesky off-diagonal complex MACs and J(J^2-1)/6 inverse complex MACs, plus diagonal norms, scaling and divisions. | O(J^2) real/complex factor temporaries; no collective. This serial cubic term does not shrink with ranks. A real/library triangular adapter is follow-up work; no numerical factor/guard change here. |
| MGS fallback / completion | Fallback: two-pass MGS, about 4 L J(J-1) dot/update FLOPs plus normalization. Completion: each attempted column against q accepted columns costs 8 L q dot/update FLOPs plus norm. | Fallback J(J+1) scalar sums; completion 2q+1 scalar sums per attempted column. Same dependence/pivot thresholds and kinetic ordering. Fallback was not reached in the recorded valid fixtures; no inference that it never occurs. |
| Atomic H application / Gram | Apply H to **K** columns in B-sized chunks, not only B: all atomic columns are needed to diagonalize in the same span. Gram **2 L K^2 FLOPs**. Nonlocal H about 4 L Q K + 2 K C2 FLOPs; FFT work is sum over chunk widths of ceil(width/2) inverse/forward complex pairs. | Q 8 L K, chunk input/output 16 L width, K-square Gram 8 K^2. W Gram sums totaling 8 K^2 bytes; if Q>0, W projector sums totaling 8 Q K. Batched H has two transpose all-to-alls/chunk plus two mirror exchanges if split; payload includes padding, not just useful G. No full HQ is retained. |
| Atomic eigensolve / final rotation | Replicated symmetric K-square EVD, backend-dependent O(K^3); selected rotation **2 L K B FLOPs** | EVD O(K^2); selected coefficient copy 16 K B read/write bytes; rotation payload 8(LK+KB+LB). No collective. Selecting only B input atomic columns changes the Ritz space and is not equivalent. |
| Mixer/metric initial transform | Driver-owned metric/history setup plus initial density transform; no new work here | Local field/history buffers and existing scalar reductions/transpose; recorded separately from start(). |

New nested setup clocks are `atomic_qr_gram`, `atomic_qr_factor_inverse`,
`atomic_qr_rotate`, `atomic_qr_mgs_fallback` (present only if taken),
`atomic_completion`, `atomic_h_apply` and `atomic_h_gram`. They are ordinary
setup-boundary clocks, with no barriers or reductions added; sum repeated
chunk/pass entries within each name. QR/completion children overlap
`atomic_orthonormalize_complete`; H children overlap `atomic_h_apply_gram`.
These and the pre-existing four atomic children all overlap
`atomic_block_total`. Do not add nested times to totals. A zero-duration
completion clock can simply mean no missing columns and clock resolution.

Mac observations, exact fixture counts, validation commands and the
**one profiled run per arm/case** dev08 request are in
[setup2-dev08-request.md](../performance/lanes/RSCF/setup2-dev08-request.md).
Derived from the initialization/shell-cache principles in the own-words
`qe-scf-cycle.md` and `qe-stage-by-stage.md`; no QE source accessed. A52's
rejected QR substitution is not repeated, and no numerical defaults change.

### P3-SETUP3: generalized atomic Ritz start and compact G lookup

Two independent switches are appended to the registry, both **default off**:
`REFSTART_GENERALIZED` (rounding class) and `REFSTART_COMPACT_G` (exact
metadata removal). They affect only `refstart`, not the legacy start or the
reference Davidson iteration. Derived from the starting-subspace principles
in the own-words notes `qe-scf-cycle.md` and `qe-stage-by-stage.md`; no QE
source accessed. Unlike the rejected A52 replacement of iterative QR, this
experiment changes only the initial Ritz representation and retains QR as a
conditioning/completion fallback.

Let X be the raw L-by-J atomic block on a rank, B the requested band count,
and K=max(J,B). All J atomic columns remain in the Ritz space. The low B
Ritz roots depend on all of them: selecting B atomic columns before applying
H changes the start and is not an equivalent saving. H remains applied in
ceil(K/B) chunks, whose widths sum to K; no K-wide FFT workspace is retained.

For J>=B the generalized path forms S=X^t X, uses the existing scaled,
guarded Cholesky inverse T with T^t S T approximately I, and forms
A=X^t H X. It solves (T^t A T) U=U e using the existing real symmetric
library adapter, then returns X (T U_B). This is A C=S C e with C=T U;
it requires no explicit orthonormal L-by-J block. One generalized EVD means
one symmetric EVD after reduction, not a newly implemented eigensolver.
Both small congruences and coefficient transforms use the existing GEMM
adapters, including REFSTART_MKL dispatch. Replicated small matrices are
**never** MPI-summed. Eigenvector signs/degenerate rotations are arbitrary;
compare Ritz energies, orthonormality and retained subspace, not raw signs.

| Step | Legacy dimensions/work | Generalized dimensions/work, J>=B |
| --- | --- | --- |
| Raw atomic fill | L x J; radial/harmonic/phase arithmetic | Identical |
| Overlap / QR | Two (J x L)(L x J) Grams and two (L x J)(J x J) rotations: 8 L J^2 FLOPs | One S Gram: 2 L J^2 FLOPs, one 8 J^2-byte reduction |
| Factor/inverse | Two J-square scaled/guarded factors | One existing J-square factor/inverse; same first-pass pivot guard |
| Metric check | Second wide QR pass | Small T^t S T, 4 J^3 FLOPs, no collective |
| H application | K columns, widths min(B,K-first) | J columns, same widths and operator capacity B |
| H Gram | (K x L)(L x width), 2 L K^2 FLOPs | (J x L)(L x width), 2 L J^2 FLOPs; ceil(J/B) reductions totaling 8 J^2 bytes |
| Reduced problem | One K-square symmetric EVD | T^t A T costs 4 J^3 FLOPs, then one J-square symmetric EVD |
| Selected coefficients | Copy first B eigenvectors | (J x J)(J x B), 2 J^2 B FLOPs |
| Final wide rotation | (L x K)(K x B), 2 L K B FLOPs | (L x J)(J x B), 2 L J B FLOPs |

Ignoring common H/fill/EVD and scalar factor work, the candidate saves
**6 L J^2 - 8 J^3 - 2 J^2 B** real GEMM FLOPs per rank for J>=B, one
J-square Gram reduction and one factor/inverse. This can be negative on a
small/empty owner; it is a measured experiment, not a guaranteed speedup.
At J=512/B=260 (Si128 if its manifest confirms those counts), wide QR work
8 L*262144 becomes overlap 2 L*262144, but the small checks/reduction add
1,073,741,824 FLOPs and coefficient rotation 136,314,880 FLOPs **on each
rank**. At J=256/B=132 the corresponding costs are 134,217,728 and
17,301,504 FLOPs. Record actual L_r from each run, not G/8 rounded up.
The legacy complex scalar factor helper is reused once; replacing it with a
library triangular adapter needs separately owned numerics work.

The first-pass pivot/scale rejection, or max-entry defect of T^t S T above
1e-12 (the adapter's existing unit-scale accuracy tolerance), sends the raw
block to unchanged CholQR2/MGS. No atomic rank is silently dropped. J<B
uses unchanged QR and kinetic-ordered completion without trying the new
path. Rank-deficient input still errors. A fallback adds the attempted S and
small-factor work; it is visible through overlap clocks followed by QR clocks.
Workspace row capacity is max(L,J) only for an eligible generalized start,
so empty/small owners can evaluate the replicated small problem too.

New setup children: `atomic_overlap_gram`, `atomic_overlap_factor_inverse`,
`atomic_overlap_metric_check`, `atomic_generalized_diagonalize_rotate`.
Successful candidates omit `atomic_qr_*`, `atomic_completion` and
`atomic_orthonormalize_complete`. The existing `atomic_h_apply_gram`,
`atomic_h_apply`, `atomic_h_gram` retain their scope. All overlap
`atomic_block_total`; do not sum nested clocks. The existing conservative
principal-array memory envelope remains; opaque GEMM workspace is excluded.

`REFSTART_COMPACT_G` omits the second FftGridSpec for the density sphere,
including its unused F-entry reverse map and D-entry forward map. It uses
the already validated grid dimensions, public canonical Miller mapping and
binary search of the FFT's sorted owned positions. This removes the S_r-entry
BTreeMap too, retaining the original density-basis order and exactly the
same owned G subset. Radials and structure factors were already evaluated
only at used owned G; no claim of removing all-rank numerical tables is made.
No interpolation, radial quadrature, FFT order or coefficient arithmetic changes.

Historical P3-SETUP3 ownership boundary (FFT part superseded by P3-SETUP4 below):
the driver constructs DensityPlan's full-grid
FFT and density basis again after start(). Both it and StartFields have
private mutable FFT buffers, with no public plan transfer/cache constructor.
DistributedPairFft also creates its own lanes/layouts; a dimensions-only cache
would omit rank, stick policy, backend and active-index identity. Sharing
immutable lane plans or a run-owned full-grid plan, and passing the existing
density metadata into DensityPlan, require changes outside OWNED_FILES
(`refscf/density.rs`, driver wiring, `distributed_fft.rs`, FFT adapters).
They are **not implemented here**. Wave/density basis enumeration and sorting
also remain in the basis crate; only redundant grid metadata is removed.
The exact follow-up boundary and one-run diagnostic request are in
[setup3-dev08-request.md](../performance/lanes/RSCF/setup3-dev08-request.md).

### P3-SETUP4: run-owned FFT kernels and atomic Gram

Three independent default-off switches, with no production default changes:

- `REFSTART_FFT_SHARE`: one constructor-scoped cache for the start and H
  plans, rustfft by length/direction and MKL by length/batch count. Identical
  kernels can serve different sparse/full layouts. DensityPlan copies the
  existing full-grid layout into fresh mutable buffers and clones shared
  kernel handles; no repeated partition or FFT plan construction. No global
  cache and no shared coefficient/transpose buffers. MKL's mutable descriptor
  API is protected by a mutex; pure rustfft retains private scratch and no
  execution lock. The unused zero-count MKL lanes allocate no descriptors.
- `REFSTART_SINGLE_GEMM`: atomic S and H Gram rectangles use the existing
  single-call MKL entry, guarded by REFSTART_MKL and the backend's MKL_DGEMM
  setting. Rounding class on MKL, unchanged pure arithmetic. Rotations already
  use the wide route. H already uses the same operator/apply_block as SCF,
  including batched pairs and adopted nonlocal projection/back-accumulation
  switches. Keep all J atomic columns and B-wide chunks; no separate H kernel.
- `REFSTART_LOCAL_SHELLS`: build exact-bit owned shell IDs/square roots once
  for every species; evaluate species-specific radial tables once per shell.
  Independent-row angle/accumulation passes preserve the original dot product,
  species/atom sums and quadrature. Both arms already operate only on owned G.
  Additional temporaries scale with owned rows/shells, never grid times atoms.

Sharing removes repeated planning, not remaining density-basis construction,
layout copying or buffer allocation. The existing field/operator/density plan
clocks retain these costs. S and all H borders each cost 2 L_r J^2 FLOPs;
single-call dispatch changes the reduction grouping, not that mathematical
work or the rank reduction. All existing nested setup clock boundaries remain.

The cache API restores its calling-thread constructor scope on exit, including
errors; returned plans retain their handles. It is not a cross-run cache.
Private constructors remain the default, and the reference driver's fork uses
the exact original rank, backend, support, axis and distribution metadata.
See [setup4-dev08-request.md](../performance/lanes/RSCF/setup4-dev08-request.md)
for the audit, correctness checks, limits and one-run/four-case dev08 screen.

### P4-START: triangular generalized start and construction reuse

Three independent **default-off** switches; existing defaults remain unchanged.
`REFSTART_CHOLESKY` factors the scaled overlap once with sequential faer LLT,
uses triangular solves for the metric check and H congruence, calls the existing
symmetric EVD once (MKL dsyevd when enabled), and back-solves only the retained
B eigenvectors before the wide rotation. No explicit inverse is formed. The
unit-column pivot guard 1e-6 and metric defect 1e-12 are preserved; rejected
metrics use unchanged QR/MGS, and J<B uses unchanged completion. Dense work
is replicated once per rank and never summed across ranks. The new factor
clock is `atomic_overlap_factor_solve`; all existing nested boundaries remain.
Faer handles the triangular operations on both builds; no new MKL FFI is added.

`REFSTART_SPECIES_FILL` forms owned radial/angular amplitudes once per species,
and direct atom/G phases once across its radial/magnetic channels. Adjacent
ORC partners reuse those values. Output ordering and expression order remain
unchanged. Temporary payload is 8 L M_s + 32 L bytes, with M_s the species'
magnetic-channel count, excluding radial maps. Local fields were already owned,
per-species and shell-cached; the P3t LOCAL_SHELLS no-gain experiment stays off.
Nonlocal diagonal and operator already share the sparse complex projector source.

`REFSTART_PARTITION_CACHE` extends the constructor-scoped PlanCache with one
exact geometry entry: shape, ordered indices, ranks, mirror and balance policies.
This removes two repeated wave partitions between owned-projector lookup and
operator construction. Private work buffers and rank selection remain separate;
geometry is backend-independent and kernels retain their existing backend keys.
Metadata copying/allocation and the first full-field plan remain. Existing
REFSTART_FFT_SHARE already covers field/operator kernels and the density fork;
full-field and sparse operator layouts cannot be exchanged. The extra O(G +
sticks + ranks) metadata is released at start return and excluded from the
principal numerical-array payload envelope.

Atomic H already uses SCF's batched pair FFT and nonlocal GEMM/back-accumulation
paths on the same retained operator. All J columns remain necessary; cold H
cannot use kept-HX without an earlier density/HX. The reported per-column gap
needs matched-width full-H counters, first-use allocation and wait attribution,
not a second kernel. CLI preflight and basis enumeration remain outside scope.
See [P4-START audit and one-run request](../performance/lanes/RSCF/p4-start-dev08-request.md)
for the five-arm si64/si128/si216 MPI8 diagnostic, provenance and validation.
No setup-at-QE or default-adoption claim is made from Mac checks.

## H-apply FFT cost accounting

P3-FFT (2026-11-24), audited against local base
`d7ec4e9565a24d2f88c5b241e4fdd597f50d2e73`. The ticket supplies a dev08
`66d45e69` Si128 MPI8 rank-max profile: inverse 5.80 s, forward 4.95 s,
local 2D transforms 2.12 + 2.14 s, alltoall 1.16 + 0.88 s, transpose
pack/unpack 1.37 + 1.10 s, ORC/grid/mirror work 1.49 s, multiply 0.43 s;
reported FFT excluding multiply 11.16 s versus QE FFTW 6.70 s (vloc_psi
7.98 s). BAANDI 14204 H columns are paired; QE reports 13592 FFTW calls.
These are supplied motivation, not newly reproduced reference results or
additive exclusive wall components: different rank maxima cannot be summed.
Raw provenance and transform-count definitions must accompany a comparison.

### Exact logical traffic, one complex pair

A pair carries two real bands (one on an odd tail). After deterministic axis
selection let `(n0,n1,n2)` be the FFT lengths, `R` ranks, `S` total occupied
sticks, `s_r` owned sticks, `Z` occupied reciprocal z planes, `L_r` local
complex coefficients, `p_r` owned first-axis real planes. Define
`A_r=s_r*n0`, `M=Z*n1`, `F=n1*n2`, `b=max(s_r)*max(p_r)`.
Complex f64 values occupy 16 bytes. The table counts value reads plus writes;
indices, write allocation, cache-line waste and FFT/MPI internal scratch are
excluded. FFT rows give **minimum interface traffic**, not measured DRAM
traffic or the butterfly traffic inside rustfft/MKL. Multiply by pair batch
width for a block; real-space rows already include all `p_r` planes.

| Current stage (fft_profile index) | Bytes per pair on rank r | Required work / avoidable part |
| --- | --- | --- |
| Coefficients to sticks (0) | `16*A_r + 32*L_r` | clear sticks, read coefficients, scatter; `16*L_r` clear writes are overwritten |
| Inverse first axis (1) | `32*A_r` | `s_r` contiguous length-n0 transforms, already pruned |
| Inverse transpose pack (2) | `32*A_r` | sticks -> destination-major send buffer |
| Inverse alltoall (3) | `16*R*b` outgoing payload | useful outgoing `16*A_r`; useful incoming `16*p_r*S` |
| Inverse middle unpack (4) | `p_r*(16*M + 32*S)` | clear dense active-z lanes and scatter received sticks; `16*p_r*S` redundant clears |
| Inverse second axis (5) | `32*p_r*M` | `p_r*Z` length-n1 transforms, already pruned |
| Middle -> slab (6) | `p_r*(16*F + 32*M)` | transpose active z lanes into contiguous z rows; `16*p_r*M` redundant clears |
| Inverse last axis (7) | `32*p_r*F` | `p_r*n1` length-n2 transforms; dense output required |
| Potential multiply (8) | `40*p_r*F` | complex read/write and real potential read |
| Forward last axis (9) | `32*p_r*F` | dense real input; same count as inverse |
| Slab -> middle (10) | `32*p_r*M` | gather only active reciprocal z planes |
| Forward second axis (11) | `32*p_r*M` | only active z planes needed |
| Forward middle pack (12) | `32*p_r*S` | gather only occupied sticks into send buffer |
| Forward alltoall (13) | `16*R*b` outgoing payload | useful outgoing `16*p_r*S`; useful incoming `16*A_r` |
| Forward transpose unpack (14) | `32*A_r` | received segments -> contiguous sticks |
| Forward first axis (15) | `32*A_r` | only owned occupied sticks |
| Gather and 1/N scale (16) | `32*L_r` | read selected G, write local coefficients |

Summing the table **excluding multiply and MPI payload** gives
`B_fft,r = 64*L_r + 144*A_r + 144*p_r*M + 80*p_r*F + 64*p_r*S`
bytes per pair. Hole-only initialization reduces this by
`16*(L_r+p_r*(S+M))` for unique G indices, the store minimum for these
same in-place buffers. Outer input fusion saves another `32*L_r` from the
combined ORC+FFT total, not from this FFT-only table alone. This is a
minimum for the retained layout/axis order, not a lower bound on every
possible out-of-place or strided FFT algorithm.

The table excludes the outer ORC packing and mirror exchanges. For a full
(two-band) pair, ORC packing reads 32 bytes per nonzero-G position and 16
bytes at G=0, then writes one complex coefficient: `48*L_r - 16*z_r`, where
`z_r` is 1 if this rank owns G=0. Let `H_r` be the number of owned nonzero
Gamma representatives (two ORC rows each). Unpacking reads plus/minus
complex values and writes four real entries: `64*H_r + 32*z_r` bytes.
The subsequent kinetic addition is separate: `32*Rorc_r` bytes per band
(read x, kinetic and y; write y), where Rorc_r is the owned ORC row count.

For split Gamma representatives let `O_r` be outgoing representative items,
`I_r` incoming items, and `u` the maximum items in a rank-to-rank block.
Per full pair the mirror-input staging clears `32*R*u` bytes and copies
`64*O_r` bytes; outgoing payload is `32*R*u`. Return staging copies
`32*I_r` bytes; outgoing payload is `16*R*u`. The returned plus/minus reads
are already in the unpack formula. These exchanges vanish when `u=0`, as
with mirror-paired ownership on symmetric supports. They are unchanged by
P3-FFT. For an odd tail, packing reads 16 bytes per nonzero G (8 at G=0)
and writes 16; unpacking reads 32 and writes 16 per nonzero representative
(16 read + 8 write at G=0). Mirror input copy/clear/payload per tail halve;
return complex payload is unchanged. Do not divide an odd block's total by
its pair count and call every pair a two-band pair.

Alltoall payload includes self blocks, **not** both send and receive traffic.
Both directions together send `32*R*b` bytes per rank. Across all ranks,
useful payload/direction is `16*S*n0`, padded payload is `16*R^2*b`;
utilization is `S*n0/(R^2*b)` when nonempty. Off-rank useful payload subtracts
`16*sum_r(s_r*p_r)`. Rank timings or Npw alone do not determine padding.

### Minimum and implemented experiments

The current code already scatters only into occupied sticks, transforms only
occupied sticks and active middle planes, uses one transpose per direction
for the **whole band-pair block**, and keeps backend plans/scratch. It never
zeros or copies a full 3D box per pair. Further generic pruning, alltoall
batching or plan reuse is therefore not a new mechanism. Transform work per
direction across ranks is `S*n0 + n0*Z*n1 + n0*n1*n2` point-transforms,
with conventional work estimate `5*sum(points_axis*log2(n_axis))` FLOPs.
The last dense pass is necessary for arbitrary local potentials; zero G
coefficients do not imply zero real-space values.

For the present in-place contiguous backend and padded collective layout,
a single scatter must initialize each input location either with its value
or with zero. Clearing **holes only** achieves that store minimum. Likewise
one transpose each way is already achieved at the collective level; packing
is still needed to reconcile destination-major MPI blocks and dense backend
lanes. Zero-copy transpose/strided out-of-place FFT would require another
backend contract, outside this ticket (especially `real.rs`). Removing
padding needs a variable-count collective, also outside owned files.

Two independent construction-cached switches, both default off:

- `BAANDI_FFT_SPARSE_CLEAR=on`: complement ranges for sticks, middle and slab
  are built once, coalescing adjacent holes into contiguous fills. It saves
  `16*(L_r + p_r*(S+M))` bytes per pair for unique coefficient indices.
  More generally replace L_r with the number of unique scatter destinations.
  The three existing profile byte counts reflect this reduction exactly.
  Every hole is re-cleared after reuse; buffers mutated by an in-place FFT
  must not be assumed to retain zeros. No per-point lookup table or new
  FFT arithmetic is introduced; range metadata/cache effects need dev08.
- `BAANDI_FFT_DIRECT_PACK=on`: the distributed H apply supplies ORC packed
  values to the first-axis stick scatter through an inlined callback. This
  removes the input `hbuf` write/read, **32*L_r bytes per pair**, including an
  odd tail. Output hbuf, mirror exchange/unpack and normalization remain.
  Old and direct paths use exactly the same scalar ORC expression. Both
  batched and diagnostic unbatched H application are supported; density and
  kept-change entry points retain their existing coefficient input path.
  This switch affects the distributed operator, including other drivers
  that explicitly opt in; it is not confined to reference SCF.

Savings add without double counting. All transforms, collectives, payloads,
axis selection and kernel choices stay unchanged. With switches off there
are no hole-range allocations; with on, metadata is allocated at plan
construction only. Transpose storage retains its existing grow-to-largest
batch behavior. H input hbuf capacity remains because forward output needs
it; this is a traffic optimization, not a peak-memory claim.

`FFT_DIRECT_PACK` moves ORC arithmetic from outer `pair_pack` into FFT stage
0; stage 0's profile still counts a logical 16-byte coefficient read rather
than callback source loads. Its abstract byte count does **not** independently
prove the 32*L_r removed traffic. Compare combined pack+scatter and full
apply totals, and derive saved traffic from L_r and the actual pair sequence.
No speedup may be claimed from a smaller renamed stage alone.

### Prior attempts, tests and limits

Derived from the principles in the own-words notes
`fft-call-structure.md`, `qe-hpsi-and-fft.md`, `qe-gamma-real-path.md` and the
engine's current code; no QE source accessed. A01/A02/A26 and M2d already
supply reuse/pruning/batching. A11/A28/A34, neutral M2_FFT, rejected async
alltoall, FFT-KEPT-FUSE-NEUTRAL and FFT-ARRIVAL-OWNERSHIP-CLOSED are not
repeated. PFFT-4's indexed transpose was worse on the old Mac path: this
candidate retains the transpose and uses contiguous complement **spans**,
not an all-element indexed transpose. The changed premise is the user's
new reference-driver Si128/Si64 dev08 profile and explicit P3-FFT work order,
not an inference that removed bytes guarantee wall improvement.

Tests in the two owned source files compare each switch alone and combined,
cubic/anisotropic/odd grids, asymmetric support, full/empty support, both
first axes, simulated ranks 1/2/3/8, dirty-buffer reuse and shrinking/growing
batches. Physics tests add cubic/skewed/elongated cells, odd/even/empty band
blocks and greedy/mirror ownership. Existing serial full-grid and pair-plan
oracles remain. These are correctness tests, not Mac performance evidence.
The [dev08 request](../performance/lanes/RSCF/fft-dev08-request.md) owns the
screen and A/A + crossed A/B protocol. Default promotion, real MPI/MKL,
full convergence and external-reference accuracy remain remote gates.

### P3-A2A: payload minimum and capped stick ownership

Audited at `ee6a0c40d7ce8c197e570681247bc6744ff45c67`. Derived from BAANDI
code and the own-words notes `fft-call-structure.md` and `qe-hpsi-and-fft.md`;
no QE source was read. Ticket-supplied ec553c93 Si128 MPI8 motivation:
206 alltoalls, 4134 MB outgoing payload/rank, 103 projector sums (58 MB),
about 2.0 s collectives, 15.9 s H apply, 28 s Davidson stage sum, 3.5 s FFT
movement and 4.9 s kernels. These rounded componentwise rank maxima are not
additive elapsed wall. The supplied QE h_psi 11.25 s / 82 calls is not new
reference evidence; this experiment supplies no QE run or speed ratio.
MB below means decimal 10^6 bytes; confirm the remote counter's formatting
before matching its rounded aggregate to the exact byte formulas.

#### Exact payload versus useful data

For `q=ceil(m/2)` packed pairs, the current collective sends
`16*q*R*s_max*p_max` bytes/rank/direction, including self. Useful inverse
send/rank is `16*q*s_r*n0`, useful forward send is `16*q*p_r*S`.
Their difference from the padded total is unused slots, not extra support.
Padding can contain stale values after reuse; receivers never read it.
Entirely empty reciprocal columns are **already absent**. Zeros within an
occupied stick are filled before the first inverse, after which all n0
spatial samples are generally nonzero. Filtering those transformed samples
by the original G-sphere would be incorrect. The forward route also gathers
only S occupied sticks, not the full n1*n2 rectangle.

For Hermitian spectra A and B of two real orbitals, C(G)=A(G)+i B(G), but
C(-G)=conj(A(G))+i conj(B(G)), not conj(C(G)). Both mirror sticks contain
independent data for the **two** real bands. The narrowest exact type for
this contract is therefore complex f64 (16 bytes/sample), not real f64 or
f32. Self-mirror sticks do not make C real. An odd tail has a zero partner
but still uses the full complex route; a separate Hermitian single-band
algorithm would need another contract and is not implemented here.

The global useful minimum for this transpose layout is `16*q*S*n0`
bytes/direction; mean/rank is `16*q*S*n0/R`. Only
`16*q*(S*n0-sum(s_r*p_r))` useful bytes cross rank boundaries. Counters
record outgoing send bytes including self, not send+receive or transport
protocol overhead. There is no second factor-of-two Gamma saving. A
variable-count collective could eliminate padding; the existing equal-block
interface instead pays the maximum stick and plane counts on every edge.

`a2a_catalog_payload_accounting` enumerates the 30 Ry = 15 Ha spheres for
the lattice/grid sizes of `benchmarks/cases/si{64,128}-gamma.toml`, then
executes the actual mirror partition at MPI8. This is local catalog
accounting without UPF/SCF, not proof that remote input hashes match it.

| Catalog, MPI8 | Grid | Npw | S | s_r, rank order | p_r | Sent bytes/pair/rank/direction | Useful mean bytes | Padding |
| --- | --- | ---: | ---: | --- | ---: | ---: | ---: | ---: |
| Si64 | 72x72x72 | 23871 | 1005 | 125,124,126,126,126,126,126,126 | 9 | 145152 | 144720 | 432 (0.2976%) |
| Si128 | 72x72x144 | 47851 | 2005 | 251,250,250,250,250,250,252,252 | 9 | 290304 | 288720 | 1584 (0.5456%) |

Si128 inverse useful bytes/rank/pair range 288000..290304; forward useful
bytes are 288720 on every rank. Across ranks/direction: 2322432 submitted,
2309760 useful, 12672 padding, 2021040 useful off-rank bytes. A full
260-band starting block (130 pairs) submits 37739520 bytes/rank/direction.
The rounded 4134 MB is consistent with 7120 pairs over 103 H calls
(`2*7120*290304 = 4133928960` bytes), but does not prove that pair count;
request actual widths. Even eliminating all padding would save only about
22.6 MB of 4134 MB. Both catalog layouts already attain the smallest
mirror-preserving equal-block capacity; the candidate leaves them unchanged.

#### Calls, latency and the earlier overlap rejection

With mirror ownership, each nonempty H block uses exactly two FFT alltoalls
and one projector sum (zero sums without projectors). Greedy ownership can
add two mirror exchanges. Unbatched diagnostics use `2*ceil(m/2)` FFT
alltoalls; density uses one inverse transpose outside these H counts. Empty
H blocks do no collectives. The low-level FFT API retains its historical
zero-width calls for nonempty support. No call order/count is changed here.

4134 MB / 206 is about **20.07 MB/call/rank**. Since 2.0 s includes sums,
alltoall mean residence is at most **9.71 ms/call**, not measured latency.
In `T = calls*alpha + payload/bandwidth + arrival/progress wait`, alpha,
bandwidth and wait cannot be separated from that aggregate. Pure latency
share is unknown: at most 12.58% of the H stage or 7.14% of the stage sum,
with the same rank-max caveat. Never call 9.71 ms startup latency or assign
all 2.0 s to alltoall. Per-rank FFT stage 3/13 clocks and existing arrival
sites 104/114 are needed to separate arrival imbalance from residence.

MPI-REMAINING-ASYNC-ALLTOALL rejected overlap at Si32 sizes. Si128 messages
are larger, but the aggregate still does not show alltoall alone above that
analysis's roughly 10%-of-solve reopening premise. `dft-mpi::Collective` and
World only expose blocking equal-count alltoall: no alltoallv, Ialltoall,
request/test/wait or progress guarantee. Overlap needs the wrapper, safe
buffer/request lifetimes, double buffers, error/unwind and communicator-order
rules, simulator tests, and forwarding through accounting/engine wrappers.
A new trait method with a padded fallback would not remove traffic through
those wrappers. No FFI is added. The old rejection is not silently reversed;
Si128 asynchronous progress and overlap benefit remain unmeasured.

#### Implemented reduction: FFT_A2A_BALANCE (default off)

The switch minimizes the largest removable part **within the existing
equal-block transport and mirror constraint**. Construction first performs
the unchanged coefficient-weighted greedy assignment. With the switch on,
form indivisible size-1/size-2 mirror units (only size-1 under greedy mode).
For K size-2 units, the smallest feasible capacity c obeys `R*c >= S` and
`R*floor(c/2) >= K`: ceil(S/R), or that value plus one. Assign pairs before
singletons, decreasing coefficient weight within each class, to the feasible
rank with least coefficient load (ties: stick count then rank). This secondary
objective limits, but does not guarantee against, coefficient-load imbalance.

Use the candidate **only if c < old s_max**; otherwise preserve old owners
exactly, including empty support and single rank. Savings per rank/direction
are `16*q*R*(old_s_max-c)*p_max` bytes, also removed from each transpose
buffer. Plane padding and indivisible-mirror padding remain. Savings can be
zero, as on both catalog Si cells. Existing collective/FFT counters record
the smaller buffer exactly. No per-call allocation, logging, new backend,
kernel, collective or batching policy is introduced; metadata is built at
partition construction only. The switch affects all partition users,
including density/operator rows, not just reference SCF.

FFT arithmetic and spatial-plane order are unchanged; assembled transforms
and local+kinetic H must be bit-identical. Projector/Gram sums and SCF can
round differently when coefficient rows change owner: this is a
**rounding-class** experiment. Check remote counts, verification and load
imbalance explicitly. Lower bytes are not a speed result. See the
[one-run-per-arm request](../performance/lanes/RSCF/a2a-dev08-request.md)
for local validation and pending remote evidence. Shared ledger/log,
README evidence and dashboard integration remain the master's work outside
this ticket's OWNED_FILES.

## Verification by energy

ADR [0044](../decisions/0044-reference-verification-energy.md) changes only the
reference driver's default verification judgment. `REFSCF_VERIFY_ENERGY=on`
requires `abs(E(rho_out -> F(rho_out)) - E(rho_in -> rho_out)) < 5e-9 Ha`,
using each map's `energy.band_form`: band sum minus the output-density integral
of the input Hartree/XC potential, plus output-density Hartree/XC energies and
Ewald energy. Neither value includes the subsequent mixed-density `descf`
correction. The delivered state stays the stop map's state. The independent
map, its tolerance, all diagnostics, failed-attempt recovery and work charging
are retained. Coulomb error and the empirical residual estimate are diagnostic
only; no density fixed-point accuracy guarantee follows from an energy pass.

Both verification predictor gates are inactive in energy mode, so every
estimator stop attempts a fresh verification (unless verification is explicitly
off). Results label the mode `energy_change`, store both variational energies
and their absolute difference per attempt, and retain the old diagnostic
fields. `refscf-check --stop-to-verified` accepts either judgment without
reinterpreting its recorded verdict. `REFSCF_VERIFY_ENERGY=off` restores the
old acceptance rule and configured gates; the legacy driver is unaffected.

## QE-matched clock table

P3-QECLOCK, diagnosis only. This map uses the own-words notes
`research/qe-study/qe-stage-by-stage.md`, `qe-scf-cycle.md`,
`qe-hpsi-and-fft.md`, and `davidson-iteration-structure-gamma.md`; implementers
have not read QE source. Scope is NC, unpolarized real Gamma, fixed occupations,
no symmetry, reference driver. Exact means the same mathematical operation
and enclosing work, not identical algorithms, number of calls, or trajectory.

With `BAANDI_DAVIDSON_STAGES=on`, `reference_scf.qe_clock` is an array of
`{name, wall_s, calls, rank, boundary}` on rank 0 only. CPU seconds are **not
measured**: the console says so, rather than labelling wall as CPU. Seconds are
rank 0's local elapsed wall, including waits inside each operation; counts
are invocations on that rank, not columns, pairs, iterations, or sums over
ranks. Existing Davidson max/min statistics remain unchanged in their original
location. No new barriers or timing collectives are introduced. The new fixed
counter scope is absent when the switch is off: no clock reads, TLS counter
access, hot-loop allocation, table allocation, serialization or console table
on that path (cached-flag/Option branches remain; literal zero instructions is
not claimed). The always-on pre-existing phase clocks are unchanged.

| Requested QE line | QE includes | BAANDI boundary and status |
|---|---|---|
| `init_run` | Setup of grids, tables, initial fields and starting-wavefunction rotation, after input/setup decisions | **Unavailable/inexact.** `baandi_setup` is driver entry through reference state, density-plan and mixer workspace setup; CLI preflight/UPF parsing and some basis work have different ownership. Do not relabel it `init_run`. |
| `electrons` | Entire electronic SCF procedure, including band driver, density, mixing, potentials, convergence and energy bookkeeping | **Unavailable/inexact.** `baandi_scf` spans first map entry through final map end, minus verification intervals. It includes direct energy audits and recovery work with no QE equivalent. It is not an `electrons` proxy. |
| `c_bands` | Band-driver setup, diagonalization/re-entry for every k, threshold/convergence handling | **Unavailable/inexact.** `baandi_solve` sums existing map diagonalization intervals (options plus solver); diagonal setup is in the potential interval and retries interleave density/audits. Calls count actual map evaluations, including retries/redo, not `c_bands` invocations. |
| `h_psi` | Complete kinetic, local-potential and nonlocal H application to one block | **Exact.** Whole distributed real operator application, including mirror exchanges, pair transforms/transposes, kinetic, projection reduction, coupling and back accumulation; one call per nonempty block. Includes atomic-start H as well as SCF H. |
| `vloc_psi` | Pair packing, inverse FFT, real-space local-potential multiply, forward FFT and unpack | **Exact, derived on the same rank.** Mirror/pack through mirror-back/unpack, subtracting only the separately measured kinetic loop inside that interval. All communication and FFTs remain inside. One call per nonempty H block. |
| `fftw` | Wavefunction FFT calls (both directions), also reached outside H, e.g. density/startup | **Unavailable/inexact.** Fused distributed transform instrumentation includes other work and batched calls have a different count unit. Neither sum of FFT kernels nor pair counts is emitted as `fftw`. Existing detailed FFT diagnostics remain available. |
| `calbec` (`h_psi:calbec` in the QE printout) | H's projector overlaps and their reduction | **Exact** for H-scoped QE label only. Local projection plus deterministic rank sum; one call per H block with projectors, even on a rank with zero local rows. An unqualified global QE `calbec` label may contain other work; the parser marks it inexact. |
| `add_vuspsi` | Nonlocal coupling of overlaps, back product and accumulation into H psi | **Exact.** Coupling, local projector back product and addition; one call per H block with projectors. Audit-only nonlocal energy calls are excluded. |
| `regterg` | One complete real Davidson solve, including starting H, borders, reduced problems, corrections, contractions and returned vectors | **Exact.** Reference solver's inclusive solve interval (including workspace/trace setup and finish rotation). Calls count solver invocations, including re-entry and first-step redo; no verification calls. |
| `rdiaghg` | Full generalized reduced eigensolve: overlap Cholesky, standard-form reduction, eigenproblem, back transform and root broadcast | **Unavailable/inexact.** Factor updates/pivot repair occur separately from the standard reduced eigensolve. `baandi_reduced_standard` is the latter interval, including transforms and broadcasts, with logical solve calls; it must not stand in for `rdiaghg`. |
| Davidson remainder: `regterg - h_psi - rdiaghg` | Only children nested within those Davidson invocations | **Unavailable/inexact.** The matching generalized child boundary is absent; moreover QE's printed global H and reduced clocks include startup. Never subtract global clocks from a parent. `baandi_davidson_remainder` subtracts solver-internal H and standard-reduced intervals from its own total, so it still includes factor work; denominator is Davidson invocations. |
| `sum_band` | Occupations/band sum and charge accumulation, transforms to reciprocal space and applicable symmetrization | **Unavailable/inexact.** `baandi_density` is real-space band accumulation and field-layout exchange only. Reciprocal transform and energy/occupation work are elsewhere; no sum of mixed phase buckets is labelled exact. |
| `v_of_rho` | Hartree/XC potential and associated energies, including final output-density evaluation | **Unavailable/inexact.** `baandi_input_potential` additionally assembles ionic potential, Hamiltonian diagonal and field redistribution, while output-density Hartree/XC is buried in energy audit and mixed-density fields in correction. It is not a complete potential clock. |
| `mix_rho` | Residual/error metric, history solve/update and density mixing | **Unavailable/inexact.** `baandi_mix` measures history update plus inverse transform/handoff; the stopping metric is evaluated separately before verification and on redo. It is not a QE mixing proxy. |
| `verification` | No counterpart | **BAANDI only.** All attempted independent checks, including failures; skipped checks have no call. Failed checks reused as later maps are charged exactly once here. |
| `output` | No counterpart with this boundary | **BAANDI only.** Materialization plus CLI output preparation/writes through the boundary just before final results serialization. Excludes final JSON write, console delivery and process/MPI exit. External delivered wall includes them. |

Unavailable QE rows have `wall_s: null, calls: null, boundary: "unavailable"`,
never fabricated zeroes. BAANDI-only measured alternatives have their own names.
Measured exact rows use `boundary: "exact"`. All H/solver rows cover setup plus
SCF, excluding verification snapshots on the **same rank**; verification is
reported separately, never removed from production delivered wall. Rows overlap
(parent/child), so adding them does not reconstruct total wall. The solver's
H interval includes its operator dispatch wrapper; operator `h_psi` excludes
that wrapper, which remains in the solver remainder.

The parser `scripts/qe_clock_compare.py` reads one full QE stdout and rank-zero
BAANDI results, emits seconds, calls and seconds/call side by side, preserves
missing rows and QE's printed precision, and rejects conflicting duplicate QE
rows (concatenated runs). It never computes speed ratios or fabricates a
remainder from global clocks. Scope/condition identity requires the manifest
review described in [the dev08 request](../performance/lanes/RSCF/qeclock-dev08-request.md).
Headline comparisons keep detailed profiling **off**; no performance default
or optimization is adopted by adding these clocks.
### P3-FFT2: pass ledger and row-fused experiment

Source base `8d1f3584eaa634f24b2bd2d11a1c1c341489919f`. The ticket's
**relayed**, single-profile motivation at `9e99bd3`, Si128 MPI8, rank max:
FFT forward + inverse 10.2 s, of which approximately 4.9 s is 1D kernels;
alltoall 0.86 + 0.98 s, middle pack/unpack 0.69 + 0.91 s, slab gather/scatter
0.53 + 0.70 s, transpose pack/unpack 0.47 + 0.40 s, multiply 0.56 s.
These maxima overlap in rank attribution and are not an additive wall model.
QE FFTW's supplied 6.7 s is motivation only, not a newly executed comparison.
P3-FFT direct packing and sparse clearing had no measurable dev08 gain;
M2_FFT=mkl was not adopted (small lengths lose). No backend default changes.

#### Every buffer traversal, one pair, one direction

Use the symbols above. `R(X,q)` / `W(X,q)` denote reading/writing q complex
values of buffer X. Rows for a real-space plane repeat p_r times; `P` is the
**same allocation as slab**, exposed by `plane()`, not a second buffer.
The potential is real (8 bytes/value), all other value streams are 16 bytes.
The MPI row counts the application-visible send read and receive write;
internal MPI copies, staging, cache-line effects and FFT scratch remain excluded.
Padding is sent/read even though the transform never consumes it.

| Direction / operation | Explicit streams on this rank, one pair |
| --- | --- |
| G -> r: initialize sticks | W(sticks,A_r) |
| coefficient scatter | R(coefficients,L_r), W(sticks,L_r) |
| first-axis inverse | R(sticks,A_r), W(sticks,A_r), plus backend scratch |
| transpose pack | R(sticks,A_r), W(send,A_r) |
| one inverse alltoall | R(send,R*b), W(recv,R*b) |
| initialize middle, per plane | W(middle,M) |
| middle unpack, per plane | R(recv,S), W(middle,S) |
| second-axis inverse, per plane | R(middle,M), W(middle,M), plus backend scratch |
| initialize slab, per plane | W(P,F) |
| slab scatter, per plane | R(middle,M), W(P,M) |
| last-axis inverse, per plane | R(P,F), W(P,F), plus backend scratch |
| potential multiply (once per round trip), per plane | R(P,F), R(V,F real values), W(P,F) |
| r -> G: last-axis forward, per plane | R(P,F), W(P,F), plus backend scratch |
| slab gather, per plane | R(P,M), W(middle,M) |
| second-axis forward, per plane | R(middle,M), W(middle,M), plus backend scratch |
| middle pack, per plane | R(middle,S), W(send,S) |
| one forward alltoall | R(send,R*b), W(recv,R*b) |
| transpose unpack | R(recv,A_r), W(sticks,A_r) |
| first-axis forward | R(sticks,A_r), W(sticks,A_r), plus backend scratch |
| gather + 1/N | R(sticks,L_r), W(coefficients,L_r) |

A **standalone transfer pass** below means a copy between two named work
buffers, separate from a transform; sparse accesses count as one pass over
the selected elements, not as full dense scans. Coefficient scatter/gather
is separately counted, so it cannot disappear through a naming change.

| Pass class | Current G -> r | Current r -> G | Ideal one-transpose, fused interfaces (each direction) |
| --- | --- | --- | --- |
| MPI collectives | 1 | 1 | 1; same useful payload, padding is a separate issue |
| Standalone transfers between FFT/MPI buffers | 3: sticks->send, recv->middle, middle->P | 3: P->middle, middle->send, recv->sticks | 0 standalone: transform output stores to send, recv unpack supplies plane transform input, adjacent local transforms share scratch |
| Separate initialization sweeps | 3: sticks, middle, P | 0 | 0 separate sweeps if every transform input is defined by fused value/zero fill or immutable zero holes |
| Coefficient input/output passes | 1 scatter | 1 gather/scale | 1 at each boundary, fused with first-axis load/store where supported |
| Transform stages | 3 | 3 | 3, same transform points; fusion does not remove arithmetic |
| Separate full-plane multiply sweep | 1 at turnaround | 0 | 0 standalone: V is read and multiplied while final-axis data are resident |

The ideal column is an **interface target**, not a claim of zero reads/writes
inside transforms. It still needs useful send writes A_r / p_r*S and useful
receive reads p_r*S / A_r in the two directions, MPI send reads + receive
writes of R*b each, all three axes' arithmetic/scratch, and F real-potential
reads per owned plane. Arbitrary sparse/rank-owned y positions cannot be
transformed in place in the present receive layout: source-rank-major padded
blocks are not contiguous y lanes. Achieving that ideal needs strided or
fused gather/scatter backend interfaces or a new communication layout; the
existing `BackendLanes` takes whole contiguous buffers. More band batching
cannot remove any collective here: M2d already batches the whole H block.

#### Implemented subset: `BAANDI_FFT_FUSED_ROWS=on` (default off)

Ordinary distributed `apply_local`, `apply_local_batch` and the direct-input
variant now optionally perform the last-axis inverse, V multiply and
forward together for **one row**. A construction-sized spectral input row
receives the active-z middle entries. RustFFT's immutable-input inverse
writes a second row, the unchanged scalar multiply updates it, the forward
runs in place, and only active z values are written back to middle. Rows are
independent: overwriting middle[k*n1+y] cannot affect a later y. The first
and second axes, their operation order, ORC input/output, final 1/N and both
collectives are unchanged. No tolerance, grid, stick partition or batch
width changes. No hand-written FFT or new dependency.

Inactive spectral-row entries remain exactly zero because the inverse never
modifies its input. Thus the entire repeated slab initialization disappears,
not just the overwritten stores removed by P3-FFT. The full slab is never
read or written by this apply route. Its allocation remains for the public
plane API, density and kept-plane entry points, which retain the old path.
The row pair plus RustFFT scratch is allocated only at construction. This
is a working-set/traffic change, **not a peak-allocation reduction**.

Exact abstract interface-byte saving per pair/rank is `16*p_r*F` with sparse
clear off, or `16*p_r*(F-M)` relative to sparse clear on. Other interface
byte counts do not fall: middle->row and row->middle still copy M elements
per plane, and inverse/multiply/forward still read/write F values, now in
reused O(n2) row storage instead of distinct rows of an O(n1*n2) slab.
No six-pass or DRAM-byte speedup is claimed from this cache-local schedule.
It targets the last-axis turnaround working set (slab scatter/gather and
multiply), not the larger but layout-constrained middle/MPI pack passes.
Immutable FFT scratch traffic and per-row dispatch can offset the benefit;
dev08 decides. For example, an F=36*36 slab formerly has 20,736 repeated
zero-store bytes per plane; a 72*36 slab has 41,472. These are derived counts,
not measured machine traffic or speed.

The switch is resolved once at plan creation through RunConfig, with normal
input key `fft_fused_rows`, manifest consultation and MPI agreement hashing.
With MKL lanes explicitly requested it records a `slab` fallback (immutable
row support is not exposed by that adapter), retaining MKL arithmetic.
An MKL **build** with the default rustfft M2 lanes uses the experiment.
Off selects the unchanged slab operations; middle helpers merely share the
existing unpack/FFT and FFT/pack code. The immutable and in-place RustFFT
interfaces do not promise identical rounding on every architecture/length,
so the candidate is conservatively **rounding class**, even where bits match.

The existing FFT profile retains all stage names and point/byte meanings.
Stages 6 through 10 are now called n1 times per plane, each for one row;
point totals stay equal, stage 6 bytes drop by the formula, and stages 7/9
flop estimates can differ by integer rounding of per-call estimates only.
Calls are instrumentation regions, not 3D transform counts. More timer reads
can perturb a profiled row arm; keep the outer H total/unattributed remainder
and report the effect as a single-run development observation. No added
barrier, collective, allocation or clock occurs with profiling disabled.

Derived from the own-words notes `fft-call-structure.md`, `qe-hpsi-and-fft.md`,
`qe-gamma-real-path.md` and this engine's code. A26/PFFT-3/P39A already fuse
at plane scale; this experiment narrows the final turnaround to one row and
uses immutable zero holes. P39A-4 resolved P39A-3's accuracy objection, not
a proof of speed for this change. A11/PFFT-4 per-lane overhead and P3-FFT's
neutral byte-removal screen remain counterevidence, explicitly retained.
The new premise is the supplied reference-driver movement profile and the
P3-FFT2 work order, not a repeat of lane batching, indexed transposes or MKL.

The [single-run dev08 request](../performance/lanes/RSCF/fft2-dev08-request.md)
contains validation scope, results and remote gates. Shared ledger/log,
README evidence and dashboard updates belong to the master integration;
this ticket edits only its owned accounting section and request.

## Per-map energy cost

P3-ENERGY is derived from the own-words `qe-scf-cycle.md`, sections 7–9,
not QE source. It does not repeat ENERGY-IDENTITY (which merely replaced a
nonlocal term and did not improve wall time): it defers the entire independent
output audit until needed. `BAANDI_REFSCF_LEAN_ENERGY=on`, or
`[performance] refscf_lean_energy=true`, is experimental and **default off**.

Here R is the owned real grid size, G the owned orbital rows, B the band count,
and P the projector count. A field evaluation consists of one forward and one
inverse distributed full-box FFT, a reciprocal Poisson pass, pointwise LDA and
real-grid integrals (plus rank reductions). These are logical costs, not timings.

| Computation | Cost per evaluation | Consumer | Lean policy |
| --- | --- | --- | --- |
| Input Hartree/XC fields | 2 FFTs, O(R) Poisson/XC/integral passes | Hamiltonian, diagonal mean, double counting, Harris–Foulkes | Reuse preceding mixed fields |
| Orbital density | Inverse FFT per occupied orbital pair, O(R B) accumulation | Stop dr2, mixing, delivered density | Unchanged |
| Reciprocal output/residual | 1 density forward transform, reciprocal metric pass | Stop, ethr/redo, Broyden | Unchanged |
| True orbital residuals | O(G B) from retained HX, band reduction | Verification residual guard, records | Unchanged |
| Kinetic energy | O(G B) pass | Direct audit, records, delivery | Deferred |
| Nonlocal energy | Second projection O(G P B), coupling, reductions | Direct audit, records, delivery | Deferred |
| Ionic local and zero-mode terms | O(R) integrals | Direct audit, records, delivery | Deferred |
| Band sum | O(B) | Variational energy, Harris–Foulkes, records | Retained |
| Output density times input HXC | O(R) integral and reduction | Variational double counting | Retained |
| Output Hartree/XC | 2 FFTs and field passes | Direct/band-form audit, energy verification, delivery | Deferred |
| Functional-difference bound | 2 O(R) residual/input-potential integrals and O(B) residual sum | Audit records only | Deferred |
| Mixed-density fields | 2 FFTs and field passes | descf, variational Hartree/XC, next Hamiltonian | Retained and cached |
| descf | O(R) integral and reduction | Variational trace energy | Unchanged |
| Independent verification | Fresh solve and full map audit | Coulomb or ADR 0044 energy verdict | Full audit |

No direct audit energy enters the stopping estimator, ethr ratchet, redo, or
Broyden update. Candidate maps at an estimator stop are audited before any
verification attempt, including failed attempts. Predictor-skipped candidates
remain lean unless they must be delivered at the iteration limit. All verification maps
are fully audited. Maximum-iteration and mixer-failure deliveries are audited
without another solve or density construction. Thus every delivered ScfResult
has a complete finite energy audit; deferred work is charged to map/energy time.

Intermediate records retain their existing keys. In lean mode the unavailable
direct total, kinetic, ionic local/zero mode, nonlocal, output Hartree/XC,
band-form and functional-difference-bound slots are JSON **null** (not zero,
not the mixed-density variational energy). The shared in-memory EnergyAudit
uses NaN for these unavailable slots; no numerical decision reads them. A
history energy difference involving an unavailable direct total is null too.
Band sum, double counting, Ewald, Harris–Foulkes and the separately labelled
mixed-density variational energy remain available. Trace uses the latter.

A call-owned, one-entry exact-density cache is active in both arms. It stores
only owned density/field planes, never gathers a grid. All ranks agree on a
cache hit before skipping collectives; density bits, including signed zero,
must match. Functional, cell and layout are fixed for the cache lifetime.
This removes mixed-fields → next-input and candidate-output → verification-input
duplicates. A miss replaces the entry; nonconsecutive duplicates can still miss.
The cache adds one scalar agreement reduction and owned-plane comparisons/copies;
its net speed is for dev08 to measure. No wall-time gain is claimed on Mac.

For an ordinary continuing map without retries, the old path had six field
FFTs (input, output and mixed density pairs). Exact mixed-input reuse reduces
this to four with lean off; lean on reduces it to two after the first map.
These counts exclude orbital-density transforms, the stopping transform,
verification and the final audit. They are a work model, not measured speed.

### P3-ENERGY Mac validation

The synthetic `lean_energy_preserves_delivery_and_counts` test compares 12
configurations (cubic, skewed, elongated; verification off, Coulomb verification,
energy verification, and one-iteration exhaustion), each on two simulated ranks.
Delivered total energy difference is exactly 0 Ha and density difference exactly
0 electron/Bohr^3; mixed variational energy bits, dr2, first stop, verification
call count, SCF count, Davidson iterations and H columns agree. The three-rank
field-cache test compares fresh/cached Hartree, XC and potential integrals by
bits and both potential arrays, including a miss on just one rank. CLI tests
check that unavailable direct audit terms serialize as null in Ha and Ry.

Mac gates: refscf 29 passed / 1 ignored (external Si matrix, dev08 only); CLI
reference and numerics switches pass; cache test, fmt and clippy pass. The
additional `scripts/check.sh fast dft-engine dft-cli dft-numerics` reports
build/fmt/clippy PASS, test FAIL: 552 passed, 3 failed, 148 ignored. All three
failures are existing CLI pseudo-cache tests rejecting the sandbox scratch path
inside the checkout (`pseudopotential cache must be outside repositories`). No
cache policy or test tolerance was changed. This is synthetic correctness
validation, not real-material/dev08 performance acceptance. Request:
`docs/performance/lanes/RSCF/energy-dev08-request.md`.
## RMM-Davidson

P3-RMM adds a default-off, standalone `dft_solvers::ref_rmm` module, derived
only from the own-words [RMM note](../research/qe-study/qe-rmm-davidson.md).
`[performance] refscf_rmm="on"` (`BAANDI_REFSCF_RMM=on`) selects Davidson
at SCF iterations 1, 2 and multiples of five, and RMM on other iterations.
The first-step redo retains index 1. Verification always calls Davidson;
failed-verification reuse, when enabled, retains the already solved Davidson
map and charges its work once. Existing ethr, stopping, energy and mixing
rules are unchanged. There is no convergence retry or Davidson fallback
inside an RMM map.

Each call applies H to the n starting orbitals, measures both X^T X and
X^T H X and performs the generalized Rayleigh-Ritz rotation through the
existing reference Davidson workspace and reduced-EVD adapters. Both X
and H X are rotated and retained. One sweep follows, at most four steps
per band, with frozen bands excluded from subsequent H batches. Every
history entry stores its own normalized state, H image and own-Rayleigh
residual. DIIS minimizes residual norm with `v^T B v=1`; the canonical
metric EVD keeps beta>1e-14 and falls back to the newest entry when at most
one direction survives. Reduced problems are replicated on all ranks.

The preconditioner uses the bare local-row kinetic energy, not Davidson's
full Hamiltonian diagonal. With T=sum(t_G*x_G^2) in Ha and y=t_G/(1.5T),
`K=-(4/3)/T*(27+18y+12y²+8y³)/(27+18y+12y²+8y³+16y⁴)`.
For y>1 the reciprocal polynomial avoids overflow. At T=0 the limit is
-1/t_G on nonzero kinetic rows and zero on G=0 (whose normalized constant
state has zero residual component). ORC already encodes Gamma weights;
no extra pair/G=0 factor is introduced.

The normalized update uses the parabola through E(0), its exact slope and
E(0.5), with step clamped to [0.05,1]. The curvature zero test is 5e-17 Ha
(the note's 1e-16 Ry). Each band freezes when its eigenvalue change is
strictly below ethr (occupied), or max(5*ethr,5e-6 Ha) (empty). Freeze flags
only save work; an unconverged sweep is accepted by the SCF loop.

Final band-order block Gram-Schmidt carries H X through every operation.
A panel overlap and the existing guarded inverse-Cholesky adapter implement
the positive-diagonal QR equivalent of classical GS, avoiding per-band
MPI reductions in accepted panels. A rejected panel uses literal one-pass
classical GS with measured norms; norm<1e-16 fails explicitly, without a
solver fallback. Finished blocks project all later bands with one block
overlap and two GEMMs. The minimum pre-normalization norm is recorded.
Final eigenvalues are diagonal Rayleigh quotients; eigenvalues, X and H X
are sorted together, without a final Ritz rotation or fresh H application.

Exactly `n + sum_k u_k <= 5n` H columns are spent at default depth four.
Scratch and histories are reserved before the sweep. Matrix products use
the existing GEMM workspaces/routes and small EVD backend; no library,
logging, file I/O or retained cross-SCF history is introduced. Options
`refscf_rmm_steps` (default 4, 1–16) and `refscf_rmm_block` (default 16, 1–64)
are explicit diagnostic parameters, not material-dependent rules.

Full diagnostics retain per-SCF-iteration RMM records with exclusive stage
seconds, per-stage GEMM FLOPs, H columns, per-step active counts, DIIS
fallback/clamp counts and minimum GS norm. GEMM totals count nominal 2mnk
only, excluding scalar band work and backend EVD internals. Operator H
profiling uses the same scope as Davidson. Standard diagnostics retain no
RMM stage records and execute identical numerical arithmetic. Existing
Davidson stage tables exclude RMM; add the separate RMM totals for run
work. Legacy iteration-count fields retain their names and count sweep
steps on records tagged `solver="rmm"`. They are not Davidson iterations
or a residual certificate.

Local evidence, limitations and the one-run MPI8 request are in
[the P3-RMM request](../performance/lanes/RSCF/rmm-dev08-request.md).

The strict local fixtures (common 1e-12 Ry estimator) meet the 1e-10 Ha
D/M energy gate. At the default 1e-8 Ry estimator all four synthetic
fixtures converge and verify but fail that stricter comparison (maximum
6.115339346e-10 Ha). This is a recorded acceptance blocker, not a relaxed
test or a basis for default adoption; explicit ignored failure reproducers
remain in `tests/refscf.rs` and the request gives their command and results.

## Cutoff ladder

P3-LADDER adds the run-scoped `BAANDI_REFSCF_LADDER=on|off` switch
(`[performance] refscf_ladder="on"`), **default off**. The legacy ladder
switch and its behavior are unchanged. The new switch acts only at the
reference driver's entry, including `run_with` and simulated collectives.
It does not recursively invoke the ladder or construct a legacy session.

The coarse cutoff is **0.5 times the requested cutoff**, C1's pre-registered
fraction, with an automatically generated smaller density grid. The fine
input grid override, if present, remains the fine grid. Coarse Broyden beta,
history length, preconditioning, Davidson rules and stopping estimator are
those of the reference driver. Its estimator threshold is **100 times the
fine threshold in Ry**: a fixed two-order separation, not fitted to any
material or cell. Its budget is `min(30, max_iterations)`, following the C1
cost bound. The coarse stage has verification off because it supplies only
a start. On budget exhaustion it exports its last output and records
`coarse_stopped=false`; numerical failures are errors, not silent fallback.
There is no change to the fine iteration budget. An insufficient coarse
basis (for example fewer plane waves than requested bands) fails explicitly;
no fraction search or system-specific cutoff adjustment is implemented.

The coarse output density and orbitals belong to the same last map. The
owned real Gamma rows are injected by Miller index, using the shared C1
mapping helper, and redistributed to their fine FFT-stick owners. Row
values are copied with no coefficient arithmetic, normalization or QR.
Both members of each Gamma pair keep their orientation and their real/imaginary
row weights. New rows are zero. Thus overlaps and Hermitian symmetry are
preserved (overlap sums themselves can differ by floating-point association).
The density uses the shared C1 Fourier-padding convention: forward-normalized
coefficients, signed frequencies, and removal of even-grid Nyquist planes,
which lie outside the coarse density sphere. Fine inverse FFT roundoff is
not clipped or charge-rescaled. LDA on the initial prolonged density uses
the same absolute-density rule as the reference map. Transfer is one
padded equal-block all-to-all for orbitals and one for density coefficients;
no full orbital or density value array is gathered. Rank-local ownership
can change with grid and rank count, including empty owners.

Fine setup builds its own Hamiltonian, ionic potential and fine-grid fields,
but skips the atomic orbital solve. It consumes the prolonged density and
owned orbitals and constructs an empty Broyden history. S2 secant seeding is
not enabled. The first fine call, first-step redo, ethr reset/ratchet,
stopping estimator, fine energy audit and verification/recovery logic remain
unchanged. The public result and verification are exclusively fine-stage
outputs. `verify="off"` still returns the fine estimator stop without a
passing certificate. Fine atomic-start H applications are recorded as zero
when the orbitals were prolonged.

Every record has `stage="fine"` or `stage="coarse"`, with stage-local
iteration numbers. To preserve existing consumers, `Run.records`,
`reference_scf.iterations`, `convergence.history` and verification attempts
remain **fine only**. Coarse rows are in
`reference_scf.cutoff_ladder.iterations`; that block also reports fraction,
coarse threshold, convergence flag, setup phases, coarse seconds, coarse/fine
map counts and separate solve-call counts (redo/retry calls are not extra
outer maps). No coarse verification or public-result assembly is performed.
Davidson stage totals include both stages. Coarse setup/maps/export are in
the run's wall; fine setup includes the redistribution/inverse transform.
All map timestamps share the original driver-entry origin. `wall_s.start`
is inclusive through the first fine map and contains the separately exposed
`coarse_in_start` subtotal; do not add that subtotal twice. With the ladder
enabled, `peak_bytes_per_rank` describes the fine-stage
arrays only and explicitly excludes the coarse/transfer workspace; use the
whole-process RSS measurement for a whole-run peak. The QE-clock
proxy table remains the fine stage's; coarse phase evidence is in the ladder
block. `refscf-check --stop-to-verified` accepts historical rows without a
stage and current fine rows, and rejects coarse rows in the fine history.
Its interval is still first **fine** estimator stop to passing verification.

Local tests and the one-run-per-arm dev08 request are recorded in
[ladder-dev08-request.md](../performance/lanes/RSCF/ladder-dev08-request.md).
No performance promotion follows from Mac correctness or map counts.

## Kept HX warm start

P3-KEPTHX adds `BAANDI_REFSCF_KEPT_HX=on` (or `[performance]
refscf_kept_hx="on"`), default **off**. An ordinary reference Davidson map
retains its exact returned Ritz block X, HX, and the owned input local
potential in the Hamiltonian's plane layout. The next eligible map supplies
`HX_old + (V_new - V_old) X` as its starting product. The fixed kinetic and
nonlocal operators cancel in this difference. Generalized overlap, rank,
Ritz-change convergence and independent fixed-point verification are unchanged.
The solver's supplied-product entry validates shape and finite values and
rejects disagreement about presence across ranks in its existing input check.

The cache holds a unique immutable operator-instance token (basis, kinetic,
projectors, row and plane layout), a bitwise copy of X, and a monotonic owned
potential epoch. It compares the entire retained X, not its address or a hash.
The actual old/new owned potentials form the delta; epochs alone never stand
in for potential data. One fixed-tree collective agrees a cache hit before any
rank enters the delta kernel. An entry is consumed on lookup, including on a
miss. The first map, first-iteration redo, bounded eigensolver retries,
independent verification maps, and rebuilt cutoff/layout stages use fresh H.
A new fine-stage cache is constructed at a ladder handoff. RMM maps do not
consume or publish kept products. Verification results reused as ordinary
maps carry no retained HX cache into subsequent entries.

`band_density_sum_keeping` retains leading complete real-space pair slabs
while constructing the usual density. Storage is capped at **512 MiB per rank**;
a replicated upper bound on plane ownership selects the same retained pair
count everywhere. `apply_local_change` reuses these slabs and transforms any
remaining pairs, including an odd final band. Thus a B-band hit replaces
exactly B full-H starting columns regardless of the slab cap. Retained X, HX
and owned potential storage are additional to the slab limit. The existing
kernel requires full-grid indexing for a transient delta buffer; only owned
planes are populated and read, and other slots are NaN.

`results.json.reference_scf.kept_hx` records `cache_hits`,
`full_h_starting_columns_replaced`, `delta_v_columns`, and
`retained_slab_bytes_per_rank`. Counts include the coarse and fine stages;
slab bytes are the maximum retained payload on each rank across those stages,
not their sum and not RSS. Counts are available with lightweight diagnostics.
No dev08 speedup or default adoption is claimed by the Mac correctness tests.

Local screening command:
`cargo test -p dft-engine --test refscf fast_order_ -- --nocapture`. The nonlocal two-rank fixture compares the delta-updated product with
a direct H application at relative Frobenius error <= 1e-13, checks a B-column
saving, and invalidates the orbital block on one rank. The two-rank skew-cell
ladder screen checks bitwise raw injection and zero fill, Gamma pairs, charge
error <= 1e-12 e, orbital Gram change <= 1e-12, and verified fine/direct energy
agreement <= 1e-10 Ha at a tight threshold. These are synthetic internal
oracles, not QE references. The mechanism follows the already adopted older
MPI starting-product reuse (ledger M2-STICKS-START-REUSE-DEV08-ADOPT); this new
reference path still requires its own dev08 common-input and full-convergence
measurements before promotion.

Local Mac screen (2026-11-25, synthetic fixtures, the command above): **2 passed,
0 failed, 0 ignored**. B = 3: relative product error
4.39708409238569948e-17; fresh/hit/invalidated starting full-H columns = 3/0/3;
one hit, three replaced columns, three delta-V columns; retained slab payload
[4096, 4096] bytes. The skew-cell injection preserves copied and zero-filled
ORC bits without arithmetic in the test gather, preserves Gamma conjugacy,
has charge error 6.21724893790087663e-15 e and maximum Gram-entry change 0.
Verified ladder versus direct fresh-H fine delivery differs by
1.26565424807267846e-14 Ha; continuing the exact inspected seed through the
same fine-stage driver differs by 1.48769885299770976e-14 Ha. The direct fine
energy is -1.10503701509437602 Ha. Independent verification remains enabled.
The screen also rejects a rebuilt equal-shaped operator, single-rank cache
invalidation, rank-deficient starts, nonfinite supplied HX on one rank, and
rank disagreement about supplied HX. These numbers establish correctness only;
no performance adoption or external-code comparison follows from them.


### P3-KEPTHX2: density-pass accounting and storage lifetime

The pre-change kept path already did **ceil(B/2) inverse pair transforms and
zero forward transforms**, exactly like the plain density sum. Both paths
pack the same ORC columns, exchange the same mirror rows and do the same
batched inverse transpose. No second FFT was responsible for the work-order
density regression. At each pair/plane the old kept path made three scans of
the inverse plane buffer: a complex slab copy, a real-band density update,
and an imaginary-band density update. The plain path made the last two scans.
In addition, `KeptHx::density` allocated and zero-initialized K*P complex values
on every map, and `prepare` dropped them after each lookup. K is the retained
complete-pair count and P is the rank's owned plane-point count. For a retained
pair, the extra copy read and wrote 16*P bytes each; allocation also initialized
16*P bytes. These are logical traffic counts, not measured DRAM traffic.

The kept branch now loads each inverse-buffer complex value once, stores it
in the retained slab and applies the two density updates at that point in
ascending band order, retaining `weight * value * value / volume` exactly.
It adds no separate full-buffer copy pass. The non-kept branch is unchanged.
Slabs remain contiguous `[pair][owned plane][row][z]`, complex re/im for bands
2t/2t+1, unscaled. The pair stride is P complex values, the plane stride is
rows*nz, and z has unit stride. `apply_local_change` reads this same layout
without repacking; no change to its transforms or arithmetic was needed.
Zero-occupation buffer bands are stored too: the next starting HX requires
them. An odd tail and pairs beyond the cap still take the existing fallback.

Slab allocation belongs to the stage-owned `KeptHx`, separately from entry
eligibility. The first density call sizes it; subsequent equal-shape maps
reuse its length/capacity with no zeroing or allocation. Invalidation and a
consumed or missed lookup discard only the entry, not this storage; a new
entry is published only after density has overwritten the retained range.
A new stage owns a new buffer, so operator changes and ladder handoffs cannot
reuse stale values. Retained payload and the 512 MiB cap are unchanged; a lower
RSS is not promised. Existing X/HX/potential snapshot allocations and the delta
workspace remain; they are outside the removed slab allocation/pass. The
reference map also allocates its plane-density output and redistributes it
back to the field layout in both arms; neither operation changes here.

Tests compare plain/kept inverse and forward counters, every density/slab bit,
and delta-HX bits against slabs constructed by the unchanged `band_values`
path (1/2/5 simulated ranks, odd/even bands, zero occupation, no/partial/full
retention). A private cache test poisons the buffer and checks its pointer,
capacity, complete overwrite and density bits across hits and invalidation.
The delivery fingerprint test includes a verified two-rank skew/nonlocal SCF
and repeatability. Mac evidence is correctness only; the exact four-run dev08
request and evidence limits are in
[the P3-KEPTHX2 request](../performance/lanes/RSCF/kepthx2-dev08-request.md).


### P4-DENSKEEP: optional slab recomputation

`BAANDI_REFSCF_KEPT_RECOMPUTE=on` (`[performance]
refscf_kept_recompute="on"`), default **off**, sets the retained pair count to
zero for the lifetime of a new `KeptHx` stage. It requires kept HX to be
enabled. X, HX, the old potential, identity/epoch checks, invalidation and
verification preservation are unchanged. The existing f64 zero-kept-pair
fallback reconstructs the inverse planes only when applying the local change;
no compressed arithmetic or new transform implementation is introduced.
Density takes the existing no-store traversal and allocates no slab. The
reported retained slab high water is zero (X/HX/potential are still retained).

For B bands and K previously retained complete pairs, density still executes
ceil(B/2) inverse pairs and zero forward pairs. Each warm hit still replaces B
full-H columns, but now executes K more inverse pairs, plus the existing
mirror-row exchange and inverse transpose when the stored arm needed none.
The forward pair count is unchanged. Thus the memory saving is exact; retention
of the previously observed solve-time gain is **not established** without a
dev08 measurement. No material-dependent cutoff or auto-selection is installed.

The storage is already contiguous `[pair][owned plane][row][z]`: 16-byte
complex stride, 16*P-byte pair stride. All K*P values are consumed by the
local-change kernel, including zero-occupation buffer bands; there is no
unused local subregion to omit generally. Storage costs 16*K*P bytes written
per ordinary map and the same payload initialized once per stage. Mac replay
and process page-fault observations confirm first touch rather than recurring
faults; they do not measure hardware cache misses or DRAM transactions.

`fast_order_denskeep_*` checks 1/2/9 simulated ranks (including zero-plane
ranks), 3/4 bands, bit-identical starting HX and density versus stored f64,
and <=1e-13 relative Frobenius error versus direct H. Four two-rank synthetic
full SCFs (cubic, skew/nonlocal, elongated, slab) have exactly zero energy
difference, identical density bits and Davidson/H-column counts, and pass
independent verification. The strict 1e-10 Ha difference gate is unchanged.
Two ignored release diagnostics isolate actual common-input density work and
large memory-only traversal/first-touch work; commands, numbers, rejected
alternatives and the one-run-per-arm remote screen are in the
[P4-DENSKEEP request](../performance/lanes/RSCF/p4-denskeep-dev08-request.md).

## Warm verification maps

P3-VERWARM adds `BAANDI_REFSCF_VERIFY_WARM=on|off` (input
`[performance] refscf_verify_warm="on"`), **default off**. It requires
`REFSCF_KEPT_HX=on` (the adopted default); disabling the latter gives a cold
fallback. This section supersedes the verification-boundary cache exclusion
in the historical P3-KEPTHX section only when the new switch is on.

### Independence and state contract

The check remains a new evaluation of `F(y)` at `y = rho_out`, with a newly
installed `H[y]`, a new Davidson solve and a newly accumulated density. The
ordinary candidate's density is never substituted for that result. Both arms
already start Davidson from the candidate's returned Ritz block X. W supplies
`H_old X + (V[y] - V_old) X` instead of recomputing the initial full-H product.
Kinetic energy and fixed nonlocal projectors cancel exactly in this identity;
its floating-point evaluation is rounding class. The operator-instance token,
bitwise X comparison, finite-value checks, potential epoch and collective
eligibility agreement remain active. Every correction is applied with the
current operator. The generalized overlap/Ritz criteria, retry policy and
all residual diagnostics are unchanged.

Independence here means a separate fixed-point map, not a statistically
independent starting guess. Supplying an equivalent starting product does
not turn the check into a replay of the old density. It also does not add a
global uniqueness or spectral-completeness guarantee absent from the cold
check. The check still measures `E_H[F(y)-y] < 5e-9 Ha` under the default
profile, with the existing residual-guard policy. The energy criterion stays
off under user decision C. No pass is inferred from a cache hit.

The verification tolerance remains **exactly**
`next_ethr_ry(ethr, dr2, electrons) / 2` Ha, including the existing per-band
thresholds and retry semantics. The starting block is already warm in D;
the fresh full-H starting application was an implementation choice. The
threshold rule is a contract, not a cost knob changed by this experiment.

During a warm check the cache retains the candidate X, original HX, original
potential and its bounded slabs. The check clones only the supplied HX for
the delta update; its density uses the ordinary accumulation kernel and
neither reads cached density nor overwrites the candidate slabs. On failure,
Broyden receives the same ordinary pair and the next ordinary map starts from
the same candidate X. It can then form HX at the mixed recovery potential
from the preserved **original** potential, not from the verification
potential. No check-state orbitals or products enter delivery. There is no
second slab cache. The extra HX copy is covered by the existing conservative
retained-X/HX + supplied-HX + delta-workspace envelope; RSS still needs dev08
measurement.

Retries invalidate the entry as before, including during verification. RMM,
changed vectors/operator/layout, first-step redo and cutoff handoff retain
their cold fallbacks. Optional VERIFY_REUSE still consumes the failed map
according to its existing contract; a preserved entry cannot hit its changed
Ritz block and falls back collectively. No default for that rejected variant
changes here.

`Record.warm_start` records an actual supplied-product hit in an ordinary
map (false for reuse of a whole map). `Record.verification_warm` is absent on skipped or
nonattempt maps and otherwise says whether any call in that attempt was
warm (retry calls can be cold). Reference verification JSON adds per-attempt
`warm_start` and `h_columns`; the latter equals the existing `h_vectors`,
including retries, and excludes delta-V columns. Aggregate `kept_hx` counts
include verification hits. The stop-to-verified reader carries the optional
flags as `verification_warm_starts`; old files remain readable and timing
closure/decisions are unchanged.

### Work accounting and scope

A B-band hit removes B starting full-H columns: kinetic and nonlocal work,
plus their starting full pair transforms. The delta still costs B local-V
columns, normally `floor(B/2)` forward pair transforms using the candidate's
retained real-space slabs; uncached pairs and an odd tail need inverses too.
All correction columns, reduced eigensolves, new-density inverse transforms,
potential construction, Coulomb metric and energy/residual audit remain.
Thus H-column reduction is not a claim that those columns cost zero.

For Si128's archived B=260, two eligible verification attempts replace **520
full-H columns** with 520 delta-V columns. One first recovery after the failed
attempt can replace a further **260 full-H columns**. Subsequent ordinary
maps, including a gated skip, already use ordinary kept HX in D. These are
conditional work estimates, not measured dev08 counts or a promised wall
saving; retries/misses or a changed stopping trajectory must be reported.

The task-relayed ac6aee9c MPI8 single run has Si128 verified/verify-off wall
30.15/23.59 s, a **6.56 s** stop-to-delivery overhead containing attempts and
recovery, not 6.56 s spent only in verification solves. It has 17 ordinary
maps, two attempts and one gated skip. Si64's corresponding 6.40/4.57 s gap
is **1.83 s**. These motivating observations are not reproduced here and
are not new QE references. This change does not revisit the rejected RMM
trajectory or failed-map reuse: it removes initial full-H work while retaining
the same candidate state, mixer and eigensolve rules.

Mac correctness evidence and the one-run-per-arm remote request are recorded
in [verwarm-dev08-request.md](../performance/lanes/RSCF/verwarm-dev08-request.md).
The Mac fixtures exercise cubic, skew, elongated, synthetic nonlocal and
actual failed-verification recovery on 1/3/8 simulated ranks. They are small
internal oracles, not a physical Si128 SCF or MPI/MKL speed evidence.

## NLCC

**Vacuum densities.** Accepted [ADR 0046](../decisions/0046-xc-small-and-negative-density.md)
adopts default-on `XC_SCREEN_NEGATIVE`: XC alone evaluates an auxiliary field
that is zero for combined channel density <=1e-10 Bohr^-3, with a per-channel
screened-charge budget of 1e-4 of total valence electrons. Physical valence and
core arrays remain signed and unchanged. PBE differentiates this auxiliary
field and masks the complete potential after divergence; the derivative is
valid away from mask crossings, not at the hard cutoff. Start, map and fresh
verification use the same model. Screening counts, integrated absolute charge
and minimum/global grid index are diagnostics, not an energy-error certificate.
`XC_SCREEN_NEGATIVE=off` retains historical rejection/continuation behavior.
The default is supported by the dev08 gen14 PBE Gamma result (7 iterations,
VERIFIED, <=1e-7 electrons screened); the inventory and qualification limits
are in the ADR.

G1-NLCC extends the implemented Start/Driver envelope above to norm-conserving
LDA-PZ and modified LDA-PW datasets with nonlinear core correction. Earlier
statements in those sections that NLCC rejects are superseded by this section.
G2B extends this support to unpolarized PBE, including the total-density
gradient (see PBE below). Spin, smearing, symmetry, density floors and
force/stress delivery remain outside these stages. Functional mismatch still
rejects unless the caller explicitly requests the existing diagnostic override;
using that override does not qualify a PBE dataset for physical LDA work.

### Core construction and ownership

`StartState.core_density` is an optional immutable field in electrons/Bohr^3,
ordered exactly like `density` on the owned real-space planes. Each species'
`PP_NLCC` stores the physical radial density, not `4 pi r^2 rho`. The existing
`spherical_density_transform_with_rule` integrates it with the selected radial
quadrature. Its coefficients are divided by cell volume and multiplied by the
same atomic structure factors as the ionic local field. Both setup variants
reuse the density-sphere ownership, shell IDs/cache and structure factors already
built for atomic charge and the ionic potential. Species without a core add
nothing. One additional owned inverse transform constructs the core; no full
core grid or reciprocal core array is gathered. LDA discards the reciprocal
core; PBE retains an owned Fourier cache as described below. The legacy global
reality/finite-value guard is preserved using scalar rank maxima. Core ringing
is neither clipped nor normalized, and G=0 is the radial integral, never the
valence electron count. A prolonged ladder seed carries only valence density;
the core is rebuilt on the target grid from the same immutable atomic data.

The LDA core adds `8 F_r` retained bytes and transient owned reciprocal/transform
storage. PBE additionally retains `32 G_r + 8 F_r` bytes for the core Fourier
cache and two combined-density scratch arrays. `StartMemory` includes those
named buffers. With no NLCC there is no
additional field allocation, transform or collective; charge, ionic and XC
arithmetic and the existing energy reduction retain their previous order.

### XC and energy accounting

For signed mixed valence density `n` and fixed core `c`, the reference LDA map
and initial potential evaluate `E_xc[|n+c|]` and `v_xc[|n+c|]`. The existing
reference absolute-density rule is applied after adding the core. On positive
total density this is the legacy `evaluate_with_core` functional exactly.
Negative-total mixed inputs retain the reference signed-density policy; they
are not asserted to have the legacy driver's negative-density behavior.

Only XC sees `c`. Hartree uses `n`, and electron count, occupations, density
output, G=0 charge, Broyden history, residual `dr2`, verification charge and
`descf` density differences remain valence-only. In particular:

- Direct energy: `T + E_local + E_nonlocal + E_H[n_out] + E_xc[n_out+c] + E_ion`.
- Double counting: `integral n_out (V_H[n_in] + v_xc[n_in+c])`.
- Harris–Foulkes: `E_band - E_H[n_in] + E_xc[n_in+c]
  - integral n_in v_xc[n_in+c] + E_ion`.
- The mixed-density variational energy uses the same XC replacement; the
  signed valence integral defines `descf`. No core-potential integral is
  subtracted a second time or added to the total energy.

The delivered `EnergyAudit.nlcc` matches the legacy output-density ledger:
core electrons, `integral c v_xc[n_out+c]` and
`integral n_out v_xc[n_out+c]`. Ledger integrals share the existing field-energy
reduction. The exact-density cache is valid because core/cell/functional/layout
are immutable for its lifetime. Lean maps defer the output ledger together
with the output audit (`nlcc=None` means unavailable there); stopped, verified,
iteration-limit and failed-mixer deliveries complete it through the existing
full-audit boundary. Verification has the same core as every ordinary map.

### Validation and limits

`nlcc_owned_start_matches_legacy_and_gaussian_fourier` checks the analytic
unit-charge Gaussian Fourier coefficients, legacy core/start density and
Hartree-plus-XC initial potential for LDA-PZ/PW, a skew two-species cell with
one core-free species, both local-shell setup modes and 1/3/8 simulated ranks.
The map unit test checks valence-only Hartree, total-density XC and both ledger
integrals, including signed mixed input. No-core Si has a before/after
fingerprint test over delivered arrays, energy and scientific decision traces.

Full delivered-state matrices are ignored in debug (ADR 0004). They compare
local cubic LDA-PZ and nonlocal skew LDA-PW to tight, certified legacy SCF,
at standard `1e-8 Ry` and tight `1e-12 Ry` estimator thresholds. Each layout is
repeated at 1/3/8 ranks. The existing `common/rank_rule.rs` checks the complete
ordinary/verification trajectory and demands threshold-margin evidence for
any first divergent discrete predicate. Energy caps are `5e-9 Ha` standard
and `1e-10 Ha` tight; charge is within `1e-10 e`. Density and occupied-root
caps are `1e-6 e/Bohr^3` (tight: `1e-7`) and `1e-5 Ha` against legacy; the independent cross-rank
density cap remains `1e-9 e/Bohr^3` unless the rank rule explains a branch.
All returned roots are checked within `1e-5 Ha`; empty roots retain the existing
loose empty-band solver threshold.

```sh
cargo test --release -p dft-engine --test refscf nlcc_delivered -- --ignored --nocapture --test-threads=1
DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf nlcc_no_core_silicon_fingerprint -- --ignored --nocapture
```

No committed real LDA NLCC UPF is available. Synthetic fixtures are not QE or
material-validation evidence. The real MPI8 and optional correctly matched
LDA dataset protocol is [the G1 dev08 request](../performance/lanes/RSCF/g1-nlcc-dev08-request.md).

#### G1 local evidence

Mac Apple Silicon, pure Rust, base `d0a93375b3b9f91ce79218e5e8dfcf7c15d9eb03`.
These are engine-to-engine synthetic comparisons, not QE references or timings.
The release matrix covers 12 system/threshold/rank combinations, each twice;
no unexplained or rounding-threshold trajectory branch occurred.

| Fixture | Estimator (Ry) | Maximum legacy energy difference (Ha) | Maximum density difference (e/Bohr^3) | Maximum occupied-root difference (Ha) |
| --- | ---: | ---: | ---: | ---: |
| Cubic local LDA-PZ | 1e-8 | 7.3443e-11 | 7.6992e-7 | 2.9851e-6 |
| Cubic local LDA-PZ | 1e-12 | 1.7764e-15 | 4.1340e-9 | 1.6235e-8 |
| Skew nonlocal LDA-PW | 1e-8 | 5.2606e-10 | 6.3892e-7 | 2.5841e-8 |
| Skew nonlocal LDA-PW | 1e-12 | 1.2413e-13 | 9.5408e-9 | 9.6301e-10 |

The initial MPI-harness probe accidentally used local **skew** LDA-PZ instead
of the matrix's local cubic fixture. At `1e-8 Ry` it passed the energy class
(`2.3783e-10 Ha`) and verification, but missed the fixed synthetic density
cap (`1.5916e-6 e/Bohr^3 > 1e-6`); occupied-root difference was `6.3477e-6 Ha`.
This is not hidden as a pass and no tolerance was relaxed. The harness now
uses the exact two matrix fixtures. Energy-class acceptance alone does not
promise the auxiliary pointwise density cap on additional systems.

Initial test compilation caught a non-Sync legacy problem captured by a
simulator closure; the test now captures only immutable arrays. The initial
start comparison also compared atomic start to the legacy default uniform
start; selecting `StartDensity::Atomic` fixed the oracle setup. Neither fix
changed a numerical tolerance or production algorithm.

No-core Si2 regression: the same release test was run with `refstart.rs` and
`refscf/map.rs` from the base commit, then with the G1 implementation. All
three FNV-1a fingerprints of the full scientific result/trajectory agree
exactly (clocks excluded): rank 1 `4d94c48aa4101f65`, rank 3
`dd3d0b80dddf9504`, rank 8 `9c03585f1f4dc0ea`. Respective delivered energies
are `-7.23375645269802980`, `-7.23375645270790546`, and
`-7.23375645263328870 Ha`. This is before/after identity at each fixed rank
count, not cross-rank bit identity. External `Si.pz-vbc.UPF` SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`;
cutoff `4 Ha`, five bands, Gamma, atomic start, default reference profile,
full diagnostics. The pseudopotential remains outside the repository.

Final local verification (all exit 0):

| Command | Result |
| --- | --- |
| `cargo check --workspace --tests` | Passed |
| `cargo test -p dft-engine --test refscf` | 65 passed, 11 ignored |
| `cargo test -p dft-engine` | 346 passed, 92 ignored across harnesses |
| `cargo fmt --check && cargo clippy -p dft-engine --all-targets -- -D warnings` | Passed |
| `DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf nlcc_ -- --include-ignored --nocapture --test-threads=1` | 9 passed, including serial World harness and Si fingerprint |
| `CHECK_SCCACHE=0 BAANDI_TEST_TMP=/private/tmp/g1-nlcc-gates scripts/check.sh` | Build, test, fmt, clippy PASS; 952 passed, 202 ignored, 95 binaries |
| `python3 -I scripts/tests/test_readme_figures.py` | 9 passed; existing generated evidence remains fresh |

The full gate's optional process-load diagnostic could not invoke `ps` inside
the sandbox; all four actual gates completed successfully (exit 0). Real MPI
is not installed or exercised here. README/figure result generation, dashboard
and session-registry updates are outside G1's explicit owned-file list and
remain the orchestrator's reporting boundary; no such files were changed.

## PBE

G2 enables unpolarized PBE in `refstart` and reference map fields. This
supersedes the historical PBE rejection in Start above. The existing CLI
unpolarized/fixed-occupation guard remains; collinear PBE is G4. Original
PW92 and explicit density floors remain outside this change. G2B adds NLCC
without changing the unpolarized/fixed-occupation guard.

For NLCC, `StartState.pbe_core` caches the normalized forward transform of
the immutable owned real core once per run (including each ladder-grid rebuild).
This differentiates exactly the same real core as the legacy evaluator, without
clipping, renormalizing or changing its G=0. Initialization and every map fill
reusable owned buffers with `n_total(r) = n_valence(r) + n_core(r)` and
`n_total(G) = n_valence(G) + n_core(G)` before the distributed GGA call. No
per-map core transform or new combined-input allocation is needed. Species
without core charge contribute zero. LDA/NLCC and PBE without NLCC retain
their previous arithmetic paths.

`dft_physics::gga::distributed` uses the run-owned `StartFields` plan through
an owned Fourier-field interface. Both initialization and each uncached field
evaluation reuse their existing normalized rho(G), before its Hartree
multiplier. For each Cartesian component the helper multiplies by iG and
inverts into owned real slabs. As in the serial GGA evaluator, each even-axis
Nyquist Miller component is zeroed **before** conversion with the reciprocal
cell matrix. The derivative is real and skew-adjoint on skew as well as
orthogonal grids; the Hartree Nyquist convention is unchanged.

The existing xcx GgaXPbe/GgaCPbe library evaluates rho and sigma=|grad rho|^2
on the owned points. Flux is 2*(vsigma_x+vsigma_c)*grad rho and
Vxc=vrho_x+vrho_c-div(flux). Each flux component takes a forward transform,
iG multiplication into a reciprocal divergence sum, followed by one shared
inverse. This is algebraically equivalent to serial componentwise inverses,
with floating-point rounding differences. Thus, given rho(G), the gradient
costs three inverse transforms and divergence costs three forward transforms
plus one inverse (seven transforms total). The density forward and Hartree inverse are shared
with the ordinary field evaluation. No full-grid payload gather occurs in
this helper or its start/map callers. Only scalar error/energy reductions and
the existing FFT transposes communicate. Even empty slab owners participate
in every collective; local library failures are agreed before later FFTs.

The density policy is the existing serial PBE policy: negative/nonfinite
total XC input fails collectively; n <= 1e-10 Bohr^-3 has zero local library inputs.
Gradients use the unscreened density and the divergence can be nonzero at a
screened point. No LDA absolute-density rule is imposed on PBE. Start density
clipping/normalization remains in the pre-existing atomic start. Prolonged
or mixed PBE density is not silently repaired. Field-sized buffers and xcx
outputs are allocated at the field-evaluation boundary, not per orbital/H
application; all payloads scale with owned slabs/sticks. Existing StartMemory
counters do not include these transient GGA arrays and are not PBE RSS bounds.

Energy is dv*sum(n_total*exc), summed locally then reduced together with Hartree
and the valence potential integral. Only XC sees the core: Hartree, occupations,
mixing, delivered density and charge verification remain valence-only. Double
counting uses dv*sum(n_valence*Vxc) with the **complete**
Vxc including divergence. The existing band form, Harris–Foulkes expression
and exact discrete XC functional remainder therefore apply unchanged; no
local-density Hessian approximation is introduced. LEAN_ENERGY cache,
continuing-map deferral, candidate/delivery audit and fresh verification
semantics are unchanged. The NLCC ledger records the core charge and the
separate core/valence potential integrals; the core integral is diagnostic,
never an extra energy term or part of the band double-counting subtraction.

G2B synthetic validation: `cargo test -p dft-engine --test refscf
pbe_nlcc_total_gradient_and_start_match_serial -- --nocapture` compares the
total-density gradient and starting potential with serial FFT/XC on three
odd/even skew grids at 1/3/8 simulated ranks (including empty slab owners).
`cargo test --release -p dft-engine --test refscf pbe_nlcc_delivered --
--ignored --nocapture` covers cubic/local and skew/nonlocal PBE NC fixtures
with core charge, thresholds 1e-8/1e-12 Ry, 1/3/8 simulated ranks, full/lean
energy audits and fresh verification. Gates are 5e-9/1e-10 Ha against the
certified legacy result, plus serial delivered XC energy/double counting.
These are synthetic Mac oracles, not real-MPI or QE/material qualification.

Measured G2B release results (same commands, full/lean, all three rank counts):
maximum legacy energy differences are 4.0734438045e-10 Ha at 1e-8 Ry and
2.5313084961e-14 Ha at 1e-12 Ry. Maximum total-gradient difference is
1.3877787808e-16 Bohr^-4 and starting-potential difference is
1.1241008124e-15 Ha. Full/lean delivered density and energy are bit-identical.
The intentionally under-resolved 1 Ha skew fixture has minimum total density
-7.1689028855e-4 electrons/Bohr^3 and is rejected by both serial PBE and
reference initialization at 1/3/8 ranks; no core clipping or tolerance change
is used to make this fixture run. Successful comparisons use explicit atomic
starts on both paths and the resolved fixtures encoded in the tests.

Validation and pending MPI8 measurements are recorded in
[the G2 request](../performance/lanes/RSCF/g2-pbe-dev08-request.md). Synthetic
Si diamond geometry is explicitly not a physical Si PBE pseudopotential
validation; no Si.pz-vbc functional override is used here.

## K points

G3-KPOINTS adds `refscf::k_driver::{run,run_with}` for validated weighted k
lists. The CLI supplies its existing mesh/list resolution and time-reversal
reduction, including summed partner weights. A single Gamma point with real
Gamma enabled dispatches directly to the unchanged real driver. A mesh that
contains Gamma uses the complex solver for all its points. No crystal symmetry
or density projection is introduced. Spin and Fermi–Dirac requests remain G4;
crystal symmetry reduction remains G6. The optional real-only RMM and cutoff
ladder experiments are explicitly rejected on the complex path.

### Complex Davidson and numerical contracts

`ref_davidson_complex` shares Options, occupied/empty thresholds, the smooth
preconditioner, guarded sqrt variant, restart guards, and accounting types
with `ref_davidson`, including its shared workspace capacity contract. It uses
conjugate products, fixed-tree real/imaginary
reductions, one common ethr supplied by the density loop, and the same strict
Ritz-change comparison. Empty bands retain max(5*ethr_Ha, 5e-6 Ha); this is
not an eigen-residual certificate. Reduced EVD may be replicated or root-only.
Every starting H column is applied unless an exact initial HX is supplied.

The overlap factor uses the existing real `PIVOT_RATIO_MIN`:
Cholesky pivot squared / original overlap diagonal >= eps/1e-12. Library
Cholesky, triangular solves and Hermitian EVD live in the approved
`dft-numerics::linalg::ComplexOverlapFactor` adapter. There is no handwritten
factorization, triangular inverse, new dependency, or new unsafe code. The
complex factor/EVD uses sequential faer even in an MKL build; complex GEMM
retains its existing backend dispatch. Native MKL execution has not been tested.

A rejected border preserves the old factor. Repair whitens the old basis and
HX together, preserves the physical Ritz block, projects corrections twice,
and removes directions at the existing 1e-12 norm threshold. Contraction
retains the lowest requested Ritz states. Measured restart refresh is the
baseline; the shared ideal-reset switch uses the same norm guard, four-reset
refresh interval and permanent fallback after a guard failure. Supplied-HX,
sqrt, restart-correction and buffer-swap variants are tested.

Row blocks, Gram/rotation matrices, factor/EVD scratch and packed collective
buffers are reserved before iteration. As with the existing adapters, opaque
backend packing and collective internals are outside that storage guarantee.
The real-only setup-fusion, single/wide GEMM, triangular-border, fused-append,
stacked-rotation, lean-clear and cross-call row-cache switches retain their
settings but report an explicit complex-adapter fallback in the run manifest.
They do not select a different mathematical rule or silently claim the real
optimized route. Complex per-call work and clocks describe the actual complex
solver; unsupported QE clock breakdowns are left unavailable.

### Per-k state, map, and distribution

Each k point has its own Bloch basis, owned complex orbital block, operator,
identity/epoch and kept HX. Kinetic energy and projectors use k+G; FFT support
uses integer G. The new distributed operator applies local FFT, kinetic and
nonlocal terms without gathering orbital blocks. All ranks participate in
one k point at a time. K-parallel pools are a later optimization.

All k points share one density grid. Atomic trial columns use the Fourier
transforms at k+G and atom phases, with owned-row orthogonalization and
lowest-kinetic coordinate complements. A projected Hermitian solve selects
the requested initial states. Setup reuses the common atomic-field initializer while skipping Gamma
orbital/projector construction and Gamma band-rank constraints. A small unused
Gamma FFT adapter remains in the common field-state type. Complex setup/output
memory estimates include
all retained k blocks; named-array estimates are not RSS measurements.

For each map, the driver installs the same potential in each operator layout,
solves each k with the same ethr, and forms the density in ascending k order:
occupation * normalized k weight * |u_nk|^2 / Omega. Time-reversal multiplicity
is already in the weight and is applied exactly once. Band, kinetic and
nonlocal energies use the same weights. The density functional, G=0 handling,
Coulomb dr2, Broyden update, ethr ratchet/first-step redo and independent
verification follow the real density loop. Verification records all k
residuals and weighted error estimates, not a first-k proxy.

Kept state compares exact X bits, operator identity and potential epoch per k.
A hit reconstructs delta-V times X on owned planes and adds it to HX; complex
slabs are not retained. Warm verification may borrow candidate state without
publishing its own state. Redo invalidates all entries. Inner exhaustion is
recorded per k and uses the same bounded six-call retry policy. Final results
contain all k blocks; the first-k fields remain compatibility views. CLI
ordinary and verification solve records include a zero-based `k_index` only
for the complex path. All-rank G distribution is reported as distributed,
not legacy replicated multi-k execution; outer k concurrency is one.

### Local validation (Mac, pure Rust, 2026-11-25)

No QE source was read and no QE run was performed. Legacy SCF is the independent
oracle here; its earlier QE comparisons are not new G3 QE evidence.

- Complex factor contracts: sizes 1/2/5/9 and reuse down to size 1, scales
  1e-3/1/1e3, independent generalized residual and metric, invalid metrics,
  and acceptance/rejection around the real pivot rule. Maximum relative
  residual 1.2374729395027181e-14; metric defect 1.7763568394002505e-15,
  dimensionless (both limits 1e-12).
- Complex solver: Hermitian dense oracle, dependent starts, root/replicated
  EVD, supplied-HX identity, switch variants, simulated ranks 1/3/8.
- Operator: 96 cubic/skew, k, projector, width and rank combinations. Maximum
  errors: H 6.206335383118183e-17 Ha; density 1.214306433183765e-17 e/Bohr^3;
  nonlocal expectation 5.204170427930421e-18 Ha. Additional tests cover empty
  owners, zero bands, delta-HX and collective rejection.
- End-to-end synthetic k list: kept-HX off/on, verified convergence, charge,
  cache hits and per-k records. CLI test exercises a shifted mesh with four
  time-reversal representatives and rejects G4/G6 combinations.

The external Si.pz-vbc.UPF SHA-256 is
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
It is not stored in the repository. Si2 and a skew deformation of the same
cell use 4 Ha cutoff, 5 bands, fixed occupations and identical legacy inputs.
Each row below covers simulated ranks 1/3/8 at 1e-8 Ry (standard) and 1e-12 Ry
(tight). Energy limits were 5e-9 and 1e-10 Ha respectively; all independent
verification attempts required for delivery passed.

| Cell / k list | Max standard energy error (Ha) | Max tight energy error (Ha) |
|---|---:|---:|
| Si2 / unshifted 2x2x2 | 1.3182477331952214e-9 | 1.0125233984581428e-13 |
| Si2 / shifted 2x2x2 | 3.1607605421868357e-11 | 2.7533531010703882e-14 |
| Si2 / explicit weighted list | 1.858868614590392e-10 | 6.750155989720952e-14 |
| Skew Si2 / unshifted 2x2x2 | 6.189271317680323e-11 | 2.4868995751603507e-14 |
| Skew Si2 / shifted 2x2x2 | 1.1889333961789816e-10 | 1.056932319443149e-13 |
| Skew Si2 / explicit weighted list | 4.0162007053368143e-10 | 7.105427357601002e-14 |

Per-k occupied eigenvalue differences are at most 5.329501798226333e-6 Ha
standard and 4.483715052364112e-8 Ha tight. Across all bands the maxima are
1.851278131737022e-5 and 1.3706304452687501e-5 Ha: the unchanged empty-band
Ritz threshold floor does not promise 1e-10 Ha empty eigenvalues. The largest
rank-count delivered-energy spread is 4.154099286779456e-11 Ha standard
(unshifted Si2), 3.730349362740526e-14 Ha tight. This is an accuracy comparison,
not a claim of identical complex trajectories across rank counts.

Reproduce the long/external-data check (ignored in the default gate):

```sh
DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine \
  --test refscf g3_k_mesh_legacy_matrix -- --ignored --nocapture
```

Gamma dispatch is checked bit-for-bit at ranks 1/3/8, including delivered
orbitals and scientific trajectories. The separate archived-base/main audit
uses `g3_gamma_base_fingerprint`, normalizing only the new empty Gamma k-index
metadata. The archived main commit
`58260cbe2387fe39a73a2839ab1efd2ebb728bc9` and this implementation produce
identical full fingerprint byte streams (not just matching hashes):

| Ranks | Bytes | SHA-256 |
| --- | ---: | --- |
| 1 | 100801 | `9afed278586f803a836e19d46470789de0c6e337c7ad59474bef8f2cdaef729d` |
| 3 | 100763 | `6f329ff3360630d73610dfbb0931de7c3f9d7626159b5315ec8aadc63a49ba5e` |
| 8 | 100729 | `0cc15ffcf4d25a1ddfff91c9c4740a0a900d8b4bc463f276c5d9da6f11cee57f` |

The final workspace `scripts/check.sh` gate exits 0: build, tests, format and
Clippy pass; 97 test binaries report 972 passed, 0 failed and 205 ignored.
The requested standalone workspace check, solver tests, reference SCF tests,
CLI reference tests and Clippy checks also exit 0. The sandbox denies the
gate's optional `ps` load probe; this does not skip any acceptance gate.
Evidence logs and fingerprint streams are under the ignored
`target/tmp/g3-validation/` directory. No changes are committed.

Native MPI, MKL and dev08 speed measurements remain unperformed. No performance
default changed. README generation, dashboard/session registry and the
performance ledger remain outside this ticket's approved ownership.

## Forces, stress and relaxation

G5 uses `refscf::delivered::force_state` (also exported as `stress_state`)
with the original `StartInput` and a converged `ScfOutcome`. This borrows the
stop map's density, orbitals, occupations, k weights and energy audit; the
independent verification map is never substituted. The independently recorded
force model identity still checks the cell, ordered atoms, pseudopotentials,
functional, cutoff, grid and occupation policy. Local/Hartree/XC derivative
fields are reconstructed from the delivered density and the original model.
The existing spin/finite-temperature views already carry per-spin/per-k
occupations and the E/F distinction, providing the extension boundary for G4.

A passed reference verification (ADR 0043 or the opt-in ADR 0044 assessment)
is accepted. This is a convergence precondition, not a rigorous derivative
error bound: the energy-only rule does not establish a density/force bound.
With an explicit `Verification::Off`, a converged delivery can instead grant
an internal unverified-delivery marker. Forces and stress then report
`verification = "unverified"`; neither a failed verification nor an exhausted
SCF grants that marker. The delivered verification report must be absent for
the opt-out; an empty report is not evidence of explicit opt-out. Changing a
failed result's requested verification mode
to Off does not make its failed attempts acceptable. Legacy derivative callers retain their existing certificate rule. G5B wires
`scf.verify = "off"` for reference `run` forces/stress only; the JSON derivative
blocks say `unverified` and retain false verification preconditions. Ionic
relaxation still requires verification on every step.

The current derivative kernels require replicated complex orbitals. The
adapter calls the existing collective `materialize_orbitals` boundary once;
subsequent force/stress calls reuse those arrays. Both Gamma and complex k
blocks are supported. Every rank participates, and no root-only full-orbital
bottleneck is introduced. This is **replicated derivative evaluation**, not a
claim of distributed force/stress kernels. For B bands, wave counts G_k and F
real grid points, delivered orbital storage is
`16 B (sum_k G_k + G_0)` bytes per rank (the common result retains a first-k
compatibility copy), plus `8 F` density bytes. Derivative potentials and FFT
workspaces are additional replicated grid arrays. These are logical payloads,
not RSS measurements; lazy-to-full conversion can temporarily coexist with
owned storage and collective buffers. No repeated orbital gather occurs in
the nonlocal term or in a stress call after a force call.

`BAANDI_FORCE_BAND_FIRST=on` (input `[performance] force_band_first = true`)
stacks all atomic projectors and their three Cartesian position derivatives,
then evaluates the overlaps for all bands in one adjoint GEMM per k block.
The remaining small contraction includes full off-diagonal projector coupling,
occupations and k weights. Extra scratch is `64 G_k P + 64 P B` bytes for P
projector columns at that k point. This changes floating-point reduction order;
the scalar implementation remains the default. No timing or speed adoption is
claimed. Synthetic identity coverage includes complex k orbitals, multiple
atoms, s/s/p channels, off-diagonal coupling and fractional/zero occupations.

Reference fixed-cell relaxation runs each ionic evaluation through the
reference driver. An in-memory owned `ladder::Seed` transfers both the preceding
step's orbitals and density through the existing injection/reorthogonalization
mechanism. Seed export checks the original model, orbital shape, finite density
and charge (relative tolerance 1e-8), before selecting the owned rows. Retained
seed payload is `8 B G_owned + 16 F_owned` bytes per rank, excluding metadata;
export additionally uses replicated density FFT work arrays. Ionic potentials,
start fields, mixer and verification are fresh.
Per-step record and restart formats are retained. A disk restart contains only
the existing density record, so the first resumed reference step reports a
fresh atomic start with an explicit missing-orbital-seed fallback; subsequent
Gamma and k-mesh steps are warm. G5B adds `IonicSeed::K`: the preceding
stop map's full complex blocks and real density are retained, then scattered
by the next start's actual owned G rows and field-plane offsets. Export checks
the source model and charge; injection rejects changed cell, cutoff, FFT grid,
k list/weights, band shape or charge. Davidson's generalized starting Gram
factorization handles the transferred block. Ionic fields, mixer and
verification are rebuilt; there is no exact-replay claim. K-seed storage is
`16 B sum_k G_k + 8 F` bytes per rank, excluding metadata. Cold initialization
is currently constructed before injection, so this is a correctness boundary,
not a claim to eliminate setup cost. Reference stress cutoff scans
use independently converged reference SCFs, sequentially so MPI collectives
cannot interleave between scans. G5B `stress_prepared` forwards the checked symmetry preparation to every
cutoff scan, including scan-grid compatibility and the irreducible list. The
CLI now reports the Pulay assessment on symmetric runs. Calling the plain
energy-scan entry on a symmetric reference state fails rather than silently
dropping its projector. Reference CLI stress requests retain the
existing stress default estimator threshold (5e-11 Ha); an explicit threshold
still wins.

R7a enables `ScfOptions::strain_reference` for real Gamma and complex k starts.
The reference cell selects the ordered wavefunction Miller lists and the
four-times-cutoff density support; `ScfOptions::basis` evaluates them with the
current cell's reciprocal metric, fractional k and volume, exactly as the
strict oracle does. Atomic positions, structure factors, radial transforms,
nonlocal projectors, kinetic energies and real-space quadrature use the current
cell. An explicit frozen `grid_override` is required at both entry points,
before the k initializer computes its shared grid. No cutoff-sphere rebuild or
resorting is performed at the strained cell. The None-reference branch and
SCF iteration arithmetic are unchanged.

`delivered::stress_state` retains the original options and current cell. Its
existing model identity includes the reference cell, so removing or replacing
that reference on delivery rejects. The analytic tensor differentiates the
frozen wavefunction and density supports at fixed fractional atoms/k and fixed
FFT dimensions. The separately reported Pulay assessment keeps the legacy
protocol: independent nominal-cutoff SCFs at the current cell with
`strain_reference = None`, a common enlarged grid and +5%/+10% cutoffs; the
estimate is not silently added to the frozen tensor. Complex ionic seeds also
retain the reference-cell identity and reject a changed reference at injection.
R7b routes cell and ions_cell evaluations through `api::solve`, with the
strained cell and frozen reference/grid. `IonicSeed::in_cell` carries the
unchanged ordered Miller coefficients (including complex k blocks) and scales
the fractional-grid density by Omega_old/Omega_new, conserving electron count.
Gamma seeds use the existing owned-row injection; complex seeds validate the
source segment and target mesh/grid before injection. Ionic fields, mixer and
verification are rebuilt. Rebuild/polish evaluations are fresh; disk resume
has the same explicit missing-orbital fallback as ionic relaxation.

D5: port C3, retaining `start_extrapolation` and `RELAX_EXTRAPOLATION` as
opt-in, default none. The C3-EXTRAPOLATION ledger measured six-step Si8
SCF counts 84/76/70 and Mac wall 5.98/5.54/5.27 s (none/first/second order),
with gains on a synthetic slab, mixed results on fcc and a regression on a
mismatched Gaussian atomic charge. These are local observations, not a dev08
speed adoption or a new default. The existing atomic-superposition builder
supplies n_at without running a legacy SCF; the predictor replaces only the
seed density with n_at(R_new) plus extrapolated bonding density, clips negative
samples and restores charge. Delivered orbitals are retained and the next
reference solve reorthogonalizes them. Existing history/restart and fallback
labels remain. As before, C3 applies only to ions mode; cell modes use affine
transfer. L4 and L5 no longer select the legacy driver.

The ignored `cell_relaxation::reference_cell_protocol` test reads the stored
P23A synthetic starts and pressures from `validation.json`, adds atomic start
data to the same Gaussian potential, asserts actual reference/legacy routing,
and compares trajectories at a tightened 5e-15 Ha SCF estimator threshold:
per-step energy <= 5e-9 Ha, final lattice <= 1e-5 Bohr, pressure <= 1e-6
Ha/Bohr^3. `reference_cell_seed_and_predictor` covers combined ion/cell,
complex k, periodic rebuilds and ionic C3 modes. Run both with
`cargo test --release -p dft-cli --test all cell_relaxation::reference_ -- --ignored --nocapture`.
Local validation: `reference_cell_seed_and_predictor` passes for the
unsymmetrized combined, complex-k and rebuild trajectories. Maximum per-step
energy differences are 7.03e-11, 1.94e-11 and 6.07e-12 Ha respectively; maximum
final lattice difference is 1.43e-9 Bohr and pressure difference 3.43e-12
Ha/Bohr^3. Both C3 orders and the none arm pass their energy/provenance checks.
The existing external Si2 k222 ionic test also passes (maximum dE 6.50e-12 Ha).

R7b continuation after R8a: `api::run_request` now forwards both the K seed
and prepared symmetry into the existing `k_driver::run_continuation`.
No stepper loop, per-map allocation or arithmetic order changes. Prepared
Gamma solves retain their complex-k representation. `Observables::ionic_seed`
exports using the delivered symmetry identity; it supports both deferred CLI
export and immediate `save_ionic_seed`. Plain export still rejects a symmetric
state, and prepared export rejects stale geometry. A direct test compares
seeded `solve()` with explicit stepper advancement: energy bits, density,
orbitals and forces are identical.

The protocol uses an explicit one-vector mixing history in both arms at the
unchanged 5e-15 Ha estimator threshold. The compressed synthetic case's default
history otherwise becomes singular at iteration 5 (condition 4.28e15, eigenvalues
2.89e-19 and 1.24e-3 Ry) before the threshold is met. No numerical default or
accuracy tolerance changes. Tetragonal-to-cubic group changes can invalidate
warm density/model identity in both drivers; the test permits only those
recorded first-attempt failures followed by the existing fresh half-step retry.
The complete stored synthetic protocol converges in both arms, with seven
frames per isotropic case and eleven for tetragonal relaxation. Same-build Mac
release comparisons (reference versus explicit legacy pin):

| Case | Maximum per-step dE (Ha) | Final lattice difference (Bohr) | Final pressure difference (Ha/Bohr^3) |
| --- | ---: | ---: | ---: |
| P=0, compressed start | 3.03e-12 | 2.71e-14 | 9.95e-17 |
| P=0, expanded start | 1.61e-11 | 1.36e-13 | 4.24e-16 |
| P=0.001, compressed start | 3.19e-11 | 7.99e-15 | 1.04e-17 |
| P=0.001, expanded start | 1.47e-10 | 2.89e-14 | 3.37e-16 |
| Tetragonal to cubic | 1.37e-9 | 1.84e-9 | 2.42e-13 |

`cargo test --release -p dft-cli --test all cell_relaxation -- --ignored --skip real_files --nocapture`
passes all eight tests (exit 0), including the existing restart, independent
energy/derivative, rebuild/FIRE and standard-versus-strict checks. The full
`CHECK_SCCACHE=0 BAANDI_TEST_TMP=/private/tmp/r7b-approved-gates scripts/check.sh`
gate passes build, test, fmt and clippy (exit 0; 1030 passed, 230 ignored).
All ticket TEST_COMMANDS also exit 0. No tolerance was changed. These are
correctness results, not performance or QE evidence.

Reference cutoff scans now preserve the caller's mixer coefficient and history
length (plain mixing retains zero history) as well as the estimator threshold.
Previously `stress::energy_at_prepared` discarded those mixer settings and
constructed a default eight-vector history. Standalone final-cell scans at
12.0/12.6/13.2 Ry, with the same 10x10x10 grid and 5e-15 Ha threshold, reproduced
singular histories (condition 2.42e15) with history eight; all three converged
and passed independent verification with history one. This fix is confined to
reference cutoff-scan setup; the SCF iteration and legacy scan paths are unchanged.
Real Si/MgO cell protocols and dev08 validation have not been rerun.

R7a tests append to `refscf_forces`: missing-grid rejection, explicit ordered
Miller/kinetic/density-support checks across cutoff-shell crossings, standard
(5e-9 Ha) and tight (5e-13 Ha estimator) comparisons against the strict oracle,
and xx/symmetric-xy energy finite differences at strain steps 2e-4 and 1e-4.
The tight acceptance limits are 1e-10 Ha in energy and 1e-8 Ha/Bohr^3 in stress;
finite differences and the separate Pulay comparison use 1e-8 Ha/Bohr^3.
The external Si2 matrix uses 4 Ha and 5 bands at Gamma and full k222 and is
ignored in debug. Fast synthetic tests cover skew Gamma and a time-reversed
two-point complex k list (the synthetic full k222 fixture was noninsulating
and is not an insulating oracle). Run
`DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf_forces frozen_strain -- --include-ignored --nocapture`.
`unstrained_delivery_fingerprints` prints exact same-build fingerprints of
energy, density, occupations, eigenvalues, orbitals and stress for comparison
before/after this change; these are not portable floating-point goldens.

R7a local evidence (2026-11-25, Mac pure Rust, strict engine oracle; no QE
run): the release frozen-strain command passes all four tests. Si2 uses
`A = (I + eta) A_ref`, with diagonal eta = (0.09, 0.08, 0.07),
eta_xy = eta_yx = 0.012, and displaced fractional second atom
(0.253, 0.249, 0.252). Both wavefunction and density cutoff-shell crossings
are asserted. External Si UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.

| Si2 sampling | Standard absolute dE (Ha) | Tight absolute dE (Ha) | Tight max tensor difference (Ha/Bohr^3) | Max FD error, xx/xy and both steps (Ha/Bohr^3) | Pulay estimate difference (Ha/Bohr^3) |
| --- | ---: | ---: | ---: | ---: | ---: |
| Gamma real | 7.885e-10 | 1.075e-13 | 4.657e-11 | 1.133e-10 | 7.601e-13 |
| Full k222 complex | 2.282e-10 | 8.260e-14 | 2.020e-10 | 3.970e-10 | 6.167e-12 |

The standard-threshold stress differences reach 2.343e-8 Ha/Bohr^3; the
1e-8 derivative gate applies to the tight runs, not the standard energy
threshold. Frozen-basis differentiation does not establish cutoff convergence.
Same-build unstrained synthetic Gamma and two-point complex-k fingerprints
are unchanged from the pre-edit production code at
`60498b5ff8b0e44eb208b07e13873d1fae4a6107`: respectively
`dbbc517c35af3115` and `56d403ae498a5177` (test-profile build).
README evidence ingestion, dashboard/session publication and dev08 validation
remain with their owners; no tracked files outside this ticket's owned set were edited.

R7a verification: `cargo check --workspace --tests`, the requested
`refscf_forces` (10 passed, 3 ignored) and `refscf` (84 passed, 24 ignored)
targets, `cargo fmt --check`, and engine all-target clippy with `-D warnings`
pass. Full `scripts/check.sh` passes build/test/fmt/clippy: 106 binaries,
1018 passed, 0 failed, 220 ignored. The sandbox run sets
`BAANDI_TEST_TMP=/private/tmp/r7a-check CHECK_NICE=0 CHECK_SCCACHE=0`; its
optional process-load probe cannot invoke `ps`, but no gate is skipped.
README freshness tests pass (9 tests).

Validation commands: `cargo test -p dft-engine --test refscf_forces` includes
synthetic finite differences, ranks 1/3/8, verification rejection/opt-out,
Gamma ionic seeds and lazy materialization (idempotent, bit-identical
force/stress against eager delivery). External-file Si2 tests are ignored by
default and run with
`DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine --test
refscf_forces -- --include-ignored --nocapture`. The real Si2 three-step ionic
trajectory oracle is `DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release
-p dft-cli reference_relax_tests -- --ignored --nocapture`. UPF files stay
external. These compare against the existing engine and reference-energy
finite differences, not a new QE run. Force and stress comparison tolerances
are respectively 1e-6 Ha/Bohr and 1e-6 Ha/Bohr^3; the trajectory also checks
energy differences below 5e-9 Ha and positions below 1e-6 Bohr.

G5B closes the previously recorded complex-k seed, CLI opt-out and symmetry
cutoff-scan gaps. The stop-to-verified reader recognizes `reference-spin-fd-1`
using its top-level verification and threshold. Those historical outputs lack
monotonic map boundaries: the reader reports the recorded verdict and attempt
durations, with a null interval and `parts_sum_checked=false`; it never invents
clocks. Gamma and fixed-occupation k outputs retain full interval validation.
The executable hash sidecar uses schema v2 with its recording time: mtime or
ctime within two seconds of recording or lookup is unverified and rehashed.
An entry recorded during the timestamp collision window stays suspect even
when later read after that window. Deterministic tests force 50 stale-digest
metadata collisions and also test the safe hit path with a controlled clock.

G5B additional commands (external files stay outside the checkout):
- `cargo test -p dft-cli --test reference_gaps`
- `DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-cli --test reference_gaps reference_k_relax -- --ignored --nocapture`
- `DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf_forces symmetric_si2_pulay -- --ignored --nocapture`

Local validation (Mac, pure Rust; continuation review, 2026-11-25): all five
G5 `TEST_COMMANDS` pass. `scripts/check.sh` also passes build/test/fmt/clippy
(979 passed, 209 ignored across 98 binaries). Its optional process-load
inspection could not run in the sandbox; this did not skip any gate. README
freshness checks pass (9 tests); adding this evidence to the generated README
and dashboard remains with the owners of those files.

The external release derivative command above passes all 6 tests. The
comparison matrix covers cubic/skew synthetic nonlocal cells at Gamma and a
weighted two-point k list, displaced Si2 at Gamma and that k list, and skew
Si2 at Gamma; each is compared at simulated ranks 1/3/8. Finite differences
hold the FFT grid fixed and assert unchanged wavefunction and density Miller
sets for strain. They check Cartesian atom-1 x force and xx stress; they do
not establish cutoff convergence or every tensor component independently.

| Quantity | Maximum observed absolute difference | Acceptance |
| --- | ---: | ---: |
| Reference vs legacy energy | 3.781e-13 Ha | 5e-9 Ha |
| Reference vs legacy force component | 8.926e-8 Ha/Bohr | 1e-6 Ha/Bohr |
| Reference vs legacy stress component | 6.373e-10 Ha/Bohr^3 | 1e-6 Ha/Bohr^3 |
| Reference-energy finite-difference force | 3.558e-7 Ha/Bohr | 1e-6 Ha/Bohr |
| Reference-energy finite-difference stress | 5.879e-10 Ha/Bohr^3 | 1e-6 Ha/Bohr^3 |
| Band-first vs scalar nonlocal fixture | 8.674e-17 Ha/Bohr | 1e-12 Ha/Bohr |

Derivative oracles deliberately tighten the standard estimator to 5e-13 Ha
(5e-12 Ha in the certificate/opt-out fixture); these are not a guarantee of
1e-6 derivative accuracy at the default energy estimator. No production
threshold is tuned. Si2 uses 4 Ha wavefunction cutoff and 5 bands; synthetic
fixtures use 2 Ha and 3 bands. The external Si UPF SHA-256 is
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
There is no new QE reference, real-MPI validation, timing or switch adoption.
The additional CLI test
`cargo test -p dft-cli --lib reference_k_mesh_records_and_unsupported_combinations`
passes verified complex-k force/stress JSON and the three reference cutoff
scans at the default stress threshold.

The external three-step Si2 release relaxation test also passes (4 frames,
including the initial frame). Every subsequent reference frame records
`warm_transfer`. Maximum differences against the legacy trajectory are
7.177e-13 Ha in energy, 7.697e-8 Ha/Bohr in force and 2.027e-7 Bohr in Cartesian
position, under unchanged gates of 5e-9 Ha, 1e-6 Ha/Bohr and 1e-6 Bohr.
This trajectory is serial Gamma; rank coverage applies to the derivative and
ionic-seed tests, not a real MPI CLI relaxation campaign.

G5B local evidence (Mac, pure Rust; 2026-11-25): the eight-point 2x2x2 Si2
trajectory completes its prescribed three ionic steps (four frames), with
`warm_transfer` on every subsequent reference frame. Both optimizers stop at
the step budget, not at the deliberately tight ionic convergence tolerance.
At the 5e-13 Ha electronic threshold, maximum reference/legacy differences
are 6.497e-12 Ha in energy, 6.789e-8 Ha/Bohr in forces and 3.599e-7 Bohr in
positions (gates 5e-9 Ha, 1e-6 Ha/Bohr and 1e-6 Bohr). The unchanged Gamma
trajectory oracle also passes. These are engine comparisons, not QE evidence.

For symmetric undisplaced Si2, the eight-point mesh reduces to three points.
The frozen tensor differs from the full-mesh reference by 1.559e-10
Ha/Bohr^3 and the Pulay correction by 5.947e-10 Ha/Bohr^3; the corrected
diagonal components also pass the 1e-8 Ha/Bohr^3 gate. The correction remains
an isotropic cutoff estimate, not a cutoff-convergence certificate.

The synthetic k-seed test passes at simulated ranks 1/3/8: maximum warm/cold
energy difference 4.863e-13 Ha and force difference 9.420e-8 Ha/Bohr.
The CLI fixture delivers identical energy, density cube, all k eigenvalues
and occupations, forces and stress with verification on/off: 22 ordinary
maps in each run, four per-k verification solves versus zero. Actual k, FD
and fixed-sector spin results pass the offline reader, including the CLI
subcommand for FD/spin. The six cache tests include 50 forced same-metadata
rewrites and concurrent publication. No real MPI, dev08 or speed adoption is
claimed; generated README and dashboard updates remain outside G5B ownership.
The full `scripts/check.sh` gate exits 0: build, test, fmt and clippy pass;
992 tests pass, 217 are ignored across 99 binaries. All six ticket command
groups pass, as do the explicit release trajectory and symmetry tests above.
README freshness checks pass (9 tests); no generated evidence files were edited.
The sandbox denies optional `ps`/`nice` load inspection; verification still ran.


## Spin and Fermi-Dirac

G4 adds the independent `refscf::electronic::{run,run_with}` entry. Its policy
is scalar Fermi–Dirac, fixed collinear sector, or collinear Fermi–Dirac with a
common chemical potential or fixed moment. Gamma uses complex orbitals for
these policies, as required by ONEPATH-3. Existing fixed unpolarized Gamma
and k entry points retain their arithmetic and even-electron restriction.

Stage 1 solves every retained band at every k before calling the existing
`occupations::FermiDirac` number/entropy solver. Every retained FD band is
mandatory in reference Davidson, including numerically zero tails. Density,
kinetic, nonlocal, band energy and the residual-square guard all use the
same occupation times normalized k weight exactly once. The explicit FD cap
is prepared before construction (CLI default 32, identical to legacy).
Positive theta, interior sector capacities, number/moment residuals and the
legacy top-band coverage limits are enforced without clipping occupations.
Internal energy E, positive TS and Mermin F=E-TS are separate records.

Stage 2 uses separate channel/k operators, orbitals and kept HX entries.
Hartree sees total valence charge; XC sees both channels and half the immutable
NLCC core in each. Only LDA-PW has an available legacy polarized LDA evaluator.
Polarized LDA-PZ is absent in legacy and rejects explicitly. G4C adds
distributed two-channel PBE, with or without NLCC, for fixed sectors and
both FD spin policies. Scalar FD PBE retains the existing distributed GGA.
No new XC parameterization is introduced.

Broyden receives one stacked (n,m) vector and one shared history. Both nonzero
blocks use the existing half-sphere Coulomb metric; the magnetization zero
mode adds 4*pi*Omega/G_min^2 in Ry, exactly twice ONEPATH-3's E_0. Fixed-sector
and fixed-moment zero modes are imposed from their constraints; free moment
is never projected away. The estimator is the sum of charge, spin and uniform
spin terms. The standard occupation-weighted residual-square guard and number
conditions also apply. Fixed sectors require each nonempty channel's sampled
indirect gap, including on verification. Empty channels have no gap gate.

Verification repeats the entire map at the candidate output, including fresh
chemical potentials and occupations. It uses the next reference ratchet and
ONEPATH-3 Coulomb/guard judgment. A failed check preserves the original
candidate and Broyden history; only an ordinary update advances the run.
Delivery contains the candidate's own densities, eigenpairs, occupations,
chemical potentials, entropy and total/absolute magnetization. Verification
results never replace it. Energy-only verification, failed-map reuse, RMM and
cutoff ladder currently reject; predictor gates and deferred energy audits
are explicitly inactive on this entry.

### G4 continuation review and acceptance status (2026-11-25)

**Scalar FD and fixed collinear sectors are qualified; collinear FD moment policies are routed to legacy.**

| Policy | Reference-driver status |
| --- | --- |
| Scalar FD (including Si2) | Qualified and enabled |
| Fixed collinear sector | Qualified and enabled |
| Collinear FD, fixed moment | Routed to legacy; G4B pending |
| Collinear FD, free moment | Routed to legacy; G4B pending |

`BAANDI_REFSCF_SPIN_FD_MOMENT=on` (or
`[performance] refscf_spin_fd_moment='on'`) enables the two unqualified
policies for G4B development only; its default is off. Default fallback is
recorded as `reference driver: collinear Fermi-Dirac not yet qualified, G4B`
in the result and switch manifest. No acceptance bound is relaxed.
The passing ignored `g4_synthetic_matrix` covers the 24 qualified comparisons.
The separate ignored `g4b_moment_acceptance` retains all 24 moment-policy
comparisons and is expected to fail until G4B qualifies them.
The approved stopping estimator and the Broyden singular-history rejection
remain unchanged. In particular, the tests retain the 5e-9 Ha standard,
1e-10 Ha tight, 1e-4 mu_B moment, 1e-4 Ha chemical-potential and 1e-4 electron
weighted-occupation bounds. Passing the estimator alone is not evidence that
these delivered-error bounds hold.

Review corrections to the interrupted implementation:

- Require every mandatory FD/gap-edge Ritz flag to converge before accepting
  an inner call; the fixed driver's allowance for pending buffer roots cannot
  apply to mandatory roots. Six calls still bound retries.
- Prepare the polarized XC adapter, its spectral scratch, channel arrays and
  immutable core once per owned slab; reuse them across ordinary and fresh
  verification maps. Empty slab owners join the reductions.
- Write the CLI run-completion/binary-identity manifest, which was missing
  from the early-return branch. Tests check spin/k/band array sizes, weighted
  charge, common/per-spin chemical potentials, fixed moments and E-TS=F.
- Test verification delivery independence (density, orbitals, eigenvalues,
  occupations, chemical potential and free energy are unchanged by a passing
  verification), and reject noncollinear/Gaussian/MP/cold requests.

The synthetic fixture is a one-electron, nonlocal Gaussian atom in the skew
cell `[[5,0,0],[0.2,5,0],[0,0.1,5]]` Bohr at `[0.13,0.21,0.32]`, 2.5 Ha
cutoff. FD uses theta=0.01 Ha and 12 bands; fixed spin uses 4 bands and sector
(1,0). The free-moment seed is 0.3 mu_B; the constrained FD moment is 0.3 mu_B.
Both Gamma and an unshifted 2x2x2 mesh run at simulated ranks 1/3/8 and
thresholds 1e-8/1e-12 Ry. The reference profile is the production default
(beta=0.7, history=8), with verification on. The strict legacy oracle uses
eigen/density tolerances 1e-10 and energy tolerance 1e-11 Ha. Its free-moment
case uses an explicitly configured full linear step: legacy Pulay oscillated,
and beta=0.2 linear mixing failed positivity because its mean correction
imposes the new channel count on a damped density. No legacy code changed.

Each numeric entry below is the maximum over both k lists and all three ranks.
Fixed-spin F is E (TS=0).

| Policy | Standard | Tight | Max standard abs(delta M) (mu_B) | Max standard abs(delta mu) (Ha) |
| --- | --- | --- | ---: | ---: |
| Scalar FD | PASS, abs(delta F)=4.5800058190e-11 Ha | PASS, 6.6752159356e-15 Ha | 0 | 8.5941571258e-7 |
| Fixed spin | PASS, abs(delta E)=1.9797810613e-10 Ha | PASS, 1.6253665081e-13 Ha | 1.9984014443e-15 | not applicable |
| Fixed-moment FD | Energy PASS, abs(delta F)=6.4546469558e-10 Ha; occupation FAIL | PASS, 1.6063539388e-14 Ha | 5.3623772089e-14 | 3.4210860895e-6 |
| Free-moment FD | FAIL, abs(delta F)=1.1167172215e-8 Ha | FAIL, singular histories | 3.3922507799e-4 | 2.1352597911e-4 |

Of the 48 comparisons, 33 pass all gates and 15 fail. Scalar FD and fixed
spin pass all 24 comparisons. Fixed-moment FD passes its energy, moment and
chemical-potential bounds, but all three standard mesh runs miss the weighted
occupation bound (maximum 1.4082986520e-4 electrons against 1e-4); its other
nine comparisons pass. All 12 free-moment comparisons fail acceptance.
Gamma standard passes its single verification
but misses the energy, moment and chemical-potential bounds. Mesh standard
meets the energy class (maximum 2.4829952330e-9 Ha) but misses the moment
and weighted-occupation bounds (maximum 3.3922507799e-4 electrons).
Tight free-moment runs reject at mixing iteration 10 (Gamma) or 13
(mesh), with history condition numbers about 1.54e16 and 8.86e14. Replacing
these errors by silent regularization would contradict the existing Broyden
contract. The interrupted beta=0.5/history=4 test profile also missed the
standard mesh moment bound (1.9081328586e-4 mu_B); it is not a default fix.
The failing default-profile regression is retained and collects errors across
the full matrix, rather than terminating at the first inaccurate result.

Si2 with the external Si.pz-vbc.UPF (SHA-256
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`), 4 Ha
cutoff, 8 bands and theta=1e-4 Ha passes all 12 Gamma/mesh/rank/threshold
comparisons. Maximum free-energy differences from strict legacy FD are
3.2998670463e-10 Ha standard and 4.9737991503e-14 Ha tight; differences from
strict integer occupations are 3.2995650656e-10 and 4.7073456244e-14 Ha.
The UPF remains external.

Closed-shell Gamma is byte-identical to the archived ticket-base commit
`d235143fcf96d2e18771af4b965c77d60343cfa5`, whose crate sources were checked
against Git blob hashes. Full scientific-history/orbital fingerprints at
ranks 1/3/8 have respectively 100801/100763/100729 bytes and the same SHA-256
values recorded in the G3 table above. This is a base comparison, not only
an in-process dispatch comparison.

Reproduction (`g4b_moment_acceptance` intentionally fails):

```sh
cargo test --release -p dft-engine --test refscf g4_synthetic_matrix -- --ignored --nocapture
cargo test --release -p dft-engine --test refscf g4b_moment_acceptance -- --ignored --nocapture
DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf g4_si_gap -- --ignored --nocapture
cargo test --release -p dft-engine --test refscf g4_free_moment_tight_short_history -- --ignored --nocapture
G3_GAMMA_DUMP_DIR=/existing/output/directory cargo test --release -p dft-engine --test refscf g3_gamma_base_fingerprint -- --nocapture
```

Raw logs and fingerprint streams are in ignored `target/tmp/g4-validation/`.
All five ticket TEST_COMMANDS exit 0: workspace check, solver tests (194
passed, 25 ignored), reference tests (78 passed, 19 ignored), CLI reference
tests (15 passed, 9 ignored), and format/Clippy with warnings denied. The full
`CHECK_SCCACHE=0 BAANDI_TEST_TMP=/private/tmp/g4-spin-fd CHECK_TEST_JOBS=2
scripts/check.sh` exits 0: build/test/fmt/Clippy pass, with 977 passing tests
across 97 binaries (209 ignored in that invocation). The later standalone
reference test invocation also compiles the added ignored short-history
diagnostic. The sandbox denies the optional `ps` load probe; no gate is
skipped. External scratch is used for cache tests that reject checkout paths.
Before the architect routing decision, the release
`g4_ --include-ignored --test-threads=1` run, with the Si
UPF environment set, exits 101: five test functions pass and the synthetic
acceptance matrix fails. Regular gates exclude that documented long test,
so their success does not establish G4 acceptance.
The short-history command is a diagnostic, not a substitute for acceptance.
It passes all six free-moment Gamma/mesh/rank comparisons at 1e-12 Ry with
beta=0.5/history=4: maximum abs(delta F)=1.2207179712e-12 Ha,
abs(delta M)=1.5505818965e-6 mu_B, abs(delta mu)=2.2515379179e-6 Ha and
weighted occupation difference=1.5505818992e-6 electrons. This separates the
map's ability to reach the legacy solution from qualification of the default
stop and history policies; it does not remove their recorded failures.
Resolution requires a reviewed collinear-FD stopping/verification policy and
a reviewed singular-history recovery policy, followed by the unchanged gates.
No QE source was read, no QE run was performed, and simulated ranks establish
no real-MPI, MKL or speed claim. Dashboard/registry, README evidence generation
and performance ledgers are outside this ticket's explicit owned files.


After the architect routing decision, the release `g4_ --include-ignored
--nocapture` run with the external Si UPF exits 0: all six test functions pass.
This includes the 24 qualified synthetic comparisons, 12 Si2 comparisons,
and the short-history diagnostic; it excludes `g4b_moment_acceptance`.
CLI regression coverage checks default legacy routing and the exact recorded
reason for both free and fixed moments, then enables the registry switch and
checks reference records. Fixed-sector and scalar FD routing remain enabled.

The post-routing full `scripts/check.sh` run exits 0 (build, test, fmt,
Clippy): 977 passed, 0 failed and 211 ignored across 97 binaries. The
standalone reference suite has 78 passed and 20 ignored, including the new
G4B acceptance test. Sandbox `ps`/priority diagnostics do not skip any gate.
All five ticket TEST_COMMANDS were rerun after routing changes and exit 0.
Post-decision logs and exit codes are saved under
`target/tmp/g4-validation/architect-*`.

### G4B diagnosis before candidate changes (2026-11-25)

The unchanged arithmetic was instrumented with an optional, preallocated,
in-memory `electronic::Trace`. It survives a rejected mixing update. Adapters
write it only after return. `G4B_TRACE_DIR` on the G4 acceptance tests writes
one CSV per policy/mesh/rank/threshold, including every ordinary, first-map
redo and fresh verification map. Each row records charge/nonuniform-spin/
uniform-spin terms (Ry), absolute errors against the independently converged
strict legacy oracle (F in Ha, M in mu_B, mu in Ha, k-weighted occupation L1
in electrons), number/moment residuals, history condition, individual Gram
diagonals and maximum absolute pair cosine. Zero condition means no history
solve at that map, not a perfectly conditioned history.

Baseline command (exit 101, expected scientific acceptance failure):
`G4B_TRACE_DIR="$PWD/target/tmp/g4b" cargo test --release -p dft-engine --test refscf g4b_moment_acceptance -- --ignored --nocapture`.
The first invocation used a relative trace path and failed to write from the
crate test working directory; the absolute-path rerun produced all 24 traces.

Representative one-rank baseline rows (all bounds unchanged):

| Policy/map | charge Ry | spin Ry | uniform Ry | abs(dF) Ha | abs(dM) mu_B | abs(dmu) Ha | occupation L1 electrons |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| Free Gamma standard, delivered 7 | 7.20513483e-12 | 1.24654425e-9 | 5.51814786e-9 | 1.11671672e-8 | 1.48839490e-4 | 2.13525967e-4 | 1.48839490e-4 |
| Free Gamma verification 7 | 8.00092706e-13 | 1.38348208e-10 | 6.16972293e-10 | 1.28744398e-9 | 5.03946931e-5 | 7.29032180e-5 | 5.03946931e-5 |
| Free mesh standard, delivered 11 | 3.20209006e-10 | 9.70878042e-10 | 5.15123959e-9 | 2.48298660e-9 | 3.39224954e-4 | 1.02320082e-5 | 3.39224954e-4 |
| Fixed moment mesh standard, delivered 4 | 2.01173880e-9 | 4.58689136e-9 | 0 | 6.45450964e-10 | 4.82947016e-15 | 1.19614857e-6 | 1.40829865e-4 |

The standard estimator and its fresh-map check accept before the delivered
state meets the oracle bounds. A verification map's better energy does not
repair the candidate it certifies. Free-moment error is predominantly in the
uniform term: this is not evidence that short-wavelength under-weighting is
the only problem. Local XC stiffness is therefore a candidate, not a promised
bound or a substitute for kinetic/occupation response.

Tight Gamma rejects at iteration 10 with condition 1.54081746e16, pair squared
norms spanning 3.71200062e-10 to 2.61055520e-3 Ry, and maximum pair cosine
0.9999957651. Its true F error decreased from 6.60563657e-9 (map 8) through
1.20780239e-9 (map 9) to 8.94843921e-11 Ha (map 10). Tight mesh rejects at
13 with condition 8.85736543e14 and cosine 0.9999959382; true F error decreases
2.48298660e-9 -> 6.39018824e-10 -> 4.06599754e-12 Ha at maps 11–13.
These are near-dependent, widely scaled secants on a converging tail, not
an observed persistent stall. Both failures are free-moment runs: the uniform
spin mode is retained, so fixed-M/zero-mode projection cannot explain them.
The same metric/projection succeeds for the constrained tight cases.

Both paths call `FermiDirac::new(theta)` with electron tolerance 1e-13 times
the requested sector count and the same normalized k weights. Both occupation
residual definitions are sum_ski w_k abs(f_new-f_old), with no second capacity
factor. Fixed-M mesh map 4 has number residual -1.50990331e-14 electrons and
moment residual -2.22044605e-15 electrons. Its occupation-oracle discrepancy
is redistribution between states, not failed channel-number solves. No
occupation-root tolerance/definition fix is justified.

### G4B candidate contracts

Two independent switches are appended to the registry, both default off:
`REFSCF_SPIN_STIFFNESS` and `REFSCF_HISTORY_CONDITIONING`. They apply only to
`Policy::SpinFermiDirac`; scalar FD, fixed sectors and closed-shell reference
drivers ignore both. The existing development-routing switch is separate.

The stiffness candidate **adds**, in Ha,
`E_xc,m = (1/2) integral |f_mm[n,m]| (m_out-m_in)^2 dr`, where
`f_mm = d[(v_up-v_down)/2]/dm` at fixed total n. The signed XC curvature alone
is frequently negative and cannot be a positive stop metric; its absolute
value prevents cancellation of positive residual errors. This is not the
full thermodynamic spin stiffness (kinetic and occupation response are absent),
nor an energy-error bound. Unlike Coulomb weights it does not attenuate short
wavelengths by 1/G^2, as motivated in ONEPATH-3 section 2. The old uniform and
nonuniform spin metric remain present; the new term also sees uniform errors.

The same existing polarized LDA-PW evaluator supplies potentials. Because its
public adapter exposes first derivatives only, the candidate uses two point-local
potential evaluations at fixed n with relative displacement cbrt(f64 epsilon).
This balances central-difference truncation and cancellation without a fitted
material scale. At polarization endpoints it uses the feasible one-sided span;
zero total density contributes zero, with the evaluator's existing vacuum
screening unchanged. Spin-neutral NLCC is included before differentiation.
Scratch channels and finite-difference buffers are constructed once. No new
XC library, fitted coefficient, or parameterization is introduced. This has
not been timed on dev08. A separate energy second-difference test checks the
curvature including NLCC and the exact zero-residual limit.

Ordinary maps, first-map redo, eigenvalue-tolerance ratchet, and independent
verification all use the augmented energy. The coefficient 2 converts this
Ha term to the driver's Ry estimator. Verification still delivers the original
candidate, never its fresh-map output. The history metric itself is unchanged.
The existing `dr2_ry` and verification `coulomb_error_ha` fields carry the
augmented estimator when stiffness is selected; opt-in traces separate all
three Coulomb terms from `stiffness_ha`.
History conditioning follows proposed ADR 0045: repeatedly remove the oldest
pair until rcond > sqrt(p*epsilon), with every removal's pre-drop condition
recorded. Nonfinite/backend errors still reject. The off branch retains the
original ring/cache, numerical singularity rejection and arithmetic. Candidate
drop metadata is kept separate from legacy `Diagnostics` to preserve existing
closed-shell scientific-record bytes.

<details>
<summary>G4B baseline per-map oracle trace (one rank)</summary>

All energies below are Ha except the three explicitly labelled Ry estimator terms.
`v` denotes a verification map and `r` a first-map redo. Number residuals
are signed electrons; all oracle differences are absolute. `none` means no
history solve. Rank 3/8 raw traces use the same CSV schema and reproduce the
same failure classifications.

`free-false-1-1e-8` (false=Gamma, true=2x2x2):

| Map | charge Ry | spin Ry | uniform Ry | dF Ha | dM mu_B | dmu Ha | occupation L1 e | number residual e | moment residual e | condition |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.3235082e-02 | 3.9155794e-03 | 3.7876821e-03 | 6.1965369e-03 | 4.4956251e-01 | 2.4030228e-02 | 4.4956251e-01 | 3.9968029e-14 | 0.0000000e+00 | none |
| 2 | 1.2369928e-03 | 2.2095545e-03 | 6.1471134e-03 | 2.7909261e-03 | 2.1199955e-01 | 2.1989185e-02 | 2.1199955e-01 | -2.5535130e-14 | 0.0000000e+00 | 1.0000000e+00 |
| 3 | 1.2310136e-05 | 5.1810333e-04 | 2.1606033e-03 | 4.3232237e-04 | 5.1120059e-02 | 1.5739322e-02 | 5.1120059e-02 | 3.2418512e-14 | 0.0000000e+00 | 5.1001306e+00 |
| 4 | 4.1632326e-07 | 3.8087689e-06 | 1.2713068e-05 | 9.4647901e-06 | 4.8388337e-03 | 5.0663546e-03 | 4.8388337e-03 | 9.8365760e-14 | 0.0000000e+00 | 7.6807899e+04 |
| 5 | 3.6799854e-08 | 2.0729494e-07 | 5.5836008e-07 | 8.0756648e-07 | 1.3057703e-03 | 1.7186156e-03 | 1.3057703e-03 | 5.9952043e-15 | 0.0000000e+00 | 4.9028961e+05 |
| 6 | 5.7346577e-10 | 1.6044258e-08 | 5.9976612e-08 | 1.0988022e-07 | 4.7113130e-04 | 6.5897632e-04 | 4.7113130e-04 | 5.2402527e-14 | 0.0000000e+00 | 1.5326664e+09 |
| 7 | 7.2051348e-12 | 1.2465443e-09 | 5.5181479e-09 | 1.1167167e-08 | 1.4883949e-04 | 2.1352597e-04 | 1.4883949e-04 | -7.4051876e-14 | 0.0000000e+00 | none |
| 7v | 8.0009271e-13 | 1.3834821e-10 | 6.1697229e-10 | 1.2874440e-09 | 5.0394693e-05 | 7.2903218e-05 | 5.0394693e-05 | -4.2077453e-14 | 0.0000000e+00 | none |

`free-true-1-1e-8` (false=Gamma, true=2x2x2):

| Map | charge Ry | spin Ry | uniform Ry | dF Ha | dM mu_B | dmu Ha | occupation L1 e | number residual e | moment residual e | condition |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.6931249e-02 | 1.3319731e-03 | 9.7308869e-05 | 6.3406180e-04 | 2.8944067e-01 | 7.1549080e-03 | 2.8944067e-01 | 3.3306691e-14 | 0.0000000e+00 | none |
| 2 | 1.4756016e-03 | 2.5103014e-04 | 1.6369230e-04 | 5.2263677e-04 | 2.5046186e-01 | 5.5643635e-03 | 2.5046186e-01 | -4.4630966e-14 | 0.0000000e+00 | 1.0000000e+00 |
| 3 | 7.3666932e-06 | 4.0388411e-05 | 1.9066726e-04 | 3.7799832e-04 | 1.9565727e-01 | 2.8708069e-03 | 1.9565727e-01 | 8.7707619e-14 | 0.0000000e+00 | 1.5610154e+04 |
| 4 | 1.5650241e-06 | 5.7685184e-06 | 1.6415290e-05 | 5.9312926e-06 | 1.5835623e-02 | 6.7887990e-04 | 1.5835623e-02 | -5.7398530e-14 | 0.0000000e+00 | 1.9161523e+04 |
| 5 | 1.9305577e-06 | 9.7098373e-06 | 4.8011577e-05 | 3.2492398e-05 | 4.2860517e-02 | 1.0714152e-03 | 4.2860517e-02 | -2.6090241e-14 | 0.0000000e+00 | 3.2068882e+05 |
| 6 | 1.8885315e-05 | 7.1726964e-05 | 3.7346945e-04 | 1.0071455e-04 | 5.7792276e-02 | 2.7012407e-03 | 5.7792276e-02 | -9.7699626e-14 | 0.0000000e+00 | 1.0048662e+06 |
| 7 | 4.6422470e-08 | 6.7169390e-08 | 3.6418560e-07 | 1.8628746e-07 | 2.9541355e-03 | 8.3593298e-05 | 2.9541355e-03 | -6.6613381e-15 | 0.0000000e+00 | 1.7341828e+08 |
| 8 | 2.0532277e-06 | 6.3517323e-06 | 3.3673278e-05 | 1.3112401e-05 | 2.3173170e-02 | 8.0972468e-04 | 2.3173170e-02 | -9.0927266e-14 | 0.0000000e+00 | 7.9409653e+08 |
| 9 | 8.1952641e-09 | 2.8696223e-08 | 1.5082269e-07 | 7.3313500e-08 | 1.8504234e-03 | 5.5734743e-05 | 1.8504234e-03 | 2.3980817e-14 | 0.0000000e+00 | 3.1317227e+10 |
| 10 | 6.4513042e-09 | 2.2457072e-08 | 1.1845131e-07 | 5.7447603e-08 | 1.6371697e-03 | 4.9355691e-05 | 1.6371697e-03 | 6.5503158e-14 | 0.0000000e+00 | 1.6208509e+10 |
| 11 | 3.2020901e-10 | 9.7087804e-10 | 5.1512396e-09 | 2.4829866e-09 | 3.3922495e-04 | 1.0232008e-05 | 3.3922495e-04 | 2.3314684e-14 | 0.0000000e+00 | none |
| 11v | 1.2096351e-10 | 2.9274760e-10 | 1.5135261e-09 | 7.4004038e-10 | 1.8503542e-04 | 5.5215213e-06 | 1.8503542e-04 | 9.1260333e-14 | 0.0000000e+00 | none |

`free-false-1-1e-12` (false=Gamma, true=2x2x2):

| Map | charge Ry | spin Ry | uniform Ry | dF Ha | dM mu_B | dmu Ha | occupation L1 e | number residual e | moment residual e | condition |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.3235082e-02 | 3.9155794e-03 | 3.7876821e-03 | 6.1965369e-03 | 4.4956251e-01 | 2.4030228e-02 | 4.4956251e-01 | 3.9968029e-14 | 0.0000000e+00 | none |
| 2 | 1.2369928e-03 | 2.2095545e-03 | 6.1471134e-03 | 2.7909261e-03 | 2.1199955e-01 | 2.1989185e-02 | 2.1199955e-01 | -2.5535130e-14 | 0.0000000e+00 | 1.0000000e+00 |
| 3 | 1.2310136e-05 | 5.1810333e-04 | 2.1606033e-03 | 4.3232237e-04 | 5.1120059e-02 | 1.5739322e-02 | 5.1120059e-02 | 3.2418512e-14 | 0.0000000e+00 | 5.1001306e+00 |
| 4 | 4.1632326e-07 | 3.8087689e-06 | 1.2713068e-05 | 9.4647901e-06 | 4.8388337e-03 | 5.0663546e-03 | 4.8388337e-03 | 9.8365760e-14 | 0.0000000e+00 | 7.6807899e+04 |
| 5 | 3.6799854e-08 | 2.0729494e-07 | 5.5836008e-07 | 8.0756648e-07 | 1.3057703e-03 | 1.7186156e-03 | 1.3057703e-03 | 5.9952043e-15 | 0.0000000e+00 | 4.9028961e+05 |
| 6 | 5.7346577e-10 | 1.6044258e-08 | 5.9976612e-08 | 1.0988022e-07 | 4.7113130e-04 | 6.5897632e-04 | 4.7113130e-04 | 5.2402527e-14 | 0.0000000e+00 | 1.5326664e+09 |
| 7 | 7.2051348e-12 | 1.2465443e-09 | 5.5181479e-09 | 1.1167167e-08 | 1.4883949e-04 | 2.1352597e-04 | 1.4883949e-04 | -7.4051876e-14 | 0.0000000e+00 | 6.3616801e+10 |
| 8 | 3.5912092e-13 | 6.7525053e-10 | 3.2240916e-09 | 6.6056366e-09 | 1.1436010e-04 | 1.6445106e-04 | 1.1436010e-04 | -5.9952043e-15 | 0.0000000e+00 | 2.7464826e+14 |
| 9 | 2.0368932e-14 | 1.2134984e-10 | 5.7664756e-10 | 1.2078024e-09 | 4.8807055e-05 | 7.0571067e-05 | 4.8807055e-05 | -1.7430501e-14 | 0.0000000e+00 | 3.3797604e+14 |
| 10 | 6.4070256e-15 | 9.0505600e-12 | 4.2129715e-11 | 8.9484392e-11 | 1.3270633e-05 | 1.9245216e-05 | 1.3270633e-05 | 2.8865799e-15 | 0.0000000e+00 | 1.5408175e+16 |

`free-true-1-1e-12` (false=Gamma, true=2x2x2):

| Map | charge Ry | spin Ry | uniform Ry | dF Ha | dM mu_B | dmu Ha | occupation L1 e | number residual e | moment residual e | condition |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.6931249e-02 | 1.3319731e-03 | 9.7308869e-05 | 6.3406180e-04 | 2.8944067e-01 | 7.1549080e-03 | 2.8944067e-01 | 3.3306691e-14 | 0.0000000e+00 | none |
| 2 | 1.4756016e-03 | 2.5103014e-04 | 1.6369230e-04 | 5.2263677e-04 | 2.5046186e-01 | 5.5643635e-03 | 2.5046186e-01 | -4.4630966e-14 | 0.0000000e+00 | 1.0000000e+00 |
| 3 | 7.3666932e-06 | 4.0388411e-05 | 1.9066726e-04 | 3.7799832e-04 | 1.9565727e-01 | 2.8708069e-03 | 1.9565727e-01 | 8.7707619e-14 | 0.0000000e+00 | 1.5610154e+04 |
| 4 | 1.5650241e-06 | 5.7685184e-06 | 1.6415290e-05 | 5.9312926e-06 | 1.5835623e-02 | 6.7887990e-04 | 1.5835623e-02 | -5.7398530e-14 | 0.0000000e+00 | 1.9161523e+04 |
| 5 | 1.9305577e-06 | 9.7098373e-06 | 4.8011577e-05 | 3.2492398e-05 | 4.2860517e-02 | 1.0714152e-03 | 4.2860517e-02 | -2.6090241e-14 | 0.0000000e+00 | 3.2068882e+05 |
| 6 | 1.8885315e-05 | 7.1726964e-05 | 3.7346945e-04 | 1.0071455e-04 | 5.7792276e-02 | 2.7012407e-03 | 5.7792276e-02 | -9.7699626e-14 | 0.0000000e+00 | 1.0048662e+06 |
| 7 | 4.6422470e-08 | 6.7169390e-08 | 3.6418560e-07 | 1.8628746e-07 | 2.9541355e-03 | 8.3593298e-05 | 2.9541355e-03 | -6.6613381e-15 | 0.0000000e+00 | 1.7341828e+08 |
| 8 | 2.0532277e-06 | 6.3517323e-06 | 3.3673278e-05 | 1.3112401e-05 | 2.3173170e-02 | 8.0972468e-04 | 2.3173170e-02 | -9.0927266e-14 | 0.0000000e+00 | 7.9409653e+08 |
| 9 | 8.1952641e-09 | 2.8696223e-08 | 1.5082269e-07 | 7.3313500e-08 | 1.8504234e-03 | 5.5734743e-05 | 1.8504234e-03 | 2.3980817e-14 | 0.0000000e+00 | 3.1317227e+10 |
| 10 | 6.4513042e-09 | 2.2457072e-08 | 1.1845131e-07 | 5.7447603e-08 | 1.6371697e-03 | 4.9355691e-05 | 1.6371697e-03 | 6.5503158e-14 | 0.0000000e+00 | 1.6208509e+10 |
| 11 | 3.2020901e-10 | 9.7087804e-10 | 5.1512396e-09 | 2.4829866e-09 | 3.3922495e-04 | 1.0232008e-05 | 3.3922495e-04 | 2.3314684e-14 | 0.0000000e+00 | 6.3329716e+12 |
| 12 | 8.4194946e-11 | 2.4951534e-10 | 1.3267481e-09 | 6.3901882e-10 | 1.7201978e-04 | 5.1870260e-06 | 1.7201978e-04 | -4.3520743e-14 | 0.0000000e+00 | 2.9077805e+13 |
| 13 | 5.7279610e-13 | 1.5843602e-12 | 8.4368893e-12 | 4.0659975e-12 | 1.3721568e-05 | 4.1260821e-07 | 1.3721568e-05 | 7.5939255e-14 | 0.0000000e+00 | 8.8573654e+14 |

`moment-true-1-1e-8` (false=Gamma, true=2x2x2):

| Map | charge Ry | spin Ry | uniform Ry | dF Ha | dM mu_B | dmu Ha | occupation L1 e | number residual e | moment residual e | condition |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.7312773e-02 | 1.0425507e-03 | 0.0000000e+00 | 2.7145670e-05 | 8.9928065e-14 | 4.3788586e-03 | 1.2750801e-02 | -2.4646951e-14 | -8.4932061e-14 | none |
| 2 | 1.6466514e-03 | 1.0518008e-04 | 0.0000000e+00 | 7.0911379e-06 | 4.5963233e-14 | 2.6235463e-03 | 4.4693890e-03 | 5.8175686e-14 | 5.1236793e-14 | 1.0000000e+00 |
| 3 | 6.3719273e-08 | 9.0788841e-08 | 0.0000000e+00 | 1.0360689e-08 | 3.8635761e-14 | 2.6819822e-05 | 4.5406926e-04 | -4.7517545e-14 | -3.6137759e-14 | 2.8666642e+04 |
| 4 | 2.0117388e-09 | 4.5868914e-09 | 0.0000000e+00 | 6.4545096e-10 | 4.8294702e-15 | 1.1961486e-06 | 1.4082987e-04 | -1.5099033e-14 | -2.2204460e-15 | none |
| 4v | 1.5940804e-12 | 3.0261548e-10 | 0.0000000e+00 | 4.7009618e-11 | 5.0515148e-14 | 7.6496600e-07 | 4.5284302e-05 | 2.0650148e-14 | -4.8849813e-14 | none |

</details>

### G4B candidate measurements (Mac, 2026-11-25)

The complete 48-comparison matrix was rerun in each arm, with the production
beta=0.7/history=8 and unchanged acceptance bounds. Each arm also reran all
12 Si2 FD comparisons with the same external UPF/hash stated above.

| Stiffness | Conditioning | Pass / 48 | Rejected runs | Remaining failures |
| --- | --- | ---: | ---: | --- |
| off | off | 33 | 6 | All standard free-moment cases; fixed-M standard mesh occupations; tight free-moment singular histories |
| on | off | 39 | 6 | Standard free-moment mesh M/occupations; tight free-moment singular histories |
| off | on | 45 | 0 | Fixed-M standard mesh occupations |
| on | on | **48** | **0** | **None in this matrix** |

The stiffness-only mesh maximum moment error is 1.7145382127e-4 mu_B;
conditioning-only fixed-M mesh occupation error remains 1.4082986520e-4
electrons. No third number-root switch was added, because both paths already
use the same root tolerance and occupation definition. The combination, not
either individual change, meets every registered bound.

Maxima over Gamma/mesh and ranks 1/3/8 with both candidates:

| Policy | threshold Ry | abs(dF) Ha | abs(dM) mu_B | abs(dmu) Ha | occupation L1 electrons |
| --- | ---: | ---: | ---: | ---: | ---: |
| scalar | 1e-08 | 4.5800058190e-11 | 0.0000000000e+00 | 8.5941571258e-07 | 7.8777769115e-07 |
| scalar | 1e-12 | 6.6752159356e-15 | 0.0000000000e+00 | 9.1950014824e-08 | 7.5100252097e-07 |
| fixed | 1e-08 | 1.9797810613e-10 | 1.9984014443e-15 | 0.0000000000e+00 | 0.0000000000e+00 |
| fixed | 1e-12 | 1.6253665081e-13 | 1.7763568394e-15 | 0.0000000000e+00 | 0.0000000000e+00 |
| free | 1e-08 | 1.8569937255e-09 | 8.6988030765e-05 | 8.7483935304e-05 | 8.6988030765e-05 |
| free | 1e-12 | 7.1997963147e-14 | 3.8134513081e-07 | 5.5365723445e-07 | 3.8134513122e-07 |
| moment | 1e-08 | 6.2009231083e-11 | 4.9626969201e-14 | 3.4210137784e-06 | 4.3497456097e-05 |
| moment | 1e-12 | 1.9394208461e-14 | 3.1363800446e-14 | 8.2458214146e-09 | 3.0991102276e-07 |

All Si2 FD arms pass 12/12, with maximum abs(dF) against strict legacy FD
3.2998670463e-10 Ha standard and 4.9737991503e-14 Ha tight. These are internal
legacy oracles; no new QE reference or speed claim is made.

With both candidates, the one-rank free-moment Gamma run drops 2 pairs at
standard and 3 at tight; mesh drops 4 and 7. Fixed-M Gamma/mesh each drop zero
pairs at standard and 2 at tight. For example Gamma map 6 removes two oldest
pairs at pre-drop conditions 1.48102953e9 and 2.70393441e8. Every event is
retained in `Iteration.history_conditioning` and optional trace CSV; there
is no silent restart. Largest retained one-rank free-Gamma tight condition
is 2.82427134e7. Standard Gamma's accepted-map stiffness term is
3.6385331062e-9 Ha; fixed-M mesh's is 8.6041613751e-11 Ha after one extra map.

The qualified paths are bit-identical with candidates off/on: all 24
scalar-FD/fixed-sector synthetic fingerprints and all 12 Si2 fingerprints
match byte for byte. The stream includes density/orbital f64 bits, eigenvalues,
occupations, energy ledger, chemical potentials, every scientific iteration,
inner-solver record, and verification record, excluding wall clocks only.
Closed-shell Gamma at ranks 1/3/8 matches the pre-candidate diagnostic build
and both final switch arms byte for byte, including scientific history and
orbitals. SHA-256 values:

| Ranks | SHA-256 |
| --- | --- |
| 1 | `9afed278586f803a836e19d46470789de0c6e337c7ad59474bef8f2cdaef729d` |
| 3 | `6f329ff3360630d73610dfbb0931de7c3f9d7626159b5315ec8aadc63a49ba5e` |
| 8 | `0cc15ffcf4d25a1ddfff91c9c4740a0a900d8b4bc463f276c5d9da6f11cee57f` |

Reproduce the passing arm (absolute trace/output paths must already exist):

```sh
BAANDI_REFSCF_SPIN_STIFFNESS=on BAANDI_REFSCF_HISTORY_CONDITIONING=on \
DFT_SI_UPF=/absolute/Si.pz-vbc.UPF \
cargo test --release -p dft-engine --test refscf g4 -- --include-ignored --nocapture --test-threads=1
```

This exits 0 (8 test functions), including the short-history diagnostic and
manufactured oldest-pair/zero-history test in addition to the full matrix.
Either single candidate and the both-off arm exit 101 on the explicitly
ignored acceptance test; each still runs all 48 comparisons and all Si2 cases.
`G4B_TRACE_DIR` additionally exports map-wise CSVs after return, and
`G4B_FINGERPRINT_DIR` exports the scientific streams. Logs, all four-arm traces,
command exit records and the machine-readable summary are under ignored
`target/tmp/g4b/`. Representative baseline traces above persist in this document.

The scope is one synthetic spin system plus scalar Si2 FD, simulated ranks,
pure Rust on Mac. Real MPI/MKL and multi-material qualification on dev08
remain the orchestrator's gate. Defaults remain off; ADR 0045 is proposed.
Dashboard/registry, performance ledger and generated README evidence are
outside this ticket's explicit owned files and are left for the orchestrator.

### Routing handoff and ownership boundary

The measured combination requires the CLI predicate in
`crates/dft-cli/src/run.rs::reference_spin_fd_fallback` to require all three
flags: `REFSCF_SPIN_FD_MOMENT`, `REFSCF_SPIN_STIFFNESS`, and
`REFSCF_HISTORY_CONDITIONING`. The ticket owns only the switch-count edit in
CLI tests, not `src/run.rs`. A concrete unified patch, including a focused
8-combination routing test and the existing end-to-end test's three-switch
input, is prepared at `target/tmp/g4b/routing.patch`. An explicit ownership
extension was requested; this out-of-scope patch has not been applied or
claimed as tested. Until it is applied, the existing development selector
alone still bypasses legacy routing. All three flags remain default off.


Verification actually executed: `cargo check --workspace --tests` exits 0;
`cargo test -p dft-engine --test refscf` exits 0 (79 passed, 20 ignored);
the LDA-PW energy-second-difference kernel unit test exits 0; the release
passing-arm G4 command above exits 0. `cargo fmt --check` and
`cargo clippy -p dft-engine -p dft-numerics --all-targets -- -D warnings`
exit 0 after fixing two lint-only iterator suggestions. The full
`CHECK_SCCACHE=0 BAANDI_TEST_TMP=/private/tmp/g4b-spin-moment CHECK_TEST_JOBS=2 scripts/check.sh`
exits 0: build/test/fmt/Clippy all pass, 984 passed, 0 failed, 213 ignored
across 98 binaries. The sandbox denies the optional `ps` load probe; the
script continues and runs every gate. No git metadata was written and all
changes are uncommitted. The routing patch remains a separate, unapplied
ownership handoff, so these checks do not qualify that CLI change.

## Symmetry

G6 adds `refscf::symmetry::Prepared` and the `run_prepared` /
`run_prepared_with` entries in `k_driver` and `electronic`. Preparation uses
`scf::PreparedSymmetry` and its existing `KOrbitMap`; representatives, ordering,
multiplicities, normalized weights and time reversal are unchanged. The CLI
reuses its cached preparation through `from_prepared`, which checks the exact
geometry, full mesh, tolerance and time-reversal setting. The reference wrapper
retains the geometry and rejects reuse after an atom/cell change. Explicit FFT
grids must preserve the selected density group and are used exactly, without
FFT-friendly rounding. `Config::prepare` owns this selection for both CLI
drivers; the direct reference API independently validates its supplied grid.
Automatic grids retain the existing `compatible_fft_grid` choice, selecting
2/3/5-smooth dimensions that preserve density translations and rotations.
Preflight records why the automatic dimensions were selected. The CLI result
`basis.fft_grid_overridden` is true only for user-specified dimensions, not
automatic symmetry padding. `crystal_density_group` selects the
full crystal group; otherwise the mesh-preserving unitary subgroup is used.
Antiunitary partners contribute only to k weights, never to the scalar density
projector. Collinear magnetization is a scalar; spatial operations do not swap
channels or rotate its quantization axis.

The full-stick distributed projector implements the serial expression

```
(P rho)(G) = sum_{(R,tau)} exp(-2 pi i G.tau) rho(R^T G) / |group|.
```

The density sphere is closed under every operation; missing rotated Miller
indices reject instead of being aliased. Coefficients outside the sphere are
zeroed, as in the serial projector. Ownership is the existing full-grid FFT
stick partition. A deterministic plan is built once per run: each destination
requests each needed source coefficient once, sorted by global FFT index within
each source rank. Each projection packs all channels into **one** equal-block
all-to-all for symmetry, including self blocks and empty owners. Each map now
performs two projections: the physical projection and an invariance diagnostic
on a reusable scratch copy. The diagnostic adds a second symmetry exchange;
it never replaces the physical density with the twice-projected result. FFT transpose
collectives are separate and unchanged. No full-valued density gather is added.
The phase/accumulation order is density-basis order then operation order,
independent of rank count; receive order never controls a floating-point sum.
Plan buffers and terms are reused. The engine wrapper uses the existing owned
field-transform API, which returns owned vectors; this does not claim that the
whole SCF map or its FFT boundary is allocation-free.

Projection follows the complete k-weighted band-density sum and precedes the
output functional audit, stopping residual and Broyden mixing. Fresh verification
maps take the identical path. Both spin channels share the same exchange.
Initial atomic or uniform channel densities are checked against their projected
values at the legacy pointwise allowance of 1e-12 electrons/Bohr^3, without
replacing the initial values. There is no additional normalization, clipping,
magnetization constraint or change to any convergence threshold.

`Run.symmetry` / `results.json.reference_scf.symmetry` records density group
order, irreducible k count, channel count and chronological map residuals.
With `rho` denoting the irreducible-k band sum before projection, each row has:

- `pre_projection_norm`: `||rho - P rho||`, in electrons/Bohr^(3/2). The existing
  `residual_norm` remains an alias for this value for compatibility.
- `pre_projection_relative_norm`: `||rho - P rho|| / ||rho||`, dimensionless
  (zero for zero density). This can remain nonzero at SCF convergence because
  the irreducible-k sum need not itself be invariant.
- `post_projection_invariance_norm`: `||P(P rho) - P rho||`, in
  electrons/Bohr^(3/2), expected at rounding for the invariant projected density.

Norms use `sqrt(Omega sum_G,spin |rho(G)|^2)` with forward-normalized coefficients
and Parseval, including discarded out-of-sphere components in the pre-projection
difference and denominator. The diagnostic uses a scratch copy of `P rho` and
a batched reduction of the three squared norms. Physical density, energy,
mixing, stopping criteria and verification are unchanged. These are projection
diagnostics, not convergence gates. Rows identify the SCF iteration and verification
phase; first-step redo produces another row with the same iteration. Startup checks
are excluded. Off-path records omit the JSON block.

G5 delivery records the same operations in the force model identity.
`delivered::force_state_prepared` validates that identity and exposes the group
to the unchanged force and frozen-stress projectors (atom permutations and
Cartesian vector/tensor rotations). A mismatched/omitted group fails closed.
Symmetry-enabled Gamma uses the complex prepared entry; symmetry-off Gamma
still dispatches to the existing real driver with unchanged arithmetic.
Real-only ladder/RMM restrictions and the existing collinear-FD qualification
switch remain unchanged. No QE source was read, no QE reference was generated,
and no speed/default promotion follows from this feature.

Validation commands (external UPF is never committed):

```sh
cargo test -p dft-physics symmetry
cargo test -p dft-engine --test refscf
DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf g6_si444 -- --ignored --nocapture
cargo test --release -p dft-engine --test refscf g6_spin_fd -- --ignored --nocapture
```

The projector oracle covers two channels, fractional translations, idempotence,
serial agreement and empty owners. The SCF matrix compares diamond Si2 4x4x4
and a skew deformation, full mesh versus reduced reference versus reduced
legacy, at standard 1e-8 Ry and tight 1e-12 Ry estimator thresholds and simulated
ranks 1/3/8. Acceptance energy differences remain 5e-9 / 1e-10 Ha. Separate
synthetic fixed-spin and scalar-FD comparisons exercise channel projection.
These are internal and legacy-code oracles, not new QE or native-MPI evidence.

### G6 local evidence and limits

Mac pure Rust, base `841682a6afa0ca668ecbcc71663ea6bd92b13a86`, simulated
collectives. Si2 uses 4 Ha cutoff, 5 bands, an unshifted 4x4x4 mesh and the
external `Si.pz-vbc.UPF` SHA-256
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
The diamond group has 48 unitary operations and 8 irreducible k points; its
skew deformation has 2 operations and 36 points. Each row covers ranks 1/3/8.
These are new engine-to-engine measurements, not QE references or timings.

| Comparison | Standard maximum difference (Ha) | Tight maximum difference (Ha) |
| --- | ---: | ---: |
| Si diamond reduced vs full reference | 4.2944e-10 | 3.1974e-14 |
| Si diamond reduced vs strict reduced legacy | 4.5784e-10 | 6.7502e-14 |
| Skew Si reduced vs full reference | 7.9936e-15 | 5.3291e-15 |
| Skew Si reduced vs strict reduced legacy | 1.0373e-10 | 1.8652e-14 |
| Synthetic scalar FD vs strict reduced legacy, free energy | 1.2418e-11 | 7.5079e-14 |
| Synthetic fixed spin vs strict reduced legacy, energy | 2.2880e-10 | 4.9127e-14 |

Maximum diamond reduced/full density differences are 2.0124e-7 and 1.1513e-8
electrons/Bohr^3 at the standard/tight thresholds. Tight symmetrized derivative
differences against legacy are at most 4.0960e-8 Ha/Bohr (force component)
and 3.6403e-10 Ha/Bohr^3 (stress component), below the unchanged G5 1e-6
limits. All delivered comparisons converge and pass fresh verification.

The synthetic spin/FD fixture is a one-electron Gaussian nonlocal potential,
width 1.8 Bohr, skew 5-Bohr cell, atom at fractional (1/4,1/4,1/4), 2.5 Ha
cutoff, grid 16^3 and shifted 2^3 mesh. Scalar FD uses LDA-PZ, theta=0.01 Ha
and 12 bands; fixed spin uses LDA-PW, sector (1,0) and 4 bands. Collinear
preparation rejects time reversal, exactly as legacy. Both channels are also
checked with nonzero data against the serial projector in the physics test.
This does not qualify the previously unqualified collinear-FD moment policies.

An initial narrow synthetic atomic-start fixture failed invariant validation:
clipping at insufficient cutoff generates coefficients outside the density
sphere. The fast test now explicitly proves that both legacy and reference
reject it at the same 1e-12 allowance. Passing fixtures use a resolved broad
atomic density; no numerical tolerance or production start policy changed.
Initial oracle setup also needed the existing strict-profile `Verification::Off`
setting and a compatible explicit FFT grid; these are fixture corrections.

The G5 optional cutoff-correction scan in `stress::energy_at` does not forward
symmetry and lies outside G6's owned files. Symmetry-enabled CLI stress therefore
returns the symmetrized analytic frozen-basis tensor, `pulay=null`, and an explicit
`pulay_unavailable_reason`; it never evaluates an irreducible mesh without its
projector. Symmetry-off stress retains its three cutoff scans. Extending those
scans requires the separate symmetry-forwarding change in `stress.rs`.
Native MPI/MKL, dev08 performance and independent external symmetry oracles
remain unmeasured. Shared README evidence generation, dashboard/session registry
and performance ledgers remain outside the ticket's ownership.

Symmetry-off Gamma preservation was checked against the actual base source,
with all modified tracked Rust files restored from `git show HEAD:<path>` for
the baseline build and then restored byte-for-byte to this implementation.
Both builds ran `g3_gamma_base_fingerprint` in the same debug profile. Complete
scientific fingerprint byte streams (density, energies, eigenvalues, histories,
verification and orbitals; clocks excluded) are identical, not merely their
hashes. Rank 1/3/8 lengths are 100801/100763/100729 bytes and SHA-256 values are
respectively `9afed278586f803a836e19d46470789de0c6e337c7ad59474bef8f2cdaef729d`,
`6f329ff3360630d73610dfbb0931de7c3f9d7626159b5315ec8aadc63a49ba5e`, and
`0cc15ffcf4d25a1ddfff91c9c4740a0a900d8b4bc463f276c5d9da6f11cee57f`.
Logs and fingerprint streams are in the ignored `target/g6-*.log` and
`target/tmp/g6-{base,current}/` directories. Changes are uncommitted.
The existing memory census is not extended to count symmetry exchange buffers
or group-term metadata; it must not be treated as a complete symmetry-run peak
or RSS measurement.

Final local gates all exited 0: `cargo check --workspace --tests`, physics
`symmetry` tests (2 passed), engine `--test refscf` (79 passed, 22 ignored),
CLI `reference` tests (15 passed, 10 ignored), `cargo fmt --check`, and Clippy
for engine/physics/CLI with `--all-targets -- -D warnings`. The release `g6_`
matrix with `--include-ignored --nocapture --test-threads=1` passed all 3 tests,
including the external-UPF cases above. The full `scripts/check.sh` also
exited 0: PASS build/test/fmt/clippy, 98 test binaries, 984 passed, 0 failed,
215 ignored. Its environment was `CHECK_SCCACHE=0 CHECK_TEST_JOBS=2
CARGO_NET_OFFLINE=true BAANDI_TEST_TMP=/private/tmp/g6-symmetry-gates`;
sandbox-denied optional `ps`/priority probes did not prevent any gate.
`python3 -I scripts/tests/test_readme_figures.py` passed all 9 checks, and
`git diff --check` passed. These checks do not add G6 results to the shared
README generator; that integration remains with the owner of those files.



### G4C distributed spin PBE

`gga::distributed::SpinPbe` owns the polarized xcx functionals, interleaved
rho/sigma, Cartesian gradient multipliers, six gradient fields, a channel and
flux buffer, derivative coefficients, and two response potentials. All are
prepared once per run with the owned slab/stick layout and reused by input,
output-audit and verification maps. Existing FFT `Fields` and xcx return-vector
allocations remain at the field boundary; there are no per-orbital allocations
or full-grid payload gathers in this adapter. Empty owners join all transforms
and scalar error reductions, skipping only the local xcx call.

The actual combined samples `n_up+c/2` and `n_down+c/2` are forward transformed
separately. This includes the same immutable real core as serial `SpinXc`,
without clipping or normalization. G2B's cached core Fourier field remains
used by the scalar initialization; spin maps forward the combined samples to
match the serial input arithmetic. Cartesian derivatives remove each even-axis
Nyquist Miller component before applying the reciprocal cell matrix. Three
invariants (up-up, up-down, down-down) feed polarized GgaXPbe/GgaCPbe. Channel
flux is `2*e_sigma_ss*grad(n_s+c/2) + e_sigma_ud*grad(n_other+c/2)`;
axis-by-axis divergence subtraction follows the serial adapter's ordering.
Each combined channel must be finite and nonnegative. The serial spin total
threshold (1e-10 electrons/Bohr^3) screens both channels and all three sigma
inputs together; gradients use unscreened samples. A fully polarized nonvacuum
point is not screened, and divergence can remain nonzero at a screened point.

Only XC sees the core. Hartree, charge, occupations, mixing and delivered
channel densities remain valence-only. The complete potentials, including
cross-channel divergence, enter the valence double-counting integral; the
separate `integral c/2*(V_up+V_down)` is diagnostic in the NLCC ledger. The
existing scalar Hartree assembly uses its local LDA-PW branch for spin PBE,
whose scalar XC result is discarded, avoiding redundant unpolarized PBE FFTs.
The G4B opt-in stiffness remains the documented local LDA-PW curvature proxy
also for PBE; it is not asserted to be the PBE Hessian or an error bound.
All existing scalar and LDA dispatch/arithmetic remains unchanged.

Reproduction: `cargo test -p dft-physics --test spin_pbe_distributed --
--nocapture` checks random positive channels, optional signed core, odd/even
skew grids, reusable storage, and ranks 1/3/8 (including empty owners), against
serial spin PBE with a 1e-13 relative gate. The release-only engine matrix is
`cargo test --release -p dft-engine --test refscf g4c_spin_pbe_delivered --
--ignored --nocapture`: synthetic one-electron open shell, core on/off, fixed
sector and fixed-moment FD, all three G4B switches on, ranks 1/3/8,
thresholds 1e-8/1e-12 Ry and legacy energy gates 5e-9/1e-10 Ha. These are Mac
synthetic/serial comparisons, not real-MPI, real-material or QE qualification.

G4C local results (2026-11-25, Mac, no tolerance changes): the random-grid
potential/energy/integral comparison passes with maximum relative difference
6.952422915745373e-16. Vacuum screening, exactly polarized channels, errors on
one rank, empty owners and recovery after a rejected call also pass. All 24
pre/post scientific fingerprints of `g4_synthetic_matrix` are byte-identical:
scalar LDA-PZ FD and fixed-spin LDA-PW, Gamma/k222, ranks 1/3/8, both thresholds.

Architect review qualifies fixed-sector and fixed-moment FD PBE with and
without NLCC. `g4c_spin_pbe_delivered` covers these 24 comparisons. The separate
ignored `g4d_free_moment_pbe_acceptance` retains all 12 free-moment comparisons
and the original accuracy assertions; three standard no-core cases are expected
to fail. Run it with `cargo test --release -p dft-engine --test refscf
g4d_free_moment_pbe_acceptance -- --ignored --nocapture`. The explicit engine
entry remains available for this diagnostic; production CLI routing below
prevents claiming this policy is qualified.

| Policy (maximum over core on/off and ranks 1/3/8) | Standard abs(delta F), Ha | Tight abs(delta F), Ha |
| --- | ---: | ---: |
| Fixed sector (F=E) | 1.0069276246138514e-9 | 2.353672812205332e-14 |
| Fixed-moment FD | 1.2032827512165056e-10 | 8.926193117986259e-14 |
| Free-moment FD | **9.206498157443566e-7 (FAIL)** | 1.7263968032921184e-14 |

The no-core free-moment legacy oracle has M=1.0000000000000473 mu_B and
F=-0.23459382746784757 Ha. At standard threshold, rank one delivers
M=0.9999446771843783 mu_B and F=-0.234592906818031827 Ha despite passing the
existing fixed-point verification. Tight threshold delivers
F=-0.234593827467864835 Ha. With NLCC, free-moment FD passes both thresholds
(maximum standard difference 5.537099667662915e-10 Ha). This exposes a remaining
near-full-polarization stopping/verification qualification issue with the
G4B local LDA curvature proxy; it does not establish a proven cause or a safe
replacement estimator. No stopping-rule adjustment, forced occupation or
material-specific branch is introduced in G4C.

The CLI routes free-moment collinear FD + PBE to legacy regardless of all
three G4B switches, with the exact reason
`reference driver: free-moment collinear FD with PBE not yet qualified, G4D`.
The reason is recorded at `reference_driver_fallback_reason`,
`execution.reference_driver_fallback_reason`, and the routing switch manifest.
Fixed-sector PBE remains reference-capable; fixed-moment PBE retains the G4B
three-switch gate used for qualification. Free-moment LDA-PW routing is unchanged.
Tests cover all eight G4B switch combinations and emitted legacy result metadata.
The obsolete PBE rejection assertion now checks convergence, while polarized
LDA-PZ still rejects. CLI limitation metadata reflects the G4D boundary.

### G4D starting diagnosis: no-core free-moment PBE

Measured with the shared G4C fixture above, Gamma, one simulated rank,
threshold 1e-8 Ry, theta=0.01 Ha, all G4B switches on. The strict legacy oracle
is F=-0.234593827467847571 Ha, M=1.00000000000004730 mu_B and minority occupation
N_down=3.09057801998514300e-79 electrons. Here N_down is the sum of the 12
minority-band occupations (the sole k point has weight one), not a single band.
Errors below are signed map output minus the oracle. F is the Mermin energy;
charge/spin/uniform terms are Ry, stiffness is Ha. The augmented estimator is
charge + spin + uniform + 2*stiffness, in Ry. Trace capture is opt-in and prints
only after the run, including verification; it does not change stopping rules.

| Map | Charge (Ry) | Spin (Ry) | Uniform (Ry) | Stiffness (Ha) | delta F (Ha) | delta M (mu_B) | delta N_down (e) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 1 | 1.3235081547e-2 | 3.9155855220e-3 | 3.7876952238e-3 | 2.5780599862e-3 | 7.3564795323e-3 | -4.5607997736e-1 | 2.2803998868e-1 |
| 2 | 1.2378305580e-3 | 2.1753968769e-3 | 6.0785868720e-3 | 2.5833494784e-3 | 3.6190665559e-3 | -2.2025401117e-1 | 1.1012700558e-1 |
| 3 | 1.3745852296e-5 | 5.7210948066e-4 | 2.3538749687e-3 | 9.9004650211e-4 | 6.3437969940e-4 | -4.9470134861e-2 | 2.4735067430e-2 |
| 4 | 1.6989383255e-6 | 1.1034603071e-5 | 2.3618040028e-5 | 4.0683161479e-5 | 4.3988174758e-5 | -6.7121104264e-3 | 3.3560552132e-3 |
| 5 | 6.2287939400e-7 | 3.1041094226e-6 | 4.5203564108e-6 | 1.4922775423e-5 | 2.1626730087e-5 | -4.1051943499e-3 | 2.0525971749e-3 |
| 6 | 3.7666646453e-10 | 1.6253431336e-7 | 9.0139350192e-7 | 3.7532555987e-6 | 1.1152301916e-5 | -2.5995944151e-3 | 1.2997972075e-3 |
| 7 | 2.3822277013e-11 | 4.0522157406e-8 | 2.1264347720e-7 | 1.3248913100e-6 | 6.6845210760e-6 | -1.7868719499e-3 | 8.9343597488e-4 |
| 8 | 4.3225903715e-11 | 1.7095746321e-9 | 4.8577228066e-9 | 7.9809844720e-8 | 2.4407404897e-6 | -6.9607938778e-4 | 3.4803969391e-4 |
| 9 | 4.2573259839e-12 | 1.0322776393e-10 | 3.0049343553e-10 | 9.3671881529e-9 | 1.5706219113e-6 | -3.1307660644e-4 | 1.5653830322e-4 |
| 10 | 5.0680109715e-12 | 2.1071153004e-11 | 4.1315659311e-11 | 3.9989240082e-9 | 9.2064981574e-7 | -5.5322815669e-5 | 2.7661407799e-5 |
| 10 verify | 1.7855355597e-12 | 1.7402515194e-11 | 3.4841595080e-11 | 4.2709047420e-9 | 7.5776810851e-7 | -3.1928575017e-5 | 1.5964287490e-5 |

Map 10's augmented estimator is 8.0653028397e-09 Ry, below 1e-8 Ry, while
its true free-energy error is 9.2064981574e-7 Ha. The fresh verification map
has augmented estimator 8.5958391298e-09 Ry and true error 7.5776810851e-7 Ha;
it passes, but delivery correctly retains map 10's candidate. The tail reaches
near-full polarization with a small residual-based estimate while a minority
occupation remains. This is diagnostic evidence of an underestimated delivered
error, not proof of a specific cause or permission to replace the estimator.
G4D must address qualification; G4C changes no tolerances, occupation floors,
stopping rules or verification criteria.

After the architect's qualification/routing review, all four requested checks
pass: workspace/test compilation; dft-physics (101 passed, 18 ignored);
dft-engine `refscf` (79 passed, 22 ignored); and fmt/three-crate clippy with
`-D warnings`. The qualified release matrix passes all 24 comparisons. The
separate G4D release acceptance test exits 101 as expected, retaining its three
standard no-core accuracy failures. With `BAANDI_TEST_TMP` and `TMPDIR` set to
`/private/tmp/g4c-check`, `CHECK_NICE=0 CHECK_SCCACHE=0 scripts/check.sh fast
dft-physics dft-engine dft-cli` exits 0: build/test/fmt/clippy all pass, 59
binaries, 632 passed, 0 failed, 167 ignored. The sandbox still denies the
non-gating `ps` load probe; cache tests now pass outside the checkout.
No stopping rules or tolerances were changed. Changes remain uncommitted.
### G6B grid and diagnostic regression

The explicit-grid replacement was in `Config::prepare`: after validating an
explicit grid it unconditionally called the automatic chooser. That rounded
28^3 to 32^3 because 28 is not 2/3/5-smooth. Selection now branches on input
provenance before invoking that chooser. Automatic dimensions retain their
previous selection and numerical mappings, while their override flag is false.
The reference driver's direct prepared API already honored compatible explicit
grids and retains its independent validation.

Local synthetic diamond regression (1 Ha cutoff, shifted 2^3 mesh, three maps):
pre-projection relative norms are 0.3414538869678375, 0.2510191207624144 and
0.2418444925004731; post-projection invariance norms are respectively
5.16519e-17, 4.08951e-17 and 2.49180e-17 electrons/Bohr^(3/2).
The automatic-grid final energy is unchanged at -0.08240083504027451 Ha;
the before/after density-bit FNV fingerprint is `4a2263a7c90dba78`.
These are internal synthetic checks, not Si.pz-vbc or QE reference results.
The direct API regression also checks explicit 28^3 and rejects 25^3 with
three simulated ranks. CLI coverage checks both drivers, mesh/crystal density
groups and Gamma/shifted meshes at the 30 Ry diamond grid-selection boundary;
actual output tests use a lower-cutoff synthetic potential.

Reproduce with `cargo test -p dft-cli --test grid_symmetry` and
`cargo test -p dft-engine --test refscf symmetry -- --nocapture`.
The added diagnostic uses one preallocated local reciprocal scratch vector per
channel and one extra projector application per map; native MPI timing is not
measured here. README evidence generation and dashboard/session registries
remain outside this ticket's owned files.

G6B final local gates exited 0: workspace check with tests, CLI tests (167
passed, 45 ignored), engine symmetry tests (3 passed, 2 ignored), formatting,
and Clippy for CLI/engine/basis with warnings denied. Full `scripts/check.sh`
passed build/test/fmt/clippy: 100 binaries, 992 passed, 0 failed, 217 ignored.
The full gate used `CARGO_NET_OFFLINE=true CHECK_SCCACHE=0 CHECK_TEST_JOBS=2
BAANDI_TEST_TMP=/private/tmp/g6b-tests`; its optional sandbox-denied `ps` probe
did not prevent verification. README generation made no tracked changes and
all 9 freshness checks passed. Logs and the before/after automatic-grid
fingerprints are retained locally under `target/tmp/g6b-validation/`.

## Non-self-consistent bands

R4A-REF-NSCF supports scalar and collinear fixed-density spectra. The user
approved the visibility-only extension to `electronic_map.rs`: `Fields`,
`Fields::hxc` and `fields` are visible to sibling reference modules. Their
arithmetic and the SCF callers are unchanged.

R4A-REF-NSCF adds `refscf::nscf::{nscf, nscf_with}` and
`refstart::start_kpoints_from_density`. The engine API accepts `StartInput`,
replicated density arrays, an ordered k list, an explicit total band count and
a positive residual tolerance in Hartree. One array denotes total unpolarized
density; two arrays denote up/down collinear densities, not total/magnetization.
The explicit-communicator variant permits the same deterministic collective
harness as the reference SCF driver. All inputs must be identical across ranks.

The supplied real-space grid is retained. There is no interpolation, clipping,
charge normalization, density update, occupation solve or mixing. The checks
retain the legacy finite-density lower bound of -1e-14 electrons/Bohr^3 and
charge tolerance of 1e-8 times the valence electron count. Unpolarized inputs
require positive even integer electron count and at least N/2 + 4 bands. For
collinear input the minimum is four above the ceiling of the larger integrated
channel charge (charge-check-sized roundoff is removed before the ceiling).
PZ has no polarized adapter and rejects collinear input. Too few plane waves or
an incompatible supplied grid reject, rather than changing the basis or grid.

Initialization builds ionic/core fields and uses the supplied valence density.
It shares the existing atomic-orbital Ritz construction and kinetic completion
with the k SCF start; it does not call an atomic-density SCF initializer. Atomic
orbitals are still required by this reference start contract. The scalar
Hartree/XC potential uses `electronic_map::fields` (delegating to `map::fields`
for scalar input), including the existing NLCC and spectral PBE core-gradient
treatment. Spin input uses the same per-channel LDA-PW/PBE fields as SCF, with
the Hartree potential from the total valence density and the nonmagnetic core
split equally between channels. The immutable channel potentials are prepared
once, then
redistributed to each k Hamiltonian. Gamma also uses complex orbitals. All ranks
process each k point in sequence; there are no k-point pools.

Every band is marked strict in `ref_davidson_complex`, bypassing the SCF loose
empty-band criterion. Ritz-value changes are not the acceptance certificate.
The inner Ritz-change threshold is the minimum positive normal `f64`, so
only exactly unchanged Ritz values can end that inner check; a zero threshold
is invalid in the shared solver. This is not an energy-accuracy threshold.
After each solver call NSCF freshly applies H to the returned orbitals, reduces
squared residual norms over owned rows, and accepts only when every finite
residual is at most the requested tolerance. Otherwise it restarts all bands
from those vectors, for at most 64 calls. Exhaustion reports channel/k context
and the maximum residual in Hartree; it never returns unconverged spectra.
The result retains original k coordinates/weights, eigenvalues and residuals
in Ha, and distributed eigenvectors with their global row indices. It is not
an SCF result or a fixed-point certificate. CLI routing is R4b's responsibility.

Validation commands:

```sh
cargo test -p dft-engine --test refscf_nscf
cargo test --release -p dft-engine --test refscf_nscf synthetic_rank_matrix -- --nocapture
DFT_SI_UPF=/absolute/path/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf_nscf silicon_path -- --ignored --nocapture
cargo test --release -p dft-engine --test refscf_nscf spin_channels_match_serial_dense_oracle -- --ignored --nocapture
```

The scalar matrix compares identical nonuniform fixed densities with legacy
NSCF at 1e-8 Ha residual certification on both sides and a 1e-8 Ha eigenvalue
bound. The legacy wrapper retains SCF buffer-band stopping and rejects the
ordinary low-band-count fixture with an unconverged empty-band residual. The
oracle therefore requests its complete plane-wave basis separately at each k
and compares the requested leading bands. Each legacy call still passes its
own all-band 1e-8 Ha guard. Its single-point weight is one; k weights do not
enter the fixed-density Hamiltonian. It includes Gamma/X/L, skew two-species local/nonlocal projectors, PBE
NLCC, and external Si2, using simulated communicator sizes 1/3/8. The spin
oracle uses the serial collinear Hartree/XC adapters and full dense Hermitian
diagonalization because legacy `nscf` has no collinear entry. Synthetic fixtures
are not material or QE evidence. Simulated ranks do not establish real MPI
transport/scaling. No SCF thresholds, defaults or arithmetic are changed.

Local evidence (2026-11-25, simulated ranks, not QE or MPI timing):

| Fixture | Ranks | Maximum eigenvalue difference (Ha) | Maximum NSCF residual (Ha) |
| --- | --- | ---: | ---: |
| Cubic synthetic, no atomic charge | 1 | 2.232e-14 | 5.002e-9 |
| Skew two-species | 1/3/8 | 2.109e-15 | 6.418e-9 |
| PBE NLCC synthetic | 1/3/8 | 2.498e-15 | 9.897e-9 |
| Si2 Gamma/X/L, 3 Ha cutoff | 1/3/8 | 2.942e-15 | 3.236e-9 |
| Spin LDA-PW NLCC, serial spin-field/dense oracle | 1/3/8 | 2.387e-15 | 9.080e-9 |
| Spin PBE NLCC, serial spin-field/dense oracle | 1/3/8 | 2.998e-15 | 7.503e-9 |

The Si file is external `Si.pz-vbc.UPF`, SHA-256
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
The spin matrix covers two unequal, nonuniform channel densities integrating
to 1.25 and 0.75 electrons, a skew cell, nonlocal projectors, NLCC and three
k points. The serial dense oracle independently certifies its retained bands
at 1e-8 Ha; the reference solve certifies every requested band at the same
bound. The input channel arrays remain unchanged. These are correctness
observations, not a performance or material-accuracy qualification.

Final gate results after the approved spin integration: `cargo check
--workspace --tests`, the default NSCF suite (3 passed, 2 ignored), the existing
`refscf` suite (84 passed, 24 ignored), formatting and engine clippy all exited
0. The complete release NSCF matrix, including external Si2 and both spin
functionals, passed all 5 tests with no ignored tests. The full gate
`BAANDI_TEST_TMP=/private/tmp/r4a-check CHECK_SCCACHE=0 CHECK_NICE=0
scripts/check.sh` exited 0: 107 binaries, 1,017 passed, 0 failed, 221 ignored;
build/fmt/clippy also passed. External scratch resolves the earlier three
CLI cache-fixture failures caused by worktree-local scratch. The optional
sandbox-denied `ps` load probe did not affect the gate. All nine README
freshness checks passed. Logs are retained under `target/` as
`nscf-all-approved.log`, `check-full-approved.log`,
`check-workspace-approved.log`, `nscf-approved.log`, `refscf-approved.log` and
`clippy-approved.log`.

## Stepper (R8a)

`refscf::stepper::Stepper` is a synchronous, single-owner continuation of the
Gamma, complex-k, and electronic (fixed spin / Fermi–Dirac) drivers. Construct it
with `driver::stepper[_with]`, `k_driver::stepper_with`,
`electronic::stepper_with`, or the adapter-facing `api::stepper(input, request)`.
Call `initialize()` once, then `step()` until `is_done()`, then `finish()`.
`result()` borrows a delivered result. `run_to_end()` is precisely that sequence;
all existing run entry points and `api::solve` use it. A caller may also call
`run_to_end()` after any successful partial sequence.

Initialization validates and builds the original numerical workspace, stopping
before the first SCF map. Calls before initialization, repeated initialization,
and calls after completion fail without advancing numerics. A numerical error
is terminal and releases the continuation; no partial result is fabricated.
`finish()` rejects incomplete or failed runs and never executes numerical work.
An exhausted run retains its original delivered nonconverged result.

Each report describes one logical SCF map, including its inner eigensolver
retries. First-iteration SCF redo is a separate `Retry` map. Verification is a
separate `Verification` map, including a failed verification. Its iteration
index names the ordinary candidate being checked. Reports carry the estimator
in Ha, dr2 and ethr in Ry, optional fully computed energy in Ha, and H columns
for that map only. Lean unaudited maps carry `None`, never a fabricated energy.
A cached failed verification consumed as the next ordinary map is marked
`reused = true` and has zero new H work. Coarse ladder maps have stage `coarse`;
fine maps have stage `fine`, with the original per-stage iteration numbering.

Mixing order is unchanged: an ordinary candidate that needs verification yields
with `safe_point = false`; the following verification step includes the pass
or failure decision and any recovery mixing of that candidate. A discarded
first map also yields an unsafe boundary before its retry. Other successful
boundaries are safe only after the original mixing/state handoff, or after
terminal delivery when no next eigensolve exists. A mixer failure is not a safe
point. Observers may read scalar reports at every return; interventions must
wait for a safe point. No map is repeated to implement stepping, and verification
never replaces the delivered candidate.

Implementation uses Rust's compiler-generated pinned continuation to retain
loop locals and their existing ownership and drop order. There is no async
executor, background thread, polling loop, wakeup, callback or transport in the
engine. Two allocations at construction hold the continuation and a single
scalar report slot. Yielding and polling that slot allocate nothing. Existing
map, audit, history and delivery allocations remain; this is not a claim that
the pre-existing SCF maps are allocation-free. No density, orbital or HX clone
was added, and no collective was inserted or reordered.

The continuation is thread-affine because existing configuration, memory and
solver diagnostic scopes are thread-local. Keep its configuration scope alive,
and finish or drop one initialized stepper before initializing another on the
same thread; overlapping initialization is rejected without advancing either
run. Independent threads/ranks can own independent steppers. MPI callers
must initialize, step, finish, and abandon a run collectively in the same order;
there is no asynchronous rank cancellation. Caller pauses are included in the
existing wall clocks, so timing gates use uninterrupted `run_to_end()`.

Checkpoint encoding/restoration and policy mutation are R8b/R8c work, not R8a
APIs. The compiler continuation itself is not a serialization format: those
stages must expose explicit versioned numerical state rather than serialize
opaque future bytes. Python/MCP adapters must drive the continuation on its
owning worker thread. Derivative and ionic-seed delivery in `api::stepper`
occurs at the terminal step just as it does in `api::solve`.

Validation: `cargo test -p dft-engine --test refscf_stepper` exercises the short
lifecycle/exhaustion case. The complete synthetic Gamma, shifted 2x2x2 k mesh,
FD and fixed-spin matrix, verify off/on, simulated ranks 1/3/8, is release-only:

```sh
cargo test --release -p dft-engine --test refscf_stepper -- --include-ignored --nocapture
```

The test compares every density/orbital bit as well as energy bits, prints
FNV-1a fingerprints over little-endian f64 bytes, and reconciles per-step H
columns and verification counts with the original run records. Simulated ranks
establish collective ordering and same-rank identity, not real-MPI timing or
cross-rank bit identity. The separate `reference_identity` test and release
`bench_scf.py --require-identical` retain their pre-refactor reference baselines.

Local validation (2026-11-25, pure Rust macOS/aarch64, Rust 1.99.1): the
run-only fixture was independently compiled from an archive of starting commit
`81763a54f574b4a63b625303a91667626f89013b`. All 24 converged
fixture/rank/verification combinations and the one exhausted Gamma case retain
its energy/free-energy bits and density/orbital fingerprints. The parent
fingerprints are stored in `tests/refscf_stepper.rs`; the stepwise comparison
also checks all array bits, internal energy and entropy, and H-column totals.
The release stepper suite passes all five tests, including failed verification,
whole-map reuse, coarse/fine boundaries, partial facade execution, rejection of
interleaved sessions, and error/drop cleanup.

Final local gates all exit 0: workspace check; debug stepper (3 passed,
2 release-only ignored); existing refscf (84 passed, 24 ignored);
`reference_identity` (1 passed); full `scripts/check.sh` (109 binaries,
1,025 passed, 226 ignored, no failures; build/fmt/clippy pass); engine clippy;
and README generation/freshness (9 tests, generated output unchanged).
The release benchmark with the default `baseline-reference.json`, one thread,
three repetitions each of Gamma/k222s/k222 reports seven iterations per case,
zero energy difference in Ha and identical density hashes. These are identity
checks, not QE comparisons or Mac performance acceptance.

Two validation setup errors were corrected without changing numerical gates:
the first release test compile selected stale parent-archive library artifacts
from the shared target (exit 101; only the engine release artifacts were cleaned
before rebuilding), and the first benchmark invocation inherited test-only
`BAANDI_TEST_TMP` (exit 3; removed for the benchmark's switch-manifest check).
Final invocations passed. The sandbox also denies the gate's optional `ps`
load probe; all four verification gates still ran and passed. Real MPI8 identity
and wall acceptance remain with the master session's
[`r8a-dev08-request.md`](../performance/lanes/RSCF/r8a-dev08-request.md).

## Runtime control

R8C-RUNTIME currently provides a lifecycle foundation in `refscf::runtime`;
CLI dispatch and L7 migration are **not complete**. The adapter owns the R8a
continuation and exposes scalar snapshots with completed-map count, last map,
state version and safe-point status. Trusted host pause/resume preserves that
continuation. Cancel releases it at a safe point and retains the scalar partial
record, without fabricating a delivered density or converged result. Controls
requested between a candidate and its verification, or before a retry, wait for
a safe boundary. Terminal numerical delivery wins over pending controls.
Transport authentication and approval are not supplied by this lifecycle API.
Collective callers must provide identical controls to every rank.

`telemetry::reference_step_json` uses the existing version 1.0 event envelope,
with one `iteration_finished` event per logical map, including verification
and retry maps. Its payload names the kind, stage, estimator in Ha, dr2 and
ethr in Ry, optional energy in Ha, H columns, reuse and safe-point flags.
Unavailable legacy density-residual and charge diagnostics are null; dr2 is
not a real-space density norm. Formatting occurs after a step returns.

The remaining policy transaction requires a mailbox at the stepper boundary
and application hooks in `driver.rs`, `k_driver.rs` and `electronic.rs`.
Those four existing files are outside R8C's current OWNED_FILES. The R8a public
API cannot update the profile, Broyden mixer or immutable verification option
held inside its pinned continuation. A wrapper must not claim to apply a policy
by updating only metadata, restarting from scratch or repeating maps.
Approval to extend ownership is pending. The intended transaction validates
beta/history/threshold/verification, resets the mixer for beta/history changes,
retains density/orbitals, and applies only at a safe boundary. Unsupported
legacy knobs must reject explicitly. L7 stays until this integration is ready.

The historical self-healing evidence in
`tests/reference/runtime/generated/results.json` uses a `MaxIterations` budget
extension. That scenario is not applicable to R8c's closed policy allow-list;
this does not establish reference-driver numerical-instability recovery.
The new `reference-results.json` records only the tested lifecycle and event
serialization scope, separately from the historical README evidence.

Reproduce the synthetic one-/three-simulated-rank identity and pause/resume
checks (verification off/on), event capture, and partial cancellation checks:

```sh
R8C_WRITE_EVIDENCE=1 cargo test -p dft-cli --test reference_runtime -- --nocapture
```

This is same-build BAANDI identity, not a QE comparison or native MPI timing.

Local validation of this partial implementation: the three engine lifecycle
unit tests pass. The new CLI integration target cannot compile because the
starting HEAD calls nonexistent `l4` and `l5` functions in `driver_select.rs`
(lines 75 and 78); those diagnostic blocks are outside the L7-only allowance.
Workspace check and CLI tests exit 101 for that pre-existing error.
An isolated Rust test harness compiled the exact integration-test source and
`reference_step_json` function against the engine library, excluding the broken
CLI crate. Both tests pass: 16 rank/mode/verification records converge with
identical same-rank energy bits and density SHA-256, plus cancellation retains
its partial scalar record. The new evidence includes the actual isolated
compile command, source hashes and this limitation. It is not a passing CLI
acceptance gate. The local reproducer is `target/tmp/r8c-isolated.py`; once the
routing build is repaired, use the Cargo command above to regenerate evidence
through the intended integration target.

Final local gates for the partial patch: engine-only Clippy with warnings
denied and formatting pass (exit 0). Full `scripts/check.sh`, with
`CHECK_SCCACHE=0` and `BAANDI_TEST_TMP=/private/tmp/r8c-runtime-tests`, exits 1:
build/test/clippy fail on the same missing routing predicates; fmt passes.
The combined CLI/engine Clippy command exits 101. No CLI acceptance or L7
migration is claimed. All changes remain uncommitted pending scope resolution.
