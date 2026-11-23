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
| Restart | Rotate V and HV by C: 4rmn FLOPs. Then fresh S/H: 4rn^2 FLOPs; factor S; two small refresh-transform GEMMs: 4n^3 FLOPs. No H call or EVD. | Rotation payload `16*(rm+mn+rn)` bytes; copies into V/HV: 32rn bytes. Fresh Gram payload `16*(2rn+n^2)` bytes and two sums totaling 16n^2 bytes, plus 48n^2 reset-store bytes, 16n^2 HR-to-H copy bytes and adapter-owned small-matrix/factor work. | Required physical rotations: two, not three. The note's ideal diagonal-H/unit-S reset eliminates the measured refresh, but that was already falsified by a 1.0154e-11 Ha degenerate-fixture drift. Preserve measured refresh and its factor; do not repeat that rejection. |
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
