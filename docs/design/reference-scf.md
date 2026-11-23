# Reference SCF path (issue #25): design

Status: design draft by the orchestrator, 2026-11-23. User decisions of the same day: all lane/A-B work stopped; a NEW SCF driver replaces the legacy one; numerical kernels are reused only where they are good, and bold rewrites are allowed wherever they promise speed (the old code stays in git history). Implementation goes through Codex tickets; QE conventions reach implementers only through the Claude note `docs/research/qe-study/qe-density-mixing-spec.md` and the published literature.

## Why a new driver
The dev08 decomposition (docs/performance/requests/si64-gap-attribution-2026-11-20.md, architecture brief addendum) shows the ratio to QE as per-map cost x map count x non-SCF overhead: si32 1.51 x 0.83 x 1.66, si64 1.56 x 1.00 x 1.49, si128 1.35 x 1.93 x 1.31 (A/D 3.586 measured), si216 1.21 x 1.50 x 1.47. The per-map cost is closing with size; the map count (si128: 29 maps vs QE 15, maps 3-13 diverge with the standard error growing 272x while no control reacts) and the overhead outside the SCF loop (coarse ladder 18-24 % of wall at every size, verification 9-11 %) are the levers. The legacy driver (`ScfSession::advance`, about 1100 lines, plus the cutoff ladder and verification scheduling, more than 130 switches) is where both live. A new driver that follows a published, QE-equivalent fixed-point scheme lets us (1) compare every iteration with QE term by term and (2) add BAANDI's own accelerations back one at a time, each with a measured gain.

## Scope of the new path
Real-Gamma, integer occupations, spin-unpolarized, norm-conserving, LDA/PBE: the cases of the benchmark suite (Si8...Si216, skew, slab, elongated, MgO Gamma, synthetic PPs). Everything else (k points, spin, smearing, DFT+U, forces/stress outputs) keeps using the legacy driver until the new one covers it; the CLI chooses by profile.

## Architecture
A new module `crates/dft-engine/src/refscf/` (no dependency on `session.rs`, `cutoff_ladder.rs` or the switch registry beyond one profile switch):

1. `profile`: default-off execution profile `BAANDI_REFERENCE_SCF=on` (input `[scf] driver = "reference"`), resolved once, recorded in results.json (`scf.driver`). The CLI branch sits where `run.rs` chooses between `run_ladder` and `run_outcome`; output schema stays the same (`ScfOutcome`), plus a `reference_scf` block per iteration.
2. `density`: the mixed quantity in reciprocal space exactly as defined in the note (half-sphere storage for real Gamma, G=0 handling, which G shells are mixed and what happens to the rest, normalization, units), distributed by G vectors under MPI (owned G slabs from the existing grid plan; no full-grid replication).
3. `metric`: the Hartree-like inner product used for the residual norm (the SCF accuracy estimate) and for the history overlap matrix, with its exact weight in |G|^2, the G=0 rule and the half-sphere weights; one reduction per inner product, batched.
4. `broyden`: the modified Broyden/Anderson update (Johnson, PRB 38, 12807 (1988)) with the startup step, beta, history length, insertion/eviction order and the linear solve, as specified in the note; deterministic reductions (fixed order) so results are rank-count independent within rounding.
5. `driver`: the loop. Input n_1 = atomic superposition at the FINAL cutoff (no coarse ladder in the reference arm). For each iteration k: inner threshold from the QE-equivalent rule (initial value, update from the previous accuracy estimate and the electron count, floor), one fixed-point map `ScfProblem::step_session(n_k, warm, threshold)` (the existing distributed Hamiltonian/Davidson/density kernels), accuracy estimate dr2 = metric(n_out - n_in), the first-iteration redo rule, stop when dr2 is below the threshold that the standard accuracy class maps to (decision 3 of 2026-11-11: match QE's conv_thr 1e-8 Ry equivalent), otherwise mix. After the stop: BAANDI's independent fixed-point verification (certificate), counted in the delivered wall.
6. `record`: per iteration: dr2, threshold, Davidson iterations and H columns, beta and history occupancy, step norm, energy terms, wall per phase; plus run totals. This is the comparison surface with QE (estimated scf accuracy, ethr, avg # of iterations).

## Patches (issue #25)
- Patch 1 (runnable density update): modules 1-4 and the driver skeleton; a reproducible test that calls the production update routine on (a) synthetic arrays with closed-form answers and (b) QE reference arrays (input density, output density, history) generated on dev08 from a recorded QE 7.6 Si128 input; compare the next density, total charge and low-G components in common units against a documented tolerance; report the first failing term if any. Deliver: commit, CLI command, comparison table.
- Patch 2 (common start convergence): same atomic start at the final cutoff for both codes; compare output densities first (a mismatch points at potential, occupations or the inner solve), then next inputs (a mismatch points at mixing); show whether Si128's early excursion disappears and the map count at equal accuracy (QE's 15 is a comparator, never a stopping rule). Record the reference-style stop and BAANDI's certificate separately, both charged to the delivered wall.
- Patch 3 (speed): dev08 pure MPI8 matched runs, Si128 with Si64 and Si216 controls, A/A floor first, crossed A/B, effect >= 3x floor and |mean|/SE >= 3, uninstrumented external wall including startup, certificate and output; skewed/slab and non-Si PP coverage; strict 1e-10 / standard 5e-9 Ha before any default change. Then add BAANDI accelerations back one at a time (kept orbitals first, then the ladder only if it pays at equal accuracy), and rewrite kernels where the per-iteration profile shows waste (per-iteration non-column cost was 99 ms vs QE 69 ms at si64; Gram/rotate/update 2.53 s vs QE 0.78 s of non-eigen algebra).

## Rules
No QE source to implementers; equations and conventions from the note and the literature only, the note's INFERRED items verified against QE reference data. Same parallel model on both sides. No constant tuned to silicon. The legacy driver stays the default until the new path passes the gates on the generality set; then the old driver and its switches are removed in one step.

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
