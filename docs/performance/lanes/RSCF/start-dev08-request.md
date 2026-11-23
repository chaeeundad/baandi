# R-SCF-START: dev08 first-Davidson boundary and RSS request

Status: pending master-session integration/measurement. Base d97d864; obtain
reviewed candidate commit from master (this worktree leaves changes uncommitted).
No remote action, QE run or Mac timing claim was made by this ticket.
Implementation/audit: [Start](../../../design/reference-scf.md#start).

## Boundary and prerequisites

Use the R-SCF-1 host adapter to call `refstart::start` from parsed input, before
constructing any legacy `ScfProblem`. Return state is already atomic-density,
first-potential and initial-Ritz ready at the final cutoff. The first Davidson
must consume its operator/owned rows directly, with H block width <= bands.
There is no CLI switch for this function yet. Do not benchmark a wrapper that
first constructs a legacy problem, nor call a missing switch and label it B.
The master owns the runnable adapter/instrumentation outside this ticket.

A is the legacy start **at the same final cutoff**, atomic starting density and
atomic wavefunctions, coarse ladder disabled. Preserve its merged defaults,
including M2_OWNED_INPUT=off, in both arms. If an A2 run with that switch on is
wanted, keep it separately labelled; never fold it into A/A noise. B is the new
start. Match atoms/coordinates, PP bytes and SHA-256, 15 Ha (30 Ry) cutoff,
explicit FFT dimensions, bands (2*N_atoms+4), XC, radial conventions, threshold,
seed, workers, binding, binary/compiler/MKL settings. Certificate stays on for
any completed SCF runs. No BAANDI-thread/QE-MPI comparisons, even informational.

Instrument host boundaries (exclusive, maximum over ranks):

1. external process launch -> parsed input ready, with MPI initialization and
   binary load separately identified if the host already measures them;
2. parsed input ready -> **entry to first Davidson** (includes all problem,
   grids, projectors, Ewald, atomic density, potential, atomic Q/HQ/Ritz and
   distributed operator/input/workspace setup);
3. external launch -> first-Davidson boundary; report host/test-harness entry
   separately, never pretend it includes process-load time;
4. first Davidson and eventual SCF/certificate/output costs separately if run.

`start` does not construct the Davidson workspace: include equivalent workspace
construction in both timed host boundaries. First-potential/atomic Ritz are
already included. A library-only `Instant(start)` is not the launch boundary.
Collect sampled /proc/<pid>/status VmRSS and VmHWM for **every rank**, maximum
and sum, at parsed input, after fields, atomic span peak, first Davidson and
process end. State sampling interval. VmHWM includes parser/binary/allocator
memory; principal-array counters do not. Retain rank-local raw logs outside the
numerical hot loops. A root-only `/usr/bin/time` launcher RSS is insufficient.

## Jobs and statistics

One case per job, sequentially, no builds or other benchmarks overlapping:

| job | cell | geometry replication | ranks | bands | purpose |
| --- | --- | --- | --- | --- | --- |
| 1 | Si64 | conventional 2x2x2 | 8 | 132 | time baseline and A/A floor |
| 2 | Si128 | conventional 4x2x2 | 8 | 260 | memory and scaling control |
| 3 | Si216 | conventional 3x3x3 | 8 | 436 | largest startup/RSS control |

Reuse master-approved geometry/input manifests and explicit FFT grids; do not
substitute guessed dimensions. The Mac correctness fixture displacement is
not automatically the headline geometry. Record all input and UPF hashes.
Use the approved MPI8 one-thread/rank binding on physical cores 0-7, same socket,
OMP_NUM_THREADS=1, MKL_NUM_THREADS=1; record mpirun version and binding output.
Pure and MKL are distinct result tables if both are requested.

For each job first run A/A (at least five crossed pairs) to measure the floor;
then at least ten crossed A/B pairs, alternating AB/BA, extending for a close
result. Warm binary/page cache equally; retain cold-launch observations as a
separate table. Report raw pairs, medians, spread, paired difference, SE,
|mean|/SE, and comparison to 3x A/A floor. The adoption rule has no arbitrary
minimum percentage: resolved positive gains may be adopted, but correctness or
resolved regressions prohibit adoption. If noise prevents resolution, report
inconclusive and extend, without switching metrics after seeing results.

Record source/binary hashes, Rust and backend versions, compilation commands
and exit codes, load/uptime before and after each job, every environment switch,
rank allocation and source files for all raw evidence. No remote installs.

## Correctness and expected counter movement

Before timing: run the release Mac-equivalent correctness matrix on dev08 and
one real-MPI start check at each rank count 1/2/3/8, same final-cutoff problem.
Compare normalized charge and exact N_el/Omega zero mode; max density/potential
error; global Gram of the owned start; Ritz energies and degenerate subspace
agreement; subsequent standard and strict energies/certificates. Keep strict
1e-10 Ha / standard 5e-9 Ha adoption classes and the 1e-12 Ha internal
unique-solution comparison. Add a skew/slab and non-Si synthetic PP control
before any default adoption. The new API currently rejects PBE and NLCC.

B must show: no coarse calls; no serial Hamiltonian/BandDensity construction;
no full Q/HQ/initial complex or full ORC start materialization; one final basis
and density basis; one distributed operator; owned potential epoch 0; retained
orbital bytes exactly 8*L_r*B; dense-field bytes/indices explicitly distinguished.
Report `StartMemory` per rank plus the true RSS observations. The full-indexed
8F potential, O(F) metadata, communication buffers and GEMM panel memory still
exist. Do not report the principal envelope as a measured peak allocation.

Required return: per-case A/A floor, all raw A/B pairs, first-Davidson phase
closure, every rank's RSS/HWM and counters, accuracy/certificate results,
resolved gain/regression verdict with limits, commands/version/hash manifest.
Only the master may then decide driver/default adoption. This patch changes
neither defaults nor any QE speed-comparison headline.
