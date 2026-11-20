> **Reopened, 2026-11-20:** the orchestrator supplied one-core STREAM
> 12.97 GB/s and explicitly requested serial fusion. The historical closing
> decision below is superseded by [dev08-request.md](dev08-request.md).
> The conservative 11.466 GB model now predicts 0.884 s, above the design
> screen. This is a testable bandwidth hypothesis, not adoption evidence.
> DGEMM 64–108 GFLOP/s does not establish FFT throughput. No MPI lever or
> acceptance rule changed. Implementation preserves ascending band additions
> and is default-off, standard/unregistered-process/one-core only.

# LANE-GLUE Phase 2 diagnosis

Decision: **recommend closing this lane under the requested minimum gain and
current ownership**. No reviewed candidate has a defensible expected saving
of at least 0.2 s on Si64 MPI8 or 0.8 s on Si64 one core. This is not a proof
that density cannot be improved. In particular, serial density fusion remains
an unpriced hypothesis; the supplied aggregates do not establish its gain.
No optimization, numerical default, or convergence contract is changed.

## Provenance and denominator

All measured seconds below are **user-relayed dev08** values: source `e947a35`,
instrumented, socket 0, CPUs 0–7 for MPI8. The supplied one-core values remain
in a separate model. No Mac profiling or timing supports this diagnosis.
Binary/UPF/input hashes, raw rank results, repeat spread, and the new GLUE
subspan/work counters were not included in this relay. They are not invented.
The user supplied 13 SCF maps plus two verification maps for the Si64 serial
density normalization. Whether coarse-stage density is additionally included
in that aggregate must be checked against the original bucket data.

Local analysis base is `8de1711`. `git diff e947a35 HEAD` is empty for
`grid_dist.rs`, `band_density.rs`, `scf.rs`, `energy.rs`,
`distributed_real_gamma.rs` and `switches.rs`; these inspected paths therefore
match the profile source. No QE source was accessed.

MPI8, each number an envelope from the same relayed run family:

| cell | GLUE s | share of external wall | density s | density/energy/potential/mixing share of GLUE |
| --- | ---: | ---: | ---: | --- |
| Si64 | 1.1690 | 13.70% | 0.7331 | 62.7 / 24.5 / 10.5 / 2.3% |
| Si32 | 0.2801 | 14.36% | 0.1336 | 47.7 / 28.7 / 19.5 / 4.1% |
| Si16 | 0.1279 | 19.00% | 0.0428 | 33.5 / 34.1 / 26.1 / 6.3% |

One core:

| cell | GLUE s | share of external wall | density s | density/energy/potential/mixing share of GLUE |
| --- | ---: | ---: | ---: | --- |
| Si64 | 8.2797 | 20.32% | 6.3032 | 76.1 / 12.0 / 9.4 / 2.6% |
| Si32 | 1.8018 | 22.91% | 1.1107 | 61.6 / 18.1 / 16.3 / 3.9% |

Shares are rounded and need not sum to exactly 100%. `m2.density` aliases the
density span: count it once. `m2.xc` includes work already inside potential
and energy. `m2.nonlocal_energy` is inside energy. `m2.solve_prepare` 0.1529 s
and IO 0.2689 s (Si64 MPI8) are context, not added to GLUE. Per-component rank
maxima do not necessarily occur on one rank; their sum is not an exclusive
critical-path decomposition.

The user relayed the same-window MPI8 QE Davidson `sum_band=0.34`,
`v_of_rho=0.03`, `mix_rho=0.03` s. The directional differences are about
+0.39, +0.09 and -0.003 s. Energy has no directly matching QE row; its proposed
+0.24 s is a residual-accounting estimate, not a measured equivalent routine.
The roughly +0.7 s GLUE-type gap out of the stated 6.1 s gap is context only.
No new QE reference, precision claim for 0.03 s timers, or cross-model ratio
is asserted without the original inputs/thresholds/hashes.

## Density: required transforms, accumulation, and MPI ownership

Source chain: `scf.rs::m2_density` / `band_density.rs::compute_keep` ->
`distributed_real_gamma.rs::band_planes` / `fused.rs::SplitPairDensity`.
At integer occupations, zero-weight tails are skipped. Two real Gamma
orbitals occupy the real/imaginary parts of one complex inverse transform.
For an occupied contiguous prefix B, there are ceil(B/2) inverse transforms
per density map, not B transforms and not a forward/inverse pair. General
complex k-point orbitals cannot use this Gamma pairing.

Si64's documented 30 Ry fixture (`benchmarks/cases/si64-gamma.toml`, historical
size in `docs/performance/mpi-scaling-gap.md`) has grid 72^3 = 373,248 points,
256 valence electrons, 128 occupied bands and normally 132 solved bands.
Thus the fixture model has 64 inverse transforms per map, 960 over 15 maps.
These sizes are a **fixture model**, not extracted e947a35 counters. A coarse
ladder grid or a changed grid override must be accounted separately.

Already adopted and not proposed again:

- Ledger A26 fused/pruned density, zero occupations skipped.
- DENSITY_PAIR and DENSITY_REDUCE defaults on: Gamma pairing, ordered band sum.
- MPI-FS-DENSITY-FUSED and subsequent adoption: standard M2 accumulates
  directly on owned planes; no bands-by-plane contribution materialization.
- M2-OWNED-DEFAULT: owned density/potential persist through maps; full density
  gathers are consumer boundaries, not the ordinary-map default.
- Lazy orbitals avoid the full orbital allgather/unpack before density/audit;
  audit kinetic energy is obtained from local rows. REUSE-PSI retains the new
  orbitals' real-space grids to accelerate the **next Hamiltonian apply**.

MPI8 is a plane/stick distribution, not a band partition followed by a
sum of eight density arrays. Each rank accumulates all occupied bands at its
owned points in ascending band order. The inverse FFT has a batched transpose
and possible mirror-row exchange (mirror ownership can make that exchange
empty). There is **no cross-rank density sum** on the eligible owned path.
If operator and grid-plan plane axes differ, PlaneMap redistributes planes.
The gathered fallback uses an allgatherv permutation, not a sum. Later audit
and charge norms gather/reduce canonical partials. No communication-volume
lever is proposed in this task.

### Arithmetic and traffic floor inputs

For a full complex 3-D FFT the conventional work proxy is
`5 N (log2 nx + log2 ny + log2 nz)` real flops per inverse transform.
For 72^3 this is 34.54 million per pair, 33.16 billion over 960 pairs.
This is an algorithmic estimate, not a hardware flop counter or a rigorous
lower bound on every FFT implementation. For a pruned stick transform the
corresponding proxy is `5[S nx log2 nx + Z nx ny log2 ny + N log2 nz]`,
with actual active-stick count S and active-z count Z from the plan. Never
assign the full-grid estimate to the pruned MPI path without its sizes.

A relevant backend distinction: `fused.rs::FusedPlan` selects
`GridFusedPlan<FftPlan3d>` under **mkl**, whereas the pure build selects
`OrderedFusedPlan`. The MKL serial adapter uses full-grid transforms and
scratch; the word "fused" does not imply that it has the pure backend's
pruned stages. The serial build/features must be confirmed from the raw
profile, although MKL is the requested dev08 configuration. FFT-lane work
already owns investigation of the FFT implementation itself.

At each grid point and occupied band, density needs the real orbital value,
its square, occupation/volume weighting and ordered accumulation (roughly
four scalar operations with the current expression including division).
It must write at least N output doubles. The current serial path additionally
clears and writes the per-band contribution arrays, later reads them for the
ordered fold, and may write 64 complex kept-psi grids. In the fixture model:

| logical work, 15 fine-grid maps | amount |
| --- | ---: |
| three FFT-axis read/write passes, 32 bytes/complex point/axis | 34.399 GB |
| kept-psi writes, if all 64 pairs fit the memory gate | 5.733 GB |
| serial split contribution clear/read-modify-write/fold read | 22.932 GB |
| estimated net reduction from replacing split storage by direct ordered rho accumulation | 11.466 GB |

GB is decimal. These are logical traffic models, **not DRAM counters**:
cache reuse, scratch copies, write allocation and vendor internal passes alter
actual traffic. The net reduction assumes 32 bytes/band/point for split
(clear 8, += read/write 16, later read 8) versus 16 for direct rho +=.
It excludes an additional final density pass and cannot be priced as an
end-to-end saving simply by dividing by nominal socket bandwidth.

With the supplied denominator, 6.3032/15 = **0.420213 s/map**, or
**6.565833 ms/paired inverse transform** including the non-FFT work. Assigning
all that time to the dense FFT proxy would give 5.261 GFLOP/s; this is a
normalization, not measured FFT throughput. To reach 0.8 s requires removing
12.692% of serial density. Serial fusion's 11.466 GB traffic model reaches
0.8 s only at an effective removed-pass throughput of 14.333 GB/s or lower,
assuming the reduction translates entirely to critical-path savings. The new
relay supplies one-core STREAM 12.97 GB/s, pricing the conservative saving at
0.884054 s. This qualifies the candidate for implementation; removed-pass
throughput and actual wall gain remain to be measured by crossed A/B.
For MPI8 the density target requires 27.281% of density, and this split-array
traffic has already been removed by the adopted fused path.

A roofline time lower-bound model is `max(F/R_peak, Q/B_peak)` for each
non-overlapping kernel, with actual backend work F and relevant traffic Q.
No measured/qualified FFT throughput, cache traffic or detailed GLUE subspan
distribution was returned. STREAM bandwidth is now available, but is not a
measurement of this kernel's removed DRAM traffic. **The numerical excess above a
hardware time floor is undetermined**, not 6.30 s minus an invented FFT time.
The large serial share is consistent with 64 transforms/map plus band-grid
traffic while most non-density glue is O(N), and with fast vendor solver
algebra shrinking its denominator. It is not evidence of duplicate transforms.

## Energy: what is already shared and what can be skipped

`energy.rs::prepare_shared` performs one output-density forward FFT using
`hartree::energy_and_density_reciprocal`; it does **not** inverse-transform
an unused output Hartree potential. It reuses those coefficients for local
ionic energy, except when clamping changes the input (raw-density fallback).
The coefficients are also retained for the standard Coulomb stop estimator.
`grid_dist::audit` likewise performs one forward transform, per-stick Hartree
contraction, output XC energy/integral, and canonical partial sums.

Every map currently records kinetic, ionic local and its zero mode, nonlocal,
Hartree, XC, constant Ewald, band sum, double counting, Harris-Foulkes,
band-form energy and functional-difference bound. Ewald and ionic potential
construction are setup work; EWALD/IONIC fast paths are already adopted,
not recalculated wholesale in each energy map. The per-map ionic local
*integral* changes with output density and is not that construction cost.

The standard stopping estimator primarily needs the input/output Coulomb
residual; not every explicit audit scalar is necessary to that estimator
alone. Nevertheless the current session consumes finite explicit energy,
records its delta/history, and uses fresh verification energy/functionals.
Harris-Foulkes depends on input fields; direct total and Bregman remainders
also depend on output fields. Reports and strict-reference checks carry real
values, not stale/zero placeholders. Deferring audit therefore requires a
session/history contract and a final/verification materialization policy; it
is not a safe branch in `grid_dist::audit` alone.

The existing ENERGY-IDENTITY already tests substituting the band-sum identity
for explicit nonlocal work and auditing delivery. ENERGY-IDENTITY-DEV08
(2026-11-13) was not adopted: +0.3–1.2% wall, unchanged density/SCF. The new
relay does not demonstrate why its former overhead disappears. Projection
reuse is also already implemented as LEAN_PROJ_CACHE; it is not a new GLUE
candidate and is outside this task's GEMM-shape scope.

Si64 MPI8 nonlocal energy is 0.1161 s: even eliminating it entirely cannot
meet 0.2 s. Whole energy is 0.2865 s; the target would remove 69.808% of it.
Serial energy is 0.9897 s; the target would remove 80.833%. No supplied
subspan assigns that much time to dispensable reporting. One must not
subtract two independently maximized rank timers to claim an exact exclusive
remainder. Required output transforms/verification cannot be discarded to
manufacture that saving.

## Potential and the 37.4% m2.xc asymmetry

Per map the input potential needs a forward density FFT, reciprocal Coulomb
kernel/contraction, inverse Hartree FFT, pointwise XC and effective-potential
assembly. Output audit needs a different forward FFT and XC evaluation.
The density transform of orbitals cannot replace a transform of their summed
squares: `FFT(|psi|^2)` is nonlinear in psi (a reciprocal convolution), not the
orbital transform already computed. Reusing output coefficients for the
*next mixed input* is invalid except on exact identity paths; normally
rho_next differs from rho_out. Existing stop-coefficient reuse is already on.

In `grid_dist.rs`, **both `potential_field` and `audit` charge Stage::Xc**.
That timer contains FFTs, packing, collectives, canonical reductions and XC;
`m2.rs::xc` additionally uses the label on its fallback path. It is not a pure
pointwise XC timer. The correct pure computation counter is
`run.counters["glue_lane.grid.xc.ns"]`, with `grid.xc_points` and rank identity.
Those counters were not supplied. The 0.1559 s/37.4% figure cannot determine
whether skew is local XC work, FFT ownership, collective arrival, cold/fallback
calls or differing bucket coverage. The precise cause remains unresolved.

For an eligible 72-plane grid and 8 ranks, each rank owns 9 planes: simple
point-count imbalance is zero on that fine grid. Coarse/other grids may have
ceil/floor plane counts, PZ branches depend on local density, and ranks can
arrive at collectives at different times. Equal point count does not imply
equal elapsed time. These are hypotheses, not measured causes. Even interpreting
37.4% as `(max-min)/max` (definition absent from relay) gives only 0.058307 s
spread in the 0.1559 s envelope, not a removable 37% of GLUE. If the reported
0.1559 s is the full same-window M2 Xc envelope, deleting all its work would
still be below 0.2 s. If it is a narrower bucket, that ceiling must be revised
from raw bucket data before claiming a gain.

Pure potential alone is 0.1226 s MPI8 and 0.7749 s serial, each below the
requested threshold even at an impossible 100% elimination. XC imbalance
fixes cannot plausibly remove the whole potential plus energy audit without
also removing required physics. No load-balancing change is selected.

## Candidate decision

| candidate | predicted scope / ceiling, not measured gain | decision |
| --- | --- | --- |
| Share audit output FFT/XC with verification input potential | optimistic uniform-map estimate (not a measured bound), two potentials out of 15: 0.01635 s MPI8 / 0.10332 s serial; actual saving smaller because inverse Hartree remains | below threshold; mixed ordinary inputs defeat general reuse |
| Skip verification density build using prior density/kept psi | nominal two-map serial ceiling 0.84043 s, barely above target | invalid: verification solves new orbitals; reusing prior output would suppress the fixed-point residual being tested |
| Serial direct density accumulation instead of per-band arrays | 11.466 GB modeled removed traffic/15 maps; >=0.8 s needs <=14.333 GB/s effective removed-pass rate | chosen after new STREAM relay: expected 0.884 s; explicitly authorized serial band_density/FFT scope; default-off implementation, MPI already fused |
| Defer explicit energy audit except final/verification | whole-phase ceiling 0.2865 s MPI8 / 0.9897 s serial; needs 70%/81% removal, no demonstrated dispensable subspan that large | no qualified expectation; ENERGY_IDENTITY previously lost; current history/verification contracts must remain truthful |
| Reduce pointwise XC imbalance | proxy spread 0.05831 s if stated asymmetry is spread/max | below MPI target, cause not established, no serial imbalance to remove |

Mixing is deliberately excluded (0.0268 s MPI8, 0.2119 s serial), consistent
with the user and `mixing-survey.md`: no universal mixer default emerged from
2,244 historical runs and trajectory changes can regress accuracy/convergence.
The survey's old standard-path certification status is historical, not the
current ADR 0043 contract. No parameter tuning or guard experiment is reopened.

## Result and limits

The STREAM follow-up supersedes the initial closing recommendation. Implemented
`SERIAL_DENSITY_FUSED` as a conservative ordered sink, preserving FFT plans,
Gamma pairing, each band's arithmetic, kept psi and band addition order.
Strict reference and registered/compiled MPI execution remain on the original
path; only a standard one-core unregistered process can enter the new sink.
A small MPI-backend compile-time constant prevents the MPI rank-one CLI (which
does not register its communicator) from accidentally entering this lane.
No numerical allocation or logging was introduced inside transform loops.

The 0.8 s design screen is not an adoption minimum. Any resolved positive
same-model gain counts after correctness, with A/A first and |mean|/SE>=3.
See [dev08-request.md](dev08-request.md) for commands and withdrawal criteria,
and [fusion-validation.md](fusion-validation.md) for actual Mac bit/counter
evidence. The prior [closing-note.md](closing-note.md) remains historical.
[phase2-floors.json](phase2-floors.json) includes the STREAM-priced analytic
expectation, not an observed performance result. No Mac timings support it.
