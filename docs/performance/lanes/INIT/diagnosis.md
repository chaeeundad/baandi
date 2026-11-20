# LANE-INIT Phase 2 diagnosis and decision

2026-11-20. Source base `c0f2143769f0876931538e0e9d24c5e0c8202de0`.
Measurements below are **dev08, relayed by the user**, instrumented runs,
physical CPUs 0–7 on socket 0. Raw files, input/binary hashes, repeat count
and distributions were not supplied; they remain on dev08. Inclusive spans
are not summed. No Mac time is used here.

| model and cell | init / external wall, s | init share |
| --- | ---: | ---: |
| MPI8 si64 | 0.8398 / 9.3551 | 8.98% |
| one-core si64 | 2.5504 / 42.9552 | 5.94% |
| MPI8 si32 | 0.2909 / 1.9130 | 15.20% |
| one-core si32 | 0.5704 / 8.0192 | 7.11% |
| MPI8 si16 | 0.1520 / 0.6601 | 23.03% |

The `problem` timer is the constructor, not all startup. Its reported shares
of init are 31% / 12.5% for si64 MPI8 / one-core and 45% / 26% for si32.
Atomic span/Ritz creation is deferred until the first cold map's potential
exists. Atomic start dominates: si64 MPI8 0.52558 s (63% of init), one-core
2.19395 s (86%). By contrast projector assembly is 0.08779 s one-core and
rank-max 0.00442 s MPI8. This phase targets the deferred start, not preflight.

`init.states`, basis, ionic, densities and Ewald are **inside** problem
construction, despite their placement after “outside problem” in the relay.
They cannot be added to `problem`. CLI/preflight is outside it. `atomic.span`
contains radial/fill/orthonormalization; `atomic.ritz` contains Gram/eigh/rotate.
At MPI, `atomic.hq` includes owned projector-overlap preparation; the actual
complex H application is inside `split_gram`/`atomic.gram`. Thus the HQ label
does not mean the same exclusive kernel on the two paths. Likewise the
old MPI `atomic.orthonormalize` wraps block fill and gather, whereas serial
MGS excludes fill. Compare whole `atomic.span` or the new component timers,
not these differently scoped orthonormalize values.

## Code and work floor

`scf.rs::atomic_guesses` -> `atomic_start::span` -> HQ -> Ritz. One-core uses
two-pass complex MGS (`orthonormalize_pooled`, same arithmetic as serial).
MPI with adopted M2_ATOMIC_BLOCK fills a contiguous row partition, performs
guarded CholQR2, gathers the complex Q, and applies H by column partitions.
Both routes use all atomic orbitals, then kinetic-ordered completion if needed.

For the repository si64 30 Ry / half-cutoff ladder input and Si.pz-vbc:
start cutoff 7.5 Ha, m=8409 full complex G rows, n=256 atomic columns
(64 atoms x s+3p), b=132 requested bands, p=256 nonlocal channels, 54^3 FFT
grid. These shapes are **derived from the repository protocol**, not recovered
from the absent dev08 input files. New INIT_PROFILE counters must confirm
them on the actual A/B input; a different ladder fraction changes m and grid.
The initial Q payload is 34,443,264 bytes; the complex reduced matrix is
1,048,576 bytes. No H application or FFT is performed by orthonormalization.

- MGS: n(n-1)=65,280 projection calls (two passes). Each is one complex dot
  and one complex axpy: 2 m n(n-1)=1,097,879,040 complex MACs, approximately 8.78 GFLOP at 8 real FLOP/complex MAC, plus O(mn)
  norms/scales. No GEMM or Cholesky. The logical dot/update streams are
  80 m n(n-1) bytes, about 43.9 GB, before cache reuse; not DRAM traffic.
- Guarded Gamma CholQR2: two real Gram GEMMs and two real update GEMMs on
  interleaved 2m x n data, 16 m n^2 real FLOP (about 8.82 GFLOP), two
  n x n Cholesky factorizations and two triangular inverses. It changes
  arithmetic intensity rather than removing the leading operation count.
  Its existing pivot guards (1e-6 first, 0.5 second) remain unchanged;
  failure returns to the original MGS and rank rejection policy.
- One-core HQ: ATOMIC_REAL packs Q to m x n real ORC, applies the real H,
  then unpacks. There are n/2=128 pair applications, each with one full-grid
  inverse and one forward complex 3-D FFT: 256 full-grid transforms, one
  transform per band in aggregate (plus an extra pair for an odd tail).
  Two nonlocal real GEMMs cost about 4mnp=2.20 GFLOP. Grid products cost
  O(grid_points*n). FFT arithmetic depends on backend/radices; a conventional
  estimate is 5*grid_points*log2(grid_points) per complex transform, not a
  measured hardware count. Real-state setup/copies are inside atomic.hq too.
- Ritz currently uses a full complex Q^H HQ GEMM (m n^2 MACs), complex
  Q C (m n b MACs), and an n x n Hermitian eigensolve. Hermitian completion
  checks, Gram symmetrization and eigenvector selection remain unchanged.

The relay's one-core times are atomic.span 1.42627 s, orthonormalize
1.21804 s, HQ 0.64629 s, Ritz 0.12756 s, Gram 0.05376 s, rotate 0.03706 s,
eigh 0.03673 s, fill 0.20073 s, radial 0.00644 s. The MGS estimate corresponds
to about 7.21 effective GFLOP/s and 36.1 GB/s **logical** stream rate at
1.21804 s. These are algorithmic work divided by an inclusive measured
timer, not peak machine performance or proof of a memory bottleneck. No
measured GEMM throughput/cache traffic was supplied, so an absolute hardware
minimum time cannot be inferred honestly.

MPI8 already uses block algebra: roughly 2 tree sums of n^2+1 doubles at
Gamma, replicated n x n factor/inverse work, one Q allgatherv (34.4 MB
delivered per rank), and a row-partitioned Ritz rotation plus another gather.
HQ is column partitioned: 32 complex bands per rank, about 64 full-grid FFTs
per rank (complex inverse/forward), plus owned projector overlaps. The full
Q is replicated; radial tables, factors and eigensolve remain replicated.
The 0.52558 s start is not proof of optimal distribution. Per-rank Gram,
factor, update and gather timers are added to distinguish these residuals.
No inference compares a one-core timing ratio against an MPI timing ratio.

## Candidates and pre-implementation choice

1. **Local guarded block orthonormalization**, default-off ATOMIC_LOCAL_BLOCK:
   reuse the already-adopted atomic CholQR2 algebra without MPI collectives,
   only when no multi-rank communicator is active. Same column order, span,
   positive factor diagonal, guards, completion and Ritz calculation. On
   decline, run the existing MGS. MPI is unchanged. Conditional estimate for
   si64 one-core: if sequential block products sustain 20–40 GFLOP/s, their
   8.82 GFLOP costs 0.22–0.44 s, before factors/copies/fill; an expected net
   saving of roughly 0.5–0.9 s (1.2–2.1% of the supplied one-core wall) is a
   hypothesis for dev08, not a measured result. No MPI gain is predicted.
2. **Gamma ORC-native span/HQ/Ritz**: generate/orthogonalize the independent
   real coefficients and keep them through HQ and real reduced algebra.
   This replaces 2m interleaved rows by m ORC rows for block products and
   removes complex pack/unpack/full-Gram work. Work ceiling: halve the
   Gamma block Gram/update FLOP; real reduced Gram/rotation replace complex
   work. Potential one-core saving is of order 0.5–1.0 s in the reported
   span/HQ region; MPI benefit is bounded by its as-yet-unseparated local
   algebra and representation-copy time, not 50% of the whole 0.526 s.
   It requires a broader representation, parity and fallback change and is
   deferred until detailed MPI spans return.

Choose candidate 1. A52 rejected CholQR2 in **Davidson correction panels**
after a skew-system iteration regression. New evidence is the adopted,
guarded INIT-ATOMIC-BLOCK on the **fixed atomic input span**, plus the dev08
measurement identifying untouched serial MGS as 1.21804 s. This is not a
retry on Davidson corrections or an attempt to re-adopt MPI ATOMIC_BLOCK.
Still, rounding can alter the SCF trajectory: run complete SCFs and retain
strict 1e-10 Ha / standard 5e-9 Ha gates, reporting the tighter 1e-11 Ha
cross-arm check as well. No default changes, Mac performance claims, or
unconditional adoption from silicon alone.
