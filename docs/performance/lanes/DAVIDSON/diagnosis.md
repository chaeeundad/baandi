# Phase 2: Davidson upper-border experiment

## Evidence and limits

Source of measured seconds: the user's dev08 relay, source 2733278, physical
CPUs 0–7 on socket 0, instrumented runs. All 15 OFF/ON integrity pairs had
zero energy difference, identical scientific history/verification/density and
zero dropped records. No Mac timing is used. Working base is 45670dd633af6573b68578fd2e3d5700b982c889.

| Model/cell | External wall s | Solve envelope s | Solve/wall | Rank max/min |
|---|---:|---:|---:|---:|
| MPI8 Si64 | 9.7715 | 6.4323 | 65.8% | 1.004 |
| MPI8 Si32 | 1.9254 | 1.1224 | 58.3% | 1.005 |
| MPI8 Si16 | 0.6625 | 0.2996 | 45.2% | 1.002 |
| One core Si64 | 45.5618 | 27.2416 | 59.8% | n/a |
| One core Si32 | 8.1225 | 4.5788 | 56.4% | n/a |

Different models are separate experiments; these rows are not speedup ratios.
The MPI8 Si64 envelope agrees with the existing solver_total 6.4384 s to
0.1%. The serial solver_total is zero/unusable; use the lane envelope.

| Inclusive segment, Si64 | MPI8 s | One-core s |
|---|---:|---:|
| reduced_gram | 0.9263 | 2.1793 |
| rotate_residual | 0.8877 | 5.8842 |
| overlap_border | not separately relayed | 2.1248 |
| reduced_eigen_step | 0.7076 | 1.2804 |
| reduced_eigen | 0.7076 | 1.1714 |
| tree_sum | 0.3666 | n/a |

Seconds above are derived from the rounded percentages in the relay, not
additional measurements. Do not add inclusive rows. The relayed exclusive
split percentages (45.9 apply, 15.2 Gram, 15.0 eigen, 14.6 rotate, 9.5 kept,
3.9 residual) sum to 104.1%; their aggregation/denominator needs clarification.
They cannot safely be turned into an additive budget. Treat apply as owned by
FFT/NONLOCAL, including operator communication; this experiment changes no H
application. Rank imbalance <=0.5% does not justify a load-balancing change.

**Missing input:** the relay contains no actual per-call n/m/k records,
counts, payloads or hardware throughput/bandwidth measurements. A local path
to results.json / tables.json was requested. Numerical accumulated FLOP/byte
floors and measured-time/floor gaps cannot be inferred from percentages.
Historical strict-path sizes in davidson-algebra-gap-analysis.md are not the
current standard trajectory; substituting those would fabricate a diagnosis.
The following exact shape formulas specify what to compute when data arrives.

## Arithmetic, data movement and communication floors

Let n be local PW rows (global on serial), m basis order, k output width,
p=m-k retained prefix. All counts refer to FP64, FLOP = multiply plus add.
For each rank sum work over calls; use the critical rank, not the sum over
ranks, for a wall floor. For rate ceilings P FLOP/s and B byte/s, a simple
optimistic roofline floor is max(F/P,D/B). The gap is measured inclusive
seconds divided by that floor (or subtract the floor for the headroom in s).
Rates must be for the relevant memory level and simultaneous-core placement;
minimum array accesses are not necessarily DRAM traffic. Cache residency can
invalidate a DRAM interpretation. No rate ceiling was supplied in the relay.

* Rotate: two NN GEMMs (n x m)(m x k) -> n x k, F=4nmk.
  Minimum distinct-array bytes = 8(2nm+mk+2nk), before residual reads/writes.
  This floor excludes residual formation/norm work inside the serial timer.
  Serial: gap 5.8842 / floor_rotate_serial. MPI: 0.8877 /
  floor_rotate_critical_rank. A rank's n is NOT global n/8 by assumption.
* Gram: TN product (m x n)(n x k) -> m x k, F=2nmk.
  S=V^T V border has distinct-input union nm, plus mk output doubles.
  H=V^T HV border needs nm+nk input and mk output doubles. For a combined
  distributed S/H record, trace k=ks+kh: F=2nm(ks+kh); conservative distinct
  array floor is 8(nm+m(ks+kh)); accurate H traffic needs kh separately.
  A fresh S block is SYRK-shaped, n x m -> symmetric m x m; unique arithmetic
  n*m*(m+1), versus 2*n*m*m for a full GEMM. The H block is symmetric only
  as a projected self-adjoint operator, not a same-operand SYRK; preserve
  the existing upper-entry definition rather than assuming bitwise symmetry.
  Serial reduced/overlap gaps are 2.1793/floor_H and 2.1248/floor_S.
  MPI combined gap is 0.9263/floor_combined, including reductions.
* EVD: replicated full order-m dsyevd (unless overridden), not distributed
  by PW rows. Record call count, order distribution and sum(m^3) per rank;
  **do not divide by 8**. Dense Householder tridiagonalization alone has the
  conventional leading work estimate (4/3)m^3; full eigenvectors add work.
  This is an algorithmic leading-term proxy, not an exact FLOP count or a
  universal eigenproblem lower bound. At least the input and output m x m
  arrays must be accessed; the BLAS2 passes require more traffic. Compute
  and bandwidth lower estimates therefore need an explicit dsyevd model or
  a same-size dev08 calibration. Gaps 1.1714/floor_EVD_serial and
  0.7076/floor_EVD_MPI remain unquantified. reduced_eigen_step contains the
  EVD, B5 attempts/cache and surrounding work; never add the pair.
* MPI border tree sum: current code sends m*(ks+kh+1) doubles, including a
  full m-element residual-norm tail even when only part is populated. Each
  border invocation has one reduction; S and H normally occur on separate
  iterations of the append/apply sequence. Payload is 8*m*(ks+kh+1) bytes
  per rank; network bytes depend on tree implementation. Actual calls and
  payload sum must come from tree_sum records with site reduced_gram or
  lean_border. Exclude apply site (another lane). A latency/bandwidth model
  requires measured collective latency L and effective Bnet, e.g. a
  recursive-doubling bound involving log2(P) stages, not payload/B alone.
  Tree_sum 0.3666 s includes other sites and must not be added to Gram.

The parser already retains the needed shape/count inputs. Its minimum-array
model is deliberately conservative, not a memory-controller measurement.

## Candidates and pre-registered hypotheses

1. **Chosen: upper-border GEMM tiles (default OFF).** Code computes full
   m x k rectangles then discards rows i>j. Partition new columns into
   fixed 64-column tiles, each ending at e; compute only e x width output.
   Actual F=2n sum(e*width); retained mathematical minimum is
   2n sum(j+1). No inner-axis split, new allocation, MPI collective, changed
   tolerance, restart rule or H column. Serial overlap and reduced matrix,
   distributed lean S/H border use this helper. Classic distributed,
   return Gram and optional fused-border path are unchanged.
   Example shapes (illustrations, NOT trace observations): fresh m=k=132
   saves 26.4% of rectangle work; m=396,k=132 saves 8.8%; width<=64 saves
   zero. Extra calls may reread/repack V. This is deliberately a GEMM-only
   experiment with existing adapters, not a new SYRK FFI.
   Conditional dev08 hypothesis if affected time falls by 10–25%: MPI8
   Si64 0.093–0.232 s (0.95–2.37% wall) from the 0.9263 s envelope; serial
   Si64 0.430–1.076 s (0.94–2.36% wall) from H+S's 4.3041 s. These are
   Amdahl scenarios, not measured or trace-weighted predictions; extra
   dispatch/packing could erase the saving. No positive numeric prediction
   for Si32/Si16 without their Gram shares. Rotation/EVD savings: zero.
2. **Deferred: triangular reduced-coordinate products.** T is explicitly
   lower triangular; T^T Y and T P currently multiply stored zeros. Fixed
   row/column panels could omit structurally zero regions while retaining
   library GEMMs. No basis/stopping change. Each such product can approach
   half its arithmetic (less with blocking), but it affects coefficients
   and transformed reduced matrix, not full dsyevd. Serial/MPI expected
   wall gain is 10–30% of those separately measured product times, which
   were not relayed. Optimistic upper limits if ALL non-EVD eigen_step time
   were removable: serial 0.1090 s (0.24% wall); MPI rounded step/EVD shares
   give no resolvable budget. These are not precise bounds on coefficients
   outside that envelope. Low priority until product trace data arrives.

Why not the largest serial rotation row? Selected-root rotation and fused
Ritz/residual routes already exist; the new relay does not identify an unused
width subset or redundant rotation. Removing a residual check without that
information would change the trajectory. No capacity/tolerance tuning here.

## Prior decisions respected

Read ledger REL-FINDING and REL-SLICES-CLOSED: sliced MRRR 0.43–0.53 x EVD
failed orthogonality (up to 6.9e-13); inverse-iteration slices preserved
orthogonality but cost 1.42–1.64 x. No new EVD evidence reopens it. Advance
increased H by 37–38% on two Si16 seeds; restart-keep increased rotations
and contractions with no consistent SCF win. Async alltoall was only 2–3%
of solve and rejected; PR-CG/direct-minimization H work lost on larger cells.
The backlog's fresh-basis symmetric-border opportunity is distinct from
those rejected algorithm changes. No QE source was read.
