# Real GEMM dispatch on MKL builds: where the pure fallback runs (F11, 2026-11-16)

Review item F11 (BAANDI code performance review 2026-11-16) asked which real GEMM
call sites of an MKL build silently run the pure-Rust kernel because of a shape
condition, and to make the effective route observable. This note documents the
counters added for that and the result for si16, si32 and si64. Defaults are
unchanged; the fixes are proposals for a dev08 A/B.

Evidence labels: **O** observed on the Mac (exact counts; counts do not depend on
the machine), **D** derived by applying the dispatch rule of an MKL build to the
observed call shapes (the same code path, `rows_as_if_mkl`), **E** estimate.
No timing in this note is evidence for dev08.

## 1. What was added (count-only, no behaviour change)

`crates/dft-numerics/src/gemm_dispatch.rs` counts, per (site, route), the calls,
2 m n k flops and the largest m, n, k (output m x n, reduction k) of every real
GEMM wrapper of `RealGemmWorkspace` (`multiply_slices_into`,
`multiply_slices_wide_into`, `rotate_pair_slices`, `rotate_slices`,
`rotate_residual_slices`, `multiply_column_groups`, the direct `multiply_into` /
`rotate_pair` entries and the NARROW streaming kernels). Routes:

| route | meaning |
|---|---|
| `mkl` | per-panel dgemm (fixed 256-row / 256-term panels) |
| `mkl_wide` | dgemm only because the call used the wide entry (whole output m n > short^2) |
| `pure_build` | build without the `mkl` feature |
| `pure_switch` | MKL build with `BAANDI_MKL_DGEMM=off` |
| `pure_output_bound` | MKL build, fallback: m n > short^2 on a non-wide (or transposed) call |
| `pure_panel_bound` | MKL build, fallback: one 256-row output panel exceeds a tile |
| `pure_long` / `pure_short` | MKL build, fallback: A longer than the workspace rows / wider than `short` |
| `pure_layout` | MKL build, fallback: empty shape, leading dimension or output view mismatch |
| `pure_direct` | direct entry into the pure kernel (never tries dgemm) |
| `narrow` | streaming dot/update kernels (NARROW switches; not a GEMM) |

The route is computed by `slices_route` (and the analogous one-line conditions of
the rotation wrappers), a clause-by-clause decomposition of the original dgemm
condition; the MKL-only unit test `real::tests::slices_route_matches_original_condition`
checks over 7,200 shape/flag combinations that dgemm is chosen exactly when the
original condition held (type-checked on the Mac with `--features mkl`; it runs in
any MKL build). Pure builds additionally count the route an MKL build with dgemm
on would take for the same call (`rows_as_if_mkl`), which is how this note gets
MKL routing from a Mac run. Results go to results.json
`convergence.work.gemm_dispatch` (`rows`, `rows_as_if_mkl`, `wide_requests`,
`mkl_flop_fraction`); the source-only accounting build (`P14B_GEMM_ACCOUNTING=1`)
also records every call shape as `real.<site>.<route>` with its time and writes
them to results.json `gemm_accounting`. Overhead: five relaxed atomic adds per
wrapper call (no allocation, no clock outside the accounting build); si8 energy
and density SHA are bit-identical to the build without counters (O).

Limits of the counters: sites are the six `GemmSite` labels plus `unlabelled`;
`other` mixes several call sites (lean factor products, MGS-B panels, nonlocal
products, B5), separated below by their shapes. Complex GEMMs (atomic start,
k points, `linalg::gemm`, `GemmWorkspace`) are not in this table.

## 2. Runs

Mac release build at the commit of this note, Si.pz-vbc, 30 Ry, Gamma, the
benchmark cells `benchmarks/cases/si{16,32,64}-gamma.toml` with the last atom
displaced by 1e-3 (fractional z), threads 8, default capacity 6 per band, both
profiles (standard = lean Davidson, strict = classic residual solves; the
standard profile includes the cutoff ladder's coarse stage). The MKL build's
orthogonalization defaults differ from the pure build's (`BAANDI_MGS=blocked`,
`BAANDI_ORTHO_BLOCK=off` on MKL builds), so the runs were repeated with those two
variables set to reproduce the MKL call shapes (table below); the strict si64 run
then differs by four rotations (rounding class of the Gram-Schmidt variant).

| cell | bands | n_pw | capacity |
|---|---|---|---|
| si16 | 36 | 5915 | 216 |
| si32 | 68 | 11877 | 408 |
| si64 | 132 | 23871 | 792 |

## 3. Result (D on O shapes; MKL defaults for the orthogonalization)

| cell, profile | real GEMM total (GFLOP) | pure fallback (GFLOP) | share | where |
|---|---|---|---|---|
| si16 standard | 8.22 | 0.055 | 0.7 % | MGS-B cross-panel apply (other, n x 16, k <= 16), 18 calls |
| si32 standard | 45.53 | 0.584 | 1.3 % | MGS-B cross-panel apply (other, n x 16, k <= 48), 48 calls |
| si64 standard | 453.48 | 0 | 0 % | none |
| si16 strict | 21.66 | 2.204 | 10.2 % | projection apply X -= Q S (79 of 156 calls, 1.898 GF) + MGS-B cross-panel (106 calls, 0.305 GF) |
| si32 strict | 162.24 | 19.228 | 11.9 % | projection apply (73 of 116 calls, 15.260 GF) + MGS-B cross-panel (338 calls, 3.968 GF) |
| si64 strict | 1725.52 | 178.792 | 10.4 % | projection apply (105 of 249 calls, 178.792 GF) |

With the pure-build orthogonalization (`ORTHO_BLOCK=on`, which already takes the
wide route) only the projection apply falls back: 8.7 / 9.3 / 10.3 % of the
strict totals (O shapes, D route).

Every other site runs dgemm: rotations (V C, H V C; Ritz return, contraction),
the overlap and reduced-matrix borders (A^T B, whole output <= capacity^2 by
construction), the lean factor and reduced-border products, B5, the nonlocal
P^T X (A^T B) and back product P W (already wide: `mkl_wide`, 55 to 110 calls,
1.0 to 132 GF).

### Why: the whole-output bound

`multiply_slices_into` takes dgemm only if the whole output fits the square
workspace, m n <= short^2 with short = capacity = 6 x bands. For the long apply
X -= Q S the output is n_pw x c (c = block width), so it falls back when
c > 36 bands^2 / n_pw: c >= 8 at si16, 15 at si32, 27 at si64 (D; the P49A row
has the same thresholds). A full-width block (c = bands) falls back whenever
n_pw > 36 x bands. For Si at 30 Ry n_pw / bands is about 165 to 181, and for
plane-wave calculations in general the number of plane waves per band is far
above 36, so **full-width projections of the classic path fall back in every
realistic cell**, more so at higher cutoff (more plane waves per band) and with
smaller capacity (e.g. the MPI rule of 3 columns per band, although the
distributed solvers take the wide route by default). The MGS-B cross-panel apply
(n_pw x 16) falls back when 16 n_pw > 36 bands^2, i.e. when bands < (2/3)
sqrt(n_pw): the small cells (si16: 36 < 51, si32: 68 < 73; not si64: 132 > 103).
The bound is not needed for A B products: `dgemm_panels` writes each 256-row
output panel directly and uses one tile per panel, which is what the wide entry
already relies on (`multiply_slices_wide_into`). It is needed for A^T B, whose
partial sums are tiles of the whole output.

### Strict and certificate

The strict profile and the fixed-point verification (certificate solves, ON by
default in the single path) use the classic residual Davidson, whose basis
projection is the site above; so the fallback also sits inside every timed
single-path run through the certificate. The lean standard solves themselves
append unprojected corrections and have no fallback at si64.

### Atomic start

The atomic start (`dft-engine/src/atomic_start.rs`) uses the complex
`linalg::gemm` (MKL builds: one sequential zgemm per product; pure: faer
`Par::Seq`): Q^H (H Q) and Q Y, n_pw x 4 n_atoms columns. It is outside the real
dispatch (no shape fallback) but single-threaded; the threaded fixed-column-group
form is already a backlog item. With `ATOMIC_REAL=on` (default) H Q itself is
applied by the real Hamiltonian, whose two projector products are `mkl` / `mkl_wide`.

## 4. Cost (E)

Extra time of a fallback = F_fallback x (1 / r_pure - 1 / r_dgemm) per core. With
r_dgemm about 40 GFLOP/s and the MKL build's pure kernel (`PURE_GEMM=nano`) about
8 to 12 GFLOP/s per core (assumed; to be measured), si64 strict at T=1 carries
178.8 GF x (1/10 - 1/40) = about 13 s of extra single-core work (E), si32 strict
about 1.1 s, si16 strict about 0.16 s. Relayed dev08 data that points the same
way: arm C (`PROJECT_WIDE=on`) of the lean A/B was -1.0 to -4.3 % in every cell
(B-LEAN-DAVIDSON-RECOUNT); that arm ran the classic standard path, whose
projections are this site.

## 5. Proposed fixes (not applied; defaults unchanged)

1. **Measure the existing switch on the strict and single-path runs.**
   `BAANDI_PROJECT_WIDE=on` takes both fallback sites (projection apply and the
   MGS-B cross-panel apply) through the per-panel dgemm. dev08 A/B: si16, si32,
   si64, strict and S8 (single path with verification), T=1 and T=8, seeds
   601-610 (si64 601-603), med_by_seed; with this commit's counters the MKL
   `rows` must show zero `pure_output_bound` calls in the on arm.
2. **Replace the switch by a shape rule** (if 1 is positive): in
   `multiply_slices_impl` apply the whole-output bound only to A^T B products;
   A B products keep only the per-panel bound (the rule the wide entry already
   uses). This makes every A B call site (projection, MGS-B, B5 coefficient
   product, future sites) use dgemm without a per-site switch; pure builds stay
   bit-identical (both routes forward to the same kernel); MKL builds are
   rounding class (same fixed panels, so thread-count independent). Gate: the
   P49A unit tests (wide vs square within 2 k eps sum |a||b|, bitwise across
   T=1,2,4,8), `bench_scf.py` energy within 1e-12 Ha, strict 1e-10 / standard
   5e-9 Ha classes.
3. **Keep the counters in the dev08 runs**: on an MKL build `rows` are the real
   routes, so the A/B report states the dgemm flop share (`mkl_flop_fraction`)
   per cell with no inference.
4. Not proposed: raising `short` (more tile memory, and the A^T B bound is
   structural), or a pure-kernel speedup for these shapes (the vendor route
   exists).

## 6. Reproduction

```
cargo build --release -p dft-cli
python3 -I scripts/gemm_dispatch_survey.py target/release/dft-cli benchmarks/cases <Si.pz-vbc.UPF> <outdir> standard|strict
BAANDI_MGS=blocked BAANDI_ORTHO_BLOCK=off python3 -I scripts/gemm_dispatch_survey.py ...   # MKL orthogonalization defaults
```
(the script writes the input of section 2 and prints `gemm_dispatch.rows` and
`rows_as_if_mkl`). Unit tests:
`cargo test -p dft-numerics --test gemm_dispatch`;
`cargo test -p dft-numerics --features mkl real::tests::slices_route_matches_original_condition`
(MKL machines).
