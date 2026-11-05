# Study note 4: adopted incremental Davidson Gram construction

Date: 2026-11-05. Only rust-dft and existing prose notes were read; no QE source
was read. The architect explicitly adopted incremental Gram construction after
reviewing the initial experiment. It is now the default in both pure Rust and
MKL, with no compile-time policy switch. MKL execution and timing await dev08.

## Decision and acceptance rule

The architect superseded this ticket's original all-system identical-SCF-count
rule for this change: realistic fixtures must retain SCF count and total energy
within 1e-12 Ha; converged toy fixtures must retain the existing suite energy
tolerance (1e-12 Ha), but reported iteration-count changes are allowed. No
previously converged case may become nonconvergent. Final certificates and their
tolerances remain unchanged. Existing nonconvergence is explicitly reported,
not reclassified as success. The historical fixed-point multiplicity note is
narrative-only evidence; this adoption relies on the explicit architect decision
and the newly stored numerical comparisons, not an independently proved cause
of each rounding-sensitive trajectory change.

[Current acceptance tables](../../../tests/reference/perf/p14b/adopted/acceptance.md)
list **all 19 cases**, every iteration count change, and T=1,2,4,8 checks.
[Initial accounting tables](../../../tests/reference/perf/p14b/gemm-tables.md)
retain the original ordinary/certificate phase attribution and initial trial.
Raw inputs, executable/UPF hashes, counters and outputs accompany both tables.
The original patch is retained only as historical experimental evidence; applying
it to the adopted source is neither necessary nor supported.

The four toy changes remain syn-fft-primes 14→16, syn-p3 8→10, syn-s1 13→14,
syn-s2 17→16. All 16 formerly converged suite cases still converge and meet the
suite energy tolerance. Three existing failures remain: si-elongated, syn-s3,
syn-z14. Therefore this is **not a claim that all 19 converge**.

The four requested fixtures retain SCF counts 15/10/6/9 and energy differences
at most 2.842171e-13 Ha. Added 20-atom Gamma cubic coverage includes both an
8-Bohr dense cube (100-iteration nonconvergence before and after) and a 12-Bohr
cube (13 SCF iterations before and after; ΔE=9.379164e-13 Ha). Both are stored;
neither changes tolerances. The 12-Bohr cube's H-vector work increases
1996→2107, so GEMM savings do not imply every kind of work always decreases.

## Algorithm and backend implementation

If an unchanged basis grows from p to q columns at N plane waves, retain both
triangles of the p×p reduced Hamiltonian and compute two rectangles: new rows
against all columns, then old rows against new columns. Work becomes
2N(q²-p²) instead of 2Nq² in the historical 2mnk convention. No triangle is
inferred by conjugation, preserving the independent Hermiticity check. The
solver invalidates the cache on restart or untrusted/replaced H products.

Pure Rust keeps the fixed 256-plane-wave panels and ordered assembly, independent
of thread count. Changed rectangle shapes change backend rounding; baseline
and adopted density files are not bit-identical. Within the adopted pure-Rust
backend T=1,2,4,8 energies, eigenvalues, densities and H-work counts are identical
on all 19 suite cases and the additional fixtures.

MKL borrows contiguous column slices of the large inputs and calls
`dft_mkl::gemm` once per nonempty rectangle, sequentially at every thread budget.
A workspace buffer allocated once holds the contiguous small output required by
the existing adapter. Outputs are scattered into the reduced matrix. No large
operand packing/copy is introduced. The small copies are explicit: 16p² bytes
for the preserved prefix plus 16(q²-p²) bytes for MKL output scatter per
extension, totaling 16q² bytes. This is a small-copy increase over the former
MKL direct-output full GEMM, **not** fulfillment of the original no-copy-increase
aspiration. Counts are recorded through the existing worker-analysis hook;
wall-time impact and MKL numerical/determinism acceptance remain unverified.
The pure-Rust path copies only the preserved 16p²-byte block and sums tile outputs.

Reduced result allocation follows the pre-existing allocating Davidson API;
tile scratch and MKL output scratch are reserved at workspace creation. This is
not a claim that the whole legacy Davidson loop is allocation-free.

Independent matrix tests compare every entry with scalar dot products, preserve
old entries exactly, cover non-Hermitian products so a conjugate-fill shortcut
cannot pass, cover partial panels and T=1,2,4,8, and check restart/full-product
identity and invalid dimensions. MKL compiles the same tests, but they must be
executed remotely with the library linked.

## Accounting and next reduction

Source-only `P14B_GEMM_ACCOUNTING=1` instrumentation reserves a bounded collector
before the first numerical loop. It records solve ID, site, shape, calls, 2mnk,
conservative minimum and logical-call wall time; default builds omit clocks and
updates. The CLI drains it through the P14A reporting boundary after SCF, never
inside numerical loops. It covers explicit Davidson workspace products, not
GEMMs internal to EVD or Hamiltonian/nonlocal kernels. Times include scheduling,
assembly and validation; they are not isolated microkernel or packing CPU times.
Concurrent completion order is not used to infer certificate identity.

Adopted work: si8 21,636,961,594→9,489,226,838 (56.14% reduction); k444s
3,300,404,692→2,258,007,134 (31.58%). The synthetic skewed and elongated cases
also decrease. Exact current site totals and second-pass ceilings are generated
in the acceptance table. Shared-host timings overlapped other work and are
observations only, not accepted speedup measurements.

The next largest site is all-root Ritz/residual formation. Current roots are
unlocked: a root below tolerance in an earlier step may fail after the next
subspace rotation. All current residuals are needed to determine convergence
and select corrections. Consequently the **guaranteed safe removable bound is
zero** under the current stopping algorithm. The full site cost is only an
upper ceiling for a future locking algorithm, not a realizable saving.

Correction columns already include only unconverged roots. Buffer thresholds
and certificates cannot be dropped. There is no separate generalized overlap
matrix. Projection against only the newest block is invalid for arbitrary
corrections. The two outer projection passes are separately counted; internal
MGS is scalar, not GEMM. An adaptive second-pass criterion still needs the
second overlap after normalization, and must handle near-dependence. Skipping
its subtraction/MGS changes arithmetic and lacks a validated universal gate.
No second reduction is adopted. Restart H contractions are already reused in
the measured runs. No reusable backend packing API or packing speedup is claimed.

## Reproduction and gates

```sh
PATH=/opt/homebrew/opt/rustup/bin:$PATH
cargo build --release -p dft-cli --offline
python3 -I tests/p14a_work_accounting.py --cases si8,k444s,syn20-skew,syn20-cubic12,syn4-elongated --threads 1 --output /tmp/p14b-t1.json
# Repeat with 2, 4, 8 and distinct output paths.
P14A_WORK_ACCOUNTING=1 P14B_GEMM_ACCOUNTING=1 cargo build --release -p dft-cli --offline
python3 -I tests/p14a_work_accounting.py --output /tmp/p14b-accounting.json
P14B_GEMM_ACCOUNTING=1 cargo test --release -p dft-numerics accounting_shapes -- --ignored
python3 -I tests/p14b_adoption.py
```

The additional cubic baseline was produced by temporarily ignoring `previous`
in `reduced_hamiltonian`, then building `/tmp/p14b-full-gram-cli` and restoring
the default immediately. This invokes the same original full product, without
changing solver policy. Its binary hash and exact inputs are in the cubic
baseline JSON. Current source hashes are stored with the adoption report.

The unmodified robustness harness still exits 5 for the four now-authorized toy
iteration changes. Cross-thread harness comparisons exit 5 because thread count
is part of input identity. The separate `p14b_adoption.py` validator enforces the
architect's exact numerical rule and scientific cross-thread identity without
rewriting baselines or concealing these harness exits. Three pre-existing
nonconverged cases and the dense cubic challenge are reported separately.

Full build/test/fmt/clippy, MKL all-target type check, SCF historical-baseline
comparison, diagnostic tests and README regeneration/tests are recorded in the
adoption report. MKLROOT is unavailable locally. No dev08 run, new QE reference,
or accepted wall-time speedup is claimed. The original report and raw snapshots
are historical; the adoption report supersedes their blocked status. Changes
remain uncommitted. Dashboard uploads and the canonical delegate registry remain
outside this ticket's owned files; no recurring report timer is active.
