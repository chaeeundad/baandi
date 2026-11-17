# GitHub issue #19: Reuse distributed error-agreement scratch and measure the new hardening overhead

Reviewed at `1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4` (2026-11-17). This is a static review of the implementation and committed experiment records; no new BAANDI/QE timing run was performed in this review. Related tracking issue: #13.

**Suggested priority:** P2 — measurable implementation cost.

### Concrete cost added by the hardening

The new failure agreement and returned-vector checks are valuable. Keep them. Their implementation currently adds avoidable allocation/copy work that existing “same call counts” tests do not quantify.

[distributed_davidson.rs:1600](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-solvers/src/distributed_davidson.rs#L1600) allocates two zero-filled vectors on **every solve**:

```text
cap = max(max_subspace, bands * capacity_per_band).min(n)
payload = 2 * cap * cap + 2 * cap + 1
scratch = (vec![0.0; payload], vec![0.0; payload])
```

Calculated allocation sizes, not measured timings:

| bands | capacity | two scratch buffers per rank per solve |
|---:|---:|---:|
| 68 | 204 | 1,338,256 bytes (1.28 MiB) |
| 132 | 396 | 5,030,800 bytes (4.80 MiB) |
| 132 | 792 | 20,097,808 bytes (19.17 MiB) |

[ApplyAgreement::tree_sum_f64:469](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-solvers/src/distributed_davidson.rs#L469) copies the complete reduction input into scratch and the result out again to append one failure-status value. The buffers disappear after the solve despite a persistent solver workspace.

Also, [workspace/lean.rs:607](https://github.com/chaeeundad/baandi/blob/1c52998f9ad11c22e4fbdbaf17cab5c2bc169ab4/crates/dft-solvers/src/distributed_davidson/workspace/lean.rs#L607) now verifies the actual returned block with a local Gram GEMM and packs the BxB Gram with the residual reduction. The payload grows from B to B(B+1) doubles: 68→4,692 or 132→17,556. Collective **count** stays unchanged while bytes and local arithmetic increase.

These observations do not establish a wall-time regression; that must be measured on the hardened current head.

### Proposed work

1. Give agreement scratch a persistent owner with the same lifetime as the distributed solver workspace. Grow on capacity change; reuse across SCF maps. Avoid an unkeyed global cache or stale failed-status state.
2. Where reduction buffers are already owned, reserve a status slot and reduce directly into the caller's buffer, eliminating full input/output copies. Preserve deterministic reduction order and failure propagation.
3. Instrument allocation bytes/count, copied bytes, Gram flops/time and reduction bytes separately. Time profiling-off production runs; do not infer performance from counter equality.
4. Keep the actual returned-vector `X^T X` check and collective error agreement. The pre-apply output clear/nonfinite scan also has a cost; change either only after establishing a complete-write operator contract and adversarial tests.

### Acceptance

Adversarial rank-local operator failures, empty ranks, nonfinite output, unstable overlap factors and returned-vector defects remain detected without a hang. Healthy-path energies, work and collective ordering retain the required identity. After warmup at fixed shape, agreement scratch performs no new heap allocations. A/B uses the current hardened head, not the earlier unhardened lean-solver benchmark.

### Validation and adoption

Use the repository's current two-lane fast-loop protocol. Start with one short case and an identity/correctness smoke; use full-SCF Stage 1 for any trajectory change. Final acceptance uses the designated MPI benchmark host, the current default stack, fixed-point verification included, the existing accuracy classes and same-model/same-rank QE comparisons. Measurable small improvements remain eligible under the adoption rule; do not impose a new minimum percentage. Report failed or negative experiments as well as wins.
