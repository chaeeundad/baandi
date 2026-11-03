# ADR 0001: isolate native numerical backends

- Date: 2026-11-02
- Task: P0-BOOTSTRAP
- Status: proposed for reviewer acceptance; implemented in P0

## Context

Numerical conventions and reproducible validation take priority over speed.
Network access is unavailable; Apple Silicon provides a system BLAS/LAPACK.

## Decision

Reserve dft-numerics as the sole backend boundary. Deny unsafe by default and
allow it only in a private, documented FFI adapter. Use the SDK's legacy LP64
CBLAS ddot and CLAPACK dsyev ABI for fixed-size P0 smoke tests. Do not expose
this probe as the P1 API. Keep all nine requested crate boundaries and zero
external crate dependencies. Pin Rust 1.99.1 and document the offline alias.

## Alternatives

A modern Accelerate LAPACK interface would avoid deprecated entry points but
requires separate symbol/ABI validation. A third-party BLAS wrapper requires
a dependency decision. Pure Rust tests alone would not verify native linkage.

## Consequences

Real-valued ABI evidence is available now. Legacy entry points are deprecated
in the installed SDK; P1 must decide whether to migrate. Complex layout,
ILP64, FFT and performance remain unvalidated. Tests deliberately reject
non-aarch64-macOS targets so absence of native coverage cannot look like success.

## Validation

Known exact dot product and analytic symmetric eigenvalues, plus eigenvector
residual and orthonormality checks; absolute tolerance 1e-12. See environment.md.

## Future ADR template

Use: title/number, date, task, status, context, decision, alternatives,
consequences, validation, and superseding ADR if applicable. Record measured
evidence separately from assumptions and outstanding questions.

## P7P6 notes

SCF orchestration may distribute disjoint mutable k states using
`std::thread::scope`. Hamiltonians remain intentionally non-Sync; ownership
transfer needs only their existing Send property. No numerical backend bounds,
unsafe code, dependencies or nested backend parallelism are introduced.
