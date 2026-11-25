Scope extension approved for G3-KPOINTS (orchestrator decision): you may also edit
  - crates/dft-numerics/src/linalg.rs   # complex Cholesky / triangular solve / Hermitian EVD entries with the same pivot and rank rules as the real ones
  - crates/dft-numerics/tests/   # contract tests for those entries only
  - crates/dft-mkl/src/adapter.rs and crates/dft-mkl/tests/   # complex (zgemm/zpotrf/ztrsm/zheevd) entries only, if the MKL route needs them; unsafe only in the adapter, preconditions documented and tested
Keep the real Gamma path bit-identical. Continue with the remaining stages of the ticket: complex reference Davidson, per-k state and driver, CLI records, SCF agreement tests (Si2 and skew cell, shifted and unshifted 2x2x2 meshes and a k list, ranks 1/3/8, standard class and 1e-10 Ha tight vs legacy, eigenvalues per k), Gamma bit identity. Write the report in English (project rule).
