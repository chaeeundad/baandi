# G1-NLCC: real MPI correctness request

Status: prepared for the orchestrator / baandi-dev08; no remote run performed
by the implementer. This is a full-convergence correctness comparison, not a
speed adoption or QE timing claim. G1 changes no numerical default. Record the
reviewed commit, source/binary hashes, compiler/backend, all commands and exit
codes, and raw stdout/stderr. Do not install anything or fetch new datasets.

## Synthetic MPI8 comparison

Use one 600 s fast job per exact test invocation. Bind eight ranks to eight
physical cores on socket 0 with one thread per rank; report the actual binding
and MPI version. Build the integration binary with the existing toolchain:

```sh
cargo test --release -p dft-engine --features dft-mpi/mpi --test refscf --no-run --message-format=json > g1-build.json
```

Resolve `executable` from the `compiler-artifact` row whose target name is
`refscf` and profile has `test=true`, then run that exact binary (not cargo)
under the site's established MPI launcher:

```sh
OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 \
mpirun -np 8 --bind-to core --map-by ppr:8:socket \
  /absolute/refscf-binary --exact nlcc_real_mpi_legacy_comparison \
  --ignored --nocapture --test-threads=1
```

The test constructs the in-memory Gaussian-core NC fixtures, including the same local
cubic LDA-PZ and nonlocal skew LDA-PW variants as the simulated matrix, runs legacy with the registered MPI8
communicator and reference with that same communicator, and checks delivered
energy/density/occupied roots and independent verification. Standard and tight
estimator thresholds are both tested; caps are `5e-9 Ha` and `1e-10 Ha`.
No UPF file is written. Run the exact invocation twice for repeatability, plus
the same binary at MPI1 and MPI3 as rank-count controls. Retain `G1_NLCC_MPI`
rows. No wall ratio is requested; any later timing comparison must keep the
same parallel model and rank count on both sides.

Also run the simulated full-trace acceptance matrix outside mpirun:

```sh
cargo test --release -p dft-engine --test refscf nlcc_delivered -- --ignored --nocapture --test-threads=1
```

It covers LDA-PZ/PW, local/nonlocal, cubic/skew, twice per 1/3/8 layout and
prints `G1_NLCC` / `G1_NLCC_BRANCH` evidence. A rank-count difference is accepted
only by the existing rank-rule predicate and its stored margins, never by
loosening a density or eigenvalue tolerance. The real-MPI harness checks the
accuracy class and verification but does not itself compare traces across
separate mpirun invocations; retain those records if extending that gate.

## Optional real LDA NLCC cell

Inventory only the already available dev08 files. If a genuinely LDA,
norm-conserving, scalar/nonrelativistic NLCC dataset with atomic charge and
l<=3 orbitals exists, record its source/license, exact header and SHA-256 and
select one documented closed-shell Gamma cell. Run legacy and reference at
MPI8 with identical cell/positions, cutoff, FFT dimensions, file, band count,
start, accuracy threshold and binding, and independent verification enabled.
Use an explicit LDA-PZ or modified LDA-PW match; **no functional override**.
Record all energy terms and NLCC ledger, charge, density maximum difference,
occupied eigenvalues, iteration/verification counts and convergence outcomes.
Compare the reference standard result to a tight certified legacy result
within `5e-9 Ha` (and tight-to-tight within `1e-10 Ha`). Preserve inputs and
raw output/hashes on dev08; never commit the UPF.

If no such file exists, report this arm as unavailable. The available
PseudoDojo NC-SR PBE files must **not** be run under LDA and described as a
physical validation. Their physical reference-driver comparison waits for G2
(PBE). No new QE reference is claimed or requested by the synthetic tests.
