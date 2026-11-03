# Building and validating the optional oneMKL backend

Default: `cargo build --release -p dft-cli`. No MKL installation, headers or
runtime are needed, including for `scripts/check.sh` on Apple Silicon.

Declaration check on the Mac (does not link or run MKL):

```sh
cargo check -p dft-numerics --features mkl
cargo check --workspace --all-targets --features dft-numerics/mkl
```

Execution is supported on Linux x86-64 with a user-installed LP64 oneMKL:

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export LD_LIBRARY_PATH="$MKLROOT/lib/intel64:$MKLROOT/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
cargo build --release -p dft-cli -p dft-numerics --features dft-numerics/mkl
```

Select both CLI and numerics packages for an MKL build: the CLI has no direct
MKL feature and Cargo does not accept a transitive feature on the CLI alone.

The build script selects `lib/intel64` or `lib` under MKLROOT, checks for
`libmkl_core.so`, and links sequential MKL (no OpenMP layer, no `mkl_rt` dispatcher).
Unset MKLROOT allows cargo check; an executable using MKL cannot link. An invalid
MKLROOT or unsupported execution target fails explicitly. Changing MKLROOT reruns
the build script. Do not set MKLROOT on the Mac or interpret cargo check as an ABI,
link or numerical test. Do not enable ILP64. The header probe below detects an
incompatible library interface before tests.

## Exact dev08 protocol for Claude

Run from the reviewed checkout on an idle node. Use the existing remote UPF;
never transfer the local Mac UPF. Set its absolute path and the available physical
CPU set (default 6-13), plus a single core for kernel timings (default 6). This performs no install or remote configuration change:

```sh
export MKLROOT=/opt/intel/oneapi/mkl/2026.1
export P10B_UPF=/absolute/existing/dev08/path/Si.pz-vbc.UPF
export P10B_CPUS=6-13
export P10B_CORE=6
export P10B_OUT="$PWD/target/p10b-evidence"
export MKL_BENCH_MAX_MIB=2048
export MKL_BENCH_REPS=10
bash crates/dft-mkl/tests/dev08.sh
```

The UPF path is the only machine-specific missing input. `MKL_BENCH_MAX_MIB=2048`
allows all requested cube edges through 192; reduce it on a memory-constrained
node. The harness prints skipped sizes instead of allocating beyond its budget.
The estimate includes concurrent outer-worker buffers, but MKL's private workspace
is opaque: leave RAM headroom. Threads are created outside timed kernel loops.
All commands use fail-fast behavior, retain full logs, and never relax thresholds.
Run under a shell with enough disk space for separate Rust/MKL release targets.
The script may take substantially longer than the normal debug gate.

Success requires all of the following, not just the last printed line:

- ABI probe: `PASS MKL LP64 ABI constants`, long/int/complex sizes 8/4/16 bytes.
- Default check: four PASS gates. MKL release tests and ignored release tests:
  zero failures; FFI tests are enabled only by feature `mkl` and require linking.
  Existing synthetic Hamiltonian, solvers and SCF tests run against the alias;
  implementation-specific Rust lane-count tests still test Rust explicitly.
- `fft-{1,2,3}.log`: `PASS FFT accuracy and thread determinism`. Every transformed
  coefficient is compared with rustfft using complex relative error <=1e-13
  (no absolute floor; exactly matching zero has error zero); backward and forward
  are both checked in FFI tests. Round-trip maximum absolute error and relative
  Parseval defect <=1e-13. Small noncubic direct DFT checks cover every coefficient;
  larger grids sample Kahan-summed coefficients, absolute error <=1e-13. The
  Kahan oracle still uses f64 trigonometry and is not arbitrary precision.
- MKL results bit-identical across 1/2/4 FFT workers and 1/2/4/8 SCF budgets.
  Engine tests also check eigenvectors, eigenvalues, convergence records and
  residuals, and Hamiltonian tests include varied band counts and projectors.
- SCF cross-backend delta <=1e-12 Ha in every case; fresh MKL outputs pass unchanged
  real QE Gamma/k-point comparators (the stored QE provenance is preserved).
  Archive total energies, charge, eigensolver residuals/orthogonality and complete
  convergence records in `qe-*/results.json`, not just timing summaries.
- Report median and range of the three FFT/Hamiltonian runs and SCF wall-time
  samples, and hardware/compiler/MKL metadata. Compare same-host Rust and MKL,
  not Mac timings. Kernel GEMM/EVD timings are diagnostic; complete SCF wins govern
  any recommendation. `scf-rust.json` is the contemporaneous backend baseline;
  also compare the existing baseline below without changing it.

```sh
python3 -I scripts/bench_scf.py --binary "$PWD/target/p10b-mkl/release/dft-cli" \
  --upf "$P10B_UPF" --cases gamma,k222 --threads 1 --repeats 3 \
  --baseline docs/performance/baseline.json
```

Near-zero Fourier coefficients can expose failure of the strict per-element
relative gate even when normwise error is small. Preserve such failures and
investigate; do not substitute a looser normwise assertion. A skipped large grid
is missing evidence, not a pass for that grid. Do not adopt automatic selection
or a default until the P10G multi-system set has also passed (cell shapes,
projectors, gaps, atom and band counts). Silicon-only timings cannot establish
general performance. Batched/fused stage 2 awaits P10C.

After review, Claude stores measured evidence under `tests/reference`, extends
`scripts/make_readme_figures.py` for any new result type, regenerates the README,
and records an accepted optimization in the log. Until then the comparison
document explicitly marks runtime results pending. The local dashboard generator
can run here; publishing the private Claude Artifact and updating its datasets
requires the orchestrator's dashboard access.
