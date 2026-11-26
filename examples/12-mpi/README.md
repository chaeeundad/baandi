# MPI with the same Si input

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, from the [QE pseudopotential library](https://pseudopotentials.quantum-espresso.org/legacy_tables/original-qe-pp-library/si). Use a supported UPF v2 file; inspect it before running. The historical local validation file has unresolved provenance/license, so no redistribution or identity with a fresh download is implied.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
cargo build --release -p dft-cli --features mpi
python3 examples/prepare.py 12-mpi
input=target/tmp/examples/12-mpi/input.toml
target/release/dft-cli validate "$input"
mpirun -np 2 target/release/dft-cli run "$input" --output-dir results/12-mpi
```

This input is identical to example 01. Install/configure a supported MPI toolchain first; see [MPI guide](../../docs/wiki/mpi.md). Rank zero owns output. Keep one thread per rank here and avoid oversubscription. Results need not be bit-identical across rank counts; exact checkpoints require the same ranks/partition. Compare performance only at equal accuracy and the same parallel model/rank count (BAANDI MPI versus QE MPI). No scaling claim is made for this tiny cell.

Laptop runtime budget: **seconds to a few minutes; startup may dominate**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
