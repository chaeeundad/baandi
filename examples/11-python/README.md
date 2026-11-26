# Python batch, forces and stepwise API

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, from the [QE pseudopotential library](https://pseudopotentials.quantum-espresso.org/legacy_tables/original-qe-pp-library/si). Use a supported UPF v2 file; inspect it before running. The historical local validation file has unresolved provenance/license, so no redistribution or identity with a fresh download is implied.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 11-python
input=target/tmp/examples/11-python/input.toml
target/release/dft-cli validate "$input"
python3 -m venv target/tmp/examples/python-env
. target/tmp/examples/python-env/bin/activate
python -m pip install ./python/baandi tomli
python examples/11-python/run.py "$input"
```

`run.py` reads the prepared input (Python 3.9+, using `tomli` before 3.11) and passes UPF **text**, Bohr lattice vectors, fractional positions and `Options` to `Session`. It prints energy in Ha, forces in Ha/Bohr, then creates a fresh session and advances one reference map per `step()`. No output files are written by this script. Check convergence before querying derivatives. The pip command builds the native binding in release mode using maturin and the pinned Rust toolchain. See [Python API](../../python/baandi/) and [binding validation](../../crates/baandi-py/VALIDATION.md).

Laptop runtime budget: **1–10 minutes for both solves**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
