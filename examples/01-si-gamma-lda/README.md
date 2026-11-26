# Si Gamma LDA energy and eigenvalues

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, from the [QE pseudopotential library](https://pseudopotentials.quantum-espresso.org/legacy_tables/original-qe-pp-library/si). Use a supported UPF v2 file; inspect it before running. The historical local validation file has unresolved provenance/license, so no redistribution or identity with a fresh download is implied.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 01-si-gamma-lda
input=target/tmp/examples/01-si-gamma-lda/input.toml
target/release/baandi validate "$input"
target/release/baandi run "$input" --output-dir results/01-si-gamma-lda
```

Look at `results.json`: `total_energy` (Ha/Ry/eV), eigenvalues, convergence, `scf.driver = "reference"` and verification. `density.bin` and `density.json` contain total density and its grid/checksum. Gamma and 30 Ry are teaching settings; converge the cutoff and k mesh before using energies scientifically.

Laptop runtime budget: **seconds to a few minutes**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
