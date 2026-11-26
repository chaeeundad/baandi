# Si Monkhorst–Pack mesh with symmetry

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, from the [QE pseudopotential library](https://pseudopotentials.quantum-espresso.org/legacy_tables/original-qe-pp-library/si). Use a supported UPF v2 file; inspect it before running. The historical local validation file has unresolved provenance/license, so no redistribution or identity with a fresh download is implied.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 02-si-kmesh
input=target/tmp/examples/02-si-kmesh/input.toml
target/release/baandi validate "$input"
target/release/baandi run "$input" --output-dir results/02-si-kmesh
```

The half-shifted 2×2×2 Monkhorst–Pack grid enables compatible spatial symmetry. Inspect `results.json.symmetry` for operations, irreducible weights, multiplicities and the full-to-reduced map. Automatic FFT sizing respects symmetry; compare symmetry on/off at an identical FFT grid.

Laptop runtime budget: **1–10 minutes**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
