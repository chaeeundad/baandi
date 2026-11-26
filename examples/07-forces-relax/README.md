# Displaced Si8 forces, stress and ion relaxation

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, from the [QE pseudopotential library](https://pseudopotentials.quantum-espresso.org/legacy_tables/original-qe-pp-library/si). Use a supported UPF v2 file; inspect it before running. The historical local validation file has unresolved provenance/license, so no redistribution or identity with a fresh download is implied.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 07-forces-relax
input=target/tmp/examples/07-forces-relax/input.toml
target/release/dft-cli validate "$input"
target/release/dft-cli run "$input" --output-dir results/07-forces-relax
target/release/dft-cli relax "$input" --output-dir results/07-forces-relax/relaxed
```

One of eight Si atoms is displaced by 0.002 fractional x. `run` reports forces (Ha/Bohr) and tensile-positive stress (Ha/Bohr^3); pressure is −trace(stress)/3. Verification stays on. `relax` moves ions in a fixed cell; inspect `relax.log`, `trajectory.extxyz` (eV, Angstrom), `POSCAR`, `final.toml` and the relaxation block in `results.json`. Exit 3 means the step budget was exhausted, not convergence. Only unpolarized fixed occupations are supported by the CLI optimizer. See [relaxation](../../docs/user/relaxation.md).

Laptop runtime budget: **10–60 minutes or more for relaxation**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
