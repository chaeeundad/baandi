# Si variable-cell relaxation

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, from the [QE pseudopotential library](https://pseudopotentials.quantum-espresso.org/legacy_tables/original-qe-pp-library/si). Use a supported UPF v2 file; inspect it before running. The historical local validation file has unresolved provenance/license, so no redistribution or identity with a fresh download is implied.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 08-cell-relax
input=target/tmp/examples/08-cell-relax/input.toml
target/release/baandi validate "$input"
target/release/baandi relax "$input" --output-dir results/08-cell-relax
```

This slightly expanded Si primitive cell relaxes isotropically at zero external pressure with fractional ions fixed. Use `mode = "ions_cell"` to also move ions. Inspect per-frame cell/stress/enthalpy, `final.cif`, `trajectory.extxyz` and `results.json`. Frozen Miller segments rebuild every five moves and at the end: inspect rebuilt pressure and Pulay warnings; a frozen-basis minimum is not a cutoff-converged equilibrium. Verification remains on; default stress threshold is 5e-11 Ha. See [cell guide](../../docs/user/cell-relaxation.md).

Laptop runtime budget: **10–60 minutes or more**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
