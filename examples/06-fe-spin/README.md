# Fe collinear spin with a fixed moment

UPF: `Fe.upf` from [PseudoDojo](https://www.pseudo-dojo.org/), NC-SR v0.4 PBE **standard**, UPF format ([archive](https://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard_upf.tgz)). Extract outside the repository and place this file directly in `PSEUDO_DIR`; retain its license and attribution.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 06-fe-spin
input=target/tmp/examples/06-fe-spin/input.toml
target/release/dft-cli validate "$input"
target/release/dft-cli run "$input" --output-dir results/06-fe-spin
```

FD avoids imposing an insulating gap on metallic Fe. The total moment is constrained to 2 mu_B per primitive cell, not a predicted ground-state moment. Inspect spin/k/band axes, channel counts, moment, free energy and verification in `results.json`, and `density-up.cube`, `density-down.cube`, `density-total.cube`, `magnetization.cube`. Magnetic time reversal is disabled. **G4D:** removing the moment constraint creates free-moment PBE, which is unqualified and routes to the strict oracle; that fallback does not guarantee convergence. See [spin guide](../../docs/user/collinear-spin.md).

Laptop runtime budget: **5–60 minutes or more**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
