# Si PBE and NLCC

UPF: `Si.upf` from [PseudoDojo](https://www.pseudo-dojo.org/), NC-SR v0.4 PBE **standard**, UPF format ([archive](https://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard_upf.tgz)). Extract outside the repository and place this file directly in `PSEUDO_DIR`; retain its license and attribution.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 03-si-pbe
input=target/tmp/examples/03-si-pbe/input.toml
target/release/baandi validate "$input"
target/release/baandi run "$input" --output-dir results/03-si-pbe
```

PseudoDojo files may contain nonlinear core correction (NLCC); current NC support includes it. Inspect the header rather than inferring NLCC from a filename. The core contributes to XC, not the valence electron count. Match the PBE functional; do not enable a functional override merely to silence a mismatch. Inspect energy, verification and executed atomic-start provenance. The 40 Ry cutoff needs a convergence study.

Laptop runtime budget: **1–10 minutes**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
