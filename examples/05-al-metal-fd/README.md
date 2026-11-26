# Al metal with Fermi–Dirac occupations

UPF: `Al.upf` from [PseudoDojo](https://www.pseudo-dojo.org/), NC-SR v0.4 PBE **standard**, UPF format ([archive](https://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard_upf.tgz)). Extract outside the repository and place this file directly in `PSEUDO_DIR`; retain its license and attribution.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 05-al-metal-fd
input=target/tmp/examples/05-al-metal-fd/input.toml
target/release/baandi validate "$input"
target/release/baandi run "$input" --output-dir results/05-al-metal-fd
```

One Al atom has three valence electrons with the specified file, so use FD rather than fixed insulating occupations. Inspect internal energy E, entropy contribution TS, free energy F=E−TS, chemical potential, occupations and top-band coverage in `results.json`. Repeat with a denser mesh, smaller width and at least four more bands before drawing physical conclusions. This is limited-validation metallic functionality, not an accuracy-qualified Al reference.

Laptop runtime budget: **2–20 minutes**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
