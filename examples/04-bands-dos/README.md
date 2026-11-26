# Band path and density of states

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, from the [QE pseudopotential library](https://pseudopotentials.quantum-espresso.org/legacy_tables/original-qe-pp-library/si). Use a supported UPF v2 file; inspect it before running. The historical local validation file has unresolved provenance/license, so no redistribution or identity with a fresh download is implied.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 04-bands-dos
input=target/tmp/examples/04-bands-dos/input.toml
target/release/baandi validate "$input"
target/release/baandi bands "$input" --output-dir results/04-bands-dos/bands
target/release/baandi dos "$input" --output-dir results/04-bands-dos/dos
```

The explicit Gamma–X–L–Gamma path uses fractional reciprocal coordinates of this fcc primitive cell. Each command performs its own SCF then fixed-density NSCF. Inspect `bands.svg`/`dos.svg`, CSV and JSON; energies are relative to the sampled VBM, with its absolute Ha offset in JSON. Other exports: `.dat`, `.agr`, QE-style `bands.dat.gnu`, VASP-style `EIGENVAL`/`DOSCAR`. DOS includes both spins; DOSCAR uses a documented VBM-zero placeholder, not a metallic Fermi energy. These exports do not come from QE/VASP runs. See [format limits](../../docs/howto-bands-dos.md).

Laptop runtime budget: **5–30 minutes**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
