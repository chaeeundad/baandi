# Structure file to BAANDI input (CIF)

`baandi-convert` turns a structure file (anything ASE reads: CIF, POSCAR,
extended XYZ, XSF, ...) or a Materials Project id into the explicit `[system]`
block BAANDI reads (`--compact` writes the one-atom-per-line `[structure]` form).
It runs on your workstation; the `baandi`/`dft-cli` executable reads only POSCAR
and extended XYZ files itself ([example 15](../15-poscar-input/)) and never
contacts a database, so compute nodes need no Python, ASE or network. See
[Converting structures](../../docs/cli.md#converting-structures).

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, as in [example 01](../01-si-gamma-lda/).
From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).

Step 1, on a workstation with ASE (`pip install ase`, or
`pip install ./python/baandi-convert[ase]` for the `baandi-convert` command):

```sh
mkdir -p target/tmp/examples/14-structure-file
python3 python/baandi-convert/baandi_convert.py examples/14-structure-file/Si.cif \
  --skeleton -o target/tmp/examples/14-structure-file/converted.toml
```

[`Si.cif`](Si.cif) is a hand-written P 1 cell of conventional diamond Si
(a = 5.43 Å, eight sites, occupancy 1). The converter writes the lattice rows
and fractional positions with 17 significant digits, a `[system.provenance]`
table (source path, SHA-256, format, reader and version, conversion time) and,
with `--skeleton`, placeholder `[pseudopotentials]`, `[basis]`, `[xc]` and
`[scf]` sections. [`input.toml`](input.toml) is that result with the UPF name
filled in (its header says how the committed copy was made). Partial
occupancies, overlapping sites, missing cells and non-periodic files are
rejected with a message instead of being guessed.

Step 2, on any machine with the BAANDI executable:

```sh
python3 examples/prepare.py 14-structure-file
input=target/tmp/examples/14-structure-file/input.toml
target/release/dft-cli validate "$input"
target/release/dft-cli run "$input" --output-dir results/14-structure-file
```

`results.json` carries the provenance table unchanged as `structure_provenance`;
it is never used for the physics. The run is an ordinary Gamma-only Si8 LDA
calculation at the 30 Ry placeholder cutoff: a teaching setting, not a converged
result.

Laptop runtime budget: **under a minute to a few minutes**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
