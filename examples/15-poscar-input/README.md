# POSCAR input with selective dynamics

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, as in [example 01](../01-si-gamma-lda/).
From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).

```sh
python3 examples/prepare.py 15-poscar-input
input=target/tmp/examples/15-poscar-input/input.toml
target/release/baandi validate "$input"
target/release/baandi relax "$input" --output-dir results/15-poscar-input
```

[`input.toml`](input.toml) names the structure with `[structure] file = "POSCAR"`
(relative to the input file; `prepare.py` copies it next to the prepared input).
The executable reads the [`POSCAR`](POSCAR) natively: VASP 5 layout, scale
factor 5.43 Å, Direct coordinates, and `Selective dynamics`. The first atom is
`F F F`, so `relax` keeps it in place; `T`/`F` along a lattice direction becomes
a Cartesian constraint (here the cell is cubic, so each direction is an axis).
Extended XYZ files (`Lattice=`, `Properties=`, `move_mask`, `initial_magmoms`)
are read the same way; CIF and other formats go through `baandi-convert`
([example 14](../14-structure-file/)).

`results.json` has a `structure` block with the lattice in Bohr, the
fractional positions, the per-atom `fixed` flags and the source
(`kind = "file"`, path, format and SHA-256 of the POSCAR), and
`relaxation.constraints` lists the frozen atom. See
[Structure input](../../docs/cli.md#structure-input).

Laptop runtime budget: **10–60 minutes or more for relaxation**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
