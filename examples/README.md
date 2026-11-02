# Input example

Run from the repository root:

```sh
cargo run --offline -p dft-cli -- validate examples/si_lda.toml
cargo run --offline -p dft-cli -- grid examples/si_lda.toml
```

`si_lda.toml` describes the two-atom diamond Si primitive cell with
conventional lattice parameter 5.43 Angstrom. Each `lattice_vectors` row
is one lattice vector, not one matrix row of the internal column-vector
matrix A. The CLI converts Angstrom to Bohr and transposes these rows into
columns through dft-geometry. `lattice_unit = "bohr"` is also accepted.
Atomic positions are fractional, need exactly three finite components,
and are not wrapped. The energy input is in Ry: 30 Ry = 15 Ha.

Only `kpoints = [[0.0, 0.0, 0.0]]` and `functional = "lda-pz"` are accepted.
All structured sections reject unknown keys. Pseudopotential keys are
species labels and their values are paths; every atom needs a matching
entry. The UPF path is an intentional placeholder. These commands parse it without
checking file existence or reading UPF contents. Neither command runs SCF.

The grid report prints the full complex Gamma sphere (both G and -G),
the four-times-energy density sphere, and independently chosen FFT grid.
To override the automatic grid, add `fft_grid = [27, 27, 27]` in `[basis]`.
Each length must be at least 2*M_i+1 for the density sphere; validate also
checks this requirement. Override lengths need not have only factors 2,3,5.
These grid dimensions and counts are not QE reference results.

The project spec's optional SCF/output schema is:

```toml
[scf]
max_iterations = 80
mixing = "linear"
mixing_beta = 0.3
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-6

[output]
directory = "results"
```

Both sections and every field within them are optional. The run command uses
these controls; `mixing = "pulay"` is also supported and is the default.
`nbands`, `extra_bands`, `eigen_tolerance`, and `seed` are optional controls.
See [CLI documentation](../docs/cli.md) for defaults and the complete output schema.

```sh
cargo run --offline -p dft-cli -- run input.toml --output-dir results
```

Supply your own supported NC SLA PZ input path; no real UPF is included here.
Relative UPF paths resolve against the TOML directory. Output paths resolve
against the working directory. The CLI writes results.json, density.bin and
density.json, with checksums and reproducibility metadata. Only engine-certified
convergence returns exit code 0; exhausted budgets write results and return 3.
Synthetic end-to-end verification commands are documented in docs/cli.md.
