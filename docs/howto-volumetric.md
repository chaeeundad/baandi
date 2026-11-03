# View density, ELF and potentials

## Request the files

Add this table to your usual input (the output directory remains the normal
`[output] directory` or `--output-dir`):

```toml
[output.volumetric]
quantities = ["density", "elf", "hartree", "local_potential", "partial_charge"]
formats = ["cube", "xsf", "vasp", "vti"]
kpoint = 1
bands = [1, 1]
```

Choose only the quantities and formats you need. `bands = [2, 4]` sums bands
2, 3 and 4. Band and k indices start at 1. The k index follows the retained list
in results.json, after any time-reversal reduction. Omitted selections mean
first k and first band. Selection options require `partial_charge`; duplicate,
empty or unknown quantities/formats fail. Out-of-range selections fail rather
than silently selecting another state. Existing volumetric files are never
overwritten. Use a fresh directory to repeat an export.

The cube filenames are `density.cube`, `elf.cube`, `hartree.cube`,
`local_potential.cube` and `partial_charge.cube`, with corresponding `.xsf`
files. VASP filenames use `.CHGCAR` for density/partial_charge, `.LOCPOT` for
Hartree/local_potential, and `.ELFCAR` for ELF. VTI uses `<quantity>.vti`.
The existing density.bin output is unchanged. If convergence fails,
fields describe the final returned state, not a converged solution; always check
results.json before interpreting them.

## Open and inspect

In VESTA, open the XSF or cube file using File > Open. XSF includes the periodic
cell as well as atoms. To import another field into an existing structure, use
Edit Data > Volumetric data > Import and choose replacement rather than adding
it to the current field. Use the isosurface controls to select a value, or a
slice/contour view to examine a plane. ELF values such as 0.5 and 0.8 make useful
comparisons. Potentials may need both positive and negative contour levels.
VESTA documents both formats and its import controls in its
[format guide](https://jp-minerals.org/vesta/en/doc/VESTAch18.html) and
[volumetric import guide](https://jp-minerals.org/vesta/en/doc/VESTAch6.html).
XSF can also be used with XCrySDen. Viewer behavior is not tested here.

Cube geometry is Bohr; XSF geometry is Angstrom. Both store scalar density values
in electrons/Bohr^3, potentials in Hartree and ELF without units. Do not apply an
automatic density conversion merely because XSF geometry is Angstrom. To compare
with electrons/Angstrom^3, divide values by (0.529177210903)^3. The units printed
by a viewer may be generic: consult this guide rather than that label.

## What the fields mean

| File stem | Meaning | Units |
| --- | --- | --- |
| density | Total valence electron density; cell integral is the valence electron count | electrons/Bohr^3 |
| elf | Localization indicator: uniform electron gas 0.5, single real orbital 1 above the vacuum screen | dimensionless |
| hartree | Electrostatic repulsion potential from valence density | Ha |
| local_potential | Ionic local + Hartree + selected XC potential | Ha |
| partial_charge | Sum of squared selected orbitals at one k | Bohr^-3 |

Partial charge is not occupation weighted. Its integral is one for one band,
three for a three-band interval, even if bands are empty; the k integration
weight is also omitted. It shows the spatial shape of the selected states.

ELF uses the spin-summed density and occupations, already containing the usual
factor of two. tau=1/2 sum(w f |grad psi|^2), tau_W=|grad n|^2/(8n),
tau_h=(3/10)(3 pi^2)^(2/3)n^(5/3), ELF=1/[1+((tau-tau_W)/tau_h)^2].
At n<=1e-10 electrons/Bohr^3 it is set to zero to avoid meaningless vacuum
ratios. Constant orbital phases and the common Bloch phase cancel; this is not
a current-corrected ELF for arbitrary electromagnetic gauge transformations.

**Pseudopotential caveat:** density and ELF are valence-only. A norm-conserving
pseudopotential smooths the core-region orbitals and omits core electrons. These
plots differ from all-electron density and all-electron ELF near nuclei. Do not
interpret missing core peaks as a physical depletion of core charge.

The Hartree cell average is zero by convention (G=0 removed). Total local
potential retains the ionic finite zero mode and XC offset. Absolute offsets
are convention dependent; compare differences with matched conventions.
Nonlocal pseudopotential operators have no scalar potential representation and
are absent from local_potential.

Cube has N samples per periodic axis and voxel a/N, last index fastest. XSF
has N+1 nodes, repeats periodic boundary planes, and writes the first index
fastest. Do not count the repeated XSF endpoints twice when integrating. General
skew and left-handed cells are preserved; some viewers support fewer cell types.

## Independent check and real-Si example

```sh
python3 -I scripts/read_cube.py results/density.cube
python3 -I scripts/read_cube.py results/elf.cube
```

The JSON includes dimensions, cell vectors, atoms, minimum, maximum and the
scalar integral. Only the density integral is an electron count; ELF's integral
is not. The reader checks the emitted scalar Bohr cube variant, not every cube
extension. The test suite verifies sample round trips at 17 significant digits
and density integral within 1e-10 electrons, including left-handed skew cells.

For an opt-in real-Si run, use your own supported NC LDA Si UPF:

```sh
DFT_SI_UPF=/absolute/path/to/Si.UPF \
DFT_VOLUMETRIC_OUTPUT=/absolute/path/to/new-output \
cargo test --offline --release -p dft-cli --test volumetric real_si_example -- --ignored --nocapture
```

It writes all fields and prints the independently read density integral and ELF
range. It requires a converged run and never commits the pseudopotential. This is an engine/example run, not QE evidence. See the measured run in
[ADR 0021](decisions/0021-volumetric-output.md).

## VASP and ParaView formats

`formats = ["vasp"]` selects a POSCAR-style header followed by a single scalar
grid in first-index-fastest order. Geometry is Angstrom with scale factor 1 and
fractional atom positions, grouped by element. The complete general cell is
preserved, including its handedness.

- `density.CHGCAR` and `partial_charge.CHGCAR` store density times **cell** volume:
  `n[Bohr^-3] * Omega[Bohr^3]`. This is also `n[Angstrom^-3] * Omega[Angstrom^3]`.
  There is no additional grid-size factor. The mean of the stored numbers is the
  electron/probability integral. A CHGCAR reader divides by the cell volume in
  Angstrom^3 to recover density in Angstrom^-3.
- `hartree.LOCPOT` and `local_potential.LOCPOT` store eV, converted from Hartree;
  no volume scaling is applied. The former is Hartree-only, as its filename says.
- `elf.ELFCAR` stores dimensionless ELF, with no volume scaling.

Open these with VESTA's VASP format reader if the filename filter does not select
them automatically. These are visualization exports, **not VASP restart files**:
there are no PAW augmentation occupancies, POTCAR data, spin channels or core
reconstruction. The header cannot carry each atom's valence charge or the k/band
selection. Keep the input TOML with your results. See the public
[CHGCAR](https://vasp.at/wiki/CHGCAR), [LOCPOT](https://vasp.at/wiki/LOCPOT), and
[ELFCAR](https://vasp.at/wiki/ELFCAR) format descriptions. Implementation uses
public format documentation only, not source code from those programs.

`formats = ["vti"]` creates ASCII VTK ImageData for ParaView. Open a `.vti`, press
Apply, and select the named scalar array. Use a Slice or Contour filter for a
plane or an isosurface. Geometry is Angstrom, but array names explicitly specify
atomic scalar units (for example `density_electrons_per_bohr3`). The full-cell
point grid repeats periodic endpoints: there are N+1 points on each axis.

For an axis i, Spacing_i = |a_i|/N_i in Angstrom. The row-major Direction matrix
has columns a_i/|a_i|. Thus `r = Direction * diag(Spacing) * index`; skew cells
are not orthogonalized and left-handed cells are not reflected. Modern,
direction-aware VTK readers/filters are required. Some filters assume orthogonal
or axis-aligned images and may distort a skew cell; the direction matrix must
survive the pipeline. This limitation is documented by the
[vtkImageData API](https://vtk.org/doc/nightly/html/classvtkImageData.html).
The file uses the [VTK XML image layout](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html).

VTI exports only the scalar image: it carries no atom geometry, species, charges,
or electronic-state metadata. Cube includes atoms and their valence charges
but no explicit periodic-cell object; viewers may treat it as a finite box.
XSF includes periodic geometry and atoms but no valence-charge column. None of
these scalar exports carries nonlocal projectors, orbital phases, all-electron
core reconstruction or a restartable SCF state.

The independent test reader in `crates/dft-cli/tests/validate_volumetric.py`
checks VASP geometry, grouped atoms, scaling, ordering and integrated charge,
and parses VTI with Python's XML reader to check the full affine map, scalar
units, sample count, ordering and periodic endpoints. It compares every value
against an independently parsed cube for all five quantities. Skew and
left-handed mixed-element fixtures are tested. These are documented format
checks, not claims that a GUI session was performed.
