# ADR 0021: Post-SCF volumetric scalar output

Status: Implemented. The orchestrator approved moving the existing dft-physics
CLI dev-dependency to dependencies; no new package or lockfile change is needed.
Normal `run` invokes preflight and exports after the final SCF result, using four
additive lines in run.rs (module declaration and two calls).

## Decision

Compute visualization fields once, after SCF, using its final density and final
per-k eigenvectors. Reconstruct basis ordering from the same cutoff, cell and
canonical k coordinates. Never solve the SCF again. The physical implementation
lives in dft-physics/src/volumetric.rs; serialization, selection and the complete
118-element symbol table live in dft-cli/src/volumetric_out.rs. The symbol table
uses [IUPAC symbols](https://iupac.org/what-we-do/periodic-table-of-elements/).
Use the UPF element symbol, not an arbitrary species label. Unknown symbols fail.

All orbital reductions run serially in final k order, band order, Cartesian
component order. Occupations are spin-summed and applied once; weights are
applied once. There are no new dependencies, unsafe blocks or hot-path changes.

With psi=exp(ik.r) sum_G c_G exp(iG.r)/sqrt(Omega), transform i(k+G)c_G for each
Cartesian gradient. The omitted common exp(ik.r) cancels in squared norms and
products. tau=1/2 sum(w f |grad psi|^2). For ELF, compute grad n through the
orbital product rule 2 Re sum(w f psi* grad psi). This is analytically identical
to differentiating density and retains relative accuracy in tiny density tails.
A spectral derivative of sampled density was rejected after the polynomial
Gaussian fixture showed an ELF error 6.515712578617494e-4 near vacuum: FFT absolute
roundoff is amplified by n^(5/3). Product-rule evaluation gives zero measured ELF
error in that fixture without changing the 1e-10 tolerance or density screen.

ELF=1/(1+chi^2), chi=(tau-|grad n|^2/(8n))/tau_h and
tau_h=(3/10)(3 pi^2)^(2/3)n^(5/3). At n<=1e-10 electrons/Bohr^3, ELF=0.
Constant orbital phase changes leave these quantities invariant. Arbitrary
position-dependent electromagnetic gauge changes do not: the requested formula
has no current correction. In particular a current-carrying single complex
orbital need not have ELF=1. A uniform electron gas means tau=tau_h, not merely a
constant density with an arbitrary orbital population.

Hartree is evaluated from final output density with its G=0 coefficient zero.
The total scalar potential is ionic + Hartree + selected LDA/PBE XC evaluated
on that density. Ionic finite alpha_z G=0 and the XC constant are retained;
nonlocal projectors cannot be exported as a scalar. These are output-density
potentials, not the previous mixed input potential when SCF is unconverged.
The Hartree routine's canonical skew-grid complex Nyquist field is projected to
its real part only at the scalar output boundary, consistent with the engine.

Partial charge sums normalized |psi|^2 over a one-based inclusive band interval
at one retained k point. No occupation or k weight is applied: each band
integrates to one, even an empty band. No renormalization masks coefficient errors.

## Formats and validation

Cube: two comments, positive atom and axis counts, origin zero, Bohr coordinates,
voxel vectors a_i/N_i, atomic number and UPF valence charge, last index fastest.
Each scalar is written with 17 significant decimal digits (stronger than the
required 1e-12 relative precision); rows are whitespace separated with at most
six values per line. General skew and left-handed cells retain their orientation.
See the [cube layout specification](https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html).

XSF: CRYSTAL, PRIMVEC and PRIMCOORD geometry in Angstrom, plus one named
BEGIN_BLOCK_DATAGRID_3D block. The general grid has N_i+1 nodes and full-cell
spanning vectors, repeating first samples at last faces/edges/corners. X varies
fastest. Scalar units stay atomic, including density in electrons/Bohr^3.
[VESTA's format documentation](https://jp-minerals.org/vesta/en/doc/VESTAch18.html)
specifies geometry in Angstrom and equal primitive and grid spanning vectors.

The independent standard-library Python reader validates the scalar Bohr cube
subset, counts, finite numbers, atomic numbers and determinant-based integral.
Tests compare every parsed value, atom, cell vector and electron integral for
both handednesses; XSF tests check endpoint duplication and first-index-fastest
ordering. A two-k synthetic SCF verifies all twenty files byte-identically for
thread budgets 1 and 4. It is not a physical Si or QE reference.

The opt-in release real-Si test requires a user-supplied NC LDA UPF and a fresh
output directory. No UPF is stored in the repository. GUI viewing has not been
performed; automated reader/structure validation is the interoperability evidence.

## Verification record (2026-11-03)

Commands use the existing offline toolchain at RUSTUP_HOME=/private/tmp/p0c-rustup
and /opt/homebrew/opt/rustup/bin on PATH.

- `cargo test --offline -p dft-physics --test volumetric -- --nocapture`: exit 0,
  four tests. Uniform-gas ELF exactly 0.5; spin-paired normalized periodic
  Gaussian-like orbital max ELF error 0, charge 1.9999999999999922 electrons;
  Bloch plane-wave tau max error 0 Ha/Bohr^3, partial charge integral
  1.0000000000000007. ELF tolerance remains 1e-10 above the vacuum threshold.
- `cargo test --offline -p dft-cli --test volumetric -- --nocapture`: exit 0,
  three tests passed, one real-Si test ignored. Both handednesses of the skew
  cube gave 7.9999999999999991 electrons against expected 8; all samples survived
  the independent Python read bitwise. XSF ordering and periodic faces passed.
  Synthetic two-k outputs: ten files byte-identical for thread budgets 1 and 4.
- `RUSTUP_HOME=/private/tmp/p0c-rustup CARGO_NET_OFFLINE=true scripts/check.sh`:
  final exit 0, build/test/fmt/clippy all PASS. An earlier run exited 1 for two
  writer range-loop Clippy warnings; these were corrected before the final run.
- `python3 -I -m unittest discover -s scripts/tests`: exit 0, 115 tests.
  Existing unrelated ResourceWarnings remain; no Python test failed.

At that initial stage no GUI session or real-Si run had been performed.
No new QE comparison was produced.
The initial verification above preceded runtime integration and the added VASP/VTI
formats; see the follow-up record below. README, artifacts, workspace manifest,
lockfile and tickets were not edited. Changes remain uncommitted for review.

## Additional standard viewer formats

Per the orchestrator's output-format instruction, `vasp` and `vti` join cube and
XSF in the same table, without dependencies. Implementations derive from public
format documentation only; no source from VASP, QE or other simulation codes was
consulted or copied.

VASP scalar layout: [POSCAR](https://vasp.at/wiki/POSCAR) header, Angstrom cell,
scale 1, exact element symbols and grouped fractional positions, then a blank
line, N1 N2 N3, and first-index-fastest scalar data. Density/partial charge use
CHGCAR scaling n*Omega, so the sample mean equals the integral (no extra Ngrid
factor). Potentials use LOCPOT eV units and ELF uses unscaled ELFCAR values.
These choices follow the public [CHGCAR](https://vasp.at/wiki/CHGCAR),
[LOCPOT](https://vasp.at/wiki/LOCPOT), and [ELFCAR](https://vasp.at/wiki/ELFCAR)
format definitions. The five filenames distinguish normalization. No PAW
augmentation or spin blocks are invented; files are visualization-only.

VTK ImageData: XML ASCII Float64 PointData, first index fastest, extents 0..Ni
including repeated periodic endpoints. Geometry is Angstrom. Spacing is
|a_i|/Ni, Direction is row-major with normalized lattice columns. Consequently
Direction*diag(Spacing) maps integer sample indices to actual Cartesian voxel
positions, including shear and negative determinant. Scalar arrays retain
atomic units in their names. No atom objects are embedded. See the public
[XML layout](https://docs.vtk.org/en/latest/vtk_file_formats/vtkxml_file_format.html)
and [vtkImageData transform API](https://vtk.org/doc/nightly/html/classvtkImageData.html).
Readers/filters that ignore direction or assume orthogonality cannot correctly
process skew exports; that limitation is explicit in the viewer guide.

Independent standard-library readers in tests/validate_volumetric.py parse
VASP and XML VTI, compare all values with the Python cube reader, reconstruct
lattice vectors, test handedness and mixed-species atom grouping, and check
unit scaling and periodic endpoint duplication. A normal CLI subprocess test
checks that all formats are emitted, existing outputs are protected, and the
cube density is bit-identical to the unchanged density.bin. The real-Si ignored
test also invokes the normal CLI rather than a private exporter.

## Runtime integration and real-Si verification (2026-11-03)

`cargo test --offline -p dft-cli --test volumetric -- --nocapture` after runtime
integration and the new formats: exit 0, five passed and the real-Si test ignored.
The synthetic two-k SCF now produces twenty byte-identical files for thread
budgets 1 and 4. Normal CLI output is checked independently for all five fields;
its density.cube samples equal density.bin bitwise. VASP/VTI checks also cover
mixed elements and both handednesses of a skew cell. No dependency was added;
only the existing dft-physics manifest line moved into runtime dependencies.

Executed real-Si command (exit 0, one ignored test explicitly run and passed):

```sh
PATH=/opt/homebrew/opt/rustup/bin:$PATH \
RUSTUP_HOME=/private/tmp/p0c-rustup \
DFT_SI_UPF=/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
DFT_VOLUMETRIC_OUTPUT=/private/tmp/p8b-real-si-volumetric \
cargo test --offline --release -p dft-cli --test volumetric real_si_example -- --ignored --nocapture
```

The test invoked normal `dft-cli run` on its generated input. Two-atom diamond
Si, a=5.43 Angstrom, Gamma, LDA-PZ, cutoff 30 Ry, 749 plane waves, grid 25^3,
12 SCF iterations, `converged=true`. UPF was read in place, never copied:
SHA-256 `e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Input SHA-256 `c2cbd52c2a2d87af1afbada2e4bd25bc27b9c9a6d56446c2ab6eb6b6228ac024`.
The input is examples/si_lda.toml with the given UPF path and all five quantities
and all four formats requested (default band=1, k=1).

Independent-reader measurements:

| Field | Measurement | Units |
| --- | --- | --- |
| density | integral 8.000000000000004, expected 8 | electrons |
| density | min 2.9606900036306047e-5, max 0.12247174657531103 | electrons/Bohr^3 |
| ELF | min 2.2842676730263424e-12, max 0.999999208442029 | dimensionless |
| partial_charge | integral 1.0000000000000029, expected 1 | probability |
| Hartree | min -0.43259383508650245, max 0.528974811045778 | Ha |
| total local | min -6.8818338630188105, max 0.1477505201343282 | Ha |

VASP density sample mean was 8.000000000000002 and partial-charge mean
1.0000000000000027. Every VASP and VTI scalar was independently compared with
its cube counterpart; geometry, unit conversion, ordering and endpoint checks
passed. ELF's spatial integral is not an electron count and is not used as one.

Files remain outside the repository under `/private/tmp/p8b-real-si-volumetric`:

- `density.cube`, `density.xsf`, `density.CHGCAR`, `density.vti`
- `elf.cube`, `elf.xsf`, `elf.ELFCAR`, `elf.vti`
- `hartree.cube`, `hartree.xsf`, `hartree.LOCPOT`, `hartree.vti`
- `local_potential.cube`, `local_potential.xsf`, `local_potential.LOCPOT`, `local_potential.vti`
- `partial_charge.cube`, `partial_charge.xsf`, `partial_charge.CHGCAR`, `partial_charge.vti`
- Existing `density.bin`, `density.json`, `results.json`, plus the test's
  `input.toml`, `run.stdout`, `run.stderr` for reproduction.

This validates the export pipeline on a physical Si input; it is not an
all-electron reconstruction or a new QE comparison. GUI opening and ParaView
filter behavior were not exercised. External temporary output is not durable
repository evidence; the command, input definition and measured results above
permit reproduction with the same read-only UPF.

Final full-workspace verification after integration:

- `RUSTUP_HOME=/private/tmp/p0c-rustup CARGO_NET_OFFLINE=true scripts/check.sh`:
  exit 0; build, test, fmt and clippy all PASS. The preceding run exited 1 on a
  test-only constant chunks_exact Clippy lint; changing it to as_chunks fixed
  the lint without changing the numerical calculation or real-Si output.
- `python3 -I -m unittest discover -s scripts/tests`: exit 0, 115 tests passed.
- `git diff --check`: exit 0.

All implementation files remain uncommitted. An unrelated existing `.DS_Store`
modification was left untouched. No UPF or generated real-Si export was added to
this repository; Cargo.lock is unchanged.
