# ADR 0024: Fixed-density spectra and portable interchange

- Date: 2026-11-03
- Task: P8A-BANDS-DOS
- Status: implemented, proposed for reviewer acceptance

## Numerical boundary

`dft_engine::nscf` accepts the original cell, atoms, species potentials, converged
grid density, requested weighted k points and `ScfOptions`. It creates the same
`ScfProblem` as SCF and evaluates exactly one `step_kpoints` at the supplied
density. The output density is discarded and never mixed or used to rebuild a
second potential. This reuses ionic transforms, Hartree, XC (including the
existing functional selection), per-k bases and nonlocal Hamiltonians, Davidson,
and deterministic ordered workers without copying any potential formula.
The input density fixes the FFT dimensions and must have finite nonnegative
values and the expected electron count. Callers remain responsible for passing
the original physical system and cutoff. The CLI does this in one process.

No changes to `scf.rs`, its private fields, solver code, or existing hot paths
are required. The tradeoff is that the existing fixed-point map also computes
an output density and energy audit that NSCF does not need. This costs time and
memory; no speedup or optimal NSCF memory claim is made. Extracting the shared
potential/solve stages can be a later coordinated SCF refactor, with profiling
and equivalence evidence under the performance policy. New orchestration does
not add allocations, timers or I/O inside Hamiltonian/solver iterations.

The result contains original fractional coordinates and weights, eigenpairs,
per-band residual norms in Hartree, the requested tolerance and a convergence
flag. Backend calculations may wrap k internally; original path coordinates
are preserved for distance and output. Default tolerance is the ScfOptions
1e-8 Ha residual norm. Default extra bands are at least four. The CLI explicitly
requires occupied plus at least four empty bands and increases Davidson's
subspace capacity to at least four times the requested band count.

At implementation start P7PC's per-band tolerance API was not present. Every
requested band is therefore strict, including buffer bands. A degenerate top
multiplet may require more bands or iterations. Errors never produce a success
flag or a partial accepted spectrum. Existing worker errors encode a `k index
N:` prefix in `ScfError::Numerical`; the new adapter converts that boundary
prefix to typed `NscfError::KPoint { index, cause }`. Setup/potential errors are
separate. Tests exercise both construction and solve failures. This conversion
is intentionally confined to one adapter; a future typed SCF worker error
should replace it when that shared API can be changed.

## Geometry, reference and DOS

Paths use the reciprocal metric B=2*pi*A^(-T). Explicit vertices retain their
fractional basis meaning. Named fcc/bcc/sc paths use the ordered standard
primitive conventions of Setyawan and Curtarolo; normalized Gram checks and
positive determinant permit a proper rotation and positive scale, but reject
incompatible primitive/conventional cells. A discontinuous branch has no
accumulated jump distance and is not joined by plotting lines. Subdivision is
an interval count or a maximum reciprocal spacing.

For the spin-unpolarized insulator, occupied bands are 0..N_e/2. The reference
is the maximum occupied eigenvalue over the computed sampling. JSON retains
the absolute offset in Hartree and both edge locations. Overlapping sampled
bands are rejected. There is no Fermi level, smearing, symmetry band tracking,
or physical gap convergence claim.

Gaussian DOS lives in the engine, not CLI numerical code. It uses normalized
k weights and spin multiplicity two, standard-deviation sigma, and an analytic
Gaussian CDF for integrated states. The CLI converts eigenvalues to eV before
calling it, constructs a uniform energy grid with a requested maximum spacing,
and serializes results. Tests use an independent trapezoidal integral to avoid
validating the Gaussian only against its own CDF. Default DOS k grid doubles
each SCF grid dimension and retains the shift; explicit non-Gamma lists require
a DOS grid. Only time reversal is applied, using the existing generator/reducer.

## Interchange and usability

CSV and JSON are the authoritative portable bundle. The CLI also writes
whitespace columns, QE-style gnuplot blocks/total-DOS columns, Grace XY projects,
and VASP-style non-spin-polarized EIGENVAL/total-only DOSCAR layouts. These
layouts were implemented from public format documentation, without QE source
access. They are exports, not claims of QE/VASP execution or universal parser
compatibility. DOSCAR's Fermi field is explicitly a zero VBM-reference placeholder;
EIGENVAL path weights are zero because a line path is not a BZ quadrature.
The how-to links each public specification and describes information lost by
each format. An independent standard-library reader validates the file structure
and every numeric row against CSV, plus counts and occupations.

Band SVG generation embeds the standalone standard-library Python script at
build time and runs `python3 -I -c` without shell interpolation. Python 3 is a
runtime requirement for `bands`; DOS also writes a native SVG. The standalone
script optionally aligns a DOS panel to the bands' VBM using JSON offsets and
verifies the same input, pseudopotentials, cell, functional and density hash.
Both panels share an energy axis. Branches stay disconnected and XML text is
escaped. Missing JSON means CSVs must already share an energy reference.

All output destinations are checked against inputs before computing, and
symlink/nonregular/multiply-linked output files are refused. Multi-file writes
are not transactional; failure can leave old or partial files. Exit status is
authoritative. No new dependency, manifest change, unsafe, UPF copy or git commit
is introduced. Ticket-specific prohibitions leave README, dashboards, artifacts
and tickets untouched; the orchestrator owns their updates.

## Validation

See [how-to evidence](../howto-bands-dos.md#measured-demonstration-and-reproduction)
for exact commands and Si observations. Synthetic tests compare NSCF and SCF
on the SCF grid at 1e-8 Ha, exercise thread counts 0/1/2/4/9 with exact equality,
verify charge/input failures and typed k failures, integrate Gaussian DOS within
1e-6 states/cell and check peak energies, and test skew metric/preset geometry.
CLI tests read CSV, DAT, GNU, AGR, EIGENVAL, DOSCAR independently and parse SVG XML.
A separate ignored release test uses only `DFT_SI_UPF`; it verifies Si Gamma
NSCF/SCF agreement and exact thread 1/4 eigenvalue equality. Real Si reference
files and the figure use the validation UPF of unverified provenance. The
post-merge QE 7.6 validation is reserved for Claude on dev08.

Final verification on 2026-11-03: `scripts/check.sh` exited 0 (build, test, fmt,
clippy all PASS); `python3 -I -m unittest discover -s scripts/tests` exited 0
(116 tests). The ignored Si release command, both real-Si CLI subcommands,
the combined SVG command and the independent reader on the real-Si output
all exited 0. The four checked-in example files total 183876 bytes. Changes
are left uncommitted. An existing unrelated `.DS_Store` change is untouched.
