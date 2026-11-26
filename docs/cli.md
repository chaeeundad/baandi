# Command-line SCF

> The command is `baandi` (`target/release/baandi`). The Cargo package is still `dft-cli` (`cargo build -p dft-cli`), and a `dft-cli` binary of the same program remains as a deprecated alias for existing scripts; it will be removed in a later release. Result files keep `"tool": {"name": "dft-cli"}` for schema stability.

```sh
cargo run --offline -p dft-cli -- run input.toml [--output-dir DIR]
```

`run` uses the production **reference driver**: real Gamma or complex k points,
NC LDA-PZ/LDA-PW/PBE, optional NLCC, fixed occupations or eligible collinear/FD
models. The sole automatic oracle exception is free-moment collinear FD PBE
(G4D); `BAANDI_REFERENCE_STRICT=on` explicitly selects StrictOracle with precedence.
See [ADR 0047](decisions/0047-reference-driver-single-path.md) for the accepted
contract and [UPF support](upf-support.md) for parser/functional restrictions.
PAW and ultrasoft inputs reject. No pseudopotentials are shipped in the repository.

Use the schema in [si_lda.toml](../examples/si_lda.toml). Each input
`system.lattice_vectors` row is a vector; it is transposed to matrix columns
internally. `lattice_unit` is `bohr` or `angstrom`. Atoms have a species label
and fractional coordinates. Relative pseudopotential paths resolve against the
input file's directory. All structured tables reject unknown fields.

Optional SCF fields and defaults:

| Key | Default | Meaning |
| --- | --- | --- |
| `max_iterations` | 200 | positive integer budget |
| `mixing` | `"broyden"` | `"plain"`, `"broyden"`, or alias `"pulay"`; history 8 |
| `mixing_beta` | 0.7 | finite coefficient in (0, 1] |
| `nbands` | omitted | occupied + `extra_bands`; explicit value takes precedence |
| `extra_bands` | 4 | unoccupied bands when `nbands` is omitted |
| `seed` | 2026 | unsigned 64-bit deterministic initialization seed |

Reference accuracy is `scf.standard_threshold_ha` (default 5e-9 Ha;
5e-11 Ha for stress/variable-cell inputs), with `scf.verify = "fixed_point"`
by default. The reference mixer is Broyden, beta 0.7, history 8 (`1..=25`);
`plain` selects its plain update and `pulay` is an alias for Broyden. Fixed
unpolarized occupations require an even valence electron count; spin sectors
and FD occupations have their own number/moment contracts. The stop estimator
and independent verification are distinct from the frozen oracle's density,
energy and eigen-residual certificate. A threshold is not a proof of true error.

Removed keys/values and the oracle ladder compatibility pin are listed under
[single-path controls](#single-scf-path-and-fixed-point-verification-adr-0043-default-since-2026-11-16).

`[output] directory = "results"` is the default. Relative output directories
resolve against the working directory, including `--output-dir`, which overrides
the TOML value. Directories are created. Previous `results.json`, `density.json`
and `density.bin` are overwritten; other files are untouched. Before computing,
the CLI rejects any output destination that aliases the input TOML or any UPF,
including symlinks and (on Unix) hard links. Do not modify input files or output
links concurrently with a run. A failed run before output leaves prior outputs
intact, so always check the process exit code. Writes are not a multi-file atomic
transaction; an I/O failure may leave partial output.

`scf.threads` sets the thread budget (default 1); MPI is an optional build
feature. Reference Gamma and complex-k/spin/FD paths use distributed numerical
work under MPI. Inspect recorded driver/backend/parallel metadata rather than
inferring the model from rank count. Strict-oracle replication limits do not
describe production reference execution. Formatting and I/O occur outside maps.

## Current result records

`run` writes `results.json`, density metadata/binary and requested standard-format
outputs. Inspect `scf.driver`, `driver_reason`, `driver_source`, `single_path`,
`execution`, verification and the `reference_scf` block. Scalar reference records
include per-map dr2/ethr (Ry), estimator/energy (Ha), solve work and map-boundary
clocks. Electronic spin/FD output uses schema version 2,
`output_profile = "reference_spin_fd_v1"` and `reference_scf.schema =
"reference-spin-fd-1"`; it lacks the scalar map-boundary clock fields. Missing
clocks mean unavailable measurements. The old scalar schema below is retained
as a legacy/oracle field reference, not a complete reference-output schema.
Preflight and executed start-fallback reason text can disagree (dev08 gen27);
prefer the executed start record when explaining what ran.

## Checkpoint/restart v2

```sh
baandi checkpoint input.toml boundary.cp --iteration 3 --kind exact
baandi restart input.toml boundary.cp result.cp
dft-cli status result.cp
baandi checkpoint input.toml evaluated.cp --kind result
```

These default to the reference stepper and v2. `restart` writes a ResultExport
**container**, not a results directory; `<output>.results.json` carries routing
metadata only. Export summaries contain energy, convergence, counts and hashes,
not normal scientific output arrays. Exact restart requires identical input bytes,
UPFs, build/backend/platform, threads, ranks and partition. Warm v2 import and
terminal exact capture reject. v1 is legacy/oracle only, with no automatic v2
conversion. See [format, safe points and limits](checkpoint-format.md).

## Exit status

| Code | Meaning |
| --- | --- |
| 0 | engine-certified convergence and successful output writes |
| 1 | invalid TOML, configuration, geometry, electron/band count |
| 2 | malformed/unsupported UPF, invalid UPF encoding, incompatible functional |
| 3 | SCF budget exhausted; final evaluated state and history written, `converged=false` |
| 4 | diagonalization/numerical error or divergence |
| 5 | filesystem/I/O failure, including missing UPF or protected-input collision |

All run failures emit their cause on stderr. Missing UPF is code 5 because it is
an access failure, not a judgment about file content. Numerical failures do not
produce new result files. These codes describe `run`; legacy commands retain
success/failure status. Clap's usage errors follow clap conventions.

## Legacy/oracle results.json schema version 1

This section describes the historical scalar serializer. Its option fields are
output provenance, not accepted current TOML keys. Ladder/Pulay diagnostic
fields below describe legacy artifacts; the legacy experiments were removed
in R10/R11, see git history. Use the current result records above for production.

Pretty JSON uses the following top-level keys:

| Key | Contents |
| --- | --- |
| `schema_version` | integer 1 |
| `tool` | `name`, `version`, `git_commit`, `rustc_version`; std-only build script captures the last two, with `unknown` fallback |
| `timestamp` | UTC RFC 3339, whole seconds |
| `input_sha256` | lowercase hex SHA-256 of original TOML bytes |
| `species` | ordered by name: `name`, canonical absolute `path`, byte `sha256`, `z_valence`, original `functional`, `n_projectors` |
| `cell` | `lattice_columns_bohr`, `volume_bohr3`, `reciprocal_columns_inverse_bohr` |
| `atoms` | `species`, original `fractional`, corresponding `cartesian_bohr` |
| `cutoff` | `ecutwfc_ry`, `ecutwfc_ha` |
| `basis` | `n_plane_waves`, `n_density_g`, `fft_grid`, `fft_grid_overridden` |
| `electrons`, `bands`, `occupations`, `xc` | electron count, actual band count, electron occupations, `lda-pz` |
| `scf` | all actual engine option values, described below |
| `convergence` | `converged`, `iterations`, `last_iteration`, full `history`; each mixed history entry has `pulay` (SCF-TRAJECTORY diagnostic, never feeds back): `differences`, `kept`, `gram_condition` (null = singular), `coefficients` (oldest difference first), `residual_l2_sum`, `predicted_ratio` (fit-metric norm of R - dR gamma over that of R), `step_l2_sum`, `coulomb_predicted_ratio` and `coulomb_refit` (`BAANDI_MIX_COULOMB_GUARD`); `work.cutoff_ladder.coarse_history` holds the coarse-stage maps in the same format (energies at the coarse cutoff) |
| `effective_start` | the start each stage actually used: with a used cutoff ladder `[{"stage":"coarse","density":"atomic"/"uniform","wavefunctions",...,"atomic_orbitals_used","atomic_start_h_applications","atomic_fallback_reason"}, {"stage":"fine","density":"prolonged coarse density","wavefunctions":"prolonged coarse orbitals"}]`, otherwise one `"stage":"single"` entry. The requested `start_density`/`start_wavefunctions` act on the FIRST stage; the top-level `atomic_start` block is the full-cutoff session's and reads "supplied initial vectors retained" under the ladder by design |
| `energies` | all audit terms, each `{ha, ry}` |
| `eigenvalues` | arrays `ha`, `ev` in band order |
| `timings` | `phases`, `setup_s`, `engine_setup_s`, `engine_total_s`, `total_s`, `fft_hpsi` availability note; P48F additive blocks: `certificate` (blocks, wall_s, screen_s, phases_s already inside `phases`, outside_phases_s), `atomic_start_s` (inside `phases.potential_s`), `steps_s` (wall of the SCF steps), `ladder_coarse` (null without a used cutoff ladder; setup_s, session_s, prolong_s, total_s, coarse phases, atomic_start_s, steps_s, certificate_s) and `closure` (total = setup + phases + certificate outside phases + step overhead + ladder coarse + remainder) |
| `threading` | `threads=1`, backend description |
| `memory_estimate` | array payload bytes and unavailable peak-RSS note |
| `units` | units summary |

Matrix arrays are row-indexed matrices with vectors in **columns**,
`A[row][column]`; reciprocal vectors include 2*pi. Atomic positions in output
are unwrapped input positions; the engine uses periodically equivalent wrapped
positions. `scf` records `cutoff_hartree`, `nbands` (null means automatic),
`extra_bands`, `grid_override`, `mixing` (`method`, `beta`, and Pulay `history`),
`density_residual_tolerance`, `energy_tolerance_ha`, `eigen_tolerance`,
`initial_eigen_tolerance` (1e-5 Ha), `eigen_tolerance_decay` (0.5),
`max_iterations`, `davidson_max_iterations` (400), `davidson_subspace` (48),
`seed`, `divergence_factor` (20), `divergence_patience` (3). The top-level
`bands` is the resolved band count. The scheduled eigen tolerance has a floor
of 0.1 times `eigen_tolerance`.

`scf.switches` (also in the failure record of a diverged or unconverged run)
is the switch manifest: one object per `BAANDI_*` environment switch with
`value` (effective value: string, number, boolean, or null), `source`
(`default`, `env`, `invalid -> default`, `invalid -> off` for `BAANDI_ORTHO`,
`invalid -> error` for input overrides), `default`, `allowed` and `meaning`.
The table lives in `crates/dft-numerics/src/switches.rs`; every switch is read
once per process there and nowhere else.

Each history/last-iteration record has `iteration`, `energies`,
`delta_energy_ha` (null initially), `density_residual`,
`max_eigen_residual_ha`, `requested_eigen_tolerance_ha`, `davidson_iterations`,
`relative_charge_error`, `relative_input_charge_error`,
`min_density_electrons_per_bohr3`, `mixing_fraction`, `mixer_reset`, `timings`,
and `fixed_point_certificate` (null if not attempted). Certificates contain
`density_residual`, `occupied_eigenvalue_change_ha`, `davidson_iterations`,
and `timings`. A certificate can be present on an unsuccessful candidate;
only the top-level convergence Boolean indicates certified completion.

Energy keys are `kinetic`, `local`, `nonlocal`, `hartree`, `xc`, `ewald`,
`total`, and diagnostics `local_zero_mode`, `band_sum`, `double_counting`,
`band_form`, `harris_foulkes`, `functional_difference_bound`. Only the first
six terms sum to `total`; diagnostics must not be added again.

All phase objects have `potential_s`, `diagonalization_s`, `density_s`,
`total_energy_s`, `mixing_s`. Aggregate phases include certificate work.
FFT and H-psi are included in enclosing phases and not separately measured.
CLI setup includes input reading and engine construction; engine setup is a
subset, not another additive term. Total is measured through report preparation,
excluding console/file writes. Engine total includes engine setup and SCF.

Memory keys are `basis_bytes`, `density_basis_bytes` (vector payload sizes),
`density_grid_bytes` (8*Ngrid), `one_complex_grid_bytes` (16*Ngrid),
`orbital_bytes` (16*Npw*Nbands), `peak_rss_bytes` (null), and `note`.
They exclude scratch arrays, projectors and allocator overhead, and are **not**
a measured working-set peak. Use the measurement helper below for peak RSS.

Units not explicit in key names: `density_residual` and its tolerance are
Bohr^(-3/2); `eigen_tolerance`, `initial_eigen_tolerance`, `cutoff_hartree` are
Ha; `z_valence`, `electrons`, occupations and normalization integrals are electron
counts; counts/indices/fractions, seed, relative errors, beta and decay factors
are dimensionless. Timings are seconds; memory is bytes. Ha-to-eV uses the
uom-derived dft-units CODATA 2018 factor. Fractional coordinates are dimensionless.

Determinism comparisons remove `timestamp` and every `timings` object recursively
(including last-iteration and certificate copies), then compare pretty-serialized
JSON bytes. No numerical or convergence keys are excluded. Density bytes match
exactly. Cross-build/platform bitwise equality is not promised. Git metadata
identifies HEAD, not whether the source tree had uncommitted changes.

## Density interchange

`density.bin` contains headerless IEEE-754 f64, little-endian, C order:
`flat=(i*n2+j)*n3+k`. **k, the last index, varies fastest**. There are exactly
`n1*n2*n3` samples in electrons/Bohr^3, representing the final evaluated output
density, even on budget exhaustion.

`density.json` contains `schema_version=1`, `fft_grid`, `lattice_columns_bohr`,
`volume_bohr3`, `units`, `ordering`, `endianness`, `dtype=f64`, `normalization`
(`expected_electrons`, `integral_electrons`), and SHA-256 `sha256` of the binary.
The integral is `sum(samples)*volume/Ngrid` and equals N_e to numerical roundoff.

## Verification and measurements

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export RUSTUP_HOME=/private/tmp/p0c-rustup  # this sandbox's offline toolchain
scripts/check.sh
cargo test --offline -p dft-cli --test all run:: -- --nocapture
cargo test --offline --release -p dft-cli --test all run:: -- --ignored --nocapture
```

The fast suite covers a small synthetic cell, option/date conversion, all exit
classes, path overrides, overwrite protection, checksums, independent density
reading, direct-engine comparison and deterministic outputs. The ignored release
test covers a two-atom fcc cell, synthetic Z=4 Gaussian UPF with two s and one p projector,
8 Ry. It performs two CLI runs and a direct engine run. Energy sum and direct
comparison tolerance: absolute 1e-12 Ha. Relative charge error: <1e-8. Independent
density integral tolerance: relative 1e-10. Selected density samples and hashes
are exact comparisons. Tests print measured errors; no tolerances are relaxed.
The synthetic fixture is not a physical Si pseudopotential or a QE reference.

```sh
scripts/benchmark.sh input.toml [output-directory]
```

This macOS helper wraps the release cargo invocation with `/usr/bin/time -l`,
records stderr, exit code, wall time and peak RSS in a unique `benchmark.*` text
report next to outputs. It explicitly overrides the TOML output directory with
the second argument (default `results`). It includes cargo startup/build overhead;
use a warm build for runtime measurements. Peak RSS is bytes on macOS. It is a
measurement helper, not a gate, and is not a portable Linux timing script.
In this restricted sandbox, `/usr/bin/time -l` fails with
`sysctl kern.clockrate: Operation not permitted` after running the command;
its report contains elapsed time but no peak RSS, and the helper returns 1.
Run it in a normal macOS terminal to obtain the missing RSS measurement.

The additive engine `run_outcome()` API returns `{converged, result}` on either
certified completion or exhausted budget. Existing `run()` retains its original
success-only result/error behavior and existing mixing checks. The final state
on budget exhaustion is returned before the unused final mixing proposal.
No existing numerical algorithms, tolerances or public signatures were changed.

## Explicit boundary units (P7U)

Retained numeric TOML fields accept the following units. Removed D4 keys and
their unit aliases stay invalid. The unit-aware `Config::parse` entry
point additionally accepts these alternatives (direct serde deserialization
continues to use the legacy structs):

| Section | Unit-string form | Legacy form it replaces |
| --- | --- | --- |
| `basis` | `cutoff = "30 Ry"` | `ecutwfc_ry = 30.0` |
| `system` | `lattice_vectors_with_units` (3×3 strings, one vector per row) | both `lattice_vectors` and `lattice_unit` |

For example, a general cell can use:

```toml
[system]
name = "example"
lattice_vectors_with_units = [
  ["5.43 Angstrom", "0 Bohr", "0 Bohr"],
  ["0 Bohr", "0.543 nm", "0 Bohr"],
  ["0 Bohr", "0 Bohr", "543 pm"],
]
# Supply atoms and the other required tables as before.
```

Both forms together are rejected, even if equal; the lattice string matrix
conflicts with either legacy lattice key. Cutoff and lattice remain required.
Omitted tolerances retain their existing defaults; finite-positive validation
is unchanged. Density residual tolerance keeps its existing Bohr^(-3/2) numeric
form: it is neither an energy nor a length. Fractional coordinates are unchanged.

The closed alias lists below are ASCII case-insensitive. No plurals or other
spellings are inferred. In particular `meV` (including `MeV`) always means
millielectronvolt; megaelectronvolt is unsupported.

| Dimension | Accepted aliases, grouped by unit |
| --- | --- |
| Energy | `ha`, `hartree`; `ry`, `rydberg`; `ev`, `electronvolt`; `mev`, `millielectronvolt`; `j`, `joule` |
| Length | `bohr`, `a0`; `angstrom`, `ang`, `a`; `nm`, `nanometer`; `pm`, `picometer`; `m`, `meter` |
| Pressure | `ha/bohr^3`; `ry/bohr^3`; `ev/angstrom^3`; `gpa`; `kbar`; `pa` |

Strings accept decimal/scientific numbers with or without whitespace before the
unit (`30 Ry`, `30Ry`, `408.2 eV`). Empty input, missing/unknown units, NaN,
infinity and extra tokens are errors that quote the offending input. Pressure
is available for explicit conversion. See P22A below for supported stress calculations.

```sh
cargo run --offline -p dft-cli -- convert 30 Ry Ha
# 1.50000000000000e1
cargo run --offline -p dft-cli -- convert 1 Ha/Bohr^3 GPa
```

`convert` prints one value in scientific notation with 15 significant digits.
Success is exit 0; usage/conversion errors (including incompatible dimensions
and nonfinite results) are exit 2. Unit-string TOML errors use the existing
configuration error status 1. Tested normal-range round trips have relative
error <1e-14; f64 underflow/overflow is outside that bound.

The additive results.json schema-version-1 key `total_energy` contains
`{ "ha": value, "ry": value, "ev": value }`, with all three derived from the
same final Hartree total through dft-units. `energies.total` and every other old
field remain. Energy conversion now uses uom's CODATA 2018 definitions; final
roundoff of eV output can differ from the historical rounded constant.

## K-point sampling (P7A2)

Under `[basis]`, legacy `kpoints = [[0,0,0]]` keeps its existing engine behavior
and exact JSON structure. New inputs use either:

```toml
kgrid = [2, 2, 2]
kshift = [0.5, 0.5, 0.5] # default [0,0,0]; each 0 or 0.5 in grid steps
ktime_reversal = true     # grid default
```

or:

```toml
kpoints = [[0.25, 0, 0, 0.5], [-0.25, 0, 0, 0.5]]
ktime_reversal = false    # explicit-list default
```

Grid dimensions are positive integers. The Gamma-centered convention is
`(m+s)/n` wrapped to `[-1/2,1/2)`; a half-shifted grid has no Gamma point.
Weights must be positive, finite and sum to one within 1e-14; no normalization
is performed. Explicit true applies time reversal to an explicit list. No space
group is consulted. The only allowed list/grid combination is an exact legacy
Gamma placeholder alongside a grid, which takes precedence. `kshift` requires a
grid. Invalid shapes, nonfinite coordinates and ambiguous combinations fail.

New forms require an empty band and the engine's sampled insulating gap
certificate (minimum 1e-7 Ha and greater than twice eigen tolerance). Failure
reports the engine cause with exit 4. Legacy Gamma keeps occupied-only support.

New forms add a `kpoints` JSON block: `grid`, `shift` (null for explicit lists),
`explicit`, `count`, `time_reversal`, and `points`. Each point contains
`fractional`, `weight`, `eigenvalues: {ha, ev}`, and `occupations`. The old
flat eigenvalue fields and basis plane-wave count describe the first k point.
Memory estimates cover one orbital block. Density files are unchanged.
Console tables show all k points up to 10, otherwise the first and last three;
results.json always contains all points. Direct and indirect HOMO/LUMO gaps
include k labels (one-based JSON order). Unconverged gaps are diagnostics.

See [ADR 0013](decisions/0013-kpoint-input.md) for the full schema, compatibility
rules and QE symmetry flags. `make_qe_input.py` emits automatic grid/shift cards
or explicit crystal lists, with `expected_kpoints` in its metadata and an
importable `compare_kpoints` helper for weighted-set verification. The existing
Gamma-only verifier remains limited to Gamma. See
[si_lda_kgrid.toml](../examples/si_lda_kgrid.toml) for an input template.

## Expression inputs (P7E2)

`Config::parse`, `validate`, `grid` and `run` accept expression strings in the
following explicit fields. Arrays may mix TOML numbers and strings. Plain
numbers retain their existing behavior; all existing physical constraints apply.

| TOML key path (indices are zero-based) | String resolution / value unit |
| --- | --- |
| `system.lattice_vectors[i][j]` | float expression, `lattice_unit` |
| `system.lattice_vectors_with_units[i][j]` | expression plus length unit, resolved Bohr; every entry requires a unit string |
| `system.atoms[i].fractional[j]` | float expression, dimensionless |
| `system.atoms[i].cartesian[j]` | float expression, `lattice_unit` (Bohr with unit-bearing lattice) |
| `basis.ecutwfc_ry` | float expression, Ry |
| `basis.cutoff` | expression plus energy unit, resolved Ry |
| `basis.kshift[i]` | float expression, grid steps; still restricted to 0 or 0.5 |
| `basis.kpoints[i][j]` | float expression, fractional coordinate or dimensionless weight |
| `basis.kgrid[i]`, `basis.fft_grid[i]` | integer expression, count |
| `scf.max_iterations`, `scf.nbands`, `scf.extra_bands`, `scf.seed` | integer expression, count/seed |
| `scf.mixing_beta` | float expression, dimensionless |

`cartesian` is an alternative to `fractional`; specifying both is an error.
Names, species labels, paths, enums and booleans are not evaluated. Their existing
validation remains in force; paths and descriptive names are literal strings.
Integer expressions use `eval_expr_integer` from ADR 0014 (integrality tolerance
`1e-12*abs(value)`, exact zero, i64 range), then existing consumer constraints;
`"8/2"` is accepted, `"7/2"` is rejected. Numeric TOML integer ranges are unchanged.

Expressions support `+ - * / ^`, parentheses, unary minus, constants `pi`, `e`,
and `sqrt abs sin cos tan asin acos atan exp ln log10`. Trigonometric arguments
are radians. Names such as `a` are unsupported: write `"5.43/2"`, not `"a/2"`.
There is no implicit multiplication. Powers associate left (`2^3^2 = 64`),
and unary minus binds tighter than power (`-2^2 = 4`). Nonfinite results fail.
For quantity expressions, separate the final unit token with whitespace:
`"0.5*30 Ry"`, `"5.43/2 Angstrom"`. Legacy `"30Ry"` still works.
See [ADR 0014](decisions/0014-expression-inputs.md) for the evaluator contract and
[si_lda_expr.toml](../examples/si_lda_expr.toml) for a complete template.

```sh
dft-cli eval 'atan(0.5*pi)'
dft-cli eval '1/4'                  # 2.50000000000000e-1
dft-cli convert '1/4*30' Ry Ha      # 3.75000000000000e0
```

`eval` requires no units and rejects unit tokens. `convert` requires source and
target units even for expressions. Both print 15 significant digits in scientific
notation, return 0 on success and 2 on expression/usage/conversion errors.
Configuration evaluation errors retain status 1 and include the key path and
quoted original string, for example `system.atoms[1].fractional[2]` and `"1/0"`.

`validate` echoes every resolved string as `key_path: "original text" = value`.
`results.json` optionally adds `input_expressions`, an array of
`{"key_path":"system.atoms[1].fractional[0]","text":"1/4","value":0.25}`.
All numeric-field strings are recorded, including legacy quantity strings and
plain numeric strings; numeric TOML values are not recorded. Units follow the
table above. Cartesian records precede conversion to fractional coordinates.
Integer record values are f64 representations of accepted integers; typed
integer controls retain their integer precision. The block is absent for
plain-number inputs; schema version and all legacy output fields stay unchanged.
No expression evaluation or provenance recording occurs inside numerical loops.

### PBE selection and explicit pseudopotential override (P7B)

The required `[xc] functional` accepts exactly `"lda-pz"` or `"pbe"`.
PBE accepts UPF labels `PBE`, `SLA PW PBX PBC` and the older `SLA PW PBE PBE` (case/whitespace normalized).
A mismatching label is rejected unless `[xc] allow_functional_override = true`.
The default override is false. The original UPF label remains in species
provenance; PBE and override-enabled results.json records `xc`, `allow_functional_override`, and
`functional_override_applied`. Console output records the same selection and
both override flags. Legacy LDA without override retains its existing JSON schema
(absent flags mean false). This permits explicit code-to-code tests using an LDA
pseudopotential with PBE; it does not establish physical transferability.
Spin polarization remains unsupported: `nspin=2` is rejected as an unknown
input field, as are all other spin controls. No functional is inferred from
an unknown UPF label and no missing functional gets a silent TOML default.

The QE input generator accepts PBE results with boolean
`allow_functional_override` and `functional_override_applied` provenance.
It checks these against the stored UPF label, emits `input_dft = 'pbe'`, and
records the label and both flags in QE input comments. It reuses the same
pseudopotential basename and SHA-256 without reading or modifying the UPF.
Missing override permission for a mismatched label is an error. Legacy LDA
QE input remains byte-identical to the stored golden file.

### Historical oracle SCF diagnostics (P7P4)

`results.json` convergence history now includes `operator_applications` (number
of Hamiltonian vector columns, not block calls), `davidson_iterations_per_k`
and `operator_applications_per_k`, in the same order as `kpoints.points`.
These counters cover the main solve only. Each `fixed_point_certificate` has
its own counters with the same names; sum both to obtain total work, including
failed certification attempts. Existing `davidson_iterations` is unchanged.

The recorded engine options additionally include `residual_eigen_schedule` and
`free_electron_initial_guess`. Both default to true. These are Rust `ScfOptions`
controls, not new TOML input keys. Setting `residual_eigen_schedule=false`
selects the legacy `initial_eigen_tolerance`/`eigen_tolerance_decay` schedule;
those controls retain their original meanings in that mode. Setting
`free_electron_initial_guess=false` restores the seeded random cold start.
The new schedule requests `max(0.1*eigen_tolerance,
min(1e-2, 0.1*previous_density_residual))` Ha, initially 1e-2 Ha. The density
residual is the existing normalized density L2 norm, numerically in atomic
units. All convergence criteria and the tighter independent certificate remain
unchanged. Cold initial vectors select the lowest kinetic energies, with basis
index breaking ties; Davidson performs its existing Rayleigh–Ritz procedure
and seeded rank completion. Warm starts retain previous eigenvectors.

### P7P6: deterministic k-point threads

`[scf] threads = N` selects standard-library scoped workers for independent
k-point Davidson solves. The default is 1; 0 also means 1. The key accepts a
nonnegative TOML integer; negative, fractional, string, boolean and overflowing
values are rejected. Actual workers are capped at the number of retained k
points; Gamma always uses the original serial loop. No automatic core detection
or nested parallelism is enabled: each worker uses sequential faer and serial FFT.

Only independent solves run concurrently, including the fixed-point certificate
solves. Density, energy, weights, mixing and certificate reductions retain their
original sequential order. Numerical results are bit-identical across thread
counts; elapsed times are naturally different. `diagonalization_s` measures wall
time around the complete solve region, including worker startup and joins,
rather than summed worker CPU time.

`results.json.threading` retains `threads` (now the actual worker count) and
`backend`, and adds `threads_requested` (raw input, default 1) and `threads_used`
(`min(max(threads_requested, 1), number_of_k_points)`). `scf.threads` also records
the requested value. The backend description identifies scoped k-point workers
and sequential numerical backends.

Every scoped worker is joined before a solve error is returned. The lowest
zero-based failing k index wins, independent of completion order. Unwinding solve
panics become indexed numerical errors; the standard panic hook may still print
a diagnostic. Abort-mode panics and allocation/process failures cannot be caught.
The serial path also converts Davidson panics into indexed errors.

Benchmark inputs can opt in with `scripts/bench_scf.py --threads N`; omitting
that option leaves the generated input unchanged.

### P7PB: total thread budget within and across k points

`[scf] threads = N` remains a total budget (default 1; zero normalizes to 1).
For `nk` points, `outer = min(max(N, 1), nk)` and
`inner = max(1, max(N, 1) / outer)` using integer division. Thus
`outer * inner <= max(N, 1)`; remainders are unused. For example, N=8 and
nk=3 gives outer=3, inner=2; Gamma with N=4 gives outer=1, inner=4.

The `results.json` `threading` object adds `outer_threads` and `inner_threads`.
The existing `threads` and `threads_used` fields retain their historical
meaning (outer k workers), and `threads_requested` retains the input value.
Inner threads are a budget: a smaller block may use fewer workers. Timings
include thread creation and joins. The backend description now mentions bands.

Within a k point, Hamiltonian block vectors and density band transforms use
scoped standard-library threads with persistent numerical scratch. Density
contributions are reduced in ascending band order, then ascending k order.
Each vector uses the existing FFT arithmetic. Single-thread runs spawn no inner
threads; small matrix operations and individual FFT transforms remain serial.
No numerical tolerance, occupation, default band count or SCF algorithm changes.

`bench_scf.py --cases si8` adds the 8-atom conventional cubic Si cell at
5.43 Angstrom and 30 Ry, with default bands and the supplied read-only UPF.
The starting binary fails to converge its Davidson solve for this input;
therefore this case currently reports failure, not an accepted SCF speedup.
Benchmark records add `phase_medians_s` and `engine_total_median_s`; the older
last-run fields and default case selection retain their meanings.

## Pseudopotential registry (P9B)

Offline commands use the embedded metadata-only registry:

```sh
dft-cli pseudo list [--element Si] [--supported]
dft-cli pseudo show ID
dft-cli pseudo verify /outside/file.upf [--id ID]
dft-cli pseudo path
dft-cli pseudo suggest SiO2
dft-cli pseudo suggest Si O
dft-cli pseudo fetch ID [--yes] [--offline]
baandi inspect-upf /outside/file.upf --json
```

`pseudo --registry /outside/registry.json COMMAND` selects an explicitly trusted
registry; validate it first with `python3 -I scripts/pseudo_registry_check.py`.
Unknown/mismatching identities and fetch failures return nonzero. JSON inspection
returns exit 0 for a completed report, with `supported_by_engine: false` and
reasons for unsupported/invalid contents; I/O failures return nonzero. Existing
text inspection behavior is unchanged. Suggestions require verified identity and
engine support, rank ascending source wavefunction cutoff (Ry), unknown last, then
ID; they print XC and reasons and do not select a potential for a calculation.

Fetch accepts only verified, supported entries with verified licenses. It rehashes
cache hits and checks provenance. Offline misses fail without network. Online
misses require `--yes` and then fail with the deliberately unavailable P9A transport;
the hash-verifying publication path is implemented and tested with synthetic data.
`BAANDI_PSEUDO_DIR` overrides the platform cache; it must be absolute and fetch
rejects repositories. The initial registry has no verified files. Input `registry:`
resolution is a future P8G contract, not implemented config syntax.

See [the user guide](howto-pseudopotentials.md) for platform cache paths and
[source evidence, schema and fetch contracts](architecture/pseudopotential-registry.md).
## Volumetric output (P8B)

The optional `[output.volumetric]` table accepts required nonempty `quantities`
(from `density`, `elf`, `hartree`, `local_potential`, `partial_charge`) and `formats`
(from `cube`, `xsf`, `vasp`, `vti`). Duplicates and unknown keys/values are rejected. Optional
`kpoint` (default 1) and inclusive `bands = [first,last]` (default `[1,1]`) are
one-based and apply only to `partial_charge`. Partial charge is unweighted,
with one unit of integrated probability per selected band. Files are named
`<quantity>.<format>` in the normal output directory and cannot overwrite
existing files/symlinks. Raw density output is unchanged.

Normal `run` emits the requested files after the final SCF state is returned.
The `vasp` format uses `<quantity>.CHGCAR` for density and partial charge,
`<quantity>.LOCPOT` for the potentials, and `elf.ELFCAR` for ELF; other suffixes
match the format name. VASP density is multiplied by cell volume, potentials are
converted to eV, and ELF is dimensionless. VTI contains an Angstrom affine grid
and explicitly named atomic-unit scalar arrays. Export failures use exit 5.
The existing dft-physics dependency is now a runtime dependency (approved
manifest exception); no new package was added.
See [the viewer guide](howto-volumetric.md) for units, pseudopotential caveats,
periodic endpoints and the independent Python validator, and
[ADR 0021](decisions/0021-volumetric-output.md) for numerical and format contracts.
## Fixed-density bands and DOS (P8A)

```sh
baandi bands input.toml [--output-dir DIR]
baandi dos input.toml [--output-dir DIR]
python3 -I scripts/plot_bands_dos.py DIR/bands.csv DIR/dos.csv -o DIR/figure.svg
```

Both commands run verified reference SCF then reference fixed-density NSCF
by default (strict requests use the oracle). `[scf] threads`
controls both stages. No change is made to `run` output behavior. The new
commands produce their spectral files, not `run`'s density/restart outputs.
`bands` requires Python 3 for its ready SVG; the plot script has no third-party
dependencies. See the complete, non-programmer [guide and measured Si example](howto-bands-dos.md).

Optional `[bands]` keys: exactly one of `preset = "fcc" | "bcc" | "sc"` or
`points = [{label="Gamma", fractional=[0,0,0]}, ...]`; points may set
`break_before=true`. Sampling is `points_per_segment` (default 20 **intervals**)
or `spacing_inv_bohr`, mutually exclusive. Primitive preset geometry must
match the documented ordered convention up to a proper rigid rotation.
`nbands`, `eigen_tolerance_ha` (default 1e-8), and `max_iterations` (default 400)
control the NSCF Davidson solve. At least occupied + 4 bands are required.
All bands use the strict residual tolerance; degenerate top multiplets can
need more bands/iterations. Nonconvergence identifies the k index.

Optional `[dos]` keys: `dos_kgrid` (default twice SCF grid dimensions, or 2×2×2
for legacy Gamma), `sigma_ev` (default 0.1 standard deviation), `spacing_ev`
(default maximum spacing 0.01), `window_ev = [min,max]` relative to the sampled
VBM (default spectrum extent plus 8 sigma tails), and the same `nbands`,
`eigen_tolerance_ha`, `max_iterations` controls. Explicit non-Gamma SCF point
lists require an explicit DOS grid. Time reversal is always used for DOS;
weights sum to one and spin multiplicity two is included. New table numeric
fields currently accept literal TOML numbers, not expression strings.

Both commands also write routing metadata to `results.json`.
Band outputs: `bands.csv`, `bands_labels.csv`, `bands.json`, `bands.dat`,
`bands.dat.gnu`, `bands.agr`, `EIGENVAL`, `bands.svg`. DOS outputs: `dos.csv`,
`dos.json`, `dos.dat`, `dos.agr`, `DOSCAR`, `dos.svg`. Units, column layouts,
public format references, and format-specific limitations are in the guide.
The VBM is the maximum occupied eigenvalue over the computed sampling; JSON
stores its absolute Hartree offset. No Fermi energy is computed. DOSCAR's
nominal Fermi slot is a documented VBM-zero placeholder. These are rust-dft
exports, not calculations by the other codes. Output input-alias/symlink/hardlink
checks occur before computing. Exit codes are 1 input, 2 UPF, 3 SCF exhaustion,
4 numerical/noninsulating/NSCF failure, 5 I/O/plot failure. Writes are not a
multi-file transaction; check exit status before using files.
## Full input preflight (P8G)

`baandi validate input.toml [--json]` now opens and validates the referenced
UPF files, resolved relative to the input file. This supersedes earlier descriptions
of `validate` as a geometry-only check. `grid` remains a geometry/grid inspection.
`run` invokes the identical preflight before constructing an SCF problem or writing
outputs. No density/checkpoint data is read and no SCF is performed by preflight;
UPF radial arrays are parsed to establish that the potential is supported.

Existing string entries remain supported. Custom species labels require an explicit
chemical symbol; symbols are case-sensitive:

```toml
[pseudopotentials]
Si = "Si.UPF"
Si_up = { path = "Si.UPF", element = "Si" }
# Deliberate identity-mismatch experiment only:
Fe = { path = "Si.UPF", allow_element_mismatch = true }
```

The declared element (or chemical-symbol label) must match the UPF header element.
An override does not change the potential or electron count. Overrides are printed
and recorded, including permission enabled without a mismatch. Custom tags still
need a valid explicit element even with the override. Functional compatibility uses
the engine's existing label matcher and `xc.allow_functional_override` permission.

Checks include cell/coordinates, k-point quadrature, supported UPF data, finite
positive tolerances, positive iteration and thread counts, mixer name and beta in
(0,1], positive integral even valence-electron count, occupied bands, and the
basis/subspace capacity at every resolved k point. A successful preflight does not
guarantee SCF convergence, a gap, numerical stability, or adequate physical accuracy.

Cutoff rule: warn when a finite positive UPF `wfc_cutoff` (Ry) exceeds
`basis.ecutwfc_ry`; zero/nonpositive/nonfinite metadata gives no recommendation.
This is a minimal wavefunction-cutoff advisory, not a density-cutoff or convergence
certification. The FFT grid and plane-wave counts are enumerated with the same
basis APIs as the run. Memory is an array-payload estimate for orbitals at all k
points plus one real and one complex grid, in bytes; scratch, projectors, basis
storage and allocator overhead are excluded. It is not a peak-memory bound.

Text output includes the UPF composition and labels, atoms, electrons, XC, requested
k mesh/shift (or Gamma/explicit list) and resolved point count, cutoff in Ry, bands
and threads. Overrides identify the actual UPF header material.
Validation errors retain exit 1 (command-line syntax errors use clap exit 2);
warnings do not change exit 0. JSON mode writes one JSON object to stdout on
success; errors retain the ordinary stderr diagnostic and nonzero exit.

The version-1 preflight JSON object has `schema_version`, `valid`, string arrays
`warnings` and `overrides`, `summary`, integer `electrons` and `bands`, integer
array `plane_waves_per_k`, three-integer `fft_grid`, integer
`estimated_array_bytes`, and `memory_note`. Successful runs add this object under
`results.json.preflight`, including runs that exhaust their iteration budget.

Compatibility note: example inputs containing `/absolute/path/to/Si.UPF` require
substitution of a local supported Si UPF. Tests use generated synthetic Si headers;
the numerical fixture arrays are unchanged. Missing/unsupported UPFs and functional
mismatches now fail run's preflight with exit 1, rather than the previous late-stage
exit 5/2. Zero threads and mixing beta above one are rejected before solver setup.
Output I/O failures still use exit 5; numerical failures still use exit 4.
### Density mixing

Production accepts `scf.mixing = "plain" | "broyden" | "pulay"` (the last is
an alias), `mixing_beta`, `mixing_history`, and `reference_preconditioning =
"plain" | "tf"`. The legacy `linear` and `kerker-pulay` values and the
`mixing_q0_inverse_bohr` and `mixing_metric` keys are no longer accepted.

### Oracle cutoff ladder compatibility

`scf.cutoff_ladder = "off"` disables the strict oracle's historical coarse
initialization. Omission retains that frozen oracle default. The production
reference path accepts `off` as a no-op. `on` is an unknown value and
`cutoff_ladder_fraction` is an unknown key.

### Oracle band accuracy fields (P7PC)

`results.json.band_accuracy` describes the legacy first-k eigenvalue list.
For explicit k-point output, each `kpoints.points[]` also has `band_accuracy`.
Each entry contains `residual_ha` (residual vector 2-norm), `buffer` and
`converged` (whether that band's requested residual target was met). Existing
energy, eigenvalue and convergence fields retain their names and units.

Occupied bands and the lowest empty band use the requested strict tolerance.
Higher empty bands use `max(requested_tolerance, 1e-5)` Ha and are printed with
`(buffer; less accurate)`. Their printed decimal places do not imply strict
accuracy. At the solver iteration limit buffers can be returned with
`converged=false`; all strict bands must still pass. A residual bounds distance
to the spectrum, not the error in an individually indexed member of a cluster.

Fixed-point history certificates add `strict_residuals_per_k_ha` and
`buffer_residuals_per_k_ha`. Only the former belong to the eigenpair certificate;
the latter are separate informational measurements. The sampled-gap test is
unchanged and uses the strict lowest empty band. Eight bands for four occupied
states means five strict bands and three buffers, not four buffers. Default
`extra_bands=4` is unchanged. The engine uses at least six subspace columns per
requested band; the configured `davidson_subspace` is a minimum working capacity.

Intermediate SCF history eigen residuals may be recurrence residuals while the
requested tolerance is above the final `eigen_tolerance`. Final-result band
residuals and all fixed-point certificate arrays remain true operator residuals;
final-accuracy maps always perform that check. This removes duplicate operator
work without weakening the final acceptance tests.

## Optional agent preflight findings (P11K, partial implementation)

`baandi validate INPUT --json` adds `findings_schema_version: 1`, typed `findings`
and a read-only `structure` summary to successful reports. Codes currently cover
cutoff suggestions and explicit element/functional overrides; each carries a
severity, TOML key path, typed evidence and fix hint. Structure vectors are rows
in Bohr and atom positions are fractional. Legacy text warnings, verdicts, exit
codes and `results.json.preflight` persistence are preserved; no agent is required
and the CLI makes no network call. Failed validation still returns stderr text
and a nonzero exit, not a findings JSON envelope.

Use the existing `explain INPUT` for human explanation. `baandi check --explain`
is not currently a supported command. The external advisory/patch workflow,
restricted-field approval rules, proposed additional structure checks and additive
provenance contract are specified in [agent preflight](architecture/agent-preflight.md).
Proposal validation, applying approved diffs and writing decision sidecars are
not implemented. Agent advice cannot change a deterministic verdict or modify
input automatically.

## Completed optional preflight review (P11K)

This section supersedes the earlier partial-implementation description above.

`dft-cli check INPUT --json` is an alias for `validate` and emits deterministic
findings on success and failure (exit 0 valid, 1 invalid). `check INPUT --explain`
prints the findings and authority explanation without an agent. Successful reports
retain legacy fields and add `findings_schema_version: 1`, typed findings and a
read-only structure summary. Failure reports identify stopped coverage rather
than inventing unavailable metrics. Runs use the same deterministic rules.

Optional `[preflight]` thresholds enable warnings for minimum periodic distance,
cell volume/skew, explicit slab vacuum and reciprocal k-mesh spacing. Overlapping
atoms and invalid cells fail. Cutoff suggestions, electron parity and species
identity remain deterministic. All units, codes and policy fields are documented
in [agent preflight](architecture/agent-preflight.md).

`dft-cli preflight-review INPUT PROPOSAL.json` validates and displays a proposed
diff without changing input or recording approval. Add `--decision approved`
(or `rejected`), `--actor NAME`, `--reason TEXT`, and `--output-dir RESULTS` to
record an explicit decision beside `results.json`. Approval requires a valid
candidate and fresh source, UPF and findings hashes. The immutable review sidecar
contains proposal, decision, candidate and recheck; existing results are untouched.
Even approved proposals are never applied by this CLI. There are no network calls,
UPF downloads or mandatory agent dependencies. See the linked contract for the
strict full-file unified-diff subset and versioned proposal fields.

## Collinear fixed integer spin (P13B / MF03)

See [the collinear guide](user/collinear-spin.md) for `[spin]` keys, full-mesh
requirements, fixed-sector interpretation, certificates, reproducible evidence
and limits. Supported polarized XC choices are `pbe` and `lda-pw`; the latter
also supports unpolarized runs. `lda-pz` remains unpolarized only. A collinear
`run` writes the version 2 spin result profile and channel/total/magnetization
cubes; omitted spin settings preserve version 1 output. `validate` and `run`
share sector/parity/capacity, axis, unsupported-functional and symmetry checks.

### P13C NC datasets and NLCC

`baandi inspect-upf FILE --json [--manifest docs/data/nc-datasets.json]`
reports SHA-256, ordered one-based projector channels, wavefunction/mesh counts,
original XC and relativistic labels, NLCC, and cutoff hints in Ry (absent hints
are null). The optional manifest is matched by hash, never filename; missing or
ambiguous identity and element mismatches fail. License provenance comes from
that manifest. Text mode prints facts and blockers before numerical validation.
JSON inspection exits zero when the report can be produced even for unsupported
content: callers must examine `supported_by_engine` and `blockers`. Inspection
and a successful parse are not material qualification.

The five PseudoDojo files now parse without modifying the vendor or dependencies.
The exact third-party failure was `missing field @wfc_cutoff`; the old adapter
also rejected `PP_BETA.4`. Compatibility normalization in dft-upf validates original
numbered suffixes, indices, order, channel count and D before use. Extra numbered
nodes map to repeatable third-party enum variants without dropping or reordering
any array, index, occupation or D entry. Missing informational header/mesh
attributes receive neutral reader-only defaults; the mesh count derives from the
header. Only the explicit PP_R and PP_RAB arrays enter radial numerics. Original
missing metadata remains absent in inspection; `raw()` is explicitly a normalized
compatibility view, not an exact representation of the source XML. A finite
projector sample at r=0 is preserved: its r*u*j_l quadrature weight is exactly
zero, with a bit-exact synthetic test rather than an acceptance tolerance.

PP_NLCC is rho_core(r), electrons/Bohr^3, not 4*pi*r^2*rho. Its spherical Fourier
transform uses the existing Simpson/3/8 index quadrature, shell caching at setup,
1/volume normalization and atomic phases on the density sphere/grid. The prepared
core grid is immutable and never enters electron count, Hartree or orbital
density. LDA-PZ, LDA-PW and PBE evaluate XC on valence+core; collinear LDA-PW/PBE
split the core equally between channels, including GGA gradients. Band double
counting uses valence alone. `energies.nlcc`, present only for NLCC, exposes core
charge, core-potential integral (audit only) and valence-potential integral.
The XC energy already includes core; never add the core audit integral to total E.
Local-potential volumetric output follows the same XC convention.

Negative combined-density rules are unchanged: LDA-PZ and unpolarized LDA-PW
allow their existing <=1e-14 Bohr^-3 roundoff clamp; unpolarized PBE and both spin
functionals reject any negative combined density. No separate clipping or
renormalization of core samples occurs. Fourier truncation can make a coarse
core grid infeasible, and the error names the combined NLCC density. The density
threshold remains 1e-10 Bohr^-3. Derivative oracles avoid threshold crossings.
At that stage forces and stress were unavailable; P20A/P22A below add strict certified derivatives. NLCC checkpoint v1 export/import rejects
explicitly because that format cannot represent the new ledger; ordinary v1
encoding/decoding is unchanged. A spin model identity records active NLCC.

Reproduction (all real files remain external):

```sh
export DFT_NC_DIR=/absolute/external/pseudodojo/ex
python3 -I scripts/check_upf_headers.py --binary target/debug/baandi \
  --output tests/reference/datasets/pseudodojo-headers.json
cargo test --release -p dft-cli --test all inspect_upf::real_nc_headers -- --ignored --nocapture
cargo test --release -p dft-engine --test all nlcc:: -- --include-ignored --nocapture
cargo test --release -p dft-engine --test all real_nc:: -- --ignored --nocapture --test-threads=1
# Additional discretization diagnostic; does not change engine defaults:
DFT_NC_CUTOFF=30 cargo test --release -p dft-engine --test all real_nc::isolated_atoms_pbe_smoke -- --ignored --nocapture
```

Frozen synthetic caps: radial/atomic integrals 1e-10 relative, XC directional
finite difference 1e-7 Ha per dimensionless perturbation with three h steps,
translation energy 1e-7 Ha. The translation cap is this ticket's explicit
pre-measurement choice at the MF06 derivative energy scale; the source table
specified no separate translation-energy cap. Cubic, skew and elongated cells,
even/odd radial meshes, spin swap and zero-spin reduction are covered. Five
coupled synthetic projectors are checked against analytic Gaussian radial
transforms and the spherical-harmonic addition theorem. Real-file D assembly
is independent of the production matrix assembly, but reuses its radial
transforms; it is not an independent scattering/ghost oracle.

Evidence and commands are in `tests/reference/datasets/p13c-status.json` and
`pseudodojo-headers.json`. Real atom runs use a 10 Bohr cubic cell, Gamma, fixed
spin sectors, 20/30 Ha and at most 60 iterations. Failures and unconverged states
are retained as smoke diagnostics, never labeled successful SCFs. Real finite
atomic-density tables have small charge deficits; they are reported without
renormalization and are not certified against the analytic 1e-10 normalization
cap. All five datasets remain unqualified: cutoff scans, ghost/scattering tests
and matched external-code validation are absent. No QE agreement is claimed.

README regeneration and its specification tests retain the existing evidence.
Adding a new README result renderer and dashboard/session publication remains
outside this ticket's explicit owned-file list.

The two-atom skew-cell total-SCF translation test additionally holds the
wavefunction cutoff at 6 Ha and refines only the FFT grid. PBE fails at
18x15x20 (3.0781221274e-6 Ha), passes at 36x30x40
(6.7368192358e-8 Ha), and improves to 5.1714472704e-9 Ha at 54x45x60.
The cap stays 1e-7 Ha throughout; coarse-grid failure is retained in evidence.
This is a resolved-grid invariance check, not a promise that default grids give
that tolerance for every system. Production grid/cutoff defaults are unchanged.

## P13D bulk NC XC floor

`xc.density_floor` is no longer an input key. The strict oracle retains its
historical NLCC-only floor internally; production uses its reference XC policy.
Historical floor measurements remain in [bulk smoke evidence](validation/materials-smoke.md).

### Standard convergence profile (ADR 0035)

The input key `scf.convergence_profile` was **removed** on 2026-11-16 (user
decision, ADR 0043): an input that still contains it (any value) is an input
error that names the replacement. Subject to the explicit strict/G4D routes, every command (`run` in every execution mode,
`relax`, `bands`, `dos`, `checkpoint`, `restart`, `explain`) uses the standard
profile with the fixed-point verification (next section; single path complete
since ONEPATH-4). `standard_threshold_ha`
(positive, finite, default **5e-9 Ha**, frozen P17A E1) sets its accuracy;
`verify = "off"` disables the reference verification; it does not restore
the removed legacy standard implementation.
The strict path is a reference/test tool only: `BAANDI_REFERENCE_STRICT=on`
(environment) runs the same input on the strict path with its certificate;
strict arithmetic, certificate, output schema version and omitted-extension
fields are unchanged there. Collinear spin and Fermi-Dirac inputs run the
standard profile too since ONEPATH-3 (see below).

Reference stops on the total reciprocal-density residual metric (without
per-electron division), represented as dr2 in Ry. The eigenvalue-change
threshold starts at 1e-2 Ry = 5e-3 Ha and ratchets as
`max(1e-13, min(previous_ethr, 0.1*dr2/max(1,N_e)))` in Ry. On the first map
all bands use that threshold; subsequently empty bands use
`max(5*threshold_ha, 5e-6 Ha)`. This is Ritz-change stopping, not the oracle's
occupied/first-empty residual certificate. Retry and verification work are
charged separately; see [driver/solver contract](design/reference-scf.md).

The independent fixed-point verification is standard-class evidence, not a
strict certificate or a proof of spectral completeness. `scf.spectral_guard`
is removed. Historical P17A profile data below retain their original driver.

`python3 -I scripts/equal_accuracy.py --output-dir /fresh/path` runs strict,
two tight references, three standard thresholds and strict/Hartree E4.
Exit 3 means incomplete evidence, including expected unavailable references;
work is still reported when available. Use `--qe-output /path/pw.out` (repeatable),
`--qe-reference /path/tight.out`, and optionally `--qe-density` plus
`--qe-reference-density` for explicit QE files. Set `--cases ''` for QE-only.
Raw QE parsing was tested on all four existing k222/k222s/si30 outputs.
The validation-claim scripts reject non-strict or explicitly uncertified results.
See [measured accuracy, limitations and commands](performance/p17a-standard-profile.md).

### Single SCF path and fixed-point verification (ADR 0043, default since 2026-11-16)

Production Reference uses the standard profile with `standard_threshold_ha` from the
input (default 5e-9 Ha) and reference Broyden beta 0.7/history 8, followed by a **fixed-point
verification** after the stop rule fires: one more ordinary map at the
delivered output density (warm start, eigen tolerance = the next ratchet value)
passes when the Coulomb fixed-point residual `E_H[F(rho_out) - rho_out]` is
below the threshold (and, with `BAANDI_STANDARD_RESIDUAL_GUARD=on`, the default,
the residual estimate of that solve; this guard is an empirical estimator, not
an error bound, see ADR 0043). An adopted predictor can defer an attempt;
verified success always requires an actual independent map. A failed verification continues SCF and verifies
again at the next stop; the delivered state is always the converged step's own
(verification never changes a number).

`[scf] verify = "fixed_point" | "off"`: default `fixed_point`; `off` skips the
verification; it does not select a different driver.

`scf.driver = "reference"` is accepted as a no-op. Production always uses
Reference unless `BAANDI_REFERENCE_STRICT=on` requests StrictOracle, or a
free-moment collinear Fermi-Dirac PBE input requires the recorded
`G4D_FALLBACK` strict oracle. The latter is still unqualified on the reference
path. `scf.driver = "legacy"` is an unknown value. `REFERENCE_SCF` is no longer
a registered switch. There is no legacy standard driver and no retry on engine errors.

Legacy-only keys are ordinary unknown-key errors, including `xc.density_floor`,
`scf.mixing_q0_inverse_bohr`, `mixing_metric`, `hartree_mixing`, `spectral_guard`, `locking`,
`preconditioner`, `davidson_capacity_per_band`, `energy_tolerance_ha`,
`density_residual_tolerance`, `eigen_tolerance`, `standard_stall_guard`,
`standard_density_forcing`, and `cutoff_ladder_fraction`. The energy-unit
aliases of the removed tolerances are also unknown. `start_density` and
`start_wavefunctions` accept only `"atomic"`; omitted oracle settings preserve
its historical defaults. These errors also apply with the strict switch on.

Every selection prints one JSON-bearing stderr line. Result metadata includes
`scf.driver` (`reference`, `strict_oracle`), `scf.driver_reason` (nullable)
and `scf.driver_source` (`default`, `env`, `input`). Source describes the winning
request even when a feature requires fallback. `reference_driver_fallback_reason`
is retained. Bands/DOS also write `results.json`; checkpoint/restart write a
`<checkpoint-output>.results.json` companion without changing checkpoint v1.
The stdout-only `explain` and `refscf-map` commands include the same `scf` object
in their JSON output. The diagnostic's `--solver` is an explicit input request;
STRICT takes precedence and retains residual certification.

`BAANDI_REFERENCE_STRICT=on|off` (default `off`) forces the strict reference path
with its fixed-point certificate for accuracy gates and tests; the standard-only
keys of the input (`standard_threshold_ha`, `verify`) are then ignored, which the result records
(`single_path.status = "reference_strict"`).

Every result carries `single_path` (`status`: `applied` or `reference_strict`
with `reason`; `not_applicable` is no longer produced since ONEPATH-4;
`effective_profile`; `verification`:
`fixed_point`, `off` or `strict_certificate`) and `execution` (review F01:
`requested` controls, `effective_profile`, `standard_threshold_ha`,
`verification`, `effective_solver` (`real_gamma_lean_davidson`,
`real_gamma_classic_davidson`, `complex_davidson` or an MPI distributed solver),
`real_gamma` and its fallback reason, `davidson_mode`, `lean_solves`, `backend`
(`linear_algebra` mkl or pure, `mpi_feature`, `ranks`, `threads`),
`parallel_model` and `non_default_switches`; the complete switch table stays in
`scf.switches`). Verified runs carry `verified`, `verification` (every attempt:
`coulomb_error_ha`, `threshold_ha`, `eigen_tolerance_ha`, `residual_bound_ha`
(the residual estimate; name kept for compatibility), `density_residual`,
`occupied_eigenvalue_change_ha`, information-only `strict_criteria_met`, work and
`wall_s`) and `timings.verification_s`; verification work is also in each history
entry's `fixed_point_certificate` and in `timings.certificate`.

Forces accept a standard state with a passed verification; `forces.preconditions`
records `standard_threshold_ha`, `fixed_point_verification_passed`,
`verification_attempts` and the last attempt's `verification_coulomb_error_ha`
and `verification_residual_bound_ha` (ONEPATH-1). The forces are computed with
the same formulas as on the strict path, from the converged step's own state.
`output.forces` with `verify = "off"` is an input error.

Fixed-cell relaxation (`baandi relax`, `relax.mode = "ions"`, the default) runs
every ionic step on the standard path with the verification (ONEPATH-1, user
decision 2026-11-16): a failed verification continues that step's SCF, warm
starts and `start_extrapolation` are unchanged, and the relaxation
`results.json` carries `single_path`, `scf.standard_threshold_ha` and a
`verification` record per frame (`passed`, `threshold_ha`, `attempts`,
`coulomb_error_ha`, `residual_bound_ha`, `strict_criteria_met`, `h_vectors`;
no wall time in this payload). Reference optimizer resume lacks an orbital
seed and is not guaranteed bit-identical; see [relaxation](user/relaxation.md). `verify = "off"`
with `[relax]` is an input error. The threshold is `scf.standard_threshold_ha`
(default 5e-9 Ha). **For force-sensitive runs (tight `force_tol_ha_bohr`, small
`energy_tol_ha`, phonon-like finite differences) tighten the threshold**: the
SCF energy is resolved to about the threshold, so an ionic energy test at
1e-12 Ha needs a threshold of the same order (the synthetic oracle in
`tests/relaxation.rs` uses 5e-13 Ha), and the measured force differences
against the strict reference shrink with the threshold (ledger ONEPATH-1;
dev08 request `performance/requests/onepath1-forces-relax.md`).

Stress (`output.stress`) and the variable-cell relaxation modes `cell` and
`ions_cell` run on the standard path with the verification since ONEPATH-2
(user decision 2026-11-16; the reasons `stress` and `stress (variable-cell
relaxation)` are gone). The stress needs a passed verification exactly like the
forces; `stress.preconditions` carries the same standard keys as
`forces.preconditions`; the three cutoff-scan SCFs of the isotropic Pulay
estimate run the standard path with the verification at
`min(1e-9 Ha, scf.standard_threshold_ha)` (repeated once at the caller's
threshold when that is looser; no strict fallback), recorded in
`stress.pulay_assessment.scan_profile` (`standard_verified`),
`scan_thresholds_ha` and `scan_energy_resolution_ha`. Every variable-cell step
(including basis rebuilds and the final nominal-cutoff point) verifies and
records `frames[].verification`; the frozen Miller chart and the Pulay handling
are unchanged. `verify = "off"` with `output.stress` or with any `[relax]` is an
input error. **Stress default threshold (provisional):** an input with
`output.stress = true` or `relax.mode = "cell"`/`"ions_cell"` and no
`scf.standard_threshold_ha` runs at 5e-11 Ha instead of 5e-9 Ha, because 5e-9
missed the 1e-7 Ha/Bohr^3 stress class on a synthetic nonlocal fixture
(1.3e-7) while 5e-11 met it on both fixtures (margin 1.6x). Forces,
`relax.mode = "ions"` and plain runs keep 5e-9 Ha; an explicit
`scf.standard_threshold_ha` always wins. `single_path.standard_threshold_ha` and
`single_path.threshold_source` (`user`, `stress_default`, `default`) record
what applied (also in the relaxation `scf` block). The rule awaits
confirmation by the dev08 request `performance/requests/onepath2-stress.md`;
measurements: ledger rows ONEPATH-2 and ONEPATH-2-STRESS-DEFAULT.

Except for G4D free-moment PBE, collinear spin and Fermi-Dirac occupations run
on the reference standard path with the verification since ONEPATH-3 (user decision
2026-11-16; the reasons `collinear spin` and `Fermi-Dirac occupations` are gone;
design `design/onepath3-spin-metals.md`). The stop estimate of a collinear run is
`E_H[d_n] + E_H[d_m] + E_0[d_m]` (charge residual, magnetization residual in the
same Coulomb metric, and a uniform-moment term weighted like the longest
wavelength of the cell; an estimator, not a bound); the residual guard sums
`f |r|^2` over both spin channels with the actual (fractional under
Fermi-Dirac) occupations. The verification repeats the whole map including the
Fermi-level solve and passes on the same estimate below the threshold plus the
number, moment and band-edge checks of the strict certificate; Fermi level,
moment and occupation changes are recorded as information
(`verification.attempts[].chemical_potential_change_ha`, `moment_change_mu_b`,
`occupation_change_electrons`, `spin_density_residual`). The delivered state
(energies, F, mu, moment, occupations) is the stop step's own. Spin and FD
results carry `convergence_profile`, `standard_threshold_ha`,
`convergence.standard_error_ha` (per step), `verified`, `verification`,
`timings.verification_s`, and in each spin history entry `standard_error_ha` and
`requested_eigen_tolerance_ha`; the strict `convergence.certificate` is null on
the standard path. These runs use the complex Davidson
(`execution.effective_solver = "complex_davidson"`; real Gamma excludes spin and
FD as before). Forces and stress of a verified spin or FD state are available
with the standard preconditions. CLI spin/FD relaxation remains unsupported. Reference spin/FD observation and
safe-point control use the same electronic stepper; oracle spin/FD (including
G4D) and oracle symmetry runtime combinations still reject.

The strict oracle is selected by G4D fallback or through `BAANDI_REFERENCE_STRICT=on`.
Python uses `Options.reference_strict = True` with `Options.driver = "reference"`;
the old Python `"legacy"` value is also removed.

- `run --observe` and `run --agent` drive `refscf::runtime::Runtime`, with one
  report per ordinary/retry/verification map. Iteration numbers can repeat;
  interventions require a safe point. Outputs include `events.jsonl`,
  `runtime.json` and the normal results; controlled runs also write the audit.
  Policies are beta, history, tightening `threshold_ry` and verification.
  Legacy eigensolver/tolerance/schedule/budget controls reject. See
  [execution modes](user/cli-modes.md) for transport and approval boundaries.
- `bands` and `dos` use reference SCF and `refscf::nscf` at fixed delivered density.
  The NSCF stage uses the requested residual tolerance and reports its driver
  separately. Strict requests retain the oracle diagnostic path.
- `checkpoint` uses v2 on Reference and v1 on StrictOracle; restart dispatches by
  checked version and compatibility. The exact v2 continuation retains kept HX,
  mixer and ratchet state, reproducing uninterrupted same-rank numerics without
  the old v1 REUSE_PSI caveat. See [checkpoint v2](checkpoint-format.md).

Reference k/spin/FD numerical work is distributed under MPI; historical F02
`mpi_replicated` descriptions apply to oracle solves that actually replicate.
Always compare equal parallel models, core counts and conditions.

## Finite-temperature occupations (P17B)

Omitting `[occupations]` preserves fixed insulating occupations and all existing
rejections. Explicit fixed mode is `[occupations] mode = "fixed"` (written as a
normal TOML table with `mode` on the next line). Fixed mode accepts no smearing.

```toml
[occupations]
mode = "fermi_dirac"
smearing_ha = 0.01
band_cap = 24
```

Supply exactly one of `smearing_ha`, `smearing_ev`, `smearing_kelvin`, finite and
positive. Kelvin uses the exact SI Boltzmann constant at the CLI boundary;
the engine uses theta=k_B T in Ha. Gaussian, MP and cold schemes reject.
Odd unpolarized electron counts are permitted only in FD mode.

The FD band count (per spin channel and k point) is `scf.nbands` when given,
otherwise `band_cap` (default 32). `band_cap` stays the documented cap: an
`scf.nbands` above an explicit `band_cap` is rejected at preflight
(`scf.nbands = N exceeds occupations.band_cap = C`). Before FD-FORCES
(2026-11-26) `scf.nbands` was ignored in FD mode and the cap was always the band
count; inputs that set a smaller `nbands` now solve fewer bands. The policy
prepares and solves every one of these bands. The band count must fit every k basis. Insufficient
empty-state coverage fails with an instruction to increase the cap. The
reported top-band diagnostics are not bounds on infinitely many omitted bands;
repeat with at least four more bands and compare F/density.

When to use it for an insulator-like cell: if fixed occupations stall (density residual
floors around 1e-3, the SCF wanders) and the printed HOMO/LUMO gap is below about
1e-3 Ha, the integer occupation cut lies inside a (near-)degenerate level cluster and
the fixed-occupation problem is ill-conditioned (elongated and slab-like cells; the
same cells defeat Quantum ESPRESSO, see ledger rows ELONGATED-GAP, ELONGATED-SMEARING and
SLAB-SMEARING). Use FD with a width above the cluster splitting (0.02 to 0.1 eV converged
both measured cells in 20 to 35 iterations) with the strict profile, and `band_cap` at least
occupied bands + 5. These are historical strict-profile observations; the current
reference standard path supports FD, except the G4D oracle fallback. A non-converged
fixed-occupation run prints a hint when the final gap is below 1e-3 Ha.

For collinear FD, omit channel counts and `total_moment` for one shared mu and
free moment, or supply `total_moment` (possibly fractional) or matching
`n_up`/`n_down` for two roots. Positive-temperature sector endpoints reject.
Use a full mesh (`basis.ktime_reversal = false`) and `lda-pw` or `pbe`.

Results contain occupations in spin/k/band order, mu, S/k_B, E, positive TS,
F=E-TS, separate E/-TS/F ledger lines, a fresh finite-T certificate, and the
proposed `state_contracts` extension. Band-gap output becomes explicit
metallic/Fermi-level information. Total density is written as Gaussian cube;
collinear channel and magnetization cubes are also written, and all channel grids
are retained in JSON. P20A/P22A add strict certified forces/stress; unsupported derivatives remain unavailable, not zero.
Fixed-policy results retain their original output path.

Reference FD supports runtime control and exact checkpoint v2. Legacy
band/DOS/volumetric options and checkpoint v1 FD restart remain unavailable;
XC density floor and spectral-guard input keys were removed.
A fixed-occupation checkpoint is never silently converted to FD.
See [occupation validation](validation/occupations.md) for measured scope,
root algorithm, thresholds, ignored Al command and remaining integration work.

## P19A crystal symmetry (opt-in)

```toml
[symmetry]
mode = "auto"       # default "off"
tolerance_bohr = 1e-7
```

`auto` detects species-aware spatial operations, uses the subgroup preserving
`basis.kgrid`/`kshift`, and reports irreducible eigenvalues, weights and integer
multiplicities in an additive `results.json.symmetry` block. It includes the
operation list and complete full-to-irreducible orbit map. `off` skips spatial
reduction, retaining the separate `basis.ktime_reversal` setting.

Density is symmetrized on the full FFT grid, with fractional-translation phases.
Explicit `basis.fft_grid` dimensions are used exactly when every selected density
operation permutes the grid; incompatible grids reject with `symmetry.fft_grid`.
FFT-friendly rounding applies only to automatic selection: the existing chooser
selects the smallest compatible 2/3/5-smooth dimensions along each axis, using a
common dimension if rotations require it. Preflight records the dimensions,
alias-free lower bounds and symmetry reason. `basis.fft_grid_overridden` means
user-specified, including when the explicit dimensions equal the automatic choice;
it is false for automatic symmetry padding. Both CLI drivers use this preparation.
Compare full and reduced calculations with identical FFT dimensions (Si diamond
30 Ry uses 32^3 with auto, versus the symmetry-off 25^3 default; explicit 28^3 is
used exactly and 25^3 rejects). Extremely skew cells or incommensurate origins
can reject with an actionable diagnostic; no invalid reduction is silently accepted.

Unpolarized fixed occupations and Fermi–Dirac are supported. Collinear input
requires `basis.ktime_reversal=false`; only spatial identity-spin operations
act on each channel. Spin swaps and magnetic antiunitary operations are disabled
without a model-invariance proof. Fields, local constraints and Hubbard input
are unsupported and rejected. Reference auto symmetry supports observation/control and checkpoint v2, whose
state retains prepared symmetry. Oracle auto symmetry still requires Classic
execution and its v1 limitations remain. Band and density file formats
are unchanged. See [validation and limits](validation/symmetry.md).

## Atomic forces (P20A)

Set `[output] forces = true` to request Cartesian forces in input atom order,
Ha/Bohr. A passed fixed-point verification (the default single path, ADR 0043)
or the strict certificate (`BAANDI_REFERENCE_STRICT=on`) is required;
`scf.verify = "off"` rejects at preflight.
`results.json.forces` contains `status`, `quantity` (E or F), `units`, `frame`,
`values`, `per_term` (ewald/local/nonlocal/core), `net_force`,
`maximum_absolute_component`, `elapsed_seconds`, `band_passes` and
`unavailable_terms` and a named `derivative_ledger`. Missing derivatives never become zero or a partial total.
The fixed unpolarized path prints the table; the strict-oracle spin and FD
serializers also write `forces.txt`.
With symmetry auto the same SCF spatial subgroup averages each term.
For FD these are derivatives of free energy at fixed temperature and charge.

FD-FORCES (2026-11-26): the default reference driver computes forces and stress
for Fermi-Dirac (scalar and collinear) and collinear fixed-sector runs too.
Before this change the reference spin/FD result (`output_profile =
"reference_spin_fd_v1"`) silently omitted both blocks. `forces.quantity` and
`stress.quantity` are `F` for FD (Mermin free energy F = E - TS; the
occupation derivatives add no term because F is variational) and `E` for a
fixed spin sector. The nonlocal and kinetic terms sum f_nk w_k over both spin
channels. As for fixed occupations, the passed reference verification is the
certificate (`verification = "verified"`, standard `preconditions`); with
`scf.verify = "off"` the blocks say `unverified`. A non-converged run writes
`{"status":"unavailable","reason":"scf_not_converged"}` for each requested
block. Both console tables are printed. The stress block includes the Pulay
assessment from three reference spin/FD cutoff-scan SCFs with the same
occupation policy.
State-contract extension version 2 includes the available-force branch and requires
strict/converged/certificate-passed preconditions. Version 1 unavailable documents
remain valid, and unavailable forces never carry numerical arrays. The outer
results schema version is unchanged. See [validation](validation/forces.md).

## Fixed-cell relaxation (P21A)

`baandi relax input.toml --output-dir DIR [--resume]` accepts the optional
`[relax]` table, forces from a verified standard SCF on every ionic step (the
strict certificate with `BAANDI_REFERENCE_STRICT=on`; see the single SCF path
section), capped BFGS/FIRE steps and Cartesian constraints. Exit 0 is converged, 3 is ionic step exhaustion, 4 is failure after
one fresh half-step retry. Density transfer is labelled `warm_transfer`; the
optimizer restart is a separate same-build ionic continuation record.
The `calculation_type="relaxation"` results profile contains the `relaxation`
block and input/executable/UPF provenance. See [input, units and output-format
losses](user/relaxation.md) and [validation scope](validation/relaxation.md).


## Stress and pressure (P22A)

Set `[output] stress = true` on a converged NC SCF: since ONEPATH-2 the default
standard path with a passed fixed-point verification (see the single SCF path
section above), or the strict certified reference path
(`BAANDI_REFERENCE_STRICT=on`); `scf.verify = "off"` rejects this request.
Spin and Fermi-Dirac stress runs use the standard path with the verification
since ONEPATH-3. The
Cartesian tensile-positive 3x3 tensor is in Ha/Bohr^3;
pressure is `-trace(stress)/3`. The JSON `stress` block contains per-term values,
symmetric/antisymmetric parts, raw antisymmetric norm, derivative ledger and
preconditions. LDA-PZ, LDA-PW, PBE, collinear spin, NLCC and FD free-energy
stress are supported in the corresponding existing engine/CLI models.

Only frozen-Miller stress is reported as the stress. The separate isotropic
Pulay estimate performs three additional SCFs at 100%, 105% and 110% cutoff,
reports a two-step uncertainty diagnostic, and warns above 1e-5 Ha/Bohr^3 or
when no shell is sampled. It never silently corrects the tensor and never
asserts cutoff convergence. This can cost more than the analytic pass.

Fixed unpolarized runs print a stress table; spin/FD write `stress.txt` as well
as JSON. ASE `get_stress` uses tensile-positive eV/Angstrom^3 and order
xx,yy,zz,yz,xz,xy for the existing Gamma/unpolarized adapter. See
[conventions, equations, protocols and evidence](validation/stress.md).

## Variable-cell relaxation (P23A)

`baandi relax` accepts `relax.mode="ions"|"cell"|"ions_cell"` (default ions),
`cell_dofs="full"|"isotropic"|"xyz_diagonal"`, Cartesian `fixed_axes`, external
`pressure_ha_bohr3`, `stress_tol_ha_bohr3`, `rebuild_every`, and
`pulay_warning_ha_bohr3`. Cell modes minimize E+PV on frozen Miller segments,
then report a separate final nominal-cutoff single point. See
[full conventions, restart and output semantics](user/cell-relaxation.md) and
[measured validation and limits](validation/cell-relaxation.md).

For nominal-cutoff material equilibrium, follow the periodic-rebuild and
fresh-reference polish recipe in [cell relaxation](user/cell-relaxation.md#recommended-material-equilibrium-protocol).
`relaxation.rebuilt_pressure_check` warns when the final rebuilt pressure misses
the requested external pressure by at least `stress_tol_ha_bohr3`. This check is
separate from the Pulay tensor-difference threshold; frozen-only convergence
alone is not a nominal-cutoff equilibrium claim. The default `rebuild_every=0`
remains unchanged.
