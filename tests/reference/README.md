# QE reference interchange

Real QE reference bundles now exist (`si30_ours`, `si30_auto`, `scan_e*`); this opening paragraph describes the original tooling change, which shipped before any real output existed. `../fixtures_format/` contains
invented FORMAT FIXTURES, not calculations. Never move those files here or
remove their `is_format_fixture: true` marker to call them references.

Use the remote runner and original MPI build described in
[reference-environment.md](../../docs/reference-environment.md). No local QE
installation, QE execution, real UPF inspection, or network access is needed by
these tools. Python 3.9+ standard library is sufficient; use `python3 -I`.

## Generate inputs

```sh
python3 -I scripts/make_qe_input.py --results runs/ours/results.json \
  --pseudo-dir /remote/pseudo --upf-name Si.pz-vbc.UPF --out qe.in
```

The default `--fft-grid ours` explicitly requests `nr1,nr2,nr3`; `--fft-grid auto`
omits all three. Produce both to study QE's grid choice. `--ecut-ry` overrides
both wavefunction cutoff and the density cutoff (always four times it). An
override differing from the source result will intentionally fail comparison;
run our engine at that cutoff first for a matched pair. Other options are
`--conv-thr` (default 1e-12 Ry), `--mixing-beta` (0.7, range (0,1]),
`--prefix` (`rust_dft`), and `--outdir` (`./qe-work`).

The single-species schema-1 results supply the resolved bands, electrons,
occupations, cutoff, grid, species digest, and geometry. Our row-indexed matrix
`A[row][column]` stores lattice vectors as columns; the generator transposes
these into three QE `CELL_PARAMETERS bohr` rows. All geometry uses 17 significant
digits, round-tripping binary64 exactly. Fractional positions stay unwrapped.
No Angstrom conversion is involved. Atomic mass 1.0 is an SCF-only placeholder.
`K_POINTS automatic` with `1 1 1 0 0 0` selects the full complex Gamma path;
`K_POINTS gamma` uses a different real-wavefunction representation and is refused.

The generated input uses `ibrav=0`, fixed occupations, `nosym`, `noinv`, Davidson
with full accuracy, 200 electron steps, low disk I/O, high verbosity, and disabled
forces/stress. `input_dft` is deliberately absent: the UPF functional must be
reported as LDA-PZ. QE's `conv_thr` has different semantics from our three SCF
criteria; the generator does not pretend to translate them. Our actual criteria
remain in the source result. Mixing options are convergence controls, not
claims of identical SCF trajectories.

`qe.in.meta.json` records the input and source-results SHA-256, expected UPF
SHA-256, UPF basename, and every generator option (including null cutoff override).
The generator never opens a UPF. The remote runner must hash the actual file
used; the local digest comes only from `results.json`. Paths containing quotes,
newlines, `!`, or `&` are rejected rather than ambiguously escaped.

## Returned bundle

Return a directory containing:

- `qe.in`: exact input bytes actually executed, without post-run editing.
- `qe.out`: verbatim stdout of one completed, converged SCF calculation.
- `provenance.json`: all fields below. Preserve the input sidecar separately when available.

| Required field | Type / rule |
| --- | --- |
| `command` | Nonempty exact executed command string, including launcher/environment setup; record `OMP_NUM_THREADS=1 MKL_NUM_THREADS=1` |
| `qe_version` | Full trimmed `Program PWSCF v.X starts ...` line, exactly matching stdout |
| `host` | Nonempty free text identifying the execution host |
| `mpi_processes` | Positive integer |
| `omp_threads` | Integer 1, per remote protocol |
| `date_utc` | ISO-8601 timestamp with UTC offset (`Z` or `+00:00`) |
| `input_sha256` | Lowercase SHA-256 of returned `qe.in` bytes |
| `upf_sha256` | Lowercase SHA-256 of the UPF actually used |
| `wall_time_s` | Finite nonnegative measured seconds; wrapper timing may differ from QE timing |
| `returned_by` | Nonempty runner identity |
| `notes` | String; describe deviations and provenance caveats; empty allowed |
| `is_format_fixture` | Boolean, false only for actual reference runs (temporary synthetic self-tests are the sole test exception) |

Additional provenance fields such as `maxrss` are retained in reports. Record
the original MPI build and exact executable path in `command`; hashes and text
are attestations by the runner, not cryptographic proof of execution. The tool
cannot authenticate the executable or verify an inaccessible remote UPF.
Duplicate JSON keys, nonfinite values, missing fields, input hash changes, UPF
hash differences, and version disagreement are refused. A true fixture marker
immediately refuses numerical work with `FORMAT FIXTURE: not a reference`.

## Verify and report

```sh
python3 -I scripts/verify_qe.py --ours runs/ours --bundle tests/reference/returned \
  --json comparison.json --markdown comparison.md
```

Defaults are `qe-comparison.json` and `qe-comparison.md` in the current directory.
Both reports are written only after parsing/provenance/conditions succeed.
Refusals do not overwrite old reports: use fresh report paths and always check
the exit code. Evidence paths cannot be used as report destinations.

| Exit | Meaning |
| --- | --- |
| 0 | Comparison completed; all noninformational tolerances passed |
| 1 | Usage, malformed our-results data, or I/O/report error |
| 2 | Provenance incomplete, malformed, or mismatched |
| 3 | Required output/input pattern absent/unsupported, or conditions mismatch |
| 4 | FORMAT FIXTURE: not a reference |
| 5 | Comparison completed; at least one noninformational tolerance exceeded |

Before energy comparisons, require matching wavefunction/density cutoffs in Ry,
electron/band counts, plane-wave count, species labels, complex Gamma sampling,
fixed occupations, and LDA-PZ. Lattice components and Cartesian positional
separations must agree within 1e-12 Bohr (no wrapping or fitted translation).
When the input contains grid dimensions, require all three to equal ours and
require stdout to report them. When all are absent, record both grids and an
automatic-grid warning; grid equality is deliberately not a gate. Our SCF must
be certified converged. All discovered condition failures are listed together;
malformed syntax fails at the offending pattern before semantic checks.

The input parser accepts the generator's three namelists in order, one supported
assignment per line, followed by its four cards in order. Unknown/duplicate
assignments, extra cards, alternative units, gamma-only tricks, multiple species,
spin settings, or `input_dft` overrides are rejected. Remote runners should use
the generated inputs, not hand-reformat them. Text metadata and SCF thresholds
may be changed only before execution and must be reflected in the returned hash.

### Exact stdout subset

This parser is independently written for text interchange, without copying or
translating QE source. It has not yet been tested against a real returned output.
Recognized records must be complete and unique; unrelated banner/iteration text
is ignored. Unsupported spelling or layout fails with the expected regex; there
is no inferred/default energy, cutoff, grid, plane-wave count, or eigenvalue.
Required records are:

- One `Program PWSCF v.<version> starts ...` line.
- `number of electrons = N`, `number of Kohn-Sham states = N`,
  `kinetic-energy cutoff = X Ry`, `charge density cutoff = X Ry`, and
  `number of k points = 1`, each on its own line.
- `Dense grid: ... FFT dimensions: ( N, N, N )` on one line. The dense
  (charge-density) dimensions are used, not a smooth-grid record.
- `Exchange-correlation = SLA PZ NOGX NOGC` (or `SLA PZ`), optionally followed
  by `( 1 1 0 0 )`; whitespace may vary. The complete printed string is reported.
- `convergence has been achieved in N iterations` and
  `End of self-consistent calculation`.
- After that end marker, exactly one
  `k = X X X ( N PWs) bands (ev):` line, with zero k components, followed by
  exactly `nbnd` sorted numeric eigenvalues over one or more lines. A blank line
  or following nonnumeric record terminates the numeric block.
- Exactly one line starting `! total energy = X Ry`, followed later by the
  four `one-electron contribution = X Ry`, `hartree contribution = X Ry`,
  `xc contribution = X Ry`, and `ewald contribution = X Ry` records.
- `JOB DONE.`. Multiple concatenated runs, partial output, and unconverged runs
  are refused.

Optional `highest occupied level (ev): X` or
`highest occupied, lowest unoccupied levels (ev): X Y` and
`PWSCF : <duration> CPU <duration> WALL` are retained if present. A present but
malformed optional record fails; absent optional values are JSON null. Wall
units are `d`, `h`, `m`, `s`, in order (seconds may be fractional). The comparison
gap is computed from the band eigenvalues; optional printed levels are reported
without substituting them. Missing empty bands explicitly mark the gap unavailable.
Numbers permit E/e/D/d exponents and must be finite.

### Numerical rules

| Comparison | Default absolute tolerance | Option |
| --- | --- | --- |
| Total energy per atom | 1e-4 Ha/atom | `--tol-energy` |
| Each energy decomposition term | 1e-4 Ha | `--tol-term` |
| Each eigenvalue, all bands | 1e-3 Ha | `--tol-eigenvalue` |
| HOMO-LUMO gap | 1e-3 Ha | `--tol-gap` |
| Electron count | 1e-8 relative (absolute bound = 1e-8 times electrons) | `--tol-electrons` |

1 Ha = 2 Ry exactly. Eigenvalue conversion uses CODATA 2018:
1 Ha = 27.211386245988 eV. QE's own constant may differ in its last digits;
eigenvalue differences at the 1e-9 relative level are irrelevant to the default
tolerances. Decimal stdout rounding also limits available precision.

Map QE one-electron to our `kinetic + local + nonlocal`, and Hartree, XC, Ewald
to identically named `energies.*.ha`. Do not add diagnostic band sums or
zero-mode terms again. Compare raw reported terms with no fitted shifts or
constant corrections. If independent convention analysis establishes a different
decomposition, supply `--decomposition-note 'specific reason'`: all four terms
remain visible as PASS/FAIL but explicitly informational; they do not determine
overall acceptance. This never exempts total energy or eigenvalues. Conventions
must not be declared different merely because values disagree.

Reports include both values, absolute differences, units, tolerances, PASS/FAIL,
maximum eigenvalue/term differences with attaining items, and the largest
difference/tolerance-ratio item in JSON. Every nondefault tolerance is flagged.
Tolerances must never be loosened to make a comparison pass.

QE's electron record is the configured count, not an independent density
integral. Our final relative integrated-charge error is also tested at 1e-8;
our occupied eigen residual is checked at 1e-7 Ha. The project orthogonality
target is 1e-8, but neither input data source exports the required overlap
matrix: this tool explicitly reports it as not assessed. No density comparison
or all-band residual validation is claimed.

For the candidate UPF digest
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`, every comparison
header contains **Pseudopotential provenance is unverified**, regardless of the
comparison outcome. See [candidate notes](../../docs/candidate-si-upf.md).

## Tests

```sh
python3 -I -m unittest discover -s scripts/tests -v
```

All stored fixtures use invented conspicuous values and `is_format_fixture: true`.
A passing comparison is exercised only inside temporary tests which build
self-consistent pairs from our own synthetic numbers, temporarily set the bundle
marker false, and label reports **synthetic self-test, not a reference**.
The real-binary scan test generates its synthetic UPF in a temporary directory;
no UPF is committed. It uses `target/debug/dft-cli` or the `DFT_CLI` environment
variable and passes it through the driver's `--dft-cli` argument. If absent,
the test skips with instructions to build via `scripts/check.sh`.

## Real bundle regression

The earlier fixture-only notes above describe the initial tooling release.
`si30_ours/` now contains the real returned QE 7.6 bundle. The regression in
`scripts/tests/test_qe_real.py` protects its printed values, full comparison,
unverified pseudopotential warning, singular/plural level labels, and refusals
for altered provenance, missing energy, and mismatched conditions. The bundle
must never be edited: all mutation tests use temporary copies, and all reports
are written to temporary directories. No external UPF or QE executable is needed.

## Fourier density comparison (P6D)

Run the standard-library tool in isolated Python, with **fresh** report paths:

```sh
python3 -I scripts/compare_rhog.py \
  --ours-dir tests/reference/si30_ours/ours \
  --qe-summary tests/reference/si30_ours/density/qe_rhoG_summary.json \
  --json-report /tmp/rhog-report.json --markdown-report /tmp/rhog-report.md
```

For all-coefficient comparison replace `--qe-summary ...` with
`--qe-file /path/to/charge-density.dat`. The strict reader supports only the
specified little-endian, four-record, non-HDF5, full-complex (`gamma_only=0`),
single-spin format. No QE installation or external Python packages are needed.

Our directory must contain `density.bin`, `density.json`, and `results.json`.
The tool checks the density digest and integrated charge (absolute error at most
1e-10 electrons), then independently computes a normalized separable DFT.
It reconstructs our density sphere from the lattice and `ecutwfc_ha` in the
results (density cutoff = 4 times the wavefunction cutoff, relative boundary
allowance 1e-12), checking its reported count. It does not select our G set from
QE's list. Each axis requires `-(n//2) <= m <= (n-1)//2`, with negative even
Nyquist and distinct wrapped indices. Full-file mode reports both set differences;
a mismatch exits 3. Summary mode checks the count and listed indices, but cannot
certify equality of the complete QE set because that list is absent from the
stored summary. Summary sums and maxima are compared over our density sphere.

Default gates: coefficient and nonzero maximum magnitude differences <= 1e-7
 e/Bohr^3, rho(0) difference <= 1e-12 e/Bohr^3, relative sum-of-squares difference
(using QE as denominator) <= 1e-6, and full-file relative L2 difference (using our
norm as denominator) <= 1e-6. Override with `--coefficient-tol`, `--rho0-tol`,
`--parseval-tol`, or `--relative-l2-tol`; reports identify non-default values and
explicit override flags. Different maximizing Miller triplets do not fail the
magnitude comparison. The maximum pairwise complex distance among the available
(1,0,0), (0,1,0), (0,0,1) samples is an informational symmetry noise estimate.
Exit codes are 0 (pass), 5 (tolerance failure), 3 (inconsistent input/format),
and 1 (usage/I/O). Existing output paths are never overwritten.

`python3 -I -m unittest discover -s scripts/tests -v` includes analytic DFT,
malformed-record, CLI and real-summary regression tests. The latter prints the
measured differences and uses only repository evidence. The original summary's
provenance limitations, including unverified pseudopotential provenance, remain
in effect; a passing density comparison does not resolve them.
