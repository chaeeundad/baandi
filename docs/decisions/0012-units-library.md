# ADR 0012: Boundary unit conversions through uom

Status: Implemented for P7U-UNITS; pending review (2026-11-03).

## Decision

Use the already vendored uom 0.38.0 (Apache-2.0 OR MIT), with defaults disabled
and `f64`, `si`, `std` enabled. uom owns all measured constants and SI scaling.
The implementation inspected its energy, length, pressure and ratio definitions.
The engine continues to use plain f64 in Hartree atomic units. No numerical
kernel, convergence tolerance, temperature unit or magnetic unit changes.

Keep Bohr, Angstrom, Hartree, Rydberg and ElectronVolt and their existing From
implementations. Preserve the public `const f64` conversion-factor API using a
small build script: uom quantities produce Rust constant declarations in OUT_DIR.
This uses uom as both a build and runtime dependency; no copied CODATA literals
enter workspace production sources. Keeping const semantics avoids breaking
consumers in concurrent tickets. The generated constants are not committed.

EnergyUnit, LengthUnit and PressureUnit select the boundary units. Conversion
functions take a value, source unit and destination unit. `to_hartree`,
`from_hartree`, `to_bohr`, `from_bohr` are convenience functions. Pressure is
derived dimensionally from energy divided by three lengths. uom has neither a
Rydberg energy nor a kilobar unit nor a millielectronvolt unit: Ry = Ha/2 is an
**exact definition**, the narrow ADR 0003 physics-glue exception, not a measured
constant. Milli/kilo scaling comes from uom's `per_mille` ratio quantity and
its electronvolt/bar units, without hand-typed decimal conversion factors.

The CODATA 2018 compatibility test checks the validated MVP Bohr/Angstrom and
Hartree/eV values within one last printed digit. SI-based Hartree/eV division
can differ from the previously rounded eV constant at f64 roundoff; bitwise
identity of every eV output is not promised. Bohr/Angstrom retains the MVP
value. The independent pressure test uses published CODATA 2018 SI numbers
only in an explicitly allowlisted integration-test file.

## Parsing and CLI compatibility

Strict typed errors retain the offending input and distinguish empty input,
invalid number, nonfinite number, missing unit, unknown unit and extra tokens.
Decimal and scientific notation may adjoin the unit (`30Ry`, `1eV`). Unit aliases
are case-insensitive and form a closed list documented in docs/cli.md. No unit
is inferred from a value. In particular `meV` means millielectronvolt under the
case-insensitive policy; megaelectronvolt is unsupported.

Config::parse normalizes the additive string forms into the unchanged public
configuration structs before strict deserialization and validation. Existing
numeric fields and all output fields remain. Supplying both forms is an error,
even if equivalent. Existing defaults apply only when neither optional tolerance
form is given. Direct serde deserialization retains the legacy schema; the
unit-aware public entry point is Config::parse. The lattice alternative is an
explicit matrix of unit strings, preserving arbitrary cells without inventing
an implicit lattice shape from a scalar lattice constant.

The convert subcommand supports energy, length and pressure within the same
dimension, prints scientific notation with 15 significant digits, and returns
usage code 2 for malformed numbers/units, dimension mismatch or nonfinite output.
Results schema version 1 gains `total_energy: {ha, ry, ev}`, all from the same
engine Hartree total through dft-units. Existing energy audit fields remain.
The schema documentation lives in docs/cli.md; the agent report schema is
unrelated and unchanged.

## Precision and enforcement

All conversions use f64 uom quantities. Every unit-pair round trip is tested
with relative error <1e-14 for representative finite normal values, including
negative values and widely separated scales. The original wrapper round-trip
bound remains <1e-15. These bounds exclude underflow/overflow and subnormals;
raw numeric conversion functions follow IEEE-754, while string boundary helpers
and the CLI reject nonfinite inputs/results. Zero is tested separately.

The standard-library Python guard scans workspace .rs sources, including
build.rs and dft-units/src. It detects digit prefixes independent of underscore,
decimal-point and exponent placement. Only named reference files under tests/
are allowlisted. Rust inline tests are scanned too; references belong in the
allowlisted test file. Documentation is outside the .rs scan. Generated target
files, vendor sources and metadata are excluded. A temporary source copy is
used to plant literals and assert guard exit 1 with path/line diagnostics.
This guard detects literals, not deliberately disguised arithmetic expressions.

## Dependency review

Claude supplied the pinned, licensed, vendored crate before implementation.
uom is a general-purpose pure-Rust dimensional-analysis library, avoiding FFI
and target-specific numerical backends. The offline host build/tests exercise
Apple Silicon compatibility. Transitive dependencies are num-traits and typenum,
with autocfg as num-traits' build dependency. This ticket does not independently
claim current download counts or release activity beyond the supplied 0.38.0
source and lockfile; future updates must rerun the CODATA compatibility tests.
