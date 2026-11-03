# ADR 0026: Shared CLI input preflight

Status: proposed implementation, pending review

## Context

The feature-gap survey T2 demonstrated that species labels could silently select
another element's potential, and validate did not reject invalid SCF settings.
The CLI must reject those inputs before expensive SCF construction.

## Decision

Use one CLI preflight entry point for validate and run. Keep numerical kernels
unchanged. It resolves the engine defaults through run::options, uses the existing
functional matcher and basis/grid constructors, and mirrors the public SCF
construction constraints without constructing a Hamiltonian or running SCF.
The CLI additionally rejects zero threads (the engine otherwise normalizes them).

Retain string pseudopotential paths; accept strict inline entries with path,
optional element, and default-false allow_element_mismatch. Normalize these in the
input boundary while preserving the existing path map used by parallel tickets.
Explicit elements must be recognized chemical symbols. Whitespace is trimmed from
the UPF element; symbol case remains significant. Overrides never modify physics.

Expose the module under config::preflight to avoid editing the unowned library
root. Config wiring adds the module, skipped identity map, and normalization calls.
Run wiring adds the shared check/print before preparation and one additive result
field. Main wiring adds --json and routes validate through the shared check.
No UPF changes are needed: header values and suggested wfc_cutoff are already public.

Warnings and overrides are structured string lists in a versioned report, printed
for human use and included in results.json. The minimum cutoff warning compares
positive finite wfc_cutoff in Ry to the input wavefunction cutoff; it does not infer
convergence or recommend from rho_cutoff. Memory is explicitly a partial payload
estimate, not peak RSS. UPF radial data is parsed; dense orbital/density files are
never read by validation.

## Consequences and verification

The previous permissive behavior is intentionally incompatible. Existing tests
asserting that nonexistent files validate or that Si may silently select synthetic
Xx were migrated after explicit orchestrator approval to edit tests/commands.rs
and tests/run.rs. Example templates still need real local potential paths.

New process tests use the existing generated fixture with a synthetic Si header;
no real UPF is committed. They exercise error key paths, custom labels, overrides,
SCF controls, odd electrons, basis capacity, warnings, report shape and persistence.
The engine remains the final validation authority; future public SCF option changes
must update preflight parity tests. Preflight enumerates basis/grid mappings and is
not an allocation-free resource-budget guard. No claim of QE numerical validation
is made by these input tests.

## Measured review evidence (2026-11-03)

- `cargo test --offline -p dft-cli --test preflight -- --nocapture`: exit 0,
  six tests passed, 0.62 seconds in debug. All three example templates validate
  after replacing their placeholder path with the generated synthetic Si file.
- One-step synthetic CLI versus direct ScfProblem: energy difference 0 Ha
  (assertion <1e-12 Ha), density.bin exactly equal to direct density bytes.
  This is a wiring equivalence test, not a converged material result or QE reference.
## Migration of existing tests

The orchestrator approved the two additional test files during this task. Each
existing test retains its numerical or process-boundary purpose:

- `commands::example_commands_and_failures` supplies a temporary synthetic Si
  potential, checks the material/electron summary, then removes it and asserts
  an early missing-potential error. Grid inspection and missing-input/usage checks
  remain in place.
- Run fixtures change only the Xx header to Si; radial arrays and numerical
  assertions are unchanged. The ignored two-species gap test uses H headers and
  explicit H elements for the custom A/B labels, preserving both distinct radial
  potentials and the expected physical gap failure.
- `verify`, `kpoint_scf`, and `expressions_preserve_energy_and_record_resolutions`
  include the additive preflight field in exact result-schema checks.
- `failures_and_paths` expects exit 1 for beta=100, unsupported PAW/USPP, XC
  mismatch and missing UPFs: all are now rejected before the former solver/parser/
  I/O stages. Budget exhaustion remains exit 3 and output I/O errors exit 5.
- `functional_compatibility_and_original_label` and
  `pbe_selection_override_and_spin_rejection` expect exit 1 for incompatible XC,
  and check the new key-path/remedy message. Original label preservation and
  successful override computations remain tested.
- `threading_metadata_and_identical_outputs` now asserts that zero threads fail
  without creating outputs. Positive-thread runs still require identical numerical
  JSON and density bytes; only the descriptive thread-count summary is removed
  after checking that it accurately describes the requested count.

Initial checks failed on the old contracts (one command test and nine run tests).
After migration, `cargo test --offline -p dft-cli --test commands --test run --
--nocapture` exits 0: 5 command tests and 10 run tests pass, 2 release-only tests
remain excluded from debug. The small converged synthetic case gives
-0.7275518620715111 Ha in 13 iterations, direct-engine difference 0 Ha, and
charge error 3.553e-15 electrons. Expression and quantity-expression energy
differences are both 0 Ha. These are synthetic tests, not QE references.

README, artifacts and tickets remain untouched under the explicit ticket scope.
Changes are uncommitted as requested.

Final verification after the authorized migration:

- `scripts/check.sh`: exit 0; build, test, fmt and clippy all PASS.
- `cargo test --offline --release -p dft-cli --test run -- --ignored --nocapture`:
  exit 0; both release-only tests pass (0.48 seconds test execution).
  FCC synthetic energy -4.6830234053637119 Ha, direct-engine difference 0 Ha,
  charge error 6.217e-15 electrons, 20 iterations. All three legacy/new Gamma
  comparisons have energy difference 0 Ha. The four k-point forms agree within
  2.2e-14 Ha (required <1e-10 Ha). The deliberately gapless case still fails
  numerically with gap -0.0669426243347313 Ha (required 1e-7 Ha).
- `git diff --check`: exit 0.

Shared-file integration locations: config.rs module declaration, Config's skipped
species_identity field and Config::parse normalization (six added lines); run.rs
execute's preflight check/print and the additive results.json field (two added
lines and one changed line); main.rs Validate's --json option and validate dispatch,
with expression reporting preserved for text output. No manifests, dependencies,
UPF parser, numerical kernels, or public engine APIs changed.
