# ADR 0015: Resolve CLI expressions before typed configuration

Status: Implemented for P7E2-EXPR-CLI; pending review (2026-11-03).

`Config::parse` walks an explicit numeric-field allowlist in the TOML value tree
before existing typed deserialization and validation. Numeric TOML values are
left untouched. Strings use only the dft-units evaluator from ADR 0014; integer
strings use `eval_expr_integer`, including its documented integrality tolerance
and i64 range. Unsigned consumers reject negative resolved integers. No parser,
dependency or numerical-kernel changes are introduced. Array shapes, unknown
keys and nonnumeric fields remain checked by serde and existing validators.

Unit alternatives use `parse_quantity_expr` and existing conversions. Their
resolved records use the destination field's unit: cutoff in Ry, energy
tolerances in Ha and unit-bearing lattice entries in Bohr. Other records use
the declared input field's units. All numeric-field strings, including numeric
strings and legacy quantity strings, are recorded without syntax heuristics.
Each record contains `key_path`, original `text`, and resolved f64 `value`.
Integer records store the accepted integer as f64 (large integers may lose
precision in this provenance representation; the typed integer is retained).
The schema-version-1 `input_expressions` array is omitted when empty, preserving
plain-input output key sets. Validate prints the same records using round-trip
f64 display. Original input hashing remains unchanged.

Cartesian atom coordinates are an additive `cartesian` alternative to
`fractional`, in `system.lattice_unit` (Bohr with `lattice_vectors_with_units`).
They accept mixed numbers and expression strings and convert through Cell's
existing inverse transformation. Both coordinate forms together are rejected.
Provenance retains original Cartesian components, before conversion to fractional
coordinates; output atoms continue to contain fractional and Cartesian Bohr values.

`eval` evaluates a unitless expression. `convert` evaluates its value argument
then requires both source and target units. Both print scientific notation with
15 significant digits and return 2 on invalid input, 0 on success. Existing
validate/run error statuses remain unchanged. Evaluation failures identify the
original key path (including zero-based array indices) and quote the input.
No expression variables other than pi/e or units within arithmetic are added.

Verification uses only synthetic UPF fixtures, comparing resolved fields and
SCF energies at absolute tolerance 1e-12 Ha, checking plain-output key sets,
integer rejection, indexed diagnostics, mixed arrays and command exit codes.
