# ADR 0014: Arithmetic expressions at the units boundary

Status: Implemented for P7E-EXPR; pending review (2026-11-03).

## Library and scope

Use the supplied, vendored exmex 0.21.0 (MIT OR Apache-2.0), with no default
features and no optional features. evalexpr is rejected because of AGPL-3.0.
Exmex owns expression tokenization, grammar, precedence and evaluation. Its
operator factory is filtered to a closed vocabulary. Its provided
`literal_matcher_from_pattern!` regex extension recognizes decimal and scientific
literals, including `1e-13` (not covered by its default numeric matcher).
There is no handwritten expression parser. Whitespace is mapped to ASCII spaces
without removing token boundaries. A character allowlist excludes alternate
variable syntax, assignments, commas and unsupported punctuation.

This is an additive dft-units API only. Existing APIs and CLI behavior are
unchanged. P7E2 will connect TOML after P7A2; no CLI, kernel, serialization or
file-I/O work is included here. No units inside expressions, user variables,
user functions, symbolic algebra or implicit multiplication are supported.

## Vocabulary and evaluation

Case-sensitive operators are `+ - * / ^`, unary minus and parentheses.
Functions are exactly `sqrt`, `abs`, `sin`, `cos`, `tan`, `asin`, `acos`, `atan`,
`exp`, `ln`, `log10`. Arguments and results of trigonometric functions use radians.
Exmex also accepts unary function application without parentheses, such as
`sin pi`; parentheses are recommended for clarity. Constants are only `pi` and
`e`, bound through exmex's variable list to `std::f64::consts::PI` and `E`.
The default exmex constants and all other functions are excluded. Every remaining
variable is an error, including one in an algebraically cancellable expression.
Plain numeric strings use Rust f64 parsing first, preserving scientific notation,
signed zero and signed numeric literals such as `+2`; unary plus on expressions
is otherwise excluded. No locale-dependent number formats are accepted.

`eval_expr(&str) -> Result<f64, ExprError>` rejects nonfinite final values,
including `1/0`, `sqrt(-1)`, `asin(2)` and literal NaN/infinity. Finite underflow
follows f64 behavior. The contract checks the final result, not every intermediate.
Parsing uses `parse_wo_compile` to avoid algebraic compilation/reassociation.
Tests establish that powers associate left: `2^3^2 = 64`; unary minus binds
more tightly than power: `-2^2 = 4`, while `-(2^2) = -4`. Multiplication/division
precede addition/subtraction; equal-precedence arithmetic associates left.
These are exmex rules, not assumed mathematical conventions.

The evaluated f64 is the value a future consumer must use, without reparsing or
rounding the displayed value. `resolve_expr` returns `ResolvedExpression { text,
value }` with the original text and that same f64 for future results.json
provenance (`expression text -> value`). No serialization dependency is added.
Transcendentals follow std's host math; cross-platform bitwise equality is not
promised. The documented tests allow at most one ULP against std on the host.

`eval_expr_integer` returns i64 after checking
`abs(value-round(value)) <= 1e-12*abs(value)`; zero is exact. Out-of-range rounded
values are rejected before casting, including the exclusive upper bound that
`i64::MAX as f64` rounds to. Negative integers are allowed; grid positivity and
application-specific limits belong to the future consumer.

## Quantity boundary and errors

`parse_quantity_expr::<U>(&str) -> Result<(f64, U), ExprError>` has the same
shape and unit types as `parse_quantity`. It first delegates to that existing
parser, preserving successful legacy strings bit-for-bit, including `30 Ry`,
`30Ry`, `1eV`, signed zero and scientific notation. Bare values have no inferred
unit. Legacy nonfinite numeric inputs retain the wrapped unit error.

For expression forms, trim the input and split at its last whitespace character.
The trailing unit token starts with an ASCII letter and contains only ASCII
letters, digits, slash or caret. Digits are allowed to preserve `a0` and pressure
aliases such as `Ha/Bohr^3`; the existing unit enum validates the actual token.
Everything before that whitespace is the expression, which may contain spaces
and parentheses. Thus `2 * 5.43 Angstrom` works. Adjacent units are supported
only for legacy numeric forms; `2/4Ry` is rejected. Unit aliases retain their
existing case-insensitive policy. The return value is in the specified unit;
existing conversions produce atomic units (e.g. `0.5*30 Ry = 15 Ry = 7.5 Ha`).

`ExprError` and `ExprErrorKind` are non-exhaustive. Every error retains complete
original input in `text`, quotes it in Display, and supplies a diagnostic detail.
Kinds are Empty, Syntax (including unknown function and trailing operator),
UnknownName (unbound variable), NonFinite, NotInteger, IntegerOutOfRange, and
Unit(UnitErrorKind). Malformed inputs to the new quantity API may have different
error categories from the legacy API; the legacy API itself is untouched.

## Validation and dependency review

Integration tests print evaluated values and compare arithmetic bit-for-bit;
all eleven functions are checked at three valid points within one ULP, with
additional negative arguments and requested constant expressions. Quantity tests
exercise length, energy, pressure, whitespace, legacy behavior and errors. Integer
tests cover tolerance and range boundaries. A source guard derives forbidden
pi/e digit prefixes from std and scans dft-units/src, avoiding typed constants.
The existing repository hardcoded-unit guard remains unchanged.

Claude supplied the dependency and license review before this ticket. The
pure-Rust library requires no FFI or system installation. Offline host build,
test and clippy gates exercise Apple Silicon compatibility. No independent
claim about current download counts or upstream release activity is made.
The dependency inventory records metadata licenses and existing transitives.
