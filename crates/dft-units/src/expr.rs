//! Closed, finite f64 arithmetic at the input boundary; see ADR 0014.
use std::{fmt, str::FromStr};

use exmex::{Express, FlatEx, FloatOpsFactory, MakeOperators, MatchLiteral, Operator};

use crate::{UnitError, parse_quantity};

// exmex's default matcher does not cover all scientific f64 literals.
// Use its library-provided regex matcher; exmex still owns expression parsing.
exmex::literal_matcher_from_pattern!(
    InputNumber,
    r"^(?:[0-9]+(?:\.[0-9]*)?|\.[0-9]+)(?:[eE][+-]?[0-9]+)?"
);

#[derive(Clone, Debug)]
struct InputOperators;
impl MakeOperators<f64> for InputOperators {
    fn make<'a>() -> Vec<Operator<'a, f64>> {
        FloatOpsFactory::<f64>::make()
            .into_iter()
            .filter(|op| {
                matches!(
                    op.repr(),
                    "+" | "-"
                        | "*"
                        | "/"
                        | "^"
                        | "sqrt"
                        | "abs"
                        | "sin"
                        | "cos"
                        | "tan"
                        | "asin"
                        | "acos"
                        | "atan"
                        | "exp"
                        | "ln"
                        | "log10"
                )
            })
            .map(|op| {
                if op.repr() == "+" {
                    Operator::make_bin("+", op.bin().expect("exmex addition is binary"))
                } else {
                    op
                }
            })
            .collect()
    }
}

/// Stable error categories; callers must allow future variants.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum ExprErrorKind {
    Empty,
    Syntax,
    UnknownName,
    NonFinite,
    NotInteger,
    IntegerOutOfRange,
    Unit(crate::UnitErrorKind),
}

/// Every error retains and quotes the complete original input.
/// Syntax includes unknown functions, unsupported operators and malformed input;
/// UnknownName identifies unbound variables in a successfully parsed expression.
#[derive(Clone, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub struct ExprError {
    pub kind: ExprErrorKind,
    pub text: String,
    pub detail: String,
}
impl fmt::Display for ExprError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{:?} in expression {:?}: {}",
            self.kind, self.text, self.detail
        )
    }
}
impl std::error::Error for ExprError {}
fn error(kind: ExprErrorKind, text: &str, detail: impl ToString) -> ExprError {
    ExprError {
        kind,
        text: text.into(),
        detail: detail.to_string(),
    }
}

/// Evaluate + - * / ^, unary minus, parentheses, pi/e, and
/// sqrt/abs/sin/cos/tan/asin/acos/atan/exp/ln/log10. No user variables.
/// Plain f64 literals retain Rust parsing behavior, including signed zero.
/// Only finite results are accepted. Transcendentals follow host std f64 math.
pub fn eval_expr(text: &str) -> Result<f64, ExprError> {
    let s = text.trim();
    if s.is_empty() {
        return Err(error(ExprErrorKind::Empty, text, "empty input"));
    }
    let value = if let Ok(value) = s.parse::<f64>() {
        value
    } else {
        // Character policy only: all grammar, tokenization and evaluation belong to exmex.
        if !s
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || c.is_whitespace() || "+-*/^().".contains(c))
        {
            return Err(error(ExprErrorKind::Syntax, text, "unsupported character"));
        }
        // Avoid algebraic compilation/reassociation; retain the parsed operation order.
        // exmex recognizes ASCII spaces; preserve token boundaries for other whitespace.
        let normalized: String = s
            .chars()
            .map(|c| if c.is_whitespace() { ' ' } else { c })
            .collect();
        let expr = FlatEx::<f64, InputOperators, InputNumber>::parse_wo_compile(&normalized)
            .map_err(|e| error(ExprErrorKind::Syntax, text, e))?;
        let values = expr
            .var_names()
            .iter()
            .map(|name| match name.as_str() {
                "pi" => Ok(std::f64::consts::PI),
                "e" => Ok(std::f64::consts::E),
                _ => Err(error(
                    ExprErrorKind::UnknownName,
                    text,
                    format!("unbound name {name:?}"),
                )),
            })
            .collect::<Result<Vec<_>, _>>()?;
        expr.eval(&values)
            .map_err(|e| error(ExprErrorKind::Syntax, text, e))?
    };
    if !value.is_finite() {
        return Err(error(
            ExprErrorKind::NonFinite,
            text,
            "result is not finite",
        ));
    }
    Ok(value)
}

/// Input provenance for future CLI result serialization (no I/O or serde here).
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedExpression {
    pub text: String,
    pub value: f64,
}
/// Return the exact f64 used by the evaluator together with unmodified input.
pub fn resolve_expr(text: &str) -> Result<ResolvedExpression, ExprError> {
    Ok(ResolvedExpression {
        text: text.into(),
        value: eval_expr(text)?,
    })
}

/// Round to i64 only when |value - round(value)| <= 1e-12 * |value|.
/// Zero is exact; values outside the i64 range are rejected before casting.
pub fn eval_expr_integer(text: &str) -> Result<i64, ExprError> {
    let value = eval_expr(text)?;
    let rounded = value.round();
    if (value - rounded).abs() > 1e-12 * value.abs() {
        return Err(error(
            ExprErrorKind::NotInteger,
            text,
            "not an integer within relative 1e-12",
        ));
    }
    // i64::MAX rounds up to the exclusive upper bound when converted to f64.
    if rounded < i64::MIN as f64 || rounded >= -(i64::MIN as f64) {
        return Err(error(
            ExprErrorKind::IntegerOutOfRange,
            text,
            "outside i64 range",
        ));
    }
    Ok(rounded as i64)
}

/// Parse an expression followed by a whitespace-separated unit, returning its
/// value in that unit. Legacy adjacent numeric units still use parse_quantity.
/// A trailing token starts with an ASCII letter and contains only ASCII letters,
/// digits, '/' or '^' (digits preserve a0 and pressure aliases). No default unit.
pub fn parse_quantity_expr<U: FromStr<Err = UnitError>>(text: &str) -> Result<(f64, U), ExprError> {
    let legacy_error = match parse_quantity::<U>(text) {
        Ok(result) => return Ok(result),
        Err(e) => e,
    };
    if matches!(legacy_error.kind, crate::UnitErrorKind::NonFinite) {
        return Err(error(
            ExprErrorKind::Unit(legacy_error.kind),
            text,
            legacy_error,
        ));
    }
    let s = text.trim();
    if let Some(i) = s.rfind(char::is_whitespace) {
        let token = &s[i..].trim();
        if token.starts_with(|c: char| c.is_ascii_alphabetic())
            && token
                .chars()
                .all(|c| c.is_ascii_alphanumeric() || "/^".contains(c))
        {
            let unit = token
                .parse::<U>()
                .map_err(|e| error(ExprErrorKind::Unit(e.kind), text, e))?;
            let value = eval_expr(s[..i].trim()).map_err(|mut e| {
                e.text = text.into();
                e
            })?;
            return Ok((value, unit));
        }
    }
    Err(error(
        ExprErrorKind::Unit(legacy_error.kind),
        text,
        legacy_error,
    ))
}
