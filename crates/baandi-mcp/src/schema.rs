//! Shared input declarations and discovery schemas. No schema dependency.
use serde::Deserialize;
use serde_json::{Value, json};

#[derive(Debug, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Mode {
    Observe,
    Supervised,
}
macro_rules! input {
    ($name:ident { $($field:ident : $ty:ty => $schema:expr),* $(,)? }) => {
        #[derive(Debug, Deserialize)]
        #[serde(deny_unknown_fields)]
        pub struct $name { $(pub $field: $ty),* }
        impl $name {
            pub fn schema() -> Value {
                json!({"type":"object", "additionalProperties":false,
                    "properties": {$(stringify!($field): $schema),*},
                    "required":[$(stringify!($field)),*]})
            }
        }
    };
}
fn label() -> Value {
    json!({"type":"string","minLength":1,"maxLength":128,"pattern":"^[!-~]+$"})
}
fn integer() -> Value {
    json!({"type":"integer","minimum":0,"maximum":u64::MAX})
}
input!(Create { input_id: String => label(), calc_id: String => label(), mode: Mode => json!({"type":"string","enum":["observe","supervised"]}) });
input!(Calculation { calc_id: String => label() });
input!(Run { calc_id: String => label(), steps: usize => json!({"type":"integer","minimum":1,"maximum":100}) });
input!(Events { calc_id: String => label(), limit: usize => json!({"type":"integer","minimum":1,"maximum":256}) });
input!(Control {
    calc_id: String => label(), request_id: String => label(),
    observed_state_version: u64 => integer(), observed_policy_revision: u64 => integer(),
    max_iteration: usize => integer(), reason: String => json!({"type":"string","maxLength":512})
});
input!(Policy {
    calc_id: String => label(), request_id: String => label(),
    observed_state_version: u64 => integer(), observed_policy_revision: u64 => integer(),
    max_iteration: usize => integer(), reason: String => json!({"type":"string","maxLength":512}),
    mixing_beta: f64 => json!({"type":"number"})
});
impl From<Policy> for (Control, f64) {
    fn from(p: Policy) -> Self {
        (
            Control {
                calc_id: p.calc_id,
                request_id: p.request_id,
                observed_state_version: p.observed_state_version,
                observed_policy_revision: p.observed_policy_revision,
                max_iteration: p.max_iteration,
                reason: p.reason,
            },
            p.mixing_beta,
        )
    }
}
pub fn tools() -> Value {
    let specs = [
        (
            "create_calculation",
            "Consume a host-validated input ID",
            Create::schema(),
        ),
        (
            "step",
            "Advance one complete SCF iteration",
            Calculation::schema(),
        ),
        (
            "run",
            "Advance a bounded number of complete iterations",
            Run::schema(),
        ),
        (
            "status",
            "Read lifecycle and boundary versions",
            Calculation::schema(),
        ),
        (
            "get_snapshot",
            "Read scalar scientific state with units",
            Calculation::schema(),
        ),
        (
            "get_events",
            "Drain bounded telemetry and read the control audit",
            Events::schema(),
        ),
        (
            "propose_policy_change",
            "Propose mixing beta; requires separate host approval",
            Policy::schema(),
        ),
        (
            "pause",
            "Propose pause; requires host approval",
            Control::schema(),
        ),
        (
            "resume",
            "Propose resume; requires host approval",
            Control::schema(),
        ),
        (
            "cancel",
            "Propose cancellation; requires host approval",
            Control::schema(),
        ),
        (
            "export_result",
            "Return terminal scalar result as JSON and CSV, without writing files",
            Calculation::schema(),
        ),
    ];
    json!(specs.into_iter().map(|(name, description, schema)| json!({"name":name,"description":description,"inputSchema":schema})).collect::<Vec<_>>())
}
/// Validate exactly the small JSON Schema vocabulary generated above. Scientific
/// constraints deliberately remain in Runtime, including rejection auditing.
pub fn validate(schema: &Value, value: &Value) -> bool {
    match schema["type"].as_str() {
        Some("object") => {
            let Some(object) = value.as_object() else {
                return false;
            };
            let props = schema["properties"].as_object().unwrap();
            object.len() == props.len()
                && props
                    .iter()
                    .all(|(key, s)| object.get(key).is_some_and(|v| validate(s, v)))
        }
        Some("string") => value.as_str().is_some_and(|s| {
            s.len() >= schema["minLength"].as_u64().unwrap_or(0) as usize
                && s.len() <= schema["maxLength"].as_u64().unwrap_or(u64::MAX) as usize
                && (schema.get("pattern").is_none() || s.bytes().all(|b| b.is_ascii_graphic()))
                && schema
                    .get("enum")
                    .is_none_or(|e| e.as_array().unwrap().contains(value))
        }),
        Some("integer") => value.as_u64().is_some_and(|n| {
            n >= schema["minimum"].as_u64().unwrap_or(0)
                && n <= schema["maximum"].as_u64().unwrap_or(u64::MAX)
        }),
        Some("number") => value.as_f64().is_some_and(f64::is_finite),
        _ => false,
    }
}
