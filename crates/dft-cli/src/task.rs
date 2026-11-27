//! UX1 single entry (ADR 0053): the `[task]` input table selects the job that a
//! bare `baandi` (or `baandi -i FILE`) runs. Each task maps onto one existing
//! subcommand implementation; the subcommands stay as aliases.
//!
//! Rule: `task.type` is the only selector. Without a `[task]` table the task
//! is `scf`. A `[relax]` table never selects relaxation by itself: a bare
//! `baandi` on an input with `[relax]` but without `task.type = "relax"` or
//! `"cell-relax"` is an input error (nothing is silently ignored).
use serde::Deserialize;
use std::path::{Path, PathBuf};

/// Default input file of the single entry, in the working directory.
pub const DEFAULT_INPUT: &str = "input.toml";

#[derive(Clone, Copy, Debug, Default, Deserialize, PartialEq, Eq)]
#[serde(rename_all = "kebab-case")]
pub enum TaskType {
    /// Self-consistent ground state (`baandi run`).
    #[default]
    Scf,
    /// Fixed-cell ionic relaxation (`baandi relax`, `relax.mode = "ions"`).
    Relax,
    /// Variable-cell relaxation (`baandi relax`, `relax.mode = "cell"` or `"ions_cell"`).
    CellRelax,
    /// SCF followed by a band path (`baandi bands`).
    Bands,
    /// SCF followed by a density of states (`baandi dos`).
    Dos,
    /// SCF followed by a many-band fixed-density solve (`baandi nscf`).
    Nscf,
    /// SCF followed by the optical-limit dielectric tensor (`baandi dielectric`, R2).
    Dielectric,
}
impl TaskType {
    pub fn name(self) -> &'static str {
        match self {
            Self::Scf => "scf",
            Self::Relax => "relax",
            Self::CellRelax => "cell-relax",
            Self::Bands => "bands",
            Self::Dos => "dos",
            Self::Nscf => "nscf",
            Self::Dielectric => "dielectric",
        }
    }
    /// The alias subcommand that runs this task.
    pub fn subcommand(self) -> &'static str {
        match self {
            Self::Scf => "run",
            Self::Relax | Self::CellRelax => "relax",
            Self::Bands => "bands",
            Self::Dos => "dos",
            Self::Nscf => "nscf",
            Self::Dielectric => "dielectric",
        }
    }
}

/// The `[task]` input table.
#[derive(Clone, Debug, Default, Deserialize, PartialEq)]
#[serde(default, deny_unknown_fields)]
pub struct TaskTable {
    #[serde(rename = "type")]
    pub kind: TaskType,
    /// relax / cell-relax: continue from `restart.json` in the output directory.
    pub resume: bool,
    /// nscf: per-k checkpoint directory (relative to the working directory).
    pub checkpoint_dir: Option<PathBuf>,
}

/// A bare invocation resolved to one existing subcommand.
#[derive(Clone, Debug, PartialEq)]
pub struct Resolved {
    pub task: TaskType,
    pub input: PathBuf,
    /// `[output] directory`, default `results`, relative to the working directory.
    pub output_dir: PathBuf,
    pub resume: bool,
    pub checkpoint_dir: Option<PathBuf>,
}

/// Output directory of an input, as every driver resolves it.
pub fn output_dir(value: &toml::Value) -> PathBuf {
    PathBuf::from(
        value
            .get("output")
            .and_then(|o| o.get("directory"))
            .and_then(|d| d.as_str())
            .unwrap_or("results"),
    )
}

/// Resolve the task of an input text. Only the `[task]`, `[relax]` and
/// `[output]` tables are read here; the selected driver parses the full input.
pub fn resolve_text(input: &Path, text: &str) -> Result<Resolved, String> {
    let value: toml::Value =
        toml::from_str(text).map_err(|e| format!("{}: {e}", input.display()))?;
    let table: TaskTable = match value.get("task") {
        Some(t) => t
            .clone()
            .try_into()
            .map_err(|e| format!("{}: [task]: {e}", input.display()))?,
        None => TaskTable::default(),
    };
    let explicit = value.get("task").is_some_and(|t| t.get("type").is_some());
    let relax = value.get("relax");
    let mode = relax
        .and_then(|r| r.get("mode"))
        .and_then(|m| m.as_str())
        .unwrap_or("ions");
    match table.kind {
        TaskType::Relax | TaskType::CellRelax if relax.is_none() => {
            return Err(format!(
                "task.type = \"{}\" requires a [relax] table",
                table.kind.name()
            ));
        }
        TaskType::Relax if mode != "ions" => {
            return Err(format!(
                "task.type = \"relax\" is fixed-cell; relax.mode = \"{mode}\" needs task.type = \"cell-relax\""
            ));
        }
        TaskType::CellRelax if mode == "ions" => {
            return Err(
                "task.type = \"cell-relax\" needs relax.mode = \"cell\" or \"ions_cell\"".into(),
            );
        }
        TaskType::Bands if value.get("bands").is_none() => {
            return Err("task.type = \"bands\" requires a [bands] table".into());
        }
        TaskType::Nscf if value.get("nscf").is_none() => {
            return Err("task.type = \"nscf\" requires an [nscf] table".into());
        }
        TaskType::Dielectric if value.get("response").is_none() => {
            return Err("task.type = \"dielectric\" requires a [response] table".into());
        }
        TaskType::Scf if relax.is_some() => {
            return Err(if explicit {
                "task.type = \"scf\" with a [relax] table: remove [relax] or set task.type = \"relax\" / \"cell-relax\"".into()
            } else {
                "a [relax] table does not select relaxation by itself: set [task] type = \"relax\" (or \"cell-relax\")".into()
            });
        }
        _ => {}
    }
    if table.resume && !matches!(table.kind, TaskType::Relax | TaskType::CellRelax) {
        return Err("task.resume applies to relax and cell-relax only".into());
    }
    if table.checkpoint_dir.is_some() && table.kind != TaskType::Nscf {
        return Err("task.checkpoint_dir applies to nscf only".into());
    }
    Ok(Resolved {
        task: table.kind,
        input: input.to_owned(),
        output_dir: output_dir(&value),
        resume: table.resume,
        checkpoint_dir: table.checkpoint_dir,
    })
}

pub fn resolve(input: &Path) -> Result<Resolved, String> {
    let text = std::fs::read_to_string(input).map_err(|e| format!("{}: {e}", input.display()))?;
    resolve_text(input, &text)
}

#[cfg(test)]
mod tests {
    use super::*;
    const BASE: &str = "[system]\nname='x'\n";
    fn resolve(extra: &str) -> Result<Resolved, String> {
        resolve_text(Path::new("input.toml"), &format!("{BASE}{extra}"))
    }

    #[test]
    fn every_task_type_maps_onto_its_alias() {
        let cases = [
            ("", TaskType::Scf, "run"),
            ("[task]\ntype='scf'\n", TaskType::Scf, "run"),
            ("[task]\ntype='relax'\n[relax]\n", TaskType::Relax, "relax"),
            (
                "[task]\ntype='relax'\n[relax]\nmode='ions'\n",
                TaskType::Relax,
                "relax",
            ),
            (
                "[task]\ntype='cell-relax'\n[relax]\nmode='cell'\n",
                TaskType::CellRelax,
                "relax",
            ),
            (
                "[task]\ntype='cell-relax'\n[relax]\nmode='ions_cell'\n",
                TaskType::CellRelax,
                "relax",
            ),
            ("[task]\ntype='bands'\n[bands]\n", TaskType::Bands, "bands"),
            ("[task]\ntype='dos'\n", TaskType::Dos, "dos"),
            ("[task]\ntype='nscf'\n[nscf]\n", TaskType::Nscf, "nscf"),
            (
                "[task]\ntype='dielectric'\n[response]\n",
                TaskType::Dielectric,
                "dielectric",
            ),
        ];
        for (extra, task, alias) in cases {
            let r = resolve(extra).unwrap_or_else(|e| panic!("{extra}: {e}"));
            assert_eq!(r.task, task, "{extra}");
            assert_eq!(r.task.subcommand(), alias);
            assert_eq!(r.output_dir, PathBuf::from("results"));
        }
        let r = resolve("[task]\ntype='relax'\nresume=true\n[relax]\n[output]\ndirectory='out'\n")
            .unwrap();
        assert!(r.resume);
        assert_eq!(r.output_dir, PathBuf::from("out"));
        let r = resolve("[task]\ntype='nscf'\ncheckpoint_dir='ck'\n[nscf]\n").unwrap();
        assert_eq!(r.checkpoint_dir, Some(PathBuf::from("ck")));
    }

    #[test]
    fn inconsistent_task_tables_reject() {
        for (extra, needle) in [
            ("[relax]\n", "does not select relaxation"),
            ("[task]\ntype='scf'\n[relax]\n", "remove [relax]"),
            ("[task]\ntype='relax'\n", "requires a [relax] table"),
            ("[task]\ntype='relax'\n[relax]\nmode='cell'\n", "cell-relax"),
            ("[task]\ntype='cell-relax'\n[relax]\n", "needs relax.mode"),
            ("[task]\ntype='bands'\n", "[bands]"),
            ("[task]\ntype='nscf'\n", "[nscf]"),
            ("[task]\ntype='dielectric'\n", "[response]"),
            ("[task]\ntype='phonon'\n", "unknown variant"),
            (
                "[task]\ntype='scf'\nresume=true\n",
                "relax and cell-relax only",
            ),
            ("[task]\ntype='dos'\ncheckpoint_dir='x'\n", "nscf only"),
            ("[task]\nkind='scf'\n", "unknown field"),
        ] {
            let e = resolve(extra).unwrap_err();
            assert!(e.contains(needle), "{extra}: {e}");
        }
    }
}
