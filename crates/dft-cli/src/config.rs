//! Strict TOML input; conversion to atomic units happens only at this boundary.
#[path = "preflight.rs"]
pub mod preflight;
use dft_basis::{BasisError, FftGridSpec, PlaneWaveBasis};
use dft_geometry::{Cell, CellError};
use dft_units::{Angstrom, Bohr, Hartree, Rydberg};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ConfigError {
    #[error("{0}")]
    Preflight(Box<preflight::Finding>),
    #[error("invalid TOML input: {0}")]
    Parse(#[from] toml::de::Error),
    #[error("invalid configuration: {0}")]
    Invalid(String),
    #[error("invalid cell: {0}")]
    Cell(#[from] CellError),
    #[error("invalid basis: {0}")]
    Basis(#[from] BasisError),
}

#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Config {
    #[serde(default)]
    pub preflight: preflight::structure::Policy,
    #[serde(skip)]
    pub species_identity: BTreeMap<String, preflight::Identity>,
    #[serde(skip)]
    pub input_expressions: Vec<InputExpression>,
    pub system: System,
    pub pseudopotentials: BTreeMap<String, PathBuf>,
    pub basis: Basis,
    pub xc: Xc,
    pub scf: Option<Scf>,
    pub output: Option<Output>,
    pub bands: Option<Bands>,
    pub dos: Option<Dos>,
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct System {
    pub name: String,
    pub lattice_unit: LatticeUnit,
    /// Each row is one lattice vector, transposed into internal columns.
    pub lattice_vectors: [[f64; 3]; 3],
    pub atoms: Vec<Atom>,
}
#[derive(Debug, Deserialize)]
#[serde(rename_all = "lowercase")]
pub enum LatticeUnit {
    Angstrom,
    Bohr,
}
impl LatticeUnit {
    fn to_bohr(&self, value: f64) -> f64 {
        match self {
            Self::Angstrom => Bohr::from(Angstrom(value)).0,
            Self::Bohr => value,
        }
    }
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Atom {
    pub species: String,
    pub fractional: [f64; 3],
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Basis {
    pub ecutwfc_ry: f64,
    pub kpoints: Option<Vec<Vec<f64>>>,
    pub kgrid: Option<[usize; 3]>,
    pub kshift: Option<[f64; 3]>,
    pub ktime_reversal: Option<bool>,
    pub fft_grid: Option<[usize; 3]>,
}
impl Basis {
    pub fn legacy_gamma(&self) -> bool {
        self.kgrid.is_none()
            && self.ktime_reversal.is_none()
            && self.kpoints.as_deref() == Some(&[vec![0.0, 0.0, 0.0]])
    }
    pub fn time_reversal(&self) -> bool {
        self.ktime_reversal.unwrap_or(self.kgrid.is_some())
    }
    pub fn resolve_kpoints(&self) -> Result<Vec<dft_geometry::kpoints::KPoint>, ConfigError> {
        use dft_geometry::kpoints::{KPoint, monkhorst_pack, reduce_time_reversal, validate};
        let invalid = |s: String| ConfigError::Invalid(format!("basis k points: {s}"));
        let gamma = self.kpoints.as_deref() == Some(&[vec![0.0, 0.0, 0.0]]);
        let points = if let Some(grid) = self.kgrid {
            if self.kpoints.is_some() && !gamma {
                return Err(invalid(
                    "cannot specify both kgrid and kpoints (except legacy Gamma placeholder)"
                        .into(),
                ));
            }
            monkhorst_pack(grid, self.kshift.unwrap_or([0.0; 3]))
                .map_err(|e| invalid(e.to_string()))?
        } else {
            if self.kshift.is_some() {
                return Err(invalid("kshift requires kgrid".into()));
            }
            if gamma {
                vec![KPoint {
                    fractional: [0.0; 3],
                    weight: 1.0,
                }]
            } else {
                self.kpoints.as_ref().ok_or_else(|| invalid("supply kgrid or kpoints".into()))?.iter().map(|row| {
                    if row.len() != 4 {
                        return Err(invalid("explicit kpoints rows require [kx, ky, kz, weight]; only [[0,0,0]] is legacy Gamma".into()));
                    }
                    Ok(KPoint { fractional: [row[0], row[1], row[2]], weight: row[3] })
                }).collect::<Result<Vec<_>, _>>()?
            }
        };
        if points.is_empty() {
            return Err(invalid("explicit kpoints list must not be empty".into()));
        }
        validate(&points).map_err(|e| invalid(e.to_string()))?;
        if self.time_reversal() {
            reduce_time_reversal(&points).map_err(|e| invalid(e.to_string()))
        } else {
            Ok(points)
        }
    }
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Xc {
    pub functional: Functional,
    #[serde(default)]
    pub allow_functional_override: bool,
}
#[derive(Debug, Deserialize)]
pub enum Functional {
    #[serde(rename = "lda-pz")]
    LdaPz,
    #[serde(rename = "pbe")]
    Pbe,
}
impl Functional {
    pub fn engine(&self) -> dft_engine::Functional {
        match self {
            Self::LdaPz => dft_engine::Functional::LdaPz,
            Self::Pbe => dft_engine::Functional::Pbe,
        }
    }
}
/// Optional SCF controls; omitted fields use engine defaults.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Scf {
    pub threads: Option<usize>,
    pub max_iterations: Option<i64>,
    pub mixing: Option<String>,
    pub energy_tolerance_ha: Option<f64>,
    pub density_residual_tolerance: Option<f64>,
    pub mixing_beta: Option<f64>,
    pub mixing_history: Option<usize>,
    pub mixing_q0_inverse_bohr: Option<f64>,
    pub nbands: Option<usize>,
    pub extra_bands: Option<usize>,
    pub eigen_tolerance: Option<f64>,
    pub seed: Option<u64>,
}
impl Scf {
    pub fn mixing_strategy(&self) -> Result<dft_engine::Mixing, ConfigError> {
        let invalid = |message: &str| ConfigError::Invalid(message.into());
        let beta = self.mixing_beta.unwrap_or(0.3);
        let history = self.mixing_history.unwrap_or(6);
        if !beta.is_finite() || beta <= 0.0 {
            return Err(invalid("mixing_beta must be finite and positive"));
        }
        if history == 0 || history.checked_add(1).is_none() {
            return Err(invalid("mixing_history must be positive and representable"));
        }
        let method = self.mixing.as_deref().unwrap_or("pulay");
        if method != "kerker-pulay" && self.mixing_q0_inverse_bohr.is_some() {
            return Err(invalid(
                "mixing_q0_inverse_bohr requires mixing = 'kerker-pulay'",
            ));
        }
        match method {
            "linear" if self.mixing_history.is_none() => Ok(dft_engine::Mixing::Linear { beta }),
            "linear" => Err(invalid("mixing_history does not apply to linear mixing")),
            "pulay" => Ok(dft_engine::Mixing::Pulay { beta, history }),
            "kerker-pulay" => {
                let q0 = self.mixing_q0_inverse_bohr.ok_or_else(|| {
                    invalid("kerker-pulay requires mixing_q0_inverse_bohr in Bohr^-1")
                })?;
                if !q0.is_finite() || q0 <= 0.0 || !(q0 * q0).is_finite() || q0 * q0 == 0.0 {
                    return Err(invalid(
                        "mixing_q0_inverse_bohr must be positive with finite nonzero square",
                    ));
                }
                Ok(dft_engine::Mixing::KerkerPulay {
                    beta,
                    history,
                    q0_inverse_bohr: q0,
                })
            }
            _ => Err(invalid("mixing must be linear, pulay or kerker-pulay")),
        }
    }
}
/// Output paths resolve against the process working directory.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Output {
    pub directory: Option<String>,
    pub volumetric: Option<Volumetric>,
}

/// Optional post-SCF scalar exports; one-based inclusive band selection.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Volumetric {
    pub quantities: Vec<String>,
    pub formats: Vec<String>,
    pub kpoint: Option<usize>,
    pub bands: Option<[usize; 2]>,
}
impl Volumetric {
    pub fn validate(&self) -> Result<(), ConfigError> {
        let fail = |s: &str| ConfigError::Invalid(format!("output.volumetric: {s}"));
        for (values, allowed) in [
            (
                &self.quantities,
                &[
                    "density",
                    "elf",
                    "hartree",
                    "local_potential",
                    "partial_charge",
                ][..],
            ),
            (&self.formats, &["cube", "xsf", "vasp", "vti"][..]),
        ] {
            if values.is_empty() {
                return Err(fail("quantities and formats must be nonempty"));
            }
            for (i, value) in values.iter().enumerate() {
                if !allowed.contains(&value.as_str()) || values[..i].contains(value) {
                    return Err(fail("unknown or duplicate quantity/format"));
                }
            }
        }
        if self.kpoint == Some(0) || self.bands.is_some_and(|[a, b]| a == 0 || a > b) {
            return Err(fail(
                "kpoint and inclusive bands must be positive and ordered",
            ));
        }
        if (self.kpoint.is_some() || self.bands.is_some())
            && !self.quantities.iter().any(|q| q == "partial_charge")
        {
            return Err(fail("band/k selection requires partial_charge"));
        }
        Ok(())
    }
}

/// Original string and resolved boundary value (units follow the documented key).
#[derive(Debug, Serialize)]
pub struct InputExpression {
    pub key_path: String,
    pub text: String,
    pub value: f64,
}

#[derive(Clone, Copy)]
enum Scalar {
    Float,
    Integer,
    Energy(dft_units::EnergyUnit),
    Length,
}

// Traverse only explicitly selected numeric fields; serde still checks shapes/types.
fn resolve(
    value: &mut toml::Value,
    path: &str,
    kind: Scalar,
    records: &mut Vec<InputExpression>,
) -> Result<(), ConfigError> {
    if let Some(array) = value.as_array_mut() {
        for (i, item) in array.iter_mut().enumerate() {
            resolve(item, &format!("{path}[{i}]"), kind, records)?;
        }
    } else if let Some(text) = value.as_str() {
        let fail = |e| ConfigError::Invalid(format!("{path}: {e}"));
        let (resolved, replacement) = match kind {
            Scalar::Integer => {
                let n = dft_units::eval_expr_integer(text).map_err(fail)?;
                if n < 0 && path != "scf.max_iterations" {
                    return Err(ConfigError::Invalid(format!(
                        "{path}: expression {text:?} must be nonnegative"
                    )));
                }
                (n as f64, toml::Value::Integer(n))
            }
            _ => {
                let n = match kind {
                    Scalar::Float => dft_units::eval_expr(text).map_err(fail)?,
                    Scalar::Energy(target) => {
                        let (n, unit) =
                            dft_units::parse_quantity_expr::<dft_units::EnergyUnit>(text)
                                .map_err(fail)?;
                        dft_units::from_hartree(dft_units::to_hartree(n, unit), target)
                    }
                    Scalar::Length => {
                        let (n, unit) =
                            dft_units::parse_quantity_expr::<dft_units::LengthUnit>(text)
                                .map_err(fail)?;
                        dft_units::to_bohr(n, unit)
                    }
                    Scalar::Integer => unreachable!(),
                };
                if !n.is_finite() {
                    return Err(ConfigError::Invalid(format!(
                        "{path}: expression {text:?} is not finite after conversion"
                    )));
                }
                (n, toml::Value::Float(n))
            }
        };
        records.push(InputExpression {
            key_path: path.into(),
            text: text.into(),
            value: resolved,
        });
        *value = replacement;
    } else if matches!(kind, Scalar::Energy(_) | Scalar::Length) {
        return Err(ConfigError::Invalid(format!(
            "{path} requires a unit string: {value}"
        )));
    }
    Ok(())
}

fn normalize(
    value: &mut toml::Value,
    records: &mut Vec<InputExpression>,
) -> Result<(), ConfigError> {
    for (section, floats, integers) in [
        ("system", &["lattice_vectors"][..], &[][..]),
        (
            "basis",
            &["ecutwfc_ry", "kpoints", "kshift"][..],
            &["kgrid", "fft_grid"][..],
        ),
        (
            "scf",
            &[
                "energy_tolerance_ha",
                "density_residual_tolerance",
                "mixing_beta",
                "mixing_q0_inverse_bohr",
                "eigen_tolerance",
            ][..],
            &[
                "max_iterations",
                "nbands",
                "extra_bands",
                "seed",
                "mixing_history",
            ][..],
        ),
    ] {
        if let Some(table) = value.get_mut(section).and_then(toml::Value::as_table_mut) {
            for (keys, kind) in [(floats, Scalar::Float), (integers, Scalar::Integer)] {
                for key in keys {
                    if let Some(v) = table.get_mut(*key) {
                        resolve(v, &format!("{section}.{key}"), kind, records)?;
                    }
                }
            }
        }
    }
    for (section, alternate, legacy, kind) in [
        (
            "basis",
            "cutoff",
            "ecutwfc_ry",
            Scalar::Energy(dft_units::EnergyUnit::Rydberg),
        ),
        (
            "scf",
            "energy_tolerance",
            "energy_tolerance_ha",
            Scalar::Energy(dft_units::EnergyUnit::Hartree),
        ),
        (
            "scf",
            "eigen_tolerance_energy",
            "eigen_tolerance",
            Scalar::Energy(dft_units::EnergyUnit::Hartree),
        ),
        (
            "system",
            "lattice_vectors_with_units",
            "lattice_vectors",
            Scalar::Length,
        ),
    ] {
        if let Some(table) = value.get_mut(section).and_then(toml::Value::as_table_mut)
            && let Some(mut input) = table.remove(alternate)
        {
            if table.contains_key(legacy)
                || (section == "system" && table.contains_key("lattice_unit"))
            {
                return Err(ConfigError::Invalid(format!(
                    "{section}.{alternate}: conflicts with {legacy} or lattice_unit"
                )));
            }
            resolve(&mut input, &format!("{section}.{alternate}"), kind, records)?;
            table.insert(legacy.into(), input);
            if section == "system" {
                table.insert("lattice_unit".into(), toml::Value::String("bohr".into()));
            }
        }
    }
    if let Some(system) = value.get_mut("system").and_then(toml::Value::as_table_mut) {
        // Construct the cell lazily so legacy missing-field errors remain serde errors.
        let cartesian_cell = if system
            .get("atoms")
            .and_then(toml::Value::as_array)
            .is_some_and(|atoms| atoms.iter().any(|a| a.get("cartesian").is_some()))
        {
            let get = |key: &str| {
                system
                    .get(key)
                    .cloned()
                    .ok_or_else(|| ConfigError::Invalid(format!("system.{key} is required")))
            };
            let rows: [[f64; 3]; 3] = get("lattice_vectors")?.try_into()?;
            let unit: LatticeUnit = get("lattice_unit")?.try_into()?;
            Some((
                Cell::from_user_rows_bohr(rows.map(|r| r.map(|x| unit.to_bohr(x))))?,
                unit,
            ))
        } else {
            None
        };
        if let Some(atoms) = system.get_mut("atoms").and_then(toml::Value::as_array_mut) {
            for (i, atom) in atoms.iter_mut().enumerate() {
                for key in ["fractional", "cartesian"] {
                    if let Some(v) = atom.get_mut(key) {
                        resolve(
                            v,
                            &format!("system.atoms[{i}].{key}"),
                            Scalar::Float,
                            records,
                        )?;
                    }
                }
                if let Some(table) = atom.as_table_mut()
                    && let Some(v) = table.remove("cartesian")
                {
                    if table.contains_key("fractional") {
                        return Err(ConfigError::Invalid(format!(
                            "system.atoms[{i}].cartesian: conflicts with fractional"
                        )));
                    }
                    let xyz: [f64; 3] = v.try_into()?;
                    if !xyz.iter().all(|x| x.is_finite()) {
                        return Err(ConfigError::Invalid(format!(
                            "system.atoms[{i}].cartesian must be finite"
                        )));
                    }
                    let (cell, unit) = cartesian_cell.as_ref().unwrap();
                    let fractional =
                        cell.cartesian_bohr_to_fractional(xyz.map(|x| unit.to_bohr(x)));
                    table.insert(
                        "fractional".into(),
                        toml::Value::Array(
                            fractional.into_iter().map(toml::Value::Float).collect(),
                        ),
                    );
                }
            }
        }
    }
    Ok(())
}

pub struct PreparedInput {
    pub cell: Cell,
    pub wavefunctions: PlaneWaveBasis,
    pub density: PlaneWaveBasis,
    pub grid: FftGridSpec,
}
impl Config {
    pub fn parse(text: &str) -> Result<Self, ConfigError> {
        let config = Self::parse_unvalidated(text)?;
        config.validate()?;
        Ok(config)
    }
    /// Parse for diagnostic reporting before geometry validation. Never run unchecked.
    pub fn parse_unvalidated(text: &str) -> Result<Self, ConfigError> {
        let mut value: toml::Value = toml::from_str(text)?;
        let mut records = Vec::new();
        normalize(&mut value, &mut records)?;
        let identities = preflight::normalize(&mut value)?;
        let mut config: Self = value.try_into()?;
        config.input_expressions = records;
        config.species_identity = identities;
        Ok(config)
    }
    /// Validate active physical settings and cell geometry. SCF/output are type-checked only.
    pub fn validate(&self) -> Result<(), ConfigError> {
        let invalid = |message: &str| ConfigError::Invalid(message.to_owned());
        if !self.basis.ecutwfc_ry.is_finite() || self.basis.ecutwfc_ry <= 0. {
            return Err(invalid("basis.ecutwfc_ry must be finite and positive"));
        }
        self.basis.resolve_kpoints()?;
        for atom in &self.system.atoms {
            if !atom.fractional.iter().all(|x| x.is_finite()) {
                return Err(invalid("atomic fractional coordinates must be finite"));
            }
            if !self.pseudopotentials.contains_key(&atom.species) {
                return Err(ConfigError::Invalid(format!(
                    "species {:?} has no pseudopotential entry",
                    atom.species
                )));
            }
        }
        self.cell()?;
        if let Some(v) = self.output.as_ref().and_then(|o| o.volumetric.as_ref()) {
            v.validate()?;
        }
        Ok(())
    }
    pub fn cell(&self) -> Result<Cell, ConfigError> {
        let rows = self
            .system
            .lattice_vectors
            .map(|row| row.map(|x| self.system.lattice_unit.to_bohr(x)));
        Ok(Cell::from_user_rows_bohr(rows)?)
    }
    /// Both CLI commands prepare the grid, so validate also rejects small overrides.
    pub fn prepare(&self) -> Result<PreparedInput, ConfigError> {
        self.validate()?;
        let cell = self.cell()?;
        let cutoff = Hartree::from(Rydberg(self.basis.ecutwfc_ry));
        let wavefunctions = PlaneWaveBasis::new(&cell, cutoff)?;
        let density = wavefunctions.density_basis(&cell)?;
        let grid = FftGridSpec::new(&wavefunctions, &density, self.basis.fft_grid)?;
        let first = self.basis.resolve_kpoints()?[0].fractional;
        let wavefunctions = PlaneWaveBasis::at_k(&cell, cutoff, first)?;
        Ok(PreparedInput {
            cell,
            wavefunctions,
            density,
            grid,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    const EXAMPLE: &str = include_str!("../../../examples/si_lda.toml");
    #[test]
    fn example_and_units_and_rows() {
        let input = Config::parse(EXAMPLE).unwrap();
        let prepared = input.prepare().unwrap();
        assert_eq!(prepared.wavefunctions.cutoff(), Hartree(15.));
        assert!(!prepared.grid.overridden());
        let rows = "[[2.0, 1.0, 0.0], [0.0, 3.0, 1.0], [1.0, 0.0, 4.0]]";
        let changed = EXAMPLE.replace("angstrom", "bohr").replace(
            "[[0.0, 2.715, 2.715], [2.715, 0.0, 2.715], [2.715, 2.715, 0.0]]",
            rows,
        );
        let c = Config::parse(&changed).unwrap();
        assert_eq!(
            c.cell().unwrap().fractional_to_cartesian_bohr([1., 0., 0.]),
            [2., 1., 0.]
        );
        let changed = EXAMPLE.replace("kpoints =", "fft_grid = [27, 27, 27]\nkpoints =");
        assert!(
            Config::parse(&changed)
                .unwrap()
                .prepare()
                .unwrap()
                .grid
                .overridden()
        );
    }
    #[test]
    fn rejects_invalid_inputs() {
        for (from, to) in [
            ("ecutwfc_ry = 30.0", "ecutwfc_ry = 0.0"),
            ("ecutwfc_ry = 30.0", "ecutwfc_ry = -1.0"),
            ("ecutwfc_ry = 30.0", "ecutwfc_ry = nan"),
            ("ecutwfc_ry = 30.0", "ecutwfc_ry = inf"),
            ("kpoints = [[0.0, 0.0, 0.0]]", "kpoints = [[0.1, 0.0, 0.0]]"),
            ("kpoints = [[0.0, 0.0, 0.0]]", "kpoints = []"),
            (
                "kpoints = [[0.0, 0.0, 0.0]]",
                "kpoints = [[0.0, 0.0, 0.0], [0.0, 0.0, 0.0]]",
            ),
            ("kpoints = [[0.0, 0.0, 0.0]]", "kpoints = [[nan, 0.0, 0.0]]"),
            ("angstrom", "meter"),
            ("lda-pz", "unsupported-xc"),
            ("species = \"Si\"", "species = \"C\""),
            (
                "fractional = [0.25, 0.25, 0.25]",
                "fractional = [0.25, 0.25]",
            ),
            (
                "fractional = [0.25, 0.25, 0.25]",
                "fractional = [inf, 0.25, 0.25]",
            ),
            ("2.715", "nan"),
        ] {
            assert!(EXAMPLE.contains(from));
            assert!(
                Config::parse(&EXAMPLE.replace(from, to)).is_err(),
                "accepted {to}"
            );
        }
        let small = EXAMPLE.replace("kpoints =", "fft_grid = [1, 1, 1]\nkpoints =");
        assert!(matches!(
            Config::parse(&small).unwrap().prepare(),
            Err(ConfigError::Basis(BasisError::GridTooSmall { .. }))
        ));
    }
    #[test]
    fn unknown_fields_rejected_at_every_level() {
        assert!(Config::parse(&format!("typo = 1\n{EXAMPLE}")).is_err());
        for header in ["[system]", "[[system.atoms]]", "[basis]", "[xc]"] {
            assert!(
                Config::parse(&EXAMPLE.replace(header, &format!("{header}\ntypo = 1"))).is_err()
            );
        }
        for section in ["scf", "output"] {
            assert!(Config::parse(&format!("{EXAMPLE}\n[{section}]\ntypo = 1\n")).is_err());
        }
        // Pseudopotentials intentionally permits arbitrary species keys, but values are paths.
        assert!(
            Config::parse(&EXAMPLE.replace("Si = \"/absolute/path/to/Si.UPF\"", "Si = 1")).is_err()
        );
    }
    #[test]
    fn mixing_strategy_parameters_and_expressions() {
        let parse = |keys: &str| {
            Config::parse(&format!("{EXAMPLE}\n[scf]\n{keys}\n"))
                .unwrap()
                .scf
                .unwrap()
        };
        assert!(matches!(
            parse("").mixing_strategy().unwrap(),
            dft_engine::Mixing::Pulay {
                beta: 0.3,
                history: 6
            }
        ));
        let keys = "mixing='kerker-pulay'\nmixing_beta='1/2'\nmixing_history='2*3'\nmixing_q0_inverse_bohr='1/2'";
        assert!(matches!(
            parse(keys).mixing_strategy().unwrap(),
            dft_engine::Mixing::KerkerPulay {
                beta: 0.5,
                history: 6,
                q0_inverse_bohr: 0.5
            }
        ));
        for keys in [
            "mixing='kerker-pulay'",
            "mixing='unknown'",
            "mixing_history=0",
            "mixing='linear'\nmixing_history=2",
            "mixing_q0_inverse_bohr=1.0",
            "mixing_beta=nan",
            "mixing='kerker-pulay'\nmixing_q0_inverse_bohr=0.0",
            "mixing='kerker-pulay'\nmixing_q0_inverse_bohr=nan",
            "mixing='kerker-pulay'\nmixing_q0_inverse_bohr=-1.0",
        ] {
            assert!(parse(keys).mixing_strategy().is_err(), "{keys}");
        }
    }
    #[test]
    fn typed_optional_sections() {
        let optional = "\n[scf]\nmax_iterations = 80\nmixing = 'linear'\nenergy_tolerance_ha = 1e-8\ndensity_residual_tolerance = 1e-6\nmixing_beta = 0.3\n[output]\ndirectory = 'results'\n";
        assert!(Config::parse(&format!("{EXAMPLE}{optional}")).is_ok());
        assert!(Config::parse(&format!("{EXAMPLE}\n[scf]\n[output]\n")).is_ok());
        for (from, to) in [
            ("max_iterations = 80", "max_iterations = '7/2'"),
            ("mixing = 'linear'", "mixing = 1"),
            (
                "energy_tolerance_ha = 1e-8",
                "energy_tolerance_ha = 'small'",
            ),
            (
                "density_residual_tolerance = 1e-6",
                "density_residual_tolerance = false",
            ),
            ("mixing_beta = 0.3", "mixing_beta = '1/0'"),
            ("directory = 'results'", "directory = 1"),
        ] {
            assert!(Config::parse(&format!("{EXAMPLE}{}", optional.replace(from, to))).is_err());
        }
        // Reject all superseded provisional keys, even when their types are valid.
        for section in [
            "[scf]\nmax_iter = 80",
            "[scf]\nenergy_tol_ha = 1e-8",
            "[scf]\ndensity_tol = 1e-6",
            "[output]\nprefix = 'si'",
            "[output]\nverbosity = 1",
        ] {
            assert!(Config::parse(&format!("{EXAMPLE}\n{section}\n")).is_err());
        }
    }
    #[test]
    fn threading_input() {
        for threads in [0, 1, 2, 8] {
            let config =
                Config::parse(&format!("{EXAMPLE}\n[scf]\nthreads = {threads}\n")).unwrap();
            assert_eq!(config.scf.unwrap().threads, Some(threads));
        }
        for invalid in ["-1", "1.5", "true", "'two'", "[2]", "9223372036854775808"] {
            assert!(Config::parse(&format!("{EXAMPLE}\n[scf]\nthreads = {invalid}\n")).is_err());
        }
        assert!(Config::parse(EXAMPLE).unwrap().scf.is_none());
    }
}

/// Optional fixed-density band path, with raw numeric values in the stated units.
#[derive(Debug, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct Bands {
    pub preset: Option<String>,
    pub points: Option<Vec<PathPoint>>,
    pub points_per_segment: Option<usize>,
    pub spacing_inv_bohr: Option<f64>,
    pub nbands: Option<usize>,
    pub eigen_tolerance_ha: Option<f64>,
    pub max_iterations: Option<usize>,
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct PathPoint {
    pub label: String,
    pub fractional: [f64; 3],
    #[serde(default)]
    pub break_before: bool,
}
#[derive(Debug, Deserialize, Default)]
#[serde(deny_unknown_fields)]
pub struct Dos {
    pub dos_kgrid: Option<[usize; 3]>,
    pub sigma_ev: Option<f64>,
    pub spacing_ev: Option<f64>,
    pub window_ev: Option<[f64; 2]>,
    pub nbands: Option<usize>,
    pub eigen_tolerance_ha: Option<f64>,
    pub max_iterations: Option<usize>,
}
