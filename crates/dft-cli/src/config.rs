//! Strict TOML input; conversion to atomic units happens only at this boundary.
use dft_basis::{BasisError, FftGridSpec, PlaneWaveBasis};
use dft_geometry::{Cell, CellError};
use dft_units::{Angstrom, Bohr, Hartree, Rydberg};
use serde::Deserialize;
use std::collections::BTreeMap;
use std::path::PathBuf;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum ConfigError {
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
    pub system: System,
    pub pseudopotentials: BTreeMap<String, PathBuf>,
    pub basis: Basis,
    pub xc: Xc,
    pub scf: Option<Scf>,
    pub output: Option<Output>,
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
}
#[derive(Debug, Deserialize)]
pub enum Functional {
    #[serde(rename = "lda-pz")]
    LdaPz,
}
/// Optional SCF controls; omitted fields use engine defaults.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Scf {
    pub max_iterations: Option<i64>,
    pub mixing: Option<String>,
    pub energy_tolerance_ha: Option<f64>,
    pub density_residual_tolerance: Option<f64>,
    pub mixing_beta: Option<f64>,
    pub nbands: Option<usize>,
    pub extra_bands: Option<usize>,
    pub eigen_tolerance: Option<f64>,
    pub seed: Option<u64>,
}
/// Output paths resolve against the process working directory.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Output {
    pub directory: Option<String>,
}

// Resolve additive unit forms before deserializing the typed tables.
fn normalize_units(value: &mut toml::Value) -> Result<(), ConfigError> {
    let invalid = |message: String| ConfigError::Invalid(message);
    for (section, alternate, legacy, target) in [
        (
            "basis",
            "cutoff",
            "ecutwfc_ry",
            dft_units::EnergyUnit::Rydberg,
        ),
        (
            "scf",
            "energy_tolerance",
            "energy_tolerance_ha",
            dft_units::EnergyUnit::Hartree,
        ),
        (
            "scf",
            "eigen_tolerance_energy",
            "eigen_tolerance",
            dft_units::EnergyUnit::Hartree,
        ),
    ] {
        if let Some(table) = value.get_mut(section).and_then(toml::Value::as_table_mut)
            && let Some(input) = table.remove(alternate)
        {
            if table.contains_key(legacy) {
                return Err(invalid(format!(
                    "{section}: cannot specify both {alternate} and {legacy}"
                )));
            }
            let text = input.as_str().ok_or_else(|| {
                invalid(format!(
                    "{section}.{alternate} requires a unit string: {input}"
                ))
            })?;
            let ha = dft_units::parse_energy(text).map_err(|e| invalid(e.to_string()))?;
            let resolved = dft_units::from_hartree(ha, target);
            if !resolved.is_finite() {
                return Err(invalid(format!("nonfinite converted value: {text:?}")));
            }
            table.insert(legacy.into(), toml::Value::Float(resolved));
        }
    }
    if let Some(table) = value.get_mut("system").and_then(toml::Value::as_table_mut)
        && let Some(input) = table.remove("lattice_vectors_with_units")
    {
        if table.contains_key("lattice_vectors") || table.contains_key("lattice_unit") {
            return Err(invalid(
                "lattice_vectors_with_units conflicts with lattice_vectors and lattice_unit".into(),
            ));
        }
        let rows: [[String; 3]; 3] = input.try_into()?;
        let mut converted = Vec::new();
        for row in rows {
            let mut values = Vec::new();
            for text in row {
                values.push(toml::Value::Float(
                    dft_units::parse_length(&text).map_err(|e| invalid(e.to_string()))?,
                ));
            }
            converted.push(toml::Value::Array(values));
        }
        table.insert("lattice_vectors".into(), toml::Value::Array(converted));
        table.insert("lattice_unit".into(), toml::Value::String("bohr".into()));
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
        let mut value: toml::Value = toml::from_str(text)?;
        normalize_units(&mut value)?;
        let config: Self = value.try_into()?;
        config.validate()?;
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
        Ok(())
    }
    pub fn cell(&self) -> Result<Cell, ConfigError> {
        let rows = self.system.lattice_vectors.map(|row| {
            row.map(|x| match self.system.lattice_unit {
                LatticeUnit::Angstrom => Bohr::from(Angstrom(x)).0,
                LatticeUnit::Bohr => x,
            })
        });
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
            ("lda-pz", "pbe"),
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
    fn typed_optional_sections() {
        let optional = "\n[scf]\nmax_iterations = 80\nmixing = 'linear'\nenergy_tolerance_ha = 1e-8\ndensity_residual_tolerance = 1e-6\nmixing_beta = 0.3\n[output]\ndirectory = 'results'\n";
        assert!(Config::parse(&format!("{EXAMPLE}{optional}")).is_ok());
        assert!(Config::parse(&format!("{EXAMPLE}\n[scf]\n[output]\n")).is_ok());
        for (from, to) in [
            ("max_iterations = 80", "max_iterations = '80'"),
            ("mixing = 'linear'", "mixing = 1"),
            (
                "energy_tolerance_ha = 1e-8",
                "energy_tolerance_ha = 'small'",
            ),
            (
                "density_residual_tolerance = 1e-6",
                "density_residual_tolerance = false",
            ),
            ("mixing_beta = 0.3", "mixing_beta = '0.3'"),
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
}
