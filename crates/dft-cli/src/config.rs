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
    #[serde(skip)]
    symmetry_cache: std::cell::RefCell<Option<std::sync::Arc<dft_engine::PreparedSymmetry>>>,
    /// ADR 0049: record of the input projection onto the detected group
    /// (symmetry.mode = auto only); the projected geometry replaced the input.
    #[serde(skip)]
    pub symmetrization: Option<dft_geometry::symmetry::Symmetrized>,
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
    #[serde(default)]
    pub compat: Compat,
    #[serde(default)]
    pub pseudo: Pseudo,
    #[serde(default)]
    pub symmetry: Symmetry,
    pub spin: Option<SpinConfig>,
    pub occupations: Option<Occupations>,
    pub scf: Option<Scf>,
    pub output: Option<Output>,
    pub bands: Option<Bands>,
    pub dos: Option<Dos>,
    pub relax: Option<crate::relax::Settings>,
    /// D3-DISPERSION (ADR 0048): Grimme D3 correction evaluated once per
    /// geometry after the SCF; absent = no code path runs.
    pub dispersion: Option<Dispersion>,
    /// RUN-CONFIG: switch values by registry name (`m2_grid = "on"`); see
    /// `crate::run_config` (environment > input > default).
    #[serde(default)]
    pub performance: BTreeMap<String, toml::Value>,
    /// RUN-CONFIG: measurement switches (`phase_profile`, `stage_timers`,
    /// `davidson_stages`, `snapshot_dir`, `phase_wait`, `first_block`,
    /// `dgemm_timers`).
    #[serde(default)]
    pub profile: BTreeMap<String, toml::Value>,
}
#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Compat {
    pub mode: Option<String>,
}
impl Compat {
    pub fn reference(&self) -> bool {
        self.mode.as_deref() == Some("reference")
    }
}
pub const REFERENCE_SCOPE: &str = "compat.mode = \"reference\" is for unpolarized PBE energies and densities; spin-polarized PBE, LDA, analytic forces, stress and relaxation are unsupported";

#[derive(Debug, Default, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum SymmetryMode {
    #[default]
    Off,
    Auto,
}
#[derive(Debug, Default, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum DensityGroup {
    #[default]
    Mesh,
    Crystal,
}
#[derive(Debug, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Symmetry {
    pub mode: SymmetryMode,
    pub density_group: DensityGroup,
    pub tolerance_bohr: f64,
}
impl Default for Symmetry {
    fn default() -> Self {
        Self {
            mode: SymmetryMode::Off,
            density_group: DensityGroup::Mesh,
            tolerance_bohr: dft_geometry::symmetry::DEFAULT_TOLERANCE_BOHR,
        }
    }
}
impl Config {
    pub fn density_operations(
        &self,
    ) -> Result<Vec<dft_geometry::symmetry::SymmetryOperation>, ConfigError> {
        Ok(self.prepared_symmetry()?.map_or_else(Vec::new, |p| {
            p.density_operations(self.symmetry.density_group == DensityGroup::Crystal)
        }))
    }
    pub fn symmetry_map(&self) -> Result<Option<dft_geometry::symmetry::KOrbitMap>, ConfigError> {
        Ok(self.prepared_symmetry()?.map(|p| p.map().clone()))
    }
    /// Shared by validation, preflight, result metadata and engine setup.
    pub fn prepared_symmetry(
        &self,
    ) -> Result<Option<std::sync::Arc<dft_engine::PreparedSymmetry>>, ConfigError> {
        if self.symmetry.density_group == DensityGroup::Crystal
            && self.symmetry.mode != SymmetryMode::Auto
        {
            return Err(ConfigError::Invalid(
                "symmetry.density_group=crystal requires symmetry.mode=auto".into(),
            ));
        }
        if self.symmetry.mode == SymmetryMode::Off {
            return Ok(None);
        }
        let fail = |e: String| {
            ConfigError::Invalid(format!(
                "symmetry.operations: {e}; cannot safely reduce; use symmetry.mode='off' or correct the structure/mesh"
            ))
        };
        let full = self.symmetry_full_mesh()?;
        let atoms: Vec<_> = self
            .system
            .atoms
            .iter()
            .map(|a| dft_geometry::Atom::new(&a.species, a.fractional))
            .collect();
        let cell = self.cell()?;
        let tr = self.basis.time_reversal();
        let mut cache = self.symmetry_cache.borrow_mut();
        if let Some(prepared) = cache.as_ref()
            && prepared.matches_for_grid(
                &cell,
                &atoms,
                &full,
                self.detection_tolerance(),
                tr,
                self.basis.fft_grid,
            )
        {
            return Ok(Some(prepared.clone()));
        }
        let prepared = std::sync::Arc::new(
            dft_engine::PreparedSymmetry::new_for_grid(
                &cell,
                &atoms,
                &full,
                self.detection_tolerance(),
                tr,
                self.basis.fft_grid,
            )
            .map_err(|e| fail(e.to_string()))?,
        );
        *cache = Some(prepared.clone());
        Ok(Some(prepared))
    }
    pub fn symmetry_full_mesh(&self) -> Result<Vec<dft_geometry::kpoints::KPoint>, ConfigError> {
        self.basis.resolve_kpoints()?; // Preserve all input consistency checks.
        if let Some(grid) = self.basis.kgrid {
            dft_geometry::kpoints::monkhorst_pack(grid, self.basis.kshift.unwrap_or([0.; 3]))
                .map_err(|e| ConfigError::Invalid(e.to_string()))
        } else if self.basis.kpoints.as_deref() == Some(&[vec![0.0, 0.0, 0.0]]) {
            self.basis.resolve_kpoints()
        } else {
            Err(ConfigError::Invalid("symmetry.mesh: auto requires kgrid or legacy Gamma; use a full Monkhorst-Pack mesh".into()))
        }
    }
    pub fn resolved_points(&self) -> Result<Vec<dft_geometry::kpoints::KPoint>, ConfigError> {
        Ok(if let Some(map) = self.symmetry_map()? {
            map.points
        } else {
            self.basis.resolve_kpoints()?
        })
    }
    pub fn symmetry_json(&self) -> Result<Option<serde_json::Value>, ConfigError> {
        let density_ops = self.density_operations()?;
        Ok(self.symmetry_map()?.map(|m| serde_json::json!({
            "density_group": if self.symmetry.density_group == DensityGroup::Crystal { "crystal" } else { "mesh" },
            "density_operations_count": density_ops.len(),
            "density_operations": density_ops.iter().map(|o| serde_json::json!({"rotation":o.rotation,"translation":o.translation,"atom_permutation":o.atom_permutation,"antiunitary":o.antiunitary})).collect::<Vec<_>>(),
            "mode":"auto", "operations_count":m.operations.iter().filter(|o|!o.antiunitary).count(),
            "irreducible_k_count":m.points.len(),"full_k_count":m.entries.len(),"tolerance_bohr":self.symmetry.tolerance_bohr,"detection_tolerance_bohr":self.detection_tolerance(),
            "crystal_operations_count": self.symmetrization.as_ref().map(|r| r.operations_count),
            "detected_operations_count": self.symmetrization.as_ref().map(|r| r.detected_operations_count),
            "dropped_operations": self.symmetrization.as_ref().filter(|r| r.dropped_for_grid > 0).map(|r| serde_json::json!({
                "count": r.dropped_for_grid,
                "reason": format!("fractional translations/rotations do not permute the explicit basis.fft_grid {:?}; the largest grid-preserving subgroup is used (ADR 0049, GS1C)", self.basis.fft_grid.unwrap_or_default()),
            })),
            "input_symmetrization": self.symmetrization.as_ref().map(|r| serde_json::json!({
                "policy": "ADR 0049: lattice and positions projected onto the group detected at tolerance_bohr; refused above it",
                "positions_changed": r.positions_changed, "lattice_changed": r.lattice_changed,
                "translations_snapped": r.translations_snapped,
                "max_atom_displacement_bohr": r.max_atom_displacement_bohr,
                "max_lattice_change_bohr": r.max_lattice_change_bohr,
            })),
            "phase_convention":m.phase_convention,"magnetic_policy":"identity spin only; spin swap and magnetic antiunitary operations disabled without model-invariance proof",
            "points":m.points.iter().enumerate().map(|(i,p)|serde_json::json!({"fractional":p.fractional,"weight":p.weight,"multiplicity":m.entries.iter().filter(|e|e.representative_id==i).count()})).collect::<Vec<_>>(),
            "operations":m.operations.iter().map(|o|serde_json::json!({"rotation":o.rotation,"translation":o.translation,"atom_permutation":o.atom_permutation,"spin_action":"identity","antiunitary":o.antiunitary,"model_invariance_proof":"nc_scalar_periodic_identity_spin_no_external_terms_v1"})).collect::<Vec<_>>(),
            "orbit_map":m.entries.iter().map(|e|serde_json::json!({"full_id":e.full_id,"representative_id":e.representative_id,"operation_id":e.operation_id,"reciprocal_shift":e.reciprocal_shift,"multiplicity":e.multiplicity})).collect::<Vec<_>>()
        })))
    }
}
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct System {
    pub name: String,
    pub lattice_unit: LatticeUnit,
    /// Each row is one lattice vector, transposed into internal columns.
    pub lattice_vectors: [[f64; 3]; 3],
    pub atoms: Vec<Atom>,
    /// Optional `[system.provenance]` table written by the workstation-side
    /// structure converter (`baandi-convert`). It never enters the physics;
    /// it is checked for plain values and copied into results.json.
    #[serde(default)]
    pub provenance: Option<BTreeMap<String, toml::Value>>,
}
impl System {
    /// `[system.provenance]` as JSON for results.json; `None` when absent.
    pub fn provenance_json(&self) -> Option<serde_json::Value> {
        self.provenance.as_ref().map(|table| {
            serde_json::Value::Object(
                table
                    .iter()
                    .map(|(k, v)| (k.clone(), provenance_value(v)))
                    .collect(),
            )
        })
    }
    fn validate_provenance(&self) -> Result<(), ConfigError> {
        fn plain(v: &toml::Value) -> bool {
            match v {
                toml::Value::String(_)
                | toml::Value::Integer(_)
                | toml::Value::Float(_)
                | toml::Value::Boolean(_) => true,
                toml::Value::Array(items) => items.iter().all(plain),
                toml::Value::Datetime(_) | toml::Value::Table(_) => false,
            }
        }
        for (key, value) in self.provenance.iter().flatten() {
            if !plain(value) {
                return Err(ConfigError::Invalid(format!(
                    "system.provenance.{key}: only strings, numbers, booleans and arrays of them are allowed (write times as strings)"
                )));
            }
        }
        Ok(())
    }
}
fn provenance_value(v: &toml::Value) -> serde_json::Value {
    match v {
        toml::Value::String(s) => serde_json::Value::from(s.as_str()),
        toml::Value::Integer(i) => serde_json::Value::from(*i),
        toml::Value::Float(f) => serde_json::Value::from(*f),
        toml::Value::Boolean(b) => serde_json::Value::from(*b),
        toml::Value::Array(items) => items.iter().map(provenance_value).collect(),
        // Rejected by validation; kept total for unvalidated diagnostics.
        other => serde_json::Value::from(other.to_string()),
    }
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
/// GT1 (ADR 0050): whether the k sampling of a run used time reversal
/// (k ~ -k; within each channel for collinear spin) and why.
pub fn time_reversal_json(config: &Config) -> serde_json::Value {
    let used = config.basis.time_reversal();
    serde_json::json!({
        "used": used,
        "source": if config.basis.ktime_reversal.is_some() { "basis.ktime_reversal" } else if config.basis.kgrid.is_some() { "kgrid default" } else { "explicit-list default" },
        "scope": if config.spin.as_ref().is_some_and(|s| s.mode != SpinMode::Unpolarized) { "within each collinear spin channel; no spin swap" } else { "spin-unpolarized" },
        "policy": "ADR 0050",
    })
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
    pub pw92_prefactor: Option<String>,
    pub functional: Functional,
    #[serde(default)]
    pub allow_functional_override: bool,
}
impl Xc {
    pub fn prefactor(&self) -> Result<dft_engine::Pw92Prefactor, ConfigError> {
        match self.pw92_prefactor.as_deref().unwrap_or("modified") {
            "modified" => Ok(dft_engine::Pw92Prefactor::Modified),
            "original" => Ok(dft_engine::Pw92Prefactor::Original),
            _ => Err(ConfigError::Invalid(
                "xc.pw92_prefactor must be modified or original".into(),
            )),
        }
    }
}
#[derive(Debug, Deserialize)]
pub enum Functional {
    #[serde(rename = "lda-pz")]
    LdaPz,
    #[serde(rename = "lda-pw")]
    LdaPw,
    #[serde(rename = "pbe")]
    Pbe,
}
impl Functional {
    pub fn engine(&self) -> dft_engine::Functional {
        match self {
            Self::LdaPz => dft_engine::Functional::LdaPz,
            Self::LdaPw => dft_engine::Functional::LdaPw,
            Self::Pbe => dft_engine::Functional::Pbe,
        }
    }
}
/// Fermi-Dirac band count when neither `scf.nbands` nor `occupations.band_cap` is given.
pub const DEFAULT_FD_BAND_CAP: usize = 32;
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Occupations {
    pub mode: String,
    pub smearing_ha: Option<f64>,
    pub smearing_ev: Option<f64>,
    pub smearing_kelvin: Option<f64>,
    pub band_cap: Option<usize>,
}
impl Config {
    pub fn finite_temperature(&self) -> Result<Option<f64>, ConfigError> {
        let Some(o) = &self.occupations else {
            return Ok(None);
        };
        let values = [
            o.smearing_ha,
            o.smearing_ev
                .map(|v| dft_units::to_hartree(v, dft_units::EnergyUnit::ElectronVolt)),
            o.smearing_kelvin
                .map(|v| dft_units::to_hartree(v * 1.380649e-23, dft_units::EnergyUnit::Joule)),
        ];
        let count = values.iter().flatten().count();
        match o.mode.as_str() {
            "fixed" if count == 0 && o.band_cap.is_none() => Ok(None),
            "fermi_dirac" if count == 1 => {
                let theta = values.into_iter().flatten().next().unwrap();
                if !theta.is_finite() || theta <= 0.0 || o.band_cap == Some(0) {
                    return Err(ConfigError::Invalid("occupations.smearing: positive finite width and band cap required; use fixed mode for zero temperature".into()));
                }
                Ok(Some(theta))
            }
            "fixed" => Err(ConfigError::Invalid("occupations: fixed mode cannot carry smearing or band cap; remove them".into())),
            "fermi_dirac" => Err(ConfigError::Invalid("occupations.smearing: supply exactly one of smearing_ha, smearing_ev, smearing_kelvin".into())),
            _ => Err(ConfigError::Invalid("occupations.mode: unsupported smearing scheme; select fixed or fermi_dirac".into())),
        }
    }
    /// Band count of a Fermi-Dirac run (per spin channel and k point):
    /// `scf.nbands` when given, otherwise `occupations.band_cap` (default
    /// `DEFAULT_FD_BAND_CAP`). An explicit `band_cap` stays the documented cap:
    /// `scf.nbands` above it is rejected instead of silently replaced.
    pub fn fd_band_count(&self) -> Result<usize, ConfigError> {
        let cap = self.occupations.as_ref().and_then(|o| o.band_cap);
        match (self.scf.as_ref().and_then(|s| s.nbands), cap) {
            (Some(0), _) => Err(ConfigError::Invalid(
                "scf.nbands: Fermi-Dirac occupations need a positive band count".into(),
            )),
            (Some(n), Some(c)) if n > c => Err(ConfigError::Invalid(format!(
                "scf.nbands = {n} exceeds occupations.band_cap = {c}; raise band_cap to at least {n} or lower scf.nbands"
            ))),
            (Some(n), _) => Ok(n),
            (None, c) => Ok(c.unwrap_or(DEFAULT_FD_BAND_CAP)),
        }
    }
    /// Reference FD band count: automatic growth from the engine's start rule
    /// up to `band_cap` (`fd_band_count`) unless `scf.nbands` is explicit,
    /// which keeps the fixed count `scf.nbands`.
    pub fn fd_band_growth(
        &self,
        finite: bool,
    ) -> Option<dft_engine::refscf::electronic::BandGrowth> {
        (finite && self.scf.as_ref().and_then(|s| s.nbands).is_none())
            .then_some(dft_engine::refscf::electronic::BandGrowth::Automatic)
    }
    pub fn finite_constraint(&self) -> dft_physics::occupations::Constraint {
        use dft_physics::occupations::Constraint;
        self.spin
            .as_ref()
            .and_then(|s| {
                s.total_moment
                    .or_else(|| s.n_up.zip(s.n_down).map(|(u, d)| u as f64 - d as f64))
            })
            .map_or(Constraint::Total, Constraint::FixedMoment)
    }
}
/// Additive fixed-sector input; magnetic seeds are provenance only.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct SpinConfig {
    pub mode: SpinMode,
    pub n_up: Option<usize>,
    pub n_down: Option<usize>,
    pub total_moment: Option<f64>,
    pub axis: Option<String>,
    pub axis_cartesian: Option<[f64; 3]>,
    pub initial_moment_seed: Option<f64>,
}
#[derive(Debug, Deserialize, PartialEq)]
#[serde(rename_all = "lowercase")]
pub enum SpinMode {
    Unpolarized,
    Collinear,
}
/// ADR 0043: the SCF path of one run and why.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ScfPath {
    pub standard: Option<dft_engine::StandardProfile>,
    pub verification: dft_engine::Verification,
    /// `applied` (the single standard path) or `reference_strict`
    /// (`BAANDI_REFERENCE_STRICT=on`). Since ONEPATH-4 there is no strict
    /// fallback (`not_applicable` is no longer produced).
    pub status: &'static str,
    pub reason: Option<&'static str>,
    /// Where the standard threshold came from (None off the standard path):
    /// `user` (`scf.standard_threshold_ha`), `stress_default`
    /// (`STRESS_DEFAULT_THRESHOLD_HA`, stress or variable-cell input) or
    /// `default` (`StandardProfile::default()`).
    pub threshold_source: Option<&'static str>,
}
/// Default standard threshold, Ha, of inputs that ask for stress
/// (`output.stress`, `relax.mode = "cell"` or `"ions_cell"`) without
/// `scf.standard_threshold_ha` (ONEPATH-2, orchestrator decision 2026-11-16
/// under the adoption rule "the default must meet the accuracy class"): the
/// default 5e-9 Ha missed the 1e-7 Ha/Bohr^3 stress class on a synthetic
/// nonlocal fixture (1.3e-7), 5e-11 Ha met it on both fixtures (margin 1.6x).
/// A rule by observable, not by material; provisional until the dev08 list
/// `docs/performance/requests/onepath2-stress.md` confirms it.
pub const STRESS_DEFAULT_THRESHOLD_HA: f64 = dft_engine::StandardProfile::STRESS_THRESHOLD_HA;
/// Reason recorded when the reference switch forces the strict path.
pub const REFERENCE_STRICT_REASON: &str = crate::driver_select::STRICT;
/// Error for the removed `scf.convergence_profile` key (user decision 2026-11-16).
pub const CONVERGENCE_PROFILE_REMOVED: &str = "scf.convergence_profile was removed (user decision 2026-11-16, ADR 0043): every run uses the single standard SCF path with the fixed-point verification after convergence; set the accuracy with scf.standard_threshold_ha (default 5e-9 Ha) and the verification with scf.verify = \"fixed_point\" (default) or \"off\"; the strict path is a reference tool only (environment BAANDI_REFERENCE_STRICT=on); delete the key";
/// Error for `scf.locking = "auto"` on the standard path (ONEPATH-4).
impl Config {
    /// Path selection (ADR 0043, default flipped 2026-11-16; single path
    /// complete since ONEPATH-4). Every command (`run` in every execution mode,
    /// `relax`, `bands`, `dos`, `checkpoint`, `restart`, `explain`) uses the
    /// standard path with the fixed-point verification. The strict path is
    /// reached only through the reference switch `BAANDI_REFERENCE_STRICT=on`.
    pub fn scf_path(&self) -> Result<ScfPath, ConfigError> {
        self.scf_path_selected(crate::driver_select::preliminary(self))
    }
    /// `scf_path` with the reference switch given explicitly (`true` = the
    /// strict reference path, as under `BAANDI_REFERENCE_STRICT=on`), for
    /// in-process callers that must match a child process with the switch set.
    pub fn scf_path_with(&self, reference_strict: bool) -> Result<ScfPath, ConfigError> {
        let mut selection = crate::driver_select::preliminary(self);
        selection.driver = if reference_strict {
            crate::driver_select::Driver::StrictOracle
        } else {
            crate::driver_select::Driver::Reference
        };
        self.scf_path_selected(selection)
    }
    pub fn scf_path_selected(
        &self,
        selection: crate::driver_select::Selection,
    ) -> Result<ScfPath, ConfigError> {
        selection.validate(self)?;
        let invalid = |m: String| ConfigError::Invalid(m);
        let s = self.scf.as_ref();
        if let Some(s) = s {
            s.reject_removed_keys()?;
        }
        let verify = s.map(Scf::verify).transpose()?.flatten();
        let threshold = s.and_then(|s| s.standard_threshold_ha);
        if threshold.is_some_and(|t| !t.is_finite() || t <= 0.0) {
            return Err(invalid(
                "scf.standard_threshold_ha must be finite and positive".into(),
            ));
        }
        if selection.strict() {
            // Reference tool: the same input runs on the strict path; the
            // standard-only keys are recorded as ignored by the status.
            return Ok(ScfPath {
                standard: None,
                verification: dft_engine::Verification::Off,
                status: "reference_strict",
                reason: selection.reason,
                threshold_source: None,
            });
        }
        self.check_single_path_input()?;
        let verification = if verify == Some(false) {
            dft_engine::Verification::Off
        } else {
            dft_engine::Verification::FixedPoint
        };
        if verification == dft_engine::Verification::Off && self.relax.is_some() {
            return Err(invalid(
                "relax needs forces from a passed fixed-point verification on every ionic step (ADR 0043 item 5, ONEPATH-1); remove scf.verify = \"off\"".into(),
            ));
        }
        // An explicit threshold always wins; otherwise inputs that ask for
        // stress get the stress default, everything else the class default.
        let (threshold_ha, source) = match threshold {
            Some(t) => (t, "user"),
            None if self.wants_stress() => (STRESS_DEFAULT_THRESHOLD_HA, "stress_default"),
            None => (
                dft_engine::StandardProfile::default().threshold_ha,
                "default",
            ),
        };
        Ok(ScfPath {
            standard: Some(dft_engine::StandardProfile { threshold_ha }),
            verification,
            status: "applied",
            reason: None,
            threshold_source: Some(source),
        })
    }
    /// The input asks for stress: `output.stress` or a variable-cell
    /// relaxation (`relax.mode = "cell"` or `"ions_cell"`).
    /// Pulay assessment mode of every stress evaluation of this input.
    pub fn pulay_assessment(&self) -> dft_engine::PulayAssessmentMode {
        if self
            .output
            .as_ref()
            .is_some_and(|o| o.stress_pulay_assessment)
        {
            dft_engine::PulayAssessmentMode::CutoffScan
        } else {
            dft_engine::PulayAssessmentMode::Off
        }
    }
    pub fn wants_stress(&self) -> bool {
        self.output.as_ref().is_some_and(|o| o.stress)
            || self
                .relax
                .as_ref()
                .is_some_and(|r| r.mode != crate::relax::cell::Mode::Ions)
    }
    /// Inputs the single standard path cannot run (ADR 0043 item 5). Stages
    /// ONEPATH-1 to 3 (user decisions 2026-11-16) moved forces, every
    /// relaxation mode, stress, collinear spin and Fermi-Dirac occupations onto
    /// the standard path with the verification; ONEPATH-4 turned the last
    /// input fallback, `scf.locking = "auto"`, into an input error (the
    /// locking eigensolver exists only for the strict residual Davidson).
    fn check_single_path_input(&self) -> Result<(), ConfigError> {
        self.finite_temperature()?;
        Ok(())
    }
    pub fn collinear(&self) -> bool {
        self.spin
            .as_ref()
            .is_some_and(|s| s.mode == SpinMode::Collinear)
    }
    pub fn spin_sector(
        &self,
        electrons: f64,
    ) -> Result<Option<dft_engine::collinear::SpinSector>, ConfigError> {
        let fail = |s: &str| ConfigError::Invalid(s.into());
        let Some(s) = &self.spin else { return Ok(None) };
        if s.mode == SpinMode::Unpolarized {
            if s.n_up.is_some()
                || s.n_down.is_some()
                || s.total_moment.is_some()
                || s.axis.is_some()
                || s.axis_cartesian.is_some()
                || s.initial_moment_seed.is_some()
            {
                return Err(fail(
                    "spin.mode: unpolarized input cannot carry spin-channel controls; remove them or select collinear",
                ));
            }
            return Ok(None);
        }
        if self.xc.functional.engine() == dft_engine::Functional::LdaPz {
            return Err(fail(
                "xc.spin: required polarized capability unavailable for lda-pz; select pbe or lda-pw",
            ));
        }
        // GT1 (ADR 0050): k ~ -k holds within each collinear channel; the
        // former rejection stays available as COLLINEAR_TIME_REVERSAL=off.
        if self.basis.time_reversal()
            && !dft_numerics::switches::on(dft_numerics::switches::Id::COLLINEAR_TIME_REVERSAL)
        {
            return Err(fail(
                "symmetry.operations: reduction is not valid for this model; use the full mesh with basis.ktime_reversal = false",
            ));
        }
        let finite = self.finite_temperature()?.is_some();
        let (n_up, n_down) = if finite {
            if s.n_up.is_some() != s.n_down.is_some()
                || (s.total_moment.is_some() && s.n_up.is_some())
            {
                return Err(fail(
                    "spin.constraint: competing or incomplete constraint; choose total_moment or both channel counts",
                ));
            }
            if let (Some(u), Some(d)) = (s.n_up, s.n_down)
                && u.checked_add(d)
                    .is_none_or(|n| (n as f64 - electrons).abs() > 1e-10)
            {
                return Err(fail(
                    "spin.constraint: channel counts must sum to total electrons; correct n_up/n_down",
                ));
            }
            if let dft_physics::occupations::Constraint::FixedMoment(m) = self.finite_constraint()
                && (!m.is_finite() || m.abs() >= electrons)
            {
                return Err(fail(
                    "occupations.electron_target: finite temperature requires interior spin-sector targets",
                ));
            }
            let n = electrons.round() as usize;
            (n.div_ceil(2), n / 2)
        } else {
            match (s.n_up, s.n_down, s.total_moment) {
                (Some(u), Some(d), None) => (u, d),
                (None, None, Some(m)) => {
                    let up = (electrons + m) / 2.0;
                    let down = (electrons - m) / 2.0;
                    if !m.is_finite()
                        || up < 0.0
                        || down < 0.0
                        || up >= usize::MAX as f64
                        || down >= usize::MAX as f64
                        || (up - up.round()).abs() > 1e-10
                        || (down - down.round()).abs() > 1e-10
                    {
                        return Err(fail(
                            "spin.constraint: fixed occupations require integer channel counts within capacity; correct total_moment",
                        ));
                    }
                    (up.round() as usize, down.round() as usize)
                }
                _ => {
                    return Err(fail(
                        "spin.constraint: select exactly one fixed sector; supply both n_up/n_down or total_moment",
                    ));
                }
            }
        };
        if s.axis.is_some() && s.axis_cartesian.is_some() {
            return Err(fail(
                "spin.axis: competing axis definitions; supply axis or axis_cartesian",
            ));
        }
        let axis_cartesian = match s.axis.as_deref() {
            None => s.axis_cartesian.unwrap_or([0.0, 0.0, 1.0]),
            Some("x") => [1.0, 0.0, 0.0],
            Some("y") => [0.0, 1.0, 0.0],
            Some("z") => [0.0, 0.0, 1.0],
            _ => return Err(fail("spin.axis: unsupported label; select x, y or z")),
        };
        let sector = dft_engine::collinear::SpinSector {
            n_up,
            n_down,
            axis_cartesian,
            initial_moment_seed: s.initial_moment_seed,
        };
        sector
            .validate(electrons)
            .map_err(|e| fail(&e.to_string()))?;
        Ok(Some(sector))
    }
}
/// Optional SCF controls; omitted fields use engine defaults.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Scf {
    /// Optional reference selector, accepted as a no-op.
    pub driver: Option<String>,
    /// Reference-driver-only screening: plain (default) or tf.
    pub reference_preconditioning: Option<String>,
    pub start_wavefunctions: Option<String>,
    pub start_density: Option<String>,
    pub gamma_real: Option<String>,
    /// Removed key (user decision 2026-11-16): kept only so that an input that
    /// still names it fails with `CONVERGENCE_PROFILE_REMOVED` instead of a
    /// generic unknown-field error.
    #[serde(rename = "convergence_profile")]
    pub removed_convergence_profile: Option<toml::Value>,
    /// Threshold of the standard stop rule in Ha (default 5e-9 Ha, ADR 0035).
    pub standard_threshold_ha: Option<f64>,
    /// ADR 0043: "fixed_point" (default) or "off".
    pub verify: Option<String>,
    pub threads: Option<usize>,
    pub max_iterations: Option<i64>,
    pub mixing: Option<String>,
    /// Original reference-only spelling before the shared preflight adapter.
    #[serde(skip)]
    pub reference_mixing_input: Option<String>,
    pub mixing_beta: Option<f64>,
    pub mixing_history: Option<usize>,
    pub nbands: Option<usize>,
    pub extra_bands: Option<usize>,
    pub seed: Option<u64>,
    /// Oracle compatibility pin: "off" keeps the full-cutoff cold start;
    /// omission keeps the historical fixed-default ladder. Other values reject.
    pub cutoff_ladder: Option<String>,
}
impl Scf {
    pub fn start_density(&self) -> Result<dft_engine::StartDensity, ConfigError> {
        match self.start_density.as_deref().unwrap_or("uniform") {
            "uniform" if self.start_density.is_none() => Ok(dft_engine::StartDensity::Uniform),
            "atomic" => Ok(dft_engine::StartDensity::Atomic),
            _ => Err(ConfigError::Invalid(
                "unknown scf.start_density value; expected atomic".into(),
            )),
        }
    }
    pub fn start_wavefunctions(&self) -> Result<dft_engine::StartWavefunctions, ConfigError> {
        match self.start_wavefunctions.as_deref().unwrap_or("atomic") {
            "atomic" => Ok(dft_engine::StartWavefunctions::Atomic),
            _ => Err(ConfigError::Invalid(
                "unknown scf.start_wavefunctions value; expected atomic".into(),
            )),
        }
    }
    pub fn gamma_real(&self) -> Result<dft_engine::GammaReal, ConfigError> {
        match self.gamma_real.as_deref().unwrap_or("auto") {
            "off" => Ok(dft_engine::GammaReal::Off),
            "auto" => Ok(dft_engine::GammaReal::Auto),
            "on" => Ok(dft_engine::GammaReal::On),
            _ => Err(ConfigError::Invalid(
                "scf.gamma_real must be off, auto or on".into(),
            )),
        }
    }

    /// Reject removed keys and check the standard-path controls' syntax. The
    /// path itself is chosen by `Config::scf_path`.
    pub fn reject_removed_keys(&self) -> Result<(), ConfigError> {
        if self.removed_convergence_profile.is_some() {
            return Err(ConfigError::Invalid(CONVERGENCE_PROFILE_REMOVED.into()));
        }
        Ok(())
    }
    /// `scf.verify`: Some(true) = fixed_point, Some(false) = off, None = default.
    pub fn verify(&self) -> Result<Option<bool>, ConfigError> {
        match self.verify.as_deref() {
            None => Ok(None),
            Some("fixed_point") => Ok(Some(true)),
            Some("off") => Ok(Some(false)),
            Some(_) => Err(ConfigError::Invalid(
                "scf.verify must be fixed_point or off".into(),
            )),
        }
    }
    /// Mixing strategy for the strict oracle; reference mixing is resolved separately.
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
        match method {
            "pulay" | "broyden" | "plain" => Ok(dft_engine::Mixing::Pulay { beta, history }),
            _ => Err(invalid(
                "unknown scf.mixing value; expected plain, broyden or pulay",
            )),
        }
    }
}
/// `[dispersion]` (ADR 0048). Damping parameters are dimensionless except
/// `a2` (Bohr); absent parameters come from the library's set for the run's
/// functional. `library` is used as given (absolute, or relative to the working
/// directory); without it `BAANDI_DFTD3_LIB`, then the platform default name.
#[derive(Debug, Clone, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Dispersion {
    pub method: String,
    pub s6: Option<f64>,
    pub s8: Option<f64>,
    pub s9: Option<f64>,
    pub a1: Option<f64>,
    pub a2: Option<f64>,
    pub rs6: Option<f64>,
    pub rs8: Option<f64>,
    pub alpha6: Option<f64>,
    pub library: Option<PathBuf>,
}

/// Output paths resolve against the process working directory.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Output {
    #[serde(default)]
    pub forces: bool,
    #[serde(default)]
    pub stress: bool,
    /// PULAY-OPTIN (2026-11-26): run the isotropic cutoff-derivative Pulay
    /// assessment of the stress (three extra SCFs); default off. Applies to
    /// `output.stress` and to the nominal-cutoff frames of cell relaxation.
    #[serde(default)]
    pub stress_pulay_assessment: bool,
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
            &["standard_threshold_ha", "mixing_beta"][..],
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

#[derive(Debug)]
pub struct PreparedInput {
    pub cell: Cell,
    pub wavefunctions: PlaneWaveBasis,
    pub density: PlaneWaveBasis,
    pub grid: FftGridSpec,
}
impl Config {
    pub fn parse(text: &str) -> Result<Self, ConfigError> {
        let mut config = Self::parse_unvalidated(text)?;
        config.symmetrize_input()?;
        config.validate()?;
        Ok(config)
    }
    /// ADR 0049: with symmetry.mode = auto, detect the group at the configured
    /// tolerance and replace the lattice and positions by their projection onto
    /// that group (refused above the tolerance). Changes at roundoff level keep
    /// the input bits. Runs once at parse time so every driver, output and
    /// relaxation start sees one consistent, exactly symmetric geometry.
    /// ADR 0049 (amendment GS1B): tolerance (Bohr) for every group detection
    /// after the input was projected onto its group at `symmetry.tolerance_bohr`.
    /// The projected input carries that group exactly and relaxations keep it to
    /// roundoff, so the tight tolerance re-finds exactly the start group. The
    /// loose input tolerance applied to a relaxed geometry would add operations
    /// once the trajectory passes within it of a higher symmetry, and the
    /// symmetrized stress/forces would then lose the component driving it there.
    /// Without a projection record (configs built without `parse`) the
    /// configured tolerance is used unchanged.
    pub fn detection_tolerance(&self) -> f64 {
        if self.symmetrization.is_some() {
            self.symmetry
                .tolerance_bohr
                .min(dft_geometry::symmetry::SYMMETRIZATION_VERIFY_BOHR)
        } else {
            self.symmetry.tolerance_bohr
        }
    }
    pub fn symmetrize_input(&mut self) -> Result<(), ConfigError> {
        if self.symmetry.mode != SymmetryMode::Auto || self.symmetrization.is_some() {
            return Ok(());
        }
        let cell = self.cell()?;
        let atoms: Vec<_> = self
            .system
            .atoms
            .iter()
            .map(|a| dft_geometry::Atom::new(&a.species, a.fractional))
            .collect();
        let record =
            dft_geometry::symmetry::symmetrize_for_grid(
                &cell,
                &atoms,
                self.symmetry.tolerance_bohr,
                self.basis.fft_grid,
            )
                .map_err(|e| {
                    ConfigError::Invalid(format!(
                        "symmetry.operations: {e}; cannot safely reduce; use symmetry.mode='off', a different symmetry.tolerance_bohr or correct the structure"
                    ))
                })?;
        if record.lattice_changed {
            self.system.lattice_unit = LatticeUnit::Bohr;
            self.system.lattice_vectors =
                std::array::from_fn(|v| std::array::from_fn(|c| record.lattice_columns_bohr[c][v]));
        }
        if record.positions_changed {
            for (atom, p) in self.system.atoms.iter_mut().zip(&record.fractional) {
                atom.fractional = *p;
            }
        }
        self.symmetrization = Some(record);
        Ok(())
    }
    /// Parse for diagnostic reporting before geometry validation. Never run unchecked.
    pub fn parse_unvalidated(text: &str) -> Result<Self, ConfigError> {
        let mut value: toml::Value = toml::from_str(text)?;
        let mut records = Vec::new();
        normalize(&mut value, &mut records)?;
        let identities = preflight::normalize(&mut value)?;
        let explicit_density_group = value
            .get("symmetry")
            .and_then(|s| s.get("density_group"))
            .cloned();
        let mut config: Self = value.try_into()?;
        // Both remaining drivers (reference and strict oracle) accept the
        // reference mixer names; broyden/plain are aliases of pulay here.
        if let Some(s) = &mut config.scf {
            s.reference_mixing_input = s.mixing.clone();
            if matches!(s.mixing.as_deref(), Some("broyden" | "plain")) {
                s.mixing = Some("pulay".into());
            }
        }
        match config.compat.mode.as_deref().unwrap_or("none") {
            "none" => {}
            "reference" => {
                config.validate_reference_scope()?;
                for (key, slot, expected) in [
                    (
                        "xc.pw92_prefactor",
                        &mut config.xc.pw92_prefactor,
                        "original",
                    ),
                    (
                        "pseudo.local_mesh",
                        &mut config.pseudo.local_mesh,
                        "qe_like",
                    ),
                    (
                        "pseudo.radial_even_rule",
                        &mut config.pseudo.radial_even_rule,
                        "simpson_quadratic_end",
                    ),
                ] {
                    if let Some(actual) = slot.as_deref()
                        && actual != expected
                    {
                        return Err(ConfigError::Invalid(format!(
                            "compat.mode = \"reference\" conflicts with {key} = {actual:?}; requires {expected:?}"
                        )));
                    }
                    *slot = Some(expected.into());
                }
                let expected = if config.symmetry.mode == SymmetryMode::Auto {
                    "crystal"
                } else {
                    "mesh"
                };
                if let Some(actual) = explicit_density_group
                    && actual.as_str() != Some(expected)
                {
                    return Err(ConfigError::Invalid(format!(
                        "compat.mode = \"reference\" conflicts with symmetry.density_group = {actual}; requires {expected:?} for this symmetry.mode"
                    )));
                }
                config.symmetry.density_group = if config.symmetry.mode == SymmetryMode::Auto {
                    DensityGroup::Crystal
                } else {
                    DensityGroup::Mesh
                };
            }
            _ => {
                return Err(ConfigError::Invalid(
                    "compat.mode must be none or reference".into(),
                ));
            }
        }
        config.input_expressions = records;
        config.species_identity = identities;
        Ok(config)
    }
    fn validate_reference_scope(&self) -> Result<(), ConfigError> {
        if self.compat.reference()
            && (self.xc.functional.engine() != dft_engine::Functional::Pbe
                || self.collinear()
                || self.relax.is_some()
                || self.output.as_ref().is_some_and(|o| o.forces || o.stress))
        {
            return Err(ConfigError::Invalid(REFERENCE_SCOPE.into()));
        }
        Ok(())
    }
    /// Validate active physical settings and cell geometry. SCF/output are type-checked only.
    pub fn validate(&self) -> Result<(), ConfigError> {
        crate::driver_select::preliminary(self).validate(self)?;
        self.validate_reference_scope()?;
        self.system.validate_provenance()?;
        if let Some(s) = &self.scf {
            s.start_density()?;
            s.start_wavefunctions()?;
            s.mixing_strategy()?;
            if !matches!(s.cutoff_ladder.as_deref(), None | Some("off")) {
                return Err(ConfigError::Invalid(
                    "unknown scf.cutoff_ladder value; expected off".into(),
                ));
            }
            if !matches!(s.driver.as_deref(), None | Some("reference")) {
                return Err(ConfigError::Invalid(
                    "unknown scf.driver value; expected reference".into(),
                ));
            }
            if !matches!(
                s.reference_preconditioning.as_deref(),
                None | Some("plain" | "tf")
            ) {
                return Err(ConfigError::Invalid(
                    "scf.reference_preconditioning must be plain or tf".into(),
                ));
            }
        }
        if self.symmetry.density_group == DensityGroup::Crystal {
            if self.symmetry.mode != SymmetryMode::Auto {
                return Err(ConfigError::Invalid(
                    "symmetry.density_group=crystal requires symmetry.mode=auto".into(),
                ));
            }
            if self.collinear()
                || self.relax.is_some()
                || self.output.as_ref().is_some_and(|o| o.forces || o.stress)
            {
                return Err(ConfigError::Invalid("symmetry.density_group=crystal does not support spin, forces, stress or relaxation".into()));
            }
        }
        self.pseudo.local_mesh()?;
        self.pseudo.radial_even_rule()?;
        if self.xc.prefactor()? == dft_engine::Pw92Prefactor::Original
            && (self.xc.functional.engine() != dft_engine::Functional::Pbe || self.collinear())
        {
            return Err(ConfigError::Invalid(
                "xc.pw92_prefactor=original requires unpolarized PBE".into(),
            ));
        }
        if self.finite_temperature()?.is_some() {
            self.fd_band_count()?;
        }
        if let Some(relax) = &self.relax {
            relax
                .validate(self.system.atoms.len())
                .map_err(|e| ConfigError::Invalid(e.to_string()))?;
            // Fixed-cell relaxation runs the standard path with verification
            // (ONEPATH-1); variable-cell relaxation keeps the strict path and
            // `scf_path` rejects standard-only controls there.
        }
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
        if let Some(s) = &self.scf {
            s.gamma_real()?;
            s.start_wavefunctions()?;
            s.start_density()?;
            s.reject_removed_keys()?;
            s.verify()?;
        }
        self.scf_path()?;
        self.cell()?;
        self.symmetry_map()?;
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
        let mut grid = FftGridSpec::new(&wavefunctions, &density, self.basis.fft_grid)?;
        if self.symmetry_map()?.is_some() {
            let ops = self.density_operations()?;
            if self.basis.fft_grid.is_some()
                && !dft_geometry::symmetry::preserves_fft_grid(grid.dims(), &ops)
            {
                return Err(ConfigError::Invalid("symmetry.fft_grid: explicit grid does not preserve translations/rotations; omit basis.fft_grid for automatic compatible dimensions".into()));
            }
            if self.basis.fft_grid.is_none() {
                let dims = dft_geometry::symmetry::compatible_fft_grid(grid.dims(), &ops)
                    .map_err(|e| ConfigError::Invalid(e.to_string()))?;
                grid = FftGridSpec::automatic_with_dims(&wavefunctions, &density, dims)?;
            }
        }
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
    fn density_group_contract_and_metadata() {
        let base = EXAMPLE.replace(
            "kpoints = [[0.0, 0.0, 0.0]]",
            "kgrid = [4,4,4]\nkshift = [0.5,0.5,0.5]",
        );
        let parse = |tail: &str| Config::parse(&format!("{base}\n{tail}\n")).unwrap();
        assert_eq!(parse("").symmetry.density_group, DensityGroup::Mesh);
        assert!(Config::parse(&format!("{base}\n[symmetry]\ndensity_group='typo'\n")).is_err());
        for tail in [
            "[symmetry]\ndensity_group='crystal'",
            "[symmetry]\nmode='off'\ndensity_group='crystal'",
            "[symmetry]\nmode='auto'\ndensity_group='crystal'\n[output]\nforces=true",
            "[symmetry]\nmode='auto'\ndensity_group='crystal'\n[output]\nstress=true",
            "[symmetry]\nmode='auto'\ndensity_group='crystal'\n[spin]\nmode='collinear'",
        ] {
            assert!(
                Config::parse(&format!("{base}\n{tail}\n"))
                    .unwrap_err()
                    .to_string()
                    .contains("density_group")
            );
        }
        let mesh = parse("[symmetry]\nmode='auto'");
        let crystal = parse("[symmetry]\nmode='auto'\ndensity_group='crystal'");
        assert_eq!(
            mesh.resolved_points().unwrap(),
            crystal.resolved_points().unwrap()
        );
        let m = mesh.symmetry_json().unwrap().unwrap();
        let c = crystal.symmetry_json().unwrap().unwrap();
        assert_eq!(m["operations"], c["operations"]);
        assert_eq!(m["orbit_map"], c["orbit_map"]);
        assert_eq!(m["density_operations_count"], 12);
        assert_eq!(c["density_operations_count"], 48);
        assert_eq!(c["density_operations"].as_array().unwrap().len(), 48);
        assert_eq!(c["density_group"], "crystal");
        assert_eq!(crystal.prepare().unwrap().grid.dims(), [32; 3]);
        let bad = base.replace("[basis]", "[basis]\nfft_grid=[25,25,25]");
        let bad = Config::parse(&format!(
            "{bad}\n[symmetry]\nmode='auto'\ndensity_group='crystal'\n"
        ))
        .unwrap();
        // GS1C: an explicit grid that cannot carry the diamond quarter
        // translations keeps the grid-preserving subgroup instead of rejecting.
        assert_eq!(bad.prepare().unwrap().grid.dims(), [25; 3]);
        let kept = bad.density_operations().unwrap();
        assert_eq!(kept.len(), 24); // Td: the quarter-translation half is dropped.
        assert!(dft_geometry::symmetry::preserves_fft_grid([25; 3], &kept));
        assert_eq!(
            bad.symmetrization.as_ref().unwrap().dropped_for_grid,
            48 - kept.len()
        );
    }

    /// User decision 2026-11-16: the profile key is removed; every value of it
    /// (including the old "strict" and "standard") fails with a message that
    /// names the replacement controls.
    #[test]
    fn removed_convergence_profile_key_is_a_clear_error() {
        for value in ["'strict'", "'standard'", "'typo'", "1"] {
            let err = Config::parse(&format!(
                "{EXAMPLE}\n[scf]\nconvergence_profile = {value}\n"
            ))
            .unwrap_err()
            .to_string();
            for needle in [
                "convergence_profile was removed",
                "scf.standard_threshold_ha",
                "scf.verify",
                "BAANDI_REFERENCE_STRICT",
            ] {
                assert!(err.contains(needle), "{value}: {err}");
            }
        }
    }
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
        let keys = "mixing='pulay'\nmixing_beta='1/2'\nmixing_history='2*3'";
        assert!(matches!(
            parse(keys).mixing_strategy().unwrap(),
            dft_engine::Mixing::Pulay {
                beta: 0.5,
                history: 6
            }
        ));
        for keys in [
            "mixing='kerker-pulay'",
            "mixing='linear'",
            "mixing='unknown'",
            "mixing_q0_inverse_bohr=1.0",
        ] {
            assert!(
                Config::parse(&format!("{EXAMPLE}\n[scf]\n{keys}\n")).is_err(),
                "{keys}"
            );
        }
    }
    #[test]
    fn typed_optional_sections() {
        let optional = "\n[scf]\nmax_iterations = 80\nmixing = 'pulay'\nstandard_threshold_ha = 1e-8\nmixing_beta = 0.3\n[output]\ndirectory = 'results'\n";
        assert!(Config::parse(&format!("{EXAMPLE}{optional}")).is_ok());
        assert!(Config::parse(&format!("{EXAMPLE}\n[scf]\n[output]\n")).is_ok());
        for (from, to) in [
            ("max_iterations = 80", "max_iterations = '7/2'"),
            ("mixing = 'pulay'", "mixing = 1"),
            (
                "standard_threshold_ha = 1e-8",
                "standard_threshold_ha = 'small'",
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

/// Independent research options; omission preserves the existing conventions.
#[derive(Debug, Default, Deserialize)]
#[serde(default, deny_unknown_fields)]
pub struct Pseudo {
    pub local_mesh: Option<String>,
    pub radial_even_rule: Option<String>,
}
impl Pseudo {
    pub fn local_mesh(&self) -> Result<dft_engine::LocalMesh, ConfigError> {
        match self.local_mesh.as_deref().unwrap_or("baandi") {
            "baandi" => Ok(dft_engine::LocalMesh::Baandi),
            "qe_like" => Ok(dft_engine::LocalMesh::QeLike),
            "qe_like_g0" => Ok(dft_engine::LocalMesh::QeLikeG0),
            _ => Err(ConfigError::Invalid(
                "pseudo.local_mesh must be baandi, qe_like or qe_like_g0".into(),
            )),
        }
    }
    pub fn radial_even_rule(&self) -> Result<dft_engine::RadialEvenRule, ConfigError> {
        match self.radial_even_rule.as_deref().unwrap_or("simpson_3_8") {
            "simpson_3_8" => Ok(dft_engine::RadialEvenRule::Simpson3_8),
            "simpson_quadratic_end" => Ok(dft_engine::RadialEvenRule::SimpsonQuadraticEnd),
            _ => Err(ConfigError::Invalid(
                "pseudo.radial_even_rule must be simpson_3_8 or simpson_quadratic_end".into(),
            )),
        }
    }
}
