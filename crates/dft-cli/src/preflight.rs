//! Shared CLI preflight. No SCF problem, orbital arrays or density files are opened.
use super::{Config, ConfigError};
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_units::Hartree;
use serde::{Deserialize, Serialize};
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

#[derive(Debug, Default, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Identity {
    pub path: PathBuf,
    pub element: Option<String>,
    #[serde(default)]
    pub allow_element_mismatch: bool,
}

pub fn normalize(value: &mut toml::Value) -> Result<BTreeMap<String, Identity>, ConfigError> {
    let mut identities = BTreeMap::new();
    if let Some(table) = value
        .get_mut("pseudopotentials")
        .and_then(toml::Value::as_table_mut)
    {
        for (label, entry) in table {
            if entry.is_table() {
                let identity: Identity = entry.clone().try_into()?;
                *entry = toml::Value::String(identity.path.to_string_lossy().into_owned());
                identities.insert(label.clone(), identity);
            }
        }
    }
    Ok(identities)
}
fn fail(key: &str, problem: impl std::fmt::Display, fix: &str) -> ConfigError {
    ConfigError::Invalid(format!(
        "{key}: {problem}; cannot safely prepare the requested calculation; {fix}"
    ))
}
fn chemical(s: &str) -> bool {
    "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og".split_whitespace().any(|x| x == s)
}

/// Version of the additive findings contract, independent of the legacy report.
pub const FINDINGS_SCHEMA_VERSION: u32 = 1;

#[derive(Debug, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum Severity {
    Error,
    Warning,
    Info,
}

/// Code and evidence are one typed variant, so they cannot disagree.
#[derive(Debug, Serialize, Deserialize, PartialEq)]
#[serde(tag = "code", content = "evidence", deny_unknown_fields)]
pub enum FindingData {
    #[serde(rename = "PF_INPUT_INVALID")]
    InputInvalid { message: String },
    #[serde(rename = "PF_CELL_VOLUME")]
    CellVolume {
        volume_bohr3: Option<f64>,
        minimum_bohr3: Option<f64>,
        reason: String,
    },
    #[serde(rename = "PF_CELL_SKEW")]
    CellSkew {
        volume_ratio: f64,
        minimum_ratio: f64,
    },
    #[serde(rename = "PF_MIN_DISTANCE")]
    MinDistance {
        atoms: [usize; 2],
        distance_bohr: f64,
        minimum_bohr: f64,
    },
    #[serde(rename = "PF_ATOM_OVERLAP")]
    AtomOverlap {
        atoms: [usize; 2],
        distance_bohr: f64,
        tolerance_bohr: f64,
    },
    #[serde(rename = "PF_SLAB_VACUUM")]
    SlabVacuum {
        axis: usize,
        gap_bohr: f64,
        minimum_bohr: f64,
    },
    #[serde(rename = "PF_KMESH_CELL")]
    KmeshCell {
        lengths_bohr: [f64; 3],
        spacing_inverse_bohr: [f64; 3],
        maximum_inverse_bohr: f64,
    },
    #[serde(rename = "PF_ELECTRON_PARITY")]
    ElectronParity { electrons: usize },
    #[serde(rename = "PF_SPECIES_IDENTITY")]
    SpeciesIdentity {
        species: String,
        declared_element: String,
        upf_element: String,
        upf_file: String,
    },
    #[serde(rename = "PF_CUTOFF_BELOW_SUGGESTED")]
    CutoffBelowSuggested {
        species: String,
        cutoff_ry: f64,
        suggested_cutoff_ry: f64,
    },
    #[serde(rename = "PF_ELEMENT_OVERRIDE")]
    ElementOverride {
        species: String,
        declared_element: String,
        upf_element: String,
        mismatch: bool,
    },
    #[serde(rename = "PF_FUNCTIONAL_OVERRIDE")]
    FunctionalOverride {
        requested_functional: String,
        upf_functional: String,
        species: String,
    },
}

#[derive(Debug, Serialize, Deserialize, PartialEq)]
pub struct Finding {
    pub severity: Severity,
    /// TOML path; dynamic species labels are quoted, including labels with dots.
    pub key_path: String,
    #[serde(flatten)]
    pub data: FindingData,
    pub fix_hint: String,
}

/// Read-only input geometry, with unit conversion only at the CLI boundary.
#[derive(Debug, Serialize)]
pub struct StructureSummary {
    /// Each row is a lattice vector, matching the input file convention.
    pub lattice_vectors_bohr: [[f64; 3]; 3],
    pub atoms: Vec<StructureAtom>,
    pub review: structure::Metrics,
}

#[derive(Debug, Serialize)]
pub struct StructureAtom {
    pub species: String,
    pub fractional: [f64; 3],
}

#[derive(Debug, Serialize)]
pub struct Report {
    pub findings_schema_version: u32,
    pub findings: Vec<Finding>,
    pub structure: StructureSummary,
    pub schema_version: u32,
    pub valid: bool,
    pub warnings: Vec<String>,
    pub overrides: Vec<String>,
    pub summary: String,
    pub electrons: usize,
    pub bands: usize,
    pub plane_waves_per_k: Vec<usize>,
    pub fft_grid: [usize; 3],
    pub estimated_array_bytes: u128,
    pub memory_note: &'static str,
}
impl Report {
    pub fn print(&self) {
        for finding in &self.findings {
            eprintln!("{finding}");
        }
        for warning in &self.warnings {
            eprintln!("warning: {warning}");
        }
        for item in &self.overrides {
            eprintln!("override: {item}");
        }
        println!("{}", self.summary);
        println!(
            "Plane waves per k: {:?}; FFT grid: {:?}; estimated array payload: {} bytes ({})",
            self.plane_waves_per_k, self.fft_grid, self.estimated_array_bytes, self.memory_note
        );
    }
}

pub fn check(config: &Config, input: &Path) -> Result<Report, ConfigError> {
    let (mut findings, metrics) = structure::review(config)?;
    config.validate()?;
    if config.system.atoms.is_empty() {
        return Err(fail(
            "system.atoms",
            "empty system",
            "add at least one atom",
        ));
    }
    if let Some(s) = &config.scf {
        for (key, value) in [
            ("energy_tolerance_ha", s.energy_tolerance_ha),
            ("density_residual_tolerance", s.density_residual_tolerance),
            ("eigen_tolerance", s.eigen_tolerance),
            ("mixing_beta", s.mixing_beta),
        ] {
            if value.is_some_and(|x| !x.is_finite() || x <= 0. || (key == "mixing_beta" && x > 1.))
            {
                return Err(fail(
                    &format!("scf.{key}"),
                    "invalid value",
                    if key == "mixing_beta" {
                        "set a finite value in (0, 1]"
                    } else {
                        "set a positive finite tolerance"
                    },
                ));
            }
        }
        if s.max_iterations.is_some_and(|n| n <= 0) {
            return Err(fail(
                "scf.max_iterations",
                "nonpositive iteration budget",
                "set a positive integer",
            ));
        }
        if s.threads == Some(0) {
            return Err(fail(
                "scf.threads",
                "zero threads",
                "set a positive integer",
            ));
        }
        if s.mixing
            .as_deref()
            .is_some_and(|s| s != "linear" && s != "pulay")
        {
            return Err(fail("scf.mixing", "unknown mixer", "use linear or pulay"));
        }
    }
    let options =
        crate::run::options(config).map_err(|e| fail("scf", e, "correct the SCF controls"))?;
    let mut warnings = Vec::new();
    let mut overrides = Vec::new();
    let mut charges = BTreeMap::new();
    let mut elements = BTreeMap::new();
    for (label, path) in &config.pseudopotentials {
        let segment = if label
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || b == b'_' || b == b'-')
        {
            label.clone()
        } else {
            toml::Value::String(label.clone()).to_string()
        };
        let key = format!("pseudopotentials.{segment}");
        let parent = input
            .parent()
            .filter(|p| !p.as_os_str().is_empty())
            .unwrap_or(Path::new("."));
        let parent = parent
            .canonicalize()
            .map_err(|e| fail("input", e, "use a readable input directory"))?;
        let path = parent.join(path);
        let p = dft_upf::from_file(&path).map_err(|e| {
            fail(
                &key,
                format!("{}: {e}", path.display()),
                "supply a readable supported NC UPF file",
            )
        })?;
        let identity = config.species_identity.get(label);
        let element = identity.and_then(|i| i.element.as_deref()).unwrap_or(label);
        if !chemical(element) {
            return Err(fail(
                &format!("{key}.element"),
                format!("label {label:?} needs an explicit chemical element"),
                "use { path = 'file.UPF', element = 'Si' } with the intended element",
            ));
        }
        if element != p.element().trim() {
            let message = format!(
                "{key}.element: label {label:?}, declared element {element:?}, file {}, UPF header element {:?}",
                path.display(),
                p.element()
            );
            if identity.is_some_and(|i| i.allow_element_mismatch) {
                overrides.push(format!("{message}; allow_element_mismatch = true"));
            } else {
                return Err(Finding::error(
                    &key,
                    FindingData::SpeciesIdentity {
                        species: label.clone(),
                        declared_element: element.into(),
                        upf_element: p.element().trim().into(),
                        upf_file: path.display().to_string(),
                    },
                    "correct the label/element or select the matching UPF; deliberate tests may set allow_element_mismatch = true",
                ));
            }
        } else if identity.is_some_and(|i| i.allow_element_mismatch) {
            overrides.push(format!(
                "{key}.allow_element_mismatch = true (permission enabled; no mismatch)"
            ));
        }
        if identity.is_some_and(|i| i.allow_element_mismatch) {
            findings.push(Finding {
                severity: Severity::Info,
                key_path: format!("pseudopotentials.{}.allow_element_mismatch", toml::Value::String(label.clone())),
                data: FindingData::ElementOverride {
                    species: label.clone(),
                    declared_element: element.to_owned(),
                    upf_element: p.element().trim().to_owned(),
                    mismatch: element != p.element().trim(),
                },
                fix_hint: "Review the intended species and UPF identity; remove the override unless deliberate.".into(),
            });
        }
        if !options.functional.matches_upf(p.functional()) {
            let message = format!(
                "xc.functional: requested {}, {} has UPF functional {:?}",
                options.functional.name(),
                path.display(),
                p.functional()
            );
            if options.allow_functional_override {
                findings.push(Finding {
                    severity: Severity::Info,
                    key_path: "xc.allow_functional_override".into(),
                    data: FindingData::FunctionalOverride {
                        requested_functional: options.functional.name().to_owned(),
                        upf_functional: p.functional().to_owned(),
                        species: label.clone(),
                    },
                    fix_hint: "Review the functional and UPF compatibility; select matching inputs unless the override is deliberate.".into(),
                });
                overrides.push(format!("{message}; xc.allow_functional_override = true"));
            } else {
                return Err(fail(
                    "xc.functional",
                    message,
                    "select a matching UPF/functional or explicitly set xc.allow_functional_override = true",
                ));
            }
        }
        let suggested = p.raw().header.wfc_cutoff;
        if suggested.is_finite() && suggested > 0. && config.basis.ecutwfc_ry < suggested {
            findings.push(Finding {
                severity: Severity::Warning,
                key_path: "basis.ecutwfc_ry".into(),
                data: FindingData::CutoffBelowSuggested {
                    species: label.clone(),
                    cutoff_ry: config.basis.ecutwfc_ry,
                    suggested_cutoff_ry: suggested,
                },
                fix_hint: "Increase cutoff and test convergence; the UPF suggestion is not a convergence certificate.".into(),
            });
            warnings.push(format!("basis.ecutwfc_ry: {} Ry is below {} suggested wfc_cutoff = {suggested} Ry; increase cutoff and test convergence", config.basis.ecutwfc_ry, path.display()));
        }
        elements.insert(label, p.element().trim().to_owned());
        charges.insert(label, p.z_valence());
    }
    let electrons: f64 = config
        .system
        .atoms
        .iter()
        .map(|a| charges[&a.species])
        .sum();
    if !electrons.is_finite()
        || electrons <= 0.
        || electrons >= usize::MAX as f64
        || (electrons - electrons.round()).abs() > 1e-10
    {
        return Err(fail(
            "system.atoms / pseudopotentials",
            "nonpositive or nonintegral valence electron count",
            "use potentials and composition giving a positive integral count",
        ));
    }
    let electrons = electrons.round() as usize;
    if !electrons.is_multiple_of(2) {
        return Err(Finding::error(
            "system.atoms / pseudopotentials",
            FindingData::ElectronParity { electrons },
            &format!(
                "odd electron count {electrons}: spin-unpolarized occupations require pairs; use an even-electron composition"
            ),
        ));
    }
    let default_bands = (electrons / 2)
        .checked_add(options.extra_bands)
        .ok_or_else(|| {
            fail(
                "scf.extra_bands",
                "band count overflow",
                "reduce extra_bands",
            )
        })?;
    let bands = options.bands.unwrap_or(default_bands);
    if bands < electrons / 2 {
        return Err(fail(
            "scf.nbands",
            "too few bands for occupied electron pairs",
            "set nbands >= electrons/2",
        ));
    }
    if !config.basis.legacy_gamma() && bands <= electrons / 2 {
        return Err(fail(
            "scf.nbands / scf.extra_bands",
            "k-point gap certification requires an empty band",
            "set nbands > electrons/2 or extra_bands >= 1",
        ));
    }
    let prepared = config.prepare().map_err(|e| {
        fail(
            "basis.ecutwfc_ry / basis.fft_grid",
            e,
            "use a representable cutoff and a sufficiently large FFT grid",
        )
    })?;
    let points = config.basis.resolve_kpoints()?;
    let mut counts = Vec::new();
    for point in &points {
        let basis = PlaneWaveBasis::at_k(
            &prepared.cell,
            Hartree(options.cutoff_hartree),
            point.fractional,
        )
        .map_err(|e| fail("basis.ecutwfc_ry", e, "correct the cutoff"))?;
        let n = basis.vectors().len();
        if bands > n || (options.davidson_subspace.min(n) <= bands && bands != n) {
            return Err(fail(
                "scf.nbands / basis.ecutwfc_ry",
                format!(
                    "{bands} bands incompatible with {n} plane waves and Davidson subspace {}",
                    options.davidson_subspace
                ),
                "increase cutoff or reduce bands within occupation and subspace limits",
            ));
        }
        FftGridSpec::new(&basis, &prepared.density, Some(prepared.grid.dims()))
            .map_err(|e| fail("basis.fft_grid", e, "increase FFT dimensions"))?;
        counts.push(n);
    }
    let mut composition = BTreeMap::<&str, usize>::new();
    for atom in &config.system.atoms {
        *composition.entry(&atom.species).or_default() += 1;
    }
    let composition = composition
        .iter()
        .map(|(s, n)| format!("{s}:{n}"))
        .collect::<Vec<_>>()
        .join(", ");
    let mut material = BTreeMap::<&str, usize>::new();
    for atom in &config.system.atoms {
        *material
            .entry(elements[&atom.species].as_str())
            .or_default() += 1;
    }
    let material = material
        .iter()
        .map(|(s, n)| format!("{s}:{n}"))
        .collect::<Vec<_>>()
        .join(", ");
    let sampling = if let Some(grid) = config.basis.kgrid {
        format!(
            "{grid:?}, shift {:?}",
            config.basis.kshift.unwrap_or([0.; 3])
        )
    } else if config.basis.legacy_gamma() {
        "Gamma only".into()
    } else {
        format!("explicit {:?}", config.basis.kpoints)
    };
    let summary = format!(
        "Valid: {}. Will compute spin-unpolarized {} SCF for UPF composition {} (labels {}; {} atoms, {electrons} electrons), k mesh {} ({} resolved points), cutoff {} Ry, {bands} bands, {} threads.",
        config.system.name,
        options.functional.name(),
        material,
        composition,
        config.system.atoms.len(),
        sampling,
        points.len(),
        config.basis.ecutwfc_ry,
        options.threads
    );
    let fft_grid = prepared.grid.dims();
    let estimated_array_bytes = counts
        .iter()
        .map(|&n| n as u128 * bands as u128 * 16)
        .sum::<u128>()
        + fft_grid.iter().map(|&n| n as u128).product::<u128>() * 24;
    Ok(Report {
        findings_schema_version: FINDINGS_SCHEMA_VERSION,
        findings,
        structure: StructureSummary {
            review: metrics,
            lattice_vectors_bohr: config
                .system
                .lattice_vectors
                .map(|v| v.map(|x| config.system.lattice_unit.to_bohr(x))),
            atoms: config
                .system
                .atoms
                .iter()
                .map(|atom| StructureAtom {
                    species: atom.species.clone(),
                    fractional: atom.fractional,
                })
                .collect(),
        },
        schema_version: 1,
        valid: true,
        warnings,
        overrides,
        summary,
        electrons,
        bands,
        plane_waves_per_k: counts,
        fft_grid,
        estimated_array_bytes,
        memory_note: "orbital plus one real and one complex grid payload only; excludes scratch, projectors, allocator overhead; not peak RSS",
    })
}

impl std::fmt::Display for Finding {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{}: {}; {}",
            self.key_path,
            serde_json::to_string(&self.data).map_err(|_| std::fmt::Error)?,
            self.fix_hint
        )
    }
}
impl Finding {
    pub fn error(key: &str, data: FindingData, hint: &str) -> ConfigError {
        ConfigError::Preflight(Box::new(Self {
            severity: Severity::Error,
            key_path: key.into(),
            data,
            fix_hint: hint.into(),
        }))
    }
    pub fn from_error(error: ConfigError) -> Self {
        match error {
            ConfigError::Preflight(f) => *f,
            other => Self { severity: Severity::Error, key_path: "$".into(), data: FindingData::InputInvalid { message: other.to_string() }, fix_hint: "Correct the reported input error and repeat preflight; dependent checks have not completed.".into() },
        }
    }
}

/// Success preserves the legacy report; failure explicitly stops dependent checks.
pub fn inspect(input: &Path) -> serde_json::Value {
    let result = std::fs::read_to_string(input)
        .map_err(|e| ConfigError::Invalid(format!("{}: {e}", input.display())))
        .and_then(|text| Config::parse_unvalidated(&text))
        .and_then(|config| check(&config, input));
    match result {
        Ok(report) => serde_json::to_value(report).expect("finite preflight report"),
        Err(error) => failure(error),
    }
}

fn failure(error: ConfigError) -> serde_json::Value {
    serde_json::json!({
        "schema_version": 1, "findings_schema_version": FINDINGS_SCHEMA_VERSION,
        "valid": false, "coverage": "stopped_at_first_error",
        "findings": [Finding::from_error(error)]
    })
}

pub fn cli_check(
    input: &Path,
    json: bool,
    explain: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    let report = inspect(input);
    if json {
        println!("{}", serde_json::to_string_pretty(&report)?);
    } else if report["valid"] == true {
        let text = std::fs::read_to_string(input)?;
        let config = Config::parse(&text)?;
        for record in &config.input_expressions {
            println!("{}: {:?} = {}", record.key_path, record.text, record.value);
        }
        check(&config, input)?.print();
        if explain {
            println!(
                "Only deterministic rules decide validity. Findings are advisory unless severity is error; no edits are applied."
            );
        }
    } else {
        for finding in report["findings"].as_array().expect("findings array") {
            let finding: Finding = serde_json::from_value(finding.clone())?;
            eprintln!("{finding}");
        }
    }
    if report["valid"] == true {
        Ok(())
    } else {
        Err("preflight failed; dependent checks stopped".into())
    }
}

pub mod structure {
    //! Deterministic CLI geometry review. No thresholds affect numerical settings.
    use super::{Config, ConfigError, Finding, FindingData, Severity};
    use dft_geometry::Cell;
    use serde::{Deserialize, Serialize};
    use std::f64::consts::TAU;

    /// Optional advisory policies, in atomic units. No material-specific defaults.
    #[derive(Debug, Default, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct Policy {
        pub minimum_distance_bohr: Option<f64>,
        pub minimum_volume_bohr3: Option<f64>,
        pub minimum_volume_ratio: Option<f64>,
        pub slab_axis: Option<usize>,
        pub minimum_vacuum_bohr: Option<f64>,
        pub maximum_k_spacing_inverse_bohr: Option<f64>,
    }
    #[derive(Debug, Serialize)]
    pub struct Metrics {
        pub volume_bohr3: f64,
        pub volume_ratio: f64,
        pub minimum_distance_bohr: f64,
        /// Equal indices denote a nonzero periodic self image.
        pub nearest_atoms: [usize; 2],
        pub vacuum_gap_bohr: Option<f64>,
        /// None for arbitrary explicit k points; no mesh is inferred.
        pub k_spacing_inverse_bohr: Option<[f64; 3]>,
    }
    fn invalid(message: &str) -> ConfigError {
        ConfigError::Invalid(message.into())
    }
    fn norm(v: [f64; 3]) -> f64 {
        v[0].hypot(v[1]).hypot(v[2])
    }
    fn warning(key: &str, data: FindingData, hint: &str) -> Finding {
        Finding {
            severity: Severity::Warning,
            key_path: key.into(),
            data,
            fix_hint: hint.into(),
        }
    }

    pub fn review(config: &Config) -> Result<(Vec<Finding>, Metrics), ConfigError> {
        let policy = &config.preflight;
        for (key, value) in [
            ("minimum_distance_bohr", policy.minimum_distance_bohr),
            ("minimum_volume_bohr3", policy.minimum_volume_bohr3),
            ("minimum_volume_ratio", policy.minimum_volume_ratio),
            ("minimum_vacuum_bohr", policy.minimum_vacuum_bohr),
            (
                "maximum_k_spacing_inverse_bohr",
                policy.maximum_k_spacing_inverse_bohr,
            ),
        ] {
            if value.is_some_and(|v| !v.is_finite() || v <= 0.) {
                return Err(invalid(&format!(
                    "preflight.{key} must be finite and positive"
                )));
            }
        }
        if policy.minimum_volume_ratio.is_some_and(|v| v > 1.)
            || policy.slab_axis.is_some_and(|axis| axis > 2)
            || policy.minimum_vacuum_bohr.is_some() && policy.slab_axis.is_none()
        {
            return Err(invalid(
                "preflight: volume ratio must be <= 1; slab axis must be 0..2 and is required for a vacuum threshold",
            ));
        }
        let cell = config.cell().map_err(|e| Finding::error(
        "system.lattice_vectors",
        FindingData::CellVolume { volume_bohr3: None, minimum_bohr3: None, reason: e.to_string() },
        "Correct lattice vectors and units; the cell must have finite, representable nonzero volume.",
    ))?;
        if config.system.atoms.is_empty()
            || config
                .system
                .atoms
                .iter()
                .any(|a| a.fractional.iter().any(|x| !x.is_finite()))
        {
            return Err(invalid(
                "system.atoms must be nonempty with finite fractional coordinates",
            ));
        }
        let a = cell.lattice_columns_bohr();
        let lengths = std::array::from_fn(|i| norm([a[0][i], a[1][i], a[2][i]]));
        let b = cell.reciprocal_columns_inverse_bohr();
        let reciprocal: [f64; 3] = std::array::from_fn(|i| norm([b[0][i], b[1][i], b[2][i]]));
        let volume = cell.volume_bohr3();
        let ratio = lengths.iter().fold(volume, |v, length| v / length).min(1.);
        if !ratio.is_finite()
            || ratio <= 0.
            || lengths
                .iter()
                .chain(reciprocal.iter())
                .any(|v| !v.is_finite() || *v <= 0.)
        {
            return Err(invalid(
                "system.lattice_vectors: geometry metrics are not representable",
            ));
        }
        let mut findings = Vec::new();
        if let Some(minimum) = policy.minimum_volume_bohr3.filter(|min| volume < *min) {
            findings.push(warning(
                "system.lattice_vectors",
                FindingData::CellVolume {
                    volume_bohr3: Some(volume),
                    minimum_bohr3: Some(minimum),
                    reason: "below explicit review policy".into(),
                },
                "Review lattice units and intended cell volume.",
            ));
        }
        if let Some(minimum) = policy.minimum_volume_ratio.filter(|min| ratio < *min) {
            findings.push(warning(
                "system.lattice_vectors",
                FindingData::CellSkew {
                    volume_ratio: ratio,
                    minimum_ratio: minimum,
                },
                "Review cell conditioning; a basis change must preserve the structure.",
            ));
        }
        // Budget bounds work, not distance accuracy. Exhaustion fails explicitly.
        let mut budget = 1_000_000_u64;
        let mut nearest = [0, 0];
        let mut distance = closest(&cell, [0.; 3], true, &mut budget)?;
        for (i, atom) in config.system.atoms.iter().enumerate() {
            for (j, other) in config.system.atoms.iter().enumerate().skip(i + 1) {
                let delta = std::array::from_fn(|k| {
                    atom.fractional[k].rem_euclid(1.) - other.fractional[k].rem_euclid(1.)
                });
                let d = closest(&cell, delta, false, &mut budget)?;
                if d <= 1e-10 {
                    return Err(Finding::error(
                        "system.atoms",
                        FindingData::AtomOverlap {
                            atoms: [i, j],
                            distance_bohr: d,
                            tolerance_bohr: 1e-10,
                        },
                        "Remove unintended periodic duplicates or correct coordinates; do not automatically move atoms.",
                    ));
                }
                if d < distance {
                    distance = d;
                    nearest = [i, j];
                }
            }
        }
        if let Some(minimum) = policy.minimum_distance_bohr.filter(|min| distance < *min) {
            findings.push(warning(
                "system.atoms",
                FindingData::MinDistance {
                    atoms: nearest,
                    distance_bohr: distance,
                    minimum_bohr: minimum,
                },
                "Review coordinates, units and bonding; no atom is moved automatically.",
            ));
        }
        let gap = policy.slab_axis.map(|axis| {
            let mut coordinates: Vec<_> = config
                .system
                .atoms
                .iter()
                .map(|a| a.fractional[axis].rem_euclid(1.))
                .collect();
            coordinates.sort_by(f64::total_cmp);
            let cyclic = coordinates[0] + 1. - coordinates[coordinates.len() - 1];
            let gap = coordinates
                .windows(2)
                .map(|p| p[1] - p[0])
                .fold(cyclic, f64::max);
            gap * (TAU / reciprocal[axis])
        });
        if let (Some(axis), Some(gap), Some(minimum)) =
            (policy.slab_axis, gap, policy.minimum_vacuum_bohr)
            && gap < minimum
        {
            findings.push(warning("system.lattice_vectors", FindingData::SlabVacuum { axis, gap_bohr: gap, minimum_bohr: minimum }, "Consider more vacuum and test convergence; this is an empty nuclear-plane gap, not electronic vacuum thickness."));
        }
        let spacing = config
            .basis
            .kgrid
            .filter(|grid| grid.iter().all(|n| *n > 0))
            .map(|grid| std::array::from_fn(|i| reciprocal[i] / grid[i] as f64));
        if let (Some(spacing), Some(maximum)) = (spacing, policy.maximum_k_spacing_inverse_bohr)
            && spacing.iter().any(|s| *s > maximum)
        {
            findings.push(warning(
                "basis.kgrid",
                FindingData::KmeshCell {
                    lengths_bohr: lengths,
                    spacing_inverse_bohr: spacing,
                    maximum_inverse_bohr: maximum,
                },
                "Review reciprocal spacing and perform a k-mesh convergence study.",
            ));
        }
        Ok((
            findings,
            Metrics {
                volume_bohr3: volume,
                volume_ratio: ratio,
                minimum_distance_bohr: distance,
                nearest_atoms: nearest,
                vacuum_gap_bohr: gap,
                k_spacing_inverse_bohr: spacing,
            },
        ))
    }

    /// Exhaustive closest-vector search with rigorous inverse-row bounds.
    /// If |A(d+n)| <= R, then |d_i+n_i| <= R |row_i(A^-1)|.
    /// Unlike componentwise rounding, this also covers highly skewed cells.
    fn closest(
        cell: &Cell,
        delta: [f64; 3],
        self_image: bool,
        budget: &mut u64,
    ) -> Result<f64, ConfigError> {
        let d = delta.map(|x| x - x.round());
        let a = cell.lattice_columns_bohr();
        let radius = if self_image {
            (0..3)
                .map(|i| norm([a[0][i], a[1][i], a[2][i]]))
                .fold(f64::INFINITY, f64::min)
        } else {
            norm(cell.fractional_to_cartesian_bohr(d))
        };
        if !radius.is_finite() {
            return Err(invalid("structure distance is not representable"));
        }
        let b = cell.reciprocal_columns_inverse_bohr();
        let mut bounds = [(0_i64, 0_i64); 3];
        let mut count = 1_u64;
        for i in 0..3 {
            let bound = radius * (norm([b[0][i], b[1][i], b[2][i]]) / TAU);
            if !bound.is_finite() || bound > 1e6 {
                return Err(invalid(
                    "structure review search budget exceeded; reduce cell skew or use an equivalent reduced basis",
                ));
            }
            let lower = (-d[i] - bound - 1e-12).ceil() as i64;
            let upper = (-d[i] + bound + 1e-12).floor() as i64;
            count = count.saturating_mul((upper - lower + 1).max(0) as u64);
            bounds[i] = (lower, upper);
        }
        if count > *budget {
            return Err(invalid(
                "structure review search budget exceeded; use an equivalent reduced lattice or a separately reviewed larger search budget",
            ));
        }
        *budget -= count;
        let mut best = radius;
        for i in bounds[0].0..=bounds[0].1 {
            for j in bounds[1].0..=bounds[1].1 {
                for k in bounds[2].0..=bounds[2].1 {
                    if self_image && i == 0 && j == 0 && k == 0 {
                        continue;
                    }
                    best = best.min(norm(cell.fractional_to_cartesian_bohr([
                        d[0] + i as f64,
                        d[1] + j as f64,
                        d[2] + k as f64,
                    ])));
                }
            }
        }
        Ok(best)
    }
}

pub mod proposal {
    //! Offline proposal validation and immutable review records. Never edits inputs.
    use super::{Config, inspect};
    use serde::{Deserialize, Serialize};
    use serde_json::Value;
    use sha2::{Digest, Sha256};
    use std::{
        collections::BTreeMap,
        fs,
        io::Write,
        path::Path,
        time::{SystemTime, UNIX_EPOCH},
    };

    type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;
    pub fn digest(bytes: &[u8]) -> String {
        format!("{:x}", Sha256::digest(bytes))
    }

    #[derive(Debug, Serialize, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct Edit {
        /// Literal TOML table segments. Arrays are replaced as a whole.
        pub key_path: Vec<String>,
        pub before: Value,
        pub after: Value,
        pub reason: String,
        pub finding_codes: Vec<String>,
    }
    #[derive(Debug, Serialize, Deserialize)]
    #[serde(deny_unknown_fields)]
    pub struct Proposal {
        pub schema_version: u32,
        pub input_sha256: String,
        pub findings_sha256: String,
        pub pseudopotentials: BTreeMap<String, String>,
        pub candidate_pseudopotentials: BTreeMap<String, String>,
        pub advisory: String,
        pub diff: String,
        pub edits: Vec<Edit>,
        pub candidate_input_sha256: String,
        pub agent: String,
    }
    #[derive(Debug, Clone, Copy, Serialize, Deserialize, clap::ValueEnum)]
    #[serde(rename_all = "snake_case")]
    pub enum Decision {
        Approved,
        Rejected,
    }

    /// Exactly the bytes printed by check/validate --json, including trailing newline.
    pub fn findings_bytes(input: &Path) -> Result<Vec<u8>> {
        Ok(format!("{}\n", serde_json::to_string_pretty(&inspect(input))?).into_bytes())
    }
    pub fn potential_hashes(config: &Config, input: &Path) -> Result<BTreeMap<String, String>> {
        config
            .pseudopotentials
            .iter()
            .map(|(species, path)| {
                let path = input.parent().unwrap_or(Path::new(".")).join(path);
                Ok((species.clone(), digest(&fs::read(path)?)))
            })
            .collect()
    }
    fn nonempty(s: &str) -> bool {
        !s.trim().is_empty()
    }
    fn is_hash(s: &str) -> bool {
        s.len() == 64
            && s.bytes()
                .all(|b| b.is_ascii_digit() || (b'a'..=b'f').contains(&b))
    }

    /// A deliberately strict unified diff subset: one full-file hunk, exact context,
    /// deletions then additions, UTF-8 files with final newline. No paths are executed.
    fn candidate(source: &str, diff: &str) -> Result<String> {
        if !source.ends_with('\n') || !diff.ends_with('\n') {
            return Err("input and diff require final newlines".into());
        }
        let mut lines = diff.split_terminator('\n');
        if lines.next() != Some("--- a/input.toml") || lines.next() != Some("+++ b/input.toml") {
            return Err("diff must target only a/input.toml and b/input.toml".into());
        }
        let header = lines.next().ok_or("missing unified hunk")?;
        let mut old = String::new();
        let mut new = String::new();
        let mut adding = false;
        for line in lines {
            if let Some(line) = line.strip_prefix('-') {
                if adding {
                    return Err("deletions must precede additions".into());
                }
                old.push_str(line);
                old.push('\n');
            } else if let Some(line) = line.strip_prefix('+') {
                adding = true;
                new.push_str(line);
                new.push('\n');
            } else {
                return Err("only a full-file replacement hunk is supported".into());
            }
        }
        let expected = format!(
            "@@ -1,{} +1,{} @@",
            old.lines().count(),
            new.lines().count()
        );
        if header != expected || old != source || new.is_empty() || old == new {
            return Err(
                "diff hunk counts, exact original bytes or candidate content mismatch".into(),
            );
        }
        Ok(new)
    }

    fn changes(
        before: &Value,
        after: &Value,
        path: &mut Vec<String>,
        out: &mut BTreeMap<Vec<String>, (Value, Value)>,
    ) {
        if before == after {
            return;
        }
        if let (Some(a), Some(b)) = (before.as_object(), after.as_object()) {
            let keys: std::collections::BTreeSet<_> = a.keys().chain(b.keys()).collect();
            for key in keys {
                path.push(key.clone());
                changes(
                    a.get(key).unwrap_or(&Value::Null),
                    b.get(key).unwrap_or(&Value::Null),
                    path,
                    out,
                );
                path.pop();
            }
        } else {
            out.insert(path.clone(), (before.clone(), after.clone()));
        }
    }

    /// Validates syntax, freshness and the complete declared TOML delta. A rejected
    /// proposal may contain a candidate that still fails deterministic preflight.
    pub fn validate(input: &Path, proposal: &Proposal) -> Result<(String, Value)> {
        let source = fs::read_to_string(input)?;
        let config = Config::parse_unvalidated(&source)?;
        if proposal.schema_version != 1
            || !nonempty(&proposal.advisory)
            || !nonempty(&proposal.agent)
            || ![
                &proposal.input_sha256,
                &proposal.findings_sha256,
                &proposal.candidate_input_sha256,
            ]
            .iter()
            .all(|s| is_hash(s))
            || proposal.edits.is_empty()
        {
            return Err("invalid proposal schema, metadata or hash encoding".into());
        }
        if digest(source.as_bytes()) != proposal.input_sha256
            || potential_hashes(&config, input)? != proposal.pseudopotentials
            || digest(&findings_bytes(input)?) != proposal.findings_sha256
        {
            return Err("stale input, pseudopotentials or findings".into());
        }
        let text = candidate(&source, &proposal.diff)?;
        if digest(text.as_bytes()) != proposal.candidate_input_sha256 {
            return Err("candidate digest mismatch".into());
        }
        let candidate_config = Config::parse_unvalidated(&text)?;
        if potential_hashes(&candidate_config, input)? != proposal.candidate_pseudopotentials {
            return Err("candidate pseudopotential hashes mismatch".into());
        }
        let before: toml::Value = toml::from_str(&source)?;
        let after: toml::Value = toml::from_str(&text)?;
        let mut delta = BTreeMap::new();
        changes(
            &serde_json::to_value(before)?,
            &serde_json::to_value(after)?,
            &mut Vec::new(),
            &mut delta,
        );
        if delta.len() != proposal.edits.len() {
            return Err("undeclared or duplicate TOML edits".into());
        }
        let report = inspect(input);
        let codes: Vec<_> = report["findings"]
            .as_array()
            .ok_or("missing findings")?
            .iter()
            .filter_map(|v| v["code"].as_str())
            .collect();
        for edit in &proposal.edits {
            if edit.key_path.is_empty()
                || !nonempty(&edit.reason)
                || edit.finding_codes.is_empty()
                || edit
                    .finding_codes
                    .iter()
                    .any(|code| !codes.contains(&code.as_str()))
                || delta.remove(&edit.key_path) != Some((edit.before.clone(), edit.after.clone()))
            {
                return Err(
                    "edit must match the actual delta and cite findings with a reason".into(),
                );
            }
        }
        // Path resolution stays relative to the original input; nothing is staged or run.
        let recheck = match super::check(&candidate_config, input) {
            Ok(report) => serde_json::to_value(report)?,
            Err(error) => super::failure(error),
        };
        Ok((text, recheck))
    }

    /// Records an explicit operator decision, without applying even approved edits.
    /// The proposal and decision are one atomic, no-clobber sidecar beside results.json.
    pub fn record(
        input: &Path,
        proposal_path: &Path,
        output: &Path,
        decision: Decision,
        actor: &str,
        reason: &str,
    ) -> Result<std::path::PathBuf> {
        if !nonempty(actor) || !nonempty(reason) {
            return Err("decision requires actor and reason".into());
        }
        let bytes = fs::read(proposal_path)?;
        if bytes.len() > 4 * 1024 * 1024 {
            return Err("proposal exceeds 4 MiB".into());
        }
        let proposal: Proposal = serde_json::from_slice(&bytes)?;
        let (candidate, report) = validate(input, &proposal)?;
        if matches!(decision, Decision::Approved) && report["valid"] != true {
            return Err("cannot approve a candidate that fails deterministic preflight".into());
        }
        println!("{}\n{}", proposal.advisory, proposal.diff);
        let timestamp = SystemTime::now().duration_since(UNIX_EPOCH)?;
        let record = serde_json::json!({
            "schema_version":1,
            "proposal_sha256":digest(&bytes),
            "proposal":proposal,
            "proposal_source":std::str::from_utf8(&bytes)?,
            "decision":{"value":decision,"actor":actor,"reason":reason,"unix_time_ns":timestamp.as_nanos().to_string()},
            "candidate_input":candidate,
            "candidate_preflight":report,
            "applied":false
        });
        // Refuse symlink destinations and ancestors. No symlink or overwrite follows.
        for component in output.ancestors() {
            if component.as_os_str().is_empty() {
                continue;
            }
            if fs::symlink_metadata(component).is_ok_and(|m| m.file_type().is_symlink()) {
                return Err("audit destination must not contain symlinks".into());
            }
        }
        fs::create_dir_all(output)?;
        let target = output.join(format!("preflight-review-{}.json", digest(&bytes)));
        let temporary = output.join(format!(
            ".preflight-review-{}-{}.tmp",
            std::process::id(),
            timestamp.as_nanos()
        ));
        // Only clean up a temporary file that this call actually created.
        let mut file = fs::OpenOptions::new()
            .write(true)
            .create_new(true)
            .open(&temporary)?;
        let result = (|| -> Result<()> {
            file.write_all(&serde_json::to_vec_pretty(&record)?)?;
            file.write_all(b"\n")?;
            file.sync_all()?;
            // Detect source or UPF changes while the review was being prepared.
            validate(input, &proposal)?;
            fs::hard_link(&temporary, &target)?;
            fs::File::open(output)?.sync_all()?;
            Ok(())
        })();
        let _ = fs::remove_file(&temporary);
        result?;
        Ok(target)
    }

    /// Display the exact validated diff before an operator supplies a decision.
    pub fn preview(input: &Path, proposal_path: &Path) -> Result<()> {
        let bytes = fs::read(proposal_path)?;
        if bytes.len() > 4 * 1024 * 1024 {
            return Err("proposal exceeds 4 MiB".into());
        }
        let proposal: Proposal = serde_json::from_slice(&bytes)?;
        let (_, recheck) = validate(input, &proposal)?;
        println!(
            "Proposal SHA-256: {}\n{}\n{}\nCandidate valid: {}; no decision recorded, input unchanged",
            digest(&bytes),
            proposal.advisory,
            proposal.diff,
            recheck["valid"]
        );
        Ok(())
    }
}
