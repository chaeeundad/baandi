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

#[derive(Debug, Serialize)]
pub struct Report {
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
        let key = format!("pseudopotentials.{label}");
        let path = input.parent().unwrap_or(Path::new(".")).join(path);
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
                return Err(fail(
                    &key,
                    message,
                    "correct the label/element or select the matching UPF; deliberate tests may set allow_element_mismatch = true",
                ));
            }
        } else if identity.is_some_and(|i| i.allow_element_mismatch) {
            overrides.push(format!(
                "{key}.allow_element_mismatch = true (permission enabled; no mismatch)"
            ));
        }
        if !options.functional.matches_upf(p.functional()) {
            let message = format!(
                "xc.functional: requested {}, {} has UPF functional {:?}",
                options.functional.name(),
                path.display(),
                p.functional()
            );
            if options.allow_functional_override {
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
        return Err(fail(
            "system.atoms / pseudopotentials",
            format!("odd electron count {electrons}; spin-unpolarized occupations require pairs"),
            "use an even-electron composition",
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
