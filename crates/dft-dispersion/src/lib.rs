//! Grimme DFT-D3 dispersion correction through the official simple-dftd3
//! library (ADR 0048).
//!
//! The library (LGPL-3.0-or-later) is never linked: with feature `dftd3` it is
//! opened at run time (`dlopen` through `libloading`) by the thin adapter in
//! `ffi.rs`, which declares only the C functions used here. Without the feature
//! the types still exist (so callers can parse and validate input), and
//! [`compute`] returns [`Error::FeatureDisabled`].
//!
//! Units are Hartree atomic units throughout: positions and lattice in Bohr,
//! energy in Hartree, gradient dE/dR in Ha/Bohr, virial in Hartree.
//!
//! The virial returned by the library is the strain derivative of the energy,
//! `virial[a][b] = dE/d eps_ab` for the homogeneous deformation
//! `r -> (1 + eps) r` of positions and lattice vectors (verified by central
//! finite differences in `tests/finite_differences.rs`). The engine's
//! tensile-positive stress is therefore `virial / Omega`
//! (see [`Dispersion::stress_ha_bohr3`]).

#[cfg(feature = "dftd3")]
mod ffi;

use std::fmt;
use std::path::{Path, PathBuf};

/// Environment variable naming the shared library when the input does not.
pub const LIBRARY_ENV: &str = "BAANDI_DFTD3_LIB";

/// Platform default file name, resolved by the dynamic loader's search path.
#[cfg(target_os = "macos")]
pub const DEFAULT_LIBRARY: &str = "libs-dftd3.dylib";
/// Platform default file name, resolved by the dynamic loader's search path.
#[cfg(not(target_os = "macos"))]
pub const DEFAULT_LIBRARY: &str = "libs-dftd3.so";

/// Damping function of the two-body term.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Damping {
    /// Rational (Becke-Johnson) damping, input name `d3bj`.
    Rational,
    /// Zero (Chai-Head-Gordon) damping, input name `d3zero`.
    Zero,
}

impl Damping {
    /// Parse the input `method` value.
    pub fn from_method(method: &str) -> Option<Self> {
        match method {
            "d3bj" => Some(Self::Rational),
            "d3zero" => Some(Self::Zero),
            _ => None,
        }
    }
    /// Input name of the method.
    pub fn method(self) -> &'static str {
        match self {
            Self::Rational => "d3bj",
            Self::Zero => "d3zero",
        }
    }
    /// Names of the parameters that define the damping completely.
    pub fn core_parameters(self) -> &'static [&'static str] {
        match self {
            Self::Rational => &["s6", "s8", "a1", "a2"],
            Self::Zero => &["s6", "s8", "rs6"],
        }
    }
}

/// Optional user parameters. When none of the core parameters of the damping
/// ([`Damping::core_parameters`]) is given, the library's own parameter set for
/// the functional is loaded; then only `s9` may be given and must be 0 (no
/// three-body term) or 1 (three-body term on, the default). When all core
/// parameters are given they define the damping completely and `s9`
/// (default 1), `alpha6` (default 14) and, for zero damping, `rs8` (default 1)
/// may be given too. Anything in between is rejected: the C API cannot read
/// back a loaded parameter set, so a partial override would be guesswork.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Overrides {
    pub s6: Option<f64>,
    pub s8: Option<f64>,
    pub s9: Option<f64>,
    pub a1: Option<f64>,
    pub a2: Option<f64>,
    pub rs6: Option<f64>,
    pub rs8: Option<f64>,
    pub alpha6: Option<f64>,
}

/// Resolved parameter request passed to the library.
#[derive(Clone, Debug, PartialEq)]
pub enum Parameters {
    /// Library parameter set for `functional` (lower case), with or without the
    /// three-body (ATM) term.
    Library {
        damping: Damping,
        functional: String,
        three_body: bool,
    },
    /// Rational damping: s6, s8, s9, a1, a2 (Bohr), alpha.
    Rational {
        s6: f64,
        s8: f64,
        s9: f64,
        a1: f64,
        a2: f64,
        alpha: f64,
    },
    /// Zero damping: s6, s8, s9, rs6, rs8, alpha.
    Zero {
        s6: f64,
        s8: f64,
        s9: f64,
        rs6: f64,
        rs8: f64,
        alpha: f64,
    },
}

impl Parameters {
    /// Resolve the request. `functional` is `None` for a functional that D3
    /// does not parameterize (for example LDA); then every core parameter must
    /// be given explicitly.
    pub fn resolve(
        damping: Damping,
        functional: Option<&str>,
        o: &Overrides,
    ) -> Result<Self, Error> {
        let core: Vec<Option<f64>> = match damping {
            Damping::Rational => vec![o.s6, o.s8, o.a1, o.a2],
            Damping::Zero => vec![o.s6, o.s8, o.rs6],
        };
        let foreign = match damping {
            Damping::Rational => [("rs6", o.rs6), ("rs8", o.rs8)],
            Damping::Zero => [("a1", o.a1), ("a2", o.a2)],
        };
        if let Some((name, _)) = foreign.iter().find(|(_, v)| v.is_some()) {
            return Err(Error::Parameters(format!(
                "`{name}` does not belong to method `{}`",
                damping.method()
            )));
        }
        for (name, v) in [
            ("s6", o.s6),
            ("s8", o.s8),
            ("s9", o.s9),
            ("a1", o.a1),
            ("a2", o.a2),
            ("rs6", o.rs6),
            ("rs8", o.rs8),
            ("alpha6", o.alpha6),
        ] {
            if v.is_some_and(|v| !v.is_finite()) {
                return Err(Error::Parameters(format!("`{name}` must be finite")));
            }
        }
        let given = core.iter().filter(|v| v.is_some()).count();
        if given == 0 {
            let Some(functional) = functional else {
                return Err(Error::Parameters(format!(
                    "the functional has no D3 parameter set; give all of {} explicitly",
                    damping.core_parameters().join(", ")
                )));
            };
            if o.alpha6.is_some() || o.rs8.is_some() {
                return Err(Error::Parameters(
                    "`alpha6`/`rs8` can only be set together with all core parameters".into(),
                ));
            }
            let three_body = match o.s9 {
                None => true,
                Some(1.0) => true,
                Some(0.0) => false,
                Some(_) => {
                    return Err(Error::Parameters(
                        "with the library parameter set `s9` must be 0 or 1; give all core parameters for another value".into(),
                    ));
                }
            };
            return Ok(Self::Library {
                damping,
                functional: functional.to_ascii_lowercase(),
                three_body,
            });
        }
        if given != core.len() {
            return Err(Error::Parameters(format!(
                "partial damping parameters: give none (library set) or all of {}",
                damping.core_parameters().join(", ")
            )));
        }
        let s9 = o.s9.unwrap_or(1.0);
        let alpha = o.alpha6.unwrap_or(14.0);
        Ok(match damping {
            Damping::Rational => Self::Rational {
                s6: o.s6.unwrap(),
                s8: o.s8.unwrap(),
                s9,
                a1: o.a1.unwrap(),
                a2: o.a2.unwrap(),
                alpha,
            },
            Damping::Zero => Self::Zero {
                s6: o.s6.unwrap(),
                s8: o.s8.unwrap(),
                s9,
                rs6: o.rs6.unwrap(),
                rs8: o.rs8.unwrap_or(1.0),
                alpha,
            },
        })
    }
    /// Short description for reports.
    pub fn source(&self) -> &'static str {
        match self {
            Self::Library { .. } => "library_parameter_set",
            _ => "explicit_input",
        }
    }
}

/// Periodic structure in Bohr.
#[derive(Clone, Debug)]
pub struct Structure<'a> {
    /// Atomic numbers, one per atom.
    pub numbers: &'a [i32],
    /// Cartesian positions, Bohr.
    pub positions_bohr: &'a [[f64; 3]],
    /// Lattice vectors as rows (`lattice[i]` is vector i), Bohr. Callers with
    /// lattice vectors as matrix columns transpose at this boundary.
    pub lattice_rows_bohr: [[f64; 3]; 3],
}

/// Result of one evaluation.
#[derive(Clone, Debug, PartialEq)]
pub struct Dispersion {
    pub energy_ha: f64,
    /// dE/dR per atom, Ha/Bohr (the force is its negative).
    pub gradient_ha_bohr: Vec<[f64; 3]>,
    /// Strain derivative dE/d eps, Hartree (symmetric).
    pub virial_ha: [[f64; 3]; 3],
    /// Library version as reported by `dftd3_get_version`, e.g. "1.7.0".
    pub library_version: String,
    /// Library path that was opened.
    pub library: PathBuf,
}

impl Dispersion {
    /// Forces -dE/dR, Ha/Bohr.
    pub fn forces_ha_bohr(&self) -> Vec<[f64; 3]> {
        self.gradient_ha_bohr
            .iter()
            .map(|g| [-g[0], -g[1], -g[2]])
            .collect()
    }
    /// Tensile-positive stress (1/Omega) dE/d eps in Ha/Bohr^3, the engine's
    /// convention.
    pub fn stress_ha_bohr3(&self, volume_bohr3: f64) -> [[f64; 3]; 3] {
        std::array::from_fn(|a| std::array::from_fn(|b| self.virial_ha[a][b] / volume_bohr3))
    }
}

/// Errors of the adapter; every library message is surfaced verbatim.
#[derive(Clone, Debug, PartialEq)]
pub enum Error {
    /// Built without feature `dftd3`.
    FeatureDisabled,
    /// The shared library or a symbol could not be loaded.
    Load(String),
    /// Invalid parameter request.
    Parameters(String),
    /// Invalid structure (sizes, element numbers, non-finite values).
    Structure(String),
    /// Error reported by the library.
    Library(String),
}

impl fmt::Display for Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::FeatureDisabled => write!(
                f,
                "this build has no D3 dispersion support: rebuild with `--features dftd3`"
            ),
            Self::Load(m) => write!(f, "cannot load the simple-dftd3 library: {m}"),
            Self::Parameters(m) => write!(f, "dispersion parameters: {m}"),
            Self::Structure(m) => write!(f, "dispersion structure: {m}"),
            Self::Library(m) => write!(f, "simple-dftd3: {m}"),
        }
    }
}

impl std::error::Error for Error {}

/// Library path: explicit input, else `BAANDI_DFTD3_LIB`, else the platform
/// default name (searched by the dynamic loader).
pub fn library_path(explicit: Option<&Path>) -> PathBuf {
    if let Some(p) = explicit {
        return p.to_path_buf();
    }
    match std::env::var_os(LIBRARY_ENV) {
        Some(v) if !v.is_empty() => PathBuf::from(v),
        _ => PathBuf::from(DEFAULT_LIBRARY),
    }
}

/// Whether the build carries the adapter.
pub const fn available() -> bool {
    cfg!(feature = "dftd3")
}

/// Open the library and return its version string ("major.minor.patch").
pub fn library_version(path: &Path) -> Result<String, Error> {
    #[cfg(feature = "dftd3")]
    {
        ffi::Api::open(path).map(|api| api.version())
    }
    #[cfg(not(feature = "dftd3"))]
    {
        let _ = path;
        Err(Error::FeatureDisabled)
    }
}

/// Evaluate energy, gradient and virial of `structure` (periodic in all three
/// directions) with `parameters`, using the library at `path`.
pub fn compute(
    path: &Path,
    structure: &Structure<'_>,
    parameters: &Parameters,
) -> Result<Dispersion, Error> {
    validate(structure)?;
    #[cfg(feature = "dftd3")]
    {
        let api = ffi::Api::open(path)?;
        let (energy_ha, gradient_ha_bohr, virial_ha) = api.dispersion(structure, parameters)?;
        Ok(Dispersion {
            energy_ha,
            gradient_ha_bohr,
            virial_ha,
            library_version: api.version(),
            library: path.to_path_buf(),
        })
    }
    #[cfg(not(feature = "dftd3"))]
    {
        let _ = (path, parameters);
        Err(Error::FeatureDisabled)
    }
}

fn validate(s: &Structure<'_>) -> Result<(), Error> {
    if s.numbers.is_empty() || s.numbers.len() != s.positions_bohr.len() {
        return Err(Error::Structure(format!(
            "{} atomic numbers for {} positions",
            s.numbers.len(),
            s.positions_bohr.len()
        )));
    }
    if i32::try_from(s.numbers.len()).is_err() {
        return Err(Error::Structure("too many atoms".into()));
    }
    if let Some(z) = s.numbers.iter().find(|z| !(1..=MAX_ELEMENT).contains(*z)) {
        return Err(Error::Structure(format!(
            "atomic number {z} outside 1..={MAX_ELEMENT}"
        )));
    }
    if s.positions_bohr
        .iter()
        .flatten()
        .chain(s.lattice_rows_bohr.iter().flatten())
        .any(|v| !v.is_finite())
    {
        return Err(Error::Structure("non-finite position or lattice".into()));
    }
    Ok(())
}

/// Highest atomic number accepted here; simple-dftd3 reference data cover
/// hydrogen to plutonium (Z = 94).
pub const MAX_ELEMENT: i32 = 94;

const SYMBOLS: [&str; 94] = [
    "H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl",
    "Ar", "K", "Ca", "Sc", "Ti", "V", "Cr", "Mn", "Fe", "Co", "Ni", "Cu", "Zn", "Ga", "Ge", "As",
    "Se", "Br", "Kr", "Rb", "Sr", "Y", "Zr", "Nb", "Mo", "Tc", "Ru", "Rh", "Pd", "Ag", "Cd", "In",
    "Sn", "Sb", "Te", "I", "Xe", "Cs", "Ba", "La", "Ce", "Pr", "Nd", "Pm", "Sm", "Eu", "Gd", "Tb",
    "Dy", "Ho", "Er", "Tm", "Yb", "Lu", "Hf", "Ta", "W", "Re", "Os", "Ir", "Pt", "Au", "Hg", "Tl",
    "Pb", "Bi", "Po", "At", "Rn", "Fr", "Ra", "Ac", "Th", "Pa", "U", "Np", "Pu",
];

/// Atomic number of an element symbol (case-insensitive, surrounding blanks
/// ignored), e.g. from a UPF header.
pub fn atomic_number(symbol: &str) -> Option<i32> {
    let s = symbol.trim();
    SYMBOLS
        .iter()
        .position(|e| e.eq_ignore_ascii_case(s))
        .map(|i| i as i32 + 1)
}

/// Element symbol of an atomic number.
pub fn symbol(z: i32) -> Option<&'static str> {
    usize::try_from(z - 1)
        .ok()
        .and_then(|i| SYMBOLS.get(i).copied())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn element_table_round_trips() {
        assert_eq!(atomic_number("C"), Some(6));
        assert_eq!(atomic_number(" si "), Some(14));
        assert_eq!(atomic_number("Pu"), Some(94));
        assert_eq!(atomic_number("Xx"), None);
        for z in 1..=MAX_ELEMENT {
            assert_eq!(atomic_number(symbol(z).unwrap()), Some(z));
        }
    }

    #[test]
    fn parameter_resolution_rules() {
        let none = Overrides::default();
        assert_eq!(
            Parameters::resolve(Damping::Rational, Some("PBE"), &none).unwrap(),
            Parameters::Library {
                damping: Damping::Rational,
                functional: "pbe".into(),
                three_body: true
            }
        );
        assert!(Parameters::resolve(Damping::Rational, None, &none).is_err());
        let s9 = Overrides {
            s9: Some(0.0),
            ..none
        };
        assert!(matches!(
            Parameters::resolve(Damping::Zero, Some("pbe"), &s9).unwrap(),
            Parameters::Library {
                three_body: false,
                ..
            }
        ));
        let bad_s9 = Overrides {
            s9: Some(0.5),
            ..none
        };
        assert!(Parameters::resolve(Damping::Zero, Some("pbe"), &bad_s9).is_err());
        let partial = Overrides {
            s6: Some(1.0),
            ..none
        };
        assert!(Parameters::resolve(Damping::Rational, Some("pbe"), &partial).is_err());
        let full = Overrides {
            s6: Some(1.0),
            s8: Some(0.7875),
            a1: Some(0.4289),
            a2: Some(4.4407),
            ..none
        };
        assert_eq!(
            Parameters::resolve(Damping::Rational, None, &full).unwrap(),
            Parameters::Rational {
                s6: 1.0,
                s8: 0.7875,
                s9: 1.0,
                a1: 0.4289,
                a2: 4.4407,
                alpha: 14.0
            }
        );
        let foreign = Overrides {
            rs6: Some(1.2),
            ..full
        };
        assert!(Parameters::resolve(Damping::Rational, None, &foreign).is_err());
        let nan = Overrides {
            s8: Some(f64::NAN),
            ..full
        };
        assert!(Parameters::resolve(Damping::Rational, None, &nan).is_err());
    }

    #[test]
    fn structure_validation() {
        let s = Structure {
            numbers: &[6, 0],
            positions_bohr: &[[0.; 3], [1.; 3]],
            lattice_rows_bohr: [[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]],
        };
        assert!(matches!(validate(&s), Err(Error::Structure(_))));
        let s = Structure {
            numbers: &[6],
            positions_bohr: &[[0.; 3], [1.; 3]],
            ..s
        };
        assert!(validate(&s).is_err());
    }

    #[cfg(not(feature = "dftd3"))]
    #[test]
    fn without_feature_compute_reports_it() {
        let s = Structure {
            numbers: &[6],
            positions_bohr: &[[0.; 3]],
            lattice_rows_bohr: [[10., 0., 0.], [0., 10., 0.], [0., 0., 10.]],
        };
        let p = Parameters::resolve(Damping::Rational, Some("pbe"), &Overrides::default()).unwrap();
        assert_eq!(compute(Path::new("x"), &s, &p), Err(Error::FeatureDisabled));
    }
}
