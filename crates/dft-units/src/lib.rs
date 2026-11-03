//! Explicit unit wrappers; internal lengths and energies use Bohr and Hartree.

include!(concat!(env!("OUT_DIR"), "/factors.rs"));
/// Exact definition of the Rydberg energy unit.
pub const RYDBERG_PER_HARTREE: f64 = 2.0;

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Bohr(pub f64);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Angstrom(pub f64);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Hartree(pub f64);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Rydberg(pub f64);
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ElectronVolt(pub f64);

impl From<Bohr> for Angstrom {
    fn from(x: Bohr) -> Self {
        Self(x.0 * BOHR_IN_ANGSTROM)
    }
}
impl From<Angstrom> for Bohr {
    fn from(x: Angstrom) -> Self {
        Self(x.0 / BOHR_IN_ANGSTROM)
    }
}
impl From<Hartree> for Rydberg {
    fn from(x: Hartree) -> Self {
        Self(x.0 * RYDBERG_PER_HARTREE)
    }
}
impl From<Rydberg> for Hartree {
    fn from(x: Rydberg) -> Self {
        Self(x.0 / RYDBERG_PER_HARTREE)
    }
}
impl From<Hartree> for ElectronVolt {
    fn from(x: Hartree) -> Self {
        Self(x.0 * HARTREE_IN_EV)
    }
}
impl From<ElectronVolt> for Hartree {
    fn from(x: ElectronVolt) -> Self {
        Self(x.0 / HARTREE_IN_EV)
    }
}

use std::{fmt, str::FromStr};
use uom::si::{
    energy,
    f64::{Energy, Length, Pressure, Ratio},
    length, pressure, ratio,
};

/// Errors retain the complete offending input; no units are inferred.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct UnitError {
    pub kind: UnitErrorKind,
    pub text: String,
}
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum UnitErrorKind {
    Empty,
    MissingUnit,
    InvalidNumber,
    NonFinite,
    UnknownUnit,
    ExtraTokens,
}
impl fmt::Display for UnitError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "{:?} in unit input {:?}", self.kind, self.text)
    }
}
impl std::error::Error for UnitError {}
fn error(kind: UnitErrorKind, text: &str) -> UnitError {
    UnitError {
        kind,
        text: text.into(),
    }
}

macro_rules! units {
    ($name:ident { $($variant:ident => [$($alias:literal),+]),+ $(,)? }) => {
        #[derive(Debug, Clone, Copy, PartialEq, Eq)]
        pub enum $name { $($variant),+ }
        impl $name { pub const ALL: &'static [Self] = &[$(Self::$variant),+]; }
        impl FromStr for $name {
            type Err = UnitError;
            fn from_str(s: &str) -> Result<Self, Self::Err> {
                match s.trim().to_ascii_lowercase().as_str() {
                    $($($alias)|+ => Ok(Self::$variant),)+
                    _ => Err(error(UnitErrorKind::UnknownUnit, s)),
                }
            }
        }
    }
}
units!(EnergyUnit {
    Hartree => ["ha", "hartree"], Rydberg => ["ry", "rydberg"],
    ElectronVolt => ["ev", "electronvolt"], MilliElectronVolt => ["mev", "millielectronvolt"], Joule => ["j", "joule"]
});
units!(LengthUnit {
    Bohr => ["bohr", "a0"], Angstrom => ["angstrom", "ang", "a"],
    Nanometer => ["nm", "nanometer"], Picometer => ["pm", "picometer"], Meter => ["m", "meter"]
});
units!(PressureUnit {
    HartreePerBohrCubed => ["ha/bohr^3"], RydbergPerBohrCubed => ["ry/bohr^3"],
    ElectronVoltPerAngstromCubed => ["ev/angstrom^3"], Gigapascal => ["gpa"], Kilobar => ["kbar"], Pascal => ["pa"]
});
fn energy_quantity(unit: EnergyUnit) -> Energy {
    match unit {
        EnergyUnit::Hartree => Energy::new::<energy::hartree>(1.),
        EnergyUnit::Rydberg => Energy::new::<energy::hartree>(1. / RYDBERG_PER_HARTREE),
        EnergyUnit::ElectronVolt => Energy::new::<energy::electronvolt>(1.),
        EnergyUnit::MilliElectronVolt => Energy::new::<energy::electronvolt>(
            Ratio::new::<ratio::per_mille>(1.).get::<ratio::ratio>(),
        ),
        EnergyUnit::Joule => Energy::new::<energy::joule>(1.),
    }
}
fn length_quantity(unit: LengthUnit) -> Length {
    match unit {
        LengthUnit::Bohr => Length::new::<length::bohr_radius>(1.),
        LengthUnit::Angstrom => Length::new::<length::angstrom>(1.),
        LengthUnit::Nanometer => Length::new::<length::nanometer>(1.),
        LengthUnit::Picometer => Length::new::<length::picometer>(1.),
        LengthUnit::Meter => Length::new::<length::meter>(1.),
    }
}
fn pressure_quantity(unit: PressureUnit) -> Pressure {
    let b = length_quantity(LengthUnit::Bohr);
    let a = length_quantity(LengthUnit::Angstrom);
    match unit {
        PressureUnit::HartreePerBohrCubed => energy_quantity(EnergyUnit::Hartree) / (b * b * b),
        PressureUnit::RydbergPerBohrCubed => energy_quantity(EnergyUnit::Rydberg) / (b * b * b),
        PressureUnit::ElectronVoltPerAngstromCubed => {
            energy_quantity(EnergyUnit::ElectronVolt) / (a * a * a)
        }
        PressureUnit::Gigapascal => Pressure::new::<pressure::gigapascal>(1.),
        PressureUnit::Kilobar => {
            Pressure::new::<pressure::bar>(Ratio::new::<ratio::ratio>(1.).get::<ratio::per_mille>())
        }
        PressureUnit::Pascal => Pressure::new::<pressure::pascal>(1.),
    }
}
/// f64 conversion; finite representable round trips are tested to relative 1e-14.
pub fn convert_energy(value: f64, from: EnergyUnit, to: EnergyUnit) -> f64 {
    value * (energy_quantity(from) / energy_quantity(to)).get::<ratio::ratio>()
}
pub fn convert_length(value: f64, from: LengthUnit, to: LengthUnit) -> f64 {
    value * (length_quantity(from) / length_quantity(to)).get::<ratio::ratio>()
}
pub fn convert_pressure(value: f64, from: PressureUnit, to: PressureUnit) -> f64 {
    value * (pressure_quantity(from) / pressure_quantity(to)).get::<ratio::ratio>()
}
pub fn to_hartree(value: f64, from: EnergyUnit) -> f64 {
    convert_energy(value, from, EnergyUnit::Hartree)
}
pub fn from_hartree(value: f64, to: EnergyUnit) -> f64 {
    convert_energy(value, EnergyUnit::Hartree, to)
}
pub fn to_bohr(value: f64, from: LengthUnit) -> f64 {
    convert_length(value, from, LengthUnit::Bohr)
}
pub fn from_bohr(value: f64, to: LengthUnit) -> f64 {
    convert_length(value, LengthUnit::Bohr, to)
}

/// Parse a decimal/scientific number followed by exactly one unit token.
/// The longest valid numeric prefix handles adjacent units, including `1eV`.
pub fn parse_quantity<U: FromStr<Err = UnitError>>(text: &str) -> Result<(f64, U), UnitError> {
    let s = text.trim();
    if s.is_empty() {
        return Err(error(UnitErrorKind::Empty, text));
    }
    let (end, value) = s
        .char_indices()
        .map(|(i, _)| i)
        .chain(std::iter::once(s.len()))
        .filter(|&i| i > 0)
        .filter_map(|i| s[..i].parse::<f64>().ok().map(|v| (i, v)))
        .next_back()
        .ok_or_else(|| error(UnitErrorKind::InvalidNumber, text))?;
    if !value.is_finite() {
        return Err(error(UnitErrorKind::NonFinite, text));
    }
    let unit = s[end..].trim();
    if unit.is_empty() {
        return Err(error(UnitErrorKind::MissingUnit, text));
    }
    if unit.split_whitespace().count() != 1 {
        return Err(error(UnitErrorKind::ExtraTokens, text));
    }
    Ok((
        value,
        unit.parse().map_err(|e: UnitError| error(e.kind, text))?,
    ))
}
pub fn parse_energy(text: &str) -> Result<f64, UnitError> {
    let (v, u) = parse_quantity::<EnergyUnit>(text)?;
    finite_result(to_hartree(v, u), text)
}
pub fn parse_length(text: &str) -> Result<f64, UnitError> {
    let (v, u) = parse_quantity::<LengthUnit>(text)?;
    finite_result(to_bohr(v, u), text)
}
fn finite_result(value: f64, text: &str) -> Result<f64, UnitError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(error(UnitErrorKind::NonFinite, text))
    }
}
