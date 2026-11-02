//! Explicit unit wrappers; internal lengths and energies use Bohr and Hartree.

/// CODATA 2018 Bohr radius, in Angstrom.
/// Source: <https://physics.nist.gov/cuu/pdf/all_2018.pdf>.
/// QE's constant must be checked independently in P6.
pub const BOHR_IN_ANGSTROM: f64 = 0.529_177_210_903;
/// CODATA 2018 Hartree energy in eV.
/// Source: <https://physics.nist.gov/cuu/pdf/all_2018.pdf>.
pub const HARTREE_IN_EV: f64 = 27.211_386_245_988;
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

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn conversions() {
        assert_eq!(Rydberg::from(Hartree(1.0)), Rydberg(2.0));
        assert_eq!(Hartree::from(Rydberg(2.0)), Hartree(1.0));
        assert_eq!(Angstrom::from(Bohr(1.0)).0, BOHR_IN_ANGSTROM);
        assert_eq!(ElectronVolt::from(Hartree(1.0)).0, HARTREE_IN_EV);
        let mut worst = 0.0_f64;
        for x in [-1e6, -3.7, 1e-9, 0.1, 1.0, 1e6] {
            for y in [
                Bohr::from(Angstrom::from(Bohr(x))).0,
                Angstrom::from(Bohr::from(Angstrom(x))).0,
                Hartree::from(Rydberg::from(Hartree(x))).0,
                Rydberg::from(Hartree::from(Rydberg(x))).0,
                Hartree::from(ElectronVolt::from(Hartree(x))).0,
                ElectronVolt::from(Hartree::from(ElectronVolt(x))).0,
            ] {
                worst = worst.max(((y - x) / x).abs());
            }
        }
        println!("unit round-trip max relative error={worst:.16e}");
        assert!(worst < 1e-15); // Dimensionless relative tolerance.
        assert_eq!(Bohr::from(Angstrom(0.0)), Bohr(0.0));
        assert_eq!(Hartree::from(ElectronVolt(0.0)), Hartree(0.0));
    }
}
