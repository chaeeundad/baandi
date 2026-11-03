use std::{env, fs, path::PathBuf};
use uom::si::{energy, f64::Energy, f64::Length, length};
fn main() {
    let bohr = Length::new::<length::bohr_radius>(1.).get::<length::angstrom>();
    let hartree = Energy::new::<energy::hartree>(1.).get::<energy::electronvolt>();
    fs::write(PathBuf::from(env::var_os("OUT_DIR").unwrap()).join("factors.rs"),
        format!("/// Bohr in Angstrom, derived from uom.\npub const BOHR_IN_ANGSTROM: f64 = {bohr:?};\n/// Hartree in eV, derived from uom.\npub const HARTREE_IN_EV: f64 = {hartree:?};\n")).unwrap();
    println!("cargo:rerun-if-changed=build.rs");
}
