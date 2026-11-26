//! Central finite differences of E_disp against the library's analytic
//! gradient and virial (ticket D3-DISPERSION): h = 1e-4 Bohr for positions,
//! strain step 1e-5 for the cell. Establishes the sign conventions used by the
//! engine: force = -gradient, and tensile-positive stress
//! (1/Omega) dE/d eps = virial / Omega for r -> (1 + eps) r applied to positions
//! and lattice vectors.
//!
//! About 40 s (some 300 library calls, dominated by the library itself, so
//! release mode does not shorten it much), hence ignored in the gate. Run:
//! `BAANDI_DFTD3_LIB=... cargo test --release -p dft-dispersion --features dftd3
//! --test finite_differences -- --ignored --nocapture`.
#![cfg(feature = "dftd3")]
// Index loops mirror the (atom, axis) and (a, b) tensor notation.
#![allow(clippy::needless_range_loop)]

mod common;

use common::Fixture;
use dft_dispersion::{Damping, Overrides, Parameters, compute};

const H: f64 = 1e-4;
const STRAIN: f64 = 1e-5;
const FORCE_TOL: f64 = 1e-6;
const STRESS_TOL: f64 = 1e-6;

fn energy(lib: &std::path::Path, f: &Fixture, p: &Parameters) -> f64 {
    compute(lib, &f.structure(), p).unwrap().energy_ha
}

fn strained(f: &Fixture, a: usize, b: usize, e: f64) -> Fixture {
    let apply = |r: [f64; 3]| -> [f64; 3] {
        let mut out = r;
        out[a] += e * r[b];
        out
    };
    Fixture {
        name: f.name,
        numbers: f.numbers.clone(),
        positions: f.positions.iter().map(|r| apply(*r)).collect(),
        lattice_rows: f.lattice_rows.map(apply),
    }
}

#[test]
#[ignore = "about 40 s; run with --ignored (see the module documentation)"]
fn analytic_derivatives_match_central_differences() {
    let Some(lib) = common::library() else { return };
    let mut worst = (0.0_f64, 0.0_f64);
    for damping in [Damping::Rational, Damping::Zero] {
        let p = Parameters::resolve(damping, Some("pbe"), &Overrides::default()).unwrap();
        for f in common::all() {
            let d = compute(&lib, &f.structure(), &p).unwrap();
            let forces = d.forces_ha_bohr();
            let mut df = 0.0_f64;
            for i in 0..f.numbers.len() {
                for c in 0..3 {
                    let mut plus = Fixture {
                        name: f.name,
                        numbers: f.numbers.clone(),
                        positions: f.positions.clone(),
                        lattice_rows: f.lattice_rows,
                    };
                    plus.positions[i][c] += H;
                    let mut minus = Fixture {
                        positions: f.positions.clone(),
                        numbers: f.numbers.clone(),
                        ..plus
                    };
                    minus.positions[i][c] -= H;
                    let fd = -(energy(&lib, &plus, &p) - energy(&lib, &minus, &p)) / (2. * H);
                    df = df.max((fd - forces[i][c]).abs());
                }
            }
            let volume = f.volume();
            let stress = d.stress_ha_bohr3(volume);
            let mut ds = 0.0_f64;
            for a in 0..3 {
                for b in 0..3 {
                    let fd = (energy(&lib, &strained(&f, a, b, STRAIN), &p)
                        - energy(&lib, &strained(&f, a, b, -STRAIN), &p))
                        / (2. * STRAIN)
                        / volume;
                    ds = ds.max((fd - stress[a][b]).abs());
                    // The opposite sign convention would fail by 2|stress|.
                }
            }
            println!(
                "{:6} {:18} max|F_fd-F| {:9.2e} Ha/Bohr  max|sigma_fd-sigma| {:9.2e} Ha/Bohr^3  (|sigma|max {:9.2e})",
                damping.method(),
                f.name,
                df,
                ds,
                stress.iter().flatten().fold(0.0_f64, |m, v| m.max(v.abs()))
            );
            assert!(df <= FORCE_TOL, "{} forces: {df:e}", f.name);
            assert!(ds <= STRESS_TOL, "{} stress: {ds:e}", f.name);
            worst = (worst.0.max(df), worst.1.max(ds));
        }
    }
    println!(
        "worst: force {:.2e} Ha/Bohr, stress {:.2e} Ha/Bohr^3",
        worst.0, worst.1
    );
}
