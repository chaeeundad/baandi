//! `ionic_local` micro-benchmark (INIT-IONIC-FAST): si64 cell (64 Si,
//! conventional 2x2x2, a = 20.5224 Bohr), wavefunction cutoff 15 Ha (density
//! sphere 60 Ha), direct phases (default) against the phase tables
//! (`BAANDI_IONIC_FAST`), single thread. Observations only.
//!
//! `cargo run --release -p dft-physics --example ionic_bench -- <Si UPF path> [reps]`
//! (the UPF is read at run time and never stored in the repository).
use dft_basis::PlaneWaveBasis;
use dft_geometry::Cell;
use dft_physics::hamiltonian::{self, Species};
use dft_units::Hartree;
use std::time::Instant;

fn main() {
    let mut args = std::env::args().skip(1);
    let Some(path) = args.next() else {
        eprintln!("usage: ionic_bench <Si UPF path> [reps]");
        return;
    };
    let reps: usize = args.next().and_then(|s| s.parse().ok()).unwrap_or(3);
    let pp = dft_upf::parse_upf(&std::fs::read_to_string(path).unwrap()).unwrap();
    let a = 20.522425713435865;
    let cell = Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap();
    let basis = [
        [0.0, 0.0, 0.0],
        [0.0, 0.5, 0.5],
        [0.5, 0.0, 0.5],
        [0.5, 0.5, 0.0],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ];
    let mut positions = Vec::new();
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                for b in basis {
                    positions.push(cell.fractional_to_cartesian_bohr([
                        (b[0] + f64::from(i)) / 2.0,
                        (b[1] + f64::from(j)) / 2.0,
                        (b[2] + f64::from(k)) / 2.0,
                    ]));
                }
            }
        }
    }
    let density = PlaneWaveBasis::new(&cell, Hartree(15.0))
        .unwrap()
        .density_basis(&cell)
        .unwrap();
    let species = [Species {
        potential: &pp,
        positions_bohr: &positions,
    }];
    let run = |fast: bool| {
        hamiltonian::set_ionic_fast_override(Some(fast));
        let mut best = f64::INFINITY;
        let mut out = Vec::new();
        for _ in 0..reps {
            let t = Instant::now();
            out = hamiltonian::ionic_local(&density, cell.volume_bohr3(), &species).unwrap();
            best = best.min(t.elapsed().as_secs_f64());
        }
        hamiltonian::set_ionic_fast_override(None);
        (best, out)
    };
    let (t_direct, direct) = run(false);
    let (t_fast, fast) = run(true);
    let scale = direct.iter().map(|z| z.norm()).fold(0.0, f64::max);
    let diff = direct
        .iter()
        .zip(&fast)
        .map(|(a, b)| (a - b).norm())
        .fold(0.0, f64::max);
    println!(
        "si64 ionic_local ({} G): direct {t_direct:.4} s | tables {t_fast:.4} s | max |dV| {diff:.3e} Ha (max |V| {scale:.3e})",
        density.vectors().len()
    );
}
