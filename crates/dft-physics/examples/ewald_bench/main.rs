//! Ewald micro-benchmark (INIT-EWALD-FAST). `old.rs` is the evaluator of
//! commit eecd57f verbatim (before INIT-SCALE-3); `ewald::energy` is the
//! current one. Cells: si64 (conventional 2x2x2 Si, a = 20.5224 Bohr, 64 ions)
//! by default. Single thread; observations only (acceptance numbers come from
//! dev08).
//!
//! `cargo run --release -p dft-physics --example ewald_bench [reps]`
mod old;

use dft_geometry::Cell;
use dft_physics::ewald::{self, EwaldOptions, Ion, Scheme};
use std::time::Instant;

fn si64() -> (Cell, Vec<Ion>) {
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
    let mut ions = Vec::new();
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                for b in basis {
                    ions.push(Ion {
                        charge: 4.0,
                        fractional_position: [
                            (b[0] + i as f64) / 2.0,
                            (b[1] + j as f64) / 2.0,
                            (b[2] + k as f64) / 2.0,
                        ],
                    });
                }
            }
        }
    }
    (cell, ions)
}

fn time<T>(reps: usize, mut f: impl FnMut() -> T) -> (f64, T) {
    let mut best = f64::INFINITY;
    let mut out = None;
    for _ in 0..reps {
        let t = Instant::now();
        let v = f();
        best = best.min(t.elapsed().as_secs_f64());
        out = Some(v);
    }
    (best, out.unwrap())
}

/// Cells of the generality set (rows in Bohr, scaled for the ion count) with
/// `n` ions at deterministic spread-out positions and mixed charges.
fn shaped(shape: usize, n: usize) -> (Cell, Vec<Ion>) {
    let rows = [
        [[1., 0., 0.], [0., 1., 0.], [0., 0., 1.]],
        [[1., 0.07, 0.03], [0.15, 1.1, 0.05], [0.08, 0.1, 1.2]],
        [[0.7, 0., 0.], [0., 0.8, 0.], [0., 0., 2.2]],
        [[1.2, 0., 0.], [0., 1.3, 0.], [0.05, 0., 3.0]],
    ][shape];
    // About 20 Bohr^3 per ion (diamond Si density).
    let det = {
        let m: nalgebra::Matrix3<f64> = nalgebra::Matrix3::from_fn(|i, j| rows[i][j]);
        m.determinant().abs()
    };
    let s = (20.0 * n as f64 / det).cbrt();
    let cell = Cell::from_user_rows_bohr(rows.map(|r| r.map(|x| x * s))).unwrap();
    let ions = (0..n)
        .map(|i| Ion {
            charge: [4.0, 6.0, 1.0, 3.5][i % 4],
            fractional_position: [
                (i as f64 * 0.618034).fract(),
                (i as f64 * 0.754878 + 0.1).fract(),
                (i as f64 * 0.569840 + 0.2).fract(),
            ],
        })
        .collect();
    (cell, ions)
}

fn main() {
    let reps: usize = std::env::args()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .unwrap_or(3);
    let (cell, ions) = si64();
    let options = EwaldOptions::default();
    let (t_old, e_old) = time(reps, || old::old_energy(&cell, &ions, options).unwrap());
    let (t_ref, e_ref) = time(reps, || {
        ewald::energy_with(&cell, &ions, options, Scheme::Reference).unwrap()
    });
    let (t_bal, e_bal) = time(reps, || {
        ewald::energy_with(&cell, &ions, options, Scheme::Balanced).unwrap()
    });
    let (t_force_ref, _) = time(reps, || {
        ewald::set_scheme_override(Some(Scheme::Reference));
        let f = ewald::forces(&cell, &ions, options).unwrap();
        ewald::set_scheme_override(None);
        f
    });
    let (t_force_bal, _) = time(reps, || {
        ewald::set_scheme_override(Some(Scheme::Balanced));
        let f = ewald::forces(&cell, &ions, options).unwrap();
        ewald::set_scheme_override(None);
        f
    });
    println!(
        "si64 energy: eecd57f {t_old:.4} s | reference (pruned) {t_ref:.4} s, bits equal {} | balanced {t_bal:.5} s, E diff {:.3e} Ha",
        e_old.energy.0.to_bits() == e_ref.energy.0.to_bits(),
        e_bal.energy.0 - e_ref.energy.0
    );
    println!("si64 forces: reference {t_force_ref:.4} s | balanced {t_force_bal:.5} s");
    for scheme in [Scheme::Reference, Scheme::Balanced] {
        let plan = ewald::EwaldPlan::with_scheme(&cell, &ions, options, scheme).unwrap();
        let (mut real, mut reciprocal) = (Vec::new(), Vec::new());
        plan.real_terms(0..plan.real_candidates(), &mut real)
            .unwrap();
        plan.reciprocal_terms(0..plan.reciprocal_candidates(), &mut reciprocal);
        println!(
            "  {scheme:?}: real candidates {} terms {} | reciprocal candidates {} terms {}",
            plan.real_candidates(),
            real.len(),
            plan.reciprocal_candidates(),
            reciprocal.len()
        );
    }
    println!(
        "  reference eta {:.4} r_cut {:.2} g_cut {:.3} | balanced eta {:.4} r_cut {:.2} g_cut {:.3} tail {:.2e}",
        e_ref.eta,
        e_ref.real_cutoff_bohr,
        e_ref.reciprocal_cutoff_inverse_bohr,
        e_bal.eta,
        e_bal.real_cutoff_bohr,
        e_bal.reciprocal_cutoff_inverse_bohr,
        e_bal.tail_bound.0
    );
    // Cost-ratio scan: balanced time with eta from ratio c on the generality set.
    let ratios = [1.0, 4.0, 8.0, 12.0, 16.0, 32.0, 64.0];
    println!("cost ratio scan (balanced seconds per ratio {ratios:?}):");
    for shape in 0..4 {
        for n in [1usize, 8, 64, 216] {
            let (cell, ions) = shaped(shape, n);
            let length = cell.volume_bohr3().cbrt();
            let mut row = Vec::new();
            for c in ratios {
                let eta = std::f64::consts::PI.sqrt() * (n as f64 * c).powf(1.0 / 6.0) / length;
                let o = EwaldOptions {
                    eta: Some(eta),
                    ..options
                };
                let (t, _) = time(reps, || {
                    ewald::energy_with(&cell, &ions, o, Scheme::Balanced).unwrap()
                });
                row.push(format!("{t:.5}"));
            }
            let (t_ref, _) = time(1, || {
                ewald::energy_with(&cell, &ions, options, Scheme::Reference).unwrap()
            });
            println!(
                "  shape {shape} n {n:>3}: {} | reference {t_ref:.4}",
                row.join(" ")
            );
        }
    }
}
