//! GK18 timing harness (ignored): per-k cost of the k basis and its FFT grid
//! map on a GeTe-shaped skewed cell (80 Ry, 48^3). Observation only.
//! `cargo test --release -p dft-basis --test gk18_basis_timing -- --ignored --nocapture`
use dft_basis::{FftGridSpec, KBasisSummary, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_units::Hartree;
use std::time::Instant;

fn gete() -> Cell {
    let a = 1.889_726_124_565_062;
    let rows = [
        [4.28 * a, 0., 0.],
        [2.249_018_666_589_306_4 * a, 3.641_471_548_334_939_4 * a, 0.],
        [
            2.249_018_666_589_306_4 * a,
            1.254_359_637_225_141_8 * a,
            3.418_610_381_110_031 * a,
        ],
    ];
    Cell::from_user_rows_bohr(rows).unwrap()
}

#[test]
#[ignore]
fn timing() {
    let cell = gete();
    let cutoff = Hartree(40.);
    let density = PlaneWaveBasis::new(&cell, Hartree(160.)).unwrap();
    let dims = FftGridSpec::new(&density, &density, None).unwrap().dims();
    let ks: Vec<[f64; 3]> = (0..64)
        .map(|i| {
            let x = i as f64;
            [
                ((x * 0.137).sin()) * 0.49,
                ((x * 0.291).cos()) * 0.49,
                ((x * 0.533).sin()) * 0.49,
            ]
        })
        .collect();
    let t = Instant::now();
    let bases: Vec<_> = ks
        .iter()
        .map(|k| PlaneWaveBasis::at_k(&cell, cutoff, *k).unwrap())
        .collect();
    let basis_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    for b in &bases {
        std::hint::black_box(FftGridSpec::new(b, &density, Some(dims)).unwrap());
    }
    let grid_s = t.elapsed().as_secs_f64();
    let t = Instant::now();
    for k in &ks {
        let s = KBasisSummary::at_k(&cell, cutoff, *k).unwrap();
        std::hint::black_box(
            FftGridSpec::dims_for_summary(&s, &density, Some(dims), || {
                PlaneWaveBasis::at_k(&cell, cutoff, *k)
            })
            .unwrap(),
        );
    }
    let summary_s = t.elapsed().as_secs_f64();
    println!(
        "summary + dims_for_summary {:.3} ms/k",
        1e3 * summary_s / ks.len() as f64
    );
    println!(
        "dims {dims:?} npw {} density {}: at_k {:.3} ms/k, FftGridSpec {:.3} ms/k",
        bases[0].vectors().len(),
        density.vectors().len(),
        1e3 * basis_s / ks.len() as f64,
        1e3 * grid_s / ks.len() as f64
    );
}
