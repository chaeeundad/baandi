#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential};
use dft_geometry::{Atom, Cell};

fn model(kind: usize) -> (Cell, Vec<Atom>, dft_upf::Pseudopotential) {
    let rows = match kind {
        0 => [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]],
        1 => [[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]],
        2 => [[4., 0., 0.], [0., 4., 0.], [0., 0., 8.]],
        _ => [[5., 0., 0.], [0., 5., 0.], [0., 0., 10.]],
    };
    let atoms = if kind < 2 {
        vec![Atom::new("Xx", [0.13, 0.21, 0.32])]
    } else {
        vec![
            Atom::new("Xx", [0.13, 0.21, 0.17]),
            Atom::new("Xx", [0.6, 0.7, 0.68]),
        ]
    };
    (
        Cell::from_user_rows_bohr(rows).unwrap(),
        atoms,
        common::potential(2., 1., kind != 2, 401),
    )
}
fn options(mode: GammaReal, threads: usize, seed: u64) -> ScfOptions {
    ScfOptions {
        gamma_real: mode,
        threads,
        seed,
        cutoff_hartree: 2.5,
        bands: Some(5),
        ..Default::default()
    }
}
fn problem(kind: usize, opts: ScfOptions) -> ScfProblem {
    let (cell, atoms, pp) = model(kind);
    ScfProblem::new(
        &cell,
        &atoms,
        &[SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        opts,
    )
    .unwrap()
}
#[allow(dead_code)]
fn run(kind: usize, mode: GammaReal, threads: usize, seed: u64) -> ScfResult {
    problem(kind, options(mode, threads, seed)).run().unwrap()
}

fn save(name: &str, threads: usize, r: ScfResult) {
    let dir = std::env::var("P33A_SNAPSHOTS").expect("set snapshot directory");
    std::fs::create_dir_all(&dir).unwrap();
    let mut bytes = Vec::new();
    for v in std::iter::once(r.total_hartree())
        .chain(r.density.iter().copied())
        .chain(r.eigenvalues.iter().copied())
    {
        bytes.extend(v.to_bits().to_le_bytes());
    }
    bytes.extend((r.history.len() as u64).to_le_bytes());
    for h in &r.history {
        bytes.extend((h.davidson_iterations as u64).to_le_bytes());
        bytes.extend((h.operator_applications as u64).to_le_bytes());
    }
    // All solver fields, including certificate counts and eigenvectors; no wall clocks.
    for k in &r.k_points {
        bytes.extend(format!("{:?}", k.solution).as_bytes());
    }
    std::fs::write(format!("{dir}/{name}-T{threads}.bin"), bytes).unwrap();
}
#[test]
#[ignore = "P33A before/after snapshot plumbing; release and external DFT_SI_UPF required"]
fn snapshots() {
    let pp =
        dft_upf::parse_upf(&std::fs::read_to_string(std::env::var("DFT_SI_UPF").unwrap()).unwrap())
            .unwrap();
    for standard in [false, true] {
        for threads in [1, 2, 4, 8, 12] {
            for kind in 0..4 {
                save(&format!("synthetic-{kind}-standard-{standard}"), threads, {
                    let mut o = options(GammaReal::On, threads, 2026);
                    o.standard = standard.then(dft_engine::StandardProfile::default);
                    problem(kind, o).run().unwrap()
                });
            }
            for size in [2, 8, 16] {
                let si8 = size != 2;
                let a = 10.261_212_856_717_933;
                let cell = if si8 {
                    Cell::from_user_rows_bohr([
                        [a, 0., 0.],
                        [0., a, 0.],
                        [0., 0., if size == 16 { 2. * a } else { a }],
                    ])
                    .unwrap()
                } else {
                    common::cell()
                };
                let positions = if si8 {
                    vec![
                        [0., 0., 0.],
                        [0., 0.5, 0.5],
                        [0.5, 0., 0.5],
                        [0.5, 0.5, 0.],
                        [0.25, 0.25, 0.25],
                        [0.25, 0.75, 0.75],
                        [0.75, 0.25, 0.75],
                        [0.75, 0.75, 0.25],
                    ]
                } else {
                    vec![[0.; 3], [0.25; 3]]
                };
                let positions = if size == 16 {
                    positions
                        .into_iter()
                        .flat_map(|p| [[p[0], p[1], p[2] / 2.], [p[0], p[1], (p[2] + 1.) / 2.]])
                        .collect()
                } else {
                    positions
                };
                let atoms: Vec<_> = positions.into_iter().map(|p| Atom::new("Si", p)).collect();
                let r = ScfProblem::new(
                    &cell,
                    &atoms,
                    &[SpeciesPotential {
                        label: "Si",
                        potential: &pp,
                    }],
                    ScfOptions {
                        gamma_real: GammaReal::On,
                        standard: standard.then(dft_engine::StandardProfile::default),
                        threads,
                        cutoff_hartree: 15.,
                        max_iterations: 100,
                        mixing: dft_engine::Mixing::Pulay {
                            beta: 0.5,
                            history: 6,
                        },
                        ..Default::default()
                    },
                )
                .unwrap()
                .run()
                .unwrap();
                save(&format!("si{size}-standard-{standard}"), threads, r);
            }
        }
    }
}
