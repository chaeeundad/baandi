//! R1 many-empty-band solve (`docs/roadmap/rpa-crpa-plan.md`, ticket R1):
//! dense complete EVD versus Chebyshev subspace iteration versus reference
//! Davidson at a fixed density. Synthetic data are not QE reference evidence;
//! the fixed density is a nonuniform supplied field or the superposition of
//! atomic charges, never an SCF fixed point (solver agreement does not need one).
//!
//! Fast gate (debug): `cargo test -p dft-engine --test refscf_nscf_many`
//! Real UPFs (release, ignored):
//! DFT_SI_UPF=/abs/Si.pz-vbc.UPF DFT_NC_DIR=/abs/pseudodojo/ex \
//!   cargo test --release -p dft-engine --test refscf_nscf_many -- --ignored --nocapture
#[allow(dead_code)]
mod common;
use dft_engine::refscf::nscf::{
    BandOptions, BandSolver, Bands, NoStore, NscfResult, nscf_bands_with,
};
use dft_engine::{refstart::StartInput, *};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_physics::data::{Array3, Complex64};

struct System<'a> {
    name: &'static str,
    cell: Cell,
    atoms: Vec<Atom>,
    potentials: Vec<SpeciesPotential<'a>>,
    options: ScfOptions,
    density: Array3<f64>,
    points: Vec<KPoint>,
}
impl System<'_> {
    fn input(&self) -> StartInput<'_> {
        StartInput {
            cell: &self.cell,
            atoms: &self.atoms,
            potentials: &self.potentials,
            options: &self.options,
        }
    }
    fn run(&self, options: BandOptions, ranks: usize) -> Vec<NscfResult> {
        dft_mpi::simulate(ranks, |comm| {
            nscf_bands_with(
                self.input(),
                std::slice::from_ref(&self.density),
                &self.points,
                options,
                comm,
                &mut NoStore,
            )
            .unwrap()
        })
    }
    fn serial(&self, options: BandOptions) -> NscfResult {
        nscf_bands_with(
            self.input(),
            std::slice::from_ref(&self.density),
            &self.points,
            options,
            &refscf::metric::Serial,
            &mut NoStore,
        )
        .unwrap()
    }
}

fn options(functional: Functional, cutoff: f64) -> ScfOptions {
    ScfOptions {
        cutoff_hartree: cutoff,
        bands: Some(8),
        eigen_tolerance: 1e-8,
        functional,
        ..ScfOptions::strict_reference()
    }
}

/// Exactly charge-normalized nonuniform density on the problem's grid.
fn fixed_density(
    cell: &Cell,
    atoms: &[Atom],
    pp: &[SpeciesPotential<'_>],
    options: &ScfOptions,
    points: &[KPoint],
    atomic: bool,
) -> Array3<f64> {
    let problem = ScfProblem::with_kpoints(
        cell,
        atoms,
        pp,
        options.clone(),
        points,
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap();
    let mut density = problem.starting_density();
    if atomic {
        return density;
    }
    let shape = density.shape().to_vec();
    let scale = problem.electron_count() / cell.volume_bohr3();
    for ((i, j, k), n) in density.indexed_iter_mut() {
        *n = scale
            * (1.
                + 0.1 * (std::f64::consts::TAU * i as f64 / shape[0] as f64).cos()
                + 0.04
                    * (std::f64::consts::TAU
                        * (j as f64 / shape[1] as f64 + k as f64 / shape[2] as f64))
                        .sin());
    }
    density
}

fn gamma_mesh(n: [usize; 3]) -> Vec<KPoint> {
    let total = (n[0] * n[1] * n[2]) as f64;
    let mut points = Vec::new();
    for i in 0..n[0] {
        for j in 0..n[1] {
            for k in 0..n[2] {
                points.push(KPoint {
                    fractional: [
                        i as f64 / n[0] as f64,
                        j as f64 / n[1] as f64,
                        k as f64 / n[2] as f64,
                    ],
                    weight: 1. / total,
                });
            }
        }
    }
    points
}

fn max_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0., f64::max)
}

fn bits(b: &Bands) -> Vec<u64> {
    b.eigenvalues
        .iter()
        .chain(&b.residuals)
        .map(|v| v.to_bits())
        .chain(
            b.eigenvectors
                .as_slice()
                .iter()
                .flat_map(|z| [z.re.to_bits(), z.im.to_bits()]),
        )
        .collect()
}

/// max |X X^H - I| and max |X^H X - I| of one-rank complete eigenvectors.
fn completeness(x: &dft_numerics::linalg::ComplexMatrix) -> (f64, f64) {
    use dft_numerics::linalg::{MatrixOp, gemm};
    let n = x.rows();
    assert_eq!(x.cols(), n, "complete basis requires bands = basis size");
    let mut worst = [0_f64; 2];
    for (slot, (a, b)) in [
        (MatrixOp::Identity, MatrixOp::Adjoint),
        (MatrixOp::Adjoint, MatrixOp::Identity),
    ]
    .into_iter()
    .enumerate()
    {
        let p = gemm(x, a, x, b).unwrap();
        for j in 0..n {
            for i in 0..n {
                let expected = Complex64::new(f64::from(i == j), 0.);
                worst[slot] = worst[slot].max((p[(i, j)] - expected).norm());
            }
        }
    }
    (worst[0], worst[1])
}

/// Checks the R1 gates on one system; returns a one-line summary.
fn gates(
    system: &System<'_>,
    low: usize,
    chebyshev_bands: usize,
    ranks: &[usize],
    complete: bool,
) -> String {
    let basis: Vec<usize> = system
        .serial(BandOptions::new(low, 1e-6, BandSolver::Dense))
        .stats[0]
        .iter()
        .map(|s| s.basis_size)
        .collect();
    let smallest = *basis.iter().min().unwrap();
    // Dense: every eigenpair of the smallest basis (complete there).
    let all = if complete {
        smallest
    } else {
        chebyshev_bands + 8
    };
    let dense = system.serial(BandOptions::new(all, 1e-6, BandSolver::Dense));
    let mut dense_residual = 0_f64;
    let mut defect = (0_f64, 0_f64);
    for (k, b) in dense.channels[0].iter().enumerate() {
        dense_residual = b.residuals.iter().copied().fold(dense_residual, f64::max);
        if complete && basis[k] == all {
            let (left, right) = completeness(&b.eigenvectors);
            defect = (defect.0.max(left), defect.1.max(right));
        }
    }
    assert!(
        dense_residual <= 1e-9,
        "{}: dense residual {dense_residual:e}",
        system.name
    );
    if complete {
        assert!(
            defect.0 <= 1e-10 && defect.1 <= 1e-10,
            "{}: completeness {defect:?}",
            system.name
        );
    }
    // Reference Davidson (the existing strict path) on the low block.
    let davidson = system.serial(BandOptions::new(low, 1e-8, BandSolver::Davidson));
    let mut dav = 0_f64;
    for (d, e) in davidson.channels[0].iter().zip(&dense.channels[0]) {
        dav = dav.max(max_diff(&d.eigenvalues, &e.eigenvalues[..low]));
    }
    assert!(dav <= 1e-8, "{}: Davidson vs dense {dav:e}", system.name);
    // Chebyshev with buffer bands and locking against the dense EVD.
    let cheb_options = BandOptions::new(chebyshev_bands, 1e-6, BandSolver::Chebyshev);
    let cheb = system.serial(cheb_options);
    let mut iterative = 0_f64;
    let mut cheb_residual = 0_f64;
    for (c, e) in cheb.channels[0].iter().zip(&dense.channels[0]) {
        iterative = iterative.max(max_diff(&c.eigenvalues, &e.eigenvalues[..chebyshev_bands]));
        cheb_residual = c.residuals.iter().copied().fold(cheb_residual, f64::max);
    }
    assert!(
        iterative <= 1e-8,
        "{}: Chebyshev vs dense {iterative:e}",
        system.name
    );
    assert!(
        cheb_residual <= 1e-6,
        "{}: Chebyshev residual {cheb_residual:e}",
        system.name
    );
    // Rank counts and the one-rank simulator against the serial adapter.
    let mut rank_dense = 0_f64;
    let mut rank_cheb = 0_f64;
    let dense_ranked = BandOptions::new(chebyshev_bands, 1e-6, BandSolver::Dense);
    let dense_serial = system.serial(dense_ranked);
    for &r in ranks {
        for (answers, serial, worst) in [
            (system.run(dense_ranked, r), &dense_serial, &mut rank_dense),
            (system.run(cheb_options, r), &cheb, &mut rank_cheb),
        ] {
            for answer in &answers {
                for (a, b) in answer.channels[0].iter().zip(&serial.channels[0]) {
                    *worst = worst.max(max_diff(&a.eigenvalues, &b.eigenvalues));
                    if r == 1 {
                        assert!(
                            bits(a) == bits(b),
                            "{}: np=1 not bit-identical",
                            system.name
                        );
                    }
                }
            }
        }
    }
    assert!(
        rank_dense <= 1e-10,
        "{}: dense rank spread {rank_dense:e}",
        system.name
    );
    assert!(
        rank_cheb <= 1e-10,
        "{}: Chebyshev rank spread {rank_cheb:e}",
        system.name
    );
    let sweeps: usize = cheb.stats[0].iter().map(|s| s.iterations).sum();
    let columns: u64 = cheb.stats[0].iter().map(|s| s.h_columns).sum();
    format!(
        "{}: k={} N_PW {}..{}; dense bands {all}, max residual {dense_residual:.2e} Ha, completeness |XX^H-I| {:.2e} |X^HX-I| {:.2e}; Davidson({low}) vs dense {dav:.2e} Ha; Chebyshev {chebyshev_bands} bands vs dense {iterative:.2e} Ha, residual {cheb_residual:.2e} Ha, {sweeps} sweeps, {columns} H columns; rank spread {ranks:?} dense {rank_dense:.2e} Chebyshev {rank_cheb:.2e} Ha; np=1 bit-identical",
        system.name,
        system.points.len(),
        smallest,
        basis.iter().max().unwrap(),
        defect.0,
        defect.1,
    )
}

fn path() -> Vec<KPoint> {
    [[0.; 3], [0.5, 0., 0.5], [0.13, -0.21, 0.37]]
        .into_iter()
        .map(|fractional| KPoint {
            fractional,
            weight: 1. / 3.,
        })
        .collect()
}

#[test]
fn synthetic_skewed_two_species_gates() {
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0.7, 5.5, 0.], [0.3, 0.4, 6.5]]).unwrap();
    let p = common::potential_with_atomic_start(2., 1.5, true, 101, "SLA PZ");
    let q = common::potential_with_atomic_start(2., 1.8, false, 101, "SLA PZ");
    let atoms = vec![
        Atom::new("Xx", [0.13, 0.21, 0.37]),
        Atom::new("Yy", [0.61, 0.52, 0.74]),
    ];
    let potentials = vec![
        SpeciesPotential {
            label: "Xx",
            potential: &p,
        },
        SpeciesPotential {
            label: "Yy",
            potential: &q,
        },
    ];
    let o = options(Functional::LdaPz, 4.0);
    let points = path();
    let density = fixed_density(&cell, &atoms, &potentials, &o, &points, false);
    let system = System {
        name: "synthetic-skew-two-species",
        cell,
        atoms,
        potentials,
        options: o,
        density,
        points,
    };
    println!("{}", gates(&system, 8, 40, &[1, 2, 4], true));
}

#[test]
fn rejects_band_requests_beyond_the_basis() {
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
    let p = common::potential_with_atomic_start(2., 1.5, true, 101, "SLA PZ");
    let atoms = vec![Atom::new("Xx", [0.1, 0.2, 0.3])];
    let potentials = vec![SpeciesPotential {
        label: "Xx",
        potential: &p,
    }];
    let o = options(Functional::LdaPz, 1.5);
    let points = path();
    let density = fixed_density(&cell, &atoms, &potentials, &o, &points, false);
    let system = System {
        name: "cubic",
        cell,
        atoms,
        potentials,
        options: o,
        density,
        points,
    };
    let sizes: Vec<usize> = system
        .serial(BandOptions::new(5, 1e-6, BandSolver::Dense))
        .stats[0]
        .iter()
        .map(|s| s.basis_size)
        .collect();
    let smallest = *sizes.iter().min().unwrap();
    for solver in [BandSolver::Dense, BandSolver::Chebyshev] {
        assert!(
            nscf_bands_with(
                system.input(),
                std::slice::from_ref(&system.density),
                &system.points,
                BandOptions::new(smallest + 1, 1e-6, solver),
                &refscf::metric::Serial,
                &mut NoStore,
            )
            .is_err()
        );
        // Below occupied + 4 bands is rejected by the start, as for Davidson.
        assert!(
            nscf_bands_with(
                system.input(),
                std::slice::from_ref(&system.density),
                &system.points,
                BandOptions::new(4, 1e-6, solver),
                &refscf::metric::Serial,
                &mut NoStore,
            )
            .is_err()
        );
    }
    // The complete basis is accepted at the smallest k.
    let complete = nscf_bands_with(
        system.input(),
        std::slice::from_ref(&system.density),
        &[KPoint {
            fractional: system.points[0].fractional,
            weight: 1.,
        }],
        BandOptions::new(sizes[0], 1e-6, BandSolver::Dense),
        &refscf::metric::Serial,
        &mut NoStore,
    )
    .unwrap();
    assert_eq!(complete.channels[0][0].eigenvalues.len(), sizes[0]);
}

fn real(path: &str) -> dft_upf::Pseudopotential {
    dft_upf::parse_upf(&std::fs::read_to_string(path).unwrap()).unwrap()
}

#[test]
#[ignore = "external Si UPF; DFT_SI_UPF=/abs/Si.pz-vbc.UPF cargo test --release -p dft-engine --test refscf_nscf_many silicon -- --ignored --nocapture"]
fn silicon_mesh_and_skewed_cell() {
    let pp = real(&std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF"));
    let potentials = vec![SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let o = options(Functional::LdaPz, 5.0);
    // Si 2-atom fcc, Gamma-centered 4x4x4 (64 k, full mesh, no reduction).
    let points = gamma_mesh([4, 4, 4]);
    let cell = common::cell();
    let atoms = common::atoms();
    let density = fixed_density(&cell, &atoms, &potentials, &o, &points, true);
    let system = System {
        name: "Si2-fcc-4x4x4",
        cell,
        atoms,
        potentials: vec![SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }],
        options: o.clone(),
        density,
        points,
    };
    println!("{}", gates(&system, 8, 60, &[1, 2, 4], true));
    // Skewed/strained Si cell: monoclinic shear of the fcc vectors, displaced atom.
    let a = 10.261_212_856_717_933;
    let cell = Cell::from_user_rows_bohr([
        [0.3, a / 2.0, a / 2.0 + 0.2],
        [a / 2.0 - 0.1, 0.4, a / 2.0],
        [a / 2.0, a / 2.0 + 0.3, -0.2],
    ])
    .unwrap();
    let atoms = vec![
        Atom::new("Xx", [0.01, -0.02, 0.0]),
        Atom::new("Xx", [0.26, 0.24, 0.255]),
    ];
    let points = gamma_mesh([2, 2, 2]);
    let density = fixed_density(&cell, &atoms, &potentials, &o, &points, true);
    let system = System {
        name: "Si2-skewed-2x2x2",
        cell,
        atoms,
        potentials,
        options: o,
        density,
        points,
    };
    println!("{}", gates(&system, 8, 60, &[1, 2, 4], true));
}

#[test]
#[ignore = "external PseudoDojo PBE UPFs; DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-engine --test refscf_nscf_many dojo -- --ignored --nocapture"]
fn dojo_lif_and_aln() {
    let dir = std::path::PathBuf::from(std::env::var("DFT_NC_DIR").expect("DFT_NC_DIR"));
    let load = |s: &str| real(dir.join(format!("{s}.upf")).to_str().unwrap());
    let (li, f, al, n) = (load("Li"), load("F"), load("Al"), load("N"));
    let mut o = options(Functional::Pbe, 8.0);
    o.allow_functional_override = false;
    o.bands = Some(20);
    // LiF rock salt (a = 4.03 Angstrom = 7.6157 Bohr), Gamma-centered 2x2x2.
    let a = 7.6157;
    let cell = Cell::from_user_rows_bohr([
        [0.0, a / 2.0, a / 2.0],
        [a / 2.0, 0.0, a / 2.0],
        [a / 2.0, a / 2.0, 0.0],
    ])
    .unwrap();
    let atoms = vec![Atom::new("Li", [0.0; 3]), Atom::new("F", [0.5; 3])];
    let potentials = vec![
        SpeciesPotential {
            label: "Li",
            potential: &li,
        },
        SpeciesPotential {
            label: "F",
            potential: &f,
        },
    ];
    let points = gamma_mesh([2, 2, 2]);
    let density = fixed_density(&cell, &atoms, &potentials, &o, &points, true);
    let system = System {
        name: "LiF-rocksalt-2x2x2",
        cell,
        atoms,
        potentials,
        options: o.clone(),
        density,
        points,
    };
    println!("{}", gates(&system, 12, 40, &[1, 2, 4], true));
    // Wurtzite AlN (a = 3.11 A, c/a = 1.601, u = 0.382), Gamma + (0, 0, 1/2).
    let (a, c, u) = (5.8770, 5.8770 * 1.601, 0.382);
    let cell = Cell::from_user_rows_bohr([
        [a, 0., 0.],
        [-a / 2., a * 3f64.sqrt() / 2., 0.],
        [0., 0., c],
    ])
    .unwrap();
    let atoms = vec![
        Atom::new("Al", [1. / 3., 2. / 3., 0.]),
        Atom::new("Al", [2. / 3., 1. / 3., 0.5]),
        Atom::new("N", [1. / 3., 2. / 3., u]),
        Atom::new("N", [2. / 3., 1. / 3., 0.5 + u]),
    ];
    let potentials = vec![
        SpeciesPotential {
            label: "Al",
            potential: &al,
        },
        SpeciesPotential {
            label: "N",
            potential: &n,
        },
    ];
    let points = vec![
        KPoint {
            fractional: [0.; 3],
            weight: 0.5,
        },
        KPoint {
            fractional: [0., 0., 0.5],
            weight: 0.5,
        },
    ];
    let density = fixed_density(&cell, &atoms, &potentials, &o, &points, true);
    let system = System {
        name: "AlN-wurtzite",
        cell,
        atoms,
        potentials,
        options: o,
        density,
        points,
    };
    println!("{}", gates(&system, 16, 80, &[1, 2, 4], false));
}
