#[allow(dead_code)]
mod common;
use common::potential;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_engine::{OccupationPolicy, ScfError, ScfOptions, ScfProblem, SpeciesPotential};
use dft_geometry::{
    Atom, Cell,
    kpoints::{KPoint, monkhorst_pack, reduce_time_reversal},
};
use dft_physics::{
    data::*,
    hamiltonian::{Hamiltonian, NonlocalOperator, Species},
};
use dft_solvers::{DenseEigensolver, LinearOperator};
use dft_upf::radial::real_spherical_harmonic;

fn cell() -> Cell {
    Cell::from_user_rows_bohr([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]).unwrap()
}
fn policy() -> OccupationPolicy {
    OccupationPolicy::Insulator {
        minimum_gap_hartree: 1e-7,
    }
}
fn points(k: [f64; 3]) -> Vec<KPoint> {
    vec![KPoint {
        fractional: k,
        weight: 1.0,
    }]
}

#[test]
fn free_electron_skew_and_nonlocal_dense() {
    let cell =
        Cell::from_user_rows_bohr([[4.0, 0.0, 0.0], [1.1, 3.7, 0.0], [0.4, 0.8, 4.1]]).unwrap();
    let k = [0.21, -0.17, 0.13];
    let basis = PlaneWaveBasis::at_k(&cell, Hartree(2.0), k).unwrap();
    let db = basis.density_basis(&cell).unwrap();
    let grid = FftGridSpec::new(&basis, &db, None).unwrap();
    let zero = Array3::zeros(grid.dims());
    let h = Hamiltonian::new(
        &cell,
        &basis,
        zero.clone(),
        NonlocalOperator::new(&basis, cell.volume_bohr3(), &[]).unwrap(),
    )
    .unwrap();
    let free = DenseEigensolver.solve(&h, basis.vectors().len()).unwrap();
    // Independent wide integer enumeration; never derive expected kinetic from momenta().
    let b = cell.reciprocal_columns_inverse_bohr();
    let mut expected = Vec::new();
    for i in -8..=8 {
        for j in -8..=8 {
            for l in -8..=8 {
                let m = [i, j, l];
                let q: [f64; 3] = std::array::from_fn(|a| {
                    (0..3).map(|v| b[a][v] * (f64::from(m[v]) + k[v])).sum()
                });
                let e = q.iter().map(|x| x * x / 2.0).sum::<f64>();
                if e <= 2.0 {
                    expected.push(e);
                }
            }
        }
    }
    expected.sort_by(f64::total_cmp);
    assert_eq!(expected.len(), free.eigenvalues.len());
    let error = expected
        .iter()
        .zip(&free.eigenvalues)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    println!(
        "free electron skew: {} waves, max eigenvalue error={error:.16e} Ha",
        expected.len()
    );
    assert!(error < 1e-12);
    let pp = potential(2.0, 0.8, true, 401);
    let positions = [[0.23, 0.61, 1.13], [1.41, 0.37, 2.1]];
    let species = [Species {
        potential: &pp,
        positions_bohr: &positions,
    }];
    let nl = NonlocalOperator::new(&basis, cell.volume_bohr3(), &species).unwrap();
    let h = Hamiltonian::new(&cell, &basis, zero, nl).unwrap();
    let n = basis.vectors().len();
    let mut matrix = vec![vec![Complex64::default(); n]; n];
    for j in 0..n {
        let mut unit = vec![Complex64::default(); n];
        unit[j].re = 1.0;
        let mut y = unit.clone();
        h.apply(&unit, &mut y).unwrap();
        for i in 0..n {
            matrix[i][j] = y[i];
        }
    }
    // Build beta from fractional k, integer G, radial transforms and angular
    // functions independently of the production projector columns.
    let mut dense = vec![vec![Complex64::default(); n]; n];
    for position in positions {
        let mut columns = Vec::new();
        let mut channels = Vec::new();
        for (p, projector) in pp.projectors().iter().enumerate() {
            let l = projector.l();
            for m in -(l as isize)..=l as isize {
                channels.push((p, l, m));
                columns.push(
                    basis
                        .vectors()
                        .iter()
                        .map(|g| {
                            let q: [f64; 3] = std::array::from_fn(|a| {
                                (0..3)
                                    .map(|v| b[a][v] * (f64::from(g.miller[v]) + k[v]))
                                    .sum()
                            });
                            let norm = q.iter().map(|v| v * v).sum::<f64>().sqrt();
                            let angle = q.iter().zip(position).map(|(q, r)| q * r).sum::<f64>();
                            Complex64::new(0.0, -1.0).powu(l as u32)
                                * real_spherical_harmonic(l, m, q).unwrap()
                                * pp.projector_transform(p, norm).unwrap()
                                / cell.volume_bohr3().sqrt()
                                * Complex64::from_polar(1.0, -angle)
                        })
                        .collect::<Vec<_>>(),
                );
            }
        }
        for (a, &(p, l, m)) in channels.iter().enumerate() {
            for (c, &(r, ll, mm)) in channels.iter().enumerate() {
                if l == ll && m == mm {
                    for i in 0..n {
                        for j in 0..n {
                            dense[i][j] += columns[a][i]
                                * pp.dij_hartree()[p * pp.projectors().len() + r]
                                * columns[c][j].conj();
                        }
                    }
                }
            }
        }
    }
    let mut hermitian: f64 = 0.0;
    let mut error: f64 = 0.0;
    for i in 0..n {
        for j in 0..n {
            if i == j {
                dense[i][j].re += basis.momenta()[i].norm_squared / 2.0;
            }
            error = error.max((matrix[i][j] - dense[i][j]).norm());
            hermitian = hermitian.max((matrix[i][j] - matrix[j][i].conj()).norm());
        }
    }
    println!("nonlocal k: Hermiticity={hermitian:.16e} Ha, dense error={error:.16e} Ha");
    assert!(hermitian < 1e-12);
    assert!(error < 1e-10);
}

#[test]
fn gamma_pre_p7_regression() {
    let pp = potential(2.0, 1.0, false, 401);
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let options = ScfOptions {
        cutoff_hartree: 1.5,
        ..Default::default()
    };
    let legacy = ScfProblem::new(&cell(), &atoms, &species, options.clone())
        .unwrap()
        .run()
        .unwrap();
    let k = ScfProblem::with_kpoints(
        &cell(),
        &atoms,
        &species,
        options,
        &points([0.0; 3]),
        policy(),
    )
    .unwrap()
    .run()
    .unwrap();
    let reference = -1.848_228_969_297_197_4; // Measured on the unmodified pre-P7A functional-label test.
    println!(
        "Gamma: pre-P7={reference:.16e} Ha, legacy={:.16e} Ha, k API={:.16e} Ha, gap={:?} Ha",
        legacy.total_hartree(),
        k.total_hartree(),
        k.gap_hartree
    );
    assert!((legacy.total_hartree() - reference).abs() < 1e-10);
    assert!((k.total_hartree() - reference).abs() < 1e-10);
    assert_eq!(
        legacy.total_hartree().to_bits(),
        k.total_hartree().to_bits()
    );
}

#[test]
#[ignore = "multiple tight SCF solves; cargo test --release -p dft-engine --test kpoints -- --ignored --nocapture"]
fn scf_translation_and_time_reversal() {
    let a = potential(1.0, 0.55, true, 401);
    let b = potential(1.0, 1.15, true, 401);
    let species = [
        SpeciesPotential {
            label: "A",
            potential: &a,
        },
        SpeciesPotential {
            label: "B",
            potential: &b,
        },
    ];
    // Different species at generic positions: no inversion center for this motif.
    let atoms = [
        Atom::new("A", [0.11, 0.21, 0.17]),
        Atom::new("B", [0.39, 0.57, 0.68]),
    ];
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(3),
        density_tolerance: 1e-12,
        energy_tolerance: 1e-13,
        eigen_tolerance: 1e-12,
        initial_eigen_tolerance: 1e-12,
        max_iterations: 300,
        ..Default::default()
    };
    let run = |p: &[KPoint]| {
        ScfProblem::with_kpoints(&cell(), &atoms, &species, options.clone(), p, policy())
            .unwrap()
            .run()
            .unwrap()
    };
    let k = run(&points([0.125, -0.25, 0.375]));
    let translated = run(&points([1.125, -2.25, 3.375]));
    let de = (k.total_hartree() - translated.total_hartree()).abs();
    let dn = k
        .density
        .iter()
        .zip(&translated.density)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let ev = k
        .eigenvalues
        .iter()
        .zip(&translated.eigenvalues)
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    println!(
        "k+G0: energy={de:.16e} Ha eigenvalues={ev:.16e} Ha density={dn:.16e} electron/Bohr^3"
    );
    assert!(de < 1e-10 && ev < 1e-10 && dn < 1e-10);
    for unequal in [false, true] {
        let mut full = monkhorst_pack([2; 3], [0.5; 3]).unwrap();
        if unequal {
            let total = (1..=full.len()).map(|i| (i * i) as f64).sum::<f64>();
            for (i, p) in full.iter_mut().enumerate() {
                p.weight = ((i + 1) * (i + 1)) as f64 / total;
            }
        }
        let reduced = reduce_time_reversal(&full).unwrap();
        let f = run(&full);
        let r = run(&reduced);
        let de = (f.total_hartree() - r.total_hartree()).abs();
        let dn = f
            .density
            .iter()
            .zip(&r.density)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, f64::max);
        println!(
            "time reversal: full E={:.16e} Ha reduced E={:.16e} Ha delta={de:.16e} Ha density={dn:.16e} electron/Bohr^3 gap={:?} Ha",
            f.total_hartree(),
            r.total_hartree(),
            r.gap_hartree
        );
        assert!(de < 1e-10 && dn < 1e-10);
        for result in [&f, &r] {
            assert!((result.charge - 2.0).abs() < 1e-12);
            assert!((result.energy.total_hartree() - result.energy.band_form).abs() < 1e-10);
            assert_eq!(
                result.k_points.len(),
                if std::ptr::eq(result, &f) {
                    full.len()
                } else {
                    reduced.len()
                }
            );
        }
    }
    // This sampled system has overlapping bands even though integer occupations
    // converge: it must be rejected rather than silently reported as insulating.
    let metallic = monkhorst_pack([3, 1, 1], [0.0; 3]).unwrap();
    let error = ScfProblem::with_kpoints(&cell(), &atoms, &species, options, &metallic, policy())
        .unwrap()
        .run()
        .unwrap_err();
    if let ScfError::NonInsulating { gap_hartree, .. } = error {
        println!("metallic grid rejected: sampled indirect gap={gap_hartree:.16e} Ha");
        assert!(gap_hartree < 0.0);
    } else {
        panic!("unexpected error: {error}");
    }
}

#[test]
fn occupation_policy_and_weight_errors() {
    let pp = potential(2.0, 1.0, false, 401);
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let options = ScfOptions {
        cutoff_hartree: 1.5,
        ..Default::default()
    };
    assert!(
        ScfProblem::with_kpoints(
            &cell(),
            &atoms,
            &species,
            options.clone(),
            &points([0.0; 3]),
            OccupationPolicy::PartialOccupations
        )
        .is_err()
    );
    assert!(
        ScfProblem::with_kpoints(
            &cell(),
            &atoms,
            &species,
            ScfOptions {
                bands: Some(1),
                ..options.clone()
            },
            &points([0.0; 3]),
            policy()
        )
        .is_err()
    );
    assert!(
        ScfProblem::with_kpoints(
            &cell(),
            &atoms,
            &species,
            options.clone(),
            &[KPoint {
                fractional: [0.0; 3],
                weight: 0.9
            }],
            policy()
        )
        .is_err()
    );
    let problem = ScfProblem::with_kpoints(
        &cell(),
        &atoms,
        &species,
        options,
        &points([0.0; 3]),
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 100.0,
        },
    )
    .unwrap();
    assert!(matches!(problem.run(), Err(ScfError::NonInsulating { .. })));
}

#[test]
#[ignore = "nonlocal Gamma regression; cargo test --release -p dft-engine --test kpoints -- --ignored --nocapture"]
fn gamma_nonlocal_pre_p7_regression() {
    let pp = potential(4.0, 1.0, true, 1201);
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let result = ScfProblem::with_kpoints(
        &common::cell(),
        &common::atoms(),
        &species,
        ScfOptions::default(),
        &points([0.0; 3]),
        policy(),
    )
    .unwrap()
    .run()
    .unwrap();
    // ADR 0009 measured the same synthetic fixture at the default 5 Ha / 15^3 grid.
    let reference = -4.684_561_184_277_859;
    let error = (result.total_hartree() - reference).abs();
    println!(
        "Gamma nonlocal: pre-P7={reference:.16e} Ha, k API={:.16e} Ha, error={error:.16e} Ha",
        result.total_hartree()
    );
    assert!(error < 1e-10);
}
