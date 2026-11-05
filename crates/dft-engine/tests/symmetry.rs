#[allow(dead_code)]
mod common;
use dft_engine::*;
use dft_geometry::{
    Atom, Cell,
    kpoints::monkhorst_pack,
    symmetry::{detect, reduce_mesh},
};
use dft_upf::parse_upf;

fn options() -> ScfOptions {
    ScfOptions {
        cutoff_hartree: 3.,
        density_tolerance: 1e-11,
        energy_tolerance: 1e-12,
        eigen_tolerance: 1e-11,
        max_iterations: 250,
        ..Default::default()
    }
}
fn compare(
    cell: &Cell,
    atoms: &[Atom],
    potentials: &[SpeciesPotential<'_>],
    grid: [usize; 3],
    shift: [f64; 3],
    mut opts: ScfOptions,
) {
    let full = monkhorst_pack(grid, shift).unwrap();
    let map = reduce_mesh(&full, &detect(cell, atoms, 1e-7).unwrap(), true, false).unwrap();
    let basis =
        dft_basis::PlaneWaveBasis::new(cell, dft_physics::data::Hartree(opts.cutoff_hartree))
            .unwrap();
    let density = basis.density_basis(cell).unwrap();
    let grid = dft_basis::FftGridSpec::new(&basis, &density, None).unwrap();
    opts.grid_override =
        Some(dft_geometry::symmetry::compatible_fft_grid(grid.dims(), &map.operations).unwrap());
    let run = |reduce, threads| {
        let mut opts = opts.clone();
        opts.threads = threads;
        let mut p = ScfProblem::with_kpoints(
            cell,
            atoms,
            potentials,
            opts,
            if reduce { &map.points } else { &full },
            OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-7,
            },
        )
        .unwrap();
        if p.core_density().is_some() {
            p.set_density_floor(1e-10).unwrap();
        }
        if reduce {
            p.enable_symmetry(&full, 1e-7, true).unwrap();
        }
        p.run().unwrap()
    };
    let a = run(false, 1);
    let b = run(true, 1);
    let c = run(true, 4);
    println!("FFT shape={:?}", a.density.shape());
    let de = (a.energy.total_hartree() - b.energy.total_hartree()).abs();
    let dn = density_distance(&a.density, &b.density, cell.volume_bohr3(), a.charge);
    println!(
        "atoms={} full={} irreducible={} dE={de:.16e} Ha dn={dn:.16e} Bohr^-3/2",
        atoms.len(),
        full.len(),
        map.points.len()
    );
    assert!(de <= 1e-10);
    assert!(dn <= 1e-9);
    assert_eq!(
        b.energy.total_hartree().to_bits(),
        c.energy.total_hartree().to_bits()
    );
    assert!(
        b.density
            .iter()
            .zip(&c.density)
            .all(|(a, b)| a.to_bits() == b.to_bits())
    );
}
#[test]
#[ignore = "release-mode SCF matrix: cargo test --release -p dft-engine --test symmetry -- --ignored --nocapture"]
fn synthetic_shapes() {
    let pp = common::potential(2., 1., true, 401);
    let refs = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    for rows in [
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]],
        [[4., 0., 0.], [0., 4., 0.], [0., 0., 12.]],
        [[5., 0.2, 0.1], [0., 5.3, 0.4], [0.3, 0., 5.7]],
    ] {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        compare(
            &cell,
            &[Atom::new("Xx", [0.; 3])],
            &refs,
            [2; 3],
            [0.5; 3],
            options(),
        );
    }
}
#[test]
#[ignore = "external Si UPF via DFT_SI_UPF; release mode"]
fn diamond_and_cubic_si() {
    let path = std::env::var("DFT_SI_UPF").expect("set DFT_SI_UPF to external Si.pz-vbc.UPF");
    let pp = parse_upf(&std::fs::read_to_string(path).unwrap()).unwrap();
    let refs = [SpeciesPotential {
        label: "Si",
        potential: &pp,
    }];
    let cell = common::cell();
    let atoms = [Atom::new("Si", [0.; 3]), Atom::new("Si", [0.25; 3])];
    let opts = ScfOptions {
        cutoff_hartree: 15.,
        ..options()
    };
    for shift in [0., 0.5] {
        compare(&cell, &atoms, &refs, [4; 3], [shift; 3], opts.clone());
    }
    let a = 10.261212856717933;
    let cell = Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap();
    let mut atoms = Vec::new();
    for f in [[0., 0., 0.], [0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]] {
        for d in [0., 0.25] {
            atoms.push(Atom::new("Si", f.map(|v| v + d)));
        }
    }
    compare(&cell, &atoms, &refs, [2; 3], [0.5; 3], opts);
}
#[test]
#[ignore = "external PseudoDojo Mg/O via DFT_NC_DIR; release mode"]
fn rocksalt_mgo() {
    let root = std::path::PathBuf::from(std::env::var("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let mg = parse_upf(&std::fs::read_to_string(root.join("Mg.upf")).unwrap()).unwrap();
    let o = parse_upf(&std::fs::read_to_string(root.join("O.upf")).unwrap()).unwrap();
    let a = 7.96;
    let cell = Cell::from_user_rows_bohr([
        [0., a / 2., a / 2.],
        [a / 2., 0., a / 2.],
        [a / 2., a / 2., 0.],
    ])
    .unwrap();
    let atoms = [Atom::new("Mg", [0.; 3]), Atom::new("O", [0.5; 3])];
    compare(
        &cell,
        &atoms,
        &[
            SpeciesPotential {
                label: "Mg",
                potential: &mg,
            },
            SpeciesPotential {
                label: "O",
                potential: &o,
            },
        ],
        [2; 3],
        [0.5; 3],
        ScfOptions {
            functional: Functional::Pbe,
            cutoff_hartree: 20.,
            ..options()
        },
    );
}
#[test]
#[ignore = "release-mode spin and finite-T SCF matrix"]
fn spin_and_fermi_dirac() {
    use dft_engine::{
        collinear::{CollinearProblem, SpinSector},
        finite_temperature::FiniteTProblem,
    };
    use dft_physics::spin::SpinDensity;
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.; 3])];
    let full = monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    let ops = detect(&cell, &atoms, 1e-7).unwrap();
    let map = reduce_mesh(&full, &ops, false, false).unwrap();
    let pp = common::potential_with_functional(1., 1., false, 401, "SLA PW");
    let refs = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let run = |reduced| {
        let mut p = CollinearProblem::new(
            &cell,
            &atoms,
            &refs,
            ScfOptions {
                functional: Functional::LdaPw,
                bands: Some(4),
                ..options()
            },
            if reduced { &map.points } else { &full },
            SpinSector {
                n_up: 1,
                n_down: 0,
                axis_cartesian: [0., 0., 1.],
                initial_moment_seed: None,
            },
        )
        .unwrap();
        if reduced {
            p.enable_symmetry(&full, 1e-7).unwrap();
        }
        p.run().unwrap()
    };
    let a = run(false);
    let b = run(true);
    assert!(a.converged && b.converged);
    let SpinDensity::Collinear { channels: ac, .. } = a.density else {
        panic!()
    };
    let SpinDensity::Collinear { channels: bc, .. } = b.density else {
        panic!()
    };
    let de = (a.energy.total_hartree() - b.energy.total_hartree()).abs();
    let dn = (0..2)
        .map(|s| density_distance(&ac[s], &bc[s], 125., 1.))
        .fold(0., f64::max);
    let dm = (a.moment_mu_b - b.moment_mu_b).abs();
    println!(
        "spin: full={} reduced={} dE={de:e} Ha dn={dn:e} Bohr^-3/2 dM={dm:e} mu_B",
        full.len(),
        map.points.len()
    );
    assert!(de <= 1e-10 && dn <= 1e-9 && dm <= 1e-9);
    let pp = common::potential(1., 1., false, 401);
    let refs = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let run = |reduced| {
        let mut p = FiniteTProblem::new(
            &cell,
            &atoms,
            &refs,
            ScfOptions {
                cutoff_hartree: 4.5,
                ..options()
            },
            if reduced { &map.points } else { &full },
            0.01,
            12,
        )
        .unwrap();
        if reduced {
            p.enable_symmetry(&full, 1e-7, true).unwrap();
        }
        p.run().unwrap()
    };
    let a = run(false);
    let b = run(true);
    assert!(a.converged && b.converged);
    let de = (a.thermodynamics.internal_energy_ha - b.thermodynamics.internal_energy_ha).abs();
    let df = (a.thermodynamics.free_energy_ha - b.thermodynamics.free_energy_ha).abs();
    let dn = density_distance(&a.state.density, &b.state.density, 125., 1.);
    println!("FD: dE={de:e} Ha dF={df:e} Ha dn={dn:e} Bohr^-3/2");
    assert!(de <= 1e-10 && df <= 1e-10 && dn <= 1e-9);
}
#[test]
fn rejects_noninvariant_density_and_mismatched_mesh() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.; 3])];
    let full = monkhorst_pack([2; 3], [0.5; 3]).unwrap();
    let map = reduce_mesh(&full, &detect(&cell, &atoms, 1e-7).unwrap(), true, false).unwrap();
    let pp = common::potential(2., 1., false, 401);
    let refs = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let mut p = ScfProblem::with_kpoints(
        &cell,
        &atoms,
        &refs,
        options(),
        &map.points,
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap();
    assert!(
        p.enable_symmetry(&monkhorst_pack([3; 3], [0.; 3]).unwrap(), 1e-7, true)
            .is_err()
    );
    p.enable_symmetry(&full, 1e-7, true).unwrap();
    let mut n = p.uniform_density();
    n[[0, 0, 1]] += 1e-4;
    assert!(p.step_kpoints(&n, None, 1e-8).is_err());
}
