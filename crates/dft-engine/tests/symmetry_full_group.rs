#[allow(dead_code)]
mod common;
use dft_engine::*;
use dft_geometry::{
    Atom,
    kpoints::monkhorst_pack,
    symmetry::{detect, reduce_mesh},
};

#[test]
fn crystal_requires_enabled_symmetry_and_commensurate_grid() {
    let cell = common::cell();
    let atoms = [Atom::new("Xx", [0.; 3]), Atom::new("Xx", [0.25; 3])];
    let pp = common::potential(2., 1., false, 401);
    let refs = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let full = monkhorst_pack([4; 3], [0.5; 3]).unwrap();
    let ops = detect(&cell, &atoms, 1e-7).unwrap();
    let map = reduce_mesh(&full, &ops, true, false).unwrap();
    assert_eq!(ops.len(), 48);
    assert_eq!(map.operations.iter().filter(|o| !o.antiunitary).count(), 12);
    for dims in [[25; 3], [32; 3]] {
        let mut p = ScfProblem::with_kpoints(
            &cell,
            &atoms,
            &refs,
            ScfOptions {
                cutoff_hartree: 3.,
                grid_override: Some(dims),
                crystal_density_group: true,
                ..Default::default()
            },
            &map.points,
            OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-7,
            },
        )
        .unwrap();
        let mut n = p.uniform_density();
        assert!(
            p.step_kpoints(&n, None, 1e-8)
                .err()
                .unwrap()
                .to_string()
                .contains("requires enabled auto")
        );
        let enabled = p.enable_symmetry(&full, 1e-7, true);
        if dims == [25; 3] {
            assert!(enabled.err().unwrap().to_string().contains("fft_grid"));
        } else {
            enabled.unwrap();
            // A mesh-invariant Fourier field that is not crystal-invariant.
            let basis = dft_basis::PlaneWaveBasis::new(&cell, dft_units::Hartree(12.)).unwrap();
            let grid = dft_basis::FftGridSpec::new(&basis, &basis, Some(dims)).unwrap();
            let mesh_ops: Vec<_> = map
                .operations
                .iter()
                .filter(|o| !o.antiunitary)
                .cloned()
                .collect();
            for ((i, j, k), v) in n.indexed_iter_mut() {
                *v += 0.001 * (std::f64::consts::TAU * (i + j + 2 * k) as f64 / 32.).cos();
            }
            let mut fft = dft_numerics::FftPlan3d::new(dims).unwrap();
            let mut mesh =
                dft_physics::symmetry::DensitySymmetrizer::new(&basis, &grid, &mesh_ops).unwrap();
            mesh.apply(&mut n, &mut fft).unwrap();
            assert!(mesh.check_invariant(&n, &mut fft).unwrap());
            assert!(
                p.step_kpoints(&n, None, 1e-8)
                    .err()
                    .unwrap()
                    .to_string()
                    .contains("input breaks")
            );
        }
    }
}

#[test]
#[ignore = "external DFT_SI_UPF; run cargo test --release -p dft-engine --test symmetry_full_group -- --ignored --nocapture"]
fn silicon_full_group_oracles() {
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF")).unwrap(),
    )
    .unwrap();
    let refs = [SpeciesPotential {
        label: "Si",
        potential: &pp,
    }];
    let cell = common::cell();
    let atoms = [Atom::new("Si", [0.; 3]), Atom::new("Si", [0.25; 3])];
    let ops = detect(&cell, &atoms, 1e-7).unwrap();
    for shift in [0., 0.5] {
        let full = monkhorst_pack([4; 3], [shift; 3]).unwrap();
        let map = reduce_mesh(&full, &ops, true, false).unwrap();
        let run = |crystal| {
            let mut p = ScfProblem::with_kpoints(
                &cell,
                &atoms,
                &refs,
                ScfOptions {
                    cutoff_hartree: 15.,
                    grid_override: Some([32; 3]),
                    crystal_density_group: crystal,
                    density_tolerance: 1e-11,
                    energy_tolerance: 1e-12,
                    eigen_tolerance: 1e-11,
                    max_iterations: 100,
                    ..Default::default()
                },
                &map.points,
                OccupationPolicy::Insulator {
                    minimum_gap_hartree: 1e-7,
                },
            )
            .unwrap();
            p.enable_symmetry(&full, 1e-7, true).unwrap();
            p.run().unwrap()
        };
        let mesh = run(false);
        let crystal = run(true);
        let mut outside = Vec::new();
        let mut max_crystal = 0.0_f64;
        for op in &ops {
            let mut err_mesh = 0.0_f64;
            let mut err_crystal = 0.0_f64;
            // Independent real-space permutation x -> R*x+tau, no FFT projector.
            for i in 0..32 {
                for j in 0..32 {
                    for k in 0..32 {
                        let x = [i, j, k];
                        let y: [usize; 3] = std::array::from_fn(|a| {
                            let v = (0..3).map(|b| op.rotation[a][b] * x[b] as i32).sum::<i32>()
                                + (32. * op.translation[a]).round() as i32;
                            v.rem_euclid(32) as usize
                        });
                        err_mesh = err_mesh.max((mesh.density[x] - mesh.density[y]).abs());
                        err_crystal =
                            err_crystal.max((crystal.density[x] - crystal.density[y]).abs());
                    }
                }
            }
            max_crystal = max_crystal.max(err_crystal);
            assert!(err_crystal < 1e-13, "{err_crystal:e}");
            if !map.operations.iter().any(|o| {
                !o.antiunitary && o.rotation == op.rotation && o.translation == op.translation
            }) {
                outside.push(err_mesh);
                assert!(err_mesh > 1e-10);
            }
        }
        let de = crystal.energy.total_hartree() - mesh.energy.total_hartree();
        let dn = mesh
            .density
            .iter()
            .zip(&crystal.density)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        println!(
            "shift={shift} delta_E={de:.17e} Ha max_density_difference={dn:e} electrons/Bohr^3 crystal_invariance={max_crystal:e} nonmesh_errors={outside:?}"
        );
        if shift == 0. {
            assert!(de.abs() < 1e-12 && dn < 1e-12);
        } else {
            assert_eq!(outside.len(), 36);
        }
    }
}
