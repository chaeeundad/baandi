#[allow(dead_code)]
#[path = "../../dft-engine/tests/common/mod.rs"]
mod common;
use dft_basis::PlaneWaveBasis;
use dft_physics::{
    data::Hartree,
    hamiltonian::{Species, ionic_local, ionic_local_with_conventions},
};
use dft_upf::radial::{
    LocalMesh as L, RadialEvenRule as R, RadialMesh, projector_transform_with_rule,
    spherical_density_transform_with_rule,
};

#[test]
fn local_count_and_ionic_reconstruction() {
    let cell = common::cell();
    let basis = PlaneWaveBasis::new(&cell, Hartree(1.0)).unwrap();
    for n in [10, 11, 20, 21, 100, 101] {
        let pp = common::potential(2., 8., true, n);
        let positions = [[0.1, 0.2, 0.3]];
        let species = [Species {
            potential: &pp,
            positions_bohr: &positions,
        }];
        let default = ionic_local(&basis, cell.volume_bohr3(), &species).unwrap();
        let explicit = ionic_local_with_conventions(
            &basis,
            cell.volume_bohr3(),
            &species,
            L::Baandi,
            R::Simpson3_8,
        )
        .unwrap();
        for (a, b) in default.iter().zip(&explicit) {
            assert_eq!(a.re.to_bits(), b.re.to_bits());
            assert_eq!(a.im.to_bits(), b.im.to_bits());
        }
        let count = pp
            .mesh()
            .r()
            .iter()
            .position(|&r| r > 10.)
            .map_or(n, |i| i + 1);
        let count = if count % 2 == 0 { count - 1 } else { count };
        for mesh in [L::Baandi, L::QeLike] {
            for rule in [R::Simpson3_8, R::SimpsonQuadraticEnd] {
                let local = pp
                    .local_transform_with_conventions(cell.volume_bohr3(), mesh, rule)
                    .unwrap();
                if mesh == L::QeLike {
                    assert_eq!(local.effective_r_max(), pp.mesh().r()[count - 1]);
                }
                let ionic =
                    ionic_local_with_conventions(&basis, cell.volume_bohr3(), &species, mesh, rule)
                        .unwrap();
                for (g, v) in basis.vectors().iter().zip(ionic) {
                    let radial = if g.norm_squared == 0. {
                        local.alpha_z().unwrap()
                    } else {
                        local.at_g(g.norm_squared.sqrt()).unwrap()
                    };
                    let expected = radial
                        * dft_physics::hamiltonian::structure_factor(g.cartesian, &positions)
                            .unwrap();
                    assert!((v - expected).norm() < 1e-13 * expected.norm().max(1.));
                }
            }
        }
    }
}

#[test]
fn core_and_projector_use_even_rule_and_keep_prefix() {
    let mesh = RadialMesh::new((0..8).map(|i| i as f64).collect(), vec![1.; 8]).unwrap();
    // At q=0 the r*u and r²*rho integrands both equal r³, whose
    // quadratic-end error is +1/4. Projector prefix has only four samples.
    let beta: Vec<_> = (0..4).map(|i| (i * i) as f64).collect();
    let rho: Vec<_> = (0..8).map(|i| i as f64).collect();
    for (a, b) in [
        (
            projector_transform_with_rule(&mesh, 0, &beta, 0., R::Simpson3_8).unwrap(),
            projector_transform_with_rule(&mesh, 0, &beta, 0., R::SimpsonQuadraticEnd).unwrap(),
        ),
        (
            spherical_density_transform_with_rule(&mesh, &rho, 0., R::Simpson3_8).unwrap(),
            spherical_density_transform_with_rule(&mesh, &rho, 0., R::SimpsonQuadraticEnd).unwrap(),
        ),
    ] {
        assert!((b - a - std::f64::consts::PI).abs() < 1e-11);
    }
    let expected = 4. * std::f64::consts::PI * 81. / 4.;
    assert!(
        (projector_transform_with_rule(&mesh, 0, &beta, 0., R::Simpson3_8).unwrap() - expected)
            .abs()
            < 1e-12
    );
}
