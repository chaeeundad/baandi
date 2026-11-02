//! Additive SCF helper contracts; existing P2B tests remain unchanged.
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_physics::{
    data::*,
    hamiltonian::{Hamiltonian, NonlocalOperator},
};
use dft_solvers::LinearOperator;
#[test]
fn potential_update_is_transactional_and_matches_fresh_operator() {
    let c = Cell::from_user_rows_bohr([[4.0, 0.0, 0.0], [0.0, 4.0, 0.0], [0.0, 0.0, 4.0]]).unwrap();
    let b = PlaneWaveBasis::new(&c, Hartree(1.5)).unwrap();
    let g = FftGridSpec::new(&b, &b.density_basis(&c).unwrap(), None).unwrap();
    let nl = || NonlocalOperator::new(&b, c.volume_bohr3(), &[]).unwrap();
    let v = Array3::from_shape_fn(g.dims(), |(i, _, _)| {
        0.1 * (std::f64::consts::TAU * i as f64 / g.dims()[0] as f64).cos()
    });
    let mut h = Hamiltonian::new(&c, &b, Array3::zeros(g.dims()), nl()).unwrap();
    h.set_local_potential(v.clone()).unwrap();
    let fresh = Hamiltonian::new(&c, &b, v, nl()).unwrap();
    let x = vec![Complex64::new(0.3, 0.2); b.vectors().len()];
    let mut a = vec![Complex64::default(); x.len()];
    let mut r = a.clone();
    h.apply(&x, &mut a).unwrap();
    fresh.apply(&x, &mut r).unwrap();
    assert_eq!(a, r);
    assert!(h.set_local_potential(Array3::zeros([1; 3])).is_err());
    assert!(
        h.set_local_potential(Array3::from_elem(g.dims(), f64::NAN))
            .is_err()
    );
    h.apply(&x, &mut r).unwrap();
    assert_eq!(a, r);
}
