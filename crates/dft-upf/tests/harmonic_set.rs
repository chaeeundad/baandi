//! GK18: `real_spherical_harmonics_upto` equals `real_spherical_harmonic` bit
//! for bit for every (l, m) up to 3, on generic, axis, near-pole and
//! large/small-magnitude directions, and reports the same errors.
use dft_upf::radial::{real_spherical_harmonic, real_spherical_harmonics_upto};

#[test]
fn harmonic_set_is_bit_identical() {
    let mut dirs: Vec<[f64; 3]> = vec![
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
        [0., 0., -1.],
        [-1., -1e-300, 0.],
        [1e-17, 0., 1.],
        [3e150, -2e150, 1e150],
        [1e-200, 2e-200, -3e-200],
        [-0.0, 0.0, 2.5],
    ];
    let mut x = 0.123_456_789_f64;
    for _ in 0..20000 {
        let mut d = [0.; 3];
        for v in &mut d {
            x = (x * 3.987_654_321 + 0.271_828).fract();
            *v = (x - 0.5) * 17.;
        }
        dirs.push(d);
    }
    let mut out = [0.; 16];
    for d in dirs {
        for lmax in 0..=3usize {
            real_spherical_harmonics_upto(lmax, d, &mut out).unwrap();
            for l in 0..=lmax {
                for m in -(l as isize)..=l as isize {
                    let single = real_spherical_harmonic(l, m, d).unwrap();
                    assert_eq!(
                        single.to_bits(),
                        out[(l * l + l).wrapping_add_signed(m)].to_bits(),
                        "l={l} m={m} d={d:?}"
                    );
                }
            }
        }
    }
    assert!(real_spherical_harmonics_upto(2, [0.; 3], &mut out).is_err());
    assert!(real_spherical_harmonic(2, 0, [0.; 3]).is_err());
    assert!(real_spherical_harmonics_upto(4, [1., 0., 0.], &mut [0.; 25]).is_err());
    assert!(real_spherical_harmonics_upto(2, [f64::NAN, 0., 0.], &mut out).is_err());
    assert!(real_spherical_harmonics_upto(3, [1., 0., 0.], &mut [0.; 15]).is_err());
}
