use dft_numerics::{Fft3d, FftError, RustFftPlan3d};
use ndarray::{Array3, ShapeBuilder};
use num_complex::Complex64 as C;

fn check_backend(mut plan: impl Fft3d) {
    let shape = plan.shape();
    let source = Array3::from_shape_fn(shape, |(i, j, k)| {
        C::new((i * 7 + j * 3 + k) as f64 / 13., (i + j + k) as f64 / -17.)
    });
    let mut scratch = Array3::zeros(shape);
    for inverse in [false, true] {
        let mut expected = Array3::zeros(shape);
        if inverse {
            plan.to_real(&source, &mut expected).unwrap();
        } else {
            plan.to_reciprocal(&source, &mut expected).unwrap();
        }
        let mut actual = source.clone();
        let pointers = (actual.as_ptr(), scratch.as_ptr());
        plan.transform_in_place(&mut actual, &mut scratch, inverse)
            .unwrap();
        assert_eq!(pointers, (actual.as_ptr(), scratch.as_ptr()));
        for (a, b) in actual.iter().zip(&expected) {
            assert_eq!(
                (a.re.to_bits(), a.im.to_bits()),
                (b.re.to_bits(), b.im.to_bits())
            );
        }
        for mut bad in [Array3::zeros([1, 1, 1]), Array3::zeros((3, 4, 5).f())] {
            assert_eq!(
                plan.transform_in_place(&mut bad, &mut scratch, inverse),
                Err(FftError::IncompatibleArray)
            );
            assert_eq!(
                plan.transform_in_place(&mut actual, &mut bad, inverse),
                Err(FftError::IncompatibleArray)
            );
        }
    }
}

#[test]
fn default_in_place_preserves_bits_and_validates_dimensions() {
    check_backend(RustFftPlan3d::new([3, 4, 5]).unwrap());
}

#[cfg(feature = "mkl")]
#[test]
fn mkl_in_place_preserves_bits_and_validates_dimensions() {
    check_backend(dft_numerics::mkl_backend::MklFftPlan3d::new([3, 4, 5]).unwrap());
}

fn check_grid(mut legacy: impl Fft3d, backend: impl Fft3d) {
    use dft_numerics::{
        PairDensityOperator,
        fused::{FusedOperator, GridFusedPlan},
    };
    let shape = legacy.shape();
    let n: usize = shape.iter().product();
    let indices = [
        [0, 0, 0],
        [1, 2, 3],
        [shape[0] - 1, shape[1] - 1, shape[2] - 1],
    ];
    let c = [C::new(0.25, -0.75), C::new(-0.3, 0.2), C::new(0.6, 0.1)];
    let potential: Vec<_> = (0..n).map(|i| (i % 19) as f64 / 13. - 0.4).collect();
    let mut fused = GridFusedPlan::from_fft(backend, &indices).unwrap();
    let mut reciprocal = Array3::zeros(shape);
    let mut real = Array3::zeros(shape);
    for (index, value) in indices.iter().zip(c) {
        reciprocal[*index] = value;
    }
    legacy.to_real(&reciprocal, &mut real).unwrap();
    let mut density = vec![0.25; n];
    let mut pair_density = density.clone();
    let mut expected_density = density.clone();
    let mut expected_pair_density = density.clone();
    for ((d, p), z) in expected_density
        .iter_mut()
        .zip(&mut expected_pair_density)
        .zip(&real)
    {
        *d += 0.7 * z.norm_sqr() / 20.;
        *p += (0.3 * z.re * z.re + 0.9 * z.im * z.im) / 20.;
    }
    for (z, v) in real.iter_mut().zip(&potential) {
        *z *= *v;
    }
    legacy.to_reciprocal(&real, &mut reciprocal).unwrap();
    let mut actual = [C::default(); 3];
    // Repeated mixed calls must not depend on old scratch contents.
    for _ in 0..3 {
        fused.apply_local(&c, &potential, &mut actual).unwrap();
        for (a, index) in actual.iter().zip(indices) {
            let b = reciprocal[index];
            assert_eq!(
                (a.re.to_bits(), a.im.to_bits()),
                (b.re.to_bits(), b.im.to_bits())
            );
        }
        density.fill(0.25);
        pair_density.fill(0.25);
        fused
            .accumulate_density(&c, 0.7, 20., &mut density)
            .unwrap();
        fused
            .accumulate_density_pair(&c, [0.3, 0.9], 20., &mut pair_density)
            .unwrap();
        for (a, b) in density
            .iter()
            .zip(&expected_density)
            .chain(pair_density.iter().zip(&expected_pair_density))
        {
            assert_eq!(a.to_bits(), b.to_bits());
        }
    }
}

#[test]
fn rust_grid_preserves_legacy_local_and_density_bits() {
    check_grid(
        RustFftPlan3d::new([3, 4, 5]).unwrap(),
        RustFftPlan3d::new([3, 4, 5]).unwrap(),
    );
}

#[cfg(feature = "mkl")]
#[test]
fn mkl_grid_preserves_legacy_pruned_pattern_bits() {
    use dft_numerics::mkl_backend::MklFftPlan3d;
    for shape in [[3, 4, 5], [25; 3], [36; 3], [72, 36, 36], [7, 9, 11]] {
        check_grid(
            MklFftPlan3d::new(shape).unwrap(),
            MklFftPlan3d::new(shape).unwrap(),
        );
    }
}
