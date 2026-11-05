//! Capture-only P15A diagnostic; called only in a source-only measurement build.
use super::*;

pub(super) fn density_error_ry(
    input: &Array3<f64>,
    output: &Array3<f64>,
    cell: &Cell,
    cutoff_hartree: f64,
    fft: &mut FftPlan3d,
) -> Result<f64, ScfError> {
    let field = Array3::from_shape_fn(input.raw_dim(), |index| {
        Complex64::new(output[index] - input[index], 0.0)
    });
    let mut reciprocal = Array3::zeros(input.raw_dim());
    fft.to_reciprocal(&field, &mut reciprocal)?;
    let sphere = PlaneWaveBasis::new(cell, Hartree(cutoff_hartree))?.density_basis(cell)?;
    let mut energy = 0.0;
    for g in sphere.vectors() {
        if g.norm_squared > 0.0 {
            let index = FftGridSpec::map_miller(g.miller, fft.shape())?;
            energy += 4.0 * std::f64::consts::PI * reciprocal[index].norm_sqr() / g.norm_squared;
        }
    }
    // 2 * E_H[delta n] (Ha) = dr2 (Ry). The FFT already divides by grid size.
    Ok(cell.volume_bohr3() * energy)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn signed_mode_metric_has_analytic_units_and_sphere_filter() {
        let side = 8.0;
        let cell =
            Cell::from_user_rows_bohr([[side, 0., 0.], [0., side, 0.], [0., 0., side]]).unwrap();
        let mut fft = FftPlan3d::new([12, 12, 12]).unwrap();
        let input = Array3::zeros([12, 12, 12]);
        let amplitude = 0.03;
        let pi = std::f64::consts::PI;
        let output = Array3::from_shape_fn([12, 12, 12], |(i, _, _)| {
            0.2 + amplitude * (2.0 * pi * i as f64 / 12.0).cos()
                + 0.07 * (8.0 * pi * i as f64 / 12.0).cos()
        });
        // Density sphere energy cutoff 4*0.2 Ha includes |m|=1, excludes |m|=4.
        let measured = density_error_ry(&input, &output, &cell, 0.2, &mut fft).unwrap();
        let g2 = (2.0 * pi / side).powi(2);
        let expected = side.powi(3) * 2.0 * pi * amplitude * amplitude / g2;
        assert!((measured - expected).abs() < 1e-12);
        let reverse = density_error_ry(&output, &input, &cell, 0.2, &mut fft).unwrap();
        assert_eq!(measured, reverse);
        assert_eq!(
            density_error_ry(&output, &output, &cell, 0.2, &mut fft).unwrap(),
            0.0
        );
    }
}
