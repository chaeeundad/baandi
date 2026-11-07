//! Full MKL transforms implementing the same support contract as the Rust path.
use crate::{Fft3d, FftError, SpherePattern};
use ndarray::Array3;
use num_complex::Complex64;

pub struct MklFftPlan3d {
    shape: [usize; 3],
    plan: dft_mkl::Plan3d,
}
impl MklFftPlan3d {
    pub fn new(shape: [usize; 3]) -> Result<Self, FftError> {
        let plan = dft_mkl::Plan3d::new(shape).map_err(|e| match e {
            dft_mkl::Error::Dimension => FftError::InvalidGrid,
            e => FftError::Mkl(e),
        })?;
        Ok(Self { shape, plan })
    }
    fn validate(
        &self,
        input: &Array3<Complex64>,
        output: &Array3<Complex64>,
    ) -> Result<(), FftError> {
        if input.shape() != self.shape
            || output.shape() != self.shape
            || !input.is_standard_layout()
            || !output.is_standard_layout()
        {
            return Err(FftError::IncompatibleArray);
        }
        Ok(())
    }
    fn transform(
        &mut self,
        input: &Array3<Complex64>,
        output: &mut Array3<Complex64>,
        inverse: bool,
    ) -> Result<(), FftError> {
        self.validate(input, output)?;
        self.plan.transform(
            input.as_slice().unwrap(),
            output.as_slice_mut().unwrap(),
            inverse,
        )?;
        Ok(())
    }
    /// Input outside reciprocal support must be canonical zero. Output is the
    /// FULL real-space field; zeroing it would violate the inverse contract.
    pub fn to_real_pruned(
        &mut self,
        input: &Array3<Complex64>,
        pattern: &SpherePattern,
        output: &mut Array3<Complex64>,
    ) -> Result<(), FftError> {
        self.validate(input, output)?;
        if pattern.shape != self.shape {
            return Err(FftError::IncompatibleArray);
        }
        debug_assert!(
            input
                .iter()
                .zip(&pattern.present)
                .all(|(z, p)| *p || (z.re.to_bits() == 0 && z.im.to_bits() == 0)),
            "inverse input outside pattern must be canonical zero"
        );
        self.to_real(input, output)
    }
    /// Full forward transform followed by exact zero outside reciprocal support.
    pub fn to_reciprocal_pruned(
        &mut self,
        input: &Array3<Complex64>,
        pattern: &SpherePattern,
        output: &mut Array3<Complex64>,
    ) -> Result<(), FftError> {
        if pattern.shape != self.shape {
            return Err(FftError::IncompatibleArray);
        }
        self.to_reciprocal(input, output)?;
        for (z, p) in output.iter_mut().zip(&pattern.present) {
            if !p {
                *z = Complex64::default();
            }
        }
        Ok(())
    }
}
impl Fft3d for MklFftPlan3d {
    fn transform_in_place(
        &mut self,
        data: &mut Array3<Complex64>,
        scratch: &mut Array3<Complex64>,
        inverse: bool,
    ) -> Result<(), FftError> {
        self.validate(data, scratch)?;
        self.plan
            .transform_in_place(data.as_slice_mut().unwrap(), inverse)?;
        Ok(())
    }
    fn shape(&self) -> [usize; 3] {
        self.shape
    }
    fn to_real(
        &mut self,
        input: &Array3<Complex64>,
        output: &mut Array3<Complex64>,
    ) -> Result<(), FftError> {
        self.transform(input, output, true)
    }
    fn to_reciprocal(
        &mut self,
        input: &Array3<Complex64>,
        output: &mut Array3<Complex64>,
    ) -> Result<(), FftError> {
        self.transform(input, output, false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::RustFftPlan3d;
    #[test]
    fn support_contract_layout_and_empty_support() {
        let shape = [3, 4, 5];
        let mut mkl = MklFftPlan3d::new(shape).unwrap();
        let mut rust = RustFftPlan3d::new(shape).unwrap();
        let mut input = Array3::zeros(shape);
        input[[0, 0, 0]] = Complex64::new(0.5, 0.2);
        input[[2, 1, 3]] = Complex64::new(-0.3, 0.4);
        let mut actual = input.clone();
        let mut expected = input.clone();
        let pattern =
            SpherePattern::from_indices(shape, &[[0, 0, 0], [2, 1, 3], [2, 1, 3]]).unwrap();
        mkl.to_real_pruned(&input, &pattern, &mut actual).unwrap();
        rust.to_real_pruned(&input, &pattern, &mut expected)
            .unwrap();
        assert!(
            actual
                .iter()
                .zip(&expected)
                .all(|(a, b)| (*a - *b).norm() < 1e-13)
        );
        let real = actual.clone();
        mkl.to_reciprocal_pruned(&real, &pattern, &mut actual)
            .unwrap();
        for ((a, b), present) in actual.iter().zip(&input).zip(&pattern.present) {
            assert!((*a - *b).norm() < 1e-13);
            if !present {
                assert_eq!((a.re.to_bits(), a.im.to_bits()), (0, 0));
            }
        }
        let empty = SpherePattern::from_indices(shape, &[]).unwrap();
        mkl.to_reciprocal_pruned(&real, &empty, &mut actual)
            .unwrap();
        assert!(
            actual
                .iter()
                .all(|z| z.re.to_bits() == 0 && z.im.to_bits() == 0)
        );
        mkl.to_real_pruned(&Array3::zeros(shape), &empty, &mut actual)
            .unwrap();
        assert!(actual.iter().all(|z| z.norm() == 0.0));
        assert_eq!(
            mkl.to_real(&Array3::zeros((3, 5, 4)), &mut actual),
            Err(FftError::IncompatibleArray)
        );
        let wrong = SpherePattern::from_indices([1, 1, 1], &[]).unwrap();
        assert_eq!(
            mkl.to_real_pruned(&input, &wrong, &mut actual),
            Err(FftError::IncompatibleArray)
        );
        use ndarray::ShapeBuilder;
        let fortran = Array3::zeros((3, 4, 5).f());
        assert_eq!(
            mkl.to_real(&fortran, &mut actual),
            Err(FftError::IncompatibleArray)
        );
        assert_eq!(
            mkl.to_real(&input, &mut Array3::zeros((3, 4, 5).f())),
            Err(FftError::IncompatibleArray)
        );
    }
}
