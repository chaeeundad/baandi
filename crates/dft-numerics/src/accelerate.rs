//! Private, fixed-size ABI probe, not a production numerical API.

#[cfg(test)]
#[link(name = "Accelerate", kind = "framework")]
unsafe extern "C" {
    fn cblas_ddot(n: i32, x: *const f64, inc_x: i32, y: *const f64, inc_y: i32) -> f64;
    fn dsyev_(
        jobz: *mut std::ffi::c_char,
        uplo: *mut std::ffi::c_char,
        n: *mut i32,
        a: *mut f64,
        lda: *mut i32,
        w: *mut f64,
        work: *mut f64,
        lwork: *mut i32,
        info: *mut i32,
    ) -> i32;
}

// The probe is test-only so no unused production API is introduced in P0.
#[cfg(test)]
mod tests {
    use super::{cblas_ddot, dsyev_};

    fn dot3(x: &[f64; 3], y: &[f64; 3]) -> f64 {
        // SAFETY: LP64 CBLAS uses 32-bit integers. Both arrays contain three
        // initialized f64 values, remain alive for the call, and use stride 1.
        unsafe { cblas_ddot(3, x.as_ptr(), 1, y.as_ptr(), 1) }
    }

    fn eigen2(mut a: [f64; 4]) -> Result<([f64; 2], [f64; 4]), i32> {
        let mut jobz = b'V' as std::ffi::c_char;
        let mut uplo = b'U' as std::ffi::c_char;
        let mut n = 2_i32;
        let mut lda = 2_i32;
        let mut w = [0.0; 2];
        let mut work = [0.0; 5];
        let mut lwork = 5_i32;
        let mut info = 0_i32;
        // SAFETY: Matches the installed legacy CLAPACK C declaration (including
        // int return and no hidden character lengths). Scalars are LP64 i32.
        // A is writable column-major 2x2 storage, lda=2; only U is read.
        // W has n entries; WORK has max(1,3*n-1)=5 entries. All pointers are
        // valid, aligned, non-overlapping, and live throughout the call.
        unsafe {
            dsyev_(
                &mut jobz,
                &mut uplo,
                &mut n,
                a.as_mut_ptr(),
                &mut lda,
                w.as_mut_ptr(),
                work.as_mut_ptr(),
                &mut lwork,
                &mut info,
            );
        }
        if info == 0 { Ok((w, a)) } else { Err(info) }
    }

    #[test]
    fn accelerate_ddot() {
        let actual = dot3(&[1.0, 2.0, 3.0], &[4.0, 5.0, 6.0]);
        println!("CBLAS ddot = {actual:.16}; expected = 32 (dimensionless)");
        assert_eq!(actual, 32.0);
    }

    #[test]
    fn accelerate_dsyev() {
        // Symmetric (therefore Hermitian) matrix [[2,1],[1,2]]. Poison the
        // unused lower triangle to verify UPLO and column-major interpretation.
        let (values, vectors) = eigen2([2.0, f64::NAN, 1.0, 2.0]).expect("DSYEV INFO=0");
        let mut residual = 0.0_f64;
        let mut orthogonality = 0.0_f64;
        for j in 0..2 {
            assert!((values[j] - [1.0, 3.0][j]).abs() < 1e-12);
            let x = vectors[2 * j];
            let y = vectors[2 * j + 1];
            residual = residual.max((2.0 * x + y - values[j] * x).abs());
            residual = residual.max((x + 2.0 * y - values[j] * y).abs());
            for k in 0..2 {
                let dot = x * vectors[2 * k] + y * vectors[2 * k + 1];
                orthogonality = orthogonality.max((dot - f64::from(j == k)).abs());
            }
        }
        println!(
            "LAPACK dsyev INFO=0; eigenvalues={values:?}; max residual={residual:.16e}; max orthogonality error={orthogonality:.16e} (all dimensionless)"
        );
        assert!(residual < 1e-12);
        assert!(orthogonality < 1e-12);
    }
}
