use num_complex::Complex64;
use std::{
    cell::Cell,
    ffi::{CStr, c_char, c_int, c_long, c_void},
    marker::PhantomData,
    ptr,
};

// MKL_LONG is C long; MKL_INT is int for the explicitly linked LP64 interface.
// Enum arguments use C int, including values in variadic argument positions.
unsafe extern "C" {
    fn DftiCreateDescriptor(
        handle: *mut *mut c_void,
        precision: c_int,
        domain: c_int,
        rank: c_long,
        ...
    ) -> c_long;
    fn DftiSetValue(handle: *mut c_void, parameter: c_int, ...) -> c_long;
    fn DftiCommitDescriptor(handle: *mut c_void) -> c_long;
    fn DftiComputeForward(handle: *mut c_void, data: *mut c_void, ...) -> c_long;
    fn DftiComputeBackward(handle: *mut c_void, data: *mut c_void, ...) -> c_long;
    fn DftiFreeDescriptor(handle: *mut *mut c_void) -> c_long;
    fn DftiErrorMessage(status: c_long) -> *const c_char;
    fn DftiErrorClass(status: c_long, class: c_long) -> c_long;
    fn cblas_zgemm(
        layout: c_int,
        ta: c_int,
        tb: c_int,
        m: c_int,
        n: c_int,
        k: c_int,
        alpha: *const c_void,
        a: *const c_void,
        lda: c_int,
        b: *const c_void,
        ldb: c_int,
        beta: *const c_void,
        c: *mut c_void,
        ldc: c_int,
    );
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Error {
    Dimension,
    Dfti(c_long),
}
impl std::fmt::Display for Error {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::Dimension => write!(f, "invalid MKL dimensions or buffer lengths"),
            Self::Dfti(status) => {
                // SAFETY: DFTI returns a static NUL-terminated diagnostic for a status.
                let message = unsafe { DftiErrorMessage(*status) };
                if message.is_null() {
                    return write!(f, "DFTI status {status}");
                }
                // SAFETY: non-null library diagnostic has static lifetime and a NUL terminator.
                let message = unsafe { CStr::from_ptr(message) };
                write!(f, "DFTI status {status}: {}", message.to_string_lossy())
            }
        }
    }
}
impl std::error::Error for Error {}
fn status(code: c_long) -> Result<(), Error> {
    // SAFETY: integer-only status query; zero is DFTI_NO_ERROR.
    if unsafe { DftiErrorClass(code, 0) } != 0 {
        Ok(())
    } else {
        Err(Error::Dfti(code))
    }
}

/// Owns one descriptor and one C-order complex buffer. Transform calls allocate
/// no Rust storage. Send but not Sync: move to a worker, never share a descriptor
/// between simultaneous calls. The linked sequential layer owns no worker pool.
pub struct Plan3d {
    handle: *mut c_void,
    buffer: Vec<Complex64>,
    not_sync: PhantomData<Cell<()>>,
}
// SAFETY: descriptor and buffer have unique ownership, no thread affinity; all
// operations require &mut self. Moving cannot cause concurrent DFTI access.
unsafe impl Send for Plan3d {}
impl Plan3d {
    pub fn new(shape: [usize; 3]) -> Result<Self, Error> {
        let len = shape
            .iter()
            .try_fold(1usize, |n, &d| n.checked_mul(d))
            .filter(|&n| n > 0 && n <= isize::MAX as usize / size_of::<Complex64>())
            .ok_or(Error::Dimension)?;
        let mut dims = [0 as c_long; 3];
        for (dst, src) in dims.iter_mut().zip(shape) {
            *dst = src.try_into().map_err(|_| Error::Dimension)?;
        }
        let mut plan = Self {
            handle: ptr::null_mut(),
            buffer: vec![Complex64::default(); len],
            not_sync: PhantomData,
        };
        // SAFETY: writable handle, rank=3 with three C-long lengths; DOUBLE=36,
        // COMPLEX=32. Dimensions and byte count validated above. Drop frees any
        // descriptor created even when a later configuration call fails.
        unsafe {
            status(DftiCreateDescriptor(
                &mut plan.handle,
                36,
                32,
                3,
                dims.as_ptr(),
            ))?;
            // PLACEMENT=11, INPLACE=43; FORWARD_SCALE=4 accepts double.
            status(DftiSetValue(plan.handle, 11, 43 as c_int))?;
            status(DftiSetValue(plan.handle, 4, 1.0 / len as f64))?;
            status(DftiCommitDescriptor(plan.handle))?;
        }
        Ok(plan)
    }
    /// Input/output lengths must equal the plan size. Safe borrows prohibit
    /// overlap; Complex64 has the C layout of two consecutive f64 values.
    pub fn transform(
        &mut self,
        input: &[Complex64],
        output: &mut [Complex64],
        inverse: bool,
    ) -> Result<(), Error> {
        if input.len() != self.buffer.len() || output.len() != self.buffer.len() {
            return Err(Error::Dimension);
        }
        self.buffer.copy_from_slice(input);
        // SAFETY: committed in-place double-complex descriptor and uniquely
        // borrowed contiguous buffer of exactly its validated size.
        unsafe {
            status(if inverse {
                DftiComputeBackward(self.handle, self.buffer.as_mut_ptr().cast())
            } else {
                DftiComputeForward(self.handle, self.buffer.as_mut_ptr().cast())
            })?;
        }
        output.copy_from_slice(&self.buffer);
        Ok(())
    }
}
impl Drop for Plan3d {
    fn drop(&mut self) {
        if !self.handle.is_null() {
            // SAFETY: descriptor owned exclusively here; freed exactly once.
            // Drop cannot propagate errors; all computational errors are returned.
            unsafe {
                DftiFreeDescriptor(&mut self.handle);
            }
        }
    }
}

/// Borrowed contiguous column-major matrix, with checked storage dimensions.
#[derive(Clone, Copy)]
pub struct Matrix<'a> {
    data: &'a [Complex64],
    rows: usize,
    cols: usize,
}
impl<'a> Matrix<'a> {
    pub fn new(data: &'a [Complex64], rows: usize, cols: usize) -> Result<Self, Error> {
        if rows.checked_mul(cols) != Some(data.len())
            || rows > c_int::MAX as usize
            || cols > c_int::MAX as usize
        {
            return Err(Error::Dimension);
        }
        Ok(Self { data, rows, cols })
    }
}
/// Computes op(A) op(B) into C, without allocation. `adjoint` selects conjugate
/// transpose. Empty products return zeros without passing dangling pointers to C.
pub fn gemm(
    a: Matrix<'_>,
    adjoint_a: bool,
    b: Matrix<'_>,
    adjoint_b: bool,
    c: &mut [Complex64],
) -> Result<(), Error> {
    let (m, k) = if adjoint_a {
        (a.cols, a.rows)
    } else {
        (a.rows, a.cols)
    };
    let (l, n) = if adjoint_b {
        (b.cols, b.rows)
    } else {
        (b.rows, b.cols)
    };
    if k != l || m.checked_mul(n) != Some(c.len()) {
        return Err(Error::Dimension);
    }
    if m == 0 || n == 0 || k == 0 {
        c.fill(Complex64::default());
        return Ok(());
    }
    let one = Complex64::new(1.0, 0.0);
    let zero = Complex64::default();
    // SAFETY: checked LP64 dimensions, exact contiguous column-major lengths,
    // positive leading dimensions, interleaved Complex64 C layout; immutable
    // inputs cannot overlap uniquely borrowed output. Scalars live through call.
    unsafe {
        cblas_zgemm(
            102,
            if adjoint_a { 113 } else { 111 },
            if adjoint_b { 113 } else { 111 },
            m as c_int,
            n as c_int,
            k as c_int,
            (&one as *const Complex64).cast(),
            a.data.as_ptr().cast(),
            a.rows as c_int,
            b.data.as_ptr().cast(),
            b.rows as c_int,
            (&zero as *const Complex64).cast(),
            c.as_mut_ptr().cast(),
            m as c_int,
        );
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn dfti_error_is_typed_and_has_diagnostic() {
        let plan = Plan3d::new([2, 3, 4]).unwrap();
        // SAFETY: valid descriptor, PLACEMENT expects an enum promoted to int.
        // Deliberately invalid enum value tests DFTI's configuration error path;
        // no array pointer or inconsistent buffer is ever passed to computation.
        let code = unsafe { DftiSetValue(plan.handle, 11, 999 as c_int) };
        let error = status(code).unwrap_err();
        assert!(matches!(error, Error::Dfti(_)));
        assert!(error.to_string().contains("DFTI status"));
    }
}
