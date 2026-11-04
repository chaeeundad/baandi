//! Keeps `libmkl_sequential` and `libmkl_core` in DT_NEEDED under the default
//! `--as-needed` linking. Each static below takes the address of one exported
//! symbol of that library, which creates an undefined-symbol reference from this
//! crate, so the linker must keep the defining library. Nothing here is ever
//! called: only the addresses are stored, signatures are irrelevant.

unsafe extern "C" {
    // Public service function exported by libmkl_core.
    fn MKL_Get_Version_String();
    // BLAS driver of the sequential threading layer (the symbol the loader reported
    // missing when the layer was dropped).
    fn mkl_blas_dgemm();
}

type Anchor = unsafe extern "C" fn();

#[used]
static KEEP_MKL_CORE: Anchor = MKL_Get_Version_String;
#[used]
static KEEP_MKL_SEQUENTIAL: Anchor = mkl_blas_dgemm;

/// Reads both anchors so the object file that holds them is always linked.
#[inline(never)]
pub(crate) fn anchor() -> usize {
    std::hint::black_box(KEEP_MKL_CORE as usize)
        ^ std::hint::black_box(KEEP_MKL_SEQUENTIAL as usize)
}
