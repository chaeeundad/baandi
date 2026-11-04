//! Keeps `libmkl_sequential` and `libmkl_core` in DT_NEEDED under the default
//! `--as-needed` linking. Each static below takes the address of one exported
//! symbol of that library, which creates an undefined-symbol reference from this
//! crate, so the linker must keep the defining library. Nothing here is ever
//! called: only the addresses are stored, signatures are irrelevant.

unsafe extern "C" {
    // Exported by libmkl_core only (checked with `nm -D` on MKL 2026.1: absent from
    // libmkl_intel_lp64 and libmkl_sequential). Not public API: if a future MKL
    // renames it the link fails loudly with an undefined symbol, never silently.
    fn mkl_serv_allocate();
    // BLAS driver of the sequential threading layer. Also defined in libmkl_core; the
    // linker binds the first library listed (mkl_sequential), which is the one we need.
    fn mkl_blas_dgemm();
}

type Anchor = unsafe extern "C" fn();

#[used]
static KEEP_MKL_CORE: Anchor = mkl_serv_allocate;
#[used]
static KEEP_MKL_SEQUENTIAL: Anchor = mkl_blas_dgemm;

/// Reads both anchors so the object file that holds them is always linked.
#[inline(never)]
pub(crate) fn anchor() -> usize {
    std::hint::black_box(KEEP_MKL_CORE as usize)
        ^ std::hint::black_box(KEEP_MKL_SEQUENTIAL as usize)
}
