//! Thin FFI adapter for the simple-dftd3 C API (`dftd3.h`, library version
//! 1.x), opened at run time. Only the functions used by this crate are
//! declared. The declarations follow the public header:
//!
//! - handles (`dftd3_error`, `dftd3_structure`, `dftd3_model`, `dftd3_param`)
//!   are opaque pointers, created by `dftd3_new_*`/`dftd3_load_*` and released by
//!   `dftd3_delete_*`, which take a pointer to the handle;
//! - every operation reports failure through the error handle, queried with
//!   `dftd3_check_error` (nonzero means error) and `dftd3_get_error`;
//! - arrays are C row-major: positions `[n][3]`, lattice `[3][3]` with one
//!   lattice vector per row, gradient `[n][3]`, virial `[3][3]`, all in Bohr and
//!   Hartree.
//!
//! Preconditions common to every unsafe call below: the function pointer was
//! resolved from a library that stays loaded for the whole process (the
//! `Library` is leaked on purpose, so the pointer never dangles), its signature
//! matches the header, handles passed in are live (created here and not yet
//! deleted), and array pointers point to buffers of the documented length.
//! Calls are serialized by `CALL_LOCK` because the library keeps no
//! documented thread-safety guarantee for concurrent calls.
#![allow(unsafe_code)]

use std::ffi::{CString, c_char, c_int, c_void};
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};

use crate::{Damping, Error, Parameters, Structure};

type Handle = *mut c_void;
type GetVersion = unsafe extern "C" fn() -> c_int;
type NewError = unsafe extern "C" fn() -> Handle;
type CheckError = unsafe extern "C" fn(Handle) -> c_int;
type GetError = unsafe extern "C" fn(Handle, *mut c_char, *const c_int);
type Delete = unsafe extern "C" fn(*mut Handle);
type NewStructure = unsafe extern "C" fn(
    Handle,
    c_int,
    *const c_int,
    *const f64,
    *const f64,
    *const bool,
) -> Handle;
type NewModel = unsafe extern "C" fn(Handle, Handle) -> Handle;
type LoadDamping = unsafe extern "C" fn(Handle, *mut c_char, bool) -> Handle;
type NewDamping = unsafe extern "C" fn(Handle, f64, f64, f64, f64, f64, f64) -> Handle;
type GetDispersion =
    unsafe extern "C" fn(Handle, Handle, Handle, Handle, *mut f64, *mut f64, *mut f64);

/// Energy (Ha), gradient (Ha/Bohr, per atom) and virial (Ha).
type Evaluation = (f64, Vec<[f64; 3]>, [[f64; 3]; 3]);

/// Resolved entry points of one loaded library.
pub(crate) struct Api {
    get_version: GetVersion,
    new_error: NewError,
    check_error: CheckError,
    get_error: GetError,
    delete_error: Delete,
    new_structure: NewStructure,
    delete_structure: Delete,
    new_d3_model: NewModel,
    delete_model: Delete,
    load_rational_damping: LoadDamping,
    load_zero_damping: LoadDamping,
    new_rational_damping: NewDamping,
    new_zero_damping: NewDamping,
    delete_param: Delete,
    get_dispersion: GetDispersion,
}

static LIBRARIES: OnceLock<Mutex<Vec<(PathBuf, &'static Api)>>> = OnceLock::new();
static CALL_LOCK: Mutex<()> = Mutex::new(());

fn symbol<T: Copy>(lib: &'static libloading::Library, name: &str) -> Result<T, Error> {
    // SAFETY: `T` is one of the function pointer types above, each matching
    // the C declaration of `name` in dftd3.h; the library is leaked, so the
    // copied pointer stays valid for the rest of the process.
    unsafe {
        lib.get::<T>(name.as_bytes())
            .map(|s| *s)
            .map_err(|e| Error::Load(format!("symbol {name}: {e}")))
    }
}

impl Api {
    /// Open (once per path and process) and resolve all entry points.
    pub(crate) fn open(path: &Path) -> Result<&'static Api, Error> {
        let registry = LIBRARIES.get_or_init(|| Mutex::new(Vec::new()));
        let mut list = registry
            .lock()
            .map_err(|_| Error::Load("library registry poisoned".into()))?;
        if let Some((_, api)) = list.iter().find(|(p, _)| p == path) {
            return Ok(api);
        }
        // SAFETY: loading runs the library's initializers (gfortran runtime,
        // OpenMP); this is the documented way to use the shared library and
        // has no further precondition. The handle is leaked below and never
        // unloaded, because unloading a Fortran runtime is not supported.
        let lib = unsafe { libloading::Library::new(path) }
            .map_err(|e| Error::Load(format!("{}: {e}", path.display())))?;
        let lib: &'static libloading::Library = Box::leak(Box::new(lib));
        let api = Api {
            get_version: symbol(lib, "dftd3_get_version")?,
            new_error: symbol(lib, "dftd3_new_error")?,
            check_error: symbol(lib, "dftd3_check_error")?,
            get_error: symbol(lib, "dftd3_get_error")?,
            delete_error: symbol(lib, "dftd3_delete_error")?,
            new_structure: symbol(lib, "dftd3_new_structure")?,
            delete_structure: symbol(lib, "dftd3_delete_structure")?,
            new_d3_model: symbol(lib, "dftd3_new_d3_model")?,
            delete_model: symbol(lib, "dftd3_delete_model")?,
            load_rational_damping: symbol(lib, "dftd3_load_rational_damping")?,
            load_zero_damping: symbol(lib, "dftd3_load_zero_damping")?,
            new_rational_damping: symbol(lib, "dftd3_new_rational_damping")?,
            new_zero_damping: symbol(lib, "dftd3_new_zero_damping")?,
            delete_param: symbol(lib, "dftd3_delete_param")?,
            get_dispersion: symbol(lib, "dftd3_get_dispersion")?,
        };
        let api: &'static Api = Box::leak(Box::new(api));
        let major = api.version_number() / 10000;
        if major != 1 {
            return Err(Error::Load(format!(
                "{}: library version {} is not 1.x (adapter written for the 1.x C API)",
                path.display(),
                api.version()
            )));
        }
        list.push((path.to_path_buf(), api));
        Ok(api)
    }

    fn version_number(&self) -> c_int {
        // SAFETY: no arguments; see the module preconditions.
        unsafe { (self.get_version)() }
    }

    /// "major.minor.patch".
    pub(crate) fn version(&self) -> String {
        let v = self.version_number();
        format!("{}.{}.{}", v / 10000, (v / 100) % 100, v % 100)
    }

    pub(crate) fn dispersion(
        &self,
        s: &Structure<'_>,
        p: &Parameters,
    ) -> Result<Evaluation, Error> {
        let _guard = CALL_LOCK
            .lock()
            .map_err(|_| Error::Library("call lock poisoned".into()))?;
        let n = s.numbers.len();
        let natoms = c_int::try_from(n).map_err(|_| Error::Structure("too many atoms".into()))?;
        let err = Owned::new(
            // SAFETY: no arguments; returns a fresh handle or null.
            unsafe { (self.new_error)() },
            self.delete_error,
        )
        .ok_or_else(|| Error::Library("dftd3_new_error returned null".into()))?;
        let positions: Vec<f64> = s.positions_bohr.iter().flatten().copied().collect();
        let lattice: Vec<f64> = s.lattice_rows_bohr.iter().flatten().copied().collect();
        let periodic = [true; 3];
        // SAFETY: `err` is live; `numbers` has `natoms` entries, `positions`
        // 3*natoms, `lattice` 9 and `periodic` 3 (validated by the caller and
        // built above); the library copies them.
        let mol = unsafe {
            (self.new_structure)(
                err.0,
                natoms,
                s.numbers.as_ptr(),
                positions.as_ptr(),
                lattice.as_ptr(),
                periodic.as_ptr(),
            )
        };
        let mol = self.checked(&err, mol, self.delete_structure, "dftd3_new_structure")?;
        // SAFETY: `err` and `mol` are live handles.
        let model = unsafe { (self.new_d3_model)(err.0, mol.0) };
        let model = self.checked(&err, model, self.delete_model, "dftd3_new_d3_model")?;
        let param = match p {
            Parameters::Library {
                damping,
                functional,
                three_body,
            } => {
                let name = CString::new(functional.as_str())
                    .map_err(|_| Error::Parameters("functional name contains NUL".into()))?;
                let load = match damping {
                    Damping::Rational => self.load_rational_damping,
                    Damping::Zero => self.load_zero_damping,
                };
                // SAFETY: `err` is live and `name` is a NUL-terminated string
                // that outlives the call; the library only reads it (the
                // header's `char*` is not const but the string is not written).
                unsafe { load(err.0, name.as_ptr().cast_mut(), *three_body) }
            }
            Parameters::Rational {
                s6,
                s8,
                s9,
                a1,
                a2,
                alpha,
            } => {
                // SAFETY: `err` is live; plain double arguments in header order
                // (s6, s8, s9, a1, a2, alp).
                unsafe { (self.new_rational_damping)(err.0, *s6, *s8, *s9, *a1, *a2, *alpha) }
            }
            Parameters::Zero {
                s6,
                s8,
                s9,
                rs6,
                rs8,
                alpha,
            } => {
                // SAFETY: `err` is live; plain double arguments in header order
                // (s6, s8, s9, rs6, rs8, alp).
                unsafe { (self.new_zero_damping)(err.0, *s6, *s8, *s9, *rs6, *rs8, *alpha) }
            }
        };
        let param = self.checked(&err, param, self.delete_param, "damping parameters")?;
        let mut energy = 0.0_f64;
        let mut gradient = vec![0.0_f64; 3 * n];
        let mut sigma = [0.0_f64; 9];
        // SAFETY: all handles are live; `energy` is one double, `gradient`
        // holds 3*natoms doubles and `sigma` 9, matching the structure the
        // model was created for.
        unsafe {
            (self.get_dispersion)(
                err.0,
                mol.0,
                model.0,
                param.0,
                &mut energy,
                gradient.as_mut_ptr(),
                sigma.as_mut_ptr(),
            )
        };
        self.status(&err, "dftd3_get_dispersion")?;
        drop((param, model, mol));
        if !energy.is_finite() || gradient.iter().chain(&sigma).any(|v| !v.is_finite()) {
            return Err(Error::Library("non-finite dispersion result".into()));
        }
        let gradient = gradient.as_chunks::<3>().0.to_vec();
        let virial = std::array::from_fn(|a| std::array::from_fn(|b| sigma[3 * a + b]));
        Ok((energy, gradient, virial))
    }

    fn status(&self, err: &Owned, what: &str) -> Result<(), Error> {
        // SAFETY: `err` is live.
        if unsafe { (self.check_error)(err.0) } == 0 {
            return Ok(());
        }
        let mut buffer = vec![0 as c_char; 512];
        let size: c_int = buffer.len() as c_int;
        // SAFETY: `buffer` holds `size` bytes; the library writes at most
        // `size` bytes including the terminating NUL.
        unsafe { (self.get_error)(err.0, buffer.as_mut_ptr(), &size) };
        let bytes: Vec<u8> = buffer
            .iter()
            .take_while(|c| **c != 0)
            .map(|c| *c as u8)
            .collect();
        Err(Error::Library(format!(
            "{what}: {}",
            String::from_utf8_lossy(&bytes).trim()
        )))
    }

    fn checked(&self, err: &Owned, h: Handle, delete: Delete, what: &str) -> Result<Owned, Error> {
        let owned = Owned::new(h, delete);
        self.status(err, what)?;
        owned.ok_or_else(|| Error::Library(format!("{what} returned null")))
    }
}

/// Owned library handle, deleted on drop with its own `dftd3_delete_*`.
struct Owned(Handle, Delete);

impl Owned {
    fn new(h: Handle, delete: Delete) -> Option<Self> {
        (!h.is_null()).then_some(Self(h, delete))
    }
}

impl Drop for Owned {
    fn drop(&mut self) {
        // SAFETY: the handle was created by the library and is deleted exactly
        // once, here, by the matching delete function (which takes a pointer
        // to the handle and nulls it).
        unsafe { (self.1)(&mut self.0) };
    }
}
