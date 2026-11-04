use std::{env, path::Path};

fn main() {
    println!("cargo:rerun-if-env-changed=MKLROOT");
    println!("cargo:rustc-check-cfg=cfg(mkl_linked)");
    if env::var_os("CARGO_FEATURE_MKL").is_none() {
        return;
    }
    // Missing MKL is permitted for cargo check. A real build then fails at link,
    // never silently substituting the Rust backend.
    let Some(root) = env::var_os("MKLROOT") else {
        println!("cargo:warning=MKLROOT unset: type checking only; MKL executables cannot link");
        return;
    };
    assert_eq!(env::var("CARGO_CFG_TARGET_ARCH").unwrap(), "x86_64");
    assert_eq!(env::var("CARGO_CFG_TARGET_OS").unwrap(), "linux");
    let root = Path::new(&root);
    let lib = [root.join("lib/intel64"), root.join("lib")]
        .into_iter()
        .find(|p| p.join("libmkl_core.so").exists())
        .expect("MKLROOT must contain lib[/intel64]/libmkl_core.so");
    println!("cargo:rustc-link-search=native={}", lib.display());
    // The MKL layers resolve each other's symbols at load time, so Rust objects
    // never reference mkl_sequential/mkl_core directly. The default linker
    // --as-needed would drop them and the loader fails with "undefined symbol:
    // mkl_blas_dgemm". Keep them with the `-as-needed` link modifier.
    for name in ["mkl_intel_lp64", "mkl_sequential", "mkl_core"] {
        println!("cargo:rustc-link-lib=dylib:-as-needed={name}");
    }
    for name in ["pthread", "m", "dl"] {
        println!("cargo:rustc-link-lib=dylib={name}");
    }
    println!("cargo:rustc-cfg=mkl_linked");
}
