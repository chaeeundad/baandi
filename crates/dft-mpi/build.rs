use std::{env, process::Command};

fn main() {
    println!("cargo:rerun-if-env-changed=MPICC");
    println!("cargo:rerun-if-env-changed=MPI_ROOT");
    println!("cargo:rustc-check-cfg=cfg(mpi_linked)");
    if env::var_os("CARGO_FEATURE_MPI").is_none() {
        return;
    }
    // Missing MPI is permitted for cargo check (type checking only); a real build
    // then fails at link, never silently falling back to the stand-in.
    let mpicc = env::var("MPICC").unwrap_or_else(|_| "mpicc".into());
    let Ok(out) = Command::new(&mpicc).arg("-show").output() else {
        println!(
            "cargo:warning=`{mpicc} -show` unavailable: type checking only; MPI executables cannot link"
        );
        return;
    };
    if !out.status.success() {
        println!("cargo:warning=`{mpicc} -show` failed: type checking only");
        return;
    }
    let text = String::from_utf8_lossy(&out.stdout);
    // The FFI declarations follow the Open MPI ABI (pointer handles, ompi_mpi_*
    // predefined objects). Refuse anything else instead of mis-linking.
    let open_mpi = Command::new(&mpicc)
        .arg("--showme:version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    assert!(
        open_mpi,
        "dft-mpi supports the Open MPI ABI only (`{mpicc} --showme:version` failed)"
    );
    for token in text.split_whitespace() {
        if let Some(path) = token.strip_prefix("-L") {
            println!("cargo:rustc-link-search=native={path}");
        } else if let Some(lib) = token.strip_prefix("-l") {
            println!("cargo:rustc-link-lib=dylib={lib}");
        } else if let Some(rest) = token.strip_prefix("-Wl,") {
            // -Wl,-rpath,<dir> and friends are passed through.
            println!("cargo:rustc-link-arg=-Wl,{rest}");
        }
    }
    println!("cargo:rustc-cfg=mpi_linked");
}
