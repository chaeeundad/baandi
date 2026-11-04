// Equivalent to pyo3_build_config::add_extension_module_link_args on macOS,
// without a build dependency that would probe Python during workspace gates.
// Only the extension cdylib resolves Python symbols from its loading process;
// Rust test executables must never embed an interpreter.
fn main() {
    if std::env::var_os("CARGO_FEATURE_PYTHON_EXTENSION").is_some()
        && std::env::var("CARGO_CFG_TARGET_OS").as_deref() == Ok("macos")
    {
        println!("cargo:rustc-cdylib-link-arg=-undefined");
        println!("cargo:rustc-cdylib-link-arg=dynamic_lookup");
    }
}
