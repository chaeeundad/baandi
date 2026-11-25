//! R4b adapter gates. External material runs are opt-in; no UPF is committed.
#[allow(dead_code)]
mod common;
use serde_json::Value;
use std::{
    fs,
    path::PathBuf,
    process::Command,
    sync::atomic::{AtomicUsize, Ordering},
};
static SERIAL: AtomicUsize = AtomicUsize::new(0);
struct Case(PathBuf);
impl Case {
    fn new(upf: &str, material: bool) -> Self {
        let dir = std::env::temp_dir().join(format!(
            "r4b-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::SeqCst)
        ));
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("ion.upf"), upf).unwrap();
        let cell = if material {
            "[[0.0,5.13,5.13],[5.13,0.0,5.13],[5.13,5.13,0.0]]"
        } else {
            "[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]"
        };
        let extra = if material {
            "[[system.atoms]]\nspecies='Si'\nfractional=[0.25,0.25,0.25]\n"
        } else {
            ""
        };
        let text = format!(
            "[system]\nname='R4b spectra'\nlattice_unit='bohr'\nlattice_vectors={cell}\n[[system.atoms]]\nspecies='Si'\nfractional=[0.0,0.0,0.0]\n{extra}[pseudopotentials]\nSi='ion.upf'\n[basis]\necutwfc_ry={}\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='{}'\n[scf]\nmax_iterations=150\n[bands]\npoints_per_segment=2\npoints=[{{label='Gamma',fractional=[0.0,0.0,0.0]}},{{label='A',fractional=[0.1,0.0,0.0]}}]\n[dos]\ndos_kgrid=[2,1,1]\nsigma_ev=0.2\nspacing_ev=0.05\n",
            if material { 6.0 } else { 4.0 },
            if material { "lda-pz" } else { "lda-pw" }
        );
        let text = if material {
            text.replace(
                "{label='A',fractional=[0.1,0.0,0.0]}",
                "{label='X',fractional=[0.5,0.0,0.5]}",
            )
        } else {
            text
        };
        // Legacy NSCF still stops buffer bands loosely. A nearly full basis
        // makes its residual guard certify all exported bands without weakening
        // the 1e-8 Ha residual or 1e-7 Ha cross-driver eigenvalue gates.
        let config = dft_cli::config::Config::parse(&text).unwrap();
        let cell = config.cell().unwrap();
        let mut count = usize::MAX;
        for x in -2..=2 {
            for y in -2..=2 {
                for z in -2..=2 {
                    let point = [x as f64 / 4., y as f64 / 4., z as f64 / 4.];
                    count = count.min(
                        dft_basis::PlaneWaveBasis::at_k(
                            &cell,
                            dft_units::Hartree(if material { 3. } else { 2. }),
                            point,
                        )
                        .unwrap()
                        .vectors()
                        .len(),
                    );
                }
            }
        }
        let text = text
            .replace("[bands]", &format!("[bands]\nnbands={count}"))
            .replace("[dos]", &format!("[dos]\nnbands={count}"));
        fs::write(dir.join("base.toml"), text).unwrap();
        Self(dir)
    }
    fn run(
        &self,
        command: &str,
        label: &str,
        controls: &str,
        spin: &str,
        strict: bool,
    ) -> (std::process::Output, PathBuf) {
        let text = fs::read_to_string(self.0.join("base.toml"))
            .unwrap()
            .replace("[scf]", &format!("[scf]\n{controls}"))
            + spin;
        let input = self.0.join(format!("{label}.toml"));
        let out = self.0.join(format!("{command}-{label}"));
        fs::write(&input, text).unwrap();
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        for (key, _) in std::env::vars().filter(|(k, _)| k.starts_with("BAANDI_")) {
            cmd.env_remove(key);
        }
        if strict {
            cmd.env("BAANDI_REFERENCE_STRICT", "on");
        }
        let result = cmd
            .arg(command)
            .arg(input)
            .arg("--output-dir")
            .arg(&out)
            .output()
            .unwrap();
        (result, out)
    }
    fn success(
        &self,
        command: &str,
        label: &str,
        controls: &str,
        spin: &str,
        strict: bool,
    ) -> (Value, PathBuf) {
        let (output, out) = self.run(command, label, controls, spin, strict);
        assert!(
            output.status.success(),
            "{command}/{label}: {}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let json = serde_json::from_slice(&fs::read(out.join(format!("{command}.json"))).unwrap())
            .unwrap();
        (json, out)
    }
}
impl Drop for Case {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}
fn atomic(z: f64) -> String {
    let mut pp = common::gaussian_with(z, 1.5, true, 101).replace("SLA PZ", "SLA PW");
    let p = dft_upf::parse_upf(&pp).unwrap();
    let array = |f: &dyn Fn(f64) -> f64| {
        p.mesh()
            .r()
            .iter()
            .map(|&r| format!("{:.17e}", f(r)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let rho = array(&|r| {
        4. * z
            * std::f64::consts::PI
            * r
            * r
            * (1. / (2.25 * std::f64::consts::PI)).powf(1.5)
            * (-r * r / 2.25).exp()
    });
    let start = pp.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
    let end = pp.find("</PP_RHOATOM>").unwrap();
    pp.replace_range(start..end, &rho);
    pp.replace("number_of_wfc=\"0\"", "number_of_wfc=\"2\"").replace("</UPF>", &format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1><PP_CHI.2 l=\"1\" label=\"p\" occupation=\"0\">{}</PP_CHI.2></PP_PSWFC></UPF>", array(&|r| r*(-0.6*r*r).exp()), array(&|r| r*r*(-0.5*r*r).exp())))
}
fn route(v: &Value, expected: &str) {
    assert_eq!(v["scf"]["driver"], expected);
    assert_eq!(v["nscf"]["driver"], expected);
    for channel in v["band_data"]["residuals_ha"].as_array().unwrap() {
        for point in channel.as_array().unwrap() {
            for value in point.as_array().unwrap() {
                assert!(value.as_f64().unwrap() <= v["eigen_tolerance_hartree"].as_f64().unwrap());
            }
        }
    }
}
fn read_spin(out: &std::path::Path, command: &str) {
    let output = Command::new("python3")
        .args([
            "-I",
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/read_spin_spectra.py"),
        ])
        .arg(out)
        .arg(command)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
}
fn compare(a: &Value, b: &Value, channel: usize) -> f64 {
    let a = a["band_data"]["eigenvalues_ha"][channel]
        .as_array()
        .unwrap();
    let b = b["band_data"]["eigenvalues_ha"][0].as_array().unwrap();
    assert_eq!(a.len(), b.len());
    let mut max: f64 = 0.;
    for (a, b) in a.iter().zip(b) {
        let (a, b) = (a.as_array().unwrap(), b.as_array().unwrap());
        assert_eq!(a.len(), b.len());
        for (a, b) in a.iter().zip(b) {
            max = max.max((a.as_f64().unwrap() - b.as_f64().unwrap()).abs());
        }
    }
    assert!(
        max <= 1e-7,
        "max eigenvalue difference {max:e} Ha exceeds 1e-7 Ha"
    );
    max
}
#[test]
fn reference_formats_and_unconverged_scf() {
    let case = Case::new(&atomic(2.), false);
    let (v, bands) = case.success("bands", "default", "", "", false);
    route(&v, "reference");
    let (v, dos) = case.success("dos", "default", "", "", false);
    route(&v, "reference");
    for name in [
        "dos.csv", "dos.dat", "dos.agr", "dos.json", "DOSCAR", "dos.svg",
    ] {
        fs::copy(dos.join(name), bands.join(name)).unwrap();
    }
    let output = Command::new("python3")
        .args([
            "-I",
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/read_spectra.py"),
        ])
        .arg(bands)
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    let base = fs::read_to_string(case.0.join("base.toml"))
        .unwrap()
        .replace("max_iterations=150", "max_iterations=1");
    fs::write(case.0.join("base.toml"), base).unwrap();
    let (output, out) = case.run("bands", "exhausted", "", "", false);
    assert_eq!(output.status.code(), Some(3));
    assert!(!out.join("bands.csv").exists());
}
// Cross-driver eigenvalue agreement is checked at a tightened SCF threshold
// (5e-13 Ha): at the default standard class the SCF densities of the two
// drivers differ by design and propagate to ~3e-6 Ha in eigenvalues.
#[test]
#[ignore = "tight-input diagnostic; DFT_SI_UPF=... cargo test --release -p dft-cli --test reference_bands_dos silicon_tight_band_path -- --ignored --nocapture"]
fn silicon_tight_band_path() {
    silicon_comparison("standard_threshold_ha=5e-13\n", &["bands"]);
}
fn silicon_comparison(controls: &str, commands: &[&str]) {
    let pp =
        fs::read_to_string(std::env::var("DFT_SI_UPF").expect("external Si UPF required")).unwrap();
    let case = Case::new(&pp, true);
    for &command in commands {
        let (reference, _) = case.success(command, "reference", controls, "", false);
        let (legacy, _) = case.success(
            command,
            "legacy",
            &format!("driver='legacy'\n{controls}"),
            "",
            false,
        );
        route(&reference, "reference");
        route(&legacy, "legacy");
        println!(
            "Si2 {command}: max eigenvalue difference {:.12e} Ha",
            compare(&reference, &legacy, 0)
        );
    }
}
#[test]
fn spin_channels_and_legacy_pin() {
    let case = Case::new(&atomic(2.), false);
    let input = case.0.join("base.toml");
    fs::write(
        &input,
        fs::read_to_string(&input)
            .unwrap()
            .replace("[scf]", "[scf]\nstandard_threshold_ha=5e-13"),
    )
    .unwrap();
    for command in ["bands", "dos"] {
        let (scalar, _) = case.success(command, "legacy", "driver='legacy'\n", "", false);
        let (spin, out) = case.success(
            command,
            "spin",
            "",
            "[spin]\nmode='collinear'\nn_up=1\nn_down=1\n",
            false,
        );
        route(&scalar, "legacy");
        route(&spin, "reference");
        read_spin(&out, command);
        assert_eq!(spin["band_data"]["shape"][0], 2);
        for channel in 0..2 {
            println!(
                "balanced spin {command} channel {channel}: max difference {:.12e} Ha",
                compare(&spin, &scalar, channel)
            );
        }
        let numbers = |s: &str| {
            s.split_whitespace()
                .map(|v| v.parse::<f64>().unwrap())
                .collect::<Vec<_>>()
        };
        if command == "bands" {
            let text = fs::read_to_string(out.join("EIGENVAL")).unwrap();
            let lines: Vec<_> = text.lines().collect();
            assert_eq!(numbers(lines[0])[3], 2.);
            assert_eq!(numbers(lines[8]).len(), 5);
            assert_eq!(&numbers(lines[8])[3..], &[1., 1.]);
        } else {
            let text = fs::read_to_string(out.join("DOSCAR")).unwrap();
            let last = numbers(text.lines().last().unwrap());
            assert_eq!(last.len(), 5);
            assert!((last[3] - spin["nbands"].as_f64().unwrap()).abs() < 1e-10);
            assert!((last[4] - spin["nbands"].as_f64().unwrap()).abs() < 1e-10);
        }
    }
    let (strict, _) = case.success("bands", "strict", "driver='reference'\n", "", true);
    route(&strict, "strict_oracle");
}

#[test]
fn polarized_density_is_not_collapsed() {
    let case = Case::new(&atomic(1.), false);
    // The single-electron fixture is insulating at Gamma, not over the full BZ.
    let path = case.0.join("base.toml");
    fs::write(
        &path,
        fs::read_to_string(&path)
            .unwrap()
            .replace("dos_kgrid=[2,1,1]", "dos_kgrid=[1,1,1]"),
    )
    .unwrap();
    for command in ["bands", "dos"] {
        let (spin, out) = case.success(
            command,
            "polarized",
            "",
            "[spin]\nmode='collinear'\nn_up=1\nn_down=0\n",
            false,
        );
        route(&spin, "reference");
        read_spin(&out, command);
        assert_eq!(spin["spin"]["occupied_bands"], serde_json::json!([1, 0]));
        let up = spin["band_data"]["eigenvalues_ha"][0][0][0]
            .as_f64()
            .unwrap();
        let down = spin["band_data"]["eigenvalues_ha"][1][0][0]
            .as_f64()
            .unwrap();
        assert!(
            (up - down).abs() > 1e-4,
            "distinct polarized potentials must reach NSCF"
        );
        println!(
            "polarized {command}: first-band splitting {:.12e} Ha",
            down - up
        );
    }
    let (failure, _) = case.run(
        "bands",
        "unsupported-legacy-spin",
        "driver='legacy'\n",
        "[spin]\nmode='collinear'\nn_up=1\nn_down=0\n",
        false,
    );
    assert_eq!(failure.status.code(), Some(1));
    assert!(
        String::from_utf8_lossy(&failure.stderr)
            .contains("legacy NSCF supports only total density")
    );
}

#[test]
#[ignore = "tight-input per-k legacy oracle; DFT_SI_UPF=... cargo test --release -p dft-cli --test reference_bands_dos silicon_tight_dos_per_k_oracle -- --ignored --nocapture"]
fn silicon_tight_dos_per_k_oracle() {
    let pp =
        fs::read_to_string(std::env::var("DFT_SI_UPF").expect("external Si UPF required")).unwrap();
    let case = Case::new(&pp, true);
    let base_path = case.0.join("base.toml");
    let base = fs::read_to_string(&base_path)
        .unwrap()
        .replace("dos_kgrid=[2,1,1]", "dos_kgrid=[2,2,2]");
    fs::write(&base_path, &base).unwrap();
    let controls = "standard_threshold_ha=5e-13\n";
    let (reference, _) = case.success("dos", "reference", controls, "", false);
    route(&reference, "reference");
    let config = dft_cli::config::Config::parse(&base).unwrap();
    let cell = config.cell().unwrap();
    let mut max = 0_f64;
    for (k, point) in reference["kpoints"].as_array().unwrap().iter().enumerate() {
        let p: [f64; 3] = std::array::from_fn(|i| point["fractional"][i].as_f64().unwrap());
        let mut equivalent = p;
        equivalent[0] += 1.;
        let count = dft_basis::PlaneWaveBasis::at_k(&cell, dft_units::Hartree(3.), p)
            .unwrap()
            .vectors()
            .len();
        // R4a's full-basis legacy oracle, exported through the CLI. The two
        // path endpoints differ by a reciprocal vector and have the same H.
        // Only the NSCF request changes; the SCF input/density is identical.
        let text = format!(
            "{}[bands]\nnbands={count}\npoints_per_segment=1\npoints=[{{label='k',fractional={p:?}}},{{label='k+G',fractional={equivalent:?}}}]\n",
            base.split("[bands]").next().unwrap()
        );
        fs::write(&base_path, text).unwrap();
        let (legacy, _) = case.success(
            "bands",
            &format!("legacy-{k}"),
            &format!("driver='legacy'\n{controls}"),
            "",
            false,
        );
        route(&legacy, "legacy");
        for (band, value) in reference["band_data"]["eigenvalues_ha"][0][k]
            .as_array()
            .unwrap()
            .iter()
            .enumerate()
        {
            let delta = (value.as_f64().unwrap()
                - legacy["band_data"]["eigenvalues_ha"][0][0][band]
                    .as_f64()
                    .unwrap())
            .abs();
            max = max.max(delta);
        }
    }
    println!(
        "Si2 DOS k222 per-k full-basis legacy oracle: max eigenvalue difference {max:.12e} Ha"
    );
    assert!(max <= 1e-7, "max difference {max:e} Ha");
}
