#[allow(dead_code)]
mod common;
use std::sync::atomic::{AtomicUsize, Ordering};
use std::{fs, path::PathBuf, process::Command};
static SERIAL: AtomicUsize = AtomicUsize::new(0);
struct Case(PathBuf);
impl Case {
    fn new() -> Self {
        let dir = std::env::temp_dir().join(format!(
            "dft-spectrum-{}-{}",
            std::process::id(),
            SERIAL.fetch_add(1, Ordering::SeqCst)
        ));
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("synthetic.xml"), common::fixture()).unwrap();
        let text = "[system]\nname='synthetic spectrum'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.0,4.0,0.0],[0.0,0.0,4.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nSi='synthetic.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[bands]\npoints_per_segment=2\npoints=[{label='Gamma',fractional=[0.0,0.0,0.0]},{label='A',fractional=[0.1,0.0,0.0]},{label='B',fractional=[0.0,0.1,0.0],break_before=true},{label='Gamma',fractional=[0.0,0.0,0.0]}]\n[dos]\ndos_kgrid=[1,1,1]\n";
        fs::write(dir.join("input.toml"), text).unwrap();
        Self(dir)
    }
    fn run(&self, command: &str) -> std::process::Output {
        Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&self.0)
            .args([command, "input.toml", "--output-dir", "out"])
            .output()
            .unwrap()
    }
}
impl Drop for Case {
    fn drop(&mut self) {
        fs::remove_dir_all(&self.0).unwrap();
    }
}
#[test]
fn portable_formats_independently_read() {
    let c = Case::new();
    for command in ["bands", "dos"] {
        let o = c.run(command);
        assert!(
            o.status.success(),
            "{}\n{}",
            String::from_utf8_lossy(&o.stdout),
            String::from_utf8_lossy(&o.stderr)
        );
    }
    let o = Command::new("python3")
        .args([
            "-I",
            concat!(env!("CARGO_MANIFEST_DIR"), "/tests/read_spectra.py"),
        ])
        .arg(c.0.join("out"))
        .output()
        .unwrap();
    assert!(o.status.success(), "{}", String::from_utf8_lossy(&o.stderr));
}
#[test]
fn refuses_wrong_preset_and_protected_output() {
    let c = Case::new();
    let input = c.0.join("input.toml");
    let text = fs::read_to_string(&input).unwrap();
    fs::write(
        &input,
        text.split("[bands]").next().unwrap().to_owned() + "[bands]\npreset='fcc'\n",
    )
    .unwrap();
    let o = c.run("bands");
    assert!(!o.status.success());
    assert!(String::from_utf8_lossy(&o.stderr).contains("primitive"));
    fs::write(&input, text).unwrap();
    fs::create_dir_all(c.0.join("out")).unwrap();
    fs::hard_link(&input, c.0.join("out/bands.csv")).unwrap();
    let before = fs::read(&input).unwrap();
    let o = c.run("bands");
    assert!(!o.status.success());
    assert_eq!(before, fs::read(&input).unwrap());
}
