//! R1 `baandi nscf` adapter gates on a synthetic UPF (never reference
//! evidence): dense and Chebyshev agree, outputs are written, and per-k
//! `nscf_bands_v1` units resume bit-identically, while corrupt or
//! mismatching units reject.
#[allow(dead_code)]
mod common;
use serde_json::Value;
use std::{fs, path::Path, path::PathBuf, process::Command};

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

struct Case(PathBuf);
impl Drop for Case {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
impl Case {
    fn new() -> Self {
        let dir = std::env::temp_dir().join(format!("r1-nscf-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("ion.upf"), atomic(2.)).unwrap();
        Self(dir)
    }
    fn input(&self, label: &str, solver: &str) -> PathBuf {
        let text = format!(
            "[system]\nname='R1 nscf'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nSi='ion.upf'\n[basis]\necutwfc_ry=4.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pw'\n[scf]\nmax_iterations=150\n[nscf]\nbands=12\nband_solver='{solver}'\nk_mesh=[2,1,1]\n"
        );
        let path = self.0.join(format!("{label}.toml"));
        fs::write(&path, text).unwrap();
        path
    }
    fn run(&self, input: &Path, out: &str, checkpoint: Option<&str>) -> std::process::Output {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        for (key, _) in std::env::vars().filter(|(k, _)| k.starts_with("BAANDI_")) {
            cmd.env_remove(key);
        }
        cmd.arg("nscf")
            .arg(input)
            .arg("--output-dir")
            .arg(self.0.join(out));
        if let Some(c) = checkpoint {
            cmd.arg("--checkpoint-dir").arg(self.0.join(c));
        }
        cmd.output().unwrap()
    }
    fn ok(&self, input: &Path, out: &str, checkpoint: Option<&str>) -> Value {
        let o = self.run(input, out, checkpoint);
        assert!(
            o.status.success(),
            "{}\n{}",
            String::from_utf8_lossy(&o.stdout),
            String::from_utf8_lossy(&o.stderr)
        );
        for name in ["nscf_eigenvalues.csv", "EIGENVAL", "results.json"] {
            assert!(self.0.join(out).join(name).is_file(), "{name}");
        }
        serde_json::from_slice(&fs::read(self.0.join(out).join("nscf.json")).unwrap()).unwrap()
    }
}

fn eigenvalues(v: &Value) -> Vec<f64> {
    v["band_data"]["eigenvalues_ha"]
        .as_array()
        .unwrap()
        .iter()
        .flat_map(|c| c.as_array().unwrap().iter())
        .flat_map(|k| k.as_array().unwrap().iter())
        .map(|e| e.as_f64().unwrap())
        .collect()
}

#[test]
fn nscf_dense_chebyshev_and_checkpoint_resume() {
    let case = Case::new();
    let dense = case.input("dense", "dense");
    let first = case.ok(&dense, "dense-1", Some("cp"));
    assert_eq!(first["band_solver"], "dense");
    assert!(first["max_residual_hartree"].as_f64().unwrap() <= 1e-6);
    assert_eq!(
        first["checkpoint"]["written_units"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert!(case.0.join("cp").join("nscf-manifest.json").is_file());
    // Resume: every unit loads, is re-certified, and the bands are unchanged.
    let second = case.ok(&dense, "dense-2", Some("cp"));
    assert_eq!(
        second["checkpoint"]["resumed_units"]
            .as_array()
            .unwrap()
            .len(),
        2
    );
    assert_eq!(
        second["checkpoint"]["written_units"]
            .as_array()
            .unwrap()
            .len(),
        0
    );
    let (a, b) = (eigenvalues(&first), eigenvalues(&second));
    assert_eq!(
        a.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
        b.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
    );
    // Chebyshev against dense.
    let cheb = case.ok(&case.input("cheb", "chebyshev"), "cheb", None);
    let error = eigenvalues(&cheb)
        .iter()
        .zip(&a)
        .map(|(x, y)| (x - y).abs())
        .fold(0., f64::max);
    assert!(error <= 1e-8, "Chebyshev vs dense {error:e} Ha");
    // A different solver is a different fingerprint: the stored units reject.
    let mismatch = case.run(&case.input("cheb2", "chebyshev"), "cheb2", Some("cp"));
    assert!(!mismatch.status.success());
    assert!(String::from_utf8_lossy(&mismatch.stderr).contains("fingerprint mismatch"));
    // A corrupted unit fails its checksum and rejects.
    let unit = case.0.join("cp").join("nscf-c0-k000001.cp");
    let mut bytes = fs::read(&unit).unwrap();
    let at = bytes.len() / 2;
    bytes[at] ^= 1;
    fs::remove_file(&unit).unwrap();
    fs::write(&unit, bytes).unwrap();
    let corrupt = case.run(&dense, "dense-3", Some("cp"));
    assert!(!corrupt.status.success());
    assert!(String::from_utf8_lossy(&corrupt.stderr).contains("checksum"));
}
