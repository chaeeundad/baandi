//! R2 `baandi dielectric` adapter gates on a synthetic UPF (never reference
//! evidence): the outputs are written, the CSV and gnuplot views carry the
//! JSON record's numbers (independent parse of the documented columns), and
//! unsupported keys reject.
#[allow(dead_code)]
mod common;
use serde_json::Value;
use std::{fs, path::Path, path::PathBuf, process::Command};

fn atomic(z: f64) -> String {
    let mut pp = common::gaussian_with(z, 0.8, true, 401)
        .replace("SLA PZ", "SLA PW")
        .replace(
            "<PP_DIJ>2 0.4 0 0.4 4 0 0 0 6</PP_DIJ>",
            "<PP_DIJ>-1 0.1 0 0.1 -0.5 0 0 0 0.8</PP_DIJ>",
        );
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
        let dir = std::env::temp_dir().join(format!("r2-dielectric-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("ion.upf"), atomic(2.)).unwrap();
        Self(dir)
    }
    fn input(&self, label: &str, response: &str) -> PathBuf {
        let text = format!(
            "[system]\nname='R2 dielectric'\nlattice_unit='bohr'\nlattice_vectors=[[6.0,0.0,0.0],[0.7,5.5,0.0],[0.3,0.4,6.5]]\n[[system.atoms]]\nspecies='Xx'\nfractional=[0.13,0.21,0.37]\n[pseudopotentials]\nXx='ion.upf'\n[basis]\necutwfc_ry=8.0\nkgrid=[2,2,2]\n[xc]\nfunctional='lda-pw'\n[scf]\nmax_iterations=150\n[response]\n{response}\n"
        );
        let path = self.0.join(format!("{label}.toml"));
        fs::write(&path, text).unwrap();
        path
    }
    fn run(&self, input: &Path, out: &str) -> std::process::Output {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        for (key, _) in std::env::vars().filter(|(k, _)| k.starts_with("BAANDI_")) {
            cmd.env_remove(key);
        }
        cmd.arg("dielectric")
            .arg(input)
            .arg("--output-dir")
            .arg(self.0.join(out));
        cmd.output().unwrap()
    }
}

const BASE: &str = "ecut_chi_ha=3.0\nbands=24\nk_mesh=[2,2,2]\n[response.frequency]\npoints=4\n";

#[test]
fn dielectric_outputs_and_guards() {
    let case = Case::new();
    let o = case.run(&case.input("ok", BASE), "out");
    assert!(
        o.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
    let dir = case.0.join("out");
    for name in [
        "epsilon.json",
        "epsilon_M.csv",
        "epsilon_M.dat",
        "chi0_head.dat",
        "results.json",
    ] {
        assert!(dir.join(name).is_file(), "{name}");
    }
    let meta: Value = serde_json::from_slice(&fs::read(dir.join("epsilon.json")).unwrap()).unwrap();
    assert_eq!(meta["kind"], "baandi_dielectric_r2");
    assert_eq!(meta["settings"]["bands"], 24);
    let freqs: Vec<f64> = meta["frequencies_hartree"]
        .as_array()
        .unwrap()
        .iter()
        .map(|v| v.as_f64().unwrap())
        .collect();
    assert_eq!(freqs.len(), 5);
    assert_eq!(freqs[0], 0.);
    let lf = &meta["epsilon_with_local_fields"];
    let nlf = &meta["epsilon_without_local_fields"];
    // CSV: header names, then one row per frequency.
    let csv = fs::read_to_string(dir.join("epsilon_M.csv")).unwrap();
    let mut lines = csv.lines().filter(|l| !l.starts_with('#'));
    let header: Vec<&str> = lines.next().unwrap().split(',').collect();
    let col = |name: &str| header.iter().position(|h| *h == name).unwrap();
    let rows: Vec<Vec<f64>> = lines
        .map(|l| l.split(',').map(|x| x.parse().unwrap()).collect())
        .collect();
    assert_eq!(rows.len(), freqs.len());
    for (i, row) in rows.iter().enumerate() {
        assert_eq!(row[col("omega_ha")], freqs[i]);
        for (name, a, b) in [("lf_xx", 0, 0), ("lf_yz", 1, 2), ("nlf_zz", 2, 2)] {
            let source = if name.starts_with("lf") { lf } else { nlf };
            assert_eq!(row[col(name)], source[i][a][b].as_f64().unwrap(), "{name}");
        }
        let avg = (row[col("lf_xx")] + row[col("lf_yy")] + row[col("lf_zz")]) / 3.;
        assert!((avg - row[col("lf_avg")]).abs() <= 1e-14 * avg);
        assert!(row[col("lf_xx")] > 1. && row[col("nlf_xx")] > 1.);
    }
    assert_eq!(
        rows[0][col("lf_avg")],
        meta["epsilon_infinity"]["average_with_local_fields"]
            .as_f64()
            .unwrap()
    );
    let dat = fs::read_to_string(dir.join("epsilon_M.dat")).unwrap();
    assert_eq!(
        dat.lines().filter(|l| !l.starts_with('#')).count(),
        freqs.len()
    );
    // Unsupported settings reject before any work.
    for (label, text) in [
        (
            "head",
            BASE.replace("bands=24", "bands=24\nhead='finite-q'"),
        ),
        ("grid", format!("{BASE}grid='minimax'\n")),
        (
            "solver",
            BASE.replace("bands=24", "bands=24\nband_solver='davidson'"),
        ),
        ("unknown", BASE.replace("bands=24", "bands=24\nq_points=3")),
    ] {
        let o = case.run(&case.input(label, &text), label);
        assert!(!o.status.success(), "{label} must reject");
    }
}
