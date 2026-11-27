//! R3 `baandi rpa` adapter gates on a synthetic UPF (never reference
//! evidence): outputs written, the CSV views carry the JSON record's numbers
//! (independent parse of the documented columns), the q-weighted sum and the
//! frequency integral reproduce E_c, the per-q checkpoint resumes
//! bit-identically, `[task] type = "rpa"` selects the job, bad keys reject.
//! Fast gate: `cargo test -p dft-cli --test rpa_r3`
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
        let dir = std::env::temp_dir().join(format!("r3-rpa-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("ion.upf"), atomic(2.)).unwrap();
        Self(dir)
    }
    fn input(&self, label: &str, tables: &str) -> PathBuf {
        let text = format!(
            "[system]\nname='R3 rpa'\nlattice_unit='bohr'\nlattice_vectors=[[6.0,0.0,0.0],[0.7,5.5,0.0],[0.3,0.4,6.5]]\n[[system.atoms]]\nspecies='Xx'\nfractional=[0.13,0.21,0.37]\n[pseudopotentials]\nXx='ion.upf'\n[basis]\necutwfc_ry=8.0\nkgrid=[2,2,2]\n[xc]\nfunctional='lda-pw'\n[scf]\nmax_iterations=150\n{tables}\n"
        );
        let path = self.0.join(format!("{label}.toml"));
        fs::write(&path, text).unwrap();
        path
    }
    fn command(&self) -> Command {
        let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
        for (key, _) in std::env::vars().filter(|(k, _)| k.starts_with("BAANDI_")) {
            cmd.env_remove(key);
        }
        cmd.current_dir(&self.0);
        cmd
    }
    fn run(&self, input: &Path, out: &str, checkpoint: Option<&str>) -> std::process::Output {
        let mut cmd = self.command();
        cmd.arg("rpa")
            .arg(input)
            .arg("--output-dir")
            .arg(self.0.join(out));
        if let Some(c) = checkpoint {
            cmd.arg("--checkpoint-dir").arg(self.0.join(c));
        }
        cmd.output().unwrap()
    }
}

const BASE: &str = "[response]\necut_chi_ha=2.5\nk_mesh=[2,2,2]\n[response.frequency]\npoints=8\n[rpa]\necut_chi_ha_series=[1.5,2.0,2.5]\nq0_directions=14\nq_symmetry='time-reversal'\n";

fn ok(o: &std::process::Output) {
    assert!(
        o.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
}
fn floats(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .map(|x| x.as_f64().unwrap())
        .collect()
}
fn table(path: &Path) -> (Vec<String>, Vec<Vec<String>>) {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines().filter(|l| !l.starts_with('#'));
    let header = lines
        .next()
        .unwrap()
        .split(',')
        .map(str::to_owned)
        .collect();
    let rows = lines
        .map(|l| l.split(',').map(str::to_owned).collect())
        .collect();
    (header, rows)
}

#[test]
fn rpa_outputs_resume_and_guards() {
    let case = Case::new();
    let input = case.input("ok", BASE);
    let o = case.run(&input, "out", Some("ck"));
    ok(&o);
    let dir = case.0.join("out");
    for name in [
        "rpa.json",
        "rpa_q.csv",
        "rpa_extrapolation.csv",
        "rpa_frequency.csv",
        "results.json",
    ] {
        assert!(dir.join(name).is_file(), "{name}");
    }
    let meta: Value = serde_json::from_slice(&fs::read(dir.join("rpa.json")).unwrap()).unwrap();
    assert_eq!(meta["kind"], "baandi_rpa_r3");
    let energies = floats(&meta["e_c_hartree_per_cell"]);
    assert_eq!(energies.len(), 3);
    assert!(energies.iter().all(|e| *e < 0.));
    assert!(energies[2] < energies[0]);
    assert!(meta["extrapolation"]["e_c_inf_hartree_per_cell"].is_f64());
    assert_eq!(meta["settings"]["q0_directions"], 14);
    // rpa_q.csv: weighted sum over q reproduces E_c per cutoff.
    let (header, rows) = table(&dir.join("rpa_q.csv"));
    let col = |n: &str| header.iter().position(|h| h == n).unwrap();
    let cutoffs = floats(&meta["settings"]["ecut_chi_series_hartree"]);
    for (c, e) in cutoffs.iter().zip(&energies) {
        let sum: f64 = rows
            .iter()
            .filter(|r| r[col("ecut_chi_ha")].parse::<f64>().unwrap() == *c)
            .map(|r| r[col("weighted_ha")].parse::<f64>().unwrap())
            .sum();
        assert!((sum - e).abs() <= 1e-13, "{sum} vs {e}");
    }
    // rpa_frequency.csv: (1/2pi) sum weight * term = E_c.
    let (header, rows) = table(&dir.join("rpa_frequency.csv"));
    let col = |n: &str| header.iter().position(|h| h == n).unwrap();
    let integral: f64 = rows
        .iter()
        .map(|r| {
            r[col("weight_ha")].parse::<f64>().unwrap()
                * r[col("term_ecut_2.5")].parse::<f64>().unwrap()
        })
        .sum::<f64>()
        / (2. * std::f64::consts::PI);
    assert!((integral - energies[2]).abs() <= 1e-12 * energies[2].abs());
    // rpa_extrapolation.csv matches the JSON.
    let (header, rows) = table(&dir.join("rpa_extrapolation.csv"));
    let col = |n: &str| header.iter().position(|h| h == n).unwrap();
    for (r, e) in rows.iter().zip(&energies) {
        assert_eq!(r[col("e_c_ha_per_cell")].parse::<f64>().unwrap(), *e);
    }
    // Resume from the per-q records: identical numbers, all q resumed.
    let records = fs::read_dir(case.0.join("ck")).unwrap().count();
    assert_eq!(
        records,
        meta["settings"]["irreducible_q"].as_u64().unwrap() as usize
    );
    let o = case.run(&input, "out2", Some("ck"));
    ok(&o);
    let again: Value =
        serde_json::from_slice(&fs::read(case.0.join("out2/rpa.json")).unwrap()).unwrap();
    assert_eq!(
        again["settings"]["resumed_q"],
        meta["settings"]["irreducible_q"]
    );
    assert_eq!(again["e_c_hartree_per_cell"], meta["e_c_hartree_per_cell"]);
    // A changed input invalidates the records.
    let changed = case.input("changed", &BASE.replace("points=8", "points=6"));
    let o = case.run(&changed, "out3", Some("ck"));
    assert!(!o.status.success(), "stale checkpoint must reject");
    // Single entry: [task] type = "rpa".
    let task = case.input(
        "task",
        &format!("[task]\ntype='rpa'\n[output]\ndirectory='single'\n{BASE}"),
    );
    let o = case.command().arg("-i").arg(&task).output().unwrap();
    ok(&o);
    assert!(case.0.join("single/rpa.json").is_file());
    // Unsupported settings reject.
    for (label, text) in [
        ("sym", BASE.replace("'time-reversal'", "'spacegroup'")),
        ("dirs", BASE.replace("q0_directions=14", "q0_directions=7")),
        ("order", BASE.replace("[1.5,2.0,2.5]", "[2.5,1.5]")),
        ("unknown", format!("{BASE}extra=1\n")),
    ] {
        let o = case.run(&case.input(label, &text), label, None);
        assert!(!o.status.success(), "{label} must reject");
    }
}
