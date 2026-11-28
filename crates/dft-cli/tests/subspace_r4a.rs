//! R4a `baandi subspace` adapter gates on a synthetic UPF (never reference
//! evidence): outputs written, independent parser check of `subspace_hr.dat`
//! against the public Wannier90 format specification (WS degeneracies sum rule,
//! Hermiticity, exact on-mesh eigenvalue recovery in eV), CSV/JSON consistency,
//! `[task] type = "subspace"` single entry, and input guards.
//!
//! Fast gate: `cargo test -p dft-cli --test subspace_r4a`
#[allow(dead_code)]
mod common;
use dft_physics::data::Complex64;
use dft_units::{EnergyUnit, from_hartree};
use serde_json::Value;
use std::{
    collections::BTreeMap, f64::consts::TAU, fs, path::Path, path::PathBuf, process::Command,
};

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
    pp.replace("number_of_wfc=\"0\"", "number_of_wfc=\"2\"")
        .replace(
            "</UPF>",
            &format!(
                "<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1><PP_CHI.2 l=\"1\" label=\"p\" occupation=\"0\">{}</PP_CHI.2></PP_PSWFC></UPF>",
                array(&|r| r * (-0.6 * r * r).exp()),
                array(&|r| r * r * (-0.5 * r * r).exp())
            ),
        )
}

struct Case(PathBuf);
impl Drop for Case {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}
impl Case {
    fn new(tag: &str) -> Self {
        let dir = std::env::temp_dir().join(format!("r4a-subspace-{tag}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("ion.upf"), atomic(2.)).unwrap();
        Self(dir)
    }
    fn input(&self, label: &str, tables: &str) -> PathBuf {
        let text = format!(
            "[system]\nname='R4a subspace'\nlattice_unit='bohr'\nlattice_vectors=[[6.0,0.0,0.0],[0.7,5.5,0.0],[0.3,0.4,6.5]]\n[[system.atoms]]\nspecies='Xx'\nfractional=[0.13,0.21,0.37]\n[pseudopotentials]\nXx='ion.upf'\n[basis]\necutwfc_ry=8.0\nkgrid=[2,2,2]\n[xc]\nfunctional='lda-pw'\n[scf]\nmax_iterations=150\n{tables}\n"
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
}

fn ok(o: &std::process::Output) {
    assert!(
        o.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&o.stdout),
        String::from_utf8_lossy(&o.stderr)
    );
}

struct ParsedHr {
    num_wann: usize,
    nrpts: usize,
    degs: Vec<usize>,
    mats: BTreeMap<[i32; 3], Vec<Complex64>>,
}

/// Independent parser of Wannier90 `seedname_hr.dat` per public format spec.
fn parse_hr_dat(path: &Path) -> ParsedHr {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let _header = lines.next().expect("header line");
    let num_wann: usize = lines.next().unwrap().trim().parse().unwrap();
    let nrpts: usize = lines.next().unwrap().trim().parse().unwrap();
    let mut degs = Vec::with_capacity(nrpts);
    while degs.len() < nrpts {
        let line = lines.next().expect("degeneracy line");
        for tok in line.split_whitespace() {
            degs.push(tok.parse::<usize>().unwrap());
        }
    }
    assert_eq!(degs.len(), nrpts);
    let mut matrices: BTreeMap<[i32; 3], Vec<Complex64>> = BTreeMap::new();
    let mut row_count = 0usize;
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split_whitespace().collect();
        assert_eq!(cols.len(), 7, "bad hr.dat row: {line}");
        let r = [
            cols[0].parse::<i32>().unwrap(),
            cols[1].parse::<i32>().unwrap(),
            cols[2].parse::<i32>().unwrap(),
        ];
        let m: usize = cols[3].parse::<usize>().unwrap();
        let n: usize = cols[4].parse::<usize>().unwrap();
        assert!((1..=num_wann).contains(&m) && (1..=num_wann).contains(&n));
        let re: f64 = cols[5].parse().unwrap();
        let im: f64 = cols[6].parse().unwrap();
        let entry = matrices
            .entry(r)
            .or_insert_with(|| vec![Complex64::default(); num_wann * num_wann]);
        entry[(n - 1) * num_wann + (m - 1)] = Complex64::new(re, im);
        row_count += 1;
    }
    assert_eq!(row_count, nrpts * num_wann * num_wann);
    assert_eq!(matrices.len(), nrpts);
    ParsedHr {
        num_wann,
        nrpts,
        degs,
        mats: matrices,
    }
}

#[test]
fn subspace_outputs_hr_dat_single_entry_and_guards() {
    let case = Case::new("main");
    let sub_table = "[subspace]\nk_mesh=[2,2,2]\nbands=8\nband_range=[1,1]\nspecies='Xx'\nl=0\n";
    let input = case.input("isolated", sub_table);
    let out_dir = case.0.join("out");
    let out = case
        .command()
        .arg("subspace")
        .arg(&input)
        .arg("--output-dir")
        .arg(&out_dir)
        .output()
        .unwrap();
    ok(&out);

    for f in [
        "subspace.json",
        "subspace_hr.dat",
        "subspace_populations.csv",
        "subspace_weights.csv",
        "results.json",
    ] {
        assert!(out_dir.join(f).is_file(), "missing {f}");
    }

    let meta: Value =
        serde_json::from_slice(&fs::read(out_dir.join("subspace.json")).unwrap()).unwrap();
    assert_eq!(meta["kind"], "baandi_subspace_r4a");
    assert_eq!(meta["summary"]["num_wann"], 1);
    assert!(
        meta["summary"]["max_orthonormality_defect"]
            .as_f64()
            .unwrap()
            < 1e-12
    );
    assert_eq!(meta["projector_fingerprint"].as_str().unwrap().len(), 64);

    // Independent parse of subspace_hr.dat: WS partition of unity and on-site value vs JSON.
    let ParsedHr {
        num_wann,
        nrpts,
        degs,
        mats,
    } = parse_hr_dat(&out_dir.join("subspace_hr.dat"));
    assert_eq!(num_wann, 1);
    let sum_inv_g: f64 = degs.iter().map(|&g| 1.0 / g as f64).sum();
    assert!((sum_inv_g - 8.0).abs() < 1e-12, "sum 1/g = {sum_inv_g}");
    let h0_json_ha = meta["on_site_hamiltonian_ha"][0][0][0].as_f64().unwrap();
    let h0_hr_ev = mats[&[0, 0, 0]][0].re;
    assert!(
        (h0_hr_ev - from_hartree(h0_json_ha, EnergyUnit::ElectronVolt)).abs() < 1e-7,
        "hr.dat on-site {h0_hr_ev} eV vs json {}",
        from_hartree(h0_json_ha, EnergyUnit::ElectronVolt)
    );
    // Check H(R) = H(-R)^* in parsed hr.dat and Gamma interpolation equals band 1 energy.
    let mut deg_by_r: BTreeMap<[i32; 3], usize> = BTreeMap::new();
    for (r, &g) in mats.keys().zip(&degs) {
        deg_by_r.insert(*r, g);
    }
    let mut h_gamma = Complex64::default();
    for (r, mat) in &mats {
        let neg = [-r[0], -r[1], -r[2]];
        assert!((mat[0] - mats[&neg][0].conj()).norm() < 1e-7);
        let phase = Complex64::from_polar(1.0 / deg_by_r[r] as f64, TAU * 0.0);
        h_gamma += mat[0] * phase;
    }
    // Compare with Gamma band 1 in subspace_weights.csv.
    let weights_text = fs::read_to_string(out_dir.join("subspace_weights.csv")).unwrap();
    let gamma_band1_ha: f64 = weights_text
        .lines()
        .filter(|l| !l.starts_with('#'))
        .nth(1)
        .unwrap()
        .split(',')
        .nth(5)
        .unwrap()
        .parse()
        .unwrap();
    assert!(
        (h_gamma.re - from_hartree(gamma_band1_ha, EnergyUnit::ElectronVolt)).abs() < 1e-6,
        "hr.dat Gamma {} eV vs KS {} eV (nrpts={nrpts})",
        h_gamma.re,
        from_hartree(gamma_band1_ha, EnergyUnit::ElectronVolt)
    );

    // Single entry with [task] type = "subspace" + multi-projection sp3 + [bands] interpolation.
    let single_input = case.input(
        "single",
        "[task]\ntype='subspace'\n[output]\ndirectory='single'\n[subspace]\nk_mesh=[2,2,2]\nbands=8\n[[subspace.projections]]\nspecies='Xx'\norbitals=['sp3']\n[bands]\npoints_per_segment=4\n[[bands.points]]\nlabel='G'\nfractional=[0.0,0.0,0.0]\n[[bands.points]]\nlabel='X'\nfractional=[0.5,0.0,0.0]\n",
    );
    let single_out = case
        .command()
        .arg("-i")
        .arg(&single_input)
        .output()
        .unwrap();
    ok(&single_out);
    assert!(case.0.join("single/subspace.json").is_file());
    assert!(case.0.join("single/subspace_hr.dat").is_file());
    assert!(case.0.join("single/subspace_bands.csv").is_file());
    let single_meta: Value =
        serde_json::from_slice(&fs::read(case.0.join("single/subspace.json")).unwrap()).unwrap();
    assert_eq!(single_meta["summary"]["num_wann"], 4);
    assert_eq!(single_meta["summary"]["interpolated_band_kpoints"], 5);

    // Guard: unknown method rejects cleanly.
    let bad = case.input(
        "bad",
        "[subspace]\nmethod='bogus'\nk_mesh=[2,2,2]\nspecies='Xx'\nl=0\n",
    );
    let bad_out = case
        .command()
        .arg("subspace")
        .arg(&bad)
        .arg("--output-dir")
        .arg(case.0.join("bad_out"))
        .output()
        .unwrap();
    assert!(!bad_out.status.success());
}
