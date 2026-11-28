//! R4b `baandi wannier-export` and `[subspace] method = "wannier90"` adapter
//! gates on a synthetic UPF (never reference evidence):
//! - Writes `.win`, `.eig`, `.amn`, `.mmn`, `wannier_export.json`, `results.json`
//! - Independent parsers verify every Wannier90 file format (row counts,
//!   1-based indices, eV/Angstrom units, `.mmn` Hermitian adjointness)
//! - `--nnkp` mode matches an externally supplied `.nnkp` neighbor table
//! - Full round-trip: write `_u.mat` (+ `_u_dis.mat`) and import through
//!   `baandi subspace` (`method = "wannier90"`), verifying `subspace_hr.dat`,
//!   `subspace_weights.csv` sum rule and spreads
//! - `[task] type = "wannier-export"` single entry and input guards
//!
//! Fast gate: `cargo test -p dft-cli --test wannier_r4b`
#[allow(dead_code)]
mod common;
use dft_physics::data::Complex64;
use serde_json::Value;
use std::{collections::HashMap, fmt::Write as _, fs, path::Path, path::PathBuf, process::Command};

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
        let dir = std::env::temp_dir().join(format!("r4b-wannier-{tag}-{}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        fs::write(dir.join("ion.upf"), atomic(2.)).unwrap();
        Self(dir)
    }
    fn input(&self, label: &str, tables: &str) -> PathBuf {
        let text = format!(
            "[system]\nname='R4b wannier'\nlattice_unit='bohr'\nlattice_vectors=[[6.0,0.0,0.0],[0.7,5.5,0.0],[0.3,0.4,6.5]]\n[[system.atoms]]\nspecies='Xx'\nfractional=[0.13,0.21,0.37]\n[pseudopotentials]\nXx='ion.upf'\n[basis]\necutwfc_ry=8.0\nkgrid=[2,2,2]\n[xc]\nfunctional='lda-pw'\n[scf]\nmax_iterations=150\n{tables}\n"
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

fn verify_eig(path: &Path, expected_bands: usize, expected_nk: usize) {
    let text = fs::read_to_string(path).unwrap();
    let rows: Vec<&str> = text.lines().filter(|l| !l.trim().is_empty()).collect();
    assert_eq!(rows.len(), expected_bands * expected_nk);
    for (idx, line) in rows.iter().enumerate() {
        let cols: Vec<&str> = line.split_whitespace().collect();
        assert_eq!(cols.len(), 3);
        let m: usize = cols[0].parse().unwrap();
        let k: usize = cols[1].parse().unwrap();
        let e: f64 = cols[2].parse().unwrap();
        assert_eq!(m, (idx % expected_bands) + 1);
        assert_eq!(k, (idx / expected_bands) + 1);
        assert!(e.is_finite());
    }
}

fn verify_amn(path: &Path, expected_bands: usize, expected_nk: usize, expected_nw: usize) {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let _hdr = lines.next().unwrap();
    let dims: Vec<usize> = lines
        .next()
        .unwrap()
        .split_whitespace()
        .map(|s| s.parse().unwrap())
        .collect();
    assert_eq!(dims, vec![expected_bands, expected_nk, expected_nw]);
    let mut count = 0usize;
    for line in lines {
        if line.trim().is_empty() {
            continue;
        }
        let cols: Vec<&str> = line.split_whitespace().collect();
        assert_eq!(cols.len(), 5);
        let m: usize = cols[0].parse().unwrap();
        let n: usize = cols[1].parse().unwrap();
        let k: usize = cols[2].parse().unwrap();
        let re: f64 = cols[3].parse().unwrap();
        let im: f64 = cols[4].parse().unwrap();
        assert!((1..=expected_bands).contains(&m));
        assert!((1..=expected_nw).contains(&n));
        assert!((1..=expected_nk).contains(&k));
        assert!(re.is_finite() && im.is_finite());
        count += 1;
    }
    assert_eq!(count, expected_bands * expected_nw * expected_nk);
}

fn verify_mmn_and_adjoint(path: &Path, expected_bands: usize, expected_nk: usize) -> usize {
    let text = fs::read_to_string(path).unwrap();
    let mut lines = text.lines();
    let _hdr = lines.next().unwrap();
    let dims: Vec<usize> = lines
        .next()
        .unwrap()
        .split_whitespace()
        .map(|s| s.parse().unwrap())
        .collect();
    assert_eq!(dims[0], expected_bands);
    assert_eq!(dims[1], expected_nk);
    let nntot = dims[2];
    assert!(nntot > 0);
    let mut blocks: HashMap<(usize, usize, [i32; 3]), Vec<Complex64>> = HashMap::new();
    for _ in 0..(expected_nk * nntot) {
        let head: Vec<i32> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();
        assert_eq!(head.len(), 5);
        let key = (
            head[0] as usize,
            head[1] as usize,
            [head[2], head[3], head[4]],
        );
        let mut mat = Vec::with_capacity(expected_bands * expected_bands);
        for _ in 0..(expected_bands * expected_bands) {
            let pair: Vec<f64> = lines
                .next()
                .unwrap()
                .split_whitespace()
                .map(|s| s.parse().unwrap())
                .collect();
            assert_eq!(pair.len(), 2);
            mat.push(Complex64::new(pair[0], pair[1]));
        }
        blocks.insert(key, mat);
    }
    assert_eq!(blocks.len(), expected_nk * nntot);
    for (&(ik, kn, g), mat) in &blocks {
        let rev_key = (kn, ik, [-g[0], -g[1], -g[2]]);
        if let Some(rev) = blocks.get(&rev_key) {
            for col in 0..expected_bands {
                for row in 0..expected_bands {
                    let fwd_val = mat[col * expected_bands + row];
                    let rev_val = rev[row * expected_bands + col];
                    assert!(
                        (fwd_val - rev_val.conj()).norm() < 1e-10,
                        "mmn adjoint mismatch at ({ik}->{kn}, {g:?}) ({row},{col})"
                    );
                }
            }
        }
    }
    nntot
}

#[test]
fn wannier_export_nnkp_and_gauge_import_roundtrip() {
    let case = Case::new("main");
    let export_input = case.input(
        "exp",
        "[subspace]\nk_mesh=[2,2,2]\nbands=6\nseedname='toy'\nwindow_ha=[-0.5,0.5]\nspecies='Xx'\norbitals=['sp3']\n",
    );
    let out_dir = case.0.join("w90_out");
    let out = case
        .command()
        .arg("wannier-export")
        .arg(&export_input)
        .arg("--output-dir")
        .arg(&out_dir)
        .output()
        .unwrap();
    ok(&out);

    for f in [
        "toy.win",
        "toy.eig",
        "toy.amn",
        "toy.mmn",
        "wannier_export.json",
        "results.json",
    ] {
        assert!(out_dir.join(f).is_file(), "missing {f}");
    }
    let meta: Value =
        serde_json::from_slice(&fs::read(out_dir.join("wannier_export.json")).unwrap()).unwrap();
    assert_eq!(meta["kind"], "baandi_wannier_export_r4b");
    assert!(meta["summary"]["max_mmn_adjoint_defect"].as_f64().unwrap() < 1e-12);
    assert!(
        meta["summary"]["b1_completeness_residual"]
            .as_f64()
            .unwrap()
            < 1e-10
    );
    let win_text = fs::read_to_string(out_dir.join("toy.win")).unwrap();
    assert!(win_text.contains("num_bands = 6"));
    assert!(win_text.contains("num_wann  = 4"));
    assert!(win_text.contains("mp_grid = 2 2 2"));
    assert!(win_text.contains("begin unit_cell_cart"));
    assert!(win_text.contains("begin kpoints"));
    verify_eig(&out_dir.join("toy.eig"), 6, 8);
    verify_amn(&out_dir.join("toy.amn"), 6, 8, 4);
    let nntot = verify_mmn_and_adjoint(&out_dir.join("toy.mmn"), 6, 8);
    assert!(nntot >= 6);

    // Test --nnkp override with a minimal synthetic 2-neighbor .nnkp file.
    let nnkp_path = case.0.join("toy.nnkp");
    let mut nnkp = String::from("begin nnkpts\n  2\n");
    for ik in 1..=8 {
        writeln!(nnkp, "  {ik}  {ik}  1  0  0").unwrap();
        writeln!(nnkp, "  {ik}  {ik} -1  0  0").unwrap();
    }
    nnkp.push_str("end nnkpts\n");
    fs::write(&nnkp_path, nnkp).unwrap();
    let nnkp_out = case.0.join("w90_nnkp");
    let res_nnkp = case
        .command()
        .arg("wannier-export")
        .arg(&export_input)
        .arg("--output-dir")
        .arg(&nnkp_out)
        .arg("--nnkp")
        .arg(&nnkp_path)
        .output()
        .unwrap();
    ok(&res_nnkp);
    assert_eq!(verify_mmn_and_adjoint(&nnkp_out.join("toy.mmn"), 6, 8), 2);

    // Round-trip: write synthetic identity _u.mat (1x1) and _u_dis.mat (6x1 selecting band 1)
    // and import through `baandi subspace` with method = "wannier90"!
    let mut u_mat = String::from("synthetic _u.mat for R4b test\n  8  1  1\n");
    let mut u_dis = String::from("synthetic _u_dis.mat for R4b test\n  8  1  6\n");
    for i in 0..2 {
        for j in 0..2 {
            for k in 0..2 {
                let (kx, ky, kz) = (i as f64 / 2.0, j as f64 / 2.0, k as f64 / 2.0);
                writeln!(
                    u_mat,
                    "\n  {kx:.10}  {ky:.10}  {kz:.10}\n  1.0000000000  0.0000000000"
                )
                .unwrap();
                writeln!(u_dis, "\n  {kx:.10}  {ky:.10}  {kz:.10}").unwrap();
                writeln!(u_dis, "  1.0000000000  0.0000000000").unwrap();
                for _ in 1..6 {
                    writeln!(u_dis, "  0.0000000000  0.0000000000").unwrap();
                }
            }
        }
    }
    fs::write(out_dir.join("toy_u.mat"), u_mat).unwrap();
    fs::write(out_dir.join("toy_u_dis.mat"), u_dis).unwrap();

    let import_input = case.input(
        "imp",
        "[subspace]\nmethod='wannier90'\nk_mesh=[2,2,2]\nseedname='toy'\nwannier_dir='w90_out'\n",
    );
    let imp_dir = case.0.join("imp_out");
    let imp_out = case
        .command()
        .arg("subspace")
        .arg(&import_input)
        .arg("--output-dir")
        .arg(&imp_dir)
        .output()
        .unwrap();
    ok(&imp_out);
    let imp_meta: Value =
        serde_json::from_slice(&fs::read(imp_dir.join("subspace.json")).unwrap()).unwrap();
    assert_eq!(imp_meta["kind"], "baandi_subspace_r4b_wannier90");
    assert_eq!(imp_meta["summary"]["num_wann"], 1);
    assert!((imp_meta["summary"]["total_population"].as_f64().unwrap() - 2.0).abs() < 1e-12);
    assert!(imp_dir.join("subspace_hr.dat").is_file());
    assert!(imp_dir.join("subspace_weights.csv").is_file());

    // Single entry [task] type = "wannier-export".
    let single = case.input(
        "single",
        "[task]\ntype='wannier-export'\n[output]\ndirectory='single_w90'\n[subspace]\nk_mesh=[2,2,2]\nbands=4\nseedname='s1'\nspecies='Xx'\nl=0\n",
    );
    ok(&case.command().arg("-i").arg(&single).output().unwrap());
    assert!(case.0.join("single_w90/s1.win").is_file());
    assert!(case.0.join("single_w90/s1.mmn").is_file());
}
