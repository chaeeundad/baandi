//! C3 ionic-step start extrapolation: relaxed path identity against the plain
//! warm transfer, SCF iterations per step, resume identity with the history.
//! Multi-step SCF, release only:
//! `cargo test --release -p dft-cli --test relax_extrapolation -- --ignored --nocapture`
#[allow(dead_code)]
mod common;
use std::{
    fs,
    io::{BufRead, BufReader},
    process::{Command, Stdio},
};

/// Synthetic Gaussian ion with a Gaussian atomic charge (normalised to z), so
/// that the atomic superposition n_at exists.
fn potential_with_charge(z: f64) -> String {
    let text = common::gaussian_with(z, 1.5, false, 801).replace("element=\"Xx\"", "element=\"H\"");
    let r: Vec<f64> = text
        .split("<PP_R>")
        .nth(1)
        .unwrap()
        .split("</PP_R>")
        .next()
        .unwrap()
        .split_whitespace()
        .map(|v| v.parse().unwrap())
        .collect();
    let a = 1.2f64;
    let rho: Vec<String> = r
        .iter()
        .map(|r| {
            let v = z
                * 4.
                * std::f64::consts::PI
                * r
                * r
                * (a / std::f64::consts::PI).powf(1.5)
                * (-a * r * r).exp();
            format!("{v:.17e}")
        })
        .collect();
    let start = text.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
    let end = text.find("</PP_RHOATOM>").unwrap();
    format!("{}{}{}", &text[..start], rho.join(" "), &text[end..])
}

fn input(extra: &str) -> String {
    format!(
        "[system]\nname='c3'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,9.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.2,0.17,0.15]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=4\nmax_iterations=300\nenergy_tolerance_ha=1e-12\ndensity_residual_tolerance=1e-10\neigen_tolerance=1e-11\n[relax]\nmax_steps=5\nforce_tol_ha_bohr=1e-9\n{extra}"
    )
}

fn frames(dir: &std::path::Path, out: &str) -> Vec<serde_json::Value> {
    let v: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join(out).join("results.json")).unwrap()).unwrap();
    v["relaxation"]["frames"].as_array().unwrap().clone()
}

#[test]
#[ignore = "multi-step SCF relaxations; run in release"]
fn extrapolated_starts_keep_the_relaxation_path_and_resume_exactly() {
    let dir = std::env::temp_dir().join(format!("c3-relax-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(dir.join("syn.xml"), potential_with_charge(1.)).unwrap();
    let mut runs = Vec::new();
    for mode in ["none", "first_order", "second_order"] {
        let name = format!("{mode}.toml");
        fs::write(
            dir.join(&name),
            input(&format!("start_extrapolation='{mode}'\n")),
        )
        .unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .env_remove("BAANDI_RELAX_EXTRAPOLATION")
            .args(["relax", &name, "--output-dir", mode])
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(3),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        runs.push((mode, frames(&dir, mode)));
    }
    let (_, reference) = &runs[0];
    for (mode, f) in &runs {
        assert_eq!(f.len(), reference.len(), "{mode}");
        let scf: Vec<_> = f
            .iter()
            .map(|f| f["scf_iterations"].as_u64().unwrap())
            .collect();
        let h: Vec<_> = f
            .iter()
            .map(|f| f["scf_h_vectors"].as_u64().unwrap())
            .collect();
        let kinds: Vec<_> = f
            .iter()
            .map(|f| f["provenance"].as_str().unwrap())
            .collect();
        let mut max_de = 0f64;
        let mut max_df = 0f64;
        for (a, b) in f.iter().zip(reference) {
            max_de = max_de
                .max((a["energy_ha"].as_f64().unwrap() - b["energy_ha"].as_f64().unwrap()).abs());
            for (x, y) in a["forces_ha_bohr"]
                .as_array()
                .unwrap()
                .iter()
                .flat_map(|v| v.as_array().unwrap())
                .zip(
                    b["forces_ha_bohr"]
                        .as_array()
                        .unwrap()
                        .iter()
                        .flat_map(|v| v.as_array().unwrap()),
                )
            {
                max_df = max_df.max((x.as_f64().unwrap() - y.as_f64().unwrap()).abs());
            }
        }
        let last = |f: &[serde_json::Value]| -> Vec<f64> {
            f.last().unwrap()["positions"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect()
        };
        let max_dx = last(f)
            .iter()
            .zip(last(reference))
            .fold(0f64, |m, (a, b)| m.max((a - b).abs()));
        println!(
            "{mode}: SCF per step {scf:?}, H vectors {h:?}, starts {kinds:?}; max |dE| {max_de:e} Ha, max |dF| {max_df:e} Ha/Bohr, final max |dx| {max_dx:e} Bohr"
        );
        // Strict-class identity of the physics along the path (design note C3).
        assert!(max_de < 1e-8, "{mode}: {max_de:e}");
        assert!(max_df < 1e-6, "{mode}: {max_df:e}");
        assert!(max_dx < 1e-5, "{mode}: {max_dx:e}");
        if *mode != "none" {
            assert!(kinds[1..].iter().all(|k| k.starts_with("extrapolated_")));
        }
    }
    assert_eq!(runs[2].1[2]["provenance"], "extrapolated_second_order");

    // Resume after step 1 reproduces the uninterrupted second-order run exactly
    // (the predictor history is part of restart.json).
    let mut child = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&dir)
        .env_remove("BAANDI_RELAX_EXTRAPOLATION")
        .args(["relax", "second_order.toml", "--output-dir", "full"])
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut reader = BufReader::new(child.stdout.take().unwrap());
    let mut line = String::new();
    loop {
        line.clear();
        assert!(reader.read_line(&mut line).unwrap() > 0, "no boundary");
        if line.starts_with("step 1 ") {
            break;
        }
    }
    fs::create_dir_all(dir.join("resumed")).unwrap();
    fs::copy(
        dir.join("full/restart.json"),
        dir.join("resumed/restart.json"),
    )
    .unwrap();
    for line in reader.lines() {
        line.unwrap();
    }
    assert_eq!(child.wait().unwrap().code(), Some(3));
    let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&dir)
        .env_remove("BAANDI_RELAX_EXTRAPOLATION")
        .args([
            "relax",
            "second_order.toml",
            "--output-dir",
            "resumed",
            "--resume",
        ])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(3));
    for name in ["restart.json", "results.json", "relax.log"] {
        assert_eq!(
            fs::read(dir.join("full").join(name)).unwrap(),
            fs::read(dir.join("resumed").join(name)).unwrap(),
            "{name}"
        );
    }
    fs::remove_dir_all(dir).unwrap();
}
