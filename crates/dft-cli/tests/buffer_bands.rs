//! External read-only validation UPF; never copied into the fixture or repository.
use serde_json::Value;
use std::{fs, process::Command};

/// DFT_SI_UPF=/path/Si.pz-vbc.UPF cargo test --release -p dft-cli --test buffer_bands -- --ignored --nocapture
#[test]
#[ignore = "real-UPF 30 Ry si8 SCF; run in release with DFT_SI_UPF"]
fn si8_strict_certificate_and_deterministic_threads() {
    let upf = std::env::var("DFT_SI_UPF").expect("set DFT_SI_UPF to the read-only validation UPF");
    let dir = std::env::temp_dir().join(format!("dft-si8-buffer-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let mut energies = Vec::new();
    for seed in [1, 7] {
        let mut serial_density = None;
        for threads in [1, 4] {
            let mut text = String::from(
                "[system]\nname='si8 buffer regression'\nlattice_unit='angstrom'\nlattice_vectors=[[5.43,0.0,0.0],[0.0,5.43,0.0],[0.0,0.0,5.43]]\n",
            );
            for site in [
                [0., 0., 0.],
                [0., 0.5, 0.5],
                [0.5, 0., 0.5],
                [0.5, 0.5, 0.],
                [0.25, 0.25, 0.25],
                [0.25, 0.75, 0.75],
                [0.75, 0.25, 0.75],
                [0.75, 0.75, 0.25],
            ] {
                text.push_str(&format!(
                    "[[system.atoms]]\nspecies='Si'\nfractional={site:?}\n"
                ));
            }
            text.push_str(&format!("[pseudopotentials]\nSi={}\n[basis]\necutwfc_ry=30.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nmax_iterations=100\nmixing='pulay'\nmixing_beta=0.5\nenergy_tolerance_ha=1e-8\ndensity_residual_tolerance=1e-7\nseed={seed}\nthreads={threads}\n", serde_json::to_string(&upf).unwrap()));
            let input = dir.join(format!("{seed}-{threads}.toml"));
            let out = dir.join(format!("{seed}-{threads}"));
            fs::write(&input, text).unwrap();
            let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
                .arg("run")
                .arg(input)
                .arg("--output-dir")
                .arg(&out)
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{}",
                String::from_utf8_lossy(&output.stderr)
            );
            let r: Value =
                serde_json::from_slice(&fs::read(out.join("results.json")).unwrap()).unwrap();
            let energy = r["total_energy"]["ha"].as_f64().unwrap();
            assert!((energy - (-31.3218036010316)).abs() < 1e-9);
            assert_eq!(r["bands"], 20);
            assert_eq!(r["convergence"]["converged"], true);
            let accuracy = r["band_accuracy"].as_array().unwrap();
            assert_eq!(accuracy.len(), 20);
            for (j, b) in accuracy.iter().enumerate() {
                assert_eq!(b["buffer"], j >= 17);
                if j < 17 {
                    assert_eq!(b["converged"], true);
                    assert!(b["residual_ha"].as_f64().unwrap() <= 1e-8);
                }
            }
            let history = r["convergence"]["history"].as_array().unwrap();
            let cert = &history.last().unwrap()["fixed_point_certificate"];
            let strict = cert["strict_residuals_per_k_ha"][0].as_array().unwrap();
            assert_eq!(strict.len(), 17);
            assert!(strict.iter().all(|r| r.as_f64().unwrap() <= 1e-9));
            assert_eq!(
                cert["buffer_residuals_per_k_ha"][0]
                    .as_array()
                    .unwrap()
                    .len(),
                3
            );
            let density = fs::read(out.join("density.bin")).unwrap();
            if let Some(serial) = &serial_density {
                assert_eq!(&density, serial);
            } else {
                serial_density = Some(density);
            }
            println!(
                "si8 seed={seed} threads={threads}: E={energy:.17} Ha; SCF={}",
                history.len()
            );
            energies.push(energy);
        }
    }
    assert!(energies.iter().all(|e| (e - energies[0]).abs() < 1e-9));
    fs::remove_dir_all(dir).unwrap();
}
