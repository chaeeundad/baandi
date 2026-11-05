#[allow(dead_code)]
mod common;
use std::{
    fs,
    io::{BufRead, BufReader},
    process::{Command, Stdio},
};
fn input() -> String {
    "[system]\nname='relaxation oracle'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=4\nmax_iterations=300\nenergy_tolerance_ha=1e-12\ndensity_residual_tolerance=1e-10\neigen_tolerance=1e-11\n".into()
}
#[test]
#[ignore = "multi-step SCF, run cargo test --release -p dft-cli --test relaxation -- --ignored --nocapture"]
fn resume_fixed_atoms_formats_and_failure() {
    let dir = std::env::temp_dir().join(format!("p21-relax-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 801).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    fs::write(
        dir.join("input.toml"),
        input() + "\n[relax]\nmax_steps=4\nforce_tol_ha_bohr=1e-8\nfixed=[{atom=0},{atom=1,components=[false,true,false]}]\n",
    )
    .unwrap();
    let mut child = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&dir)
        .args(["relax", "input.toml", "--output-dir", "full"])
        .stdout(Stdio::piped())
        .spawn()
        .unwrap();
    let mut reader = BufReader::new(child.stdout.take().unwrap());
    let mut line = String::new();
    loop {
        line.clear();
        assert!(
            reader.read_line(&mut line).unwrap() > 0,
            "no boundary: {line}"
        );
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
    // Keep draining stdout while the uninterrupted reference finishes.
    for line in reader.lines() {
        println!("{}", line.unwrap());
    }
    assert_eq!(child.wait().unwrap().code(), Some(3));
    let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&dir)
        .args(["relax", "input.toml", "--output-dir", "resumed", "--resume"])
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(3),
        "{}",
        String::from_utf8_lossy(&output.stderr)
    );
    for name in [
        "restart.json",
        "trajectory.extxyz",
        "POSCAR",
        "results.json",
        "final.toml",
        "relax.log",
    ] {
        assert_eq!(
            fs::read(dir.join("full").join(name)).unwrap(),
            fs::read(dir.join("resumed").join(name)).unwrap(),
            "{name}"
        );
    }
    let result: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("full/results.json")).unwrap()).unwrap();
    let frames = result["relaxation"]["frames"].as_array().unwrap();
    for f in frames {
        assert_eq!(f["positions"][0], 0.);
        assert_eq!(f["positions"][1], 0.);
        assert_eq!(f["positions"][2], 0.);
        assert_eq!(f["positions"][4], frames[0]["positions"][4]);
    }
    assert_eq!(frames[1]["provenance"], "warm_transfer");
    let poscar = fs::read_to_string(dir.join("full/POSCAR")).unwrap();
    let lines: Vec<_> = poscar.lines().collect();
    assert_eq!(lines[1], "1.0");
    assert_eq!(lines[5], "H");
    assert_eq!(lines[6], "2");
    assert_eq!(lines[7], "Cartesian");
    let xyz = fs::read_to_string(dir.join("full/trajectory.extxyz")).unwrap();
    assert_eq!(xyz.lines().count(), 4 * frames.len());
    assert!(xyz.contains("Properties=species:S:1:pos:R:3:forces:R:3"));
    fs::write(
        dir.join("failed.toml"),
        input().replace("max_iterations=300", "max_iterations=1") + "\n[relax]\nmax_steps=4\n",
    )
    .unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&dir)
        .args(["relax", "failed.toml", "--output-dir", "failed"])
        .output()
        .unwrap();
    assert_eq!(output.status.code(), Some(4));
    let failure: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("failed/results.json")).unwrap()).unwrap();
    assert_eq!(failure["relaxation"]["status"], "failed");
    assert!(failure["relaxation"]["final_max_force_ha_bohr"].is_null());
    println!(
        "relaxation restart: complete files bit-identical; 4 steps; frozen atom exact; initial SCF failure exit=4"
    );
    fs::remove_dir_all(dir).unwrap();
}

#[test]
#[ignore = "energy scan and two converged trajectories; run in release"]
fn independent_energy_scan_two_starts_and_directional_derivative() {
    use dft_engine::{ScfOptions, ScfProblem, SpeciesPotential};
    use dft_geometry::{Atom, Cell};
    let dir = std::env::temp_dir().join(format!("p21-scan-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let pp_text =
        common::gaussian_with(1., 1.5, false, 801).replace("element=\"Xx\"", "element=\"H\"");
    fs::write(dir.join("syn.xml"), &pp_text).unwrap();
    let pp = dft_upf::parse_upf(&pp_text).unwrap();
    let refs = [SpeciesPotential {
        label: "H",
        potential: &pp,
    }];
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
    let opts = ScfOptions {
        cutoff_hartree: 1.5,
        bands: Some(4),
        density_tolerance: 1e-10,
        energy_tolerance: 1e-12,
        eigen_tolerance: 1e-11,
        max_iterations: 300,
        ..Default::default()
    };
    let energy = |position: [f64; 3]| {
        ScfProblem::new(
            &cell,
            &[Atom::new("H", [0.; 3]), Atom::new("H", position)],
            &refs,
            opts.clone(),
        )
        .unwrap()
        .run()
        .unwrap()
        .total_hartree()
    };
    // Scalar energy samples alone determine the vertex; no force/optimizer call.
    let h = 0.0001;
    let e0 = energy([0.5 - h; 3]);
    let e1 = energy([0.5; 3]);
    let e2 = energy([0.5 + h; 3]);
    assert!(e0 > e1 && e2 > e1, "scan is not a minimum {e0} {e1} {e2}");
    let vertex = 0.5 + 0.5 * h * (e0 - e2) / (e0 - 2. * e1 + e2);
    let scan_distance = 8. * 3f64.sqrt() * vertex;
    println!(
        "scan energies Ha: {e0:.17e} {e1:.17e} {e2:.17e}; fitted distance={scan_distance:.17e} Bohr"
    );
    for (name, t) in [("compressed", 0.4), ("stretched", 0.6)] {
        let text = input().replace("[0.27,0.23,0.19]", &format!("[{t},{t},{t}]"))
            + "\n[relax]\nmax_steps=200\nforce_tol_ha_bohr=1e-8\nenergy_tol_ha=1e-12\nfixed=[{atom=0}]\n";
        fs::write(dir.join(format!("{name}.toml")), text).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .args(["relax", &format!("{name}.toml"), "--output-dir", name])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}\n{}",
            String::from_utf8_lossy(&output.stdout),
            String::from_utf8_lossy(&output.stderr)
        );
        let result: serde_json::Value =
            serde_json::from_slice(&fs::read(dir.join(name).join("results.json")).unwrap())
                .unwrap();
        let frames = result["relaxation"]["frames"].as_array().unwrap();
        let last = frames.last().unwrap();
        let x = last["positions"].as_array().unwrap();
        let distance = x[3..]
            .iter()
            .map(|v| v.as_f64().unwrap().powi(2))
            .sum::<f64>()
            .sqrt();
        println!(
            "{name}: steps={} distance={distance:.17e} Bohr error={:.17e} Bohr maxF={} Ha/Bohr",
            frames.len() - 1,
            (distance - scan_distance).abs(),
            last["max_force_ha_bohr"]
        );
        assert!((distance - scan_distance).abs() < 1e-4);
        for pair in frames.windows(2).take(4) {
            let x: Vec<_> = pair[0]["positions"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect();
            let y: Vec<_> = pair[1]["positions"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect();
            let mut direction: [f64; 3] = std::array::from_fn(|i| y[3 + i] - x[3 + i]);
            let norm = direction.iter().map(|v| v * v).sum::<f64>().sqrt();
            if norm < 1e-10 {
                continue;
            }
            for v in &mut direction {
                *v /= norm;
            }
            let force = &pair[0]["forces_ha_bohr"][1];
            let projection: f64 = (0..3)
                .map(|i| force[i].as_f64().unwrap() * direction[i])
                .sum();
            let mut errors = Vec::new();
            for dh in [0.01, 0.005, 0.0025] {
                let plus = std::array::from_fn(|i| (x[3 + i] + dh * direction[i]) / 8.);
                let minus = std::array::from_fn(|i| (x[3 + i] - dh * direction[i]) / 8.);
                errors.push((-(energy(plus) - energy(minus)) / (2. * dh) - projection).abs());
            }
            println!("{name}: directional FD errors={errors:?} Ha/Bohr");
            assert!(errors[2] < 1e-6);
            assert!(errors[2] < 0.3 * errors[0] || errors.iter().all(|e| *e < 1e-9));
        }
    }
    fs::remove_dir_all(dir).unwrap();
}

#[test]
#[ignore = "multiple symmetry-projected SCFs; run in release"]
fn symmetry_preserves_axis_and_space_group() {
    let dir = std::env::temp_dir().join(format!("p21-sym-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 801).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    fs::write(
        dir.join("input.toml"),
        input().replace("[0.27,0.23,0.19]", "[0.25,0.0,0.0]")
            + "\n[symmetry]\nmode='auto'\n[relax]\nmax_steps=4\n",
    )
    .unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&dir)
        .args(["relax", "input.toml", "--output-dir", "out"])
        .output()
        .unwrap();
    assert_eq!(
        output.status.code(),
        Some(3),
        "{} {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let result: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
    let cell = dft_geometry::Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]])
        .unwrap();
    let mut counts = Vec::new();
    for frame in result["relaxation"]["frames"].as_array().unwrap() {
        let xyz = frame["positions"].as_array().unwrap();
        for j in [1, 2, 4, 5] {
            assert!(xyz[j].as_f64().unwrap().abs() < 1e-14); /* Cartesian averaging roundoff, Bohr. */
        }
        let atoms: Vec<_> = frame["fractional"]
            .as_array()
            .unwrap()
            .iter()
            .map(|a| dft_geometry::Atom::new("H", std::array::from_fn(|i| a[i].as_f64().unwrap())))
            .collect();
        counts.push(
            dft_geometry::symmetry::detect(&cell, &atoms, 1e-7)
                .unwrap()
                .len(),
        );
    }
    assert!(counts.iter().all(|n| *n == counts[0]));
    println!("symmetry operation counts over trajectory: {counts:?}");
    fs::remove_dir_all(dir).unwrap();
}

#[test]
#[ignore = "external PseudoDojo Si/MgO relaxations; set DFT_NC_DIR; run release"]
fn real_si_and_mgo() {
    let ppdir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    for material in ["Si", "MgO"] {
        let dir = std::env::temp_dir().join(format!("p21-real-{material}-{}", std::process::id()));
        fs::create_dir_all(&dir).unwrap();
        let (a, t, cutoff, labels) = if material == "Si" {
            (10.2612, 0.25, 60., ["Si", "Si"])
        } else {
            (7.96, 0.5, 120., ["Mg", "O"])
        };
        let p0 = ppdir.join(format!("{}.upf", labels[0]));
        let p1 = ppdir.join(format!("{}.upf", labels[1]));
        let mut pp = format!("{}={:?}\n", labels[0], p0);
        if material != "Si" {
            pp += &format!("{}={:?}\n", labels[1], p1);
        }
        let half = a / 2.;
        let displaced = t + 0.03 / a;
        let text = format!(
            "[system]\nname='{material}'\nlattice_unit='bohr'\nlattice_vectors=[[0.0,{half},{half}],[{half},0.0,{half}],[{half},{half},0.0]]\n[[system.atoms]]\nspecies='{}'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='{}'\nfractional=[{displaced},{displaced},{displaced}]\n[pseudopotentials]\n{pp}[basis]\necutwfc_ry={cutoff}\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='pbe'\ndensity_floor=1e-10\n[scf]\nmax_iterations=300\ndensity_residual_tolerance=1e-9\nenergy_tolerance_ha=1e-11\neigen_tolerance=1e-10\n[relax]\nmax_steps=80\nforce_tol_ha_bohr=1e-4\nenergy_tol_ha=1e-9\nfixed=[{{atom=0}}]\n",
            labels[0], labels[1]
        );
        fs::write(dir.join("input.toml"), text).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .args(["relax", "input.toml", "--output-dir", "out"])
            .output()
            .unwrap();
        println!("{material}: {}", String::from_utf8_lossy(&output.stdout));
        assert!(
            output.status.success(),
            "{material}: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let result: serde_json::Value =
            serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
        let r = &result["relaxation"];
        assert!(r["final_max_force_ha_bohr"].as_f64().unwrap() < 1e-4);
        let energies = r["energy_history_ha"].as_array().unwrap();
        for pair in energies.windows(2) {
            assert!(
                pair[1].as_f64().unwrap() <= pair[0].as_f64().unwrap() + 1e-9,
                "nonmonotone {material}"
            );
        }
        let last = r["frames"].as_array().unwrap().last().unwrap();
        println!(
            "{material} final fractional coordinates: {}",
            last["fractional"]
        );
        // Si's ideal-bond oracle is required by the ticket. Gamma-only MgO
        // has an unstable ideal structure in this finite model; its oracle is
        // convergence and monotonic lowering, not bulk structural qualification.
        if material == "Si" {
            for i in 0..3 {
                assert!(
                    (last["fractional"][1][i].as_f64().unwrap() - t).abs() * a < 0.01,
                    "ideal-position recovery Si"
                );
            }
        }
        if material == "Si" {
            let failure_input = fs::read_to_string(dir.join("input.toml"))
                .unwrap()
                .replace("max_iterations=300", "max_iterations=17");
            fs::write(dir.join("failure.toml"), failure_input).unwrap();
            let failed = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
                .current_dir(&dir)
                .args(["relax", "failure.toml", "--output-dir", "failed"])
                .output()
                .unwrap();
            assert_eq!(
                failed.status.code(),
                Some(4),
                "{} {}",
                String::from_utf8_lossy(&failed.stdout),
                String::from_utf8_lossy(&failed.stderr)
            );
            let failed: serde_json::Value =
                serde_json::from_slice(&fs::read(dir.join("failed/results.json")).unwrap())
                    .unwrap();
            let r = &failed["relaxation"];
            assert_eq!(r["status"], "failed");
            assert!(r["steps"].as_u64().unwrap() > 0);
            assert_eq!(r["failures"].as_array().unwrap().len(), 2);
            let f = r["frames"].as_array().unwrap().last().unwrap();
            assert_eq!(r["final_max_force_ha_bohr"], f["max_force_ha_bohr"]);
            println!(
                "Si budget-17 failure retained step {} with two rejected attempts",
                r["steps"]
            );
        }
        fs::remove_dir_all(dir).unwrap();
    }
}

#[test]
#[ignore = "FIRE trajectory with many small ionic steps; run release"]
fn fire_converges_from_displaced_structure() {
    let dir = std::env::temp_dir().join(format!("p21-fire-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 801).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    fs::write(dir.join("input.toml"),input().replace("[0.27,0.23,0.19]","[0.4,0.4,0.4]")+"\n[relax]\nmethod='fire'\nmax_steps=2000\nforce_tol_ha_bohr=1e-7\nenergy_tol_ha=1e-10\nfixed=[{atom=0}]\n").unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&dir)
        .args(["relax", "input.toml", "--output-dir", "out"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{} {}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let result: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
    let r = &result["relaxation"];
    println!(
        "FIRE: steps={} maxF={} Ha/Bohr",
        r["steps"], r["final_max_force_ha_bohr"]
    );
    assert!(r["final_max_force_ha_bohr"].as_f64().unwrap() < 1e-7);
    fs::remove_dir_all(dir).unwrap();
}

#[test]
fn rejects_standard_profile_and_invalid_constraints_without_scf() {
    let dir = std::env::temp_dir().join(format!("p21-invalid-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    for (i, text) in [
        input().replace("[scf]", "[scf]\nconvergence_profile='standard'") + "\n[relax]\n",
        input() + "\n[relax]\nfixed=[{atom=2}]\n",
        input() + "\n[relax]\nforce_tol_ha_bohr=nan\n",
        input().replace("[0.27,0.23,0.19]", "[0.25,0.0,0.0]")
            + "\n[symmetry]\nmode='auto'\n[relax]\nfixed=[{atom=0}]\n",
    ]
    .iter()
    .enumerate()
    {
        fs::write(dir.join("input.toml"), text).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .args(["relax", "input.toml", "--output-dir", &format!("out{i}")])
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(1),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        assert!(!dir.join(format!("out{i}/results.json")).exists());
    }
    fs::remove_dir_all(dir).unwrap();
}
