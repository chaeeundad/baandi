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
#[ignore = "multi-step SCF, run cargo test --release -p dft-cli --test cell_relaxation -- --ignored --nocapture"]
fn cell_resume_monotonicity_and_protocols() {
    let dir = std::env::temp_dir().join(format!("p23-cell-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 801).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    fs::write(
        dir.join("input.toml"),
        input() + "\n[relax]\nmode='ions_cell'\ncell_dofs='isotropic'\nmax_steps=4\nforce_tol_ha_bohr=1e-8\n",
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
        "final.cif",
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
    let enthalpies = result["relaxation"]["enthalpy_history_ha"]
        .as_array()
        .unwrap();
    for pair in enthalpies.windows(2) {
        assert!(pair[1].as_f64().unwrap() <= pair[0].as_f64().unwrap() + 1e-8);
    }
    assert!(result["relaxation"]["rebuilt_final"]["pressure_ha_bohr3"].is_number());
    assert!(
        result["relaxation"]["basis_rebuild_difference"]["maximum_stress_difference_ha_bohr3"]
            .is_number()
    );
    assert_ne!(
        frames[0]["cell_columns_bohr"],
        frames[4]["cell_columns_bohr"]
    );
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
}

#[test]
#[ignore = "independent frozen-basis volume scans; run in release mode"]
fn cubic_scan_and_pressure() {
    use dft_engine::{ScfProblem, SpeciesPotential};
    use dft_geometry::{Atom, Cell};
    let dir = std::env::temp_dir().join(format!("p23-scan-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let xml = common::gaussian_with(1., 1.5, false, 801).replace("element=\"Xx\"", "element=\"H\"");
    fs::write(dir.join("syn.xml"), &xml).unwrap();
    let pp = dft_upf::parse_upf(&xml).unwrap();
    let atoms = [Atom::new("H", [0.; 3]), Atom::new("H", [0.5; 3])];
    let refs = [SpeciesPotential {
        label: "H",
        potential: &pp,
    }];
    let base = input()
        .replace("[0.27,0.23,0.19]", "[0.5,0.5,0.5]")
        .replace("ecutwfc_ry=3.0", "ecutwfc_ry=12.0");
    let config = dft_cli::config::Config::parse(&base).unwrap();
    let mut opts = dft_cli::run::options(&config).unwrap();
    opts.initial_eigen_tolerance = opts.eigen_tolerance;
    opts.residual_eigen_schedule = false;
    let cell = |a: f64| Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap();
    // Broad scalar scan, independent of derivatives and the optimizer.
    let coarse: Vec<_> = (0..17)
        .map(|i| {
            let a = 3. + i as f64 * 0.25;
            let e = ScfProblem::new(&cell(a), &atoms, &refs, opts.clone())
                .unwrap()
                .run()
                .unwrap()
                .total_hartree();
            println!("coarse a={a} Bohr E={e} Ha");
            (a, e)
        })
        .collect();
    let index = (0..coarse.len())
        .min_by(|&i, &j| coarse[i].1.total_cmp(&coarse[j].1))
        .unwrap();
    assert!(
        index > 0 && index + 1 < coarse.len(),
        "minimum not bracketed"
    );
    let mut pressure_minima = Vec::new();
    for p in [0., 0.001] {
        let k = (1..coarse.len() - 1)
            .min_by(|&i, &j| {
                (coarse[i].1 + p * coarse[i].0.powi(3))
                    .total_cmp(&(coarse[j].1 + p * coarse[j].0.powi(3)))
            })
            .unwrap();
        let hs: Vec<_> = coarse[k - 1..=k + 1]
            .iter()
            .map(|(a, e)| e + p * a.powi(3))
            .collect();
        let center = coarse[k].0 + 0.125 * (hs[0] - hs[2]) / (hs[0] - 2. * hs[1] + hs[2]);
        for factor in [0.98, 1.02] {
            let start = center * factor;
            let text = base.replace(
                "[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]",
                &format!("[[{start},0.0,0.0],[0.0,{start},0.0],[0.0,0.0,{start}]]"),
            ) + &format!(
                "\n[relax]\nmode='cell'\ncell_dofs='isotropic'\nmax_steps=150\nstress_tol_ha_bohr3=1e-10\nenergy_tol_ha=1e-12\npressure_ha_bohr3={p}\n[symmetry]\nmode='auto'\n"
            );
            let name = format!("p{p}-f{factor}");
            fs::write(dir.join(format!("{name}.toml")), text).unwrap();
            let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
                .current_dir(&dir)
                .args(["relax", &format!("{name}.toml"), "--output-dir", &name])
                .output()
                .unwrap();
            assert!(
                output.status.success(),
                "{}\n{}",
                String::from_utf8_lossy(&output.stdout),
                String::from_utf8_lossy(&output.stderr)
            );
            let result: serde_json::Value =
                serde_json::from_slice(&fs::read(dir.join(&name).join("results.json")).unwrap())
                    .unwrap();
            let a = result["relaxation"]["frozen_final"]["frame"]["cell_columns_bohr"][0][0]
                .as_f64()
                .unwrap();
            let reference = cell(start);
            let grid = ScfProblem::new(&reference, &atoms, &refs, opts.clone())
                .unwrap()
                .grid_dimensions();
            let mut frozen = opts.clone();
            frozen.grid_override = Some(grid);
            frozen.strain_reference = Some(Box::new(reference));
            let energy = |a: f64| {
                ScfProblem::new(&cell(a), &atoms, &refs, frozen.clone())
                    .unwrap()
                    .run()
                    .unwrap()
                    .total_hartree()
                    + p * a.powi(3)
            };
            let h = 1e-3;
            let samples: Vec<_> = (-4..=4)
                .map(|i| (i as f64 * h, energy(a + i as f64 * h)))
                .collect();
            // Least-squares parabola on nine equally spaced independent energies.
            let sx2 = samples.iter().map(|(x, _)| x * x).sum::<f64>();
            let sx4 = samples.iter().map(|(x, _)| x.powi(4)).sum::<f64>();
            let sy = samples.iter().map(|(_, y)| y).sum::<f64>();
            let sxy = samples.iter().map(|(x, y)| x * y).sum::<f64>();
            let sx2y = samples.iter().map(|(x, y)| x * x * y).sum::<f64>();
            let b = sxy / sx2;
            let c = (sx2y - sx2 * sy / 9.) / (sx4 - sx2 * sx2 / 9.);
            let minimum = a - b / (2. * c);
            let relative = (minimum - a).abs() / minimum;
            assert!(
                if factor < 1. { start < a } else { start > a },
                "starts must bracket relaxed minimum"
            );
            pressure_minima.push(a);
            println!(
                "P={p} Ha/Bohr^3 start={start} Bohr relaxed={a} Bohr scan={minimum} Bohr relative={relative:e}"
            );
            assert!(c > 0. && relative < 1e-5);
            let lattice = &result["relaxation"]["frozen_final"]["frame"]["cell_columns_bohr"];
            for i in 0..3 {
                for j in 0..3 {
                    let x = lattice[i][j].as_f64().unwrap();
                    assert!((x - if i == j { a } else { 0. }).abs() < 1e-10);
                }
            }
        }
    }
    assert!(pressure_minima[2] < pressure_minima[0] && pressure_minima[3] < pressure_minima[1]);
}

#[test]
#[ignore = "combined strain/ionic SCF finite differences; release only"]
fn combined_scf_enthalpy_derivative() {
    use dft_engine::{ElectronicState, ForceState, ScfProblem, SpeciesPotential};
    use dft_geometry::{Atom, Cell};
    let xml = common::gaussian_with(1., 1.5, false, 801);
    let pp = dft_upf::parse_upf(&xml).unwrap();
    let refs = [SpeciesPotential {
        label: "H",
        potential: &pp,
    }];
    let cfg = dft_cli::config::Config::parse(&input()).unwrap();
    let a = cfg.cell().unwrap();
    let atoms = [Atom::new("H", [0.; 3]), Atom::new("H", [0.27, 0.23, 0.19])];
    let mut opts = dft_cli::run::options(&cfg).unwrap();
    opts.grid_override = Some(
        ScfProblem::new(&a, &atoms, &refs, opts.clone())
            .unwrap()
            .grid_dimensions(),
    );
    opts.strain_reference = Some(Box::new(a.clone()));
    let r = ScfProblem::new(&a, &atoms, &refs, opts.clone())
        .unwrap()
        .run()
        .unwrap();
    let view = ForceState {
        cell: &a,
        atoms: &atoms,
        potentials: &refs,
        options: &opts,
        electronic: ElectronicState::Fixed(&r),
        density_floor: None,
        symmetry: &[],
    };
    let forces = dft_engine::forces(&view).unwrap();
    let stress = dft_engine::stress::frozen_stress(&view).unwrap();
    let direction = [
        [0.17, 0.13, -0.07],
        [0.13, -0.09, 0.11],
        [-0.07, 0.11, 0.21],
    ];
    let ds = [[0.09, -0.12, 0.06], [-0.04, 0.03, 0.08]];
    let pressure = 0.0003;
    let mut derivative = 0.;
    for (i, row) in direction.iter().enumerate() {
        for (j, value) in row.iter().enumerate() {
            derivative += a.volume_bohr3()
                * (stress.values[i][j] + if i == j { pressure } else { 0. })
                * value;
        }
    }
    for (f, s) in forces.values.iter().zip(ds) {
        let dr = a.fractional_to_cartesian_bohr(s);
        derivative -= f.iter().zip(dr).map(|(f, x)| f * x).sum::<f64>();
    }
    let energy = |h: f64| {
        let matrix = a.lattice_columns_bohr();
        let b = Cell::from_matrix_columns_bohr(std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                matrix[i][j] + h * (0..3).map(|k| direction[i][k] * matrix[k][j]).sum::<f64>()
            })
        }))
        .unwrap();
        let atoms: Vec<_> = atoms
            .iter()
            .zip(ds)
            .map(|(a, s)| {
                Atom::new(
                    "H",
                    std::array::from_fn(|i| a.fractional_position[i] + h * s[i]),
                )
            })
            .collect();
        ScfProblem::new(&b, &atoms, &refs, opts.clone())
            .unwrap()
            .run()
            .unwrap()
            .total_hartree()
            + pressure * b.volume_bohr3()
    };
    let h = 5e-5;
    let fd = (energy(h) - energy(-h)) / (2. * h);
    let relative = (fd - derivative).abs() / derivative.abs();
    println!("combined SCF dH/dt analytic={derivative:e} Ha FD={fd:e} Ha relative={relative:e}");
    assert!(relative < 1e-7);
}

#[test]
#[ignore = "shape relaxation SCFs; release only"]
fn tetragonal_shape_returns_to_cubic() {
    let dir = std::env::temp_dir().join(format!("p23-shape-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 801).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    let text = input()
        .replace("[0.27,0.23,0.19]", "[0.5,0.5,0.5]")
        .replace("ecutwfc_ry=3.0", "ecutwfc_ry=12.0")
        .replace(
            "[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]",
            "[[3.85,0.0,0.0],[0.0,3.85,0.0],[0.0,0.0,3.95]]",
        )
        + "\n[relax]\nmode='cell'\ncell_dofs='full'\nmax_steps=150\nstress_tol_ha_bohr3=1e-9\nenergy_tol_ha=1e-12\n[symmetry]\nmode='auto'\n";
    fs::write(dir.join("input.toml"), text).unwrap();
    let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
        .current_dir(&dir)
        .args(["relax", "input.toml", "--output-dir", "out"])
        .output()
        .unwrap();
    assert!(
        output.status.success(),
        "{}\n{}",
        String::from_utf8_lossy(&output.stdout),
        String::from_utf8_lossy(&output.stderr)
    );
    let result: serde_json::Value =
        serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
    let a = &result["relaxation"]["frozen_final"]["frame"]["cell_columns_bohr"];
    let x = a[0][0].as_f64().unwrap();
    let z = a[2][2].as_f64().unwrap();
    println!(
        "tetragonal relaxed x={x:e} Bohr z={z:e} Bohr relative={:e}",
        (z - x).abs() / x
    );
    assert!((z - x).abs() / x < 1e-5);
}

#[test]
#[ignore = "external PseudoDojo Si 30/40 and MgO 50/60 Ha shifted-mesh trajectories; release only, set DFT_NC_DIR"]
fn real_files_cutoff_report() {
    let nc = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let dir = std::env::var_os("DFT_CELL_RUN_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            std::env::temp_dir().join(format!("p23-materials-{}", std::process::id()))
        });
    fs::create_dir_all(&dir).unwrap();
    println!("real-file artifacts: {}", dir.display());
    let templates =
        std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/reference/materials");
    let mut records = Vec::new();
    let mut failures = Vec::new();
    for (material, cutoffs) in [("Si", [30, 40]), ("MgO", [50, 60])] {
        for cutoff in cutoffs {
            for factor in [0.98, 1.02] {
                let name = format!("{material}-{cutoff}-{factor}");
                if std::env::var("DFT_CELL_CASE").is_ok_and(|v| v != name) {
                    continue;
                }
                let template = templates.join(format!("{}-{cutoff}.toml", material.to_lowercase()));
                let mut input: toml::Value =
                    toml::from_str(&fs::read_to_string(&template).unwrap()).unwrap();
                let nominal = 2. * input["system"]["lattice_vectors"][0][1].as_float().unwrap();
                for row in input["system"]["lattice_vectors"].as_array_mut().unwrap() {
                    for x in row.as_array_mut().unwrap() {
                        let value = x
                            .as_float()
                            .unwrap_or_else(|| x.as_integer().unwrap() as f64);
                        *x = toml::Value::Float(value * factor);
                    }
                }
                for (_, path) in input["pseudopotentials"].as_table_mut().unwrap().iter_mut() {
                    *path = toml::Value::String(
                        nc.join(path.as_str().unwrap()).to_str().unwrap().into(),
                    );
                }
                input["scf"]["max_iterations"] = toml::Value::Integer(300);
                input["scf"]["energy_tolerance_ha"] = toml::Value::Float(1e-11);
                input["scf"]["density_residual_tolerance"] = toml::Value::Float(1e-9);
                input["scf"]
                    .as_table_mut()
                    .unwrap()
                    .insert("eigen_tolerance".into(), toml::Value::Float(1e-10));
                let text = toml::to_string(&input).unwrap()
                    + "\n[relax]\nmode='cell'\ncell_dofs='isotropic'\nmax_steps=100\nstress_tol_ha_bohr3=1e-7\nenergy_tol_ha=1e-9\n[symmetry]\nmode='auto'\n";
                // Validate generated input before launching an expensive calculation.
                dft_cli::config::Config::parse(&text).unwrap();
                fs::write(dir.join(format!("{name}.toml")), text).unwrap();
                println!(
                    "START {name}, nominal={nominal} Bohr, start={} Bohr",
                    nominal * factor
                );
                let started = std::time::Instant::now();
                let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
                    .current_dir(&dir)
                    .args(["relax", &format!("{name}.toml"), "--output-dir", &name])
                    .output()
                    .unwrap();
                let elapsed = started.elapsed().as_secs_f64();
                fs::write(dir.join(format!("{name}.stdout")), &output.stdout).unwrap();
                fs::write(dir.join(format!("{name}.stderr")), &output.stderr).unwrap();
                let mut record = serde_json::json!({"case":name,"material":material,"cutoff_ha":cutoff,"start_factor":factor,"nominal_lattice_bohr":nominal,"start_lattice_bohr":nominal*factor,"elapsed_s":elapsed,"exit_code":output.status.code(),"input":input});
                if let Ok(bytes) = fs::read(dir.join(&name).join("results.json")) {
                    let result: serde_json::Value = serde_json::from_slice(&bytes).unwrap();
                    record["provenance"] = result["provenance"].clone();
                    let r = &result["relaxation"];
                    record["status"] = r["status"].clone();
                    record["failures"] = r["failures"].clone();
                    if let Some(half) =
                        r["frozen_final"]["frame"]["cell_columns_bohr"][0][1].as_f64()
                    {
                        let lattice = 2. * half;
                        record["lattice_bohr"] = serde_json::json!(lattice);
                        for key in [
                            "steps",
                            "frozen_final",
                            "rebuilt_final",
                            "basis_rebuild_difference",
                            "space_group_operations",
                        ] {
                            record[key] = r[key].clone();
                        }
                        println!(
                            "RESULT {name}: a={lattice:.12} Bohr frozen_pressure={} rebuilt_pressure={} Ha/Bohr^3 elapsed={elapsed:.3} s",
                            r["frozen_final"]["pressure_ha_bohr3"],
                            r["rebuilt_final"]["pressure_ha_bohr3"]
                        );
                        let a = &r["frozen_final"]["frame"]["cell_columns_bohr"];
                        for i in 0..3 {
                            for j in 0..3 {
                                if (a[i][j].as_f64().unwrap() - if i == j { 0. } else { half })
                                    .abs()
                                    > 1e-8
                                {
                                    failures.push(format!("{name}: cubic symmetry lost"));
                                }
                            }
                        }
                    }
                }
                if !output.status.success() {
                    println!("FAILED {name}: {}", String::from_utf8_lossy(&output.stderr));
                    failures.push(format!("{name}: exit {:?}", output.status.code()));
                }
                records.push(record);
                fs::write(
                    dir.join("summary.json"),
                    serde_json::to_vec_pretty(&records).unwrap(),
                )
                .unwrap();
            }
        }
    }
    assert!(
        failures.is_empty(),
        "{}; artifacts {}",
        failures.join("; "),
        dir.display()
    );
}

#[test]
#[ignore = "FIRE/BFGS cell constraints and periodic basis rebuild; release only"]
fn periodic_rebuild_fixed_axes_and_fire() {
    let dir = std::env::temp_dir().join(format!("p23-periodic-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 801).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    for (method, dofs) in [("bfgs", "full"), ("fire", "xyz_diagonal")] {
        let text = input()
            + &format!(
                "\n[relax]\nmode='cell'\nmethod='{method}'\ncell_dofs='{dofs}'\nfixed_axes=[false,false,true]\nrebuild_every=1\nmax_steps=2\n"
            );
        fs::write(dir.join(format!("{method}.toml")), text).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .current_dir(&dir)
            .args(["relax", &format!("{method}.toml"), "--output-dir", method])
            .output()
            .unwrap();
        assert_eq!(
            output.status.code(),
            Some(3),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let result: serde_json::Value =
            serde_json::from_slice(&fs::read(dir.join(method).join("results.json")).unwrap())
                .unwrap();
        assert_eq!(result["relaxation"]["steps"], 2);
        let frames = result["relaxation"]["frames"].as_array().unwrap();
        assert!(frames.iter().any(|f| f["provenance"] == "basis_rebuild"));
        for frame in frames {
            assert_eq!(frame["cell_columns_bohr"][2][2], 8.);
            assert_eq!(frame["fractional"], frames[0]["fractional"]);
        }
        // Independent CIF specification checks: lengths/angles and explicit sites.
        let cif = fs::read_to_string(dir.join(method).join("final.cif")).unwrap();
        let c = cif
            .lines()
            .find(|l| l.starts_with("_cell_length_c "))
            .unwrap()
            .split_whitespace()
            .nth(1)
            .unwrap()
            .parse::<f64>()
            .unwrap();
        assert!((c - 8. * dft_units::Angstrom::from(dft_units::Bohr(1.)).0).abs() < 1e-12);
        assert_eq!(
            cif.lines()
                .filter(|l| l.starts_with("H1 ") || l.starts_with("H2 "))
                .count(),
            2
        );
        assert!(cif.contains("_atom_site_fract_x"));
    }
}

#[test]
#[ignore = "external MgO 50 Ha and Si 30 Ha N=3/5 rebuild plus 2-3 polish cycles; release only, set DFT_NC_DIR"]
fn real_files_rebuild_start_agreement() {
    let nc = std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR");
    let work = std::env::var_os("DFT_CELL_REBUILD_WORK_DIR")
        .map(std::path::PathBuf::from)
        .unwrap_or_else(|| {
            std::env::temp_dir().join(format!("p23-rebuild-{}", std::process::id()))
        });
    let script = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/reference/relaxation/cell/rebuild_protocol.py");
    // The runner asserts the REBUILT pressure and agreement after polish.
    // Frozen-only start dependence remains a reported finite-basis observation.
    let status = Command::new("python3")
        .arg("-I")
        .arg(script)
        .args(["--binary", env!("CARGO_BIN_EXE_dft-cli")])
        .arg("--nc-dir")
        .arg(nc)
        .arg("--work-dir")
        .arg(&work)
        .status()
        .unwrap();
    assert!(
        status.success(),
        "rebuild comparison inconclusive; see {}",
        work.display()
    );
}
