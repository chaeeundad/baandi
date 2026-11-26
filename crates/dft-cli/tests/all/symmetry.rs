use config::{Config, SymmetryMode};
use dft_cli::config;
const INPUT: &str = include_str!("../../../../examples/si_lda.toml");
#[test]
fn defaults_mesh_weights_and_invalid_models() {
    let off = Config::parse(INPUT).unwrap();
    assert_eq!(off.symmetry.mode, SymmetryMode::Off);
    assert!(off.symmetry_map().unwrap().is_none());
    let base = INPUT.replace(
        "kpoints = [[0.0, 0.0, 0.0]]",
        "kgrid = [4,4,4]\nkshift = [0.5,0.5,0.5]",
    );
    let auto = Config::parse(&(base.clone() + "\n[symmetry]\nmode='auto'\n")).unwrap();
    let map = auto.symmetry_map().unwrap().unwrap();
    assert_eq!(map.points.len(), 10);
    assert_eq!(map.entries.len(), 64);
    assert_eq!(auto.prepare().unwrap().grid.dims(), [32; 3]);
    for block in [
        "[external_field]\nactive=true",
        "[constraints]\nactive=true",
        "[hubbard]\nactive=true",
        "[symmetry]\nmode='invalid'",
        "[symmetry]\nmode='auto'\ntolerance_bohr=-1.0",
    ] {
        assert!(Config::parse(&(base.clone() + "\n" + block)).is_err());
    }
    let invalid = base.replace(
        "ecutwfc_ry = 30.0",
        "ecutwfc_ry = 30.0\nfft_grid=[25,25,25]",
    ) + "\n[symmetry]\nmode='auto'\n";
    assert!(Config::parse(&invalid).unwrap().prepare().is_err());
}

#[test]
fn preparation_is_shared_and_invalidated_by_physical_inputs() {
    use std::sync::Arc;
    let input = INPUT.replace(
        "kpoints = [[0.0, 0.0, 0.0]]",
        "kgrid = [4,4,4]\nkshift = [0.5,0.5,0.5]",
    ) + "\n[symmetry]\nmode='auto'\n";
    let mut config = Config::parse(&input).unwrap();
    let first = config.prepared_symmetry().unwrap().unwrap();
    config.validate().unwrap();
    config.prepare().unwrap();
    config.resolved_points().unwrap();
    config.symmetry_json().unwrap();
    config.density_operations().unwrap();
    assert!(Arc::ptr_eq(
        &first,
        &config.prepared_symmetry().unwrap().unwrap()
    ));
    // Selecting the crystal density group reuses detection and k reduction.
    config.symmetry.density_group = config::DensityGroup::Crystal;
    assert_eq!(config.density_operations().unwrap().len(), 48);
    assert!(Arc::ptr_eq(
        &first,
        &config.prepared_symmetry().unwrap().unwrap()
    ));
    // After the parse-time projection, detection runs at min(tolerance, 1e-9)
    // Bohr (ADR 0049 amendment GS1B): only a tolerance below 1e-9 changes it.
    config.symmetry.tolerance_bohr = 5e-10;
    let tolerance = config.prepared_symmetry().unwrap().unwrap();
    assert!(!Arc::ptr_eq(&first, &tolerance));
    config.basis.kgrid = Some([2, 2, 2]);
    let mesh = config.prepared_symmetry().unwrap().unwrap();
    assert!(!Arc::ptr_eq(&tolerance, &mesh));
    assert_eq!(mesh.map().entries.len(), 8);
    config.basis.ktime_reversal = Some(false);
    let reversal = config.prepared_symmetry().unwrap().unwrap();
    assert!(!Arc::ptr_eq(&mesh, &reversal));
    config.system.atoms[1].fractional[0] += 0.01;
    let atoms = config.prepared_symmetry().unwrap().unwrap();
    assert!(!Arc::ptr_eq(&mesh, &atoms));
    config.system.lattice_vectors[0][1] += 0.01;
    let lattice = config.prepared_symmetry().unwrap().unwrap();
    assert!(!Arc::ptr_eq(&atoms, &lattice));
    config.symmetry.density_group = config::DensityGroup::Mesh;
    config.symmetry.mode = SymmetryMode::Off;
    assert!(config.prepared_symmetry().unwrap().is_none());
}

#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;

/// ADR 0049 end to end: an hcp cell with six-decimal fractional coordinates.
/// The default tolerance finds the full 24-operation group, the input is
/// projected onto exact sites, and the energy agrees with a symmetry-off run
/// on exact coordinates within the strict class (1e-10 Ha).
#[test]
fn six_decimal_hcp_full_group_matches_symmetry_off() {
    use std::{fs, process::Command};
    let dir = std::env::temp_dir().join(format!("gs1-sym-tol-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("syn.xml"),
        common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    let (a, c) = (4.6_f64, 7.4_f64); // Bohr; lattice written at full precision.
    let h = a * 3f64.sqrt() / 2.;
    let input = |p: [&str; 2], symmetry: &str| {
        format!(
            "[system]\nname='hcp'\nlattice_unit='bohr'\nlattice_vectors=[[{:?},{:?},0.0],[{:?},{:?},0.0],[0.0,0.0,{c:?}]]\n[[system.atoms]]\nspecies='H'\nfractional=[{},0.25]\n[[system.atoms]]\nspecies='H'\nfractional=[{},0.75]\n[pseudopotentials]\nH='syn.xml'\n[basis]\necutwfc_ry=4.0\nkgrid=[3,3,2]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=4\nmax_iterations=200\n[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\n{symmetry}",
            a / 2.,
            -h,
            a / 2.,
            h,
            p[0],
            p[1]
        )
    };
    let six = ["0.333333,0.666667", "0.666667,0.333333"];
    let exact = [
        format!("{:?},{:?}", 1. / 3., 2. / 3.),
        format!("{:?},{:?}", 2. / 3., 1. / 3.),
    ];
    let exact = [exact[0].as_str(), exact[1].as_str()];
    let parsed = Config::parse(&input(six, "[symmetry]\nmode='auto'\n")).unwrap();
    let record = parsed.symmetrization.as_ref().unwrap();
    assert_eq!(record.operations_count, 24);
    assert!(record.positions_changed && !record.lattice_changed);
    assert!(record.max_atom_displacement_bohr < 1e-5);
    let json = parsed.symmetry_json().unwrap().unwrap();
    assert_eq!(json["crystal_operations_count"], 24);
    assert_eq!(
        json["tolerance_bohr"],
        dft_geometry::symmetry::DEFAULT_TOLERANCE_BOHR
    );
    // The old 1e-7 Bohr tolerance keeps the six-decimal input at a lower group.
    let old = Config::parse(&input(
        six,
        "[symmetry]\nmode='auto'\ntolerance_bohr=1e-7\n",
    ))
    .unwrap();
    assert!(old.symmetrization.as_ref().unwrap().operations_count < 24);
    assert!(
        old.symmetry_map().unwrap().unwrap().points.len()
            > parsed.symmetry_map().unwrap().unwrap().points.len()
    );
    let mut energies = Vec::new();
    for (name, text) in [
        ("auto", input(six, "[symmetry]\nmode='auto'\n")),
        ("off", input(exact, "")),
    ] {
        fs::write(dir.join(format!("{name}.toml")), text).unwrap();
        let output = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .env("BAANDI_REFERENCE_STRICT", "on")
            .current_dir(&dir)
            .args(["run", &format!("{name}.toml"), "--output-dir", name])
            .output()
            .unwrap();
        assert!(
            output.status.success(),
            "{}",
            String::from_utf8_lossy(&output.stderr)
        );
        let result: serde_json::Value =
            serde_json::from_slice(&fs::read(dir.join(name).join("results.json")).unwrap())
                .unwrap();
        if name == "auto" {
            assert_eq!(result["symmetry"]["crystal_operations_count"], 24);
            assert_eq!(
                result["symmetry"]["input_symmetrization"]["positions_changed"],
                true
            );
        }
        energies.push(result["energies"]["total"]["ha"].as_f64().unwrap());
    }
    println!(
        "hcp six-decimal auto vs exact off: {:.17e} {:.17e} diff {:.3e} Ha",
        energies[0],
        energies[1],
        energies[0] - energies[1]
    );
    assert!((energies[0] - energies[1]).abs() <= 1e-10);
    fs::remove_dir_all(dir).unwrap();
}
