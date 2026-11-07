use config::{Config, SymmetryMode};
use dft_cli::config;
const INPUT: &str = include_str!("../../../examples/si_lda.toml");
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
    config.symmetry.tolerance_bohr *= 2.;
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
