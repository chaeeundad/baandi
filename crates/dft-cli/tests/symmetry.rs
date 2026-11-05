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
