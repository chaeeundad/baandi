#[test]
fn relax_without_table_rejects_reference_before_scf() {
    let dir = std::env::temp_dir().join(format!("reference-relax-{}", std::process::id()));
    std::fs::create_dir_all(&dir).unwrap();
    for profile in ["strict", "standard"] {
        let path = dir.join("input.toml");
        std::fs::write(
            &path,
            input(&format!("{MODE}[scf]\nconvergence_profile='{profile}'\n")),
        )
        .unwrap();
        let result = std::process::Command::new(env!("CARGO_BIN_EXE_dft-cli"))
            .arg("relax")
            .arg(&path)
            .arg("--output-dir")
            .arg(dir.join(profile))
            .output()
            .unwrap();
        assert_eq!(result.status.code(), Some(1));
        assert!(
            String::from_utf8_lossy(&result.stderr).contains(
                "compat.mode = \"reference\" is for unpolarized PBE energies and densities"
            )
        );
        assert!(!dir.join(profile).exists());
    }
    std::fs::remove_dir_all(dir).unwrap();
}

use dft_cli::{
    config::{Config, DensityGroup},
    run::options,
};
use dft_engine::CompatMode;

fn input(extra: &str) -> String {
    format!(
        "[system]\nname='reference test'\nlattice_unit='bohr'\nlattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]\n[[system.atoms]]\nspecies='Si'\nfractional=[0.0,0.0,0.0]\n[pseudopotentials]\nSi='unused.upf'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='pbe'\n{extra}"
    )
}
const MODE: &str = "\n[compat]\nmode='reference'\n";

#[test]
fn defaults_resolution_and_explicit_agreement() {
    let absent = Config::parse(&input("")).unwrap();
    let none = Config::parse(&input("[compat]\nmode='none'")).unwrap();
    for c in [&absent, &none] {
        let o = options(c).unwrap();
        assert_eq!(o.compat_mode, CompatMode::None);
        assert_eq!(o.pw92_prefactor.name(), "modified");
        assert_eq!(o.local_mesh.name(), "baandi");
        assert_eq!(o.radial_even_rule.name(), "simpson_3_8");
        assert!(!o.crystal_density_group);
    }
    for symmetry in ["", "[symmetry]\nmode='off'", "[symmetry]\nmode='auto'"] {
        let c = Config::parse(&input(&format!("{MODE}{symmetry}"))).unwrap();
        let o = options(&c).unwrap();
        assert_eq!(o.compat_mode, CompatMode::Reference);
        assert_eq!(o.pw92_prefactor.name(), "original");
        assert_eq!(o.local_mesh.name(), "qe_like");
        assert_eq!(o.radial_even_rule.name(), "simpson_quadratic_end");
        assert_eq!(
            c.symmetry.density_group == DensityGroup::Crystal,
            symmetry.contains("auto")
        );
    }
    Config::parse(&input(&format!("pw92_prefactor='original'\n{MODE}[pseudo]\nlocal_mesh='qe_like'\nradial_even_rule='simpson_quadratic_end'\n[symmetry]\nmode='auto'\ndensity_group='crystal'"))).unwrap();
}

#[test]
fn conflicts_name_the_mode_and_explicit_key() {
    for (extra, key) in [
        ("pw92_prefactor='modified'", "xc.pw92_prefactor"),
        ("[pseudo]\nlocal_mesh='baandi'", "pseudo.local_mesh"),
        (
            "[pseudo]\nradial_even_rule='simpson_3_8'",
            "pseudo.radial_even_rule",
        ),
        (
            "[symmetry]\nmode='auto'\ndensity_group='mesh'",
            "symmetry.density_group",
        ),
        (
            "[symmetry]\nmode='off'\ndensity_group='crystal'",
            "symmetry.density_group",
        ),
    ] {
        let err = Config::parse(&input(&format!("{extra}{MODE}")))
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("compat.mode = \"reference\"") && err.contains(key),
            "{err}"
        );
    }
    for extra in [
        "[compat]\nmode='other'",
        "[compat]\nmode='reference'\nunknown=true",
    ] {
        assert!(Config::parse(&input(extra)).is_err());
    }
}

#[test]
fn unsupported_scope_has_one_mode_message() {
    for extra in [
        "[spin]\nmode='collinear'",
        "[output]\nforces=true",
        "[output]\nstress=true",
        "[relax]",
    ] {
        let err = Config::parse(&input(&format!("{MODE}{extra}")))
            .unwrap_err()
            .to_string();
        assert!(
            err.contains(
                "compat.mode = \"reference\" is for unpolarized PBE energies and densities"
            ),
            "{err}"
        );
    }
    for lda in ["lda-pz", "lda-pw"] {
        let err = Config::parse(&input(MODE).replace("'pbe'", &format!("'{lda}'")))
            .unwrap_err()
            .to_string();
        assert!(
            err.contains(
                "compat.mode = \"reference\" is for unpolarized PBE energies and densities"
            ),
            "{err}"
        );
    }
}
