use dft_cli::{config::Config, run};

const DIAMOND: &str = include_str!("../../../examples/si_lda.toml");

#[test]
fn symmetry_grid_explicit_and_automatic_both_drivers() {
    for driver in ["legacy", "reference"] {
        for group in ["mesh", "crystal"] {
            for mesh in [
                "kpoints = [[0.0, 0.0, 0.0]]",
                "kgrid=[4,4,4]\nkshift=[0.5,0.5,0.5]",
            ] {
                let input = format!(
                    "{}\n[symmetry]\nmode='auto'\ndensity_group='{group}'\n[scf]\ndriver='{driver}'\n",
                    DIAMOND.replace("kpoints = [[0.0, 0.0, 0.0]]", mesh)
                );
                let auto = Config::parse(&input).unwrap();
                let grid = auto.prepare().unwrap().grid;
                assert_eq!(grid.dims(), [32; 3]);
                assert!(
                    !grid.overridden(),
                    "automatic symmetry padding is not a user override"
                );
                assert_eq!(run::options(&auto).unwrap().grid_override, Some([32; 3]));
                for n in [28, 32] {
                    let explicit = Config::parse(
                        &input.replace("[basis]", &format!("[basis]\nfft_grid=[{n},{n},{n}]")),
                    )
                    .unwrap();
                    let grid = explicit.prepare().unwrap().grid;
                    assert_eq!(grid.dims(), [n; 3]);
                    assert!(grid.overridden());
                    assert_eq!(run::options(&explicit).unwrap().grid_override, Some([n; 3]));
                }
                let bad = Config::parse(&input.replace("[basis]", "[basis]\nfft_grid=[25,25,25]"))
                    .unwrap();
                assert!(bad.prepare().unwrap_err().to_string().contains(
                    "symmetry.fft_grid: explicit grid does not preserve translations/rotations"
                ));
            }
        }
    }
}

#[allow(dead_code)]
mod common;

/// Actual output boundary on synthetic diamond, including recorded provenance.
#[test]
fn symmetry_grid_run_records_both_drivers() {
    use serde_json::Value;
    use std::{fs, process::Command};
    let dir = std::env::temp_dir().join(format!("g6b-grid-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    let mut pp =
        common::gaussian_with(2., 1.8, false, 161).replace("element=\"Xx\"", "element=\"Si\"");
    let parsed = dft_upf::parse_upf(&pp).unwrap();
    let array = |f: &dyn Fn(f64) -> f64| {
        parsed
            .mesh()
            .r()
            .iter()
            .map(|&r| format!("{:.17e}", f(r)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let start = pp.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
    let end = pp.find("</PP_RHOATOM>").unwrap();
    pp.replace_range(
        start..end,
        &array(&|r| {
            8. / (std::f64::consts::PI.sqrt() * 1.8_f64.powi(3))
                * r
                * r
                * (-(r / 1.8).powi(2)).exp()
        }),
    );
    pp = pp.replace("number_of_wfc=\"0\"", "number_of_wfc=\"1\"").replace("</UPF>", &format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1></PP_PSWFC></UPF>", array(&|r| r * (-0.6 * r * r).exp())));
    fs::write(dir.join("synthetic.xml"), pp).unwrap();
    for driver in ["legacy", "reference"] {
        for explicit in [false, true] {
            let base = DIAMOND
                .replace("/absolute/path/to/Si.UPF", "synthetic.xml")
                .replace("30.0", "2.0")
                .replace(
                    "kpoints = [[0.0, 0.0, 0.0]]",
                    "kgrid=[2,2,2]\nkshift=[0.5,0.5,0.5]",
                );
            let base = if explicit {
                base.replace("[basis]", "[basis]\nfft_grid=[28,28,28]")
            } else {
                base
            };
            let text = format!(
                "{base}\n[symmetry]\nmode='auto'\n[scf]\ndriver='{driver}'\nnbands=3\nmax_iterations=2\nverify='off'\n"
            );
            let config = Config::parse(&text).unwrap();
            let expected = config.prepare().unwrap().grid.dims();
            fs::write(dir.join("in.toml"), text).unwrap();
            let out = format!("{driver}-{explicit}");
            let result = Command::new(env!("CARGO_BIN_EXE_dft-cli"))
                .current_dir(&dir)
                .env("BAANDI_REFERENCE_SCF", "off")
                .env("BAANDI_REFERENCE_STRICT", "off")
                .env("BAANDI_CUTOFF_LADDER", "off")
                .args(["run", "in.toml", "--output-dir", &out])
                .output()
                .unwrap();
            assert!(
                matches!(result.status.code(), Some(0 | 3)),
                "{}",
                String::from_utf8_lossy(&result.stderr)
            );
            let json: Value =
                serde_json::from_slice(&fs::read(dir.join(&out).join("results.json")).unwrap())
                    .unwrap();
            assert_eq!(json["basis"]["fft_grid"], serde_json::json!(expected));
            assert_eq!(json["basis"]["fft_grid_overridden"], explicit);
            let density: Value =
                serde_json::from_slice(&fs::read(dir.join(&out).join("density.json")).unwrap())
                    .unwrap();
            assert_eq!(density["fft_grid"], serde_json::json!(expected));
            assert_eq!(
                fs::metadata(dir.join(&out).join("density.bin"))
                    .unwrap()
                    .len(),
                (expected.iter().product::<usize>() * 8) as u64
            );
            let console = String::from_utf8_lossy(&result.stderr);
            assert!(
                console.contains(if explicit {
                    "using explicit dimensions"
                } else {
                    "automatic dimensions"
                }),
                "{console}"
            );
            if driver == "reference" {
                let sym = &json["reference_scf"]["symmetry"];
                assert_eq!(sym["pre_projection_relative_units"], "dimensionless");
                for row in sym["maps"].as_array().unwrap() {
                    assert!(row["pre_projection_relative_norm"].as_f64().unwrap() > 1e-6);
                    assert!(row["post_projection_invariance_norm"].as_f64().unwrap() < 1e-13);
                    assert_eq!(row["residual_norm"], row["pre_projection_norm"]);
                }
            }
        }
    }
    fs::remove_dir_all(dir).unwrap();
}
