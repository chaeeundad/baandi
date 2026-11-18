//! ONEPATH-2 (user decision 2026-11-16): stress (`output.stress`) and
//! variable-cell relaxation (`relax.mode = "cell"`, `"ions_cell"`) run on the
//! single standard path with the fixed-point verification (ADR 0043 item 5,
//! amendment ONEPATH-2). The strict path (`BAANDI_REFERENCE_STRICT=on`) is the
//! reference. Synthetic fixtures only (no real pseudopotential data).
//!
//! Bounds were fixed in `tickets/ONEPATH-2-STRESS.md` BEFORE any measurement
//! and are never widened:
//! - stress class: the synthetic stress oracle cap of `docs/validation/stress.md`
//!   (1e-7 Ha/Bohr^3 per component), applied to the frozen-Miller tensor
//!   difference standard against strict, to the pressure difference and to the
//!   isotropic Pulay correction difference, at 5e-9 Ha and at 5e-11 Ha
//!   (5e-9 failed on the skewed fixture; stress inputs now default to
//!   5e-11 Ha, so the default-threshold test runs at 5e-11 with the same bound;
//!   the 5e-9 numbers stay in `stress_error_decomposition_information`);
//! - energy: |E_standard - E_strict| <= the standard threshold;
//! - variable-cell relaxation: final lattice columns and Cartesian positions
//!   within 1e-3 Bohr per component of the strict reference (the ONEPATH-1
//!   position bound), every frame verified.
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use serde_json::Value;
use std::{fs, path::Path, process::Command};

/// Two-atom H fixture in a cubic 8 Bohr cell (local-only Gaussian ions).
const CUBIC: &str = "[system]\nname='onepath2 cubic'\nlattice_unit='bohr'\nlattice_vectors=[[8.0,0.0,0.0],[0.0,8.0,0.0],[0.0,0.0,8.0]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.27,0.23,0.19]\n[pseudopotentials]\nH='loc.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n";
/// Four-atom fixture (two species, six electrons) in a skewed cell; the He
/// ions carry s, s and p projectors (nonlocal Gaussian ions).
const SKEW: &str = "[system]\nname='onepath2 skew'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[1.4,6.6,0.0],[0.9,1.2,7.4]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='He'\nfractional=[0.31,0.22,0.17]\n[[system.atoms]]\nspecies='H'\nfractional=[0.62,0.55,0.48]\n[[system.atoms]]\nspecies='H'\nfractional=[0.71,0.83,0.66]\n[pseudopotentials]\nHe='nl.xml'\nH='loc.xml'\n[basis]\necutwfc_ry=3.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n";

/// Stress class (docs/validation/stress.md, synthetic oracle cap), Ha/Bohr^3.
const STRESS_CLASS: f64 = 1e-7;

fn setup(tag: &str) -> std::path::PathBuf {
    let dir = std::env::temp_dir().join(format!("onepath2-{tag}-{}", std::process::id()));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("loc.xml"),
        common::gaussian_with(1., 1.5, false, 401).replace("element=\"Xx\"", "element=\"H\""),
    )
    .unwrap();
    fs::write(
        dir.join("nl.xml"),
        common::gaussian_with(2., 1.2, true, 401).replace("element=\"Xx\"", "element=\"He\""),
    )
    .unwrap();
    dir
}

struct Out {
    code: Option<i32>,
    text: String,
    result: Option<Value>,
    wall: f64,
}

/// `command`: run or relax; `strict`: BAANDI_REFERENCE_STRICT=on.
fn cli(dir: &Path, command: &str, name: &str, input: &str, strict: bool) -> Out {
    fs::write(dir.join(format!("{name}.toml")), input).unwrap();
    let mut cmd = Command::new(env!("CARGO_BIN_EXE_dft-cli"));
    cmd.current_dir(dir);
    if command == "run" {
        cmd.args(["run", &format!("{name}.toml")]);
    } else {
        cmd.args([command, &format!("{name}.toml"), "--output-dir", name]);
    }
    if strict {
        cmd.env("BAANDI_REFERENCE_STRICT", "on");
    } else {
        cmd.env_remove("BAANDI_REFERENCE_STRICT");
    }
    cmd.env_remove("BAANDI_RELAX_EXTRAPOLATION");
    let started = std::time::Instant::now();
    let out = cmd.output().unwrap();
    Out {
        code: out.status.code(),
        text: format!(
            "{}\n{}",
            String::from_utf8_lossy(&out.stdout),
            String::from_utf8_lossy(&out.stderr)
        ),
        result: fs::read(dir.join(name).join("results.json"))
            .ok()
            .map(|b| serde_json::from_slice(&b).unwrap()),
        wall: started.elapsed().as_secs_f64(),
    }
}

fn ok(o: &Out) -> &Value {
    assert_eq!(o.code, Some(0), "{}", o.text);
    o.result.as_ref().unwrap()
}

fn flat(v: &Value) -> Vec<f64> {
    v.as_array()
        .unwrap()
        .iter()
        .flat_map(|row| {
            row.as_array()
                .unwrap()
                .iter()
                .map(|x| x.as_f64().unwrap())
                .collect::<Vec<_>>()
        })
        .collect()
}

fn max_abs_diff(a: &[f64], b: &[f64]) -> f64 {
    assert_eq!(a.len(), b.len());
    a.iter()
        .zip(b)
        .map(|(x, y)| (x - y).abs())
        .fold(0., f64::max)
}

fn schema_ok(path: &Path) -> bool {
    Command::new("python3")
        .arg("-I")
        .arg(Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/validate_stress_schema.py"))
        .arg(path)
        .status()
        .unwrap()
        .success()
}

/// Stress at 5e-11 Ha against the strict reference on both fixtures.
/// Bounds: module header (fixed before measuring). Measured values printed.
#[test]
#[ignore = "stress runs with three cutoff-scan SCFs each; run cargo test --release -p dft-cli --test all onepath_stress:: -- --ignored --nocapture"]
fn standard_stress_at_5e_11_agrees_with_strict_reference() {
    compare_stress("t5e-11", Some(5e-11));
}

/// Stress at the default threshold of a stress input against the strict
/// reference, bound unchanged (criterion a). History: with the class default
/// 5e-9 Ha this test FAILED on the skewed nonlocal fixture (1.305e-7 >
/// 1e-7 Ha/Bohr^3, relative 2.7e-5; cubic 2.1e-9, ledger ONEPATH-2). Since the
/// orchestrator decision of 2026-11-16 (adoption rule: the default must meet the
/// accuracy class) a stress input without `scf.standard_threshold_ha` runs at
/// `STRESS_DEFAULT_THRESHOLD_HA` = 5e-11 Ha (`threshold_source =
/// "stress_default"`); this test checks that rule and the bound (provisional
/// until the dev08 list confirms; margin 1.6x on the skewed fixture).
#[test]
#[ignore = "stress runs with three cutoff-scan SCFs each; run cargo test --release -p dft-cli --test all onepath_stress:: -- --ignored --nocapture"]
fn standard_stress_at_default_threshold_agrees_with_strict_reference() {
    compare_stress("tdefault", None);
}

/// Stress of the standard path with verification against the strict
/// reference on two fixtures (cubic local-only; skewed two-species with
/// s, s, p projectors), at the user threshold `user` or, when None, at the
/// stress default (5e-11 Ha). Every row is measured and printed before the
/// bounds are asserted.
fn compare_stress(tag: &str, user: Option<f64>) {
    let threshold = user.unwrap_or(5e-11);
    let dir = setup(tag);
    let mut violations = Vec::new();
    for (fixture, base) in [("cubic", CUBIC), ("skew", SKEW)] {
        let input = |name: &str, scf: &str| {
            format!("{base}[scf]\nnbands=4\n{scf}[output]\ndirectory='{name}'\nstress=true\n")
        };
        let strict_name = format!("{fixture}-strict");
        let strict = cli(&dir, "run", &strict_name, &input(&strict_name, ""), true);
        let s = ok(&strict);
        assert_eq!(s["single_path"]["status"], "reference_strict");
        assert_eq!(
            s["stress"]["preconditions"]["convergence_profile"],
            "strict"
        );
        assert!(schema_ok(&dir.join(&strict_name).join("results.json")));
        let sigma_strict = flat(&s["stress"]["values"]);
        let scale = sigma_strict.iter().map(|v| v.abs()).fold(0., f64::max);
        assert!(scale > 1e-6, "fixture has no stress to compare: {scale:e}");
        let p_strict = s["stress"]["pressure_ha_bohr3"].as_f64().unwrap();
        let pulay_strict = s["stress"]["pulay_assessment"]["isotropic_correction_ha_bohr3"]
            .as_f64()
            .unwrap();
        println!(
            "ONEPATH-2 stress {fixture} strict: max|sigma| {scale:.6e} Ha/Bohr^3, pressure {p_strict:.6e}, Pulay {pulay_strict:.6e}, wall {:.2} s",
            strict.wall
        );
        {
            let name = format!("{fixture}-{tag}");
            let scf = user.map_or(String::new(), |t| format!("standard_threshold_ha={t:e}\n"));
            let standard = cli(&dir, "run", &name, &input(&name, &scf), false);
            let v = ok(&standard);
            assert_eq!(v["single_path"]["status"], "applied");
            assert_eq!(
                v["single_path"]["threshold_source"],
                if user.is_some() {
                    "user"
                } else {
                    "stress_default"
                }
            );
            assert_eq!(v["single_path"]["standard_threshold_ha"], threshold);
            assert_eq!(v["single_path"]["verification"], "fixed_point");
            assert_eq!(v["verified"], true);
            let p = &v["stress"]["preconditions"];
            assert_eq!(p["convergence_profile"], "standard");
            assert_eq!(p["fixed_point_verification_passed"], true);
            assert_eq!(p["standard_threshold_ha"].as_f64().unwrap(), threshold);
            assert!(p["verification_coulomb_error_ha"].as_f64().unwrap() < threshold);
            let pulay = &v["stress"]["pulay_assessment"];
            assert_eq!(pulay["scan_profile"], "standard_verified");
            assert!(schema_ok(&dir.join(&name).join("results.json")));
            let diff = max_abs_diff(&flat(&v["stress"]["values"]), &sigma_strict);
            let dp = v["stress"]["pressure_ha_bohr3"].as_f64().unwrap() - p_strict;
            let dpulay = pulay["isotropic_correction_ha_bohr3"].as_f64().unwrap() - pulay_strict;
            let de = v["total_energy"]["ha"].as_f64().unwrap()
                - s["total_energy"]["ha"].as_f64().unwrap();
            println!(
                "ONEPATH-2 stress {fixture} threshold {threshold:e} Ha: max|sigma_std - sigma_strict| {diff:.3e} Ha/Bohr^3 (relative {:.3e}), dP {dp:.3e}, dPulay {dpulay:.3e} Ha/Bohr^3, dE {de:.3e} Ha, scans {}, verification attempts {}, coulomb {:.3e} Ha, wall {:.2} s",
                diff / scale,
                pulay["scan_thresholds_ha"],
                p["verification_attempts"],
                p["verification_coulomb_error_ha"].as_f64().unwrap(),
                standard.wall
            );
            for (what, value, bound) in [
                ("max|dsigma|", diff, STRESS_CLASS),
                ("|dP|", dp.abs(), STRESS_CLASS),
                ("|dPulay|", dpulay.abs(), STRESS_CLASS),
                ("|dE|", de.abs(), threshold),
            ] {
                if value > bound {
                    violations.push(format!(
                        "{fixture} {threshold:e}: {what} {value:e} > {bound:e}"
                    ));
                }
            }
        }
    }
    fs::remove_dir_all(dir).unwrap();
    assert!(violations.is_empty(), "{violations:?}");
}

/// Information only (no bound, nothing asserted beyond successful runs): the
/// error of each path against a tightened strict reference (strict path with
/// energy 1e-12 Ha, density residual 1e-10, eigen 1e-11 Ha), to separate the
/// standard path's own stress error from that of the default strict reference.
#[test]
#[ignore = "information-only decomposition; run cargo test --release -p dft-cli --test all onepath_stress::stress_error_decomposition -- --ignored --nocapture"]
fn stress_error_decomposition_information() {
    let dir = setup("decomp");
    for (fixture, base) in [("cubic", CUBIC), ("skew", SKEW)] {
        let input = |name: &str, scf: &str| {
            format!("{base}[scf]\nnbands=4\n{scf}[output]\ndirectory='{name}'\nstress=true\n")
        };
        let tight_name = format!("{fixture}-tight");
        let tight = cli(
            &dir,
            "run",
            &tight_name,
            &input(
                &tight_name,
                "max_iterations=300\nenergy_tolerance_ha=1e-12\ndensity_residual_tolerance=1e-10\neigen_tolerance=1e-11\n",
            ),
            true,
        );
        let t = ok(&tight);
        let sigma = flat(&t["stress"]["values"]);
        let e = t["total_energy"]["ha"].as_f64().unwrap();
        for (label, scf, strict) in [
            ("strict default", String::new(), true),
            (
                "standard 5e-9",
                "standard_threshold_ha=5e-9\n".into(),
                false,
            ),
            (
                "standard 5e-10",
                "standard_threshold_ha=5e-10\n".into(),
                false,
            ),
            (
                "standard 5e-11",
                "standard_threshold_ha=5e-11\n".into(),
                false,
            ),
            (
                "standard 5e-13",
                "standard_threshold_ha=5e-13\n".into(),
                false,
            ),
        ] {
            let name = format!("{fixture}-{}", label.replace(' ', "-"));
            let out = cli(&dir, "run", &name, &input(&name, &scf), strict);
            let v = ok(&out);
            println!(
                "ONEPATH-2 decomposition {fixture} {label}: max|sigma - sigma_tight| {:.3e} Ha/Bohr^3, dE {:.3e} Ha, SCF iterations {}",
                max_abs_diff(&flat(&v["stress"]["values"]), &sigma),
                v["total_energy"]["ha"].as_f64().unwrap() - e,
                v["convergence"]["iterations"]
            );
        }
    }
    fs::remove_dir_all(dir).unwrap();
}

/// Two H ions in a body-centred cubic arrangement, compressed cube (3.7 Bohr):
/// the isotropic cell relaxation expands it toward the stress-free volume.
const BCC: &str = "[system]\nname='onepath2 bcc'\nlattice_unit='bohr'\nlattice_vectors=[[3.7,0.0,0.0],[0.0,3.7,0.0],[0.0,0.0,3.7]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.5,0.5,0.5]\n[pseudopotentials]\nH='loc.xml'\n[basis]\necutwfc_ry=12.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n";
/// Two H ions in a compact skewed cell (the ONEPATH-1 relaxation fixture, its
/// stiff ionic minimum), with ions and full cell free.
const COMPACT: &str = "[system]\nname='onepath2 compact skew'\nlattice_unit='bohr'\nlattice_vectors=[[4.0,0.0,0.0],[0.4,4.0,0.0],[0.3,0.2,4.2]]\n[[system.atoms]]\nspecies='H'\nfractional=[0.0,0.0,0.0]\n[[system.atoms]]\nspecies='H'\nfractional=[0.4,0.42,0.37]\n[pseudopotentials]\nH='loc.xml'\n[basis]\necutwfc_ry=6.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='lda-pz'\n";

/// Short variable-cell relaxations (cell isotropic on BCC; ions_cell with the
/// full cell on COMPACT) on the standard path with verification reach the
/// same cell and structure as the strict reference. Bounds: module header.
#[test]
#[ignore = "four variable-cell relaxations; run cargo test --release -p dft-cli --test all onepath_stress:: -- --ignored --nocapture"]
fn standard_variable_cell_relaxation_matches_strict_reference() {
    let dir = setup("cell");
    for (fixture, base, relax) in [
        (
            "bcc-cell",
            BCC,
            "[relax]\nmode='cell'\ncell_dofs='isotropic'\nmax_steps=60\nstress_tol_ha_bohr3=1e-7\n",
        ),
        (
            "compact-ions-cell",
            COMPACT,
            "[relax]\nmode='ions_cell'\ncell_dofs='full'\nmax_steps=80\nforce_tol_ha_bohr=1e-4\nstress_tol_ha_bohr3=1e-7\n",
        ),
    ] {
        let input = format!("{base}[scf]\nnbands=4\n{relax}");
        let mut finals = Vec::new();
        for strict in [true, false] {
            let name = format!("{fixture}-{}", if strict { "strict" } else { "standard" });
            let out = cli(&dir, "relax", &name, &input, strict);
            let v = ok(&out);
            let r = &v["relaxation"];
            assert_eq!(r["status"], "converged", "{name}");
            let frames = r["frames"].as_array().unwrap();
            if strict {
                assert_eq!(v["single_path"]["status"], "reference_strict");
                assert!(frames.iter().all(|f| f.get("verification").is_none()));
            } else {
                assert_eq!(v["single_path"]["status"], "applied");
                assert_eq!(v["single_path"]["verification"], "fixed_point");
                // Variable-cell input: the stress default applies.
                assert_eq!(v["scf"]["standard_threshold_ha"], 5e-11);
                assert_eq!(v["scf"]["threshold_source"], "stress_default");
                for f in frames {
                    let ver = &f["verification"];
                    assert_eq!(ver["passed"], true, "{f}");
                    assert!(ver["coulomb_error_ha"].as_f64().unwrap() < 5e-11);
                }
            }
            let last = frames.last().unwrap();
            let cell = flat(&last["cell_columns_bohr"]);
            let x: Vec<f64> = last["positions"]
                .as_array()
                .unwrap()
                .iter()
                .map(|v| v.as_f64().unwrap())
                .collect();
            let scf: u64 = frames
                .iter()
                .map(|f| f["scf_iterations"].as_u64().unwrap())
                .sum();
            let h: u64 = frames
                .iter()
                .map(|f| f["scf_h_vectors"].as_u64().unwrap())
                .sum();
            let stress = flat(&last["stress_ha_bohr3"]);
            println!(
                "ONEPATH-2 relax {name}: frames {}, total SCF iterations {scf}, H vectors {h}, final cell {:?} Bohr, final max force {:.3e} Ha/Bohr, final max |sigma| {:.3e} Ha/Bohr^3, wall {:.2} s",
                frames.len(),
                cell,
                last["max_force_ha_bohr"].as_f64().unwrap(),
                stress.iter().map(|v| v.abs()).fold(0., f64::max),
                out.wall
            );
            finals.push((cell, x));
        }
        let dcell = max_abs_diff(&finals[0].0, &finals[1].0);
        let dx = max_abs_diff(&finals[0].1, &finals[1].1);
        println!(
            "ONEPATH-2 relax {fixture}: max |a_standard - a_strict| {dcell:.3e} Bohr, max |x_standard - x_strict| {dx:.3e} Bohr"
        );
        assert!(dcell < 1e-3, "{fixture}: cell {dcell:e} Bohr");
        assert!(dx < 1e-3, "{fixture}: positions {dx:e} Bohr");
    }
    fs::remove_dir_all(dir).unwrap();
}

/// ONEPATH-2B release check on a real pseudopotential (never committed):
/// `DFT_NC_DIR/Si.upf` (PseudoDojo Si, PBE), conventional Si8 at Gamma, 30 Ry
/// (the dev08 si8-gamma case, FFT 36^3). The standard arms failed at step 3
/// with "Miller index [0, 0, 18] is outside the canonical frequencies for
/// [36, 36, 36]" before the Coulomb metric used the frozen density list.
/// Checks: no Miller-index failure in any arm; the standard arm ends with the
/// strict reference's status; the mode-cell case (the dev08 case) converges;
/// where both converge, final cell and positions within 1e-3 Bohr (the
/// synthetic test's bound).
#[test]
#[ignore = "real PseudoDojo Si PBE; release only, set DFT_NC_DIR: cargo test --release -p dft-cli --test all onepath_stress::real_si8 -- --ignored --nocapture"]
fn real_si8_variable_cell_matches_strict_reference() {
    let nc = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let upf = nc.join("Si.upf");
    assert!(upf.is_file(), "missing {}", upf.display());
    let dir = setup("real-si8");
    let a = 10.261212856717933;
    let atoms = [
        [0., 0., 0.],
        [0., 0.5, 0.5],
        [0.5, 0., 0.5],
        [0.5, 0.5, 0.],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ];
    let input = |scale: f64, shift: f64, relax: &str| {
        let mut s = format!(
            "[system]\nname='si8 cell'\nlattice_unit='bohr'\nlattice_vectors=[[{x},0.0,0.0],[0.0,{x},0.0],[0.0,0.0,{x}]]\n",
            x = a * scale
        );
        for (i, f) in atoms.iter().enumerate() {
            let d = if i == 4 { shift } else { 0. };
            s += &format!(
                "[[system.atoms]]\nspecies='Si'\nfractional=[{},{},{}]\n",
                f[0] + d,
                f[1] - d,
                f[2] + 0.5 * d
            );
        }
        s + &format!(
            "[pseudopotentials]\nSi='{}'\n[basis]\necutwfc_ry=30.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional='pbe'\n{relax}",
            upf.display()
        )
    };
    for (case, text) in [
        (
            "cell",
            input(
                0.98,
                0.,
                "[relax]\nmode='cell'\ncell_dofs='isotropic'\nforce_tol_ha_bohr=1e-3\nmax_steps=30\n",
            ),
        ),
        (
            "ions-cell",
            // ions_cell, isotropic cell, one atom displaced. On this input the
            // STRICT reference itself does not converge (the trajectory reaches
            // structures where the SCF diverges near step 23; also with the full
            // cell, lattice x 1.01 and displacement 0.01: step 20 "SCF
            // unconverged"), and the standard arm follows the same trajectory
            // to the same failure (ledger ONEPATH-2B). That is a separate,
            // path-independent finding; here the check is that the standard arm
            // never fails with the Miller-index error and ends like the
            // reference.
            input(
                1.01,
                0.005,
                "[relax]\nmode='ions_cell'\ncell_dofs='isotropic'\nforce_tol_ha_bohr=1e-3\nmax_steps=60\n",
            ),
        ),
    ] {
        let mut finals = Vec::new();
        let mut statuses = Vec::new();
        for strict in [true, false] {
            let name = format!("{case}-{}", if strict { "strict" } else { "standard" });
            let out = cli(&dir, "relax", &name, &text, strict);
            let v = out
                .result
                .as_ref()
                .unwrap_or_else(|| panic!("{}", out.text));
            let r = &v["relaxation"];
            for f in r["failures"].as_array().into_iter().flatten() {
                assert!(
                    !f["reason"].as_str().unwrap_or("").contains("Miller index"),
                    "{name}: {f}"
                );
            }
            statuses.push(r["status"].clone());
            let frames = r["frames"].as_array().unwrap();
            let last = frames.last().unwrap();
            println!(
                "ONEPATH-2B real Si8 {name}: status {}, steps {}, SCF iterations {}, final a00 {} Bohr, max force {:.3e}, wall {:.1} s",
                r["status"],
                frames.len() - 1,
                frames
                    .iter()
                    .map(|f| f["scf_iterations"].as_u64().unwrap())
                    .sum::<u64>(),
                last["cell_columns_bohr"][0][0],
                last["max_force_ha_bohr"].as_f64().unwrap(),
                out.wall
            );
            finals.push((
                flat(&last["cell_columns_bohr"]),
                last["positions"]
                    .as_array()
                    .unwrap()
                    .iter()
                    .map(|v| v.as_f64().unwrap())
                    .collect::<Vec<_>>(),
            ));
        }
        // The standard arm ends like the reference; the mode-cell case is the
        // failing dev08 case and must converge.
        assert_eq!(statuses[1], statuses[0], "{case}");
        if case == "cell" {
            assert_eq!(statuses[0], "converged");
        }
        let dcell = max_abs_diff(&finals[0].0, &finals[1].0);
        let dx = max_abs_diff(&finals[0].1, &finals[1].1);
        println!(
            "ONEPATH-2B real Si8 {case}: status {}, max |a_std - a_strict| {dcell:.3e} Bohr, max |x_std - x_strict| {dx:.3e} Bohr",
            statuses[0]
        );
        if statuses[0] == "converged" {
            assert!(dcell < 1e-3, "{case}: cell {dcell:e} Bohr");
            assert!(dx < 1e-3, "{case}: positions {dx:e} Bohr");
        }
    }
    fs::remove_dir_all(dir).unwrap();
}
