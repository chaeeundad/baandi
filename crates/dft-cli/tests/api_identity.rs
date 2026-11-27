//! R3 adapter identity. Compare against a pre-refactor executable built from
//! the same checkout/toolchain. No QE data or pseudopotentials are committed.
//! Release matrix (includes stress cutoff scans and ionic continuation):
//! BAANDI_API_BASELINE_BINARY=/path/to/pre-r3/dft-cli cargo test --release
//!   -p dft-cli --test api_identity -- --ignored --nocapture
//! Baseline policy (I-API-GAMMA-LASTBITS): build the baseline with the same
//! features and RUSTFLAGS as this test (a production binary built with
//! `-C target-cpu=native` or `--features mkl` is not a valid baseline for a
//! plain `cargo test --release` build), and from a commit whose results schema
//! matches (RELAX-STOP-RULE added relaxation fields, so relax cases need a
//! baseline at or after 0fd3ac96). Commit ids and switch help text are masked.
#[allow(dead_code)]
mod common;
use serde_json::{Value, json};
use std::{fs, path::Path, process::Command};

// Mask wall clocks only, retaining schemas, counts and other numeric fields in
// timing containers. RSS is separately normalized: it is process/allocator
// telemetry, not deterministic scientific or work-accounting evidence.
fn mask_wall_clocks(value: &mut Value, path: &str) {
    match value {
        Value::Object(map) => {
            for (key, value) in map {
                let child = format!("{path}/{key}");
                let clock = key.ends_with("_s")
                    || key.ends_with("_seconds")
                    || matches!(
                        key.as_str(),
                        "timestamp" | "elapsed_seconds" | "wall_seconds"
                    )
                    || path.starts_with("/timings")
                    || path.starts_with("/reference_scf/setup_s")
                    || path.starts_with("/reference_scf/output_s")
                    || path.starts_with("/reference_scf/wall_s")
                    || path.ends_with("/stage_s");
                if clock {
                    mask_numeric_leaves(value);
                } else {
                    mask_wall_clocks(value, &child);
                }
            }
        }
        Value::Array(items) => {
            for item in items {
                mask_wall_clocks(item, path);
            }
        }
        _ => {}
    }
}
fn mask_numeric_leaves(value: &mut Value) {
    match value {
        Value::Object(map) => {
            for value in map.values_mut() {
                mask_numeric_leaves(value);
            }
        }
        Value::Array(items) => {
            for value in items {
                mask_numeric_leaves(value);
            }
        }
        Value::Number(_) => *value = json!(0),
        Value::String(_) => {}
        _ => {}
    }
}
fn normalize(value: &mut Value) {
    mask_wall_clocks(value, "");
    // Explicit process provenance, not scientific/work fields. Different
    // executables necessarily have different paths/digests; their digests are
    // independently checked by execute() before normalization.
    for pointer in [
        "/provenance/executable_sha256",
        "/timestamp",
        "/run/run_id",
        "/run/source_commit",
        "/tool/git_commit",
        "/run/binary/path",
        "/run/binary/sha256",
        "/run/binary/hash_cache/status",
        "/run/binary/hash_cache/hashed_bytes",
        "/run/binary/hash_cache/published",
    ] {
        if let Some(field) = value.pointer_mut(pointer) {
            *field = Value::Null;
        }
    }
    // The switch catalog's help text (allowed choices, meaning) is
    // documentation that grows with diagnostic choices; the effective
    // value, source, default and request of every switch stay compared.
    if let Some(Value::Object(switches)) = value.pointer_mut("/scf/switches") {
        for entry in switches.values_mut() {
            for key in ["allowed", "meaning"] {
                if let Some(field) = entry.get_mut(key) {
                    *field = Value::Null;
                }
            }
        }
    }
    if let Some(rss) = value.pointer_mut("/reference_scf/run/peak_rss_bytes_per_rank") {
        mask_numeric_leaves(rss);
    }
}

fn first_difference(a: &Value, b: &Value, path: &str) -> Option<String> {
    match (a, b) {
        (Value::Object(a), Value::Object(b)) if a.keys().eq(b.keys()) => a
            .iter()
            .find_map(|(k, v)| first_difference(v, &b[k], &format!("{path}/{k}"))),
        (Value::Array(a), Value::Array(b)) if a.len() == b.len() => a
            .iter()
            .zip(b)
            .enumerate()
            .find_map(|(i, (a, b))| first_difference(a, b, &format!("{path}/{i}"))),
        _ if a == b => None,
        _ => Some(format!("{path}: {a} != {b}")),
    }
}
#[test]
fn mask_retains_science_work_and_verification() {
    let mut a = json!({"timings":{"total_s":1.},"energy_ha":-2.,
        "verification":{"passed":true,"wall_seconds":2.,"operator_applications":4},
        "reference_scf":{"iterations":[{"map_start_s":1.,"dr2_ry":0.01,"calls":[{"iterations":3}]}]}});
    let mut b = a.clone();
    b["timings"]["total_s"] = json!(7.);
    b["verification"]["wall_seconds"] = json!(9.);
    normalize(&mut a);
    normalize(&mut b);
    assert_eq!(a, b);
    for pointer in [
        "/energy_ha",
        "/verification/operator_applications",
        "/verification/passed",
        "/reference_scf/iterations/0/dr2_ry",
        "/reference_scf/iterations/0/calls/0/iterations",
    ] {
        let mut changed = b.clone();
        *changed.pointer_mut(pointer).unwrap() = json!(99);
        normalize(&mut changed);
        assert!(first_difference(&a, &changed, "").is_some(), "{pointer}");
    }
}
fn synthetic(root: &Path) -> String {
    let mut pp =
        common::gaussian_with(2., 1.8, true, 161).replace("element=\"Xx\"", "element=\"He\"");
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
    pp = pp.replace("number_of_wfc=\"0\"", "number_of_wfc=\"1\"").replace("</UPF>",
        &format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1></PP_PSWFC></UPF>", array(&|r| r * (-0.6 * r * r).exp())));
    fs::write(root.join("syn.xml"), pp).unwrap();
    "[system]\nname='R3 identity'\nlattice_unit='bohr'\nlattice_vectors=[[7.0,0.0,0.0],[0.2,7.0,0.0],[0.0,0.1,7.0]]\n[[system.atoms]]\nspecies='He'\nfractional=[0.13,0.21,0.32]\n[pseudopotentials]\nHe='syn.xml'\n[basis]\necutwfc_ry=2.0\nkgrid=[2,2,2]\nkshift=[0.5,0.5,0.5]\n[xc]\nfunctional='lda-pz'\n[scf]\nnbands=3\nmax_iterations=150\nstandard_threshold_ha=5e-13\ndriver='reference'\n".into()
}

fn execute(
    binary: &Path,
    command: &str,
    input: &Path,
    out: &Path,
) -> (i32, Value, Option<Vec<u8>>) {
    let output = Command::new(binary)
        .arg(command)
        .arg(input)
        .arg("--output-dir")
        .arg(out)
        .output()
        .unwrap();
    let code = output.status.code().unwrap();
    assert!(
        matches!(code, 0 | 3),
        "{command}: {code}: {}",
        String::from_utf8_lossy(&output.stderr)
    );
    let report = fs::read(out.join("results.json")).unwrap();
    let value: Value = serde_json::from_slice(&report).unwrap();
    for digest in [
        value.pointer("/run/binary/sha256"),
        value.pointer("/provenance/executable_sha256"),
    ]
    .into_iter()
    .flatten()
    {
        use sha2::{Digest, Sha256};
        assert_eq!(
            digest,
            &json!(format!("{:x}", Sha256::digest(fs::read(binary).unwrap())))
        );
    }
    (code, value, fs::read(out.join("density.bin")).ok())
}
/// Size of a density mismatch (density.bin is little-endian f64): differing
/// values, first index, largest absolute and relative difference.
fn density_difference(old: Option<&[u8]>, new: Option<&[u8]>) -> String {
    let (Some(old), Some(new)) = (old, new) else {
        return format!("present old={} new={}", old.is_some(), new.is_some());
    };
    if old.len() != new.len() || old.len() % 8 != 0 {
        return format!("lengths old={} new={}", old.len(), new.len());
    }
    let values = |b: &[u8]| -> Vec<f64> {
        b.as_chunks::<8>()
            .0
            .iter()
            .map(|c| f64::from_le_bytes(*c))
            .collect()
    };
    let (old, new) = (values(old), values(new));
    let (mut count, mut first, mut abs, mut rel) = (0, None, 0f64, 0f64);
    for (i, (a, b)) in old.iter().zip(&new).enumerate() {
        if a.to_bits() != b.to_bits() {
            count += 1;
            first.get_or_insert(i);
            abs = abs.max((a - b).abs());
            rel = rel.max((a - b).abs() / a.abs().max(b.abs()));
        }
    }
    format!(
        "{count} of {} values differ, first index {first:?}, max |diff| {abs:e}, max relative {rel:e}",
        old.len()
    )
}
#[test]
fn density_difference_reports_size() {
    let bytes = |v: &[f64]| v.iter().flat_map(|x| x.to_le_bytes()).collect::<Vec<_>>();
    let a = bytes(&[1., 2., 4.]);
    let b = bytes(&[1., 2. * (1. + f64::EPSILON), 4.]);
    let text = density_difference(Some(&a), Some(&b));
    assert!(
        text.starts_with("1 of 3 values differ, first index Some(1)"),
        "{text}"
    );
    assert_eq!(
        density_difference(Some(&a), None),
        "present old=true new=false"
    );
}
#[test]
#[ignore = "release pre-refactor binary matrix; see module header"]
fn cli_results_match_pre_facade_binary() {
    let baseline = std::env::var_os("BAANDI_API_BASELINE_BINARY")
        .expect("set BAANDI_API_BASELINE_BINARY to the pre-refactor executable");
    let root = std::env::temp_dir().join(format!("r3-api-identity-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let base = synthetic(&root).replace("driver='reference'\n", "threads=1\n");
    let gamma = base.replace(
        "kgrid=[2,2,2]\nkshift=[0.5,0.5,0.5]",
        "kpoints=[[0.0,0.0,0.0]]",
    );
    let symmetry = base
        .replace("7.0", "5.0")
        .replace("fractional=[0.13,0.21,0.32]", "fractional=[0.25,0.25,0.25]")
        + "[symmetry]\nmode='auto'\n";
    let electronic = base
        .replace(
            "functional='lda-pz'",
            "functional='lda-pw'\nallow_functional_override=true",
        )
        .replace(
            "kshift=[0.5,0.5,0.5]",
            "kshift=[0.5,0.5,0.5]\nktime_reversal=false",
        );
    let derivatives = "[output]\nforces=true\nstress=true\n";
    for (name, command, text) in [
        ("gamma", "run", gamma.clone()),
        ("k222s", "run", base.clone()),
        (
            "fd",
            "run",
            // FD-FORCES: scf.nbands sets the FD band count; keep the cap (12).
            electronic.replace("nbands=3\n", "")
                + "[occupations]\nmode='fermi_dirac'\nsmearing_ha=0.01\nband_cap=12\n",
        ),
        (
            "fixed-spin",
            "run",
            electronic + "[spin]\nmode='collinear'\nn_up=1\nn_down=1\n",
        ),
        ("symmetry", "run", symmetry.clone()),
        ("gamma-derivatives", "run", gamma.clone() + derivatives),
        ("k-derivatives", "run", base.clone() + derivatives),
        ("symmetry-derivatives", "run", symmetry + derivatives),
        (
            "unverified-derivatives",
            "run",
            base.clone() + "verify='off'\n" + derivatives,
        ),
        (
            "relax-gamma",
            "relax",
            gamma + "[relax]\nmax_steps=2\nforce_tol_ha_bohr=1e-16\nenergy_tol_ha=1e-16\n",
        ),
        (
            "relax-k",
            "relax",
            base + "[relax]\nmax_steps=2\nforce_tol_ha_bohr=1e-16\nenergy_tol_ha=1e-16\n",
        ),
    ] {
        let input = root.join(format!("{name}.toml"));
        let out = root.join(name);
        fs::write(&input, text).unwrap();
        let (old_code, mut old, old_density) = execute(Path::new(&baseline), command, &input, &out);
        fs::remove_dir_all(&out).unwrap();
        let (new_code, mut new, new_density) = execute(
            Path::new(env!("CARGO_BIN_EXE_dft-cli")),
            command,
            &input,
            &out,
        );
        assert_eq!(old_code, new_code, "{name}: exit status");
        // A baseline built with other features (mkl) takes other GEMM/EVD
        // routes and differs in the last bits by construction; name that
        // instead of reporting it as a density change. RUSTFLAGS are not
        // recorded: build both executables with the same RUSTFLAGS.
        if let (Some(a), Some(b)) = (
            old.pointer("/execution/backend"),
            new.pointer("/execution/backend"),
        ) {
            assert_eq!(
                a, b,
                "{name}: build backend differs; build the baseline with the same features and RUSTFLAGS as this test"
            );
        }
        if old_density != new_density {
            panic!(
                "{name}: density bytes: {}",
                density_difference(old_density.as_deref(), new_density.as_deref())
            );
        }
        if command == "relax" {
            let frames = new["relaxation"]["frames"].as_array().unwrap();
            assert_eq!(frames.len(), 3, "{name}: exercise two ionic continuations");
            for frame in &frames[1..] {
                assert_eq!(frame["provenance"], "warm_transfer", "{name}");
                assert!(frame["start_fallback"].is_null(), "{name}");
            }
        }
        normalize(&mut old);
        normalize(&mut new);
        if let Some(diff) = first_difference(&old, &new, "") {
            fs::write(
                root.join(format!("{name}-old.json")),
                serde_json::to_vec_pretty(&old).unwrap(),
            )
            .unwrap();
            fs::write(
                root.join(format!("{name}-new.json")),
                serde_json::to_vec_pretty(&new).unwrap(),
            )
            .unwrap();
            panic!("{name}: {diff}; evidence {}", root.display());
        }
        assert_eq!(
            serde_json::to_vec(&old).unwrap(),
            serde_json::to_vec(&new).unwrap(),
            "{name}: masked JSON bytes"
        );
        println!("{name}: identical masked JSON and density bytes (exit {new_code})");
    }
    fs::remove_dir_all(root).unwrap();
}

#[test]
#[ignore = "release eager/deferred derivative and seed identity; cargo test --release -p dft-cli --test api_identity facade_delivery -- --ignored --nocapture"]
fn facade_delivery_matches_existing_dispatch() {
    use dft_engine::{
        SpeciesPotential,
        refscf::{self, api},
        refstart::StartInput,
    };
    let root = std::env::temp_dir().join(format!("r3-api-delivery-{}", std::process::id()));
    fs::create_dir_all(&root).unwrap();
    let base = synthetic(&root);
    let pp = dft_upf::parse_upf(&fs::read_to_string(root.join("syn.xml")).unwrap()).unwrap();
    let potentials = [SpeciesPotential {
        label: "He",
        potential: &pp,
    }];
    for gamma in [true, false] {
        let text = if gamma {
            base.replace(
                "kgrid=[2,2,2]\nkshift=[0.5,0.5,0.5]",
                "kpoints=[[0.0,0.0,0.0]]",
            )
        } else {
            base.clone()
        };
        let config = dft_cli::config::Config::parse(
            &(text + "threads=1\n[output]\nforces=true\nstress=true\n"),
        )
        .unwrap();
        let options = dft_cli::run::options(&config).unwrap();
        let cell = config.cell().unwrap();
        let atoms: Vec<_> = config
            .system
            .atoms
            .iter()
            .map(|a| dft_geometry::Atom::new(&a.species, a.fractional))
            .collect();
        let input = || StartInput {
            cell: &cell,
            atoms: &atoms,
            potentials: &potentials,
            options: &options,
        };
        let points = config.resolved_points().unwrap();
        let profile = dft_cli::run::reference_profile(&config).unwrap();
        let mut request = api::Request::new(profile, points.clone());
        request.derivatives.forces = true;
        request.derivatives.stress = true;
        request.derivatives.pulay_assessment = dft_engine::PulayAssessmentMode::CutoffScan;
        request.save_ionic_seed = true;
        let eager = api::solve(input(), request).unwrap();
        let eager_run = eager.run.into_scalar().unwrap();
        let mut direct = refscf::delivered::run(input(), profile, &points, None).unwrap();
        let state = refscf::delivered::force_state(input(), &mut direct.outcome).unwrap();
        let forces = dft_engine::forces(&state).unwrap();
        let stress = dft_engine::stress(&state).unwrap();
        assert_eq!(
            eager_run.outcome.result.total_hartree().to_bits(),
            direct.outcome.result.total_hartree().to_bits()
        );
        assert_eq!(
            eager_run.outcome.result.density,
            direct.outcome.result.density
        );
        assert_eq!(eager.derivatives.forces.unwrap().values, forces.values);
        let eager_stress = eager.derivatives.stress.unwrap();
        assert_eq!(eager_stress.values, stress.values);
        assert_eq!(
            eager_stress.pulay.unwrap().correction,
            stress.pulay.unwrap().correction
        );
        let seed = eager.ionic_seed.unwrap();
        // FD-FORCES (2026-11-26): the facade delivers Fermi-Dirac forces.
        let mut fd = api::Request::new(profile, points.clone());
        fd.policy = Some(refscf::electronic::Policy::FermiDirac {
            theta_ha: 0.01,
            band_cap: 12,
        });
        fd.derivatives.forces = true;
        let fd = api::solve(input(), fd).unwrap();
        assert!(matches!(fd.run, api::Run::Electronic(_)));
        let fd_forces = fd.derivatives.forces.expect("FD forces delivered");
        assert!(fd_forces.values.iter().flatten().all(|v| v.is_finite()));
        // A real-Gamma seed with prepared symmetry stays rejected (prepared
        // solves use complex k blocks). API-SEED-GRID: complex-k seeds with
        // prepared symmetry are supported since R7B-CELL-RELAX, but only on
        // the prepared grid. The k seed here comes from the solve without
        // symmetry on its natural grid ([8,8,8]); the off-origin atom gives an
        // inversion with fractional translation 2*tau = (0.26,0.42,0.64), so
        // the symmetric grid is [50,50,25] and the seed must be rejected,
        // never resampled. The supported warm start (seed exported by a
        // symmetric solve on the prepared grid) is covered by
        // `refscf::api::seeded_symmetry_tests::prepared_seed_preserves_stepper_delivery`.
        let mut symmetric = api::Request::new(profile, points.clone());
        symmetric.seed = Some(seed.clone());
        symmetric.symmetry =
            Some(refscf::symmetry::Prepared::new(&cell, &atoms, &points, 1e-5, false).unwrap());
        let expected = if gamma {
            assert!(matches!(seed, api::IonicSeed::Gamma(_)));
            "reference ionic seed: Gamma/k mismatch"
        } else {
            assert!(matches!(seed, api::IonicSeed::K(_)));
            "reference k ionic seed: incompatible cell, basis, grid, mesh or charge"
        };
        match api::solve(input(), symmetric) {
            Err(dft_engine::ScfError::InvalidInput(message)) => assert_eq!(message, expected),
            other => panic!("gamma={gamma}: expected rejection, got {:?}", other.err()),
        }
        let mut request = api::Request::new(profile, points.clone());
        request.seed = Some(seed.clone());
        request.derivatives = api::Derivatives {
            forces: true,
            stress: true,
            deferred: true,
            // PULAY-OPTIN default: no assessment; the tensor below is compared
            // bit for bit with `dft_engine::stress` (assessment on).
            pulay_assessment: dft_engine::PulayAssessmentMode::Off,
        };
        let mut deferred = api::solve(input(), request).unwrap();
        assert!(deferred.derivatives.forces.is_none());
        assert!(deferred.derivatives.stress.is_none());
        let mut warm = deferred.run.into_scalar().unwrap();
        assert!(
            deferred
                .derivatives
                .evaluate_forces(input(), &warm.outcome.result)
                .is_err()
        );
        deferred
            .derivatives
            .prepare(input(), &mut warm.outcome)
            .unwrap();
        assert!(
            deferred
                .derivatives
                .evaluate_stress(input(), &warm.outcome.result)
                .is_err()
        );
        deferred
            .derivatives
            .evaluate_forces(input(), &warm.outcome.result)
            .unwrap();
        deferred
            .derivatives
            .evaluate_stress(input(), &warm.outcome.result)
            .unwrap();
        let mut direct_warm =
            refscf::delivered::run(input(), profile, &points, Some(seed)).unwrap();
        let state = refscf::delivered::force_state(input(), &mut direct_warm.outcome).unwrap();
        assert_eq!(
            deferred.derivatives.forces.take().unwrap().values,
            dft_engine::forces(&state).unwrap().values
        );
        let deferred_stress = deferred.derivatives.stress.take().unwrap();
        assert!(deferred_stress.pulay.is_none());
        assert_eq!(
            deferred_stress.values,
            dft_engine::stress(&state).unwrap().values
        );
        assert_eq!(
            warm.outcome.result.total_hartree().to_bits(),
            direct_warm.outcome.result.total_hartree().to_bits()
        );
        assert_eq!(
            warm.outcome.result.density,
            direct_warm.outcome.result.density
        );
        deferred
            .derivatives
            .evaluate_forces(input(), &warm.outcome.result)
            .unwrap();
        deferred
            .derivatives
            .evaluate_stress(input(), &warm.outcome.result)
            .unwrap();
        assert!(
            deferred.derivatives.forces.is_none() && deferred.derivatives.stress.is_none(),
            "do not repeat completed derivatives"
        );
        println!("gamma={gamma}: eager/deferred force/stress and warm seed identity");
    }
    fs::remove_dir_all(root).unwrap();
}
