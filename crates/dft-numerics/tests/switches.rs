//! Switch manifest gates (work order A1).
use dft_numerics::switches::{self, Id, Kind, Source, TABLE, Value};
use std::path::Path;

/// No source file outside the switches module reads a `BAANDI_*` variable.
#[test]
fn no_baandi_env_reads_outside_switches_module() {
    let crates = Path::new(env!("CARGO_MANIFEST_DIR")).parent().unwrap();
    let module = Path::new(env!("CARGO_MANIFEST_DIR")).join("src/switches.rs");
    let patterns = [
        concat!("env::var", "(\"BAANDI_"),
        concat!("env::var_os", "(\"BAANDI_"),
        concat!("env::var", "(\\\"BAANDI_"),
    ];
    let mut offenders = Vec::new();
    let mut stack: Vec<_> = std::fs::read_dir(crates)
        .unwrap()
        .map(|e| e.unwrap().path().join("src"))
        .filter(|p| p.is_dir())
        .collect();
    let mut scanned = 0;
    while let Some(dir) = stack.pop() {
        for entry in std::fs::read_dir(&dir).unwrap() {
            let path = entry.unwrap().path();
            if path.is_dir() {
                stack.push(path);
            } else if path.extension().is_some_and(|e| e == "rs") && path != module {
                scanned += 1;
                let text = std::fs::read_to_string(&path).unwrap();
                for (i, line) in text.lines().enumerate() {
                    if patterns.iter().any(|p| line.contains(p)) {
                        offenders.push(format!("{}:{}", path.display(), i + 1));
                    }
                }
            }
        }
    }
    assert!(scanned > 50, "scan found only {scanned} files");
    assert!(
        offenders.is_empty(),
        "read BAANDI_* switches through dft_numerics::switches: {offenders:?}"
    );
}

fn rendered(v: &Value) -> String {
    match v {
        Value::Choice(s) => s.to_string(),
        Value::Number(Some(x)) => format!("{x:e}"),
        Value::Number(None) => "none".into(),
        Value::Raw(None) => "unset".into(),
        Value::Raw(Some(s)) => s.clone(),
        Value::Present(b) => b.to_string(),
    }
}

/// Every table default equals the parsed value with the variable unset, and the
/// process accessor agrees where the test environment leaves the variable unset.
#[test]
fn table_defaults_match_accessors_with_environment_unset() {
    for s in TABLE {
        let (value, source) = switches::resolve(s, None);
        assert_eq!(source, Source::Default, "{}", s.name);
        let expected = match s.kind {
            Kind::Choice { default, .. } => {
                assert_eq!(s.default, default, "{}", s.name);
                default.to_string()
            }
            Kind::Number { default, .. } => {
                assert_eq!(s.default.parse::<f64>().ok(), default, "{}", s.name);
                default.map_or("none".into(), |x| format!("{x:e}"))
            }
            Kind::Ortho => {
                assert_eq!(s.default, "dgks:0.1");
                format!("{:e}", 0.1)
            }
            Kind::Override { .. } | Kind::Path => "unset".into(),
            Kind::Presence => "false".into(),
        };
        assert_eq!(rendered(&value), expected, "{}", s.name);
        if std::env::var_os(s.name).is_none() {
            let entry = &switches::manifest()[s.id as usize];
            assert_eq!(entry.value, value, "{}", s.name);
            assert_eq!(entry.source, Source::Default, "{}", s.name);
        }
    }
}

/// The defaults in force before the manifest, typed independently of the table:
/// a table edit that flips one of them fails here.
#[test]
fn legacy_defaults_are_unchanged() {
    let on = |id: Id| switches::resolve(switches::switch(id), None).0 == Value::Choice("on");
    for id in [
        Id::CERT_SCREEN,
        Id::POTENTIAL_REUSE,
        Id::EMPTY_LOOSE,
        Id::DENSITY_PAIR,
        Id::DENSITY_REDUCE,
        Id::NL_ENERGY_REAL,
        Id::ATOMIC_REAL,
        Id::MPI_M2,
        Id::M2_BATCH,
        Id::B1,
        Id::B5,
        Id::B6,
        Id::WIDE_DOT,
        Id::MKL_NONLOCAL,
        Id::APPLY_LEAN,
        Id::MKL_DGEMM,
        Id::MKL_EVD,
        Id::CUTOFF_LADDER_STANDARD,
        Id::STANDARD_RESIDUAL_GUARD,
        Id::PHASE_THREADS,
        Id::STAGE_TIMERS,
        Id::M2_WORKSPACE,
    ] {
        assert!(on(id), "{id:?} must default on");
    }
    for id in [
        Id::ENERGY_IDENTITY,
        Id::RITZ_FUSED,
        Id::NARROW_PROJ,
        Id::NARROW_REDUCED,
        Id::DAVIDSON_STAGES,
        Id::HAPPLY_OVERLAP,
        Id::HAPPLY_FOLD,
        Id::SECOND_VECTOR,
        Id::PROJECT_WIDE,
        Id::B1_NARROW,
        Id::REFERENCE_STRICT,
        Id::M2_B5,
        Id::LEAN_PROJ_CACHE,
        Id::PACK_FIXED,
    ] {
        assert!(!on(id), "{id:?} must default off");
    }
    assert_eq!(switches::choice(Id::LEAN_BORDER), "separate");
    assert!(on(Id::SPIN_SHARE), "SPIN_SHARE defaults on");
    assert!(on(Id::M2_DGEMM) && on(Id::M2_LEAN) && on(Id::REUSE_PSI) && on(Id::M2_SELECTED_ROOTS));
    assert_eq!(
        switches::resolve(switches::switch(Id::VERIFY_REUSE), None).0,
        Value::Choice("half")
    );
    assert_eq!(
        switches::resolve(switches::switch(Id::M2_CAPACITY), None).0,
        Value::Choice("auto")
    );
    assert_eq!(
        switches::resolve(switches::switch(Id::M2_EIGEN), None).0,
        Value::Choice("replicated")
    );
    assert_eq!(
        switches::resolve(switches::switch(Id::M2_ORTHO), None).0,
        Value::Choice("serial")
    );
    assert_eq!(
        switches::resolve(switches::switch(Id::M2_FFT), None).0,
        Value::Choice("rust")
    );
    assert_eq!(
        switches::resolve(switches::switch(Id::M2_SCREEN), None).0,
        Value::Choice("dist")
    );
    assert_eq!(
        switches::resolve(switches::switch(Id::M2_ATOMIC), None).0,
        Value::Choice("split")
    );
    let value = |id: Id| switches::resolve(switches::switch(id), None).0;
    assert_eq!(value(Id::NL_ENERGY), Value::Choice("fast"));
    assert_eq!(value(Id::B1_SKETCH), Value::Choice("count"));
    assert_eq!(value(Id::LADDER_STANDARD_SEED), Value::Choice("coarse"));
    assert_eq!(value(Id::LADDER_SECANT), Value::Choice("off"));
    assert_eq!(value(Id::LADDER_COARSE_STOP), Value::Number(Some(1e-4)));
    assert_eq!(value(Id::LADDER_COARSE_RULE), Value::Choice("fixed"));
    assert_eq!(value(Id::LADDER_FIRST_TOL), Value::Choice("auto"));
    assert_eq!(value(Id::STANDARD_STOP), Value::Choice("plain"));
    assert_eq!(value(Id::ROTATE_KERNEL), Value::Choice("legacy"));
    assert_eq!(
        value(Id::PURE_GEMM),
        Value::Choice(if dft_numerics::MKL_BACKEND {
            "nano"
        } else {
            "packed"
        })
    );
    assert_eq!(
        value(Id::ORTHO_BLOCK),
        Value::Choice(if dft_numerics::MKL_BACKEND {
            "off"
        } else {
            "on"
        })
    );
    assert_eq!(
        value(Id::PAIR_FFT),
        Value::Choice(if dft_numerics::MKL_BACKEND {
            "pruned"
        } else {
            "grid"
        })
    );
    assert_eq!(value(Id::DAVIDSON), Value::Choice("lean"));
    assert_eq!(value(Id::NONLOCAL_BLOCK), Value::Choice("auto"));
    assert_eq!(
        value(Id::MGS),
        Value::Choice(if dft_numerics::MKL_BACKEND {
            "blocked"
        } else {
            "classic"
        })
    );
    assert_eq!(value(Id::ORTHO), Value::Number(Some(0.1)));
    assert_eq!(value(Id::BUFFER_TOL), Value::Number(Some(1e-3)));
    assert_eq!(value(Id::EIGEN_FACTOR), Value::Number(Some(0.1)));
    assert_eq!(value(Id::MIX_GROWTH_RESET), Value::Number(None));
    assert_eq!(value(Id::B1_MARGIN), Value::Number(Some(0.5)));
    assert_eq!(value(Id::B1_ROWS), Value::Number(None));
    assert_eq!(value(Id::CUTOFF_LADDER), Value::Raw(None));
    assert_eq!(value(Id::RELAX_EXTRAPOLATION), Value::Raw(None));
    assert_eq!(value(Id::DGEMM_TIMERS), Value::Present(false));
    // MPI-KEPT-START stays opt-in until dev08.
    assert_eq!(value(Id::M2_REUSE), Value::Choice("on"));
    // MPI-SCALING-GAP switches stay opt-in until dev08.
    assert_eq!(value(Id::M2_STICKS), Value::Choice("mirror"));
    assert_eq!(value(Id::M2_START), Value::Choice("gram"));
    // MPI-OVERLAP-DAVIDSON stays opt-in until dev08.
    assert_eq!(value(Id::M2_LEAN_ALGEBRA), Value::Choice("on"));
    // Adopted 2026-11-17 on dev08 evidence (ledger M2-LAZY-DEV08-STAGE0/-STAGE1, MPI-FS-PANEL-GS).
    assert_eq!(value(Id::M2_LAZY), Value::Choice("on"));
    assert_eq!(value(Id::M2_PANEL_GS), Value::Choice("fused"));
    // MPI-REMAINING grid switches stay opt-in until dev08.
    assert_eq!(value(Id::MIX_INCREMENTAL), Value::Choice("on"));
    assert_eq!(value(Id::STOP_METRIC), Value::Choice("kept"));
    // BREAKTHROUGH-1-RESTART R1-b stays off (REJECT at stage 0 and 1).
    assert_eq!(value(Id::RESTART_KEEP), Value::Choice("off"));
    // BREAKTHROUGH-3 bounded subspace SCF stays off (rejected at stage 1, ledger BT3-*).
    assert_eq!(value(Id::SCF_ADVANCE), Value::Choice("off"));
    assert_eq!(value(Id::SCF_ADVANCE_BLOCKS), Value::Number(Some(2.0)));
    // BREAKTHROUGH-0-INSTRUMENT measurement switches stay off (results unchanged either way).
    assert_eq!(value(Id::PHASE_PROFILE), Value::Choice("off"));
    assert_eq!(value(Id::FIRST_BLOCK), Value::Choice("off"));
    assert_eq!(value(Id::SNAPSHOT_DIR), Value::Raw(None));
    // BREAKTHROUGH-4-OTHER-AND-INIT: wait probe and the three removal switches
    // stay off until dev08 (all bit-identical; the probe changes timing only).
    assert_eq!(value(Id::PHASE_WAIT), Value::Choice("off"));
    assert_eq!(value(Id::EWALD_CACHE), Value::Choice("on"));
    assert_eq!(value(Id::M2_LIGHT_STATE), Value::Choice("on"));
    assert_eq!(value(Id::M2_DENSITY), Value::Choice("fused"));
    assert_eq!(value(Id::M2_LAYOUT_CACHE), Value::Choice("off"));
    // ISSUE-20: lazy serial operator of M2 states, on since dev08 (si64 -1.22 %, bit-identical).
    assert_eq!(value(Id::M2_LAZY_SERIAL), Value::Choice("on"));
    // LADDER-IMPORT-SPLIT: split import check (same verdict, bit-identical results).
    assert_eq!(value(Id::IMPORT_CHECK), Value::Choice("split"));
    // LAZY-RESULT: lazy delivered orbitals when the caller reads none.
    assert_eq!(value(Id::M2_LAZY_RESULT), Value::Choice("on"));
    assert_eq!(TABLE.len(), 101);
}
