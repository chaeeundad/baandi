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
/// a table edit that flips one of them fails here. R2 explicitly adopts the
/// reference production driver; R1B adopts ADR 0046 XC screening.
#[test]
fn legacy_defaults_are_unchanged() {
    let on = |id: Id| switches::resolve(switches::switch(id), None).0 == Value::Choice("on");
    for id in [
        Id::M2_DGEMM,
        Id::MKL_DGEMM,
        Id::MKL_EVD,
        Id::STANDARD_RESIDUAL_GUARD,
        Id::STAGE_TIMERS,
    ] {
        assert!(on(id), "{id:?} must default on");
    }
    for id in [
        Id::ENERGY_IDENTITY,
        Id::PACK_FIXED,
        Id::DAVIDSON_STAGES,
        Id::REFERENCE_STRICT,
        Id::DAVIDSON_PROFILE,
    ] {
        assert!(!on(id), "{id:?} must default off");
    }
    assert_eq!(
        switches::table_of(Id::DAVIDSON_PROFILE),
        switches::Table::Profile
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
        value(Id::PAIR_FFT),
        Value::Choice(if dft_numerics::MKL_BACKEND {
            "pruned"
        } else {
            "grid"
        })
    );

    assert_eq!(value(Id::NONLOCAL_BLOCK), Value::Choice("auto"));

    assert_eq!(value(Id::REFSCF_LADDER_RATIO), Value::Number(Some(100.0)));

    assert_eq!(value(Id::RELAX_EXTRAPOLATION), Value::Raw(None));
    assert_eq!(value(Id::DGEMM_TIMERS), Value::Present(false));
    // MPI-KEPT-START stays opt-in until dev08.
    // MPI-SCALING-GAP switches stay opt-in until dev08.
    assert_eq!(value(Id::M2_STICKS), Value::Choice("mirror"));

    // MPI-OVERLAP-DAVIDSON stays opt-in until dev08.

    // Adopted 2026-11-17 on dev08 evidence (ledger M2-LAZY-DEV08-STAGE0/-STAGE1, MPI-FS-PANEL-GS).

    // MPI-REMAINING grid switches stay opt-in until dev08.

    // BREAKTHROUGH-1-RESTART R1-b stays off (REJECT at stage 0 and 1).

    // BREAKTHROUGH-3 bounded subspace SCF stays off (rejected at stage 1, ledger BT3-*).
    // BREAKTHROUGH-0-INSTRUMENT measurement switches stay off (results unchanged either way).
    assert_eq!(value(Id::PHASE_PROFILE), Value::Choice("off"));
    assert_eq!(value(Id::SNAPSHOT_DIR), Value::Raw(None));
    // BREAKTHROUGH-4-OTHER-AND-INIT: wait probe and the three removal switches
    // stay off until dev08 (all bit-identical; the probe changes timing only).
    assert_eq!(value(Id::PHASE_WAIT), Value::Choice("off"));

    // ISSUE-20: lazy serial operator of M2 states, on since dev08 (si64 -1.22 %, bit-identical).

    // LADDER-IMPORT-SPLIT: split import check (same verdict, bit-identical results).

    // LAZY-RESULT: lazy delivered orbitals when the caller reads none.
    assert_eq!(value(Id::M2_LAZY_RESULT), Value::Choice("on"));
    // M2-GRID: distributed full-grid stages of the M2 path, rounding class, off until dev08.
    assert_eq!(value(Id::M2_GRID), Value::Choice("on"));

    // INIT-EWALD-FAST / INIT-IONIC-FAST: rounding-class init schemes, on since 2026-11-18 (independent reference: the default Ewald sum is 1.794e-9 Ha low at si64).
    assert_eq!(value(Id::EWALD_FAST), Value::Choice("on"));
    assert_eq!(value(Id::IONIC_FAST), Value::Choice("on"));
    // SCF-TRAJECTORY: Coulomb guard of the L2 Pulay fit, off until dev08.

    // M2-OWNED: owned potential and density prototype, off until dev08.
    assert_eq!(value(Id::M2_LOCAL_FIELDS), Value::Choice("on"));
    assert_eq!(value(Id::M2_LOCAL_DENSITY), Value::Choice("on"));
    assert_eq!(value(Id::PHASE_VALIDATE), Value::Choice("off"));
    assert_eq!(value(Id::INIT_PROFILE), Value::Choice("off"));
    assert_eq!(value(Id::FFT_PROFILE), Value::Choice("off"));

    assert_eq!(value(Id::FFT_KEPT_FUSE), Value::Choice("off"));

    assert_eq!(value(Id::COLLECTIVE_TRACE), Value::Choice("off"));
    assert_eq!(value(Id::BINARY_HASH_CACHE), Value::Choice("on"));
    assert_eq!(value(Id::XC_SCREEN_NEGATIVE), Value::Choice("on"));
    assert_eq!(value(Id::REFSCF_REAL_START), Value::Choice("off"));
    assert_eq!(value(Id::REFSCF_CANONICAL_CLUSTER), Value::Choice("off"));
    assert_eq!(value(Id::REFSCF_SPIN_PULAY), Value::Choice("off"));
    assert_eq!(value(Id::K_LEAN_ENERGY), Value::Choice("on"));
    assert_eq!(value(Id::K_START_MEMO), Value::Choice("on"));
    assert_eq!(value(Id::CDAV_REUSE), Value::Choice("on"));
    assert_eq!(value(Id::K_NL_TABLE_MEMO), Value::Choice("on"));
    assert_eq!(value(Id::K_START_SPIN_SHARE), Value::Choice("on"));
    assert_eq!(value(Id::CDAV_BORDER_GRAM), Value::Choice("on"));
    assert_eq!(value(Id::RADIAL_TABLE_START), Value::Choice("on"));
    assert_eq!(value(Id::RADIAL_TABLE_NL), Value::Choice("on"));
    assert_eq!(value(Id::CDAV_BORDER_FACTOR), Value::Choice("on"));
    assert_eq!(value(Id::CDAV_STREAM), Value::Choice("on"));
    assert_eq!(value(Id::K_FIRST_EMPTY_RULE), Value::Choice("off"));
    assert_eq!(value(Id::K_REDO_KEPT_HX), Value::Choice("on"));
    assert_eq!(value(Id::K_START_KEPT_HX), Value::Choice("on"));
    assert_eq!(value(Id::REFSCF_K_POOLS), Value::Choice("auto"));
    assert_eq!(value(Id::REFSCF_K_NPOOL), Value::Number(None));
    assert_eq!(value(Id::RADIAL_TABLE_FIELDS), Value::Choice("on"));
    assert_eq!(value(Id::FFT_ONE_RANK_DIRECT), Value::Choice("on"));
    assert_eq!(value(Id::K_KEPT_SLABS), Value::Choice("off"));
    assert_eq!(value(Id::K_KEPT_SLAB_MB), Value::Number(Some(1024.)));
    assert_eq!(value(Id::SPIN_FD_START_KEPT_HX), Value::Choice("on"));
    assert_eq!(value(Id::SPIN_FD_REDO_KEPT_HX), Value::Choice("on"));
    assert_eq!(TABLE.len(), 83);
}

#[test]
fn init_profile_off_count_only_and_scoped() {
    use dft_numerics::switches::{InitProfile, RunConfig, RunCounter as C};
    let off = RunConfig::from_pairs(&[]).unwrap();
    let _scope = off.enter();
    let before = off.counters();
    assert!(InitProfile::current().is_none());
    assert_eq!(off.counters(), before);
    let counts = RunConfig::from_pairs(&[("init_profile", "on"), ("stage_timers", "off")]).unwrap();
    {
        let _scope = counts.enter();
        let p = InitProfile::current().unwrap();
        let _span = p.span(C::InitProjectorsCalls, C::InitProjectorsNanos);
        p.count(C::InitFactorRows, 37);
    }
    assert_eq!(counts.counters()[C::InitProjectorsCalls as usize], 1);
    assert_eq!(counts.counters()[C::InitFactorRows as usize], 37);
    assert_eq!(counts.counters()[C::InitProjectorsNanos as usize], 0);
    assert_eq!(off.counters(), before);
    assert!(InitProfile::current().is_none());
}

#[test]
fn glue_profile_requires_both_switches_and_scopes_its_counters() {
    use dft_numerics::switches::{GlueProfile, RunConfig, RunCounter as C};
    for pairs in [
        vec![],
        vec![("glue_profile", "on")],
        vec![("phase_profile", "on")],
    ] {
        let config = RunConfig::from_pairs(&pairs).unwrap();
        let _scope = config.enter();
        assert!(GlueProfile::current().is_none());
        assert_eq!(config.counters()[C::GlueDistMixCalls as usize], 0);
    }
    let cfg = RunConfig::from_pairs(&[
        ("glue_profile", "on"),
        ("phase_profile", "on"),
        ("stage_timers", "off"),
    ])
    .unwrap();
    {
        let _scope = cfg.enter();
        let p = GlueProfile::current().unwrap();
        let _span = p.span(C::GlueDistMixCalls, C::GlueDistMixNanos);
        p.count(C::GlueDistPoints, 37);
    }
    assert_eq!(cfg.counters()[C::GlueDistMixCalls as usize], 1);
    assert_eq!(cfg.counters()[C::GlueDistPoints as usize], 37);
    assert_eq!(cfg.counters()[C::GlueDistMixNanos as usize], 0);
}

// R6a adopts the measured G4B combination while retaining explicit A/B overrides.
#[test]
fn moment_fd_defaults_on_and_off_remains_available() {
    {
        let id = Id::REFSCF_SPIN_FD_MOMENT;
        let switch = switches::switch(id);
        assert_eq!(
            switches::resolve(switch, None),
            (Value::Choice("on"), Source::Default)
        );
        assert_eq!(
            switches::resolve(switch, Some(std::ffi::OsStr::new("off"))).0,
            Value::Choice("off")
        );
    }
}
