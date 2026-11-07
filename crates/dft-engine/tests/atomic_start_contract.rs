#[allow(dead_code)]
#[path = "../../dft-cli/tests/common/mod.rs"]
mod common;
use dft_engine::*;
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_physics::data::ComplexMatrix;

fn potential(orbitals: bool, duplicate: bool) -> dft_upf::Pseudopotential {
    let mut text = common::gaussian_with(2., 1., false, 101);
    if orbitals {
        // Use the fixture's actual mesh, independent of its mesh generation.
        let pp = dft_upf::parse_upf(&text).unwrap();
        let values = pp
            .mesh()
            .r()
            .iter()
            .map(|r| format!("{}", r * (-0.7 * r * r).exp()))
            .collect::<Vec<_>>()
            .join(" ");
        let count = if duplicate { 2 } else { 1 };
        text = text.replace("number_of_wfc=\"0\"", &format!("number_of_wfc=\"{count}\""));
        let channels = (1..=count)
            .map(|i| {
                format!("<PP_CHI.{i} l=\"0\" label=\"s\" occupation=\"0\">{values}</PP_CHI.{i}>")
            })
            .collect::<String>();
        text = text.replace("</UPF>", &format!("<PP_PSWFC>{channels}</PP_PSWFC></UPF>"));
    }
    dft_upf::parse_upf(&text).unwrap()
}
fn problem(
    orbitals: bool,
    duplicate: bool,
    mode: StartWavefunctions,
    real: GammaReal,
    k: [f64; 3],
    threads: usize,
) -> ScfProblem {
    let p = potential(orbitals, duplicate);
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    ScfProblem::with_kpoints(
        &cell,
        &[Atom::new("Xx", [0.13, 0.21, 0.32])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &p,
        }],
        ScfOptions {
            start_wavefunctions: mode,
            gamma_real: real,
            threads,
            cutoff_hartree: 2.,
            bands: Some(5),
            ..Default::default()
        },
        &[KPoint {
            fractional: k,
            weight: 1.,
        }],
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-6,
        },
    )
    .unwrap()
}

#[test]
fn missing_section_falls_back_with_identical_map() {
    let mut free = problem(
        false,
        false,
        StartWavefunctions::FreeElectron,
        GammaReal::Off,
        [0.; 3],
        1,
    );
    let mut atomic = problem(
        false,
        false,
        StartWavefunctions::Atomic,
        GammaReal::Off,
        [0.; 3],
        1,
    );
    let n = free.uniform_density();
    let a = free.step_kpoints(&n, None, 1e-5).unwrap();
    let b = atomic.step_kpoints(&n, None, 1e-5).unwrap();
    assert_eq!(a.density, b.density);
    assert_eq!(
        a.energy.total_hartree().to_bits(),
        b.energy.total_hartree().to_bits()
    );
    assert!(
        atomic
            .atomic_start_report()
            .fallback_reason
            .as_deref()
            .unwrap()
            .contains("PP_PSWFC")
    );
    assert_eq!(atomic.atomic_start_report().initial_h_applications, 0);
}

#[test]
fn rank_failure_and_supplied_vectors_take_precedence() {
    let mut bad = problem(
        true,
        true,
        StartWavefunctions::Atomic,
        GammaReal::Off,
        [0.; 3],
        1,
    );
    let n = bad.uniform_density();
    assert!(
        bad.step_kpoints(&n, None, 1e-5)
            .err()
            .unwrap()
            .to_string()
            .contains("rank deficient")
    );
    let mut free = problem(
        false,
        false,
        StartWavefunctions::FreeElectron,
        GammaReal::Off,
        [0.; 3],
        1,
    );
    let first = free.step_kpoints(&n, None, 1e-5).unwrap();
    let warm: Vec<ComplexMatrix> = first
        .k_points
        .iter()
        .map(|k| k.solution.eigenvectors.clone())
        .collect();
    let a = free.step_kpoints(&n, Some(&warm), 1e-5).unwrap();
    let b = bad.step_kpoints(&n, Some(&warm), 1e-5).unwrap();
    assert_eq!(a.density, b.density);
    assert_eq!(bad.atomic_start_report().initial_h_applications, 0);
    assert!(
        bad.atomic_start_report()
            .fallback_reason
            .as_deref()
            .unwrap()
            .contains("supplied")
    );
}

#[test]
fn completion_gamma_orc_and_shifted_k_are_deterministic() {
    for (real, k) in [
        (GammaReal::Off, [0.; 3]),
        (GammaReal::On, [0.; 3]),
        (GammaReal::Off, [0.25, 0.125, 0.]),
    ] {
        let mut reference = None;
        for threads in [1, 2, 4, 8, 12] {
            let mut p = problem(true, false, StartWavefunctions::Atomic, real, k, threads);
            let step = p.step_kpoints(&p.uniform_density(), None, 1e-6).unwrap();
            assert_eq!(p.atomic_start_report().atomic_orbitals, 1);
            assert_eq!(p.atomic_start_report().initial_h_applications, 5);
            assert_eq!(p.atomic_start_report().evd_orders, [5]);
            assert!(p.atomic_start_report().mgs_projection_macs > 0);
            let bits: Vec<_> = step.density.iter().map(|v| v.to_bits()).collect();
            if let Some(expected) = &reference {
                assert_eq!(&bits, expected);
            } else {
                reference = Some(bits);
            }
        }
    }
}

#[test]
fn exact_restart_compares_start_option_and_preserves_work() {
    use dft_engine::checkpoint::*;
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let model = ModelIdentity {
        cell_bohr: cell.lattice_columns_bohr(),
        atoms_fractional: vec![[0.13, 0.21, 0.32]],
        species: vec!["Xx".into()],
        pseudopotential_sha256: vec!["0".repeat(64)],
        cutoff_hartree: 2.,
        kpoints: vec![[0., 0., 0., 1.]],
        electrons: 2.,
        xc: "lda-pz".into(),
        spin: "unpolarized".into(),
        environment: "none".into(),
    };
    let execution = ExecutionIdentity {
        build_sha256: "1".repeat(64),
        platform: "test".into(),
        threads: 1,
    };
    let make = |mode| problem(true, false, mode, GammaReal::Off, [0.; 3], 1);
    let mut a = ScfSession::new(
        make(StartWavefunctions::Atomic),
        CompletionMode::RetainLast,
        ExecutionConfig::default(),
    )
    .unwrap();
    a.initialize().unwrap();
    a.step().unwrap();
    let cp = a
        .checkpoint(
            CheckpointKind::ExactContinuation,
            model.clone(),
            execution.clone(),
        )
        .unwrap();
    assert_eq!(
        cp.continuation
            .as_ref()
            .unwrap()
            .atomic_start
            .as_ref()
            .unwrap()
            .initial_h_applications,
        5
    );
    assert!(matches!(
        ScfSession::restore_exact(
            make(StartWavefunctions::FreeElectron),
            &cp,
            &model,
            &execution
        ),
        Err(Rejection::Options)
    ));
    let mut b =
        ScfSession::restore_exact(make(StartWavefunctions::Atomic), &cp, &model, &execution)
            .unwrap();
    a.step().unwrap();
    b.step().unwrap();
    assert_eq!(a.density(), b.density());
    let warm = a
        .checkpoint(CheckpointKind::WarmStart, model.clone(), execution.clone())
        .unwrap();
    let config = ExecutionConfig {
        calc_id: "atomic-warm".into(),
        ..Default::default()
    };
    let mut warm_atomic = ScfSession::warm_start(
        make(StartWavefunctions::Atomic),
        &warm,
        &model,
        &execution,
        CompletionMode::RetainLast,
        config.clone(),
    )
    .unwrap();
    let mut warm_free = ScfSession::warm_start(
        make(StartWavefunctions::FreeElectron),
        &warm,
        &model,
        &execution,
        CompletionMode::RetainLast,
        config,
    )
    .unwrap();
    warm_atomic.step().unwrap();
    warm_free.step().unwrap();
    assert_eq!(warm_atomic.density(), warm_free.density());
    let warm_cp = warm_atomic
        .checkpoint(
            CheckpointKind::ExactContinuation,
            model.clone(),
            execution.clone(),
        )
        .unwrap();
    assert_eq!(
        warm_cp
            .continuation
            .unwrap()
            .atomic_start
            .unwrap()
            .initial_h_applications,
        0
    );
    let resumed = b
        .checkpoint(CheckpointKind::ExactContinuation, model, execution)
        .unwrap();
    assert_eq!(
        resumed
            .continuation
            .unwrap()
            .atomic_start
            .unwrap()
            .initial_h_applications,
        5
    );
}
