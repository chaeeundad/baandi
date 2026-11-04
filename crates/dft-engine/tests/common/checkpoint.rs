//! Shared synthetic checkpoint fixtures; no real pseudopotential files.
use crate::common;
use dft_engine::*;
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use std::time::Duration;
pub fn normalized(mut r: ScfResult) -> ScfResult {
    r.setup_time = Duration::ZERO;
    r.total_time = Duration::ZERO;
    for h in &mut r.history {
        h.timings = PhaseTimings::default();
        if let Some(c) = &mut h.fixed_point_certificate {
            c.timings = PhaseTimings::default();
        }
    }
    r
}
pub fn assert_identical(a: &ScfResult, b: &ScfResult) {
    assert_eq!(format!("{a:?}"), format!("{b:?}"));
    assert_eq!(a.total_hartree().to_bits(), b.total_hartree().to_bits());
    for (x, y) in a.history.iter().zip(&b.history) {
        assert_eq!(
            x.energy.total_hartree().to_bits(),
            y.energy.total_hartree().to_bits()
        );
        assert_eq!(x.density_residual.to_bits(), y.density_residual.to_bits());
    }
    assert_eq!(a.density.shape(), b.density.shape());
    for (x, y) in a.density.iter().zip(&b.density) {
        assert_eq!(x.to_bits(), y.to_bits());
    }
    for (x, y) in a.k_points.iter().zip(&b.k_points) {
        assert_eq!(x.solution, y.solution);
        for j in 0..x.solution.eigenvectors.cols() {
            for (v, w) in x
                .solution
                .eigenvectors
                .column(j)
                .iter()
                .zip(y.solution.eigenvectors.column(j))
            {
                assert_eq!(v.re.to_bits(), w.re.to_bits());
                assert_eq!(v.im.to_bits(), w.im.to_bits());
            }
        }
    }
}

pub fn fixture(kind: usize, mixing: Mixing, threads: usize, max_iterations: usize) -> ScfProblem {
    let p = common::potential(2.0, 1.0, true, 401);
    let q = common::potential(2.0, 1.3, false, 401);
    let cell = Cell::from_user_rows_bohr(if kind == 2 {
        [[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]]
    } else {
        [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]
    })
    .unwrap();
    let mut atoms = vec![Atom::new("Xx", [0.13, 0.21, 0.32])];
    if kind == 2 {
        atoms.push(Atom::new("Yy", [0.6, 0.7, 0.8]));
    }
    let points: Vec<_> = (0..if kind == 1 { 4 } else { 1 })
        .map(|i| KPoint {
            fractional: [0.025 * i as f64, 0.015 * i as f64, 0.],
            weight: if kind == 1 { 0.25 } else { 1. },
        })
        .collect();
    let species = [
        SpeciesPotential {
            label: "Xx",
            potential: &p,
        },
        SpeciesPotential {
            label: "Yy",
            potential: &q,
        },
    ];
    ScfProblem::with_kpoints(
        &cell,
        &atoms,
        &species[..if kind == 2 { 2 } else { 1 }],
        ScfOptions {
            threads,
            mixing,
            max_iterations,
            cutoff_hartree: 2.5,
            bands: Some(if kind == 2 { 4 } else { 3 }),
            ..Default::default()
        },
        &points,
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap()
}

pub fn identity(
    kind: usize,
    threads: usize,
) -> (checkpoint::ModelIdentity, checkpoint::ExecutionIdentity) {
    let cell = Cell::from_user_rows_bohr(if kind == 2 {
        [[5., 0., 0.], [0.5, 4., 0.], [0.2, 0.4, 6.]]
    } else {
        [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]
    })
    .unwrap();
    let model = checkpoint::ModelIdentity {
        cell_bohr: cell.lattice_columns_bohr(),
        atoms_fractional: if kind == 2 {
            vec![[0.13, 0.21, 0.32], [0.6, 0.7, 0.8]]
        } else {
            vec![[0.13, 0.21, 0.32]]
        },
        species: if kind == 2 {
            vec!["Xx".into(), "Yy".into()]
        } else {
            vec!["Xx".into()]
        },
        pseudopotential_sha256: vec!["0".repeat(64); if kind == 2 { 2 } else { 1 }],
        cutoff_hartree: 2.5,
        kpoints: (0..if kind == 1 { 4 } else { 1 })
            .map(|i| {
                [
                    0.025 * i as f64,
                    0.015 * i as f64,
                    0.,
                    if kind == 1 { 0.25 } else { 1. },
                ]
            })
            .collect(),
        electrons: if kind == 2 { 4. } else { 2. },
        xc: "lda-pz".into(),
        spin: "unpolarized".into(),
        environment: "none".into(),
    };
    (
        model,
        checkpoint::ExecutionIdentity {
            build_sha256: "1".repeat(64),
            platform: "synthetic-same-test-build".into(),
            threads,
        },
    )
}
pub fn drain(s: &mut ScfSession) {
    if let Some(q) = s.events() {
        while q.pop().is_some() {}
    }
}
pub fn finish(mut s: ScfSession) -> ScfResult {
    while s.snapshot().status == StepStatus::Ready {
        s.step().unwrap();
        drain(&mut s);
    }
    let out = s.into_result().unwrap();
    assert!(out.converged);
    normalized(out.result)
}
