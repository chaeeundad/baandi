//! GK14: error of the uniform radial transform table (`RadialTransformTable`)
//! against the direct quadrature it tabulates, as a function of the node
//! spacing and the Lagrange stencil.
//!
//! Metric per channel: max over q in [0, q_max] of |table - quadrature|
//! divided by max over q of |quadrature| (relative to the channel's peak).
//! q is sampled on an irregular grid (step 0.0013 Bohr^-1, not commensurate
//! with any table spacing).
//!
//! Real files (ignored, prints the scan):
//! DFT_NC_DIR=/abs/pseudodojo/ex [DFT_SI_UPF=/abs/Si.pz-vbc.UPF] \
//!   cargo test --release -p dft-upf --test all radial_table::scan -- --ignored --nocapture
#[allow(dead_code)]
#[path = "../common/mod.rs"]
mod common;
use dft_upf::{
    Pseudopotential, parse_upf,
    radial::{RadialEvenRule, RadialTransformTable, projector_transform_with_rule},
};

const RULE: RadialEvenRule = RadialEvenRule::SimpsonQuadraticEnd;

/// (label, l, radial_times_r) of every projector and atomic orbital.
fn channels(pp: &Pseudopotential) -> Vec<(String, usize, Vec<f64>)> {
    let mut out: Vec<_> = pp
        .projectors()
        .iter()
        .enumerate()
        .map(|(i, p)| (format!("beta{i}"), p.l(), p.radial_times_r().to_vec()))
        .collect();
    for (i, c) in pp.atomic_wavefunctions().unwrap_or(&[]).iter().enumerate() {
        if c.l <= 3 {
            out.push((format!("chi{i}"), c.l, c.radial_times_r.clone()));
        }
    }
    out
}

/// Worst relative error of one channel's table.
fn relative_error(
    pp: &Pseudopotential,
    l: usize,
    f: &[f64],
    q_max: f64,
    step: f64,
    points: usize,
) -> f64 {
    let t = RadialTransformTable::new(pp.mesh(), l, f, RULE, q_max, step, points).unwrap();
    let (mut worst, mut peak) = (0.0_f64, 0.0_f64);
    let samples = (q_max / 0.0013) as usize;
    for i in 0..=samples {
        let q = (i as f64 * 0.0013).min(q_max);
        let exact = projector_transform_with_rule(pp.mesh(), l, f, q, RULE).unwrap();
        peak = peak.max(exact.abs());
        worst = worst.max((t.eval(q).unwrap() - exact).abs());
    }
    worst / peak
}

/// The default table (step 0.005, 8 points) on the synthetic fixture: every
/// projector and orbital within 1e-12 of its peak up to q = 10 Bohr^-1
/// (cutoff 50 Ha). Frozen cap; the measured value is printed.
#[test]
fn default_table_on_synthetic_projectors() {
    let pp = parse_upf(&common::fixture_with_options(
        1201,
        &common::Options {
            wavefunctions: true,
            ..Default::default()
        },
    ))
    .unwrap();
    for (name, l, f) in channels(&pp) {
        let e = relative_error(&pp, l, &f, 10.0, 0.005, 8);
        println!("GK14 synthetic {name} l={l} step=0.005 points=8 rel_error={e:.3e}");
        assert!(e < 1e-12, "{name}: {e:e}");
    }
}

#[test]
#[ignore = "requires DFT_NC_DIR; prints the error-vs-spacing scan"]
fn scan() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let mut files: Vec<_> = ["Ge", "Te", "Si", "O", "Mg", "Fe", "Ir", "I", "B"]
        .iter()
        .map(|e| dir.join(format!("{e}.upf")))
        .collect();
    if let Some(si) = std::env::var_os("DFT_SI_UPF") {
        files.push(si.into());
    }
    let q_max = 11.0; // |k+G| sphere of a 60 Ha cutoff
    for path in files {
        let Ok(pp) = dft_upf::from_file(&path) else {
            println!("GK14 skip {}", path.display());
            continue;
        };
        let name = path.file_stem().unwrap().to_string_lossy().to_string();
        for points in [4, 6, 8] {
            for step in [0.02, 0.01, 0.005, 0.0025] {
                let worst = channels(&pp)
                    .iter()
                    .map(|(_, l, f)| relative_error(&pp, *l, f, q_max, step, points))
                    .fold(0.0_f64, f64::max);
                println!(
                    "GK14_SCAN file={name} points={points} step={step} worst_rel_error={worst:.3e}"
                );
            }
        }
    }
}
