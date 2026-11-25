//! Strict Gamma oracle identity, captured from pre-R10A main 28041772.
//! Pure-Rust macOS/aarch64, Rust 1.99.1, registry defaults, one thread.
//!
//! The reported dev08 Si8/Si.pz-vbc/30 Ry/36^3 regression was
//! c03f5261b8862aca -> c03f5261b8862ad4 (Ha binary64), SCF 11 -> 9,
//! final-map Davidson 16 -> 7. Those are relayed numbers, not this synthetic
//! fixture's goldens. Local Si8 pre-R10A and HEAD both give c03f5261b8862aca,
//! SCF 10, final-map Davidson 6 with explicit Pulay beta=0.5/history=6.
//! This test protects the retained default ladder; it does not reproduce or
//! resolve the reported dev08 discrepancy.

#[allow(dead_code)]
mod common;

use dft_engine::{ScfOptions, ScfProblem, SpeciesPotential, cutoff_ladder};
use dft_geometry::{Atom, Cell};
use dft_numerics::switches::RunConfig;

#[test]
fn strict_gamma_default_ladder_identity() {
    // Backend-specific goldens must come from the same historical build, not
    // from a tolerance relaxation or recapture of a changed implementation.
    if (std::env::consts::OS, std::env::consts::ARCH) != ("macos", "aarch64")
        || dft_numerics::MKL_BACKEND
    {
        eprintln!("No pre-R10A strict-oracle golden for this platform/backend");
        return;
    }
    let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
    let pp = common::potential(2., 1., true, 161);
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let potentials = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let options = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(5),
        ..ScfOptions::strict_reference()
    };
    let problem = ScfProblem::new(&cell, &atoms, &potentials, options).unwrap();
    // Historical run_ladder reported used=true and four coarse iterations.
    // Bypassing it gives bfe3f626a6467a16 and 12 fine maps on this fixture,
    // so these goldens detect deletion of the default coarse initialization.
    let outcome = cutoff_ladder::run_oracle(problem, &potentials).unwrap();
    assert!(outcome.converged);
    let result = outcome.result;
    assert!(result.gamma_real_used);
    assert_eq!(result.total_hartree().to_bits(), 0xbfe3_f626_a646_7a14);
    assert_eq!(result.history.len(), 9);
    assert_eq!(
        result
            .history
            .iter()
            .map(|h| h.davidson_iterations)
            .collect::<Vec<_>>(),
        [2, 2, 4, 6, 2, 3, 1, 3, 3]
    );
    // Ordered binary64 fingerprint: ndarray iteration order, wrapping multiply
    // then XOR each density sample's bits (not a cryptographic digest).
    let density_bits = result.density.iter().fold(0_u64, |hash, value| {
        hash.wrapping_mul(0x100_0000_01b3) ^ value.to_bits()
    });
    assert_eq!(density_bits, 0xe3c2_28c4_058c_b267);
}
