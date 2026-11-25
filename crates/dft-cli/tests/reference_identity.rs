//! R0 identity capture, not a physical Si or cross-driver accuracy oracle.
//! Pure-Rust macOS/aarch64, pinned toolchain, registry defaults, one thread.
//! Captured at source 181275fbad3856e0eeac60e74bfd1d34cda37939, Rust 1.99.1.
//! Density hash: ndarray iteration order, each f64 encoded little-endian.
#[allow(dead_code)]
mod common;

use dft_cli::{config::Config, run};
use dft_engine::{ScfProblem, SpeciesPotential, refscf, refstart::StartInput};
use dft_numerics::switches::RunConfig;
use sha2::{Digest, Sha256};

// Same logarithmic mesh, Gaussian local potential and three projectors as
// scripts/gen_synthetic_upf.py (the existing Rust counterpart avoids Python I/O).
// Add analytic atomic charge and s/p starting orbitals; no UPF is committed.
fn synthetic() -> String {
    let mut pp =
        common::gaussian_with(4., 1.0, true, 161).replace("element=\"Xx\"", "element=\"Si\"");
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
        &array(&|r| 16. / std::f64::consts::PI.sqrt() * r * r * (-r.powi(2)).exp()),
    );
    pp.replace("number_of_wfc=\"0\"", "number_of_wfc=\"2\"")
        .replace("</UPF>", &format!(
            "<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1><PP_CHI.2 l=\"1\" label=\"p\" occupation=\"2\">{}</PP_CHI.2></PP_PSWFC></UPF>",
            array(&|r| r * (-0.6 * r * r).exp()),
            array(&|r| r * r * (-0.6 * r * r).exp())))
}

#[test]
fn reference_and_strict_oracle_identity() {
    let expected = match (std::env::consts::OS, std::env::consts::ARCH) {
        ("macos", "aarch64") if !cfg!(feature = "mkl") => [
            (
                "reference",
                0xc012_bb05_409e_658e_u64,
                "c17a0bf53a44e578fe386b965f79f732a76d1f78b8f909702a7da663f25455d1",
            ),
            (
                "strict_oracle",
                0xc012_bb05_409f_01ad_u64,
                "ba397f136e63aa4c77072166e70381fbab4d300f91093679190e2a8708262fa0",
            ),
        ],
        target => {
            eprintln!(
                "R0 identity skipped: no pure-Rust golden for {target:?} (mkl={})",
                cfg!(feature = "mkl")
            );
            return;
        }
    };
    let pp = dft_upf::parse_upf(&synthetic()).unwrap();
    let potentials = [SpeciesPotential {
        label: "Si",
        potential: &pp,
    }];
    let mut actual = Vec::new();
    for (name, _, _) in expected {
        // Isolate the two drivers without changing process environment. STRICT
        // currently loses to REFERENCE_SCF, so explicitly disable the latter.
        let reference = name == "reference";
        let _scope = RunConfig::from_pairs(&[
            ("reference_scf", if reference { "on" } else { "off" }),
            ("reference_strict", if reference { "off" } else { "on" }),
        ])
        .unwrap()
        .enter();
        let config = Config::parse(
            "\
[system]
name = 'R0 synthetic Si2'
lattice_unit = 'bohr'
lattice_vectors = [[0.0,5.13,5.13],[5.13,0.0,5.13],[5.13,5.13,0.0]]
[[system.atoms]]
species = 'Si'
fractional = [0.0,0.0,0.0]
[[system.atoms]]
species = 'Si'
fractional = [0.25,0.25,0.25]
[pseudopotentials]
Si = 'synthetic-in-memory.upf'
[basis]
ecutwfc_ry = 8.0
kpoints = [[0.0,0.0,0.0]]
[xc]
functional = 'lda-pz'
[scf]
threads = 1
max_iterations = 100
start_density = 'atomic'
mixing = 'pulay'
mixing_beta = 0.5
mixing_history = 8
",
        )
        .unwrap();
        assert_eq!(run::reference_requested(&config), reference);
        let options = run::options(&config).unwrap();
        assert_eq!(options.standard.is_some(), reference);
        let cell = config.cell().unwrap();
        let atoms: Vec<_> = config
            .system
            .atoms
            .iter()
            .map(|a| dft_geometry::Atom::new(&a.species, a.fractional))
            .collect();
        let result = if reference {
            let output = refscf::driver::run(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                run::reference_profile(&config).unwrap(),
            )
            .unwrap();
            assert!(output.outcome.converged, "{}", output.exit_reason);
            output.outcome.result
        } else {
            ScfProblem::new(&cell, &atoms, &potentials, options)
                .unwrap()
                .run()
                .unwrap_or_else(|e| panic!("strict oracle: {e}"))
        };
        let mut hash = Sha256::new();
        for value in &result.density {
            hash.update(value.to_le_bytes());
        }
        let actual_hash = format!("{:x}", hash.finalize());
        println!(
            "{name}: E={:.17e} Ha bits={:016x} density_sha256={actual_hash} iterations={}",
            result.total_hartree(),
            result.total_hartree().to_bits(),
            result.history.len()
        );
        actual.push((result.total_hartree().to_bits(), actual_hash));
    }
    for ((name, energy_bits, density_hash), (bits, hash)) in expected.into_iter().zip(actual) {
        assert_eq!(bits, energy_bits, "{name} energy (Ha)");
        assert_eq!(hash, density_hash, "{name} density (electrons/Bohr^3)");
    }
}
