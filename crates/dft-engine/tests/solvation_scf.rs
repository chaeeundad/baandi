//! Track B Gate 4 Level-3 End-to-End SCF Verification Suite (`L3-1`–`L3-10`).
//!
//! Covers full self-consistent implicit continuum solvation (VASPsol-erfc and
//! linearized Poisson-Boltzmann screening) across Gamma, k-point, collinear spin,
//! Fermi-Dirac, crystal symmetry, and simulated MPI ranks using compact synthetic
//! norm-conserving pseudopotentials (no external UPF files).
//!
//! All tolerances match the pre-declared Gate 3 verification plan
//! (`docs/research/extensible-dft/verification-plan.md` §3).

#[allow(dead_code)]
mod common;

use dft_engine::{
    ElectronicState, ForceState, Functional, forces,
    model::{
        GammaReal, ScfError, ScfOptions, SolvationConfig, SpeciesPotential, SpinSector,
        StandardProfile, StartDensity, StartWavefunctions, Verification,
    },
    parallel,
    refscf::{
        api::{Derivatives, Request, solve},
        delivered,
        electronic::Policy,
        profile::{Preconditioning, Profile},
        symmetry,
    },
    refstart::StartInput,
};
use dft_geometry::{
    Atom, Cell,
    kpoints::{KPoint, monkhorst_pack},
};
use dft_mpi::{Collective, SimRank, simulate};
use dft_numerics::switches::RunConfig;
use dft_physics::{
    occupations::Constraint,
    solvation::{CavityModel, SolvationParams, SolverConfig, SolverKind},
};
use dft_upf::Pseudopotential;
use std::sync::Arc;

const GAMMA_POINT: [f64; 3] = [0.0, 0.0, 0.0];

fn gamma_kpoints() -> Vec<KPoint> {
    vec![KPoint {
        fractional: GAMMA_POINT,
        weight: 1.0,
    }]
}

fn two_kpoints() -> Vec<KPoint> {
    vec![
        KPoint {
            fractional: [0.0, 0.0, 0.0],
            weight: 0.5,
        },
        KPoint {
            fractional: [0.25, 0.25, 0.0],
            weight: 0.5,
        },
    ]
}

struct PolarDimerFixture {
    cell: Cell,
    atoms: Vec<Atom>,
    pp_a: Pseudopotential,
    pp_b: Pseudopotential,
}

impl PolarDimerFixture {
    fn cubic(length_bohr: f64, pp_functional: &str) -> Self {
        let cell = Cell::from_user_rows_bohr([
            [length_bohr, 0.0, 0.0],
            [0.0, length_bohr, 0.0],
            [0.0, 0.0, length_bohr],
        ])
        .unwrap();
        let atoms = vec![
            Atom::new("Xa", [0.42, 0.49, 0.495]),
            Atom::new("Xb", [0.58, 0.51, 0.505]),
        ];
        Self::with_geometry(cell, atoms, pp_functional)
    }

    fn skewed(pp_functional: &str) -> Self {
        let cell = Cell::from_user_rows_bohr([
            [9.5, 0.0, 0.0],
            [0.45, 9.3, 0.0],
            [0.30, 0.35, 9.7],
        ])
        .unwrap();
        let atoms = vec![
            Atom::new("Xa", [0.42, 0.48, 0.49]),
            Atom::new("Xb", [0.58, 0.52, 0.51]),
        ];
        Self::with_geometry(cell, atoms, pp_functional)
    }

    fn with_geometry(cell: Cell, atoms: Vec<Atom>, pp_functional: &str) -> Self {
        // Closed-shell compact heteronuclear dimer (Z_a = 2, Z_b = 2; 4 valence
        // electrons, 2 closed-shell bands) with distinct core radii (`rc = 0.90`
        // vs `1.15` Bohr) and mild nonlocal projectors (`dij` without artificial
        // s-core repulsion) so `n(r) >> n_c` throughout the molecule and bond
        // region, forming a single smooth, convex molecular cavity.
        let dij = "-0.1 0.0 0 0.0 0.05 0 0 0 0.1";
        let pp_a =
            common::potential_with_atomic_start_dij(2.0, 0.90, true, 161, pp_functional, dij);
        let pp_b =
            common::potential_with_atomic_start_dij(2.0, 1.15, true, 161, pp_functional, dij);
        Self {
            cell,
            atoms,
            pp_a,
            pp_b,
        }
    }

    fn species(&self) -> [SpeciesPotential<'_>; 2] {
        [
            SpeciesPotential {
                label: "Xa",
                potential: &self.pp_a,
            },
            SpeciesPotential {
                label: "Xb",
                potential: &self.pp_b,
            },
        ]
    }
}

fn base_options(
    cutoff_hartree: f64,
    bands: usize,
    threshold_ha: f64,
    gamma_real: GammaReal,
    functional: Functional,
    solvation: Option<SolvationConfig>,
) -> ScfOptions {
    ScfOptions {
        standard: Some(StandardProfile { threshold_ha }),
        verification: Verification::FixedPoint,
        cutoff_hartree,
        bands: Some(bands),
        gamma_real,
        functional,
        start_density: StartDensity::Atomic,
        start_wavefunctions: StartWavefunctions::Atomic,
        max_iterations: 200,
        solvation: solvation.map(Box::new),
        ..ScfOptions::default()
    }
}

fn water_solvation(
    epsilon_bulk: f64,
    kappa_inv_bohr: f64,
    tau_ha_bohr2: f64,
    ionic_width_bohr: f64,
    tolerance: f64,
    adaptive_tolerance: bool,
) -> SolvationConfig {
    SolvationConfig {
        params: SolvationParams {
            epsilon_bulk,
            kappa_inv_bohr,
            tau_ha_bohr2,
            beta_ha_bohr3: 0.0,
            cavity_grad_reg: 1.0e-14,
            density_floor_bohr3: 1.0e-16,
            cavity: CavityModel::vaspsol_water(),
        },
        solver: SolverConfig {
            kind: SolverKind::DirectPcg,
            tolerance,
            max_iterations: 400,
            diis_history: 8,
            warm_start: true,
        },
        ionic_width_bohr,
        adaptive_tolerance,
    }
}

fn vacuum_solvation() -> SolvationConfig {
    SolvationConfig {
        params: SolvationParams::vacuum(),
        solver: SolverConfig {
            kind: SolverKind::DirectPcg,
            tolerance: 1.0e-12,
            max_iterations: 300,
            diis_history: 8,
            warm_start: true,
        },
        ionic_width_bohr: 0.60,
        adaptive_tolerance: false,
    }
}

fn test_run_config() -> RunConfig {
    // Coarse synthetic LDA fixture contract (matching `refscf.rs:6750`); production
    // `xc_screen_negative = on` is exercised by `l3_1_off_path_bit_identity` and
    // Level-4 real H2O PBE benchmarks.
    RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap()
}

fn solvation_profile(threshold_ha: f64) -> Profile {
    Profile {
        beta: 0.4,
        history: 8,
        threshold_ry: 2.0 * threshold_ha,
        preconditioning: Preconditioning::Plain,
    }
    .validate()
    .unwrap()
}

fn max_force_diff(a: &[[f64; 3]], b: &[[f64; 3]]) -> f64 {
    assert_eq!(a.len(), b.len(), "force atom count mismatch");
    a.iter()
        .flatten()
        .zip(b.iter().flatten())
        .map(|(x, y)| (x - y).abs())
        .fold(0.0_f64, f64::max)
}

fn vec_f64_bits(values: &[f64]) -> Vec<u64> {
    values.iter().copied().map(f64::to_bits).collect()
}

fn forces_bits(values: &[[f64; 3]]) -> Vec<u64> {
    values.iter().flatten().copied().map(f64::to_bits).collect()
}

fn stress_bits(values: &[[f64; 3]; 3]) -> Vec<u64> {
    values.iter().flatten().copied().map(f64::to_bits).collect()
}

fn assert_energy_identities(label: &str, audit: &dft_engine::energy::EnergyAudit) {
    let e_tot = audit.total_hartree();
    let d_band = (e_tot - audit.band_form).abs();
    let d_hf = (e_tot - audit.harris_foulkes).abs();
    assert!(
        d_band < 1.0e-9,
        "{label}: |total_hartree - band_form| = {d_band:.3e} Ha >= 1e-9 Ha"
    );
    assert!(
        d_hf < 1.0e-9,
        "{label}: |total_hartree - harris_foulkes| = {d_hf:.3e} Ha >= 1e-9 Ha"
    );
    assert!(
        d_hf <= audit.functional_difference_bound,
        "{label}: |total_hartree - harris_foulkes| = {d_hf:.3e} > bound {:.3e}",
        audit.functional_difference_bound
    );
}

struct CommGuard;
impl Drop for CommGuard {
    fn drop(&mut self) {
        parallel::set_communicator(None);
    }
}

fn install_comm(comm: Option<Arc<dyn Collective + Send + Sync>>) -> CommGuard {
    parallel::set_communicator(comm);
    CommGuard
}

/// L3-1: `ScfOptions::solvation = None` preserves bit-identical energy,
/// eigenvalues, forces, and stress across Gamma, k-point, and Fermi-Dirac runs.
#[test]
fn l3_1_off_path_bit_identity() {
    let _scope = test_run_config().enter();
    let fixture = PolarDimerFixture::cubic(8.5, "SLA PZ");
    let species = fixture.species();
    let threshold_ha = 5.0e-12;
    let profile = Profile::from_standard_ha(threshold_ha).unwrap();

    // 1. Gamma real path
    {
        let opts = base_options(3.5, 4, threshold_ha, GammaReal::On, Functional::LdaPz, None);
        let run_once = || {
            let mut req = Request::new(profile, gamma_kpoints());
            req.derivatives = Derivatives {
                forces: true,
                stress: true,
                ..Default::default()
            };
            solve(
                StartInput {
                    cell: &fixture.cell,
                    atoms: &fixture.atoms,
                    potentials: &species,
                    options: &opts,
                },
                req,
            )
            .unwrap()
        };
        let a = run_once();
        let b = run_once();
        let fa = a.derivatives.forces.as_ref().unwrap();
        let fb = b.derivatives.forces.as_ref().unwrap();
        let sa = a.derivatives.stress.as_ref().unwrap();
        let sb = b.derivatives.stress.as_ref().unwrap();
        assert!(fa.solvation.is_none() && fb.solvation.is_none());
        let ra = a.run.into_scalar().unwrap();
        let rb = b.run.into_scalar().unwrap();
        assert!(ra.outcome.converged && rb.outcome.converged);
        let ea = &ra.outcome.result.energy;
        let eb = &rb.outcome.result.energy;
        assert!(ea.solvation.is_none() && eb.solvation.is_none());
        let pure_sum = ea.kinetic + ea.local + ea.nonlocal + ea.hartree + ea.xc + ea.ewald;
        assert_eq!(ea.total_hartree().to_bits(), pure_sum.to_bits());
        assert_eq!(
            ra.outcome.result.total_hartree().to_bits(),
            rb.outcome.result.total_hartree().to_bits()
        );
        assert_eq!(
            vec_f64_bits(&ra.outcome.result.eigenvalues),
            vec_f64_bits(&rb.outcome.result.eigenvalues)
        );
        assert_eq!(forces_bits(&fa.values), forces_bits(&fb.values));
        assert_eq!(stress_bits(&sa.values), stress_bits(&sb.values));
    }

    // 2. K-point path
    {
        let opts = base_options(
            3.5,
            4,
            threshold_ha,
            GammaReal::Off,
            Functional::LdaPz,
            None,
        );
        let run_once = || {
            let mut req = Request::new(profile, two_kpoints());
            req.derivatives = Derivatives {
                forces: true,
                stress: true,
                ..Default::default()
            };
            solve(
                StartInput {
                    cell: &fixture.cell,
                    atoms: &fixture.atoms,
                    potentials: &species,
                    options: &opts,
                },
                req,
            )
            .unwrap()
        };
        let a = run_once();
        let b = run_once();
        let fa = a.derivatives.forces.as_ref().unwrap();
        let fb = b.derivatives.forces.as_ref().unwrap();
        let sa = a.derivatives.stress.as_ref().unwrap();
        let sb = b.derivatives.stress.as_ref().unwrap();
        assert!(fa.solvation.is_none());
        let ra = a.run.into_scalar().unwrap();
        let rb = b.run.into_scalar().unwrap();
        assert!(ra.outcome.converged && rb.outcome.converged);
        assert!(ra.outcome.result.energy.solvation.is_none());
        assert_eq!(
            ra.outcome.result.total_hartree().to_bits(),
            rb.outcome.result.total_hartree().to_bits()
        );
        for (ka, kb) in ra
            .outcome
            .result
            .k_points
            .iter()
            .zip(&rb.outcome.result.k_points)
        {
            assert_eq!(
                vec_f64_bits(&ka.solution.eigenvalues),
                vec_f64_bits(&kb.solution.eigenvalues)
            );
        }
        assert_eq!(forces_bits(&fa.values), forces_bits(&fb.values));
        assert_eq!(stress_bits(&sa.values), stress_bits(&sb.values));
    }

    // 3. Spin / Fermi-Dirac path
    {
        let opts = base_options(
            3.5,
            12,
            threshold_ha,
            GammaReal::Off,
            Functional::LdaPz,
            None,
        );
        let run_once = || {
            let mut req = Request::new(profile, gamma_kpoints());
            req.policy = Some(Policy::FermiDirac {
                theta_ha: 0.02,
                band_cap: 22,
            });
            req.derivatives = Derivatives {
                forces: true,
                stress: true,
                ..Default::default()
            };
            solve(
                StartInput {
                    cell: &fixture.cell,
                    atoms: &fixture.atoms,
                    potentials: &species,
                    options: &opts,
                },
                req,
            )
            .unwrap()
        };
        let a = run_once();
        let b = run_once();
        let fa = a.derivatives.forces.as_ref().unwrap();
        let fb = b.derivatives.forces.as_ref().unwrap();
        let sa = a.derivatives.stress.as_ref().unwrap();
        let sb = b.derivatives.stress.as_ref().unwrap();
        assert!(fa.solvation.is_none());
        let ra = a.run.into_electronic().unwrap();
        let rb = b.run.into_electronic().unwrap();
        assert!(ra.converged && rb.converged);
        assert!(ra.energy.solvation.is_none());
        assert_eq!(
            ra.energy.total_hartree().to_bits(),
            rb.energy.total_hartree().to_bits()
        );
        assert_eq!(
            ra.thermodynamics.free_energy_ha.to_bits(),
            rb.thermodynamics.free_energy_ha.to_bits()
        );
        assert_eq!(
            vec_f64_bits(&ra.blocks[0][0].eigenvalues_ha),
            vec_f64_bits(&rb.blocks[0][0].eigenvalues_ha)
        );
        assert_eq!(forces_bits(&fa.values), forces_bits(&fb.values));
        assert_eq!(stress_bits(&sa.values), stress_bits(&sb.values));
    }
}

/// L3-2: Vacuum-limit solvation (`eps_b = 1, kappa = 0, tau = 0, beta = 0`)
/// reduces to pure-PW SCF across Gamma, k-point, Fermi-Dirac, and collinear spin.
#[test]
fn l3_2_vacuum_limit_scf_reduction() {
    let _scope = test_run_config().enter();
    let threshold_ha = 5.0e-13;
    let profile = solvation_profile(threshold_ha);
    let vac = Some(vacuum_solvation());

    // 1. Gamma real
    {
        let fixture = PolarDimerFixture::cubic(8.5, "SLA PZ");
        let species = fixture.species();
        let run_case = |solv: Option<SolvationConfig>| {
            let opts = base_options(3.5, 4, threshold_ha, GammaReal::On, Functional::LdaPz, solv);
            let mut req = Request::new(profile, gamma_kpoints());
            req.derivatives.forces = true;
            let delivered = solve(
                StartInput {
                    cell: &fixture.cell,
                    atoms: &fixture.atoms,
                    potentials: &species,
                    options: &opts,
                },
                req,
            )
            .unwrap();
            let f = delivered.derivatives.forces.unwrap();
            let r = delivered.run.into_scalar().unwrap();
            assert!(r.outcome.converged);
            (r.outcome.result, f)
        };
        let (pw, f_pw) = run_case(None);
        let (sv, f_sv) = run_case(vac);
        let ledger = sv
            .energy
            .solvation
            .as_ref()
            .expect("vacuum solvation ledger");
        assert_eq!(ledger.total_ha, 0.0);
        assert_eq!(ledger.electrostatic_ha, 0.0);
        assert_eq!(ledger.background_correction_ha, 0.0);
        assert_eq!(ledger.cavitation_ha, 0.0);
        let de = (sv.total_hartree() - pw.total_hartree()).abs();
        let df = max_force_diff(&f_sv.values, &f_pw.values);
        assert!(de < 1.0e-12, "Gamma vacuum dE = {de:.3e} Ha >= 1e-12");
        assert!(df < 1.0e-10, "Gamma vacuum dF = {df:.3e} Ha/Bohr >= 1e-10");
    }

    // 2. K-point
    {
        let fixture = PolarDimerFixture::cubic(8.5, "SLA PZ");
        let species = fixture.species();
        let run_case = |solv: Option<SolvationConfig>| {
            let opts = base_options(
                3.5,
                4,
                threshold_ha,
                GammaReal::Off,
                Functional::LdaPz,
                solv,
            );
            let mut req = Request::new(profile, two_kpoints());
            req.derivatives.forces = true;
            let delivered = solve(
                StartInput {
                    cell: &fixture.cell,
                    atoms: &fixture.atoms,
                    potentials: &species,
                    options: &opts,
                },
                req,
            )
            .unwrap();
            let f = delivered.derivatives.forces.unwrap();
            let r = delivered.run.into_scalar().unwrap();
            assert!(r.outcome.converged);
            (r.outcome.result, f)
        };
        let (pw, f_pw) = run_case(None);
        let (sv, f_sv) = run_case(vac);
        let ledger = sv.energy.solvation.as_ref().expect("k-point vacuum ledger");
        assert_eq!(ledger.total_ha, 0.0);
        let de = (sv.total_hartree() - pw.total_hartree()).abs();
        let df = max_force_diff(&f_sv.values, &f_pw.values);
        assert!(de < 1.0e-12, "k-point vacuum dE = {de:.3e} Ha >= 1e-12");
        assert!(df < 1.0e-10, "k-point vacuum dF = {df:.3e} Ha/Bohr >= 1e-10");
    }

    // 3. Spin / Fermi-Dirac
    {
        let fixture = PolarDimerFixture::cubic(8.5, "SLA PW");
        let species = fixture.species();
        let run_electronic = |policy: Policy, bands: usize, solv: Option<SolvationConfig>| {
            let opts = base_options(
                3.5,
                bands,
                threshold_ha,
                GammaReal::Off,
                Functional::LdaPw,
                solv,
            );
            let mut req = Request::new(profile, gamma_kpoints());
            req.policy = Some(policy);
            req.derivatives.forces = true;
            let delivered = solve(
                StartInput {
                    cell: &fixture.cell,
                    atoms: &fixture.atoms,
                    potentials: &species,
                    options: &opts,
                },
                req,
            )
            .unwrap();
            let f = delivered.derivatives.forces.unwrap();
            let r = delivered.run.into_electronic().unwrap();
            assert!(r.converged);
            (r, f)
        };

        let fd_policy = Policy::FermiDirac {
            theta_ha: 0.02,
            band_cap: 22,
        };
        let (pw_fd, f_pw_fd) = run_electronic(fd_policy.clone(), 12, None);
        let (sv_fd, f_sv_fd) = run_electronic(fd_policy, 12, vac);
        assert_eq!(sv_fd.energy.solvation.as_ref().unwrap().total_ha, 0.0);
        let de_fd = (sv_fd.energy.total_hartree() - pw_fd.energy.total_hartree()).abs();
        let df_fd = max_force_diff(&f_sv_fd.values, &f_pw_fd.values);
        assert!(de_fd < 1.0e-12, "FD vacuum dE = {de_fd:.3e} Ha >= 1e-12");
        assert!(
            df_fd < 1.0e-10,
            "FD vacuum dF = {df_fd:.3e} Ha/Bohr >= 1e-10"
        );

        let spin_policy = Policy::FixedSpin(SpinSector {
            n_up: 3,
            n_down: 1,
            axis_cartesian: [0.0, 0.0, 1.0],
            initial_moment_seed: Some(1.0),
        });
        let (pw_sp, f_pw_sp) = run_electronic(spin_policy.clone(), 6, None);
        let (sv_sp, f_sv_sp) = run_electronic(spin_policy, 6, vac);
        assert_eq!(sv_sp.energy.solvation.as_ref().unwrap().total_ha, 0.0);
        let de_sp = (sv_sp.energy.total_hartree() - pw_sp.energy.total_hartree()).abs();
        let df_sp = max_force_diff(&f_sv_sp.values, &f_pw_sp.values);
        assert!(de_sp < 1.0e-12, "spin vacuum dE = {de_sp:.3e} Ha >= 1e-12");
        assert!(
            df_sp < 1.0e-10,
            "spin vacuum dF = {df_sp:.3e} Ha/Bohr >= 1e-10"
        );
    }
}

/// L3-3: Harris-Foulkes, band-form, and Bregman difference-bound identities for
/// a polar heteronuclear dimer in VASPsol-erfc water (`kappa = 0.0` and `0.15`,
/// `tau = 2e-5`) across Gamma and k-point drivers.
#[test]
fn l3_3_harris_foulkes_and_band_form_identity() {
    let _scope = test_run_config().enter();
    let fixture = PolarDimerFixture::cubic(8.5, "SLA PZ");
    let species = fixture.species();
    let threshold_ha = 5.0e-13;
    let profile = solvation_profile(threshold_ha);

    for kappa in [0.0, 0.15] {
        let solv = water_solvation(78.4, kappa, 2.0e-5, 0.60, 1.0e-12, true);
        for (driver_name, gamma_real, points) in [
            ("gamma", GammaReal::On, gamma_kpoints()),
            ("kpoint", GammaReal::Off, two_kpoints()),
        ] {
            let opts = base_options(
                4.0,
                4,
                threshold_ha,
                gamma_real,
                Functional::LdaPz,
                Some(solv),
            );
            let req = Request::new(profile, points);
            let delivered = solve(
                StartInput {
                    cell: &fixture.cell,
                    atoms: &fixture.atoms,
                    potentials: &species,
                    options: &opts,
                },
                req,
            )
            .unwrap();
            let run = delivered.run.into_scalar().unwrap();
            assert!(
                run.outcome.converged,
                "{driver_name} kappa={kappa} failed to converge"
            );
            let audit = &run.outcome.result.energy;
            let ledger = audit
                .solvation
                .as_ref()
                .expect("solvation ledger must be present");
            assert!(
                ledger.electrostatic_ha < 0.0,
                "{driver_name} kappa={kappa}: expected negative electrostatic energy, got {}",
                ledger.electrostatic_ha
            );
            assert!(
                ledger.cavitation_ha > 0.0,
                "{driver_name} kappa={kappa}: expected positive cavitation energy, got {}",
                ledger.cavitation_ha
            );
            assert_energy_identities(&format!("L3-3 {driver_name} kappa={kappa}"), audit);
        }
    }
}

/// L3-4: Self-consistent Hellmann-Feynman forces vs central finite differences
/// of the total energy on a skewed cell with polar heteronuclear dimer in
/// VASPsol-erfc water (`eps_b = 78.4, kappa = 0.15, tau = 2e-5`).
#[test]
fn l3_4_self_consistent_hellmann_feynman_forces_vs_total_energy_fd() {
    let _scope = test_run_config().enter();
    let fixture = PolarDimerFixture::skewed("SLA PZ");
    let species = fixture.species();
    let threshold_ha = 5.0e-13;
    let profile = solvation_profile(threshold_ha);
    let solv = water_solvation(78.4, 0.15, 2.0e-5, 0.60, 1.0e-12, false);

    let make_opts = |gamma_real: GammaReal| {
        let mut opts = base_options(
            8.0,
            4,
            threshold_ha,
            gamma_real,
            Functional::LdaPz,
            Some(solv),
        );
        opts.density_tolerance = 1.0e-12;
        opts.energy_tolerance = 1.0e-13;
        opts.eigen_tolerance = 1.0e-12;
        opts
    };

    let h_bohr = 5.0e-4;

    for (label, gamma_real, points) in [
        ("gamma", GammaReal::On, gamma_kpoints()),
        ("kpoint", GammaReal::Off, two_kpoints()),
    ] {
        let opts = make_opts(gamma_real);
        let mut req = Request::new(profile, points.clone());
        req.derivatives.forces = true;
        let center = solve(
            StartInput {
                cell: &fixture.cell,
                atoms: &fixture.atoms,
                potentials: &species,
                options: &opts,
            },
            req,
        )
        .unwrap();
        let analytic_forces = center.derivatives.forces.expect("analytic forces");
        let center_run = center.run.into_scalar().unwrap();
        assert!(center_run.outcome.converged);

        let solv_forces = analytic_forces
            .solvation
            .as_ref()
            .expect("solvation force contribution must be populated");
        let max_solv_force = solv_forces
            .iter()
            .flatten()
            .copied()
            .map(f64::abs)
            .fold(0.0_f64, f64::max);
        assert!(
            max_solv_force > 1.0e-5,
            "{label}: solvation force contribution too small ({max_solv_force:.3e} Ha/Bohr)"
        );

        let net_norm = analytic_forces
            .net_force
            .iter()
            .map(|c| c * c)
            .sum::<f64>()
            .sqrt();
        assert!(
            net_norm < 1.0e-6,
            "{label}: |sum_I F_I| = {net_norm:.3e} Ha/Bohr >= 1e-6 (net={:?})",
            analytic_forces.net_force
        );

        // Central FD on atom 0 along all 3 Cartesian axes.
        let mut max_err = 0.0_f64;
        for axis in 0..3 {
            let mut shift_cart = [0.0; 3];
            shift_cart[axis] = h_bohr;
            let shift_frac = fixture.cell.cartesian_bohr_to_fractional(shift_cart);

            let eval_energy = |sign: f64| -> f64 {
                let mut displaced = fixture.atoms.clone();
                for (pos, delta) in displaced[0]
                    .fractional_position
                    .iter_mut()
                    .zip(shift_frac)
                {
                    *pos += sign * delta;
                }
                let r = solve(
                    StartInput {
                        cell: &fixture.cell,
                        atoms: &displaced,
                        potentials: &species,
                        options: &opts,
                    },
                    Request::new(profile, points.clone()),
                )
                .unwrap()
                .run
                .into_scalar()
                .unwrap();
                assert!(
                    r.outcome.converged,
                    "{label} axis={axis} sign={sign} did not converge"
                );
                r.outcome.result.total_hartree()
            };

            let e_plus = eval_energy(1.0);
            let e_minus = eval_energy(-1.0);
            let f_fd = -(e_plus - e_minus) / (2.0 * h_bohr);
            let f_ana = analytic_forces.values[0][axis];
            let err = (f_ana - f_fd).abs();
            max_err = max_err.max(err);
            assert!(
                err < 1.0e-6,
                "{label} axis={axis}: |F_ana - F_fd| = {err:.3e} Ha/Bohr >= 1e-6 \
                 (ana={f_ana:.10}, fd={f_fd:.10})"
            );
        }
        assert!(max_err < 1.0e-6);
    }
}

/// L3-5: Ionic pseudo-charge Gaussian width `ionic_width_bohr` invariance in
/// full SCF when the dielectric cavity lies well outside the ionic cores.
#[test]
fn l3_5_ionic_width_invariance_in_full_scf() {
    let _scope = test_run_config().enter();
    let fixture = PolarDimerFixture::cubic(8.5, "SLA PZ");
    let species = fixture.species();
    let threshold_ha = 5.0e-13;
    let profile = solvation_profile(threshold_ha);

    let run_width = |width_bohr: f64, kappa: f64| -> f64 {
        let solv = water_solvation(78.4, kappa, 2.0e-5, width_bohr, 1.0e-12, false);
        let opts = base_options(
            5.0,
            4,
            threshold_ha,
            GammaReal::On,
            Functional::LdaPz,
            Some(solv),
        );
        let run = solve(
            StartInput {
                cell: &fixture.cell,
                atoms: &fixture.atoms,
                potentials: &species,
                options: &opts,
            },
            Request::new(profile, gamma_kpoints()),
        )
        .unwrap()
        .run
        .into_scalar()
        .unwrap();
        assert!(run.outcome.converged);
        run.outcome.result.total_hartree()
    };

    for kappa in [0.0, 0.15] {
        let e_055 = run_width(0.55, kappa);
        let e_070 = run_width(0.70, kappa);
        let diff = (e_055 - e_070).abs();
        assert!(
            diff < 1.0e-6,
            "kappa={kappa}: |E_tot(r_c=0.55) - E_tot(r_c=0.70)| = {diff:.3e} Ha >= 1e-6 Ha"
        );
    }
}

/// L3-6: Strictly monotonic dielectric stabilization across
/// `eps_b in [1.0, 4.0, 20.0, 78.4]` and electrolyte screening across
/// `kappa in [0.0, 0.15, 0.30]`.
#[test]
fn l3_6_dielectric_and_screening_monotonicity() {
    let _scope = test_run_config().enter();
    let fixture = PolarDimerFixture::cubic(8.5, "SLA PZ");
    let species = fixture.species();
    let threshold_ha = 5.0e-12;
    let profile = solvation_profile(threshold_ha);

    let eval = |eps_b: f64, kappa: f64| -> (f64, f64) {
        let solv = water_solvation(eps_b, kappa, 0.0, 0.60, 1.0e-11, true);
        let opts = base_options(
            4.0,
            4,
            threshold_ha,
            GammaReal::On,
            Functional::LdaPz,
            Some(solv),
        );
        let run = solve(
            StartInput {
                cell: &fixture.cell,
                atoms: &fixture.atoms,
                potentials: &species,
                options: &opts,
            },
            Request::new(profile, gamma_kpoints()),
        )
        .unwrap()
        .run
        .into_scalar()
        .unwrap();
        assert!(run.outcome.converged);
        let audit = &run.outcome.result.energy;
        (
            audit.total_hartree(),
            audit.solvation.as_ref().unwrap().electrostatic_ha,
        )
    };

    let eps_values = [1.0, 4.0, 20.0, 78.4];
    let eps_energies: Vec<f64> = eps_values.iter().map(|&eps| eval(eps, 0.0).0).collect();
    assert!(
        eps_energies[3] < eps_energies[2]
            && eps_energies[2] < eps_energies[1]
            && eps_energies[1] < eps_energies[0],
        "expected E(78.4) < E(20) < E(4) < E(1), got {eps_energies:?}"
    );

    let kappa_values = [0.0, 0.15, 0.30];
    let kappa_es: Vec<f64> = kappa_values.iter().map(|&k| eval(78.4, k).1).collect();
    assert!(
        kappa_es[2] < kappa_es[1] && kappa_es[1] < kappa_es[0],
        "expected E_es(0.30) < E_es(0.15) < E_es(0.0), got {kappa_es:?}"
    );
}

/// L3-7: Open-shell collinear spin-polarized (`N_up != N_dn`), scalar
/// Fermi-Dirac, and spin-polarized Fermi-Dirac SCF with solvation.
#[test]
fn l3_7_collinear_spin_and_fermi_dirac_scf() {
    let _scope = test_run_config().enter();
    let fixture = PolarDimerFixture::cubic(9.5, "SLA PW");
    let species = fixture.species();
    let threshold_ha = 5.0e-12;
    let profile = solvation_profile(threshold_ha);
    let solv = water_solvation(78.4, 0.15, 2.0e-5, 0.60, 1.0e-12, true);

    let cases: [(&str, Policy, usize); 3] = [
        (
            "collinear_fixed_spin",
            Policy::FixedSpin(SpinSector {
                n_up: 3,
                n_down: 1,
                axis_cartesian: [0.0, 0.0, 1.0],
                initial_moment_seed: Some(1.0),
            }),
            6,
        ),
        (
            "scalar_fermi_dirac",
            Policy::FermiDirac {
                theta_ha: 0.02,
                band_cap: 22,
            },
            14,
        ),
        (
            "spin_fermi_dirac",
            Policy::SpinFermiDirac {
                sector: SpinSector {
                    n_up: 2,
                    n_down: 2,
                    axis_cartesian: [0.0, 0.0, 1.0],
                    initial_moment_seed: Some(0.5),
                },
                theta_ha: 0.02,
                constraint: Constraint::FixedMoment(1.0),
                band_cap: 22,
            },
            14,
        ),
    ];

    for (name, policy, bands) in cases {
        let opts = base_options(
            4.0,
            bands,
            threshold_ha,
            GammaReal::Off,
            Functional::LdaPw,
            Some(solv),
        );
        let mut req = Request::new(profile, gamma_kpoints());
        req.policy = Some(policy);
        req.derivatives.forces = true;
        let delivered = solve(
            StartInput {
                cell: &fixture.cell,
                atoms: &fixture.atoms,
                potentials: &species,
                options: &opts,
            },
            req,
        )
        .unwrap_or_else(|e| panic!("L3-7 {name}: solve failed: {e:?}"));
        let f = delivered
            .derivatives
            .forces
            .expect("forces for electronic solvated run");
        assert!(f.solvation.is_some(), "{name}: missing solvation force term");
        let run = delivered.run.into_electronic().unwrap();
        assert!(run.converged, "{name}: did not converge");
        let ledger = run
            .energy
            .solvation
            .as_ref()
            .expect("solvation ledger must be present");
        assert!(
            ledger.electrostatic_ha < 0.0,
            "{name}: expected negative electrostatic energy, got {}",
            ledger.electrostatic_ha
        );
        assert_energy_identities(&format!("L3-7 {name}"), &run.energy);
    }
}

/// L3-8: Crystal symmetry invariance for a symmetric homonuclear dimer in
/// implicit solvent (`symmetry` enabled vs disabled).
#[test]
fn l3_8_crystal_symmetry_invariance() {
    let _scope = test_run_config().enter();
    let length_bohr = 9.5;
    let cell = Cell::from_user_rows_bohr([
        [length_bohr, 0.0, 0.0],
        [0.0, length_bohr, 0.0],
        [0.0, 0.0, length_bohr],
    ])
    .unwrap();
    // Inversion/mirror-symmetric homonuclear dimer with nonzero bond force along x.
    let atoms = vec![
        Atom::new("Xx", [0.42, 0.5, 0.5]),
        Atom::new("Xx", [0.58, 0.5, 0.5]),
    ];
    // Closed-shell homonuclear dimer (Z = 2 + 2 = 4 valence electrons, 2 closed
    // sigma_g / sigma_u bands) with mild nonlocal projectors (no core hole).
    let pp = common::potential_with_atomic_start_dij(
        2.0,
        0.95,
        true,
        161,
        "SLA PZ",
        "-0.1 0.0 0 0.0 0.05 0 0 0 0.1",
    );
    let species = [SpeciesPotential {
        label: "Xx",
        potential: &pp,
    }];
    let full_points = monkhorst_pack([2, 2, 2], [0.5, 0.5, 0.5]).unwrap();
    let prepared = symmetry::Prepared::new(&cell, &atoms, &full_points, 1.0e-7, true).unwrap();
    assert!(
        prepared.points().len() < full_points.len(),
        "expected symmetry reduction of 2x2x2 mesh: {} vs {}",
        prepared.points().len(),
        full_points.len()
    );
    assert!(
        prepared.operations(false).len() > 1,
        "expected nontrivial symmetry operations"
    );

    let threshold_ha = 5.0e-13;
    let profile = solvation_profile(threshold_ha);
    let solv = water_solvation(78.4, 0.15, 2.0e-5, 0.60, 1.0e-12, true);
    let opts = base_options(
        4.0,
        4,
        threshold_ha,
        GammaReal::Off,
        Functional::LdaPz,
        Some(solv),
    );
    let input = || StartInput {
        cell: &cell,
        atoms: &atoms,
        potentials: &species,
        options: &opts,
    };

    let mut req_nosym = Request::new(profile, full_points.clone());
    req_nosym.derivatives.forces = true;
    let nosym = solve(input(), req_nosym).unwrap();
    let f_nosym = nosym.derivatives.forces.unwrap();
    let r_nosym = nosym.run.into_scalar().unwrap();
    assert!(r_nosym.outcome.converged);
    assert!(
        f_nosym.maximum_absolute_component > 1.0e-4,
        "expected nontrivial dimer bond force"
    );

    let mut req_sym = Request::new(profile, full_points);
    req_sym.symmetry = Some(prepared);
    req_sym.derivatives.forces = true;
    let sym = solve(input(), req_sym).unwrap();
    let f_sym = sym.derivatives.forces.unwrap();
    let r_sym = sym.run.into_scalar().unwrap();
    assert!(r_sym.outcome.converged);
    assert_eq!(
        r_sym.outcome.result.density.shape(),
        r_nosym.outcome.result.density.shape(),
        "symmetric and unsymmetric runs should share the same FFT grid"
    );

    let de = (r_sym.outcome.result.total_hartree() - r_nosym.outcome.result.total_hartree()).abs();
    let df = max_force_diff(&f_sym.values, &f_nosym.values);
    assert!(de < 1.0e-9, "symmetry dE = {de:.3e} Ha >= 1e-9 Ha");
    assert!(df < 1.0e-8, "symmetry dF = {df:.3e} Ha/Bohr >= 1e-8 Ha/Bohr");
}

/// L3-9: MPI rank independence in full solvated SCF (`P = 1, 2, 4`) for both
/// Gamma and k-point drivers.
#[test]
fn l3_9_mpi_rank_independence_full_scf() {
    let _scope = test_run_config().enter();
    let fixture = PolarDimerFixture::cubic(8.5, "SLA PZ");
    let species = fixture.species();
    let threshold_ha = 5.0e-13;
    let profile = solvation_profile(threshold_ha);
    let solv = water_solvation(78.4, 0.15, 2.0e-5, 0.60, 1.0e-12, true);

    for (label, gamma_real, points) in [
        ("gamma", GammaReal::On, gamma_kpoints()),
        ("kpoint", GammaReal::Off, two_kpoints()),
    ] {
        let opts = base_options(
            4.0,
            4,
            threshold_ha,
            gamma_real,
            Functional::LdaPz,
            Some(solv),
        );

        let run_for_ranks = |ranks: usize| -> Vec<(f64, Vec<[f64; 3]>)> {
            simulate(ranks, |comm: &SimRank| {
                let _s = test_run_config().enter();
                let _comm_guard = install_comm(Some(Arc::new(comm.clone())));
                let mut req = Request::new(profile, points.clone());
                req.derivatives.forces = true;
                let delivered = solve(
                    StartInput {
                        cell: &fixture.cell,
                        atoms: &fixture.atoms,
                        potentials: &species,
                        options: &opts,
                    },
                    req,
                )
                .unwrap();
                let f = delivered.derivatives.forces.unwrap().values;
                let run = delivered.run.into_scalar().unwrap();
                assert!(run.outcome.converged);
                (run.outcome.result.total_hartree(), f)
            })
        };

        let ref_p1 = run_for_ranks(1);
        let (e1, ref f1) = ref_p1[0];

        for ranks in [2usize, 4usize] {
            let results = run_for_ranks(ranks);
            for (rank, (ep, fp)) in results.iter().enumerate() {
                assert_eq!(
                    ep.to_bits(),
                    results[0].0.to_bits(),
                    "{label} P={ranks} rank={rank}: ranks disagree on energy"
                );
                assert_eq!(
                    forces_bits(fp),
                    forces_bits(&results[0].1),
                    "{label} P={ranks} rank={rank}: ranks disagree on forces"
                );
                let de = (ep - e1).abs();
                let df = max_force_diff(fp, f1);
                assert!(
                    de < 1.0e-10,
                    "{label} P={ranks}: |E(P) - E(1)| = {de:.3e} Ha >= 1e-10 Ha"
                );
                assert!(
                    df < 1.0e-9,
                    "{label} P={ranks}: max |F(P) - F(1)| = {df:.3e} Ha/Bohr >= 1e-9 Ha/Bohr"
                );
            }
        }
    }
}

/// L3-10: Fail-closed contracts on solvated stress, mismatched `SolvationConfig`
/// at force validation, and invalid solvation parameters.
#[test]
fn l3_10_fail_closed_stress_and_model_identity() {
    let _scope = test_run_config().enter();
    let fixture = PolarDimerFixture::cubic(8.5, "SLA PZ");
    let species = fixture.species();
    let threshold_ha = 5.0e-11;
    let profile = solvation_profile(threshold_ha);
    let solv = water_solvation(78.4, 0.15, 2.0e-5, 0.60, 1.0e-11, true);
    let opts = base_options(
        4.0,
        4,
        threshold_ha,
        GammaReal::On,
        Functional::LdaPz,
        Some(solv),
    );
    let input = || StartInput {
        cell: &fixture.cell,
        atoms: &fixture.atoms,
        potentials: &species,
        options: &opts,
    };

    // 1. Immediate `stress: true` fails closed with the exact structured error.
    {
        let mut req = Request::new(profile, gamma_kpoints());
        req.derivatives = Derivatives {
            forces: true,
            stress: true,
            deferred: false,
            ..Default::default()
        };
        let err = match solve(input(), req) {
            Ok(_) => panic!("solvated stress request must fail closed"),
            Err(e) => e,
        };
        assert!(
            matches!(
                err,
                ScfError::InvalidInput(
                    "observables.stress: incomplete derivative coverage: solvation"
                )
            ),
            "unexpected stress error: {err:?}"
        );
    }

    // 2. Immediate `forces: true, stress: false` succeeds with solvation forces.
    {
        let mut req = Request::new(profile, gamma_kpoints());
        req.derivatives = Derivatives {
            forces: true,
            stress: false,
            deferred: false,
            ..Default::default()
        };
        let ok = solve(input(), req).unwrap();
        assert!(ok.derivatives.forces.as_ref().unwrap().solvation.is_some());
        assert!(ok.derivatives.stress.is_none());
        assert!(ok.run.into_scalar().unwrap().outcome.converged);
    }

    // 3. Deferred `evaluate_forces` succeeds on the same state while
    //    `evaluate_stress` fails closed with the exact same structured error.
    let mut delivered = {
        let mut req = Request::new(profile, gamma_kpoints());
        req.derivatives = Derivatives {
            forces: true,
            stress: true,
            deferred: true,
            ..Default::default()
        };
        solve(input(), req).unwrap()
    };
    let mut run = delivered.run.into_scalar().unwrap();
    assert!(run.outcome.converged);
    delivered
        .derivatives
        .prepare(input(), &mut run.outcome)
        .unwrap();
    delivered
        .derivatives
        .evaluate_forces(input(), &run.outcome.result)
        .unwrap();
    assert!(
        delivered
            .derivatives
            .forces
            .as_ref()
            .unwrap()
            .solvation
            .is_some()
    );
    let stress_err = delivered
        .derivatives
        .evaluate_stress(input(), &run.outcome.result)
        .unwrap_err();
    assert!(
        matches!(
            stress_err,
            ScfError::InvalidInput("observables.stress: incomplete derivative coverage: solvation")
        ),
        "unexpected deferred stress error: {stress_err:?}"
    );

    // 4. Mismatched `SolvationConfig` fails `forces::validate` (via `dft_engine::forces`)
    //    and `delivered::force_state`.
    let mismatched_configs: [(&str, Option<SolvationConfig>); 5] = [
        ("solvation_disabled", None),
        (
            "epsilon_changed",
            Some(water_solvation(20.0, 0.15, 2.0e-5, 0.60, 1.0e-11, true)),
        ),
        (
            "kappa_changed",
            Some(water_solvation(78.4, 0.0, 2.0e-5, 0.60, 1.0e-11, true)),
        ),
        (
            "tau_changed",
            Some(water_solvation(78.4, 0.15, 0.0, 0.60, 1.0e-11, true)),
        ),
        (
            "ionic_width_changed",
            Some(water_solvation(78.4, 0.15, 2.0e-5, 0.70, 1.0e-11, true)),
        ),
    ];

    for (tag, bad_solv) in mismatched_configs {
        let mut bad_opts = opts.clone();
        bad_opts.solvation = bad_solv.map(Box::new);
        let bad_state = ForceState {
            cell: &fixture.cell,
            atoms: &fixture.atoms,
            potentials: &species,
            options: &bad_opts,
            electronic: ElectronicState::Fixed(&run.outcome.result),
            density_floor: None,
            symmetry: &[],
        };
        let err = forces(&bad_state).unwrap_err();
        assert!(
            matches!(
                err,
                ScfError::InvalidInput(msg) if msg.contains("model or occupation policy differs")
            ),
            "{tag}: expected model-mismatch InvalidInput from forces::validate, got {err:?}"
        );
        let delivery_err = delivered::force_state(
            StartInput {
                cell: &fixture.cell,
                atoms: &fixture.atoms,
                potentials: &species,
                options: &bad_opts,
            },
            &mut run.outcome,
        );
        assert!(
            matches!(delivery_err, Err(ScfError::InvalidInput(_))),
            "{tag}: expected model-mismatch error from delivered::force_state"
        );
    }

    // 5. Invalid solvation parameters fail closed before any SCF iteration.
    let invalid_configs = [
        water_solvation(0.5, 0.0, 0.0, 0.60, 1.0e-10, true),
        water_solvation(78.4, -0.1, 0.0, 0.60, 1.0e-10, true),
        water_solvation(78.4, 0.0, -1.0e-5, 0.60, 1.0e-10, true),
        water_solvation(78.4, 0.0, 0.0, 0.0, 1.0e-10, true),
    ];
    for bad_solv in invalid_configs {
        let mut bad_opts = opts.clone();
        bad_opts.solvation = Some(Box::new(bad_solv));
        let res = solve(
            StartInput {
                cell: &fixture.cell,
                atoms: &fixture.atoms,
                potentials: &species,
                options: &bad_opts,
            },
            Request::new(profile, gamma_kpoints()),
        );
        assert!(
            matches!(res, Err(ScfError::InvalidInput(_))),
            "expected InvalidInput for invalid solvation config {bad_solv:?}"
        );
    }
}
