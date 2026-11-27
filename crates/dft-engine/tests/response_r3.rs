//! R3 RPA correlation-energy driver (`docs/roadmap/rpa-crpa-plan.md`, ticket
//! R3; ADR 0055) on a synthetic in-memory pseudopotential with nonlocal
//! projectors in a skewed one-atom cell and a fixed nonuniform density (not
//! an SCF fixed point and not reference evidence). Internal gates: rank-count
//! independence (np = 1 bit-identical, others within 1e-9 Ha), q-reduction
//! exactness, frequency-grid and direction-set convergence, band-count
//! monotonicity, Hermiticity and positivity of 1 - v chi0, spin channels
//! equal to the unpolarized path for an unpolarized density, checkpoint
//! resume bit-identical.
//!
//! Fast gate (debug): `cargo test -p dft-engine --test response_r3 -- --nocapture`
#[allow(dead_code)]
mod common;
use dft_engine::response::rpa::{
    NoQStore, QContribution, QStore, QSymmetry, RpaOptions, RpaResult, gamma_mesh, rpa_with,
};
use dft_engine::{refstart::StartInput, *};
use dft_geometry::{Atom, Cell};
use dft_physics::data::Array3;

/// Ha -> meV.
const MEV: f64 = 27_211.386_245_988;

struct System<'a> {
    cell: Cell,
    atoms: Vec<Atom>,
    potentials: Vec<SpeciesPotential<'a>>,
    options: ScfOptions,
    density: Array3<f64>,
}
impl System<'_> {
    fn input(&self) -> StartInput<'_> {
        StartInput {
            cell: &self.cell,
            atoms: &self.atoms,
            potentials: &self.potentials,
            options: &self.options,
        }
    }
    fn run(&self, o: &RpaOptions) -> RpaResult {
        rpa_with(
            self.input(),
            std::slice::from_ref(&self.density),
            o,
            &refscf::metric::Serial,
            &mut NoQStore,
        )
        .unwrap()
    }
}

fn with_system(test: impl FnOnce(&System<'_>)) {
    with_functional(Functional::LdaPz, test)
}
fn with_functional(functional: Functional, test: impl FnOnce(&System<'_>)) {
    let p = common::potential_with_atomic_start_dij(
        2.,
        0.8,
        true,
        401,
        if functional == Functional::LdaPw {
            "SLA PW"
        } else {
            "SLA PZ"
        },
        "-1 0.1 0 0.1 -0.5 0 0 0 0.8",
    );
    let atoms = vec![Atom::new("Xx", [0.13, 0.21, 0.37])];
    let potentials = vec![SpeciesPotential {
        label: "Xx",
        potential: &p,
    }];
    let options = ScfOptions {
        cutoff_hartree: 4.0,
        bands: Some(8),
        eigen_tolerance: 1e-8,
        functional,
        ..ScfOptions::strict_reference()
    };
    let cell = Cell::from_user_rows_bohr([[6.0, 0., 0.], [0.7, 5.5, 0.], [0.3, 0.4, 6.5]]).unwrap();
    let points = gamma_mesh([2, 2, 2]);
    let problem = ScfProblem::with_kpoints(
        &cell,
        &atoms,
        &potentials,
        options.clone(),
        &points,
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap();
    let density = problem.starting_density();
    test(&System {
        cell,
        atoms,
        potentials,
        options,
        density,
    });
}

fn options() -> RpaOptions {
    RpaOptions {
        frequency_points: 12,
        band_block: 7,
        ..RpaOptions::new(vec![1.5, 2.5], [2, 2, 2])
    }
}

fn bits(r: &RpaResult) -> Vec<u64> {
    r.energies
        .iter()
        .chain(r.frequency_terms.iter().flatten())
        .map(|v| v.to_bits())
        .collect()
}

#[test]
fn rpa_energy_contracts() {
    with_system(|system| {
        let base = system.run(&options());
        println!(
            "synthetic skewed Xx, 2x2x2: N_G {:?}, bands {:?}, occupied {:?}, q points {} ({}), omega0 {:.4} Ha, transitions {}",
            base.response_sizes,
            base.bands,
            base.occupied,
            base.q_points.len(),
            base.q_reduction,
            base.omega0,
            base.transitions
        );
        for q in &base.q_points {
            println!(
                "q {:?} w {:.4} optical {}: E_c(q) {:?} Ha",
                q.fractional,
                q.weight,
                q.optical,
                q.cutoffs.iter().map(|c| c.energy).collect::<Vec<_>>()
            );
        }
        println!(
            "E_c per cutoff {:?} Ha; min pivot {:.6}; hermiticity {:.2e}; fit {:?}",
            base.energies, base.min_pivot, base.max_hermiticity, base.extrapolation
        );
        assert!(base.energies.iter().all(|e| *e < 0.));
        assert!(base.min_pivot >= 1. - 1e-12, "{}", base.min_pivot);
        assert!(base.max_hermiticity <= 1e-12, "{}", base.max_hermiticity);
        assert!(base.max_residual <= 1e-6);
        // Larger response basis (bands tied to N_G) lowers E_c.
        assert!(base.energies[1] < base.energies[0]);
        assert_eq!(base.q_points.len(), 8); // 2x2x2: every point is its own -q
        assert!(base.extrapolation.is_some());

        // q reduction: none vs time reversal (identical orbits on 2x2x2) and a
        // 3x3x3 mesh where time reversal merges q and -q.
        let full = system.run(&RpaOptions {
            q_symmetry: QSymmetry::Full,
            ..options()
        });
        let d = (full.energies[1] - base.energies[1]).abs();
        println!("q reduction 2x2x2: |none - time reversal| {d:.3e} Ha");
        assert!(d <= 1e-12);

        // Frequency-grid convergence (plan 2.4: doubling until < 0.05 meV/atom).
        let mut last: Option<f64> = None;
        for points in [6, 12, 24] {
            let r = system.run(&RpaOptions {
                frequency_points: points,
                ecut_chi_series: vec![2.5],
                ..options()
            });
            if let Some(prev) = last {
                let diff = (r.energies[0] - prev) * MEV / r.atoms as f64;
                println!(
                    "N_w {points}: E_c {:.12} Ha, change {diff:.3e} meV/atom",
                    r.energies[0]
                );
                if points == 24 {
                    assert!(diff.abs() <= 0.05, "frequency grid not converged");
                }
            }
            last = Some(r.energies[0]);
        }

        // Direction sets at q = 0 (all from one run).
        let q0 = base.q_points.iter().find(|q| q.optical).unwrap();
        let sets = &q0.cutoffs[1].direction_energies;
        let reference = sets.last().unwrap().1;
        for (label, e) in sets {
            println!(
                "q0 directions {label}: E_c(q=0) {e:.12} Ha, vs product {:.3e} meV",
                (e - reference) * MEV
            );
        }
        let lebedev26 = sets.iter().find(|s| s.0 == "lebedev-26").unwrap().1;
        assert!((lebedev26 - reference).abs() * MEV <= 0.01);
        println!(
            "q0 body only (G = 0 dropped): {:.12} Ha",
            q0.cutoffs[1].body_only_energy.unwrap()
        );

        // Band count: more empty bands at a fixed response basis lower E_c.
        let mut previous = 0.;
        for bands in [6, 12, 24] {
            let r = system.run(&RpaOptions {
                bands: Some(bands),
                ecut_chi_series: vec![2.5],
                ..options()
            });
            println!("bands {bands}: E_c {:.12} Ha", r.energies[0]);
            assert!(r.energies[0] < previous);
            previous = r.energies[0];
        }
    });
}

#[test]
fn rpa_rank_spin_and_resume_contracts() {
    with_system(|system| {
        let o = RpaOptions {
            k_mesh: [3, 1, 1],
            ..options()
        };
        let serial = system.run(&o);
        println!(
            "3x1x1: q points {} weights {:?}, E_c {:?} Ha",
            serial.q_points.len(),
            serial.q_points.iter().map(|q| q.weight).collect::<Vec<_>>(),
            serial.energies
        );
        assert_eq!(serial.q_points.len(), 2, "time reversal merges q and -q");
        let full = system.run(&RpaOptions {
            q_symmetry: QSymmetry::Full,
            ..o.clone()
        });
        let d = (full.energies[1] - serial.energies[1]).abs();
        println!("3x1x1: |none - time reversal| {d:.3e} Ha");
        assert!(d <= 1e-12);

        let reference = bits(&serial);
        for ranks in [1, 2, 3] {
            let runs = dft_mpi::simulate(ranks, |comm| {
                rpa_with(
                    system.input(),
                    std::slice::from_ref(&system.density),
                    &o,
                    comm,
                    &mut NoQStore,
                )
                .unwrap()
            });
            let first = bits(&runs[0]);
            for run in &runs {
                assert_eq!(bits(run), first, "ranks disagree at np={ranks}");
            }
            let spread = runs[0]
                .energies
                .iter()
                .zip(&serial.energies)
                .map(|(a, b)| (a - b).abs())
                .fold(0., f64::max);
            println!("np={ranks}: max |E_c - serial| {spread:.3e} Ha");
            if ranks == 1 {
                assert_eq!(first, reference, "np=1 not bit-identical to serial");
            }
            assert!(spread <= 1e-9);
        }

        // Checkpoint per q: a second run resumes every q, bit-identical.
        #[derive(Default)]
        struct Memory(Vec<QContribution>);
        impl QStore<refscf::metric::Serial> for Memory {
            fn load(
                &mut self,
                _: &refscf::metric::Serial,
                index: usize,
            ) -> Result<Option<QContribution>, ScfError> {
                Ok(self.0.iter().find(|q| q.index == index).cloned())
            }
            fn save(
                &mut self,
                _: &refscf::metric::Serial,
                record: &QContribution,
            ) -> Result<(), ScfError> {
                self.0.push(record.clone());
                Ok(())
            }
        }
        let mut store = Memory::default();
        let run = |store: &mut Memory, o: &RpaOptions| {
            rpa_with(
                system.input(),
                std::slice::from_ref(&system.density),
                o,
                &refscf::metric::Serial,
                store,
            )
        };
        let first = run(&mut store, &o).unwrap();
        assert_eq!(first.resumed, 0);
        // Drop one q: it is recomputed, the other resumed.
        store.0.pop();
        let second = run(&mut store, &o).unwrap();
        assert_eq!(second.resumed, 1);
        assert_eq!(bits(&second), bits(&first));
        // A record from a different frequency grid rejects.
        let other = RpaOptions {
            frequency_points: 8,
            ..o.clone()
        };
        assert!(run(&mut store, &other).is_err());
    });
}

#[test]
fn rpa_spin_channels_equal_unpolarized() {
    // LSDA needs PW92; at zero polarization it equals LDA PW92.
    with_functional(Functional::LdaPw, |system| {
        let o = RpaOptions {
            k_mesh: [3, 1, 1],
            ..options()
        };
        let unpolarized = system.run(&o);
        let half = system.density.mapv(|x| 0.5 * x);
        let spin = rpa_with(
            system.input(),
            &[half.clone(), half],
            &o,
            &refscf::metric::Serial,
            &mut NoQStore,
        )
        .unwrap();
        let d = spin
            .energies
            .iter()
            .zip(&unpolarized.energies)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        println!(
            "spin channels (unpolarized density): occupied {:?}, E_c {:?}, |E_c - unpolarized| {d:.3e} Ha",
            spin.occupied, spin.energies
        );
        assert!(d <= 1e-10);
        // Fully polarized fixed moment (2 up, 0 down): an empty channel adds
        // nothing; a different (LSDA) spectrum gives a different finite E_c.
        let zero = system.density.mapv(|_| 0.);
        let polarized = rpa_with(
            system.input(),
            &[system.density.clone(), zero],
            &o,
            &refscf::metric::Serial,
            &mut NoQStore,
        )
        .unwrap();
        println!(
            "fully polarized: occupied {:?}, gap {:.4} Ha, E_c {:?} Ha",
            polarized.occupied,
            polarized.conduction_minimum - polarized.valence_maximum,
            polarized.energies
        );
        assert_eq!(polarized.occupied, vec![2, 0]);
        assert!(polarized.energies.iter().all(|e| *e < 0. && e.is_finite()));
        assert!(polarized.min_pivot >= 1. - 1e-12);
    });
}
