//! Reference-start correctness and payload evidence (no timing assertions).
//! Extended matrix: DFT_SI_UPF=/absolute/Si.pz-vbc.UPF cargo test --release
//! -p dft-engine --test refstart silicon_matrix -- --ignored --nocapture
#[allow(dead_code)]
#[path = "../../dft-cli/tests/common/mod.rs"]
mod fixture;
use dft_engine::{
    refstart::{self, StartInput},
    *,
};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use dft_numerics::RealMatrix;
use dft_physics::data::*;

fn synthetic(nonlocal: bool, orbitals: usize) -> dft_upf::Pseudopotential {
    let mut text = fixture::gaussian_with(2., 1., nonlocal, 101);
    let pp = dft_upf::parse_upf(&text).unwrap();
    let array = |f: &dyn Fn(f64) -> f64| {
        pp.mesh()
            .r()
            .iter()
            .map(|&r| format!("{:.17e}", f(r)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let a = text.find("<PP_RHOATOM>").unwrap() + "<PP_RHOATOM>".len();
    let b = text.find("</PP_RHOATOM>").unwrap();
    text.replace_range(
        a..b,
        &array(&|r| 8. / std::f64::consts::PI.sqrt() * r * r * (-r * r).exp()),
    );
    text = text.replace(
        "number_of_wfc=\"0\"",
        &format!("number_of_wfc=\"{orbitals}\""),
    );
    let channels = (0..orbitals)
        .map(|i| {
            format!(
                "<PP_CHI.{} l=\"{}\" label=\"test\" occupation=\"0\">{}</PP_CHI.{}>",
                i + 1,
                i % 2,
                array(&|r| r.powi((i % 2 + 1) as i32) * (-0.6 * r * r).exp()),
                i + 1
            )
        })
        .collect::<String>();
    text = text.replace("</UPF>", &format!("<PP_PSWFC>{channels}</PP_PSWFC></UPF>"));
    dft_upf::parse_upf(&text).unwrap()
}
fn problem(
    cell: &Cell,
    atoms: &[Atom],
    pp: &dft_upf::Pseudopotential,
    o: ScfOptions,
) -> ScfProblem {
    ScfProblem::new(
        cell,
        atoms,
        &[SpeciesPotential {
            label: "Xx",
            potential: pp,
        }],
        o,
    )
    .unwrap()
}
fn options(cutoff: f64, bands: usize) -> ScfOptions {
    ScfOptions {
        cutoff_hartree: cutoff,
        bands: Some(bands),
        start_density: StartDensity::Atomic,
        ..ScfOptions::strict_reference()
    }
}
fn full_orbitals<K: Collective>(state: &refstart::StartState, comm: &K) -> ComplexMatrix {
    let n = state.basis.vectors().len();
    let bands = state.orbitals.cols();
    let mut local = vec![0.; n * bands];
    for j in 0..bands {
        for (i, &r) in state.operator.local_rows().iter().enumerate() {
            local[j * n + r] = state.orbitals[(i, j)];
        }
    }
    let mut full = vec![0.; local.len()];
    comm.tree_sum_f64(&local, &mut full).unwrap();
    let mut c = ComplexMatrix::zeros(n, bands).unwrap();
    for j in 0..bands {
        state
            .half
            .unpack(&full[j * n..(j + 1) * n], c.column_mut(j))
            .unwrap();
    }
    c
}
fn compare(
    name: &str,
    cell: &Cell,
    atoms: &[Atom],
    pp: &dft_upf::Pseudopotential,
    o: &ScfOptions,
    ranks: usize,
) {
    let mut old = problem(cell, atoms, pp, o.clone());
    let density = old.starting_density();
    let guess = old
        .atomic_start_vectors(&density)
        .unwrap()
        .unwrap()
        .remove(0);
    let mut fft = FftPlan3d::new(old.grid_dimensions()).unwrap();
    let hartree = dft_physics::hartree::evaluate(&density, cell, &mut fft).unwrap();
    let xc = o.functional.evaluate(&density, cell, &mut fft).unwrap();
    let potential: Vec<_> = old
        .ionic_potential()
        .iter()
        .zip(hartree.potential.iter())
        .zip(xc.potential.iter())
        .map(|((&ion, h), &x)| ion + (h.re + x))
        .collect();
    let basis = dft_basis::PlaneWaveBasis::new(cell, Hartree(o.cutoff_hartree)).unwrap();
    let db = basis.density_basis(cell).unwrap();
    let grid = dft_basis::FftGridSpec::new(&basis, &db, o.grid_override).unwrap();
    let half = dft_basis::HalfSphere::new(&basis, &grid).unwrap();
    let mut packed = RealMatrix::zeros(basis.vectors().len(), guess.cols()).unwrap();
    for j in 0..guess.cols() {
        let n = packed.rows();
        half.pack(
            guess.column(j),
            &mut packed.as_mut_slice()[j * n..(j + 1) * n],
        )
        .unwrap();
    }
    let ewald_energy = old.ewald_terms().energy.0;
    let results = simulate(ranks, |comm: &SimRank| {
        let species = [SpeciesPotential {
            label: "Xx",
            potential: pp,
        }];
        let state = refstart::start(
            StartInput {
                cell,
                atoms,
                potentials: &species,
                options: o,
            },
            comm,
        )
        .unwrap();
        assert_eq!(state.orbitals.rows(), state.operator.local_len());
        assert_eq!(
            state.memory.start_block_bytes,
            8 * state.orbitals.as_slice().len()
        );
        assert_eq!(state.memory.real_projector_bytes % 8, 0);
        let n = state.fields.gather(comm, &state.density).unwrap();
        let v = state.fields.gather(comm, &state.potential).unwrap();
        let c = full_orbitals(&state, comm);
        let drho = n
            .iter()
            .zip(&density)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        let dv = v
            .iter()
            .zip(&potential)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        let charge = n.sum() * cell.volume_bohr3() / n.len() as f64;
        assert!((charge - state.electrons).abs() < 2e-12 * state.electrons);
        for (z, &g) in state
            .density_g
            .iter()
            .zip(state.fields.reciprocal_indices())
        {
            if g == 0 {
                assert_eq!(
                    z.re.to_bits(),
                    (state.electrons / cell.volume_bohr3()).to_bits()
                );
                assert_eq!(z.im, 0.);
            }
        }
        assert!(drho < 2e-12, "{name} P={ranks} density {drho:e}");
        assert!(dv < 2e-11, "{name} P={ranks} potential {dv:e}");
        let b = c.cols();
        let mut orth: f64 = 0.;
        let mut subspace: f64 = 0.;
        for j in 0..b {
            let mut residual = c.column(j).to_vec();
            for k in 0..b {
                let dot: Complex64 = guess
                    .column(k)
                    .iter()
                    .zip(c.column(j))
                    .map(|(a, b)| a.conj() * b)
                    .sum();
                for (r, &x) in residual.iter_mut().zip(guess.column(k)) {
                    *r -= x * dot;
                }
                let dot: Complex64 = c
                    .column(k)
                    .iter()
                    .zip(c.column(j))
                    .map(|(a, b)| a.conj() * b)
                    .sum();
                orth = orth.max((dot - Complex64::new(f64::from(j == k), 0.)).norm());
            }
            subspace = subspace.max(residual.iter().map(|z| z.norm_sqr()).sum::<f64>().sqrt());
        }
        assert!(orth < 2e-12, "orthogonality {orth:e}");
        assert!(subspace < 2e-8, "{name} P={ranks} subspace {subspace:e}");
        assert!((state.ewald.energy.0 - ewald_energy).abs() < 1e-12);
        // Check actual H images too, not only the span (projectors and layouts).
        let local = state.operator.extract_local(&packed).unwrap();
        let mut op = state.operator;
        let mut image = RealMatrix::zeros(local.rows(), b).unwrap();
        op.apply_block(comm, &local, &mut image).unwrap();
        let mut image_full = vec![0.; packed.rows() * b];
        for j in 0..b {
            for (i, &r) in op.local_rows().iter().enumerate() {
                image_full[j * packed.rows() + r] = image[(i, j)];
            }
        }
        let mut total = vec![0.; image_full.len()];
        comm.tree_sum_f64(&image_full, &mut total).unwrap();
        let eigenvalues: Vec<_> = (0..b)
            .map(|j| {
                (0..packed.rows())
                    .map(|i| packed[(i, j)] * total[j * packed.rows() + i])
                    .sum::<f64>()
            })
            .collect();
        let de = state
            .ritz_values_hartree
            .iter()
            .zip(&eigenvalues)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        assert!(de < 2e-11, "{name} Ritz energies {de:e} Ha");
        println!(
            "START case={name} ranks={ranks} rank={} density_max={drho:.6e} e/Bohr^3 potential_max={dv:.6e} Ha orth={orth:.6e} span={subspace:.6e} ritz_max={de:.6e} Ha charge={charge:.15} memory={:?}",
            comm.rank(),
            state.memory
        );
        (n, c, state.memory)
    });
    assert_eq!(
        results.iter().map(|r| r.2.owned_rows).sum::<usize>(),
        basis.vectors().len()
    );
    assert_eq!(
        results
            .iter()
            .map(|r| r.2.owned_grid_coefficients)
            .sum::<usize>(),
        density.len()
    );
    assert_eq!(
        results
            .iter()
            .map(|r| r.2.owned_plane_points)
            .sum::<usize>(),
        density.len()
    );
}

#[test]
fn synthetic_owned_start() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.4, 4.7, 0.], [0.2, 0.3, 5.2]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    for nonlocal in [false, true] {
        let pp = synthetic(nonlocal, 2);
        for p in [1, 2, 3, 8] {
            compare("synthetic", &cell, &atoms, &pp, &options(2.5, 3), p);
        }
    }
}
#[test]
fn completion_and_empty_owners() {
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let pp = synthetic(false, 1);
    for (cutoff, bands) in [(2.5, 4), (0.02, 1)] {
        for p in [1, 2, 3, 8] {
            compare("completion", &cell, &atoms, &pp, &options(cutoff, bands), p);
        }
    }
}
#[test]
fn unsupported_and_dependent_input_fails_collectively() {
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.; 3])];
    let pp = synthetic(false, 3);
    for ranks in [1, 3, 8] {
        simulate(ranks, |comm: &SimRank| {
            let species = [SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }];
            let o = options(2.5, 3);
            assert!(
                refstart::start(
                    StartInput {
                        cell: &cell,
                        atoms: &atoms,
                        potentials: &species,
                        options: &o
                    },
                    comm
                )
                .is_err()
            );
        });
    }
}

#[test]
#[ignore = "real external UPF and sixteen final-cutoff starts; release command in module docs"]
fn silicon_matrix() {
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("DFT_SI_UPF").expect("DFT_SI_UPF required"))
            .unwrap(),
    )
    .unwrap();
    let base = [
        [0., 0., 0.],
        [0., 0.5, 0.5],
        [0.5, 0., 0.5],
        [0.5, 0.5, 0.],
        [0.25, 0.25, 0.25],
        [0.25, 0.75, 0.75],
        [0.75, 0.25, 0.75],
        [0.75, 0.75, 0.25],
    ];
    for (rx, ry, skew) in [(1, 1, false), (2, 1, false), (2, 2, false), (1, 1, true)] {
        let a = 10.261_212_856_717_933;
        let cell = Cell::from_user_rows_bohr([
            [a * rx as f64, 0., 0.],
            [if skew { 0.7 } else { 0. }, a * ry as f64, 0.],
            [0., if skew { 0.3 } else { 0. }, a],
        ])
        .unwrap();
        let mut atoms = Vec::new();
        for x in 0..rx {
            for y in 0..ry {
                for p in base {
                    let i = atoms.len();
                    atoms.push(Atom::new(
                        "Xx",
                        [
                            (p[0] + x as f64) / rx as f64 + 1e-4 * ((i * 7 + 1) as f64).sin(),
                            (p[1] + y as f64) / ry as f64 + 1e-4 * ((i * 11 + 3) as f64).sin(),
                            p[2] + 1e-4 * ((i * 13 + 5) as f64).sin(),
                        ],
                    ));
                }
            }
        }
        let name = format!("si{}{}", atoms.len(), if skew { "-skew" } else { "" });
        for p in [1, 2, 3, 8] {
            compare(
                &name,
                &cell,
                &atoms,
                &pp,
                &options(15., 2 * atoms.len() + 4),
                p,
            );
        }
    }
}

#[test]
#[ignore = "two verified SCF continuations; run release with --ignored --nocapture"]
fn legacy_driver_from_owned_start() {
    use dft_engine::checkpoint::*;
    let pp = synthetic(true, 2);
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let o = ScfOptions {
        cutoff_hartree: 2.5,
        bands: Some(3),
        start_density: StartDensity::Atomic,
        max_iterations: 200,
        ..ScfOptions::default()
    };
    let mut old = problem(&cell, &atoms, &pp, o.clone());
    let n = old.starting_density();
    let c = old.atomic_start_vectors(&n).unwrap().unwrap().remove(0);
    let owned = simulate(3, |comm: &SimRank| {
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let s = refstart::start(
            StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &o,
            },
            comm,
        )
        .unwrap();
        (
            s.fields.gather(comm, &s.density).unwrap(),
            full_orbitals(&s, comm),
        )
    })
    .remove(0);
    let model = ModelIdentity {
        cell_bohr: cell.lattice_columns_bohr(),
        atoms_fractional: atoms.iter().map(|a| a.fractional_position).collect(),
        species: vec!["Xx".into()],
        pseudopotential_sha256: vec!["0".repeat(64)],
        cutoff_hartree: o.cutoff_hartree,
        kpoints: vec![[0., 0., 0., 1.]],
        electrons: 2.,
        xc: "lda-pz".into(),
        spin: "unpolarized".into(),
        environment: "none".into(),
    };
    let execution = ExecutionIdentity {
        build_sha256: "1".repeat(64),
        platform: "same-process synthetic test".into(),
        threads: 1,
    };
    let run = |(n, c): (Array3<f64>, ComplexMatrix)| {
        let mut template = ScfSession::new(
            problem(&cell, &atoms, &pp, o.clone()),
            CompletionMode::Strict,
            ExecutionConfig::default(),
        )
        .unwrap();
        template.initialize().unwrap();
        let mut cp = template
            .checkpoint(CheckpointKind::WarmStart, model.clone(), execution.clone())
            .unwrap();
        cp.electronic_state.density = n.iter().copied().collect();
        cp.electronic_state.orbitals = vec![Orbitals {
            rows: c.rows(),
            cols: c.cols(),
            coefficients: c.as_slice().iter().flat_map(|z| [z.re, z.im]).collect(),
        }];
        let mut s = ScfSession::warm_start(
            problem(&cell, &atoms, &pp, o.clone()),
            &cp,
            &model,
            &execution,
            CompletionMode::Strict,
            ExecutionConfig {
                calc_id: "owned-start-continuation".into(),
                ..ExecutionConfig::default()
            },
        )
        .unwrap();
        while s.step().unwrap().status == StepStatus::Ready {}
        s.into_result().unwrap()
    };
    let a = run((n, c));
    let b = run(owned);
    assert!(a.converged && b.converged);
    let delta = (a.result.total_hartree() - b.result.total_hartree()).abs();
    assert!(delta < 1e-12, "end-to-end {delta:e} Ha");
    println!(
        "END_TO_END legacy={:.17e} Ha owned={:.17e} Ha delta={delta:.6e} Ha steps={}/{} verification={:?}/{:?}",
        a.result.total_hartree(),
        b.result.total_hartree(),
        a.result.history.len(),
        b.result.history.len(),
        a.result.verification,
        b.result.verification
    );
}

#[test]
fn deterministic_seed_and_library_pw() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.4, 4.7, 0.], [0.2, 0.3, 5.2]]).unwrap();
    let atoms = [Atom::new("Xx", [0.13, 0.21, 0.32])];
    let pp = synthetic(true, 2);
    let mut o = options(2.5, 3);
    o.functional = Functional::LdaPw;
    o.allow_functional_override = true;
    for p in [1, 2, 3, 8] {
        compare("synthetic-pw", &cell, &atoms, &pp, &o, p);
    }
    simulate(3, |comm: &SimRank| {
        let species = [SpeciesPotential {
            label: "Xx",
            potential: &pp,
        }];
        let a = refstart::start(
            StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &o,
            },
            comm,
        )
        .unwrap();
        let mut changed = o.clone();
        changed.seed = o.seed.wrapping_add(97);
        let b = refstart::start(
            StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &species,
                options: &changed,
            },
            comm,
        )
        .unwrap();
        assert_eq!(a.orbitals.as_slice(), b.orbitals.as_slice());
        assert_eq!(a.density, b.density);
        assert_eq!(a.density_g, b.density_g);
        assert_eq!(a.potential, b.potential);
        changed.functional = Functional::Pbe;
        assert!(
            refstart::start(
                StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &species,
                    options: &changed
                },
                comm
            )
            .is_err()
        );
    });
}
