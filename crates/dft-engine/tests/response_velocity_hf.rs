//! R2-VEL-CHECK (engine level): Hellmann-Feynman check of the dielectric
//! velocity against the PRODUCTION fixed-density Hamiltonian. The bands come
//! from `nscf_bands_with` (dense EVD, the `baandi dielectric` NSCF path, k
//! projectors from the GK14 radial tables of `refstart`); the velocity is the
//! one `response::dielectric` forms (`dft_response::velocity::velocity` with
//! `dft_physics::stress::projector_k_gradient`, exact radial quadrature). For
//! every non-degenerate band, `d e_n / d k_b = <u_n| dH_k/dk_b |u_n>`; the
//! left side is a central difference of the NSCF eigenvalues at k +- h e_b
//! (Richardson over h and h/2) on the same Miller set (checked). The density
//! is the superposition of atomic charges (the identity holds for any fixed
//! density; this is not an SCF fixed point and not QE reference evidence).
//!
//! `DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-engine --features strict-oracle \
//!   --test response_velocity_hf -- --ignored --nocapture`
use dft_engine::refscf::nscf::{BandOptions, BandSolver, NoStore, nscf_bands_with};
use dft_engine::{refstart::StartInput, *};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::hamiltonian::Species;
use std::f64::consts::PI;

fn real(path: &std::path::Path) -> dft_upf::Pseudopotential {
    dft_upf::parse_upf(&std::fs::read_to_string(path).unwrap()).unwrap()
}

struct Case<'a> {
    name: &'static str,
    cell: Cell,
    atoms: Vec<Atom>,
    potentials: Vec<SpeciesPotential<'a>>,
    bands: usize,
    ks: Vec<[f64; 3]>,
}

/// Central-difference steps (Bohr^-1); Richardson pairs (h, h/2) of consecutive entries.
const STEPS: [f64; 5] = [8e-4, 4e-4, 2e-4, 1e-4, 5e-5];
/// Bands closer than this to a neighbour (at k) are skipped (Ha).
const SEPARATION: f64 = 5e-3;

/// Per step: max |central FD - p_nn| / max_n |p_nn| and the Richardson value
/// (h, h/2) of the same, over k, direction and non-degenerate bands.
struct Scan {
    plain: Vec<f64>,
    richardson: Vec<f64>,
    checked: usize,
    skipped: usize,
}

fn hellmann_feynman(case: &Case<'_>) -> Scan {
    let options = ScfOptions {
        cutoff_hartree: 30.0,
        bands: Some(12),
        eigen_tolerance: 1e-8,
        functional: Functional::Pbe,
        allow_functional_override: false,
        ..ScfOptions::strict_reference()
    };
    let problem = ScfProblem::with_kpoints(
        &case.cell,
        &case.atoms,
        &case.potentials,
        options.clone(),
        &[KPoint {
            fractional: [0.; 3],
            weight: 1.,
        }],
        OccupationPolicy::Insulator {
            minimum_gap_hartree: 1e-7,
        },
    )
    .unwrap();
    let density = problem.starting_density();
    let a = case.cell.lattice_columns_bohr();
    let volume = case.cell.volume_bohr3();
    let positions: Vec<Vec<[f64; 3]>> = case
        .potentials
        .iter()
        .map(|p| {
            case.atoms
                .iter()
                .filter(|x| x.species == p.label)
                .map(|x| {
                    case.cell
                        .fractional_to_cartesian_bohr(x.fractional_position)
                })
                .collect()
        })
        .collect();
    let species: Vec<Species<'_>> = case
        .potentials
        .iter()
        .zip(&positions)
        .map(|(p, r)| Species {
            potential: p.potential,
            positions_bohr: r,
        })
        .collect();
    let mut scan = Scan {
        plain: vec![0.; STEPS.len()],
        richardson: vec![0.; STEPS.len() - 1],
        checked: 0,
        skipped: 0,
    };
    let per_axis = 2 * STEPS.len();
    for &k in &case.ks {
        // k, then (+h, -h) for every step and Cartesian direction.
        let mut points = vec![k];
        // a[b][j]: Cartesian component b of lattice vector j.
        for row in &a {
            for h in STEPS {
                for s in [h, -h] {
                    points.push(std::array::from_fn(|j| k[j] + s * row[j] / (2. * PI)));
                }
            }
        }
        let count = points.len() as f64;
        let points: Vec<KPoint> = points
            .into_iter()
            .map(|fractional| KPoint {
                fractional,
                weight: 1. / count,
            })
            .collect();
        let input = StartInput {
            cell: &case.cell,
            atoms: &case.atoms,
            potentials: &case.potentials,
            options: &options,
        };
        let r = nscf_bands_with(
            input,
            std::slice::from_ref(&density),
            &points,
            BandOptions::new(case.bands, 1e-9, BandSolver::Dense),
            &refscf::metric::Serial,
            &mut NoStore,
        )
        .unwrap();
        let miller =
            |i: usize| -> Vec<[i32; 3]> { r.bases[i].vectors().iter().map(|g| g.miller).collect() };
        let reference = miller(0);
        for i in 1..points.len() {
            assert_eq!(
                miller(i),
                reference,
                "{}: Miller set changed at k +- h",
                case.name
            );
        }
        let bands = &r.channels[0];
        let basis = &r.bases[0];
        let n = basis.vectors().len();
        let mut full = ComplexMatrix::zeros(n, case.bands).unwrap();
        for (i, &g) in bands[0].rows.iter().enumerate() {
            for b in 0..case.bands {
                full[(g, b)] = bands[0].eigenvectors[(i, b)];
            }
        }
        let momenta: Vec<[f64; 3]> = basis.momenta().iter().map(|g| g.cartesian).collect();
        let grad = dft_physics::stress::projector_k_gradient(basis, volume, &species).unwrap();
        let p = dft_response::velocity::velocity(&momenta, Some(&grad), &full, &full).unwrap();
        let e = &bands[0].eigenvalues;
        let max_residual = bands
            .iter()
            .flat_map(|b| b.residuals.iter())
            .fold(0_f64, |m, &x| m.max(x));
        for (b, pb) in p.iter().enumerate() {
            let scale = (0..case.bands)
                .map(|m| pb[(m, m)].norm())
                .fold(0_f64, f64::max)
                .max(1e-3);
            let ev = |i: usize| &bands[1 + per_axis * b + i].eigenvalues;
            for m in 0..case.bands {
                let lo = if m == 0 {
                    f64::INFINITY
                } else {
                    e[m] - e[m - 1]
                };
                let hi = if m + 1 == case.bands {
                    f64::INFINITY
                } else {
                    e[m + 1] - e[m]
                };
                // The last requested band has no computed upper neighbour:
                // skip it as well (its partner may be degenerate).
                if lo.min(hi) < SEPARATION || m + 1 == case.bands {
                    scan.skipped += 1;
                    continue;
                }
                let pnn = pb[(m, m)];
                assert!(pnn.im.abs() < 1e-10 * scale.max(1.), "p_nn must be real");
                let d: Vec<f64> = STEPS
                    .iter()
                    .enumerate()
                    .map(|(i, h)| (ev(2 * i)[m] - ev(2 * i + 1)[m]) / (2. * h))
                    .collect();
                for (i, di) in d.iter().enumerate() {
                    scan.plain[i] = scan.plain[i].max((di - pnn.re).abs() / scale);
                }
                for i in 0..STEPS.len() - 1 {
                    let r = (4. * d[i + 1] - d[i]) / 3.;
                    scan.richardson[i] = scan.richardson[i].max((r - pnn.re).abs() / scale);
                }
                scan.checked += 1;
            }
        }
        println!(
            "{} k {k:?}: N_pw {n}, max NSCF residual {max_residual:.2e} Ha",
            case.name
        );
    }
    scan
}

#[test]
#[ignore = "external PseudoDojo PBE UPFs; DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-engine --features strict-oracle --test response_velocity_hf -- --ignored --nocapture"]
fn dielectric_velocity_is_the_k_derivative_of_the_nscf_hamiltonian() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("DFT_NC_DIR"));
    let load = |s: &str| real(&dir.join(format!("{s}.upf")));
    let (li, f, al, nn) = (load("Li"), load("F"), load("Al"), load("N"));
    // Geometry of tests/reference/response/r2/{lif,aln}-444.toml, 60 Ry.
    let lif = Case {
        name: "LiF (Dojo PBE), 60 bands",
        cell: Cell::from_user_rows_bohr([
            [0.0, 3.80785, 3.80785],
            [3.80785, 0.0, 3.80785],
            [3.80785, 3.80785, 0.0],
        ])
        .unwrap(),
        atoms: vec![Atom::new("Li", [0.0; 3]), Atom::new("F", [0.5; 3])],
        potentials: vec![
            SpeciesPotential {
                label: "Li",
                potential: &li,
            },
            SpeciesPotential {
                label: "F",
                potential: &f,
            },
        ],
        bands: 60,
        ks: vec![[0.25, 0.0, 0.25], [0.25, -0.25, 0.0], [0.13, 0.31, -0.22]],
    };
    let aln = Case {
        name: "AlN wurtzite (Dojo PBE), 120 bands",
        cell: Cell::from_user_rows_bohr([
            [5.877, 0.0, 0.0],
            [-2.9385, 5.089631298041145, 0.0],
            [0.0, 0.0, 9.409077],
        ])
        .unwrap(),
        atoms: vec![
            Atom::new("Al", [1. / 3., 2. / 3., 0.]),
            Atom::new("Al", [2. / 3., 1. / 3., 0.5]),
            Atom::new("N", [1. / 3., 2. / 3., 0.382]),
            Atom::new("N", [2. / 3., 1. / 3., 0.882]),
        ],
        potentials: vec![
            SpeciesPotential {
                label: "Al",
                potential: &al,
            },
            SpeciesPotential {
                label: "N",
                potential: &nn,
            },
        ],
        bands: 120,
        ks: vec![[0.25, 0.0, 1. / 3.], [0.13, 0.31, -0.22]],
    };
    for case in [lif, aln] {
        let scan = hellmann_feynman(&case);
        println!(
            "{}: {} non-degenerate (band, k, direction) checked, {} skipped \
             (neighbour within {SEPARATION} Ha or last band)",
            case.name, scan.checked, scan.skipped
        );
        println!("  h (Bohr^-1) | max |central FD - p_nn| / max|p_nn| | Richardson (h, h/2)");
        for (i, h) in STEPS.iter().enumerate() {
            let r = scan
                .richardson
                .get(i)
                .map_or("-".to_string(), |r| format!("{r:.3e}"));
            println!("  {h:9.1e} | {:.3e} | {r}", scan.plain[i]);
        }
        let best = scan
            .richardson
            .iter()
            .copied()
            .fold(f64::INFINITY, f64::min);
        assert!(best <= 1e-6, "{}: {best:e}", case.name);
    }
}
