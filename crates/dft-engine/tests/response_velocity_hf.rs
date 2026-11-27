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
use dft_numerics::linalg::{ComplexMatrix, MatrixOp, gemm};
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

/// R2-GAMMA-VEL steps (Bohr^-1) of the matrix check at Gamma; Richardson
/// pairs (h, h/2) of consecutive entries.
const GAMMA_STEPS: [f64; 3] = [2e-3, 1e-3, 5e-4];

/// The cutoff closest to `centre` (0.01 Ha steps) with no |G|^2/2 within
/// `margin` Ha of it, so that k = 0 +- h e_b (h <= `GAMMA_STEPS[0]`) keeps
/// the Miller set of Gamma (|G+h|^2/2 - |G|^2/2 <= |G| h + h^2/2).
fn cutoff_with_gap(cell: &Cell, centre: f64, margin: f64) -> f64 {
    let b = cell.reciprocal_columns_inverse_bohr();
    let n = 24_i32;
    let mut energies = Vec::new();
    for m0 in -n..=n {
        for m1 in -n..=n {
            for m2 in -n..=n {
                let g: [f64; 3] = std::array::from_fn(|c| {
                    b[c][0] * m0 as f64 + b[c][1] * m1 as f64 + b[c][2] * m2 as f64
                });
                energies.push(g.iter().map(|v| v * v).sum::<f64>() / 2.);
            }
        }
    }
    for i in 0..200 {
        let step = if i % 2 == 0 { i / 2 } else { -(i + 1) / 2 };
        let cutoff = centre + 0.01 * step as f64;
        if energies.iter().all(|e| (e - cutoff).abs() > margin) {
            return cutoff;
        }
    }
    panic!("no cutoff with a {margin} Ha gap near {centre} Ha");
}

/// Result of [`gamma_matrix`]: per step, max_b max_ij |FD - p| / max |p|
/// for the production velocity and for the velocity without the k + G = 0
/// gradient rows; Richardson values of the production one.
struct GammaScan {
    plain: Vec<f64>,
    richardson: Vec<f64>,
    without_g0: Vec<f64>,
    g0_share: f64,
}

/// R2-GAMMA-VEL (engine level): at k = 0 exactly, the full velocity matrix
/// `<u_i| dH/dk_b |u_j>` (lowest `case.bands`, i == j and i != j) of the
/// dielectric path against central differences of the PRODUCTION NSCF
/// Hamiltonian at 0 +- h e_b. The NSCF solve at the shifted points keeps the
/// complete basis (bands = N_pw, dense EVD), so `H(k') = V E V^H` exactly on
/// the common Miller set and `<u_i(0)| H(k') |u_j(0)>` needs no truncation.
/// The cutoff is the one nearest 30 Ha with no G shell within reach of the
/// steps (otherwise the Miller set changes between the points).
fn gamma_matrix(case: &Case<'_>) -> GammaScan {
    let g_max = (2. * 31.0_f64).sqrt();
    let h = GAMMA_STEPS[0];
    let cutoff = cutoff_with_gap(&case.cell, 30.0, g_max * h + h * h);
    let options = ScfOptions {
        cutoff_hartree: cutoff,
        bands: Some(12),
        eigen_tolerance: 1e-8,
        functional: Functional::Pbe,
        allow_functional_override: false,
        ..ScfOptions::strict_reference()
    };
    let gamma = KPoint {
        fractional: [0.; 3],
        weight: 1.,
    };
    let problem = ScfProblem::with_kpoints(
        &case.cell,
        &case.atoms,
        &case.potentials,
        options.clone(),
        &[gamma],
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
    let input = || StartInput {
        cell: &case.cell,
        atoms: &case.atoms,
        potentials: &case.potentials,
        options: &options,
    };
    let solve = |points: &[KPoint], bands: usize| {
        nscf_bands_with(
            input(),
            std::slice::from_ref(&density),
            points,
            BandOptions::new(bands, 1e-9, BandSolver::Dense),
            &refscf::metric::Serial,
            &mut NoStore,
        )
        .unwrap()
    };
    let n = solve(&[gamma], case.bands).bases[0].vectors().len();
    let mut points = vec![gamma];
    for row in &a {
        for h in GAMMA_STEPS {
            for s in [h, -h] {
                points.push(KPoint {
                    fractional: std::array::from_fn(|j| s * row[j] / (2. * PI)),
                    weight: 1.,
                });
            }
        }
    }
    let count = points.len() as f64;
    for p in &mut points {
        p.weight = 1. / count;
    }
    let r = solve(&points, n);
    let miller =
        |i: usize| -> Vec<[i32; 3]> { r.bases[i].vectors().iter().map(|g| g.miller).collect() };
    let reference = miller(0);
    for i in 1..points.len() {
        assert_eq!(miller(i), reference, "{}: Miller set changed", case.name);
    }
    let bands = &r.channels[0];
    // Complete eigenvectors of point i in the global row order.
    let complete = |i: usize, cols: usize| {
        let mut m = ComplexMatrix::zeros(n, cols).unwrap();
        for (row, &g) in bands[i].rows.iter().enumerate() {
            for b in 0..cols {
                m[(g, b)] = bands[i].eigenvectors[(row, b)];
            }
        }
        m
    };
    let x = complete(0, case.bands);
    // <x_i| H(k_point) |x_j> = (X^H V) E (X^H V)^H.
    let projected = |i: usize| {
        let v = complete(i, n);
        let xv = gemm(&x, MatrixOp::Adjoint, &v, MatrixOp::Identity).unwrap();
        let mut scaled = xv.clone();
        for (c, e) in bands[i].eigenvalues.iter().enumerate() {
            for z in scaled.column_mut(c) {
                *z *= *e;
            }
        }
        gemm(&scaled, MatrixOp::Identity, &xv, MatrixOp::Adjoint).unwrap()
    };
    let basis = &r.bases[0];
    let momenta: Vec<[f64; 3]> = basis.momenta().iter().map(|g| g.cartesian).collect();
    let grad = dft_physics::stress::projector_k_gradient(basis, volume, &species).unwrap();
    let mut dropped = grad.clone();
    for (ig, g) in basis.momenta().iter().enumerate() {
        if g.norm_squared == 0. {
            for m in dropped.gradients.iter_mut() {
                for c in 0..m.cols() {
                    m[(ig, c)] = Default::default();
                }
            }
        }
    }
    let p = dft_response::velocity::velocity(&momenta, Some(&grad), &x, &x).unwrap();
    let p0 = dft_response::velocity::velocity(&momenta, Some(&dropped), &x, &x).unwrap();
    let max_diff = |u: &ComplexMatrix, w: &ComplexMatrix| {
        u.as_slice()
            .iter()
            .zip(w.as_slice())
            .map(|(a, b)| (a - b).norm())
            .fold(0_f64, f64::max)
    };
    let mut scan = GammaScan {
        plain: vec![0.; GAMMA_STEPS.len()],
        richardson: vec![0.; GAMMA_STEPS.len() - 1],
        without_g0: vec![0.; GAMMA_STEPS.len()],
        g0_share: 0.,
    };
    let per_axis = 2 * GAMMA_STEPS.len();
    for b in 0..3 {
        let scale = p[b]
            .as_slice()
            .iter()
            .map(|z| z.norm())
            .fold(0_f64, f64::max);
        scan.g0_share = scan.g0_share.max(max_diff(&p[b], &p0[b]) / scale);
        let fd: Vec<ComplexMatrix> = GAMMA_STEPS
            .iter()
            .enumerate()
            .map(|(s, h)| {
                let plus = projected(1 + per_axis * b + 2 * s);
                let minus = projected(2 + per_axis * b + 2 * s);
                let mut d = plus.clone();
                for (z, m) in d.as_mut_slice().iter_mut().zip(minus.as_slice()) {
                    *z = (*z - m) / (2. * h);
                }
                d
            })
            .collect();
        for (s, d) in fd.iter().enumerate() {
            scan.plain[s] = scan.plain[s].max(max_diff(d, &p[b]) / scale);
            scan.without_g0[s] = scan.without_g0[s].max(max_diff(d, &p0[b]) / scale);
        }
        for s in 0..GAMMA_STEPS.len() - 1 {
            let mut rich = fd[s + 1].clone();
            for (z, c) in rich.as_mut_slice().iter_mut().zip(fd[s].as_slice()) {
                *z = (*z * 4. - c) / 3.;
            }
            scan.richardson[s] = scan.richardson[s].max(max_diff(&rich, &p[b]) / scale);
        }
    }
    let max_residual = bands
        .iter()
        .flat_map(|b| b.residuals.iter())
        .fold(0_f64, |m, &x| m.max(x));
    println!(
        "{} at Gamma, cutoff {cutoff:.2} Ha: N_pw {n} (complete EVD at {} points), {} bands, \
         max residual {max_residual:.2e} Ha, k+G = 0 term max|dp|/max|p| {:.3e}",
        case.name,
        points.len(),
        case.bands,
        scan.g0_share
    );
    println!("  h (Bohr^-1) | max rel err full p | Richardson (h, h/2) | without k+G = 0 gradient");
    for (s, h) in GAMMA_STEPS.iter().enumerate() {
        let rich = scan
            .richardson
            .get(s)
            .map_or("-".to_string(), |r| format!("{r:.3e}"));
        println!(
            "  {h:9.1e} | {:.3e} | {rich} | {:.3e}",
            scan.plain[s], scan.without_g0[s]
        );
    }
    scan
}

fn dojo_cases<'a>(
    li: &'a dft_upf::Pseudopotential,
    f: &'a dft_upf::Pseudopotential,
    al: &'a dft_upf::Pseudopotential,
    nn: &'a dft_upf::Pseudopotential,
) -> [Case<'a>; 2] {
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
                potential: li,
            },
            SpeciesPotential {
                label: "F",
                potential: f,
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
                potential: al,
            },
            SpeciesPotential {
                label: "N",
                potential: nn,
            },
        ],
        bands: 120,
        ks: vec![[0.25, 0.0, 1. / 3.], [0.13, 0.31, -0.22]],
    };
    [lif, aln]
}

/// R2-GAMMA-VEL: the velocity of the dielectric path at k = 0 equals the k
/// derivative of the production NSCF Hamiltonian there, including the l = 1
/// projector term at k + G = 0 (whose omission the check detects).
#[test]
#[ignore = "external PseudoDojo PBE UPFs, complete dense EVD; DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-engine --features strict-oracle --test response_velocity_hf gamma -- --ignored --nocapture"]
fn dielectric_velocity_at_gamma_is_the_k_derivative_of_the_nscf_hamiltonian() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("DFT_NC_DIR"));
    let load = |s: &str| real(&dir.join(format!("{s}.upf")));
    let (li, f, al, nn) = (load("Li"), load("F"), load("Al"), load("N"));
    for mut case in dojo_cases(&li, &f, &al, &nn) {
        // The lowest bands (valence and the first conduction groups) carry
        // the transitions of the R2 discrepancy; the complete basis is solved.
        case.bands = 24;
        let scan = gamma_matrix(&case);
        let best = |v: &[f64]| v.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(scan.g0_share > 1e-4, "{}: {:e}", case.name, scan.g0_share);
        assert!(best(&scan.richardson) <= 1e-6, "{}", case.name);
        assert!(best(&scan.without_g0) > 10. * best(&scan.richardson));
    }
}

#[test]
#[ignore = "external PseudoDojo PBE UPFs; DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-engine --features strict-oracle --test response_velocity_hf -- --ignored --nocapture"]
fn dielectric_velocity_is_the_k_derivative_of_the_nscf_hamiltonian() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("DFT_NC_DIR"));
    let load = |s: &str| real(&dir.join(format!("{s}.upf")));
    let (li, f, al, nn) = (load("Li"), load("F"), load("Al"), load("N"));
    let [lif, aln] = dojo_cases(&li, &f, &al, &nn);
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
