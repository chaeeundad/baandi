//! GD1 derivative kernels: hoisted shell values (bit identity), the GEMM
//! nonlocal force/stress kernel against the former scalar loops and against
//! independent finite differences, and the radial tables against the
//! quadrature. Child module of `stress_terms` (shares its fixtures).
use super::*;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::{
    hamiltonian::{NonlocalOperator, Species},
    radial_shells::{RadialShellValues, ShellEvaluator, ShellKind},
    stress::{self, ProjectorRadialTables, nonlocal_k_derivatives, projector_transform_slope},
};
use dft_upf::radial::{RADIAL_TABLE_POINTS, RADIAL_TABLE_STEP, RadialTransformTable};

fn orbitals(rows: usize, bands: usize, seed: f64) -> ComplexMatrix {
    let mut m = ComplexMatrix::zeros(rows, bands).unwrap();
    for b in 0..bands {
        for g in 0..rows {
            m[(g, b)] = Complex64::new(
                (g as f64 * 0.71 + b as f64 * 1.3 + seed).sin(),
                (g as f64 * 0.29 - b as f64 * 0.7 + seed).cos(),
            ) / (1. + 0.01 * g as f64);
        }
        let norm = m.column(b).iter().map(|c| c.norm_sqr()).sum::<f64>().sqrt();
        for g in 0..rows {
            m[(g, b)] /= norm;
        }
    }
    m
}

fn tables(pp: &dft_upf::Pseudopotential, q_max: f64) -> ProjectorRadialTables {
    let nodes = RadialTransformTable::nodes(q_max, RADIAL_TABLE_STEP, RADIAL_TABLE_POINTS).unwrap();
    let build = |i: usize, slope: bool| {
        let l = pp.projectors()[i].l();
        let values = (0..nodes)
            .map(|n| {
                let q = n as f64 * RADIAL_TABLE_STEP;
                if slope {
                    projector_transform_slope(pp, i, q).unwrap()
                } else {
                    pp.projector_transform(i, q).unwrap()
                }
            })
            .collect();
        RadialTransformTable::from_nodes(
            if slope { l + 1 } else { l },
            q_max,
            RADIAL_TABLE_STEP,
            RADIAL_TABLE_POINTS,
            values,
        )
        .unwrap()
    };
    let n = pp.projectors().len();
    ProjectorRadialTables {
        value: (0..n).map(|i| build(i, false)).collect(),
        slope: (0..n).map(|i| build(i, true)).collect(),
    }
}

fn max_abs<'a>(v: impl Iterator<Item = &'a f64>) -> f64 {
    v.fold(0., |m: f64, x| m.max(x.abs()))
}

/// Every hoisted shell value equals the former per-shell expression bit for
/// bit, and the cached local/core forces, stress and core density equal the
/// uncached ones bit for bit.
#[test]
fn gd1_shell_values_and_cached_terms_are_bit_identical() {
    let cell = Cell::from_user_rows_bohr([[8., 0.3, 0.], [1., 7., 0.2], [0.5, 1., 9.]]).unwrap();
    let pp = synthetic(true, "PBE");
    let volume = cell.volume_bohr3();
    let db = PlaneWaveBasis::new(&cell, Hartree(3.)).unwrap();
    let positions = [
        cell.fractional_to_cartesian_bohr([0.2, 0.3, 0.4]),
        cell.fractional_to_cartesian_bohr([0.61, 0.17, 0.83]),
    ];
    let species: Vec<_> = positions
        .iter()
        .map(|p| Species {
            potential: &pp,
            positions_bohr: std::slice::from_ref(p),
        })
        .collect();
    let mut evaluator = ShellEvaluator::new(&pp, volume).unwrap();
    let local = pp.local_transform(volume, None).unwrap();
    let mut cache = RadialShellValues::default();
    for g in db.vectors().iter().filter(|g| g.miller != [0; 3]) {
        let q = g.norm_squared.sqrt();
        for kind in [
            ShellKind::LocalValue,
            ShellKind::LocalSlope,
            ShellKind::CoreValue,
            ShellKind::CoreSlope,
        ] {
            let v = evaluator.eval(kind, q).unwrap().unwrap();
            cache.insert(&pp, kind, q, v);
        }
        assert_eq!(
            cache.get(&pp, ShellKind::LocalValue, q).unwrap().to_bits(),
            local.at_g(q).unwrap().to_bits()
        );
        assert_eq!(
            cache.get(&pp, ShellKind::CoreValue, q).unwrap().to_bits(),
            dft_upf::radial::spherical_density_transform(pp.mesh(), pp.core_density().unwrap(), q)
                .unwrap()
                .to_bits()
        );
    }
    let field: Vec<_> = db
        .vectors()
        .iter()
        .map(|g| {
            Complex64::new(
                (-0.3 * g.norm_squared).exp(),
                0.1 * (g.miller[0] as f64).sin(),
            )
        })
        .collect();
    let bits3 = |v: &[[f64; 3]]| v.iter().flatten().map(|x| x.to_bits()).collect::<Vec<_>>();
    for core in [false, true] {
        let a = dft_physics::forces::local_or_core(&db, volume, &species, &field, core).unwrap();
        let b = dft_physics::forces::local_or_core_cached(
            &db,
            volume,
            &species,
            &field,
            core,
            Some(&cache),
        )
        .unwrap();
        assert_eq!(bits3(&a), bits3(&b), "forces core={core}");
        let a = stress::local_or_core(&db, volume, &species, &field, core).unwrap();
        let b = stress::local_or_core_cached(&db, volume, &species, &field, core, Some(&cache))
            .unwrap();
        assert_eq!(bits3(&a), bits3(&b), "stress core={core}");
    }
    let grid = dft_basis::FftGridSpec::new(&db, &db, None).unwrap();
    let mut fft = dft_numerics::FftPlan3d::new(grid.dims()).unwrap();
    let a = dft_physics::density::core_density(&db, &grid, volume, &species, &mut fft)
        .unwrap()
        .unwrap();
    let b = dft_physics::density::core_density_cached(
        &db,
        &grid,
        volume,
        &species,
        &mut fft,
        Some(&cache),
    )
    .unwrap()
    .unwrap();
    assert!(
        a.iter()
            .zip(b.iter())
            .all(|(x, y)| x.to_bits() == y.to_bits())
    );
}

/// The GEMM kernel (quadrature and tabulated radial values) reproduces the
/// former per-G scalar forces and stress on a skewed cell, Gamma and generic
/// k, l = 0..3 projectors with off-diagonal D, fractional and zero
/// occupations; and its tabulated stress and forces agree with independent
/// central finite differences of the nonlocal energy.
#[test]
fn gd1_nonlocal_kernel_matches_former_loops_and_finite_differences() {
    let pp = stress_support::potential(1., 1.5, true, 401);
    let cell = Cell::from_user_rows_bohr([[8., 0.2, 0.], [1., 7., 0.3], [0.5, 0.7, 9.]]).unwrap();
    let volume = cell.volume_bohr3();
    let fractional = [[0.17, 0.28, 0.39], [0.61, 0.07, 0.74]];
    let positions: Vec<_> = fractional
        .iter()
        .map(|f| cell.fractional_to_cartesian_bohr(*f))
        .collect();
    let species: Vec<_> = positions
        .iter()
        .map(|p| Species {
            potential: &pp,
            positions_bohr: std::slice::from_ref(p),
        })
        .collect();
    let cutoff: f64 = 1.5;
    let table = tables(&pp, 1.1 * (2. * cutoff).sqrt() + 0.05);
    let with_tables = [Some(&table), Some(&table)];
    let occupations = [0.9, 0.35, 0., 0.02];
    for k in [[0.; 3], [0.12, 0.04, -0.1]] {
        let basis = PlaneWaveBasis::at_k(&cell, Hartree(cutoff), k).unwrap();
        let psi = orbitals(basis.vectors().len(), occupations.len(), k[0]);
        let operator = NonlocalOperator::new(&basis, volume, &species).unwrap();
        let old_f = dft_physics::forces::nonlocal(&operator, &basis, &psi, &occupations).unwrap();
        let old_s = stress::nonlocal(&basis, volume, &species, &psi, &occupations).unwrap();
        let scale_f = max_abs(old_f.iter().flatten());
        let scale_s = max_abs(old_s.iter().flatten());
        for (name, t) in [
            ("quadrature", &[None, None][..]),
            ("tables", &with_tables[..]),
        ] {
            let (f, _) =
                nonlocal_k_derivatives(&basis, volume, &species, t, &psi, &occupations, false)
                    .unwrap();
            let (_, s) =
                nonlocal_k_derivatives(&basis, volume, &species, t, &psi, &occupations, true)
                    .unwrap();
            let df = max_abs(
                f.iter()
                    .flatten()
                    .zip(old_f.iter().flatten())
                    .map(|(a, b)| a - b)
                    .collect::<Vec<_>>()
                    .iter(),
            );
            let ds = max_abs(
                s.iter()
                    .flatten()
                    .zip(old_s.iter().flatten())
                    .map(|(a, b)| a - b)
                    .collect::<Vec<_>>()
                    .iter(),
            );
            println!(
                "GD1 kernel k={k:?} {name}: max|dF|={df:e} (scale {scale_f:e}) max|dS|={ds:e} (scale {scale_s:e})"
            );
            assert!(df <= 1e-12 * scale_f.max(1e-300) + 1e-16, "forces {name}");
            assert!(ds <= 1e-12 * scale_s.max(1e-300) + 1e-16, "stress {name}");
        }
    }
    // Independent oracle: finite differences of E_nl at fixed coefficients.
    let k = [0.12, 0.04, -0.1];
    let basis = PlaneWaveBasis::at_k(&cell, Hartree(cutoff), k).unwrap();
    let psi = orbitals(basis.vectors().len(), occupations.len(), 0.3);
    let energy = |c: &Cell, frac: &[[f64; 3]; 2]| {
        let b = basis.in_cell(c).unwrap();
        let pos: Vec<_> = frac
            .iter()
            .map(|f| c.fractional_to_cartesian_bohr(*f))
            .collect();
        let sp: Vec<_> = pos
            .iter()
            .map(|p| Species {
                potential: &pp,
                positions_bohr: std::slice::from_ref(p),
            })
            .collect();
        let op = NonlocalOperator::new(&b, c.volume_bohr3(), &sp).unwrap();
        let mut y = vec![Complex64::default(); psi.rows()];
        (0..occupations.len())
            .map(|band| {
                op.apply(psi.column(band), &mut y).unwrap();
                occupations[band]
                    * psi
                        .column(band)
                        .iter()
                        .zip(&y)
                        .map(|(a, b)| (a.conj() * b).re)
                        .sum::<f64>()
            })
            .sum::<f64>()
    };
    let (_, s) = nonlocal_k_derivatives(
        &basis,
        volume,
        &species,
        &with_tables,
        &psi,
        &occupations,
        true,
    )
    .unwrap();
    check_differences(&cell, s, |c| energy(c, &fractional));
    let (f, _) = nonlocal_k_derivatives(
        &basis,
        volume,
        &species,
        &with_tables,
        &psi,
        &occupations,
        false,
    )
    .unwrap();
    let h = 1e-4;
    let mut worst: f64 = 0.;
    for atom in 0..2 {
        for axis in 0..3 {
            let shift = |sign: f64| {
                let mut frac = fractional;
                let mut cart = cell.fractional_to_cartesian_bohr(frac[atom]);
                cart[axis] += sign * h;
                frac[atom] = cell.cartesian_bohr_to_fractional(cart);
                energy(&cell, &frac)
            };
            let fd = -(shift(1.) - shift(-1.)) / (2. * h);
            worst = worst.max((fd - f[atom][axis]).abs());
        }
    }
    println!("GD1 kernel forces vs central differences (h={h} Bohr): max error {worst:e} Ha/Bohr");
    assert!(worst < 1e-7);
}

/// Real-pseudopotential kernel timing and agreement on a Hg-like k block
/// (observation only): DFT_NC_DIR=/path/to/pseudodojo cargo test --release
/// -p dft-physics --test all derivative_kernels::gd1_real -- --ignored --nocapture
#[test]
#[ignore = "real UPF timing; set DFT_NC_DIR and use the release command"]
fn gd1_real_upf_kernel_timing() {
    let dir = std::path::PathBuf::from(std::env::var("DFT_NC_DIR").expect("DFT_NC_DIR"));
    let load = |e: &str| {
        dft_upf::parse_upf(&std::fs::read_to_string(dir.join(format!("{e}.upf"))).unwrap()).unwrap()
    };
    let (ta, zn) = (load("Ta"), load("Zn"));
    for (name, cutoff, bands) in [("small", 20., 24), ("Hg-like", 30., 48)] {
        let cell =
            Cell::from_user_rows_bohr([[9.5, 0.3, 0.], [0.4, 9.1, 0.2], [0.3, 0.5, 10.2]]).unwrap();
        let volume = cell.volume_bohr3();
        let frac = [[0.1, 0.2, 0.3], [0.55, 0.62, 0.71]];
        let pos: Vec<_> = frac
            .iter()
            .map(|f| cell.fractional_to_cartesian_bohr(*f))
            .collect();
        let species = [
            Species {
                potential: &ta,
                positions_bohr: std::slice::from_ref(&pos[0]),
            },
            Species {
                potential: &zn,
                positions_bohr: std::slice::from_ref(&pos[1]),
            },
        ];
        let basis = PlaneWaveBasis::at_k(&cell, Hartree(cutoff), [0.13, -0.21, 0.07]).unwrap();
        let psi = orbitals(basis.vectors().len(), bands, 0.1);
        let occupations: Vec<f64> = (0..bands)
            .map(|b| 2. / (1. + (b as f64 - 30.).exp()))
            .collect();
        let q_max = 1.1 * (2. * cutoff).sqrt() + 0.05;
        let clock = std::time::Instant::now();
        let (tt, tz) = (tables(&ta, q_max), tables(&zn, q_max));
        let build = clock.elapsed().as_secs_f64();
        let t = [Some(&tt), Some(&tz)];
        let clock = std::time::Instant::now();
        let operator = NonlocalOperator::new(&basis, volume, &species).unwrap();
        let old_f = dft_physics::forces::nonlocal(&operator, &basis, &psi, &occupations).unwrap();
        let old_force_s = clock.elapsed().as_secs_f64();
        let clock = std::time::Instant::now();
        let old_s = stress::nonlocal(&basis, volume, &species, &psi, &occupations).unwrap();
        let old_stress_s = clock.elapsed().as_secs_f64();
        let clock = std::time::Instant::now();
        let (f, _) =
            nonlocal_k_derivatives(&basis, volume, &species, &t, &psi, &occupations, false)
                .unwrap();
        let new_force_s = clock.elapsed().as_secs_f64();
        let clock = std::time::Instant::now();
        let (_, s) =
            nonlocal_k_derivatives(&basis, volume, &species, &t, &psi, &occupations, true).unwrap();
        let new_stress_s = clock.elapsed().as_secs_f64();
        let df = f
            .iter()
            .flatten()
            .zip(old_f.iter().flatten())
            .fold(0_f64, |m, (a, b)| m.max((a - b).abs()));
        let ds = s
            .iter()
            .flatten()
            .zip(old_s.iter().flatten())
            .fold(0_f64, |m, (a, b)| m.max((a - b).abs()));
        println!(
            "{name}: npw {} bands {bands}: forces former {old_force_s:.3}s new {new_force_s:.4}s; stress former {old_stress_s:.3}s new {new_stress_s:.4}s; table build (serial) {build:.3}s; max|dF| {df:e} (max|F| {:e}) max|dS| {ds:e} (max|S| {:e})",
            basis.vectors().len(),
            max_abs(old_f.iter().flatten()),
            max_abs(old_s.iter().flatten()),
        );
        assert!(df <= 1e-11 * max_abs(old_f.iter().flatten()) + 1e-15);
        assert!(ds <= 1e-11 * max_abs(old_s.iter().flatten()) + 1e-15);
    }
}
