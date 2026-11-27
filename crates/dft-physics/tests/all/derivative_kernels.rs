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
            let (joint_f, s) =
                nonlocal_k_derivatives(&basis, volume, &species, t, &psi, &occupations, true)
                    .unwrap();
            // GD2: the stress call returns the same forces bit for bit.
            assert_eq!(
                joint_f
                    .iter()
                    .flatten()
                    .map(|v| v.to_bits())
                    .collect::<Vec<_>>(),
                f.iter().flatten().map(|v| v.to_bits()).collect::<Vec<_>>(),
                "{name}"
            );
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

/// GD2: table of `ShellEvaluator::integral` on the uniform nodes, as the
/// engine builds it (`forces::tabulated_shell_values`).
fn shell_table(e: &mut ShellEvaluator<'_>, kind: ShellKind, q_max: f64) -> RadialTransformTable {
    let nodes = RadialTransformTable::nodes(q_max, RADIAL_TABLE_STEP, RADIAL_TABLE_POINTS).unwrap();
    let values = (0..nodes)
        .map(|i| {
            e.integral(kind, i as f64 * RADIAL_TABLE_STEP)
                .unwrap()
                .unwrap()
        })
        .collect();
    let parity = match kind {
        ShellKind::LocalValue | ShellKind::CoreValue => 0,
        _ => 1,
    };
    RadialTransformTable::from_nodes(
        parity,
        q_max,
        RADIAL_TABLE_STEP,
        RADIAL_TABLE_POINTS,
        values,
    )
    .unwrap()
}

/// Worst |table - quadrature| / max |quadrature| per kind over `qs`.
fn shell_table_errors(pp: &dft_upf::Pseudopotential, volume: f64, qs: &[f64]) -> [f64; 4] {
    let q_max = qs.iter().copied().fold(0., f64::max);
    let mut e = ShellEvaluator::new(pp, volume).unwrap();
    let kinds = [
        ShellKind::LocalValue,
        ShellKind::LocalSlope,
        ShellKind::CoreValue,
        ShellKind::CoreSlope,
    ];
    kinds.map(|kind| {
        if pp.core_density().is_none()
            && matches!(kind, ShellKind::CoreValue | ShellKind::CoreSlope)
        {
            return 0.;
        }
        let table = shell_table(&mut e, kind, q_max);
        let (mut worst, mut peak) = (0f64, 0f64);
        for &q in qs {
            let exact = e.eval(kind, q).unwrap().unwrap();
            // `finish` of the exact integral is the exact value bit for bit.
            let integral = e.integral(kind, q).unwrap().unwrap();
            assert_eq!(
                e.finish(kind, q, integral).unwrap().to_bits(),
                exact.to_bits()
            );
            // Error of the interpolated radial integral (the only tabulated
            // part) relative to that integral's own peak.
            worst = worst.max((table.eval(q).unwrap() - integral).abs());
            peak = peak.max(integral.abs());
        }
        worst / peak
    })
}

/// GD2: tabulated local/core shell values and slopes against the per-shell
/// quadrature on a synthetic NLCC potential (skewed cell, every distinct |G|
/// of a 12 Ha density sphere and off-node values).
#[test]
fn gd2_shell_tables_match_quadrature() {
    let cell = Cell::from_user_rows_bohr([[8., 0.3, 0.], [1., 7., 0.2], [0.5, 1., 9.]]).unwrap();
    let pp = synthetic(true, "PBE");
    let volume = cell.volume_bohr3();
    let db = PlaneWaveBasis::new(&cell, Hartree(12.)).unwrap();
    let mut qs: Vec<f64> = db
        .vectors()
        .iter()
        .filter(|g| g.miller != [0; 3])
        .map(|g| g.norm_squared.sqrt())
        .collect();
    qs.sort_by(f64::total_cmp);
    qs.dedup();
    let errors = shell_table_errors(&pp, volume, &qs);
    println!(
        "GD2 synthetic shell tables ({} shells): relative errors local {:e} local slope {:e} core {:e} core slope {:e}",
        qs.len(),
        errors[0],
        errors[1],
        errors[2],
        errors[3]
    );
    assert!(errors.iter().all(|e| *e < 1e-12), "{errors:?}");
}

/// GD2 table error scan on real PseudoDojo files (ignored):
/// DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-physics --test all gd2_shell_tables_real -- --ignored --nocapture
#[test]
#[ignore = "requires DFT_NC_DIR; release"]
fn gd2_shell_tables_real() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let volume = 300.;
    // 4000 off-node q in (0, 20] Bohr^-1.
    let qs: Vec<f64> = (1..=4000)
        .map(|i| 20. * (i as f64 - 0.37) / 4000.)
        .collect();
    let mut worst = [0f64; 4];
    for element in [
        "Sr", "Ta", "Zn", "Ru", "Hg", "K", "B", "Ge", "Te", "Si", "O", "Fe", "Ir", "I", "Au", "Al",
        "Ti", "Ga", "Hf",
    ] {
        let path = dir.join(format!("{element}.upf"));
        if !path.exists() {
            continue;
        }
        let pp = dft_upf::from_file(&path).unwrap();
        let e = shell_table_errors(&pp, volume, &qs);
        println!(
            "GD2 {element}: local {:.2e} local slope {:.2e} core {:.2e} core slope {:.2e}",
            e[0], e[1], e[2], e[3]
        );
        for (w, e) in worst.iter_mut().zip(e) {
            *w = w.max(e);
        }
    }
    println!("GD2 worst: {worst:?}");
    assert!(worst.iter().all(|e| *e < 1e-12), "{worst:?}");
}

/// GD2: `gga::PointSplit` over a simulated communicator (contiguous blocks,
/// one allgatherv with a status value), as the engine's adapter.
struct SimSplit<'a>(&'a dft_mpi::SimRank);
impl dft_physics::gga::PointSplit for SimSplit<'_> {
    fn block(&self, len: usize) -> std::ops::Range<usize> {
        use dft_mpi::Collective;
        let counts = dft_mpi::block_counts(len, self.0.size());
        let first = dft_mpi::block_displacements(&counts)[self.0.rank()];
        first..first + counts[self.0.rank()]
    }
    fn gather(
        &self,
        local: &[f64],
        len: usize,
        width: usize,
        failed: bool,
    ) -> Result<Vec<f64>, dft_physics::PhysicsError> {
        use dft_mpi::Collective;
        let counts = dft_mpi::block_counts(len, self.0.size());
        let mut send = vec![f64::from(u8::from(failed))];
        send.extend_from_slice(local);
        send.resize(1 + counts[self.0.rank()] * width, 0.);
        let counts: Vec<_> = counts.iter().map(|c| 1 + c * width).collect();
        let mut all = vec![0.; counts.iter().sum()];
        self.0.allgatherv_f64(&send, &counts, &mut all).unwrap();
        let mut out = Vec::new();
        let mut at = 0;
        for c in counts {
            if all[at] != 0. {
                return Err(dft_physics::PhysicsError::Xc("failed".into()));
            }
            out.extend_from_slice(&all[at + 1..at + c]);
            at += c;
        }
        Ok(out)
    }
}

/// GD2: `gga::evaluate_for_derivatives` equals `gga::evaluate` (library
/// correlation) and `gga::evaluate_pw92` (original PW92 prefactor) bit for
/// bit, its metric equals `stress_metric` / `stress_metric_pw92`, for 1, 2,
/// 3 and 7 simulated ranks, with screening off and on; with a screened point
/// the metric is withheld (None).
#[test]
fn gd2_pbe_for_derivatives_bit_identical() {
    use dft_numerics::switches::RunConfig;
    use dft_physics::{gga, pbe_pw92::Pw92Prefactor};
    let cell = Cell::from_user_rows_bohr([[7., 0.3, 0.], [0.6, 6.5, 0.2], [0.4, 0.5, 8.]]).unwrap();
    let dims = [15, 16, 18];
    let smooth = |i: usize, j: usize, k: usize| {
        let (x, y, z) = (i as f64 / 15., j as f64 / 16., k as f64 / 18.);
        0.02 + 0.5 * (-8. * ((x - 0.3).powi(2) + (y - 0.6).powi(2) + (z - 0.4).powi(2))).exp()
            + 0.01 * (std::f64::consts::TAU * (x + 2. * y - z)).sin().powi(2)
    };
    let dense = ndarray::Array3::from_shape_fn(dims, |(i, j, k)| smooth(i, j, k));
    let mut sparse = dense.clone();
    sparse[[3, 4, 5]] = 1e-12;
    let bits = |r: &dft_physics::xc::XcGrid| -> Vec<u64> {
        r.potential
            .iter()
            .chain([&r.energy.0, &r.potential_integral.0])
            .map(|v| v.to_bits())
            .collect()
    };
    let mbits = |m: &[[f64; 3]; 3]| m.iter().flatten().map(|v| v.to_bits()).collect::<Vec<_>>();
    for screening in ["off", "on"] {
        let scope = RunConfig::from_pairs(&[("xc_screen_negative", screening)]).unwrap();
        for (name, n) in [("dense", &dense), ("sparse", &sparse)] {
            for prefactor in [None, Some(Pw92Prefactor::Original)] {
                let _s = scope.enter();
                let mut fft = dft_numerics::FftPlan3d::new(dims).unwrap();
                let (reference, metric) = match prefactor {
                    None => (
                        gga::evaluate(n, &cell, &mut fft).unwrap(),
                        gga::stress_metric(n, &cell, &mut fft).unwrap(),
                    ),
                    Some(p) => (
                        gga::evaluate_pw92(p, n, &cell, &mut fft).unwrap(),
                        gga::stress_metric_pw92(p, n, &cell, &mut fft).unwrap(),
                    ),
                };
                let screened = screening == "on" && name == "sparse";
                let check = |r: gga::DerivativeXc, tag: &str| {
                    assert_eq!(bits(&r.grid), bits(&reference), "{tag}");
                    match r.metric {
                        Some(m) => {
                            assert!(!screened, "{tag}");
                            assert_eq!(mbits(&m), mbits(&metric), "{tag}");
                        }
                        None => assert!(screened, "{tag}"),
                    }
                };
                let tag = format!("{name} {prefactor:?} screening={screening}");
                check(
                    gga::evaluate_for_derivatives(prefactor, n, &cell, &mut fft, None, true)
                        .unwrap(),
                    &tag,
                );
                for ranks in [2, 3, 7] {
                    let results = dft_mpi::simulate(ranks, |comm| {
                        let _s = scope.enter();
                        let mut fft = dft_numerics::FftPlan3d::new(dims).unwrap();
                        let split = SimSplit(comm);
                        gga::evaluate_for_derivatives(
                            prefactor,
                            n,
                            &cell,
                            &mut fft,
                            Some(&split),
                            true,
                        )
                        .unwrap()
                    });
                    for (rank, r) in results.into_iter().enumerate() {
                        check(r, &format!("{tag} ranks={ranks} rank={rank}"));
                    }
                }
            }
        }
    }
}
