//! GD2 stage timing (observation only, ignored): the derivative-pass shell
//! values and the XC evaluation of the NLCC force / XC stress on a Ta6Zn7-like
//! cell (hexagonal lattice, V ~ 1425 Bohr^3, 13 atoms, real PseudoDojo Ta and
//! Zn PBE + NLCC, 40 Ha), legacy arm against the GD2 arm, on 1 and
//! GD2_RANKS (default 8) simulated ranks (threads on one host):
//! DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-engine --lib gd2_stage_timing -- --ignored --nocapture
use super::*;
use dft_numerics::switches::RunConfig;

fn hexagonal(a: f64, c: f64) -> Cell {
    Cell::from_user_rows_bohr([
        [a, 0., 0.],
        [-0.5 * a, 0.5 * 3f64.sqrt() * a, 0.],
        [0., 0., c],
    ])
    .unwrap()
}

#[test]
#[ignore = "requires DFT_NC_DIR (Ta.upf, Zn.upf); release; timing observation"]
fn gd2_stage_timing() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let ta = dft_upf::from_file(dir.join("Ta.upf")).unwrap();
    let zn = dft_upf::from_file(dir.join("Zn.upf")).unwrap();
    let ranks: usize = std::env::var("GD2_RANKS")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(8);
    let cutoff: f64 = std::env::var("GD2_CUTOFF_HA")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or(40.);
    let cell = hexagonal(13.3, 9.3);
    let volume = cell.volume_bohr3();
    let fractional: Vec<[f64; 3]> = (0..13)
        .map(|i| {
            let x = i as f64;
            [
                (0.137 * x + 0.05).fract(),
                (0.291 * x + 0.11).fract(),
                (0.413 * x + 0.07).fract(),
            ]
        })
        .collect();
    let positions: Vec<_> = fractional
        .iter()
        .map(|f| cell.fractional_to_cartesian_bohr(*f))
        .collect();
    let species: Vec<_> = positions
        .iter()
        .enumerate()
        .map(|(i, p)| Species {
            potential: if i < 6 { &ta } else { &zn },
            positions_bohr: std::slice::from_ref(p),
        })
        .collect();
    let options = ScfOptions {
        cutoff_hartree: cutoff,
        functional: crate::Functional::Pbe,
        ..ScfOptions::default()
    };
    let db = options.basis(&cell, [0.; 3], 4.).unwrap();
    let dg = FftGridSpec::new(&db, &db, None).unwrap();
    let dims = dg.dims();
    let mut seen = std::collections::HashSet::new();
    let shells = db
        .vectors()
        .iter()
        .filter(|g| g.miller != [0; 3] && seen.insert(g.norm_squared.sqrt().to_bits()))
        .count();
    // A smooth positive valence density with the cell's electron count.
    let electrons = 6. * ta.z_valence() + 7. * zn.z_valence();
    let mut density = Array3::from_shape_fn(dims, |(i, j, k)| {
        let (x, y, z) = (
            i as f64 / dims[0] as f64,
            j as f64 / dims[1] as f64,
            k as f64 / dims[2] as f64,
        );
        1. + 0.6 * (std::f64::consts::TAU * (x + y)).cos() * (std::f64::consts::TAU * z).sin()
            + 0.3 * (2. * std::f64::consts::TAU * x).sin().powi(2)
    });
    let total = density.sum() * volume / density.len() as f64;
    density.mapv_inplace(|v| v * electrons / total);
    println!(
        "GD2 timing: cutoff {cutoff} Ha, grid {dims:?} ({} points), {} density G, {shells} distinct |G|, ranks {ranks}",
        density.len(),
        db.vectors().len()
    );
    let legacy = RunConfig::from_pairs(&[
        ("deriv_shell_tables", "off"),
        ("deriv_xc_split", "off"),
        ("xc_screen_negative", "off"),
    ])
    .unwrap();
    let new = RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
    // One stage pass: shells (with slopes), core density, XC (+ metric).
    let pass = |scope: &RunConfig, legacy_arm: bool| -> [f64; 4] {
        let _s = scope.enter();
        let mut fft = FftPlan3d::new(dims).unwrap();
        let clock = std::time::Instant::now();
        let shells = shell_values(&species, &db, volume, true).unwrap();
        let t_shells = clock.elapsed().as_secs_f64();
        let clock = std::time::Instant::now();
        let core = dft_physics::density::core_density_cached(
            &db,
            &dg,
            volume,
            &species,
            &mut fft,
            shells.as_ref(),
        )
        .unwrap()
        .unwrap();
        let t_core = clock.elapsed().as_secs_f64();
        let mut scratch = Array3::zeros(dims);
        let clock = std::time::Instant::now();
        // Forces: the NLCC potential alone.
        let force_xc = options
            .evaluate_xc_for_derivatives(
                &density,
                Some(&core),
                &mut scratch,
                &cell,
                &mut fft,
                None,
                false,
            )
            .unwrap();
        let t_force = clock.elapsed().as_secs_f64();
        let clock = std::time::Instant::now();
        // Stress: the evaluation and the GGA metric.
        let (_, metric) = options
            .evaluate_xc_for_derivatives(
                &density,
                Some(&core),
                &mut scratch,
                &cell,
                &mut fft,
                None,
                true,
            )
            .unwrap();
        if metric.is_none() {
            let combined = &density + &core;
            dft_physics::gga::stress_metric(&combined, &cell, &mut fft).unwrap();
        }
        let t_stress = clock.elapsed().as_secs_f64();
        assert_eq!(legacy_arm, metric.is_none());
        drop(force_xc);
        [t_shells, t_core, t_force, t_stress]
    };
    for p in [1, ranks] {
        for (name, scope, legacy_arm) in [("legacy", &legacy, true), ("gd2", &new, false)] {
            let times = if p == 1 {
                vec![pass(scope, legacy_arm)]
            } else {
                dft_mpi::simulate(p, |comm| {
                    crate::parallel::set_communicator(Some(std::sync::Arc::new(comm.clone())));
                    let t = pass(scope, legacy_arm);
                    crate::parallel::set_communicator(None);
                    t
                })
            };
            let max = |i: usize| times.iter().map(|t| t[i]).fold(0., f64::max);
            println!(
                "GD2 timing P={p} {name}: shells {:.3} s, core density {:.3} s, force XC {:.3} s, stress XC+metric {:.3} s (max over ranks)",
                max(0),
                max(1),
                max(2),
                max(3)
            );
        }
    }
}
