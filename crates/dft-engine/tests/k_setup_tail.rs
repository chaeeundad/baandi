//! GK16: setup tail of the reference k path (common start + k start).
//!
//! Mac stage timing (ignored, observation only):
//!   GK16_RANKS=8 GK16_POOLS=on GK16_K=16 DFT_NC_DIR=/abs/pseudodojo/ex \
//!   cargo test --release -p dft-engine --test k_setup_tail timing -- --ignored --nocapture
//! Without DFT_NC_DIR a synthetic two-species pair stands in for Ir/I.
#[allow(dead_code)]
mod common;
use dft_engine::{
    ScfOptions, SpeciesPotential,
    refstart::{self, StartInput},
};
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_numerics::switches::RunConfig;

fn kp(f: [f64; 3], w: f64) -> KPoint {
    KPoint {
        fractional: f,
        weight: w,
    }
}

/// IrI3-shaped cell: 2 + 6 atoms, skewed (monoclinic-like), about 1350 Bohr^3.
fn iri3_atoms() -> Vec<Atom> {
    vec![
        Atom::new("A", [0.0, 0.0, 0.0]),
        Atom::new("A", [0.5, 0.5, 0.0]),
        Atom::new("B", [0.21, 0.02, 0.23]),
        Atom::new("B", [0.79, 0.98, 0.77]),
        Atom::new("B", [0.27, 0.51, 0.26]),
        Atom::new("B", [0.73, 0.49, 0.74]),
        Atom::new("B", [0.71, 0.03, 0.24]),
        Atom::new("B", [0.29, 0.97, 0.76]),
    ]
}

fn iri3_cell() -> Cell {
    Cell::from_user_rows_bohr([[12.6, 0., 0.], [-1.9, 11.4, 0.], [1.3, -0.9, 9.6]]).unwrap()
}

/// Observation only: every setup phase (rank 0 and the max over ranks) for
/// the common start plus the k start. GK16_K (default 16), GK16_ECUT
/// (default 44 Ha), GK16_RANKS (default 1), GK16_POOLS (off/on, default on),
/// GK16_BANDS (default 42), GK16_PAIRS ("key=value,key=value" switch pairs).
#[test]
#[ignore = "Mac timing observation; run in release"]
fn timing() {
    let nk: usize = std::env::var("GK16_K").map_or(16, |v| v.parse().unwrap());
    let ecut: f64 = std::env::var("GK16_ECUT").map_or(44., |v| v.parse().unwrap());
    let ranks: usize = std::env::var("GK16_RANKS").map_or(1, |v| v.parse().unwrap());
    let bands: usize = std::env::var("GK16_BANDS").map_or(42, |v| v.parse().unwrap());
    let pools = std::env::var("GK16_POOLS").unwrap_or_else(|_| "on".into());
    let extra = std::env::var("GK16_PAIRS").unwrap_or_default();
    let (a, b) = match std::env::var_os("DFT_NC_DIR") {
        Some(dir) => {
            let dir = std::path::PathBuf::from(dir);
            (
                dft_upf::from_file(dir.join("Ir.upf")).unwrap(),
                dft_upf::from_file(dir.join("I.upf")).unwrap(),
            )
        }
        None => (
            common::potential_with_atomic_start(9., 1.8, true, 1301, "SLA PZ"),
            common::potential_with_atomic_start(7., 1.6, true, 1201, "SLA PZ"),
        ),
    };
    let real = std::env::var_os("DFT_NC_DIR").is_some();
    let potentials = [
        SpeciesPotential {
            label: "A",
            potential: &a,
        },
        SpeciesPotential {
            label: "B",
            potential: &b,
        },
    ];
    let atoms = iri3_atoms();
    let cell = iri3_cell();
    let options = ScfOptions {
        cutoff_hartree: ecut,
        bands: Some(bands),
        functional: if real {
            dft_engine::Functional::Pbe
        } else {
            dft_engine::Functional::LdaPz
        },
        ..Default::default()
    };
    let generic = [0.13, -0.29, 0.21];
    let points: Vec<_> = (0..nk)
        .map(|i| {
            kp(
                generic.map(|x| x * (i as f64 + 1.) / nk as f64),
                1. / nk as f64,
            )
        })
        .collect();
    let mut pairs: Vec<(String, String)> = vec![("refscf_k_pools".into(), pools.clone())];
    for kv in extra.split(',').filter(|s| !s.is_empty()) {
        let (k, v) = kv.split_once('=').expect("key=value");
        pairs.push((k.into(), v.into()));
    }
    let report = dft_mpi::simulate(ranks, |comm| {
        let refs: Vec<(&str, &str)> = pairs
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
        let config = RunConfig::from_pairs(&refs).unwrap();
        let _scope = config.enter();
        let clock = std::time::Instant::now();
        let state = refstart::start_kpoints(
            StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            },
            &points,
            comm,
        )
        .unwrap();
        let npw = state
            .blocks
            .iter()
            .map(|b| b.basis.vectors().len())
            .max()
            .unwrap_or(0);
        (
            clock.elapsed().as_secs_f64(),
            npw,
            state.grid.dims(),
            state.setup.phases.clone(),
        )
    });
    let (wall, npw, dims, phases) = &report[0];
    eprintln!(
        "GK16_TIMING real={real} ranks={ranks} pools={pools} pairs={extra} k={nk} npw={npw} dims={dims:?} wall0={wall:.3} wallmax={:.3}",
        report.iter().map(|r| r.0).fold(0., f64::max)
    );
    let mut names: Vec<_> = phases.iter().collect();
    names.sort_by(|a, b| b.1.total_cmp(a.1));
    for (name, s) in names {
        let max = report
            .iter()
            .map(|r| r.3.get(name).copied().unwrap_or(0.))
            .fold(0., f64::max);
        if max >= 0.005 {
            eprintln!("GK16_TIMING   {name:<44} rank0 {s:8.4}  max {max:8.4}");
        }
    }
}

// The GK16 field-table A/B tests against the per-shell quadrature arm
// (`RADIAL_TABLE_FIELDS=off`: delivered energies within the strict 1e-10 Ha
// class on cubic, skewed, elongated and shifted-mesh k sets, spin and FD, and
// the real IrI3 pair) were retired with the switch in CL1 (2026-11-27); the
// evidence is in the attempts ledger (GK16) and
// docs/performance/switch-retirement.md. The table accuracy scan below needs
// no switch and stays.

/// Table error scan on real files (no SCF): for each species, the three field
/// functions (atomic charge integral, Coulomb-subtracted local integral, core
/// transform) tabulated exactly as the engine does (step 0.005, 8 points)
/// against the direct quadrature on 4000 off-node q in (0, 20] Bohr^-1.
/// Prints the worst |table - exact| relative to the function's peak |value|.
#[test]
#[ignore = "requires DFT_NC_DIR; release"]
fn fields_table_error_scan_real() {
    use dft_upf::radial::{
        LocalMesh, RADIAL_TABLE_POINTS as POINTS, RADIAL_TABLE_STEP as STEP, RadialEvenRule,
        RadialTransformTable, spherical_bessel, spherical_density_transform_with_rule,
    };
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let q_max = 20.;
    let nodes = RadialTransformTable::nodes(q_max, STEP, POINTS).unwrap();
    let mut worst_all = 0f64;
    for name in [
        "Ir", "I", "Ge", "Te", "B", "Si", "O", "Fe", "Mg", "Sn", "Br", "Au", "Ga", "As",
    ] {
        let Ok(pp) = dft_upf::from_file(dir.join(format!("{name}.upf"))) else {
            continue;
        };
        let mesh = pp.mesh();
        let rho_samples = &pp.raw().rhoatom.values;
        let rho = |q: f64| {
            let work: Vec<f64> = mesh
                .r()
                .iter()
                .zip(rho_samples)
                .map(|(&r, &f)| {
                    if q == 0. {
                        f
                    } else {
                        f * spherical_bessel(0, q * r).unwrap()
                    }
                })
                .collect();
            mesh.integrate(&work).unwrap()
        };
        let local = pp
            .local_transform_with_conventions(1., LocalMesh::Baandi, RadialEvenRule::default())
            .unwrap();
        let core = |q: f64| {
            spherical_density_transform_with_rule(
                mesh,
                pp.core_density().unwrap(),
                q,
                RadialEvenRule::default(),
            )
            .unwrap()
        };
        type Function<'a> = (&'a str, Box<dyn Fn(f64) -> f64 + 'a>);
        let mut functions: Vec<Function<'_>> = vec![
            ("rho", Box::new(rho)),
            ("local", Box::new(|q| local.short_integral_at(q))),
        ];
        if pp.core_density().is_some() {
            functions.push(("core", Box::new(core)));
        }
        for (f, func) in &functions {
            let values: Vec<f64> = (0..nodes).map(|i| func(i as f64 * STEP)).collect();
            let peak = values.iter().fold(0f64, |m, v| m.max(v.abs()));
            let table = RadialTransformTable::from_nodes(0, q_max, STEP, POINTS, values).unwrap();
            let mut worst = 0f64;
            for i in 1..=4000 {
                let q = i as f64 * (q_max / 4000.) * 0.999_731;
                worst = worst.max((table.eval(q).unwrap() - func(q)).abs());
            }
            worst_all = worst_all.max(worst / peak);
            println!(
                "GK16_SCAN {name:<3} {f:<6} peak={peak:.4e} worst_abs={worst:.3e} worst_rel={:.3e}",
                worst / peak
            );
        }
    }
    println!("GK16_SCAN worst_rel_all={worst_all:.3e}");
    assert!(worst_all < 1e-12, "{worst_all:e}");
}

#[test]
#[ignore = "Mac memory observation; requires DFT_NC_DIR (B.upf); release"]
fn b12_memory() {
    use dft_engine::{Verification, refscf::k_driver, refscf::profile::Profile};
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    let nk: usize = std::env::var("GK16_K").map_or(32, |v| v.parse().unwrap());
    let ecut: f64 = std::env::var("GK16_ECUT").map_or(30., |v| v.parse().unwrap());
    let ranks: usize = std::env::var("GK16_RANKS").map_or(8, |v| v.parse().unwrap());
    let maps: usize = std::env::var("GK16_MAPS").map_or(3, |v| v.parse().unwrap());
    let pools = std::env::var("GK16_POOLS").unwrap_or_else(|_| "off".into());
    let extra = std::env::var("GK16_PAIRS").unwrap_or_default();
    let b = dft_upf::from_file(dir.join("B.upf")).unwrap();
    let potentials = [SpeciesPotential {
        label: "B",
        potential: &b,
    }];
    let (x1, z1, x2, z2) = (0.0103, 0.654, 0.2211, 0.6289);
    let mut atoms = Vec::new();
    for (x, z) in [(x1, z1), (x2, z2)] {
        for p in [[x, x, z], [x, z, x], [z, x, x]] {
            atoms.push(Atom::new("B", p));
            atoms.push(Atom::new("B", p.map(|v| -v)));
        }
    }
    let (a, alpha) = (9.57_f64, 58.06_f64.to_radians());
    let (c, s) = (alpha.cos(), alpha.sin());
    let cy = (c - c * c) / s;
    let cz = (1. - c * c - cy * cy).sqrt();
    let cell =
        Cell::from_user_rows_bohr([[a, 0., 0.], [a * c, a * s, 0.], [a * c, a * cy, a * cz]])
            .unwrap();
    let options = ScfOptions {
        cutoff_hartree: ecut,
        bands: Some(22),
        functional: dft_engine::Functional::Pbe,
        max_iterations: maps,
        verification: Verification::Off,
        ..Default::default()
    };
    let generic = [0.13, -0.29, 0.21];
    let points: Vec<_> = (0..nk)
        .map(|i| {
            kp(
                generic.map(|x| x * (i as f64 + 1.) / nk as f64),
                1. / nk as f64,
            )
        })
        .collect();
    let mut pairs: Vec<(String, String)> = vec![("refscf_k_pools".into(), pools.clone())];
    for kv in extra.split(',').filter(|s| !s.is_empty()) {
        let (k, v) = kv.split_once('=').expect("key=value");
        pairs.push((k.into(), v.into()));
    }
    let report = dft_mpi::simulate(ranks, |comm| {
        let refs: Vec<(&str, &str)> = pairs
            .iter()
            .map(|(k, v)| (k.as_str(), v.as_str()))
            .collect();
        let config = RunConfig::from_pairs(&refs).unwrap();
        let _scope = config.enter();
        let clock = std::time::Instant::now();
        let run = k_driver::run_with(
            StartInput {
                cell: &cell,
                atoms: &atoms,
                potentials: &potentials,
                options: &options,
            },
            Profile::default(),
            &points,
            comm,
        )
        .unwrap();
        (
            clock.elapsed().as_secs_f64(),
            run.records.len(),
            run.peak_payload_bytes_per_rank.clone(),
            run.outcome.result.total_hartree(),
            run.setup.phases.clone(),
            run.start_s,
        )
    });
    let (wall, maps, payload, e, phases, start_s) = &report[0];
    let mut names: Vec<_> = phases.iter().collect();
    names.sort_by(|a, b| b.1.total_cmp(a.1));
    eprintln!("GK16_MEM start_s={start_s:.3}");
    for (name, s) in names.into_iter().take(14) {
        eprintln!("GK16_MEM   {name:<44} rank0 {s:8.4}");
    }
    eprintln!(
        "GK16_MEM ranks={ranks} pools={pools} pairs={extra} k={nk} maps={maps} wall={wall:.2} E={e:.12} payload_per_rank_MB={:?}",
        payload.iter().map(|b| b / 1_000_000).collect::<Vec<_>>()
    );
}
