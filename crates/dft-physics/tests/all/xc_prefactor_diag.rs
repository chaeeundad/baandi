//! Test-only PW92/PBE prefactor probe; no production defaults are changed.
//! Sources: PW92 PRB 45, 13244, Eq. (10), Table I; PBE PRL 77, 3865,
//! Eqs. (7), (8), (10), (14). Full-precision beta/mu: ADR 0018.
//! Real runs: build release dft-cli, then set DFT_NC_DIR and run this test
//! in release with --ignored --nocapture. See xc-prefactor-diagnostic.md.
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_numerics::{Fft3d, FftPlan3d};
use dft_physics::{density::core_density, gga, hamiltonian::Species};
use dft_units::Hartree;
use ndarray::{Array1, Array3};
use num_complex::Complex64;
use std::{f64::consts::PI, path::PathBuf, process::Command};

const MODIFIED: f64 = 0.0310907;
const ORIGINAL: f64 = 0.031091;
const CAP_HA: f64 = 1e-12;

// Energy per electron, Ha. The PBE auxiliary coefficient is distinct from PW92 A.
fn pbe(n: f64, sigma: f64, a: f64) -> f64 {
    if n <= gga::DENSITY_THRESHOLD {
        return 0.0;
    }
    let rs = (3.0 / (4.0 * PI * n)).cbrt();
    let ec = -2.0
        * a
        * (1.0 + 0.21370 * rs)
        * (1.0
            / (2.0
                * a
                * (7.5957 * rs.sqrt() + 3.5876 * rs + 1.6382 * rs.powf(1.5) + 0.49294 * rs * rs)))
            .ln_1p();
    let kf = (3.0 * PI * PI * n).cbrt();
    let s2 = sigma / (4.0 * kf * kf * n * n);
    let ex =
        -3.0 * kf / (4.0 * PI) * (1.0 + 0.804 - 0.804 / (1.0 + 0.2195149727645171 * s2 / 0.804));
    let gamma = (1.0 - 2.0_f64.ln()) / (PI * PI);
    let beta = 0.06672455060314922;
    let t2 = sigma / (4.0 * (4.0 * kf / PI) * n * n);
    let aux = (beta / gamma) / (-ec / gamma).exp_m1();
    let at = aux * t2;
    let h = gamma * ((beta / gamma) * t2 * (1.0 + at) / (1.0 + at + at * at)).ln_1p();
    ex + ec + h
}

// Same real spectral operator as gga.rs: zero Nyquist Miller components
// before multiplying by the reciprocal lattice, including for skew cells.
fn sigma(n: &Array3<f64>, cell: &Cell, fft: &mut FftPlan3d) -> Array3<f64> {
    let dims = fft.shape();
    let b = cell.reciprocal_columns_inverse_bohr();
    let mut ng = Array3::zeros(n.raw_dim());
    fft.to_reciprocal(&n.mapv(|x| Complex64::new(x, 0.0)), &mut ng)
        .unwrap();
    let mut sigma = Array3::zeros(n.raw_dim());
    for row in b {
        let mut dg = ng.clone();
        for ((i, j, k), z) in dg.indexed_iter_mut() {
            let index = [i, j, k];
            let g: f64 = (0..3)
                .map(|a| {
                    let m = if dims[a].is_multiple_of(2) && index[a] == dims[a] / 2 {
                        0.0
                    } else if index[a] <= (dims[a] - 1) / 2 {
                        index[a] as f64
                    } else {
                        index[a] as f64 - dims[a] as f64
                    };
                    row[a] * m
                })
                .sum();
            *z *= Complex64::new(0.0, g);
        }
        let mut dr = Array3::zeros(n.raw_dim());
        fft.to_real(&dg, &mut dr).unwrap();
        sigma += &dr.mapv(|z| z.re * z.re);
    }
    sigma
}

fn energy(n: &Array3<f64>, s: &Array3<f64>, cell: &Cell, a: f64) -> f64 {
    cell.volume_bohr3() / n.len() as f64
        * Array1::from_iter(n.iter().zip(s).map(|(&n, &s)| n * pbe(n, s, a))).sum()
}

#[test]
fn synthetic_grid_probe_matches_engine() {
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [1., 6., 0.], [0.5, 1., 7.]]).unwrap();
    for dims in [[8, 10, 12], [9, 7, 5]] {
        let n = Array3::from_shape_fn(dims, |(i, j, k)| {
            0.1 + 0.02 * (2.0 * PI * i as f64 / dims[0] as f64).cos()
                + 0.01 * ((j + k) as f64).sin()
        });
        let mut fft = FftPlan3d::new(dims).unwrap();
        let s = sigma(&n, &cell, &mut fft);
        let error =
            energy(&n, &s, &cell, MODIFIED) - gga::evaluate(&n, &cell, &mut fft).unwrap().energy.0;
        println!("synthetic {dims:?} evaluator_minus_engine={error:.16e} Ha");
        assert!(error.abs() <= CAP_HA, "STOP: invalid probe: {error:e} Ha");
    }
    assert_eq!(pbe(1e-10, 0.0, MODIFIED), 0.0);
}

#[test]
#[ignore = "real PseudoDojo SCFs; requires DFT_NC_DIR and release dft-cli; see validation doc"]
fn materials_prefactor_diagnostic() {
    let Some(nc) = std::env::var_os("DFT_NC_DIR").map(PathBuf::from) else {
        eprintln!("SKIP: set DFT_NC_DIR to the PseudoDojo UPF directory");
        return;
    };
    if ["Si", "Mg", "O"]
        .iter()
        .any(|s| !nc.join(format!("{s}.upf")).is_file())
    {
        eprintln!("SKIP: Si.upf, Mg.upf and O.upf are required in DFT_NC_DIR");
        return;
    }
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .canonicalize()
        .unwrap();
    let work = std::env::var_os("XC_DIAG_WORK")
        .map(PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join(format!("p27-xc-{}", std::process::id())));
    std::fs::create_dir_all(&work).unwrap();
    let cases = std::env::var("XC_DIAG_CASES").unwrap_or_else(|_| "si-20,mgo-40".into());
    let specs: Vec<_> = cases
        .split(',')
        .map(|name| match name {
            "si-20" => (name, 5.13, 20.0, [30; 3], true),
            "si-30" => (name, 5.13, 30.0, [36; 3], true),
            "si-40" => (name, 5.13, 40.0, [45; 3], true),
            "mgo-40" => (name, 3.975, 40.0, [36; 3], false),
            "mgo-50" => (name, 3.975, 50.0, [36; 3], false),
            "mgo-60" => (name, 3.975, 60.0, [40; 3], false),
            _ => panic!("unsupported material case {name}"),
        })
        .collect();
    // Existing runner checks UPF hashes, convergence, fresh certificate and charge.
    // Input TOMLs are preserved except for absolute runtime UPF paths.
    assert!(
        Command::new("python3")
            .args(["-I", "scripts/run_materials_smoke.py", "--cases", &cases])
            .arg("--work-dir")
            .arg(&work)
            .arg("--output")
            .arg(work.join("smoke.json"))
            .current_dir(&root)
            .status()
            .unwrap()
            .success(),
        "SCF runner failed; see {}",
        work.display()
    );
    println!("SCF evidence: {}", work.display());
    for (name, l, cutoff, dims, silicon) in specs {
        let cell = Cell::from_user_rows_bohr([[0., l, l], [l, 0., l], [l, l, 0.]]).unwrap();
        let basis = PlaneWaveBasis::new(&cell, Hartree(cutoff)).unwrap();
        let db = basis.density_basis(&cell).unwrap();
        let grid = FftGridSpec::new(&db, &db, Some(dims)).unwrap();
        let elements = if silicon {
            vec!["Si", "Si"]
        } else {
            vec!["Mg", "O"]
        };
        let potentials: Vec<_> = elements
            .iter()
            .map(|e| {
                dft_upf::parse_upf(&std::fs::read_to_string(nc.join(format!("{e}.upf"))).unwrap())
                    .unwrap()
            })
            .collect();
        let second = if silicon { 0.25 } else { 0.5 };
        let positions = [
            [cell.fractional_to_cartesian_bohr([0.; 3])],
            [cell.fractional_to_cartesian_bohr([second; 3])],
        ];
        let species: Vec<_> = potentials
            .iter()
            .zip(&positions)
            .map(|(p, r)| Species {
                potential: p,
                positions_bohr: r,
            })
            .collect();
        let mut fft = FftPlan3d::new(dims).unwrap();
        let core = core_density(&db, &grid, cell.volume_bohr3(), &species, &mut fft)
            .unwrap()
            .unwrap();
        // Verify the CLI file contract and provenance without a new Rust dependency.
        let metadata = Command::new("python3")
            .args([
                "-I",
                "-c",
                r#"
import hashlib, json, pathlib, sys
p = pathlib.Path(sys.argv[1])
m = json.loads((p / 'density.json').read_text())
r = json.loads((p / 'results.json').read_text())
d = int(sys.argv[2]); l = float(sys.argv[3])
assert m['fft_grid'] == [d, d, d]
assert m['lattice_columns_bohr'] == [[0,l,l],[l,0,l],[l,l,0]]
assert m['dtype'] == 'f64' and m['endianness'] == 'little'
assert m['ordering'] == 'C: (i*n2+j)*n3+k; k fastest'
assert m['units'] == 'electrons per Bohr^3'
assert hashlib.sha256((p / 'density.bin').read_bytes()).hexdigest() == m['sha256']
assert r['convergence']['converged']
print(r['energies']['xc']['ha'])
"#,
            ])
            .arg(work.join(name))
            .arg(dims[0].to_string())
            .arg(l.to_string())
            .output()
            .unwrap();
        assert!(
            metadata.status.success(),
            "{}",
            String::from_utf8_lossy(&metadata.stderr)
        );
        let stored_xc: f64 = String::from_utf8(metadata.stdout)
            .unwrap()
            .trim()
            .parse()
            .unwrap();
        let bytes = std::fs::read(work.join(name).join("density.bin")).unwrap();
        assert_eq!(bytes.len(), dims.iter().product::<usize>() * 8);
        let valence = Array3::from_shape_vec(
            dims,
            bytes
                .as_chunks::<8>()
                .0
                .iter()
                .map(|b| f64::from_le_bytes(*b))
                .collect(),
        )
        .unwrap();
        let raw = &valence + &core;
        assert!(raw.iter().all(|n| n.is_finite()));
        // Explicit density_floor in material TOMLs acts before the FFT derivative.
        let n = raw.mapv(|n| n.max(1e-10));
        let s = sigma(&n, &cell, &mut fft);
        let engine = gga::evaluate(&n, &cell, &mut fft).unwrap().energy.0;
        println!(
            "{name} reconstructed_minus_stored_xc={:.16e} Ha",
            engine - stored_xc
        );
        assert!(
            (engine - stored_xc).abs() <= CAP_HA,
            "STOP: density/core reconstruction mismatch"
        );
        let modified = energy(&n, &s, &cell, MODIFIED);
        let error = modified - engine;
        println!(
            "{name} engine_xc={engine:.16e} evaluator_xc={modified:.16e} difference={error:.16e} Ha"
        );
        assert!(error.abs() <= CAP_HA, "STOP: invalid probe: {error:e} Ha");
        let shift = energy(&n, &s, &cell, ORIGINAL) - modified;
        let mut sorted: Vec<_> = s.iter().copied().collect();
        sorted.sort_by(f64::total_cmp);
        println!("{name} shift={shift:.16e} Ha = {:.16e} Ry", 2.0 * shift);
        println!(
            "{name} points={} min_density={:.16e} e/Bohr^3 below_1e-6={} below_1e-10={} negative={} sigma_below_1e-10={} smallest_sigma={:?} Bohr^-8",
            raw.len(),
            raw.iter().copied().fold(f64::INFINITY, f64::min),
            raw.iter().filter(|&&n| n < 1e-6).count(),
            raw.iter().filter(|&&n| n < 1e-10).count(),
            raw.iter().filter(|&&n| n < 0.).count(),
            s.iter().filter(|&&s| s < 1e-10).count(),
            &sorted[..8]
        );
        let threshold_shift = cell.volume_bohr3() / n.len() as f64
            * Array1::from_iter(n.iter().zip(&s).map(|(&n, &s)| {
                if n <= 1e-6 || s <= 1e-10 {
                    n * (pbe(n, 0., MODIFIED) - pbe(n, s, MODIFIED))
                } else {
                    0.
                }
            }))
            .sum();
        println!("{name} suppress_gradient_correction_shift={threshold_shift:.16e} Ha");
    }
}
