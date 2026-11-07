//! Test-only radial convention sensitivity. See docs/validation/mgo-residual-diagnostic.md.
use dft_upf::{Pseudopotential, radial::spherical_bessel};
use std::{f64::consts::PI, path::PathBuf};

#[derive(Clone, Copy)]
struct Convention {
    limit: Option<f64>,
    odd: bool,
    trap: bool,
}
const BASE: Convention = Convention {
    limit: None,
    odd: false,
    trap: false,
};
const CASES: [(&str, Convention); 5] = [
    ("odd", Convention { odd: true, ..BASE }),
    (
        "r9",
        Convention {
            limit: Some(9.0),
            ..BASE
        },
    ),
    (
        "r9.5",
        Convention {
            limit: Some(9.5),
            ..BASE
        },
    ),
    (
        "r10.5",
        Convention {
            limit: Some(10.5),
            ..BASE
        },
    ),
    ("trap", Convention { trap: true, ..BASE }),
];

// Index-coordinate quadrature of samples already multiplied by PP_RAB.
// Independent weight assembly, including overlapping Simpson/3/8 endpoint.
fn integral(f: &[f64], trap: bool) -> f64 {
    let n = f.len();
    if n < 2 {
        return 0.0;
    }
    let mut w = vec![0.0; n];
    if trap || n == 2 {
        w.fill(1.0);
        w[0] = 0.5;
        w[n - 1] = 0.5;
    } else {
        let end = if n.is_multiple_of(2) { n - 4 } else { n - 1 };
        for i in (0..end).step_by(2) {
            for (j, v) in [1.0 / 3.0, 4.0 / 3.0, 1.0 / 3.0].iter().enumerate() {
                w[i + j] += v;
            }
        }
        if n.is_multiple_of(2) {
            for (j, v) in [3.0 / 8.0, 9.0 / 8.0, 9.0 / 8.0, 3.0 / 8.0]
                .iter()
                .enumerate()
            {
                w[end + j] += v;
            }
        }
    }
    // Compensated sum keeps tiny differences distinguishable from summation noise.
    let (mut sum, mut correction) = (0.0, 0.0);
    for (x, w) in f.iter().zip(w) {
        let y = x * w - correction;
        let t = sum + y;
        correction = (t - sum) - y;
        sum = t;
    }
    sum
}

// kind 0: volume-free local (Ha Bohr^3), kind 1: core (electrons),
// kind >=2: projector radial factor (Bohr^3/2).
fn transform(p: &Pseudopotential, kind: usize, q: f64, c: Convention) -> f64 {
    let r = p.mesh().r();
    let available = if kind >= 2 {
        p.projectors()[kind - 2].radial_times_r().len()
    } else {
        r.len()
    };
    let limit = c.limit.or(if kind == 0 { Some(10.0) } else { None });
    let mut n = limit.map_or(available, |x| r[..available].partition_point(|r| *r <= x));
    if c.odd {
        // Truncate the mesh endpoint, not the separately stored projector prefix.
        if kind >= 2 {
            let mesh_n = r.len() - usize::from(r.len().is_multiple_of(2));
            n = n.min(mesh_n);
        } else if n.is_multiple_of(2) {
            n -= 1;
        }
    }
    let samples: Vec<_> = (0..n)
        .map(|i| {
            let (f, l) = match kind {
                0 => (
                    r[i] * r[i] * p.local_hartree()[i] + p.z_valence() * r[i] * libm::erf(r[i]),
                    0,
                ),
                1 => (r[i] * r[i] * p.core_density().map_or(0.0, |rho| rho[i]), 0),
                _ => (
                    r[i] * p.projectors()[kind - 2].radial_times_r()[i],
                    p.projectors()[kind - 2].l(),
                ),
            };
            f * spherical_bessel(l, q * r[i]).unwrap() * p.mesh().rab()[i]
        })
        .collect();
    let value = 4.0 * PI * integral(&samples, c.trap);
    if kind != 0 {
        value
    } else if q == 0.0 {
        value + PI * p.z_valence()
    } else {
        value - 4.0 * PI * p.z_valence() * (-q * q / 4.0).exp() / (q * q)
    }
}

// Interpolate only the smooth short-range local part; Coulomb pole stays analytic.
fn table(p: &Pseudopotential, kind: usize, q: f64) -> f64 {
    let start = ((q / 0.01).floor() as usize).saturating_sub(1);
    let pole = |x: f64| {
        if kind == 0 && x > 0.0 {
            -4.0 * PI * p.z_valence() * (-x * x / 4.0).exp() / (x * x)
        } else {
            0.0
        }
    };
    let mut value = 0.0;
    for i in 0..4 {
        let x = (start + i) as f64 * 0.01;
        let mut w = 1.0;
        for j in 0..4 {
            if i != j {
                let y = (start + j) as f64 * 0.01;
                w *= (q - y) / (x - y);
            }
        }
        let f = transform(p, kind, x, BASE)
            - pole(x)
            - if kind == 0 && x == 0.0 {
                PI * p.z_valence()
            } else {
                0.0
            };
        value += w * f;
    }
    value + pole(q)
}

#[test]
fn quadrature_polynomials() {
    for n in [5, 6, 11, 12] {
        let end = (n - 1) as f64;
        for power in 0..=3 {
            let f: Vec<_> = (0..n).map(|i| (i as f64).powi(power)).collect();
            assert!(
                (integral(&f, false) - end.powi(power + 1) / f64::from(power + 1)).abs() < 1e-10
            );
        }
        assert!(
            (integral(&(0..n).map(|i| i as f64).collect::<Vec<_>>(), true) - end * end / 2.0).abs()
                < 1e-12
        );
    }
}

#[test]
#[ignore = "external PseudoDojo files; set DFT_NC_DIR; release command in diagnostic doc"]
fn radial_conventions() {
    let Some(dir) = std::env::var_os("DFT_NC_DIR").map(PathBuf::from) else {
        eprintln!("SKIP: set DFT_NC_DIR");
        return;
    };
    if ["Si", "Mg", "O"]
        .iter()
        .any(|s| !dir.join(format!("{s}.upf")).is_file())
    {
        eprintln!("SKIP: Si/Mg/O UPFs absent");
        return;
    }
    let mut alpha = [[0.0; 5]; 3];
    for (is, s) in ["Si", "Mg", "O"].iter().enumerate() {
        let p = dft_upf::from_file(dir.join(format!("{s}.upf"))).unwrap();
        println!(
            "SPECIES {s} mesh={} z={} alpha={:.16e} Ha Bohr^3",
            p.mesh().r().len(),
            p.z_valence(),
            transform(&p, 0, 0.0, BASE)
        );
        for kind in 0..2 + p.projectors().len() {
            for q in [0.0, 0.137, 0.513, 1.017] {
                let base = transform(&p, kind, q, BASE);
                let engine = match kind {
                    0 => {
                        let t = p.local_transform(1.0, None).unwrap();
                        if q == 0.0 {
                            t.alpha_z().unwrap()
                        } else {
                            t.at_g(q).unwrap()
                        }
                    }
                    1 => p.core_density().map_or(0.0, |rho| {
                        dft_upf::radial::spherical_density_transform(p.mesh(), rho, q).unwrap()
                    }),
                    _ => p.projector_transform(kind - 2, q).unwrap(),
                };
                assert!(
                    (base - engine).abs() < 1e-11,
                    "scalar baseline differs: {s} {kind} {q}"
                );
                for (ic, (label, c)) in CASES.iter().enumerate() {
                    let delta = transform(&p, kind, q, *c) - base;
                    assert!(delta.is_finite());
                    println!("DELTA {s} kind={kind} q={q} {label} {delta:.16e}");
                    if kind == 0 && q == 0.0 {
                        alpha[is][ic] = delta;
                    }
                }
                if q > 0.0 {
                    println!(
                        "TABLE {s} kind={kind} q={q} delta={:.16e}",
                        table(&p, kind, q) - base
                    );
                }
            }
        }
    }
    for (ic, (label, _)) in CASES.iter().enumerate() {
        let si = 2.0 * 8.0 * 2.0 * alpha[0][ic] / (2.0 * 5.13_f64.powi(3));
        let mgo = 2.0 * 16.0 * (alpha[1][ic] + alpha[2][ic]) / (2.0 * 3.975_f64.powi(3));
        println!(
            "ENERGY {label} alt-base Si={si:.16e} Ry MgO={mgo:.16e} Ry; base-alt/residual Si={:.8} MgO={:.8}",
            -si / 7.97e-8,
            -mgo / 1.279e-6
        );
    }
}

#[test]
#[ignore = "fixed converged densities: DFT_NC_DIR and RADIAL_DIAG_DENSITIES; see doc"]
fn fixed_density_local() {
    use dft_basis::{FftGridSpec, PlaneWaveBasis};
    use dft_geometry::Cell;
    use dft_physics::data::*;
    use std::{collections::BTreeMap, process::Command};
    let (Some(nc), Some(work)) = (
        std::env::var_os("DFT_NC_DIR"),
        std::env::var_os("RADIAL_DIAG_DENSITIES"),
    ) else {
        eprintln!("SKIP: set DFT_NC_DIR and RADIAL_DIAG_DENSITIES");
        return;
    };
    let (nc, work) = (PathBuf::from(nc), PathBuf::from(work));
    if ["Si", "Mg", "O"]
        .iter()
        .any(|s| !nc.join(format!("{s}.upf")).is_file())
    {
        eprintln!("SKIP: Si/Mg/O UPFs absent");
        return;
    }
    for (name, l, cutoff, d, els, second, nel, residual) in [
        ("si-20", 5.13, 20.0, 30, ["Si", "Si"], 0.25, 8.0, 7.97e-8),
        ("mgo-40", 3.975, 40.0, 36, ["Mg", "O"], 0.5, 16.0, 1.279e-6),
    ] {
        let path = work.join(name);
        // Verify provenance and the density contract; compare reconstructed local energy too.
        let metadata = Command::new("python3")
            .args([
                "-I",
                "-c",
                r#"
import hashlib,json,pathlib,sys
p=pathlib.Path(sys.argv[1]); d=int(sys.argv[2]); l=float(sys.argv[3])
m=json.loads((p/'density.json').read_text()); r=json.loads((p/'results.json').read_text())
assert m['fft_grid']==[d,d,d] and m['lattice_columns_bohr']==[[0,l,l],[l,0,l],[l,l,0]]
assert m['dtype']=='f64' and m['endianness']=='little'
assert m['ordering']=='C: (i*n2+j)*n3+k; k fastest' and m['units']=='electrons per Bohr^3'
assert hashlib.sha256((p/'density.bin').read_bytes()).hexdigest()==m['sha256']
assert r['convergence']['converged']
for pp in r['species']:
 assert hashlib.sha256((pathlib.Path(sys.argv[4])/pathlib.Path(pp['path']).name).read_bytes()).hexdigest()==pp['sha256']
print(r['energies']['local']['ha'], r['energies']['xc']['ha'])
"#,
            ])
            .arg(&path)
            .arg(d.to_string())
            .arg(l.to_string())
            .arg(&nc)
            .output()
            .unwrap();
        assert!(
            metadata.status.success(),
            "{}",
            String::from_utf8_lossy(&metadata.stderr)
        );
        let stored_values: Vec<f64> = String::from_utf8(metadata.stdout)
            .unwrap()
            .split_whitespace()
            .map(|s| s.parse().unwrap())
            .collect();
        let stored = stored_values[0];
        let bytes = std::fs::read(path.join("density.bin")).unwrap();
        assert_eq!(bytes.len(), d * d * d * 8);
        let density = Array3::from_shape_vec(
            [d; 3],
            bytes
                .as_chunks::<8>()
                .0
                .iter()
                .map(|b| Complex64::new(f64::from_le_bytes(*b), 0.0))
                .collect(),
        )
        .unwrap();
        let mut fft = FftPlan3d::new([d; 3]).unwrap();
        let mut ng = Array3::zeros([d; 3]);
        fft.to_reciprocal(&density, &mut ng).unwrap();
        let cell = Cell::from_user_rows_bohr([[0., l, l], [l, 0., l], [l, l, 0.]]).unwrap();
        let volume = cell.volume_bohr3();
        assert!((ng[[0, 0, 0]].re * volume - nel).abs() < 1e-8);
        let basis = PlaneWaveBasis::new(&cell, Hartree(cutoff)).unwrap();
        let db = basis.density_basis(&cell).unwrap();
        let grid = FftGridSpec::new(&db, &db, Some([d; 3])).unwrap();
        let ps: Vec<_> = els
            .iter()
            .map(|e| dft_upf::from_file(nc.join(format!("{e}.upf"))).unwrap())
            .collect();
        let positions = [[0.; 3], cell.fractional_to_cartesian_bohr([second; 3])];
        let mut energy = 0.0;
        let mut delta = [0.0; 6];
        for (p, pos) in ps.iter().zip(positions) {
            let mut cache = BTreeMap::new();
            for g in db.vectors() {
                let q = g.norm_squared.sqrt();
                let (base, changes) = cache.entry(q.to_bits()).or_insert_with(|| {
                    let base = transform(p, 0, q, BASE);
                    let mut changes = [0.0; 6];
                    for (i, (_, c)) in CASES.iter().enumerate() {
                        changes[i] = transform(p, 0, q, *c) - base;
                    }
                    // Alpha is independently evaluated, not interpolated across the Coulomb pole.
                    changes[5] = if q == 0.0 { 0.0 } else { table(p, 0, q) - base };
                    (base, changes)
                });
                let ix = FftGridSpec::map_miller(g.miller, grid.dims()).unwrap();
                let phase = Complex64::from_polar(
                    1.0,
                    -g.cartesian.iter().zip(pos).map(|(g, r)| g * r).sum::<f64>(),
                );
                let weight = (ng[ix].conj() * phase).re;
                energy += weight * *base;
                for (out, change) in delta.iter_mut().zip(changes.iter()) {
                    *out += weight * change;
                }
            }
        }
        assert!(
            (energy - stored).abs() < 1e-10,
            "local reconstruction: {energy} vs {stored}"
        );
        println!(
            "FIXED {name} baseline={energy:.16e} Ha stored_difference={:.16e} Ha",
            energy - stored
        );
        let mut core_g = vec![Array3::<Complex64>::zeros([d; 3]); 7];
        for (p, pos) in ps.iter().zip(positions) {
            if p.core_density().is_none() {
                continue;
            }
            let mut cache = BTreeMap::new();
            for g in db.vectors() {
                let q = g.norm_squared.sqrt();
                let values = cache.entry(q.to_bits()).or_insert_with(|| {
                    let mut values = [transform(p, 1, q, BASE); 7];
                    for (i, (_, c)) in CASES.iter().enumerate() {
                        values[i + 1] = transform(p, 1, q, *c);
                    }
                    values[6] = table(p, 1, q);
                    values
                });
                let ix = FftGridSpec::map_miller(g.miller, grid.dims()).unwrap();
                let phase = Complex64::from_polar(
                    1.0,
                    -g.cartesian.iter().zip(pos).map(|(g, r)| g * r).sum::<f64>(),
                );
                for (out, value) in core_g.iter_mut().zip(values.iter()) {
                    out[ix] += phase * *value / volume;
                }
            }
        }
        let mut xc = [0.0; 7];
        for (out, core) in xc.iter_mut().zip(core_g) {
            let mut real = Array3::zeros([d; 3]);
            fft.to_real(&core, &mut real).unwrap();
            let total = (&density + &real).mapv(|v| v.re.max(1e-10));
            *out = dft_physics::gga::evaluate(&total, &cell, &mut fft)
                .unwrap()
                .energy
                .0;
        }
        assert!(
            (xc[0] - stored_values[1]).abs() < 1e-12,
            "core/XC reconstruction"
        );
        println!("CORE {name} baseline_xc={:.16e} Ha", xc[0]);
        for (i, label) in ["odd", "r9", "r9.5", "r10.5", "trap", "table"]
            .iter()
            .enumerate()
        {
            println!(
                "CORE {name} {label} alt-base={:.16e} Ry",
                2.0 * (xc[i + 1] - xc[0])
            );
            println!(
                "FIXED {name} {label} alt-base={:.16e} Ry base-alt/residual={:.10}",
                2.0 * delta[i],
                -2.0 * delta[i] / residual
            );
        }
    }
}
