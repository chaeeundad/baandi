//! Real-grid equivalence and six changed-convention SCFs; defaults unchanged.
//! Sources: PW92 PRB 45, 13244, Eq. (10), Table I; PBE PRL 77, 3865,
//! Eqs. (7), (8), (10), (14). Full-precision beta/mu: ADR 0018.
//! Real runs: build release dft-cli, then set DFT_NC_DIR and run this test
//! in release with --ignored --nocapture. See xc-prefactor-diagnostic.md.
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_numerics::FftPlan3d;
use dft_physics::pbe_pw92::Pw92Prefactor as A;
use dft_physics::{density::core_density, gga, hamiltonian::Species};
use dft_units::Hartree;
use ndarray::Array3;
use std::{path::PathBuf, process::Command};
const CAP_HA: f64 = 1e-12;

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
    let work = std::env::var_os("PW92_WORK")
        .map(PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join(format!("p27b-xc-{}", std::process::id())));
    std::fs::create_dir_all(&work).unwrap();
    let cases = std::env::var("PW92_CASES")
        .unwrap_or_else(|_| "si-20,si-30,si-40,mgo-40,mgo-50,mgo-60".into());
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
        let old = gga::evaluate(&n, &cell, &mut fft).unwrap();
        let new = gga::evaluate_pw92(A::Modified, &n, &cell, &mut fft).unwrap();
        let verr = old
            .potential
            .iter()
            .zip(&new.potential)
            .map(|(a, b)| (a - b).abs())
            .fold(0.0_f64, f64::max);
        assert!((old.energy.0 - stored_xc).abs() <= CAP_HA);
        assert!((old.energy.0 - new.energy.0).abs() <= CAP_HA);
        assert!(verr <= CAP_HA);
        println!(
            "{name} real-grid dE={} Ha max_dV={verr:e} Ha",
            new.energy.0 - old.energy.0
        );
        let script = r#"
import hashlib, json, os, pathlib, subprocess, sys
root, work, name = pathlib.Path(sys.argv[1]), pathlib.Path(sys.argv[2]), sys.argv[3]
text = (work / (name+'.toml')).read_text().replace('[xc]', '[xc]\npw92_prefactor = "original"')
inp=work/(name+'-original.toml'); inp.write_text(text)
out=work/(name+'-original')
r=subprocess.run([str(root/'target/release/dft-cli'),'run',str(inp),'--output-dir',str(out)],capture_output=True,text=True)
(work/(name+'-original.log')).write_text(r.stdout+r.stderr)
assert r.returncode==0, r.stderr
original=json.loads((out/'results.json').read_text())
modified=json.loads((work/name/'results.json').read_text())
assert original['convergence']['converged'] and modified['convergence']['converged']
assert original['pw92_prefactor']=='original'
# Compare default PBE against the pre-change executable when supplied.
if os.environ.get('PW92_BASELINE_CLI'):
    baseline_out=work/(name+'-baseline')
    proc=subprocess.run([os.environ['PW92_BASELINE_CLI'],'run',str(work/(name+'.toml')),'--output-dir',str(baseline_out)],capture_output=True,text=True)
    assert proc.returncode==0, proc.stderr
    baseline=json.loads((baseline_out/'results.json').read_text())
    assert baseline['total_energy']==modified['total_energy']
    assert (baseline_out/'density.bin').read_bytes()==(work/name/'density.bin').read_bytes()
    print(name, 'default PBE baseline energy/density bit-identical')
# Explicit modified is equivalent to omission, including every density bit.
mod_inp=work/(name+'-modified.toml'); mod_inp.write_text(text.replace('"original"','"modified"'))
mod_out=work/(name+'-modified')
proc=subprocess.run([str(root/'target/release/dft-cli'),'run',str(mod_inp),'--output-dir',str(mod_out)],capture_output=True,text=True)
assert proc.returncode==0, proc.stderr
explicit=json.loads((mod_out/'results.json').read_text())
assert explicit['total_energy']==modified['total_energy']
assert (mod_out/'density.bin').read_bytes()==(work/name/'density.bin').read_bytes()
e0=modified['total_energy']['ha']; e1=original['total_energy']['ha']
prediction=-2.1002e-6 if name.startswith('si') else -6.8706e-6
row=dict(case=name,modified_ha=e0,original_ha=e1,delta_ry=2*(e1-e0),prediction_ry=prediction,iterations=original['convergence']['iterations'])
(work/(name+'-comparison.json')).write_text(json.dumps(row,indent=2)+'\n')
print(json.dumps(row))
"#;
        assert!(
            Command::new("python3")
                .args(["-I", "-c", script])
                .arg(&root)
                .arg(&work)
                .arg(name)
                .status()
                .unwrap()
                .success()
        );
    }
}
