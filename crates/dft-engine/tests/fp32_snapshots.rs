//! P47B measurement only. P14A_WORK_ACCOUNTING=1 cargo test --release -p dft-engine
//! --test fp32_snapshots -- --ignored --nocapture --test-threads=1
//! External P47B_UPF; output P47B-evidence. No engine precision option.
use dft_basis::{FftGridSpec, HalfSphere, PlaneWaveBasis};
use dft_engine::{
    CompletionMode, ExecutionConfig, GammaReal, ScfOptions, ScfProblem, ScfSession,
    SpeciesPotential, StandardProfile, StepStatus,
};
use dft_geometry::{Atom, Cell};
use dft_numerics::real::RealMatrix;
use dft_physics::{
    data::*,
    hamiltonian::{self, Hamiltonian, NonlocalOperator, Species},
    hartree,
    real_gamma::RealHamiltonian,
    xc,
};
use std::{io::Write, path::PathBuf, time::Instant};

fn numbers(text: &str, key: &str) -> Vec<f64> {
    text.lines()
        .find_map(|l| {
            l.split_once('=')
                .filter(|(k, _)| k.trim() == key)
                .map(|(_, v)| v)
        })
        .unwrap()
        .split(['[', ']', ',', ' ', '\t'])
        .filter(|x| !x.is_empty())
        .map(|x| x.parse().unwrap())
        .collect()
}
fn save(
    name: &str,
    cell: &Cell,
    atoms: &[Atom],
    pp: &dft_upf::Pseudopotential,
    density: &Array3<f64>,
    vectors: &ComplexMatrix,
) {
    let basis = PlaneWaveBasis::new(cell, Hartree(15.)).unwrap();
    let db = basis.density_basis(cell).unwrap();
    let dims = std::array::from_fn(|i| density.shape()[i]);
    let dg = FftGridSpec::new(&db, &db, Some(dims)).unwrap();
    let pos: Vec<_> = atoms
        .iter()
        .map(|a| cell.fractional_to_cartesian_bohr(a.fractional_position))
        .collect();
    let species = [Species {
        potential: pp,
        positions_bohr: &pos,
    }];
    let mut fft = FftPlan3d::new(dims).unwrap();
    let ig = dg
        .scatter(&hamiltonian::ionic_local(&db, cell.volume_bohr3(), &species).unwrap())
        .unwrap();
    let mut ir = Array3::zeros(dims);
    fft.to_real(&ig, &mut ir).unwrap();
    let hp = hartree::evaluate(density, cell, &mut fft).unwrap();
    let xp = xc::evaluate(density, cell.volume_bohr3()).unwrap();
    let mut potential = ir.mapv(|z| z.re);
    for ((v, h), x) in potential.iter_mut().zip(&hp.potential).zip(&xp.potential) {
        *v += h.re + x;
    }
    let h = Hamiltonian::new(
        cell,
        &basis,
        potential,
        NonlocalOperator::new(&basis, cell.volume_bohr3(), &species).unwrap(),
    )
    .unwrap();
    let half = HalfSphere::new(&basis, h.fft_grid()).unwrap();
    let n = basis.vectors().len();
    let b = vectors.cols();
    let mut x = RealMatrix::zeros(n, b).unwrap();
    for j in 0..b {
        half.pack(
            &vectors.as_slice()[j * n..(j + 1) * n],
            &mut x.as_mut_slice()[j * n..(j + 1) * n],
        )
        .unwrap();
    }
    let p: usize = h.nonlocal().atoms().iter().map(|a| a.columns().len()).sum();
    let mut proj = vec![0.; n * p];
    let mut coupling = vec![0.; p * p];
    let mut offset = 0;
    for atom in h.nonlocal().atoms() {
        let c = atom.columns().len();
        for (j, beta) in atom.columns().iter().enumerate() {
            let out = &mut proj[(offset + j) * n..(offset + j + 1) * n];
            out[0] = beta[0].re;
            for (i, &g) in half.plus().iter().enumerate().skip(1) {
                out[2 * i - 1] = std::f64::consts::SQRT_2 * beta[g as usize].re;
                out[2 * i] = std::f64::consts::SQRT_2 * beta[g as usize].im;
            }
        }
        for i in 0..c {
            for j in 0..c {
                coupling[(offset + j) * p + offset + i] = atom.coupling_hartree()[i * c + j];
            }
        }
        offset += c;
    }
    let mut kin = vec![0.; n];
    for (i, &g) in half.plus().iter().enumerate().skip(1) {
        kin[2 * i - 1] = h.kinetic_diagonal()[g as usize];
        kin[2 * i] = kin[2 * i - 1];
    }
    let mut real =
        RealHamiltonian::from_complex(&h, &basis, b, WorkerPool::new(1).unwrap()).unwrap();
    let mut y = RealMatrix::zeros(n, b).unwrap();
    real.apply_block(&x, &mut y).unwrap();
    let mut samples = Vec::new();
    for _ in 0..9 {
        let t = Instant::now();
        real.apply_block(&x, &mut y).unwrap();
        samples.push(t.elapsed().as_secs_f64());
    }
    println!("PRODUCTION_H,{name},{n},{b},{p},{dims:?},{samples:?}");
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../P47B-evidence");
    let mut f =
        std::io::BufWriter::new(std::fs::File::create(root.join(format!("{name}.bin"))).unwrap());
    for k in dims.into_iter().chain([n, b, p]) {
        f.write_all(&(k as u64).to_le_bytes()).unwrap();
    }
    for list in [half.plus(), half.minus()] {
        for &g in list {
            let ix = h.fft_grid().g_to_fft()[g as usize];
            let k = (ix[0] * dims[1] + ix[1]) * dims[2] + ix[2];
            f.write_all(&(k as u64).to_le_bytes()).unwrap();
        }
    }
    for slice in [
        h.local_potential_values(),
        kin.as_slice(),
        x.as_slice(),
        proj.as_slice(),
        coupling.as_slice(),
        y.as_slice(),
    ] {
        for v in slice {
            f.write_all(&v.to_le_bytes()).unwrap();
        }
    }
}
#[test]
#[ignore = "external Si UPF; SCFs and snapshots, release measurement only"]
fn histories_and_snapshots() {
    if cfg!(debug_assertions) {
        panic!("run SCF measurements in release mode");
    }
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    std::fs::create_dir_all(root.join("P47B-evidence")).unwrap();
    let pp = dft_upf::parse_upf(
        &std::fs::read_to_string(std::env::var("P47B_UPF").unwrap_or_else(|_| {
            root.join("../pseudo/Si.pz-vbc.UPF")
                .to_string_lossy()
                .into()
        }))
        .unwrap(),
    )
    .unwrap();
    for size in [8, 16, 32] {
        let text =
            std::fs::read_to_string(root.join(format!("benchmarks/cases/si{size}-gamma.toml")))
                .unwrap();
        let a = numbers(&text, "lattice");
        let cell = Cell::from_user_rows_bohr(std::array::from_fn(|i| {
            std::array::from_fn(|j| a[3 * i + j])
        }))
        .unwrap();
        let atoms: Vec<_> = numbers(&text, "atoms")
            .as_chunks::<3>()
            .0
            .iter()
            .map(|p| Atom::new("Si", [p[0], p[1], p[2]]))
            .collect();
        println!(
            "BASIS,si{size},{}",
            PlaneWaveBasis::new(&cell, Hartree(15.))
                .unwrap()
                .vectors()
                .len()
        );
        let species = [SpeciesPotential {
            label: "Si",
            potential: &pp,
        }];
        let options = ScfOptions {
            gamma_real: GammaReal::On,
            standard: Some(StandardProfile::default()),
            threads: 1,
            cutoff_hartree: 15.,
            ..Default::default()
        };
        let make = || ScfProblem::new(&cell, &atoms, &species, options.clone()).unwrap();
        for seed in 0..3 {
            let source = make();
            let mut density = source.uniform_density();
            if seed > 0 {
                let mut s = seed as u64;
                for x in &mut density {
                    s = s.wrapping_mul(6364136223846793005).wrapping_add(1);
                    *x *= 1. + 1e-10 * ((s >> 11) as f64 / ((1_u64 << 53) as f64) - 0.5);
                }
            }
            // Exclude the extra source model used only to seed density transfer.
            let start = Instant::now();
            let mut session = ScfSession::new(
                make(),
                CompletionMode::RetainLast,
                ExecutionConfig::default(),
            )
            .unwrap();
            session
                .initialize_density_transfer(&source, &density)
                .unwrap();
            dft_engine::work_accounting::take();
            let last_input = loop {
                let input = session.density().unwrap().clone();
                let r = session.step().unwrap();
                if r.status != StepStatus::Ready {
                    assert_eq!(r.status, StepStatus::Converged);
                    break input;
                }
            };
            let wall = start.elapsed().as_secs_f64();
            for r in session.history() {
                println!(
                    "HISTORY,si{size},{seed},{},{:.17e},{},{},{:.9},{},{}",
                    r.iteration,
                    r.requested_eigen_tolerance,
                    r.davidson_iterations,
                    r.operator_applications,
                    r.timings.diagonalization.as_secs_f64(),
                    r.standard_redo,
                    r.standard_error_ha.unwrap()
                );
            }
            for w in dft_engine::work_accounting::take() {
                println!(
                    "COUNT,si{size},{seed},{:.17e},{},{},{},{}",
                    w.tolerance, w.iterations, w.vector_applications, w.block_calls, w.returned
                );
            }
            println!(
                "SCF,si{size},{seed},{wall:.9},{:.17e},{}",
                session.result().unwrap().result.energy.total_hartree(),
                session.history().len()
            );
            if seed == 0 && std::env::var("P47B_SKIP_SNAPSHOTS").is_err() {
                save(
                    &format!("si{size}-converged"),
                    &cell,
                    &atoms,
                    &pp,
                    &last_input,
                    &session.result().unwrap().result.k_points[0]
                        .solution
                        .eigenvectors,
                );
                let mut problem = make();
                let initial = problem.uniform_density();
                let early = problem.step(&initial, None, 5e-3).unwrap();
                save(
                    &format!("si{size}-early"),
                    &cell,
                    &atoms,
                    &pp,
                    &initial,
                    &early.solution.eigenvectors,
                );
            }
        }
    }
}
