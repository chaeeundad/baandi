//! MPI entry-cost profile (work order beat-QE E2, docs/performance/mpi-entry-cost-analysis.md):
//! the same Si SCF on the threads-only path (one rank, T=1) and on the M2
//! stick-row path with simulated ranks, printing the wall, the distributed
//! Davidson stages, the M2 operator stages and the collective counts. Simulated
//! ranks are threads of one process, so every process-wide timer is the sum
//! over ranks; the report divides by the rank count. Mac observations only
//! (pure build, no MKL): this isolates the structural cost of the distributed
//! path at equal kernels, it is not dev08 evidence.
//!
//! DFT_SI_UPF=/path/Si.pz-vbc.UPF BAANDI_DAVIDSON_STAGES=on \
//!   cargo test --release -p dft-engine --test m2_entry_cost -- --ignored --nocapture
//! Optional: M2_PROFILE_CASE=si8, si16 or si32 (default si8), M2_PROFILE_RANKS=2 (comma list).
#[allow(dead_code)]
mod common;
use dft_engine::{GammaReal, ScfOptions, ScfProblem, ScfResult, SpeciesPotential, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::sync::Arc;
use std::time::Instant;

const A: f64 = 10.261_212_856_717_933;
const BASE: [[f64; 3]; 8] = [
    [0., 0., 0.],
    [0., 0.5, 0.5],
    [0.5, 0., 0.5],
    [0.5, 0.5, 0.],
    [0.25, 0.25, 0.25],
    [0.25, 0.75, 0.75],
    [0.75, 0.25, 0.75],
    [0.75, 0.75, 0.25],
];

fn structure(case: &str) -> (Cell, Vec<Atom>) {
    let (rx, ry) = match case {
        "si16" => (2, 1),
        "si32" => (2, 2),
        _ => (1, 1),
    };
    let cell = Cell::from_user_rows_bohr([
        [A * rx as f64, 0., 0.],
        [0., A * ry as f64, 0.],
        [0., 0., A],
    ])
    .unwrap();
    // Fixed fractional displacements of amplitude 1e-3 (deterministic LCG, seed 601).
    let mut s: u64 = 601;
    let mut next = || {
        s = s
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((s >> 11) as f64 / (1u64 << 53) as f64) * 2. - 1.
    };
    let mut atoms = Vec::new();
    for i in 0..rx {
        for j in 0..ry {
            for p in BASE {
                let q = [
                    (p[0] + i as f64) / rx as f64 + 1e-3 * next(),
                    (p[1] + j as f64) / ry as f64 + 1e-3 * next(),
                    p[2] + 1e-3 * next(),
                ];
                atoms.push(Atom::new("Si", q));
            }
        }
    }
    (cell, atoms)
}

fn run(case: &str, comm: Option<Arc<dyn Collective + Send + Sync>>) -> ScfResult {
    let path = std::env::var("DFT_SI_UPF").expect("set DFT_SI_UPF to the external Si UPF");
    let pp = dft_upf::parse_upf(&std::fs::read_to_string(path).unwrap()).unwrap();
    let (cell, atoms) = structure(case);
    let m2 = comm.is_some();
    if m2 {
        parallel::set_m2_override(Some(true));
    }
    parallel::set_communicator(comm);
    let result = ScfProblem::new(
        &cell,
        &atoms,
        &[SpeciesPotential {
            label: "Si",
            potential: &pp,
        }],
        ScfOptions {
            gamma_real: GammaReal::Auto,
            cutoff_hartree: 15.,
            threads: 1,
            ..Default::default()
        },
    )
    .unwrap()
    .run()
    .unwrap();
    parallel::set_communicator(None);
    parallel::set_m2_override(None);
    result
}

fn work(r: &ScfResult) -> (usize, usize, usize) {
    (
        r.history.len(),
        r.history.iter().map(|h| h.davidson_iterations).sum(),
        r.history.iter().map(|h| h.operator_applications).sum(),
    )
}

type Snap = (
    Vec<(&'static str, f64)>,
    Vec<(String, f64)>,
    Vec<(&'static str, usize, f64)>,
);
fn snapshot() -> Snap {
    (
        parallel::stage_seconds(),
        parallel::m2_stats().1,
        parallel::collective_stats(),
    )
}

#[test]
#[ignore = "requires DFT_SI_UPF; release only; Mac observation of the M2 entry cost"]
fn m2_entry_cost_profile() {
    let case = std::env::var("M2_PROFILE_CASE").unwrap_or_else(|_| "si8".into());
    let ranks: Vec<usize> = std::env::var("M2_PROFILE_RANKS")
        .unwrap_or_else(|_| "2".into())
        .split(',')
        .map(|s| s.trim().parse().unwrap())
        .collect();
    let t = Instant::now();
    let serial = run(&case, None);
    let serial_wall = t.elapsed().as_secs_f64();
    let st = dft_solvers::davidson_stages::stage_seconds();
    let wt = dft_numerics::work_timers::snapshot();
    let wt = wt[0].plus(&wt[1]).plus(&wt[2]);
    println!(
        "{case} serial T=1: wall {serial_wall:.3} s, E {:.12} Ha, (SCF, Davidson, H) {:?}",
        serial.total_hartree(),
        work(&serial)
    );
    println!(
        "  serial stages (BAANDI_DAVIDSON_STAGES=on): reduced {:.3} eigen {:.3} rotate {:.3} residual {:.3} orthonormalize {:.3} h_apply {:.3} precondition {:.3}",
        st.reduced,
        st.eigen,
        st.rotate,
        st.residual,
        st.orthonormalize,
        st.h_apply,
        st.precondition
    );
    println!(
        "  serial h_apply parts: apply {:.3} local {:.3} nonlocal_projection {:.3} nonlocal_back {:.3} columns {}",
        wt.apply, wt.local, wt.nonlocal_projection, wt.nonlocal_back, wt.apply_columns
    );
    for &p in &ranks {
        let before = snapshot();
        let t = Instant::now();
        let results = simulate(p, |comm: &SimRank| run(&case, Some(Arc::new(comm.clone()))));
        let wall = t.elapsed().as_secs_f64();
        let after = snapshot();
        let r = &results[0];
        println!(
            "{case} M2 np={p}: wall {wall:.3} s, dE {:.3e} Ha, (SCF, Davidson, H) {:?}",
            (r.total_hartree() - serial.total_hartree()).abs(),
            work(r)
        );
        let per = |a: f64, b: f64| (a - b) / p as f64;
        let line: Vec<String> = after
            .0
            .iter()
            .zip(&before.0)
            .map(|((n, a), (_, b))| format!("{n} {:.3}", per(*a, *b)))
            .collect();
        println!("  davidson stages per rank: {}", line.join(", "));
        let line: Vec<String> = after
            .1
            .iter()
            .zip(&before.1)
            .map(|((n, a), (_, b))| format!("{n} {:.3}", per(*a, *b)))
            .collect();
        println!("  m2 stages per rank: {}", line.join(", "));
        let line: Vec<String> = after
            .2
            .iter()
            .zip(&before.2)
            .map(|((n, c, a), (_, d, b))| format!("{n} {} calls {:.3} s", (c - d) / p, per(*a, *b)))
            .collect();
        println!("  collectives per rank: {}", line.join(", "));
    }
}
