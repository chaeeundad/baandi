//! Stage-0 benchmark of ONE standard-path Davidson solve at a fixed Hamiltonian
//! (ticket FASTBENCH-DAVIDSON, docs/performance/fast-loop-protocol.md, stage 0).
//!
//! The SCF convergence process (step counts, verification attempts, mixing) is
//! NOT measured here; this binary times the eigensolve of one SCF map at a fixed
//! input density, in-process and repeatably, and prints the exact work counters.
//!
//! Problem: a CLI input file (`--input`, the same TOML `dft-cli run` reads; the
//! harness `scripts/fast_screen_dav.py` writes it from a benchmarks/cases file with
//! the stage-1 displacement stream). Options are resolved by `dft_cli::run::options`
//! (the single standard path, Ritz-change solves, lean Davidson by default, the MPI
//! capacity rule), and the problem is built as `dft-cli run` builds it. The cutoff
//! ladder is not used: every solve runs at the full cutoff.
//!
//! Base density rho0: `--density FILE` (the headerless little-endian f64 C-order
//! `density.bin` a converged `dft-cli run` writes, docs/cli.md) or, without it, the
//! superposition of atomic densities of the problem (the SCF start density). The
//! JSON records which (`base_density`).
//!
//! Setup (untimed): one solve at rho0 from the atomic start at `--prev-tol`
//! (default 1e-8 Ha) gives X_prev and rho_out; the warm density is
//! rho_w = rho0 + mix (rho_out - rho0) (`--mix`, default 0.5, the default Pulay
//! beta); the warm tolerance is the standard ratchet after that map,
//! ratchet(prev_tol, E_H[rho_out - rho0], N_e) (`--warm-tol` overrides).
//!   cold: map at rho0 from the atomic start at the standard first-step tolerance
//!         `StandardProfile::CAP_HA` (5e-3 Ha) (`--cold-tol` overrides). The atomic
//!         start vectors are computed once by the SCF's routine
//!         (`ScfProblem::atomic_start_vectors`) and passed to every cold map, which
//!         then runs the identical solve (`--cold-start internal` recomputes them in
//!         every map, the SCF's own path, to check the identity; slower);
//!   warm: an untimed re-arm map at rho0 from X_prev at the cold tolerance (two
//!         Davidson iterations from converged vectors; it leaves X, H X and psi(r)
//!         kept for REUSE-PSI exactly as an SCF step does), then the timed map at
//!         rho_w from the re-arm orbitals: the typical mid-SCF solve.
//! Each kind runs one discarded warm-up and `--repeats` (default 5) timed repeats.
//! The time per solve is the diagonalization phase of the map (`PhaseTimings`),
//! i.e. the solver call only; the whole map (potential, density, audit) is also
//! reported. Counters are deltas over each timed map (process-wide counters; with
//! simulated ranks the sum over ranks divided by the rank count).
//!
//! Usage (threads, pure or `--features mkl`):
//!   cargo build --release -p dft-cli --example davidson_bench
//!   target/release/examples/davidson_bench --input in.toml [--density density.bin]
//!       [--repeats 5] [--kinds cold,warm] [--sim-ranks P] [--prev-tol 1e-8] [--mix 0.5]
//! Simulated ranks (threads of one process, the M2 stick-row path, counters only for
//! latency questions): `--sim-ranks 4`.
//! Real MPI (dev08): `cargo build --release -p dft-cli --features mpi --example
//! davidson_bench` and `mpirun -np 4 target/release/examples/davidson_bench --input in.toml`.
//! Output: human lines and one `DAVBENCH_JSON {...}` line (rank 0), read by
//! `scripts/fast_screen_dav.py`. Switches (`BAANDI_*`) are read once per process:
//! an A/B runs the binary twice. With `BAANDI_DAVIDSON_STAGES=on` the FFT pair
//! counts and stage seconds are filled (timer cost: never for the timed arms).
use dft_cli::config::Config;
use dft_engine::{ScfProblem, ScfStep, SpeciesPotential, StandardProfile, parallel};
use dft_geometry::Atom;
use dft_mpi::Collective;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::data::Array3;
use serde_json::{Value, json};
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
    sync::Arc,
    time::Instant,
};

type Shared = Arc<dyn Collective + Send + Sync>;

struct Args {
    input: PathBuf,
    density: Option<PathBuf>,
    repeats: usize,
    kinds: Vec<String>,
    sim_ranks: usize,
    prev_tol: f64,
    mix: f64,
    cold_tol: f64,
    warm_tol: Option<f64>,
    /// Recompute the atomic start inside every cold map (the SCF's own code path)
    /// instead of passing the vectors computed once (identical solve; slower).
    cold_internal: bool,
}

fn usage(msg: &str) -> ! {
    eprintln!(
        "davidson_bench: {msg}\nusage: davidson_bench --input FILE.toml [--density density.bin] [--repeats N] [--kinds cold,warm] [--sim-ranks P] [--prev-tol X] [--mix X] [--cold-tol X] [--warm-tol X] [--cold-start once|internal]"
    );
    std::process::exit(2)
}

fn parse_args() -> Args {
    let mut a = Args {
        input: PathBuf::new(),
        density: None,
        repeats: 5,
        kinds: vec!["cold".into(), "warm".into()],
        sim_ranks: 1,
        prev_tol: 1e-8,
        mix: 0.5,
        cold_tol: StandardProfile::CAP_HA,
        warm_tol: None,
        cold_internal: false,
    };
    let mut it = std::env::args().skip(1);
    while let Some(k) = it.next() {
        let mut v = || {
            it.next()
                .unwrap_or_else(|| usage(&format!("{k} needs a value")))
        };
        let num = |s: String| -> f64 { s.parse().unwrap_or_else(|_| usage("bad number")) };
        match k.as_str() {
            "--input" => a.input = v().into(),
            "--density" => a.density = Some(v().into()),
            "--repeats" => a.repeats = num(v()) as usize,
            "--kinds" => a.kinds = v().split(',').map(str::to_string).collect(),
            "--sim-ranks" => a.sim_ranks = num(v()) as usize,
            "--prev-tol" => a.prev_tol = num(v()),
            "--mix" => a.mix = num(v()),
            "--cold-tol" => a.cold_tol = num(v()),
            "--warm-tol" => a.warm_tol = Some(num(v())),
            "--cold-start" => {
                a.cold_internal = match v().as_str() {
                    "once" => false,
                    "internal" => true,
                    _ => usage("--cold-start is once or internal"),
                }
            }
            _ => usage(&format!("unknown argument {k}")),
        }
    }
    if a.input.as_os_str().is_empty() {
        usage("--input is required");
    }
    if a.repeats == 0 || a.sim_ranks == 0 || !(a.mix > 0.0 && a.mix <= 1.0) {
        usage("repeats and sim-ranks must be positive, mix in (0, 1]");
    }
    if a.kinds.iter().any(|k| k != "cold" && k != "warm") {
        usage("kinds are cold and warm");
    }
    a
}

/// Process-wide counters at one instant.
#[derive(Clone)]
struct Snap {
    gemm: BTreeMap<String, (u64, u64)>,
    flops: dft_solvers::davidson_stages::FlopCounters,
    lean: dft_solvers::davidson_stages::LeanCounters,
    work: dft_numerics::work_timers::WorkTimes,
    stages: dft_engine::StageSecondsReport,
    nl_audits: (u64, u64),
    reuse: dft_engine::ReusePsiCounters,
    coll_calls: Vec<(&'static str, usize, f64)>,
    coll_bytes: Vec<(&'static str, u64)>,
    m2_stage_calls: Vec<(&'static str, u64)>,
    m2_stage_seconds: Vec<(&'static str, f64)>,
    reductions: (u64, u64),
    reduced: dft_solvers::distributed_davidson::ReducedSolveStats,
    rotation_flops: u64,
}

fn snap() -> Snap {
    let mut gemm = BTreeMap::new();
    for r in dft_numerics::gemm_dispatch::rows() {
        gemm.insert(format!("{}/{}", r.site, r.route), (r.calls, r.flops));
    }
    let w = dft_numerics::work_timers::snapshot();
    let work = w.iter().skip(1).fold(w[0], |acc, x| acc.plus(x));
    let red = parallel::reduction_stats();
    Snap {
        gemm,
        flops: dft_solvers::davidson_stages::flop_counters(),
        lean: dft_solvers::davidson_stages::lean_counters(),
        work,
        stages: dft_engine::stage_seconds(),
        nl_audits: dft_engine::nl_audit_counters(),
        reuse: dft_engine::reuse_psi_counters(),
        coll_calls: parallel::collective_stats(),
        coll_bytes: parallel::collective_bytes(),
        m2_stage_calls: parallel::stage_calls(),
        m2_stage_seconds: parallel::stage_seconds(),
        reductions: (
            red.1.iter().map(|s| s.calls).sum(),
            red.1.iter().map(|s| s.payload).sum(),
        ),
        reduced: parallel::reduced_solve_stats(),
        rotation_flops: dft_solvers::distributed_davidson::rotation_stats().flops,
    }
}

/// Counter deltas `b - a` per rank (`p` = simulated ranks sharing the process).
/// Deterministic counters are integers; seconds go to `seconds`.
fn delta(a: &Snap, b: &Snap, p: u64) -> (BTreeMap<String, Value>, BTreeMap<String, f64>) {
    let mut c = BTreeMap::new();
    let mut s = BTreeMap::new();
    let mut put = |k: &str, v: u64| {
        c.insert(k.to_string(), json!(v / p));
    };
    let (mut calls, mut flops) = (0, 0);
    for (k, (n1, f1)) in &b.gemm {
        let (n0, f0) = a.gemm.get(k).copied().unwrap_or((0, 0));
        if n1 > &n0 {
            put(&format!("gemm.{k}.calls"), n1 - n0);
            put(&format!("gemm.{k}.flop"), f1 - f0);
            calls += n1 - n0;
            flops += f1 - f0;
        }
    }
    put("gemm_calls", calls);
    put("gemm_flop", flops);
    put("davidson_flop", b.flops.total() - a.flops.total());
    put(
        "reduced_eigen_full",
        b.flops.eigen_full - a.flops.eigen_full,
    );
    put(
        "reduced_eigen_full_m3",
        b.flops.eigen_full_m3 - a.flops.eigen_full_m3,
    );
    let (l0, l1) = (&a.lean, &b.lean);
    put("lean.iterations", l1.iterations - l0.iterations);
    put(
        "lean.residual_columns",
        l1.residual_columns - l0.residual_columns,
    );
    put(
        "lean.appended_columns",
        l1.appended_columns - l0.appended_columns,
    );
    put("lean.contractions", l1.contractions - l0.contractions);
    put("lean.fallbacks", l1.fallbacks - l0.fallbacks);
    put("lean.kept_starts", l1.kept_starts - l0.kept_starts);
    put("lean.kept_declined", l1.kept_declined - l0.kept_declined);
    put("lean.kept_columns", l1.kept_columns - l0.kept_columns);
    put(
        "lean.proj_cache_returns",
        l1.proj_cache_returns - l0.proj_cache_returns,
    );
    put(
        "lean.proj_cache_columns",
        l1.proj_cache_columns - l0.proj_cache_columns,
    );
    let (w0, w1) = (&a.work, &b.work);
    put(
        "fft.pair_transforms",
        w1.pair_transforms - w0.pair_transforms,
    );
    put(
        "fft.kept_transforms",
        w1.kept_transforms - w0.kept_transforms,
    );
    put("h_apply.calls", w1.apply_calls - w0.apply_calls);
    put("h_apply.columns", w1.apply_columns - w0.apply_columns);
    put("nl_audit.cached", b.nl_audits.0 - a.nl_audits.0);
    put("nl_audit.direct", b.nl_audits.1 - a.nl_audits.1);
    put("reuse.armed", b.reuse.armed - a.reuse.armed);
    put(
        "reuse.columns_saved",
        b.reuse.columns_saved - a.reuse.columns_saved,
    );
    for ((n, c1, s1), (_, c0, s0)) in b.coll_calls.iter().zip(&a.coll_calls) {
        put(&format!("coll.{n}.calls"), (c1 - c0) as u64);
        s.insert(format!("coll.{n}"), (s1 - s0) / p as f64);
    }
    for ((n, b1), (_, b0)) in b.coll_bytes.iter().zip(&a.coll_bytes) {
        put(&format!("coll.{n}.bytes"), b1 - b0);
    }
    for ((n, c1), (_, c0)) in b.m2_stage_calls.iter().zip(&a.m2_stage_calls) {
        put(&format!("dist.{n}.calls"), c1 - c0);
    }
    for ((n, s1), (_, s0)) in b.m2_stage_seconds.iter().zip(&a.m2_stage_seconds) {
        s.insert(format!("dist.{n}"), (s1 - s0) / p as f64);
    }
    put("dist.reductions.calls", b.reductions.0 - a.reductions.0);
    put("dist.reductions.doubles", b.reductions.1 - a.reductions.1);
    put("dist.reduced_solves", b.reduced.solves - a.reduced.solves);
    put(
        "dist.reduced_sum_m3",
        (b.reduced.sum_cubed - a.reduced.sum_cubed).max(0.0) as u64,
    );
    put("dist.rotation_flop", b.rotation_flops - a.rotation_flops);
    let (t0, t1) = (&a.stages, &b.stages);
    for (k, v) in [
        ("stage.reduced", t1.reduced - t0.reduced),
        ("stage.eigen", t1.eigen - t0.eigen),
        ("stage.rotate", t1.rotate - t0.rotate),
        ("stage.residual", t1.residual - t0.residual),
        (
            "stage.orthonormalize",
            t1.orthonormalize - t0.orthonormalize,
        ),
        ("stage.h_apply", t1.h_apply - t0.h_apply),
        ("stage.precondition", t1.precondition - t0.precondition),
        ("stage.local", w1.local - w0.local),
        (
            "stage.nonlocal",
            w1.nonlocal_projection - w0.nonlocal_projection + w1.nonlocal_back - w0.nonlocal_back,
        ),
    ] {
        s.insert(k.to_string(), v / p as f64);
    }
    (c, s)
}

fn fnv(values: &[f64]) -> u64 {
    values.iter().fold(0xcbf2_9ce4_8422_2325, |h: u64, x| {
        (h ^ x.to_bits()).wrapping_mul(0x0100_0000_01b3)
    })
}

/// E_H[a - b] (Ha) from the Hartree energies of the two nonnegative densities:
/// E[a - b] = E[a] + E[b] - integral a V_H[b].
fn hartree_difference(a: &Array3<f64>, b: &Array3<f64>, cell: &dft_geometry::Cell) -> f64 {
    let dims = [a.shape()[0], a.shape()[1], a.shape()[2]];
    let mut fft = dft_numerics::FftPlan3d::new(dims).expect("FFT plan");
    let ha = dft_physics::hartree::evaluate(a, cell, &mut fft).expect("Hartree a");
    let hb = dft_physics::hartree::evaluate(b, cell, &mut fft).expect("Hartree b");
    let cross: f64 = a
        .iter()
        .zip(&hb.potential)
        .map(|(x, v)| x * v.re)
        .sum::<f64>()
        * cell.volume_bohr3()
        / a.len() as f64;
    (ha.energy.0 + hb.energy.0 - cross).max(0.0)
}

struct Repeat {
    diag_s: f64,
    map_s: f64,
    iterations: usize,
    h_columns: usize,
    eig_fnv: u64,
    eigenvalues: Vec<f64>,
    residuals: Vec<f64>,
    converged: usize,
    counters: BTreeMap<String, Value>,
    seconds: BTreeMap<String, f64>,
}

struct KindOut {
    tolerance: f64,
    repeats: Vec<Repeat>,
}

struct RankOut {
    setup_s: f64,
    base_density: String,
    electrons: f64,
    bands: usize,
    plane_waves: usize,
    grid: [usize; 3],
    gamma_real: bool,
    warm_tol: f64,
    e_h_change: f64,
    kinds: Vec<(String, KindOut)>,
}

fn bench(args: &Args, comm: Option<&dyn Collective>, p: u64) -> RankOut {
    let t0 = Instant::now();
    let rank0 = comm.is_none_or(|c| c.rank() == 0);
    let barrier = || {
        if let Some(c) = comm {
            c.barrier().expect("barrier");
        }
    };
    let text = std::fs::read_to_string(&args.input).expect("read --input");
    let config = Config::parse(&text).expect("parse --input");
    let opts = dft_cli::run::options(&config).expect("options");
    let base_dir = args.input.parent().unwrap_or(Path::new("."));
    let mut potentials = Vec::new();
    for (name, path) in &config.pseudopotentials {
        let data = std::fs::read_to_string(base_dir.join(path)).expect("read UPF");
        potentials.push((name.clone(), dft_upf::parse_upf(&data).expect("parse UPF")));
    }
    let refs: Vec<_> = potentials
        .iter()
        .map(|(label, potential)| SpeciesPotential { label, potential })
        .collect();
    let atoms: Vec<_> = config
        .system
        .atoms
        .iter()
        .map(|a| Atom::new(&a.species, a.fractional))
        .collect();
    let prepared = config.prepare().expect("prepare");
    let points = config.resolved_points().expect("k points");
    if points.len() != 1 {
        usage("one k point (Gamma) expected");
    }
    let mut problem = if config.basis.legacy_gamma() {
        ScfProblem::new(&prepared.cell, &atoms, &refs, opts.clone())
    } else {
        ScfProblem::with_kpoints(
            &prepared.cell,
            &atoms,
            &refs,
            opts.clone(),
            &points,
            dft_engine::OccupationPolicy::Insulator {
                minimum_gap_hartree: 1e-7,
            },
        )
    }
    .expect("problem");
    if let Some(floor) = config
        .xc
        .density_floor
        .or_else(|| problem.core_density().map(|_| 1e-10))
    {
        problem.set_density_floor(floor).expect("density floor");
    }
    let dims = problem.grid_dimensions();
    let electrons = problem.electron_count();
    let (rho0, base_density) = match &args.density {
        Some(path) => {
            let bytes = std::fs::read(path).expect("read --density");
            let n = dims.iter().product::<usize>();
            if bytes.len() != 8 * n {
                usage(&format!(
                    "--density has {} bytes, grid {dims:?} needs {}",
                    bytes.len(),
                    8 * n
                ));
            }
            let v: Vec<f64> = bytes
                .as_chunks::<8>()
                .0
                .iter()
                .map(|c| f64::from_le_bytes(*c))
                .collect();
            (
                Array3::from_shape_vec((dims[0], dims[1], dims[2]), v).unwrap(),
                format!("file {}", path.display()),
            )
        }
        None => (
            problem.starting_density(),
            if problem.atomic_superposition_density().is_some() {
                "atomic superposition".to_string()
            } else {
                "uniform (no atomic charge in the UPF)".to_string()
            },
        ),
    };
    let one = |m: &ComplexMatrix| vec![m.clone()];
    let construct_s = t0.elapsed().as_secs_f64();
    // Atomic-start vectors of a cold map at rho0, computed once (the SCF's own
    // routine, `ScfProblem::atomic_start_vectors`); None without atomic orbitals
    // (then a cold map starts from the solver's own guess).
    let t_atomic = Instant::now();
    let atomic = problem.atomic_start_vectors(&rho0).expect("atomic start");
    let atomic_s = t_atomic.elapsed().as_secs_f64();
    // Setup map: X_prev and rho_out at rho0, from the atomic start.
    let setup = problem
        .step_kpoints_with_certification(&rho0, atomic.as_deref(), args.prev_tol, false)
        .expect("setup map");
    let gamma_real = problem.gamma_real_diagnostics().all(|d| d.gamma_real_used);
    let x_prev = one(&setup.k_points[0].solution.eigenvectors);
    let plane_waves = x_prev[0].rows();
    let bands = setup.k_points[0].solution.eigenvalues.len();
    let e_h_change = hartree_difference(&setup.density, &rho0, &prepared.cell);
    let warm_tol = args
        .warm_tol
        .unwrap_or_else(|| StandardProfile::ratchet(args.prev_tol, e_h_change, electrons));
    let mut rho_w = rho0.clone();
    for (w, o) in rho_w.iter_mut().zip(&setup.density) {
        *w += args.mix * (o - *w);
    }
    drop(setup);
    let setup_s = t0.elapsed().as_secs_f64();
    if rank0 {
        eprintln!(
            "davidson_bench: construct {construct_s:.3} s, atomic start {atomic_s:.3} s, setup map {:.3} s",
            setup_s - construct_s - atomic_s
        );
    }
    let mut kinds = Vec::new();
    for kind in &args.kinds {
        let tolerance = if kind == "cold" {
            args.cold_tol
        } else {
            warm_tol
        };
        let mut repeats = Vec::new();
        for r in 0..=args.repeats {
            let initial = if kind == "warm" {
                let rearm = problem
                    .step_kpoints_with_certification(&rho0, Some(&x_prev), args.cold_tol, false)
                    .expect("re-arm map");
                if r == 0 && rank0 {
                    let s = &rearm.k_points[0].solution;
                    eprintln!(
                        "davidson_bench: re-arm map {:.3} s (solve {:.3} s), Davidson {} H columns {}",
                        rearm.timings.potential.as_secs_f64()
                            + rearm.timings.diagonalization.as_secs_f64()
                            + rearm.timings.density.as_secs_f64()
                            + rearm.timings.energy.as_secs_f64(),
                        rearm.timings.diagonalization.as_secs_f64(),
                        s.iterations,
                        s.operator_applications
                    );
                }
                Some(one(&rearm.k_points[0].solution.eigenvectors))
            } else if args.cold_internal {
                None
            } else {
                atomic.clone()
            };
            let density = if kind == "warm" { &rho_w } else { &rho0 };
            barrier();
            let before = snap();
            let t = Instant::now();
            let step: ScfStep = problem
                .step_kpoints_with_certification(density, initial.as_deref(), tolerance, false)
                .expect("timed map");
            let map_s = t.elapsed().as_secs_f64();
            barrier();
            let after = snap();
            if r == 0 {
                continue; // warm-up
            }
            let s = &step.k_points[0].solution;
            let (counters, seconds) = delta(&before, &after, p);
            repeats.push(Repeat {
                diag_s: step.timings.diagonalization.as_secs_f64(),
                map_s,
                iterations: s.iterations,
                h_columns: s.operator_applications,
                eig_fnv: fnv(&s.eigenvalues),
                eigenvalues: s.eigenvalues.clone(),
                residuals: s.residuals.clone(),
                converged: s.converged.iter().filter(|&&c| c).count(),
                counters,
                seconds,
            });
        }
        kinds.push((kind.clone(), KindOut { tolerance, repeats }));
    }
    RankOut {
        setup_s,
        base_density,
        electrons,
        bands,
        plane_waves,
        grid: dims,
        gamma_real,
        warm_tol,
        e_h_change,
        kinds,
    }
}

fn median(v: &mut [f64]) -> f64 {
    v.sort_by(f64::total_cmp);
    let n = v.len();
    if n % 2 == 1 {
        v[n / 2]
    } else {
        0.5 * (v[n / 2 - 1] + v[n / 2])
    }
}

/// Rank 0 report; `diag[k][r]` / `map[k][r]` are the maxima over ranks.
fn report(
    args: &Args,
    out: &RankOut,
    ranks: usize,
    mode: &str,
    diag: &[Vec<f64>],
    map: &[Vec<f64>],
) {
    let mut kinds = Vec::new();
    println!(
        "davidson_bench {mode} ranks {ranks}: bands {} plane waves {} grid {:?} electrons {} base density {} gamma_real {}; setup {:.3} s; E_H[rho_out - rho0] {:.3e} Ha, warm tolerance {:.3e} Ha",
        out.bands,
        out.plane_waves,
        out.grid,
        out.electrons,
        out.base_density,
        out.gamma_real,
        out.setup_s,
        out.e_h_change,
        out.warm_tol
    );
    for (i, (kind, k)) in out.kinds.iter().enumerate() {
        let mut d: Vec<f64> = diag[i].iter().map(|x| x * 1e3).collect();
        let mut m: Vec<f64> = map[i].iter().map(|x| x * 1e3).collect();
        let (dmin, dmax) = (
            d.iter().cloned().fold(f64::INFINITY, f64::min),
            d.iter().cloned().fold(0.0, f64::max),
        );
        let dmed = median(&mut d);
        let mmed = median(&mut m);
        let first = &k.repeats[0];
        let deterministic = k.repeats.iter().all(|r| {
            r.counters == first.counters
                && r.iterations == first.iterations
                && r.h_columns == first.h_columns
        });
        let identical = k.repeats.iter().all(|r| r.eig_fnv == first.eig_fnv);
        let occupied = (out.electrons / 2.0).round() as usize;
        println!(
            "  {kind}: tol {:.3e} Ha, solve median {dmed:.2} ms (min {dmin:.2} max {dmax:.2}, {} repeats), map median {mmed:.2} ms; Davidson {} H columns {} converged {}/{}; eigenvalue FNV {:#018x} (repeats identical {identical}); counters deterministic {deterministic}",
            k.tolerance,
            d.len(),
            first.iterations,
            first.h_columns,
            first.converged,
            first.eigenvalues.len(),
            first.eig_fnv
        );
        kinds.push(json!({
            "kind": kind,
            "tolerance_ha": k.tolerance,
            "solve_ms": diag[i].iter().map(|x| x * 1e3).collect::<Vec<_>>(),
            "map_ms": map[i].iter().map(|x| x * 1e3).collect::<Vec<_>>(),
            "solve_ms_median": dmed,
            "map_ms_median": mmed,
            "davidson_iterations": first.iterations,
            "h_columns": first.h_columns,
            "converged": first.converged,
            "eigen_fnv": format!("{:#018x}", first.eig_fnv),
            "repeats_identical": identical,
            "counters_deterministic": deterministic,
            "eigenvalues": first.eigenvalues,
            "band_energy_occupied": first.eigenvalues.iter().take(occupied).sum::<f64>(),
            "residuals": first.residuals,
            "residual_max": first.residuals.iter().cloned().fold(0.0, f64::max),
            "counters": first.counters,
            "seconds": first.seconds,
        }));
    }
    let doc = json!({
        "schema": "davbench-1",
        "mode": mode,
        "ranks": ranks,
        "input": args.input.display().to_string(),
        "base_density": out.base_density,
        "bands": out.bands,
        "occupied": (out.electrons / 2.0).round() as usize,
        "plane_waves": out.plane_waves,
        "grid": out.grid,
        "gamma_real": out.gamma_real,
        "setup_s": out.setup_s,
        "prev_tol_ha": args.prev_tol,
        "mix": args.mix,
        "e_h_change_ha": out.e_h_change,
        "warm_tol_ha": out.warm_tol,
        "repeats": args.repeats,
        "mkl_backend": dft_numerics::MKL_BACKEND,
        "kinds": kinds,
    });
    println!("DAVBENCH_JSON {doc}");
}

fn main() {
    let args = parse_args();
    let universe = dft_mpi::Universe::init().expect("MPI init");
    let world = universe.world();
    if dft_numerics::switches::choice(dft_numerics::switches::Id::MPI_TREE) == "rd" {
        dft_mpi::set_tree_sum(dft_mpi::TreeSum::RecursiveDoubling);
    }
    let maxima = |outs: &[&RankOut], f: fn(&Repeat) -> f64| -> Vec<Vec<f64>> {
        (0..outs[0].kinds.len())
            .map(|k| {
                (0..outs[0].kinds[k].1.repeats.len())
                    .map(|r| {
                        outs.iter()
                            .map(|o| f(&o.kinds[k].1.repeats[r]))
                            .fold(0.0, f64::max)
                    })
                    .collect()
            })
            .collect()
    };
    if args.sim_ranks > 1 {
        if world.size() > 1 {
            usage("--sim-ranks needs a single-process run");
        }
        let p = args.sim_ranks;
        // The CLI capacity rule reads the rank count of the process.
        dft_cli::parallel::set(0, p);
        let outs = dft_mpi::simulate(p, |comm| {
            parallel::set_m2_override(Some(true));
            parallel::set_communicator(Some(Arc::new(comm.clone()) as Shared));
            let out = bench(&args, Some(comm as &dyn Collective), p as u64);
            parallel::set_communicator(None);
            parallel::set_m2_override(None);
            out
        });
        let refs: Vec<&RankOut> = outs.iter().collect();
        let diag = maxima(&refs, |r| r.diag_s);
        let map = maxima(&refs, |r| r.map_s);
        report(&args, &outs[0], p, "simulated", &diag, &map);
        return;
    }
    let size = world.size();
    dft_cli::parallel::set(world.rank(), size);
    let comm: Option<Shared> = (size > 1).then(|| Arc::new(world) as Shared);
    parallel::set_communicator(comm.clone());
    let out = bench(&args, comm.as_deref().map(|c| c as &dyn Collective), 1);
    // Maxima over ranks of every timed repeat (one allgatherv per series).
    let gather = |f: fn(&Repeat) -> f64| -> Vec<Vec<f64>> {
        out.kinds
            .iter()
            .map(|(_, k)| {
                let local: Vec<f64> = k.repeats.iter().map(f).collect();
                match &comm {
                    None => local,
                    Some(c) => {
                        let n = local.len();
                        let mut all = vec![0.0; n * size];
                        c.allgatherv_f64(&local, &vec![n; size], &mut all)
                            .expect("allgatherv");
                        (0..n)
                            .map(|r| (0..size).map(|q| all[q * n + r]).fold(0.0, f64::max))
                            .collect()
                    }
                }
            })
            .collect()
    };
    let diag = gather(|r| r.diag_s);
    let map = gather(|r| r.map_s);
    parallel::set_communicator(None);
    if world.rank() == 0 {
        report(
            &args,
            &out,
            size,
            if size > 1 { "mpi" } else { "threads" },
            &diag,
            &map,
        );
    }
}
