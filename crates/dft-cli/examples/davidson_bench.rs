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
use dft_engine::snapshot::{ReplayContract, ReplayMode};
use dft_engine::{ScfProblem, ScfStep, SpeciesPotential, StandardProfile, parallel};
use dft_geometry::Atom;
use dft_mpi::Collective;
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
    /// Snapshot replay (BREAKTHROUGH-0-INSTRUMENT): the manifest written by
    /// `BAANDI_SNAPSHOT_DIR`; the timed kind is then `replay` (one map at the
    /// captured density from the captured start, kept H X armed when captured).
    snapshot: Option<PathBuf>,
    /// ISSUE-18: `exact` (default: operator, rows, potential, profile and
    /// kept status must match; a mismatch exits 3, an unreproduced solve
    /// exits 4), `ab` (algorithm A/B on the same physical problem: solver
    /// settings may differ, a changed potential rounding recomputes H X),
    /// `legacy` (manifest without fingerprint, unchecked, labeled).
    replay_mode: ReplayMode,
}

/// The parts of a snapshot manifest the replay needs.
struct Manifest {
    /// ISSUE-18: the captured contract (fingerprint, profile, potential).
    contract: ReplayContract,
    dir: PathBuf,
    label: String,
    stage: String,
    kind: String,
    iteration: u64,
    ranks: usize,
    layout: String,
    bands: usize,
    plane_waves: usize,
    grid: [usize; 3],
    cutoff_hartree: f64,
    tolerance: f64,
    start: String,
    kept_armed: bool,
    density_file: Option<String>,
    result_iterations: u64,
    result_h_columns: u64,
    result_eigen_fnv: String,
    result_kept_taken: bool,
}

fn read_manifest(path: &Path, mode: ReplayMode) -> Manifest {
    let text = std::fs::read_to_string(path).unwrap_or_else(|e| usage(&format!("--snapshot: {e}")));
    let v: Value =
        serde_json::from_str(&text).unwrap_or_else(|e| usage(&format!("--snapshot: {e}")));
    if v["schema"] != "baandi-snapshot-1" {
        usage("--snapshot: not a baandi-snapshot-1 manifest");
    }
    let contract = ReplayContract::from_manifest(&text, mode).unwrap_or_else(|e| {
        eprintln!("davidson_bench: replay identity check FAILED: --snapshot: {e}");
        std::process::exit(3)
    });
    let s = |k: &str| v[k].as_str().unwrap_or("").to_string();
    let u = |k: &str| v[k].as_u64().unwrap_or(0);
    let g = v["grid"].as_array().cloned().unwrap_or_default();
    let grid = [0, 1, 2].map(|i| g.get(i).and_then(Value::as_u64).unwrap_or(0) as usize);
    Manifest {
        contract,
        dir: path.parent().unwrap_or(Path::new(".")).to_path_buf(),
        label: s("label"),
        stage: s("stage"),
        kind: s("kind"),
        iteration: u("iteration"),
        ranks: u("ranks") as usize,
        layout: s("layout"),
        bands: u("bands") as usize,
        plane_waves: u("plane_waves") as usize,
        grid,
        cutoff_hartree: v["cutoff_hartree"].as_f64().unwrap_or(0.0),
        tolerance: v["tolerance_ha"].as_f64().unwrap_or(0.0),
        start: s("start"),
        kept_armed: v["kept_armed"].as_bool().unwrap_or(false),
        density_file: v["density_file"].as_str().map(str::to_string),
        result_iterations: v["result"]["davidson_iterations"].as_u64().unwrap_or(0),
        result_h_columns: v["result"]["h_columns"].as_u64().unwrap_or(0),
        result_eigen_fnv: v["result"]["eigen_fnv"].as_str().unwrap_or("").to_string(),
        result_kept_taken: v["result"]["kept_taken"].as_bool().unwrap_or(false),
    }
}

fn usage(msg: &str) -> ! {
    eprintln!(
        "davidson_bench: {msg}\nusage: davidson_bench --input FILE.toml [--density density.bin] [--repeats N] [--kinds cold,warm] [--sim-ranks P] [--prev-tol X] [--mix X] [--cold-tol X] [--warm-tol X] [--cold-start once|internal] [--snapshot MANIFEST.json [--replay-mode exact|ab|legacy]]"
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
        snapshot: None,
        replay_mode: ReplayMode::Exact,
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
            "--snapshot" => a.snapshot = Some(v().into()),
            "--replay-mode" => {
                a.replay_mode = match v().as_str() {
                    "exact" => ReplayMode::Exact,
                    "ab" => ReplayMode::AbSamePhysics,
                    "legacy" => ReplayMode::Legacy,
                    _ => usage("--replay-mode is exact, ab or legacy"),
                }
            }
            _ => usage(&format!("unknown argument {k}")),
        }
    }
    if a.input.as_os_str().is_empty() {
        usage("--input is required");
    }
    if a.snapshot.is_some() {
        a.kinds = vec!["replay".into()];
    }
    if a.repeats == 0 || a.sim_ranks == 0 || !(a.mix > 0.0 && a.mix <= 1.0) {
        usage("repeats and sim-ranks must be positive, mix in (0, 1]");
    }
    if a.kinds
        .iter()
        .any(|k| k != "cold" && k != "warm" && k != "replay")
    {
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
    /// Engine and operator M2 stages (`parallel::m2_stats`: density, XC, solve
    /// interface, FFT transposes, nonlocal, ...).
    m2_engine_seconds: Vec<(String, f64)>,
    /// Engine and operator M2 counters (the `lazy_*` ones are reported).
    m2_counts: Vec<(String, usize)>,
    reductions: (u64, u64),
    reduced: dft_solvers::distributed_davidson::ReducedSolveStats,
    rotation_flops: u64,
    start_gram: (u64, u64),
    lean_algebra: dft_solvers::distributed_davidson::LeanAlgebraStats,
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
        m2_engine_seconds: parallel::m2_stats().1,
        m2_counts: parallel::m2_stats().0,
        reductions: (
            red.1.iter().map(|s| s.calls).sum(),
            red.1.iter().map(|s| s.payload).sum(),
        ),
        reduced: parallel::reduced_solve_stats(),
        rotation_flops: dft_solvers::distributed_davidson::rotation_stats().flops,
        start_gram: dft_solvers::distributed_davidson::start_gram_stats(),
        lean_algebra: dft_solvers::distributed_davidson::lean_algebra_stats(),
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
    put(
        "lean.keep_contractions",
        l1.keep_contractions - l0.keep_contractions,
    );
    put(
        "lean.kept_directions",
        l1.kept_directions - l0.kept_directions,
    );
    put("lean.keep_dropped", l1.keep_dropped - l0.keep_dropped);
    put("lean.keep_fallbacks", l1.keep_fallbacks - l0.keep_fallbacks);
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
    put("lean.fused_borders", l1.fused_borders - l0.fused_borders);
    put(
        "lean.fused_wasted_columns",
        l1.fused_wasted_columns - l0.fused_wasted_columns,
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
    for ((n, s1), (_, s0)) in b.m2_engine_seconds.iter().zip(&a.m2_engine_seconds) {
        s.insert(format!("m2.{n}"), (s1 - s0) / p as f64);
    }
    // MPI-LAZY-ORBITALS counters (`m2.lazy_*`).
    for ((n, c1), (_, c0)) in b.m2_counts.iter().zip(&a.m2_counts) {
        if n.starts_with("lazy_") {
            put(&format!("m2.{n}"), (c1 - c0) as u64);
        }
    }
    put("dist.reductions.calls", b.reductions.0 - a.reductions.0);
    put("dist.reductions.doubles", b.reductions.1 - a.reductions.1);
    put("dist.reduced_solves", b.reduced.solves - a.reduced.solves);
    put(
        "dist.reduced_sum_m3",
        (b.reduced.sum_cubed - a.reduced.sum_cubed).max(0.0) as u64,
    );
    put("dist.rotation_flop", b.rotation_flops - a.rotation_flops);
    put("dist.start_gram_taken", b.start_gram.0 - a.start_gram.0);
    put("dist.start_gram_declined", b.start_gram.1 - a.start_gram.1);
    // MPI-OVERLAP-DAVIDSON (`BAANDI_M2_LEAN_ALGEBRA=on`).
    let (la, lb) = (&a.lean_algebra, &b.lean_algebra);
    put("dist.lean_solves", lb.solves - la.solves);
    put(
        "dist.lean_border_reductions",
        lb.border_reductions - la.border_reductions,
    );
    put("dist.lean_appends", lb.appends - la.appends);
    put(
        "dist.lean_reprojections",
        lb.reprojections - la.reprojections,
    );
    put("dist.lean_fallbacks", lb.fallbacks - la.fallbacks);
    put("dist.lean_rebuilds", lb.rebuilds - la.rebuilds);
    put(
        "dist.lean_wasted_columns",
        lb.wasted_columns - la.wasted_columns,
    );
    put(
        "dist.lean_factored_starts",
        lb.factored_starts - la.factored_starts,
    );
    // BREAKTHROUGH-1-RESTART R1-b (`BAANDI_RESTART_KEEP`).
    put("dist.lean_contractions", lb.contractions - la.contractions);
    put(
        "dist.keep_contractions",
        lb.keep_contractions - la.keep_contractions,
    );
    put(
        "dist.kept_directions",
        lb.kept_directions - la.kept_directions,
    );
    put("dist.keep_dropped", lb.keep_dropped - la.keep_dropped);
    put("dist.keep_fallbacks", lb.keep_fallbacks - la.keep_fallbacks);
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
    /// Replay mode: the manifest summary and the captured solve's counters.
    snapshot: Option<Value>,
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
    let mut opts = dft_cli::run::options(&config).expect("options");
    let manifest = args
        .snapshot
        .as_deref()
        .map(|p| read_manifest(p, args.replay_mode));
    if let Some(m) = &manifest {
        // A coarse-stage snapshot: the ladder's coarse problem is the same
        // input at the fraction of the cutoff without a fixed grid.
        if (m.cutoff_hartree - opts.cutoff_hartree).abs() > 1e-12 * opts.cutoff_hartree.abs() {
            opts.cutoff_hartree = m.cutoff_hartree;
            opts.grid_override = None;
        }
        let ranks = comm.map_or(1, |c| c.size());
        if m.ranks != ranks {
            usage(&format!(
                "--snapshot was captured with {} ranks, this run has {ranks}",
                m.ranks
            ));
        }
    }
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
    let construct_s = t0.elapsed().as_secs_f64();
    dft_engine::profile::add_init(construct_s);
    if let Some(m) = &manifest {
        return replay(args, m, &mut problem, comm, p, dims, electrons, construct_s);
    }
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
    // Warm starts of the session maps: blocks plus their lazy rows
    // (MPI-LAZY-ORBITALS; with `BAANDI_M2_LAZY` off the blocks are full and the
    // maps are `step_kpoints_with_certification`).
    let warm_of = |step: &ScfStep| {
        let (warm, lazy) = dft_engine::split_warm(Some(&step.k_points));
        (warm.expect("one block"), lazy)
    };
    // Atomic-start vectors of a cold map at rho0, computed once (the SCF's own
    // routine, `ScfProblem::atomic_start_vectors`); None without atomic orbitals
    // (then a cold map starts from the solver's own guess).
    let t_atomic = Instant::now();
    let atomic = problem.atomic_start_vectors(&rho0).expect("atomic start");
    let atomic_s = t_atomic.elapsed().as_secs_f64();
    // Setup map: X_prev and rho_out at rho0, from the atomic start.
    let setup = problem
        .step_kpoints_lazy(&rho0, atomic.as_deref(), &[], args.prev_tol, false)
        .expect("setup map");
    let gamma_real = problem.gamma_real_diagnostics().all(|d| d.gamma_real_used);
    let (x_prev, x_prev_lazy) = warm_of(&setup);
    let plane_waves = x_prev_lazy
        .first()
        .and_then(|l| l.as_ref())
        .map_or(x_prev[0].rows(), |l| l.global_len());
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
                    .step_kpoints_lazy(&rho0, Some(&x_prev), &x_prev_lazy, args.cold_tol, false)
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
                Some(warm_of(&rearm))
            } else if args.cold_internal {
                None
            } else {
                atomic.clone().map(|a| (a, Vec::new()))
            };
            let density = if kind == "warm" { &rho_w } else { &rho0 };
            let (initial, initial_lazy) = match &initial {
                Some((w, l)) => (Some(w.as_slice()), l.as_slice()),
                None => (None, &[][..]),
            };
            barrier();
            let before = snap();
            let t = Instant::now();
            let step: ScfStep = problem
                .step_kpoints_lazy(density, initial, initial_lazy, tolerance, false)
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
        snapshot: None,
    }
}

/// Snapshot replay: `--repeats` timed maps (plus one discarded warm-up) at the
/// captured density from the captured start, each injected afresh.
#[allow(clippy::too_many_arguments)]
fn replay(
    args: &Args,
    m: &Manifest,
    problem: &mut ScfProblem,
    comm: Option<&dyn Collective>,
    p: u64,
    dims: [usize; 3],
    electrons: f64,
    construct_s: f64,
) -> RankOut {
    let rank0 = comm.is_none_or(|c| c.rank() == 0);
    let rank = comm.map_or(0, |c| c.rank());
    let barrier = || {
        if let Some(c) = comm {
            c.barrier().expect("barrier");
        }
    };
    if dims != m.grid {
        usage(&format!(
            "--snapshot grid {:?} differs from the problem grid {dims:?}",
            m.grid
        ));
    }
    let density_file = m
        .density_file
        .as_deref()
        .unwrap_or_else(|| usage("--snapshot manifest has no density file"));
    let rho = dft_engine::snapshot::read_density_file(&m.dir.join(density_file), dims)
        .unwrap_or_else(|e| usage(&e));
    let (x, hx, rows_fnv, layout) = if m.start == "none" {
        (None, None, 0, m.layout.clone())
    } else {
        let path = m.dir.join(format!("{}.rank{rank}.bin", m.label));
        let (x, hx, f, l) =
            dft_engine::snapshot::read_rank_file(&path).unwrap_or_else(|e| usage(&e));
        (Some(x), hx, f, l)
    };
    if let Some(x) = &x
        && x.cols() != m.bands
    {
        usage(&format!(
            "--snapshot start has {} columns, manifest bands {}",
            x.cols(),
            m.bands
        ));
    }
    if rank0 {
        eprintln!(
            "davidson_bench: replay {} ({} {} iteration {}): start {} (kept armed {}, H X {}), layout {layout}, rows fnv {rows_fnv:#018x}, tolerance {:.3e} Ha; captured solve: Davidson {} H columns {} kept taken {}",
            m.label,
            m.stage,
            m.kind,
            m.iteration,
            m.start,
            m.kept_armed,
            hx.is_some(),
            m.tolerance,
            m.result_iterations,
            m.result_h_columns,
            m.result_kept_taken
        );
    }
    // ISSUE-18: this rank's contract (manifest plus its rank-file header).
    let contract = m.contract.clone().with_rank_file(rows_fnv, &layout);
    let identity_failed = |e: &dyn std::fmt::Display| -> ! {
        eprintln!("davidson_bench: replay identity check FAILED (rank {rank}): {e}");
        std::process::exit(3)
    };
    // A cold snapshot has no injected start: the operator fingerprint and the
    // density bytes are still checked (agreed over ranks) before any solve.
    if x.is_none() && contract.mode != ReplayMode::Legacy {
        let current = problem.operator_fingerprint(&rho);
        let mut problems = contract
            .operator
            .as_ref()
            .map_or(vec!["operator fingerprint absent".into()], |c| {
                c.differences(&current)
            });
        let density = dft_engine::snapshot::fnv(rho.as_slice().unwrap_or(&[]));
        if contract.density_fnv != Some(density) {
            problems.push(format!(
                "density bytes {}",
                dft_engine::snapshot::hex(density)
            ));
        }
        let mut any = [f64::from(u8::from(!problems.is_empty()))];
        if let Some(c) = comm {
            let mut out = [0.0];
            c.tree_sum_f64(&any, &mut out).expect("identity agreement");
            any = out;
        }
        if any[0] > 0.0 {
            identity_failed(&if problems.is_empty() {
                "mismatch on another rank".to_string()
            } else {
                problems.join("; ")
            });
        }
    }
    let mut replay_report: Option<dft_engine::snapshot::ReplayReport> = None;
    let setup_s = construct_s;
    let mut repeats = Vec::new();
    for r in 0..=args.repeats {
        barrier();
        let before = snap();
        let t = Instant::now();
        let step: ScfStep = match &x {
            Some(x) => {
                let (step, rep) = problem
                    .replay_checked(
                        &rho,
                        dft_engine::snapshot::ReplayStart {
                            x: x.clone(),
                            hx: hx.clone(),
                        },
                        &contract,
                    )
                    .unwrap_or_else(|e| identity_failed(&e));
                replay_report = Some(rep);
                step
            }
            None => problem
                .step_kpoints_lazy(&rho, None, &[], m.tolerance, false)
                .expect("replay cold map"),
        };
        let map_s = t.elapsed().as_secs_f64();
        barrier();
        let after = snap();
        if r == 0 {
            continue;
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
    let bands = m.bands;
    RankOut {
        setup_s,
        base_density: format!("snapshot {density_file}"),
        electrons,
        bands,
        plane_waves: m.plane_waves,
        grid: dims,
        gamma_real: problem.gamma_real_diagnostics().all(|d| d.gamma_real_used),
        warm_tol: m.tolerance,
        e_h_change: f64::NAN,
        kinds: vec![(
            "replay".to_string(),
            KindOut {
                tolerance: m.tolerance,
                repeats,
            },
        )],
        snapshot: Some(json!({
            "manifest": args.snapshot.as_ref().map(|p| p.display().to_string()),
            "label": m.label, "stage": m.stage, "kind": m.kind, "iteration": m.iteration,
            "start": m.start, "kept_armed": m.kept_armed, "layout": layout,
            "replay_mode": contract.mode.name(),
            "hx_injected": replay_report.as_ref().map(|r| r.hx_injected),
            "hx_recomputed": replay_report.as_ref().map(|r| r.hx_recomputed),
            "replay_notes": replay_report.as_ref().map(|r| r.notes.clone()),
            "captured": {"davidson_iterations": m.result_iterations, "h_columns": m.result_h_columns,
                "eigen_fnv": m.result_eigen_fnv, "kept_taken": m.result_kept_taken},
        })),
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
    profile: Option<Value>,
) {
    let mut kinds = Vec::new();
    let mut unreproduced = false;
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
        if let Some(s) = &out.snapshot {
            let c = &s["captured"];
            let work_reproduced = c["davidson_iterations"].as_u64()
                == Some(first.iterations as u64)
                && c["h_columns"].as_u64() == Some(first.h_columns as u64);
            let eigen_reproduced =
                c["eigen_fnv"].as_str() == Some(&format!("{:#018x}", first.eig_fnv));
            // ISSUE-18: an exact replay that does not reproduce its captured
            // contract fails the run (after the report is printed).
            if s["replay_mode"] == "exact" && !(work_reproduced && eigen_reproduced) {
                unreproduced = true;
            }
            println!(
                "  replay of {} ({} {} iteration {}, start {}): captured Davidson {} H columns {} eigenvalue FNV {}; work reproduced {work_reproduced}, eigenvalues reproduced {eigen_reproduced}",
                s["label"].as_str().unwrap_or(""),
                s["stage"].as_str().unwrap_or(""),
                s["kind"].as_str().unwrap_or(""),
                s["iteration"],
                s["start"].as_str().unwrap_or(""),
                c["davidson_iterations"],
                c["h_columns"],
                c["eigen_fnv"].as_str().unwrap_or("")
            );
        }
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
        "snapshot": out.snapshot,
        "phase_profile": profile,
    });
    println!("DAVBENCH_JSON {doc}");
    if unreproduced {
        eprintln!(
            "davidson_bench: exact replay did not reproduce the captured work/eigenvalues: FAILED"
        );
        std::process::exit(4);
    }
}

/// Phase profile of the process (`BAANDI_PHASE_PROFILE=on`): the per-rank
/// profiles gathered (real MPI) or taken from the registry (simulated ranks).
fn profile_report(profiles: Vec<dft_engine::profile::RankProfile>) -> Option<Value> {
    if profiles.is_empty() {
        return None;
    }
    let report = dft_engine::profile::Report::from_ranks(&profiles);
    print!("{}", report.table(None));
    Some(dft_cli::run::profile_json(&report, f64::NAN))
}

fn print_first_block() {
    let records = dft_engine::first_block::take_records();
    if !records.is_empty() {
        print!("{}", dft_engine::first_block::table(&records));
    }
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
            dft_engine::profile::finish_rank(comm.rank());
            parallel::set_communicator(None);
            parallel::set_m2_override(None);
            out
        });
        let refs: Vec<&RankOut> = outs.iter().collect();
        let diag = maxima(&refs, |r| r.diag_s);
        let map = maxima(&refs, |r| r.map_s);
        let profile = profile_report(dft_engine::profile::take_registry());
        print_first_block();
        report(&args, &outs[0], p, "simulated", &diag, &map, profile);
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
    let profiles = match dft_engine::profile::finish_rank(world.rank()) {
        Some(mine) => match &comm {
            Some(c) => dft_engine::profile::gather(c.as_ref(), &mine).expect("profile gather"),
            None => vec![mine],
        },
        None => Vec::new(),
    };
    parallel::set_communicator(None);
    if world.rank() == 0 {
        let profile = profile_report(profiles);
        print_first_block();
        report(
            &args,
            &out,
            size,
            if size > 1 { "mpi" } else { "threads" },
            &diag,
            &map,
            profile,
        );
    }
}
