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
//! Optional: M2_PROFILE_CASE=si8, si16, si32 or si64 (default si8), M2_PROFILE_RANKS=2 (comma list),
//! M2_PROFILE_CAPACITY=<2..12> (Davidson capacity per band, default 6),
//! M2_PROFILE_SERIAL=off (skip the threads-only reference run).
//! M2_PROFILE_STANDARD=on runs the standard (Ritz-change) profile instead of strict.
//! M2_PROFILE_MODE=production runs the production default (`ScfOptions::default()`:
//! standard 5e-9 Ha, fixed-point verification, Pulay 0.5) with the CLI capacity rule
//! (3 per band at four or more ranks unless `BAANDI_M2_CAPACITY=off` or
//! M2_PROFILE_CAPACITY is set); combine with M2_PROFILE_LADDER=on for the CLI default;
//! M2_PROFILE_VERIFY=off drops the verification (to size its cost only).
//! Every distributed run also prints one machine-readable line `M2SCREEN_JSON {...}`
//! (per-rank counters and stages) read by `scripts/fast_screen_mpi.py`.
//! The Ritz rotation work (`rotation_stats`: pairs, columns, flops) is printed per rank.
//! The reduced-eigenproblem distribution (orders, sum of m^3, single-column
//! appends, bordered and broadcast solves) is printed per rank
//! (docs/performance/mpi-reduced-eigen-analysis.md).
//! M2_PROFILE_LADDER=on runs the C1 cutoff ladder (coarse stage report, certificate,
//! screen and atomic-start clocks; docs/performance/mpi-coarse-certificate-analysis.md);
//! the energy bits and a density FNV-1a hash are printed for A/B identity checks.
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
    let (rx, ry, rz) = match case {
        "si16" => (2, 1, 1),
        "si32" => (2, 2, 1),
        "si64" => (2, 2, 2),
        _ => (1, 1, 1),
    };
    let cell = Cell::from_user_rows_bohr([
        [A * rx as f64, 0., 0.],
        [0., A * ry as f64, 0.],
        [0., 0., A * rz as f64],
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
    for k in 0..rz {
        for i in 0..rx {
            for j in 0..ry {
                for p in BASE {
                    let q = [
                        (p[0] + i as f64) / rx as f64 + 1e-3 * next(),
                        (p[1] + j as f64) / ry as f64 + 1e-3 * next(),
                        (p[2] + k as f64) / rz as f64 + 1e-3 * next(),
                    ];
                    atoms.push(Atom::new("Si", q));
                }
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
    let ranks = comm.as_ref().map_or(1, |c| c.size());
    let production = std::env::var("M2_PROFILE_MODE").as_deref() == Ok("production");
    let capacity = std::env::var("M2_PROFILE_CAPACITY")
        .ok()
        .and_then(|v| v.parse().ok())
        .unwrap_or_else(|| {
            // The CLI rule (MPI-CAPACITY-RULE) in production mode, else the default 6.
            let auto =
                dft_numerics::switches::choice(dft_numerics::switches::Id::M2_CAPACITY) == "auto";
            if production && ranks >= 4 && auto {
                3
            } else {
                6
            }
        });
    let base = if production {
        ScfOptions {
            // M2_PROFILE_VERIFY=off: production without the fixed-point
            // verification (sizes the verification cost; never a production arm).
            verification: if std::env::var("M2_PROFILE_VERIFY").as_deref() == Ok("off") {
                dft_engine::Verification::Off
            } else {
                dft_engine::Verification::FixedPoint
            },
            ..ScfOptions::default()
        }
    } else {
        ScfOptions {
            standard: (std::env::var("M2_PROFILE_STANDARD").as_deref() == Ok("on"))
                .then(dft_engine::StandardProfile::default),
            ..ScfOptions::strict_reference()
        }
    };
    if m2 {
        parallel::set_m2_override(Some(true));
    }
    parallel::set_communicator(comm);
    let potentials = [SpeciesPotential {
        label: "Si",
        potential: &pp,
    }];
    let t_construct = Instant::now();
    let problem = ScfProblem::new(
        &cell,
        &atoms,
        &potentials,
        ScfOptions {
            gamma_real: GammaReal::Auto,
            cutoff_hartree: 15.,
            threads: 1,
            davidson_capacity_per_band: capacity,
            ..base
        },
    )
    .unwrap();
    // BREAKTHROUGH-0-INSTRUMENT: construction is the `init` bucket of the
    // phase profile (BAANDI_PHASE_PROFILE=on; otherwise a no-op).
    dft_engine::profile::add_init(t_construct.elapsed().as_secs_f64());
    let result = if std::env::var("M2_PROFILE_LADDER").as_deref() == Ok("on") {
        let out = dft_engine::cutoff_ladder::run_ladder(
            problem,
            &potentials,
            &dft_engine::cutoff_ladder::CutoffLadderOptions::default(),
            false,
        )
        .unwrap();
        let r = &out.report;
        println!(
            "  ladder: used {} coarse grid {:?} pw {} iterations {} Davidson {} setup {:.3} s session {:.3} s prolong {:.3} s, coarse phases potential {:.3} diag {:.3} density {:.3} energy {:.3} mixing {:.3}, atomic start {:.3} s",
            r.used,
            r.coarse_grid,
            r.coarse_plane_waves,
            r.coarse_iterations,
            r.coarse_davidson_iterations,
            r.coarse_setup.as_secs_f64(),
            r.coarse_session.as_secs_f64(),
            r.coarse_prolong.as_secs_f64(),
            r.coarse_phases.potential.as_secs_f64(),
            r.coarse_phases.diagonalization.as_secs_f64(),
            r.coarse_phases.density.as_secs_f64(),
            r.coarse_phases.energy.as_secs_f64(),
            r.coarse_phases.mixing.as_secs_f64(),
            r.coarse_clock.atomic_start.as_secs_f64()
        );
        out.outcome.result
    } else {
        problem.run().unwrap()
    };
    dft_engine::profile::finish_rank(parallel::SharedComm::active().map_or(0, |c| c.rank()));
    parallel::set_communicator(None);
    parallel::set_m2_override(None);
    result
}

/// Phase profile of the run just finished (all ranks of a simulated run, or the
/// serial run): the table and a JSON object (max over ranks of every entry,
/// `wall_min` the minimum). Empty without `BAANDI_PHASE_PROFILE=on`.
fn profile_block(wall: f64) -> String {
    let profiles = dft_engine::profile::take_registry();
    if profiles.is_empty() {
        return String::new();
    }
    let report = dft_engine::profile::Report::from_ranks(&profiles);
    print!("{}", report.table(Some(wall)));
    let buckets: Vec<String> = report
        .entries()
        .into_iter()
        .map(|(label, entries)| {
            let items: Vec<String> = entries
                .iter()
                .map(|(k, mx, _)| format!("\"{k}\":{mx:.6}"))
                .chain(
                    entries
                        .iter()
                        .take(1)
                        .map(|(_, _, mn)| format!("\"wall_min\":{mn:.6}")),
                )
                .collect();
            format!("\"{label}\":{{{}}}", items.join(","))
        })
        .collect();
    format!(
        ",\"profile\":{{\"ranks\":{},\"covered_wall_s\":{:.6},{}}}",
        report.ranks,
        report.covered_wall(),
        buckets.join(",")
    )
}

fn print_first_block() {
    let records = dft_engine::first_block::take_records();
    if !records.is_empty() {
        print!("{}", dft_engine::first_block::table(&records));
    }
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
    dft_engine::phase_clock::PhaseClock,
    Vec<(&'static str, u64)>,
);
fn snapshot() -> Snap {
    (
        parallel::stage_seconds(),
        parallel::m2_stats().1,
        parallel::collective_stats(),
        dft_engine::phase_clock::snapshot(),
        parallel::collective_bytes(),
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
    if std::env::var("M2_PROFILE_SERIAL").as_deref() == Ok("off") {
        for &p in &ranks {
            distributed_rows(&case, p, None);
        }
        return;
    }
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
    let _ = profile_block(serial_wall);
    print_first_block();
    for &p in &ranks {
        distributed_rows(&case, p, Some(&serial));
    }
}

fn distributed_rows(case: &str, p: usize, serial: Option<&ScfResult>) {
    {
        let solves_before = dft_solvers::distributed_davidson::reduced_solve_stats();
        let rot_before = dft_solvers::distributed_davidson::rotation_stats();
        let red_before = dft_solvers::distributed_davidson::reduction_stats();
        let lean_before = dft_solvers::distributed_davidson::lean_algebra_stats();
        let grid_before = (
            dft_solvers::mixing::mixer_counters(),
            dft_engine::stop_metric_counters(),
            dft_engine::grid_dist::counters(),
        );
        let before = snapshot();
        let t = Instant::now();
        let results = simulate(p, |comm: &SimRank| run(case, Some(Arc::new(comm.clone()))));
        let wall = t.elapsed().as_secs_f64();
        let after = snapshot();
        let r = &results[0];
        println!(
            "{case} M2 np={p}: wall {wall:.3} s, E {:.12} Ha, dE vs serial {:.3e} Ha, (SCF, Davidson, H) {:?}",
            r.total_hartree(),
            serial.map_or(f64::NAN, |s| (r.total_hartree() - s.total_hartree()).abs()),
            work(r)
        );
        let red = dft_solvers::distributed_davidson::reduction_stats();
        let solves_d = (red.0 - red_before.0).max(1) as f64;
        for (now, then) in red.1.iter().zip(&red_before.1) {
            let calls = (now.calls - then.calls) as f64 / p as f64;
            if calls > 0.0 {
                let h: Vec<u64> = (0..6)
                    .map(|b| (now.histogram[b] - then.histogram[b]) / p as u64)
                    .collect();
                println!(
                    "  reductions per rank {}: calls {calls:.0} ({:.1} per solve), payload {:.0} doubles, histogram {:?} (buckets {:?})",
                    now.site,
                    calls / (solves_d / p as f64),
                    (now.payload - then.payload) as f64 / p as f64,
                    h,
                    dft_solvers::distributed_davidson::REDUCTION_BUCKETS
                );
            }
        }
        let rot = dft_solvers::distributed_davidson::rotation_stats();
        println!(
            "  ritz rotations per rank: pairs {} columns {} flops {:.4e}",
            (rot.rotations - rot_before.rotations) / p as u64,
            (rot.columns - rot_before.columns) / p as u64,
            (rot.flops - rot_before.flops) as f64 / p as f64
        );
        let s = dft_solvers::distributed_davidson::reduced_solve_stats();
        let pp = p as u64;
        let solves = (s.solves - solves_before.solves) / pp;
        println!(
            "  reduced solves per rank: {solves}, orders mean {:.1} (cumulative over the process) p90 {} max {}, sum m^3 {:.3e}, single-column {} bordered (summed over ranks) {} broadcast {}",
            s.mean_order,
            s.p90_order,
            s.max_order,
            (s.sum_cubed - solves_before.sum_cubed) / p as f64,
            (s.single_column - solves_before.single_column) / pp,
            (s.bordered - solves_before.bordered),
            (s.broadcast - solves_before.broadcast) / pp
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
        let clock = dft_engine::phase_clock::snapshot().since(&before.3);
        // The ladder clock printed above is process-wide (summed over simulated ranks).
        println!(
            "  clocks per rank: certificate {:.3} s in {} blocks, screen {:.3} s, atomic start {:.3} s",
            clock.certificate.as_secs_f64() / p as f64,
            clock.certificate_blocks / p as u64,
            clock.screen.as_secs_f64() / p as f64,
            clock.atomic_start.as_secs_f64() / p as f64
        );
        // SCF phase seconds summed over the recorded iterations (MPI-SCALING-GAP):
        // maximum over the ranks of each phase; the remainder of the wall is work
        // outside the recorded maps (setup, verification maps, ladder coarse stage).
        let phase = |f: fn(&dft_engine::PhaseTimings) -> std::time::Duration| -> f64 {
            results
                .iter()
                .map(|r| {
                    r.history
                        .iter()
                        .map(|h| f(&h.timings).as_secs_f64())
                        .sum::<f64>()
                })
                .fold(0.0, f64::max)
        };
        let phases = [
            ("potential", phase(|t| t.potential)),
            ("diagonalization", phase(|t| t.diagonalization)),
            ("density", phase(|t| t.density)),
            ("energy", phase(|t| t.energy)),
            ("mixing", phase(|t| t.mixing)),
        ];
        println!(
            "  phases (recorded iterations, max over ranks): {}; outside {:.3} s",
            phases
                .iter()
                .map(|(n, s)| format!("{n} {s:.3}"))
                .collect::<Vec<_>>()
                .join(", "),
            wall - phases.iter().map(|(_, s)| s).sum::<f64>()
        );
        let density: u64 = r.density.iter().fold(0xcbf2_9ce4_8422_2325, |h: u64, x| {
            (h ^ x.to_bits()).wrapping_mul(0x0100_0000_01b3)
        });
        println!(
            "  energy bits {:#018x}, density FNV-1a {density:#018x}",
            r.total_hartree().to_bits()
        );
        // Machine-readable line for scripts/fast_screen_mpi.py (per-rank values).
        let (scf, dav, h) = work(r);
        let mut j = format!(
            "{{\"case\":\"{case}\",\"ranks\":{p},\"wall_s\":{wall:.6},\"energy_ha\":{:.17e},\"energy_bits\":\"{:#018x}\",\"density_fnv\":\"{density:#018x}\",\"scf\":{scf},\"davidson\":{dav},\"h_columns\":{h}",
            r.total_hartree(),
            r.total_hartree().to_bits()
        );
        j.push_str(&format!(
            ",\"rotation_flops\":{:.6e},\"reduced_solves\":{solves},\"reduced_sum_m3\":{:.6e}",
            (rot.flops - rot_before.flops) as f64 / p as f64,
            (s.sum_cubed - solves_before.sum_cubed) / p as f64
        ));
        let sites: Vec<String> = red
            .1
            .iter()
            .zip(&red_before.1)
            .map(|(now, then)| {
                format!(
                    "\"{}\":[{},{}]",
                    now.site,
                    (now.calls - then.calls) / p as u64,
                    (now.payload - then.payload) / p as u64
                )
            })
            .collect();
        j.push_str(&format!(",\"reductions\":{{{}}}", sites.join(",")));
        let bytes = parallel::collective_bytes();
        let coll: Vec<String> = after
            .2
            .iter()
            .zip(&before.2)
            .zip(bytes.iter().zip(&before.4))
            .map(|(((n, c, a), (_, d, b)), ((_, x), (_, y)))| {
                format!(
                    "\"{n}\":[{},{},{:.6}]",
                    (c - d) / p,
                    (x - y) / p as u64,
                    per(*a, *b)
                )
            })
            .collect();
        j.push_str(&format!(",\"collectives\":{{{}}}", coll.join(",")));
        let stages: Vec<String> = after
            .0
            .iter()
            .zip(&before.0)
            .map(|((n, a), (_, b))| format!("\"{n}\":{:.6}", per(*a, *b)))
            .chain(
                after
                    .1
                    .iter()
                    .zip(&before.1)
                    .map(|((n, a), (_, b))| format!("\"m2.{n}\":{:.6}", per(*a, *b))),
            )
            .collect();
        j.push_str(&format!(",\"stages\":{{{}}}", stages.join(",")));
        let ph: Vec<String> = phases
            .iter()
            .map(|(n, s)| format!("\"{n}\":{s:.6}"))
            .collect();
        j.push_str(&format!(",\"phases\":{{{}}}", ph.join(",")));
        let mut counters: Vec<String> = parallel::m2_stats()
            .0
            .iter()
            .map(|(n, v)| format!("\"{n}\":{}", v / p))
            .collect();
        // MPI-OVERLAP-DAVIDSON lean algebra counters (per rank).
        let la = dft_solvers::distributed_davidson::lean_algebra_stats();
        let lb = &lean_before;
        for (n, v) in [
            ("lean_solves", la.solves - lb.solves),
            (
                "lean_border_reductions",
                la.border_reductions - lb.border_reductions,
            ),
            ("lean_reprojections", la.reprojections - lb.reprojections),
            ("lean_fallbacks", la.fallbacks - lb.fallbacks),
            ("lean_rebuilds", la.rebuilds - lb.rebuilds),
            ("lean_wasted_columns", la.wasted_columns - lb.wasted_columns),
            (
                "lean_factored_starts",
                la.factored_starts - lb.factored_starts,
            ),
            // BREAKTHROUGH-1-RESTART R1-b (`BAANDI_RESTART_KEEP`).
            ("lean_contractions", la.contractions - lb.contractions),
            (
                "keep_contractions",
                la.keep_contractions - lb.keep_contractions,
            ),
            ("kept_directions", la.kept_directions - lb.kept_directions),
            ("keep_dropped", la.keep_dropped - lb.keep_dropped),
            ("keep_fallbacks", la.keep_fallbacks - lb.keep_fallbacks),
        ] {
            counters.push(format!("\"{n}\":{}", v / p as u64));
        }
        j.push_str(&format!(",\"m2_counters\":{{{}}}", counters.join(",")));
        // MPI-REMAINING: replicated grid-stage counters (process-wide over the
        // simulated ranks, divided by the rank count like every other entry).
        let mix = dft_solvers::mixing::mixer_counters();
        let stop = dft_engine::stop_metric_counters();
        // M2-GRID: distributed potential maps, audits, stop estimates, plans, declined.
        let dist = dft_engine::grid_dist::counters();
        j.push_str(&format!(
            ",\"grid_counters\":{{\"mix_calls\":{},\"mix_diffs\":{},\"mix_dots\":{},\"mix_metric\":{},\"stop_kept\":{},\"stop_fft\":{},\"dist_potential\":{},\"dist_audit\":{},\"dist_stop\":{},\"dist_plans\":{},\"dist_declined\":{}}}",
            (mix[0] - grid_before.0[0]) / p as u64,
            (mix[1] - grid_before.0[1]) / p as u64,
            (mix[2] - grid_before.0[2]) / p as u64,
            (mix[3] - grid_before.0[3]) / p as u64,
            (stop[0] - grid_before.1[0]) / p as u64,
            (stop[1] - grid_before.1[1]) / p as u64,
            (dist[0] - grid_before.2[0]) / p as u64,
            (dist[1] - grid_before.2[1]) / p as u64,
            (dist[2] - grid_before.2[2]) / p as u64,
            (dist[3] - grid_before.2[3]) / p as u64,
            (dist[4] - grid_before.2[4]) / p as u64,
        ));
        // BREAKTHROUGH-0-INSTRUMENT: exclusive phase profile (per rank, max
        // over ranks) and the first-block residual table (both default off).
        j.push_str(&profile_block(wall));
        print_first_block();
        j.push('}');
        println!("M2SCREEN_JSON {j}");
    }
}
