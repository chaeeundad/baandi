//! ISSUE-18 (ticket ISSUE18-SNAPSHOT-IDENTITY): a snapshot replay validates the
//! operator and row-layout identity BEFORE the captured X and H X are used.
//! Exact replay rejects (a) a moved atom, (b) a different same-valence
//! pseudopotential, (c) a changed XC functional, (d) a changed local row
//! order/ownership (2 simulated ranks; a rank-local mismatch becomes an agreed
//! error on every rank), (e) modified density bytes and (f) an incompatible
//! stopping profile, while ranks, grid and bands are unchanged; the original
//! snapshot reproduces its captured contract. An explicitly labeled A/B replay
//! with a changed potential rounding recomputes H X instead of using the
//! stored one. Manifests without a fingerprint load only in legacy mode.
//! Own test binary with a lock (process-wide forced switches).
#[allow(dead_code)]
mod common;
use dft_engine::snapshot::{self, ReplayContract, ReplayMode};
use dft_engine::{Functional, GammaReal, ScfOptions, ScfProblem, parallel};
use dft_geometry::{Atom, Cell};
use dft_mpi::{Collective, SimRank, simulate};
use std::path::{Path, PathBuf};
use std::sync::{Arc, Mutex};

static LOCK: Mutex<()> = Mutex::new(());

fn options() -> ScfOptions {
    ScfOptions {
        gamma_real: GammaReal::Auto,
        cutoff_hartree: 3.,
        bands: Some(6),
        threads: 1,
        ..ScfOptions::strict_reference()
    }
}

fn problem_with(
    pp: &dft_upf::Pseudopotential,
    second: [f64; 3],
    options: ScfOptions,
) -> ScfProblem {
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
    let atoms = vec![Atom::new("Xx", [0.13, 0.21, 0.32]), Atom::new("Xx", second)];
    ScfProblem::new(
        &cell,
        &atoms,
        &[dft_engine::SpeciesPotential {
            label: "Xx",
            potential: pp,
        }],
        options,
    )
    .unwrap()
}

const SECOND: [f64; 3] = [0.61, 0.57, 0.74];

fn problem(pp: &dft_upf::Pseudopotential) -> ScfProblem {
    problem_with(pp, SECOND, options())
}

fn scratch(name: &str) -> PathBuf {
    let dir = std::env::temp_dir().join(format!("baandi-identity-{}-{name}", std::process::id()));
    let _ = std::fs::remove_dir_all(&dir);
    std::fs::create_dir_all(&dir).unwrap();
    dir
}

fn with_ranks<R: Send>(ranks: usize, f: impl Fn(Option<&SimRank>) -> R + Sync) -> Vec<R> {
    if ranks == 1 {
        return vec![f(None)];
    }
    simulate(ranks, |comm: &SimRank| {
        parallel::set_m2_override(Some(true));
        dft_engine::gamma_real::set_m2_reuse_override(Some(true));
        let shared: Arc<dyn Collective + Send + Sync> = Arc::new(comm.clone());
        parallel::set_communicator(Some(shared));
        let r = f(Some(comm));
        parallel::set_communicator(None);
        dft_engine::gamma_real::set_m2_reuse_override(None);
        parallel::set_m2_override(None);
        r
    })
}

/// Run the SCF with capture on; return the sorted manifests.
fn capture(ranks: usize, pp: &dft_upf::Pseudopotential, dir: &Path) -> Vec<PathBuf> {
    snapshot::force_dir(Some(Some(dir.to_path_buf())));
    with_ranks(ranks, |_| problem(pp).run().unwrap());
    snapshot::force_dir(Some(None));
    let mut v: Vec<PathBuf> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|e| e == "json"))
        .collect();
    v.sort();
    v
}

/// A kept-start manifest (X and the stored H X) from the middle of the run.
fn kept_manifest(files: &[PathBuf]) -> PathBuf {
    let kept: Vec<&PathBuf> = files
        .iter()
        .filter(|f| {
            let t = std::fs::read_to_string(f).unwrap();
            t.contains("\"start\": \"kept\"") && t.contains("\"has_hx\": true")
        })
        .collect();
    assert!(!kept.is_empty(), "a kept-start snapshot with H X exists");
    kept[kept.len() / 2].clone()
}

struct Loaded {
    text: String,
    label: String,
    grid: [usize; 3],
    dir: PathBuf,
}
fn load(path: &Path) -> Loaded {
    let text = std::fs::read_to_string(path).unwrap();
    let c = ReplayContract::from_manifest(&text, ReplayMode::Exact).unwrap();
    Loaded {
        label: c.label.clone(),
        grid: c.grid,
        dir: path.parent().unwrap().to_path_buf(),
        text,
    }
}

type Outcome = Result<(usize, usize, String, Vec<f64>, snapshot::ReplayReport), String>;

/// Replay `m` on `p` (built per rank) with the given mode; `tamper` edits the
/// rank's contract and `density_edit` the loaded density.
#[allow(clippy::too_many_arguments)]
fn replay(
    m: &Loaded,
    ranks: usize,
    mode: ReplayMode,
    text: Option<&str>,
    build: &(dyn Fn() -> ScfProblem + Sync),
    rank_file: &(dyn Fn(usize) -> usize + Sync),
    tamper: &(dyn Fn(usize, &mut ReplayContract) + Sync),
    density_edit: &(dyn Fn(&mut [f64]) + Sync),
) -> Vec<Outcome> {
    let text = text.unwrap_or(&m.text);
    with_ranks(ranks, |comm| {
        let rank = comm.map_or(0, |c| c.rank());
        let mut density =
            snapshot::read_density_file(&m.dir.join(format!("{}.density.bin", m.label)), m.grid)
                .unwrap();
        density_edit(density.as_slice_mut().unwrap());
        let file = m
            .dir
            .join(format!("{}.rank{}.bin", m.label, rank_file(rank)));
        let (x, hx, rows_fnv, layout) = snapshot::read_rank_file(&file).unwrap();
        let mut contract = ReplayContract::from_manifest(text, mode)
            .map_err(|e| e.to_string())?
            .with_rank_file(rows_fnv, &layout);
        tamper(rank, &mut contract);
        let mut p = build();
        let (step, report) = p
            .replay_checked(&density, snapshot::ReplayStart { x, hx }, &contract)
            .map_err(|e| e.to_string())?;
        let s = &step.k_points[0].solution;
        Ok((
            s.iterations,
            s.operator_applications,
            format!("{:#018x}", snapshot::fnv(&s.eigenvalues)),
            s.eigenvalues.clone(),
            report,
        ))
    })
}

fn field(text: &str, key: &str) -> String {
    let at = text.find(&format!("\"{key}\":")).unwrap();
    let rest = &text[at + key.len() + 3..];
    let end = rest.find([',', '}', '\n']).unwrap();
    rest[..end].trim().trim_matches('"').to_string()
}

fn same(_: usize) -> usize {
    0
}
fn own(r: usize) -> usize {
    r
}
fn keep(_: usize, _: &mut ReplayContract) {}
fn no_edit(_: &mut [f64]) {}

fn expect_reject(outs: &[Outcome], what: &str, needle: &str) {
    for (r, o) in outs.iter().enumerate() {
        match o {
            Ok(_) => panic!("{what}: rank {r} accepted a stale snapshot"),
            Err(e) => {
                println!("{what}: rank {r} rejected: {e}");
                assert!(e.contains(needle), "{what}: rank {r}: {e} lacks {needle:?}");
            }
        }
    }
}

#[test]
fn exact_replay_rejects_every_changed_operator_and_reproduces_the_original() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    let dir = scratch("p1");
    let files = capture(1, &pp, &dir);
    let m = load(&kept_manifest(&files));
    assert!(
        m.text.contains("\"operator_fingerprint\": {"),
        "the manifest records the operator fingerprint"
    );
    // The original snapshot reproduces its captured contract.
    let build = || problem(&pp);
    let outs = replay(
        &m,
        1,
        ReplayMode::Exact,
        None,
        &build,
        &same,
        &keep,
        &no_edit,
    );
    let (iterations, h, fnv, eigen, report) = outs[0].as_ref().expect("original replays");
    assert_eq!(
        *iterations,
        field(&m.text, "davidson_iterations")
            .parse::<usize>()
            .unwrap()
    );
    assert_eq!(*h, field(&m.text, "h_columns").parse::<usize>().unwrap());
    assert_eq!(*fnv, field(&m.text, "eigen_fnv"));
    assert!(report.hx_injected && !report.hx_recomputed, "{report:?}");
    // (a) a moved atom.
    let moved = || problem_with(&pp, [0.62, 0.57, 0.74], options());
    expect_reject(
        &replay(
            &m,
            1,
            ReplayMode::Exact,
            None,
            &moved,
            &same,
            &keep,
            &no_edit,
        ),
        "moved atom",
        "positions",
    );
    // (b) a different pseudopotential with the same valence.
    let pp_b = common::potential(4., 1.1, true, 401);
    let other_pp = || problem(&pp_b);
    expect_reject(
        &replay(
            &m,
            1,
            ReplayMode::Exact,
            None,
            &other_pp,
            &same,
            &keep,
            &no_edit,
        ),
        "same-valence pseudopotential",
        "pseudo",
    );
    // (c) a changed XC functional.
    let xc = || {
        problem_with(
            &pp,
            SECOND,
            ScfOptions {
                functional: Functional::LdaPw,
                allow_functional_override: true,
                ..options()
            },
        )
    };
    expect_reject(
        &replay(&m, 1, ReplayMode::Exact, None, &xc, &same, &keep, &no_edit),
        "changed XC",
        "xc",
    );
    // (e) modified density bytes.
    let flip = |d: &mut [f64]| d[7] = f64::from_bits(d[7].to_bits() ^ 1);
    expect_reject(
        &replay(&m, 1, ReplayMode::Exact, None, &build, &same, &keep, &flip),
        "modified density",
        "density",
    );
    // (f) an incompatible stopping profile (strict reference arithmetic).
    let strict = || {
        problem_with(
            &pp,
            SECOND,
            ScfOptions {
                gamma_real: GammaReal::Auto,
                cutoff_hartree: 3.,
                bands: Some(6),
                threads: 1,
                ..ScfOptions::strict_reference()
            },
        )
    };
    expect_reject(
        &replay(
            &m,
            1,
            ReplayMode::Exact,
            None,
            &strict,
            &same,
            &keep,
            &no_edit,
        ),
        "stopping profile",
        "profile",
    );
    // A tampered captured tolerance is also an incompatible profile.
    let tol = |_: usize, c: &mut ReplayContract| c.strict_bands += 1;
    expect_reject(
        &replay(
            &m,
            1,
            ReplayMode::Exact,
            None,
            &build,
            &same,
            &tol,
            &no_edit,
        ),
        "strict band count",
        "profile",
    );
    // Legacy manifests (no fingerprint): exact replay refuses, the labeled
    // legacy mode accepts with a warning.
    let start = m.text.find(",\n  \"operator_fingerprint\"").unwrap();
    let end = start + 2 + m.text[start + 2..].find('\n').unwrap();
    let legacy = format!("{}{}", &m.text[..start], &m.text[end..]);
    assert!(ReplayContract::from_manifest(&legacy, ReplayMode::Exact).is_err());
    let outs = replay(
        &m,
        1,
        ReplayMode::Legacy,
        Some(&legacy),
        &build,
        &same,
        &keep,
        &no_edit,
    );
    let (_, _, legacy_fnv, _, report) = outs[0].as_ref().expect("legacy replays");
    assert_eq!(*legacy_fnv, *fnv);
    assert_eq!(report.mode, ReplayMode::Legacy);
    let _ = eigen;
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn ab_replay_with_changed_potential_rounding_recomputes_h_x() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    let dir = scratch("ab");
    let files = capture(1, &pp, &dir);
    let m = load(&kept_manifest(&files));
    let build = || problem(&pp);
    let exact = replay(
        &m,
        1,
        ReplayMode::Exact,
        None,
        &build,
        &same,
        &keep,
        &no_edit,
    );
    let (_, _, _, eigen, _) = exact[0].as_ref().unwrap();
    // A permitted kernel change that alters the potential's rounding: the
    // captured potential fingerprint no longer matches the current operator.
    let rounding = |_: usize, c: &mut ReplayContract| {
        c.potential_fnv = c.potential_fnv.map(|f| f ^ 1);
    };
    expect_reject(
        &replay(
            &m,
            1,
            ReplayMode::Exact,
            None,
            &build,
            &same,
            &rounding,
            &no_edit,
        ),
        "exact replay, changed potential rounding",
        "potential",
    );
    let outs = replay(
        &m,
        1,
        ReplayMode::AbSamePhysics,
        None,
        &build,
        &same,
        &rounding,
        &no_edit,
    );
    let (_, _, _, ab_eigen, report) = outs[0].as_ref().expect("A/B replays");
    assert_eq!(report.mode, ReplayMode::AbSamePhysics);
    assert!(
        !report.hx_injected && report.hx_recomputed,
        "stale H X dropped and recomputed: {report:?}"
    );
    let tolerance: f64 = field(&m.text, "tolerance_ha").parse().unwrap();
    let occupied = 4;
    for (a, b) in eigen.iter().zip(ab_eigen).take(occupied) {
        assert!(
            (a - b).abs() <= tolerance.max(1e-10),
            "A/B occupied eigenvalue {b} vs exact {a} (tolerance {tolerance:e} Ha)"
        );
    }
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn rank_local_row_mismatch_is_an_agreed_error() {
    let _lock = LOCK.lock().unwrap_or_else(|e| e.into_inner());
    let pp = common::potential(4., 1., true, 401);
    let dir = scratch("p2");
    let files = capture(2, &pp, &dir);
    let m = load(&kept_manifest(&files));
    let build = || problem(&pp);
    let outs = replay(
        &m,
        2,
        ReplayMode::Exact,
        None,
        &build,
        &own,
        &keep,
        &no_edit,
    );
    for o in &outs {
        let (iterations, h, fnv, _, report) = o.as_ref().expect("original replays on 2 ranks");
        assert_eq!(
            *iterations,
            field(&m.text, "davidson_iterations")
                .parse::<usize>()
                .unwrap()
        );
        assert_eq!(*h, field(&m.text, "h_columns").parse::<usize>().unwrap());
        assert_eq!(*fnv, field(&m.text, "eigen_fnv"));
        assert!(report.hx_injected);
    }
    // (d) changed row ownership: the ranks read each other's rows.
    let swapped = |r: usize| 1 - r;
    expect_reject(
        &replay(
            &m,
            2,
            ReplayMode::Exact,
            None,
            &build,
            &swapped,
            &keep,
            &no_edit,
        ),
        "swapped row ownership",
        "rows",
    );
    // (d) a changed local row order on rank 1 only: rank 0 agrees on the error
    // before any rank enters the solver's collectives (no hang, no solve).
    let reorder = |r: usize, c: &mut ReplayContract| {
        if r == 1 {
            c.rows_fnv = c.rows_fnv.map(|f| f.wrapping_add(1));
        }
    };
    let outs = replay(
        &m,
        2,
        ReplayMode::Exact,
        None,
        &build,
        &own,
        &reorder,
        &no_edit,
    );
    expect_reject(&outs, "rank-local row order", "rows");
    // A/B mode never accepts a different row interpretation either.
    expect_reject(
        &replay(
            &m,
            2,
            ReplayMode::AbSamePhysics,
            None,
            &build,
            &own,
            &reorder,
            &no_edit,
        ),
        "A/B rank-local row order",
        "rows",
    );
    let _ = std::fs::remove_dir_all(&dir);
}
