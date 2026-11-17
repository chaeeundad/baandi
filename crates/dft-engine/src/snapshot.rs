//! SCF snapshot capture and replay (ticket BREAKTHROUGH-0-INSTRUMENT, review
//! 2026-11-17 section 13.3).
//!
//! `BAANDI_SNAPSHOT_DIR=<absolute dir>` (default unset): every real Gamma
//! eigensolve of the run (cold / atomic, ladder coarse maps, the coarse
//! hand-off solve, fine early / middle / late, the first-step redo at the same
//! Hamiltonian, every verification map) writes, at the solver entry, the state
//! that replays exactly that solve:
//!
//! * `<label>.json` (rank 0): the manifest: context (stage, kind, SCF iteration,
//!   sequence number), problem identity (electrons, bands, occupied, plane waves,
//!   grid, cutoff), the operator version (FNV-1a of the local potential),
//!   tolerance and band policy, the start kind (none / warm / kept and whether
//!   the kept start was armed), rank count and layout, and after the solve its
//!   outcome (Davidson iterations, H columns, restarts, contractions, whether
//!   the kept start was taken, eigenvalue FNV-1a) so that a replay can check
//!   that it reproduced the work;
//! * `<label>.density.bin` (rank 0): the map's input density, headerless
//!   little-endian f64, C order (the `density.bin` layout of the CLI);
//! * `<label>.rank<r>.bin` (every rank): this rank's rows of the start X (the
//!   stick rows on the M2 path, the full ORC block on the threads path), the
//!   kept H X when one is armed (already brought to the current potential), the
//!   FNV-1a of the row list.
//!
//! Previous Ritz directions are not kept by the current solvers (review R1):
//! the manifest records `"previous_directions": 0`. The delta-V update of the
//! kept H X is not replayed (its result is stored), so a replay reproduces the
//! solver's own counters (iterations, H columns, restarts, eigenvalues), not
//! the `potential_change_columns` of the operator.
//!
//! File writing happens outside the stage clocks and is subtracted from the
//! map's `diagonalization` phase (`take_io_seconds`); the run's total wall is
//! not a timing measurement while the switch is set.
//!
//! Replay: `davidson_bench --snapshot <label>.json` loads the files, builds
//! the problem from the same input (the manifest's cutoff for coarse maps),
//! injects the start through [`inject`] and runs the session map
//! (`ScfProblem::replay_solve`).
use crate::profile;
use dft_numerics::real::RealMatrix;
use dft_physics::data::Array3;
use dft_solvers::real_davidson::RealEigenSolution;
use std::path::{Path, PathBuf};
use std::sync::{Mutex, OnceLock};
use std::time::Instant;

static FORCED_DIR: Mutex<Option<Option<PathBuf>>> = Mutex::new(None);
/// Tests: force the capture directory (`Some(None)` = capture off,
/// `None` = the switch).
pub fn force_dir(dir: Option<Option<PathBuf>>) {
    *FORCED_DIR.lock().expect("snapshot dir") = dir;
    // A new capture series: every rank's sequence numbers restart at 0 (the
    // labels must agree across ranks).
    SLOTS.lock().expect("snapshot slots").clear();
}
/// The capture directory, when the switch is set to an absolute path.
pub fn dir() -> Option<PathBuf> {
    if let Some(forced) = FORCED_DIR.lock().expect("snapshot dir").as_ref() {
        return forced.clone();
    }
    static DIR: OnceLock<Option<PathBuf>> = OnceLock::new();
    DIR.get_or_init(|| {
        let t = dft_numerics::switches::text(dft_numerics::switches::Id::SNAPSHOT_DIR)?;
        let p = PathBuf::from(t);
        if p.is_absolute() {
            Some(p)
        } else {
            eprintln!("BAANDI_SNAPSHOT_DIR must be an absolute directory; capture off");
            None
        }
    })
    .clone()
}

#[derive(Default)]
struct RankSlot {
    density: Option<([usize; 3], Vec<f64>)>,
    seq: usize,
    io_seconds: f64,
    last_manifest: Option<PathBuf>,
}
static SLOTS: Mutex<Vec<RankSlot>> = Mutex::new(Vec::new());
fn with_slot<R>(f: impl FnOnce(&mut RankSlot) -> R) -> R {
    let r = profile::rank_key();
    let mut v = SLOTS.lock().expect("snapshot slots");
    if v.len() <= r {
        v.resize_with(r + 1, RankSlot::default);
    }
    f(&mut v[r])
}

/// The input density of the map about to solve (copied only when capturing).
pub(crate) fn set_density(n: &Array3<f64>) {
    if dir().is_none() {
        return;
    }
    let dims = [n.shape()[0], n.shape()[1], n.shape()[2]];
    let v = n.iter().copied().collect();
    with_slot(|s| s.density = Some((dims, v)));
}
/// Seconds spent writing since the last call (the caller subtracts them from
/// its timed phase).
pub(crate) fn take_io_seconds() -> f64 {
    if dir().is_none() {
        return 0.0;
    }
    with_slot(|s| std::mem::take(&mut s.io_seconds))
}

/// Start layout of a rank file.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Layout {
    /// This rank's stick rows (M2 path), `rows_fnv` of the row list.
    StickRows,
    /// The full ORC block (threads path).
    OrcFull,
}
impl Layout {
    fn name(self) -> &'static str {
        match self {
            Layout::StickRows => "stick_rows",
            Layout::OrcFull => "orc_full",
        }
    }
}

/// What the solver entry hands to the capture.
pub(crate) struct Entry<'a> {
    pub rank: usize,
    pub ranks: usize,
    pub layout: Layout,
    pub x: Option<&'a RealMatrix>,
    pub hx: Option<&'a RealMatrix>,
    pub armed: bool,
    pub rows_fnv: u64,
    pub n_global: usize,
    pub bands: usize,
    pub occupied: usize,
    pub tolerance: f64,
    pub strict_bands: usize,
    pub buffer_tolerance: f64,
    pub ritz: bool,
    pub potential: &'a [f64],
    pub grid: [usize; 3],
    pub cutoff_hartree: f64,
    pub electrons: f64,
    pub capacity_per_band: usize,
}

/// A capture in progress (manifest completed by `capture_exit`).
pub(crate) struct Capture {
    label: String,
    manifest: Option<PathBuf>,
    head: String,
}

pub fn fnv(values: &[f64]) -> u64 {
    values.iter().fold(0xcbf2_9ce4_8422_2325, |h: u64, x| {
        (h ^ x.to_bits()).wrapping_mul(0x0100_0000_01b3)
    })
}
pub fn fnv_usize(values: &[usize]) -> u64 {
    values.iter().fold(0xcbf2_9ce4_8422_2325, |h: u64, x| {
        (h ^ *x as u64).wrapping_mul(0x0100_0000_01b3)
    })
}

const MAGIC: &[u8; 8] = b"BAANDSN1";

/// Write a rank file: `x` (rows x cols, column-major), optional `hx`.
pub fn write_rank_file(
    path: &Path,
    x: &RealMatrix,
    hx: Option<&RealMatrix>,
    rows_fnv: u64,
    layout: &str,
) -> std::io::Result<()> {
    let mut bytes = Vec::with_capacity(64 + 8 * x.as_slice().len() * 2);
    bytes.extend_from_slice(MAGIC);
    let layout_code: u64 = if layout == "stick_rows" { 1 } else { 2 };
    for v in [
        x.rows() as u64,
        x.cols() as u64,
        u64::from(hx.is_some()),
        rows_fnv,
        layout_code,
    ] {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    for v in x.as_slice() {
        bytes.extend_from_slice(&v.to_le_bytes());
    }
    if let Some(h) = hx {
        for v in h.as_slice() {
            bytes.extend_from_slice(&v.to_le_bytes());
        }
    }
    std::fs::write(path, bytes)
}
/// Read a rank file: (x, hx, rows_fnv, layout name).
pub fn read_rank_file(
    path: &Path,
) -> Result<(RealMatrix, Option<RealMatrix>, u64, String), String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    if bytes.len() < 48 || &bytes[..8] != MAGIC {
        return Err(format!("{}: not a snapshot rank file", path.display()));
    }
    let word = |i: usize| u64::from_le_bytes(bytes[8 + 8 * i..16 + 8 * i].try_into().unwrap());
    let (rows, cols, has_hx, rows_fnv, layout) = (
        word(0) as usize,
        word(1) as usize,
        word(2) != 0,
        word(3),
        word(4),
    );
    let n = rows * cols;
    let need = 48 + 8 * n * (1 + usize::from(has_hx));
    if bytes.len() != need {
        return Err(format!(
            "{}: {} bytes, expected {need}",
            path.display(),
            bytes.len()
        ));
    }
    let read = |from: usize| -> Vec<f64> {
        bytes[from..from + 8 * n]
            .as_chunks::<8>()
            .0
            .iter()
            .map(|c| f64::from_le_bytes(*c))
            .collect()
    };
    let x = RealMatrix::new(rows, cols, read(48)).map_err(|e| e.to_string())?;
    let hx = if has_hx {
        Some(RealMatrix::new(rows, cols, read(48 + 8 * n)).map_err(|e| e.to_string())?)
    } else {
        None
    };
    Ok((
        x,
        hx,
        rows_fnv,
        if layout == 1 {
            "stick_rows"
        } else {
            "orc_full"
        }
        .to_string(),
    ))
}
/// Read a density file (headerless little-endian f64, C order) of `dims`.
pub fn read_density_file(path: &Path, dims: [usize; 3]) -> Result<Array3<f64>, String> {
    let bytes = std::fs::read(path).map_err(|e| format!("{}: {e}", path.display()))?;
    let n: usize = dims.iter().product();
    if bytes.len() != 8 * n {
        return Err(format!(
            "{}: {} bytes, grid {dims:?} needs {}",
            path.display(),
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
    Array3::from_shape_vec((dims[0], dims[1], dims[2]), v).map_err(|e| e.to_string())
}

/// Capture the solver entry (every rank writes its rows; rank 0 the density
/// and the manifest head). `None` when capture is off or the start is absent
/// on this rank.
pub(crate) fn capture_entry(e: Entry<'_>) -> Option<Capture> {
    let dir = dir()?;
    let clock = Instant::now();
    let ctx = profile::context();
    let (seq, density) = with_slot(|s| {
        let seq = s.seq;
        s.seq += 1;
        (seq, s.density.take())
    });
    let label = format!(
        "s{seq:03}-{}-{}-it{:02}",
        ctx.stage, ctx.kind, ctx.iteration
    );
    let mut ok = true;
    if let Some(x) = e.x {
        let path = dir.join(format!("{label}.rank{}.bin", e.rank));
        if let Err(err) = write_rank_file(
            &path,
            x,
            e.hx.filter(|_| e.armed),
            e.rows_fnv,
            e.layout.name(),
        ) {
            eprintln!("snapshot: {}: {err}", path.display());
            ok = false;
        }
    }
    let mut manifest = None;
    let mut head = String::new();
    if e.rank == 0 {
        let mut density_fnv = 0u64;
        if let Some((dims, v)) = &density {
            density_fnv = fnv(v);
            let path = dir.join(format!("{label}.density.bin"));
            let bytes: Vec<u8> = v.iter().flat_map(|x| x.to_le_bytes()).collect();
            if let Err(err) = std::fs::write(&path, bytes) {
                eprintln!("snapshot: {}: {err}", path.display());
                ok = false;
            }
            if *dims != e.grid {
                eprintln!("snapshot: density grid {dims:?} differs from the operator grid");
            }
        }
        let start = match (e.x, e.armed) {
            (None, _) => "none",
            (Some(_), true) => "kept",
            (Some(_), false) => "warm",
        };
        head = format!(
            concat!(
                "{{\n  \"schema\": \"baandi-snapshot-1\",\n  \"label\": \"{label}\",\n",
                "  \"seq\": {seq},\n  \"stage\": \"{stage}\",\n  \"kind\": \"{kind}\",\n  \"iteration\": {it},\n",
                "  \"ranks\": {ranks},\n  \"layout\": \"{layout}\",\n  \"bands\": {bands},\n  \"occupied\": {occ},\n",
                "  \"plane_waves\": {n},\n  \"grid\": [{g0}, {g1}, {g2}],\n  \"cutoff_hartree\": {cut:.17e},\n",
                "  \"electrons\": {el:.17e},\n  \"capacity_per_band\": {cap},\n  \"tolerance_ha\": {tol:.17e},\n",
                "  \"strict_bands\": {sb},\n  \"buffer_tolerance_ha\": {bt:.17e},\n  \"ritz\": {ritz},\n",
                "  \"start\": \"{start}\",\n  \"kept_armed\": {armed},\n  \"has_hx\": {hashx},\n",
                "  \"previous_directions\": 0,\n  \"potential_fnv\": \"{pfnv:#018x}\",\n",
                "  \"density_fnv\": \"{dfnv:#018x}\",\n  \"density_file\": {dfile},\n  \"rank_files\": \"{label}.rank<r>.bin\",\n",
                "  \"switches\": {{\"M2_REUSE\": \"{sw_reuse}\", \"M2_LEAN_ALGEBRA\": \"{sw_lean}\", \"M2_LAZY\": \"{sw_lazy}\", \"REUSE_PSI\": \"{sw_psi}\", \"DAVIDSON\": \"{sw_dav}\"}},\n",
                "  \"verification_passed\": null"
            ),
            label = label,
            seq = seq,
            stage = ctx.stage,
            kind = ctx.kind,
            it = ctx.iteration,
            ranks = e.ranks,
            layout = e.layout.name(),
            bands = e.bands,
            occ = e.occupied,
            n = e.n_global,
            g0 = e.grid[0],
            g1 = e.grid[1],
            g2 = e.grid[2],
            cut = e.cutoff_hartree,
            el = e.electrons,
            cap = e.capacity_per_band,
            tol = e.tolerance,
            sb = e.strict_bands,
            bt = e.buffer_tolerance,
            ritz = e.ritz,
            start = start,
            armed = e.armed,
            hashx = e.hx.is_some() && e.armed,
            pfnv = fnv(e.potential),
            dfnv = density_fnv,
            dfile = if density.is_some() {
                format!("\"{label}.density.bin\"")
            } else {
                "null".to_string()
            },
            sw_reuse = dft_numerics::switches::choice(dft_numerics::switches::Id::M2_REUSE),
            sw_lean = dft_numerics::switches::choice(dft_numerics::switches::Id::M2_LEAN_ALGEBRA),
            sw_lazy = dft_numerics::switches::choice(dft_numerics::switches::Id::M2_LAZY),
            sw_psi = dft_numerics::switches::choice(dft_numerics::switches::Id::REUSE_PSI),
            sw_dav = dft_numerics::switches::choice(dft_numerics::switches::Id::DAVIDSON),
        );
        manifest = Some(dir.join(format!("{label}.json")));
    }
    with_slot(|s| s.io_seconds += clock.elapsed().as_secs_f64());
    ok.then_some(Capture {
        label,
        manifest,
        head,
    })
}

/// Complete the manifest with the solve's outcome (rank 0).
pub(crate) fn capture_exit(c: Capture, result: &RealEigenSolution, kept_taken: bool) {
    let clock = Instant::now();
    if let Some(path) = &c.manifest {
        let text = format!(
            "{},\n  \"result\": {{\"davidson_iterations\": {}, \"h_columns\": {}, \"restarts\": {}, \"contractions\": {}, \"kept_taken\": {}, \"converged\": {}, \"eigen_fnv\": \"{:#018x}\", \"residual_max\": {:.17e}}}\n}}\n",
            c.head,
            result.iterations,
            result.operator_applications,
            result.restarts,
            result.contracted_restarts,
            kept_taken,
            result.converged.iter().filter(|&&b| b).count(),
            fnv(&result.eigenvalues),
            result.residuals.iter().copied().fold(0.0, f64::max)
        );
        if let Err(err) = std::fs::write(path, text) {
            eprintln!("snapshot: {}: {err}", path.display());
        }
        with_slot(|s| s.last_manifest = Some(path.clone()));
    }
    let _ = c.label;
    with_slot(|s| s.io_seconds += clock.elapsed().as_secs_f64());
}

/// Record the outcome of the verification whose map was captured last
/// (rank 0 rewrites the manifest; outside any timed phase).
pub(crate) fn annotate_verification(passed: bool) {
    if dir().is_none() {
        return;
    }
    let path = with_slot(|s| s.last_manifest.take());
    if let Some(path) = path
        && let Ok(text) = std::fs::read_to_string(&path)
    {
        let text = text.replacen(
            "\"verification_passed\": null",
            &format!("\"verification_passed\": {passed}"),
            1,
        );
        let _ = std::fs::write(&path, text);
    }
}

/// A start injected into the next solve of a rank (replay).
pub struct ReplayStart {
    /// This rank's rows of X (M2 path) or the full ORC block (threads path).
    pub x: RealMatrix,
    /// H X in the same layout, when the kept start was armed.
    pub hx: Option<RealMatrix>,
}
static INJECT: Mutex<Vec<Option<ReplayStart>>> = Mutex::new(Vec::new());
static INJECT_PENDING: std::sync::atomic::AtomicBool = std::sync::atomic::AtomicBool::new(false);
/// Inject a start for the next real Gamma solve on `rank`.
pub fn inject(rank: usize, start: ReplayStart) {
    let mut v = INJECT.lock().expect("replay injection");
    if v.len() <= rank {
        v.resize_with(rank + 1, || None);
    }
    v[rank] = Some(start);
    INJECT_PENDING.store(true, std::sync::atomic::Ordering::Relaxed);
}
pub(crate) fn take_injection() -> Option<ReplayStart> {
    if !INJECT_PENDING.load(std::sync::atomic::Ordering::Relaxed) {
        return None;
    }
    let r = profile::rank_key();
    let mut v = INJECT.lock().expect("replay injection");
    let out = v.get_mut(r).and_then(Option::take);
    if v.iter().all(Option::is_none) {
        INJECT_PENDING.store(false, std::sync::atomic::Ordering::Relaxed);
    }
    out
}
