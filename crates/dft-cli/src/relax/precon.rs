//! Preconditioned L-BFGS for fixed-cell ionic relaxation (PRELAX-P2, ADR 0058).
//!
//! All coordinates, steps and gradients are flat Cartesian Bohr vectors of
//! length 3N (atom-major); fixed Cartesian components are zero in every vector
//! and eliminated from the preconditioner (rows and columns removed).
//!
//! Preconditioner (model Hessian) P = mu P1 with
//!   P1 = sum over neighbour pairs (i, j, image) of w_ij B_ij (e_i - e_j)(e_i - e_j)^T
//!        + c_stab I,
//!   B_ij = (1 - sigma) I + sigma n n^T  (n the unit bond vector, stretch-dominant),
//!   w_ij = exp(-A (d_ij / r0_ij - 1)),  r0_ij = r_cov,i + r_cov,j.
//! P1 is dimensionless; mu (Ha/Bohr^2) is estimated from the observed secant
//! pairs (no extra DFT evaluation). The L-BFGS two-loop recursion uses
//! H0 = P^-1 = (mu P1)^-1. Solves use block-Jacobi preconditioned conjugate
//! gradients on the sparse 3x3-block matrix (no new dependency).
//!
//! Literature: Packwood et al., J. Chem. Phys. 144, 164109 (2016) (sparse
//! neighbour-graph preconditioner for L-BFGS with Armijo line search);
//! Nocedal and Wright, Numerical Optimization, ch. 3 and 7. The bond-spring
//! shape (stretch-dominant blocks, covalent-radius graph) is the project's own
//! A3 model of PRELAX-P1 (docs/roadmap/precon-relax-design.md section 4). No
//! code of other packages was used.
use serde::{Deserialize, Serialize};

/// Decay of the bond stiffness with the bond length, w = exp(-A (d/r0 - 1)).
/// A3 value (PRELAX-P1); Packwood et al. use 3 for the same functional form.
pub const STIFFNESS_DECAY: f64 = 2.5;
/// Fraction of the bond block along the bond: B = (1 - s) I + s n n^T. A
/// central pair force has a Hessian n n^T phi'' plus a transverse part
/// phi'/r (I - n n^T) that is small near equilibrium; 0.7 keeps a transverse
/// share so that bending (angular) stiffness of covalent networks is not zero.
pub const STRETCH_FRACTION: f64 = 0.7;
/// Isotropic stabiliser relative to a unit bond weight: removes the
/// translation null space of the graph Laplacian and pins isolated atoms and
/// disconnected fragments (A3 ratio 0.05/5 eV/A^2).
pub const STABILISER: f64 = 0.01;
/// Bond cutoff as a multiple of the covalent-radius sum.
pub const BOND_FACTOR: f64 = 1.3;
/// Cutoff multiple for a molecular element (H, C, N, O) bonded to a
/// non-molecular one: adsorbate-surface contacts are longer than covalent
/// bonds and would otherwise leave the molecule coupled only by the stabiliser.
pub const MIXED_BOND_FACTOR: f64 = 1.8;
/// L-BFGS memory (secant pairs kept).
pub const MEMORY: usize = 10;
/// Rebuild P when an atom moved more than this fraction of the shortest bond
/// at the last build (Packwood et al. use 0.5 of the nearest-neighbour distance).
pub const REBUILD_FRACTION: f64 = 0.1;
/// Armijo sufficient-decrease constant (Nocedal and Wright, section 3.1).
pub const ARMIJO_C1: f64 = 1e-4;
/// Line-search trials per direction before the fallback (or failure).
pub const MAX_TRIALS: usize = 3;
/// Initial spring of the native BFGS (k = 1 Ha/Bohr^2 on every Cartesian
/// coordinate, `optimizer::Optimizer::new`): the first P is scaled so that its
/// mean diagonal over free coordinates equals this spring until a secant pair
/// exists. Stiff by design (Si: about 3 to 10 times the measured scale), so
/// the first step is short and its pair fixes mu.
pub const PRIOR_SPRING_HA_BOHR2: f64 = 1.;

/// Covalent radii, Angstrom, Z = 1..96: Cordero et al., Dalton Trans. 2008,
/// 2832 (Table 2; C sp3, Mn/Fe/Co low spin).
const COVALENT_RADII_ANGSTROM: [f64; 96] = [
    0.31, 0.28, 1.28, 0.96, 0.84, 0.76, 0.71, 0.66, 0.57, 0.58, 1.66, 1.41, 1.21, 1.11, 1.07, 1.05,
    1.02, 1.06, 2.03, 1.76, 1.70, 1.60, 1.53, 1.39, 1.39, 1.32, 1.26, 1.24, 1.32, 1.22, 1.22, 1.20,
    1.19, 1.20, 1.20, 1.16, 2.20, 1.95, 1.90, 1.75, 1.64, 1.54, 1.47, 1.46, 1.42, 1.39, 1.45, 1.44,
    1.42, 1.39, 1.39, 1.38, 1.39, 1.40, 2.44, 2.15, 2.07, 2.04, 2.03, 2.01, 1.99, 1.98, 1.98, 1.96,
    1.94, 1.92, 1.92, 1.89, 1.90, 1.87, 1.87, 1.75, 1.70, 1.62, 1.51, 1.44, 1.41, 1.36, 1.36, 1.32,
    1.45, 1.46, 1.48, 1.40, 1.50, 1.50, 2.60, 2.21, 2.15, 2.06, 2.00, 1.96, 1.90, 1.87, 1.80, 1.69,
];

/// Per-atom data of the model: covalent radius (Bohr) and whether the element
/// is one of the molecular elements H, C, N, O.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AtomKind {
    pub radius_bohr: f64,
    pub molecular: bool,
}
impl AtomKind {
    /// From an element symbol (exact case).
    pub fn from_element(symbol: &str) -> Result<Self, String> {
        let z =
            crate::run::volumetric_out::atomic_number(symbol.trim()).map_err(|e| e.to_string())?;
        let radius = COVALENT_RADII_ANGSTROM.get(z - 1).ok_or_else(|| {
            format!("precon_lbfgs: no covalent radius for {symbol} (Z = {z} > 96)")
        })?;
        Ok(Self {
            radius_bohr: dft_units::to_bohr(*radius, dft_units::LengthUnit::Angstrom),
            molecular: matches!(z, 1 | 6 | 7 | 8),
        })
    }
}

/// Geometry-independent model data: cell (columns, Bohr), atom kinds and the
/// free-coordinate mask (3N).
#[derive(Clone, Debug)]
pub struct Model {
    pub lattice: [[f64; 3]; 3],
    pub kinds: Vec<AtomKind>,
    pub free: Vec<bool>,
}
impl Model {
    fn cutoff(&self, i: usize, j: usize) -> f64 {
        let (a, b) = (self.kinds[i], self.kinds[j]);
        let factor = if a.molecular != b.molecular {
            MIXED_BOND_FACTOR
        } else {
            BOND_FACTOR
        };
        factor * (a.radius_bohr + b.radius_bohr)
    }
    fn max_cutoff(&self) -> f64 {
        let r = self.kinds.iter().map(|k| k.radius_bohr).fold(0., f64::max);
        MIXED_BOND_FACTOR.max(BOND_FACTOR) * 2. * r
    }
}

fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn inverse3(m: &[[f64; 3]; 3]) -> Option<[[f64; 3]; 3]> {
    let c = |r: usize, s: usize| {
        let (r1, r2) = ((r + 1) % 3, (r + 2) % 3);
        let (s1, s2) = ((s + 1) % 3, (s + 2) % 3);
        m[r1][s1] * m[r2][s2] - m[r1][s2] * m[r2][s1]
    };
    let det = m[0][0] * c(0, 0) + m[0][1] * c(0, 1) + m[0][2] * c(0, 2);
    (det.is_finite() && det != 0.)
        .then(|| std::array::from_fn(|i| std::array::from_fn(|j| c(j, i) / det)))
}

/// Directed neighbour pairs (i, j, d) with d = r_j + T - r_i the Cartesian
/// vector to atom j in lattice image T, |d| < cutoff(i, j), excluding (i, i, 0).
/// Every unordered bond appears in both directions. Cell lists on wrapped
/// fractional coordinates: bins at least `max_cut` wide along each plane
/// normal, so neighbours lie in adjacent bins; cells thinner than the cutoff
/// scan as many images as needed. O(N) for a fixed density.
pub fn neighbours(
    lattice: &[[f64; 3]; 3],
    positions: &[f64],
    max_cut: f64,
    cutoff: impl Fn(usize, usize) -> f64,
) -> Vec<(usize, usize, [f64; 3])> {
    let n = positions.len() / 3;
    let inv = inverse3(lattice).expect("nonsingular lattice");
    let col = |k: usize| [lattice[0][k], lattice[1][k], lattice[2][k]];
    let volume = dot(&col(0), &cross(col(1), col(2))).abs();
    let frac: Vec<[f64; 3]> = positions
        .chunks(3)
        .map(|r| {
            std::array::from_fn(|k| {
                let f = inv[k][0] * r[0] + inv[k][1] * r[1] + inv[k][2] * r[2];
                f - f.floor()
            })
        })
        .collect();
    let per_axis_cap = ((n as f64).cbrt().ceil() as usize + 2).max(1);
    let mut bins = [1usize; 3];
    let mut reach = [1i64; 3];
    for k in 0..3 {
        let c = cross(col((k + 1) % 3), col((k + 2) % 3));
        let height = volume / dot(&c, &c).sqrt();
        bins[k] = ((height / max_cut).floor() as usize).clamp(1, per_axis_cap);
        reach[k] = ((max_cut / (height / bins[k] as f64)).ceil() as i64).max(1);
    }
    let bin_of = |f: &[f64; 3]| -> [usize; 3] {
        std::array::from_fn(|k| ((f[k] * bins[k] as f64).floor() as usize).min(bins[k] - 1))
    };
    let linear = |b: [usize; 3]| (b[0] * bins[1] + b[1]) * bins[2] + b[2];
    let mut members = vec![Vec::new(); bins[0] * bins[1] * bins[2]];
    for (i, f) in frac.iter().enumerate() {
        members[linear(bin_of(f))].push(i);
    }
    let mut pairs = Vec::new();
    for (i, fi) in frac.iter().enumerate() {
        let bi = bin_of(fi);
        for ox in -reach[0]..=reach[0] {
            for oy in -reach[1]..=reach[1] {
                for oz in -reach[2]..=reach[2] {
                    let o = [ox, oy, oz];
                    let mut target = [0usize; 3];
                    let mut image = [0i64; 3];
                    for k in 0..3 {
                        let c = bi[k] as i64 + o[k];
                        target[k] = c.rem_euclid(bins[k] as i64) as usize;
                        image[k] = c.div_euclid(bins[k] as i64);
                    }
                    for &j in &members[linear(target)] {
                        if i == j && image == [0; 3] {
                            continue;
                        }
                        let df: [f64; 3] =
                            std::array::from_fn(|k| frac[j][k] + image[k] as f64 - fi[k]);
                        let d: [f64; 3] = std::array::from_fn(|a| {
                            lattice[a][0] * df[0] + lattice[a][1] * df[1] + lattice[a][2] * df[2]
                        });
                        let r = dot(&d, &d).sqrt();
                        if r < cutoff(i, j) {
                            pairs.push((i, j, d));
                        }
                    }
                }
            }
        }
    }
    pairs
}

type Block = [[f64; 3]; 3];

/// Sparse symmetric P1 (3x3 blocks per atom pair) restricted to free coordinates.
#[derive(Clone, Debug)]
pub struct Matrix {
    row_start: Vec<usize>,
    cols: Vec<usize>,
    blocks: Vec<Block>,
    free: Vec<bool>,
    jacobi: Vec<Block>,
    /// Undirected bonds (image pairs) in the graph.
    pub bonds: usize,
    /// Shortest bond, Bohr (None: no bond).
    pub min_bond_bohr: Option<f64>,
    /// Positions at which this matrix was built.
    pub positions: Vec<f64>,
}
impl Matrix {
    pub fn build(model: &Model, positions: &[f64]) -> Self {
        let n = model.kinds.len();
        assert_eq!(positions.len(), 3 * n);
        assert_eq!(model.free.len(), 3 * n);
        let pairs = neighbours(&model.lattice, positions, model.max_cutoff(), |i, j| {
            model.cutoff(i, j)
        });
        let mut rows: Vec<Vec<(usize, Block)>> = vec![Vec::new(); n];
        let mut min_bond = f64::INFINITY;
        for &(i, j, d) in &pairs {
            let r = dot(&d, &d).sqrt();
            min_bond = min_bond.min(r);
            let r0 = model.kinds[i].radius_bohr + model.kinds[j].radius_bohr;
            let w = (-STIFFNESS_DECAY * (r / r0 - 1.)).exp();
            let u = [d[0] / r, d[1] / r, d[2] / r];
            let b: Block = std::array::from_fn(|a| {
                std::array::from_fn(|c| {
                    w * ((1. - STRETCH_FRACTION) * f64::from(u8::from(a == c))
                        + STRETCH_FRACTION * u[a] * u[c])
                })
            });
            let neg: Block = std::array::from_fn(|a| std::array::from_fn(|c| -b[a][c]));
            rows[i].push((i, b));
            rows[i].push((j, neg));
        }
        let mut stab = [[0.; 3]; 3];
        for (a, row) in stab.iter_mut().enumerate() {
            row[a] = STABILISER;
        }
        let (mut row_start, mut cols, mut blocks) = (vec![0], Vec::new(), Vec::new());
        let mut jacobi = Vec::with_capacity(n);
        for (i, mut row) in rows.into_iter().enumerate() {
            row.push((i, stab));
            // Deterministic merge: stable sort by column, then sum in order.
            row.sort_by_key(|(j, _)| *j);
            let start = cols.len();
            for (j, b) in row {
                if cols.len() > start && *cols.last().unwrap() == j {
                    let last: &mut Block = blocks.last_mut().unwrap();
                    for a in 0..3 {
                        for c in 0..3 {
                            last[a][c] += b[a][c];
                        }
                    }
                } else {
                    cols.push(j);
                    blocks.push(b);
                }
            }
            row_start.push(cols.len());
            let diag = (start..cols.len())
                .find(|&k| cols[k] == i)
                .map(|k| blocks[k])
                .expect("diagonal block");
            let masked: Block = std::array::from_fn(|a| {
                std::array::from_fn(|c| {
                    let (fa, fc) = (model.free[3 * i + a], model.free[3 * i + c]);
                    if fa && fc {
                        diag[a][c]
                    } else if a == c {
                        1.
                    } else {
                        0.
                    }
                })
            });
            jacobi.push(inverse3(&masked).expect("SPD diagonal block"));
        }
        Self {
            row_start,
            cols,
            blocks,
            free: model.free.clone(),
            jacobi,
            bonds: pairs.len() / 2,
            min_bond_bohr: min_bond.is_finite().then_some(min_bond),
            positions: positions.to_vec(),
        }
    }
    fn mask(&self, v: &mut [f64]) {
        for (x, f) in v.iter_mut().zip(&self.free) {
            if !f {
                *x = 0.;
            }
        }
    }
    /// P1 v on free coordinates (fixed components of v ignored, of the result zero).
    pub fn apply(&self, v: &[f64]) -> Vec<f64> {
        let mut out = vec![0.; v.len()];
        for i in 0..self.row_start.len() - 1 {
            let mut acc = [0.; 3];
            for k in self.row_start[i]..self.row_start[i + 1] {
                let (j, b) = (self.cols[k], &self.blocks[k]);
                for (c, fc) in (0..3).map(|c| (c, self.free[3 * j + c])) {
                    if fc {
                        let x = v[3 * j + c];
                        for (a, acc) in acc.iter_mut().enumerate() {
                            *acc += b[a][c] * x;
                        }
                    }
                }
            }
            out[3 * i..3 * i + 3].copy_from_slice(&acc);
        }
        self.mask(&mut out);
        out
    }
    /// Mean diagonal of P1 over free coordinates.
    pub fn mean_free_diagonal(&self) -> f64 {
        let (mut sum, mut count) = (0., 0usize);
        for i in 0..self.row_start.len() - 1 {
            for k in self.row_start[i]..self.row_start[i + 1] {
                if self.cols[k] == i {
                    for a in 0..3 {
                        if self.free[3 * i + a] {
                            sum += self.blocks[k][a][a];
                            count += 1;
                        }
                    }
                }
            }
        }
        if count == 0 { 1. } else { sum / count as f64 }
    }
    /// P1^-1 b on free coordinates: block-Jacobi preconditioned conjugate
    /// gradients to a relative residual of 1e-13 (or 10 n_free + 100
    /// iterations). Deterministic (fixed summation order).
    pub fn solve(&self, b: &[f64]) -> Vec<f64> {
        let mut r = b.to_vec();
        self.mask(&mut r);
        let mut x = vec![0.; b.len()];
        let norm_b = dot(&r, &r).sqrt();
        if norm_b == 0. || !norm_b.is_finite() {
            return x;
        }
        let jacobi = |r: &[f64]| -> Vec<f64> {
            let mut z = vec![0.; r.len()];
            for (i, m) in self.jacobi.iter().enumerate() {
                for a in 0..3 {
                    z[3 * i + a] = (0..3).map(|c| m[a][c] * r[3 * i + c]).sum();
                }
            }
            self.mask(&mut z);
            z
        };
        let n_free = self.free.iter().filter(|f| **f).count();
        let mut z = jacobi(&r);
        let mut p = z.clone();
        let mut rz = dot(&r, &z);
        for _ in 0..10 * n_free + 100 {
            let q = self.apply(&p);
            let pq = dot(&p, &q);
            if pq <= 0. || !pq.is_finite() {
                break;
            }
            let alpha = rz / pq;
            for k in 0..x.len() {
                x[k] += alpha * p[k];
                r[k] -= alpha * q[k];
            }
            if dot(&r, &r).sqrt() <= 1e-13 * norm_b {
                break;
            }
            z = jacobi(&r);
            let rz_new = dot(&r, &z);
            let beta = rz_new / rz;
            rz = rz_new;
            for k in 0..p.len() {
                p[k] = z[k] + beta * p[k];
            }
        }
        x
    }
    /// Whether an atom moved more than REBUILD_FRACTION of the shortest bond
    /// (or, without bonds, of the smallest covalent-radius sum) since the build.
    pub fn needs_rebuild(&self, model: &Model, positions: &[f64]) -> bool {
        let scale = self.min_bond_bohr.unwrap_or_else(|| {
            2. * model
                .kinds
                .iter()
                .map(|k| k.radius_bohr)
                .fold(f64::INFINITY, f64::min)
        });
        let moved = positions
            .chunks(3)
            .zip(self.positions.chunks(3))
            .map(|(a, b)| {
                let d: Vec<f64> = a.iter().zip(b).map(|(x, y)| x - y).collect();
                dot(&d, &d).sqrt()
            })
            .fold(0., f64::max);
        moved > REBUILD_FRACTION * scale
    }
}

/// Optimizer memory persisted in restart.json (positions of the last P build,
/// scale, secant pairs and counters). P itself is rebuilt from
/// `build_positions` on resume (deterministic).
#[derive(Clone, Debug, Serialize, Deserialize, PartialEq)]
#[serde(deny_unknown_fields)]
pub struct State {
    /// Scale of P, Ha/Bohr^2.
    pub mu: f64,
    /// "prior" until the first stored secant pair, then "secant".
    pub mu_source: String,
    pub s: Vec<Vec<f64>>,
    pub y: Vec<Vec<f64>>,
    pub build_positions: Vec<f64>,
    pub rebuilds: usize,
    pub history_resets: usize,
    /// Every attempted SCF evaluation (accepted, rejected, failed), initial included.
    pub evaluations: usize,
    pub rejected_trials: usize,
}
impl State {
    pub fn new(matrix: &Matrix) -> Self {
        Self {
            mu: PRIOR_SPRING_HA_BOHR2 / matrix.mean_free_diagonal(),
            mu_source: "prior".into(),
            s: vec![],
            y: vec![],
            build_positions: matrix.positions.clone(),
            rebuilds: 0,
            history_resets: 0,
            evaluations: 0,
            rejected_trials: 0,
        }
    }
    /// H g with the two-loop recursion and H0 = (mu P1)^-1.
    pub fn inverse_hessian_apply(&self, matrix: &Matrix, g: &[f64]) -> Vec<f64> {
        let mut q = g.to_vec();
        let mut alphas = Vec::with_capacity(self.s.len());
        for (s, y) in self.s.iter().zip(&self.y).rev() {
            let a = dot(s, &q) / dot(y, s);
            for (q, y) in q.iter_mut().zip(y) {
                *q -= a * y;
            }
            alphas.push(a);
        }
        let mut z = matrix.solve(&q);
        for v in &mut z {
            *v /= self.mu;
        }
        for ((s, y), a) in self.s.iter().zip(&self.y).zip(alphas.iter().rev()) {
            let b = dot(y, &z) / dot(y, s);
            for (z, s) in z.iter_mut().zip(s) {
                *z += (a - b) * s;
            }
        }
        z
    }
    /// Store a secant pair (s = x_new - x, y = g_new - g) when the curvature
    /// s.y > 1e-12 |s||y| holds, keep the last MEMORY pairs, and set
    /// mu = y.P1^-1 y / s.y (the P-metric form of the L-BFGS scaling
    /// gamma = s.y / y.y). Returns whether the pair was stored.
    pub fn store(&mut self, matrix: &Matrix, s: Vec<f64>, y: Vec<f64>) -> bool {
        let sy = dot(&s, &y);
        if !(sy > 1e-12 * (dot(&s, &s) * dot(&y, &y)).sqrt() && sy > 0. && sy.is_finite()) {
            return false;
        }
        let mu = dot(&y, &matrix.solve(&y)) / sy;
        if mu.is_finite() && mu > 0. {
            self.mu = mu;
            self.mu_source = "secant".into();
        }
        self.s.push(s);
        self.y.push(y);
        if self.s.len() > MEMORY {
            self.s.remove(0);
            self.y.remove(0);
        }
        true
    }
    fn reset_history(&mut self) {
        self.s.clear();
        self.y.clear();
        self.history_resets += 1;
    }
}

/// Largest per-atom Euclidean norm of a flat 3N vector.
pub fn max_atom_norm(v: &[f64]) -> f64 {
    v.chunks(3).map(|a| dot(a, a).sqrt()).fold(0., f64::max)
}

/// One ionic step's line search along a fixed direction `p` (already capped):
/// trial positions x0 + t p.
#[derive(Clone, Debug)]
pub struct LineSearch {
    pub x0: Vec<f64>,
    pub e0: f64,
    pub g0: Vec<f64>,
    pub p: Vec<f64>,
    pub gp: f64,
    pub t: f64,
    /// Successful trials judged on this direction.
    pub trials: usize,
    /// Trials on all directions of this step (the step's line-search count).
    pub total_trials: usize,
    pub fallback: bool,
    cap: f64,
}
/// Decision after a trial evaluation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Judgement {
    Accept,
    /// Try again at the new `t` (same or fallback direction).
    Retry,
    /// No acceptable point after the fallback direction.
    Fail,
}
fn capped(mut p: Vec<f64>, cap: f64) -> Vec<f64> {
    let m = max_atom_norm(&p);
    if m > cap {
        for v in &mut p {
            *v *= cap / m;
        }
    }
    p
}
impl LineSearch {
    /// Start a step at x0 with energy e0 and projected gradient g0 (= -force):
    /// p = -H g0 capped to `cap` per atom; a non-descent or non-finite
    /// direction resets the history and uses -P^-1 g0 (preconditioned
    /// steepest descent).
    pub fn begin(
        state: &mut State,
        matrix: &Matrix,
        x0: &[f64],
        e0: f64,
        g0: &[f64],
        cap: f64,
    ) -> Self {
        let mut p: Vec<f64> = state
            .inverse_hessian_apply(matrix, g0)
            .iter()
            .map(|v| -v)
            .collect();
        let mut fallback = false;
        if dot(g0, &p).is_nan() || dot(g0, &p) >= 0. || p.iter().any(|v| !v.is_finite()) {
            state.reset_history();
            p = steepest(state, matrix, g0);
            fallback = true;
        }
        let p = capped(p, cap);
        Self {
            x0: x0.to_vec(),
            e0,
            g0: g0.to_vec(),
            gp: dot(g0, &p),
            p,
            t: 1.,
            trials: 0,
            total_trials: 0,
            fallback,
            cap,
        }
    }
    /// Trial positions at the current t.
    pub fn trial(&self) -> Vec<f64> {
        self.x0
            .iter()
            .zip(&self.p)
            .map(|(x, p)| if *p == 0. { *x } else { x + self.t * p })
            .collect()
    }
    /// The current step t p.
    pub fn step(&self) -> Vec<f64> {
        self.p.iter().map(|p| self.t * p).collect()
    }
    /// An SCF failure at the current trial: halve t (the retry runs from a
    /// fresh start, as the BFGS driver does).
    pub fn failed(&mut self) {
        self.t *= 0.5;
        self.total_trials += 1;
    }
    /// Judge a successful trial with energy `e` and projected gradient `g`.
    /// `stop` is true when the trial already meets the stop rule (accepted
    /// regardless of the Armijo test). `noise` is the energy tolerance (Ha)
    /// of the Armijo test (SCF energy accuracy).
    pub fn judge(
        &mut self,
        state: &mut State,
        matrix: &Matrix,
        e: f64,
        g: &[f64],
        stop: bool,
        noise: f64,
    ) -> Judgement {
        self.trials += 1;
        self.total_trials += 1;
        if stop || e <= self.e0 + ARMIJO_C1 * self.t * self.gp + noise {
            return Judgement::Accept;
        }
        state.rejected_trials += 1;
        if self.trials < MAX_TRIALS {
            self.t = backtrack(self.t, self.e0, e, self.gp, Some(dot(g, &self.p)));
            return Judgement::Retry;
        }
        if self.fallback {
            return Judgement::Fail;
        }
        state.reset_history();
        let p = steepest(state, matrix, &self.g0);
        self.p = capped(p, self.cap);
        self.gp = dot(&self.g0, &self.p);
        if self.gp.is_nan() || self.gp >= 0. {
            return Judgement::Fail;
        }
        self.t = 1.;
        self.trials = 0;
        self.fallback = true;
        Judgement::Retry
    }
}
fn steepest(state: &State, matrix: &Matrix, g: &[f64]) -> Vec<f64> {
    matrix.solve(g).iter().map(|v| -v / state.mu).collect()
}
/// Next step length after a rejected trial at t (energy e0 -> et, directional
/// derivatives gp at 0 and gpt at t): the zero of the linear interpolation of
/// the directional derivative when the curvature along p is positive (exact on
/// a quadratic), else the minimiser of the quadratic through e0, gp and et;
/// safeguarded to [0.1 t, 0.5 t].
pub fn backtrack(t: f64, e0: f64, et: f64, gp: f64, gpt: Option<f64>) -> f64 {
    let curvature = gpt.map(|d| (d - gp) / t);
    let candidate = match curvature {
        Some(c) if c > 0. && c.is_finite() => -gp / c,
        _ => {
            let denom = 2. * (et - e0 - gp * t);
            if denom > 0. {
                -gp * t * t / denom
            } else {
                0.5 * t
            }
        }
    };
    let candidate = if candidate.is_finite() {
        candidate
    } else {
        0.5 * t
    };
    candidate.clamp(0.1 * t, 0.5 * t)
}

#[cfg(test)]
mod tests {
    //! Units: Bohr, Ha, Ha/Bohr^2. Tolerances are stated per assertion.
    use super::*;

    /// Deterministic generator (64-bit LCG, Knuth MMIX constants) in [0, 1).
    struct Lcg(u64);
    impl Lcg {
        fn next(&mut self) -> f64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            (self.0 >> 11) as f64 / (1u64 << 53) as f64
        }
    }
    fn si_kind() -> AtomKind {
        AtomKind::from_element("Si").unwrap()
    }
    fn random_cell(rng: &mut Lcg, n: usize, length: f64) -> ([[f64; 3]; 3], Vec<f64>) {
        // Skewed cell: columns a1, a2, a3 with random shear.
        let mut l = [[0.; 3]; 3];
        for (k, row) in l.iter_mut().enumerate() {
            for (j, x) in row.iter_mut().enumerate() {
                *x = if k == j {
                    length
                } else {
                    0.3 * length * (rng.next() - 0.5)
                };
            }
        }
        let positions: Vec<f64> = (0..n)
            .flat_map(|_| {
                let f = [rng.next(), rng.next(), rng.next()];
                (0..3)
                    .map(|a| (0..3).map(|k| l[a][k] * f[k]).sum::<f64>())
                    .collect::<Vec<_>>()
            })
            .collect();
        (l, positions)
    }
    /// Brute-force neighbour list: all images within a wide range.
    fn brute(
        lattice: &[[f64; 3]; 3],
        positions: &[f64],
        cut: f64,
    ) -> Vec<(usize, usize, [i64; 3])> {
        let n = positions.len() / 3;
        let mut out = Vec::new();
        for i in 0..n {
            for j in 0..n {
                for a in -4i64..=4 {
                    for b in -4i64..=4 {
                        for c in -4i64..=4 {
                            if i == j && [a, b, c] == [0; 3] {
                                continue;
                            }
                            let t = [a as f64, b as f64, c as f64];
                            let d: Vec<f64> = (0..3)
                                .map(|x| {
                                    positions[3 * j + x] - positions[3 * i + x]
                                        + (0..3).map(|k| lattice[x][k] * t[k]).sum::<f64>()
                                })
                                .collect();
                            if dot(&d, &d).sqrt() < cut {
                                out.push((i, j, [a, b, c]));
                            }
                        }
                    }
                }
            }
        }
        out
    }
    #[test]
    fn neighbour_list_matches_brute_force_including_thin_cells() {
        let mut rng = Lcg(7);
        for (n, length) in [(5, 4.0), (12, 9.0), (30, 14.0), (3, 2.5)] {
            let (l, x) = random_cell(&mut rng, n, length);
            let cut = 5.2;
            let fast = neighbours(&l, &x, cut, |_, _| cut);
            let slow = brute(&l, &x, cut);
            assert_eq!(fast.len(), slow.len(), "n {n} length {length}");
            // Same multiset of (i, j, |d|) to 1e-12 Bohr.
            let key = |v: &mut Vec<(usize, usize, f64)>| {
                v.sort_by(|a, b| (a.0, a.1).cmp(&(b.0, b.1)).then(a.2.total_cmp(&b.2)))
            };
            let mut a: Vec<_> = fast
                .iter()
                .map(|(i, j, d)| (*i, *j, dot(d, d).sqrt()))
                .collect();
            let mut b: Vec<_> = slow
                .iter()
                .map(|(i, j, t)| {
                    let tt = [t[0] as f64, t[1] as f64, t[2] as f64];
                    let d: Vec<f64> = (0..3)
                        .map(|c| {
                            x[3 * j + c] - x[3 * i + c]
                                + (0..3).map(|k| l[c][k] * tt[k]).sum::<f64>()
                        })
                        .collect();
                    (*i, *j, dot(&d, &d).sqrt())
                })
                .collect();
            key(&mut a);
            key(&mut b);
            for (p, q) in a.iter().zip(&b) {
                assert_eq!((p.0, p.1), (q.0, q.1));
                assert!((p.2 - q.2).abs() < 1e-12);
            }
        }
    }
    fn model(l: [[f64; 3]; 3], n: usize, fixed: &[(usize, [bool; 3])]) -> Model {
        let mut free = vec![true; 3 * n];
        for (atom, comps) in fixed {
            for a in 0..3 {
                if comps[a] {
                    free[3 * atom + a] = false;
                }
            }
        }
        // Mixed kinds: Si and O (molecular) to exercise both cutoffs.
        let o = AtomKind::from_element("O").unwrap();
        Model {
            lattice: l,
            kinds: (0..n)
                .map(|i| if i % 3 == 0 { o } else { si_kind() })
                .collect(),
            free,
        }
    }
    fn random_vector(rng: &mut Lcg, m: &Model) -> Vec<f64> {
        m.free
            .iter()
            .map(|f| if *f { rng.next() - 0.5 } else { 0. })
            .collect()
    }
    #[test]
    fn spd_symmetric_and_inverse_accurate_with_images_and_fixed_atoms() {
        let mut rng = Lcg(11);
        for (n, length, fixed) in [
            (8, 7.0, vec![]),
            (4, 4.5, vec![(0usize, [true; 3])]),
            (20, 12.0, vec![(1, [true, false, true]), (5, [true; 3])]),
            (2, 3.0, vec![(1, [false, true, false])]),
        ] {
            let (l, x) = random_cell(&mut rng, n, length);
            let m = model(l, n, &fixed);
            let p = Matrix::build(&m, &x);
            assert!(p.bonds > 0, "test cells must contain bonds");
            for _ in 0..5 {
                let u = random_vector(&mut rng, &m);
                let v = random_vector(&mut rng, &m);
                let (pu, pv) = (p.apply(&u), p.apply(&v));
                // Symmetry to rounding.
                assert!((dot(&u, &pv) - dot(&v, &pu)).abs() < 1e-12 * (1. + dot(&u, &pv).abs()));
                // Lower bound from the stabiliser: v.P1 v >= c_stab |v|^2.
                assert!(dot(&v, &pv) >= STABILISER * dot(&v, &v) * (1. - 1e-12));
                // Inverse accuracy |P1 P1^-1 v - v| < 1e-10 |v|.
                let back = p.apply(&p.solve(&v));
                let err: f64 = back
                    .iter()
                    .zip(&v)
                    .map(|(a, b)| (a - b).powi(2))
                    .sum::<f64>();
                assert!(err.sqrt() < 1e-10 * dot(&v, &v).sqrt(), "{}", err.sqrt());
                // Fixed components stay zero.
                for (k, f) in m.free.iter().enumerate() {
                    if !f {
                        assert_eq!(pv[k], 0.);
                    }
                }
            }
        }
    }
    #[test]
    fn invariant_under_rigid_translation_and_permutation() {
        let mut rng = Lcg(5);
        let n = 10;
        let (l, x) = random_cell(&mut rng, n, 8.0);
        let m = model(l, n, &[]);
        let p = Matrix::build(&m, &x);
        // Rigid translation (also by lattice vectors: wrapping changes nothing).
        let shift = [0.37, -2.1, 11.3];
        let xt: Vec<f64> = x
            .iter()
            .enumerate()
            .map(|(k, v)| v + shift[k % 3])
            .collect();
        let pt = Matrix::build(&m, &xt);
        let v = random_vector(&mut rng, &m);
        let (a, b) = (p.apply(&v), pt.apply(&v));
        for (a, b) in a.iter().zip(&b) {
            assert!((a - b).abs() < 1e-12, "{a} {b}");
        }
        // Permutation of the atoms (kinds permuted with them).
        let perm: Vec<usize> = (0..n).rev().collect();
        let mp = Model {
            kinds: perm.iter().map(|&i| m.kinds[i]).collect(),
            ..m.clone()
        };
        let xp: Vec<f64> = perm
            .iter()
            .flat_map(|&i| x[3 * i..3 * i + 3].to_vec())
            .collect();
        let vp: Vec<f64> = perm
            .iter()
            .flat_map(|&i| v[3 * i..3 * i + 3].to_vec())
            .collect();
        let pp = Matrix::build(&mp, &xp).apply(&vp);
        for (k, &i) in perm.iter().enumerate() {
            for c in 0..3 {
                assert!((pp[3 * k + c] - a[3 * i + c]).abs() < 1e-12);
            }
        }
    }
    #[test]
    fn disconnected_fragment_and_isolated_atom_stay_spd() {
        // A molecule (two O far apart from a Si pair) in a large box: three
        // components, every one pinned only by the stabiliser.
        let l = [[30., 0., 0.], [0., 30., 0.], [0., 0., 30.]];
        let x = vec![
            0., 0., 0., 4.4, 0., 0., // Si-Si bond
            15., 15., 15., 17.3, 15., 15., // O-O
            7., 22., 3., // isolated Si
        ];
        let m = Model {
            lattice: l,
            kinds: vec![
                si_kind(),
                si_kind(),
                AtomKind::from_element("O").unwrap(),
                AtomKind::from_element("O").unwrap(),
                si_kind(),
            ],
            free: vec![true; 15],
        };
        let p = Matrix::build(&m, &x);
        assert_eq!(p.bonds, 2);
        // The isolated atom is a pure stabiliser block; translation of each
        // fragment has curvature exactly c_stab.
        let mut v = vec![0.; 15];
        v[12] = 1.;
        assert!((dot(&v, &p.apply(&v)) - STABILISER).abs() < 1e-15);
        let mut t = vec![0.; 15];
        t[6] = 1.;
        t[9] = 1.;
        assert!((dot(&t, &p.apply(&t)) - 2. * STABILISER).abs() < 1e-14);
        let back = p.apply(&p.solve(&t));
        for (a, b) in back.iter().zip(&t) {
            assert!((a - b).abs() < 1e-10);
        }
    }
    #[test]
    fn exact_on_a_quadratic_whose_hessian_is_mu_p1() {
        // E = 1/2 (x - x*).(mu P1)(x - x*), mu = 0.05 Ha/Bohr^2 (unknown to the
        // optimizer, whose prior is stiffer). The first step has the exact
        // direction; its pair fixes mu, the second step lands on x*.
        let mut rng = Lcg(3);
        let n = 8;
        let (l, x_star) = random_cell(&mut rng, n, 7.5);
        let m = model(l, n, &[(2, [true; 3])]);
        let p1 = Matrix::build(&m, &x_star);
        let mu_true = 0.05;
        let energy_grad = |x: &[f64]| {
            let d: Vec<f64> = x.iter().zip(&x_star).map(|(a, b)| a - b).collect();
            let g: Vec<f64> = p1.apply(&d).iter().map(|v| mu_true * v).collect();
            (0.5 * dot(&d, &g), g)
        };
        let mut x = x_star.clone();
        let dx = random_vector(&mut rng, &m);
        for (x, d) in x.iter_mut().zip(&dx) {
            *x += 0.2 * d;
        }
        let mut state = State::new(&p1);
        let (mut e, mut g) = energy_grad(&x);
        let mut accepted = 0;
        let mut evaluations = 0;
        while max_atom_norm(&g) > 1e-10 && accepted < 5 {
            let mut ls = LineSearch::begin(&mut state, &p1, &x, e, &g, 10.);
            loop {
                let xt = ls.trial();
                let (et, gt) = energy_grad(&xt);
                evaluations += 1;
                match ls.judge(&mut state, &p1, et, &gt, false, 0.) {
                    Judgement::Accept => {
                        let s = ls.step();
                        let y: Vec<f64> = gt.iter().zip(&g).map(|(a, b)| a - b).collect();
                        state.store(&p1, s, y);
                        x = xt;
                        e = et;
                        g = gt;
                        accepted += 1;
                        break;
                    }
                    Judgement::Retry => {}
                    Judgement::Fail => panic!("line search failed"),
                }
            }
        }
        let err = x
            .iter()
            .zip(&x_star)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        assert!(
            accepted <= 2 && evaluations <= 3,
            "{accepted} {evaluations}"
        );
        assert!(err < 1e-9, "{err}");
        assert!((state.mu - mu_true).abs() < 1e-9 * mu_true, "{}", state.mu);
    }
    #[test]
    fn backtracking_is_exact_on_a_parabola_and_safeguarded() {
        // f(t) = -t + t^2 / 2 (minimum at 1): a trial at t = 4 is rejected; the
        // derivative interpolation returns exactly 1 (inside [0.4, 2]).
        let (gp, t) = (-1., 4.);
        let et = -t + t * t / 2.;
        let gpt = -1. + t;
        assert_eq!(backtrack(t, 0., et, gp, Some(gpt)), 1.);
        // Energy-only quadratic interpolation gives the same point.
        assert!((backtrack(t, 0., et, gp, None) - 1.).abs() < 1e-15);
        // Negative curvature along p: energy-only quadratic, 1/(2 (0.5 + 1)).
        assert!((backtrack(1., 0., 0.5, -1., Some(-2.)) - 1. / 3.).abs() < 1e-15);
        // Far overshoot: clamped at 0.1 t.
        assert_eq!(backtrack(100., 0., 4900., -1., Some(99.)), 10.);
    }
    /// Periodic Morse crystal (Cu-like parameters in atomic units: D = 0.0125
    /// Ha, alpha = 0.72 / Bohr, r_e = 5 Bohr) with a cosine taper between
    /// 1.15 and 1.35 r_e (nearest neighbours only, so fcc at a = sqrt(2) r_e
    /// is an unstressed minimum), so energy and forces are smooth and consistent.
    struct Morse {
        lattice: [[f64; 3]; 3],
    }
    const RE: f64 = 5.;
    impl Morse {
        fn pair(r: f64) -> (f64, f64) {
            let (d, a) = (0.0125, 0.72);
            let e = (-a * (r - RE)).exp();
            let (phi, dphi) = (d * ((1. - e).powi(2) - 1.), 2. * d * a * e * (1. - e));
            let (r1, rc) = (1.15 * RE, 1.35 * RE);
            if r < r1 {
                (phi, dphi)
            } else if r < rc {
                let u = std::f64::consts::PI * (r - r1) / (rc - r1);
                let s = 0.5 * (1. + u.cos());
                let ds = -0.5 * std::f64::consts::PI / (rc - r1) * u.sin();
                (phi * s, dphi * s + phi * ds)
            } else {
                (0., 0.)
            }
        }
        fn eval(&self, x: &[f64]) -> (f64, Vec<f64>) {
            let rc = 1.35 * RE;
            let mut e = 0.;
            let mut g = vec![0.; x.len()];
            for (i, _, d) in neighbours(&self.lattice, x, rc, |_, _| rc) {
                let r = dot(&d, &d).sqrt();
                let (phi, dphi) = Self::pair(r);
                e += 0.5 * phi;
                for a in 0..3 {
                    g[3 * i + a] -= dphi * d[a] / r;
                }
            }
            (e, g)
        }
    }
    /// fcc 2x2x2 conventional cells at the Morse nearest-neighbour distance,
    /// atom 0 removed (vacancy), every atom displaced uniformly in +-amp Bohr.
    fn morse_start(amp: f64, seed: u64) -> (Morse, Model, Vec<f64>) {
        let a = RE * std::f64::consts::SQRT_2;
        let l = [[2. * a, 0., 0.], [0., 2. * a, 0.], [0., 0., 2. * a]];
        let mut rng = Lcg(seed);
        let mut x = Vec::new();
        for c in 0..8usize {
            let o = [(c & 1) as f64, ((c >> 1) & 1) as f64, ((c >> 2) & 1) as f64];
            for b in [[0., 0., 0.], [0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]] {
                for k in 0..3 {
                    x.push(a * (o[k] + b[k]) + amp * (2. * rng.next() - 1.));
                }
            }
        }
        x.drain(0..3);
        let n = x.len() / 3;
        let kind = AtomKind {
            radius_bohr: RE / 2.,
            molecular: false,
        };
        let model = Model {
            lattice: l,
            kinds: vec![kind; n],
            free: vec![true; 3 * n],
        };
        (Morse { lattice: l }, model, x)
    }
    /// Evaluations to max atomic |g| < tol with the preconditioned L-BFGS;
    /// `rebuild` = Some(keep_pairs) rebuilds P after every accepted step.
    fn precon_evaluations(seed: u64, tol: f64, rebuild: Option<bool>) -> usize {
        let (pot, m, mut x) = morse_start(0.4, seed);
        let mut p = Matrix::build(&m, &x);
        let mut state = State::new(&p);
        let (mut e, mut g) = pot.eval(&x);
        let mut evaluations = 1;
        while max_atom_norm(&g) >= tol {
            assert!(evaluations < 500, "no convergence");
            let mut ls = LineSearch::begin(&mut state, &p, &x, e, &g, 0.3);
            loop {
                let xt = ls.trial();
                let (et, gt) = pot.eval(&xt);
                evaluations += 1;
                let stop = max_atom_norm(&gt) < tol;
                match ls.judge(&mut state, &p, et, &gt, stop, 1e-12) {
                    Judgement::Accept => {
                        let y: Vec<f64> = gt.iter().zip(&g).map(|(a, b)| a - b).collect();
                        state.store(&p, ls.step(), y);
                        (x, e, g) = (xt, et, gt);
                        break;
                    }
                    Judgement::Retry => {}
                    Judgement::Fail => panic!("line search failed"),
                }
            }
            if let Some(keep) = rebuild {
                p = Matrix::build(&m, &x);
                if !keep {
                    state.reset_history();
                }
            } else if p.needs_rebuild(&m, &x) {
                // Production rule (REBUILD_FRACTION), pairs kept.
                p = Matrix::build(&m, &x);
            }
        }
        evaluations
    }
    /// The same start with the native dense BFGS (H0 = I Bohr^2/Ha, trust
    /// radius halved on an energy rise and grown by 1.1 otherwise, as relax.rs).
    fn bfgs_evaluations(seed: u64, tol: f64) -> usize {
        let (pot, _, mut x) = morse_start(0.4, seed);
        let mut opt = crate::relax::optimizer::Optimizer::new(x.len(), 0.3);
        let (mut e, mut g) = pot.eval(&x);
        let mut evaluations = 1;
        while max_atom_norm(&g) >= tol {
            assert!(evaluations < 2000, "no convergence");
            let f: Vec<f64> = g.iter().map(|v| -v).collect();
            let s = opt.step(&x, &f, false, 0.3);
            for (x, s) in x.iter_mut().zip(&s) {
                *x += s;
            }
            let (en, gn) = pot.eval(&x);
            evaluations += 1;
            opt.radius = if en - e > 1e-12 {
                opt.radius * 0.5
            } else {
                (opt.radius * 1.1).min(0.3)
            };
            (e, g) = (en, gn);
        }
        evaluations
    }
    #[test]
    fn morse_vacancy_needs_fewer_evaluations_than_native_bfgs() {
        // Observation printed for the ADR; the assertion is only the direction.
        let tol = 1e-5;
        for seed in 1..=3 {
            let (pre, keep, reset, bfgs) = (
                precon_evaluations(seed, tol, None),
                precon_evaluations(seed, tol, Some(true)),
                precon_evaluations(seed, tol, Some(false)),
                bfgs_evaluations(seed, tol),
            );
            println!(
                "morse fcc vacancy seed {seed}: precon {pre} (rebuild every step: keep pairs {keep}, reset {reset}), native BFGS {bfgs}"
            );
            assert!(pre < bfgs, "seed {seed}: {pre} vs {bfgs}");
        }
    }
    #[test]
    fn element_data_and_units() {
        let si = AtomKind::from_element("Si").unwrap();
        // 1.11 Angstrom in Bohr.
        let expected = dft_units::to_bohr(1.11, dft_units::LengthUnit::Angstrom);
        assert_eq!(si.radius_bohr, expected);
        assert!(!si.molecular);
        assert!(AtomKind::from_element("H").unwrap().molecular);
        assert!(AtomKind::from_element("Bk").is_err());
        assert!(AtomKind::from_element("Xx").is_err());
    }
}
