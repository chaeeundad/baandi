//! Independent high-accuracy Ewald reference (INIT-EWALD-INDEPENDENT-REFERENCE,
//! 2026-11-18). Test code only, never shipped in the engine.
//!
//! It shares no code with the engine evaluators (`Scheme::Reference`,
//! `Scheme::Balanced`) or with `old_reference.rs`: its own lattice algebra
//! (volume in double-double), its own image enumeration (per-axis box from the
//! reciprocal vector norms, exact sphere predicate), its own cutoffs (eta r_c =
//! g_c / (2 eta) = t with t = 7.5 by default, analytic tails far below 1e-16 Ha)
//! and its own reciprocal half space with phases reduced modulo 1 in
//! double-double before the sine and cosine. Every term (real-space images,
//! structure factor components, |S|^2, reciprocal terms, self and background)
//! is accumulated in double-double arithmetic, so the only rounding left is the
//! relative rounding of each individual term (a few ulp, uncorrelated), never
//! the growth of a running sum. Same convention as the engine: the G = 0 term
//! is dropped, the neutralizing background is -pi Q^2 / (2 Omega eta^2), the
//! self term -eta/sqrt(pi) sum q^2.
//!
//! The evaluator is validated against analytic Madelung constants (NaCl, CsCl)
//! and against exact supercell relations of the one-component lattice with
//! background (sc, fcc, bcc), its forces and stress against central finite
//! differences of its own energy; its eta independence measures its own error.
// Cartesian and lattice axes are indexed explicitly, as in the formulas.
#![allow(clippy::needless_range_loop)]
use super::super::{EwaldOptions, EwaldPlan, energy_with, forces, set_scheme_override, stress};
use super::super::{Ion, Scheme};
use dft_geometry::Cell;
use std::f64::consts::{PI, TAU};

include!("benchmark_positions.rs");

// ---------------------------------------------------------------- double-double

/// Unevaluated sum hi + lo with |lo| <= ulp(hi)/2 (about 106 significant bits).
#[derive(Clone, Copy, Debug, Default, PartialEq)]
struct Dd {
    hi: f64,
    lo: f64,
}
fn two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    let bb = s - a;
    (s, (a - (s - bb)) + (b - bb))
}
fn quick_two_sum(a: f64, b: f64) -> (f64, f64) {
    let s = a + b;
    (s, b - (s - a))
}
/// Exact product a b as a double-double (fused multiply-add).
fn dd_prod(a: f64, b: f64) -> Dd {
    let p = a * b;
    Dd {
        hi: p,
        lo: a.mul_add(b, -p),
    }
}
impl Dd {
    const ZERO: Dd = Dd { hi: 0.0, lo: 0.0 };
    /// pi to double-double precision.
    const PI: Dd = Dd {
        hi: PI,
        lo: 1.2246467991473532e-16,
    };
    fn new(x: f64) -> Self {
        Dd { hi: x, lo: 0.0 }
    }
    fn add(self, o: Dd) -> Dd {
        let (s, e) = two_sum(self.hi, o.hi);
        let (t, f) = two_sum(self.lo, o.lo);
        let (s, e) = quick_two_sum(s, e + t);
        let (hi, lo) = quick_two_sum(s, e + f);
        Dd { hi, lo }
    }
    fn add_f(self, x: f64) -> Dd {
        let (s, e) = two_sum(self.hi, x);
        let (hi, lo) = quick_two_sum(s, e + self.lo);
        Dd { hi, lo }
    }
    fn neg(self) -> Dd {
        Dd {
            hi: -self.hi,
            lo: -self.lo,
        }
    }
    fn sub(self, o: Dd) -> Dd {
        self.add(o.neg())
    }
    fn mul(self, o: Dd) -> Dd {
        let p = dd_prod(self.hi, o.hi);
        let (hi, lo) = quick_two_sum(p.hi, p.lo + (self.hi * o.lo + self.lo * o.hi));
        Dd { hi, lo }
    }
    fn mul_f(self, x: f64) -> Dd {
        self.mul(Dd::new(x))
    }
    fn div(self, o: Dd) -> Dd {
        let q1 = self.hi / o.hi;
        let r = self.sub(o.mul_f(q1));
        let q2 = r.hi / o.hi;
        let r = r.sub(o.mul_f(q2));
        let q3 = r.hi / o.hi;
        let (hi, lo) = quick_two_sum(q1, q2);
        Dd { hi, lo }.add_f(q3)
    }
    fn sqrt(self) -> Dd {
        let x = self.hi.sqrt();
        let r = self.sub(dd_prod(x, x));
        Dd::new(x).add_f(r.hi / (2.0 * x))
    }
    fn value(self) -> f64 {
        self.hi + self.lo
    }
    /// self - x for a double x, as a double (exact up to the final rounding).
    fn minus(self, x: f64) -> f64 {
        self.add_f(-x).value()
    }
}

// --------------------------------------------------------------------- lattice

fn cross(u: [f64; 3], v: [f64; 3]) -> [f64; 3] {
    [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ]
}
fn dot(u: [f64; 3], v: [f64; 3]) -> f64 {
    u[0] * v[0] + u[1] * v[1] + u[2] * v[2]
}
fn norm(u: [f64; 3]) -> f64 {
    dot(u, u).sqrt()
}

/// Lattice vectors a_k (rows), reciprocal vectors b_k with a_i . b_k = 2 pi
/// delta_ik, volume in double-double, cell circumradii.
struct Lattice {
    a: [[f64; 3]; 3],
    b: [[f64; 3]; 3],
    omega: Dd,
    rho_r: f64,
    rho_g: f64,
}
fn circumradius(v: &[[f64; 3]; 3]) -> f64 {
    let mut radius: f64 = 0.0;
    for corner in 0..8 {
        let s = [0, 1, 2].map(|k| if corner >> k & 1 == 0 { -0.5 } else { 0.5 });
        let p = [0, 1, 2].map(|c| s[0] * v[0][c] + s[1] * v[1][c] + s[2] * v[2][c]);
        radius = radius.max(norm(p));
    }
    radius
}
impl Lattice {
    fn new(a: [[f64; 3]; 3]) -> Self {
        let c = |u: [f64; 3], v: [f64; 3], i: usize, j: usize| {
            dd_prod(u[i], v[j]).sub(dd_prod(u[j], v[i]))
        };
        let a12 = [
            c(a[1], a[2], 1, 2),
            c(a[1], a[2], 2, 0),
            c(a[1], a[2], 0, 1),
        ];
        let det = a12[0]
            .mul_f(a[0][0])
            .add(a12[1].mul_f(a[0][1]))
            .add(a12[2].mul_f(a[0][2]));
        let omega = if det.hi < 0.0 { det.neg() } else { det };
        let d = det.value();
        let b = [
            cross(a[1], a[2]).map(|x| TAU * x / d),
            cross(a[2], a[0]).map(|x| TAU * x / d),
            cross(a[0], a[1]).map(|x| TAU * x / d),
        ];
        Lattice {
            a,
            b,
            omega,
            rho_r: circumradius(&a),
            rho_g: circumradius(&b),
        }
    }
    fn cart(&self, x: [f64; 3]) -> [f64; 3] {
        [0, 1, 2].map(|c| x[0] * self.a[0][c] + x[1] * self.a[1][c] + x[2] * self.a[2][c])
    }
    fn recip(&self, m: [i32; 3]) -> [f64; 3] {
        let m = m.map(f64::from);
        [0, 1, 2].map(|c| m[0] * self.b[0][c] + m[1] * self.b[1][c] + m[2] * self.b[2][c])
    }
}

/// cos and sin of 2 pi m . s, the argument reduced modulo 1 in double-double.
fn phase(m: [i32; 3], s: &[f64; 3]) -> (f64, f64) {
    let mut x = Dd::ZERO;
    for k in 0..3 {
        x = x.add(dd_prod(f64::from(m[k]), s[k]));
    }
    let n = x.hi.round();
    let f = (x.hi - n) + x.lo;
    let (sin, cos) = (TAU * f).sin_cos();
    (cos, sin)
}

/// Real-space images of pair (i, j) inside `rc`: calls `visit(r, d)` with the
/// Cartesian separation r = R_i - R_j + A m and d = |r| (self image excluded).
fn real_images(
    l: &Lattice,
    si: &[f64; 3],
    sj: &[f64; 3],
    same: bool,
    rc: f64,
    mut visit: impl FnMut([f64; 3], f64),
) {
    let mut d = [si[0] - sj[0], si[1] - sj[1], si[2] - sj[2]];
    for x in &mut d {
        *x -= x.round();
    }
    let h = [0, 1, 2].map(|k| norm(l.b[k]) * rc / TAU);
    let lo = [0, 1, 2].map(|k| (-d[k] - h[k]).floor() as i32 - 1);
    let hi = [0, 1, 2].map(|k| (-d[k] + h[k]).ceil() as i32 + 1);
    for m0 in lo[0]..=hi[0] {
        for m1 in lo[1]..=hi[1] {
            for m2 in lo[2]..=hi[2] {
                if same && [m0, m1, m2] == [0; 3] {
                    continue;
                }
                let r = l.cart([
                    d[0] + f64::from(m0),
                    d[1] + f64::from(m1),
                    d[2] + f64::from(m2),
                ]);
                let d2 = dot(r, r);
                assert!(d2 > 0.0, "coincident ions");
                if d2 <= rc * rc {
                    visit(r, d2.sqrt());
                }
            }
        }
    }
}

/// Reciprocal half space (G and -G represented once) inside `gc`.
fn half_space(l: &Lattice, gc: f64, mut visit: impl FnMut([i32; 3], [f64; 3], f64)) {
    let n = [0, 1, 2].map(|k| (norm(l.a[k]) * gc / TAU).floor() as i32 + 1);
    for m0 in 0..=n[0] {
        for m1 in -n[1]..=n[1] {
            for m2 in -n[2]..=n[2] {
                let positive = m0 > 0 || (m0 == 0 && (m1 > 0 || (m1 == 0 && m2 > 0)));
                if !positive {
                    continue;
                }
                let g = l.recip([m0, m1, m2]);
                let g2 = dot(g, g);
                if g2 <= gc * gc {
                    visit([m0, m1, m2], g, g2);
                }
            }
        }
    }
}

/// Structure factor S(G) = sum q exp(i G . R) in double-double, with the
/// per-ion phases.
fn structure(m: [i32; 3], ions: &[Ion], phases: &mut Vec<(f64, f64)>) -> (Dd, Dd) {
    phases.clear();
    let (mut re, mut im) = (Dd::ZERO, Dd::ZERO);
    for ion in ions {
        let (c, s) = phase(m, &ion.fractional_position);
        phases.push((c, s));
        re = re.add(dd_prod(ion.charge, c));
        im = im.add(dd_prod(ion.charge, s));
    }
    (re, im)
}

#[derive(Clone, Copy, Debug)]
struct RefEwald {
    energy: Dd,
    real: Dd,
    reciprocal: Dd,
    background: Dd,
    /// Analytic truncation bound (Ha) of both sums.
    tail: f64,
    real_terms: usize,
    reciprocal_vectors: usize,
}

/// Lattice-counting tail bound of sum_{|x| > c} f(|x|) over a shifted lattice
/// of cell volume `v` and circumradius `rho`, given f(c) and an upper bound of
/// integral_c^inf f dN_max with N_max(r) = 4 pi (r + rho)^3 / (3 v).
fn tail_bound(f_c: f64, c: f64, rho: f64, v: f64, integral: f64) -> f64 {
    f_c * 4.0 * PI / (3.0 * v) * ((c + rho).powi(3) - (c - rho).max(0.0).powi(3)) + integral
}

/// Ewald energy with splitting parameter `eta` and cutoffs `rc` (Bohr), `gc`
/// (1/Bohr), every term accumulated in double-double.
fn reference_with(rows: [[f64; 3]; 3], ions: &[Ion], eta: f64, rc: f64, gc: f64) -> RefEwald {
    let l = Lattice::new(rows);
    let omega = l.omega.value();
    let mut real = Dd::ZERO;
    let mut real_terms = 0usize;
    for (i, x) in ions.iter().enumerate() {
        for (j, y) in ions.iter().enumerate() {
            let w = 0.5 * x.charge * y.charge;
            real_images(
                &l,
                &x.fractional_position,
                &y.fractional_position,
                i == j,
                rc,
                |_, d| {
                    real = real.add_f(w * libm::erfc(eta * d) / d);
                    real_terms += 1;
                },
            );
        }
    }
    let mut reciprocal = Dd::ZERO;
    let mut reciprocal_vectors = 0usize;
    let mut phases = Vec::new();
    let four_eta2 = 4.0 * eta * eta;
    half_space(&l, gc, |m, _, g2| {
        let (re, im) = structure(m, ions, &mut phases);
        let s2 = re.mul(re).add(im.mul(im));
        // G and -G: 2 (2 pi / Omega) exp(-G^2 / 4 eta^2) / G^2 |S|^2.
        let w = 4.0 * PI / omega * (-g2 / four_eta2).exp() / g2;
        reciprocal = reciprocal.add(s2.mul_f(w));
        reciprocal_vectors += 1;
    });
    let charge = ions.iter().fold(Dd::ZERO, |s, i| s.add_f(i.charge));
    let q2 = ions
        .iter()
        .fold(Dd::ZERO, |s, i| s.add(dd_prod(i.charge, i.charge)));
    let self_term = Dd::new(eta).mul(q2).div(Dd::PI.sqrt()).neg();
    let background = Dd::PI
        .mul(charge.mul(charge))
        .div(l.omega.mul(dd_prod(eta, eta)).mul_f(2.0))
        .neg();
    let energy = real.add(reciprocal).add(self_term).add(background);
    // Tails per (sum |q|)^2 (see the module notes in the ledger row).
    let abs_q: f64 = ions.iter().map(|i| i.charge.abs()).sum();
    let real_tail = 0.5
        * tail_bound(
            libm::erfc(eta * rc) / rc,
            rc,
            l.rho_r,
            omega,
            2.0 * PI * (1.0 + l.rho_r / rc).powi(2) * libm::erfc(eta * rc) / (omega * eta * eta),
        );
    let omega_g = TAU.powi(3) / omega;
    let f_gc = 2.0 * PI / omega * (-gc * gc / four_eta2).exp() / (gc * gc);
    let reciprocal_tail = tail_bound(
        f_gc,
        gc,
        l.rho_g,
        omega_g,
        2.0 * PI / omega * 4.0 * PI / omega_g
            * (1.0 + l.rho_g / gc).powi(2)
            * eta
            * PI.sqrt()
            * libm::erfc(gc / (2.0 * eta)),
    );
    RefEwald {
        energy,
        real,
        reciprocal,
        background,
        tail: abs_q * abs_q * (real_tail + reciprocal_tail),
        real_terms,
        reciprocal_vectors,
    }
}

/// Cutoffs eta r_c = g_c / (2 eta) = `t`.
fn reference(rows: [[f64; 3]; 3], ions: &[Ion], eta: f64, t: f64) -> RefEwald {
    reference_with(rows, ions, eta, t / eta, 2.0 * eta * t)
}

/// The default truncation parameter: exp(-t^2) = 4.8e-25.
const T: f64 = 7.5;

/// sqrt(pi) / Omega^(1/3), the engine's default scale (used only as a unit).
fn eta0(rows: [[f64; 3]; 3]) -> f64 {
    PI.sqrt() / Lattice::new(rows).omega.value().cbrt()
}

/// -dE/dR (Ha/Bohr), double-double accumulation.
fn reference_forces(rows: [[f64; 3]; 3], ions: &[Ion], eta: f64, t: f64) -> Vec<[f64; 3]> {
    let l = Lattice::new(rows);
    let omega = l.omega.value();
    let (rc, gc) = (t / eta, 2.0 * eta * t);
    let mut out = vec![[Dd::ZERO; 3]; ions.len()];
    let c = 2.0 * eta / PI.sqrt();
    for (i, x) in ions.iter().enumerate() {
        for (j, y) in ions.iter().enumerate() {
            let q = x.charge * y.charge;
            real_images(
                &l,
                &x.fractional_position,
                &y.fractional_position,
                i == j,
                rc,
                |r, d| {
                    let f = q
                        * (libm::erfc(eta * d) / (d * d * d)
                            + c * (-(eta * d) * (eta * d)).exp() / (d * d));
                    for k in 0..3 {
                        out[i][k] = out[i][k].add_f(f * r[k]);
                    }
                },
            );
        }
    }
    let mut phases = Vec::new();
    half_space(&l, gc, |m, g, g2| {
        let (re, im) = structure(m, ions, &mut phases);
        let (re, im) = (re.value(), im.value());
        // G and -G: 2 (4 pi / Omega) exp(-G^2 / 4 eta^2) / G^2 Im(S* q_i e^{i G.R_i}) G.
        let w = 8.0 * PI / omega * (-g2 / (4.0 * eta * eta)).exp() / g2;
        for (i, ion) in ions.iter().enumerate() {
            let (cos, sin) = phases[i];
            let imag = ion.charge * (re * sin - im * cos);
            for k in 0..3 {
                out[i][k] = out[i][k].add_f(w * imag * g[k]);
            }
        }
    });
    out.iter().map(|f| f.map(Dd::value)).collect()
}

/// (1/Omega) dE/d(strain) at fixed fractional positions (the engine's sign
/// convention, Ha/Bohr^3), double-double accumulation.
fn reference_stress(rows: [[f64; 3]; 3], ions: &[Ion], eta: f64, t: f64) -> [[f64; 3]; 3] {
    let l = Lattice::new(rows);
    let omega = l.omega.value();
    let (rc, gc) = (t / eta, 2.0 * eta * t);
    let mut out = [[Dd::ZERO; 3]; 3];
    let c = 2.0 * eta / PI.sqrt();
    for (i, x) in ions.iter().enumerate() {
        for (j, y) in ions.iter().enumerate() {
            let q = x.charge * y.charge;
            real_images(
                &l,
                &x.fractional_position,
                &y.fractional_position,
                i == j,
                rc,
                |r, d| {
                    let f = -0.5 * q / omega
                        * (libm::erfc(eta * d) / (d * d * d)
                            + c * (-(eta * d) * (eta * d)).exp() / (d * d));
                    for a in 0..3 {
                        for b in 0..3 {
                            out[a][b] = out[a][b].add_f(f * r[a] * r[b]);
                        }
                    }
                },
            );
        }
    }
    let mut phases = Vec::new();
    half_space(&l, gc, |m, g, g2| {
        let (re, im) = structure(m, ions, &mut phases);
        let s2 = re.mul(re).add(im.mul(im)).value();
        let w = 2.0 * 2.0 * PI / (omega * omega) * (-g2 / (4.0 * eta * eta)).exp() / g2 * s2;
        let k = 2.0 / g2 + 0.5 / (eta * eta);
        for a in 0..3 {
            for b in 0..3 {
                let delta = if a == b { 1.0 } else { 0.0 };
                out[a][b] = out[a][b].add_f(w * (k * g[a] * g[b] - delta));
            }
        }
    });
    let background = reference(rows, ions, eta, t).background.value();
    for (a, row) in out.iter_mut().enumerate() {
        row[a] = row[a].add_f(-background / omega);
    }
    out.map(|r| r.map(Dd::value))
}

// ------------------------------------------------------------------- fixtures

fn ion(charge: f64, s: [f64; 3]) -> Ion {
    Ion {
        charge,
        fractional_position: s,
    }
}
fn cubic(a: f64) -> [[f64; 3]; 3] {
    [[a, 0., 0.], [0., a, 0.], [0., 0., a]]
}
fn scale(rows: [[f64; 3]; 3], s: f64) -> [[f64; 3]; 3] {
    rows.map(|r| r.map(|x| x * s))
}
const FCC: [[f64; 3]; 4] = [[0., 0., 0.], [0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]];
/// NaCl Madelung constant (nearest-neighbour distance), 40 digits.
const MADELUNG_NACL: f64 = 1.747_564_594_633_182_190_636_212_035_544_397;
/// CsCl Madelung constant (nearest-neighbour distance).
const MADELUNG_CSCL: f64 = 1.762_674_773_070_988_397_935_673_320_638_644;

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs()
}

/// Reference energies of `rows`/`ions` at eta = f eta0 for each f, with the
/// relative spread printed.
fn energies(rows: [[f64; 3]; 3], ions: &[Ion], factors: &[f64]) -> Vec<RefEwald> {
    let e0 = eta0(rows);
    factors
        .iter()
        .map(|f| {
            let r = reference(rows, ions, f * e0, T);
            assert!(r.tail < 1e-18, "tail {}", r.tail);
            r
        })
        .collect()
}

#[test]
fn double_double_arithmetic_is_exact_to_about_1e_30() {
    let three = Dd::new(1.0).div(Dd::new(3.0)).mul_f(3.0);
    assert!(three.sub(Dd::new(1.0)).value().abs() < 1e-30);
    let s = Dd::new(2.0).sqrt();
    assert!(s.mul(s).sub(Dd::new(2.0)).value().abs() < 1e-30);
    // 1 + 1e-20 survives in the low word.
    let x = Dd::new(1.0).add_f(1e-20).add_f(-1.0);
    assert!((x.value() - 1e-20).abs() < 1e-35);
    let sp = Dd::PI.sqrt();
    assert!(sp.mul(sp).sub(Dd::PI).value().abs() < 1e-30);
}

/// Analytic checks: NaCl and CsCl Madelung constants (neutral cells, no
/// background) to 1e-13 relative at three eta; the one-component lattice with
/// background: sc 2x2x2 supercell = 8 x primitive, fcc conventional = 4 x
/// primitive, bcc conventional = 2 x primitive, each to 1e-13 relative at three
/// eta (different lattices, different G sets: the background and G = 0
/// convention must be right for these to hold).
#[test]
fn independent_reference_reproduces_madelung_constants_and_supercells() {
    let factors = [0.6, 1.0, 2.0];
    // NaCl: conventional cube a = 2, nearest neighbour 1; 4 formula units.
    let mut nacl = Vec::new();
    for s in FCC {
        nacl.push(ion(1.0, s));
        nacl.push(ion(-1.0, [s[0] + 0.5, s[1], s[2]]));
    }
    for r in energies(cubic(2.0), &nacl, &factors) {
        let want = -4.0 * MADELUNG_NACL;
        println!(
            "NaCl E = {:.17} (hi) {:+.3e} (lo), rel err {:.2e}",
            r.energy.hi,
            r.energy.lo,
            rel(r.energy.value(), want)
        );
        assert!(rel(r.energy.value(), want) < 1e-13);
    }
    // CsCl: cube a = 1, nearest neighbour sqrt(3)/2.
    let cscl = [ion(1.0, [0.; 3]), ion(-1.0, [0.5; 3])];
    for r in energies(cubic(1.0), &cscl, &factors) {
        let want = -MADELUNG_CSCL / (3f64.sqrt() / 2.0);
        println!(
            "CsCl E = {:.17}, rel err {:.2e}",
            r.energy.value(),
            rel(r.energy.value(), want)
        );
        assert!(rel(r.energy.value(), want) < 1e-13);
    }
    // One-component lattices with background (charge 3, a = 7 Bohr).
    let z = 3.0;
    let a = 7.0;
    let sc_super: Vec<Ion> = (0..8)
        .map(|k| ion(z, [0, 1, 2].map(|b| 0.5 * f64::from(k >> b & 1))))
        .collect();
    let fcc_prim = scale([[0., 0.5, 0.5], [0.5, 0., 0.5], [0.5, 0.5, 0.]], a);
    let bcc_prim = scale([[-0.5, 0.5, 0.5], [0.5, -0.5, 0.5], [0.5, 0.5, -0.5]], a);
    /// Name, primitive cell and ions, supercell and ions, ion-count ratio.
    type Supercell = (
        &'static str,
        [[f64; 3]; 3],
        Vec<Ion>,
        [[f64; 3]; 3],
        Vec<Ion>,
        f64,
    );
    let checks: [Supercell; 3] = [
        (
            "sc",
            cubic(a),
            vec![ion(z, [0.; 3])],
            cubic(2.0 * a),
            sc_super,
            8.0,
        ),
        (
            "fcc",
            fcc_prim,
            vec![ion(z, [0.; 3])],
            cubic(a),
            FCC.iter().map(|&s| ion(z, s)).collect(),
            4.0,
        ),
        (
            "bcc",
            bcc_prim,
            vec![ion(z, [0.; 3])],
            cubic(a),
            vec![ion(z, [0.; 3]), ion(z, [0.5; 3])],
            2.0,
        ),
    ];
    // Wigner-Seitz Madelung constants of the one-component lattices (E = -alpha
    // Z^2 / (2 r_s)), literature values quoted to 9 digits.
    let literature = [1.76011888, 1.79174723, 1.79185851];
    for ((name, small_rows, small, big_rows, big, ratio), alpha) in
        checks.into_iter().zip(literature)
    {
        let es = energies(small_rows, &small, &factors);
        let r_s = (3.0 * Lattice::new(small_rows).omega.value() / (4.0 * PI)).cbrt();
        assert!(
            rel(-2.0 * r_s * es[0].energy.value() / (z * z), alpha) < 1e-8,
            "{name}"
        );
        let eb = energies(big_rows, &big, &factors);
        for (s, b) in es.iter().zip(&eb) {
            println!(
                "{name}: E_prim = {:.17}, E_super/{ratio} = {:.17}, rel {:.2e}, alpha_WS = 2 r_s |E|/Z^2 = {:.15}",
                s.energy.value(),
                b.energy.value() / ratio,
                rel(b.energy.value() / ratio, s.energy.value()),
                -2.0 * r_s * s.energy.value() / (z * z)
            );
            assert!(rel(b.energy.value() / ratio, s.energy.value()) < 1e-13);
        }
        for s in &es {
            assert!(rel(s.energy.value(), es[0].energy.value()) < 1e-14);
        }
    }
}

/// Small skewed cell with mixed charges (net charge, background on).
fn skewed() -> ([[f64; 3]; 3], Vec<Ion>) {
    (
        [[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]],
        vec![
            ion(4.0, [0.1, 0.2, 0.3]),
            ion(6.0, [0.55, 0.47, 0.81]),
            ion(1.0, [0.9, 0.73, 0.12]),
            ion(3.5, [0.31, 0.88, 0.6]),
        ],
    )
}

/// Reference forces and stress equal fourth-order central differences of the
/// reference energy (h = 1e-3 Bohr, strain 1e-4): caps 1e-9 Ha/Bohr and
/// 1e-9 Ha (Omega sigma), set before the measurement.
#[test]
fn independent_reference_derivatives_match_finite_differences() {
    let (rows, ions) = skewed();
    let eta = 2.0 * eta0(rows);
    let e = |rows: [[f64; 3]; 3], ions: &[Ion]| reference(rows, ions, eta, T).energy;
    let stencil = |f: &dyn Fn(f64) -> Dd, h: f64| {
        let d1 = f(h).sub(f(-h)).value();
        let d2 = f(2.0 * h).sub(f(-2.0 * h)).value();
        (8.0 * d1 - d2) / (12.0 * h)
    };
    let l = Lattice::new(rows);
    let f_ref = reference_forces(rows, &ions, eta, T);
    let mut worst_f: f64 = 0.0;
    for i in 0..ions.len() {
        for c in 0..3 {
            let shifted = |h: f64| {
                let mut moved = ions.clone();
                for k in 0..3 {
                    moved[i].fractional_position[k] += h * l.b[k][c] / TAU;
                }
                e(rows, &moved)
            };
            let fd = -stencil(&shifted, 1e-3);
            worst_f = worst_f.max((fd - f_ref[i][c]).abs());
        }
    }
    let s_ref = reference_stress(rows, &ions, eta, T);
    let omega = l.omega.value();
    let mut worst_s: f64 = 0.0;
    for a in 0..3 {
        for b in 0..3 {
            let strained = |h: f64| {
                let mut eps = [[0.0; 3]; 3];
                eps[a][b] += 0.5 * h;
                eps[b][a] += 0.5 * h;
                let r = rows.map(|v| {
                    [0, 1, 2].map(|c| v[c] + eps[c][0] * v[0] + eps[c][1] * v[1] + eps[c][2] * v[2])
                });
                e(r, &ions)
            };
            let fd = stencil(&strained, 1e-4);
            worst_s = worst_s.max((fd - omega * s_ref[a][b]).abs());
        }
    }
    println!(
        "reference vs finite differences: force {worst_f:.3e} Ha/Bohr, Omega sigma {worst_s:.3e} Ha"
    );
    assert!(worst_f < 1e-9, "force {worst_f}");
    assert!(worst_s < 1e-9, "stress {worst_s}");
}

/// Both engine schemes against the independent reference on the analytic and
/// skewed cells (debug-cheap); per-atom caps 1e-13 Ha (engine tolerance 1e-13).
#[test]
fn engine_schemes_match_the_independent_reference_on_small_cells() {
    let (rows, ions) = skewed();
    let mut nacl = Vec::new();
    for s in FCC {
        nacl.push(ion(1.0, s));
        nacl.push(ion(-1.0, [s[0] + 0.5, s[1], s[2]]));
    }
    for (name, rows, ions) in [("skewed4", rows, ions), ("nacl8", cubic(2.0), nacl)] {
        let r = reference(rows, &ions, 2.0 * eta0(rows), T).energy;
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        for scheme in [Scheme::Reference, Scheme::Balanced] {
            let e = energy_with(&cell, &ions, EwaldOptions::default(), scheme).unwrap();
            let d = r.minus(e.energy.0);
            println!("{name} {scheme:?}: ref - engine = {d:.3e} Ha");
            assert!(
                d.abs() / ions.len() as f64 <= 1e-13,
                "{name} {scheme:?} {d}"
            );
        }
    }
}

fn si_ions(positions: &[[f64; 3]]) -> Vec<Ion> {
    positions.iter().map(|&s| ion(4.0, s)).collect()
}

/// si32 / si64 at the benchmark geometry (seed 601, amplitude 1e-3, charge 4):
/// the independent reference at eta = (1, 2, 4, 8) x eta0 (t = 7.5), the
/// reference at the engine default scheme's own (eta, r_cut, g_cut), the
/// engine's default-scheme term list summed in double-double, both engine
/// schemes; forces and stress of both schemes against the reference.
/// About 1.9 s in the debug gate profile, 1.3 s in release (Mac). Numbers:
/// `cargo test --release -p dft-physics --lib ewald::tests::independent -- --nocapture --test-threads=1`
#[test]
fn benchmark_cells_independent_reference() {
    let a = 20.522425713435865;
    let cases = [
        (
            "si32",
            [[a, 0., 0.], [0., a, 0.], [0., 0., 10.261212856717933]],
            si_ions(&SI32_FRACTIONAL),
        ),
        ("si64", cubic(a), si_ions(&SI64_FRACTIONAL)),
    ];
    for (name, rows, ions) in cases {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        let options = EwaldOptions::default();
        let default = energy_with(&cell, &ions, options, Scheme::Reference).unwrap();
        let fast = energy_with(&cell, &ions, options, Scheme::Balanced).unwrap();
        let refs = energies(rows, &ions, &[1.0, 2.0, 4.0, 8.0]);
        let best = refs[1].energy;
        println!("{name}: eta0 {:.6}", eta0(rows));
        let mut spread: f64 = 0.0;
        for (f, r) in [1.0, 2.0, 4.0, 8.0].iter().zip(&refs) {
            let d = r.energy.sub(best).value();
            spread = spread.max(d.abs());
            println!(
                "  ref eta {f} eta0: E = {:.16} {:+.3e} | - ref(2 eta0) {d:+.3e} | real terms {} G (half) {} tail {:.1e}",
                r.energy.hi, r.energy.lo, r.real_terms, r.reciprocal_vectors, r.tail
            );
        }
        // The default scheme's own parameters, double-double summation.
        let own = reference_with(
            rows,
            &ions,
            default.eta,
            default.real_cutoff_bohr,
            default.reciprocal_cutoff_inverse_bohr,
        );
        // The default scheme's exact term list, summed in double-double.
        let plan = EwaldPlan::with_scheme(&cell, &ions, options, Scheme::Reference).unwrap();
        let (mut rt, mut gt) = (Vec::new(), Vec::new());
        plan.real_terms(0..plan.real_candidates(), &mut rt).unwrap();
        plan.reciprocal_terms(0..plan.reciprocal_candidates(), &mut gt);
        let dd_real = rt.iter().fold(Dd::ZERO, |s, &t| s.add_f(t));
        let dd_recip = gt.iter().fold(Dd::ZERO, |s, &t| s.add_f(t));
        let dd_default = dd_real
            .add(dd_recip)
            .add_f(default.self_energy.0)
            .add_f(default.background.0);
        // Cost of a compensated (Neumaier) accumulation of the same term list
        // against the plain one (observation, best of 5).
        let time = |f: &dyn Fn() -> f64| {
            (0..5)
                .map(|_| {
                    let t = std::time::Instant::now();
                    std::hint::black_box(f());
                    t.elapsed().as_secs_f64()
                })
                .fold(f64::INFINITY, f64::min)
        };
        let t_plain = time(&|| rt.iter().fold(0.0, |s, &t| s + std::hint::black_box(t)));
        let t_neumaier = time(&|| {
            let (mut s, mut c) = (0.0f64, 0.0f64);
            for &t in &rt {
                let t = std::hint::black_box(t);
                let u = s + t;
                c += if s.abs() >= t.abs() {
                    (s - u) + t
                } else {
                    (t - u) + s
                };
                s = u;
            }
            s + c
        });
        let t_energy = time(&|| {
            energy_with(&cell, &ions, options, Scheme::Reference)
                .unwrap()
                .energy
                .0
        });
        println!(
            "    cost: default energy {t_energy:.4} s; summing its {} real terms plain {t_plain:.5} s, Neumaier {t_neumaier:.5} s",
            rt.len()
        );
        // Mechanism: terms below half an ulp of the final accumulator.
        let half_ulp = default.real.0.abs() * f64::EPSILON / 2.0;
        let small: Vec<f64> = rt.iter().copied().filter(|t| t.abs() < half_ulp).collect();
        println!(
            "    terms below half an ulp of the real sum ({half_ulp:.2e}): {} of {}, their sum {:.3e} Ha",
            small.len(),
            rt.len(),
            small.iter().fold(Dd::ZERO, |s, &t| s.add_f(t)).value()
        );
        println!(
            "    ref(own params) real {:.16} recip {:.16} | default real {:.16} recip {:.16}",
            own.real.value(),
            own.reciprocal.value(),
            default.real.0,
            default.reciprocal.0
        );
        println!(
            "  default  E = {:.16} | ref - default {:+.3e} (eta {:.4} r_cut {:.2} g_cut {:.3} tail {:.1e})",
            default.energy.0,
            best.minus(default.energy.0),
            default.eta,
            default.real_cutoff_bohr,
            default.reciprocal_cutoff_inverse_bohr,
            default.tail_bound.0
        );
        println!(
            "    default real sum naive {:.16} dd {:.16} {:+.3e} (naive - dd {:+.3e}, {} terms); recip naive - dd {:+.3e}; ref(own params) - default {:+.3e}; ref - dd-summed default {:+.3e}",
            default.real.0,
            dd_real.hi,
            dd_real.lo,
            -dd_real.minus(default.real.0),
            rt.len(),
            -dd_recip.minus(default.reciprocal.0),
            own.energy.minus(default.energy.0),
            best.sub(dd_default).value()
        );
        println!(
            "  fast     E = {:.16} | ref - fast {:+.3e} (eta {:.4} r_cut {:.2} g_cut {:.3} tail {:.1e})",
            fast.energy.0,
            best.minus(fast.energy.0),
            fast.eta,
            fast.real_cutoff_bohr,
            fast.reciprocal_cutoff_inverse_bohr,
            fast.tail_bound.0
        );
        println!(
            "  spread of the reference over eta: {spread:.3e} Ha; fast - default {:+.3e}",
            fast.energy.0 - default.energy.0
        );
        // Pre-registered caps: reference spread 1e-12 Ha (64 ions, |E| ~ 500 Ha:
        // a few ulp); fast within 1e-11 Ha of the reference (its truncation
        // budget 6.4e-12 plus rounding). The default is only bounded loosely
        // (standard class 5e-9 Ha); its deviation is the measured finding.
        assert!(spread < 1e-12, "spread {spread}");
        assert!(best.minus(fast.energy.0).abs() < 1e-11);
        assert!(best.minus(default.energy.0).abs() < 5e-9);

        // Forces and stress against the reference at eta = 2 eta0.
        let eta = 2.0 * eta0(rows);
        let f_ref = reference_forces(rows, &ions, eta, T);
        let s_ref = reference_stress(rows, &ions, eta, T);
        let f_ref4 = reference_forces(rows, &ions, 2.0 * eta, T);
        let s_ref4 = reference_stress(rows, &ions, 2.0 * eta, T);
        let max_f = |f: &[[f64; 3]], g: &[[f64; 3]]| {
            f.iter()
                .flatten()
                .zip(g.iter().flatten())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max)
        };
        let max_s = |s: &[[f64; 3]; 3], t: &[[f64; 3]; 3]| {
            s.iter()
                .flatten()
                .zip(t.iter().flatten())
                .map(|(a, b)| (a - b).abs())
                .fold(0.0, f64::max)
        };
        println!(
            "  reference self-consistency (eta 2 vs 4 eta0): force {:.3e} Ha/Bohr, stress {:.3e} Ha/Bohr^3",
            max_f(&f_ref, &f_ref4),
            max_s(&s_ref, &s_ref4)
        );
        for scheme in [Scheme::Reference, Scheme::Balanced] {
            set_scheme_override(Some(scheme));
            let f = forces(&cell, &ions, options).unwrap();
            let s = stress(&cell, &ions, options).unwrap();
            set_scheme_override(None);
            let (df, ds) = (max_f(&f, &f_ref), max_s(&s, &s_ref));
            println!("  {scheme:?}: max |dF| {df:.3e} Ha/Bohr, max |dsigma| {ds:.3e} Ha/Bohr^3");
            assert!(df < 1e-10 && ds < 1e-12, "{scheme:?} {df} {ds}");
        }
    }
}
