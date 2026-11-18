//! Point-charge Ewald sum with a uniform neutralizing background.
//!
//! Two schemes (INIT-EWALD-FAST, `BAANDI_EWALD_FAST`):
//! - [`Scheme::Reference`] (default): eta = sqrt(pi)/Omega^(1/3), the historical
//!   truncation spheres and summation order. The real-space images of every ion
//!   pair are enumerated only near the cutoff sphere (conservative pruning of
//!   the old enumeration box); every candidate the old loop accepted is still
//!   tested with the same predicate and added in the same order, so the result
//!   is bit-identical to the previous evaluator.
//! - [`Scheme::Balanced`]: eta from the real/reciprocal cost balance for the
//!   cell volume and ion count, cutoffs from rigorous lattice-counting tail
//!   bounds, ion pairs i <= j and the reciprocal half space with per-ion phase
//!   tables. Rounding class (different truncation and summation).
use crate::PhysicsError;
use dft_geometry::Cell;
use dft_units::Hartree;
use nalgebra::{Matrix3, Vector3};
use num_complex::Complex64;
use std::f64::consts::{PI, TAU};

#[derive(Clone, Copy, Debug)]
pub struct Ion {
    pub charge: f64,
    pub fractional_position: [f64; 3],
}

#[derive(Clone, Copy, Debug)]
pub struct EwaldOptions {
    /// Inverse Bohr. None chooses sqrt(pi)/Omega^(1/3) (reference scheme) or
    /// the cost-balanced value (balanced scheme).
    pub eta: Option<f64>,
    /// Absolute truncation budget = tolerance * (sum |Z|)^2/Omega^(1/3).
    /// This is a scale-relative bound, not relative to a possibly zero energy.
    /// The balanced scheme also caps it at tolerance * N Ha.
    pub tolerance: f64,
}
impl Default for EwaldOptions {
    fn default() -> Self {
        Self {
            eta: None,
            tolerance: 1e-13,
        }
    }
}

#[derive(Clone, Debug)]
pub struct EwaldEnergy {
    pub energy: Hartree,
    pub real: Hartree,
    pub reciprocal: Hartree,
    pub self_energy: Hartree,
    pub background: Hartree,
    pub eta: f64,
    pub real_cutoff_bohr: f64,
    pub reciprocal_cutoff_inverse_bohr: f64,
    /// Conservative analytic truncation bound; excludes floating-point error.
    pub tail_bound: Hartree,
}

/// Evaluation scheme of [`energy`], [`forces`], [`stress`] and [`EwaldPlan`].
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Scheme {
    /// Historical parameters and summation, bit-identical to the evaluator
    /// before INIT-EWALD-FAST.
    Reference,
    /// Cost-balanced splitting parameter and tight truncation (rounding class).
    Balanced,
}

thread_local! {
    static SCHEME_OVERRIDE: std::cell::Cell<Option<Scheme>> = const { std::cell::Cell::new(None) };
}
/// Override `BAANDI_EWALD_FAST` on the calling thread (tests); `None` restores
/// the environment.
pub fn set_scheme_override(scheme: Option<Scheme>) {
    SCHEME_OVERRIDE.with(|o| o.set(scheme));
}
/// The scheme in effect on the calling thread.
pub fn scheme() -> Scheme {
    SCHEME_OVERRIDE.with(|o| o.get()).unwrap_or_else(|| {
        if dft_numerics::switches::on(dft_numerics::switches::Id::EWALD_FAST) {
            Scheme::Balanced
        } else {
            Scheme::Reference
        }
    })
}

/// Cost of one real-space image term (erfc, square root, division, image
/// enumeration) relative to one ion term of a reciprocal vector (two complex
/// products) in the balanced scheme. A property of the implementation and the
/// machine, not of any material: scanned with `examples/ewald_bench` on
/// cubic, skewed, elongated and slab cells with 1 to 216 ions (Apple M-series,
/// release, 2026-11-18): the time is flat within about 10 % for ratios 12 to
/// 64 on every cell, ratio 1 is about 2x slower. Any value changes only the
/// speed; the truncation bound holds for every eta.
pub const BALANCE_COST_RATIO: f64 = 16.0;

fn matrix(a: [[f64; 3]; 3]) -> Matrix3<f64> {
    Matrix3::from_fn(|i, j| a[i][j])
}

/// Tensile-positive ion/background strain derivative in Ha/Bohr^3, at fixed
/// fractional ion positions. Eta is held fixed in the derivative: its response
/// cancels between the four Ewald terms in the infinite sums. The energy
/// evaluator's truncation spheres are reused; its energy tail bound is not a
/// stress error bound. This is an isolated contribution, not a total stress.
pub fn stress(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
) -> Result<[[f64; 3]; 3], PhysicsError> {
    let parameters = energy(cell, ions, options)?;
    let p = prepare(cell, ions, options, scheme())?;
    let eta = parameters.eta;
    let omega = cell.volume_bohr3();
    let a = p.a;
    let b = p.b;
    let gb = bounds(&a, parameters.reciprocal_cutoff_inverse_bohr, 0.0)?;
    let positions = &p.positions;
    let mut result = [[0.0; 3]; 3];
    for (i, ion) in ions.iter().enumerate() {
        for (j, other) in ions.iter().enumerate() {
            let delta = positions[i] - positions[j];
            p.walker.images(&delta, |m| {
                if i == j && m == [0; 3] {
                    return Ok(());
                }
                let r = a * (positions[i] - positions[j] + Vector3::from(m.map(f64::from)));
                let distance = r.norm();
                if distance <= parameters.real_cutoff_bohr {
                    let factor = -0.5 * ion.charge * other.charge / omega
                        * (libm::erfc(eta * distance) / distance.powi(3)
                            + 2.0 * eta / PI.sqrt() * (-(eta * distance).powi(2)).exp()
                                / distance.powi(2));
                    for (a, row) in result.iter_mut().enumerate() {
                        for (b, value) in row.iter_mut().enumerate() {
                            *value += factor * r[a] * r[b];
                        }
                    }
                }
                Ok(())
            })?;
        }
    }
    enumerate(
        |m| {
            if m == [0; 3] {
                return Ok(());
            }
            let mv = Vector3::from(m.map(f64::from));
            let g = b * mv;
            let g2 = g.norm_squared();
            if g2 <= parameters.reciprocal_cutoff_inverse_bohr.powi(2) {
                let structure: Complex64 = ions
                    .iter()
                    .zip(positions)
                    .map(|(ion, s)| Complex64::from_polar(ion.charge, TAU * mv.dot(s)))
                    .sum();
                let weight = 2.0 * PI / omega.powi(2) * (-g2 / (4.0 * eta * eta)).exp() / g2
                    * structure.norm_sqr();
                for (a, row) in result.iter_mut().enumerate() {
                    for (b, value) in row.iter_mut().enumerate() {
                        *value += weight
                            * ((2.0 / g2 + 0.5 / eta.powi(2)) * g[a] * g[b]
                                - if a == b { 1.0 } else { 0.0 });
                    }
                }
            }
            Ok(())
        },
        gb,
    )?;
    // The self term has no explicit fixed-eta strain derivative. Background
    // scales as 1/Omega and must not be mistaken for a zero-mode constant.
    for (a, row) in result.iter_mut().enumerate() {
        row[a] -= parameters.background.0 / omega;
    }
    if result.iter().flatten().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(result)
}
fn covering_radius(a: &Matrix3<f64>) -> f64 {
    (0..3).map(|i| a.column(i).norm()).sum::<f64>() / 2.0
}
/// Circumradius of the parallelepiped spanned by the columns of `a`, centred
/// at the origin: every lattice cell (translated to its lattice point) lies in
/// a ball of this radius around that point.
fn circumradius(a: &Matrix3<f64>) -> f64 {
    let mut radius: f64 = 0.0;
    for corner in 0..8 {
        let s = Vector3::new(
            if corner & 1 == 0 { -0.5 } else { 0.5 },
            if corner & 2 == 0 { -0.5 } else { 0.5 },
            if corner & 4 == 0 { -0.5 } else { 0.5 },
        );
        radius = radius.max((a * s).norm());
    }
    radius
}
fn bounds(dual: &Matrix3<f64>, radius: f64, shift: f64) -> Result<[i32; 3], PhysicsError> {
    let mut bounds = [0; 3];
    for (i, b) in bounds.iter_mut().enumerate() {
        let value = (dual.column(i).norm() * radius / TAU + shift).ceil() + 1.0;
        if !value.is_finite() || value >= i32::MAX as f64 {
            return Err(PhysicsError::Unrepresentable);
        }
        *b = value as i32;
    }
    Ok(bounds)
}
fn enumerate(
    mut visit: impl FnMut([i32; 3]) -> Result<(), PhysicsError>,
    b: [i32; 3],
) -> Result<(), PhysicsError> {
    for i in -b[0]..=b[0] {
        for j in -b[1]..=b[1] {
            for k in -b[2]..=b[2] {
                visit([i, j, k])?;
            }
        }
    }
    Ok(())
}

/// Conservative enumeration of the lattice translations m of the box `rb`
/// whose image A (delta + m) can lie within `radius`, in the order of
/// [`enumerate`] (INIT-EWALD-FAST). Necessary conditions, each widened by one
/// index and a relative radius slack of 1e-9: |delta_k + m_k| <= |b_k| r / 2 pi
/// for k = 0, 1 (b_k . A x = 2 pi x_k), and for every (m0, m1) row the interval
/// of m2 where the quadratic |p + a_2 m2|^2 <= r^2 holds. A superset of the
/// sphere's points inside the box, so a caller that applies its exact
/// predicate to every visited point accepts exactly the points the full box
/// loop accepted, in the same order.
struct ImageWalker {
    a: Matrix3<f64>,
    a2: Vector3<f64>,
    a2_norm2: f64,
    half: [f64; 3],
    radius2: f64,
    rb: [i32; 3],
}
impl ImageWalker {
    fn new(a: Matrix3<f64>, b: &Matrix3<f64>, radius: f64, rb: [i32; 3]) -> Self {
        let slack = radius * (1.0 + 1e-9);
        let a2: Vector3<f64> = a.column(2).into();
        Self {
            a,
            a2,
            a2_norm2: a2.norm_squared(),
            half: std::array::from_fn(|k| b.column(k).norm() * slack / TAU),
            radius2: slack * slack,
            rb,
        }
    }
    /// Integers m in [-bound, bound] that can satisfy |center + m| <= half.
    fn range(center: f64, half: f64, bound: i32) -> std::ops::RangeInclusive<i32> {
        let lo = (-center - half).floor() - 1.0;
        let hi = (-center + half).ceil() + 1.0;
        if !(lo.is_finite() && hi.is_finite()) {
            return -bound..=bound;
        }
        let lo = lo.max(-f64::from(bound)) as i32;
        let hi = hi.min(f64::from(bound)) as i32;
        lo..=hi
    }
    fn m0_range(&self, delta: &Vector3<f64>) -> std::ops::RangeInclusive<i32> {
        Self::range(delta[0], self.half[0], self.rb[0])
    }
    /// Visits the candidates of the row m0 (which must lie in the box).
    #[inline]
    fn row(
        &self,
        delta: &Vector3<f64>,
        m0: i32,
        mut visit: impl FnMut([i32; 3]) -> Result<(), PhysicsError>,
    ) -> Result<(), PhysicsError> {
        for m1 in Self::range(delta[1], self.half[1], self.rb[1]) {
            let p =
                self.a * Vector3::new(delta[0] + f64::from(m0), delta[1] + f64::from(m1), delta[2]);
            let pa = p.dot(&self.a2);
            let disc = pa * pa - self.a2_norm2 * (p.norm_squared() - self.radius2);
            let span = if disc.is_nan() {
                -self.rb[2]..=self.rb[2]
            } else if disc < 0.0 {
                continue;
            } else {
                let s = disc.sqrt();
                let lo = ((-pa - s) / self.a2_norm2).floor() - 1.0;
                let hi = ((-pa + s) / self.a2_norm2).ceil() + 1.0;
                if !(lo.is_finite() && hi.is_finite()) {
                    -self.rb[2]..=self.rb[2]
                } else {
                    let lo = lo.max(-f64::from(self.rb[2])) as i32;
                    let hi = hi.min(f64::from(self.rb[2])) as i32;
                    lo..=hi
                }
            };
            for m2 in span {
                visit([m0, m1, m2])?;
            }
        }
        Ok(())
    }
    fn images(
        &self,
        delta: &Vector3<f64>,
        mut visit: impl FnMut([i32; 3]) -> Result<(), PhysicsError>,
    ) -> Result<(), PhysicsError> {
        for m0 in self.m0_range(delta) {
            self.row(delta, m0, &mut visit)?;
        }
        Ok(())
    }
}

/// Parameters shared by the serial evaluation and the distributed plan.
struct Prepared {
    scheme: Scheme,
    eta: f64,
    a: Matrix3<f64>,
    b: Matrix3<f64>,
    omega: f64,
    r_cut: f64,
    g_cut: f64,
    tail: f64,
    rb: [i32; 3],
    gb: [i32; 3],
    positions: Vec<Vector3<f64>>,
    walker: ImageWalker,
    /// Balanced only: ion pairs i <= j in order.
    pairs: Vec<(usize, usize)>,
    /// Balanced only: reciprocal half space inside the cutoff (Miller, |G|^2).
    half_space: Vec<([i32; 3], f64)>,
    /// Balanced only: phase tables per axis, entry ((m + gb_k) * n + ion);
    /// axis 0 carries the ion charge.
    tables: [Vec<Complex64>; 3],
}

fn validate(eta: f64, options: &EwaldOptions, ions: &[Ion]) -> Result<(), PhysicsError> {
    if !eta.is_finite()
        || eta <= 0.0
        || !options.tolerance.is_finite()
        || options.tolerance <= 0.0
        || options.tolerance >= 1.0
        || ions
            .iter()
            .any(|i| !i.charge.is_finite() || i.fractional_position.iter().any(|x| !x.is_finite()))
    {
        return Err(PhysicsError::EwaldInput);
    }
    Ok(())
}

fn prepare(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
    scheme: Scheme,
) -> Result<Prepared, PhysicsError> {
    match scheme {
        Scheme::Reference => prepare_reference(cell, ions, options),
        Scheme::Balanced => prepare_balanced(cell, ions, options),
    }
}

fn positions_of(ions: &[Ion]) -> Vec<Vector3<f64>> {
    ions.iter()
        .map(|i| Vector3::from(i.fractional_position.map(|s| s.rem_euclid(1.0))))
        .collect()
}

fn prepare_reference(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
) -> Result<Prepared, PhysicsError> {
    let length = cell.volume_bohr3().cbrt();
    let eta = options.eta.unwrap_or(PI.sqrt() / length);
    validate(eta, &options, ions)?;
    let a = matrix(cell.lattice_columns_bohr());
    let b = matrix(cell.reciprocal_columns_inverse_bohr());
    let rho_r = covering_radius(&a);
    let rho_g = covering_radius(&b);
    let omega = cell.volume_bohr3();
    let omega_g = TAU.powi(3) / omega;
    // Bounds are per (sum |Z|)^2, avoiding division by zero for empty/zero charges.
    let budget = options.tolerance / length;
    let mut t = (-options.tolerance.ln()).sqrt();
    let (r_cut, g_cut, tail) = loop {
        let u = t / eta;
        let v = 2.0 * eta * t;
        let real_bound =
            PI.sqrt() / (omega * eta.powi(3) * u) * (1.0 + rho_r / u).powi(2) * (-t * t).exp();
        let reciprocal_bound = (2.0 * PI / omega)
            * (8.0 * PI * eta * eta / (omega_g * v))
            * (1.0 + rho_g / v).powi(2)
            * (-t * t).exp();
        if !u.is_finite()
            || !v.is_finite()
            || !real_bound.is_finite()
            || !reciprocal_bound.is_finite()
            || budget == 0.0
        {
            return Err(PhysicsError::Unrepresentable);
        }
        if real_bound + reciprocal_bound <= budget {
            break (
                u + 2.0 * rho_r,
                v + 2.0 * rho_g,
                real_bound + reciprocal_bound,
            );
        }
        t += 1.0;
    };
    let rb = bounds(&b, r_cut, 1.0)?;
    let gb = bounds(&a, g_cut, 0.0)?;
    Ok(Prepared {
        scheme: Scheme::Reference,
        eta,
        a,
        b,
        omega,
        r_cut,
        g_cut,
        tail,
        rb,
        gb,
        positions: positions_of(ions),
        walker: ImageWalker::new(a, &b, r_cut, rb),
        pairs: Vec::new(),
        half_space: Vec::new(),
        tables: [Vec::new(), Vec::new(), Vec::new()],
    })
}

/// Truncation bound of sum_{|x| > cut} f(|x|) over a (shifted) lattice with
/// cell volume `volume` and cell circumradius `rho`, for f positive and
/// decreasing: the counting function lies between the volumes of the balls of
/// radius r -+ rho divided by `volume`, so the sum is at most
/// f(cut) [N_max(cut) - N_min(cut)] + integral_cut^inf f dN_max; `integral`
/// is an upper bound of the latter.
fn lattice_tail(f_cut: f64, cut: f64, rho: f64, volume: f64, integral: f64) -> f64 {
    let shell = 4.0 * PI / (3.0 * volume) * ((cut + rho).powi(3) - (cut - rho).max(0.0).powi(3));
    f_cut * shell + integral
}

/// Balanced truncation bounds (per (sum |Z|)^2) at t = eta r_c = g_c / (2 eta).
fn balanced_bounds(t: f64, eta: f64, omega: f64, rho_r: f64, rho_g: f64) -> (f64, f64, f64) {
    let rc = t / eta;
    let gc = 2.0 * eta * t;
    let e = (-t * t).exp();
    // Real space: 1/2 sum_{|x| > rc} erfc(eta x)/x; with erfc(y) <= exp(-y^2)/(y sqrt(pi))
    // and (r + rho)^2 / r <= r (1 + rho/rc)^2 the integral is at most
    // 2 sqrt(pi) (1 + rho/rc)^2 exp(-t^2) / (Omega eta^3 rc).
    let real_integral =
        2.0 * PI.sqrt() * (1.0 + rho_r / rc).powi(2) * e / (omega * eta.powi(3) * rc);
    let real = 0.5 * lattice_tail(libm::erfc(t) / rc, rc, rho_r, omega, real_integral);
    // Reciprocal: sum_{|G| > gc} (2 pi/Omega) exp(-G^2/4eta^2)/G^2 (|S|^2 <= (sum |Z|)^2);
    // integral at most (2 pi/Omega)(4 pi/Omega_g)(1 + rho/gc)^2 2 eta^2 exp(-t^2)/gc.
    let omega_g = TAU.powi(3) / omega;
    let f_gc = 2.0 * PI / omega * e / (gc * gc);
    let reciprocal_integral =
        2.0 * PI / omega * 4.0 * PI / omega_g * (1.0 + rho_g / gc).powi(2) * 2.0 * eta * eta * e
            / gc;
    let reciprocal = lattice_tail(f_gc, gc, rho_g, omega_g, reciprocal_integral);
    (rc, gc, real + reciprocal)
}

fn prepare_balanced(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
) -> Result<Prepared, PhysicsError> {
    let omega = cell.volume_bohr3();
    let length = omega.cbrt();
    let n = ions.len();
    // Real-space work ~ N^2/2 (4 pi/3) rc^3/Omega, reciprocal ~ N (1/2)(4 pi/3) gc^3 Omega/(2 pi)^3
    // with rc = t/eta, gc = 2 eta t: equal for eta^6 = pi^3 N c / Omega^2.
    let eta = options
        .eta
        .unwrap_or(PI.sqrt() * ((n.max(1) as f64) * BALANCE_COST_RATIO).powf(1.0 / 6.0) / length);
    validate(eta, &options, ions)?;
    let a = matrix(cell.lattice_columns_bohr());
    let b = matrix(cell.reciprocal_columns_inverse_bohr());
    let rho_r = circumradius(&a);
    let rho_g = circumradius(&b);
    let abs_charge: f64 = ions.iter().map(|i| i.charge.abs()).sum();
    let mut budget = options.tolerance / length;
    if n > 0 && abs_charge > 0.0 {
        budget = budget.min(options.tolerance * n as f64 / (abs_charge * abs_charge));
    }
    // Positive and finite (a NaN budget is rejected too).
    if budget.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) || !budget.is_finite() {
        return Err(PhysicsError::Unrepresentable);
    }
    let mut t: f64 = 1.0;
    let (r_cut, g_cut, tail) = loop {
        let (rc, gc, bound) = balanced_bounds(t, eta, omega, rho_r, rho_g);
        if !rc.is_finite() || !gc.is_finite() || !bound.is_finite() || t > 64.0 {
            return Err(PhysicsError::Unrepresentable);
        }
        if bound <= budget {
            break (rc, gc, bound);
        }
        t += 0.01;
    };
    let rb = bounds(&b, r_cut, 1.0)?;
    let gb = bounds(&a, g_cut, 0.0)?;
    let positions = positions_of(ions);
    let pairs = (0..n).flat_map(|i| (i..n).map(move |j| (i, j))).collect();
    let mut half_space = Vec::new();
    enumerate(
        |m| {
            let positive = m[0] > 0 || (m[0] == 0 && (m[1] > 0 || (m[1] == 0 && m[2] > 0)));
            if positive {
                let g2 = (b * Vector3::from(m.map(f64::from))).norm_squared();
                if g2 <= g_cut * g_cut {
                    half_space.push((m, g2));
                }
            }
            Ok(())
        },
        gb,
    )?;
    let tables = std::array::from_fn(|k| {
        let mut table = Vec::with_capacity((2 * gb[k] as usize + 1) * n);
        for m in -gb[k]..=gb[k] {
            for (ion, s) in ions.iter().zip(&positions) {
                let r = if k == 0 { ion.charge } else { 1.0 };
                table.push(Complex64::from_polar(r, TAU * f64::from(m) * s[k]));
            }
        }
        table
    });
    Ok(Prepared {
        scheme: Scheme::Balanced,
        eta,
        a,
        b,
        omega,
        r_cut,
        g_cut,
        tail,
        rb,
        gb,
        positions,
        walker: ImageWalker::new(a, &b, r_cut, rb),
        pairs,
        half_space,
        tables,
    })
}

/// Reference real-space terms of ion pair (i, j) and box row m0 in serial order.
fn reference_real_row(
    p: &Prepared,
    ions: &[Ion],
    i: usize,
    j: usize,
    m0: i32,
    visit: &mut impl FnMut(f64),
) -> Result<(), PhysicsError> {
    let delta = p.positions[i] - p.positions[j];
    if !p.walker.m0_range(&delta).contains(&m0) {
        return Ok(());
    }
    let (qi, qj) = (ions[i].charge, ions[j].charge);
    p.walker.row(&delta, m0, |m| {
        if i == j && m == [0; 3] {
            return Ok(());
        }
        let distance = (p.a * (delta + Vector3::from(m.map(f64::from)))).norm();
        if distance == 0.0 {
            return Err(PhysicsError::CoincidentIons);
        }
        if distance <= p.r_cut {
            visit(0.5 * qi * qj * libm::erfc(p.eta * distance) / distance);
        }
        Ok(())
    })
}

/// Balanced real-space term of the ion pair (i <= j): w q_i q_j sum erfc/x.
fn balanced_real_pair(p: &Prepared, ions: &[Ion], i: usize, j: usize) -> Result<f64, PhysicsError> {
    let delta = p.positions[i] - p.positions[j];
    let mut sum = 0.0;
    p.walker.images(&delta, |m| {
        if i == j && m == [0; 3] {
            return Ok(());
        }
        let r = p.a * (delta + Vector3::from(m.map(f64::from)));
        let d2 = r.norm_squared();
        if d2 == 0.0 {
            return Err(PhysicsError::CoincidentIons);
        }
        if d2 <= p.r_cut * p.r_cut {
            let d = d2.sqrt();
            sum += libm::erfc(p.eta * d) / d;
        }
        Ok(())
    })?;
    let weight = if i == j { 0.5 } else { 1.0 };
    Ok(weight * ions[i].charge * ions[j].charge * sum)
}

/// Reciprocal term of reciprocal lattice point m; None for m = 0 or outside the
/// cutoff sphere (reference scheme).
fn reciprocal_term(p: &Prepared, ions: &[Ion], m: [i32; 3]) -> Option<f64> {
    if m == [0; 3] {
        return None;
    }
    let mv = Vector3::from(m.map(f64::from));
    let g2 = (p.b * mv).norm_squared();
    // The serial predicate exactly (a NaN g2 is skipped, as before).
    let inside = g2 <= p.g_cut * p.g_cut;
    if !inside {
        return None;
    }
    let structure: Complex64 = ions
        .iter()
        .zip(&p.positions)
        .map(|(ion, s)| Complex64::from_polar(ion.charge, TAU * mv.dot(s)))
        .sum();
    Some((2.0 * PI / p.omega) * (-g2 / (4.0 * p.eta * p.eta)).exp() / g2 * structure.norm_sqr())
}

/// Balanced reciprocal term of half-space vector `index` (G and -G together).
fn balanced_reciprocal(p: &Prepared, index: usize) -> f64 {
    let (m, g2) = p.half_space[index];
    let n = p.positions.len();
    let row = |k: usize| {
        let at = (m[k] + p.gb[k]) as usize * n;
        &p.tables[k][at..at + n]
    };
    let (t0, t1, t2) = (row(0), row(1), row(2));
    let mut s = Complex64::new(0.0, 0.0);
    for ((x, y), z) in t0.iter().zip(t1).zip(t2) {
        s += x * y * z;
    }
    (4.0 * PI / p.omega) * (-g2 / (4.0 * p.eta * p.eta)).exp() / g2 * s.norm_sqr()
}

fn finish(
    p: &Prepared,
    ions: &[Ion],
    real: f64,
    reciprocal: f64,
) -> Result<EwaldEnergy, PhysicsError> {
    let eta = p.eta;
    let omega = p.omega;
    let total_charge: f64 = ions.iter().map(|i| i.charge).sum();
    let abs_charge: f64 = ions.iter().map(|i| i.charge.abs()).sum();
    let self_energy = -eta / PI.sqrt() * ions.iter().map(|i| i.charge * i.charge).sum::<f64>();
    let background = -PI * total_charge * total_charge / (2.0 * omega * eta * eta);
    let energy = real + reciprocal + self_energy + background;
    if !energy.is_finite() || !(p.tail * abs_charge * abs_charge).is_finite() {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(EwaldEnergy {
        energy: Hartree(energy),
        real: Hartree(real),
        reciprocal: Hartree(reciprocal),
        self_energy: Hartree(self_energy),
        background: Hartree(background),
        eta,
        real_cutoff_bohr: p.r_cut,
        reciprocal_cutoff_inverse_bohr: p.g_cut,
        tail_bound: Hartree(p.tail * abs_charge * abs_charge),
    })
}

/// Ewald energy with the scheme of [`scheme`] (reference unless
/// `BAANDI_EWALD_FAST=on`).
pub fn energy(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
) -> Result<EwaldEnergy, PhysicsError> {
    energy_with(cell, ions, options, scheme())
}

/// Ewald energy with an explicit scheme.
pub fn energy_with(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
    scheme: Scheme,
) -> Result<EwaldEnergy, PhysicsError> {
    let plan = EwaldPlan::with_scheme(cell, ions, options, scheme)?;
    let mut real = 0.0;
    plan.visit_real(0..plan.real_candidates, &mut |t| real += t)?;
    let mut reciprocal = 0.0;
    plan.visit_reciprocal(0..plan.reciprocal_candidates, &mut |t| reciprocal += t);
    plan.finish_sums(real, reciprocal)
}

/// Points of the enumeration box `b` (the order of `enumerate`).
fn box_points(b: [i32; 3]) -> Result<usize, PhysicsError> {
    b.iter()
        .map(|&e| usize::try_from(2 * i64::from(e) + 1).ok())
        .try_fold(1usize, |acc, n| n.and_then(|n| acc.checked_mul(n)))
        .ok_or(PhysicsError::Unrepresentable)
}

/// Point `k` of the enumeration box `b`, in the order of `enumerate`.
fn box_point(b: [i32; 3], k: usize) -> [i32; 3] {
    let n1 = (2 * b[1] + 1) as usize;
    let n2 = (2 * b[2] + 1) as usize;
    [
        (k / (n1 * n2)) as i32 - b[0],
        ((k / n2) % n1) as i32 - b[1],
        (k % n2) as i32 - b[2],
    ]
}

/// The Ewald sum split into independent terms (INIT-SCALE-3): the candidates of
/// the real-space sum and of the reciprocal sum are numbered; `real_terms` /
/// `reciprocal_terms` return the accepted terms of a candidate range in serial
/// order; `finish` adds a complete term list sequentially from +0.0. Skipped
/// candidates add nothing in the serial loop, so the terms of contiguous
/// candidate ranges, concatenated in range order and summed by `finish`, give
/// [`energy`] bit for bit for any partition.
///
/// Reference scheme: a real-space candidate is an ion pair (i, j in order) and
/// a row m0 of the enumeration box (INIT-EWALD-FAST; the row's translations are
/// pruned to the cutoff sphere, see `ImageWalker`); a reciprocal candidate is a
/// point of the reciprocal box. Balanced scheme: a real-space candidate is an
/// ion pair i <= j (one term per pair), a reciprocal candidate a vector of the
/// half space inside the cutoff.
pub struct EwaldPlan<'a> {
    ions: &'a [Ion],
    p: Prepared,
    rows: usize,
    real_candidates: usize,
    reciprocal_candidates: usize,
}
impl<'a> EwaldPlan<'a> {
    /// Plan with the scheme of [`scheme`].
    pub fn new(cell: &Cell, ions: &'a [Ion], options: EwaldOptions) -> Result<Self, PhysicsError> {
        Self::with_scheme(cell, ions, options, scheme())
    }
    pub fn with_scheme(
        cell: &Cell,
        ions: &'a [Ion],
        options: EwaldOptions,
        scheme: Scheme,
    ) -> Result<Self, PhysicsError> {
        let p = prepare(cell, ions, options, scheme)?;
        let (rows, real_candidates, reciprocal_candidates) = match scheme {
            Scheme::Reference => {
                let rows = usize::try_from(2 * i64::from(p.rb[0]) + 1)
                    .map_err(|_| PhysicsError::Unrepresentable)?;
                let real = ions
                    .len()
                    .checked_mul(ions.len())
                    .and_then(|n| n.checked_mul(rows))
                    .ok_or(PhysicsError::Unrepresentable)?;
                (rows, real, box_points(p.gb)?)
            }
            Scheme::Balanced => (0, p.pairs.len(), p.half_space.len()),
        };
        Ok(Self {
            ions,
            p,
            rows,
            real_candidates,
            reciprocal_candidates,
        })
    }
    pub fn scheme(&self) -> Scheme {
        self.p.scheme
    }
    pub fn real_candidates(&self) -> usize {
        self.real_candidates
    }
    pub fn reciprocal_candidates(&self) -> usize {
        self.reciprocal_candidates
    }
    fn visit_real(
        &self,
        range: std::ops::Range<usize>,
        visit: &mut impl FnMut(f64),
    ) -> Result<(), PhysicsError> {
        let n = self.ions.len();
        let end = range.end.min(self.real_candidates);
        match self.p.scheme {
            Scheme::Reference => {
                for c in range.start..end {
                    let pair = c / self.rows;
                    let m0 = (c % self.rows) as i32 - self.p.rb[0];
                    reference_real_row(&self.p, self.ions, pair / n, pair % n, m0, visit)?;
                }
            }
            Scheme::Balanced => {
                for c in range.start..end {
                    let (i, j) = self.p.pairs[c];
                    visit(balanced_real_pair(&self.p, self.ions, i, j)?);
                }
            }
        }
        Ok(())
    }
    fn visit_reciprocal(&self, range: std::ops::Range<usize>, visit: &mut impl FnMut(f64)) {
        let end = range.end.min(self.reciprocal_candidates);
        match self.p.scheme {
            Scheme::Reference => {
                for c in range.start..end {
                    if let Some(term) = reciprocal_term(&self.p, self.ions, box_point(self.p.gb, c))
                    {
                        visit(term);
                    }
                }
            }
            Scheme::Balanced => {
                for c in range.start..end {
                    visit(balanced_reciprocal(&self.p, c));
                }
            }
        }
    }
    /// Appends the accepted real-space terms of candidates `range` to `out`.
    pub fn real_terms(
        &self,
        range: std::ops::Range<usize>,
        out: &mut Vec<f64>,
    ) -> Result<(), PhysicsError> {
        self.visit_real(range, &mut |t| out.push(t))
    }
    /// Appends the accepted reciprocal terms of candidates `range` to `out`.
    pub fn reciprocal_terms(&self, range: std::ops::Range<usize>, out: &mut Vec<f64>) {
        self.visit_reciprocal(range, &mut |t| out.push(t))
    }
    /// The energy from complete term lists (all candidates, serial order).
    pub fn finish(&self, real: &[f64], reciprocal: &[f64]) -> Result<EwaldEnergy, PhysicsError> {
        self.finish_sums(accumulate(0.0, real), accumulate(0.0, reciprocal))
    }
    /// The energy from the two sums accumulated in the serial order (for
    /// example by `accumulate` chained over contiguous candidate blocks).
    pub fn finish_sums(&self, real: f64, reciprocal: f64) -> Result<EwaldEnergy, PhysicsError> {
        finish(&self.p, self.ions, real, reciprocal)
    }
}

/// `start` plus `terms` added one at a time in order: the serial loop's
/// association, so chaining over the blocks of a partition from +0.0 gives the
/// serial sum bit for bit.
pub fn accumulate(start: f64, terms: &[f64]) -> f64 {
    let mut sum = start;
    for t in terms {
        sum += t;
    }
    sum
}

/// Cartesian -dE/dR in Ha/Bohr, including real and reciprocal terms.
/// Self and uniform-background energies have no position derivative at fixed cell.
/// Uses exactly the energy evaluator's truncation spheres and eta.
pub fn forces(
    cell: &Cell,
    ions: &[Ion],
    options: EwaldOptions,
) -> Result<Vec<[f64; 3]>, PhysicsError> {
    let parameters = energy(cell, ions, options)?;
    let p = prepare(cell, ions, options, scheme())?;
    let eta = parameters.eta;
    let a = p.a;
    let b = p.b;
    let gb = bounds(&a, parameters.reciprocal_cutoff_inverse_bohr, 0.0)?;
    let positions = &p.positions;
    let mut result = vec![[0.0; 3]; ions.len()];
    for (i, ion) in ions.iter().enumerate() {
        for (j, other) in ions.iter().enumerate() {
            let delta = positions[i] - positions[j];
            p.walker.images(&delta, |m| {
                if i == j && m == [0; 3] {
                    return Ok(());
                }
                let r = a * (positions[i] - positions[j] + Vector3::from(m.map(f64::from)));
                let distance = r.norm();
                if distance <= parameters.real_cutoff_bohr {
                    let factor = ion.charge
                        * other.charge
                        * (libm::erfc(eta * distance) / distance.powi(3)
                            + 2.0 * eta / PI.sqrt() * (-(eta * distance).powi(2)).exp()
                                / distance.powi(2));
                    for axis in 0..3 {
                        result[i][axis] += factor * r[axis];
                    }
                }
                Ok(())
            })?;
        }
    }
    enumerate(
        |m| {
            if m == [0; 3] {
                return Ok(());
            }
            let mv = Vector3::from(m.map(f64::from));
            let g = b * mv;
            let g2 = g.norm_squared();
            if g2 <= parameters.reciprocal_cutoff_inverse_bohr.powi(2) {
                let structure: Complex64 = ions
                    .iter()
                    .zip(positions)
                    .map(|(ion, s)| Complex64::from_polar(ion.charge, TAU * mv.dot(s)))
                    .sum();
                let prefactor =
                    4.0 * PI / cell.volume_bohr3() * (-g2 / (4.0 * eta * eta)).exp() / g2;
                for (i, ion) in ions.iter().enumerate() {
                    let phase = Complex64::from_polar(ion.charge, TAU * mv.dot(&positions[i]));
                    for axis in 0..3 {
                        result[i][axis] += prefactor * (structure.conj() * phase).im * g[axis];
                    }
                }
            }
            Ok(())
        },
        gb,
    )?;
    if result.iter().flatten().any(|v| !v.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
