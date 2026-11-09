//! Direct radial quadrature; lengths are Bohr and wave numbers Bohr^-1.
use crate::Pseudopotential;
use std::f64::consts::PI;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum RadialError {
    #[error("invalid radial mesh: {0}")]
    Mesh(&'static str),
    #[error("invalid radial argument: {0}")]
    Argument(&'static str),
}
#[derive(Debug)]
pub struct RadialMesh {
    r: Vec<f64>,
    rab: Vec<f64>,
}
impl RadialMesh {
    pub fn new(r: Vec<f64>, rab: Vec<f64>) -> Result<Self, RadialError> {
        if r.len() < 2 || r.len() != rab.len() {
            return Err(RadialError::Mesh("need equal lengths >= 2"));
        }
        if r.iter().any(|x| !x.is_finite() || *x < 0.0) || r.windows(2).any(|w| w[1] <= w[0]) {
            return Err(RadialError::Mesh(
                "r must be finite, nonnegative and strictly increasing",
            ));
        }
        if rab.iter().any(|x| !x.is_finite() || *x <= 0.0) {
            return Err(RadialError::Mesh("rab must be finite and positive"));
        }
        Ok(Self { r, rab })
    }
    pub fn r(&self) -> &[f64] {
        &self.r
    }
    pub fn rab(&self) -> &[f64] {
        &self.rab
    }
    /// Integrate sampled f(r) dr using Simpson weights in the mesh index.
    pub fn integrate(&self, values: &[f64]) -> Result<f64, RadialError> {
        if values.len() != self.r.len() || values.iter().any(|x| !x.is_finite()) {
            return Err(RadialError::Argument(
                "integrand length or non-finite values",
            ));
        }
        finite(integrate(values.len(), |i| values[i] * self.rab[i]))
    }
}
fn finite(x: f64) -> Result<f64, RadialError> {
    if x.is_finite() {
        Ok(x)
    } else {
        Err(RadialError::Argument("unrepresentable result"))
    }
}
// Composite 1/3 rule; even point counts use a final 3/8 panel.
// Two points alone use trapezoidal integration in the index coordinate.
fn integrate(n: usize, f: impl Fn(usize) -> f64) -> f64 {
    if n < 2 {
        return 0.0;
    }
    if n == 2 {
        return (f(0) + f(1)) * 0.5;
    }
    let last = if n.is_multiple_of(2) { n - 4 } else { n - 1 };
    let mut sum = 0.0;
    let mut correction = 0.0;
    let mut add = |x: f64| {
        let y = x - correction;
        let t = sum + y;
        correction = (t - sum) - y;
        sum = t;
    };
    for i in (0..last).step_by(2) {
        add((f(i) + 4.0 * f(i + 1) + f(i + 2)) / 3.0);
    }
    if n.is_multiple_of(2) {
        add(3.0 * (f(last) + 3.0 * f(last + 1) + 3.0 * f(last + 2) + f(last + 3)) / 8.0);
    }
    sum
}
/// Spherical Bessel j_l, l <= 3, x >= 0.
pub fn spherical_bessel(l: usize, x: f64) -> Result<f64, RadialError> {
    if l > 3 || !x.is_finite() || x < 0.0 {
        return Err(RadialError::Argument("j_l requires l <= 3, finite x >= 0"));
    }
    Ok(bessel(l, x))
}
fn bessel(l: usize, x: f64) -> f64 {
    if x < 4.0 {
        let mut term = 1.0;
        for k in 1..=l {
            term *= x / (2 * k + 1) as f64;
        }
        let mut sum = term;
        for k in 1..80 {
            term *= -x * x / (2 * k * (2 * l + 2 * k + 1)) as f64;
            sum += term;
            if term.abs() <= f64::EPSILON * sum.abs() {
                break;
            }
        }
        return sum;
    }
    let (sin, cos) = x.sin_cos();
    let mut prev = sin / x;
    if l == 0 {
        return prev;
    }
    let mut cur = (prev - cos) / x;
    for k in 1..l {
        let next = (2 * k + 1) as f64 / x * cur - prev;
        prev = cur;
        cur = next;
    }
    cur
}
/// Real normalized harmonics. m=-l..l; negative m are sine, positive cosine.
/// Includes Condon-Shortley P_l^m and an additional (-1)^m in real combinations.
/// Direction is normalized internally; zero/non-finite directions are rejected.
pub fn real_spherical_harmonic(
    l: usize,
    m: isize,
    direction: [f64; 3],
) -> Result<f64, RadialError> {
    if l > 3 || m.unsigned_abs() > l || direction.iter().any(|v| !v.is_finite()) {
        return Err(RadialError::Argument("Y_lm index or direction"));
    }
    let norm = direction[0].hypot(direction[1]).hypot(direction[2]);
    if norm == 0.0 || !norm.is_finite() {
        return Err(RadialError::Argument("zero or unrepresentable direction"));
    }
    let z = (direction[2] / norm).clamp(-1.0, 1.0);
    let phi = direction[1].atan2(direction[0]);
    let a = m.unsigned_abs();
    let mut p = 1.0;
    for k in 1..=a {
        p *= -((2 * k - 1) as f64) * ((1.0 - z) * (1.0 + z)).sqrt();
    }
    if l > a {
        let mut prev = p;
        p *= (2 * a + 1) as f64 * z;
        for k in a + 2..=l {
            let next = ((2 * k - 1) as f64 * z * p - (k + a - 1) as f64 * prev) / (k - a) as f64;
            prev = p;
            p = next;
        }
    }
    let ratio = ((l - a + 1)..=l + a).fold(1.0, |v, k| v / k as f64);
    let norm = ((2 * l + 1) as f64 / (4.0 * PI) * ratio).sqrt();
    let angular = if m == 0 {
        1.0
    } else {
        let sign = if a.is_multiple_of(2) { 1.0 } else { -1.0 };
        2.0_f64.sqrt()
            * sign
            * if m < 0 {
                (a as f64 * phi).sin()
            } else {
                (a as f64 * phi).cos()
            }
    };
    Ok(norm * p * angular)
}

/// Local transform view. The endpoint is the last mesh point <= requested r_max.
pub struct LocalTransform<'a> {
    mesh: &'a RadialMesh,
    v: &'a [f64],
    z: f64,
    volume: f64,
    n: usize,
}
impl<'a> LocalTransform<'a> {
    pub fn new(
        mesh: &'a RadialMesh,
        v_hartree: &'a [f64],
        z: f64,
        volume: f64,
        r_max: Option<f64>,
    ) -> Result<Self, RadialError> {
        if v_hartree.len() != mesh.r.len()
            || v_hartree.iter().any(|v| !v.is_finite())
            || !z.is_finite()
            || z <= 0.0
            || !volume.is_finite()
            || volume <= 0.0
        {
            return Err(RadialError::Argument("local samples, charge or volume"));
        }
        let end = *mesh.r.last().expect("validated mesh");
        let limit = r_max.unwrap_or(end.min(10.0));
        if !limit.is_finite() || limit <= 0.0 || limit > end {
            return Err(RadialError::Argument("r_max outside mesh"));
        }
        let n = mesh.r.partition_point(|r| *r <= limit);
        if n < 2 {
            return Err(RadialError::Argument(
                "r_max includes fewer than two points",
            ));
        }
        Ok(Self {
            mesh,
            v: v_hartree,
            z,
            volume,
            n,
        })
    }
    pub fn effective_r_max(&self) -> f64 {
        self.mesh.r[self.n - 1]
    }
    fn short_integral(&self, g: f64) -> f64 {
        integrate(self.n, |i| {
            let r = self.mesh.r[i];
            // r^2[V + Z erf(r)/r], finite even at r=0; no division by r.
            (r * r * self.v[i] + self.z * r * libm::erf(r)) * bessel(0, g * r) * self.mesh.rab[i]
        })
    }
    /// V_loc(G), in Ha, for strictly positive G in Bohr^-1.
    pub fn at_g(&self, g: f64) -> Result<f64, RadialError> {
        if !g.is_finite() || g <= 0.0 || !(g * self.effective_r_max()).is_finite() {
            return Err(RadialError::Argument("G must be finite and positive"));
        }
        finite(
            4.0 * PI / self.volume
                * (self.short_integral(g) - self.z * (-g * g / 4.0).exp() / (g * g)),
        )
    }
    /// Finite Coulomb-subtracted G=0 constant, in Ha. Adds pi*Z/Omega
    /// to the erf-subtracted integral, including its analytic tail to infinity.
    pub fn alpha_z(&self) -> Result<f64, RadialError> {
        finite(4.0 * PI / self.volume * self.short_integral(0.0) + PI * self.z / self.volume)
    }
}
impl Pseudopotential {
    pub fn local_transform(
        &self,
        volume_bohr3: f64,
        r_max: Option<f64>,
    ) -> Result<LocalTransform<'_>, RadialError> {
        LocalTransform::new(
            &self.mesh,
            &self.local,
            self.z_valence(),
            volume_bohr3,
            r_max,
        )
    }
    pub fn projector_transform(&self, index: usize, q: f64) -> Result<f64, RadialError> {
        let p = self
            .projectors()
            .get(index)
            .ok_or(RadialError::Argument("projector index"))?;
        projector_transform(&self.mesh, p.l, &p.radial_times_r, q)
    }
}
/// 4*pi*integral r*u(r)*j_l(qr) dr, u=r*beta in Bohr^-1/2.
/// Samples may be a prefix ending at the projector cutoff, with zero beyond.
/// No volume, angular, phase, or structure-factor prefactors are included.
pub fn projector_transform(
    mesh: &RadialMesh,
    l: usize,
    radial_times_r: &[f64],
    q: f64,
) -> Result<f64, RadialError> {
    let n = radial_times_r.len();
    if l > 3
        || !q.is_finite()
        || q < 0.0
        || n == 0
        || n > mesh.r.len()
        || radial_times_r.iter().any(|x| !x.is_finite())
        || !(q * mesh.r[n - 1]).is_finite()
    {
        return Err(RadialError::Argument("projector samples, l or q"));
    }
    finite(
        4.0 * PI
            * integrate(n, |i| {
                mesh.r[i] * radial_times_r[i] * bessel(l, q * mesh.r[i]) * mesh.rab[i]
            }),
    )
}

/// Fourier transform of a spherical density rho(r), in electrons/Bohr^3.
/// Returns electrons: 4*pi integral r^2 rho(r) j0(qr) dr. Unlike PP_RHOATOM,
/// PP_NLCC is a density, not a radial charge already multiplied by 4*pi*r^2.
/// Uses the same index-coordinate Simpson/3/8 quadrature as local/projectors.
pub fn spherical_density_transform(
    mesh: &RadialMesh,
    rho: &[f64],
    q: f64,
) -> Result<f64, RadialError> {
    if rho.len() != mesh.r.len()
        || rho.iter().any(|v| !v.is_finite())
        || !q.is_finite()
        || q < 0.0
        || !(q * mesh.r[mesh.r.len() - 1]).is_finite()
    {
        return Err(RadialError::Argument("density samples or q"));
    }
    finite(
        4.0 * PI
            * integrate(rho.len(), |i| {
                let r = mesh.r[i];
                r * r * rho[i] * bessel(0, q * r) * mesh.rab[i]
            }),
    )
}

/// Research-only radial endpoint selection. Defaults preserve existing numerics.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum LocalMesh {
    #[default]
    Baandi,
    QeLike,
    /// `QeLike` endpoint plus the reference code's G=0 convention: the constant
    /// is the plain mesh integral of r (r V + Z), which omits the first interval
    /// [0, r_1] of a mesh that does not start at r = 0 (BAANDI integrates the
    /// Coulomb part analytically and includes it).
    QeLikeG0,
}
impl LocalMesh {
    pub fn name(self) -> &'static str {
        match self {
            Self::Baandi => "baandi",
            Self::QeLike => "qe_like",
            Self::QeLikeG0 => "qe_like_g0",
        }
    }
}
/// Research-only even-count quadrature; odd counts always use the existing rule.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub enum RadialEvenRule {
    #[default]
    Simpson3_8,
    SimpsonQuadraticEnd,
}
impl RadialEvenRule {
    pub fn name(self) -> &'static str {
        match self {
            Self::Simpson3_8 => "simpson_3_8",
            Self::SimpsonQuadraticEnd => "simpson_quadratic_end",
        }
    }
}
fn integrate_rule(n: usize, f: impl Fn(usize) -> f64, rule: RadialEvenRule) -> f64 {
    if rule == RadialEvenRule::Simpson3_8 || !n.is_multiple_of(2) || n < 4 {
        return integrate(n, f);
    }
    integrate(n - 1, &f) + (-f(n - 3) + 8.0 * f(n - 2) + 5.0 * f(n - 1)) / 12.0
}
impl RadialMesh {
    pub fn integrate_with_rule(
        &self,
        values: &[f64],
        rule: RadialEvenRule,
    ) -> Result<f64, RadialError> {
        if rule == RadialEvenRule::Simpson3_8 {
            return self.integrate(values);
        }
        if values.len() != self.r.len() || values.iter().any(|x| !x.is_finite()) {
            return Err(RadialError::Argument(
                "integrand length or non-finite values",
            ));
        }
        finite(integrate_rule(
            values.len(),
            |i| values[i] * self.rab[i],
            rule,
        ))
    }
}
/// Additive research view; the original LocalTransform implementation is unchanged.
pub struct ConventionLocalTransform<'a> {
    base: LocalTransform<'a>,
    rule: RadialEvenRule,
    g0_plain: bool,
}
impl ConventionLocalTransform<'_> {
    pub fn effective_r_max(&self) -> f64 {
        self.base.effective_r_max()
    }
    fn short_integral(&self, g: f64) -> f64 {
        let b = &self.base;
        integrate_rule(
            b.n,
            |i| {
                let r = b.mesh.r[i];
                (r * r * b.v[i] + b.z * r * libm::erf(r)) * bessel(0, g * r) * b.mesh.rab[i]
            },
            self.rule,
        )
    }
    pub fn at_g(&self, g: f64) -> Result<f64, RadialError> {
        if self.rule == RadialEvenRule::Simpson3_8 {
            return self.base.at_g(g);
        }
        if !g.is_finite() || g <= 0.0 || !(g * self.effective_r_max()).is_finite() {
            return Err(RadialError::Argument("G must be finite and positive"));
        }
        finite(
            4.0 * PI / self.base.volume
                * (self.short_integral(g) - self.base.z * (-g * g / 4.0).exp() / (g * g)),
        )
    }
    pub fn alpha_z(&self) -> Result<f64, RadialError> {
        if self.g0_plain {
            let b = &self.base;
            // Plain mesh integral of r (r V + Z), no analytic head or tail.
            let plain = integrate_rule(
                b.n,
                |i| {
                    let r = b.mesh.r[i];
                    r * (r * b.v[i] + b.z) * b.mesh.rab[i]
                },
                self.rule,
            );
            return finite(4.0 * PI / b.volume * plain);
        }
        if self.rule == RadialEvenRule::Simpson3_8 {
            return self.base.alpha_z();
        }
        finite(
            4.0 * PI / self.base.volume * self.short_integral(0.0)
                + PI * self.base.z / self.base.volume,
        )
    }
}
impl Pseudopotential {
    pub fn local_transform_with_conventions(
        &self,
        volume: f64,
        local_mesh: LocalMesh,
        rule: RadialEvenRule,
    ) -> Result<ConventionLocalTransform<'_>, RadialError> {
        let mut base = self.local_transform(volume, None)?;
        if local_mesh != LocalMesh::Baandi {
            let n = self
                .mesh
                .r
                .iter()
                .position(|&r| r > 10.0)
                .map_or(self.mesh.r.len(), |i| i + 1);
            base.n = 2 * n.div_ceil(2) - 1;
            if base.n < 2 {
                return Err(RadialError::Argument(
                    "local mesh includes fewer than two points",
                ));
            }
        }
        Ok(ConventionLocalTransform {
            base,
            rule,
            g0_plain: local_mesh == LocalMesh::QeLikeG0,
        })
    }
    pub fn projector_transform_with_rule(
        &self,
        index: usize,
        q: f64,
        rule: RadialEvenRule,
    ) -> Result<f64, RadialError> {
        if rule == RadialEvenRule::Simpson3_8 {
            return self.projector_transform(index, q);
        }
        let p = self
            .projectors()
            .get(index)
            .ok_or(RadialError::Argument("projector index"))?;
        projector_transform_with_rule(&self.mesh, p.l, &p.radial_times_r, q, rule)
    }
}

pub fn projector_transform_with_rule(
    mesh: &RadialMesh,
    l: usize,
    radial_times_r: &[f64],
    q: f64,
    rule: RadialEvenRule,
) -> Result<f64, RadialError> {
    if rule == RadialEvenRule::Simpson3_8 {
        return projector_transform(mesh, l, radial_times_r, q);
    }
    let n = radial_times_r.len();
    if l > 3
        || !q.is_finite()
        || q < 0.0
        || n == 0
        || n > mesh.r.len()
        || radial_times_r.iter().any(|x| !x.is_finite())
        || !(q * mesh.r[n - 1]).is_finite()
    {
        return Err(RadialError::Argument("projector samples, l or q"));
    }
    finite(
        4.0 * PI
            * integrate_rule(
                n,
                |i| mesh.r[i] * radial_times_r[i] * bessel(l, q * mesh.r[i]) * mesh.rab[i],
                rule,
            ),
    )
}

pub fn spherical_density_transform_with_rule(
    mesh: &RadialMesh,
    rho: &[f64],
    q: f64,
    rule: RadialEvenRule,
) -> Result<f64, RadialError> {
    if rule == RadialEvenRule::Simpson3_8 {
        return spherical_density_transform(mesh, rho, q);
    }
    if rho.len() != mesh.r.len()
        || rho.iter().any(|v| !v.is_finite())
        || !q.is_finite()
        || q < 0.0
        || !(q * mesh.r[mesh.r.len() - 1]).is_finite()
    {
        return Err(RadialError::Argument("density samples or q"));
    }
    finite(
        4.0 * PI
            * integrate_rule(
                rho.len(),
                |i| {
                    let r = mesh.r[i];
                    r * r * rho[i] * bessel(0, q * r) * mesh.rab[i]
                },
                rule,
            ),
    )
}

#[cfg(test)]
mod core_tests {
    use super::*;
    #[test]
    fn gaussian_core_fourier_oracle() {
        // rho = Q*(a/pi)^(3/2)*exp(-a*r^2); transform Q*exp(-q^2/(4a)).
        // Frozen MF06 relative cap 1e-10; both odd and even mesh lengths.
        for n in [2001, 2002] {
            let h = 0.005;
            let mesh = RadialMesh::new((0..n).map(|i| i as f64 * h).collect(), vec![h; n]).unwrap();
            for a in [0.5, 1.0, 2.0] {
                let rho: Vec<_> = mesh
                    .r()
                    .iter()
                    .map(|r| 3.0 * (a / PI).powf(1.5) * (-a * r * r).exp())
                    .collect();
                for q in [0.0, 0.3, 1.0, 3.0] {
                    let actual = spherical_density_transform(&mesh, &rho, q).unwrap();
                    let expected = 3.0 * (-q * q / (4.0 * a)).exp();
                    let error = (actual - expected).abs() / expected;
                    println!("core Gaussian n={n} a={a} q={q} relative_error={error:e}");
                    assert!(error <= 1e-10);
                }
            }
        }
    }
}

#[cfg(test)]
mod convention_tests {
    use super::*;
    #[test]
    fn polynomial_degree_and_odd_identity() {
        for n in 3..20 {
            for h in [0.125_f64, 1.0, 2.0] {
                for degree in 0..=3 {
                    let f = |i: usize| (i as f64 * h).powi(degree) * h;
                    let exact = ((n - 1) as f64 * h).powi(degree + 1) / f64::from(degree + 1);
                    let cubic = integrate_rule(n, f, RadialEvenRule::Simpson3_8);
                    let quadratic = integrate_rule(n, f, RadialEvenRule::SimpsonQuadraticEnd);
                    // Simpson panels integrate cubics; the one-interval quadratic
                    // interpolant integrates quadratics and has +h^4/4 cubic error.
                    let expected_error = if n % 2 == 0 && degree == 3 {
                        h.powi(4) / 4.0
                    } else {
                        0.0
                    };
                    assert!((cubic - exact).abs() < 1e-12 * exact.max(1.0));
                    assert!((quadratic - exact - expected_error).abs() < 1e-12 * exact.max(1.0));
                    if n % 2 == 1 {
                        assert_eq!(cubic.to_bits(), quadratic.to_bits());
                    }
                }
            }
        }
    }
    #[test]
    fn specified_last_interval_weights() {
        // Scale impulse samples by 12 so all expected arithmetic is exact binary.
        for n in [4, 6, 10] {
            for (offset, weight) in [(3, -1.0), (2, 8.0), (1, 5.0)] {
                let impulse = |i| if i == n - offset { 12.0 } else { 0.0 };
                let panel = integrate(n - 1, impulse);
                assert_eq!(
                    integrate_rule(n, impulse, RadialEvenRule::SimpsonQuadraticEnd) - panel,
                    weight
                );
            }
        }
        assert_eq!(
            integrate_rule(2, |i| i as f64, RadialEvenRule::SimpsonQuadraticEnd),
            0.5
        );
    }
}
