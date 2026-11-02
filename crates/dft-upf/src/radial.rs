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
