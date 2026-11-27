//! Deterministic fixed-cell optimizers; all coordinates and steps are Cartesian Bohr.
use serde::{Deserialize, Serialize};

#[derive(Clone, Debug, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Optimizer {
    pub h: Vec<f64>,
    pub velocity: Vec<f64>,
    pub previous_x: Vec<f64>,
    pub previous_force: Vec<f64>,
    pub radius: f64,
    /// FIRE time step in ASE FIRE units, Angstrom/sqrt(eV) (unit masses,
    /// forces in eV/Angstrom, positions in Angstrom); see [`fire_time_scale`].
    pub dt: f64,
    /// Upper bound of `dt`, same units (ADR 0056).
    pub dt_max: f64,
    pub alpha: f64,
    pub positive: usize,
}
/// ASE FIRE defaults (ase.optimize.FIRE: dt 0.1, dtmax 1.0, in its eV/Angstrom
/// units with unit masses).
pub const FIRE_DT: f64 = 0.1;
pub const FIRE_DT_MAX: f64 = 1.0;
const FIRE_N_MIN: usize = 5;
const FIRE_F_INC: f64 = 1.1;
const FIRE_F_DEC: f64 = 0.5;
const FIRE_ALPHA_START: f64 = 0.1;
const FIRE_F_ALPHA: f64 = 0.99;
/// Factor that turns a FIRE time step in ASE units (Angstrom/sqrt(eV)) into
/// the native one (Bohr/sqrt(Ha)) such that the same dimensionless FIRE
/// parameters produce the same physical displacements (ADR 0056).
///
/// With unit masses an Euler step moves x = dt^2 F. In ASE units
/// x_A = dt_A^2 F_eV/A. With E = eV per Ha and L = Angstrom per Bohr,
/// F_eV/A = F_Ha/Bohr E/L and x_Bohr = x_A/L, so x_Bohr = dt_A^2 (E/L^2) F_Ha/Bohr:
/// the native step dt_n = dt_A sqrt(E)/L (about 9.858) reproduces the ASE
/// trajectory exactly (the velocity mixing and the power test are homogeneous
/// in velocity and force, hence unit-independent).
pub fn fire_time_scale() -> f64 {
    use dft_units::{Angstrom, Bohr, EnergyUnit, from_hartree};
    let ev_per_ha = from_hartree(1., EnergyUnit::ElectronVolt);
    let angstrom_per_bohr = Angstrom::from(Bohr(1.)).0;
    ev_per_ha.sqrt() / angstrom_per_bohr
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
impl Optimizer {
    /// BFGS spring model and FIRE state with the given FIRE time steps (ASE units).
    pub fn with_fire(n: usize, radius: f64, dt: f64, dt_max: f64) -> Self {
        Self {
            dt,
            dt_max,
            ..Self::new(n, radius)
        }
    }
    pub fn new(n: usize, radius: f64) -> Self {
        // Independent harmonic springs, k=1 Ha/Bohr^2, H^-1=I Bohr^2/Ha.
        let mut h = vec![0.; n * n];
        for i in 0..n {
            h[i * n + i] = 1.;
        }
        Self {
            h,
            velocity: vec![0.; n],
            previous_x: vec![],
            previous_force: vec![],
            radius,
            dt: FIRE_DT,
            dt_max: FIRE_DT_MAX,
            alpha: FIRE_ALPHA_START,
            positive: 0,
        }
    }
    pub fn step(&mut self, x: &[f64], f: &[f64], fire: bool, cap: f64) -> Vec<f64> {
        let n = x.len();
        let mut step = vec![0.; n];
        if fire {
            // Unit fictitious masses in ASE FIRE units (ADR 0056). The first
            // step only kicks (no power test), as ase.optimize.FIRE; the
            // driver zeroes the velocity after a retry, which the power test
            // then treats as a restart (dt halved, alpha reset).
            if !self.previous_x.is_empty() {
                let power = dot(&self.velocity, f);
                if power > 0. {
                    let ratio = (dot(&self.velocity, &self.velocity) / dot(f, f)).sqrt();
                    for (v, f) in self.velocity.iter_mut().zip(f) {
                        *v = (1. - self.alpha) * *v + self.alpha * ratio * f;
                    }
                    if self.positive > FIRE_N_MIN {
                        self.dt = (self.dt * FIRE_F_INC).min(self.dt_max);
                        self.alpha *= FIRE_F_ALPHA;
                    }
                    self.positive += 1;
                } else {
                    self.positive = 0;
                    self.dt *= FIRE_F_DEC;
                    self.alpha = FIRE_ALPHA_START;
                    self.velocity.fill(0.);
                }
            }
            let dt = self.dt * fire_time_scale();
            for (v, f) in self.velocity.iter_mut().zip(f) {
                *v += dt * f;
            }
            for (s, v) in step.iter_mut().zip(&self.velocity) {
                *s = dt * v;
            }
        } else {
            if !self.previous_x.is_empty() {
                let s: Vec<_> = x.iter().zip(&self.previous_x).map(|(a, b)| a - b).collect();
                let y: Vec<_> = self
                    .previous_force
                    .iter()
                    .zip(f)
                    .map(|(a, b)| a - b)
                    .collect();
                let sy = dot(&s, &y);
                if sy > 1e-12 * (dot(&s, &s) * dot(&y, &y)).sqrt() && sy > 1e-30 {
                    let hy: Vec<_> = self.h.chunks(n).map(|row| dot(row, &y)).collect();
                    let factor = (sy + dot(&y, &hy)) / (sy * sy);
                    for i in 0..n {
                        for j in 0..n {
                            self.h[i * n + j] +=
                                factor * s[i] * s[j] - (hy[i] * s[j] + s[i] * hy[j]) / sy;
                        }
                    }
                }
            }
            for (s, row) in step.iter_mut().zip(self.h.chunks(n)) {
                *s = dot(row, f);
            }
            if dot(&step, f) < 0. || step.iter().any(|v| !v.is_finite()) {
                *self = Self::new(n, self.radius);
                step.copy_from_slice(f);
            }
        }
        let max_atom = step.chunks(3).map(|v| dot(v, v).sqrt()).fold(0., f64::max);
        let limit = cap.min(self.radius);
        if max_atom > limit {
            for s in &mut step {
                *s *= limit / max_atom;
            }
        }
        self.previous_x = x.to_vec();
        self.previous_force = f.to_vec();
        step
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    /// Independent FIRE written in ASE units (positions in Angstrom, forces in
    /// eV/Angstrom, unit masses) following the published algorithm (Bitzek et
    /// al., PRL 97, 170201 (2006)) with the ASE parameter names and defaults;
    /// the step cap is left out (never binding in the comparison below).
    struct AseUnitsFire {
        v: Option<Vec<f64>>,
        dt: f64,
        a: f64,
        n_steps: usize,
    }
    impl AseUnitsFire {
        fn step(&mut self, f: &[f64]) -> Vec<f64> {
            match &mut self.v {
                None => self.v = Some(vec![0.; f.len()]),
                Some(v) => {
                    let vf = dot(v, f);
                    if vf > 0. {
                        let scale = (dot(v, v) / dot(f, f)).sqrt();
                        for (v, f) in v.iter_mut().zip(f) {
                            *v = (1. - self.a) * *v + self.a * f * scale;
                        }
                        if self.n_steps > 5 {
                            self.dt = (self.dt * 1.1).min(1.0);
                            self.a *= 0.99;
                        }
                        self.n_steps += 1;
                    } else {
                        v.fill(0.);
                        self.a = 0.1;
                        self.dt *= 0.5;
                        self.n_steps = 0;
                    }
                }
            }
            let v = self.v.as_mut().unwrap();
            for (v, f) in v.iter_mut().zip(f) {
                *v += self.dt * f;
            }
            v.iter().map(|v| self.dt * v).collect()
        }
    }
    #[test]
    fn fire_time_scale_matches_codata_units() {
        // sqrt(Hartree in eV) / (Bohr in Angstrom) = 9.8576...; its square is
        // the displacement ratio 97.17 of the audit (PRELAX-P0, 2.2).
        let s = fire_time_scale();
        assert!((s * s - 97.1736).abs() < 1e-3, "{}", s * s);
        // First FIRE step: a pure kick, x = (dt s)^2 F; no power test, so dt
        // is not halved on the first call (ASE behaviour).
        let mut opt = Optimizer::new(3, 10.);
        let f = [1e-3, -2e-3, 0.5e-3];
        let step = opt.step(&[0.; 3], &f, true, 10.);
        for (s, f) in step.iter().zip(f) {
            let expected = (FIRE_DT * fire_time_scale()).powi(2) * f;
            assert!(
                (s - expected).abs() <= 1e-15 * expected.abs(),
                "{s} {expected}"
            );
        }
        assert_eq!(opt.dt, FIRE_DT);
    }
    #[test]
    fn fire_trajectory_equals_ase_units_fire() {
        // Anisotropic harmonic wells (Ha/Bohr^2) with a coupling term, chosen so
        // that the run contains positive-power stretches (dt growth after
        // N_min) and at least one power reset. The cap is never binding.
        use dft_units::{Angstrom, Bohr, EnergyUnit, from_hartree};
        let ev = from_hartree(1., EnergyUnit::ElectronVolt);
        let ang = Angstrom::from(Bohr(1.)).0;
        let k = [2e-2, 1e-1, 5e-3, 3e-2, 6e-2, 1e-2];
        let force = |x: &[f64]| -> Vec<f64> {
            let mut f: Vec<f64> = x.iter().zip(&k).map(|(x, k)| -k * x).collect();
            f[0] -= 1e-2 * x[3];
            f[3] -= 1e-2 * x[0];
            f
        };
        let mut x: Vec<f64> = vec![0.3, -0.2, 0.25, 0.1, -0.15, 0.05];
        let mut x_ase: Vec<f64> = x.iter().map(|x| x * ang).collect();
        let mut native = Optimizer::new(6, 1e3);
        let mut reference = AseUnitsFire {
            v: None,
            dt: 0.1,
            a: 0.1,
            n_steps: 0,
        };
        let (mut resets, mut max_dt) = (0, 0f64);
        for _ in 0..60 {
            let f = force(&x);
            let positive_before = native.positive;
            let s = native.step(&x, &f, true, 1e3);
            if native.positive == 0 && positive_before > 0 {
                resets += 1;
            }
            let x_bohr: Vec<f64> = x_ase.iter().map(|x| x / ang).collect();
            let f_ase: Vec<f64> = force(&x_bohr).iter().map(|f| f * ev / ang).collect();
            let s_ase = reference.step(&f_ase);
            for ((x, s), (xa, sa)) in x.iter_mut().zip(&s).zip(x_ase.iter_mut().zip(&s_ase)) {
                *x += s;
                *xa += sa;
            }
            for (x, xa) in x.iter().zip(&x_ase) {
                assert!(
                    (x - xa / ang).abs() <= 1e-10 * (1. + x.abs()),
                    "native {x} Bohr vs ASE units {} Bohr",
                    xa / ang
                );
            }
            assert!((native.dt - reference.dt).abs() <= 1e-15);
            max_dt = max_dt.max(native.dt);
        }
        assert!(
            resets > 0 && max_dt > FIRE_DT,
            "path must exercise both branches: resets {resets} max dt {max_dt}"
        );
    }
    #[test]
    fn quadratic_and_rosenbrock() {
        for fire in [false, true] {
            for rosen in [false, true] {
                let mut x: Vec<f64> = vec![-1.2, 1., 0.];
                let mut opt = Optimizer::new(3, 0.1);
                for _ in 0..100000 {
                    let f = if rosen {
                        vec![
                            2. * (1. - x[0]) + 40. * x[0] * (x[1] - x[0] * x[0]),
                            -20. * (x[1] - x[0] * x[0]),
                            -x[2],
                        ]
                    } else {
                        vec![-x[0], -3. * x[1], -x[2]]
                    };
                    if f.iter().fold(0f64, |a, b| a.max(b.abs())) < 1e-8 {
                        break;
                    }
                    let s = opt.step(&x, &f, fire, 0.1);
                    assert!(dot(&s, &s).sqrt() <= 0.100000000000001);
                    for (x, s) in x.iter_mut().zip(s) {
                        *x += s;
                    }
                }
                let target = if rosen { 1. } else { 0. };
                assert!(
                    (x[0] - target).abs() < 1e-6 && (x[1] - target).abs() < 1e-6,
                    "{fire} {rosen} {x:?}"
                );
            }
        }
    }
}
