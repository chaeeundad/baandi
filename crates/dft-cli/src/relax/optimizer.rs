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
    pub dt: f64,
    pub alpha: f64,
    pub positive: usize,
}
fn dot(a: &[f64], b: &[f64]) -> f64 {
    a.iter().zip(b).map(|(a, b)| a * b).sum()
}
impl Optimizer {
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
            dt: 0.1,
            alpha: 0.1,
            positive: 0,
        }
    }
    pub fn step(&mut self, x: &[f64], f: &[f64], fire: bool, cap: f64) -> Vec<f64> {
        let n = x.len();
        let mut step = vec![0.; n];
        if fire {
            // Unit fictitious masses; inspect power before the Euler kick.
            let power = dot(&self.velocity, f);
            if power > 0. {
                self.positive += 1;
                let ratio = (dot(&self.velocity, &self.velocity) / dot(f, f)).sqrt();
                for (v, f) in self.velocity.iter_mut().zip(f) {
                    *v = (1. - self.alpha) * *v + self.alpha * ratio * f;
                }
                if self.positive > 5 {
                    self.dt = (self.dt * 1.1).min(1.);
                    self.alpha *= 0.99;
                }
            } else {
                self.positive = 0;
                self.dt *= 0.5;
                self.alpha = 0.1;
                self.velocity.fill(0.);
            }
            for (v, f) in self.velocity.iter_mut().zip(f) {
                *v += self.dt * f;
            }
            for (s, v) in step.iter_mut().zip(&self.velocity) {
                *s = self.dt * v;
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
