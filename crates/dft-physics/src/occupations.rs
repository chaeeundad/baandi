//! Finite-temperature number roots. Arrays are spin/k/band ordered; k weights
//! are applied only in sums. Callers own all output storage.
use std::fmt;

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Constraint {
    Total,
    FixedMoment(f64),
}
#[derive(Clone, Copy, Debug)]
pub struct Layout<'a> {
    pub capacities: &'a [f64],
    pub weights: &'a [f64],
    pub bands: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct FermiDirac {
    pub theta_ha: f64,
    pub electron_tolerance: f64,
    pub max_iterations: usize,
}
impl FermiDirac {
    pub fn new(theta_ha: f64) -> Self {
        Self {
            theta_ha,
            electron_tolerance: 1e-13,
            max_iterations: 256,
        }
    }
}
#[derive(Clone, Copy, Debug)]
pub struct Root {
    pub mu_ha: f64,
    pub bracket_ha: [f64; 2],
    pub iterations: usize,
}
#[derive(Clone, Copy, Debug)]
pub struct Evaluation {
    pub roots: [Option<Root>; 2],
    pub electron_count: f64,
    pub electron_residual: f64,
    pub moment_residual: Option<f64>,
    pub entropy_over_kb: f64,
    pub entropy_term_ha: f64,
    pub edge_max_fraction: f64,
    pub edge_electrons: f64,
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OccupationError(pub &'static str);
impl fmt::Display for OccupationError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(self.0)
    }
}
impl std::error::Error for OccupationError {}

// Both tails use exp(-abs(x)); no overflow and no log(0). Compute entropy from
// x rather than a rounded p=1 so the small hole entropy is retained.
fn fraction_entropy(x: f64) -> (f64, f64) {
    let a = x.abs();
    let t = (-a).exp();
    let q = t / (1.0 + t);
    let p = if x >= 0.0 { q } else { 1.0 / (1.0 + t) };
    let entropy = t.ln_1p() + if a.is_infinite() { 0.0 } else { a * q };
    (p, entropy)
}
impl FermiDirac {
    /// Expands a monotone bracket, then bisects until the number residual meets
    /// the relative target. A narrow bracket alone never certifies a root.
    /// No allocations, root hints, clipping, or endpoint approximation.
    pub fn evaluate(
        self,
        layout: Layout<'_>,
        energies: &[f64],
        target: f64,
        constraint: Constraint,
        output: &mut [f64],
    ) -> Result<Evaluation, OccupationError> {
        let fail = |s| Err(OccupationError(s));
        let ns = layout.capacities.len();
        let nk = layout.weights.len();
        let len = ns.checked_mul(nk).and_then(|v| v.checked_mul(layout.bands));
        if !matches!(layout.capacities, [2.0] | [1.0, 1.0])
            || nk == 0
            || layout.bands == 0
            || len != Some(energies.len())
            || output.len() != energies.len()
            || layout.weights.iter().any(|w| !w.is_finite() || *w <= 0.0)
            || (layout.weights.iter().sum::<f64>() - 1.0).abs() > 1e-14
            || energies.iter().any(|v| !v.is_finite())
            || energies
                .chunks(layout.bands)
                .any(|b| b.windows(2).any(|w| w[0] > w[1]))
        {
            return fail("occupations.layout: invalid sorted spin/k/band layout or weights");
        }
        if !self.theta_ha.is_finite() || self.theta_ha <= 0.0 {
            return fail(
                "occupations.smearing: require positive finite theta in Ha; use fixed occupations for zero temperature",
            );
        }
        if !target.is_finite()
            || !self.electron_tolerance.is_finite()
            || self.electron_tolerance <= 0.0
            || self.electron_tolerance > 1e-13
            || self.max_iterations == 0
        {
            return fail("occupations.root: invalid target, tolerance or iteration limit");
        }
        let sectors = match constraint {
            Constraint::Total => 1,
            Constraint::FixedMoment(m) if ns == 2 && m.is_finite() => 2,
            _ => {
                return fail(
                    "occupations.constraint: fixed moment requires two collinear channels",
                );
            }
        };
        let mut result = Evaluation {
            roots: [None; 2],
            electron_count: 0.0,
            electron_residual: 0.0,
            moment_residual: None,
            entropy_over_kb: 0.0,
            entropy_term_ha: 0.0,
            edge_max_fraction: 0.0,
            edge_electrons: 0.0,
        };
        let mut channel_counts = [0.0; 2];
        for sector in 0..sectors {
            let wanted = match constraint {
                Constraint::Total => target,
                Constraint::FixedMoment(m) => (target + if sector == 0 { m } else { -m }) / 2.0,
            };
            let first = if sectors == 1 { 0 } else { sector };
            let last = if sectors == 1 { ns } else { sector + 1 };
            let start = first * nk * layout.bands;
            let end = last * nk * layout.bands;
            let capacity = layout.capacities[first..last].iter().sum::<f64>()
                * layout.bands as f64
                * layout.weights.iter().sum::<f64>();
            if wanted <= 0.0 || wanted >= capacity {
                return fail(
                    "occupations.electron_target: no finite chemical potential at this sector endpoint; use fixed occupations or an interior target",
                );
            }
            let number = |mu: f64| {
                (start..end)
                    .map(|i| {
                        let s = i / (nk * layout.bands);
                        let k = (i / layout.bands) % nk;
                        layout.weights[k]
                            * layout.capacities[s]
                            * fraction_entropy((energies[i] - mu) / self.theta_ha).0
                    })
                    .sum::<f64>()
            };
            let low_energy = energies[start..end]
                .iter()
                .copied()
                .fold(f64::INFINITY, f64::min);
            let high_energy = energies[start..end]
                .iter()
                .copied()
                .fold(f64::NEG_INFINITY, f64::max);
            let mut width = self.theta_ha;
            let mut lo = low_energy;
            let mut hi = high_energy;
            let mut bracketed = false;
            for _ in 0..self.max_iterations {
                lo = low_energy - width;
                hi = high_energy + width;
                if !lo.is_finite() || !hi.is_finite() {
                    break;
                }
                if number(lo) <= wanted && number(hi) >= wanted {
                    bracketed = true;
                    break;
                }
                width *= 2.0;
            }
            if !bracketed {
                return fail("occupations.root: unable to bracket finite chemical potential");
            }
            let mut root = None;
            for iteration in 1..=self.max_iterations {
                let mu = 0.5 * lo + 0.5 * hi;
                let count = number(mu);
                if (count - wanted).abs() <= self.electron_tolerance * wanted {
                    root = Some(Root {
                        mu_ha: mu,
                        bracket_ha: [lo, hi],
                        iterations: iteration,
                    });
                    break;
                }
                if mu == lo || mu == hi {
                    break;
                }
                if count < wanted {
                    lo = mu;
                } else {
                    hi = mu;
                }
            }
            let root = root.ok_or(OccupationError(
                "occupations.root: number residual did not converge",
            ))?;
            result.roots[sector] = Some(root);
            for i in start..end {
                let s = i / (nk * layout.bands);
                let k = (i / layout.bands) % nk;
                let (p, entropy) = fraction_entropy((energies[i] - root.mu_ha) / self.theta_ha);
                output[i] = layout.capacities[s] * p;
                let weighted = layout.weights[k] * output[i];
                channel_counts[s] += weighted;
                result.entropy_over_kb += layout.weights[k] * layout.capacities[s] * entropy;
                if i % layout.bands == layout.bands - 1 {
                    result.edge_max_fraction = result.edge_max_fraction.max(p);
                    result.edge_electrons += weighted;
                }
            }
        }
        result.electron_count = channel_counts.iter().sum();
        result.electron_residual = result.electron_count - target;
        if let Constraint::FixedMoment(m) = constraint {
            result.moment_residual = Some(channel_counts[0] - channel_counts[1] - m);
        }
        result.entropy_term_ha = self.theta_ha * result.entropy_over_kb;
        if !result.entropy_term_ha.is_finite() {
            return fail("occupations.entropy: unrepresentable TS");
        }
        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // Independent logistic/entropy oracle, intentionally not production helpers.
    fn analytic(e: f64, mu: f64, t: f64) -> f64 {
        0.5 * (1.0 - ((e - mu) / (2.0 * t)).tanh())
    }
    fn entropy(p: f64) -> f64 {
        if p == 0.0 || p == 1.0 {
            0.0
        } else {
            -p * p.ln() - (1.0 - p) * (1.0 - p).ln()
        }
    }
    #[test]
    fn symmetric_and_degenerate_oracles() {
        for a in [0.0, 0.01, 1.0] {
            for theta in [0.001, 0.01, 0.1] {
                let mut f = [0.0; 2];
                let r = FermiDirac::new(theta)
                    .evaluate(
                        Layout {
                            capacities: &[2.0],
                            weights: &[1.0],
                            bands: 2,
                        },
                        &[-a, a],
                        2.0,
                        Constraint::Total,
                        &mut f,
                    )
                    .unwrap();
                assert_eq!(r.roots[0].unwrap().mu_ha, 0.0);
                assert!((f.iter().sum::<f64>() - 2.0).abs() < 1e-12);
                for (i, e) in [-a, a].into_iter().enumerate() {
                    assert!((f[i] - 2.0 * analytic(e, 0.0, theta)).abs() < 1e-12);
                }
                let expected = [-a, a]
                    .iter()
                    .map(|e| 2.0 * entropy(analytic(*e, 0.0, theta)))
                    .sum::<f64>();
                assert!((r.entropy_over_kb - expected).abs() < 1e-10);
                if a == 0.0 {
                    assert!((r.entropy_over_kb - 4.0 * 2.0_f64.ln()).abs() < 1e-10);
                }
            }
        }
    }
    #[test]
    fn unequal_weights_and_fixed_moment() {
        let e = [-0.3, 0.2, -0.1, 0.4, -0.2, 0.3, 0.0, 0.5];
        let layout = Layout {
            capacities: &[1.0, 1.0],
            weights: &[0.3, 0.7],
            bands: 2,
        };
        for constraint in [Constraint::Total, Constraint::FixedMoment(0.4)] {
            let mut f = [0.0; 8];
            let r = FermiDirac::new(0.1)
                .evaluate(layout, &e, 1.7, constraint, &mut f)
                .unwrap();
            assert!(r.electron_residual.abs() < 1e-12);
            assert!(r.moment_residual.unwrap_or(0.0).abs() < 1e-12);
            let mut counts = [0.0; 2];
            for i in 0..8 {
                let s = i / 4;
                let mu = r.roots[if constraint == Constraint::Total {
                    0
                } else {
                    s
                }]
                .unwrap()
                .mu_ha;
                assert!((f[i] - analytic(e[i], mu, 0.1)).abs() < 1e-12);
                counts[s] += layout.weights[(i / 2) % 2] * f[i];
            }
            assert!((counts.iter().sum::<f64>() - 1.7).abs() < 1e-12);
        }
    }
    #[test]
    fn rejects_endpoints_and_invalid_parameters() {
        let layout = Layout {
            capacities: &[2.0],
            weights: &[1.0],
            bands: 2,
        };
        for n in [0.0, 4.0, -1.0, f64::NAN] {
            assert!(
                FermiDirac::new(0.01)
                    .evaluate(layout, &[-1.0, 1.0], n, Constraint::Total, &mut [0.0; 2])
                    .is_err()
            );
        }
        for t in [0.0, -1.0, f64::INFINITY, f64::NAN] {
            assert!(
                FermiDirac::new(t)
                    .evaluate(layout, &[-1.0, 1.0], 2.0, Constraint::Total, &mut [0.0; 2])
                    .is_err()
            );
        }
        assert!(
            FermiDirac {
                max_iterations: 1,
                ..FermiDirac::new(0.1)
            }
            .evaluate(layout, &[-1.0, 1.0], 1.3, Constraint::Total, &mut [0.0; 2])
            .is_err()
        );
    }
}
