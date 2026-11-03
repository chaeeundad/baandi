// Independent direct DFT oracle; no production FFT code is used here.
use dft_solvers::{SolverError, kerker, mixing};
use kerker::{DensityTransform, Kerker, PreconditionedPulay, ResidualPreconditioner};
use mixing::{DensityMixer, PulayMixer};
use num_complex::Complex64;
use std::{f64::consts::TAU, time::Instant};
struct Dft {
    shape: [usize; 3],
    scratch: Vec<Complex64>,
}
impl Dft {
    fn new(shape: [usize; 3]) -> Self {
        Self {
            shape,
            scratch: vec![Complex64::default(); shape.iter().product()],
        }
    }
}
fn coord(i: usize, s: [usize; 3]) -> [usize; 3] {
    [i / (s[1] * s[2]), i / s[2] % s[1], i % s[2]]
}
impl DensityTransform for Dft {
    fn transform(&mut self, v: &mut [Complex64], inverse: bool) -> Result<(), SolverError> {
        for (a, out) in self.scratch.iter_mut().enumerate() {
            let m = coord(a, self.shape);
            *out = v
                .iter()
                .enumerate()
                .map(|(b, &z)| {
                    let n = coord(b, self.shape);
                    let phase = TAU
                        * (0..3)
                            .map(|d| (m[d] * n[d]) as f64 / self.shape[d] as f64)
                            .sum::<f64>();
                    z * Complex64::from_polar(1.0, if inverse { phase } else { -phase })
                })
                .sum::<Complex64>()
                / if inverse { 1.0 } else { v.len() as f64 };
        }
        v.copy_from_slice(&self.scratch);
        Ok(())
    }
}
fn cell(stretch: f64) -> [[f64; 3]; 3] {
    [
        [0.7, 0.12, 0.0],
        [0.0, 0.9, 0.13],
        [0.0, 0.0, 0.6 / stretch],
    ]
}
#[test]
fn reciprocal_modes_charge_and_nyquist_in_skew_cells() {
    for shape in [[5, 3, 7], [4, 4, 6]] {
        let b = cell(4.0);
        let size = shape.iter().product();
        let mut p = Kerker::new(shape, b, 0.8, Dft::new(shape)).unwrap();
        assert_eq!(p.weights()[0], 0.0);
        for i in 0..size {
            let c = coord(i, shape);
            let neg = [
                (shape[0] - c[0]) % shape[0],
                (shape[1] - c[1]) % shape[1],
                (shape[2] - c[2]) % shape[2],
            ];
            let j = (neg[0] * shape[1] + neg[1]) * shape[2] + neg[2];
            assert!((p.weights()[i] - p.weights()[j]).abs() < 2e-16);
        }
        let mode = [1.0, -1.0, 1.0];
        let g2: f64 = b
            .iter()
            .map(|row| {
                row.iter()
                    .zip(mode)
                    .map(|(v, m)| v * m)
                    .sum::<f64>()
                    .powi(2)
            })
            .sum();
        let factor = g2 / (g2 + 0.64);
        let r: Vec<_> = (0..size)
            .map(|i| {
                let c = coord(i, shape);
                2.0 + (TAU
                    * (0..3)
                        .map(|a| mode[a] * c[a] as f64 / shape[a] as f64)
                        .sum::<f64>())
                .cos()
            })
            .collect();
        let mut filtered = vec![0.0; size];
        p.apply(&r, &mut filtered).unwrap();
        for (a, r) in filtered.iter().zip(r) {
            assert!((a - factor * (r - 2.0)).abs() < 2e-14);
        }
        assert!(filtered.iter().sum::<f64>().abs() < 2e-12);
    }
}
#[test]
fn validation_restart_and_diagnostics() {
    let shape = [5, 1, 1];
    let b = cell(1.0);
    for q in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(Kerker::new(shape, b, q, Dft::new(shape)).is_err());
    }
    assert!(Kerker::new([0, 1, 1], b, 1.0, Dft::new(shape)).is_err());
    assert!(Kerker::new(shape, [[0.0; 3]; 3], 1.0, Dft::new(shape)).is_err());
    let create = || {
        PreconditionedPulay::new(
            0.3,
            3,
            5,
            Kerker::new(shape, b, 1.0, Dft::new(shape)).unwrap(),
        )
        .unwrap()
    };
    let mut p = create();
    let x = [1.2, 0.9, 0.9, 1.0, 1.0];
    let y = [1.1, 0.95, 0.95, 1.0, 1.0];
    let first = p.mix(&x, &y, 1.0).unwrap();
    assert!(p.diagnostics().contraction.is_none());
    p.mix(&x, &y, 1.0).unwrap();
    assert_eq!(p.diagnostics().contraction, Some(1.0));
    assert_eq!(p.diagnostics().history_used, 1);
    p.reset();
    assert_eq!(p.diagnostics().resets, 1);
    assert_eq!(p.mix(&x, &y, 1.0).unwrap(), first);
    assert!(p.mix(&[f64::NAN; 5], &y, 1.0).is_err());
    assert!(p.mix(&x[..4], &y, 1.0).is_err());
    let z = p.mix(&x, &y, 2.0).unwrap();
    assert_eq!(p.diagnostics().resets, 2);
    assert!((z.iter().sum::<f64>() - 10.0).abs() < 2e-14);
}
// R_G=-(1+1/G^2)*(x_G-target_G), with known inverse P at q0=1.
fn dielectric(x: &[f64], weights: &[f64], shape: [usize; 3]) -> Vec<f64> {
    let mut v: Vec<_> = x.iter().map(|x| Complex64::from(x - 1.0)).collect();
    let mut fft = Dft::new(shape);
    fft.transform(&mut v, false).unwrap();
    for (r, &p) in v.iter_mut().zip(weights) {
        *r = if p == 0.0 {
            Complex64::default()
        } else {
            -*r / p
        };
    }
    fft.transform(&mut v, true).unwrap();
    x.iter().zip(v).map(|(x, r)| x + r.re).collect()
}
fn solve(
    mut mixer: Box<dyn DensityMixer>,
    shape: [usize; 3],
    weights: &[f64],
) -> (usize, f64, f64, Vec<f64>) {
    let size: usize = shape.iter().product();
    let mut x: Vec<_> = (0..size)
        .map(|i| 1.0 + 1e-4 * ((i * 17 + 3) as f64).sin())
        .collect();
    let offset = x.iter().sum::<f64>() / size as f64 - 1.0;
    for v in &mut x {
        *v -= offset;
    }
    let start = Instant::now();
    let mut trace = Vec::new();
    let mut charge: f64 = 0.0;
    for n in 1..=100 {
        let y = dielectric(&x, weights, shape);
        let r = y
            .iter()
            .zip(&x)
            .map(|(y, x)| (y - x).powi(2))
            .sum::<f64>()
            .sqrt();
        trace.push(r);
        if r < 1e-10 {
            return (n, start.elapsed().as_secs_f64(), charge, trace);
        }
        if !r.is_finite() || r > 1e4 {
            return (101, start.elapsed().as_secs_f64(), charge, trace);
        }
        x = mixer.mix(&x, &y, 1.0).unwrap();
        charge = charge.max((x.iter().sum::<f64>() - size as f64).abs());
    }
    (101, start.elapsed().as_secs_f64(), charge, trace)
}
#[test]
fn analytic_sloshing_comparison() {
    let shape = [5, 3, 5];
    let size = shape.iter().product();
    for stretch in [1.0, 4.0, 12.0] {
        let b = cell(stretch);
        // Independent analytic response on an odd grid, not production weights.
        let exact: Vec<f64> = (0..size)
            .map(|i| {
                let c = coord(i, shape);
                let m: [f64; 3] = std::array::from_fn(|a| {
                    if c[a] <= shape[a] / 2 {
                        c[a] as f64
                    } else {
                        c[a] as f64 - shape[a] as f64
                    }
                });
                let g2: f64 = b
                    .iter()
                    .map(|row| row.iter().zip(m).map(|(v, n)| v * n).sum::<f64>().powi(2))
                    .sum();
                g2 / (g2 + 1.0)
            })
            .collect();
        let base = solve(Box::new(PulayMixer::default()), shape, &exact);
        println!(
            "stretch={stretch} default iterations={} seconds={:.6} charge={:.3e} residuals={:?}",
            base.0, base.1, base.2, base.3
        );
        for q0 in [0.5, 1.0, 2.0] {
            let p = Kerker::new(shape, b, q0, Dft::new(shape)).unwrap();
            let r = solve(
                Box::new(PreconditionedPulay::new(0.3, 6, size, p).unwrap()),
                shape,
                &exact,
            );
            println!(
                "stretch={stretch} q0={q0} iterations={} seconds={:.6} charge={:.3e} residuals={:?}",
                r.0, r.1, r.2, r.3
            );
            assert!(r.0 <= 100);
            assert!(r.2 < 2e-12);
            if q0 == 1.0 {
                assert!(r.0 <= 4);
                assert!(r.0 < base.0);
            }
        }
    }
}

#[test]
fn nearly_stationary_map_and_transform_failure_are_explicit() {
    let shape = [5, 1, 1];
    for beta in [0.1, 0.3, 0.8] {
        let filter = Kerker::new(shape, cell(2.0), 0.5, Dft::new(shape)).unwrap();
        let mut p = PreconditionedPulay::new(beta, 6, 5, filter).unwrap();
        let mut x = vec![1.1, 0.9, 1.0, 1.0, 1.0];
        for _ in 0..20 {
            let y: Vec<_> = x.iter().map(|v| 1.0 + 0.999 * (v - 1.0)).collect();
            x = p.mix(&x, &y, 1.0).unwrap();
            assert!((x.iter().sum::<f64>() - 5.0).abs() < 2e-14);
        }
        assert!(x.iter().all(|v| (v - 1.0).abs() < 2e-8));
    }
    let fail = |_: &mut [Complex64], _: bool| {
        Err(SolverError::Application(
            "injected transform failure".into(),
        ))
    };
    let filter = Kerker::new(shape, cell(1.0), 1.0, fail).unwrap();
    let mut p = PreconditionedPulay::new(0.3, 6, 5, filter).unwrap();
    assert!(matches!(
        p.mix(&[1.0; 5], &[1.1; 5], 1.0),
        Err(SolverError::Application(_))
    ));
    assert_eq!(p.diagnostics().history_used, 0);
}

#[test]
fn infeasible_first_filtered_step_recovers_without_clipping() {
    let shape = [5, 1, 1];
    let mut p = PreconditionedPulay::new(
        0.3,
        6,
        5,
        Kerker::new(shape, cell(1.0), 1.0, Dft::new(shape)).unwrap(),
    )
    .unwrap();
    let x = [0.0, 1.0, 1.0, 1.0, 2.0];
    let y = [0.0, 1.2, 0.8, 1.0, 2.0];
    let expected = mixing::LinearMixer { beta: 0.3 }.mix(&x, &y, 1.0).unwrap();
    let actual = p.mix(&x, &y, 1.0).unwrap();
    assert_eq!(actual, expected);
    assert!(p.diagnostics().used_raw_recovery);
    assert_eq!(p.diagnostics().raw_recoveries, 1);
    assert_eq!(p.diagnostics().history_used, 0);
    p.reset();
    assert_eq!(p.diagnostics().raw_recoveries, 1);
    assert!(!p.diagnostics().used_raw_recovery);
}
