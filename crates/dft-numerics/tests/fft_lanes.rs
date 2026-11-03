//! Rejected packed-lane experiment; never selected by the production backend.
use dft_numerics::{Fft3d, FftPlan3d};
use ndarray::Array3;
use num_complex::Complex64 as C;
use rustfft::{Fft, FftPlanner};
use std::{sync::Arc, time::Instant};
struct PackedLanes {
    n: usize,
    packed: Vec<C>,
    scratch: Vec<Vec<C>>,
    plan: Arc<dyn Fft<f64>>,
}
impl PackedLanes {
    fn new(n: usize, threads: usize) -> Self {
        let plan = FftPlanner::new().plan_fft_inverse(n);
        Self {
            n,
            packed: vec![C::default(); n * n * n],
            scratch: vec![vec![C::default(); plan.get_inplace_scratch_len()]; threads],
            plan,
        }
    }
    fn inverse(&mut self, input: &Array3<C>, output: &mut Array3<C>) {
        let n = self.n;
        let data = output.as_slice_mut().unwrap();
        data.copy_from_slice(input.as_slice().unwrap());
        for axis in 0..3 {
            let stride = n.pow((2 - axis) as u32);
            for (lane, values) in self.packed.chunks_exact_mut(n).enumerate() {
                let start = (lane / stride) * n * stride + lane % stride;
                for (i, z) in values.iter_mut().enumerate() {
                    *z = data[start + i * stride];
                }
            }
            let chunk = (n * n).div_ceil(self.scratch.len()) * n;
            let plan = &self.plan;
            std::thread::scope(|scope| {
                for (values, scratch) in self.packed.chunks_mut(chunk).zip(&mut self.scratch) {
                    scope.spawn(move || {
                        for lane in values.chunks_exact_mut(n) {
                            plan.process_with_scratch(lane, scratch);
                            for z in lane {
                                *z *= 1.0 / n as f64;
                            }
                        }
                    });
                }
            });
            for (lane, values) in self.packed.chunks_exact(n).enumerate() {
                let start = (lane / stride) * n * stride + lane % stride;
                for (i, &z) in values.iter().enumerate() {
                    data[start + i * stride] = z;
                }
            }
        }
        for z in data {
            *z *= (n * n * n) as f64;
        }
    }
}
fn median(mut samples: Vec<f64>) -> f64 {
    samples.sort_by(f64::total_cmp);
    samples[samples.len() / 2]
}
#[test]
#[ignore = "release measurement: cargo test --release -p dft-numerics --test fft_lanes -- --ignored --nocapture"]
fn packed_lane_experiment() {
    // 36^3 is the conventional Si 30 Ry grid; 48^3 is the requested estimate.
    for n in [25, 36, 48] {
        let input = Array3::from_shape_fn([n; 3], |(i, j, k)| {
            C::new(
                ((i + 3 * j + 7 * k) as f64).sin(),
                ((i + 2 * j + k) as f64).cos(),
            )
        });
        let mut output = Array3::zeros([n; 3]);
        let mut reference = output.clone();
        let mut serial = FftPlan3d::new([n; 3]).unwrap();
        serial.to_real(&input, &mut reference).unwrap();
        let baseline = median(
            (0..11)
                .map(|_| {
                    let start = Instant::now();
                    for _ in 0..10 {
                        serial.to_real(&input, &mut output).unwrap();
                    }
                    start.elapsed().as_secs_f64() / 10.0
                })
                .collect(),
        );
        for threads in [2, 4, 8] {
            let mut experiment = PackedLanes::new(n, threads);
            experiment.inverse(&input, &mut output);
            assert!(
                output
                    .iter()
                    .zip(&reference)
                    .all(|(a, b)| a.re.to_bits() == b.re.to_bits()
                        && a.im.to_bits() == b.im.to_bits())
            );
            let timing = median(
                (0..11)
                    .map(|_| {
                        let start = Instant::now();
                        for _ in 0..10 {
                            experiment.inverse(&input, &mut output);
                        }
                        start.elapsed().as_secs_f64() / 10.0
                    })
                    .collect(),
            );
            println!(
                "grid={n}^3 threads={threads} serial_s={baseline:.9} lanes_s={timing:.9} speedup={:.3} bit_identical=true",
                baseline / timing
            );
        }
    }
}
