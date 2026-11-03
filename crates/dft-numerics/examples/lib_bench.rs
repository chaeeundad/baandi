//! Micro-benchmark of the numerical kernels used on the hot path, for comparison with
//! vendor libraries (tools/lib_bench/mkl_bench.c prints the same table).
//! Run: cargo run --release -p dft-numerics --example lib_bench
use dft_numerics::linalg::{ComplexMatrix, MatrixOp, gemm, hermitian_eigh};
use dft_numerics::{Fft3d, FftPlan3d};
use ndarray::Array3;
use num_complex::Complex64;
use std::time::Instant;

struct Lcg(u64);
impl Lcg {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        ((self.0 >> 11) as f64) / ((1u64 << 53) as f64) - 0.5
    }
}

fn main() {
    let mut rng = Lcg(12345);
    println!(
        "kind size microseconds_per_op   (op = one complex 3D transform, forward or inverse, in place cost; one gemm; one eigendecomposition)"
    );
    for n in [24usize, 25, 27, 36, 48] {
        let mut plan = FftPlan3d::new([n, n, n]).unwrap();
        let mut a = Array3::<Complex64>::zeros((n, n, n));
        for z in a.iter_mut() {
            *z = Complex64::new(rng.next(), rng.next());
        }
        let mut b = a.clone();
        let reps = (3_000_000 / (n * n * n)).max(20);
        for _ in 0..3 {
            plan.to_real(&a, &mut b).unwrap();
            plan.to_reciprocal(&b, &mut a).unwrap();
        }
        let t = Instant::now();
        for _ in 0..reps {
            plan.to_real(&a, &mut b).unwrap();
            plan.to_reciprocal(&b, &mut a).unwrap();
        }
        let us = t.elapsed().as_secs_f64() * 1e6 / (2 * reps) as f64;
        println!("fft3d {n}^3 {us:.2}");
    }
    let npw = 750usize;
    for m in [8usize, 16, 32, 48] {
        let mut a = ComplexMatrix::zeros(npw, m).unwrap();
        let mut b = ComplexMatrix::zeros(npw, 8).unwrap();
        for j in 0..m {
            for z in a.column_mut(j) {
                *z = Complex64::new(rng.next(), rng.next());
            }
        }
        for j in 0..8 {
            for z in b.column_mut(j) {
                *z = Complex64::new(rng.next(), rng.next());
            }
        }
        let reps = 2000;
        for _ in 0..10 {
            let _ = gemm(&a, MatrixOp::Adjoint, &b, MatrixOp::Identity).unwrap();
        }
        let t = Instant::now();
        for _ in 0..reps {
            let _ = gemm(&a, MatrixOp::Adjoint, &b, MatrixOp::Identity).unwrap();
        }
        println!(
            "zgemm_CN {npw}x{m}^H * {npw}x8 {:.2}",
            t.elapsed().as_secs_f64() * 1e6 / reps as f64
        );
    }
    for m in [16usize, 32, 48] {
        let mut h = ComplexMatrix::zeros(m, m).unwrap();
        for i in 0..m {
            for j in 0..=i {
                let z = Complex64::new(rng.next(), if i == j { 0.0 } else { rng.next() });
                h[(i, j)] = z;
                h[(j, i)] = z.conj();
            }
        }
        let reps = 300;
        for _ in 0..3 {
            let _ = hermitian_eigh(&h).unwrap();
        }
        let t = Instant::now();
        for _ in 0..reps {
            let _ = hermitian_eigh(&h).unwrap();
        }
        println!(
            "zheev {m}x{m} {:.2}",
            t.elapsed().as_secs_f64() * 1e6 / reps as f64
        );
    }
}
