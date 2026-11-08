//! P47B: release-only, single-thread rustfft/faer precision experiment.
//! First run dft-engine --test fp32_snapshots. Then cargo test --release
//! -p dft-numerics --test fp32_feasibility -- --ignored --nocapture --test-threads=1.
//! Repeat with RUSTFLAGS='-C target-cpu=native' and a separate CARGO_TARGET_DIR.
use num_complex::Complex;
use rustfft::{Fft, FftNum, FftPlanner};
use std::{hint::black_box, path::PathBuf, sync::Arc, time::Instant};

struct Lanes<T: FftNum> {
    dims: [usize; 3],
    forward: [Arc<dyn Fft<T>>; 3],
    inverse: [Arc<dyn Fft<T>>; 3],
    grid: Vec<Complex<T>>,
    lane: Vec<Complex<T>>,
    scratch: Vec<Complex<T>>,
}
impl<T: FftNum + Default> Lanes<T> {
    fn new(dims: [usize; 3]) -> Self {
        let mut p = FftPlanner::<T>::new();
        let forward = dims.map(|n| p.plan_fft_forward(n));
        let inverse = dims.map(|n| p.plan_fft_inverse(n));
        let max = dims.into_iter().max().unwrap();
        let scratch = forward
            .iter()
            .chain(&inverse)
            .map(|p| p.get_inplace_scratch_len())
            .max()
            .unwrap();
        Self {
            dims,
            forward,
            inverse,
            grid: vec![Complex::default(); dims.iter().product()],
            lane: vec![Complex::default(); max],
            scratch: vec![Complex::default(); scratch],
        }
    }
    fn transform(&mut self, inverse: bool) {
        let plans = if inverse {
            &self.inverse
        } else {
            &self.forward
        };
        let [nx, ny, nz] = self.dims;
        for lane in self.grid.chunks_exact_mut(nz) {
            plans[2].process_with_scratch(lane, &mut self.scratch);
        }
        for x in 0..nx {
            for z in 0..nz {
                for y in 0..ny {
                    self.lane[y] = self.grid[(x * ny + y) * nz + z];
                }
                plans[1].process_with_scratch(&mut self.lane[..ny], &mut self.scratch);
                for y in 0..ny {
                    self.grid[(x * ny + y) * nz + z] = self.lane[y];
                }
            }
        }
        for y in 0..ny {
            for z in 0..nz {
                for x in 0..nx {
                    self.lane[x] = self.grid[(x * ny + y) * nz + z];
                }
                plans[0].process_with_scratch(&mut self.lane[..nx], &mut self.scratch);
                for x in 0..nx {
                    self.grid[(x * ny + y) * nz + z] = self.lane[x];
                }
            }
        }
    }
    // FftNum does not require the assignment-operator bounds used by Complex<T>.
    #[allow(clippy::assign_op_pattern)]
    fn cycle(&mut self) {
        self.transform(true);
        let scale = T::from_f64(1. / self.grid.len() as f64).unwrap();
        for v in &mut self.grid {
            *v = *v * scale;
        }
        self.transform(false);
        black_box(&self.grid);
    }
}
#[derive(Clone)]
struct Snapshot {
    dims: [usize; 3],
    n: usize,
    b: usize,
    p: usize,
    plus: Vec<usize>,
    minus: Vec<usize>,
    v: Vec<f64>,
    kin: Vec<f64>,
    x: Vec<f64>,
    proj: Vec<f64>,
    d: Vec<f64>,
    reference: Vec<f64>,
}
fn read(name: &str) -> Snapshot {
    let path =
        PathBuf::from(env!("CARGO_MANIFEST_DIR")).join(format!("../../P47B-evidence/{name}.bin"));
    let bytes = std::fs::read(path).unwrap();
    let (words, remainder) = bytes.as_chunks::<8>();
    assert!(remainder.is_empty());
    let mut chunks = words.iter();
    let mut uint = || u64::from_le_bytes(*chunks.next().unwrap()) as usize;
    let dims = [uint(), uint(), uint()];
    let n = uint();
    let b = uint();
    let p = uint();
    let plus = (0..n.div_ceil(2)).map(|_| uint()).collect();
    let minus = (0..n.div_ceil(2)).map(|_| uint()).collect();
    let mut vals = |n| {
        (0..n)
            .map(|_| f64::from_le_bytes(*chunks.next().unwrap()))
            .collect::<Vec<_>>()
    };
    let result = Snapshot {
        dims,
        n,
        b,
        p,
        plus,
        minus,
        v: vals(dims.iter().product()),
        kin: vals(n),
        x: vals(n * b),
        proj: vals(n * p),
        d: vals(p * p),
        reference: vals(n * b),
    };
    assert!(chunks.next().is_none());
    result
}
// Each expansion uses the library's native scalar type, never f64 arithmetic in the f32 kernel.
macro_rules! precision {
    ($mod:ident,$t:ty) => {
        mod $mod {
            use super::*;
            use faer::{Accum, Mat, Par, linalg::matmul::matmul};
            pub struct Kernel {
                fft: Lanes<$t>,
                x: Mat<$t>,
                p: Mat<$t>,
                d: Mat<$t>,
                q: Mat<$t>,
                w: Mat<$t>,
                nl: Mat<$t>,
                y: Vec<$t>,
                v: Vec<$t>,
                kin: Vec<$t>,
                plus: Vec<usize>,
                minus: Vec<usize>,
                n: usize,
                b: usize,
            }
            impl Kernel {
                pub fn new(s: &Snapshot) -> Self {
                    Self {
                        fft: Lanes::new(s.dims),
                        x: Mat::from_fn(s.n, s.b, |i, j| s.x[j * s.n + i] as $t),
                        p: Mat::from_fn(s.n, s.p, |i, j| s.proj[j * s.n + i] as $t),
                        d: Mat::from_fn(s.p, s.p, |i, j| s.d[j * s.p + i] as $t),
                        q: Mat::zeros(s.p, s.b),
                        w: Mat::zeros(s.p, s.b),
                        nl: Mat::zeros(s.n, s.b),
                        y: vec![0.; s.n * s.b],
                        v: s.v.iter().map(|&x| x as $t).collect(),
                        kin: s.kin.iter().map(|&x| x as $t).collect(),
                        plus: s.plus.clone(),
                        minus: s.minus.clone(),
                        n: s.n,
                        b: s.b,
                    }
                }
                pub fn nonlocal(&mut self) {
                    matmul(
                        self.q.as_mut(),
                        Accum::Replace,
                        self.p.as_ref().transpose(),
                        self.x.as_ref(),
                        1.,
                        Par::Seq,
                    );
                    matmul(
                        self.w.as_mut(),
                        Accum::Replace,
                        self.d.as_ref(),
                        self.q.as_ref(),
                        1.,
                        Par::Seq,
                    );
                    matmul(
                        self.nl.as_mut(),
                        Accum::Replace,
                        self.p.as_ref(),
                        self.w.as_ref(),
                        1.,
                        Par::Seq,
                    );
                    black_box(&self.nl);
                }
                pub fn local(&mut self) {
                    let rt = (2. as $t).sqrt();
                    let scale = 1. / self.fft.grid.len() as $t;
                    for j in (0..self.b).step_by(2) {
                        let get = |i, k| if k < self.b { self.x[(i, k)] } else { 0. };
                        self.fft.grid.fill(Complex::default());
                        self.fft.grid[0] = Complex::new(get(0, j), get(0, j + 1));
                        for i in 1..self.plus.len() {
                            let a = Complex::new(get(2 * i - 1, j), get(2 * i, j)) / rt;
                            let b = Complex::new(get(2 * i - 1, j + 1), get(2 * i, j + 1)) / rt;
                            self.fft.grid[self.plus[i]] = a + Complex::new(-b.im, b.re);
                            self.fft.grid[self.minus[i]] = a.conj() + Complex::new(b.im, b.re);
                        }
                        self.fft.transform(true);
                        for (z, &v) in self.fft.grid.iter_mut().zip(&self.v) {
                            *z *= v * scale;
                        }
                        self.fft.transform(false);
                        self.y[j * self.n] = self.fft.grid[0].re;
                        if j + 1 < self.b {
                            self.y[(j + 1) * self.n] = self.fft.grid[0].im;
                        }
                        for i in 1..self.plus.len() {
                            let z = self.fft.grid[self.plus[i]];
                            let m = self.fft.grid[self.minus[i]].conj();
                            let a = (z + m) / rt;
                            let b = (z - m) / rt;
                            self.y[j * self.n + 2 * i - 1] = a.re;
                            self.y[j * self.n + 2 * i] = a.im;
                            if j + 1 < self.b {
                                self.y[(j + 1) * self.n + 2 * i - 1] = b.im;
                                self.y[(j + 1) * self.n + 2 * i] = -b.re;
                            }
                        }
                    }
                    black_box(&self.y);
                }
                pub fn apply(&mut self) {
                    self.local();
                    self.nonlocal();
                    for j in 0..self.b {
                        for i in 0..self.n {
                            self.y[j * self.n + i] +=
                                self.kin[i] * self.x[(i, j)] + self.nl[(i, j)];
                        }
                    }
                    black_box(&self.y);
                }
                pub fn result(&self) -> Vec<f64> {
                    self.y.iter().map(|&x| x as f64).collect()
                }
            }
            pub struct Product {
                a: Mat<$t>,
                b: Mat<$t>,
                c: Mat<$t>,
                transpose: bool,
            }
            impl Product {
                pub fn new(m: usize, n: usize, k: usize, transpose: bool) -> Self {
                    let a = if transpose {
                        Mat::from_fn(k, m, |i, j| ((i * 7 + j * 13) as f64).sin() as $t)
                    } else {
                        Mat::from_fn(m, k, |i, j| ((i * 7 + j * 13) as f64).sin() as $t)
                    };
                    Self {
                        a,
                        b: Mat::from_fn(k, n, |i, j| ((i * 11 + j * 3) as f64).cos() as $t),
                        c: Mat::zeros(m, n),
                        transpose,
                    }
                }
                pub fn run(&mut self) {
                    let a = if self.transpose {
                        self.a.as_ref().transpose()
                    } else {
                        self.a.as_ref()
                    };
                    matmul(
                        self.c.as_mut(),
                        Accum::Replace,
                        a,
                        self.b.as_ref(),
                        1.,
                        Par::Seq,
                    );
                    black_box(&self.c);
                }
            }
        }
    };
}
precision!(single, f32);
precision!(double, f64);
fn paired(label: &str, repeats: usize, mut a: impl FnMut(), mut b: impl FnMut()) {
    a();
    b();
    for pair in 0..30 {
        let mut measure = |single: bool| {
            let t = Instant::now();
            for _ in 0..repeats {
                if single { a() } else { b() }
            }
            t.elapsed().as_secs_f64() / repeats as f64
        };
        let (t32, t64) = if pair % 2 == 0 {
            (measure(true), measure(false))
        } else {
            let t64 = measure(false);
            (measure(true), t64)
        };
        println!("TIME,{label},{pair},{t32:.12},{t64:.12}");
    }
}
#[test]
#[ignore = "release-only precision timings, 30 alternating pairs, external snapshots"]
fn measure() {
    if cfg!(debug_assertions) {
        panic!("run precision measurements in release mode");
    }
    for dims in [
        [25, 25, 25],
        [32, 32, 32],
        [36, 36, 36],
        [36, 36, 72],
        [48, 48, 48],
        [72, 36, 36],
        [72, 72, 36],
    ] {
        let mut a = Lanes::<f32>::new(dims);
        let mut b = Lanes::<f64>::new(dims);
        for (i, (a, b)) in a.grid.iter_mut().zip(&mut b.grid).enumerate() {
            let x = (i as f64).sin();
            *a = Complex::new(x as f32, 0.);
            *b = Complex::new(x, 0.);
        }
        paired(
            &format!("fft-{}x{}x{}", dims[0], dims[1], dims[2]),
            5,
            || a.cycle(),
            || b.cycle(),
        );
    }
    for name in [
        "si8-early",
        "si8-converged",
        "si16-early",
        "si16-converged",
        "si32-early",
        "si32-converged",
    ] {
        let s = read(name);
        let mut a = single::Kernel::new(&s);
        let mut b = double::Kernel::new(&s);
        a.apply();
        b.apply();
        let ya = a.result();
        let yb = b.result();
        for j in 0..s.b {
            let mut norm = 0.;
            let mut err = 0.;
            let mut oracle = 0.;
            let mut ray = 0.;
            for i in 0..s.n {
                let k = j * s.n + i;
                norm += yb[k] * yb[k];
                err += (ya[k] - yb[k]).powi(2);
                oracle += (yb[k] - s.reference[k]).powi(2);
                ray += s.x[k] * (ya[k] - yb[k]);
            }
            println!(
                "ERROR,{name},{j},{:.17e},{:.17e},{:.17e},{:.17e},{:.17e}",
                err.sqrt(),
                (err / norm).sqrt(),
                ray.abs(),
                (err - ray * ray).max(0.).sqrt(),
                oracle.sqrt()
            );
            assert!(
                oracle.sqrt() < 1e-11,
                "f64 experiment must agree with production H in Ha"
            );
        }
        paired(&format!("h-{name}"), 2, || a.apply(), || b.apply());
        paired(&format!("local-{name}"), 2, || a.local(), || b.local());
        paired(
            &format!("nonlocal-{name}"),
            5,
            || a.nonlocal(),
            || b.nonlocal(),
        );
        for width in [1, 2, 4, 8] {
            let mut narrow = s.clone();
            narrow.b = width;
            narrow.x.truncate(s.n * width);
            narrow.reference.truncate(s.n * width);
            let mut a = single::Kernel::new(&narrow);
            let mut b = double::Kernel::new(&narrow);
            paired(
                &format!("h-{name}-width{width}"),
                2,
                || a.apply(),
                || b.apply(),
            );
        }
    }
    let s8 = read("si8-early");
    let s16 = read("si16-early");
    // Si32 dimension is supplied by the engine basis count, not the FFT box volume.
    let n32 = std::env::var("P47B_SI32_N")
        .expect("set actual Si32 plane-wave count")
        .parse()
        .unwrap();
    for (name, n, b, p) in [
        ("si8", s8.n, s8.b, s8.p),
        ("si16", s16.n, s16.b, s16.p),
        ("si32", n32, 68, 128),
    ] {
        for width in [2, b] {
            for q in [b, 2 * b, 6 * b] {
                for (site, m, k, cols, tr) in [
                    ("gram", q, n, width, true),
                    ("rotation", n, q, b, false),
                    ("projection", n, q, width, false),
                ] {
                    let mut a = single::Product::new(m, cols, k, tr);
                    let mut d = double::Product::new(m, cols, k, tr);
                    paired(
                        &format!("gemm-{name}-{site}-{m}x{cols}x{k}"),
                        3,
                        || a.run(),
                        || d.run(),
                    );
                }
            }
        }
        for width in [2, b] {
            for (site, m, k, cols, tr) in [
                ("project", p, n, width, true),
                ("reconstruct", n, p, width, false),
            ] {
                let mut a = single::Product::new(m, cols, k, tr);
                let mut d = double::Product::new(m, cols, k, tr);
                paired(
                    &format!("nl-{name}-{site}-{m}x{cols}x{k}"),
                    5,
                    || a.run(),
                    || d.run(),
                );
            }
        }
    }
}

#[test]
#[ignore = "release-only projected Ritz accuracy on external physical snapshots"]
fn projected_ritz_accuracy() {
    use dft_numerics::real::{RealMatrix, symmetric_eigh};
    for name in [
        "si8-early",
        "si8-converged",
        "si16-early",
        "si16-converged",
        "si32-early",
        "si32-converged",
    ] {
        let s = read(name);
        let mut a = single::Kernel::new(&s);
        let mut b = double::Kernel::new(&s);
        a.apply();
        b.apply();
        let y32 = a.result();
        let y64 = b.result();
        let mut q32 = RealMatrix::zeros(s.b, s.b).unwrap();
        let mut q64 = RealMatrix::zeros(s.b, s.b).unwrap();
        for j in 0..s.b {
            for i in 0..s.b {
                q32[(i, j)] = (0..s.n).map(|k| s.x[i * s.n + k] * y32[j * s.n + k]).sum();
                q64[(i, j)] = (0..s.n).map(|k| s.x[i * s.n + k] * y64[j * s.n + k]).sum();
            }
        }
        let mut defect = 0.0_f64;
        for j in 0..s.b {
            for i in 0..j {
                defect = defect.max((q32[(i, j)] - q32[(j, i)]).abs());
                let v = (q32[(i, j)] + q32[(j, i)]) / 2.;
                q32[(i, j)] = v;
                q32[(j, i)] = v;
                let v = (q64[(i, j)] + q64[(j, i)]) / 2.;
                q64[(i, j)] = v;
                q64[(j, i)] = v;
            }
        }
        let (e32, _) = symmetric_eigh(&q32).unwrap();
        let (e64, _) = symmetric_eigh(&q64).unwrap();
        let max_error = e32
            .iter()
            .zip(&e64)
            .map(|(a, b)| (a - b).abs())
            .fold(0., f64::max);
        println!("RITZ,{name},{max_error:.17e},{defect:.17e}");
        for tau in [1e-4, 1e-5, 1e-6, 1e-7, 1e-8] {
            println!(
                "TOLERANCE,{name},{tau:.1e},{},{}",
                e32.iter()
                    .zip(&e64)
                    .filter(|(a, b)| (*a - *b).abs() > tau)
                    .count(),
                s.b
            );
        }
    }
}
