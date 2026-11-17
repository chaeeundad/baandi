//! REDUCED-EIGEN-SLICES (dev08 microbenchmark, MKL, one core): can the
//! replicated reduced eigenproblem of the MPI lean Davidson be shortened by
//! keeping the tridiagonalization replicated and splitting only the wanted
//! eigenvectors over the ranks?
//!
//! Per order m (k = ceil(m * frac) wanted lowest roots, P ranks):
//! - `dsyevd`: the production solve (all pairs, `dft_mkl::SymmetricEigen`).
//! - `trd`: `dsytrd` alone (the part that stays replicated).
//! - `vals`: `dstebz` for the k lowest eigenvalues of the tridiagonal (replicated).
//! - `stein(P)` / `mrrr(P)`: the largest rank's slice of k / P eigenvectors by
//!   inverse iteration (`dstein`) or MRRR (`dstemr`, range by index), plus
//!   `dormtr` on those columns.
//! - `sliced(P)` = trd + vals + max slice (+ the allgather of m x k doubles,
//!   NOT included: take it from the dev08 collective timers).
//!
//! Accuracy of the assembled slices (computed independently, slice bounds moved
//! to the widest eigenvalue gap within +-2 of the even split): max |Z^T Z - I|
//! and max_j |A z_j - w_j z_j| / |A|_max, against the same for dsyevd. The lean
//! solver's return Gram gate is 1e-12; a variant above it is not usable.
//!
//! Matrix models: `dense` (diagonal ramp plus a dense symmetric random part)
//! and `ritz` (Davidson-like: diagonal Ritz block with a dense border of
//! m/3 columns), both scaled like a Hamiltonian in Hartree.
//!
//!   MKLROOT=/opt/intel/oneapi/mkl/2026.1 cargo run --release -p dft-mkl \
//!       --features mkl --example evd_slices -- [frac=0.34] [reps=20]
//!   (bind to one core: taskset -c 2; optionally run 8 copies at once with
//!   taskset -c 0..7 to see the shared-node effect)
use std::ffi::{c_char, c_int};
use std::time::Instant;

unsafe extern "C" {
    fn dsytrd_(
        uplo: *const c_char,
        n: *const c_int,
        a: *mut f64,
        lda: *const c_int,
        d: *mut f64,
        e: *mut f64,
        tau: *mut f64,
        work: *mut f64,
        lwork: *const c_int,
        info: *mut c_int,
    );
    fn dormtr_(
        side: *const c_char,
        uplo: *const c_char,
        trans: *const c_char,
        m: *const c_int,
        n: *const c_int,
        a: *const f64,
        lda: *const c_int,
        tau: *const f64,
        c: *mut f64,
        ldc: *const c_int,
        work: *mut f64,
        lwork: *const c_int,
        info: *mut c_int,
    );
    fn dstebz_(
        range: *const c_char,
        order: *const c_char,
        n: *const c_int,
        vl: *const f64,
        vu: *const f64,
        il: *const c_int,
        iu: *const c_int,
        abstol: *const f64,
        d: *const f64,
        e: *const f64,
        m: *mut c_int,
        nsplit: *mut c_int,
        w: *mut f64,
        iblock: *mut c_int,
        isplit: *mut c_int,
        work: *mut f64,
        iwork: *mut c_int,
        info: *mut c_int,
    );
    fn dstein_(
        n: *const c_int,
        d: *const f64,
        e: *const f64,
        m: *const c_int,
        w: *const f64,
        iblock: *const c_int,
        isplit: *const c_int,
        z: *mut f64,
        ldz: *const c_int,
        work: *mut f64,
        iwork: *mut c_int,
        ifail: *mut c_int,
        info: *mut c_int,
    );
    fn dstemr_(
        jobz: *const c_char,
        range: *const c_char,
        n: *const c_int,
        d: *mut f64,
        e: *mut f64,
        vl: *const f64,
        vu: *const f64,
        il: *const c_int,
        iu: *const c_int,
        m: *mut c_int,
        w: *mut f64,
        z: *mut f64,
        ldz: *const c_int,
        nzc: *const c_int,
        isuppz: *mut c_int,
        tryrac: *mut c_int,
        work: *mut f64,
        lwork: *const c_int,
        iwork: *mut c_int,
        liwork: *const c_int,
        info: *mut c_int,
    );
}

fn ch(c: u8) -> c_char {
    c as c_char
}

struct Rng(u64);
impl Rng {
    fn next(&mut self) -> f64 {
        self.0 = self
            .0
            .wrapping_mul(6364136223846793005)
            .wrapping_add(1442695040888963407);
        (self.0 >> 11) as f64 / (1u64 << 53) as f64 - 0.5
    }
}

/// Column-major symmetric test matrix (both triangles filled).
fn matrix(m: usize, model: &str, seed: u64) -> Vec<f64> {
    let mut r = Rng(seed);
    let mut a = vec![0.0; m * m];
    let border = m - m / 3;
    for j in 0..m {
        for i in 0..=j {
            let v = match model {
                "dense" => {
                    let d = if i == j {
                        -0.5 + 20.0 * j as f64 / m as f64
                    } else {
                        0.0
                    };
                    d + 0.2 * r.next()
                }
                _ => {
                    // Ritz block: diagonal ascending, tiny couplings; border: dense.
                    if i == j {
                        -0.5 + 20.0 * (j as f64 / m as f64).powi(2) + 0.01 * r.next()
                    } else if j >= border {
                        0.05 * r.next()
                    } else {
                        1e-9 * r.next()
                    }
                }
            };
            a[i + j * m] = v;
            a[j + i * m] = v;
        }
    }
    a
}

/// Tridiagonal reduction result: Householder data, d, e, tau.
struct Trd {
    a: Vec<f64>,
    d: Vec<f64>,
    e: Vec<f64>,
    tau: Vec<f64>,
}

fn sytrd(m: usize, a0: &[f64], work: &mut Vec<f64>) -> Trd {
    let n = m as c_int;
    let mut t = Trd {
        a: a0.to_vec(),
        d: vec![0.0; m],
        e: vec![0.0; m.max(2) - 1],
        tau: vec![0.0; m.max(2) - 1],
    };
    let mut info = 0;
    if work.is_empty() {
        let (mut q, query) = (0.0, -1 as c_int);
        // SAFETY: workspace query (lwork = -1): only `q` and `info` are written;
        // the other buffers have the documented lengths for order m.
        unsafe {
            dsytrd_(
                &ch(b'L'),
                &n,
                t.a.as_mut_ptr(),
                &n,
                t.d.as_mut_ptr(),
                t.e.as_mut_ptr(),
                t.tau.as_mut_ptr(),
                &mut q,
                &query,
                &mut info,
            );
        }
        assert_eq!(info, 0);
        work.resize(q as usize, 0.0);
    }
    let lwork = work.len() as c_int;
    // SAFETY: A is m x m (lda = m), d has m, e and tau m - 1 elements, work has
    // the queried length; all buffers are uniquely borrowed for the call.
    unsafe {
        dsytrd_(
            &ch(b'L'),
            &n,
            t.a.as_mut_ptr(),
            &n,
            t.d.as_mut_ptr(),
            t.e.as_mut_ptr(),
            t.tau.as_mut_ptr(),
            work.as_mut_ptr(),
            &lwork,
            &mut info,
        );
    }
    assert_eq!(info, 0, "dsytrd");
    t
}

/// The k lowest eigenvalues of the tridiagonal with their block data.
struct Vals {
    w: Vec<f64>,
    iblock: Vec<c_int>,
    isplit: Vec<c_int>,
}

fn stebz(t: &Trd, k: usize) -> Vals {
    let m = t.d.len();
    let n = m as c_int;
    let (il, iu) = (1 as c_int, k as c_int);
    let (mut found, mut nsplit, mut info) = (0, 0, 0);
    let mut v = Vals {
        w: vec![0.0; m],
        iblock: vec![0; m],
        isplit: vec![0; m],
    };
    let mut work = vec![0.0; 4 * m];
    let mut iwork = vec![0 as c_int; 3 * m];
    // SAFETY: d has m and e m - 1 elements; w, iblock, isplit have m; work 4m and
    // iwork 3m as documented; vl/vu are not referenced for RANGE = 'I'.
    unsafe {
        dstebz_(
            &ch(b'I'),
            &ch(b'B'),
            &n,
            &0.0,
            &0.0,
            &il,
            &iu,
            &0.0,
            t.d.as_ptr(),
            t.e.as_ptr(),
            &mut found,
            &mut nsplit,
            v.w.as_mut_ptr(),
            v.iblock.as_mut_ptr(),
            v.isplit.as_mut_ptr(),
            work.as_mut_ptr(),
            iwork.as_mut_ptr(),
            &mut info,
        );
    }
    assert_eq!(info, 0, "dstebz");
    assert_eq!(found as usize, k);
    v.w.truncate(k);
    v.iblock.truncate(k);
    v
}

/// Eigenvectors of the tridiagonal for roots `lo..hi` (0-based), m x (hi - lo).
fn stein(t: &Trd, v: &Vals, lo: usize, hi: usize) -> Vec<f64> {
    let m = t.d.len();
    let c = hi - lo;
    let (n, nc) = (m as c_int, c as c_int);
    let mut z = vec![0.0; m * c];
    let mut work = vec![0.0; 5 * m];
    let mut iwork = vec![0 as c_int; m];
    let mut ifail = vec![0 as c_int; c];
    let mut info = 0;
    // SAFETY: d (m), e (m - 1), w/iblock (c entries of the slice, block order
    // as dstebz returned them), isplit (m), z (m x c, ldz = m), work 5m,
    // iwork m, ifail c, as documented.
    unsafe {
        dstein_(
            &n,
            t.d.as_ptr(),
            t.e.as_ptr(),
            &nc,
            v.w[lo..hi].as_ptr(),
            v.iblock[lo..hi].as_ptr(),
            v.isplit.as_ptr(),
            z.as_mut_ptr(),
            &n,
            work.as_mut_ptr(),
            iwork.as_mut_ptr(),
            ifail.as_mut_ptr(),
            &mut info,
        );
    }
    assert_eq!(info, 0, "dstein");
    z
}

/// MRRR eigenpairs of the tridiagonal for roots `lo..hi` (0-based).
fn mrrr(t: &Trd, lo: usize, hi: usize) -> (Vec<f64>, Vec<f64>) {
    let m = t.d.len();
    let c = hi - lo;
    let n = m as c_int;
    let (il, iu) = (lo as c_int + 1, hi as c_int);
    let mut d = t.d.clone();
    let mut e = t.e.clone();
    e.push(0.0); // dstemr needs e of length n (workspace in the last entry)
    let mut w = vec![0.0; m];
    let mut z = vec![0.0; m * c.max(1)];
    let mut isuppz = vec![0 as c_int; 2 * c.max(1)];
    let (mut found, mut tryrac, mut info) = (0, 1 as c_int, 0);
    let nzc = c as c_int;
    let (mut wq, mut iq, query) = (0.0, 0 as c_int, -1 as c_int);
    // SAFETY: workspace query (lwork = liwork = -1): only wq, iq and info are
    // written besides the documented scalar outputs; buffers have their sizes.
    unsafe {
        dstemr_(
            &ch(b'V'),
            &ch(b'I'),
            &n,
            d.as_mut_ptr(),
            e.as_mut_ptr(),
            &0.0,
            &0.0,
            &il,
            &iu,
            &mut found,
            w.as_mut_ptr(),
            z.as_mut_ptr(),
            &n,
            &nzc,
            isuppz.as_mut_ptr(),
            &mut tryrac,
            &mut wq,
            &query,
            &mut iq,
            &query,
            &mut info,
        );
    }
    assert_eq!(info, 0, "dstemr query");
    let mut work = vec![0.0; wq as usize];
    let mut iwork = vec![0 as c_int; iq as usize];
    let (lw, liw) = (work.len() as c_int, iwork.len() as c_int);
    tryrac = 1;
    // SAFETY: d (m) and e (m, last entry workspace) are private copies, w has m,
    // z is m x c (ldz = m, nzc = c), isuppz 2c, work/iwork the queried lengths.
    unsafe {
        dstemr_(
            &ch(b'V'),
            &ch(b'I'),
            &n,
            d.as_mut_ptr(),
            e.as_mut_ptr(),
            &0.0,
            &0.0,
            &il,
            &iu,
            &mut found,
            w.as_mut_ptr(),
            z.as_mut_ptr(),
            &n,
            &nzc,
            isuppz.as_mut_ptr(),
            &mut tryrac,
            work.as_mut_ptr(),
            &lw,
            iwork.as_mut_ptr(),
            &liw,
            &mut info,
        );
    }
    assert_eq!(info, 0, "dstemr");
    assert_eq!(found as usize, c);
    w.truncate(c);
    (w, z)
}

/// Z (m x c) := Q Z with the Householder vectors of `t`.
fn ormtr(t: &Trd, z: &mut [f64], c: usize, work: &mut Vec<f64>) {
    let m = t.d.len();
    let (n, nc) = (m as c_int, c as c_int);
    let mut info = 0;
    if work.len() < m * 64 {
        work.resize(m * 64, 0.0);
    }
    let lwork = work.len() as c_int;
    // SAFETY: A/tau from dsytrd of order m (lda = m), C is m x c (ldc = m),
    // work has at least the minimum c (we pass m * 64 >= n * nb for nb <= 64).
    unsafe {
        dormtr_(
            &ch(b'L'),
            &ch(b'L'),
            &ch(b'N'),
            &n,
            &nc,
            t.a.as_ptr(),
            &n,
            t.tau.as_ptr(),
            z.as_mut_ptr(),
            &n,
            work.as_mut_ptr(),
            &lwork,
            &mut info,
        );
    }
    assert_eq!(info, 0, "dormtr");
}

/// Slice bounds of k roots over p ranks, each bound moved to the widest gap
/// within +-2 of the even split (never cutting below 1 root per slice).
fn bounds(w: &[f64], p: usize) -> Vec<usize> {
    let k = w.len();
    let mut b = vec![0];
    for r in 1..p {
        let even = r * k / p;
        let lo = (even.saturating_sub(2)).max(b[r - 1] + 1);
        let hi = (even + 2).min(k - 1);
        let best = (lo..=hi.max(lo))
            .filter(|&i| i > 0 && i < k)
            .max_by(|&x, &y| (w[x] - w[x - 1]).total_cmp(&(w[y] - w[y - 1])))
            .unwrap_or(even);
        b.push(best);
    }
    b.push(k);
    b
}

/// (max |Z^T Z - I|, max_j |A z_j - w_j z_j| / max |A|).
fn quality(a: &[f64], m: usize, w: &[f64], z: &[f64]) -> (f64, f64) {
    let k = w.len();
    let amax = a.iter().fold(0.0_f64, |s, x| s.max(x.abs()));
    let mut ortho = 0.0_f64;
    for i in 0..k {
        for j in 0..=i {
            let dot: f64 = (0..m).map(|r| z[r + i * m] * z[r + j * m]).sum();
            ortho = ortho.max((dot - if i == j { 1.0 } else { 0.0 }).abs());
        }
    }
    let mut res = 0.0_f64;
    for j in 0..k {
        for r in 0..m {
            let az: f64 = (0..m).map(|c| a[r + c * m] * z[c + j * m]).sum();
            res = res.max((az - w[j] * z[r + j * m]).abs() / amax);
        }
    }
    (ortho, res)
}

fn best_of(reps: usize, mut f: impl FnMut()) -> f64 {
    f();
    let mut best = f64::INFINITY;
    for _ in 0..reps {
        let t = Instant::now();
        f();
        best = best.min(t.elapsed().as_secs_f64());
    }
    best * 1e6
}

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let frac: f64 = args.first().map_or(0.34, |s| s.parse().expect("frac"));
    let reps: usize = args.get(1).map_or(20, |s| s.parse().expect("reps"));
    let orders = [136usize, 204, 272, 408, 544, 816];
    println!(
        "REDUCED-EIGEN-SLICES frac={frac} reps={reps} (us, best of reps; sliced excludes the allgather)"
    );
    println!(
        "model m k | dsyevd | trd vals | stein P=1/4/8 | mrrr P=1/4/8 | ormtr k/1,k/4,k/8 | sliced-stein P=4/8 (ratio) | sliced-mrrr P=4/8 (ratio) | ortho/res dsyevd, stein8, mrrr8"
    );
    for model in ["dense", "ritz"] {
        for &m in &orders {
            let k = ((m as f64 * frac).ceil() as usize).clamp(1, m);
            let a = matrix(m, model, 7 + m as u64);
            let mut evd = dft_mkl::SymmetricEigen::new(m).expect("dsyevd workspace");
            let mut vecs = vec![0.0; m * m];
            let mut vals = vec![0.0; m];
            let t_evd = best_of(reps, || {
                vecs.copy_from_slice(&a);
                evd.solve(m, &mut vecs, &mut vals).expect("dsyevd");
            });
            let q_evd = quality(&a, m, &vals[..k], &vecs[..m * k]);
            let mut trd_work = Vec::new();
            let t_trd = best_of(reps, || {
                std::hint::black_box(sytrd(m, &a, &mut trd_work));
            });
            let t = sytrd(m, &a, &mut trd_work);
            let t_vals = best_of(reps, || {
                std::hint::black_box(stebz(&t, k));
            });
            let v = stebz(&t, k);
            let mut orm_work = Vec::new();
            let mut stein_t = [0.0; 3];
            let mut mrrr_t = [0.0; 3];
            let mut orm_t = [0.0; 3];
            for (slot, &p) in [1usize, 4, 8].iter().enumerate() {
                let b = bounds(&v.w, p.min(k));
                // The largest slice decides the critical path.
                let (lo, hi) = (0..b.len() - 1)
                    .map(|r| (b[r], b[r + 1]))
                    .max_by_key(|&(lo, hi)| hi - lo)
                    .unwrap();
                stein_t[slot] = best_of(reps, || {
                    std::hint::black_box(stein(&t, &v, lo, hi));
                });
                mrrr_t[slot] = best_of(reps, || {
                    std::hint::black_box(mrrr(&t, lo, hi));
                });
                let mut z = stein(&t, &v, lo, hi);
                let base = z.clone();
                orm_t[slot] = best_of(reps, || {
                    z.copy_from_slice(&base);
                    ormtr(&t, &mut z, hi - lo, &mut orm_work);
                });
            }
            // Assembled slices at P = 8, each slice computed independently.
            let b = bounds(&v.w, 8.min(k));
            let mut zs = Vec::with_capacity(m * k);
            let mut zm = Vec::with_capacity(m * k);
            let mut wm = Vec::with_capacity(k);
            for r in 0..b.len() - 1 {
                let (lo, hi) = (b[r], b[r + 1]);
                let mut z = stein(&t, &v, lo, hi);
                ormtr(&t, &mut z, hi - lo, &mut orm_work);
                zs.extend_from_slice(&z);
                let (w, mut z) = mrrr(&t, lo, hi);
                ormtr(&t, &mut z, hi - lo, &mut orm_work);
                zm.extend_from_slice(&z);
                wm.extend_from_slice(&w);
            }
            let q_stein = quality(&a, m, &v.w, &zs);
            let q_mrrr = quality(&a, m, &wm, &zm);
            let sliced = |vec_t: &[f64; 3], i: usize| t_trd + t_vals + vec_t[i] + orm_t[i];
            let sliced_m = |i: usize| t_trd + mrrr_t[i] + orm_t[i];
            println!(
                "{model} {m} {k} | {t_evd:.0} | {t_trd:.0} {t_vals:.0} | {:.0}/{:.0}/{:.0} | {:.0}/{:.0}/{:.0} | {:.0}/{:.0}/{:.0} | {:.0} ({:.2}) / {:.0} ({:.2}) | {:.0} ({:.2}) / {:.0} ({:.2}) | {:.1e}/{:.1e} {:.1e}/{:.1e} {:.1e}/{:.1e}",
                stein_t[0],
                stein_t[1],
                stein_t[2],
                mrrr_t[0],
                mrrr_t[1],
                mrrr_t[2],
                orm_t[0],
                orm_t[1],
                orm_t[2],
                sliced(&stein_t, 1),
                sliced(&stein_t, 1) / t_evd,
                sliced(&stein_t, 2),
                sliced(&stein_t, 2) / t_evd,
                sliced_m(1),
                sliced_m(1) / t_evd,
                sliced_m(2),
                sliced_m(2) / t_evd,
                q_evd.0,
                q_evd.1,
                q_stein.0,
                q_stein.1,
                q_mrrr.0,
                q_mrrr.1,
            );
        }
    }
}
