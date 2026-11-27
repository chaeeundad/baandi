//! chi0(q, i w) accumulation by GEMM (plan 2.1, decision D1). For every
//! frequency the transition block contributes `-(Q Q^H)` with
//! `Q = P diag(sqrt(-f_j(w)))`, `f_j(w) = s 2 D_j / (D_j^2 + w^2) < 0`, so each
//! chi0(i w) is Hermitian negative semidefinite by construction.
//!
//! Determinism: frequencies are independent; each is accumulated in the order
//! of the `add` calls with the sequential GEMM backend, whatever the thread
//! count. Partial accumulators of MPI ranks are combined by `tree_sum_f64`
//! (fixed association, ADR 0041): bit-reproducible at a fixed rank count, and
//! one rank is bit-identical to the threads-only run.
use crate::{ResponseError, input, numerical};
use dft_mpi::Collective;
use dft_numerics::linalg::{ComplexMatrix, MatrixOp, gemm_into};
use num_complex::Complex64 as C;

/// Accumulated response at one q for every frequency.
#[derive(Clone, Debug)]
pub struct Chi0 {
    /// Frequencies (Ha) on the imaginary axis, as given.
    pub frequencies: Vec<f64>,
    /// chi0_GG' per frequency (Bohr^-3 Ha^-1, response sphere order). In
    /// optical (q -> 0) mode the G = 0 row and column are zero: their limits
    /// are `head` and `wings`.
    pub body: Vec<ComplexMatrix>,
    /// Optical mode: the 3x3 tensor L with `chi0_00(q) = q^2 q^.L.q^`.
    pub head: Vec<[[C; 3]; 3]>,
    /// Optical mode: N_G x 3, `chi0_G0(q) = q sum_b q^_b W[(G, b)]` (row G = 0 zero).
    pub wings: Vec<ComplexMatrix>,
    /// Transitions added over all ranks.
    pub transitions: u64,
}

struct Workspace {
    scaled: ComplexMatrix,
    product: ComplexMatrix,
}

/// Rank-local accumulator.
pub struct Accumulator {
    n_g: usize,
    optical: bool,
    threads: usize,
    chi0: Chi0,
    workspaces: Vec<Workspace>,
}

impl Accumulator {
    /// `optical`: q -> 0 with head and wings (the G = 0 row of the pair
    /// densities is ignored and the q-linear coefficients are supplied).
    pub fn new(
        n_g: usize,
        frequencies: &[f64],
        optical: bool,
        threads: usize,
    ) -> Result<Self, ResponseError> {
        if n_g == 0
            || frequencies.is_empty()
            || frequencies.iter().any(|w| !w.is_finite() || *w < 0.)
        {
            return Err(input("chi0: empty basis or invalid frequencies"));
        }
        let nw = frequencies.len();
        let threads = threads.clamp(1, nw);
        let zeros = |r, c| ComplexMatrix::zeros(r, c).map_err(numerical);
        let mut body = Vec::with_capacity(nw);
        let mut wings = Vec::with_capacity(nw);
        for _ in 0..nw {
            body.push(zeros(n_g, n_g)?);
            wings.push(zeros(n_g, if optical { 3 } else { 0 })?);
        }
        let mut workspaces = Vec::with_capacity(threads);
        for _ in 0..threads {
            workspaces.push(Workspace {
                scaled: zeros(n_g, 0)?,
                product: zeros(n_g, n_g)?,
            });
        }
        Ok(Self {
            n_g,
            optical,
            threads,
            chi0: Chi0 {
                frequencies: frequencies.to_vec(),
                body,
                head: vec![[[C::default(); 3]; 3]; if optical { nw } else { 0 }],
                wings,
                transitions: 0,
            },
            workspaces,
        })
    }

    /// Add `pairs.cols()` transitions: `pairs[(G, j)] = rho_j(q+G)`,
    /// `delta[j] = e_occ - e_emp < 0` (Ha), `scale = 2 w_k / Omega`
    /// (spin-unpolarized; Bohr^-3), and in optical mode `slopes[j] = a_j`
    /// with `rho_j(q) = q . a_j` (Bohr).
    pub fn add(
        &mut self,
        pairs: &ComplexMatrix,
        delta: &[f64],
        scale: f64,
        slopes: Option<&[[C; 3]]>,
    ) -> Result<(), ResponseError> {
        let t = pairs.cols();
        if pairs.rows() != self.n_g
            || delta.len() != t
            || self.optical != slopes.is_some()
            || slopes.is_some_and(|s| s.len() != t)
            || !scale.is_finite()
            || scale <= 0.
        {
            return Err(input("chi0: transition block shape"));
        }
        if delta.iter().any(|d| !d.is_finite() || *d >= 0.) {
            return Err(input(
                "chi0: every transition needs e_occ < e_emp (insulator, D6)",
            ));
        }
        if t == 0 {
            return Ok(());
        }
        pairs.validate_finite().map_err(numerical)?;
        let optical = self.optical;
        let per = self.chi0.frequencies.len().div_ceil(self.threads);
        let Chi0 {
            frequencies,
            body,
            head,
            wings,
            ..
        } = &mut self.chi0;
        let mut heads: Vec<Option<&mut [[C; 3]; 3]>> = if optical {
            head.iter_mut().map(Some).collect()
        } else {
            (0..frequencies.len()).map(|_| None).collect()
        };
        let results: Vec<Result<(), ResponseError>> = std::thread::scope(|scope| {
            let mut handles = Vec::new();
            for ((((w, b), wi), h), ws) in frequencies
                .chunks(per)
                .zip(body.chunks_mut(per))
                .zip(wings.chunks_mut(per))
                .zip(heads.chunks_mut(per))
                .zip(self.workspaces.iter_mut())
            {
                handles.push(scope.spawn(move || {
                    for (((&omega, body), wing), head) in w
                        .iter()
                        .zip(b.iter_mut())
                        .zip(wi.iter_mut())
                        .zip(h.iter_mut())
                    {
                        accumulate(pairs, delta, scale, slopes, omega, body, wing, head, ws)?;
                    }
                    Ok(())
                }));
            }
            handles
                .into_iter()
                .map(|h| {
                    h.join()
                        .unwrap_or_else(|_| Err(numerical("chi0 worker panicked")))
                })
                .collect()
        });
        for r in results {
            r?;
        }
        self.chi0.transitions += t as u64;
        Ok(())
    }

    /// Combine the rank partial sums (fixed-association tree); every rank
    /// receives identical bits. Collective.
    pub fn reduce<K: Collective>(self, comm: &K) -> Result<Chi0, ResponseError> {
        let mut chi0 = self.chi0;
        let mut flat = Vec::new();
        let push = |flat: &mut Vec<f64>, z: &C| flat.extend([z.re, z.im]);
        for m in chi0.body.iter().chain(&chi0.wings) {
            m.as_slice().iter().for_each(|z| push(&mut flat, z));
        }
        for h in &chi0.head {
            h.iter().flatten().for_each(|z| push(&mut flat, z));
        }
        flat.push(chi0.transitions as f64);
        let mut sum = vec![0.; flat.len()];
        comm.tree_sum_f64(&flat, &mut sum).map_err(numerical)?;
        let mut it = sum.as_chunks::<2>().0.iter();
        let mut next = || {
            let p = it.next().expect("reduced length");
            C::new(p[0], p[1])
        };
        for m in chi0.body.iter_mut().chain(chi0.wings.iter_mut()) {
            m.as_mut_slice().iter_mut().for_each(|z| *z = next());
        }
        for h in &mut chi0.head {
            h.iter_mut().flatten().for_each(|z| *z = next());
        }
        chi0.transitions = *sum.last().expect("count") as u64;
        Ok(chi0)
    }
}

#[allow(clippy::too_many_arguments)]
fn accumulate(
    pairs: &ComplexMatrix,
    delta: &[f64],
    scale: f64,
    slopes: Option<&[[C; 3]]>,
    omega: f64,
    body: &mut ComplexMatrix,
    wing: &mut ComplexMatrix,
    head: &mut Option<&mut [[C; 3]; 3]>,
    ws: &mut Workspace,
) -> Result<(), ResponseError> {
    let t = pairs.cols();
    crate::set_columns(&mut ws.scaled, t)?;
    let factor = |d: f64| scale * 2. * d / (d * d + omega * omega);
    for (j, &d) in delta.iter().enumerate() {
        let s = (-factor(d)).sqrt();
        let column = ws.scaled.column_mut(j);
        for (q, p) in column.iter_mut().zip(pairs.column(j)) {
            *q = p * s;
        }
        if slopes.is_some() {
            column[0] = C::default();
        }
    }
    gemm_into(
        &ws.scaled,
        MatrixOp::Identity,
        &ws.scaled,
        MatrixOp::Adjoint,
        &mut ws.product,
    )
    .map_err(numerical)?;
    for (c, p) in body.as_mut_slice().iter_mut().zip(ws.product.as_slice()) {
        *c -= p;
    }
    if let (Some(slopes), Some(head)) = (slopes, head.as_deref_mut()) {
        for (j, (&d, a)) in delta.iter().zip(slopes).enumerate() {
            let f = factor(d);
            for (x, row) in head.iter_mut().enumerate() {
                for (y, h) in row.iter_mut().enumerate() {
                    *h += f * a[x] * a[y].conj();
                }
            }
            let rho = pairs.column(j);
            for (y, ay) in a.iter().enumerate() {
                let c = f * ay.conj();
                for (g, w) in wing.column_mut(y).iter_mut().enumerate().skip(1) {
                    *w += rho[g] * c;
                }
            }
        }
    }
    Ok(())
}
