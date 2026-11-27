//! Pair densities `rho_nm(q+G) = sum_G1 conj(c_n(G1)) c_m(G1 + G)` by FFT:
//! orbitals to real space, the product `conj(u_n(r)) u_m(r)`, forward to the
//! response sphere. Both directions use the batched stick FFT
//! (`DistributedPairFft`) of the plane-wave layout in its one-rank form: the
//! response driver distributes transitions over ranks (MPI level 1, plan
//! 3.4), so every pair density of a transition is formed on one rank and the
//! product is exact. Real-space fields are kept as full row-major grids
//! `[x][y][z]`, which decouples the first-axis choices of the orbital and
//! response plans.
//!
//! Exactness: with the grid of the density sphere (every Miller difference of
//! two orbital spheres representable without aliasing) and a response sphere
//! inside it (E_cut^chi <= 4 E_cut), the FFT result equals the convolution sum
//! up to roundoff.
use crate::{ResponseError, input, numerical};
use dft_mpi::Comm;
use dft_numerics::linalg::ComplexMatrix;
use dft_solvers::distributed_fft::DistributedPairFft;
use num_complex::Complex64 as C;

/// Columns per batched transform.
const BATCH: usize = 16;

fn grid_indices(dims: [usize; 3], miller: &[[i32; 3]]) -> Result<Vec<[usize; 3]>, ResponseError> {
    if dims.contains(&0) {
        return Err(input("empty FFT grid"));
    }
    let mut seen = vec![false; dims.iter().product()];
    let mut out = Vec::with_capacity(miller.len());
    for m in miller {
        let index: [usize; 3] =
            std::array::from_fn(|i| i64::from(m[i]).rem_euclid(dims[i] as i64) as usize);
        let flat = (index[0] * dims[1] + index[1]) * dims[2] + index[2];
        if seen[flat] {
            return Err(input(format!(
                "FFT grid {dims:?} cannot hold the plane-wave set (Miller {m:?} aliases)"
            )));
        }
        seen[flat] = true;
        out.push(index);
    }
    Ok(out)
}

/// One-rank batched stick plan of a coefficient set on the grid.
pub struct Plan {
    fft: DistributedPairFft,
    len: usize,
}
impl Plan {
    /// `miller` lists the coefficient rows (global order). A row with Miller
    /// index M is placed at grid point M - `shift`, i.e. the plan represents
    /// `u(r) = sum_M c(M) exp(i (M - shift).G r)`; `shift` = G0 of an umklapp
    /// `k + q = k' + G0` turns the orbitals of k' into those of k + q.
    pub fn new(
        dims: [usize; 3],
        miller: &[[i32; 3]],
        shift: [i32; 3],
    ) -> Result<Self, ResponseError> {
        let shifted: Vec<[i32; 3]> = miller
            .iter()
            .map(|m| std::array::from_fn(|i| m[i] - shift[i]))
            .collect();
        let indices = grid_indices(dims, &shifted)?;
        let mut fft = DistributedPairFft::new_direct(dims, &indices, 0, 1).map_err(numerical)?;
        // Grow the per-vector stick buffers of both directions to the batch
        // width once (the forward stage writes into them without growing).
        fft.begin_to_real_batch(
            &Comm::Single,
            &vec![C::default(); miller.len() * BATCH],
            BATCH,
        )
        .map_err(numerical)?;
        if fft.local_positions().len() != miller.len()
            || fft
                .local_positions()
                .iter()
                .enumerate()
                .any(|(i, &p)| i != p)
        {
            return Err(numerical(
                "one-rank stick plan does not own every row in order",
            ));
        }
        Ok(Self {
            fft,
            len: miller.len(),
        })
    }
    pub fn len(&self) -> usize {
        self.len
    }
    pub fn is_empty(&self) -> bool {
        self.len == 0
    }
    fn grid_len(&self) -> usize {
        self.fft.shape().iter().product()
    }

    /// Real-space values `u(r_j) = sum c e^{i G.r_j}` (no 1/sqrt(Omega)) of
    /// `columns` of `coefficients` (rows = this plan's rows), written to
    /// `out` as consecutive full grids (resized; reused across calls).
    pub fn to_real(
        &mut self,
        coefficients: &ComplexMatrix,
        columns: std::ops::Range<usize>,
        out: &mut Vec<C>,
    ) -> Result<(), ResponseError> {
        if coefficients.rows() != self.len || columns.end > coefficients.cols() {
            return Err(input("orbital block shape"));
        }
        let n = self.grid_len();
        let nz = self.fft.shape()[2];
        out.resize(columns.len() * n, C::default());
        let mut start = columns.start;
        while start < columns.end {
            let end = (start + BATCH).min(columns.end);
            let np = end - start;
            let data = &coefficients.as_slice()[start * self.len..end * self.len];
            self.fft
                .begin_to_real_batch(&Comm::Single, data, np)
                .map_err(numerical)?;
            for t in 0..np {
                let field = &mut out[(start - columns.start + t) * n..][..n];
                for i in 0..self.fft.plane_count() {
                    self.fft.to_real_plane_in(t, np, i).map_err(numerical)?;
                    for y in 0..self.fft.rows() {
                        let at = self.fft.grid_start(i, y);
                        field[at..at + nz].copy_from_slice(&self.fft.plane()[y * nz..(y + 1) * nz]);
                    }
                }
            }
            start = end;
        }
        Ok(())
    }

    /// Forward `(1/N) sum_j f(r_j) e^{-i G.r_j}` of `count` consecutive full
    /// grids into columns `col..col + count` of `out` (rows = this plan's rows).
    pub fn to_reciprocal(
        &mut self,
        fields: &[C],
        count: usize,
        out: &mut ComplexMatrix,
        col: usize,
    ) -> Result<(), ResponseError> {
        let n = self.grid_len();
        let nz = self.fft.shape()[2];
        if fields.len() < count * n || out.rows() != self.len || col + count > out.cols() {
            return Err(input("pair field batch shape"));
        }
        let mut start = 0;
        while start < count {
            let end = (start + BATCH).min(count);
            let np = end - start;
            for t in 0..np {
                let field = &fields[(start + t) * n..][..n];
                for i in 0..self.fft.plane_count() {
                    for y in 0..self.fft.rows() {
                        let at = self.fft.grid_start(i, y);
                        self.fft.plane_mut()[y * nz..(y + 1) * nz]
                            .copy_from_slice(&field[at..at + nz]);
                    }
                    self.fft
                        .to_reciprocal_plane_in(t, np, i)
                        .map_err(numerical)?;
                }
            }
            let slice = &mut out.as_mut_slice()[(col + start) * self.len..(col + end) * self.len];
            self.fft
                .finish_to_reciprocal_batch(&Comm::Single, slice, np)
                .map_err(numerical)?;
            start = end;
        }
        Ok(())
    }
}

/// Pair densities of one left orbital against a block of right orbitals.
pub struct PairFft {
    response: Plan,
    product: Vec<C>,
}
impl PairFft {
    /// `response_miller`: the response sphere {G : |q+G|^2/2 <= E_cut^chi}
    /// (Miller indices relative to the q used for the orbitals' k + q).
    pub fn new(dims: [usize; 3], response_miller: &[[i32; 3]]) -> Result<Self, ResponseError> {
        Ok(Self {
            response: Plan::new(dims, response_miller, [0; 3])?,
            product: Vec::new(),
        })
    }
    pub fn response_len(&self) -> usize {
        self.response.len()
    }
    /// `rho_nm(q+G)` for the left field `u_n` (one full grid) and `count`
    /// right fields `u_m,k+q` (consecutive full grids), into columns
    /// `col..col + count` of `out` (rows: the response sphere).
    pub fn pairs(
        &mut self,
        left: &[C],
        right: &[C],
        count: usize,
        out: &mut ComplexMatrix,
        col: usize,
    ) -> Result<(), ResponseError> {
        let n = left.len();
        if n != self.response.grid_len() || right.len() < count * n {
            return Err(input("pair field shape"));
        }
        self.product.resize(count * n, C::default());
        for (t, product) in self.product.chunks_exact_mut(n).take(count).enumerate() {
            for ((p, l), r) in product.iter_mut().zip(left).zip(&right[t * n..(t + 1) * n]) {
                *p = l.conj() * r;
            }
        }
        let product = std::mem::take(&mut self.product);
        let result = self.response.to_reciprocal(&product, count, out, col);
        self.product = product;
        result
    }
}
