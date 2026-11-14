//! Row-block distributed dense algebra for the real Davidson (ADR 0041, M1b).
//!
//! The plane-wave rows of a matrix are split into fixed contiguous blocks, one
//! per rank (`dft_mpi::block_counts`). A product that sums over the rows (V^T W,
//! Q^T X) is a local GEMM on the local rows followed by the deterministic tree
//! sum, so every rank holds the same bits and the result for a fixed rank count
//! does not depend on the library's reduction tree. Products that keep the row
//! index (rotation, projection back-application) are purely local. With one rank
//! every primitive performs exactly the serial GEMM calls (bit-identical).
use crate::SolverError;
use dft_mpi::{Collective, block_counts, block_displacements};
use dft_numerics::{
    GemmSite,
    real::{RealGemmWorkspace, RealMatrix, RealOp},
};

fn comm_error(e: dft_mpi::MpiError) -> SolverError {
    SolverError::Application(e.to_string())
}

/// Row-block partition of `total` rows over `ranks` for one rank.
#[derive(Debug, Clone)]
pub struct RowBlock {
    rank: usize,
    counts: Vec<usize>,
    displs: Vec<usize>,
}
impl RowBlock {
    pub fn new(total: usize, rank: usize, ranks: usize) -> Self {
        let counts = block_counts(total, ranks);
        let displs = block_displacements(&counts);
        Self {
            rank,
            counts,
            displs,
        }
    }
    pub fn total(&self) -> usize {
        self.counts.iter().sum()
    }
    pub fn local_rows(&self) -> usize {
        self.counts[self.rank]
    }
    pub fn offset(&self) -> usize {
        self.displs[self.rank]
    }
    /// The local rows of a full (replicated) matrix; used by tests and by the
    /// setup that distributes a replicated start.
    pub fn scatter(&self, full: &RealMatrix) -> Result<RealMatrix, SolverError> {
        let (rows, cols) = (self.local_rows(), full.cols());
        let mut data = Vec::with_capacity(rows * cols);
        for j in 0..cols {
            let column = &full.as_slice()[j * full.rows()..(j + 1) * full.rows()];
            data.extend_from_slice(&column[self.offset()..self.offset() + rows]);
        }
        Ok(RealMatrix::new(rows, cols, data)?)
    }
    /// Reassemble the full matrix on every rank (column by column allgather).
    pub fn gather<C: Collective>(
        &self,
        comm: &C,
        local: &RealMatrix,
    ) -> Result<RealMatrix, SolverError> {
        if local.rows() != self.local_rows() {
            return Err(SolverError::InvalidParameter("local row count"));
        }
        let (total, cols) = (self.total(), local.cols());
        let mut data = vec![0.0; total * cols];
        for j in 0..cols {
            let column = &local.as_slice()[j * local.rows()..(j + 1) * local.rows()];
            comm.allgatherv_f64(column, &self.counts, &mut data[j * total..(j + 1) * total])
                .map_err(comm_error)?;
        }
        Ok(RealMatrix::new(total, cols, data)?)
    }
}

/// out (ma x nb) = sum over all ranks of A_local^T B_local; identical on every
/// rank. `a` and `b` are the local row blocks with the same row count.
pub fn gram<C: Collective>(
    comm: &C,
    ws: &mut RealGemmWorkspace,
    a: &RealMatrix,
    b: &RealMatrix,
    out: &mut RealMatrix,
) -> Result<(), SolverError> {
    let rows = a.rows();
    if b.rows() != rows || out.rows() != a.cols() || out.cols() != b.cols() {
        return Err(SolverError::InvalidParameter("gram dimensions"));
    }
    let mut partial = RealMatrix::zeros(out.rows(), out.cols())?;
    if rows > 0 && a.cols() > 0 && b.cols() > 0 {
        ws.multiply_slices_into(
            GemmSite::Other,
            a.as_slice(),
            rows,
            a.cols(),
            rows,
            RealOp::Transpose,
            b.as_slice(),
            b.cols(),
            rows,
            partial.as_mut(),
        )?;
    }
    comm.tree_sum_f64(partial.as_slice(), out.as_mut_slice())
        .map_err(comm_error)
}

/// out_local (rows x cols) = A_local (rows x k) times C (k x cols), C replicated.
pub fn rotate_local(
    ws: &mut RealGemmWorkspace,
    a: &RealMatrix,
    c: &RealMatrix,
    out: &mut RealMatrix,
) -> Result<(), SolverError> {
    rotate_local_with(ws, a, c, out, false)
}

/// `rotate_local` with the GEMM route chosen: `wide` (`BAANDI_M2_DGEMM=on`) takes
/// `multiply_slices_wide_into`, so on MKL builds the product reaches the
/// per-panel dgemm for any output size (only the fixed 256-row panel bound
/// applies) instead of falling back to the pure kernel once rows x cols exceeds
/// the square workspace bound; on pure builds both routes are the same product
/// (bit-identical). Panels are fixed by the shapes, so the result does not
/// depend on the thread count.
pub fn rotate_local_with(
    ws: &mut RealGemmWorkspace,
    a: &RealMatrix,
    c: &RealMatrix,
    out: &mut RealMatrix,
    wide: bool,
) -> Result<(), SolverError> {
    if c.rows() != a.cols() || out.rows() != a.rows() || out.cols() != c.cols() {
        return Err(SolverError::InvalidParameter("rotate dimensions"));
    }
    if a.rows() == 0 || a.cols() == 0 || c.cols() == 0 {
        out.as_mut_slice().fill(0.0);
        return Ok(());
    }
    let product = if wide {
        RealGemmWorkspace::multiply_slices_wide_into
    } else {
        RealGemmWorkspace::multiply_slices_into
    };
    product(
        ws,
        GemmSite::Other,
        a.as_slice(),
        a.rows(),
        a.cols(),
        a.rows(),
        RealOp::Identity,
        c.as_slice(),
        c.cols(),
        c.rows(),
        out.as_mut(),
    )?;
    Ok(())
}

/// x_local -= Q_local (Q^T X): the block projection against the columns of Q.
/// `overlap` (k x w) is left holding Q^T X (identical on every rank), `scratch`
/// (rows x w) is workspace.
pub fn project_out<C: Collective>(
    comm: &C,
    ws: &mut RealGemmWorkspace,
    q: &RealMatrix,
    x: &mut RealMatrix,
    overlap: &mut RealMatrix,
    scratch: &mut RealMatrix,
) -> Result<(), SolverError> {
    project_out_with(comm, ws, q, x, overlap, scratch, false)
}

/// `project_out` with the rotation route of `rotate_local_with`.
pub fn project_out_with<C: Collective>(
    comm: &C,
    ws: &mut RealGemmWorkspace,
    q: &RealMatrix,
    x: &mut RealMatrix,
    overlap: &mut RealMatrix,
    scratch: &mut RealMatrix,
    wide: bool,
) -> Result<(), SolverError> {
    gram(comm, ws, q, x, overlap)?;
    rotate_local_with(ws, q, overlap, scratch, wide)?;
    for (dst, &src) in x.as_mut_slice().iter_mut().zip(scratch.as_slice()) {
        *dst -= src;
    }
    Ok(())
}

/// Squared column norms of a row-distributed block (one reduction).
pub fn column_norms_squared<C: Collective>(
    comm: &C,
    x: &RealMatrix,
) -> Result<Vec<f64>, SolverError> {
    let rows = x.rows();
    let partial: Vec<f64> = (0..x.cols())
        .map(|j| {
            x.as_slice()[j * rows..(j + 1) * rows]
                .iter()
                .map(|v| v * v)
                .sum::<f64>()
        })
        .collect();
    let mut out = vec![0.0; x.cols()];
    comm.tree_sum_f64(&partial, &mut out).map_err(comm_error)?;
    Ok(out)
}
