//! Owned-G complex Bloch Hamiltonian for the reference k-point extension.
//! Each rank participates in every k point. Coefficients are the periodic
//! Bloch part; kinetic/projectors use k+G, while the FFT uses integer G.
//! No Gamma pair weights or time-reversal multiplicities belong here.
use crate::hamiltonian::NonlocalOperator;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_mpi::Collective;
use dft_numerics::{
    linalg::{ComplexMatrix, GemmWorkspace, MatrixOp},
    pool::WorkerPool,
};
use dft_solvers::{SolverError, distributed_fft::DistributedPairFft};
use num_complex::Complex64 as C;
use std::sync::Arc;

fn error(e: impl std::fmt::Display) -> SolverError {
    SolverError::Application(e.to_string())
}

/// Kinetic plus nonlocal diagonal per owned row:
/// `T_g + sum_atom sum_ij D_ij P_gi conj(P_gj)`, added per row in the order
/// (atom, i, j). GK16: the pair loop is outside the row loop, so each pass
/// reads two contiguous projector columns, and pairs with an exactly zero
/// coupling (different l or m) are skipped. Bit-identical to the row-outer
/// form: every row receives the same terms in the same order, and a skipped
/// term is a signed zero (finite projectors, validated by the caller) whose
/// addition cannot change a row value that starts at T_g >= +0 and never
/// becomes -0 under round-to-nearest.
fn static_diagonal(
    kinetic: &[f64],
    projectors: &ComplexMatrix,
    channels: &[usize],
    couplings: &[Vec<f64>],
) -> Vec<f64> {
    let mut diagonal = kinetic.to_vec();
    let mut offset = 0;
    for (&count, coupling) in channels.iter().zip(couplings) {
        for i in 0..count {
            for j in 0..count {
                let d = coupling[i * count + j];
                if d == 0. {
                    continue;
                }
                let (pi, pj) = (projectors.column(offset + i), projectors.column(offset + j));
                for ((value, a), b) in diagonal.iter_mut().zip(pi).zip(pj) {
                    *value += (d * *a * b.conj()).re;
                }
            }
        }
        offset += count;
    }
    diagonal
}

/// The layout, k point, kinetic energy and projectors are immutable. Construct
/// a distinct operator for each k point; a potential update changes only V(r).
/// Scratch is bounded by `capacity` bands, with no full orbital gather.
pub struct DistributedComplexHamiltonian {
    fft: DistributedPairFft,
    global_len: usize,
    capacity: usize,
    rank: usize,
    ranks: usize,
    kinetic: Vec<f64>,
    diagonal: Vec<f64>,
    projectors: ComplexMatrix,
    couplings: Vec<Vec<f64>>,
    channels: Vec<usize>,
    /// GK16: V(r) on the owned planes only (`plane_points()` values in the
    /// plane-local layout of `set_owned_potential`), not a full grid per k.
    potential: Vec<f64>,
    epoch: Option<u64>,
    identity: Arc<()>,
    pool: Arc<WorkerPool>,
    gemm: GemmWorkspace,
    projection: ComplexMatrix,
    coupled: ComplexMatrix,
    back: ComplexMatrix,
    partial: Vec<f64>,
    reduced: Vec<f64>,
}

impl DistributedComplexHamiltonian {
    #[allow(clippy::too_many_arguments)]
    pub fn new_owned(
        basis: &PlaneWaveBasis,
        grid: &FftGridSpec,
        nonlocal: &NonlocalOperator,
        capacity: usize,
        rank: usize,
        ranks: usize,
        pool: Arc<WorkerPool>,
    ) -> Result<Self, SolverError> {
        if ranks == 0
            || rank >= ranks
            || capacity == 0
            || grid.g_to_fft().len() != basis.vectors().len()
        {
            return Err(SolverError::InvalidParameter("complex distributed layout"));
        }
        for (g, &index) in basis.vectors().iter().zip(grid.g_to_fft()) {
            if FftGridSpec::map_miller(g.miller, grid.dims()).map_err(error)? != index {
                return Err(SolverError::InvalidParameter(
                    "complex distributed grid mapping",
                ));
            }
        }
        let fft =
            DistributedPairFft::new(grid.dims(), grid.g_to_fft(), rank, ranks).map_err(error)?;
        let rows = fft.local_len();
        let plane_points = fft.plane_points();
        let p = nonlocal.projector_count();
        let projectors = ComplexMatrix::new(
            rows,
            p,
            nonlocal
                .projector_rows(fft.local_positions())
                .map_err(error)?,
        )?;
        projectors.validate_finite()?;
        let channels: Vec<_> = nonlocal.atoms().iter().map(|a| a.channel_count()).collect();
        let couplings: Vec<_> = nonlocal
            .atoms()
            .iter()
            .map(|a| a.coupling_hartree().to_vec())
            .collect();
        let kinetic: Vec<_> = fft
            .local_positions()
            .iter()
            .map(|&g| basis.momenta()[g].norm_squared / 2.)
            .collect();
        let diagonal = static_diagonal(&kinetic, &projectors, &channels, &couplings);
        if diagonal.iter().any(|v| !v.is_finite()) {
            return Err(SolverError::InvalidParameter("complex nonfinite diagonal"));
        }
        Ok(Self {
            fft,
            global_len: basis.vectors().len(),
            capacity,
            rank,
            ranks,
            kinetic,
            diagonal,
            projectors,
            couplings,
            channels,
            potential: vec![f64::NAN; plane_points],
            epoch: None,
            identity: Arc::new(()),
            pool: pool.clone(),
            gemm: GemmWorkspace::with_pool(
                rows.max(p),
                capacity.max(p),
                pool.workers_created().max(1),
                pool,
            )?,
            projection: ComplexMatrix::zeros(p, capacity)?,
            coupled: ComplexMatrix::zeros(p, capacity)?,
            back: ComplexMatrix::zeros(rows, capacity)?,
            partial: vec![0.; 2 * p * capacity],
            reduced: vec![0.; 2 * p * capacity],
        })
    }

    /// Named owned numerical arrays; excludes opaque FFT/GEMM scratch.
    pub fn payload_bytes(&self) -> usize {
        8 * (self.kinetic.len()
            + self.diagonal.len()
            + self.potential.len()
            + self.partial.len()
            + self.reduced.len()
            + self.couplings.iter().map(Vec::len).sum::<usize>())
            + 16 * (self.projectors.as_slice().len()
                + self.projection.as_slice().len()
                + self.coupled.as_slice().len()
                + self.back.as_slice().len())
    }
    pub fn global_len(&self) -> usize {
        self.global_len
    }
    pub fn local_rows(&self) -> &[usize] {
        self.fft.local_positions()
    }
    pub fn kinetic_diagonal(&self) -> &[f64] {
        &self.kinetic
    }
    /// Kinetic plus nonlocal diagonal (Ha). Add the mean current local V.
    pub fn diagonal_static(&self) -> &[f64] {
        &self.diagonal
    }
    pub fn identity(&self) -> Arc<()> {
        self.identity.clone()
    }
    pub fn owned_epoch(&self) -> Option<u64> {
        self.epoch
    }
    pub fn worker_pool(&self) -> Arc<WorkerPool> {
        self.pool.clone()
    }
    /// Scratch band capacity (columns per apply).
    pub fn capacity(&self) -> usize {
        self.capacity
    }
    /// Grow the band scratch to at least `capacity` columns. A no-op when the
    /// scratch is already large enough; the potential, epoch and identity are
    /// unchanged, so cached H X entries keyed on this operator stay valid.
    pub fn ensure_capacity(&mut self, capacity: usize) -> Result<(), SolverError> {
        if capacity <= self.capacity {
            return Ok(());
        }
        let rows = self.local_rows().len();
        let p = self.projectors.cols();
        self.gemm = GemmWorkspace::with_pool(
            rows.max(p),
            capacity.max(p),
            self.pool.workers_created().max(1),
            self.pool.clone(),
        )?;
        self.projection = ComplexMatrix::zeros(p, capacity)?;
        self.coupled = ComplexMatrix::zeros(p, capacity)?;
        self.back = ComplexMatrix::zeros(rows, capacity)?;
        self.partial = vec![0.; 2 * p * capacity];
        self.reduced = vec![0.; 2 * p * capacity];
        self.capacity = capacity;
        Ok(())
    }
    pub fn first_axis(&self) -> usize {
        self.fft.first_axis()
    }
    pub fn plane_points(&self) -> usize {
        self.fft.plane_points()
    }
    pub fn plane_offsets(&self) -> Vec<usize> {
        let nz = self.fft.shape()[2];
        (0..self.fft.plane_count())
            .flat_map(|i| {
                (0..self.fft.rows()).flat_map(move |j| {
                    let start = self.fft.grid_start(i, j);
                    start..start + nz
                })
            })
            .collect()
    }

    /// Callers agree errors before any subsequent collective. Only the owned
    /// planes are stored (GK16; formerly a full-grid array per operator whose
    /// unowned part stayed poisoned).
    pub fn set_owned_potential(&mut self, values: &[f64], epoch: u64) -> Result<(), SolverError> {
        if values.len() != self.plane_points()
            || values.iter().any(|v| !v.is_finite())
            || self.epoch.is_some_and(|old| epoch <= old)
        {
            return Err(SolverError::InvalidParameter(
                "complex owned potential/epoch",
            ));
        }
        self.potential.copy_from_slice(values);
        self.epoch = Some(epoch);
        Ok(())
    }

    fn agree<K: Collective>(&self, comm: &K, bad: bool) -> Result<(), SolverError> {
        let mut result = [0.];
        comm.tree_sum_f64(
            &[f64::from(
                bad || comm.rank() != self.rank || comm.size() != self.ranks,
            )],
            &mut result,
        )
        .map_err(error)?;
        if result[0] != 0. {
            Err(SolverError::InvalidParameter(
                "complex distributed input on a rank",
            ))
        } else {
            Ok(())
        }
    }

    fn project<K: Collective>(&mut self, comm: &K, x: &ComplexMatrix) -> Result<(), SolverError> {
        use MatrixOp::{Adjoint, Identity};
        let p = self.projectors.cols();
        self.projection.resize_columns(x.cols())?;
        self.coupled.resize_columns(x.cols())?;
        if p == 0 {
            return Ok(());
        }
        let product =
            self.gemm
                .multiply_into(&self.projectors, Adjoint, x, Identity, &mut self.projection);
        self.agree(comm, product.is_err())?;
        product?;
        let len = 2 * p * x.cols();
        for (pair, z) in self.partial[..len]
            .as_chunks_mut::<2>()
            .0
            .iter_mut()
            .zip(self.projection.as_slice())
        {
            pair.copy_from_slice(&[z.re, z.im]);
        }
        comm.tree_sum_f64(&self.partial[..len], &mut self.reduced[..len])
            .map_err(error)?;
        for (z, pair) in self
            .projection
            .as_mut_slice()
            .iter_mut()
            .zip(self.reduced[..len].as_chunks::<2>().0)
        {
            *z = C::new(pair[0], pair[1]);
        }
        for band in 0..x.cols() {
            let mut offset = 0;
            for (&count, d) in self.channels.iter().zip(&self.couplings) {
                for i in 0..count {
                    self.coupled[(offset + i, band)] = (0..count)
                        .map(|j| d[i * count + j] * self.projection[(offset + j, band)])
                        .sum();
                }
                offset += count;
            }
        }
        Ok(())
    }

    pub fn apply_block<K: Collective>(
        &mut self,
        comm: &K,
        x: &ComplexMatrix,
        out: &mut ComplexMatrix,
    ) -> Result<(), SolverError> {
        self.agree(
            comm,
            self.epoch.is_none()
                || x.rows() != self.local_rows().len()
                || x.cols() > self.capacity
                || out.rows() != x.rows()
                || out.cols() != x.cols()
                || x.validate_finite().is_err(),
        )?;
        if x.cols() == 0 {
            return Ok(());
        }
        out.as_mut_slice().copy_from_slice(x.as_slice());
        self.fft
            .apply_local_batch_planes(comm, out.as_mut_slice(), &self.potential, x.cols())
            .map_err(error)?;
        for band in 0..x.cols() {
            for (g, &t) in self.kinetic.iter().enumerate() {
                out[(g, band)] += t * x[(g, band)];
            }
        }
        self.project(comm, x)?;
        if self.projectors.cols() > 0 {
            self.back.resize_columns(x.cols())?;
            let product = self.gemm.multiply_into(
                &self.projectors,
                MatrixOp::Identity,
                &self.coupled,
                MatrixOp::Identity,
                &mut self.back,
            );
            self.agree(comm, product.is_err())?;
            product?;
            for (y, z) in out.as_mut_slice().iter_mut().zip(self.back.as_slice()) {
                *y += z;
            }
        }
        self.agree(comm, out.validate_finite().is_err())
    }

    /// Add (V_new - V_old) X to caller-owned HX. The caller owns the exact
    /// X/HX/potential identity check; this kernel never decides cache validity.
    /// Plane values are reconstructed, with no retained real-space slab cache.
    pub fn add_local_change<K: Collective>(
        &mut self,
        comm: &K,
        x: &ComplexMatrix,
        change: &[f64],
        hx: &mut ComplexMatrix,
    ) -> Result<(), SolverError> {
        self.agree(
            comm,
            x.rows() != self.local_rows().len()
                || x.cols() > self.capacity
                || hx.rows() != x.rows()
                || hx.cols() != x.cols()
                || change.len() != self.plane_points()
                || change.iter().any(|v| !v.is_finite())
                || x.validate_finite().is_err()
                || hx.validate_finite().is_err(),
        )?;
        if x.cols() == 0 {
            return Ok(());
        }
        self.back.resize_columns(x.cols())?;
        self.back.as_mut_slice().copy_from_slice(x.as_slice());
        // GK16: the plane-local change is used in place (no full-grid copy).
        self.fft
            .apply_local_batch_planes(comm, self.back.as_mut_slice(), change, x.cols())
            .map_err(error)?;
        for (h, d) in hx.as_mut_slice().iter_mut().zip(self.back.as_slice()) {
            *h += d;
        }
        self.agree(comm, hx.validate_finite().is_err())
    }

    /// Sum f_b |u_b(r)|² / volume in band order on owned planes. The driver
    /// supplies weights f_b = occupation_b * k_weight, including TR reduction.
    pub fn band_density_sum<K: Collective>(
        &mut self,
        comm: &K,
        x: &ComplexMatrix,
        weights: &[f64],
        volume: f64,
        density: &mut [f64],
    ) -> Result<(), SolverError> {
        self.agree(
            comm,
            x.rows() != self.local_rows().len()
                || x.cols() > self.capacity
                || weights.len() != x.cols()
                || density.len() != self.plane_points()
                || !volume.is_finite()
                || volume <= 0.
                || x.validate_finite().is_err()
                || weights.iter().any(|v| !v.is_finite() || *v < 0.),
        )?;
        density.fill(0.);
        if x.cols() == 0 {
            return Ok(());
        }
        self.fft
            .begin_to_real_batch(comm, x.as_slice(), x.cols())
            .map_err(error)?;
        let plane_len = self.fft.rows() * self.fft.shape()[2];
        for (band, &weight) in weights.iter().enumerate() {
            for i in 0..self.fft.plane_count() {
                self.fft
                    .to_real_plane_in(band, x.cols(), i)
                    .map_err(error)?;
                for (n, z) in density[i * plane_len..(i + 1) * plane_len]
                    .iter_mut()
                    .zip(self.fft.plane())
                {
                    *n += weight * z.norm_sqr() / volume;
                }
            }
        }
        self.agree(comm, density.iter().any(|v| !v.is_finite()))
    }

    /// Per-band nonlocal expectation in Ha, already reduced across G owners.
    pub fn nonlocal_band_energies<K: Collective>(
        &mut self,
        comm: &K,
        x: &ComplexMatrix,
        energies: &mut [f64],
    ) -> Result<(), SolverError> {
        self.agree(
            comm,
            x.rows() != self.local_rows().len()
                || x.cols() > self.capacity
                || energies.len() != x.cols()
                || x.validate_finite().is_err(),
        )?;
        self.project(comm, x)?;
        for (band, e) in energies.iter_mut().enumerate() {
            *e = self
                .projection
                .column(band)
                .iter()
                .zip(self.coupled.column(band))
                .map(|(a, b)| (a.conj() * b).re)
                .sum();
        }
        self.agree(comm, energies.iter().any(|v| !v.is_finite()))
    }
}

impl dft_solvers::ref_davidson_complex::DistributedComplexOperator
    for DistributedComplexHamiltonian
{
    fn global_len(&self) -> usize {
        self.global_len()
    }
    fn local_rows(&self) -> &[usize] {
        self.local_rows()
    }
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.worker_pool())
    }
    fn apply_block<K: Collective>(
        &mut self,
        comm: &K,
        x: &ComplexMatrix,
        out: &mut ComplexMatrix,
    ) -> Result<(), SolverError> {
        self.apply_block(comm, x, out)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{density, hamiltonian::Hamiltonian};
    use dft_geometry::Cell;
    use dft_numerics::{FftPlan3d, linalg::norm};
    use dft_solvers::LinearOperator;
    use dft_units::Hartree;
    use ndarray::Array3;

    fn fixture(
        cell: &Cell,
        k: [f64; 3],
        projectors: bool,
    ) -> (PlaneWaveBasis, FftGridSpec, Hamiltonian) {
        let basis = PlaneWaveBasis::at_k(cell, Hartree(0.9), k).unwrap();
        let db = PlaneWaveBasis::new(cell, Hartree(3.6)).unwrap();
        let grid = FftGridSpec::new(&basis, &db, None).unwrap();
        let atoms = if projectors {
            [2, 1]
                .into_iter()
                .map(|count| {
                    let columns = (0..count)
                        .map(|j| {
                            basis
                                .momenta()
                                .iter()
                                .map(|g| {
                                    C::from_polar(
                                        0.2 / (1. + g.norm_squared),
                                        g.cartesian[0] * (j + 1) as f64 * 0.37
                                            + g.cartesian[1] * 0.23,
                                    )
                                })
                                .collect()
                        })
                        .collect();
                    let coupling = (0..count * count)
                        .map(|i| if i / count == i % count { 0.7 } else { -0.13 })
                        .collect();
                    (columns, coupling)
                })
                .collect()
        } else {
            Vec::new()
        };
        let nonlocal = NonlocalOperator::from_columns(basis.momenta().to_vec(), atoms, 2);
        let dims = grid.dims();
        let potential = Array3::from_shape_fn(dims, |(i, j, l)| {
            0.2 * (std::f64::consts::TAU * i as f64 / dims[0] as f64).cos()
                + 0.1 * (std::f64::consts::TAU * j as f64 / dims[1] as f64).sin()
                - 0.07 * (std::f64::consts::TAU * l as f64 / dims[2] as f64).cos()
        });
        let h = Hamiltonian::new_with_grid(&basis, grid.clone(), potential, nonlocal).unwrap();
        (basis, grid, h)
    }

    /// GK16: the column-pass diagonal equals the former row-outer loop bit for
    /// bit, including skipped zero couplings and a G = 0 row with T = +0.
    #[test]
    fn static_diagonal_matches_row_outer_loop_bitwise() {
        let rows = 37;
        let channels = [4usize, 1, 3];
        let p: usize = channels.iter().sum();
        let mut projectors = ComplexMatrix::zeros(rows, p).unwrap();
        for c in 0..p {
            for g in 0..rows {
                projectors[(g, c)] = if g == 0 && c % 2 == 0 {
                    C::new(0., -0.)
                } else {
                    C::new(
                        ((g * 5 + c * 3 + 1) as f64).sin() * 0.3,
                        ((g * 2 + c * 7) as f64).cos() * 0.2,
                    )
                };
            }
        }
        let couplings: Vec<Vec<f64>> = channels
            .iter()
            .enumerate()
            .map(|(a, &n)| {
                (0..n * n)
                    .map(|k| match (k + a) % 3 {
                        0 => 0.,
                        1 => -0.,
                        _ => 0.4 - 0.11 * k as f64,
                    })
                    .collect()
            })
            .collect();
        let kinetic: Vec<f64> = (0..rows).map(|g| 0.5 * (g * g) as f64 / 7.).collect();
        let mut reference = kinetic.clone();
        let mut offset = 0;
        for (&count, coupling) in channels.iter().zip(&couplings) {
            for (g, value) in reference.iter_mut().enumerate() {
                for i in 0..count {
                    for j in 0..count {
                        *value += (coupling[i * count + j]
                            * projectors[(g, offset + i)]
                            * projectors[(g, offset + j)].conj())
                        .re;
                    }
                }
            }
            offset += count;
        }
        let fast = static_diagonal(&kinetic, &projectors, &channels, &couplings);
        let bits = |v: &[f64]| v.iter().map(|x| x.to_bits()).collect::<Vec<_>>();
        assert_eq!(bits(&fast), bits(&reference));
    }

    fn waves(rows: usize, bands: usize) -> ComplexMatrix {
        let mut x = ComplexMatrix::zeros(rows, bands).unwrap();
        for j in 0..bands {
            for i in 0..rows {
                x[(i, j)] = C::new(
                    ((i * 7 + j * 13 + 1) as f64).sin(),
                    ((i * 11 + j * 3 + 2) as f64).cos(),
                );
            }
            let length = norm(x.column(j));
            for z in x.column_mut(j) {
                *z /= length;
            }
        }
        x
    }

    #[test]
    fn complex_owned_h_density_and_nonlocal_match_serial() {
        let mut worst_h = 0_f64;
        let mut worst_density = 0_f64;
        let mut worst_nl = 0_f64;
        for skew in [false, true] {
            let cell = Cell::from_user_rows_bohr(if skew {
                [[6., 0., 0.], [0.7, 5.5, 0.], [0.3, -0.2, 7.]]
            } else {
                [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]
            })
            .unwrap();
            for k in [[0.; 3], [0.25; 3], [-0.5, 0., -0.5], [0.13, -0.21, 0.07]] {
                for projectors in [false, true] {
                    let (basis, grid, h) = fixture(&cell, k, projectors);
                    for bands in [1, 3] {
                        let x = waves(basis.vectors().len(), bands);
                        let mut hx = ComplexMatrix::zeros(x.rows(), bands).unwrap();
                        h.apply_block(&x, &mut hx).unwrap();
                        let weights: Vec<_> =
                            (0..bands).map(|j| if j == 0 { 0.7 } else { 0. }).collect();
                        let rho = density::from_orbitals(
                            &(0..bands).map(|j| x.column(j).to_vec()).collect::<Vec<_>>(),
                            &weights,
                            &grid,
                            cell.volume_bohr3(),
                            &mut FftPlan3d::new(grid.dims()).unwrap(),
                        )
                        .unwrap();
                        let expected_nl: Vec<_> = (0..bands)
                            .map(|j| {
                                let mut y = vec![C::default(); x.rows()];
                                h.nonlocal().apply(x.column(j), &mut y).unwrap();
                                x.column(j)
                                    .iter()
                                    .zip(y)
                                    .map(|(a, b)| (a.conj() * b).re)
                                    .sum::<f64>()
                            })
                            .collect();
                        for ranks in [1, 3, 8] {
                            // Hamiltonian owns mutable scratch and is not Sync; only
                            // immutable input/oracle arrays cross the simulator boundary.
                            let outputs = dft_mpi::simulate(ranks, |comm| {
                                let (_, _, h) = fixture(&cell, k, projectors);
                                let mut op = DistributedComplexHamiltonian::new_owned(
                                    &basis,
                                    &grid,
                                    h.nonlocal(),
                                    bands,
                                    comm.rank(),
                                    ranks,
                                    WorkerPool::new(1).unwrap(),
                                )
                                .unwrap();
                                let offsets = op.plane_offsets();
                                let potential: Vec<_> = offsets
                                    .iter()
                                    .map(|&i| h.local_potential_values()[i])
                                    .collect();
                                op.set_owned_potential(&potential, 0).unwrap();
                                let mut local =
                                    ComplexMatrix::zeros(op.local_rows().len(), bands).unwrap();
                                for j in 0..bands {
                                    for (i, &g) in op.local_rows().iter().enumerate() {
                                        local[(i, j)] = x[(g, j)];
                                    }
                                }
                                let mut out = ComplexMatrix::zeros(local.rows(), bands).unwrap();
                                op.apply_block(comm, &local, &mut out).unwrap();
                                let dh = (0..bands)
                                    .flat_map(|j| {
                                        op.local_rows()
                                            .iter()
                                            .enumerate()
                                            .map(move |(i, &g)| (i, g, j))
                                    })
                                    .map(|(i, g, j)| (out[(i, j)] - hx[(g, j)]).norm())
                                    .fold(0., f64::max);
                                let mut density = vec![0.; offsets.len()];
                                op.band_density_sum(
                                    comm,
                                    &local,
                                    &weights,
                                    cell.volume_bohr3(),
                                    &mut density,
                                )
                                .unwrap();
                                let dn = offsets
                                    .iter()
                                    .zip(&density)
                                    .map(|(&i, n)| (n - rho.as_slice().unwrap()[i]).abs())
                                    .fold(0., f64::max);
                                let mut energies = vec![0.; bands];
                                op.nonlocal_band_energies(comm, &local, &mut energies)
                                    .unwrap();
                                let de = energies
                                    .iter()
                                    .zip(&expected_nl)
                                    .map(|(a, b)| (a - b).abs())
                                    .fold(0., f64::max);
                                // Reuse the buffers after their density and projection mutations.
                                let old = out.clone();
                                op.apply_block(comm, &local, &mut out).unwrap();
                                assert_eq!(old, out);
                                assert!(op.set_owned_potential(&potential, 0).is_err());
                                (dh, dn, de)
                            });
                            for (dh, dn, de) in outputs {
                                worst_h = worst_h.max(dh);
                                worst_density = worst_density.max(dn);
                                worst_nl = worst_nl.max(de);
                                assert!(dh < 1e-13, "H error {dh:e} Ha");
                                assert!(dn < 1e-13, "density error {dn:e} e/Bohr^3");
                                assert!(de < 1e-13, "nonlocal error {de:e} Ha");
                            }
                        }
                    }
                }
            }
        }
        eprintln!(
            "G3 complex operator: max H={worst_h:e} Ha, density={worst_density:e} e/Bohr^3, nonlocal={worst_nl:e} Ha"
        );
    }

    #[test]
    fn complex_empty_owners_delta_and_collective_rejection() {
        let cell = Cell::from_user_rows_bohr([[2., 0., 0.], [0.1, 2., 0.], [0., 0., 2.]]).unwrap();
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let (basis, grid, h) = fixture(&cell, [0.1, -0.05, 0.], true);
                let mut op = DistributedComplexHamiltonian::new_owned(
                    &basis,
                    &grid,
                    h.nonlocal(),
                    2,
                    comm.rank(),
                    ranks,
                    WorkerPool::new(1).unwrap(),
                )
                .unwrap();
                let mut x = ComplexMatrix::zeros(op.local_rows().len(), 2).unwrap();
                for j in 0..2 {
                    for z in x.column_mut(j) {
                        *z = C::new(0.3, 0.4);
                    }
                }
                let potential: Vec<_> = op
                    .plane_offsets()
                    .iter()
                    .map(|&i| h.local_potential_values()[i])
                    .collect();
                let mut hx = ComplexMatrix::zeros(x.rows(), 2).unwrap();
                // An unset potential on even an empty owner fails collectively.
                if comm.rank() + 1 != ranks {
                    op.set_owned_potential(&potential, 0).unwrap();
                }
                assert!(op.apply_block(comm, &x, &mut hx).is_err());
                if comm.rank() + 1 == ranks {
                    op.set_owned_potential(&potential, 0).unwrap();
                }
                op.apply_block(comm, &x, &mut hx).unwrap();
                let change: Vec<_> = potential.iter().map(|v| 0.2 - 0.13 * v).collect();
                let mut updated = hx.clone();
                op.add_local_change(comm, &x, &change, &mut updated)
                    .unwrap();
                let next: Vec<_> = potential.iter().zip(&change).map(|(v, d)| v + d).collect();
                op.set_owned_potential(&next, 1).unwrap();
                op.apply_block(comm, &x, &mut hx).unwrap();
                let delta = hx
                    .as_slice()
                    .iter()
                    .zip(updated.as_slice())
                    .map(|(a, b)| (*a - *b).norm())
                    .fold(0., f64::max);
                assert!(delta < 1e-13, "delta-HX error {delta:e} Ha");
                let empty = ComplexMatrix::zeros(x.rows(), 0).unwrap();
                let mut out = empty.clone();
                op.apply_block(comm, &empty, &mut out).unwrap();
                let mut density = vec![1.; op.plane_points()];
                op.band_density_sum(comm, &empty, &[], cell.volume_bohr3(), &mut density)
                    .unwrap();
                assert!(density.iter().all(|&n| n == 0.));
                // Invalid shape on one rank must be agreed before FFT work.
                let bad =
                    ComplexMatrix::zeros(x.rows() + usize::from(comm.rank() == 0), 2).unwrap();
                assert!(op.apply_block(comm, &bad, &mut hx).is_err());
            });
        }
    }
}
