//! Allocation-free application of a prepared density group projector.
use crate::PhysicsError;
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::symmetry::{SpinAction, SymmetryOperation};
use dft_numerics::Fft3d;
use ndarray::Array3;
use num_complex::Complex64;
use std::{collections::BTreeMap, f64::consts::TAU};

pub struct DensitySymmetrizer {
    /// For each output G: input R^T G and exp(-2*pi*i*G.tau)/|group|.
    terms: Vec<(usize, usize, Complex64)>,
    real: Array3<Complex64>,
    input: Array3<Complex64>,
    output: Array3<Complex64>,
}
impl DensitySymmetrizer {
    /// The density sphere must be closed under every operation. Grid padding
    /// outside that sphere is zero; no modulo alias of a rotated G is allowed.
    pub fn new(
        basis: &PlaneWaveBasis,
        grid: &FftGridSpec,
        operations: &[SymmetryOperation],
    ) -> Result<Self, PhysicsError> {
        if operations.is_empty()
            || operations
                .iter()
                .any(|o| o.spin_action != SpinAction::Identity || o.antiunitary)
        {
            return Err(PhysicsError::Shape);
        }
        let dims = grid.dims();
        let index =
            |m| FftGridSpec::map_miller(m, dims).map(|p| (p[0] * dims[1] + p[1]) * dims[2] + p[2]);
        let lookup: BTreeMap<_, _> = basis
            .vectors()
            .iter()
            .map(|g| Ok((g.miller, index(g.miller)?)))
            .collect::<Result<_, dft_basis::BasisError>>()?;
        let mut terms = Vec::new();
        for g in basis.vectors() {
            for op in operations {
                let source =
                    std::array::from_fn(|i| (0..3).map(|j| op.rotation[j][i] * g.miller[j]).sum());
                let src = *lookup.get(&source).ok_or(PhysicsError::Shape)?;
                let phase = -TAU
                    * (0..3)
                        .map(|i| g.miller[i] as f64 * op.translation[i])
                        .sum::<f64>();
                terms.push((
                    lookup[&g.miller],
                    src,
                    Complex64::from_polar(1.0 / operations.len() as f64, phase),
                ));
            }
        }
        Ok(Self {
            terms,
            real: Array3::zeros(dims),
            input: Array3::zeros(dims),
            output: Array3::zeros(dims),
        })
    }
    /// Fail closed for externally supplied densities outside the invariant
    /// subspace. Absolute allowance is 1e-12 electrons/Bohr^3 for FFT roundoff.
    pub fn check_invariant(
        &mut self,
        density: &Array3<f64>,
        fft: &mut impl Fft3d,
    ) -> Result<bool, PhysicsError> {
        self.transform(density, fft)?;
        Ok(density
            .iter()
            .zip(&self.real)
            .all(|(n, z)| n.is_finite() && (n - z.re).abs() <= 1e-12 && z.im.abs() <= 1e-12))
    }
    /// Each collinear channel uses the same unitary identity-spin projector.
    pub fn apply(
        &mut self,
        density: &mut Array3<f64>,
        fft: &mut impl Fft3d,
    ) -> Result<(), PhysicsError> {
        self.transform(density, fft)?;
        for (n, z) in density.iter_mut().zip(&self.real) {
            *n = z.re;
        }
        Ok(())
    }
    fn transform(
        &mut self,
        density: &Array3<f64>,
        fft: &mut impl Fft3d,
    ) -> Result<(), PhysicsError> {
        if density.shape() != self.real.shape() {
            return Err(PhysicsError::Shape);
        }
        for (z, &n) in self.real.iter_mut().zip(density.iter()) {
            *z = Complex64::new(n, 0.0);
        }
        fft.to_reciprocal(&self.real, &mut self.input)?;
        self.output.fill(Complex64::default());
        let input = self.input.as_slice().unwrap();
        let output = self.output.as_slice_mut().unwrap();
        for &(dst, src, phase) in &self.terms {
            output[dst] += phase * input[src];
        }
        fft.to_real(&self.output, &mut self.real)?;
        if self
            .real
            .iter()
            .any(|z| !z.re.is_finite() || z.im.abs() > 1e-12 * z.re.abs().max(1.0))
        {
            return Err(PhysicsError::Unrepresentable);
        }
        Ok(())
    }
}

/// Prepared full-stick projector. Each destination requests a source coefficient
/// once, even when several operations need it. Channels share one padded exchange.
/// Terms retain the serial G/operation order; no value array is replicated.
pub struct DistributedDensitySymmetrizer {
    local_indices: Vec<usize>,
    send_terms: Vec<(usize, usize)>,
    terms: Vec<(usize, usize, Complex64)>,
    send: Vec<f64>,
    recv: Vec<f64>,
    output: Vec<Complex64>,
    block: usize,
    channels: usize,
    rank: usize,
    ranks: usize,
}
impl DistributedDensitySymmetrizer {
    pub fn new(
        basis: &PlaneWaveBasis,
        dims: [usize; 3],
        operations: &[SymmetryOperation],
        rank: usize,
        ranks: usize,
        channels: usize,
    ) -> Result<Self, PhysicsError> {
        use std::collections::BTreeSet;
        if ranks == 0
            || dims.contains(&0)
            || rank >= ranks
            || channels == 0
            || operations.is_empty()
            || operations
                .iter()
                .any(|o| o.spin_action != SpinAction::Identity || o.antiunitary)
            // GS1D-B (ADR 0049): with GRID_INCOMPATIBLE_OPS=reciprocal a
            // rotation need not permute the grid (the projection is in G space
            // on the density sphere); translations must still be commensurate.
            || !operations.iter().all(|o| {
                dft_geometry::symmetry::grid_keeps(
                    dims,
                    o,
                    dft_numerics::switches::choice(
                        dft_numerics::switches::Id::GRID_INCOMPATIBLE_OPS,
                    ) == "reciprocal",
                )
            })
        {
            return Err(PhysicsError::Shape);
        }
        let positions: Vec<_> = (0..dims[0])
            .flat_map(|x| (0..dims[1]).flat_map(move |y| (0..dims[2]).map(move |z| [x, y, z])))
            .collect();
        let owners = dft_solvers::distributed_fft::position_owners(dims, &positions, ranks)
            .map_err(|_| PhysicsError::Shape)?;
        let index =
            |m| FftGridSpec::map_miller(m, dims).map(|p| (p[0] * dims[1] + p[1]) * dims[2] + p[2]);
        let lookup: BTreeMap<_, _> = basis
            .vectors()
            .iter()
            .map(|g| Ok((g.miller, index(g.miller)?)))
            .collect::<Result<_, dft_basis::BasisError>>()?;
        let mut requests = vec![vec![BTreeSet::new(); ranks]; ranks];
        let mut raw = Vec::new();
        for g in basis.vectors() {
            let dst = lookup[&g.miller];
            for op in operations {
                let source =
                    std::array::from_fn(|i| (0..3).map(|j| op.rotation[j][i] * g.miller[j]).sum());
                let src = *lookup.get(&source).ok_or(PhysicsError::Shape)?;
                requests[owners[dst]][owners[src]].insert(src);
                if owners[dst] == rank {
                    let phase = -TAU
                        * (0..3)
                            .map(|i| g.miller[i] as f64 * op.translation[i])
                            .sum::<f64>();
                    raw.push((
                        dst,
                        src,
                        Complex64::from_polar(1. / operations.len() as f64, phase),
                    ));
                }
            }
        }
        let width = requests
            .iter()
            .flatten()
            .map(BTreeSet::len)
            .max()
            .unwrap_or(0);
        let block = 2 * channels * width;
        let local_indices: Vec<_> = owners
            .iter()
            .enumerate()
            .filter_map(|(i, &o)| (o == rank).then_some(i))
            .collect();
        let local: BTreeMap<_, _> = local_indices
            .iter()
            .enumerate()
            .map(|(i, &g)| (g, i))
            .collect();
        let mut send_terms = Vec::new();
        let mut received = BTreeMap::new();
        for (peer, outgoing) in requests.iter().enumerate() {
            for (slot, &g) in outgoing[rank].iter().enumerate() {
                send_terms.push((peer * block + 2 * channels * slot, local[&g]));
            }
            for (slot, &g) in requests[rank][peer].iter().enumerate() {
                received.insert(g, peer * block + 2 * channels * slot);
            }
        }
        let terms = raw
            .into_iter()
            .map(|(dst, src, phase)| (local[&dst], received[&src], phase))
            .collect();
        Ok(Self {
            output: vec![Complex64::default(); channels * local_indices.len()],
            local_indices,
            send_terms,
            terms,
            send: vec![0.; ranks * block],
            recv: vec![0.; ranks * block],
            block,
            channels,
            rank,
            ranks,
        })
    }
    pub fn local_indices(&self) -> &[usize] {
        &self.local_indices
    }
    /// Forward-normalized coefficients. Returns the local squared Fourier L2
    /// difference summed over channels; multiply the global sum by volume for
    /// the real-space integral. Padding outside the density sphere is zeroed.
    pub fn apply<K: dft_mpi::Collective>(
        &mut self,
        comm: &K,
        density: &mut [Vec<Complex64>],
    ) -> Result<f64, PhysicsError> {
        let n = self.local_indices.len();
        let invalid = comm.rank() != self.rank
            || comm.size() != self.ranks
            || density.len() != self.channels
            || density
                .iter()
                .any(|v| v.len() != n || v.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()));
        let mut agreed = [0.];
        comm.tree_sum_f64(&[f64::from(invalid)], &mut agreed)
            .map_err(|_| PhysicsError::Unrepresentable)?;
        if agreed[0] != 0. {
            return Err(PhysicsError::Shape);
        }
        for &(slot, src) in &self.send_terms {
            for (c, values) in density.iter().enumerate() {
                self.send[slot + 2 * c] = values[src].re;
                self.send[slot + 2 * c + 1] = values[src].im;
            }
        }
        comm.alltoall_f64(&self.send, self.block, &mut self.recv)
            .map_err(|_| PhysicsError::Unrepresentable)?;
        self.output.fill(Complex64::default());
        for &(dst, slot, phase) in &self.terms {
            for c in 0..self.channels {
                self.output[c * n + dst] +=
                    phase * Complex64::new(self.recv[slot + 2 * c], self.recv[slot + 2 * c + 1]);
            }
        }
        let mut residual = 0.;
        for (c, values) in density.iter_mut().enumerate() {
            for (i, z) in values.iter_mut().enumerate() {
                let projected = self.output[c * n + i];
                residual += (*z - projected).norm_sqr();
                *z = projected;
            }
        }
        Ok(residual)
    }
}

#[cfg(test)]
mod distributed_tests {
    use super::*;
    use dft_geometry::{Atom, Cell};
    use dft_mpi::Collective;
    use dft_numerics::FftPlan3d;
    use dft_units::Hartree;

    #[test]
    fn symmetry_distributed_translations_channels_and_empty_owners() {
        let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
        let ops =
            dft_geometry::symmetry::detect(&cell, &[Atom::new("X", [0.25; 3])], 1e-7).unwrap();
        assert!(ops.iter().any(|o| o.translation.iter().any(|t| *t != 0.)));
        let basis = PlaneWaveBasis::new(&cell, Hartree(2.)).unwrap();
        let grid = FftGridSpec::new(&basis, &basis, Some([8; 3])).unwrap();
        let dims = grid.dims();
        let initial = Array3::from_shape_fn(dims, |(i, j, k)| {
            Complex64::new(
                2. + 0.2 * (TAU * (i + j) as f64 / 8.).cos() + 0.1 * (TAU * k as f64 / 8.).sin(),
                0.,
            )
        });
        let mut fft = FftPlan3d::new(dims).unwrap();
        let mut reciprocal = Array3::zeros(dims);
        fft.to_reciprocal(&initial, &mut reciprocal).unwrap();
        let mut expected = initial.mapv(|z| z.re);
        DensitySymmetrizer::new(&basis, &grid, &ops)
            .unwrap()
            .apply(&mut expected, &mut fft)
            .unwrap();
        for ranks in [1, 3, 8, 70] {
            dft_mpi::simulate(ranks, |comm| {
                let mut plan = DistributedDensitySymmetrizer::new(
                    &basis,
                    dims,
                    &ops,
                    comm.rank(),
                    comm.size(),
                    2,
                )
                .unwrap();
                let indices = plan.local_indices().to_vec();
                let mut channels = vec![
                    indices
                        .iter()
                        .map(|&g| reciprocal.as_slice().unwrap()[g])
                        .collect::<Vec<_>>(),
                    indices
                        .iter()
                        .map(|&g| 0.3 * reciprocal.as_slice().unwrap()[g])
                        .collect(),
                ];
                let residual = plan.apply(comm, &mut channels).unwrap();
                assert!(residual.is_finite());
                let once = channels.clone();
                plan.apply(comm, &mut channels).unwrap();
                assert!(
                    once.iter()
                        .flatten()
                        .zip(channels.iter().flatten())
                        .all(|(a, b)| (*a - *b).norm() < 1e-14)
                );
                let mut flat = vec![0.; 2 * expected.len()];
                for (i, &g) in indices.iter().enumerate() {
                    flat[2 * g] = once[0][i].re;
                    flat[2 * g + 1] = once[0][i].im;
                }
                let mut all = flat.clone();
                comm.tree_sum_f64(&flat, &mut all).unwrap();
                let array = Array3::from_shape_vec(
                    dims,
                    all.as_chunks::<2>()
                        .0
                        .iter()
                        .map(|z| Complex64::new(z[0], z[1]))
                        .collect(),
                )
                .unwrap();
                let mut real = Array3::zeros(dims);
                FftPlan3d::new(dims)
                    .unwrap()
                    .to_real(&array, &mut real)
                    .unwrap();
                assert!(
                    real.iter()
                        .zip(&expected)
                        .all(|(a, b)| (a.re - b).abs() < 1e-13 && a.im.abs() < 1e-13)
                );
                assert!(
                    once[0]
                        .iter()
                        .zip(&once[1])
                        .all(|(a, b)| (*b - 0.3 * a).norm() < 1e-14)
                );
            });
        }
    }
}
