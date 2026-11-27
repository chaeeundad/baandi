//! Transition blocks of one (k, k + q) pair: occupied orbitals at k (left)
//! against a block of empty orbitals at k + q (right), pair densities by
//! [`crate::pair`], optional q -> 0 slopes by [`crate::velocity`], added to a
//! [`crate::chi0::Accumulator`]. The unit of work the response driver
//! distributes over MPI ranks is one call of [`Transitions::add`].
use crate::{
    ResponseError,
    chi0::Accumulator,
    input, numerical,
    pair::{PairFft, Plan},
    velocity::velocity,
};
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::stress::ProjectorKGradient;
use num_complex::Complex64 as C;

/// Orbitals of one Bloch vector: rows in plane-wave order (`miller`),
/// columns are bands in ascending energy (Ha), every row present (a rank
/// holding a transition holds the whole vectors).
#[derive(Clone, Copy)]
pub struct KOrbitals<'a> {
    pub miller: &'a [[i32; 3]],
    pub coefficients: &'a ComplexMatrix,
    pub energies: &'a [f64],
}

/// Data of the velocity operator at k (optical limit only).
#[derive(Clone, Copy)]
pub struct Velocity<'a> {
    /// k+G (Cartesian, Bohr^-1) per row.
    pub momenta: &'a [[f64; 3]],
    /// KB projectors and k gradients on the same rows; `None` = local H.
    pub projectors: Option<&'a ProjectorKGradient>,
}

pub struct Transitions {
    left_fields: Vec<C>,
    grid: usize,
    left_energies: Vec<f64>,
    occupied: ComplexMatrix,
    right_plan: Plan,
    right_fields: Vec<C>,
    pairs: ComplexMatrix,
    delta: Vec<f64>,
    slopes: Vec<[C; 3]>,
}

impl Transitions {
    /// `occupied` lowest left bands; right orbitals at k' with
    /// `k + q = k' + shift` (`shift` in reciprocal-lattice units).
    pub fn new(
        dims: [usize; 3],
        left: KOrbitals<'_>,
        occupied: usize,
        right_miller: &[[i32; 3]],
        shift: [i32; 3],
    ) -> Result<Self, ResponseError> {
        let rows = left.miller.len();
        if left.coefficients.rows() != rows
            || occupied == 0
            || occupied > left.coefficients.cols()
            || left.energies.len() != left.coefficients.cols()
        {
            return Err(input("transitions: left orbital shapes"));
        }
        let mut left_plan = Plan::new(dims, left.miller, [0; 3])?;
        let mut left_fields = Vec::new();
        left_plan.to_real(left.coefficients, 0..occupied, &mut left_fields)?;
        let mut occ = ComplexMatrix::zeros(rows, occupied).map_err(numerical)?;
        for n in 0..occupied {
            occ.column_mut(n)
                .copy_from_slice(left.coefficients.column(n));
        }
        Ok(Self {
            grid: dims.iter().product(),
            left_fields,
            left_energies: left.energies[..occupied].to_vec(),
            occupied: occ,
            right_plan: Plan::new(dims, right_miller, shift)?,
            right_fields: Vec::new(),
            pairs: ComplexMatrix::zeros(0, 0).map_err(numerical)?,
            delta: Vec::new(),
            slopes: Vec::new(),
        })
    }

    /// Transitions (n occupied) x (m in `block` of `right`) with weight
    /// `scale` (= 2 w_k / Omega). `velocity` (optical limit, left = right = k)
    /// supplies the q-linear pair densities `a_nm = p_nm / (e_m - e_n)`.
    pub fn add(
        &mut self,
        accumulator: &mut Accumulator,
        fft: &mut PairFft,
        right: KOrbitals<'_>,
        block: std::ops::Range<usize>,
        scale: f64,
        velocity_data: Option<Velocity<'_>>,
    ) -> Result<(), ResponseError> {
        let b = self.block(fft, right, block, velocity_data)?;
        if b.delta.is_empty() {
            return Ok(());
        }
        accumulator.add(b.pairs, b.delta, scale, b.slopes)
    }

    /// The transition block without accumulation (R3: one block feeds the
    /// accumulators of several response cutoffs). Columns are `n * count + t`
    /// for occupied n and block member t; empty when `block` is empty.
    pub fn block(
        &mut self,
        fft: &mut PairFft,
        right: KOrbitals<'_>,
        block: std::ops::Range<usize>,
        velocity_data: Option<Velocity<'_>>,
    ) -> Result<TransitionBlock<'_>, ResponseError> {
        let count = block.len();
        let occ = self.left_energies.len();
        if right.coefficients.rows() != self.right_plan.len()
            || block.end > right.coefficients.cols()
            || right.energies.len() != right.coefficients.cols()
        {
            return Err(input("transitions: right orbital shapes"));
        }
        if count == 0 {
            self.delta.clear();
            return Ok(TransitionBlock {
                pairs: &self.pairs,
                delta: &self.delta,
                slopes: None,
            });
        }
        self.right_plan
            .to_real(right.coefficients, block.clone(), &mut self.right_fields)?;
        let n_g = fft.response_len();
        if self.pairs.rows() != n_g {
            self.pairs = ComplexMatrix::zeros(n_g, occ * count).map_err(numerical)?;
        } else {
            crate::set_columns(&mut self.pairs, occ * count)?;
        }
        self.delta.clear();
        for n in 0..occ {
            let left = &self.left_fields[n * self.grid..(n + 1) * self.grid];
            fft.pairs(left, &self.right_fields, count, &mut self.pairs, n * count)?;
            for m in block.clone() {
                self.delta.push(self.left_energies[n] - right.energies[m]);
            }
        }
        let slopes = if let Some(v) = velocity_data {
            let mut block_vectors =
                ComplexMatrix::zeros(right.coefficients.rows(), count).map_err(numerical)?;
            for (t, m) in block.clone().enumerate() {
                block_vectors
                    .column_mut(t)
                    .copy_from_slice(right.coefficients.column(m));
            }
            let p = velocity(v.momenta, v.projectors, &self.occupied, &block_vectors)?;
            self.slopes.clear();
            for n in 0..occ {
                for t in 0..count {
                    // rho_nm(q) = q . p_nm / (e_m - e_n), e_m - e_n = -delta.
                    let d = -self.delta[n * count + t];
                    self.slopes.push(std::array::from_fn(|b| p[b][(n, t)] / d));
                }
            }
            true
        } else {
            false
        };
        Ok(TransitionBlock {
            pairs: &self.pairs,
            delta: &self.delta,
            slopes: slopes.then_some(self.slopes.as_slice()),
        })
    }
}

/// One formed transition block: `pairs[(G, j)] = rho_j(q+G)` on the
/// response sphere of the [`PairFft`], `delta[j] = e_occ - e_emp` (Ha) and,
/// in the optical limit, the q-linear slopes (Bohr).
pub struct TransitionBlock<'a> {
    pub pairs: &'a ComplexMatrix,
    pub delta: &'a [f64],
    pub slopes: Option<&'a [[C; 3]]>,
}
