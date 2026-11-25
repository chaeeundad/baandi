//! Read-only finite-angle orbital-rotation density probe. Band selection is
//! explicit diagnostic input, never inferred to be a physical zero mode.
use super::{density::DensityPlan, electronic_map::Map, metric::Metric, numerical};
use crate::{ScfError, grid_dist::PlaneMap, refstart::KStartState};
use dft_mpi::Collective;
use dft_numerics::{
    linalg::ComplexMatrix,
    real::{RealEighWorkspace, RealMatrix},
};
use dft_physics::data::Complex64 as C;

pub(super) struct Probe {
    pairs: Vec<[usize; 4]>,
    angle: f64,
    volume: f64,
    orbitals: Vec<ComplexMatrix>,
    planes: Vec<PlaneMap>,
    plus: Vec<Vec<f64>>,
    minus: Vec<Vec<f64>>,
    local: Vec<f64>,
    tangent: Vec<Vec<C>>,
    gram: RealMatrix,
    vectors: RealMatrix,
    eigenvalues: Vec<f64>,
    rhs: Vec<f64>,
    solver: RealEighWorkspace,
}
impl Probe {
    pub fn new<K: Collective>(
        pairs: &[[usize; 4]],
        angle: f64,
        states: &[KStartState],
        plan: &DensityPlan,
        comm: &K,
        volume: f64,
    ) -> Result<Self, ScfError> {
        if pairs.is_empty() || !angle.is_finite() || angle <= 0. || angle >= 0.1 {
            return Err(numerical("rotation probe needs pairs and 0 < angle < 0.1 radians").into());
        }
        let mut orbitals = Vec::new();
        let mut planes = Vec::new();
        let mut plus = Vec::new();
        let mut minus = Vec::new();
        for &[s, k, i, j] in pairs {
            let state = states
                .get(s)
                .ok_or_else(|| numerical("rotation probe channel"))?;
            let block = state
                .blocks
                .get(k)
                .ok_or_else(|| numerical("rotation probe k point"))?;
            if i == j || i.max(j) >= state.orbitals.0[k].cols() {
                return Err(numerical("rotation probe band pair").into());
            }
            orbitals.push(ComplexMatrix::zeros(state.orbitals.0[k].rows(), 2).map_err(numerical)?);
            planes.push(PlaneMap::new(
                state.grid.dims(),
                comm.rank(),
                comm.size(),
                state.fields.first_axis(),
                block.operator.first_axis(),
            ));
            plus.push(vec![0.; block.operator.plane_points()]);
            minus.push(vec![0.; block.operator.plane_points()]);
        }
        let n = 2 * pairs.len();
        Ok(Self {
            pairs: pairs.to_vec(),
            angle,
            volume,
            orbitals,
            planes,
            plus,
            minus,
            local: vec![0.; states[0].density.len()],
            tangent: vec![vec![C::default(); states.len() * plan.miller.len()]; n],
            gram: RealMatrix::zeros(n, n).map_err(numerical)?,
            vectors: RealMatrix::zeros(n, n).map_err(numerical)?,
            eigenvalues: vec![0.; n],
            rhs: vec![0.; n],
            solver: RealEighWorkspace::new(n).map_err(numerical)?,
        })
    }
    pub fn measure<K: Collective>(
        &mut self,
        states: &mut [KStartState],
        step: &Map,
        plan: &mut DensityPlan,
        metric: &mut Metric,
        comm: &K,
        residual: &[C],
    ) -> Result<(f64, usize), ScfError> {
        let size = plan.miller.len();
        let channels = states.len();
        for (p, &[s, k, i, j]) in self.pairs.iter().enumerate() {
            let state = &mut states[s];
            let block = &step.blocks[s][k];
            for quadrature in 0..2 {
                let phase = if quadrature == 0 {
                    C::new(1., 0.)
                } else {
                    C::new(0., 1.)
                };
                for r in 0..block.orbitals.rows() {
                    let a = block.orbitals[(r, i)] * self.angle.cos();
                    let b = block.orbitals[(r, j)] * phase * self.angle.sin();
                    self.orbitals[p][(r, 0)] = a + b;
                    self.orbitals[p][(r, 1)] = a - b;
                }
                let op = &mut state.blocks[k].operator;
                op.band_density_sum(
                    comm,
                    &self.orbitals[p],
                    &[1., 0.],
                    self.volume,
                    &mut self.plus[p],
                )
                .map_err(numerical)?;
                op.band_density_sum(
                    comm,
                    &self.orbitals[p],
                    &[0., 1.],
                    self.volume,
                    &mut self.minus[p],
                )
                .map_err(numerical)?;
                let weight = state.blocks[k].point.weight
                    * (block.occupations[i] - block.occupations[j])
                    / (2. * self.angle);
                for (a, b) in self.plus[p].iter_mut().zip(&self.minus[p]) {
                    *a = (*a - b) * weight;
                }
                self.planes[p].b_to_a(comm, &self.plus[p], &mut self.local)?;
                let t = &mut self.tangent[2 * p + quadrature];
                plan.forward(comm, &self.local, &mut t[..size])?;
                for (g, m) in plan.miller.iter().enumerate() {
                    if *m == [0, 0, 0] {
                        t[g] = C::default();
                    }
                }
                if channels == 2 {
                    for g in 0..size {
                        t[size + g] = t[g] * if s == 0 { 1. } else { -1. };
                    }
                }
            }
        }
        for i in 0..self.tangent.len() {
            let mut value = [0.];
            metric.batch(comm, &[(&self.tangent[i], residual)], &mut value)?;
            self.rhs[i] = value[0];
            for j in 0..=i {
                metric.batch(comm, &[(&self.tangent[i], &self.tangent[j])], &mut value)?;
                self.gram[(i, j)] = value[0];
                self.gram[(j, i)] = value[0];
            }
        }
        self.solver
            .solve_into(&self.gram, &mut self.eigenvalues, &mut self.vectors)
            .map_err(numerical)?;
        let cutoff = 64. * f64::EPSILON * self.eigenvalues.last().copied().unwrap_or(0.);
        let mut projected = 0.;
        let mut rank = 0;
        for (j, &lambda) in self.eigenvalues.iter().enumerate() {
            if lambda > cutoff {
                let a = (0..self.rhs.len())
                    .map(|i| self.vectors[(i, j)] * self.rhs[i])
                    .sum::<f64>();
                projected += a * a / lambda;
                rank += 1;
            }
        }
        let norm = metric.norm(comm, residual)?;
        Ok((
            if norm > 0. {
                (projected / norm).clamp(0., 1.)
            } else {
                0.
            },
            rank,
        ))
    }
}
