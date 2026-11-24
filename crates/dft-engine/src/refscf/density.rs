//! Valence rho(G) in e/bohr^3, including G=0=N_el/Omega. FFT forward carries
//! 1/N_r; inverse is the unscaled sum. No NLCC core, clipping or renormalization.
//! The half-sphere is i>0 || (i=0,j>0) || (i=j=0,k>=0). BAANDI's order is
//! retained; external entries are matched by Miller triple, never position.
//! History/metric use the smooth sphere; higher dense G use simple mixing.
use super::{Error, invalid, numerical};
use dft_basis::PlaneWaveBasis;
use dft_geometry::Cell;
use dft_mpi::Collective;
use dft_physics::data::{Array3, Complex64 as C, Hartree};
use dft_solvers::distributed_fft::DistributedPairFft;

pub fn half(m: [i32; 3]) -> bool {
    m[0] > 0 || (m[0] == 0 && (m[1] > 0 || (m[1] == 0 && m[2] >= 0)))
}
pub fn g_squared(b: [[f64; 3]; 3], m: [i32; 3]) -> f64 {
    (0..3)
        .map(|i| {
            (0..3)
                .map(|j| b[i][j] * f64::from(m[j]))
                .sum::<f64>()
                .powi(2)
        })
        .sum()
}
/// Metadata is replicated; all numerical density/history arrays are local G
/// rows or local real-space planes. Uses the same full-grid stick plan as
/// grid_dist, including its first-axis choice and local_positions ordering.
pub struct DensityPlan {
    fft: DistributedPairFft,
    pub miller: Vec<[i32; 3]>,
    pub g2: Vec<f64>,
    pub smooth: Vec<bool>,
    selected: Vec<usize>,
    coefficients: Vec<C>,
    offsets: Vec<usize>,
}
impl DensityPlan {
    pub fn new<K: Collective>(
        comm: &K,
        cell: &Cell,
        dims: [usize; 3],
        dense_cutoff_ha: f64,
        smooth_cutoff_ha: f64,
    ) -> Result<Self, Error> {
        Self::with_fft(comm, cell, dims, dense_cutoff_ha, smooth_cutoff_ha, None)
    }
    pub(crate) fn with_fft<K: Collective>(
        comm: &K,
        cell: &Cell,
        dims: [usize; 3],
        dense_cutoff_ha: f64,
        smooth_cutoff_ha: f64,
        shared: Option<DistributedPairFft>,
    ) -> Result<Self, Error> {
        if smooth_cutoff_ha <= 0.0 || smooth_cutoff_ha > dense_cutoff_ha {
            return Err(invalid("density cutoffs"));
        }
        let basis = PlaneWaveBasis::new(cell, Hartree(dense_cutoff_ha)).map_err(numerical)?;
        let fft = if let Some(fft) = shared {
            fft
        } else {
            let indices: Vec<_> = (0..dims[0])
                .flat_map(|i| (0..dims[1]).flat_map(move |j| (0..dims[2]).map(move |k| [i, j, k])))
                .collect();
            DistributedPairFft::new(dims, &indices, comm.rank(), comm.size()).map_err(numerical)?
        };
        let local: std::collections::HashMap<_, _> = fft
            .local_positions()
            .iter()
            .enumerate()
            .map(|(i, &g)| (g, i))
            .collect();
        let mut result = Self {
            coefficients: vec![C::default(); fft.local_len()],
            fft,
            miller: Vec::new(),
            g2: Vec::new(),
            smooth: Vec::new(),
            selected: Vec::new(),
            offsets: Vec::new(),
        };
        for g in basis.vectors().iter().filter(|g| half(g.miller)) {
            // Reject an aliased / missing sphere rather than wrapping a cutoff
            // beyond the grid's canonical frequency range.
            let ijk = g.miller.map(|m| m as i64);
            if (0..3).any(|d| ijk[d].unsigned_abs() as usize * 2 >= dims[d] && ijk[d] != 0) {
                return Err(invalid(
                    "density sphere reaches FFT Nyquist; enlarge fft_grid",
                ));
            }
            let idx: [usize; 3] =
                std::array::from_fn(|d| ijk[d].rem_euclid(dims[d] as i64) as usize);
            let offset = (idx[0] * dims[1] + idx[1]) * dims[2] + idx[2];
            if let Some(&p) = local.get(&offset) {
                result.miller.push(g.miller);
                result.g2.push(g.norm_squared);
                result.smooth.push(
                    g.norm_squared <= 2.0 * smooth_cutoff_ha * (1.0 + dft_basis::CUTOFF_TOLERANCE),
                );
                result.selected.push(p);
            }
        }
        for i in 0..result.fft.plane_count() {
            for y in 0..result.fft.rows() {
                let start = result.fft.grid_start(i, y);
                result.offsets.extend(start..start + dims[2]);
            }
        }
        Ok(result)
    }
    pub fn plane_len(&self) -> usize {
        self.offsets.len()
    }
    pub fn first_axis(&self) -> usize {
        self.fft.first_axis()
    }
    pub fn planes_from_full(&self, full: &Array3<f64>, planes: &mut [f64]) -> Result<(), Error> {
        if full.shape() != self.fft.shape() || planes.len() != self.offsets.len() {
            return Err(invalid("density grid shape"));
        }
        let full = full
            .as_slice()
            .ok_or_else(|| invalid("noncontiguous density"))?;
        for (v, &i) in planes.iter_mut().zip(&self.offsets) {
            *v = full[i];
        }
        Ok(())
    }
    pub fn forward<K: Collective>(
        &mut self,
        comm: &K,
        planes: &[f64],
        out: &mut [C],
    ) -> Result<(), Error> {
        if planes.len() != self.offsets.len() || out.len() != self.selected.len() {
            return Err(invalid("density forward shape"));
        }
        let size = self.fft.rows() * self.fft.shape()[2];
        for i in 0..self.fft.plane_count() {
            for (z, &v) in self
                .fft
                .plane_mut()
                .iter_mut()
                .zip(&planes[i * size..(i + 1) * size])
            {
                *z = C::new(v, 0.);
            }
            self.fft.to_reciprocal_plane(i).map_err(numerical)?;
        }
        self.fft
            .finish_to_reciprocal(comm, &mut self.coefficients)
            .map_err(numerical)?;
        for (o, &i) in out.iter_mut().zip(&self.selected) {
            *o = self.coefficients[i];
        }
        Ok(())
    }
    pub fn inverse<K: Collective>(
        &mut self,
        comm: &K,
        coefficients: &[C],
        planes: &mut [f64],
    ) -> Result<(), Error> {
        if coefficients.len() != self.selected.len() || planes.len() != self.offsets.len() {
            return Err(invalid("density inverse shape"));
        }
        self.coefficients.fill(C::default());
        let mut zero = [0.0];
        for ((&v, &i), &m) in coefficients.iter().zip(&self.selected).zip(&self.miller) {
            self.coefficients[i] = v;
            if m == [0; 3] {
                zero[0] = v.re;
            }
        }
        let mut g0 = [0.0];
        comm.tree_sum_f64(&zero, &mut g0).map_err(numerical)?;
        // Inverse of the positive half alone: 2 Re(sum_half)-rho(0).
        // No conjugate-grid replication or mirror communication is necessary.
        self.fft
            .begin_to_real(comm, &self.coefficients)
            .map_err(numerical)?;
        let size = self.fft.rows() * self.fft.shape()[2];
        for i in 0..self.fft.plane_count() {
            self.fft.to_real_plane(i).map_err(numerical)?;
            for (v, z) in planes[i * size..(i + 1) * size]
                .iter_mut()
                .zip(self.fft.plane())
            {
                *v = 2.0 * z.re - g0[0];
            }
        }
        Ok(())
    }
    pub fn write_serial_full(&self, planes: &[f64], full: &mut Array3<f64>) -> Result<(), Error> {
        if self.offsets.len() != self.fft.shape().iter().product::<usize>()
            || planes.len() != self.offsets.len()
            || full.shape() != self.fft.shape()
        {
            return Err(invalid("full density requested from distributed planes"));
        }
        let values = full
            .as_slice_mut()
            .ok_or_else(|| invalid("noncontiguous full density"))?;
        for (&v, &i) in planes.iter().zip(&self.offsets) {
            values[i] = v;
        }
        Ok(())
    }
}
