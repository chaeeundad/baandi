//! Serial SpinXc is the oracle; simulated ranks run the production FFT transposes.
use dft_geometry::Cell;
use dft_mpi::Collective;
use dft_physics::{
    PhysicsError,
    data::{Array3, Complex64 as C},
    gga::distributed::{Fields, SpinPbe},
    spin::{SpinFunctional, SpinXc, SpinXcResponse},
};
use dft_solvers::distributed_fft::DistributedPairFft;

struct Slabs(DistributedPairFft);
impl Slabs {
    fn new<K: Collective>(dims: [usize; 3], comm: &K) -> Self {
        let indices: Vec<_> = (0..dims[0])
            .flat_map(|i| (0..dims[1]).flat_map(move |j| (0..dims[2]).map(move |k| [i, j, k])))
            .collect();
        Self(DistributedPairFft::new(dims, &indices, comm.rank(), comm.size()).unwrap())
    }
    fn offsets(&self) -> Vec<usize> {
        (0..self.0.plane_count())
            .flat_map(|i| {
                (0..self.0.rows()).flat_map(move |j| {
                    let start = self.0.grid_start(i, j);
                    start..start + self.0.shape()[2]
                })
            })
            .collect()
    }
}
impl<K: Collective> Fields<K> for Slabs {
    type Error = PhysicsError;
    fn dims(&self) -> [usize; 3] {
        self.0.shape()
    }
    fn reciprocal_indices(&self) -> &[usize] {
        self.0.local_positions()
    }
    fn real_indices(&self) -> Option<Vec<usize>> {
        Some(self.offsets())
    }
    fn plane_points(&self) -> usize {
        self.0.plane_points()
    }
    fn forward(&mut self, comm: &K, n: &[f64]) -> Result<Vec<C>, PhysicsError> {
        let len = self.0.rows() * self.0.shape()[2];
        for i in 0..self.0.plane_count() {
            for (z, &n) in self
                .0
                .plane_mut()
                .iter_mut()
                .zip(&n[i * len..(i + 1) * len])
            {
                *z = C::new(n, 0.);
            }
            self.0.to_reciprocal_plane(i).unwrap();
        }
        let mut out = vec![C::default(); self.0.local_len()];
        self.0.finish_to_reciprocal(comm, &mut out).unwrap();
        Ok(out)
    }
    fn inverse(&mut self, comm: &K, n: &[C]) -> Result<Vec<f64>, PhysicsError> {
        self.0.begin_to_real(comm, n).unwrap();
        let mut out = Vec::new();
        for i in 0..self.0.plane_count() {
            self.0.to_real_plane(i).unwrap();
            out.extend(self.0.plane().iter().map(|z| z.re));
        }
        Ok(out)
    }
}

#[test]
fn screened_spin_pbe_matches_serial_at_one_three_eight_ranks() {
    use dft_numerics::switches::RunConfig;
    use dft_physics::xc;
    let dims = [4, 3, 5];
    let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0.3, 7., 0.], [0.2, 0.1, 6.]]).unwrap();
    let mut density = [Array3::from_elem(dims, 0.06), Array3::from_elem(dims, 0.04)];
    let mut core = Array3::from_elem(dims, 0.01);
    density[0][[1, 2, 3]] = 1e-8;
    density[1][[1, 2, 3]] = 2e-8;
    core[[1, 2, 3]] = -4.2e-7;
    let config = RunConfig::from_pairs(&[("xc_screen_negative", "on")]).unwrap();
    let _s = config.enter();
    let mut serial = SpinXc::new(SpinFunctional::Pbe, &cell, dims).unwrap();
    let mut response = SpinXcResponse::new(dims);
    serial
        .evaluate_with_core(&density, Some(&core), &mut response)
        .unwrap();
    let mut reference_records = None;
    for ranks in [1, 3, 8] {
        let runs = dft_mpi::simulate(ranks, |comm| {
            let _s = config.enter();
            xc::clear_screening_records();
            let mut fields = Slabs::new(dims, comm);
            let offsets = fields.offsets();
            let local: [Vec<_>; 2] = std::array::from_fn(|s| {
                offsets
                    .iter()
                    .map(|&i| density[s].as_slice().unwrap()[i])
                    .collect()
            });
            let local_core: Vec<_> = offsets
                .iter()
                .map(|&i| core.as_slice().unwrap()[i])
                .collect();
            let mut xc = SpinPbe::new::<dft_mpi::SimRank, Slabs>(&fields, &cell).unwrap();
            let out = xc
                .evaluate(&mut fields, comm, [&local[0], &local[1]], Some(&local_core))
                .unwrap();
            let mut energy = [0.];
            comm.tree_sum_f64(&[out.energy], &mut energy).unwrap();
            assert!((energy[0] - response.energy_ha).abs() < 1e-12);
            for s in 0..2 {
                for (j, &i) in offsets.iter().enumerate() {
                    assert!(
                        (out.potentials[s][j] - response.potentials[s].as_slice().unwrap()[i])
                            .abs()
                            < 1e-12
                    );
                }
            }
            xc::screening_records()
        });
        for records in runs {
            if let Some(expected) = &reference_records {
                assert_eq!(&records, expected);
            } else {
                reference_records = Some(records);
            }
        }
    }
}
