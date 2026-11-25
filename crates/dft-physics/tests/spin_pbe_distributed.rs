//! Serial SpinXc is the oracle; simulated ranks run the production FFT transposes.
use dft_geometry::Cell;
use dft_mpi::{Collective, MpiError};
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
/// No payload allgathers are allowed inside the evaluator or FFT.
struct NoGather<'a, K>(&'a K);
impl<K: Collective> Collective for NoGather<'_, K> {
    fn rank(&self) -> usize {
        self.0.rank()
    }
    fn size(&self) -> usize {
        self.0.size()
    }
    fn barrier(&self) -> Result<(), MpiError> {
        self.0.barrier()
    }
    fn tree_sum_f64(&self, a: &[f64], b: &mut [f64]) -> Result<(), MpiError> {
        self.0.tree_sum_f64(a, b)
    }
    fn bcast_f64(&self, a: &mut [f64], root: usize) -> Result<(), MpiError> {
        self.0.bcast_f64(a, root)
    }
    fn allgatherv_f64(&self, _: &[f64], _: &[usize], _: &mut [f64]) -> Result<(), MpiError> {
        panic!("field gather")
    }
    fn alltoall_f64(&self, a: &[f64], n: usize, b: &mut [f64]) -> Result<(), MpiError> {
        self.0.alltoall_f64(a, n, b)
    }
}

#[test]
fn spin_pbe_matches_serial_with_core_and_empty_owners() {
    let cell = Cell::from_user_rows_bohr([[7., 0., 0.], [0.7, 5., 0.], [0.3, 0.2, 9.]]).unwrap();
    for dims in [[4, 6, 5], [5, 3, 7]] {
        let mut seed = 1301_u64;
        let mut random = || {
            seed = seed.wrapping_mul(6364136223846793005).wrapping_add(1);
            (seed >> 11) as f64 / ((1_u64 << 53) as f64)
        };
        let density: [Array3<f64>; 2] = std::array::from_fn(|s| {
            Array3::from_shape_fn(dims, |_| 0.02 + 0.01 * (s + 1) as f64 * random())
        });
        // Includes signed core samples, permitted whenever combined channels are positive.
        let core = Array3::from_shape_fn(dims, |_| 0.02 * (random() - 0.2));
        for with_core in [false, true] {
            let core = with_core.then_some(&core);
            for ranks in [1, 3, 8] {
                let errors = dft_mpi::simulate(ranks, |rank| {
                    let comm = NoGather(rank);
                    let mut fields = Slabs::new(dims, &comm);
                    let offsets = fields.offsets();
                    let local: [Vec<_>; 2] = std::array::from_fn(|s| {
                        offsets
                            .iter()
                            .map(|&i| density[s].as_slice().unwrap()[i])
                            .collect()
                    });
                    let local_core: Option<Vec<_>> =
                        core.map(|c| offsets.iter().map(|&i| c.as_slice().unwrap()[i]).collect());
                    let mut work =
                        SpinPbe::new::<NoGather<'_, dft_mpi::SimRank>, _>(&fields, &cell).unwrap();
                    let mut error = 0_f64;
                    let mut pointers = None;
                    for repeat in 0..2 {
                        // Reuse with a second density, catching stale gradient/flux state.
                        let scaled: [Vec<_>; 2] = std::array::from_fn(|s| {
                            local[s]
                                .iter()
                                .map(|n| n * (1. + 0.1 * repeat as f64))
                                .collect()
                        });
                        let actual = work
                            .evaluate(
                                &mut fields,
                                &comm,
                                [&scaled[0], &scaled[1]],
                                local_core.as_deref(),
                            )
                            .unwrap();
                        let now = actual.potentials.each_ref().map(|p| p.as_ptr());
                        if let Some(old) = pointers {
                            assert_eq!(old, now);
                        }
                        pointers = Some(now);
                        // Comparison below is filled against an independent serial response.
                        let mut reference = SpinXcResponse::new(dims);
                        let n = density
                            .each_ref()
                            .map(|n| n.mapv(|v| v * (1. + 0.1 * repeat as f64)));
                        SpinXc::new(SpinFunctional::Pbe, &cell, dims)
                            .unwrap()
                            .evaluate_with_core(&n, core, &mut reference)
                            .unwrap();
                        for s in 0..2 {
                            for (&v, &i) in actual.potentials[s].iter().zip(&offsets) {
                                let r = reference.potentials[s].as_slice().unwrap()[i];
                                error = error.max((v - r).abs() / r.abs());
                            }
                        }
                        let mut totals = [0.; 3];
                        comm.tree_sum_f64(
                            &[
                                actual.energy,
                                actual.valence_potential_integral,
                                actual.core_potential_integral,
                            ],
                            &mut totals,
                        )
                        .unwrap();
                        for (a, b) in totals.into_iter().zip([
                            reference.energy_ha,
                            reference.valence_potential_integral_ha,
                            reference.core_potential_integral_ha,
                        ]) {
                            if b != 0. {
                                error = error.max((a - b).abs() / b.abs());
                            } else {
                                assert_eq!(a, b);
                            }
                        }
                    }
                    assert!(error < 1e-13, "relative error {error:e}");
                    error
                });
                eprintln!(
                    "G4C grid={dims:?} core={with_core} ranks={ranks} max_relative={:e}",
                    errors.into_iter().fold(0., f64::max)
                );
            }
        }
    }
}

#[test]
fn spin_pbe_screening_polarization_and_collective_rejection() {
    let dims = [3, 4, 2];
    let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 6., 0.], [0., 0.1, 7.]]).unwrap();
    for channels in [[0., 0.], [4e-11, 5e-11], [0.03, 0.], [0., 0.03]] {
        let density = channels.map(|n| Array3::from_elem(dims, n));
        let mut reference = SpinXcResponse::new(dims);
        SpinXc::new(SpinFunctional::Pbe, &cell, dims)
            .unwrap()
            .evaluate(&density, &mut reference)
            .unwrap();
        dft_mpi::simulate(8, |rank| {
            let comm = NoGather(rank);
            let mut fields = Slabs::new(dims, &comm);
            let offsets = fields.offsets();
            let local = channels.map(|n| vec![n; offsets.len()]);
            let mut work =
                SpinPbe::new::<NoGather<'_, dft_mpi::SimRank>, _>(&fields, &cell).unwrap();
            let out = work
                .evaluate(&mut fields, &comm, [&local[0], &local[1]], None)
                .unwrap();
            for s in 0..2 {
                for (&v, &i) in out.potentials[s].iter().zip(&offsets) {
                    assert!((v - reference.potentials[s].as_slice().unwrap()[i]).abs() < 1e-14);
                }
            }
            for bad in [-1e-20, f64::NAN, f64::INFINITY] {
                let mut invalid = local.clone();
                if offsets.contains(&0) {
                    invalid[0][offsets.iter().position(|&i| i == 0).unwrap()] = bad;
                }
                assert!(
                    work.evaluate(&mut fields, &comm, [&invalid[0], &invalid[1]], None)
                        .is_err()
                );
            }
            let mut invalid = local.clone();
            if comm.rank() == 7 {
                invalid[1].push(0.);
            }
            assert!(
                work.evaluate(&mut fields, &comm, [&invalid[0], &invalid[1]], None)
                    .is_err()
            );
            let mut core = vec![0.; offsets.len()];
            if offsets.contains(&0) {
                core[offsets.iter().position(|&i| i == 0).unwrap()] = f64::NAN;
            }
            assert!(
                work.evaluate(&mut fields, &comm, [&local[0], &local[1]], Some(&core))
                    .is_err()
            );
            // A rejected evaluation must not contaminate the reusable workspace.
            let recovered = work
                .evaluate(&mut fields, &comm, [&local[0], &local[1]], None)
                .unwrap();
            assert!(recovered.potentials.iter().flatten().all(|v| v.is_finite()));
        });
    }
}
