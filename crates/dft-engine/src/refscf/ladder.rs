//! Reference cutoff ladder: the legacy C1 fraction, owned exact injection,
//! and a fresh fine driver. No legacy session or mixer state is transferred.
use super::{
    numerical,
    record::{Record, SetupTimings},
};
use crate::{model::ScfError, refstart::StartFields};
use dft_basis::{FftGridSpec, HalfSphere, PlaneWaveBasis};
use dft_mpi::Collective;
use dft_numerics::RealMatrix;
use dft_physics::data::Complex64;
use dft_solvers::distributed_fft::position_owners;

pub const FRACTION: f64 = 0.5;
/// A fixed estimator ratio, independent of material, volume and band count.
pub const THRESHOLD_RATIO: f64 = 100.;
pub const MAX_COARSE_MAPS: usize = 30;

pub struct Report {
    pub kept_hx: super::kept_hx::Report,
    pub coarse_records: Vec<Record>,
    pub coarse_setup: SetupTimings,
    /// Setup + maps + export + teardown; no coarse verification or public result assembly.
    pub coarse_seconds: f64,
    pub coarse_threshold_ry: f64,
    pub coarse_stopped: bool,
}

/// Owned coarse state for exact injection before fine-stage reorthogonalization.
#[derive(Clone)]
pub struct Seed {
    pub basis: PlaneWaveBasis,
    pub grid: FftGridSpec,
    pub half: HalfSphere,
    pub orbitals: RealMatrix,
    pub rows: Vec<usize>,
    pub density_g: Vec<Complex64>,
    pub indices: Vec<usize>,
}
impl Seed {
    pub fn inject<K: Collective>(
        self,
        basis: &PlaneWaveBasis,
        grid: &FftGridSpec,
        half: &HalfSphere,
        rows: &[usize],
        fields: &StartFields,
        comm: &K,
    ) -> Result<(Vec<Complex64>, RealMatrix), ScfError> {
        let complex_map = crate::refscf::prolong::prolong_indices(&self.basis, basis)?;
        let mut fine_pair = vec![usize::MAX; basis.vectors().len()];
        for (p, &g) in half.plus().iter().enumerate() {
            fine_pair[g as usize] = p;
        }
        let mut mapping = vec![None; self.basis.vectors().len()];
        mapping[0] = Some(0);
        for p in 1..self.half.plus().len() {
            let f = fine_pair[complex_map[self.half.plus()[p] as usize]];
            if f == usize::MAX
                || complex_map[self.half.minus()[p] as usize] != half.minus()[f] as usize
            {
                return Err(ScfError::InvalidInput(
                    "reference ladder: Gamma pair orientation",
                ));
            }
            mapping[2 * p - 1] = Some(2 * f - 1);
            mapping[2 * p] = Some(2 * f);
        }
        let owners =
            position_owners(grid.dims(), grid.g_to_fft(), comm.size()).map_err(numerical)?;
        let target: Vec<_> = (0..basis.vectors().len())
            .map(|r| owners[half.plus()[r.div_ceil(2)] as usize])
            .collect();
        let orbitals = inject_rows(comm, &mapping, &self.rows, &self.orbitals, &target, rows)?;

        let cd = self.grid.dims();
        let fd = grid.dims();
        let mut mapping = Vec::with_capacity(cd.iter().product());
        for i in 0..cd[0] {
            for j in 0..cd[1] {
                for k in 0..cd[2] {
                    let signed = crate::refscf::prolong::prolong_frequency;
                    let next = match (signed(i, cd[0]), signed(j, cd[1]), signed(k, cd[2])) {
                        (Some(a), Some(b), Some(c)) => {
                            let [x, y, z] =
                                FftGridSpec::map_miller([a, b, c], fd).map_err(numerical)?;
                            Some((x * fd[1] + y) * fd[2] + z)
                        }
                        _ => None,
                    };
                    mapping.push(next);
                }
            }
        }
        let indices: Vec<_> = (0..fd[0])
            .flat_map(|i| (0..fd[1]).flat_map(move |j| (0..fd[2]).map(move |k| [i, j, k])))
            .collect();
        let target = position_owners(fd, &indices, comm.size()).map_err(numerical)?;
        let mut values = RealMatrix::zeros(self.density_g.len(), 2).map_err(numerical)?;
        for (i, z) in self.density_g.iter().enumerate() {
            values.as_mut_slice()[i] = z.re;
            values.as_mut_slice()[self.density_g.len() + i] = z.im;
        }
        let out = inject_rows(
            comm,
            &mapping,
            &self.indices,
            &values,
            &target,
            fields.reciprocal_indices(),
        )?;
        let density = (0..out.rows())
            .map(|i| Complex64::new(out.as_slice()[i], out.as_slice()[out.rows() + i]))
            .collect();
        Ok((density, orbitals))
    }
}

/// Copy each owned row to its new owner, with no floating-point arithmetic.
/// Only destination rows travel; equal-block transport pads to the largest edge.
fn inject_rows<K: Collective>(
    comm: &K,
    mapping: &[Option<usize>],
    rows: &[usize],
    values: &RealMatrix,
    target: &[usize],
    fine_rows: &[usize],
) -> Result<RealMatrix, ScfError> {
    let p = comm.size();
    let mut edges = vec![Vec::new(); p];
    for (l, &r) in rows.iter().enumerate() {
        if let Some(f) = mapping[r] {
            edges[target[f]].push((l, f));
        }
    }
    let mut maxima = vec![0.; p];
    maxima[comm.rank()] = edges.iter().map(Vec::len).max().unwrap_or(0) as f64;
    let mut all = maxima.clone();
    comm.tree_sum_f64(&maxima, &mut all).map_err(numerical)?;
    let edge = all.into_iter().fold(0., f64::max) as usize;
    let width = values.cols() + 1;
    let block = edge * width;
    let mut send = vec![0.; p * block];
    // Zero means empty; indices are transmitted as exact integers plus one.
    for (dest, entries) in edges.iter().enumerate() {
        for (i, &(l, f)) in entries.iter().enumerate() {
            let at = dest * block + i * width;
            send[at] = (f + 1) as f64;
            for b in 0..values.cols() {
                send[at + 1 + b] = values.as_slice()[b * values.rows() + l];
            }
        }
    }
    let mut recv = vec![0.; send.len()];
    comm.alltoall_f64(&send, block, &mut recv)
        .map_err(numerical)?;
    let mut slot = vec![usize::MAX; target.len()];
    for (i, &r) in fine_rows.iter().enumerate() {
        slot[r] = i;
    }
    let mut out = RealMatrix::zeros(fine_rows.len(), values.cols()).map_err(numerical)?;
    for packet in recv.chunks_exact(width) {
        if packet[0] == 0. {
            continue;
        }
        let l = slot[packet[0] as usize - 1];
        for b in 0..values.cols() {
            out.as_mut_slice()[b * fine_rows.len() + l] = packet[1 + b];
        }
    }
    Ok(out)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{
        SpeciesPotential,
        model::GammaReal,
        model::ScfOptions,
        refstart::{self, StartInput},
    };
    use dft_numerics::linalg::ComplexMatrix;
    use dft_numerics::switches::RunConfig;

    #[test]
    fn cutoff_ladder_owned_prolongation_exactness() {
        for ranks in [1, 3, 8] {
            dft_mpi::simulate(ranks, |comm| {
                let _scope = RunConfig::from_pairs(&[]).unwrap().enter();
                let pp = crate::scf::pool_fixture::potential_with_atomic_start(
                    2., 1.8, true, 161, "SLA PZ",
                );
                let cell = dft_geometry::Cell::from_user_rows_bohr([
                    [5., 0., 0.],
                    [0.2, 6., 0.],
                    [0., 0.1, 9.],
                ])
                .unwrap();
                let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
                let potentials = [SpeciesPotential {
                    label: "Xx",
                    potential: &pp,
                }];
                let options = ScfOptions {
                    cutoff_hartree: 2.,
                    bands: Some(3),
                    gamma_real: GammaReal::On,
                    ..Default::default()
                };
                let coarse_options = ScfOptions {
                    cutoff_hartree: 1.,
                    ..options.clone()
                };
                let input = |options| StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options,
                };
                let mut coarse = refstart::start(input(&coarse_options), comm).unwrap();
                let mut fine = refstart::start(input(&options), comm).unwrap();
                // Use an analytic band-limited positive density, independent of
                // the writer and of the atomic charge's high-G clipping.
                let cd = coarse.grid.dims();
                let density: Vec<_> = coarse
                    .fields
                    .plane_offsets()
                    .iter()
                    .map(|&g| {
                        let x = g / (cd[1] * cd[2]);
                        (2. + 0.1 * (std::f64::consts::TAU * x as f64 / cd[0] as f64).cos())
                            / cell.volume_bohr3()
                    })
                    .collect();
                let density_g = coarse.fields.forward(comm, &density).unwrap();
                let coarse_density = coarse.fields.gather(comm, &density).unwrap();
                let gather = |n: usize, rows: &[usize], x: &RealMatrix| {
                    let mut data = vec![0.; n * x.cols()];
                    for b in 0..x.cols() {
                        for (i, &r) in rows.iter().enumerate() {
                            data[b * n + r] = x.as_slice()[b * x.rows() + i];
                        }
                    }
                    let mut out = data.clone();
                    comm.tree_sum_f64(&data, &mut out).unwrap();
                    out
                };
                let c = gather(
                    coarse.basis.vectors().len(),
                    coarse.operator.local_rows(),
                    &coarse.orbitals,
                );
                let mut complex = ComplexMatrix::zeros(coarse.basis.vectors().len(), 3).unwrap();
                for b in 0..3 {
                    coarse
                        .half
                        .unpack(
                            &c[b * complex.rows()..(b + 1) * complex.rows()],
                            complex.column_mut(b),
                        )
                        .unwrap();
                }
                let expected =
                    crate::refscf::prolong::prolong_orbitals(&coarse.basis, &fine.basis, &complex)
                        .unwrap();
                let map =
                    crate::refscf::prolong::prolong_indices(&coarse.basis, &fine.basis).unwrap();
                let seed = Seed {
                    basis: coarse.basis,
                    grid: coarse.grid,
                    half: coarse.half,
                    rows: coarse.operator.local_rows().to_vec(),
                    orbitals: coarse.orbitals,
                    indices: coarse.fields.reciprocal_indices().to_vec(),
                    density_g,
                };
                let (g, x) = seed
                    .inject(
                        &fine.basis,
                        &fine.grid,
                        &fine.half,
                        fine.operator.local_rows(),
                        &fine.fields,
                        comm,
                    )
                    .unwrap();
                let real = gather(fine.basis.vectors().len(), fine.operator.local_rows(), &x);
                let mut actual = ComplexMatrix::zeros(fine.basis.vectors().len(), 3).unwrap();
                for b in 0..3 {
                    fine.half
                        .unpack(
                            &real[b * actual.rows()..(b + 1) * actual.rows()],
                            actual.column_mut(b),
                        )
                        .unwrap();
                    for &r in &map {
                        assert_eq!(
                            expected.column(b)[r].re.to_bits(),
                            actual.column(b)[r].re.to_bits()
                        );
                        assert_eq!(
                            expected.column(b)[r].im.to_bits(),
                            actual.column(b)[r].im.to_bits()
                        );
                    }
                    for (&p, &m) in fine.half.plus().iter().zip(fine.half.minus()).skip(1) {
                        let (a, z) = (actual.column(b)[p as usize], actual.column(b)[m as usize]);
                        assert_eq!(a.re.to_bits(), z.re.to_bits());
                        assert_eq!((-a.im).to_bits(), z.im.to_bits());
                    }
                    for i in 0..=b {
                        let dot: Complex64 = actual
                            .column(i)
                            .iter()
                            .zip(actual.column(b))
                            .map(|(a, b)| a.conj() * b)
                            .sum();
                        assert!((dot - Complex64::new(f64::from(i == b), 0.)).norm() < 1e-12);
                    }
                }
                let planes = fine.fields.inverse(comm, &g).unwrap();
                let density = fine.fields.gather(comm, &planes).unwrap();
                let oracle =
                    crate::refscf::prolong::prolong_density(&coarse_density, fine.grid.dims())
                        .unwrap();
                for (a, b) in density.iter().zip(&oracle) {
                    assert!((a - b).abs() < 1e-15);
                }
                assert!(
                    (density.sum() * cell.volume_bohr3() / density.len() as f64 - 2.).abs() < 1e-13
                );
                for (&i, z) in fine.fields.reciprocal_indices().iter().zip(&g) {
                    if i == 0 {
                        assert!((z.re * cell.volume_bohr3() - 2.).abs() < 1e-14);
                    }
                }
            });
        }
    }
}
