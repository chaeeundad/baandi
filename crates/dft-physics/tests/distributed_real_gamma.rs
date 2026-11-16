use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_mpi::{Collective, MpiError, SimRank, simulate};
use dft_numerics::{RealMatrix, pool::WorkerPool};
use dft_physics::{
    RealHamiltonian,
    distributed_real_gamma::DistributedRealHamiltonian,
    hamiltonian::{Hamiltonian, NonlocalOperator, Species},
};
use dft_units::Hartree;
use ndarray::Array3;
use std::cell::Cell as Counter;

fn pseudo(ls: &[usize]) -> dft_upf::Pseudopotential {
    let n = 81;
    let a = |f: &dyn Fn(f64) -> f64| {
        (0..n)
            .map(|i| format!("{:.17e}", f(i as f64 * 0.15)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let mut beta = String::new();
    for (j, &l) in ls.iter().enumerate() {
        let i = j + 1;
        beta.push_str(&format!("<PP_BETA.{i} index=\"{i}\" angular_momentum=\"{l}\" cutoff_radius=\"12\" cutoff_radius_index=\"{n}\">{}</PP_BETA.{i}>",a(&|r| r.powi(l as i32+1)*(-r*r/(1.+j as f64*0.1)).exp())));
    }
    let dij = (0..ls.len())
        .flat_map(|i| {
            (0..ls.len()).map(move |j| {
                if i == j {
                    1. + i as f64 * 0.1
                } else if ls[i] == ls[j] {
                    0.2
                } else {
                    0.
                }
            })
        })
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    dft_upf::parse_upf(&format!(r#"<UPF version="2.0.1"><PP_INFO>Synthetic R3 fixture</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ" z_valence="4" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="3" l_max_rho="0" l_local="-1" mesh_size="{n}" number_of_wfc="0" number_of_proj="{}"/>
<PP_MESH dx="0.15" mesh="{n}" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH><PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>{dij}</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,ls.len(),a(&|r|r),a(&|_|0.15),a(&|_|0.),a(&|_|0.))).unwrap()
}
fn fixture(rows: [[f64; 3]; 3], atoms: usize, cutoff: f64) -> (PlaneWaveBasis, Hamiltonian) {
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
    let basis = PlaneWaveBasis::new(&cell, Hartree(cutoff)).unwrap();
    let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
    let positions: Vec<_> = (0..atoms)
        .map(|i| {
            cell.fractional_to_cartesian_bohr([
                ((i * 7 + 1) % 23) as f64 / 23.,
                ((i * 11 + 2) % 29) as f64 / 29.,
                ((i * 13 + 3) % 31) as f64 / 31.,
            ])
        })
        .collect();
    let p = pseudo(&[0, 0, 1, 1, 2, 2, 3, 3]);
    let q = pseudo(&[0, 1]);
    let split = atoms.div_ceil(2);
    let nl = NonlocalOperator::new(
        &basis,
        cell.volume_bohr3(),
        &[
            Species {
                potential: &p,
                positions_bohr: &positions[..split],
            },
            Species {
                potential: &q,
                positions_bohr: &positions[split..],
            },
        ],
    )
    .unwrap();
    let v = Array3::from_shape_fn(grid.dims(), |(i, j, k)| {
        0.3 * ((i * 13 + j * 7 + k * 3) as f64).sin()
    });
    let h = Hamiltonian::new(&cell, &basis, v, nl).unwrap();
    (basis, h)
}

fn block(n: usize, m: usize, seed: usize) -> RealMatrix {
    // Deterministic pseudo-random columns (not normalized; norms stay comparable).
    let mut s = (seed as u64)
        .wrapping_mul(2862933555777941757)
        .wrapping_add(3037000493);
    let data = (0..n * m)
        .map(|_| {
            s = s
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            ((s >> 11) as f64 / (1u64 << 53) as f64) - 0.5
        })
        .collect();
    RealMatrix::new(n, m, data).unwrap()
}

struct Counting<'a> {
    inner: &'a SimRank,
    alltoall: Counter<usize>,
    tree: Counter<usize>,
    gather: Counter<usize>,
}
impl Collective for Counting<'_> {
    fn rank(&self) -> usize {
        self.inner.rank()
    }
    fn size(&self) -> usize {
        self.inner.size()
    }
    fn barrier(&self) -> Result<(), MpiError> {
        self.inner.barrier()
    }
    fn allgatherv_f64(&self, l: &[f64], c: &[usize], o: &mut [f64]) -> Result<(), MpiError> {
        self.gather.set(self.gather.get() + 1);
        self.inner.allgatherv_f64(l, c, o)
    }
    fn tree_sum_f64(&self, l: &[f64], o: &mut [f64]) -> Result<(), MpiError> {
        self.tree.set(self.tree.get() + 1);
        self.inner.tree_sum_f64(l, o)
    }
    fn bcast_f64(&self, d: &mut [f64], r: usize) -> Result<(), MpiError> {
        self.inner.bcast_f64(d, r)
    }
    fn alltoall_f64(&self, s: &[f64], b: usize, r: &mut [f64]) -> Result<(), MpiError> {
        self.alltoall.set(self.alltoall.get() + 1);
        self.inner.alltoall_f64(s, b, r)
    }
}

type Counts = Vec<(usize, usize)>;

/// Distributed H of the global block on `ranks` simulated ranks, gathered to ORC order,
/// with per-rank measured and predicted (alltoall, tree_sum) calls of one apply.
fn distributed(
    rows: [[f64; 3]; 3],
    atoms: usize,
    x: &RealMatrix,
    ranks: usize,
) -> (RealMatrix, Counts, Counts) {
    // The collective counts asserted by the callers are those of the greedy
    // stick partition (two mirror exchanges per apply on more than one rank);
    // the mirror partition, default since 0d87d6e, has its own test below.
    distributed_sticks(rows, atoms, x, ranks, Some(false))
}

/// `distributed` with the stick rule forced on every rank (`Some(true)` =
/// `BAANDI_M2_STICKS=mirror`); `None` keeps the switch.
fn distributed_sticks(
    rows: [[f64; 3]; 3],
    atoms: usize,
    x: &RealMatrix,
    ranks: usize,
    mirror: Option<bool>,
) -> (RealMatrix, Counts, Counts) {
    let parts = simulate(ranks, |comm| {
        dft_solvers::distributed_fft::set_mirror_sticks_override(mirror);
        // Hamiltonian is not Sync: every simulated rank builds its own snapshot.
        let (basis, h) = fixture(rows, atoms, 1.2);
        let pool = WorkerPool::new(1).unwrap();
        let mut d =
            DistributedRealHamiltonian::new(&h, &basis, 16, comm.rank(), ranks, pool).unwrap();
        let counting = Counting {
            inner: comm,
            alltoall: Counter::new(0),
            tree: Counter::new(0),
            gather: Counter::new(0),
        };
        let xl = d.extract_local(x).unwrap();
        let mut yl = RealMatrix::zeros(xl.rows(), xl.cols()).unwrap();
        d.apply_block(&counting, &xl, &mut yl).unwrap();
        let measured = (counting.alltoall.get(), counting.tree.get());
        let mut global = RealMatrix::zeros(x.rows(), x.cols()).unwrap();
        d.scatter_local(&yl, &mut global).unwrap();
        (
            d.local_rows().to_vec(),
            global,
            measured,
            d.collectives_per_apply(x.cols()),
        )
    });
    let mut y = RealMatrix::zeros(x.rows(), x.cols()).unwrap();
    let n = x.rows();
    for (rows, global, _, _) in &parts {
        for c in 0..x.cols() {
            for &r in rows {
                y.as_mut_slice()[c * n + r] = global.as_slice()[c * n + r];
            }
        }
    }
    let measured = parts.iter().map(|p| p.2).collect();
    let predicted = parts.iter().map(|p| p.3).collect();
    (y, measured, predicted)
}

#[test]
fn distributed_h_matches_serial_on_1_to_8_ranks() {
    // (cell rows, atoms): cubic, skewed, elongated slab-like. The fixture has two
    // species with 8 and 2 projectors per atom (multi-projector, mixed l).
    let cases = [
        ([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], 2),
        ([[6.4, 0.4, 0.2], [1., 7.3, 0.3], [0.5, 0.7, 8.4]], 4),
        ([[4.2, 0., 0.], [0., 5.2, 0.], [0., 0., 14.2]], 3),
    ];
    let mut worst_by_ranks = [0.0_f64; 8];
    let mut bit_one_rank = true;
    for (case, (rows, atoms)) in cases.into_iter().enumerate() {
        let (basis, h) = fixture(rows, atoms, 1.2);
        let n = basis.vectors().len();
        let mut serial =
            RealHamiltonian::from_complex(&h, &basis, 16, WorkerPool::new(1).unwrap()).unwrap();
        for m in [1, 2, 3, 16] {
            let x = block(n, m, 5 + m);
            let mut want = RealMatrix::zeros(n, m).unwrap();
            serial.apply_block(&x, &mut want).unwrap();
            let norm = want.as_slice().iter().map(|v| v * v).sum::<f64>().sqrt();
            for ranks in 1..=8 {
                let (got, measured, predicted) = distributed(rows, atoms, &x, ranks);
                let err = got
                    .as_slice()
                    .iter()
                    .zip(want.as_slice())
                    .map(|(a, b)| (a - b) * (a - b))
                    .sum::<f64>()
                    .sqrt()
                    / norm;
                worst_by_ranks[ranks - 1] = worst_by_ranks[ranks - 1].max(err);
                assert!(err < 1e-12, "case {case} m {m} ranks {ranks}: {err:e}");
                if ranks == 1 {
                    bit_one_rank &= got.as_slice() == want.as_slice();
                }
                // Fixed collective counts (M2d batched transposes): 2 FFT transposes
                // for the whole block plus 2 mirror exchanges (none on one rank)
                // alltoalls, one tree_sum per apply. Unbatched: 2*ceil(m/2) + 2.
                let a2a = 2 + if ranks > 1 { 2 } else { 0 };
                for (meas, pred) in measured.iter().zip(&predicted) {
                    assert_eq!(*pred, (a2a, 1));
                    assert_eq!(*meas, (a2a, 1), "case {case} m {m} ranks {ranks}");
                }
                if m == 3 || m == 16 {
                    let (again, _, _) = distributed(rows, atoms, &x, ranks);
                    assert_eq!(again.as_slice(), got.as_slice(), "bit reproducibility");
                }
            }
        }
    }
    println!("max relative error by rank count 1..8: {worst_by_ranks:?}");
    println!("one rank bit-identical to serial: {bit_one_rank}");
}

#[test]
fn mirror_sticks_match_serial_without_mirror_exchanges() {
    // BAANDI_M2_STICKS=mirror (MPI-SCALING-GAP): with +G and -G of every pair
    // representative on one rank the two mirror alltoalls vanish (2 alltoalls and
    // one tree_sum per apply at every rank count); H X agrees with the serial
    // operator to rounding and is bit-reproducible at a fixed rank count.
    let cases = [
        ([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], 2),
        ([[6.4, 0.4, 0.2], [1., 7.3, 0.3], [0.5, 0.7, 8.4]], 4),
        ([[4.2, 0., 0.], [0., 5.2, 0.], [0., 0., 14.2]], 3),
    ];
    for (case, (rows, atoms)) in cases.into_iter().enumerate() {
        let (basis, h) = fixture(rows, atoms, 1.2);
        let n = basis.vectors().len();
        let mut serial =
            RealHamiltonian::from_complex(&h, &basis, 16, WorkerPool::new(1).unwrap()).unwrap();
        for m in [1, 3, 16] {
            let x = block(n, m, 11 + m);
            let mut want = RealMatrix::zeros(n, m).unwrap();
            serial.apply_block(&x, &mut want).unwrap();
            let norm = want.as_slice().iter().map(|v| v * v).sum::<f64>().sqrt();
            for ranks in 1..=6 {
                let (got, measured, predicted) =
                    distributed_sticks(rows, atoms, &x, ranks, Some(true));
                let err = got
                    .as_slice()
                    .iter()
                    .zip(want.as_slice())
                    .map(|(a, b)| (a - b) * (a - b))
                    .sum::<f64>()
                    .sqrt()
                    / norm;
                assert!(err < 1e-12, "case {case} m {m} ranks {ranks}: {err:e}");
                for (meas, pred) in measured.iter().zip(&predicted) {
                    assert_eq!(*pred, (2, 1), "case {case} m {m} ranks {ranks}");
                    assert_eq!(*meas, (2, 1), "case {case} m {m} ranks {ranks}");
                }
                let (again, _, _) = distributed_sticks(rows, atoms, &x, ranks, Some(true));
                assert_eq!(again.as_slice(), got.as_slice(), "bit reproducibility");
            }
        }
    }
}

#[test]
fn row_partition_is_a_permutation_and_apply_does_not_reallocate() {
    let (basis, h) = fixture([[6.4, 0.4, 0.2], [1., 7.3, 0.3], [0.5, 0.7, 8.4]], 4, 1.2);
    let n = basis.vectors().len();
    for ranks in 1..=8 {
        let mut seen = vec![0; n];
        for r in 0..ranks {
            let d = DistributedRealHamiltonian::new(
                &h,
                &basis,
                4,
                r,
                ranks,
                WorkerPool::new(1).unwrap(),
            )
            .unwrap();
            for &row in d.local_rows() {
                seen[row] += 1;
            }
            assert!(d.local_rows().windows(2).all(|w| w[0] < w[1]));
            if ranks == 1 {
                assert!(d.local_rows().iter().copied().eq(0..n));
            }
        }
        assert!(seen.iter().all(|&s| s == 1), "ranks {ranks}");
    }
    simulate(2, |comm| {
        let (basis, h) = fixture([[6.4, 0.4, 0.2], [1., 7.3, 0.3], [0.5, 0.7, 8.4]], 4, 1.2);
        let mut d = DistributedRealHamiltonian::new(
            &h,
            &basis,
            4,
            comm.rank(),
            2,
            WorkerPool::new(1).unwrap(),
        )
        .unwrap();
        let x = d.extract_local(&block(n, 3, 1)).unwrap();
        let mut y = RealMatrix::zeros(x.rows(), 3).unwrap();
        let p = y.as_slice().as_ptr();
        for _ in 0..3 {
            d.apply_block(comm, &x, &mut y).unwrap();
        }
        assert_eq!(p, y.as_slice().as_ptr());
    });
}

// ---- M2c: Davidson on the stick partition with the distributed H ----
use dft_solvers::{BandTolerances, RealDavidson, RealLinearOperator, SolverError};

struct SerialOp {
    h: std::cell::RefCell<RealHamiltonian>,
    n: usize,
    pool: std::sync::Arc<WorkerPool>,
}
impl RealLinearOperator for SerialOp {
    fn dimension(&self) -> usize {
        self.n
    }
    fn worker_pool(&self) -> Option<std::sync::Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) -> Result<(), SolverError> {
        let mut xi = RealMatrix::zeros(self.n, 1).unwrap();
        xi.as_mut_slice().copy_from_slice(x);
        let mut yo = RealMatrix::zeros(self.n, 1).unwrap();
        self.apply_block(&xi, &mut yo)?;
        y.copy_from_slice(yo.as_slice());
        Ok(())
    }
    fn apply_block(&self, x: &RealMatrix, y: &mut RealMatrix) -> Result<(), SolverError> {
        self.h
            .borrow_mut()
            .apply_block(x, y)
            .map_err(|e| SolverError::Application(e.to_string()))
    }
}

fn orc_diagonal(h: &Hamiltonian, plus: &[u32]) -> Vec<f64> {
    let n = 2 * plus.len() - 1;
    let mut d = vec![0.; n];
    d[0] = h.kinetic_diagonal()[0];
    for (i, &p) in plus.iter().enumerate().skip(1) {
        d[2 * i - 1] = h.kinetic_diagonal()[p as usize];
        d[2 * i] = d[2 * i - 1];
    }
    d
}

#[test]
fn davidson_on_stick_partition_matches_serial_and_reports_collectives() {
    let cases = [
        ([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], 2),
        ([[6.4, 0.4, 0.2], [1., 7.3, 0.3], [0.5, 0.7, 8.4]], 3),
    ];
    let bands = 6;
    let solver = RealDavidson {
        bands,
        tolerance: 1e-10,
        max_iterations: 300,
        max_subspace: 6 * bands,
        seed: 42,
    };
    let pol = BandTolerances {
        strict_bands: bands,
        buffer_tolerance: 1e-10,
    };
    for (rows, atoms) in cases {
        let (basis, h) = fixture(rows, atoms, 4.0);
        let pool = WorkerPool::new(1).unwrap();
        let real = RealHamiltonian::from_complex(&h, &basis, 6 * bands, pool.clone()).unwrap();
        let n = real.half_sphere().full_len();
        let diag = orc_diagonal(&h, real.half_sphere().plus());
        let op = SerialOp {
            h: std::cell::RefCell::new(real),
            n,
            pool,
        };
        let serial = solver
            .solve_with_band_tolerances(&op, None, &diag, pol)
            .unwrap();
        for ranks in 2..=5usize {
            let run = || {
                simulate(ranks, |comm| {
                    let (basis, h) = fixture(rows, atoms, 4.0);
                    let mut d = DistributedRealHamiltonian::new(
                        &h,
                        &basis,
                        6 * bands,
                        comm.rank(),
                        ranks,
                        WorkerPool::new(1).unwrap(),
                    )
                    .unwrap();
                    let dl: Vec<f64> = d.local_rows().iter().map(|&r| diag[r]).collect();
                    let counting = Counting {
                        inner: comm,
                        alltoall: Counter::new(0),
                        tree: Counter::new(0),
                        gather: Counter::new(0),
                    };
                    let r = solver
                        .solve_distributed_rows(&counting, &mut d, None, &dl, pol, false, true, 6)
                        .unwrap();
                    let rows_counts = (
                        counting.gather.get(),
                        counting.tree.get(),
                        counting.alltoall.get(),
                    );
                    // The M1 path (gather the block, apply H to owned pairs, gather products).
                    let counting1 = Counting {
                        inner: comm,
                        alltoall: Counter::new(0),
                        tree: Counter::new(0),
                        gather: Counter::new(0),
                    };
                    let real = RealHamiltonian::from_complex(
                        &h,
                        &basis,
                        6 * bands,
                        WorkerPool::new(1).unwrap(),
                    )
                    .unwrap();
                    let op1 = SerialOp {
                        h: std::cell::RefCell::new(real),
                        n,
                        pool: WorkerPool::new(1).unwrap(),
                    };
                    let r1 = solver
                        .solve_distributed(&counting1, &op1, None, &diag, pol, false, true, 6)
                        .unwrap();
                    let m1_counts = (
                        counting1.gather.get(),
                        counting1.tree.get(),
                        counting1.alltoall.get(),
                    );
                    (r, rows_counts, r1.iterations, m1_counts)
                })
            };
            let results = run();
            let bits =
                |m: &RealMatrix| m.as_slice().iter().map(|x| x.to_bits()).collect::<Vec<_>>();
            for (r, rc, it1, m1c) in &results {
                for (a, b) in r.eigenvalues.iter().zip(&serial.eigenvalues) {
                    assert!((a - b).abs() < 1e-9, "P={ranks}: {a} vs {b}");
                }
                for j in 0..bands {
                    let col = &r.eigenvectors.as_slice()[j * n..(j + 1) * n];
                    let mut y = vec![0.; n];
                    op.apply(col, &mut y).unwrap();
                    let res = y
                        .iter()
                        .zip(col)
                        .map(|(y, c)| (y - r.eigenvalues[j] * c).powi(2))
                        .sum::<f64>()
                        .sqrt();
                    assert!(res < 1e-8, "P={ranks} band {j} residual {res}");
                }
                assert_eq!(bits(&r.eigenvectors), bits(&results[0].0.eigenvectors));
                // Gathers: setup (counts, row lists) plus the final eigenvector gather,
                // independent of the iteration count.
                assert!(rc.0 <= 4, "allgatherv calls {}", rc.0);
                let it = r.iterations as f64;
                if std::ptr::eq(r, &results[0].0) {
                    println!(
                        "M2c n={n} P={ranks}: it={} rows (allgatherv,tree_sum,alltoall)/iteration = ({:.2},{:.2},{:.2}) total {:?}; M1 it={} per iteration = ({:.2},{:.2},{:.2}) total {:?}",
                        r.iterations,
                        rc.0 as f64 / it,
                        rc.1 as f64 / it,
                        rc.2 as f64 / it,
                        rc,
                        it1,
                        m1c.0 as f64 / *it1 as f64,
                        m1c.1 as f64 / *it1 as f64,
                        m1c.2 as f64 / *it1 as f64,
                        m1c
                    );
                }
            }
            let again = run();
            assert_eq!(
                bits(&results[0].0.eigenvectors),
                bits(&again[0].0.eigenvectors),
                "not reproducible at P={ranks}"
            );
        }
    }
}

#[test]
fn distributed_update_from_complex_follows_the_potential() {
    let (basis, h) = fixture([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], 2, 1.2);
    let mut d =
        DistributedRealHamiltonian::new(&h, &basis, 4, 0, 1, WorkerPool::new(1).unwrap()).unwrap();
    assert!(d.update_from_complex(&h, &basis));
    let (_, other) = fixture([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], 2, 1.2);
    assert!(
        !d.update_from_complex(&other, &basis),
        "a different operator must be rebuilt"
    );
}

// ---- M2d: batched transposes, distributed density and nonlocal energies ----

/// H x on `ranks` simulated ranks with batched or per-pair transposes:
/// (global result, alltoall calls per rank).
fn apply_mode(
    rows: [[f64; 3]; 3],
    atoms: usize,
    x: &RealMatrix,
    ranks: usize,
    batched: bool,
) -> (Vec<f64>, Vec<usize>) {
    apply_mode_with(rows, atoms, x, ranks, batched, false)
}

/// `apply_mode` with the nonlocal dgemm route (`BAANDI_M2_DGEMM`) forced.
fn apply_mode_with(
    rows: [[f64; 3]; 3],
    atoms: usize,
    x: &RealMatrix,
    ranks: usize,
    batched: bool,
    dgemm: bool,
) -> (Vec<f64>, Vec<usize>) {
    let parts = simulate(ranks, |comm| {
        // Greedy stick partition: the callers count its mirror exchanges.
        dft_solvers::distributed_fft::set_mirror_sticks_override(Some(false));
        let (basis, h) = fixture(rows, atoms, 1.2);
        let mut d = DistributedRealHamiltonian::new(
            &h,
            &basis,
            16,
            comm.rank(),
            ranks,
            WorkerPool::new(1).unwrap(),
        )
        .unwrap();
        d.set_batched(batched);
        d.set_dgemm(dgemm);
        let counting = Counting {
            inner: comm,
            alltoall: Counter::new(0),
            tree: Counter::new(0),
            gather: Counter::new(0),
        };
        let xl = d.extract_local(x).unwrap();
        let mut yl = RealMatrix::zeros(xl.rows(), xl.cols()).unwrap();
        d.apply_block(&counting, &xl, &mut yl).unwrap();
        assert_eq!(counting.alltoall.get(), d.collectives_per_apply(x.cols()).0);
        let mut global = RealMatrix::zeros(x.rows(), x.cols()).unwrap();
        d.scatter_local(&yl, &mut global).unwrap();
        (d.local_rows().to_vec(), global, counting.alltoall.get())
    });
    let n = x.rows();
    let mut y = vec![0.; n * x.cols()];
    for (rows, global, _) in &parts {
        for c in 0..x.cols() {
            for &r in rows {
                y[c * n + r] = global.as_slice()[c * n + r];
            }
        }
    }
    (y, parts.iter().map(|p| p.2).collect())
}

#[test]
fn batched_transposes_equal_per_pair_transposes_bitwise() {
    let cases = [
        ([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], 2),
        ([[6.4, 0.4, 0.2], [1., 7.3, 0.3], [0.5, 0.7, 8.4]], 4),
        ([[4.2, 0., 0.], [0., 5.2, 0.], [0., 0., 14.2]], 3),
    ];
    for (case, (rows, atoms)) in cases.into_iter().enumerate() {
        let (basis, _) = fixture(rows, atoms, 1.2);
        let n = basis.vectors().len();
        for m in [1, 3, 16] {
            let x = block(n, m, 31 + m);
            for ranks in [1, 2, 3, 5, 8] {
                let (yb, cb) = apply_mode(rows, atoms, &x, ranks, true);
                let (yu, cu) = apply_mode(rows, atoms, &x, ranks, false);
                let mirror = if ranks > 1 { 2 } else { 0 };
                assert!(cb.iter().all(|&c| c == 2 + mirror));
                assert!(cu.iter().all(|&c| c == 2 * m.div_ceil(2) + mirror));
                let bits = yb.iter().zip(&yu).all(|(a, b)| a.to_bits() == b.to_bits());
                let den = yu.iter().map(|v| v * v).sum::<f64>().sqrt();
                let err = yb
                    .iter()
                    .zip(&yu)
                    .map(|(a, b)| (a - b) * (a - b))
                    .sum::<f64>()
                    .sqrt()
                    / den;
                assert!(err <= 1e-14, "case {case} m {m} ranks {ranks}: {err:e}");
                println!(
                    "case {case} m {m} P={ranks}: alltoall per apply batched {} unbatched {}, bit-identical {bits}",
                    cb[0], cu[0]
                );
                assert!(bits, "batching must not change any bit");
            }
        }
    }
}

/// `BAANDI_M2_DGEMM=on` (nonlocal products through `multiply_slices_wide_into`):
/// bit-identical on pure builds (the same product), within 1e-14 relative on MKL
/// builds (per-panel dgemm), for 1..5 ranks, cubic, skewed and elongated cells.
#[test]
fn nonlocal_dgemm_route_matches_the_default_route() {
    let cases = [
        ([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], 2),
        ([[6.4, 0.4, 0.2], [1., 7.3, 0.3], [0.5, 0.7, 8.4]], 4),
        ([[4.2, 0., 0.], [0., 5.2, 0.], [0., 0., 14.2]], 3),
    ];
    for (case, (rows, atoms)) in cases.into_iter().enumerate() {
        let (basis, _) = fixture(rows, atoms, 1.2);
        let n = basis.vectors().len();
        for m in [1, 3, 16] {
            let x = block(n, m, 57 + m);
            for ranks in [1, 2, 3, 5] {
                let (y0, _) = apply_mode_with(rows, atoms, &x, ranks, true, false);
                let (y1, _) = apply_mode_with(rows, atoms, &x, ranks, true, true);
                let den = y0.iter().map(|v| v * v).sum::<f64>().sqrt();
                let err = y1
                    .iter()
                    .zip(&y0)
                    .map(|(a, b)| (a - b) * (a - b))
                    .sum::<f64>()
                    .sqrt()
                    / den;
                assert!(err <= 1e-14, "case {case} m {m} ranks {ranks}: {err:e}");
                if !dft_numerics::MKL_BACKEND {
                    assert!(
                        y0.iter().zip(&y1).all(|(a, b)| a.to_bits() == b.to_bits()),
                        "pure build: the dgemm route is the same product"
                    );
                }
            }
        }
    }
}

#[test]
fn distributed_band_densities_and_nonlocal_energies_match_serial() {
    use dft_basis::HalfSphere;
    use dft_numerics::{Fft3d, FftPlan3d};
    use num_complex::Complex64 as C;
    let cases = [
        ([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], 2),
        ([[6.4, 0.4, 0.2], [1., 7.3, 0.3], [0.5, 0.7, 8.4]], 4),
        ([[4.2, 0., 0.], [0., 5.2, 0.], [0., 0., 14.2]], 3),
    ];
    for (case, (rows, atoms)) in cases.into_iter().enumerate() {
        let (basis, h) = fixture(rows, atoms, 1.2);
        let grid = h.fft_grid().clone();
        let half = HalfSphere::new(&basis, &grid).unwrap();
        let n = half.full_len();
        let volume = Cell::from_user_rows_bohr(rows).unwrap().volume_bohr3();
        let m = 5;
        let x = block(n, m, 77);
        let weights = [2., 2., 1.5, 0.5, 2.];
        // Serial reference: complex coefficients, full-grid FFT, |psi|^2.
        let mut fft = FftPlan3d::new(grid.dims()).unwrap();
        let mut want = vec![0.; grid.dims().iter().product()];
        let mut want_nl = vec![0.; m];
        for b in 0..m {
            let mut c = vec![C::default(); n];
            half.unpack(&x.as_slice()[b * n..(b + 1) * n], &mut c)
                .unwrap();
            let full = grid.scatter(&c).unwrap();
            let mut field = Array3::<C>::zeros(grid.dims());
            fft.to_real(&full, &mut field).unwrap();
            for (w, v) in want.iter_mut().zip(field.iter()) {
                *w += weights[b] * v.norm_sqr() / volume;
            }
            h.nonlocal()
                .expectation_values(&[&c], &mut want_nl[b..b + 1])
                .unwrap();
        }
        let scale = want.iter().fold(0f64, |a, v| a.max(v.abs()));
        let mut previous: Option<Vec<f64>> = None;
        for ranks in 1..=6 {
            let run = || {
                simulate(ranks, |comm| {
                    let (basis, h) = fixture(rows, atoms, 1.2);
                    let mut d = DistributedRealHamiltonian::new(
                        &h,
                        &basis,
                        8,
                        comm.rank(),
                        ranks,
                        WorkerPool::new(1).unwrap(),
                    )
                    .unwrap();
                    let xl = d.extract_local(&x).unwrap();
                    let points = d.plane_points();
                    let mut contrib = vec![0.; m * points];
                    d.band_densities(comm, &xl, &weights, volume, &mut contrib)
                        .unwrap();
                    let mut local = vec![0.; points];
                    for b in 0..m {
                        for (l, c) in local.iter_mut().zip(&contrib[b * points..]) {
                            *l += c;
                        }
                    }
                    let mut full = vec![0.; grid.dims().iter().product()];
                    d.gather_planes(comm, &local, &mut full).unwrap();
                    let offsets = d.plane_grid_offsets();
                    for (p, &o) in offsets.iter().enumerate() {
                        assert_eq!(full[o].to_bits(), local[p].to_bits());
                    }
                    let mut nl = vec![0.; m];
                    d.nonlocal_band_energies(comm, &xl, &mut nl).unwrap();
                    (full, nl)
                })
            };
            let parts = run();
            let again = run();
            for (r, (full, nl)) in parts.iter().enumerate() {
                assert_eq!(full, &parts[0].0, "replicated density on rank {r}");
                assert_eq!(nl, &parts[0].1, "replicated energies on rank {r}");
                assert_eq!(full, &again[r].0, "bit reproducibility");
                assert_eq!(nl, &again[r].1, "bit reproducibility");
            }
            let (full, nl) = &parts[0];
            let err = full
                .iter()
                .zip(&want)
                .fold(0f64, |a, (g, w)| a.max((g - w).abs()))
                / scale;
            let nl_err = nl
                .iter()
                .zip(&want_nl)
                .fold(0f64, |a, (g, w)| a.max((g - w).abs() / w.abs().max(1e-300)));
            println!("case {case} P={ranks}: density rel {err:e}, nonlocal rel {nl_err:e}");
            assert!(err < 1e-13, "case {case} P={ranks}: density {err:e}");
            assert!(nl_err < 1e-12, "case {case} P={ranks}: nonlocal {nl_err:e}");
            // The plane split is a pure permutation: same bits for any rank count.
            if let Some(p) = &previous {
                assert_eq!(p, full, "density independent of the rank count");
            }
            previous = Some(full.clone());
        }
    }
}

/// `band_values` (distributed certificate screen): plane-local psi_b(r) equal to
/// the real part of the complex full-grid transform of the Hermitian-completed
/// column within 1e-13 relative, weight * psi^2 / volume equal to
/// `band_densities` bit for bit, and the values independent of the rank count
/// (the plane split is a permutation). Three cell shapes, ranks 1 to 5.
#[test]
fn distributed_band_values_match_serial_and_band_densities() {
    use dft_basis::HalfSphere;
    use dft_numerics::{Fft3d, FftPlan3d};
    use num_complex::Complex64 as C;
    let cases = [
        ([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], 2),
        ([[6.4, 0.4, 0.2], [1., 7.3, 0.3], [0.5, 0.7, 8.4]], 4),
        ([[4.2, 0., 0.], [0., 5.2, 0.], [0., 0., 14.2]], 3),
    ];
    for (case, (rows, atoms)) in cases.into_iter().enumerate() {
        let (basis, h) = fixture(rows, atoms, 1.2);
        let grid = h.fft_grid().clone();
        let half = HalfSphere::new(&basis, &grid).unwrap();
        let n = half.full_len();
        let volume = Cell::from_user_rows_bohr(rows).unwrap().volume_bohr3();
        let m = 5;
        let x = block(n, m, 91);
        let weights = [2., 2., 1.5, 0.5, 2.];
        let points_all: usize = grid.dims().iter().product();
        let mut fft = FftPlan3d::new(grid.dims()).unwrap();
        let mut want = vec![0.; m * points_all];
        for b in 0..m {
            let mut c = vec![C::default(); n];
            half.unpack(&x.as_slice()[b * n..(b + 1) * n], &mut c)
                .unwrap();
            let mut field = Array3::<C>::zeros(grid.dims());
            fft.to_real(&grid.scatter(&c).unwrap(), &mut field).unwrap();
            for (w, v) in want[b * points_all..].iter_mut().zip(field.iter()) {
                *w = v.re;
            }
        }
        let scale = want.iter().fold(0f64, |a, v| a.max(v.abs()));
        let mut previous: Option<Vec<f64>> = None;
        for ranks in 1..=5 {
            let parts = simulate(ranks, |comm| {
                let (basis, h) = fixture(rows, atoms, 1.2);
                let mut d = DistributedRealHamiltonian::new(
                    &h,
                    &basis,
                    8,
                    comm.rank(),
                    ranks,
                    WorkerPool::new(1).unwrap(),
                )
                .unwrap();
                let xl = d.extract_local(&x).unwrap();
                let points = d.plane_points();
                let mut values = vec![0.; m * points];
                d.band_values(comm, &xl, &mut values).unwrap();
                let mut density = vec![0.; m * points];
                d.band_densities(comm, &xl, &weights, volume, &mut density)
                    .unwrap();
                for b in 0..m {
                    for p in 0..points {
                        let v = values[b * points + p];
                        assert_eq!(
                            (weights[b] * v * v / volume).to_bits(),
                            density[b * points + p].to_bits()
                        );
                    }
                }
                // BAANDI_M2_DENSITY=fused: the accumulated plane density equals
                // the ascending-band sum of the per-band densities bit for bit.
                let mut sum = vec![0.; points];
                for b in 0..m {
                    for (s, c) in sum.iter_mut().zip(&density[b * points..(b + 1) * points]) {
                        *s += c;
                    }
                }
                let mut fused = vec![f64::NAN; points];
                d.band_density_sum(comm, &xl, &weights, volume, &mut fused)
                    .unwrap();
                assert_eq!(
                    fused.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    sum.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                    "fused plane density"
                );
                // Scatter this rank's planes into the full grid (zeros elsewhere).
                let offsets = d.plane_grid_offsets();
                let mut full = vec![0.; m * points_all];
                for b in 0..m {
                    for (p, &o) in offsets.iter().enumerate() {
                        full[b * points_all + o] = values[b * points + p];
                    }
                }
                full
            });
            let mut full = vec![0.; m * points_all];
            for part in &parts {
                for (f, v) in full.iter_mut().zip(part) {
                    if *v != 0. {
                        *f = *v;
                    }
                }
            }
            let err = full
                .iter()
                .zip(&want)
                .fold(0f64, |a, (g, w)| a.max((g - w).abs()))
                / scale;
            println!("case {case} P={ranks}: band values rel {err:e}");
            assert!(err < 1e-13, "case {case} P={ranks}: {err:e}");
            if let Some(p) = &previous {
                assert_eq!(p, &full, "values independent of the rank count");
            }
            previous = Some(full);
        }
    }
}

/// MPI-KEPT-START: H_new X = H_old X + (V_new - V_old) X. The potential change
/// from kept slabs (`band_densities_keeping` of the leading pairs) has the bits
/// of the change formed by transforming every pair, the kept sum equals a plain
/// apply with the new potential to rounding (gate 1e-13 relative), and the
/// density bits do not depend on keeping. 1 to 4 ranks; even and odd column
/// counts; two cells.
#[test]
fn potential_change_from_kept_planes_matches_plain_apply() {
    let cases = [
        ([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]], 2),
        ([[6.4, 0.4, 0.2], [1., 7.3, 0.3], [0.5, 0.7, 8.4]], 4),
    ];
    for (rows, atoms) in cases {
        let (basis, h) = fixture(rows, atoms, 1.2);
        let half = dft_basis::HalfSphere::new(&basis, h.fft_grid()).unwrap();
        let n = half.full_len();
        let volume = Cell::from_user_rows_bohr(rows).unwrap().volume_bohr3();
        for m in [5usize, 6] {
            let x = block(n, m, 400 + m);
            let occupied = 4;
            let weights = vec![2.; occupied];
            for ranks in 1..=4 {
                let parts = simulate(ranks, |comm| {
                    let (basis, mut h) = fixture(rows, atoms, 1.2);
                    let mut d = DistributedRealHamiltonian::new(
                        &h,
                        &basis,
                        8,
                        comm.rank(),
                        ranks,
                        WorkerPool::new(1).unwrap(),
                    )
                    .unwrap();
                    let xl = d.extract_local(&x).unwrap();
                    let nl = xl.rows();
                    let mut old = RealMatrix::zeros(nl, m).unwrap();
                    d.apply_block(comm, &xl, &mut old).unwrap();
                    // Density of the occupied prefix, keeping both full pairs.
                    let mut xm = RealMatrix::zeros(nl, occupied).unwrap();
                    xm.as_mut_slice()
                        .copy_from_slice(&xl.as_slice()[..nl * occupied]);
                    let pp = d.plane_points();
                    let mut plain_density = vec![0.; occupied * pp];
                    d.band_densities(comm, &xm, &weights, volume, &mut plain_density)
                        .unwrap();
                    let mut kept_density = vec![0.; occupied * pp];
                    let mut planes = vec![num_complex::Complex64::default(); 2 * pp];
                    d.band_densities_keeping(
                        comm,
                        &xm,
                        &weights,
                        volume,
                        &mut kept_density,
                        &mut planes,
                    )
                    .unwrap();
                    assert_eq!(plain_density, kept_density, "keeping changed the density");
                    // New potential of the same operator.
                    let v_old = d.potential().to_vec();
                    let dims = h.fft_grid().dims();
                    let v_new = Array3::from_shape_fn(dims, |(i, j, k)| {
                        0.3 * ((i * 13 + j * 7 + k * 3) as f64).sin()
                            + 0.05 * ((i + 2 * j + 5 * k) as f64).cos()
                    });
                    h.set_local_potential(v_new).unwrap();
                    assert!(d.update_from_complex(&h, &basis));
                    let delta: Vec<f64> = d
                        .potential()
                        .iter()
                        .zip(&v_old)
                        .map(|(a, b)| a - b)
                        .collect();
                    let mut change = Vec::new();
                    for kept in 0..=2 {
                        let mut y = RealMatrix::zeros(nl, m).unwrap();
                        d.apply_local_change(comm, &xl, &planes, kept, &delta, &mut y)
                            .unwrap();
                        change.push(y);
                    }
                    let mut fresh = RealMatrix::zeros(nl, m).unwrap();
                    d.apply_block(comm, &xl, &mut fresh).unwrap();
                    let mut sum = old.clone();
                    for (s, c) in sum.as_mut_slice().iter_mut().zip(change[2].as_slice()) {
                        *s += c;
                    }
                    let identical = change.iter().all(|c| c.as_slice() == change[0].as_slice());
                    let mut g_sum = RealMatrix::zeros(n, m).unwrap();
                    let mut g_fresh = RealMatrix::zeros(n, m).unwrap();
                    d.scatter_local(&sum, &mut g_sum).unwrap();
                    d.scatter_local(&fresh, &mut g_fresh).unwrap();
                    (d.local_rows().to_vec(), g_sum, g_fresh, identical)
                });
                let (mut sum, mut fresh) = (vec![0.; n * m], vec![0.; n * m]);
                for (rows_l, gs, gf, identical) in &parts {
                    assert!(identical, "kept planes changed the bits (P={ranks})");
                    for c in 0..m {
                        for &r in rows_l {
                            sum[c * n + r] = gs.as_slice()[c * n + r];
                            fresh[c * n + r] = gf.as_slice()[c * n + r];
                        }
                    }
                }
                let scale = fresh.iter().fold(0f64, |a, v| a.max(v.abs()));
                let err = sum
                    .iter()
                    .zip(&fresh)
                    .fold(0f64, |a, (s, f)| a.max((s - f).abs()))
                    / scale;
                println!("atoms {atoms} m {m} P={ranks}: H_old X + dV X vs H_new X rel {err:e}");
                assert!(err < 1e-13, "atoms {atoms} m {m} P={ranks}: {err:e}");
            }
        }
    }
}
