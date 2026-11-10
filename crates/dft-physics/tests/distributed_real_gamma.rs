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
    let parts = simulate(ranks, |comm| {
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
                // Fixed collective counts: 2*ceil(m/2) FFT transposes plus 2 mirror
                // exchanges (none on one rank) alltoalls, one tree_sum per apply.
                let pairs = m.div_ceil(2);
                let a2a = 2 * pairs + if ranks > 1 { 2 } else { 0 };
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
