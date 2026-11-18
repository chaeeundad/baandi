//! M1b: row-block distributed primitives against the serial GEMM, on simulated
//! ranks (threads), 1 to 8 ranks.
use dft_mpi::{Collective, simulate};
use dft_numerics::{
    pool::WorkerPool,
    real::{RealGemmWorkspace, RealMatrix, RealOp},
};
use dft_solvers::distributed::{RowBlock, column_norms_squared, gram, project_out, rotate_local};

fn matrix(n: usize, m: usize, phase: usize) -> RealMatrix {
    RealMatrix::new(
        n,
        m,
        (0..n * m)
            .map(|i| ((i * 7 + phase) as f64).sin() * (1.0 + (i % 11) as f64 * 0.013))
            .collect(),
    )
    .unwrap()
}
fn workspace(rows: usize, cols: usize) -> RealGemmWorkspace {
    RealGemmWorkspace::with_pool(rows.max(1), cols, WorkerPool::new(1).unwrap()).unwrap()
}
fn serial_gram(a: &RealMatrix, b: &RealMatrix) -> RealMatrix {
    let mut ws = workspace(a.rows(), a.cols().max(b.cols()));
    let mut out = RealMatrix::zeros(a.cols(), b.cols()).unwrap();
    ws.multiply_into(a.as_ref(), RealOp::Transpose, b.as_ref(), out.as_mut())
        .unwrap();
    out
}
fn bits(m: &RealMatrix) -> Vec<u64> {
    m.as_slice().iter().map(|x| x.to_bits()).collect()
}

#[test]
fn gram_matches_serial_within_rounding_and_is_identical_on_every_rank() {
    let (n, ma, mb) = (1013usize, 24usize, 9usize);
    let (a, b) = (matrix(n, ma, 1), matrix(n, mb, 5));
    let serial = serial_gram(&a, &b);
    let scale = serial
        .as_slice()
        .iter()
        .fold(1.0_f64, |x, y| x.max(y.abs()));
    for ranks in [1usize, 2, 3, 4, 8] {
        let results = simulate(ranks, |comm| {
            let block = RowBlock::new(n, comm.rank(), ranks);
            let (al, bl) = (block.scatter(&a).unwrap(), block.scatter(&b).unwrap());
            let mut ws = workspace(al.rows(), ma.max(mb));
            let mut out = RealMatrix::zeros(ma, mb).unwrap();
            gram(comm, &mut ws, &al, &bl, &mut out).unwrap();
            out
        });
        for out in &results {
            assert_eq!(bits(out), bits(&results[0]), "ranks differ at P={ranks}");
            for (got, want) in out.as_slice().iter().zip(serial.as_slice()) {
                assert!(
                    (got - want).abs() <= 8.0 * n as f64 * f64::EPSILON * scale,
                    "P={ranks}"
                );
            }
        }
        if ranks == 1 {
            assert_eq!(
                bits(&results[0]),
                bits(&serial),
                "one rank is the serial GEMM"
            );
        }
        // Reproducible for a fixed rank count.
        let again = simulate(ranks, |comm| {
            let block = RowBlock::new(n, comm.rank(), ranks);
            let (al, bl) = (block.scatter(&a).unwrap(), block.scatter(&b).unwrap());
            let mut ws = workspace(al.rows(), ma.max(mb));
            let mut out = RealMatrix::zeros(ma, mb).unwrap();
            gram(comm, &mut ws, &al, &bl, &mut out).unwrap();
            out
        });
        assert_eq!(bits(&again[0]), bits(&results[0]));
    }
}

#[test]
fn projection_rotation_and_norms_match_serial() {
    let (n, k, w) = (777usize, 12usize, 5usize);
    let (q, x) = (matrix(n, k, 3), matrix(n, w, 9));
    // Serial reference: X - Q (Q^T X).
    let mut ws = workspace(n, k.max(w));
    let overlap = serial_gram(&q, &x);
    let mut qo = RealMatrix::zeros(n, w).unwrap();
    ws.multiply_into(q.as_ref(), RealOp::Identity, overlap.as_ref(), qo.as_mut())
        .unwrap();
    let want: Vec<f64> = x
        .as_slice()
        .iter()
        .zip(qo.as_slice())
        .map(|(a, b)| a - b)
        .collect();
    let c = matrix(k, 4, 2);
    let mut rot = RealMatrix::zeros(n, 4).unwrap();
    ws.multiply_into(q.as_ref(), RealOp::Identity, c.as_ref(), rot.as_mut())
        .unwrap();
    for ranks in [1usize, 2, 4, 8] {
        let results = simulate(ranks, |comm| {
            let block = RowBlock::new(n, comm.rank(), ranks);
            let ql = block.scatter(&q).unwrap();
            let mut xl = block.scatter(&x).unwrap();
            let mut ws = workspace(ql.rows(), k.max(w).max(4));
            let mut ov = RealMatrix::zeros(k, w).unwrap();
            let mut scratch = RealMatrix::zeros(ql.rows(), w).unwrap();
            project_out(comm, &mut ws, &ql, &mut xl, &mut ov, &mut scratch).unwrap();
            let mut rl = RealMatrix::zeros(ql.rows(), 4).unwrap();
            rotate_local(&mut ws, &ql, &c, &mut rl).unwrap();
            let norms = column_norms_squared(comm, &xl).unwrap();
            (
                block.gather(comm, &xl).unwrap(),
                block.gather(comm, &rl).unwrap(),
                norms,
            )
        });
        for (xfull, rfull, norms) in &results {
            for (got, w2) in xfull.as_slice().iter().zip(&want) {
                assert!((got - w2).abs() <= 1e-11, "projection P={ranks}");
            }
            for (got, w2) in rfull.as_slice().iter().zip(rot.as_slice()) {
                assert!((got - w2).abs() <= 1e-11, "rotation P={ranks}");
            }
            for (j, nn) in norms.iter().enumerate() {
                let serial: f64 = xfull.as_slice()[j * n..(j + 1) * n]
                    .iter()
                    .map(|v| v * v)
                    .sum();
                assert!(
                    (nn - serial).abs() <= 1e-10 * serial.max(1.0),
                    "norm P={ranks}"
                );
            }
            assert_eq!(
                norms.iter().map(|x| x.to_bits()).collect::<Vec<_>>(),
                results[0].2.iter().map(|x| x.to_bits()).collect::<Vec<_>>()
            );
        }
    }
}
