//! dev08 microbenchmark of the M2 collectives at the payloads of one Davidson
//! solve (ticket MPI-SCALING-GAP, docs/performance/mpi-scaling-gap.md):
//! `mpirun -np P target/release/examples/mpi_transpose_bench` (build with
//! `--features mpi`; np 2, 4 and 8 with the binding of the timed runs).
//!
//! Question: does the per-call time of the FFT transposes (one `alltoall_f64` per
//! direction per apply) and of the N x B return gather (`allgatherv_f64`) shrink
//! with the rank count? Both move a node-level volume that does not depend on P
//! (transposes: every occupied stick of length n0 once per pair; gather: N x B
//! doubles to every rank, i.e. P times N x B in total), so on one node their time
//! per call is expected to stay flat or grow, unlike the compute stages.
//!
//! Payloads (C, from the stick partition of the stage-0 inputs, `sticks.py` in the
//! analysis): si32 499 sticks of length 72, N = 11877, B = 68; si64 1005 sticks of
//! length 72, N = 23871, B = 132. Transpose block per rank pair:
//! 2 ceil(S/P) ceil(n0/P) pairs doubles (the padded block of `DistributedPairFft`),
//! for 1, 6 and B/2 pairs per call (single column, a typical warm-iteration
//! append, the initial apply). Rank 0 prints one line per payload: worst-rank
//! microseconds per call and the node volume per call (useful doubles x 8 bytes).
use dft_mpi::Universe;

fn main() {
    let universe = Universe::init().expect("MPI init");
    let world = universe.world();
    let (rank, size) = (world.rank(), world.size());
    let worst = |us: f64| -> f64 {
        let mut all = vec![0.0; size];
        world
            .allgatherv_f64(&[us], &vec![1; size], &mut all)
            .expect("gather");
        all.iter().cloned().fold(0.0, f64::max)
    };
    for (case, sticks, n0, n, bands) in [
        ("si32", 499usize, 72usize, 11877usize, 68usize),
        ("si64", 1005, 72, 23871, 132),
    ] {
        for pairs in [1usize, 6, bands / 2] {
            let block = 2 * sticks.div_ceil(size) * n0.div_ceil(size) * pairs;
            let send: Vec<f64> = (0..block * size).map(|i| (i + rank) as f64).collect();
            let mut recv = vec![0.0; block * size];
            let reps = (200_000_000 / (block * size * size).max(1)).clamp(5, 2000);
            for _ in 0..5 {
                world.alltoall_f64(&send, block, &mut recv).expect("warm");
            }
            world.barrier().expect("barrier");
            let t = world.wtime();
            for _ in 0..reps {
                world
                    .alltoall_f64(&send, block, &mut recv)
                    .expect("alltoall");
            }
            let us = worst((world.wtime() - t) / reps as f64 * 1e6);
            if rank == 0 {
                let useful = 2.0 * (sticks * n0 * pairs) as f64 * 8.0;
                println!(
                    "np={size} {case} alltoall pairs {pairs}: block {block} doubles, {us:.1} us per call (worst rank), node volume {:.2} MB, {:.2} GB/s",
                    useful / 1e6,
                    useful / (us * 1e3)
                );
            }
        }
        // Return gather of the eigenvectors: N x B doubles assembled on every rank.
        let counts: Vec<usize> = (0..size)
            .map(|r| (n / size + usize::from(r < n % size)) * bands)
            .collect();
        let local = vec![1.0 + rank as f64; counts[rank]];
        let mut full = vec![0.0; n * bands];
        let reps = 20;
        for _ in 0..3 {
            world
                .allgatherv_f64(&local, &counts, &mut full)
                .expect("warm");
        }
        world.barrier().expect("barrier");
        let t = world.wtime();
        for _ in 0..reps {
            world
                .allgatherv_f64(&local, &counts, &mut full)
                .expect("allgatherv");
        }
        let us = worst((world.wtime() - t) / reps as f64 * 1e6);
        if rank == 0 {
            let total = (n * bands * size) as f64 * 8.0;
            println!(
                "np={size} {case} allgatherv N x B ({n} x {bands}): {us:.1} us per call (worst rank), node receive volume {:.1} MB, {:.2} GB/s",
                total / 1e6,
                total / (us * 1e3)
            );
        }
    }
    if rank == 0 {
        println!("mpi_transpose_bench OK on {size} rank(s)");
    }
}
