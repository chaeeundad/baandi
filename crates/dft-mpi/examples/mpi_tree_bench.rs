//! dev08 first run of the recursive-doubling tree sum (`BAANDI_MPI_TREE=rd`,
//! docs/performance/requests/mpi-fast-screen-2026-11-16.md):
//! `mpirun -np P target/release/examples/mpi_tree_bench` (build with `--features mpi`).
//! 1. Bit check: for several payloads the binomial and the recursive-doubling
//!    patterns return the same bits on every rank (P a power of two).
//! 2. Latency: microseconds per call of both patterns for payloads of 1, 8, 64,
//!    408 and 27744 doubles (the in-panel dot/norm sizes and a 408 x 68 Gram),
//!    rank 0 prints one line per payload.
use dft_mpi::{TreeSum, Universe, set_tree_sum};

fn main() {
    let universe = Universe::init().expect("MPI init");
    let world = universe.world();
    let (rank, size) = (world.rank(), world.size());
    let payloads = [1usize, 8, 64, 408, 408 * 68];
    for &n in &payloads {
        let local: Vec<f64> = (0..n)
            .map(|i| {
                0.1 * (rank as f64 + 1.0) + 1e-9 * (i * (rank + 3)) as f64
                    - 1e15 * ((rank % 2) as f64)
            })
            .collect();
        let mut a = vec![0.0; n];
        let mut b = vec![0.0; n];
        set_tree_sum(TreeSum::Binomial);
        world.tree_sum_f64(&local, &mut a).expect("binomial");
        set_tree_sum(TreeSum::RecursiveDoubling);
        world.tree_sum_f64(&local, &mut b).expect("doubling");
        let same = a.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits());
        if size.is_power_of_two() {
            assert!(same, "rank {rank}: patterns differ at payload {n}");
        }
        // Every rank holds the same bits as rank 0.
        let mut root = b.clone();
        world.bcast_f64(&mut root, 0).expect("bcast");
        assert!(root.iter().zip(&b).all(|(x, y)| x.to_bits() == y.to_bits()));
    }
    let reps = |n: usize| if n > 1000 { 200 } else { 5000 };
    for &n in &payloads {
        let local = vec![1.0 + rank as f64; n];
        let mut out = vec![0.0; n];
        let mut line = format!("np={size} payload {n} doubles:");
        for (name, pattern) in [
            ("binomial", TreeSum::Binomial),
            ("rd", TreeSum::RecursiveDoubling),
        ] {
            set_tree_sum(pattern);
            for _ in 0..50 {
                world.tree_sum_f64(&local, &mut out).expect("warm");
            }
            world.barrier().expect("barrier");
            let t = world.wtime();
            for _ in 0..reps(n) {
                world.tree_sum_f64(&local, &mut out).expect("tree");
            }
            let us = (world.wtime() - t) / reps(n) as f64 * 1e6;
            let mut worst = [us];
            let mut all = vec![0.0; size];
            world
                .allgatherv_f64(&worst, &vec![1; size], &mut all)
                .expect("gather");
            worst[0] = all.iter().cloned().fold(0.0, f64::max);
            line.push_str(&format!(" {name} {:.2} us", worst[0]));
        }
        if rank == 0 {
            println!("{line}");
        }
    }
    set_tree_sum(TreeSum::Binomial);
    if rank == 0 {
        println!(
            "mpi_tree_bench OK on {size} rank(s) (bit check {})",
            if size.is_power_of_two() {
                "enforced"
            } else {
                "skipped: not a power of two"
            }
        );
    }
}
