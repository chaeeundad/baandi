//! M0 smoke test: `mpirun -np N target/release/examples/mpi_smoke`
//! (build with `--features mpi`; runs as a single rank without it).
use dft_mpi::{Universe, block_counts};

fn main() {
    let universe = Universe::init().expect("MPI init");
    let world = universe.world();
    let (rank, size) = (world.rank(), world.size());
    // Allgatherv of uneven blocks: rank r contributes block_counts(total)[r] items.
    let total = 11877usize;
    let counts = block_counts(total, size);
    let local: Vec<f64> = (0..counts[rank])
        .map(|i| (rank * 1_000_000 + i) as f64)
        .collect();
    let mut all = vec![0.0; total];
    world
        .allgatherv_f64(&local, &counts, &mut all)
        .expect("allgatherv");
    let mut at = 0;
    for (r, &c) in counts.iter().enumerate() {
        for i in 0..c {
            assert_eq!(
                all[at + i],
                (r * 1_000_000 + i) as f64,
                "allgatherv content"
            );
        }
        at += c;
    }
    // Ordered sum equals the sequential sum in rank order, bitwise.
    let mine = [0.1 * (rank as f64 + 1.0), 1e-17 * (rank as f64 + 1.0)];
    let mut sum = [0.0; 2];
    world.ordered_sum_f64(&mine, &mut sum).expect("ordered_sum");
    let mut expect = [0.0; 2];
    for r in 0..size {
        expect[0] += 0.1 * (r as f64 + 1.0);
        expect[1] += 1e-17 * (r as f64 + 1.0);
    }
    assert_eq!(sum[0].to_bits(), expect[0].to_bits());
    assert_eq!(sum[1].to_bits(), expect[1].to_bits());
    // Tree sum: the same bits on every rank, equal to the locally simulated tree.
    let vals: Vec<f64> = (0..size)
        .map(|r| 0.1 * (r as f64 + 1.0) + 1e-9 * (r * r) as f64)
        .collect();
    let mut tree = [0.0];
    world
        .tree_sum_f64(&[vals[rank]], &mut tree)
        .expect("tree_sum");
    let mut acc = vals.clone();
    let mut step = 1;
    while step < size {
        for r in (0..size).step_by(2 * step) {
            if r + step < size {
                acc[r] += acc[r + step];
            }
        }
        step <<= 1;
    }
    assert_eq!(tree[0].to_bits(), acc[0].to_bits(), "tree order");
    let mine_bits = [f64::from_bits(tree[0].to_bits())];
    let mut all_bits = vec![0.0; size];
    world
        .allgatherv_f64(&mine_bits, &vec![1; size], &mut all_bits)
        .expect("allgather bits");
    for b in &all_bits {
        assert_eq!(
            b.to_bits(),
            tree[0].to_bits(),
            "tree result differs between ranks"
        );
    }
    // Broadcast.
    let mut b = [0.0; 4];
    if rank == 0 {
        b = [1.0, 2.0, 3.0, 4.0];
    }
    world.bcast_f64(&mut b, 0).expect("bcast");
    assert_eq!(b, [1.0, 2.0, 3.0, 4.0]);
    // Alltoall with 3 items per pair.
    let send: Vec<f64> = (0..3 * size).map(|i| (rank * 100 + i) as f64).collect();
    let mut recv = vec![0.0; 3 * size];
    world.alltoall_f64(&send, 3, &mut recv).expect("alltoall");
    for r in 0..size {
        for k in 0..3 {
            assert_eq!(recv[3 * r + k], (r * 100 + 3 * rank + k) as f64);
        }
    }
    world.barrier().expect("barrier");
    if rank == 0 {
        println!("mpi_smoke OK on {size} rank(s)");
    }
}
