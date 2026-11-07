#![cfg(feature = "mkl")]
use dft_mkl::{Error, Plan3d};
use num_complex::Complex64 as C;
use std::{hint::black_box, time::Instant};

fn input(shape: [usize; 3], pruned: bool) -> Vec<C> {
    let mut state = 0x123456789abcdefu64;
    (0..shape.iter().product())
        .map(|i| {
            let mut next = || {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                (state >> 11) as f64 / (1u64 << 53) as f64 - 0.5
            };
            let z = C::new(next(), next());
            if pruned && i % 13 != 0 {
                C::default()
            } else {
                z
            }
        })
        .collect()
}

#[test]
fn in_place_matches_wrapper_bits() {
    for shape in [[3, 4, 5], [25; 3], [36; 3], [72, 36, 36], [7, 9, 11]] {
        let mut plan = Plan3d::new(shape).unwrap();
        for pruned in [false, true] {
            let source = input(shape, pruned);
            for inverse in [false, true] {
                let mut expected = source.clone();
                plan.transform(&source, &mut expected, inverse).unwrap();
                // Offset by one complex element to exercise caller alignment too.
                for offset in [0, 1] {
                    let mut storage = vec![C::default(); source.len() + offset];
                    let actual = &mut storage[offset..];
                    actual.copy_from_slice(&source);
                    plan.transform_in_place(actual, inverse).unwrap();
                    for (i, (a, b)) in actual.iter().zip(&expected).enumerate() {
                        assert_eq!(
                            (a.re.to_bits(), a.im.to_bits()),
                            (b.re.to_bits(), b.im.to_bits()),
                            "shape={shape:?} pruned={pruned} inverse={inverse} offset={offset} i={i}"
                        );
                    }
                }
            }
        }
        for len in [
            0,
            shape.iter().product::<usize>() - 1,
            shape.iter().product::<usize>() + 1,
        ] {
            assert_eq!(
                plan.transform_in_place(&mut vec![C::default(); len], true),
                Err(Error::Dimension)
            );
        }
    }
}

#[test]
#[ignore = "dev08 release timing; see docs/performance/p39b-results.md"]
fn in_place_microbenchmark() {
    for shape in [[25; 3], [36; 3], [72, 36, 36]] {
        let source = input(shape, false);
        let mut plan = Plan3d::new(shape).unwrap();
        let mut old = source.clone();
        let mut new = source.clone();
        for inverse in [true, false] {
            for _ in 0..8 {
                plan.transform(&source, &mut old, inverse).unwrap();
                new.copy_from_slice(&source);
                plan.transform_in_place(&mut new, inverse).unwrap();
            }
            let mut samples = [Vec::with_capacity(32), Vec::with_capacity(32)];
            for round in 0..32 {
                for path in if round % 2 == 0 { [0, 1] } else { [1, 0] } {
                    // Reset outside the timed interval: repeated inverse calls otherwise overflow.
                    new.copy_from_slice(&source);
                    let start = Instant::now();
                    if path == 0 {
                        plan.transform(black_box(&source), black_box(&mut old), inverse)
                            .unwrap();
                        black_box(&old);
                    } else {
                        plan.transform_in_place(black_box(&mut new), inverse)
                            .unwrap();
                        black_box(&new);
                    }
                    samples[path].push(start.elapsed().as_secs_f64() * 1e6);
                }
            }
            for (path, mut times) in samples.into_iter().enumerate() {
                println!("shape={shape:?} inverse={inverse} path={path} samples_us={times:?}");
                times.sort_by(f64::total_cmp);
                println!(
                    "shape={shape:?} inverse={inverse} path={} n=32 median_us={:.6} min_us={:.6} max_us={:.6}",
                    ["wrapper", "in_place"][path],
                    (times[15] + times[16]) / 2.0,
                    times[0],
                    times[31]
                );
            }
        }
    }
}
