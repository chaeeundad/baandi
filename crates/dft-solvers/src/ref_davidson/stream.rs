//! P3-STREAM: contiguous, ordered streaming kernels. No fast-math or reassociation.
use super::{
    Iteration, error, preconditioner_denominator,
    work::{self, Stage, span},
};
use crate::SolverError;
use dft_mpi::Collective;
use dft_numerics::RealMatrix;

// A column-level guard keeps z*z well below overflow. The ordinary physical
// range uses sqrt, including underflow of z*z (1 + z*z then rounds to 1).
const SQRT_LIMIT: f64 = 3.273390607896142e150; // 2^500
const CHUNK: usize = 32;

#[inline]
pub(crate) fn sqrt_denominator(h: f64, e: f64) -> f64 {
    let x = 2.0 * (h - e);
    let z = (x - 1.0).abs();
    // Same cancellation-free rationalization as the baseline; max is a
    // branchless select on finite inputs, not a change of preconditioner.
    0.25 * (2.0 * x.max(1.0) + 1.0 / ((1.0 + z * z).sqrt() + z))
}

pub(crate) fn sqrt_safe(bounds: (f64, f64), e: f64) -> bool {
    [bounds.0, bounds.1].into_iter().all(|h| {
        let z = (2.0 * (h - e) - 1.0).abs();
        z <= SQRT_LIMIT
    })
}

/// Each chunk computes independent quotients/squares (auto-vectorizable), then
/// accumulates those squares in exactly ascending row order. Stack scratch
/// breaks the loop-carried dependency without changing the reduction tree.
pub(super) fn column<const SQRT: bool>(
    h: &[f64],
    x: &[f64],
    hx: &[f64],
    out: &mut [f64],
    e: f64,
) -> f64 {
    assert_eq!(h.len(), out.len());
    assert_eq!(x.len(), out.len());
    assert_eq!(hx.len(), out.len());
    let mut norm = 0.0;
    let mut squares = [0.0; CHUNK];
    for (((h, x), hx), out) in h
        .chunks(CHUNK)
        .zip(x.chunks(CHUNK))
        .zip(hx.chunks(CHUNK))
        .zip(out.chunks_mut(CHUNK))
    {
        for i in 0..out.len() {
            let d = if SQRT {
                sqrt_denominator(h[i], e)
            } else {
                preconditioner_denominator(h[i], e)
            };
            let z = (hx[i] - e * x[i]) / d;
            out[i] = z;
            squares[i] = z * z;
        }
        for &square in &squares[..out.len()] {
            norm += square;
        }
    }
    norm
}

pub(super) fn count(profile: bool, rows: usize, sqrt_columns: usize, hypot_columns: usize) {
    work::add(profile, |s| {
        let sqrt = (rows * sqrt_columns) as u64;
        let hypot = (rows * hypot_columns) as u64;
        s.scalar_flops += 15 * sqrt + 13 * hypot;
        s.sqrt_calls += sqrt;
        s.hypot_calls += hypot;
        s.streamed_bytes += 32 * (sqrt + hypot) + 8 * (sqrt_columns + hypot_columns) as u64;
    });
}

#[allow(clippy::too_many_arguments)]
pub(super) fn normalize<C: Collective>(
    c: &C,
    e: &mut RealMatrix,
    v: &mut RealMatrix,
    norms: &[f64],
    sums: &mut [f64],
    m: usize,
    cap: usize,
    r: &mut Iteration,
    profile: bool,
) -> Result<Option<usize>, SolverError> {
    let _stage = span(profile, Stage::Normalize);
    let rows = e.rows();
    let k = e.cols();
    c.tree_sum_f64(&norms[..k], &mut sums[..k]).map_err(error)?;
    // Capacity-exhaustion repair must see the original basis. Only stage an
    // append when the whole uncompressed block fits (the ordinary path).
    let stage = m + k <= cap;
    if stage {
        v.reshape_within_storage(rows, m + k)?;
    }
    let mut kept = 0;
    for (j, &sum) in sums[..k].iter().enumerate() {
        let norm = sum.sqrt();
        if !norm.is_finite() {
            return Err(SolverError::Application(
                "nonfinite reference correction".into(),
            ));
        }
        if norm <= 32.0 * f64::EPSILON {
            r.dropped_corrections += 1;
            continue;
        }
        let data = e.as_mut_slice();
        if stage {
            let dst = &mut v.as_mut_slice()[(m + kept) * rows..(m + kept + 1) * rows];
            for (i, dst) in dst.iter_mut().enumerate() {
                let value = data[j * rows + i] / norm;
                data[kept * rows + i] = value;
                *dst = value;
            }
        } else {
            for i in 0..rows {
                data[kept * rows + i] = data[j * rows + i] / norm;
            }
        }
        kept += 1;
    }
    e.reshape_within_storage(rows, kept)?;
    if stage {
        v.reshape_within_storage(rows, m + kept)?;
    }
    work::add(profile, |s| {
        s.sqrt_calls += k as u64;
        s.scalar_flops += (rows * kept) as u64;
        s.streamed_bytes += (if stage { 24 } else { 16 }) * (rows * kept) as u64;
        // The second store duplicates the scaled value, with no reread of E.
        if stage {
            s.copied_bytes += 8 * (rows * kept) as u64;
        }
    });
    Ok(stage.then_some(m))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn sqrt_bound_and_ordered_chunk_norm() {
        let mut max_relative = 0.0_f64;
        let mut differing = 0;
        // Dense ordinary range and both signs over 1000 binary exponents.
        let values = (-10000..=10000)
            .map(|i| i as f64 / 37.0)
            .chain((-500..=499).flat_map(|i| [2.0_f64.powi(i), -2.0_f64.powi(i)]));
        let random = (0..65536).scan(0x1234_5678_9abc_def0_u64, |state, i| {
            *state = state.wrapping_mul(6364136223846793005).wrapping_add(1);
            let exponent = (1013 + i % 21) as u64;
            Some(
                f64::from_bits((*state & ((1 << 52) - 1)) | (exponent << 52))
                    * if i % 2 == 0 { 1.0 } else { -1.0 },
            )
        });
        for x in values.chain(random) {
            let h = x * 0.5;
            assert!(sqrt_safe((h, h), 0.0));
            let a = preconditioner_denominator(h, 0.0);
            let b = sqrt_denominator(h, 0.0);
            max_relative = max_relative.max((a - b).abs() / a);
            differing += usize::from(a.to_bits() != b.to_bits());
            assert!((a - b).abs() <= 8.0 * f64::EPSILON * a, "x={x} a={a} b={b}");
        }
        assert!(!sqrt_safe((-f64::MAX / 4.0, 1.0), 0.0));
        assert!(!sqrt_safe((0.0, 1.0), f64::INFINITY));
        for rows in [0, 1, 31, 32, 33, 65, 257] {
            let h: Vec<_> = (0..rows).map(|i| i as f64 * 0.17 - 5.0).collect();
            let x: Vec<_> = (0..rows).map(|i| (i as f64).sin()).collect();
            let hx: Vec<_> = (0..rows).map(|i| (i as f64).cos()).collect();
            let mut out = vec![0.0; rows];
            let norm = column::<true>(&h, &x, &hx, &mut out, 0.3);
            let mut expected = 0.0;
            for i in 0..rows {
                let z = (hx[i] - 0.3 * x[i]) / sqrt_denominator(h[i], 0.3);
                assert_eq!(z.to_bits(), out[i].to_bits());
                expected += z * z;
            }
            assert_eq!(norm.to_bits(), expected.to_bits());
            let norm = column::<false>(&h, &x, &hx, &mut out, 0.3);
            let mut expected = 0.0;
            for i in 0..rows {
                let z = (hx[i] - 0.3 * x[i]) / preconditioner_denominator(h[i], 0.3);
                assert_eq!(z.to_bits(), out[i].to_bits());
                expected += z * z;
            }
            assert_eq!(norm.to_bits(), expected.to_bits());
        }
        eprintln!(
            "STREAM denominator samples=87537 changed={differing} max_relative={max_relative:e} (dimensionless)"
        );
    }
    #[test]
    fn normalized_compaction_stages_only_retained_columns() {
        use crate::ref_davidson::{Record, ReducedPlacement};
        dft_mpi::simulate(1, |c| {
            for cap in [4, 5] {
                let mut e = RealMatrix::new(2, 3, vec![0., 0., 3., 4., 0., 2.]).unwrap();
                let mut v = RealMatrix::new(2, 5, vec![123.; 10]).unwrap();
                let mut sums = [0.; 3];
                let mut r = Record::new(1, 5, false, ReducedPlacement::Root);
                let staged = normalize(
                    c,
                    &mut e,
                    &mut v,
                    &[0., 25., 4.],
                    &mut sums,
                    2,
                    cap,
                    &mut r.setup,
                    false,
                )
                .unwrap();
                assert_eq!(e.as_slice(), &[0.6, 0.8, 0., 1.]);
                assert_eq!(r.setup.dropped_corrections, 1);
                assert_eq!(staged, (cap == 5).then_some(2));
                assert_eq!(&v.as_slice()[..4], &[123.; 4]);
                if cap == 5 {
                    assert_eq!(&v.as_slice()[4..], e.as_slice());
                }
            }
        });
    }
}
