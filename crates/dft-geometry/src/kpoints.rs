//! Fractional reciprocal sampling; no space-group symmetry is assumed.
use thiserror::Error;

/// A fractional reciprocal point and its normalized Brillouin-zone weight.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct KPoint {
    /// Dimensionless reciprocal coordinates; constructors wrap to [-1/2,1/2).
    pub fractional: [f64; 3],
    /// Positive normalized integration weight, dimensionless.
    pub weight: f64,
}

/// Invalid sampling data; weights are never silently normalized.
#[derive(Clone, Debug, Error, PartialEq)]
#[error(
    "k points require finite coordinates, positive weights summing to one, and a nonzero grid with shifts 0 or 0.5"
)]
pub struct KPointError;

/// Wrap into [-1/2,1/2), preserving coordinates already in that interval.
pub fn wrap_fractional(k: [f64; 3]) -> Result<[f64; 3], KPointError> {
    if k.iter().any(|x| !x.is_finite()) {
        return Err(KPointError);
    }
    Ok(k.map(|x| {
        if (-0.5..0.5).contains(&x) {
            x
        } else {
            x - (x + 0.5).floor()
        }
    }))
}

// Pairwise summation keeps validation insensitive to linear accumulation error
// on larger grids; it does not change or renormalize any supplied weight.
fn weight_sum(points: &[KPoint]) -> f64 {
    match points.len() {
        0 => 0.0,
        1 => points[0].weight,
        n => {
            let (left, right) = points.split_at(n / 2);
            weight_sum(left) + weight_sum(right)
        }
    }
}

/// Validate a nonempty quadrature (absolute pairwise weight-sum allowance 1e-14).
pub fn validate(points: &[KPoint]) -> Result<(), KPointError> {
    if points.is_empty()
        || points.iter().any(|p| {
            p.fractional.iter().any(|x| !x.is_finite()) || !p.weight.is_finite() || p.weight <= 0.0
        })
        || (weight_sum(points) - 1.0).abs() > 1e-14
    {
        return Err(KPointError);
    }
    Ok(())
}

/// Uniform grid k_i=(m_i+s_i)/n_i, m_i=0..n_i-1, wrapped to [-1/2,1/2).
/// s_i is 0 or 1/2 in grid-step units. In particular the unshifted even grid
/// includes Gamma and boundary points (not the original centered MP convention).
pub fn monkhorst_pack(grid: [usize; 3], shift: [f64; 3]) -> Result<Vec<KPoint>, KPointError> {
    if grid.contains(&0) || shift.iter().any(|s| *s != 0.0 && *s != 0.5) {
        return Err(KPointError);
    }
    let count = grid
        .iter()
        .try_fold(1usize, |n, m| n.checked_mul(*m))
        .ok_or(KPointError)?;
    let mut points = Vec::new();
    points.try_reserve_exact(count).map_err(|_| KPointError)?;
    for i in 0..grid[0] {
        for j in 0..grid[1] {
            for k in 0..grid[2] {
                let m = [i, j, k];
                points.push(KPoint {
                    fractional: wrap_fractional(std::array::from_fn(|a| {
                        (m[a] as f64 + shift[a]) / grid[a] as f64
                    }))?,
                    weight: 1.0 / count as f64,
                });
            }
        }
    }
    validate(&points)?;
    Ok(points)
}

/// Merge only equivalent k and -k modulo reciprocal integers, summing weights.
/// Keeps the first representative. Also accepts an already reduced quadrature.
/// Coordinate equivalence uses an absolute 1e-12 fractional allowance.
pub fn reduce_time_reversal(points: &[KPoint]) -> Result<Vec<KPoint>, KPointError> {
    validate(points)?;
    let mut result: Vec<KPoint> = Vec::new();
    for point in points {
        let k = wrap_fractional(point.fractional)?;
        let found = result.iter_mut().find(|p| {
            [1.0, -1.0].iter().any(|s| {
                (0..3).all(|a| {
                    let d = k[a] - s * p.fractional[a];
                    (d - d.round()).abs() <= 1e-12
                })
            })
        });
        if let Some(p) = found {
            p.weight += point.weight;
        } else {
            result.push(KPoint {
                fractional: k,
                weight: point.weight,
            });
        }
    }
    validate(&result)?;
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn independent_orbits() {
        for n in [2, 3, 4] {
            for shifted in [false, true] {
                let shift = if shifted { 0.5 } else { 0.0 };
                let full = monkhorst_pack([n; 3], [shift; 3]).unwrap();
                let reduced = reduce_time_reversal(&full).unwrap();
                // Independent integer labels: inversion m -> -m-2s modulo n.
                let mut orbits = std::collections::BTreeSet::new();
                for i in 0..n {
                    for j in 0..n {
                        for k in 0..n {
                            let m = [i, j, k];
                            let partner = m.map(|x| (2 * n - x - usize::from(shifted)) % n);
                            orbits.insert(m.min(partner));
                        }
                    }
                }
                assert_eq!(reduced.len(), orbits.len());
                if n == 2 {
                    assert_eq!(reduced.len(), if shifted { 4 } else { 8 });
                }
                if n == 4 && !shifted {
                    assert_eq!(reduced.len(), 36);
                }
                let error = (reduced.iter().map(|p| p.weight).sum::<f64>() - 1.0).abs();
                println!(
                    "grid={n}^3 shift={shift}: full={} reduced={} weight error={error:.16e}",
                    full.len(),
                    reduced.len()
                );
                assert!(error < 1e-14);
                assert!(
                    full.iter()
                        .all(|p| p.fractional.iter().all(|x| (-0.5..0.5).contains(x)))
                );
                assert_eq!(reduce_time_reversal(&reduced).unwrap(), reduced);
            }
        }
        let large = monkhorst_pack([12; 3], [0.0; 3]).unwrap();
        let sum = large
            .chunks(32)
            .map(|c| c.iter().map(|p| p.weight).sum::<f64>())
            .sum::<f64>();
        println!(
            "12^3 weight error (independent chunk sum)={:.16e}",
            (sum - 1.0).abs()
        );
        assert!((sum - 1.0).abs() < 1e-14);
        assert!(monkhorst_pack([0, 2, 2], [0.0; 3]).is_err());
        assert!(monkhorst_pack([2; 3], [0.25; 3]).is_err());
        assert!(wrap_fractional([f64::NAN; 3]).is_err());
        assert!(reduce_time_reversal(&[]).is_err());
    }
}
