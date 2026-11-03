//! Piecewise reciprocal-space paths in the supplied primitive basis.
use crate::Cell;
use thiserror::Error;

#[derive(Clone, Debug)]
pub struct Vertex {
    pub label: String,
    pub fractional: [f64; 3],
    /// Start a disconnected branch here; no distance is added across the break.
    pub break_before: bool,
}
#[derive(Clone, Debug)]
pub struct PathLabel {
    pub label: String,
    pub index: usize,
    pub distance: f64,
}
#[derive(Clone, Debug)]
pub struct KPath {
    pub fractional: Vec<[f64; 3]>,
    pub distances: Vec<f64>,
    pub labels: Vec<PathLabel>,
    pub breaks: Vec<usize>,
}
#[derive(Debug, Error)]
#[error("invalid k path: {0}")]
pub struct PathError(pub &'static str);

fn norm(cell: &Cell, delta: [f64; 3]) -> f64 {
    let b = cell.reciprocal_columns_inverse_bohr();
    (0..3)
        .map(|i| (0..3).map(|j| b[i][j] * delta[j]).sum::<f64>().powi(2))
        .sum::<f64>()
        .sqrt()
}
/// `intervals` subdivisions per connected segment, or a maximum target spacing
/// in inverse Bohr. Segment endpoints are included exactly once.
pub fn generate(
    cell: &Cell,
    vertices: &[Vertex],
    intervals: usize,
    spacing: Option<f64>,
) -> Result<KPath, PathError> {
    if vertices.len() < 2
        || vertices.len() > 1_000_000
        || vertices[0].break_before
        || intervals == 0
        || spacing.is_some_and(|s| !s.is_finite() || s <= 0.0)
        || vertices.iter().any(|v| {
            v.label.is_empty()
                || v.label.chars().any(char::is_control)
                || v.fractional.iter().any(|v| !v.is_finite())
        })
    {
        return Err(PathError(
            "finite labeled vertices and positive sampling required",
        ));
    }
    let mut result = KPath {
        fractional: vec![vertices[0].fractional],
        distances: vec![0.0],
        labels: vec![PathLabel {
            label: vertices[0].label.clone(),
            index: 0,
            distance: 0.0,
        }],
        breaks: Vec::new(),
    };
    let mut distance = 0.0;
    for pair in vertices.windows(2) {
        let a = &pair[0];
        let b = &pair[1];
        if b.break_before {
            result.breaks.push(result.fractional.len());
            result.fractional.push(b.fractional);
            result.distances.push(distance);
        } else {
            let delta = std::array::from_fn(|i| b.fractional[i] - a.fractional[i]);
            let length = norm(cell, delta);
            if !length.is_finite() || length <= 0.0 {
                return Err(PathError("zero or nonfinite segment length"));
            }
            let count = spacing.map_or(intervals, |s| (length / s).ceil().max(1.0) as usize);
            if count > 1_000_000 || result.fractional.len().saturating_add(count) > 1_000_000 {
                return Err(PathError("path exceeds one million points"));
            }
            for step in 1..=count {
                let t = step as f64 / count as f64;
                result.fractional.push(if step == count {
                    b.fractional
                } else {
                    std::array::from_fn(|i| a.fractional[i] + t * delta[i])
                });
                result.distances.push(distance + length * t);
            }
            distance += length;
        }
        result.labels.push(PathLabel {
            label: b.label.clone(),
            index: result.fractional.len() - 1,
            distance,
        });
    }
    if !distance.is_finite() || distance <= 0.0 {
        return Err(PathError(
            "path requires a finite positive connected length",
        ));
    }
    Ok(result)
}

/// Standard ordered primitive vectors: fcc (0,a/2,a/2), (a/2,0,a/2),
/// (a/2,a/2,0); bcc (-a/2,a/2,a/2), (a/2,-a/2,a/2),
/// (a/2,a/2,-a/2); sc Cartesian a axes. Gram matrices must agree
/// to relative 1e-10 and handedness must be positive (proper rotation only).
pub fn preset(cell: &Cell, name: &str) -> Result<Vec<Vertex>, PathError> {
    let cosine = match name {
        "fcc" => 0.5,
        "bcc" => -1.0 / 3.0,
        "sc" => 0.0,
        _ => return Err(PathError("preset must be fcc, bcc or sc")),
    };
    let a = cell.lattice_columns_bohr();
    let gram: [[f64; 3]; 3] =
        std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|r| a[r][i] * a[r][j]).sum()));
    let scale = gram[0][0];
    let determinant = a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0]);
    if determinant <= 0.0
        || (0..3).any(|i| {
            (0..3).any(|j| (gram[i][j] / scale - if i == j { 1.0 } else { cosine }).abs() > 1e-10)
        })
    {
        return Err(PathError(
            "cell does not match the ordered standard primitive preset up to a rigid rotation; use explicit fractional vertices",
        ));
    }
    let g = ("Gamma", [0.0; 3], false);
    let list = match name {
        "fcc" => vec![
            g,
            ("X", [0.5, 0.0, 0.5], false),
            ("W", [0.5, 0.25, 0.75], false),
            ("K", [0.375, 0.375, 0.75], false),
            g,
            ("L", [0.5; 3], false),
            ("U", [0.625, 0.25, 0.625], false),
            ("W", [0.5, 0.25, 0.75], false),
            ("L", [0.5; 3], false),
            ("K", [0.375, 0.375, 0.75], false),
            ("U", [0.625, 0.25, 0.625], true),
            ("X", [0.5, 0.0, 0.5], false),
        ],
        "bcc" => vec![
            g,
            ("H", [0.5, -0.5, 0.5], false),
            ("N", [0.0, 0.0, 0.5], false),
            g,
            ("P", [0.25; 3], false),
            ("H", [0.5, -0.5, 0.5], false),
            ("P", [0.25; 3], true),
            ("N", [0.0, 0.0, 0.5], false),
        ],
        _ => vec![
            g,
            ("X", [0.0, 0.5, 0.0], false),
            ("M", [0.5, 0.5, 0.0], false),
            g,
            ("R", [0.5; 3], false),
            ("X", [0.0, 0.5, 0.0], false),
            ("M", [0.5, 0.5, 0.0], true),
            ("R", [0.5; 3], false),
        ],
    };
    Ok(list
        .into_iter()
        .map(|(label, fractional, break_before)| Vertex {
            label: label.into(),
            fractional,
            break_before,
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn skew_metric_and_breaks() {
        let cell =
            Cell::from_user_rows_bohr([[2.0, 0.0, 0.0], [1.0, 3.0, 0.0], [0.0, 0.0, 4.0]]).unwrap();
        let v = |label: &str, fractional, break_before| Vertex {
            label: label.into(),
            fractional,
            break_before,
        };
        let p = generate(
            &cell,
            &[
                v("A", [0.0; 3], false),
                v("B", [1.0, 0.0, 0.0], false),
                v("C", [0.0; 3], true),
                v("D", [0.0, 0.0, 1.0], false),
            ],
            2,
            None,
        )
        .unwrap();
        let length = std::f64::consts::PI * (10.0_f64 / 9.0).sqrt();
        assert!((p.distances[2] - length).abs() < 1e-12);
        assert!((p.distances[5] - length - std::f64::consts::FRAC_PI_2).abs() < 1e-12);
        assert_eq!(p.breaks, vec![3]);
        assert_eq!(p.distances[2], p.distances[3]);
        assert!(preset(&cell, "fcc").is_err());
    }
    #[test]
    fn geometry_checks_and_rotation() {
        for (name, rows) in [
            ("fcc", [[0., 1., 1.], [1., 0., 1.], [1., 1., 0.]]),
            ("bcc", [[-1., 1., 1.], [1., -1., 1.], [1., 1., -1.]]),
            ("sc", [[2., 0., 0.], [0., 2., 0.], [0., 0., 2.]]),
        ] {
            let cell = Cell::from_user_rows_bohr(rows).unwrap();
            assert!(preset(&cell, name).is_ok());
            let rotated = rows.map(|v| [-v[1], v[0], v[2]]);
            assert!(preset(&Cell::from_user_rows_bohr(rotated).unwrap(), name).is_ok());
            let mut distorted = rows;
            distorted[0][0] += 0.1;
            assert!(preset(&Cell::from_user_rows_bohr(distorted).unwrap(), name).is_err());
        }
    }
}
