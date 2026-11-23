//! In-memory reader for the four-record QE charge-density.dat format specified
//! in the own-words density-mixing note, section 6.2. No QE source was used.
//! Fortran sequential records: 4-byte markers, int32 logical/ngm/nspin,
//! float64 reciprocal columns, interleaved int32 Miller triples, complex128 rho.
//! Both byte orders are recognized; HDF5, 8-byte markers, spin and non-Gamma
//! files fail explicitly. The CLI owns filesystem access.
use super::{
    Error,
    density::{g_squared, half},
    invalid,
};
use dft_physics::data::Complex64 as C;
use std::collections::HashMap;
#[derive(Clone, Debug)]
pub struct QeDensity {
    /// Matrix with reciprocal vectors as columns, including 2 pi, bohr^-1.
    pub reciprocal: [[f64; 3]; 3],
    pub miller: Vec<[i32; 3]>,
    pub rho: Vec<C>,
}
impl QeDensity {
    pub fn decode(bytes: &[u8]) -> Result<Self, Error> {
        let first: [u8; 4] = bytes
            .get(..4)
            .ok_or_else(|| invalid("truncated QE header"))?
            .try_into()
            .expect("four");
        let little = if u32::from_le_bytes(first) == 12 {
            true
        } else if u32::from_be_bytes(first) == 12 {
            false
        } else {
            return Err(invalid(
                "QE record 1 must have 12 bytes (four-byte Fortran markers; HDF5 unsupported)",
            ));
        };
        let int = |b: &[u8]| {
            let a = b.try_into().expect("int32 field");
            if little {
                i32::from_le_bytes(a)
            } else {
                i32::from_be_bytes(a)
            }
        };
        let float = |b: &[u8]| {
            let a = b.try_into().expect("float64 field");
            if little {
                f64::from_le_bytes(a)
            } else {
                f64::from_be_bytes(a)
            }
        };
        let mut offset = 0usize;
        let mut record = |expected: usize| -> Result<&[u8], Error> {
            let marker = bytes
                .get(offset..offset + 4)
                .ok_or_else(|| invalid("truncated QE record marker"))?;
            if int(marker) < 0 || int(marker) as usize != expected {
                return Err(invalid(format!(
                    "QE record length: expected {expected}, found {}",
                    int(marker)
                )));
            }
            let start = offset + 4;
            let end = start
                .checked_add(expected)
                .ok_or_else(|| invalid("QE length overflow"))?;
            let payload = bytes
                .get(start..end)
                .ok_or_else(|| invalid("truncated QE record"))?;
            if bytes.get(end..end + 4) != Some(marker) {
                return Err(invalid("QE trailing record marker mismatch"));
            }
            offset = end + 4;
            Ok(payload)
        };
        let h = record(12)?;
        if ![1, -1].contains(&int(&h[..4])) || int(&h[8..]) != 1 {
            return Err(invalid("QE reader requires gamma_only=true, nspin=1"));
        }
        let n = usize::try_from(int(&h[4..8])).map_err(|_| invalid("QE ngm is negative"))?;
        if n == 0 || n > bytes.len() / 28 {
            return Err(invalid("QE ngm exceeds file length"));
        }
        let b = record(72)?;
        let reciprocal = std::array::from_fn(|i| {
            std::array::from_fn(|j| float(&b[8 * (3 * j + i)..8 * (3 * j + i + 1)]))
        });
        if reciprocal.iter().flatten().any(|v| !v.is_finite()) {
            return Err(invalid("nonfinite QE reciprocal basis"));
        }
        let m = record(12 * n)?;
        let miller: Vec<[i32; 3]> = m
            .as_chunks::<12>()
            .0
            .iter()
            .map(|r| std::array::from_fn(|i| int(&r[4 * i..4 * i + 4])))
            .collect();
        let r = record(16 * n)?;
        let rho: Vec<_> = r
            .as_chunks::<16>()
            .0
            .iter()
            .map(|v| C::new(float(&v[..8]), float(&v[8..])))
            .collect();
        if offset != bytes.len() {
            return Err(invalid("unexpected trailing QE records"));
        }
        let mut seen = std::collections::HashSet::new();
        if miller.iter().any(|&m| !half(m) || !seen.insert(m)) || !seen.contains(&[0; 3]) {
            return Err(invalid(
                "QE Miller list is not a unique Gamma half sphere including G=0",
            ));
        }
        if rho.iter().any(|v| !v.re.is_finite() || !v.im.is_finite()) {
            return Err(invalid("nonfinite QE density"));
        }
        Ok(Self {
            reciprocal,
            miller,
            rho,
        })
    }
    /// Match every coefficient by Miller triple, independently of shell ties.
    pub fn aligned(&self, order: &[[i32; 3]]) -> Result<Vec<C>, Error> {
        if order.len() != self.miller.len() {
            return Err(invalid("QE Miller membership/count mismatch"));
        }
        let map: HashMap<_, _> = self
            .miller
            .iter()
            .copied()
            .zip(self.rho.iter().copied())
            .collect();
        let mut seen = std::collections::HashSet::new();
        order
            .iter()
            .map(|m| {
                if !seen.insert(m) {
                    return Err(invalid("duplicate target Miller index"));
                }
                map.get(m)
                    .copied()
                    .ok_or_else(|| invalid(format!("missing QE Miller index {m:?}")))
            })
            .collect()
    }
    pub fn validate_cell(&self, reciprocal: [[f64; 3]; 3]) -> Result<(), Error> {
        for (a, b) in self
            .reciprocal
            .iter()
            .flatten()
            .zip(reciprocal.iter().flatten())
        {
            if (a - b).abs() > 1e-11 * a.abs().max(b.abs()).max(1.0) {
                return Err(invalid(
                    "QE/input reciprocal cell mismatch (bohr^-1, including 2 pi)",
                ));
            }
        }
        Ok(())
    }
    /// Use the file's membership (cutoff boundary rounding is external data),
    /// sorted in BAANDI's norm/Miller order, rather than QE shell tie order.
    pub fn order(&self) -> Vec<[i32; 3]> {
        let mut m = self.miller.clone();
        m.sort_by(|a, b| {
            g_squared(self.reciprocal, *a)
                .total_cmp(&g_squared(self.reciprocal, *b))
                .then(a.cmp(b))
        });
        m
    }
}
#[derive(Clone, Copy, Debug)]
pub struct PrintedAccuracy {
    pub ry: f64,
    pub rounding_ry: f64,
    pub upper_bound: bool,
}
/// Final printed estimate in a k-step stdout. The '<' form is a bound, not an
/// exact rounded value. Missing estimates are reported unavailable by callers.
pub fn printed_accuracy(stdout: &str) -> Option<PrintedAccuracy> {
    stdout
        .lines()
        .filter_map(|line| {
            let (_, tail) = line.split_once("estimated scf accuracy")?;
            let number = tail.split_whitespace().find_map(|s| {
                s.replace(['D', 'd'], "e")
                    .parse::<f64>()
                    .ok()
                    .map(|v| (s, v))
            })?;
            if !number.1.is_finite() || number.1 < 0.0 {
                return None;
            }
            let text = number.0.to_lowercase().replace('d', "e");
            let (mantissa, exp) = text.split_once('e').map_or((text.as_str(), 0), |(m, e)| {
                (m, e.parse::<i32>().ok().unwrap_or(0))
            });
            let digits = mantissa.split_once('.').map_or(0, |(_, d)| d.len() as i32);
            Some(PrintedAccuracy {
                ry: number.1,
                rounding_ry: 0.5 * 10f64.powi(exp - digits),
                upper_bound: tail.contains('<'),
            })
        })
        .next_back()
}
#[derive(Clone, Copy, Debug)]
pub struct Difference {
    pub max: f64,
    pub rms: f64,
    pub relative_max: f64,
    pub relative_rms: f64,
    pub scale: f64,
}
pub fn difference(a: &[C], b: &[C], indices: &[usize]) -> Difference {
    let mut max = 0.0_f64;
    let mut sum = 0.0;
    let mut scale = 0.0_f64;
    let mut norm = 0.0;
    for &i in indices {
        let d = (a[i] - b[i]).norm();
        max = max.max(d);
        sum += d * d;
        scale = scale.max(b[i].norm());
        norm += b[i].norm_sqr();
    }
    Difference {
        max,
        rms: (sum / indices.len().max(1) as f64).sqrt(),
        relative_max: if scale > 0.0 { max / scale } else { max },
        relative_rms: if norm > 0.0 {
            (sum / norm).sqrt()
        } else {
            sum.sqrt()
        },
        scale,
    }
}
