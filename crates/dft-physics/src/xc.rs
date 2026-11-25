//! Slater exchange (alpha=2/3) and unpolarized PZ81 correlation.
use crate::{PhysicsError, volume};
use dft_units::Hartree;
use ndarray::Array3;
use std::f64::consts::PI;

pub const DENSITY_THRESHOLD: f64 = 1e-10;
pub const NEGATIVE_ROUNDOFF: f64 = 1e-14;

pub(crate) fn checked_density(n: f64) -> Result<f64, PhysicsError> {
    if !n.is_finite() || n < -NEGATIVE_ROUNDOFF {
        Err(PhysicsError::Density(n))
    } else {
        Ok(n.max(0.0))
    }
}

#[derive(Clone, Copy, Debug, Default)]
pub struct XcPoint {
    /// Exchange energy per electron (Ha).
    pub eps_x: f64,
    /// Correlation energy per electron (Ha).
    pub eps_c: f64,
    /// Exchange potential (Ha).
    pub v_x: f64,
    /// Correlation potential (Ha).
    pub v_c: f64,
}
impl XcPoint {
    pub fn eps_xc(self) -> f64 {
        self.eps_x + self.eps_c
    }
    pub fn v_xc(self) -> f64 {
        self.v_x + self.v_c
    }
}

/// Allocation-free point evaluation. At/below threshold all four values vanish.
pub fn lda_pz(n: f64) -> Result<XcPoint, PhysicsError> {
    let n = if screen_negative() && n.is_finite() {
        screened_argument(n)
    } else {
        checked_density(n)?
    };
    lda_pz_checked(n)
}

fn lda_pz_checked(n: f64) -> Result<XcPoint, PhysicsError> {
    if n <= DENSITY_THRESHOLD {
        return Ok(XcPoint::default());
    }
    let rs = (3.0 / (4.0 * PI)).cbrt() / n.cbrt();
    let eps_x = -0.75 * (3.0 / PI).cbrt() * n.cbrt();
    let (eps_c, derivative) = if rs >= 1.0 {
        let denominator = 1.0 + 1.0529 * rs.sqrt() + 0.3334 * rs;
        (
            -0.1423 / denominator,
            0.1423 * (1.0529 / (2.0 * rs.sqrt()) + 0.3334) / denominator.powi(2),
        )
    } else {
        (
            0.0311 * rs.ln() - 0.048 + 0.0020 * rs * rs.ln() - 0.0116 * rs,
            0.0311 / rs + 0.0020 * (rs.ln() + 1.0) - 0.0116,
        )
    };
    // drs/dn = -rs/(3n); d(n eps_c)/dn = eps_c - rs eps_c'/3.
    Ok(XcPoint {
        eps_x,
        eps_c,
        v_x: 4.0 * eps_x / 3.0,
        v_c: eps_c - rs * derivative / 3.0,
    })
}

pub struct XcGrid {
    pub potential: Array3<f64>,
    pub energy: Hartree,
    /// Integral n(r) v_xc(r) d^3r (Ha), for energy double counting.
    pub potential_integral: Hartree,
}

pub fn evaluate(n: &Array3<f64>, volume_bohr3: f64) -> Result<XcGrid, PhysicsError> {
    volume(volume_bohr3)?;
    if n.is_empty() {
        return Err(PhysicsError::Shape);
    }
    audit_density(n, volume_bohr3, "total", -NEGATIVE_ROUNDOFF)?;
    let mut potential = Array3::zeros(n.raw_dim());
    let mut energy = 0.0;
    let mut integral = 0.0;
    let screening = screen_negative();
    for (v, &n) in potential.iter_mut().zip(n) {
        let n = if screening {
            screened_argument(n)
        } else {
            checked_density(n)?
        };
        let p = lda_pz_checked(n)?;
        *v = p.v_xc();
        energy += n * p.eps_xc();
        integral += n * *v;
    }
    let weight = volume_bohr3 / n.len() as f64;
    energy *= weight;
    integral *= weight;
    if !energy.is_finite() || !integral.is_finite() {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(XcGrid {
        potential,
        energy: Hartree(energy),
        potential_integral: Hartree(integral),
    })
}

/// Proposed ADR 0046. This bounds a model modification, not an energy error.
pub const SCREENED_CHARGE_FRACTION: f64 = 1e-4;

pub fn screen_negative() -> bool {
    dft_numerics::switches::on(dft_numerics::switches::Id::XC_SCREEN_NEGATIVE)
}

/// XC argument only; never store this value in an electronic density.
pub fn screened_argument(n: f64) -> f64 {
    if n <= DENSITY_THRESHOLD { 0.0 } else { n }
}

/// In-memory field audit. The caller owns persistence; no I/O in evaluators.
#[derive(Clone, Debug, PartialEq)]
pub struct ScreeningRecord {
    pub stage: &'static str,
    pub channel: &'static str,
    pub points: usize,
    pub charge_electrons: f64,
    pub total_electrons: f64,
    pub minimum_bohr3: f64,
    pub minimum_grid_index: usize,
}

thread_local! {
    static AUDIT: std::cell::RefCell<(u64, Vec<ScreeningRecord>)> = const { std::cell::RefCell::new((0, Vec::new())) };
    static CONTEXT: std::cell::Cell<(&'static str, bool, Option<f64>)> = const { std::cell::Cell::new(("map/verification", false, None)) };
}

/// Dynamic context for a single XC boundary, restored even on an error.
pub struct ScreeningScope((&'static str, bool, Option<f64>));
impl Drop for ScreeningScope {
    fn drop(&mut self) {
        CONTEXT.set(self.0);
    }
}
/// `prechecked` is for distributed callers that already checked the GLOBAL
/// charge budget; never apply a per-rank charge-fraction test.
pub fn screening_scope(
    stage: &'static str,
    prechecked: bool,
    electrons: Option<f64>,
) -> ScreeningScope {
    ScreeningScope(CONTEXT.replace((stage, prechecked, electrons)))
}
pub fn screening_stage() -> &'static str {
    CONTEXT.get().0
}

pub fn screening_records() -> Vec<ScreeningRecord> {
    let id = dft_numerics::switches::RunConfig::current().run_id();
    AUDIT.with_borrow(|a| if a.0 == id { a.1.clone() } else { Vec::new() })
}
pub fn clear_screening_records() {
    AUDIT.with_borrow_mut(|a| {
        a.0 = dft_numerics::switches::RunConfig::current().run_id();
        a.1.clear();
    });
}

/// A fixed z-row partition makes the diagnostic sums independent of MPI rank
/// count. Each complete row must belong to exactly one rank. These are small
/// row summaries, never a gathered density field. Five columns per row:
/// screened count, screened |n|, signed n, minimum n, minimum global index.
pub fn screening_rows(n: &[f64], indices: &[usize], dims: [usize; 3]) -> Vec<f64> {
    let mut rows = vec![0.; 5 * dims[0] * dims[1]];
    for (&n, &i) in n.iter().zip(indices) {
        let row = &mut rows[5 * (i / dims[2])..5 * (i / dims[2]) + 5];
        if n <= DENSITY_THRESHOLD {
            row[0] += 1.;
            row[1] += n.abs();
        }
        row[2] += n;
        if i % dims[2] == 0 || n < row[3] || !n.is_finite() {
            row[3] = n;
            row[4] = i as f64;
        }
    }
    rows
}

pub fn finish_screening(
    rows: &[f64],
    weight: f64,
    channel: &'static str,
    electrons: Option<f64>,
) -> Result<ScreeningRecord, PhysicsError> {
    let mut record = ScreeningRecord {
        stage: screening_stage(),
        channel,
        points: 0,
        charge_electrons: 0.,
        total_electrons: 0.,
        minimum_bohr3: f64::INFINITY,
        minimum_grid_index: 0,
    };
    for row in rows.as_chunks::<5>().0.iter() {
        record.points += row[0] as usize;
        record.charge_electrons += row[1];
        record.total_electrons += row[2];
        if row[3] < record.minimum_bohr3 || !row[3].is_finite() {
            record.minimum_bohr3 = row[3];
            record.minimum_grid_index = row[4] as usize;
        }
    }
    record.charge_electrons *= weight;
    record.total_electrons = electrons.unwrap_or(record.total_electrons * weight);
    let enabled = screen_negative();
    let bad = !record.minimum_bohr3.is_finite()
        || !record.total_electrons.is_finite()
        || (enabled
            && record.charge_electrons > SCREENED_CHARGE_FRACTION * record.total_electrons.max(0.));
    if enabled {
        let id = dft_numerics::switches::RunConfig::current().run_id();
        AUDIT.with_borrow_mut(|a| {
            if a.0 != id {
                a.0 = id;
                a.1.clear();
            }
            a.1.push(record.clone());
        });
    }
    if bad {
        return Err(PhysicsError::Xc(format!(
            "{}: XC channel {} minimum {:.17e} Bohr^-3 at grid point {}; screened charge {:.17e} e exceeds fraction {} of {:.17e} e, or nonfinite input",
            record.stage,
            channel,
            record.minimum_bohr3,
            record.minimum_grid_index,
            record.charge_electrons,
            SCREENED_CHARGE_FRACTION,
            record.total_electrons,
        )));
    }
    Ok(record)
}

/// Check an ordinary complete field. Legacy off-mode arithmetic remains intact.
pub fn audit_density(
    n: &Array3<f64>,
    volume: f64,
    channel: &'static str,
    negative_limit: f64,
) -> Result<(), PhysicsError> {
    if CONTEXT.get().1 {
        return Ok(());
    }
    if !screen_negative() && n.iter().all(|v| v.is_finite() && *v >= negative_limit) {
        return Ok(());
    }
    let dims = [n.shape()[0], n.shape()[1], n.shape()[2]];
    let values: Vec<_> = n.iter().copied().collect();
    let rows = screening_rows(&values, &(0..n.len()).collect::<Vec<_>>(), dims);
    let record = finish_screening(&rows, volume / n.len() as f64, channel, CONTEXT.get().2)?;
    if !screen_negative() && record.minimum_bohr3 < negative_limit {
        return Err(PhysicsError::Xc(format!(
            "{}: XC channel {} minimum {:.17e} Bohr^-3 at grid point {} (screening off)",
            record.stage, channel, record.minimum_bohr3, record.minimum_grid_index
        )));
    }
    Ok(())
}

pub fn screening_prechecked() -> bool {
    CONTEXT.get().1
}

/// Preserve the first rank's detailed failure on every rank before continuing.
pub fn agree_error<K: dft_mpi::Collective>(
    comm: &K,
    error: Option<String>,
) -> Result<(), PhysicsError> {
    let mut failures = [0.];
    comm.tree_sum_f64(&[f64::from(error.is_some())], &mut failures)
        .map_err(|e| PhysicsError::Xc(e.to_string()))?;
    if failures[0] == 0. {
        return Ok(());
    }
    for root in 0..comm.size() {
        let mut length = [if root == comm.rank() {
            error.as_ref().map_or(0, String::len) as f64
        } else {
            0.
        }];
        comm.bcast_f64(&mut length, root)
            .map_err(|e| PhysicsError::Xc(e.to_string()))?;
        if length[0] > 0. {
            let mut bytes = if root == comm.rank() {
                error
                    .as_ref()
                    .unwrap()
                    .bytes()
                    .map(f64::from)
                    .collect::<Vec<_>>()
            } else {
                vec![0.; length[0] as usize]
            };
            comm.bcast_f64(&mut bytes, root)
                .map_err(|e| PhysicsError::Xc(e.to_string()))?;
            return Err(PhysicsError::Xc(
                String::from_utf8_lossy(&bytes.iter().map(|&v| v as u8).collect::<Vec<_>>())
                    .into_owned(),
            ));
        }
    }
    Ok(())
}
