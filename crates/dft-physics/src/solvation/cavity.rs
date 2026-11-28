//! Density-dependent implicit-solvation cavity fields, skew-adjoint Nyquist-zeroed
//! reciprocal vectors, Gaussian ionic pseudo-charges, and exact discrete ionic
//! solvation forces in Hartree atomic units.

use crate::solvation::{CavityModel, IonicPseudoCharge, SolvationError, SolvationParams};
use dft_geometry::Cell;
use dft_mpi::Collective;
use num_complex::Complex64;
use std::f64::consts::{PI, SQRT_2, TAU};

/// Prevalidated Gaussian ionic pseudo-charge converted to Cartesian coordinates.
struct PreparedIon {
    q_over_omega: f64,
    neg_quarter_rc2: f64,
    tau_bohr: [f64; 3],
    alpha_z_sm: f64,
}

/// Returns `true` iff any axis `a` of the C-order `(i, j, k)` index lies on an
/// even-axis Nyquist plane (`dims[a].is_multiple_of(2) && ijk[a] == dims[a] / 2`).
pub fn is_even_nyquist(dims: [usize; 3], global_index: usize) -> bool {
    if dims.contains(&0) {
        return false;
    }
    let ijk = [
        global_index / (dims[1] * dims[2]),
        (global_index / dims[2]) % dims[1],
        global_index % dims[2],
    ];
    (0..3).any(|a| dims[a].is_multiple_of(2) && ijk[a] == dims[a] / 2)
}

/// Cartesian reciprocal vector (Bohr^-1) for BAANDI's real, skew-adjoint spectral
/// derivative: zero each even-axis Nyquist Miller component BEFORE multiplying by
/// `recip_cols` (`G[axis] = sum_a recip_cols[axis][a] * m[a]`).
pub fn skew_adjoint_g_cart(
    dims: [usize; 3],
    global_index: usize,
    recip_cols: [[f64; 3]; 3],
) -> [f64; 3] {
    if dims.contains(&0) {
        return [0.0; 3];
    }
    let ijk = [
        global_index / (dims[1] * dims[2]),
        (global_index / dims[2]) % dims[1],
        global_index % dims[2],
    ];
    let miller: [f64; 3] = std::array::from_fn(|a| {
        if dims[a].is_multiple_of(2) && ijk[a] == dims[a] / 2 {
            0.0
        } else if ijk[a] <= (dims[a] - 1) / 2 {
            ijk[a] as f64
        } else {
            ijk[a] as f64 - dims[a] as f64
        }
    });
    std::array::from_fn(|axis| (0..3).map(|a| recip_cols[axis][a] * miller[a]).sum())
}

/// Evaluate switching function `S(n)`, relative permittivity `eps(n)`, LPB
/// screening coefficient `eps_b * kappa^2 * S(n)`, and their exact derivatives
/// with respect to `n` on caller-owned plane slices without heap allocation.
#[allow(clippy::too_many_arguments)]
pub fn evaluate_cavity_fields(
    density_planes: &[f64],
    params: &SolvationParams,
    s: &mut [f64],
    ds_dn: &mut [f64],
    eps: &mut [f64],
    deps_dn: &mut [f64],
    kappa2: &mut [f64],
    dkappa2_dn: &mut [f64],
) -> Result<(), SolvationError> {
    params.validate()?;
    let len = density_planes.len();
    if s.len() != len
        || ds_dn.len() != len
        || eps.len() != len
        || deps_dn.len() != len
        || kappa2.len() != len
        || dkappa2_dn.len() != len
    {
        return Err(SolvationError::Shape);
    }
    if density_planes.iter().any(|n| !n.is_finite()) {
        return Err(SolvationError::NonFinite);
    }

    let eps_b = params.epsilon_bulk;
    let eps_delta = eps_b - 1.0;
    let k2_bulk = eps_b * params.kappa_inv_bohr * params.kappa_inv_bohr;
    if !k2_bulk.is_finite() {
        return Err(SolvationError::NonFinite);
    }
    let floor = params.density_floor_bohr3;

    match params.cavity {
        CavityModel::VaspSol { nc_bohr3, sigma } => {
            let sigma_sqrt2 = sigma * SQRT_2;
            let sigma_sqrt2pi = sigma * (2.0 * PI).sqrt();
            for (idx, &n) in density_planes.iter().enumerate() {
                let (sv, dsv) = if n <= floor {
                    (1.0, 0.0)
                } else {
                    let x = (n / nc_bohr3).ln() / sigma_sqrt2;
                    let sv = 0.5 * libm::erfc(x);
                    let dsv = -(-x * x).exp() / (n * sigma_sqrt2pi);
                    (sv, dsv)
                };
                s[idx] = sv;
                ds_dn[idx] = dsv;
                eps[idx] = 1.0 + eps_delta * sv;
                deps_dn[idx] = eps_delta * dsv;
                kappa2[idx] = k2_bulk * sv;
                dkappa2_dn[idx] = k2_bulk * dsv;
            }
        }
        CavityModel::Sccs {
            rho_min_bohr3,
            rho_max_bohr3,
            log_epsilon,
        } => {
            let denom = (rho_max_bohr3 / rho_min_bohr3).ln();
            if !denom.is_finite() || denom <= 0.0 {
                return Err(SolvationError::NonFinite);
            }
            let ln_eb = eps_b.ln();
            for (idx, &n) in density_planes.iter().enumerate() {
                let (sv, dsv) = if n <= floor || n <= rho_min_bohr3 {
                    (1.0, 0.0)
                } else if n >= rho_max_bohr3 {
                    (0.0, 0.0)
                } else {
                    let t = (rho_max_bohr3 / n).ln() / denom;
                    let sv = t - (TAU * t).sin() / TAU;
                    let dsv = -(1.0 - (TAU * t).cos()) / (n * denom);
                    (sv, dsv)
                };
                let (ev, dev) = if log_epsilon {
                    let ev = (ln_eb * sv).exp();
                    (ev, ev * ln_eb * dsv)
                } else {
                    (1.0 + eps_delta * sv, eps_delta * dsv)
                };
                s[idx] = sv;
                ds_dn[idx] = dsv;
                eps[idx] = ev;
                deps_dn[idx] = dev;
                kappa2[idx] = k2_bulk * sv;
                dkappa2_dn[idx] = k2_bulk * dsv;
            }
        }
    }

    Ok(())
}

fn validate_grid_positions(
    dims: [usize; 3],
    local_positions: &[usize],
    slice_len: usize,
) -> Result<(), SolvationError> {
    if dims.contains(&0) || slice_len != local_positions.len() {
        return Err(SolvationError::Shape);
    }
    let total_points = dims[0]
        .checked_mul(dims[1])
        .and_then(|n| n.checked_mul(dims[2]))
        .ok_or(SolvationError::Shape)?;
    if local_positions.iter().any(|&idx| idx >= total_points) {
        return Err(SolvationError::Shape);
    }
    Ok(())
}

fn prepare_ions(
    cell: &Cell,
    ions: &[IonicPseudoCharge],
) -> Result<Vec<PreparedIon>, SolvationError> {
    let a = cell.lattice_columns_bohr();
    let omega = cell.volume_bohr3();
    let mut prepared = Vec::with_capacity(ions.len());
    for ion in ions {
        if !ion.width_bohr.is_finite() || ion.width_bohr <= 0.0 {
            return Err(SolvationError::InvalidParameter(
                "ionic pseudo-charge width_bohr must be finite and > 0.0",
            ));
        }
        if !ion.charge.is_finite() || ion.fractional_position.iter().any(|v| !v.is_finite()) {
            return Err(SolvationError::NonFinite);
        }
        let f = ion.fractional_position;
        let tau_bohr = [
            a[0][0] * f[0] + a[0][1] * f[1] + a[0][2] * f[2],
            a[1][0] * f[0] + a[1][1] * f[1] + a[1][2] * f[2],
            a[2][0] * f[0] + a[2][1] * f[1] + a[2][2] * f[2],
        ];
        let q_over_omega = ion.charge / omega;
        let neg_quarter_rc2 = -0.25 * ion.width_bohr * ion.width_bohr;
        let alpha_z_sm = PI * ion.charge * ion.width_bohr * ion.width_bohr / omega;
        if tau_bohr.iter().any(|v| !v.is_finite())
            || !q_over_omega.is_finite()
            || !neg_quarter_rc2.is_finite()
            || !alpha_z_sm.is_finite()
        {
            return Err(SolvationError::NonFinite);
        }
        prepared.push(PreparedIon {
            q_over_omega,
            neg_quarter_rc2,
            tau_bohr,
            alpha_z_sm,
        });
    }
    Ok(prepared)
}

/// Assemble the Gaussian ionic pseudo-charge Fourier coefficients `rho_ion(G)`
/// on `local_positions` (zeroing even-axis Nyquist planes to preserve real density
/// and exact discrete ionic forces) and return `sum_I alpha_z,I^sm`.
pub fn build_ionic_pseudo_charge_g(
    cell: &Cell,
    dims: [usize; 3],
    local_positions: &[usize],
    ions: &[IonicPseudoCharge],
    out_g: &mut [Complex64],
) -> Result<f64, SolvationError> {
    validate_grid_positions(dims, local_positions, out_g.len())?;
    let prepared = prepare_ions(cell, ions)?;
    let sum_alpha_z_sm: f64 = prepared.iter().map(|ion| ion.alpha_z_sm).sum();
    if !sum_alpha_z_sm.is_finite() {
        return Err(SolvationError::NonFinite);
    }
    if prepared.is_empty() {
        out_g.fill(Complex64::default());
        return Ok(0.0);
    }

    let b = cell.reciprocal_columns_inverse_bohr();
    for (out, &g_idx) in out_g.iter_mut().zip(local_positions) {
        if is_even_nyquist(dims, g_idx) {
            *out = Complex64::default();
            continue;
        }
        let g = skew_adjoint_g_cart(dims, g_idx, b);
        let g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
        let mut sum = Complex64::default();
        for ion in &prepared {
            let g_dot_tau =
                g[0] * ion.tau_bohr[0] + g[1] * ion.tau_bohr[1] + g[2] * ion.tau_bohr[2];
            let amp = ion.q_over_omega * (ion.neg_quarter_rc2 * g2).exp();
            sum += Complex64::from_polar(amp, -g_dot_tau);
        }
        if !sum.re.is_finite() || !sum.im.is_finite() {
            return Err(SolvationError::NonFinite);
        }
        *out = sum;
    }

    Ok(sum_alpha_z_sm)
}

/// Exact discrete derivative of the electrostatic solvation energy with respect to
/// Gaussian ionic pseudo-charge positions (Ha/Bohr), reduced across ranks via
/// `comm.tree_sum_f64`.
pub fn ionic_solvation_forces<K: Collective>(
    comm: &K,
    cell: &Cell,
    dims: [usize; 3],
    local_positions: &[usize],
    ions: &[IonicPseudoCharge],
    phi_pol_g: &[Complex64],
) -> Result<Vec<[f64; 3]>, SolvationError> {
    validate_grid_positions(dims, local_positions, phi_pol_g.len())?;
    if phi_pol_g
        .iter()
        .any(|z| !z.re.is_finite() || !z.im.is_finite())
    {
        return Err(SolvationError::NonFinite);
    }
    let prepared = prepare_ions(cell, ions)?;
    if prepared.is_empty() {
        return Ok(Vec::new());
    }

    let b = cell.reciprocal_columns_inverse_bohr();
    let omega = cell.volume_bohr3();
    let mut local_forces = vec![0.0; 3 * prepared.len()];

    for (&g_idx, &phi_g) in local_positions.iter().zip(phi_pol_g) {
        if is_even_nyquist(dims, g_idx) {
            continue;
        }
        let g = skew_adjoint_g_cart(dims, g_idx, b);
        let g2 = g[0] * g[0] + g[1] * g[1] + g[2] * g[2];
        let phi_conj = phi_g.conj();
        for (ion_idx, ion) in prepared.iter().enumerate() {
            let g_dot_tau =
                g[0] * ion.tau_bohr[0] + g[1] * ion.tau_bohr[1] + g[2] * ion.tau_bohr[2];
            let amp = ion.q_over_omega * (ion.neg_quarter_rc2 * g2).exp();
            let product = phi_conj * Complex64::from_polar(amp, -g_dot_tau);
            let base = 3 * ion_idx;
            local_forces[base] -= omega * g[0] * product.im;
            local_forces[base + 1] -= omega * g[1] * product.im;
            local_forces[base + 2] -= omega * g[2] * product.im;
        }
    }

    let mut global_forces = vec![0.0; local_forces.len()];
    comm.tree_sum_f64(&local_forces, &mut global_forces)?;
    if global_forces.iter().any(|v| !v.is_finite()) {
        return Err(SolvationError::NonFinite);
    }

    Ok(global_forces.as_chunks::<3>().0.to_vec())
}

#[cfg(test)]
mod tests {
    use super::*;
    use dft_mpi::{Comm, simulate};

    fn skewed_cell() -> Cell {
        Cell::from_user_rows_bohr([[8.5, 0.7, 0.0], [-0.4, 9.2, 0.5], [0.3, -0.6, 10.1]]).unwrap()
    }

    #[test]
    fn nyquist_and_skew_adjoint_g_vectors_match_rules() {
        let cell = skewed_cell();
        let b = cell.reciprocal_columns_inverse_bohr();
        let dims = [4, 5, 6];
        let n_tot = dims[0] * dims[1] * dims[2];

        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    let idx = (i * dims[1] + j) * dims[2] + k;
                    let expected_nyq = i == 2 || k == 3;
                    assert_eq!(is_even_nyquist(dims, idx), expected_nyq);

                    let g = skew_adjoint_g_cart(dims, idx, b);
                    let neg_i = (dims[0] - i) % dims[0];
                    let neg_j = (dims[1] - j) % dims[1];
                    let neg_k = (dims[2] - k) % dims[2];
                    let neg_idx = (neg_i * dims[1] + neg_j) * dims[2] + neg_k;
                    let g_neg = skew_adjoint_g_cart(dims, neg_idx, b);
                    for axis in 0..3 {
                        assert!((g[axis] + g_neg[axis]).abs() < 1e-13);
                    }
                }
            }
        }
        assert_eq!(n_tot, 120);

        // Even Nyquist on axis 0 (i=2) with nonzero j=1, k=0 must zero m[0] before B*m.
        let idx = (2 * dims[1] + 1) * dims[2];
        let g = skew_adjoint_g_cart(dims, idx, b);
        for axis in 0..3 {
            assert!((g[axis] - b[axis][1]).abs() < 1e-14);
        }
    }

    #[test]
    fn cavity_fields_match_finite_differences_and_limits() {
        let models = [
            SolvationParams {
                kappa_inv_bohr: 0.15,
                ..SolvationParams::vaspsol_water()
            },
            SolvationParams {
                kappa_inv_bohr: 0.12,
                ..SolvationParams::sccs_water()
            },
            SolvationParams {
                kappa_inv_bohr: 0.12,
                cavity: CavityModel::Sccs {
                    rho_min_bohr3: 1.0e-4,
                    rho_max_bohr3: 1.5e-3,
                    log_epsilon: false,
                },
                ..SolvationParams::sccs_water()
            },
        ];

        for params in models {
            let densities = [-1.0e-6, 0.0, 1.0e-18, 2.0e-4, 6.0e-4, 1.2e-3, 3.0e-3, 0.05];
            let len = densities.len();
            let mut s = vec![0.0; len];
            let mut ds = vec![0.0; len];
            let mut eps = vec![0.0; len];
            let mut deps = vec![0.0; len];
            let mut k2 = vec![0.0; len];
            let mut dk2 = vec![0.0; len];
            evaluate_cavity_fields(
                &densities, &params, &mut s, &mut ds, &mut eps, &mut deps, &mut k2, &mut dk2,
            )
            .unwrap();

            // Below density floor: smooth S = 1, S' = 0 without NaN.
            for idx in 0..3 {
                assert_eq!(s[idx], 1.0);
                assert_eq!(ds[idx], 0.0);
                assert_eq!(deps[idx], 0.0);
                assert_eq!(dk2[idx], 0.0);
            }

            // Central finite-difference check on interior points.
            for n0 in [2.5e-4, 7.0e-4, 1.1e-3] {
                let h = n0 * 1.0e-6;
                let probe = [n0 - h, n0, n0 + h];
                let mut sp = [0.0; 3];
                let mut dsp = [0.0; 3];
                let mut ep = [0.0; 3];
                let mut dep = [0.0; 3];
                let mut kp = [0.0; 3];
                let mut dkp = [0.0; 3];
                evaluate_cavity_fields(
                    &probe, &params, &mut sp, &mut dsp, &mut ep, &mut dep, &mut kp, &mut dkp,
                )
                .unwrap();

                let fd_s = (sp[2] - sp[0]) / (2.0 * h);
                let fd_e = (ep[2] - ep[0]) / (2.0 * h);
                let fd_k = (kp[2] - kp[0]) / (2.0 * h);
                assert!((dsp[1] - fd_s).abs() < 1.0e-7 * fd_s.abs().max(1.0));
                assert!((dep[1] - fd_e).abs() < 1.0e-7 * fd_e.abs().max(1.0));
                assert!((dkp[1] - fd_k).abs() < 1.0e-7 * fd_k.abs().max(1.0));
            }
        }

        // Exact vacuum limit.
        let vac = SolvationParams::vacuum();
        let n = [0.0, 1.0e-3, 1.0e-1];
        let mut s = [0.0; 3];
        let mut ds = [0.0; 3];
        let mut eps = [0.0; 3];
        let mut deps = [0.0; 3];
        let mut k2 = [0.0; 3];
        let mut dk2 = [0.0; 3];
        evaluate_cavity_fields(
            &n, &vac, &mut s, &mut ds, &mut eps, &mut deps, &mut k2, &mut dk2,
        )
        .unwrap();
        assert_eq!(eps, [1.0; 3]);
        assert_eq!(deps, [0.0; 3]);
        assert_eq!(k2, [0.0; 3]);
        assert_eq!(dk2, [0.0; 3]);
    }

    #[test]
    fn ionic_forces_match_finite_difference_and_mpi_distribution() {
        let cell = skewed_cell();
        let dims = [4, 5, 6];
        let n_tot = dims[0] * dims[1] * dims[2];
        let all_positions: Vec<usize> = (0..n_tot).collect();

        let ions = vec![
            IonicPseudoCharge {
                charge: 4.0,
                fractional_position: [0.17, 0.31, 0.53],
                width_bohr: 0.75,
            },
            IonicPseudoCharge {
                charge: 1.0,
                fractional_position: [0.62, 0.19, 0.81],
                width_bohr: 0.5,
            },
        ];

        // Construct a Hermitian-symmetric polarization potential phi_pol(G).
        let mut phi_pol_g = vec![Complex64::default(); n_tot];
        for i in 0..dims[0] {
            for j in 0..dims[1] {
                for k in 0..dims[2] {
                    let idx = (i * dims[1] + j) * dims[2] + k;
                    let neg_i = (dims[0] - i) % dims[0];
                    let neg_j = (dims[1] - j) % dims[1];
                    let neg_k = (dims[2] - k) % dims[2];
                    let neg_idx = (neg_i * dims[1] + neg_j) * dims[2] + neg_k;
                    if idx == neg_idx {
                        phi_pol_g[idx] = Complex64::new(0.13 * (idx as f64 + 1.0).sin(), 0.0);
                    } else if idx < neg_idx {
                        let re = 0.07 * (0.31 * idx as f64).cos();
                        let im = -0.05 * (0.23 * idx as f64).sin();
                        phi_pol_g[idx] = Complex64::new(re, im);
                        phi_pol_g[neg_idx] = Complex64::new(re, -im);
                    }
                }
            }
        }

        let mut rho_g = vec![Complex64::default(); n_tot];
        let alpha_sum =
            build_ionic_pseudo_charge_g(&cell, dims, &all_positions, &ions, &mut rho_g).unwrap();
        let omega = cell.volume_bohr3();
        let expected_alpha = PI * (4.0 * 0.75 * 0.75 + 1.0 * 0.5 * 0.5) / omega;
        assert!((alpha_sum - expected_alpha).abs() < 1e-14);
        assert!((rho_g[0].re - 5.0 / omega).abs() < 1e-14);
        assert!(rho_g[0].im.abs() < 1e-14);
        for &idx in &all_positions {
            if is_even_nyquist(dims, idx) {
                assert_eq!(rho_g[idx], Complex64::default());
            }
        }

        let analytic =
            ionic_solvation_forces(&Comm::Single, &cell, dims, &all_positions, &ions, &phi_pol_g)
                .unwrap();

        // Central finite-difference check of E = Omega * sum_G Re(conj(phi_pol(G)) * rho_ion(G))
        // with respect to Cartesian displacements of each ion.
        let h = 1.0e-6;
        let energy_for = |shifted: &[IonicPseudoCharge]| -> f64 {
            let mut rg = vec![Complex64::default(); n_tot];
            build_ionic_pseudo_charge_g(&cell, dims, &all_positions, shifted, &mut rg).unwrap();
            omega
                * phi_pol_g
                    .iter()
                    .zip(&rg)
                    .map(|(p, r)| (p.conj() * r).re)
                    .sum::<f64>()
        };

        for ion_idx in 0..ions.len() {
            let tau0 = cell.fractional_to_cartesian_bohr(ions[ion_idx].fractional_position);
            for axis in 0..3 {
                let mut plus = ions.clone();
                let mut minus = ions.clone();
                let mut tau_p = tau0;
                let mut tau_m = tau0;
                tau_p[axis] += h;
                tau_m[axis] -= h;
                plus[ion_idx].fractional_position = cell.cartesian_bohr_to_fractional(tau_p);
                minus[ion_idx].fractional_position = cell.cartesian_bohr_to_fractional(tau_m);
                let fd_force = -(energy_for(&plus) - energy_for(&minus)) / (2.0 * h);
                assert!((analytic[ion_idx][axis] - fd_force).abs() < 1.0e-8);
            }
        }

        // Verify 4-rank MPI distribution reproduces the single-rank forces.
        let distributed = simulate(4, |comm| {
            let rank = comm.rank();
            let size = comm.size();
            let local_pos: Vec<usize> = all_positions
                .iter()
                .copied()
                .filter(|&idx| idx % size == rank)
                .collect();
            let local_phi: Vec<Complex64> = local_pos.iter().map(|&idx| phi_pol_g[idx]).collect();
            ionic_solvation_forces(comm, &cell, dims, &local_pos, &ions, &local_phi).unwrap()
        });
        for rank_forces in distributed {
            for (a, b) in analytic.iter().zip(&rank_forces) {
                for axis in 0..3 {
                    assert!((a[axis] - b[axis]).abs() < 1.0e-13);
                }
            }
        }
    }
}
