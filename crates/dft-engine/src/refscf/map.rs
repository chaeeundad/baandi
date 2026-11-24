//! Owned-row map for the reference loop. Explicit kernels and call-owned
//! options only: no ScfProblem, session, retained-H or legacy solver dispatcher.
use super::numerical;
use crate::{Functional, PhaseTimings, ScfError, energy::EnergyAudit, refstart::StartState};
use dft_geometry::Cell;
use dft_mpi::Collective;
use dft_physics::data::{Array3, Complex64 as C};
use dft_solvers::ref_davidson;
use std::{f64::consts::PI, time::Instant};

pub(super) struct Fields {
    pub hartree: f64,
    pub xc: f64,
    pub xc_integral: f64,
    pub hxc: Vec<f64>,
    pub vh: Vec<f64>,
}
pub(super) struct Map {
    pub solve: ref_davidson::Solution,
    pub density: Vec<f64>,
    pub residuals: Vec<f64>,
    pub energy: EnergyAudit,
    pub timings: PhaseTimings,
}
fn sum<K: Collective>(comm: &K, values: &[f64]) -> Result<Vec<f64>, ScfError> {
    let mut out = vec![0.; values.len()];
    comm.tree_sum_f64(values, &mut out).map_err(numerical)?;
    Ok(out)
}
/// Full-box Hartree and plane-local LDA, in the same convention as refstart.
/// All local errors are agreed before the next transform/collective.
pub(super) fn fields<K: Collective>(
    state: &mut StartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    density: &[f64],
) -> Result<Fields, ScfError> {
    let dims = state.grid.dims();
    let dv = cell.volume_bohr3() / dims.iter().product::<usize>() as f64;
    let mut coeff = state.fields.forward(comm, density)?;
    let b = cell.reciprocal_columns_inverse_bohr();
    for (v, &g) in coeff.iter_mut().zip(state.fields.reciprocal_indices()) {
        let ijk = [g / (dims[1] * dims[2]), g / dims[2] % dims[1], g % dims[2]];
        let m = std::array::from_fn(|a| {
            if ijk[a] <= (dims[a] - 1) / 2 {
                ijk[a] as i32
            } else {
                ijk[a] as i32 - dims[a] as i32
            }
        });
        let g2 = super::density::g_squared(b, m);
        *v = if g == 0 {
            C::default()
        } else {
            *v * (4. * PI / g2)
        };
    }
    let mut hxc = state.fields.inverse(comm, &coeff)?;
    let eh = 0.5 * dv * density.iter().zip(&hxc).map(|(n, v)| n * v).sum::<f64>();
    let xc = if density.is_empty() {
        Ok((Vec::new(), 0., 0.))
    } else {
        let n = Array3::from_shape_vec(
            [density.len(), 1, 1],
            density.iter().map(|n| n.abs()).collect(),
        )
        .map_err(numerical)?;
        let x = match functional {
            Functional::LdaPz => dft_physics::xc::evaluate(&n, dv * n.len() as f64),
            Functional::LdaPw => dft_physics::spin::unpolarized_pw(&n, dv * n.len() as f64),
            Functional::Pbe => {
                return Err(ScfError::InvalidInput("reference map: PBE unsupported"));
            }
        };
        x.map(|x| {
            (
                x.potential.into_raw_vec_and_offset().0,
                x.energy.0,
                x.potential_integral.0,
            )
        })
    };
    if sum(comm, &[f64::from(xc.is_err())])?[0] != 0. {
        return Err(ScfError::InvalidInput(
            "reference map: invalid XC density on a rank",
        ));
    }
    let (vxc, ex, _) = xc.map_err(numerical)?;
    // The SCF-cycle note specifies |rho| for LDA evaluation; the signed
    // mixed density itself is never clipped, damped or renormalized.
    let vx = dv * density.iter().zip(&vxc).map(|(n, v)| n * v).sum::<f64>();
    let vh = hxc.clone();
    for (h, x) in hxc.iter_mut().zip(vxc) {
        *h += x;
    }
    let e = sum(comm, &[eh, ex, vx])?;
    Ok(Fields {
        hartree: e[0],
        xc: e[1],
        xc_integral: e[2],
        hxc,
        vh,
    })
}
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate<K: Collective>(
    state: &mut StartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    input: &[f64],
    occupations: &[f64],
    ethr_ry: f64,
    first_step: bool,
    diagnostics: bool,
    rmm: Option<dft_solvers::ref_rmm::Options>,
) -> Result<Map, ScfError> {
    let mut timings = PhaseTimings::default();
    let clock = Instant::now();
    let input_fields = fields(state, comm, cell, functional, input)?;
    let potential: Vec<_> = state
        .ionic
        .iter()
        .zip(&input_fields.hxc)
        .map(|(i, v)| i + v)
        .collect();
    let points = state.grid.dims().iter().product::<usize>();
    let v0 = sum(comm, &[potential.iter().sum::<f64>()])?[0] / points as f64;
    let diagonal: Vec<_> = state.diagonal_static.iter().map(|d| d + v0).collect();
    let mut planes = crate::grid_dist::PlaneMap::new(
        state.grid.dims(),
        comm.rank(),
        comm.size(),
        state.fields.first_axis(),
        state.operator.plane_layout().0,
    );
    let mut op_potential = Vec::new();
    planes.a_to_b(comm, &potential, &mut op_potential)?;
    state
        .operator
        .set_owned_potential(&op_potential, 0)
        .map_err(numerical)?;
    timings.potential = clock.elapsed();
    let clock = Instant::now();
    let mut options = ref_davidson::Options::new(
        occupations.len(),
        ethr_ry / 2.,
        state.electrons as usize / 2,
    );
    options.first_step = first_step;
    options.return_hx = true;
    options.profile = diagnostics;
    let solve = if let Some(mut rmm) = rmm {
        // The RMM preconditioner needs bare t_G, not Davidson's full H diagonal.
        let kinetic: Vec<_> = state
            .operator
            .local_rows()
            .iter()
            .map(|&row| {
                let g = if row == 0 {
                    0
                } else {
                    state.half.plus()[row.div_ceil(2)] as usize
                };
                state.basis.vectors()[g].norm_squared / 2.
            })
            .collect();
        rmm.bands = options.bands;
        rmm.occupied = options.occupied;
        rmm.threshold_ha = options.threshold;
        rmm.profile = options.profile;
        dft_solvers::ref_rmm::solve(comm, &mut state.operator, &state.orbitals, &kinetic, rmm)?
    } else {
        ref_davidson::solve(
            comm,
            &mut state.operator,
            &state.orbitals,
            &diagonal,
            options,
        )?
    };
    // An exhausted inner solve must never be silently labelled converged.
    // Re-entry policy belongs to the caller and is recorded as another call.
    timings.diagonalization = clock.elapsed();
    let clock = Instant::now();
    let mut op_density = vec![0.; state.operator.plane_points()];
    state
        .operator
        .band_density_sum(
            comm,
            &solve.eigenvectors,
            occupations,
            cell.volume_bohr3(),
            &mut op_density,
        )
        .map_err(numerical)?;
    let mut density = Vec::new();
    planes.b_to_a(comm, &op_density, &mut density)?;
    timings.density = clock.elapsed();
    let clock = Instant::now();
    let rows = state.operator.local_rows();
    let hx = solve.hx.as_ref().expect("requested HX");
    let local_r: Vec<_> = solve
        .eigenvalues
        .iter()
        .enumerate()
        .map(|(j, &e)| {
            (0..rows.len())
                .map(|i| (hx[(i, j)] - e * solve.eigenvectors[(i, j)]).powi(2))
                .sum::<f64>()
        })
        .collect();
    let residuals: Vec<_> = sum(comm, &local_r)?.into_iter().map(f64::sqrt).collect();
    let mut kinetic = 0.;
    for (j, &f) in occupations.iter().enumerate() {
        for (i, &row) in rows.iter().enumerate() {
            let g = if row == 0 {
                0
            } else {
                state.half.plus()[row.div_ceil(2)] as usize
            };
            kinetic +=
                f * solve.eigenvectors[(i, j)].powi(2) * state.basis.vectors()[g].norm_squared / 2.;
        }
    }
    let mut nl = vec![0.; occupations.len()];
    state
        .operator
        .nonlocal_band_energies(comm, &solve.eigenvectors, &mut nl)
        .map_err(numerical)?;
    let nonlocal: f64 = nl.iter().zip(occupations).map(|(v, f)| v * f).sum();
    let dv = cell.volume_bohr3() / points as f64;
    let integrals = sum(
        comm,
        &[
            kinetic,
            dv * density
                .iter()
                .zip(&state.ionic)
                .map(|(n, v)| n * v)
                .sum::<f64>(),
            dv * density
                .iter()
                .zip(&input_fields.hxc)
                .map(|(n, v)| n * v)
                .sum::<f64>(),
            state.ionic.iter().sum::<f64>() / points as f64,
            dv * density
                .iter()
                .zip(input)
                .zip(&input_fields.hxc)
                .map(|((o, i), v)| (o - i) * v)
                .sum::<f64>(),
            dv * density
                .iter()
                .zip(input)
                .zip(&input_fields.vh)
                .map(|((o, i), v)| (o - i) * v)
                .sum::<f64>(),
        ],
    )?;
    let output_fields = fields(state, comm, cell, functional, &density)?;
    let band_sum: f64 = solve
        .eigenvalues
        .iter()
        .zip(occupations)
        .map(|(e, f)| e * f)
        .sum();
    let ewald = state.ewald.energy.0;
    let energy = EnergyAudit {
        xc_floor: None,
        nlcc: None,
        kinetic: integrals[0],
        local: integrals[1],
        nonlocal,
        local_zero_mode: state.electrons * integrals[3],
        hartree: output_fields.hartree,
        xc: output_fields.xc,
        ewald,
        band_sum,
        double_counting: integrals[2],
        band_form: band_sum - integrals[2] + output_fields.hartree + output_fields.xc + ewald,
        harris_foulkes: band_sum - input_fields.hartree + input_fields.xc
            - input_fields.xc_integral
            + ewald,
        functional_difference_bound: (output_fields.hartree - input_fields.hartree - integrals[5])
            .abs()
            + (output_fields.xc - input_fields.xc - (integrals[4] - integrals[5])).abs()
            + residuals
                .iter()
                .zip(occupations)
                .map(|(r, f)| r * f)
                .sum::<f64>(),
    };
    timings.energy = clock.elapsed();
    if !energy.total_hartree().is_finite() || residuals.iter().any(|r| !r.is_finite()) {
        return Err(ScfError::NonFinite);
    }
    Ok(Map {
        solve,
        density,
        residuals,
        energy,
        timings,
    })
}
