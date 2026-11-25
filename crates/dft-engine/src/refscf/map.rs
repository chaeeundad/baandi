//! Owned-row map for the reference loop. Explicit kernels and call-owned
//! options only: no ScfProblem, session or legacy solver dispatcher.
use super::numerical;
use crate::{Functional, PhaseTimings, ScfError, energy::EnergyAudit, refstart::StartState};
use dft_geometry::Cell;
use dft_mpi::Collective;
use dft_physics::data::{Array3, Complex64 as C};
use dft_solvers::ref_davidson;
use std::{f64::consts::PI, time::Instant};

#[derive(Clone)]
pub(super) struct Fields {
    pub hartree: f64,
    pub xc: f64,
    pub xc_integral: f64,
    pub hxc: Vec<f64>,
    pub vh: Vec<f64>,
}
pub(super) struct Map {
    pub solve: ref_davidson::Solution,
    pub warm_start: bool,
    pub density: Vec<f64>,
    pub residuals: Vec<f64>,
    pub energy: EnergyAudit,
    pub timings: PhaseTimings,
    input_fields: Fields,
}
fn sum<K: Collective>(comm: &K, values: &[f64]) -> Result<Vec<f64>, ScfError> {
    let mut out = vec![0.; values.len()];
    comm.tree_sum_f64(values, &mut out).map_err(numerical)?;
    Ok(out)
}
/// Full-box Hartree and LDA / distributed PBE, in the same convention as refstart.
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
    // Reuse rho(G) before the Poisson multiplier overwrites it. GGA owns only
    // this rank's slabs/sticks and follows the serial signed-density policy.
    let pbe = if functional == Functional::Pbe {
        Some(dft_physics::gga::distributed::evaluate(
            &mut state.fields,
            comm,
            cell,
            density,
            &coeff,
        )?)
    } else {
        None
    };
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
    let xc = if let Some(x) = pbe {
        Ok((x.potential, x.energy, 0.))
    } else if density.is_empty() {
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
                return Err(ScfError::InvalidInput(
                    "reference map: PBE must use the distributed GGA evaluator",
                ));
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
    cache: &mut Option<(Vec<f64>, Fields)>,
    lean: bool,
    rmm: Option<dft_solvers::ref_rmm::Options>,
    kept: &mut super::kept_hx::KeptHx,
) -> Result<Map, ScfError> {
    let mut timings = PhaseTimings::default();
    let clock = Instant::now();
    let input_fields = cached_fields(state, comm, cell, functional, input, cache)?;
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
    let epoch = state
        .operator
        .owned_epoch()
        .unwrap_or(0)
        .checked_add(1)
        .ok_or(ScfError::InvalidInput("reference potential epoch overflow"))?;
    state
        .operator
        .set_owned_potential(&op_potential, epoch)
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
    if rmm.is_some() {
        kept.invalidate();
    }
    let initial_hx = if rmm.is_none() {
        kept.prepare(comm, &mut state.operator, &state.orbitals, &op_potential)?
    } else {
        None
    };
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
        ref_davidson::solve_with_initial_hx(
            comm,
            &mut state.operator,
            &state.orbitals,
            &diagonal,
            options,
            initial_hx.as_ref(),
        )?
    };
    // An exhausted inner solve must never be silently labelled converged.
    // Re-entry policy belongs to the caller and is recorded as another call.
    timings.diagonalization = clock.elapsed();
    let clock = Instant::now();
    let mut op_density = vec![0.; state.operator.plane_points()];
    kept.density(
        comm,
        &mut state.operator,
        &solve.eigenvectors,
        solve.hx.as_ref().expect("requested HX"),
        &op_potential,
        occupations,
        cell.volume_bohr3(),
        &mut op_density,
    )?;
    if rmm.is_some() {
        kept.invalidate();
    }
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
    let energy = audit(
        state,
        comm,
        cell,
        functional,
        input,
        occupations,
        &solve,
        &density,
        &residuals,
        &input_fields,
        cache,
        !lean,
    )?;
    timings.energy = clock.elapsed();
    if (!lean && !energy.total_hartree().is_finite())
        || !energy.band_sum.is_finite()
        || residuals.iter().any(|r| !r.is_finite())
    {
        return Err(ScfError::NonFinite);
    }
    Ok(Map {
        warm_start: initial_hx.is_some(),
        solve,
        density,
        residuals,
        energy,
        timings,
        input_fields,
    })
}

/// One-entry cache. Density equality is agreed globally before skipping any
/// collective: local equality alone could deadlock ranks with empty planes.
pub(super) fn cached_fields<K: Collective>(
    state: &mut StartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    density: &[f64],
    cache: &mut Option<(Vec<f64>, Fields)>,
) -> Result<Fields, ScfError> {
    let hit = cache.as_ref().is_some_and(|(n, _)| {
        n.len() == density.len()
            && n.iter()
                .zip(density)
                .all(|(a, b)| a.to_bits() == b.to_bits())
    });
    if sum(comm, &[f64::from(!hit)])?[0] == 0. {
        return Ok(cache.as_ref().expect("global cache hit").1.clone());
    }
    let value = fields(state, comm, cell, functional, density)?;
    *cache = Some((density.to_vec(), value.clone()));
    Ok(value)
}

#[allow(clippy::too_many_arguments)]
fn audit<K: Collective>(
    state: &mut StartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    input: &[f64],
    occupations: &[f64],
    solve: &ref_davidson::Solution,
    density: &[f64],
    residuals: &[f64],
    input_fields: &Fields,
    cache: &mut Option<(Vec<f64>, Fields)>,
    full: bool,
) -> Result<EnergyAudit, ScfError> {
    let rows = state.operator.local_rows();
    let points = state.grid.dims().iter().product::<usize>();
    if !full {
        // NaN is an unavailable audit slot, never a convergence input. The
        // existing JSON float serializer writes it as null. Delivery fills all
        // slots before exposing a ScfResult.
        let dv = cell.volume_bohr3() / points as f64;
        let dc = sum(
            comm,
            &[dv * density
                .iter()
                .zip(&input_fields.hxc)
                .map(|(n, v)| n * v)
                .sum::<f64>()],
        )?[0];
        let band_sum = solve
            .eigenvalues
            .iter()
            .zip(occupations)
            .map(|(e, f)| e * f)
            .sum();
        let ewald = state.ewald.energy.0;
        return Ok(EnergyAudit {
            xc_floor: None,
            nlcc: None,
            kinetic: f64::NAN,
            local: f64::NAN,
            nonlocal: f64::NAN,
            local_zero_mode: f64::NAN,
            hartree: f64::NAN,
            xc: f64::NAN,
            ewald,
            band_sum,
            double_counting: dc,
            band_form: f64::NAN,
            harris_foulkes: band_sum - input_fields.hartree + input_fields.xc
                - input_fields.xc_integral
                + ewald,
            functional_difference_bound: f64::NAN,
        });
    }
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
    let nonlocal = nl.iter().zip(occupations).map(|(v, f)| v * f).sum();

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
    let output_fields = cached_fields(state, comm, cell, functional, density, cache)?;
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
    Ok(energy)
}

#[allow(clippy::too_many_arguments)]
pub(super) fn complete_audit<K: Collective>(
    step: &mut Map,
    state: &mut StartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    input: &[f64],
    occupations: &[f64],
    cache: &mut Option<(Vec<f64>, Fields)>,
) -> Result<(), ScfError> {
    if step.energy.kinetic.is_nan() {
        step.energy = audit(
            state,
            comm,
            cell,
            functional,
            input,
            occupations,
            &step.solve,
            &step.density,
            &step.residuals,
            &step.input_fields,
            cache,
            true,
        )?;
        if !step.energy.total_hartree().is_finite() {
            return Err(ScfError::NonFinite);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn cached_fields_equal_fresh_bits_and_agree_on_local_misses() {
        dft_mpi::simulate(3, |comm| {
            let pp =
                crate::scf::pool_fixture::potential_with_atomic_start(2., 1.8, true, 161, "SLA PZ");
            let cell =
                Cell::from_user_rows_bohr([[7., 0., 0.], [0.2, 7., 0.], [0., 0.1, 7.]]).unwrap();
            let atoms = [dft_geometry::Atom::new("Xx", [0.13, 0.21, 0.32])];
            let potentials = [crate::SpeciesPotential {
                label: "Xx",
                potential: &pp,
            }];
            let options = crate::ScfOptions {
                gamma_real: crate::GammaReal::On,
                cutoff_hartree: 1.,
                bands: Some(3),
                ..Default::default()
            };
            let mut state = crate::refstart::start(
                crate::refstart::StartInput {
                    cell: &cell,
                    atoms: &atoms,
                    potentials: &potentials,
                    options: &options,
                },
                comm,
            )
            .unwrap();
            let mut density = state.density.clone();
            let mut cache = None;
            for change in [false, false, true, false] {
                if change && comm.rank() == 0 {
                    density[0] += 1e-8;
                }
                let cached = cached_fields(
                    &mut state,
                    comm,
                    &cell,
                    Functional::LdaPz,
                    &density,
                    &mut cache,
                )
                .unwrap();
                let fresh = fields(&mut state, comm, &cell, Functional::LdaPz, &density).unwrap();
                assert_eq!(cached.hartree.to_bits(), fresh.hartree.to_bits());
                assert_eq!(cached.xc.to_bits(), fresh.xc.to_bits());
                assert_eq!(cached.xc_integral.to_bits(), fresh.xc_integral.to_bits());
                assert_eq!(cached.hxc, fresh.hxc);
                assert_eq!(cached.vh, fresh.vh);
            }
        });
    }
}
