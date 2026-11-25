//! Spin/k/band map: all eigenproblems precede the global occupation solve.
use super::{electronic::Policy, k_kept::KeptHx, numerical, record::SolveRecord};
use crate::{Functional, ScfError, energy::EnergyAudit, refstart::KStartState};
use dft_geometry::Cell;
use dft_mpi::Collective;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::{
    data::Array3,
    occupations::{Evaluation, FermiDirac, Layout},
    spin::{SpinFunctional, SpinXc, SpinXcResponse},
};
use dft_solvers::{ref_davidson, ref_davidson_complex};

/// The local LDA adapter owns spectral scratch even though LDA needs no
/// derivatives. Construct it once, never once per SCF map or verification.
pub(super) struct Workspace {
    spin: Option<SpinFields>,
}
struct SpinFields {
    xc: SpinXc,
    response: SpinXcResponse,
    channels: [Array3<f64>; 2],
    core: Option<Array3<f64>>,
    kernel_difference: Vec<f64>,
    kernel_span: Vec<f64>,
}
impl Workspace {
    pub fn new<K: Collective>(
        state: &KStartState,
        comm: &K,
        cell: &Cell,
        spin: bool,
    ) -> Result<Self, ScfError> {
        let local = (|| -> Result<_, ScfError> {
            let len = state.density.len();
            if !spin || len == 0 {
                return Ok(None);
            }
            let dims = [len, 1, 1];
            Ok(Some(SpinFields {
                xc: SpinXc::new(SpinFunctional::LdaPw, cell, dims)?,
                response: SpinXcResponse::new(dims),
                kernel_difference: vec![0.; len],
                kernel_span: vec![0.; len],
                channels: std::array::from_fn(|_| Array3::zeros(dims)),
                core: state
                    .core_density
                    .as_ref()
                    .map(|c| Array3::from_shape_vec(dims, c.clone()).map_err(numerical))
                    .transpose()?,
            }))
        })();
        if sum(comm, &[f64::from(local.is_err())])?[0] != 0. {
            return Err(ScfError::InvalidInput(
                "reference spin XC setup failed on a rank",
            ));
        }
        Ok(Self { spin: local? })
    }
    /// Positive residual measure 1/2 integral |d²e_xc/dm²| d_m², Ha.
    /// Finite differences of the existing polarized LDA-PW potentials hold n
    /// fixed. h/n=cbrt(epsilon) balances second-order truncation and roundoff.
    /// At full polarization use the feasible one-sided span, never negative rho.
    pub fn spin_stiffness<K: Collective>(
        &mut self,
        comm: &K,
        input: &[Vec<f64>],
        output: &[Vec<f64>],
        dv: f64,
    ) -> Result<f64, ScfError> {
        let local = (|| -> Result<f64, ScfError> {
            let Some(spin) = self.spin.as_mut() else {
                return Ok(0.);
            };
            for sign in [-1., 1.] {
                for (i, (&up, &down)) in input[0].iter().zip(&input[1]).enumerate() {
                    let core = spin.core.as_ref().map_or(0., |c| c.as_slice().unwrap()[i]);
                    let u = up + 0.5 * core;
                    let d = down + 0.5 * core;
                    let n = u + d;
                    let h = n * f64::EPSILON.cbrt();
                    let shift = (sign * h).max(-2. * u).min(2. * d);
                    spin.channels[0].as_slice_mut().unwrap()[i] = u + 0.5 * shift;
                    spin.channels[1].as_slice_mut().unwrap()[i] = d - 0.5 * shift;
                    if sign < 0. {
                        spin.kernel_span[i] = -shift;
                    } else {
                        spin.kernel_span[i] += shift;
                    }
                }
                // The scratch channels already include the immutable core.
                spin.xc.evaluate(&spin.channels, &mut spin.response)?;
                for i in 0..input[0].len() {
                    let vm = 0.5
                        * (spin.response.potentials[0].as_slice().unwrap()[i]
                            - spin.response.potentials[1].as_slice().unwrap()[i]);
                    if sign < 0. {
                        spin.kernel_difference[i] = -vm;
                    } else {
                        spin.kernel_difference[i] += vm;
                    }
                }
            }
            let mut energy = 0.;
            for i in 0..input[0].len() {
                let span = spin.kernel_span[i];
                if span > 0. {
                    let kernel = spin.kernel_difference[i] / span;
                    let dm = (output[0][i] - output[1][i]) - (input[0][i] - input[1][i]);
                    energy += 0.5 * kernel.abs() * dm * dm;
                }
            }
            Ok(dv * energy)
        })();
        if sum(comm, &[f64::from(local.is_err())])?[0] != 0. {
            return Err(numerical("spin stiffness evaluator failed on a rank").into());
        }
        let energy = sum(comm, &[local?])?[0];
        if !energy.is_finite() {
            return Err(ScfError::NonFinite);
        }
        Ok(energy)
    }
}

pub(super) struct Block {
    pub eigenvalues: Vec<f64>,
    pub occupations: Vec<f64>,
    pub residuals: Vec<f64>,
    pub orbitals: ComplexMatrix,
    pub calls: Vec<SolveRecord>,
}
pub(super) struct Map {
    pub density: Vec<Vec<f64>>,
    pub blocks: Vec<Vec<Block>>,
    pub evaluation: Option<Evaluation>,
    pub energy: EnergyAudit,
    pub guard: f64,
    pub charges: Vec<f64>,
    pub moment: f64,
    pub absolute_moment: f64,
}
impl Map {
    pub fn entropy(&self) -> f64 {
        self.evaluation.map_or(0., |e| e.entropy_term_ha)
    }
    pub fn free_energy(&self) -> f64 {
        self.energy.total_hartree() - self.entropy()
    }
}
pub(super) fn sum<K: Collective>(comm: &K, v: &[f64]) -> Result<Vec<f64>, ScfError> {
    let mut out = vec![0.; v.len()];
    comm.tree_sum_f64(v, &mut out).map_err(numerical)?;
    Ok(out)
}
struct Fields {
    hxc: Vec<Vec<f64>>,
    hartree: f64,
    xc: f64,
    xc_integral: f64,
    vh: Vec<f64>,
    nlcc: Option<crate::energy::CoreXcLedger>,
}
fn fields<K: Collective>(
    workspace: &mut Workspace,
    state: &mut KStartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    density: &[Vec<f64>],
) -> Result<Fields, ScfError> {
    let total: Vec<_> = (0..density[0].len())
        .map(|i| density.iter().map(|n| n[i]).sum::<f64>())
        .collect();
    let scalar = super::map::fields(&mut state.common, comm, cell, functional, &total)?;
    if density.len() == 1 {
        return Ok(Fields {
            hxc: vec![scalar.hxc],
            vh: scalar.vh,
            hartree: scalar.hartree,
            xc: scalar.xc,
            xc_integral: scalar.xc_integral,
            nlcc: scalar.nlcc,
        });
    }
    // LDA is point-local. Evaluate the existing polarized adapter on owned
    // points, then convert its local quadrature to the physical full-grid dv.
    // No density payload is gathered. Empty slab owners join every reduction.
    let len = total.len();
    let points = state.grid.dims().iter().product::<usize>();
    let local = (|| -> Result<_, ScfError> {
        if len == 0 {
            return Ok((vec![Vec::new(), Vec::new()], [0.; 4]));
        }
        let spin = workspace
            .spin
            .as_mut()
            .expect("prepared nonempty spin slab");
        for (buffer, input) in spin.channels.iter_mut().zip(density) {
            buffer.as_slice_mut().unwrap().copy_from_slice(input);
        }
        spin.xc
            .evaluate_with_core(&spin.channels, spin.core.as_ref(), &mut spin.response)?;
        let out = &spin.response;
        let scale = len as f64 / points as f64;
        let values = [
            out.energy_ha * scale,
            out.valence_potential_integral_ha * scale,
            out.core_potential_integral_ha * scale,
            state
                .core_density
                .as_ref()
                .map_or(0., |c| c.iter().sum::<f64>())
                * cell.volume_bohr3()
                / points as f64,
        ];
        Ok((
            out.potentials
                .iter()
                .map(|v| v.as_slice().unwrap().to_vec())
                .collect::<Vec<_>>(),
            values,
        ))
    })();
    if sum(comm, &[f64::from(local.is_err())])?[0] != 0. {
        return Err(numerical(format!("reference spin XC: {:?}", local.err())).into());
    }
    let (mut potentials, values) = local?;
    let values = sum(comm, &values)?;
    for potential in &mut potentials {
        for (v, h) in potential.iter_mut().zip(&scalar.vh) {
            *v += h;
        }
    }
    Ok(Fields {
        hxc: potentials,
        vh: scalar.vh,
        hartree: scalar.hartree,
        xc: values[0],
        xc_integral: values[1],
        nlcc: state
            .core_density
            .as_ref()
            .map(|_| crate::energy::CoreXcLedger {
                core_electrons: values[3],
                core_potential_integral_ha: values[2],
                valence_potential_integral_ha: values[1],
            }),
    })
}
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate<K: Collective>(
    workspace: &mut Workspace,
    states: &mut [KStartState],
    kept: &mut [KeptHx],
    comm: &K,
    cell: &Cell,
    functional: Functional,
    input: &[Vec<f64>],
    policy: &Policy,
    ethr: f64,
    first: bool,
    diagnostics: bool,
) -> Result<Map, ScfError> {
    let fields = fields(workspace, &mut states[0], comm, cell, functional, input)?;
    let bands = states[0].orbitals.cols();
    let nk = states[0].blocks.len();
    let ns = states.len();
    let electrons = states[0].electrons;
    let dims = states[0].grid.dims();
    let points = dims.iter().product::<usize>();
    let dv = cell.volume_bohr3() / points as f64;
    let mut blocks = Vec::with_capacity(ns);
    for (s, state) in states.iter_mut().enumerate() {
        let potential: Vec<_> = state
            .ionic
            .iter()
            .zip(&fields.hxc[s])
            .map(|(a, b)| a + b)
            .collect();
        let v0 = sum(comm, &[potential.iter().sum::<f64>()])?[0] / points as f64;
        let mut channel = Vec::with_capacity(nk);
        for k in 0..nk {
            let mut planes = crate::grid_dist::PlaneMap::new(
                dims,
                comm.rank(),
                comm.size(),
                state.fields.first_axis(),
                state.blocks[k].operator.first_axis(),
            );
            let mut op_potential = Vec::new();
            planes.a_to_b(comm, &potential, &mut op_potential)?;
            let op = &mut state.blocks[k].operator;
            let epoch = op
                .owned_epoch()
                .unwrap_or(0)
                .checked_add(1)
                .ok_or(ScfError::InvalidInput("spin potential epoch overflow"))?;
            op.set_owned_potential(&op_potential, epoch)
                .map_err(numerical)?;
            let diagonal: Vec<_> = op.diagonal_static().iter().map(|d| d + v0).collect();
            let hx = kept[s].prepare(k, comm, op, &state.orbitals.0[k], &op_potential)?;
            // FD: all retained bands are mandatory, even underflowed tails.
            // Fixed sector: first empty is mandatory as in the legacy policy.
            let strict = policy.strict_bands(s, bands);
            let mut options = ref_davidson::Options::new(bands, ethr / 2., strict);
            options.first_step = first;
            options.return_hx = true;
            options.profile = diagnostics;
            let mut solution = ref_davidson_complex::solve_with_initial_hx(
                comm,
                op,
                &state.orbitals.0[k],
                &diagonal,
                options,
                hx.as_ref(),
            )
            .map_err(numerical)?;
            let mut calls = Vec::new();
            for retry in 0..6 {
                calls.push(super::k_driver::solve_record(
                    &mut solution,
                    k,
                    ethr,
                    first && retry == 0,
                    diagnostics,
                ));
                // The fixed-occupation reference driver permits a few pending
                // buffer roots. That allowance cannot include mandatory FD
                // roots (or the fixed sector's gap edge): any of them can
                // change the number solve and the resulting spin density.
                if solution.converged[..strict].iter().all(|&v| v)
                    && solution.converged.iter().filter(|&&v| !v).count() <= 5.max(bands / 4)
                {
                    break;
                }
                if retry == 5 {
                    return Err(ScfError::InvalidInput(
                        "electronic reference Davidson exhausted six calls",
                    ));
                }
                options.first_step = false;
                solution = ref_davidson_complex::solve(
                    comm,
                    op,
                    &solution.eigenvectors,
                    &diagonal,
                    options,
                )
                .map_err(numerical)?;
            }
            let hx = solution.hx.as_ref().expect("requested HX");
            let local: Vec<_> = (0..bands)
                .map(|j| {
                    (0..op.local_rows().len())
                        .map(|i| {
                            (hx[(i, j)] - solution.eigenvalues[j] * solution.eigenvectors[(i, j)])
                                .norm_sqr()
                        })
                        .sum()
                })
                .collect();
            let residuals = sum(comm, &local)?.into_iter().map(f64::sqrt).collect();
            kept[s].retain(k, state, &solution.eigenvectors, hx, &op_potential);
            channel.push(Block {
                eigenvalues: solution.eigenvalues,
                occupations: vec![0.; bands],
                residuals,
                orbitals: solution.eigenvectors,
                calls,
            });
        }
        blocks.push(channel);
    }
    let weights: Vec<_> = states[0].blocks.iter().map(|b| b.point.weight).collect();
    let evaluation = if let Some((theta, constraint)) = policy.finite() {
        let energies: Vec<_> = blocks
            .iter()
            .flatten()
            .flat_map(|b| b.eigenvalues.iter().copied())
            .collect();
        let mut occupations = vec![0.; energies.len()];
        let evaluation = FermiDirac::new(theta)
            .evaluate(
                Layout {
                    capacities: if ns == 1 { &[2.] } else { &[1., 1.] },
                    weights: &weights,
                    bands,
                },
                &energies,
                electrons,
                constraint,
                &mut occupations,
            )
            .map_err(numerical)?;
        if evaluation.edge_max_fraction > 1e-8
            || evaluation.edge_electrons > 1e-10 * electrons.max(1.)
        {
            return Err(ScfError::InvalidInput(
                "occupations.bands: insufficient empty-state coverage at explicit cap; increase cap and repeat the convergence study",
            ));
        }
        for (block, f) in blocks
            .iter_mut()
            .flatten()
            .zip(occupations.chunks_exact(bands))
        {
            block.occupations.copy_from_slice(f);
        }
        Some(evaluation)
    } else {
        for (s, channel) in blocks.iter_mut().enumerate() {
            for block in channel {
                block.occupations[..policy.fixed_count(s)].fill(1.);
            }
        }
        None
    };
    let mut density = vec![vec![0.; input[0].len()]; ns];
    let (mut kinetic, mut nonlocal, mut band_sum, mut guard, mut residual_bound) =
        (0., 0., 0., 0., 0.);
    for (s, channel) in blocks.iter().enumerate() {
        let state = &mut states[s];
        for (k, block) in channel.iter().enumerate() {
            let op = &mut state.blocks[k].operator;
            let f: Vec<_> = block.occupations.iter().map(|v| v * weights[k]).collect();
            for (j, &weight) in f.iter().enumerate() {
                band_sum += weight * block.eigenvalues[j];
                guard += weight * block.residuals[j].powi(2);
                residual_bound += weight * block.residuals[j];
                kinetic += weight
                    * op.kinetic_diagonal()
                        .iter()
                        .enumerate()
                        .map(|(i, t)| t * block.orbitals[(i, j)].norm_sqr())
                        .sum::<f64>();
            }
            let mut nl = vec![0.; bands];
            op.nonlocal_band_energies(comm, &block.orbitals, &mut nl)
                .map_err(numerical)?;
            nonlocal += nl.iter().zip(&f).map(|(a, b)| a * b).sum::<f64>();
            let mut plane_density = vec![0.; op.plane_points()];
            op.band_density_sum(
                comm,
                &block.orbitals,
                &f,
                cell.volume_bohr3(),
                &mut plane_density,
            )
            .map_err(numerical)?;
            let mut planes = crate::grid_dist::PlaneMap::new(
                dims,
                comm.rank(),
                comm.size(),
                state.fields.first_axis(),
                state.blocks[k].operator.first_axis(),
            );
            let mut local = Vec::new();
            planes.b_to_a(comm, &plane_density, &mut local)?;
            for (a, b) in density[s].iter_mut().zip(local) {
                *a += b;
            }
        }
    }
    let total: Vec<_> = (0..density[0].len())
        .map(|i| density.iter().map(|n| n[i]).sum::<f64>())
        .collect();
    let integral = sum(
        comm,
        &[
            kinetic,
            dv * total
                .iter()
                .zip(&states[0].ionic)
                .map(|(n, v)| n * v)
                .sum::<f64>(),
            dv * density
                .iter()
                .zip(&fields.hxc)
                .map(|(n, v)| n.iter().zip(v).map(|(n, v)| n * v).sum::<f64>())
                .sum::<f64>(),
            states[0].ionic.iter().sum::<f64>() / points as f64,
            dv * density[0].iter().sum::<f64>(),
            if ns == 2 {
                dv * density[1].iter().sum::<f64>()
            } else {
                0.
            },
            if ns == 2 {
                dv * density[0]
                    .iter()
                    .zip(&density[1])
                    .map(|(u, d)| (u - d).abs())
                    .sum::<f64>()
            } else {
                0.
            },
        ],
    )?;
    let remainder = sum(
        comm,
        &[
            dv * density
                .iter()
                .zip(input)
                .zip(&fields.hxc)
                .map(|((a, b), v)| {
                    a.iter()
                        .zip(b)
                        .zip(v)
                        .map(|((a, b), v)| (a - b) * v)
                        .sum::<f64>()
                })
                .sum::<f64>(),
            dv * total
                .iter()
                .enumerate()
                .zip(&fields.vh)
                .map(|((i, n), v)| (n - input.iter().map(|c| c[i]).sum::<f64>()) * v)
                .sum::<f64>(),
        ],
    )?;
    let output = fields_output(workspace, &mut states[0], comm, cell, functional, &density)?;
    let ewald = states[0].ewald.energy.0;
    let energy = EnergyAudit {
        xc_floor: None,
        nlcc: output.nlcc,
        kinetic: integral[0],
        local: integral[1],
        nonlocal,
        local_zero_mode: electrons * integral[3],
        hartree: output.hartree,
        xc: output.xc,
        ewald,
        band_sum,
        double_counting: integral[2],
        band_form: band_sum - integral[2] + output.hartree + output.xc + ewald,
        harris_foulkes: band_sum - fields.hartree + fields.xc - fields.xc_integral + ewald,
        functional_difference_bound: (output.hartree - fields.hartree - remainder[1]).abs()
            + (output.xc - fields.xc - (remainder[0] - remainder[1])).abs()
            + residual_bound,
    };
    if !energy.total_hartree().is_finite() {
        return Err(ScfError::NonFinite);
    }
    Ok(Map {
        density,
        blocks,
        evaluation,
        energy,
        guard,
        charges: integral[4..4 + ns].to_vec(),
        moment: if ns == 2 {
            integral[4] - integral[5]
        } else {
            0.
        },
        absolute_moment: integral[6],
    })
}
// Separate name avoids shadowing the input Fields value in evaluate.
fn fields_output<K: Collective>(
    workspace: &mut Workspace,
    state: &mut KStartState,
    comm: &K,
    cell: &Cell,
    functional: Functional,
    density: &[Vec<f64>],
) -> Result<Fields, ScfError> {
    fields(workspace, state, comm, cell, functional, density)
}

#[cfg(test)]
mod g4b_kernel_tests {
    use super::*;
    use crate::refscf::metric::Serial;

    #[test]
    fn spin_curvature_matches_energy_second_difference_with_core() {
        let cell = Cell::from_user_rows_bohr([[5., 0., 0.], [0.2, 4., 0.], [0., 0.1, 3.]]).unwrap();
        let dims = [2; 3];
        for core in [0., 0.02] {
            let mut work = Workspace {
                spin: Some(SpinFields {
                    xc: SpinXc::new(SpinFunctional::LdaPw, &cell, dims).unwrap(),
                    response: SpinXcResponse::new(dims),
                    channels: std::array::from_fn(|_| Array3::zeros(dims)),
                    core: Some(Array3::from_elem(dims, core)),
                    kernel_difference: vec![0.; 8],
                    kernel_span: vec![0.; 8],
                }),
            };
            let input = vec![vec![0.04; 8], vec![0.02; 8]];
            let dm = 1e-5;
            let output = vec![vec![0.04 + dm / 2.; 8], vec![0.02 - dm / 2.; 8]];
            let measured = work
                .spin_stiffness(&Serial, &input, &output, cell.volume_bohr3() / 8.)
                .unwrap();
            let h = (0.06 + core) * 1e-4;
            let spin = work.spin.as_mut().unwrap();
            let mut energy = |delta: f64| {
                spin.channels[0].fill(0.04 + delta / 2.);
                spin.channels[1].fill(0.02 - delta / 2.);
                spin.xc
                    .evaluate_with_core(&spin.channels, spin.core.as_ref(), &mut spin.response)
                    .unwrap();
                spin.response.energy_ha
            };
            let curvature = (energy(h) + energy(-h) - 2. * energy(0.)) / (h * h);
            let expected = 0.5 * curvature.abs() * dm * dm;
            assert!(
                (measured / expected - 1.).abs() < 1e-5,
                "{measured:e} Ha vs {expected:e} Ha"
            );
            assert_eq!(
                work.spin_stiffness(&Serial, &input, &input, cell.volume_bohr3() / 8.)
                    .unwrap(),
                0.
            );
            // Empty and exactly polarized points use a feasible span, never a
            // negative channel; no arbitrary density floor is introduced.
            for u in [0., 0.06] {
                let edge = vec![vec![u; 8], vec![0.; 8]];
                assert_eq!(
                    work.spin_stiffness(&Serial, &edge, &edge, cell.volume_bohr3() / 8.)
                        .unwrap(),
                    0.
                );
            }
        }
    }
}
