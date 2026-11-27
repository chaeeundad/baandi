//! Spin/k/band map: all eigenproblems precede the global occupation solve.
use super::{electronic::Policy, k_kept::KeptHx, numerical, record::SolveRecord};
use crate::{Functional, energy::EnergyAudit, model::ScfError, refstart::KStartState};
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
/// derivatives. It also supplies G4B's local curvature proxy for PBE. Both
/// it and the separate distributed PBE scratch are constructed once per run.
pub(super) struct Workspace {
    pub symmetry: Option<super::symmetry::Workspace>,
    spin: Option<SpinFields>,
    pbe: Option<dft_physics::gga::distributed::SpinPbe>,
}
struct SpinFields {
    xc: SpinXc,
    response: SpinXcResponse,
    channels: [Array3<f64>; 2],
    core: Option<Array3<f64>>,
    kernel_difference: Vec<f64>,
    kernel_span: Vec<f64>,
    kernel_indices: Vec<usize>,
}
impl Workspace {
    pub fn new<K: Collective>(
        state: &KStartState,
        comm: &K,
        cell: &Cell,
        spin: bool,
        functional: Functional,
    ) -> Result<Self, ScfError> {
        if spin && state.pw92_prefactor == crate::Pw92Prefactor::Original {
            return Err(ScfError::InvalidInput(
                "original PW92 prefactor requires unpolarized PBE",
            ));
        }
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
                kernel_indices: state.fields.plane_offsets(),
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
        if !spin || functional != Functional::Pbe {
            return Ok(Self {
                spin: local?,
                pbe: None,
                symmetry: None,
            });
        }
        let pbe = dft_physics::gga::distributed::SpinPbe::new::<K, _>(&state.fields, cell);
        if sum(comm, &[f64::from(pbe.is_err())])?[0] != 0. {
            return Err(ScfError::InvalidInput(
                "reference spin PBE setup failed on a rank",
            ));
        }
        Ok(Self {
            spin: local?,
            pbe: Some(pbe?),
            symmetry: None,
        })
    }
    /// Positive residual measure 1/2 integral |d²e_xc/dm²| d_m², Ha.
    /// Finite differences of the existing polarized LDA-PW potentials hold n
    /// fixed. h/n=cbrt(epsilon) balances second-order truncation and roundoff.
    /// Freeze the ADR 0046 mask of the already-audited combined input. A
    /// screened channel stays zero in both probes (zero response); fully
    /// screened points have zero span. Restrict active probes to half the
    /// distance to the threshold when needed, retaining a one-sided span.
    /// Off retains the historical feasible polarization-endpoint span.
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
            let screening = dft_physics::xc::screen_negative();
            // Main XC has audited this input with the GLOBAL charge budget.
            // Probes are derivative scratch, not new densities to audit using
            // the local adapter's artificial [len,1,1] grid and volume.
            let _scope =
                dft_physics::xc::screening_scope(dft_physics::xc::screening_stage(), true, None);
            for sign in [-1., 1.] {
                for (i, (&up, &down)) in input[0].iter().zip(&input[1]).enumerate() {
                    let core = spin.core.as_ref().map_or(0., |c| c.as_slice().unwrap()[i]);
                    let mut u = up + 0.5 * core;
                    let mut d = down + 0.5 * core;
                    if !u.is_finite() || !d.is_finite() {
                        return Err(numerical("nonfinite combined input density").into());
                    }
                    if !screening && (u < 0. || d < 0.) {
                        return Err(
                            numerical("negative combined input density (screening off)").into()
                        );
                    }
                    if screening {
                        u = dft_physics::xc::screened_argument(u);
                        d = dft_physics::xc::screened_argument(d);
                    }
                    let n = u + d;
                    let h = n * f64::EPSILON.cbrt();
                    let shift = if screening {
                        let t = dft_physics::xc::DENSITY_THRESHOLD;
                        let lo = if u > 0. { -(u - t) } else { -h };
                        let hi = if d > 0. { d - t } else { h };
                        let shift = (sign * h).max(lo).min(hi);
                        // Adjacent representable values can round a midpoint
                        // onto t. Use the unperturbed endpoint in that direction.
                        if (u > 0. && u + 0.5 * shift <= t) || (d > 0. && d - 0.5 * shift <= t) {
                            0.
                        } else {
                            shift
                        }
                    } else {
                        (sign * h).max(-2. * u).min(2. * d)
                    };
                    spin.channels[0].as_slice_mut().unwrap()[i] = if screening && u == 0. {
                        0.
                    } else {
                        u + 0.5 * shift
                    };
                    spin.channels[1].as_slice_mut().unwrap()[i] = if screening && d == 0. {
                        0.
                    } else {
                        d - 0.5 * shift
                    };
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
            if !energy.is_finite() || !(dv * energy).is_finite() {
                return Err(numerical("nonfinite spin stiffness contribution").into());
            }
            Ok(dv * energy)
        })();
        if sum(comm, &[f64::from(local.is_err())])?[0] != 0. {
            // Error-only communication: every rank receives the first failing
            // rank's actual reason and GLOBAL C-order indices, even empty owners.
            for rank in 0..comm.size() {
                let mut failed = [f64::from(comm.rank() == rank && local.is_err())];
                comm.bcast_f64(&mut failed, rank).map_err(numerical)?;
                if failed[0] == 0. {
                    continue;
                }
                let message = if comm.rank() == rank {
                    let mut detail = String::new();
                    if let Some(spin) = &self.spin {
                        use std::fmt::Write;
                        for (s, name) in ["up", "down"].into_iter().enumerate() {
                            let mut minimum = (f64::INFINITY, usize::MAX);
                            for (i, &n) in input[s].iter().enumerate() {
                                let n = n + spin
                                    .core
                                    .as_ref()
                                    .map_or(0., |c| 0.5 * c.as_slice().unwrap()[i]);
                                let index = spin.kernel_indices[i];
                                if !n.is_finite()
                                    || n < minimum.0
                                    || (n == minimum.0 && index < minimum.1)
                                {
                                    minimum = (n, index);
                                }
                            }
                            write!(
                                detail,
                                "; channel {name} minimum {:.17e} Bohr^-3 at grid point {}",
                                minimum.0, minimum.1
                            )
                            .unwrap();
                        }
                    }
                    format!(
                        "{}: spin stiffness rank {rank}{detail}; rank-local reason: {}",
                        dft_physics::xc::screening_stage(),
                        local.as_ref().unwrap_err()
                    )
                } else {
                    String::new()
                };
                let mut len = [message.len() as f64];
                comm.bcast_f64(&mut len, rank).map_err(numerical)?;
                let mut bytes = vec![0.; len[0] as usize];
                for (b, v) in bytes.iter_mut().zip(message.bytes()) {
                    *b = f64::from(v);
                }
                comm.bcast_f64(&mut bytes, rank).map_err(numerical)?;
                let message = String::from_utf8(bytes.into_iter().map(|b| b as u8).collect())
                    .map_err(numerical)?;
                return Err(numerical(message).into());
            }
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
    /// FD buffer bands: the complete solved block (retained bands followed by
    /// the buffer), kept only as the next warm start. `None` without a buffer.
    pub solved: Option<ComplexMatrix>,
    pub calls: Vec<SolveRecord>,
}
impl Block {
    /// The orbitals that warm-start the next map (with the buffer, if any).
    pub fn warm_start(&self) -> &ComplexMatrix {
        self.solved.as_ref().unwrap_or(&self.orbitals)
    }
}
/// FD buffer bands (2026-11-26): the reference Davidson solves `n + delta`
/// columns, `delta = max(4, ceil(n / 10))`, but requires convergence and
/// computes occupations, coverage, density, energies and delivery only for the
/// leading `n`. The width is bounded by the smallest k basis. Fixed-sector
/// policies solve exactly `n` (no buffer).
pub(super) fn buffer_bands(bands: usize) -> usize {
    4.max(bands.div_ceil(10))
}
pub(super) fn fd_width(bands: usize, basis: usize) -> usize {
    (bands + buffer_bands(bands)).min(basis).max(bands)
}
pub(super) fn solve_width(policy: &Policy, bands: usize, basis: usize) -> usize {
    if policy.finite().is_some() {
        fd_width(bands, basis)
    } else {
        bands
    }
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
    pub xc_response_ha: f64,
    pub xc_remainder_ha: f64,
    pub fd_number_slope_e_per_ha: f64,
    pub xc_mask_crossings: usize,
    /// FD empty-state coverage numbers when this map misses the criterion
    /// (fixed band count only; with automatic growth such a map is rejected
    /// and solved again at more bands). `None` when covered or not FD.
    pub coverage: Option<Shortfall>,
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
pub(super) struct Fields {
    pub(super) hxc: Vec<Vec<f64>>,
    hartree: f64,
    xc: f64,
    xc_integral: f64,
    vh: Vec<f64>,
    nlcc: Option<crate::energy::CoreXcLedger>,
}
pub(super) fn fields<K: Collective>(
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
    // The scalar XC is discarded for spin; use the existing Hartree assembly
    // without paying for a second, unpolarized PBE gradient/divergence.
    let scalar_functional = if density.len() == 2 && functional == Functional::Pbe {
        Functional::LdaPw
    } else {
        functional
    };
    // This scalar XC is discarded for spin. Only the actual channel XC may
    // consume the screening budget or publish its field audit.
    let scalar_scope = (density.len() == 2 && dft_physics::xc::screen_negative())
        .then(|| dft_physics::xc::screening_scope(dft_physics::xc::screening_stage(), true, None));
    let scalar = super::map::fields(&mut state.common, comm, cell, scalar_functional, &total)?;
    drop(scalar_scope);
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
    let screening = dft_physics::xc::screen_negative();
    if screening {
        let dv = cell.volume_bohr3() / state.grid.dims().iter().product::<usize>() as f64;
        for (s, n) in density.iter().enumerate() {
            let combined: Vec<_> = n
                .iter()
                .enumerate()
                .map(|(i, n)| {
                    n + state
                        .common
                        .core_density
                        .as_ref()
                        .map_or(0., |c| 0.5 * c[i])
                })
                .collect();
            dft_physics::gga::distributed::audit(
                &state.common.fields,
                comm,
                &combined,
                dv,
                if s == 0 { "up" } else { "down" },
                Some(state.electrons),
            )?;
        }
    }
    let _scope = screening.then(|| {
        dft_physics::xc::screening_scope(
            dft_physics::xc::screening_stage(),
            true,
            Some(state.electrons),
        )
    });
    if let Some(pbe) = &mut workspace.pbe {
        let core = state.common.core_density.as_deref();
        let out = pbe.evaluate(
            &mut state.common.fields,
            comm,
            [&density[0], &density[1]],
            core,
        )?;
        let dv = cell.volume_bohr3() / state.grid.dims().iter().product::<usize>() as f64;
        let values = sum(
            comm,
            &[
                out.energy,
                out.valence_potential_integral,
                out.core_potential_integral,
                dv * core.map_or(0., |c| c.iter().sum::<f64>()),
            ],
        )?;
        let potentials = out
            .potentials
            .iter()
            .map(|v| v.iter().zip(&scalar.vh).map(|(x, h)| x + h).collect())
            .collect();
        return Ok(Fields {
            hxc: potentials,
            vh: scalar.vh,
            hartree: scalar.hartree,
            xc: values[0],
            xc_integral: values[1],
            nlcc: core.map(|_| crate::energy::CoreXcLedger {
                core_electrons: values[3],
                core_potential_integral_ha: values[2],
                valence_potential_integral_ha: values[1],
            }),
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
/// The FD empty-state coverage rejection (criterion: top-band occupation
/// fraction above `EDGE_FRACTION` or edge electrons above `EDGE_ELECTRONS` N
/// at any k). GB1 (2026-11-27): with a fixed (explicit) band count it is
/// applied to the delivered state only (the stopped map, the verification
/// map, the last map of an exhausted budget); a transient SCF map that misses
/// it is recorded as a diagnostic. Automatic growth reacts to every map.
pub(super) const COVERAGE: &str = "occupations.bands: insufficient empty-state coverage at explicit cap; increase cap and repeat the convergence study";
/// Top-band occupation fraction (f / channel capacity) limit of the coverage
/// criterion: 1e-8, about 18.4 theta above the chemical potential.
pub(super) const EDGE_FRACTION: f64 = 1e-8;
/// Edge electrons (k-weighted top-band electrons) limit relative to max(N, 1).
pub(super) const EDGE_ELECTRONS: f64 = 1e-10;
/// Coverage numbers of an FD map that misses the criterion.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub(super) struct Shortfall {
    pub edge_max_fraction: f64,
    pub edge_electrons: f64,
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
    response_diagnostics: bool,
    shortfall: Option<&mut Shortfall>,
    bands: usize,
) -> Result<Map, ScfError> {
    let fields = fields(workspace, &mut states[0], comm, cell, functional, input)?;
    // Solved columns: the retained `bands` plus the FD buffer (`solve_width`).
    let width = states[0].orbitals.cols();
    if bands == 0 || bands > width {
        return Err(ScfError::InvalidInput(
            "electronic reference map: retained bands exceed the solved block",
        ));
    }
    // Local k blocks: every k without pools, this rank's own with GK9 pools.
    let nk = states[0].blocks.len();
    let ns = states.len();
    let electrons = states[0].electrons;
    let dims = states[0].grid.dims();
    let points = dims.iter().product::<usize>();
    let dv = cell.volume_bohr3() / points as f64;
    let pool = states[0].pool.clone();
    let pooled = pool.is_some();
    let layout = super::k_pool::KPool::comm_of(pool.as_ref());
    let kc = super::k_pool::KComm::new(comm, &layout);
    let field_axis = states[0].fields.first_axis();
    // Channel potentials and their mean (collectives on every rank) before
    // the rank-local solves; with pools also gathered to the full grid.
    let mut potentials = Vec::with_capacity(ns);
    for (s, state) in states.iter().enumerate() {
        let potential: Vec<_> = state
            .ionic
            .iter()
            .zip(&fields.hxc[s])
            .map(|(a, b)| a + b)
            .collect();
        let v0 = sum(comm, &[potential.iter().sum::<f64>()])?[0] / points as f64;
        let full = if pooled {
            Some(super::k_pool::FullField::gather(
                comm, dims, field_axis, &potential,
            )?)
        } else {
            None
        };
        potentials.push((potential, v0, full));
    }
    let solved = (|| -> Result<Vec<Vec<Block>>, ScfError> {
        let mut blocks = Vec::with_capacity(ns);
        for (s, (state, (potential, v0, full))) in
            states.iter_mut().zip(potentials.iter_mut()).enumerate()
        {
            let v0 = *v0;
            let mut channel = Vec::with_capacity(nk);
            for k in 0..nk {
                let mut op_potential = Vec::new();
                if let Some(full) = full.as_mut() {
                    op_potential.extend_from_slice(full.for_operator(
                        state.blocks[k].operator.first_axis(),
                        kc.rank(),
                        kc.size(),
                    )?);
                } else {
                    let mut planes = crate::grid_dist::PlaneMap::new(
                        dims,
                        comm.rank(),
                        comm.size(),
                        field_axis,
                        state.blocks[k].operator.first_axis(),
                    );
                    planes.a_to_b(comm, potential, &mut op_potential)?;
                }
                let op = &mut state.blocks[k].operator;
                let epoch = op
                    .owned_epoch()
                    .unwrap_or(0)
                    .checked_add(1)
                    .ok_or(ScfError::InvalidInput("spin potential epoch overflow"))?;
                op.set_owned_potential(&op_potential, epoch)
                    .map_err(numerical)?;
                let diagonal: Vec<_> = op.diagonal_static().iter().map(|d| d + v0).collect();
                let hx = kept[s].prepare(k, &kc, op, &state.orbitals.0[k], &op_potential)?;
                // FD: all retained bands are mandatory, even underflowed tails;
                // the buffer columns beyond them (`solve_width`) are not.
                // Fixed sector: first empty is mandatory as in the legacy policy.
                let strict = policy.strict_bands(s, bands);
                let mut options = ref_davidson::Options::new(width, ethr / 2., strict);
                options.first_step = first;
                options.return_hx = true;
                options.profile = diagnostics;
                let mut solution = ref_davidson_complex::solve_with_initial_hx(
                    &kc,
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
                    // FD buffer columns are never required: with them the pending
                    // count is at most delta <= max(5, width / 4).
                    if solution.converged[..strict].iter().all(|&v| v)
                        && solution.converged.iter().filter(|&&v| !v).count() <= 5.max(width / 4)
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
                        &kc,
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
                                (hx[(i, j)]
                                    - solution.eigenvalues[j] * solution.eigenvectors[(i, j)])
                                    .norm_sqr()
                            })
                            .sum()
                    })
                    .collect();
                let residuals = sum(&kc, &local)?.into_iter().map(f64::sqrt).collect();
                kept[s].retain(k, state, &solution.eigenvectors, hx, &op_potential);
                let mut eigenvalues = solution.eigenvalues;
                let (orbitals, solved) = if width > bands {
                    // Only the retained bands enter occupations, coverage,
                    // density, energies and delivery; the buffer is dropped.
                    eigenvalues.truncate(bands);
                    let mut retained = solution.eigenvectors.clone();
                    retained.resize_columns(bands).map_err(numerical)?;
                    (retained, Some(solution.eigenvectors))
                } else {
                    (solution.eigenvectors, None)
                };
                channel.push(Block {
                    eigenvalues,
                    occupations: vec![0.; bands],
                    residuals,
                    orbitals,
                    solved,
                    calls,
                });
            }
            blocks.push(channel);
        }
        Ok(blocks)
    })();
    let mut blocks = match &pool {
        None => solved?,
        Some(pool) => share_blocks(comm, pool, solved, ns, bands)?,
    };
    let weights: Vec<_> = super::k_map::k_weights(&states[0]);
    let mut coverage = None;
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
        if evaluation.edge_max_fraction > EDGE_FRACTION
            || evaluation.edge_electrons > EDGE_ELECTRONS * electrons.max(1.)
        {
            coverage = Some(Shortfall {
                edge_max_fraction: evaluation.edge_max_fraction,
                edge_electrons: evaluation.edge_electrons,
            });
        }
        // Automatic band growth (`shortfall` given) rejects the map and keeps
        // its solved orbitals as the warm start of the enlarged block; a fixed
        // count keeps the map and records the miss (`Map::coverage`).
        if let (Some(found), Some(shortfall)) = (coverage, shortfall) {
            *shortfall = found;
            for (state, channel) in states.iter_mut().zip(&blocks) {
                // GK9 pools: this rank's blocks are the global k `local()`.
                let first = state.pool.as_ref().map_or(0, |p| p.local().start);
                for (x, b) in state.orbitals.0.iter_mut().zip(channel.iter().skip(first)) {
                    x.clone_from(b.warm_start());
                }
            }
            return Err(ScfError::InvalidInput(COVERAGE));
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
    // Per (channel, local k): the weighted kinetic term of every band and the
    // nonlocal term; the densities go straight to the field planes (or, with
    // pools, into the one exchange of the map).
    let mut exchange = pool
        .as_ref()
        .map(|p| super::k_pool::DensityExchange::new(p, dims, field_axis, ns));
    let first = pool.as_ref().map_or(0, |p| p.local().start);
    let terms = (|| -> Result<Vec<(Vec<f64>, f64)>, ScfError> {
        let mut terms = Vec::with_capacity(ns * nk);
        for (s, channel) in blocks.iter().enumerate() {
            let state = &mut states[s];
            for (k, block) in channel.iter().skip(first).take(nk).enumerate() {
                let op = &mut state.blocks[k].operator;
                let f: Vec<_> = block
                    .occupations
                    .iter()
                    .map(|v| v * weights[first + k])
                    .collect();
                let band_kinetic: Vec<f64> = f
                    .iter()
                    .enumerate()
                    .map(|(j, &weight)| {
                        weight
                            * op.kinetic_diagonal()
                                .iter()
                                .enumerate()
                                .map(|(i, t)| t * block.orbitals[(i, j)].norm_sqr())
                                .sum::<f64>()
                    })
                    .collect();
                // GK9B: row partials of a pool of several ranks are summed
                // in the pool so the records are replicated there.
                let band_kinetic = if pooled && kc.size() > 1 {
                    super::k_map::pool_sum(&kc, Vec::new(), band_kinetic)?.1
                } else {
                    band_kinetic
                };
                let mut nl = vec![0.; bands];
                op.nonlocal_band_energies(&kc, &block.orbitals, &mut nl)
                    .map_err(numerical)?;
                terms.push((
                    band_kinetic,
                    nl.iter().zip(&f).map(|(a, b)| a * b).sum::<f64>(),
                ));
                let mut plane_density = vec![0.; op.plane_points()];
                // GF1 port of GK17 `K_KEPT_SLABS` (default off): keep the
                // planes of the leading retained bands for the next map's
                // kept-HX update of this channel and k (the budget is shared
                // by every channel's k; rank-independent band count).
                let pp = op.plane_points();
                let (slab, slab_bands) =
                    kept[s].slab_sink(k, ns * nk, bands, 16 * points.div_ceil(kc.size()), pp);
                op.band_density_sum_keeping(
                    &kc,
                    &block.orbitals,
                    &f,
                    cell.volume_bohr3(),
                    &mut plane_density,
                    slab,
                    slab_bands,
                )
                .map_err(numerical)?;
                kept[s].mark_slabs(k, slab_bands);
                if let (Some(ex), Some(pool)) = (exchange.as_mut(), pool.as_ref()) {
                    ex.add(pool, s, k, op.first_axis(), &plane_density)?;
                } else {
                    let mut planes = crate::grid_dist::PlaneMap::new(
                        dims,
                        comm.rank(),
                        comm.size(),
                        field_axis,
                        op.first_axis(),
                    );
                    let mut local = Vec::new();
                    planes.b_to_a(comm, &plane_density, &mut local)?;
                    for (a, b) in density[s].iter_mut().zip(local) {
                        *a += b;
                    }
                }
            }
        }
        Ok(terms)
    })();
    let terms = match &pool {
        None => terms?,
        Some(pool) => {
            let terms = channel_major(
                super::k_pool::allgather(comm, terms, pool.leader())?,
                pool,
                ns,
            )?;
            exchange
                .take()
                .expect("pooled exchange")
                .finish(comm, pool, &mut density)?;
            terms
        }
    };
    // The one-rank loop's accumulation in global (channel, k) order.
    let global_k = weights.len();
    for (s, channel) in blocks.iter().enumerate() {
        for (k, block) in channel.iter().enumerate() {
            let (band_kinetic, band_nonlocal) = &terms[s * global_k + k];
            for (j, &occupation) in block.occupations.iter().enumerate() {
                let weight = occupation * weights[k];
                band_sum += weight * block.eigenvalues[j];
                guard += weight * block.residuals[j].powi(2);
                residual_bound += weight * block.residuals[j];
                kinetic += band_kinetic[j];
            }
            nonlocal += band_nonlocal;
        }
    }
    // GK9 pools: the replayed kinetic term is already global on every rank and
    // enters a tree sum of rank partials below; only rank 0 contributes it
    // (x + 0 = x exactly).
    if pooled && comm.rank() != 0 {
        kinetic = 0.;
    }
    if let Some(s) = &mut workspace.symmetry {
        s.apply(&mut states[0].fields, comm, &mut density)?;
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
    // Diagnose the complete functional along the actual map displacement.
    // These potentials already include PBE's gradient/divergence response;
    // subtract Hartree before contracting each channel, without cancellation.
    let response_diagnostics = response_diagnostics
        && functional == Functional::Pbe
        && matches!(
            policy,
            Policy::SpinFermiDirac {
                constraint: dft_physics::occupations::Constraint::Total,
                ..
            }
        );
    let (xc_response_ha, xc_remainder_ha, fd_number_slope_e_per_ha, xc_mask_crossings) =
        if response_diagnostics {
            let local = density
                .iter()
                .zip(input)
                .enumerate()
                .map(|(s, (a, b))| {
                    a.iter()
                        .zip(b)
                        .enumerate()
                        .map(|(i, (a, b))| {
                            let dvxc = (output.hxc[s][i] - output.vh[i])
                                - (fields.hxc[s][i] - fields.vh[i]);
                            ((a - b) * dvxc).abs()
                        })
                        .sum::<f64>()
                })
                .sum::<f64>();
            let mut crossings = 0.;
            if dft_physics::xc::screen_negative() {
                for (a, b) in density.iter().zip(input) {
                    for (i, (a, b)) in a.iter().zip(b).enumerate() {
                        let core = states[0].core_density.as_ref().map_or(0., |c| 0.5 * c[i]);
                        let t = dft_physics::xc::DENSITY_THRESHOLD;
                        crossings += f64::from((a + core > t) != (b + core > t));
                    }
                }
            }
            let values = sum(comm, &[local, crossings])?;
            let theta = policy.finite().unwrap().0;
            let mu = evaluation.unwrap().roots[0].unwrap().mu_ha;
            // dN/dmu at this frozen spectrum. Use the exponential tail, retaining
            // holes even when a stored majority occupation rounds to one.
            let slope = blocks
                .iter()
                .map(|channel| {
                    channel
                        .iter()
                        .zip(&weights)
                        .map(|(b, w)| {
                            b.eigenvalues
                                .iter()
                                .map(|e| {
                                    let t = (-((e - mu) / theta).abs()).exp();
                                    w * t / ((1. + t).powi(2) * theta)
                                })
                                .sum::<f64>()
                        })
                        .sum::<f64>()
                })
                .sum::<f64>();
            (
                0.5 * dv * values[0],
                (output.xc - fields.xc - (remainder[0] - remainder[1])).abs(),
                slope,
                values[1] as usize,
            )
        } else {
            (0., 0., 0., 0)
        };
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
        xc_response_ha,
        xc_remainder_ha,
        fd_number_slope_e_per_ha,
        xc_mask_crossings,
        coverage,
    })
}
/// GK9 pools: records arrive rank-major (each rank's channel-major block of
/// its own k points); reorder to global channel-major (channel, k) order.
fn channel_major<T>(
    records: Vec<T>,
    pool: &super::k_pool::KPool,
    channels: usize,
) -> Result<Vec<T>, ScfError> {
    if records.len() != channels * pool.nk {
        return Err(ScfError::InvalidInput("k pools: record count"));
    }
    let mut slots: Vec<Option<T>> = records.into_iter().map(Some).collect();
    let mut out = Vec::with_capacity(slots.len());
    for s in 0..channels {
        for k in 0..pool.nk {
            let r = pool.owner(k);
            let at = channels * pool.first[r] + s * pool.counts[r] + (k - pool.first[r]);
            out.push(slots[at].take().expect("each record once"));
        }
    }
    Ok(out)
}
/// GK9 pools: allgather the solve records of every rank's own (channel, k)
/// blocks; the orbitals stay with their owner (other k get an empty block).
fn share_blocks<K: Collective>(
    comm: &K,
    pool: &super::k_pool::KPool,
    local: Result<Vec<Vec<Block>>, ScfError>,
    channels: usize,
    bands: usize,
) -> Result<Vec<Vec<Block>>, ScfError> {
    let (mut owned, records) = match local {
        Ok(blocks) => {
            let records = blocks
                .iter()
                .flatten()
                .map(|b| super::k_map::PerK {
                    eigenvalues: b.eigenvalues.clone(),
                    residuals: b.residuals.clone(),
                    converged: Vec::new(),
                    calls: b.calls.clone(),
                })
                .collect();
            (blocks, Ok(records))
        }
        Err(e) => (Vec::new(), Err(e)),
    };
    let records = channel_major(
        super::k_pool::allgather(comm, records, pool.leader())?,
        pool,
        channels,
    )?;
    let local = pool.local();
    let mut records = records.into_iter();
    let mut out = Vec::with_capacity(channels);
    owned.resize_with(channels, Vec::new);
    for mine in owned {
        let mut mine = mine.into_iter();
        let mut channel = Vec::with_capacity(pool.nk);
        for k in 0..pool.nk {
            let record = records.next().expect("channel-major records");
            channel.push(if local.contains(&k) {
                mine.next().expect("owned block")
            } else {
                Block {
                    eigenvalues: record.eigenvalues,
                    occupations: vec![0.; bands],
                    residuals: record.residuals,
                    orbitals: ComplexMatrix::zeros(0, 0).map_err(numerical)?,
                    solved: None,
                    calls: record.calls,
                }
            });
        }
        out.push(channel);
    }
    Ok(out)
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
                pbe: None,
                spin: Some(SpinFields {
                    xc: SpinXc::new(SpinFunctional::LdaPw, &cell, dims).unwrap(),
                    response: SpinXcResponse::new(dims),
                    channels: std::array::from_fn(|_| Array3::zeros(dims)),
                    core: Some(Array3::from_elem(dims, core)),
                    kernel_difference: vec![0.; 8],
                    kernel_span: vec![0.; 8],
                    kernel_indices: (0..8).collect(),
                }),
                symmetry: None,
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

#[cfg(test)]
mod g4f_stiffness_tests {
    use super::*;
    use crate::refscf::metric::Serial;
    use dft_numerics::switches::RunConfig;

    fn workspace() -> Workspace {
        let cell = Cell::from_user_rows_bohr([[2., 0., 0.], [0., 2., 0.], [0., 0., 2.]]).unwrap();
        Workspace {
            pbe: None,
            symmetry: None,
            spin: Some(SpinFields {
                xc: SpinXc::new(SpinFunctional::LdaPw, &cell, [2; 3]).unwrap(),
                response: SpinXcResponse::new([2; 3]),
                channels: std::array::from_fn(|_| Array3::zeros([2; 3])),
                core: Some(Array3::from_elem([2; 3], 0.02)),
                kernel_difference: vec![0.; 8],
                kernel_span: vec![0.; 8],
                kernel_indices: (17..25).collect(),
            }),
        }
    }

    #[test]
    fn g4f_frozen_mask_and_unscreened_bits() {
        let config = RunConfig::from_pairs(&[("xc_screen_negative", "on")]).unwrap();
        let _scope = config.enter();
        let mut work = workspace();
        let t = dft_physics::xc::DENSITY_THRESHOLD;
        // Cancel the core exactly at vacuum samples, including negative input.
        for u in [-0.03, -0.01, -0.01 + 0.5 * t] {
            let input = vec![vec![u; 8], vec![-0.02; 8]];
            let output = vec![vec![u + 1e-5; 8], vec![-0.02; 8]];
            assert_eq!(
                work.spin_stiffness(&Serial, &input, &output, 1.).unwrap(),
                0.
            );
        }
        work.spin.as_mut().unwrap().core = None;
        for down in [-1e-14, 0., t, f64::from_bits(t.to_bits() + 1), 1.01 * t] {
            let input = vec![vec![0.04; 8], vec![down; 8]];
            let output = vec![vec![0.04 + 1e-5; 8], vec![down; 8]];
            let value = work.spin_stiffness(&Serial, &input, &output, 1.).unwrap();
            assert!(value.is_finite() && value > 0.);
            let spin = work.spin.as_ref().unwrap();
            for &probe in &spin.channels[1] {
                if down <= t {
                    assert_eq!(probe, 0.);
                } else {
                    assert!(probe > t, "probe crossed the frozen mask");
                }
            }
        }
        assert!(
            dft_physics::xc::screening_records().is_empty(),
            "probes must not publish density audits"
        );
        // A screened channel contributes no potential derivative. Check the
        // surviving channel's curvature against an independent energy stencil.
        let input = vec![vec![0.04; 8], vec![0.; 8]];
        let dm = 1e-5;
        let output = vec![vec![0.04 + dm / 2.; 8], vec![-dm / 2.; 8]];
        let measured = work.spin_stiffness(&Serial, &input, &output, 1.).unwrap();
        let spin = work.spin.as_mut().unwrap();
        let mut energy = |m: f64| {
            spin.channels[0].fill(0.04 + m / 2.);
            spin.channels[1].fill(0.);
            spin.xc
                .evaluate(&spin.channels, &mut spin.response)
                .unwrap();
            spin.response.energy_ha
        };
        let h = 4e-6;
        let curvature = (energy(h) + energy(-h) - 2. * energy(0.)) / (h * h);
        let expected = 0.5 * curvature.abs() * dm * dm;
        assert!((measured / expected - 1.).abs() < 1e-5);
        let input = vec![vec![0.04; 8], vec![0.02; 8]];
        let output = vec![vec![0.04001; 8], vec![0.01999; 8]];
        let on = work.spin_stiffness(&Serial, &input, &output, 1.).unwrap();
        let off = RunConfig::from_pairs(&[("xc_screen_negative", "off")]).unwrap();
        let _scope = off.enter();
        let off = work.spin_stiffness(&Serial, &input, &output, 1.).unwrap();
        assert_eq!(on.to_bits(), off.to_bits());
    }

    #[test]
    fn g4f_rank_local_error_reaches_empty_owner_with_stage_and_global_index() {
        for stage in ["start", "map", "verification"] {
            dft_mpi::simulate(3, |comm| {
                let config = RunConfig::from_pairs(&[("xc_screen_negative", "on")]).unwrap();
                let _scope = config.enter();
                let _stage = dft_physics::xc::screening_scope(stage, false, None);
                let mut work = workspace();
                let mut input = vec![vec![0.04; 8], vec![0.02; 8]];
                if comm.rank() == 1 {
                    input[1][3] = f64::NAN;
                }
                if comm.rank() == 2 {
                    work.spin = None;
                    input = vec![vec![], vec![]];
                }
                let error = work
                    .spin_stiffness(comm, &input, &input, 1.)
                    .unwrap_err()
                    .to_string();
                for expected in [
                    stage,
                    "rank 1",
                    "channel down minimum NaN",
                    "grid point 20",
                    "rank-local reason",
                    "nonfinite combined input density",
                ] {
                    assert!(error.contains(expected), "{error}");
                }
            });
        }
    }
}
