//! Post-SCF derivatives of E, or F at fixed temperature/electron constraints.
//! No occupation-response term is added: stationarity of F cancels it together
//! with orbital/density response. The supplied model must be the state model.
use crate::{
    ScfError, ScfOptions, ScfResult, SpeciesPotential, model::collinear::CollinearResult,
    model::finite_temperature::FiniteTResult,
};
use dft_basis::FftGridSpec;
use dft_geometry::{Atom, Cell, kpoints::KPoint};
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::{
    data::*,
    ewald, forces as terms,
    hamiltonian::{NonlocalOperator, Species},
};

#[derive(Clone, Debug)]
pub struct DerivativeEntry {
    pub active: bool,
    pub energy_ha: Option<f64>,
    pub included_in: &'static str,
    pub dependencies: Vec<&'static str>,
    pub term_id: &'static str,
    pub owner: &'static str,
    pub force_available: bool,
    pub stress_available: bool,
    pub reason: &'static str,
}
/// Every physical term in the supported fixed-cell NC functional. Audit-only
/// band sums and double-counting integrals are not independent energy terms.
pub fn derivative_ledger() -> Vec<DerivativeEntry> {
    [
        (
            "kinetic",
            "physics::kinetic",
            "stationary fixed plane-wave basis",
        ),
        (
            "local",
            "physics::forces",
            "ionic structure-factor derivative",
        ),
        (
            "nonlocal",
            "physics::forces",
            "Hermitian projector derivative",
        ),
        (
            "hartree",
            "physics::hartree",
            "no explicit ionic dependence",
        ),
        (
            "xc",
            "physics::xc",
            "valence response cancels by stationarity",
        ),
        (
            "core",
            "physics::forces",
            "moving core against XC potential",
        ),
        ("ewald", "physics::ewald", "real and reciprocal derivatives"),
        (
            "entropy",
            "engine::occupations",
            "stationary free energy at fixed theta,N,M",
        ),
    ]
    .into_iter()
    .map(|(term_id, owner, reason)| DerivativeEntry {
        active: true,
        energy_ha: None,
        included_in: "E",
        dependencies: Vec::new(),
        term_id,
        owner,
        reason,
        force_available: true,
        stress_available: false,
    })
    .collect()
}
/// Check the complete active ledger before calculating any derivative. Hosts
/// adding a physical term must append its capability here, never omit it.
pub fn require_complete(ledger: &[DerivativeEntry]) -> Result<(), String> {
    let missing: Vec<_> = ledger
        .iter()
        .filter(|e| e.active && !e.force_available)
        .map(|e| e.term_id)
        .collect();
    if missing.is_empty() {
        Ok(())
    } else {
        Err(format!(
            "observables.forces: incomplete derivative coverage: {}; disable this request or select a validated combination",
            missing.join(", ")
        ))
    }
}

#[derive(Debug)]
pub struct ForceResult {
    pub verification: &'static str,
    pub ledger: Vec<DerivativeEntry>,
    pub units: &'static str,
    pub quantity: &'static str,
    pub values: Vec<[f64; 3]>,
    pub ewald: Vec<[f64; 3]>,
    pub local: Vec<[f64; 3]>,
    pub nonlocal: Vec<[f64; 3]>,
    pub core: Vec<[f64; 3]>,
    pub net_force: [f64; 3],
    pub maximum_absolute_component: f64,
    pub elapsed: std::time::Duration,
    pub band_passes: usize,
    /// GD1: wall seconds per stage of this evaluation, in execution order
    /// (observation only; per rank).
    pub timings: Vec<(&'static str, f64)>,
}
/// Explicit borrowed stationary electronic state/model view. This API deliberately
/// does not accept a density-only checkpoint as a stationary electronic state.
/// Callers must retain the exact cell, atoms, potentials and options used to
/// produce `electronic`; matching dimensions alone do not establish identity.
pub enum ElectronicState<'a> {
    Fixed(&'a ScfResult),
    Collinear(&'a CollinearResult),
    FiniteTemperature(&'a FiniteTResult),
    /// Reference spin/Fermi-Dirac delivery (`refscf::electronic`): scalar FD,
    /// collinear fixed sector and collinear FD. Occupations carry the channel
    /// capacity (2 unpolarized, 1 per spin channel) and exclude k weights.
    Reference(&'a crate::refscf::electronic::Run),
}
/// One occupied k block of a derivative: point (with weight), complete
/// replicated coefficients in the `options.basis(cell, k, 1)` order, and the
/// occupations without k weight.
pub(crate) type DerivativeBlock<'a> = (KPoint, &'a ComplexMatrix, &'a [f64]);

/// MPI-DERIV: evaluate the per-block part `eval` of a derivative (a fixed
/// number `stride` of values per block) and hand every part to `accumulate` in
/// global block order. Without an active communicator (threads-only path) the
/// blocks are evaluated and accumulated one after the other, exactly as the
/// former single loop. With two or more ranks each rank evaluates a contiguous
/// range of blocks (`dft_mpi::block_counts`), and one allgatherv delivers all
/// parts to every rank; the accumulation then runs in the same global order on
/// the same values, so the result is bit-identical to the replicated
/// evaluation and independent of the rank count. Requires the replicated
/// delivered state (every rank holds every block), which the reference
/// delivery guarantees. Each rank's segment starts with a status value so a
/// failure on one rank is reported on all ranks instead of leaving the others
/// waiting in the collective.
pub(crate) fn accumulate_blocks(
    blocks: &[DerivativeBlock<'_>],
    stride: usize,
    mut eval: impl FnMut(&DerivativeBlock<'_>, &mut [f64]) -> Result<(), ScfError>,
    mut accumulate: impl FnMut(&[f64]),
) -> Result<(), ScfError> {
    use dft_mpi::Collective;
    let Some(comm) = crate::parallel::SharedComm::active() else {
        let mut part = vec![0.0; stride];
        for block in blocks {
            eval(block, &mut part)?;
            accumulate(&part);
        }
        return Ok(());
    };
    let counts = dft_mpi::block_counts(blocks.len(), comm.size());
    let first = dft_mpi::block_displacements(&counts)[comm.rank()];
    let own = &blocks[first..first + counts[comm.rank()]];
    let mut local = vec![0.0; 1 + own.len() * stride];
    let mut failure = None;
    if stride > 0 {
        for (block, slot) in own.iter().zip(local[1..].chunks_exact_mut(stride)) {
            if let Err(e) = eval(block, slot) {
                failure = Some(e);
                local[0] = 1.0;
                break;
            }
        }
    }
    let counts: Vec<_> = counts.iter().map(|c| 1 + c * stride).collect();
    let mut all = vec![0.0; counts.iter().sum()];
    comm.allgatherv_f64(&local, &counts, &mut all)
        .map_err(|e| ScfError::Numerical(format!("derivative block gather: {e}")))?;
    if let Some(e) = failure {
        return Err(e);
    }
    let mut at = 0;
    for &count in &counts {
        if all[at] != 0.0 {
            return Err(ScfError::Numerical(
                "derivative block evaluation failed on another rank".into(),
            ));
        }
        at += count;
    }
    if stride == 0 {
        return Ok(());
    }
    let mut at = 0;
    for count in counts {
        for values in all[at + 1..at + count].chunks_exact(stride) {
            accumulate(values);
        }
        at += count;
    }
    Ok(())
}
impl<'a> ElectronicState<'a> {
    /// Channel densities (one unpolarized, up/down collinear), electrons/Bohr^3.
    pub(crate) fn densities(&self) -> Result<Vec<&'a Array3<f64>>, ScfError> {
        Ok(match *self {
            Self::Fixed(r) => vec![&r.density],
            Self::FiniteTemperature(r) => vec![&r.state.density],
            Self::Collinear(r) => {
                let dft_physics::spin::SpinDensity::Collinear { channels, .. } = &r.density else {
                    return Err(ScfError::InvalidInput(
                        "forces.spin: collinear density required",
                    ));
                };
                channels.iter().collect()
            }
            Self::Reference(r) => r.density.iter().collect(),
        })
    }
    /// Every (spin, k) block with its occupations, spin-major.
    pub(crate) fn blocks(&self) -> Vec<DerivativeBlock<'a>> {
        match *self {
            Self::Fixed(r) => r
                .k_points
                .iter()
                .map(|k| (k.point, &k.solution.eigenvectors, r.occupations.as_slice()))
                .collect(),
            Self::FiniteTemperature(r) => r
                .state
                .k_points
                .iter()
                .enumerate()
                .map(|(k, p)| (p.point, &p.solution.eigenvectors, r.occupations.block(0, k)))
                .collect(),
            Self::Collinear(r) => r
                .k_points
                .iter()
                .enumerate()
                .flat_map(|(s, ks)| {
                    ks.iter().enumerate().map(move |(k, p)| {
                        (
                            p.point,
                            &p.solution.eigenvectors,
                            r.finite_occupations
                                .as_ref()
                                .map_or(r.occupations[s].as_slice(), |f| f.block(s, k)),
                        )
                    })
                })
                .collect(),
            Self::Reference(r) => r
                .blocks
                .iter()
                .flatten()
                .map(|b| (b.point, &b.orbitals, b.occupations.as_slice()))
                .collect(),
        }
    }
    /// Sampled k points of the first channel (all channels share the mesh).
    pub(crate) fn points(&self) -> Vec<KPoint> {
        match *self {
            Self::Fixed(r) => r.k_points.iter().map(|k| k.point).collect(),
            Self::FiniteTemperature(r) => r.state.k_points.iter().map(|k| k.point).collect(),
            Self::Collinear(r) => r.k_points[0].iter().map(|k| k.point).collect(),
            Self::Reference(r) => r.blocks[0].iter().map(|b| b.point).collect(),
        }
    }
}
pub struct ForceState<'a> {
    pub cell: &'a Cell,
    pub atoms: &'a [Atom],
    pub potentials: &'a [SpeciesPotential<'a>],
    pub options: &'a ScfOptions,
    pub electronic: ElectronicState<'a>,
    pub density_floor: Option<f64>,
    /// Exact spatial operations used by the SCF, empty when symmetry is off.
    pub symmetry: &'a [dft_geometry::symmetry::SymmetryOperation],
}

/// Outcome of the shared derivative preconditions ([`validate`]).
pub(crate) struct Validated {
    pub(crate) verification: &'static str,
    pub(crate) ledger: Vec<DerivativeEntry>,
    pub(crate) quantity: &'static str,
}

pub fn forces(state: &ForceState<'_>) -> Result<ForceResult, ScfError> {
    let start = std::time::Instant::now();
    let validated = validate(state)?;
    let mut timings = vec![("validate", start.elapsed().as_secs_f64())];
    let mut clock = std::time::Instant::now();
    let mut lap = |name: &'static str, timings: &mut Vec<(&'static str, f64)>| {
        timings.push((name, clock.elapsed().as_secs_f64()));
        clock = std::time::Instant::now();
    };
    let Validated {
        verification,
        ledger,
        quantity,
    } = validated;
    let densities = state.electronic.densities()?;
    let blocks = state.electronic.blocks();
    let mut density = densities[0].clone();
    if densities.len() == 2 {
        density += densities[1];
    }
    let volume = state.cell.volume_bohr3();
    let density_basis = state.options.basis(state.cell, [0.; 3], 4.0)?;
    let dims: [usize; 3] = density
        .shape()
        .try_into()
        .map_err(|_| ScfError::InvalidInput("forces.density: invalid shape"))?;
    let dg = FftGridSpec::new(&density_basis, &density_basis, Some(dims))?;
    // One Species entry per input atom prevents grouping from permuting output.
    let positions: Vec<_> = state
        .atoms
        .iter()
        .map(|a| {
            state
                .cell
                .fractional_to_cartesian_bohr(a.fractional_position.map(|x| x.rem_euclid(1.0)))
        })
        .collect();
    let species: Vec<_> = state
        .atoms
        .iter()
        .zip(&positions)
        .map(|(a, p)| {
            let potential = state
                .potentials
                .iter()
                .find(|s| s.label == a.species)
                .ok_or(ScfError::InvalidInput("forces.species: missing potential"))?
                .potential;
            Ok(Species {
                potential,
                positions_bohr: std::slice::from_ref(p),
            })
        })
        .collect::<Result<_, ScfError>>()?;
    let mut fft = FftPlan3d::new(dims)?;
    let mut reciprocal = Array3::zeros(dims);
    fft.to_reciprocal(&density.mapv(|x| Complex64::new(x, 0.0)), &mut reciprocal)?;
    lap("setup", &mut timings);
    let shells = shell_values(&species, &density_basis, volume, false)?;
    lap("shells", &mut timings);
    let mut local = terms::local_or_core_cached(
        &density_basis,
        volume,
        &species,
        &dg.gather(&reciprocal)?,
        false,
        shells.as_ref(),
    )?;
    lap("local", &mut timings);
    let core_density = dft_physics::density::core_density_cached(
        &density_basis,
        &dg,
        volume,
        &species,
        &mut fft,
        shells.as_ref(),
    )?;
    let mut core = if let Some(core_density) = core_density {
        let mut scratch = Array3::zeros(dims);
        let potential = if densities.len() == 1 {
            state
                .options
                .evaluate_xc(
                    &density,
                    Some(&core_density),
                    &mut scratch,
                    state.cell,
                    &mut fft,
                    state.density_floor,
                )?
                .potential
        } else {
            use dft_physics::spin::{SpinFunctional, SpinXc, SpinXcResponse};
            if state.density_floor.is_some() {
                return Err(ScfError::InvalidInput(
                    "forces.spin: density floor unavailable",
                ));
            }
            let kind = match state.options.functional {
                crate::Functional::LdaPw => SpinFunctional::LdaPw,
                crate::Functional::Pbe => SpinFunctional::Pbe,
                _ => {
                    return Err(ScfError::InvalidInput(
                        "forces.spin: polarized PZ unavailable",
                    ));
                }
            };
            let mut xc = SpinXc::new(kind, state.cell, dims)?;
            let mut out = SpinXcResponse::new(dims);
            xc.evaluate_with_core(
                &[densities[0].clone(), densities[1].clone()],
                Some(&core_density),
                &mut out,
            )?;
            (&out.potentials[0] + &out.potentials[1]) * 0.5
        };
        lap("core_xc", &mut timings);
        fft.to_reciprocal(&potential.mapv(|x| Complex64::new(x, 0.0)), &mut reciprocal)?;
        terms::local_or_core_cached(
            &density_basis,
            volume,
            &species,
            &dg.gather(&reciprocal)?,
            true,
            shells.as_ref(),
        )?
    } else {
        vec![[0.0; 3]; state.atoms.len()]
    };
    lap("core", &mut timings);
    let ions: Vec<_> = state
        .atoms
        .iter()
        .zip(&species)
        .map(|(a, s)| ewald::Ion {
            charge: s.potential.z_valence(),
            fractional_position: a.fractional_position,
        })
        .collect();
    let mut ewald = ewald::forces(state.cell, &ions, ewald::EwaldOptions::default())?;
    lap("ewald", &mut timings);
    let mut nonlocal = vec![[0.0; 3]; state.atoms.len()];
    let mut band_passes = 0;
    let band_first = dft_numerics::switches::on(dft_numerics::switches::Id::FORCE_BAND_FIRST);
    let gemm_kernel = dft_numerics::switches::on(dft_numerics::switches::Id::DERIV_NL_GEMM);
    let tables = if gemm_kernel {
        projector_tables(state, &species, blocks.len(), false)?
    } else {
        Vec::new()
    };
    let per_atom = tables_per_atom(&species, &tables);
    lap("nonlocal_tables", &mut timings);
    // MPI-DERIV: per-k parts distributed over ranks, accumulated in k order.
    accumulate_blocks(
        &blocks,
        3 * state.atoms.len(),
        |(point, orbitals, occupations), out| {
            let basis = state.options.basis(state.cell, point.fractional, 1.0)?;
            let occupations: Vec<_> = occupations.iter().map(|f| f * point.weight).collect();
            let part = if gemm_kernel {
                dft_physics::stress::nonlocal_k_derivatives(
                    &basis,
                    volume,
                    &species,
                    &per_atom,
                    orbitals,
                    &occupations,
                    false,
                )?
                .0
            } else {
                let operator = NonlocalOperator::new(&basis, volume, &species)?;
                if band_first {
                    nonlocal_band_first(&operator, &basis, orbitals, &occupations)?
                } else {
                    terms::nonlocal(&operator, &basis, orbitals, &occupations)?
                }
            };
            if part.len() != out.len() / 3 {
                return Err(ScfError::InvalidInput("nonlocal forces: atom count"));
            }
            for (slot, value) in out.as_chunks_mut::<3>().0.iter_mut().zip(part) {
                slot.copy_from_slice(&value);
            }
            Ok(())
        },
        |part| {
            for (sum, part) in nonlocal.iter_mut().zip(part.as_chunks::<3>().0) {
                for axis in 0..3 {
                    sum[axis] += part[axis];
                }
            }
        },
    )?;
    band_passes += blocks.iter().map(|b| b.2.len()).sum::<usize>();
    lap("nonlocal", &mut timings);
    for term in [&mut ewald, &mut local, &mut nonlocal, &mut core] {
        symmetrize(state.cell, term, state.symmetry)?;
    }
    let values: Vec<[f64; 3]> = (0..state.atoms.len())
        .map(|i| std::array::from_fn(|a| ewald[i][a] + local[i][a] + nonlocal[i][a] + core[i][a]))
        .collect();
    if values.iter().flatten().any(|f| !f.is_finite()) {
        return Err(ScfError::NonFinite);
    }
    let net_force: [f64; 3] = std::array::from_fn(|a| values.iter().map(|f| f[a]).sum());
    if net_force.iter().any(|f| !f.is_finite()) {
        return Err(ScfError::NonFinite);
    }
    let maximum_absolute_component = values.iter().flatten().map(|v| v.abs()).fold(0.0, f64::max);
    lap("symmetrize", &mut timings);
    Ok(ForceResult {
        verification,
        ledger,
        units: "Ha/Bohr",
        quantity,
        values,
        ewald,
        local,
        nonlocal,
        core,
        net_force,
        maximum_absolute_component,
        elapsed: start.elapsed(),
        band_passes,
        timings,
    })
}

/// GD1: projector radial tables (values and, with `slopes`, q derivatives)
/// per distinct potential of `species`, in first-appearance order, when the
/// derivative pass has two or more k blocks and `DERIV_RADIAL_TABLES` is on;
/// otherwise empty (the kernel uses the quadrature). The range covers
/// |k+G| <= 1.1 sqrt(2 E_cut) + 0.05 Bohr^-1 (a strained frozen-Miller set
/// may exceed it; such rows fall back to the quadrature). Node values are the
/// quadrature itself; (channel, node) units are split into contiguous rank
/// blocks of an active communicator and shared by one allgatherv, so the
/// tables are bit-identical for any rank count.
pub(crate) fn projector_tables(
    state: &ForceState<'_>,
    species: &[Species<'_>],
    blocks: usize,
    slopes: bool,
) -> Result<Vec<(usize, dft_physics::stress::ProjectorRadialTables)>, ScfError> {
    use dft_upf::radial::{
        RADIAL_TABLE_POINTS as POINTS, RADIAL_TABLE_STEP as STEP, RadialTransformTable,
    };
    if blocks < 2 || !dft_numerics::switches::on(dft_numerics::switches::Id::DERIV_RADIAL_TABLES) {
        return Ok(Vec::new());
    }
    let q_max = 1.1 * (2.0 * state.options.cutoff_hartree).sqrt() + 0.05;
    let nodes = RadialTransformTable::nodes(q_max, STEP, POINTS)
        .map_err(|e| ScfError::Numerical(e.to_string()))?;
    let mut potentials: Vec<&dft_upf::Pseudopotential> = Vec::new();
    for s in species {
        if !potentials.iter().any(|p| std::ptr::eq(*p, s.potential)) {
            potentials.push(s.potential);
        }
    }
    // (potential, projector, slope?)
    let mut channels: Vec<(usize, usize, bool)> = Vec::new();
    for (p, pp) in potentials.iter().enumerate() {
        for i in 0..pp.projectors().len() {
            channels.push((p, i, false));
            if slopes {
                channels.push((p, i, true));
            }
        }
    }
    let total = channels.len() * nodes;
    let comm = crate::parallel::SharedComm::active();
    let (rank, size) = comm.as_ref().map_or((0, 1), |c| {
        (dft_mpi::Collective::rank(c), dft_mpi::Collective::size(c))
    });
    // Contiguous unit blocks per rank, one allgatherv: every value is the
    // quadrature of exactly one rank, copied unchanged to all ranks.
    let counts = dft_mpi::block_counts(total, size);
    let first = dft_mpi::block_displacements(&counts)[rank];
    let mut local = vec![0.0; 1 + counts[rank]];
    for (slot, u) in local[1..].iter_mut().zip(first..first + counts[rank]) {
        let ((p, i, slope), n) = (channels[u / nodes], u % nodes);
        let q = n as f64 * STEP;
        let value = if slope {
            dft_physics::stress::projector_transform_slope(potentials[p], i, q)
                .map_err(|e| e.to_string())
        } else {
            potentials[p]
                .projector_transform(i, q)
                .map_err(|e| e.to_string())
        };
        match value {
            Ok(v) => *slot = v,
            Err(_) => {
                local[0] = 1.0;
                break;
            }
        }
    }
    let values = if let Some(comm) = &comm {
        let counts: Vec<_> = counts.iter().map(|c| 1 + c).collect();
        let mut all = vec![0.0; counts.iter().sum()];
        dft_mpi::Collective::allgatherv_f64(comm, &local, &counts, &mut all)
            .map_err(|e| ScfError::Numerical(format!("derivative radial tables: {e}")))?;
        let mut values = Vec::with_capacity(total);
        let mut at = 0;
        for count in counts {
            if all[at] != 0.0 {
                local[0] = 1.0;
            }
            values.extend_from_slice(&all[at + 1..at + count]);
            at += count;
        }
        values
    } else {
        local[1..].to_vec()
    };
    if local[0] != 0.0 {
        return Err(ScfError::Numerical(
            "derivative radial tables: quadrature failed".into(),
        ));
    }
    let mut out: Vec<(usize, dft_physics::stress::ProjectorRadialTables)> = potentials
        .iter()
        .map(|pp| {
            (
                *pp as *const _ as usize,
                dft_physics::stress::ProjectorRadialTables {
                    value: Vec::new(),
                    slope: Vec::new(),
                },
            )
        })
        .collect();
    for (c, &(p, i, slope)) in channels.iter().enumerate() {
        let l = potentials[p].projectors()[i].l();
        // The slope of an l transform has the opposite parity.
        let table = RadialTransformTable::from_nodes(
            if slope { l + 1 } else { l },
            q_max,
            STEP,
            POINTS,
            values[c * nodes..(c + 1) * nodes].to_vec(),
        )
        .map_err(|e| ScfError::Numerical(e.to_string()))?;
        if slope {
            out[p].1.slope.push(table);
        } else {
            out[p].1.value.push(table);
        }
    }
    Ok(out)
}

/// GD1: the local and core radial transforms (and with `slopes` their q
/// derivatives) on every distinct nonzero |G| of the density sphere, for each
/// distinct potential, evaluated once per derivative pass with hoisted
/// integrands (`dft_physics::radial_shells::ShellEvaluator`, every value
/// bit-identical to the former per-shell expression). With an active
/// communicator the (potential, function, shell) units are split into
/// contiguous rank blocks and combined by one allgatherv, so every rank holds
/// the identical values for any rank count. `None` when `DERIV_SHELLS` is off.
pub(crate) fn shell_values(
    species: &[Species<'_>],
    density_basis: &dft_basis::PlaneWaveBasis,
    volume: f64,
    slopes: bool,
) -> Result<Option<dft_physics::radial_shells::RadialShellValues>, ScfError> {
    use dft_mpi::Collective;
    use dft_physics::radial_shells::{RadialShellValues, ShellEvaluator, ShellKind};
    if !dft_numerics::switches::on(dft_numerics::switches::Id::DERIV_SHELLS) {
        return Ok(None);
    }
    let mut potentials: Vec<&dft_upf::Pseudopotential> = Vec::new();
    for s in species {
        if !potentials.iter().any(|p| std::ptr::eq(*p, s.potential)) {
            potentials.push(s.potential);
        }
    }
    let mut seen = std::collections::HashSet::new();
    let shells: Vec<f64> = density_basis
        .vectors()
        .iter()
        .filter(|g| g.miller != [0; 3])
        .map(|g| g.norm_squared.sqrt())
        .filter(|q| seen.insert(q.to_bits()))
        .collect();
    let mut units: Vec<(usize, ShellKind)> = Vec::new();
    for (p, pp) in potentials.iter().enumerate() {
        units.push((p, ShellKind::LocalValue));
        if slopes {
            units.push((p, ShellKind::LocalSlope));
        }
        if pp.core_density().is_some() {
            units.push((p, ShellKind::CoreValue));
            if slopes {
                units.push((p, ShellKind::CoreSlope));
            }
        }
    }
    let total = units.len() * shells.len();
    let comm = crate::parallel::SharedComm::active();
    let (rank, size) = comm.as_ref().map_or((0, 1), |c| (c.rank(), c.size()));
    let counts = dft_mpi::block_counts(total, size);
    let first = dft_mpi::block_displacements(&counts)[rank];
    let mut local = vec![0.0; 1 + counts[rank]];
    let mut evaluators: Vec<Option<ShellEvaluator<'_>>> = potentials.iter().map(|_| None).collect();
    for (slot, u) in local[1..].iter_mut().zip(first..first + counts[rank]) {
        let ((p, kind), q) = (units[u / shells.len()], shells[u % shells.len()]);
        if evaluators[p].is_none() {
            // A failure is reported through the status value, never by an
            // early return that would leave the other ranks in the gather.
            match ShellEvaluator::new(potentials[p], volume) {
                Ok(e) => evaluators[p] = Some(e),
                Err(_) => {
                    local[0] = 1.0;
                    break;
                }
            }
        }
        match evaluators[p].as_mut().expect("evaluator").eval(kind, q) {
            Ok(Some(v)) => *slot = v,
            _ => {
                local[0] = 1.0;
                break;
            }
        }
    }
    let all = if let Some(comm) = &comm {
        let counts: Vec<_> = counts.iter().map(|c| 1 + c).collect();
        let mut all = vec![0.0; counts.iter().sum()];
        comm.allgatherv_f64(&local, &counts, &mut all)
            .map_err(|e| ScfError::Numerical(format!("derivative shell gather: {e}")))?;
        let mut values = Vec::with_capacity(total);
        let mut at = 0;
        for count in counts {
            if all[at] != 0.0 {
                // A failing shell is recomputed (and reported) by the kernels.
                return Ok(None);
            }
            values.extend_from_slice(&all[at + 1..at + count]);
            at += count;
        }
        values
    } else {
        if local[0] != 0.0 {
            return Ok(None);
        }
        local.split_off(1)
    };
    let mut out = RadialShellValues::default();
    for (u, v) in all.into_iter().enumerate() {
        let ((p, kind), q) = (units[u / shells.len()], shells[u % shells.len()]);
        out.insert(potentials[p], kind, q, v);
    }
    Ok(Some(out))
}

/// Per-entry view of [`projector_tables`] for a species list.
pub(crate) fn tables_per_atom<'t>(
    species: &[Species<'_>],
    tables: &'t [(usize, dft_physics::stress::ProjectorRadialTables)],
) -> Vec<Option<&'t dft_physics::stress::ProjectorRadialTables>> {
    species
        .iter()
        .map(|s| {
            let id = s.potential as *const _ as usize;
            tables.iter().find(|t| t.0 == id).map(|t| &t.1)
        })
        .collect()
}

/// Shared preconditions of every derivative (forces and frozen stress): the
/// independent setup identity, the certificate or passed verification, the
/// strict tolerance caps, complete state arrays and the complete active
/// ledger. Computes no derivative.
pub(crate) fn validate(state: &ForceState<'_>) -> Result<Validated, ScfError> {
    if let ElectronicState::Fixed(r) = &state.electronic
        && !r.orbitals_delivered()
    {
        return Err(ScfError::InvalidInput(
            "forces need the delivered orbitals: call ScfResult::materialize_orbitals",
        ));
    }
    if let ElectronicState::Reference(r) = &state.electronic
        && !r.orbitals_delivered()
    {
        return Err(ScfError::InvalidInput(
            "forces need the delivered orbitals: call electronic::Run::materialize_orbitals",
        ));
    }
    let crystal_density = match &state.electronic {
        ElectronicState::Fixed(r) => r.crystal_density_group,
        ElectronicState::FiniteTemperature(r) => r.state.crystal_density_group,
        ElectronicState::Collinear(_) | ElectronicState::Reference(_) => false,
    };
    if state.options.crystal_density_group || crystal_density {
        return Err(ScfError::InvalidInput(
            "symmetry.density_group=crystal does not support forces or stress",
        ));
    }
    if state.options.local_mesh != crate::LocalMesh::Baandi
        || state.options.radial_even_rule != crate::RadialEvenRule::Simpson3_8
    {
        return Err(ScfError::InvalidInput(
            "research pseudo conventions do not support analytic derivatives",
        ));
    }

    use crate::force_identity::{ModelDescription, OccupationIdentity};
    let (identity, policy, density, points) = match &state.electronic {
        ElectronicState::Fixed(r) => (
            &r.force_identity,
            OccupationIdentity::new(
                &vec![
                    r.occupations.as_slice();
                    if r.force_identity.reference_delivery && !r.gamma_real_used {
                        r.k_points.len()
                    } else {
                        1
                    }
                ],
                None,
                None,
            ),
            &r.density,
            r.k_points.iter().map(|k| k.point).collect(),
        ),
        ElectronicState::FiniteTemperature(r) => (
            &r.force_identity,
            OccupationIdentity::new(&[], Some(&r.occupations), None),
            &r.state.density,
            r.state.k_points.iter().map(|k| k.point).collect(),
        ),
        ElectronicState::Collinear(r) => {
            let dft_physics::spin::SpinDensity::Collinear { channels, .. } = &r.density else {
                return Err(ScfError::InvalidInput(
                    "observables.forces.model: spin model mismatch; rerun SCF with the requested model",
                ));
            };
            if r.k_points[0].iter().map(|k| k.point).collect::<Vec<_>>()
                != r.k_points[1].iter().map(|k| k.point).collect::<Vec<_>>()
            {
                return Err(ScfError::InvalidInput(
                    "observables.forces.model: spin k-point mismatch; rerun SCF with the requested model",
                ));
            }
            (
                &r.force_identity,
                OccupationIdentity::new(
                    &[&r.occupations[0], &r.occupations[1]],
                    r.finite_occupations.as_ref(),
                    Some(&r.sector),
                ),
                &channels[0],
                r.k_points[0].iter().map(|k| k.point).collect(),
            )
        }
        ElectronicState::Reference(r) => {
            if r.blocks.is_empty()
                || r.density.is_empty()
                || r.blocks.len() != r.density.len()
                || r.blocks.iter().any(|c| {
                    c.iter()
                        .map(|b| b.point)
                        .ne(r.blocks[0].iter().map(|b| b.point))
                })
            {
                return Err(ScfError::InvalidInput(
                    "observables.forces.model: spin k-point mismatch; rerun SCF with the requested model",
                ));
            }
            (
                &r.force_identity,
                r.occupation_identity(),
                &r.density[0],
                state.electronic.points(),
            )
        }
    };
    let dims: [usize; 3] = density
        .shape()
        .try_into()
        .map_err(|_| ScfError::InvalidInput("observables.forces.model: invalid grid; rerun SCF"))?;
    let mut model = ModelDescription::new(
        state.cell,
        state.atoms,
        state.potentials,
        state.options,
        dims,
        points,
    );
    model.floor = state.density_floor;
    model.symmetry = state.symmetry.to_vec();
    if !identity.matches(&model, &policy) || state.options.grid_override.is_some_and(|g| g != dims)
    {
        return Err(ScfError::InvalidInput(
            "observables.forces.model: model or occupation policy differs from the certified state; use its original model or rerun SCF with the requested model",
        ));
    }

    // ADR 0043: a standard state with a passed fixed-point verification counts
    // as certified; an unverified standard state is rejected as before.
    let verified = match &state.electronic {
        ElectronicState::Fixed(r) => {
            r.standard.is_some() && r.verification.as_ref().is_some_and(|v| v.passed())
        }
        // ONEPATH-3: spin and Fermi-Dirac states on the standard path.
        ElectronicState::FiniteTemperature(r) => {
            r.converged
                && r.standard.is_some()
                && r.verification.as_ref().is_some_and(|v| v.passed())
        }
        ElectronicState::Collinear(r) => {
            r.converged
                && r.standard.is_some()
                && r.verification.as_ref().is_some_and(|v| v.passed())
        }
        // The reference electronic driver always runs the standard profile.
        ElectronicState::Reference(r) => {
            r.converged
                && state.options.standard.is_some()
                && r.verification.as_ref().is_some_and(|v| v.passed())
        }
    };
    let (certified, standard, quantity) = match &state.electronic {
        ElectronicState::Fixed(r) => (
            verified
                || (r.certified
                    && r.history
                        .last()
                        .and_then(|h| h.fixed_point_certificate.as_ref())
                        .is_some()),
            r.standard.is_some() && !verified,
            "E",
        ),
        ElectronicState::FiniteTemperature(r) => (
            verified || (r.converged && r.certificate.as_ref().is_some_and(|c| c.passed)),
            r.standard.is_some() && !verified,
            "F",
        ),
        ElectronicState::Collinear(r) => {
            let certified = verified
                || (r.converged
                    && if r.finite_occupations.is_some() {
                        r.finite_certificate.as_ref().is_some_and(|c| c.passed)
                    } else {
                        r.certificate.as_ref().is_some_and(|c| c.passed)
                    });
            (
                certified,
                r.standard.is_some() && !verified,
                if r.finite_occupations.is_some() {
                    "F"
                } else {
                    "E"
                },
            )
        }
        // ADR 0043: the passed reference verification is the certificate.
        ElectronicState::Reference(r) => (verified, !verified, r.quantity()),
    };
    let densities = state.electronic.densities()?;
    let blocks = state.electronic.blocks();
    let unverified = match &state.electronic {
        ElectronicState::Fixed(r) => {
            r.force_identity.unverified_delivery
                && state.options.verification == crate::Verification::Off
                && r.standard.is_some()
                && r.verification.is_none()
        }
        ElectronicState::Reference(r) => {
            r.force_identity.unverified_delivery
                && state.options.verification == crate::Verification::Off
                && r.converged
                && r.verification.is_none()
        }
        _ => false,
    };
    if !unverified && ((state.options.standard.is_some() && !verified) || standard) {
        return Err(ScfError::InvalidInput(
            "observables.forces: strict convergence profile required, or a passed fixed-point verification on the standard path; rerun with scf.verify = fixed_point (default) or on the strict reference path (BAANDI_REFERENCE_STRICT=on)",
        ));
    }
    if state.options.density_tolerance > 1e-7
        || state.options.energy_tolerance > 1e-8
        || state.options.eigen_tolerance > 1e-8
    {
        return Err(ScfError::InvalidInput(
            "observables.forces: certificate tolerances exceed strict caps; tighten density<=1e-7 and energy/eigen<=1e-8 Ha",
        ));
    }
    if !certified && !unverified {
        return Err(ScfError::InvalidInput(
            "observables.forces: passing fixed-point certificate required; converge the state",
        ));
    }
    if state.atoms.is_empty()
        || blocks.is_empty()
        || densities[0].is_empty()
        || densities.iter().any(|d| d.shape() != densities[0].shape())
    {
        return Err(ScfError::InvalidInput(
            "observables.forces: incomplete state arrays; use the complete converged electronic state",
        ));
    }
    let (energy, entropy) = match &state.electronic {
        ElectronicState::Fixed(r) => (&r.energy, None),
        ElectronicState::FiniteTemperature(r) => {
            (&r.state.energy, Some(r.thermodynamics.entropy_term_ha))
        }
        ElectronicState::Collinear(r) => (&r.energy, r.thermodynamics.map(|t| t.entropy_term_ha)),
        ElectronicState::Reference(r) => (
            &r.energy,
            r.policy
                .is_finite()
                .then_some(r.thermodynamics.entropy_term_ha),
        ),
    };
    // Intentionally exhaustive: adding an energy/audit field must trigger a
    // force-coverage review instead of silently retaining a partial total.
    let crate::energy::EnergyAudit {
        xc_floor: _,
        nlcc,
        kinetic,
        local,
        local_zero_mode: _,
        nonlocal,
        hartree,
        xc,
        ewald,
        band_sum: _,
        double_counting: _,
        band_form: _,
        harris_foulkes: _,
        functional_difference_bound: _,
    } = energy;
    let mut ledger = derivative_ledger();
    for entry in &mut ledger {
        entry.energy_ha = Some(match entry.term_id {
            "kinetic" => *kinetic,
            "local" => *local,
            "nonlocal" => *nonlocal,
            "hartree" => *hartree,
            "xc" => *xc,
            "ewald" => *ewald,
            "entropy" => {
                entry.active = entropy.is_some();
                entry.included_in = "correction";
                -entropy.unwrap_or(0.0)
            }
            "core" => {
                entry.active = nlcc.is_some();
                entry.included_in = "audit_only";
                entry.dependencies.push("xc");
                nlcc.as_ref().map_or(0.0, |c| c.core_potential_integral_ha)
            }
            _ => unreachable!("closed NC ledger"),
        });
    }
    require_complete(&ledger).map_err(ScfError::Numerical)?;
    Ok(Validated {
        verification: if unverified {
            "unverified"
        } else if verified {
            "verified"
        } else {
            "strict_certificate"
        },
        ledger,
        quantity,
    })
}

/// Stack beta and its three position derivatives into one projector matrix.
/// One adjoint GEMM per k evaluates all overlaps for every band and atom.
/// This changes reduction order only; occupations already include k weights.
pub fn nonlocal_band_first(
    operator: &NonlocalOperator,
    basis: &dft_basis::PlaneWaveBasis,
    orbitals: &ComplexMatrix,
    occupations: &[f64],
) -> Result<Vec<[f64; 3]>, ScfError> {
    use dft_numerics::linalg::{MatrixOp, gemm};
    let ng = basis.momenta().len();
    if orbitals.rows() != ng || orbitals.cols() != occupations.len() {
        return Err(ScfError::InvalidInput("nonlocal forces: orbital shape"));
    }
    if occupations.iter().any(|f| !f.is_finite() || *f < 0.) {
        return Err(ScfError::NonFinite);
    }
    let count: usize = operator.atoms().iter().map(|a| a.columns().len()).sum();
    let mut forces = vec![[0.; 3]; operator.atoms().len()];
    if count == 0 {
        return Ok(forces);
    }
    let mut projectors =
        ComplexMatrix::zeros(ng, 4 * count).map_err(|e| ScfError::Numerical(e.to_string()))?;
    let mut index = 0;
    for atom in operator.atoms() {
        for column in atom.columns() {
            if column.len() != ng {
                return Err(ScfError::InvalidInput("nonlocal forces: projector shape"));
            }
            for (g, (&beta, momentum)) in column.iter().zip(basis.momenta()).enumerate() {
                projectors[(g, index)] = beta;
                for axis in 0..3 {
                    projectors[(g, (axis + 1) * count + index)] =
                        Complex64::new(0., -momentum.cartesian[axis]) * beta;
                }
            }
            index += 1;
        }
    }
    let overlaps = gemm(&projectors, MatrixOp::Adjoint, orbitals, MatrixOp::Identity)
        .map_err(|e| ScfError::Numerical(e.to_string()))?;
    let mut offset = 0;
    for (atom, force) in operator.atoms().iter().zip(&mut forces) {
        let n = atom.columns().len();
        for (band, &occupation) in occupations.iter().enumerate() {
            if occupation == 0. {
                continue;
            }
            for p in 0..n {
                let weighted: Complex64 = (0..n)
                    .map(|q| atom.coupling_hartree()[p * n + q] * overlaps[(offset + q, band)])
                    .sum();
                for (axis, component) in force.iter_mut().enumerate() {
                    *component -= 2.
                        * occupation
                        * (overlaps[((axis + 1) * count + offset + p, band)].conj() * weighted).re;
                }
            }
        }
        offset += n;
    }
    if forces.iter().flatten().any(|f| !f.is_finite()) {
        return Err(ScfError::NonFinite);
    }
    Ok(forces)
}

/// Polar-vector transformation A R A^-1 and the SCF atom permutation. No
/// optimizer constraint or subtraction of a residual net force is performed.
fn symmetrize(
    cell: &Cell,
    values: &mut Vec<[f64; 3]>,
    ops: &[dft_geometry::symmetry::SymmetryOperation],
) -> Result<(), ScfError> {
    if ops.is_empty() {
        return Ok(());
    }
    let mut result = vec![[0.0; 3]; values.len()];
    for op in ops {
        if op.atom_permutation.len() != values.len() {
            return Err(ScfError::InvalidInput(
                "forces.symmetry: invalid permutation",
            ));
        }
        let mut seen = vec![false; values.len()];
        for (i, v) in values.iter().enumerate() {
            let target = op.atom_permutation[i];
            if target >= values.len() || seen[target] {
                return Err(ScfError::InvalidInput(
                    "forces.symmetry: invalid permutation",
                ));
            }
            seen[target] = true;
            let f = cell.cartesian_bohr_to_fractional(*v);
            let rotated = cell.fractional_to_cartesian_bohr(std::array::from_fn(|a| {
                (0..3).map(|b| op.rotation[a][b] as f64 * f[b]).sum()
            }));
            for (a, component) in rotated.iter().enumerate() {
                result[target][a] += component / ops.len() as f64;
            }
        }
    }
    *values = result;
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    /// MPI-DERIV: global k order on every rank count, and a failing block on
    /// one rank is reported on every rank (no rank left in the collective).
    #[test]
    fn distributed_blocks_keep_order_and_report_failures() {
        let matrix = ComplexMatrix::zeros(1, 1).unwrap();
        let occupations = [2.0];
        let blocks: Vec<DerivativeBlock<'_>> = (0..5)
            .map(|i| {
                (
                    KPoint {
                        fractional: [i as f64, 0., 0.],
                        weight: 1.,
                    },
                    &matrix,
                    &occupations[..],
                )
            })
            .collect();
        for ranks in [2, 3, 7] {
            let results = dft_mpi::simulate(ranks, |comm| {
                crate::parallel::set_communicator(Some(std::sync::Arc::new(comm.clone())));
                let mut order = Vec::new();
                let ok = accumulate_blocks(
                    &blocks,
                    2,
                    |b, out| {
                        out.copy_from_slice(&[b.0.fractional[0], -b.0.fractional[0]]);
                        Ok(())
                    },
                    |p| order.push(p.to_vec()),
                );
                let failed = accumulate_blocks(
                    &blocks,
                    2,
                    |b, _| {
                        if b.0.fractional[0] == 3. {
                            Err(ScfError::NonFinite)
                        } else {
                            Ok(())
                        }
                    },
                    |_| {},
                );
                crate::parallel::set_communicator(None);
                (ok.is_ok(), order, failed.is_err())
            });
            for (ok, order, failed) in results {
                assert!(ok && failed);
                let expected: Vec<_> = (0..5).map(|i| vec![i as f64, -(i as f64)]).collect();
                assert_eq!(order, expected);
            }
        }
    }
    #[test]
    fn symmetry_rotates_polar_vectors_and_permutes_atoms() {
        use dft_geometry::symmetry::{SpinAction, SymmetryOperation};
        let cell = Cell::from_user_rows_bohr([[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]]).unwrap();
        let identity = SymmetryOperation {
            rotation: [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            translation: [0.; 3],
            atom_permutation: vec![0, 1],
            spin_action: SpinAction::Identity,
            antiunitary: false,
            model_invariance_proof: "test",
        };
        let inversion = SymmetryOperation {
            rotation: [[-1, 0, 0], [0, -1, 0], [0, 0, -1]],
            atom_permutation: vec![1, 0],
            ..identity.clone()
        };
        let mut values = vec![[1., 2., 3.], [2., 4., 6.]];
        symmetrize(&cell, &mut values, &[identity, inversion]).unwrap();
        assert_eq!(values, vec![[-0.5, -1., -1.5], [0.5, 1., 1.5]]);
    }
}
