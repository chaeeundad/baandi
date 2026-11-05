# ADR 0033: Grid roles, energy ownership and derivative ledger

- Status: **accepted by the user 2026-11-05** (contract approved; implementation capability still absent until its ticket delivers it); P13A-STATE-CONTRACTS / MF02.
- Basis: [master plan D5, D6, D8, D10–12 and MF12–16](../roadmap/master-plan-proposal.md), [physics conventions](../physics-conventions.md), [ADR 0017](0017-coupling-commitments.md), [ADR 0019](0019-optional-mkl-backend.md).

## Grid contract

`GridDescriptor` fields: `id`, `role:wavefunction|density|augmentation_reserved`,
`cell_columns_bohr:[3,3]`, `origin_fractional:[3]`, `dimensions:[3]` positive
integers, `reciprocal_cutoff_ha`, `ordering:C_z_fastest`,
`normalization:forward_1_over_N_inverse_unscaled`, `periodic:[true,true,true]`,
`quadrature_weight_bohr3:Omega/product(dimensions)`, `basis_order_digest`, and
`derivative_convention_id` (including Nyquist handling). Wavefunction role
transforms periodic Bloch parts; it does not include Omega^-1/2 automatically.
Density role contains total/channel density and local potentials. Different role
IDs may share dimensions/storage only under explicit aliasing ownership.
Augmentation is reserved and rejects under NC-first scope U3. The present
4*Ecut density sphere does not imply arbitrary independent real grids (D5).

`GridTransfer` fields: `source_id`, `target_id`, `method:fourier_embed_truncate`,
`common_miller_order_digest`, `origin_phase_convention`,
`adjoint_id`, `charge_error`, `discarded_norm`, `validation_version`.
Let P map orbital-grid samples to density-grid samples and W be grid quadrature.
The reverse potential action is P†=W_w^-1 P* W_d; test
<Pu,v>_d=<u,P†v>_w. Preserve the constant Fourier mode and phase for origin
changes. Truncation is lossy and is not an exact inverse. Transfer for changed
cells is not admitted by this same-cell contract. An explicit future warm-transfer
plan must preserve charge and record interpolation error (MF16).

Plans and scratch belong to `WorkerWorkspace{worker_id,grid_ids,basis_ids,
backend_id,plan_generation}`; one mutable workspace per active worker, immutable
model shared. Prepare plans outside kernels, rebuild on grid/basis/backend
change, never serialize handles. Scheduling owns the total budget; library
threads do not multiply it. No GPU/MPI buffers or runtime backend switching
are introduced (D10–11, ADR 0019). Legacy single-grid runs keep their path.

## EnergyContribution and capabilities

Each entry has `term_id`, `owner`, `energy_ha`, `included_in:E|correction|audit_only`,
`state_digest`, `force_capability`, `stress_capability`, and
`dependencies:[term_id]`. Capabilities are `unavailable{reason}` or future
`available{implementation_id,validated_model_tags,evidence_id}`. Availability
requires all active model terms and a fresh stationary certificate. Today every
force/stress capability below is unavailable; this document does not implement
one. `not_applicable` is allowed only for an inactive contribution, never an
active term with a missing derivative. This ledger applies the energy audit in
[physics conventions](../physics-conventions.md) and D6.

| Contribution | Single proposed owner | E / audit role | Required displacement derivative | Required strain derivative; delivery gate |
| --- | --- | --- | --- | --- |
| Kinetic | physics kinetic, engine sums w*f | E | Stationarity cancels orbital response in complete fixed plane-wave basis | Reciprocal metric, normalization and cutoff/basis treatment; MF14 |
| Local electron-ion, including G=0 convention | physics local pseudopotential | E | Ionic phase/radial dependence; MF13 | Volume, G metric and local transform; MF14 |
| Nonlocal NC | physics projectors | E | Projector center derivatives for every channel; MF13 | Projector radial/angular and cell/basis dependence; MF14 |
| Hartree | physics Hartree | E, factor 1/2 | No explicit ionic coordinate at fixed valence density; cancellation must be demonstrated with total variational derivative | Coulomb kernel, G=0 and volume; MF14 |
| XC LDA/PBE | physics XC | E | NLCC moving core if active; pure valence variational cancellation; MF13 | Grid/gradient metric, volume, core; MF14 |
| Ewald real, reciprocal, self, background | physics Ewald | E summed once | Full ion-ion derivative, self/background invariance checked; MF13 | All four pieces, eta dependence cancels within truncation bounds; MF14 |
| Hubbard interaction minus double counting | physics Hubbard (MF10) | E once; convention/projectors in model | Projector and occupation-matrix dependence; MF13 before U forces advertised | Strain/projector/overlap terms; MF14 |
| Entropic term -TS | engine occupation policy | correction E -> F | Fixed theta,N stationary free-energy derivative; no duplicate df term | Same at fixed theta,N; MF14 |
| External electron-energy field, ion/environment terms | separately named physics/host owners (MF19, ADR 0017) | Only approved functional; no incomplete coupled total | Unavailable until each active coordinate dependency is supplied | Unavailable until each cell dependency is supplied |
| Band sum, integral n V_H, sum_s integral n_s V_xc, Hubbard potential subtraction, Harris-Foulkes remainder | engine energy audit | audit_only, never extra physical E | Diagnostics, not independent forces | Diagnostics, not independent stress |

New functional terms must add a ledger entry and invalidate total derivative
availability until covered. PAW/USPP overlap, meta-GGA differential operators,
nonlocal exchange and charged reservoir corrections remain rejected/reserved
(D4, D7, D12), rather than represented as zero-valued NC contributions.

## Force and stress semantics

For Q=E (fixed insulating occupations) or F (FD at fixed theta,N and any explicit
fixed M), force F_Ia=-partial Q/partial R_Ia at fixed cell, Ha/Bohr, Cartesian
ordered by input atoms. `ForceResult` contains status, reason, quantity Q,
state_digest, atom_order_digest, units, frame and optional `[atom,3]` values;
values are absent when unavailable. Constraints are recorded separately: raw
physical forces must not be overwritten by projected optimizer forces.

Use tensile-positive Cauchy stress sigma_ab=(1/Omega)partial Q/partial eta_ab
for A(eta)=(I+eta)A and fixed fractional ionic coordinates. Units Ha/Bohr^3,
Cartesian tensor `[3,3]`; delta Q=Omega sum_ab sigma_ab delta eta_ab.
Pressure=-trace(sigma)/3, positive in compression. For a symmetric shear with
eta_xy=eta_yx=h, derivative dQ/dh=2 Omega sigma_xy. No implicit Voigt factor.
External-code pressure/stress signs are converted explicitly before comparison.

`StressResult` adds `basis_protocol:frozen_miller_set|fixed_cutoff_converged`,
`cutoff_ha`, `strain_steps`, and `pulay_assessment` to force-style metadata.
A fixed-cell force test has no cell/basis change. A frozen Miller-set strain test
changes reciprocal vectors and normalization but keeps the integer G list and
FFT dimensions; it tests a differentiable finite-basis functional, **not** the
production fixed-cutoff stress limit. A fixed-cutoff strain scan rebuilds the
basis, can jump at shell crossings and requires cutoff/grid convergence and
reported Pulay/discretization error. Never mix these protocols in one tolerance.

MF15/16 optimizers must request complete derivatives of the same Q and state,
stop on unavailable/unconverged results, retain failed steps, and distinguish
exact restart from cell-transfer warm starts. See [stage oracles](../validation/stage-oracles.md)
for multi-step displacement/strain tests and [migration](../architecture/state-contracts.md#migration-and-invalid-combinations)
for explicit rejections.
