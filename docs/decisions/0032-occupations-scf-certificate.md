# ADR 0032: Occupation policy, entropy and finite-temperature certificate

- Status: **accepted by the user 2026-11-05** (contract approved; implementation capability still absent until its ticket delivers it); P13A-STATE-CONTRACTS / MF02.
- Basis: [master plan D2, MF04](../roadmap/master-plan-proposal.md), [physics conventions §25–26](../physics-conventions.md), [ADR 0025](0025-preconditioned-density-mixing.md).

## Proposed policy contract

A pure `OccupationPolicy.evaluate` boundary accepts `BandLayout`, all sorted
`eigenvalues_ha:[s,k,b]`, electron target N, spin constraint, policy parameters,
`electron_tolerance`, `root_tolerance_ha`, `max_root_iterations`, and previous
chemical potentials as optional numerical hints. It returns either the record
below or a typed failure; hints cannot change the selected policy or certify a
failed root. No file access, allocation policy hidden in kernels, or temperature
inference from material identity. Workspace allocation belongs to preparation.

| Return field | Type and meaning |
| --- | --- |
| `occupations_electrons` | f64 [s,k,b], includes channel capacity, excludes k weight |
| `chemical_potential` | `not_unique{reason:fixed_insulator}` or `common{mu_ha}` or `per_spin{up_ha,down_ha}`; per-spin only for fixed M |
| `entropy_kind` | `zero`, `physical_fermi_dirac`, or future `generalized_smearing` |
| `entropy_over_kb` | f64 dimensionless for FD/fixed, absent for future generalized schemes unless a named definition is accepted |
| `entropy_term_ha` | positive TS for FD; 0 for fixed; signed generalized correction for later schemes, explicitly not physical TS |
| `electron_count`, `electron_residual` | sum w*f and signed sum w*f-N, electron count |
| `moment_residual` | signed M-M_target for fixed moment, otherwise absent; not a demand for zero spontaneous magnetization |
| `root_iterations`, `root_bracket_ha` | u64; absent/common/per-spin bracket matching chemical-potential variant |
| `edge_max_fraction`, `edge_electrons` | maximum top-band f/capacity over s,k and sum_sk w_k f_sk,last; diagnostics, not a bound on all omitted states |
| `policy_id`, `policy_version` | explicit reproducible definition |

Policy tags: `fixed_insulator{channel_electrons, gap_threshold_ha}`;
`fermi_dirac{temperature_energy_ha}`; reserved
`gaussian{width_ha}`, `methfessel_paxton{width_ha,order}`, and
`cold{width_ha,variant}`. The last three reject until separately validated.
FD uses theta=k_B T>0 in Ha and p=1/(1+exp((epsilon-mu)/theta)), f=g_s p.
S/k_B=-sum_skb w_k g_s[p log p+(1-p) log(1-p)] with endpoint limits zero.
A monotone bracketed number solve is required, including bracket expansion,
overflow-safe tails, capacity checks and bounded failure. Fixed M uses two
number roots; unconstrained collinear uses one root. Equal energies in the same
constraint sector receive equal fractions, independent of ordering.

For a finite-mu FD root, the target in each independently constrained sector must
be strictly between zero and its weighted retained-band capacity. Endpoint
targets (including an exactly empty spin channel at fixed M) have no finite root
at positive theta and reject with `occupations.electron_target: no finite
chemical potential at this sector endpoint; use fixed occupations or an interior
target`. Never serialize infinite mu or quietly approximate an endpoint.

Gaussian/MP/cold are named integration schemes, not automatically physical
finite temperatures. Some generalized occupations can violate [0,g]; they must
not inherit FD positivity or entropy assertions. Their energy correction,
root-selection rule and density feasibility need a later accepted policy.
Zero theta is rejected for FD with a remedy selecting fixed_insulator, not a
silent policy switch. Fixed policy preserves the existing even-N legacy rule,
sampled gap, empty-band requirements of k APIs and Gamma behavior.

## Energy and certification

`ThermodynamicLedger` contains `internal_energy_ha:E`, `entropy_term_ha:TS`,
`free_energy_ha:F=E-TS`, `entropy_kind`, `temperature_energy_ha`,
`variational_quantity:internal_energy|free_energy`, and
`grand_potential:unavailable{reason:fixed_electron_number}`. For fixed policy,
TS=0, F=E, theta=0. Future grand canonical Omega=F-mu*N requires a separately
approved reservoir/environment energy ledger (D12, ADR 0017); it is not a new
alias for E. Both direct and band-sum audits evaluate E from the same output
orbitals/density/occupations, then subtract the **same** TS once to audit F.
All components refer to the same evaluated state, never mixed-input energies.

FD map: input spin density -> Hamiltonians -> eigenpairs for every retained
band -> number roots -> occupations -> output density and E,S,F. At convergence,
build a fresh Hamiltonian from that output density, repeat the **entire** map,
and test the following on that one fresh result. `FiniteTCertificate` fields
include `map_version`, `policy_digest`, raw charge/spin density distances
(Bohr^-3/2, sqrt(integral delta^2)/N), maximum retained-band eigen residual
(Ha), weighted occupation change sum w*abs(delta f) (electrons), signed number
and optional moment residuals, abs(delta F) (Ha), edge diagnostics, configured
thresholds, `gap_status:not_applicable`, and `passed`. Never use a filtered
residual as a physical convergence gate. Require all retained-band residuals,
not only an integer occupied prefix; retain the existing tightened fresh solve
factor 0.1. Root convergence and F alone do not establish a fixed point.

Proposed initial finite-T gates: raw density distances use existing 1e-7
Bohr^-3/2, retained eigen tolerance uses existing 1e-8 Ha, number/moment errors
<=1e-10*max(1,N) electrons, occupation change <=1e-8*max(1,N) electrons,
abs(delta F)<=1e-8 Ha, edge_max_fraction<=1e-8 and
edge_electrons<=1e-10*max(1,N). These new gates are proposals for MF04 tests,
not relaxed replacements of legacy gates. Record root bracket width; near an
insulating gap mu can be ill-conditioned, so mu difference is diagnostic, not
a universal convergence test. Root residual is decisive.

Finite top-band occupation cannot prove an infinite-tail bound. Certificate
labels its scope `retained_basis`; material acceptance additionally repeats with
at least four more bands per block until E/F and density agree under the MF04
oracle. Insufficient bands return `occupations.bands: insufficient empty-state
coverage; increase bands and repeat the convergence study`, never convergence
by clipping fractions. No auto-growth without a recorded policy change.
T->0 for a sampled gapped system must recover fixed-policy E/density/occupations;
its fixed-policy gap rejection remains unchanged. See [stage oracles](../validation/stage-oracles.md).
