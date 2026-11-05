# MF03–MF14 independent oracle plan

Status: **proposed**, P13A-STATE-CONTRACTS / MF02, 2026-11-05. This is a test
plan, **not executed evidence**. Authority: [master plan MF03–14 and D1–D9](../roadmap/master-plan-proposal.md),
[ADRs 0031–34](../architecture/state-contracts.md), [physics conventions](../physics-conventions.md).
Existing measured claims remain in [evidence-map](evidence-map.md). Proposed caps
below are acceptance targets with units, never observed errors. Freeze them in
the implementation ticket before seeing cross-code differences.

## Common acceptance rules

Use independent analytic formulas or central differences of the scalar functional,
not the derivative implementation reused as its own oracle. For each finite
difference use three steps h,h/2,h/4; require an O(h^2) trend before cancellation
or an already sub-cap plateau with estimated SCF/roundoff noise below one third
of the cap. Tighten SCF, cutoff, k mesh and grid independently. Failure to resolve
noise is **inconclusive**, not passing. A stable discrepancy above the frozen cap,
wrong sign, wrong axis/capacity or missing contribution is failure. No cap may be
relaxed to accommodate another code. The existing overhead-only 1e-12 Ha rule
is not a general force or material tolerance (master-plan performance section).

Every stage retains the existing legacy tests and raw fixed-point gates. Fast
analytic fixtures belong in the default gate; runs above about 10 s debug are
ignored with an explicit release command recorded by their implementation ticket
([AGENTS.md](../../AGENTS.md)). Cover cubic, skew, elongated and slab-like cells,
Gamma and unequal-weight/shifted k points, small/large band counts, varied NC
projector structures and small/large gaps; real Si alone cannot set defaults.

| Stage | Independent oracle and concrete cases | Proposed tolerance and rationale | Failure condition |
| --- | --- | --- | --- |
| MF03 collinear | Uniform exchange e_x=-(3/4)(6/pi)^(1/3)(n_up^(4/3)+n_down^(4/3)); analytic v_s=-(6/pi)^(1/3)n_s^(1/3). Densities 0.01,0.1,1 Bohr^-3, polarization 0,0.3,0.9,1; polarized correlation oracle from independently transcribed published parameterization identified in MF03. Central functional perturbations away from density screens. Spin swap, zero-moment reduction and N/M sums on suite S below; real recipe R1 requested. | Exchange energy per electron <=1e-12 Ha; v relative <=1e-8 or abs<=1e-10 Ha near zero. Same-discretization zero-spin E <=1e-10 Ha, density distance <=1e-9 Bohr^-3/2; N/M <=1e-10*max(1,N) electrons. These isolate capacity/derivative bugs below SCF tolerances; use tighter solves for reduction. | Factor-of-two error, cross-spin GGA derivative mismatch, failure of reduction or polarized XC capability; a finite reference magnetic moment is not substituted for a functional oracle. |
| MF04 occupations | Two levels [-a,+a], equal capacities, half total capacity: mu=0 and p(-a)+p(a)=1. a=0,0.01,1 Ha; theta=0.001,0.01,0.1 Ha. Degenerate level p=1/2 gives S/k_B=g ln2 per weighted level. Unequal k weights and fixed-M roots; dF/dtheta=-S/k_B on stationary solutions. Gapped S fixture theta->0 sequence and bands b,b+4,b+8; real Al R2. | Analytic f/number <=1e-12 (electrons for N); F=E-TS <=1e-12 Ha; entropy identity <=1e-10 dimensionless. Root/certificate thresholds ADR 0032. Band refinement delta F<=1e-8 Ha, density<=1e-8 Bohr^-3/2. Low-T limiting E<=1e-8 Ha once exponential occupation tail is below cap. | Root not bracketed/converged, negative FD entropy, bands pass edge test but fail expansion, fresh map disagrees, insulating rejection bypassed. |
| MF05 datasets | Independent UPF header/radial metadata inspection against provider documentation, file SHA-256 and license text; O/Ni/Al plus Mg/O option. Negative fixtures for PAW/USPP/SOC and malformed channels; no real files committed. | Exact hash/type/element/valence/license provenance agreement; no floating tolerance for identity. Header cutoff is not a convergence oracle. | Missing license/hash/semicore/XC information, unsupported file accepted, parser success labeled material qualification. |
| MF06 NC coverage | Independent radial quadrature/analytic Gaussian local+multiple nonlocal projectors, normalized atomic density integrates valence; synthetic NLCC translated-core derivative against finite difference. Licensed O/Ni/MgO R1/R3 cutoff scans, atomic scattering/log-derivative ghost checks against independently produced provider/reference data. | Radial/atomic integrals relative<=1e-10; functional derivative<=1e-7 Ha per unit perturbation on dimensionless perturbation path. Real energy target 1e-4 Ha/atom and selected eigenvalues 1e-3 Ha as existing comparison scale, with discretization uncertainty<=one third; not a transferability proof alone. | Ghost state, wrong core double count, core derivative missing, cutoff nonconvergence or unsupported NC feature coerced. If independent scattering data absent, qualification stays incomplete. |
| MF07 state identity | Independent table-driven mutations of every identity field; same model direct density start vs fresh solve on S; charge integral computed by independent sum; legacy Gamma/k regression. | Exact expected compatible/incompatible decisions and unchanged legacy bits on same build. Warm converged energy <=1e-10 Ha only for demonstrably same branch; charge import <=1e-8 relative and policy-specific stricter gates. | Changed spin/theta/U/PP/environment accepted as identical; seed changes model; warm start labeled exact; different branch called restart error without classification. |
| MF08 persistence | Independent byte-layout reader extending P11C tests, golden scalar bits including signed zero, spin/k layouts; truncation at boundaries, trailing bytes, corruption, rehashed semantic errors, duplicate fields, unsupported version and crash-before-publication; all three kinds. | Payload and exact-continuation scientific history bit identity on same build/backend/topology, excluding timings; v1 untouched. Warm energy as MF07 with measured branch/noise checks. | Silent coercion, partial file published as complete, unavailable orbital fabricated, exact downgraded, result/unconverged export restores Converged. |
| MF09 multi-start | Independent replay of fixed recipes/seeds on S and qualified R3: uniform, bounded perturbation, normalized atomic and saved density. Include deliberately failed attempt. Analytic equal-energy different-density pair and same-density orbital unitary rotation. | Deterministic attempt ordering/status under same execution; calibrate density/moment noise using >=3 repeats, use unresolved region [noise,10*noise] rather than energy-only merge. Publish actual noise; no universal numeric state threshold before measurement. | Failure dropped, distinct fixture merged by E, gauge rotation declared distinct solely by orbitals, global-minimum claim from finite starts. |
| MF10 DFT+U | Independent matrix formula in state-contracts; diagonal n=(0,0.25,0.5,1), U_eff=0,0.1,0.2 Ha; Hermitian perturbations and random unitary basis transforms. R3 AFM/FM branches only with exactly matched projector and double-counting conventions. | E_U/unitary invariance<=1e-12 Ha; central energy-potential derivative<=1e-8 Ha for dimensionless perturbation path; real total E<=1e-4 Ha/atom, moment<=1e-3 mu_B after same branch/projector check. | U=0 fails base limit; non-Hermitian matrices; dependence on m gauge; comparing unlike projector populations or advertising U forces before MF13. |
| MF11 symmetry | Full-grid independent oracle for nonmagnetic and allowed collinear S, shifted meshes and skew origin; reconstruct density with explicit phases. Inverse/group composition and k multiplicity tests; field/constraint-breaking operations must reject. | Weights sum error<=1e-14; full/reduced E<=1e-10 Ha and density distance<=1e-9 Bohr^-3/2 under tightened SCF; moment difference<=1e-9 mu_B. | Magnetic time reversal applied without proof, non-invariant Hubbard/field term reduced, wrong phase/weights or density symmetrization. |
| MF12 grids/ledger | Analytic single Fourier modes and constants on 12x14x16 and 15x12x10 grids; representable common modes, changed origin, truncated modes; direct weighted dot-product adjoint identity and Parseval; enumerated active-term coverage. | Relative Parseval/adjoint<=1e-12, constant-mode charge<=1e-12 relative. Truncation error is reported, not expected zero. Every missing derivative yields unavailable. | Transfer changes charge, P† not quadrature-adjoint, aliasing not declared, shared mutable scratch, absent derivative represented by zero. |
| MF13 forces | Central displacement differences of Q on distorted S and R1/R2/R3, h=0.01,0.005,0.0025 Bohr, all Cartesian components and translational/net-force checks. Term-isolated local/nonlocal/Ewald, LDA/PBE, spin, NLCC, U and FD combinations; QE R4 only after complete term coverage. | Synthetic derivative cap 1e-6 Ha/Bohr; translation/net force <=1e-8 Ha/Bohr at resolved grid. Real matched force cap 1e-4 Ha/Bohr with SCF/cutoff/k/grid uncertainty <=one third. Several steps separate discretization from sign/term errors. | Any active derivative missing, E differentiated instead of FD F, net force/sign failure, plateau above cap; reject unsupported combination rather than partial total. |
| MF14 stress | Central strain differences at h=0.002,0.001,0.0005 for three diagonal and three symmetric shears, fixed fractional ions; frozen Miller and rebuilt fixed-cutoff protocols separate. Cubic/skew/slab cells, Ewald/background, core/U/FD terms, R4 QE stress with explicit sign conversion. | Frozen-set synthetic cap 1e-7 Ha/Bohr^3; converged real cap 1e-5 Ha/Bohr^3, uncertainty<=one third. Tensor antisymmetric norm<=1e-10 Ha/Bohr^3 for supported rotationally invariant model. Cutoff shell jumps require a converged scan, not a relaxed cap. | Factor-two shear/sign/volume error, Pulay estimate absent, frozen-set result advertised as cutoff-converged, omission of any active term. |


| MF15 fixed-cell relaxation (P21A ticket) | Independent scalar-energy scan/parabolic minimum for a synthetic pair from compressed/stretched starts; analytic quadratic/Rosenbrock BFGS/FIRE; displaced Si and MgO, space-group preservation, exact frozen coordinates, three-step directional energy differences, interrupted ionic continuation and failed SCF retention. | Synthetic distance <=1e-4 Bohr; warm/fresh E <=1e-10 Ha at three structures; transferred charge <=1e-8 relative; stopping requires force and energy tolerances for two consecutive steps. Same-build resumed accepted trajectory is bit-identical. Use MF13 directional force caps. | Uncertified forces, standard profile accepted, silent model coercion, lost failed trials, constraints violated, warm transfer labelled exact SCF continuation, missing last good structure on failure. Gamma-only material minima are not bulk structural qualification. |

The MF15 row is added from the explicit P21A-RELAXATION ticket; the original
P13A table stopped at MF14. Measured scope and negative MgO ideal-position
curvature are recorded in [relaxation](relaxation.md).

## Fixture S: reproducible synthetic coverage

Reuse the analytic fixture definitions in `crates/dft-engine/tests/scf.rs`,
`crates/dft-engine/tests/checkpoint.rs` and `crates/dft-physics/tests/pbe.rs`;
implementation tickets must record generator/input hashes so changes are visible.
Add explicit model-independent matrix/occupation fixtures above. Proposed geometry
variants use lattice columns (in Bohr) (8,0,0),(0,8,0),(0,0,8);
(8,0,0),(1,7,0),(0.5,1,9); (6,0,0),(0,6,0),(0,0,18); and
(8,0,0),(0,8,0),(0,0,24). Use fractional atoms (0,0,0),(.27,.23,.19),
then a four-atom replication; synthetic charges 1/1 and 2/2, widths and
projector coefficients from the cited fixture generators, explicitly retained
in each test manifest. Compare Gamma, full half-shifted 2x2x2, and a
three-point synthetic list with weights .2/.3/.5 for non-symmetry integral tests.
Use 8 and 16 bands, plus capacity-boundary negative tests. These are synthetic
models, not qualified element pseudopotentials or material predictions.

## Remote reference requests: exact-input freeze protocol

No remote job is sent by this documentation ticket. The orchestrator requests
runs on dev08 through dft5090 under [reference-environment](../reference-environment.md).
Use installed QE 7.6, never install locally or read its source. Exact inputs
cannot honestly be completed before MF05 supplies licensed NC file identities.
Therefore R1–R4 are **pending request recipes**, not ready-to-run reference
bundles. The following numeric geometry/settings are fixed proposals; a request
is blocked until every symbolic dataset identifier resolves to actual filename,
SHA-256, valence, XC/NLCC/semicore flags and license record. No placeholder SHA
or invented QE result is allowed. MF03 analytic gates can proceed independently.

| Recipe / stages | Fixed geometry and electronic protocol to freeze into both inputs |
| --- | --- |
| R1 spin / MF03,06 | O2 in cubic 20 Bohr cell, positions (0,0,-1.15),(0,0,1.15) Bohr, Gamma, qualified O NC PBE; collinear total moment 2 mu_B fixed sector and zero-moment control; 24 bands per channel, wavefunction cutoffs 40/60/80 Ha, density cutoffs four times each. Repeat cell 24 Bohr for isolation. |
| R2 metal / MF04 | Primitive fcc Al columns (0,3.8,3.8),(3.8,0,3.8),(3.8,3.8,0) Bohr; atom at origin; qualified Al NC PBE, unpolarized FD theta=0.001/0.005/0.01 Ha (QE widths 0.002/0.01/0.02 Ry), full half-shifted 4^3/6^3/8^3 meshes, 16/20/24 bands, Ecut 20/30/40 Ha and density 4*Ecut. |
| R3 oxide / MF06,09,10 | Rocksalt conventional NiO cube 7.9 Bohr. Ni fractional (0,0,0),(0,.5,.5),(.5,0,.5),(.5,.5,0); O (.5,0,0),(0,.5,0),(0,0,.5),(.5,.5,.5). Qualified NC Ni/O PBE with semicore/core flags frozen. FM seeds (+,+,+,+) and trial compensated (+,+,-,-), magnitudes 1 mu_B per Ni (initialization, not a claim of AFM ground order); shared-mu FD theta=0.001 Ha; U_eff=0/0.2 Ha with exact projectors matched before nonzero-U comparison. Full half-shifted 2^3/4^3/6^3, 64/68/72 bands per channel, Ecut 60/80/100 Ha, density 4*Ecut. Increase capacity only via a revised frozen manifest if qualified valence requires it. |
| R4 derivatives / MF13,14 | R1–R3 center geometries, displace first atom x by 0.03 Bohr to avoid symmetry-zero-only tests; use converged choices from their scans, then the exact displacement/strain steps in the table. Preserve N, theta, fixed-M constraint if present, k fractional coordinates and branch at +/- steps. Request forces and full stress, record QE sign/units before conversion. |

Every dispatched bundle must contain exact BAANDI input, exact `qe.in`, resolved
FFT dimensions/cutoffs, occupation/smearing and spin conventions, symmetry flags
(full meshes initially), Hubbard projector definitions if active, all displaced/
strained inputs, SCF/eigen thresholds, executable version/hash, build/library and
thread/MPI environment, UTC, exact launcher command, raw stdout/stderr and checksums.
Initial QE electronic threshold target is 1e-12 Ry; tighten until finite-difference
noise is below one third of the cap. BAANDI density/eigen/energy thresholds are
independently tightened and all certificate values retained. QE's total-energy
threshold alone is not equivalent to BAANDI's density/eigen certificate.
Use the same pseudopotential bytes and XC; any explicit mismatch override remains
a consistency-only case. No automatic assumption that two smearing-width names
or projector labels have the same definition. Full original output is required;
a relayed summary is labeled as such per [evidence-map](evidence-map.md).

These recipes are requests for future independent calculations, not matches to
published material accuracy. Missing hashes, unqualified datasets, unmatched U
projectors or unresolved branches block the respective real-system gate. A
literature value without identical model/input provenance cannot replace it.
