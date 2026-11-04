# DFT capability matrix

| Scope | Status | Companion |
| --- | --- | --- |
| 77 capability rows; eleven established codes and rust-dft; materials-driven priorities | Planning survey; documentation evidence is not a numerical certification | [Master plan](../roadmap/master-plan-proposal.md) |

This document is being assembled for P8F-MASTER-FEATURE-MATRIX. It distinguishes documented support, restricted or external support, and unverified claims. No external implementation source is used.

## Reading the matrix

Survey date: 2026-11-04; repository snapshot: `318b1c2`. Columns use
Y = documented support, P = partial/restricted support, E = support through an
external driver/tool (a partial capability), N = incompatible with the stated
method or explicitly excluded, and **? = weak documentation evidence; no support
verdict**. An unanswered cell is deliberately not called N. Y does not establish
that arbitrary combinations work, that derivatives exist, or that a particular
release/build contains the feature. Table headings link the public documentation
used; references give titles and access year. Restrictions in cells take precedence
over broad feature lists. Code-family basis labels are not interchangeable:
GPAW PAW reconstruction and CP2K GAPW are not LAPW, and a density FFT grid is not
a real-space orbital discretization.

Engine entries are strictly Y/N for the stated scope; a restricted Y is not a
claim of general support. E1–E8 below identify inspected local evidence. Priorities
mean must = required for target-material credibility, should = next useful
extension, could = research-dependent, not planned = separate engine strategy.
They are proposals, not permission to implement out-of-scope physics.

### Basis, Hamiltonian and pseudopotentials

| Capability | VASP [V] | QE [Q] | ABINIT [A] | CASTEP [C] | GPAW [G] | CP2K [P] | WIEN2k [W] | Elk [E] | FHI-aims [F] | SIESTA [S] | Octopus [O] | Engine today | Priority and reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Plane-wave orbital basis | Y | Y | Y | Y | Y | P: optional | P: LAPW | P: LAPW | N: NAO | N: NAO | N: real-space | Y: E1 | must: existing foundation |
| Localized orbital basis | N | N | N | N | Y: LCAO | Y: Gaussian | P: local orbitals | P: local orbitals | Y: NAO | Y: NAO | ? | N: E1 | not planned: avoid second solver family |
| All-electron Hamiltonian | N: PAW | N: pseudo | N: PAW | N: pseudo | P: PAW reconstruction | Y: GAPW | Y | Y | Y | N: pseudo | P: selected potentials | N: E2 | not planned: use independent all-electron oracle |
| Real-space orbital discretization | N | N | N | N | Y: FD | ? | N | N | N | N | Y | N: E1 | not planned: retain plane waves |
| Norm-conserving pseudo | ? | Y | Y | Y | ? | Y: GTH | N: all-electron | N: all-electron | N: all-electron | Y | Y | Y: restricted E2 | must: validate 3d/O semicore choices |
| USPP | P: legacy | Y | N: NC/PAW | Y | N: PAW | ? | N | N | N | N | ? | N: E2 rejects | could: augmentation cost |
| PAW ground state | Y | Y | Y | N: properties reconstruction only | Y | N: GAPW distinct | N | N | N | N | ? | N: E2 rejects | should: contingent on NC cost evidence |
| ONCV datasets | ? | Y: NC | Y: NC | ? | ? | ? | N | N | N | P: PSML | P: NC | N: no general qualification E2 | must: qualify specific files, not generator name |
| UPF versions/import | N: POTCAR | Y: v1/v2 | ? | ? | N: setups | ? | N | N | N | P: conversion | P: UPF | Y: 2.0/2.0.0/2.0.1 subset E2 | must: fail closed on unsupported metadata |
| Dataset database/provenance | Y: POTCAR | E: SSSP/PseudoDojo | E: PseudoDojo/JTH | Y: OTFG | Y: setups | Y: basis/potential sets | Y: atomic data | Y: species | Y: species defaults | E: pseudo libraries | E: pseudo libraries | N: general catalog E2 | must: per-file hashes and license review |

Dataset license is independent of code license. POTCAR is restricted; do not
redistribute it. Public availability of SSSP, PseudoDojo, JTH, GPAW setups,
CP2K potential tables, PSML or species defaults is **not** evidence of permission
to redistribute every dataset. This survey did not inspect each dataset license:
all such license verdicts are **weak/unverified**, including ONCV-derived files.
Future manifests must record provider, version, generator, valence/semicore,
relativistic and XC settings, original URL, SHA-256, license text/reference,
redistribution permission and qualification tests. No UPF is added to this repo.

### XC and many-body scope (detailed design belongs to P8D)

| Capability | VASP [V] | QE [Q] | ABINIT [A] | CASTEP [C] | GPAW [G] | CP2K [P] | WIEN2k [W] | Elk [E] | FHI-aims [F] | SIESTA [S] | Octopus [O] | Engine today | Priority and reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LDA | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y: PZ E3 | must: regression anchor |
| GGA | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y: PBE E3 | must: structural baseline |
| Meta-GGA | Y | P | Y | Y | P | Y | Y | Y | Y | ? | P | N: E3 | should: improved energetics |
| Hybrid/exact exchange | Y | Y | Y | Y | P | Y | P | P | Y | ? | P | N: E3 | should: gaps and defect localization |
| vdW/dispersion | Y | Y | Y | Y | Y | Y | ? | ? | Y | Y | P | N: E3 | should: layered cathodes |
| DFT+U | Y | Y | Y | Y | Y | Y | Y | Y | ? | Y | ? | N: E3 | must: 3d magnetic manifold |
| RPA/ACFDT | Y | P: suite | Y | ? | Y | P | ? | P | P: scope limited | ? | ? | N: E3 | could: expensive benchmark tier |
| GW | Y | E: Yambo | Y | ? | Y | P | E: ecosystem | ? | P: system dependent | E: ecosystem | ? | N: E3 | could: quasiparticle benchmarks |
| BSE | Y | E: Yambo | Y | ? | Y | ? | ? | Y | ? | ? | ? | N: E3 | could: excited-state scope |
| DMFT interface | E | E | Y | ? | ? | ? | E | E | ? | ? | ? | N: E3 | could: correlated research workflows |

The broad many-body labels are not force/stress guarantees. P8D owns functional
formulation, library choices, Hubbard projectors, cRPA and many-body validation;
this matrix only places those capabilities in the dependency order.

### Spin, occupations and fields

| Capability | VASP [V] | QE [Q] | ABINIT [A] | CASTEP [C] | GPAW [G] | CP2K [P] | WIEN2k [W] | Elk [E] | FHI-aims [F] | SIESTA [S] | Octopus [O] | Engine today | Priority and reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Collinear spin | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | N: E4 | must: transition metals |
| Noncollinear spin | Y | Y | Y | P | Y | ? | ? | Y | ? | Y | Y | N: E4 | should: frustrated magnetism |
| Spin-orbit coupling | Y | Y | Y | Y | P | ? | Y | Y | Y | Y | Y | N: E2/E4 | should: anisotropy/heavy elements |
| Constrained magnetism | Y | Y | Y | ? | P | P | P | Y | ? | ? | ? | N: E4 | should: metastable states |
| Spin spirals | Y | ? | Y | ? | P | ? | ? | Y | ? | ? | ? | N: E4 | could: noncollinear phases |
| Magnetic anisotropy | Y | P: SOC differences | P | P | P | ? | Y | Y | P | Y | ? | N: E4 | should: magnetic applications |
| Fixed integer occupations | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y: 2/0 insulator E4 | must: preserve baseline |
| Fractional occupations | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | N: E4 rejects | must: alloys/metals |
| Smearing | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | N: E4 | must: metallic SCF |
| Tetrahedron integration | Y | Y | Y | ? | Y | ? | Y | Y | ? | ? | ? | N: E4 | should: accurate DOS/integration |
| Finite electronic temperature | Y | Y | Y | Y | Y | Y | ? | Y | Y | Y | Y | N: E4 | must: consistent free energy |
| Berry-phase polarization | Y | Y | Y | Y | Y | Y | Y | ? | ? | Y | P | N: E5 | should: insulating polarization |
| Slab dipole correction | Y | Y | Y | ? | Y | P | ? | ? | Y | Y | P | N: E5 | should: surface/work-function studies |
| Finite electric field | Y | Y | Y | Y | Y | Y | ? | Y | P | Y | Y | N: E5 | should: gated coupling physics |
| Magnetic field response | P | P | P | P | ? | P | P | Y | ? | ? | Y | N: E4/E5 | could: distinguish Zeeman/orbital terms |
| Implicit solvent | E: VASPsol | E: Environ | ? | ? | Y | P | ? | ? | P | ? | ? | N: E5 | could: electrochemical extension |
| Constant-potential electrochemistry | E: extensions | P: ESM/FCP/GC-SCF | ? | ? | P: SJM | ? | ? | ? | ? | ? | ? | N: E4/E5 | could: needs grand-potential model |
| Charged cells/defect corrections | P: corrections | P: corrections | P: corrections | P | E | P | ? | ? | P | P | P | N: qualified workflow E5 | must: defect validation before production |

Finite-temperature here concerns electrons; it is not NVT ionic dynamics.
A charged-cell calculation with a neutralizing background is not by itself a
validated defect formation energy. Constant electron number in an external
potential is not constant-potential electrochemistry. The first coupling gate
must retain the neutral, fixed-N restriction until separately approved.
### Derivatives, dynamics and response

| Capability | VASP [V] | QE [Q] | ABINIT [A] | CASTEP [C] | GPAW [G] | CP2K [P] | WIEN2k [W] | Elk [E] | FHI-aims [F] | SIESTA [S] | Octopus [O] | Engine today | Priority and reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Ionic forces | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | N: E5 | must: structure/MLIP prerequisite |
| Stress | Y | Y | Y | Y | P: mode dependent | Y | ? | Y | Y | Y | ? | N: E5 | must: cell and voltage studies |
| Ionic relaxation/constraints | Y | Y | Y | Y | E: ASE | Y | Y | Y | Y | Y | Y | N: E5 | must: stable geometries |
| Cell relaxation/constraints | Y | Y | Y | Y | E: ASE | Y | P | Y | Y | Y | ? | N: E5 | must: phase energetics |
| Ab initio MD, NVE | Y | Y | Y | Y | E: ASE | Y | ? | ? | Y | Y | Y | N: E5 | should: force-drift validation |
| NVT dynamics | Y | Y | Y | Y | E: ASE | Y | ? | ? | Y | Y | ? | N: E5 | should: thermal sampling |
| NPT dynamics | Y | P: driver | Y | Y | E: ASE | Y | ? | ? | ? | Y | ? | N: E5 | could: needs stress certification |
| ML force fields/active learning | Y: MLFF | E | E | ? | E: ASE | P: integrations | ? | ? | E | E | ? | N: E5/E7 | should: data first, external trainer |
| NEB | Y | Y: neb.x | Y | Y | E: ASE | Y | E | ? | E: ASE | E | ? | N: E5 | should: ion migration |
| Dimer transition-state search | Y | ? | ? | ? | E: ASE | Y | ? | ? | E | ? | ? | N: E5 | could: unknown saddle pathways |
| DFPT phonons | Y | Y | Y | Y | ? | ? | ? | Y | ? | ? | ? | N: E5 | could: finite displacement first |
| Finite-displacement phonons | Y | E | E | Y | E: ASE | P | E: PHONON | P | E | Y | ? | N: E5 | should: dynamical stability |
| Quasi-harmonic free energy | E: phonopy | E | Y: anaddb | P | E | E | E | ? | E | E | ? | N: E5 | could: volume/temperature phase maps |
| Elastic constants | P | P | Y | P | E | P | P | P | P | P | ? | N: E5 | should: strain energetics |
| Static dielectric response | Y | Y | Y | Y | Y | P | P | Y | P | Y | P | N: E5 | should: charged-defect correction |
| Born effective charges | Y | Y | Y | Y | P | P | ? | ? | ? | P | ? | N: E5 | should: polar lattice response |
| Piezoelectricity | Y | P | Y | P | ? | ? | ? | ? | ? | ? | ? | N: E5 | could: coupled electromechanics |
| Independent-particle optics | Y | Y | Y | Y | Y | P | Y | Y | P | P | Y | N: E5 | should: baseline spectra |
| TDDFT/linear optical response | P | Y: turbo suite | P | ? | Y | Y | ? | Y | P | Y | Y | N: E5 | could: separate response operators |

### Analysis, execution and ecosystem

| Capability | VASP [V] | QE [Q] | ABINIT [A] | CASTEP [C] | GPAW [G] | CP2K [P] | WIEN2k [W] | Elk [E] | FHI-aims [F] | SIESTA [S] | Octopus [O] | Engine today | Priority and reason |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| DOS/PDOS | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | N: analysis workflow E7 | must: electronic diagnostics |
| Bands/path sampling | Y | Y | Y | Y | Y | P | Y | Y | Y | Y | Y | N: dedicated workflow E7 | must: characterize states |
| Fermi surface | P | E | P | P | E | ? | Y | Y | ? | E | ? | N: E4/E7 | should: alloys and metals |
| Wannier functions | E: Wannier90 | E: Wannier90 | E: Wannier90 | ? | E | P: localized orbitals | E: wien2wannier | E | ? | E: Wannier90 | ? | N: E7 | should: correlated subspaces |
| Charge analysis | E: Bader | P | P | Y | E | Y | Y | Y | Y | Y | P | N: beyond total charge E7 | should: oxidation-state diagnostics |
| ELF | Y | Y | Y | Y | Y | Y | Y | Y | ? | P | Y | N: E7 | could: bonding diagnostics |
| Work function | P: slab potential | P | P | P | P | P | P | P | P | P | ? | N: calibrated vacuum E5/E7 | should: surface electrochemistry |
| STM simulation | Y | Y | P | ? | E: ASE | P | Y | Y | ? | P | ? | N: E7 | could: experimental comparison |
| Space-group reduction/density symmetrization | Y | Y | Y | Y | Y | P | Y | Y | Y | P | P | N: E6 | should: efficient large k grids |
| Weighted k-point meshes | Y | Y | Y | Y | Y | P | Y | Y | Y | Y | Y | Y: MP/explicit E6 | must: periodic materials |
| MPI | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | N: E8 | should: after measured CPU limits |
| Threads | Y | Y | Y | Y | Y | Y | P | Y | P | Y | Y | Y: scoped k/band E8 | must: retain deterministic folds |
| GPU | P | P | P | ? | P | P | ? | ? | P | P | P | N: E8 | should: measured backend project |
| ASE/workflow integration | E | E | E | E | Y: ASE | E | ? | E | E | E | E | N: official adapter E7 | should: P8E owns implementation |
| AiiDA/pymatgen integration | E | E | E | ? | ? | E | ? | ? | ? | E | ? | N: E7 | should: external provenance workflows |
| Machine-readable input/output | Y | Y: XML | Y: NetCDF | Y | Y: Python/GPW | Y | P: text | P: text | P: text | Y: FDF/NetCDF | Y | Y: TOML/JSON/binary E7 | must: explicit units/schema |
| Mixing/preconditioning/restarts | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y | Y: Pulay/Davidson; no reload E4/E7 | must: robust SCF, P8E |
| ACE exchange acceleration | P: VASP6/build | Y | ? | ? | ? | ? | ? | ? | ? | ? | ? | N: E3/E8 | could: only after exact exchange oracle |
| Block/batched numerical operations | Y | Y | Y | Y | Y | Y | ? | ? | Y | Y | P | Y: apply_block/GEMM E8 | must: profile each change |
| Mixed precision | P | P | ? | ? | ? | P | ? | ? | ? | ? | ? | N: E8 | could: require residual correction |

E entries without a named adapter are **weak evidence of ecosystem availability**,
not a verified integration or compatibility promise. For NPT, transition-state,
phonon and analysis workflows, external orchestration can only use a calculator
that supplies the required observables; a generic ASE interface does not create
missing stress or forces. GPU, mixed precision, ACE, DMFT and lesser-used response
cells require release/build-specific checks before they can inform procurement or
a development commitment. Broad Y/P cells inherited from an overview have
**medium evidence**, not the strength of a tested manual example.

### Explicit weak-evidence follow-up register

The following positive-looking cells are **tentative planning leads**, not
fully verified support. The cited overview does not establish their exact scope;
read the appropriate release manual/example before promoting them to a firm Y.
This qualification overrides the shorthand cell value.

| Cells/area | Evidence gap and required follow-up |
| --- | --- |
| GPU, mixed precision and block/batched rows | Hardware/build/kernel coverage is largely unspecified; all mixed-precision positives and batching claims outside the inspected engine are weak. Verify actual precision mode and supported operators. |
| ASE and AiiDA/pymatgen rows, except GPAW's explicit ASE integration | Named frameworks are ecosystem leads, not checked plugin installations. Verify maintained plugin, supported version and result/error semantics for each code. |
| DMFT interfaces, except QE input controls and ABINIT feature inventory | Named capability is not an audited export/import contract; projector normalization, double counting and supported solvers remain unverified. |
| External NEB/dimer, finite-displacement, QHA and ML workflows | Driver availability and force/stress compatibility for each code were not exercised. E means partial/external, not a native executable or all-Hamiltonian support. |
| Constant-potential extensions outside QE controls and GPAW SJM | Extension version, thermodynamic functional and force consistency were not verified. Do not equate an implicit-solvent flag with constant potential. |
| Dataset formats outside explicit UPF/PSML/setup documentation | ONCV labels and conversions need file-level tests. Dataset license/provenance cells are all unverified for redistribution. |
| Specialized response/fields and derived-property P cells | Overview-level possibilities do not establish analytic derivatives or a packaged workflow; strain differences, slab references and experimental response need individual manual checks. |
| Octopus historical version 12 cells | Evidence is release-specific and weaker for current releases; do not infer present availability or removal. |

## Verified engine evidence map

The following files were inspected locally. Negative entries mean no implemented
end-to-end feature at this snapshot, not a theorem that extension is impossible.
All paths are relative to this document. Existing numerical results are not rerun
external references by this survey.

| ID | Inspected evidence | What it establishes and what it does not |
| --- | --- | --- |
| E1 | [basis](../../crates/dft-basis/src/lib.rs), [Hamiltonian](../../crates/dft-physics/src/hamiltonian.rs), [SCF](../../crates/dft-engine/src/scf.rs) | Complex plane-wave Bloch states; matrix-free kinetic/local/NC nonlocal apply. One real scalar density grid; no alternative orbital basis implementation. |
| E2 | [UPF parser](../../crates/dft-upf/src/parser.rs), [UPF support](../upf-support.md) | NC subset only; UPF versions explicitly enumerated. PAW, USPP, SOC, bare Coulomb, NLCC and all-electron-wavefunction flags rejected. A specific ONCV file might fit the subset; generator-wide compatibility is not established. |
| E3 | [Functional](../../crates/dft-engine/src/functional.rs), [GGA](../../crates/dft-physics/src/gga.rs), [ADR 0018](../decisions/0018-pbe-via-xcx.md) | LdaPz and Pbe exist; xcx-backed PBE has mathematical/synthetic evidence. ADR still defers the real-Si PBE/QE comparison. No other functional selector or orbital-dependent XC operator. |
| E4 | [SCF](../../crates/dft-engine/src/scf.rs), [mixing](../../crates/dft-solvers/src/mixing.rs) | Fixed 2/0 occupations; PartialOccupations explicitly rejected; scalar density and spin-unpolarized problem. Occupied residual, density and fresh fixed-point checks; Pulay/linear mixing. No entropy or chemical-potential solve. |
| E5 | [energy ledger](../../crates/dft-engine/src/energy.rs), [physics module list](../../crates/dft-physics/src/lib.rs), [SCF](../../crates/dft-engine/src/scf.rs), [coupling ADR](../decisions/0017-coupling-commitments.md) | Direct/band/Harris-Foulkes energy audit, not derivative or response infrastructure. No force/stress output or field/coupled functional implementation. Coupling ADR is proposed, not delivered. |
| E6 | [k points](../../crates/dft-geometry/src/kpoints.rs), [ADR 0011](../decisions/0011-kpoints.md) | MP and explicit weighted k points, time reversal only. Space-group and magnetic reduction not implemented. |
| E7 | [CLI config](../../crates/dft-cli/src/config.rs), [CLI run](../../crates/dft-cli/src/run.rs), [engine result](../../crates/dft-engine/src/scf.rs), [ESLE ADR](../decisions/0016-esle-layering.md) | Structured reports/density export and provenance; eigenpairs at SCF k points are not a band/DOS workflow. No supported checkpoint reload, ElectronicState implementation, or production workflow adapter. |
| E8 | [numerics](../../crates/dft-numerics/src/lib.rs), [linear algebra](../../crates/dft-numerics/src/linalg.rs), [solver](../../crates/dft-solvers/src/lib.rs), [SCF](../../crates/dft-engine/src/scf.rs), [Hamiltonian](../../crates/dft-physics/src/hamiltonian.rs), [ADR 0019](../decisions/0019-optional-mkl-backend.md) | Reusable FFT plans, block applies, GEMM and scoped k/band workers exist. Current numerical source has no MKL selection; accepted ADR 0019 is not proof of an installed backend. No GPU/MPI/mixed-precision implementation. |

The older [architecture overview](../architecture.md) still describes LDA-only
and serial execution in places. Actual source and implementation results in ADRs
0018 and 0001 override those stale statements for this survey. Likewise proposed
ESLE/coupling ADRs do not count as implemented features. P8D
`docs/research/xc-and-manybody-roadmap.md` and P8E
`docs/research/feature-gap-survey.md` were absent at inspection; reconcile with
their owners before issuing overlapping tickets.

## Target-material shortlist

| Target | Minimum credible bundle | Main dependency/risk |
| --- | --- | --- |
| Transition-metal oxides | Qualified O/3d NC datasets, PBE, spin, DFT+U, k/cutoff convergence | NLCC and semicore support may block dataset choice; PAW is an alternative project, not an automatic requirement |
| Battery cathodes | Oxide bundle, forces/stress, ionic/cell relaxation, consistent phase references | Compare magnetic branches and identical U/XC policies; add dispersion for layered phases and NEB for migration |
| Defects | Supercells, relaxation, multiple starts, charge/reference policy, dielectric corrections | Neutral defects first; charged formation energies require finite-size and potential-alignment convergence |
| Alloys | Smearing/free energies, spin when relevant, k sampling, structural ensembles | Disorder/SQS orchestration external; composition/reference consistency precedes phase ranking |
| Magnetic materials | Spin state identity, occupation control, DFT+U, constrained starts | Noncollinear/SOC/anisotropy follow validated spinors and relativistic projectors |
| Multiscale and MLIP | Typed observables, reproducible state, complete ledger, certified derivatives | Energy-only datasets first; never fabricate forces from unavailable observables |

Detailed prerequisites, effort estimates, user decisions and the ordered queue
are in the [master plan](../roadmap/master-plan-proposal.md).

## Limits of this document

This is a bounded planning survey, not an exhaustive test of eleven products.
Only public documentation was read; no QE or other external implementation source
was read, cloned or copied. No external executable was run. Web pages mix release
lines and some URLs failed; support may depend on compile options or an extension.
Question marks, unnamed ecosystem adapters and database licenses remain explicit
verification tasks. In particular do not infer that the all-electron codes lack
features merely because their public overview is less detailed.

Priorities, estimates and architectural recommendations are engineering opinion.
Cross-products (hybrid + SOC + stress, U + phonons, PAW + response, GPU + a given
functional) remain unverified even where individual features are Y. NC accuracy
for a specific oxide must be measured; Si agreement cannot establish it. Code
license, dataset license, redistributability and physical transferability are
separate questions. No vendor ranking, performance comparison or new QE numerical
reference is asserted here.

## References

All online references are public documentation, accessed 2026-11-04 (2026 survey
year; undated living pages are not assigned a fictitious publication year).
These are feature inventories, not source-code references. Column letters refer
to the following bundles; feature-specific exceptions and weak cells are marked
above. The linked inventories do not certify every feature combination.

[V]: https://vasp.at/info/about/
[Q]: https://www.quantum-espresso.org/Doc/pw_user_guide/node3.html
[A]: https://docs.abinit.org/topics/features/
[C]: https://www.castep.org/features/capabilities
[G]: https://gpaw.readthedocs.io/documentation/documentation.html
[P]: https://www.cp2k.org/features
[W]: https://susi.theochem.tuwien.ac.at/features/index.html
[E]: https://elk.sourceforge.io/
[F]: https://fhi-aims.org/learn-more
[S]: https://www.siesta-project.org/siesta/About/overview.html
[O]: https://www.octopus-code.org/documentation/12/octopus_features/

- V — VASP, *About / features*, living documentation (2026 access); [Many-body perturbation theory](https://vasp.at/wiki/index.php/Category:Many-body_perturbation_theory), [Spin spirals](https://vasp.at/wiki/Spin_spirals), [LFOCKACE](https://vasp.at/wiki/LFOCKACE), and [Machine learning force fields](https://vasp.at/wiki/Machine_learning_force_field_calculations). General methods/category URLs failed; the official About page supplied the broad inventory. ACE is build/algorithm dependent in VASP 6.
- Q — Quantum ESPRESSO, *What can PWscf do*, living user guide (2026 access); [pw.x input description](https://www.quantum-espresso.org/Doc/INPUT_PW.html), retrieved page labeled 7.5, including occupations, ACE, DMFT, constrained moments and electrochemical controls. Suite/external capabilities must not be confused with pw.x itself, and manual version is not the repository reference-run version.
- A — ABINIT, *Features*, living documentation (2026 access); [Many-Body Perturbation Theory](https://docs.abinit.org/theory/mbt/).
- C — CASTEP, *Capabilities*, living documentation (2026 access); [Vibrational Spectroscopy](https://www.castep.org/features/capabilities/vibrational-spectroscopy).
- G — GPAW, *Documentation*, [Introduction](https://gpaw.readthedocs.io/documentation/basic.html), [Advanced topics](https://gpaw.readthedocs.io/documentation/advanced.html) and [Solvated Jellium](https://gpaw.readthedocs.io/documentation/sjm/sjm.html), living documentation (2026 access); mode restrictions matter.
- P — CP2K, *Features*, living documentation (2026 access); GPW/GAPW distinctions matter.
- W — WIEN2k, *Features and Calculated Properties*, official TU Wien host, living documentation (2026 access); primary homepage timed out.
- E — Elk, *Features*, living project documentation (2026 access); no code browser followed.
- F — FHI-aims, *Learn more / Quick overview*, living documentation (2026 access); distinguishes molecular and periodic many-body support.
- S — SIESTA, *About* and [Guide to SIESTA versions](https://siesta-project.org/siesta/Developers/Guide-to-Siesta-versions.html), living documentation (2026 access); [5.4 manual, spin-orbit coupling](https://docs.siesta-project.org/projects/siesta/en/5.4/reference/siesta.html).
- O — Octopus, *Octopus Features*, version 12 documentation (2026 access; historical release coverage). The attempted main/capabilities URL was unavailable; current support must be rechecked against the maintained manual.
- Local sources: [validation](../validation.md), [P7 plan](../p7-plan.md), [performance rules](../performance/README.md), and E1–E8 above.
