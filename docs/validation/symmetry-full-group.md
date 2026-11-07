# Full crystal-group density experiment (P34A)

## Conclusion and frozen criterion

**The hypothesis is supported.** On shifted Si 4x4x4, averaging the density over
all 48 crystal operations produces a crystal-minus-off energy difference of
-1.5538989465824216e-5 Ha with tightened convergence. The specified relayed
QE 7.6 difference is -1.554e-5 Ha: absolute discrepancy 1.010534175784e-9 Ha,
relative discrepancy **0.00650279%**, with the same negative sign.

Before measurement, the classification was fixed at relative discrepancy <=5%
for supported, <=20% for partly supported, otherwise rejected (a wrong sign also
rejects). The reference is the ticket's rounded difference, not a newly executed
QE calculation. Provenance is [the relayed record](../../tests/reference/perf/dev08_materials_qe_relayed.md),
section (5) and its correction. No QE source was accessed. This supports the
proposed explanation; it does not prove the reference code's internal path.

The effect is consistent with different effective k sampling: projecting over
the crystal group corresponds to the union of rotated meshes at an invariant
potential. This is not a claim that either sampling is more correct. No constant
energy offset is fitted or used. No default or k-point reduction rule changed.

## Option and implementation

```toml
[symmetry]
mode = "auto"
density_group = "crystal" # research option; default is "mesh"
```

Both choices use exactly the same mesh-preserving reduction, representatives,
weights and order. Shifted Si uses 12 spatial mesh operations, 64 full points,
and 10 representatives. `mesh` keeps the existing density projector and its
summation order. `crystal` supplies all 48 detected spatial operations to that
same prepared Fourier projector, including its input-density invariance check.
No change to the physics projector was needed; no allocation or I/O was added
to its application. Time-reversal entries are not density operations.

Grid preparation uses the selected density group. Every operation must satisfy
integer `N_i*tau_i` (allowance 1e-8) and integer `N_i*R_ij/N_j`. An explicit
nonpermutable grid fails, including diamond at 25^3. Automatic grids use the
existing compatible-grid search, without dropping any detected operation; failure
to find a compatible grid is an error. The engine independently checks the whole
group. All prediction runs explicitly use 32^3.

Results `symmetry.density_group`, `density_operations_count` and
`density_operations` identify the policy and every rotation, translation and atom
permutation actually selected. Existing `operations` and `orbit_map` continue to
identify the k-point reduction. Preflight reports both operation counts.

`crystal` with off/missing auto, collinear spin, force/stress requests or relaxation
rejects. The engine rejects missing symmetry activation, spin construction and
analytic derivatives independently. The projection policy travels with engine
results so supplying mesh options for a crystal result cannot enable derivatives.
Stress uses the checked force boundary.

Exact option comparison includes `ScfOptions.crystal_density_group`; wire omission
means false for old checkpoints. **Symmetry checkpoints remain unsupported**:
the v1 envelope lacks an operation record. Capture/import rejects instead of
claiming exact continuation. CLI checkpoint context now explicitly rejects auto
(the old context could bypass activation). Tests cover option mismatch and wire
roundtrip/legacy omission, not a newly supported crystal restart.

## Prediction measurements

Mac Apple Silicon, pure Rust portable release, T=1, strict profile, Si.pz-vbc.UPF,
30 Ry and explicit 32^3 grid in every run. Input generation uses
`scripts/bench_scf.py::case_input(k444s)`; its structure is the same primitive cell
as `benchmarks/cases/si-k444s.toml` (the latter is a flat benchmark specification,
not a CLI TOML input). Defaults: density 1e-7 Bohr^-3/2, energy 1e-8 Ha and engine
eigen tolerance 1e-8 Ha. Tightened: density 1e-11 Bohr^-3/2, energy 1e-12 Ha,
eigen tolerance 1e-11 Ha. No acceptance tolerance was relaxed.

| Mesh / convergence | E(off), Ha | E(auto, mesh), Ha | E(auto, crystal), Ha |
|---|---:|---:|---:|
| Shifted / default | -7.925976134331574 | -7.925976134331593 | -7.9259916733210085 |
| Shifted / tight | -7.9259761343315365 | -7.925976134331535 | -7.925991673321002 |
| Unshifted / default | -7.919059567875023 | -7.919059567874930 | -7.919059567874930 |
| Unshifted / tight | -7.9190595678749816 | -7.919059567875055 | -7.919059567875055 |

| Mesh / convergence | mesh - off, Ha | crystal - off, Ha | crystal - mesh, Ha |
|---|---:|---:|---:|
| Shifted / default | -1.9539925233402755e-14 | -1.553898943473797e-5 | -1.5538989415198046e-5 |
| Shifted / tight | 1.7763568394002505e-15 | -1.5538989465824216e-5 | -1.5538989467600572e-5 |
| Unshifted / default | 9.325873406851315e-14 | 9.325873406851315e-14 | 0 |
| Unshifted / tight | -7.37188088351104e-14 | -7.37188088351104e-14 | 0 |

On unshifted 4x4x4 all 48 operations preserve the mesh; mesh and crystal energy
and density are bit-identical at both tolerances (maximum density difference 0).
The shifted prediction changes by only 3.10862447e-14 Ha upon tightening.

## Independent mechanism check

The evidence analyzer directly permutes each real-space grid point by
`x -> R*x + tau`. It does not call the Fourier projector, its phases or its
invariance helper. It checks all 48 operations separately on saved CLI binary
densities, and records every operation's maximum error in `analysis.json`.
Units below are electrons/Bohr^3; these are maximum pointwise errors.

| Shifted density | Default convergence | Tight convergence |
|---|---:|---:|
| Off: error for each of the 36 non-mesh operations, min–max | 0.002126881950608006–0.0021268819562402785 | 0.0021268811697882917–0.002126881169820849 |
| Mesh: error for each non-mesh operation, min–max | 0.002126881952876830–0.002126881952876844 | 0.0021268811698035572–0.0021268811698035572 |
| Crystal: maximum over all 48 operations | 6.938893903907228e-17 | 6.245004513516506e-17 |
| Maximum crystal-versus-mesh density difference | 0.001589446794227789 | 0.0015894455123635631 |

All 36 non-mesh operations fail invariance for the full unsymmetrized mesh density;
all 48 pass the crystal-density 1e-13 gate. The unshifted crystal errors are
6.938893903907228e-17 (default) and 5.551115123125783e-17 (tight).
An independent Rust release test also checks real-space permutations on engine
results, the unshifted equality, and each of the 36 failures. Its mixer defaults
differ from the CLI Pulay fixture, so its numbers are separately retained.

## Reproduction and gates

All raw inputs, outputs, logs, scripts and metadata are under `P34A-evidence/`
in this worktree, not under target and not committed UPF files. Preserve this
directory before worktree cleanup. Starting commit:
`5951f23b05f5feccdd080310a0f2efcc23feb0c9`.

UPF SHA-256:
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
Baseline executable SHA-256:
`e083ce669f0376b5d57cf8a20615979349b096905eb073b2a9dc6c5c2e67baf9`.
Prediction executable (preserved as `prediction-cli`) SHA-256:
`9d35279487aadda65856da3aa65c9c20e75680b54da4c0ec01298502c653e39d`.
Final executable hashes and command exits are in `P34A-evidence/manifest.json`.

```sh
cargo build --release -p dft-cli
python3 P34A-evidence/identity.py baseline  # requires saved pre-edit baseline-cli
python3 P34A-evidence/identity.py candidate
python3 P34A-evidence/prediction.py
python3 P34A-evidence/analyze.py
cargo test -p dft-engine --test symmetry_full_group
DFT_SI_UPF="$PWD/../pseudo/Si.pz-vbc.UPF" cargo test --release -p dft-engine --test symmetry_full_group -- --ignored --nocapture
cargo test -p dft-cli density_group
scripts/check.sh
```

The identity driver invokes `bench_scf.py --require-identical` in 20 candidate
commands: five cases (gamma, k222s, k222, k444s, si8), T=1,2,4,8,12,
symmetry off/auto and strict/standard: **100/100 identical energy bits and density
hashes**, with identical SCF, Davidson and main/certificate H-work counts.
Historical Si8 absolute-anchor bypass is explicit; baseline comparisons remain
active. No additional target directory was created.

Final verification: `scripts/check.sh` exited **0**, with build, workspace tests,
format and clippy all PASS (`verified-check.log`). The final release integration
command passed both tests (`verified-oracles.log`); CLI policy/metadata and wire
unit tests passed. All 20 final candidate benchmark commands exited 0, and the
12 final prediction CLI runs exited 0. The independent analyzer exited 0.
The final result-policy guards did not change any reported prediction number.

Initial validation mistakes were corrected rather than hidden: a misspelled
CLI eigen-tolerance key rejected six initial tight inputs (exit 1); use
`eigen_tolerance`. Test-only `unwrap_err` on a non-Debug success type caused a
compile failure; a config test initially expected validation later than parsing.
A provenance guard initially addressed the finite-T step as a result, caught by
the compiler and corrected. Final gate logs supersede these intermediate failures.

## Changed-file scope and limits

Engine files: `scf.rs` (selection and result policy), `session.rs` (result policy
and checkpoint boundary), `checkpoint.rs` (exact comparison), `collinear.rs`
(spin rejection), `forces.rs` (force/stress rejection). The last two are solely
fail-closed plumbing for the explicitly excluded modes. CLI files: `config.rs`,
`run.rs`, `preflight.rs`, `checkpoint_io.rs`. New integration test:
`crates/dft-engine/tests/symmetry_full_group.rs`. Documentation: this study and
append-only performance ledger entry. The existing physics projector is unchanged.

No new QE run, source inspection, dev08 timing, general performance acceptance,
cutoff convergence, material generality or preferred-sampling claim is made.
This research option does not validate spin or derivatives, or enable symmetry
restart. Artifacts and README are intentionally untouched per this ticket.
Changes remain uncommitted for review.
