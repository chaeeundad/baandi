# Reference compatibility mode: user guide and evidence

## Configuration

This supplements the [CLI configuration guide](../cli.md). For an existing
unpolarized PBE energy/density input, add:

```toml
[compat]
mode = "reference"
```

It resolves `xc.pw92_prefactor="original"`, the reference-style local mesh
`pseudo.local_mesh="qe_like"`, and
`pseudo.radial_even_rule="simpson_quadratic_end"`. If the input already has
`[symmetry] mode="auto"`, it also resolves `density_group="crystal"`.
With symmetry off or absent, `density_group` stays `mesh`. The switch does not
activate symmetry or change convergence settings, the functional, or UPFs.

Omit the table, or use `mode="none"`, to retain existing defaults. Explicit
keys agreeing with the mode are allowed; a disagreement rejects and names both
settings. LDA, spin-polarized PBE, analytic forces/stress and relaxation reject:
this mode is for unpolarized PBE energies and densities only.

`dft-cli validate input.toml --json` exposes `compat_mode`, the four values and
scope. `dft-cli explain input.toml` and `dft-cli run input.toml` print them.
`results.json` stores `compat_mode="reference"` and a `compat` record, also in
`scf`; other result profiles retain the shared `scf` record. Explain of scalar
results retains that record. Defaults omit the new fields.

Exact restart requires the same mode label, even if manually selected options
are identical; use explicit WarmStart for a label change. Old checkpoint fields
that are absent mean `none`. Existing checkpoint restrictions (including
symmetry and NLCC support) still apply. See [ADR 0038](../decisions/0038-reference-compat-mode.md),
[ADR 0036](../decisions/0036-pw92-prefactor-convention.md), and
[ADR 0037](../decisions/0037-radial-convention-options.md).

## Prior measured evidence

These are copied from [the radial score](radial-convention-score.md), using the
original PW92 prefactor, reference-style mesh and quadratic endpoint. Reference
totals are the real QE 7.6 runs in the [dev08 relayed record](../../tests/reference/perf/dev08_materials_qe_relayed.md),
not new reference runs. That record contains input hashes, binary/version,
commands and all three UPF SHA-256 identities. QE prints eight decimal places
in Ry (approximately +/-5e-9 Ry rounding uncertainty). The six material inputs
have symmetry off. A case label's number is its cutoff in Ha, not Ry.

| Case | Prior BAANDI total (Ry) | QE total (Ry) | BAANDI - QE (Ry) | SCF iterations |
| --- | ---: | ---: | ---: | ---: |
| si-20 | -16.92463781634268 | -16.92463785 | +3.3657322263e-8 | 15 |
| si-30 | -16.92466015555931 | -16.92466016 | +4.4406931465e-9 | 15 |
| si-40 | -16.92467578740471 | -16.92467579 | +2.5952893168e-9 | 15 |
| mgo-40 | -152.10176985204663 | -152.10176985 | -2.0466188744e-9 | 18 |
| mgo-50 | -152.10194662866155 | -152.10194662 | -8.6615443706e-9 | 19 |
| mgo-60 | -152.10196248930606 | -152.10196248 | -9.3060634754e-9 | 19 |

The [prefactor diagnostic](xc-prefactor-diagnostic.md) separates the original
prefactor effect from the remaining radial difference; no tolerance is tuned.
The si-20 residual is 3.37e-8 Ry: “1e-8 Ry level” is a scale, not a universal
absolute error cap of 1e-8 Ry. Neither convention set is more correct by claim.

From [the full-group study](symmetry-full-group.md), shifted Si k444s with
30 Ry cutoff, explicit 32^3 FFT and tightened convergence:

| Quantity | Value (Ha) |
| --- | ---: |
| Symmetry off total | -7.9259761343315365 |
| Auto, mesh density group total | -7.925976134331535 |
| Auto, crystal density group total | -7.925991673321002 |
| Crystal minus off | -1.5538989465824216e-5 |
| Relayed QE crystal minus off | -1.554e-5 |
| Absolute shift discrepancy | 1.010534175784e-9 |

Relative discrepancy is 0.00650279%. At default tolerances, off/mesh/crystal
are -7.925976134331574, -7.925976134331593, -7.9259916733210085 Ha.
On unshifted tight 4x4x4 they are -7.9190595678749816,
-7.919059567875055, -7.919059567875055 Ha, respectively. The full-group study
used LDA and independently validates that existing convention; it is not a
run of this PBE-only switch. No new QE run or source access occurred here.

## Commands and evidence

Existing study reproduction commands (raw scripts must be retained with their
original evidence directories):

```sh
DFT_NC_DIR=/Users/mino/Develop/rust-dft-wt/pseudo/pseudodojo/ex \
  python3 -I crates/dft-engine/tests/radial_convention_score.py \
  --evidence /tmp/P27D-evidence/materials-final
python3 P34A-evidence/prediction.py
python3 P34A-evidence/analyze.py
```

Reference-mode reproduction from the worktree root:

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
cargo build --offline --release -p dft-cli
DFT_NC_DIR=/Users/mino/Develop/rust-dft-wt/pseudo/pseudodojo/ex \
  cargo test --offline --release -p dft-cli --test reference_mode_equivalence \
  material_mode_matches_manual -- --ignored --nocapture
python3 -I crates/dft-engine/tests/reference_mode_materials.py \
  --binary target/release/dft-cli --evidence P44A-evidence/materials
cargo test --offline -p dft-cli --test reference_mode_config \
  --test reference_mode_restart --test reference_mode_equivalence
scripts/check.sh
```

The material runner uses DFT_NC_DIR, otherwise `../pseudo/pseudodojo/ex`, checks
all three hashes, and skips cleanly when files are missing. Invalid present files
or failed SCFs fail. It runs all six mode cases plus manual si-20/mgo-40 controls,
and requires exact energy, density SHA-256, SCF count and H-work equality.
Inputs, stdout/stderr, result files and exact executable commands are preserved
in `P44A-evidence/materials`; no UPFs are copied. `score.json` records the
executable hash, UPF hashes, results and residuals in Ry.

The small ordinary test uses a synthetic Si-labelled Gaussian NC potential;
it is a plumbing identity test, not a realistic Si accuracy claim. It checks
symmetry off and auto against manual keys. The external material tests are
ignored in debug and documented above. Restart tests change only the label,
check both exact rejection directions, same-label exact success and explicit
WarmStart success. Wire tests check legacy omission, both modes and invalid
values.

Default-path identity uses the saved pre-edit binary and unchanged benchmark:

```sh
python3 -I P44A-evidence/identity.py baseline
python3 -I P44A-evidence/identity.py final
```

The scripts record every expanded `scripts/bench_scf.py --require-identical`
command and exit code. They cover gamma/k222s/k222/k444s/si8, threads
1/2/4/8/12, strict/standard and gamma_real off. With gamma_real on, gamma and
si8 are eligible; multiple-k cases are rejected by the existing implementation,
so those combinations require unchanged rejection rather than an energy value.
No performance improvement is claimed from benchmark timing fields.

## P44A measurements (2026-11-08)

All six mode SCFs converged. The current totals reproduce the prior radial-score
values; no tolerance or parameter was tuned.

| Case | Mode total (Ry) | Reference total (Ry) | Difference (Ry) | SCF iterations |
| --- | ---: | ---: | ---: | ---: |
| si-20 | -16.924637816342678 | -16.92463785 | 3.365732226257023e-08 | 15 |
| si-30 | -16.924660155559305 | -16.92466016 | 4.440693146534613e-09 | 15 |
| si-40 | -16.92467578740471 | -16.92467579 | 2.595289316786875e-09 | 15 |
| mgo-40 | -152.10176985204663 | -152.10176985 | -2.046618874373962e-09 | 18 |
| mgo-50 | -152.10194662866155 | -152.10194662 | -8.661544370625052e-09 | 19 |
| mgo-60 | -152.10196248930606 | -152.10196248 | -9.306063475378323e-09 | 19 |

Manual versus mode on si-20 and mgo-40: exactly equal total energy (0 Ha
difference), density SHA-256, SCF iterations and operator applications.
The default benchmark matrix passed all 70 eligible comparisons with zero energy
change (0 Ha), identical density hashes, SCF iterations and H/Davidson work.
All 30 ineligible gamma_real=on multiple-k combinations produced exactly the
same rejection before and after (benchmark exit 3, expected).

These Mac runs establish identity, not a performance acceptance result or a
new independent external reference measurement. Material executable and UPF
hashes are in `P44A-evidence/materials/score.json`; expanded commands and logs
are retained beside the outputs.

## Relaxation entry guard (approved ownership extension)

The architect extended ownership to `crates/dft-cli/src/relax.rs` for the entry
guard only. Exactly three lines were added after parsing: reference mode returns
the shared scope diagnostic before default relaxation settings or SCF setup.
No other line in that file changed. The command-level regression checks both
strict and standard profiles without a `[relax]` table, exit 1, the mode's own
message, and absence of an output directory. The prepared patches were applied.

A pre-change versus final-executable comparison additionally runs four-step ionic
and isotropic cell trajectories without the mode. Trajectory, POSCAR, final TOML,
log and cell CIF files are byte-identical. Results and complete restart payloads
are byte-identical after replacing only the independently verified executable
SHA-256 and its derived build/input identity. Raw files necessarily differ in
that provenance because the executables differ; each restart checksum is checked.
No energy, density, force, stress, optimizer-state or numerical-history field is
excluded. The original exact restart/resume tests compare all files without
normalization within one executable.

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
cargo build --offline --release -p dft-cli
python3 -I P44A-evidence/identity.py final
python3 -I P44A-evidence/rejections.py
DFT_REFERENCE_BASELINE="$PWD/P44A-evidence/baseline-cli" \
  cargo test --offline --release -p dft-cli --test reference_mode_relaxation \
  -- --ignored --nocapture
DFT_NC_DIR="$PWD/../pseudo/pseudodojo/ex" \
  cargo test --offline --release -p dft-cli --test relaxation \
  -- --include-ignored --nocapture --test-threads=1
cargo test --offline --release -p dft-cli --test cell_relaxation \
  -- --include-ignored --skip real_files --nocapture --test-threads=1
scripts/check.sh
```

The two pre-existing external cell cutoff/rebuild studies (`real_files*`) are
not repeated: these are extended material qualification studies rather than
entry-guard regressions. All existing ionic tests, including real Si/MgO,
and the five synthetic cell tests are included above. Detailed command outcomes
are in `P44A-evidence/completion-report.json`.

The added cross-build harness initially rejected expected provenance differences,
then encountered an existing output directory and the cell restart's nested
identity field. These test-harness issues were corrected without any production
change, tolerance change or exclusion of scientific fields. The failed logs are
retained alongside the passing log.

Workspace verification: `scripts/check.sh` exited 0 (build, test, fmt, clippy).
The focused configuration/restart/equivalence tests and the legacy checkpoint
wire test passed. `scripts/tests/test_readme_figures.py` passed all nine tests;
README and artifacts were left unchanged. The initial equivalence test used
an unsupported `explain --json` flag and exited 101; the test was corrected
to the existing text explain interface and then passed. The missing-file runner
reported SKIP and exited 0. All changes remain uncommitted.

Completion verification after the approved guard: all four workspace gates passed
(exit 0), 4 configuration/entry tests passed, 6 existing ionic tests passed
(including real Si/MgO), 5 synthetic cell tests passed, and both cross-build
relaxation trajectories passed the byte comparisons described above. The
70 eligible SCF comparisons and 30 expected rejections were repeated successfully.
`relax.rs` was mechanically checked against HEAD: removing only the three-line
guard reproduces the original file byte-for-byte. No extra target directory
was created. Source-file inventory and raw-evidence inventory are recorded in
`P44A-evidence/touched-files.json`.

### Source files touched

- `crates/dft-cli/src/checkpoint_io.rs`
- `crates/dft-cli/src/config.rs`
- `crates/dft-cli/src/inspect.rs`
- `crates/dft-cli/src/preflight.rs`
- `crates/dft-cli/src/relax.rs`
- `crates/dft-cli/src/run.rs`
- `crates/dft-cli/tests/reference_mode_config.rs`
- `crates/dft-cli/tests/reference_mode_equivalence.rs`
- `crates/dft-cli/tests/reference_mode_relaxation.rs`
- `crates/dft-cli/tests/reference_mode_restart.rs`
- `crates/dft-engine/src/checkpoint.rs`
- `crates/dft-engine/src/scf.rs`
- `crates/dft-engine/src/session.rs`
- `crates/dft-engine/tests/reference_mode_materials.py`
- `docs/decisions/0038-reference-compat-mode.md`
- `docs/validation/reference-mode.md`
