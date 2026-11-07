# P27D: Radial convention score

Measured on 2026-11-07, Apple Silicon Mac, pure Rust release build. This is
an unpolarized PBE diagnostic with `pw92_prefactor = "original"`, not a
change of defaults or a claim that one convention is more correct.

## Result

**One rule set fits both materials together:** `local_mesh = "qe_like"`
and `radial_even_rule = "simpson_quadratic_end"`. It reduces the absolute
BAANDI-minus-QE residual in all six cases. Si residuals become
+3.37e-8, +4.44e-9, +2.60e-9 Ry; MgO residuals become
-2.05e-9, -8.66e-9, -9.31e-9 Ry. This is a candidate explanation of the
remaining offsets under this six-case protocol, not a universal accuracy
claim. No species-dependent selection, fitting, or tolerance change was used.

The mesh choice alone leaves Si unchanged and increases MgO's positive
residual. The even-count rule alone improves Si but overcorrects MgO to
roughly -7.7e-7 Ry. Both choices together bring MgO to the 1e-8 Ry level.
Their effects are consistent with the supplied note: Si's local count stays
1001; O's local count changes from 936 to 935, while the even projector
prefixes retain their lengths and change only their end integration rule.

## Complete SCF score

All energies are for the whole cell. Residual = 2*E_BAANDI(Ha) - E_QE(Ry).
Case names retain the existing input naming (e.g. si-20 means 20 Ha / 40 Ry
wavefunction cutoff). Every row uses the original PW92 prefactor, the same
input settings and external UPFs. All 24 runs exited 0 and converged.

| Case | local_mesh | radial_even_rule | Total (Ha) | Total (Ry) | BAANDI - QE (Ry) | SCF iterations |
| --- | --- | --- | ---: | ---: | ---: | ---: |
| si-20 | baandi | simpson_3_8 | -8.46231888703151 | -16.92463777406301 | +7.5936988253e-08 | 15 |
| si-20 | baandi | simpson_quadratic_end | -8.46231890817134 | -16.92463781634268 | +3.3657322263e-08 | 15 |
| si-20 | qe_like | simpson_3_8 | -8.46231888703151 | -16.92463777406301 | +7.5936988253e-08 | 15 |
| si-20 | qe_like | simpson_quadratic_end | -8.46231890817134 | -16.92463781634268 | +3.3657322263e-08 | 15 |
| si-30 | baandi | simpson_3_8 | -8.46233005664336 | -16.92466011328672 | +4.6713278579e-08 | 15 |
| si-30 | baandi | simpson_quadratic_end | -8.46233007777965 | -16.92466015555931 | +4.4406931465e-09 | 15 |
| si-30 | qe_like | simpson_3_8 | -8.46233005664336 | -16.92466011328672 | +4.6713278579e-08 | 15 |
| si-30 | qe_like | simpson_quadratic_end | -8.46233007777965 | -16.92466015555931 | +4.4406931465e-09 | 15 |
| si-40 | baandi | simpson_3_8 | -8.46233787255100 | -16.92467574510200 | +4.4897994655e-08 | 15 |
| si-40 | baandi | simpson_quadratic_end | -8.46233789370235 | -16.92467578740471 | +2.5952893168e-09 | 15 |
| si-40 | qe_like | simpson_3_8 | -8.46233787255100 | -16.92467574510200 | +4.4897994655e-08 | 15 |
| si-40 | qe_like | simpson_quadratic_end | -8.46233789370235 | -16.92467578740471 | +2.5952893168e-09 | 15 |
| mgo-40 | baandi | simpson_3_8 | -76.05088428700304 | -152.10176857400609 | +1.2759939239e-06 | 18 |
| mgo-40 | baandi | simpson_quadratic_end | -76.05088530944683 | -152.10177061889365 | -7.6889364209e-07 | 18 |
| mgo-40 | qe_like | simpson_3_8 | -76.05088390358867 | -152.10176780717734 | +2.0428226719e-06 | 18 |
| mgo-40 | qe_like | simpson_quadratic_end | -76.05088492602331 | -152.10176985204663 | -2.0466188744e-09 | 18 |
| mgo-50 | baandi | simpson_3_8 | -76.05097267573814 | -152.10194535147627 | +1.2685237323e-06 | 19 |
| mgo-50 | baandi | simpson_quadratic_end | -76.05097369775243 | -152.10194739550485 | -7.7550484434e-07 | 19 |
| mgo-50 | qe_like | simpson_3_8 | -76.05097229232564 | -152.10194458465128 | +2.0353487287e-06 | 19 |
| mgo-50 | qe_like | simpson_quadratic_end | -76.05097331433078 | -152.10194662866155 | -8.6615443706e-09 | 19 |
| mgo-60 | baandi | simpson_3_8 | -76.05098060606232 | -152.10196121212465 | +1.2678753478e-06 | 19 |
| mgo-60 | baandi | simpson_quadratic_end | -76.05098162807570 | -152.10196325615141 | -7.7615140981e-07 | 19 |
| mgo-60 | qe_like | simpson_3_8 | -76.05098022264967 | -152.10196044529934 | +2.0347006568e-06 | 19 |
| mgo-60 | qe_like | simpson_quadratic_end | -76.05098124465303 | -152.10196248930606 | -9.3060634754e-09 | 19 |

QE energies are the previously executed QE 7.6 results relayed in
`tests/reference/perf/dev08_materials_qe_relayed.md`, with version, input
hashes and UPF identities there. They are printed to eight decimal places
in Ry: the rounding uncertainty is approximately +/-5e-9 Ry. Residuals below
that scale should not be interpreted as resolved physical differences.
This task made no new QE run and did not read QE source. The only source
of alternate convention definitions was the architect's own-words note
`docs/research/qe-study/qe-radial-endpoint-convention.md`.

## Corrected polynomial gate (architect decision)

The original ticket's cubic-exact requirement for both even-count rules was
incorrect. The architect explicitly corrected it before implementation.
Simpson 3/8 integrates degree-three polynomials exactly; the quadratic
interpolant for the last interval integrates only through degree two.
The preceding Simpson 1/3 panels remain cubic-exact. For f(x)=x^3 on
x=0,1,2,3, the quadratic-end total is 20.5 whereas the exact integral and
Simpson 3/8 result are 20.25. On a uniform physical grid of step h, its
cubic error (computed minus exact) is +h^4/4 for a monic cubic.
The specified weights (-1,8,5)/12 were retained without adjustment.

Unit tests check degrees 0..3 for point counts 3..19 and steps
0.125, 1, 2, including the expected cubic defect for even counts;
odd-count results are bit-identical. Impulse probes scaled by 12 verify
exactly the three last-interval weights. The two-point fallback remains
trapezoidal (exact through degree one); an odd-forced local count with fewer
than two samples is rejected.

## Implementation and scope

```toml
[pseudo]
local_mesh = "qe_like"                   # default: "baandi"
radial_even_rule = "simpson_quadratic_end" # default: "simpson_3_8"
```

The knobs are independent and reject unknown strings and unknown table keys.
New additive entry points dispatch defaults to the unchanged original
functions. The local view and ionic-local reconstruction use the same
endpoint; the core mesh remains full length; stored projector prefixes are
unchanged. Alternative quadrature is selected for local, core and projector
transforms. No sampling tables, pseudopotential values or defaults change.
The engine's common setup passes these choices to all SCF modes and NSCF;
the measured evidence here is only unpolarized PBE.

The enum fields default to the legacy choices. CLI results include both
names under `pseudo` and `scf.pseudo` when either is nondefault; omission
means both defaults and preserves the legacy default JSON shape. Other
result profiles retain `scf.pseudo` via the common options serializer.
Checkpoint options serialize nondefault strings; missing fields in old
checkpoints decode as defaults. The checkpoint version is unchanged.
Exact restart compares both knobs independently in both directions;
changing them requires the existing explicit warm-start mechanism.
Tests cover roundtrip, legacy omission and unknown-value rejection.

The derivative model identity includes both choices. Analytic force/stress
routines have separate radial derivative implementations not extended by
this diagnostic. They explicitly reject nondefault research conventions;
thus a calculation cannot silently report derivatives of a different
Hamiltonian. Default derivatives are unchanged. No force/stress capability
or spin/LDA validation is claimed here.

## Reproduction and raw evidence

No UPF was copied or committed. Set DFT_NC_DIR to an external directory
containing the three PseudoDojo files and run from the repository root:

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
cargo build --offline --release -p dft-cli
DFT_NC_DIR=/Users/mino/Develop/rust-dft-wt/pseudo/pseudodojo/ex \
  python3 -I crates/dft-engine/tests/radial_convention_score.py \
  --evidence /tmp/P27D-evidence/materials-final
```

The runner skips explicitly if the variable or files are missing; invalid
present data or an unsuccessful SCF fails. It runs the 24 requested cases
plus six omission runs to compare both explicit defaults against omission.
It preserves input TOML, stdout/stderr, results, density files, hashes and
command lines, and writes `score.json`. It rejects an evidence path inside
`target/`. Original inputs are read from `tests/reference/materials`.
No extra target directory was created. `/tmp/P27D-evidence` is outside the
build tree and is temporary local storage; preserve it before OS cleanup.

Raw evidence locations:

- `/tmp/P27D-evidence/materials-final/score.json`: all rows, commands, hashes, identity result.
- `/tmp/P27D-evidence/materials-final/<case>-<mesh>-<rule>/`: complete SCF files.
- `/tmp/P27D-evidence/baseline-cli`, `baseline.json`, `baseline.log`: pre-change build and T=1 benchmark.
- `/tmp/P27D-evidence/final-default-t{1,2,4,8,12}.{json,log}`: default identity comparisons.
- `/tmp/P27D-evidence/final-identity-commands.json`: benchmark commands and exit codes.
- `/tmp/P27D-evidence/{check,upf-tests,physics-tests,config-tests,readme-tests,score}.log`: verification logs.

UPF SHA-256 (matches relayed QE evidence):

- Si: `39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282`
- Mg: `15cdd9380687a47a1f4bc4ab5fe1f0c6d9c50891350fd8534c8f292b3ce0be14`
- O: `2a543889f6522e2a1f94e2c2ff1477e61cb89f25c67b27c54940e8263ec6fb11`

## Verification

The baseline command used the unchanged executable, the five cases below,
T=1, `--repeats 1 --record /tmp/P27D-evidence/baseline.json`.
All five subsequent commands exited 0:

```sh
python3 -I scripts/bench_scf.py --binary target/release/dft-cli \
  --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF \
  --cases gamma,k222s,k222,k444s,si8 --repeats 1 --threads T \
  --baseline /tmp/P27D-evidence/baseline.json --require-identical \
  --record /tmp/P27D-evidence/final-default-tT.json
```

T = 1, 2, 4, 8, 12: **25/25** comparisons have exactly zero energy change
(0 Ha) and byte-identical density files. No timing improvement is claimed;
the host was shared with other checks.

Final verification (all commands below exited 0):

| Command | Result |
| --- | --- |
| `scripts/check.sh` | PASS build, test, fmt, clippy |
| `cargo test --offline -p dft-upf` | Includes polynomial, specified-weight, odd-identity and existing transform tests |
| `cargo test --offline -p dft-physics --test radial_convention_score` | Two tests passed: endpoint/ionic consistency and core/projector polynomial oracle |
| `cargo test --offline -p dft-cli --test pseudo_convention_config` | Two tests passed: strict config and exact/warm restart in both directions for each knob |
| Material runner above | 24 converged score SCFs plus six omission SCFs; all six explicit-default identities passed |
| Five benchmark commands above | 25/25 zero-energy-change, byte-identical density comparisons |
| `python3 -I scripts/tests/test_readme_figures.py` | Nine tests passed; existing unclosed-file ResourceWarnings |
| Material runner without DFT_NC_DIR and with a nonexistent directory | Explicit SKIP in each case, exit 0 |
| `git diff --check` | No whitespace errors |

An additional JSON inspection confirmed the convention metadata on all 24
outputs. The scored executable is preserved as
`/tmp/P27D-evidence/final-cli`, SHA-256
`4b41fd005258634bdab5567a6d329f1e4b2eab975e34d604baf818d378540ee0`.
`source-preservation.json` records that the original radial production
source and original physics Hamiltonian/data source remain unchanged
prefixes of their files. The default branch calls those original routines.

During implementation, an initial generated helper had a delimiter error,
and an overbroad edit initially inserted option validation into a volumetric
validator; both compile failures (exit 101) were corrected before successful
SCF runs. The first workspace check exited 1: its doctest could not find
the engine crate after a test-module move overlapped the build. The subsequent
full check with source edits stopped passed all four gates. The original
failure log is retained in `check-initial.log`; no test, tolerance or gate
was relaxed. README is intentionally unchanged per ticket.


The final release rebuild changed the binary hash after the test-module
relocation. The full 24-case score, six explicit-default comparisons, and
25 benchmark identities were therefore repeated with the final executable.
Final evidence uses `materials-final`, `score-final.log`,
`final-default-t*.{json,log}` and `final-identity-commands.json`.
The earlier `materials` and `scored-cli` artifacts are retained separately.

## Changed-file inventory

Only additive numerical entry points and option dispatch were introduced.
All pre-existing public function signatures are unchanged. The following
files carry the implementation or its verification:

- UPF: `crates/dft-upf/src/radial.rs` (enums, additive transforms, quadrature tests).
- Physics: `crates/dft-physics/src/hamiltonian.rs` (local/projector dispatch),
  `crates/dft-physics/src/data.rs` (core dispatch),
  `crates/dft-physics/tests/radial_convention_score.rs` (synthetic oracles).
- Engine: `crates/dft-engine/src/lib.rs` (enum facade),
  `crates/dft-engine/src/scf.rs` (defaulted options, setup and result fields),
  `crates/dft-engine/src/session.rs` (result metadata),
  `crates/dft-engine/src/checkpoint.rs` (exact option comparison),
  `crates/dft-engine/src/force_identity.rs` (model identity),
  `crates/dft-engine/src/forces.rs` and `crates/dft-engine/src/stress.rs`
  (explicit derivative rejection for nondefault research options),
  `crates/dft-engine/tests/radial_convention_score.py` (external SCF runner).
- CLI: `crates/dft-cli/src/config.rs` (strict keys),
  `crates/dft-cli/src/run.rs` (engine options and output metadata),
  `crates/dft-cli/src/checkpoint_io.rs` (backward-compatible wire fields/tests),
  `crates/dft-cli/tests/pseudo_convention_config.rs` (config/restart tests).
- Documentation: this file and `docs/decisions/0037-radial-convention-options.md`.

Changes are uncommitted on `task/P27D-RADIAL-CONVENTION-SCORE` for review.
