# Python/ASE: certified energies compared with Classic

## 1. Problem and motivation

An ASE workflow should receive the same energy as the agent-free CLI and must
not mistake unavailable forces or unconverged SCF for useful scientific labels.
This case tests a native session, stepwise execution and an ASE energy calculator.

## 2. Setup

Two deliberately synthetic systems vary cell shape, ion count and projectors:
cubic H-labelled two-ion valence-one, s projector; skew He-labelled one-ion
valence-two, s/p projectors. These labels do not make the synthetic potentials
physical H or He references. Generator radius=1.4 Bohr, 401 points; cutoff=4 Ha,
four bands, Gamma, LDA-PZ. Exact options, source hashes and UPF-byte SHA-256 values
are in [python-setup.json](data/recorded/python-setup.json) and
[provenance.json](data/recorded/provenance.json). The native comparisons use
Bohr cells; ASE comparisons use 2 Angstrom cubic and rows (2,0,0), (.25,2,0),
(0,.25,2.5) Angstrom. These are different physical cells; do not mix their energies.

The runner builds the native wheel and CLI in the debug profile, installs the
wheel into its own output directory without dependencies or network, and runs
all five installed-wheel tests. Supply an existing interpreter with ASE/maturin.
The wheel hash, compiler versions, interpreter and package versions are recorded.

## 3. Classic and external path

The CLI `run` without mode flags is Classic. For each system the workflow loads
identical generated potential text and geometry, executes native batch and
stepwise sessions and ASE `get_potential_energy`, and compares against the CLI.
A separate timed pair uses precisely the same ASE geometry in the Classic input.
Exact inputs, Classic results and CLI logs are retained under `data/recorded/ase-*`.

## 4. External actions and provenance

The scripted client constructs options and potential mappings, initializes the
session, steps to completion, inspects the certificate-backed result, and obtains
an ASE energy. Native results retain owned eigenvalues after sessions are dropped.
ASE converts units through the binding's engine wrappers. No agent is required,
and no policy intervention is submitted or implicitly approved. This binding
has no public audited control-request API; the run log is a test/workflow log,
not a runtime control transcript. All decisions are fixed by the script.
See [python.log](data/recorded/python.log), [ASE comparisons](data/recorded/ase.json),
[native comparisons](data/recorded/native.json) and [paired measurements](data/recorded/python-timings.json).

## 5. Results and timing boundaries

The table below comes from the separate paired measurements. Tests also check
bitwise energy equality in Ha/eV and eigenvalue equality, and charge to ten decimal
places. Historical numerical verification uses 1e-10 Ha and identical iteration
counts; same-run energy comparison is stronger (bitwise).

Wall times are single sequential observations: CLI includes process startup and
file output; ASE is in-process and includes calculator construction. The signed
difference is end-to-end workflow cost, **not isolated binding overhead**. These
short, warmed, sequential runs cannot support a speedup claim or timing threshold.

| System | Energy (Ha) | Iterations Classic/ASE | Delta (Ha) | Classic (s) | ASE (s) | Signed workflow difference (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| cubic_h2 | -0.9120966424598976 | 12/12 | 0.0 | 0.068530 | 0.057639 | -0.010891 |
| skew_he | -1.4971019807693309 | 8/8 | 0.0 | 0.073323 | 0.061977 | -0.011346 |

## 6. Failures

All five tests execute without skips. They reject malformed shapes and zero
iteration budgets, retain an uncertified one-step native outcome, and require
ASE to raise ConvergenceError with no published result for that outcome.
Forces, stress, free energy, dipole and charges raise NotImplementedError.
No missing force becomes a zero training label. No stale control request is
applicable because this Python binding does not expose audited interventions.

## 7. Limits

Two small synthetic insulator fixtures on one Mac; no real materials, MD, force
training, MLIP loop, PBE, k meshes, remote execution or live agent. The wheel is
not a portability/distribution claim. Negative timing differences include startup
and cache effects. Unsupported capabilities remain unavailable.

## 8. One-command reproduction

```sh
python3 -I scripts/run_case_studies.py python --verify --python /path/to/python
```

The runner fails if tests skip; a clean exit with an unavailable extension is not
accepted. Archived TOMLs contain expired temporary paths; regenerate them.
