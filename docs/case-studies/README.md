# External-control case studies

P11L is **partial, not accepted**. These documents distinguish executed evidence
from missing implementations and missing comparisons. Every number here comes
from a new P11L run, not from earlier tickets. No QE reference is claimed.

| Path | Evidence in this revision | Remaining acceptance gap |
| --- | --- | --- |
| [Python/ASE](python-ase.md) | Two synthetic systems; Classic/native/ASE equality, paired wall times, failure checks | No live agent or runtime control audit (binding does not expose that capability) |
| [Agent-assisted preflight](preflight.md) | Two systems; Classic before/after review, proposals, approved/rejected sidecars, stale rejection | Scripted review only; candidate intentionally never applied |
| [Self-Healing SCF](self-healing.md) | Repeated budget recovery, full events/audit and certified results | One system; reference uses AgentControlled, not Classic; no instability recovery |
| [MCP control](mcp.md) | Nine newly executed adapter integration checks | No matched Classic energy/timing pair, second system or persisted wire/audit transcript |
| [MLIP active learning](active-learning.md) | Not executed | Learning/deployment workflow and MLIP implementation absent |
| [Low-level coupling](low-level-coupling.md) | Not executed | Public field/operator coupling API implementation absent |

## Reproduce available evidence

From the repository root, with the pinned Rust toolchain, sibling offline vendor
directory and a Python environment already containing ASE and maturin:

```sh
python3 -I scripts/run_case_studies.py all --verify --python /path/to/python
```

The command prints a fresh temporary output directory. `--output /new/directory`
selects a persistent destination and refuses an existing directory. Individual
selectors are `python`, `preflight`, `self-healing`, and `mcp`. No network installs
or UPF downloads occur. Python wheels install only into the output directory.
The recorded interpreter path and exact nested commands are in
[data/recorded/provenance.json](data/recorded/provenance.json). Build logs record
compiler/Cargo versions. The source SHA-256 map binds code and synthetic generators;
the source commit is not a claim that these uncommitted documents were committed.

`--verify` checks energy within 1e-10 Ha and recorded iteration counts for
self-healing and Python test results; preflight checks energy within 1e-10 Ha.
Within each paired execution, stronger bitwise energy equality assertions apply.
The installed-wheel suite must execute five tests with zero skips. MCP and
preflight must execute nine and eleven tests respectively, with zero ignored
cases. Wall times are observations, not regression thresholds or speed claims.
An unavailable dependency or failed assertion returns nonzero and retains a
failed provenance record. A passing command means selected implemented checks
passed, **not** that all six case studies meet P11L acceptance.

The committed bundle retains text/JSON evidence, not wheels, installed packages,
density binaries or UPFs. Input TOML and proposal paths intentionally retain the
original temporary pseudopotential path to preserve the recorded hashes. These
paths expire; reproduce with the runner, which regenerates the fixture and all
path-dependent hashes. Do not execute archived TOMLs unchanged.

## Integration handoff

README's required link to this list is still pending: README.md is outside this
ticket's explicit OWNED_FILES. Suggested link label: "External-control case studies";
repository-root target: `docs/case-studies/README.md`.
No engine feature, default, tolerance or reference README result was changed.
Dashboard/registry updates likewise belong to the orchestrator; this checkout's
registry has no P11L running entry. No recurring session timer was installed.
