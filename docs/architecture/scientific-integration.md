# Scientific integration contracts

Status: design only, not implemented. No Python, MCP, MLIP or coupling code is authorized by P11R. Adapter calls ultimately wrap `ScfProblem::new/with_kpoints`, the proposed session, `ScfOutcome` and `EnergyAudit`. [ADR 0017](../decisions/0017-coupling-commitments.md) remains authoritative for physical roles and energy ownership.

## Ownership boundary with P11I (D1/D2/D3)

P11I-COUPLING-API-SPEC owns the canonical high-level workflow/ASE/Python contract (`docs/architecture/workflow-api.md`), low-level coupling contract (`docs/architecture/coupling-api-lowlevel.md`) and MLIP/active-learning contract (`docs/architecture/active-learning-contract.md`). Those files are not present in this checkout when this addendum is written; see the P11I ticket on main `98f0b50`. Do not treat the absent files as implemented features. The first P11R draft's duplicate structure, prediction and field request schemas are withdrawn in favor of that single owner.

P11R supplies the runtime adapter boundary: validate a P11I scientific request into `ScfProblem::new/with_kpoints`, drive the one session under a declared mode, and return qualified `ScfOutcome`/`EnergyAudit` references with provenance. Modes never change P11I's physical interpretation. P11G implements optional Python/ASE binding; P11F implements MCP over the same driver; P11H presents CLI capabilities. None creates a second SCF/control implementation.

Integration invariants retained here, not competing schemas:

- ASE units convert only at I/O boundaries; unavailable forces/stress raise explicit unsupported errors, never zeros. Certified energy is distinct from partial unconverged output. P11I owns shape/layout/cache/observable details and P11G validates binding conventions.
- MLIP uncertainty must distinguish calibrated from uncalibrated/unavailable and retain calibration reference; no silent model mixing or zero-filled DFT force labels. Producer/model identity survives runtime caching, provenance and workflow edges. P11I owns the prediction/dataset contract.
- External numerical solvers use P11I's low-level contract, not agent prompt/event text. The runtime enforces resource bounds, immutable environment per inner calculation, and separate inner certificate versus outer status. No LLM or HTTP dependency enters the tight numerical exchange path.
- E1 records P11I request/response artifact digests; E2 references P11I typed artifact identities. They do not redefine field layouts, force terms, zero-copy lifetimes or active-learning data schemas.

## Provenance E1 proposal

`ProvenanceRecord`: `schema_version:string`, `provenance_id:string` (content digest), `calc_id:string`, `request_id:string`, `parent_ids:string[]`, `input_sha256:string`, `model_digest:string`, `environment_digest:string`, `source_commit:string`, `dirty_patch_sha256:string|null`, `lockfile_sha256:string`, `compiler:string`, `backend:string`, `execution:ExecutionIdentity`, `pseudopotential_sha256_by_species:map<string,string>`, `initialization:{kind:uniform/exact/warm,source_id:string|null}`, `policy_log_sha256:string`, `output_manifests:string[]`, `status:converged/exhausted/failed/cancelled`, `certificate_present:bool`, `created_utc:string`. This fills metadata missing from `ScfResult`; values are adapter-supplied facts, never inferred from energy. Exclude timestamps/host data from content identity; canonical schema ordering and round-trip float representation are versioned. Record exact command as data with secrets redacted, never evaluate it on import. UPF bytes are not embedded or sent to agents.

## Workflow DAG E2 proposal

`Workflow`: `schema_version:string`, `workflow_id:string`, `nodes:Node[]`, `edges:Edge[]`, `budget:{max_nodes:u64,max_total_iterations:u64,max_wall_seconds:u64,max_artifact_bytes:u64}`. `Node`: `id:string`, `kind:scf/inspect/export/mlip/external_solver`, `request_digest:string`, `input_artifact_ids:string[]`, `requested_observables:string[]`, `retry_limit:u64`, `requires_certified_inputs:bool`. `Edge`: `from:string`, `to:string`, `artifact_role:string`. Node status is `pending/ready/running/succeeded/failed/cancelled/blocked`; result references a ProvenanceRecord. Validate unique IDs, acyclicity, typed artifact compatibility and aggregate quotas before scheduling. Unknown node kinds are rejected; no shell or arbitrary-code node.

Execution wraps session creation and results; independent ready nodes may run only under host resource limits. Repeated request id plus equal canonical payload deduplicates; differing payload conflicts. Retries have new attempt ids and lineage, preserve failures, and cannot silently change model/policy. Inner SCF certification is separate from workflow success and outer coupled convergence. DAG caching needs matching full model/execution requirements; an old unconverged result never satisfies a certified dependency. No scheduler currently exists.

External numerical solver D3 uses the separate [contract](../research/external-solver-contracts.md). These adapter proposals refine ADR 0017's deferred runtime without authorizing field insertion or new physical models.

## Python usage (P11G implementation)

The separate `baandi-py` crate wraps the engine's Gamma-only LDA session.
The native binding and ARM64 CPython 3.9 wheel have been exercised. **ASE energy
is pending the `dft-units` dependency approval**: the proposed bridge is in
`crates/baandi-py/units-integration.patch`. Until applied, `Baandi` rejects
unsupported properties but cannot calculate energy. This does not implement
the full P11I workflow schema or audited policy intervention API.

Build from `python/baandi` (do not pass `--manifest-path`, which changes maturin's
Python project discovery). ASE remains an optional extra:

```sh
cd python/baandi
PYO3_PYTHON=/Users/mino/Develop/rust-dft-wt/pyenv/bin/python maturin build --offline --locked
# Install the resulting wheel with pip; request baandi[ase] for the ASE adapter.
```

The current shared vendor directory lacks PyO3 despite its presence in the local
Cargo registry cache. An offline Cargo invocation from outside the worktree,
with an absolute manifest path, uses that cache. See the crate validation note
for the tested local wrapper and the Python framework link workaround. No new
packages were downloaded.

```python
from pathlib import Path
from baandi import Options, Session

options = Options()
options.cutoff_hartree = 5.0
session = Session(
    cell_rows_bohr=[[8., 0., 0.], [0., 8., 0.], [0., 0., 8.]],
    species=["He"], positions_fractional=[[0., 0., 0.]],
    pseudopotentials={"He": Path("He.UPF").read_text()}, options=options,
)
session.initialize()
while session.status == "Ready":
    session.step()  # exactly one full SCF iteration, with engine certification
result = session.result()  # terminal output, including an exhausted attempt
if not result.converged:
    raise RuntimeError("SCF exhausted; result is uncertified")
print(result.energy_hartree, result.eigenvalues_hartree)
print(session.density_summary())
```

`run()` initializes a Created session and steps until terminal; it also resumes
an initialized session. Status values are Created, Ready, Converged,
NotConverged and Failed. Results are owned copies; eigenvalues are Gamma-only
in Ha. Density summary returns `(shape, min, max, mean, integrated_electrons)`;
the density values are Bohr^-3. At an active boundary it summarizes the next
input density; at termination it summarizes the evaluated output. No full
orbital or density arrays are exported. No background SCF runs between calls.

Options are copied on construction and retain engine defaults except explicitly
set cutoff, band counts, iteration limit, tolerances, seed and Pulay settings.
Changing a Python Options object afterwards never mutates a live calculation.
No alternate policy setter is exposed. The engine's non-Send mixer remains on
one dedicated worker thread per session; only owned data crosses channels.
Waiting for setup and numerical steps releases the GIL. Dropping the session
closes the command channel and the worker releases the engine. There is no
handwritten unsafe code or Python callback inside the numerical loop.

After the proposed units bridge is approved and applied:

```python
from ase import Atoms
from baandi.ase import Baandi
atoms = Atoms("He", cell=[4., 4., 4.], pbc=True)  # Angstrom, row vectors
atoms.calc = Baandi(pseudopotentials={"He": Path("He.UPF").read_text()})
energy_ev = atoms.get_potential_energy()
```

ASE requires a fully periodic cell and rejects initial charge or magnetic
moments. It copies UPF text/settings at construction and disables energy caching
until complete model-identity caching is available. Forces, stress, free energy
and every other unsupported property raise `NotImplementedError`, including
force-consistent energy requests; they never return zero. An exhausted SCF
raises `ConvergenceError`, leaves `results` empty and retains its session for
uncertified diagnostics. The intended conversion uses exactly `dft-units`, not
ASE's different CODATA constants. The CLI comparison test checks hexadecimal
float equality in Ha and eV with the same settings and represented geometry.
