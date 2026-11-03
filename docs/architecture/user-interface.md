# Researcher-facing interface architecture

| Item | Design decision |
| --- | --- |
| Status | P9A-USER-INTERFACE-DESIGN; proposal, no feature implementation |
| User outcome | Import a known format, choose a material explicitly, review an editable calculation, then run and view links |
| Trust boundary | Every frontend produces the same validated plan; unsupported physics cannot be approved away |
| Core boundary | Existing numerical crates receive typed local data only; no network, serialization, language model, viewer or workflow logic |
| First delivery | Import translation reports plus P8G species identity/preflight; recipes and database search follow |
| Companions | [Evidence and dated references](../research/ui-landscape.md), [implementation roadmap](../roadmap/ui.md) |

## 1. Scope and current truth

This design follows ADR 0017's adapter isolation and explicit unavailable/error semantics, ADRs 0012–0015's units/expressions/k-point contracts, and the standard-output policy in AGENTS.md. New commands and schemas below are **proposed**. No new crate is authorized here: initial modules live in `dft-cli` or an external process adapter; any future adapter crate needs a new ownership/manifests ticket. The numerical dependency graph is unchanged.

The inspected CLI supports TOML, `run`, `validate`, `grid`, `inspect-upf`, `eval`, and `convert`. Configuration supports spin-unpolarized, insulating, norm-conserving LDA-PZ/PBE with Gamma, uniform grids or explicit weighted k points. Compatibility remains constrained by the UPF subset and validated material coverage; PBE support does not imply a library of suitable PBE potentials. Older help/architecture prose is stale. `validate` alone is not full preflight. P8A bands/DOS and P8B cube/XSF are pending in this checkout. Forces, stress, relaxation, spin, +U, smearing/metal occupations, PAW and USPP must not be advertised as executable.

## 2. Layers and owners

| Layer | Proposed owner | Responsibility and boundary |
| --- | --- | --- |
| L0: explicit input | Existing `dft-cli::config`; `dft-units` conversions | Keep existing TOML and command behavior. The resolved plan materializes defaults and units. Core consumes the same typed inputs as today. |
| L1: import | New CLI `import/` modules; optional external format helper | Standard structures and foreign-code input to normalized structure/settings plus translation report. No solver, potential conversion or framework dependency in core. |
| L2: recipes | CLI `request` and `recipes` modules | Closed bilingual request grammar, versioned recipes, explicit editable TOML and workflow/view manifest. Never run during compilation. |
| L3: guided mode | CLI `guide` module | Wizard/REPL over L1/L2/preflight; immutable drafts and explicit selections, no alternative hidden state machine for physics. |
| L4: agent interface | CLI JSON interface; external `adapters/mcp/` process | Capability/schema discovery, import/search/plan/preflight/run/status/view tools. General language interpretation belongs to the assistant, not engine. |
| L5: result viewing | CLI `report`/`view`; optional external local server | Read completed output bundles and standard files, render offline HTML, SVG/PNG and links. No numerical kernel imports a viewer. |
| Shared acquisition | Optional CLI `providers` feature, later separately approved adapter crate | HTTP, provider policies, cache and provenance. Produces data for L1 and cannot call SCF. |

Data flow: user text/file → request or import → structure selection → capability check → pseudopotential identity → resolved input + translation/recommendation reports → full preflight + estimate → approval bound to plan digest → existing engine → output bundle → viewer. Requests that fail early may still export a clearly blocked draft, never a runnable plan.

## 3. Import and translation contract

Every import, including geometry-only, emits human-readable `translation.md` and versioned `translation.json`. Each source field/record is accounted for exactly once, with source byte digest, format/dialect/version, location, original value/unit, destination path/value/unit, category, reason and diagnostic code. Categories:

- **mapped**: equivalent representation, including explicit conversion and atom/cell permutation;
- **defaulted**: absent source value, engine/recipe policy and version, proposed value, acceptance state;
- **ignored**: recognized nonphysical metadata with a reason, retained in provenance, visible to the user;
- **unsupported**: unrecognized or non-equivalent semantics; blocks executable translation.

No “best effort” dropping of unknown tags. Unsupported physics cannot become ignored after a generic yes. A user may separately request **structure-only** extraction; preserve the original blocked full-import report and record discarded calculation intent. Imported settings never imply equal QE/VASP numerical results. A lossless geometry mapping is distinct from a changed Hamiltonian.

| Input | First supported subset and mandatory diagnostics |
| --- | --- |
| POSCAR/CONTCAR | Cell scaling (including negative-volume convention), direct/Cartesian positions, explicit species order. Missing element symbols require supplied mapping; never infer from filenames or POTCAR. Selective-dynamics flags, velocities and predictor data are inventoried; relaxation/MD intent is unsupported. Multiple scaling forms need tests or explicit refusal. |
| CIF | Select data block; expand documented symmetry operations through a reviewed parser; record tolerance, site multiplicity and resulting cell. Reject partial occupancy, unresolved disorder, ambiguous species and unsupported magnetic CIF. Preserve uncertainty/measurement metadata; do not round occupancies to integers. |
| XYZ/extxyz | Plain XYZ needs explicit periodic cell and units; no invented vacuum. extxyz maps lattice, positions, species and PBC with documented matrix convention. Require frame selection; refuse nonperiodic/mixed-PBC execution for the current 3D engine. Inventory attached forces, charges, constraints and magnetic moments. |
| XSF | Explicit periodic cell and Cartesian species/positions; select animation frame. Unsupported dimensionality is blocking. Data grids remain named auxiliary artifacts, never confused with atoms or automatically used as a restart. Official specification check is pending [R2, F9]. |
| QE pw.x | Start with `ibrav=0`, `CELL_PARAMETERS`, `ATOMIC_POSITIONS`, `ATOMIC_SPECIES`, explicit units and `K_POINTS gamma/automatic/crystal`. Other `ibrav` and path modes need implemented/tested expansion or rejection. `ecutwfc` converts Ry at the boundary; `ecutrho` only maps if the target density-basis policy is equivalent. Do not equate `conv_thr` with engine energy tolerance. Spin, Hubbard, smearing, relaxation, PAW/USPP, unsupported functionals and fields block. Paths are inert until preflight. |
| VASP POSCAR + INCAR + KPOINTS | Map geometry, approved functional intent, explicitly supplied cutoff (ENCUT eV) as a suggestion pending NC-potential validation. Enumerate actual reciprocal points to verify Gamma/Monkhorst-Pack parity and shifts against ADR 0013; do not copy mesh labels blindly. Line mode requires P8A. EDIFF/ALGO/ISMEAR/ISPIN/LDAU/IBRION are not equivalent engine controls. Unknown/nested tags block full import. **Never open, inspect, hash, read or redistribute POTCAR**; no automatic directory scan of it. User supplies compatible NC UPFs explicitly. |
| ASE JSON | Explicit versioned ASE database/Atoms dialect and record selector; data-only cell, PBC, numbers, positions. Calculator results/constraints are inventoried, never executed or treated as current engine results. Unsupported dialect rejected. |
| pymatgen JSON | Explicit Structure/MSON dictionary dialect, lattice/sites/species/occupancy; no dynamic module loading from `@module`/`@class`, pickle or arbitrary object construction. Reject unsupported object types and disorder. |

All structures retain elements separately from species aliases. Normalize units at the I/O boundary (input lattice vectors are rows, internal vectors columns; reciprocal vectors include 2π). Preserve original cell, coordinates and atom ordering; any primitive-cell reduction or symmetry operation is an explicit recorded transformation, not an optimization hidden in import. Files have byte/atom/depth limits; reject nonfinite coordinates, singular cells and malformed compressed payloads. Tests use independently constructed skew/left-handed cells, unit conversions and specification fixtures plus independent readers where permitted. No upstream numerical source is needed [R2, F1–F8].

## 4. Deterministic requests and recipes

The convenience language is a finite grammar with documented aliases, not a claim of understanding unrestricted English/Korean. Normalize Unicode NFC and whitespace; keep element capitalization and provider IDs intact. Tokenize quoted paths, IDs, numbers and a closed bilingual vocabulary. Reject leftover text, ambiguous aliases, conflicting clauses and unknown units. No shell interpolation or arbitrary expression execution; quantity expressions reuse ADR 0014 only in delimited quantity fields.

Proposed grammar v1 (braces mean optional clauses, not literal syntax):

```text
calculation := action target {"with" xc} {"at a" mesh "mesh"} {"along" path}
action      := "scf" | "band structure of" | "dos of" | "relax"
target      := quoted-file | material-id | formula {material-id} | known-alias
mesh        := positive-int "x" positive-int "x" positive-int
xc          := "LDA-PZ" | "PBE"
path        := "the standard path" | quoted-path-id
view        := "show the density" | "show the density slice" | "plot bands"
ko-calc     := target {xc "로"} {mesh "격자로"} ("계산" | "구조 최적화")
ko-bands    := target "표준 경로로 밴드 구조"
ko-view     := "밀도 보여줘" | "밀도 단면 보여줘" | "밴드 그려줘"
```

Finite aliases map `silicon`/`실리콘` to element Si for **search**, not to an assumed polymorph. Bare formula input in the guide means search, not execute. Unrestricted paraphrases are handed to an external agent only if the user chooses that interface.

| Request | Deterministic outcome |
| --- | --- |
| `relax Fe3O4 mp-19306 with PBE at a 4x4x4 mesh` | Parse action, formula, ID, XC and grid; reject relaxation capability before any solve. If acquisition is requested, verify ID/formula agreement; mp-19306 is illustrative, not a verified Fe3O4 structure here. |
| `Fe3O4 mp-19306 PBE로 4x4x4 격자로 구조 최적화` | Same typed request and same blocking diagnostics. Korean particles shown in grammar are the accepted spelling. |
| `band structure of silicon along the standard path` / `실리콘 표준 경로로 밴드 구조` | Ask which structure/cell. P8A prerequisite; a standard-path recipe needs a named convention, geometry check, transformation and tolerance. Initially only a validated preset, not arbitrary symmetry inference. |
| `show the density` / `밀도 보여줘` | Resolve explicitly selected completed run; emit a view manifest referring to its density. Never trigger SCF. |
| `show the density slice` / `밀도 단면 보여줘` | Guide requests plane, offset and units when absent; propose fractional c=0.5 visibly, do not silently choose it. |

Compiler output is an editable `input.toml` for each executable calculation, plus `request.json` (original text, language, grammar/recipe versions, normalized intent), `plan.json` (workflow steps, resolved references, status) and reports. All defaulted SCF controls must be printed. A blocked request gets a manifest listing unresolved fields and no runnable input. A view-only request gets a printable/editable `view.json`; it must not manufacture a new calculation input. Editing TOML invalidates plan approval; reparse and preflight it. Same request + selected source bytes + recipe versions + capability snapshot yields the same plan digest (timestamps live outside digest inputs).

Recipes initially: fixed-geometry insulating SCF; cutoff/k convergence study as an explicit finite batch; bands/DOS after P8A. Recommendations are not convergence certificates. Record recipe identity/version and every chosen setting with origin `user`, `import`, `provider`, `registry`, or `recipe`.

## 5. Capability negotiation and pseudopotentials

`capabilities --json` is generated from tested build features and adapter registrations, not a language-model answer. Each capability has `supported`, `unsupported`, `experimental` or `unavailable_in_build`, constraints, evidence pointer, and roadmap pointer. Distinguish parsing an action from executing it. Planning and run dispatch both compare requirements against the same snapshot; recheck at execution so changing binaries cannot bypass validation.

Required capability dimensions include XC, spin, occupations, charge, potential type/subset, force/stress availability, k sampling, output quantity, import dialect, provider and viewer support. Unknown magnetic/metal metadata means unknown suitability, not permission to assume an insulator. A finite database band gap does not prove convergence or validity for another functional. Imported magnetic moments or +U settings must remain blockers even if the user omitted those words in the request.

Example refusal (localizable):

> Cannot run this Fe3O4 request. It requires relaxation (forces/optimizer), magnetic spin states and, if requested, Hubbard U. This build supports fixed-geometry, spin-unpolarized insulating NC LDA-PZ/PBE only. Spin, +U and metallic occupations are planned in ESLE phase E; forces/relaxation need a separate physics ticket (UI roadmap deferred gate PHYS-RELAX). You may inspect/export the selected structure. No calculation has started.

Korean equivalent: `이 Fe3O4 요청은 실행할 수 없습니다. 구조 최적화·자기 스핀·요청한 +U 기능이 필요합니다. 현재는 고정 구조의 비자성 절연체 NC LDA-PZ/PBE만 지원합니다. ESLE E 단계와 PHYS-RELAX를 참고하세요. 구조 확인·내보내기는 가능합니다. 계산은 시작하지 않았습니다.`

Do not claim all Fe3O4 IDs have identical ordering; list requirements actually requested or documented by the selected record. A separate nonmagnetic model would be a different explicit research request with a new fingerprint, not a “continue anyway” fallback. Model/state/branch/request identities follow ESLE; near-equal energies do not identify equivalent states, and a database stable flag does not certify an engine ground state.

**P8G prerequisite:** atom/species expected element is checked against `PP_HEADER` through the existing UPF parser; aliases such as `Fe_A` remain valid only with explicit element Fe. Check every species, valence/electron count, parser subset, XC compatibility, projector limits and the bytes actually opened. Do not trust file names. Hash and retain the exact file identity used at launch; detect changes after approval. No UPFs are bundled or committed.

Proposed curated registry stores element, potential family/version, NC type, XC, valence, relativistic/core treatment, parser compatibility, SHA-256, publisher URL, license/redistribution status, independent validation evidence and cutoff recommendations with units/source. A matching hash is identity, not scientific validation. Unknown license or absent multi-element evidence is visible; no silent Si or generic potential default. Registry lists candidates, user selects; downloads require the network grant and valid license. Existing expert functional overrides stay explicit and prominently recorded; recipes never enable them automatically.

Cutoff proposal: use the maximum validated recommended wavefunction cutoff over selected potentials; density cutoff uses the engine's documented supported relationship, not VASP ENCUT or an arbitrary QE ecutrho. Missing validated recommendations means require an explicit value and suggest a convergence scan. Mesh proposal policy `reciprocal-spacing-v1`: user chooses Δk in Å⁻¹; compute `n_i=max(1,ceil(|b_i|/Δk))` from reciprocal vectors including 2π, Gamma-centered shift zero unless explicitly selected. Explain that skew cells and anisotropy need convergence checks. Record rule/version, Δk, original metric, resulting full/reduced grid, cutoff sources and user edits. No universal Fe3O4 cutoff or mesh is asserted.

## 6. Acquisition, disambiguation and offline behavior

`StructureProvider` exposes `describe`, `search(query,cursor,fields)`, `fetch(id,revision)` and `citations`. It returns provider-scoped candidate IDs, typed units, nullable properties, per-property origins, license status, revision and web URL. Implement MP enrichment first and shared OPTIMADE structure retrieval second; native additions remain provider-specific [R2].

MP configuration uses HTTPS `https://api.materialsproject.org` and `X-API-KEY`. Formula, chemical system and contains-elements modes are distinct. Summary candidate fields and core/thermo origin joins, rate/version details and unverified raw route spelling are inventoried in [R2]. Summary entries show formula, material ID, space-group symbol/number, site count, volume, hull eV/atom, formation eV/atom, gap eV, magnetic ordering and calculation provenance, each nullable. Avoid fetching full structures for every row.

Stable candidate ordering: within one provider/release/thermo model, exclude deprecated entries by default (show count/filter); exact formula first; `is_stable=true`, then false, then unknown; hull energy ascending with unknown last; site count then opaque ID as deterministic tie breaks. Do not sort unlike providers/functionals on raw energy. Label the stable indicator `on provider convex hull`, **not** `ground state`. If a provider has a genuine ground-state flag, display its definition, scope and origin separately; otherwise display `ground-state flag: unavailable`. A hull point may be a mixture competitor, and magnetism/polymorph identification still needs explicit choice.

Proposed dialogue (no fabricated candidate values):

```text
researcher> Fe3O4
cli> Search Materials Project for exact formula Fe3O4?
     Sends formula and field list to api.materialsproject.org; no local files.
researcher> search
cli> [table from response: ID | space group | sites | hull | ordering | link]
     Choose an ID, inspect an ID, next page, refine filters, or cancel.
researcher> inspect <id>
cli> [record properties, license, structure preview, source calculation, limitations]
researcher> select <id>
cli> Selected <id> at release <release>, structure digest <digest>.
     [capability diagnostics; export-only if unsupported]
```

Even a single candidate requires explicit selection or an already explicit ID; never choose the first match automatically. Numbered selection binds to a frozen result-list digest, so pagination/reordering cannot select another record. An explicit ID skips list selection, but formula disagreement blocks and an unknown/deprecated alias needs resolution. Record both requested and canonical IDs. Distinguish zero matches, partial pages, authentication failure, license restriction, timeout and offline cache miss.

Network approval precedes the first provider request: show host, query, data to be sent, bounded page/byte budget and applicable terms. Explicit batch `--allow-network=materialsproject` authorizes that scope; absence in noninteractive mode returns `needs_network_permission` JSON rather than hanging. Scope grants persist for the displayed search/fetch flow, not arbitrary future hosts or uploads. Changing terms, provider or data disclosure requires a new grant.

Key storage: environment `MP_API_KEY`, or user config outside the repository with directory mode 0700/file 0600 on Unix (equivalent user-only ACL elsewhere). Reject insecure config and symlinks. Never accept key command-line arguments; never print, log, commit, put into request/result files, send to MCP clients or retain in REPL history. Redact HTTP errors/headers before formatting; do not dump environment variables. No key forwarding across redirects; restrict endpoints and TLS validation. No key/network: local imports and previously authorized cached records still work; display setup link or a precise cache-miss message, no synthetic structure and no unannounced alternate database.

Cache raw sanitized response bytes and normalized structures, query/field set, provider/schema/client versions, license/citations and UTC fetch time. Partition access-controlled datasets by account-local cache namespace without storing the secret. Cache keys include provider, endpoint, parameters, release and normalizer version. Bound size and eviction; pinned run inputs are never evicted as ordinary cache. Offline mode performs zero DNS/HTTP, marks cache age/release visibly and never claims current search completeness. Refresh creates a new snapshot, not an overwrite of provenance. Cross-release search/fetch detects and requests re-selection; REST cannot promise historical replay.

Result provenance stores provider, canonical/requested ID, database release (or explicit unknown), query date, retrieval date, query/selection digest, API/schema/helper versions, source task IDs/run types, source response SHA-256, structure SHA-256 and normalization version, source URL, license/citations and transformations. Hash a versioned canonical structure encoding (ordered cell/positions/species, units, PBC, finite round-trip numbers); also hash exact source bytes. This hash establishes input identity, not structural equivalence under symmetry. Reproduction uses pinned bytes rather than requerying a moving database.

Tests require no network: synthetic fixtures first; recorded provider responses only when redistribution is licensed, sanitized and attributed. Include multiple Fe3O4-like candidates without asserting real IDs, alphabetical IDs, formula mismatch, unknown ordering, release changes, missing properties, pagination, 401/403/429, malformed/oversized responses, offline hits/misses, key-leak canaries and terms refusal. An opt-in credentialed smoke test is separate from `scripts/check.sh` and never commits keys or restricted records.

## 7. HTTP and dependency boundary

| Option | Tradeoff and disposition |
| --- | --- |
| Optional Rust HTTP/TLS client | Recommended production path in CLI `providers` feature; reqwest/rustls candidate [R2, H1]. Typed errors, timeouts and response limits; default offline build excludes it. Audit crypto backend and transitives before calling the full stack pure Rust. |
| System curl subprocess | Possible bootstrap, no shell invocation. Version/TLS behavior varies across machines; keys must arrive via protected stdin configuration, not argv. Disable implicit curl config and sanitize proxy/environment behavior. Less predictable packaging, not preferred production contract. |
| Small external helper | Official mp-api can track server changes and return a bounded data-only JSON envelope. Pin helper/environment, use `MP_API_KEY` only in helper, no runtime pip install. Python/native dependencies and schema drift need packaging tests. Optional fallback, not required for local calculations. |

No client is fetched or added here. A later dependency ticket records SPDX licenses/notices for all transitives, recent maintenance/security review, exact version/MSRV, reproducible offline vendoring and arm64 macOS build/test/TLS certificate tests. New crates/manifests need owner authorization. Apply connect/read/total timeouts, decompressed byte limits, bounded redirects/retries, page ceilings and cancellation. Honor `Retry-After`; retry idempotent reads only. Numerical core has no network feature at all.

## 8. Guided execution, estimates and trust

Wizard and REPL are views over saved plans: `import/search → select → model → potentials → settings → preflight → review → run → inspect`. Commands include `back`, `explain <setting>`, `edit`, `save`, `cancel`, `show plan` and `quit`. Noninteractive mode never prompts; it returns unresolved choices. History excludes credentials, HTTP payloads and sensitive absolute paths. The chosen run is explicit, not whatever output directory was modified last.

Full preflight opens/hashes potentials, checks identity, occupations, every k basis/band capacity, resource budgets and output-path protection before solving. Estimate shows atom/electron/band/k counts, cutoff, FFT dimensions, memory **range** in MiB/GiB and time **range** in s/min plus calibration machine and uncertainty. Existing first-k array bytes are not peak all-k memory. If no calibration exists, say `runtime unknown`, show a user-specified wall-time/iteration budget and do not invent seconds. Estimates must bound preprocessing too. An expensive batch summarizes total jobs and maximum cost; optional paid-service money estimates appear only if a real tariff exists.

Before execution show resolved files, model/potential digests, defaults, changes, warnings and estimate. Confirmation binds plan/build/potential/output identities and resource limits; any change invalidates it. Batch automation supplies explicit digest-bound approval and budget; `--yes` alone cannot ignore physics errors. Cancellation is not convergence or restart support. Distinguish queued/running/converged/unconverged/failed/cancelled; publish a bundle completion manifest after payloads, with run ID and hashes, so viewers reject mixed/stale outputs. This completion contract is a planned CLI prerequisite, not an existing guarantee.

## 9. Agent interface and schema contract

Prefer a local stdio MCP server wrapping the same CLI operations; no inbound network listener by default. Tool inventory: `capabilities`, `import_structure`, `search_materials`, `select_material`, `compile_request`, `preflight`, `run_plan`, `run_status`, `cancel_run`, `list_results`, `render_view`. No generic shell tool. Search requires a network grant, run requires approval digest/budget, and render only reads approved output paths. MP's existing MCP server may supply data, but its output must enter the same import/provenance gate [R2, A1–A2].

Publish JSON Schemas for input, resolved plan, translation report, capabilities, result manifest and view request. This ticket does not edit frozen `docs/schemas/`; UI-01 owns the future explicit authorization. Schema version negotiation is independent of binary version. Strict request objects reject unknown keys; old result readers tolerate documented additive fields but reject incompatible major versions. Typed unavailable quantities carry a reason, never zero. Binary grids remain standard artifacts, not enormous JSON arrays.

Proposed contract fragments (field inventory, not an implemented complete schema):

| Object | Required content |
| --- | --- |
| Capability snapshot | schema version, binary/build fingerprint, supported models/actions/formats, constraints, evidence and roadmap pointers |
| Calculation input | structure reference/digest, explicit XC/occupations, potential map, cutoff with unit, k specification, SCF controls, required observables; native TOML remains executable authority |
| Plan | schema version, request ID, deterministic plan digest, resolved input hashes, capabilities digest, blockers, recommendations, cost estimate and approval state |
| Result | schema version, request/run/model IDs, status, engine result reference, energy units/reference, provenance, completion manifest and artifact hashes |
| Artifact | relative path, MIME/format/version, size, SHA-256, units, availability/reason; no key or secret URLs |
| Tool error | stable code, stage, field path, message, required choice/capability, roadmap link, retryability |

A capability example is `{"name":"spin_polarized","status":"unsupported","roadmap":"docs/roadmap/esle.md"}`; the definitive explanation points to phase E in that document. Schemas must test the current legacy Gamma and all-k/PBE variants, not assume top-level eigenvalues represent every k point. MCP protocol version is negotiated; the reviewed 2025-06-18 schema features are a baseline, not a latest-version assertion.

The following minimal JSON Schema 2020-12 envelopes illustrate the proposed machine boundary. They are design examples, not substitutes for UI-01's full native-input/results schemas. The input references a local TOML file so agents cannot accidentally invent a competing physical-input dialect; the resolved plan validates that file against the full native-input schema. The result references the engine record and typed artifact metadata; the full result schema adds all numerical fields/units and provenance described above.

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "title": "Proposed calculation input envelope v1",
  "type": "object", "additionalProperties": false,
  "required": ["schema_version", "input_path", "input_sha256", "required_capabilities"],
  "properties": {
    "schema_version": {"const": 1},
    "input_path": {"type": "string", "minLength": 1},
    "input_sha256": {"type": "string", "pattern": "^[0-9a-f]{64}$"},
    "required_capabilities": {"type": "array", "uniqueItems": true,
      "items": {"type": "string", "minLength": 1}}
  }
}
```

```json
{
  "$schema": "https://json-schema.org/draft/2020-12/schema",
  "title": "Proposed calculation result envelope v1",
  "type": "object", "additionalProperties": false,
  "required": ["schema_version", "run_id", "plan_sha256", "status", "artifacts", "diagnostics"],
  "properties": {
    "schema_version": {"const": 1},
    "run_id": {"type": "string", "minLength": 1},
    "plan_sha256": {"type": "string", "pattern": "^[0-9a-f]{64}$"},
    "status": {"enum": ["converged", "unconverged", "failed", "cancelled"]},
    "artifacts": {"type": "array", "items": {
      "type": "object", "additionalProperties": false,
      "required": ["path", "format", "sha256"],
      "properties": {
        "path": {"type": "string", "minLength": 1},
        "format": {"type": "string", "minLength": 1},
        "sha256": {"type": "string", "pattern": "^[0-9a-f]{64}$"}
      }
    }},
    "diagnostics": {"type": "array", "items": {"type": "string"}}
  }
}
```

Path-root enforcement, digest verification, completion status and required scientific metadata are semantic validators in addition to schema validation. A syntactically valid JSON object is not permission to run or proof of convergence.

LLM output is untrusted proposed input. File/database comments are data, never instructions. The server independently enforces schemas, file roots, capabilities, grants and budgets; a model cannot self-issue approval. Do not send private structures/prompts to an external LLM from the CLI. A future opt-in mode would need separate provider/key/disclosure/cost design.

## 10. Link-based results and visualization

Print material-page HTTPS links, canonical result-folder `file://` links, `report.html` file links, and an optional `http://127.0.0.1:<port>/...` viewer link. OSC 8 framing follows public terminal documentation [R2, V1–V2]; provide `--hyperlinks=auto|always|never`, TTY detection and plain URL/path fallback in logs. Strip control characters, URI-encode paths and allow only approved schemes. Never put credentials into a URL. For SSH, explain that a file link points to the remote host and provide download/port-forward instructions; do not silently publish results.

First offline report: one self-contained HTML file with inline CSS, SVG, small project-owned JavaScript/Canvas, inline sanitized data and no CDN/fonts/analytics/fetch at runtime. No embedded third-party library: repository MIT licensing applies to renderer code; source data retain their own licenses. Use SVG for bands, DOS, SCF residual/energy history, energy-component bars and a static projected cell/atom view; Canvas adds rotation and density slices. No guessed chemical bonds, no crystallographic analysis in the renderer. Show a static/table fallback when scripts are disabled. This bounded technique is chosen for inspectability, portable exports and small offline artifacts; advanced isosurfaces remain external-viewer work initially.

Report panels: structure/cell/species and source provenance; capability/model/potential summary; convergence certificate and iteration history; the six additive energy components (kinetic/local/nonlocal/Hartree/XC/Ewald) separate from diagnostics; bands/DOS only when those artifacts exist; density slices only when valid volumetric data exist. Never call an SCF k-point eigenvalue list a computed path band structure. Failed/unconverged data carry a persistent banner. Missing outputs display `unavailable` with reason. Units and energy zero (P8A insulating VBM plus absolute offset, no invented Fermi level), valence-only density, potential G=0 convention and sample/grid details are visible.

Standard files remain primary: cube/XSF for grids; additional CHGCAR/LOCPOT and VTI writers only under separate validated tickets; bands/DOS CSV, plain `.dat`, xmgrace `.agr`, QE-style `bands.dat.gnu`/`dos.dat`, VASP-style EIGENVAL/DOSCAR plus JSON metadata and SVG. Structure export prioritizes POSCAR/CIF/XYZ/extxyz/XSF. Preserve what each format cannot carry in a sidecar. Every writer needs an independent reader or specification test; the HTML renderer does not bypass P8A/P8B gates.

`plot bands` resolves bands CSV/metadata; `show the density slice` resolves cube/XSF plus metadata into the same versioned view request. Record plane/origin/offset, coordinate convention, interpolation, color range and units. A default report may embed bounded preview grids; label downsampling and retain the full standard artifact separately. The embedded preview is self-contained offline even if the original bundle is absent; large-grid full exploration requires the local bundle. Export plots as SVG and Canvas views as PNG where supported; browser print-to-PDF is optional. Reproducibility concerns plot data/view parameters, not cross-browser pixel identity.

Optional viewer binds loopback only, read-only bundle root, no directory traversal/symlink escape, Host/Origin checks, no CORS or run endpoints. Use an ephemeral access token kept out of persisted reports/logs and an explicit lifetime; avoid dumping the full home directory. Do not open a browser or bind a remote address automatically. Content Security Policy denies connections and third-party resources; escape HTML, JSON script terminators and user labels, with hashes for permitted inline scripts. Offline acceptance intercepts requests and asserts **zero external requests**; injected HTML/OSC sequences remain inert.

## 11. Limits and verification plan

This is design work. No user study establishes the grammar or ranking quality, no physics capability is added, and no performance estimate is measured here. API/license verification limits are explicit in [R2]. Contract tests should cover translation coverage, bilingual equivalence and ambiguity, blocked magnetism/+U/metal/relaxation, identity mismatches, repeated-plan determinism, changed-input approval invalidation, offline fixtures, schema evolution, interrupted manifests and malicious labels/paths. Heavy scientific tests stay ignored in debug with release commands; routine gates remain fast.

Only the three P9A documents are owned. This task therefore does not update dashboard/session JSON, README figures, ticket files, dependency manifests or schemas. Roadmap entries are proposed scopes, not created/registered tasks. Run `scripts/check.sh` and report its actual exit code; a failure in the existing build is not evidence that any proposed UI works or fails.

## References

- R1. rust-dft, *Current architecture*, *CLI*, ADRs 0012–0015/0017, *feature-gap survey*, P8A/P8B tickets and *ESLE roadmap*, 2026; local paths: [architecture](../architecture.md), [CLI](../cli.md), [ADR 0017](../decisions/0017-coupling-commitments.md), [survey](../research/feature-gap-survey.md), [ESLE](../roadmap/esle.md). Inspected 2026-11-03.
- R2. *User-interface ecosystem: evidence and choices*, 2026, [dated primary references M1–M11, D1–D6, F1–F9, V1–V2, A1–A2, H1](../research/ui-landscape.md#references), all accessed 2026-11-03. External factual claims above use these entries; architecture mechanisms and numerical suggestion rules are proposals.
