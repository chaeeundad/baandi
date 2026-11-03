# User-interface implementation roadmap

| Item | Proposal |
| --- | --- |
| Status | P9A design only; the fifteen rows below are candidate tickets, not created or approved tasks |
| Largest immediate value | Import familiar structures with a complete translation report; require P8G preflight before executing imported inputs |
| Sequence | Imports → safe preflight → recipes → MP → guide → agent interface → offline report |
| Effort | 41–60 agent-days for the fifteen scopes; planning estimates, not elapsed-time promises or measurements |
| External dependencies | P8G species identity, P8A bands/DOS, P8B volumetric output; dependency and frozen-path approvals as listed |
| User decisions | Provider scope, external-agent-only language mode, licensing policy, renderer scope; recommendations below |
| Source | [Architecture](../architecture/user-interface.md), [ecosystem evidence](../research/ui-landscape.md) |

## Delivery principles and quick wins

Finish import/preflight before automatic acquisition-to-run. A structure import saves users from private lattice tags immediately; a chatbot without species checks can make the wrong material easier to run. Capability truth is a release criterion for every frontend. All fifteen tickets implement adapters/contracts, not new physics. Large tests stay outside the fast debug gate, with documented release-mode commands.

Quick wins within the scopes below: print header element alongside species alias through P8G; correct stale Gamma-only/LDA help in UI-01; add raw material/output links with terminal fallback in UI-10/UI-14; provide `explain` text and printable defaults in UI-08. Each is reviewable in approximately 0.5–1 agent-day, already included in row estimates. Do not schedule competing edits to CLI dispatch/config/run; assign one wiring owner per merge window.

## First fifteen implementation tickets

Every row requires `scripts/check.sh` plus its listed gate, real exit codes and stored evidence. Proposed owned areas require exact file lists before delegation; mentioning a frozen file here is not authorization to edit it now. Existing tickets own their work; do not recreate P8G or claim it delivered.

| Order / proposed ID | Title and owned areas | Acceptance gate | Dependencies | Agent-days |
| --- | --- | --- | --- | --- |
| 1 / UI-01 | Capability and frontend contracts. CLI `capabilities.rs`, dispatch/tests; proposed input/plan/result/capability schemas and CLI docs, with explicit frozen-schema authorization | Machine inventory matches actual build; legacy Gamma/all-k/PBE result variants validate; unsupported spin/+U/metals/forces advertised accurately; schema version/refusal tests; stale help corrected | Architecture acceptance; schema owner review | 2–3 |
| 2 / UI-02 | Translation report and import scaffold. CLI `import/mod.rs`, `translation.rs`, import tests/docs | Every parsed source item assigned mapped/defaulted/ignored/unsupported; unknown tags block; deterministic JSON/Markdown report and source digest; geometry-only extraction preserves lost-intent report | UI-01 | 2–3 |
| 3 / UI-03 | POSCAR/CONTCAR import/export. CLI `import/poscar.rs`, writer/tests/how-to | Independent specification fixtures cover negative/positive scaling, species order, skew cells, Cartesian/direct and absent elements; POTCAR never opened; selective flags inventoried; independent export readback | UI-02; P8G before execution | 2–3 |
| 4 / UI-04 | CIF import/export. CLI/external format adapter, fixture tests, dependency inventory if needed | Symmetry-expanded atom count/geometry independently checked; block/frame selection; occupancy/disorder refusal; malicious input bounds; CIF export specification checks | UI-02; parser license/maintenance/arm64 review | 3–5 |
| 5 / UI-05 | XYZ/extxyz and XSF import/export. CLI `import/{xyz,xsf}.rs`, writers/tests/docs | Explicit cell/PBC/unit enforcement; multiple frames require choice; skew cell and lattice-order fixtures; unavailable dimensions refused; independent readers/specification checks, including XSF public-document recheck | UI-02; coordinate contracts; P8G before execution | 3–4 |
| 6 / UI-06 | Foreign calculation input translation. CLI `import/{qe,vasp}.rs`, documentation fixtures and import guide | QE ibrav=0/units/cards and VASP INCAR/KPOINTS mapped by documented semantics; actual weighted k sets checked for parity/shifts; unsupported physical tags and unknowns block; no POTCAR access or numerical-equivalence claim | UI-03, UI-02; ADR 0013 | 4–6 |
| 7 / UI-07 | Data-only ASE/pymatgen JSON adapters. External `adapters/structures/` or CLI JSON module, tests/docs | Dialect/version/record selection; no class execution/pickle; cell/species/occupancy checks; constraints and calculator outputs accounted for; parser dependency license gate | UI-02; UI-04/05 for cross-format fixtures | 2–3 |
| 8 / UI-08 | P8G integration, potential registry metadata and recommendation policy. CLI preflight presentation/registry/estimate modules and tests; no numerical parser changes in this scope | P8G full identity checks used by all frontends; Fe label with Si UPF fails; aliases with explicit correct element pass; registry hashes/licenses and XC checked; all-k memory estimate bounds stated; cutoff/Δk rules/version and unknown runtime shown | Existing P8G complete; UI-01/02; registry licenses and evidence | 3–4 |
| 9 / UI-09 | Versioned recipes and bilingual deterministic request compiler. CLI `request.rs`, `recipes.rs`, fixtures/docs | English/Korean examples produce same typed intent; unresolved polymorph/standard path prompts; unsupported relax/magnetism/+U blocked; exact editable TOML/defaults; plan digest stable; edits invalidate approval | UI-08; P8A only for executable bands/DOS recipes | 3–4 |
| 10 / UI-10 | MP provider, credentials, cache and selection. Optional CLI provider module, offline fixtures, dependency record and guide | Current authenticated OpenAPI/terms verified; formula/chemsys/elements semantics; sorted stable lists without ground-state invention; opaque IDs; release drift; safe header/key redaction; offline zero-network; licensed fixtures; no unattended first-match runs | UI-08/09; approved HTTP dependency; account access/terms verification | 4–6 |
| 11 / UI-11 | OPTIMADE provider and provider policies. CLI provider modules, provider conformance fixtures/docs | `/info` discovery/version negotiation; pagination/filter distinction; units/PBC/unknown properties and provider licenses retained; one verified COD or OQMD deployment; no cross-provider energy ranking | UI-10 common provider/cache contract; provider-specific access approval | 3–4 |
| 12 / UI-12 | Guided CLI and run approval. CLI `guide.rs`, dispatch/state-machine tests, guide documentation | Search/import/select/review/save/cancel/back flow; noninteractive structured unresolved choices; network/run grants distinct; digest-bound budget approval; no unsupported-physics override; interrupted output detection | UI-09/10; P8G; result completion contract integration from UI-01 | 3–4 |
| 13 / UI-13 | Agent JSON tools and stdio MCP. External `adapters/mcp/`, CLI JSON wrappers, schema conformance/security tests | Tools expose input/output schemas; same gates as guide; no generic shell; agent cannot self-authorize; key/path sanitization; typed unavailable results; cancellation never claims restart | UI-01/12; MCP version/dependency review | 2–4 |
| 14 / UI-14 | Offline HTML report and terminal links. CLI report renderer/templates/assets/tests, report guide | One offline file, no CDN or external requests; sanitized OSC/HTML; plots from standard artifacts; model/provenance/license/convergence visible; static fallback; failed/mixed bundle rejection; SVG export | UI-01 completion manifest; P8A/P8B for corresponding optional panels, not basic report | 3–4 |
| 15 / UI-15 | Conversational views and optional local server. CLI view dispatcher, external loopback server, browser tests/docs | Density slice/plot-bands resolve same artifact IDs; no SCF on view; units/plane/downsampling logged; SVG/PNG export; loopback/root/Origin restrictions; missing outputs explicit | UI-09/12/14; P8A/P8B; optional server review | 2–3 |

The estimates sum to **41–60 agent-days** when every row is included. This row-wise sum is authoritative; it includes fixture work and routine review corrections but excludes new physics, credential procurement, dependency acquisition by the architect, independent scientific reference runs and waiting for approvals. Calendar duration depends on sequencing and available agents; no delegation is started by this document.

## Milestones and handoff gates

| Milestone | Deliverable | Exit condition |
| --- | --- | --- |
| M1: trusted imports | UI-01–08, coordinated with P8G | Researcher can import standard geometry, see every loss/default and prove species identity before run; local-only use works without Python/network dependencies unless an optional format needs them |
| M2: reproducible setup | UI-09–11 | Same request/selection/versions yields the same plan; MP disambiguation has authenticated fixtures/terms clearance; offline replay works; provider metadata does not substitute for engine validation |
| M3: interactive/agent control | UI-12–13 | Guide and agent invoke one approval/preflight path; blocked models remain blocked; run identity and completion manifests prevent stale-result success |
| M4: understandable results | UI-14–15 with P8A/P8B | Downloadable standard outputs and one offline report; view commands consume existing results; no external requests, secret exposure or unsupported quantity fabricated |

Structure-only imports may land before UI-08, but must not expose an automatic run shortcut until P8G is merged. UI-08 integrates P8G rather than implementing a second identity checker. UI-01 defines the completion contract; UI-12 implements the CLI manifest lifecycle before any guide/agent treats bundle presence as success. UI-14 can render existing bundles only with explicit legacy status caveats until that lifecycle gate passes.

P8A/P8B are parallel prerequisites, not superseded by this roadmap. Coordinate output schema/version/units and the project's additional standard-format policy; a panel remains unavailable until its producer is validated. P8A's initial limited path presets do not authorize a universal “standard path” claim. New symmetry/path libraries require their own review. No UI ticket may silently add smearing, force expressions, occupation control or potential conversion.

## Decisions for the user before implementation (not required to finish P9A)

| Decision | Options | One recommendation and consequence |
| --- | --- | --- |
| D1: first database | MP alone; MP plus COD/OQMD; all six sources | **MP first, then one COD/OQMD OPTIMADE deployment.** Matches the requested workflow while keeping licensing and provenance tractable; defer AFLOW/NOMAD/JARVIS enrichment. |
| D2: unrestricted language | External agent interface only; an opt-in LLM service inside CLI | **External agent only.** CLI retains deterministic bilingual grammar and guide; no private structures/prompts leave the machine for inference through this CLI. |
| D3: data license policy | Permit reviewed noncommercial datasets; require commercial-compatible licenses only | **Commercial-compatible sources by default; restricted datasets explicit opt-in only after their actual terms are established.** Unknown licenses block acquisition/redistribution; MP GNoME is not blanket-authorized. |
| D4: viewer scope | Offline SVG/Canvas first; immediate embedded full 3D/isurface library | **SVG/Canvas first**, with external cube/XSF viewer links for advanced isosurfaces. No third-party JS dependency until justified and license-cleared. |
| D5: network implementation | Optional Rust HTTP/TLS; system curl; pinned mp-api helper | **Optional Rust HTTP/TLS after OpenAPI fixtures**, retaining helper as a compatibility fallback. Architect approves versions/transitive licenses and arm64 checks first. |
| D6: expected scientific target | Current nonmagnetic insulators first; prioritize Fe3O4 physics | **Deliver current-insulator workflows first.** Fe3O4 search/export remains useful, but executable magnetic/relaxation workflows await separately validated physics tickets. |

These decisions are design review items, not a request for permission to alter extra files in P9A. No implementation ticket is created in this task. At creation, the orchestrator must register the actual ID/title/ownership and dashboard/session records under its own authorized scope.

## Deferred physics and larger product gates

- **ESLE phase E:** spin polarization, metallic occupations/smearing, +U and occupation control. This is an existing roadmap phase, not a new UI ticket or delivery date. Each method needs separate physical validation before capability advertisement.
- **PHYS-RELAX (placeholder, not created):** forces, then optimizer and stress/cell relaxation under separate ADRs and tests. A parsed `relax` request remains blocked until these exist.
- **FORMAT-EXTRA (placeholder):** CHGCAR/LOCPOT, VTI, EIGENVAL/DOSCAR, AGR and other policy-required writers not delivered by the first P8 producers; independent readers/specification tests and loss documentation required. Do not advertise formats solely because this roadmap lists them.
- **REPORT-3D (placeholder):** advanced periodic bonds/isosurfaces with a reviewed library, GPU/browser fallback and explicit licenses; no CDN at runtime.
- Remote hosted reports, uploads/sharing, schedulers, billing, plugin marketplaces and LLM-provider keys require separate trust/privacy scopes. Local links do not imply public hosting authorization.

## Evidence gates, uncertainty and limits

Before UI-10 ships, an authorized maintainer must verify modern MP `X-API-KEY`, raw summary/core/thermo routes, actual field availability and query pagination against the current OpenAPI, obtain legal/attribution terms and capture licensed fixtures. The P9A investigation found MP public query documentation but received 403 on terms/OpenAPI. Full current key terms and blanket data licensing are unverified. Dataset permission is not inferred from a key that happens to work. API schema and database release are independently recorded.

XSF public documentation retrieval failed; UI-05 cannot treat a private round trip as independent format validation. JARVIS per-dataset licenses and provider live uptime were not established. Estimates and order are engineering judgment; there has been no user study, browser test or working UI prototype. P9A changes no feature code, generated evidence, dashboard, README, manifests or schemas, and claims no new physical numerical result.

Review each delivery against the [architecture](../architecture/user-interface.md): complete translation accounting, explicit model, immutable identity, reproducible plan, no key leakage, no unapproved networking, no silent approximation and standard outputs first. Report actual checks and scientific evidence separately from synthetic adapter tests. Follow the README-generation and dashboard rules only in implementation tickets that authorize those paths; P9A's three-file limit takes precedence here.

## References

- rust-dft, *Researcher-facing interface architecture*, 2026, [design and contracts](../architecture/user-interface.md), inspected 2026-11-03.
- rust-dft, *User-interface ecosystem: evidence and choices*, 2026, [primary sources with title/year/URL/access dates](../research/ui-landscape.md#references), accessed 2026-11-03. Ecosystem-dependent choices above rely on those sources; estimates are opinion.
- rust-dft, *Practical feature and workflow gap survey*, 2026, [verified preflight/identity gaps](../research/feature-gap-survey.md), inspected 2026-11-03.
- rust-dft, *ESLE roadmap*, 2026, [physics phases](esle.md); *Coupling commitments*, ADR 0017, 2026, [adapter boundary](../decisions/0017-coupling-commitments.md); inspected 2026-11-03.
- rust-dft, *P8A-BANDS-DOS* and *P8B-VOLUMETRIC*, 2026, [bands/DOS ticket](../../tickets/P8A-BANDS-DOS.md), [volumetric ticket](../../tickets/P8B-VOLUMETRIC.md), inspected 2026-11-03.
