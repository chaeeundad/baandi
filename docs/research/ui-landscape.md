# User-interface ecosystem: evidence and choices

| Item | Finding as of 2026-11-03 |
| --- | --- |
| Status | P9A design research; public documentation only, no upstream numerical source inspected |
| Recommendation | Standard-file import first; MP-specific metadata adapter plus a shared OPTIMADE structure adapter |
| Main qualification | Database structures do not establish compatibility with this engine's physics or pseudopotentials |
| Current MP | Modern API, not legacy `/rest/v2`; latest documented release v2026.09.28 |
| Verification boundary | No authenticated queries; MP OpenAPI/terms pages returned 403; some license and wire-schema details remain release gates |
| Companion documents | [Architecture](../architecture/user-interface.md), [first fifteen tickets](../roadmap/ui.md) |

## Evidence method

Access date for every external reference below: **2026-11-03**. “Undated” means the page does not provide a reliable publication year; a documentation copyright year is not treated as a release date. Descriptions are paraphrases. Search snippets from third parties are not used to establish licensing. No QE, VASP, ASE or other GPL/proprietary implementation was opened, cloned or reproduced; public format manuals and protocol documentation are the basis for future importers. This is a design, not a claim of implemented interoperability.

Repository evidence: `crates/dft-cli/src/{main,config,run}.rs`, examples, ADRs 0012–0015 and 0017, and the [feature-gap survey](feature-gap-survey.md). The checked configuration accepts LDA-PZ and PBE, and Gamma, meshes and weighted k points. Older architecture/help prose still says Gamma-only LDA. Existing `validate` does not read UPFs. P8A bands/DOS and P8B volumetric output are tickets, not available commands in this snapshot. The species-alias/UPF-element mismatch documented in the survey is a prerequisite blocker for automatic setup.

## Materials Project: documented contract and remaining checks

| Concern | Evidence and implication |
| --- | --- |
| Base and credentials | `https://api.materialsproject.org`; account API key, sent as the `X-API-KEY` HTTP header over HTTPS. Header names are case-insensitive. Modern client reference documents `MP_API_KEY` [M1, M2]. Header usage is also documented by MP's older API manual [M8]; live modern OpenAPI confirmation remains necessary. Never use legacy URL query-key authentication. |
| Routes | Modern client `materials.summary.search` and `materials.search` [M3] correspond to the intended REST `/materials/summary/` and `/materials/core/` routes. `/materials/thermo/` supplies related calculation provenance. Raw routes, pagination and query parameter spelling must be checked against OpenAPI before shipping [M4]; `/materials/` must not be assumed to be the core search route. |
| Formula search | Public client `formula="Fe3O4"` searches composition; `chemsys="Fe-O"` is an exact chemical system; `elements=["Fe","O"]` is a contains-elements query and may include other elements [M5]. Expose these as different search modes; never broaden a formula silently. |
| Fields | Summary inventory should request `material_id`, `formula_pretty`, `symmetry`, `nsites`, `volume`, `energy_above_hull`, `formation_energy_per_atom`, `is_stable`, `band_gap`, `is_metal`, `is_magnetic`, `ordering`, `structure`, `origins`, `deprecated`, and available license metadata. Validate availability against the server schema; unknown/missing is not false or zero [M3, M5]. |
| Units and origin | Hull and formation energies are eV/atom; gap eV; volume Å³; site count integer. Space group comes from `symmetry.symbol`/`number`. Summary properties can originate in different calculations. Resolve structure `origins.task_id`, then the matching thermo entry `run_type`; core `calc_types` is a task-to-calculation-type map, not a single material-wide functional [M3, M6]. |
| Versions | Latest changelog entry v2026.09.28 went live 2026-10-19; client >=0.47.0 is advised for compatibility. REST always serves the current database; client `db_version` selects historical S3 data, not a historical REST snapshot [M2, M6]. Record version at query/fetch time and preserve bytes. |
| IDs | MP introduced alphabetical IDs in 2026; historical numeric IDs remain aliases [M6]. Store opaque provider IDs and resolved aliases, never parse IDs as integers or restrict them to `mp-[0-9]+`. |
| Rate limits | Documentation says limits start at 25 requests/s [M7], not a guaranteed per-user entitlement. Batch IDs, select fields, cache results. Proposed client policy: at most 2 requests/s, bounded concurrency, honor `Retry-After` and bounded backoff on 429/temporary failures. |
| Web links | MP recommends its main domain; next-gen points to its newest site [M9]. Display `https://next-gen.materialsproject.org/materials/<encoded-id>`; keep URL construction provider-owned and configurable. No claim that the example mp-19306 was retrieved. |

The public client interface is the verified query description; it is **not** a guarantee that Python keyword arguments can be pasted into an HTTP URL. Candidate REST requests such as `/materials/summary/?formula=Fe3O4&_fields=material_id,formula_pretty` are design sketches pending the OpenAPI gate. Fixture capture must record the exact successful route, parameters, schema digest and response metadata. Avoid downloading all structures merely to display a candidate list: summary first, selected structure second.

### License, attribution and API-key terms

MP's official terms page [M10] could not be fetched (403), so the full current data license, exceptions and key-specific contractual terms were **not independently verified**. Do not advertise blanket commercial reuse or describe every MP record as CC BY 4.0. The release notes explicitly identify GNoME-originated structures as BY-NC and require explicit acceptance for access [M6]. The implementation must retrieve/verify the applicable record/dataset license and terms through an authorized browser before enabling redistribution or fixture publication. Unknown license means private local use is not automatically authorized either: show the terms link and block the adapter's download until the user can establish applicable permission. No automatic acceptance of provider terms.

The older official MP manual [M8] says keys identify the account, must remain secret, must not be shared, and the account holder is responsible for use; these are conservative handling requirements, not proof that all legacy contractual restrictions still apply. Recheck current terms at release. Never ship a shared project key.

For attribution, retain provider, record link, dataset/version, originating paper/DOI, applicable license URL, and a list of transformations. When CC BY applies, preserve attribution, link the license and indicate changes; verify the precise license version before distribution. MP's 2025 staff citation announcement requests the new Horton et al. Nature Materials paper for current website/API data, DOI `10.1038/s41563-025-02272-0`, and the Jain et al. 2013 paper for legacy data [M11]. Include property-specific citations returned by the provider. A citation is not a license grant.

## Other structure sources and why one adapter is insufficient

| Source | Verified interface/value | License/access evidence and design consequence |
| --- | --- | --- |
| COD | Public crystal structure collection, CIF download/search; experimental structures complement calculated structures [D1] | CC0/public-domain dedication; acknowledge original structural authors. Handle occupancy, disorder and measurement conditions, not fabricated hull energies. |
| OQMD | Native REST and OPTIMADE, anonymous queries; calculated formation-energy data [D2] | Official API page states CC BY 4.0 and requests citations. Prefer shared OPTIMADE geometry; optional native enrichment. |
| AFLOW | AFLUX property queries and AFLOW REST, JSON search results [D3] | Official documentation restricts data to scientific, academic and noncommercial purposes. Do not infer a permissive data license from software availability. Provider opt-in and terms review required. |
| NOMAD | Repository with rich calculation metadata and APIs [D4] | Official site states CC BY 4.0 for shared data. Private/unpublished access differs; preserve entry-level provenance and access rules. Avoid mixing calculations just because structures look alike. |
| JARVIS | NIST documents Figshare datasets through JARVIS-Tools, including DFT collections and distinct methods [D5] | Dataset-specific license and live REST authentication were not established here. Start with explicitly licensed downloaded datasets and dataset IDs; do not label every mirrored collection with one license. |
| OPTIMADE | Versioned structures API, filters, info discovery, pagination and provider registry [D6] | Common transport/geometry reduces implementation duplication. Provider extensions, missing values, authentication, licensing and properties still differ. Discover `/info` and `/info/structures`; negotiate supported versions. |

Recommendation: a common `StructureProvider` contract with native MP enrichment and an OPTIMADE implementation plus small provider policies. OPTIMADE can cover baseline structure search across compatible deployments; it does **not** standardize every hull energy, magnetic label, calculation provenance, citation, authentication policy or web link. No provider's current uptime/conformance was live-tested. Verify each deployment and advertised fields before enabling it. Filter examples: `chemical_formula_reduced="Fe3O4"` (provider-normalized spelling) versus `elements HAS ALL "Fe", "O"`; exact element system additionally constrains `nelements=2`. Preserve pagination and unknown values. Do not merge IDs across providers; geometrical deduplication is an optional annotated presentation, never provenance erasure.

## Imports and interaction ecosystem

Public VASP POSCAR/KPOINTS/INCAR documentation [F1–F3] supports independent input translation, not numerical equivalence. POTCAR supplies potential identity in VASP [F4]; this project must never open, inspect, hash, import or redistribute it. A POSCAR without explicit element names needs a user-supplied identity map. QE's publicly served `pw.x` input description currently identifies itself as version 7.5 [F5]; pin importer support to tested documentation versions, not an assumed latest executable. Read documentation only.

CIF has a formal IUCr specification [F6]. ASE documents extxyz and JSON formats [F7]; pymatgen documents JSON serialization [F8]. Neither means all JSON is a standard structure format. Specify supported dialect/version, reject arbitrary class deserialization and require an explicit record/frame selector. The XSF official page returned 502 [F9]; P8B supplies repository conventions, but external format verification is still an importer release gate.

OSC 8 hyperlinks are implemented by terminals including WezTerm and Contour [V1, V2]; support is not universal. Plain URLs remain necessary. For the first offline report, use project-owned SVG and small Canvas interaction (MIT repository license), with no third-party embedded library. This is a design preference: static vector plots plus projected cell/atom rotation are sufficient for the first report and avoid a large rendering dependency. Dedicated isosurfaces can use external viewers reading cube/XSF; 3Dmol is a future candidate, not selected or license-cleared in this task.

MCP tools expose input/output schemas [A1]. MP itself documents an MCP server since mp-api 0.46.0 [A2]. Reuse external assistant capabilities for unrestricted language; do not build a second general LLM stack inside the numerical engine. MP's MCP does not replace this engine's preflight, capability negotiation or run approval.

For HTTP, reqwest documents Rust TLS support and MIT/Apache-2.0 licensing [H1]. This is a candidate, not an approved dependency: pin a version and audit all TLS/crypto transitives, native build requirements, MSRV, security maintenance and Apple Silicon offline build before adoption. “Rust HTTP client” does not establish a wholly pure-Rust crypto dependency tree.

## Limits and opinions

No authenticated MP/other-provider query, rate-limit experiment, browser compatibility trial, importer equivalence test, or dependency installation was performed. API responses, license metadata keys, magnetic ordering completeness and record-specific ground-state meaning need fixtures. No real Fe3O4 candidate table is claimed. Current MP terms and OpenAPI are blocked; XSF documentation retrieval failed. JARVIS dataset licenses remain per-dataset checks. Ranking, interface wording, ticket estimates and the SVG/Canvas choice are recommendations, not ecosystem facts or measured usability results.

## References

All accessed **2026-11-03**. Years below are publication years when known; otherwise “undated, current documentation”.

- M1. Materials Project, *Getting Started*, undated: https://docs.materialsproject.org/downloading-data/using-the-api/getting-started
- M2. Materials Project, *MPRester API reference*, current documentation (site banner says 2022): https://materialsproject.github.io/api/_autosummary/mp_api.client.mprester.MPRester.html
- M3. Materials Project, *Querying Data*, undated: https://docs.materialsproject.org/downloading-data/using-the-api/querying-data
- M4. Materials Project, *Advanced Usage*, undated: https://docs.materialsproject.org/downloading-data/using-the-api/advanced-usage ; linked OpenAPI https://api.materialsproject.org/docs (403).
- M5. Materials Project, *SummaryRester API reference*, current documentation (site banner says 2022): https://materialsproject.github.io/api/_autosummary/mp_api.client.routes.materials.summary.SummaryRester.html
- M6. Materials Project, *Database Versions*, 2026 release entries: https://docs.materialsproject.org/changes/database-versions
- M7. Materials Project, *Tips for Large Downloads*, undated: https://docs.materialsproject.org/downloading-data/using-the-api/tips-for-large-downloads
- M8. Materials Project, *The Materials API*, legacy documentation, undated (documentation prose only): https://github.com/materialsproject/docs/blob/master/docs/open-apis/the-materials-api.md
- M9. Materials Project, *Changes and Updates*, undated: https://docs.materialsproject.org/changes
- M10. Materials Project, *Terms of Use*, undated, inaccessible (403): https://next-gen.materialsproject.org/about/terms and https://materialsproject.org/about/terms
- M11. Aaron Kaplan (MP), *Materials Project Paper / Citation*, 2025: https://matsci.org/t/materials-project-paper-citation/64422
- D1. COD, *Crystallography Open Database*, undated/current: https://www.crystallography.net/cod/
- D2. OQMD, *OQMD API*, undated: https://www.oqmd.org/api/ ; *RESTful API*, documentation v1.2.0: https://static.oqmd.org/static/docs/restful.html
- D3. AFLOW, *The Aflux Search API*, undated; cites 2017 AFLUX paper: https://www.aflow.org/documentation/
- D4. NOMAD, *Materials science data, managed and shared*, undated: https://nomad-lab.eu/nomad-lab/ ; *Upload and share with the API*, undated: https://docs.nomad-lab.eu/tutorial/upload_publish_api.html
- D5. NIST, *Databases — JARVIS-Tools*, undated: https://pages.nist.gov/jarvis/databases/
- D6. OPTIMADE, *API specification v1.3.0*, undated release page: https://www.optimade.org/specification/latest/
- F1. VASP, *POSCAR*, undated/current wiki: https://vasp.at/wiki/index.php/POSCAR
- F2. VASP, *KPOINTS*, undated/current wiki: https://vasp.at/wiki/index.php/KPOINTS
- F3. VASP, *INCAR*, undated/current wiki: https://vasp.at/wiki/index.php/INCAR
- F4. VASP, *POTCAR*, undated/current public format description: https://vasp.at/wiki/POTCAR
- F5. Quantum ESPRESSO, *pw.x Input File Description*, version 7.5, year not stated: https://www.quantum-espresso.org/Doc/INPUT_PW.html
- F6. IUCr, *CIF Version 1.1 Working specification*, undated page: https://www.iucr.org/what-we-do/digital-standards/cif/cif1
- F7. ASE, *Format Specific Options*, undated/current documentation: https://docs.ase-lib.org/ase/io/formatoptions.html
- F8. pymatgen, *Usage*, undated/current documentation: https://pymatgen.org/usage.html
- F9. XCrySDen, *XSF format*, undated, retrieval failed (502): https://www.xcrysden.org/doc/XSF.html
- V1. WezTerm, *Use hyperlinks directly in the terminal*, undated: https://wezterm.org/recipes/hyperlinks.html
- V2. Contour, *Clickable Links*, undated: https://contour-terminal.org/vt-extensions/clickable-links/
- A1. MCP, *Tools*, specification 2025-06-18 (explicit baseline, not a latest-version claim): https://modelcontextprotocol.io/specification/2025-06-18/server/tools
- A2. Materials Project, *MP + LLMs (MCP)*, undated: https://docs.materialsproject.org/services/ml-and-ai-applications/mcp
- H1. reqwest, *Crate documentation*, 0.13.5/current page: https://docs.rs/crate/reqwest/latest
