# ADR 0023: Metadata-only pseudopotential registry and explicit fetch boundary

Status: Proposed for review (P9B, 2026-11-03).

Users need reproducible file identity, legal attribution and actionable engine
compatibility before obtaining pseudopotentials. Family names alone cannot
establish any of these. The present reader has a three-projector total ceiling;
NC does not imply that every ONCV UPF can be loaded.

Store versioned JSON metadata and exact-file SHA-256 only. Preserve unknown facts
as null and unknown identities as unverified. Use the strict existing reader and
engine XC classifier to compute support, plus independent consistency checks.
Expose offline list/show/verify/path/suggest and cached fetch in the CLI adapter.
No registry I/O belongs in numerical code. No new dependency is required: use
existing serde_json, sha2 and clap, and Python's standard library for maintainers.

Downloads require explicit --yes, license verification and an exact verified
identity. Hash and size are checked before atomic, non-overwriting publication
outside repositories; provenance is recorded beside the file. Offline misses
fail. The production Fetcher remains unavailable until P9A decides HTTP policy;
redirects, TLS errors and hash mismatches must be fatal. Test the full acceptance
path with a synthetic in-memory transport, without shipping any real UPF bytes.

Input `registry:` resolution is a P8G contract, not new config behavior in P9B.
Source license uncertainties block automatic fetch, not inspection of user files.
SSSP original member licenses survive collection selection. PseudoDojo's official
data statement verifies CC BY 4.0 independently of its software license. SG15
license/bytes still require original-site rechecking in a network-enabled run.

Consequences: the initial registry may truthfully recommend zero files. Maintainers
must review sources and build exact-byte evidence before enabling entries. Parser
expansion and NLCC/SOC/PAW/USPP are separate roadmap work. No blanket family support
or presumed scientific accuracy is advertised. See the [design and source evidence](../architecture/pseudopotential-registry.md)
for schema semantics, cache failure handling, ranking and integration contracts.
