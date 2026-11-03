# Pseudopotential registry (P9B)

Status: implemented offline adapter; HTTP transport intentionally unavailable pending
P9A. Input registry resolution is a P8G integration contract, not implemented here.
No numerical hot path changes, dependencies, real UPF files or reference results.

## Identity and evidence

`registry/pseudopotentials.json` is embedded in the CLI. `pseudo --registry PATH`
selects a maintainer-reviewed alternative; it is a trust input, not a download URL.
Run `python3 -I scripts/pseudo_registry_check.py PATH` before using an alternative.
The CLI also recomputes capability exclusions and refuses contradictory support.
The JSON Schema is `docs/schemas/pseudo-registry.schema.json` (2020-12).
The standard-library checker implements only the keywords used by this schema;
unknown keywords fail closed. It checks types, required fields, unknown fields,
patterns, numeric bounds, duplicate IDs, dates and consistency invariants.

An entry describes one exact file, never an archive. IDs are immutable, lowercase
ASCII identifiers; version changes or changed bytes require new IDs. SHA-256 is
lowercase hexadecimal over original bytes, size is bytes, date_checked is ISO date.
Null means unknown, never false or zero. Entries without SHA-256 must have
`status: unverified`. Verified means exact bytes were inspected and hashed, not
that a potential is accurate for every calculation. Deprecated entries keep their
hash and provenance for identification but cannot be recommended or fetched.
`license.verified` is independent of file identity and physics support. SPDX IDs
are preferred; a descriptive license identifier plus authoritative URL is allowed.

Source metadata includes family, version, exact file URL (null if unknown), and
project URL. Attribution and citations must name both collection and original
file authors. Formats list available/inspected formats; only an inspected `UPF2`
file can be supported. Functional is the actual header string, valence is electrons,
projector_counts maps angular momentum l to count. Relativistic, NLCC, SOC, PAW and
ultrasoft facts are explicit. Suggested wavefunction and density cutoffs are in Ry
(1 Ha = 2 Ry), with the exact recommendation source; null is no verified advice.
Source recommendations still require convergence testing for the user's system.

## Actual engine compatibility

The code was inspected, not inferred from the ticket: `parse_upf` accepts UPF
2.0/2.0.0/2.0.1, NC, scalar/nonrelativistic/`no`, and **at most three projectors
TOTAL per file**, not three per channel. l must be <=3. No NLCC, SOC, PAW/USPP,
GIPAW, Coulomb datasets, non-nominal units, or unsupported sections. Numerical
array and mesh validation must pass. Current engine XC accepts LDA-PZ and PBE
labels via `dft_engine::Functional::matches_upf`; PBE is already implemented in
this checkout. Support is not a statement that an LDA input may use a PBE file.

`inspect-upf FILE --json` returns facts and `supported_by_engine` with reasons.
A valid inspection result has exit 0 even when support is false. Malformed input
produces a false support result; metadata may be absent if it cannot be scanned.
I/O failures return nonzero. The ordinary text inspection keeps its earlier
behavior. An additive UPF header accessor permits reporting unsupported metadata;
it does not change strict parse acceptance or authorize numerical use.

| Family / inspected elements | Present result | Capability needed or remaining check |
| --- | --- | --- |
| SG15 (Si candidate; no exact file inspected) | Unverified, none recommended | Check total projector count, units/metadata and XC; NC alone is insufficient |
| PseudoDojo (Si candidate; no exact file inspected) | Unverified, none recommended | Choose UPF2 scalar variant; inspect total projectors, NLCC and XC; fully relativistic needs SOC |
| SSSP efficiency / precision (Si candidates) | Unverified, none recommended | Inspect original member; NC may qualify, USPP/PAW cannot; original license required |
| pslibrary USPP / PAW | Not usable | Augmentation density, generalized overlap, PAW/USPP physics |
| GBRV UPF ultrasoft | Not usable | Ultrasoft augmentation and generalized overlap |
| Hamann ONCVPSP outputs | Unverified, none recommended | Generator is not a licensed file collection; inspect each output and its license |

No transition-metal or other element is claimed compatible without file evidence.
If an inspected transition-metal file has NLCC, `suggest` reports nonlinear core
correction as the missing capability; it never substitutes an unsupported entry.
Reader expansion belongs to feature-gap survey T4; NLCC, SOC and PAW/USPP physics
belong to P8D. See [feature-gap survey](../research/feature-gap-survey.md).
P8D's planned `docs/research/xc-and-manybody-roadmap.md` was absent in this checkout;
link/coordinate its capability milestones when it lands, without presuming dates.

## Offline commands and ranking

`list [--element Si] [--supported]` prints status and reasons. `show ID` prints
complete JSON including license, citations and links. `verify FILE` hashes bytes,
reports all matching IDs and their support/status; an unknown hash fails nonzero.
`verify FILE --id ID` distinguishes an expected-identity mismatch from an unknown
file. A match is identity evidence, not a new engine validation. `path` performs
no writes. `suggest SiO2` or `suggest Si O` accepts positive integer stoichiometry,
space/comma separated elements and returns unique elements in alphabetical order.
Parentheses, charges, isotopes and fractional stoichiometry are not supported.

Suggestions include only verified, nondeprecated, supported entries, sorted by
ascending recommended wavefunction cutoff, unknown cutoff last, then ID. This is
a reproducible cost-oriented listing, not a scientific quality ranking or an
automatic XC choice. Each result prints XC, cutoff units/source and ranking rule.
No candidates produces an explicit compatibility/coverage message. There is no
fallback to unsafe files. License clearance is additionally required by fetch.

## Cache and fetch boundary

Cache resolution: absolute `BAANDI_PSEUDO_DIR` overrides everything. Otherwise:
macOS `$HOME/Library/Caches/baandi/pseudopotentials`; Linux/other Unix
`$XDG_CACHE_HOME/baandi/pseudopotentials`, falling back to
`$HOME/.cache/baandi/pseudopotentials`; Windows
`%LOCALAPPDATA%/baandi/pseudopotentials`. Empty/relative configured paths fail.
Fetch rejects repository paths (including existing symlink parents), `..`, and
nonregular cached files. Use a private user-owned cache, not a shared writable
directory: this std-only implementation is not hardened against hostile concurrent
replacement of directory symlinks. Platform defaults are compile-time selected;
macOS was exercised here, Linux/Windows branches need native CI.

`pseudo fetch ID [--yes] [--offline]` prints license, attribution, citations and
project URL first. Unknown licenses require reading the source page and maintainer
verification; `--yes` does not override that block. Unverified, unsupported and
deprecated entries are refused. A cached hit is rehashed and its manifest checked
before success. Existing different bytes are never overwritten or deleted.
Offline mode never calls transport and a miss fails. A miss without `--yes` fails
before any network access; this release chooses explicit confirmation only, not
an interactive prompt (safe for automation). There is no cache creation on misses
with the unavailable production transport.

`Fetcher::download(url, max_bytes)` in `crates/dft-cli/src/pseudo.rs` is the marked
P9A seam; `NetworkUnavailable` is the only production implementation. The future
client must enforce HTTPS, normal TLS certificate/hostname checks, no redirects,
status 200 only, timeout and size limit. TLS errors and redirects must fail, never
retry with weaker security. The controller verifies SHA-256 and exact size before
publication. Current test transport returns bytes in memory: a mismatch drops
them and creates no temporary file. A streaming P9A replacement must remove partial
downloads on any failure and preserve this verification-before-publication rule.

After validation, a create-new staging file is synced and hard-linked atomically
to `<id>-<sha256>.upf`, which cannot replace another entry. Staging is removed on
success/error. `<id>-<sha256>.json` is created exclusively with id, URL, SHA-256
and `downloaded_at_unix_seconds` (UTC Unix seconds, the download date/time).
Manifest creation failure fails the operation; a verified payload can remain
orphaned. Subsequent hits fail closed until the user repairs provenance or removes
the orphan. A crash can leave staging debris: do not interpret it as an accepted
file. Hard-link publication requires filesystem support. Cache reads require the
manifest and recheck identity even with `--offline`.

## Input and preflight integration (design only)

```toml
[pseudopotentials]
Si = "registry:sg15-si-v1.2"
O = "/absolute/path/to/O.upf"
```

P8G resolves `registry:` IDs using this registry/cache without implicit network.
A miss gives the explicit fetch command. It checks file hash/size, element against
the species key, strict parser support, calculation XC, and provenance. Plain
paths preserve existing semantics, but unknown registry identity should be
reported honestly. Numerical code receives parsed potentials only. Current input
parser does not implement the example's registry syntax; this is an integration
contract and must not be advertised as working SCF input yet.

## Maintainer workflow

Start with the registry entry fields as a complete metadata template; map each
local filename (or exact URL) to its template in a JSON object. Fill immutable ID,
original direct file URL, version, legal evidence, author attribution and citations.
For SSSP preserve the original member's license/citation and the collection's
version/variant. Supply recommendation metadata only after reading the source.
Run outside any repository for local downloaded files:

```sh
python3 -I scripts/pseudo_registry_build.py /outside/downloads \
  --metadata /outside/metadata.json --cli target/debug/dft-cli \
  --output /outside/candidate-registry.json
python3 -I scripts/pseudo_registry_check.py /outside/candidate-registry.json
```

The builder uses the standard library, computes SHA-256 of exact bytes, invokes
the strict Rust JSON inspection and derives support conservatively. Unsupported
files retain header facts and parser reasons. Local inspection can record an
unverified license; it cannot authorize fetch. HTTPS URL inputs require explicit
`--allow-network`, verified license metadata, exact matching URL, no redirects,
default TLS validation, 30-second timeout and 64 MiB cap. No archive extraction.
Temporary inspection bytes are deleted outside the repository. Output is create-new
and never replaces an existing registry. Review the candidate's source claims
before copying metadata into the repository: computed hash alone cannot establish
that a manually supplied file came from its claimed URL.

## Source and legal evidence (checked 2026-11-03)

Only metadata/hashes are stored here. We neither redistribute nor modify UPF files;
users obtain original bytes from the original site. CC BY requires appropriate
credit, license link and change indication; CC BY-SA additionally governs sharing
adapted material. Keep notices/citations when using or sharing files. These terms
do not turn our metadata catalog into a blanket relicense of the files. See the
[CC BY deed](https://creativecommons.org/licenses/by/4.0/) and
[CC BY-SA deed](https://creativecommons.org/licenses/by-sa/4.0/).

- [SG15 project](https://www.quantum-simulation.org/potentials/sg15_oncv/):
  Schlipf and Gygi, CPC 196, 36 (2015), doi:10.1016/j.cpc.2015.05.011.
  CC BY-SA 4.0 and v1.0/v1.2 were reported in the task; the original site timed
  out on recheck. The 2020-02-06 archive date was not independently verified.
  Initial license verification remains false, SHA-256/size/file URL null.
- [PseudoDojo official repository license statement](https://github.com/abinit/pseudo_dojo#license)
  explicitly assigns CC BY 4.0 to potentials/input data, distinct from its code
  license; UPF2, psp8 and PSML are listed. van Setten et al., CPC 226, 39 (2018),
  doi:10.1016/j.cpc.2018.01.012. Data license uncertainty is resolved at family
  level; file identity and compatibility remain unverified.
- [SSSP Discover](https://legacy.materialscloud.org/discover/sssp) describes
  efficiency/precision variants and wavefunction/density recommendations in Ry.
  It explicitly preserves each original member's license and attribution. Do
  not apply a collection-level CC BY claim to every member. Prandini et al.,
  npj Comput. Mater. 4, 72 (2018), doi:10.1038/s41524-018-0127-2.
  [Materials Cloud record 2023.65](https://archive.materialscloud.org/record/2023.65)
  lists v1.3.0 PBE/PBEsol cutoff JSON and verification data. JSON/license downloads
  redirected to inaccessible object storage, so no cutoff values or file hashes
  were copied. Published MD5 checksums are not SHA-256 and cannot be substituted.
- [pslibrary](https://dalcorso.github.io/pslibrary/) describes USPP/PAW and GPL;
  Dal Corso, Comput. Mater. Sci. 95, 337 (2014). Exact generated file licensing
  needs checking; no specific GPL version is invented in the registry.
- [GBRV](https://www.physics.rutgers.edu/gbrv/) describes ultrasoft UPF and PAW
  variants and GNU public licensing; Garrity et al., Comput. Mater. Sci. 81, 446
  (2014). Exact per-file version/license remains unverified; UPF USPP is unusable.
- [Hamann ONCVPSP](http://www.mat-simresearch.com/) timed out on recheck.
  Generator outputs have no blanket verified license or engine compatibility here.

Shell retrieval of the SG15 page failed with DNS resolution error in the sandbox.
Web page retrieval is not a byte-preserving UPF download. Consequently all seven
initial candidates have null SHA-256 and unverified status; no real UPF was
inspected and the matrix is explicitly unverified where appropriate.

The validation `Si.pz-vbc.UPF` has unresolved provenance and license, as documented
in [candidate Si facts](../candidate-si-upf.md). Its historical digest is not
recast as newly verified evidence. Replace it with a registered, verifiably
licensed Si NC LDA-PZ file once one passes the reader and reference validation;
do not silently substitute SG15 PBE for an LDA reference calculation.

## Validation evidence for this implementation

On 2026-11-03, `scripts/check.sh` exited 0 with PASS build, test, fmt and clippy.
`python3 -I -m unittest discover -s scripts/tests` exited 0: 119 tests in 48.690 s
(existing ResourceWarnings in unrelated tests). The new Python module has 6 tests;
the CLI integration suite has 5 tests (0.62 s), and the injected fetch-boundary
unit test has 1 test (0.01 s). The registry checker exited 0 for 7 entries, all
unverified with null hashes. Tests use generated synthetic bytes, not real UPFs.
No new physical accuracy or QE comparison result is claimed.

Shared-file changes are limited to `main.rs`: the `mod pseudo` declaration,
`InspectUpf` JSON flag, `Pseudo` variant, two early dispatch paths and two exhaustive
match arms (22 added lines, 2 removed after formatting). No changes to config.rs,
run.rs, manifests, README, tickets, baselines or artifacts. The UPF change adds a
header-facts accessor using the existing scanner with a separate inspection mode;
the strict parser still invokes its original rejection policy. Changes are left
uncommitted for review as required by the task.
