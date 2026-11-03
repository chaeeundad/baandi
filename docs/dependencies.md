# Dependency inventory

Versions are resolved package versions from [Cargo.lock](../Cargo.lock), checked
against workspace manifests and `cargo metadata --offline --format-version 1`.
License expressions below are metadata declarations, not a legal determination.
The toolchain is Rust 1.99.1, edition 2024, workspace resolver 3.

## Direct production dependencies

| Package | Locked version | Declared license | Consumers |
| --- | --- | --- | --- |
| clap | 4.6.7 | MIT OR Apache-2.0 | dft-cli |
| exmex | 0.21.0 | MIT OR Apache-2.0 | dft-units (all features off) |
| faer | 0.24.4 | MIT | dft-numerics |
| libm | 0.2.16 | MIT | dft-physics, dft-upf |
| nalgebra | 0.35.0 | Apache-2.0 | dft-basis, dft-geometry, dft-physics |
| ndarray | 0.17.2 | MIT OR Apache-2.0 | dft-basis, dft-numerics, dft-physics |
| ndrustfft | 0.6.2 | MIT | dft-numerics |
| num-complex | 0.4.6 | MIT OR Apache-2.0 | dft-basis, dft-numerics, dft-physics, dft-solvers |
| rustfft | 6.4.1 | MIT OR Apache-2.0 | dft-numerics |
| serde | 1.0.229 | MIT OR Apache-2.0 | dft-cli |
| serde_json | 1.0.154 | MIT OR Apache-2.0 | dft-cli |
| sha2 | 0.10.9 | MIT OR Apache-2.0 | dft-cli |
| thiserror | 2.0.21 | MIT OR Apache-2.0 | dft-basis, dft-cli, dft-geometry, dft-numerics, dft-physics, dft-solvers, dft-upf |
| toml | 1.1.8+spec-1.1.0 | MIT OR Apache-2.0 | dft-cli |
| uom | 0.38.0 | Apache-2.0 OR MIT | dft-units (runtime and build) |
| upf | 0.1.1 | MIT OR Apache-2.0 | dft-upf |
| xcx | 0.3.1 | (MIT OR Apache-2.0) AND MPL-2.0 | dft-physics |

`xcx` 0.3.1 is a user-approved normal dependency of `dft-physics` for PBE.
It also cross-checks LDA exchange and provides a loose PW92 correlation sanity
check in tests, not a PZ81 oracle. Its MPL-2.0 component is file-level copyleft
derived from libxc and is now linked into released binaries.
[ADR 0018](decisions/0018-pbe-via-xcx.md) supersedes the dev-only restriction
of ADR 0007. The production LDA-PZ path remains in-tree in `dft-physics`.
There are no direct external dev-only packages.

dft-units also uses uom as a build-dependency to preserve its public constant API.
The CLI build script uses std only.
Workspace dev edges additionally include `dft-basis -> dft-numerics` and
`dft-cli -> dft-physics`; all workspace edges are in
[architecture](architecture.md). `dft-units` uses uom with default features disabled and `f64`, `si`, `std` enabled.

## Transitive inventory and licenses

The unfiltered offline metadata graph contains **343 external package/version
entries** plus nine workspace crates: 17 direct production packages, no direct
external dev-only package, and **326 additional transitive entries**. This includes
development and target-specific resolution, not just libraries linked into a
host release binary. Multiple versions of one crate count separately. In
particular, direct nalgebra is 0.35.0 and direct thiserror is 2.0.21; older
versions in the graph are transitive.

Declared license families are MIT, Apache-2.0, BSD-2-Clause, Zlib, Unlicense,
Unicode-3.0, MPL-2.0 and LGPL-2.1-or-later, with LLVM-exception on some Apache
alternatives. Metadata also contains legacy slash-separated MIT/Apache and
Unlicense/MIT spellings. Unicode appears both alone and in AND expressions.
There are no missing license fields in this graph. The complete list of
copyleft-containing expressions is:

| Package | Expression | Introduced through |
| --- | --- | --- |
| xcx 0.3.1 | `(MIT OR Apache-2.0) AND MPL-2.0` | dft-physics normal dependency |
| r-efi 5.3.0 | `MIT OR Apache-2.0 OR LGPL-2.1-or-later` | dft-numerics -> faer 0.24.4 -> rand 0.9.5 -> rand_core 0.9.5 -> getrandom 0.3.4 -> r-efi |

The r-efi expression offers permissive alternatives; it is not LGPL-only.
The xcx expression includes MPL with AND, rather than offering it as an OR
choice. No GPL-only or LGPL-only expression appears in this metadata snapshot.

## Adoption and offline sources

[ADR 0003](decisions/0003-prefer-general-purpose-libraries.md) requires license,
maintenance and Apple Silicon checks before adopting a numerical library.
[ADR 0004](decisions/0004-backend-libraries.md) selects faer and ndrustfft/rustfft;
[ADR 0005](decisions/0005-small-utility-crates.md) covers geometry and utilities,
[ADR 0007](decisions/0007-xc-and-special-functions.md) XC/special functions,
[ADR 0008](decisions/0008-upf-and-radial-math.md) UPF, and
[ADR 0010](decisions/0010-cli-output-crates.md) checksums and JSON.
The exact metadata license for libm 0.2.16 is MIT; ADR 0007's broader
MIT OR Apache-2.0 description is historical and differs from this package.
The former native adapter was removed; current workspace manifests declare no
LAPACK or Accelerate source-selector crates. The
[backend spike](backend-spike-results.md) is historical evidence, not a current
production dependency list or a general performance guarantee.

Per ADR 0004, Claude fetched and vendored the crates for offline Codex use.
The shared vendored directory is untracked; `vendor/` is ignored in
[.gitignore](../.gitignore). This worktree's local `.cargo/config.toml` replaces
crates.io with `../rust-dft-vendor`, and
[scripts/delegate.sh](../scripts/delegate.sh) configures that shared source for
new worktrees when available. Existing records establish this workflow, but do
not give the original exact vendoring command; none is reconstructed here.
No dependencies or vendor contents were changed for this documentation refresh.

The Python tooling under `scripts/` uses only the standard library and sibling
modules; comparison tools support Python 3.9+ and are invoked with `python3 -I`.
No pip dependency is required. Shell verification uses Cargo; the macOS benchmark
helper uses `/usr/bin/time -l` as documented in [CLI documentation](cli.md).

The project's own license is MIT (see `LICENSE`; selected by the user on 2026-11-03).

P7A-KPOINTS (ADR 0011): no new dependencies or manifest changes. Bloch sampling
and weighted SCF use the existing geometry, FFT, complex-array and Davidson APIs.


P7U adds uom 0.38.0, already fetched and license-reviewed by Claude. Its enabled
features are `f64`, `si`, `std` (no default features); no serialization or other
storage types are enabled. Its transitive crates are num-traits 0.2.19 and
typenum 1.20.1 (both MIT OR Apache-2.0), plus num-traits' build dependency
autocfg 1.5.1 (MIT OR Apache-2.0). These transitive packages were already present.
The pure-Rust implementation needs no system library or FFI; offline host gates
exercise Apple Silicon support. See ADR 0012 for adoption and precision checks.


P7E adds exmex 0.21.0 with default features disabled and no optional features;
`cargo metadata --offline --format-version 1` confirms its resolved feature list
is empty. The regex dependencies retain their own features. License declarations
and versions below come from that metadata and agree with Cargo.lock:

| Transitive package | Version | License | Present before exmex adoption |
| --- | --- | --- | --- |
| regex | 1.13.1 | MIT OR Apache-2.0 | yes |
| regex-automata | 0.4.18 | MIT OR Apache-2.0 | yes |
| regex-syntax | 0.8.11 | MIT OR Apache-2.0 | yes |
| aho-corasick | 1.1.5 | Unlicense OR MIT | yes |
| memchr | 2.8.3 | Unlicense OR MIT | yes |
| smallvec | 1.16.4 | MIT OR Apache-2.0 | yes |

The exmex direct edges are regex and smallvec; the remaining packages are regex
transitives. The older regex-syntax 0.6.29 elsewhere in the graph is not an
exmex transitive. All six versions above were already in the parent lockfile
before adoption commit 7ba4c04; only exmex itself was added by that commit.
The supplied pure-Rust source builds offline on Apple Silicon without FFI.
ADR 0014 records API restrictions, maintenance-review limits and validation.

P7P1-FFT-ALLOC uses rustfft directly for production 3D axis passes with persistent
plans, scratch and a strided lane buffer. ndrustfft remains in the unchanged
manifest and is used only by the historical FFT test oracle; moving it to a
dev-dependency is a separate manifest decision for Claude. No dependency versions
or features changed.

### P7B linked PBE component

Per accepted ADR 0018, xcx 0.3.1 is a normal dft-physics dependency, used in production
behind the spin-unpolarized PBE grid adapter. License: `(MIT OR Apache-2.0)
AND MPL-2.0`; linked releases must retain the MPL notices and source availability
for the covered library files. No xcx source is modified. Pure Rust, edition
2021, MSRV 1.87; no native/FFI dependency added by this ticket. The independent
closed-form oracle and discrete variational tests constrain replacement of this
young, single-author dependency. Existing LDA-PZ code is unchanged.
