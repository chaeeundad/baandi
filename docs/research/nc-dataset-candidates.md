# Norm-conserving dataset candidates for the target materials (MF05 groundwork)

Date: 2026-11-05. Written by the orchestrator (Claude Code) after fetching and inspecting the files outside the repository (`/Users/mino/Develop/rust-dft-wt/pseudo/pseudodojo/`, never committed). Scope: facts needed to plan MF05 and MF06. Not a qualification: no cutoff scan, ghost check or transferability test has been run.

## Chosen first candidate set: PseudoDojo NC-SR v0.4, PBE, "standard" accuracy, UPF2

- Download: `http://www.pseudo-dojo.org/pseudos/nc-sr-04_pbe_standard_upf.tgz` (5,230,385 bytes), archive sha256 `603d7f07d364b60b321aa5c8eb0ae809aab339c3e884e4b62f31c35370991067`, 72 UPF files, fetched 2026-11-05.
- License: the PseudoDojo repository README states that the pseudopotential files and input files are released under CC BY 4.0 (credit the PseudoDojo project; citation list at http://www.pseudo-dojo.org/faq.html). The Python code is LGPL/GPL, not relevant for the data. The ONCVPSP generator header in each file asks (not requires) a citation of D. R. Hamann, Phys. Rev. B 88, 085117 (2013). The UPF files carry no license text inside, so the license reference must be recorded in the manifest, not read from the file. UPF files are still never committed to this repository (project rule); the archive hash and per-file hashes identify them.
- Generator: ONCVPSP 3.3.0 (08/16/2017) scalar-relativistic, PBE, UPF version 2.0.1.

| Element | file sha256 | z_valence | projectors | l_max | NLCC | note |
| --- | --- | --- | --- | --- | --- | --- |
| O | `2a543889f6522e2a1f94e2c2ff1477e61cb89f25c67b27c54940e8263ec6fb11` | 6 | 5 (l=0:2, l=1:2, l=2:1) | 2 | yes | 2 valence shells |
| Ni | `f83c9ee8062152d26e736c0ff9c0af7c86936f0642f3ab66ffb2eda42bd46a64` | 18 | 6 | 2 | yes | semicore 3s3p plus 3d4s |
| Al | `b02eaa07c5d98f5eeae1bec4155854a7a94ea3635ca6b42e914f3f7ebd6912fb` | 3 | 6 | 2 | yes | |
| Mg | `15cdd9380687a47a1f4bc4ab5fe1f0c6d9c50891350fd8534c8f292b3ce0be14` | 10 | 4 | 1 | no | semicore 2s2p |
| Si | `39822757f53f36e3bf3bfb779356152a8d3f21199c7db9dd5a931e5d18c45282` | 4 | 6 | 2 | yes | the repository's validation Si is a different, older file (Si.pz-vbc, LDA) |

All five: `pseudo_type=NC`, `relativistic=scalar`, `functional=PBE`. Header cutoffs are hints, not convergence evidence.

## Consequences for the plan (to be confirmed by tickets, not claims of support)

1. **NLCC is required** for O, Ni, Al and Si of this set (`core_correction=T`); the engine rejects it today (feature matrix E2). MF06 must implement a model core charge in the XC energy and potential before any of these files can run, with the translated-core derivative oracle in stage-oracles.md.
2. **Several projectors per angular momentum** (two for s and p in O) with a full coupling matrix: the engine supports multiple coupled projectors (synthetic `syn-s2` style) but must be checked on the real file layout (radial grid 936 to 1876 points, projector normalization, `PP_DIJ`).
3. **Semicore and many valence electrons**: Ni 18 and Mg 10 electrons per atom make cells expensive; the dev08 runs and the performance weakness of single-k-point larger cells matter for MF18.
4. **These are PBE files.** LDA-PZ results with them are not meaningful; spin-polarized PBE (MF03) is the relevant functional path, and QE comparisons must use the same files on the same grids.
5. Alternatives noted but not fetched: SG15 ONCV (CC BY-SA 4.0 per the Materials Cloud archive; share-alike obligations would bind any redistributed derived data) and pslibrary (not norm-conserving in general). The PseudoDojo choice avoids the share-alike condition.

## Next steps (owner: tickets after MF03)
- MF05: a manifest schema (hash, element, valence, XC, semicore, projector channel list, license reference, qualification id) and a read-only inspection tool in the CLI layer that reports exactly what the file needs and which engine capability is missing, with explicit blockers; independent check of the header values above against a second reader.
- MF06: NLCC and real-file coverage as above, cutoff scans and ghost/log-derivative checks against independently produced reference data (requested through the remote session).
