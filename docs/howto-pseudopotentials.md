# Finding and checking pseudopotentials

A pseudopotential describes the atomic core for one element. A file for the right
element can still use physics this program does not support. The registry records
where files come from, their licenses, citations and whether the engine can read
and use them.

This release provides the offline tools. The initial catalog contains seven
unverified candidates and **no approved download yet**. Actual internet downloads
will be enabled by P9A after the network client is selected. Nothing is downloaded
when you list, inspect or ask for suggestions.

```sh
dft-cli pseudo list --element Si
dft-cli pseudo list --supported
dft-cli pseudo show sg15-si-v1.2
dft-cli pseudo suggest SiO2
```

`show` includes the source website, license, citation text and any known cutoff.
A null value means it has not been checked. `suggest` explains why no supported
file is available; it does not choose an unsupported substitute. When verified
entries become available, suggestions list the lowest recommended wavefunction
cutoff first, followed by entries with no recommendation, then sort ties by ID.
Cutoffs are in Rydberg (Ry); 2 Ry = 1 Hartree. Check the XC label against your input
and test convergence: lower cost does not guarantee better accuracy.

If you already have a file outside the project directory:

```sh
dft-cli pseudo verify /path/to/Si.upf
dft-cli pseudo verify /path/to/Si.upf --id sg15-si-v1.2
dft-cli inspect-upf /path/to/Si.upf
```

A matching SHA-256 means the bytes match a registry record. `UNKNOWN` means the
file is not identified by the registry. A mismatch with `--id` means it is not the
expected file. Neither unknown nor mismatch is approval to use the file. Inspection
separately checks the format; calculation settings must also match the functional.
The machine-readable form is `inspect-upf /path/to/Si.upf --json`.

The current reader accepts only its norm-conserving UPF version 2 subset, with no
more than **three projectors in the entire file**. Nonlinear core correction,
spin-orbit, PAW and ultrasoft files are not supported. ONCV and PseudoDojo labels
alone do not guarantee compatibility. Missing support requires engine development,
not removing a flag from a downloaded file.

```sh
dft-cli pseudo path
dft-cli pseudo fetch VERIFIED_ID --offline
# After the P9A transport is installed and a verified entry is available:
dft-cli pseudo fetch VERIFIED_ID --yes
```

`VERIFIED_ID` is a placeholder, not an initial catalog entry. The offline command
uses only an existing cache and checks both file hash and provenance. The online
form currently reports the unavailable transport. Future downloads display the
license and citation and require your explicit `--yes` before connecting. A bad
hash is rejected. Existing different files are preserved. An unknown license
cannot be overridden with `--yes`; read the source website and obtain a verified
registry update. Cite original authors as well as the collection you used.

Files belong outside the repository. To choose another cache, set
`BAANDI_PSEUDO_DIR` to an absolute folder path outside all repositories. Defaults:

| Platform | Cache |
| --- | --- |
| macOS | `$HOME/Library/Caches/baandi/pseudopotentials` |
| Linux / other Unix | `$XDG_CACHE_HOME/baandi/pseudopotentials`, or `$HOME/.cache/baandi/pseudopotentials` |
| Windows | `%LOCALAPPDATA%/baandi/pseudopotentials` |

Do not commit UPF files. Keep their accompanying provenance JSON. An interrupted
cache write may leave an orphan file; a missing or incorrect manifest causes an
error instead of silently trusting it. Use a private cache directory.

Use ordinary pseudopotential file paths in current calculation inputs. Future
P8G integration will accept `Si = "registry:sg15-si-v1.2"` inside a
`[pseudopotentials]` table, but that syntax is **design only** today and will never
cause implicit downloading during a calculation.

The older validation file `Si.pz-vbc.UPF` has unverified provenance and license.
It is not an endorsed downloadable registry file. We recommend replacing it with
a verifiably licensed Si NC LDA file once it is registered and validated.
See [technical design and source links](architecture/pseudopotential-registry.md)
for details or give that page to the person maintaining your installation.
