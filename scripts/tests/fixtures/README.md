# QE interchange fixtures

The `qe_legacy` inputs were captured with the pre-P16A generator from the
repository results listed in `qe_legacy/index.json`, using `/legacy/pseudo` and
`Si.pz-vbc.UPF`. These are byte-regression inputs, not new QE measurements.

The `materials` JSON files retain literal fields from the full P13D BAANDI
outputs below. Only unused top-level fields and convergence history/last-iteration
details were omitted. No numbers, commits, paths or hashes were reconstructed.
These are BAANDI output metadata, not QE references and contain no UPF content.

| Case | Original full result location | Full result SHA-256 |
| --- | --- | --- |
| si-20 | `/tmp/p13d-materials/si-20/results.json` | `65f5523b9bc00a53854a825d9cda321fdf8e5489248f8e7da202d30b7d179c41` |
| mgo-40 | `/tmp/p13d-materials/mgo-40/results.json` | `dc90f7a0636c6a93ed6321ffeca93800c04e33e4b9b255292a6cf58285a5f496` |

The executed-input hashes differ from the summary template hashes because the
P13D runner substituted absolute pseudopotential paths before invoking BAANDI.
The stored source commit is the commit reported by that executable, not a claim
that it identifies every uncommitted working-tree change in the original run.

Goldens in `tests/reference/materials/qe_inputs` are generated from these
projections. Their comment records the literal generation invocation. Tests
render to temporary paths while passing that same provenance string for stable
byte comparison; isolated CLI tests separately verify invocation recording.
