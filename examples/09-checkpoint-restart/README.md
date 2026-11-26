# Map-boundary checkpoint and restart

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, from the [QE pseudopotential library](https://pseudopotentials.quantum-espresso.org/legacy_tables/original-qe-pp-library/si). Use a supported UPF v2 file; inspect it before running. The historical local validation file has unresolved provenance/license, so no redistribution or identity with a fresh download is implied.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 09-checkpoint-restart
input=target/tmp/examples/09-checkpoint-restart/input.toml
target/release/dft-cli validate "$input"
mkdir -p results/09-checkpoint-restart
target/release/dft-cli checkpoint "$input" results/09-checkpoint-restart/boundary.cp --iteration 3 --kind exact
target/release/dft-cli restart "$input" results/09-checkpoint-restart/boundary.cp results/09-checkpoint-restart/result.cp
target/release/dft-cli status results/09-checkpoint-restart/result.cp
```

Capture at a nonterminal safe map boundary once iteration 3 is reached; map count and iteration count can differ. `cutoff_ladder = "off"` gives supported fixed-cutoff continuation. Exact restart requires identical input bytes, UPFs, build/backend/platform, threads, ranks and partition. `result.cp` is a **ResultExport container**, not a results directory or another restart boundary. `result.cp.results.json` carries routing metadata only; use `status` for the export summary. No warm v2 import or v1 conversion. See [checkpoint format](../../docs/checkpoint-format.md).

Laptop runtime budget: **seconds to a few minutes**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
