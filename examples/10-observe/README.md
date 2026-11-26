# Observe the SCF event stream

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, from the [QE pseudopotential library](https://pseudopotentials.quantum-espresso.org/legacy_tables/original-qe-pp-library/si). Use a supported UPF v2 file; inspect it before running. The historical local validation file has unresolved provenance/license, so no redistribution or identity with a fresh download is implied.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 10-observe
input=target/tmp/examples/10-observe/input.toml
target/release/baandi validate "$input"
mkdir -p results/10-observe
target/release/baandi run "$input" --observe --events - --output-dir results/10-observe > results/10-observe/events-mirror.jsonl
```

`--events -` mirrors JSONL to stdout; the command above saves it separately from the engine's durable `events.jsonl` and `runtime.json`. Inspect map kinds (ordinary, verification, retry), safe points and policy revisions. A map is not always a new ordinary iteration. Observe is read-only. For supervised proposals use `--agent`, the local control socket and the [CLI modes guide](../../docs/user/cli-modes.md); controls need host approval and produce an audit. [MCP](../../crates/baandi-mcp/README.md) is a host-embedded adapter, not a standalone CLI server.

Laptop runtime budget: **seconds to a few minutes**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
