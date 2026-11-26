# Strict oracle for validation

UPF: `Si.pz-vbc.UPF`, NC LDA-PZ, from the [QE pseudopotential library](https://pseudopotentials.quantum-espresso.org/legacy_tables/original-qe-pp-library/si). Use a supported UPF v2 file; inspect it before running. The historical local validation file has unresolved provenance/license, so no redistribution or identity with a fresh download is implied.

From the repository root, first [build and set PSEUDO_DIR](../README.md#setup).
No UPF or generated output belongs in this directory.

```sh
python3 examples/prepare.py 13-strict-oracle
input=target/tmp/examples/13-strict-oracle/input.toml
target/release/dft-cli validate "$input"
BAANDI_REFERENCE_STRICT=on target/release/dft-cli run "$input" --output-dir results/13-strict-oracle/strict
BAANDI_REFERENCE_STRICT=off target/release/dft-cli run "$input" --output-dir results/13-strict-oracle/production
```

`BAANDI_REFERENCE_STRICT=on` takes precedence and selects the frozen validation oracle. Inspect `scf.driver = "strict_oracle"` and the certificate before treating this as an oracle result. The production run is separate and should record `reference`. Compare energies in Ha only when both converge with identical physical settings; this tutorial supplies no measured difference or QE reference. See [ADR 0047](../../docs/decisions/0047-reference-driver-single-path.md).

Laptop runtime budget: **1–10 minutes for both solves**, release build, after compilation.
This is an unmeasured planning estimate, dependent on hardware, UPF and convergence.
**Validate-only in debug tests**; no SCF or optimizer is run by the example test.
