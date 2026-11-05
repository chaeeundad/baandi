# P19A symmetry evidence

No pseudopotential files or new QE outputs are stored here.

- `provenance.json`: executed commands, external UPF hashes, baseline commit and limits.
- `scf-oracles.log`: full/reduced Si2 (shifted/unshifted 4^3), Si8, MgO,
  cubic/skew/slab synthetic, spin and FD, with tightened caps. Reduced T=1/T=4
  energy and density bits are asserted by the shape/material helper.
- `check.log`: complete workspace build/test/fmt/clippy gate result.
- `legacy-baseline.json`, `legacy-t*.json`, `legacy-t*.log`, `legacy-gates.json`:
  five legacy cases, symmetry off, exact energy/density comparison to the
  pre-integration binary, T=1/2/4/8/12. Wall times in these regression files are
  not controlled performance measurements.
- `mac-k444s.json`: separate three-run alternating timing/work observation,
  including full input text and binary/UPF hashes. Compare at identical 32^3
  FFT discretization, not the legacy off-only 25^3 grid.

Reproduce from the repository root, with external datasets:

```sh
scripts/check.sh
DFT_SI_UPF=/external/Si.pz-vbc.UPF DFT_NC_DIR=/external/pseudodojo \
  cargo test --release -p dft-engine --test symmetry \
  -- --include-ignored --nocapture --test-threads=1
cargo test -p dft-geometry --test symmetry -- --nocapture
cargo test -p dft-physics --test symmetry -- --nocapture
cargo test -p dft-cli --test symmetry
python3 -I tests/reference/symmetry/measure.py \
  --binary target/release/dft-cli --upf /external/Si.pz-vbc.UPF \
  --output /tmp/symmetry-timings.json
python3 -I scripts/tests/test_readme_figures.py
```

See `docs/validation/symmetry.md` for the initial failed quadrature test, its
resolution, the unchanged NLCC floor policy, limitations and required owner edits.
The README generator does not yet consume this new result type; freshness of
its existing block is a separate gate from publishing these new measurements.
