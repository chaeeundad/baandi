# k-point reference bundles (R6)

Three Si runs (30 Ry, LDA-PZ) on Monkhorst-Pack grids: `k222_c` (2x2x2, unshifted), `k222s_c` (2x2x2, half-step shift), `k444s_c` (4x4x4, half-step shift). Each directory holds the QE input produced by `scripts/make_qe_input.py` (`qe.in`, `qe.in.meta.json` with the input and UPF sha256), our `ours_results.json` (rust-dft `results.json`), and the QE evidence: `qe.out` (QE 7.6, dev08) for `k222_c` and `k222s_c`; `qe_summary.json` for `k444s_c` (relayed plain text, see its `_note`).
The UPF file is not distributed (see `docs/candidate-si-upf.md`, `docs/mvp-review.md`).
Check: `python3 -I scripts/compare_kpoints_qe.py --ours <dir>/ours_results.json --qe-out <dir>/qe.out` (or `--qe-summary`), exit 0 means all default tolerances are met. Tests: `python3 -I -m unittest discover -s scripts/tests`. Results and caveats: `docs/validation.md` R6.
