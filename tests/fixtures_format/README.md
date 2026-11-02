# FORMAT FIXTURES ONLY

These files are invented schema examples, not QE output or physical results.
`results.json` and `bundle/provenance.json` carry `is_format_fixture: true`.
The verifier must return 4, `FORMAT FIXTURE: not a reference`, for this bundle.
`golden.in` was written from the documented input contract for the synthetic
non-symmetric lattice. `bundle/qe.out` was constructed from invented values;
no pw.x command was executed. No pseudopotential is stored here.
