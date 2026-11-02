# si30_ours: first real QE reference bundle

- `qe.in`, `qe.out`, `provenance.json`: returned by the remote runner (dft5090 session, dev08), verified on receipt by SHA-256 (input, gzip, stdout). QE 7.6, MPI x4, `nr1=nr2=nr3=25` as requested.
- `ours/results.json`: the rust-dft run this bundle is compared against (Si diamond, a = 5.43 Angstrom, Gamma only, 30 Ry; see docs/validation.md R1). It records absolute paths and a timestamp of the machine that ran it.
- Pseudopotential Si.pz-vbc.UPF: provenance unverified (see docs/validation.md). The UPF itself is not in the repository.
- Run B (`qe_si30_auto`, QE's automatic grid) was reported by the runner as giving the same grid 25x25x25 and the same energy; its stdout is not stored here until its transfer is verified.
