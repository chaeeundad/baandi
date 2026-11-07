# Band-count scan (orchestrator, Mac, 2026-11-07): work versus extra bands

Question raised by the dev08 QE sweep (tests/reference/perf/dev08_materials_qe_relayed.md, Q1 section; QE with nbnd equal to the occupied bands was 20 to 26 percent faster on si16 and si32, no gain on si8): how does BAANDI's work depend on the number of extra bands above the occupied ones? Method: `scripts/bench_scf.py` with `--extra-toml "nbands = N"`, one run per row, T=1, pure Rust release at main 2f1c, strict default profile, gamma_real off. Davidson steps and H vectors are deterministic counts (wall times were not used: the Mac was loaded by three concurrent builds, load average 8 to 12). Energies agree to about 1e-13 Ha across all rows of a case.

| Case (occupied bands) | extra bands: Davidson steps / H vectors |
| --- | --- |
| si8 (16) | 1: 720 / 2151; 2: 650 / 2579; 3: 413 / 2321; **4 (default): 283 / 2063**; 5: 212 / 1866; 6: 87 / 1253; 8: 67 / 1042; 10: 59 / 1092; 12: 63 / 1150 |
| Si primitive Gamma (4) | 1: 59 / 236; 2: 54 / 245; 3: 49 / 261; **4: 41 / 281**; 5: 44 / 308; 6: 42 / 334; 8: 41 / 396 |
| k222 (4) | 1: 463 / 1766; 2: 401 / 1964; **4: 408 / 2543**; 6: 412 / 2990; 8: 403 / 3457 |
| k444s (4) | 1: 1957 / 7574; 2: 1782 / 8349; **4: 1686 / 10290**; 6: 1727 / 12183; 8: 1800 / 14226 |

Reading: (1) for the three small cells the work grows with the number of bands, so one extra band needs 16 to 26 percent fewer H vectors than the default of four (k444s 7574 versus 10290); (2) si8 is the opposite: its top bands sit in a degenerate multiplet, a cut inside the multiplet (extra 1 to 5) makes the Davidson stagnate (the P7PC top-band regression, benchmark `si8-gamma`), and completing the multiplet (extra 6 to 10) cuts the work by 40 to 50 percent (1042 H vectors at extra 8 versus 2063 at the default 4); (3) the QE sweep shows the same two effects (nbnd = occupied is faster on si16 and si32, not on si8). So the band-count policy is a real lever with a degeneracy-dependent optimum; a single fixed default cannot be best for both groups. Not yet established: wall-time effect on dev08, other cells and other pseudopotentials, the standard profile, and whether completing a degenerate shell can be detected cheaply (for example from the first SCF iteration's eigenvalues).
