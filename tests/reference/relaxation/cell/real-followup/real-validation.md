# P23A actual material runs

Pressure units: Ha/Bohr^3; lattice: Bohr; time: seconds.

| Case | Lattice | Frozen P | Rebuilt P | Rebuilt - frozen P | Time | Exit |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| Si-30-0.98 | 10.334460599298 | 1.830861273e-11 | 3.715320317e-08 | 3.713489456e-08 | 115.548 | 0 |
| Si-30-1.02 | 10.334514952422 | 6.341755246e-12 | -1.028274771e-08 | -1.028908946e-08 | 113.395 | 0 |
| Si-40-0.98 | 10.334354861162 | 1.836427307e-11 | -1.129938760e-08 | -1.131775187e-08 | 171.004 | 0 |
| Si-40-1.02 | 10.334352638539 | 6.435082791e-12 | -9.359662457e-09 | -9.366097540e-09 | 180.640 | 0 |
| MgO-50-0.98 | 8.036842438247 | 1.236542350e-11 | 1.758221561e-07 | 1.758097907e-07 | 73.002 | 0 |
| MgO-50-1.02 | 8.036996494869 | 3.875580412e-12 | -1.157683548e-07 | -1.157722304e-07 | 71.241 | 0 |
| MgO-60-0.98 | 8.037045644138 | 9.664112103e-12 | -2.223937761e-07 | -2.224034402e-07 | 95.135 | 0 |
| MgO-60-1.02 | 8.036914908800 | 2.671853479e-12 | 2.479184669e-08 | 2.478917484e-08 | 91.944 | 0 |

## Start differences

Optimizer lattice tolerances below are local slope estimates from the original stress tolerance; they are not an independently specified lattice tolerance.

- Si 30 Ha: expanded minus compressed = 5.435312430e-05 Bohr; combined estimated tolerance = 2.310411789e-04 Bohr; within estimate = True.
- Si 40 Ha: expanded minus compressed = -2.222622646e-06 Bohr; combined estimated tolerance = 2.310502669e-04 Bohr; within estimate = True.
- MgO 50 Ha: expanded minus compressed = 1.540566228e-04 Bohr; combined estimated tolerance = 1.069353092e-04 Bohr; within estimate = False.
- MgO 60 Ha: expanded minus compressed = -1.307353382e-04 Bohr; combined estimated tolerance = 1.069962277e-04 Bohr; within estimate = False.

## Cutoff dependence

- Si start factor 0.98: higher minus lower cutoff lattice = -1.057381356e-04 Bohr.
- Si start factor 1.02: higher minus lower cutoff lattice = -1.623138826e-04 Bohr.
- MgO start factor 0.98: higher minus lower cutoff lattice = 2.032058914e-04 Bohr.
- MgO start factor 1.02: higher minus lower cutoff lattice = -8.158606961e-05 Bohr.

Finite frozen-Miller minima and rebuilt nominal-cutoff diagnostics; no experiment or external-code comparison.
