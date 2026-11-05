# Work versus delivered accuracy

Errors are against each code's own tight reference. Missing values are not zero.

| Case | Code / profile | Status | SCF | Davidson | H vectors | Wall (s) | dE (Ha) | max dEigen (Ha) | density L2 | sqrt(E_H) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| syn4-elongated | BAANDI/tight | ok | 11 | 41 | 215 | 0.0107815 | 0 | 0 | 0 | 0 |
| syn4-elongated | BAANDI/tighter | ok | 11 | 42 | 220 | 0.009984458 | 0 | 0 | 0 | 0 |
| syn4-elongated | BAANDI/strict | ok | 9 | 31 | 171 | 0.009599875 | 5.32907052e-15 | 2.81200396e-10 | 1.37593162e-09 | 2.11611921e-09 |
| syn4-elongated | BAANDI/standard-5e-08 | failed | unavailable | unavailable | unavailable | 0.002401333 | unavailable | unavailable | unavailable | unavailable |
| syn4-elongated | BAANDI/standard-5e-09 | failed | unavailable | unavailable | unavailable | 0.002542875 | unavailable | unavailable | unavailable | unavailable |
| syn4-elongated | BAANDI/standard-5e-10 | failed | unavailable | unavailable | unavailable | 0.002805166 | unavailable | unavailable | unavailable | unavailable |
| syn4-elongated | BAANDI/strict-hartree | failed | unavailable | unavailable | unavailable | 0.002429292 | unavailable | unavailable | unavailable | unavailable |
