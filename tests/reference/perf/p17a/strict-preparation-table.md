# Work versus delivered accuracy

Errors are against each code's own tight reference. Missing values are not zero.

| Case | Code / profile | Status | SCF | Davidson | H vectors | Wall (s) | dE (Ha) | max dEigen (Ha) | density L2 | sqrt(E_H) |
| --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| si8 | BAANDI/tight | ok | 20 | 594 | 2969 | 3.95168021 | -2.13162821e-14 | 2.57561937e-10 | 1.244743e-09 | 3.42005164e-09 |
| si8 | BAANDI/tighter | ok | 23 | 895 | 3730 | 4.91073554 | 0 | 0 | 0 | 0 |
| si8 | BAANDI/strict | ok | 15 | 345 | 2265 | 4.12324871 | 4.97379915e-13 | 4.08629737e-08 | 2.30261707e-07 | 6.76875764e-07 |
| k444s | BAANDI/tight | ok | 12 | 2150 | 12700 | 5.13885571 | -1.15463195e-14 | 9.81931758e-11 | 1.37393872e-10 | 2.59087783e-10 |
| k444s | BAANDI/tighter | ok | 14 | 2402 | 13954 | 4.63958804 | 0 | 0 | 0 | 0 |
| k444s | BAANDI/strict | ok | 10 | 1827 | 11181 | 2.68155133 | -1.77635684e-15 | 3.97484995e-09 | 5.31728585e-09 | 9.66318509e-09 |
| si-20 | BAANDI/tight | ok | 21 | 3856 | 20853 | 9.58906967 | 5.32907052e-15 | 4.8005902e-10 | 4.01771375e-10 | 6.40089326e-10 |
| si-20 | BAANDI/tighter | ok | 24 | 4290 | 22853 | 12.1944857 | 0 | 0 | 0 | 0 |
| si-20 | BAANDI/strict | ok | 15 | 2865 | 16526 | 6.50328283 | 3.37507799e-14 | 9.25526275e-09 | 1.91718776e-08 | 3.01537344e-08 |
| mgo-40 | BAANDI/tight | ok | 20 | 1371 | 12356 | 8.44968283 | -7.10542736e-14 | 2.24418262e-10 | 5.77213127e-10 | 6.61217263e-10 |
| mgo-40 | BAANDI/tighter | ok | 22 | 1504 | 13318 | 9.02840021 | 0 | 0 | 0 | 0 |
| mgo-40 | BAANDI/strict | ok | 18 | 1306 | 12202 | 8.27986083 | 3.83693077e-13 | 3.10354853e-09 | 8.57850237e-09 | 1.17943746e-08 |
| syn20-skew | BAANDI/tight | failed | unavailable | unavailable | unavailable | 0.030384708 | unavailable | unavailable | unavailable | unavailable |
| syn20-skew | BAANDI/tighter | failed | unavailable | unavailable | unavailable | 0.027204417 | unavailable | unavailable | unavailable | unavailable |
| syn20-skew | BAANDI/strict | ok | 6 | 15 | 255 | 0.022656833 | unavailable | unavailable | unavailable | unavailable |
| syn4-elongated | BAANDI/tight | ok | 11 | 41 | 215 | 0.009845 | 0 | 0 | 0 | 0 |
| syn4-elongated | BAANDI/tighter | ok | 11 | 42 | 220 | 0.010053791 | 0 | 0 | 0 | 0 |
| syn4-elongated | BAANDI/strict | ok | 9 | 31 | 171 | 0.00951975 | 5.32907052e-15 | 2.81200396e-10 | 1.37593162e-09 | 2.11611921e-09 |
| syn4-slab | BAANDI/tight | failed | unavailable | unavailable | unavailable | 0.094681125 | unavailable | unavailable | unavailable | unavailable |
| syn4-slab | BAANDI/tighter | failed | unavailable | unavailable | unavailable | 0.094506209 | unavailable | unavailable | unavailable | unavailable |
| syn4-slab | BAANDI/strict | failed | unavailable | unavailable | unavailable | 0.092868916 | unavailable | unavailable | unavailable | unavailable |
| syn16-gamma | BAANDI/tight | ok | 12 | 98 | 773 | 0.042154375 | 1.42108547e-14 | 3.97878619e-10 | 1.91598674e-09 | 6.62559953e-09 |
| syn16-gamma | BAANDI/tighter | ok | 14 | 106 | 893 | 0.044703542 | 0 | 0 | 0 | 0 |
| syn16-gamma | BAANDI/strict | ok | 10 | 69 | 593 | 0.035973875 | 7.10542736e-15 | 9.579191e-09 | 1.44077448e-08 | 4.37666869e-08 |
| syn16-k211 | BAANDI/tight | failed | unavailable | unavailable | unavailable | 0.18483675 | unavailable | unavailable | unavailable | unavailable |
| syn16-k211 | BAANDI/tighter | failed | unavailable | unavailable | unavailable | 0.180514542 | unavailable | unavailable | unavailable | unavailable |
| syn16-k211 | BAANDI/strict | failed | unavailable | unavailable | unavailable | 0.2001165 | unavailable | unavailable | unavailable | unavailable |
