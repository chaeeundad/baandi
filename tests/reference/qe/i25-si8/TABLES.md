## si8

| k | mode | stop line | step-k estimated scf accuracy (Ry) | step-k ethr | step-k avg # Davidson its | step-k total energy (Ry) | '!' energy | wall s |
|---|---|---|---|---|---|---|---|---|
| 1 | mix | convergence NOT achieved after   1 iterations: stopping | 0.43855533 | 1.00E-02 | 3.0 | -62.61942841 |  | 0.16 |
| 1 | out | convergence has been achieved in   1 iterations | 0.43855533 | 1.00E-02 | 3.0 | -62.60391409 | -62.60391409 | 0.16 |
| 2 | mix | convergence NOT achieved after   2 iterations: stopping | 0.01956404 | 1.37E-03 | 1.0 | -62.64131413 |  | 0.20 |
| 2 | out | convergence has been achieved in   2 iterations | 0.01956404 | 1.37E-03 | 1.0 | -62.64125463 | -62.64125463 | 0.16 |
| 3 | mix | convergence NOT achieved after   3 iterations: stopping | 0.00085399 | 6.11E-05 | 4.0 | -62.64319427 |  | 0.18 |
| 3 | out | convergence has been achieved in   3 iterations | 0.00085399 | 6.11E-05 | 4.0 | -62.64316560 | -62.64316560 | 0.17 |
| 4 | mix | convergence NOT achieved after   4 iterations: stopping | 0.00120905 | 2.67E-06 | 3.0 | -62.64337763 |  | 0.23 |
| 4 | out | convergence has been achieved in   4 iterations | 0.00120905 | 2.67E-06 | 3.0 | -62.64290026 | -62.64290026 | 0.17 |
| 5 | mix | convergence NOT achieved after   5 iterations: stopping | 0.00002339 | 2.67E-06 | 2.0 | -62.64356086 |  | 0.22 |
| 5 | out | convergence has been achieved in   5 iterations | 0.00002339 | 2.67E-06 | 2.0 | -62.64355385 | -62.64355385 | 0.22 |
| 6 | mix | convergence NOT achieved after   6 iterations: stopping | 0.00000495 | 7.31E-08 | 6.0 | -62.64357264 |  | 0.19 |
| 6 | out | convergence has been achieved in   6 iterations | 0.00000495 | 7.31E-08 | 6.0 | -62.64357081 | -62.64357081 | 0.21 |

Trajectory prefix identical across all runs (acc, energy, ethr per iteration): yes

Full trajectory (from the k=6 mix run):

| iteration | estimated scf accuracy (Ry) | ethr | avg # Davidson its | total energy (Ry) |
|---|---|---|---|---|
| 1 | 0.43855533 | 1.00E-02 | 3.0 | -62.61942841 |
| 2 | 0.01956404 | 1.37E-03 | 1.0 | -62.64131413 |
| 3 | 0.00085399 | 6.11E-05 | 4.0 | -62.64319427 |
| 4 | 0.00120905 | 2.67E-06 | 3.0 | -62.64337763 |
| 5 | 0.00002339 | 2.67E-06 | 2.0 | -62.64356086 |
| 6 | 0.00000495 | 7.31E-08 | 6.0 | -62.64357264 |

Density files: gamma_only=1 ngm_g=11936 nspin=1 G0=0.0296178744919
Miller index list identical in every file: yes

### Relative differences ||a-b|| / ||b|| (2-norm over the stored half-sphere coefficients)

| check | value | meaning |
|---|---|---|
| k1 mix beta=1.0 vs k1 out | 3.006e-18 | first-step mix at beta=1 equals rho_out(1) up to rounding |
| rho_in(1) from beta 0.7 vs from beta 0.3 | 1.723e-16 | reconstruction rho_in(1) = (rho_in(2) - beta rho_out(1))/(1-beta) is consistent |
| ||rho_out(1) - rho_in(1)|| / ||rho_in(1)|| | 1.968e-01 | size of the first update |
| rho_in(1) G0 | 0.0296178744919 | electron count / Omega |
| k1: mix(k) vs out(k) | 5.578e-02 | rho_in(k+1) vs rho_out(k) |
| k2: mix(k) vs out(k) | 7.079e-03 | rho_in(k+1) vs rho_out(k) |
| k3: mix(k) vs out(k) | 1.799e-03 | rho_in(k+1) vs rho_out(k) |
| k4: mix(k) vs out(k) | 4.035e-03 | rho_in(k+1) vs rho_out(k) |
| k5: mix(k) vs out(k) | 5.460e-04 | rho_in(k+1) vs rho_out(k) |
| k6: mix(k) vs out(k) | 2.856e-04 | rho_in(k+1) vs rho_out(k) |
| k2: mix(k) vs linear rho_in(k)+0.7(rho_out(k)-rho_in(k)) | 7.042e-03 | 0 would mean plain linear mixing; nonzero = Broyden history |
| k3: mix(k) vs linear rho_in(k)+0.7(rho_out(k)-rho_in(k)) | 3.565e-04 | 0 would mean plain linear mixing; nonzero = Broyden history |
| k4: mix(k) vs linear rho_in(k)+0.7(rho_out(k)-rho_in(k)) | 2.314e-03 | 0 would mean plain linear mixing; nonzero = Broyden history |
| k5: mix(k) vs linear rho_in(k)+0.7(rho_out(k)-rho_in(k)) | 2.823e-04 | 0 would mean plain linear mixing; nonzero = Broyden history |
| k6: mix(k) vs linear rho_in(k)+0.7(rho_out(k)-rho_in(k)) | 1.538e-04 | 0 would mean plain linear mixing; nonzero = Broyden history |

### Recomputed 'estimated scf accuracy' from the stored files (decides what each file holds)

dr2(k) recomputed with QE's own metric from rho_in(k) and the k-out file; rho_in(1) = reconstruction above, rho_in(k>=2) = the (k-1)-mix file.
Agreement to every printed digit means: k-out holds rho_out(k) and (k-1)-mix holds rho_in(k), for every k.

| k | recomputed dr2 (Ry) | QE printed (Ry) | agree to printed digits |
|---|---|---|---|
| 1 | 0.438555332 | 0.43855533 | yes |
| 2 | 0.01956403971 | 0.01956404 | yes |
| 3 | 0.0008539880469 | 0.00085399 | yes |
| 4 | 0.001209052555 | 0.00120905 | yes |
| 5 | 2.339192519e-05 | 0.00002339 | yes |
| 6 | 4.953395323e-06 | 0.00000495 | yes |

| check | value | meaning |
|---|---|---|
| k2: out(k) vs out(k-1) | 2.236e-02 | successive rho_out |
| k3: out(k) vs out(k-1) | 5.519e-03 | successive rho_out |
| k4: out(k) vs out(k-1) | 6.550e-03 | successive rho_out |
| k5: out(k) vs out(k-1) | 4.275e-03 | successive rho_out |
| k6: out(k) vs out(k-1) | 7.719e-04 | successive rho_out |

## si8-np1

| k | mode | stop line | step-k estimated scf accuracy (Ry) | step-k ethr | step-k avg # Davidson its | step-k total energy (Ry) | '!' energy | wall s |
|---|---|---|---|---|---|---|---|---|
| 1 | mix | convergence NOT achieved after   1 iterations: stopping | 0.43855533 | 1.00E-02 | 3.0 | -62.61942841 |  | 0.12 |
| 1 | out | convergence has been achieved in   1 iterations | 0.43855533 | 1.00E-02 | 3.0 | -62.60391409 | -62.60391409 | 0.12 |
| 2 | mix | convergence NOT achieved after   2 iterations: stopping | 0.01956404 | 1.37E-03 | 1.0 | -62.64131413 |  | 0.13 |
| 2 | out | convergence has been achieved in   2 iterations | 0.01956404 | 1.37E-03 | 1.0 | -62.64125463 | -62.64125463 | 0.13 |
| 3 | mix | convergence NOT achieved after   3 iterations: stopping | 0.00085399 | 6.11E-05 | 4.0 | -62.64319427 |  | 0.17 |
| 3 | out | convergence has been achieved in   3 iterations | 0.00085399 | 6.11E-05 | 4.0 | -62.64316560 | -62.64316560 | 0.15 |
| 4 | mix | convergence NOT achieved after   4 iterations: stopping | 0.00120905 | 2.67E-06 | 3.0 | -62.64337763 |  | 0.20 |
| 4 | out | convergence has been achieved in   4 iterations | 0.00120905 | 2.67E-06 | 3.0 | -62.64290026 | -62.64290026 | 0.18 |
| 5 | mix | convergence NOT achieved after   5 iterations: stopping | 0.00002339 | 2.67E-06 | 2.0 | -62.64356086 |  | 0.20 |
| 5 | out | convergence has been achieved in   5 iterations | 0.00002339 | 2.67E-06 | 2.0 | -62.64355385 | -62.64355385 | 0.20 |
| 6 | mix | convergence NOT achieved after   6 iterations: stopping | 0.00000495 | 7.31E-08 | 6.0 | -62.64357264 |  | 0.22 |
| 6 | out | convergence has been achieved in   6 iterations | 0.00000495 | 7.31E-08 | 6.0 | -62.64357081 | -62.64357081 | 0.23 |

Trajectory prefix identical across all runs (acc, energy, ethr per iteration): yes

Full trajectory (from the k=6 mix run):

| iteration | estimated scf accuracy (Ry) | ethr | avg # Davidson its | total energy (Ry) |
|---|---|---|---|---|
| 1 | 0.43855533 | 1.00E-02 | 3.0 | -62.61942841 |
| 2 | 0.01956404 | 1.37E-03 | 1.0 | -62.64131413 |
| 3 | 0.00085399 | 6.11E-05 | 4.0 | -62.64319427 |
| 4 | 0.00120905 | 2.67E-06 | 3.0 | -62.64337763 |
| 5 | 0.00002339 | 2.67E-06 | 2.0 | -62.64356086 |
| 6 | 0.00000495 | 7.31E-08 | 6.0 | -62.64357264 |

Density files: gamma_only=1 ngm_g=11936 nspin=1 G0=0.0296178744919
Miller index list identical in every file: yes

### Relative differences ||a-b|| / ||b|| (2-norm over the stored half-sphere coefficients)

| check | value | meaning |
|---|---|---|
| k1 mix beta=1.0 vs k1 out | 2.429e-18 | first-step mix at beta=1 equals rho_out(1) up to rounding |
| rho_in(1) from beta 0.7 vs from beta 0.3 | 1.766e-16 | reconstruction rho_in(1) = (rho_in(2) - beta rho_out(1))/(1-beta) is consistent |
| ||rho_out(1) - rho_in(1)|| / ||rho_in(1)|| | 1.968e-01 | size of the first update |
| rho_in(1) G0 | 0.0296178744919 | electron count / Omega |
| k1: mix(k) vs out(k) | 5.578e-02 | rho_in(k+1) vs rho_out(k) |
| k2: mix(k) vs out(k) | 7.079e-03 | rho_in(k+1) vs rho_out(k) |
| k3: mix(k) vs out(k) | 1.799e-03 | rho_in(k+1) vs rho_out(k) |
| k4: mix(k) vs out(k) | 4.035e-03 | rho_in(k+1) vs rho_out(k) |
| k5: mix(k) vs out(k) | 5.460e-04 | rho_in(k+1) vs rho_out(k) |
| k6: mix(k) vs out(k) | 2.856e-04 | rho_in(k+1) vs rho_out(k) |
| k2: mix(k) vs linear rho_in(k)+0.7(rho_out(k)-rho_in(k)) | 7.042e-03 | 0 would mean plain linear mixing; nonzero = Broyden history |
| k3: mix(k) vs linear rho_in(k)+0.7(rho_out(k)-rho_in(k)) | 3.565e-04 | 0 would mean plain linear mixing; nonzero = Broyden history |
| k4: mix(k) vs linear rho_in(k)+0.7(rho_out(k)-rho_in(k)) | 2.314e-03 | 0 would mean plain linear mixing; nonzero = Broyden history |
| k5: mix(k) vs linear rho_in(k)+0.7(rho_out(k)-rho_in(k)) | 2.823e-04 | 0 would mean plain linear mixing; nonzero = Broyden history |
| k6: mix(k) vs linear rho_in(k)+0.7(rho_out(k)-rho_in(k)) | 1.538e-04 | 0 would mean plain linear mixing; nonzero = Broyden history |

### Recomputed 'estimated scf accuracy' from the stored files (decides what each file holds)

dr2(k) recomputed with QE's own metric from rho_in(k) and the k-out file; rho_in(1) = reconstruction above, rho_in(k>=2) = the (k-1)-mix file.
Agreement to every printed digit means: k-out holds rho_out(k) and (k-1)-mix holds rho_in(k), for every k.

| k | recomputed dr2 (Ry) | QE printed (Ry) | agree to printed digits |
|---|---|---|---|
| 1 | 0.438555332 | 0.43855533 | yes |
| 2 | 0.01956403971 | 0.01956404 | yes |
| 3 | 0.0008539880469 | 0.00085399 | yes |
| 4 | 0.001209052555 | 0.00120905 | yes |
| 5 | 2.339192519e-05 | 0.00002339 | yes |
| 6 | 4.953395323e-06 | 0.00000495 | yes |

| check | value | meaning |
|---|---|---|
| k2: out(k) vs out(k-1) | 2.236e-02 | successive rho_out |
| k3: out(k) vs out(k-1) | 5.519e-03 | successive rho_out |
| k4: out(k) vs out(k-1) | 6.550e-03 | successive rho_out |
| k5: out(k) vs out(k-1) | 4.275e-03 | successive rho_out |
| k6: out(k) vs out(k-1) | 7.719e-04 | successive rho_out |
