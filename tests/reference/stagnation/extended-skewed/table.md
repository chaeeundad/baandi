# Measured sweep

All energies are Ha; residuals are Bohr^-3/2. Full windows/audit in JSON.

| Model | Policy | Arm | State | Iterations | Final residual | Certified energy | P11D triggered rules | Recovery ΔE | Sampled gap | Validated intervention |
|---|---|---|---|---:|---:|---:|---|---:|---:|---|
| skewed | reference | reference | Converged | 19 | 5.796688160226345e-08 | -0.6441418148101661 | none | None | 0.024375866177302696 | None |
| skewed | linear-1.0-davidson | baseline | Exhausted | 1000 | 0.08623068394761953 | None | Oscillation, Stagnation | None | 0.0656825395304168 | None |
| skewed | linear-1.0-davidson | intervention | Converged | 39 | 6.667116061439711e-08 | -0.6441418148101419 | Stagnation | 2.4202861936828413e-14 | 0.02437591106872128 | True |
