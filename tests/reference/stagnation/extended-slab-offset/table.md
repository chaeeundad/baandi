# Measured sweep

All energies are Ha; residuals are Bohr^-3/2. Full windows/audit in JSON.

| Model | Policy | Arm | State | Iterations | Final residual | Certified energy | P11D triggered rules | Recovery ΔE | Sampled gap | Validated intervention |
|---|---|---|---|---:|---:|---:|---|---:|---:|---|
| slab-offset | reference | reference | Converged | 32 | 4.786518210819649e-08 | -0.721558424816071 | none | None | 0.03815640571276868 | None |
| slab-offset | linear-1.0-davidson | baseline | Exhausted | 1000 | 0.11227601418830685 | None | Oscillation, Stagnation | None | 0.0798943111243297 | None |
| slab-offset | linear-1.0-davidson | intervention | Converged | 31 | 2.1375552115903586e-08 | -0.7215584248160892 | Stagnation | 1.8207657603852567e-14 | 0.03815640113315731 | True |
