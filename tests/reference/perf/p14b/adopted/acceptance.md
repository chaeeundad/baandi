# P14B architect-adopted incremental Gram

Realistic fixtures retain SCF counts and |ΔE| ≤ 1e-12 Ha. Converged toy fixtures retain the suite energy tolerance (1e-12 Ha); iteration changes are reported. Existing failures are retained, not counted as converged. No new convergence failures are permitted.

| Case | Before SCF | Adopted SCF | Δ iterations | ΔE (Ha) | Status |
| --- | ---: | ---: | ---: | ---: | --- |
| si-elongated | 100 | 100 | +0 | 4.603995e-02 | nonconverged |
| si-gamma | 10 | 10 | +0 | 4.440892e-15 | pass |
| si-k222 | 10 | 10 | +0 | 1.776357e-15 | pass |
| si-k222s | 10 | 10 | +0 | -3.552714e-15 | pass |
| si-low | 11 | 11 | +0 | -6.217249e-15 | pass |
| si-skew | 13 | 13 | +0 | -1.776357e-15 | pass |
| si8-gamma | 15 | 15 | +0 | -1.421085e-14 | pass |
| syn-cubic-degenerate | 10 | 10 | +0 | 2.664535e-15 | pass |
| syn-d3 | 8 | 8 | +0 | 1.776357e-15 | pass |
| syn-fft-primes | 14 | 16 | +2 | 0.000000e+00 | pass |
| syn-gaussian-fcc | 17 | 17 | +0 | -8.881784e-16 | pass |
| syn-odd-shift | 9 | 9 | +0 | 2.664535e-15 | pass |
| syn-p3 | 8 | 10 | +2 | -7.105427e-15 | pass |
| syn-s1 | 13 | 14 | +1 | 5.329071e-15 | pass |
| syn-s2 | 17 | 16 | -1 | -8.881784e-16 | pass |
| syn-s3 | 100 | 100 | +0 | 4.944725e-02 | nonconverged |
| syn-spd | 11 | 11 | +0 | 0.000000e+00 | pass |
| syn-z1 | 7 | 7 | +0 | 1.110223e-16 | pass |
| syn-z14 | 100 | 100 | +0 | -1.027704e-02 | nonconverged |

All 19 suite cases have bit-identical energy, density hashes, eigenvalues and H counts at T=1,2,4,8, including retained nonconverged outputs. The legacy harness exits 5 for accepted toy iteration changes at T=1; its cross-thread comparisons also exit 5 because thread counts are part of input identity. This independent test checks scientific determinism directly.

## Additional realistic fixtures

| Case | Before / adopted SCF | ΔE (Ha) | T=1,2,4,8 identical |
| --- | --- | ---: | --- |
| si8 | 15 / 15 | -1.421085e-14 | yes |
| k444s | 10 / 10 | 0.000000e+00 | yes |
| syn20-skew | 6 / 6 | -2.842171e-13 | yes |
| syn20-cubic | 100 / 100 | 1.446232e-04 | yes |
| syn4-elongated | 9 / 9 | 1.776357e-15 | yes |
| syn20-cubic12 | 13 / 13 | 9.379164e-13 | yes |

The dense 8-Bohr cubic 20-atom fixture is nonconverged before and after (100 iterations), so its energy is not an accepted converged comparison. The 12-Bohr cubic fixture converges in 13 iterations in both builds; its H-vector work increases 1996→2107, despite unchanged SCF count. Both input definitions and both outcomes are retained.

## Adopted GEMM accounting and next-site bound

| Case | Before 2mnk | Adopted 2mnk | Reduction % | Ritz/residual 2mnk | Second projection pass 2mnk |
| --- | ---: | ---: | ---: | ---: | ---: |
| si8 | 21636961594 | 9489226838 | 56.14 | 5386241040 | 1296609804 |
| k444s | 3300404692 | 2258007134 | 31.58 | 946782528 | 400012164 |
| syn20-skew | 2846384 | 2426462 | 14.75 | 1605888 | 156456 |
| syn4-elongated | 612150 | 496818 | 18.84 | 258720 | 65352 |

Ritz/residual products are the next-largest site. For the current unlocked convergence test, every root residual must be recomputed: guaranteed removable 2mnk is zero without changing the algorithm. The listed full site totals are ceilings on potential savings, not attainable savings. Second-pass projection totals are similarly ceilings; testing the normalized block still needs its second overlap, and omitting the subtract/MGS would require a new orthogonality gate. No second optimization is adopted.

MKL executes the same two rectangular products with borrowed contiguous input columns and one sequential cblas call per rectangle. It adds only reduced-matrix copies: 16*p*p bytes to preserve the prefix plus 16*(q*q-p*p) bytes to scatter outputs, i.e. 16*q*q bytes per extension. This is an explicit small-copy increase versus full direct-output MKL GEMM, not a claim of reduced copied bytes. MKL runtime/determinism/timings require dev08.
