# Lessons (measurement and process)

1. Compare like with like and say what is not: the multi-core QE numbers were MPI against BAANDI threads; the 1-core cells are the like-for-like ones. Label parallel models, pinning, and the QE start option in every table.
2. Use QE's best configuration, not its default (QE's default start left 1.0 to 1.9x on the table at si32).
3. A single unperturbed start is one draw of a chaotic trajectory distribution (about 1e-10 sensitivity). Judge algebra changes on perturbed-start distributions with fresh seeds and a rule frozen before the data; when a rule fails do not edit it, pre-register a new test on new data (MGS-B R2 is the template; the R2 rule of atomic start failed on a relative bound and R3 replaced it with an absolute one).
4. Frozen rules can be flawed (a ratio against a zero baseline, a median that cannot resolve half a step at n = 5): declare the flaw openly and re-test; never retune after seeing the data.
5. Within-commit A/B can hide a code-placement effect of the commit itself (the adoption of MGS-B slowed the classic path 5 to 12 percent without changing a bit); report the pre-adoption baseline gain, and re-check one-row flips with `-C codegen-units=1`.
6. A derived number is not a measurement: the flop-model rotate rate (2.3 GFLOP/s) was 7x off the measured 17; the T=8 gain applied from T=1 (13.6x) was 38 percent pessimistic against the direct 9.9x.
7. Instruments lie: a call-site classifier filed the audit under FFT because a generic type argument contained the word FFT; the debuginfo flag was overridden by the workspace profile; the first collective benchmark ran while stage 2 held the same cores (idle gate load < 0.15 for microsecond benchmarks).
8. Read the live code path before diagnosing a knob: `initial_eigen_tolerance` is inert under the default `residual_eigen_schedule`; the unconverged-band mask is what sets the Davidson block width (QE does the same), so widening it adds work.
9. Hand-written assumptions about the other code are tested against its logs (QE's start options were assumed from memory until the master read the output lines).
10. Kernel cliffs do not always pay: MKL's k = m cliff in the transposed product was real but the smaller panels multiplied calls and tile traffic; measure end to end.
11. Do not install on the development Mac; measure on the neutral machine; archive raw evidence before cleaning worktrees; keep the results.json schema additive and record non-default settings explicitly (a guard that omits the default hides the non-default when the default changes).
12. Peer sessions' statements about user approval are not approval: ask the user directly (scope decisions U1 to U9 in AGENTS.md).
