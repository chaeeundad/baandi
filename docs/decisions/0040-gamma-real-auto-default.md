# ADR 0040: scf.gamma_real defaults to auto

Status: Accepted (user, 2026-11-09: "오케이" to the proposal "change the default of the real Gamma path from off to auto, recommendation: change")

## Context
The real Gamma-only path (half-sphere ORC basis, paired FFT, real Davidson; tickets P26 to P39) was implemented and validated against the complex path but shipped with `gamma_real = off`. Every dev08 speed figure in recent reports was measured with the real path switched on; a user running defaults got the complex path, which is 1.5 to 2 times slower on Gamma cells (dev08 88cb469: si8 standard 4.15x QE complex against 2.94x real, si16 5.30x against 2.84x, si32 strict 17.87x against 7.29x). An external review (2026-11-09) pointed at the same default.

## Decision
`ScfOptions::default().gamma_real` and the CLI key `scf.gamma_real` default to `auto`. `auto` selects the real path when the problem qualifies (single Gamma point, fixed occupations, unpolarized, no unsupported solver options) and otherwise runs the complex path and records the reason in `results.json` (`gamma_real_used`, `gamma_real_fallback_reason`). `off` and `on` keep their meaning (`on` fails closed). Old checkpoints without the key are read as `off`, so exact continuation of old runs keeps its path.

## Consequences
- Default results move at rounding level (energies within the validated 1e-10 class; density hashes and SCF counts differ from the complex path); strict certification is unchanged.
- results.json carries two extra keys by default (26 instead of 24 top-level keys on the small fixtures).
- Test expectations that implicitly assumed the complex path were updated or pinned to `off` (spin and spectral-guard comparisons, legacy checkpoint header); baselines for "default path identical" comparisons now use `gamma_real = "off"` explicitly.
- Cases that fall back to the complex path (k points, spin, Fermi-Dirac, locking, RMM, spectral guard) behave as before.

## Evidence
P26 to P39 validation series and `docs/performance/real-gamma-*.md`; dev08 relays in `tests/reference/perf/` (88cb469 re-measurement, stage decomposition). A dev08 comparison of the new default against `off` was requested after the change.
