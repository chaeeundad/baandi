# ADR 0008: UPF reading and radial mathematics

Status: Accepted (2026-11-02, Claude Code). Applies ADR 0003.

## Decisions
- **UPF text -> typed structs:** the `upf` crate 0.1.1 (MIT OR Apache-2.0, repository ShapelessCat/upf, released 2026-04). It is new and little used (71 downloads when checked), so it is used for reading only. Our `dft-upf` layer validates every field the engine relies on and rejects anything outside the supported subset with explicit errors (no silent skipping). Its dependency tree (quick-xml, serde, serde_with, validify, env_logger, anyhow...) was license-checked: no GPL-only crate. Re-evaluate if a field is mis-parsed; the fallback is a small `quick-xml` reader of our own. Authoritative format reference: the UPF specification shipped inside the crate (`docs/reference/upf-spec.html`) and the public UPF format page; no QE (GPL) source is read or translated.
- **Spherical Bessel functions j_l and real spherical harmonics Y_lm (l <= 3):** hand-written in `dft-upf` (small, elementary, testable by identities). Libraries checked: `scilib` is GPL-3.0 (rejected); `scirs2-special` is Apache-2.0 but belongs to a large ecosystem and was not trialed. Verification does not depend on a library: series and closed-form agreement, recurrence relations, and the addition theorem sum_m Y_lm(a) Y_lm(b) = (2l+1)/(4*pi) P_l(a.b).
- **Radial integration:** Simpson-type quadrature on the UPF mesh using the supplied `rab` weights, in `dft-upf`; the Coulomb tail of the local potential is handled analytically (erf subtraction) rather than by truncation.

## Consequences
UPF units (Rydberg for potentials and D_ij in UPF v2, projector storage as r*beta(r)) must be converted exactly once, inside `dft-upf`, and covered by tests. Cross-checks against the real Si file wait for the user's approval of that file.
