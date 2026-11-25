# ADR 0046: XC small and negative densities

- Status: proposed; opt-in experiment, not a production default.
- Task: G4E-VACUUM-DENSITY.
- Date: 2026-11-25.

## Inventory before this change

All densities are electrons/Bohr^3; energies and potentials are Hartree.

| Evaluator | Negative input | Small positive input | Gradient / potential |
| --- | --- | --- | --- |
| Serial unpolarized LDA-PZ | Reject below -1e-14; roundoff maps to zero inside XC | Zero through 1e-10 | Point-local |
| Serial unpolarized LDA-PW | Same checked-density rule | Zero through 1e-10 | Point-local, xcx |
| Serial unpolarized PBE, both PW92 conventions | Reject every negative sample | Zero local energy/derivatives through 1e-10 | Differentiate original samples; divergence may survive at screened points |
| Serial polarized LDA-PW | Each channel must be nonnegative | Screen when sum of channels <=1e-10 | Does not independently screen an empty channel at a nonempty point |
| Serial polarized PBE | Each channel must be nonnegative | Screen sum of channels and all sigma components together | Original-channel gradients; divergence may survive screening |
| Reference scalar LDA start/maps | Maps use abs(n), including combined core; start uses abs for core/warm starts | Same local LDA screen | Historical signed-density continuation differs from legacy |
| Reference distributed scalar PBE | Reject every negative sample collectively | Same as serial PBE | Original combined-density Fourier gradient |
| Reference distributed polarized PBE | Reject negative combined channels collectively | Same total-density screen as serial spin | Original combined channels forwarded; no channel floor |
| Reference polarized LDA-PW | Serial spin adapter on owned points | Same as serial spin | Point-local |

With NLCC, scalar XC sees n+c and spin XC sees n_up+c/2, n_down+c/2.
The immutable core is never normalized or clipped. The potential integral used
for band double counting contracts **valence**, not valence plus core. Core
potential integrals are diagnostic. Without NLCC the same rules apply with c=0.
The optional historical density-floor model is a separate functional; this ADR
neither silently replaces it nor qualifies its combination with screening.

## Proposed common rule

`BAANDI_XC_SCREEN_NEGATIVE=on` (registry key `xc_screen_negative`) selects an
XC-only projection q_s(n_s)=n_s for n_s>t and zero otherwise, t=1e-10 Bohr^-3.
The default is off and retains historical scientific behavior.

The threshold is an absolute density scale, not a relative threshold or a
material-dependent tolerance. It preserves the adapters' established numerical
cutoff: xcx LDA exchange/PW and PBE exchange have density limits 1e-15 and PBE
correlation 1e-12; 1e-10 gives a two-decade margin over the largest component
limit before division by fractional density powers. This is an explicit
numerical safety margin, not a derived physical density scale or a rigorous
bound on XC error. Negative values of any magnitude fall below t; their
acceptability is decided by integrated charge, **not** by widening t to fit O.
Nonfinite inputs always fail.

Per channel and XC evaluation, record the count of n_s<=t, integral of |n_s|
over these points (electrons), signed total charge used by the gate, and the
minimum input with its global C-order grid index. With core, the denominator
is total valence electrons in engine calls. Reject if a channel's screened
charge exceeds 1e-4 of that total. This is an explicit proposed model-error
budget (0.01% per channel), not an energy-accuracy guarantee or a proof of
valence infeasibility: a large core truncation error can also fail it. For an
empty system, nonzero screened charge fails. The same gate applies at start,
ordinary maps, and verification. Fix cutoff/grid/start or density feasibility
when the budget fails; do not silently renormalize.

No stored density, core, electron count, kinetic term or Hartree input is
changed. q is private to the XC functional. A fully screened point contributes
zero energy density and zero potential. In spin correlation the energy is a
joint function of both channels, so a separate correlation energy assigned to
one channel is not uniquely defined: use e(q_up,q_down), with zero derivative
with respect to a screened input channel.

## Variational consistency and limitations

Define E_screen[n]=dv sum e(q(n), D q(n)), where D is the existing skew-adjoint
spectral derivative. On a region where the mask is unchanged,
V_screen=M*(e_q-D dot e_gradient), M=1[n>t]. Thus **both** gradients and local
energy use q, and the final mask is applied after the complete divergence.
Masking only local derivatives while differentiating the unscreened density
would leave a nonzero nonlocal derivative at screened points; masking that
potential afterwards would not differentiate the stated energy.

This hard cutoff is discontinuous at t. The formula is a functional derivative
only away from mask crossings; identical masks for energy and potential do
not make a hard threshold globally differentiable. Derivative tests must keep
the mask fixed. Verification evaluates the same screened functional with the
same charge gate and includes its actual XC remainder. Passing that certificate
does not certify the unscreened functional or an energy error due to screening.
No new force/stress qualification is claimed. Positive fields with every
channel above t preserve the previous arithmetic; fields near a historical
screen boundary need not be bit-identical.

## Diagnostics and validation scope

Diagnostic sums use fixed global z rows, then a global row-order sum, independent
of rank partition. Reductions communicate row summaries, never a density
payload. Each complete z row must have exactly one owner. Minimum ties choose
the smallest grid index. Records are in-memory and run-keyed; only the CLI
serializes them, outside numerical loops. They describe XC evaluations, which
can include both input and output fields, not just one record per SCF iteration.
The discarded scalar XC used to assemble spin Hartree does not consume a second
screening budget or emit records.

Synthetic tests and real run results must be reported separately. Simulated
ranks exercise production collectives without establishing real MPI execution
on dev08. No QE source was read and no QE comparison is claimed.

## Local evidence and remaining integration boundary

Mac synthetic isolated one-electron ion, 10 Bohr cubic cell, 2.5 Ha wavefunction
cutoff, synthetic Gaussian NLCC core, fixed up=1/down=0, uniform valence start:
both drivers reject negative XC input with screening off and converge with it
on. Simulated ranks 1/3/8 at reference thresholds 1e-8 and 1e-12 Ry give maximum
legacy energy differences 2.3192092690749178e-10 and 1.5654144647214707e-13 Ha,
respectively. The tight reference energies span -0.466715817218312123 to
-0.466715817218203211 Ha. These are synthetic code comparisons, not QE data.
Reproduce: `cargo test -p dft-engine --test refscf g4e_vacuum -- --nocapture`.

The screened PBE directional derivative test gives absolute error
6.915330530432584e-10 Ha at step 1e-5; the threshold mask is fixed. The distributed
spin PBE test compares serial fields and energy to 1e-12 Ha, and requires exact
screening-record equality for identical samples at simulated ranks 1/3/8.
Reproduce: `cargo test -p dft-physics --test vacuum_distributed --test
vacuum_screening -- --nocapture`.

The current allowed evaluator boundaries identify `start` explicitly, but the
shared map functions are called for both ordinary and verification maps without
a stage argument. Their fallback label is honestly `map/verification`. Exact
stage labels and iteration association require caller hooks in the reference
loops and legacy session/collinear drivers, outside this ticket's OWNED_FILES;
these hooks have been requested, not silently added. Run-keyed diagnostics are
available through `xc::screening_records()` and serialized as `xc_screening` in
CLI output. This stage-label limitation prevents full ticket acceptance.

Final Mac gates: workspace check, physics tests (108 passed, 18 ignored),
reference SCF tests (82 passed, 24 ignored), fmt and targeted clippy passed.
`CHECK_NICE=0 CHECK_SCCACHE=0 BAANDI_TEST_TMP=/private/tmp/g4e-gates
scripts/check.sh fast dft-physics dft-engine dft-cli dft-numerics` exited 0:
743 passed, 0 failed, 193 ignored across 69 binaries; build/fmt/clippy passed.
The sandbox denied the informational `ps` load probe. The first gate attempt
also exited 0, despite a `nice` priority warning; the repeat disabled that
priority adjustment explicitly. No numerical gate failure remains. Positive-density scalar LDA-PZ/PW/PBE SCFs preserve every checked
energy/density bit in both drivers; positive spin LDA-PW/PBE core fields also
preserve energy, potential and both integrals exactly.
