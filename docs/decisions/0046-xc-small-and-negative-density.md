# ADR 0046: XC small and negative densities

- Status: accepted; default on (R1B-ROUTING-AND-SCREEN, 2026-11-25).
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

## Accepted common rule

`BAANDI_XC_SCREEN_NEGATIVE=on` (registry key `xc_screen_negative`) selects an
XC-only projection q_s(n_s)=n_s for n_s>t and zero otherwise, t=1e-10 Bohr^-3.
The default is on. `BAANDI_XC_SCREEN_NEGATIVE=off` retains historical
rejection/continuation behavior for regression checks and rollback.

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
charge exceeds 1e-4 of that total. This is an explicit model-error
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
on dev08. No QE source was read. The original local experiments made no QE comparison;
the subsequent dev08 adoption evidence is recorded below.

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
CLI output. This stage-label limitation remains an integration follow-up; default adoption
does not claim it is resolved.

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


## Default adoption: dev08 gen14 (R1B)

The user-supplied dev08 record `/data/sig-dft/gen14/RESULT.md` establishes the
README PBE Gamma case: Si2, PBE override on Si.pz-vbc, 30 Ry wavefunction cutoff,
25^3 FFT grid, 8 bands. With screening off, a Broyden map transiently produces
-5.6e-6 Bohr^-3 at one atom-site grid point and the reference driver rejects it.
With screening on, the run converges in 7 iterations and is VERIFIED; 2 points
are screened across the run, with integrated screened charge <=1e-7 electrons.
E_reference - E_legacy = -6.5e-11 Ha; the difference from QE is unchanged from
the existing README value. This is the supplied remote evidence, not a new
local QE run or a new independent QE reference value.

Adopt the existing general density threshold and charge gate as the default;
no material-specific cutoff or tolerance changes. The off switch remains
available. Positive densities above the threshold retain their arithmetic;
local reference/oracle identity checks are required for this default change.
The hard-mask, force/stress, and diagnostic stage-label limitations above stand.

### Local R1B baseline check

On macOS/aarch64 with the pure-Rust release build, the three strict-oracle
cases (Gamma, k222s, k222) reproduce `baseline-oracle.json` exactly in energy
and density hash. With the new default, reference k222s and k222 also reproduce
`baseline-reference.json` exactly. Reference Gamma does not: E changes from
-7.295530466369128 to -7.295530466374794 Ha (delta -5.666578317686799e-12 Ha),
and the density SHA-256 changes from
`35295c543cfde3e7303c2edfa97a47d4d6999f52607750f355b8145c88e471d7` to
`a0c76c8c7741125ad2cf59dee1534a0e4c29bc306002e6de22e56b92c99573e2`.
This is a signed-density trajectory: evaluation index 5 has minimum
-8.4477852131698e-5 Bohr^-3 at grid index 0, one screened point and
1.4603566607452055e-6 electrons screened. All other evaluations screen zero
points. Explicit `threads=1` reproduces the same difference;
`BAANDI_XC_SCREEN_NEGATIVE=off` restores the old Gamma energy and density hash
exactly. The architect approved one screening rule for LDA and GGA alike, without a
signed-density exception. In particular, the default replaces the reference
scalar LDA map's historical evaluation at |n| with ADR 0046 screening. The
Gamma energy difference is an intentional policy change, not a relaxation of
the 1e-12 Ha regression tolerance. Only the Gamma entry of
`docs/performance/baseline-reference.json` is re-recorded; k222s/k222 and
`baseline-oracle.json` remain unchanged. Screening off reproduces the old
Gamma baseline exactly in energy and density hash.

The Gamma entry was recorded with the pure-Rust release build using
`python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf <UPF>
--cases gamma --threads 1 --repeats 3 --record <temporary-record.json>`;
only `cases.gamma` was copied into the baseline. The recording command exits 5
against the old baseline, as expected for the measured policy change.
The external Si.pz-vbc.UPF SHA-256 is
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`.
These are local arithmetic checks, not speed-adoption evidence or new QE
measurements. The parent baseline's provenance describes the untouched cases;
this paragraph records the replacement Gamma provenance.

The approved source-local test migrations pin the old negative-density
rejection contract to `xc_screen_negative=off` and make the spin/FD fallback
fixture's driver request implicit. They change no production code.

After those migrations, `cargo check --workspace --tests` and
`cargo fmt --check && cargo clippy --workspace --all-targets -- -D warnings`
exit 0. `CHECK_NICE=0 CHECK_SCCACHE=0
BAANDI_TEST_TMP=/private/tmp/r1b-approved-gates scripts/check.sh` exits 0:
108 binaries, 1024 passed, 0 failed, 224 ignored; build/test/fmt/clippy all pass.
The sandbox denies the informational `ps` probe, without preventing any gate.
The explicit routing/corpus/reference_identity targets pass (8 + 1 + 1 tests).
The benchmark command with `--cases gamma,k222s,k222 --require-identical`
exits 0 against the updated reference baseline; adding
`BAANDI_REFERENCE_STRICT=on` and `--baseline docs/performance/baseline-oracle.json`
also exits 0. Every case has zero energy difference and identical density SHA-256.
The Gamma-only `BAANDI_XC_SCREEN_NEGATIVE=off` check against the saved original
baseline exits 0. No oracle or non-Gamma baseline entry was changed.
