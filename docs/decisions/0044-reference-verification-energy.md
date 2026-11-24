# ADR 0044: Reference-driver verification by energy

Accepted: user decision 2026-11-25, option B; ticket I25-VERIFY-ENERGY.

**Status update 2026-11-25 (user decision C): implemented but default OFF.** dev08 MPI8 one-run screen at
7b2c624 (P3f): the verification map's energy change tracks the Coulomb fixed-point error (ratio about
0.85 at every attempt), so the energy rule rarely passes earlier: si64 7.62 -> 7.28 s (one map saved),
si128 38.39 -> 40.64 s (three attempts instead of two because VERIFY_GATE is inactive in energy mode),
elongated 1.51 -> 1.69 s, rskew equal; delivered energies within 5e-11 Ha of the strict reference in
both modes. The reference driver keeps the Coulomb fixed-point rule with REFSCF_VERIFY_GATE; the energy
rule stays available as BAANDI_REFSCF_VERIFY_ENERGY=on. Reporting rule (decision C): every QE comparison
of the reference driver shows two BAANDI numbers, verification off (same conditions as QE, which runs no
independent check) and verification on (default), see docs/performance/comparison-rules.md.
This supersedes the reference driver's Coulomb acceptance rule only. ADR 0043
and the legacy SCF driver retain their existing behavior.

## Decision and energy definition

`BAANDI_REFSCF_VERIFY_ENERGY=on` (default; input
`[performance] refscf_verify_energy = true`) judges the existing independent
verification map by `abs(E_check - E_stop) < 5e-9 Ha`. The threshold is fixed,
not inherited from a custom stopping-estimator threshold. With verification
enabled, a fresh map always runs at every estimator stop; no predictor skips
it. Explicit `verification = Off` remains an opt-out, not a verified result.

For a map with input density x, eigenvalues epsilon_i and output density y,
both energies use the same existing variational band form:

    E(x -> y) = sum_i f_i epsilon_i[H(x)]
                - integral y(r) (V_H[x](r) + V_xc[x](r)) dr
                + E_H[y] + E_xc[y] + E_ion-ion.

The band sum includes kinetic, local ionic (including its zero mode), nonlocal
pseudopotential, and input-density Hartree/XC potentials. Subtracting the
output-density integral removes the input Hartree/XC contribution; Hartree
and XC energies are then evaluated at the output density. Ion-ion energy is
the same Ewald term in both maps. Occupations are the driver's fixed integer
occupations. The existing Hartree G=0 convention and LDA evaluation at |rho|
are unchanged. These are `Map.energy.band_form` values, not Harris-Foulkes
energies or the separately audited kinetic/local/nonlocal sum.

At the stop, `E_stop = E(rho_in -> rho_out)`. The verification uses
`E_check = E(rho_out -> F(rho_out))`. No mixed-next-density `descf` correction
is included in either value. At successful delivery this is exactly the
reported `variational_energy_ha` (`descf = 0`). A rejected attempt records its
candidate energy separately: the later ordinary-iteration variational report
may include the existing mixed-density correction.

The verification tolerance stays `next_ethr_ry(ethr, dr2, electrons) / 2` Ha,
with the existing warm orbitals and inner-solver retry handling. Nonfinite
energies fail the energy judgment. Coulomb fixed-point error, the empirical
residual estimate (`residual_bound_ha`), density residual and strict diagnostic
criteria are still computed and recorded; none is an additional energy-mode
acceptance gate. Failed attempts retain the existing Broyden recovery and
optional failed-map reuse. All attempts and recovery remain charged. The
delivered density, orbitals, eigenvalues and energy audit remain the stop
map's own state, never the verification map's state.

`REFSCF_VERIFY_GATE` and `REFSCF_VERIFY_GATE_FIRST` are inactive in energy
mode, including the first stop. Their configuration values remain visible,
with an explicit inactive path note; verification JSON records the energy
criterion and inactive gate. `REFSCF_VERIFY_REUSE` still applies after failure.

## Evidence motivating the decision

User-supplied dev08 issue-25 evidence (not reproduced on this Mac): MPI8 dr2
tracks QE's estimated accuracy map by map, stopping at si128 map 15 and si64
map 11. The old fixed-point check gives 1.51e-7 Ha and 5.80e-8 Ha respectively,
12–30 times the 5e-9 Ha threshold. Two additional maps plus attempts cost about
25% of si128 wall. Stop-map energies differ from strict references by
1.0e-9 to 3.3e-9 Ha across si64/si128/si216 np8 and si32/si64 np1. This ADR
records the supplied rationale, not a newly generated QE reference bundle or
an independently reproduced timing result.

## Guarantees, limits and rollback

With verification enabled, successful delivery requires an independent map
and agreement of the delivered candidate's variational energy with that fresh
map to less than 5e-9 Ha. This is an energy-change check, not a mathematical
bound on the error against the exact self-consistent solution. The Coulomb
fixed-point error remains observable but is no longer enforced to 5e-9 Ha.
Density fixed-point accuracy, forces and stress are not certified by this rule.

Set `BAANDI_REFSCF_VERIFY_ENERGY=off` (or the input value to false) to restore
the old Coulomb-plus-residual-estimate rule and the configured predictor gates.
The off branch preserves the old numerical operations, state and work counts.
The legacy driver does not consult this switch. The pre-existing console banner
retains its historical `fixed_point` label (this ticket owns verification JSON
only); the JSON `mode` and `criterion` identify the actual judgment. Local checks are restricted to
the ticket's debug checks/tests; real MPI/MKL and delivered-wall measurements
remain work for the master dev08 session.

## Local validation and observed limits (Mac)

The existing synthetic skew-cell fixtures use cutoff 2 Ha, one atom at
(0.13, 0.21, 0.32), Gaussian atomic start, SLA/PZ, and default beta 0.7.
The broad local fixture has Z=2, rc=1.8 Bohr, length 7 Bohr, two bands,
history 8. The projector fixture has Z=4, rc=1.2 Bohr, length 11 Bohr,
three bands, history 4. Both have row shears 0.2 and 0.1 Bohr.

| Fixture | Old / energy final map | Old / energy attempts | First energy change (Ha) | Final energy change (Ha) | Delivered energy difference from old (Ha) |
| --- | --- | --- | --- | --- | --- |
| Broad local | 6 / 6 | 1 / 1 | 2.264195497758692e-10 | 2.264195497758692e-10 | 0 |
| Projector | 28 / 27 | 2 / 2 | 7.085897646419426e-9 | 2.6875612846311014e-10 | 1.576103692002562e-10 |

The projector's first Coulomb diagnostic is 1.783487664218233e-8 Ha.
**The requested first-attempt improvement is not demonstrated by these
existing fixtures:** the first energy difference also exceeds 5e-9 Ha and
correctly fails. The energy mode disables the predictor's skipped attempt and
converges one ordinary map earlier, with two genuine attempts. Tests retain
this failure/recovery case rather than relaxing the energy threshold.

A temporary before/after fingerprint check covered all 16 old-rule fixture
and gate/reuse combinations (two fixtures, eight variants). Scientific
records, energy audits, density/orbital bits, eigenvalues, occupations and
work counts matched exactly, excluding clocks/optional inner diagnostics
and the new absent energy-record field. The permanent old-rule matrix keeps
its pre-existing counts and comparisons; it explicitly switches energy off.

Exploratory variations of these small fixtures (history, beta and cell length
or cubic shape; no production settings changed) did not establish the missing
first-attempt improvement either. They also show why this is not an error
bound against the old rule: on the projector fixture with history 3, the old
rule accepts map 18 (Coulomb 2.042367798543282e-9 Ha), but its energy change is
8.178304700123817e-9 Ha. Energy verification instead accepts map 19, and the
delivered energies differ by 8.221326286417252e-9 Ha, exceeding 5e-9 Ha. No
strict oracle was run to decide which delivered energy is closer. Some broader
Gaussian variants exhausted 100 iterations under both rules; identical
unconverged outputs are not accuracy evidence. These exploratory outcomes
are limitations, not reasons to change a tolerance or adopt a mixing setting.
The original two fixtures satisfy the requested old/new energy agreement;
universal old/new agreement and the dev08 first-attempt claim remain unproven.

Permanent tests cover first-pass delivery, real failure/recovery, explicit
verification off, unchanged candidate state, inactive gates including the
first gate, fixed threshold under a deliberately loose estimator, strict
inequality/nonfinite rejection, diagnostic levels, and actual CLI JSON plus
`refscf-check --stop-to-verified` in both judgment modes. Run the ticket's
`cargo check --workspace --tests`, `cargo test -p dft-engine --test refscf`,
`cargo test -p dft-cli stop_to_verified`, `cargo test -p dft-numerics switches`,
and fmt/clippy commands. No release matrices, full check.sh, remote run or
performance acceptance was performed by this implementer.
