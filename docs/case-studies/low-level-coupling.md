# Low-level external-solver or QM/MM coupling — implementation blocked

## 1. Problem

Scientific coupling must distinguish electron potential energy from electrostatic
potential, and inner SCF convergence from outer coupled convergence.

## 2. Setup

[The low-level API](../architecture/coupling-api-lowlevel.md) and
[ADR 0027](../decisions/0027-coupling-api-split.md) specify a design, not a merged
public field/operator exchange implementation. No external solver version,
physical coupling model, UPF hash or executable input can be recorded yet.

## 3. Baseline and external path

Not executed. An internal Hamiltonian setter or Hartree solve is not the public
coupling API. A future case needs uncoupled Classic, zero-coupling equivalence
and a defined nonzero coupling with independent outer convergence checks.

## 4. Decisions and audit

No executed exchange or audit exists. Future evidence must retain sign/units,
cell/grid/origin identities, immutable environment identity, operator/buffer
lifetimes, request versions, accepted/rejected exchanges and energy ownership.
Optional agent advice must not replace numerical convergence certification.

## 5. Results

Unavailable: no coupled energies, outer residuals, iteration counts, wall times
or transport overhead were measured. No QM/MM capability is claimed.

## 6. Failure cases

Required future checks include stale environment/grid identities, wrong units,
unsupported boundary conditions, buffer expiry and incomplete energy ledgers.
These are pending requirements, not successful test results.

## 7. Limits

No API implementation, external solver driver, partition/link-atom model or
validated coupling derivatives exists for this case. Forces/stress and coupled
total-energy claims remain separately gated. New coupling features are outside
P11L ownership and scope.

## 8. Reproduction gate

No one-command reproduction is currently possible. Resume after the coupling API
and a supported physical model are implemented under a separate ticket. The
available-evidence runner does not count this design page as a passing case.
