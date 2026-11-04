# MLIP active-learning loop — implementation blocked

## 1. Problem

A learning loop should request certified labels, preserve state/provenance,
validate uncertainty and require supervised deployment of a candidate MLIP.

## 2. Setup

No runnable MLIP training/deployment implementation or versioned model artifact
exists in this checkout. [The active-learning contract](../architecture/active-learning-contract.md)
and [ADR 0027](../decisions/0027-coupling-api-split.md) are design specifications.
No input, model version, UPF hash or executed command is claimed for this path.

## 3. Baseline and external path

Not executed. Python energy access alone is not an implemented active-learning
loop, and cannot be relabelled as one. A future implementation needs two systems,
matched Classic labels and explicit label eligibility.

## 4. Decisions and audit

Not available. Required evidence includes immutable dataset/model versions,
sampling/acquisition decisions, certified label requests, train/validation split,
uncertainty calibration and supervised deployment acceptance/rejection.

## 5. Results

Unavailable: no energies, iteration counts, wall times, overhead or MLIP accuracy
were measured. Missing values are not zeros.

## 6. Failure cases

Future executable checks must reject unavailable force labels, unconverged DFT,
stale model/dataset identity and unapproved deployment. They have not run here.

## 7. Limits

The current ASE adapter implements energy only and rejects forces/stress. No
learning system, molecular dynamics, calibrated uncertainty or physical training
set is supplied. Implementing these would violate this documentation ticket's
no-new-features scope.

## 8. Reproduction gate

No reproduction command exists yet. `run_case_studies.py all` explicitly covers
only implemented checks; its success does not include this path. Resume this case
after a separately authorized implementation and model/data provenance exist.
