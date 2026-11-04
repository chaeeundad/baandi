# Self-Healing SCF: bounded iteration-budget repair

Status: executed **partial case**. This is budget recovery, not numerical-instability healing.

## 1. Problem and motivation

An SCF budget of one iteration cannot establish the energy-change certificate.
A host can permit a bounded budget repair without letting an optional agent
change the physical model or convergence thresholds.

## 2. Setup and exact execution

The public synthetic fixture is `common::potential(4.0, 1.0, true, 801)` in
`crates/dft-engine/tests/common/mod.rs`: two Xx ions at fractional (0,0,0)
and (1/4,1/4,1/4), primitive cell with a=10.261212856717933 Bohr, Gamma,
5 Ha cutoff, LDA-PZ, Pulay beta=0.3/history=6, one worker, seed=2026.
The fixture is generated in memory, not a real-element potential. No UPF file
SHA-256 was exported by this harness; generator/source hashes and expanded
options are recorded instead. That provenance gap is explicit.
Build/date/commands: [provenance](data/recorded/provenance.json),
[build](data/recorded/self-healing/build.txt).

## 3. Baseline and external path

Four arms, each executed twice from a fresh identical state: one-step baseline;
one-step baseline repaired by a deterministic rule; the same repair proposed by
a scripted recorded agent; and a fresh 1000-step reference. All four use
AgentControlled. **There is no Classic comparator in this harness.**
It therefore cannot measure overhead versus Classic or satisfy that acceptance
item. The rules/agent grant allows one budget extension to at most 200.

## 4. Decisions and evidence

Initialization supplies state_version=1, policy_revision=0. At iteration zero the
rule or recorded decision submits `extend_budget` to 200 with exact observation
and expiry. The host-installed grant supplies authority; the agent cannot grant
itself authority. The agent arm first probes invalid beta=1.1, which is rejected,
then applies the budget request. Neither tolerance nor model changes.
See [agent audit](data/recorded/self-healing/agent-0/audit.jsonl),
[events](data/recorded/self-healing/agent-0/events.jsonl) and
[summary](data/recorded/self-healing/agent-0/summary.json).
The actual audit distinguishes accepted/applied/rejected outcomes and retains
request debug text, authority, before/after policy and a hash chain. This is a
supervised scripted replay, not a live LLM or asynchronous human review.

## 5. Results

See the measured table below and [aggregate](data/recorded/self-healing/results.json).
The harness checks equal initial density hashes, replay equality, full normalized
result/array identity and certified energy agreement within 1e-10 Ha. Baseline
energy is retained but uncertified; its certified energy is null. Wall times
include harness observation, and no isolated observer overhead is available.

| Arm (repeat 0) | Certified energy (Ha) | Iterations | Applied / rejected | Wall (s) |
| --- | ---: | ---: | ---: | ---: |
| baseline | unavailable | 1 | 0 / 0 | 0.001724 |
| rules | -4.684561184277879 | 30 | 1 / 0 | 0.044914 |
| agent | -4.684561184277879 | 30 | 1 / 1 | 0.041499 |
| reference | -4.684561184277879 | 30 | 0 / 0 | 0.035148 |

Rules/agent/reference energy differences are 0 Ha. Each converged arm records 3 certificate solves and 951 operator applications including certificates.

## 6. Failures

The baseline exhausts its budget. The agent's invalid beta is rejected without
changing the model. Separately executed tests reject stale/future observations
and unknown authority fields; see [negative-control log](data/recorded/self-healing-rejections.log).
Those tests do not emit a persistent rejection transcript. No stagnation or
oscillation rescue is demonstrated.

## 7. Limits

One two-ion synthetic model, one Mac, one backend/thread setting; no generality,
physical silicon, QE, distributed execution or ground-state uniqueness claim.
A second model and a matched Classic arm remain necessary. Existing synthetic
constructors could support these, but the current harness does not expose them.

## 8. One-command reproduction

```sh
python3 -I scripts/run_case_studies.py self-healing --verify
```
