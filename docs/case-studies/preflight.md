# Agent-assisted preflight: review without silently changing a calculation

## 1. Problem and motivation

An optional adviser can suggest a correction, but scientific inputs must remain
under operator control. This case records both operator decisions and a stale
proposal failure, then checks that review alone leaves the Classic calculation
unchanged.

## 2. Setup

Two generated synthetic fixtures: H-labelled cubic two-ion valence-one/s-projector
and He-labelled skew one-ion valence-two/s,p-projector cells. The cubic cell is
4 Bohr on each side; the skew rows are (4,0,0), (.5,4,0), (0,.5,5) Bohr.
Positions are (0,0,0), (.5,.5,.5) for the first and origin for the second.
Both use Gamma, LDA-PZ, 8 Ry cutoff and four bands. Generator radius=1.4 Bohr,
401 points. These are nonphysical test pseudopotentials, not material references.
Exact options appear in the stored input/results; UPF SHA-256 values are in
[preflight-results.json](data/recorded/preflight-results.json). Compiler, build,
source hashes, UTC time, every command and exit code are in
[provenance.json](data/recorded/provenance.json).

## 3. Baseline and external path

Run Classic on the original input. The deliberately oversized advisory
`minimum_distance_bohr=100.0` generates PF_MIN_DISTANCE. An external scripted
adviser proposes changing only this review threshold to 0.1 Bohr. The proposal
is previewed; explicit approved and rejected decisions are recorded separately.
Neither decision applies the candidate. Run Classic again with identical input
bytes and require convergence, equal iteration count and bit-identical energy.
This tests review-only invariance, not an agent choosing a physically correct
bond length or improving an unconverged model.

## 4. Actions, decisions and audit

For each system the runner preserves original input, raw findings stdout,
full-file unified diff, source/findings/candidate/UPF hashes, exact TOML edit,
and declared finding code. The actor is explicitly `P11L-scripted-operator`.
The frozen script records both decision branches; no live human or LLM approval
is claimed. A real operator should review the diff before recording a decision.
See [cubic proposal](data/recorded/cubic_h2/proposal.json) and the `approved/` and
`rejected/` sidecars under [recorded data](data/recorded/).
The sidecars retain proposal source, decision reason, candidate findings and
`applied: false`. The original input is checked byte-for-byte afterward.

## 5. Results

The table below records newly executed pairs. Both paths retain certified energies
and identical iterations. Timing records separate Classic run, review sequence
(check + preview + both decisions), and reviewed run. Workflow overhead is their
signed difference `(review + reviewed run) - Classic`; the stale negative control
is excluded. Single sequential observations include process/file I/O and cache
variation; negative overhead is not evidence of acceleration.

| System | Energy (Ha) | Iterations, each run | Delta (Ha) | Classic (s) | Review (s) | Reviewed run (s) | Signed workflow difference (s) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| cubic_h2 | -0.9291990923945251 | 12 | 0.0 | 0.520104 | 0.090478 | 0.073767 | -0.355858 |
| skew_he | -1.4687372786010826 | 9 | 0.0 | 0.079492 | 0.100643 | 0.079536 | 0.100687 |

## 6. Failures and stale requests

A separate proposal replaces the original input SHA-256 with zeros. Both systems
reject that proposal with exit code 1; logs are `cubic_h2-stale.log` and
`skew_he-stale.log`. No input is modified. Eleven additional integration tests
exercise malformed/undeclared changes, changed source or UPF bytes, invalid
candidate approval, duplicate decisions, symlink destinations, periodic geometry
and bounded searches. See [preflight.log](data/recorded/preflight.log).

## 7. Limits

Two small synthetic cells, one machine; no actual input edit applied, no live
agent, no optimization of scientific settings, no proof of physical convergence
versus cutoff or cell size. The proposal reduces an intentionally artificial
review threshold; it is not advice for a real material. Changing cutoff/geometry
would change the model and require a separate comparison. Original UPF paths
expire after execution; hashes are preserved and the runner regenerates files.

## 8. One-command reproduction

```sh
python3 -I scripts/run_case_studies.py preflight --verify
```

No Python extension or ASE is needed for this case. Expected stale failures are
checked as exit 1, not treated as a runner failure. Numerical rerun tolerance is
1e-10 Ha; each current pair must match bitwise and converge.
