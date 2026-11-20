# LANE-VERIFY first-attempt schedule — dev08 request

Candidate: `BAANDI_VERIFY_SCHEDULE=first`; control: `on` (current default).
**Conditional diagnostic request, not an adoption recommendation.** The final
Mac Si64 work check retained two attempts and increased maps 13->14 and total
H 6,644->6,661. The 48-run generality matrix reduces aggregate verification
work but adds ordinary maps and has wrong-skip controls. Thus the intended
Si64 mechanism is unproven; do not present the optimistic failed-attempt cost
as an expected gain. MPI8 and hardware timing remain unmeasured for this patch.
Same binary, all other switches identical. `STOP_PREDICTOR=off` is pinned in
both arms: this experiment must not change the acceptance rule. No adoption
or performance conclusion comes from Mac durations. See
[first-schedule.md](first-schedule.md) for the model and its false-positive
limitations. Run through the orchestrator, one case/model per job, no competing
build or benchmark. No QE comparison is requested.

## Preconditions and exact commands

Use the accepted dev08 build configuration and the **same archived Phase 0
inputs and UPF** (including geometry, seed, threshold, bands, grid, ladder,
threads=1). Do not regenerate a slightly different geometry for comparison
against the relay. Record full commit/dirty diff, build command, compiler,
MPI/MKL versions, host, affinity, binary/input/UPF SHA-256 and effective switches.
The driver preserves the three hashes, command, inherited BAANDI/MKL/OMP
settings, load observations, scientific JSON and density hash. Record version
and binding evidence separately with the returned archive.

Set these to the existing absolute paths in the orchestrator's dev08 job:

```sh
export BINARY=/absolute/path/to/candidate/dft-cli
export SI_UPF=/absolute/path/to/Si.pz-vbc.UPF
export SI64_INPUT=/absolute/path/to/phase0/si64/input.toml
export SI32_INPUT=/absolute/path/to/phase0/si32/input.toml
export VERIFY_OUT=/absolute/path/to/new/verify-first-results
mkdir -p "$VERIFY_OUT"
```

`dev08_schedule.py` refuses Mac timing. It performs one unretained warmup per
arm, then crossed AB/BA pairs, and never discards an outlier. Run each A/A
first; inspect its floor before starting the corresponding A/B. CPU bindings
below assume the approved same-socket CPUs 0–7/0; if the host's reserved CPUs
differ, replace both arms with the identical approved binding and report it.

MPI8 Si64 (same MPI binary and rank model throughout):

```sh
python3 -I docs/performance/lanes/VERIFY/dev08_schedule.py \
  --binary "$BINARY" --input "$SI64_INPUT" --upf "$SI_UPF" \
  --launcher 'taskset -c 0-7 mpirun -np 8 --bind-to core --map-by core' \
  --mode aa --pairs 6 --out "$VERIFY_OUT/si64-np8-aa"
python3 -I docs/performance/lanes/VERIFY/dev08_schedule.py \
  --binary "$BINARY" --input "$SI64_INPUT" --upf "$SI_UPF" \
  --launcher 'taskset -c 0-7 mpirun -np 8 --bind-to core --map-by core' \
  --mode ab --pairs 12 --out "$VERIFY_OUT/si64-np8-ab"
```

One core Si64 and Si32, each a separate sequential job (use the accepted
one-core binary/configuration, not an eight-rank wall comparison):

```sh
python3 -I docs/performance/lanes/VERIFY/dev08_schedule.py \
  --binary "$BINARY" --input "$SI64_INPUT" --upf "$SI_UPF" \
  --launcher 'taskset -c 0' --mode aa --pairs 6 --out "$VERIFY_OUT/si64-one-aa"
python3 -I docs/performance/lanes/VERIFY/dev08_schedule.py \
  --binary "$BINARY" --input "$SI64_INPUT" --upf "$SI_UPF" \
  --launcher 'taskset -c 0' --mode ab --pairs 12 --out "$VERIFY_OUT/si64-one-ab"
python3 -I docs/performance/lanes/VERIFY/dev08_schedule.py \
  --binary "$BINARY" --input "$SI32_INPUT" --upf "$SI_UPF" \
  --launcher 'taskset -c 0' --mode aa --pairs 6 --out "$VERIFY_OUT/si32-one-aa"
python3 -I docs/performance/lanes/VERIFY/dev08_schedule.py \
  --binary "$BINARY" --input "$SI32_INPUT" --upf "$SI_UPF" \
  --launcher 'taskset -c 0' --mode ab --pairs 12 --out "$VERIFY_OUT/si32-one-ab"
```

For an isolated count/decision diagnostic, run the identical CLI/launcher with
`BAANDI_VERIFY_SCHEDULE=first-trace`, then `first`; retain those outputs outside
the timing pairs. First-trace keeps baseline numerical decisions and supplies
the actual pass/fail label at the proposed first skip. With profiles desired,
use Phase 0 switches only for these separate diagnostics. No profile clocks or
per-rank gathers are enabled in the timed pairs.

## Pre-registered return counters and gates

For **each** retained run return:

- `verification.attempts[]`: iteration, passed, coulomb_error_ha, threshold_ha,
  residual_bound_ha, eigen_tolerance_ha, h_vectors, davidson_iterations.
- `verification.schedule.{mode,prediction_evaluations,records[]}`: the first
  predicted delay, `map_index`, actual wait/attempt `action`, and null `pass` outcomes for skipped
  attempts. A wait is never reported as a pass.
- `convergence.iterations`; `convergence.work.{h_vectors_recorded,
  certificate_h_vectors,davidson_iterations_recorded,certificate_davidson_iterations}`;
  `convergence.work.cutoff_ladder.{coarse_h_vectors,coarse_start_h_vectors}`.
  Full H total is ordinary recorded + certificate + coarse + coarse start;
  the ladder's fine aliases must **not** be added a second time.
- Total failed-attempt H and Davidson iterations (sum only failed attempts),
  total attempts, final energy, verification status and density hash. The
  driver's manifest includes these reductions and external wall per run.
- `scf.switches`, `execution`, `run`, and complete JSON/history, not just medians.
  Check the MPI model is distributed and ranks agree. Return both A/A and A/B
  archives with failed/incomplete runs included.

Primary hypothesis for the final geometric-tail model on Si64: attempts **2 -> 1**, first failed-attempt H
**380 -> 0** on MPI8 / **395 -> 0** on one core, and failed-attempt Davidson
iterations decrease. Total H must decrease, not merely move into ordinary
maps. This is a falsifiable target, not a promised count: HALF omission can
alter the trajectory. The lag-only precursor failed to remove the first failure in the Mac Si64 check and is not this candidate. Si32 one core is a wrong-skip/regression control; do not
hide an added map or its cost. Inspect both cutoff levels because the policy
also applies to eligible coarse sessions.

For each model/case separately compute paired d_i = wall_B - wall_A, mean,
sample SE, median ratio and all raw deltas. Require negative mean with
`abs(mean)/SE >= 3`, `-mean_AB > max(abs(mean_AA), 3*SE_AA)` in seconds; an A/A resolved
bias or unstable load invalidates the block (rerun the entire block with a
reason, retain the first). No wins are inferred from comparing MPI and one
core. No rejection merely for a small resolved gain.

Every run must converge and pass the unchanged certificate. Compare to a
same-input strict reference with the same launcher/model: standard energy
error <=5e-9 Ha; the tighter-class generality run <=1e-10 Ha. A/B equality is
not expected for a trajectory change, and cannot replace the strict-reference
check. Withdraw the candidate on correctness failure, a resolved regression
on a control. Lack of failed-work/total-work reduction refutes the proposed
mechanism and must be reported; if wall nevertheless improves significantly,
identify the actual work change rather than reject a small valid speedup solely
for missing a counter target. If unresolved, report inconclusive and request more complete
pairs, not a default flip. Defaults remain unchanged pending review.
