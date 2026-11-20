# LANE-GLUE Phase 2 closing recommendation

**Recommend closing this lane for the current work order.** No candidate has
a substantiated expected saving >=0.2 s on Si64 MPI8 or >=0.8 s on Si64 one
core. [diagnosis.md](diagnosis.md) records the relayed e947a35 profile, source
inspection, work/traffic models, five candidates and reasons not to implement.
No code/default change, switch, A/B request, performance claim or remote job.

Key distinctions:

- Density is 0.7331 s MPI8 / 6.3032 s one core. Gamma pairing and MPI fused
  accumulation are already adopted. The owned path has no final cross-rank
  density sum to remove. A serial density fusion could save logical traffic,
  but the supplied data does not price those passes or establish >=0.8 s.
- `m2.xc` is an inclusive potential/audit label, not the pure XC kernel.
  37.4% rank asymmetry is not evidence of equivalent removable work.
- Energy terms and FFTs already share substantial work. Omitting verification
  density would invalidate verification. Earlier ENERGY_IDENTITY was slower.
- The original owned implementation files cannot change serial orbital
  density/FFT or session audit scheduling. No out-of-scope edit was made.

Reopen only with a new premise: the existing raw rank-local GLUE counters and
FFT/work records identify removable work above the requested budget, or a
separately owned density/FFT task establishes a qualifying serial mechanism.
Useful evidence is `glue_lane.grid.{forward,inverse,xc,gather,reduce}` times and
work, actual per-phase grids/occupied pairs, and isolated serial contribution
clear/store/fold costs on dev08. These are reopening criteria, not a pending
request to the user or a request to run further jobs now. The precise XC
asymmetry cause and seconds above a hardware FFT/bandwidth floor remain unknown.

If reopened for a concrete candidate, pre-register its counter, require an
A/A noise floor followed by crossed order-alternated A/B pairs on Si64/Si32
MPI8 and Si64/Si32 one core, identical numerical inputs/binding within each
model, `abs(mean saving)/SE >= 3`, and the requested <=1e-11 Ha equivalence
and rank-count independence. No current candidate is authorized by this note.

The implementation-free closing branch does not require a dummy
`dev08-request.md`: that filename is reserved for an actual implementable A/B.
Changes stay uncommitted for review. Stop after reporting the verification
command results; do not start another optimization or measurement phase.

## Verification performed

- `python3 -I docs/performance/lanes/GLUE/phase2_floors.py`: exit 0;
  regenerates analytic JSON without running DFT. Counts/unit assertions:
  960 paired transforms, 6.3032/15 normalization and break-even units pass.
- `TMPDIR_990="$PWD/target/scratch" CHECK_SCCACHE=0 scripts/check.sh fast dft-engine`:
  exit **0**, build/test/fmt/clippy PASS; 33 binaries, **260 passed, 0 failed,
  74 ignored**. `target/scratch/glue-phase2-check.log` preserves the output.
  The sandbox-denied `ps` load probe is ancillary; all four gates executed.
- `git diff --check`: exit 0. No Rust source changed in this phase; no new
  equivalence or speed claim is made. Phase 0's tests remain historical.
