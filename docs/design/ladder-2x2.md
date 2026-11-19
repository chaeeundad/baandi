# F-LADDER-2X2: captured state × first-fine tolerance

Mac pure-Rust release, one core; source base `d6be1d0cab87cb0654d6b6de2c45f8c9a1fe6cfb`.
This is a causal diagnostic, not a dev08 performance acceptance or a QE comparison.
No production source, switch, default, solver, mixer, threshold or certificate is changed.

## Decision: interaction and geometry dependence; keep the default

The si32 seed-601 state × seed-source cells have fine map counts
**10 / 11 / 11 / 13**, and fine-plus-verification H **1918 / 1987 / 2455 /
2233**, in order `(0.50, seed50)`, `(0.50, seed40)`, `(0.40, seed50)`,
`(0.40, seed40)`. The tighter seed removes two maps from the 0.40 state but
**adds 222 H columns (+9.94%)**. Loosening the 0.50 state adds one map and
69 H columns (+3.60%). The map difference-in-differences is one map; the H
interaction is −291 columns. A map-count recovery alone is not a work recovery. At this geometry, state40 also costs
more fine H than state50 at **either matched seed**, so a state-quality/history
penalty remains after the scalar exchange. The cross-geometry evidence below
makes interaction the overall decision, rather than a universal one-factor rule.

Seed 602 has seven coarse maps at both fractions and all four cells take ten
fine maps. Its 0.40 natural seed is already tighter than 0.50's. Seed 603 has
six coarse maps at both fractions; **loosening** the 0.40 state's seed removes
a map (12 → 11) and reduces fine-plus-verification H (2211 → 2064).
Consequently the experiment does not support a transferable "tighten the
0.40 handoff" rule. No production candidate or switch is introduced; the
conditional seven-system × three-seed rule-validation stage is not triggered.
This does not exclude a different future handoff policy with independent evidence.

Si64/601 is a held-out size, not another independent material family. Both
natural coarse stages take six maps; their seeds (4.296279381e-6 and
3.856824574e-6 Ha) differ by only 11.4%. The 0.40 state takes 12 fine maps
and 4061 fine-plus-verification H with either seed. Tightening the 0.50
state leaves 13 maps and raises H from 4621 to 4698. This control also
provides no evidence for spending extra inner-solve work to recover maps.

Si64/603 does have a discrete coarse-stop crossing (7 → 6 maps) and a
40.3-fold seed change (1.525041715e-7 → 6.152031800e-6 Ha). Tightening its
0.40 state **does not remove a fine map** (12 at either seed), and raises
fine-plus-verification H from 4154 to 4861 (+17.02%). Conversely, loosening
the 0.50 state adds a map (11 → 12) while reducing H (4348 → 4156). The
si32 cliff is therefore not a reliable argument for a universal tight-seed rule.
Si64/602 also worsens when its 0.50 state is tightened: 12 → 13 maps,
4556 → 4690 H. All four cells for every geometry are reported below.

The common-low-G density and occupied-subspace diagnostics below quantify
representation differences without assigning all excess work to them. In
si32/601 the common-low-G density error grows from 1.523828513e-4 to
4.094508017e-4 Bohr^-3/2 and the maximum occupied angle from 0.01995961 to
0.04223398 rad. The same size of representation degradation exists for seed
602 without the map cliff. Every arm starts with an empty fine mixer; its
subsequent nonlinear history and verification retries remain part of the measured
state/seed interaction. This is the plan's third branch, not a proof of a
representation-only or mixer-only mechanism.

## Intervention and accounting

The earlier [fraction study](ladder-fraction.md) exchanged seeds between six- and
seven-map states **at fraction 0.40**. It did not test the natural 0.50 and 0.40
handoffs. This experiment implements section 3 of the
[follow-up plan](../performance/plan-2026-11-19-followup.md) literally.

For each geometry (fractional displacements ±0.001, seeds 601–603), capture the
natural 0.50, 0.40 and 0.45 handoffs. A capture holds the mixed **next input
density** and the **last solved warm orbitals**, not a density reconstructed
from those orbitals. The first two immutable in-memory captures are each
continued with both observed scalar seeds. 0.45 remains a natural cross-check.
The seed uses the last coarse standard error in Ha through the existing ratchet;
it does not use the coarse density residual. Scalars remain full binary64 values.

Every fine continuation constructs the same full-cutoff problem and uses
`initialize_ladder_electronic` and `seed_standard_tolerance`, the exact import
and setter used by `finish_from`. No coarse mixer memory is transferred in any
arm. The fine Hamiltonian is identical within a captured-state pair. Between
states its density-dependent potential legitimately differs. No counterfactual
literal or new rule enters production. The extra driver is wholly `#[cfg(test)]`
and ignored by fast tests.

Each geometry also runs the unsplit production `run_ladder` at all three natural
fractions. Diagonal continuations must byte-match those controls in final energy,
eigenvalues, density, complex orbitals and recorded scientific history (effective tolerance,
redo, H work, Pulay diagnostics and charge errors; timings excluded). Repeated
arms must also match. A mismatch aborts the job, rather than weakening to 1e-12 Ha.

The fine timer includes import, session setup, all fine maps, every failed and
passed verification attempt, and result delivery. Fine problem construction is
separately charged. State construction includes coarse problem construction,
coarse session initialization/atomic start, all coarse maps, export and prolongation.
The reported **full price** is their sum. Because a capture is reused, this is a
sum of measured components, not an uninterrupted production wall; the natural
controls also retain their actual uninterrupted wall. Capture I/O, strict reference
and post-run diagnostics are outside both timers. A separate full-construction
pair driver below rebuilds every starting state inside each timed arm.

H counts are unweighted full-H columns, with atomic-start work included in Hc,
first-step redo included in Hf, and **all** verification attempts in Hv. The
per-map tables also give full-H block calls. Delta-potential cache updates are
not full-H applications. Coulomb residual means `E_H[rho_out-rho_in]` in Ha,
the actual standard estimator; it is distinct from the L2 density residual.
Charge is the absolute electron-number error `|Q-Ne|`, derived from the recorded
relative error and the known electron count. `maps.csv` retains the input charge
error too. A rejected attempt means a failed fixed-point verification; redo and
mixer reset/refit are separate columns. No attempts are silently dropped.

Two captured-state rounds reverse the entire execution order; every round has
an A/A repeat of the natural 0.50 cell. These small, loaded-Mac samples diagnose
work/trajectory, not a small speedup. Full-price dev08 pairs must independently
measure their noise floor, discard the first pair and retain raw ordered samples.

## Harness corrections

- `si64` already maps to `si64-gamma.toml` at this task's base commit. The new
  regression test resolves the real shipped fixture and checks that the alias
  and explicit name yield identical 64-atom inputs, without reading a UPF.
  It covers 8/16/32 too. No fictitious `si64.toml` is created.
- The legacy ignored `minimal-grid` arm previously constructed and discarded
  a full coarse problem before constructing another at its selected grid. It
  now enumerates the density basis directly, plans the same dimensions, and
  constructs exactly one coarse problem. The full timer includes this planning.
  This preserves its old **non-smooth** minimal-grid experiment; it does not
  revive or change the closed production smooth-grid policy.
- The legacy driver now accepts `fraction:stop:first-seed:<Ha>`, an explicit
  arbitrary seed override applied after `c.standard_next` is computed. This
  test-only option enables uninterrupted, full-construction crossed pairs.
  The natural-control equality assertion applies to natural arms, not to an
  intentionally counterfactual seed.

## Reproduction on the Mac

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export TMPDIR_990="$PWD/target/scratch" TMPDIR="$PWD/target/scratch"
mkdir -p "$TMPDIR_990"
cargo test --offline --release -p dft-engine --lib --no-run
python3 -I scripts/ladder_2x2_tests.py
python3 -I scripts/ladder_2x2_run.py \
  --binary target/release/deps/dft_engine-5a2965bdd204a9df \
  --out target/scratch/2x2 --cases si32,si64 --seeds 601,602,603 --pairs 2
python3 -I scripts/ladder_2x2_analyze.py target/scratch/2x2
scripts/check.sh fast dft-engine dft-numerics dft-cli
```

Use the executable printed by the local build; its suffix is not portable.
Raw immutable executable/provenance, fixture and UPF hashes, captured arrays,
strict/natural final arrays, per-map histories, exit records and paired times
are under `target/scratch/2x2`. Binary arrays are internal experiment captures:
little-endian f64, interleaved real/imaginary orbital entries in column-major
order, density in the existing ndarray traversal order. Shape/grid metadata
and complete principal-angle arrays are in each `state-*/state.json`; the
analysis also records the band count inferred from file size. They are not a
portable production checkpoint. `rows.json`, `maps.csv` and the Markdown tables
are the readable scientific record. Do not remove that scratch before the master
archives the evidence. The main result tables below are retained in this document.

## dev08 commands (master session only; not executed here)

Use a quiet reservation, the same checked-out sources, backend and bound core,
with `MKL_NUM_THREADS=1` and `OMP_NUM_THREADS=1`. Select the lib-test executable
printed by the release build. Replace the binary path below with that executable.
No MPI or QE speed comparison is inferred from this one-core experiment.

```sh
export TMPDIR_990="$PWD/target/scratch" TMPDIR="$PWD/target/scratch"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$TMPDIR_990"
cargo test --offline --release -p dft-engine --lib --no-run --message-format=json \
  > "$TMPDIR_990/2x2-build.jsonl"
ENGINE_TEST=$(python3 -c 'import json,sys; rows=[json.loads(s) for s in open(sys.argv[1])]; print(next(r["executable"] for r in rows if r.get("executable") and r.get("target",{}).get("name")=="dft_engine"))' "$TMPDIR_990/2x2-build.jsonl")
# Bind the Python parent too, so every child inherits the same allowed core.
taskset -c 0 python3 -I scripts/ladder_2x2_run.py \
  --binary "$ENGINE_TEST" --out target/scratch/dev08-capture \
  --cases si32,si64 --seeds 601,602,603 --pairs 2
python3 -I scripts/ladder_2x2_analyze.py target/scratch/dev08-capture
# Repeat separately for each case and seed: one geometry per job.
taskset -c 0 python3 -I scripts/ladder_2x2_pairs.py \
  --capture-root target/scratch/dev08-capture --binary "$ENGINE_TEST" \
  --case si32 --seed 601 --pairs 9 --out target/scratch/dev08-si32-601-full
```

`ladder_2x2_pairs.py` checks the capture/binary SHA-256 match, reads the two seeds
from that host's natural trajectories, and constructs explicit `first-seed` arms.
Within each round it runs an A/A floor plus both tolerance pairs; odd rounds
reverse order. It marks round 0 for exclusion and saves the full manifest and
per-pair raw wall/percentage differences. It never imports Mac literal seeds as
a dev08 policy. Repeat with `--case si64` and seeds 602/603 in separate jobs.

These harnesses enable existing diagnostic timers to count work. Use a separate
profiling reservation and report the instrumentation explicitly; do not promote
these diagnostic walls into an uninstrumented headline. Require a moving work
counter, an effect at least three times the A/A floor, separated uncertainty
intervals and stable raw spread before attributing a full-wall effect. No remote
job, reservation, system change or acceptance decision was made in this task.

## Complete measured results

Measured 19 November 2026 UTC (19–20 November KST). **72/72 continuations**,
**18/18 unsplit natural controls** and six separately tightened strict references
converged. All standard continuations passed the unchanged verification.
**48 diagonal/repeated-control byte comparisons passed**. The 30 unique cells
(24 factorial cells plus six 0.45 cross-checks) produced 824 timed map records
including repeats. The largest absolute strict-reference difference is
**3.0014746243978152e-10 Ha**, si32/602 natural 0.40, below 5e-9 Ha.
Only **17/30** unique cells are also within 1e-10 Ha; this is standard-class
evidence, not a claim that every standard result meets the strict class.

All recorded mixer resets and Coulomb refits are zero. The maximum output
charge error is **2.2737367544323206e-12 electrons**. There is one unique
first-map redo: si64/603, state 0.50 with seed40. Its injected seed is
6.152031800255166e-6 Ha, but the existing redo reduces the effective requested
tolerance to **5.426896967241076e-7 Ha**. Its first map costs **706 H columns
and 6 full-H block calls**, including both attempts. Both order repetitions
reproduce this redo. All other first maps use their injected seed without redo.

Binary SHA-256: `59d44797595dc0063bca9bdcb01796c050236d6728eeeadc8a23a018621c9382`. Measured engine source SHA-256:
`744da831c4b340ecf673efc15166e33e7d92bfaea45400002b1960d817434ed7`. Platform: `macOS-27.0.1-arm64-arm-64bit`.
The initial host load was [3.5771484375, 3.5107421875, 2.42578125]. No build ran concurrently
with the matrix. Shared-host activity remained uncontrolled. Raw arrays and
hashes are retained in the worktree scratch; no pseudopotential is committed.

### Natural seeds and state quality

The low-G domain is the common 0.40 **wavefunction** sphere for all three states
of a geometry. Density errors use sqrt(volume × mean(delta-density²))/Ne
(Bohr^-3/2). Angles are the full occupied-subspace principal angles against
the strict reference: acos(singular values of X_handoff^H X_strict), radians.
This is invariant to rotations within occupied degeneracies. The table gives
the maximum; complete spectra are retained in `state-*/state.json`.

| System/seed | Fraction | Coarse maps | Last L2 residual | Last E_H (Ha) | Injected seed (Ha) | Common low-G error | Max occupied angle (rad) |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| si32/601 | 0.50 | 7 | 2.6783792145e-05 | 1.2276468362e-04 | 9.5909909074236446e-08 | 1.5238285130e-04 | 1.9959605477e-02 |
| si32/601 | 0.40 | 6 | 9.9494722265e-05 | 3.3434762494e-03 | 2.6120908198666445e-06 | 4.0945080174e-04 | 4.2233977588e-02 |
| si32/601 | 0.45 | 7 | 2.5777905581e-05 | 1.1254779206e-04 | 8.7927962548893628e-08 | 2.6052269733e-04 | 2.6777330425e-02 |
| si32/602 | 0.50 | 7 | 2.7077188687e-05 | 1.3961132929e-04 | 1.0907135101095130e-07 | 1.5253132501e-04 | 1.9964820857e-02 |
| si32/602 | 0.40 | 7 | 2.4713350898e-05 | 1.1897740821e-04 | 9.2951100161367309e-08 | 4.1098258595e-04 | 4.1635753368e-02 |
| si32/602 | 0.45 | 7 | 2.5645794840e-05 | 1.2360376311e-04 | 9.6565439928614723e-08 | 2.6073529365e-04 | 2.6771049506e-02 |
| si32/603 | 0.50 | 6 | 8.7724406735e-05 | 2.5671228414e-03 | 2.0055647198471947e-06 | 1.5154516970e-04 | 2.0153694945e-02 |
| si32/603 | 0.40 | 6 | 7.9495783752e-05 | 2.1101045847e-03 | 1.6485192068268693e-06 | 4.0988835620e-04 | 4.1949794843e-02 |
| si32/603 | 0.45 | 6 | 8.4055362679e-05 | 2.3527802462e-03 | 1.8381095673349779e-06 | 2.5952126068e-04 | 2.6834803785e-02 |
| si64/601 | 0.50 | 6 | 8.7557413902e-05 | 1.0998475216e-02 | 4.2962793811650328e-06 | 1.0379383588e-04 | 2.0211598343e-02 |
| si64/601 | 0.40 | 6 | 8.3516322005e-05 | 9.8734709099e-03 | 3.8568245741957998e-06 | 2.9033013040e-04 | 4.2582520947e-02 |
| si64/601 | 0.45 | 6 | 8.6523149509e-05 | 1.0641621698e-02 | 4.1568834757117239e-06 | 1.8280627585e-04 | 2.6970717138e-02 |
| si64/602 | 0.50 | 6 | 8.2296729263e-05 | 9.8491978373e-03 | 3.8473429052012039e-06 | 1.0340723675e-04 | 2.0178739279e-02 |
| si64/602 | 0.40 | 6 | 7.6507612423e-05 | 8.3887391386e-03 | 3.2768512260015537e-06 | 2.9017671682e-04 | 4.2698754661e-02 |
| si64/602 | 0.45 | 6 | 8.0419852642e-05 | 9.3393222628e-03 | 3.6481727589237726e-06 | 1.8262227291e-04 | 2.6976658651e-02 |
| si64/603 | 0.50 | 7 | 2.4056609086e-05 | 3.9041067913e-04 | 1.5250417153458297e-07 | 1.0477060324e-04 | 1.9974361033e-02 |
| si64/603 | 0.40 | 6 | 9.8772486128e-05 | 1.5749201409e-02 | 6.1520318002551660e-06 | 2.9015951778e-04 | 4.3038862797e-02 |
| si64/603 | 0.45 | 6 | 9.9578594966e-05 | 1.5980531827e-02 | 6.2423952450781534e-06 | 1.8256853663e-04 | 2.7787022143e-02 |

### All four cells and 0.45 cross-checks

Time columns are medians of the observed component prices in seconds; they
are not acceptance timings. The A/A floor below shows why. Tighter/looser
labels refer to the initial injected seed, before any ordinary first-map redo.
Hf excludes Hv; the full price includes every verification attempt.

| System / seed | State / seed source | Fine maps | Hc / Hf / Hv | Total H | State cost (s) | Fine setup (s) | Fine incl. verify (s) | Full price (s) | Abs. error (Ha) |
| --- | --- | ---: | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| si32/601 | 0.50 / 0.50 (tight) | 10 | 1146 / 1717 / 201 | 3064 | 1.0040 | 0.0216 | 3.8296 | 4.8551 | 1.058e-10 |
| si32/601 | 0.50 / 0.40 (loose) | 11 | 1146 / 1596 / 391 | 3133 | 1.0040 | 0.0206 | 3.9677 | 4.9924 | 7.806e-11 |
| si32/601 | 0.40 / 0.50 (tight) | 11 | 1022 / 2066 / 389 | 3477 | 0.5979 | 0.0203 | 5.0804 | 5.6986 | 1.325e-10 |
| si32/601 | 0.40 / 0.40 (loose) | 13 | 1022 / 1872 / 361 | 3255 | 0.5979 | 0.0187 | 4.1223 | 4.7389 | 2.068e-11 |
| si32/601 | 0.45 / 0.45 (natural) | 10 | 1147 / 1782 / 207 | 3136 | 0.7833 | 0.0213 | 3.8090 | 4.6135 | 1.702e-10 |
| si32/602 | 0.50 / 0.50 (loose) | 10 | 1145 / 1735 / 197 | 3077 | 1.0560 | 0.0258 | 6.8928 | 7.9746 | 1.349e-10 |
| si32/602 | 0.50 / 0.40 (tight) | 10 | 1145 / 1753 / 198 | 3096 | 1.0560 | 0.0317 | 6.7063 | 7.7940 | 1.207e-10 |
| si32/602 | 0.40 / 0.50 (loose) | 10 | 1138 / 1842 / 191 | 3171 | 0.5886 | 0.0359 | 7.5528 | 8.1772 | 2.461e-10 |
| si32/602 | 0.40 / 0.40 (tight) | 10 | 1138 / 1843 / 194 | 3175 | 0.5886 | 0.0420 | 6.5348 | 7.1654 | 3.001e-10 |
| si32/602 | 0.45 / 0.45 (natural) | 10 | 1143 / 1779 / 201 | 3123 | 0.8614 | 0.1321 | 7.8429 | 8.8364 | 1.568e-10 |
| si32/603 | 0.50 / 0.50 (loose) | 11 | 1030 / 1809 / 378 | 3217 | 2.2090 | 0.0195 | 4.4990 | 6.7424 | 6.393e-11 |
| si32/603 | 0.50 / 0.40 (tight) | 11 | 1030 / 1810 / 383 | 3223 | 2.2090 | 0.0227 | 4.9080 | 7.1397 | 5.515e-11 |
| si32/603 | 0.40 / 0.50 (loose) | 11 | 1027 / 1672 / 392 | 3091 | 0.9367 | 0.0194 | 3.9989 | 4.9550 | 1.427e-10 |
| si32/603 | 0.40 / 0.40 (tight) | 12 | 1027 / 1825 / 386 | 3238 | 0.9367 | 0.0186 | 4.4438 | 5.3991 | 5.407e-11 |
| si32/603 | 0.45 / 0.45 (natural) | 12 | 1032 / 1832 / 356 | 3220 | 1.6919 | 0.0202 | 4.2263 | 5.9384 | 2.082e-11 |
| si64/601 | 0.50 / 0.50 (loose) | 13 | 2023 / 3871 / 750 | 6644 | 4.5752 | 0.0535 | 23.7518 | 28.3773 | 5.312e-11 |
| si64/601 | 0.50 / 0.40 (tight) | 13 | 2023 / 3933 / 765 | 6721 | 4.5752 | 0.0527 | 23.2122 | 27.8401 | 6.381e-11 |
| si64/601 | 0.40 / 0.50 (loose) | 12 | 2013 / 3334 / 727 | 6074 | 3.3159 | 0.0538 | 20.6975 | 24.0672 | 3.632e-11 |
| si64/601 | 0.40 / 0.40 (tight) | 12 | 2013 / 3334 / 727 | 6074 | 3.3159 | 0.0505 | 20.2126 | 23.5790 | 3.632e-11 |
| si64/601 | 0.45 / 0.45 (natural) | 13 | 2016 / 3944 / 703 | 6663 | 3.8948 | 0.0528 | 23.1262 | 27.0738 | 2.754e-11 |
| si64/602 | 0.50 / 0.50 (loose) | 12 | 2017 / 3778 / 778 | 6573 | 4.7027 | 0.0547 | 23.8749 | 28.6298 | 1.373e-10 |
| si64/602 | 0.50 / 0.40 (tight) | 13 | 2017 / 4021 / 669 | 6707 | 4.7027 | 0.0516 | 24.0779 | 28.8323 | 1.398e-11 |
| si64/602 | 0.40 / 0.50 (loose) | 12 | 2013 / 3355 / 795 | 6163 | 3.1188 | 0.0535 | 21.5874 | 24.7597 | 1.064e-10 |
| si64/602 | 0.40 / 0.40 (tight) | 12 | 2013 / 3381 / 787 | 6181 | 3.1188 | 0.0578 | 21.6680 | 24.8446 | 9.936e-11 |
| si64/602 | 0.45 / 0.45 (natural) | 11 | 2010 / 3558 / 775 | 6343 | 3.9551 | 0.0540 | 22.3073 | 26.3163 | 1.051e-10 |
| si64/603 | 0.50 / 0.50 (tight) | 11 | 2251 / 3574 / 774 | 6599 | 5.4409 | 0.0536 | 20.6673 | 26.1618 | 7.765e-11 |
| si64/603 | 0.50 / 0.40 (loose) | 12 | 2251 / 3387 / 769 | 6407 | 5.4409 | 0.0496 | 19.6630 | 25.1535 | 7.134e-11 |
| si64/603 | 0.40 / 0.50 (tight) | 12 | 2022 / 4158 / 703 | 6883 | 3.0075 | 0.0491 | 22.9884 | 26.0451 | 2.981e-11 |
| si64/603 | 0.40 / 0.40 (loose) | 12 | 2022 / 3365 / 789 | 6176 | 3.0075 | 0.0508 | 20.1354 | 23.1937 | 1.170e-10 |
| si64/603 | 0.45 / 0.45 (natural) | 13 | 2032 / 3931 / 769 | 6732 | 4.1069 | 0.0491 | 22.2684 | 26.4244 | 9.288e-11 |

### Actual uninterrupted natural-control times and references

| System/seed | Fraction | Natural full verified time (s) | Strict reference (Ha) | Absolute error (Ha) |
| --- | ---: | ---: | ---: | ---: |
| si32/601 | 0.50 | 4.930085 | -126.7034320500943 | 1.0579981335e-10 |
| si32/601 | 0.40 | 4.771001 | -126.7034320500943 | 2.0676793611e-11 |
| si32/601 | 0.45 | 4.737711 | -126.7034320500943 | 1.7020340692e-10 |
| si32/602 | 0.50 | 4.984593 | -126.703445499512 | 1.3490364381e-10 |
| si32/602 | 0.40 | 5.016451 | -126.703445499512 | 3.0014746244e-10 |
| si32/602 | 0.45 | 4.844969 | -126.703445499512 | 1.5680257093e-10 |
| si32/603 | 0.50 | 11.081906 | -126.7032917972904 | 6.3934635364e-11 |
| si32/603 | 0.40 | 15.459466 | -126.7032917972904 | 5.4072302191e-11 |
| si32/603 | 0.45 | 10.124806 | -126.7032917972904 | 2.0818902158e-11 |
| si64/601 | 0.50 | 26.579574 | -253.3781677212694 | 5.3120174925e-11 |
| si64/601 | 0.40 | 25.001598 | -253.3781677212694 | 3.6322944652e-11 |
| si64/601 | 0.45 | 27.221710 | -253.3781677212694 | 2.7540636438e-11 |
| si64/602 | 0.50 | 28.359713 | -253.3779038983341 | 1.3733369997e-10 |
| si64/602 | 0.40 | 25.534367 | -253.3779038983341 | 9.9362296169e-11 |
| si64/602 | 0.45 | 26.802242 | -253.3779038983341 | 1.0510348147e-10 |
| si64/603 | 0.50 | 26.859833 | -253.3779189682817 | 7.7648110164e-11 |
| si64/603 | 0.40 | 24.628667 | -253.3779189682817 | 1.1701217772e-10 |
| si64/603 | 0.45 | 27.623074 | -253.3779189682817 | 9.2882146419e-11 |

### Interaction and observed A/A floor

Difference-in-differences is cell(40,seed40) − cell(40,seed50) −
cell(50,seed40) + cell(50,seed50). The column is indexed by seed **source**,
not by a fixed assumption about which seed is tighter. H includes verification.
A/A values are ordered signed fine-time changes between identical 0.50 arms,
one per round. In particular si32/602 has a 39.8% A/A drift: no small wall
effect is resolved by this matrix, and times must not be pooled across geometries.

| System/seed | Map interaction | Fine+verify H interaction | A/A round 0 (%) | A/A round 1 (%) |
| --- | ---: | ---: | ---: | ---: |
| si32/601 | 1 | -291 | 1.1738 | -2.8234 |
| si32/602 | 0 | -15 | 14.3539 | 39.8226 |
| si32/603 | 1 | 141 | 1.8486 | 2.8144 |
| si64/601 | 0 | -77 | 7.0934 | -0.3773 |
| si64/602 | -1 | -116 | -5.3018 | -2.7457 |
| si64/603 | -1 | -515 | 0.7488 | -0.3731 |

### Per-map record for every unique cell

State/seed-source IDs: **0 = 0.50, 1 = 0.40, 2 = 0.45**. `H / blocks`
includes a first-step redo; `Hv / blocks` and rejected verifies are separate.
Every row below is the first occurrence of that cell; all 824 map observations
and map seconds, input charge error, density L2 residual, cumulative resets and
the full Pulay record are in `target/scratch/2x2/maps.csv` (Pulay is a JSON cell)
and the per-arm JSON histories. Zero verification H means no attempt, including
an ordinary scheduled wait; it does not mean an unrecorded successful certificate.

| System/seed | State/seed source | Map | Effective tol (Ha) | Redo | H / blocks | Hv / blocks | E_H residual (Ha) | Charge error (e) | Reset / refit | Rejected verifies |
| --- | --- | ---: | ---: | --- | --- | --- | ---: | ---: | --- | ---: |
| si32/601 | 0/0 | 1 | 9.59099091e-08 | False | 403/7 | 0/0 | 7.36224151e-04 | 1.990e-13 | False/False | 0 |
| si32/601 | 0/0 | 2 | 9.59099091e-08 | False | 151/3 | 0/0 | 3.17122228e-04 | 5.826e-13 | False/False | 0 |
| si32/601 | 0/0 | 3 | 9.59099091e-08 | False | 239/8 | 0/0 | 3.95033995e-03 | 4.832e-13 | False/False | 0 |
| si32/601 | 0/0 | 4 | 9.59099091e-08 | False | 210/5 | 0/0 | 7.55680994e-06 | 3.411e-13 | False/False | 0 |
| si32/601 | 0/0 | 5 | 5.90375777e-09 | False | 122/3 | 0/0 | 1.58709646e-06 | 0.000e+00 | False/False | 0 |
| si32/601 | 0/0 | 6 | 1.23991911e-09 | False | 112/3 | 0/0 | 1.74871499e-07 | 5.684e-14 | False/False | 0 |
| si32/601 | 0/0 | 7 | 1.36618358e-10 | False | 105/2 | 0/0 | 2.69935524e-08 | 1.563e-13 | False/False | 0 |
| si32/601 | 0/0 | 8 | 2.10887128e-11 | False | 149/3 | 0/0 | 5.61931386e-08 | 1.563e-13 | False/False | 0 |
| si32/601 | 0/0 | 9 | 2.10887128e-11 | False | 128/3 | 0/0 | 7.54677709e-09 | 7.105e-14 | False/False | 0 |
| si32/601 | 0/0 | 10 | 5.89591960e-12 | False | 98/2 | 201/4 | 1.07032003e-10 | 7.105e-14 | False/False | 0 |
| si32/601 | 0/1 | 1 | 2.61209082e-06 | False | 314/5 | 0/0 | 7.33110001e-04 | 2.416e-13 | False/False | 0 |
| si32/601 | 0/1 | 2 | 5.72742188e-07 | False | 81/2 | 0/0 | 2.84761382e-04 | 1.563e-13 | False/False | 0 |
| si32/601 | 0/1 | 3 | 2.22469830e-07 | False | 201/8 | 0/0 | 2.86802894e-03 | 3.411e-13 | False/False | 0 |
| si32/601 | 0/1 | 4 | 2.22469830e-07 | False | 172/5 | 0/0 | 1.30777265e-05 | 4.974e-13 | False/False | 0 |
| si32/601 | 0/1 | 5 | 1.02169738e-08 | False | 104/4 | 0/0 | 1.83071027e-06 | 1.137e-13 | False/False | 0 |
| si32/601 | 0/1 | 6 | 1.43024240e-09 | False | 102/3 | 0/0 | 4.11221209e-07 | 3.979e-13 | False/False | 0 |
| si32/601 | 0/1 | 7 | 3.21266569e-10 | False | 128/3 | 0/0 | 4.03551823e-07 | 5.684e-14 | False/False | 0 |
| si32/601 | 0/1 | 8 | 3.15274862e-10 | False | 122/3 | 0/0 | 1.70906387e-07 | 1.990e-13 | False/False | 0 |
| si32/601 | 0/1 | 9 | 1.33520615e-10 | False | 97/2 | 189/4 | 3.14851234e-09 | 3.268e-13 | False/False | 1 |
| si32/601 | 0/1 | 10 | 2.45977527e-12 | False | 184/5 | 0/0 | 1.26610252e-09 | 3.979e-13 | False/False | 0 |
| si32/601 | 0/1 | 11 | 9.89142594e-13 | False | 91/2 | 202/5 | 9.21801766e-11 | 1.705e-13 | False/False | 0 |
| si32/601 | 1/1 | 1 | 2.61209082e-06 | False | 339/6 | 0/0 | 3.34578879e-03 | 0.000e+00 | False/False | 0 |
| si32/601 | 1/1 | 2 | 2.61209082e-06 | False | 68/1 | 0/0 | 9.91610797e-04 | 1.421e-14 | False/False | 0 |
| si32/601 | 1/1 | 3 | 7.74695935e-07 | False | 196/9 | 0/0 | 4.04396741e-03 | 1.137e-13 | False/False | 0 |
| si32/601 | 1/1 | 4 | 7.74695935e-07 | False | 154/5 | 0/0 | 1.46807577e-04 | 5.684e-14 | False/False | 0 |
| si32/601 | 1/1 | 5 | 1.14693419e-07 | False | 97/2 | 0/0 | 1.01648994e-05 | 1.421e-13 | False/False | 0 |
| si32/601 | 1/1 | 6 | 7.94132762e-09 | False | 115/3 | 0/0 | 1.20781077e-06 | 7.105e-14 | False/False | 0 |
| si32/601 | 1/1 | 7 | 9.43602167e-10 | False | 122/3 | 0/0 | 4.43842203e-07 | 1.137e-13 | False/False | 0 |
| si32/601 | 1/1 | 8 | 3.46751721e-10 | False | 135/3 | 0/0 | 6.02224422e-07 | 1.421e-14 | False/False | 0 |
| si32/601 | 1/1 | 9 | 3.46751721e-10 | False | 110/3 | 0/0 | 3.35570037e-08 | 1.705e-13 | False/False | 0 |
| si32/601 | 1/1 | 10 | 2.62164092e-11 | False | 103/2 | 197/4 | 4.18839327e-09 | 1.563e-13 | False/False | 1 |
| si32/601 | 1/1 | 11 | 3.27218224e-12 | False | 185/4 | 0/0 | 3.53586371e-10 | 1.705e-13 | False/False | 0 |
| si32/601 | 1/1 | 12 | 2.76239353e-13 | False | 141/3 | 0/0 | 2.23829641e-10 | 1.990e-13 | False/False | 0 |
| si32/601 | 1/1 | 13 | 1.74866907e-13 | False | 107/3 | 164/4 | 2.73343001e-11 | 2.842e-13 | False/False | 0 |
| si32/601 | 1/0 | 1 | 9.59099091e-08 | False | 426/7 | 0/0 | 3.36709926e-03 | 0.000e+00 | False/False | 0 |
| si32/601 | 1/0 | 2 | 9.59099091e-08 | False | 223/10 | 0/0 | 1.21273543e-03 | 5.400e-13 | False/False | 0 |
| si32/601 | 1/0 | 3 | 9.59099091e-08 | False | 259/8 | 0/0 | 1.07776847e-02 | 8.527e-14 | False/False | 0 |
| si32/601 | 1/0 | 4 | 9.59099091e-08 | False | 238/7 | 0/0 | 7.97243208e-05 | 2.274e-13 | False/False | 0 |
| si32/601 | 1/0 | 5 | 6.22846256e-08 | False | 114/3 | 0/0 | 1.31922458e-05 | 2.842e-13 | False/False | 0 |
| si32/601 | 1/0 | 6 | 1.03064420e-08 | False | 108/3 | 0/0 | 5.04608611e-07 | 2.274e-13 | False/False | 0 |
| si32/601 | 1/0 | 7 | 3.94225477e-10 | False | 117/2 | 0/0 | 2.71618995e-08 | 2.842e-13 | False/False | 0 |
| si32/601 | 1/0 | 8 | 2.12202340e-11 | False | 156/3 | 0/0 | 7.97296119e-08 | 1.137e-13 | False/False | 0 |
| si32/601 | 1/0 | 9 | 2.12202340e-11 | False | 149/3 | 0/0 | 2.93923918e-08 | 1.990e-13 | False/False | 0 |
| si32/601 | 1/0 | 10 | 2.12202340e-11 | False | 100/2 | 191/4 | 4.19660632e-10 | 0.000e+00 | False/False | 1 |
| si32/601 | 1/0 | 11 | 3.27859869e-13 | False | 176/4 | 198/4 | 1.67051523e-10 | 5.684e-14 | False/False | 0 |
| si32/601 | 2/2 | 1 | 8.79279625e-08 | False | 412/7 | 0/0 | 1.59369207e-03 | 1.705e-13 | False/False | 0 |
| si32/601 | 2/2 | 2 | 8.79279625e-08 | False | 161/3 | 0/0 | 5.38000739e-04 | 1.705e-13 | False/False | 0 |
| si32/601 | 2/2 | 3 | 8.79279625e-08 | False | 245/8 | 0/0 | 4.02207769e-03 | 2.274e-13 | False/False | 0 |
| si32/601 | 2/2 | 4 | 8.79279625e-08 | False | 214/5 | 0/0 | 7.61077122e-06 | 1.563e-13 | False/False | 0 |
| si32/601 | 2/2 | 5 | 5.94591502e-09 | False | 138/3 | 0/0 | 5.55538399e-06 | 1.421e-13 | False/False | 0 |
| si32/601 | 2/2 | 6 | 4.34014374e-09 | False | 108/3 | 0/0 | 2.11724952e-07 | 1.563e-13 | False/False | 0 |
| si32/601 | 2/2 | 7 | 1.65410119e-10 | False | 119/3 | 0/0 | 2.51420403e-08 | 5.684e-14 | False/False | 0 |
| si32/601 | 2/2 | 8 | 1.96422190e-11 | False | 151/4 | 0/0 | 7.41260332e-08 | 2.842e-14 | False/False | 0 |
| si32/601 | 2/2 | 9 | 1.96422190e-11 | False | 142/3 | 0/0 | 1.77353015e-08 | 8.527e-14 | False/False | 0 |
| si32/601 | 2/2 | 10 | 1.38557043e-11 | False | 92/2 | 207/4 | 1.02933295e-10 | 5.684e-14 | False/False | 0 |
| si32/602 | 0/0 | 1 | 1.09071351e-07 | False | 403/7 | 0/0 | 7.48846641e-04 | 2.416e-13 | False/False | 0 |
| si32/602 | 0/0 | 2 | 1.09071351e-07 | False | 176/12 | 0/0 | 3.62057453e-04 | 1.563e-13 | False/False | 0 |
| si32/602 | 0/0 | 3 | 1.09071351e-07 | False | 230/5 | 0/0 | 5.21002937e-03 | 3.411e-13 | False/False | 0 |
| si32/602 | 0/0 | 4 | 1.09071351e-07 | False | 208/4 | 0/0 | 1.10355814e-05 | 8.527e-14 | False/False | 0 |
| si32/602 | 0/0 | 5 | 8.62154798e-09 | False | 118/3 | 0/0 | 1.52374209e-06 | 2.842e-14 | False/False | 0 |
| si32/602 | 0/0 | 6 | 1.19042351e-09 | False | 110/3 | 0/0 | 1.46438264e-07 | 4.121e-13 | False/False | 0 |
| si32/602 | 0/0 | 7 | 1.14404894e-10 | False | 107/2 | 0/0 | 2.28222991e-08 | 1.137e-13 | False/False | 0 |
| si32/602 | 0/0 | 8 | 1.78299211e-11 | False | 152/3 | 0/0 | 4.38773616e-08 | 5.684e-14 | False/False | 0 |
| si32/602 | 0/0 | 9 | 1.78299211e-11 | False | 135/3 | 0/0 | 1.06906117e-08 | 2.842e-14 | False/False | 0 |
| si32/602 | 0/0 | 10 | 8.35204042e-12 | False | 96/2 | 197/4 | 1.17247899e-10 | 2.842e-13 | False/False | 0 |
| si32/602 | 0/1 | 1 | 9.29511002e-08 | False | 405/7 | 0/0 | 7.48721722e-04 | 7.105e-14 | False/False | 0 |
| si32/602 | 0/1 | 2 | 9.29511002e-08 | False | 180/12 | 0/0 | 3.62597051e-04 | 3.411e-13 | False/False | 0 |
| si32/602 | 0/1 | 3 | 9.29511002e-08 | False | 235/5 | 0/0 | 5.23065897e-03 | 5.684e-14 | False/False | 0 |
| si32/602 | 0/1 | 4 | 9.29511002e-08 | False | 211/4 | 0/0 | 1.09287705e-05 | 1.563e-13 | False/False | 0 |
| si32/602 | 0/1 | 5 | 8.53810194e-09 | False | 118/3 | 0/0 | 1.52932740e-06 | 1.421e-13 | False/False | 0 |
| si32/602 | 0/1 | 6 | 1.19478703e-09 | False | 110/3 | 0/0 | 1.34820389e-07 | 1.990e-13 | False/False | 0 |
| si32/602 | 0/1 | 7 | 1.05328429e-10 | False | 113/2 | 0/0 | 1.89030339e-08 | 2.416e-13 | False/False | 0 |
| si32/602 | 0/1 | 8 | 1.47679952e-11 | False | 147/3 | 0/0 | 3.75227456e-08 | 2.842e-14 | False/False | 0 |
| si32/602 | 0/1 | 9 | 1.47679952e-11 | False | 138/3 | 0/0 | 1.25669397e-08 | 0.000e+00 | False/False | 0 |
| si32/602 | 0/1 | 10 | 9.81792160e-12 | False | 96/2 | 198/4 | 8.72910098e-11 | 1.137e-13 | False/False | 0 |
| si32/602 | 1/1 | 1 | 9.29511002e-08 | False | 418/7 | 0/0 | 3.37823987e-03 | 1.137e-13 | False/False | 0 |
| si32/602 | 1/1 | 2 | 9.29511002e-08 | False | 205/17 | 0/0 | 9.69381359e-04 | 2.274e-13 | False/False | 0 |
| si32/602 | 1/1 | 3 | 9.29511002e-08 | False | 237/5 | 0/0 | 3.50974979e-03 | 5.684e-14 | False/False | 0 |
| si32/602 | 1/1 | 4 | 9.29511002e-08 | False | 206/4 | 0/0 | 3.14164412e-05 | 8.527e-14 | False/False | 0 |
| si32/602 | 1/1 | 5 | 2.45440947e-08 | False | 125/3 | 0/0 | 1.97794805e-05 | 0.000e+00 | False/False | 0 |
| si32/602 | 1/1 | 6 | 1.54527192e-08 | False | 108/3 | 0/0 | 1.13615951e-07 | 2.274e-13 | False/False | 0 |
| si32/602 | 1/1 | 7 | 8.87624619e-11 | False | 163/4 | 0/0 | 2.49954157e-08 | 1.563e-13 | False/False | 0 |
| si32/602 | 1/1 | 8 | 1.95276685e-11 | False | 144/3 | 0/0 | 6.04311542e-08 | 4.121e-13 | False/False | 0 |
| si32/602 | 1/1 | 9 | 1.95276685e-11 | False | 140/3 | 0/0 | 2.50122512e-08 | 7.105e-14 | False/False | 0 |
| si32/602 | 1/1 | 10 | 1.95276685e-11 | False | 97/2 | 194/4 | 2.98439133e-10 | 8.527e-14 | False/False | 0 |
| si32/602 | 1/0 | 1 | 1.09071351e-07 | False | 417/7 | 0/0 | 3.37826595e-03 | 1.563e-13 | False/False | 0 |
| si32/602 | 1/0 | 2 | 1.09071351e-07 | False | 197/17 | 0/0 | 9.69377511e-04 | 3.979e-13 | False/False | 0 |
| si32/602 | 1/0 | 3 | 1.09071351e-07 | False | 235/5 | 0/0 | 3.50413606e-03 | 3.268e-13 | False/False | 0 |
| si32/602 | 1/0 | 4 | 1.09071351e-07 | False | 197/4 | 0/0 | 3.13402849e-05 | 1.421e-13 | False/False | 0 |
| si32/602 | 1/0 | 5 | 2.44845976e-08 | False | 126/3 | 0/0 | 1.98795240e-05 | 1.137e-13 | False/False | 0 |
| si32/602 | 1/0 | 6 | 1.55308781e-08 | False | 108/3 | 0/0 | 1.42436855e-07 | 5.684e-14 | False/False | 0 |
| si32/602 | 1/0 | 7 | 1.11278793e-10 | False | 155/3 | 0/0 | 2.06239048e-08 | 5.684e-14 | False/False | 0 |
| si32/602 | 1/0 | 8 | 1.61124256e-11 | False | 147/3 | 0/0 | 4.52334252e-08 | 5.684e-14 | False/False | 0 |
| si32/602 | 1/0 | 9 | 1.61124256e-11 | False | 148/3 | 0/0 | 3.41124433e-08 | 1.421e-14 | False/False | 0 |
| si32/602 | 1/0 | 10 | 1.61124256e-11 | False | 112/3 | 191/4 | 2.59543302e-10 | 8.527e-14 | False/False | 0 |
| si32/602 | 2/2 | 1 | 9.65654399e-08 | False | 412/7 | 0/0 | 1.60195336e-03 | 1.421e-13 | False/False | 0 |
| si32/602 | 2/2 | 2 | 9.65654399e-08 | False | 177/14 | 0/0 | 5.45616096e-04 | 1.705e-13 | False/False | 0 |
| si32/602 | 2/2 | 3 | 9.65654399e-08 | False | 231/5 | 0/0 | 4.17862781e-03 | 8.527e-14 | False/False | 0 |
| si32/602 | 2/2 | 4 | 9.65654399e-08 | False | 207/4 | 0/0 | 1.33175921e-05 | 0.000e+00 | False/False | 0 |
| si32/602 | 2/2 | 5 | 1.04043688e-08 | False | 129/3 | 0/0 | 6.65910111e-06 | 2.558e-13 | False/False | 0 |
| si32/602 | 2/2 | 6 | 5.20242274e-09 | False | 110/3 | 0/0 | 1.78300121e-07 | 8.527e-14 | False/False | 0 |
| si32/602 | 2/2 | 7 | 1.39296970e-10 | False | 121/3 | 0/0 | 2.26801098e-08 | 2.416e-13 | False/False | 0 |
| si32/602 | 2/2 | 8 | 1.77188358e-11 | False | 151/3 | 0/0 | 6.33606722e-08 | 2.274e-13 | False/False | 0 |
| si32/602 | 2/2 | 9 | 1.77188358e-11 | False | 145/3 | 0/0 | 2.06158995e-08 | 9.948e-14 | False/False | 0 |
| si32/602 | 2/2 | 10 | 1.61061715e-11 | False | 96/2 | 201/5 | 9.70167987e-11 | 3.268e-13 | False/False | 0 |
| si32/603 | 0/0 | 1 | 2.00556472e-06 | False | 327/7 | 0/0 | 7.94868150e-04 | 5.684e-14 | False/False | 0 |
| si32/603 | 0/0 | 2 | 6.20990742e-07 | False | 138/12 | 0/0 | 5.49488551e-04 | 0.000e+00 | False/False | 0 |
| si32/603 | 0/0 | 3 | 4.29287930e-07 | False | 206/6 | 0/0 | 9.89931345e-03 | 0.000e+00 | False/False | 0 |
| si32/603 | 0/0 | 4 | 4.29287930e-07 | False | 193/6 | 0/0 | 5.97230534e-05 | 2.558e-13 | False/False | 0 |
| si32/603 | 0/0 | 5 | 4.66586354e-08 | False | 106/3 | 0/0 | 1.76190859e-06 | 1.137e-13 | False/False | 0 |
| si32/603 | 0/0 | 6 | 1.37649108e-09 | False | 134/3 | 0/0 | 6.16790845e-08 | 5.684e-14 | False/False | 0 |
| si32/603 | 0/0 | 7 | 4.81867848e-11 | False | 147/3 | 0/0 | 3.82839808e-08 | 1.137e-13 | False/False | 0 |
| si32/603 | 0/0 | 8 | 2.99093600e-11 | False | 151/3 | 0/0 | 9.64014149e-08 | 1.847e-13 | False/False | 0 |
| si32/603 | 0/0 | 9 | 2.99093600e-11 | False | 139/3 | 0/0 | 2.37882693e-08 | 5.684e-14 | False/False | 0 |
| si32/603 | 0/0 | 10 | 1.85845854e-11 | False | 99/2 | 188/4 | 6.45837520e-10 | 1.847e-13 | False/False | 1 |
| si32/603 | 0/0 | 11 | 5.04560562e-13 | False | 169/4 | 190/5 | 8.31556651e-11 | 2.842e-13 | False/False | 0 |
| si32/603 | 0/1 | 1 | 1.64851921e-06 | False | 356/18 | 0/0 | 7.99866400e-04 | 0.000e+00 | False/False | 0 |
| si32/603 | 0/1 | 2 | 6.24895625e-07 | False | 110/3 | 0/0 | 5.57660660e-04 | 1.705e-13 | False/False | 0 |
| si32/603 | 0/1 | 3 | 4.35672391e-07 | False | 205/6 | 0/0 | 9.86108650e-03 | 8.527e-14 | False/False | 0 |
| si32/603 | 0/1 | 4 | 4.35672391e-07 | False | 192/6 | 0/0 | 8.45749251e-05 | 5.684e-14 | False/False | 0 |
| si32/603 | 0/1 | 5 | 6.60741603e-08 | False | 104/2 | 0/0 | 1.73775237e-06 | 1.421e-14 | False/False | 0 |
| si32/603 | 0/1 | 6 | 1.35761904e-09 | False | 148/4 | 0/0 | 6.61329664e-08 | 2.842e-14 | False/False | 0 |
| si32/603 | 0/1 | 7 | 5.16663800e-11 | False | 135/3 | 0/0 | 4.27898160e-08 | 1.705e-13 | False/False | 0 |
| si32/603 | 0/1 | 8 | 3.34295438e-11 | False | 153/3 | 0/0 | 1.24174195e-07 | 1.421e-13 | False/False | 0 |
| si32/603 | 0/1 | 9 | 3.34295438e-11 | False | 140/3 | 0/0 | 2.20031394e-08 | 5.684e-14 | False/False | 0 |
| si32/603 | 0/1 | 10 | 1.71899527e-11 | False | 100/2 | 184/4 | 4.78473811e-10 | 0.000e+00 | False/False | 1 |
| si32/603 | 0/1 | 11 | 3.73807665e-13 | False | 167/4 | 199/4 | 7.14932781e-11 | 1.705e-13 | False/False | 0 |
| si32/603 | 1/1 | 1 | 1.64851921e-06 | False | 367/8 | 0/0 | 3.39229940e-03 | 8.527e-14 | False/False | 0 |
| si32/603 | 1/1 | 2 | 1.64851921e-06 | False | 80/3 | 0/0 | 1.05803423e-03 | 1.421e-14 | False/False | 0 |
| si32/603 | 1/1 | 3 | 8.26589245e-07 | False | 196/10 | 0/0 | 5.88015287e-03 | 1.137e-13 | False/False | 0 |
| si32/603 | 1/1 | 4 | 8.26589245e-07 | False | 159/5 | 0/0 | 1.64624309e-04 | 5.684e-14 | False/False | 0 |
| si32/603 | 1/1 | 5 | 1.28612741e-07 | False | 102/3 | 0/0 | 1.46271034e-05 | 5.684e-14 | False/False | 0 |
| si32/603 | 1/1 | 6 | 1.14274245e-08 | False | 115/3 | 0/0 | 3.28601300e-07 | 7.105e-14 | False/False | 0 |
| si32/603 | 1/1 | 7 | 2.56719766e-10 | False | 152/3 | 0/0 | 3.27953820e-07 | 2.842e-13 | False/False | 0 |
| si32/603 | 1/1 | 8 | 2.56213922e-10 | False | 138/3 | 0/0 | 5.85778570e-07 | 1.705e-13 | False/False | 0 |
| si32/603 | 1/1 | 9 | 2.56213922e-10 | False | 119/3 | 0/0 | 3.32986062e-08 | 3.126e-13 | False/False | 0 |
| si32/603 | 1/1 | 10 | 2.60145361e-11 | False | 103/2 | 196/4 | 1.27230256e-09 | 1.563e-13 | False/False | 1 |
| si32/603 | 1/1 | 11 | 9.93986374e-13 | False | 179/4 | 0/0 | 3.37825330e-10 | 1.421e-14 | False/False | 0 |
| si32/603 | 1/1 | 12 | 2.63926039e-13 | False | 115/3 | 190/4 | 6.96634616e-11 | 1.563e-13 | False/False | 0 |
| si32/603 | 1/0 | 1 | 2.00556472e-06 | False | 355/7 | 0/0 | 3.39122964e-03 | 7.105e-14 | False/False | 0 |
| si32/603 | 1/0 | 2 | 2.00556472e-06 | False | 69/2 | 0/0 | 1.00634443e-03 | 1.705e-13 | False/False | 0 |
| si32/603 | 1/0 | 3 | 7.86206587e-07 | False | 197/10 | 0/0 | 4.25582093e-03 | 2.842e-13 | False/False | 0 |
| si32/603 | 1/0 | 4 | 7.86206587e-07 | False | 157/5 | 0/0 | 2.25464816e-04 | 4.121e-13 | False/False | 0 |
| si32/603 | 1/0 | 5 | 1.76144388e-07 | False | 99/2 | 0/0 | 1.44735463e-05 | 2.842e-14 | False/False | 0 |
| si32/603 | 1/0 | 6 | 1.13074580e-08 | False | 114/3 | 0/0 | 7.45946594e-07 | 5.684e-14 | False/False | 0 |
| si32/603 | 1/0 | 7 | 5.82770776e-10 | False | 139/3 | 0/0 | 4.71681260e-07 | 2.842e-14 | False/False | 0 |
| si32/603 | 1/0 | 8 | 3.68500984e-10 | False | 137/3 | 0/0 | 6.98098607e-07 | 3.837e-13 | False/False | 0 |
| si32/603 | 1/0 | 9 | 3.68500984e-10 | False | 114/3 | 0/0 | 3.22856269e-08 | 0.000e+00 | False/False | 0 |
| si32/603 | 1/0 | 10 | 2.52231460e-11 | False | 106/2 | 197/4 | 1.49427547e-09 | 2.842e-13 | False/False | 1 |
| si32/603 | 1/0 | 11 | 1.16740271e-12 | False | 185/4 | 195/4 | 1.83939756e-10 | 1.990e-13 | False/False | 0 |
| si32/603 | 2/2 | 1 | 1.83810957e-06 | False | 346/7 | 0/0 | 1.62623625e-03 | 1.137e-13 | False/False | 0 |
| si32/603 | 2/2 | 2 | 1.27049707e-06 | False | 99/9 | 0/0 | 6.86940000e-04 | 7.105e-14 | False/False | 0 |
| si32/603 | 2/2 | 3 | 5.36671875e-07 | False | 204/10 | 0/0 | 8.00831108e-03 | 5.684e-14 | False/False | 0 |
| si32/603 | 2/2 | 4 | 5.36671875e-07 | False | 178/6 | 0/0 | 6.61262864e-05 | 4.832e-13 | False/False | 0 |
| si32/603 | 2/2 | 5 | 5.16611613e-08 | False | 107/3 | 0/0 | 4.62834675e-06 | 1.847e-13 | False/False | 0 |
| si32/603 | 2/2 | 6 | 3.61589590e-09 | False | 115/4 | 0/0 | 5.15270707e-07 | 1.137e-13 | False/False | 0 |
| si32/603 | 2/2 | 7 | 4.02555240e-10 | False | 114/3 | 0/0 | 7.48324122e-08 | 3.268e-13 | False/False | 0 |
| si32/603 | 2/2 | 8 | 5.84628220e-11 | False | 148/3 | 0/0 | 1.00249951e-07 | 2.842e-13 | False/False | 0 |
| si32/603 | 2/2 | 9 | 5.84628220e-11 | False | 127/3 | 0/0 | 4.80533552e-08 | 2.842e-14 | False/False | 0 |
| si32/603 | 2/2 | 10 | 3.75416838e-11 | False | 102/3 | 191/4 | 1.41328698e-09 | 2.558e-13 | False/False | 1 |
| si32/603 | 2/2 | 11 | 1.10413045e-12 | False | 178/4 | 0/0 | 2.85468293e-10 | 1.421e-13 | False/False | 0 |
| si32/603 | 2/2 | 12 | 2.23022104e-13 | False | 114/3 | 165/4 | 2.69254275e-11 | 2.274e-13 | False/False | 0 |
| si64/601 | 0/0 | 1 | 4.29627938e-06 | False | 593/5 | 0/0 | 1.57219274e-03 | 4.547e-13 | False/False | 0 |
| si64/601 | 0/0 | 2 | 6.14137788e-07 | False | 232/5 | 0/0 | 1.17775296e-03 | 1.137e-13 | False/False | 0 |
| si64/601 | 0/0 | 3 | 4.60059749e-07 | False | 416/13 | 0/0 | 1.93784673e-02 | 5.116e-13 | False/False | 0 |
| si64/601 | 0/0 | 4 | 4.60059749e-07 | False | 380/5 | 0/0 | 2.31077406e-04 | 5.116e-13 | False/False | 0 |
| si64/601 | 0/0 | 5 | 9.02646118e-08 | False | 207/3 | 0/0 | 2.89141056e-06 | 6.821e-13 | False/False | 0 |
| si64/601 | 0/0 | 6 | 1.12945725e-09 | False | 302/4 | 0/0 | 4.93155951e-07 | 6.253e-13 | False/False | 0 |
| si64/601 | 0/0 | 7 | 1.92639043e-10 | False | 217/3 | 0/0 | 2.08882447e-07 | 1.023e-12 | False/False | 0 |
| si64/601 | 0/0 | 8 | 8.15947060e-11 | False | 289/3 | 0/0 | 2.18297878e-07 | 6.253e-13 | False/False | 0 |
| si64/601 | 0/0 | 9 | 8.15947060e-11 | False | 204/2 | 0/0 | 3.96945754e-08 | 1.251e-12 | False/False | 0 |
| si64/601 | 0/0 | 10 | 1.55056935e-11 | False | 210/3 | 0/0 | 7.79579010e-09 | 2.160e-12 | False/False | 0 |
| si64/601 | 0/0 | 11 | 3.04523051e-12 | False | 208/3 | 395/5 | 2.37027148e-09 | 5.684e-14 | False/False | 1 |
| si64/601 | 0/0 | 12 | 9.25887295e-13 | False | 378/5 | 0/0 | 3.18791262e-10 | 5.684e-14 | False/False | 0 |
| si64/601 | 0/0 | 13 | 1.24527837e-13 | False | 235/3 | 355/4 | 7.73816764e-11 | 5.684e-14 | False/False | 0 |
| si64/601 | 0/1 | 1 | 3.85682457e-06 | False | 600/5 | 0/0 | 1.57452910e-03 | 9.663e-13 | False/False | 0 |
| si64/601 | 0/1 | 2 | 6.15050429e-07 | False | 231/5 | 0/0 | 1.16297339e-03 | 2.842e-13 | False/False | 0 |
| si64/601 | 0/1 | 3 | 4.54286479e-07 | False | 416/13 | 0/0 | 1.88121671e-02 | 5.684e-14 | False/False | 0 |
| si64/601 | 0/1 | 4 | 4.54286479e-07 | False | 379/5 | 0/0 | 2.56495601e-04 | 6.821e-13 | False/False | 0 |
| si64/601 | 0/1 | 5 | 1.00193594e-07 | False | 204/2 | 0/0 | 2.52797436e-06 | 3.126e-13 | False/False | 0 |
| si64/601 | 0/1 | 6 | 9.87489986e-10 | False | 331/4 | 0/0 | 4.87694548e-07 | 5.116e-13 | False/False | 0 |
| si64/601 | 0/1 | 7 | 1.90505683e-10 | False | 217/3 | 0/0 | 2.03871645e-07 | 3.126e-13 | False/False | 0 |
| si64/601 | 0/1 | 8 | 7.96373615e-11 | False | 292/3 | 0/0 | 2.60532036e-07 | 7.390e-13 | False/False | 0 |
| si64/601 | 0/1 | 9 | 7.96373615e-11 | False | 229/3 | 0/0 | 7.05319222e-08 | 5.116e-13 | False/False | 0 |
| si64/601 | 0/1 | 10 | 2.75515321e-11 | False | 201/3 | 0/0 | 6.90084190e-09 | 1.705e-13 | False/False | 0 |
| si64/601 | 0/1 | 11 | 2.69564137e-12 | False | 212/2 | 395/5 | 3.14123737e-09 | 1.165e-12 | False/False | 1 |
| si64/601 | 0/1 | 12 | 1.22704585e-12 | False | 378/5 | 0/0 | 3.56319296e-10 | 2.274e-13 | False/False | 0 |
| si64/601 | 0/1 | 13 | 1.39187225e-13 | False | 243/3 | 370/4 | 9.34942840e-11 | 1.421e-13 | False/False | 0 |
| si64/601 | 1/1 | 1 | 3.85682457e-06 | False | 652/6 | 0/0 | 6.76003957e-03 | 5.116e-13 | False/False | 0 |
| si64/601 | 1/1 | 2 | 2.64064046e-06 | False | 132/1 | 0/0 | 2.02730191e-03 | 7.674e-13 | False/False | 0 |
| si64/601 | 1/1 | 3 | 7.91914810e-07 | False | 371/14 | 0/0 | 7.61330050e-03 | 1.705e-13 | False/False | 0 |
| si64/601 | 1/1 | 4 | 7.91914810e-07 | False | 297/5 | 0/0 | 4.13939733e-04 | 5.684e-14 | False/False | 0 |
| si64/601 | 1/1 | 5 | 1.61695208e-07 | False | 196/2 | 0/0 | 2.60275368e-05 | 5.116e-13 | False/False | 0 |
| si64/601 | 1/1 | 6 | 1.01670066e-08 | False | 225/3 | 0/0 | 2.59936646e-06 | 8.527e-13 | False/False | 0 |
| si64/601 | 1/1 | 7 | 1.01537752e-09 | False | 219/3 | 0/0 | 6.69297488e-07 | 2.842e-14 | False/False | 0 |
| si64/601 | 1/1 | 8 | 2.61444331e-10 | False | 248/3 | 0/0 | 3.32044228e-07 | 4.547e-13 | False/False | 0 |
| si64/601 | 1/1 | 9 | 1.29704777e-10 | False | 234/3 | 0/0 | 1.19989976e-07 | 1.137e-13 | False/False | 0 |
| si64/601 | 1/1 | 10 | 4.68710844e-11 | False | 216/3 | 0/0 | 1.26978043e-08 | 9.095e-13 | False/False | 0 |
| si64/601 | 1/1 | 11 | 4.96007980e-12 | False | 209/2 | 395/4 | 5.16713404e-10 | 1.137e-13 | False/False | 1 |
| si64/601 | 1/1 | 12 | 2.01841174e-13 | False | 335/4 | 332/4 | 5.27643585e-11 | 7.390e-13 | False/False | 0 |
| si64/601 | 1/0 | 1 | 4.29627938e-06 | False | 652/6 | 0/0 | 6.76003957e-03 | 5.116e-13 | False/False | 0 |
| si64/601 | 1/0 | 2 | 2.64064046e-06 | False | 132/1 | 0/0 | 2.02730191e-03 | 7.674e-13 | False/False | 0 |
| si64/601 | 1/0 | 3 | 7.91914810e-07 | False | 371/14 | 0/0 | 7.61330050e-03 | 1.705e-13 | False/False | 0 |
| si64/601 | 1/0 | 4 | 7.91914810e-07 | False | 297/5 | 0/0 | 4.13939733e-04 | 5.684e-14 | False/False | 0 |
| si64/601 | 1/0 | 5 | 1.61695208e-07 | False | 196/2 | 0/0 | 2.60275368e-05 | 5.116e-13 | False/False | 0 |
| si64/601 | 1/0 | 6 | 1.01670066e-08 | False | 225/3 | 0/0 | 2.59936646e-06 | 8.527e-13 | False/False | 0 |
| si64/601 | 1/0 | 7 | 1.01537752e-09 | False | 219/3 | 0/0 | 6.69297488e-07 | 2.842e-14 | False/False | 0 |
| si64/601 | 1/0 | 8 | 2.61444331e-10 | False | 248/3 | 0/0 | 3.32044228e-07 | 4.547e-13 | False/False | 0 |
| si64/601 | 1/0 | 9 | 1.29704777e-10 | False | 234/3 | 0/0 | 1.19989976e-07 | 1.137e-13 | False/False | 0 |
| si64/601 | 1/0 | 10 | 4.68710844e-11 | False | 216/3 | 0/0 | 1.26978043e-08 | 9.095e-13 | False/False | 0 |
| si64/601 | 1/0 | 11 | 4.96007980e-12 | False | 209/2 | 395/4 | 5.16713404e-10 | 1.137e-13 | False/False | 1 |
| si64/601 | 1/0 | 12 | 2.01841174e-13 | False | 335/4 | 332/4 | 5.27643585e-11 | 7.390e-13 | False/False | 0 |
| si64/601 | 2/2 | 1 | 4.15688348e-06 | False | 630/5 | 0/0 | 3.31625661e-03 | 3.411e-13 | False/False | 0 |
| si64/601 | 2/2 | 2 | 1.29541274e-06 | False | 194/5 | 0/0 | 1.48045970e-03 | 3.411e-13 | False/False | 0 |
| si64/601 | 2/2 | 3 | 5.78304572e-07 | False | 397/13 | 0/0 | 1.63482749e-02 | 5.684e-14 | False/False | 0 |
| si64/601 | 2/2 | 4 | 5.78304572e-07 | False | 352/5 | 0/0 | 2.83055096e-04 | 1.194e-12 | False/False | 0 |
| si64/601 | 2/2 | 5 | 1.10568397e-07 | False | 203/2 | 0/0 | 4.32213476e-06 | 4.547e-13 | False/False | 0 |
| si64/601 | 2/2 | 6 | 1.68833389e-09 | False | 296/4 | 0/0 | 3.91990832e-07 | 5.684e-14 | False/False | 0 |
| si64/601 | 2/2 | 7 | 1.53121419e-10 | False | 308/4 | 0/0 | 5.93964034e-07 | 2.842e-13 | False/False | 0 |
| si64/601 | 2/2 | 8 | 1.53121419e-10 | False | 279/3 | 0/0 | 3.50229946e-07 | 6.821e-13 | False/False | 0 |
| si64/601 | 2/2 | 9 | 1.36808573e-10 | False | 233/3 | 0/0 | 1.15498500e-07 | 1.251e-12 | False/False | 0 |
| si64/601 | 2/2 | 10 | 4.51166015e-11 | False | 203/3 | 393/5 | 3.36262572e-09 | 7.390e-13 | False/False | 1 |
| si64/601 | 2/2 | 11 | 1.31352567e-12 | False | 373/5 | 0/0 | 4.88549883e-10 | 5.684e-14 | False/False | 0 |
| si64/601 | 2/2 | 12 | 1.90839798e-13 | False | 258/3 | 0/0 | 2.79422876e-10 | 0.000e+00 | False/False | 0 |
| si64/601 | 2/2 | 13 | 1.09149561e-13 | False | 218/3 | 310/3 | 3.75603348e-11 | 3.695e-13 | False/False | 0 |
| si64/602 | 0/0 | 1 | 3.84734291e-06 | False | 597/5 | 0/0 | 1.60527513e-03 | 5.684e-13 | False/False | 0 |
| si64/602 | 0/0 | 2 | 6.27060598e-07 | False | 251/5 | 0/0 | 1.49790499e-03 | 1.137e-13 | False/False | 0 |
| si64/602 | 0/0 | 3 | 5.85119138e-07 | False | 420/13 | 0/0 | 2.76027183e-02 | 2.274e-12 | False/False | 0 |
| si64/602 | 0/0 | 4 | 5.85119138e-07 | False | 389/10 | 0/0 | 9.56255709e-05 | 2.274e-13 | False/False | 0 |
| si64/602 | 0/0 | 5 | 3.73537386e-08 | False | 232/3 | 0/0 | 3.49247611e-06 | 9.379e-13 | False/False | 0 |
| si64/602 | 0/0 | 6 | 1.36424848e-09 | False | 283/4 | 0/0 | 2.77290211e-07 | 1.705e-13 | False/False | 0 |
| si64/602 | 0/0 | 7 | 1.08316489e-10 | False | 250/3 | 0/0 | 1.71940606e-07 | 5.116e-13 | False/False | 0 |
| si64/602 | 0/0 | 8 | 6.71642993e-11 | False | 284/3 | 0/0 | 1.78639320e-07 | 1.705e-13 | False/False | 0 |
| si64/602 | 0/0 | 9 | 6.71642993e-11 | False | 239/3 | 0/0 | 4.99286704e-08 | 1.506e-12 | False/False | 0 |
| si64/602 | 0/0 | 10 | 1.95033869e-11 | False | 202/2 | 391/4 | 1.35025429e-09 | 1.023e-12 | False/False | 1 |
| si64/602 | 0/0 | 11 | 5.27443081e-13 | False | 376/4 | 0/0 | 3.06676243e-10 | 5.684e-13 | False/False | 0 |
| si64/602 | 0/0 | 12 | 1.19795407e-13 | False | 255/3 | 387/4 | 1.83711734e-10 | 2.842e-13 | False/False | 0 |
| si64/602 | 0/1 | 1 | 3.27685123e-06 | False | 616/5 | 0/0 | 1.60386007e-03 | 0.000e+00 | False/False | 0 |
| si64/602 | 0/1 | 2 | 6.26507840e-07 | False | 246/5 | 0/0 | 1.45175716e-03 | 0.000e+00 | False/False | 0 |
| si64/602 | 0/1 | 3 | 5.67092640e-07 | False | 417/12 | 0/0 | 2.62988438e-02 | 1.194e-12 | False/False | 0 |
| si64/602 | 0/1 | 4 | 5.67092640e-07 | False | 388/10 | 0/0 | 1.23859917e-04 | 1.194e-12 | False/False | 0 |
| si64/602 | 0/1 | 5 | 4.83827799e-08 | False | 218/3 | 0/0 | 2.12334130e-06 | 1.194e-12 | False/False | 0 |
| si64/602 | 0/1 | 6 | 8.29430196e-10 | False | 318/4 | 0/0 | 3.34737153e-07 | 4.547e-13 | False/False | 0 |
| si64/602 | 0/1 | 7 | 1.30756701e-10 | False | 224/3 | 0/0 | 1.66416633e-07 | 4.547e-13 | False/False | 0 |
| si64/602 | 0/1 | 8 | 6.50064973e-11 | False | 296/3 | 0/0 | 2.66280412e-07 | 3.411e-13 | False/False | 0 |
| si64/602 | 0/1 | 9 | 6.50064973e-11 | False | 257/3 | 0/0 | 6.03186424e-08 | 8.527e-13 | False/False | 0 |
| si64/602 | 0/1 | 10 | 2.35619697e-11 | False | 200/2 | 392/4 | 1.67096478e-09 | 9.095e-13 | False/False | 1 |
| si64/602 | 0/1 | 11 | 6.52720619e-13 | False | 379/4 | 0/0 | 6.95607611e-10 | 1.137e-13 | False/False | 0 |
| si64/602 | 0/1 | 12 | 2.71721723e-13 | False | 255/3 | 0/0 | 4.06140214e-10 | 1.421e-13 | False/False | 0 |
| si64/602 | 0/1 | 13 | 1.58648521e-13 | False | 207/2 | 277/3 | 1.83434885e-11 | 2.842e-14 | False/False | 0 |
| si64/602 | 1/1 | 1 | 3.27685123e-06 | False | 664/6 | 0/0 | 6.74893237e-03 | 3.411e-13 | False/False | 0 |
| si64/602 | 1/1 | 2 | 2.63630171e-06 | False | 137/4 | 0/0 | 2.11222513e-03 | 5.116e-13 | False/False | 0 |
| si64/602 | 1/1 | 3 | 8.25087940e-07 | False | 371/11 | 0/0 | 1.02825166e-02 | 9.948e-13 | False/False | 0 |
| si64/602 | 1/1 | 4 | 8.25087940e-07 | False | 312/5 | 0/0 | 2.81213064e-04 | 3.979e-13 | False/False | 0 |
| si64/602 | 1/1 | 5 | 1.09848853e-07 | False | 193/3 | 0/0 | 2.75242188e-05 | 3.695e-13 | False/False | 0 |
| si64/602 | 1/1 | 6 | 1.07516480e-08 | False | 218/3 | 0/0 | 4.21130321e-06 | 4.547e-13 | False/False | 0 |
| si64/602 | 1/1 | 7 | 1.64504032e-09 | False | 216/3 | 0/0 | 6.27527487e-07 | 5.969e-13 | False/False | 0 |
| si64/602 | 1/1 | 8 | 2.45127925e-10 | False | 282/3 | 0/0 | 6.26782636e-07 | 2.274e-13 | False/False | 0 |
| si64/602 | 1/1 | 9 | 2.44836967e-10 | False | 207/3 | 0/0 | 4.88890156e-08 | 1.137e-13 | False/False | 0 |
| si64/602 | 1/1 | 10 | 1.90972717e-11 | False | 217/3 | 0/0 | 1.00635842e-08 | 4.547e-13 | False/False | 0 |
| si64/602 | 1/1 | 11 | 3.93108758e-12 | False | 216/3 | 395/4 | 5.90963638e-10 | 2.274e-13 | False/False | 1 |
| si64/602 | 1/1 | 12 | 2.30845171e-13 | False | 348/4 | 392/4 | 1.30973594e-10 | 7.958e-13 | False/False | 0 |
| si64/602 | 1/0 | 1 | 3.84734291e-06 | False | 652/6 | 0/0 | 6.74622636e-03 | 2.274e-13 | False/False | 0 |
| si64/602 | 1/0 | 2 | 2.63524467e-06 | False | 137/4 | 0/0 | 2.11462114e-03 | 6.821e-13 | False/False | 0 |
| si64/602 | 1/0 | 3 | 8.26023883e-07 | False | 372/11 | 0/0 | 1.03479834e-02 | 6.537e-13 | False/False | 0 |
| si64/602 | 1/0 | 4 | 8.26023883e-07 | False | 312/5 | 0/0 | 2.83071440e-04 | 6.821e-13 | False/False | 0 |
| si64/602 | 1/0 | 5 | 1.10574781e-07 | False | 191/3 | 0/0 | 2.75341941e-05 | 7.390e-13 | False/False | 0 |
| si64/602 | 1/0 | 6 | 1.07555446e-08 | False | 219/3 | 0/0 | 4.21984366e-06 | 5.116e-13 | False/False | 0 |
| si64/602 | 1/0 | 7 | 1.64837643e-09 | False | 215/3 | 0/0 | 6.74192692e-07 | 2.842e-13 | False/False | 0 |
| si64/602 | 1/0 | 8 | 2.63356520e-10 | False | 275/3 | 0/0 | 6.54772334e-07 | 3.126e-13 | False/False | 0 |
| si64/602 | 1/0 | 9 | 2.55770443e-10 | False | 207/3 | 0/0 | 4.90521432e-08 | 5.684e-14 | False/False | 0 |
| si64/602 | 1/0 | 10 | 1.91609934e-11 | False | 216/3 | 0/0 | 9.82698051e-09 | 1.137e-13 | False/False | 0 |
| si64/602 | 1/0 | 11 | 3.83866426e-12 | False | 216/3 | 400/4 | 5.64275140e-10 | 0.000e+00 | False/False | 1 |
| si64/602 | 1/0 | 12 | 2.20419976e-13 | False | 343/4 | 395/4 | 1.35954191e-10 | 1.108e-12 | False/False | 0 |
| si64/602 | 2/2 | 1 | 3.64817276e-06 | False | 630/5 | 0/0 | 3.32228931e-03 | 0.000e+00 | False/False | 0 |
| si64/602 | 2/2 | 2 | 1.29776926e-06 | False | 205/5 | 0/0 | 1.72325804e-03 | 1.421e-13 | False/False | 0 |
| si64/602 | 2/2 | 3 | 6.73147672e-07 | False | 399/17 | 0/0 | 2.40276860e-02 | 1.023e-12 | False/False | 0 |
| si64/602 | 2/2 | 4 | 6.73147672e-07 | False | 372/10 | 0/0 | 1.21675503e-04 | 6.821e-13 | False/False | 0 |
| si64/602 | 2/2 | 5 | 4.75294934e-08 | False | 216/3 | 0/0 | 3.58565832e-06 | 4.832e-13 | False/False | 0 |
| si64/602 | 2/2 | 6 | 1.40064778e-09 | False | 292/4 | 0/0 | 2.32459782e-07 | 3.979e-13 | False/False | 0 |
| si64/602 | 2/2 | 7 | 9.08046024e-11 | False | 318/3 | 0/0 | 5.20630137e-07 | 6.537e-13 | False/False | 0 |
| si64/602 | 2/2 | 8 | 9.08046024e-11 | False | 306/3 | 0/0 | 4.23101530e-07 | 8.242e-13 | False/False | 0 |
| si64/602 | 2/2 | 9 | 9.08046024e-11 | False | 259/3 | 0/0 | 6.42759102e-08 | 3.411e-13 | False/False | 0 |
| si64/602 | 2/2 | 10 | 2.51077774e-11 | False | 205/2 | 380/4 | 1.60326008e-09 | 5.969e-13 | False/False | 1 |
| si64/602 | 2/2 | 11 | 6.26273470e-13 | False | 356/4 | 395/4 | 1.35878025e-10 | 1.421e-13 | False/False | 0 |
| si64/603 | 0/0 | 1 | 1.52504172e-07 | False | 775/6 | 0/0 | 1.39964425e-03 | 5.684e-14 | False/False | 0 |
| si64/603 | 0/0 | 2 | 1.52504172e-07 | False | 247/2 | 0/0 | 4.51136741e-04 | 2.842e-13 | False/False | 0 |
| si64/603 | 0/0 | 3 | 1.52504172e-07 | False | 370/14 | 0/0 | 1.89552572e-03 | 1.535e-12 | False/False | 0 |
| si64/603 | 0/0 | 4 | 1.52504172e-07 | False | 314/4 | 0/0 | 1.10726805e-04 | 1.990e-13 | False/False | 0 |
| si64/603 | 0/0 | 5 | 4.32526582e-08 | False | 199/2 | 0/0 | 1.36506408e-06 | 2.274e-13 | False/False | 0 |
| si64/603 | 0/0 | 6 | 5.33228155e-10 | False | 312/4 | 0/0 | 8.83484175e-07 | 4.832e-13 | False/False | 0 |
| si64/603 | 0/0 | 7 | 3.45111006e-10 | False | 217/3 | 0/0 | 2.24585986e-07 | 1.705e-13 | False/False | 0 |
| si64/603 | 0/0 | 8 | 8.77289009e-11 | False | 219/3 | 0/0 | 2.79455500e-08 | 2.274e-13 | False/False | 0 |
| si64/603 | 0/0 | 9 | 1.09162305e-11 | False | 283/3 | 0/0 | 6.32554500e-08 | 2.842e-13 | False/False | 0 |
| si64/603 | 0/0 | 10 | 1.09162305e-11 | False | 261/3 | 398/4 | 2.04295129e-09 | 2.274e-13 | False/False | 1 |
| si64/603 | 0/0 | 11 | 7.98027849e-13 | False | 377/5 | 376/4 | 9.44655560e-11 | 1.990e-13 | False/False | 0 |
| si64/603 | 0/1 | 1 | 5.42689697e-07 | True | 706/6 | 0/0 | 1.40185968e-03 | 9.095e-13 | False/False | 0 |
| si64/603 | 0/1 | 2 | 5.42689697e-07 | False | 132/1 | 0/0 | 4.11698881e-04 | 1.137e-13 | False/False | 0 |
| si64/603 | 0/1 | 3 | 1.60819876e-07 | False | 337/3 | 0/0 | 1.13261260e-03 | 1.137e-13 | False/False | 0 |
| si64/603 | 0/1 | 4 | 1.60819876e-07 | False | 302/10 | 0/0 | 1.16048660e-04 | 1.137e-13 | False/False | 0 |
| si64/603 | 0/1 | 5 | 4.53315078e-08 | False | 194/2 | 0/0 | 4.23240384e-06 | 3.979e-13 | False/False | 0 |
| si64/603 | 0/1 | 6 | 1.65328275e-09 | False | 219/3 | 0/0 | 8.77136640e-07 | 1.990e-13 | False/False | 0 |
| si64/603 | 0/1 | 7 | 3.42631500e-10 | False | 218/3 | 0/0 | 3.23708464e-07 | 1.137e-13 | False/False | 0 |
| si64/603 | 0/1 | 8 | 1.26448619e-10 | False | 231/3 | 0/0 | 6.78193412e-08 | 9.948e-13 | False/False | 0 |
| si64/603 | 0/1 | 9 | 2.64919302e-11 | False | 223/3 | 0/0 | 1.88870066e-08 | 5.684e-13 | False/False | 0 |
| si64/603 | 0/1 | 10 | 7.37773697e-12 | False | 246/3 | 0/0 | 6.58601951e-09 | 7.958e-13 | False/False | 0 |
| si64/603 | 0/1 | 11 | 2.57266387e-12 | False | 205/3 | 394/4 | 3.14004750e-10 | 1.990e-13 | False/False | 1 |
| si64/603 | 0/1 | 12 | 1.22658106e-13 | False | 374/4 | 375/4 | 8.87545698e-11 | 4.832e-13 | False/False | 0 |
| si64/603 | 1/1 | 1 | 6.15203180e-06 | False | 648/6 | 0/0 | 6.77642004e-03 | 7.674e-13 | False/False | 0 |
| si64/603 | 1/1 | 2 | 2.64703908e-06 | False | 141/4 | 0/0 | 2.04663078e-03 | 1.137e-13 | False/False | 0 |
| si64/603 | 1/1 | 3 | 7.99465150e-07 | False | 360/10 | 0/0 | 7.60328976e-03 | 5.684e-14 | False/False | 0 |
| si64/603 | 1/1 | 4 | 7.99465150e-07 | False | 313/10 | 0/0 | 4.93430981e-04 | 9.095e-13 | False/False | 0 |
| si64/603 | 1/1 | 5 | 1.92746477e-07 | False | 192/2 | 0/0 | 2.60496942e-05 | 1.421e-13 | False/False | 0 |
| si64/603 | 1/1 | 6 | 1.01756618e-08 | False | 221/3 | 0/0 | 2.75581741e-06 | 8.527e-13 | False/False | 0 |
| si64/603 | 1/1 | 7 | 1.07649117e-09 | False | 213/3 | 0/0 | 1.05568076e-06 | 6.253e-13 | False/False | 0 |
| si64/603 | 1/1 | 8 | 4.12375298e-10 | False | 253/3 | 0/0 | 6.13835135e-07 | 3.126e-13 | False/False | 0 |
| si64/603 | 1/1 | 9 | 2.39779350e-10 | False | 231/3 | 0/0 | 1.43152465e-07 | 1.137e-13 | False/False | 0 |
| si64/603 | 1/1 | 10 | 5.59189318e-11 | False | 207/3 | 0/0 | 7.24813217e-09 | 5.684e-13 | False/False | 0 |
| si64/603 | 1/1 | 11 | 2.83130163e-12 | False | 238/3 | 391/4 | 1.03749273e-09 | 3.979e-13 | False/False | 1 |
| si64/603 | 1/1 | 12 | 4.05270597e-13 | False | 348/4 | 398/5 | 1.48388330e-10 | 1.137e-13 | False/False | 0 |
| si64/603 | 1/0 | 1 | 1.52504172e-07 | False | 822/8 | 0/0 | 6.81551810e-03 | 3.979e-13 | False/False | 0 |
| si64/603 | 1/0 | 2 | 1.52504172e-07 | False | 383/11 | 0/0 | 2.38889156e-03 | 1.705e-13 | False/False | 0 |
| si64/603 | 1/0 | 3 | 1.52504172e-07 | False | 466/5 | 0/0 | 1.50418624e-02 | 1.080e-12 | False/False | 0 |
| si64/603 | 1/0 | 4 | 1.52504172e-07 | False | 435/5 | 0/0 | 5.40350078e-04 | 4.547e-13 | False/False | 0 |
| si64/603 | 1/0 | 5 | 1.52504172e-07 | False | 231/3 | 0/0 | 1.66827418e-05 | 8.527e-13 | False/False | 0 |
| si64/603 | 1/0 | 6 | 6.51669600e-09 | False | 219/3 | 0/0 | 1.38004314e-06 | 1.336e-12 | False/False | 0 |
| si64/603 | 1/0 | 7 | 5.39079353e-10 | False | 221/3 | 0/0 | 1.62277766e-07 | 3.979e-13 | False/False | 0 |
| si64/603 | 1/0 | 8 | 6.33897522e-11 | False | 274/3 | 0/0 | 7.89772724e-08 | 4.832e-13 | False/False | 0 |
| si64/603 | 1/0 | 9 | 3.08504970e-11 | False | 272/3 | 0/0 | 1.11057233e-07 | 1.990e-13 | False/False | 0 |
| si64/603 | 1/0 | 10 | 3.08504970e-11 | False | 229/3 | 391/4 | 3.52414614e-09 | 9.379e-13 | False/False | 1 |
| si64/603 | 1/0 | 11 | 1.37661958e-12 | False | 373/4 | 0/0 | 3.03595929e-10 | 3.411e-13 | False/False | 0 |
| si64/603 | 1/0 | 12 | 1.18592160e-13 | False | 233/3 | 312/3 | 3.44010996e-11 | 4.547e-13 | False/False | 0 |
| si64/603 | 2/2 | 1 | 6.24239525e-06 | False | 628/6 | 0/0 | 3.35285158e-03 | 4.832e-13 | False/False | 0 |
| si64/603 | 2/2 | 2 | 1.30970765e-06 | False | 208/5 | 0/0 | 1.58513093e-03 | 3.126e-13 | False/False | 0 |
| si64/603 | 2/2 | 3 | 6.19191768e-07 | False | 400/14 | 0/0 | 1.80759604e-02 | 1.307e-12 | False/False | 0 |
| si64/603 | 2/2 | 4 | 6.19191768e-07 | False | 354/5 | 0/0 | 4.03248026e-04 | 7.390e-13 | False/False | 0 |
| si64/603 | 2/2 | 5 | 1.57518760e-07 | False | 203/2 | 0/0 | 3.24659745e-06 | 2.842e-14 | False/False | 0 |
| si64/603 | 2/2 | 6 | 1.26820213e-09 | False | 329/4 | 0/0 | 7.08703365e-07 | 1.080e-12 | False/False | 0 |
| si64/603 | 2/2 | 7 | 2.76837252e-10 | False | 261/3 | 0/0 | 5.80289629e-07 | 4.547e-13 | False/False | 0 |
| si64/603 | 2/2 | 8 | 2.26675636e-10 | False | 268/3 | 0/0 | 6.72015239e-07 | 1.137e-13 | False/False | 0 |
| si64/603 | 2/2 | 9 | 2.26675636e-10 | False | 238/3 | 0/0 | 1.54047683e-07 | 7.390e-13 | False/False | 0 |
| si64/603 | 2/2 | 10 | 6.01748760e-11 | False | 201/2 | 385/4 | 3.63397838e-09 | 0.000e+00 | False/False | 1 |
| si64/603 | 2/2 | 11 | 1.41952281e-12 | False | 370/4 | 0/0 | 2.46898972e-09 | 2.842e-13 | False/False | 0 |
| si64/603 | 2/2 | 12 | 9.64449110e-13 | False | 265/3 | 0/0 | 1.78481688e-09 | 1.364e-12 | False/False | 0 |
| si64/603 | 2/2 | 13 | 6.97194094e-13 | False | 206/2 | 384/4 | 1.15983404e-10 | 1.421e-13 | False/False | 0 |


## Additional harness validation

The pristine pre-edit release executable is retained as
`target/scratch/baseline-engine`. Three si32/601 natural trajectories at
0.50/0.40/0.45 match it in energy bits, handoff density/orbital bytes and every
shared coarse/fine history field. The only additional history fields in the
new writer are input/output charge errors. Separately, the entire production
prefix of `cutoff_ladder.rs` before `#[cfg(test)] mod handoff_experiment` is
byte-identical to HEAD; no production numerical path was edited.

The minimal-grid regression uses si16/601, with default coarse grid
25×25×54 and the old experimental non-smooth minimal grid 25×25×51. Before
and after the setup fix, **each arm** preserves all four coarse/handoff binary
arrays, both histories, energy/reference values and H counters. Both report
646 ordinary coarse H, 916 fine H and 197 verification H. The test does not
assert that the default and minimal grids are mutually identical. Constructor
count falls from two to one only in the minimal arm; single samples of elapsed
time do not establish a speedup. Raw evidence: `target/scratch/minimal-before`
and `target/scratch/minimal-after`.

The full-construction pair harness also ran locally on si32/601: **14/14 arms**,
including two order-alternated rounds, two A/A runs per round, the natural 0.50
control and the 0.45 cross-check. Every rebuilt handoff matches the original
capture byte for byte, and fixture/UPF hashes match. After discarding the first
round, the observed full verified times were:

| Starting state | seed50 (s) | seed40 (s) |
| --- | ---: | ---: |
| 0.50 | 4.391110875 | 4.624032375 |
| 0.40 | 4.779347375 | 4.470593083 |

The retained A/A difference is +1.83597%. This is a **functional smoke test**
with only one retained pair, not a resolved timing effect or dev08 acceptance.
It confirms that the full-construction command path accepts the literal seeds,
charges all setup and reproduces the captured states. Raw ordered measurements,
A/A samples and warmup flags are in `target/scratch/full-pairs-smoke/pairs.json`.

`python3 -I scripts/ladder_2x2_tests.py`: **3 passed**, exit 0 (real alias
resolution/counts and pair ordering/A/A separation/seed round-trip precision).
`python3 -I scripts/ladder_2x2_analyze.py target/scratch/2x2`: exit 0; validates
all six geometries, same fine basis/grid, independent ratchet arithmetic,
work/history reconciliation, charge, verification and byte-identical controls.
`python3 target/scratch/extra_checks.py`: exit 0; reproduces the pristine-natural,
minimal-grid and full-construction comparisons described above. Exact child
commands/exit codes are saved in `target/scratch/extra-command-results.json`.

The independent README check ran 9 tests: **8 passed, 1 failed**, exit 1.
`test_structural_scope_and_sources` expects the obsolete phrase
`standard profile is unsupported`. This pre-existing assertion is outside
OWNED_FILES; README, its generator and stored reference artifacts are unchanged.
The task adds no adopted production result to the README. Dashboard/registry
publication and evidence archival remain with the master session; the worktree's
copied registry is stale, while the master's registry identifies this run.

### Final code and required gate

The first required gate reported a Clippy `useless_vec` in the ignored driver.
The fixed six-arm execution list now uses an array; its construction is outside
all measured SCF intervals. After that lint-only edit, the final release build
and a fresh si32/601 six-arm confirmation pass. **49 captured/reference/control/
continuation binary and history files are byte-identical** to the main matrix;
work and all four cells reproduce. `target/scratch/final-identity.json` records
this check. The main 72-arm evidence keeps its original binary/source hashes
above, rather than relabelling it as a different build.

Final binary SHA-256:
`4645466c0b971c02b6cbc420935e301465a0a806006e7a5d9eb4ae8a0f7b4c8f`.
Final engine-source SHA-256:
`b81d23a978f71d8685f08587959d15ee75c70b1854502cc5ace5fa002767ea3f`.
The entire production prefix remains byte-identical to HEAD.

`CHECK_SCCACHE=0 scripts/check.sh fast dft-engine dft-numerics dft-cli`:
**exit 1**, **build/fmt/clippy PASS**, **470 passed / 3 failed / 135 ignored**
across 51 binaries. The final gate took 17 s; the first cold gate took 152 s
and also exposed the now-fixed Clippy issue. All three remaining failures are
the existing in-repository TMPDIR pseudopotential-cache restriction:

- `pseudo::cache_resolution_and_no_network_side_effects`
- `pseudo::cached_fetch_checks_hash_and_provenance`
- `pseudo::tests::fetch_boundary_verifies_before_publication_and_never_overwrites`

The sandbox also denies the gate's optional `ps` load-inventory probe; build,
test, fmt and clippy still run and report their real results. No cache test,
scientific threshold or ownership boundary was bypassed. Logs are retained in
`target/scratch/check.log` and `target/scratch/check-final.log`. The README
failure above is separate from these three CLI failures. No dev08 gate is
claimed. Changes are uncommitted for master review.
