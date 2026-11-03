# ADR 0013: CLI k-point input and reporting

- Date: 2026-11-03
- Task: P7A2-KPOINTS-CLI
- Status: implemented; proposed for reviewer acceptance

## Input and compatibility

Under `[basis]`, choose one of:

- Legacy `kpoints = [[0,0,0]]`: preserve the legacy engine constructor and
  results.json structure, including occupied-only band compatibility.
- `kgrid = [n1,n2,n3]`, positive integers, with optional
  `kshift = [s1,s2,s3]` (default zeros; each exactly 0 or 0.5), and
  `ktime_reversal = true|false` (default true).
- `kpoints = [[kx,ky,kz,weight], ...]`, finite fractional reciprocal coordinates
  and positive finite weights summing to one within the geometry validator's
  1e-14 allowance. No normalization. `ktime_reversal` defaults false here;
  explicit true requests the geometry library's reduction.

A grid may accompany the exact legacy Gamma placeholder, which it supersedes.
Any other simultaneous list/grid is rejected. A shift without a grid is rejected.
Empty lists, three-column non-Gamma lists, mixed row lengths, nonfinite values,
zero grid dimensions and other shifts fail explicitly. A legacy Gamma list with
an explicit `ktime_reversal` option opts into the new API and output block.

Generation and reduction delegate to dft-geometry (ADR 0011):
`k_i=(m_i+s_i)/n_i`, wrapped to `[-1/2,1/2)`. There is no space-group reduction.
New forms use the engine's insulator policy with minimum gap 1e-7 Ha; the engine
also requires the gap to exceed twice the eigen tolerance. Gap failures retain
the engine explanation and exit 4. No occupations or tolerances are weakened.

## Output

Legacy inputs omit the new block entirely. New forms add `kpoints` with `grid`
(integer triple or null), `shift` (triple or null), `explicit` (Boolean), `count`,
`time_reversal` (whether reduction was requested) and `points`. Each point has
`fractional`, `weight`, `eigenvalues: {ha, ev}`, and `occupations`. Coordinates
are the engine's canonical coordinates; eV conversion uses dft-units.
Existing eigenvalue/occupation fields refer to the first k point. The basis
plane-wave count likewise describes the first k point; memory remains an
estimate for one orbital block, not total multi-k storage. Density formats stay
unchanged. Schema version remains 1 for the additive block.

Console eigenvalue tables include all points through count 10, otherwise the
first and last three with a results.json pointer. Labels k1, k2, ... follow JSON
order. The minimum direct gap includes its k label; the indirect gap identifies
both HOMO and LUMO k labels. An unconverged gap is diagnostic only.

## QE interchange

The generator consumes results.json, whose grid and shift retain the TOML
request. `K_POINTS automatic` uses the full grid dimensions and flags `2*s_i`:
0 means Gamma-centered unshifted, 1 means a half-grid-step shift. This is the
same formula as ADR 0011, not the original centered Monkhorst-Pack formula.
`nosym=.true.` prevents space-group reduction. With time reversal requested,
`noinv` is omitted: QE may merge k/-k, with multiplicity weights. With reduction
disabled, `noinv=.true.` retains the full grid. Representatives and ordering may
differ, so compare modulo reciprocal translations and, when enabled, sign.

Explicit inputs emit `K_POINTS crystal` using the resolved CLI list and weights.
Both flags are true to preserve that list, which already includes any explicitly
requested reduction. Legacy Gamma output is byte-for-byte unchanged, including
its automatic 1x1x1 card and symmetry flags.

The `.meta.json` sidecar includes `expected_kpoints`. The importable
`make_qe_input.compare_kpoints` checks coordinates and weights independent of
ordering; callers must convert QE Cartesian coordinates to fractional reciprocal
coordinates and spin weights summing to 2 into integration weights summing to 1.
The existing Gamma-only `verify_qe.py` protocol is not extended by this ticket
(it is outside owned paths); multi-k callers use this comparison helper.
No QE execution or new reference numbers are claimed.

## Verification

Run `scripts/check.sh` and `python3 -I -m unittest discover -s scripts/tests`.
Multiple tight SCF comparisons are ignored in debug to preserve the fast gate:
`cargo test --offline --release -p dft-cli --test run kpoint_scf -- --ignored --nocapture`.
All pseudopotentials in tests are synthetic and temporary. No dependencies added.

Measured synthetic CLI evidence (not QE references): full 2x1x1 half-shifted
x grid, equivalent explicit list, time-reversal-reduced grid and explicitly
reduced list all give -39.680996009155528 Ha (maximum energy difference 0 Ha,
required bound 1e-10 Ha). Legacy Gamma versus 1x1x1 grid, weighted Gamma list,
and grid with legacy Gamma placeholder gives bitwise equal energies (0 Ha
difference). The two-species gap-less 3x1x1 case exits 4 and reports sampled
indirect gap -0.06694262433472553 Ha, required gap 1e-7 Ha. The unchanged
large Gamma fixture gives -4.6830234053637234 Ha with direct-engine difference
0 Ha and charge error 1.243e-14 electrons. A 16-point budget-exhausted run keeps
all 16 JSON records and prints exactly six console tables.
