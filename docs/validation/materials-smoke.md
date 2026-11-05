# P13D real NC bulk smoke

Scope: PseudoDojo NC-SR v0.4 PBE standard Si and Mg/O, fixed full occupations,
unpolarized PBE through xcx. No QE calculation, accuracy claim, or cutoff
convergence claim is made. The external files are identified by the existing
manifest and individual SHA-256 comments in `tests/reference/materials/*.toml`.
No UPF is copied or committed.

## Reproduction

```sh
cargo build --offline --release -p dft-cli
DFT_NC_DIR=/absolute/path/to/pseudodojo/ex python3 -I scripts/run_materials_smoke.py
DFT_NC_DIR=/absolute/path/to/pseudodojo/ex cargo test --offline --release -p dft-cli --test materials -- --ignored --nocapture
```

The Python runner uses only the standard library, including on Python 3.9.
It resolves relative UPF names into temporary input files, verifies both manifest
and individual UPF hashes, and stores numerical summaries (not pseudopotential
contents). `--work-dir /tmp/materials` retains full CLI logs and output.
`--implicit-spin` removes only the explicit unpolarized spin block for a paired
reduction check; use `--output` to keep its summary separately.

The fcc primitive cells contain two atoms. Silicon uses a=10.26 Bohr and
fractional positions (0,0,0), (1/4,1/4,1/4); MgO uses a=7.95 Bohr and
(0,0,0), (1/2,1/2,1/2). Shifted meshes are 4x4x4 and 3x3x3 respectively,
with shift (1/2,1/2,1/2) and existing time reversal reduction. Inputs select
Kerker-Pulay beta=0.3, q0=1 Bohr^-1 and four extra bands, four threads,
energy tolerance 1e-8 Ha and density tolerance 1e-7. No mixing/solver default
was changed. The CLI preflight previously rejected unpolarized Kerker-Pulay
although its option resolver and engine support it; that rejection is fixed.

## XC density policy

For the unpolarized run path, an explicit `[xc] density_floor` is positive and
finite in electrons/Bohr^3. NLCC inputs default to 1e-10; non-NLCC inputs
without the key retain their original path. The engine's explicit
`ScfProblem::set_density_floor` permits the same opt-in at the Rust boundary.

Define t=n_valence+n_core and u=max(t,floor). Only XC evaluates u; Hartree,
orbital density, stored density, and immutable core are unchanged. PBE takes
spectral gradients of u, evaluates the library functional and assembles its
usual divergence, then multiplies the completed potential by du/dt (zero
below the floor). LDA uses the same chain rule. This is the derivative of the
discrete clamped functional away from the floor and the library's screening
boundary. Existing functional screening remains in place. Equality with the
floor is a nondifferentiable boundary, not a finite-difference oracle point.

The reported clipped negative charge is integral max(-t,0) dV, in electrons;
it does not count positive density raised to the floor. Input and output values
are recorded per iteration in `energies.xc_floor` and printed to the log. A
converged candidate and its fresh certification output must each satisfy
Q_minus <= max(1e-3, 1e-4*N_valence). Excess fails with the key path, measured
charge, limit and a cutoff/grid/dataset remedy. Divergence histories available
from the ordinary run path are written to failure results and a log as well.
No numerical kernel performs logging or JSON serialization.

Density thresholds are commonly used to avoid ill-conditioned low-density XC
and gradients. This policy explicitly distinguishes a functional floor from
changing the conserved density, and audits negative charge instead of silently
accepting it. No external source or implementation was read for this choice.

Synthetic grid tests exercise LDA-PZ, LDA-PW and PBE with negative combined
points, positive neighboring points, three finite-difference steps and unchanged
input arrays. The strictest-step maximum derivative error was
2.2271806621176893e-10 Ha Bohr^3 (cap 1e-7); clipped-point derivatives were zero.
The existing P13C NLCC oracles remain required.

## Initial density and remaining scope

This checkout starts SCF from `ScfProblem::uniform_density`, not PP_RHOATOM.
It already uses N_valence/volume and is positive. Provenance explicitly records
`method=uniform`, `atomic_table_used=false`. Therefore the radial atomic-table
charge deficits do not affect these bulk starts. An atomic-superposition start
and its rescaling have **not** been implemented. That requested extension needs
a validated atomic-density accessor/transform contract in dft-upf (outside this
ticket's owned source paths); its current `raw()` interface explicitly forbids
use as numerical-kernel input. No raw UPF table bypass was introduced.

The floor does not add polarized magnetic, metallic, force or stress support.
Explicit polarized floor inputs are rejected. Checkpoint v1 rejects floor-ledger
results because its model does not encode this new functional policy. The floor
is currently a `run`/ScfProblem feature; explicit floor inputs to bands/DOS and checkpoint/restart fail early with a
capability message rather than silently ignoring the policy.

README generation code and the dashboard registry are outside OWNED_FILES.
The architect must extend the README generator for this evidence type, generate
its output, and reconcile/upload the dashboard. The existing README freshness
test is run without hand-editing the generated block.

## Measured bulk results

Four requested CPU threads; wall times include CLI startup and I/O and were
measured on a shared host. No performance conclusion is drawn.

| System | Cutoff (Ha) | Iterations | Total energy (Ha) | Indirect gap (Ha) | Electrons | Wall (s) | Exit |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| si-20 | 20 | 15 | -8.46231783690604 | 0.043356604690 | 8.00000000000006 | 2.714 | 0 |
| si-30 | 30 | 15 | -8.46232900651845 | 0.043355757674 | 8.00000000000004 | 6.340 | 0 |
| si-40 | 40 | 15 | -8.46233682242691 | 0.043356644338 | 7.99999999999999 | 7.142 | 0 |
| mgo-40 | 40 | 18 | -76.05088085171062 | 0.234595439332 | 16.00000000000016 | 3.360 | 0 |
| mgo-50 | 50 | 19 | -76.05096924043031 | 0.234593031640 | 16.00000000000003 | 4.274 | 0 |
| mgo-60 | 60 | 19 | -76.05097717077302 | 0.234592354608 | 16.00000000000018 | 5.069 | 0 |

All six pass the standard fresh fixed-point certificate. Final clipped negative
charge is zero electrons in all cases. Total energies decrease over these
three cutoffs; this trend is not evidence of an asymptotic convergence tolerance.
The Mg potential has ten valence electrons, so MgO contains sixteen, not eight.
Both explicit unpolarized/floor and implicit unpolarized/default-NLCC-floor runs
produce identical density SHA-256, energies, and iteration counts. The ignored
real-file test verifies all six paired reductions (twelve CLI runs).

## Legacy and synthetic validation

The historical `docs/performance/baseline.json` comparison returned exit 5:
energy differences were at most 1.67e-13 Ha but stored density hashes differed.
A fresh build of the unmodified parent a3947c6e6499437837baf848d991c217ce533dd6
was therefore compared on the same toolchain/backend. Gamma, k222s, k222,
k444s and si8 each had exactly 0 Ha energy difference and bit-identical density
(`--require-identical`, exit 0). Neither the historical baseline nor its tolerance
was edited. This isolates this ticket's changes from pre-existing baseline drift.
Details and hashes are in `tests/reference/materials/legacy-equivalence.json`.
Single samples do not establish a speedup; the benchmark's printed ratios are
retained as raw output only.

The P13C release NLCC suite, including its ignored SCF/translation oracles,
passed all six tests. The new floor finite-difference and charge-budget tests
are included in the ordinary workspace gate. The cutoff, occupation and gap
rules remain unchanged.
