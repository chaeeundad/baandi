# ADR 0039: opt-in atomic starting wavefunctions

Status: Proposed.
Ticket: P46A-ATOMIC-START-OPTION.

P45A reduced standard SCF counts to 13 for Si8 and 9 for Si16 on both
arithmetic paths, but failed its frozen skew-cell spread and accuracy gates.
Those results do not justify changing defaults. The architect authorizes an
opt-in implementation and independent evaluation under preregistered R2.

Options are retaining only free-electron starts or adding atomic starts with
explicit fallback and initialization work metadata. The latter must preserve
supplied restart/warm vectors, default arithmetic and exact-restart identity.
Atomic radial functions are chi(r)=r R(r), with all magnetic components and
all supplied occupations represented; occupations do not weight the basis.
The initialization uses the first-density Hamiltonian, with unchanged density,
mixing and stopping rules. Gamma completion must preserve real-function
conjugate symmetry before ORC conversion.

Implementation and evidence are tracked in
[atomic-start-option.md](../performance/atomic-start-option.md). Recommendation
is conditional on R2. No recommendation or acceptance is recorded yet.

Local verification now reproduces P45A's 13/9 standard counts on both paths,
passes 70 default identity comparisons and 48 atomic thread fingerprint
comparisons, and passes workspace build/test/fmt/clippy. This qualifies the
local implementation contracts only. At that initial stage, new-system coverage and full-cost acceptance were
pending the clarifications subsequently resolved below.
Status remains Proposed and defaults remain free_electron.

## R2 evaluation and recommendation

The architect clarified before new data that si32 uses five standard starts
and one strict reference per arm; other systems use ten and two. Total cost
is median single-thread whole-SCF wall time with alternating paired arms;
H applications and MAC counts are diagnostics. See the frozen protocol snapshot.

R2 fails on the third system, si-k222: atomic maximum standard energy error
8.8816420884541e-10 Ha is 15.682825463 times the default error, above the
factor-3 limit, despite SCF 10 -> 9 and local median wall 0.3089391875 ->
0.2814508955 s. The earlier si32 and k444s results do not override this failure.
The mandatory stop leaves MgO, PseudoDojo Si and the low-symmetry cell unmeasured.

Recommendation: do not recommend atomic starts under R2 and do not change the
free_electron default. Retain the separately authorized opt-in implementation
for review with this negative evidence. Status remains Proposed; no universal
accuracy or performance improvement is established. No tuning or new remote
measurement follows from this ADR.
