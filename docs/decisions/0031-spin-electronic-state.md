# ADR 0031: Tagged spin and ElectronicState

- Status: **proposed**, 2026-11-05; P13A-STATE-CONTRACTS / MF02. User approval is required before MF03.
- Authority: [master plan D1, D3, D4, MF03, MF07](../roadmap/master-plan-proposal.md), scope B in [AGENTS.md](../../AGENTS.md). This is a contract, not implemented capability.
- Dependencies: [ADR 0016](0016-esle-layering.md) remains proposed; [P11C](../checkpoint-format.md) is the concrete persistence baseline. [Evidence map](../validation/evidence-map.md) supersedes older survey statements about missing restart support.

## Proposed representation

The pure engine owns the following logical types; serialization stays in the CLI.
All floating values are finite binary64; lengths are checked unsigned integers.
These are field specifications, not Rust source or permission to alter existing structs.
The shared records and fingerprint rules are detailed in [state-contracts](../architecture/state-contracts.md#shared-records).

| Type | Complete fields / invariant |
| --- | --- |
| `SpinMode` | Tagged `unpolarized`, `collinear{axis_cartesian:[3], constraint}`, or `noncollinear_reserved`. Axis has unit Euclidean norm within 1e-12; channel order is up then down along that axis. Reserved mode fails construction. |
| `SpinConstraint` | Tagged `none` or `fixed_total_moment{electron_difference:f64}`. Fixed moment sets N_up=(N+M)/2 and N_down=(N-M)/2; abs(M)<=N. No local constraint fields accepted. |
| `BandLayout` | `spin_channels:u64`, `k_count:u64`, `band_count:u64`, `channel_capacity_electrons:[spin_channels]`, ordered `k_points:[k_count]`. Same band count across spin/k initially. Unpolarized: 1 channel, capacity 2; collinear: 2 channels, capacity 1 each. |
| `BandData` | `layout`, `eigenvalues_ha:[s,k,b]`, `occupations_electrons:[s,k,b]`, `eigen_residuals_ha:[s,k,b]`. Last axis contiguous. Occupations already include spin capacity, never k weight. No hidden factor of two. |
| `OrbitalBlock` | `spin_index`, `k_index`, `basis_order_digest`, `g_count`, `band_count`, `coefficients:[g,b]` complex128, column-major within a block. Each column norm 1. Blocks ordered spin-major then k; G count may differ by k. Optional as a complete payload, not fictitious zero orbitals. |
| `SpinDensity` | Tagged `unpolarized{grid_id, total:[x,y,z]}` or `collinear{grid_id, channels:[2,x,y,z]}`. C order, z contiguous, electron number/Bohr^3. Total and magnetization are derived views, not separately mutable authoritative arrays. |
| `MagnetizationDensity` | Collinear signed scalar m=n_up-n_down along the recorded axis, electron difference/Bohr^3; total moment M=integral m, reported numerically in mu_B using the conventional spin-polarization sign. This is not signed electron angular magnetic moment or an orbital moment. Unpolarized m=0 by definition. Future spinor n,m_x,m_y,m_z requires a separate accepted variant and positivity test n>=abs(m). |

For each channel, n_s(r)=sum_kb w_k f_skb abs(psi_skb(r))^2;
N_s=Omega/N_grid sum_r n_s(r). Total N=sum_s N_s=sum_skb w_k f_skb;
M=N_up-N_down. Weights are positive and sum to one within 1e-14, unchanged
from [physics conventions §26](../physics-conventions.md). In unpolarized mode
the single density is total n, not n/2. In collinear mode both channel densities
must be feasible under the selected XC's existing negative-density rule; never
clip magnetization or independently rescale channels after eigensolving.

## Proposed XC boundary

`SpinXcRequest` contains `functional_id`, `functional_version`, `spin_mode`,
`density`, `grid_descriptor`, and `core_density: absent | scalar_grid`.
Core density is immutable, electron/Bohr^3, split equally between channels for
nonmagnetic NLCC only after MF06 validates the functional/core convention.
Requesting it earlier rejects. `SpinXcResponse` contains `energy_xc_ha`,
`potential_ha:[s,x,y,z]`, `valence_potential_integral_ha` (sum_s integral n_s V_s),
and `derivative_capabilities`. Energy includes the chosen core convention; the
band double-counting integral uses valence density only.

For GGA the library adapter consumes n_up,n_down and three invariants
sigma_uu, sigma_ud, sigma_dd (Bohr^-8). Its outputs are energy per volume
(Ha/Bohr^3), partial derivatives with respect to both densities (Ha), and
three sigma derivatives (Ha Bohr^5). The adapter explicitly converts a library's
per-electron energy convention. With the skew-adjoint D of
[ADR 0018](0018-pbe-via-xcx.md), V_up = e_nup - div(2 e_sigma_uu grad n_up +
e_sigma_ud grad n_down), and analogously down. Cross derivative multiplicity
must be tested; it is not guessed from an unpolarized API. MF03 must authorize
and independently verify the polarized XC backend; this ADR does not assert
that the existing xcx wrapper supplies it. LDA omits gradients.

## Identity, compatibility and proof obligations

Spin mode, axis, constraints, XC convention and occupation policy belong to the
model fingerprint. Actual densities, occupations, moments and optional orbital
hashes belong to state identity. Ordering and capacities belong to basis identity.
Initial magnetic seeds are provenance, not a new Hamiltonian. Opposite moments
remain distinct state records even if their energies match. A global spin flip
may be reported as a symmetry relation only under MF11, never silently merged.

Unpolarized legacy Gamma and k APIs keep their arithmetic and first-k views.
A new spin result view is additive; accessing a legacy scalar orbital view on a
collinear result rejects rather than selecting up implicitly. Initial MF03 fixed
occupations require explicit integer channel electron counts (1/0 per channel);
this selects a fixed-moment sector. Unconstrained moment equilibration requires
MF04's shared chemical potential. No odd-electron exception is added to legacy
unpolarized constructors. Mixing acts on both channels with one feasibility step;
raw charge **and** raw spin residuals certify convergence, independently of any
Kerker filter ([ADR 0025](0025-preconditioned-density-mixing.md)).

For new fixed-collinear runs, each nonempty occupied channel needs at least one
empty band and its own sampled indirect gap above max(requested gap, twice the
eigen tolerance). A zero-electron channel has no occupied-edge gap; record that
channel's gap as not applicable, still solve its retained spectrum. Cross-channel
level ordering does not invalidate a deliberately fixed-M sector. This is not
the unconstrained ground-state certificate; the sector constraint is reported.

Acceptance: polarized uniform exchange analytic derivatives; channel swap;
zero-magnetization collinear versus unpolarized E, density and spectrum; N and M
integrals against weighted occupations; unequal k weights; nonzero/zero moment;
odd collinear electron count; cubic, skew and elongated cells. Exact proposed
limits and failures are in [MF03 oracles](../validation/stage-oracles.md).

Rejected alternative: duplicated spin solvers or a spinor flag implemented by
extra bands (D1). Consequence: explicit shapes and capability dispatch add boundary
validation, but do not change numerical dependencies or authorize SOC.
