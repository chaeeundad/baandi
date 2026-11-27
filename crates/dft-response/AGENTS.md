# dft-response

**Purpose:** Response core (R2 onward of `docs/roadmap/rpa-crpa-plan.md`):
pair densities, chi0(q, i w) by the sum over states (decision D1), head and
wings at q -> 0 with the velocity d H_k / d k including the nonlocal
commutator, the symmetrized dielectric matrix and the macroscopic tensor
eps_M with and without local fields. No I/O, serde, logging or transport.

**Invariants to preserve and test:** pair FFT equals the convolution sum
(umklapp included); chi0(i w) Hermitian negative semidefinite; static
columns equal supercell finite differences of the occupied density
response; high-frequency closure over every empty state; k.p head and wings
equal finite-q extrapolation (nonlocal velocity included); one rank
bit-identical to the threads-only run; frequencies independent of the
thread count. Conventions: module docs of `src/lib.rs` and ADR 0054.

**Before editing:** read ../../AGENTS.md, ../../ownership.toml (allowed dependencies),
the plan and ADR 0054 in ../../docs/decisions/.
**Verify:** `scripts/check.sh dft-response`, then `scripts/check.sh` for the whole workspace.
Oracles: `cargo test -p dft-response --test oracles -- --nocapture` (debug, under 2 s);
larger free-electron mesh: `cargo test --release -p dft-response --test oracles -- --ignored --nocapture`.
Velocity vs central differences of H_k at fixed Miller set (R2-VEL-CHECK, production
`NonlocalOperator` rows, multi-projector s/p/d/f fixture): `cargo test -p dft-response --test velocity_fd`;
PseudoDojo LiF/AlN: `DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-response --test velocity_fd -- --ignored --nocapture`.
Diagnostic switch `BAANDI_RESPONSE_VELOCITY=full|local|sp` (engine; default full).
**Rules:** keep public API changes out of unrelated tickets; every behavior change adds or updates a test with stated tolerance and units; prefer general-purpose crates (ADR 0003); record design decisions in an ADR.

## Layout

- `pair`: one-rank batched stick FFT plans (`DistributedPairFft`) per
  coefficient set; real-space fields as full row-major grids (decouples the
  first-axis choices of orbital and response plans).
- `transitions`: (k, k + q) unit, occupied x block of empty bands; the MPI
  work unit of `dft-engine::response`.
- `chi0`: per-frequency GEMM accumulation `-(Q Q^H)`; frequencies split over
  threads (each frequency sequential); `reduce` = fixed-association tree sum.
- `velocity`: `<u_n| dH/dk |u_m>`; projector gradients from
  `dft_physics::stress::projector_k_gradient`.
- `dielectric`: optical tensors (Schur complement over G != 0) and finite q.
- `frequency`: mapped Gauss-Legendre grid (decision D5), `omega0` rule.

Pair fields, products, the scaled GEMM operand and the N_G^2 product are
reused scratch; the velocity matrices and the frequency worker threads are
allocated once per transition unit (never per transition or per G). No speed target applies (decision D10); costs are
measured and ledgered when tuned.
