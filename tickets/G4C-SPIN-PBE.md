TASK_ID: G4C-SPIN-PBE
GOAL: |
  Spin-polarized PBE (with and without NLCC) in the reference driver. The reference spin entry (crates/dft-engine/src/refscf/electronic.rs) rejects it: "distributed PBE has no two-channel evaluator". Real pseudopotential sets available to us (PseudoDojo) are PBE, so spin on real materials in the reference driver needs this. The legacy path has a serial spin PBE (crates/dft-physics/src/spin.rs SpinFunctional::Pbe, collinear.rs) that is the oracle.
READ FIRST: docs/design/reference-scf.md ("PBE", "NLCC", "Spin and Fermi-Dirac"); crates/dft-physics/src/gga/distributed.rs (unpolarized distributed GGA: gradients of the owned field, divergence term) and gga.rs/pbe_pw92.rs; crates/dft-physics/src/spin.rs (serial polarized PBE, core split rule: half the core in each channel); refscf/electronic.rs and electronic_map.rs (how the spin XC is called, LDA-PW case); G2B PbeCore in refstart.rs. Never read QE source.
TASK: |
  1. Distributed two-channel PBE: gradients of n_up+c/2 and n_down+c/2 (and whatever the serial spin PBE needs: total gradient, per-channel gradients) on the owned slabs with the existing distributed FFT; energy density and potentials per channel including the divergence term; same signed-density and threshold policy as the serial spin PBE. Reuse scratch across maps (prepared once per run).
  2. Wire it into the reference spin entry for fixed sector and Fermi-Dirac (scalar is unaffected), including NLCC; remove the rejection for PBE (keep polarized LDA-PZ rejected, as legacy).
  3. Tests (Mac): distributed vs serial spin PBE potentials and energy on random positive densities with and without core (1e-13 relative), ranks 1/3/8; reference vs legacy delivered energy on a synthetic open-shell PBE fixture (fixed sector, and FD with the G4B switches on) within the standard class and 1e-10 Ha tight; unpolarized and LDA paths bit-identical.
OWNED_FILES:
  - crates/dft-physics/src/gga/distributed.rs
  - crates/dft-physics/src/gga.rs   # shared kernels only if needed
  - crates/dft-engine/src/refscf/electronic.rs
  - crates/dft-engine/src/refscf/electronic_map.rs
  - crates/dft-engine/tests/refscf.rs   # new tests appended at the end only
  - crates/dft-physics/tests/   # new test file only
  - crates/dft-cli/src/run.rs   # remove the CLI-side PBE spin rejection only
  - docs/design/reference-scf.md   # "Spin and Fermi-Dirac" section only
TEST_COMMANDS:
  - cargo check --workspace --tests
  - cargo test -p dft-physics
  - cargo test -p dft-engine --test refscf
  - cargo fmt --check && cargo clippy -p dft-physics -p dft-engine -p dft-cli --all-targets -- -D warnings
ACCEPTANCE_CRITERIA: spin PBE (with NLCC) runs in the reference driver with legacy agreement; other paths bit-identical.
GIT: leave changes uncommitted.
RETURN_FORMAT: docs/schemas/report.schema.json (report in English)
