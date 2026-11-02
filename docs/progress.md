# Progress log (상태 원천)

## 2026-11-02 환경 확인
- arm64 / macOS 27.0.1 / Xcode CLT / Accelerate.framework 확인.
- Rust 1.99.1 (rustup via Homebrew, keg-only; PATH=/opt/homebrew/opt/rustup/bin) 설치 완료, clippy/rustfmt 포함.
- pw.x(QE) 없음 → P6는 BLOCKED: QE reference unavailable.
- Codex CLI 0.166.0 (ChatGPT 로그인). `codex mcp-server`는 제거됨(≥0.154 추정, 공식 deprecated 2026-08-24).
  결정: `codex exec` + scripts/delegate.sh 방식 (사용자 승인). 컨텍스트는 AGENTS.md + docs + fork/resume.

## P0-BOOTSTRAP — 2026-11-02
- task/P0-BOOTSTRAP: nine std-only skeleton crates, Rust 1.99.1 pin, resolver 3, shared lints.
- Private Accelerate LP64 ddot/dsyev smoke tests: dot=32, eigenvalues=[1,3], residual=0, orthogonality error=2.2204460492503131e-16 (dimensionless).
- Required build/test/fmt/clippy commands executed successfully (exit 0); 10 unit tests. Initial clippy lint priority failure (101) corrected.
- Environment, dependency inventory, P1 API proposal and ADR 0001 recorded. Offline tests use a temporary rustup home referencing verified installed 1.99.1; see docs/environment.md.
- BLOCKED on commit: git add exited 128 because sandbox denied creation of /Users/mino/Develop/rust-dft/.git/worktrees/P0-BOOTSTRAP/index.lock. Implementation and required checks complete; no commit created. Complex ABI, FFT, QE reference and license selection remain outside P0.

## 2026-11-02 P0 verified
- Claude re-ran build, test, fmt, clippy in the P0 worktree: all exit 0 (10 unit tests). Merged to main as 7ce7209.
- Accelerate FFI verified for real ddot and dsyev only. Complex zheev/zgemm and LP64/ILP64 choice remain open (issue I8).
- Codex cannot write git metadata for a worktree; Claude commits after verification.
- Status dashboard data lives in artifacts/ (see artifacts/README.md).

## 2026-11-02 QE reference policy
- User decision: do NOT install pw.x on this Mac. Reference runs execute on dev08 through the dft5090 session; results are returned with command, QE version, input/UPF SHA-256 and wall time.
- Readiness check sent to dft5090 (read-only). P6 stays blocked until the reply and a valid reference setup.

## 2026-11-02 Policies added
- ADR 0003: prefer general-purpose libraries; P0B (hand-written complex FFI) stopped and discarded; spike P0C will choose faer vs LAPACK and ndrustfft.
- AI-first structure: AGENTS.md (read order), docs/INDEX.md, ownership.toml, tickets/, scripts/check.sh, per-crate AGENTS.md.
- Codex run registry: artifacts/sessions.json written by scripts/delegate.sh (PID, role, session, status, exit code); `delegate.sh stop` kills only the recorded PID.
- dft5090 reported it cannot reach dev08; waiting for user to supply a connection path.

## 2026-11-02 QE reference route ready (reported)
- dft5090 reported dev08 access, pw.x 7.6 (orig-mpi-dfti), idle node, scp exchange, and the same Si.pz-vbc.UPF digest on both of its hosts. Details and run protocol: docs/reference-environment.md. UPF provenance remains unverified (blank author/license).

## 2026-11-02 P0 closed, backend decided
- Spike P0C (Codex, session 01a19bfa) re-run by Claude in release: all tolerances held. ADR 0004 adopts faer (complex Hermitian eig, GEMM) and ndrustfft (3D FFT). Spike code kept at tag spike/p0c, results doc merged.
- P0 gates 5/5. Next: P1 (units, lattice, G-vectors, FFT on ndrustfft), which also deletes the hand-written FFI.

## 2026-11-02 P1A merged
- Codex P1A (session 01a19c0c): units, Cell/reciprocal lattice, FftPlan3d on ndrustfft, physics-conventions.md. Claude review found value-dependent shortcuts in the FFT that made the constant test vacuous; removed via resume (run #3), tests now tolerance-based.
- delegate.sh fixes found by real use: resume rejects -C/-s, stale report was read as result, relative prompt path after cd.
- Remaining for P1: G-vector cutoff, density/FFT grid rule, plane-wave<->grid mapping, CLI (ticket P1B).

## 2026-11-02 P1 closed
- P1B (Codex session 01a19c17) merged: PlaneWaveBasis, density sphere, FftGridSpec with 2,3,5 rule, scatter/gather, TOML config with deny_unknown_fields, `validate` and `grid` CLI. Run #1 stopped BLOCKED on manifest ownership (Claude fixed ndarray and Cargo.lock), run #2 finished.
- Si fcc a=5.43 A, 30 Ry: 749 plane waves, 5985 density G, FFT 25^3 (independent rule; compare with QE in P6).
- Next: P4 math pieces (eigensolver) can run now; P2 waits for the user's UPF decision.

## 2026-11-02 P4A and P3A merged, P2A started
- P4A (Codex session 01a19c21): faer adapter, LinearOperator, dense baseline, block Davidson; ADR 0006. Claude read the Davidson loop: convergence only from recomputed residuals.
- P3A (Codex session 01a19c2e): density, kinetic, Hartree, LDA-PZ, Ewald; ADR 0007 (xcx dev-only, no PZ81 in xcx so PZ is in-tree). Si Ewald -8.39947 Ha informational.
- P2A (UPF parser via upf crate, radial transforms, inspect-upf) delegated on synthetic fixtures; ADR 0008.

## 2026-11-02 P2A merged, Si candidate facts
- P2A (Codex session 01a19c3a) merged; Claude checked the UPF units against the specification (Rydberg atomic units; vloc and D in Ry; r*beta storage).
- dft5090 reported header facts of the candidate Si file (docs/candidate-si-upf.md): 2 projectors, D diagonal, tail -2Z/r matches, provenance blank. P2B (Hamiltonian assembly on synthetic pseudopotentials) delegated.

## 2026-11-02 P2B merged
- P2B (Codex session 01a19c4c) merged. Review question on projector norms (p stuck at 9.0e-3) settled by measurement: isolated-cell norm error 3e-13, Si image overlap predicted to 1e-13.
- P3 gates complete on synthetic pseudopotentials. Next: P5A (SCF, mixing, total energy audit, electrostatic cross-check) on synthetic PPs; real Si file waits for the user's decision.

## 2026-11-02 P5A merged
- P5A (Codex session 01a19c5d; run #1 hit a transient 'model at capacity' error, retried once via resume): SCF driver, LinearMixer and PulayMixer, term-by-term energy, ADR 0009.
- Gaussian electrostatic cross-check passes at 3.3e-14 over 48 cases; Claude re-derived the self-energy constant. This closes the main G=0 / Ewald convention risk on synthetic systems; QE energy-component comparison remains for P6.
- Next: P5B (CLI `run`, results.json, density file, exit codes), then the real-Si run once the user decides on the UPF.

## 2026-11-02 Real Si run (R1), P2 and P5 closed
- Si.pz-vbc.UPF transferred from dft5090 in three base64 chunks and verified (gz and original SHA-256). Real-file parsing needed P2C (relativistic="no", newline-separated arrays) and P5C (functional label with NOGX NOGC); the file was not given to Codex.
- Real Si SCF (Gamma only, 30 Ry, a=5.43 A): converged in 11 iterations, -7.295530468210927 Ha; degenerate levels preserved to 1e-12. Recorded in docs/validation.md (R1). Not yet compared with QE.
- Next: P6 (QE input generator, comparison script, cutoff/FFT/tolerance scans, reference runs through dft5090).

## 2026-11-02 First QE comparison (R3)
- Two real QE 7.6 runs (dev08 via dft5090) verified by chunk/gz/stdout SHA-256 and stored as tests/reference bundles. Strict tooling initially refused the real output (singular HOMO/LUMO label); fixed in P6B with a regression test on the real bundle.
- Comparison passes: total energy 1.28e-6 Ha/atom, Hartree/XC/Ewald <= 2e-8 Ha, one-electron 2.6e-6 Ha, eigenvalues <= 1.3e-6 Ha. Cause of the 2.6e-6 Ha one-electron residual unknown (issue I16).
- Requested from dft5090: QE cutoff scan (12 inputs; full stdout for the 6 same-grid runs, grid and energy lines for the 6 auto-grid runs).

## 2026-11-02 QE cutoff scan (R4)
- Six QE bundles (20 to 80 Ry) verified and compared: all pass; difference constant at 1.28e-6 Ha/atom. Diagnostic on the real pseudopotential traced the offset to the mesh-origin segment of the G=0 constant (predicted 2.548e-6 Ha, observed 2.556e-6 Ha); stated as consistent with QE, not proven.
- QE's automatic FFT grid equals ours at all scanned cutoffs. Remaining P6 item: density comparison; then final report and approval.

## 2026-11-03 Density comparison (R5), P6 gates complete
- pp.x is not installed on the runner; QE's native charge-density.dat was parsed read-only instead (Fourier coefficients, no grid dependence). Comparison tool P6D merged (90 Python tests). Result: rho(0) identical, sum |rho|^2 relative 1.1e-8, largest sampled coefficient difference 4.0e-9 e/Bohr^3.
- P6 comparison gates are all met. Remaining: documentation refresh (P6E), final review, and the user's decisions (license, MVP approval, UPF provenance statement).
