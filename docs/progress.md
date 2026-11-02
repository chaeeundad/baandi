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
