# AGENTS.md — rust-dft (Codex/구현 에이전트용 상시 지침)

Rust-first plane-wave DFT 엔진. 첫 성과는 빠른 계산이 아니라 **검증 가능한 계산**이다.

## 역할
- Claude Code = 아키텍트/검수자. Codex(당신) = 구현/테스트/수정 책임자. 사용자 = 최종 승인.
- 완료 선언은 신뢰되지 않는다. 실제 명령을 실행해 종료 코드와 수치 출력을 보고하라. 미실행 테스트를 성공으로 쓰지 마라.

## 환경
- Apple Silicon macOS. Rust stable 1.99 (rustup, `/opt/homebrew/opt/rustup/bin`이 PATH에 없으면 추가).
- 필수 통과 명령: `cargo build --workspace`, `cargo test --workspace`, `cargo fmt --check`,
  `cargo clippy --workspace --all-targets -- -D warnings`.
- 네트워크가 막힐 수 있다. 새 crate 추가가 필요하면 임의로 추가하지 말고 BLOCKED로 보고하라.

## 물리/수치 규약 (변경 금지, 상세는 docs/physics-conventions.md)
- 내부 단위: Hartree atomic units (Bohr, Hartree). 1 Ha = 2 Ry. 입출력에서만 Å/Ry/eV 변환.
- 역격자 B = 2π(A⁻¹)ᵀ, 격자벡터는 A의 열. ψ(r)=Ω^-1/2 Σ c(G)e^{iG·r}, Σ|c|²=1.
- Γ-only 복소수 계산. spin-unpolarized, norm-conserving, LDA만.
- 물리적 정합성 > 성능. 허용오차를 임의로 완화하지 마라. 불일치 원인 추적 순서:
  단위 → FFT 정규화 → G=0 → UPF projector → double counting → XC 파라미터화.

## 금지
- CUDA/Metal/MPI/PAW/USPP/PBE/힘/응력/smearing/다중 k-point (MVP 범위 밖).
- QE(GPL) 소스 복사·번역. UPF 파서는 공개 명세로 독립 구현.
- UPF 파일을 Git에 커밋. 실제 계산으로 얻지 않은 수치를 'QE reference'로 기록.
- 소유 파일(OWNED_FILES) 밖 수정. 공통 타입/Cargo workspace manifest/public API는 티켓이 허락한 경우에만.
- 핵심 수치 루프에 로깅/JSON/파일 I/O. unsafe는 FFI 어댑터에만, 전제조건 주석+테스트 필수.
- 원격 서버, 대량 설치, 파괴적 Git 작업(reset --hard, force push 등).

## 상태의 원천
`docs/progress.md`, `docs/decisions/`(ADR), 티켓 ID, Git 커밋만이 상태 원천이다. 설계 결정은 ADR에 남겨라.

## 완료 보고
최종 메시지는 `docs/schemas/report.schema.json` 형식을 따른다 (status, changed files, 테스트 명령+종료 코드,
수치 출력(단위 포함), 한계/위험, 열린 질문, 커밋 해시). 작업 브랜치에 커밋하고 해시를 보고하라.
