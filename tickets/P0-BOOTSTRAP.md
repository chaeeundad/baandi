TASK_ID: P0-BOOTSTRAP
GOAL: Apple Silicon Mac mini에서 Rust-first plane-wave DFT 엔진의 재현 가능한 CPU 개발 환경과 Cargo workspace를 구성한다.
OWNED_FILES: Cargo.toml, rust-toolchain.toml, README.md, crates/* (최소 skeleton), docs/architecture.md, docs/environment.md, docs/dependencies.md, docs/decisions/*, docs/progress.md(P0 항목 추가만), .github/ 는 만들지 말 것.
READ_ONLY_CONTEXT: AGENTS.md, docs/progress.md, docs/schemas/report.schema.json
CONSTRAINTS:
- 샌드박스에서 네트워크가 막혀 있다. P0에서는 외부 crate 의존성 0개(std만)로 구성하라. crate 추가가 필요해 보이면 추가하지 말고 open_questions에 사유와 후보를 적어라.
- workspace crate 구조: dft-units, dft-geometry, dft-numerics, dft-basis, dft-upf, dft-physics, dft-solvers, dft-engine, dft-cli (빈 skeleton + 컴파일되는 최소 단위 테스트). 과도하다고 판단되면 축소안을 open_questions로 제안.
IMPLEMENTATION:
1. 환경 조사 보고(docs/environment.md): arch, rustc/cargo, Xcode CLT, Accelerate, pw.x/fftw/openblas 유무 (실제 실행 결과만 기록).
2. rust-toolchain.toml로 버전 고정(현재 1.99.1), workspace Cargo.toml(resolver, 공통 lints), .gitignore는 이미 있음.
3. dft-numerics에 Apple Accelerate ABI smoke test: build.rs 또는 #[link(name="Accelerate", kind="framework")]로 CBLAS(zgemm 또는 ddot)와 LAPACK(zheev 또는 dsyev) 호출을 FFI로 수행. 작은 알려진 행렬(예: 2x2 Hermitian, 고유값 해석해)로 결과 검증. LP64/ILP64 정수 폭, column-major, 복소수 호출 규약, 스레딩을 docs/environment.md에 기록. unsafe는 이 FFI 어댑터에만 격리하고 전제조건을 주석+테스트로 명시.
4. 다음 단계(P1)용 가장 작은 수치 백엔드 API 스케치를 docs/architecture.md에 제안 (구현 아님, 트레이트 시그니처 수준). FFT 후보(FFTW, rustfft 등)의 라이선스/유지보수/Apple Silicon 이슈는 네트워크 없이 아는 범위만 쓰고 "미검증"으로 표기.
5. docs/decisions/0001-*.md: ADR 템플릿과 첫 결정(수치 백엔드 격리 방침). docs/dependencies.md: 현재 의존성 0개 보고.
NON_GOALS: DFT 물리, FFT 구현, CUDA, CI 워크플로, LICENSE 확정.
TEST_COMMANDS (모두 실행하고 종료 코드 보고; PATH에 /opt/homebrew/opt/rustup/bin 추가):
cargo build --workspace; cargo test --workspace; cargo fmt --check; cargo clippy --workspace --all-targets -- -D warnings
ACCEPTANCE: 위 4개 성공, Accelerate FFI 테스트 통과(수치 출력 포함), 변경은 task 브랜치에 커밋.
RETURN_FORMAT: report schema. 미실행 테스트를 성공으로 보고하지 말 것.
