TASK_ID: P5C-FUNCTIONAL-STRING
GOAL: Fix the UPF functional-string compatibility check so that real-world norm-conserving LDA files are accepted. Real files written by older generators declare the functional as `SLA  PZ   NOGX NOGC` (Slater exchange, Perdew-Zunger correlation, and explicit "no gradient correction" markers for exchange and correlation, with leading/double spaces). The current check in two places requires the whitespace-split tokens to equal exactly ["SLA", "PZ"], so the real file is rejected with `unsupported UPF functional ... requires SLA PZ`. The functionals are the same LDA-PZ the engine implements.
OWNED_FILES: crates/dft-upf/** (add one shared helper and its tests), crates/dft-engine/** (use the helper, tests), crates/dft-cli/** (use the helper, tests), docs/upf-support.md and docs/cli.md (update the functional paragraph). Do not touch manifests, ADRs, AGENTS.md or other crates.
READ_ONLY_CONTEXT: AGENTS.md, crates/dft-cli/src/run.rs (the check near the UPF parsing, error code 2), crates/dft-engine/src/scf.rs (the check "only explicitly labelled SLA PZ potentials are supported"), crates/dft-upf/src, docs/upf-support.md, docs/cli.md.
REQUIRED BEHAVIOR:
1. One shared function in dft-upf (for example `functional_is_lda_pz(&str) -> bool`, or a method on the validated pseudopotential type) that tokenizes on ASCII whitespace (so leading, trailing and repeated spaces do not matter) and accepts exactly two token sequences, case-sensitive as written in UPF files: ["SLA", "PZ"] and ["SLA", "PZ", "NOGX", "NOGC"]. Everything else is rejected: other correlation (["SLA","PW"]), gradient corrections (["SLA","PZ","PBX","PBC"], ["PBE"]), partial markers (["SLA","PZ","NOGX"]), reordered tokens, empty string. Do not add fuzzy matching.
2. dft-cli (run pipeline) and dft-engine (SCF input validation) both call this helper; remove the duplicated token comparisons. Error messages stay explicit: name the offending functional string and list the two accepted forms. Exit code 2 behavior in the CLI is unchanged.
3. The stored functional string stays exactly as in the file in results.json (do not normalize it there), so the recorded value matches the UPF file.
TESTS (no real UPF; all synthetic fixtures):
a. Unit tests of the helper over accepted forms (including " SLA  PZ   NOGX NOGC", tabs, CRLF remnants if tokenization sees them) and the rejected forms above.
b. CLI end-to-end: a synthetic UPF with functional="SLA  PZ   NOGX NOGC" runs and converges; results.json records the original string; a synthetic UPF with functional="SLA PW" (or "PBE") fails with exit code 2 and the cause on stderr.
c. Engine-level test with both accepted forms giving bitwise-identical total energies for the same synthetic data.
d. All existing tests still pass unchanged (the existing negative functional test in the CLI run tests must keep passing, adapt only if its expected message text changes and say so).
TEST_COMMANDS: scripts/check.sh; `cargo test --offline -p dft-upf -p dft-engine -p dft-cli -- --nocapture`.
ACCEPTANCE_CRITERIA: all gates exit 0; no unsafe; no new dependencies; no value-dependent shortcuts; docs match the code.
NON_GOALS: other functionals, GGA, any physics change, touching the real UPF.
DEPENDENCIES: P2C merged (done).
GIT: The sandbox cannot write git metadata for the worktree. Leave changes uncommitted and say so.
RETURN_FORMAT: docs/schemas/report.schema.json. Real exit codes only.
