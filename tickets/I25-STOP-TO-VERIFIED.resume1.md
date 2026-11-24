Approved (orchestrator): add exactly these to OWNED_FILES and apply your integration patch (target/tmp/stop-to-verified-integration.patch):
  - crates/dft-cli/src/main.rs   # option wiring for the stop-to-verified reader only
  - crates/dft-cli/src/run.rs    # JSON serialization of the new per-map timestamps only (reference-driver output section)
  - crates/dft-engine/tests/common/rank_rule.rs   # clear the new timestamp fields in fingerprint() (they are wall clocks, like inner_clocks)
The timestamps must be present at every diagnostic level (lightweight default included). Then rerun the Mac commands of the ticket (cargo check, the two test commands, fmt, clippy -D warnings) and report. The full gate and release runs stay with the orchestrator on dev08. Leave changes uncommitted.
