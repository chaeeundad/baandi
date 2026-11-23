# Request: run the full gate on dev08 for a pushed commit (GATE-DEV08, 2026-11-19)

Purpose: the Mac is shared by many sessions (load 30 to 100, 16 GB RAM with swap nearly full), so a
`scripts/check.sh full` there takes 5 to 50 minutes depending on load. The Mac runs only
`scripts/check.sh fast` while iterating (docs/performance/fast-loop-protocol.md, "Two lanes"); the
pre-push `full` gate of a pushed commit can run on dev08 instead, through the master session.
This is a correctness gate, not a performance measurement: no QE, no timing table, no load gate
below 0.30 (the gate's own load line is reported). Job class: `gate-<short hash>`; it must finish
within the 10-minute stage limit (a warm persistent target directory makes it about 1 to 3
minutes; a cold first build on dev08 is the only part that may come near the limit, see step 2).

## Inputs from the requester
- `HASH`: the full commit hash, already pushed to `origin` (`https://github.com/chaeeundad/baandi.git`).
- Optional `SLOW=1`: also run the release-mode ignored tests listed under "Tests ignored for runtime"
  in docs/performance/slow-tests.md (currently none; the step is then a no-op).

## Commands (dev08, master session)
```
set -euo pipefail
HASH=<full hash>
GATE=$HOME/baandi-gate                         # one clean checkout reused across gates
CACHE=$HOME/baandi-gate-target/1.99.1-dev      # persistent target: toolchain + profile
SIG=<shared NFS dir>/signals                   # the usual signal directory
echo "$(date -Is) gate-${HASH:0:8} started" > $SIG/gate-${HASH:0:8}.started

# 1. clean checkout of exactly HASH (no local changes survive: reset + clean)
[ -d $GATE/.git ] || git clone https://github.com/chaeeundad/baandi.git $GATE
git -C $GATE fetch origin
git -C $GATE checkout --force --detach $HASH
git -C $GATE clean -fdx -e target
test "$(git -C $GATE rev-parse HEAD)" = "$HASH"

# 2. pinned toolchain (rust-toolchain.toml: 1.99.1, rustfmt, clippy); rustup installs it once
cd $GATE && rustup show active-toolchain && cargo --version
#    no network on the build step? vendor once on the Mac and sync it:
#      (Mac) cargo vendor --locked ../rust-dft-vendor && rsync -a ../rust-dft-vendor dev08:baandi-vendor
#      (dev08) export CARGO_HOME=$HOME/baandi-gate-cargo; mkdir -p $CARGO_HOME
#              printf '[source.crates-io]\nreplace-with="v"\n[source.v]\ndirectory="%s"\n' \
#                     $HOME/baandi-vendor > $CARGO_HOME/config.toml
#      check.sh calls cargo plainly; the config in CARGO_HOME redirects crates.io to the
#      vendor directory, so the build needs no network.

# 3. the gate itself (same script as the Mac; a persistent target keeps it warm)
export CARGO_TARGET_DIR=$CACHE TMPDIR_990=$HOME/baandi-gate-tmp
export CHECK_TEST_JOBS=8 CHECK_NICE=0
t0=$(date +%s); uptime
scripts/check.sh full 2>&1 | tee $HOME/gate-${HASH:0:8}.log; rc=${PIPESTATUS[0]}
echo "GATE_WALL $(( $(date +%s) - t0 )) s rc=$rc"; uptime

# 4. optional: the slow ignored tests in release (only with SLOW=1; list in slow-tests.md)
#    cargo test --release -p <crate> --test <target> <name> -- --ignored --exact

R=$([ $rc = 0 ] && echo done || echo failed)
echo "$(date -Is) rc=$rc log=$HOME/gate-${HASH:0:8}.log" > $SIG/gate-${HASH:0:8}.$R
```
`TMPDIR_990` only names the scratch directory for the gate's logs (on the Mac it points at the
990-2TB volume); on dev08 any local directory works. Equivalent without the script, for a
cross-check: `cargo build --workspace && cargo test --workspace --no-fail-fast && cargo fmt --check
&& cargo clippy --workspace --all-targets -- -D warnings` (same coverage; check.sh only runs the
test binaries concurrently).

## Report back (to the orchestrator)
- `HASH`, `git rev-parse HEAD` of the checkout, `rustc --version`, `cargo --version`.
- The gate output in full: the `LOAD` line, the four `PASS`/`FAIL` lines, the runner summary
  (`N binaries, P passed, F failed, I ignored, W s wall`), the `WALL` line, `GATE_WALL` and rc.
- On a failure: the failing binaries' tails as printed by check.sh (they are in the log).
- Whether the target directory was cold (first gate or after a toolchain change) or warm.
- Signal files: `gate-<hash>.started`, then `.done` or `.failed`; `dev08.idle` as usual when the
  queue is empty.

## Rules
- The checkout is reset to `HASH` before every gate; nothing is edited on dev08.
- Never mixed with performance jobs on the same cores at the same time: the gate takes the cores
  the queue gives it, and a timed job never runs beside a gate.
- A gate result is valid for `HASH` only. A merge commit is a new hash and gets its own gate.


## Update 2026-11-23: who runs the gate
Since 2026-11-23 the measurement queue and the gate are operated by the `baandi-dev08` session (the dft5090 relay is retired). The gate runs on socket 1 from `/scr/mino/gate/gate.sh` (not `$HOME/baandi-gate`); measurements run on socket 0. A gate job in the queue is named `final-a0-gate-<hash>` and calls that script, so no measurement runs next to a gate. Paths for artifacts stay under `/scr` (the root filesystem under `$HOME` is small).
