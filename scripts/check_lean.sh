#!/usr/bin/env bash
# Reproducible CI verification for BAANDI Lean 4 + Mathlib formal proofs.
#
# Checks:
#   1. Static scan of all project .lean sources (excluding .lake/ and Lean comments)
#      for forbidden `sorry`, `admit`, or custom `axiom` declarations.
#   2. `lake build` in the target Lean project directory.
#   3. Re-elaboration of every module containing `#print axioms` via `lake env lean`
#      so cached `.olean` builds cannot skip the axiom audit, verifying that every
#      expected theorem is audited and depends ONLY on standard Lean/Mathlib axioms:
#        {propext, Classical.choice, Quot.sound}.
#
# Usage:
#   bash scripts/check_lean.sh [LEAN_PROJECT_DIR]
#   LEAN_PROJECT_DIR=/scr/mino/lean/BaandiMath bash scripts/check_lean.sh

set -uo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "${SCRIPT_DIR}/.." && pwd)"
PROJECT_ARG="${1:-${LEAN_PROJECT_DIR:-formal/BaandiMath}}"

if [[ "${PROJECT_ARG}" = /* ]]; then
  PROJECT_DIR="${PROJECT_ARG}"
else
  PROJECT_DIR="${REPO_ROOT}/${PROJECT_ARG}"
fi

# Make elan-managed lean/lake visible in non-interactive shells (e.g. dev08).
if [[ -d "${HOME}/.elan/bin" ]]; then
  export PATH="${HOME}/.elan/bin:${PATH}"
fi

TOOLCHAIN="unknown"
if [[ -f "${PROJECT_DIR}/lean-toolchain" ]]; then
  TOOLCHAIN="$(tr -d '[:space:]' < "${PROJECT_DIR}/lean-toolchain")"
fi

fail() {
  local reason="$1"
  local code="${2:-1}"
  echo "LEAN_CHECK_FAIL reason=${reason} toolchain=${TOOLCHAIN} project=${PROJECT_DIR} exit_code=${code}"
  exit "${code}"
}

if [[ ! -d "${PROJECT_DIR}" ]]; then
  echo "ERROR: Lean project directory not found: ${PROJECT_DIR}" >&2
  fail "project_dir_missing" 2
fi

if ! command -v python3 >/dev/null 2>&1; then
  echo "ERROR: python3 is required for comment-aware Lean source and axiom auditing." >&2
  fail "python3_missing" 2
fi

# Step 1: Comment-aware static scan for forbidden sorry / admit / custom axiom declarations,
# and collection of all modules containing `#print axioms`.
TMP_DIR="$(mktemp -d)"
trap 'rm -rf "${TMP_DIR}"' EXIT
EXPECTED_AXIOMS_FILE="${TMP_DIR}/expected_axioms.txt"
AUDIT_MODULES_FILE="${TMP_DIR}/audit_modules.txt"

python3 - "${PROJECT_DIR}" "${EXPECTED_AXIOMS_FILE}" "${AUDIT_MODULES_FILE}" <<'PY'
import pathlib
import re
import sys

project_dir = pathlib.Path(sys.argv[1])
expected_out = pathlib.Path(sys.argv[2])
modules_out = pathlib.Path(sys.argv[3])

def strip_lean_comments_and_strings(src: str) -> str:
    """Replace Lean line comments (-- ...), nested block comments (/- ... -/),
    and string literals ("...") with spaces while preserving line numbers."""
    out = []
    i = 0
    n = len(src)
    block_depth = 0
    in_line_comment = False
    in_string = False

    while i < n:
        ch = src[i]
        nxt = src[i + 1] if i + 1 < n else ""

        if in_line_comment:
            if ch == "\n":
                in_line_comment = False
                out.append("\n")
            else:
                out.append(" ")
            i += 1
            continue

        if block_depth > 0:
            if ch == "/" and nxt == "-":
                block_depth += 1
                out.append("  ")
                i += 2
            elif ch == "-" and nxt == "/":
                block_depth -= 1
                out.append("  ")
                i += 2
            elif ch == "\n":
                out.append("\n")
                i += 1
            else:
                out.append(" ")
                i += 1
            continue

        if in_string:
            if ch == "\\" and i + 1 < n:
                out.append("  ")
                i += 2
            elif ch == '"':
                in_string = False
                out.append(" ")
                i += 1
            elif ch == "\n":
                out.append("\n")
                i += 1
            else:
                out.append(" ")
                i += 1
            continue

        if ch == "-" and nxt == "-":
            in_line_comment = True
            out.append("  ")
            i += 2
        elif ch == "/" and nxt == "-":
            block_depth = 1
            out.append("  ")
            i += 2
        elif ch == '"':
            in_string = True
            out.append(" ")
            i += 1
        else:
            out.append(ch)
            i += 1

    return "".join(out)

lean_files = sorted(
    p for p in project_dir.rglob("*.lean")
    if ".lake" not in p.parts and "build" not in p.parts
)
if not lean_files:
    print(f"ERROR: No .lean files found under {project_dir}", file=sys.stderr)
    sys.exit(2)

violations = []
expected_axioms = []
audit_modules = []

sorry_re = re.compile(r"\b(sorry|admit)\b")
axiom_decl_re = re.compile(r"^\s*(?:(?:private|protected|noncomputable|unsafe)\s+)*axiom\b")
print_axioms_re = re.compile(r"^\s*#print\s+axioms\s+([A-Za-z0-9_.']+)")

for path in lean_files:
    rel = path.relative_to(project_dir)
    text = path.read_text(encoding="utf-8")
    stripped = strip_lean_comments_and_strings(text)
    has_print_axioms = False
    for lineno, line in enumerate(stripped.splitlines(), start=1):
        m_sorry = sorry_re.search(line)
        if m_sorry:
            violations.append(f"{rel}:{lineno}: forbidden token `{m_sorry.group(1)}`")
        if axiom_decl_re.match(line):
            violations.append(f"{rel}:{lineno}: forbidden custom `axiom` declaration")
        m_pa = print_axioms_re.match(line)
        if m_pa:
            has_print_axioms = True
            expected_axioms.append(f"{rel}\t{m_pa.group(1)}")
    if has_print_axioms:
        audit_modules.append(str(rel))

if violations:
    print("ERROR: Forbidden Lean constructs detected:", file=sys.stderr)
    for v in violations:
        print(f"  - {v}", file=sys.stderr)
    sys.exit(3)

if not expected_axioms:
    print("ERROR: No `#print axioms` audit statements found in project.", file=sys.stderr)
    sys.exit(4)

expected_out.write_text("\n".join(expected_axioms) + "\n", encoding="utf-8")
modules_out.write_text("\n".join(audit_modules) + "\n", encoding="utf-8")
print(f"[static-scan] OK ({len(lean_files)} .lean files scanned, {len(expected_axioms)} `#print axioms` entries found)")
PY
STATIC_RC=$?
if [[ "${STATIC_RC}" -ne 0 ]]; then
  fail "static_scan_failed" "${STATIC_RC}"
fi

if ! command -v lake >/dev/null 2>&1; then
  echo "ERROR: 'lake' not found in PATH (checked ${HOME}/.elan/bin). Run on dev08 or install elan." >&2
  fail "lake_not_found" 127
fi

LEAN_VER="$(cd "${PROJECT_DIR}" && lean --version 2>/dev/null | head -n 1 || true)"
echo "[toolchain] ${TOOLCHAIN} (${LEAN_VER:-unknown})"

# Step 2: Build the Lean project with `lake build`.
BUILD_LOG="${TMP_DIR}/lake_build.log"
echo "[lake-build] Running lake build in ${PROJECT_DIR}..."
(cd "${PROJECT_DIR}" && lake build) 2>&1 | tee "${BUILD_LOG}"
BUILD_RC="${PIPESTATUS[0]}"
if [[ "${BUILD_RC}" -ne 0 ]]; then
  fail "lake_build_failed" "${BUILD_RC}"
fi

# Step 3: Run `lake env lean` on every module with `#print axioms` to guarantee fresh
# axiom output even when `lake build` reuses cached `.olean` files.
AXIOM_LOG="${TMP_DIR}/axioms.log"
: > "${AXIOM_LOG}"
while IFS= read -r rel_mod; do
  [[ -z "${rel_mod}" ]] && continue
  echo "[axiom-audit] Elaborating ${rel_mod} via lake env lean..."
  (cd "${PROJECT_DIR}" && lake env lean "${rel_mod}") 2>&1 | tee -a "${AXIOM_LOG}"
  MOD_RC="${PIPESTATUS[0]}"
  if [[ "${MOD_RC}" -ne 0 ]]; then
    fail "lean_elab_failed" "${MOD_RC}"
  fi
done < "${AUDIT_MODULES_FILE}"

# Step 4: Parse every `#print axioms` line and enforce strict standard-axiom allowlist.
AUDIT_SUMMARY="$(
python3 - "${EXPECTED_AXIOMS_FILE}" "${AXIOM_LOG}" <<'PY'
import pathlib
import re
import sys

expected_lines = [
    line.strip().split("\t")
    for line in pathlib.Path(sys.argv[1]).read_text(encoding="utf-8").splitlines()
    if line.strip()
]
log_text = pathlib.Path(sys.argv[2]).read_text(encoding="utf-8")

ALLOWED = {"propext", "Classical.choice", "Quot.sound"}

dep_re = re.compile(r"^'([^']+)' depends on axioms: \[([^\]]*)\]\s*$")
none_re = re.compile(r"^'([^']+)' does not depend on any axioms\s*$")

observed = {}
bad_axioms = []

for line in log_text.splitlines():
    s = line.strip()
    m_dep = dep_re.match(s)
    if m_dep:
        name = m_dep.group(1)
        raw = m_dep.group(2).strip()
        axioms = [a.strip() for a in raw.split(",") if a.strip()]
        observed[name] = axioms
        disallowed = [a for a in axioms if a not in ALLOWED]
        if disallowed:
            bad_axioms.append((name, disallowed))
        continue
    m_none = none_re.match(s)
    if m_none:
        observed[m_none.group(1)] = []

if bad_axioms:
    print("ERROR: Disallowed axioms detected:", file=sys.stderr)
    for name, dis in bad_axioms:
        print(f"  - {name}: {dis}", file=sys.stderr)
    sys.exit(5)

missing = []
for rel, ident in expected_lines:
    matched = any(
        obs == ident or obs.endswith("." + ident)
        for obs in observed
    )
    if not matched:
        missing.append(f"{rel}: {ident}")

if missing:
    print("ERROR: Expected `#print axioms` output missing for:", file=sys.stderr)
    for item in missing:
        print(f"  - {item}", file=sys.stderr)
    sys.exit(6)

print(len(observed))
PY
)"
AUDIT_RC=$?
if [[ "${AUDIT_RC}" -ne 0 ]]; then
  fail "axiom_audit_failed" "${AUDIT_RC}"
fi

echo "LEAN_CHECK_PASS toolchain=${TOOLCHAIN} project=${PROJECT_DIR} axioms_checked=${AUDIT_SUMMARY} exit_code=0"
exit 0
