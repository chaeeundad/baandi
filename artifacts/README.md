# artifacts/

Human-facing status views. Source of truth for the status dashboard.

- `milestones.json`, `issues.json`, `tests.json`: dashboard datasets. Update these after every verified change.
- `dashboard.html`: page source of the dashboard (published as a Claude Artifact).

Rules: record only results that were actually run in this environment. A gate counts as done only after Claude re-runs it.
Dashboard: https://claude.ai/artifact/CEv2B7eMBsafD6f4JcKPvQ (private to the owner unless shared)
