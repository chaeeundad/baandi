# BAANDI project wiki (synthesised knowledge, for people and agents)

What this is: the distilled layer over the raw records. Raw records stay where they are and remain the evidence: `docs/performance/attempts-ledger.md` (every optimisation tried, with numbers), `docs/decisions/` (ADRs), `tests/reference/perf/*relayed*.md` (dev08 measurements), `docs/research/` (QE and literature notes), `artifacts/*.json` (state). The wiki pages below say what is true now, why, and where the evidence is. Read order for a performance task: [state-of-play](state-of-play.md), [decision-table](decision-table.md), [lessons](lessons.md), then the raw ledger rows you need.

Pages
- [state-of-play.md](state-of-play.md): current numbers against QE, defaults, open work.
- [decision-table.md](decision-table.md): every optimisation idea with status (adopted, opt-in, rejected, open), the reason, and the ledger/ADR row.
- [qe-knowledge.md](qe-knowledge.md): what we learned about QE 7.6 (start, thresholds, timers, hard cells, the energy offset).
- [lessons.md](lessons.md): measurement and process lessons (what went wrong before, how to avoid it).
- [dev08-playbook.md](dev08-playbook.md): how to measure on dev08 and what each switch/environment variable does.
- [mpi.md](mpi.md): the MPI layer in one page.

Maintenance rules
1. Update the wiki in the same commit as a ledger row that changes a status or a headline number; keep one fact in one place and link to the evidence instead of copying tables.
2. A statement in the wiki carries a date and a commit or ledger id; a number without a source is not written.
3. When a page contradicts the ledger, the ledger wins and the page is fixed.
4. User decisions are quoted with their date (AGENTS.md scope decisions U1 to U9).
