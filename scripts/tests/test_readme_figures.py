import importlib.util, os, unittest
from pathlib import Path
from unittest.mock import patch

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


def load():
    spec = importlib.util.spec_from_file_location("mrf", os.path.join(ROOT, "scripts", "make_readme_figures.py"))
    m = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(m)
    return m


class ReadmeIsCurrent(unittest.TestCase):
    def test_readme_validation_block_matches_the_evidence(self):
        m = load()
        rows, kruns = m.cutoff_data(), m.kpoint_data()
        g = next(r for r in rows if r["ecut"] == 30.0)
        gamma = {"e": abs(g["ours"] - g["qe"]) / 2, "eig": g["eig"], "iters": g["iters"]}
        expected = m.block(rows, kruns, gamma, m.perf_data())
        text = Path(ROOT, "README.md").read_text()
        self.assertIn(m.BEGIN, text)
        got = text[text.index(m.BEGIN): text.index(m.END) + len(m.END)]
        self.assertEqual(got, expected, "README validation block is stale: run python3 -I scripts/make_readme_figures.py")

    def test_new_evidence_changes_make_the_readme_stale(self):
        m = load()
        original_load = m.load
        baseline = "\n".join(m.extended_results())
        changes = [
            ("crates/dft-engine/tests/collinear-evidence.json",
             lambda d: d["matrix"][0].update(energy_difference_ha=0.123)),
            ("tests/reference/datasets/p13c-status.json",
             lambda d: d.update(gaussian_core_max_relative_error=0.123)),
            ("tests/reference/materials/results.json",
             lambda d: d["cases"][0].update(total_energy_ha=-123.0)),
            ("tests/reference/spectral/summary.json",
             lambda d: next(r for r in d if r["case"] == "si-20").update(occupied_max_error_ha=0.123)),
            ("tests/reference/perf/p17a/final-report.json",
             lambda d: next(r for r in d["rows"] if r["label"] == "strict").update(davidson_steps=123456)),
            ("tests/reference/occupations/summary.json",
             lambda d: d["synthetic"].update(free_energy_ha=-123.0)),
        ]
        for relative, mutate in changes:
            def changed_load(path):
                data = original_load(path)
                if os.path.relpath(path, ROOT) == relative:
                    mutate(data)
                return data
            with self.subTest(source=relative), patch.object(m, "load", side_effect=changed_load):
                self.assertNotEqual(baseline, "\n".join(m.extended_results()))

    def test_new_evidence_references_and_limits_are_visible(self):
        m = load()
        text = "\n".join(m.extended_results())
        for warning in ("comparison with QE pending for Si and MgO",
                        "Standard results are not certified", "relayed QE",
                        "not a universal completeness proof", "unqualified",
                        "finite-temperature checkpoint persistence is unavailable"):
            self.assertIn(warning, text)
        # Every generated evidence link must resolve, including the spin exception.
        import re
        for target in re.findall(r"\]\(([^)]+)\)", text):
            self.assertTrue(Path(ROOT, target).is_file(), target)
        rows = m.standard_rows()
        failed = [(a, b) for a, b in rows if b["status"] != "ok"]
        self.assertTrue(failed, "Do not hide the failed stress cases")
        for strict, standard in failed:
            row = next(line for line in text.splitlines() if line.startswith("| " + strict["case"] + " |"))
            self.assertIn("unavailable", row)
            self.assertNotIn("0.00e+00", row)

    def test_status_block_is_current(self):
        m = load()
        text = Path(ROOT, "README.md").read_text()
        got = text[text.index(m.STATUS_BEGIN): text.index(m.STATUS_END) + len(m.STATUS_END)]
        self.assertEqual(got, m.status_block(), "README status block is stale: run python3 -I scripts/make_readme_figures.py")

    def test_figures_match_the_evidence_and_all_runs_pass(self):
        m = load()
        rows, kruns = m.cutoff_data(), m.kpoint_data()
        g = next(r for r in rows if r["ecut"] == 30.0)
        gamma = {"e": abs(g["ours"] - g["qe"]) / 2, "eig": g["eig"], "iters": g["iters"]}
        for name, svg in (("cutoff-scan.svg", m.fig_cutoff(rows)), ("kpoint-comparison.svg", m.fig_kpoints(gamma, kruns)),
                          ("performance.svg", m.fig_perf(m.perf_data())), ("thread-scaling.svg", m.fig_threads(m.threads_data()))):
            self.assertEqual(Path(ROOT, "docs", "img", name).read_text(), svg, name + " is stale")
        self.assertTrue(all(k["ok"] for k in kruns))


if __name__ == "__main__":
    unittest.main()
