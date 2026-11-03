import importlib.util, os, unittest

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
        text = open(os.path.join(ROOT, "README.md")).read()
        self.assertIn(m.BEGIN, text)
        got = text[text.index(m.BEGIN): text.index(m.END) + len(m.END)]
        self.assertEqual(got, expected, "README validation block is stale: run python3 -I scripts/make_readme_figures.py")

    def test_figures_match_the_evidence_and_all_runs_pass(self):
        m = load()
        rows, kruns = m.cutoff_data(), m.kpoint_data()
        g = next(r for r in rows if r["ecut"] == 30.0)
        gamma = {"e": abs(g["ours"] - g["qe"]) / 2, "eig": g["eig"], "iters": g["iters"]}
        for name, svg in (("cutoff-scan.svg", m.fig_cutoff(rows)), ("kpoint-comparison.svg", m.fig_kpoints(gamma, kruns)),
                          ("performance.svg", m.fig_perf(m.perf_data())), ("thread-scaling.svg", m.fig_threads(m.threads_data()))):
            self.assertEqual(open(os.path.join(ROOT, "docs", "img", name)).read(), svg, name + " is stale")
        self.assertTrue(all(k["ok"] for k in kruns))


if __name__ == "__main__":
    unittest.main()
