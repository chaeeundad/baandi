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
    def test_structural_evidence_changes_invalidate_generated_results(self):
        m = load()
        baseline = m.structural_results()
        original = m.load
        changes = [
            ("symmetry/evidence.json", lambda d: d["spin"].update(delta_energy_ha=0.123)),
            ("symmetry/mac-k444s.json", lambda d: [r.update(wall_s=123) for r in d["cases"]["crystal"]]),
            ("forces/manifest.json", lambda d: d["observations"].update(si_uncertainty_sum_ha_per_bohr=0.123)),
            ("forces/cost.json", lambda d: d[0]["forces"].update(elapsed_seconds=123)),
            ("relaxation/manifest.json", lambda d: d["ase_comparison"].update(energy_error_ha=0.123)),
            ("relaxation/mgo-curvature.json", lambda d: d["curvatures"][-1].update(directional_curvature_ha_bohr2=0.123)),
            ("stress/validation.json", lambda d: d["real_materials"]["Si-a10.6"].update(uncertainty_ha_bohr3=0.123)),
            ("stress/cost.json", lambda d: d["cases"]["si8"]["runs"][1]["stress"].update(elapsed_seconds=123)),
        ]
        for relative, mutate in changes:
            def changed(path):
                data = original(path)
                if os.path.relpath(path, m.REF) == relative:
                    mutate(data)
                return data
            with self.subTest(source=relative), patch.object(m, "load", side_effect=changed):
                self.assertNotEqual(baseline, m.structural_results())
        original_text = m.evidence_text
        for path, old, new in (
            ("tests/reference/stress/physics-terms.log", "2.3386910873168798e-8", "0.123"),
            ("tests/reference/relaxation/logs/cli-all.log", "0.2499999999839776", "0.123"),
        ):
            def changed_text(target):
                text = original_text(target)
                return text.replace(old, new) if target == path else text
            with self.subTest(source=path), patch.object(m, "evidence_text", side_effect=changed_text):
                self.assertNotEqual(baseline, m.structural_results())

    def test_structural_scope_and_sources(self):
        import re
        m = load()
        text = "\n".join(m.structural_results())
        for phrase in ("32 time-reversal representatives to 10", "inconclusive",
                       "DFT+U forces are unsupported", "standard profile is unsupported",
                       "CLI spin-polarized or finite-temperature relaxation is unsupported",
                       "not an error bound", "Second scalar energy reimplementation",
                       "byte-identical", "Cutoff convergence is not claimed",
                       "not full-tensor relaxation", "no comparison with experiment or other codes"):
            self.assertIn(phrase, text)
        self.assertEqual(text.count("comparison with QE pending"), 5)
        for target in re.findall(r"\]\(([^)]+)\)", text):
            self.assertTrue(Path(ROOT, target).is_file(), target)
        self.assertEqual(len(m.stress_term_errors()), 7)

    def test_cell_evidence_changes_invalidate_generated_results(self):
        m = load()
        original = m.load
        baseline = m.cell_results()
        changes = [
            ("validation.json", lambda d: d['synthetic_scans'][0].update(relative_error=0.123)),
            ("validation.json", lambda d: d['tetragonal'].update(relative_difference=0.123)),
            ("validation.json", lambda d: d['combined_scf_directional_derivative'].update(relative_error=0.123)),
            ("validation.json", lambda d: d.update(restart="changed restart evidence")),
            ("real-followup/real-validation.json", lambda d: d['cases'][0].update(lattice_bohr=12.3)),
            ("real-followup/real-validation.json", lambda d: d['cases'][0]['rebuilt_final'].update(pressure_ha_bohr3=0.123)),
            ("real-followup/real-validation.json", lambda d: d['cutoff_comparisons'][0].update(higher_minus_lower_cutoff_bohr=0.123)),
            ("real-followup/real-validation.json", lambda d: d['start_comparisons'][0].update(expanded_minus_compressed_bohr=0.123)),
            ("rebuild-followup/validation.json", lambda d: d['cases'][0]['stages'][-1].update(lattice_bohr=12.3)),
            ("rebuild-followup/validation.json", lambda d: d['cases'][0]['stages'][-1].update(rebuilt_pressure_ha_bohr3=0.123)),
            ("rebuild-followup/validation.json", lambda d: d['comparisons'][0].update(final_difference_bohr=0.123)),
        ]
        for relative, mutate in changes:
            def changed(path):
                data = original(path)
                if os.path.relpath(path, m.REF) == "relaxation/cell/" + relative:
                    mutate(data)
                return data
            with self.subTest(source=relative), patch.object(m, "load", side_effect=changed):
                self.assertNotEqual(baseline, m.cell_results())

    def test_cell_values_and_protocol_limits_are_visible(self):
        m = load()
        text = "\n".join(m.cell_results())
        for phrase in ("9.24e-07", "1.46e-07", "5.85e-10", "parabola fit",
                       "not electronic exact continuation", "5.435e-05", "2.223e-06",
                       "1.541e-04", "1.307e-04", "1.257e-08", "2.919e-09",
                       "6.825e-08", "4.036e-08", "9.212e-11", "1e-07 Ha/Bohr^3",
                       "frozen-only MgO start comparisons failed", "different effective cutoffs"):
            self.assertIn(phrase, text)
        status = m.load(os.path.join(ROOT, "docs/feature-status.json"))['items']
        self.assertEqual(next(r['status'] for r in status if r['feature'] == 'Cell relaxation'), 'done')
        self.assertEqual(next(r['status'] for r in status if r['feature'] == 'Molecular dynamics and phonons'), 'planned')

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
                          ("relaxation-forces.svg", m.fig_relaxation()),
                          ("performance.svg", m.fig_perf(m.perf_data())), ("thread-scaling.svg", m.fig_threads(m.threads_data()))):
            self.assertEqual(Path(ROOT, "docs", "img", name).read_text(), svg, name + " is stale")
        self.assertTrue(all(k["ok"] for k in kruns))


if __name__ == "__main__":
    unittest.main()
