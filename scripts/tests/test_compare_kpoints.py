import copy, json, os, subprocess, sys, tempfile, unittest

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
SCRIPT = os.path.join(ROOT, "scripts", "compare_kpoints_qe.py")
REF = os.path.join(ROOT, "tests", "reference", "kpoints")


def run(*args):
    p = subprocess.run([sys.executable, "-I", SCRIPT, *args], capture_output=True, text=True)
    return p.returncode, p.stdout


class RealEvidence(unittest.TestCase):
    def test_full_qe_outputs(self):
        for t in ("k222_c", "k222s_c"):
            code, out = run("--ours", f"{REF}/{t}/ours_results.json", "--qe-out", f"{REF}/{t}/qe.out")
            self.assertEqual(code, 0, out)
            res = json.loads(out[: out.rindex("}") + 1])
            print(t, "eig diff Ha", res["max_eigenvalue_diff_ha"], "energy diff Ha/atom", res["total_energy_diff_ha_per_atom"])
            self.assertLess(res["max_eigenvalue_diff_ha"], 1e-5)
            self.assertLess(abs(res["total_energy_diff_ha_per_atom"]), 2e-6)

    def test_relayed_summary(self):
        code, out = run("--ours", f"{REF}/k444s_c/ours_results.json", "--qe-summary", f"{REF}/k444s_c/qe_summary.json")
        self.assertEqual(code, 0, out)


class Failures(unittest.TestCase):
    def setUp(self):
        self.d = tempfile.TemporaryDirectory()
        self.addCleanup(self.d.cleanup)

    def write(self, name, obj):
        p = os.path.join(self.d.name, name)
        json.dump(obj, open(p, "w"))
        return p

    def test_perturbed_eigenvalue_exits_5(self):
        ours = json.load(open(f"{REF}/k444s_c/ours_results.json"))
        ours["kpoints"]["points"][3]["eigenvalues"]["ev"][2] += 0.5  # 18 mHa, above the 1e-3 Ha tolerance
        code, _ = run("--ours", self.write("o.json", ours), "--qe-summary", f"{REF}/k444s_c/qe_summary.json")
        self.assertEqual(code, 5)

    def test_missing_k_point_exits_3(self):
        qe = json.load(open(f"{REF}/k444s_c/qe_summary.json"))
        qe["kpoints_cryst"][0][0] = 0.0625  # no ours k matches
        code, out = run("--ours", f"{REF}/k444s_c/ours_results.json", "--qe-summary", self.write("q.json", qe))
        self.assertEqual(code, 3, out)

    def test_count_mismatch_exits_3(self):
        ours = json.load(open(f"{REF}/k222_c/ours_results.json"))
        code, out = run("--ours", self.write("o.json", ours), "--qe-summary", f"{REF}/k444s_c/qe_summary.json")
        self.assertEqual(code, 3, out)

    def test_nondefault_tolerance_is_flagged(self):
        code, out = run("--ours", f"{REF}/k222_c/ours_results.json", "--qe-out", f"{REF}/k222_c/qe.out", "--eig-tol", "1e-2")
        self.assertIn("NON-DEFAULT", out)


if __name__ == "__main__":
    unittest.main()
