#!/usr/bin/env python3
"""Re-run implemented P11L evidence without modifying historical reference bundles.

Use --output for a new evidence directory; default output is a temporary directory.
Incomplete coverage is reported explicitly, never counted as a completed case.
"""
import argparse
import datetime
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import platform
import re
import subprocess
import sys
import tempfile
import time
import unittest

ROOT = Path(__file__).resolve().parents[1]
REFERENCE = ROOT / "docs/case-studies/data/recorded"


def write(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n")


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def preflight_cases(out, run):
    """Review-only paired Classic runs; proposed advisory policy is never applied."""
    spec = importlib.util.spec_from_file_location("synthetic", ROOT / "scripts/gen_synthetic_upf.py")
    generator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(generator)
    cli = str(ROOT / "target/debug/dft-cli")
    results = []
    for name, symbol, z, channels, cell, positions in (
        ("cubic_h2", "H", 1, (0,), [[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]], [[0., 0., 0.], [.5, .5, .5]]),
        ("skew_he", "He", 2, (0, 1), [[4., 0., 0.], [.5, 4., 0.], [0., .5, 5.]], [[0., 0., 0.]]),
    ):
        case = out / name
        case.mkdir()
        with tempfile.TemporaryDirectory(prefix="p11l-pseudo-") as folder:
            pp = generator.generate(z, channels, radius=1.4, points=401).replace('element="Xx"', f'element="{symbol}"')
            pseudo = Path(folder) / "synthetic.upf"
            pseudo.write_text(pp)
            text = f'[system]\nname="{name}"\nlattice_unit="bohr"\nlattice_vectors={json.dumps(cell)}\n'
            for position in positions:
                text += f'[[system.atoms]]\nspecies="{symbol}"\nfractional={json.dumps(position)}\n'
            text += f'[pseudopotentials]\n{symbol}={json.dumps(str(pseudo))}\n[basis]\necutwfc_ry=8.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional="lda-pz"\n[scf]\nnbands=4\n[preflight]\nminimum_distance_bohr=100.0\n'
            input_path = case / "input.toml"
            input_path.write_text(text)
            start = time.perf_counter()
            run(name + "-classic", [cli, "run", str(input_path), "--output-dir", str(case / "classic")])
            classic_seconds = time.perf_counter() - start
            start = time.perf_counter()
            findings = run(name + "-findings", [cli, "check", str(input_path), "--json"])
            (case / "findings.json").write_text(findings)
            candidate = text.replace("minimum_distance_bohr=100.0", "minimum_distance_bohr=0.1")
            old, new = text.splitlines(), candidate.splitlines()
            diff = f'--- a/input.toml\n+++ b/input.toml\n@@ -1,{len(old)} +1,{len(new)} @@\n'
            diff += ''.join('-' + line + '\n' for line in old) + ''.join('+' + line + '\n' for line in new)
            sha = lambda value: hashlib.sha256(value.encode()).hexdigest()
            proposal = dict(schema_version=1, input_sha256=sha(text), findings_sha256=sha(findings),
                pseudopotentials={symbol: sha(pp)}, candidate_pseudopotentials={symbol: sha(pp)},
                advisory="Review the deliberately oversized distance warning threshold; no physical convergence claim.",
                agent="P11L scripted optional adviser; no LLM", diff=diff,
                edits=[dict(key_path=["preflight", "minimum_distance_bohr"], before=100., after=.1,
                            reason="Explicit fixture review policy only", finding_codes=["PF_MIN_DISTANCE"])],
                candidate_input_sha256=sha(candidate))
            proposal_path = case / "proposal.json"
            write(proposal_path, proposal)
            run(name + "-preview", [cli, "preflight-review", str(input_path), str(proposal_path)])
            for decision in ("approved", "rejected"):
                run(name + "-" + decision, [cli, "preflight-review", str(input_path), str(proposal_path),
                    "--decision", decision, "--actor", "P11L-scripted-operator",
                    "--reason", "Frozen review-only fixture; never apply candidate",
                    "--output-dir", str(case / decision)])
            review_seconds = time.perf_counter() - start
            stale = dict(proposal, input_sha256="0" * 64)
            write(case / "stale-proposal.json", stale)
            run(name + "-stale", [cli, "preflight-review", str(input_path), str(case / "stale-proposal.json")], expected=1)
            if input_path.read_text() != text:
                raise RuntimeError("Preflight review changed original input")
            start = time.perf_counter()
            run(name + "-reviewed", [cli, "run", str(input_path), "--output-dir", str(case / "reviewed")])
            reviewed_seconds = time.perf_counter() - start
            a = json.loads((case / "classic/results.json").read_text())
            b = json.loads((case / "reviewed/results.json").read_text())
            energy = a["total_energy"]["ha"]
            if not a["convergence"]["converged"] or not b["convergence"]["converged"]:
                raise RuntimeError("Preflight paired run did not converge")
            if energy.hex() != b["total_energy"]["ha"].hex() or len(a["convergence"]["history"]) != len(b["convergence"]["history"]):
                raise RuntimeError("Review-only path changed Classic energy")
            results.append(dict(system=name, upf_sha256=sha(pp), energy_hartree=energy,
                iterations=len(a["convergence"]["history"]),
                energy_difference_hartree=b["total_energy"]["ha"] - energy,
                classic_wall_seconds=classic_seconds, reviewed_run_wall_seconds=reviewed_seconds,
                review_wall_seconds=review_seconds,
                workflow_overhead_seconds=review_seconds + reviewed_seconds - classic_seconds,
                timing_scope="single sequential CLI process runs, includes I/O; review includes preview and both decisions; excludes stale negative control",
                reference="same-build Classic, identical input bytes; candidate never applied"))
    write(out / "preflight-results.json", results)


def python_worker(out):
    spec = importlib.util.spec_from_file_location(
        "binding_cases", ROOT / "python/baandi/tests/test_binding.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    suite = unittest.defaultTestLoader.loadTestsFromModule(module)
    result = unittest.TextTestRunner(verbosity=2).run(suite)
    if not result.wasSuccessful() or result.skipped or result.testsRun != 5:
        raise RuntimeError("All five installed-wheel tests must execute; skips fail")
    import ase
    import baandi
    from ase import Atoms
    from baandi.ase import Baandi
    fixtures = []
    timings = []
    for name, symbol, z, channels, cell, positions in module.CASES:
        pp = module.potential(symbol, z, channels)
        fixtures.append(dict(system=name, symbol=symbol, valence=z,
                             projectors=channels, native_cell_bohr=cell,
                             positions_fractional=positions,
                             upf_sha256=hashlib.sha256(pp.encode()).hexdigest()))
        # This separate paired run measures end-to-end paths after the regression suite.
        case = out / ("ase-" + name)
        case.mkdir()
        ase_cell = [[2., 0., 0.], [0., 2., 0.], [0., 0., 2.]] if name == "cubic_h2" else [[2., 0., 0.], [.25, 2., 0.], [0., .25, 2.5]]
        atoms = Atoms([symbol] * len(positions), scaled_positions=positions, cell=ase_cell, pbc=True)
        fractions = atoms.get_scaled_positions(wrap=False).tolist()
        with tempfile.TemporaryDirectory(prefix="p11l-ase-") as folder:
            pseudo = Path(folder) / "synthetic.upf"
            pseudo.write_text(pp)
            text = f'[system]\nname="{name}"\nlattice_unit="angstrom"\nlattice_vectors={json.dumps(ase_cell)}\n'
            for pos in fractions:
                text += f'[[system.atoms]]\nspecies="{symbol}"\nfractional={json.dumps(pos)}\n'
            text += f'[pseudopotentials]\n{symbol}={json.dumps(str(pseudo))}\n[basis]\necutwfc_ry=8.0\nkpoints=[[0.0,0.0,0.0]]\n[xc]\nfunctional="lda-pz"\n[scf]\nnbands=4\n'
            (case / "input.toml").write_text(text)
            command = [os.environ["BAANDI_CLI"], "run", str(case / "input.toml"), "--output-dir", str(case / "classic")]
            start = time.perf_counter()
            completed = subprocess.run(command, capture_output=True, text=True, check=True)
            cli_seconds = time.perf_counter() - start
            (case / "classic.log").write_text(completed.stdout + completed.stderr)
            classic = json.loads((case / "classic/results.json").read_text())
            start = time.perf_counter()
            atoms.calc = Baandi(pseudopotentials={symbol: pp}, options=module.settings())
            ev = atoms.get_potential_energy()
            ase_seconds = time.perf_counter() - start
            native = atoms.calc.session.result()
            if not classic["convergence"]["converged"] or not native.converged or ev.hex() != classic["total_energy"]["ev"].hex():
                raise RuntimeError("Timed ASE/Classic mismatch")
            timings.append(dict(system=name, command=command, exit_code=completed.returncode,
                energy_hartree=native.energy_hartree, energy_ev=ev, iterations=native.iterations,
                classic_iterations=len(classic["convergence"]["history"]),
                energy_difference_hartree=native.energy_hartree - classic["total_energy"]["ha"],
                classic_wall_seconds=cli_seconds, ase_wall_seconds=ase_seconds,
                workflow_overhead_seconds=ase_seconds-cli_seconds,
                timing_scope="one sequential pair; CLI process includes file I/O, ASE in-process includes calculator construction; no isolated adapter-overhead claim"))
    write(out / "python-timings.json", timings)
    options = module.settings()
    write(out / "python-setup.json", dict(
        python=sys.version, ase=ase.__version__, binding=baandi.__file__,
        fixtures=fixtures, options={k: getattr(options, k) for k in (
            "cutoff_hartree", "bands", "extra_bands", "max_iterations",
            "density_tolerance", "energy_tolerance", "eigen_tolerance", "seed",
            "mixing_beta", "mixing_history")}))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", choices=["all", "self-healing", "mcp", "preflight", "python", "python-worker"])
    parser.add_argument("--output", type=Path)
    parser.add_argument("--python", default=sys.executable,
                        help="Interpreter with ASE and maturin installed (no network install)")
    parser.add_argument("--verify", action="store_true", help="Compare to recorded numerical evidence")
    args = parser.parse_args()
    if args.case == "python-worker":
        python_worker(args.output)
        return
    out = args.output.resolve() if args.output else Path(tempfile.mkdtemp(prefix="p11l-"))
    if args.output:
        out.mkdir(parents=True, exist_ok=False)
    print("Evidence:", out, flush=True)
    env = dict(os.environ, PATH="/opt/homebrew/opt/rustup/bin:" + os.environ["PATH"],
               CARGO_NET_OFFLINE="true")
    report = dict(schema_version=1, task_id="P11L-CASE-STUDIES",
                  started_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                  platform=platform.platform(), commands=[], status="running",
                  coverage="partial: no MLIP or low-level coupling implementation; see case list")
    def run(label, command, cwd=ROOT, extra=None, expected=0):
        start = time.perf_counter()
        p = subprocess.run(command, cwd=cwd, env={**env, **(extra or {})},
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
        (out / (label + ".log")).write_text(p.stdout)
        report["commands"].append(dict(label=label, command=command, cwd=str(cwd),
                                       environment=extra or {}, exit_code=p.returncode,
                                       wall_seconds=time.perf_counter() - start))
        write(out / "provenance.json", report)
        if p.returncode != expected:
            raise RuntimeError(f"{label} failed: see {out / (label + '.log')}")
        return p.stdout

    try:
        report["revision"] = run("revision", ["git", "rev-parse", "HEAD"]).strip()
        run("worktree", ["git", "status", "--short"])
        run("rustc", ["rustc", "-Vv"])
        run("cargo", ["cargo", "-V"])
        sources = [ROOT / "Cargo.lock", ROOT / "rust-toolchain.toml", Path(__file__)]
        for directory in ("crates", "python", "tests/reference/runtime"):
            sources.extend(p for p in (ROOT / directory).rglob("*")
                           if p.is_file() and p.suffix in (".rs", ".toml", ".py", ".json"))
        sources.append(ROOT / "scripts/gen_synthetic_upf.py")
        report["sources_sha256"] = {str(p.relative_to(ROOT)): digest(p) for p in sorted(set(sources))}
        selected = ["self-healing", "mcp", "preflight", "python"] if args.case == "all" else [args.case]
        if "self-healing" in selected:
            run("self-healing", ["bash", "scripts/demo_self_healing.sh"],
                extra={"DFT_HEAL_OUTPUT": str(out / "self-healing")})
            run("self-healing-rejections", ["cargo", "test", "--offline", "--locked", "--release",
                "-p", "dft-cli", "--test", "self_healing", "--", "--skip", "three_arm_demo"])
            if args.verify:
                current = json.loads((out / "self-healing/results.json").read_text())
                old = json.loads((REFERENCE / "self-healing/results.json").read_text())
                if len(current["runs"]) != len(old["runs"]):
                    raise RuntimeError("Changed self-healing run count")
                for a, b in zip(current["runs"], old["runs"]):
                    for key in ("state", "iterations", "converged", "policy_changes", "rejected"):
                        if a[key] != b[key]:
                            raise RuntimeError(f"Changed self-healing {key}")
                    if abs(a["energy_hartree"] - b["energy_hartree"]) > 1e-10:
                        raise RuntimeError("Self-healing energy differs by >1e-10 Ha")
        for case, package, test, count in (("mcp", "baandi-mcp", "server", 9),
                                           ("preflight", "dft-cli", "preflight_findings", 11)):
            if case in selected:
                log = run(case, ["cargo", "test", "--offline", "--locked", "--release",
                    "-p", package, "--test", test, "--", "--nocapture", "--test-threads=1"])
                if not re.search(rf"test result: ok\. {count} passed; 0 failed; 0 ignored", log):
                    raise RuntimeError(f"Expected all {count} {case} checks")
        if "preflight" in selected:
            run("preflight-cli-build", ["cargo", "build", "--offline", "--locked", "-p", "dft-cli"])
            preflight_cases(out, run)
            if args.verify:
                current = json.loads((out / "preflight-results.json").read_text())
                old = json.loads((REFERENCE / "preflight-results.json").read_text())
                if len(current) != len(old):
                    raise RuntimeError("Changed preflight system count")
                for a, b in zip(current, old):
                    if a["system"] != b["system"] or abs(a["energy_hartree"] - b["energy_hartree"]) > 1e-10:
                        raise RuntimeError("Preflight reference mismatch (>1e-10 Ha)")
        if "python" in selected:
            run("cli-build", ["cargo", "build", "--offline", "--locked", "-p", "dft-cli"])
            wheels = out / "wheels"
            run("wheel-build", [args.python, "-m", "maturin", "build", "--offline", "--locked",
                "--out", str(wheels)], cwd=ROOT / "python/baandi", extra={"PYO3_PYTHON": args.python})
            wheel, = wheels.glob("*.whl")
            report["wheel_sha256"] = digest(wheel)
            install = out / "installed"
            run("wheel-install", [args.python, "-m", "pip", "install", "--no-deps", "--no-index",
                                  "--target", str(install), str(wheel)])
            run("python", [args.python, str(Path(__file__).resolve()), "python-worker", "--output", str(out)],
                extra={"PYTHONPATH": str(install), "BAANDI_CLI": str(ROOT / "target/debug/dft-cli"),
                       "BAANDI_EVIDENCE": str(out / "ase.json"),
                       "BAANDI_NATIVE_EVIDENCE": str(out / "native.json")})
            if args.verify:
                for name in ("ase.json", "native.json"):
                    current = json.loads((out / name).read_text())
                    old = json.loads((REFERENCE / name).read_text())
                    if len(current) != len(old):
                        raise RuntimeError("Changed Python system count")
                    for a, b in zip(current, old):
                        if a["system"] != b["system"] or a["iterations"] != b["iterations"] or abs(a["energy_ha"] - b["energy_ha"]) > 1e-10:
                            raise RuntimeError("Python reference mismatch (>1e-10 Ha or iteration change)")
        report["status"] = "passed_selected_checks"
    except Exception as exc:
        report["status"] = "failed"
        report["error"] = str(exc)
        raise
    finally:
        report["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat()
        write(out / "provenance.json", report)


if __name__ == "__main__":
    main()
