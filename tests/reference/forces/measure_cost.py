"""Mac observation only; exact input, UPF hash, force pass counters and times."""
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile

root = Path(__file__).resolve().parents[3]
spec = importlib.util.spec_from_file_location("bench", root / "scripts/bench_scf.py")
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)
binary, upf, destination = map(Path, sys.argv[1:])
records = []
for case in ("si8", "k444s"):
    with tempfile.TemporaryDirectory() as folder:
        folder = Path(folder)
        text = bench.case_input(str(upf.resolve()), case) + "\n[output]\nforces=true\n"
        path = folder / "input.toml"
        path.write_text(text)
        command = [str(binary.resolve()), "run", str(path), "--output-dir", str(folder / "out")]
        subprocess.run(command, check=True, capture_output=True)
        result = json.loads((folder / "out/results.json").read_text())
        records.append({"case": case, "input": text, "upf_sha256": bench.sha256(upf),
                        "binary_sha256": bench.sha256(binary),
                        "forces": result["forces"], "timings": result.get("timings"),
                        "force_share_scf_plus_force": result["forces"]["elapsed_seconds"] / (result["forces"]["elapsed_seconds"] + result["timings"]["engine_total_s"]),
                        "scope": "one shared-host Mac observation; not an optimization acceptance benchmark"})
destination.write_text(json.dumps(records, indent=2) + "\n")
