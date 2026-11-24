"""Standard-library-only controller fixture; launched with python3 -I."""
import json
import sys

mode = sys.argv[1] if len(sys.argv) > 1 else "echo"
for line in sys.stdin:
    row = json.loads(line)
    assert row["schema"] == "exp-scf-trace-1"
    assert row["phase"] == "observation"
    assert row["beta_used"] is None
    if mode == "eof":
        break
    if mode == "garbage":
        print("not-json", flush=True)
    elif mode == "range":
        print('{"beta":1.5001,"reset_history":true}', flush=True)
    elif mode == "missing":
        print('{"beta":0.7}', flush=True)
    elif mode == "reset":
        print(json.dumps({"beta": 1.2 if row["iteration"] == 2 else 0.7,
                          "reset_history": row["iteration"] == 2}), flush=True)
    else:
        print('{"beta":0.7,"reset_history":false}', flush=True)
