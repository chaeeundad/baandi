#!/usr/bin/env python3
"""Prepare a tutorial input without editing the checked-in template (Python 3.9+)."""
import argparse
import json
import os
from pathlib import Path
import re
import shutil


def main():
    examples = Path(__file__).resolve().parent
    workflows = sorted(p.parent.name for p in examples.glob("*/input.toml"))
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("workflow", choices=workflows)
    args = parser.parse_args()
    if not os.environ.get("PSEUDO_DIR"):
        parser.error("set PSEUDO_DIR to the directory containing your external UPF files")
    pseudo = Path(os.environ["PSEUDO_DIR"]).expanduser().resolve(strict=True)
    source = (examples / args.workflow / "input.toml").read_text()

    def substitute(match):
        path = pseudo / match[1]
        if not path.is_file():
            parser.error(f"missing pseudopotential: {path}")
        return json.dumps(str(path), ensure_ascii=False)

    prepared = re.sub(r'"/absolute/path/to/pseudo/([^"/]+)"', substitute, source)
    destination = examples.parent / "target/tmp/examples" / args.workflow / "input.toml"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(prepared)
    # Structure files named by `[structure] file` resolve against the input
    # file, so the other files of the workflow are copied next to it.
    for extra in sorted((examples / args.workflow).iterdir()):
        if extra.is_file() and extra.name not in ("input.toml", "README.md"):
            shutil.copyfile(extra, destination.parent / extra.name)
    print(destination)


if __name__ == "__main__":
    main()
