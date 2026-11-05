#!/usr/bin/env python3
"""Link the existing P13B fixture generator after a release CLI build."""
import json
import os
from pathlib import Path
import subprocess
ROOT=Path(__file__).resolve().parents[3]
lib=max((ROOT/'target/release/deps').glob('libdft_physics-*.rlib'),key=lambda p:p.stat().st_mtime)
env=dict(os.environ,PATH='/opt/homebrew/opt/rustup/bin:'+os.environ['PATH'])
command=['rustc','--edition=2024',str(Path(__file__).with_name('spin_fixture.rs')),
         '-L',f'dependency={ROOT}/target/release/deps','--extern',f'dft_physics={lib}',
         '-o','/tmp/p15c-spin-fixture']
proc=subprocess.run(command,env=env)
Path(__file__).with_name('spin-generator-build.json').write_text(json.dumps(dict(command=command,exit_code=proc.returncode),indent=2)+'\n')
raise SystemExit(proc.returncode)
