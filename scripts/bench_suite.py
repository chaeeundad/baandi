#!/usr/bin/env python3
"""Reproducible generality/robustness suite (Python >=3.9, standard library).

Case files use a deliberately restricted flat TOML subset: JSON-compatible
values on one line, no tables or inline comments. This avoids a TOML dependency
on the Mac's Python 3.9. Generated CLI inputs are ordinary TOML.
Exit 0: completed (findings allowed in robustness mode); 1: setup error;
3: unexpected run failure; 5: baseline regression/incomparable input.
"""
import argparse
import copy
import datetime
import hashlib
import importlib.util
import json
import math
import os
from pathlib import Path
import platform
import random
import re
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
_spec = importlib.util.spec_from_file_location("suite_upf", ROOT / "scripts/gen_synthetic_upf.py")
_upf = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_upf)
ENERGY_TOL_HA = 1e-12
KEYS = set("id purpose system_class tier potential lattice atoms cutoff_ry kgrid kshift fft_grid nbands valence channels radius shape timeout_s expected_error deterministic max_iterations".split())


def digest(data):
    return hashlib.sha256(data).hexdigest()


def file_hash(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as f:
        for block in iter(lambda: f.read(1 << 20), b''):
            h.update(block)
    return h.hexdigest()


def load_case(path):
    case = {}
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.lstrip().startswith('#'):
            continue
        key, sep, value = line.partition('=')
        key = key.strip()
        if not sep or key not in KEYS or key in case:
            raise ValueError(f"{path}: duplicate/unknown key or invalid line: {line}")
        case[key] = json.loads(value, parse_constant=lambda s: (_ for _ in ()).throw(ValueError(s)))
    for key in ('id', 'purpose', 'system_class', 'tier', 'potential', 'lattice', 'atoms', 'cutoff_ry', 'timeout_s'):
        if key not in case:
            raise ValueError(f"{path}: missing {key}")
    if not re.fullmatch(r'[a-z0-9][a-z0-9_-]*', case['id']):
        raise ValueError("unsafe case id")
    if case['tier'] not in ('default', 'long') or case['potential'] not in ('si', 'synthetic'):
        raise ValueError("invalid tier/potential")
    if not isinstance(case['timeout_s'], (int, float)) or not 0 < case['timeout_s'] < math.inf:
        raise ValueError("timeout must be positive and finite")
    if len(case['lattice']) != 3 or not case['atoms']:
        raise ValueError("three lattice vectors and at least one atom required")
    for row in case['lattice'] + case['atoms']:
        if len(row) != 3 or any(type(x) not in (int, float) or not math.isfinite(x) for x in row):
            raise ValueError("coordinates must be finite three-vectors")
    return case


def synthetic(case):
    return _upf.generate(case.get('valence', 4), case.get('channels', [0, 0, 1]),
                         case.get('radius', 1.0), case.get('shape', 'gaussian'))


def variant(case, seed, index, displacement):
    c = copy.deepcopy(case)
    # Stable independent streams: selection/order changes do not alter a variant.
    rng = random.Random(digest(f"{case['id']}:{seed}:{index}".encode()))
    c['id'] += f'-variant-{index}'
    c.pop('expected_error', None)
    c['atoms'] = [[(x + rng.uniform(-displacement, displacement)) % 1.0 for x in site]
                  for site in c['atoms']]
    c['purpose'] += f'; fractional displacement +/-{displacement}, variant {index}, seed {seed}'
    return c


def case_input(case, upf, seed, threads):
    species = 'Si' if case['potential'] == 'si' else 'Xx'
    lines = ['[system]', f'name = {json.dumps(case["id"])}', 'lattice_unit = "bohr"',
             f'lattice_vectors = {json.dumps(case["lattice"])}']
    for site in case['atoms']:
        lines += ['[[system.atoms]]', f'species = "{species}"', f'fractional = {json.dumps(site)}']
    lines += ['[pseudopotentials]', f'{species} = {json.dumps(str(upf))}', '[basis]',
              f'ecutwfc_ry = {case["cutoff_ry"]}']
    if 'kgrid' in case:
        lines += [f'kgrid = {json.dumps(case["kgrid"])}', f'kshift = {json.dumps(case.get("kshift", [0,0,0]))}']
    else:
        lines += ['kpoints = [[0.0, 0.0, 0.0]]']
    if 'fft_grid' in case:
        lines += [f'fft_grid = {json.dumps(case["fft_grid"])}']
    lines += ['[xc]', 'functional = "lda-pz"', '[scf]', f'seed = {seed}', f'threads = {threads}',
              f'max_iterations = {case.get("max_iterations", 100)}', 'mixing = "pulay"', 'mixing_beta = 0.5',
              'energy_tolerance_ha = 1e-8', 'density_residual_tolerance = 1e-7']
    if 'nbands' in case:
        lines += [f'nbands = {case["nbands"]}']
    return '\n'.join(lines) + '\n'


def measurements(report):
    if type(report['convergence']['converged']) is not bool:
        raise ValueError('converged must be boolean')
    if not math.isfinite(report['total_energy']['ha']):
        raise ValueError('nonfinite total energy')
    history = report['convergence']['history']
    # Certificate work occurs in history; do not count last_iteration again.
    counts = []
    for step in history:
        counts.append(step.get('operator_applications'))
        cert = step.get('fixed_point_certificate')
        if cert:
            counts.append(cert.get('operator_applications'))
    points = report.get('kpoints', {}).get('points') or [report]
    gaps = []
    homos, lumos = [], []
    for point in points:
        ev = point.get('eigenvalues', {}).get('ha', [])
        occ = point.get('occupations', [])
        occupied = [v for v, n in zip(ev, occ) if n > 0]
        empty = [v for v, n in zip(ev, occ) if n == 0]
        if occupied and empty:
            homos.append(max(occupied)); lumos.append(min(empty))
            gaps.append(min(empty)-max(occupied))
    return dict(converged=report['convergence']['converged'], iterations=len(history),
                total_energy_ha=report['total_energy']['ha'], phases_s=report['timings']['phases'],
                engine_total_s=report['timings']['engine_total_s'],
                operator_applications=sum(counts) if counts and all(x is not None for x in counts) else None,
                fft_grid=report['basis']['fft_grid'], bands=report['bands'],
                tool=report.get('tool'), threading=report.get('threading'),
                eigenvalues_ha=[p.get('eigenvalues', {}).get('ha', []) for p in points],
                last_iteration_diagnostics={key: history[-1].get(key) for key in
                    ('density_residual', 'max_eigen_residual_ha', 'relative_charge_error',
                     'fixed_point_certificate')} if history else None,
                minimum_direct_gap_ha=min(gaps) if gaps else None,
                indirect_gap_ha=min(lumos)-max(homos) if gaps else None)


def run_case(binary, case, upf, directory, seed, threads, repeats, timeout):
    directory.mkdir(parents=True)
    if case['potential'] == 'synthetic':
        upf = directory / 'synthetic.UPF'
        upf.write_text(synthetic(case), encoding='utf-8')
    input_text = case_input(case, upf, seed, threads)
    (directory / 'input.toml').write_text(input_text, encoding='utf-8')
    # Portable reproducer embedded in JSON; real potential is never embedded/copied.
    portable = case_input(case, 'SI_UPF_PATH' if case['potential'] == 'si' else 'synthetic.UPF', seed, threads)
    fingerprint = digest((portable + file_hash(upf)).encode())
    record = dict(case=case, input_toml=portable, input_sha256=digest(input_text.encode()),
                  fingerprint=fingerprint, upf_sha256=file_hash(upf), runs=[])
    for i in range(repeats):
        output = directory / f'out-{i}'
        command = [str(binary), 'run', str(directory / 'input.toml'), '--output-dir', str(output)]
        start = time.perf_counter()
        run = dict(command=command, converged=False, iterations=None, total_energy_ha=None,
                   density_sha256=None, phases_s=None, operator_applications=None)
        try:
            p = subprocess.run(command, capture_output=True, text=True, timeout=timeout or case['timeout_s'])
            stdout, stderr = p.stdout, p.stderr
            run.update(exit_code=p.returncode, status='error', error=stderr.strip())
            result_path = output / 'results.json'
            if result_path.exists():
                try:
                    r = json.loads(result_path.read_text())
                    run.update(measurements(r))
                    run['density_sha256'] = file_hash(output / 'density.bin')
                    run['status'] = 'pass' if p.returncode == 0 and run['converged'] else 'nonconverged'
                    if run['status'] != 'pass' and not run['error']:
                        run['error'] = 'SCF did not converge'
                except (KeyError, TypeError, ValueError, OSError) as e:
                    run.update(status='invalid_output', converged=False, error=f'{type(e).__name__}: {e}')
            elif p.returncode == 0:
                run.update(status='invalid_output', error='exit 0 without results.json')
        except subprocess.TimeoutExpired as e:
            def decoded(s):
                return s.decode(errors='replace') if isinstance(s, bytes) else (s or '')
            stdout, stderr = decoded(e.stdout), decoded(e.stderr)
            run.update(exit_code=None, status='timeout', error=f'timeout after {timeout or case["timeout_s"]} s')
        run['wall_s'] = time.perf_counter() - start
        (directory / f'stdout-{i}.txt').write_text(stdout)
        (directory / f'stderr-{i}.txt').write_text(stderr)
        record['runs'].append(run)
    runs = record['runs']
    identities = {(r['status'], r['total_energy_ha'], r['density_sha256'], r['iterations']) for r in runs}
    record.update(status=runs[-1]['status'], converged=all(r['converged'] for r in runs),
                  repeat_identical=len(identities) == 1, wall_median_s=statistics.median(r['wall_s'] for r in runs),
                  wall_min_s=min(r['wall_s'] for r in runs), wall_max_s=max(r['wall_s'] for r in runs),
                  timeout_s=timeout or case['timeout_s'])
    if len(identities) != 1:
        record['status'] = 'nondeterministic'
    for key in ('total_energy_ha', 'iterations', 'density_sha256', 'operator_applications'):
        record[key] = runs[-1][key]
    expected = case.get('expected_error')
    record['expected_failure'] = bool(expected and all(r['status'] == 'error' and expected in r['error'] for r in runs))
    record['resolved_expected_failure'] = bool(expected and record['status'] == 'pass')
    return record


def compare(current, baseline, require_identical=False):
    findings = []
    for key, now in current.items():
        old = baseline.get(key)
        issues, improvements = [], []
        if old is None:
            issues.append('missing baseline case')
        elif now['fingerprint'] != old['fingerprint']:
            issues.append('incomparable input or UPF (includes seed and threads)')
        else:
            if old['status'] != 'pass' and now['status'] == 'pass':
                improvements.append('previous failure now passes')
            elif now['status'] != old['status'] or now['converged'] != old['converged']:
                issues.append('status/convergence changed')
            if old['status'] == now['status'] == 'pass':
                de = abs(now['total_energy_ha'] - old['total_energy_ha'])
                if now['case'].get('deterministic', True) and (not math.isfinite(de) or de > ENERGY_TOL_HA):
                    issues.append(f'energy changed by {de:.17g} Ha (>1e-12 Ha)')
                if now['iterations'] != old['iterations']:
                    issues.append('iteration count changed')
                if require_identical and now['density_sha256'] != old['density_sha256']:
                    issues.append('density hash changed')
            elif not improvements and {r['error'] for r in now['runs']} != {r['error'] for r in old['runs']}:
                issues.append('failure message changed')
            if not now['repeat_identical']:
                issues.append('repeats differ')
        findings.append(dict(id=key, regressions=issues, improvements=improvements))
    return findings


def markdown(data):
    lines = ['# Benchmark suite results', '', 'Synthetic results are non-physical. Timings are indicative on a shared host.', '',
             '| Case | Class | Status | Seconds (median) | Iterations | Energy (Ha) |', '|---|---|---|---:|---:|---:|']
    groups = {}
    for key, r in data['cases'].items():
        lines.append(f"| {key} | {r['case']['system_class']} | {r['status']} | {r['wall_median_s']:.3f} | {r['iterations']} | {r['total_energy_ha']} |")
        g = groups.setdefault(r['case']['system_class'], [0,0,0.0])
        g[0] += 1; g[1] += r['status'] == 'pass'; g[2] += r['wall_median_s']
    lines += ['', '## By system class', '', '| Class | Cases | Passed | Sum median seconds |', '|---|---:|---:|---:|']
    for key, g in sorted(groups.items()):
        lines.append(f'| {key} | {g[0]} | {g[1]} | {g[2]:.3f} |')
    lines += ['', '## Baseline comparison', '', 'Only selected cases are compared; missing selected entries fail comparison.']
    lines += ['- ' + x['id'] + ': ' + '; '.join(x['regressions'] + x['improvements'] or ['unchanged']) for x in data['comparison']]
    lines += ['', '## Findings and exact reproducing inputs', '',
              'Replace SI_UPF_PATH with the read-only validation file. For synthetic cases regenerate synthetic.UPF using the case parameters in JSON. Runtime directories contain executable inputs and logs.']
    for key, r in data['cases'].items():
        if r['status'] == 'pass':
            continue
        lines += ['', f'### {key}', '', f"Expected failure matched: {r['expected_failure']}", '', '```text']
        lines += [run['error'] for run in r['runs']]
        lines += ['```', '', '```toml', r['input_toml'].rstrip(), '```']
    return '\n'.join(lines) + '\n'


def metadata(binary):
    def command(args):
        p = subprocess.run(args, cwd=ROOT, text=True, capture_output=True)
        return p.stdout.strip() if p.returncode == 0 else None
    return dict(timestamp_utc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
                platform=platform.platform(), machine=platform.machine(), processor=platform.processor(),
                cpu_count=os.cpu_count(), python=platform.python_version(),
                commit=command(['git', 'rev-parse', 'HEAD']), dirty=command(['git', 'status', '--porcelain']),
                binary=str(binary), binary_sha256=file_hash(binary),
                generator_sha256=file_hash(ROOT / 'scripts/gen_synthetic_upf.py'),
                runner_sha256=file_hash(Path(__file__)),
                toolchain=(ROOT / 'rust-toolchain.toml').read_text(),
                cargo_lock_sha256=file_hash(ROOT / 'Cargo.lock'))


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--binary', type=Path, default=ROOT / 'target/release/dft-cli')
    p.add_argument('--upf', type=Path)
    p.add_argument('--case-dir', type=Path, default=ROOT / 'benchmarks/cases')
    p.add_argument('--cases', help='comma-separated ids; explicit selection includes long cases')
    p.add_argument('--include-long', action='store_true')
    p.add_argument('--robustness', action='store_true')
    p.add_argument('--seed', type=int, default=1)
    p.add_argument('--threads', type=int, default=1)
    p.add_argument('--repeats', type=int, default=1)
    p.add_argument('--variants', type=int, default=0, help='random displaced variants per selected case')
    p.add_argument('--displacement', type=float, default=0.01, help='fractional displacement amplitude')
    p.add_argument('--timeout', type=float, help='override each per-run timeout, seconds')
    p.add_argument('--baseline', type=Path)
    p.add_argument('--require-identical', action='store_true')
    p.add_argument('--output', type=Path, required=True, help='new JSON report (never overwritten)')
    p.add_argument('--work-dir', type=Path, help='new directory for exact inputs, logs and CLI output')
    a = p.parse_args(argv)
    if a.seed < 0 or a.threads < 0 or a.repeats < 1 or a.variants < 0 or not 0 <= a.displacement < 0.5:
        p.error('invalid seed/threads/repeats/variants/displacement')
    if a.timeout is not None and not 0 < a.timeout < math.inf:
        p.error('timeout must be positive and finite')
    try:
        cases = [load_case(path) for path in sorted(a.case_dir.glob('*.toml'))]
        if len({c['id'] for c in cases}) != len(cases):
            raise ValueError('duplicate case id')
        selected = a.cases.split(',') if a.cases else None
        if selected and set(selected) - {c['id'] for c in cases}:
            raise ValueError('unknown selected case')
        cases = [c for c in cases if (c['id'] in selected if selected else a.include_long or c['tier'] == 'default')]
        if not cases:
            raise ValueError('empty case selection')
        binary = a.binary.resolve(strict=True)
        upf = a.upf.resolve(strict=True) if a.upf else None
        if any(c['potential'] == 'si' for c in cases) and upf is None:
            raise ValueError('--upf required for Si cases')
        baseline = json.loads(a.baseline.read_text())['cases'] if a.baseline else {}
        md_path = a.output.with_suffix('.md')
        if a.output.exists() or md_path.exists() or a.output == md_path:
            raise ValueError('output JSON/Markdown paths must be distinct and new')
        if a.work_dir:
            work = a.work_dir.resolve()
            work.mkdir(parents=True, exist_ok=False)
        else:
            work = Path(tempfile.mkdtemp(prefix='dft-bench-suite-')).resolve()
        expanded = []
        for c in cases:
            expanded.append(c)
            expanded.extend(variant(c, a.seed, i, a.displacement) for i in range(a.variants))
        if len({c['id'] for c in expanded}) != len(expanded):
            raise ValueError('variant IDs collide with case IDs')
        data = dict(schema_version=1, metadata=metadata(binary), seed=a.seed, threads=a.threads,
                    repeats=a.repeats, work_dir=str(work), energy_tolerance_ha=ENERGY_TOL_HA,
                    invocation=sys.argv if argv is None else argv,
                    cases={}, comparison=[], aborted=False)
        status = 0
        for c in expanded:
            r = run_case(binary, c, upf, work / c['id'], a.seed, a.threads, a.repeats, a.timeout)
            data['cases'][c['id']] = r
            print(f"{c['id']}: {r['status']}, {r['wall_median_s']:.3f} s", flush=True)
            if r['status'] != 'pass' and not r['expected_failure'] and not a.robustness:
                data['aborted'] = True
                status = 3
                break
        if a.baseline:
            data['comparison'] = compare(data['cases'], baseline, a.require_identical)
            if any(c['regressions'] for c in data['comparison']):
                status = 5
        a.output.parent.mkdir(parents=True, exist_ok=True)
        a.output.write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')
        md_path.write_text(markdown(data))
        print(f'Report: {a.output}; reproducers: {work}')
        return status
    except (OSError, ValueError, KeyError) as e:
        print(f'suite setup/I/O error: {e}', file=sys.stderr)
        return 1


if __name__ == '__main__':
    sys.exit(main())
