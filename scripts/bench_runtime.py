#!/usr/bin/env python3
"""Matched release runtime benchmark. Standard library only; see runtime/README.md.

Exit 0: measurements completed (NOT performance acceptance); 5: bit mismatch.
Builds an independent archived pre-extraction revision without modifying git.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import statistics
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parents[1]
DRIVER = ROOT / 'tests/reference/runtime/bench_driver.rs'
MODES = {'classic': 0.5, 'observed': 1.0, 'agent': 1.0}


def sha(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def command(args, cwd=ROOT, **kwargs):
    result = subprocess.run(args, cwd=cwd, text=True, capture_output=True, **kwargs)
    if result.returncode:
        raise RuntimeError(f'{args}: exit {result.returncode}: {result.stderr[-12000:]}')
    return result.stdout.strip()


def optional_command(args):
    try:
        return command(args)
    except (OSError, RuntimeError) as error:
        return 'unavailable: ' + str(error)


def distribution(values):
    ordered = sorted(values)
    def percentile(p):
        pos = (len(ordered) - 1) * p
        lo = int(pos)
        return ordered[lo] + (ordered[min(lo + 1, len(ordered) - 1)] - ordered[lo]) * (pos - lo)
    return {'n': len(values), 'median': statistics.median(values), 'p5': percentile(.05), 'p95': percentile(.95)}


def compare(ref, candidate, target, seed, sufficient):
    overhead = [100 * (c / r - 1) for r, c in zip(ref, candidate)]
    rng = random.Random(seed)
    samples = sorted(statistics.median(rng.choices(overhead, k=len(overhead))) for _ in range(10000))
    ci = [samples[249], samples[9749]]
    status = ('met' if ci[1] < target else 'not_met' if ci[0] >= target else 'inconclusive') if sufficient else 'insufficient_samples'
    return {'target_percent': target, 'overhead_percent': distribution(overhead), 'paired_overhead_percent': overhead,
            'median_overhead_ci95_percent': ci, 'status': status}


def build(source, binary, pre_refactor, env, target):
    args = ['cargo', 'build', '--offline', '--release', '-p', 'dft-cli', '--message-format=json',
            '--target-dir', str(target), '--config',
            'source.vendored-sources.directory=' + json.dumps(str((ROOT / '../rust-dft-vendor').resolve()))]
    output = command(args, cwd=source, env=env)
    libs = {}
    for line in output.splitlines():
        item = json.loads(line)
        if item.get('reason') == 'compiler-artifact':
            for name in item['filenames']:
                if name.endswith('.rlib'):
                    libs[item['target']['name']] = name
    needed = ['dft_engine', 'dft_geometry', 'dft_physics', 'dft_upf', 'serde_json', 'sha2']
    if not pre_refactor:
        needed.append('dft_cli')
    args = ['rustc', '--edition=2024', '-C', 'opt-level=3', str(DRIVER), '-o', str(binary),
            '-L', 'dependency=' + str(target / 'release/deps')]
    for name in needed:
        args += ['--extern', name + '=' + libs[name]]
    if pre_refactor:
        args += ['--cfg', 'pre_refactor']
    command(args, env=env)
    return {'binary_sha256': sha(binary), 'cargo_command': 'cargo build --offline --release -p dft-cli',
            'rustc_command': args}


def run(binary, mode, kind, budget, threads, probe='normal'):
    args = [str(binary), mode, str(kind), str(budget), str(threads), probe]
    # wait4 gives per-child RSS (not cumulative RUSAGE_CHILDREN high water).
    with tempfile.TemporaryFile() as output, tempfile.TemporaryFile() as error:
        start = time.perf_counter_ns()
        child = subprocess.Popen(args, stdout=output, stderr=error)
        _, status, usage = os.wait4(child.pid, 0)
        wall = time.perf_counter_ns() - start
        child.returncode = os.waitstatus_to_exitcode(status)
        output.seek(0)
        error.seek(0)
        if child.returncode:
            raise RuntimeError(f'{args}: exit {child.returncode}: {error.read().decode()}')
        result = json.load(output)
    result.update(process_wall_ns=wall, process_cpu_s=usage.ru_utime + usage.ru_stime,
                  peak_rss_bytes=usage.ru_maxrss * (1 if sys.platform == 'darwin' else 1024), command=args, exit_code=0)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--reference-revision', default='9495d9fbaaf0e0bf6bf141ea61675e1d803eeeec')
    parser.add_argument('--output', type=Path, default=ROOT / 'tests/reference/runtime/benchmark.json')
    parser.add_argument('--repeats', type=int, default=100)
    parser.add_argument('--warmups', type=int, default=3)
    parser.add_argument('--threads', type=int, default=1)
    parser.add_argument('--seed', type=int, default=20261104)
    args = parser.parse_args()
    if args.repeats < 1 or args.warmups < 0 or args.threads < 1:
        parser.error('positive repeats/threads and nonnegative warmups required')
    env = dict(os.environ, PATH='/opt/homebrew/opt/rustup/bin:' + os.environ['PATH'], CARGO_NET_OFFLINE='true')
    revision = command(['git', 'rev-parse', args.reference_revision])
    evidence = {'schema_version': 1, 'task_id': 'P11J-RUNTIME-BENCH', 'status': 'running',
                'command': [sys.executable, *sys.argv], 'reference_revision': revision,
                'candidate_revision': command(['git', 'rev-parse', 'HEAD']),
                'git_status': command(['git', 'status', '--short']), 'platform': platform.platform(),
                'cpu': optional_command(['sysctl', '-n', 'machdep.cpu.brand_string']) if sys.platform == 'darwin' else platform.processor(),
                'rustc': command(['rustc', '-Vv'], env=env), 'load_start': os.getloadavg(),
                'requested_threads': args.threads, 'warmups': args.warmups, 'repeats': args.repeats,
                'bootstrap_seed': args.seed, 'bootstrap_resamples': 10000, 'exclusions': [],
                'timer_resolution_s': time.get_clock_info('perf_counter').resolution,
                'baselines': {str(p.relative_to(ROOT)): {'sha256': sha(p), 'data': json.loads(p.read_text())}
                              for p in sorted((ROOT / 'docs/performance').glob('baseline*.json'))},
                'sources': {str(p.relative_to(ROOT)): sha(p) for p in
                            [DRIVER, Path(__file__), ROOT / 'Cargo.lock', ROOT / 'crates/dft-engine/tests/common/mod.rs']},
                'cases': [], 'limitations': [
                    'Shared local host: quiet dev08 acceptance must be repeated by Claude.',
                    'Historical Si baselines are provenance, not matched synthetic timing denominators.',
                    'No allocation trace or isolated internal event construction/enqueue timers.',
                    'No async persistence adapter, slow-subscriber reconnect journal, or network control transport measured.',
                    'No slab/elongated, large-band, systematically varied gap or real Si fixture in this matrix.',
                    'Process wall/RSS include fingerprint/output instrumentation identically in both builds.',
                    'Checkpoint model hashes use existing synthetic test identity; not production provenance.' ]}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    def save():
        args.output.write_text(json.dumps(evidence, indent=2) + '\n')
    save()  # Replace a stale successful report before any build/run.
    try:
        with tempfile.TemporaryDirectory(prefix='dft-runtime-') as temporary:
            tmp = Path(temporary)
            archive = tmp / 'reference.tar'
            with archive.open('wb') as stream:
                subprocess.run(['git', 'archive', revision], cwd=ROOT, stdout=stream, check=True)
            reference = tmp / 'reference'
            reference.mkdir()
            command(['tar', '-xf', str(archive), '-C', str(reference)])
            target = ROOT / 'target'
            print('Building candidate and independent pre-extraction driver', flush=True)
            evidence['candidate_build'] = build(ROOT, tmp / 'candidate', False, env, target)
            # The reference revision must not share the candidate's target directory: both
            # build package `dft-cli` and overwrite each other's artifacts across runs.
            evidence['reference_build'] = build(reference, tmp / 'reference-driver', True, env, ROOT / 'target' / 'bench-reference')
            all_identical = True
            for kind in range(3):
                for budget in [1, 200]:
                    case = {'kind': kind, 'budget': budget, 'pairs': [], 'probes': []}
                    for mode, threshold in MODES.items():
                        pairs = []
                        for index in range(args.warmups + args.repeats):
                            pair = {}
                            order = ['reference', mode] if index % 2 == 0 else [mode, 'reference']
                            for arm in order:
                                pair[arm] = run(tmp / ('reference-driver' if arm == 'reference' else 'candidate'),
                                                arm, kind, budget, args.threads)
                            if index >= args.warmups:
                                pair['order'] = order
                                pairs.append(pair)
                        identical = all(p['reference']['identity'] == p[mode]['identity'] for p in pairs)
                        # Compare repeat trajectories as well as each matched pair.
                        identical &= all(p['reference']['identity'] == pairs[0]['reference']['identity'] for p in pairs)
                        all_identical &= identical
                        timing = compare([p['reference']['total_ns'] for p in pairs], [p[mode]['total_ns'] for p in pairs],
                                         threshold, args.seed, args.repeats >= 20 and args.warmups >= 3)
                        if not identical:
                            timing['status'] = 'equivalence_failed'
                        case['pairs'].append({'mode': mode, 'identity': identical, 'in_process_total': timing,
                                              'process_wall': compare([p['reference']['process_wall_ns'] for p in pairs],
                                                                      [p[mode]['process_wall_ns'] for p in pairs], threshold,
                                                                      args.seed, args.repeats >= 20 and args.warmups >= 3),
                                              'rss_delta_bytes': distribution([p[mode]['peak_rss_bytes'] - p['reference']['peak_rss_bytes'] for p in pairs]),
                                              'samples': pairs})
                        print(f'kind={kind} budget={budget} {mode}: {timing["status"]} {timing["median_overhead_ci95_percent"]}, identical={identical}', flush=True)
                    if budget > 1:
                        expected = case['pairs'][0]['samples'][0]['reference']['identity']
                        for probe in ['saturation', 'control', 'checkpoint', 'events']:
                            for _ in range(args.repeats):
                                value = run(tmp / 'candidate', 'agent', kind, budget, args.threads, probe)
                                value['probe_name'] = probe
                                value['identity_matches_reference'] = value['identity'] == expected
                                all_identical &= value['identity_matches_reference']
                                case['probes'].append(value)
                    evidence['cases'].append(case)
                    save()
            evidence['all_identity_checks_pass'] = all_identical
            evidence['status'] = 'measured' if all_identical else 'equivalence_failed'
            evidence['load_end'] = os.getloadavg()
            save()
            return 0 if all_identical else 5
    except Exception as error:
        evidence['status'] = 'failed'
        evidence['error'] = str(error)
        save()
        raise


if __name__ == '__main__':
    sys.exit(main())
