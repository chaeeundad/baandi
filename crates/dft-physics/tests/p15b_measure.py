#!/usr/bin/env python3
"""Temporary source instrumentation and alternating-pair measurement for P15B.

Run from repository root. `build` restores the production file in a finally
block; it never edits git metadata. Output paths must be outside version control.
Do not run `build` concurrently with another Cargo invocation in this worktree.
Instrumentation counts logical column-element reads (not hardware cache misses)
and times nonlocal/H calls; it is absent from the delivered numerical hot path.
"""
import argparse
import csv
import json
import os
from pathlib import Path
import statistics
import subprocess
import tempfile
import time
import importlib.util

SOURCE = Path('crates/dft-physics/src/hamiltonian.rs')
INSTRUMENT = r'''
mod p15b_probe {
    use std::sync::atomic::{AtomicU64, Ordering::Relaxed};
    use std::time::Instant;
    static H_NS: AtomicU64 = AtomicU64::new(0);
    static H_VECTORS: AtomicU64 = AtomicU64::new(0);
    static NL_NS: AtomicU64 = AtomicU64::new(0);
    static VECTORS: AtomicU64 = AtomicU64::new(0);
    static BYTES: AtomicU64 = AtomicU64::new(0);
    static OLD_BYTES: AtomicU64 = AtomicU64::new(0);
    pub struct Timer(Instant, bool);
    impl Timer {
        pub fn h(vectors: usize) -> Self { H_VECTORS.fetch_add(vectors as u64, Relaxed); Self(Instant::now(), true) }
        pub fn nl(vectors: usize, bytes: usize) -> Self {
            VECTORS.fetch_add(vectors as u64, Relaxed);
            BYTES.fetch_add(bytes as u64, Relaxed);
            OLD_BYTES.fetch_add((bytes * vectors) as u64, Relaxed);
            Self(Instant::now(), false)
        }
    }
    impl Drop for Timer {
        fn drop(&mut self) {
            let counter = if self.1 { &H_NS } else { &NL_NS };
            counter.fetch_add(self.0.elapsed().as_nanos() as u64, Relaxed);
        }
    }
    pub fn dump() {
        eprintln!("P15B {},{},{},{},{},{}", H_NS.load(Relaxed), NL_NS.load(Relaxed), VECTORS.load(Relaxed), BYTES.load(Relaxed), OLD_BYTES.load(Relaxed), H_VECTORS.load(Relaxed));
    }
}
impl Drop for NonlocalOperator {
    fn drop(&mut self) {
        let columns: usize = self.atoms.iter().map(|a| a.columns.len()).sum();
        eprintln!("P15B_SHAPE {},{},{},{}", self.basis.len(), self.atoms.len(), columns, self.basis.len() * columns * 16);
        p15b_probe::dump();
    }
}
'''


def instrument(source):
    source += INSTRUMENT
    for signature in ['fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {',
                      'fn apply_block(&self, x: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {']:
        assert signature in source
        source = source.replace(signature, signature + ('\nlet _timer = p15b_probe::Timer::h(x.cols());' if 'apply_block' in signature else '\nlet _timer = p15b_probe::Timer::h(1);'))
    signature = 'fn add(&self, x: &[Complex64], y: &mut [Complex64], scratch: &mut [Complex64]) {'
    source = source.replace(signature, signature + '\nlet _timer = p15b_probe::Timer::nl(1, self.atoms.iter().map(|a| a.columns.len() * x.len() * 32).sum());')
    if 'fn add_block<const G:' in source:
        pos = source.index('        debug_assert!(bands <= G);')
        source = source[:pos] + 'let _timer = (bands > 1).then(|| p15b_probe::Timer::nl(bands, self.atoms.iter().map(|a| a.columns.len() * rows * 32).sum()));\n' + source[pos:]
    return source


def build(out, features=None):
    current = SOURCE.read_text()
    try:
        for name, source in [('baseline', subprocess.check_output(['git', 'show', 'HEAD:' + str(SOURCE)], text=True)), ('candidate', current)]:
            SOURCE.write_text(instrument(source))
            command = ['cargo', 'build', '--release', '-p', 'dft-cli']
            if features:
                command += ['--features', features]
            subprocess.run(command, check=True)
            dest = out / (name + '-probe')
            dest.write_bytes(Path('target/release/dft-cli').read_bytes())
            dest.chmod(0o755)
    finally:
        SOURCE.write_text(current)


def measure(out, upf, pairs):
    spec = importlib.util.spec_from_file_location('bench', 'scripts/bench_scf.py')
    bench = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(bench)
    binaries = [out / 'baseline-probe', out / 'candidate-probe']
    result = {'head': subprocess.check_output(['git', 'rev-parse', 'HEAD'], text=True).strip(),
              'upf_sha256': bench.sha256(upf), 'binaries': {str(b): bench.sha256(b) for b in binaries},
              'load_start': os.getloadavg(), 'samples': []}
    for threads in [8, 1]:
        for case in ['si8', 'k444s']:
            for pair in range(pairs):
                for mode in ([0, 1] if pair % 2 == 0 else [1, 0]):
                    with tempfile.TemporaryDirectory() as td:
                        td = Path(td)
                        inp = td / 'input.toml'
                        inp.write_text(bench.case_input(upf, case, threads))
                        start = time.perf_counter()
                        p = subprocess.run([str(binaries[mode]), 'run', str(inp), '--output-dir', str(td / 'out')], capture_output=True, text=True, check=True)
                        wall = time.perf_counter() - start
                        report = json.loads((td / 'out/results.json').read_text())
                        counters = [list(map(int, line[5:].split(','))) for line in p.stderr.splitlines() if line.startswith('P15B ')]
                        assert counters
                        h_ns, nl_ns, vectors, actual, old, h_vectors = counters[-1]
                        history = report['convergence']['history']
                        row = {'case': case, 'threads': threads, 'pair': pair, 'mode': mode,
                               'wall_s': wall, 'h_us_per_vector': h_ns / h_vectors / 1000,
                               'nl_us_per_vector': nl_ns / vectors / 1000, 'vectors': vectors, 'h_vectors': h_vectors,
                               'projector_bytes': actual, 'scalar_projector_bytes': old,
                               'energy_ha': report['total_energy']['ha'],
                               'density_sha256': bench.sha256(td / 'out/density.bin'),
                               'scf_iterations': len(history), 'load': os.getloadavg(),
                               'projector_shapes': sorted(set(line for line in p.stderr.splitlines() if line.startswith('P15B_SHAPE ')))}
                        result['samples'].append(row)
                        (out / 'pairs.json').write_text(json.dumps(result, indent=2) + '\n')
                print(case, threads, pair + 1, flush=True)
    fields = ['case', 'threads', 'pair', 'mode', 'wall_s', 'h_us_per_vector',
              'nl_us_per_vector', 'h_vectors', 'vectors', 'projector_bytes',
              'scalar_projector_bytes', 'energy_ha', 'density_sha256', 'scf_iterations']
    with (out / 'pairs.csv').open('w', newline='') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, extrasaction='ignore')
        writer.writeheader()
        writer.writerows(result['samples'])
    for case in ['si8', 'k444s']:
        for threads in [1, 8]:
            rows = [r for r in result['samples'] if r['case'] == case and r['threads'] == threads]
            assert len({r['energy_ha'] for r in rows}) == 1
            assert len({r['density_sha256'] for r in rows}) == 1
            for mode in [0, 1]:
                subset = [r for r in rows if r['mode'] == mode]
                print(case, threads, mode, {k: statistics.median(r[k] for r in subset) for k in ['h_us_per_vector', 'nl_us_per_vector', 'wall_s', 'vectors', 'projector_bytes', 'scalar_projector_bytes']})

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=['build', 'measure'])
    parser.add_argument('--out', type=Path, required=True)
    parser.add_argument('--upf', type=str)
    parser.add_argument('--pairs', type=int, default=30)
    parser.add_argument('--features', help='Cargo features for build, e.g. mkl on dev08')
    args = parser.parse_args()
    args.out.mkdir(parents=True, exist_ok=True)
    if args.action == 'build':
        build(args.out, args.features)
    else:
        measure(args.out, str(Path(args.upf).resolve()), args.pairs)
