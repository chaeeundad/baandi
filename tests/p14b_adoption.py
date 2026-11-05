#!/usr/bin/env python3
"""Validate the architect's P14B acceptance rule without changing legacy baselines."""
import hashlib
import json
from pathlib import Path
import runpy

HERE = Path(__file__).resolve().parent
ROOT = HERE / 'reference/perf/p14b'
OUT = ROOT / 'adopted'
records = runpy.run_path(str(HERE / 'p14b_summarize.py'))['records']


def read(path):
    return json.loads(path.read_text())['cases']


def main():
    old, new = read(ROOT/'suite-final.json'), read(OUT/'suite1.json')
    lines = ['# P14B architect-adopted incremental Gram', '',
             'Realistic fixtures retain SCF counts and |ΔE| ≤ 1e-12 Ha. Converged toy fixtures retain the suite energy tolerance (1e-12 Ha); iteration changes are reported. Existing failures are retained, not counted as converged. No new convergence failures are permitted.', '',
             '| Case | Before SCF | Adopted SCF | Δ iterations | ΔE (Ha) | Status |',
             '| --- | ---: | ---: | ---: | ---: | --- |']
    results = []
    assert len(new) == 19
    for name, a in new.items():
        b = old[name]
        assert a['status'] == b['status'], name
        de = a['total_energy_ha'] - b['total_energy_ha']
        if a['converged']:
            assert abs(de) <= 1e-12, (name, de)
            if name.startswith('si'):
                assert a['iterations'] == b['iterations'], name
        lines.append(f"| {name} | {b['iterations']} | {a['iterations']} | {a['iterations']-b['iterations']:+d} | {de:.6e} | {a['status']} |")
        for threads in (2,4,8):
            t = read(OUT/f'suite{threads}.json')[name]
            for key in ('status','iterations','total_energy_ha','density_sha256','operator_applications'):
                assert a[key] == t[key], (name, threads, key)
            assert a['runs'][0]['eigenvalues_ha'] == t['runs'][0]['eigenvalues_ha']
        results.append(dict(case=name, status=a['status'], iterations_before=b['iterations'],
                            iterations_after=a['iterations'], delta_energy_ha=de, threads_identical=True))
    lines += ['', 'All 19 suite cases have bit-identical energy, density hashes, eigenvalues and H counts at T=1,2,4,8, including retained nonconverged outputs. The legacy harness exits 5 for accepted toy iteration changes at T=1; its cross-thread comparisons also exit 5 because thread counts are part of input identity. This independent test checks scientific determinism directly.', '',
              '## Additional realistic fixtures', '',
              '| Case | Before / adopted SCF | ΔE (Ha) | T=1,2,4,8 identical |', '| --- | --- | ---: | --- |']
    original = read(ROOT/'default-final.json')
    current = read(OUT/'threads1.json')
    original.update(read(OUT/'cubic-full-gram.json'))
    current.update(read(OUT/'cubic12-threads1.json'))
    for name, a in current.items():
        b = original[name]
        ma, mb = a['measurements'], b['measurements']
        assert ma['converged'] == mb['converged']
        de = ma['total_energy_ha'] - mb['total_energy_ha']
        for threads in (2,4,8):
            f = f'cubic12-threads{threads}.json' if name == 'syn20-cubic12' else f'threads{threads}.json'
            t = read(OUT/f)[name]
            assert a['density_sha256'] == t['density_sha256']
            for key in ('total_energy_ha','iterations','operator_applications','eigenvalues_ha'):
                assert ma[key] == t['measurements'][key]
        if ma['converged']:
            assert ma['iterations'] == mb['iterations'], name
            assert abs(de) <= 1e-12, (name,de)
        lines.append(f"| {name} | {mb['iterations']} / {ma['iterations']} | {de:.6e} | yes |")
    lines += ['', 'The dense 8-Bohr cubic 20-atom fixture is nonconverged before and after (100 iterations), so its energy is not an accepted converged comparison. The 12-Bohr cubic fixture converges in 13 iterations in both builds; its H-vector work increases 1996→2107, despite unchanged SCF count. Both input definitions and both outcomes are retained.', '', '## Adopted GEMM accounting and next-site bound', '',
              '| Case | Before 2mnk | Adopted 2mnk | Reduction % | Ritz/residual 2mnk | Second projection pass 2mnk |',
              '| --- | ---: | ---: | ---: | ---: | ---: |']
    measured = read(OUT/'accounting.json')
    for name, a in measured.items():
        if name not in read(ROOT/'accounting-final.json'):
            continue
        before = sum(r['flops'] for r in records(read(ROOT/'accounting-final.json')[name]))
        rs = list(records(a)); total = sum(r['flops'] for r in rs)
        ritz = sum(r['flops'] for r in rs if r['site']=='ritz_and_residual')
        second = sum(r['flops'] for r in rs if r['site'].endswith('_2'))
        assert a['density_sha256'] == current[name]['density_sha256']
        lines.append(f'| {name} | {before} | {total} | {100*(1-total/before):.2f} | {ritz} | {second} |')
    lines += ['', 'Ritz/residual products are the next-largest site. For the current unlocked convergence test, every root residual must be recomputed: guaranteed removable 2mnk is zero without changing the algorithm. The listed full site totals are ceilings on potential savings, not attainable savings. Second-pass projection totals are similarly ceilings; testing the normalized block still needs its second overlap, and omitting the subtract/MGS would require a new orthogonality gate. No second optimization is adopted.', '',
              'MKL executes the same two rectangular products with borrowed contiguous input columns and one sequential cblas call per rectangle. It adds only reduced-matrix copies: 16*p*p bytes to preserve the prefix plus 16*(q*q-p*p) bytes to scatter outputs, i.e. 16*q*q bytes per extension. This is an explicit small-copy increase versus full direct-output MKL GEMM, not a claim of reduced copied bytes. MKL runtime/determinism/timings require dev08.', '']
    (OUT/'acceptance.md').write_text('\n'.join(lines))
    (OUT/'acceptance.json').write_text(json.dumps(dict(rule='Architect decision, 2026-11-05',cases=results,
        all_previously_converged_cases_converge=True, existing_nonconvergence=[n for n,a in new.items() if not a['converged']],
        pure_rust_threads=[1,2,4,8], mkl_runtime_verified=False),indent=2)+'\n')
    source_files=['crates/dft-numerics/src/linalg/tiled.rs','crates/dft-solvers/src/lib.rs','tests/p14a_work_accounting.py']
    repo=HERE.parent
    (OUT/'source-sha256.json').write_text(json.dumps({p:hashlib.sha256((repo/p).read_bytes()).hexdigest() for p in source_files},indent=2)+'\n')
    print('Architect acceptance: PASS; 16 converged suite cases, 3 retained failures; all 19 thread-identical')


if __name__ == '__main__':
    main()
