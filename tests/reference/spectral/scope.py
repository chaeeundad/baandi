#!/usr/bin/env python3
"""Frozen-H scope probe; compile with P15C_ORACLE=1. Never modifies SCF states.
Caps fixed before measurement: spectrum and occupied band sum each 1e-9 Ha.
"""
import argparse
import gzip
import importlib.util
import json
import os
from pathlib import Path
import re
import subprocess
import tempfile

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location('replay', ROOT / 'tests/p15a_replay.py')
replay = importlib.util.module_from_spec(spec)
spec.loader.exec_module(replay)
CAP = 1e-9


def compare(values, dense, occupations):
    errors = [abs(a-b) for a,b in zip(values,dense)]
    occupied = sum(f > 0 for f in occupations)
    # Match multiplicities, not merely the closest eigenvalue for each root.
    available = set(range(len(dense)))
    matched = []
    for value in values:
        within = [j for j in available if abs(dense[j]-value) <= CAP]
        i = min(within) if within else min(available, key=lambda j: abs(dense[j]-value))
        if abs(dense[i]-value) <= CAP:
            available.remove(i)
            matched.append(i)
        else:
            matched.append(None)
    clusters = []
    for i, value in enumerate(dense):
        if not clusters or abs(value-dense[clusters[-1][0]]) > CAP:
            clusters.append([i])
        else:
            clusters[-1].append(i)
    missing = sorted(i for i in available if i < len(values))
    return dict(occupied_max_error_ha=max(errors[:occupied], default=0),
                first_empty_error_ha=errors[occupied] if occupied < len(errors) else None,
                buffer_max_error_ha=max(errors[occupied+1:], default=0),
                band_sum_excess_ha=sum(f*(a-b) for f,a,b in zip(occupations,values,dense)),
                missing_dense_indices_zero_based=missing if all(i is not None for i in matched) else None,
                unmatched_dense_indices_zero_based=missing,
                missing_indices_certified=all(i is not None for i in matched),
                missing_clusters=[c for c in clusters if set(c) & set(missing)],
                matched_dense_indices_zero_based=matched,
                all_roots_within_cap=max(errors, default=0) <= CAP)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--binary', type=Path, required=True)
    ap.add_argument('--baseline', type=Path, default=Path('/tmp/p15c-baseline-cli'))
    ap.add_argument('--output', type=Path, default=OUT)
    ap.add_argument('--cases', default='si8,k444s,gamma,k222s,k222,syn20-skew,syn4-elongated,syn4-slab,syn16-gamma,syn16-k211,si-20,mgo-40')
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    for name in args.cases.split(','):
        with tempfile.TemporaryDirectory(prefix='p15c-') as td:
            work = Path(td)
            upf = ROOT.parent / 'pseudo/Si.pz-vbc.UPF'
            hashes = {}
            if name.startswith(('si-', 'mgo-')):
                text = (ROOT / 'tests/reference/materials' / (name+'.toml')).read_text().replace('threads = 4', 'threads = 1')
                nc = Path(os.environ.get('DFT_NC_DIR', ROOT.parent / 'pseudo/pseudodojo/ex'))
                for element in set(re.findall(r'species = "([A-Za-z]+)"', text)):
                    path = (nc / (element+'.upf')).resolve()
                    hashes[element] = replay.sha(path)
                    assert hashes[element] in text
                    text = text.replace('"'+element+'.upf"', json.dumps(str(path)))
            elif name in ('gamma','k222','k222s'):
                text = replay.load_module('bench_scf').case_input(str(upf), name, threads=1)
                hashes['Si'] = replay.sha(upf)
            else:
                text, pseudo = replay.inputs(name, work, upf)
                hashes['pseudo'] = replay.sha(pseudo)
            inp = work / 'input.toml'
            inp.write_text(text)
            base, br, bd, bc = replay.run(args.baseline, inp, work/'base')
            assert br is not None, base.stderr
            nk = br.get('kpoints', {}).get('count', 1)
            hist = br['convergence']['history']
            ids, locations, sid = [], {}, 0
            for row in hist:
                for k in range(nk):
                    locations[sid+k] = dict(iteration=row['iteration'], k=k, phase='ordinary')
                if row['iteration'] in (1,2,4,(len(hist)+1)//2,len(hist)):
                    ids += [sid+k for k in sorted(set((0,nk//2,nk-1)))]
                sid += nk
                if row.get('fixed_point_certificate'):
                    for k in range(nk):
                        locations[sid+k] = dict(iteration=row['iteration'], k=k, phase='certificate')
                    ids += [sid+k for k in sorted(set((0,nk//2,nk-1)))]
                    sid += nk
            os.environ['P15C_DENSE_IDS'] = ','.join(map(str, ids))
            proc, report, density, command = replay.run(args.binary, inp, work/'oracle')
            solves, _ = replay.parse(proc.stderr)
            assert report is not None and solves, proc.stderr[-4000:]
            summaries = []
            for i, solve in enumerate(solves):
                if 'dense_eigenvalues_ha' not in solve:
                    continue
                native = solve[0]['steps'][-1]
                variants = {}
                for variant, trace in solve.items():
                    if isinstance(variant, int):
                        last = trace['steps'][-1]
                        variants[variant] = dict(success=trace['success'],h_vectors=last['h_vectors'],
                            **compare(last['eigenvalues_ha'],solve['dense_eigenvalues_ha'],br['occupations']))
                summaries.append(dict(solve=i, **locations[i], residuals_ha=native['residuals_ha'],variants=variants,
                    **compare(native['eigenvalues_ha'],solve['dense_eigenvalues_ha'],br['occupations'])))
            result = dict(case=name, eigenvalue_cap_ha=CAP, band_sum_cap_ha=CAP,
                commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
                input=text, upf_sha256=hashes, binary_sha256=replay.sha(args.binary),
                baseline_binary_sha256=replay.sha(args.baseline), baseline_command=bc, command=command,
                dense_ids=ids, exit_code=proc.returncode, baseline_exit_code=base.returncode,
                converged=report['convergence']['converged'], iterations=len(hist),
                total_energy_ha=report['total_energy']['ha'],
                energy_identical=report['total_energy']==br['total_energy'],density_identical=density==bd,
                history_identical=replay.scientific(report['convergence'])==replay.scientific(br['convergence']),
                comparisons=summaries, stderr_tail='\n'.join(l for l in proc.stderr.splitlines() if not l.startswith('P15A_'))[-1000:])
            (args.output/(name+'-traces.json.gz')).write_bytes(gzip.compress(json.dumps(solves,separators=(',',':')).encode(),mtime=0))
            (args.output/(name+'.json')).write_text(json.dumps(result,indent=2)+'\n')
            final = [r for r in summaries if r['iteration']==len(hist)]
            print(name, 'converged',result['converged'],'identity',result['energy_identical'],result['density_identical'], 'final',[(r['phase'],r['k'],r['occupied_max_error_ha'],r['buffer_max_error_ha']) for r in final],flush=True)
            assert result['energy_identical'] and result['density_identical'] and result['history_identical']

if __name__ == '__main__':
    main()
