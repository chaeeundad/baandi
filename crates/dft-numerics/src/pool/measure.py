#!/usr/bin/env python3
"""Run from the repository root. Analysis binary only; raw CSV timestamps are ns.
Usage: python3 -I crates/dft-numerics/src/pool/measure.py BINARY OUTPUT_DIR UPF
Each child performs one SCF, preserving input, trace, result and density evidence.
Phases are inclusive; their times must not be added. Busy time covers instrumented
leaf kernels only, not all main-thread work. Imbalance includes idle pool workers.
"""
import importlib.util,pathlib,subprocess,json,hashlib,time,collections,sys,csv
root=pathlib.Path.cwd(); spec=importlib.util.spec_from_file_location('bench',root/'scripts/bench_scf.py'); b=importlib.util.module_from_spec(spec);spec.loader.exec_module(b)
binary,tag,upf=sys.argv[1:];binary=str(pathlib.Path(binary).resolve());upf=str(pathlib.Path(upf).resolve());base=pathlib.Path(tag).resolve();base.mkdir(exist_ok=True)
def union(spans):
    end=-1; total=0
    for a,z in sorted(spans):
        total+=max(0,z-max(a,end));end=max(end,z)
    return total/1e9
def merged(spans):
    result=[]
    for a,z in sorted(spans):
        if result and a <= result[-1][1]: result[-1][1]=max(result[-1][1],z)
        else: result.append([a,z])
    return result

def overlap(a,b):
    i=j=total=0
    while i<len(a) and j<len(b):
        total+=max(0,min(a[i][1],b[j][1])-max(a[i][0],b[j][0]))
        if a[i][1]<b[j][1]: i+=1
        else: j+=1
    return total/1e9

rows=[]
for case in ['si8','skew4']:
 for t in [1,2,4,8,12]:
    d=base/f'{case}-{t}';d.mkdir(exist_ok=True)
    text=b.case_input(upf,'si8' if case=='si8' else 'gamma',t)
    if case=='skew4':
        text=text.replace('[[0.0, 2.715, 2.715],','[[0.0, 5.43, 5.43],')
        a=text.index('[[system.atoms]]');z=text.index('[pseudopotentials]')
        text=text[:a]+''.join(f'[[system.atoms]]\nspecies = "Si"\nfractional = {site}\n\n' for site in [[0,0,0],[.125,.25,.25],[.5,0,0],[.625,.25,.25]])+text[z:]
    inp=d/'input.toml';inp.write_text(text)
    start=time.monotonic()
    with (d/'trace.csv').open('w') as err:
        p=subprocess.run([binary,'run',str(inp),'--output-dir',str(d/'out')],stdout=subprocess.DEVNULL,stderr=err)
    assert p.returncode==0,(case,t,(d/'trace.csv').read_text()[-1000:])
    events=[]
    for line in (d/'trace.csv').read_text().splitlines():
        if not line.startswith('P10L,'):continue
        fields=line.split(',');_,name,a,z,w,n,busy=fields[:7];events.append(dict(name=name,a=int(a),z=int(z),w=int(w),n=int(n),busy=list(map(int,busy.split(';'))) if busy else [],nested_busy=len(fields)>7 and fields[7]=='1'))
    result=json.loads((d/'out/results.json').read_text());span=next(e for e in events if e['name']=='scf');regions=[e for e in events if e['w']]
    # Preserve individual-region imbalance, including idle workers. Serial
    # fallback occupies the external slot (255); do not silently discard it.
    with (d/'regions.csv').open('w') as f:
        writer=csv.writer(f)
        writer.writerow(['phase','start_ns','end_ns','scheduled_slots','pieces','max_over_mean','busy_fraction_by_pool_worker','external_busy_fraction'])
        for e in regions:
            duration=e['z']-e['a'];total=sum(e['busy'])
            writer.writerow([e['name'],e['a'],e['z'],e['w'],e['n'],max(e['busy'])*t/total if total else '',
                ';'.join(str(v/duration) if duration else '' for v in e['busy'][:t]),e['busy'][255]/duration if duration else ''])
    phases={}
    parallel=merged((e["a"],e["z"]) for e in regions if e["w"]>1)
    logical=merged((e["a"],e["z"]) for e in regions)
    for name in sorted(set(e['name'] for e in events)):
        ev=[e for e in events if e['name']==name];wall=union((e['a'],e['z']) for e in ev)
        busy=[sum(e['busy'][i] for e in ev if e['busy'])/1e9 for i in range(256)]
        pieces=[e['n'] for e in ev if e['w']]
        phases[name]=dict(wall_s=wall,count=len(ev),busy_s=sum(busy),max_mean=max(busy)*t/sum(busy) if sum(busy) else None,pieces_min=min(pieces) if pieces else None,pieces_max=max(pieces) if pieces else None,worker_s=busy if any(busy) else [])
        phases[name]['parallel_regions_s']=overlap(merged((e['a'],e['z']) for e in ev),parallel)
        phases[name]['between_logical_regions_s']=wall-overlap(merged((e['a'],e['z']) for e in ev),logical)
    wall=(span['z']-span['a'])/1e9;inside=union((max(e['a'],span['a']),min(e['z'],span['z'])) for e in regions if e['z']>span['a'] and e['a']<span['z'])
    row=dict(case=case,threads=t,scf_s=wall,regions_s=inside,between_s=wall-inside,busy_s=sum(sum(e['busy']) for e in regions)/1e9,phases=phases,energy=result['total_energy']['ha'],iterations=len(result['convergence']['history']),density=hashlib.sha256((d/'out/density.bin').read_bytes()).hexdigest())
    row['parallel_regions_s']=union(parallel)
    row['observed_busy_fraction_by_worker']=[sum(e['busy'][i] for e in events if e['busy'] and not e['nested_busy'])/1e9/wall for i in range(t)]
    row['observed_external_busy_fraction']=sum(e['busy'][255] for e in events if e['busy'] and not e['nested_busy'])/1e9/wall
    rows.append(row);(base/'summary.json').write_text(json.dumps(rows,indent=2));print(case,t,round(wall,4),round(inside,4),round(row['busy_s']/wall/t,3),flush=True)
