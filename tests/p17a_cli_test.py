"""Release CLI/profile/checkpoint contracts; generated UPF remains temporary."""
import json
from pathlib import Path
import struct
import subprocess
import sys
import tempfile

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'))
import equal_accuracy as ea
BINARY=ROOT/'target/release/dft-cli'


def run(*args, ok=True):
    p=subprocess.run([str(BINARY),*map(str,args)],capture_output=True,text=True)
    assert (p.returncode==0)==ok, (args,p.returncode,p.stderr[:1000])
    return p


def read_cp(path):
    data=path.read_bytes();size=struct.unpack('<Q',data[8:16])[0]
    return json.loads(data[24:24+size]),data[24+size:-32]


def scientific(value):
    if isinstance(value,dict):return {k:scientific(v) for k,v in value.items() if k not in ('timings','timestamp')}
    if isinstance(value,list):return list(map(scientific,value))
    return value


with tempfile.TemporaryDirectory(prefix='p17a-cli-') as d:
    d=Path(d)
    text,_=ea.case_input('syn4-elongated',d,Path('/unused'),1)
    for metric in ['l2','hartree']:
        inp=d/(metric+'.toml')
        inp.write_text(ea.scf_settings(text,dict(convergence_profile='standard',mixing_metric=metric)))
        output=d/(metric+'-run')
        p=run('run',inp,'--output-dir',output)
        assert 'convergence_profile: standard; certified: false' in p.stdout
        r=json.loads((output/'results.json').read_text())
        assert r['schema_version']==1 and r['certified'] is False
        assert r['convergence_profile']=='standard' and r['standard_threshold_ha']==5e-9
        assert r['convergence']['converged']
        h=r['convergence']['history']
        assert all(x['fixed_point_certificate'] is None for x in h)
        assert h[-1]['standard_error_ha']<5e-9
        assert all(a['requested_eigen_tolerance_ha']>=b['requested_eigen_tolerance_ha'] for a,b in zip(h,h[1:]))
        assert all(x['converged'] for x in r['band_accuracy'][:sum(x>0 for x in r['occupations'])+1])
        assert r['work']['h_vectors']==sum(x['operator_applications'] for x in h)
        exact,direct,resumed=[d/(metric+x+'.cp') for x in ('-exact','-direct','-resumed')]
        run('checkpoint',inp,exact,'--iteration','2')
        run('checkpoint',inp,direct,'--kind','result')
        run('restart',inp,exact,resumed)
        a,ab=read_cp(direct);b,bb=read_cp(resumed)
        assert a['standard']['convergence_profile']=='standard' and a['standard']['certified'] is False
        assert scientific(a['result'])==scientific(b['result'])
        assert ab==bb, 'exact restart must preserve density and orbital bits'
        warm=d/(metric+'-warm.cp');warm_result=d/(metric+'-warm-result.cp')
        run('checkpoint',inp,warm,'--kind','warm','--iteration',str(min(5,len(h)-1)))
        run('restart',inp,warm,warm_result,'--warm')
        wh,_=read_cp(warm_result)
        assert wh['result']['history'][0]['standard_redo'], 'near-fixed-point warm start must trigger the first redo'
        strict=d/'strict.toml';strict.write_text(text)
        run('restart',strict,exact,d/'mismatch.cp',ok=False)
    for settings in [dict(convergence_profile='typo'),dict(convergence_profile='standard',standard_threshold_ha=0.),
                     dict(mixing_metric='typo'),dict(standard_threshold_ha=5e-9)]:
        inp=d/'bad.toml';inp.write_text(ea.scf_settings(text,settings))
        run('run',inp,'--output-dir',d/'bad',ok=False)
print('PASS standard labels, redo, ratchet, roots, counters, both metrics, exact checkpoint bits, incompatible options and invalid inputs')
