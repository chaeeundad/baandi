//! Cross-layer characterization through the CLI process, without introducing
//! engine/CLI dependencies into dft-solvers. Synthetic UPFs reproduce the
//! Gaussian-ion formulas in dft-engine/tests/common/mod.rs (Python math.erf).
//! Release command and retained evidence are documented in the comparison note.
use std::{path::Path, process::Command};

#[test]
#[ignore = "full SCF matrix; build release dft-cli, then run the documented release command"]
fn material_mixing_matrix() {
    let root = Path::new(env!("CARGO_MANIFEST_DIR")).join("../..");
    let status = Command::new("python3")
        .args(["-I", "-c", DRIVER])
        .current_dir(root)
        .status()
        .expect("Python 3 is required by the cross-layer characterization harness");
    assert!(
        status.success(),
        "SCF characterization or its invariants failed"
    );
}

const DRIVER: &str = r#"
import hashlib, importlib.util, json, math, os, pathlib, statistics, subprocess, tempfile, time
root = pathlib.Path.cwd()
binary = pathlib.Path(os.environ.get('DFT_CLI_BINARY', str(root/'target/release/dft-cli'))).resolve()
real_upf = pathlib.Path(os.environ.get('DFT_SI_UPF', '/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF'))
assert binary.is_file(), 'build release dft-cli first'
assert real_upf.is_file(), 'set DFT_SI_UPF to the read-only validation UPF'
out = pathlib.Path(os.environ.get('DFT_MIXING_EVIDENCE') or tempfile.mkdtemp(prefix='p10d-scf-')).resolve()
out.mkdir(parents=True, exist_ok=True)
repeats = int(os.environ.get('DFT_MIXING_REPEATS', '3'))
assert repeats > 0

def synthetic(z, nonlocal_):
    points, end, r0, rc = 1201, 12., .02, 1.
    dx = math.log1p(end/r0)/(points-1)
    r = [r0*math.expm1(i*dx) for i in range(points)]
    arr = lambda v: ' '.join(format(x,'.17e') for x in v)
    beta = ''
    if nonlocal_:
        for i,l in enumerate([0,0,1]):
            beta += '<PP_BETA.{0} index="{0}" angular_momentum="{1}" cutoff_radius="12" cutoff_radius_index="{2}">{3}</PP_BETA.{0}>'.format(i+1,l,points,arr([x**(l+1)*math.exp(-x*x/(i+1)) for x in r]))
        beta += '<PP_DIJ>2 0.4 0 0.4 4 0 0 0 6</PP_DIJ>'
    xml = '''<UPF version="2.0.1"><PP_INFO>SYNTHETIC Gaussian ion, no real-element reference.</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ" z_valence="{z}" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="1" l_max_rho="0" l_local="-1" mesh_size="{points}" number_of_wfc="0" number_of_proj="{np}" />
<PP_MESH dx="{dx}" mesh="{points}" xmin="0" rmax="12" zmesh="1"><PP_R>{r}</PP_R><PP_RAB>{rab}</PP_RAB></PP_MESH><PP_LOCAL>{local}</PP_LOCAL><PP_NONLOCAL>{beta}</PP_NONLOCAL><PP_RHOATOM>{rho}</PP_RHOATOM></UPF>'''.format(z=z, points=points, np=3 if nonlocal_ else 0, dx=dx, r=arr(r), rab=arr([(x+r0)*dx for x in r]), local=arr([-4*z/(rc*math.sqrt(math.pi)) if x==0 else -2*z*math.erf(x/rc)/x for x in r]), beta=beta, rho=arr([0.]*points))
    path = out/('synthetic-z%d-p%d.UPF'%(z,3 if nonlocal_ else 0))
    path.write_text(xml)
    return path

def input_text(name, lattice, sites, upf, label='Xx', ecut=10., mesh=None, bands=None, grid=None):
    text = '[system]\nname='+json.dumps(name)+'\nlattice_unit="bohr"\nlattice_vectors='+json.dumps(lattice)+'\n'
    for site in sites:
        text += '\n[[system.atoms]]\nspecies='+json.dumps(label)+'\nfractional='+json.dumps(site)+'\n'
    text += '\n[pseudopotentials]\n'+label+'='+json.dumps(str(upf))+'\n[basis]\necutwfc_ry='+str(ecut)+'\n'
    text += 'kpoints=[[0.0,0.0,0.0]]\n' if mesh is None else 'kgrid='+json.dumps(mesh)+'\nkshift=[0.5,0.5,0.5]\n'
    if grid is not None: text += 'fft_grid='+json.dumps(grid)+'\n'
    text += '\n[xc]\nfunctional="lda-pz"\n[scf]\nthreads=1\nmax_iterations=150\nmixing="pulay"\nmixing_beta=0.3\nenergy_tolerance_ha=1e-8\ndensity_residual_tolerance=1e-7\n'
    if bands is not None: text += 'nbands='+str(bands)+'\n'
    return text

z2, z2nl, z4nl = synthetic(2,False), synthetic(2,True), synthetic(4,True)
a=10.261212856717933
cases = {
 'synthetic-cubic-z2-p0': input_text('cubic Z2 local', [[8.,0.,0.],[0.,8.,0.],[0.,0.,8.]], [[0.,0.,0.]], z2, bands=5),
 'synthetic-skew-z2-p3-k211': input_text('skew elongated Z2 projectors', [[6.,0.,0.],[1.,7.,0.],[.5,.7,16.]], [[0.,0.,0.],[.25,.25,.25]], z2nl, mesh=[2,1,1], bands=8),
 'synthetic-slab-z2-p0': input_text('synthetic slab vacuum', [[6.,0.,0.],[0.,6.,0.],[0.,0.,24.]], [[0.,0.,.45],[.5,.5,.55]], z2, bands=8),
 'synthetic-fcc-z4-p3': input_text('fcc Z4 projectors', [[0.,a/2,a/2],[a/2,0.,a/2],[a/2,a/2,0.]], [[0.,0.,0.],[.25,.25,.25]], z4nl, bands=8),
 'synthetic-even-skew-z2-p0-k221': input_text('even grid skew Z2', [[6.,0.,0.],[1.,6.,0.],[.5,.7,10.]], [[0.,0.,0.]], z2, mesh=[2,2,1], bands=12, grid=[24,24,32]),
}
spec=importlib.util.spec_from_file_location('bench',root/'scripts/bench_scf.py'); bench=importlib.util.module_from_spec(spec); spec.loader.exec_module(bench)
for key in ['gamma','k222s','k222','k444s']:
    cases['si-'+key]=bench.case_input(str(real_upf),key,1)
# Replicate the primitive cell along its third lattice vector, then stretch it.
prim=[[0.,a/2,a/2],[a/2,0.,a/2],[a/2,a/2,0.]]
long=[prim[0],prim[1],[v*2*1.15 for v in prim[2]]]
sites=[[x,y,(z+t)/2] for t in [0,1] for x,y,z in [[0.,0.,0.],[.25,.25,.25]]]
cases['si-supercell-stretch']=input_text('Si 1x1x2 primitive, stretch 1.15',long,sites,real_upf,'Si',ecut=20.,bands=12)
cases['si-supercell-stretch-wide']=input_text('Si stretched supercell, 20 bands',long,sites,real_upf,'Si',ecut=20.,bands=20)
conv=[[0.,0.,0.],[.25,.25,.25],[0.,.5,.5],[.5,0.,.5],[.5,.5,0.],[.25,.75,.75],[.75,.25,.75],[.75,.75,.25]]
cases['si-slab-vacuum']=input_text('Si conventional cell plus 2 cell lengths vacuum',[[a,0.,0.],[0.,a,0.],[0.,0.,3*a]],[[x,y,(z+1)/3] for x,y,z in conv],real_upf,'Si',ecut=12.,bands=20)
beta_override=os.environ.get('DFT_MIXING_BETA')
if beta_override:
    assert math.isfinite(float(beta_override)) and float(beta_override)>0
    cases={k:v.replace('mixing_beta=0.3','mixing_beta='+beta_override).replace('mixing_beta = 0.5','mixing_beta = '+beta_override) for k,v in cases.items()}
selected=os.environ.get('DFT_MIXING_CASES')
if selected:
    names=selected.split(','); assert all(k in cases for k in names); cases={k:cases[k] for k in names}
records=[]

def run(name, text, method, repeat):
    folder=out/(name+'-'+method+'-'+str(repeat));folder.mkdir(exist_ok=True)
    inp=folder/'input.toml';inp.write_text(text)
    t0=time.perf_counter()
    result=subprocess.run([str(binary),'run',str(inp),'--output-dir',str(folder/'out')],capture_output=True,text=True,timeout=600)
    wall=time.perf_counter()-t0
    (folder/'stdout.txt').write_text(result.stdout);(folder/'stderr.txt').write_text(result.stderr)
    row={'case':name,'strategy':method,'repeat':repeat,'exit':result.returncode,'wall_s':wall,'input':str(inp),'stderr':result.stderr[:400]+'\n...\n'+result.stderr[-400:]}
    assert result.returncode in [0,3,4], row
    report=folder/'out/results.json'
    if report.exists():
        r=json.loads(report.read_text());h=r['convergence']['history'];n=r['electrons']
        row.update(iterations=len(h),energy_ha=r['total_energy']['ha'],raw_residual=h[-1]['density_residual'],max_charge_error_electrons=max(max(v['relative_charge_error'],v['relative_input_charge_error'])*n for v in h),min_density=min(v['min_density_electrons_per_bohr3'] for v in h),resets=h[-1]['mixing_resets'],raw_recoveries=sum(v['mixing_raw_recovery'] for v in h),engine_s=r['timings']['engine_total_s'],mixing_s=r['timings']['phases']['mixing_s'],certificate=h[-1]['fixed_point_certificate'],density_sha256=hashlib.sha256((folder/'out/density.bin').read_bytes()).hexdigest(),history=h,parameters=r['scf']['mixing'])
        assert row['max_charge_error_electrons'] < 1e-10, row['max_charge_error_electrons']
        assert row['min_density'] >= -1e-14
        count=0;resets=0
        for i,v in enumerate(h):
            ratio=None if i==0 or h[i-1]['density_residual']==0 else v['density_residual']/h[i-1]['density_residual']
            assert v['density_contraction']==ratio
            if v['mixing_history_used'] is not None:
                expected=0 if r['scf']['mixing']['method']=='linear' else min(count,r['scf']['mixing']['history'])
                assert v['mixing_history_used']==expected
                count+=1
            if v['mixing_raw_recovery']:
                assert r['scf']['mixing']['method']=='kerker-pulay' and v['mixing_history_used']==0
            if v['mixer_reset']: resets+=1;count=0
            assert v['mixing_resets']==resets
        if result.returncode==0:
            last=h[-1];cert=last['fixed_point_certificate'];assert cert is not None
            assert last['density_residual']<=1e-7 and abs(last['delta_energy_ha'])<=1e-8 and last['max_eigen_residual_ha']<=1e-8
            assert cert['density_residual']<=1e-7 and cert['occupied_eigenvalue_change_ha']<=1e-7
    else:
        assert result.returncode==4 and any(word in result.stderr for word in ['Diverged','NonInsulating','Numerical']), row
    records.append(row)
    (out/'scf-matrix.json').write_text(json.dumps({'binary':str(binary),'real_upf_sha256':hashlib.sha256(real_upf.read_bytes()).hexdigest(),'repeats':repeats,'records':records},indent=2))
    print(name,method,'repeat',repeat,'exit',row['exit'],'iterations',row.get('iterations'),'energy',row.get('energy_ha'),'wall_s',round(wall,4),flush=True)
    return row

for name,text in cases.items():
    for method,q in [('pulay',None),('kerker-0.5',.5),('kerker-1',1.),('kerker-2',2.)]:
        effective=text if q is None else text.replace('mixing="pulay"','mixing="kerker-pulay"\nmixing_q0_inverse_bohr='+str(q)).replace('mixing = "pulay"','mixing = "kerker-pulay"\nmixing_q0_inverse_bohr='+str(q))
        group=[run(name,effective,method,i) for i in range(repeats)]
        assert len({v['exit'] for v in group})==1
        if 'energy_ha' in group[0]:
            assert len({(v['energy_ha'],v['density_sha256'],v['iterations']) for v in group})==1
    success=[v for v in records if v['case']==name and v['exit']==0]
    base=next((v for v in success if v['strategy']=='pulay'),None)
    if base:
        assert all(abs(v['energy_ha']-base['energy_ha'])<=1e-12 for v in success)
    if name in ['synthetic-cubic-z2-p0','si-gamma']:

        assert all(v['exit']==0 for v in records if v['case']==name)

# Small-filtered-residual stagnation must not be mistaken for convergence.
text=bench.case_input(str(real_upf),'gamma',1)
text=text.replace('max_iterations = 100','max_iterations = 3').replace('mixing = "pulay"','mixing = "kerker-pulay"\nmixing_q0_inverse_bohr=1e100')
r=run('certificate-raw-residual',text,'kerker-large-q',0)
assert r['exit']==3 and r['raw_residual']>1e-7
# Explicit legacy divergence remains a nonzero, descriptive failure.
text=bench.case_input(str(real_upf),'gamma',1).replace('mixing = "pulay"','mixing = "linear"').replace('mixing_beta = 0.5','mixing_beta = 100.0')
r=run('divergence',text,'linear-100',0)
assert r['exit']==4 and 'Diverged' in r['stderr']
text=input_text('explicit new strategy divergence',[[0.,a/2,a/2],[a/2,0.,a/2],[a/2,a/2,0.]],[[0.,0.,0.],[.25,.25,.25]],z4nl,bands=8)
text=text.replace('mixing="pulay"','mixing="kerker-pulay"\nmixing_q0_inverse_bohr=1.0').replace('mixing_beta=0.3','mixing_beta=1e308')
r=run('divergence',text,'kerker-extreme-beta',0)
assert r['exit']==4 and 'Diverged' in r['stderr']
print('Retained complete parameters, residual/energy histories and certificates:',out/'scf-matrix.json',flush=True)
"#;
