//! Isolated-source measurement runner; production sources and SCF results stay unchanged.
//! P37A_SI_UPF=/absolute/Si.pz-vbc.UPF DFT_NC_DIR=/absolute/pseudodojo/ex
//! cargo test --release -p dft-engine --test generalized_gram_stage0 -- --ignored --nocapture
//! P37A_CASES, P37A_PROFILE and P37A_START_RUN/P37A_END_RUN select a bounded batch.
//! Evidence and the instrumented source copy are retained in P37A-evidence/.
use std::{path::PathBuf, process::Command};

#[path = "../../dft-solvers/tests/generalized_gram_stage0.rs"]
mod payload;

#[test]
#[ignore = "release-only SCF histories; creates isolated measurement sources and invokes cargo"]
fn histories() {
    if cfg!(debug_assertions) {
        panic!("run this measurement in release mode");
    }
    let root = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../..");
    let evidence = root.join("P37A-evidence");
    std::fs::create_dir_all(&evidence).unwrap();
    std::fs::write(evidence.join("probe.rs"), payload::PROBE).unwrap();
    std::fs::write(evidence.join("run.py"), DRIVER).unwrap();
    std::fs::write(evidence.join("analyze.py"), ANALYSIS).unwrap();
    assert!(
        Command::new("python3")
            .args(["-I", "P37A-evidence/run.py"])
            .current_dir(root)
            .status()
            .unwrap()
            .success()
    );
}

const DRIVER: &str = r###"
import pathlib, shutil, os, subprocess, json, hashlib
root=pathlib.Path.cwd().resolve(); evidence=root/'P37A-evidence'; dst=evidence/'source'
# Never copy target, .git, artifacts, external datasets or UPFs.
if not dst.exists():
    dst.mkdir()
    for name in ['Cargo.toml','Cargo.lock','rust-toolchain.toml']:
        shutil.copy2(root/name,dst/name)
    shutil.copytree(root/'crates',dst/'crates',ignore=shutil.ignore_patterns('target','*.upf','*.UPF','__pycache__'))
    (dst/'.cargo').mkdir()
    (dst/'.cargo/config.toml').write_text('[source.crates-io]\nreplace-with="vendored-sources"\n[source.vendored-sources]\ndirectory='+json.dumps(str(root.parent/'rust-dft-vendor'))+'\n')
    for name in ['benchmarks','tests']:
        (dst/name).symlink_to(root/name,target_is_directory=True)
# Refuse reuse after production sources or manifests have changed.
identity_files=sorted(set(root.glob('crates/*/src/**/*.rs'))|set(root.glob('crates/*/Cargo.toml'))|{root/'Cargo.toml',root/'Cargo.lock'})
identity={str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in identity_files}
identity_file=dst/'production-identity.json'
if identity_file.exists():assert json.loads(identity_file.read_text())==identity,'measurement copy belongs to different production sources'
else:identity_file.write_text(json.dumps(identity,indent=2))
# Reset only instrumented copies from production on every invocation.
sol=dst/'crates/dft-solvers'; eng=dst/'crates/dft-engine'
manifest=(root/'crates/dft-solvers/Cargo.toml').read_text()
(sol/'Cargo.toml').write_text(manifest+'\nfaer.workspace = true\n')
shutil.copy2(evidence/'probe.rs',sol/'src/p37a.rs')
s=(root/'crates/dft-solvers/src/lib.rs').read_text()
needle='        let mut history = Vec::new();\n        let mut operator_applications = 0;'
assert s.count(needle)==1
s=s.replace(needle,'        let _p37a_id = crate::p37a::complex(operator, &basis, preconditioner, self.bands, policy, self.tolerance, ritz, self.max_iterations);\n'+needle)
s+='\nmod p37a;\n'
(sol/'src/lib.rs').write_text(s)
s=(root/'crates/dft-solvers/src/real_davidson.rs').read_text()
needle='        let mut history = Vec::with_capacity(self.max_iterations);'
assert s.count(needle)==1
s=s.replace(needle,'        let p37a_id = crate::p37a::real(op, &basis, diagonal, self.bands, policy, self.tolerance, ritz, self.max_iterations);\n        let _p37a_guard = P37aCount::begin(p37a_id);\n'+needle)
s=s.replace('            max_basis_dimension = max_basis_dimension.max(basis.cols());','            P37aCount::step(basis.cols());\n            max_basis_dimension = max_basis_dimension.max(basis.cols());')
s=s.replace('ws.multiply_into(', 'p37a_multiply(ws, ')
s=s.replace('                        let p = dot(&self.column, column(&self.output, k));','                        P37aCount::scalar(2 * n);\n                        let p = dot(&self.column, column(&self.output, k));')
s+=r'''
// Measurement-copy instrumentation of each real solver GEMM call and MGS dot/update.
thread_local! { static P37A_COUNT: std::cell::RefCell<(usize,usize,u64,u64,u64)> = const {std::cell::RefCell::new((0,0,0,0,0))}; }
struct P37aCount;
impl P37aCount {
    fn begin(id:usize)->Self{P37A_COUNT.with(|x|*x.borrow_mut()=(id,0,0,0,0)); Self}
    fn step(m:usize){P37A_COUNT.with(|x|{let mut x=x.borrow_mut(); x.1+=1;x.4+=2*(m as u64).pow(3);});}
    fn scalar(n:usize){P37A_COUNT.with(|x|x.borrow_mut().3+=n as u64);}
}
impl Drop for P37aCount {fn drop(&mut self){P37A_COUNT.with(|x|{let x=x.borrow();println!("P37A_BASE,{},{},{},{},{}",x.0,x.1,x.2,x.3,x.4);});}}
fn p37a_multiply(ws:&mut RealGemmWorkspace,a:RealMatrixRef<'_>,op:RealOp,b:RealMatrixRef<'_>,out:dft_numerics::real::RealMatrixMut<'_>)->Result<(),LinalgError>{
    let (m,k)=if op==RealOp::Identity{(a.nrows(),a.ncols())}else{(a.ncols(),a.nrows())};
    P37A_COUNT.with(|x|x.borrow_mut().2+=(m*k*b.ncols()) as u64);
    ws.multiply_into(a,op,b,out)
}
'''
(sol/'src/real_davidson.rs').write_text(s)
# Reuse the established P31A geometry/start driver, preserving the solver defaults.
s=(root/'crates/dft-engine/tests/capacity_trajectory.rs').read_text().replace('P31A','P37A')
s=s.replace('unwrap_or("7".into())','unwrap_or("3".into())')
s=s.replace('let real = std::env::var("P37A_REAL").is_ok_and(|v| v == "on");','let real = mesh.len() == 1 && mesh[0].fractional == [0.;3];')
# The A52 fixture is generated in-memory from the existing synthetic generator.
s=s.replace('let text = std::fs::read_to_string(root.join(file)).unwrap();', '''let text = if case == "syn20-skew" {
            "lattice = [14,0,0,1.5,14,0,0.5,0.7,14]\\ncutoff_ry = 20\\natoms = [".to_owned()
            + &(0..20).flat_map(|i| [(i%4) as f64/4.+0.07,(i/4) as f64/5.+0.03,if i%2==0 {0.35}else{0.65}]).map(|v|v.to_string()).collect::<Vec<_>>().join(",") + "]\\n"
        } else {std::fs::read_to_string(root.join(file)).unwrap()};
        let synthetic = if case == "syn20-skew" {
            let output=std::process::Command::new("python3").args(["-c", "import importlib.util; s=importlib.util.spec_from_file_location('g','scripts/gen_synthetic_upf.py');g=importlib.util.module_from_spec(s);s.loader.exec_module(g);print(g.generate(2,(0,0,1),points=101))"]).current_dir(std::env::var("P37A_ROOT").unwrap()).output().unwrap();
            assert!(output.status.success());Some(dft_upf::parse_upf(std::str::from_utf8(&output.stdout).unwrap()).unwrap())
        }else{None};''')
s=s.replace('potential: &pp,','potential: synthetic.as_ref().unwrap_or(&pp),')
s=s.replace('        for run in first..last {', '        for run in first..if case == "si8-gamma" {last}else{first+1} {')
s=s.replace('            if material {\n                opts.functional', '            if case == "syn20-skew" {opts.max_iterations=6; opts.bands=Some(24);}\n            if material {\n                opts.functional')
# Runs 2..7 = five seeded starts; strict uses the original strict defaults.
s=s.replace('            if run >= 2 {', '            if run >= 2 && std::env::var("P37A_UNPERTURBED").is_err() {')
(eng/'tests/p37a_driver.rs').write_text(s)
if os.environ.get('P37A_CONTROL')=='1':
    for name in ['lib.rs','real_davidson.rs']:
        shutil.copy2(root/'crates/dft-solvers/src'/name,sol/'src'/name)
    shutil.copy2(root/'crates/dft-solvers/Cargo.toml',sol/'Cargo.toml')
env=os.environ.copy();env['CARGO_TARGET_DIR']=str(root/'target');env['P37A_ROOT']=str(root)
env.setdefault('P37A_SI_UPF',str(root.parent/'pseudo/Si.pz-vbc.UPF'))
env.setdefault('DFT_NC_DIR',str(root.parent/'pseudo/pseudodojo/ex'))
env.setdefault('P37A_CASES','si8-gamma');env.setdefault('P37A_START_RUN','2');env.setdefault('P37A_END_RUN','3')
cmd=['cargo','test','--offline','--release','--manifest-path',str(dst/'Cargo.toml'),'-p','dft-engine','--test','p37a_driver','sweep','--','--ignored','--nocapture','--test-threads=1']
label=env['P37A_CASES'].replace(',','_')+'-'+env.get('P37A_PROFILE','standard')+'-'+env['P37A_START_RUN']+'-'+env['P37A_END_RUN']
if env.get('P37A_UNPERTURBED')=='1':label+='-plain'
if env.get('P37A_CONTROL')=='1':label='control-'+label
log=evidence/(label+'.log')
with log.open('w') as f:
    result=subprocess.run(cmd,cwd=dst,env=env,stdout=f,stderr=subprocess.STDOUT)
(evidence/(label+'.command.json')).write_text(json.dumps({'command':cmd,'exit_code':result.returncode,'environment':{k:v for k,v in env.items() if k.startswith(('P37A_','DFT_NC','CARGO_TARGET'))},'production_sources':{str(p.relative_to(root)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [root/'crates/dft-solvers/src/lib.rs',root/'crates/dft-solvers/src/real_davidson.rs',root/'crates/dft-engine/src/gamma_real.rs']}},indent=2))
(evidence/(label+'.provenance.json')).write_text(json.dumps({'head':subprocess.check_output(['git','rev-parse','HEAD'],cwd=root,text=True).strip(),'upf_sha256':{str(p):hashlib.sha256(p.read_bytes()).hexdigest() for p in [pathlib.Path(env['P37A_SI_UPF'])]+[pathlib.Path(env['DFT_NC_DIR'])/(x+'.upf') for x in ['Si','Mg','O']] if p.exists()}},indent=2))
print(log, result.returncode, flush=True)
if result.returncode: print(log.read_text()[-10000:]);raise SystemExit(result.returncode)
"###;

const ANALYSIS: &str = r###"
import pathlib,json,statistics,collections,math,csv
root=pathlib.Path(__file__).resolve().parent
groups=collections.defaultdict(list); outcomes=[]; commands=[]
for manifest in sorted(root.glob('*.command.json')):
    if manifest.name.startswith(('pilot-','control-')):continue
    meta=json.loads(manifest.read_text()); commands.append(meta)
    if meta['exit_code']!=0:continue
    log=manifest.with_name(manifest.name.replace('.command.json','.log'))
    pending={};base={}
    for line in log.read_text().splitlines():
        if 'P37A_GRAM,' in line:
            x=line[line.index('P37A_GRAM,'):].split(','); sid=int(x[1]);cap=int(x[4]);
            pending.setdefault((sid,cap),[]).append(dict(step=int(x[5]),n=int(x[6]),m=int(x[7]),rank=int(x[8]),smallest=float(x[9]),relative=float(x[10]),condition=float(x[11]),cholesky=x[12]=='true',rectangle=int(x[13]),rotation=int(x[14]),contraction=int(x[15]),small_model=int(x[16]),state=x[17],gram_drift=float(x[18]) if len(x)>18 else 0.))
        elif line.startswith('P37A_BASE,'):
            x=line.split(',');base[int(x[1])]=dict(steps=int(x[2]),gemm=int(x[3]),mgs=int(x[4]),small_model=int(x[5]))
        elif line.startswith('SUMMARY,'):
            x=line.split(',');case,seed,profile=x[1],(0 if meta['environment'].get('P37A_UNPERTURBED')=='1' else int(x[2])-1),x[4]
            for (sid,cap),rows in pending.items(): groups[(case,profile,cap)].append(dict(seed=seed,solve=sid,rows=rows,base=base.get(sid)))
            pending={};base={}
        elif line.startswith('END,'):
            x=line.split(',');outcomes.append(dict(case=x[1],seed=(0 if meta['environment'].get('P37A_UNPERTURBED')=='1' else int(x[2])-1),profile=meta['environment'].get('P37A_PROFILE','standard'),status=','.join(x[3:])))
assert groups,'no completed measurements'
def q(a,p):
    a=sorted(a);x=(len(a)-1)*p;i=int(x);return a[i]+(a[min(i+1,len(a)-1)]-a[i])*(x-i)
summary=[]
for (case,profile,cap),solves in sorted(groups.items()):
    rows=[r for s in solves for r in s['rows']];rel=[r['relative'] for r in rows];val=[r['smallest'] for r in rows];n=len(rows)
    r=dict(case=case,profile=profile,capacity=cap,seeds=sorted(set(s['seed'] for s in solves)),solves=len(solves),steps=n,min_eigenvalue=min(val),p01_eigenvalue=q(val,.01),p10_eigenvalue=q(val,.1),p50_eigenvalue=q(val,.5),min_relative=min(rel),p01_relative=q(rel,.01),p10_relative=q(rel,.1),p50_relative=q(rel,.5),cholesky_fail_fraction=sum(not x['cholesky'] for x in rows)/n,fallback_1e6_fraction=sum(x['relative']<1e-6 or not x['cholesky'] for x in rows)/n,fallback_1e8_fraction=sum(x['relative']<1e-8 or not x['cholesky'] for x in rows)/n,max_rank_deficiency=max(x['m']-x['rank'] for x in rows),diagnostic_unfinished=sum(s['rows'][-1]['state']!='converged' for s in solves),rank_exhausted=sum(s['rows'][-1]['state']=='rank_exhausted' for s in solves),max_gram_drift=max(x['gram_drift'] for x in rows))
    for field in ['rectangle','rotation','contraction','small_model']:r[field+'_per_step']=sum(x[field] for x in rows)/n
    r['generalized_mac_per_step']=sum(r[f+'_per_step'] for f in ['rectangle','rotation','contraction'])
    if all(s['base'] for s in solves):
        bsteps=sum(s['base']['steps'] for s in solves); bmac=sum(s['base']['gemm']+s['base']['mgs'] for s in solves);bsmall=sum(s['base']['small_model'] for s in solves)
        r.update(baseline_steps=bsteps,baseline_mac_per_step=bmac/bsteps,baseline_small_model_per_step=bsmall/bsteps,mac_saving_fraction=1-r['generalized_mac_per_step']/(bmac/bsteps),model_inclusive_saving_fraction=1-(r['generalized_mac_per_step']+r['small_model_per_step'])/((bmac+bsmall)/bsteps),total_mac_ratio=r['generalized_mac_per_step']*n/bmac)
        savings=[1-(sum(x['rectangle']+x['rotation']+x['contraction'] for x in s['rows'])/len(s['rows']))/((s['base']['gemm']+s['base']['mgs'])/s['base']['steps']) for s in solves]
        r['worst_solve_mac_saving_fraction']=min(savings);r['solves_below_40pct_saving']=sum(x<.4 for x in savings)
    summary.append(r)
(root/'summary.json').write_text(json.dumps(dict(rows=summary,outcomes=outcomes),indent=2,allow_nan=False))
with (root/'summary.csv').open('w') as f:
    keys=list(dict.fromkeys(k for r in summary for k in r));w=csv.DictWriter(f,keys);w.writeheader();w.writerows(summary)
lines=['## Measured overlap statistics','', 'Eigenvalues and relative ratios are dimensionless. Percentiles use linear interpolation across steps; all five Si8 seeds are pooled; A52 pools the uniform start and seed 1 when both are present. Fallback includes Cholesky failure or a relative eigenvalue below the threshold.','', '| Case | Profile | Capacity | Solves / steps | min lambda | p1 | p10 | p50 | min relative | Cholesky fail % | fallback 1e-6 % | fallback 1e-8 % | max rank deficit | unfinished solves |','|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
for r in summary:lines.append(f"| {r['case']} | {r['profile']} | {r['capacity']} | {r['solves']} / {r['steps']} | {r['min_eigenvalue']:.3e} | {r['p01_eigenvalue']:.3e} | {r['p10_eigenvalue']:.3e} | {r['p50_eigenvalue']:.3e} | {r['min_relative']:.3e} | {100*r['cholesky_fail_fraction']:.2f} | {100*r['fallback_1e6_fraction']:.2f} | {100*r['fallback_1e8_fraction']:.2f} | {r['max_rank_deficiency']} | {r['diagnostic_unfinished']} |")
lines+=['','## Measured algebra dimensions','', 'MAC/step columns are millions of real MACs for Gamma; complex cases are millions of complex MACs and have no real baseline. The small-solve columns are explicitly modeled MAC-equivalents. Total MAC ratio includes the change in diagnostic step count. These are work counts, not timings.','', '| Case | Profile | Capacity | H+S rectangle | Packed Ritz | Contraction | Small model | Current real MAC/step | Current small model | MAC saving % | Model-inclusive saving % | Worst solve MAC saving % | Total MAC ratio |','|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|']
for r in summary:
    a=f"| {r['case']} | {r['profile']} | {r['capacity']} | {r['rectangle_per_step']/1e6:.4f} | {r['rotation_per_step']/1e6:.4f} | {r['contraction_per_step']/1e6:.4f} | {r['small_model_per_step']/1e6:.4f} |"
    if 'baseline_steps' in r:a+=f" {r['baseline_mac_per_step']/1e6:.4f} | {r['baseline_small_model_per_step']/1e6:.4f} | {100*r['mac_saving_fraction']:.2f} | {100*r['model_inclusive_saving_fraction']:.2f} | {100*r['worst_solve_mac_saving_fraction']:.2f} | {r['total_mac_ratio']:.3f} |"
    else:a+=' n/a | n/a | n/a | n/a | n/a | n/a |'
    lines.append(a)
lines+=['','## SCF outcomes','', '| Case | Profile | Seed | Original SCF outcome |','|---|---|---:|---|']
for r in outcomes:lines.append(f"| {r['case']} | {r['profile']} | {r['seed']} | {r['status']} |")
(root/'tables.md').write_text('\n'.join(lines)+'\n')
print('groups',len(summary),'steps',sum(r['steps'] for r in summary),'original solves',sum(r['solves'] for r in summary)//3)
"###;
