# P10A RMM-DIIS benchmark status and protocol

Status: implemented experimental option; final measurements are recorded below.
Base: `0142ebf`, 2026-11-03. Default stays Davidson. See the
[design](../research/rmm-diis-design.md) and ADR 0020.

## Evidence boundaries and rejected approaches

The initial adapter ownership blocker was explicitly cleared by the orchestrator.
Reusable EVD, GEMM output and MGS output adapters now support the RMM loops.
Existing Davidson code, numerical adapter functions/tests and the QE comparison
script remain unchanged. No dependency or unsafe code was added.

Existing optimization-log profiles motivate investigating band-count scaling,
but are historical evidence, not fresh profiles of this worktree. The P7P7
outcomes argue against changing Davidson capacity or its convergence criterion;
this ticket leaves both alone. No performance default is adopted from RMM.
Optimization evidence is recorded here because the ticket excludes edits to the
shared optimization-log, baselines, README and artifacts.

An initial hybrid skipped warm-up whenever the starting block was orthonormal.
Engine measurements showed excessive RMM work with coordinate starts (Si Gamma:
2750 applications vs Davidson 470; four-cell Si: 46480 vs 2679). This variant
was rejected: the final policy screens true Rayleigh residuals and warms up
poor starts. These preliminary numbers are not mixed with final timings.

An initial benchmark attempted to generate moderate guesses using a separate
Davidson run at capacity m+4. At n=264, m=128, wide spectrum, that setup itself
failed after 800 iterations (max residual 162.89542282050255, dimensionless),
exit 101. This was not an RMM failure and is not hidden: the final harness uses
a reproducible perturbation of the dense eigenspace for the warm-start category,
without a separate iterative convergence dependency. Guess/oracle construction
is outside timings, identically for both solvers. This is controlled synthetic
initialization, not evidence that real SCF can obtain these guesses for free.

## Experiment matrix and measurement rules

Synthetic Hermitian operators must include both diagonal and complex rotated
operators, with n at least 2m and preferably 4m; m = 8, 16, 32, 64, 128, and
larger only if resources allow. Dense oracles are limited to n <= 1024. Use fixed
seeds, exact repeated eigenvalues, separations down to 1e-10, tightly clustered
groups separated by gaps, and wide spectra. Include a wanted boundary inside a
cluster. Do not benchmark only preconditioner-exact diagonal problems.

Compare identical Hamiltonians, guesses and requested tolerances for Davidson,
RMM alone and the hybrid. Initial guesses: random, moderately contaminated low
eigenspace, accurate prior eigenspace, deliberately duplicated columns and an
invariant block missing lower eigenstates. Record convergence failures alongside
successful timings; never remove failed runs from the crossover analysis.

The engine matrix must vary cubic, skewed, elongated and slab-like cells, atom
counts, Gamma and k meshes, projector multiplicities and angular channels, gaps
and band counts. Use existing synthetic NC pseudopotential fixtures plus the
read-only Si UPF at `/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF`.
Do not commit/copy the UPF. Si supercell measurements alone cannot adopt a default.
Synthetic fixed-H tests do not substitute for complete SCF runs.

For each case record:

- Configuration, seed, hardware, compiler/backend version, source revision,
  thread count and initial-guess quality; UPF SHA-256 for real systems.
- All vector-column operator applications, including trial directions, warm-up,
  certification and failed RMM work before Davidson fallback.
- Residual history and final per-band residuals in eigenvalue units, maximum
  entrywise Gram defect (dimensionless), dense eigenvalue error and complete
  cluster projector errors. Synthetic eigenvalue units are dimensionless.
- Wall time (s), peak resident memory (bytes), estimated persistent solver
  workspace separately, and phases taken. Process RSS includes backend/operator
  and runtime storage, not just the solver arrays.
- End-to-end energy difference in Ha, density and charge diagnostics, final
  SCF fixed-point certificate, gap and SCF iterations.

Use at least five alternating A/B release runs after build/warm-up, report every
sample and median/min/max. Use isolated processes for RSS. Thread counts 1, 2
and 4 must be compared for deterministic outputs; this finite test set must not
be described as proof for every possible thread count. Count arithmetic order
and independent k-point scheduling separately.

Candidate crossover: the first band range where the hybrid repeatedly wins on
complete SCF wall time at unchanged correctness, across several system families.
Report variability and subsequent reversals; there may be no single crossover.
A cheap synthetic operator can make orthogonalization dominant and predict a
crossover that disappears when Hamiltonian applications become expensive.

## Reproduction commands available now

These run the unchanged Davidson numerical baseline as well as workspace gates.
Long numerical tests remain ignored in the debug gate.

```sh
export PATH=/opt/homebrew/opt/rustup/bin:$PATH
export RUSTUP_HOME=/private/tmp/p0c-rustup
export CARGO_NET_OFFLINE=true
scripts/check.sh
cargo test --offline --release -p dft-solvers -- --ignored --nocapture
```

The synthetic driver is the ignored `benchmark` test in
`crates/dft-solvers/tests/rmm.rs`:

```sh
RMM_REPEATS=5 cargo test --offline --release -p dft-solvers --test rmm benchmark -- --ignored --nocapture
# Optional: RMM_BANDS=8,16; RMM_SOLVER=davidson or hybrid; RMM_HISTORY=1.
```

It uses complex Fourier-rotated spectra (linear, threefold degenerate, 1e-7
near-degenerate and quadratic/wide), n=2m+8, 8–128 bands, tolerance 1e-9,
800 Davidson iterations with capacity min(n,6m), and seed 2026. RMM controls are
history 4, at most 80 sweeps, 10 stale sweeps, warm-up residual 1e-3. All are
experimental settings, not adopted production defaults. Every successful row
checks eigenvalues against the dense oracle within 1e-8, true residuals within
the requested 1e-9, and Gram defect below 1e-10. Fallback reason is printed.
`RMM_HISTORY=1` prints complete per-sweep residual histories; the ordinary rows
report final residuals, orthogonality, error, operator count and RMM sweeps.

## Acceptance and recommendations

The final command results and measurements below supersede the design-only
ledger. RMM is an experimental opt-in, not a recommended replacement for
Davidson. A successful hybrid row may contain Davidson warm-up or fallback and
must not be labeled pure RMM convergence. Near-degenerate boundary stagnation is
reported explicitly and tested through recovery. No measured universal crossover
is assumed. Tests cover no false duplicate-vector convergence, per-band
certificates, output shape/finiteness failures, safe stagnation and deterministic
repeated solves. The engine retains its full fixed-point certificate.

## Measured unchanged Davidson baseline

The release command above completed successfully on this Apple Silicon worktree.
These are baseline correctness results, not RMM results or timing comparisons.
All numerical quantities in this table are dimensionless. The test executable
reported 4.93 s for the two tests together; that includes oracle and validation
work and is not a solver benchmark.

| Case | Maximum eigenvalue error vs dense | Maximum residual | Gram defect | Iterations |
| --- | ---: | ---: | ---: | ---: |
| N=200, 8 bands | 8.526512829121202e-14 | 8.877270171161885e-11 | 8.437695117484581e-15 | 92 |
| N=300, 4 bands, preconditioned | 1.173661168252238e-11 | 1.3656221937221573e-11 | 2.4426435631164226e-15 | 51 |
| N=300, 4 bands, unpreconditioned | 1.3728129744094986e-11 | 9.989769242729403e-10 | 3.885780632587054e-14 | 2215 |

N=200 relative operator Hermiticity defect: 1.7236847321817372e-15.
The existing assertions passed without edits or relaxed tolerances.

## Reproducible engine driver

The driver is embedded here to keep the engine source change limited to its
selection field and two call sites. It reuses the existing synthetic fixture
module without editing it. Build in a temporary standalone Cargo project using
only existing workspace crates, offline. Run from the worktree root:

```sh
python3 - <<'PYDRIVER'
from pathlib import Path
import re
root = Path.cwd()
d = Path('/private/tmp/p10a-driver')
(d/'src').mkdir(parents=True, exist_ok=True)
(d/'.cargo').mkdir(exist_ok=True)
(d/'Cargo.toml').write_text('[package]\nname="p10a-driver"\nversion="0.0.0"\nedition="2024"\n[dependencies]\n' + ''.join(f'{c}={{path="{root}/crates/{c}"}}\n' for c in ['dft-engine','dft-geometry','dft-physics','dft-upf','dft-solvers']))
(d/'.cargo/config.toml').write_text('[source.crates-io]\nreplace-with="vendored-sources"\n[source.vendored-sources]\ndirectory="'+str(root.parent/'rust-dft-vendor')+'"\n')
text=(root/'docs/performance/rmm-diis-benchmark.md').read_text()
code=re.search(r'```rust\n(.*?)\n```',text,re.S).group(1)
prefix=f'#[allow(dead_code)]\n#[path="{root}/crates/dft-engine/tests/common/mod.rs"]\nmod fixture;\n'
(d/'src/main.rs').write_text(prefix+code)
PYDRIVER
cargo build --offline --release --manifest-path /private/tmp/p10a-driver/Cargo.toml
RMM_REPEATS=3 /private/tmp/p10a-driver/target/release/p10a-driver synthetic
# Other modes: gamma, k222, k222s, k444s, supercell.
# RMM_CELLS=2, 4 or 8 selects 16, 32 or 64 bands in supercell mode.
```

The recorded runs use the repository's pinned Rust toolchain and
`RUSTUP_HOME=/private/tmp/p0c-rustup`. The harness writes reduced comparison JSON
to `/private/tmp/p10a-<mode>-<hybrid>.json`; these are inputs to the unchanged QE
comparison script, not new QE reference outputs. It asserts <=1e-9 Ha energy
agreement for every pair. Synthetic repeat 1 uses random rather than free
electron initialization for both algorithms, recording initial-guess variability.
Thread comparisons require exactly identical energies, densities, final
eigensolutions and SCF iteration counts for 0/1/2/3/4/8 threads.

```rust
use dft_engine::*;
use dft_geometry::{Atom,Cell,kpoints::{KPoint,monkhorst_pack,reduce_time_reversal}};
use dft_solvers::{SolverChoice,HybridRmm};
use std::time::Instant;
fn run(cell:&Cell,atoms:&[Atom],pp:&dft_upf::Pseudopotential,points:&[KPoint],bands:usize,cutoff:f64,hybrid:bool,threads:usize,free_guess:bool)->ScfResult {
    ScfProblem::with_kpoints(cell,atoms,&[SpeciesPotential{label:&atoms[0].species,potential:pp}],ScfOptions {
        eigensolver:if hybrid {SolverChoice::Hybrid(HybridRmm::default())} else {SolverChoice::Davidson},
        bands:Some(bands), cutoff_hartree:cutoff, threads,free_electron_initial_guess:free_guess,
        davidson_subspace:6*bands,eigen_tolerance:1e-9,density_tolerance:1e-8,energy_tolerance:1e-10,
        ..Default::default()
    },points,OccupationPolicy::Insulator{minimum_gap_hartree:1e-7}).unwrap().run().unwrap()
}
fn record(case:&str,bands:usize,solver:&str,rep:usize,seconds:f64,r:&ScfResult) {
    let last=r.history.last().unwrap();
    assert!(last.fixed_point_certificate.is_some());
    let calls:usize=r.history.iter().map(|h|h.operator_applications).sum();
    println!("SCF,{case},{bands},{solver},{rep},{seconds:.9},{:.17},{},{calls},{:e},{:e},{:?}",r.total_hartree(),r.history.len(),last.max_eigen_residual,last.charge_error,r.gap_hartree);
}
fn output(path:&str,r:&ScfResult) {
    let points=r.k_points.iter().map(|k|format!("{{\"fractional\":{:?},\"weight\":{},\"eigenvalues\":{{\"ev\":{:?}}}}}",k.point.fractional,k.point.weight,k.solution.eigenvalues.iter().map(|v|v*27.211386245988).collect::<Vec<_>>())).collect::<Vec<_>>().join(",");
    std::fs::write(path,format!("{{\"total_energy\":{{\"ha\":{}}},\"kpoints\":{{\"count\":{},\"points\":[{}]}}}}",r.total_hartree(),r.k_points.len(),points)).unwrap();
}
fn main(){
    let mode=std::env::args().nth(1).unwrap_or("synthetic".into());
    let repeats:usize=std::env::var("RMM_REPEATS").unwrap_or("3".into()).parse().unwrap();
    if mode=="synthetic" {
        for (index,rows) in [ [[4.,0.,0.],[0.,4.,0.],[0.,0.,4.]], [[4.,0.,0.],[0.6,4.2,0.],[0.2,0.3,4.5]], [[4.,0.,0.],[0.,4.,0.],[0.,0.,7.]], [[5.,0.,0.],[0.,5.,0.],[0.,0.,3.5]] ].into_iter().enumerate() {
            let cell=Cell::from_user_rows_bohr(rows).unwrap();
            let pp=fixture::potential(if index==2 {4.0}else{2.0},1.0,index%2==1,401);
            let atoms=[Atom::new("Xx",[0.13,0.21,0.32])];
            let nk=[1,2,3,5][index];
            let points:Vec<_>=(0..nk).map(|i|KPoint{fractional:[0.025*i as f64,0.015*i as f64,0.0],weight:1.0/nk as f64}).collect();
            let bands=if index==0 {8} else {16};
            let mut energies=[0.0;2];
            let mut reference=None;
            for rep in 0..repeats {for hybrid in if rep%2==0 {[false,true]}else{[true,false]} {
                let start=Instant::now();let r=run(&cell,&atoms,&pp,&points,bands,5.0,hybrid,1,rep%2==0);
                record(&format!("synthetic-{index}"),bands,if hybrid{"hybrid"}else{"davidson"},rep,start.elapsed().as_secs_f64(),&r);
                energies[usize::from(hybrid)]=r.total_hartree();
                if hybrid && rep==0 {reference=Some(r);}
            } assert!((energies[0]-energies[1]).abs()<1e-9);}
            let reference=reference.unwrap();
            for threads in [0,2,3,4,8] {
                let r=run(&cell,&atoms,&pp,&points,bands,5.0,true,threads,true);
                assert_eq!(reference.total_hartree().to_bits(),r.total_hartree().to_bits());
                assert_eq!(reference.density,r.density);
                assert_eq!(reference.history.len(),r.history.len());
                for (a,b) in reference.k_points.iter().zip(&r.k_points){assert_eq!(a.solution,b.solution);}
            }
            println!("DETERMINISTIC,synthetic-{index},threads=0/1/2/3/4/8");
        }
    } else {
        let pp=dft_upf::parse_upf(&std::fs::read_to_string("/Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF").unwrap()).unwrap();
        let base=fixture::cell();
        let points=match mode.as_str(){"k222"=>monkhorst_pack([2;3],[0.0;3]).unwrap(),"k222s"=>reduce_time_reversal(&monkhorst_pack([2;3],[0.5;3]).unwrap()).unwrap(),"k444s"=>reduce_time_reversal(&monkhorst_pack([4;3],[0.5;3]).unwrap()).unwrap(),_=>vec![KPoint{fractional:[0.0;3],weight:1.0}]};
        let count=if mode=="supercell" {std::env::var("RMM_CELLS").unwrap_or("2".into()).parse::<usize>().unwrap()}else{1};
        let a=10.261_212_856_717_933/2.0;
        let cell=if count==1 {base}else{Cell::from_user_rows_bohr([[0.,a*count as f64,a*count as f64],[a,0.,a],[a,a,0.]]).unwrap()};
        let atoms:Vec<_>=(0..count).flat_map(|i|[Atom::new("Si",[i as f64/count as f64,0.,0.]),Atom::new("Si",[(i as f64+0.25)/count as f64,0.25,0.25])]).collect();
        let bands=8*count;
        let mut energies=[0.0;2];
        for rep in 0..repeats{for hybrid in if rep%2==0 {[false,true]}else{[true,false]}{
            let start=Instant::now();let r=run(&cell,&atoms,&pp,&points,bands,15.,hybrid,1,true);
            record(&mode,bands,if hybrid{"hybrid"}else{"davidson"},rep,start.elapsed().as_secs_f64(),&r);
            energies[usize::from(hybrid)]=r.total_hartree();
            output(&format!("/private/tmp/p10a-{mode}-{hybrid}.json"),&r);
        }assert!((energies[0]-energies[1]).abs()<1e-9,"energies {energies:?}");}
    }
}

```

## Final synthetic results (shared Apple Silicon Mac, release)

Rust 1.99.1 (58cf8f9c3, 2026-10-13), arm64, faer 0.24.4. The sandbox denies
`sysctl machdep.cpu.brand_string`, so the exact chip model was not independently
queried. Timings are indicative: other tasks and builds share this Mac; these
are not quiet-dev08 speedup claims. Five alternating repetitions per cell.
All 400 successful rows passed the dense/residual/orthogonality assertions.
Maximum observed eigenvalue error was 4.911271389573812e-10, true residual
9.975623882843934e-10, and Gram defect 5.173639295056041e-14 (dimensionless).

There is **no monotone crossover** from 8 through 128 bands. The wide-spectrum
cold cases at 16 and 64 bands favor the hybrid, but the 128-band version reverses
that result. One degenerate warm 64-band case has a small timing advantage despite
more applications; this is not convincing evidence of a general improvement.
Keep Davidson as default; switch only for measured fixed-system workloads with
these same correctness checks. There is no band-count-only switching rule.

Table entries are median seconds [min, max], then deterministic application counts
for Davidson/hybrid. Spectra: 0 linear, 1 threefold degenerate, 2 near-degenerate
(1e-7 spacing), 3 quadratic. Synthetic energies and residuals are dimensionless.
The fallback column records the hybrid result (identical across repetitions).

| Spectrum | Bands | Guess | Davidson s [min,max] | Hybrid s [min,max] | H columns D/H | Fallback |
| --- | ---: | --- | --- | --- | --- | --- |
| 0 | 8 | cold | 8.9208e-05 [8.6625e-05, 0.000213333] | 0.000105667 [0.000103709, 0.000119417] | 32/48 | None |
| 0 | 8 | warm | 4.1292e-05 [4.1166e-05, 4.3667e-05] | 9.9333e-05 [9.8291e-05, 0.000103208] | 24/96 | None |
| 1 | 8 | cold | 8.175e-05 [8.1583e-05, 8.4833e-05] | 9.8084e-05 [9.6542e-05, 0.000100417] | 32/48 | None |
| 1 | 8 | warm | 5.85e-05 [5.4375e-05, 0.000254] | 9.825e-05 [9.7708e-05, 0.000102417] | 25/96 | None |
| 2 | 8 | cold | 8.4416e-05 [8.1584e-05, 9.1708e-05] | 9.9417e-05 [9.6834e-05, 0.000100958] | 32/48 | None |
| 2 | 8 | warm | 3.9708e-05 [3.8833e-05, 4.4791e-05] | 0.000101416 [9.7041e-05, 0.000104041] | 24/96 | None |
| 3 | 8 | cold | 8.425e-05 [8.3875e-05, 8.85e-05] | 0.000100875 [9.9916e-05, 0.000103666] | 32/48 | None |
| 3 | 8 | warm | 4.0875e-05 [4.0625e-05, 4.1292e-05] | 9.8333e-05 [9.7584e-05, 9.875e-05] | 24/96 | None |
| 0 | 16 | cold | 0.00213029 [0.00201467, 0.00245375] | 0.00592325 [0.00573783, 0.00783679] | 248/1889 | None |
| 0 | 16 | warm | 0.000216875 [0.000211958, 0.000221] | 0.000494542 [0.000470792, 0.000513125] | 48/192 | None |
| 1 | 16 | cold | 0.00289908 [0.00269804, 0.00362567] | 0.003498 [0.00324983, 0.00382171] | 305/1016 | None |
| 1 | 16 | warm | 0.000355041 [0.00035275, 0.000363042] | 0.000480416 [0.000455667, 0.000484125] | 50/192 | None |
| 2 | 16 | cold | 0.00298125 [0.00296037, 0.00346708] | 0.0069765 [0.00680021, 0.0089795] | 327/2216 | Some(Stagnation) |
| 2 | 16 | warm | 0.000266042 [0.000265333, 0.00026925] | 0.000638584 [0.000638167, 0.000641125] | 48/192 | None |
| 3 | 16 | cold | 0.00572479 [0.00489675, 0.00641725] | 0.00182088 [0.00181463, 0.00183167] | 510/238 | None |
| 3 | 16 | warm | 0.000200792 [0.000199708, 0.000201458] | 0.000457625 [0.000456, 0.000459208] | 48/192 | None |
| 0 | 32 | cold | 0.0184805 [0.0154711, 0.0199278] | 0.0416813 [0.0366792, 0.0449933] | 575/4352 | None |
| 0 | 32 | warm | 0.00124463 [0.00117967, 0.0016725] | 0.00256708 [0.00251217, 0.00298458] | 96/384 | None |
| 1 | 32 | cold | 0.0116068 [0.0103361, 0.012835] | 0.0239793 [0.0229661, 0.024261] | 451/2573 | None |
| 1 | 32 | warm | 0.00211892 [0.00193462, 0.00291938] | 0.00274262 [0.00252792, 0.00296625] | 99/384 | None |
| 2 | 32 | cold | 0.027169 [0.0268725, 0.0274619] | 0.0465606 [0.0455509, 0.0480414] | 644/5261 | Some(Stagnation) |
| 2 | 32 | warm | 0.00126637 [0.00117483, 0.00158171] | 0.00254671 [0.00251675, 0.00312558] | 96/384 | None |
| 3 | 32 | cold | 0.0243206 [0.0219008, 0.0243936] | 0.0471422 [0.045734, 0.0486937] | 810/4548 | None |
| 3 | 32 | warm | 0.00126687 [0.00118533, 0.00161787] | 0.00261471 [0.00253979, 0.00317883] | 96/384 | None |
| 0 | 64 | cold | 0.0897013 [0.0869895, 0.0906174] | 0.328534 [0.322909, 0.342899] | 950/11545 | None |
| 0 | 64 | warm | 0.00900938 [0.00824612, 0.010587] | 0.0188221 [0.0174813, 0.0201008] | 192/768 | None |
| 1 | 64 | cold | 0.094167 [0.0862883, 0.0960044] | 0.210289 [0.202558, 0.211899] | 974/6936 | None |
| 1 | 64 | warm | 0.0182651 [0.0180516, 0.020076] | 0.0173056 [0.0171673, 0.0177144] | 196/768 | None |
| 2 | 64 | cold | 0.10287 [0.0902597, 0.118753] | 0.42325 [0.356543, 0.431982] | 989/13414 | Some(Stagnation) |
| 2 | 64 | warm | 0.00980796 [0.00828237, 0.0119956] | 0.0209378 [0.0204406, 0.0225216] | 192/768 | None |
| 3 | 64 | cold | 0.220465 [0.207273, 0.237708] | 0.0941519 [0.0900942, 0.0983433] | 1539/931 | None |
| 3 | 64 | warm | 0.0140155 [0.0125736, 0.0232012] | 0.0415162 [0.025639, 0.107397] | 192/768 | None |
| 0 | 128 | cold | 1.23794 [0.963929, 2.15594] | 4.29739 [3.96433, 4.46057] | 2048/42765 | None |
| 0 | 128 | warm | 0.0665571 [0.0641268, 0.0695142] | 0.143319 [0.142005, 0.146049] | 384/1536 | None |
| 1 | 128 | cold | 0.709257 [0.701586, 0.780859] | 1.93337 [1.87099, 2.08834] | 1811/19617 | None |
| 1 | 128 | warm | 0.0929042 [0.092313, 0.0938717] | 0.126878 [0.126452, 0.127381] | 388/1536 | None |
| 2 | 128 | cold | 1.11623 [1.10669, 1.15673] | 2.77983 [2.77279, 3.0146] | 2092/29142 | Some(Stagnation) |
| 2 | 128 | warm | 0.0742877 [0.0704583, 0.078549] | 0.158918 [0.146514, 0.211537] | 384/1536 | None |
| 3 | 128 | cold | 1.03403 [0.97809, 1.17786] | 7.89255 [7.51078, 8.39859] | 2643/79163 | None |
| 3 | 128 | warm | 0.110008 [0.069888, 0.133998] | 0.19866 [0.17132, 0.514069] | 384/1536 | None |

### Pure RMM and memory

A separate pure-RMM warm-start sweep at 128 bands converged all four spectra in
two outer sweeps and 1408 applications (no Davidson, no fallback). The maximum
eigenvalue error was 1.0913936421275139e-10 and maximum residual
1.216434539976613e-10. This uses the accurately contaminated dense guess described
above; it is not a cold-start capability claim. Near-degenerate partial-boundary
stagnation with a less accurate guess remains covered by the unit test.

Peak RSS was measured with `os.wait4(child_pid, 0).ru_maxrss` (bytes on Darwin),
with a separate process for each solver, `RMM_BANDS=128 RMM_REPEATS=1` and
`RMM_SOLVER` selecting the solver. These include fixture/oracle/backend/runtime
memory, not just solver storage. Davidson/hybrid processes include both cold and
warm cases; pure RMM includes warm only, so this is not a matched memory-speedup
claim. All three child exit codes were 0.

| Process | Peak RSS (bytes) |
| --- | ---: |
| davidson | 70631424 |
| hybrid | 77070336 |
| rmm | 28721152 |

The explicit RMM complex arrays occupy approximately
`16 * ((9 + 2h)nm + 2m² + 2h² + n)` bytes, plus scalar histories and backend EVD
scratch; h is the history size. For n=264, m=128, h=4 this is 9,720,448 bytes
before scalar arrays/backend scratch. History storage is O(nmh), not free, and
hybrid peak memory can exceed Davidson's. No general memory reduction is claimed.

`/usr/bin/time -l` ran the pure-RMM test successfully but itself exited 1 after
`sysctl kern.clockrate` was denied; its incomplete output was not used for RSS.
`sample <owned-driver-pid> 1` exited 255 because the sandbox could not examine
the child. No fresh sampling profile or quiet-node performance claim is made.
The driver itself completed that profiled-attempt run with exit 0.

Example measured pure-RMM maximum-residual histories for 8-band warm starts:

```text
HISTORY,0,8,warm,true,[1.3888888889439229e-5, 7.86297319844503e-15]
HISTORY,1,8,warm,true,[5.208333333086196e-6, 3.890120528938143e-14]
HISTORY,2,8,warm,true,[5.2083331601912725e-6, 2.683549374998111e-15]
HISTORY,3,8,warm,true,[0.0003055555555466844, 1.008622800532327e-13]
```

## Reference comparison and baseline preservation

The unchanged `scripts/compare_kpoints_qe.py` passed all eight combinations below
(exit 0 each), using freshly computed driver output and the existing QE 7.6
bundles. No new QE run was performed. QE provenance, inputs and commands remain
in their existing reference directories; the real UPF SHA-256 is
`e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217`, matching the
stored k-point metadata. The UPF was read in place and not copied.

| Case | Solver | Energy difference vs QE (Ha/atom) | Max eigenvalue difference (Ha) | Weight difference |
| --- | --- | ---: | ---: | ---: |
| gamma | Davidson | 1.280894553535e-06 | 1.270854584781e-06 | 0 |
| gamma | hybrid | 1.280894549538e-06 | 1.269410814410e-06 | 0 |
| k222 | Davidson | 1.281439412359e-06 | 1.987271588685e-06 | 0 |
| k222 | hybrid | 1.281439443002e-06 | 1.986245993844e-06 | 0 |
| k222s | Davidson | 1.281265468389e-06 | 1.773158624518e-06 | 0 |
| k222s | hybrid | 1.281265459507e-06 | 1.773070007204e-06 | 0 |
| k444s | Davidson | 1.281675959586e-06 | 2.100188964359e-06 | 0 |
| k444s | hybrid | 1.281675957809e-06 | 2.100285347795e-06 | 0 |

Comparison thresholds were unchanged: 1e-4 Ha/atom energy, 1e-3 Ha eigenvalues,
and 1e-6 weights. These QE tolerances are separate from the much tighter
**1e-9 Ha hybrid-versus-Davidson SCF gate**. Reproduce for each emitted JSON:

```sh
python3 -I scripts/compare_kpoints_qe.py --ours /private/tmp/p10a-gamma-false.json --qe-out tests/reference/si30_ours/qe.out
python3 -I scripts/compare_kpoints_qe.py --ours /private/tmp/p10a-k222-false.json --qe-out tests/reference/kpoints/k222_c/qe.out
python3 -I scripts/compare_kpoints_qe.py --ours /private/tmp/p10a-k222s-false.json --qe-out tests/reference/kpoints/k222s_c/qe.out
python3 -I scripts/compare_kpoints_qe.py --ours /private/tmp/p10a-k444s-false.json --qe-summary tests/reference/kpoints/k444s_c/qe_summary.json
# Repeat with true.json for the hybrid runs.
```

Default-solver equivalence command (exit 0):

```sh
python3 -I scripts/bench_scf.py --binary target/release/dft-cli --upf /Users/mino/Develop/rust-dft-wt/pseudo/Si.pz-vbc.UPF --cases gamma,k222s,k222,k444s --repeats 1 --baseline docs/performance/baseline.json
```

Differences versus the original pre-optimization baseline: Gamma 1.55e-13 Ha,
k222s 3.55e-15 Ha, k222 6.84e-14 Ha, all below 1e-12 Ha. That historical baseline
has no k444s entry; the tool ran it without a baseline comparison. Density hashes
differ from that old baseline (other already-merged algorithms changed arithmetic);
no bit-identity claim or new default-solver speedup is inferred. The current
change leaves the default solver implementation and its tests unchanged.

## Validation commands and scope

- `scripts/check.sh`: build, test, fmt and clippy all pass. Use
  `CARGO_NET_OFFLINE=true RUSTUP_HOME=/private/tmp/p0c-rustup` on this host.
- `cargo test --offline --release -p dft-solvers --test rmm -- --nocapture`:
  8 passed, 1 ignored benchmark, including explicit duplicate-vector fallback,
  strict/loose buffer certificate and terminal warm-up/fallback diagnostics.
- `cargo test --offline --release -p dft-solvers --lib -- --ignored --nocapture`:
  existing two Davidson numerical tests pass; their numbers remain unchanged.
- The 400-row synthetic sweep, pure-RMM sweep and engine driver assert every
  numerical gate before reporting their successful measurements.
- Existing `linalg.rs` and `linalg/tests.rs` content was checked byte-for-byte
  against `git show HEAD:<path>`: both original files are unchanged prefixes of
  the new files. Only appended functions/tests were added.

No timing or numerical data here is written from a solver hot loop. Benchmark
stdout, JSON comparison files and timers are in the test/driver layer. Raw local
logs are in `/private/tmp/p10a-*.log`; the tables, histories and executable driver
above preserve the reviewable evidence without editing forbidden artifacts.

Adapter additions are `gemm_into`, `HermitianWorkspace::{new,solve}` and
`orthonormalize_into`; the appended adapter tests are
`reusable_evd_gemm_and_rank_filter` and
`reusable_adapters_rectangular_and_empty_shapes`. The engine change consists of the selection
field/default and the two solve dispatches (seven added lines, two replaced).
The new RMM module and tests are isolated; `lib.rs` adds only module/re-exports.
Changes remain uncommitted for the orchestrator's review.

## Complete SCF results

The driver compares complete setup plus SCF, with density tolerance 1e-8,
energy tolerance 1e-10 Ha and eigen tolerance 1e-9 Ha for both choices. These are
tighter than the default benchmark CLI settings, so do not compare their wall
times directly to the historical baseline timings. Three alternating repetitions
per pair; time entries are median [minimum, maximum] seconds. Every run passed
the engine fixed-point certificate; every pair agreed within 1e-9 Ha.

Synthetic cells are cubic (4,4,4 Bohr), skewed, elongated (4,4,7 Bohr), and flattened
(5,5,3.5 Bohr). They use valence 2 or 4, zero projectors or three projectors with
two s channels and one p channel, and 1/2/3/5 weighted k points. Initial guesses
alternate free-electron and seeded random starts. All four synthetic cases
produce bit-identical energies, densities, final eigensolutions and iteration
counts across thread requests 0/1/2/3/4/8 (workers capped at the k-point count).
The Si supercells repeat the first primitive vector 2/4/8 times: 4/8/16 atoms,
16/32/64 bands. No new real-element physical reference is claimed for synthetics.

The H-column column below is explicitly the **sum over ordinary SCF iteration
solves**, as emitted by `IterationRecord::operator_applications`; it includes
hybrid screening/warm-up/refinement/fallback and eigensolver residual certification
inside those solves. The separate final SCF fixed-point map's applications are
not included in this column. Whole-run times include that map. The standalone
solver benchmark counts all applications including its certification. This scope
must be retained when comparing the engine counts to other tools.

| Case | Bands | Davidson seconds [min,max] | Hybrid seconds [min,max] | Max pair dE (Ha) | Iteration H columns D/H |
| --- | ---: | --- | --- | ---: | --- |
| synthetic-0 | 8 | 0.00479521 [0.00458713, 0.00595604] | 0.00943396 [0.00902321, 0.0128742] | 8.881784e-16 | 272/319 ; 1077/1877 |
| synthetic-1 | 16 | 0.017031 [0.0160792, 0.0171546] | 0.175782 [0.172052, 0.234865] | 2.664535e-15 | 1156/1192 ; 28309/32307 |
| synthetic-2 | 16 | 0.0491024 [0.0466597, 0.0518571] | 0.945448 [0.923243, 1.00494] | 1.065814e-14 | 2269/2379 ; 87604/99807 |
| synthetic-3 | 16 | 0.0542549 [0.0523834, 0.0543345] | 0.717378 [0.711964, 0.741388] | 8.881784e-16 | 3138/3323 ; 90277/93408 |
| gamma | 8 | 0.207651 [0.193776, 0.207854] | 0.374344 [0.352888, 0.409849] | 7.993606e-15 | 470 ; 1219 |
| k222 | 8 | 1.61985 [1.50527, 1.94325] | 4.92797 [3.86948, 5.15265] | 6.128431e-14 | 4317 ; 14982 |
| k222s | 8 | 0.822105 [0.773646, 0.829058] | 2.29401 [2.26772, 2.31628] | 1.776357e-14 | 2308 ; 9721 |
| k444s | 8 | 4.83729 [4.76701, 5.67217] | 28.8532 [25.1885, 32.6134] | 3.552714e-15 | 17628 ; 136022 |
| super2 | 16 | 0.842601 [0.767645, 0.848667] | 2.68118 [2.13158, 3.30887] | 4.796163e-14 | 1086 ; 4480 |
| super4 | 32 | 11.8954 [9.6965, 12.9466] | 32.5729 [31.2413, 46.9605] | 1.385558e-13 | 2679 ; 33073 |
| super8 | 64 | 61.6889 [52.8596, 69.3041] | 229.351 [227.564, 242.809] | 5.897505e-13 | 6099 ; 105331 |

Maximum hybrid/Davidson energy difference: 5.89750470680883154e-13 Ha. Maximum reported
SCF iteration eigen residual: 5.46865691280404891e-10 Ha; charge error: 5.55111512312578270e-15
(relative). Sampled insulating gaps span 0.0172138258 to 1.15839147 Ha.

Every measured complete SCF workload favors Davidson. Even the 64-band Si case
shows no end-to-end crossover. The synthetic fixed-H wins at selected sizes do
not transfer to these SCF cases. Keep the module available for larger or different
workloads, but do not recommend enabling it for these systems. Quiet dev08
remeasurement should use the same driver and data (including failure and warm-up
counts) before any performance/default decision; no remote change was made here.

Timing runs used the final numerical algorithm. Subsequent changes only enriched
terminal-error diagnostics and handled zero-dimensional workspace construction,
which these successful positive-dimensional runs do not exercise. Final workspace
gates and the RMM acceptance tests cover those changes as well.

The design note distinguishes this affine-coefficient residual-minimization
variant from the normalized-history formulation in Kresse–Furthmüller (1996).
Negative performance here is not a general result about every RMM-DIIS variant.
