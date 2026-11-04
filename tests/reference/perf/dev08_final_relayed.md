# dev08 final measurement, relayed summary (2026-11-04)

Provenance: measured on dev08 by the remote measurement session (not by the Claude Code session that holds this repository), relayed to this repository as plain text in cross-session messages. The full source file `dev08_final.json` (63176 B, sha256 `8fb74fb0115b00d73f388e3a1663dca655253683c0e414d4802da0b62278ed53`) stayed on dev08 (`~/rustdft_p6/perf7/dev08_final.json`); attempts to transfer it as base64 chunks failed because of character corruption on the relay path (chunk sha256 mismatch), so it has not been verified here. Numbers below were copied from the messages without rounding. When the original file arrives, compare and correct this document. The architect did not re-run anything on dev08.

## Part 1: host, builds, gates, BAANDI wall times

```
HOST
 cpu     Intel(R) Xeon(R) Gold 6230 CPU @ 2.10GHz
 topology CPU(s): 80  Core(s) per socket: 20  Socket(s): 2
 os      Rocky Linux release 9.7 (Blue Onyx)
 kernel  5.14.0-611.24.1.el9_7.x86_64
 stage1  before 00:10:18 up 270 days, load 0.28 1.84 1.90
         after  00:13:42 up 270 days, load 2.21 1.79 1.85
 stage2  before 00:16:44 up 270 days, load 0.25 1.06 1.55
         after  00:21:30 up 270 days, load 1.69 1.54 1.62
 qe      before 00:27:13 up 270 days, load 0.28 1.06 1.39
         after  00:28:09 up 270 days, load 1.90 1.41 1.50
 qe_redo before 00:30:35 up 270 days, load 0.17 0.86 1.27
         after  00:31:17 up 270 days, load 1.15 1.02 1.31

CODE
 commit   fe0492b8b8865e8029214861f666dfc816045800
 rustc    rustc 1.99.1 (58cf8f9c3 2026-10-13)
 mklroot  /opt/intel/oneapi/mkl/2026.1
 sha pure 630e058c8fd0188882bad8257490c8e842df93ff209d96200d829336f8049e6a
 sha mkl  07c49b490a8330554dfeb7370b71f8ec116be9bd850936a699f9eb2c9699d2a6
 pinning  taskset -c 0-(T-1), one socket; MKL_NUM_THREADS=1 (sequential layer)

QE
 7.6  ~/fq/builds/orig-mpi-dfti/bin/pw.x (MKL 2026.1)
 env  OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
 run  mpirun -np N --map-by core --bind-to core
 grid ours (make_qe_input.py --fft-grid ours)
 upf  e8d933754cd51c6bb4b2a809151f89e0647e53d878bab88d26e1b5a5d68d5217

GATE energy pure vs mkl (tol 1e-12 Ha)
 gamma  0.0                      pass
 k222s  5.329070518200751e-15    pass
 k222   7.993605777301127e-15    pass
 k444s  0.0                      pass
 si8    4.583000645652646e-13    pass

GATE density sha256 identical across threads 1,4,8,16
 pure   gamma T  k222s T  k222 T  k444s T  si8 T
 mkl    gamma T  k222s T  k222 T  k444s T  si8 T
 across backends: all False (gamma, k222s, k222, k444s, si8)

QE spread >5% remaining after 9 repeats (6 of 80)
 gamma conv_thr 1e-12 np1
 k222  conv_thr 1e-8  np8
 k444s conv_thr 1e-12 np16
 k444s conv_thr 1e-8  np4
 k444s conv_thr 1e-8  np8
 si8   conv_thr 1e-12 np8

BAANDI wall median s   case T pure mkl
 gamma  1  0.34542716294527054   0.30094488710165024
 gamma  4  0.3121382109820843    0.29170583561062813
 gamma  8  0.3244693912565708    0.31313569098711014
 gamma 16  0.361531101167202     0.3395419120788574
 k222s  1  1.0033023953437805    0.8073595240712166
 k222s  4  0.4279301092028618    0.36500227078795433
 k222s  8  0.42478300631046295   0.37121379002928734
 k222s 16  0.4577600732445717    0.40640338510274887
 k222   1  1.8497117385268211    1.452968668192625
 k222   4  0.6710673943161964    0.5747056677937508
 k222   8  0.4982187934219837    0.4237099513411522
 k222  16  0.528404988348484     0.44845157489180565
 k444s  1  7.115919779986143     5.493066567927599
 k444s  4  2.2357924543321133    1.7936346009373665
 k444s  8  1.4030840508639812    1.1598914787173271
 k444s 16  1.0025044530630112    0.8527095355093479
 si8    1  9.784100040793419     6.334394369274378
 si8    4  9.024114556610584     5.571705963462591
 si8    8  9.10754743590951      5.765270806849003
 si8   16  9.401152145117521     5.947770871222019
```

Iterations: gamma, k222s, k222, k444s: 10 for both backends; si8: 16 (pure) and 15 (mkl). BAANDI 3-run spread was within 5 percent for all entries (largest 4.9 percent, k444s pure T4).

## Part 2: QE 7.6 comparison

`qe_int_s` is QE's internal WALL average (mpirun start-up excluded); the ratios are BAANDI over QE, so a value above 1 means BAANDI is slower. Entries with the 9-repeat re-measurement use its median. QE internal timer resolution is 0.01 s: for gamma and k222s (0.05 to 0.27 s) one tick is up to 20 percent, so only k444s and si8 are valid comparisons (the relaying session's and the architect's judgement).

```
conv_thr 1e-8    case  N  qe_int_s  mkl/qe  pure/qe
 gamma   1   0.09     3.344    3.838
 gamma   4   0.06     4.862    5.202
 gamma   8   0.05     6.263    6.489
 gamma  16   0.06     5.659    6.026
 k222s   1   0.2      4.037    5.017
 k222s   4   0.11     3.318    3.89
 k222s   8   0.1      3.712    4.248
 k222s  16   0.09     4.516    5.086
 k222    1   0.34     4.273    5.44
 k222    4   0.17     3.381    3.947
 k222    8   0.16     2.648    3.114
 k222   16   0.13     3.45     4.065
 k444s   1   1.19     4.616    5.98
 k444s   4   0.54     3.322    4.14
 k444s   8   0.53     2.188    2.647
 k444s  16   0.39     2.186    2.571
 si8     1   0.68     9.315   14.388
 si8     4   0.27    20.636   33.423
 si8     8   0.21    27.454   43.369
 si8    16   0.23    25.86    40.875

conv_thr 1e-12   case  N  qe_int_s  mkl/qe  pure/qe
 gamma   1   0.12     2.508    2.879
 gamma   4   0.07     4.167    4.459
 gamma   8   0.06     5.219    5.408
 gamma  16   0.07     4.851    5.165
 k222s   1   0.27     2.99     3.716
 k222s   4   0.14     2.607    3.057
 k222s   8   0.13     2.855    3.268
 k222s  16   0.11     3.695    4.161
 k222    1   0.47     3.091    3.936
 k222    4   0.2267   2.535    2.96
 k222    8   0.2133   1.986    2.336
 k222   16   0.17     2.638    3.108
 k444s   1   1.7033   3.225    4.178
 k444s   4   0.76     2.36     2.942
 k444s   8   0.7367   1.574    1.905
 k444s  16   0.54     1.579    1.856
 si8     1   0.91     6.961   10.752
 si8     4   0.3333  16.717   27.075
 si8     8   0.29    19.88    31.405
 si8    16   0.3     19.826   31.337
```

QE total energy (Ry): gamma -14.59106606, k222s -15.83885195, k222 -15.66246524, k444s -15.85195738, si8 -62.64362757. conv_thr 1e-8 and 1e-12 give identical energies for all cases (8 printed digits); 1e-12 only adds iterations. QE iterations (conv_thr 1e-8 to 1e-12): gamma 7 to 11, k222s 6 to 9, k222 6 to 9, k444s 5 to 8, si8 14 to 20 (np1).

## Part 3: P11J three-mode overhead (relayed from dev08)

Columns: kind, budget, mode, target percent, median overhead percent, 95 percent bootstrap CI low and high, status. All five runs: `status=measured`, `all_identity_checks_pass=true`, `exclusions=[]`, 3 warm-ups; targets are classic 0.5, observed 1.0, agent 1.0.

Builds. Runs a, b, c and p100: candidate commit `fe0492b8b8865e8029214861f666dfc816045800`, candidate binary sha256 `05a37652db0d851930070c8fddebc8540fa050261ada12ae1f5d5452b9411aac`. Run p300: candidate commit `e10cfa3a11fcc4aac78a435d9afd9ed18e334c8a`, candidate binary sha256 `592efc7f9525343fa536c22b98d7aa482c2c937436ab95640e4503554cbb149b`. Reference (all five): sha256 `f4d6b054b7ac443fb300e9f30e3162b3a82365d766dc1b7d69b2d64189adbe38`, revision `9495d9fbaaf0e0bf6bf141ea61675e1d803eeeec`.

Why the candidate hash differs between fe0492b and e10cfa3 although `crates/`, the workspace Cargo.toml, `scripts/bench_runtime.py` and the bench driver are unchanged between them (architect finding): `crates/dft-cli/build.rs` embeds the git commit (`DFT_GIT_COMMIT`, written to `results.json` as `tool.git_commit`) into the binary, so any new commit gives a different binary hash with identical code. The four recorded source hashes (Cargo.lock, bench_driver.rs, dft-engine/tests/common/mod.rs, bench_runtime.py) are identical for p100 and p300, which is consistent with this. The five runs therefore measure the same code.

Load at start and end (1, 5, 15 minute): a [0.04, 0.41, 0.92] to [4.72, 2.52, 1.67]; b [0.27, 1.42, 1.38] to [0.52, 1.38, 1.37]; c [0.27, 1.21, 1.31] to [0.56, 1.19, 1.3]; p100 [0.28, 1.04, 1.24] to [0.95, 1.07, 1.23]; p300 [0.2, 0.51, 0.91] to [1.0, 0.89, 0.97].

p100 (100 pairs):
```
kind budget mode     target median                 ci_lo                  ci_hi                 status
0 1   agent    1.0  0.2877624812578694   -0.3071088319988813   0.5750247083821103  met
0 1   classic  0.5 -0.09923764640045185  -0.2808621182230797   0.21426464515368338 met
0 1   observed 1.0  0.24466471211861407  -0.12643158256590015  0.4831289285755602  met
0 200 agent    1.0  0.40910828237816776  -0.12138989263300037  0.7142055025266769  met
0 200 classic  0.5 -0.13172103031531823  -0.482907559607898    0.1719983874921449  met
0 200 observed 1.0  0.7057393229685571    0.39366332440960994  1.095080987537922   inconclusive
1 1   agent    1.0  0.6522224910375352    0.10218385561955867  0.9764623081810586  met
1 1   classic  0.5  0.12352132141866257 -0.48161868183403955  0.5157830910465511  inconclusive
1 1   observed 1.0  0.5692530435910559    0.07179795912866105  0.8244330163990732  met
1 200 agent    1.0  0.5523407901072308    0.280370755672088    0.8437162990057478  met
1 200 classic  0.5  0.16467021832186068  -0.006526051060035121 0.3262155630206709 met
1 200 observed 1.0  0.646859194597238     0.4621824263997154   0.9451939265883391  met
2 1   agent    1.0  0.18434321591482838 -0.042618178243658544 0.5486098319361421  met
2 1   classic  0.5 -0.43086921731666505 -0.7328865912903604  -0.10063620736803225 met
2 1   observed 1.0  0.28440085650422997 -0.04201564347652287  0.6941670829762825  met
2 200 agent    1.0  0.18030981092059806 -0.06931182171123762  0.6817674314168576  met
2 200 classic  0.5  0.2911636497672787   -7.209334719204818e-05 0.5287365588602211 inconclusive
2 200 observed 1.0  0.5290153209548953    0.09367619828819507  0.7981775269542246  met
```

p300 (300 pairs):
```
kind budget mode     target median                 ci_lo                  ci_hi                 status
0 1   agent    1.0  0.011829617946346183 -0.1996581943103981   0.30600430165304005 met
0 1   classic  0.5  0.2141466594521857    0.02837888898559182  0.44120629061665806 met
0 1   observed 1.0  0.17024189237898613  -0.10553210928053725  0.46533634895264964 met
0 200 agent    1.0  0.20063109175513238   0.04715605444655102  0.378170675522449   met
0 200 classic  0.5  0.015338995845559644 -0.1514729694024941   0.18567154761635063 met
0 200 observed 1.0 -0.048171889576498206 -0.22532683146019306  0.17042075144781954 met
1 1   agent    1.0  0.3797496890170948    0.21654862757777282  0.5358956397643722  met
1 1   classic  0.5  0.09186571058459325   0.01572338218677194  0.25267540652158527 met
1 1   observed 1.0  0.31304455580977164   0.15800816485719293  0.45693864394558137 met
1 200 agent    1.0  0.5520492880559891    0.45557341984461264  0.6209005789942879  met
1 200 classic  0.5  0.19110832387991827   0.10198827668008725  0.2938955779429042  met
1 200 observed 1.0  0.38814320907571975   0.2669877652833774   0.4875214652443116  met
2 1   agent    1.0  0.38360462051338917   0.1705747288311521   0.49986302367437396 met
2 1   classic  0.5 -0.20701028556839307 -0.3692754234355389  -0.07137729565656059 met
2 1   observed 1.0  0.12453670017366081 -0.1596756824545409   0.2683927963428645  met
2 200 agent    1.0  0.2896015993043566    0.11203370615403552  0.50300438118277    met
2 200 classic  0.5  0.018443958845337782 -0.14390641181237784  0.1798794965757633  met
2 200 observed 1.0  0.3225982244043246    0.09014495429795488  0.5013047593399822  met
```

Runs a, b and c (20 pairs each, same fe0492b build, load at start [0.04, 0.41, 0.92], [0.27, 1.42, 1.38], [0.27, 1.21, 1.31]) were also received as full per-combination tables in a later message but are not transcribed here: they carry no verdict weight, because with 20 pairs the 95 percent intervals reach 3.9 percentage points and the null control (classic) was inconclusive in 11 of 18 combinations in every one of the three runs. Verdict counts per run (met / inconclusive over 18 combinations): a 8 / 10, b 10 / 8, c 9 / 9. Their numbers remain on dev08 (`~/rustdft_p6/perf8/`, raw files sha256 given in the transfer messages: run_a 82205e29de6ef7f4..., run_b 94e5ea1bb967c69b..., run_c 164b944eb6e88768...).

## Reading
- Overhead verdict (P11J, provisional because relayed): with 300 pairs all 18 combinations meet their targets (classic below 0.5 percent, observed and agent below 1.0 percent); medians span classic -0.207 to +0.214, observed -0.048 to +0.388, agent +0.012 to +0.552 percent. With 100 pairs 15 of 18 met, with 20 pairs the interval is too wide to judge. The default of `scripts/bench_runtime.py` is now 100 pairs.
- The classic null control is not exactly zero in some combinations (CI excluding zero by up to 0.2 percent) with signs that differ between kinds; no systematic regression is claimed.
- MKL and QE: see parts 1 and 2. Syn-s2 backend differences are explained in `docs/validation/scf-fixed-point-multiplicity.md`.

## Part 4: fresh BAANDI and QE measurement on main fdce016 (after P10M), relayed 2026-11-05

Same session, same idle node, same pw.x and input hashes as part 2. BAANDI 5-run median, QE 5-run mean. Ratio = BAANDI over QE internal time (above 1: BAANDI slower). Builds: pure sha256 `f5803c7a1004c7aa465ab8f6f9849d808d415d57df7cd09f00bf7f7c7d550e19`, mkl sha256 `41da71e9dbbd9df470f4e1458d4a5358739794211511d4262901bec621d7db6e`. Uptime/load: BAANDI before 08:34:12 load 0.29, after 08:38:00 load 1.44; QE before 08:39:36 load 0.29, after 08:40:24 load 2.80. BAANDI spread at most 2.7 percent, QE spread at most 5.9 percent (k444s 1e-8 N1). QE total energies unchanged (k444s -15.85195738 Ry, si8 -62.64362757 Ry). `*`: QE internal time below 0.3 s, not valid for the README (resolution 0.01 s).

```
case   conv     N   QE int(s)  spr%   pure(s)   mkl(s)  pure/QE  mkl/QE  QE iters
k444s  1e-8     1     1.2040    5.9     6.127    5.034     5.09    4.18     5
k444s  1e-8     8     0.5300    0.0     1.173    1.053     2.21    1.99     5
k444s  1e-8    16     0.3840    2.6     0.834    0.788     2.17    2.05     5
k444s  1e-12    1     1.7000    0.0     6.127    5.034     3.60    2.96     8
k444s  1e-12    8     0.7300    0.0     1.173    1.053     1.61    1.44     8
k444s  1e-12   16     0.5320    1.9     0.834    0.788     1.57    1.48     8
si8    1e-8     1     0.6840    2.9     5.801    6.016     8.48    8.80    14
si8    1e-8     8     0.2100*   0.0     4.151    4.809    19.77   22.90    13
si8    1e-8    16     0.2300*   0.0     4.568    4.998    19.86   21.73    13
si8    1e-12    1     0.9120    1.1     5.801    6.016     6.36    6.60    20
si8    1e-12    8     0.2900*   0.0     4.151    4.809    14.31   16.58    21
si8    1e-12   16     0.3000    1.0     4.568    4.998    15.23   16.66    19
```

Worker utilization si8 pure T=8 (12 runs, 250 samples): busy workers 0: 9.6 percent, 1: 29.6, 2: 34.8, 3: 20.4, 4: 4.4, 5: 1.2, 6 or more: 0.0; mean 1.84 of 8 (23.0 percent; P10L 22.5 percent); max/mean imbalance 3.01 (P10L 3.29). Amdahl: serial 9.2 percent (0.382 s), parallel region 3.769 s with worker utilization 26.2 percent; at 100 percent utilization 1.369 s (4.24x); actual 4.151 s (1.40x). si8 pure T=1 5.801 s to T=8 4.151 s is 1.40x.

