# Benchmark suite results

Synthetic results are non-physical. Timings are indicative on a shared host.

| Case | Class | Status | Seconds (median) | Iterations | Energy (Ha) |
|---|---|---|---:|---:|---:|
| si-elongated | elongated | nonconverged | 0.445 | 100 | -15.128518095558501 |
| si-gamma | si-primitive | pass | 0.112 | 10 | -7.295530468210765 |
| si-k222 | si-primitive | pass | 0.656 | 10 | -7.831230057121137 |
| si-k222s | si-primitive | pass | 0.352 | 10 | -7.919423412469083 |
| si-low | si-cutoff | pass | 0.019 | 11 | -7.233756452824798 |
| si-skew | skewed | pass | 0.130 | 13 | -7.250555672960788 |
| si8-gamma | si-conventional | pass | 3.673 | 15 | -31.321803601031256 |
| syn-cubic-degenerate | synthetic-gap | pass | 0.009 | 10 | -4.2612772341353455 |
| syn-d3 | synthetic-projectors | pass | 0.005 | 8 | -0.9465076331236242 |
| syn-fft-primes | awkward-fft | pass | 0.006 | 14 | -0.7346843876893185 |
| syn-gaussian-fcc | synthetic-gap | pass | 0.041 | 17 | -4.68302340536371 |
| syn-odd-shift | synthetic-mesh | pass | 0.010 | 9 | -1.590066401990696 |
| syn-p3 | synthetic-projectors | pass | 0.005 | 8 | -0.9465076597173674 |
| syn-s1 | synthetic-projectors | pass | 0.005 | 13 | -0.7346754280328689 |
| syn-s2 | synthetic-projectors | pass | 0.005 | 17 | 0.6304095283333497 |
| syn-s3 | synthetic-projectors | nonconverged | 0.009 | 100 | 1.7660002310813194 |
| syn-spd | synthetic-projectors | pass | 0.005 | 11 | -0.7346751923893944 |
| syn-z1 | synthetic-valence | pass | 0.008 | 7 | -0.9291990923944965 |
| syn-z14 | synthetic-valence | nonconverged | 0.043 | 100 | -71.13959832030042 |

## By system class

| Class | Cases | Passed | Sum median seconds |
|---|---:|---:|---:|
| awkward-fft | 1 | 1 | 0.006 |
| elongated | 1 | 0 | 0.445 |
| si-conventional | 1 | 1 | 3.673 |
| si-cutoff | 1 | 1 | 0.019 |
| si-primitive | 3 | 3 | 1.119 |
| skewed | 1 | 1 | 0.130 |
| synthetic-gap | 2 | 2 | 0.050 |
| synthetic-mesh | 1 | 1 | 0.010 |
| synthetic-projectors | 6 | 5 | 0.035 |
| synthetic-valence | 2 | 1 | 0.051 |

## Baseline comparison

Only selected cases are compared; missing selected entries fail comparison.
- si-elongated: failure message changed
- si-gamma: unchanged
- si-k222: unchanged
- si-k222s: unchanged
- si-low: unchanged
- si-skew: unchanged
- si8-gamma: previous failure now passes
- syn-cubic-degenerate: incomparable input or UPF (includes seed and threads)
- syn-d3: incomparable input or UPF (includes seed and threads)
- syn-fft-primes: incomparable input or UPF (includes seed and threads)
- syn-gaussian-fcc: incomparable input or UPF (includes seed and threads)
- syn-odd-shift: incomparable input or UPF (includes seed and threads)
- syn-p3: incomparable input or UPF (includes seed and threads)
- syn-s1: incomparable input or UPF (includes seed and threads)
- syn-s2: incomparable input or UPF (includes seed and threads)
- syn-s3: incomparable input or UPF (includes seed and threads)
- syn-spd: incomparable input or UPF (includes seed and threads)
- syn-z1: incomparable input or UPF (includes seed and threads)
- syn-z14: incomparable input or UPF (includes seed and threads)

## Findings and exact reproducing inputs

Replace SI_UPF_PATH with the read-only validation file. For synthetic cases regenerate synthetic.UPF using the case parameters in JSON. Runtime directories contain executable inputs and logs.

### si-elongated

Expected failure matched: False

```text
P14B_GEMM,reduced_hamiltonian,12,12,569,100,16387200,1716796
P14B_GEMM,ritz_and_residual,569,12,12,200,32774400,2859711
P14B_GEMM,projection_overlap_1,12,12,569,49,8029728,830084
P14B_GEMM,projection_subtract_1,569,12,12,49,8029728,693045
P14B_GEMM,projection_overlap_2,12,12,569,49,8029728,808878
P14B_GEMM,projection_subtract_2,569,12,12,49,8029728,688041
P14B_GEMM,reduced_hamiltonian,24,24,569,55,36051840,3205418
P14B_GEMM,ritz_and_residual,569,12,24,110,36051840,2868917
P14B_GEMM,projection_overlap_1,24,12,569,18,5899392,558624
P14B_GEMM,projection_subtract_1,569,12,24,18,5899392,468498
P14B_GEMM,projection_overlap_2,24,12,569,18,5899392,555541
P14B_GEMM,projection_subtract_2,569,12,24,18,5899392,469790
P14B_GEMM,reduced_hamiltonian,36,36,569,25,36871200,3110791
P14B_GEMM,ritz_and_residual,569,12,36,50,24580800,1927466
P14B_GEMM,projection_overlap_1,36,12,569,2,983232,92709
P14B_GEMM,projection_subtract_1,569,12,36,2,983232,77167
P14B_GEMM,projection_overlap_2,36,12,569,2,983232,90041
P14B_GEMM,projection_subtract_2,569,12,36,2,983232,77000
P14B_GEMM,reduced_hamiltonian,48,48,569,7,18353664,1507625
P14B_GEMM,ritz_and_residual,569,12,48,14,9176832,715709
P14B_GEMM,projection_overlap_1,48,12,569,1,655488,61417
P14B_GEMM,projection_subtract_1,569,12,48,1,655488,51166
P14B_GEMM,projection_overlap_2,48,12,569,1,655488,59334
P14B_GEMM,projection_subtract_2,569,12,48,1,655488,51084
P14B_GEMM,reduced_hamiltonian,60,60,569,1,4096800,339125
P14B_GEMM,ritz_and_residual,569,12,60,2,1638720,127334
P14B_GEMM,projection_overlap_1,60,11,569,1,751080,73833
P14B_GEMM,projection_subtract_1,569,11,60,1,751080,62458
P14B_GEMM,projection_overlap_2,60,11,569,1,751080,73500
P14B_GEMM,projection_subtract_2,569,11,60,1,751080,62209
P14B_GEMM,reduced_hamiltonian,71,71,569,1,5736658,476000
P14B_GEMM,ritz_and_residual,569,12,71,2,1939152,150084
P14B_GEMM,projection_overlap_1,12,7,569,2,191184,23500
P14B_GEMM,projection_subtract_1,569,7,12,2,191184,17500
P14B_GEMM,projection_overlap_2,12,7,569,2,191184,23958
P14B_GEMM,projection_subtract_2,569,7,12,2,191184,17125
P14B_GEMM,reduced_hamiltonian,19,19,569,5,2054090,211791
P14B_GEMM,ritz_and_residual,569,12,19,10,2594640,211124
P14B_GEMM,projection_overlap_1,19,4,569,1,86488,11292
P14B_GEMM,projection_subtract_1,569,4,19,1,86488,7500
P14B_GEMM,projection_overlap_2,19,4,569,1,86488,11125
P14B_GEMM,projection_subtract_2,569,4,19,1,86488,7291
P14B_GEMM,reduced_hamiltonian,23,23,569,9,5418018,516999
P14B_GEMM,ritz_and_residual,569,12,23,18,5653584,450459
P14B_GEMM,projection_overlap_1,23,3,569,1,78522,13125
P14B_GEMM,projection_subtract_1,569,3,23,1,78522,7459
P14B_GEMM,projection_overlap_2,23,3,569,1,78522,12834
P14B_GEMM,projection_subtract_2,569,3,23,1,78522,7334
P14B_GEMM,reduced_hamiltonian,26,26,569,12,9231456,879041
P14B_GEMM,ritz_and_residual,569,12,26,24,8521344,671916
P14B_GEMM,projection_overlap_1,26,3,569,1,88764,14125
P14B_GEMM,projection_subtract_1,569,3,26,1,88764,9125
P14B_GEMM,projection_overlap_2,26,3,569,1,88764,14125
P14B_GEMM,projection_subtract_2,569,3,26,1,88764,8125
P14B_GEMM,reduced_hamiltonian,29,29,569,6,5742348,545291
P14B_GEMM,ritz_and_residual,569,12,29,12,4752288,374458
P14B_GEMM,projection_overlap_1,29,3,569,1,99006,15625
P14B_GEMM,projection_subtract_1,569,3,29,1,99006,9125
P14B_GEMM,projection_overlap_2,29,3,569,1,99006,15583
P14B_GEMM,projection_subtract_2,569,3,29,1,99006,8958
P14B_GEMM,reduced_hamiltonian,32,32,569,12,13983744,1221999
P14B_GEMM,ritz_and_residual,569,12,32,24,10487808,825502
P14B_GEMM,projection_overlap_1,32,2,569,3,218496,49916
P14B_GEMM,projection_subtract_1,569,2,32,3,218496,25292
P14B_GEMM,projection_overlap_2,32,2,569,3,218496,49417
P14B_GEMM,projection_subtract_2,569,2,32,3,218496,24709
P14B_GEMM,reduced_hamiltonian,34,34,569,9,11839752,1097583
P14B_GEMM,ritz_and_residual,569,12,34,18,8357472,655920
P14B_GEMM,projection_overlap_1,34,1,569,2,77384,15084
P14B_GEMM,projection_subtract_1,569,1,34,2,77384,10750
P14B_GEMM,projection_overlap_2,34,1,569,2,77384,14875
P14B_GEMM,projection_subtract_2,569,1,34,2,77384,10417
P14B_GEMM,reduced_hamiltonian,35,35,569,10,13940500,1233001
P14B_GEMM,ritz_and_residual,569,12,35,20,9559200,748791
P14B_GEMM,projection_overlap_1,35,1,569,6,238980,45707
P14B_GEMM,projection_subtract_1,569,1,35,6,238980,32832
P14B_GEMM,projection_overlap_2,35,1,569,6,238980,45543
P14B_GEMM,projection_subtract_2,569,1,35,6,238980,32084
P14B_GEMM,projection_overlap_1,48,8,569,1,436992,42500
P14B_GEMM,projection_subtract_1,569,8,48,1,436992,34167
P14B_GEMM,projection_overlap_2,48,8,569,1,436992,42250
P14B_GEMM,projection_subtract_2,569,8,48,1,436992,34292
P14B_GEMM,reduced_hamiltonian,56,56,569,1,3568768,293708
P14B_GEMM,ritz_and_residual,569,12,56,2,1529472,119416
P14B_GEMM,projection_overlap_1,56,2,569,1,127456,27875
P14B_GEMM,projection_subtract_1,569,2,56,1,127456,15625
P14B_GEMM,projection_overlap_2,56,2,569,1,127456,27416
P14B_GEMM,projection_subtract_2,569,2,56,1,127456,15084
P14B_GEMM,reduced_hamiltonian,58,58,569,1,3828232,327375
P14B_GEMM,ritz_and_residual,569,12,58,2,1584096,123125
P14B_GEMM,projection_overlap_1,36,6,569,3,737424,92917
P14B_GEMM,projection_subtract_1,569,6,36,3,737424,65625
P14B_GEMM,projection_overlap_2,36,6,569,3,737424,92584
P14B_GEMM,projection_subtract_2,569,6,36,3,737424,65291
P14B_GEMM,reduced_hamiltonian,42,42,569,7,14052024,1210626
P14B_GEMM,ritz_and_residual,569,12,42,14,8029728,627458
P14B_GEMM,projection_overlap_1,42,2,569,2,191184,43167
P14B_GEMM,projection_subtract_1,569,2,42,2,191184,21542
P14B_GEMM,projection_overlap_2,42,2,569,2,191184,41333
P14B_GEMM,projection_subtract_2,569,2,42,2,191184,21250
P14B_GEMM,reduced_hamiltonian,44,44,569,5,11015840,927750
P14B_GEMM,ritz_and_residual,569,12,44,10,6008640,469417
P14B_GEMM,projection_overlap_1,36,10,569,3,1229040,132374
P14B_GEMM,projection_subtract_1,569,10,36,3,1229040,103875
P14B_GEMM,projection_overlap_2,36,10,569,3,1229040,132291
P14B_GEMM,projection_subtract_2,569,10,36,3,1229040,103792
P14B_GEMM,reduced_hamiltonian,46,46,569,7,16856056,1487793
P14B_GEMM,ritz_and_residual,569,12,46,14,8794464,686210
P14B_GEMM,projection_overlap_1,46,1,569,3,157044,29918
P14B_GEMM,projection_subtract_1,569,1,46,3,157044,21708
P14B_GEMM,projection_overlap_2,46,1,569,3,157044,29876
P14B_GEMM,projection_subtract_2,569,1,46,3,157044,21584
P14B_GEMM,reduced_hamiltonian,47,47,569,9,22624578,1924875
P14B_GEMM,ritz_and_residual,569,12,47,18,11552976,900378
P14B_GEMM,projection_overlap_1,24,2,569,5,273120,62376
P14B_GEMM,projection_subtract_1,569,2,24,5,273120,29167
P14B_GEMM,projection_overlap_2,24,2,569,5,273120,61792
P14B_GEMM,projection_subtract_2,569,2,24,5,273120,28833
P14B_GEMM,projection_overlap_1,24,11,569,3,901296,92583
P14B_GEMM,projection_subtract_1,569,11,24,3,901296,74624
P14B_GEMM,projection_overlap_2,24,11,569,3,901296,92167
P14B_GEMM,projection_subtract_2,569,11,24,3,901296,74167
P14B_GEMM,projection_overlap_1,35,4,569,1,159320,19084
P14B_GEMM,projection_subtract_1,569,4,35,1,159320,12875
P14B_GEMM,projection_overlap_2,35,4,569,1,159320,19042
P14B_GEMM,projection_subtract_2,569,4,35,1,159320,12875
P14B_GEMM,reduced_hamiltonian,39,39,569,6,10385388,887208
P14B_GEMM,ritz_and_residual,569,12,39,12,6391008,499876
P14B_GEMM,projection_overlap_1,39,2,569,4,355056,77333
P14B_GEMM,projection_subtract_1,569,2,39,4,355056,39666
P14B_GEMM,projection_overlap_2,39,2,569,4,355056,77000
P14B_GEMM,projection_subtract_2,569,2,39,4,355056,39583
P14B_GEMM,reduced_hamiltonian,41,41,569,7,13390846,1195833
P14B_GEMM,ritz_and_residual,569,12,41,14,7838544,613749
P14B_GEMM,projection_overlap_1,41,2,569,4,373264,82833
P14B_GEMM,projection_subtract_1,569,2,41,4,373264,41916
P14B_GEMM,projection_overlap_2,41,2,569,4,373264,82249
P14B_GEMM,projection_subtract_2,569,2,41,4,373264,41668
P14B_GEMM,reduced_hamiltonian,43,43,569,8,16833296,1476333
P14B_GEMM,ritz_and_residual,569,12,43,16,9395328,733706
P14B_GEMM,projection_overlap_1,43,1,569,2,97868,18666
P14B_GEMM,projection_subtract_1,569,1,43,2,97868,13501
P14B_GEMM,projection_overlap_2,43,1,569,2,97868,18542
P14B_GEMM,projection_subtract_2,569,1,43,2,97868,13375
P14B_GEMM,projection_overlap_1,44,1,569,2,100144,19041
P14B_GEMM,projection_subtract_1,569,1,44,2,100144,13834
P14B_GEMM,projection_overlap_2,44,1,569,2,100144,19042
P14B_GEMM,projection_subtract_2,569,1,44,2,100144,13749
P14B_GEMM,reduced_hamiltonian,45,45,569,8,18435600,1596582
P14B_GEMM,ritz_and_residual,569,12,45,16,9832320,767335
P14B_GEMM,projection_overlap_1,45,1,569,2,102420,19500
P14B_GEMM,projection_subtract_1,569,1,45,2,102420,14208
P14B_GEMM,projection_overlap_2,45,1,569,2,102420,19584
P14B_GEMM,projection_subtract_2,569,1,45,2,102420,13791
P14B_GEMM,projection_overlap_1,47,1,569,3,160458,30625
P14B_GEMM,projection_subtract_1,569,1,47,3,160458,22083
P14B_GEMM,projection_overlap_2,47,1,569,3,160458,30624
P14B_GEMM,projection_subtract_2,569,1,47,3,160458,21584
P14B_GEMM,projection_overlap_1,48,1,569,3,163872,31250
P14B_GEMM,projection_subtract_1,569,1,48,3,163872,22416
P14B_GEMM,projection_overlap_2,48,1,569,3,163872,31208
P14B_GEMM,projection_subtract_2,569,1,48,3,163872,22167
P14B_GEMM,reduced_hamiltonian,49,49,569,6,16394028,1449458
P14B_GEMM,ritz_and_residual,569,12,49,12,8029728,627626
P14B_GEMM,projection_overlap_1,24,1,569,11,300432,57374
P14B_GEMM,projection_subtract_1,569,1,24,11,300432,41539
P14B_GEMM,projection_overlap_2,24,1,569,11,300432,57163
P14B_GEMM,projection_subtract_2,569,1,24,11,300432,40125
P14B_GEMM,reduced_hamiltonian,25,25,569,11,7823750,790709
P14B_GEMM,ritz_and_residual,569,12,25,22,7510800,594333
P14B_GEMM,projection_overlap_1,36,5,569,2,409680,59499
P14B_GEMM,projection_subtract_1,569,5,36,2,409680,38292
P14B_GEMM,projection_overlap_2,36,5,569,2,409680,59291
P14B_GEMM,projection_subtract_2,569,5,36,2,409680,38083
P14B_GEMM,projection_overlap_1,24,6,569,3,491616,63334
P14B_GEMM,projection_subtract_1,569,6,24,3,491616,43209
P14B_GEMM,projection_overlap_2,24,6,569,3,491616,63125
P14B_GEMM,projection_subtract_2,569,6,24,3,491616,42833
P14B_GEMM,reduced_hamiltonian,30,30,569,8,8193600,740126
P14B_GEMM,ritz_and_residual,569,12,30,16,6554880,516168
P14B_GEMM,projection_overlap_1,30,3,569,1,102420,16333
P14B_GEMM,projection_subtract_1,569,3,30,1,102420,9833
P14B_GEMM,projection_overlap_2,30,3,569,1,102420,16000
P14B_GEMM,projection_subtract_2,569,3,30,1,102420,9458
P14B_GEMM,reduced_hamiltonian,33,33,569,8,9914256,896835
P14B_GEMM,ritz_and_residual,569,12,33,16,7210368,566999
P14B_GEMM,projection_overlap_1,33,2,569,3,225324,50834
P14B_GEMM,projection_subtract_1,569,2,33,3,225324,25917
P14B_GEMM,projection_overlap_2,33,2,569,3,225324,50041
P14B_GEMM,projection_subtract_2,569,2,33,3,225324,25708
P14B_GEMM,projection_overlap_1,36,1,569,3,122904,23459
P14B_GEMM,projection_subtract_1,569,1,36,3,122904,16583
P14B_GEMM,projection_overlap_2,36,1,569,3,122904,23375
P14B_GEMM,projection_subtract_2,569,1,36,3,122904,16375
P14B_GEMM,reduced_hamiltonian,37,37,569,8,12463376,1154166
P14B_GEMM,ritz_and_residual,569,12,37,16,8084352,632665
P14B_GEMM,projection_overlap_1,12,9,569,1,122904,15917
P14B_GEMM,projection_subtract_1,569,9,12,1,122904,11500
P14B_GEMM,projection_overlap_2,12,9,569,1,122904,15542
P14B_GEMM,projection_subtract_2,569,9,12,1,122904,11292
P14B_GEMM,reduced_hamiltonian,21,21,569,2,1003716,99583
P14B_GEMM,ritz_and_residual,569,12,21,4,1147104,92377
P14B_GEMM,projection_overlap_1,12,3,569,3,122904,24000
P14B_GEMM,projection_subtract_1,569,3,12,3,122904,13000
P14B_GEMM,projection_overlap_2,12,3,569,3,122904,20541
P14B_GEMM,projection_subtract_2,569,3,12,3,122904,11958
P14B_GEMM,reduced_hamiltonian,15,15,569,3,768150,79167
P14B_GEMM,ritz_and_residual,569,12,15,6,1229040,103623
P14B_GEMM,projection_overlap_1,25,1,569,6,170700,32625
P14B_GEMM,projection_subtract_1,569,1,25,6,170700,23334
P14B_GEMM,projection_overlap_2,25,1,569,6,170700,32499
P14B_GEMM,projection_subtract_2,569,1,25,6,170700,22833
P14B_GEMM,projection_overlap_1,26,1,569,7,207116,39666
P14B_GEMM,projection_subtract_1,569,1,26,7,207116,28167
P14B_GEMM,projection_overlap_2,26,1,569,7,207116,39457
P14B_GEMM,projection_subtract_2,569,1,26,7,207116,27624
P14B_GEMM,reduced_hamiltonian,27,27,569,8,6636816,595082
P14B_GEMM,ritz_and_residual,569,12,27,16,5899392,466249
P14B_GEMM,projection_overlap_1,27,1,569,4,122904,23583
P14B_GEMM,projection_subtract_1,569,1,27,4,122904,16959
P14B_GEMM,projection_overlap_2,27,1,569,4,122904,23416
P14B_GEMM,projection_subtract_2,569,1,27,4,122904,16415
P14B_GEMM,reduced_hamiltonian,28,28,569,9,8029728,734543
P14B_GEMM,ritz_and_residual,569,12,28,18,6882624,542582
P14B_GEMM,projection_overlap_1,12,11,569,4,600864,67416
P14B_GEMM,projection_subtract_1,569,11,12,4,600864,53041
P14B_GEMM,projection_overlap_2,12,11,569,4,600864,65749
P14B_GEMM,projection_subtract_2,569,11,12,4,600864,52374
P14B_GEMM,projection_overlap_1,23,1,569,6,157044,30042
P14B_GEMM,projection_subtract_1,569,1,23,6,157044,21626
P14B_GEMM,projection_overlap_2,23,1,569,6,157044,29916
P14B_GEMM,projection_subtract_2,569,1,23,6,157044,21042
P14B_GEMM,projection_overlap_1,12,6,569,1,81936,11500
P14B_GEMM,projection_subtract_1,569,6,12,1,81936,8083
P14B_GEMM,projection_overlap_2,12,6,569,1,81936,11292
P14B_GEMM,projection_subtract_2,569,6,12,1,81936,7584
P14B_GEMM,reduced_hamiltonian,18,18,569,3,1106136,110958
P14B_GEMM,ritz_and_residual,569,12,18,6,1474848,122042
P14B_GEMM,projection_overlap_1,18,2,569,1,40968,9500
P14B_GEMM,projection_subtract_1,569,2,18,1,40968,4375
P14B_GEMM,projection_overlap_2,18,2,569,1,40968,9375
P14B_GEMM,projection_subtract_2,569,2,18,1,40968,4292
P14B_GEMM,reduced_hamiltonian,20,20,569,4,1820800,172375
P14B_GEMM,ritz_and_residual,569,12,20,8,2184960,176957
P14B_GEMM,projection_overlap_1,20,2,569,1,45520,10792
P14B_GEMM,projection_subtract_1,569,2,20,1,45520,4750
P14B_GEMM,projection_overlap_2,20,2,569,1,45520,10625
P14B_GEMM,projection_subtract_2,569,2,20,1,45520,4708
P14B_GEMM,reduced_hamiltonian,22,22,569,6,3304752,337749
P14B_GEMM,ritz_and_residual,569,12,22,12,3605184,288750
P14B_GEMM,projection_overlap_1,22,1,569,4,100144,19123
P14B_GEMM,projection_subtract_1,569,1,22,4,100144,13916
P14B_GEMM,projection_overlap_2,22,1,569,4,100144,19000
P14B_GEMM,projection_subtract_2,569,1,22,4,100144,13375
P14B_GEMM,projection_overlap_1,12,4,569,2,109248,14292
P14B_GEMM,projection_subtract_1,569,4,12,2,109248,10374
P14B_GEMM,projection_overlap_2,12,4,569,2,109248,14083
P14B_GEMM,projection_subtract_2,569,4,12,2,109248,9959
P14B_GEMM,reduced_hamiltonian,16,16,569,2,582656,61000
P14B_GEMM,ritz_and_residual,569,12,16,4,873984,73833
P14B_GEMM,projection_overlap_1,15,2,569,1,34140,8000
P14B_GEMM,projection_subtract_1,569,2,15,1,34140,4209
P14B_GEMM,projection_overlap_2,15,2,569,1,34140,7833
P14B_GEMM,projection_subtract_2,569,2,15,1,34140,3750
P14B_GEMM,reduced_hamiltonian,17,17,569,5,1644410,180373
P14B_GEMM,ritz_and_residual,569,12,17,10,2321520,194875
P14B_GEMM,projection_overlap_1,17,2,569,1,38692,9292
P14B_GEMM,projection_subtract_1,569,2,17,1,38692,4208
P14B_GEMM,projection_overlap_2,17,2,569,1,38692,9208
P14B_GEMM,projection_subtract_2,569,2,17,1,38692,4084
P14B_GEMM,projection_overlap_1,19,2,569,1,43244,10667
P14B_GEMM,projection_subtract_1,569,2,19,1,43244,4542
P14B_GEMM,projection_overlap_2,19,2,569,1,43244,10542
P14B_GEMM,projection_subtract_2,569,2,19,1,43244,4500
P14B_GEMM,projection_overlap_1,21,1,569,1,23898,4584
P14B_GEMM,projection_subtract_1,569,1,21,1,23898,3333
P14B_GEMM,projection_overlap_2,21,1,569,1,23898,4541
P14B_GEMM,projection_subtract_2,569,1,21,1,23898,3209
P14B_GEMM,projection_overlap_1,12,5,569,3,204840,32833
P14B_GEMM,projection_subtract_1,569,5,12,3,204840,20875
P14B_GEMM,projection_overlap_2,12,5,569,3,204840,32500
P14B_GEMM,projection_subtract_2,569,5,12,3,204840,19833
P14B_GEMM,projection_overlap_1,17,1,569,2,38692,7416
P14B_GEMM,projection_subtract_1,569,1,17,2,38692,5375
P14B_GEMM,projection_overlap_2,17,1,569,2,38692,7333
P14B_GEMM,projection_subtract_2,569,1,17,2,38692,5333
P14B_GEMM,projection_overlap_1,18,1,569,2,40968,7875
P14B_GEMM,projection_subtract_1,569,1,18,2,40968,5834
P14B_GEMM,projection_overlap_2,18,1,569,2,40968,7833
P14B_GEMM,projection_subtract_2,569,1,18,2,40968,5626
P14B_GEMM,projection_overlap_1,12,10,569,4,546240,65166
P14B_GEMM,projection_subtract_1,569,10,12,4,546240,48959
P14B_GEMM,projection_overlap_2,12,10,569,4,546240,63667
P14B_GEMM,projection_subtract_2,569,10,12,4,546240,48250
P14B_GEMM,projection_overlap_1,22,7,569,1,175252,21458
P14B_GEMM,projection_subtract_1,569,7,22,1,175252,15042
P14B_GEMM,projection_overlap_2,22,7,569,1,175252,21292
P14B_GEMM,projection_subtract_2,569,7,22,1,175252,14917
P14B_GEMM,projection_overlap_1,29,4,569,1,132008,15959
P14B_GEMM,projection_subtract_1,569,4,29,1,132008,11041
P14B_GEMM,projection_overlap_2,29,4,569,1,132008,16042
P14B_GEMM,projection_subtract_2,569,4,29,1,132008,10750
P14B_GEMM,projection_overlap_1,37,1,569,3,126318,24168
P14B_GEMM,projection_subtract_1,569,1,37,3,126318,17207
P14B_GEMM,projection_overlap_2,37,1,569,3,126318,24042
P14B_GEMM,projection_subtract_2,569,1,37,3,126318,16875
P14B_GEMM,reduced_hamiltonian,38,38,569,3,4929816,438375
P14B_GEMM,ritz_and_residual,569,12,38,6,3113568,243959
P14B_GEMM,projection_overlap_1,26,2,569,3,177528,41291
P14B_GEMM,projection_subtract_1,569,2,26,3,177528,19041
P14B_GEMM,projection_overlap_2,26,2,569,3,177528,40832
P14B_GEMM,projection_subtract_2,569,2,26,3,177528,18876
P14B_GEMM,projection_overlap_1,28,2,569,3,191184,45083
P14B_GEMM,projection_subtract_1,569,2,28,3,191184,20917
P14B_GEMM,projection_overlap_2,28,2,569,3,191184,44708
P14B_GEMM,projection_subtract_2,569,2,28,3,191184,20792
P14B_GEMM,projection_overlap_1,30,2,569,4,273120,61542
P14B_GEMM,projection_subtract_1,569,2,30,4,273120,31332
P14B_GEMM,projection_overlap_2,30,2,569,4,273120,60708
P14B_GEMM,projection_subtract_2,569,2,30,4,273120,31084
P14B_GEMM,projection_overlap_1,32,1,569,5,182080,34708
P14B_GEMM,projection_subtract_1,569,1,32,5,182080,25083
P14B_GEMM,projection_overlap_2,32,1,569,5,182080,34625
P14B_GEMM,projection_subtract_2,569,1,32,5,182080,24374
P14B_GEMM,projection_overlap_1,46,4,569,2,418784,50083
P14B_GEMM,projection_subtract_1,569,4,46,2,418784,33416
P14B_GEMM,projection_overlap_2,46,4,569,2,418784,49958
P14B_GEMM,projection_subtract_2,569,4,46,2,418784,33375
P14B_GEMM,reduced_hamiltonian,50,50,569,4,11380000,977083
P14B_GEMM,ritz_and_residual,569,12,50,8,5462400,425666
P14B_GEMM,projection_overlap_1,50,2,569,1,113800,24958
P14B_GEMM,projection_subtract_1,569,2,50,1,113800,13208
P14B_GEMM,projection_overlap_2,50,2,569,1,113800,24709
P14B_GEMM,projection_subtract_2,569,2,50,1,113800,12959
P14B_GEMM,reduced_hamiltonian,52,52,569,1,3077152,259000
P14B_GEMM,ritz_and_residual,569,12,52,2,1420224,110708
P14B_GEMM,projection_overlap_1,52,2,569,1,118352,26208
P14B_GEMM,projection_subtract_1,569,2,52,1,118352,13875
P14B_GEMM,projection_overlap_2,52,2,569,1,118352,25959
P14B_GEMM,projection_subtract_2,569,2,52,1,118352,13791
P14B_GEMM,reduced_hamiltonian,54,54,569,2,6636816,556166
P14B_GEMM,ritz_and_residual,569,12,54,4,2949696,229917
P14B_GEMM,projection_overlap_1,36,3,569,3,368712,57041
P14B_GEMM,projection_subtract_1,569,3,36,3,368712,34749
P14B_GEMM,projection_overlap_2,36,3,569,3,368712,56750
P14B_GEMM,projection_subtract_2,569,3,36,3,368712,34542
P14B_GEMM,projection_overlap_1,41,1,569,2,93316,17750
P14B_GEMM,projection_subtract_1,569,1,41,2,93316,12708
P14B_GEMM,projection_overlap_2,41,1,569,2,93316,17750
P14B_GEMM,projection_subtract_2,569,1,41,2,93316,12418
P14B_GEMM,projection_overlap_1,35,2,569,3,238980,54624
P14B_GEMM,projection_subtract_1,569,2,35,3,238980,27583
P14B_GEMM,projection_overlap_2,35,2,569,3,238980,54083
P14B_GEMM,projection_subtract_2,569,2,35,3,238980,27291
P14B_GEMM,projection_overlap_1,37,2,569,2,168424,38083
P14B_GEMM,projection_subtract_1,569,2,37,2,168424,19083
P14B_GEMM,projection_overlap_2,37,2,569,2,168424,37917
P14B_GEMM,projection_subtract_2,569,2,37,2,168424,18917
P14B_GEMM,projection_overlap_1,39,3,569,1,133146,20208
P14B_GEMM,projection_subtract_1,569,3,39,1,133146,12541
P14B_GEMM,projection_overlap_2,39,3,569,1,133146,20250
P14B_GEMM,projection_subtract_2,569,3,39,1,133146,12333
P14B_GEMM,projection_overlap_1,42,1,569,1,47796,9084
P14B_GEMM,projection_subtract_1,569,1,42,1,47796,6542
P14B_GEMM,projection_overlap_2,42,1,569,1,47796,9125
P14B_GEMM,projection_subtract_2,569,1,42,1,47796,6417
P14B_GEMM,projection_overlap_1,24,3,569,1,81936,12916
P14B_GEMM,projection_subtract_1,569,3,24,1,81936,7917
P14B_GEMM,projection_overlap_2,24,3,569,1,81936,12750
P14B_GEMM,projection_subtract_2,569,3,24,1,81936,7708
P14B_GEMM,projection_overlap_1,27,2,569,1,61452,14042
P14B_GEMM,projection_subtract_1,569,2,27,1,61452,6625
P14B_GEMM,projection_overlap_2,27,2,569,1,61452,13833
P14B_GEMM,projection_subtract_2,569,2,27,1,61452,6667
P14B_GEMM,projection_overlap_1,29,2,569,1,66004,15167
P14B_GEMM,projection_subtract_1,569,2,29,1,66004,7334
P14B_GEMM,projection_overlap_2,29,2,569,1,66004,15084
P14B_GEMM,projection_subtract_2,569,2,29,1,66004,7375
P14B_GEMM,reduced_hamiltonian,31,31,569,7,7655326,705958
P14B_GEMM,ritz_and_residual,569,12,31,14,5926704,466166
P14B_GEMM,projection_overlap_1,31,1,569,4,141112,27042
P14B_GEMM,projection_subtract_1,569,1,31,4,141112,19291
P14B_GEMM,projection_overlap_2,31,1,569,4,141112,26959
P14B_GEMM,projection_subtract_2,569,1,31,4,141112,18876
P14B_GEMM,projection_overlap_1,43,2,569,4,391472,87875
P14B_GEMM,projection_subtract_1,569,2,43,4,391472,43708
P14B_GEMM,projection_overlap_2,43,2,569,4,391472,87084
P14B_GEMM,projection_subtract_2,569,2,43,4,391472,43542
P14B_GEMM,projection_overlap_1,45,2,569,5,512100,110875
P14B_GEMM,projection_subtract_1,569,2,45,5,512100,56752
P14B_GEMM,projection_overlap_2,45,2,569,5,512100,110666
P14B_GEMM,projection_subtract_2,569,2,45,5,512100,56584
P14B_GEMM,projection_overlap_1,30,1,569,2,68280,13042
P14B_GEMM,projection_subtract_1,569,1,30,2,68280,9376
P14B_GEMM,projection_overlap_2,30,1,569,2,68280,13042
P14B_GEMM,projection_subtract_2,569,1,30,2,68280,9082
P14B_GEMM,projection_overlap_1,33,1,569,3,112662,21457
P14B_GEMM,projection_subtract_1,569,1,33,3,112662,15291
P14B_GEMM,projection_overlap_2,33,1,569,3,112662,21375
P14B_GEMM,projection_subtract_2,569,1,33,3,112662,15041
P14B_GEMM,projection_overlap_1,24,8,569,3,655488,66124
P14B_GEMM,projection_subtract_1,569,8,24,3,655488,52833
P14B_GEMM,projection_overlap_2,24,8,569,3,655488,65793
P14B_GEMM,projection_subtract_2,569,8,24,3,655488,52542
P14B_GEMM,projection_overlap_1,34,2,569,1,77384,17958
P14B_GEMM,projection_subtract_1,569,2,34,1,77384,8875
P14B_GEMM,projection_overlap_2,34,2,569,1,77384,17750
P14B_GEMM,projection_subtract_2,569,2,34,1,77384,8833
P14B_GEMM,projection_overlap_1,12,1,569,1,13656,2625
P14B_GEMM,projection_subtract_1,569,1,12,1,13656,2083
P14B_GEMM,projection_overlap_2,12,1,569,1,13656,2584
P14B_GEMM,projection_subtract_2,569,1,12,1,13656,1958
P14B_GEMM,reduced_hamiltonian,13,13,569,1,192322,24583
P14B_GEMM,ritz_and_residual,569,12,13,2,355056,30625
P14B_GEMM,projection_overlap_1,36,7,569,3,860328,95833
P14B_GEMM,projection_subtract_1,569,7,36,3,860328,72374
P14B_GEMM,projection_overlap_2,36,7,569,3,860328,95708
P14B_GEMM,projection_subtract_2,569,7,36,3,860328,72041
P14B_GEMM,projection_overlap_1,47,2,569,3,320916,70417
P14B_GEMM,projection_subtract_1,569,2,47,3,320916,36125
P14B_GEMM,projection_overlap_2,47,2,569,3,320916,69957
P14B_GEMM,projection_subtract_2,569,2,47,3,320916,36000
P14B_GEMM,projection_overlap_1,49,1,569,2,111524,21124
P14B_GEMM,projection_subtract_1,569,1,49,2,111524,14959
P14B_GEMM,projection_overlap_2,49,1,569,2,111524,21209
P14B_GEMM,projection_subtract_2,569,1,49,2,111524,14667
P14B_GEMM,projection_overlap_1,50,1,569,1,56900,10875
P14B_GEMM,projection_subtract_1,569,1,50,1,56900,7583
P14B_GEMM,projection_overlap_2,50,1,569,1,56900,10833
P14B_GEMM,projection_subtract_2,569,1,50,1,56900,7459
P14B_GEMM,reduced_hamiltonian,51,51,569,1,2959938,245625
P14B_GEMM,ritz_and_residual,569,12,51,2,1392912,108542
P14B_GEMM,projection_overlap_1,28,1,569,2,63728,12209
P14B_GEMM,projection_subtract_1,569,1,28,2,63728,8751
P14B_GEMM,projection_overlap_2,28,1,569,2,63728,12167
P14B_GEMM,projection_subtract_2,569,1,28,2,63728,8500
P14B_GEMM,projection_overlap_1,32,3,569,2,218496,35041
P14B_GEMM,projection_subtract_1,569,3,32,2,218496,21125
P14B_GEMM,projection_overlap_2,32,3,569,2,218496,34543
P14B_GEMM,projection_subtract_2,569,3,32,2,218496,20333
P14B_GEMM,projection_overlap_1,19,1,569,2,43244,8333
P14B_GEMM,projection_subtract_1,569,1,19,2,43244,6125
P14B_GEMM,projection_overlap_2,19,1,569,2,43244,8250
P14B_GEMM,projection_subtract_2,569,1,19,2,43244,5917
P14B_GEMM,projection_overlap_1,50,3,569,1,170700,26333
P14B_GEMM,projection_subtract_1,569,3,50,1,170700,16500
P14B_GEMM,projection_overlap_2,50,3,569,1,170700,26083
P14B_GEMM,projection_subtract_2,569,3,50,1,170700,16208
P14B_GEMM,reduced_hamiltonian,53,53,569,1,3196642,275333
P14B_GEMM,ritz_and_residual,569,12,53,2,1447536,112624
P14B_GEMM,projection_overlap_1,53,1,569,1,60314,11542
P14B_GEMM,projection_subtract_1,569,1,53,1,60314,8125
P14B_GEMM,projection_overlap_2,53,1,569,1,60314,11458
P14B_GEMM,projection_subtract_2,569,1,53,1,60314,7958
P14B_GEMM,projection_overlap_1,16,1,569,1,18208,3500
P14B_GEMM,projection_subtract_1,569,1,16,1,18208,2584
P14B_GEMM,projection_overlap_2,16,1,569,1,18208,3500
P14B_GEMM,projection_subtract_2,569,1,16,1,18208,2542
P14B_GEMM,projection_overlap_1,24,7,569,3,573552,65500
P14B_GEMM,projection_subtract_1,569,7,24,3,573552,48624
P14B_GEMM,projection_overlap_2,24,7,569,3,573552,65333
P14B_GEMM,projection_subtract_2,569,7,24,3,573552,48292
P14B_GEMM,projection_overlap_1,31,3,569,2,211668,34084
P14B_GEMM,projection_subtract_1,569,3,31,2,211668,19250
P14B_GEMM,projection_overlap_2,31,3,569,2,211668,33875
P14B_GEMM,projection_subtract_2,569,3,31,2,211668,19124
P14B_GEMM,projection_overlap_1,34,3,569,2,232152,36917
P14B_GEMM,projection_subtract_1,569,3,34,2,232152,21167
P14B_GEMM,projection_overlap_2,34,3,569,2,232152,36750
P14B_GEMM,projection_subtract_2,569,3,34,2,232152,21209
P14B_GEMM,projection_overlap_1,24,5,569,1,136560,20209
P14B_GEMM,projection_subtract_1,569,5,24,1,136560,12709
P14B_GEMM,projection_overlap_2,24,5,569,1,136560,20167
P14B_GEMM,projection_subtract_2,569,5,24,1,136560,12625
P14B_GEMM,projection_overlap_1,29,1,569,2,66004,12750
P14B_GEMM,projection_subtract_1,569,1,29,2,66004,9041
P14B_GEMM,projection_overlap_2,29,1,569,2,66004,12542
P14B_GEMM,projection_subtract_2,569,1,29,2,66004,8832
P14B_GEMM,projection_overlap_1,12,8,569,1,109248,11875
P14B_GEMM,projection_subtract_1,569,8,12,1,109248,9708
P14B_GEMM,projection_overlap_2,12,8,569,1,109248,11750
P14B_GEMM,projection_subtract_2,569,8,12,1,109248,9625
P14B_GEMM,projection_overlap_1,44,2,569,1,100144,22083
P14B_GEMM,projection_subtract_1,569,2,44,1,100144,11166
P14B_GEMM,projection_overlap_2,44,2,569,1,100144,21875
P14B_GEMM,projection_subtract_2,569,2,44,1,100144,11084
P14B_GEMM,projection_overlap_1,46,2,569,2,209392,46750
P14B_GEMM,projection_subtract_1,569,2,46,2,209392,23584
P14B_GEMM,projection_overlap_2,46,2,569,2,209392,46416
P14B_GEMM,projection_subtract_2,569,2,46,2,209392,23375
P14B_GEMM,projection_overlap_1,23,5,569,1,130870,20250
P14B_GEMM,projection_subtract_1,569,5,23,1,130870,12417
P14B_GEMM,projection_overlap_2,23,5,569,1,130870,19958
P14B_GEMM,projection_subtract_2,569,5,23,1,130870,12125
P14B_GEMM,projection_overlap_1,42,3,569,1,143388,22083
P14B_GEMM,projection_subtract_1,569,3,42,1,143388,13750
P14B_GEMM,projection_overlap_2,42,3,569,1,143388,21917
P14B_GEMM,projection_subtract_2,569,3,42,1,143388,13416
P14B_GEMM,projection_overlap_1,36,8,569,1,327744,32125
P14B_GEMM,projection_subtract_1,569,8,36,1,327744,25833
P14B_GEMM,projection_overlap_2,36,8,569,1,327744,31958
P14B_GEMM,projection_subtract_2,569,8,36,1,327744,25875
P14B_GEMM,projection_overlap_1,44,3,569,1,150216,23458
P14B_GEMM,projection_subtract_1,569,3,44,1,150216,14416
P14B_GEMM,projection_overlap_2,44,3,569,1,150216,23042
P14B_GEMM,projection_subtract_2,569,3,44,1,150216,14166
P14B_GEMM,projection_overlap_1,24,10,569,1,273120,30250
P14B_GEMM,projection_subtract_1,569,10,24,1,273120,22875
P14B_GEMM,projection_overlap_2,24,10,569,1,273120,29875
P14B_GEMM,projection_subtract_2,569,10,24,1,273120,22750
P14B_GEMM,projection_overlap_1,34,6,569,1,232152,30541
P14B_GEMM,projection_subtract_1,569,6,34,1,232152,20792
P14B_GEMM,projection_overlap_2,34,6,569,1,232152,30541
P14B_GEMM,projection_subtract_2,569,6,34,1,232152,20708
P14B_GEMM,reduced_hamiltonian,40,40,569,2,3641600,317958
P14B_GEMM,ritz_and_residual,569,12,40,4,2184960,171376
P14B_GEMM,projection_overlap_1,40,1,569,1,45520,8750
P14B_GEMM,projection_subtract_1,569,1,40,1,45520,6208
P14B_GEMM,projection_overlap_2,40,1,569,1,45520,8708
P14B_GEMM,projection_subtract_2,569,1,40,1,45520,6083
P14B_GEMM,projection_overlap_1,36,4,569,1,163872,20083
P14B_GEMM,projection_subtract_1,569,4,36,1,163872,13417
P14B_GEMM,projection_overlap_2,36,4,569,1,163872,19333
P14B_GEMM,projection_subtract_2,569,4,36,1,163872,13291
P14B_GEMM,projection_overlap_1,40,2,569,1,91040,20542
P14B_GEMM,projection_subtract_1,569,2,40,1,91040,10166
P14B_GEMM,projection_overlap_2,40,2,569,1,91040,20459
P14B_GEMM,projection_subtract_2,569,2,40,1,91040,10083
P14B_GEMM,projection_overlap_1,24,4,569,1,109248,13250
P14B_GEMM,projection_subtract_1,569,4,24,1,109248,9291
P14B_GEMM,projection_overlap_2,24,4,569,1,109248,13083
P14B_GEMM,projection_subtract_2,569,4,24,1,109248,9042
P14B_GEMM,projection_overlap_1,28,3,569,1,95592,15791
P14B_GEMM,projection_subtract_1,569,3,28,1,95592,9167
P14B_GEMM,projection_overlap_2,28,3,569,1,95592,15542
P14B_GEMM,projection_subtract_2,569,3,28,1,95592,8709
P14B_GEMM,projection_overlap_1,31,2,569,1,70556,16417
P14B_GEMM,projection_subtract_1,569,2,31,1,70556,8167
P14B_GEMM,projection_overlap_2,31,2,569,1,70556,16292
P14B_GEMM,projection_subtract_2,569,2,31,1,70556,8083
P14B_GEMM,projection_overlap_1,43,3,569,1,146802,22917
P14B_GEMM,projection_subtract_1,569,3,43,1,146802,13833
P14B_GEMM,projection_overlap_2,43,3,569,1,146802,22875
P14B_GEMM,projection_subtract_2,569,3,43,1,146802,13583
P14B_GEMM,projection_overlap_1,36,9,569,1,368712,42583
P14B_GEMM,projection_subtract_1,569,9,36,1,368712,31750
P14B_GEMM,projection_overlap_2,36,9,569,1,368712,42542
P14B_GEMM,projection_subtract_2,569,9,36,1,368712,31542
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "si-elongated"
lattice_unit = "bohr"
lattice_vectors = [[8, 0, 0], [0, 8, 0], [0, 0, 24]]
[[system.atoms]]
species = "Si"
fractional = [0, 0, 0]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.25, 0.125]
[[system.atoms]]
species = "Si"
fractional = [0, 0, 0.5]
[[system.atoms]]
species = "Si"
fractional = [0.25, 0.25, 0.625]
[pseudopotentials]
Si = "SI_UPF_PATH"
[basis]
ecutwfc_ry = 8
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### syn-s3

Expected failure matched: False

```text
P14B_GEMM,reduced_hamiltonian,6,6,7,100,50400,19126
P14B_GEMM,ritz_and_residual,7,6,6,200,100800,17622
P14B_GEMM,projection_overlap_1,6,1,7,23,1932,1127
P14B_GEMM,projection_subtract_1,7,1,6,23,1932,996
P14B_GEMM,projection_overlap_2,6,1,7,23,1932,878
P14B_GEMM,projection_subtract_2,7,1,6,23,1932,1369
P14B_GEMM,reduced_hamiltonian,7,7,7,23,15778,3083
P14B_GEMM,ritz_and_residual,7,6,7,46,27048,4374
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "syn-s3"
lattice_unit = "bohr"
lattice_vectors = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
[[system.atoms]]
species = "Si"
fractional = [0.13, 0.21, 0.32]
[pseudopotentials]
Si = "synthetic.UPF"
[basis]
ecutwfc_ry = 3
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
```

### syn-z14

Expected failure matched: False

```text
P14B_GEMM,reduced_hamiltonian,12,12,33,100,950400,120333
P14B_GEMM,ritz_and_residual,33,12,12,200,1900800,182379
P14B_GEMM,projection_overlap_1,12,12,33,6,57024,6331
P14B_GEMM,projection_subtract_1,33,12,12,6,57024,5584
P14B_GEMM,projection_overlap_2,12,12,33,6,57024,6000
P14B_GEMM,projection_subtract_2,33,12,12,6,57024,5458
P14B_GEMM,reduced_hamiltonian,24,24,33,20,760320,73665
P14B_GEMM,ritz_and_residual,33,12,24,40,760320,65456
P14B_GEMM,projection_overlap_1,12,10,33,12,95040,11791
P14B_GEMM,projection_subtract_1,33,10,12,12,95040,10417
P14B_GEMM,projection_overlap_2,12,10,33,12,95040,10834
P14B_GEMM,projection_subtract_2,33,10,12,12,95040,8624
P14B_GEMM,reduced_hamiltonian,22,22,33,24,766656,81503
P14B_GEMM,ritz_and_residual,33,12,22,48,836352,71417
P14B_GEMM,projection_overlap_1,22,9,33,1,13068,2459
P14B_GEMM,projection_subtract_1,33,9,22,1,13068,1375
P14B_GEMM,projection_overlap_2,22,9,33,1,13068,1541
P14B_GEMM,projection_subtract_2,33,9,22,1,13068,1250
P14B_GEMM,reduced_hamiltonian,31,31,33,6,380556,43124
P14B_GEMM,ritz_and_residual,33,12,31,12,294624,24710
P14B_GEMM,projection_overlap_1,12,3,33,4,9504,1792
P14B_GEMM,projection_subtract_1,33,3,12,4,9504,1542
P14B_GEMM,projection_overlap_2,12,3,33,4,9504,1542
P14B_GEMM,projection_subtract_2,33,3,12,4,9504,917
P14B_GEMM,reduced_hamiltonian,15,15,33,5,74250,8875
P14B_GEMM,ritz_and_residual,33,12,15,10,118800,10959
P14B_GEMM,projection_overlap_1,15,1,33,2,1980,500
P14B_GEMM,projection_subtract_1,33,1,15,2,1980,542
P14B_GEMM,projection_overlap_2,15,1,33,2,1980,374
P14B_GEMM,projection_subtract_2,33,1,15,2,1980,374
P14B_GEMM,reduced_hamiltonian,16,16,33,10,168960,18253
P14B_GEMM,ritz_and_residual,33,12,16,20,253440,23373
P14B_GEMM,projection_overlap_1,16,1,33,2,2112,458
P14B_GEMM,projection_subtract_1,33,1,16,2,2112,375
P14B_GEMM,projection_overlap_2,16,1,33,2,2112,416
P14B_GEMM,projection_subtract_2,33,1,16,2,2112,332
P14B_GEMM,reduced_hamiltonian,17,17,33,9,171666,19167
P14B_GEMM,ritz_and_residual,33,12,17,18,242352,25415
P14B_GEMM,projection_overlap_1,12,11,33,8,69696,8165
P14B_GEMM,projection_subtract_1,33,11,12,8,69696,7083
P14B_GEMM,projection_overlap_2,12,11,33,8,69696,7751
P14B_GEMM,projection_subtract_2,33,11,12,8,69696,6792
P14B_GEMM,reduced_hamiltonian,23,23,33,15,523710,52749
P14B_GEMM,ritz_and_residual,33,12,23,30,546480,47247
P14B_GEMM,projection_overlap_1,23,3,33,5,22770,4000
P14B_GEMM,projection_subtract_1,33,3,23,5,22770,2333
P14B_GEMM,projection_overlap_2,23,3,33,5,22770,3501
P14B_GEMM,projection_subtract_2,33,3,23,5,22770,2209
P14B_GEMM,reduced_hamiltonian,26,26,33,23,1026168,103582
P14B_GEMM,ritz_and_residual,33,12,26,46,947232,80297
P14B_GEMM,projection_overlap_1,26,2,33,4,13728,2875
P14B_GEMM,projection_subtract_1,33,2,26,4,13728,1834
P14B_GEMM,projection_overlap_2,26,2,33,4,13728,2584
P14B_GEMM,projection_subtract_2,33,2,26,4,13728,1624
P14B_GEMM,reduced_hamiltonian,28,28,33,16,827904,79294
P14B_GEMM,ritz_and_residual,33,12,28,32,709632,59752
P14B_GEMM,projection_overlap_1,28,2,33,3,11088,2333
P14B_GEMM,projection_subtract_1,33,2,28,3,11088,1333
P14B_GEMM,projection_overlap_2,28,2,33,3,11088,2040
P14B_GEMM,projection_subtract_2,33,2,28,3,11088,1374
P14B_GEMM,reduced_hamiltonian,30,30,33,8,475200,45082
P14B_GEMM,ritz_and_residual,33,12,30,16,380160,31544
P14B_GEMM,projection_overlap_1,12,5,33,6,23760,4208
P14B_GEMM,projection_subtract_1,33,5,12,6,23760,2750
P14B_GEMM,projection_overlap_2,12,5,33,6,23760,3709
P14B_GEMM,projection_subtract_2,33,5,12,6,23760,2376
P14B_GEMM,projection_overlap_1,17,4,33,4,17952,2584
P14B_GEMM,projection_subtract_1,33,4,17,4,17952,1708
P14B_GEMM,projection_overlap_2,17,4,33,4,17952,2458
P14B_GEMM,projection_subtract_2,33,4,17,4,17952,1667
P14B_GEMM,reduced_hamiltonian,21,21,33,21,611226,61418
P14B_GEMM,ritz_and_residual,33,12,21,42,698544,60707
P14B_GEMM,projection_overlap_1,21,4,33,4,22176,3084
P14B_GEMM,projection_subtract_1,33,4,21,4,22176,2000
P14B_GEMM,projection_overlap_2,21,4,33,4,22176,2876
P14B_GEMM,projection_subtract_2,33,4,21,4,22176,2001
P14B_GEMM,reduced_hamiltonian,25,25,33,21,866250,87209
P14B_GEMM,ritz_and_residual,33,12,25,42,831600,70042
P14B_GEMM,projection_overlap_1,12,9,33,11,78408,10291
P14B_GEMM,projection_subtract_1,33,9,12,11,78408,7751
P14B_GEMM,projection_overlap_2,12,9,33,11,78408,9792
P14B_GEMM,projection_subtract_2,33,9,12,11,78408,7416
P14B_GEMM,projection_overlap_1,21,5,33,5,34650,5751
P14B_GEMM,projection_subtract_1,33,5,21,5,34650,3709
P14B_GEMM,projection_overlap_2,21,5,33,5,34650,5209
P14B_GEMM,projection_subtract_2,33,5,21,5,34650,3584
P14B_GEMM,projection_overlap_1,26,1,33,7,12012,2332
P14B_GEMM,projection_subtract_1,33,1,26,7,12012,2541
P14B_GEMM,projection_overlap_2,26,1,33,7,12012,2376
P14B_GEMM,projection_subtract_2,33,1,26,7,12012,2166
P14B_GEMM,reduced_hamiltonian,27,27,33,14,673596,63957
P14B_GEMM,ritz_and_residual,33,12,27,28,598752,49707
P14B_GEMM,projection_overlap_1,22,5,33,4,29040,4833
P14B_GEMM,projection_subtract_1,33,5,22,4,29040,3166
P14B_GEMM,projection_overlap_2,22,5,33,4,29040,4416
P14B_GEMM,projection_subtract_2,33,5,22,4,29040,3167
P14B_GEMM,projection_overlap_1,27,2,33,2,7128,1584
P14B_GEMM,projection_subtract_1,33,2,27,2,7128,876
P14B_GEMM,projection_overlap_2,27,2,33,2,7128,1374
P14B_GEMM,projection_subtract_2,33,2,27,2,7128,916
P14B_GEMM,reduced_hamiltonian,29,29,33,9,499554,48625
P14B_GEMM,ritz_and_residual,33,12,29,18,413424,34331
P14B_GEMM,projection_overlap_1,12,8,33,20,126720,14999
P14B_GEMM,projection_subtract_1,33,8,12,20,126720,11290
P14B_GEMM,projection_overlap_2,12,8,33,20,126720,13918
P14B_GEMM,projection_subtract_2,33,8,12,20,126720,11085
P14B_GEMM,reduced_hamiltonian,20,20,33,23,607200,60375
P14B_GEMM,ritz_and_residual,33,12,20,46,728640,63626
P14B_GEMM,projection_overlap_1,20,3,33,1,3960,708
P14B_GEMM,projection_subtract_1,33,3,20,1,3960,417
P14B_GEMM,projection_overlap_2,20,3,33,1,3960,625
P14B_GEMM,projection_subtract_2,33,3,20,1,3960,417
P14B_GEMM,projection_overlap_1,23,1,33,3,4554,917
P14B_GEMM,projection_subtract_1,33,1,23,3,4554,917
P14B_GEMM,projection_overlap_2,23,1,33,3,4554,875
P14B_GEMM,projection_subtract_2,33,1,23,3,4554,709
P14B_GEMM,projection_overlap_1,12,1,33,2,1584,375
P14B_GEMM,projection_subtract_1,33,1,12,2,1584,292
P14B_GEMM,projection_overlap_2,12,1,33,2,1584,333
P14B_GEMM,projection_subtract_2,33,1,12,2,1584,292
P14B_GEMM,reduced_hamiltonian,13,13,33,2,22308,2792
P14B_GEMM,ritz_and_residual,33,12,13,4,41184,3958
P14B_GEMM,projection_overlap_1,25,3,33,2,9900,1709
P14B_GEMM,projection_subtract_1,33,3,25,2,9900,1083
P14B_GEMM,projection_overlap_2,25,3,33,2,9900,1542
P14B_GEMM,projection_subtract_2,33,3,25,2,9900,1042
P14B_GEMM,projection_overlap_1,28,1,33,4,7392,1375
P14B_GEMM,projection_subtract_1,33,1,28,4,7392,1374
P14B_GEMM,projection_overlap_2,28,1,33,4,7392,1459
P14B_GEMM,projection_subtract_2,33,1,28,4,7392,1167
P14B_GEMM,projection_overlap_1,29,1,33,2,3828,750
P14B_GEMM,projection_subtract_1,33,1,29,2,3828,625
P14B_GEMM,projection_overlap_2,29,1,33,2,3828,750
P14B_GEMM,projection_subtract_2,33,1,29,2,3828,708
P14B_GEMM,projection_overlap_1,30,1,33,2,3960,792
P14B_GEMM,projection_subtract_1,33,1,30,2,3960,834
P14B_GEMM,projection_overlap_2,30,1,33,2,3960,750
P14B_GEMM,projection_subtract_2,33,1,30,2,3960,792
P14B_GEMM,projection_overlap_1,27,3,33,2,10692,1792
P14B_GEMM,projection_subtract_1,33,3,27,2,10692,1167
P14B_GEMM,projection_overlap_2,27,3,33,2,10692,1626
P14B_GEMM,projection_subtract_2,33,3,27,2,10692,1082
P14B_GEMM,projection_overlap_1,12,7,33,13,72072,9376
P14B_GEMM,projection_subtract_1,33,7,12,13,72072,6916
P14B_GEMM,projection_overlap_2,12,7,33,13,72072,8625
P14B_GEMM,projection_subtract_2,33,7,12,13,72072,6748
P14B_GEMM,reduced_hamiltonian,19,19,33,15,357390,37833
P14B_GEMM,ritz_and_residual,33,12,19,30,451440,39498
P14B_GEMM,projection_overlap_1,19,4,33,3,15048,2167
P14B_GEMM,projection_subtract_1,33,4,19,3,15048,1375
P14B_GEMM,projection_overlap_2,19,4,33,3,15048,2042
P14B_GEMM,projection_subtract_2,33,4,19,3,15048,1334
P14B_GEMM,projection_overlap_1,12,2,33,2,3168,750
P14B_GEMM,projection_subtract_1,33,2,12,2,3168,376
P14B_GEMM,projection_overlap_2,12,2,33,2,3168,666
P14B_GEMM,projection_subtract_2,33,2,12,2,3168,375
P14B_GEMM,reduced_hamiltonian,14,14,33,3,38808,4542
P14B_GEMM,ritz_and_residual,33,12,14,6,66528,6832
P14B_GEMM,projection_overlap_1,14,2,33,2,3696,876
P14B_GEMM,projection_subtract_1,33,2,14,2,3696,999
P14B_GEMM,projection_overlap_2,14,2,33,2,3696,750
P14B_GEMM,projection_subtract_2,33,2,14,2,3696,458
P14B_GEMM,projection_overlap_1,20,5,33,6,39600,6416
P14B_GEMM,projection_subtract_1,33,5,20,6,39600,4167
P14B_GEMM,projection_overlap_2,20,5,33,6,39600,6042
P14B_GEMM,projection_subtract_2,33,5,20,6,39600,4041
P14B_GEMM,projection_overlap_1,24,3,33,1,4752,791
P14B_GEMM,projection_subtract_1,33,3,24,1,4752,458
P14B_GEMM,projection_overlap_2,24,3,33,1,4752,708
P14B_GEMM,projection_subtract_2,33,3,24,1,4752,458
P14B_GEMM,projection_overlap_1,22,4,33,1,5808,833
P14B_GEMM,projection_subtract_1,33,4,22,1,5808,583
P14B_GEMM,projection_overlap_2,22,4,33,1,5808,750
P14B_GEMM,projection_subtract_2,33,4,22,1,5808,542
P14B_GEMM,projection_overlap_1,23,5,33,2,15180,2500
P14B_GEMM,projection_subtract_1,33,5,23,2,15180,1708
P14B_GEMM,projection_overlap_2,23,5,33,2,15180,2292
P14B_GEMM,projection_subtract_2,33,5,23,2,15180,1666
P14B_GEMM,projection_overlap_1,28,3,33,1,5544,916
P14B_GEMM,projection_subtract_1,33,3,28,1,5544,584
P14B_GEMM,projection_overlap_2,28,3,33,1,5544,792
P14B_GEMM,projection_subtract_2,33,3,28,1,5544,584
P14B_GEMM,projection_overlap_1,22,6,33,3,26136,3459
P14B_GEMM,projection_subtract_1,33,6,22,3,26136,2793
P14B_GEMM,projection_overlap_2,22,6,33,3,26136,3208
P14B_GEMM,projection_subtract_2,33,6,22,3,26136,2541
P14B_GEMM,projection_overlap_1,16,2,33,2,4224,958
P14B_GEMM,projection_subtract_1,33,2,16,2,4224,916
P14B_GEMM,projection_overlap_2,16,2,33,2,4224,833
P14B_GEMM,projection_subtract_2,33,2,16,2,4224,500
P14B_GEMM,reduced_hamiltonian,18,18,33,11,235224,24208
P14B_GEMM,ritz_and_residual,33,12,18,22,313632,28122
P14B_GEMM,projection_overlap_1,12,4,33,6,19008,2707
P14B_GEMM,projection_subtract_1,33,4,12,6,19008,1832
P14B_GEMM,projection_overlap_2,12,4,33,6,19008,2501
P14B_GEMM,projection_subtract_2,33,4,12,6,19008,1832
P14B_GEMM,projection_overlap_1,16,4,33,2,8448,1208
P14B_GEMM,projection_subtract_1,33,4,16,2,8448,916
P14B_GEMM,projection_overlap_2,16,4,33,2,8448,1126
P14B_GEMM,projection_subtract_2,33,4,16,2,8448,791
P14B_GEMM,projection_overlap_1,21,1,33,2,2772,583
P14B_GEMM,projection_subtract_1,33,1,21,2,2772,1083
P14B_GEMM,projection_overlap_2,21,1,33,2,2772,541
P14B_GEMM,projection_subtract_2,33,1,21,2,2772,625
P14B_GEMM,projection_overlap_1,24,7,33,2,22176,2583
P14B_GEMM,projection_subtract_1,33,7,24,2,22176,2166
P14B_GEMM,projection_overlap_2,24,7,33,2,22176,2417
P14B_GEMM,projection_subtract_2,33,7,24,2,22176,2041
P14B_GEMM,projection_overlap_1,31,1,33,1,2046,375
P14B_GEMM,projection_subtract_1,33,1,31,1,2046,291
P14B_GEMM,projection_overlap_2,31,1,33,1,2046,417
P14B_GEMM,projection_subtract_2,33,1,31,1,2046,291
P14B_GEMM,reduced_hamiltonian,32,32,33,3,202752,19166
P14B_GEMM,ritz_and_residual,33,12,32,6,152064,12541
P14B_GEMM,projection_overlap_1,21,7,33,2,19404,2375
P14B_GEMM,projection_subtract_1,33,7,21,2,19404,1958
P14B_GEMM,projection_overlap_2,21,7,33,2,19404,2167
P14B_GEMM,projection_subtract_2,33,7,21,2,19404,1876
P14B_GEMM,projection_overlap_1,23,4,33,1,6072,875
P14B_GEMM,projection_subtract_1,33,4,23,1,6072,583
P14B_GEMM,projection_overlap_2,23,4,33,1,6072,791
P14B_GEMM,projection_subtract_2,33,4,23,1,6072,583
P14B_GEMM,projection_overlap_1,27,1,33,3,5346,1042
P14B_GEMM,projection_subtract_1,33,1,27,3,5346,833
P14B_GEMM,projection_overlap_2,27,1,33,3,5346,959
P14B_GEMM,projection_subtract_2,33,1,27,3,5346,792
P14B_GEMM,projection_overlap_1,24,8,33,2,25344,2708
P14B_GEMM,projection_subtract_1,33,8,24,2,25344,2292
P14B_GEMM,projection_overlap_2,24,8,33,2,25344,2500
P14B_GEMM,projection_subtract_2,33,8,24,2,25344,2208
P14B_GEMM,projection_overlap_1,23,6,33,1,9108,1167
P14B_GEMM,projection_subtract_1,33,6,23,1,9108,959
P14B_GEMM,projection_overlap_2,23,6,33,1,9108,1125
P14B_GEMM,projection_subtract_2,33,6,23,1,9108,917
P14B_GEMM,projection_overlap_1,12,6,33,6,28512,4126
P14B_GEMM,projection_subtract_1,33,6,12,6,28512,2834
P14B_GEMM,projection_overlap_2,12,6,33,6,28512,3875
P14B_GEMM,projection_subtract_2,33,6,12,6,28512,2707
P14B_GEMM,projection_overlap_1,18,6,33,1,7128,1000
P14B_GEMM,projection_subtract_1,33,6,18,1,7128,708
P14B_GEMM,projection_overlap_2,18,6,33,1,7128,917
P14B_GEMM,projection_subtract_2,33,6,18,1,7128,667
P14B_GEMM,projection_overlap_1,24,1,33,7,11088,2208
P14B_GEMM,projection_subtract_1,33,1,24,7,11088,2208
P14B_GEMM,projection_overlap_2,24,1,33,7,11088,2251
P14B_GEMM,projection_subtract_2,33,1,24,7,11088,2459
P14B_GEMM,projection_overlap_1,25,1,33,6,9900,1877
P14B_GEMM,projection_subtract_1,33,1,25,6,9900,1960
P14B_GEMM,projection_overlap_2,25,1,33,6,9900,1915
P14B_GEMM,projection_subtract_2,33,1,25,6,9900,2291
P14B_GEMM,projection_overlap_1,21,3,33,1,4158,791
P14B_GEMM,projection_subtract_1,33,3,21,1,4158,416
P14B_GEMM,projection_overlap_2,21,3,33,1,4158,667
P14B_GEMM,projection_subtract_2,33,3,21,1,4158,417
P14B_GEMM,projection_overlap_1,20,2,33,4,10560,2332
P14B_GEMM,projection_subtract_1,33,2,20,4,10560,1334
P14B_GEMM,projection_overlap_2,20,2,33,4,10560,2126
P14B_GEMM,projection_subtract_2,33,2,20,4,10560,1250
P14B_GEMM,projection_overlap_1,22,2,33,2,5808,1292
P14B_GEMM,projection_subtract_1,33,2,22,2,5808,750
P14B_GEMM,projection_overlap_2,22,2,33,2,5808,1209
P14B_GEMM,projection_subtract_2,33,2,22,2,5808,666
P14B_GEMM,projection_overlap_1,16,3,33,1,3168,583
P14B_GEMM,projection_subtract_1,33,3,16,1,3168,834
P14B_GEMM,projection_overlap_2,16,3,33,1,3168,500
P14B_GEMM,projection_subtract_2,33,3,16,1,3168,333
P14B_GEMM,projection_overlap_1,19,1,33,1,1254,292
P14B_GEMM,projection_subtract_1,33,1,19,1,1254,375
P14B_GEMM,projection_overlap_2,19,1,33,1,1254,250
P14B_GEMM,projection_subtract_2,33,1,19,1,1254,500
P14B_GEMM,projection_overlap_1,22,8,33,1,11616,1333
P14B_GEMM,projection_subtract_1,33,8,22,1,11616,1042
P14B_GEMM,projection_overlap_2,22,8,33,1,11616,1250
P14B_GEMM,projection_subtract_2,33,8,22,1,11616,1042
P14B_GEMM,projection_overlap_1,19,7,33,2,17556,2209
P14B_GEMM,projection_subtract_1,33,7,19,2,17556,1791
P14B_GEMM,projection_overlap_2,19,7,33,2,17556,2042
P14B_GEMM,projection_subtract_2,33,7,19,2,17556,1708
P14B_GEMM,projection_overlap_1,15,2,33,1,1980,458
P14B_GEMM,projection_subtract_1,33,2,15,1,1980,583
P14B_GEMM,projection_overlap_2,15,2,33,1,1980,417
P14B_GEMM,projection_subtract_2,33,2,15,1,1980,250
P14B_GEMM,projection_overlap_1,17,1,33,1,1122,208
P14B_GEMM,projection_subtract_1,33,1,17,1,1122,208
P14B_GEMM,projection_overlap_2,17,1,33,1,1122,208
P14B_GEMM,projection_subtract_2,33,1,17,1,1122,167
P14B_GEMM,projection_overlap_1,20,6,33,3,23760,3251
P14B_GEMM,projection_subtract_1,33,6,20,3,23760,2333
P14B_GEMM,projection_overlap_2,20,6,33,3,23760,3083
P14B_GEMM,projection_subtract_2,33,6,20,3,23760,2250
P14B_GEMM,projection_overlap_1,17,5,33,2,11220,1958
P14B_GEMM,projection_subtract_1,33,5,17,2,11220,1333
P14B_GEMM,projection_overlap_2,17,5,33,2,11220,1750
P14B_GEMM,projection_subtract_2,33,5,17,2,11220,1125
P14B_GEMM,projection_overlap_1,22,1,33,3,4356,915
P14B_GEMM,projection_subtract_1,33,1,22,3,4356,750
P14B_GEMM,projection_overlap_2,22,1,33,3,4356,834
P14B_GEMM,projection_subtract_2,33,1,22,3,4356,709
P14B_GEMM,projection_overlap_1,15,3,33,2,5940,1167
P14B_GEMM,projection_subtract_1,33,3,15,2,5940,875
P14B_GEMM,projection_overlap_2,15,3,33,2,5940,1000
P14B_GEMM,projection_subtract_2,33,3,15,2,5940,625
P14B_GEMM,projection_overlap_1,18,3,33,3,10692,1875
P14B_GEMM,projection_subtract_1,33,3,18,3,10692,1125
P14B_GEMM,projection_overlap_2,18,3,33,3,10692,1666
P14B_GEMM,projection_subtract_2,33,3,18,3,10692,1083
P14B_GEMM,projection_overlap_1,19,2,33,1,2508,625
P14B_GEMM,projection_subtract_1,33,2,19,1,2508,292
P14B_GEMM,projection_overlap_2,19,2,33,1,2508,541
P14B_GEMM,projection_subtract_2,33,2,19,1,2508,291
P14B_GEMM,projection_overlap_1,22,7,33,2,20328,2500
P14B_GEMM,projection_subtract_1,33,7,22,2,20328,2083
P14B_GEMM,projection_overlap_2,22,7,33,2,20328,2333
P14B_GEMM,projection_subtract_2,33,7,22,2,20328,1959
P14B_GEMM,projection_overlap_1,20,4,33,5,26400,3708
P14B_GEMM,projection_subtract_1,33,4,20,5,26400,2458
P14B_GEMM,projection_overlap_2,20,4,33,5,26400,3418
P14B_GEMM,projection_subtract_2,33,4,20,5,26400,2375
P14B_GEMM,projection_overlap_1,13,1,33,1,858,208
P14B_GEMM,projection_subtract_1,33,1,13,1,858,167
P14B_GEMM,projection_overlap_2,13,1,33,1,858,208
P14B_GEMM,projection_subtract_2,33,1,13,1,858,333
P14B_GEMM,projection_overlap_1,14,1,33,1,924,208
P14B_GEMM,projection_subtract_1,33,1,14,1,924,167
P14B_GEMM,projection_overlap_2,14,1,33,1,924,209
P14B_GEMM,projection_subtract_2,33,1,14,1,924,166
P14B_GEMM,projection_overlap_1,19,5,33,2,12540,2042
P14B_GEMM,projection_subtract_1,33,5,19,2,12540,1333
P14B_GEMM,projection_overlap_2,19,5,33,2,12540,2000
P14B_GEMM,projection_subtract_2,33,5,19,2,12540,1250
P14B_GEMM,projection_overlap_1,18,7,33,2,16632,2084
P14B_GEMM,projection_subtract_1,33,7,18,2,16632,1832
P14B_GEMM,projection_overlap_2,18,7,33,2,16632,1958
P14B_GEMM,projection_subtract_2,33,7,18,2,16632,1625
P14B_GEMM,projection_overlap_1,19,3,33,3,11286,1959
P14B_GEMM,projection_subtract_1,33,3,19,3,11286,1167
P14B_GEMM,projection_overlap_2,19,3,33,3,11286,1792
P14B_GEMM,projection_subtract_2,33,3,19,3,11286,1167
P14B_GEMM,projection_overlap_1,18,1,33,1,1188,250
P14B_GEMM,projection_subtract_1,33,1,18,1,1188,208
P14B_GEMM,projection_overlap_2,18,1,33,1,1188,250
P14B_GEMM,projection_subtract_2,33,1,18,1,1188,208
P14B_GEMM,projection_overlap_1,19,6,33,2,15048,2167
P14B_GEMM,projection_subtract_1,33,6,19,2,15048,1500
P14B_GEMM,projection_overlap_2,19,6,33,2,15048,2083
P14B_GEMM,projection_subtract_2,33,6,19,2,15048,1417
P14B_GEMM,projection_overlap_1,25,2,33,1,3300,791
P14B_GEMM,projection_subtract_1,33,2,25,1,3300,500
P14B_GEMM,projection_overlap_2,25,2,33,1,3300,667
P14B_GEMM,projection_subtract_2,33,2,25,1,3300,417
P14B_GEMM,projection_overlap_1,20,1,33,1,1320,250
P14B_GEMM,projection_subtract_1,33,1,20,1,1320,250
P14B_GEMM,projection_overlap_2,20,1,33,1,1320,250
P14B_GEMM,projection_subtract_2,33,1,20,1,1320,250
P14B_GEMM,projection_overlap_1,16,5,33,1,5280,875
P14B_GEMM,projection_subtract_1,33,5,16,1,5280,750
P14B_GEMM,projection_overlap_2,16,5,33,1,5280,834
P14B_GEMM,projection_subtract_2,33,5,16,1,5280,541
P14B_GEMM,projection_overlap_1,18,4,33,1,4752,708
P14B_GEMM,projection_subtract_1,33,4,18,1,4752,459
P14B_GEMM,projection_overlap_2,18,4,33,1,4752,625
P14B_GEMM,projection_subtract_2,33,4,18,1,4752,459
P14B_GEMM,projection_overlap_1,24,2,33,1,3168,667
P14B_GEMM,projection_subtract_1,33,2,24,1,3168,416
P14B_GEMM,projection_overlap_2,24,2,33,1,3168,625
P14B_GEMM,projection_subtract_2,33,2,24,1,3168,416
error: NOT CONVERGED: iteration budget exhausted
```

```toml
[system]
name = "syn-z14"
lattice_unit = "bohr"
lattice_vectors = [[4, 0, 0], [0, 4, 0], [0, 0, 4]]
[[system.atoms]]
species = "Si"
fractional = [0.13, 0.21, 0.32]
[pseudopotentials]
Si = "synthetic.UPF"
[basis]
ecutwfc_ry = 12
kpoints = [[0.0, 0.0, 0.0]]
[xc]
functional = "lda-pz"
[scf]
seed = 1
threads = 1
max_iterations = 100
mixing = "pulay"
mixing_beta = 0.5
energy_tolerance_ha = 1e-8
density_residual_tolerance = 1e-7
nbands = 12
```
