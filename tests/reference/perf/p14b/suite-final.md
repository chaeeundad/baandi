# Benchmark suite results

Synthetic results are non-physical. Timings are indicative on a shared host.

| Case | Class | Status | Seconds (median) | Iterations | Energy (Ha) |
|---|---|---|---:|---:|---:|
| si-elongated | elongated | nonconverged | 0.573 | 100 | -15.128518095558501 |
| si-gamma | si-primitive | pass | 0.159 | 10 | -7.295530468210765 |
| si-k222 | si-primitive | pass | 0.920 | 10 | -7.831230057121137 |
| si-k222s | si-primitive | pass | 0.446 | 10 | -7.919423412469083 |
| si-low | si-cutoff | pass | 0.027 | 11 | -7.233756452824798 |
| si-skew | skewed | pass | 0.157 | 13 | -7.250555672960788 |
| si8-gamma | si-conventional | pass | 4.287 | 15 | -31.321803601031256 |
| syn-cubic-degenerate | synthetic-gap | pass | 0.012 | 10 | -4.2612772341353455 |
| syn-d3 | synthetic-projectors | pass | 0.007 | 8 | -0.9465076331236242 |
| syn-fft-primes | awkward-fft | pass | 0.007 | 14 | -0.7346843876893185 |
| syn-gaussian-fcc | synthetic-gap | pass | 0.047 | 17 | -4.68302340536371 |
| syn-odd-shift | synthetic-mesh | pass | 0.011 | 9 | -1.590066401990696 |
| syn-p3 | synthetic-projectors | pass | 0.007 | 8 | -0.9465076597173674 |
| syn-s1 | synthetic-projectors | pass | 0.005 | 13 | -0.7346754280328689 |
| syn-s2 | synthetic-projectors | pass | 0.006 | 17 | 0.6304095283333497 |
| syn-s3 | synthetic-projectors | nonconverged | 0.011 | 100 | 1.7660002310813194 |
| syn-spd | synthetic-projectors | pass | 0.006 | 11 | -0.7346751923893944 |
| syn-z1 | synthetic-valence | pass | 0.009 | 7 | -0.9291990923944965 |
| syn-z14 | synthetic-valence | nonconverged | 0.049 | 100 | -71.13959832030042 |

## By system class

| Class | Cases | Passed | Sum median seconds |
|---|---:|---:|---:|
| awkward-fft | 1 | 1 | 0.007 |
| elongated | 1 | 0 | 0.573 |
| si-conventional | 1 | 1 | 4.287 |
| si-cutoff | 1 | 1 | 0.027 |
| si-primitive | 3 | 3 | 1.525 |
| skewed | 1 | 1 | 0.157 |
| synthetic-gap | 2 | 2 | 0.059 |
| synthetic-mesh | 1 | 1 | 0.011 |
| synthetic-projectors | 6 | 5 | 0.042 |
| synthetic-valence | 2 | 1 | 0.058 |

## Baseline comparison

Only selected cases are compared; missing selected entries fail comparison.
- si-elongated: failure message changed
- si-gamma: unchanged
- si-k222: unchanged
- si-k222s: unchanged
- si-low: unchanged
- si-skew: unchanged
- si8-gamma: unchanged
- syn-cubic-degenerate: unchanged
- syn-d3: unchanged
- syn-fft-primes: unchanged
- syn-gaussian-fcc: unchanged
- syn-odd-shift: unchanged
- syn-p3: unchanged
- syn-s1: unchanged
- syn-s2: unchanged
- syn-s3: failure message changed
- syn-spd: unchanged
- syn-z1: unchanged
- syn-z14: failure message changed

## Findings and exact reproducing inputs

Replace SI_UPF_PATH with the read-only validation file. For synthetic cases regenerate synthetic.UPF using the case parameters in JSON. Runtime directories contain executable inputs and logs.

### si-elongated

Expected failure matched: False

```text
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
