# Benchmark suite results

Synthetic results are non-physical. Timings are indicative on a shared host.

| Case | Class | Status | Seconds (median) | Iterations | Energy (Ha) |
|---|---|---|---:|---:|---:|
| si-elongated | elongated | nonconverged | 0.513 | 100 | -15.08247814257598 |
| si-gamma | si-primitive | pass | 0.130 | 10 | -7.29553046821076 |
| si-k222 | si-primitive | pass | 0.779 | 10 | -7.831230057121135 |
| si-k222s | si-primitive | pass | 0.413 | 10 | -7.9194234124690865 |
| si-low | si-cutoff | pass | 0.022 | 11 | -7.233756452824804 |
| si-skew | skewed | pass | 0.152 | 13 | -7.25055567296079 |
| si8-gamma | si-conventional | pass | 2.887 | 15 | -31.32180360103127 |
| syn-cubic-degenerate | synthetic-gap | pass | 0.011 | 10 | -4.261277234135343 |
| syn-d3 | synthetic-projectors | pass | 0.006 | 8 | -0.9465076331236224 |
| syn-fft-primes | awkward-fft | pass | 0.008 | 16 | -0.7346843876893185 |
| syn-gaussian-fcc | synthetic-gap | pass | 0.048 | 17 | -4.683023405363711 |
| syn-odd-shift | synthetic-mesh | pass | 0.011 | 9 | -1.5900664019906934 |
| syn-p3 | synthetic-projectors | pass | 0.007 | 10 | -0.9465076597173745 |
| syn-s1 | synthetic-projectors | pass | 0.006 | 14 | -0.7346754280328636 |
| syn-s2 | synthetic-projectors | pass | 0.006 | 16 | 0.6304095283333488 |
| syn-s3 | synthetic-projectors | nonconverged | 0.011 | 100 | 1.8154474817076043 |
| syn-spd | synthetic-projectors | pass | 0.007 | 11 | -0.7346751923893944 |
| syn-z1 | synthetic-valence | pass | 0.009 | 7 | -0.9291990923944964 |
| syn-z14 | synthetic-valence | nonconverged | 0.047 | 100 | -71.14987536499353 |

## By system class

| Class | Cases | Passed | Sum median seconds |
|---|---:|---:|---:|
| awkward-fft | 1 | 1 | 0.008 |
| elongated | 1 | 0 | 0.513 |
| si-conventional | 1 | 1 | 2.887 |
| si-cutoff | 1 | 1 | 0.022 |
| si-primitive | 3 | 3 | 1.323 |
| skewed | 1 | 1 | 0.152 |
| synthetic-gap | 2 | 2 | 0.059 |
| synthetic-mesh | 1 | 1 | 0.011 |
| synthetic-projectors | 6 | 5 | 0.043 |
| synthetic-valence | 2 | 1 | 0.056 |

## Baseline comparison

Only selected cases are compared; missing selected entries fail comparison.

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
