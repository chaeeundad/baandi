#!/usr/bin/env bash
# Run by Claude on the quiet dev08 node, from the repository root.
set -euo pipefail
: "${MKLROOT:?Set MKLROOT=/opt/intel/oneapi/mkl/2026.1}"
: "${P10B_UPF:?Absolute path to the existing dev08 Si.pz-vbc.UPF; do not transfer it}"
P10B_OUT="${P10B_OUT:-$PWD/target/p10b-evidence}"
P10B_CPUS="${P10B_CPUS:-6-13}"
P10B_CORE="${P10B_CORE:-6}"
mkdir -p "$P10B_OUT"
export LD_LIBRARY_PATH="$MKLROOT/lib/intel64:$MKLROOT/lib${LD_LIBRARY_PATH:+:$LD_LIBRARY_PATH}"
export MKL_NUM_THREADS=1 OMP_NUM_THREADS=1
{
    date -u
    printf "CPU set=%s kernel core=%s memory budget MiB=%s repetitions=%s\n" "$P10B_CPUS" "$P10B_CORE" "${MKL_BENCH_MAX_MIB:-1024}" "${MKL_BENCH_REPS:-10}"
    git rev-parse HEAD
    git diff --stat
    rustc -Vv
    lscpu
    sha256sum "$P10B_UPF"
    cat "$MKLROOT/include/mkl_version.h"
} > "$P10B_OUT/environment.txt"
cc -std=c11 -Wall -Wextra -Werror -I"$MKLROOT/include" \
    crates/dft-mkl/tests/abi_probe.c -o "$P10B_OUT/abi_probe"
"$P10B_OUT/abi_probe" | tee "$P10B_OUT/abi.log"
# Default gates never enable the proprietary runtime, even with MKLROOT set.
scripts/check.sh | tee "$P10B_OUT/default-check.log"
cargo test --release --workspace --features dft-numerics/mkl -- --nocapture \
    > "$P10B_OUT/mkl-tests.log" 2>&1
cargo test --release --workspace --features dft-numerics/mkl -- --ignored --nocapture \
    > "$P10B_OUT/mkl-ignored.log" 2>&1
cargo clippy --workspace --all-targets --features dft-numerics/mkl -- -D warnings \
    > "$P10B_OUT/mkl-clippy.log" 2>&1
cargo build --release -p dft-mkl --features mkl --example fft_compare
for trial in 1 2 3; do
    taskset -c "$P10B_CORE" target/release/examples/fft_compare \
        > "$P10B_OUT/fft-$trial.log" 2>&1
done
for backend in rust mkl; do
    feature=()
    if [[ "$backend" == mkl ]]; then feature=(--features dft-numerics/mkl); fi
    CARGO_TARGET_DIR="target/p10b-$backend" cargo build --release -p dft-cli -p dft-numerics "${feature[@]}"
    for trial in 1 2 3; do
        CARGO_TARGET_DIR="target/p10b-$backend" taskset -c "$P10B_CORE" cargo run --release -p dft-numerics \
            --example lib_bench "${feature[@]}" > "$P10B_OUT/kernels-$backend-$trial.log" 2>&1
        CARGO_TARGET_DIR="target/p10b-$backend" taskset -c "$P10B_CPUS" cargo test \
            --release -p dft-physics -p dft-numerics "${feature[@]}" --test hamiltonian conventional_cell_block_timing \
            -- --ignored --nocapture > "$P10B_OUT/hamiltonian-$backend-$trial.log" 2>&1
    done
    compare=()
    if [[ "$backend" == mkl ]]; then compare=(--baseline "$P10B_OUT/scf-rust.json"); fi
    taskset -c "$P10B_CORE" python3 -I scripts/bench_scf.py \
        --binary "$PWD/target/p10b-$backend/release/dft-cli" --upf "$P10B_UPF" \
        --cases gamma,k222,k222s,k444s,si8 --threads 1 --repeats 3 \
        --record "$P10B_OUT/scf-$backend.json" "${compare[@]}" \
        | tee "$P10B_OUT/scf-$backend.log"
done
for threads in 2 4 8; do
    taskset -c "$P10B_CPUS" python3 -I scripts/bench_scf.py \
        --binary "$PWD/target/p10b-mkl/release/dft-cli" --upf "$P10B_UPF" \
        --cases gamma,k222,k222s,k444s,si8 --threads "$threads" --repeats 3 \
        --baseline "$P10B_OUT/scf-mkl.json" --require-identical \
        --record "$P10B_OUT/scf-mkl-t$threads.json" | tee "$P10B_OUT/scf-mkl-t$threads.log"
done
# Retain fresh MKL CLI results for real QE comparisons. The benchmark's temporary
# outputs are deleted, so comparing old stored Rust results would prove nothing.
python3 -I - "$P10B_UPF" "$P10B_OUT" <<'PY'
import importlib.util, pathlib, sys
spec = importlib.util.spec_from_file_location('bench', 'scripts/bench_scf.py')
bench = importlib.util.module_from_spec(spec)
spec.loader.exec_module(bench)
for case in ['gamma', 'k222', 'k222s', 'k444s']:
    pathlib.Path(sys.argv[2], case+'.toml').write_text(bench.case_input(sys.argv[1], case, 1))
PY
for case in gamma k222 k222s k444s; do
    taskset -c "$P10B_CPUS" target/p10b-mkl/release/dft-cli run \
        "$P10B_OUT/$case.toml" --output-dir "$P10B_OUT/qe-$case"
done
python3 -I scripts/verify_qe.py --ours "$P10B_OUT/qe-gamma/results.json" \
    --bundle tests/reference/si30_ours --json "$P10B_OUT/qe-gamma.json" --markdown "$P10B_OUT/qe-gamma.md"
for case in k222 k222s; do
    python3 -I scripts/compare_kpoints_qe.py --ours "$P10B_OUT/qe-$case/results.json" \
        --qe-out "tests/reference/kpoints/${case}_c/qe.out" | tee "$P10B_OUT/qe-$case.log"
done
python3 -I scripts/compare_kpoints_qe.py --ours "$P10B_OUT/qe-k444s/results.json" \
    --qe-summary tests/reference/kpoints/k444s_c/qe_summary.json | tee "$P10B_OUT/qe-k444s.log"
# bench_scf enforces the energy tolerance and density hash, but not energy bit
# identity across thread counts; check that extra requirement explicitly.
python3 -I - "$P10B_OUT" <<'PY'
import json, pathlib, struct, sys
root = pathlib.Path(sys.argv[1])
base = json.loads((root/'scf-mkl.json').read_text())['cases']
for threads in [2,4,8]:
    other = json.loads((root/f'scf-mkl-t{threads}.json').read_text())['cases']
    for case in base:
        assert struct.pack('d',base[case]['total_energy_ha']) == struct.pack('d',other[case]['total_energy_ha']), (case,threads)
print('PASS MKL SCF energy bits across 1/2/4/8 threads')
PY
printf '%s\n' 'PASS P10B dev08 protocol; retain all evidence, review generality before recommendation'
