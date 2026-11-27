#!/bin/bash
# R2 oracle comparison, one case (written by the dev08 session 2026-11-27, MIT; paths parameterized by the orchestrator).
# usage: r2_case.sh <case>   case: si2-444 si2-888 c2-444 lif-444 aln-444 si2-skewed-444
# BAANDI `dielectric` -> QE 7.6 pw.x scf + nscf (make_qe_r2.py: full Gamma mesh, nosym/noinv, nbnd = N_b + 8, same nr/ecut/UPF) ->
# Yambo 5.3.0 p2y + setup + optics (Chimod HARTREE gives eps with LF and eps_o without LF) along x, y, z at omega = 0,
# NGsBlkXd = BAANDI response_size in RL, BndsRnXd = DipBands = 1..N_b, QpntsRXd 1..1, damping 1e-4 eV (yambo_lf_inputs.py).
# Accuracy only (decision D10). Environment (defaults are the dev08 layout of 2026-11-27):
#   R2_WORK  work root (case dir is $R2_WORK/<case>)     R2_SRC  directory with the case tomls (this directory)
#   R2_INPUTS optional directory with variant tomls       BAANDI  baandi binary (a BAANDI.sha256 file next to it is read if present)
#   PW  pw.x    YAMBO_BIN  directory with yambo, p2y      PY311  Python >= 3.11 (tomllib) for make_qe_r2.py
#   PSEUDO_DIRS  space-separated directories searched for the UPF files     CORES8 / CORES4  core lists for binding
set -u
CASE=$1
HERE=$(cd "$(dirname "$0")" && pwd)
R2_WORK=${R2_WORK:-/scr/mino/rustdft/r2-yambo}; SRC=${R2_SRC:-$HERE}; R2_INPUTS=${R2_INPUTS:-$R2_WORK/inputs}
B=${BAANDI:?set BAANDI to the baandi binary}; PW=${PW:-/scr/mino/oracles/qe-7.6/bin/pw.x}; Y=${YAMBO_BIN:-/scr/mino/oracles/yambo/bin}
PY=${PY311:-python3}; PSEUDO_DIRS=${PSEUDO_DIRS:-"$HOME/fq/pseudo /scr/mino/rustdft/pseudo"}
CORES8=${CORES8:-16,17,18,19,36,37,38,39}; CORES4=${CORES4:-16,17,18,19}
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
W=$R2_WORK/$CASE
[ -e $W ] && { echo "ABORT: $W exists"; exit 1; }; mkdir -p $W; cd $W
[ -f $R2_INPUTS/$CASE.toml ] && SRC=$R2_INPUTS
cp $SRC/$CASE.toml in.toml
for u in $(python3 -c 'import re,sys; s=open("in.toml").read(); sec=s.split("[pseudopotentials]")[1].split("[")[0]; print(" ".join(re.findall(r"=\s*\"([^\"]+)\"",sec)))'); do
  for d in $PSEUDO_DIRS; do [ -f $d/$u ] && { cp $d/$u .; break; }; done
  [ -f $u ] || { echo "missing UPF $u"; exit 1; }; done
{ echo "r2 $CASE host $(hostname) start $(date -Is)"; echo "script $0 sha256 $(sha256sum $0 | cut -c1-64)"; echo "baandi $B sha256 $(sha256sum $B | cut -c1-64)";
  echo "input $SRC/$CASE.toml sha256 $(sha256sum $SRC/$CASE.toml | cut -c1-64)"; for u in *.upf *.UPF; do [ -f $u ] && echo "UPF $u sha256 $(sha256sum $u | cut -c1-64)"; done
  echo "pw.x $PW sha256 $(sha256sum $(readlink -f $PW) | cut -c1-64)"; echo "yambo $Y/yambo sha256 $(sha256sum $Y/yambo | cut -c1-64); p2y sha256 $(sha256sum $Y/p2y | cut -c1-64)"; } > MANIFEST.txt
tm() { local n=$1; shift; local t0=$(date +%s%N); "$@"; local rc=$?; echo "$n rc=$rc wall_ms=$(( ($(date +%s%N)-t0)/1000000 ))" | tee -a MANIFEST.txt; return $rc; }
NR=$($B grid in.toml 2>&1 | grep -o 'FFT dims[^0-9]*[0-9]* x [0-9]* x [0-9]*' | grep -o '[0-9]*' | tr '\n' ' '); echo "baandi grid: $NR" >> MANIFEST.txt
tm baandi_dielectric taskset -c $CORES8 env OMP_NUM_THREADS=8 RAYON_NUM_THREADS=8 $B dielectric in.toml --output-dir baandi > baandi.out 2> baandi.err || { echo "BAANDI failed"; tail -5 baandi.err | tee -a MANIFEST.txt; exit 2; }
NG=$(python3 -c 'import json; print(json.load(open("baandi/epsilon.json"))["settings"]["response_size"])'); NB=$(python3 -c 'import json; print(json.load(open("baandi/epsilon.json"))["settings"]["bands"])')
echo "BAANDI response_size N_G $NG bands $NB" >> MANIFEST.txt
$PY -I $HERE/make_qe_r2.py in.toml --nr $NR --prefix qe --dir . | tee -a MANIFEST.txt
MP="mpirun -np 8 --map-by pe-list=$CORES8:ordered --bind-to core"
tm qe_scf $MP $PW -in qe-scf.in > qe-scf.out 2> qe-scf.err || exit 3
tm qe_nscf $MP $PW -nk 8 -in qe-nscf.in > qe-nscf.out 2> qe-nscf.err || exit 3
echo "qe nscf not-converged lines: $(grep -c 'not converged' qe-nscf.out)" >> MANIFEST.txt; grep -m1 'Program PWSCF' qe-scf.out >> MANIFEST.txt
( cd qe-work/qe.save && $Y/p2y > ../../p2y.out 2>&1 ); echo "p2y rc=$?" >> MANIFEST.txt; mkdir -p yambo; mv qe-work/qe.save/SAVE yambo/; cd yambo
$Y/yambo > setup.out 2>&1; $Y/yambo -o c -k hartree -V all -F gen.in -Q > gen.out 2>&1
python3 -I $HERE/yambo_lf_inputs.py $NG $NB
for d in x y z; do tm yambo_$d mpirun -np 4 --map-by pe-list=$CORES4:ordered --bind-to core $Y/yambo -F lf_$d.in -J lf$d > run_$d.out 2>&1; grep -h 'Hartree size' o-lf$d.eps_q1_inv_rpa_dyson >> ../MANIFEST.txt; done
grep -h -m1 'Version' o-lfx.eps_q1_inv_rpa_dyson >> ../MANIFEST.txt; grep -h 'break level degeneracy' r-lf*_optics* | sort -u >> ../MANIFEST.txt; cd ..
python3 - <<'PYEOF' | tee compare.txt
import json
d = json.load(open("baandi/epsilon.json")); e = d["epsilon_infinity"]; out = {"case": d["title"], "N_G": d["settings"]["response_size"], "bands": d["settings"]["bands"], "dirs": {}}
for i, a in enumerate("xyz"):
    row = [l.split() for l in open(f"yambo/o-lf{a}.eps_q1_inv_rpa_dyson") if not l.startswith("#") and l.strip()][0]
    w, lf, nolf = float(row[0]), float(row[2]), float(row[4])
    bl, bn = e["with_local_fields"][i][i], e["without_local_fields"][i][i]
    out["dirs"][a] = dict(omega_eV=w, yambo_LF=lf, baandi_LF=bl, rel_LF=(bl - lf) / lf, yambo_noLF=nolf, baandi_noLF=bn, rel_noLF=(bn - nolf) / nolf)
mx = max(max(abs(v["rel_LF"]), abs(v["rel_noLF"])) for v in out["dirs"].values()); out["max_rel"] = mx; out["gate_0p5pct"] = "PASS" if mx <= 5e-3 else "FAIL"
json.dump(out, open("r2-compare.json", "w"), indent=1); print(json.dumps(out, indent=1))
PYEOF
echo "end $(date -Is)" >> MANIFEST.txt; echo R2_CASE_DONE $CASE
