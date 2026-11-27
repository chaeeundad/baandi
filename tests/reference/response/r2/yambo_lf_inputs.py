#!/usr/bin/env python3
"""R2 oracle: turn the Yambo-generated gen.in (`yambo -o c -k hartree -V all -F gen.in -Q`) into lf_x/y/z.in matching BAANDI.
NG = BAANDI settings.response_size (RL units), NB = settings.bands. Static point omega = 0, damping 1e-4 eV, q = 1 only.
Then per direction: `mpirun -np 4 yambo -F lf_<d>.in -J lf<d>`; read the omega = 0 row of o-lf<d>.eps_q1_inv_rpa_dyson:
column 3 = Re eps (with local fields), column 5 = Re eps_o (without).
Written by the dev08 session on 2026-11-27 (MIT, no Yambo source).
usage: yambo_lf_inputs.py NG NB   (in the directory holding gen.in)"""
import re, sys
ng, nb = sys.argv[1], sys.argv[2]; s = open("gen.in").read()
def sub(pat, rep):
    global s; n = len(re.findall(pat, s, flags=re.M)); assert n == 1, (pat, n); s = re.sub(pat, rep, s, flags=re.M)
sub(r"^NGsBlkXd=.*$", f"NGsBlkXd= {ng} RL    # = BAANDI response_size")
sub(r"^(\s*)1 \|\s*\d+ \|(\s*# \[Xd\] Polarization function bands)", rf"   1 | {nb} |\2")
sub(r"^(\s*)1 \|\s*\d+ \|(\s*# \[DIP\] Bands range for dipoles)", rf"   1 | {nb} |\2")
sub(r"^(\s*)1 \|\s*\d+ \|(\s*# \[Xd\] Transferred momenta)", r"   1 |  1 |\2")
sub(r"^\s*[0-9.]+ \|\s*[0-9.]+ \|(\s*eV\s*# \[Xd\] Energy range)", r"  0.00000 | 0.01000 |\1")
sub(r"^\s*[0-9.]+ \|\s*[0-9.]+ \|(\s*eV\s*# \[Xd\] Damping range)", r" 0.000100 | 0.000100 |\1")
sub(r"^ETStpsXd=.*$", "ETStpsXd= 2")
for d, v in (("x", "1.0 | 0.0 | 0.0"), ("y", "0.0 | 1.0 | 0.0"), ("z", "0.0 | 0.0 | 1.0")):
    t = re.sub(r"^\s*[-0-9.]+ \|\s*[-0-9.]+ \|\s*[-0-9.]+ \|(\s*# \[Xd\] \[cc\] Electric Field)", f" {v} |\\1", s, flags=re.M)
    open(f"lf_{d}.in", "w").write(t)
