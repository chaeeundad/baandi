#!/usr/bin/env python3
"""Generate NON-PHYSICAL NC-format UPF v2 fixtures; standard library only.

No atomic reference solution or norm-conservation/transferability claim is made.
The NC header exercises the supported parser/operator contract, as in crate tests.
Radial coordinates and beta=r*projector use Bohr; local and DIJ use Ry.
"""
import argparse
import math
from pathlib import Path


def generate(valence=4, channels=(0, 0, 1), radius=1.0, shape="gaussian", points=1201):
    if isinstance(valence, bool) or not isinstance(valence, int) or not 1 <= valence <= 14:
        raise ValueError("valence must be an integer from 1 to 14")
    if not 1 <= len(channels) <= 3 or any(type(l) is not int or l not in (0, 1, 2) for l in channels):
        raise ValueError("one to three total projectors with l=0,1,2 are supported")
    if not math.isfinite(radius) or radius <= 0 or shape not in ("gaussian", "soft"):
        raise ValueError("positive finite radius and gaussian/soft shape required")
    if type(points) is not int or points < 101 or points % 2 != 1:
        raise ValueError("points must be odd and at least 101")
    dx = math.log1p(12.0 / 0.02) / (points - 1)
    r = [0.02 * math.expm1(i * dx) for i in range(points)]
    rab = [dx * (x + 0.02) for x in r]
    if shape == "gaussian":
        local = [-4 * valence / (radius * math.sqrt(math.pi)) if x == 0
                 else -2 * valence * math.erf(x / radius) / x for x in r]
    else:
        local = [-2 * valence / math.sqrt(x*x + radius*radius) for x in r]
    def array(xs):
        return " ".join(format(x, ".17e") for x in xs)
    beta = []
    for i, l in enumerate(channels):
        values = [x**(l+1) * math.exp(-x*x/(i+1)) for x in r]
        beta.append(f'<PP_BETA.{i+1} index="{i+1}" angular_momentum="{l}" cutoff_radius="12" cutoff_radius_index="{points}">{array(values)}</PP_BETA.{i+1}>')
    dij = [2.0*(i+1) if i == j else (0.4 if channels[i] == channels[j] else 0.0)
           for i in range(len(channels)) for j in range(len(channels))]
    return f'''<UPF version="2.0.1"><PP_INFO>SYNTHETIC NON-PHYSICAL benchmark fixture; not an atomic pseudopotential.</PP_INFO>
<PP_HEADER generated="synthetic" author="rust-dft benchmark" date="2026" comment="NON-PHYSICAL" element="Xx"
pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F"
has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ"
z_valence="{valence}" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="{max(channels)}" l_max_rho="0"
l_local="-1" mesh_size="{points}" number_of_wfc="0" number_of_proj="{len(channels)}"/>
<PP_MESH dx="{dx:.17e}" mesh="{points}" xmin="0" rmax="12" zmesh="1"><PP_R>{array(r)}</PP_R><PP_RAB>{array(rab)}</PP_RAB></PP_MESH>
<PP_LOCAL>{array(local)}</PP_LOCAL><PP_NONLOCAL>{''.join(beta)}<PP_DIJ>{array(dij)}</PP_DIJ></PP_NONLOCAL>
<PP_RHOATOM>{array([0.0]*points)}</PP_RHOATOM></UPF>
'''


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("output", type=Path)
    p.add_argument("--valence", type=int, default=4)
    p.add_argument("--channels", default="0,0,1", help="angular momenta, e.g. 0,1,2 or 2,2,2")
    p.add_argument("--radius", type=float, default=1.0)
    p.add_argument("--shape", choices=["gaussian", "soft"], default="gaussian")
    p.add_argument("--points", type=int, default=1201)
    a = p.parse_args()
    try:
        data = generate(a.valence, tuple(int(l) for l in a.channels.split(',')), a.radius, a.shape, a.points)
    except ValueError as e:
        p.error(str(e))
    a.output.write_text(data, encoding="utf-8")


if __name__ == "__main__":
    main()
