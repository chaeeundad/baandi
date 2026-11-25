"""Synthetic G4G Gaussian-projector atom; prints UPF, never uses real UPF data.

Standard-library generation keeps this CLI fixture independent of the engine.
Units: local potential and projector coefficients in Ry, radial mesh in Bohr.
"""
import math


def array(values):
    return " ".join(format(value, ".17e") for value in values)


def potential():
    n = 401
    z = 6.0
    dx = math.log(601.0) / (n - 1)
    r = [0.02 * math.expm1(i * dx) for i in range(n)]
    local = [
        (-4 * z / math.sqrt(math.pi) if x == 0 else -2 * z * math.erf(x) / x)
        - 20 * math.exp(-x * x)
        for x in r
    ]
    beta = "".join(
        f'<PP_BETA.{i} index="{i}" angular_momentum="{l}" '
        f'cutoff_radius="12" cutoff_radius_index="401">'
        f'{array(x ** (l + 1) * math.exp(-x * x / i) for x in r)}</PP_BETA.{i}>'
        for i, l in enumerate([0, 0, 1], 1)
    )
    return f'''<UPF version="2.0.1"><PP_INFO>SYNTHETIC G4G</PP_INFO><PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="O" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="PBE" z_valence="6" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="1" l_max_rho="0" l_local="-1" mesh_size="401" number_of_wfc="2" number_of_proj="3"/><PP_MESH mesh="401"><PP_R>{array(r)}</PP_R><PP_RAB>{array(dx * (x + 0.02) for x in r)}</PP_RAB></PP_MESH><PP_LOCAL>{array(local)}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>2 0.4 0 0.4 4 0 0 0 6</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{array(z * 4 / math.sqrt(math.pi) * x * x * math.exp(-x * x) for x in r)}</PP_RHOATOM><PP_PSWFC><PP_CHI.1 l="0" label="s" occupation="1">{array(x * math.exp(-0.6 * x * x) for x in r)}</PP_CHI.1><PP_CHI.2 l="1" label="p" occupation="0">{array(x * x * math.exp(-0.5 * x * x) for x in r)}</PP_CHI.2></PP_PSWFC></UPF>'''


if __name__ == "__main__":
    print(potential())
