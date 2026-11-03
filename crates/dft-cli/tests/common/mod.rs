//! Synthetic soft-Coulomb fixture, generated only in temporary test directories.
pub fn fixture() -> String {
    let n = 401;
    let r = (0..n)
        .map(|i| format!("{}", i as f64 * 0.025))
        .collect::<Vec<_>>()
        .join(" ");
    let rab = vec!["0.025"; n].join(" ");
    // Smooth Coulomb-tailed synthetic potential, in Ry.
    let local = (0..n)
        .map(|i| format!("{}", -8.0 / ((i as f64 * 0.025).powi(2) + 1.0).sqrt()))
        .collect::<Vec<_>>()
        .join(" ");
    let beta = (0..n)
        .map(|i| {
            let r = i as f64 * 0.025;
            format!("{}", r * (-r * r).exp())
        })
        .collect::<Vec<_>>()
        .join(" ");
    let zero = vec!["0"; n].join(" ");
    format!(
        r#"<UPF version="2.0.1"><PP_INFO>SYNTHETIC CLI TEST ONLY</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx"
pseudo_type="NC" relativistic="nonrelativistic" is_ultrasoft="F" is_paw="F" is_coulomb="F"
has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ"
z_valence="4" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="0" l_max_rho="0"
l_local="-1" mesh_size="{n}" number_of_wfc="0" number_of_proj="1"/>
<PP_MESH dx="0.025" mesh="{n}" xmin="0" rmax="10" zmesh="1"><PP_R>{r}</PP_R><PP_RAB>{rab}</PP_RAB></PP_MESH>
<PP_LOCAL>{local}</PP_LOCAL><PP_NONLOCAL><PP_BETA.1 index="1" angular_momentum="0" cutoff_radius="10" cutoff_radius_index="401">{beta}</PP_BETA.1><PP_DIJ>2</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{zero}</PP_RHOATOM></UPF>"#
    )
}

use dft_physics::data::erf;
use std::f64::consts::PI;
fn array(v: &[f64]) -> String {
    v.iter()
        .map(|v| format!("{v:.17e}"))
        .collect::<Vec<_>>()
        .join(" ")
}
pub fn gaussian() -> String {
    gaussian_with(4.0, 1.0, true, 1201)
}

#[allow(dead_code)]
pub fn gaussian_with(z: f64, rc: f64, nonlocal: bool, points: usize) -> String {
    let end = 12.0;
    let r0 = 0.02;
    let dx = (1.0_f64 + end / r0).ln() / (points - 1) as f64;
    let r: Vec<_> = (0..points).map(|i| r0 * (i as f64 * dx).exp_m1()).collect();
    let rab: Vec<_> = r.iter().map(|r| dx * (r + r0)).collect();
    let local: Vec<_> = r
        .iter()
        .map(|r| {
            if *r == 0.0 {
                -4.0 * z / (rc * PI.sqrt())
            } else {
                -2.0 * z * erf(r / rc) / r
            }
        })
        .collect();
    let mut beta = String::new();
    let (np, dij) = if nonlocal {
        for (i, l) in [0, 0, 1].into_iter().enumerate() {
            let values: Vec<_> = r
                .iter()
                .map(|r| r.powi(l + 1) * (-r * r / (i + 1) as f64).exp())
                .collect();
            beta.push_str(&format!("<PP_BETA.{} index=\"{}\" angular_momentum=\"{l}\" cutoff_radius=\"12\" cutoff_radius_index=\"{points}\">{}</PP_BETA.{}>",i+1,i+1,array(&values),i+1));
        }
        (3, "2 0.4 0 0.4 4 0 0 0 6")
    } else {
        (0, "")
    };
    let xml = format!(
        r#"<UPF version="2.0.1"><PP_INFO>SYNTHETIC Gaussian ion, no real data.</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ" z_valence="{z}" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="1" l_max_rho="0" l_local="-1" mesh_size="{points}" number_of_wfc="0" number_of_proj="{np}" />
<PP_MESH dx="{dx}" mesh="{points}" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH><PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>{dij}</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,
        array(&r),
        array(&rab),
        array(&local),
        array(&vec![0.0; points])
    );
    if nonlocal {
        xml
    } else {
        xml.replace("<PP_DIJ></PP_DIJ>", "")
    }
}
