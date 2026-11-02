//! Synthetic fixtures only; no real element or pseudopotential data.
pub fn mesh(n: usize, end: f64) -> (Vec<f64>, Vec<f64>) {
    let r0 = 0.02;
    let dx = (1.0 + end / r0).ln() / (n - 1) as f64;
    let r: Vec<_> = (0..n).map(|i| r0 * (i as f64 * dx).exp_m1()).collect();
    let rab = r.iter().map(|r| dx * (r + r0)).collect();
    (r, rab)
}
pub fn array(x: &[f64]) -> String {
    x.iter()
        .map(|v| format!("{v:.17e}"))
        .collect::<Vec<_>>()
        .join(" ")
}
pub fn fixture() -> String {
    let n = 601;
    let (r, rab) = mesh(n, 12.0);
    let local: Vec<_> = r
        .iter()
        .map(|r| {
            if *r == 0.0 {
                -16.0 / std::f64::consts::PI.sqrt()
            } else {
                -8.0 * libm::erf(*r) / r
            }
        })
        .collect();
    let mut beta = String::new();
    for (i, l) in [0, 0, 1].into_iter().enumerate() {
        let values: Vec<_> = r
            .iter()
            .map(|r| r.powi(l + 1) * (-r * r / (i + 1) as f64).exp())
            .collect();
        beta.push_str(&format!("<PP_BETA.{} index=\"{}\" angular_momentum=\"{l}\" cutoff_radius=\"12\" cutoff_radius_index=\"{n}\">{}</PP_BETA.{}>",i+1,i+1,array(&values),i+1));
    }
    format!(
        r#"<UPF version="2.0.1">
<PP_INFO>SYNTHETIC Xx fixture; not a real pseudopotential.</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic"
 element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F"
 has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ"
 z_valence="4" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="1" l_max_rho="0" l_local="-1"
 mesh_size="{n}" number_of_wfc="0" number_of_proj="3" />
<PP_MESH dx="0.01" mesh="{n}" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH>
<PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>2 0.4 0 0.4 4 0 0 0 6</PP_DIJ></PP_NONLOCAL>
<PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,
        array(&r),
        array(&rab),
        array(&local),
        array(&vec![0.0; n])
    )
}
