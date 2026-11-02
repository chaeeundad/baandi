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
    fixture_with_points(601)
}

/// Generate synthetic UPF text at a selected radial resolution.
pub fn fixture_with_points(n: usize) -> String {
    fixture_with_options(n, &Options::default())
}

pub struct Options<'a> {
    pub relativistic: &'a str,
    pub wrap: bool,
    pub separator: &'a str,
    pub crlf: bool,
    pub wavefunctions: bool,
}
impl Default for Options<'_> {
    fn default() -> Self {
        Self {
            relativistic: "scalar",
            wrap: false,
            separator: " ",
            crlf: false,
            wavefunctions: false,
        }
    }
}

pub fn fixture_with_options(n: usize, options: &Options<'_>) -> String {
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
    let mut text = format!(
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
    );
    text = text.replace(
        "relativistic=\"scalar\"",
        &format!("relativistic=\"{}\"", options.relativistic),
    );
    if options.wavefunctions {
        text = text.replace("number_of_wfc=\"0\"", "number_of_wfc=\"3\"");
        let mut chi = String::from("<PP_PSWFC>");
        for i in 1..=3 {
            chi.push_str(&format!(
                "<PP_CHI.{i} label=\"synthetic{i}\" l=\"0\" occupation=\"1\">{}</PP_CHI.{i}>",
                array(&r)
            ));
        }
        chi.push_str("</PP_PSWFC></UPF>");
        text = text.replace("</UPF>", &chi);
    }
    for tag in [
        "PP_R",
        "PP_RAB",
        "PP_LOCAL",
        "PP_BETA.1",
        "PP_BETA.2",
        "PP_BETA.3",
        "PP_DIJ",
        "PP_CHI.1",
        "PP_CHI.2",
        "PP_CHI.3",
        "PP_RHOATOM",
    ] {
        let Some(start) = text
            .find(&format!("<{tag}>"))
            .or_else(|| text.find(&format!("<{tag} ")))
        else {
            continue;
        };
        let body = start + text[start..].find('>').unwrap() + 1;
        let end = text.find(&format!("</{tag}>")).unwrap();
        let tokens: Vec<_> = text[body..end].split_ascii_whitespace().collect();
        let mut formatted = String::new();
        for (i, token) in tokens.iter().enumerate() {
            if i > 0 {
                formatted.push_str(if options.wrap && i % 4 == 0 {
                    "\n"
                } else {
                    options.separator
                });
            }
            formatted.push_str(token);
        }
        text.replace_range(body..end, &formatted);
        if options.wrap {
            text.insert_str(body - 1, " columns=\"4\"");
        }
    }
    if options.crlf {
        text = text.replace('\n', "\r\n");
    }
    text
}
