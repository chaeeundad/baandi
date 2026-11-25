//! Purely synthetic in-memory UPF data. Never a real-element reference.
use dft_geometry::{Atom, Cell};
use dft_physics::data::erf;
use dft_upf::{Pseudopotential, parse_upf};
use std::f64::consts::PI;
pub fn cell() -> Cell {
    let a = 10.261_212_856_717_933;
    Cell::from_user_rows_bohr([
        [0.0, a / 2.0, a / 2.0],
        [a / 2.0, 0.0, a / 2.0],
        [a / 2.0, a / 2.0, 0.0],
    ])
    .unwrap()
}
pub fn atoms() -> Vec<Atom> {
    vec![Atom::new("Xx", [0.0; 3]), Atom::new("Xx", [0.25; 3])]
}
fn array(v: &[f64]) -> String {
    v.iter()
        .map(|v| format!("{v:.17e}"))
        .collect::<Vec<_>>()
        .join(" ")
}
pub fn potential(z: f64, rc: f64, nonlocal: bool, points: usize) -> Pseudopotential {
    potential_with_functional(z, rc, nonlocal, points, "SLA PZ")
}
pub fn potential_with_functional(
    z: f64,
    rc: f64,
    nonlocal: bool,
    points: usize,
    functional: &str,
) -> Pseudopotential {
    build(z, rc, nonlocal, points, functional, false)
}
/// As `potential`, with a Gaussian atomic charge in PP_RHOATOM (4 pi r^2 rho(r),
/// width rc, integrating to z) so that an atomic start density exists.
#[allow(dead_code)]
pub fn potential_with_atomic_charge(
    z: f64,
    rc: f64,
    nonlocal: bool,
    points: usize,
) -> Pseudopotential {
    build(z, rc, nonlocal, points, "SLA PZ", true)
}
/// As `potential_with_functional`, with synthetic s and p atomic
/// orbitals in PP_PSWFC (r^(l+1) exp(-a r^2), not normalized; the atomic start
/// orthonormalizes its span) so that the atomic start applies (ONEPATH-3B).
#[allow(dead_code)]
pub fn potential_with_orbitals(
    z: f64,
    rc: f64,
    nonlocal: bool,
    points: usize,
    functional: &str,
) -> Pseudopotential {
    let text = build_text(z, rc, nonlocal, points, functional, false);
    let pp = parse_upf(&text).unwrap();
    let chi = |l: i32, a: f64| {
        array(
            &pp.mesh()
                .r()
                .iter()
                .map(|r| r.powi(l + 1) * (-a * r * r).exp())
                .collect::<Vec<_>>(),
        )
    };
    let text = text
        .replace("number_of_wfc=\"0\"", "number_of_wfc=\"2\"")
        .replace(
            "</UPF>",
            &format!(
                "<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1><PP_CHI.2 l=\"1\" label=\"p\" occupation=\"2\">{}</PP_CHI.2></PP_PSWFC></UPF>",
                chi(0, 0.6),
                chi(1, 0.5)
            ),
        );
    parse_upf(&text).unwrap()
}
/// Synthetic PP with both charge and orbitals for the reference driver.
#[allow(dead_code)]
pub fn potential_with_atomic_start(
    z: f64,
    rc: f64,
    nonlocal: bool,
    points: usize,
    functional: &str,
) -> Pseudopotential {
    let text = build_text(z, rc, nonlocal, points, functional, true);
    let pp = parse_upf(&text).unwrap();
    let chi = |l: i32, a: f64| {
        array(
            &pp.mesh()
                .r()
                .iter()
                .map(|r| r.powi(l + 1) * (-a * r * r).exp())
                .collect::<Vec<_>>(),
        )
    };
    let text = text
        .replace("number_of_wfc=\"0\"", "number_of_wfc=\"2\"")
        .replace(
            "</UPF>",
            &format!(
                "<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1><PP_CHI.2 l=\"1\" label=\"p\" occupation=\"2\">{}</PP_CHI.2></PP_PSWFC></UPF>",
                chi(0, 0.6),
                chi(1, 0.5)
            ),
        );
    parse_upf(&text).unwrap()
}
fn build(
    z: f64,
    rc: f64,
    nonlocal: bool,
    points: usize,
    functional: &str,
    atomic_charge: bool,
) -> Pseudopotential {
    parse_upf(&build_text(
        z,
        rc,
        nonlocal,
        points,
        functional,
        atomic_charge,
    ))
    .unwrap()
}
fn build_text(
    z: f64,
    rc: f64,
    nonlocal: bool,
    points: usize,
    functional: &str,
    atomic_charge: bool,
) -> String {
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
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="{functional}" z_valence="{z}" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="1" l_max_rho="0" l_local="-1" mesh_size="{points}" number_of_wfc="0" number_of_proj="{np}" />
<PP_MESH dx="{dx}" mesh="{points}" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH><PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>{dij}</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,
        array(&r),
        array(&rab),
        array(&local),
        array(&if atomic_charge {
            let alpha = 1.0 / (rc * rc);
            r.iter()
                .map(|r| z * 4.0 * PI * r * r * (alpha / PI).powf(1.5) * (-alpha * r * r).exp())
                .collect()
        } else {
            vec![0.0; points]
        })
    );
    if nonlocal {
        xml
    } else {
        xml.replace("<PP_DIJ></PP_DIJ>", "")
    }
}

/// G1 synthetic NC dataset: analytic Gaussian core with unit charge and width
/// parameter a=0.5 Bohr^-2, plus the existing atomic charge and s/p span.
#[allow(dead_code)]
pub fn potential_with_nlcc(nonlocal: bool, functional: &str) -> Pseudopotential {
    let text = build_text(2., 1.8, nonlocal, 401, functional, true);
    let pp = parse_upf(&text).unwrap();
    let samples =
        |f: &dyn Fn(f64) -> f64| array(&pp.mesh().r().iter().map(|&r| f(r)).collect::<Vec<_>>());
    let text = text
        .replace("core_correction=\"F\"", "core_correction=\"T\"")
        .replace("number_of_wfc=\"0\"", "number_of_wfc=\"2\"")
        .replace("</UPF>", &format!(
            "<PP_NLCC>{}</PP_NLCC><PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1><PP_CHI.2 l=\"1\" label=\"p\" occupation=\"0\">{}</PP_CHI.2></PP_PSWFC></UPF>",
            samples(&|r| (0.5 / PI).powf(1.5) * (-0.5 * r * r).exp()),
            samples(&|r| r * (-0.6 * r * r).exp()),
            samples(&|r| r * r * (-0.5 * r * r).exp()),
        ));
    parse_upf(&text).unwrap()
}
