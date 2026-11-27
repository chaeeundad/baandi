//! R2-VEL-CHECK: the full velocity matrix `<psi_i| dH_k/dk_b |psi_j>`
//! (kinetic + nonlocal, `dft_response::velocity::velocity`, as used by the
//! dielectric head and wings) against a central finite difference of the
//! Hamiltonian at k +- h e_b with the SAME plane-wave set (Miller indices
//! fixed, only k + G changes). The finite-difference side never calls
//! `projector_k_gradient`: its nonlocal part is `sum_atom R D R^H` with the
//! rows R of the production `NonlocalOperator::new` (the operator the SCF and
//! NSCF Hamiltonians are built from) at the shifted k, so it checks the
//! projector gradients, the couplings and the assembly of the velocity
//! independently. The local potential does not depend on k and drops out of
//! the difference; it only shapes the test vectors.
//!
//! Fast gate (debug, synthetic multi-projector s/p/d/f dataset with coupled D,
//! skewed cell, two atoms): `cargo test -p dft-response --test velocity_fd -- --nocapture`
//! Real PseudoDojo datasets (LiF and wurtzite AlN of the r2 cases, 30 Ha, release):
//! `DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-response --test velocity_fd -- --ignored --nocapture`
mod common;
use common::array;
use dft_basis::PlaneWaveBasis;
use dft_geometry::Cell;
use dft_numerics::linalg::{ComplexMatrix, hermitian_eigh};
use dft_physics::hamiltonian::{NonlocalOperator, Species};
use dft_physics::stress::projector_k_gradient;
use dft_response::{Complex64 as C, velocity::velocity};
use dft_units::Hartree;
use std::collections::HashMap;
use std::f64::consts::PI;

/// Synthetic NC dataset (not a real element): two projectors per channel for
/// l = 0, 1, 2 (the second with a radial node) and one f projector, with
/// off-diagonal D inside each l block (the multi-projector structure of
/// PseudoDojo/ONCVPSP files, here also with the coupling they may carry).
fn multi_projector_upf() -> dft_upf::Pseudopotential {
    let points = 601;
    let end = 12.0;
    let r0 = 0.01;
    let dx = (1.0_f64 + end / r0).ln() / (points - 1) as f64;
    let r: Vec<_> = (0..points).map(|i| r0 * (i as f64 * dx).exp_m1()).collect();
    let rab: Vec<_> = r.iter().map(|r| dx * (r + r0)).collect();
    let mut beta = String::new();
    let shapes: [(i32, f64, bool); 7] = [
        (0, 1.0, false),
        (0, 1.6, true),
        (1, 1.2, false),
        (1, 1.9, true),
        (2, 1.4, false),
        (2, 2.2, true),
        (3, 1.7, false),
    ];
    for (i, &(l, width, node)) in shapes.iter().enumerate() {
        let values: Vec<_> = r
            .iter()
            .map(|r| {
                let radial = if node {
                    1.0 - r * r / (1.5 * width)
                } else {
                    1.0
                };
                r.powi(l + 1) * radial * (-r * r / width).exp()
            })
            .collect();
        beta.push_str(&format!(
            "<PP_BETA.{n} index=\"{n}\" angular_momentum=\"{l}\" cutoff_radius=\"12\" cutoff_radius_index=\"{points}\">{}</PP_BETA.{n}>",
            array(&values),
            n = i + 1
        ));
    }
    let mut d = [[0.0_f64; 7]; 7];
    for (block, values) in [
        (0, [0.9, 0.3, -0.4]),
        (2, [-0.7, 0.2, 0.6]),
        (4, [0.5, -0.25, 0.35]),
    ] {
        d[block][block] = values[0];
        d[block][block + 1] = values[1];
        d[block + 1][block] = values[1];
        d[block + 1][block + 1] = values[2];
    }
    d[6][6] = 0.2;
    let dij: Vec<f64> = d.iter().flatten().copied().collect();
    let text = format!(
        r#"<UPF version="2.0.1"><PP_INFO>SYNTHETIC, no real data.</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ" z_valence="2" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="3" l_max_rho="0" l_local="-1" mesh_size="{points}" number_of_wfc="0" number_of_proj="7" />
<PP_MESH dx="{dx}" mesh="{points}" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH><PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>{}</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,
        array(&r),
        array(&rab),
        array(&vec![0.0; points]),
        array(&dij),
        array(&vec![0.0; points]),
    );
    dft_upf::parse_upf(&text).unwrap()
}

/// `<x_i| V_NL |y_j>` from the production operator restricted to `rows`.
fn nonlocal_matrix(
    op: &NonlocalOperator,
    rows: &[usize],
    x: &ComplexMatrix,
) -> Result<ComplexMatrix, String> {
    let n = rows.len();
    let p = op.projector_count();
    let table = op.projector_rows(rows).map_err(|e| format!("{e:?}"))?;
    let bands = x.cols();
    // projections[band * p + c] = sum_g conj(R_c(g)) x_band(g)
    let mut projections = vec![C::default(); bands * p];
    for band in 0..bands {
        let column = x.column(band);
        for c in 0..p {
            let r = &table[c * n..(c + 1) * n];
            projections[band * p + c] = r.iter().zip(column).map(|(a, b)| a.conj() * b).sum();
        }
    }
    let mut coupled = vec![C::default(); bands * p];
    op.couple_projections(&projections, &mut coupled)
        .map_err(|e| format!("{e:?}"))?;
    let mut out = ComplexMatrix::zeros(bands, bands).map_err(|e| format!("{e:?}"))?;
    for i in 0..bands {
        for j in 0..bands {
            out[(i, j)] = (0..p)
                .map(|c| projections[i * p + c].conj() * coupled[j * p + c])
                .sum();
        }
    }
    Ok(out)
}

/// One crystal of the check: cell, species, cutoff and a local potential
/// (Fourier components by Miller difference; shapes the test vectors only).
struct Crystal<'a> {
    name: &'a str,
    cell: Cell,
    species: Vec<(&'a dft_upf::Pseudopotential, Vec<[f64; 3]>)>,
    cutoff: f64,
    local: HashMap<[i32; 3], C>,
    bands: usize,
}

impl Crystal<'_> {
    fn species(&self) -> Vec<Species<'_>> {
        self.species
            .iter()
            .map(|(p, r)| Species {
                potential: p,
                positions_bohr: r,
            })
            .collect()
    }
    /// Lowest `bands` eigenvectors of the dense H_k (production operator).
    fn vectors(&self, basis: &PlaneWaveBasis) -> ComplexMatrix {
        let n = basis.vectors().len();
        let species = self.species();
        let op = NonlocalOperator::new(basis, self.cell.volume_bohr3(), &species).unwrap();
        let all: Vec<usize> = (0..n).collect();
        let mut h = ComplexMatrix::zeros(n, n).unwrap();
        let identity = {
            let mut m = ComplexMatrix::zeros(n, n).unwrap();
            for g in 0..n {
                m[(g, g)] = C::new(1., 0.);
            }
            m
        };
        let v = nonlocal_matrix(&op, &all, &identity).unwrap();
        let miller: Vec<_> = basis.vectors().iter().map(|g| g.miller).collect();
        for j in 0..n {
            for i in 0..n {
                let d = std::array::from_fn(|a| miller[i][a] - miller[j][a]);
                h[(i, j)] = v[(i, j)] + self.local.get(&d).copied().unwrap_or_default();
            }
            h[(j, j)] += basis.momenta()[j].norm_squared / 2.;
        }
        for j in 0..n {
            h[(j, j)] = C::new(h[(j, j)].re, 0.);
            for i in 0..j {
                let a = (h[(i, j)] + h[(j, i)].conj()) * 0.5;
                h[(i, j)] = a;
                h[(j, i)] = a.conj();
            }
        }
        let (_, vectors) = hermitian_eigh(&h).unwrap();
        let mut out = ComplexMatrix::zeros(n, self.bands).unwrap();
        for b in 0..self.bands {
            out.column_mut(b).copy_from_slice(vectors.column(b));
        }
        out
    }
    /// `<x_i| H_{k'} |x_j>` without the local potential (k-independent) at
    /// the fractional k' (Cartesian shift `s` along `axis`), on the Miller
    /// set of `basis`: kinetic |k'+G|^2/2 and the production nonlocal rows.
    fn shifted(
        &self,
        basis: &PlaneWaveBasis,
        x: &ComplexMatrix,
        axis: usize,
        s: f64,
    ) -> ComplexMatrix {
        let a = self.cell.lattice_columns_bohr();
        let k = basis.k_fractional();
        let kf: [f64; 3] = std::array::from_fn(|j| k[j] + s * a[axis][j] / (2. * PI));
        // A larger sphere contains every Miller index of `basis` at k'.
        let wide = PlaneWaveBasis::at_k(&self.cell, Hartree(self.cutoff * 1.5 + 1.), kf).unwrap();
        assert!(
            (0..3).all(|j| (wide.k_fractional()[j] - kf[j]).abs() < 1e-12),
            "k' must not wrap (the Miller labels would shift)"
        );
        let index: HashMap<[i32; 3], usize> = wide
            .vectors()
            .iter()
            .enumerate()
            .map(|(i, g)| (g.miller, i))
            .collect();
        let rows: Vec<usize> = basis.vectors().iter().map(|g| index[&g.miller]).collect();
        let species = self.species();
        let op = NonlocalOperator::new(&wide, self.cell.volume_bohr3(), &species).unwrap();
        let mut m = nonlocal_matrix(&op, &rows, x).unwrap();
        for i in 0..x.cols() {
            for j in 0..x.cols() {
                let t: C = rows
                    .iter()
                    .enumerate()
                    .map(|(g, &r)| {
                        x[(g, i)].conj() * x[(g, j)] * (wide.momenta()[r].norm_squared / 2.)
                    })
                    .sum();
                m[(i, j)] += t;
            }
        }
        m
    }
}

fn max_abs(m: &ComplexMatrix) -> f64 {
    m.as_slice().iter().map(|z| z.norm()).fold(0., f64::max)
}
fn max_diff(a: &ComplexMatrix, b: &ComplexMatrix) -> f64 {
    a.as_slice()
        .iter()
        .zip(b.as_slice())
        .map(|(x, y)| (x - y).norm())
        .fold(0., f64::max)
}
fn diff(a: &ComplexMatrix, b: &ComplexMatrix) -> ComplexMatrix {
    let mut out = a.clone();
    for (z, y) in out.as_mut_slice().iter_mut().zip(b.as_slice()) {
        *z -= y;
    }
    out
}

struct Errors {
    /// Per h: max over k and b of max_ij |FD - p| / max_ij |p| (full velocity).
    full: Vec<f64>,
    /// Same for the nonlocal part (FD of V_NL vs velocity minus its kinetic part).
    nonlocal: Vec<f64>,
    /// Nonlocal part against the Richardson value (4 D(h/2) - D(h)) / 3, O(h^4).
    richardson: Vec<f64>,
    /// max over k, b of max |p_NL| / max |p| (size of the nonlocal term).
    nonlocal_share: f64,
}

fn check(crystal: &Crystal<'_>, ks: &[[f64; 3]], steps: &[f64]) -> Errors {
    let mut errors = Errors {
        full: vec![0.; steps.len()],
        nonlocal: vec![0.; steps.len()],
        richardson: vec![0.; steps.len()],
        nonlocal_share: 0.,
    };
    let species = crystal.species();
    for &k in ks {
        let basis = PlaneWaveBasis::at_k(&crystal.cell, Hartree(crystal.cutoff), k).unwrap();
        let x = crystal.vectors(&basis);
        let momenta: Vec<[f64; 3]> = basis.momenta().iter().map(|g| g.cartesian).collect();
        let grad = projector_k_gradient(&basis, crystal.cell.volume_bohr3(), &species).unwrap();
        let full = velocity(&momenta, Some(&grad), &x, &x).unwrap();
        let kinetic = velocity(&momenta, None, &x, &x).unwrap();
        for b in 0..3 {
            let scale = max_abs(&full[b]);
            let nl = diff(&full[b], &kinetic[b]);
            let nl_scale = max_abs(&nl);
            errors.nonlocal_share = errors.nonlocal_share.max(nl_scale / scale);
            let central = |h: f64| {
                let mut fd = diff(
                    &crystal.shifted(&basis, &x, b, h),
                    &crystal.shifted(&basis, &x, b, -h),
                );
                for z in fd.as_mut_slice() {
                    *z /= 2. * h;
                }
                fd
            };
            for (s, &h) in steps.iter().enumerate() {
                let fd = central(h);
                let e_full = max_diff(&fd, &full[b]) / scale;
                let fd_nl = diff(&fd, &kinetic[b]);
                let e_nl = max_diff(&fd_nl, &nl) / nl_scale;
                errors.full[s] = errors.full[s].max(e_full);
                errors.nonlocal[s] = errors.nonlocal[s].max(e_nl);
                let half = central(h / 2.);
                let mut r = half.clone();
                for (z, (a, c)) in r
                    .as_mut_slice()
                    .iter_mut()
                    .zip(half.as_slice().iter().zip(fd.as_slice()))
                {
                    *z = (a * 4. - c) / 3.;
                }
                let e_r = max_diff(&diff(&r, &kinetic[b]), &nl) / nl_scale;
                errors.richardson[s] = errors.richardson[s].max(e_r);
            }
        }
    }
    println!(
        "{}: {} bands (i == j and i != j), k {:?}, |p_NL|/|p| up to {:.3e}",
        crystal.name, crystal.bands, ks, errors.nonlocal_share
    );
    println!(
        "  h (Bohr^-1) | max rel err full p | max rel err nonlocal part | nonlocal, Richardson (h, h/2)"
    );
    for (s, &h) in steps.iter().enumerate() {
        println!(
            "  {h:9.1e} | {:.3e} | {:.3e} | {:.3e}",
            errors.full[s], errors.nonlocal[s], errors.richardson[s]
        );
    }
    errors
}

fn local_components() -> HashMap<[i32; 3], C> {
    let mut map = HashMap::new();
    for (h, v) in [
        ([1, 0, 0], C::new(-0.6, 0.)),
        ([0, 1, 0], C::new(-0.5, 0.1)),
        ([0, 0, 1], C::new(-0.55, 0.)),
        ([1, 1, 0], C::new(-0.2, 0.1)),
        ([0, 1, -1], C::new(-0.15, -0.05)),
    ] {
        *map.entry(h).or_insert(C::default()) += v;
        *map.entry([-h[0], -h[1], -h[2]]).or_insert(C::default()) += v.conj();
    }
    map
}

const STEPS: [f64; 4] = [2e-2, 1e-2, 1e-3, 1e-4];

#[test]
fn velocity_matches_central_differences_of_the_hamiltonian_multi_projector() {
    let pp = multi_projector_upf();
    let cell =
        Cell::from_user_rows_bohr([[5.0, 0.0, 0.0], [0.6, 4.8, 0.0], [0.3, 0.5, 5.2]]).unwrap();
    let positions = vec![
        cell.fractional_to_cartesian_bohr([0.13, 0.27, 0.61]),
        cell.fractional_to_cartesian_bohr([0.58, 0.71, 0.09]),
    ];
    let crystal = Crystal {
        name: "synthetic s/p/d/f, two projectors per l <= 2, coupled D, skewed cell",
        cell,
        species: vec![(&pp, positions)],
        cutoff: 3.0,
        local: local_components(),
        bands: 8,
    };
    let e = check(
        &crystal,
        &[[0., 0., 0.], [0.21, -0.34, 0.12], [-0.4, 0.15, 0.37]],
        &STEPS,
    );
    // The fixture must exercise the nonlocal velocity.
    assert!(e.nonlocal_share > 1e-2, "{}", e.nonlocal_share);
    // Central differences: error h^2 down to the rounding floor eps / h; the
    // Richardson value removes the h^2 term (an O(1) error in the velocity
    // would stay at that level for every h).
    let best_full = e.full.iter().copied().fold(f64::INFINITY, f64::min);
    let best_nl = e.nonlocal.iter().copied().fold(f64::INFINITY, f64::min);
    let best_r = e.richardson.iter().copied().fold(f64::INFINITY, f64::min);
    assert!(best_full <= 1e-7, "full velocity vs FD {best_full:e}");
    assert!(best_nl <= 1e-7, "nonlocal velocity vs FD {best_nl:e}");
    assert!(
        best_r <= 1e-9,
        "nonlocal velocity vs Richardson FD {best_r:e}"
    );
    // Second order: h = 1e-2 to 1e-3 lowers the truncation error ~100x.
    assert!(e.nonlocal[2] < 0.05 * e.nonlocal[1], "{:?}", e.nonlocal);
}

fn dojo(dir: &std::path::Path, element: &str) -> dft_upf::Pseudopotential {
    let path = dir.join(format!("{element}.upf"));
    let text = std::fs::read_to_string(&path).unwrap_or_else(|e| panic!("{path:?}: {e}"));
    dft_upf::parse_upf(&text).unwrap()
}

#[test]
#[ignore = "external PseudoDojo PBE UPFs; DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-response --test velocity_fd -- --ignored --nocapture"]
fn velocity_matches_central_differences_on_pseudodojo_lif_and_aln() {
    let dir = std::path::PathBuf::from(
        std::env::var_os("DFT_NC_DIR").expect("DFT_NC_DIR=/abs/pseudodojo/ex"),
    );
    let (li, f, al, n) = (
        dojo(&dir, "Li"),
        dojo(&dir, "F"),
        dojo(&dir, "Al"),
        dojo(&dir, "N"),
    );
    // Geometry of tests/reference/response/r2/{lif,aln}-444.toml; 30 Ha = 60 Ry.
    let lif = Cell::from_user_rows_bohr([
        [0.0, 3.80785, 3.80785],
        [3.80785, 0.0, 3.80785],
        [3.80785, 3.80785, 0.0],
    ])
    .unwrap();
    let lif = Crystal {
        name: "PseudoDojo LiF (Li s,p x2; F s,p x2 + d), 30 Ha",
        species: vec![
            (&li, vec![lif.fractional_to_cartesian_bohr([0., 0., 0.])]),
            (&f, vec![lif.fractional_to_cartesian_bohr([0.5, 0.5, 0.5])]),
        ],
        cell: lif,
        cutoff: 30.0,
        local: local_components(),
        bands: 12,
    };
    let ks = [[0., 0., 0.], [0.25, 0., 0.25], [0.25, -0.25, 0.]];
    let e1 = check(&lif, &ks, &STEPS);
    let aln = Cell::from_user_rows_bohr([
        [5.877, 0.0, 0.0],
        [-2.9385, 5.089631298041145, 0.0],
        [0.0, 0.0, 9.409077],
    ])
    .unwrap();
    let frac = |f: [f64; 3]| aln.fractional_to_cartesian_bohr(f);
    let aln = Crystal {
        name: "PseudoDojo AlN wurtzite (Al s,p,d x2; N s,p x2), 30 Ha",
        species: vec![
            (
                &al,
                vec![frac([1. / 3., 2. / 3., 0.]), frac([2. / 3., 1. / 3., 0.5])],
            ),
            (
                &n,
                vec![
                    frac([1. / 3., 2. / 3., 0.382]),
                    frac([2. / 3., 1. / 3., 0.882]),
                ],
            ),
        ],
        cell: aln,
        cutoff: 30.0,
        local: local_components(),
        bands: 16,
    };
    let e2 = check(&aln, &[[0., 0., 0.], [0.25, 0., 1. / 3.]], &STEPS);
    for e in [e1, e2] {
        let best = e.richardson.iter().copied().fold(f64::INFINITY, f64::min);
        assert!(best <= 1e-7, "nonlocal velocity vs Richardson FD {best:e}");
    }
}

/// R2-VEL-CHECK (b): the NSCF Hamiltonian's k projectors interpolate GK14
/// radial tables (`refstart::k_radial_tables`: node values = the quadrature
/// with the run's even-mesh rule, default Simpson 3/8, 8-point Lagrange, step
/// 0.005 Bohr^-1); the velocity takes `projector_transform` and
/// `projector_transform_slope` (exact quadrature, no table). Compare the
/// table value and its q derivative (central difference of the interpolant,
/// i.e. the derivative of the H actually solved) with the exact quadrature
/// on 0 <= q <= q_max (q_max = 8 Bohr^-1 covers |k+G| at 30 Ha).
#[test]
#[ignore = "external PseudoDojo PBE UPFs; DFT_NC_DIR=/abs/pseudodojo/ex cargo test --release -p dft-response --test velocity_fd tables -- --ignored --nocapture"]
fn nscf_radial_tables_match_the_velocity_quadrature_on_pseudodojo() {
    use dft_upf::radial::{
        RADIAL_TABLE_POINTS, RADIAL_TABLE_STEP, RadialEvenRule, RadialTransformTable,
    };
    let dir = std::path::PathBuf::from(
        std::env::var_os("DFT_NC_DIR").expect("DFT_NC_DIR=/abs/pseudodojo/ex"),
    );
    let q_max = 8.0;
    let samples = 11_000;
    let mut worst_value = 0_f64;
    let mut worst_slope = 0_f64;
    for element in ["Li", "F", "Al", "N"] {
        let pp = dojo(&dir, element);
        for (i, p) in pp.projectors().iter().enumerate() {
            let table = RadialTransformTable::new(
                pp.mesh(),
                p.l(),
                p.radial_times_r(),
                RadialEvenRule::Simpson3_8,
                q_max,
                RADIAL_TABLE_STEP,
                RADIAL_TABLE_POINTS,
            )
            .unwrap();
            let (mut peak, mut peak_slope) = (0_f64, 0_f64);
            let (mut err_value, mut err_slope) = (0_f64, 0_f64);
            let d = 1e-4;
            for s in 1..samples {
                // Irrational offset: never on a node.
                let q = (s as f64 + 0.377_219_431) * (q_max - 2. * d) / samples as f64 + d;
                let exact = pp.projector_transform(i, q).unwrap();
                let slope = dft_physics::stress::projector_transform_slope(&pp, i, q).unwrap();
                let tabulated = table.eval(q).unwrap();
                let table_slope =
                    (table.eval(q + d).unwrap() - table.eval(q - d).unwrap()) / (2. * d);
                peak = peak.max(exact.abs());
                peak_slope = peak_slope.max(slope.abs());
                err_value = err_value.max((tabulated - exact).abs());
                err_slope = err_slope.max((table_slope - slope).abs());
            }
            println!(
                "{element} projector {i} (l = {}): |beta|max {peak:.4e}, table abs err {err_value:.2e} \
                 (rel {:.2e}); |dbeta/dq|max {peak_slope:.4e}, d(table)/dq vs exact slope abs err \
                 {err_slope:.2e} (rel {:.2e})",
                p.l(),
                err_value / peak,
                err_slope / peak_slope
            );
            worst_value = worst_value.max(err_value / peak);
            worst_slope = worst_slope.max(err_slope / peak_slope);
        }
    }
    println!("worst relative: value {worst_value:.2e}, slope {worst_slope:.2e}");
    assert!(worst_value <= 1e-12, "{worst_value:e}");
    // The interpolant's derivative carries the central-difference error
    // (d = 1e-4: d^2 / 6 x beta''' ~ 1e-9 relative) and rounding (eps / d).
    assert!(worst_slope <= 1e-7, "{worst_slope:e}");
}
