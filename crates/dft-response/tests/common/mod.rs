//! Synthetic band problems for the R2 response oracles: explicit dense
//! H_k = |k+G|^2/2 + V(G - G') + sum_atom B D B^H in a complete plane-wave
//! basis, diagonalized with the faer Hermitian EVD. No UPF file, no engine
//! code path: the local potential is a handful of Fourier components, the
//! nonlocal part (optional) a synthetic in-memory KB dataset.
#![allow(dead_code)]
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_numerics::linalg::{ComplexMatrix, hermitian_eigh};
use dft_physics::hamiltonian::Species;
use dft_physics::stress::{ProjectorKGradient, projector_k_gradient};
use dft_response::Complex64 as C;
use dft_units::Hartree;
use std::collections::HashMap;

pub fn array(v: &[f64]) -> String {
    v.iter()
        .map(|v| format!("{v:.17e}"))
        .collect::<Vec<_>>()
        .join(" ")
}

/// Synthetic NC dataset with s, p and d projectors (not a real element).
pub fn synthetic_upf() -> dft_upf::Pseudopotential {
    let points = 401;
    let end = 12.0;
    let r0 = 0.02;
    let dx = (1.0_f64 + end / r0).ln() / (points - 1) as f64;
    let r: Vec<_> = (0..points).map(|i| r0 * (i as f64 * dx).exp_m1()).collect();
    let rab: Vec<_> = r.iter().map(|r| dx * (r + r0)).collect();
    let local: Vec<_> = r.iter().map(|_| 0.0).collect();
    let mut beta = String::new();
    for (i, l) in [0, 1, 2].into_iter().enumerate() {
        let values: Vec<_> = r
            .iter()
            .map(|r| r.powi(l + 1) * (-r * r / (1.0 + 0.4 * i as f64)).exp())
            .collect();
        beta.push_str(&format!(
            "<PP_BETA.{n} index=\"{n}\" angular_momentum=\"{l}\" cutoff_radius=\"12\" cutoff_radius_index=\"{points}\">{}</PP_BETA.{n}>",
            array(&values),
            n = i + 1
        ));
    }
    let text = format!(
        r#"<UPF version="2.0.1"><PP_INFO>SYNTHETIC, no real data.</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ" z_valence="2" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="2" l_max_rho="0" l_local="-1" mesh_size="{points}" number_of_wfc="0" number_of_proj="3" />
<PP_MESH dx="{dx}" mesh="{points}" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH><PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>0.9 0 0 0 -0.7 0 0 0 0.5</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,
        array(&r),
        array(&rab),
        array(&local),
        array(&vec![0.0; points]),
    );
    dft_upf::parse_upf(&text).unwrap()
}

/// Fourier components V(h) (Ha) of a real periodic potential; the map holds
/// both h and -h (V(-h) = conj V(h)).
pub fn potential(components: &[([i32; 3], C)]) -> HashMap<[i32; 3], C> {
    let mut map = HashMap::new();
    for &(h, v) in components {
        *map.entry(h).or_insert(C::default()) += v;
        *map.entry([-h[0], -h[1], -h[2]]).or_insert(C::default()) += v.conj();
    }
    map
}

pub struct Problem<'a> {
    pub cell: Cell,
    pub cutoff: f64,
    pub local: HashMap<[i32; 3], C>,
    /// Nonlocal atoms: dataset and Cartesian positions (Bohr).
    pub nonlocal: Option<(&'a dft_upf::Pseudopotential, Vec<[f64; 3]>)>,
}

pub struct Solved {
    pub basis: PlaneWaveBasis,
    pub miller: Vec<[i32; 3]>,
    pub momenta: Vec<[f64; 3]>,
    pub energies: Vec<f64>,
    pub vectors: ComplexMatrix,
    pub projectors: Option<ProjectorKGradient>,
}

impl Problem<'_> {
    pub fn dims(&self) -> [usize; 3] {
        let basis = PlaneWaveBasis::new(&self.cell, Hartree(self.cutoff)).unwrap();
        let density = basis.density_basis(&self.cell).unwrap();
        FftGridSpec::new(&basis, &density, None).unwrap().dims()
    }
    pub fn volume(&self) -> f64 {
        self.cell.volume_bohr3()
    }
    pub fn solve(&self, k: [f64; 3]) -> Solved {
        self.solve_with(k, &HashMap::new())
    }
    /// Complete EVD with an extra local potential (Miller map of this cell).
    pub fn solve_with(&self, k: [f64; 3], extra: &HashMap<[i32; 3], C>) -> Solved {
        let (mut solved, h) = self.hamiltonian(k, extra);
        let (energies, vectors) = hermitian_eigh(&h).unwrap();
        solved.energies = energies;
        solved.vectors = vectors;
        solved
    }
    /// The explicit Hermitian H_k (rows/columns in basis order) and the
    /// basis data (energies and vectors left empty).
    pub fn hamiltonian(
        &self,
        k: [f64; 3],
        extra: &HashMap<[i32; 3], C>,
    ) -> (Solved, ComplexMatrix) {
        let basis = PlaneWaveBasis::at_k(&self.cell, Hartree(self.cutoff), k).unwrap();
        let miller: Vec<[i32; 3]> = basis.vectors().iter().map(|g| g.miller).collect();
        let momenta: Vec<[f64; 3]> = basis.momenta().iter().map(|g| g.cartesian).collect();
        let n = miller.len();
        let mut h = ComplexMatrix::zeros(n, n).unwrap();
        for j in 0..n {
            for i in 0..n {
                let d = [
                    miller[i][0] - miller[j][0],
                    miller[i][1] - miller[j][1],
                    miller[i][2] - miller[j][2],
                ];
                let mut v = self.local.get(&d).copied().unwrap_or_default();
                v += extra.get(&d).copied().unwrap_or_default();
                h[(i, j)] = v;
            }
            h[(j, j)] += basis.momenta()[j].norm_squared / 2.;
        }
        let projectors = self.nonlocal.as_ref().map(|(pp, positions)| {
            let p = projector_k_gradient(
                &basis,
                self.volume(),
                &[Species {
                    potential: pp,
                    positions_bohr: positions,
                }],
            )
            .unwrap();
            let b = &p.values;
            let mut offset = 0;
            for (&count, d) in p.channels.iter().zip(&p.couplings) {
                for a in 0..count {
                    for c in 0..count {
                        let coupling = d[a * count + c];
                        if coupling == 0. {
                            continue;
                        }
                        for j in 0..n {
                            let right = b[(j, offset + c)].conj() * coupling;
                            for i in 0..n {
                                h[(i, j)] += b[(i, offset + a)] * right;
                            }
                        }
                    }
                }
                offset += count;
            }
            p
        });
        for j in 0..n {
            h[(j, j)] = C::new(h[(j, j)].re, 0.);
            for i in 0..j {
                let a = (h[(i, j)] + h[(j, i)].conj()) * 0.5;
                h[(i, j)] = a;
                h[(j, i)] = a.conj();
            }
        }
        (
            Solved {
                basis,
                miller,
                momenta,
                energies: Vec::new(),
                vectors: ComplexMatrix::zeros(0, 0).unwrap(),
                projectors,
            },
            h,
        )
    }
}

pub fn gamma_mesh(n: [usize; 3]) -> Vec<[f64; 3]> {
    let mut out = Vec::new();
    for i in 0..n[0] {
        for j in 0..n[1] {
            for k in 0..n[2] {
                out.push([
                    i as f64 / n[0] as f64,
                    j as f64 / n[1] as f64,
                    k as f64 / n[2] as f64,
                ]);
            }
        }
    }
    out
}

/// Integer shift G0 with `left + q = right + G0` (fractional, within 1e-9).
pub fn shift(left: [f64; 3], q: [f64; 3], right: [f64; 3]) -> [i32; 3] {
    std::array::from_fn(|i| {
        let s = left[i] + q[i] - right[i];
        let r = s.round();
        assert!((s - r).abs() < 1e-9, "k + q is not k' + G: {s}");
        r as i32
    })
}

pub fn max_abs(m: &ComplexMatrix) -> f64 {
    m.as_slice().iter().map(|z| z.norm()).fold(0., f64::max)
}
