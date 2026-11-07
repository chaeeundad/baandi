//! Bitwise setup regression against the pre-P28A implementation (534595c).
use dft_basis::{GVector, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_physics::hamiltonian::{self, HamiltonianError, NonlocalOperator};
use dft_units::Hartree;
use dft_upf::{Pseudopotential, radial::real_spherical_harmonic};
use num_complex::Complex64;
fn finite(z: &Complex64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn volume(v: f64) -> Result<(), HamiltonianError> {
    if v.is_finite() && v > 0.0 {
        Ok(())
    } else {
        Err(HamiltonianError::Volume)
    }
}
fn phase(g: [f64; 3], position: [f64; 3]) -> Result<Complex64, HamiltonianError> {
    let angle: f64 = g.iter().zip(position).map(|(g, r)| g * r).sum();
    if !angle.is_finite() {
        return Err(HamiltonianError::NonFinite);
    }
    Ok(Complex64::from_polar(1.0, -angle))
}
/// Positions are Cartesian, in Bohr. Multiple entries may use the same potential.
struct Species<'a> {
    pub potential: &'a Pseudopotential,
    pub positions_bohr: &'a [[f64; 3]],
}
/// S(G) = sum_I exp(-i G.tau_I).
pub fn structure_factor(
    g: [f64; 3],
    positions_bohr: &[[f64; 3]],
) -> Result<Complex64, HamiltonianError> {
    positions_bohr.iter().map(|&r| phase(g, r)).sum()
}
/// Coefficients in the supplied density-sphere order; zero outside that sphere.
/// The caller supplies the density basis belonging to its wavefunction basis/cell.
fn old_ionic_local(
    density: &PlaneWaveBasis,
    volume_bohr3: f64,
    species: &[Species<'_>],
) -> Result<Vec<Complex64>, HamiltonianError> {
    volume(volume_bohr3)?;
    let mut result = vec![Complex64::default(); density.vectors().len()];
    for s in species {
        let local = s.potential.local_transform(volume_bohr3, None)?;
        for (out, g) in result.iter_mut().zip(density.vectors()) {
            let v = if g.miller == [0; 3] {
                local.alpha_z()?
            } else {
                local.at_g(g.norm_squared.sqrt())?
            };
            *out += v * structure_factor(g.cartesian, s.positions_bohr)?;
        }
    }
    if !result.iter().all(finite) {
        return Err(HamiltonianError::NonFinite);
    }
    Ok(result)
}
/// A single atom's projector columns, ordered by radial index then m=-l..l.
/// D is expanded only over channels with equal l and m.
pub struct AtomProjectors {
    columns: Vec<Vec<Complex64>>,
    coupling: Vec<f64>,
}
impl AtomProjectors {
    pub fn columns(&self) -> &[Vec<Complex64>] {
        &self.columns
    }
    pub fn coupling_hartree(&self) -> &[f64] {
        &self.coupling
    }
}
/// O(N_G * N_projectors) storage, never an N_G by N_G matrix.
struct OldNonlocalOperator {
    basis: Vec<GVector>,
    atoms: Vec<AtomProjectors>,
    max_channels: usize,
}
impl OldNonlocalOperator {
    /// Build projectors at q=k+G from the basis momenta, including exp(-i q.tau).
    /// Gamma uses exactly the same radial/angular and structure-factor path.
    pub fn new(
        basis: &PlaneWaveBasis,
        volume_bohr3: f64,
        species: &[Species<'_>],
    ) -> Result<Self, HamiltonianError> {
        volume(volume_bohr3)?;
        let mut atoms = Vec::new();
        let mut max_channels = 0;
        for s in species {
            let pp = s.potential;
            let np = pp.projectors().len();
            for i in 0..np {
                for j in 0..np {
                    if pp.projectors()[i].l() != pp.projectors()[j].l()
                        && pp.dij_hartree()[i * np + j] != 0.0
                    {
                        return Err(HamiltonianError::AngularCoupling);
                    }
                }
            }
            let channels: Vec<_> = pp
                .projectors()
                .iter()
                .enumerate()
                .flat_map(|(i, p)| {
                    let l = p.l() as isize;
                    (-l..=l).map(move |m| (i, p.l(), m))
                })
                .collect();
            let nc = channels.len();
            max_channels = max_channels.max(nc);
            let mut table = vec![vec![Complex64::default(); basis.momenta().len()]; nc];
            // Radial transforms are evaluated once per radial projector and G,
            // independent of the number of atoms and magnetic channels.
            for (ig, g) in basis.momenta().iter().enumerate() {
                let radial: Vec<_> = (0..np)
                    .map(|i| pp.projector_transform(i, g.norm_squared.sqrt()))
                    .collect::<Result<_, _>>()?;
                // Direction at zero is arbitrary: all l>0 radial transforms vanish.
                let direction = if g.norm_squared == 0.0 {
                    [0.0, 0.0, 1.0]
                } else {
                    g.cartesian
                };
                for (c, &(i, l, m)) in channels.iter().enumerate() {
                    table[c][ig] = Complex64::new(0.0, -1.0).powu(l as u32)
                        * real_spherical_harmonic(l, m, direction)?
                        * radial[i]
                        / volume_bohr3.sqrt();
                }
            }
            let mut coupling = vec![0.0; nc * nc];
            for (a, &(i, l, m)) in channels.iter().enumerate() {
                for (b, &(j, ll, mm)) in channels.iter().enumerate() {
                    if l == ll && m == mm {
                        coupling[a * nc + b] = pp.dij_hartree()[i * np + j];
                    }
                }
            }
            for &r in s.positions_bohr {
                let phases = basis
                    .momenta()
                    .iter()
                    .map(|g| phase(g.cartesian, r))
                    .collect::<Result<Vec<_>, _>>()?;
                let columns: Vec<Vec<_>> = table
                    .iter()
                    .map(|col| {
                        col.iter()
                            .zip(&phases)
                            .map(|(p, phase)| p * phase)
                            .collect()
                    })
                    .collect();
                if !columns.iter().flatten().all(finite) {
                    return Err(HamiltonianError::NonFinite);
                }
                atoms.push(AtomProjectors {
                    columns,
                    coupling: coupling.clone(),
                });
            }
        }
        Ok(Self {
            basis: basis.momenta().to_vec(),
            atoms,
            max_channels,
        })
    }
}
fn pseudo(ls: &[usize], n: usize) -> dft_upf::Pseudopotential {
    let a = |f: &dyn Fn(f64) -> f64| {
        (0..n)
            .map(|i| format!("{:.17e}", f(i as f64 * 0.15)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let mut beta = String::new();
    for (j, &l) in ls.iter().enumerate() {
        let i = j + 1;
        beta.push_str(&format!("<PP_BETA.{i} index=\"{i}\" angular_momentum=\"{l}\" cutoff_radius=\"12\" cutoff_radius_index=\"{n}\">{}</PP_BETA.{i}>",a(&|r| r.powi(l as i32+1)*(-r*r/(1.+j as f64*0.1)).exp())));
    }
    let dij = (0..ls.len())
        .flat_map(|i| {
            (0..ls.len()).map(move |j| {
                if i == j {
                    1. + i as f64 * 0.1
                } else if ls[i] == ls[j] {
                    0.2
                } else {
                    0.
                }
            })
        })
        .map(|v| v.to_string())
        .collect::<Vec<_>>()
        .join(" ");
    dft_upf::parse_upf(&format!(r#"<UPF version="2.0.1"><PP_INFO>Synthetic R3 fixture</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ" z_valence="4" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="3" l_max_rho="0" l_local="-1" mesh_size="{n}" number_of_wfc="0" number_of_proj="{}"/>
<PP_MESH dx="0.15" mesh="{n}" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH><PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>{dij}</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,ls.len(),a(&|r|r),a(&|_|0.15),a(&|r| -2.0 / (1.0 + r)),a(&|_|0.)).replace("<PP_NONLOCAL><PP_DIJ></PP_DIJ></PP_NONLOCAL>", "<PP_NONLOCAL/>")).unwrap()
}

fn assert_bits(a: &[Complex64], b: &[Complex64]) {
    assert_eq!(a.len(), b.len());
    for (a, b) in a.iter().zip(b) {
        assert_eq!(
            (a.re.to_bits(), a.im.to_bits()),
            (b.re.to_bits(), b.im.to_bits())
        );
    }
}

#[test]
fn setup_matches_old_bits() {
    let shapes = [
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]],
        [[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]],
        [[4., 0., 0.], [0., 5., 0.], [0., 0., 14.]],
        [[7., 0., 0.], [0., 8., 0.], [0.3, 0., 18.]],
    ];
    let p = pseudo(&[0, 0, 1, 1, 2, 2, 3, 3], 81);
    let q = pseudo(&[0, 1], 66);
    let mut cases = 0;
    for rows in shapes {
        let cell = Cell::from_user_rows_bohr(rows).unwrap();
        for atoms in 1usize..=8 {
            let positions: Vec<_> = (0..atoms)
                .map(|i| {
                    cell.fractional_to_cartesian_bohr([
                        ((i * 7 + 1) % 23) as f64 / 23.,
                        ((i * 11 + 2) % 29) as f64 / 29.,
                        ((i * 13 + 3) % 31) as f64 / 31.,
                    ])
                })
                .collect();
            let split = atoms.div_ceil(2);
            let species = [
                Species {
                    potential: &p,
                    positions_bohr: &positions[..split],
                },
                Species {
                    potential: &q,
                    positions_bohr: &positions[split..],
                },
            ];
            let actual_species: Vec<_> = species
                .iter()
                .map(|s| hamiltonian::Species {
                    potential: s.potential,
                    positions_bohr: s.positions_bohr,
                })
                .collect();
            for k in [[0.; 3], [0.25; 3], [0.13, -0.21, 0.07]] {
                let basis = PlaneWaveBasis::at_k(&cell, Hartree(0.6), k).unwrap();
                let density = basis.density_basis(&cell).unwrap();
                assert_bits(
                    &hamiltonian::ionic_local(&density, cell.volume_bohr3(), &actual_species)
                        .unwrap(),
                    &old_ionic_local(&density, cell.volume_bohr3(), &species).unwrap(),
                );
                let old = OldNonlocalOperator::new(&basis, cell.volume_bohr3(), &species).unwrap();
                let new =
                    NonlocalOperator::new(&basis, cell.volume_bohr3(), &actual_species).unwrap();
                assert_eq!(old.basis, basis.momenta());
                assert_eq!(new.scratch_len(), old.max_channels * 2);
                assert_eq!(new.atoms().len(), old.atoms.len());
                for (a, b) in new.atoms().iter().zip(&old.atoms) {
                    assert_eq!(a.columns().len(), b.columns.len());
                    for (a, b) in a.columns().iter().zip(b.columns()) {
                        assert_bits(a, b);
                    }
                    assert_eq!(
                        a.coupling_hartree()
                            .iter()
                            .map(|x| x.to_bits())
                            .collect::<Vec<_>>(),
                        b.coupling_hartree()
                            .iter()
                            .map(|x| x.to_bits())
                            .collect::<Vec<_>>()
                    );
                }
                cases += 1;
            }
        }
    }
    println!("{cases} shape/atom/k cases: local and all nonlocal components bit-identical");
}

#[test]
fn local_mesh_panels_and_empty_projectors() {
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6.000000000000001, 0.], [0., 0., 6.]])
        .unwrap();
    let basis = PlaneWaveBasis::new(&cell, Hartree(1.)).unwrap();
    let density = basis.density_basis(&cell).unwrap();
    for n in [2, 3, 4, 5, 6, 65, 66, 67, 80, 81] {
        let pp = pseudo(&[], n);
        let positions = [[0.; 3], [0.1, -0.3, 0.7]];
        let old = [Species {
            potential: &pp,
            positions_bohr: &positions,
        }];
        let new = [hamiltonian::Species {
            potential: &pp,
            positions_bohr: &positions,
        }];
        assert_bits(
            &old_ionic_local(&density, cell.volume_bohr3(), &old).unwrap(),
            &hamiltonian::ionic_local(&density, cell.volume_bohr3(), &new).unwrap(),
        );
        let nl = NonlocalOperator::new(&basis, cell.volume_bohr3(), &new).unwrap();
        assert_eq!(nl.scratch_len(), 0);
        assert_eq!(nl.atoms().len(), 2);
        assert!(
            nl.atoms()
                .iter()
                .all(|a| a.columns().is_empty() && a.coupling_hartree().is_empty())
        );
    }
}
