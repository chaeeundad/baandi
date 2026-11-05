//! Independent setup-time model record. Exact parsed PP content is compared
//! because the engine never receives original-file SHA-256 values. No hashing,
//! I/O or work in SCF numerical kernels; PP storage is shared by result clones.
use crate::collinear::SpinSector;
use crate::finite_temperature::OccupationWorkspace;
use crate::{Functional, ScfOptions, SpeciesPotential};
use dft_geometry::{Atom, Cell, kpoints::KPoint, symmetry::SymmetryOperation};
use dft_physics::occupations::Constraint;
use std::sync::Arc;
fn bits(values: &[f64]) -> Vec<u64> {
    values.iter().map(|x| x.to_bits()).collect()
}
#[derive(Clone, PartialEq)]
struct PseudoIdentity {
    label: String,
    element: String,
    functional: String,
    charge: u64,
    radial: Vec<Vec<u64>>,
    core: Option<Vec<u64>>,
    projectors: Vec<(usize, usize, Vec<u64>)>,
}
fn species(potentials: &[SpeciesPotential<'_>]) -> Arc<[PseudoIdentity]> {
    let mut result: Vec<_> = potentials
        .iter()
        .map(|s| {
            let p = s.potential;
            PseudoIdentity {
                label: s.label.into(),
                element: p.element().into(),
                functional: p.functional().into(),
                charge: p.z_valence().to_bits(),
                radial: vec![
                    bits(p.mesh().r()),
                    bits(p.mesh().rab()),
                    bits(p.local_hartree()),
                    bits(p.dij_hartree()),
                ],
                core: p.core_density().map(bits),
                projectors: p
                    .projectors()
                    .iter()
                    .map(|p| (p.l(), p.cutoff_index(), bits(p.radial_times_r())))
                    .collect(),
            }
        })
        .collect();
    result.sort_by(|a, b| a.label.cmp(&b.label));
    result.into()
}
#[derive(Clone, PartialEq)]
pub(crate) struct ModelDescription {
    cell: [[u64; 3]; 3],
    atoms: Vec<(String, [u64; 3])>,
    species: Arc<[PseudoIdentity]>,
    cutoff: u64,
    strain_reference: Option<[[u64; 3]; 3]>,
    functional: Functional,
    pub(crate) grid: [usize; 3],
    pub(crate) floor: Option<f64>,
    pub(crate) symmetry: Vec<SymmetryOperation>,
    pub(crate) points: Vec<KPoint>,
}
impl ModelDescription {
    /// Fixed-cell density transfer allows only ordered atomic positions to differ.
    pub(crate) fn transfer_matches(&self, other: &Self) -> bool {
        let mut source = other.clone();
        if self.atoms.len() != source.atoms.len()
            || self
                .atoms
                .iter()
                .zip(&source.atoms)
                .any(|(a, b)| a.0 != b.0)
        {
            return false;
        }
        source.atoms.clone_from(&self.atoms);
        *self == source
    }
    pub(crate) fn new(
        cell: &Cell,
        atoms: &[Atom],
        potentials: &[SpeciesPotential<'_>],
        options: &ScfOptions,
        grid: [usize; 3],
        points: Vec<KPoint>,
    ) -> Self {
        Self {
            cell: cell.lattice_columns_bohr().map(|row| row.map(f64::to_bits)),
            atoms: atoms
                .iter()
                .map(|a| (a.species.clone(), a.fractional_position.map(f64::to_bits)))
                .collect(),
            species: species(potentials),
            cutoff: options.cutoff_hartree.to_bits(),
            strain_reference: options
                .strain_reference
                .as_ref()
                .map(|c| c.lattice_columns_bohr().map(|r| r.map(f64::to_bits))),
            functional: options.functional,
            grid,
            floor: None,
            symmetry: Vec::new(),
            points,
        }
    }
    pub(crate) fn with_policy(&self, policy: OccupationIdentity) -> ForceIdentity {
        ForceIdentity {
            model: self.clone(),
            policy,
        }
    }
}
#[derive(Clone, PartialEq)]
pub(crate) struct OccupationIdentity {
    fixed: Vec<Vec<u64>>,
    finite: Option<(u64, Constraint, usize)>,
    sector: Option<(usize, usize, [u64; 3])>,
}
impl OccupationIdentity {
    pub(crate) fn new(
        fixed: &[&[f64]],
        finite: Option<&OccupationWorkspace>,
        sector: Option<&SpinSector>,
    ) -> Self {
        Self {
            fixed: if finite.is_some() {
                Vec::new()
            } else {
                fixed.iter().map(|f| bits(f)).collect()
            },
            finite: finite.map(|f| (f.theta_ha().to_bits(), f.constraint(), f.bands)),
            sector: sector.map(|s| (s.n_up, s.n_down, s.axis_cartesian.map(f64::to_bits))),
        }
    }
}
#[derive(Clone)]
pub(crate) struct ForceIdentity {
    model: ModelDescription,
    policy: OccupationIdentity,
}
impl std::fmt::Debug for ForceIdentity {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str("ForceIdentity { independent_setup_record: true }")
    }
}
impl ForceIdentity {
    pub(crate) fn matches(&self, model: &ModelDescription, policy: &OccupationIdentity) -> bool {
        self.model == *model && self.policy == *policy
    }
}
