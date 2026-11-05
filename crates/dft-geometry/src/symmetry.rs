//! Crystal operations in the input fractional basis (no cell standardization).
use crate::{
    Atom, Cell,
    kpoints::{self, KPoint},
};
use nalgebra::{Matrix3, Vector3};
use thiserror::Error;

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SpinAction {
    Identity,
    Swap,
}

#[derive(Clone, Debug, PartialEq)]
pub struct SymmetryOperation {
    /// Matrix rows; x' = R x + tau, modulo direct lattice integers.
    pub rotation: [[i32; 3]; 3],
    pub translation: [f64; 3],
    pub atom_permutation: Vec<usize>,
    pub spin_action: SpinAction,
    pub antiunitary: bool,
    pub model_invariance_proof: &'static str,
}
impl SymmetryOperation {
    pub fn reciprocal_action(&self, k: [f64; 3]) -> [f64; 3] {
        let r = Matrix3::from_fn(|i, j| f64::from(self.rotation[i][j]));
        let v = r
            .try_inverse()
            .expect("validated unimodular rotation")
            .transpose()
            * Vector3::from(k);
        let sign = if self.antiunitary { -1.0 } else { 1.0 };
        [sign * v[0], sign * v[1], sign * v[2]]
    }
}
#[derive(Debug, Error, PartialEq)]
#[error("symmetry: {0}")]
pub struct SymmetryError(pub &'static str);

fn matrix(r: [[i32; 3]; 3]) -> Matrix3<f64> {
    Matrix3::from_fn(|i, j| f64::from(r[i][j]))
}
fn periodic_match(cell: &Cell, a: [f64; 3], b: [f64; 3], tol: f64) -> bool {
    // A match within tol must have each fractional residual bounded by
    // tol * ||row_i(A^-1)||. Enumerate all possible images, including skew cells.
    let d = Vector3::from(a) - Vector3::from(b);
    let bounds: [f64; 3] = std::array::from_fn(|i| tol * cell.inverse.row(i).norm());
    let lo: [i32; 3] = std::array::from_fn(|i| (d[i] - bounds[i]).ceil() as i32);
    let hi: [i32; 3] = std::array::from_fn(|i| (d[i] + bounds[i]).floor() as i32);
    for x in lo[0]..=hi[0] {
        for y in lo[1]..=hi[1] {
            for z in lo[2]..=hi[2] {
                if (cell.lattice * (d - Vector3::new(x as f64, y as f64, z as f64))).norm() <= tol {
                    return true;
                }
            }
        }
    }
    false
}

/// Exhaustive metric-bounded integer column search. Atomic tolerance is a
/// Cartesian distance in Bohr; fractional component bounds are tol*||A^-1_i||.
/// Metric products allow tol*(|a_i|+|a_j|)+tol^2, in Bohr^2.
/// Fail closed if the search budget is exceeded or atom matching is ambiguous.
pub fn detect(
    cell: &Cell,
    atoms: &[Atom],
    tolerance_bohr: f64,
) -> Result<Vec<SymmetryOperation>, SymmetryError> {
    let tol = tolerance_bohr;
    if !tol.is_finite()
        || tol <= 0.0
        || atoms.is_empty()
        || atoms
            .iter()
            .any(|a| a.fractional_position.iter().any(|v| !v.is_finite()))
    {
        return Err(SymmetryError(
            "finite positive tolerance and nonempty finite structure required",
        ));
    }
    let lengths: [f64; 3] = std::array::from_fn(|i| cell.lattice.column(i).norm());
    if (0..3).any(|i| tol * cell.inverse.row(i).norm() >= 0.1) {
        return Err(SymmetryError("tolerance too large for cell"));
    }
    let metric = cell.lattice.transpose() * cell.lattice;
    let mut columns: [Vec<Vector3<i32>>; 3] = std::array::from_fn(|_| Vec::new());
    for j in 0..3 {
        let bounds: [i32; 3] = std::array::from_fn(|i| {
            ((lengths[j] + tol) * cell.inverse.row(i).norm()).ceil() as i32
        });
        let budget = bounds
            .iter()
            .try_fold(1_u64, |n, &b| n.checked_mul(2 * b as u64 + 1));
        if budget.is_none_or(|n| n > 2_000_000) {
            return Err(SymmetryError(
                "integer search budget exceeded; supply a better conditioned lattice basis",
            ));
        }
        for x in -bounds[0]..=bounds[0] {
            for y in -bounds[1]..=bounds[1] {
                for z in -bounds[2]..=bounds[2] {
                    let v = Vector3::new(x, y, z);
                    if ((cell.lattice * v.cast::<f64>()).norm() - lengths[j]).abs() <= tol {
                        columns[j].push(v);
                    }
                }
            }
        }
    }
    let mut result = Vec::new();
    for a in &columns[0] {
        for b in &columns[1] {
            for c in &columns[2] {
                let ri = Matrix3::from_columns(&[*a, *b, *c]);
                let r = ri.cast::<f64>();
                if (r.determinant().abs() - 1.0).abs() > 1e-8 {
                    continue;
                }
                let delta = r.transpose() * metric * r - metric;
                if (0..3).any(|i| {
                    (0..3)
                        .any(|j| delta[(i, j)].abs() > tol * (lengths[i] + lengths[j]) + tol * tol)
                }) {
                    continue;
                }
                let origin = r * Vector3::from(atoms[0].fractional_position);
                for target in atoms.iter().filter(|a| a.species == atoms[0].species) {
                    let translation: [f64; 3] = std::array::from_fn(|i| {
                        (target.fractional_position[i] - origin[i]).rem_euclid(1.0)
                    });
                    let mut permutation = Vec::new();
                    for atom in atoms {
                        let position: [f64; 3] = (r * Vector3::from(atom.fractional_position)
                            + Vector3::from(translation))
                        .into();
                        let matches: Vec<_> = atoms
                            .iter()
                            .enumerate()
                            .filter(|(_, a)| {
                                a.species == atom.species
                                    && periodic_match(cell, position, a.fractional_position, tol)
                            })
                            .map(|(i, _)| i)
                            .collect();
                        if matches.len() > 1 {
                            return Err(SymmetryError(
                                "ambiguous atom permutation at requested tolerance",
                            ));
                        }
                        if matches.is_empty() || permutation.contains(&matches[0]) {
                            break;
                        }
                        permutation.push(matches[0]);
                    }
                    if permutation.len() == atoms.len() {
                        let operation = SymmetryOperation {
                            rotation: std::array::from_fn(|i| std::array::from_fn(|j| ri[(i, j)])),
                            translation,
                            atom_permutation: permutation,
                            spin_action: SpinAction::Identity,
                            antiunitary: false,
                            model_invariance_proof: "species-aware geometry only; engine must validate active model",
                        };
                        if !result.iter().any(|o: &SymmetryOperation| {
                            o.rotation == operation.rotation
                                && periodic_match(cell, o.translation, operation.translation, tol)
                        }) {
                            result.push(operation);
                        }
                    }
                }
            }
        }
    }
    result.sort_by(|a, b| {
        a.rotation
            .cmp(&b.rotation)
            .then_with(|| a.translation.partial_cmp(&b.translation).unwrap())
    });
    validate_group(cell, &result, tol)?;
    Ok(result)
}

/// Closure includes nonsymmorphic translations and the atom action. Approximate
/// positions that yield a non-group at the tolerance boundary are rejected.
pub fn validate_group(
    cell: &Cell,
    operations: &[SymmetryOperation],
    tol: f64,
) -> Result<(), SymmetryError> {
    if operations.is_empty() {
        return Err(SymmetryError("empty operation group"));
    }
    for a in operations {
        for b in operations {
            let r = matrix(a.rotation) * matrix(b.rotation);
            let t: [f64; 3] = (matrix(a.rotation) * Vector3::from(b.translation)
                + Vector3::from(a.translation))
            .into();
            if !operations.iter().any(|c| {
                matrix(c.rotation) == r
                    && periodic_match(cell, c.translation, t, tol)
                    && c.atom_permutation.len() == b.atom_permutation.len()
                    && b.atom_permutation
                        .iter()
                        .enumerate()
                        .all(|(i, &j)| a.atom_permutation.get(j) == c.atom_permutation.get(i))
            }) {
                return Err(SymmetryError(
                    "operations do not close at requested tolerance",
                ));
            }
        }
    }
    Ok(())
}

#[derive(Clone, Debug)]
pub struct KOrbitEntry {
    pub full_id: usize,
    /// Index in the reduced list.
    pub representative_id: usize,
    pub operation_id: usize,
    /// k_full = (+/- R^-T) k_representative + shift.
    pub reciprocal_shift: [i32; 3],
    pub multiplicity: usize,
}
#[derive(Clone, Debug)]
pub struct KOrbitMap {
    pub points: Vec<KPoint>,
    /// Mesh-preserving operations, including explicit antiunitary partners.
    pub operations: Vec<SymmetryOperation>,
    pub entries: Vec<KOrbitEntry>,
    pub phase_convention: &'static str,
}
fn equivalent(a: [f64; 3], b: [f64; 3]) -> bool {
    (0..3).all(|i| {
        let d = a[i] - b[i];
        (d - d.round()).abs() <= 1e-12
    })
}
/// Select the mesh-preserving subgroup, or reject any excluded operation when
/// strict=true. Only uniform, duplicate-free full meshes are accepted.
pub fn reduce_mesh(
    full: &[KPoint],
    operations: &[SymmetryOperation],
    time_reversal: bool,
    strict: bool,
) -> Result<KOrbitMap, SymmetryError> {
    kpoints::validate(full).map_err(|_| SymmetryError("invalid full mesh"))?;
    if full
        .iter()
        .any(|p| (p.weight - 1.0 / full.len() as f64).abs() > 1e-14)
        || full.iter().enumerate().any(|(i, p)| {
            full[..i]
                .iter()
                .any(|q| equivalent(p.fractional, q.fractional))
        })
    {
        return Err(SymmetryError("uniform duplicate-free full mesh required"));
    }
    let mut subgroup = Vec::new();
    for op in operations {
        if op.spin_action != SpinAction::Identity
            || op.antiunitary
            || (matrix(op.rotation).determinant().abs() - 1.0).abs() > 1e-8
        {
            return Err(SymmetryError(
                "unitary unimodular identity-spin operations required",
            ));
        }
        let invariant = full.iter().all(|p| {
            full.iter()
                .any(|q| equivalent(op.reciprocal_action(p.fractional), q.fractional))
        });
        if invariant {
            subgroup.push(op.clone());
        } else if strict {
            return Err(SymmetryError("operation does not preserve shifted mesh"));
        }
    }
    if subgroup.is_empty() {
        return Err(SymmetryError("no mesh-preserving operations"));
    }
    if time_reversal {
        let partners: Vec<_> = subgroup
            .iter()
            .cloned()
            .map(|mut o| {
                o.antiunitary = true;
                o
            })
            .collect();
        if partners.iter().any(|o| {
            full.iter().any(|p| {
                !full
                    .iter()
                    .any(|q| equivalent(o.reciprocal_action(p.fractional), q.fractional))
            })
        }) {
            return Err(SymmetryError("time reversal does not preserve mesh"));
        }
        subgroup.extend(partners);
    }
    let mut points = Vec::new();
    let mut entries: Vec<Option<KOrbitEntry>> = vec![None; full.len()];
    for i in 0..full.len() {
        if entries[i].is_some() {
            continue;
        }
        let rep = points.len();
        let mut count = 0;
        for (j, p) in full.iter().enumerate() {
            if let Some((op_id, op)) = subgroup
                .iter()
                .enumerate()
                .find(|(_, o)| equivalent(o.reciprocal_action(full[i].fractional), p.fractional))
            {
                if entries[j].is_some() {
                    return Err(SymmetryError(
                        "overlapping orbits: operation set is not a group",
                    ));
                }
                let transformed = op.reciprocal_action(full[i].fractional);
                entries[j] = Some(KOrbitEntry {
                    full_id: j,
                    representative_id: rep,
                    operation_id: op_id,
                    reciprocal_shift: std::array::from_fn(|a| {
                        (p.fractional[a] - transformed[a]).round() as i32
                    }),
                    multiplicity: 0,
                });
                count += 1;
            }
        }
        if count == 0 {
            return Err(SymmetryError("operation group lacks identity"));
        }
        for e in entries
            .iter_mut()
            .flatten()
            .filter(|e| e.representative_id == rep)
        {
            e.multiplicity = count;
        }
        points.push(KPoint {
            fractional: full[i].fractional,
            weight: count as f64 / full.len() as f64,
        });
    }
    kpoints::validate(&points).map_err(|_| SymmetryError("reduced weight sum"))?;
    Ok(KOrbitMap {
        points,
        operations: subgroup,
        entries: entries.into_iter().map(Option::unwrap).collect(),
        phase_convention: "exp(-i(k+G).tau), physical reciprocal vectors; fractional dot includes 2*pi",
    })
}

/// Real-space quadrature must be permuted by x -> R*x+tau. Fourier phases
/// alone do not make a nonlinear XC quadrature invariant on an arbitrary grid.
pub fn preserves_fft_grid(dims: [usize; 3], operations: &[SymmetryOperation]) -> bool {
    !dims.contains(&0)
        && operations.iter().all(|o| {
            (0..3).all(|i| {
                let t = dims[i] as f64 * o.translation[i];
                (t - t.round()).abs() <= 1e-8
                    && (0..3).all(|j| {
                        (dims[i] as i128 * o.rotation[i][j] as i128) % dims[j] as i128 == 0
                    })
            })
        })
}
/// Choose smooth dimensions at least as large as the alias-free density grid.
/// First preserve anisotropic dimensions; if rotations mix incompatible axes,
/// search a common dimension. Arbitrary incommensurate origins fail explicitly.
pub fn compatible_fft_grid(
    required: [usize; 3],
    operations: &[SymmetryOperation],
) -> Result<[usize; 3], SymmetryError> {
    fn smooth(mut n: usize) -> bool {
        for p in [2, 3, 5] {
            while n > 1 && n.is_multiple_of(p) {
                n /= p;
            }
        }
        n == 1
    }
    let mut dims = required;
    for i in 0..3 {
        dims[i]=(required[i]..=512).find(|&n| smooth(n) && operations.iter().all(|o| {let t=n as f64*o.translation[i];(t-t.round()).abs()<=1e-8})).ok_or(SymmetryError("no commensurate FFT dimension <=512; choose a symmetry-compatible origin or use symmetry off"))?;
    }
    if preserves_fft_grid(dims, operations) {
        return Ok(dims);
    }
    (*dims.iter().max().unwrap()..=512)
        .find(|&n| smooth(n) && preserves_fft_grid([n; 3], operations))
        .map(|n| [n; 3])
        .ok_or(SymmetryError("no symmetry-compatible FFT grid <=512"))
}
