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

/// Default Cartesian atom/lattice matching tolerance, Bohr (ADR 0049).
/// About 5.3e-5 Angstrom: between VASP's SYMPREC default (1e-5 Angstrom) and a
/// 1e-5 crystal-coordinate acceptance on a 10 Bohr cell. It accepts structures
/// whose fractional coordinates are rounded to six decimals in cells up to
/// roughly 50 Bohr and is still far below any physically meaningful distortion.
pub const DEFAULT_TOLERANCE_BOHR: f64 = 1e-4;
/// Changes at or below this Cartesian size (Bohr) are treated as roundoff:
/// the input bits are kept, so exactly symmetric inputs stay bit-identical.
pub const SYMMETRIZATION_ROUNDOFF_BOHR: f64 = 1e-12;
/// Tight tolerance (Bohr) used to verify that a symmetrized structure carries
/// the detected group exactly (up to floating-point roundoff).
pub const SYMMETRIZATION_VERIFY_BOHR: f64 = 1e-9;
/// Largest denominator tried when snapping fractional translations to rationals.
const SNAP_MAX_DENOMINATOR: i64 = 24;

/// Input geometry projected onto the group detected at a tolerance (ADR 0049).
#[derive(Clone, Debug, PartialEq)]
pub struct Symmetrized {
    /// Lattice vectors as matrix columns [Cartesian][vector], Bohr.
    pub lattice_columns_bohr: [[f64; 3]; 3],
    /// Fractional positions in the input order (not wrapped).
    pub fractional: Vec<[f64; 3]>,
    /// Number of affine operations (rotation + translation) kept: the
    /// detected group, or its explicit-FFT-grid subgroup (GS1C).
    pub operations_count: usize,
    /// Number of affine operations detected at the tolerance.
    pub detected_operations_count: usize,
    /// Detected operations dropped because they do not permute the explicit
    /// FFT grid (GS1C); zero without an explicit grid.
    pub dropped_for_grid: usize,
    /// Of those, operations whose rotation does not map the grid (GS1D).
    pub dropped_for_grid_rotation: usize,
    /// Of those, operations with a commensurate rotation but a translation
    /// that is not a multiple of the grid spacing (GS1D).
    pub dropped_for_grid_translation: usize,
    /// Largest Cartesian atom displacement, Bohr, measured in the input lattice.
    pub max_atom_displacement_bohr: f64,
    /// Largest lattice-vector change |a_i' - a_i|, Bohr.
    pub max_lattice_change_bohr: f64,
    /// True when positions were replaced (displacement above roundoff).
    pub positions_changed: bool,
    /// True when the lattice was replaced (change above roundoff).
    pub lattice_changed: bool,
    /// True when translations were snapped to small-denominator rationals.
    pub translations_snapped: bool,
}

fn wrap(v: f64) -> f64 {
    v - v.round()
}

/// Snap one fractional component to the nearest p/q (q <= SNAP_MAX_DENOMINATOR,
/// smallest q first) within `bound`; None if no such rational is close enough.
fn snap_component(t: f64, bound: f64) -> Option<f64> {
    (1..=SNAP_MAX_DENOMINATOR).find_map(|q| {
        let p = (t * q as f64).round();
        let s = p / q as f64;
        ((t - s).abs() <= bound).then_some(s)
    })
}

/// Largest subgroup of `operations` whose members each permute the FFT grid
/// `dims` (rotation and fractional translation; GS1C, ADR 0049). The
/// grid-preserving affine maps are closed under composition and inversion, so
/// this subset of a group is a group; closure is re-checked and the
/// identity-only group (always grid-preserving) is the fallback.
pub fn grid_subgroup(
    cell: &Cell,
    operations: &[SymmetryOperation],
    dims: [usize; 3],
    tolerance_bohr: f64,
) -> Vec<SymmetryOperation> {
    grid_subgroup_with(cell, operations, dims, tolerance_bohr, false)
}
/// GS1D-B: with `rotations_free` an operation is kept when its fractional
/// translation is a multiple of the grid spacing, whatever its rotation does to
/// the grid (the caller then projects the local potential in G space).
/// Translation-compatible operations are closed under composition, so the kept
/// subset is again a group.
pub fn grid_subgroup_with(
    cell: &Cell,
    operations: &[SymmetryOperation],
    dims: [usize; 3],
    tolerance_bohr: f64,
    rotations_free: bool,
) -> Vec<SymmetryOperation> {
    let kept: Vec<_> = operations
        .iter()
        .filter(|o| grid_keeps(dims, o, rotations_free))
        .cloned()
        .collect();
    if !kept.is_empty() && validate_group(cell, &kept, tolerance_bohr).is_ok() {
        return kept;
    }
    operations
        .iter()
        .filter(|o| {
            o.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
                && o.translation.iter().all(|t| wrap(*t).abs() <= 1e-12)
        })
        .take(1)
        .cloned()
        .collect()
}

/// `symmetrize_for_grid` without an explicit FFT grid (full detected group).
pub fn symmetrize(
    cell: &Cell,
    atoms: &[Atom],
    tolerance_bohr: f64,
) -> Result<Symmetrized, SymmetryError> {
    symmetrize_for_grid(cell, atoms, tolerance_bohr, None)
}

/// Group-averaged (Reynolds) projection of atoms and lattice onto the group
/// detected at `tolerance_bohr`. Fractional translations are first snapped to
/// small-denominator rationals when every snapped translation stays within the
/// tolerance and the snapped set closes exactly; otherwise translations are the
/// least-squares estimate from all atoms (consistent with an origin shift).
/// With an explicit FFT grid (GS1C) only the operations that permute that grid
/// (with the chosen translations) are kept, and the structure is projected onto
/// that subgroup only. The lattice metric is averaged over the kept rotations
/// and the lattice vectors are updated by the symmetric stretch
/// A' = A M^-1/2 M'^1/2 (no rigid rotation). Refuses (error) when any atom or
/// lattice vector would move by more than the tolerance, or when the result
/// does not carry the kept group at a tight tolerance.
pub fn symmetrize_for_grid(
    cell: &Cell,
    atoms: &[Atom],
    tolerance_bohr: f64,
    grid: Option<[usize; 3]>,
) -> Result<Symmetrized, SymmetryError> {
    symmetrize_for_grid_with(cell, atoms, tolerance_bohr, grid, false)
}
/// `symmetrize_for_grid` with the GS1D-B policy (`grid_subgroup_with`).
pub fn symmetrize_for_grid_with(
    cell: &Cell,
    atoms: &[Atom],
    tolerance_bohr: f64,
    grid: Option<[usize; 3]>,
    rotations_free: bool,
) -> Result<Symmetrized, SymmetryError> {
    let detected = detect(cell, atoms, tolerance_bohr)?;
    let tol = tolerance_bohr;
    let a = cell.lattice;
    let x: Vec<Vector3<f64>> = atoms
        .iter()
        .map(|a| Vector3::from(a.fractional_position))
        .collect();
    // Least-squares translations from all atoms (origin-shift consistent).
    let estimated: Vec<Vector3<f64>> = detected
        .iter()
        .map(|o| {
            let r = matrix(o.rotation);
            let t = Vector3::from(o.translation);
            let mean = o
                .atom_permutation
                .iter()
                .enumerate()
                .fold(Vector3::zeros(), |s, (j, &i)| {
                    s + (x[i] - r * x[j] - t).map(wrap)
                })
                / x.len() as f64;
            t + mean
        })
        .collect();
    let bounds: [f64; 3] = std::array::from_fn(|i| tol * cell.inverse.row(i).norm());
    let snapped: Option<Vec<Vector3<f64>>> = estimated
        .iter()
        .map(|t| {
            let s = Vector3::new(
                snap_component(t[0], bounds[0])?,
                snap_component(t[1], bounds[1])?,
                snap_component(t[2], bounds[2])?,
            );
            ((a * (s - t)).norm() <= tol).then_some(s)
        })
        .collect();
    // Exact closure of a candidate operation set (translations included).
    let closes = |ops: &[SymmetryOperation]| {
        ops.iter().all(|oa| {
            ops.iter().all(|ob| {
                let r = matrix(oa.rotation) * matrix(ob.rotation);
                let t = matrix(oa.rotation) * Vector3::from(ob.translation)
                    + Vector3::from(oa.translation);
                ops.iter().any(|oc| {
                    matrix(oc.rotation) == r
                        && ob
                            .atom_permutation
                            .iter()
                            .enumerate()
                            .all(|(i, &j)| oa.atom_permutation[j] == oc.atom_permutation[i])
                        && (0..3).all(|k| wrap(oc.translation[k] - t[k]).abs() <= 1e-12)
                })
            })
        })
    };
    let with_taus = |taus: &[Vector3<f64>]| -> Vec<SymmetryOperation> {
        detected
            .iter()
            .zip(taus)
            .map(|(o, t)| SymmetryOperation {
                translation: (*t).into(),
                ..o.clone()
            })
            .collect()
    };
    let restrict = |ops: Vec<SymmetryOperation>| -> Vec<SymmetryOperation> {
        match grid {
            None => ops,
            Some(dims) => grid_subgroup_with(cell, &ops, dims, tol, rotations_free),
        }
    };
    let project = |ops: &[SymmetryOperation]| -> (Vec<Vector3<f64>>, f64) {
        let mut sums = vec![Vector3::zeros(); x.len()];
        for o in ops {
            let r = matrix(o.rotation);
            let t = Vector3::from(o.translation);
            for (j, &i) in o.atom_permutation.iter().enumerate() {
                sums[i] += (r * x[j] + t - x[i]).map(wrap);
            }
        }
        let n = ops.len() as f64;
        let delta: Vec<_> = sums.into_iter().map(|s| s / n).collect();
        let worst = delta.iter().map(|d| (a * d).norm()).fold(0.0, f64::max);
        (delta, worst)
    };
    let snapped_ops = snapped.map(|s| with_taus(&s)).filter(|ops| closes(ops));
    let mut candidates = Vec::new();
    if let Some(ops) = &snapped_ops {
        candidates.push((true, restrict(ops.clone())));
    }
    candidates.push((false, restrict(with_taus(&estimated))));
    let (translations_snapped, kept, delta, displacement) = candidates
        .into_iter()
        .map(|(snap, ops)| {
            let (d, w) = project(&ops);
            (snap, ops, d, w)
        })
        .find(|(_, _, _, w)| *w <= tol)
        .ok_or(SymmetryError(
            "atom symmetrization would exceed the tolerance",
        ))?;
    let dropped_reasons = match grid {
        Some(dims) if kept.len() < detected.len() => {
            // Reasons for the operations not kept, with the chosen translations.
            let chosen = if translations_snapped {
                snapped_ops.clone().unwrap_or_default()
            } else {
                with_taus(&estimated)
            };
            let dropped: Vec<_> = chosen
                .into_iter()
                .filter(|o| {
                    !kept.iter().any(|k| {
                        k.rotation == o.rotation
                            && (0..3)
                                .all(|i| wrap(k.translation[i] - o.translation[i]).abs() <= 1e-12)
                    })
                })
                .collect();
            if rotations_free {
                GridDropReasons {
                    rotation: 0,
                    translation: dropped.len(),
                }
            } else {
                grid_drop_reasons(&dropped, dims)
            }
        }
        _ => GridDropReasons::default(),
    };
    // Lattice: Reynolds average of the metric over the kept rotations.
    let metric = a.transpose() * a;
    let averaged = kept.iter().fold(Matrix3::zeros(), |s, o| {
        let r = matrix(o.rotation);
        s + r.transpose() * metric * r
    }) / kept.len() as f64;
    let averaged = (averaged + averaged.transpose()) * 0.5;
    let power = |m: Matrix3<f64>, p: f64| -> Option<Matrix3<f64>> {
        let e = m.symmetric_eigen();
        if e.eigenvalues.iter().any(|&v| v <= 0.0 || !v.is_finite()) {
            return None;
        }
        let d = Matrix3::from_diagonal(&e.eigenvalues.map(|v| v.powf(p)));
        Some(e.eigenvectors * d * e.eigenvectors.transpose())
    };
    let stretched = power(metric, -0.5)
        .zip(power(averaged, 0.5))
        .map(|(mi, ms)| a * mi * ms)
        .ok_or(SymmetryError("lattice metric is not positive definite"))?;
    let lattice_change = (0..3)
        .map(|j| (stretched.column(j) - a.column(j)).norm())
        .fold(0.0, f64::max);
    if lattice_change > tol {
        return Err(SymmetryError(
            "lattice symmetrization would exceed the tolerance",
        ));
    }
    let lattice_changed = lattice_change > SYMMETRIZATION_ROUNDOFF_BOHR;
    let new_lattice = if lattice_changed { stretched } else { a };
    let positions_changed = displacement > SYMMETRIZATION_ROUNDOFF_BOHR;
    let fractional: Vec<[f64; 3]> = if positions_changed {
        x.iter().zip(&delta).map(|(p, d)| (p + d).into()).collect()
    } else {
        atoms.iter().map(|a| a.fractional_position).collect()
    };
    let result = Symmetrized {
        lattice_columns_bohr: std::array::from_fn(|i| std::array::from_fn(|j| new_lattice[(i, j)])),
        fractional,
        operations_count: kept.len(),
        detected_operations_count: detected.len(),
        dropped_for_grid: detected.len() - kept.len(),
        dropped_for_grid_rotation: dropped_reasons.rotation,
        dropped_for_grid_translation: dropped_reasons.translation,
        max_atom_displacement_bohr: displacement,
        max_lattice_change_bohr: lattice_change,
        positions_changed,
        lattice_changed,
        translations_snapped: translations_snapped && positions_changed,
    };
    if positions_changed || lattice_changed {
        let new_cell = Cell::from_matrix_columns_bohr(result.lattice_columns_bohr)
            .map_err(|_| SymmetryError("symmetrized lattice is not a valid cell"))?;
        let new_atoms: Vec<Atom> = atoms
            .iter()
            .zip(&result.fractional)
            .map(|(a, p)| Atom::new(a.species.clone(), *p))
            .collect();
        let exact = detect(&new_cell, &new_atoms, SYMMETRIZATION_VERIFY_BOHR)?;
        let exact = match grid {
            None => exact,
            Some(dims) => grid_subgroup_with(
                &new_cell,
                &exact,
                dims,
                SYMMETRIZATION_VERIFY_BOHR,
                rotations_free,
            ),
        };
        if exact.len() != kept.len() {
            return Err(SymmetryError(
                "symmetrized structure does not carry the kept group exactly",
            ));
        }
    }
    Ok(result)
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
/// How a full mesh is reduced (GT2, ADR 0051).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum KReduction {
    /// Orbits under the subgroup that maps the whole mesh onto itself. The
    /// weighted irreducible sum equals the full-mesh sum exactly.
    MeshSubgroup,
    /// Point-wise orbits under every crystal operation: k and S k are merged
    /// whenever S k lies on the mesh (mod G). The weighted irreducible sum
    /// equals the full-mesh sum projected onto the crystal group.
    Pointwise,
    /// No reduction (every mesh point is its own orbit) with the crystal group
    /// as the operation list: the oracle of `Pointwise` (full mesh, crystal
    /// projection of the density, crystal symmetrization of derivatives).
    Unreduced,
}
impl KReduction {
    pub fn label(self) -> &'static str {
        match self {
            Self::MeshSubgroup => "mesh_subgroup",
            Self::Pointwise => "pointwise_crystal",
            Self::Unreduced => "unreduced_crystal",
        }
    }
}
#[derive(Clone, Debug)]
pub struct KOrbitMap {
    pub points: Vec<KPoint>,
    /// Operations of the reduction, including explicit antiunitary partners:
    /// the mesh-preserving subgroup (`MeshSubgroup`) or the whole crystal
    /// group (`Pointwise`, `Unreduced`). Their unitary part is the group the
    /// density and the derivatives are symmetrized with.
    pub operations: Vec<SymmetryOperation>,
    pub entries: Vec<KOrbitEntry>,
    pub phase_convention: &'static str,
    /// Reduction actually applied. `Pointwise` on a mesh that every operation
    /// (and time reversal, when used) preserves is the `MeshSubgroup` result
    /// bit for bit and is recorded as `MeshSubgroup`.
    pub reduction: KReduction,
    /// Unitary crystal operations that map the whole mesh onto itself.
    pub mesh_preserving_operations: usize,
}
fn equivalent(a: [f64; 3], b: [f64; 3]) -> bool {
    (0..3).all(|i| {
        let d = a[i] - b[i];
        (d - d.round()).abs() <= 1e-12
    })
}
/// Hash index of mesh points modulo reciprocal lattice vectors (GT2): replaces
/// the former all-pairs searches (N^2 |G| equivalence tests) by O(1) lookups.
/// Buckets are the wrapped coordinates rounded to 1e-6; a lookup scans the 27
/// neighbouring buckets and accepts with the same 1e-12 test as before, so the
/// result does not depend on bucket boundaries.
struct MeshIndex<'a> {
    points: &'a [KPoint],
    buckets: std::collections::HashMap<[i64; 3], Vec<usize>>,
}
const MESH_KEY_SCALE: f64 = 1e6;
impl<'a> MeshIndex<'a> {
    fn key(k: [f64; 3]) -> [i64; 3] {
        k.map(|x| ((x - x.floor()) * MESH_KEY_SCALE).round() as i64 % MESH_KEY_SCALE as i64)
    }
    fn new(points: &'a [KPoint]) -> Self {
        let mut buckets: std::collections::HashMap<[i64; 3], Vec<usize>> =
            std::collections::HashMap::with_capacity(points.len());
        for (i, p) in points.iter().enumerate() {
            buckets.entry(Self::key(p.fractional)).or_default().push(i);
        }
        Self { points, buckets }
    }
    /// Smallest index of a mesh point equivalent to k, if any.
    fn find(&self, k: [f64; 3]) -> Option<usize> {
        let base = Self::key(k);
        let n = MESH_KEY_SCALE as i64;
        let mut best: Option<usize> = None;
        for d0 in -1..=1 {
            for d1 in -1..=1 {
                for d2 in -1..=1 {
                    let key = [
                        (base[0] + d0).rem_euclid(n),
                        (base[1] + d1).rem_euclid(n),
                        (base[2] + d2).rem_euclid(n),
                    ];
                    for &j in self.buckets.get(&key).into_iter().flatten() {
                        if equivalent(k, self.points[j].fractional) && best.is_none_or(|b| j < b) {
                            best = Some(j);
                        }
                    }
                }
            }
        }
        best
    }
    fn preserved_by(&self, o: &SymmetryOperation) -> bool {
        self.points
            .iter()
            .all(|p| self.find(o.reciprocal_action(p.fractional)).is_some())
    }
}
/// Uniform, duplicate-free full mesh and unitary unimodular identity-spin
/// operations (shared input contract of both reductions).
fn check_mesh_input<'a>(
    full: &'a [KPoint],
    operations: &[SymmetryOperation],
) -> Result<MeshIndex<'a>, SymmetryError> {
    kpoints::validate(full).map_err(|_| SymmetryError("invalid full mesh"))?;
    let index = MeshIndex::new(full);
    if full
        .iter()
        .any(|p| (p.weight - 1.0 / full.len() as f64).abs() > 1e-14)
        || full
            .iter()
            .enumerate()
            .any(|(i, p)| index.find(p.fractional) != Some(i))
    {
        return Err(SymmetryError("uniform duplicate-free full mesh required"));
    }
    for op in operations {
        if op.spin_action != SpinAction::Identity
            || op.antiunitary
            || (matrix(op.rotation).determinant().abs() - 1.0).abs() > 1e-8
        {
            return Err(SymmetryError(
                "unitary unimodular identity-spin operations required",
            ));
        }
    }
    Ok(index)
}
/// Orbits of the mesh under `group` (images off the mesh ignored). The
/// representative is the first unassigned mesh index; each member records the
/// first operation of `group` mapping the representative onto it.
fn orbits(
    full: &[KPoint],
    index: &MeshIndex<'_>,
    group: &[SymmetryOperation],
) -> Result<(Vec<KPoint>, Vec<KOrbitEntry>), SymmetryError> {
    let mut points = Vec::new();
    let mut entries: Vec<Option<KOrbitEntry>> = vec![None; full.len()];
    let mut members: Vec<(usize, usize)> = Vec::with_capacity(group.len());
    for i in 0..full.len() {
        if entries[i].is_some() {
            continue;
        }
        let rep = points.len();
        members.clear();
        for (op_id, op) in group.iter().enumerate() {
            if let Some(j) = index.find(op.reciprocal_action(full[i].fractional))
                && !members.iter().any(|&(m, _)| m == j)
            {
                members.push((j, op_id));
            }
        }
        members.sort_unstable();
        if members.is_empty() {
            return Err(SymmetryError("operation group lacks identity"));
        }
        let count = members.len();
        for &(j, op_id) in &members {
            if entries[j].is_some() {
                return Err(SymmetryError(
                    "overlapping orbits: operation set is not a group",
                ));
            }
            let transformed = group[op_id].reciprocal_action(full[i].fractional);
            entries[j] = Some(KOrbitEntry {
                full_id: j,
                representative_id: rep,
                operation_id: op_id,
                reciprocal_shift: std::array::from_fn(|a| {
                    (full[j].fractional[a] - transformed[a]).round() as i32
                }),
                multiplicity: count,
            });
        }
        points.push(KPoint {
            fractional: full[i].fractional,
            weight: count as f64 / full.len() as f64,
        });
    }
    kpoints::validate(&points).map_err(|_| SymmetryError("reduced weight sum"))?;
    Ok((points, entries.into_iter().map(Option::unwrap).collect()))
}
/// Select the mesh-preserving subgroup, or reject any excluded operation when
/// strict=true. Only uniform, duplicate-free full meshes are accepted.
pub fn reduce_mesh(
    full: &[KPoint],
    operations: &[SymmetryOperation],
    time_reversal: bool,
    strict: bool,
) -> Result<KOrbitMap, SymmetryError> {
    let index = check_mesh_input(full, operations)?;
    let mut subgroup = Vec::new();
    for op in operations {
        if index.preserved_by(op) {
            subgroup.push(op.clone());
        } else if strict {
            return Err(SymmetryError("operation does not preserve shifted mesh"));
        }
    }
    if subgroup.is_empty() {
        return Err(SymmetryError("no mesh-preserving operations"));
    }
    let preserving = subgroup.len();
    if time_reversal {
        let partners: Vec<_> = subgroup
            .iter()
            .cloned()
            .map(|mut o| {
                o.antiunitary = true;
                o
            })
            .collect();
        if partners.iter().any(|o| !index.preserved_by(o)) {
            return Err(SymmetryError("time reversal does not preserve mesh"));
        }
        subgroup.extend(partners);
    }
    let (points, entries) = orbits(full, &index, &subgroup)?;
    Ok(KOrbitMap {
        points,
        operations: subgroup,
        entries,
        phase_convention: PHASE_CONVENTION,
        reduction: KReduction::MeshSubgroup,
        mesh_preserving_operations: preserving,
    })
}
const PHASE_CONVENTION: &str =
    "exp(-i(k+G).tau), physical reciprocal vectors; fractional dot includes 2*pi";

/// Reduce a full mesh with the selected policy (GT2, ADR 0051).
///
/// `Pointwise`: orbits are star(k) intersected with the mesh, star(k) = {S k}
/// over all crystal operations (and -S k with time reversal); images that
/// leave the mesh are ignored. Because the operations form a group, k ~ S k
/// is an equivalence on the mesh and each member of an orbit is the image of
/// its representative under one operation. The weight of an orbit is its size
/// over the mesh size. For a fixed crystal-invariant potential rho_{S k}(r) =
/// rho_k(S^-1 r), so with P the crystal projector
/// P(sum_mesh rho_k / N) = sum_orbits (|orbit| / N) P rho_rep: the reduced set
/// with the crystal projection reproduces the full mesh with the crystal
/// projection exactly. It equals the unprojected full-mesh sum only when the
/// mesh is invariant; otherwise the k quadrature becomes the group average of
/// the mesh (a k-sampling change, ADR 0051). When every operation (and time
/// reversal) preserves the mesh the result is `reduce_mesh` bit for bit.
pub fn reduce_mesh_with(
    full: &[KPoint],
    operations: &[SymmetryOperation],
    time_reversal: bool,
    reduction: KReduction,
) -> Result<KOrbitMap, SymmetryError> {
    if reduction == KReduction::MeshSubgroup {
        return reduce_mesh(full, operations, time_reversal, false);
    }
    let index = check_mesh_input(full, operations)?;
    let preserving = operations.iter().filter(|o| index.preserved_by(o)).count();
    if preserving == 0 {
        return Err(SymmetryError("no mesh-preserving operations"));
    }
    let partner = |o: &SymmetryOperation| {
        let mut o = o.clone();
        o.antiunitary = true;
        o
    };
    let tr_preserves = !time_reversal || operations.iter().all(|o| index.preserved_by(&partner(o)));
    if reduction == KReduction::Pointwise && preserving == operations.len() && tr_preserves {
        return reduce_mesh(full, operations, time_reversal, false);
    }
    let mut group = operations.to_vec();
    if reduction == KReduction::Pointwise && time_reversal {
        group.extend(operations.iter().map(partner));
    }
    let (points, entries) = if reduction == KReduction::Unreduced {
        let identity = group
            .iter()
            .position(|o| o.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
            .ok_or(SymmetryError("operation group lacks identity"))?;
        let entries: Vec<KOrbitEntry> = (0..full.len())
            .map(|j| KOrbitEntry {
                full_id: j,
                representative_id: j,
                operation_id: identity,
                reciprocal_shift: [0; 3],
                multiplicity: 1,
            })
            .collect();
        (full.to_vec(), entries)
    } else {
        orbits(full, &index, &group)?
    };
    Ok(KOrbitMap {
        points,
        operations: group,
        entries,
        phase_convention: PHASE_CONVENTION,
        reduction,
        mesh_preserving_operations: preserving,
    })
}

/// Rotation part of grid preservation: the grid point m_j/n_j e_j maps to
/// (R_ij m_j / n_j)_i, which lies on the grid for every m iff n_i R_ij / n_j
/// is an integer for all i, j. Necessary for any rotation that mixes axes of
/// unequal dimensions (GS1D).
pub fn grid_rotation_compatible(dims: [usize; 3], o: &SymmetryOperation) -> bool {
    !dims.contains(&0)
        && (0..3).all(|i| {
            (0..3).all(|j| (dims[i] as i128 * o.rotation[i][j] as i128) % dims[j] as i128 == 0)
        })
}
/// Translation part: each n_i tau_i is an integer (to 1e-8).
pub fn grid_translation_compatible(dims: [usize; 3], o: &SymmetryOperation) -> bool {
    (0..3).all(|i| {
        let t = dims[i] as f64 * o.translation[i];
        (t - t.round()).abs() <= 1e-8
    })
}
/// Number of operations that fail grid preservation, by first failing part
/// (rotation checked first).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct GridDropReasons {
    pub rotation: usize,
    pub translation: usize,
}
pub fn grid_drop_reasons(operations: &[SymmetryOperation], dims: [usize; 3]) -> GridDropReasons {
    let mut r = GridDropReasons::default();
    for o in operations {
        if !grid_rotation_compatible(dims, o) {
            r.rotation += 1;
        } else if !grid_translation_compatible(dims, o) {
            r.translation += 1;
        }
    }
    r
}
/// Whether `o` is kept for an explicit grid: it permutes the grid, or (GS1D-B,
/// `rotations_free`) at least its translation is a multiple of the spacing.
pub fn grid_keeps(dims: [usize; 3], o: &SymmetryOperation, rotations_free: bool) -> bool {
    !dims.contains(&0)
        && (rotations_free || grid_rotation_compatible(dims, o))
        && grid_translation_compatible(dims, o)
}
/// Real-space quadrature must be permuted by x -> R*x+tau. Fourier phases
/// alone do not make a nonlinear XC quadrature invariant on an arbitrary grid.
pub fn preserves_fft_grid(dims: [usize; 3], operations: &[SymmetryOperation]) -> bool {
    !dims.contains(&0)
        && operations
            .iter()
            .all(|o| grid_rotation_compatible(dims, o) && grid_translation_compatible(dims, o))
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
