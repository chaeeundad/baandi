use dft_geometry::{Atom, Cell, kpoints::monkhorst_pack, symmetry::*};
#[test]
fn diamond_counts_and_full_mesh_reconstruction() {
    let c = Cell::from_user_rows_bohr([[0., 5., 5.], [5., 0., 5.], [5., 5., 0.]]).unwrap();
    let a = [Atom::new("Si", [0.; 3]), Atom::new("Si", [0.25; 3])];
    let ops = detect(&c, &a, 1e-7).unwrap();
    assert_eq!(ops.len(), 48);
    for shift in [0., 0.5] {
        let full = monkhorst_pack([4; 3], [shift; 3]).unwrap();
        let map = reduce_mesh(&full, &ops, true, false).unwrap();
        println!(
            "diamond shift={shift}: operations={} reduced={} weight_error={:e}",
            map.operations.len() / 2,
            map.points.len(),
            (map.points.iter().map(|p| p.weight).sum::<f64>() - 1.).abs()
        );
        assert_eq!(map.points.len(), if shift == 0. { 8 } else { 10 });
        for e in &map.entries {
            let op = &map.operations[e.operation_id];
            // Independent reconstruction: R^T*(k_full-shift) = +/- k_rep.
            for i in 0..3 {
                let lhs: f64 = (0..3)
                    .map(|j| {
                        op.rotation[j][i] as f64
                            * (full[e.full_id].fractional[j] - e.reciprocal_shift[j] as f64)
                    })
                    .sum();
                let rhs = map.points[e.representative_id].fractional[i]
                    * if op.antiunitary { -1. } else { 1. };
                assert!((lhs - rhs).abs() < 1e-14);
            }
            assert_eq!(
                map.entries
                    .iter()
                    .filter(|q| q.representative_id == e.representative_id)
                    .count(),
                e.multiplicity
            );
        }
        if shift == 0.5 {
            assert!(reduce_mesh(&full, &ops, true, true).is_err());
        }
    }
}
#[test]
fn cubic_slab_skew_species_and_origin() {
    for rows in [
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]],
        [[5., 0., 0.], [0., 5., 0.], [0., 0., 20.]],
        [[5., 0.3, 0.1], [0.4, 6., 0.2], [0.6, 0.7, 7.]],
    ] {
        let c = Cell::from_user_rows_bohr(rows).unwrap();
        let ops = detect(&c, &[Atom::new("X", [0.13, 0.23, 0.37])], 1e-7).unwrap();
        assert_eq!(
            ops.len(),
            if rows[0][1] != 0. {
                2
            } else if rows[2][2] == 20. {
                16
            } else {
                48
            }
        );
        validate_group(&c, &ops, 1e-7).unwrap();
        let a = [
            Atom::new("A", [0.; 3]),
            Atom::new("B", [0.173, 0.237, 0.319]),
        ];
        assert_eq!(detect(&c, &a, 1e-7).unwrap().len(), 1);
    }
}
#[test]
fn ambiguous_atoms_and_invalid_operations_reject() {
    let c = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    assert!(
        detect(
            &c,
            &[Atom::new("X", [0.; 3]), Atom::new("X", [1e-9; 3])],
            1e-7
        )
        .is_err()
    );
    let mut ops = detect(&c, &[Atom::new("X", [0.; 3])], 1e-7).unwrap();
    ops[0].spin_action = SpinAction::Swap;
    assert!(reduce_mesh(&monkhorst_pack([2; 3], [0.; 3]).unwrap(), &ops, true, false).is_err());
}
#[test]
fn bloch_phase_reconstruction_and_fft_quadrature() {
    let cell = Cell::from_user_rows_bohr([[0., 5., 5.], [5., 0., 5.], [5., 5., 0.]]).unwrap();
    let ops = detect(
        &cell,
        &[Atom::new("Si", [0.; 3]), Atom::new("Si", [0.25; 3])],
        1e-7,
    )
    .unwrap();
    assert!(!preserves_fft_grid([25; 3], &ops));
    assert!(preserves_fft_grid([32; 3], &ops));
    let full = monkhorst_pack([4; 3], [0.5; 3]).unwrap();
    let map = reduce_mesh(&full, &ops, true, false).unwrap();
    for entry in &map.entries {
        let op = &map.operations[entry.operation_id];
        let inverse = ops
            .iter()
            .find(|b| {
                (0..3).all(|i| {
                    (0..3).all(|j| {
                        (0..3)
                            .map(|k| op.rotation[i][k] * b.rotation[k][j])
                            .sum::<i32>()
                            == i32::from(i == j)
                    })
                })
            })
            .unwrap();
        for g in [[0, 0, 0], [1, -2, 3], [-2, 1, -1]] {
            let q = std::array::from_fn(|i| {
                map.points[entry.representative_id].fractional[i] + g[i] as f64
            });
            let transformed = op.reciprocal_action(q);
            let x = [0.173, 0.271, 0.413];
            // Direct wave at R^-1(x-tau), conjugated for an antiunitary op.
            let preimage: [f64; 3] = std::array::from_fn(|i| {
                (0..3)
                    .map(|j| inverse.rotation[i][j] as f64 * (x[j] - op.translation[j]))
                    .sum()
            });
            let direct = (0..3).map(|i| q[i] * preimage[i]).sum::<f64>()
                * if op.antiunitary { -1. } else { 1. };
            let phase = -(0..3)
                .map(|i| transformed[i] * op.translation[i])
                .sum::<f64>();
            let rebuilt = phase + (0..3).map(|i| transformed[i] * x[i]).sum::<f64>();
            let tau = std::f64::consts::TAU;
            assert!(((tau * direct).cos() - (tau * rebuilt).cos()).abs() < 1e-13);
            assert!(((tau * direct).sin() - (tau * rebuilt).sin()).abs() < 1e-13);
        }
    }
    let cubic = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let edge = 1e-7 * (1. - 1e-8) / 5.;
    assert!(
        detect(
            &cubic,
            &[Atom::new("X", [0.; 3]), Atom::new("X", [edge, 0., 0.])],
            1e-7
        )
        .is_err()
    );
}

// ADR 0049: default tolerance and input symmetrization.
const ANGSTROM_TO_BOHR: f64 = 1.0 / 0.529_177_210_903;
/// Materials Project style hcp record: lattice in Angstrom and fractional
/// coordinates rounded to six decimals (Ru-like a=2.714, c=4.282 Angstrom).
fn hcp_six_decimals(a: f64, c: f64) -> (Cell, Vec<Atom>) {
    let r6 = |v: f64| (v * 1e6).round() / 1e6;
    let h = r6(a * 3f64.sqrt() / 2.);
    let rows = [[r6(a / 2.), -h, 0.], [r6(a / 2.), h, 0.], [0., 0., r6(c)]]
        .map(|v| v.map(|x| x * ANGSTROM_TO_BOHR));
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
    let atoms = vec![
        Atom::new("Ru", [0.333333, 0.666667, 0.25]),
        Atom::new("Ru", [0.666667, 0.333333, 0.75]),
    ];
    (cell, atoms)
}
fn exact_hcp(a: f64, c: f64) -> (Cell, Vec<Atom>) {
    let h = a * 3f64.sqrt() / 2.;
    let rows =
        [[a / 2., -h, 0.], [a / 2., h, 0.], [0., 0., c]].map(|v| v.map(|x| x * ANGSTROM_TO_BOHR));
    let cell = Cell::from_user_rows_bohr(rows).unwrap();
    let atoms = vec![
        Atom::new("Ru", [1. / 3., 2. / 3., 0.25]),
        Atom::new("Ru", [2. / 3., 1. / 3., 0.75]),
    ];
    (cell, atoms)
}
fn irreducible(cell: &Cell, atoms: &[Atom], tol: f64, grid: [usize; 3]) -> usize {
    let ops = detect(cell, atoms, tol).unwrap();
    let full = monkhorst_pack(grid, [0.; 3]).unwrap();
    reduce_mesh(&full, &ops, true, false).unwrap().points.len()
}

#[test]
fn hcp_six_decimal_coordinates_find_full_group_and_symmetrize() {
    for (a, c) in [(2.714, 4.282), (3.209, 5.211)] {
        // Ru-like and Mg-like c/a.
        let (cell, atoms) = hcp_six_decimals(a, c);
        // The old 1e-7 Bohr default misses most of D6h.
        let strict = detect(&cell, &atoms, 1e-7).unwrap().len();
        assert!(strict < 24, "strict tolerance found {strict}");
        let ops = detect(&cell, &atoms, DEFAULT_TOLERANCE_BOHR).unwrap();
        assert_eq!(ops.len(), 24); // P6_3/mmc, includes inversion.
        assert!(
            ops.iter()
                .any(|o| o.rotation == [[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
        );
        let s = symmetrize(&cell, &atoms, DEFAULT_TOLERANCE_BOHR).unwrap();
        assert_eq!(s.operations_count, 24);
        assert!(s.positions_changed && s.translations_snapped);
        assert!(s.max_atom_displacement_bohr > 1e-7 && s.max_atom_displacement_bohr < 1e-5);
        assert!(s.max_lattice_change_bohr < 1e-5);
        // Rounding errors cancel by pairs here, so the projection lands on the
        // exact Wyckoff 2c/2d sites to roundoff (fractional, dimensionless).
        let exact = exact_hcp(a, c).1;
        for (p, e) in s.fractional.iter().zip(&exact) {
            for (x, y) in p.iter().zip(e.fractional_position) {
                assert!((x - y).abs() < 1e-14);
            }
        }
        let new_cell = Cell::from_matrix_columns_bohr(s.lattice_columns_bohr).unwrap();
        let new_atoms: Vec<_> = s.fractional.iter().map(|p| Atom::new("Ru", *p)).collect();
        // Exact to roundoff: the tight tolerance finds the same group.
        let tight = detect(&new_cell, &new_atoms, 1e-10).unwrap();
        assert_eq!(tight.len(), 24);
        // Translations are commensurate with an even FFT grid along c.
        assert!(compatible_fft_grid([24, 24, 40], &tight).is_ok());
        // Same irreducible k set as the exact structure on hcp-compatible meshes.
        let (ec, ea) = exact_hcp(a, c);
        for grid in [[6, 6, 4], [12, 12, 8], [9, 9, 6]] {
            let reduced = irreducible(&new_cell, &new_atoms, DEFAULT_TOLERANCE_BOHR, grid);
            assert_eq!(reduced, irreducible(&ec, &ea, 1e-7, grid));
            println!(
                "hcp a={a} c={c} grid={grid:?}: tol 1e-7 -> {} k, default -> {reduced} k",
                irreducible(&cell, &atoms, 1e-7, grid)
            );
        }
    }
}

#[test]
fn exact_inputs_are_left_bit_identical() {
    let diamond = Cell::from_user_rows_bohr([[0., 5., 5.], [5., 0., 5.], [5., 5., 0.]]).unwrap();
    let si = [Atom::new("Si", [0.; 3]), Atom::new("Si", [0.25; 3])];
    let (hcp, ru) = exact_hcp(2.714, 4.282);
    let skew = Cell::from_user_rows_bohr([[5., 0.3, 0.1], [0.4, 6., 0.2], [0.6, 0.7, 7.]]).unwrap();
    let low = [
        Atom::new("X", [0.13, 0.23, 0.37]),
        Atom::new("Y", [0.41, 0.07, 0.83]),
    ];
    for (cell, atoms) in [(&diamond, &si[..]), (&hcp, &ru[..]), (&skew, &low[..])] {
        let s = symmetrize(cell, atoms, DEFAULT_TOLERANCE_BOHR).unwrap();
        assert!(!s.positions_changed && !s.lattice_changed, "{s:?}");
        assert!(s.max_atom_displacement_bohr <= SYMMETRIZATION_ROUNDOFF_BOHR);
        assert_eq!(s.lattice_columns_bohr, cell.lattice_columns_bohr());
        for (p, a) in s.fractional.iter().zip(atoms) {
            assert_eq!(p.map(f64::to_bits), a.fractional_position.map(f64::to_bits));
        }
        assert_eq!(
            s.operations_count,
            detect(cell, atoms, 1e-7).unwrap().len(),
            "exact inputs: same group as the old default"
        );
    }
}

#[test]
fn distortion_above_tolerance_keeps_lower_group() {
    let (cell, mut atoms) = exact_hcp(2.714, 4.282);
    // 1e-3 in fractional z = 8.1e-3 Bohr, far above the 1e-4 Bohr default.
    atoms[1].fractional_position[2] += 1e-3;
    let low = detect(&cell, &atoms, DEFAULT_TOLERANCE_BOHR).unwrap().len();
    assert!(low < 24);
    let s = symmetrize(&cell, &atoms, DEFAULT_TOLERANCE_BOHR).unwrap();
    assert_eq!(s.operations_count, low);
    // The projection only removes the in-tolerance noise of the lower group: the
    // physical distortion stays (fractional z difference, dimensionless).
    assert!(((s.fractional[1][2] - s.fractional[0][2]) - (0.5 + 1e-3)).abs() < 1e-12);
    // Just above the tolerance (2e-4 Bohr) also keeps the lower group.
    let (cell, mut atoms) = exact_hcp(2.714, 4.282);
    atoms[1].fractional_position[2] += 2e-4 / (4.282 * ANGSTROM_TO_BOHR);
    assert!(detect(&cell, &atoms, DEFAULT_TOLERANCE_BOHR).unwrap().len() < 24);
}

#[test]
fn origin_shifts_snap_or_stay_consistent() {
    let (cell, exact) = exact_hcp(2.714, 4.282);
    // Within tolerance: a uniform 3e-6 shift is removed by snapping translations.
    let small: Vec<_> = exact
        .iter()
        .map(|a| Atom::new("Ru", a.fractional_position.map(|x| x + 3e-6)))
        .collect();
    let s = symmetrize(&cell, &small, DEFAULT_TOLERANCE_BOHR).unwrap();
    assert!(s.translations_snapped && s.operations_count == 24);
    for (p, e) in s.fractional.iter().zip(&exact) {
        for (x, y) in p.iter().zip(e.fractional_position) {
            assert!((x - y).abs() < 1e-14);
        }
    }
    // An arbitrary origin is not snapped; the group stays exact and the
    // structure is untouched (no change beyond roundoff).
    let far: Vec<_> = exact
        .iter()
        .map(|a| Atom::new("Ru", a.fractional_position.map(|x| x + 0.0123)))
        .collect();
    let s = symmetrize(&cell, &far, DEFAULT_TOLERANCE_BOHR).unwrap();
    assert_eq!(s.operations_count, 24);
    assert!(!s.positions_changed);
    // A user-chosen loose tolerance absorbs a real 8.1e-3 Bohr distortion; the
    // applied displacement is reported and stays within that tolerance.
    let mut distorted = exact.clone();
    distorted[1].fractional_position[2] += 1e-3;
    let s = symmetrize(&cell, &distorted, 2e-2).unwrap();
    assert_eq!(s.operations_count, 24);
    assert!(s.max_atom_displacement_bohr > 1e-3 && s.max_atom_displacement_bohr <= 2e-2);
}

#[test]
#[ignore = "diagnostic scan of hcp k counts; prints only"]
fn hcp_k_count_scan() {
    let (cell, atoms) = hcp_six_decimals(2.714, 4.282);
    println!(
        "ops at 1e-7: {}",
        detect(&cell, &atoms, 1e-7).unwrap().len()
    );
    for n in 6..=16 {
        for m in 3..=10 {
            let old = irreducible(&cell, &atoms, 1e-7, [n, n, m]);
            let new = irreducible(&cell, &atoms, DEFAULT_TOLERANCE_BOHR, [n, n, m]);
            println!("grid {n}x{n}x{m}: {old} -> {new}");
        }
    }
}

// GS1C: an explicit FFT grid keeps the largest subgroup that permutes it.
#[test]
fn explicit_grid_keeps_grid_preserving_subgroup() {
    let (cell, atoms) = hcp_six_decimals(2.714, 4.282);
    // An odd grid along c cannot carry the 6_3 screw translation (0,0,1/2).
    let odd = [24, 24, 25];
    let s = symmetrize_for_grid(&cell, &atoms, DEFAULT_TOLERANCE_BOHR, Some(odd)).unwrap();
    assert_eq!(s.detected_operations_count, 24);
    assert!(s.operations_count < 24 && s.operations_count > 1);
    assert_eq!(s.dropped_for_grid, 24 - s.operations_count);
    assert!(s.max_atom_displacement_bohr <= DEFAULT_TOLERANCE_BOHR);
    let new_cell = Cell::from_matrix_columns_bohr(s.lattice_columns_bohr).unwrap();
    let new_atoms: Vec<_> = s.fractional.iter().map(|p| Atom::new("Ru", *p)).collect();
    let kept = grid_subgroup(
        &new_cell,
        &detect(&new_cell, &new_atoms, 1e-9).unwrap(),
        odd,
        1e-9,
    );
    assert_eq!(kept.len(), s.operations_count);
    validate_group(&new_cell, &kept, 1e-9).unwrap();
    assert!(preserves_fft_grid(odd, &kept));
    assert!(
        kept.iter()
            .all(|o| o.translation[2].abs() < 1e-12 || (o.translation[2] - 1.).abs() < 1e-12)
    );
    // The subgroup still reduces a compatible k mesh.
    let full = monkhorst_pack([6, 6, 4], [0.; 3]).unwrap();
    let reduced = reduce_mesh(&full, &kept, true, false).unwrap();
    println!(
        "hcp odd grid {odd:?}: kept {} of 24 operations, 6x6x4 -> {} k",
        kept.len(),
        reduced.points.len()
    );
    // A grid that cannot carry the 3-fold axis (n1 != n2) keeps fewer still,
    // and never less than the identity.
    let skew =
        symmetrize_for_grid(&cell, &atoms, DEFAULT_TOLERANCE_BOHR, Some([24, 25, 25])).unwrap();
    assert!(skew.operations_count >= 1 && skew.operations_count < s.operations_count);
}

#[test]
fn compatible_grids_and_exact_inputs_are_bit_identical() {
    for (cell, atoms) in [exact_hcp(2.714, 4.282), hcp_six_decimals(2.714, 4.282)] {
        let free = symmetrize(&cell, &atoms, DEFAULT_TOLERANCE_BOHR).unwrap();
        let gridded =
            symmetrize_for_grid(&cell, &atoms, DEFAULT_TOLERANCE_BOHR, Some([24, 24, 40])).unwrap();
        assert_eq!(gridded.dropped_for_grid, 0);
        assert_eq!(free, gridded);
        for (a, b) in free.fractional.iter().zip(&gridded.fractional) {
            assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
        }
    }
    // Exact input with an incompatible grid: projected onto the subgroup it is
    // already exact, so the input bits are kept.
    let (cell, atoms) = exact_hcp(2.714, 4.282);
    let s = symmetrize_for_grid(&cell, &atoms, DEFAULT_TOLERANCE_BOHR, Some([24, 24, 25])).unwrap();
    assert!(s.dropped_for_grid > 0 && !s.positions_changed && !s.lattice_changed);
}

// GS1D: Cl2-like Cmca (No. 64) in a primitive cell with unequal first two
// vectors, a1 = b_conv, a2 = (a+b)/2, a3 = c (Angstrom a=6.145, b=4.395,
// c=8.154; Cl on 8f (0,y,z), y=0.1015, z=0.0464), six-decimal coordinates.
fn cl2_like() -> (Cell, Vec<Atom>) {
    let (a, b, c) = (
        6.145 * ANGSTROM_TO_BOHR,
        4.395 * ANGSTROM_TO_BOHR,
        8.154 * ANGSTROM_TO_BOHR,
    );
    let cell = Cell::from_user_rows_bohr([[0., b, 0.], [a / 2., b / 2., 0.], [0., 0., c]]).unwrap();
    let (y, z) = (0.1015, 0.0464);
    let r6 = |v: f64| (v * 1e6).round() / 1e6;
    let atoms = [
        [0., y, z],
        [0., 0.5 - y, z + 0.5],
        [0., y + 0.5, 0.5 - z],
        [0., -y, -z],
    ]
    .iter()
    .map(|p| {
        let f = cell.cartesian_bohr_to_fractional([p[0] * a, p[1] * b, p[2] * c]);
        Atom::new("Cl", f.map(|v| r6(v.rem_euclid(1.0))))
    })
    .collect();
    (cell, atoms)
}

#[test]
fn cl2_like_explicit_grid_drops_only_operations_that_do_not_map_the_grid() {
    let (cell, atoms) = cl2_like();
    let ops = detect(&cell, &atoms, DEFAULT_TOLERANCE_BOHR).unwrap();
    assert_eq!(ops.len(), 8);
    let s = symmetrize(&cell, &atoms, DEFAULT_TOLERANCE_BOHR).unwrap();
    let exact_cell = Cell::from_matrix_columns_bohr(s.lattice_columns_bohr).unwrap();
    let exact_atoms: Vec<_> = s.fractional.iter().map(|p| Atom::new("Cl", *p)).collect();
    let ops = detect(&exact_cell, &exact_atoms, 1e-9).unwrap();
    let translated = |o: &SymmetryOperation| o.translation.iter().any(|t| t.abs() > 1e-12);
    assert_eq!(ops.iter().filter(|o| translated(o)).count(), 4);
    let dims = [40, 36, 80];
    for o in &ops {
        let rotation_ok = grid_rotation_compatible(dims, o);
        let translation_ok = grid_translation_compatible(dims, o);
        // Independent check: map every unit grid step and test grid membership.
        let maps_grid = (0..3).all(|j| {
            (0..3).all(|i| {
                let image = o.rotation[i][j] as f64 / dims[j] as f64 * dims[i] as f64;
                (image - image.round()).abs() < 1e-12
            })
        });
        assert_eq!(rotation_ok, maps_grid, "{o:?}");
        // Every Cmca translation here is a half: commensurate with even dims.
        assert!(translation_ok, "{o:?}");
        println!(
            "rotation {:?} translation {:?}: rotation maps grid {rotation_ok}",
            o.rotation, o.translation
        );
    }
    let kept = grid_subgroup(&exact_cell, &ops, dims, 1e-9);
    assert_eq!(kept.len(), 4);
    let reasons = grid_drop_reasons(&ops, dims);
    assert_eq!((reasons.rotation, reasons.translation), (4, 0));
    // The dropped rotations mix a1 and a2 (x -> -x maps a2 to a1 - a2): a grid
    // step 1/36 along a2 becomes 1/36 along a1, not a multiple of 1/40. Any
    // grid with n2 dividing n1 keeps all eight.
    for dims in [[40, 40, 80], [72, 36, 80]] {
        assert_eq!(grid_subgroup(&exact_cell, &ops, dims, 1e-9).len(), 8);
    }
    let gridded =
        symmetrize_for_grid(&cell, &atoms, DEFAULT_TOLERANCE_BOHR, Some([40, 36, 80])).unwrap();
    assert_eq!((gridded.operations_count, gridded.dropped_for_grid), (4, 4));
    assert_eq!(gridded.dropped_for_grid_rotation, 4);
}

/// GT1 (ADR 0050): time reversal adds the k ~ -k partners. Without inversion
/// (regular tetrahedron, -43m) it shrinks the irreducible set; with inversion
/// in the mesh-preserving group the reduced points and weights are identical
/// with and without it, so such runs are unchanged bit for bit.
#[test]
fn time_reversal_folds_only_without_inversion() {
    let c = Cell::from_user_rows_bohr([[5., 0., 0.], [0., 5., 0.], [0., 0., 5.]]).unwrap();
    let a = 0.1;
    let tetra: Vec<_> = [[a, a, a], [a, -a, -a], [-a, a, -a], [-a, -a, a]]
        .into_iter()
        .map(|p| Atom::new("H", p))
        .collect();
    let ops = detect(&c, &tetra, 1e-7).unwrap();
    assert_eq!(ops.len(), 24);
    assert!(
        !ops.iter()
            .any(|o| o.rotation == [[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
    );
    let full = monkhorst_pack([3; 3], [0.; 3]).unwrap();
    let identity: Vec<_> = ops
        .iter()
        .filter(|o| o.rotation == [[1, 0, 0], [0, 1, 0], [0, 0, 1]])
        .cloned()
        .collect();
    assert_eq!(identity.len(), 1);
    let count =
        |ops: &[SymmetryOperation], tr| reduce_mesh(&full, ops, tr, false).unwrap().points.len();
    let (none, tr_only, crystal, crystal_tr) = (
        count(&identity, false),
        count(&identity, true),
        count(&ops, false),
        count(&ops, true),
    );
    println!(
        "tetrahedron 3x3x3: full={none} tr={tr_only} crystal={crystal} crystal+tr={crystal_tr}"
    );
    assert_eq!((none, tr_only), (27, 14));
    assert!(crystal_tr < crystal, "{crystal_tr} {crystal}");
    // Centrosymmetric: diamond 4x4x4 shifted and unshifted, hcp-like skew mesh.
    let d = Cell::from_user_rows_bohr([[0., 5., 5.], [5., 0., 5.], [5., 5., 0.]]).unwrap();
    let si = [Atom::new("Si", [0.; 3]), Atom::new("Si", [0.25; 3])];
    let dops = detect(&d, &si, 1e-7).unwrap();
    for (grid, shift) in [([4; 3], 0.), ([4; 3], 0.5), ([3, 4, 5], 0.)] {
        let full = monkhorst_pack(grid, [shift; 3]).unwrap();
        let with = reduce_mesh(&full, &dops, true, false).unwrap();
        let without = reduce_mesh(&full, &dops, false, false).unwrap();
        assert!(
            without
                .operations
                .iter()
                .any(|o| o.rotation == [[-1, 0, 0], [0, -1, 0], [0, 0, -1]])
        );
        assert_eq!(with.points.len(), without.points.len(), "{grid:?} {shift}");
        for (p, q) in with.points.iter().zip(&without.points) {
            assert_eq!(
                p.fractional.map(f64::to_bits),
                q.fractional.map(f64::to_bits)
            );
            assert_eq!(p.weight.to_bits(), q.weight.to_bits());
        }
    }
}

/// GT1 cause record (dev08 MP campaign, Ta6Zn7 mp-aaackwro, 6x6x3 Gamma
/// mesh): QE 7.6 reports 12 operations with inversion and 26 k; BAANDI 34.
/// Time reversal is on and adds nothing (inversion). The difference is the
/// mesh-preserving subgroup: in this cell setting only 4 of the 12
/// operations map the 6x6x3 mesh onto itself, and reduce_mesh uses only
/// those, while a point-wise reduction with all 12 (QE-like) gives 26.
#[test]
fn ta6zn7_mesh_subgroup_explains_k_count() {
    let c = Cell::from_user_rows_bohr([
        [-4.745669651409705, -8.219736886585833, 0.0],
        [4.745669651409705, -8.219736886585833, 0.0],
        [0.0, -5.4798258697719, 17.22086864710133],
    ])
    .unwrap();
    let sites: [(&str, [f64; 3]); 13] = [
        ("Ta", [0.16470429, 0.16470429, 0.50588812]),
        ("Ta", [0.83529571, 0.83529571, 0.49411188]),
        ("Ta", [0.35067452, 0.35067452, 0.94797944]),
        ("Ta", [0.64932548, 0.64932548, 0.05202056]),
        ("Ta", [0.45447874, 0.45447874, 0.63656577]),
        ("Ta", [0.54552126, 0.54552126, 0.36343423]),
        ("Zn", [0.0, 0.0, 0.0]),
        ("Zn", [0.58629293, 0.09410151, 0.22550305]),
        ("Zn", [0.09410151, 0.09410151, 0.22550305]),
        ("Zn", [0.09410151, 0.58629293, 0.22550305]),
        ("Zn", [0.41370707, 0.90589849, 0.77449695]),
        ("Zn", [0.90589849, 0.90589849, 0.77449695]),
        ("Zn", [0.90589849, 0.41370707, 0.77449695]),
    ];
    let atoms: Vec<_> = sites.iter().map(|(s, p)| Atom::new(*s, *p)).collect();
    let ops = detect(&c, &atoms, 1e-4).unwrap();
    assert_eq!(ops.len(), 12);
    let full = monkhorst_pack([6, 6, 3], [0.; 3]).unwrap();
    let with = reduce_mesh(&full, &ops, true, false).unwrap();
    let without = reduce_mesh(&full, &ops, false, false).unwrap();
    let mesh_ops = with.operations.iter().filter(|o| !o.antiunitary).count();
    println!(
        "Ta6Zn7 6x6x3: crystal={} mesh-preserving={mesh_ops} irreducible tr={} no-tr={}",
        ops.len(),
        with.points.len(),
        without.points.len()
    );
    assert_eq!(mesh_ops, 4);
    assert_eq!((with.points.len(), without.points.len()), (34, 34));
    // Point-wise orbits under all 12 operations and k ~ -k (QE-like count).
    let key = |k: [f64; 3]| k.map(|x| ((x - x.floor()) * 36.).round() as i64 % 36);
    let index: std::collections::HashMap<_, _> = full
        .iter()
        .enumerate()
        .map(|(i, p)| (key(p.fractional), i))
        .collect();
    let mut parent: Vec<usize> = (0..full.len()).collect();
    fn root(p: &mut [usize], mut i: usize) -> usize {
        while p[i] != i {
            i = p[i];
        }
        i
    }
    for (i, p) in full.iter().enumerate() {
        for op in &ops {
            for sign in [1., -1.] {
                let k = op.reciprocal_action(p.fractional).map(|x| sign * x);
                if let Some(&j) = index.get(&key(k)) {
                    let (a, b) = (root(&mut parent, i), root(&mut parent, j));
                    parent[a] = b;
                }
            }
        }
    }
    let orbits = (0..full.len())
        .filter(|&i| root(&mut parent, i) == i)
        .count();
    assert_eq!(orbits, 26);
}
