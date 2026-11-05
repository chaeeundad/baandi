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
