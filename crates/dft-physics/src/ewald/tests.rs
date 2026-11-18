use super::*;

mod old_reference;
use old_reference::{old_energy, old_forces, old_stress};

/// Cubic, skewed, elongated and slab-like cells (rows in Bohr), scaled by `s`.
fn shapes(s: f64) -> [[[f64; 3]; 3]; 4] {
    [
        [[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]],
        [[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]],
        [[4., 0., 0.], [0., 5., 0.], [0., 0., 14.]],
        [[7., 0., 0.], [0., 8., 0.], [0.3, 0., 18.]],
    ]
    .map(|rows| rows.map(|r| r.map(|x| x * s)))
}

/// `atoms` ions with mixed charges (4, 6, 1, 3.5), deterministic spread-out
/// positions (some outside [0, 1)).
fn ions(atoms: usize) -> Vec<Ion> {
    (0..atoms)
        .map(|i| Ion {
            charge: [4.0, 6.0, 1.0, 3.5][i % 4],
            fractional_position: [
                ((i * 7 + 1) % 23) as f64 / 23. - 0.5 + 0.013 * (i / 23) as f64,
                ((i * 11 + 2) % 29) as f64 / 29. + 0.017 * (i / 29) as f64,
                ((i * 13 + 3) % 31) as f64 / 31. + 1.0 + 0.011 * (i / 31) as f64,
            ],
        })
        .collect()
}

fn cases() -> Vec<(Cell, Vec<Ion>)> {
    let mut out = Vec::new();
    for rows in shapes(1.0) {
        for atoms in [1usize, 2, 5, 9] {
            let mut list = ions(atoms);
            for (k, ion) in list.iter_mut().enumerate() {
                ion.charge = if k % 3 == 2 { 6.0 } else { 4.0 };
            }
            out.push((Cell::from_user_rows_bohr(rows).unwrap(), list));
        }
    }
    out
}

fn bits(e: &EwaldEnergy) -> [u64; 9] {
    [
        e.energy.0.to_bits(),
        e.real.0.to_bits(),
        e.reciprocal.0.to_bits(),
        e.self_energy.0.to_bits(),
        e.background.0.to_bits(),
        e.eta.to_bits(),
        e.real_cutoff_bohr.to_bits(),
        e.reciprocal_cutoff_inverse_bohr.to_bits(),
        e.tail_bound.0.to_bits(),
    ]
}

fn split(n: usize, parts: usize, k: usize) -> std::ops::Range<usize> {
    let first = (0..k).map(|q| n / parts + usize::from(q < n % parts)).sum();
    first..first + n / parts + usize::from(k < n % parts)
}

/// Every contiguous partition of the plan's candidates into 1 to 7 parts
/// reproduces `want` bit for bit, concatenated and chained.
fn assert_partitions(plan: &EwaldPlan<'_>, want: [u64; 9]) {
    for parts in 1..=7usize {
        let (mut real, mut reciprocal) = (Vec::new(), Vec::new());
        let (mut r, mut g) = (0.0, 0.0);
        for k in 0..parts {
            let (mut rk, mut gk) = (Vec::new(), Vec::new());
            plan.real_terms(split(plan.real_candidates(), parts, k), &mut rk)
                .unwrap();
            plan.reciprocal_terms(split(plan.reciprocal_candidates(), parts, k), &mut gk);
            r = accumulate(r, &rk);
            g = accumulate(g, &gk);
            real.extend(rk);
            reciprocal.extend(gk);
        }
        assert_eq!(
            bits(&plan.finish(&real, &reciprocal).unwrap()),
            want,
            "{parts}"
        );
        assert_eq!(bits(&plan.finish_sums(r, g).unwrap()), want, "{parts}");
    }
}

/// Reference scheme: the serial evaluator (pruned image enumeration) and every
/// contiguous partition of the plan's candidates into 1 to 7 parts reproduce
/// the eecd57f evaluator bit for bit (cubic, skewed, elongated, slab; 1 to 9
/// ions; mixed charges; positions outside [0, 1)); forces and stress too.
#[test]
fn reference_scheme_matches_the_previous_evaluator_bit_for_bit() {
    set_scheme_override(Some(Scheme::Reference));
    let options = EwaldOptions::default();
    for (cell, ions) in cases() {
        let want = bits(&old_energy(&cell, &ions, options).unwrap());
        assert_eq!(bits(&energy(&cell, &ions, options).unwrap()), want);
        let plan = EwaldPlan::with_scheme(&cell, &ions, options, Scheme::Reference).unwrap();
        assert_partitions(&plan, want);
        let f_old = old_forces(&cell, &ions, options).unwrap();
        let f_new = forces(&cell, &ions, options).unwrap();
        for (a, b) in f_old.iter().zip(&f_new) {
            assert_eq!(a.map(f64::to_bits), b.map(f64::to_bits));
        }
        let s_old = old_stress(&cell, &ions, options).unwrap();
        let s_new = stress(&cell, &ions, options).unwrap();
        assert_eq!(
            s_old.map(|r| r.map(f64::to_bits)),
            s_new.map(|r| r.map(f64::to_bits))
        );
    }
    set_scheme_override(None);
}

/// Balanced scheme: the serial evaluator equals every partition of its plan
/// bit for bit.
#[test]
fn balanced_plan_partitions_match_the_serial_sum_bit_for_bit() {
    let options = EwaldOptions::default();
    for rows in shapes(1.0) {
        for atoms in [1usize, 2, 8] {
            let cell = Cell::from_user_rows_bohr(rows).unwrap();
            let ions = ions(atoms);
            let want = bits(&energy_with(&cell, &ions, options, Scheme::Balanced).unwrap());
            let plan = EwaldPlan::with_scheme(&cell, &ions, options, Scheme::Balanced).unwrap();
            assert_eq!(plan.scheme(), Scheme::Balanced);
            assert_partitions(&plan, want);
        }
    }
}

/// Coincident ions fail in the serial evaluator and in the plan, both schemes.
#[test]
fn coincident_ions_fail_in_both_evaluators() {
    let cell = Cell::from_user_rows_bohr([[6., 0., 0.], [0., 6., 0.], [0., 0., 6.]]).unwrap();
    let ion = Ion {
        charge: 4.0,
        fractional_position: [0.1, 0.2, 0.3],
    };
    let ions = [ion, ion];
    let options = EwaldOptions::default();
    for scheme in [Scheme::Reference, Scheme::Balanced] {
        assert!(matches!(
            energy_with(&cell, &ions, options, scheme),
            Err(PhysicsError::CoincidentIons)
        ));
        let plan = EwaldPlan::with_scheme(&cell, &ions, options, scheme).unwrap();
        let mut out = Vec::new();
        assert!(matches!(
            plan.real_terms(0..plan.real_candidates(), &mut out),
            Err(PhysicsError::CoincidentIons)
        ));
    }
}

/// Neumaier-compensated sum.
#[derive(Default)]
struct Compensated {
    sum: f64,
    c: f64,
}
impl Compensated {
    fn add(&mut self, x: f64) {
        let t = self.sum + x;
        if self.sum.abs() >= x.abs() {
            self.c += (self.sum - t) + x;
        } else {
            self.c += (x - t) + self.sum;
        }
        self.sum = t;
    }
    fn value(&self) -> f64 {
        self.sum + self.c
    }
}

/// Independent accuracy reference of the energy: the eecd57f truncation
/// protocol at `tolerance` with splitting parameter `eta`, full enumeration
/// boxes, every term (including the structure factors) accumulated with
/// compensated summation, so that its rounding stays far below the cap even
/// for 10^7 terms.
fn compensated_energy(cell: &Cell, ions: &[Ion], eta: f64, tolerance: f64) -> f64 {
    let options = EwaldOptions {
        eta: Some(eta),
        tolerance,
    };
    let e = old_energy(cell, ions, options).unwrap();
    assert!(e.tail_bound.0 < 1e-14 * ions.len() as f64);
    let a = matrix(cell.lattice_columns_bohr());
    let b = matrix(cell.reciprocal_columns_inverse_bohr());
    let rb = bounds(&b, e.real_cutoff_bohr, 1.0).unwrap();
    let gb = bounds(&a, e.reciprocal_cutoff_inverse_bohr, 0.0).unwrap();
    let positions = positions_of(ions);
    let mut total = Compensated::default();
    for (i, x) in ions.iter().enumerate() {
        for (j, y) in ions.iter().enumerate() {
            enumerate(
                |m| {
                    if i == j && m == [0; 3] {
                        return Ok(());
                    }
                    let d = (a * (positions[i] - positions[j] + Vector3::from(m.map(f64::from))))
                        .norm();
                    if d <= e.real_cutoff_bohr {
                        total.add(0.5 * x.charge * y.charge * libm::erfc(eta * d) / d);
                    }
                    Ok(())
                },
                rb,
            )
            .unwrap();
        }
    }
    let omega = cell.volume_bohr3();
    enumerate(
        |m| {
            if m == [0; 3] {
                return Ok(());
            }
            let mv = Vector3::from(m.map(f64::from));
            let g2 = (b * mv).norm_squared();
            if g2 <= e.reciprocal_cutoff_inverse_bohr.powi(2) {
                let (mut re, mut im) = (Compensated::default(), Compensated::default());
                for (ion, s) in ions.iter().zip(&positions) {
                    let z = Complex64::from_polar(ion.charge, TAU * mv.dot(s));
                    re.add(z.re);
                    im.add(z.im);
                }
                let s2 = re.value().powi(2) + im.value().powi(2);
                total.add((2.0 * PI / omega) * (-g2 / (4.0 * eta * eta)).exp() / g2 * s2);
            }
            Ok(())
        },
        gb,
    )
    .unwrap();
    total.add(e.self_energy.0);
    total.add(e.background.0);
    total.value()
}

/// Largest deviations of the balanced scheme from tight references: energy
/// per atom (Ha) against the compensated reference (at tolerance 1e-16, eta =
/// sqrt(pi)/Omega^(1/3); the same reference at twice that eta must agree with
/// it to 1e-13 Ha/atom, which bounds the reference's own error); force
/// component (Ha/Bohr) and stress component (Ha/Bohr^3) against the eecd57f
/// evaluator at tolerance 1e-16.
fn balanced_errors(cell: &Cell, ions: &[Ion]) -> (f64, f64, f64) {
    let tight = EwaldOptions {
        eta: None,
        tolerance: 1e-16,
    };
    let options = EwaldOptions::default();
    let eta0 = PI.sqrt() / cell.volume_bohr3().cbrt();
    let reference = compensated_energy(cell, ions, eta0, 1e-16);
    let other = compensated_energy(cell, ions, 2.0 * eta0, 1e-16);
    let n = ions.len() as f64;
    assert!(
        (reference - other).abs() / n < 1e-13,
        "reference eta dependence {}",
        (reference - other).abs() / n
    );
    set_scheme_override(Some(Scheme::Balanced));
    let fast = energy(cell, ions, options).unwrap();
    let f_fast = forces(cell, ions, options).unwrap();
    let s_fast = stress(cell, ions, options).unwrap();
    set_scheme_override(None);
    assert!(fast.tail_bound.0 <= 1e-13 * ions.len() as f64);
    let f_ref = old_forces(cell, ions, tight).unwrap();
    let s_ref = old_stress(cell, ions, tight).unwrap();
    let e = (fast.energy.0 - reference).abs() / n;
    let f = f_fast
        .iter()
        .flatten()
        .zip(f_ref.iter().flatten())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    let s = s_fast
        .iter()
        .flatten()
        .zip(s_ref.iter().flatten())
        .map(|(a, b)| (a - b).abs())
        .fold(0.0, f64::max);
    (e, f, s)
}

/// Pre-registered caps (INIT-EWALD-FAST): energy 1e-12 Ha/atom, force
/// 1e-10 Ha/Bohr, stress 1e-12 Ha/Bohr^3, against the tight reference.
fn check_balanced(scale: f64, counts: &[usize]) {
    for (shape, rows) in shapes(scale).into_iter().enumerate() {
        for &atoms in counts {
            let cell = Cell::from_user_rows_bohr(rows).unwrap();
            let (e, f, s) = balanced_errors(&cell, &ions(atoms));
            println!(
                "shape {shape} scale {scale} atoms {atoms}: |dE|/N {e:.3e} Ha, |dF| {f:.3e} Ha/Bohr, |dS| {s:.3e} Ha/Bohr^3"
            );
            assert!(e < 1e-12, "energy {e}");
            assert!(f < 1e-10, "force {f}");
            assert!(s < 1e-12, "stress {s}");
        }
    }
}

#[test]
fn balanced_scheme_matches_a_tight_reference_small_cells() {
    check_balanced(1.0, &[1, 2, 8]);
}

/// 64 ions in cells twice as large; about 10 s in debug, so release only:
/// `cargo test --release -p dft-physics --lib ewald::tests::balanced_scheme_matches_a_tight_reference_64 -- --ignored --nocapture`
#[test]
#[ignore]
fn balanced_scheme_matches_a_tight_reference_64() {
    check_balanced(2.0, &[64]);
}
