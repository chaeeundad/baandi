//! P0.2 work timers in the real Gamma `h_apply`: timed applications are
//! bit-identical to untimed ones (plain, cached projections, overlap) and the
//! counts of calls, columns and pair transforms are exact at every pool budget.
//! Own binary: the counters are process-global. Synthetic NC data.
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_numerics::{RealMatrix, pool::WorkerPool, work_timers};
use dft_physics::{
    RealHamiltonian,
    hamiltonian::{Hamiltonian, NonlocalOperator, Species},
};
use dft_units::Hartree;
use ndarray::Array3;

fn pseudo() -> dft_upf::Pseudopotential {
    let n = 81;
    let a = |f: &dyn Fn(f64) -> f64| {
        (0..n)
            .map(|i| format!("{:.17e}", f(i as f64 * 0.15)))
            .collect::<Vec<_>>()
            .join(" ")
    };
    let beta = format!(
        "<PP_BETA.1 index=\"1\" angular_momentum=\"1\" cutoff_radius=\"12\" cutoff_radius_index=\"{n}\">{}</PP_BETA.1>",
        a(&|r| r * r * (-r * r).exp())
    );
    dft_upf::parse_upf(&format!(r#"<UPF version="2.0.1"><PP_INFO>Synthetic P0.2 fixture</PP_INFO>
<PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_wfc="F" has_gipaw="F" core_correction="F" functional="SLA PZ" z_valence="4" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max="1" l_max_rho="0" l_local="-1" mesh_size="{n}" number_of_wfc="0" number_of_proj="1"/>
<PP_MESH dx="0.15" mesh="{n}" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH><PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL>{beta}<PP_DIJ>1.3</PP_DIJ></PP_NONLOCAL><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,a(&|r|r),a(&|_|0.15),a(&|_|0.),a(&|_|0.))).unwrap()
}
fn bits(x: &RealMatrix) -> Vec<u64> {
    x.as_slice().iter().map(|v| v.to_bits()).collect()
}

#[test]
fn h_apply_timers_change_no_bit_and_count_exactly() {
    let cell = Cell::from_user_rows_bohr([[6., 0.4, 0.2], [1., 7., 0.3], [0.5, 0.7, 8.]]).unwrap();
    let basis = PlaneWaveBasis::new(&cell, Hartree(1.2)).unwrap();
    let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
    let p = pseudo();
    let positions = [
        cell.fractional_to_cartesian_bohr([0.1, 0.2, 0.3]),
        cell.fractional_to_cartesian_bohr([0.6, 0.5, 0.8]),
    ];
    let nl = NonlocalOperator::new(
        &basis,
        cell.volume_bohr3(),
        &[Species {
            potential: &p,
            positions_bohr: &positions,
        }],
    )
    .unwrap();
    let v = Array3::from_shape_fn(grid.dims(), |(i, j, k)| {
        0.3 * ((i * 13 + j * 7 + k * 3) as f64).sin()
    });
    let h = Hamiltonian::new(&cell, &basis, v, nl).unwrap();
    let n = basis.vectors().len();
    let m = 5;
    let x = RealMatrix::new(
        n,
        m,
        (0..n * m)
            .map(|i| ((i * 17 + 3) as f64).sin() / (n as f64).sqrt())
            .collect(),
    )
    .unwrap();
    let fine = |s: [work_timers::WorkTimes; work_timers::RUN_STAGES]| s[0];

    // Reference without timers; nothing is recorded.
    work_timers::set_enabled(false);
    let before = work_timers::snapshot();
    let mut real =
        RealHamiltonian::from_complex(&h, &basis, m, WorkerPool::new(1).unwrap()).unwrap();
    let mut plain = RealMatrix::zeros(n, m).unwrap();
    real.apply_block(&x, &mut plain).unwrap();
    let projections = real.last_projections(m).unwrap().to_vec();
    let mut cached = RealMatrix::zeros(n, m).unwrap();
    real.apply_block_with_projections(&x, &projections, &mut cached)
        .unwrap();
    assert_eq!(work_timers::snapshot(), before, "off must record nothing");

    work_timers::set_enabled(true);
    for threads in [1, 2, 4] {
        for overlap in [false, true] {
            let mut real =
                RealHamiltonian::from_complex(&h, &basis, m, WorkerPool::new(threads).unwrap())
                    .unwrap();
            real.set_apply_overlap(overlap);
            let start = fine(work_timers::snapshot());
            let mut y = RealMatrix::zeros(n, m).unwrap();
            real.apply_block(&x, &mut y).unwrap();
            assert_eq!(bits(&y), bits(&plain), "T={threads} overlap={overlap}");
            let mut yc = RealMatrix::zeros(n, m).unwrap();
            real.apply_block_with_projections(&x, &projections, &mut yc)
                .unwrap();
            assert_eq!(bits(&yc), bits(&cached), "T={threads} overlap={overlap}");
            let w = fine(work_timers::snapshot());
            assert_eq!(w.apply_calls - start.apply_calls, 2);
            assert_eq!(w.apply_columns - start.apply_columns, 2 * m as u64);
            // Two applications of five columns: three pairs each.
            assert_eq!(w.pair_transforms - start.pair_transforms, 6);
            assert_eq!(w.pair_columns - start.pair_columns, 2 * m as u64);
            assert!(w.local > start.local);
            assert!(w.nonlocal_projection > start.nonlocal_projection);
            assert!(w.nonlocal_back > start.nonlocal_back);
            assert!(w.apply - start.apply >= w.nonlocal_projection - start.nonlocal_projection);
        }
    }
    work_timers::set_enabled(false);
}
