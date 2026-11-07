use super::*;
use dft_geometry::Cell;
use dft_physics::data::Array3;
use dft_physics::hamiltonian::NonlocalOperator;

#[test]
fn operator_options_bands_and_executor_invalidate_cached_state() {
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    let mut options = ScfOptions {
        gamma_real: GammaReal::On,
        cutoff_hartree: 2.5,
        ..Default::default()
    };
    let basis = options.basis(&cell, [0.; 3], 1.).unwrap();
    let grid =
        dft_basis::FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
    let make = || {
        Hamiltonian::new(
            &cell,
            &basis,
            Array3::zeros(grid.dims()),
            NonlocalOperator::new(&basis, cell.volume_bohr3(), &[]).unwrap(),
        )
        .unwrap()
    };
    let mut h = make();
    let mut state = None;
    let mut diagnostics = GammaRealDiagnostics::default();
    let mut pool = WorkerPool::new(1).unwrap();
    let mut previous_storage = None;
    for step in 0..7 {
        if step == 1 {
            h.set_local_potential(Array3::from_elem(grid.dims(), 0.25))
                .unwrap();
        }
        if step == 2 {
            options.seed += 1;
        }
        if step == 4 {
            pool = WorkerPool::new(2).unwrap();
        }
        if step == 5 {
            options.davidson_capacity_per_band = 2;
        }
        if step == 6 {
            h = make();
        }
        let bands = if step < 3 { 3 } else { 4 };
        let policy = BandTolerances {
            strict_bands: bands,
            buffer_tolerance: 1e-9,
        };
        let run = |state: &mut Option<RealSolveState>, diagnostics: &mut GammaRealDiagnostics| {
            solve_real_gamma(
                &h,
                &basis,
                None,
                h.kinetic_diagonal(),
                bands,
                1e-9,
                policy,
                true,
                &options,
                pool.clone(),
                None,
                diagnostics,
                state,
            )
            .unwrap()
            .unwrap()
        };
        let result = run(&mut state, &mut diagnostics);
        let current_storage = state
            .as_ref()
            .unwrap()
            .op
            .h
            .borrow()
            .half_sphere()
            .plus()
            .as_ptr();
        if let Some(old) = previous_storage {
            if step == 1 {
                assert_eq!(old, current_storage);
            } else {
                assert_ne!(old, current_storage);
            }
        }
        previous_storage = Some(current_storage);
        let fresh = run(&mut None, &mut GammaRealDiagnostics::default());
        assert_eq!(format!("{result:?}"), format!("{fresh:?}"));
    }
}
