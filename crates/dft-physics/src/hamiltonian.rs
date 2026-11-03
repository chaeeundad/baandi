//! Ionic and Kleinman–Bylander terms, and a reusable Bloch Hamiltonian.
use dft_basis::{FftGridSpec, GVector, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_numerics::{
    Fft3d, FftPlan3d,
    fused::{FusedOperator, FusedPlan},
    linalg::ComplexMatrix,
};
use dft_solvers::{LinearOperator, SolverError};
use dft_upf::{Pseudopotential, radial::real_spherical_harmonic};
use ndarray::Array3;
use num_complex::Complex64;
use std::cell::RefCell;
use thiserror::Error;

#[derive(Debug, Error)]
pub enum HamiltonianError {
    #[error("incompatible basis, grid, or vector dimensions")]
    Shape,
    #[error("nonfinite input or result")]
    NonFinite,
    #[error("volume must be finite and positive")]
    Volume,
    #[error("D couples different angular momenta")]
    AngularCoupling,
    #[error("Hamiltonian scratch is already borrowed")]
    Busy,
    #[error(transparent)]
    Basis(#[from] dft_basis::BasisError),
    #[error(transparent)]
    Fft(#[from] dft_numerics::FftError),
    #[error(transparent)]
    Radial(#[from] dft_upf::radial::RadialError),
}
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
pub struct Species<'a> {
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
pub fn ionic_local(
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
pub struct NonlocalOperator {
    basis: Vec<GVector>,
    atoms: Vec<AtomProjectors>,
    max_channels: usize,
}
impl NonlocalOperator {
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
    pub fn atoms(&self) -> &[AtomProjectors] {
        &self.atoms
    }
    fn add(&self, x: &[Complex64], y: &mut [Complex64], projections: &mut [Complex64]) {
        for atom in &self.atoms {
            let nc = atom.columns.len();
            for (p, col) in projections.iter_mut().zip(&atom.columns) {
                *p = col.iter().zip(x).map(|(b, c)| b.conj() * c).sum();
            }
            for (i, col) in atom.columns.iter().enumerate() {
                let weight: Complex64 = (0..nc)
                    .map(|j| atom.coupling[i * nc + j] * projections[j])
                    .sum();
                for (out, b) in y.iter_mut().zip(col) {
                    *out += b * weight;
                }
            }
        }
    }
    pub fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), HamiltonianError> {
        if x.len() != self.basis.len() || y.len() != x.len() {
            return Err(HamiltonianError::Shape);
        }
        if !x.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        y.fill(Complex64::default());
        self.add(x, y, &mut vec![Complex64::default(); self.max_channels]);
        if !y.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        Ok(())
    }
}

/// Filter represented Fourier coefficients to the density sphere, then invert.
/// Already aliased samples cannot be recovered by filtering.
pub fn truncate_to_density_sphere(
    values: &Array3<f64>,
    density: &PlaneWaveBasis,
) -> Result<Array3<f64>, HamiltonianError> {
    let dims = [values.shape()[0], values.shape()[1], values.shape()[2]];
    let grid = FftGridSpec::new(density, density, Some(dims))?;
    if !values.iter().all(|x| x.is_finite()) {
        return Err(HamiltonianError::NonFinite);
    }
    let mut fft = FftPlan3d::new(dims)?;
    let mut real = Array3::from_shape_fn(dims, |idx| Complex64::new(values[idx], 0.0));
    let mut reciprocal = Array3::zeros(dims);
    fft.to_reciprocal(&real, &mut reciprocal)?;
    for (z, index) in reciprocal.iter_mut().zip(grid.fft_to_g()) {
        if index.is_none() {
            *z = Complex64::default();
        }
    }
    fft.to_real(&reciprocal, &mut real)?;
    Ok(real.mapv(|z| z.re))
}
struct Worker {
    scratch: Scratch,
    output: Vec<Complex64>,
}
struct Scratch {
    fft: FusedPlan,
    projections: Vec<Complex64>,
}
/// Immutable physical data with interior scratch reuse; intentionally not Sync.
/// Local grid values are in Ha and must use the internally derived grid dimensions.
pub struct Hamiltonian {
    data: HamiltonianData,
    scratch: RefCell<Vec<Worker>>,
}
struct HamiltonianData {
    grid: FftGridSpec,
    kinetic: Vec<f64>,
    potential: Array3<f64>,
    nonlocal: NonlocalOperator,
}
impl Hamiltonian {
    /// Bloch Hamiltonian for the supplied basis: |k+G|^2/2 + V_local + V_NL(k).
    /// FFTs act on the periodic part; local samples are on the common density grid.
    pub fn new(
        cell: &Cell,
        basis: &PlaneWaveBasis,
        potential: Array3<f64>,
        nonlocal: NonlocalOperator,
    ) -> Result<Self, HamiltonianError> {
        if basis.momenta() != nonlocal.basis {
            return Err(HamiltonianError::Shape);
        }
        let dims = [
            potential.shape()[0],
            potential.shape()[1],
            potential.shape()[2],
        ];
        let grid = FftGridSpec::new(basis, &basis.density_basis(cell)?, Some(dims))?;
        if !potential.iter().all(|x| x.is_finite()) {
            return Err(HamiltonianError::NonFinite);
        }
        let scratch = RefCell::new(vec![Worker {
            scratch: Scratch {
                fft: FusedPlan::new(dims, grid.g_to_fft())?,
                projections: vec![Complex64::default(); nonlocal.max_channels],
            },
            output: vec![Complex64::default(); basis.vectors().len()],
        }]);
        Ok(Self {
            data: HamiltonianData {
                grid,
                kinetic: basis
                    .momenta()
                    .iter()
                    .map(|g| 0.5 * g.norm_squared)
                    .collect(),
                potential: if potential.is_standard_layout() {
                    potential
                } else {
                    potential.as_standard_layout().into_owned()
                },
                nonlocal,
            },
            scratch,
        })
    }
    /// Set the intra-k budget once during construction. Zero means one.
    /// Numerical buffers are reused; only scoped OS-thread bookkeeping is transient.
    pub fn with_band_threads(mut self, threads: usize) -> Result<Self, HamiltonianError> {
        let dims = self.data.grid.dims();
        let pool = self.scratch.get_mut();
        for _ in pool.len()..threads.max(1) {
            pool.push(Worker {
                scratch: Scratch {
                    fft: FusedPlan::new(dims, self.data.grid.g_to_fft())?,
                    projections: vec![Complex64::default(); self.data.nonlocal.max_channels],
                },
                output: vec![Complex64::default(); self.data.kinetic.len()],
            });
        }
        pool.truncate(threads.max(1));
        Ok(self)
    }
    /// Replace only the effective local samples, retaining FFT and projector plans.
    /// Invalid input leaves the previous potential intact.
    pub fn set_local_potential(&mut self, potential: Array3<f64>) -> Result<(), HamiltonianError> {
        if potential.shape() != self.data.potential.shape() {
            return Err(HamiltonianError::Shape);
        }
        if !potential.iter().all(|v| v.is_finite()) {
            return Err(HamiltonianError::NonFinite);
        }
        self.data.potential = if potential.is_standard_layout() {
            potential
        } else {
            potential.as_standard_layout().into_owned()
        };
        Ok(())
    }
    /// Read-only access for term-by-term orbital energy audits.
    pub fn nonlocal(&self) -> &NonlocalOperator {
        &self.data.nonlocal
    }
    pub fn kinetic_diagonal(&self) -> &[f64] {
        &self.data.kinetic
    }
}
impl HamiltonianData {
    fn apply_with(
        &self,
        x: &[Complex64],
        y: &mut [Complex64],
        scratch: &mut Scratch,
    ) -> Result<(), HamiltonianError> {
        if x.len() != self.kinetic.len() || y.len() != x.len() {
            return Err(HamiltonianError::Shape);
        }
        if !x.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        scratch.fft.apply_local(
            x,
            self.potential.as_slice().ok_or(HamiltonianError::Shape)?,
            y,
        )?;
        for (i, value) in y.iter_mut().enumerate() {
            *value += self.kinetic[i] * x[i];
        }
        self.nonlocal.add(x, y, &mut scratch.projections);
        if !y.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        Ok(())
    }
}
fn application(error: HamiltonianError) -> SolverError {
    SolverError::Application(error.to_string())
}
impl LinearOperator for Hamiltonian {
    fn dimension(&self) -> usize {
        self.data.kinetic.len()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        let mut scratch = self
            .scratch
            .try_borrow_mut()
            .map_err(|_| application(HamiltonianError::Busy))?;
        self.data
            .apply_with(x, y, &mut scratch[0].scratch)
            .map_err(application)
    }
    fn apply_block(&self, x: &ComplexMatrix, y: &mut ComplexMatrix) -> Result<(), SolverError> {
        if x.rows() != self.dimension() || y.rows() != x.rows() || y.cols() != x.cols() {
            return Err(application(HamiltonianError::Shape));
        }
        let mut scratch = self
            .scratch
            .try_borrow_mut()
            .map_err(|_| application(HamiltonianError::Busy))?;
        if scratch.len() == 1 || x.cols() <= 1 {
            for j in 0..x.cols() {
                self.data
                    .apply_with(x.column(j), y.column_mut(j), &mut scratch[0].scratch)
                    .map_err(|e| SolverError::Application(format!("band index {j}: {e}")))?;
            }
            return Ok(());
        }
        let workers = scratch.len().min(x.cols());
        let output = std::sync::Mutex::new(y);
        apply_workers(&self.data, &mut scratch[..workers], x, &output, 0, x.cols())
    }
}

// Each recursion uses its caller as one worker, so at most pool.len() threads
// execute numerical work. Joins precede error propagation in ascending band order.
fn apply_workers(
    data: &HamiltonianData,
    workers: &mut [Worker],
    x: &ComplexMatrix,
    y: &std::sync::Mutex<&mut ComplexMatrix>,
    first: usize,
    end: usize,
) -> Result<(), SolverError> {
    if workers.len() == 1 {
        let Worker { scratch, output } = &mut workers[0];
        for j in first..end {
            data.apply_with(x.column(j), output, scratch)
                .map_err(|e| SolverError::Application(format!("band index {j}: {e}")))?;
            // The lock covers only a copy into a disjoint column, never arithmetic.
            y.lock()
                .map_err(|_| SolverError::Application("block output poisoned".into()))?
                .column_mut(j)
                .copy_from_slice(output);
        }
        return Ok(());
    }
    let left_count = workers.len() / 2;
    let middle = first + (end - first) * left_count / workers.len();
    let (left, right) = workers.split_at_mut(left_count);
    std::thread::scope(|scope| {
        let handle = scope.spawn(|| apply_workers(data, left, x, y, first, middle));
        let right = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
            apply_workers(data, right, x, y, middle, end)
        }));
        let left = handle.join();
        left.map_err(|_| {
            SolverError::Application(format!("band index {first}: worker panicked"))
        })??;
        right.map_err(|_| {
            SolverError::Application(format!("band index {middle}: worker panicked"))
        })?
    })
}

#[cfg(test)]
mod pool_tests {
    use super::*;
    #[test]
    fn scratch_addresses_survive_blocks_and_errors() {
        let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
        let basis = PlaneWaveBasis::new(&cell, dft_units::Hartree(2.5)).unwrap();
        let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
        let h = Hamiltonian::new(
            &cell,
            &basis,
            Array3::zeros(grid.dims()),
            NonlocalOperator::new(&basis, cell.volume_bohr3(), &[]).unwrap(),
        )
        .unwrap()
        .with_band_threads(3)
        .unwrap();
        let addresses = || {
            h.scratch
                .borrow()
                .iter()
                .map(|w| {
                    (
                        &w.scratch.fft as *const FusedPlan as usize,
                        w.scratch.projections.as_ptr() as usize,
                        w.output.as_ptr() as usize,
                    )
                })
                .collect::<Vec<_>>()
        };
        let before = addresses();
        let mut x = ComplexMatrix::zeros(basis.vectors().len(), 5).unwrap();
        let mut y = x.clone();
        for _ in 0..3 {
            h.apply_block(&x, &mut y).unwrap();
            assert_eq!(before, addresses());
        }
        x.column_mut(1)[0].re = f64::NAN;
        assert!(h.apply_block(&x, &mut y).is_err());
        assert_eq!(before, addresses());
    }
}

#[cfg(test)]
mod fused_tests {
    use super::*;
    use dft_numerics::SpherePattern;
    use std::{hint::black_box, time::Instant};
    struct Legacy {
        fft: FftPlan3d,
        pattern: SpherePattern,
        reciprocal: Array3<Complex64>,
        real: Array3<Complex64>,
        projections: Vec<Complex64>,
    }
    impl Legacy {
        fn new(h: &Hamiltonian) -> Self {
            let dims = h.data.grid.dims();
            Self {
                fft: FftPlan3d::new(dims).unwrap(),
                pattern: SpherePattern::from_indices(dims, h.data.grid.g_to_fft()).unwrap(),
                reciprocal: Array3::zeros(dims),
                real: Array3::zeros(dims),
                projections: vec![Complex64::default(); h.data.nonlocal.max_channels],
            }
        }
        fn apply(&mut self, h: &Hamiltonian, x: &[Complex64], y: &mut [Complex64]) {
            assert!(x.iter().all(finite));
            self.reciprocal.fill(Complex64::default());
            for (&i, &z) in h.data.grid.g_to_fft().iter().zip(x) {
                self.reciprocal[i] = z;
            }
            self.fft
                .to_real_pruned(&self.reciprocal, &self.pattern, &mut self.real)
                .unwrap();
            for (z, &v) in self.real.iter_mut().zip(&h.data.potential) {
                *z *= v;
            }
            self.fft
                .to_reciprocal_pruned(&self.real, &self.pattern, &mut self.reciprocal)
                .unwrap();
            for (i, &idx) in h.data.grid.g_to_fft().iter().enumerate() {
                y[i] = self.reciprocal[idx] + h.data.kinetic[i] * x[i];
            }
            h.data.nonlocal.add(x, y, &mut self.projections);
            assert!(y.iter().all(finite));
        }
    }
    fn vector(n: usize, seed: usize) -> Vec<Complex64> {
        (0..n)
            .map(|i| {
                Complex64::new(
                    ((i * 17 + seed) as f64).sin(),
                    ((i * 31 + seed) as f64).cos(),
                )
            })
            .collect()
    }
    fn bits(v: &[Complex64]) -> Vec<(u64, u64)> {
        v.iter().map(|z| (z.re.to_bits(), z.im.to_bits())).collect()
    }
    fn synthetic(shape: [usize; 3], k: [f64; 3], channels: usize) -> Hamiltonian {
        let cell =
            Cell::from_user_rows_bohr([[5., 0., 0.], [0.7, 6., 0.], [0.2, 0.3, 7.]]).unwrap();
        let basis = PlaneWaveBasis::at_k(&cell, dft_units::Hartree(2.), k).unwrap();
        let nonlocal = NonlocalOperator {
            basis: basis.momenta().to_vec(),
            max_channels: channels,
            atoms: vec![AtomProjectors {
                columns: (0..channels)
                    .map(|i| vector(basis.vectors().len(), i))
                    .collect(),
                coupling: (0..channels * channels)
                    .map(|i| {
                        if i / channels == i % channels {
                            0.03
                        } else {
                            0.01
                        }
                    })
                    .collect(),
            }],
        };
        let potential = Array3::from_shape_vec(
            shape,
            vector(shape.iter().product(), 4)
                .iter()
                .map(|z| z.re)
                .collect(),
        )
        .unwrap();
        Hamiltonian::new(&cell, &basis, potential, nonlocal).unwrap()
    }
    #[test]
    fn ordered_complete_hamiltonian_shapes_and_projectors() {
        for shape in [
            [24, 24, 24],
            [25, 25, 25],
            [27, 27, 27],
            [36, 36, 36],
            [15, 20, 24],
            [15, 17, 19],
        ] {
            for (k, channels) in [([0.; 3], 0), ([0.25, -0.25, 0.125], 3), ([0.125; 3], 5)] {
                let h = synthetic(shape, k, channels);
                let mut old = Legacy::new(&h);
                let x = vector(h.dimension(), 12);
                let mut expected = x.clone();
                let mut actual = x.clone();
                old.apply(&h, &x, &mut expected);
                h.apply(&x, &mut actual).unwrap();
                assert_eq!(bits(&actual), bits(&expected));
            }
        }
    }
    #[test]
    #[ignore = "P10C_SI_UPF=/external/Si.upf cargo test --release -p dft-physics complete_hamiltonian_timing -- --ignored --nocapture"]
    fn complete_hamiltonian_timing() {
        let pp =
            dft_upf::from_file(std::env::var("P10C_SI_UPF").expect("external Si UPF required"))
                .unwrap();
        let a = 10.261212856717933;
        let primitive = Cell::from_user_rows_bohr([
            [0., a / 2., a / 2.],
            [a / 2., 0., a / 2.],
            [a / 2., a / 2., 0.],
        ])
        .unwrap();
        let conventional =
            Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap();
        let mut cases = Vec::new();
        for (label, cell, sites, k) in [
            (
                "Si-Gamma",
                primitive.clone(),
                vec![[0.; 3], [0.25; 3]],
                [0.; 3],
            ),
            ("Si-shifted", primitive, vec![[0.; 3], [0.25; 3]], [0.25; 3]),
            (
                "Si8",
                conventional,
                vec![
                    [0.; 3],
                    [0.25; 3],
                    [0., 0.5, 0.5],
                    [0.5, 0., 0.5],
                    [0.5, 0.5, 0.],
                    [0.25, 0.75, 0.75],
                    [0.75, 0.25, 0.75],
                    [0.75, 0.75, 0.25],
                ],
                [0.; 3],
            ),
        ] {
            let basis = PlaneWaveBasis::at_k(&cell, dft_units::Hartree(15.), k).unwrap();
            let grid =
                FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
            let positions: Vec<_> = sites
                .iter()
                .map(|&r| cell.fractional_to_cartesian_bohr(r))
                .collect();
            let nonlocal = NonlocalOperator::new(
                &basis,
                cell.volume_bohr3(),
                &[Species {
                    potential: &pp,
                    positions_bohr: &positions,
                }],
            )
            .unwrap();
            let potential = Array3::from_shape_vec(
                grid.dims(),
                vector(grid.dims().iter().product(), 4)
                    .iter()
                    .map(|z| z.re)
                    .collect(),
            )
            .unwrap();
            cases.push((
                label,
                Hamiltonian::new(&cell, &basis, potential, nonlocal).unwrap(),
            ));
        }
        cases.push((
            "synthetic-skew-3-projectors",
            synthetic([15, 20, 24], [0.25, -0.25, 0.125], 3),
        ));
        for (label, h) in cases {
            let mut old = Legacy::new(&h);
            let x = vector(h.dimension(), 17);
            let mut out = x.clone();
            let mut expected = x.clone();
            old.apply(&h, &x, &mut expected);
            h.apply(&x, &mut out).unwrap();
            assert_eq!(bits(&expected), bits(&out));
            let mut before = Vec::new();
            let mut after = Vec::new();
            for _ in 0..5 {
                let t = Instant::now();
                for _ in 0..128 {
                    old.apply(&h, black_box(&x), &mut out);
                    black_box(&out);
                }
                before.push(t.elapsed().as_secs_f64() * 1e6 / 128.);
                let t = Instant::now();
                for _ in 0..128 {
                    h.apply(black_box(&x), &mut out).unwrap();
                    black_box(&out);
                }
                after.push(t.elapsed().as_secs_f64() * 1e6 / 128.);
            }
            println!(
                "{label} shape={:?} npw={} before_us={before:?} after_us={after:?}",
                h.data.grid.dims(),
                h.dimension()
            );
        }
    }
}
