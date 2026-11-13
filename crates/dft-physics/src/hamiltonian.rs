//! Ionic and Kleinman–Bylander terms, and a reusable Bloch Hamiltonian.
use dft_basis::{FftGridSpec, GVector, PlaneWaveBasis};
use dft_geometry::Cell;
use dft_numerics::pool::WorkerPool;
use dft_numerics::pool::analysis::{Piece, Span};
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
use std::sync::Arc;
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
        // Keep LocalTransform's endpoint and quadrature, but hoist only the
        // G-independent parenthesized expression. In particular, rab must
        // still multiply AFTER j0: reassociation would change the bits.
        let pp = s.potential;
        let n = pp
            .mesh()
            .r()
            .partition_point(|r| *r <= local.effective_r_max());
        let mesh = dft_upf::radial::RadialMesh::new(
            pp.mesh().r()[..n].to_vec(),
            pp.mesh().rab()[..n].to_vec(),
        )?;
        let short: Vec<_> = mesh
            .r()
            .iter()
            .zip(pp.local_hartree())
            .map(|(&r, &v)| r * r * v + pp.z_valence() * r * libm::erf(r))
            .collect();
        let mut integrand = vec![0.0; n];
        // At most one entry per G; reserve up front so misses never allocate.
        let mut cache = std::collections::HashMap::with_capacity(density.vectors().len());
        for (out, g) in result.iter_mut().zip(density.vectors()) {
            let v = if g.miller == [0; 3] {
                local.alpha_z()?
            } else {
                let key = g.norm_squared.to_bits();
                if let Some(&v) = cache.get(&key) {
                    v
                } else {
                    let q = g.norm_squared.sqrt();
                    if !q.is_finite() || q <= 0.0 || !(q * local.effective_r_max()).is_finite() {
                        return Err(dft_upf::radial::RadialError::Argument(
                            "G must be finite and positive",
                        )
                        .into());
                    }
                    for ((value, &base), &r) in integrand.iter_mut().zip(&short).zip(mesh.r()) {
                        *value = base * dft_upf::radial::spherical_bessel(0, q * r)?;
                    }
                    let v = 4.0 * std::f64::consts::PI / volume_bohr3
                        * (mesh.integrate(&integrand)?
                            - pp.z_valence() * (-q * q / 4.0).exp() / (q * q));
                    if !v.is_finite() {
                        return Err(dft_upf::radial::RadialError::Argument(
                            "unrepresentable result",
                        )
                        .into());
                    }
                    cache.insert(key, v);
                    v
                }
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
/// Band vectors sharing one pass over a projector column in `expectation_values`.
const EXPECTATION_BLOCK: usize = 4;
/// P48C (`BAANDI_NONLOCAL_BLOCK=on`): fixed column panel of the block
/// nonlocal. A band block of width m is applied in ceil(m / 32) panels of at
/// most 32 bands (the last one narrower); each panel is two GEMMs whose shapes
/// depend only on the basis size, the projector count and the panel width.
/// A data-layout constant like the 256-row GEMM panel: it never depends on a
/// material, a thread budget or a run, so the bits of a band depend only on
/// the band block, not on how panels are distributed over workers.
const NONLOCAL_BAND_PANEL: usize = 32;
/// P48C shape gate of `BAANDI_NONLOCAL_BLOCK=auto` (the default): the block
/// path is used when the packed projector matrix of one k point,
/// n_pw x n_projectors x 16 bytes, is at least a quarter of a fixed 1 MiB L2
/// (256 KiB). Below it the fixed GEMM call and packing cost outweighs the
/// BLAS-1 per-band loops (dev08 2026-11-13: tiny cells lost 1 to 1.5 % at
/// T=8). A data-size rule with a fixed hardware-shape constant, the same style
/// as P48D's `ORTHO_PANEL_DOUBLES` (1 MiB is the smallest per-core L2 of the
/// measurement hosts); it never depends on a material, a thread budget or the
/// host, so the path choice (and every rounded bit) is the same everywhere.
pub const NONLOCAL_BLOCK_GATE_BYTES: usize = (1 << 20) / 4;

/// P48C auto-gate decision from shapes only: `true` when the packed projector
/// matrix of `basis_len` plane waves by `projectors` channels reaches
/// `NONLOCAL_BLOCK_GATE_BYTES`.
pub fn nonlocal_block_gate(basis_len: usize, projectors: usize) -> bool {
    basis_len
        .saturating_mul(projectors)
        .saturating_mul(size_of::<Complex64>())
        >= NONLOCAL_BLOCK_GATE_BYTES
}

/// P48C packed projector matrix B (n x p, column-major, channels in atom
/// order, each atom's channels as in `AtomProjectors::columns`). A copy of
/// the existing columns made once at construction; no radial recomputation.
struct BlockProjectors {
    packed: Vec<Complex64>,
    count: usize,
}

/// O(N_G * N_projectors) storage, never an N_G by N_G matrix.
pub struct NonlocalOperator {
    basis: Vec<GVector>,
    atoms: Vec<AtomProjectors>,
    max_channels: usize,
    block: Option<BlockProjectors>,
}
impl NonlocalOperator {
    /// Shared tail of the constructors: the packed block matrix is built only
    /// when the `BAANDI_NONLOCAL_BLOCK` mode selects the block path, so the
    /// per-band path carries no copy.
    fn assemble(basis: Vec<GVector>, atoms: Vec<AtomProjectors>, max_channels: usize) -> Self {
        let mut op = Self {
            basis,
            atoms,
            max_channels,
            block: None,
        };
        op.set_block_mode(dft_numerics::switches::choice(
            dft_numerics::switches::Id::NONLOCAL_BLOCK,
        ));
        op
    }
    /// Total projector channels (columns of the packed matrix).
    pub fn projector_count(&self) -> usize {
        self.atoms.iter().map(|a| a.columns.len()).sum()
    }
    /// Apply a `BAANDI_NONLOCAL_BLOCK` mode: `on` always uses the block path,
    /// `off` never, `auto` (any other value) uses it when
    /// `nonlocal_block_gate(n_pw, n_projectors)` holds. Depends on shapes only.
    pub fn set_block_mode(&mut self, mode: &str) {
        let on = match mode {
            "on" => true,
            "off" => false,
            _ => nonlocal_block_gate(self.basis.len(), self.projector_count()),
        };
        self.set_block(on);
    }
    /// Enable or disable the P48C block path (builds or drops the packed copy
    /// of the projector columns). Test and experiment hook; the engine uses the
    /// `BAANDI_NONLOCAL_BLOCK` mode read at construction.
    pub fn set_block(&mut self, on: bool) {
        if !on {
            self.block = None;
            return;
        }
        if self.block.is_some() {
            return;
        }
        let n = self.basis.len();
        let count = self.projector_count();
        let mut packed = Vec::with_capacity(n * count);
        for atom in &self.atoms {
            for col in &atom.columns {
                packed.extend_from_slice(&col[..n]);
            }
        }
        self.block = Some(BlockProjectors { packed, count });
    }
    /// Whether band blocks of width >= 2 use the P48C two-GEMM path.
    pub fn block_enabled(&self) -> bool {
        self.block.is_some()
    }
    /// Bytes of the packed projector copy (zero when the block path is off).
    pub fn block_bytes(&self) -> usize {
        self.block
            .as_ref()
            .map_or(0, |b| b.packed.len() * size_of::<Complex64>())
    }
    /// Per-worker scratch of one band panel: projections and weights
    /// (p x 32 each) and the image (n x 32); zero when the block path is off.
    fn block_scratch_len(&self) -> usize {
        self.block.as_ref().map_or(0, |b| {
            (2 * b.count + self.basis.len()) * NONLOCAL_BAND_PANEL
        })
    }
    /// Y += B D B^H X for one panel of `w` <= 32 bands: `x` and `y` are the
    /// contiguous n x w column-major panels. Pk = B^H X and the image B W are
    /// single sequential GEMMs of fixed shape; W = D Pk uses each atom's
    /// coupling block with the per-band order of `add` (sum over j ascending).
    /// Rounding class against `add` (different summation order).
    fn add_block_panel(
        &self,
        x: &[Complex64],
        y: &mut [Complex64],
        w: usize,
        scratch: &mut [Complex64],
    ) -> Result<(), HamiltonianError> {
        let Some(block) = self.block.as_ref() else {
            return Err(HamiltonianError::Shape);
        };
        let n = self.basis.len();
        let p = block.count;
        if x.len() != n * w || y.len() != n * w || w > NONLOCAL_BAND_PANEL {
            return Err(HamiltonianError::Shape);
        }
        if p == 0 {
            return Ok(());
        }
        let (projections, rest) = scratch.split_at_mut(p * NONLOCAL_BAND_PANEL);
        let (weights, image) = rest.split_at_mut(p * NONLOCAL_BAND_PANEL);
        let (projections, weights, image) = (
            &mut projections[..p * w],
            &mut weights[..p * w],
            &mut image[..n * w],
        );
        dft_numerics::linalg::complex_gemm_slices_into(
            &block.packed,
            n,
            p,
            true,
            x,
            w,
            projections,
        )
        .map_err(|_| HamiltonianError::Shape)?;
        for band in 0..w {
            let pk = &projections[band * p..(band + 1) * p];
            let wk = &mut weights[band * p..(band + 1) * p];
            let mut start = 0;
            for atom in &self.atoms {
                let nc = atom.columns.len();
                for i in 0..nc {
                    wk[start + i] = (0..nc)
                        .map(|j| atom.coupling[i * nc + j] * pk[start + j])
                        .sum();
                }
                start += nc;
            }
        }
        dft_numerics::linalg::complex_gemm_slices_into(
            &block.packed,
            n,
            p,
            false,
            weights,
            w,
            image,
        )
        .map_err(|_| HamiltonianError::Shape)?;
        for (out, v) in y.iter_mut().zip(image.iter()) {
            *out += v;
        }
        Ok(())
    }
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
            // Exact squared-momentum keys only; no shell rounding or interpolation.
            // Flat preallocated storage avoids allocation on cache misses, including
            // generic k where nearly every momentum can have a different norm.
            let ng = basis.momenta().len();
            let mut radial_cache = std::collections::HashMap::with_capacity(ng);
            let mut radial_values = vec![0.0; ng * np];
            for (ig, g) in basis.momenta().iter().enumerate() {
                let next = radial_cache.len();
                let slot = match radial_cache.entry(g.norm_squared.to_bits()) {
                    std::collections::hash_map::Entry::Occupied(entry) => *entry.get(),
                    std::collections::hash_map::Entry::Vacant(entry) => {
                        for i in 0..np {
                            radial_values[next * np + i] =
                                pp.projector_transform(i, g.norm_squared.sqrt())?;
                        }
                        entry.insert(next);
                        next
                    }
                };
                let radial = &radial_values[slot * np..(slot + 1) * np];
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
        Ok(Self::assemble(
            basis.momenta().to_vec(),
            atoms,
            max_channels,
        ))
    }
    pub fn atoms(&self) -> &[AtomProjectors] {
        &self.atoms
    }
    fn add(&self, x: &[Complex64], y: &mut [Complex64], scratch: &mut [Complex64]) {
        let (projections, weights) = scratch.split_at_mut(self.max_channels);
        for atom in &self.atoms {
            let nc = atom.columns.len();
            if nc == 0 {
                continue;
            }
            for (p, col) in projections.iter_mut().zip(&atom.columns) {
                *p = col.iter().zip(x).map(|(b, c)| b.conj() * c).sum();
            }
            for (i, weight) in weights[..nc].iter_mut().enumerate() {
                *weight = (0..nc)
                    .map(|j| atom.coupling[i * nc + j] * projections[j])
                    .sum();
            }
            // Original loop structure (one pass over y per projector channel).
            // A 256-element blocked traversal (P10I) was measured 2.9x slower in
            // this function on dev08 (0.929 s versus 0.323 s, MKL T=1 si8).
            for (col, &weight) in atom.columns.iter().zip(&weights[..nc]) {
                for (out, b) in y.iter_mut().zip(col) {
                    *out += b * weight;
                }
            }
        }
    }
    /// Re <x|V_nl|x> for every vector of `bands`, from the projections alone (no
    /// back-application, no output vector). Up to `EXPECTATION_BLOCK` vectors share
    /// each pass over a projector column, and each vector keeps its own two
    /// accumulators per column, so a value never depends on how the bands are
    /// grouped or on the thread count. Rounding differs from `apply` followed by a
    /// dot product (different summation order); the difference is a few ulp.
    pub fn expectation_values(
        &self,
        bands: &[&[Complex64]],
        out: &mut [f64],
    ) -> Result<(), HamiltonianError> {
        if bands.len() != out.len() || bands.iter().any(|x| x.len() != self.basis.len()) {
            return Err(HamiltonianError::Shape);
        }
        let mut done = 0;
        while done < bands.len() {
            let take = (bands.len() - done).min(EXPECTATION_BLOCK);
            let (group, result) = (&bands[done..done + take], &mut out[done..done + take]);
            match take {
                4 => self.expectation_block::<4>(group, result),
                3 => self.expectation_block::<3>(group, result),
                2 => self.expectation_block::<2>(group, result),
                _ => self.expectation_block::<1>(group, result),
            }
            done += take;
        }
        if out.iter().any(|v| !v.is_finite()) {
            return Err(HamiltonianError::NonFinite);
        }
        Ok(())
    }
    fn expectation_block<const B: usize>(&self, bands: &[&[Complex64]], out: &mut [f64]) {
        let mut energy = [0.0_f64; B];
        let mut projections = vec![[Complex64::default(); B]; self.max_channels];
        for atom in &self.atoms {
            let nc = atom.columns.len();
            for (p, col) in projections.iter_mut().zip(&atom.columns) {
                let mut re = [0.0_f64; B];
                let mut im = [0.0_f64; B];
                let n = col.len();
                let xs: [&[Complex64]; B] = std::array::from_fn(|k| &bands[k][..n]);
                for (g, b) in col.iter().enumerate() {
                    for k in 0..B {
                        let c = xs[k][g];
                        re[k] += b.re * c.re + b.im * c.im;
                        im[k] += b.re * c.im - b.im * c.re;
                    }
                }
                for k in 0..B {
                    p[k] = Complex64::new(re[k], im[k]);
                }
            }
            // E_k += Re sum_ij conj(p_i) D_ij p_j (D is real).
            for i in 0..nc {
                for j in 0..nc {
                    let d = atom.coupling[i * nc + j];
                    if d == 0.0 {
                        continue;
                    }
                    for k in 0..B {
                        energy[k] += d * (projections[i][k].conj() * projections[j][k]).re;
                    }
                }
            }
        }
        out.copy_from_slice(&energy);
    }
    /// Projection storage required by `apply_into`.
    pub fn scratch_len(&self) -> usize {
        2 * self.max_channels
    }
    /// Caller-owned projection scratch; no allocation per orbital.
    pub fn apply_into(
        &self,
        x: &[Complex64],
        y: &mut [Complex64],
        scratch: &mut [Complex64],
    ) -> Result<(), HamiltonianError> {
        if x.len() != self.basis.len() || y.len() != x.len() || scratch.len() != self.scratch_len()
        {
            return Err(HamiltonianError::Shape);
        }
        if !x.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        y.fill(Complex64::default());
        self.add(x, y, scratch);
        if !y.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        Ok(())
    }
    pub fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), HamiltonianError> {
        if x.len() != self.basis.len() || y.len() != x.len() {
            return Err(HamiltonianError::Shape);
        }
        if !x.iter().all(finite) {
            return Err(HamiltonianError::NonFinite);
        }
        y.fill(Complex64::default());
        self.add(x, y, &mut vec![Complex64::default(); 2 * self.max_channels]);
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
    scratch: std::sync::Mutex<Scratch>,
}
struct Scratch {
    fft: FusedPlan,
    projections: Vec<Complex64>,
    /// P48C band-panel scratch (empty when the block path is off).
    block: Vec<Complex64>,
}
/// Immutable physical data with interior scratch reuse; intentionally not Sync.
/// Local grid values are in Ha and must use the internally derived grid dimensions.
pub struct Hamiltonian {
    identity: Arc<()>,
    data: HamiltonianData,
    scratch: RefCell<Vec<Worker>>,
    pool: Arc<WorkerPool>,
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
            scratch: std::sync::Mutex::new(Scratch {
                fft: FusedPlan::new(dims, grid.g_to_fft())?,
                projections: vec![Complex64::default(); 2 * nonlocal.max_channels],
                block: vec![Complex64::default(); nonlocal.block_scratch_len()],
            }),
        }]);
        Ok(Self {
            identity: Arc::new(()),
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
            pool: WorkerPool::new(1).map_err(|_| HamiltonianError::Busy)?,
        })
    }
    /// Set the intra-k budget once during construction. Zero means one.
    /// Numerical buffers and persistent workers are reused.
    pub fn with_band_threads(self, threads: usize) -> Result<Self, HamiltonianError> {
        let pool = WorkerPool::new(threads).map_err(|_| HamiltonianError::Busy)?;
        self.with_band_pool(threads, pool)
    }
    /// Share the SCF executor; the band budget limits tasks, not pool size.
    pub fn with_band_pool(
        mut self,
        threads: usize,
        executor: Arc<WorkerPool>,
    ) -> Result<Self, HamiltonianError> {
        self.pool = executor;
        let dims = self.data.grid.dims();
        let pool = self.scratch.get_mut();
        for _ in pool.len()..threads.max(1) {
            pool.push(Worker {
                scratch: std::sync::Mutex::new(Scratch {
                    fft: FusedPlan::new(dims, self.data.grid.g_to_fft())?,
                    projections: vec![Complex64::default(); 2 * self.data.nonlocal.max_channels],
                    block: vec![Complex64::default(); self.data.nonlocal.block_scratch_len()],
                }),
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
    /// Copy local samples into retained storage, without replacing its allocation.
    pub fn copy_local_potential(
        &mut self,
        potential: &Array3<f64>,
    ) -> Result<(), HamiltonianError> {
        if potential.shape() != self.data.potential.shape() {
            return Err(HamiltonianError::Shape);
        }
        if potential.iter().any(|v| !v.is_finite()) {
            return Err(HamiltonianError::NonFinite);
        }
        self.data.potential.assign(potential);
        Ok(())
    }
    /// Stable identity of the immutable cell/basis/kinetic/projector data.
    /// Potential updates preserve it; constructing another operator creates a
    /// distinct token. Retaining the token prevents address-reuse ambiguity.
    pub fn operator_identity(&self) -> Arc<()> {
        self.identity.clone()
    }
    /// Read-only basis identity for isolated operator adapters.
    pub fn basis_momenta(&self) -> &[GVector] {
        &self.data.nonlocal.basis
    }
    /// Read-only FFT mapping for isolated operator adapters.
    pub fn fft_grid(&self) -> &FftGridSpec {
        &self.data.grid
    }
    /// Real potential in standard contiguous grid order, in Hartree.
    pub fn local_potential_values(&self) -> &[f64] {
        self.data
            .potential
            .as_slice()
            .expect("Hamiltonian stores standard layout")
    }
    /// P48C test and experiment hook: enable or disable the block nonlocal
    /// path regardless of `BAANDI_NONLOCAL_BLOCK` (and of its auto gate), resizing every worker's
    /// panel scratch (allocation here, never during application).
    pub fn set_nonlocal_block(&mut self, on: bool) -> Result<(), HamiltonianError> {
        self.data.nonlocal.set_block(on);
        let len = self.data.nonlocal.block_scratch_len();
        for worker in self.scratch.get_mut().iter_mut() {
            let scratch = worker
                .scratch
                .get_mut()
                .map_err(|_| HamiltonianError::Busy)?;
            scratch.block = vec![Complex64::default(); len];
        }
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
    /// `nonlocal = false` leaves the nonlocal term to the P48C block phase.
    fn apply_with(
        &self,
        x: &[Complex64],
        y: &mut [Complex64],
        scratch: &mut Scratch,
        nonlocal: bool,
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
        if nonlocal {
            self.nonlocal.add(x, y, &mut scratch.projections);
        }
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
    fn worker_pool(&self) -> Option<Arc<WorkerPool>> {
        Some(self.pool.clone())
    }
    fn inner_threads(&self) -> usize {
        self.scratch.borrow().len()
    }
    fn dimension(&self) -> usize {
        self.data.kinetic.len()
    }
    fn apply(&self, x: &[Complex64], y: &mut [Complex64]) -> Result<(), SolverError> {
        let mut scratch = self
            .scratch
            .try_borrow_mut()
            .map_err(|_| application(HamiltonianError::Busy))?;
        self.data
            .apply_with(
                x,
                y,
                scratch[0]
                    .scratch
                    .get_mut()
                    .map_err(|_| application(HamiltonianError::Busy))?,
                true,
            )
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
        // P48C: blocks of two or more bands take the nonlocal term out of the
        // per-band phase and apply it once per fixed band panel afterwards.
        // Single-column calls keep the per-band path bit for bit.
        let block = self.data.nonlocal.block_enabled() && x.cols() >= 2;
        {
            let _region = Span::region("hamiltonian", scratch.len().min(x.cols()).max(1));
            self.local_phase(&mut scratch, x, y, !block)?;
        }
        if block {
            self.nonlocal_block(&mut scratch, x, y)?;
        }
        Ok(())
    }
}

impl Hamiltonian {
    /// Kinetic + local (and, with `nonlocal`, the per-band nonlocal) term of
    /// every column, distributed over the worker slots exactly as before P48C.
    fn local_phase(
        &self,
        scratch: &mut [Worker],
        x: &ComplexMatrix,
        y: &mut ComplexMatrix,
        nonlocal: bool,
    ) -> Result<(), SolverError> {
        if scratch.len() == 1 || x.cols() <= 1 {
            for j in 0..x.cols() {
                let _piece = Piece::new();
                self.data
                    .apply_with(
                        x.column(j),
                        y.column_mut(j),
                        scratch[0]
                            .scratch
                            .get_mut()
                            .map_err(|_| application(HamiltonianError::Busy))?,
                        nonlocal,
                    )
                    .map_err(|e| SolverError::Application(format!("band index {j}: {e}")))?;
            }
            return Ok(());
        }
        let output = y.as_mut_slice();
        if scratch.len() < self.pool.workers_created() || x.cols() <= scratch.len() {
            // With outer k parallelism, use disjoint static scratch. A modulo
            // mapping of the entire pool onto a smaller inner budget causes
            // lock collisions even when another scratch slot is free.
            let count = scratch.len().min(x.cols());
            let data = &self.data;
            let scratch = &mut scratch[..count];
            return self
                .pool
                .install(|| apply_static_workers(data, scratch, x, output, 0, x.cols(), nonlocal));
        }
        let scratch = &scratch[..];
        let data = &self.data;
        self.pool
            .install(|| apply_workers(data, scratch, x, output, 0, x.cols(), nonlocal))
    }

    /// P48C block phase: Y += B D B^H X in fixed panels of
    /// `NONLOCAL_BAND_PANEL` bands. Panels are split statically over the worker
    /// slots; a panel's arithmetic never depends on which slot runs it, so the
    /// result is identical at every thread budget.
    fn nonlocal_block(
        &self,
        scratch: &mut [Worker],
        x: &ComplexMatrix,
        y: &mut ComplexMatrix,
    ) -> Result<(), SolverError> {
        let (n, m) = (x.rows(), x.cols());
        let panels = m.div_ceil(NONLOCAL_BAND_PANEL);
        let count = scratch.len().min(panels).max(1);
        let _region = Span::region("nonlocal_block", count);
        let data = &self.data;
        let workers = &mut scratch[..count];
        let input = x.as_slice();
        let output = y.as_mut_slice();
        if count == 1 {
            nonlocal_panels(data, workers, input, output, n, m, 0, panels)
        } else {
            self.pool
                .install(|| nonlocal_panels(data, workers, input, output, n, m, 0, panels))
        }
        .map_err(application)?;
        for j in 0..m {
            if !y.column(j).iter().all(finite) {
                return Err(SolverError::Application(format!(
                    "band index {j}: {}",
                    HamiltonianError::NonFinite
                )));
            }
        }
        Ok(())
    }
}

// P48C: panels [first, end) of the band block; `y` holds exactly their columns.
#[allow(clippy::too_many_arguments)] // block operands plus a fixed panel range
fn nonlocal_panels(
    data: &HamiltonianData,
    workers: &mut [Worker],
    x: &[Complex64],
    y: &mut [Complex64],
    n: usize,
    m: usize,
    first: usize,
    end: usize,
) -> Result<(), HamiltonianError> {
    if workers.len() == 1 {
        let scratch = workers[0]
            .scratch
            .get_mut()
            .map_err(|_| HamiltonianError::Busy)?;
        let mut offset = 0;
        for panel in first..end {
            let _piece = Piece::new();
            let column = panel * NONLOCAL_BAND_PANEL;
            let w = NONLOCAL_BAND_PANEL.min(m - column);
            data.nonlocal.add_block_panel(
                &x[column * n..(column + w) * n],
                &mut y[offset..offset + w * n],
                w,
                &mut scratch.block,
            )?;
            offset += w * n;
        }
        return Ok(());
    }
    let left_count = workers.len() / 2;
    let middle = first + (end - first) * left_count / workers.len();
    let (left, right) = workers.split_at_mut(left_count);
    let (left_output, right_output) = y.split_at_mut((middle - first) * NONLOCAL_BAND_PANEL * n);
    let (left, right) = WorkerPool::join(
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                nonlocal_panels(data, left, x, left_output, n, m, first, middle)
            }))
        },
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                nonlocal_panels(data, right, x, right_output, n, m, middle, end)
            }))
        },
    );
    left.map_err(|_| HamiltonianError::Busy)??;
    right.map_err(|_| HamiltonianError::Busy)?
}

// Disjoint column-major output slices need neither a mutex nor a copy.
// Joins precede error propagation in ascending band order.
fn apply_workers(
    data: &HamiltonianData,
    workers: &[Worker],
    x: &ComplexMatrix,
    y: &mut [Complex64],
    first: usize,
    end: usize,
    nonlocal: bool,
) -> Result<(), SolverError> {
    // One fixed band per piece. Scratch belongs to a worker slot and is
    // allocated at Hamiltonian construction, never per band. No nested work
    // runs while holding a slot. Modulo bounds concurrency under outer-k work.
    if end - first == 1 {
        let slot = WorkerPool::current_worker_index().unwrap_or(0) % workers.len();
        let mut scratch = workers[slot]
            .scratch
            .lock()
            .map_err(|_| application(HamiltonianError::Busy))?;
        let _piece = Piece::new();
        return data
            .apply_with(x.column(first), y, &mut scratch, nonlocal)
            .map_err(|e| SolverError::Application(format!("band index {first}: {e}")));
    }
    let middle = first + (end - first) / 2;
    let (left_output, right_output) = y.split_at_mut((middle - first) * x.rows());
    let (left, right) = WorkerPool::join(
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                apply_workers(data, workers, x, left_output, first, middle, nonlocal)
            }))
        },
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                apply_workers(data, workers, x, right_output, middle, end, nonlocal)
            }))
        },
    );
    left.map_err(|_| SolverError::Application(format!("band index {first}: worker panicked")))??;
    right.map_err(|_| SolverError::Application(format!("band index {middle}: worker panicked")))?
}

fn apply_static_workers(
    data: &HamiltonianData,
    workers: &mut [Worker],
    x: &ComplexMatrix,
    y: &mut [Complex64],
    first: usize,
    end: usize,
    nonlocal: bool,
) -> Result<(), SolverError> {
    if workers.len() == 1 {
        let _piece = Piece::new();
        let scratch = workers[0]
            .scratch
            .get_mut()
            .map_err(|_| application(HamiltonianError::Busy))?;
        for j in first..end {
            let offset = (j - first) * x.rows();
            data.apply_with(
                x.column(j),
                &mut y[offset..offset + x.rows()],
                scratch,
                nonlocal,
            )
            .map_err(|e| SolverError::Application(format!("band index {j}: {e}")))?;
        }
        return Ok(());
    }
    let left_count = workers.len() / 2;
    let middle = first + (end - first) * left_count / workers.len();
    let (left, right) = workers.split_at_mut(left_count);
    let (left_output, right_output) = y.split_at_mut((middle - first) * x.rows());
    let (left, right) = WorkerPool::join(
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                apply_static_workers(data, left, x, left_output, first, middle, nonlocal)
            }))
        },
        || {
            std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                apply_static_workers(data, right, x, right_output, middle, end, nonlocal)
            }))
        },
    );
    left.map_err(|_| SolverError::Application(format!("band index {first}: worker panicked")))??;
    right.map_err(|_| SolverError::Application(format!("band index {middle}: worker panicked")))?
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
                    let scratch = w.scratch.lock().unwrap();
                    (
                        &scratch.fft as *const FusedPlan as usize,
                        scratch.projections.as_ptr() as usize,
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
                projections: vec![Complex64::default(); 2 * h.data.nonlocal.max_channels],
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
            block: None,
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
    // Independent pre-P10J traversal oracle, retained for bitwise comparison.
    fn legacy_nonlocal(
        op: &NonlocalOperator,
        x: &[Complex64],
        y: &mut [Complex64],
        p: &mut [Complex64],
    ) {
        for atom in &op.atoms {
            let nc = atom.columns.len();
            for (p, col) in p.iter_mut().zip(&atom.columns) {
                *p = col.iter().zip(x).map(|(b, c)| b.conj() * c).sum();
            }
            for (i, col) in atom.columns.iter().enumerate() {
                let weight: Complex64 = (0..nc).map(|j| atom.coupling[i * nc + j] * p[j]).sum();
                for (out, b) in y.iter_mut().zip(col) {
                    *out += b * weight;
                }
            }
        }
    }

    fn nonlocal_fixture(atoms: usize, rows: usize, channels: usize) -> NonlocalOperator {
        NonlocalOperator {
            basis: Vec::new(),
            max_channels: channels,
            block: None,
            atoms: (0..atoms)
                .map(|a| AtomProjectors {
                    columns: (0..channels)
                        .map(|c| vector(rows, a * channels + c))
                        .collect(),
                    coupling: (0..channels * channels)
                        .map(|i| if i % (channels + 1) == 0 { 0.7 } else { 0.03 })
                        .collect(),
                })
                .collect(),
        }
    }

    #[test]
    fn nonlocal_traversal_preserves_bits() {
        for (atoms, rows, channels) in [
            (4, 257, 1),
            (8, 513, 3),
            (8, 1537, 4),
            (20, 769, 9),
            (3, 259, 0),
        ] {
            let op = nonlocal_fixture(atoms, rows, channels);
            let x = vector(rows, 71);
            let mut a = vector(rows, 31);
            let mut b = a.clone();
            let mut scratch = vec![Complex64::default(); 2 * channels];
            legacy_nonlocal(&op, &x, &mut a, &mut scratch);
            op.add(&x, &mut b, &mut scratch);
            assert_eq!(bits(&a), bits(&b));
        }
    }

    #[test]
    fn expectation_blocks_match_apply_and_ignore_grouping() {
        for (atoms, rows, channels) in [(4, 257, 1), (8, 513, 3), (5, 1537, 4), (3, 259, 0)] {
            let op = nonlocal_fixture(atoms, rows, channels);
            let owned: Vec<_> = (0..7).map(|b| vector(rows, 100 + b)).collect();
            let bands: Vec<&[Complex64]> = owned.iter().map(|v| v.as_slice()).collect();
            // Reference: back-application followed by the dot product.
            let mut scratch = vec![Complex64::default(); 2 * channels];
            let reference: Vec<f64> = bands
                .iter()
                .map(|x| {
                    let mut y = vec![Complex64::default(); rows];
                    op.add(x, &mut y, &mut scratch);
                    x.iter().zip(&y).map(|(c, v)| (c.conj() * v).re).sum()
                })
                .collect();
            let mut alone = [0.0; 7];
            for b in 0..7 {
                op.expectation_block::<1>(&bands[b..b + 1], &mut alone[b..b + 1]);
            }
            let mut grouped = [0.0; 7];
            op.expectation_block::<4>(&bands[..4], &mut grouped[..4]);
            op.expectation_block::<3>(&bands[4..], &mut grouped[4..]);
            for b in 0..7 {
                // Grouping never changes a band's bits.
                assert_eq!(alone[b].to_bits(), grouped[b].to_bits());
                let scale = reference[b].abs().max(1.0);
                assert!(
                    (alone[b] - reference[b]).abs() <= 1e-12 * scale * (rows as f64).sqrt(),
                    "{atoms} {rows} {channels} band {b}: {} vs {}",
                    alone[b],
                    reference[b]
                );
            }
        }
    }

    #[test]
    #[ignore] // cargo test --release -p dft-physics nonlocal_traversal -- --ignored --nocapture
    fn nonlocal_traversal() {
        println!("atoms,rows,channels,old_us,new_us");
        for (atoms, rows, channels) in [
            (4, 257, 1),
            (2, 769, 4),
            (4, 1025, 4),
            (4, 4097, 4),
            (8, 2969, 4),
            (8, 4097, 3),
            (20, 16385, 9),
        ] {
            let op = nonlocal_fixture(atoms, rows, channels);
            let x = vector(rows, 71);
            let mut a = vector(rows, 31);
            let mut b = a.clone();
            let mut scratch = vec![Complex64::default(); 2 * channels];
            legacy_nonlocal(&op, &x, &mut a, &mut scratch);
            op.add(&x, &mut b, &mut scratch);
            assert_eq!(bits(&a), bits(&b));
            let mut times = [0.0; 2];
            for (mode, time) in times.iter_mut().enumerate() {
                let mut samples = Vec::new();
                for _ in 0..5 {
                    let start = Instant::now();
                    for _ in 0..50 {
                        b.fill(Complex64::default());
                        if mode == 0 {
                            legacy_nonlocal(black_box(&op), &x, &mut b, &mut scratch);
                        } else {
                            black_box(&op).add(&x, &mut b, &mut scratch);
                        }
                        black_box(&b);
                    }
                    samples.push(start.elapsed().as_secs_f64() * 1e6 / 50.0);
                }
                samples.sort_by(f64::total_cmp);
                *time = samples[2];
            }
            println!("{atoms},{rows},{channels},{:.3},{:.3}", times[0], times[1]);
        }
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

pub fn ionic_local_with_conventions(
    density: &PlaneWaveBasis,
    volume_bohr3: f64,
    species: &[Species<'_>],
    local_mesh: dft_upf::radial::LocalMesh,
    rule: dft_upf::radial::RadialEvenRule,
) -> Result<Vec<Complex64>, HamiltonianError> {
    if local_mesh == dft_upf::radial::LocalMesh::Baandi
        && rule == dft_upf::radial::RadialEvenRule::Simpson3_8
    {
        return ionic_local(density, volume_bohr3, species);
    }

    volume(volume_bohr3)?;
    let mut result = vec![Complex64::default(); density.vectors().len()];
    for s in species {
        let local = s
            .potential
            .local_transform_with_conventions(volume_bohr3, local_mesh, rule)?;
        // Keep LocalTransform's endpoint and quadrature, but hoist only the
        // G-independent parenthesized expression. In particular, rab must
        // still multiply AFTER j0: reassociation would change the bits.
        let pp = s.potential;
        let n = pp
            .mesh()
            .r()
            .partition_point(|r| *r <= local.effective_r_max());
        let mesh = dft_upf::radial::RadialMesh::new(
            pp.mesh().r()[..n].to_vec(),
            pp.mesh().rab()[..n].to_vec(),
        )?;
        let short: Vec<_> = mesh
            .r()
            .iter()
            .zip(pp.local_hartree())
            .map(|(&r, &v)| r * r * v + pp.z_valence() * r * libm::erf(r))
            .collect();
        let mut integrand = vec![0.0; n];
        // At most one entry per G; reserve up front so misses never allocate.
        let mut cache = std::collections::HashMap::with_capacity(density.vectors().len());
        for (out, g) in result.iter_mut().zip(density.vectors()) {
            let v = if g.miller == [0; 3] {
                local.alpha_z()?
            } else {
                let key = g.norm_squared.to_bits();
                if let Some(&v) = cache.get(&key) {
                    v
                } else {
                    let q = g.norm_squared.sqrt();
                    if !q.is_finite() || q <= 0.0 || !(q * local.effective_r_max()).is_finite() {
                        return Err(dft_upf::radial::RadialError::Argument(
                            "G must be finite and positive",
                        )
                        .into());
                    }
                    for ((value, &base), &r) in integrand.iter_mut().zip(&short).zip(mesh.r()) {
                        *value = base * dft_upf::radial::spherical_bessel(0, q * r)?;
                    }
                    let v = 4.0 * std::f64::consts::PI / volume_bohr3
                        * (mesh.integrate_with_rule(&integrand, rule)?
                            - pp.z_valence() * (-q * q / 4.0).exp() / (q * q));
                    if !v.is_finite() {
                        return Err(dft_upf::radial::RadialError::Argument(
                            "unrepresentable result",
                        )
                        .into());
                    }
                    cache.insert(key, v);
                    v
                }
            };
            *out += v * structure_factor(g.cartesian, s.positions_bohr)?;
        }
    }
    if !result.iter().all(finite) {
        return Err(HamiltonianError::NonFinite);
    }
    Ok(result)
}

impl NonlocalOperator {
    pub fn new_with_rule(
        basis: &PlaneWaveBasis,
        volume_bohr3: f64,
        species: &[Species<'_>],
        rule: dft_upf::radial::RadialEvenRule,
    ) -> Result<Self, HamiltonianError> {
        if rule == dft_upf::radial::RadialEvenRule::Simpson3_8 {
            return Self::new(basis, volume_bohr3, species);
        }

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
            // Exact squared-momentum keys only; no shell rounding or interpolation.
            // Flat preallocated storage avoids allocation on cache misses, including
            // generic k where nearly every momentum can have a different norm.
            let ng = basis.momenta().len();
            let mut radial_cache = std::collections::HashMap::with_capacity(ng);
            let mut radial_values = vec![0.0; ng * np];
            for (ig, g) in basis.momenta().iter().enumerate() {
                let next = radial_cache.len();
                let slot = match radial_cache.entry(g.norm_squared.to_bits()) {
                    std::collections::hash_map::Entry::Occupied(entry) => *entry.get(),
                    std::collections::hash_map::Entry::Vacant(entry) => {
                        for i in 0..np {
                            radial_values[next * np + i] =
                                pp.projector_transform_with_rule(i, g.norm_squared.sqrt(), rule)?;
                        }
                        entry.insert(next);
                        next
                    }
                };
                let radial = &radial_values[slot * np..(slot + 1) * np];
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
        Ok(Self::assemble(
            basis.momenta().to_vec(),
            atoms,
            max_channels,
        ))
    }
}

#[cfg(test)]
mod block_tests {
    //! P48C: block nonlocal (two GEMMs per fixed band panel) against the
    //! per-band path; thread-count independence; width-1 bit identity.
    use super::*;

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
    /// Skewed cell at a general k; atoms with 4, 0, 9 and 1 channels and
    /// off-diagonal couplings inside each atom (mixed projector structure).
    fn hamiltonian(threads: usize, block: bool) -> Hamiltonian {
        let cell =
            Cell::from_user_rows_bohr([[5., 0., 0.], [0.7, 6., 0.], [0.2, 0.3, 7.]]).unwrap();
        let basis =
            PlaneWaveBasis::at_k(&cell, dft_units::Hartree(4.), [0.25, -0.25, 0.125]).unwrap();
        let n = basis.vectors().len();
        let channels = [4usize, 0, 9, 1];
        let atoms = channels
            .iter()
            .enumerate()
            .map(|(a, &nc)| AtomProjectors {
                columns: (0..nc)
                    .map(|c| vector(n, 7 * a + c).into_iter().map(|z| z * 0.05).collect())
                    .collect(),
                coupling: (0..nc * nc)
                    .map(|i| {
                        if i / nc == i % nc {
                            0.7 - 0.1 * a as f64
                        } else if (i / nc + i % nc) % 3 == 0 {
                            0.03
                        } else {
                            0.0
                        }
                    })
                    .collect(),
            })
            .collect();
        let nonlocal = NonlocalOperator {
            basis: basis.momenta().to_vec(),
            atoms,
            max_channels: 9,
            block: None,
        };
        let grid = FftGridSpec::new(&basis, &basis.density_basis(&cell).unwrap(), None).unwrap();
        let potential = Array3::from_shape_vec(
            grid.dims(),
            vector(grid.dims().iter().product(), 4)
                .iter()
                .map(|z| 0.1 * z.re)
                .collect(),
        )
        .unwrap();
        let mut h = Hamiltonian::new(&cell, &basis, potential, nonlocal)
            .unwrap()
            .with_band_threads(threads)
            .unwrap();
        h.set_nonlocal_block(block).unwrap();
        h
    }
    fn block(n: usize, width: usize) -> ComplexMatrix {
        let mut x = ComplexMatrix::zeros(n, width).unwrap();
        for j in 0..width {
            x.column_mut(j).copy_from_slice(&vector(n, 1000 + 13 * j));
        }
        x
    }

    #[test]
    fn block_nonlocal_matches_per_band_and_is_thread_independent() {
        let reference = hamiltonian(1, false);
        let n = reference.dimension();
        let mut worst = 0.0_f64;
        for width in [1, 2, 7, 33, 64] {
            let x = block(n, width);
            let mut expected = ComplexMatrix::zeros(n, width).unwrap();
            reference.apply_block(&x, &mut expected).unwrap();
            // Nonlocal image alone (per band), the scale of the rounding change.
            let mut image = vec![Complex64::default(); n];
            let mut first: Option<Vec<(u64, u64)>> = None;
            for threads in [1, 2, 4, 8] {
                let h = hamiltonian(threads, true);
                assert!(h.nonlocal().block_enabled());
                let mut y = ComplexMatrix::zeros(n, width).unwrap();
                h.apply_block(&x, &mut y).unwrap();
                let b = bits(y.as_slice());
                match &first {
                    None => first = Some(b),
                    Some(f) => assert_eq!(f, &b, "width {width} T={threads}"),
                }
                if width == 1 {
                    assert_eq!(bits(expected.as_slice()), bits(y.as_slice()));
                }
                for j in 0..width {
                    reference.nonlocal().apply(x.column(j), &mut image).unwrap();
                    let scale = image.iter().map(|z| z.norm()).fold(0.0, f64::max);
                    for (a, e) in y.column(j).iter().zip(expected.column(j)) {
                        let relative = (a - e).norm() / scale;
                        worst = worst.max(relative);
                        assert!(
                            relative <= 1e-14,
                            "width {width} T={threads} band {j}: {relative:e}"
                        );
                    }
                }
            }
        }
        println!("P48C block vs per-band: max |dy| / max|V_nl x| = {worst:e}");
    }

    #[test]
    fn block_nonlocal_is_relative_per_entry() {
        // Per entry against |y_i| as well (y = T + V_loc + V_nl applied).
        let reference = hamiltonian(1, false);
        let h = hamiltonian(4, true);
        let n = reference.dimension();
        let x = block(n, 7);
        let mut expected = ComplexMatrix::zeros(n, 7).unwrap();
        let mut y = ComplexMatrix::zeros(n, 7).unwrap();
        reference.apply_block(&x, &mut expected).unwrap();
        h.apply_block(&x, &mut y).unwrap();
        let mut worst = 0.0_f64;
        for (a, e) in y.as_slice().iter().zip(expected.as_slice()) {
            if e.norm() > 0.0 {
                worst = worst.max((a - e).norm() / e.norm());
            }
        }
        println!("P48C per-entry relative max = {worst:e}");
        assert!(worst <= 1e-14, "{worst:e}");
    }

    #[test]
    fn block_nonlocal_rejects_nonfinite_and_reports_memory() {
        let h = hamiltonian(2, true);
        let n = h.dimension();
        let mut x = block(n, 5);
        let mut y = ComplexMatrix::zeros(n, 5).unwrap();
        x.column_mut(3)[2].re = f64::NAN;
        assert!(h.apply_block(&x, &mut y).is_err());
        let off = hamiltonian(2, false);
        assert!(!off.nonlocal().block_enabled());
        assert_eq!(off.nonlocal().block_bytes(), 0);
        let p: usize = h.nonlocal().atoms().iter().map(|a| a.columns().len()).sum();
        assert_eq!(h.nonlocal().block_bytes(), n * p * 16);
    }

    /// Operator with a dummy basis of `n` vectors and `atoms` x `channels`
    /// projector columns, packed copy not built.
    fn shaped(n: usize, atoms: usize, channels: usize) -> NonlocalOperator {
        NonlocalOperator {
            basis: vec![
                GVector {
                    miller: [0; 3],
                    cartesian: [0.0; 3],
                    norm_squared: 0.0,
                };
                n
            ],
            max_channels: channels,
            block: None,
            atoms: (0..atoms)
                .map(|a| AtomProjectors {
                    columns: (0..channels).map(|c| vector(n, a * channels + c)).collect(),
                    coupling: vec![0.0; channels * channels],
                })
                .collect(),
        }
    }

    #[test]
    fn auto_gate_depends_only_on_shapes() {
        // 256 KiB = 16384 complex entries: n p = 16384 is the first gated shape.
        assert_eq!(NONLOCAL_BLOCK_GATE_BYTES, 262_144);
        assert!(nonlocal_block_gate(2048, 8));
        assert!(!nonlocal_block_gate(2047, 8));
        assert!(!nonlocal_block_gate(749, 8)); // si2, 30 Ry: 95,872 B
        assert!(nonlocal_block_gate(1495, 21)); // MgO, 80 Ry: 502,320 B
        assert!(!nonlocal_block_gate(0, 0));
        for (n, atoms, channels, expected) in [(2048, 2, 4, true), (2047, 2, 4, false)] {
            let mut op = shaped(n, atoms, channels);
            for mode in ["auto", "off", "on", "auto"] {
                op.set_block_mode(mode);
                let want = match mode {
                    "on" => true,
                    "off" => false,
                    _ => expected,
                };
                assert_eq!(op.block_enabled(), want, "n {n} mode {mode}");
            }
        }
        // The decision is made on the operator before any thread budget
        // exists and is unchanged by the Hamiltonian's band threads.
        for threads in [1, 2, 4, 8] {
            let mut h = hamiltonian(threads, false);
            let n = h.dimension();
            let p = h.nonlocal().projector_count();
            let expected = nonlocal_block_gate(n, p);
            assert!(!expected, "fixture must be below the gate (n {n}, p {p})");
            h.data.nonlocal.set_block_mode("auto");
            assert_eq!(h.nonlocal().block_enabled(), expected, "T={threads}");
        }
    }

    #[test]
    fn below_gate_auto_is_bit_identical_to_off() {
        for threads in [1, 2, 4, 8] {
            let off = hamiltonian(threads, false);
            let mut auto = hamiltonian(threads, false);
            auto.data.nonlocal.set_block_mode("auto");
            assert!(!auto.nonlocal().block_enabled());
            assert_eq!(auto.nonlocal().block_bytes(), 0);
            let n = off.dimension();
            for width in [1, 2, 7, 33, 64] {
                let x = block(n, width);
                let mut expected = ComplexMatrix::zeros(n, width).unwrap();
                let mut y = ComplexMatrix::zeros(n, width).unwrap();
                off.apply_block(&x, &mut expected).unwrap();
                auto.apply_block(&x, &mut y).unwrap();
                assert_eq!(
                    bits(expected.as_slice()),
                    bits(y.as_slice()),
                    "width {width} T={threads}"
                );
            }
        }
    }

    /// Mac observation only (dev08 decides): nonlocal of m bands, per-band
    /// `add` against the block panels; GFLOP/s of the block path counts
    /// 16 n p m real flops (two complex GEMMs).
    #[test]
    #[ignore = "cargo test --release -p dft-physics --lib nonlocal_block_microbench -- --ignored --nocapture"]
    fn nonlocal_block_microbench() {
        use std::{hint::black_box, time::Instant};
        println!("label,n,p,m,per_band_us,block_us,ratio,block_gflops,per_band_gflops");
        for (label, atoms, channels, n) in [
            ("si2-k30Ry", 2, 4, 749),
            ("mgo-40-class", 2, 13, 1495),
            ("si16-class", 16, 4, 5990),
            ("si64-class", 64, 4, 23960),
        ] {
            let mut op = NonlocalOperator {
                basis: vec![
                    GVector {
                        miller: [0; 3],
                        cartesian: [0.0; 3],
                        norm_squared: 0.0,
                    };
                    n
                ],
                max_channels: channels,
                block: None,
                atoms: (0..atoms)
                    .map(|a| AtomProjectors {
                        columns: (0..channels).map(|c| vector(n, a * channels + c)).collect(),
                        coupling: (0..channels * channels)
                            .map(|i| if i % (channels + 1) == 0 { 0.7 } else { 0.03 })
                            .collect(),
                    })
                    .collect(),
            };
            op.set_block(true);
            let p = atoms * channels;
            let mut scratch = vec![Complex64::default(); op.block_scratch_len()];
            let mut small = vec![Complex64::default(); 2 * channels];
            for m in [4usize, 8, 16, 32] {
                let x: Vec<Complex64> = (0..m).flat_map(|j| vector(n, 500 + j)).collect();
                let mut y = vec![Complex64::default(); n * m];
                let reps = (2e8 / (16.0 * (n * p * m) as f64)).ceil() as usize;
                let mut times = [0.0; 2];
                for (mode, time) in times.iter_mut().enumerate() {
                    let mut samples = Vec::new();
                    for _ in 0..5 {
                        let start = Instant::now();
                        for _ in 0..reps {
                            if mode == 0 {
                                for j in 0..m {
                                    black_box(&op).add(
                                        &x[j * n..(j + 1) * n],
                                        &mut y[j * n..(j + 1) * n],
                                        &mut small,
                                    );
                                }
                            } else {
                                black_box(&op)
                                    .add_block_panel(&x, &mut y, m, &mut scratch)
                                    .unwrap();
                            }
                            black_box(&y);
                        }
                        samples.push(start.elapsed().as_secs_f64() * 1e6 / reps as f64);
                    }
                    samples.sort_by(f64::total_cmp);
                    *time = samples[2];
                }
                let flops = 16.0 * (n * p * m) as f64;
                println!(
                    "{label},{n},{p},{m},{:.2},{:.2},{:.2},{:.2},{:.2}",
                    times[0],
                    times[1],
                    times[0] / times[1],
                    flops / times[1] / 1e3,
                    flops / times[0] / 1e3
                );
            }
        }
    }
}
