//! Post-SCF scalar fields. Atomic units; no SCF iteration or parallel reduction.
use crate::{PhysicsError, data::*, orbitals, volume};
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_geometry::Cell;
use std::f64::consts::PI;

/// Vacuum rule: ELF is zero at n <= 1e-10 electrons/Bohr^3.
pub const ELF_DENSITY_THRESHOLD: f64 = 1e-10;

/// One k point, supplied in the final SCF's canonical order. Occupations already
/// include spin (usually 2); weight is applied exactly once.
pub struct OrbitalBlock<'a> {
    pub basis: &'a PlaneWaveBasis,
    pub coefficients: &'a [Vec<Complex64>],
    pub occupations: &'a [f64],
    pub weight: f64,
}

/// tau = 1/2 sum_k,n w_k f_n |grad psi_nk|^2, in Ha/Bohr^3.
/// psi = exp(ik.r) sum_G c_G exp(iG.r)/sqrt(volume). The common Bloch
/// phase drops out of the norm, so transform i(k+G)c, without exp(ik.r).
/// Invariant under constant orbital phases, not arbitrary position-dependent
/// electromagnetic gauge transformations (no current correction is included).
pub fn kinetic_density(
    blocks: &[OrbitalBlock<'_>],
    cell: &Cell,
    dims: [usize; 3],
) -> Result<Array3<f64>, PhysicsError> {
    Ok(orbital_fields(blocks, cell, dims)?.0)
}

fn orbital_fields(
    blocks: &[OrbitalBlock<'_>],
    cell: &Cell,
    dims: [usize; 3],
) -> Result<(Array3<f64>, Array3<f64>), PhysicsError> {
    let mut fft = FftPlan3d::new(dims)?;
    let mut psi = Array3::zeros(dims);
    let mut gradients: [Array3<f64>; 3] = std::array::from_fn(|_| Array3::zeros(dims));
    let mut tau = Array3::zeros(dims);
    let mut real = Array3::zeros(dims);
    let density_basis = PlaneWaveBasis::new(
        cell,
        blocks.first().ok_or(PhysicsError::Shape)?.basis.cutoff(),
    )?
    .density_basis(cell)?;
    for block in blocks {
        if !block.weight.is_finite() || block.weight < 0.0 {
            return Err(PhysicsError::Orbitals);
        }
        orbitals(
            block.coefficients,
            block.occupations,
            block.basis.vectors().len(),
        )?;
        let grid = FftGridSpec::new(block.basis, &density_basis, Some(dims))?;
        for (c, &f) in block.coefficients.iter().zip(block.occupations) {
            fft.to_real(&grid.scatter(c)?, &mut psi)?;
            for (axis, gradient) in gradients.iter_mut().enumerate() {
                let derivative: Vec<_> = c
                    .iter()
                    .zip(block.basis.momenta())
                    .map(|(c, g)| *c * Complex64::new(0.0, g.cartesian[axis]))
                    .collect();
                fft.to_real(&grid.scatter(&derivative)?, &mut real)?;
                for ((g, p), d) in gradient.iter_mut().zip(&psi).zip(&real) {
                    *g += 2.0 * block.weight * f * (p.conj() * d).re / cell.volume_bohr3();
                }
                tau.zip_mut_with(&real, |t, z| {
                    *t += 0.5 * block.weight * f * z.norm_sqr() / cell.volume_bohr3()
                });
            }
        }
    }
    if tau.iter().any(|x: &f64| !x.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    let mut gradient_squared = Array3::zeros(dims);
    for g in gradients {
        gradient_squared.zip_mut_with(&g, |s, g| *s += g * g);
    }
    Ok((tau, gradient_squared))
}

/// Unweighted, unoccupied band sum: integral equals the number of normalized
/// bands supplied, including empty bands. No renormalization hides input errors.
pub fn partial_charge(
    coefficients: &[Vec<Complex64>],
    grid: &FftGridSpec,
    volume_bohr3: f64,
) -> Result<Array3<f64>, PhysicsError> {
    volume(volume_bohr3)?;
    crate::density::from_orbitals(
        coefficients,
        &vec![1.0; coefficients.len()],
        grid,
        volume_bohr3,
        &mut FftPlan3d::new(grid.dims())?,
    )
}

/// Spin-summed ELF from density, tau and |grad n|^2. The uniform gas has
/// tau=tau_h, gradient=0 and exactly ELF=1/2. Vacuum is explicitly zero.
pub fn elf_from_fields(
    n: &Array3<f64>,
    tau: &Array3<f64>,
    gradient_squared: &Array3<f64>,
) -> Result<Array3<f64>, PhysicsError> {
    if n.shape() != tau.shape() || n.shape() != gradient_squared.shape() || n.is_empty() {
        return Err(PhysicsError::Shape);
    }
    let mut elf = Array3::zeros(n.raw_dim());
    for (((e, &n), &t), &g) in elf.iter_mut().zip(n).zip(tau).zip(gradient_squared) {
        if [n, t, g].iter().any(|x| !x.is_finite() || *x < 0.0) {
            return Err(PhysicsError::Orbitals);
        }
        if n > ELF_DENSITY_THRESHOLD {
            let chi =
                (t - g / (8.0 * n)) / (0.3 * (3.0 * PI * PI).powf(2.0 / 3.0) * n.powf(5.0 / 3.0));
            *e = 1.0 / (1.0 + chi * chi);
        }
    }
    Ok(elf)
}

/// ELF using the orbital product rule grad n = 2 Re sum(w f psi* grad psi).
/// This avoids loss of relative accuracy from differentiating tiny density tails.
/// The Bloch k term cancels in Re(psi* i k psi). Uses the final SCF density,
/// whose occupations and weights must match the supplied orbital blocks.
pub fn elf_from_orbitals(
    n: &Array3<f64>,
    blocks: &[OrbitalBlock<'_>],
    cell: &Cell,
) -> Result<Array3<f64>, PhysicsError> {
    let dims = [n.shape()[0], n.shape()[1], n.shape()[2]];
    let (tau, g2) = orbital_fields(blocks, cell, dims)?;
    elf_from_fields(n, &tau, &g2)
}

/// Real part of the canonical Hartree field, G=0 exactly zero. The engine's
/// complex skew-grid Nyquist convention is retained until this scalar boundary.
pub fn hartree_potential(n: &Array3<f64>, cell: &Cell) -> Result<Array3<f64>, PhysicsError> {
    let dims = [n.shape()[0], n.shape()[1], n.shape()[2]];
    Ok(
        crate::hartree::evaluate(n, cell, &mut FftPlan3d::new(dims)?)?
            .potential
            .mapv(|z| z.re),
    )
}

/// Local scalar potential in Ha, preserving ionic finite G=0 and XC offset.
/// Nonlocal projectors are not representable as a scalar field.
pub fn total_local(
    ionic: &Array3<f64>,
    hartree: &Array3<f64>,
    xc: &Array3<f64>,
) -> Result<Array3<f64>, PhysicsError> {
    if ionic.shape() != hartree.shape() || ionic.shape() != xc.shape() {
        return Err(PhysicsError::Shape);
    }
    let v = ionic + hartree + xc;
    if v.iter().any(|x| !x.is_finite()) {
        return Err(PhysicsError::Unrepresentable);
    }
    Ok(v)
}
