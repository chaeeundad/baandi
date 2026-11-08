//! Explicit ionic derivatives at fixed cell/basis and fixed electronic state.
//! Arrays follow species order, then each species' position order. The engine
//! must restore input atom order before reporting a total.
use crate::hamiltonian::{HamiltonianError, NonlocalOperator, Species, structure_factor};
use dft_basis::PlaneWaveBasis;
use dft_numerics::linalg::ComplexMatrix;
use num_complex::Complex64;

/// Local or core structure-factor derivative against a conjugate Fourier field.
/// For local terms field is valence density; for core it is V_xc (the mean of
/// the two spin potentials). Fourier fields use forward 1/N normalization.
pub fn local_or_core(
    basis: &PlaneWaveBasis,
    volume: f64,
    species: &[Species<'_>],
    field: &[Complex64],
    core: bool,
) -> Result<Vec<[f64; 3]>, HamiltonianError> {
    if field.len() != basis.vectors().len() {
        return Err(HamiltonianError::Shape);
    }
    if !volume.is_finite() || volume <= 0.0 {
        return Err(HamiltonianError::Volume);
    }
    if field.iter().any(|z| !z.re.is_finite() || !z.im.is_finite()) {
        return Err(HamiltonianError::NonFinite);
    }
    let mut output = Vec::new();
    // Radial values depend only on (potential, |G|): reuse them across atoms of
    // one species and across G vectors with the identical |G| bit pattern.
    let mut per_potential: std::collections::HashMap<usize, Vec<f64>> =
        std::collections::HashMap::new();
    for s in species {
        let potential_id = s.potential as *const _ as *const () as usize;
        // `?` inside the computation prevents the entry API.
        #[allow(clippy::map_entry)]
        if !per_potential.contains_key(&potential_id) {
            let local = s.potential.local_transform(volume, None)?;
            let mut shells: std::collections::HashMap<u64, f64> = std::collections::HashMap::new();
            let radial: Vec<f64> = basis
                .vectors()
                .iter()
                .map(|g| {
                    if g.miller == [0; 3] {
                        return Ok(0.0);
                    }
                    let q = g.norm_squared.sqrt();
                    if let Some(v) = shells.get(&q.to_bits()) {
                        return Ok(*v);
                    }
                    let v = if core {
                        s.potential.core_density().map_or(Ok(0.0), |n| {
                            dft_upf::radial::spherical_density_transform(s.potential.mesh(), n, q)
                                .map(|v| v / volume)
                        })
                    } else {
                        local.at_g(q)
                    }?;
                    shells.insert(q.to_bits(), v);
                    Ok(v)
                })
                .collect::<Result<_, HamiltonianError>>()?;
            per_potential.insert(potential_id, radial);
        }
        let radial = &per_potential[&potential_id];
        for &position in s.positions_bohr {
            let mut force = [0.0; 3];
            for ((g, &v), &field) in basis.vectors().iter().zip(radial).zip(field) {
                let product = field.conj() * v * structure_factor(g.cartesian, &[position])?;
                for (axis, value) in force.iter_mut().enumerate() {
                    *value -= volume * g.cartesian[axis] * product.im;
                }
            }
            output.push(force);
        }
    }
    if output.iter().flatten().any(|v| !v.is_finite()) {
        return Err(HamiltonianError::NonFinite);
    }
    Ok(output)
}

/// One pass through each band's atom/projector coefficients. Occupations already
/// include the normalized k weight and the spin channel capacity. Both sides of
/// the Hermitian KB expectation are differentiated, including off-diagonal D.
pub fn nonlocal(
    operator: &NonlocalOperator,
    basis: &PlaneWaveBasis,
    orbitals: &ComplexMatrix,
    occupations: &[f64],
) -> Result<Vec<[f64; 3]>, HamiltonianError> {
    if orbitals.rows() != basis.momenta().len() || orbitals.cols() != occupations.len() {
        return Err(HamiltonianError::Shape);
    }
    if occupations.iter().any(|f| !f.is_finite() || *f < 0.0) {
        return Err(HamiltonianError::NonFinite);
    }
    let mut forces = vec![[0.0; 3]; operator.atoms().len()];
    for (atom, force) in operator.atoms().iter().zip(&mut forces) {
        let nc = atom.columns().len();
        let mut overlaps = vec![Complex64::default(); nc];
        let mut derivatives = vec![[Complex64::default(); 3]; nc];
        for (band, &occupation) in occupations.iter().enumerate() {
            if occupation == 0.0 {
                continue;
            }
            overlaps.fill(Complex64::default());
            derivatives.fill([Complex64::default(); 3]);
            for (p, column) in atom.columns().iter().enumerate() {
                if column.len() != basis.momenta().len() {
                    return Err(HamiltonianError::Shape);
                }
                for (ig, (beta, g)) in column.iter().zip(basis.momenta()).enumerate() {
                    let z = beta.conj() * orbitals[(ig, band)];
                    overlaps[p] += z;
                    for (axis, derivative) in derivatives[p].iter_mut().enumerate() {
                        *derivative += Complex64::new(0.0, g.cartesian[axis]) * z;
                    }
                }
            }
            for (p, derivative) in derivatives.iter().enumerate() {
                let weighted: Complex64 = (0..nc)
                    .map(|q| atom.coupling_hartree()[p * nc + q] * overlaps[q])
                    .sum();
                for axis in 0..3 {
                    force[axis] -= 2.0 * occupation * (derivative[axis].conj() * weighted).re;
                }
            }
        }
    }
    if forces.iter().flatten().any(|v| !v.is_finite()) {
        return Err(HamiltonianError::NonFinite);
    }
    Ok(forces)
}
