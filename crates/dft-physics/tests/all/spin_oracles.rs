//! MF03 capability gate. Oracle transcribed from the architect-supplied
//! equations of Perdew and Wang, Phys. Rev. B 45, 13244 (1992), in Hartree
//! atomic units. No library evaluation is used inside the oracle.
use std::f64::consts::PI;
use xcx::{Functional, FunctionalId, Spin, XcInput};

// Frozen before measurement: MF03 stage-oracles caps, also applied to PW92
// per the architect's instruction. Energy density cap is n * ENERGY_CAP.
const ENERGY_CAP: f64 = 1e-12;
const POTENTIAL_REL_CAP: f64 = 1e-8;
const POTENTIAL_ABS_CAP: f64 = 1e-10;

// Three central steps as required by stage-oracles. Accept second-order
// convergence or a sub-cap plateau whose step variation is below cap / 3.
fn central_check(energy: impl Fn(f64) -> f64, h: f64, got: f64) -> f64 {
    let slopes = [h, h / 2.0, h / 4.0].map(|step| (energy(step) - energy(-step)) / (2.0 * step));
    let cap = POTENTIAL_ABS_CAP.max(POTENTIAL_REL_CAP * got.abs());
    let errors = slopes.map(|slope| (slope - got).abs());
    let d0 = (slopes[0] - slopes[1]).abs();
    let d1 = (slopes[1] - slopes[2]).abs();
    let trend = d1 <= 0.4 * d0;
    let plateau = errors.iter().all(|&e| e < cap) && d0.max(d1) < cap / 3.0;
    assert!(
        trend || plateau,
        "unresolved finite difference: slopes={slopes:?}, derivative={got:e}, cap={cap:e}"
    );
    let extrapolated = (4.0 * slopes[2] - slopes[1]) / 3.0;
    let error = (extrapolated - got).abs();
    assert!(error <= cap, "derivative error={error:e}, cap={cap:e}");
    error
}

fn exchange(rho: [f64; 2]) -> f64 {
    -0.75 * (6.0 / PI).cbrt() * (rho[0].powf(4.0 / 3.0) + rho[1].powf(4.0 / 3.0))
}

fn correlation(rho: [f64; 2]) -> f64 {
    let n = rho[0] + rho[1];
    let rs = (3.0 / (4.0 * PI * n)).cbrt();
    let z = (rho[0] - rho[1]) / n;
    let g = |[a, a1, b1, b2, b3, b4]: [f64; 6]| {
        -2.0 * a
            * (1.0 + a1 * rs)
            * (1.0
                + 1.0 / (2.0 * a * (b1 * rs.sqrt() + b2 * rs + b3 * rs.powf(1.5) + b4 * rs * rs)))
                .ln()
    };
    let c0 = g([0.031091, 0.21370, 7.5957, 3.5876, 1.6382, 0.49294]);
    let c1 = g([0.015545, 0.20548, 14.1189, 6.1977, 3.3662, 0.62517]);
    let alpha = -g([0.016887, 0.11125, 10.357, 3.6231, 0.88026, 0.49671]);
    let f = ((1.0 + z).powf(4.0 / 3.0) + (1.0 - z).powf(4.0 / 3.0) - 2.0)
        / (2.0_f64.powf(4.0 / 3.0) - 2.0);
    n * (c0 + alpha * f / 1.709921 * (1.0 - z.powi(4)) + (c1 - c0) * f * z.powi(4))
}

#[test]
fn uniform_exchange_and_independent_pw92() {
    let x = Functional::new(FunctionalId::LdaX, Spin::Polarized).unwrap();
    let c = Functional::new(FunctionalId::LdaCPw, Spin::Polarized).unwrap();
    let mut failures = Vec::new();
    let mut max_energy_error: f64 = 0.0;
    let mut max_potential_error: f64 = 0.0;
    for n in [0.01_f64, 0.1, 1.0] {
        for z in [0.0, 0.3, 0.9, 1.0] {
            let rho = [n * (1.0 + z) / 2.0, n * (1.0 - z) / 2.0];
            let xo = x.eval(1, &XcInput::lda(&rho)).unwrap();
            let co = c.eval(1, &XcInput::lda(&rho)).unwrap();
            let error = (n * (xo.exc[0] + co.exc[0]) - exchange(rho) - correlation(rho)).abs() / n;
            max_energy_error = max_energy_error.max(error);
            if error > ENERGY_CAP {
                failures.push(format!("n={n} z={z} XC energy error={error:e} Ha/electron"));
            }
            assert!((n * xo.exc[0] - exchange(rho)).abs() <= n * ENERGY_CAP);
            for s in 0..2 {
                let exact_x = -(6.0 / PI).cbrt() * rho[s].cbrt();
                let error_x = (xo.vrho[s] - exact_x).abs();
                if error_x > POTENTIAL_ABS_CAP.max(POTENTIAL_REL_CAP * exact_x.abs()) {
                    failures.push(format!(
                        "n={n} z={z} spin={s} exchange potential error={error_x:e} Ha"
                    ));
                }
                // Central perturbations are only defined in the interior of the
                // nonnegative-density domain, away from density screens.
                if rho[s] == 0.0 {
                    continue;
                }
                let h = rho[s] * 1e-4;
                let got = xo.vrho[s] + co.vrho[s];
                let error_v = central_check(
                    |step| {
                        let mut perturbed = rho;
                        perturbed[s] += step;
                        exchange(perturbed) + correlation(perturbed)
                    },
                    h,
                    got,
                );
                max_potential_error = max_potential_error.max(error_v);
            }
        }
    }
    println!(
        "PW92+X maximum energy error={max_energy_error:e} Ha/electron; interior potential error={max_potential_error:e} Ha"
    );
    assert!(failures.is_empty(), "{}", failures.join("\n"));
}

#[test]
fn spin_swap_and_unpolarized_reduction() {
    let mut max_error: f64 = 0.0;
    for id in [
        FunctionalId::LdaX,
        FunctionalId::LdaCPw,
        FunctionalId::GgaXPbe,
        FunctionalId::GgaCPbe,
    ] {
        let pol = Functional::new(id, Spin::Polarized).unwrap();
        let unpol = Functional::new(id, Spin::Unpolarized).unwrap();
        for n in [0.01_f64, 0.1, 1.0] {
            for z in [0.0, 0.3, 0.9, 1.0] {
                let rho = [n * (1.0 + z) / 2.0, n * (1.0 - z) / 2.0];
                // Parallel spin gradients of a spatially constant polarization.
                let sigma = [
                    0.04 * rho[0].powi(2),
                    0.04 * rho[0] * rho[1],
                    0.04 * rho[1].powi(2),
                ];
                let a = pol.eval(1, &XcInput::gga(&rho, &sigma)).unwrap();
                let b = pol
                    .eval(
                        1,
                        &XcInput::gga(&[rho[1], rho[0]], &[sigma[2], sigma[1], sigma[0]]),
                    )
                    .unwrap();
                assert!((a.exc[0] - b.exc[0]).abs() <= ENERGY_CAP);
                for s in 0..2 {
                    let error = (a.vrho[s] - b.vrho[1 - s]).abs();
                    max_error = max_error.max(error);
                    assert!(error < POTENTIAL_ABS_CAP);
                }
                for s in 0..a.vsigma.len() {
                    assert!(
                        (a.vsigma[s] - b.vsigma[2 - s]).abs()
                            <= POTENTIAL_ABS_CAP.max(POTENTIAL_REL_CAP * a.vsigma[s].abs())
                    );
                }
                if z == 0.0 {
                    let u = unpol.eval(1, &XcInput::gga(&[n], &[0.04 * n * n])).unwrap();
                    assert!((a.exc[0] - u.exc[0]).abs() <= ENERGY_CAP);
                    assert!((a.vrho[0] - u.vrho[0]).abs() < POTENTIAL_ABS_CAP);
                    if !a.vsigma.is_empty() {
                        let reduced = (a.vsigma[0] + a.vsigma[1] + a.vsigma[2]) / 4.0;
                        assert!(
                            (reduced - u.vsigma[0]).abs()
                                <= POTENTIAL_ABS_CAP.max(POTENTIAL_REL_CAP * u.vsigma[0].abs())
                        );
                    }
                }
            }
        }
    }
    println!(
        "spin swap maximum potential error={max_error:e} Ha; LDA/PBE uniform-polarization zero-spin reduction passed"
    );
}

#[test]
fn pbe_exchange_scaling_and_polarized_derivatives() {
    // Absolute polarized PBE correlation is NOT independently transcribed.
    // This gate checks library energy derivatives and exchange spin scaling.
    let px = Functional::new(FunctionalId::GgaXPbe, Spin::Polarized).unwrap();
    let ux = Functional::new(FunctionalId::GgaXPbe, Spin::Unpolarized).unwrap();
    let mut max_scaling_error: f64 = 0.0;
    let mut max_rho_error: f64 = 0.0;
    let mut max_sigma_error: f64 = 0.0;
    for n in [0.01_f64, 0.1, 1.0] {
        for z in [0.0, 0.3, 0.9, 1.0] {
            let rho = [n * (1.0 + z) / 2.0, n * (1.0 - z) / 2.0];
            let sigma = [
                0.04 * rho[0].powi(2),
                0.01 * rho[0] * rho[1],
                0.09 * rho[1].powi(2),
            ];
            let x = px.eval(1, &XcInput::gga(&rho, &sigma)).unwrap();
            let mut scaled = 0.0;
            for s in 0..2 {
                let u = ux
                    .eval(1, &XcInput::gga(&[2.0 * rho[s]], &[4.0 * sigma[2 * s]]))
                    .unwrap();
                scaled += rho[s] * u.exc[0];
            }
            let error = (n * x.exc[0] - scaled).abs() / n;
            max_scaling_error = max_scaling_error.max(error);
            assert!(error <= ENERGY_CAP);
            if z == 1.0 {
                continue; // Central derivatives require interior channel densities.
            }
            for id in [FunctionalId::GgaXPbe, FunctionalId::GgaCPbe] {
                let f = Functional::new(id, Spin::Polarized).unwrap();
                let out = f.eval(1, &XcInput::gga(&rho, &sigma)).unwrap();
                for s in 0..2 {
                    let error = central_check(
                        |delta| {
                            let mut r = rho;
                            r[s] += delta;
                            let energy = f.eval(1, &XcInput::gga(&r, &sigma)).unwrap().exc[0];
                            (r[0] + r[1]) * energy
                        },
                        1e-4 * rho[s],
                        out.vrho[s],
                    );
                    max_rho_error = max_rho_error.max(error);
                }
                for s in 0..3 {
                    let error = central_check(
                        |delta| {
                            let mut sig = sigma;
                            sig[s] += delta;
                            n * f.eval(1, &XcInput::gga(&rho, &sig)).unwrap().exc[0]
                        },
                        1e-2 * sigma[s],
                        out.vsigma[s],
                    );
                    max_sigma_error = max_sigma_error.max(error);
                }
                // V_up flux = 2 e_uu grad(n_up) + e_ud grad(n_down).
                // Perturb the up gradient's x component: the cross invariant
                // varies once, unlike sigma_uu. Nonparallel gradients ensure
                // this detects a missing or doubled cross term.
                let gu = [0.2 * rho[0], 0.0];
                let gd = [0.05 * rho[1], 0.0875_f64.sqrt() * rho[1]];
                let flux = 2.0 * out.vsigma[0] * gu[0] + out.vsigma[1] * gd[0];
                central_check(
                    |delta| {
                        let g = gu[0] + delta;
                        let sig = [g * g, g * gd[0], gd[0] * gd[0] + gd[1] * gd[1]];
                        n * f.eval(1, &XcInput::gga(&rho, &sig)).unwrap().exc[0]
                    },
                    1e-3 * gu[0],
                    flux,
                );
            }
        }
    }
    println!(
        "PBE exchange scaling error={max_scaling_error:e} Ha/electron; density derivative error={max_rho_error:e} Ha; sigma derivative error={max_sigma_error:e} Ha Bohr^5"
    );
}
