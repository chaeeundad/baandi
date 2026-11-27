//! R2 internal oracles (`docs/roadmap/rpa-crpa-plan.md`, ticket R2) on
//! synthetic band problems (no UPF file, no engine path):
//! (a) Lindhard function of the free-electron gas on the imaginary axis,
//!     against the exact discrete sum and the analytic continuum function;
//! (b) static chi0 columns at q != 0 against finite differences of the
//!     occupied density response of a supercell to a small periodic potential;
//! (c) Hermiticity and negative semidefiniteness on the imaginary axis;
//! (d) the high-frequency limit -(q+G).(q+G') n(G-G') / w^2;
//! (e) k.p head and wings (with the nonlocal velocity) against finite-q
//!     extrapolation; the projector k-gradient against finite differences;
//! plus band-convergence monotonicity and the MPI reduction contract.
//!
//! Fast gate (debug): `cargo test -p dft-response --test oracles`
//! Larger mesh (release, ignored):
//! `cargo test --release -p dft-response --test oracles -- --ignored --nocapture`
mod common;
use common::*;
use dft_basis::PlaneWaveBasis;
use dft_geometry::Cell;
use dft_numerics::linalg::{ComplexMatrix, hermitian_eigh};
use dft_response::{
    Complex64 as C,
    chi0::{Accumulator, Chi0},
    dielectric,
    pair::{PairFft, Plan},
    transitions::{KOrbitals, Transitions, Velocity},
};
use dft_units::Hartree;
use std::collections::HashMap;
use std::f64::consts::PI;

fn rel(a: f64, b: f64) -> f64 {
    (a - b).abs() / b.abs().max(1e-300)
}

/// Response basis at fractional q: Miller, |q+G|, wrapped q.
fn response_basis(cell: &Cell, ecut: f64, q: [f64; 3]) -> PlaneWaveBasis {
    PlaneWaveBasis::at_k(cell, Hartree(ecut), q).unwrap()
}
fn millers(b: &PlaneWaveBasis) -> Vec<[i32; 3]> {
    b.vectors().iter().map(|g| g.miller).collect()
}

/// Options of a chi0 run on a synthetic problem.
struct Run<'a> {
    problem: &'a Problem<'a>,
    mesh: Vec<[f64; 3]>,
    q: [f64; 3],
    ecut_chi: f64,
    omegas: Vec<f64>,
    /// Occupied bands at a k (Bloch vector) from its energies.
    occupied: &'a (dyn Fn(&[f64]) -> usize + Sync),
    /// Highest band (exclusive) used as an empty state; None = complete.
    bands: Option<usize>,
    optical: bool,
    block: usize,
}

/// chi0 by the library path with transitions distributed round robin over
/// `comm` ranks in blocks of `block` empty bands.
fn chi0_run<K: dft_mpi::Collective>(run: &Run<'_>, comm: &K) -> (Chi0, PlaneWaveBasis) {
    let p = run.problem;
    let dims = p.dims();
    let resp = response_basis(&p.cell, run.ecut_chi, run.q);
    let q = resp.k_fractional();
    let mut fft = PairFft::new(dims, &millers(&resp)).unwrap();
    let mut acc = Accumulator::new(resp.vectors().len(), &run.omegas, run.optical, 2).unwrap();
    let scale = 2. / (run.mesh.len() as f64 * p.volume());
    let mut unit = 0;
    for &k in &run.mesh {
        let left = p.solve(k);
        let right = if run.optical {
            None
        } else {
            Some(p.solve([k[0] + q[0], k[1] + q[1], k[2] + q[2]]))
        };
        let r = right.as_ref().unwrap_or(&left);
        let shift = shift(left.basis.k_fractional(), q, r.basis.k_fractional());
        let occ = (run.occupied)(&left.energies);
        let occ_right = (run.occupied)(&r.energies);
        let top = run.bands.unwrap_or(r.energies.len()).min(r.energies.len());
        let lo = KOrbitals {
            miller: &left.miller,
            coefficients: &left.vectors,
            energies: &left.energies,
        };
        let ro = KOrbitals {
            miller: &r.miller,
            coefficients: &r.vectors,
            energies: &r.energies,
        };
        if occ == 0 {
            continue;
        }
        let mut t = Transitions::new(dims, lo, occ, &r.miller, shift).unwrap();
        let mut m = occ_right;
        while m < top {
            let end = (m + run.block).min(top);
            if unit % comm.size() == comm.rank() {
                let v = run.optical.then(|| Velocity {
                    momenta: &left.momenta,
                    projectors: left.projectors.as_ref(),
                });
                t.add(&mut acc, &mut fft, ro, m..end, scale, v).unwrap();
            }
            unit += 1;
            m = end;
        }
    }
    (acc.reduce(comm).unwrap(), resp)
}

// ---------------------------------------------------------------- pair FFT

#[test]
fn pair_densities_equal_the_convolution_sum_with_umklapp() {
    let cell =
        Cell::from_user_rows_bohr([[5.0, 0.0, 0.0], [0.7, 4.6, 0.0], [0.3, 0.4, 5.3]]).unwrap();
    let p = Problem {
        cell,
        cutoff: 3.0,
        local: potential(&[]),
        nonlocal: None,
    };
    let dims = p.dims();
    let k = [0.25, -0.5, 1. / 3.];
    let q = [0.75, 0.5, 1. / 3.]; // k + q wraps in every direction
    let resp = response_basis(&p.cell, 6.0, q);
    let left = PlaneWaveBasis::at_k(&p.cell, Hartree(p.cutoff), k).unwrap();
    let qw = resp.k_fractional();
    let right = PlaneWaveBasis::at_k(
        &p.cell,
        Hartree(p.cutoff),
        [k[0] + qw[0], k[1] + qw[1], k[2] + qw[2]],
    )
    .unwrap();
    let g0 = shift(left.k_fractional(), qw, right.k_fractional());
    assert_ne!(g0, [0; 3], "the fixture must exercise an umklapp");
    let rand = |seed: usize, n: usize, cols: usize| {
        let mut m = ComplexMatrix::zeros(n, cols).unwrap();
        for (i, z) in m.as_mut_slice().iter_mut().enumerate() {
            let x = ((i * 7919 + seed * 104_729) % 1000) as f64 / 1000.;
            let y = ((i * 6271 + seed * 15_485_863) % 997) as f64 / 997.;
            *z = C::new(x - 0.5, y - 0.5);
        }
        m
    };
    let (lm, rm) = (millers(&left), millers(&right));
    let a = rand(1, lm.len(), 2);
    let b = rand(2, rm.len(), 5);
    let mut plan_l = Plan::new(dims, &lm, [0; 3]).unwrap();
    let mut fields = Vec::new();
    plan_l.to_real(&a, 0..2, &mut fields).unwrap();
    let mut plan_r = Plan::new(dims, &rm, g0).unwrap();
    let mut rf = Vec::new();
    plan_r.to_real(&b, 0..5, &mut rf).unwrap();
    let mut fft = PairFft::new(dims, &millers(&resp)).unwrap();
    let n = fields.len() / 2;
    let mut out = ComplexMatrix::zeros(resp.vectors().len(), 10).unwrap();
    for i in 0..2 {
        fft.pairs(&fields[i * n..(i + 1) * n], &rf, 5, &mut out, 5 * i)
            .unwrap();
    }
    // Direct: rho(q+G) = sum_G1 conj(a(G1)) b'(G1 + G), b' the k+q
    // representation: b'(M) = b(M + G0).
    let index: HashMap<[i32; 3], usize> = rm.iter().enumerate().map(|(i, m)| (*m, i)).collect();
    let mut worst = 0_f64;
    let mut scale = 0_f64;
    for i in 0..2 {
        for j in 0..5 {
            for (g, gm) in millers(&resp).iter().enumerate() {
                let mut sum = C::default();
                for (r1, m1) in lm.iter().enumerate() {
                    let target = [
                        m1[0] + gm[0] + g0[0],
                        m1[1] + gm[1] + g0[1],
                        m1[2] + gm[2] + g0[2],
                    ];
                    if let Some(&r2) = index.get(&target) {
                        sum += a[(r1, i)].conj() * b[(r2, j)];
                    }
                }
                worst = worst.max((out[(g, 5 * i + j)] - sum).norm());
                scale = scale.max(sum.norm());
            }
        }
    }
    println!("pair FFT vs convolution: max |diff| {worst:.3e}, max |rho| {scale:.3e}");
    assert!(worst <= 1e-12 * scale, "{worst:e}");
}

// ------------------------------------------------------- (a) Lindhard / FEG

fn lindhard(q: f64, w: f64, kf: f64) -> f64 {
    let z = q / (2. * kf);
    let u = w / (q * kf);
    let n = kf / (PI * PI);
    -n / 2.
        * (1.
            + (1. - z * z + u * u) / (4. * z)
                * (((z + 1.).powi(2) + u * u) / ((z - 1.).powi(2) + u * u)).ln()
            - u * (((1. + z) / u).atan() + ((1. - z) / u).atan()))
}

/// Exact discrete free-electron chi0_GG(q, i w) on the same Bloch vectors.
fn feg_direct(
    p: &Problem<'_>,
    mesh: &[[f64; 3]],
    resp: &PlaneWaveBasis,
    kf: f64,
    w: f64,
    g: usize,
) -> f64 {
    let b = p.cell.reciprocal_columns_inverse_bohr();
    let cart =
        |f: [f64; 3]| -> [f64; 3] { std::array::from_fn(|i| (0..3).map(|j| b[i][j] * f[j]).sum()) };
    let qg = resp.momenta()[g].cartesian;
    let mut sum = 0.;
    for &k in mesh {
        let basis = PlaneWaveBasis::at_k(&p.cell, Hartree(p.cutoff), k).unwrap();
        let kc = cart(basis.k_fractional());
        for g1 in basis.vectors() {
            let pk: [f64; 3] = std::array::from_fn(|i| kc[i] + g1.cartesian[i]);
            let e1: f64 = pk.iter().map(|x| x * x).sum::<f64>() / 2.;
            if e1 >= kf * kf / 2. {
                continue;
            }
            let e2: f64 = (0..3).map(|i| (pk[i] + qg[i]).powi(2)).sum::<f64>() / 2.;
            if e2 <= kf * kf / 2. || e2 > p.cutoff {
                continue;
            }
            let d = e1 - e2;
            sum += 2. * d / (d * d + w * w);
        }
    }
    sum * 2. / (mesh.len() as f64 * p.volume())
}

fn feg_case(n: usize) -> (f64, f64) {
    let a = 8.0;
    let kf = 0.9;
    let p = Problem {
        cell: Cell::from_user_rows_bohr([[a, 0., 0.], [0., a, 0.], [0., 0., a]]).unwrap(),
        // Every k+q+G1+G of an occupied G1 and |q+G| <= 1 stays inside.
        cutoff: 2.0,
        local: potential(&[]),
        nonlocal: None,
    };
    let fermi = kf * kf / 2.;
    let occupied = move |e: &[f64]| e.iter().filter(|&&x| x < fermi).count();
    let omegas = vec![0.1, 0.4];
    let run = Run {
        problem: &p,
        mesh: gamma_mesh([n; 3]),
        q: [2. / n as f64, 0., 0.],
        ecut_chi: 0.5,
        omegas: omegas.clone(),
        occupied: &occupied,
        bands: None,
        optical: false,
        block: 32,
    };
    let (chi0, resp) = chi0_run(&run, &dft_mpi::Comm::Single);
    let mut discrete = 0_f64;
    let mut continuum = 0_f64;
    let mut off = 0_f64;
    for (w, &omega) in omegas.iter().enumerate() {
        let body = &chi0.body[w];
        for g in 0..resp.vectors().len() {
            let exact = feg_direct(&p, &run.mesh, &resp, kf, omega, g);
            discrete = discrete.max(rel(body[(g, g)].re, exact));
            for h in 0..resp.vectors().len() {
                if h != g {
                    off = off.max(body[(g, h)].norm() / body[(g, g)].norm());
                }
            }
        }
        let q = resp.momenta()[0].norm_squared.sqrt();
        let l = lindhard(q, omega, kf);
        continuum = continuum.max(rel(body[(0, 0)].re, l));
        println!(
            "FEG n={n}: |q| {q:.4} Bohr^-1, w {omega} Ha: chi0_00 {:.10e}, Lindhard {l:.10e}, rel {:.3e}",
            body[(0, 0)].re,
            rel(body[(0, 0)].re, l)
        );
    }
    println!("FEG n={n}: discrete oracle max rel {discrete:.3e}, off-diagonal/diagonal {off:.3e}");
    assert!(discrete <= 1e-12, "discrete Lindhard sum {discrete:e}");
    assert!(off <= 1e-12, "plane-wave pairs couple only one G: {off:e}");
    (discrete, continuum)
}

#[test]
fn a_free_electron_gas_matches_the_lindhard_function() {
    // Pre-registered from the plain lattice-sum convergence of the
    // continuum integral (n = 6: 0.3 to 0.7 percent): 2 percent.
    let (_, continuum) = feg_case(6);
    assert!(continuum <= 2e-2, "continuum Lindhard {continuum:e}");
}

#[test]
#[ignore = "release: 8x8x8 free-electron mesh"]
fn a_free_electron_gas_converges_to_the_lindhard_function() {
    let (_, c6) = feg_case(6);
    let (_, c8) = feg_case(8);
    println!("continuum error n=6 {c6:.3e}, n=8 {c8:.3e}");
    assert!(c8 <= 1e-2);
}

// ---------------------------------------- insulating synthetic local crystal

fn skew_rows() -> [[f64; 3]; 3] {
    [[5.0, 0.0, 0.0], [0.6, 4.8, 0.0], [0.3, 0.5, 5.2]]
}
fn crystal_components() -> Vec<([i32; 3], C)> {
    vec![
        ([1, 0, 0], C::new(-1.3, 0.)),
        ([0, 1, 0], C::new(-1.1, 0.)),
        ([0, 0, 1], C::new(-1.2, 0.)),
        ([1, 1, 0], C::new(-0.3, 0.15)),
        ([0, 1, -1], C::new(-0.2, -0.1)),
    ]
}
fn crystal(cutoff: f64) -> Problem<'static> {
    Problem {
        cell: Cell::from_user_rows_bohr(skew_rows()).unwrap(),
        cutoff,
        local: potential(&crystal_components()),
        nonlocal: None,
    }
}
fn one_band(_: &[f64]) -> usize {
    1
}

/// Direct density coefficient n(Q) (Bohr^-3) of the lowest `occ` bands,
/// Q given by supercell Miller indices.
fn density_coefficient(solved: &[Solved], occ: usize, volume: f64, q: [i32; 3]) -> C {
    let mut sum = C::default();
    for s in solved {
        let index: HashMap<[i32; 3], usize> =
            s.miller.iter().enumerate().map(|(i, m)| (*m, i)).collect();
        for n in 0..occ {
            for (r1, m1) in s.miller.iter().enumerate() {
                if let Some(&r2) = index.get(&[m1[0] + q[0], m1[1] + q[1], m1[2] + q[2]]) {
                    sum += s.vectors[(r1, n)].conj() * s.vectors[(r2, n)];
                }
            }
        }
    }
    sum * 2. / (solved.len() as f64 * volume)
}

#[test]
fn b_static_columns_match_supercell_finite_differences() {
    let p = crystal(3.0);
    let mesh = gamma_mesh([3, 2, 2]);
    let q = [1. / 3., 0., 0.];
    // Gap check: one occupied band per cell.
    let mut vbm = f64::NEG_INFINITY;
    let mut cbm = f64::INFINITY;
    for &k in &mesh {
        let s = p.solve(k);
        vbm = vbm.max(s.energies[0]);
        cbm = cbm.min(s.energies[1]);
    }
    assert!(
        cbm - vbm > 0.05,
        "fixture must be an insulator: gap {}",
        cbm - vbm
    );
    let run = Run {
        problem: &p,
        mesh: mesh.clone(),
        q,
        ecut_chi: 2.0,
        omegas: vec![0.],
        occupied: &one_band,
        bands: None,
        optical: false,
        block: 16,
    };
    let (chi0, resp) = chi0_run(&run, &dft_mpi::Comm::Single);
    // Supercell 3 a1, a2, a3; primitive Miller h -> (3 h1, h2, h3).
    let rows = skew_rows();
    let sc = Problem {
        cell: Cell::from_user_rows_bohr([
            [3. * rows[0][0], 3. * rows[0][1], 3. * rows[0][2]],
            rows[1],
            rows[2],
        ])
        .unwrap(),
        cutoff: p.cutoff,
        local: p
            .local
            .iter()
            .map(|(h, v)| ([3 * h[0], h[1], h[2]], *v))
            .collect(),
        nonlocal: None,
    };
    let sc_mesh = gamma_mesh([1, 2, 2]);
    let to_sc = |g: [i32; 3]| [1 + 3 * g[0], g[1], g[2]];
    let resp_m = millers(&resp);
    let lambda = 1e-4;
    let mut worst = 0_f64;
    for g0 in [[0, 0, 0], [0, 1, 0], [-1, 0, 1]] {
        let column = resp_m.iter().position(|m| *m == g0).expect("G0 in sphere");
        let qq = to_sc(g0);
        let mut response = Vec::new();
        for sign in [1., -1.] {
            let extra = potential(&[(qq, C::new(sign * lambda / 2., 0.))]);
            let solved: Vec<Solved> = sc_mesh.iter().map(|&k| sc.solve_with(k, &extra)).collect();
            response.push(
                resp_m
                    .iter()
                    .map(|&g| density_coefficient(&solved, 3, sc.volume(), to_sc(g)))
                    .collect::<Vec<_>>(),
            );
        }
        let mut col_worst = 0_f64;
        let mut col_scale = 0_f64;
        for (g, _) in resp_m.iter().enumerate() {
            let fd = (response[0][g] - response[1][g]) / lambda;
            let sos = chi0.body[0][(g, column)];
            col_worst = col_worst.max((fd - sos).norm());
            col_scale = col_scale.max(sos.norm());
        }
        let r = col_worst / col_scale;
        println!(
            "(b) column G0={g0:?}: max |FD - SOS| {col_worst:.3e}, max |chi0| {col_scale:.3e} Bohr^-3/Ha, rel {r:.3e}"
        );
        worst = worst.max(r);
    }
    assert!(worst <= 1e-6, "finite-difference columns {worst:e}");
}

#[test]
fn c_d_hermitian_negative_and_high_frequency_limit() {
    let p = crystal(6.0);
    let mesh = gamma_mesh([2, 2, 1]);
    let big = 1e6;
    let run = Run {
        problem: &p,
        mesh: mesh.clone(),
        q: [0.5, 0., 0.],
        ecut_chi: 1.5,
        omegas: vec![0., 0.7, big],
        occupied: &one_band,
        bands: None,
        optical: false,
        block: 64,
    };
    let (chi0, resp) = chi0_run(&run, &dft_mpi::Comm::Single);
    for (w, body) in chi0.body.iter().enumerate() {
        let scale = max_abs(body);
        let mut herm = 0_f64;
        for j in 0..body.cols() {
            for i in 0..body.rows() {
                herm = herm.max((body[(i, j)] - body[(j, i)].conj()).norm());
            }
        }
        let (lambda, _) = hermitian_eigh(body).unwrap();
        let top = *lambda.last().unwrap();
        println!(
            "(c) w={}: max|chi0 - chi0^H| {herm:.3e}, largest eigenvalue {top:.3e}, scale {scale:.3e}",
            chi0.frequencies[w]
        );
        assert!(herm <= 1e-14 * scale);
        assert!(
            top <= 1e-13 * scale,
            "chi0 not negative semidefinite: {top:e}"
        );
    }
    // (d) chi0_GG'(i W) W^2 -> -2 s sum_k sum_n sum_{m emp} (e_m - e_n)
    // rho rho'^*: the closure over every empty state of the k + q basis,
    // evaluated without eigenvectors at k + q as (Y_G)^H (H' - e_n) Y_G' minus
    // the occupied-occupied part, Y_G = P' exp(i(q+G).r) psi_n by direct
    // convolution (no FFT). Exact up to roundoff and O((D_max / W)^2).
    let s = 2. / (mesh.len() as f64 * p.volume());
    let qw = resp.k_fractional();
    let m = millers(&resp);
    let ng = m.len();
    let mut closure = vec![C::default(); ng * ng];
    let mut occ_occ = vec![C::default(); ng * ng];
    for &k in &mesh {
        let left = p.solve(k);
        let kq = [k[0] + qw[0], k[1] + qw[1], k[2] + qw[2]];
        let right = p.solve(kq);
        let (_, h) = p.hamiltonian(kq, &HashMap::new());
        let g0 = shift(left.basis.k_fractional(), qw, right.basis.k_fractional());
        let index: HashMap<[i32; 3], usize> = right
            .miller
            .iter()
            .enumerate()
            .map(|(i, x)| (*x, i))
            .collect();
        let nr = right.miller.len();
        let y: Vec<Vec<C>> = m
            .iter()
            .map(|g| {
                let mut v = vec![C::default(); nr];
                for (r1, m1) in left.miller.iter().enumerate() {
                    let t = [
                        m1[0] + g[0] + g0[0],
                        m1[1] + g[1] + g0[1],
                        m1[2] + g[2] + g0[2],
                    ];
                    if let Some(&r2) = index.get(&t) {
                        v[r2] = left.vectors[(r1, 0)];
                    }
                }
                v
            })
            .collect();
        let en = left.energies[0];
        let hy: Vec<Vec<C>> = y
            .iter()
            .map(|v| {
                (0..nr)
                    .map(|i| (0..nr).map(|j| h[(i, j)] * v[j]).sum::<C>() - en * v[i])
                    .collect()
            })
            .collect();
        let proj: Vec<C> = y
            .iter()
            .map(|v| (0..nr).map(|i| right.vectors[(i, 0)].conj() * v[i]).sum())
            .collect();
        for a in 0..ng {
            for b in 0..ng {
                let full: C = (0..nr).map(|i| y[a][i].conj() * hy[b][i]).sum();
                closure[a + b * ng] += -2. * s * full;
                occ_occ[a + b * ng] +=
                    -2. * s * (right.energies[0] - en) * proj[a].conj() * proj[b];
            }
        }
    }
    let mut worst = 0_f64;
    let mut scale = 0_f64;
    let mut tr = 0_f64;
    for a in 0..ng {
        for b in 0..ng {
            let expected = closure[a + b * ng] - occ_occ[a + b * ng];
            let got = chi0.body[2][(a, b)] * big * big;
            worst = worst.max((got - expected).norm());
            scale = scale.max(expected.norm());
            tr = tr.max(occ_occ[a + b * ng].norm());
        }
    }
    println!(
        "(d) high-frequency closure: max |chi0 W^2 - closure| {worst:.3e}, scale {scale:.3e}, rel {:.3e}; occupied-occupied part (cancels by time reversal) {tr:.3e}",
        worst / scale
    );
    assert!(worst <= 1e-9 * scale, "closure {:e}", worst / scale);
    assert!(tr <= 1e-12 * scale, "time-reversal cancellation {tr:e}");
    // Physical f-sum -(q+G).(q+G') n(G-G'): holds for a complete basis only
    // (exp(i(q+G).r) psi_n leaves the k + q sphere); recorded, and its
    // truncation error must fall with the cutoff.
    let fsum = |p: &Problem<'_>| {
        let run = Run {
            problem: p,
            mesh: mesh.clone(),
            q: [0.5, 0., 0.],
            ecut_chi: 1.5,
            omegas: vec![big],
            occupied: &one_band,
            bands: None,
            optical: false,
            block: 64,
        };
        let (chi0, resp) = chi0_run(&run, &dft_mpi::Comm::Single);
        let solved: Vec<Solved> = mesh.iter().map(|&k| p.solve(k)).collect();
        let m = millers(&resp);
        let mut worst = 0_f64;
        let mut scale = 0_f64;
        for (i, gi) in m.iter().enumerate() {
            for (j, gj) in m.iter().enumerate() {
                let h = [gi[0] - gj[0], gi[1] - gj[1], gi[2] - gj[2]];
                let n = density_coefficient(&solved, 1, p.volume(), h);
                let a = resp.momenta()[i].cartesian;
                let b = resp.momenta()[j].cartesian;
                let dot: f64 = (0..3).map(|x| a[x] * b[x]).sum();
                worst = worst.max((chi0.body[0][(i, j)] * big * big + dot * n).norm());
                scale = scale.max((dot * n).norm());
            }
        }
        worst / scale
    };
    let (f6, f12) = (fsum(&p), fsum(&crystal(12.0)));
    println!("(d) physical f-sum relative deviation: E_cut 6 Ha {f6:.3e}, 12 Ha {f12:.3e}");
    assert!(
        f12 < 0.5 * f6,
        "f-sum truncation error must fall with the cutoff"
    );
}

// ------------------------------------------- (e) k.p versus finite q, nonlocal

fn nonlocal_crystal(pp: &dft_upf::Pseudopotential) -> Problem<'_> {
    let cell = Cell::from_user_rows_bohr(skew_rows()).unwrap();
    let tau = cell.fractional_to_cartesian_bohr([0.13, 0.27, 0.61]);
    Problem {
        cell,
        cutoff: 3.0,
        local: potential(&crystal_components()),
        nonlocal: Some((pp, vec![tau])),
    }
}

#[test]
fn e_kp_head_and_wings_match_finite_q_with_the_nonlocal_velocity() {
    let pp = synthetic_upf();
    let p = nonlocal_crystal(&pp);
    let mesh = gamma_mesh([2, 2, 2]);
    let optical = Run {
        problem: &p,
        mesh: mesh.clone(),
        q: [0.; 3],
        ecut_chi: 2.0,
        omegas: vec![0., 0.5],
        occupied: &one_band,
        bands: None,
        optical: true,
        block: 16,
    };
    let (kp, resp) = chi0_run(&optical, &dft_mpi::Comm::Single);
    let a = p.cell.lattice_columns_bohr();
    let mut worst_head = 0_f64;
    let mut worst_wing = 0_f64;
    let mut worst_local = 0_f64;
    for dir in [[1.0_f64, 0., 0.], [0.3, -0.5, 0.81]] {
        let norm: f64 = (dir[0] * dir[0] + dir[1] * dir[1] + dir[2] * dir[2]).sqrt();
        let hat: [f64; 3] = std::array::from_fn(|i| dir[i] / norm);
        let finite = |delta: f64| {
            // fractional q = A^T q / (2 pi)
            let qf: [f64; 3] = std::array::from_fn(|j| {
                (0..3).map(|i| a[i][j] * hat[i] * delta).sum::<f64>() / (2. * PI)
            });
            let run = Run {
                problem: &p,
                mesh: mesh.clone(),
                q: qf,
                ecut_chi: 2.0,
                omegas: vec![0., 0.5],
                occupied: &one_band,
                bands: None,
                optical: false,
                block: 16,
            };
            chi0_run(&run, &dft_mpi::Comm::Single).0
        };
        let d = 2e-3;
        let (c1, c2) = (finite(d), finite(2. * d));
        for w in 0..2 {
            let head: f64 = (0..3)
                .flat_map(|x| (0..3).map(move |y| (x, y)))
                .map(|(x, y)| (kp.head[w][x][y] * hat[x] * hat[y]).re)
                .sum();
            let f1 = c1.body[w][(0, 0)].re / (d * d);
            let f2 = c2.body[w][(0, 0)].re / (4. * d * d);
            let extrapolated = (4. * f1 - f2) / 3.;
            let r = rel(extrapolated, head);
            println!(
                "(e) q^={hat:.3?} w={}: k.p q^.L.q^ {head:.10e}, finite-q extrapolated {extrapolated:.10e}, rel {r:.3e}",
                kp.frequencies[w]
            );
            worst_head = worst_head.max(r);
            for g in 1..resp.vectors().len() {
                let wing: C = (0..3).map(|b| kp.wings[w][(g, b)] * hat[b]).sum();
                let e1 = c1.body[w][(g, 0)] / d;
                let e2 = c2.body[w][(g, 0)] / (2. * d);
                let ex = e1 * 2. - e2;
                worst_wing =
                    worst_wing.max((ex - wing).norm() / wing.norm().max(1e-3 * head.abs()));
            }
        }
        // Local-only velocity must miss the nonlocal commutator: the
        // nonlocal term is measurable on this fixture.
        let local_only = {
            let mut acc = Accumulator::new(resp.vectors().len(), &[0.], true, 1).unwrap();
            let dims = p.dims();
            let mut fft = PairFft::new(dims, &millers(&resp)).unwrap();
            for &k in &mesh {
                let s = p.solve(k);
                let o = KOrbitals {
                    miller: &s.miller,
                    coefficients: &s.vectors,
                    energies: &s.energies,
                };
                let mut t = Transitions::new(dims, o, 1, &s.miller, [0; 3]).unwrap();
                let v = Velocity {
                    momenta: &s.momenta,
                    projectors: None,
                };
                t.add(
                    &mut acc,
                    &mut fft,
                    o,
                    1..s.energies.len(),
                    2. / (mesh.len() as f64 * p.volume()),
                    Some(v),
                )
                .unwrap();
            }
            let c = acc.reduce(&dft_mpi::Comm::Single).unwrap();
            (0..3)
                .flat_map(|x| (0..3).map(move |y| (x, y)))
                .map(|(x, y)| (c.head[0][x][y] * hat[x] * hat[y]).re)
                .sum::<f64>()
        };
        let full: f64 = (0..3)
            .flat_map(|x| (0..3).map(move |y| (x, y)))
            .map(|(x, y)| (kp.head[0][x][y] * hat[x] * hat[y]).re)
            .sum();
        worst_local = worst_local.max(rel(local_only, full));
    }
    println!(
        "(e) max rel head {worst_head:.3e}, wings {worst_wing:.3e}; local-only velocity differs by {worst_local:.3e}"
    );
    assert!(worst_head <= 1e-6, "k.p head vs finite q {worst_head:e}");
    assert!(worst_wing <= 1e-5, "k.p wings vs finite q {worst_wing:e}");
    assert!(
        worst_local >= 1e-3,
        "the fixture must exercise the nonlocal velocity"
    );
    // Macroscopic tensor: symmetric, > 1, local fields reduce it.
    let eps = dielectric::optical(
        &kp,
        0,
        &resp
            .vectors()
            .iter()
            .map(|g| g.norm_squared.sqrt())
            .collect::<Vec<_>>(),
    )
    .unwrap();
    println!(
        "eps_M without LF {:.6?}\neps_M with LF {:.6?}",
        eps.without_local_fields, eps.with_local_fields
    );
    assert!(eps.max_asymmetry < 1e-10);
    for a in 0..3 {
        assert!(eps.without_local_fields[a][a] > 1.);
        assert!(eps.with_local_fields[a][a] > 1.);
        assert!(eps.with_local_fields[a][a] < eps.without_local_fields[a][a]);
    }
}

#[test]
fn projector_gradients_match_the_operator_and_finite_differences() {
    use dft_physics::hamiltonian::{NonlocalOperator, Species};
    use dft_physics::stress::projector_k_gradient;
    let pp = synthetic_upf();
    let cell = Cell::from_user_rows_bohr(skew_rows()).unwrap();
    let tau = [
        cell.fractional_to_cartesian_bohr([0.13, 0.27, 0.61]),
        cell.fractional_to_cartesian_bohr([0.5, 0.5, 0.5]),
    ];
    let species = [Species {
        potential: &pp,
        positions_bohr: &tau,
    }];
    let a = cell.lattice_columns_bohr();
    let mut worst_value = 0_f64;
    let mut worst_fd = 0_f64;
    for k in [[0., 0., 0.], [0.2, -0.35, 0.1]] {
        let basis = PlaneWaveBasis::at_k(&cell, Hartree(4.0), k).unwrap();
        let grad = projector_k_gradient(&basis, cell.volume_bohr3(), &species).unwrap();
        let op = NonlocalOperator::new(&basis, cell.volume_bohr3(), &species).unwrap();
        let n = basis.vectors().len();
        let rows = op.projector_rows(&(0..n).collect::<Vec<_>>()).unwrap();
        for c in 0..grad.values.cols() {
            for g in 0..n {
                worst_value = worst_value.max((rows[c * n + g] - grad.values[(g, c)]).norm());
            }
        }
        // d/dk_b of the operator's rows = gradient - i tau_b beta.
        let h = 1e-5;
        let index: HashMap<[i32; 3], usize> = basis
            .vectors()
            .iter()
            .enumerate()
            .map(|(i, g)| (g.miller, i))
            .collect();
        let atom_of = |c: usize| if c < grad.channels[0] { 0 } else { 1 };
        for b in 0..3 {
            let shifted = |s: f64| {
                let kf: [f64; 3] = std::array::from_fn(|j| k[j] + s * a[b][j] / (2. * PI));
                let bs = PlaneWaveBasis::at_k(&cell, Hartree(4.0), kf).unwrap();
                let op = NonlocalOperator::new(&bs, cell.volume_bohr3(), &species).unwrap();
                let m = bs.vectors().len();
                let r = op.projector_rows(&(0..m).collect::<Vec<_>>()).unwrap();
                (
                    bs.vectors().iter().map(|g| g.miller).collect::<Vec<_>>(),
                    r,
                    m,
                )
            };
            let (mp, rp, np_) = shifted(h);
            let (mm, rm, nm_) = shifted(-h);
            let ip: HashMap<[i32; 3], usize> =
                mp.iter().enumerate().map(|(i, m)| (*m, i)).collect();
            let im: HashMap<[i32; 3], usize> =
                mm.iter().enumerate().map(|(i, m)| (*m, i)).collect();
            for (miller, &g) in &index {
                // The oracle's own harmonics lose digits within h of the
                // polar axis (cos theta cancellation, as in the stress
                // harmonic test); the gradient there is covered by that test.
                let q = basis.momenta()[g].cartesian;
                if q[0].hypot(q[1]) < 1e-2 * basis.momenta()[g].norm_squared.sqrt() {
                    continue;
                }
                let (Some(&gp), Some(&gm)) = (ip.get(miller), im.get(miller)) else {
                    continue;
                };
                for c in 0..grad.values.cols() {
                    let fd = (rp[c * np_ + gp] - rm[c * nm_ + gm]) / (2. * h);
                    let t = tau[atom_of(c)][b];
                    let expected = grad.gradients[b][(g, c)] + C::new(0., -t) * grad.values[(g, c)];
                    worst_fd = worst_fd.max((fd - expected).norm());
                }
            }
        }
    }
    println!(
        "projector values vs NonlocalOperator {worst_value:.3e}; k gradient vs FD {worst_fd:.3e}"
    );
    assert!(worst_value <= 1e-14);
    assert!(worst_fd <= 1e-7, "{worst_fd:e}");
}

#[test]
fn optical_head_converges_monotonically_in_the_band_count() {
    let pp = synthetic_upf();
    let p = nonlocal_crystal(&pp);
    let mesh = gamma_mesh([2, 1, 1]);
    let mut previous = f64::INFINITY;
    let mut line = Vec::new();
    for bands in [Some(4), Some(8), Some(16), Some(24), None] {
        let run = Run {
            problem: &p,
            mesh: mesh.clone(),
            q: [0.; 3],
            ecut_chi: 1.0,
            omegas: vec![0.],
            occupied: &one_band,
            bands,
            optical: true,
            block: 8,
        };
        let (c, _) = chi0_run(&run, &dft_mpi::Comm::Single);
        let eps = 1. - 4. * PI * c.head[0][0][0].re;
        line.push(format!("{bands:?}: {eps:.8}"));
        // Every added transition adds a negative term to L_xx.
        assert!(eps >= previous - 1e-12 || previous.is_infinite());
        previous = eps;
    }
    println!("eps_xx (no local fields) versus bands: {}", line.join(", "));
}

#[test]
fn rank_reduction_is_deterministic_and_one_rank_is_bit_identical() {
    let p = crystal(3.0);
    let run = Run {
        problem: &p,
        mesh: gamma_mesh([3, 1, 1]),
        q: [1. / 3., 0., 0.],
        ecut_chi: 1.5,
        omegas: vec![0., 0.3],
        occupied: &one_band,
        bands: None,
        optical: false,
        block: 5,
    };
    let (serial, _) = chi0_run(&run, &dft_mpi::Comm::Single);
    let bits = |c: &Chi0| -> Vec<u64> {
        c.body
            .iter()
            .flat_map(|m| {
                m.as_slice()
                    .iter()
                    .flat_map(|z| [z.re.to_bits(), z.im.to_bits()])
            })
            .collect()
    };
    for ranks in [1, 2, 3] {
        let results = dft_mpi::simulate(ranks, |comm| chi0_run(&run, comm).0);
        let first = bits(&results[0]);
        for r in &results {
            assert_eq!(bits(r), first, "ranks disagree at np={ranks}");
            assert_eq!(r.transitions, serial.transitions);
        }
        if ranks == 1 {
            assert_eq!(first, bits(&serial), "np=1 not bit-identical");
        }
        let mut worst = 0_f64;
        for (a, b) in results[0].body.iter().zip(&serial.body) {
            for (x, y) in a.as_slice().iter().zip(b.as_slice()) {
                worst = worst.max((x - y).norm());
            }
        }
        println!("np={ranks}: max |chi0 - serial| {worst:.3e}");
        assert!(worst <= 1e-14 * max_abs(&serial.body[0]));
    }
}
