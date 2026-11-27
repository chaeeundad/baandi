//! `baandi dielectric` (R2, `docs/roadmap/rpa-crpa-plan.md` 2.11, 3.6):
//! reference SCF, then the optical-limit dielectric response of
//! `dft_engine::response` on the `[response]` Gamma-centred full k mesh.
//! Outputs: `epsilon.json` (complete record, units and conventions),
//! `epsilon_M.csv` (one row per imaginary frequency, both tensors),
//! `epsilon_M.dat` (gnuplot columns), `chi0_head.dat` (head tensor L versus
//! i w), `results.json` (copy of `epsilon.json`).
use crate::bands::{Context, err, ev};
use dft_cli::run::{RunError, sha256};
use dft_engine::{
    parallel::SharedComm,
    refscf::nscf::BandSolver,
    refstart::StartInput,
    response::{DielectricOptions, DipoleDump, dielectric},
};
use dft_geometry::kpoints::KPoint;
use serde_json::json;
use std::{fmt::Write as _, path::Path};

const OUTPUTS: [&str; 6] = [
    "epsilon.json",
    "dipoles.csv",
    "epsilon_M.csv",
    "epsilon_M.dat",
    "chi0_head.dat",
    "results.json",
];

fn gamma_mesh(n: [usize; 3]) -> Result<Vec<KPoint>, RunError> {
    if n.contains(&0) {
        return Err(err(1, "response.k_mesh entries must be positive"));
    }
    let total = (n[0] * n[1] * n[2]) as f64;
    let mut points = Vec::new();
    for i in 0..n[0] {
        for j in 0..n[1] {
            for k in 0..n[2] {
                points.push(KPoint {
                    fractional: [
                        i as f64 / n[0] as f64,
                        j as f64 / n[1] as f64,
                        k as f64 / n[2] as f64,
                    ],
                    weight: 1. / total,
                });
            }
        }
    }
    Ok(points)
}

/// Six independent components (xx, yy, zz, xy, xz, yz) and the trace / 3.
fn components(t: &[[f64; 3]; 3]) -> [f64; 7] {
    [
        t[0][0],
        t[1][1],
        t[2][2],
        t[0][1],
        t[0][2],
        t[1][2],
        (t[0][0] + t[1][1] + t[2][2]) / 3.,
    ]
}

/// R2-VEL-CHECK dipole record (format: tests/reference/response/r2/README.md).
fn dipole_csv(records: &[dft_engine::response::KDipoles]) -> String {
    let mut out = String::from(
        "# baandi dielectric dipoles (R2-VEL-CHECK diagnostic): gauge-invariant q -> 0 dipole sums per pair of degenerate groups (V occupied, C empty; a group = consecutive bands with eigenvalue spacing < 1e-5 Ha); bands 1-based inclusive in ascending energy per k; k1..k3 fractional (reciprocal-lattice basis), kx..kz Cartesian Bohr^-1 (input lattice frame); e_v_ha, e_c_ha group mean eigenvalues (Ha); orbitals cell-normalized, spin-unpolarized, one spin channel (no factor 2); p_a = <u_v| dH_k/dk_a |u_c> (Hartree atomic units, Ha Bohr = a.u. of velocity) = <v| -i[r_a, H] |c>, including the nonlocal commutator; s_a = sum_{v,c} |p_a|^2 (a.u.^2); r_a = sum_{v,c} |p_a|^2 / (e_c - e_v)^2 = sum |<v| r_a |c>|^2 (Bohr^2, per-band energies); *_local the same with the kinetic velocity k + G only (no [V_NL, r]); c_at_band_limit = 1 when C ends at the last computed band (possibly incomplete group)\nk_index,k1,k2,k3,kx,ky,kz,v_first,v_last,c_first,c_last,e_v_ha,e_c_ha,c_at_band_limit,s_x,s_y,s_z,r_x,r_y,r_z,s_local_x,s_local_y,s_local_z,r_local_x,r_local_y,r_local_z\n",
    );
    for k in records {
        for p in &k.pairs {
            write!(
                out,
                "{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{},{},{},{},{:.12e},{:.12e},{}",
                k.k_index,
                k.k_fractional[0],
                k.k_fractional[1],
                k.k_fractional[2],
                k.k_cartesian[0],
                k.k_cartesian[1],
                k.k_cartesian[2],
                p.v.0 + 1,
                p.v.1 + 1,
                p.c.0 + 1,
                p.c.1 + 1,
                p.e_v,
                p.e_c,
                u8::from(p.c_at_band_limit)
            )
            .unwrap();
            for v in p
                .velocity
                .iter()
                .chain(&p.position)
                .chain(&p.velocity_local)
                .chain(&p.position_local)
            {
                write!(out, ",{v:.12e}").unwrap();
            }
            out.push('\n');
        }
    }
    out
}

pub fn execute(input: &Path, output: Option<&Path>) -> Result<(), RunError> {
    let ctx = Context::load(
        input,
        output,
        dft_cli::driver_select::Command::Bands,
        &OUTPUTS,
    )?;
    let section = ctx
        .config
        .response
        .as_ref()
        .ok_or_else(|| err(1, "dielectric requires a [response] table"))?;
    if !ctx.selection.reference() {
        return Err(err(1, "dielectric requires the reference driver"));
    }
    if ctx.config.collinear() {
        return Err(err(
            1,
            "dielectric: spin-polarized response is not implemented (R3)",
        ));
    }
    let head = section.head.as_deref().unwrap_or("kp");
    if head != "kp" {
        return Err(err(
            1,
            format!(
                "response.head = {head:?}: only \"kp\" is available (finite-q is an internal check)"
            ),
        ));
    }
    let name = section.band_solver.as_deref().unwrap_or("dense");
    let solver = match BandSolver::parse(name) {
        Some(s @ (BandSolver::Dense | BandSolver::Chebyshev)) => s,
        _ => {
            return Err(err(
                1,
                format!("response.band_solver = {name:?}: use dense or chebyshev"),
            ));
        }
    };
    let bands = match &section.bands {
        None => None,
        Some(toml::Value::String(s)) if s == "auto" => None,
        Some(toml::Value::Integer(n)) if *n > 0 => Some(*n as usize),
        Some(v) => {
            return Err(err(
                1,
                format!("response.bands = {v}: use \"auto\" or a positive integer"),
            ));
        }
    };
    let mut options = DielectricOptions::new(section.ecut_chi_ha);
    options.bands = bands;
    options.solver = solver;
    if let Some(t) = section.eigen_tolerance_ha {
        options.eigen_tolerance = t;
    }
    if let Some(f) = &section.frequency {
        let grid = f.grid.as_deref().unwrap_or("gauss-legendre");
        if grid != "gauss-legendre" {
            return Err(err(
                1,
                format!("response.frequency.grid = {grid:?}: only \"gauss-legendre\" is available"),
            ));
        }
        if let Some(p) = f.points {
            options.frequency_points = p;
        }
        options.omega0_hartree = match &f.omega0_ha {
            None => None,
            Some(toml::Value::String(s)) if s == "auto" => None,
            Some(toml::Value::Float(w)) if *w > 0. => Some(*w),
            Some(toml::Value::Integer(w)) if *w > 0 => Some(*w as f64),
            Some(v) => {
                return Err(err(
                    1,
                    format!(
                        "response.frequency.omega0_ha = {v}: use \"auto\" or a positive number"
                    ),
                ));
            }
        };
    }
    // R2-VEL-CHECK diagnostic (switch RESPONSE_DIPOLE_DUMP, default off).
    let dipole_dump =
        dft_numerics::switches::choice(dft_numerics::switches::Id::RESPONSE_DIPOLE_DUMP);
    options.dipole_dump = match dipole_dump {
        "gamma" => DipoleDump::GammaAndFirst,
        "all" => DipoleDump::All,
        _ => DipoleDump::Off,
    };
    let headline_lf = section.local_fields.unwrap_or(true);
    let points = gamma_mesh(section.k_mesh)?;
    let scf = ctx.scf()?;
    let density_bytes: Vec<u8> = scf
        .density
        .iter()
        .flat_map(|s| s.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<_>>())
        .collect();
    let ranks = SharedComm::active().map_or(1, |c| dft_mpi::Collective::size(&c));
    println!(
        "dielectric: {} k points, ecut_chi {} Ha, bands {}, solver {}, {} imaginary frequencies + static",
        points.len(),
        section.ecut_chi_ha,
        bands.map_or("auto".to_owned(), |b| b.to_string()),
        solver.name(),
        options.frequency_points
    );
    let potentials: Vec<_> = ctx
        .potentials
        .iter()
        .map(|(label, potential)| dft_engine::SpeciesPotential { label, potential })
        .collect();
    let result = dielectric(
        StartInput {
            cell: &ctx.cell,
            atoms: &ctx.atoms,
            potentials: &potentials,
            options: &ctx.options,
        },
        &scf.density,
        &points,
        options,
    )
    .map_err(|e| err(4, e))?;
    let mut csv = String::from(
        "# baandi dielectric (R2): macroscopic dielectric tensor eps_M(i omega) in the optical limit q -> 0 (RPA, no xc kernel), dimensionless; omega on the imaginary axis (omega_ha Hartree, omega_ev eV); weight_ha = Gauss-Legendre weight of int_0^inf d omega (0 for the static point); nlf_* without local fields (head only), lf_* with local fields; Cartesian axes of the input lattice frame\nindex,omega_ha,omega_ev,weight_ha,nlf_xx,nlf_yy,nlf_zz,nlf_xy,nlf_xz,nlf_yz,nlf_avg,lf_xx,lf_yy,lf_zz,lf_xy,lf_xz,lf_yz,lf_avg\n",
    );
    let mut dat = String::from(
        "# baandi dielectric (R2): eps_M(i omega), q -> 0, RPA; columns: omega_ev, avg without local fields, avg with local fields, lf_xx, lf_yy, lf_zz\n",
    );
    let mut head_dat = String::from(
        "# baandi dielectric (R2): head of chi0 in the optical limit, chi0_00(q, i omega) = q^2 q^.L.q^; columns: omega_ha, Re L_xx, Re L_yy, Re L_zz, Re L_xy, Re L_xz, Re L_yz (Bohr^-1 Ha^-1)\n",
    );
    for (i, (w, e)) in result.frequencies.iter().zip(&result.epsilon).enumerate() {
        let n = components(&e.without_local_fields);
        let l = components(&e.with_local_fields);
        write!(
            csv,
            "{i},{w:.16e},{:.16e},{:.16e}",
            ev(*w),
            result.weights[i]
        )
        .unwrap();
        for v in n.iter().chain(&l) {
            write!(csv, ",{v:.16e}").unwrap();
        }
        csv.push('\n');
        writeln!(
            dat,
            "{:.10e} {:.12e} {:.12e} {:.12e} {:.12e} {:.12e}",
            ev(*w),
            n[6],
            l[6],
            l[0],
            l[1],
            l[2]
        )
        .unwrap();
        let h = &result.head[i];
        writeln!(
            head_dat,
            "{w:.10e} {:.12e} {:.12e} {:.12e} {:.12e} {:.12e} {:.12e}",
            h[0][0].re, h[1][1].re, h[2][2].re, h[0][1].re, h[0][2].re, h[1][2].re
        )
        .unwrap();
    }
    let e0 = &result.epsilon[0];
    let headline = if headline_lf {
        e0.with_local_fields
    } else {
        e0.without_local_fields
    };
    let mut meta = json!({
        "schema_version": 1,
        "kind": "baandi_dielectric_r2",
        "code": "rust-dft",
        "code_version": env!("CARGO_PKG_VERSION"),
        "title": ctx.config.system.name,
        "input_sha256": ctx.input_hash,
        "pseudopotentials": ctx.pp_metadata,
        "scf_density_sha256": sha256(&density_bytes),
        "scf_energy_hartree": scf.energy,
        "scf_iterations": scf.iterations,
        "scf_path": scf.path,
        "functional": ctx.options.functional.name(),
        "cutoff_hartree": ctx.options.cutoff_hartree,
        "threads": ctx.options.threads,
        "ranks": ranks,
        "theory": "RPA (Hartree kernel only, no xc kernel), independent-particle chi0 by the sum over states (Adler-Wiser, decision D1); optical limit q -> 0 with k.p head and wings, velocity d H_k / d k including the nonlocal commutator",
        "units": {"epsilon": "dimensionless", "frequency": "Hartree on the imaginary axis", "head_tensor": "Bohr^-1 Ha^-1 (chi0_00 = q^2 q^.L.q^)", "energies": "Hartree"},
        "conventions": "eps~ = 1 - v^1/2 chi0 v^1/2, v = 4 pi / |q+G|^2; eps_M(q^) = q^.E.q^ = 1 / [eps~^-1]_00; Cartesian axes of the input lattice frame; tensors are the real symmetric parts (max_asymmetry gives the dropped part)",
        "settings": {
            "ecut_chi_hartree": section.ecut_chi_ha,
            "response_size": result.response_size,
            "bands": result.bands,
            "bands_rule": if bands.is_some() { "explicit" } else { "auto = N_G(ecut_chi) clamped to the smallest k basis" },
            "band_solver": solver.name(),
            "eigen_tolerance_hartree": options.eigen_tolerance,
            "k_mesh": section.k_mesh,
            "k_sampling": "Gamma-centred full mesh, no symmetry reduction",
            "head": "kp",
            "headline_local_fields": headline_lf,
            "frequency_grid": "gauss-legendre mapped w = w0 t / (1 - t), plus w = 0",
            "frequency_points": options.frequency_points,
            "omega0_hartree": result.omega0,
            "omega0_rule": if options.omega0_hartree.is_some() { "explicit" } else { "auto = mean(lowest direct transition, largest transition)" },
            "band_block": options.band_block,
            "dipole_dump": dipole_dump,
        },
        "spectrum": {
            "electrons": result.electrons,
            "occupied_bands": result.occupied,
            "valence_maximum_hartree": result.valence_maximum,
            "conduction_minimum_hartree": result.conduction_minimum,
            "indirect_gap_ev": ev(result.conduction_minimum - result.valence_maximum),
            "lowest_direct_gap_ev": ev(result.lowest_direct_gap),
            "largest_transition_hartree": result.largest_transition,
            "max_residual_hartree": result.max_residual,
            "basis_size_min": result.basis_sizes.0,
            "basis_size_max": result.basis_sizes.1,
        },
        "transitions": result.transitions,
        "transition_units": result.units,
        "epsilon_infinity": {
            "headline": headline,
            "headline_average": components(&headline)[6],
            "with_local_fields": e0.with_local_fields,
            "without_local_fields": e0.without_local_fields,
            "average_with_local_fields": components(&e0.with_local_fields)[6],
            "average_without_local_fields": components(&e0.without_local_fields)[6],
            "min_body_eigenvalue": e0.min_body_eigenvalue,
            "max_asymmetry": e0.max_asymmetry,
        },
        "frequencies_hartree": result.frequencies,
        "weights_hartree": result.weights,
        "epsilon_with_local_fields": result.epsilon.iter().map(|e| e.with_local_fields).collect::<Vec<_>>(),
        "epsilon_without_local_fields": result.epsilon.iter().map(|e| e.without_local_fields).collect::<Vec<_>>(),
        "head_tensor_real": result.head.iter().map(|h| h.map(|r| r.map(|z| z.re))).collect::<Vec<_>>(),
        "head_tensor_imag": result.head.iter().map(|h| h.map(|r| r.map(|z| z.im))).collect::<Vec<_>>(),
        "limitations": "insulators only (metals reject, D6); spin-unpolarized; q -> 0 only; no symmetry reduction; no checkpoint; RPA overestimates experiment (LDA gap), never compare to experiment as a gate",
    });
    ctx.selection.record(&mut meta);
    if !result.dipoles.is_empty() {
        ctx.write("dipoles.csv", &dipole_csv(&result.dipoles))?;
    }
    ctx.write("epsilon_M.csv", &csv)?;
    ctx.write("epsilon_M.dat", &dat)?;
    ctx.write("chi0_head.dat", &head_dat)?;
    ctx.json("epsilon.json", &meta)?;
    ctx.json("results.json", &meta)?;
    println!(
        "dielectric done: epsilon_infinity (average, {} local fields) {:.6}; outputs in {}",
        if headline_lf { "with" } else { "without" },
        components(&headline)[6],
        ctx.out.display()
    );
    Ok(())
}
