//! `baandi subspace` (R4a, `docs/roadmap/rpa-crpa-plan.md` 2.8, 2.11, 3.6;
//! ADR 0059): reference SCF, then atomic-orbital subspace projection with
//! Loewdin orthonormalization and Wigner-Seitz real-space hoppings `H(R)`.
//!
//! Outputs:
//! - `subspace.json` (complete record, units, conventions, `projector_fingerprint`,
//!   per-orbital populations, per-k overlap/orthogonality/spilling metrics, `H(R)` summary)
//! - `subspace_hr.dat` (Wannier90 `seedname_hr.dat` layout, hoppings in eV via `dft_units`)
//! - `subspace_populations.csv` (per-orbital raw and Loewdin populations)
//! - `subspace_weights.csv` (per-k, per-band energies and subspace weights `p_nk`)
//! - `subspace_bands.csv` (optional, when `[bands]` path is present in input)
//! - `results.json` (copy of `subspace.json`)
use crate::bands::{Context, err, ev};
use dft_cli::run::{RunError, sha256};
use dft_engine::{
    refscf::nscf::BandSolver,
    refstart::StartInput,
    response::subspace::{
        BandWindow, ProjectionSelector, SubspaceOptions, SubspaceResult, subspace,
    },
};
use dft_geometry::kpath;
use serde_json::json;
use std::{fmt::Write as _, path::Path};

const OUTPUTS: [&str; 6] = [
    "subspace.json",
    "subspace_hr.dat",
    "subspace_populations.csv",
    "subspace_weights.csv",
    "subspace_bands.csv",
    "results.json",
];

pub(crate) fn parse_solver(name: Option<&str>) -> Result<BandSolver, RunError> {
    let label = name.unwrap_or("dense");
    BandSolver::parse(label).ok_or_else(|| {
        err(
            1,
            format!("subspace.band_solver '{label}' must be 'dense', 'davidson' or 'chebyshev'"),
        )
    })
}

fn parse_window(sub: &dft_cli::config::Subspace) -> Result<BandWindow, RunError> {
    match (sub.window_ha, sub.band_range) {
        (Some(_), Some(_)) => Err(err(
            1,
            "subspace: set either window_ha or band_range, not both",
        )),
        (Some([min_ha, max_ha]), None) => {
            if !min_ha.is_finite() || !max_ha.is_finite() || min_ha >= max_ha {
                return Err(err(1, "subspace.window_ha must satisfy min < max"));
            }
            Ok(BandWindow::EnergyRelative { min_ha, max_ha })
        }
        (None, Some([first, last])) => {
            if first == 0 || first > last {
                return Err(err(
                    1,
                    "subspace.band_range uses 1-based inclusive indices with first <= last",
                ));
            }
            Ok(BandWindow::BandRange {
                first: first - 1,
                last: last - 1,
            })
        }
        (None, None) => Ok(BandWindow::All),
    }
}

pub(crate) fn parse_projections(
    sub: &dft_cli::config::Subspace,
) -> Result<Vec<ProjectionSelector>, RunError> {
    let has_shorthand = sub.species.is_some()
        || sub.l.is_some()
        || sub.orbitals.is_some()
        || sub.atoms.is_some()
        || sub.channel.is_some();
    let has_list = !sub.projections.is_empty();
    if has_shorthand && has_list {
        return Err(err(
            1,
            "subspace: use either shorthand keys (species/l/orbitals) or [[subspace.projections]], not both",
        ));
    }
    if has_list {
        return Ok(sub
            .projections
            .iter()
            .map(|p| ProjectionSelector {
                species: p.species.clone(),
                atoms: p.atoms.clone(),
                l: p.l,
                channel: p.channel.clone(),
                orbitals: p.orbitals.clone().unwrap_or_default(),
            })
            .collect());
    }
    let species = sub.species.clone().ok_or_else(|| {
        err(
            1,
            "subspace: provide either 'species' (with 'l' or 'orbitals') or [[subspace.projections]]",
        )
    })?;
    Ok(vec![ProjectionSelector {
        species,
        atoms: sub.atoms.clone(),
        l: sub.l,
        channel: sub.channel.clone(),
        orbitals: sub.orbitals.clone().unwrap_or_default(),
    }])
}

fn projector_fingerprint(ctx: &Context, res: &SubspaceResult, window: BandWindow) -> String {
    let window_json = match window {
        BandWindow::All => json!({"kind": "all"}),
        BandWindow::BandRange { first, last } => {
            json!({"kind": "band_range_1based", "first": first + 1, "last": last + 1})
        }
        BandWindow::EnergyRelative { min_ha, max_ha } => {
            json!({"kind": "window_ha_relative_vbm", "min_bits": min_ha.to_bits(), "max_bits": max_ha.to_bits()})
        }
    };
    let payload = json!({
        "version": "baandi_subspace_projector_v1",
        "method": "atomic-projection",
        "orthonormalization": "lowdin",
        "harmonic_convention": "dft_upf::radial::real_spherical_harmonic(Condon-Shortley,(-1)^m)",
        "k_mesh": res.k_mesh,
        "window": window_json,
        "pseudopotentials": ctx.pp_metadata,
        "orbitals": res.orbitals.iter().map(|o| json!({
            "index": o.index,
            "atom_index": o.atom_index,
            "species": o.species,
            "orbital_name": o.orbital_name,
            "channel_labels": o.channel_labels,
            "l_values": o.l_values,
            "terms": o.terms.iter().map(|t| json!([t.channel_index, t.l, t.m, t.coeff.to_bits()])).collect::<Vec<_>>(),
        })).collect::<Vec<_>>(),
    });
    let canonical = serde_json::to_string(&payload).expect("serializable fingerprint payload");
    sha256(canonical.as_bytes())
}

/// Wannier90 `seedname_hr.dat` text representation (public Wannier90 format
/// specification): 1-based orbital indices, WS degeneracies 15 per line,
/// hoppings in eV via `dft_units`.
fn format_hr_dat_from_real_space(
    hr: &dft_engine::response::subspace::RealSpaceHamiltonian,
    tag: &str,
) -> String {
    let n = hr.num_wann;
    let nrpts = hr.points.len();
    let mut out = String::new();
    writeln!(
        out,
        "written by BAANDI subspace ({tag}); energies in eV, k_mesh={}x{}x{}",
        hr.mesh[0], hr.mesh[1], hr.mesh[2]
    )
    .unwrap();
    writeln!(out, "{n:>12}").unwrap();
    writeln!(out, "{nrpts:>12}").unwrap();
    for (idx, pt) in hr.points.iter().enumerate() {
        write!(out, "{:>5}", pt.degeneracy).unwrap();
        if (idx + 1) % 15 == 0 || idx + 1 == nrpts {
            out.push('\n');
        }
    }
    for (pt, mat) in hr.points.iter().zip(&hr.hoppings) {
        for j in 0..n {
            for i in 0..n {
                let z = mat[(i, j)];
                writeln!(
                    out,
                    "{:>5}{:>5}{:>5}{:>5}{:>5}{:>16.8}{:>16.8}",
                    pt.r[0],
                    pt.r[1],
                    pt.r[2],
                    i + 1,
                    j + 1,
                    ev(z.re),
                    ev(z.im)
                )
                .unwrap();
            }
        }
    }
    out
}

fn format_hr_dat(res: &SubspaceResult) -> String {
    format_hr_dat_from_real_space(&res.real_space, "R4a, ADR 0059")
}

fn format_populations_csv(res: &SubspaceResult) -> String {
    let mut out = String::from(
        "# baandi subspace populations (R4a, ADR 0059); spin-unpolarized (factor 2); 1-based orbital_index, 0-based atom_index\norbital_index,atom_index,species,orbital,l_values,channels,raw_all_bands,raw_in_window,lowdin_in_window\n",
    );
    for (i, orb) in res.orbitals.iter().enumerate() {
        let ls = orb
            .l_values
            .iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join("+");
        let chs = orb.channel_labels.join("+");
        writeln!(
            out,
            "{},{},{},{},{},{},{:.12e},{:.12e},{:.12e}",
            i + 1,
            orb.atom_index,
            orb.species,
            orb.orbital_name,
            ls,
            chs,
            res.raw_populations_all[i],
            res.raw_populations_window[i],
            res.lowdin_populations[i]
        )
        .unwrap();
    }
    out
}

fn format_weights_csv(res: &SubspaceResult) -> String {
    let mut out = String::from(
        "# baandi subspace band weights p_nk = sum_i |U_ni(k)|^2 (R4a, ADR 0059); 1-based band_index; energies in Ha and eV relative to reference_energy_ha\nk_index,k1,k2,k3,band_index,energy_ha,energy_rel_ev,in_window,weight_p_nk\n",
    );
    for (k, (pt, proj)) in res.k_points.iter().zip(&res.projections).enumerate() {
        for (b, &e_ha) in res.eigenvalues[k].iter().enumerate() {
            let in_win = u8::from(proj.window_bands.binary_search(&b).is_ok());
            writeln!(
                out,
                "{},{:.12e},{:.12e},{:.12e},{},{:.12e},{:.12e},{},{:.12e}",
                k,
                pt.fractional[0],
                pt.fractional[1],
                pt.fractional[2],
                b + 1,
                e_ha,
                ev(e_ha - res.reference_energy_ha),
                in_win,
                proj.band_weights[b]
            )
            .unwrap();
        }
    }
    out
}

fn write_interpolated_bands_from_hr(
    ctx: &Context,
    hr: &dft_engine::response::subspace::RealSpaceHamiltonian,
    reference_energy_ha: f64,
    tag: &str,
) -> Result<Option<usize>, RunError> {
    let Some(bands_cfg) = &ctx.config.bands else {
        return Ok(None);
    };
    let vertices = match (&bands_cfg.preset, &bands_cfg.points) {
        (Some(name), None) => kpath::preset(&ctx.cell, name).map_err(|e| err(1, e))?,
        (None, Some(points)) => points
            .iter()
            .map(|p| kpath::Vertex {
                label: p.label.clone(),
                fractional: p.fractional,
                break_before: p.break_before,
            })
            .collect(),
        _ => {
            return Err(err(
                1,
                "supply exactly one of bands.preset and bands.points",
            ));
        }
    };
    if bands_cfg.points_per_segment.is_some() && bands_cfg.spacing_inv_bohr.is_some() {
        return Err(err(1, "choose points_per_segment or spacing_inv_bohr"));
    }
    let path = kpath::generate(
        &ctx.cell,
        &vertices,
        bands_cfg.points_per_segment.unwrap_or(20),
        bands_cfg.spacing_inv_bohr,
    )
    .map_err(|e| err(1, format!("bands path: {e}")))?;
    let mut csv = format!(
        "# baandi subspace Wannier-interpolated bands ({tag}); energies in Ha and eV relative to reference_energy_ha\nk_index,distance_inv_bohr,k1,k2,k3,band_index,energy_ha,energy_rel_ev\n"
    );
    for (k_idx, (&frac, &dist)) in path.fractional.iter().zip(&path.distances).enumerate() {
        let evals = hr.interpolate_eigenvalues(frac).map_err(|e| err(4, e))?;
        for (b, &e_ha) in evals.iter().enumerate() {
            writeln!(
                csv,
                "{},{:.12e},{:.12e},{:.12e},{:.12e},{},{:.12e},{:.12e}",
                k_idx,
                dist,
                frac[0],
                frac[1],
                frac[2],
                b + 1,
                e_ha,
                ev(e_ha - reference_energy_ha)
            )
            .unwrap();
        }
    }
    if dft_cli::parallel::rank() == 0 {
        std::fs::write(ctx.out.join("subspace_bands.csv"), csv)?;
    }
    Ok(Some(path.fractional.len()))
}

fn write_interpolated_bands(
    ctx: &Context,
    res: &SubspaceResult,
) -> Result<Option<usize>, RunError> {
    write_interpolated_bands_from_hr(
        ctx,
        &res.real_space,
        res.reference_energy_ha,
        "R4a, ADR 0059",
    )
}

fn execute_wannier_import(
    input: &Path,
    ctx: &Context,
    section: &dft_cli::config::Subspace,
) -> Result<(), RunError> {
    use dft_engine::response::wannier::{
        DEFAULT_GAUGE_TOLERANCE, WannierImportOptions, subspace_from_wannier,
    };
    use dft_units::{LengthUnit, from_bohr};
    let ang = |bohr: f64| from_bohr(bohr, LengthUnit::Angstrom);
    let ang2 = |bohr2: f64| {
        let s = from_bohr(1.0, LengthUnit::Angstrom);
        bohr2 * s * s
    };

    let k_mesh = section
        .k_mesh
        .or_else(|| ctx.config.response.as_ref().map(|r| r.k_mesh))
        .ok_or_else(|| {
            err(
                1,
                "subspace(wannier90) requires k_mesh (or [response].k_mesh)",
            )
        })?;
    let solver = parse_solver(section.band_solver.as_deref())?;
    let eigen_tolerance = section.eigen_tolerance_ha.unwrap_or(1e-6);
    if !eigen_tolerance.is_finite() || eigen_tolerance <= 0.0 {
        return Err(err(1, "subspace.eigen_tolerance_ha must be positive"));
    }
    let gauge_tolerance = section.gauge_tolerance.unwrap_or(DEFAULT_GAUGE_TOLERANCE);
    if !gauge_tolerance.is_finite() || gauge_tolerance <= 0.0 {
        return Err(err(1, "subspace.gauge_tolerance must be positive"));
    }
    let seedname = section.seedname.as_deref().unwrap_or("baandi");
    let wdir_rel = section.wannier_dir.as_ref().ok_or_else(|| {
        err(
            1,
            "subspace.method = 'wannier90' requires subspace.wannier_dir",
        )
    })?;
    let wdir = input.parent().unwrap_or(Path::new(".")).join(wdir_rel);
    let u_path = wdir.join(format!("{seedname}_u.mat"));
    let u_dis_path = wdir.join(format!("{seedname}_u_dis.mat"));

    let (k_fractional, u_matrices) = crate::wannier::parse_wannier_mat(&u_path)?;
    let u_dis_matrices = if u_dis_path.is_file() {
        let (k_dis, mats) = crate::wannier::parse_wannier_mat(&u_dis_path)?;
        if k_dis.len() != k_fractional.len() {
            return Err(err(
                1,
                "subspace(wannier90): _u_dis.mat k-point count does not match _u.mat",
            ));
        }
        Some(mats)
    } else {
        None
    };

    let scf = ctx.scf()?;
    let potentials: Vec<_> = ctx
        .potentials
        .iter()
        .map(|(label, potential)| dft_engine::SpeciesPotential { label, potential })
        .collect();
    let start = StartInput {
        cell: &ctx.cell,
        atoms: &ctx.atoms,
        potentials: &potentials,
        options: &ctx.options,
    };
    let opts = WannierImportOptions {
        k_mesh,
        solver,
        eigen_tolerance,
        k_fractional,
        u_matrices,
        u_dis_matrices,
        gauge_tolerance,
    };
    let res = subspace_from_wannier(start, &scf.density, &opts).map_err(|e| err(4, e))?;
    let interp_k_points = write_interpolated_bands_from_hr(
        ctx,
        &res.real_space,
        res.reference_energy_ha,
        "R4b Wannier90, ADR 0060",
    )?;

    if dft_cli::parallel::rank() == 0 {
        let u_bytes = std::fs::read(&u_path)?;
        let u_dis_sha = if u_dis_path.is_file() {
            Some(sha256(&std::fs::read(&u_dis_path)?))
        } else {
            None
        };
        let fp_payload = json!({
            "version": "baandi_subspace_projector_v1",
            "method": "wannier90",
            "k_mesh": res.k_mesh,
            "seedname": seedname,
            "u_mat_sha256": sha256(&u_bytes),
            "u_dis_mat_sha256": u_dis_sha,
            "pseudopotentials": ctx.pp_metadata,
        });
        let fingerprint = sha256(serde_json::to_string(&fp_payload).unwrap().as_bytes());

        let mut pop_csv = String::from(
            "# baandi subspace Wannier90 orbital populations (R4b, ADR 0060); spin-unpolarized (factor 2)\norbital_index,population,spread_bohr2,spread_ang2,center_x_ang,center_y_ang,center_z_ang\n",
        );
        for j in 0..res.num_wann {
            let sp = res.spreads.orbital_spreads_bohr2[j];
            let c = res.spreads.centers_bohr[j];
            writeln!(
                pop_csv,
                "{},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",
                j + 1,
                res.populations[j],
                sp,
                ang2(sp),
                ang(c[0]),
                ang(c[1]),
                ang(c[2])
            )
            .unwrap();
        }
        let mut weights_csv = String::from(
            "# baandi subspace Wannier90 band weights p_nk = sum_j |U_tot_nj(k)|^2 (R4b, ADR 0060)\nk_index,k1,k2,k3,band_index,energy_ha,energy_rel_ev,weight_p_nk\n",
        );
        for (k, (pt, proj)) in res.k_points.iter().zip(&res.projections).enumerate() {
            for (b, &e_ha) in res.eigenvalues[k].iter().enumerate() {
                writeln!(
                    weights_csv,
                    "{},{:.12e},{:.12e},{:.12e},{},{:.12e},{:.12e},{:.12e}",
                    k,
                    pt.fractional[0],
                    pt.fractional[1],
                    pt.fractional[2],
                    b + 1,
                    e_ha,
                    ev(e_ha - res.reference_energy_ha),
                    proj.band_weights[b]
                )
                .unwrap();
            }
        }

        let record = json!({
            "kind": "baandi_subspace_r4b_wannier90",
            "schema_version": 1,
            "adr": "0060",
            "method": "wannier90",
            "seedname": seedname,
            "u_mat": u_path.display().to_string(),
            "u_dis_mat": u_dis_path.is_file().then(|| u_dis_path.display().to_string()),
            "projector_fingerprint": fingerprint,
            "input_sha256": ctx.input_hash,
            "pseudopotentials": ctx.pp_metadata,
            "scf": {
                "total_energy_ha": scf.energy,
                "iterations": scf.iterations,
                "path": scf.path,
            },
            "summary": {
                "num_wann": res.num_wann,
                "bands": res.bands,
                "electrons": res.electrons,
                "occupied_bands": res.occupied,
                "reference_energy_ha": res.reference_energy_ha,
                "max_nscf_residual_ha": res.max_residual_ha,
                "max_orthonormality_defect": res.max_orthonormality_defect,
                "wigner_seitz_points": res.real_space.points.len(),
                "hr_hermiticity_defect_ha": res.hr_hermiticity_defect,
                "max_mmn_adjoint_defect": res.max_mmn_adjoint_defect,
                "interpolated_band_kpoints": interp_k_points,
                "total_population": res.populations.iter().sum::<f64>(),
                "spreads": {
                    "omega_i_bohr2": res.spreads.omega_i_bohr2,
                    "omega_tilde_bohr2": res.spreads.omega_tilde_bohr2,
                    "omega_total_bohr2": res.spreads.omega_total_bohr2,
                    "omega_i_ang2": ang2(res.spreads.omega_i_bohr2),
                    "omega_tilde_ang2": ang2(res.spreads.omega_tilde_bohr2),
                    "omega_total_ang2": ang2(res.spreads.omega_total_bohr2),
                },
            },
            "outputs": {
                "json": "subspace.json",
                "hr_dat": "subspace_hr.dat",
                "populations_csv": "subspace_populations.csv",
                "weights_csv": "subspace_weights.csv",
                "bands_csv": interp_k_points.map(|_| "subspace_bands.csv"),
            },
        });
        let text = serde_json::to_string_pretty(&record).expect("json");
        std::fs::write(ctx.out.join("subspace.json"), &text)?;
        std::fs::write(ctx.out.join("results.json"), &text)?;
        std::fs::write(
            ctx.out.join("subspace_hr.dat"),
            format_hr_dat_from_real_space(&res.real_space, "R4b Wannier90, ADR 0060"),
        )?;
        std::fs::write(ctx.out.join("subspace_populations.csv"), pop_csv)?;
        std::fs::write(ctx.out.join("subspace_weights.csv"), weights_csv)?;
        println!(
            "Subspace R4b (wannier90) complete: {} wann, Omega_I={:.4} Ang^2, Omega={:.4} Ang^2, max orth defect {:.3e}",
            res.num_wann,
            ang2(res.spreads.omega_i_bohr2),
            ang2(res.spreads.omega_total_bohr2),
            res.max_orthonormality_defect
        );
    }
    Ok(())
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
        .subspace
        .as_ref()
        .ok_or_else(|| err(1, "subspace requires a [subspace] table"))?;
    if !ctx.selection.reference() {
        return Err(err(1, "subspace requires the reference driver"));
    }
    let method = section.method.as_deref().unwrap_or("atomic-projection");
    if method == "wannier90" {
        return execute_wannier_import(input, &ctx, section);
    }
    if method != "atomic-projection" {
        return Err(err(
            1,
            format!("subspace.method '{method}' must be 'atomic-projection' or 'wannier90'"),
        ));
    }
    let k_mesh = section
        .k_mesh
        .or_else(|| ctx.config.response.as_ref().map(|r| r.k_mesh))
        .ok_or_else(|| err(1, "subspace requires k_mesh (or [response].k_mesh)"))?;
    let solver = parse_solver(section.band_solver.as_deref())?;
    let window = parse_window(section)?;
    let projections = parse_projections(section)?;
    let eigen_tolerance = section.eigen_tolerance_ha.unwrap_or(1e-6);
    if !eigen_tolerance.is_finite() || eigen_tolerance <= 0.0 {
        return Err(err(1, "subspace.eigen_tolerance_ha must be positive"));
    }
    let singular_tolerance = section
        .singular_tolerance
        .unwrap_or(dft_engine::response::subspace::DEFAULT_SINGULAR_TOLERANCE);

    let scf = ctx.scf()?;
    let potentials: Vec<_> = ctx
        .potentials
        .iter()
        .map(|(label, potential)| dft_engine::SpeciesPotential { label, potential })
        .collect();
    let start = StartInput {
        cell: &ctx.cell,
        atoms: &ctx.atoms,
        potentials: &potentials,
        options: &ctx.options,
    };
    let mut opts = SubspaceOptions::new(k_mesh, projections);
    opts.bands = section.bands;
    opts.solver = solver;
    opts.eigen_tolerance = eigen_tolerance;
    opts.window = window;
    opts.singular_tolerance = singular_tolerance;

    let res = subspace(start, &scf.density, &opts).map_err(|e| err(4, e))?;
    let fingerprint = projector_fingerprint(&ctx, &res, window);
    let interp_k_points = write_interpolated_bands(&ctx, &res)?;

    if dft_cli::parallel::rank() == 0 {
        let on_site_ha: Vec<Vec<[f64; 2]>> = res
            .real_space
            .points
            .iter()
            .position(|p| p.r == [0, 0, 0])
            .map(|idx| {
                let h0 = &res.real_space.hoppings[idx];
                (0..res.real_space.num_wann)
                    .map(|i| {
                        (0..res.real_space.num_wann)
                            .map(|j| [h0[(i, j)].re, h0[(i, j)].im])
                            .collect()
                    })
                    .collect()
            })
            .unwrap_or_default();

        let record = json!({
            "kind": "baandi_subspace_r4a",
            "schema_version": 1,
            "adr": "0059",
            "method": "atomic-projection",
            "orthonormalization": "lowdin",
            "harmonic_convention": "dft_upf::radial::real_spherical_harmonic(Condon-Shortley,(-1)^m)",
            "projector_fingerprint": fingerprint,
            "input_sha256": ctx.input_hash,
            "pseudopotentials": ctx.pp_metadata,
            "scf": {
                "total_energy_ha": scf.energy,
                "iterations": scf.iterations,
                "path": scf.path,
            },
            "settings": {
                "k_mesh": res.k_mesh,
                "k_points": res.k_points.len(),
                "bands": res.bands,
                "band_solver": solver.name(),
                "eigen_tolerance_ha": eigen_tolerance,
                "singular_tolerance": singular_tolerance,
                "window": match window {
                    BandWindow::All => json!({"kind": "all"}),
                    BandWindow::BandRange { first, last } => {
                        json!({"kind": "band_range", "first_1based": first + 1, "last_1based": last + 1})
                    }
                    BandWindow::EnergyRelative { min_ha, max_ha } => {
                        json!({"kind": "window_ha", "min_ha": min_ha, "max_ha": max_ha})
                    }
                },
            },
            "summary": {
                "num_wann": res.orbitals.len(),
                "electrons": res.electrons,
                "occupied_bands": res.occupied,
                "reference_energy_ha": res.reference_energy_ha,
                "valence_maximum_ha": res.valence_maximum_ha,
                "conduction_minimum_ha": res.conduction_minimum_ha,
                "max_nscf_residual_ha": res.max_residual_ha,
                "basis_sizes": [res.basis_sizes.0, res.basis_sizes.1],
                "max_orthonormality_defect": res.max_orthonormality_defect,
                "min_overlap_eigenvalue": res.min_overlap_eigenvalue,
                "mean_spilling": res.mean_spilling,
                "wigner_seitz_points": res.real_space.points.len(),
                "hr_hermiticity_defect_ha": res.hr_hermiticity_defect,
                "interpolated_band_kpoints": interp_k_points,
                "total_lowdin_population": res.lowdin_populations.iter().sum::<f64>(),
                "total_raw_population_window": res.raw_populations_window.iter().sum::<f64>(),
                "total_raw_population_all_bands": res.raw_populations_all.iter().sum::<f64>(),
            },
            "orbitals": res.orbitals.iter().enumerate().map(|(i, o)| json!({
                "index_1based": i + 1,
                "atom_index": o.atom_index,
                "species": o.species,
                "position_fractional": o.position_fractional,
                "position_cartesian_bohr": o.position_cartesian_bohr,
                "orbital_name": o.orbital_name,
                "l_values": o.l_values,
                "channel_labels": o.channel_labels,
                "raw_population_all_bands": res.raw_populations_all[i],
                "raw_population_window": res.raw_populations_window[i],
                "lowdin_population": res.lowdin_populations[i],
            })).collect::<Vec<_>>(),
            "k_projections": res.k_points.iter().zip(&res.projections).enumerate().map(|(k, (pt, p))| json!({
                "k_index": k,
                "k_fractional": pt.fractional,
                "weight": pt.weight,
                "window_bands_1based": p.window_bands.iter().map(|b| b + 1).collect::<Vec<_>>(),
                "overlap_eigenvalues": p.overlap_eigenvalues,
                "orthonormality_defect": p.orthonormality_defect,
                "spilling": p.spilling,
                "hermiticity_defect_ha": p.hermiticity_defect,
            })).collect::<Vec<_>>(),
            "on_site_hamiltonian_ha": on_site_ha,
            "outputs": {
                "json": "subspace.json",
                "hr_dat": "subspace_hr.dat",
                "populations_csv": "subspace_populations.csv",
                "weights_csv": "subspace_weights.csv",
                "bands_csv": interp_k_points.map(|_| "subspace_bands.csv"),
            },
        });
        let text = serde_json::to_string_pretty(&record).expect("json");
        std::fs::write(ctx.out.join("subspace.json"), &text)?;
        std::fs::write(ctx.out.join("results.json"), &text)?;
        std::fs::write(ctx.out.join("subspace_hr.dat"), format_hr_dat(&res))?;
        std::fs::write(
            ctx.out.join("subspace_populations.csv"),
            format_populations_csv(&res),
        )?;
        std::fs::write(
            ctx.out.join("subspace_weights.csv"),
            format_weights_csv(&res),
        )?;
        println!(
            "Subspace R4a complete: {} orbitals, {} WS points, max orth defect {:.3e}, min S eig {:.3e}, mean spilling {:.4}",
            res.orbitals.len(),
            res.real_space.points.len(),
            res.max_orthonormality_defect,
            res.min_overlap_eigenvalue,
            res.mean_spilling
        );
    }
    Ok(())
}
