//! `baandi wannier-export` (R4b, `docs/roadmap/rpa-crpa-plan.md` 2.8, 2.11, 3.6;
//! ADR 0060) plus parsers for Wannier90 `.nnkp`, `_u.mat` and `_u_dis.mat`.
//!
//! File layouts follow the public Wannier90 3.x format documentation (no GPL
//! source read or copied; MIT/GPL boundary per decision D2). All unit
//! conversions go through `dft_units`.
use crate::bands::{Context, err, ev};
use crate::subspace::{parse_projections, parse_solver};
use dft_cli::run::RunError;
use dft_engine::{
    refstart::StartInput,
    response::wannier::{KNeighbor, WannierExportOptions, WannierExportResult, wannier_export},
};
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::data::Complex64;
use dft_units::{EnergyUnit, LengthUnit, from_bohr, from_hartree};
use serde_json::json;
use std::{
    fmt::Write as _,
    fs,
    path::{Path, PathBuf},
};

fn ang(bohr: f64) -> f64 {
    from_bohr(bohr, LengthUnit::Angstrom)
}

fn ang2(bohr2: f64) -> f64 {
    let s = from_bohr(1.0, LengthUnit::Angstrom);
    bohr2 * s * s
}

/// Parse the `begin nnkpts` ... `end nnkpts` block of a Wannier90 `.nnkp` file.
pub fn parse_nnkp(path: &Path, expected_nk: usize) -> Result<Vec<Vec<KNeighbor>>, RunError> {
    let text =
        fs::read_to_string(path).map_err(|e| err(1, format!("nnkp {}: {e}", path.display())))?;
    let mut in_block = false;
    let mut nntot: Option<usize> = None;
    let mut by_k: Vec<Vec<KNeighbor>> = vec![Vec::new(); expected_nk];
    let mut found_block = false;

    for raw_line in text.lines() {
        let line = raw_line.split('!').next().unwrap_or("").trim();
        if line.is_empty() {
            continue;
        }
        let lower = line.to_ascii_lowercase();
        if lower == "begin nnkpts" {
            in_block = true;
            found_block = true;
            continue;
        }
        if lower == "end nnkpts" {
            break;
        }
        if !in_block {
            continue;
        }
        if nntot.is_none() {
            let n: usize = line
                .parse()
                .map_err(|_| err(1, format!("nnkp: invalid nntot line '{line}'")))?;
            if n == 0 {
                return Err(err(1, "nnkp: nntot must be positive"));
            }
            nntot = Some(n);
            continue;
        }
        let toks: Vec<&str> = line.split_whitespace().collect();
        if toks.len() != 5 {
            return Err(err(
                1,
                format!("nnkp: expected 5 integers per neighbor row, got '{line}'"),
            ));
        }
        let ik: usize = toks[0]
            .parse()
            .map_err(|_| err(1, format!("nnkp: bad ik '{line}'")))?;
        let kn: usize = toks[1]
            .parse()
            .map_err(|_| err(1, format!("nnkp: bad neighbor ik '{line}'")))?;
        let g0: i32 = toks[2]
            .parse()
            .map_err(|_| err(1, format!("nnkp: bad g1 '{line}'")))?;
        let g1: i32 = toks[3]
            .parse()
            .map_err(|_| err(1, format!("nnkp: bad g2 '{line}'")))?;
        let g2: i32 = toks[4]
            .parse()
            .map_err(|_| err(1, format!("nnkp: bad g3 '{line}'")))?;
        if ik == 0 || ik > expected_nk || kn == 0 || kn > expected_nk {
            return Err(err(
                1,
                format!("nnkp: k-point index out of 1..={expected_nk} in '{line}'"),
            ));
        }
        by_k[ik - 1].push(KNeighbor {
            k_index: ik - 1,
            neighbor_k_index: kn - 1,
            g_shift: [g0, g1, g2],
        });
    }
    if !found_block {
        return Err(err(
            1,
            format!("nnkp {}: missing 'begin nnkpts' block", path.display()),
        ));
    }
    let expected_nn = nntot.ok_or_else(|| err(1, "nnkp: empty nnkpts block"))?;
    for (k, row) in by_k.iter().enumerate() {
        if row.len() != expected_nn {
            return Err(err(
                1,
                format!(
                    "nnkp: k-point {} has {} neighbors, expected {expected_nn}",
                    k + 1,
                    row.len()
                ),
            ));
        }
    }
    Ok(by_k)
}

/// Parse a Wannier90 matrix file (`_u.mat` or `_u_dis.mat`):
/// - line 1: comment header
/// - line 2: `num_kpts cols rows` (`num_kpts num_wann num_wann` for `_u.mat`,
///   `num_kpts num_wann num_bands` for `_u_dis.mat`)
/// - per k-point: fractional k-point line (`kx ky kz`), followed by `rows * cols`
///   lines of `Re Im` in column-major order (`row` fastest, `col` slower).
pub fn parse_wannier_mat(path: &Path) -> Result<(Vec<[f64; 3]>, Vec<ComplexMatrix>), RunError> {
    let text = fs::read_to_string(path)
        .map_err(|e| err(1, format!("wannier mat {}: {e}", path.display())))?;
    let mut lines = text.lines();
    let _header = lines
        .next()
        .ok_or_else(|| err(1, format!("wannier mat {}: empty file", path.display())))?;
    let dim_line = lines.find(|l| !l.trim().is_empty()).ok_or_else(|| {
        err(
            1,
            format!("wannier mat {}: missing dimensions", path.display()),
        )
    })?;
    let dims: Vec<usize> = dim_line
        .split_whitespace()
        .map(|t| {
            t.parse::<usize>()
                .map_err(|_| err(1, format!("wannier mat: invalid dimension '{t}'")))
        })
        .collect::<Result<_, _>>()?;
    if dims.len() != 3 || dims.contains(&0) {
        return Err(err(
            1,
            format!("wannier mat: expected 3 positive integers, got '{dim_line}'"),
        ));
    }
    let (nk, cols, rows) = (dims[0], dims[1], dims[2]);
    let mut ks = Vec::with_capacity(nk);
    let mut mats = Vec::with_capacity(nk);
    for k in 0..nk {
        let k_line = lines.find(|l| !l.trim().is_empty()).ok_or_else(|| {
            err(
                1,
                format!(
                    "wannier mat {}: truncated at k-point {}",
                    path.display(),
                    k + 1
                ),
            )
        })?;
        let kt: Vec<f64> = k_line
            .split_whitespace()
            .map(|t| {
                t.parse::<f64>()
                    .map_err(|_| err(1, format!("wannier mat: bad k coord '{t}'")))
            })
            .collect::<Result<_, _>>()?;
        if kt.len() != 3 || kt.iter().any(|x| !x.is_finite()) {
            return Err(err(
                1,
                format!("wannier mat: expected 3 finite k coords, got '{k_line}'"),
            ));
        }
        ks.push([kt[0], kt[1], kt[2]]);
        let mut mat = ComplexMatrix::zeros(rows, cols).map_err(|e| err(4, e))?;
        for col in 0..cols {
            for row in 0..rows {
                let entry = lines.find(|l| !l.trim().is_empty()).ok_or_else(|| {
                    err(
                        1,
                        format!(
                            "wannier mat {}: truncated inside k-point {}",
                            path.display(),
                            k + 1
                        ),
                    )
                })?;
                let vals: Vec<f64> = entry
                    .split_whitespace()
                    .map(|t| {
                        t.parse::<f64>()
                            .map_err(|_| err(1, format!("wannier mat: bad float '{t}'")))
                    })
                    .collect::<Result<_, _>>()?;
                if vals.len() != 2 || !vals[0].is_finite() || !vals[1].is_finite() {
                    return Err(err(
                        1,
                        format!("wannier mat: expected 're im', got '{entry}'"),
                    ));
                }
                mat[(row, col)] = Complex64::new(vals[0], vals[1]);
            }
        }
        mats.push(mat);
    }
    Ok((ks, mats))
}

fn format_win(ctx: &Context, res: &WannierExportResult, sub: &dft_cli::config::Subspace) -> String {
    let mut out = String::new();
    writeln!(out, "! Generated by BAANDI wannier-export (R4b, ADR 0060)").unwrap();
    writeln!(out, "num_bands = {}", res.bands).unwrap();
    writeln!(out, "num_wann  = {}", res.orbitals.len()).unwrap();
    if let Some([emin, emax]) = sub.window_ha {
        writeln!(
            out,
            "dis_win_min = {:.10}",
            from_hartree(res.valence_maximum_ha + emin, EnergyUnit::ElectronVolt)
        )
        .unwrap();
        writeln!(
            out,
            "dis_win_max = {:.10}",
            from_hartree(res.valence_maximum_ha + emax, EnergyUnit::ElectronVolt)
        )
        .unwrap();
    }
    writeln!(
        out,
        "mp_grid = {} {} {}",
        res.k_mesh[0], res.k_mesh[1], res.k_mesh[2]
    )
    .unwrap();
    writeln!(out, "write_u_matrices = .true.").unwrap();
    writeln!(out, "write_hr = .true.").unwrap();
    writeln!(out).unwrap();
    writeln!(out, "begin unit_cell_cart").unwrap();
    writeln!(out, "Ang").unwrap();
    let cols = ctx.cell.lattice_columns_bohr();
    let vecs: [[f64; 3]; 3] = std::array::from_fn(|j| [cols[0][j], cols[1][j], cols[2][j]]);
    for v in vecs {
        let [x, y, z] = v.map(ang);
        writeln!(out, "  {x:>16.10} {y:>16.10} {z:>16.10}").unwrap();
    }
    writeln!(out, "end unit_cell_cart").unwrap();
    writeln!(out).unwrap();
    writeln!(out, "begin atoms_frac").unwrap();
    for atom in &ctx.atoms {
        writeln!(
            out,
            "  {:<4} {:>16.10} {:>16.10} {:>16.10}",
            atom.species,
            atom.fractional_position[0],
            atom.fractional_position[1],
            atom.fractional_position[2]
        )
        .unwrap();
    }
    writeln!(out, "end atoms_frac").unwrap();
    writeln!(out).unwrap();
    writeln!(out, "begin projections").unwrap();
    for orb in &res.orbitals {
        writeln!(
            out,
            "  f={:.8},{:.8},{:.8}:{} ! atom {} ({})",
            orb.position_fractional[0],
            orb.position_fractional[1],
            orb.position_fractional[2],
            orb.orbital_name,
            orb.atom_index + 1,
            orb.species
        )
        .unwrap();
    }
    writeln!(out, "end projections").unwrap();
    writeln!(out).unwrap();
    writeln!(out, "begin kpoints").unwrap();
    for pt in &res.k_points {
        writeln!(
            out,
            "  {:>14.10} {:>14.10} {:>14.10}",
            pt.fractional[0], pt.fractional[1], pt.fractional[2]
        )
        .unwrap();
    }
    writeln!(out, "end kpoints").unwrap();
    out
}

fn format_eig(res: &WannierExportResult) -> String {
    let mut out = String::new();
    for (k, evals) in res.eigenvalues.iter().enumerate() {
        for (m, &e_ha) in evals.iter().enumerate() {
            writeln!(out, "{:>5}{:>5}{:>18.12}", m + 1, k + 1, ev(e_ha)).unwrap();
        }
    }
    out
}

fn format_amn(res: &WannierExportResult) -> String {
    let mut out = String::new();
    let nb = res.bands;
    let nk = res.k_points.len();
    let nw = res.orbitals.len();
    writeln!(out, "Created by BAANDI wannier-export (R4b, ADR 0060)").unwrap();
    writeln!(out, "{nb:>12}{nk:>12}{nw:>12}").unwrap();
    for (k, a) in res.amn.iter().enumerate() {
        for n in 0..nw {
            for m in 0..nb {
                let z = a[(m, n)];
                writeln!(
                    out,
                    "{:>5}{:>5}{:>5}{:>18.12}{:>18.12}",
                    m + 1,
                    n + 1,
                    k + 1,
                    z.re,
                    z.im
                )
                .unwrap();
            }
        }
    }
    out
}

fn format_mmn(res: &WannierExportResult) -> String {
    let mut out = String::new();
    let nb = res.bands;
    let nk = res.k_points.len();
    let nntot = res.neighbors.first().map_or(0, Vec::len);
    writeln!(out, "Created by BAANDI wannier-export (R4b, ADR 0060)").unwrap();
    writeln!(out, "{nb:>12}{nk:>12}{nntot:>12}").unwrap();
    for (k, (neighs, mats)) in res.neighbors.iter().zip(&res.mmn).enumerate() {
        for (nb_entry, mat) in neighs.iter().zip(mats) {
            writeln!(
                out,
                "{:>5}{:>5}{:>5}{:>5}{:>5}",
                k + 1,
                nb_entry.neighbor_k_index + 1,
                nb_entry.g_shift[0],
                nb_entry.g_shift[1],
                nb_entry.g_shift[2]
            )
            .unwrap();
            for n in 0..nb {
                for m in 0..nb {
                    let z = mat[(m, n)];
                    writeln!(out, "{:>18.12}{:>18.12}", z.re, z.im).unwrap();
                }
            }
        }
    }
    out
}

pub fn execute(
    input: &Path,
    output: Option<&Path>,
    seedname_override: Option<&str>,
    nnkp_override: Option<&Path>,
) -> Result<(), RunError> {
    let base_outputs = ["wannier_export.json", "results.json"];
    let ctx = Context::load(
        input,
        output,
        dft_cli::driver_select::Command::Bands,
        &base_outputs,
    )?;
    let section = ctx
        .config
        .subspace
        .as_ref()
        .ok_or_else(|| err(1, "wannier-export requires a [subspace] table"))?;
    if !ctx.selection.reference() {
        return Err(err(1, "wannier-export requires the reference driver"));
    }
    let seedname = seedname_override
        .or(section.seedname.as_deref())
        .unwrap_or("baandi");
    if seedname.is_empty()
        || seedname.contains('/')
        || seedname.contains('\\')
        || seedname == "."
        || seedname == ".."
    {
        return Err(err(
            1,
            format!("wannier-export: invalid seedname '{seedname}'"),
        ));
    }

    let k_mesh = section
        .k_mesh
        .or_else(|| ctx.config.response.as_ref().map(|r| r.k_mesh))
        .ok_or_else(|| err(1, "wannier-export requires k_mesh (or [response].k_mesh)"))?;
    let nk = k_mesh[0] * k_mesh[1] * k_mesh[2];
    if nk == 0 {
        return Err(err(1, "wannier-export: k_mesh entries must be positive"));
    }
    let solver = parse_solver(section.band_solver.as_deref())?;
    let projections = parse_projections(section)?;
    let eigen_tolerance = section.eigen_tolerance_ha.unwrap_or(1e-6);
    if !eigen_tolerance.is_finite() || eigen_tolerance <= 0.0 {
        return Err(err(1, "subspace.eigen_tolerance_ha must be positive"));
    }

    let nnkp_path: Option<PathBuf> = nnkp_override.map(Path::to_owned).or_else(|| {
        section
            .nnkp
            .as_ref()
            .map(|p| input.parent().unwrap_or(Path::new(".")).join(p))
    });
    let explicit_neighbors = match &nnkp_path {
        Some(p) => Some(parse_nnkp(p, nk)?),
        None => None,
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
    let mut opts = WannierExportOptions::new(k_mesh, projections);
    opts.bands = section.bands;
    opts.solver = solver;
    opts.eigen_tolerance = eigen_tolerance;
    opts.explicit_neighbors = explicit_neighbors;

    let res = wannier_export(start, &scf.density, &opts).map_err(|e| err(4, e))?;

    if dft_cli::parallel::rank() == 0 {
        let win_name = format!("{seedname}.win");
        let eig_name = format!("{seedname}.eig");
        let amn_name = format!("{seedname}.amn");
        let mmn_name = format!("{seedname}.mmn");
        fs::write(ctx.out.join(&win_name), format_win(&ctx, &res, section))?;
        fs::write(ctx.out.join(&eig_name), format_eig(&res))?;
        fs::write(ctx.out.join(&amn_name), format_amn(&res))?;
        fs::write(ctx.out.join(&mmn_name), format_mmn(&res))?;

        let record = json!({
            "kind": "baandi_wannier_export_r4b",
            "schema_version": 1,
            "adr": "0060",
            "seedname": seedname,
            "nnkp_source": nnkp_path.as_ref().map(|p| p.display().to_string()),
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
                "num_bands": res.bands,
                "num_wann": res.orbitals.len(),
                "nntot": res.neighbors.first().map_or(0, Vec::len),
                "band_solver": solver.name(),
                "eigen_tolerance_ha": eigen_tolerance,
            },
            "summary": {
                "electrons": res.electrons,
                "occupied_bands": res.occupied,
                "valence_maximum_ha": res.valence_maximum_ha,
                "conduction_minimum_ha": res.conduction_minimum_ha,
                "max_nscf_residual_ha": res.max_residual_ha,
                "max_mmn_adjoint_defect": res.max_adjoint_defect,
                "b1_shells": res.bset.as_ref().map(|b| b.shells),
                "b1_completeness_residual": res.bset.as_ref().map(|b| b.completeness_residual),
                "initial_spreads": res.initial_spreads.as_ref().map(|s| json!({
                    "omega_i_bohr2": s.omega_i_bohr2,
                    "omega_tilde_bohr2": s.omega_tilde_bohr2,
                    "omega_total_bohr2": s.omega_total_bohr2,
                    "omega_i_ang2": ang2(s.omega_i_bohr2),
                    "omega_tilde_ang2": ang2(s.omega_tilde_bohr2),
                    "omega_total_ang2": ang2(s.omega_total_bohr2),
                    "orbital_spreads_ang2": s.orbital_spreads_bohr2.iter().map(|&v| ang2(v)).collect::<Vec<_>>(),
                    "centers_ang": s.centers_bohr.iter().map(|c| [ang(c[0]), ang(c[1]), ang(c[2])]).collect::<Vec<_>>(),
                })),
            },
            "outputs": {
                "win": win_name,
                "eig": eig_name,
                "amn": amn_name,
                "mmn": mmn_name,
                "json": "wannier_export.json",
            },
        });
        let text = serde_json::to_string_pretty(&record).expect("json");
        fs::write(ctx.out.join("wannier_export.json"), &text)?;
        fs::write(ctx.out.join("results.json"), &text)?;
        println!(
            "Wannier90 export R4b complete: seedname={seedname}, {} bands, {} wann, {} neighbors/k, max M_mn adjoint defect {:.3e}",
            res.bands,
            res.orbitals.len(),
            res.neighbors.first().map_or(0, Vec::len),
            res.max_adjoint_defect
        );
    }
    Ok(())
}
