//! `baandi nscf` (R1, `docs/roadmap/rpa-crpa-plan.md`): reference SCF, then a
//! many-band fixed-density solve on the `[nscf]` k points with the selected
//! solver. Outputs (Hartree record, eV views): `nscf.json`,
//! `nscf_eigenvalues.csv`, VASP-layout `EIGENVAL` (absolute eV, not VBM
//! referenced). Optional per-k `nscf_bands_v1` checkpoint units with resume.
use crate::bands::{Context, err, ev};
use dft_cli::{nscf_store::NscfStore, run::RunError, run::sha256};
use dft_engine::{
    parallel::SharedComm,
    refscf::{
        self,
        nscf::{BandOptions, BandSolver, NoStore},
    },
    refstart::StartInput,
};
use dft_geometry::kpoints::KPoint;
use serde_json::json;
use std::{collections::BTreeMap, fmt::Write as _, path::Path};

const OUTPUTS: [&str; 4] = [
    "nscf.json",
    "nscf_eigenvalues.csv",
    "EIGENVAL",
    "results.json",
];

fn gamma_mesh(n: [usize; 3]) -> Result<Vec<KPoint>, RunError> {
    if n.contains(&0) {
        return Err(err(1, "nscf.k_mesh entries must be positive"));
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

pub fn execute(
    input: &Path,
    output: Option<&Path>,
    checkpoint_dir: Option<&Path>,
) -> Result<(), RunError> {
    let ctx = Context::load(
        input,
        output,
        dft_cli::driver_select::Command::Bands,
        &OUTPUTS,
    )?;
    let section = ctx
        .config
        .nscf
        .as_ref()
        .ok_or_else(|| err(1, "nscf requires an [nscf] table"))?;
    if !ctx.selection.reference() {
        return Err(err(1, "nscf requires the reference driver"));
    }
    let name = section.band_solver.as_deref().unwrap_or("dense");
    let solver = BandSolver::parse(name).ok_or_else(|| {
        err(
            1,
            format!("nscf.band_solver = {name:?}: use davidson, dense or chebyshev"),
        )
    })?;
    let tolerance = section.eigen_tolerance_ha.unwrap_or(1e-6);
    let mut band_options = BandOptions::new(section.bands, tolerance, solver);
    band_options.buffer = section.buffer_bands;
    if let Some(d) = section.filter_degree {
        band_options.degree = d;
    }
    if let Some(m) = section.max_iterations {
        band_options.max_iterations = m;
    }
    if solver != BandSolver::Chebyshev
        && (section.buffer_bands.is_some()
            || section.filter_degree.is_some()
            || section.max_iterations.is_some())
    {
        return Err(err(
            1,
            "nscf.buffer_bands/filter_degree/max_iterations apply only to band_solver = \"chebyshev\"",
        ));
    }
    let points = match section.k_mesh {
        Some(n) => gamma_mesh(n)?,
        None => ctx.config.basis.resolve_kpoints().map_err(|e| err(1, e))?,
    };
    let scf = ctx.scf()?;
    let density_bytes: Vec<u8> = scf
        .density
        .iter()
        .flat_map(|s| s.iter().flat_map(|v| v.to_le_bytes()).collect::<Vec<_>>())
        .collect();
    let density_sha = sha256(&density_bytes);
    let ranks = SharedComm::active().map_or(1, |c| dft_mpi::Collective::size(&c));
    println!(
        "NSCF: {} k points, {} bands per channel, solver {}, residual certificate {tolerance} Ha",
        points.len(),
        section.bands,
        solver.name()
    );
    let potentials: Vec<_> = ctx
        .potentials
        .iter()
        .map(|(label, potential)| dft_engine::SpeciesPotential { label, potential })
        .collect();
    let input_state = StartInput {
        cell: &ctx.cell,
        atoms: &ctx.atoms,
        potentials: &potentials,
        options: &ctx.options,
    };
    let mut store = None;
    let result = if let Some(dir) = checkpoint_dir {
        let points_bits: Vec<u8> = points
            .iter()
            .flat_map(|p| {
                p.fractional
                    .iter()
                    .chain(std::iter::once(&p.weight))
                    .flat_map(|v| v.to_le_bytes())
                    .collect::<Vec<_>>()
            })
            .collect();
        let pp = ctx
            .pp_metadata
            .as_array()
            .map(|a| {
                a.iter()
                    .map(|p| p["sha256"].as_str().unwrap_or("").to_owned())
                    .collect::<Vec<_>>()
                    .join(",")
            })
            .unwrap_or_default();
        let fingerprint: BTreeMap<String, String> = [
            ("input_sha256", ctx.input_hash.clone()),
            ("upf_sha256", pp),
            ("density_sha256", density_sha.clone()),
            ("points_sha256", sha256(&points_bits)),
            ("solver", solver.name().to_owned()),
            ("bands", section.bands.to_string()),
            ("tolerance_ha_bits", tolerance.to_bits().to_string()),
            ("buffer_bands", format!("{:?}", band_options.buffer)),
            ("filter_degree", band_options.degree.to_string()),
            ("ranks", ranks.to_string()),
            ("threads", ctx.options.threads.to_string()),
            (
                "platform",
                format!("{}-{}", std::env::consts::ARCH, std::env::consts::OS),
            ),
            (
                "backend",
                if dft_numerics::MKL_BACKEND {
                    "mkl"
                } else {
                    "faer"
                }
                .to_owned(),
            ),
            ("code_version", env!("CARGO_PKG_VERSION").to_owned()),
        ]
        .into_iter()
        .map(|(k, v)| (k.to_owned(), v))
        .collect();
        let mut s = NscfStore::new(dir, fingerprint).map_err(|e| err(5, e))?;
        let r = refscf::nscf::nscf_bands(input_state, &scf.density, &points, band_options, &mut s);
        store = Some(s);
        r
    } else {
        refscf::nscf::nscf_bands(
            input_state,
            &scf.density,
            &points,
            band_options,
            &mut NoStore,
        )
    }
    .map_err(|e| err(4, e))?;
    if let Some(s) = &store {
        let units: Vec<(usize, usize)> = (0..result.channels.len())
            .flat_map(|c| (0..points.len()).map(move |k| (c, k)))
            .collect();
        match SharedComm::active() {
            Some(c) => s.write_manifest(&c, &units),
            None => s.write_manifest(&refscf::metric::Serial, &units),
        }
        .map_err(|e| err(5, e))?;
        println!(
            "NSCF checkpoint: {} units resumed (re-certified), {} written",
            s.resumed.len(),
            s.written.len()
        );
    }
    let eigenvalues: Vec<Vec<Vec<f64>>> = result
        .channels
        .iter()
        .map(|c| c.iter().map(|b| b.eigenvalues.clone()).collect())
        .collect();
    let residuals: Vec<Vec<Vec<f64>>> = result
        .channels
        .iter()
        .map(|c| c.iter().map(|b| b.residuals.clone()).collect())
        .collect();
    let max_residual = residuals
        .iter()
        .flatten()
        .flatten()
        .copied()
        .fold(0., f64::max);
    let polarized = eigenvalues.len() == 2;
    let mut csv = String::from(
        "# baandi nscf (R1): absolute Kohn-Sham eigenvalues of the fixed SCF density; energy_ha Hartree, energy_ev eV, residual_ha = fresh ||H psi - e psi|| (Ha)\nchannel,k_index,kx,ky,kz,weight,band,energy_ha,energy_ev,residual_ha\n",
    );
    for (c, channel) in eigenvalues.iter().enumerate() {
        for (k, values) in channel.iter().enumerate() {
            let p = &points[k];
            for (b, e) in values.iter().enumerate() {
                writeln!(
                    csv,
                    "{c},{k},{:.16e},{:.16e},{:.16e},{:.16e},{},{e:.16e},{:.16e},{:.6e}",
                    p.fractional[0],
                    p.fractional[1],
                    p.fractional[2],
                    p.weight,
                    b + 1,
                    ev(*e),
                    residuals[c][k][b]
                )
                .unwrap();
            }
        }
    }
    let mut eigenval = ctx.vasp_header(true, polarized);
    writeln!(
        eigenval,
        "{:.0} {} {}",
        result.electrons,
        points.len(),
        section.bands
    )
    .unwrap();
    for (k, p) in points.iter().enumerate() {
        writeln!(
            eigenval,
            "\n{:.16e} {:.16e} {:.16e} {:.16e}",
            p.fractional[0], p.fractional[1], p.fractional[2], p.weight
        )
        .unwrap();
        #[allow(clippy::needless_range_loop)]
        for b in 0..section.bands {
            if polarized {
                writeln!(
                    eigenval,
                    "{} {:.16e} {:.16e} {:.1} {:.1}",
                    b + 1,
                    ev(eigenvalues[0][k][b]),
                    ev(eigenvalues[1][k][b]),
                    if b < scf.occupied[0] { 1.0 } else { 0.0 },
                    if b < scf.occupied[1] { 1.0 } else { 0.0 }
                )
                .unwrap();
            } else {
                writeln!(
                    eigenval,
                    "{} {:.16e} {:.1}",
                    b + 1,
                    ev(eigenvalues[0][k][b]),
                    if b < scf.occupied[0] { 2.0 } else { 0.0 }
                )
                .unwrap();
            }
        }
    }
    let mut meta = json!({
        "schema_version": 1,
        "kind": "baandi_nscf_r1",
        "code": "rust-dft",
        "code_version": env!("CARGO_PKG_VERSION"),
        "title": ctx.config.system.name,
        "input_sha256": ctx.input_hash,
        "pseudopotentials": ctx.pp_metadata,
        "scf_density_sha256": density_sha,
        "scf_energy_hartree": scf.energy,
        "scf_iterations": scf.iterations,
        "scf_path": scf.path,
        "functional": ctx.options.functional.name(),
        "cutoff_hartree": ctx.options.cutoff_hartree,
        "threads": ctx.options.threads,
        "ranks": ranks,
        "electrons": result.electrons,
        "occupied_bands": scf.occupied,
        "nbands": section.bands,
        "band_solver": solver.name(),
        "eigen_tolerance_hartree": tolerance,
        "max_residual_hartree": max_residual,
        "energy_reference": "absolute Kohn-Sham eigenvalues (Hartree), BAANDI G=0 convention; not VBM referenced",
        "hartree_to_ev": ev(1.0),
        "k_mesh": section.k_mesh,
        "k_sampling": if section.k_mesh.is_some() { "Gamma-centered full mesh, no symmetry reduction" } else { "[basis] k points" },
        "lattice_columns_bohr": ctx.cell.lattice_columns_bohr(),
        "atoms": ctx.atoms.iter().map(|a| json!({"species": a.species, "fractional": a.fractional_position})).collect::<Vec<_>>(),
        "kpoints": points.iter().map(|p| json!({"fractional": p.fractional, "weight": p.weight})).collect::<Vec<_>>(),
        "stats": result.stats.iter().map(|c| c.iter().map(|s| json!({
            "basis_size": s.basis_size, "resumed": s.resumed,
            "iterations": s.iterations, "h_columns": s.h_columns})).collect::<Vec<_>>()).collect::<Vec<_>>(),
        "band_data": {"ordering": "channel_k_band",
            "shape": [eigenvalues.len(), points.len(), section.bands],
            "eigenvalues_ha": eigenvalues, "residuals_ha": residuals},
        "checkpoint": store.as_ref().map(|s| json!({
            "directory": checkpoint_dir.map(|d| d.display().to_string()),
            "format": "nscf_bands_v1 (docs/checkpoint-format.md)",
            "resumed_units": s.resumed, "written_units": s.written})),
    });
    ctx.selection.record(&mut meta);
    ctx.write("nscf_eigenvalues.csv", &csv)?;
    ctx.write("EIGENVAL", &eigenval)?;
    ctx.json("nscf.json", &meta)?;
    ctx.json("results.json", &meta)?;
    println!(
        "NSCF done: max certified residual {max_residual:.3e} Ha; outputs in {}",
        ctx.out.display()
    );
    Ok(())
}
