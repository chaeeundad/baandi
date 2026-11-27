//! `baandi rpa` (R3, `docs/roadmap/rpa-crpa-plan.md` 2.5, 2.11, 3.6; ADR
//! 0055): reference SCF, then the RPA correlation energy of
//! `dft_engine::response::rpa` on the `[response]` Gamma-centred k mesh with
//! the `[rpa]` cutoff series. Outputs: `rpa.json` (complete record, units,
//! settings, per-q and per-cutoff values, direction sets, extrapolation),
//! `rpa_q.csv` (one row per q and cutoff), `rpa_extrapolation.csv` (one row
//! per cutoff with the fit), `rpa_frequency.csv` (q-weighted integrand per
//! frequency and cutoff), `results.json` (copy of `rpa.json`). Optional per-q
//! checkpoint records `rpa_q_v1` (`--checkpoint-dir` / `task.checkpoint_dir`).
use crate::bands::{Context, err, ev};
use dft_cli::run::{RunError, sha256};
use dft_engine::{
    ScfError,
    parallel::SharedComm,
    refscf::{metric::Serial, nscf::BandSolver},
    refstart::StartInput,
    response::rpa::{QContribution, QCutoff, QStore, QSymmetry, RpaOptions, rpa},
};
use serde_json::{Value, json};
use std::{
    fmt::Write as _,
    path::{Path, PathBuf},
};

const OUTPUTS: [&str; 5] = [
    "rpa.json",
    "rpa_q.csv",
    "rpa_extrapolation.csv",
    "rpa_frequency.csv",
    "results.json",
];
/// Per-q JSON records in a directory; floats stored as IEEE bits (exact
/// resume) next to readable values. Rank 0 writes; every rank reads.
struct RpaStore {
    dir: Option<PathBuf>,
    fingerprint: String,
}
impl RpaStore {
    fn path(&self, index: usize) -> Option<PathBuf> {
        self.dir
            .as_ref()
            .map(|d| d.join(format!("rpa_q_{index:05}.json")))
    }
    fn load_record(&self, index: usize) -> Result<Option<QContribution>, ScfError> {
        let Some(path) = self.path(index) else {
            return Ok(None);
        };
        let Ok(text) = std::fs::read_to_string(&path) else {
            return Ok(None);
        };
        let bad =
            |what: &str| ScfError::Numerical(format!("rpa checkpoint {}: {what}", path.display()));
        let v: Value = serde_json::from_str(&text).map_err(|_| bad("not JSON"))?;
        if v["kind"] != "rpa_q_v1" {
            return Err(bad("kind is not rpa_q_v1"));
        }
        if v["fingerprint"] != self.fingerprint.as_str() {
            return Err(bad("fingerprint does not match this input and density"));
        }
        let f = |x: &Value| {
            x.as_u64()
                .map(f64::from_bits)
                .ok_or_else(|| bad("float bits"))
        };
        let fs = |x: &Value| -> Result<Vec<f64>, ScfError> {
            x.as_array()
                .ok_or_else(|| bad("array"))?
                .iter()
                .map(f)
                .collect()
        };
        let u = |x: &Value| x.as_u64().ok_or_else(|| bad("integer"));
        let bits = &v["bits"];
        let q = fs(&bits["fractional"])?;
        let mut cutoffs = Vec::new();
        for c in bits["cutoffs"].as_array().ok_or_else(|| bad("cutoffs"))? {
            let mut direction_energies = Vec::new();
            for d in c["direction_energies"].as_array().unwrap_or(&Vec::new()) {
                direction_energies.push((
                    d[0].as_str()
                        .ok_or_else(|| bad("direction label"))?
                        .to_owned(),
                    f(&d[1])?,
                ));
            }
            cutoffs.push(QCutoff {
                response_size: u(&c["response_size"])? as usize,
                terms: fs(&c["terms"])?,
                energy: f(&c["energy"])?,
                min_pivot: f(&c["min_pivot"])?,
                hermiticity: f(&c["hermiticity"])?,
                direction_energies,
                body_only_energy: match &c["body_only_energy"] {
                    Value::Null => None,
                    x => Some(f(x)?),
                },
            });
        }
        if q.len() != 3 {
            return Err(bad("q"));
        }
        Ok(Some(QContribution {
            index: u(&bits["index"])? as usize,
            fractional: [q[0], q[1], q[2]],
            weight: f(&bits["weight"])?,
            multiplicity: u(&bits["multiplicity"])? as usize,
            optical: bits["optical"].as_bool().ok_or_else(|| bad("optical"))?,
            transitions: u(&bits["transitions"])?,
            omega0: f(&bits["omega0"])?,
            cutoffs,
        }))
    }
    fn save_record(&self, record: &QContribution) -> Result<(), ScfError> {
        let Some(path) = self.path(record.index) else {
            return Ok(());
        };
        if dft_cli::parallel::rank() != 0 {
            return Ok(());
        }
        let b = |x: f64| json!(x.to_bits());
        let value = json!({
            "kind": "rpa_q_v1",
            "fingerprint": self.fingerprint,
            "note": "per-q RPA correlation contribution (ADR 0055); 'bits' holds IEEE-754 bit patterns for exact resume, 'readable' the same values in Hartree",
            "readable": {
                "q_fractional": record.fractional,
                "weight": record.weight,
                "optical": record.optical,
                "e_c_q_hartree_per_cutoff": record.cutoffs.iter().map(|c| c.energy).collect::<Vec<_>>(),
            },
            "bits": {
                "index": record.index,
                "fractional": record.fractional.map(f64::to_bits),
                "weight": b(record.weight),
                "multiplicity": record.multiplicity,
                "optical": record.optical,
                "transitions": record.transitions,
                "omega0": b(record.omega0),
                "cutoffs": record.cutoffs.iter().map(|c| json!({
                    "response_size": c.response_size,
                    "terms": c.terms.iter().map(|t| t.to_bits()).collect::<Vec<_>>(),
                    "energy": b(c.energy),
                    "min_pivot": b(c.min_pivot),
                    "hermiticity": b(c.hermiticity),
                    "direction_energies": c.direction_energies.iter().map(|(l, e)| json!([l, e.to_bits()])).collect::<Vec<_>>(),
                    "body_only_energy": c.body_only_energy.map(f64::to_bits),
                })).collect::<Vec<_>>(),
            },
        });
        let io = |e: std::io::Error| ScfError::Numerical(format!("rpa checkpoint write: {e}"));
        let dir = path.parent().expect("record in a directory");
        std::fs::create_dir_all(dir).map_err(io)?;
        let tmp = path.with_extension("json.tmp");
        std::fs::write(
            &tmp,
            serde_json::to_string_pretty(&value).expect("serializable"),
        )
        .map_err(io)?;
        std::fs::rename(&tmp, &path).map_err(io)
    }
}
impl QStore<SharedComm> for RpaStore {
    fn load(&mut self, _: &SharedComm, index: usize) -> Result<Option<QContribution>, ScfError> {
        self.load_record(index)
    }
    fn save(&mut self, _: &SharedComm, record: &QContribution) -> Result<(), ScfError> {
        self.save_record(record)
    }
}
impl QStore<Serial> for RpaStore {
    fn load(&mut self, _: &Serial, index: usize) -> Result<Option<QContribution>, ScfError> {
        self.load_record(index)
    }
    fn save(&mut self, _: &Serial, record: &QContribution) -> Result<(), ScfError> {
        self.save_record(record)
    }
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
        .response
        .as_ref()
        .ok_or_else(|| err(1, "rpa requires a [response] table"))?;
    let default_rpa = dft_cli::config::Rpa::default();
    let rpa_section = ctx.config.rpa.as_ref().unwrap_or(&default_rpa);
    if !ctx.selection.reference() {
        return Err(err(1, "rpa requires the reference driver"));
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
    let series = rpa_section
        .ecut_chi_ha_series
        .clone()
        .unwrap_or_else(|| vec![section.ecut_chi_ha]);
    let mut options = RpaOptions::new(series.clone(), section.k_mesh);
    options.bands = bands;
    options.solver = solver;
    if let Some(t) = section.eigen_tolerance_ha {
        options.eigen_tolerance = t;
    }
    if let Some(n) = rpa_section.q0_directions {
        options.q0_directions = n;
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
    let symmetry_name = rpa_section.q_symmetry.as_deref().unwrap_or("crystal");
    options.q_symmetry = match symmetry_name {
        "none" => QSymmetry::Full,
        "time-reversal" => QSymmetry::TimeReversal,
        "crystal" => QSymmetry::Crystal(
            dft_geometry::symmetry::detect(&ctx.cell, &ctx.atoms, ctx.config.detection_tolerance())
                .map_err(|e| err(4, e))?,
        ),
        other => {
            return Err(err(
                1,
                format!("rpa.q_symmetry = {other:?}: use crystal, time-reversal or none"),
            ));
        }
    };
    let sphere_name = rpa_section.response_sphere.as_deref().unwrap_or("q+G");
    options.sphere_at_gamma = match sphere_name {
        "q+G" => false,
        "G" => true,
        other => {
            return Err(err(
                1,
                format!("rpa.response_sphere = {other:?}: use \"q+G\" or \"G\""),
            ));
        }
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
        "rpa: k mesh {:?}, ecut_chi series {:?} Ha, bands {}, solver {}, {} imaginary frequencies, q symmetry {}, q0 directions {}",
        section.k_mesh,
        series,
        bands.map_or("auto (N_G per cutoff)".to_owned(), |b| b.to_string()),
        solver.name(),
        options.frequency_points,
        symmetry_name,
        options.q0_directions
    );
    let potentials: Vec<_> = ctx
        .potentials
        .iter()
        .map(|(label, potential)| dft_engine::SpeciesPotential { label, potential })
        .collect();
    let mut store = RpaStore {
        dir: checkpoint_dir.map(Path::to_path_buf),
        fingerprint: sha256(format!("rpa_q_v1|{}|{density_sha}", ctx.input_hash).as_bytes()),
    };
    let result = rpa(
        StartInput {
            cell: &ctx.cell,
            atoms: &ctx.atoms,
            potentials: &potentials,
            options: &ctx.options,
        },
        &scf.density,
        &options,
        &mut store,
    )
    .map_err(|e| err(4, e))?;
    let natoms = result.atoms as f64;
    let per_atom_mev =
        |e: f64| dft_units::from_hartree(e, dft_units::EnergyUnit::MilliElectronVolt) / natoms;

    let mut q_csv = String::from(
        "# baandi rpa (R3): per-q RPA correlation contribution E_c(q) = (1/2pi) int dw Tr[ln(1 - v chi0) + v chi0], Hartree per cell (not q-weighted); weighted_ha = weight * E_c(q); q in fractional reciprocal coordinates; q = 0 uses the direction average of the k.p head/wings limit\nq_index,q1,q2,q3,weight,multiplicity,optical,ecut_chi_ha,response_size,e_c_q_ha,weighted_ha,min_pivot,hermiticity\n",
    );
    for q in &result.q_points {
        for (c, qc) in q.cutoffs.iter().enumerate() {
            writeln!(
                q_csv,
                "{},{:.12},{:.12},{:.12},{:.16e},{},{},{},{},{:.16e},{:.16e},{:.16e},{:.3e}",
                q.index,
                q.fractional[0],
                q.fractional[1],
                q.fractional[2],
                q.weight,
                q.multiplicity,
                q.optical,
                result.cutoffs[c],
                qc.response_size,
                qc.energy,
                q.weight * qc.energy,
                qc.min_pivot,
                qc.hermiticity
            )
            .unwrap();
        }
    }
    let fit = result.extrapolation.as_ref();
    let mut ext_csv = String::from(
        "# baandi rpa (R3): E_c per response cutoff (bands rule in rpa.json); x = ecut_chi_ha^-1.5 (Ha^-1.5); fit E_c = E_inf + A x by least squares (plan 2.5); residual = E_c - fit\necut_chi_ha,x,response_size_q0,bands,e_c_ha_per_cell,e_c_mev_per_atom,fit_ha_per_cell,residual_mev_per_atom\n",
    );
    for c in 0..result.cutoffs.len() {
        let x = result.cutoffs[c].powf(-1.5);
        let (fit_value, residual) = fit.map_or((f64::NAN, f64::NAN), |f| {
            (f.value + f.slope * x, per_atom_mev(f.residuals[c]))
        });
        writeln!(
            ext_csv,
            "{},{:.16e},{},{},{:.16e},{:.10},{:.16e},{:.6e}",
            result.cutoffs[c],
            x,
            result.response_sizes[c],
            result.bands[c],
            result.energies[c],
            per_atom_mev(result.energies[c]),
            fit_value,
            residual
        )
        .unwrap();
    }
    let mut freq_csv = String::from(
        "# baandi rpa (R3): q-weighted integrand sum_q w_q Tr[ln(1 - v chi0(q, i omega)) + v chi0] (dimensionless, per cell) per imaginary frequency; E_c = (1/2pi) sum weight_ha * term\nindex,omega_ha,omega_ev,weight_ha",
    );
    for c in &result.cutoffs {
        write!(freq_csv, ",term_ecut_{c}").unwrap();
    }
    freq_csv.push('\n');
    for (i, (w, g)) in result.frequencies.iter().zip(&result.weights).enumerate() {
        write!(freq_csv, "{i},{w:.16e},{:.16e},{g:.16e}", ev(*w)).unwrap();
        for terms in &result.frequency_terms {
            write!(freq_csv, ",{:.16e}", terms[i]).unwrap();
        }
        freq_csv.push('\n');
    }
    let q0 = result.q_points.iter().find(|q| q.optical);
    let direction_check = q0.map(|q| {
        q.cutoffs
            .iter()
            .zip(&result.cutoffs)
            .map(|(c, e)| {
                json!({
                    "ecut_chi_hartree": e,
                    "e_c_q0_hartree": c.direction_energies.iter().map(|(l, v)| json!({"set": l, "value": v})).collect::<Vec<_>>(),
                    "weighted_spread_vs_product_mev_per_atom": c.direction_energies.iter().map(|(l, v)| {
                        let reference = c.direction_energies.last().map_or(*v, |r| r.1);
                        json!({"set": l, "value": per_atom_mev(q.weight * (v - reference))})
                    }).collect::<Vec<_>>(),
                    "body_only_e_c_q0_hartree": c.body_only_energy,
                })
            })
            .collect::<Vec<_>>()
    });
    let settings = json!({
        "ecut_chi_series_hartree": result.cutoffs,
        "response_size_q0": result.response_sizes,
        "bands_per_cutoff": result.bands,
        "bands_rule": if bands.is_some() { "explicit, same for every cutoff" } else { "auto = N_G(ecut_chi) at q = 0 per cutoff, clamped to the smallest k basis (plan 2.5)" },
        "nscf_bands": result.nscf_bands,
        "band_solver": solver.name(),
        "eigen_tolerance_hartree": options.eigen_tolerance,
        "k_mesh": section.k_mesh,
        "k_sampling": "Gamma-centred full mesh NSCF (no wavefunction rotation)",
        "q_mesh": section.k_mesh,
        "q_symmetry": symmetry_name,
        "q_reduction": result.q_reduction,
        "irreducible_q": result.q_points.len(),
        "response_sphere": if options.sphere_at_gamma { "G: |G|^2/2 <= ecut_chi for every q (v uses |q+G|); q = first-Brillouin-zone image of each mesh q, mean over tied images (ADR 0055 item 2)" } else { "q+G: |q+G|^2/2 <= ecut_chi" },
        "fft_grid": result.grid_dims,
        "q0_treatment": "direction average of the q -> 0 limit (k.p head and wings), Lebedev set",
        "q0_directions": options.q0_directions,
        "frequency_grid": "gauss-legendre mapped w = w0 t / (1 - t)",
        "frequency_points": options.frequency_points,
        "omega0_hartree": result.omega0,
        "omega0_rule": if options.omega0_hartree.is_some() { "explicit" } else { "auto = mean(lowest direct transition, largest transition)" },
        "band_block": options.band_block,
        "checkpoint_dir": checkpoint_dir.map(|p| p.display().to_string()),
        "resumed_q": result.resumed,
    });
    let mut meta = json!({
        "settings": settings,
        "schema_version": 1,
        "kind": "baandi_rpa_r3",
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
        "theory": "RPA correlation energy from the ACFDT on the imaginary axis, E_c = sum_q w_q (1/2pi) int_0^inf dw Tr[ln(1 - v chi0(q, i w)) + v chi0(q, i w)], chi0 by the sum over states (Adler-Wiser, D1) on KS orbitals of the reference SCF density; q = 0 term as the Lebedev direction average of the k.p head/wings limit (velocity with the nonlocal commutator); no xc kernel, no EXX",
        "units": {"energy": "Hartree per cell unless stated; meV per atom where named", "frequency": "Hartree on the imaginary axis", "q": "fractional reciprocal coordinates"},
        "conventions": "v = 4 pi / |q+G|^2; X = v^1/2 chi0 v^1/2 Hermitized; ln det(1 - X) by Cholesky; q weights sum to one; spin-unpolarized chi0 carries the factor 2, collinear chi0 = chi0_up + chi0_down",
        "spectrum": {
            "electrons": result.electrons,
            "occupied_bands_per_channel": result.occupied,
            "valence_maximum_hartree": result.valence_maximum,
            "conduction_minimum_hartree": result.conduction_minimum,
            "indirect_gap_ev": ev(result.conduction_minimum - result.valence_maximum),
            "lowest_direct_gap_ev": ev(result.lowest_direct_gap),
            "largest_transition_hartree": result.largest_transition,
            "max_residual_hartree": result.max_residual,
            "basis_size_min": result.basis_sizes.0,
            "basis_size_max": result.basis_sizes.1,
        },
        "atoms": result.atoms,
        "lattice_columns_bohr": ctx.cell.lattice_columns_bohr(),
        "structure": ctx.atoms.iter().map(|a| json!({"species": a.species, "fractional": a.fractional_position})).collect::<Vec<_>>(),
        "e_c_hartree_per_cell": result.energies,
        "e_c_mev_per_atom": result.energies.iter().map(|e| per_atom_mev(*e)).collect::<Vec<_>>(),
        "extrapolation": fit.map(|f| json!({
            "law": "E_c(ecut_chi) = E_inf + A ecut_chi^-1.5 (Harl and Kresse), least squares over the series",
            "e_c_inf_hartree_per_cell": f.value,
            "e_c_inf_mev_per_atom": per_atom_mev(f.value),
            "slope_hartree_ha1_5": f.slope,
            "residuals_mev_per_atom": f.residuals.iter().map(|r| per_atom_mev(*r)).collect::<Vec<_>>(),
            "max_residual_last3_mev_per_atom": f.max_residual_last3.map(per_atom_mev),
            "linearity_check": match f.max_residual_last3 {
                None => "not available (fewer than three cutoffs)",
                Some(r) if per_atom_mev(r).abs() <= 0.5 => "passed (last three residuals <= 0.5 meV/atom)",
                Some(_) => "failed (add cutoffs, plan 2.5)",
            },
            "points": result.cutoffs.len(),
            "protocol_note": if result.cutoffs.len() < 4 { "plan 2.5 asks for at least four cutoffs spanning a factor >= 2" } else { "four or more cutoffs" },
        })),
        "extrapolation_inverse_ng": result.extrapolation_inverse_ng.as_ref().map(|f| json!({
            "law": "E_c = E_inf + B / N_G(q = 0) (discrete-shell form of the same law; diagnostic, not the headline)",
            "e_c_inf_hartree_per_cell": f.value,
            "e_c_inf_mev_per_atom": per_atom_mev(f.value),
            "slope_hartree_times_ng": f.slope,
            "residuals_mev_per_atom": f.residuals.iter().map(|r| per_atom_mev(*r)).collect::<Vec<_>>(),
            "max_residual_last3_mev_per_atom": f.max_residual_last3.map(per_atom_mev),
        })),
        "q_points": result.q_points.iter().map(|q| json!({
            "index": q.index,
            "fractional": q.fractional,
            "weight": q.weight,
            "multiplicity": q.multiplicity,
            "optical": q.optical,
            "transitions": q.transitions,
            "e_c_q_hartree": q.cutoffs.iter().map(|c| c.energy).collect::<Vec<_>>(),
            "response_size": q.cutoffs.iter().map(|c| c.response_size).collect::<Vec<_>>(),
        })).collect::<Vec<_>>(),
        "q0_direction_check": direction_check,
        "frequencies_hartree": result.frequencies,
        "weights_hartree": result.weights,
        "frequency_terms": result.frequency_terms,
        "diagnostics": {
            "min_cholesky_pivot": result.min_pivot,
            "max_relative_hermiticity_defect": result.max_hermiticity,
            "transitions": result.transitions,
        },
        "limitations": "insulators only (metals reject, D6); collinear spin for integer-moment insulators (functional lda-pw or pbe); full-mesh NSCF (no IBZ wavefunction rotation); q = 0 by the direction average of the q -> 0 limit (|q| dependence in the mini zone not integrated: O(1/N_q) error, measure with the k-mesh series); absolute NC-pseudopotential E_c lacks core-valence correlation: compare differences with literature, absolute values only with an oracle at matched settings",
    });
    ctx.selection.record(&mut meta);
    ctx.write("rpa_q.csv", &q_csv)?;
    ctx.write("rpa_extrapolation.csv", &ext_csv)?;
    ctx.write("rpa_frequency.csv", &freq_csv)?;
    ctx.json("rpa.json", &meta)?;
    ctx.json("results.json", &meta)?;
    let last = result.energies.len() - 1;
    println!(
        "rpa done: E_c(ecut_chi {} Ha) = {:.10} Ha/cell = {:.4} meV/atom{}; {} irreducible q ({} resumed); outputs in {}",
        result.cutoffs[last],
        result.energies[last],
        per_atom_mev(result.energies[last]),
        fit.map_or(String::new(), |f| format!(
            "; extrapolated {:.10} Ha/cell = {:.4} meV/atom",
            f.value,
            per_atom_mev(f.value)
        )),
        result.q_points.len(),
        result.resumed,
        ctx.out.display()
    );
    Ok(())
}
