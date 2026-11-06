//! Finite-temperature output boundary; all energies are explicitly labelled.
use super::*;
use dft_engine::finite_temperature::{FiniteTCertificate, FiniteTProblem};

fn certificate(c: &FiniteTCertificate, opts: &ScfOptions) -> Value {
    json!({"passed":c.passed,"scope":"retained_basis","map_version":1,"gap_status":"not_applicable","charge_density_residual":c.charge_density_residual,"spin_density_residual":c.spin_density_residual,"retained_eigen_residual_ha":c.retained_eigen_residual_ha,"occupation_change_electrons":c.occupation_change_electrons,"electron_residual":c.electron_residual,"moment_residual":c.moment_residual,"free_energy_change_ha":c.free_energy_change_ha,"edge_max_fraction":c.edge_max_fraction,"edge_electrons":c.edge_electrons,"thresholds":{"density":opts.density_tolerance,"eigen_ha":opts.eigen_tolerance,"free_energy_ha":opts.energy_tolerance,"occupation_relative":1e-8,"number_relative":1e-10,"edge_fraction":1e-8,"edge_electrons_relative":1e-10}})
}
#[allow(clippy::too_many_arguments)]
pub(super) fn execute(
    config: &Config,
    prepared: &crate::config::PreparedInput,
    opts: ScfOptions,
    atoms: &[Atom],
    potentials: &[SpeciesPotential<'_>],
    points: &[dft_geometry::kpoints::KPoint],
    out: &Path,
    theta: f64,
    target: f64,
    preflight: crate::config::preflight::Report,
    species: Vec<Value>,
    input_bytes: &[u8],
) -> Result<RunReport, RunError> {
    let cap = config
        .occupations
        .as_ref()
        .and_then(|o| o.band_cap)
        .unwrap_or(32);
    let sector = config.spin_sector(target).map_err(|e| error(1, e))?;
    let gamma_real_fallback_reason;
    let force_result;
    let stress_result;
    let (
        converged,
        iterations,
        densities,
        energy,
        thermo,
        occupation,
        cert,
        eigenvalues,
        residuals,
        spin,
    ) = if let Some(sector) = sector {
        let mut problem = dft_engine::collinear::CollinearProblem::new_fermi_dirac(
            &prepared.cell,
            atoms,
            potentials,
            opts.clone(),
            points,
            sector,
            theta,
            config.finite_constraint(),
            cap,
        )
        .map_err(engine_error)?;
        if config.symmetry.mode == crate::config::SymmetryMode::Auto {
            problem
                .enable_symmetry(
                    &config.symmetry_full_mesh().map_err(|e| error(1, e))?,
                    config.symmetry.tolerance_bohr,
                )
                .map_err(engine_error)?;
        }
        let r = problem.run().map_err(engine_error)?;
        gamma_real_fallback_reason = r.gamma_real_fallback_reason.clone();
        force_result = requested_forces(
            config,
            &prepared.cell,
            atoms,
            potentials,
            &opts,
            dft_engine::forces::ElectronicState::Collinear(&r),
            None,
        )?;
        stress_result = requested_stress(
            config,
            &prepared.cell,
            atoms,
            potentials,
            &opts,
            dft_engine::forces::ElectronicState::Collinear(&r),
            None,
        )?;
        let eigenvalues = r
            .k_points
            .iter()
            .flatten()
            .flat_map(|k| k.solution.eigenvalues.iter().copied())
            .collect::<Vec<_>>();
        let residuals = r
            .k_points
            .iter()
            .flatten()
            .flat_map(|k| k.solution.residuals.iter().copied())
            .collect::<Vec<_>>();
        let dft_physics::spin::SpinDensity::Collinear { channels, .. } = r.density else {
            unreachable!()
        };
        let constraint = match config.finite_constraint() {
            dft_physics::occupations::Constraint::Total => json!({"kind":"none"}),
            dft_physics::occupations::Constraint::FixedMoment(m) => {
                json!({"kind":"fixed_total_moment","electron_difference":m})
            }
        };
        let spin = json!({"mode":"collinear","channel_order":["up","down"],"channel_capacity_electrons":[1,1],"axis_cartesian":r.sector.axis_cartesian,"moment_mu_b":r.moment_mu_b,"constraint":constraint,"total_electrons":r.channel_electrons.iter().sum::<f64>()});
        (
            r.converged,
            r.history.len(),
            channels.to_vec(),
            r.energy,
            r.thermodynamics.unwrap(),
            r.finite_occupations.unwrap(),
            r.finite_certificate,
            eigenvalues,
            residuals,
            spin,
        )
    } else {
        let mut problem = FiniteTProblem::new(
            &prepared.cell,
            atoms,
            potentials,
            opts.clone(),
            points,
            theta,
            cap,
        )
        .map_err(engine_error)?;
        if config.symmetry.mode == crate::config::SymmetryMode::Auto {
            problem
                .enable_symmetry(
                    &config.symmetry_full_mesh().map_err(|e| error(1, e))?,
                    config.symmetry.tolerance_bohr,
                    config.basis.time_reversal(),
                )
                .map_err(engine_error)?;
        }
        let r = problem.run().map_err(engine_error)?;
        gamma_real_fallback_reason = r.gamma_real_fallback_reason.clone();
        force_result = requested_forces(
            config,
            &prepared.cell,
            atoms,
            potentials,
            &opts,
            dft_engine::forces::ElectronicState::FiniteTemperature(&r),
            None,
        )?;
        stress_result = requested_stress(
            config,
            &prepared.cell,
            atoms,
            potentials,
            &opts,
            dft_engine::forces::ElectronicState::FiniteTemperature(&r),
            None,
        )?;
        let eigenvalues = r
            .state
            .k_points
            .iter()
            .flat_map(|k| k.solution.eigenvalues.iter().copied())
            .collect::<Vec<_>>();
        let residuals = r
            .state
            .k_points
            .iter()
            .flat_map(|k| k.solution.residuals.iter().copied())
            .collect::<Vec<_>>();
        (
            r.converged,
            r.iterations,
            vec![r.state.density],
            r.state.energy,
            r.thermodynamics,
            r.occupations,
            r.certificate,
            eigenvalues,
            residuals,
            json!({"mode":"unpolarized","channel_order":["total"],"channel_capacity_electrons":[2],"axis_cartesian":null,"constraint":{"kind":"none"},"total_electrons":target,"moment_mu_b":0}),
        )
    };
    let e = occupation.evaluation.unwrap();
    let roots: Vec<_> = e.roots.into_iter().flatten().collect();
    let mu: Vec<_> = roots.iter().map(|r| r.mu_ha).collect();
    let block = json!({"policy_version":1,"policy":"fermi_dirac","temperature_energy_ha":theta,"occupations_electrons":{"axes":["spin","k","band"],"shape":[densities.len(),points.len(),cap],"ordering":"C_band_fastest","values":occupation.occupations()},"k_weights":points.iter().map(|k|k.weight).collect::<Vec<_>>(),"electron_target":target,"electron_count":e.electron_count,"electron_residual":e.electron_residual,"moment_residual":e.moment_residual,"chemical_potential":{"kind":if roots.len()==1 {"common"} else {"per_spin"},"values_ha":mu},"chemical_potential_reason":null,"entropy_kind":"physical_fermi_dirac","entropy_over_kb":e.entropy_over_kb,"entropy_term_ha":e.entropy_term_ha,"root_iterations":roots.iter().map(|r|r.iterations).sum::<usize>(),"root_brackets_ha":roots.iter().map(|r|r.bracket_ha).collect::<Vec<_>>(),"edge_max_fraction":e.edge_max_fraction,"edge_electrons":e.edge_electrons});
    let thermodynamics = json!({"internal_energy_ha":thermo.internal_energy_ha,"entropy_term_ha":thermo.entropy_term_ha,"free_energy_ha":thermo.free_energy_ha,"entropy_kind":"physical_fermi_dirac","entropy_over_kb":e.entropy_over_kb,"temperature_energy_ha":theta,"variational_quantity":"free_energy","grand_potential":{"status":"unavailable","reason":"fixed_electron_number"}});
    let mut total = densities[0].clone();
    if densities.len() == 2 {
        total += &densities[1];
    }
    let cube_atoms: Option<Vec<_>> = atoms
        .iter()
        .map(|a| {
            let p = potentials
                .iter()
                .find(|p| p.label == a.species)
                .unwrap()
                .potential;
            volumetric_out::atomic_number(p.element())
                .ok()
                .map(|atomic_number| volumetric_out::OutputAtom {
                    atomic_number,
                    valence_charge: p.z_valence(),
                    position_bohr: prepared
                        .cell
                        .fractional_to_cartesian_bohr(a.fractional_position),
                })
        })
        .collect();
    volumetric_out::write_cube(
        std::io::BufWriter::new(fs::File::create(out.join("density-total.cube"))?),
        &prepared.cell,
        cube_atoms.as_deref().unwrap_or(&[]),
        &total,
        "Finite-T total density; electrons/Bohr^3; synthetic atomic numbers omitted",
    )
    .map_err(|e| error(5, e))?;
    if densities.len() == 2 {
        let moment = &densities[0] - &densities[1];
        for (name, data) in [
            ("density-up", &densities[0]),
            ("density-down", &densities[1]),
            ("magnetization", &moment),
        ] {
            volumetric_out::write_cube(
                std::io::BufWriter::new(fs::File::create(out.join(format!("{name}.cube")))?),
                &prepared.cell,
                cube_atoms.as_deref().unwrap_or(&[]),
                data,
                &format!("Finite-T {name}; electrons/Bohr^3"),
            )
            .map_err(|e| error(5, e))?;
        }
    }
    let bytes: Vec<_> = total.iter().flat_map(|v| v.to_le_bytes()).collect();
    fs::write(out.join("density.bin"), &bytes)?;
    write_json(
        &out.join("density.json"),
        &json!({"fft_grid":prepared.grid.dims(),"lattice_columns_bohr":prepared.cell.lattice_columns_bohr(),"units":"electrons/Bohr^3","ordering":"C_z_fastest","dtype":"f64","endianness":"little","sha256":sha256(&bytes)}),
    )?;
    let physical_spin = json!({"mode":spin["mode"],"axis_cartesian":spin["axis_cartesian"],"constraint":spin["constraint"],"channel_capacity_electrons":spin["channel_capacity_electrons"]});
    let policy_record = json!({"policy":"fermi_dirac","version":1,"theta_ha":theta,"target":target,"spin":physical_spin,"weights":points.iter().map(|p|p.weight).collect::<Vec<_>>(),"bands":cap,"number_relative_tolerance":1e-13,"max_root_iterations":256});
    let policy_digest = super::spin_run::digest(&policy_record);
    let cert = cert.as_ref().map(|c| {
        let mut value = certificate(c, &opts);
        value["policy_digest"] = json!(policy_digest);
        value
    });
    let mut k_basis = Vec::with_capacity(points.len());
    for point in points {
        let basis = dft_basis::PlaneWaveBasis::at_k(
            &prepared.cell,
            dft_units::Hartree(opts.cutoff_hartree),
            point.fractional,
        )
        .map_err(|e| error(4, e))?;
        let miller: Vec<_> = basis.vectors().iter().map(|g| g.miller).collect();
        k_basis.push(json!({"fractional":point.fractional,"weight":point.weight,"g_count":miller.len(),"basis_order_digest":super::spin_run::digest(&json!(miller))}));
    }
    let basis_record = json!({"canonicalization_version":"baandi-mf03-record-v1","cell_columns_bohr":prepared.cell.lattice_columns_bohr(),"grid":prepared.grid.dims(),"origin_fractional":[0,0,0],"periodic":[true,true,true],"cutoff_ha":opts.cutoff_hartree,"density_cutoff_ha":4.0*opts.cutoff_hartree,"ordering":"C_z_fastest","normalization":"forward_1_over_N_inverse_unscaled","derivative":"cartesian_zero_axis_nyquist_v1","band_count":cap,"channel_capacity":spin["channel_capacity_electrons"],"k_points":k_basis,"grid_id":"density","grid_role":"density","wavefunction_grid_alias":"same_dimensions_density_grid","quadrature_weight_bohr3":prepared.cell.volume_bohr3()/total.len() as f64});
    let basis_digest = super::spin_run::digest(&basis_record);
    let mut model_record = json!({"canonicalization_version":"baandi-mf03-record-v1","basis_digest":basis_digest,"functional":opts.functional.name(),"xc_version":"xcx-0.3.1_or_internal_pz","hubbard":"none","core_density":if energy.nlcc.is_some(){"immutable_spin_neutral_nlcc"}else{"absent"},"relativity":"scalar_or_nonrelativistic_no_SOC","allow_functional_override":opts.allow_functional_override,"spin":physical_spin,"occupation_policy":policy_record,"symmetry":"resolved_mesh","environment":"periodic_neutral_valence_G0_zero","atoms":atoms.iter().map(|a|json!({"species":a.species,"fractional":a.fractional_position})).collect::<Vec<_>>(),"datasets":species.iter().map(|p|json!({"name":p["name"],"sha256":p["sha256"],"functional":p["functional"],"valence":p["z_valence"]})).collect::<Vec<_>>()});
    if let Some(symmetry) = config.symmetry_json().map_err(|e| error(1, e))? {
        model_record["symmetry"] = symmetry;
    }
    let model_digest = super::spin_run::digest(&model_record);
    let channel_bytes: Vec<_> = densities
        .iter()
        .flat_map(|n| {
            n.iter()
                .flat_map(|v| (if *v == 0.0 { 0.0 } else { *v }).to_le_bytes())
        })
        .collect();
    let state_digest = super::spin_run::digest(
        &json!({"model_digest":model_digest,"density_sha256":sha256(&channel_bytes),"band_data":{"eigenvalues_ha":eigenvalues,"residuals_ha":residuals,"occupations":block},"thermodynamics":thermodynamics,"converged":converged,"certificate":cert}),
    );
    let forces = json!({"status":"unavailable","reason":"not_requested","units":"Ha/Bohr","frame":"cartesian_input_atom_order"});
    let stress = json!({"status":"unavailable","reason":"not_requested","units":"Ha/Bohr^3","frame":"cartesian_tensor","sign_convention":"tensile_positive"});
    let extension = json!({"proposal_status":"proposed","extension_version":2,"model_digest":model_digest,"basis_digest":basis_digest,"state_digest":state_digest,"spin":spin,"occupations":block,"thermodynamics":thermodynamics,"forces":forces,"stress":stress});
    let mut result = json!({"schema_version":2,"output_profile":"fermi_dirac_v1","spin":spin,"occupations":block,"thermodynamics":thermodynamics,"energies":energy_json(&energy),"energy_ledger":[{"term":"E","energy_ha":thermo.internal_energy_ha},{"term":"-TS","energy_ha":-thermo.entropy_term_ha},{"term":"F","energy_ha":thermo.free_energy_ha}],"band_gap":{"status":"not_applicable","reason":"metallic/Fermi level","chemical_potential_ha":mu},"band_data":{"shape":[densities.len(),points.len(),cap],"eigenvalues_ha":eigenvalues,"residuals_ha":residuals},"band_policy":{"kind":"prepare_explicit_cap","cap":cap,"retained_bands":cap,"scope":"retained_basis","omitted_tail_bound":null},"convergence":{"converged":converged,"iterations":iterations,"certificate":cert},"forces":{"status":"unavailable","reason":"not_requested"},"stress":{"status":"unavailable","reason":"not_requested"},"density":{"channel_values":densities.iter().map(|d|d.iter().copied().collect::<Vec<_>>()).collect::<Vec<_>>()},"limitations":["synthetic validation only; no real-metal QE comparison","finite edge occupancy is not an infinite omitted-state bound","finite-T restart and runtime control unavailable"]});
    if opts.gamma_real != dft_engine::GammaReal::Off {
        result["gamma_real_used"] = json!(false);
        result["gamma_real_fallback_reason"] = json!(gamma_real_fallback_reason);
    }
    if let Some(symmetry) = config.symmetry_json().map_err(|e| error(1, e))? {
        result["symmetry"] = symmetry;
    }
    result["state_contracts"] = extension;
    result["identity"] = json!({"model":model_record,"basis":basis_record});
    result["preflight"] = json!(preflight);
    result["species"] = json!(species);
    result["scf"] = options_json(&opts);
    result["xc"] = json!(opts.functional.name());
    result["provenance"] = json!({"initial_moment_seed":config.spin.as_ref().and_then(|s|s.initial_moment_seed),"seed_role":"initial_density_only_not_constraint"});
    result["input_sha256"] = json!(sha256(input_bytes));
    result["tool"] = json!({"name":"dft-cli","version":env!("CARGO_PKG_VERSION"),"git_commit":env!("DFT_GIT_COMMIT")});
    result["energy_audit"] = json!({"free_energy_band_ha":energy.band_form-e.entropy_term_ha,"free_energy_direct_ha":thermo.free_energy_ha});
    result["forces"] = forces;
    result["stress"] = stress;
    if let Some(force) = force_result {
        result["forces"] = force_json(&force);
        result["state_contracts"]["forces"] = force_json(&force);
        fs::write(out.join("forces.txt"), force_table(&force)).map_err(|e| error(1, e))?;
    }
    if let Some(stress) = stress_result {
        result["stress"] = stress_json(&stress);
        result["state_contracts"]["stress"] = stress_json(&stress);
        fs::write(out.join("stress.txt"), stress_table(&stress)).map_err(|e| error(1, e))?;
    }
    write_json(&out.join("results.json"), &result)?;
    Ok(RunReport {
        converged,
        console: format!(
            "Finite-T metallic/Fermi level mu={mu:?} Ha\nE={:.16e} Ha; -TS={:.16e} Ha; F={:.16e} Ha\niterations={iterations}; converged={converged}\n",
            thermo.internal_energy_ha, -thermo.entropy_term_ha, thermo.free_energy_ha
        ),
    })
}

#[cfg(test)]
use crate::test_fixtures as fixtures;
#[test]
fn finite_cli_and_fixed_rejection() {
    let dir = std::env::temp_dir().join(format!(
        "p17b-cli-{}-{}",
        std::process::id(),
        std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos()
    ));
    fs::create_dir_all(&dir).unwrap();
    fs::write(
        dir.join("X.upf"),
        fixtures::gaussian_with(1.0, 1.0, false, 401).replace("element=\"Xx\"", "element=\"Si\""),
    )
    .unwrap();
    let text = r#"
[system]
name="synthetic FD"
lattice_unit="bohr"
lattice_vectors=[[5.0,0.0,0.0],[0.0,5.0,0.0],[0.0,0.0,5.0]]
atoms=[{species="Si",fractional=[0.0,0.0,0.0]}]
[pseudopotentials]
Si="X.upf"
[basis]
ecutwfc_ry=5.0
kpoints=[[0.0,0.0,0.0]]
[xc]
functional="lda-pz"
[occupations]
mode="fermi_dirac"
smearing_ha=0.01
band_cap=12
"#;
    let input = dir.join("input.toml");
    fs::write(&input, text).unwrap();
    let report = super::execute(&input, Some(&dir.join("out"))).unwrap();
    assert!(report.converged);
    let result: Value =
        serde_json::from_slice(&fs::read(dir.join("out/results.json")).unwrap()).unwrap();
    assert_eq!(result["occupations"]["policy"], "fermi_dirac");
    assert_eq!(result["convergence"]["certificate"]["passed"], true);
    assert_eq!(result["forces"]["status"], "unavailable");
    let extension = &result["state_contracts"];
    let schema: Value = serde_json::from_str(include_str!(
        "../../../docs/schemas/results-state.proposal.json"
    ))
    .unwrap();
    let required = schema["required"].as_array().unwrap();
    assert_eq!(extension.as_object().unwrap().len(), required.len());
    for key in required {
        assert!(extension.get(key.as_str().unwrap()).is_some());
    }
    for key in ["model_digest", "basis_digest", "state_digest"] {
        assert_eq!(extension[key].as_str().unwrap().len(), 64);
    }
    if let Ok(path) = std::env::var("DFT_OCCUPATIONS_EVIDENCE_DIR") {
        fs::create_dir_all(&path).unwrap();
        fs::write(
            Path::new(&path).join("cli-unpolarized.json"),
            serde_json::to_vec_pretty(&result).unwrap(),
        )
        .unwrap();
    }
    let e = result["thermodynamics"]["internal_energy_ha"]
        .as_f64()
        .unwrap();
    let ts = result["thermodynamics"]["entropy_term_ha"]
        .as_f64()
        .unwrap();
    let f = result["thermodynamics"]["free_energy_ha"].as_f64().unwrap();
    assert!((f - e + ts).abs() < 1e-12);
    for bad in [
        text.replace("smearing_ha=0.01", "smearing_ha=0.01\nsmearing_ev=0.1"),
        text.replace("fermi_dirac", "gaussian"),
        text.replace("smearing_ha=0.01", "smearing_ha=0.0"),
    ] {
        assert!(Config::parse(&bad).is_err());
    }
    // Exercise the emitted v2 available envelope through the independent checker.
    fs::write(&input, format!("{text}\n[output]\nforces=true\n")).unwrap();
    let force_destination = dir.join("forces");
    assert!(
        super::execute(&input, Some(&force_destination))
            .unwrap()
            .converged
    );
    let checker =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("tests/validate_force_schema.py");
    assert!(
        std::process::Command::new("python3")
            .arg("-I")
            .arg(&checker)
            .arg(force_destination.join("results.json"))
            .status()
            .unwrap()
            .success()
    );
    let fixed = text.replace(
        "mode=\"fermi_dirac\"\nsmearing_ha=0.01\nband_cap=12",
        "mode=\"fixed\"",
    );
    fs::write(&input, fixed).unwrap();
    assert!(super::execute(&input, Some(&dir.join("fixed"))).is_err());
    // Exercise both CLI spin dispatches and their explicit schema shapes.
    fs::write(
        dir.join("X.upf"),
        fixtures::gaussian_with(1.0, 1.0, false, 401)
            .replace("element=\"Xx\"", "element=\"Si\"")
            .replace("SLA PZ", "SLA PW"),
    )
    .unwrap();
    for (name, constraint, roots) in [("shared", "", 1), ("fixed-moment", "total_moment=0.3", 2)] {
        let spin_text = format!(
            "{}\n[spin]\nmode=\"collinear\"\n{constraint}\n[output]\nforces=true\n",
            text.replace("lda-pz", "lda-pw")
        );
        fs::write(&input, spin_text).unwrap();
        let destination = dir.join(name);
        assert!(
            super::execute(&input, Some(&destination))
                .unwrap()
                .converged
        );
        let spin: Value =
            serde_json::from_slice(&fs::read(destination.join("results.json")).unwrap()).unwrap();
        assert_eq!(
            spin["state_contracts"]["occupations"]["occupations_electrons"]["shape"],
            json!([2, 1, 12])
        );
        assert_eq!(
            spin["occupations"]["chemical_potential"]["values_ha"]
                .as_array()
                .unwrap()
                .len(),
            roots
        );
        assert!(spin["convergence"]["certificate"]["spin_density_residual"].is_number());
        assert!(
            std::process::Command::new("python3")
                .arg("-I")
                .arg(&checker)
                .arg(destination.join("results.json"))
                .status()
                .unwrap()
                .success()
        );
        if roots == 2 {
            assert!((spin["spin"]["moment_mu_b"].as_f64().unwrap() - 0.3).abs() < 1e-10);
        }
        if let Ok(path) = std::env::var("DFT_OCCUPATIONS_EVIDENCE_DIR") {
            fs::write(
                Path::new(&path).join(format!("cli-{name}.json")),
                serde_json::to_vec_pretty(&spin).unwrap(),
            )
            .unwrap();
        }
    }
    fs::write(
        &input,
        text.replace(
            "mode=\"fermi_dirac\"\nsmearing_ha=0.01\nband_cap=12",
            "mode=\"fixed\"",
        ),
    )
    .unwrap();
    // A real v1 fixed checkpoint must not become an FD state by input changes.
    fs::write(
        dir.join("X.upf"),
        fixtures::gaussian_with(2.0, 1.0, false, 401).replace("element=\"Xx\"", "element=\"Si\""),
    )
    .unwrap();
    let checkpoint = dir.join("fixed.cp");
    crate::checkpoint_io::create(&input, &checkpoint, 1, crate::checkpoint_io::Kind::Exact)
        .unwrap();
    fs::write(&input, text).unwrap();
    let error =
        crate::checkpoint_io::restart(&input, &checkpoint, &dir.join("bad.cp"), false).unwrap_err();
    assert!(error.to_string().contains("checkpoint.occupation_policy"));
    assert!(!dir.join("bad.cp").exists());
    fs::remove_dir_all(dir).unwrap();
}
