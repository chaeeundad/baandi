//! Collinear result boundary. No scalar legacy orbital view is fabricated.
use super::*;
use crate::config::PreparedInput;
use dft_engine::collinear::{CollinearProblem, SpinSector};
use dft_geometry::kpoints::KPoint;
use dft_physics::spin::SpinDensity;

/// Versioned canonical binary encoding of finite JSON-shaped records for this
/// result profile. Length-prefixed UTF-8, sorted object keys, little-endian
/// integers/f64 and normalized signed zero; never hash decimal float text.
pub(super) fn digest(value: &Value) -> String {
    fn put(v: &Value, b: &mut Vec<u8>) {
        match v {
            Value::Null => b.push(0),
            Value::Bool(v) => {
                b.push(1);
                b.push(u8::from(*v));
            }
            Value::Number(v) => {
                if let Some(n) = v.as_u64() {
                    b.push(2);
                    b.extend_from_slice(&n.to_le_bytes());
                } else if let Some(n) = v.as_i64() {
                    b.push(6);
                    b.extend_from_slice(&n.to_le_bytes());
                } else {
                    b.push(7);
                    let f = v.as_f64().expect("finite model number");
                    b.extend_from_slice(&(if f == 0.0 { 0.0 } else { f }).to_le_bytes());
                }
            }
            Value::String(v) => {
                b.push(3);
                b.extend_from_slice(&(v.len() as u64).to_le_bytes());
                b.extend_from_slice(v.as_bytes());
            }
            Value::Array(v) => {
                b.push(4);
                b.extend_from_slice(&(v.len() as u64).to_le_bytes());
                for x in v {
                    put(x, b);
                }
            }
            Value::Object(v) => {
                b.push(5);
                b.extend_from_slice(&(v.len() as u64).to_le_bytes());
                let mut keys: Vec<_> = v.keys().collect();
                keys.sort();
                for k in keys {
                    put(&Value::String(k.clone()), b);
                    put(&v[k], b);
                }
            }
        }
    }
    let mut bytes = b"baandi-mf03-record-v1".to_vec();
    put(value, &mut bytes);
    sha256(&bytes)
}
#[allow(clippy::too_many_arguments)]
pub(super) fn execute(
    prepared: &PreparedInput,
    opts: ScfOptions,
    atoms: &[Atom],
    potentials: &[SpeciesPotential<'_>],
    points: &[KPoint],
    sector: SpinSector,
    out: &Path,
    preflight: crate::config::preflight::Report,
    species: Vec<Value>,
    input_bytes: &[u8],
) -> Result<RunReport, RunError> {
    let r = CollinearProblem::new(
        &prepared.cell,
        atoms,
        potentials,
        opts.clone(),
        points,
        sector,
    )
    .map_err(engine_error)?
    .run()
    .map_err(engine_error)?;
    let SpinDensity::Collinear { channels, .. } = &r.density else {
        unreachable!()
    };
    let total = &channels[0] + &channels[1];
    let moment = &channels[0] - &channels[1];
    let n = (r.sector.n_up + r.sector.n_down) as f64;
    let nb = r.occupations[0].len();
    let nk = points.len();
    let spin = json!({"mode":"collinear","channel_order":["up","down"],"channel_capacity_electrons":[1,1],"axis_cartesian":r.sector.axis_cartesian,"constraint":{"kind":"fixed_total_moment","electron_difference":r.sector.moment()},"total_electrons":r.channel_electrons.iter().sum::<f64>(),"moment_mu_b":r.moment_mu_b});
    let values: Vec<f64> = (0..2)
        .flat_map(|s| {
            (0..nk).flat_map({
                let occupations = &r.occupations;
                move |_| occupations[s].iter().copied()
            })
        })
        .collect();
    let occupations = json!({"policy_version":1,"policy":"fixed_insulator","occupations_electrons":{"axes":["spin","k","band"],"shape":[2,nk,nb],"ordering":"C_band_fastest","values":values},"k_weights":points.iter().map(|k|k.weight).collect::<Vec<_>>(),"electron_target":n,"electron_count":r.channel_electrons.iter().sum::<f64>(),"electron_residual":r.channel_electrons.iter().sum::<f64>()-n,"moment_residual":r.moment_mu_b-r.sector.moment(),"temperature_energy_ha":0,"channel_electrons":[r.sector.n_up,r.sector.n_down],"gap_threshold_ha":1e-7_f64.max(2.0*opts.eigen_tolerance),"chemical_potential":null,"chemical_potential_reason":"fixed_insulator","entropy_kind":"zero","entropy_over_kb":0,"entropy_term_ha":0});
    let thermodynamics = json!({"internal_energy_ha":r.energy.total_hartree(),"free_energy_ha":r.energy.total_hartree(),"entropy_term_ha":0,"entropy_over_kb":0,"entropy_kind":"zero","temperature_energy_ha":0,"variational_quantity":"internal_energy","grand_potential":{"status":"unavailable","reason":"fixed_electron_number"}});
    let eigenvalues: Vec<f64> = r
        .k_points
        .iter()
        .flatten()
        .flat_map(|k| k.solution.eigenvalues.iter().copied())
        .collect();
    let residuals: Vec<f64> = r
        .k_points
        .iter()
        .flatten()
        .flat_map(|k| k.solution.residuals.iter().copied())
        .collect();
    let mut k_basis_records = Vec::with_capacity(nk);
    for k in &r.k_points[0] {
        let basis = dft_basis::PlaneWaveBasis::at_k(
            &prepared.cell,
            dft_units::Hartree(opts.cutoff_hartree),
            k.point.fractional,
        )
        .map_err(|e| error(4, e))?;
        let miller: Vec<_> = basis.vectors().iter().map(|g| g.miller).collect();
        k_basis_records.push(json!({"fractional":k.point.fractional,"weight":k.point.weight,"g_count":miller.len(),"basis_order_digest":digest(&json!(miller))}));
    }
    let basis_record = json!({"canonicalization_version":"baandi-mf03-record-v1","cell_columns_bohr":prepared.cell.lattice_columns_bohr(),"grid":prepared.grid.dims(),"origin_fractional":[0,0,0],"cutoff_ha":opts.cutoff_hartree,"density_cutoff_ha":4.0*opts.cutoff_hartree,"ordering":"C_z_fastest","normalization":"forward_1_over_N_inverse_unscaled","derivative":"cartesian_zero_axis_nyquist_v1","channel_capacity":[1,1],"band_count":nb,"k_points":k_basis_records,"grid_id":"density","grid_role":"density","wavefunction_grid_alias":"same_dimensions_density_grid","periodic":[true,true,true],"quadrature_weight_bohr3":prepared.cell.volume_bohr3()/total.len() as f64});
    let basis_digest = digest(&basis_record);
    let model_record = json!({"canonicalization_version":"baandi-mf03-record-v1","basis_digest":basis_digest,"functional":opts.functional.name(),"xc_version":"xcx-0.3.1","spin_mode":"collinear","hubbard":"none","core_density":if r.energy.nlcc.is_some() {"immutable_spin_neutral_nlcc"} else {"absent"},"relativity":"scalar_or_nonrelativistic_no_SOC","allow_functional_override":opts.allow_functional_override,"axis_cartesian":r.sector.axis_cartesian,"channel_electrons":[r.sector.n_up,r.sector.n_down],"occupation_policy":"fixed_insulator_v1","gap_threshold_ha":1e-7_f64.max(2.0*opts.eigen_tolerance),"symmetry":"full_mesh","environment":"periodic_neutral_valence_G0_zero","atoms":atoms.iter().map(|a|json!({"species":a.species,"fractional":a.fractional_position})).collect::<Vec<_>>(),"datasets":species.iter().map(|p|json!({"name":p["name"],"sha256":p["sha256"],"functional":p["functional"],"valence":p["z_valence"]})).collect::<Vec<_>>()});
    let model_digest = digest(&model_record);
    let density_bits: Vec<u8> = channels
        .iter()
        .flat_map(|n| {
            n.iter()
                .flat_map(|v| (if *v == 0.0 { 0.0 } else { *v }).to_le_bytes())
        })
        .collect();
    let certificate=r.certificate.as_ref().map(|c|json!({"passed":c.passed,"scope":"fixed_moment_sector_retained_basis","map_version":1,"charge_density_residual":c.charge_density_residual,"spin_density_residual":c.spin_density_residual,"energy_change_ha":c.energy_change_ha,"occupied_eigenvalue_change_ha":c.occupied_eigenvalue_change_ha,"strict_eigen_residual_ha":c.strict_eigen_residual_ha,"channel_electrons":c.channel_electrons,"gap_ha":c.gap_ha,"density_tolerance":c.density_tolerance,"energy_tolerance_ha":c.energy_tolerance_ha,"eigen_tolerance_ha":c.eigen_tolerance_ha,"occupied_eigenvalue_change_tolerance_ha":1e-7}));
    let state_digest = digest(
        &json!({"model_digest":model_digest,"density_sha256":sha256(&density_bits),"eigenvalues_ha":eigenvalues,"residuals_ha":residuals,"occupations":values,"thermodynamics":thermodynamics,"converged":r.converged,"certificate":certificate}),
    );
    let extension = json!({"proposal_status":"proposed","extension_version":1,"model_digest":model_digest,"basis_digest":basis_digest,"state_digest":state_digest,"spin":spin,"occupations":occupations,"thermodynamics":thermodynamics,"forces":{"status":"unavailable","reason":"not_implemented","units":"Ha/Bohr","frame":"cartesian_input_atom_order"},"stress":{"status":"unavailable","reason":"not_implemented","units":"Ha/Bohr^3","frame":"cartesian_tensor","sign_convention":"tensile_positive"}});

    let history:Vec<_>=r.history.iter().map(|h|json!({"iteration":h.iteration,"energy_ha":h.energy_ha,"charge_density_residual":h.charge_density_residual,"spin_density_residual":h.spin_density_residual,"mixing_fraction":h.mixing_fraction})).collect();
    // Real atoms use the established independently tested cube writer. Synthetic
    // symbols have no atomic number: retain their full structure in JSON and
    // explicitly omit atoms from the grid-only cubes instead of inventing Z.
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
    for (name, data) in [
        ("density-up", &channels[0]),
        ("density-down", &channels[1]),
        ("density-total", &total),
        ("magnetization", &moment),
    ] {
        let file = fs::File::create(out.join(format!("{name}.cube")))?;
        let label = format!(
            "{name}; electrons/Bohr^3; {}",
            if cube_atoms.is_some() {
                "input atoms"
            } else {
                "grid only: synthetic atomic numbers unavailable"
            }
        );
        volumetric_out::write_cube(
            std::io::BufWriter::new(file),
            &prepared.cell,
            cube_atoms.as_deref().unwrap_or(&[]),
            data,
            &label,
        )
        .map_err(|e| error(5, e))?;
    }
    let mut result = json!({"schema_version":2,"output_profile":"collinear_fixed_v1","preflight":preflight,"tool":{"name":"dft-cli","version":env!("CARGO_PKG_VERSION"),"git_commit":env!("DFT_GIT_COMMIT")},"input_sha256":sha256(input_bytes),"species":species,"xc":opts.functional.name(),"scf":options_json(&opts),"spin":spin,"state_contracts":extension,"identity":{"model":model_record,"basis":basis_record},"provenance":{"initial_moment_seed":r.sector.initial_moment_seed,"seed_role":"provenance_only_uniform_fixed_sector_initialization"},"energies":energy_json(&r.energy),"total_energy":{"ha":r.energy.total_hartree(),"ry":2.0*r.energy.total_hartree()},"band_data":{"axes":["spin","k","band"],"shape":[2,nk,nb],"ordering":"C_band_fastest","eigenvalues_ha":eigenvalues,"residuals_ha":residuals,"occupations_electrons":values},"convergence":{"converged":r.converged,"iterations":r.history.len(),"certificate":certificate,"history":history},"density":{"shape":[2,prepared.grid.dims()[0],prepared.grid.dims()[1],prepared.grid.dims()[2]],"ordering":"C_spin_xyz_z_fastest","units":"electrons/Bohr^3","channel_values":channels.iter().map(|n|n.iter().copied().collect::<Vec<_>>()).collect::<Vec<_>>(),"cube_atoms_included":cube_atoms.is_some()},"limitations":["fixed moment sector, not unconstrained magnetic ground state","full meshes only; no spin restart, runtime control or band/DOS postprocessing"]});
    if let Some(report) = &r.spectral_completeness {
        result["spectral_completeness"] = super::spectral_json(report);
    }
    let _ = dft_engine::work_accounting::take();
    write_json(&out.join("results.json"), &result)?;
    // Conventional total-density binary remains a documented additional format.
    let total_bytes: Vec<u8> = total.iter().flat_map(|v| v.to_le_bytes()).collect();
    fs::write(out.join("density.bin"), &total_bytes)?;
    write_json(
        &out.join("density.json"),
        &json!({"schema_version":1,"fft_grid":prepared.grid.dims(),"lattice_columns_bohr":prepared.cell.lattice_columns_bohr(),"volume_bohr3":prepared.cell.volume_bohr3(),"units":"electrons per Bohr^3","ordering":"C: (i*n2+j)*n3+k; k fastest","endianness":"little","dtype":"f64","normalization":{"expected_electrons":n,"integral_electrons":r.channel_electrons.iter().sum::<f64>()},"sha256":sha256(&total_bytes)}),
    )?;
    Ok(RunReport {
        converged: r.converged,
        console: format!(
            "Spectral completeness: {}\nCollinear fixed sector: N_up={}, N_down={}, M={:.12} mu_B\nTotal energy: {:.16} Ha\n{}\n",
            super::spectral_status(r.spectral_completeness.as_ref()),
            r.sector.n_up,
            r.sector.n_down,
            r.moment_mu_b,
            r.energy.total_hartree(),
            if r.converged {
                "converged"
            } else {
                "NOT CONVERGED: iteration budget exhausted"
            }
        ),
    })
}
