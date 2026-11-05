//! Variable-cell enthalpy driver. Coordinates are L*s and L*(D-I), A=D*A0.
use super::*;
type Matrix = [[f64; 3]; 3];
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum Mode {
    #[default]
    Ions,
    Cell,
    IonsCell,
}
#[derive(Clone, Debug, Default, Serialize, Deserialize, PartialEq)]
#[serde(rename_all = "snake_case")]
pub enum Dofs {
    #[default]
    Full,
    Isotropic,
    XyzDiagonal,
}
fn mul(a: Matrix, b: Matrix) -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| (0..3).map(|k| a[i][k] * b[k][j]).sum()))
}
fn transpose(a: Matrix) -> Matrix {
    std::array::from_fn(|i| std::array::from_fn(|j| a[j][i]))
}
fn inverse(a: Matrix) -> Result<Matrix, RunError> {
    let c = Cell::from_matrix_columns_bohr(a).map_err(|e| error(4, e))?;
    Ok(transpose(std::array::from_fn(|i| {
        c.cartesian_bohr_to_fractional(std::array::from_fn(|j| if i == j { 1. } else { 0. }))
    })))
}
fn deformation(q: &[f64], l: f64) -> Matrix {
    let n = q.len() - 9;
    std::array::from_fn(|i| {
        std::array::from_fn(|j| q[n + 3 * i + j] / l + if i == j { 1. } else { 0. })
    })
}
fn set_cell(ctx: &mut Context, a: Matrix) -> Result<(), RunError> {
    ctx.cell = Cell::from_matrix_columns_bohr(a).map_err(|e| error(4, e))?;
    ctx.config.system.lattice_unit = crate::config::LatticeUnit::Bohr;
    ctx.config.system.lattice_vectors = transpose(a);
    Ok(())
}
fn place(ctx: &mut Context, q: &[f64], a0: Matrix, l: f64) -> Result<Vec<f64>, RunError> {
    let d = deformation(q, l);
    // Symmetric positive definite deformation excludes collapsed/inverted cells.
    if d[0][0] <= 0. || d[0][0] * d[1][1] - d[0][1] * d[1][0] <= 0. || determinant(d) <= 0. {
        return Err(error(4, "cell deformation is not positive definite"));
    }
    set_cell(ctx, mul(d, a0))?;
    let frac: Vec<_> = q[..q.len() - 9]
        .chunks(3)
        .map(|s| [s[0] / l, s[1] / l, s[2] / l])
        .collect();
    ctx.position(&frac);
    Ok(frac
        .iter()
        .flat_map(|f| ctx.cell.fractional_to_cartesian_bohr(*f))
        .collect())
}
fn determinant(a: Matrix) -> f64 {
    a[0][0] * (a[1][1] * a[2][2] - a[1][2] * a[2][1])
        - a[0][1] * (a[1][0] * a[2][2] - a[1][2] * a[2][0])
        + a[0][2] * (a[1][0] * a[2][1] - a[1][1] * a[2][0])
}
fn project_tensor(mut t: Matrix, settings: &Settings, rotations: &[Matrix]) -> Matrix {
    t = std::array::from_fn(|i| std::array::from_fn(|j| (t[i][j] + t[j][i]) * 0.5));
    if !rotations.is_empty() {
        let mut sum = [[0.; 3]; 3];
        for r in rotations {
            let v = mul(mul(*r, t), transpose(*r));
            for i in 0..3 {
                for j in 0..3 {
                    sum[i][j] += v[i][j] / rotations.len() as f64;
                }
            }
        }
        t = sum;
    }
    let mean = (0..3).map(|i| t[i][i]).sum::<f64>() / 3.;
    std::array::from_fn(|i| {
        std::array::from_fn(|j| {
            if settings.fixed_axes[i] || settings.fixed_axes[j] {
                0.
            } else {
                match settings.cell_dofs {
                    Dofs::Isotropic => {
                        if i == j {
                            mean
                        } else {
                            0.
                        }
                    }
                    Dofs::XyzDiagonal => {
                        if i == j {
                            t[i][j]
                        } else {
                            0.
                        }
                    }
                    Dofs::Full => t[i][j],
                }
            }
        })
    })
}
fn project(v: &mut [f64], settings: &Settings, rotations: &[Matrix]) {
    let n = v.len() - 9;
    if settings.mode == Mode::Cell {
        v[..n].fill(0.);
    }
    let t = project_tensor(
        std::array::from_fn(|i| std::array::from_fn(|j| v[n + 3 * i + j])),
        settings,
        rotations,
    );
    for i in 0..3 {
        for j in 0..3 {
            v[n + 3 * i + j] = t[i][j];
        }
    }
}
fn force(
    frame: &Frame,
    q: &[f64],
    l: f64,
    settings: &Settings,
    rotations: &[Matrix],
) -> Result<Vec<f64>, RunError> {
    let a = frame.cell_columns_bohr.unwrap();
    let mut f: Vec<f64> = frame
        .forces_ha_bohr
        .iter()
        .flat_map(|f| (0..3).map(move |i| (0..3).map(|j| a[j][i] * f[j] / l).sum::<f64>()))
        .collect();
    let volume = determinant(a).abs();
    let mut sigma = frame.stress_ha_bohr3.unwrap();
    for (i, row) in sigma.iter_mut().enumerate() {
        row[i] += settings.pressure_ha_bohr3;
    }
    let g = mul(sigma, transpose(inverse(deformation(q, l))?));
    f.extend(g.iter().flatten().map(|v| -volume * v / l));
    project(&mut f, settings, rotations);
    Ok(f)
}
fn enthalpy(f: &Frame, p: f64) -> f64 {
    f.energy_ha + p * determinant(f.cell_columns_bohr.unwrap()).abs()
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct State {
    version: u32,
    ionic: Restart,
    q: Vec<f64>,
    initial: Matrix,
    reference: Matrix,
    grid: [usize; 3],
    rebuilt: Option<Frame>,
    accepted_steps: usize,
}
fn persist(path: &Path, s: &State) -> Result<(), RunError> {
    let mut v = serde_json::to_value(s).map_err(|e| error(5, e))?;
    bits(&mut v, true)?;
    let payload = serde_json::to_string(&v).map_err(|e| error(5, e))?;
    atomic(
        path,
        &serde_json::to_vec(&json!({"payload":payload,"sha256":sha256(payload.as_bytes())}))
            .map_err(|e| error(5, e))?,
    )
}
fn restore(path: &Path, id: &str, n: usize) -> Result<State, RunError> {
    if fs::metadata(path)?.len() > 512 * 1024 * 1024 {
        return Err(error(1, "cell restart too large"));
    }
    let v: Value = serde_json::from_slice(&fs::read(path)?).map_err(|e| error(1, e))?;
    let p = v["payload"]
        .as_str()
        .ok_or_else(|| error(1, "missing cell restart payload"))?;
    if v["sha256"].as_str() != Some(&sha256(p.as_bytes())) {
        return Err(error(1, "cell restart checksum mismatch"));
    }
    let mut v: Value = serde_json::from_str(p).map_err(|e| error(1, e))?;
    bits(&mut v, false)?;
    let s: State = serde_json::from_value(v).map_err(|e| error(1, e))?;
    if s.version != 1
        || s.ionic.identity != id
        || s.q.len() != n
        || s.ionic.frames.is_empty()
        || s.ionic.optimizer.h.len() != n * n
        || s.ionic.optimizer.velocity.len() != n
        || s.ionic.frames.iter().any(|f| {
            f.positions.len() != n - 9
                || f.fractional.len() * 3 != n - 9
                || f.forces_ha_bohr.len() * 3 != n - 9
                || f.cell_columns_bohr.is_none()
                || f.stress_ha_bohr3.is_none()
        })
        || s.grid.contains(&0)
        || s.ionic.density.len() != s.grid.iter().product::<usize>()
        || (!s.ionic.optimizer.previous_x.is_empty()
            && (s.ionic.optimizer.previous_x.len() != n
                || s.ionic.optimizer.previous_force.len() != n))
    {
        return Err(error(1, "cell restart identity/version/shape mismatch"));
    }
    Ok(s)
}
// Orthogonal projector onto invariant fractional displacements. Fractional
// crystallographic rotations need not be orthogonal in the Euclidean metric.
fn invariant_vectors(ctx: &Context) -> Result<Vec<Vec<f64>>, RunError> {
    let Some(map) = ctx.config.symmetry_map().map_err(|e| error(1, e))? else {
        return Ok(vec![]);
    };
    let ops: Vec<_> = map.operations.iter().filter(|o| !o.antiunitary).collect();
    let n = ctx.config.system.atoms.len() * 3;
    let mut basis: Vec<Vec<f64>> = Vec::new();
    for col in 0..n {
        let mut v = vec![0.; n];
        for op in &ops {
            for i in 0..3 {
                v[3 * op.atom_permutation[col / 3] + i] +=
                    op.rotation[i][col % 3] as f64 / ops.len() as f64;
            }
        }
        for b in &basis {
            let dot: f64 = v.iter().zip(b).map(|(x, y)| x * y).sum();
            for (x, y) in v.iter_mut().zip(b) {
                *x -= dot * y;
            }
        }
        let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
        if norm > 1e-10 {
            for x in &mut v {
                *x /= norm;
            }
            basis.push(v);
        }
    }
    Ok(basis)
}
fn project_ions(v: &mut [f64], basis: &[Vec<f64>], enabled: bool) {
    if !enabled {
        return;
    }
    let n = v.len() - 9;
    let mut p = vec![0.; n];
    for b in basis {
        let dot: f64 = v[..n].iter().zip(b).map(|(x, y)| x * y).sum();
        for (x, y) in p.iter_mut().zip(b) {
            *x += dot * y;
        }
    }
    v[..n].copy_from_slice(&p);
}
pub(super) fn execute(
    mut ctx: Context,
    settings: Settings,
    id: String,
    out: &Path,
    resume: bool,
) -> Result<(), RunError> {
    if !settings.fixed.is_empty() {
        return Err(error(
            1,
            "variable-cell relaxation does not yet support Cartesian atom constraints",
        ));
    }
    if settings.cell_dofs == Dofs::Isotropic && settings.fixed_axes.iter().any(|x| *x) {
        return Err(error(
            1,
            "isotropic strain cannot be combined with fixed axes",
        ));
    }
    // Resolve the initial electronic solve at the requested strict tolerance.
    // Otherwise a uniform first iterate can have exactly zero density residual
    // while its loose Ritz solve is nonstationary, triggering the growth guard.
    ctx.opts.initial_eigen_tolerance = ctx.opts.eigen_tolerance;
    ctx.opts.residual_eigen_schedule = false;
    let ionic_basis = invariant_vectors(&ctx)?;
    let symmetry_enabled = ctx.config.symmetry.mode == crate::config::SymmetryMode::Auto;
    let a0 = ctx.cell.lattice_columns_bohr();
    let l = ctx.cell.volume_bohr3().cbrt();
    let rotations: Vec<Matrix> = ctx
        .config
        .symmetry_map()
        .map_err(|e| error(1, e))?
        .map_or_else(Vec::new, |m| {
            m.operations
                .into_iter()
                .filter(|o| !o.antiunitary)
                .map(|o| {
                    mul(
                        mul(a0, o.rotation.map(|r| r.map(|v| v as f64))),
                        inverse(a0).unwrap(),
                    )
                })
                .collect()
        });
    // Require commuting projectors; sequential projections then equal their intersection.
    for i in 0..9 {
        let mut t = [[0.; 3]; 3];
        t[i / 3][i % 3] = 1.;
        let projected = project_tensor(t, &settings, &rotations);
        let symmetric = project_tensor(projected, &Settings::default(), &rotations);
        if projected
            .iter()
            .flatten()
            .zip(symmetric.iter().flatten())
            .any(|(a, b)| (a - b).abs() > 1e-10)
        {
            return Err(error(1, "cell constraints break auto symmetry"));
        }
    }
    let grid = ctx.problem()?.grid_dimensions();
    ctx.opts.grid_override = Some(grid);
    ctx.opts.strain_reference = Some(Box::new(ctx.cell.clone()));
    let path = out.join("restart.json");
    let n = ctx.config.system.atoms.len() * 3 + 9;
    let mut s = if resume {
        restore(&path, &id, n)?
    } else {
        let mut q: Vec<_> = ctx
            .atoms()
            .iter()
            .flat_map(|a| a.fractional_position.map(|x| x * l))
            .collect();
        q.extend([0.; 9]);
        let x = place(&mut ctx, &q, a0, l)?;
        let mut initial_failures = Vec::new();
        let initial = ctx.evaluate(x.clone(), 0., None, &settings).or_else(|e| {
            initial_failures
                .push(json!({"stage":"initial_scf","attempt":1,"reason":e.to_string()}));
            ctx.evaluate(x, 0., None, &settings)
        });
        let (frame, density) = match initial {
            Ok(v) => v,
            Err(e) => {
                initial_failures
                    .push(json!({"stage":"initial_scf","attempt":2,"reason":e.to_string()}));
                let report = json!({"schema_version":1,"calculation_type":"relaxation","provenance":ctx.provenance,"relaxation":{"status":"failed","last_good_structure":null,"failures":initial_failures}});
                atomic(
                    &out.join("results.json"),
                    &serde_json::to_vec_pretty(&report).map_err(|e| error(5, e))?,
                )?;
                return Err(e);
            }
        };
        State {
            version: 1,
            initial: a0,
            reference: a0,
            grid,
            q,
            rebuilt: None,
            accepted_steps: 0,
            ionic: Restart {
                version: 1,
                identity: id,
                optimizer: optimizer::Optimizer::new(n, settings.trust_radius_bohr),
                frames: vec![frame],
                density,
                consecutive: 0,
                status: "running".into(),
                failures: initial_failures,
            },
        }
    };
    if s.initial != a0 {
        return Err(error(1, "cell restart initial cell mismatch"));
    }
    ctx.opts.strain_reference = Some(Box::new(
        Cell::from_matrix_columns_bohr(s.reference).map_err(|e| error(1, e))?,
    ));
    ctx.opts.grid_override = Some(s.grid);
    place(&mut ctx, &s.q, a0, l)?;
    persist(&path, &s)?;
    emit(&ctx, &settings, &s, out)?;
    while s.ionic.status == "running" && s.accepted_steps < settings.max_steps {
        let last = s.ionic.frames.last().unwrap().clone();
        let source = ctx.problem()?;
        let mut f = force(&last, &s.q, l, &settings, &rotations)?;
        project_ions(&mut f, &ionic_basis, symmetry_enabled);
        let mut opt = s.ionic.optimizer.clone();
        let mut step = opt.step(
            &s.q,
            &f,
            settings.method == Method::Fire,
            settings.max_displacement_bohr,
        );
        project(&mut step, &settings, &rotations);
        project_ions(&mut step, &ionic_basis, symmetry_enabled);
        // Cap actual Cartesian atom movement, including affine cell motion.
        let candidate: Vec<_> = s.q.iter().zip(&step).map(|(x, d)| x + d).collect();
        if let Ok(x) = place(&mut ctx, &candidate, a0, l) {
            let maximum = x
                .chunks(3)
                .zip(last.positions.chunks(3))
                .map(|(a, b)| {
                    a.iter()
                        .zip(b)
                        .map(|(x, y)| (x - y).powi(2))
                        .sum::<f64>()
                        .sqrt()
                })
                .fold(0., f64::max);
            if maximum > settings.max_displacement_bohr {
                for d in &mut step {
                    *d *= settings.max_displacement_bohr / maximum;
                }
            }
        }
        place(&mut ctx, &s.q, a0, l)?;
        let mut accepted = None;
        let mut failed_scf_attempts = 0;
        // Backtracking enforces enthalpy monotonicity on each smooth basis segment.
        for attempt in 0..16 {
            let scale = 0.5_f64.powi(attempt);
            let q: Vec<_> = s.q.iter().zip(&step).map(|(x, d)| x + scale * d).collect();
            let trial = place(&mut ctx, &q, a0, l).and_then(|x| {
                let length = x
                    .chunks(3)
                    .zip(last.positions.chunks(3))
                    .map(|(a, b)| {
                        a.iter()
                            .zip(b)
                            .map(|(x, y)| (x - y).powi(2))
                            .sum::<f64>()
                            .sqrt()
                    })
                    .fold(0., f64::max);
                if length > settings.max_displacement_bohr * (1. + 1e-12) {
                    return Err(error(4, "candidate exceeds Cartesian displacement cap"));
                }
                ctx.evaluate(
                    x,
                    length,
                    (attempt == 0 && settings.warm_start)
                        .then_some((&source, s.ionic.density.as_slice())),
                    &settings,
                )
            });
            match trial {
                Ok((frame,density)) if enthalpy(&frame,settings.pressure_ha_bohr3)<=enthalpy(&last,settings.pressure_ha_bohr3)+settings.energy_tol_ha=>{accepted=Some((frame,density,q,attempt));break;}
                Ok((frame,_))=>s.ionic.failures.push(json!({"step":s.ionic.frames.len(),"attempt":attempt+1,"reason":"enthalpy increase","enthalpy_ha":enthalpy(&frame,settings.pressure_ha_bohr3)})),
                Err(e) => {
                    s.ionic.failures.push(json!({"step":s.ionic.frames.len(),"attempt":attempt+1,"reason":e.to_string(),"coordinates":q}));
                    failed_scf_attempts+=1;
                    if failed_scf_attempts>=2 {break;}
                }
            }
        }
        if let Some((frame, density, q, attempt)) = accepted {
            let delta = enthalpy(&frame, settings.pressure_ha_bohr3)
                - enthalpy(&last, settings.pressure_ha_bohr3);
            let gf = force(&frame, &q, l, &settings, &rotations)?;
            let stress = gf[n - 9..]
                .iter()
                .map(|v| v.abs() * l / determinant(frame.cell_columns_bohr.unwrap()).abs())
                .fold(0., f64::max);
            let ions_ok =
                settings.mode == Mode::Cell || frame.max_force_ha_bohr < settings.force_tol_ha_bohr;
            s.ionic.consecutive = if ions_ok
                && stress < settings.stress_tol_ha_bohr3
                && delta.abs() < settings.energy_tol_ha
            {
                s.ionic.consecutive + 1
            } else {
                0
            };
            if attempt > 0 {
                opt.radius *= 0.5;
                opt.velocity.fill(0.);
            } else {
                opt.radius = (opt.radius * 1.1).min(settings.trust_radius_bohr);
            }
            s.accepted_steps += 1;
            s.q = q;
            s.ionic.optimizer = opt;
            s.ionic.density = density;
            s.ionic.frames.push(frame);
            if s.ionic.consecutive >= 2 {
                s.ionic.status = "converged".into();
            }
            if settings.rebuild_every > 0
                && s.accepted_steps.is_multiple_of(settings.rebuild_every)
                && s.ionic.status == "running"
            {
                let old_options = ctx.opts.clone();
                let rebuild = (|| {
                    ctx.opts.strain_reference = None;
                    ctx.opts.grid_override = options(&ctx.config)?.grid_override;
                    let grid = ctx.problem()?.grid_dimensions();
                    ctx.opts.grid_override = Some(grid);
                    ctx.opts.strain_reference = Some(Box::new(ctx.cell.clone()));
                    let x = place(&mut ctx, &s.q, a0, l)?;
                    ctx.evaluate(x, 0., None, &settings)
                })();
                match rebuild {
                    Ok((mut frame, density)) => {
                        s.reference = ctx.cell.lattice_columns_bohr();
                        s.grid = ctx.opts.grid_override.unwrap();
                        frame.provenance = "basis_rebuild".into();
                        s.ionic.frames.push(frame);
                        s.ionic.density = density;
                        s.ionic.consecutive = 0;
                        s.ionic.optimizer =
                            optimizer::Optimizer::new(n, settings.trust_radius_bohr);
                    }
                    Err(e) => {
                        ctx.opts = old_options;
                        s.ionic.status = "failed".into();
                        s.ionic
                            .failures
                            .push(json!({"stage":"periodic_rebuild","reason":e.to_string()}));
                    }
                }
            }
        } else {
            s.ionic.status = "failed".into();
            place(&mut ctx, &s.q, a0, l)?;
        }
        if s.ionic.status == "running" && s.accepted_steps >= settings.max_steps {
            s.ionic.status = "unconverged".into();
        }
        persist(&path, &s)?;
        emit(&ctx, &settings, &s, out)?;
    }
    if s.rebuilt.is_none() {
        ctx.opts.strain_reference = None;
        ctx.opts.grid_override = if ctx.config.symmetry.mode == crate::config::SymmetryMode::Auto {
            Some(ctx.config.prepare().map_err(|e| error(4, e))?.grid.dims())
        } else {
            ctx.config.basis.fft_grid
        };
        let x = place(&mut ctx, &s.q, a0, l)?;
        match ctx.evaluate(x, 0., None, &settings) {
            Ok((rebuilt, _)) => s.rebuilt = Some(rebuilt),
            Err(e) => {
                s.ionic.status = "failed".into();
                s.ionic
                    .failures
                    .push(json!({"stage":"final_rebuild","reason":e.to_string()}));
                persist(&path, &s)?;
                emit(&ctx, &settings, &s, out)?;
                return Err(e);
            }
        }
        persist(&path, &s)?;
    }
    emit(&ctx, &settings, &s, out)?;
    match s.ionic.status.as_str() {
        "converged" => Ok(()),
        "failed" => Err(error(
            4,
            "variable-cell relaxation failed; last good structure retained",
        )),
        _ => Err(error(3, "variable-cell relaxation max_steps reached")),
    }
}
fn pressure(f: &Frame) -> f64 {
    let t = f.stress_ha_bohr3.unwrap();
    -(t[0][0] + t[1][1] + t[2][2]) / 3.
}
fn emit(ctx: &Context, settings: &Settings, s: &State, out: &Path) -> Result<(), RunError> {
    outputs(ctx, settings, &s.ionic, out)?;
    let last = s.ionic.frames.last().unwrap();
    let mut result: Value =
        serde_json::from_slice(&fs::read(out.join("results.json"))?).map_err(|e| error(5, e))?;
    let block = &mut result["relaxation"];
    block["steps"] = json!(s.accepted_steps);
    block["enthalpy_history_ha"] = json!(
        s.ionic
            .frames
            .iter()
            .map(|f| enthalpy(f, settings.pressure_ha_bohr3))
            .collect::<Vec<_>>()
    );
    block["frozen_final"] = json!({"frame":last,"pressure_ha_bohr3":pressure(last),"basis_protocol":"frozen_miller_set"});
    block["rebuilt_final"]=s.rebuilt.as_ref().map_or(Value::Null,|r|json!({"frame":r,"pressure_ha_bohr3":pressure(r),"basis_protocol":"nominal_cutoff_rebuilt"}));
    if let Some(r) = &s.rebuilt {
        let diff: Matrix = std::array::from_fn(|i| {
            std::array::from_fn(|j| {
                r.stress_ha_bohr3.unwrap()[i][j] - last.stress_ha_bohr3.unwrap()[i][j]
            })
        });
        let max = diff.iter().flatten().map(|v| v.abs()).fold(0., f64::max);
        let force_difference: Vec<[f64; 3]> = r
            .forces_ha_bohr
            .iter()
            .zip(&last.forces_ha_bohr)
            .map(|(a, b)| std::array::from_fn(|i| a[i] - b[i]))
            .collect();
        block["basis_rebuild_difference"] = json!({"forces_ha_bohr":force_difference,"energy_ha":r.energy_ha-last.energy_ha,"stress_ha_bohr3":diff,"pressure_ha_bohr3":pressure(r)-pressure(last),"maximum_stress_difference_ha_bohr3":max,"warning":max>settings.pulay_warning_ha_bohr3,"warning_threshold_ha_bohr3":settings.pulay_warning_ha_bohr3,"scope":"finite-basis diagnostic, not cutoff convergence"});
        let pressure_residual = (pressure(r) - settings.pressure_ha_bohr3).abs();
        let pressure_warning = pressure_residual >= settings.stress_tol_ha_bohr3;
        block["rebuilt_pressure_check"] = json!({
            "residual_ha_bohr3": pressure_residual,
            "target_pressure_ha_bohr3": settings.pressure_ha_bohr3,
            "threshold_ha_bohr3": settings.stress_tol_ha_bohr3,
            "warning": pressure_warning,
            "message": if pressure_warning {
                "Rebuilt-basis pressure exceeds the requested stress tolerance; frozen-basis convergence does not establish nominal-cutoff equilibrium. Use periodic rebuilds and a fresh-basis polish."
            } else { "Rebuilt-basis pressure meets the requested stress tolerance." },
        });
        if pressure_warning {
            eprintln!(
                "warning: final rebuilt-basis pressure residual {pressure_residual:.6e} Ha/Bohr^3 exceeds stress tolerance {:.6e}; use periodic rebuilds and a fresh-basis polish",
                settings.stress_tol_ha_bohr3
            );
        }
        if max > settings.pulay_warning_ha_bohr3 {
            eprintln!("warning: final rebuilt stress differs by {max:.6e} Ha/Bohr^3");
        }
    }
    block["space_group_operations"] = json!(
        dft_geometry::symmetry::detect(&ctx.cell, &ctx.atoms(), ctx.config.symmetry.tolerance_bohr)
            .map_err(|e| error(4, e))?
            .len()
    );
    block["restart_semantics"] = json!(
        "same-build cell optimizer continuation; affine fractional-grid density warm transfer, not exact electronic replay"
    );
    use std::fmt::Write as _;
    let mut log = String::from(
        "frame energy_Ha enthalpy_Ha volume_Bohr3 pressure_Ha_Bohr3 max_force_Ha_Bohr displacement_Bohr provenance\n",
    );
    for (i, f) in s.ionic.frames.iter().enumerate() {
        writeln!(
            log,
            "{i} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {}",
            f.energy_ha,
            enthalpy(f, settings.pressure_ha_bohr3),
            determinant(f.cell_columns_bohr.unwrap()).abs(),
            pressure(f),
            f.max_force_ha_bohr,
            f.step_length_bohr,
            f.provenance
        )
        .unwrap();
    }
    atomic(&out.join("relax.log"), log.as_bytes())?;
    structure_files(ctx, last, out)?;
    atomic(
        &out.join("results.json"),
        &serde_json::to_vec_pretty(&result).map_err(|e| error(5, e))?,
    )
}

// CIF 1.1 core dictionary: lengths Angstrom, angles degrees, fractional sites.
// P1 explicitly carries all atoms; detected group count is metadata in JSON.
fn structure_files(ctx: &Context, frame: &Frame, out: &Path) -> Result<(), RunError> {
    use std::fmt::Write as _;
    let rows = transpose(frame.cell_columns_bohr.unwrap());
    let ang = dft_units::Angstrom::from(dft_units::Bohr(1.)).0;
    let lengths = rows.map(|r| r.iter().map(|x| x * x).sum::<f64>().sqrt());
    let angle = |i: usize, j: usize| {
        (rows[i].iter().zip(rows[j]).map(|(a, b)| a * b).sum::<f64>() / (lengths[i] * lengths[j]))
            .clamp(-1., 1.)
            .acos()
            .to_degrees()
    };
    let mut cif = format!(
        "data_baandi\n_cell_length_a {:.17e}\n_cell_length_b {:.17e}\n_cell_length_c {:.17e}\n_cell_angle_alpha {:.17e}\n_cell_angle_beta {:.17e}\n_cell_angle_gamma {:.17e}\n_symmetry_space_group_name_H-M 'P 1'\n_symmetry_Int_Tables_number 1\nloop_\n_atom_site_label\n_atom_site_type_symbol\n_atom_site_fract_x\n_atom_site_fract_y\n_atom_site_fract_z\n",
        lengths[0] * ang,
        lengths[1] * ang,
        lengths[2] * ang,
        angle(1, 2),
        angle(0, 2),
        angle(0, 1)
    );
    let mut toml = format!(
        "# Replace system and system.atoms; retain other input sections.\n[system]\nname = {:?}\nlattice_unit = \"bohr\"\nlattice_vectors = {:?}\n",
        ctx.config.system.name, rows
    );
    for (i, (atom, f)) in ctx
        .config
        .system
        .atoms
        .iter()
        .zip(&frame.fractional)
        .enumerate()
    {
        let element = ctx
            .potentials
            .iter()
            .find(|(label, _)| label == &atom.species)
            .unwrap()
            .1
            .element();
        writeln!(
            cif,
            "{element}{} {element} {:.17e} {:.17e} {:.17e}",
            i + 1,
            f[0],
            f[1],
            f[2]
        )
        .unwrap();
        writeln!(
            toml,
            "[[system.atoms]]\nspecies = {:?}\nfractional = {:?}",
            atom.species, f
        )
        .unwrap();
    }
    atomic(&out.join("final.cif"), cif.as_bytes())?;
    atomic(&out.join("final.toml"), toml.as_bytes())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn combined_enthalpy_conjugacy_and_shear() {
        // Independent Cartesian harmonic energy and quadratic volume energy.
        // Nonzero shear, skew reference and external P exercise D^-T and A^T.
        let a0 = [[4., 0.7, 0.2], [0., 5., 0.4], [0., 0., 6.]];
        let l = determinant(a0).cbrt();
        let q = vec![
            0.2 * l,
            0.3 * l,
            0.4 * l,
            0.1 * l,
            0.03 * l,
            0.,
            0.03 * l,
            -0.05 * l,
            0.02 * l,
            0.,
            0.02 * l,
            0.04 * l,
        ];
        let settings = Settings {
            mode: Mode::IonsCell,
            pressure_ha_bohr3: 0.013,
            ..Default::default()
        };
        let evaluate = |q: &[f64]| {
            let a = mul(deformation(q, l), a0);
            let v = determinant(a);
            let s = [q[0] / l, q[1] / l, q[2] / l];
            let r: [f64; 3] = std::array::from_fn(|i| (0..3).map(|j| a[i][j] * s[j]).sum());
            let e = 0.5 * r.iter().map(|x| x * x).sum::<f64>() + 0.002 * (v - 100.).powi(2);
            let stress = std::array::from_fn(|i| {
                std::array::from_fn(|j| {
                    r[i] * r[j] / v + if i == j { 0.004 * (v - 100.) } else { 0. }
                })
            });
            Frame {
                positions: r.to_vec(),
                fractional: vec![s],
                energy_ha: e,
                forces_ha_bohr: vec![r.map(|x| -x)],
                stress_ha_bohr3: Some(stress),
                cell_columns_bohr: Some(a),
                pulay_assessment: None,
                max_force_ha_bohr: 0.,
                step_length_bohr: 0.,
                scf_iterations: 0,
                provenance: "analytic oracle".into(),
            }
        };
        for dofs in [Dofs::Full, Dofs::Isotropic, Dofs::XyzDiagonal] {
            let settings = Settings {
                cell_dofs: dofs,
                ..settings.clone()
            };
            let mut direction: Vec<_> = (0..q.len()).map(|i| ((i * 17 + 3) as f64).sin()).collect();
            project(&mut direction, &settings, &[]);
            let h = 1e-5;
            let plus: Vec<_> = q.iter().zip(&direction).map(|(x, d)| x + h * d).collect();
            let minus: Vec<_> = q.iter().zip(&direction).map(|(x, d)| x - h * d).collect();
            let fd = (enthalpy(&evaluate(&plus), settings.pressure_ha_bohr3)
                - enthalpy(&evaluate(&minus), settings.pressure_ha_bohr3))
                / (2. * h);
            let analytic = -force(&evaluate(&q), &q, l, &settings, &[])
                .unwrap()
                .iter()
                .zip(&direction)
                .map(|(f, d)| f * d)
                .sum::<f64>();
            let relative = (fd - analytic).abs() / analytic.abs();
            assert!(
                relative < 1e-7,
                "fd={fd:e} analytic={analytic:e} relative={relative:e}"
            );
        }
    }
}
