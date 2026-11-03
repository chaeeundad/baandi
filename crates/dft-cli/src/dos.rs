//! Weighted Gaussian DOS serialization in portable and familiar layouts.
use crate::bands::{Context, agr_header, agr_world, edges, err, ev};
use dft_cli::run::RunError;
use dft_geometry::kpoints::{monkhorst_pack, reduce_time_reversal};
use serde_json::json;
use std::{fmt::Write as _, path::Path};

pub fn execute(input: &Path, output: Option<&Path>) -> Result<(), RunError> {
    let ctx = Context::load(
        input,
        output,
        &[
            "dos.csv", "dos.dat", "dos.agr", "dos.json", "DOSCAR", "dos.svg",
        ],
    )?;
    let default = dft_cli::config::Dos::default();
    let d = ctx.config.dos.as_ref().unwrap_or(&default);
    if d.dos_kgrid.is_none() && ctx.config.basis.kgrid.is_none() && !ctx.config.basis.legacy_gamma()
    {
        return Err(err(
            1,
            "dos_kgrid is required when the SCF uses an explicit non-Gamma k list",
        ));
    }
    let grid = if let Some(grid) = d.dos_kgrid {
        grid
    } else {
        let source = ctx.config.basis.kgrid.unwrap_or([1; 3]);
        let mut grid = [0; 3];
        for i in 0..3 {
            grid[i] = source[i]
                .checked_mul(2)
                .ok_or_else(|| err(1, "DOS grid overflow"))?;
        }
        grid
    };
    let sigma = d.sigma_ev.unwrap_or(0.1);
    let spacing = d.spacing_ev.unwrap_or(0.01);
    if !sigma.is_finite() || sigma <= 0.0 || !spacing.is_finite() || spacing <= 0.0 {
        return Err(err(
            1,
            "sigma_ev and spacing_ev must be finite and positive",
        ));
    }
    if d.window_ev
        .is_some_and(|w| !w.iter().all(|e| e.is_finite()) || w[1] <= w[0])
    {
        return Err(err(1, "window_ev must be a finite increasing pair"));
    }
    let points = reduce_time_reversal(
        &monkhorst_pack(grid, ctx.config.basis.kshift.unwrap_or([0.0; 3]))
            .map_err(|e| err(1, e))?,
    )
    .map_err(|e| err(1, e))?;
    let scf = ctx.scf()?;
    let result = ctx.spectrum(
        &scf,
        &points,
        d.nbands,
        d.eigen_tolerance_ha,
        d.max_iterations,
    )?;
    let edge = edges(&result)?;
    let spectra: Vec<Vec<f64>> = result
        .solutions
        .iter()
        .map(|k| {
            k.solution
                .eigenvalues
                .iter()
                .map(|e| ev(e - edge.vbm))
                .collect()
        })
        .collect();
    let min = spectra
        .iter()
        .flatten()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let max = spectra
        .iter()
        .flatten()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    let window = d
        .window_ev
        .unwrap_or([min - 8.0 * sigma, max + 8.0 * sigma]);
    let intervals = ((window[1] - window[0]) / spacing).ceil();
    if !intervals.is_finite() || intervals > 1_000_000.0 || intervals < 1.0 {
        return Err(err(1, "DOS energy grid requires 1..1000000 intervals"));
    }
    let intervals = intervals as usize;
    let energy: Vec<_> = (0..=intervals)
        .map(|i| window[0] + (window[1] - window[0]) * i as f64 / intervals as f64)
        .collect();
    let weights: Vec<_> = points.iter().map(|p| p.weight).collect();
    let rows =
        dft_engine::gaussian_dos(&spectra, &weights, &energy, sigma).map_err(|e| err(1, e))?;
    let mut csv = String::from(
        "# energy eV relative to sampled VBM; DOS states/eV/cell including spin factor 2; integral states/cell\nenergy_eV,dos_states_per_eV_per_cell,integrated_dos\n",
    );
    let mut dat = String::from(
        "# energy_eV_relative_to_VBM DOS_states_per_eV_per_cell integrated_states_per_cell (spin factor 2)\n",
    );
    let mut agr = agr_header("Energy relative to VBM (eV)", "DOS (states/eV/cell)");
    let upper_dos = rows.iter().map(|r| r[1]).fold(0.0, f64::max).max(1e-12);
    agr.push_str(&agr_world(window[0], window[1], 0.0, 1.05 * upper_dos));
    agr.push_str("@target G0.S0\n@type xy\n");
    let mut doscar = ctx.vasp_header(false);
    writeln!(
        doscar,
        "{:.16e} {:.16e} {} 0.0 1.0",
        window[1],
        window[0],
        rows.len()
    )
    .unwrap();
    for r in &rows {
        writeln!(csv, "{:.16e},{:.16e},{:.16e}", r[0], r[1], r[2]).unwrap();
        writeln!(dat, "{:.16e} {:.16e} {:.16e}", r[0], r[1], r[2]).unwrap();
        writeln!(doscar, "{:.16e} {:.16e} {:.16e}", r[0], r[1], r[2]).unwrap();
        writeln!(agr, "{:.16e} {:.16e}", r[0], r[1]).unwrap();
    }
    agr.push_str("&\n");
    let mut meta = ctx.metadata(&result, &edge, &scf);
    meta["gaussian"] = json!({"sigma_ev":sigma,"formula":"sum_kn 2*w_k/(sigma*sqrt(2*pi))*exp(-(E-e_kn)^2/(2*sigma^2))","integrated_dos":"analytic Gaussian CDF from minus infinity","spin_degeneracy":2,"window_ev":window,"spacing_ev":(window[1]-window[0])/intervals as f64,"total_states":2*spectra[0].len()});
    meta["dos_kgrid"] = json!(grid);
    meta["time_reversal"] = json!(true);
    meta["kshift"] = json!(ctx.config.basis.kshift.unwrap_or([0.0; 3]));
    meta["DOSCAR_fermi_slot"] =
        json!("0 is the VBM reference placeholder, NOT a computed Fermi energy");
    for (name, text) in [
        ("dos.csv", csv),
        ("dos.dat", dat),
        ("dos.agr", agr),
        ("DOSCAR", doscar),
    ] {
        ctx.write(name, &text)?;
    }
    ctx.json("dos.json", &meta)?;
    let maxdos = rows.iter().map(|r| r[1]).fold(0.0, f64::max);
    let mut svg = String::from(
        "<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"900\" height=\"600\" viewBox=\"0 0 900 600\"><rect width=\"900\" height=\"600\" fill=\"white\"/><g font-family=\"sans-serif\" font-size=\"16\"><text x=\"80\" y=\"30\">Gaussian DOS (states/eV/cell; spin included)</text><text x=\"300\" y=\"580\">Energy relative to VBM (eV)</text></g><path d=\"M80 60 V520 H840\" fill=\"none\" stroke=\"black\"/>",
    );
    for i in 0..=5 {
        let x = 80.0 + 760.0 * i as f64 / 5.0;
        let e = window[0] + (window[1] - window[0]) * i as f64 / 5.0;
        write!(
            svg,
            "<text x=\"{x}\" y=\"545\" text-anchor=\"middle\">{e:.2}</text>"
        )
        .unwrap();
        let y = 520.0 - 460.0 * i as f64 / 5.0;
        write!(
            svg,
            "<text x=\"70\" y=\"{y}\" text-anchor=\"end\">{:.2}</text>",
            maxdos * i as f64 / 5.0
        )
        .unwrap();
    }
    svg.push_str("<polyline fill=\"none\" stroke=\"#b44623\" points=\"");
    for r in &rows {
        write!(
            svg,
            "{:.5},{:.5} ",
            80.0 + 760.0 * (r[0] - window[0]) / (window[1] - window[0]),
            520.0 - 460.0 * r[1] / maxdos.max(f64::MIN_POSITIVE)
        )
        .unwrap();
    }
    svg.push_str("\"/></svg>");
    ctx.write("dos.svg", &svg)?;
    println!(
        "DOS: {} weighted k points, {} energies; sampled gap {:.9} eV; integrated window states {:.9}",
        points.len(),
        rows.len(),
        ev(edge.cbm - edge.vbm),
        rows.last().unwrap()[2] - rows[0][2]
    );
    Ok(())
}
