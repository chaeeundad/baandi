//! Portable spectral output; all numerical work remains in the engine.
use dft_cli::{
    config::Config,
    run::{RunError, options, sha256},
};
use dft_engine::{
    NscfResult, OccupationPolicy, ScfOptions, ScfProblem, ScfResult, SpeciesPotential,
};
use dft_geometry::{
    Atom, Cell,
    kpath::{self, KPath},
    kpoints::KPoint,
};
use dft_units::{Angstrom, Bohr, EnergyUnit, from_hartree};
use serde_json::{Value, json};
use std::{
    fmt::Write as _,
    fs,
    path::{Path, PathBuf},
};

pub(crate) fn err(code: u8, cause: impl std::fmt::Display) -> RunError {
    RunError {
        code,
        cause: cause.to_string(),
    }
}
pub(crate) fn ev(ha: f64) -> f64 {
    from_hartree(ha, EnergyUnit::ElectronVolt)
}
pub(crate) struct Context {
    pub config: Config,
    pub cell: Cell,
    pub atoms: Vec<Atom>,
    pub potentials: Vec<(String, dft_upf::Pseudopotential)>,
    pub options: ScfOptions,
    pub out: PathBuf,
    pub input_hash: String,
    pub pp_metadata: Value,
}
impl Context {
    pub fn load(input: &Path, output: Option<&Path>, names: &[&str]) -> Result<Self, RunError> {
        let bytes = fs::read(input)?;
        let config = Config::parse(std::str::from_utf8(&bytes).map_err(|e| err(1, e))?)
            .map_err(|e| err(1, e))?;
        let cell = config.cell().map_err(|e| err(1, e))?;
        let options = options(&config)?;
        let atoms = config
            .system
            .atoms
            .iter()
            .map(|a| Atom::new(&a.species, a.fractional))
            .collect();
        let mut potentials = Vec::new();
        let mut inputs = vec![fs::canonicalize(input)?];
        let mut pp = Vec::new();
        for (label, path) in &config.pseudopotentials {
            let path = input.parent().unwrap_or(Path::new(".")).join(path);
            let data = fs::read(&path)?;
            let p = dft_upf::parse_upf(std::str::from_utf8(&data).map_err(|e| err(2, e))?)
                .map_err(|e| err(2, e))?;
            inputs.push(fs::canonicalize(&path)?);
            pp.push(json!({"species":label,"path":path,"sha256":sha256(&data),"functional":p.functional()}));
            potentials.push((label.clone(), p));
        }
        let out = output.map(Path::to_owned).unwrap_or_else(|| {
            PathBuf::from(
                config
                    .output
                    .as_ref()
                    .and_then(|o| o.directory.as_deref())
                    .unwrap_or("results"),
            )
        });
        fs::create_dir_all(&out)?;
        // Refuse symlinks and multiply linked outputs, including output/output
        // aliases. Check every destination before any computation or writing.
        for name in names {
            let path = out.join(name);
            if let Ok(meta) = fs::symlink_metadata(&path) {
                if !meta.file_type().is_file() {
                    return Err(err(
                        5,
                        format!("refusing nonregular output {}", path.display()),
                    ));
                }
                #[cfg(unix)]
                {
                    use std::os::unix::fs::MetadataExt;
                    if meta.nlink() > 1 {
                        return Err(err(5, "refusing hard-linked output"));
                    }
                }
                if inputs.contains(&fs::canonicalize(&path)?) {
                    return Err(err(5, "refusing to overwrite an input"));
                }
            }
        }
        Ok(Self {
            config,
            cell,
            atoms,
            potentials,
            options,
            out,
            input_hash: sha256(&bytes),
            pp_metadata: json!(pp),
        })
    }
    fn refs(&self) -> Vec<SpeciesPotential<'_>> {
        self.potentials
            .iter()
            .map(|(label, potential)| SpeciesPotential { label, potential })
            .collect()
    }
    pub fn scf(&self) -> Result<ScfResult, RunError> {
        let points = self.config.basis.resolve_kpoints().map_err(|e| err(1, e))?;
        let problem = if self.config.basis.legacy_gamma() {
            ScfProblem::new(&self.cell, &self.atoms, &self.refs(), self.options.clone())
        } else {
            ScfProblem::with_kpoints(
                &self.cell,
                &self.atoms,
                &self.refs(),
                self.options.clone(),
                &points,
                OccupationPolicy::Insulator {
                    minimum_gap_hartree: 1e-7,
                },
            )
        }
        .map_err(|e| err(4, e))?;
        let result = problem.run().map_err(|e| {
            err(
                if matches!(e, dft_engine::ScfError::NotConverged { .. }) {
                    3
                } else {
                    4
                },
                e,
            )
        })?;
        println!(
            "SCF converged: {:.12} Ha, {} iterations",
            result.total_hartree(),
            result.history.len()
        );
        Ok(result)
    }
    pub fn spectrum(
        &self,
        scf: &ScfResult,
        points: &[KPoint],
        bands: Option<usize>,
        tolerance: Option<f64>,
        iterations: Option<usize>,
    ) -> Result<NscfResult, RunError> {
        let occupied = scf.occupations.iter().filter(|&&f| f > 0.0).count();
        let mut opts = self.options.clone();
        let count = bands.unwrap_or(occupied + opts.extra_bands.max(4));
        if count < occupied + 4 {
            return Err(err(1, "NSCF requires occupied + at least four empty bands"));
        }
        opts.bands = Some(count);
        opts.eigen_tolerance = tolerance.unwrap_or(1e-8);
        opts.davidson_max_iterations = iterations.unwrap_or(opts.davidson_max_iterations);
        opts.davidson_subspace = opts.davidson_subspace.max(count.saturating_mul(4));
        println!(
            "NSCF: {} k points, {count} bands, residual tolerance {} Ha; all bands strict (degenerate top bands may need more bands/iterations)",
            points.len(),
            opts.eigen_tolerance
        );
        dft_engine::nscf(
            &self.cell,
            &self.atoms,
            &self.refs(),
            &scf.density,
            points,
            opts,
        )
        .map_err(|e| err(4, e))
    }
    pub fn metadata(&self, result: &NscfResult, edges: &Edges, scf: &ScfResult) -> Value {
        let density_bytes: Vec<_> = scf.density.iter().flat_map(|v| v.to_le_bytes()).collect();
        json!({"schema_version":1,"scf_density_sha256":sha256(&density_bytes),"scf_energy_hartree":scf.total_hartree(),"scf_iterations":scf.history.len(),"eigen_tolerance_hartree":result.tolerance_hartree,"title":self.config.system.name,"code":"rust-dft","code_version":env!("CARGO_PKG_VERSION"),"input_sha256":self.input_hash,"pseudopotentials":self.pp_metadata,
        "energy_reference":"maximum occupied eigenvalue over this computed sampling (VBM); not a Fermi level","hartree_to_ev":ev(1.0),"energy_reference_hartree":edges.vbm,"vbm":{"hartree":edges.vbm,"k_index":edges.vi,"fractional":result.points[edges.vi].fractional},"cbm":{"hartree":edges.cbm,"k_index":edges.ci,"fractional":result.points[edges.ci].fractional},"gap_ev":ev(edges.cbm-edges.vbm),"gap_hartree":edges.cbm-edges.vbm,"electrons":result.electrons,"nbands":result.solutions[0].solution.eigenvalues.len(),"cutoff_hartree":self.options.cutoff_hartree,"threads":self.options.threads,"converged":result.converged,"functional":self.options.functional.name(),"lattice_columns_bohr":self.cell.lattice_columns_bohr(),"kpoints":result.points.iter().map(|p| json!({"fractional":p.fractional,"weight":p.weight})).collect::<Vec<_>>(),"residuals_hartree":result.solutions.iter().map(|k| &k.solution.residuals).collect::<Vec<_>>()})
    }
    pub fn write(&self, name: &str, text: &str) -> Result<(), RunError> {
        fs::write(self.out.join(name), text)?;
        Ok(())
    }
    pub fn json(&self, name: &str, value: &Value) -> Result<(), RunError> {
        self.write(
            name,
            &serde_json::to_string_pretty(value).map_err(|e| err(5, e))?,
        )
    }
    pub fn vasp_header(&self, eigenval: bool) -> String {
        // EIGENVAL: NBLOCK*KBLOCK=1, ISPIN=1. DOSCAR: no partial DOS,
        // unused NCDIJ=0. Remaining header fields have the same layout.
        let flags = if eigenval { "1 1" } else { "0 0" };
        let n = self.atoms.len();
        let a = self.cell.lattice_columns_bohr();
        let lengths: [f64; 3] = std::array::from_fn(|j| {
            Angstrom::from(Bohr((0..3).map(|i| a[i][j] * a[i][j]).sum::<f64>().sqrt())).0 * 1e-10
        });
        let volume = self.cell.volume_bohr3() * Angstrom::from(Bohr(1.0)).0.powi(3);
        format!(
            "{n} {n} {flags}\n{volume:.16e} {:.16e} {:.16e} {:.16e} 0.0\n0.0\nCAR\nrust-dft VBM-referenced export (not VASP)\n",
            lengths[0], lengths[1], lengths[2]
        )
    }
}
pub(crate) struct Edges {
    pub vbm: f64,
    pub cbm: f64,
    pub vi: usize,
    pub ci: usize,
}
pub(crate) fn edges(result: &NscfResult) -> Result<Edges, RunError> {
    let occupied = result.electrons.round() as usize / 2;
    let mut e = Edges {
        vbm: f64::NEG_INFINITY,
        cbm: f64::INFINITY,
        vi: 0,
        ci: 0,
    };
    for (i, k) in result.solutions.iter().enumerate() {
        if k.solution.eigenvalues[occupied - 1] > e.vbm {
            e.vbm = k.solution.eigenvalues[occupied - 1];
            e.vi = i;
        }
        if k.solution.eigenvalues[occupied] < e.cbm {
            e.cbm = k.solution.eigenvalues[occupied];
            e.ci = i;
        }
    }
    if e.cbm <= e.vbm {
        return Err(err(
            4,
            "computed sampling has overlapping bands; insulating VBM reference is unavailable (smearing unsupported)",
        ));
    }
    Ok(e)
}
pub(crate) fn agr_world(xmin: f64, xmax: f64, ymin: f64, ymax: f64) -> String {
    format!(
        "@ world xmin {xmin:.16e}\n@ world xmax {xmax:.16e}\n@ world ymin {ymin:.16e}\n@ world ymax {ymax:.16e}\n@ xaxis tick major {:.16e}\n@ yaxis tick major {:.16e}\n",
        (xmax - xmin) / 5.0,
        (ymax - ymin) / 5.0
    )
}

pub(crate) fn agr_header(x: &str, y: &str) -> String {
    format!(
        "# rust-dft portable Grace project\n@version 50122\n@with g0\n@ xaxis label \"{x}\"\n@ yaxis label \"{y}\"\n"
    )
}

fn path(config: &Config, cell: &Cell) -> Result<KPath, RunError> {
    let b = config
        .bands
        .as_ref()
        .ok_or_else(|| err(1, "bands requires a [bands] table"))?;
    let vertices = match (&b.preset, &b.points) {
        (Some(name), None) => kpath::preset(cell, name).map_err(|e| err(1, e))?,
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
    if b.points_per_segment.is_some() && b.spacing_inv_bohr.is_some() {
        return Err(err(1, "choose points_per_segment or spacing_inv_bohr"));
    }
    kpath::generate(
        cell,
        &vertices,
        b.points_per_segment.unwrap_or(20),
        b.spacing_inv_bohr,
    )
    .map_err(|e| err(1, e))
}
fn csv_label(s: &str) -> String {
    format!("\"{}\"", s.replace('"', "\"\""))
}

pub fn execute(input: &Path, output: Option<&Path>) -> Result<(), RunError> {
    let ctx = Context::load(
        input,
        output,
        &[
            "bands.csv",
            "bands_labels.csv",
            "bands.json",
            "bands.dat",
            "bands.dat.gnu",
            "bands.agr",
            "EIGENVAL",
            "bands.svg",
        ],
    )?;
    let path = path(&ctx.config, &ctx.cell)?;
    let b = ctx.config.bands.as_ref().unwrap();
    let points: Vec<_> = path
        .fractional
        .iter()
        .map(|&fractional| KPoint {
            fractional,
            weight: 1.0 / path.fractional.len() as f64,
        })
        .collect();
    let scf = ctx.scf()?;
    let result = ctx.spectrum(
        &scf,
        &points,
        b.nbands,
        b.eigen_tolerance_ha,
        b.max_iterations,
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
    let nb = spectra[0].len();
    let mut csv = String::from(
        "# k fractional; distance 1/Bohr; energies eV relative to sampled VBM\nk_index,distance_inv_bohr,kx,ky,kz",
    );
    for band in 1..=nb {
        write!(csv, ",band_{band}").unwrap();
    }
    csv.push('\n');
    let mut dat = String::from(
        "# distance_inv_bohr band_1..band_N_eV_relative_to_VBM; blank lines mark path breaks\n",
    );
    let mut eigenval = ctx.vasp_header(true);
    writeln!(eigenval, "{:.0} {} {nb}", result.electrons, points.len()).unwrap();
    for (i, (p, energies)) in points.iter().zip(&spectra).enumerate() {
        if path.breaks.contains(&i) {
            dat.push('\n');
        }
        write!(
            csv,
            "{i},{:.16e},{:.16e},{:.16e},{:.16e}",
            path.distances[i], p.fractional[0], p.fractional[1], p.fractional[2]
        )
        .unwrap();
        write!(dat, "{:.16e}", path.distances[i]).unwrap();
        writeln!(
            eigenval,
            "\n{:.16e} {:.16e} {:.16e} 0.0",
            p.fractional[0], p.fractional[1], p.fractional[2]
        )
        .unwrap();
        for (band, &energy) in energies.iter().enumerate() {
            write!(csv, ",{energy:.16e}").unwrap();
            write!(dat, " {energy:.16e}").unwrap();
            writeln!(
                eigenval,
                "{} {energy:.16e} {:.1}",
                band + 1,
                if band < result.electrons as usize / 2 {
                    2.0
                } else {
                    0.0
                }
            )
            .unwrap();
        }
        csv.push('\n');
        dat.push('\n');
    }
    let mut labels = String::from("label,k_index,distance\n");
    for l in &path.labels {
        writeln!(
            labels,
            "{},{},{:.16e}",
            csv_label(&l.label),
            l.index,
            l.distance
        )
        .unwrap();
    }
    let mut gnu = String::from(
        "# rust-dft: x in 1/Bohr; energy in eV relative to VBM; one band/branch per block\n",
    );
    let mut agr = agr_header("Path (1/Bohr)", "Energy relative to VBM (eV)");
    let min_energy = spectra
        .iter()
        .flatten()
        .copied()
        .fold(f64::INFINITY, f64::min);
    let max_energy = spectra
        .iter()
        .flatten()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);
    agr.push_str(&agr_world(
        0.0,
        *path.distances.last().unwrap(),
        min_energy - 0.1,
        max_energy + 0.1,
    ));
    let mut set = 0;
    let mut starts = vec![0];
    starts.extend(&path.breaks);
    starts.push(points.len());
    for band in 0..nb {
        for range in starts.windows(2) {
            writeln!(agr, "@target G0.S{set}\n@type xy").unwrap();
            set += 1;
            for (i, energies) in spectra.iter().enumerate().take(range[1]).skip(range[0]) {
                writeln!(gnu, "{:.16e} {:.16e}", path.distances[i], energies[band]).unwrap();
                writeln!(agr, "{:.16e} {:.16e}", path.distances[i], energies[band]).unwrap();
            }
            gnu.push('\n');
            agr.push_str("&\n");
        }
    }
    let mut meta = ctx.metadata(&result, &edge, &scf);
    meta["path_break_indices"] = json!(path.breaks);
    meta["path_labels"] = json!(
        path.labels
            .iter()
            .map(|l| json!({"label":l.label,"k_index":l.index,"distance_inv_bohr":l.distance}))
            .collect::<Vec<_>>()
    );
    meta["sampling_role"] =
        json!("path points are not a Brillouin-zone quadrature; EIGENVAL weights are zero");
    for (name, text) in [
        ("bands.csv", csv),
        ("bands_labels.csv", labels),
        ("bands.dat", dat),
        ("bands.dat.gnu", gnu),
        ("bands.agr", agr),
        ("EIGENVAL", eigenval),
    ] {
        ctx.write(name, &text)?;
    }
    ctx.json("bands.json", &meta)?;
    let status = std::process::Command::new("python3")
        .args([
            "-I",
            "-c",
            include_str!("../../../scripts/plot_bands_dos.py"),
        ])
        .arg(ctx.out.join("bands.csv"))
        .arg("-o")
        .arg(ctx.out.join("bands.svg"))
        .status()?;
    if !status.success() {
        return Err(err(5, "SVG generation failed; Python 3 required"));
    }
    println!(
        "VBM {:.12} Ha at {:?}; CBM {:.12} Ha at {:?}; sampled gap {:.9} eV",
        edge.vbm,
        points[edge.vi].fractional,
        edge.cbm,
        points[edge.ci].fractional,
        ev(edge.cbm - edge.vbm)
    );
    Ok(())
}
