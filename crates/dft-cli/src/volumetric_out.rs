//! Post-SCF visualization I/O; cube is Bohr/z-fastest, XSF Angstrom/x-fastest.
use crate::config::{Config, PreparedInput, Volumetric};
use dft_basis::{FftGridSpec, PlaneWaveBasis};
use dft_engine::{ScfResult, SpeciesPotential};
use dft_geometry::Cell;
use dft_physics::{data::*, hamiltonian, volumetric};
use std::{
    fs::OpenOptions,
    io::{self, BufWriter, Write},
    path::Path,
};

type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

const SYMBOLS: &str = "H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og";

/// IUPAC symbols in atomic-number order (all 118); exact case, no label guessing.
/// Source: https://iupac.org/what-we-do/periodic-table-of-elements/
pub fn atomic_number(symbol: &str) -> Result<usize> {
    SYMBOLS
        .split_whitespace()
        .position(|s| s == symbol)
        .map(|i| i + 1)
        .ok_or_else(|| format!("unknown element symbol {symbol:?}").into())
}

pub struct OutputAtom {
    pub atomic_number: usize,
    pub valence_charge: f64,
    pub position_bohr: [f64; 3],
}
fn validate(data: &Array3<f64>, atoms: &[OutputAtom]) -> Result<()> {
    if data.is_empty()
        || data.iter().any(|x| !x.is_finite())
        || atoms.iter().any(|a| {
            !(1..=118).contains(&a.atomic_number)
                || !a.valence_charge.is_finite()
                || a.position_bohr.iter().any(|x| !x.is_finite())
        })
    {
        return Err("invalid volumetric samples or atoms".into());
    }
    Ok(())
}

/// 17 significant decimal digits retain f64 samples on a decimal round trip.
pub fn write_cube(
    mut w: impl Write,
    cell: &Cell,
    atoms: &[OutputAtom],
    data: &Array3<f64>,
    label: &str,
) -> Result<()> {
    validate(data, atoms)?;
    if label.contains(['\n', '\r']) {
        return Err("multiline label".into());
    }
    writeln!(
        w,
        "rust-dft {label}\nBohr coordinates; scalar atomic units; z fastest"
    )?;
    writeln!(w, "{} 0.0 0.0 0.0", atoms.len())?;
    let a = cell.lattice_columns_bohr();
    for (axis, &n) in data.shape().iter().enumerate() {
        writeln!(
            w,
            "{n} {:.16e} {:.16e} {:.16e}",
            a[0][axis] / n as f64,
            a[1][axis] / n as f64,
            a[2][axis] / n as f64
        )?;
    }
    for atom in atoms {
        writeln!(
            w,
            "{} {:.16e} {:.16e} {:.16e} {:.16e}",
            atom.atomic_number,
            atom.valence_charge,
            atom.position_bohr[0],
            atom.position_bohr[1],
            atom.position_bohr[2]
        )?;
    }
    for (i, v) in data.iter().enumerate() {
        write!(w, " {v:.16e}")?;
        if (i + 1) % 6 == 0 {
            writeln!(w)?;
        }
    }
    writeln!(w)?;
    Ok(())
}

/// XSF general grid: N+1 nodes along each axis, duplicated periodic endpoints.
/// Geometry is Angstrom; field values retain atomic units (no density rescale).
pub fn write_xsf(
    mut w: impl Write,
    cell: &Cell,
    atoms: &[OutputAtom],
    data: &Array3<f64>,
    label: &str,
) -> Result<()> {
    validate(data, atoms)?;
    if label.is_empty() || !label.chars().all(|c| c.is_ascii_alphanumeric() || c == '_') {
        return Err("invalid XSF label".into());
    }
    let a = cell.lattice_columns_bohr();
    let ang = dft_units::BOHR_IN_ANGSTROM;
    writeln!(w, "# {label}: scalar atomic units\nCRYSTAL\nPRIMVEC")?;
    for ((x, y), z) in a[0].iter().zip(&a[1]).zip(&a[2]) {
        writeln!(w, "{:.16e} {:.16e} {:.16e}", x * ang, y * ang, z * ang)?;
    }
    writeln!(w, "PRIMCOORD\n{} 1", atoms.len())?;
    for atom in atoms {
        writeln!(
            w,
            "{} {:.16e} {:.16e} {:.16e}",
            atom.atomic_number,
            atom.position_bohr[0] * ang,
            atom.position_bohr[1] * ang,
            atom.position_bohr[2] * ang
        )?;
    }
    let d = data.shape();
    writeln!(
        w,
        "BEGIN_BLOCK_DATAGRID_3D\n{label}\nBEGIN_DATAGRID_3D_{label}\n{} {} {}\n0.0 0.0 0.0",
        d[0] + 1,
        d[1] + 1,
        d[2] + 1
    )?;
    for ((x, y), z) in a[0].iter().zip(&a[1]).zip(&a[2]) {
        writeln!(w, "{:.16e} {:.16e} {:.16e}", x * ang, y * ang, z * ang)?;
    }
    for k in 0..=d[2] {
        for j in 0..=d[1] {
            for i in 0..=d[0] {
                writeln!(w, "{:.16e}", data[[i % d[0], j % d[1], k % d[2]]])?;
            }
        }
    }
    writeln!(w, "END_DATAGRID_3D\nEND_BLOCK_DATAGRID_3D")?;
    Ok(())
}

fn selection(config: &Config) -> Option<&Volumetric> {
    config.output.as_ref()?.volumetric.as_ref()
}

/// Reject existing files/symlinks before SCF and never overwrite them at write.
pub fn preflight(config: &Config, out: &Path) -> Result<()> {
    if let Some(v) = selection(config) {
        v.validate()?;
        for q in &v.quantities {
            for f in &v.formats {
                let path = out.join(output_name(q, f));
                match path.symlink_metadata() {
                    Ok(_) => return Err(format!("refusing to overwrite {}", path.display()).into()),
                    Err(e) if e.kind() == io::ErrorKind::NotFound => {}
                    Err(e) => return Err(e.into()),
                }
            }
        }
    }
    Ok(())
}

pub fn write_outputs(
    config: &Config,
    prepared: &PreparedInput,
    potentials: &[SpeciesPotential<'_>],
    result: &ScfResult,
    out: &Path,
) -> Result<()> {
    let Some(v) = selection(config) else {
        return Ok(());
    };
    v.validate()?;
    let cell = &prepared.cell;
    let dims = [
        result.density.shape()[0],
        result.density.shape()[1],
        result.density.shape()[2],
    ];
    let atoms = config
        .system
        .atoms
        .iter()
        .map(|a| {
            let p = potentials
                .iter()
                .find(|p| p.label == a.species)
                .ok_or("missing species")?
                .potential;
            Ok(OutputAtom {
                atomic_number: atomic_number(p.element().trim())?,
                valence_charge: p.z_valence(),
                position_bohr: cell.fractional_to_cartesian_bohr(a.fractional),
            })
        })
        .collect::<Result<Vec<_>>>()?;
    let mut fields = Vec::new();
    for q in &v.quantities {
        let data = match q.as_str() {
            "density" => result.density.clone(),
            "hartree" => volumetric::hartree_potential(&result.density, cell)?,
            "local_potential" => {
                let positions: Vec<Vec<_>> = potentials
                    .iter()
                    .map(|p| {
                        config
                            .system
                            .atoms
                            .iter()
                            .filter(|a| a.species == p.label)
                            .map(|a| cell.fractional_to_cartesian_bohr(a.fractional))
                            .collect()
                    })
                    .collect();
                let species: Vec<_> = potentials
                    .iter()
                    .zip(&positions)
                    .map(|(p, r)| hamiltonian::Species {
                        potential: p.potential,
                        positions_bohr: r,
                    })
                    .collect();
                let grid = FftGridSpec::new(&prepared.density, &prepared.density, Some(dims))?;
                let reciprocal = grid.scatter(&hamiltonian::ionic_local(
                    &prepared.density,
                    cell.volume_bohr3(),
                    &species,
                )?)?;
                let mut fft = FftPlan3d::new(dims)?;
                let mut ionic = Array3::zeros(dims);
                fft.to_real(&reciprocal, &mut ionic)?;
                let core = dft_physics::density::core_density(
                    &prepared.density,
                    &grid,
                    cell.volume_bohr3(),
                    &species,
                    &mut fft,
                )?;
                let mut xc_density = Array3::zeros(dims);
                let xc = config.xc.functional.engine().evaluate_with_core(
                    &result.density,
                    core.as_ref(),
                    &mut xc_density,
                    cell,
                    &mut fft,
                )?;
                volumetric::total_local(
                    &ionic.mapv(|z| z.re),
                    &volumetric::hartree_potential(&result.density, cell)?,
                    &xc.potential,
                )?
            }
            "elf" => {
                let bases = result
                    .k_points
                    .iter()
                    .map(|k| {
                        PlaneWaveBasis::at_k(
                            cell,
                            prepared.wavefunctions.cutoff(),
                            k.point.fractional,
                        )
                    })
                    .collect::<std::result::Result<Vec<_>, _>>()?;
                let coefficients: Vec<Vec<Vec<_>>> = result
                    .k_points
                    .iter()
                    .map(|k| {
                        let m = &k.solution.eigenvectors;
                        (0..m.cols())
                            .map(|j| (0..m.rows()).map(|i| m[(i, j)]).collect())
                            .collect()
                    })
                    .collect();
                let blocks: Vec<_> = bases
                    .iter()
                    .zip(&coefficients)
                    .zip(&result.k_points)
                    .map(|((b, c), k)| volumetric::OrbitalBlock {
                        basis: b,
                        coefficients: c,
                        occupations: &result.occupations,
                        weight: k.point.weight,
                    })
                    .collect();
                volumetric::elf_from_orbitals(&result.density, &blocks, cell)?
            }
            "partial_charge" => {
                let k = result
                    .k_points
                    .get(v.kpoint.unwrap_or(1) - 1)
                    .ok_or("selected k point out of range")?;
                let m = &k.solution.eigenvectors;
                let [first, last] = v.bands.unwrap_or([1, 1]);
                if last > m.cols() {
                    return Err("selected band out of range".into());
                }
                let coefficients: Vec<Vec<_>> = (first - 1..last)
                    .map(|j| (0..m.rows()).map(|i| m[(i, j)]).collect())
                    .collect();
                let basis = PlaneWaveBasis::at_k(
                    cell,
                    prepared.wavefunctions.cutoff(),
                    k.point.fractional,
                )?;
                let grid = FftGridSpec::new(&basis, &prepared.density, Some(dims))?;
                volumetric::partial_charge(&coefficients, &grid, cell.volume_bohr3())?
            }
            _ => return Err("unknown volumetric quantity".into()),
        };
        fields.push((q, data));
    }
    for (q, data) in fields {
        for f in &v.formats {
            let mut file = BufWriter::new(
                OpenOptions::new()
                    .write(true)
                    .create_new(true)
                    .open(out.join(output_name(q, f)))?,
            );
            match f.as_str() {
                "cube" => write_cube(&mut file, cell, &atoms, &data, q)?,
                "xsf" => write_xsf(&mut file, cell, &atoms, &data, q)?,
                "vasp" => write_vasp(&mut file, cell, &atoms, &data, q)?,
                "vti" => write_vti(&mut file, cell, &data, q)?,
                _ => return Err("unknown volumetric format".into()),
            }
            file.flush()?;
        }
    }
    Ok(())
}

/// Viewer filename; VASP suffix records the different scalar normalization.
pub fn output_name(quantity: &str, format: &str) -> String {
    let extension = if format == "vasp" {
        match quantity {
            "density" | "partial_charge" => "CHGCAR",
            "elf" => "ELFCAR",
            _ => "LOCPOT",
        }
    } else {
        format
    };
    format!("{quantity}.{extension}")
}

/// VASP visualization layout, from public POSCAR/CHGCAR/LOCPOT/ELFCAR docs.
/// No VASP implementation was consulted. No PAW augmentation/restart payload.
/// Density samples are n[Bohr^-3]*Omega[Bohr^3], potentials eV, ELF dimensionless.
pub fn write_vasp(
    mut w: impl Write,
    cell: &Cell,
    atoms: &[OutputAtom],
    data: &Array3<f64>,
    quantity: &str,
) -> Result<()> {
    validate(data, atoms)?;
    let (scale, units) = match quantity {
        "density" | "partial_charge" => (cell.volume_bohr3(), "density times cell volume"),
        "hartree" | "local_potential" => (
            dft_units::from_hartree(1.0, dft_units::EnergyUnit::ElectronVolt),
            "eV",
        ),
        "elf" => (1.0, "dimensionless"),
        _ => return Err("unknown VASP scalar quantity".into()),
    };
    if atoms.is_empty() || data.iter().any(|x| !(x * scale).is_finite()) {
        return Err("VASP output requires atoms and representable scaled values".into());
    }
    writeln!(
        w,
        "rust-dft {quantity}; {units}; valence-only visualization\n1.0"
    )?;
    let a = cell.lattice_columns_bohr();
    let ang = dft_units::BOHR_IN_ANGSTROM;
    for ((x, y), z) in a[0].iter().zip(&a[1]).zip(&a[2]) {
        writeln!(w, "{:.16e} {:.16e} {:.16e}", x * ang, y * ang, z * ang)?;
    }
    // Group atoms by element in first-appearance order, retaining order per group.
    let mut elements = Vec::new();
    for atom in atoms {
        if !elements.contains(&atom.atomic_number) {
            elements.push(atom.atomic_number);
        }
    }
    for &z in &elements {
        write!(w, "{} ", SYMBOLS.split_whitespace().nth(z - 1).unwrap())?;
    }
    writeln!(w)?;
    for &z in &elements {
        write!(
            w,
            "{} ",
            atoms.iter().filter(|a| a.atomic_number == z).count()
        )?;
    }
    writeln!(w, "\nDirect")?;
    for &z in &elements {
        for atom in atoms.iter().filter(|a| a.atomic_number == z) {
            let r = cell.cartesian_bohr_to_fractional(atom.position_bohr);
            writeln!(w, "{:.16e} {:.16e} {:.16e}", r[0], r[1], r[2])?;
        }
    }
    let d = data.shape();
    writeln!(w, "\n{} {} {}", d[0], d[1], d[2])?;
    let mut count = 0;
    for k in 0..d[2] {
        for j in 0..d[1] {
            for i in 0..d[0] {
                write!(w, " {:.16e}", data[[i, j, k]] * scale)?;
                count += 1;
                if count % 5 == 0 {
                    writeln!(w)?;
                }
            }
        }
    }
    writeln!(w)?;
    Ok(())
}

/// ASCII VTK ImageData, point samples with periodic endpoints (N+1 nodes).
/// r = Origin + Direction * diag(Spacing) * index, in Angstrom. Direction is
/// row-major, with columns a_i/|a_i| (not orthogonalized, even for skew cells).
/// Scalar values stay atomic units; no atom geometry in this image-only file.
pub fn write_vti(mut w: impl Write, cell: &Cell, data: &Array3<f64>, quantity: &str) -> Result<()> {
    validate(data, &[])?;
    let units = match quantity {
        "density" => "electrons_per_bohr3",
        "partial_charge" => "per_bohr3",
        "hartree" | "local_potential" => "hartree",
        "elf" => "dimensionless",
        _ => return Err("unknown VTI scalar quantity".into()),
    };
    let d = data.shape();
    let a = cell.lattice_columns_bohr();
    let lengths: [f64; 3] =
        std::array::from_fn(|axis| a[0][axis].hypot(a[1][axis]).hypot(a[2][axis]));
    let spacing: [f64; 3] =
        std::array::from_fn(|i| lengths[i] * dft_units::BOHR_IN_ANGSTROM / d[i] as f64);
    if lengths
        .iter()
        .chain(&spacing)
        .any(|x| !x.is_finite() || *x <= 0.0)
    {
        return Err("unrepresentable VTI spacing".into());
    }
    let extent = format!("0 {} 0 {} 0 {}", d[0], d[1], d[2]);
    writeln!(
        w,
        "<?xml version=\"1.0\"?>\n<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\">"
    )?;
    write!(
        w,
        "<ImageData WholeExtent=\"{extent}\" Origin=\"0 0 0\" Spacing=\"{:.16e} {:.16e} {:.16e}\" Direction=\"",
        spacing[0], spacing[1], spacing[2]
    )?;
    for row in a {
        for (value, length) in row.iter().zip(lengths) {
            write!(w, "{:.16e} ", value / length)?;
        }
    }
    let name = format!("{quantity}_{units}");
    writeln!(
        w,
        "\">\n<Piece Extent=\"{extent}\">\n<PointData Scalars=\"{name}\">\n<DataArray type=\"Float64\" Name=\"{name}\" NumberOfComponents=\"1\" format=\"ascii\">"
    )?;
    for k in 0..=d[2] {
        for j in 0..=d[1] {
            for i in 0..=d[0] {
                writeln!(w, "{:.16e}", data[[i % d[0], j % d[1], k % d[2]]])?;
            }
        }
    }
    writeln!(
        w,
        "</DataArray>\n</PointData>\n<CellData/>\n</Piece>\n</ImageData>\n</VTKFile>"
    )?;
    Ok(())
}
