//! SI1 (user decision 2026-11-27): compact `[structure]` block, native
//! structure-file readers (VASP POSCAR/CONTCAR, extended XYZ) and the
//! `[pseudopotentials] dir` lookup.
//!
//! Everything here rewrites the parsed TOML value into the explicit
//! `[system]` form (lattice rows, fractional or Cartesian atoms) before
//! deserialization, so the rest of the program, and the numbers it computes,
//! are those of the equivalent explicit input. No numerics live here.
//! Formats are implemented from their public descriptions (VASP wiki POSCAR
//! page; the extended XYZ specification of libAtoms/ASE); no source of other
//! codes is copied. Readers are validated against ASE 3.26 on committed
//! fixtures (`tests/fixtures/structure_input/`).
use super::{ConfigError, InputExpression};
use serde::{Deserialize, Serialize};
use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

fn invalid(message: impl Into<String>) -> ConfigError {
    ConfigError::Invalid(message.into())
}

/// How the atomic positions of a [`Structure`] are given.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Coordinates {
    /// In the basis of the lattice vectors.
    Fractional,
    /// Cartesian, in the structure's length unit.
    Cartesian,
}

/// One site with its optional per-atom flags.
#[derive(Clone, Debug, PartialEq)]
pub struct Site {
    pub species: String,
    pub position: [f64; 3],
    /// Cartesian components frozen during relaxation (`true` = frozen);
    /// `None` = free.
    pub fixed: Option<[bool; 3]>,
    /// Initial collinear moment in Bohr magnetons.
    pub magmom: Option<f64>,
}

/// A structure as read, before conversion to the `[system]` form.
#[derive(Clone, Debug)]
pub struct Structure {
    /// Rows are lattice vectors, in `unit`.
    pub lattice_rows: [[f64; 3]; 3],
    /// `"angstrom"` or `"bohr"`.
    pub unit: String,
    pub coordinates: Coordinates,
    pub sites: Vec<Site>,
}

/// Where the structure of a run came from (results.json `structure.source`).
#[derive(Clone, Debug, Default, Serialize)]
pub struct Source {
    /// `system` (explicit `[system]` table), `compact` (`[structure]` with
    /// `atoms`) or `file` (`[structure] file`).
    pub kind: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub file: Option<FileSource>,
    /// `fractional` or `cartesian` as given (compact and file forms).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub coordinates: Option<String>,
    /// Per-atom Cartesian relax constraints as given (empty: none given).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub fixed: Vec<Option<[bool; 3]>>,
    /// Per-atom initial moments in Bohr magnetons (empty: none given).
    #[serde(skip_serializing_if = "Vec::is_empty")]
    pub magmoms_mu_b: Vec<Option<f64>>,
    /// Sum of the moments written to `spin.initial_moment_seed` (None when
    /// no moment was nonzero).
    #[serde(skip_serializing_if = "Option::is_none")]
    pub initial_moment_seed_from_magmoms: Option<f64>,
}

#[derive(Clone, Debug, Serialize)]
pub struct FileSource {
    /// Path as written in the input.
    pub path: String,
    /// Path the file was read from (absolute when the input path is known).
    pub resolved_path: String,
    pub sha256: String,
    /// `poscar` or `extxyz`.
    pub format: String,
    /// Zero-based frame read (extended XYZ only) and the frame count.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub frame: Option<usize>,
    #[serde(skip_serializing_if = "Option::is_none")]
    pub frames: Option<usize>,
}

/// Result of [`expand`]: the record and the per-atom relax constraints.
#[derive(Debug, Default)]
pub struct Expanded {
    pub source: Source,
    pub fixed: Vec<crate::relax::Fixed>,
}

/// Record of the `[pseudopotentials] dir` lookup (results.json
/// `pseudopotential_lookup`).
#[derive(Clone, Debug, Serialize)]
pub struct PseudoLookup {
    pub dir: String,
    /// species -> path written into the pseudopotential table (relative paths
    /// resolve against the input file's directory, as every pseudopotential path).
    pub resolved: BTreeMap<String, String>,
    /// species -> `explicit` (per-species entry) or the matching rule.
    pub rule: BTreeMap<String, String>,
}

#[derive(Deserialize)]
#[serde(deny_unknown_fields)]
struct StructureTable {
    name: Option<String>,
    lattice_unit: Option<String>,
    lattice: Option<toml::Value>,
    coordinates: Option<String>,
    atoms: Option<String>,
    file: Option<String>,
    format: Option<String>,
    frame: Option<i64>,
    provenance: Option<toml::Value>,
}

fn parse_number(
    token: &str,
    path: &str,
    records: &mut Vec<InputExpression>,
) -> Result<f64, ConfigError> {
    if let Ok(v) = token.parse::<f64>() {
        if !v.is_finite() {
            return Err(invalid(format!("{path}: {token:?} is not finite")));
        }
        return Ok(v);
    }
    let v = dft_units::eval_expr(token).map_err(|e| invalid(format!("{path}: {e}")))?;
    if !v.is_finite() {
        return Err(invalid(format!("{path}: {token:?} is not finite")));
    }
    records.push(InputExpression {
        key_path: path.into(),
        text: token.into(),
        value: v,
    });
    Ok(v)
}

/// `fixed=` value: any nonempty subset of `x`, `y`, `z` (any order, no
/// repeats), `true` (= xyz) or `false` (= free).
pub fn parse_fixed(value: &str) -> Result<Option<[bool; 3]>, String> {
    match value {
        "true" => return Ok(Some([true; 3])),
        "false" => return Ok(None),
        _ => {}
    }
    let mut mask = [false; 3];
    if value.is_empty() {
        return Err("fixed: empty value; use x, y, z, xy, xz, yz, xyz, true or false".into());
    }
    for c in value.chars() {
        let axis = match c {
            'x' => 0,
            'y' => 1,
            'z' => 2,
            _ => {
                return Err(format!(
                    "fixed={value}: unknown component {c:?}; use x, y, z, xy, xz, yz, xyz, true or false"
                ));
            }
        };
        if mask[axis] {
            return Err(format!("fixed={value}: component {c} repeated"));
        }
        mask[axis] = true;
    }
    Ok(Some(mask))
}

fn check_species(species: &str, where_: &str) -> Result<(), String> {
    if species.is_empty()
        || !species.chars().next().unwrap().is_ascii_alphabetic()
        || species.contains('=')
    {
        return Err(format!(
            "{where_}: species label {species:?} must start with a letter (one label per line, then three coordinates)"
        ));
    }
    Ok(())
}

/// Parse the `atoms` string of `[structure]`: one atom per line,
/// `SPECIES c1 c2 c3 [fixed=...] [magmom=...]`; `#` starts a comment; blank
/// lines are skipped. Coordinates are numbers or expressions (`1/3`).
pub fn parse_atoms_block(
    text: &str,
    records: &mut Vec<InputExpression>,
) -> Result<Vec<Site>, ConfigError> {
    let mut sites = Vec::new();
    for (line_no, raw) in text.lines().enumerate() {
        let line = raw.split('#').next().unwrap().trim();
        if line.is_empty() {
            continue;
        }
        let where_ = format!("structure.atoms line {}", line_no + 1);
        let tokens: Vec<&str> = line.split_whitespace().collect();
        if tokens.len() < 4 {
            return Err(invalid(format!(
                "{where_}: expected `SPECIES c1 c2 c3 [fixed=...] [magmom=...]`, found {line:?}"
            )));
        }
        check_species(tokens[0], &where_).map_err(invalid)?;
        let mut position = [0.0; 3];
        for (j, token) in tokens[1..4].iter().enumerate() {
            if token.contains('=') {
                return Err(invalid(format!(
                    "{where_}: expected three coordinates before the flags, found {token:?}"
                )));
            }
            position[j] = parse_number(token, &format!("{where_} coordinate {}", j + 1), records)?;
        }
        let mut site = Site {
            species: tokens[0].into(),
            position,
            fixed: None,
            magmom: None,
        };
        let mut seen = Vec::new();
        for flag in &tokens[4..] {
            let Some((key, value)) = flag.split_once('=') else {
                return Err(invalid(format!(
                    "{where_}: unexpected token {flag:?}; flags are key=value (fixed=..., magmom=...)"
                )));
            };
            if seen.contains(&key) {
                return Err(invalid(format!("{where_}: flag {key} given twice")));
            }
            seen.push(key);
            match key {
                "fixed" => {
                    site.fixed =
                        parse_fixed(value).map_err(|e| invalid(format!("{where_}: {e}")))?
                }
                "magmom" => {
                    site.magmom = Some(parse_number(value, &format!("{where_} magmom"), records)?)
                }
                _ => {
                    return Err(invalid(format!(
                        "{where_}: unknown flag {key:?}; supported flags are fixed and magmom"
                    )));
                }
            }
        }
        sites.push(site);
    }
    if sites.is_empty() {
        return Err(invalid("structure.atoms: no atoms given"));
    }
    Ok(sites)
}

fn det(m: &[[f64; 3]; 3]) -> f64 {
    m[0][0] * (m[1][1] * m[2][2] - m[1][2] * m[2][1])
        - m[0][1] * (m[1][0] * m[2][2] - m[1][2] * m[2][0])
        + m[0][2] * (m[1][0] * m[2][1] - m[1][1] * m[2][0])
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

/// Convert constraints along lattice directions (VASP selective dynamics:
/// fractional coordinate i frozen) into Cartesian components. Freezing
/// fractional coordinate i restricts displacements d to b_i . d = 0 with b_i
/// the reciprocal direction a_j x a_k; this equals freezing Cartesian
/// component c exactly when b_i is parallel to the axis e_c. All three frozen
/// is always representable.
pub fn lattice_constraint_to_cartesian(
    rows: &[[f64; 3]; 3],
    frozen: [bool; 3],
) -> Result<Option<[bool; 3]>, String> {
    if frozen == [false; 3] {
        return Ok(None);
    }
    if frozen == [true; 3] {
        return Ok(Some([true; 3]));
    }
    let mut mask = [false; 3];
    for i in 0..3 {
        if !frozen[i] {
            continue;
        }
        let b = cross(rows[(i + 1) % 3], rows[(i + 2) % 3]);
        let norm = b.iter().map(|x| x * x).sum::<f64>().sqrt();
        let axis = (0..3).find(|&c| (0..3).all(|o| o == c || b[o].abs() <= 1e-12 * norm));
        match axis {
            Some(c) => mask[c] = true,
            None => {
                return Err(format!(
                    "selective dynamics freezes lattice direction {} of an atom, which is not a Cartesian axis in this cell; BAANDI relax constraints are Cartesian components; freeze all three directions of that atom or orient the cell so that the frozen reciprocal direction is a Cartesian axis",
                    i + 1
                ));
            }
        }
    }
    Ok(Some(mask))
}

/// VASP 5 POSCAR/CONTCAR (species line required). Lengths in Angstrom.
/// Scale: one positive factor, a negative number (target cell volume in
/// Angstrom^3), or three positive factors of the Cartesian x, y, z components.
/// Selective dynamics flags (T = free, F = frozen along that lattice
/// direction) become Cartesian `fixed` flags (see
/// [`lattice_constraint_to_cartesian`]). Velocities and predictor blocks after
/// the positions are ignored.
pub fn read_poscar(text: &str) -> Result<Structure, String> {
    let mut lines = text.lines().enumerate();
    let mut next = |what: &str| -> Result<(usize, Vec<String>), String> {
        lines
            .next()
            .map(|(i, l)| (i + 1, l.split_whitespace().map(str::to_owned).collect()))
            .ok_or_else(|| format!("POSCAR: file ends before the {what}"))
    };
    let float = |line: usize, t: &str| -> Result<f64, String> {
        t.parse::<f64>()
            .ok()
            .filter(|v| v.is_finite())
            .ok_or_else(|| format!("POSCAR line {line}: {t:?} is not a finite number"))
    };
    next("comment line")?;
    let (l, scale_tokens) = next("scale factor")?;
    let scale: Vec<f64> = match scale_tokens.len() {
        1 | 3 => scale_tokens
            .iter()
            .map(|t| float(l, t))
            .collect::<Result<_, _>>()?,
        n => {
            return Err(format!(
                "POSCAR line {l}: expected one scale factor or three, found {n} values"
            ));
        }
    };
    let mut rows = [[0.0; 3]; 3];
    for row in &mut rows {
        let (l, t) = next("lattice vectors")?;
        if t.len() < 3 {
            return Err(format!(
                "POSCAR line {l}: a lattice vector needs three numbers"
            ));
        }
        for c in 0..3 {
            row[c] = float(l, &t[c])?;
        }
    }
    let scale3 = if scale.len() == 3 {
        if scale.iter().any(|&s| s <= 0.0) {
            return Err("POSCAR line 2: three scale factors must all be positive".into());
        }
        [scale[0], scale[1], scale[2]]
    } else if scale[0] > 0.0 {
        [scale[0]; 3]
    } else if scale[0] < 0.0 {
        let v0 = det(&rows).abs();
        if v0 == 0.0 {
            return Err("POSCAR: the lattice vectors are linearly dependent".into());
        }
        [(-scale[0] / v0).cbrt(); 3]
    } else {
        return Err("POSCAR line 2: the scale factor must not be zero".into());
    };
    for row in &mut rows {
        for c in 0..3 {
            row[c] *= scale3[c];
        }
    }
    let (l, species_line) = next("species line")?;
    if species_line.is_empty() {
        return Err(format!("POSCAR line {l}: empty species line"));
    }
    if species_line[0].parse::<u64>().is_ok() {
        return Err(format!(
            "POSCAR line {l}: no species line (VASP 4 layout); add a line with the element symbols above the counts (VASP 5 layout)"
        ));
    }
    // VASP 6 may append "/hash" of the POTCAR to a species name.
    let species: Vec<String> = species_line
        .iter()
        .map(|s| s.split('/').next().unwrap().to_owned())
        .collect();
    for s in &species {
        check_species(s, &format!("POSCAR line {l}"))?;
    }
    let (l, count_tokens) = next("atom counts")?;
    if count_tokens.len() != species.len() {
        return Err(format!(
            "POSCAR line {l}: {} counts for {} species",
            count_tokens.len(),
            species.len()
        ));
    }
    let counts: Vec<usize> = count_tokens
        .iter()
        .map(|t| {
            t.parse::<usize>()
                .ok()
                .filter(|&n| n > 0)
                .ok_or_else(|| format!("POSCAR line {l}: {t:?} is not a positive atom count"))
        })
        .collect::<Result<_, _>>()?;
    let (mut l, mut mode) = next("coordinate mode line")?;
    let selective = mode
        .first()
        .is_some_and(|t| t.starts_with('s') || t.starts_with('S'));
    if selective {
        (l, mode) = next("coordinate mode line")?;
    }
    let cartesian = match mode.first().and_then(|t| t.chars().next()) {
        Some('c' | 'C' | 'k' | 'K') => true,
        Some('d' | 'D') => false,
        _ => {
            return Err(format!(
                "POSCAR line {l}: expected Direct or Cartesian, found {:?}",
                mode.join(" ")
            ));
        }
    };
    let mut sites = Vec::new();
    for (name, &count) in species.iter().zip(&counts) {
        for _ in 0..count {
            let (l, t) = next("atomic positions")?;
            let need = if selective { 6 } else { 3 };
            if t.len() < need {
                return Err(format!(
                    "POSCAR line {l}: expected {need} values (positions{})",
                    if selective {
                        " and three T/F flags"
                    } else {
                        ""
                    }
                ));
            }
            let mut p = [0.0; 3];
            for c in 0..3 {
                p[c] = float(l, &t[c])?;
                if cartesian {
                    p[c] *= scale3[c];
                }
            }
            let fixed = if selective {
                let mut frozen = [false; 3];
                for c in 0..3 {
                    frozen[c] = match t[3 + c].as_str() {
                        "T" | "t" => false,
                        "F" | "f" => true,
                        other => {
                            return Err(format!(
                                "POSCAR line {l}: selective dynamics flag {other:?} is not T or F"
                            ));
                        }
                    };
                }
                lattice_constraint_to_cartesian(&rows, frozen)
                    .map_err(|e| format!("POSCAR line {l}: {e}"))?
            } else {
                None
            };
            sites.push(Site {
                species: name.clone(),
                position: p,
                fixed,
                magmom: None,
            });
        }
    }
    Ok(Structure {
        lattice_rows: rows,
        unit: "angstrom".into(),
        coordinates: if cartesian {
            Coordinates::Cartesian
        } else {
            Coordinates::Fractional
        },
        sites,
    })
}

/// Split an extended XYZ comment line into key/value pairs. Values may be
/// bare words, double-quoted strings (backslash escapes `\"` and `\\`) or
/// `{...}` lists; a key without `=` has the value `T`. Keys are returned in
/// lower case (the specification treats them case-insensitively).
fn comment_pairs(line: &str) -> Result<Vec<(String, String)>, String> {
    let chars: Vec<char> = line.chars().collect();
    let mut i = 0;
    let mut pairs = Vec::new();
    let skip_ws = |i: &mut usize| {
        while *i < chars.len() && chars[*i].is_whitespace() {
            *i += 1;
        }
    };
    loop {
        skip_ws(&mut i);
        if i >= chars.len() {
            break;
        }
        let mut key = String::new();
        if chars[i] == '"' {
            i += 1;
            while i < chars.len() && chars[i] != '"' {
                key.push(chars[i]);
                i += 1;
            }
            i += 1;
        } else {
            while i < chars.len() && chars[i] != '=' && !chars[i].is_whitespace() {
                key.push(chars[i]);
                i += 1;
            }
        }
        skip_ws(&mut i);
        if i < chars.len() && chars[i] == '=' {
            i += 1;
            skip_ws(&mut i);
            let mut value = String::new();
            if i < chars.len() && chars[i] == '"' {
                i += 1;
                let mut closed = false;
                while i < chars.len() {
                    match chars[i] {
                        '\\' if i + 1 < chars.len() => {
                            value.push(chars[i + 1]);
                            i += 2;
                        }
                        '"' => {
                            closed = true;
                            i += 1;
                            break;
                        }
                        c => {
                            value.push(c);
                            i += 1;
                        }
                    }
                }
                if !closed {
                    return Err(format!("unterminated quoted value of key {key:?}"));
                }
            } else if i < chars.len() && chars[i] == '{' {
                i += 1;
                while i < chars.len() && chars[i] != '}' {
                    value.push(chars[i]);
                    i += 1;
                }
                if i >= chars.len() {
                    return Err(format!("unterminated {{...}} value of key {key:?}"));
                }
                i += 1;
            } else {
                while i < chars.len() && !chars[i].is_whitespace() {
                    value.push(chars[i]);
                    i += 1;
                }
            }
            pairs.push((key.to_lowercase(), value));
        } else {
            pairs.push((key.to_lowercase(), "T".into()));
        }
    }
    Ok(pairs)
}

fn logical(t: &str) -> Option<bool> {
    match t {
        "T" | "t" | "True" | "true" | "TRUE" => Some(true),
        "F" | "f" | "False" | "false" | "FALSE" => Some(false),
        _ => None,
    }
}

/// Extended XYZ (one frame): `Lattice="a1 a2 a3"` (nine numbers, rows are
/// lattice vectors, Angstrom), `Properties=name:type:ncols:...` (default
/// `species:S:1:pos:R:3`), optional `pbc` (default `T T T` with a lattice;
/// every direction must be periodic). Columns used: `species`, `pos`
/// (Cartesian Angstrom), optional `move_mask` (L:1 for the whole atom or L:3
/// per Cartesian component, `T` = free) and `initial_magmoms` or `magmoms`
/// (R:1, Bohr magnetons). Other columns are skipped. `frame` selects the frame
/// (zero-based, negative counts from the end; default the last).
/// Returns the structure, the zero-based frame index and the frame count.
pub fn read_extxyz(text: &str, frame: Option<i64>) -> Result<(Structure, usize, usize), String> {
    let lines: Vec<&str> = text.lines().collect();
    let mut starts = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        if lines[i].trim().is_empty() {
            i += 1;
            continue;
        }
        let n: usize = lines[i].trim().parse().map_err(|_| {
            format!(
                "extended XYZ line {}: expected the atom count, found {:?}",
                i + 1,
                lines[i]
            )
        })?;
        if n == 0 {
            return Err(format!("extended XYZ line {}: a frame has no atoms", i + 1));
        }
        // Count line, comment line, n atom lines.
        if i + 2 + n > lines.len() {
            return Err(format!(
                "extended XYZ line {}: the file ends inside a frame of {n} atoms",
                i + 1
            ));
        }
        starts.push((i, n));
        i += 2 + n;
    }
    if starts.is_empty() {
        return Err("extended XYZ: no frame in the file".into());
    }
    let count = starts.len();
    let wanted = frame.unwrap_or(-1);
    let index = if wanted < 0 {
        count as i64 + wanted
    } else {
        wanted
    };
    if index < 0 || index >= count as i64 {
        return Err(format!(
            "structure.frame = {wanted}: the file has {count} frame(s)"
        ));
    }
    let (start, n) = starts[index as usize];
    let comment_line = start + 2;
    let pairs = comment_pairs(lines[start + 1])
        .map_err(|e| format!("extended XYZ line {comment_line}: {e}"))?;
    let get = |key: &str| {
        pairs
            .iter()
            .find(|(k, _)| k == key)
            .map(|(_, v)| v.as_str())
    };
    let lattice = get("lattice").ok_or_else(|| {
        format!(
            "extended XYZ line {comment_line}: no Lattice=\"...\"; BAANDI needs a periodic cell"
        )
    })?;
    let numbers: Vec<f64> = lattice
        .split(|c: char| c.is_whitespace() || c == ',')
        .filter(|t| !t.is_empty())
        .map(|t| t.parse::<f64>().ok().filter(|v| v.is_finite()))
        .collect::<Option<_>>()
        .ok_or_else(|| format!("extended XYZ line {comment_line}: Lattice is not nine numbers"))?;
    if numbers.len() != 9 {
        return Err(format!(
            "extended XYZ line {comment_line}: Lattice has {} numbers, expected nine",
            numbers.len()
        ));
    }
    let rows: [[f64; 3]; 3] = std::array::from_fn(|r| std::array::from_fn(|c| numbers[3 * r + c]));
    if let Some(pbc) = get("pbc") {
        let flags: Vec<Option<bool>> = pbc
            .split(|c: char| c.is_whitespace() || c == ',')
            .filter(|t| !t.is_empty())
            .map(logical)
            .collect();
        if flags.len() != 3 || flags.iter().any(Option::is_none) {
            return Err(format!(
                "extended XYZ line {comment_line}: pbc must be three logical values"
            ));
        }
        if flags.contains(&Some(false)) {
            return Err(format!(
                "extended XYZ line {comment_line}: pbc=\"{pbc}\" is not periodic in every direction; BAANDI treats every cell as periodic: give the frame a full cell with enough vacuum and pbc=\"T T T\""
            ));
        }
    }
    let properties = get("properties").unwrap_or("species:S:1:pos:R:3");
    let parts: Vec<&str> = properties.split(':').collect();
    if !parts.len().is_multiple_of(3) {
        return Err(format!(
            "extended XYZ line {comment_line}: Properties must be name:type:count triples"
        ));
    }
    struct Column {
        name: String,
        kind: char,
        offset: usize,
        width: usize,
    }
    let mut columns = Vec::new();
    let mut offset = 0;
    for triple in parts.chunks(3) {
        let kind = match triple[1] {
            "S" | "R" | "I" | "L" => triple[1].chars().next().unwrap(),
            other => {
                return Err(format!(
                    "extended XYZ line {comment_line}: unknown property type {other:?}"
                ));
            }
        };
        let width: usize = triple[2].parse().ok().filter(|&w| w > 0).ok_or_else(|| {
            format!(
                "extended XYZ line {comment_line}: property {} has an invalid column count",
                triple[0]
            )
        })?;
        columns.push(Column {
            name: triple[0].to_lowercase(),
            kind,
            offset,
            width,
        });
        offset += width;
    }
    let find = |name: &str| columns.iter().find(|c| c.name == name);
    let species_col = find("species")
        .filter(|c| c.kind == 'S' && c.width == 1)
        .ok_or_else(|| format!("extended XYZ line {comment_line}: Properties needs species:S:1"))?;
    let pos_col = find("pos")
        .filter(|c| c.kind == 'R' && c.width == 3)
        .ok_or_else(|| format!("extended XYZ line {comment_line}: Properties needs pos:R:3"))?;
    let mask_col = find("move_mask");
    if let Some(c) = mask_col
        && (c.kind != 'L' || !(c.width == 1 || c.width == 3))
    {
        return Err(format!(
            "extended XYZ line {comment_line}: move_mask must be L:1 or L:3"
        ));
    }
    let magmom_col = match (find("initial_magmoms"), find("magmoms")) {
        (Some(_), Some(_)) => {
            return Err(format!(
                "extended XYZ line {comment_line}: both initial_magmoms and magmoms; keep one"
            ));
        }
        (a, b) => a.or(b),
    };
    if let Some(c) = magmom_col
        && (c.kind != 'R' || c.width != 1)
    {
        return Err(format!(
            "extended XYZ line {comment_line}: {} must be R:1 (collinear); non-collinear moments are not supported",
            c.name
        ));
    }
    let mut sites = Vec::with_capacity(n);
    for a in 0..n {
        let line_no = start + 3 + a;
        let t: Vec<&str> = lines[start + 2 + a].split_whitespace().collect();
        if t.len() != offset {
            return Err(format!(
                "extended XYZ line {line_no}: {} values, Properties declares {offset}",
                t.len()
            ));
        }
        let species = t[species_col.offset].to_owned();
        check_species(&species, &format!("extended XYZ line {line_no}"))?;
        let mut position = [0.0; 3];
        for c in 0..3 {
            let token = t[pos_col.offset + c];
            position[c] = token
                .parse::<f64>()
                .ok()
                .filter(|v| v.is_finite())
                .ok_or_else(|| {
                    format!("extended XYZ line {line_no}: {token:?} is not a finite number")
                })?;
        }
        let fixed = match mask_col {
            None => None,
            Some(c) => {
                let mut frozen = [false; 3];
                for k in 0..c.width {
                    let token = t[c.offset + k];
                    let movable = logical(token).ok_or_else(|| {
                        format!(
                            "extended XYZ line {line_no}: move_mask value {token:?} is not logical"
                        )
                    })?;
                    if c.width == 1 {
                        frozen = [!movable; 3];
                    } else {
                        frozen[k] = !movable;
                    }
                }
                (frozen != [false; 3]).then_some(frozen)
            }
        };
        let magmom = match magmom_col {
            None => None,
            Some(c) => {
                let token = t[c.offset];
                Some(
                    token
                        .parse::<f64>()
                        .ok()
                        .filter(|v| v.is_finite())
                        .ok_or_else(|| {
                            format!("extended XYZ line {line_no}: moment {token:?} is not a finite number")
                        })?,
                )
            }
        };
        sites.push(Site {
            species,
            position,
            fixed,
            magmom,
        });
    }
    Ok((
        Structure {
            lattice_rows: rows,
            unit: "angstrom".into(),
            coordinates: Coordinates::Cartesian,
            sites,
        },
        index as usize,
        count,
    ))
}

/// Format of a structure file from `format` or, when absent, its name:
/// `POSCAR*`, `CONTCAR*`, `*.vasp`, `*.poscar` are POSCAR; `*.extxyz`, `*.xyz`
/// extended XYZ (case-insensitive).
pub fn detect_format(path: &Path, format: Option<&str>) -> Result<&'static str, String> {
    match format {
        Some("poscar" | "vasp") => return Ok("poscar"),
        Some("extxyz") => return Ok("extxyz"),
        Some(other) => {
            return Err(format!(
                "structure.format = {other:?}: supported formats are poscar and extxyz (convert CIF and others with baandi-convert)"
            ));
        }
        None => {}
    }
    let name = path
        .file_name()
        .map(|n| n.to_string_lossy().to_lowercase())
        .unwrap_or_default();
    if name.starts_with("poscar")
        || name.starts_with("contcar")
        || name.ends_with(".vasp")
        || name.ends_with(".poscar")
    {
        Ok("poscar")
    } else if name.ends_with(".extxyz") || name.ends_with(".xyz") {
        Ok("extxyz")
    } else {
        Err(format!(
            "structure.file {:?}: cannot tell the format from the name; set structure.format = \"poscar\" or \"extxyz\" (CIF and other formats: baandi-convert)",
            path.display()
        ))
    }
}

/// Read a structure file (`poscar` or `extxyz`) and describe its source.
pub fn read_structure_file(
    given: &str,
    base: Option<&Path>,
    format: Option<&str>,
    frame: Option<i64>,
) -> Result<(Structure, FileSource), ConfigError> {
    let path = match base {
        Some(base) => base.join(given),
        None => PathBuf::from(given),
    };
    let format = detect_format(&path, format).map_err(invalid)?;
    if frame.is_some() && format != "extxyz" {
        return Err(invalid(
            "structure.frame applies to extended XYZ files only",
        ));
    }
    let bytes = std::fs::read(&path)
        .map_err(|e| invalid(format!("structure.file {}: {e}", path.display())))?;
    let text = std::str::from_utf8(&bytes).map_err(|e| {
        invalid(format!(
            "structure.file {}: not UTF-8 text ({e})",
            path.display()
        ))
    })?;
    let (structure, frame, frames) = match format {
        "poscar" => (read_poscar(text), None, None),
        _ => match read_extxyz(text, frame) {
            Ok((s, f, n)) => (Ok(s), Some(f), Some(n)),
            Err(e) => (Err(e), None, None),
        },
    };
    let structure =
        structure.map_err(|e| invalid(format!("structure.file {}: {e}", path.display())))?;
    let resolved = std::fs::canonicalize(&path).unwrap_or(path);
    Ok((
        structure,
        FileSource {
            path: given.into(),
            resolved_path: resolved.to_string_lossy().into_owned(),
            sha256: crate::run::sha256(&bytes),
            format: format.into(),
            frame,
            frames,
        },
    ))
}

fn formula(sites: &[Site]) -> String {
    let mut order: Vec<(&str, usize)> = Vec::new();
    for s in sites {
        match order.iter_mut().find(|(n, _)| *n == s.species) {
            Some((_, c)) => *c += 1,
            None => order.push((&s.species, 1)),
        }
    }
    order.iter().map(|(n, c)| format!("{n}{c}")).collect()
}

/// Rewrite `[structure]` into the explicit `[system]` table, map nonzero
/// per-atom moments onto `spin.initial_moment_seed` and collect the per-atom
/// relax constraints. Without `[structure]` the value is unchanged and the
/// source kind is `system`. `base` is the input file's directory (relative
/// structure files resolve against it; `None` = the working directory).
pub fn expand(
    value: &mut toml::Value,
    base: Option<&Path>,
    records: &mut Vec<InputExpression>,
) -> Result<Expanded, ConfigError> {
    let Some(root) = value.as_table_mut() else {
        return Ok(Expanded::default());
    };
    let Some(table) = root.remove("structure") else {
        return Ok(Expanded {
            source: Source {
                kind: "system".into(),
                ..Source::default()
            },
            fixed: vec![],
        });
    };
    if root.contains_key("system") {
        return Err(invalid(
            "[structure] and [system] both given: describe the structure once, either as [structure] (compact form or file) or as [system] with [[system.atoms]]",
        ));
    }
    let table: StructureTable = table
        .try_into()
        .map_err(|e: toml::de::Error| invalid(format!("[structure]: {e}")))?;
    let (structure, mut source) = if let Some(file) = &table.file {
        if table.lattice.is_some()
            || table.lattice_unit.is_some()
            || table.coordinates.is_some()
            || table.atoms.is_some()
        {
            return Err(invalid(
                "structure.file conflicts with lattice, lattice_unit, coordinates and atoms: the file supplies them; remove one of the two descriptions",
            ));
        }
        let (structure, file_source) =
            read_structure_file(file, base, table.format.as_deref(), table.frame)?;
        let coordinates = match structure.coordinates {
            Coordinates::Fractional => "fractional",
            Coordinates::Cartesian => "cartesian",
        };
        (
            structure,
            Source {
                kind: "file".into(),
                file: Some(file_source),
                coordinates: Some(coordinates.into()),
                ..Source::default()
            },
        )
    } else {
        if table.format.is_some() || table.frame.is_some() {
            return Err(invalid(
                "structure.format and structure.frame require structure.file",
            ));
        }
        let missing = |key: &str| {
            invalid(format!(
                "structure.{key} is required (compact form: lattice_unit, lattice, coordinates and atoms; or give structure.file)"
            ))
        };
        let unit = table
            .lattice_unit
            .clone()
            .ok_or_else(|| missing("lattice_unit"))?;
        if !matches!(unit.as_str(), "angstrom" | "bohr") {
            return Err(invalid(format!(
                "structure.lattice_unit = {unit:?}: use \"angstrom\" or \"bohr\""
            )));
        }
        let mut lattice = table.lattice.clone().ok_or_else(|| missing("lattice"))?;
        if let Some(rows) = lattice.as_array_mut() {
            for (r, row) in rows.iter_mut().enumerate() {
                if let Some(row) = row.as_array_mut() {
                    for (c, x) in row.iter_mut().enumerate() {
                        if let Some(text) = x.as_str() {
                            let v = parse_number(
                                text,
                                &format!("structure.lattice[{r}][{c}]"),
                                records,
                            )?;
                            *x = toml::Value::Float(v);
                        }
                    }
                }
            }
        }
        let rows: [[f64; 3]; 3] = lattice.try_into().map_err(|e: toml::de::Error| {
            invalid(format!(
                "structure.lattice: three rows of three numbers required ({e})"
            ))
        })?;
        let coordinates = match table.coordinates.as_deref() {
            Some("fractional") => Coordinates::Fractional,
            Some("cartesian") => Coordinates::Cartesian,
            Some(other) => {
                return Err(invalid(format!(
                    "structure.coordinates = {other:?}: use \"fractional\" or \"cartesian\" (Cartesian in lattice_unit)"
                )));
            }
            None => return Err(missing("coordinates")),
        };
        let atoms = table.atoms.as_deref().ok_or_else(|| missing("atoms"))?;
        let sites = parse_atoms_block(atoms, records)?;
        (
            Structure {
                lattice_rows: rows,
                unit,
                coordinates,
                sites,
            },
            Source {
                kind: "compact".into(),
                coordinates: table.coordinates.clone(),
                ..Source::default()
            },
        )
    };
    // [system] in the explicit form; Cartesian positions are converted by the
    // same code as explicit `cartesian` atoms (config::normalize).
    let mut system = toml::map::Map::new();
    system.insert(
        "name".into(),
        toml::Value::String(
            table
                .name
                .clone()
                .unwrap_or_else(|| formula(&structure.sites)),
        ),
    );
    system.insert(
        "lattice_unit".into(),
        toml::Value::String(structure.unit.clone()),
    );
    system.insert(
        "lattice_vectors".into(),
        toml::Value::Array(
            structure
                .lattice_rows
                .iter()
                .map(|r| toml::Value::Array(r.iter().map(|&x| toml::Value::Float(x)).collect()))
                .collect(),
        ),
    );
    let key = match structure.coordinates {
        Coordinates::Fractional => "fractional",
        Coordinates::Cartesian => "cartesian",
    };
    system.insert(
        "atoms".into(),
        toml::Value::Array(
            structure
                .sites
                .iter()
                .map(|s| {
                    let mut atom = toml::map::Map::new();
                    atom.insert("species".into(), toml::Value::String(s.species.clone()));
                    atom.insert(
                        key.into(),
                        toml::Value::Array(
                            s.position.iter().map(|&x| toml::Value::Float(x)).collect(),
                        ),
                    );
                    toml::Value::Table(atom)
                })
                .collect(),
        ),
    );
    if let Some(p) = table.provenance {
        if !p.is_table() {
            return Err(invalid("structure.provenance must be a table"));
        }
        system.insert("provenance".into(), p);
    }
    root.insert("system".into(), toml::Value::Table(system));

    // Relax constraints (Cartesian components, as [relax] fixed).
    let fixed: Vec<crate::relax::Fixed> = structure
        .sites
        .iter()
        .enumerate()
        .filter_map(|(atom, s)| {
            s.fixed
                .map(|components| crate::relax::Fixed { atom, components })
        })
        .collect();
    if !fixed.is_empty() && root.get("relax").and_then(|r| r.get("fixed")).is_some() {
        return Err(invalid(
            "relax.fixed and per-atom fixed flags of [structure] both given: give the constraints once (remove relax.fixed or the fixed flags / selective dynamics)",
        ));
    }
    if structure.sites.iter().any(|s| s.fixed.is_some()) {
        source.fixed = structure.sites.iter().map(|s| s.fixed).collect();
    }

    // Initial collinear moments: the engine starts from a uniform
    // polarization, so only the sum enters (spin.initial_moment_seed).
    if structure.sites.iter().any(|s| s.magmom.is_some()) {
        source.magmoms_mu_b = structure.sites.iter().map(|s| s.magmom).collect();
    }
    if structure
        .sites
        .iter()
        .any(|s| s.magmom.is_some_and(|m| m != 0.0))
    {
        let seed: f64 = structure
            .sites
            .iter()
            .map(|s| s.magmom.unwrap_or(0.0))
            .sum();
        let spin = root.get_mut("spin").and_then(toml::Value::as_table_mut);
        let Some(spin) =
            spin.filter(|s| s.get("mode").and_then(toml::Value::as_str) == Some("collinear"))
        else {
            return Err(invalid(
                "magmom: nonzero initial moments need [spin] mode = \"collinear\" (with its electron sector); remove the moments for an unpolarized run",
            ));
        };
        if spin.contains_key("initial_moment_seed") {
            return Err(invalid(
                "spin.initial_moment_seed and per-atom magmom both given: give the initial moment once",
            ));
        }
        spin.insert("initial_moment_seed".into(), toml::Value::Float(seed));
        source.initial_moment_seed_from_magmoms = Some(seed);
    }
    Ok(Expanded { source, fixed })
}

/// Files of `dir` whose names satisfy `matches`, sorted.
fn matching(names: &[String], matches: impl Fn(&str) -> bool) -> Vec<String> {
    names.iter().filter(|n| matches(n)).cloned().collect()
}

/// Lookup of one stem in the directory listing; `Ok(None)` = no match.
fn lookup_stem(names: &[String], stem: &str) -> Result<Option<(String, &'static str)>, String> {
    for exact in [format!("{stem}.upf"), format!("{stem}.UPF")] {
        if names.contains(&exact) {
            return Ok(Some((exact, "exact")));
        }
    }
    let lower = stem.to_lowercase();
    let folded = matching(names, |n| n.to_lowercase() == format!("{lower}.upf"));
    match folded.len() {
        0 => {}
        1 => return Ok(Some((folded[0].clone(), "case_insensitive"))),
        _ => {
            return Err(format!(
                "several files match {stem}.upf ignoring case: {folded:?}"
            ));
        }
    }
    let prefixed = matching(names, |n| {
        let n = n.to_lowercase();
        n.ends_with(".upf")
            && [".", "_", "-"]
                .iter()
                .any(|sep| n.starts_with(&format!("{lower}{sep}")))
    });
    match prefixed.len() {
        0 => Ok(None),
        1 => Ok(Some((prefixed[0].clone(), "prefix"))),
        _ => Err(format!(
            "several files start with {stem}. / {stem}_ / {stem}- : {prefixed:?}"
        )),
    }
}

/// Element symbol at the start of a species label (`Fe1`, `Fe_up` -> `Fe`).
fn element_prefix(label: &str) -> Option<String> {
    let mut chars = label.chars();
    let first = chars.next().filter(char::is_ascii_uppercase)?;
    let mut symbol = first.to_string();
    if let Some(c) = chars.next().filter(char::is_ascii_lowercase) {
        symbol.push(c);
    }
    Some(symbol)
}

/// `[pseudopotentials] dir = "..."`: fill the pseudopotential entry of every
/// species used by `system.atoms` that has no explicit entry. Rules, in
/// order, first for the species label and then for its leading element
/// symbol: exact `<stem>.upf`, exact `<stem>.UPF`, the unique file named
/// `<stem>.upf` ignoring case, the unique `*.upf` file (ignoring case) whose
/// name starts with `<stem>.`, `<stem>_` or `<stem>-` ignoring case. Several
/// candidates at one rule are an error; so is no match.
pub fn resolve_pseudo_dir(
    value: &mut toml::Value,
    base: Option<&Path>,
) -> Result<Option<PseudoLookup>, ConfigError> {
    let species: Vec<String> = value
        .get("system")
        .and_then(|s| s.get("atoms"))
        .and_then(toml::Value::as_array)
        .map(|atoms| {
            let mut v: Vec<String> = Vec::new();
            for a in atoms {
                if let Some(s) = a.get("species").and_then(toml::Value::as_str)
                    && !v.iter().any(|x| x == s)
                {
                    v.push(s.to_owned());
                }
            }
            v
        })
        .unwrap_or_default();
    let Some(table) = value
        .get_mut("pseudopotentials")
        .and_then(toml::Value::as_table_mut)
    else {
        return Ok(None);
    };
    let Some(dir) = table.remove("dir") else {
        return Ok(None);
    };
    let dir = dir
        .as_str()
        .map(str::to_owned)
        .ok_or_else(|| invalid("pseudopotentials.dir must be a string (a directory path)"))?;
    let on_disk = match base {
        Some(base) => base.join(&dir),
        None => PathBuf::from(&dir),
    };
    let names: Vec<String> = std::fs::read_dir(&on_disk)
        .map_err(|e| invalid(format!("pseudopotentials.dir {}: {e}", on_disk.display())))?
        .filter_map(|e| e.ok())
        .filter(|e| e.file_type().is_ok_and(|t| t.is_file() || t.is_symlink()))
        .map(|e| e.file_name().to_string_lossy().into_owned())
        .collect::<Vec<_>>();
    let mut names = names;
    names.sort();
    let mut lookup = PseudoLookup {
        dir: dir.clone(),
        resolved: BTreeMap::new(),
        rule: BTreeMap::new(),
    };
    for s in &species {
        if table.contains_key(s) {
            lookup.rule.insert(s.clone(), "explicit".into());
            continue;
        }
        let mut stems = vec![s.clone()];
        if let Some(e) = element_prefix(s)
            && e != *s
        {
            stems.push(e);
        }
        let mut found = None;
        for stem in &stems {
            if let Some((file, rule)) = lookup_stem(&names, stem)
                .map_err(|e| invalid(format!("pseudopotentials.dir {dir:?}, species {s:?}: {e}; add an explicit entry {s} = \"...\"")))?
            {
                found = Some((file, format!("{rule}:{stem}")));
                break;
            }
        }
        let Some((file, rule)) = found else {
            return Err(invalid(format!(
                "pseudopotentials.dir {dir:?}: no pseudopotential for species {s:?} (tried {} .upf/.UPF, ignoring case, and prefixed names); add an explicit entry {s} = \"...\"",
                stems.join(" and ")
            )));
        };
        let path = Path::new(&dir).join(&file).to_string_lossy().into_owned();
        table.insert(s.clone(), toml::Value::String(path.clone()));
        lookup.resolved.insert(s.clone(), path);
        lookup.rule.insert(s.clone(), rule);
    }
    Ok(Some(lookup))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn atoms(text: &str) -> Result<Vec<Site>, String> {
        parse_atoms_block(text, &mut Vec::new()).map_err(|e| e.to_string())
    }

    #[test]
    fn fixed_flag_values() {
        let all = [
            ("x", [true, false, false]),
            ("y", [false, true, false]),
            ("z", [false, false, true]),
            ("xy", [true, true, false]),
            ("yx", [true, true, false]),
            ("xz", [true, false, true]),
            ("yz", [false, true, true]),
            ("xyz", [true; 3]),
            ("zyx", [true; 3]),
            ("true", [true; 3]),
        ];
        for (v, mask) in all {
            assert_eq!(parse_fixed(v).unwrap(), Some(mask), "{v}");
        }
        assert_eq!(parse_fixed("false").unwrap(), None);
        for bad in ["", "xx", "w", "X", "x y", "1", "TRUE"] {
            assert!(parse_fixed(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn atoms_block_lines_flags_and_expressions() {
        let mut records = Vec::new();
        let sites = parse_atoms_block(
            "\n# comment\nSi 0 0 0\nSi 1/4 0.25 2.5e-1 fixed=xz magmom=-0.5  # tail\n\nFe1 0.5 0.5 0.5 magmom=2 fixed=true\nO 0.1 0.2 0.3 fixed=false\n",
            &mut records,
        )
        .unwrap();
        assert_eq!(sites.len(), 4);
        assert_eq!(sites[0].fixed, None);
        assert_eq!(sites[0].magmom, None);
        assert_eq!(sites[1].position, [0.25; 3]);
        assert_eq!(sites[1].fixed, Some([true, false, true]));
        assert_eq!(sites[1].magmom, Some(-0.5));
        assert_eq!(sites[2].species, "Fe1");
        assert_eq!(sites[2].fixed, Some([true; 3]));
        assert_eq!(sites[2].magmom, Some(2.0));
        assert_eq!(sites[3].fixed, None);
        // Only the expression is recorded; plain numbers keep their bits.
        assert_eq!(records.len(), 1);
        assert_eq!(records[0].text, "1/4");
        assert!(records[0].key_path.contains("line 4 coordinate 1"));
        assert_eq!(
            atoms("Si 0.1 0.2 0.3").unwrap()[0].position[0].to_bits(),
            0.1f64.to_bits()
        );
    }

    #[test]
    fn atoms_block_errors() {
        for (text, needle) in [
            ("", "no atoms"),
            ("# only a comment", "no atoms"),
            ("Si 0 0", "expected `SPECIES"),
            ("Si 0 0 fixed=x", "three coordinates before the flags"),
            ("Si 0 0 x=1 2", "three coordinates"),
            ("1 0 0 0", "must start with a letter"),
            ("Si 0 0 abc", "coordinate 3"),
            ("Si 0 0 inf", "not finite"),
            ("Si 0 0 0 fixed", "flags are key=value"),
            ("Si 0 0 0 fixed=x fixed=y", "given twice"),
            ("Si 0 0 0 magmom=a", "magmom"),
            ("Si 0 0 0 label=a", "unknown flag"),
        ] {
            let err = atoms(text).unwrap_err();
            assert!(err.contains(needle), "{text:?}: {err}");
        }
    }

    const SI: &str = "Si\n1.0\n0 2.715 2.715\n2.715 0 2.715\n2.715 2.715 0\nSi\n2\nDirect\n0 0 0\n0.25 0.25 0.25\n";

    #[test]
    fn poscar_layouts() {
        let s = read_poscar(SI).unwrap();
        assert_eq!(s.coordinates, Coordinates::Fractional);
        assert_eq!(s.sites[1].position, [0.25; 3]);
        assert_eq!(s.lattice_rows[0], [0.0, 2.715, 2.715]);
        // Scale factor, Cartesian (also scaled), VASP 6 species hash, comments.
        let c = SI
            .replace("1.0\n", "2.0\n")
            .replace("Si\n2\n", "Si/abc123\n2\n")
            .replace(
                "Direct\n0 0 0\n0.25 0.25 0.25\n",
                "cart\n0 0 0 ! one\n0.5 0.5 0.5 Si\n",
            );
        let s = read_poscar(&c).unwrap();
        assert_eq!(s.coordinates, Coordinates::Cartesian);
        assert_eq!(s.sites[0].species, "Si");
        assert_eq!(s.sites[1].position, [1.0; 3]);
        assert_eq!(s.lattice_rows[0], [0.0, 5.43, 5.43]);
        // Negative scale = volume.
        let v = SI.replace("1.0\n", "-80.0\n");
        let s = read_poscar(&v).unwrap();
        assert!((det(&s.lattice_rows).abs() - 80.0).abs() < 1e-12);
        // Selective dynamics in a cubic cell: lattice direction = Cartesian axis.
        let sd = "x\n1\n4 0 0\n0 5 0\n0 0 6\nH O\n1 1\nSel\nD\n0 0 0 F F F\n0.5 0.5 0.5 T F T\n";
        let s = read_poscar(sd).unwrap();
        assert_eq!(s.sites[0].fixed, Some([true; 3]));
        assert_eq!(s.sites[1].fixed, Some([false, true, false]));
        assert_eq!(s.sites[1].species, "O");
    }

    #[test]
    fn poscar_errors() {
        for (text, needle) in [
            (SI.replace("Si\n2\n", "2\n"), "VASP 4"),
            (SI.replace("1.0\n", "0\n"), "must not be zero"),
            (SI.replace("1.0\n", "1 2\n"), "one scale factor or three"),
            (SI.replace("1.0\n", "1 -1 1\n"), "all be positive"),
            (SI.replace("Si\n2\n", "Si\n0\n"), "positive atom count"),
            (SI.replace("Si\n2\n", "Si O\n2\n"), "1 counts for 2 species"),
            (SI.replace("Direct", "Fractional"), "Direct or Cartesian"),
            (SI.replace("0.25 0.25 0.25\n", ""), "file ends before"),
            (SI.replace("0 2.715 2.715\n", "0 2.715\n"), "three numbers"),
            (SI.replace("Direct\n", "Selective\nDirect\n"), "T/F flags"),
            (
                SI.replace("Direct\n0 0 0\n", "S\nDirect\n0 0 0 T X T\n")
                    .replace("0.25 0.25 0.25\n", "0.25 0.25 0.25 T T T\n"),
                "not T or F",
            ),
            // fcc cell: freezing one lattice direction is not a Cartesian axis.
            (
                SI.replace("Direct\n0 0 0\n", "S\nDirect\n0 0 0 F T T\n")
                    .replace("0.25 0.25 0.25\n", "0.25 0.25 0.25 T T T\n"),
                "not a Cartesian axis",
            ),
        ] {
            let err = read_poscar(&text).unwrap_err();
            assert!(err.contains(needle), "{needle}: {err}");
        }
    }

    #[test]
    fn lattice_constraints_to_cartesian() {
        let hex = [
            [3.0, 0.0, 0.0],
            [-1.5, 3.0 * 0.75f64.sqrt(), 0.0],
            [0.0, 0.0, 5.0],
        ];
        // c frozen: the reciprocal direction of a3 is z.
        assert_eq!(
            lattice_constraint_to_cartesian(&hex, [false, false, true]).unwrap(),
            Some([false, false, true])
        );
        assert!(lattice_constraint_to_cartesian(&hex, [true, false, false]).is_err());
        assert_eq!(
            lattice_constraint_to_cartesian(&hex, [true; 3]).unwrap(),
            Some([true; 3])
        );
        assert_eq!(
            lattice_constraint_to_cartesian(&hex, [false; 3]).unwrap(),
            None
        );
        // Monoclinic: b axis along y, a and c in xz; a2 x a3 is not along x.
        let mono = [[4.0, 0.0, 0.0], [0.0, 5.0, 0.0], [1.0, 0.0, 6.0]];
        assert_eq!(
            lattice_constraint_to_cartesian(&mono, [false, true, false]).unwrap(),
            Some([false, true, false])
        );
        assert_eq!(
            lattice_constraint_to_cartesian(&mono, [false, false, true]).unwrap(),
            Some([false, false, true])
        );
        assert!(lattice_constraint_to_cartesian(&mono, [true, false, false]).is_err());
    }

    const XYZ: &str = "2\nLattice=\"5 0 0 0 5 0 0 0 5\" Properties=species:S:1:pos:R:3 pbc=\"T T T\"\nH 0 0 0\nH 0.7 0 0\n";

    #[test]
    fn extxyz_frames_columns_and_flags() {
        let (s, frame, frames) = read_extxyz(XYZ, None).unwrap();
        assert_eq!((frame, frames), (0, 1));
        assert_eq!(s.coordinates, Coordinates::Cartesian);
        assert_eq!(s.sites[1].position, [0.7, 0.0, 0.0]);
        // Without Properties the default columns apply; without pbc periodic.
        let bare = "1\nLattice=\"5 0 0 0 5 0 0 0 5\"\nO 1 2 3\n";
        assert_eq!(
            read_extxyz(bare, None).unwrap().0.sites[0].position,
            [1.0, 2.0, 3.0]
        );
        // Several frames, blank lines between frames; negative and positive index.
        let multi = format!("{XYZ}\n{}", XYZ.replace("0.7", "0.8"));
        assert_eq!(
            read_extxyz(&multi, None).unwrap().0.sites[1].position[0],
            0.8
        );
        assert_eq!(
            read_extxyz(&multi, Some(0)).unwrap().0.sites[1].position[0],
            0.7
        );
        assert_eq!(read_extxyz(&multi, Some(-2)).unwrap().1, 0);
        // move_mask L:1 and L:3, magmoms, extra columns, quoted info.
        let flags = "2\nenergy=-1 info=\"a \\\"b\\\" c\" Lattice=\"5 0 0 0 5 0 0 0 5\" Properties=species:S:1:pos:R:3:forces:R:3:move_mask:L:3:magmoms:R:1\nFe 0 0 0 0 0 0 F T F 2.0\nFe 1 1 1 0 0 0 T T T -2.0\n";
        let s = read_extxyz(flags, None).unwrap().0;
        assert_eq!(s.sites[0].fixed, Some([true, false, true]));
        assert_eq!(s.sites[1].fixed, None);
        assert_eq!(s.sites[0].magmom, Some(2.0));
        let l1 = XYZ
            .replace("pos:R:3", "pos:R:3:move_mask:L:1")
            .replace("H 0 0 0", "H 0 0 0 F")
            .replace("H 0.7 0 0", "H 0.7 0 0 True");
        let s = read_extxyz(&l1, None).unwrap().0;
        assert_eq!(s.sites[0].fixed, Some([true; 3]));
        assert_eq!(s.sites[1].fixed, None);
    }

    #[test]
    fn extxyz_errors() {
        for (text, needle) in [
            (String::new(), "no frame"),
            (XYZ.replace("2\n", "x\n"), "atom count"),
            (XYZ.replace("H 0.7 0 0\n", ""), "ends inside a frame"),
            (
                XYZ.replace("Lattice=\"5 0 0 0 5 0 0 0 5\" ", ""),
                "no Lattice",
            ),
            (
                XYZ.replace("5 0 0 0 5 0 0 0 5", "5 0 0 0 5 0 0 0"),
                "expected nine",
            ),
            (XYZ.replace("T T T", "T F T"), "not periodic"),
            (XYZ.replace("T T T", "T T"), "three logical"),
            (XYZ.replace("pos:R:3", "pos:R:2"), "pos:R:3"),
            (XYZ.replace("species:S:1:", ""), "species:S:1"),
            (
                XYZ.replace("pos:R:3", "pos:R:3:magmoms:R:3"),
                "non-collinear",
            ),
            (
                XYZ.replace("pos:R:3", "pos:R:3:move_mask:L:2"),
                "L:1 or L:3",
            ),
            (
                XYZ.replace("pos:R:3", "pos:R:3:magmoms:R:1:initial_magmoms:R:1"),
                "keep one",
            ),
            (XYZ.replace("pos:R:3", "pos:X:3"), "unknown property type"),
            (XYZ.replace("H 0 0 0", "H 0 0 0 1"), "Properties declares 4"),
            (XYZ.replace("H 0 0 0", "H 0 0 nan"), "finite"),
            (XYZ.replace("Lattice=\"5", "Lattice=\"x"), "nine numbers"),
            (XYZ.replace("pbc=\"T T T\"", "pbc=\"T T T"), "unterminated"),
        ] {
            let err = read_extxyz(&text, None).unwrap_err();
            assert!(err.contains(needle), "{needle}: {err}");
        }
        assert!(
            read_extxyz(XYZ, Some(1))
                .unwrap_err()
                .contains("1 frame(s)")
        );
        assert!(read_extxyz(XYZ, Some(-2)).is_err());
    }

    #[test]
    fn format_detection() {
        for (name, f) in [
            ("POSCAR", "poscar"),
            ("CONTCAR", "poscar"),
            ("poscar_relaxed", "poscar"),
            ("cell.vasp", "poscar"),
            ("x.POSCAR", "poscar"),
            ("a.extxyz", "extxyz"),
            ("a.XYZ", "extxyz"),
        ] {
            assert_eq!(detect_format(Path::new(name), None).unwrap(), f, "{name}");
        }
        assert!(
            detect_format(Path::new("a.cif"), None)
                .unwrap_err()
                .contains("baandi-convert")
        );
        assert_eq!(
            detect_format(Path::new("a.cif"), Some("extxyz")).unwrap(),
            "extxyz"
        );
        assert!(detect_format(Path::new("POSCAR"), Some("cif")).is_err());
    }

    #[test]
    fn formula_and_element_prefix() {
        let site = |s: &str| Site {
            species: s.into(),
            position: [0.0; 3],
            fixed: None,
            magmom: None,
        };
        assert_eq!(formula(&[site("Si"), site("O"), site("Si")]), "Si2O1");
        assert_eq!(element_prefix("Fe1").as_deref(), Some("Fe"));
        assert_eq!(element_prefix("O_pbe").as_deref(), Some("O"));
        assert_eq!(element_prefix("fe"), None);
    }
}
