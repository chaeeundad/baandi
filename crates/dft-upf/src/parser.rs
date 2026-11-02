use dft_units::{Hartree, Rydberg};
use std::{collections::BTreeMap, path::Path};
use thiserror::Error;
use upf::model::PpBetaNode;

#[derive(Debug, Error)]
pub enum UpfError {
    #[error("UPF I/O: {0}")]
    Io(#[from] std::io::Error),
    #[error("malformed UPF: {0}")]
    Parse(String),
    #[error("unsupported UPF version: {0}")]
    Version(String),
    #[error("unsupported feature: {0}")]
    Unsupported(String),
    #[error("missing field: {0}")]
    Missing(String),
    #[error("invalid field {field}: {reason}")]
    Invalid { field: String, reason: String },
    #[error("PP_DIJ is not symmetric within 1e-12 Ha")]
    NonSymmetricD,
}
fn invalid(field: &str, reason: &str) -> UpfError {
    UpfError::Invalid {
        field: field.into(),
        reason: reason.into(),
    }
}

/// Immutable after validation. The stored function is u(r)=r*beta(r), Bohr^-1/2.
#[derive(Debug)]
pub struct Projector {
    pub(crate) l: usize,
    pub(crate) radial_times_r: Vec<f64>,
    pub(crate) cutoff: usize,
}
impl Projector {
    pub fn l(&self) -> usize {
        self.l
    }
    pub fn radial_times_r(&self) -> &[f64] {
        &self.radial_times_r
    }
    /// Number of valid points, equal to the one-based final UPF index.
    pub fn cutoff_index(&self) -> usize {
        self.cutoff
    }
}
#[derive(Debug)]
pub struct Pseudopotential {
    element: String,
    z_valence: f64,
    functional: String,
    pub(crate) mesh: crate::radial::RadialMesh,
    pub(crate) local: Vec<f64>,
    projectors: Vec<Projector>,
    dij: Vec<f64>,
    raw: upf::UpfData,
}
impl Pseudopotential {
    /// Original reader model, in UPF units; never use directly in numerical kernels.
    pub fn raw(&self) -> &upf::UpfData {
        &self.raw
    }
    pub fn element(&self) -> &str {
        &self.element
    }
    pub fn z_valence(&self) -> f64 {
        self.z_valence
    }
    pub fn functional(&self) -> &str {
        &self.functional
    }
    pub fn mesh(&self) -> &crate::radial::RadialMesh {
        &self.mesh
    }
    pub fn local_hartree(&self) -> &[f64] {
        &self.local
    }
    pub fn projectors(&self) -> &[Projector] {
        &self.projectors
    }
    /// Dense row-major real matrix, in Ha.
    pub fn dij_hartree(&self) -> &[f64] {
        &self.dij
    }
}

// Supplemental XML metadata only: the upf crate reads all numerical arrays.
// Check balanced markup so truncated documents cannot be accepted by serde.
#[derive(Debug)]
struct Tag {
    name: String,
    attrs: BTreeMap<String, String>,
}
fn decode_attribute(mut text: &str) -> Result<String, UpfError> {
    let mut result = String::new();
    while let Some(i) = text.find('&') {
        result.push_str(&text[..i]);
        text = &text[i + 1..];
        let end = text
            .find(';')
            .ok_or_else(|| UpfError::Parse("unterminated XML entity".into()))?;
        let entity = &text[..end];
        let c = match entity {
            "amp" => '&',
            "lt" => '<',
            "gt" => '>',
            "quot" => '"',
            "apos" => '\'',
            _ => {
                let code = if let Some(hex) = entity.strip_prefix("#x") {
                    u32::from_str_radix(hex, 16).ok()
                } else {
                    entity.strip_prefix('#').and_then(|s| s.parse::<u32>().ok())
                };
                code.and_then(char::from_u32)
                    .ok_or_else(|| UpfError::Parse("invalid XML entity".into()))?
            }
        };
        result.push(c);
        text = &text[end + 1..];
    }
    result.push_str(text);
    Ok(result)
}
fn tags(text: &str) -> Result<Vec<Tag>, UpfError> {
    let mut rest = text;
    let mut stack: Vec<String> = Vec::new();
    let mut out = Vec::new();
    let mut closed_root = false;
    while let Some(pos) = rest.find('<') {
        if stack.is_empty() && !rest[..pos].trim().is_empty() {
            return Err(UpfError::Parse("text outside root".into()));
        }
        rest = &rest[pos..];
        let special = if rest.starts_with("<!--") {
            Some(("-->", 4))
        } else if rest.starts_with("<?") {
            Some(("?>", 2))
        } else if rest.starts_with("<![CDATA[") && !stack.is_empty() {
            Some(("]]>", 9))
        } else {
            None
        };
        if let Some((end, start)) = special {
            let i = rest[start..]
                .find(end)
                .ok_or_else(|| UpfError::Parse("unterminated XML metadata".into()))?
                + start;
            rest = &rest[i + end.len()..];
            continue;
        }
        let mut quote = None;
        let end = rest
            .char_indices()
            .skip(1)
            .find_map(|(i, c)| {
                if let Some(q) = quote {
                    if c == q {
                        quote = None;
                    }
                } else if c == '\'' || c == '"' {
                    quote = Some(c);
                } else if c == '>' {
                    return Some(i);
                }
                None
            })
            .ok_or_else(|| UpfError::Parse("unterminated tag".into()))?;
        let body = rest[1..end].trim();
        rest = &rest[end + 1..];
        if let Some(name) = body.strip_prefix('/') {
            if stack.pop().as_deref() != Some(name.trim()) {
                return Err(UpfError::Parse("unbalanced tags".into()));
            }
            if stack.is_empty() {
                closed_root = true;
            }
            continue;
        }
        if closed_root {
            return Err(UpfError::Parse("multiple roots".into()));
        }
        let empty = body.ends_with('/');
        let body = body.trim_end_matches('/').trim();
        let n = body.find(char::is_whitespace).unwrap_or(body.len());
        let name = &body[..n];
        if name.starts_with('!') || name.is_empty() {
            return Err(UpfError::Parse("unsupported XML declaration".into()));
        }
        if stack.is_empty() && name != "UPF" {
            return Err(UpfError::Parse("expected UPF root".into()));
        }
        let mut attrs = BTreeMap::new();
        let mut a = body[n..].trim();
        while !a.is_empty() {
            let eq = a
                .find('=')
                .ok_or_else(|| UpfError::Parse("attribute without value".into()))?;
            let key = a[..eq].trim();
            a = a[eq + 1..].trim_start();
            let q = a
                .chars()
                .next()
                .filter(|c| *c == '\'' || *c == '"')
                .ok_or_else(|| UpfError::Parse("unquoted attribute".into()))?;
            a = &a[1..];
            let i = a
                .find(q)
                .ok_or_else(|| UpfError::Parse("unclosed attribute".into()))?;
            let value = decode_attribute(&a[..i])?;
            if attrs.insert(key.into(), value).is_some() {
                return Err(invalid(key, "duplicate attribute"));
            }
            a = a[i + 1..].trim_start();
        }
        if !stack.iter().any(|s| s == "PP_INFO") {
            let expected = match name {
                "UPF" => None,
                "PP_R" | "PP_RAB" => Some("PP_MESH"),
                "PP_DIJ" | "PP_BETA.1" | "PP_BETA.2" | "PP_BETA.3" => Some("PP_NONLOCAL"),
                "PP_CHI.1" | "PP_CHI.2" | "PP_CHI.3" => Some("PP_PSWFC"),
                "PP_INFO" | "PP_HEADER" | "PP_MESH" | "PP_LOCAL" | "PP_NONLOCAL" | "PP_RHOATOM"
                | "PP_PSWFC" => Some("UPF"),
                _ => return Err(UpfError::Unsupported(name.into())),
            };
            if stack.last().map(String::as_str) != expected {
                return Err(invalid(name, "unexpected section nesting"));
            }
        }
        // Informational input may contain arbitrary tags; never interpret them as data.
        if !stack.iter().any(|s| s == "PP_INFO") {
            out.push(Tag {
                name: name.into(),
                attrs,
            });
        }
        if !empty {
            stack.push(name.into());
        } else if stack.is_empty() {
            closed_root = true;
        }
    }
    if !stack.is_empty() || !closed_root || !rest.trim().is_empty() {
        return Err(UpfError::Parse(
            "truncated document or trailing text".into(),
        ));
    }
    Ok(out)
}
fn attr<'a>(tag: &'a Tag, key: &str) -> Result<&'a str, UpfError> {
    tag.attrs
        .get(key)
        .map(String::as_str)
        .ok_or_else(|| UpfError::Missing(format!("{}/{key}", tag.name)))
}
fn flag(s: &str) -> Option<bool> {
    match s.trim().to_ascii_lowercase().as_str() {
        "t" | ".t." | "true" | ".true." | "1" => Some(true),
        "f" | ".f." | "false" | ".false." | "0" => Some(false),
        _ => None,
    }
}
pub fn from_file(path: impl AsRef<Path>) -> Result<Pseudopotential, UpfError> {
    parse_upf(&std::fs::read_to_string(path)?)
}
pub fn parse_upf(text: &str) -> Result<Pseudopotential, UpfError> {
    let metadata = tags(text)?;
    let find = |name: &str| -> Result<&Tag, UpfError> {
        let mut matching = metadata.iter().filter(|t| t.name == name);
        let first = matching
            .next()
            .ok_or_else(|| UpfError::Missing(name.into()))?;
        if matching.next().is_some() {
            return Err(invalid(name, "duplicate section"));
        }
        Ok(first)
    };
    let version = attr(find("UPF")?, "version")?;
    if !matches!(version, "2.0" | "2.0.0" | "2.0.1") {
        return Err(UpfError::Version(version.into()));
    }
    let h = find("PP_HEADER")?;
    for (key, feature) in [
        ("is_paw", "PAW"),
        ("is_ultrasoft", "ultrasoft"),
        ("has_so", "spin-orbit"),
        ("core_correction", "core correction"),
        ("is_coulomb", "bare Coulomb"),
        ("has_gipaw", "GIPAW"),
        ("has_wfc", "all-electron wavefunctions"),
    ] {
        match flag(attr(h, key)?) {
            Some(true) => return Err(UpfError::Unsupported(feature.into())),
            Some(false) => (),
            None => return Err(invalid(key, "invalid boolean")),
        }
    }
    if let Some(s) = h.attrs.get("paw_as_gipaw")
        && flag(s) != Some(false)
    {
        return Err(UpfError::Unsupported("paw_as_gipaw".into()));
    }
    if attr(h, "pseudo_type")? != "NC" {
        return Err(UpfError::Unsupported(format!(
            "pseudo_type={}",
            attr(h, "pseudo_type")?
        )));
    }
    if !matches!(
        attr(h, "relativistic")?,
        "scalar" | "nonrelativistic" | "no"
    ) {
        return Err(UpfError::Unsupported(format!(
            "relativistic={}",
            attr(h, "relativistic")?
        )));
    }
    let functional = attr(h, "functional")?.to_string();
    if functional.trim().is_empty() {
        return Err(invalid("functional", "empty label"));
    }
    for t in &metadata {
        if !matches!(
            t.name.as_str(),
            "UPF"
                | "PP_INFO"
                | "PP_HEADER"
                | "PP_MESH"
                | "PP_R"
                | "PP_RAB"
                | "PP_LOCAL"
                | "PP_NONLOCAL"
                | "PP_DIJ"
                | "PP_BETA.1"
                | "PP_BETA.2"
                | "PP_BETA.3"
                | "PP_RHOATOM"
                | "PP_PSWFC"
                | "PP_CHI.1"
                | "PP_CHI.2"
                | "PP_CHI.3"
        ) {
            return Err(UpfError::Unsupported(t.name.clone()));
        }
        if let Some(unit) = t.attrs.get("units") {
            let expected = match t.name.as_str() {
                "PP_LOCAL" | "PP_DIJ" => "Ry",
                "PP_R" | "PP_RAB" => "Bohr",
                n if n.starts_with("PP_BETA.") => "Bohr^-1/2",
                _ => "",
            };
            if unit != expected {
                return Err(UpfError::Unsupported(format!("{} units={unit}", t.name)));
            }
        }
    }
    find("PP_LOCAL")?;
    let np: usize = attr(h, "number_of_proj")?
        .parse()
        .map_err(|_| invalid("number_of_proj", "not an integer"))?;
    if np > 3 {
        return Err(UpfError::Unsupported(
            "more than three projectors (upf reader limit)".into(),
        ));
    }
    for i in 1..=np {
        find(&format!("PP_BETA.{i}"))?;
    }
    // The vendor's serde array decoder separates numbers on spaces only.
    // Replace one ASCII whitespace character with one space: numeric spelling,
    // non-ASCII text and byte offsets are preserved, but line numbers are not.
    let normalized: String = text
        .chars()
        .map(|c| {
            if c.is_ascii_whitespace() || c == '\x0b' {
                ' '
            } else {
                c
            }
        })
        .collect();
    let raw = upf::from_str(&normalized).map_err(|e| match e {
        upf::UpfError::Validation(s) => {
            let field = ["PP_RAB", "PP_LOCAL", "PP_RHOATOM", "PP_R"]
                .into_iter()
                .find(|f| s.contains(f))
                .unwrap_or("UPF");
            invalid(field, &s)
        }
        other => UpfError::Parse(other.to_string()),
    })?;
    for t in &metadata {
        if let Some(size) = t.attrs.get("size") {
            let declared: usize = size
                .parse()
                .map_err(|_| invalid(&t.name, "invalid size attribute"))?;
            let actual = match t.name.as_str() {
                "PP_R" => raw.mesh.r.values.len(),
                "PP_RAB" => raw.mesh.rab.values.len(),
                "PP_LOCAL" => raw.local.values.len(),
                "PP_RHOATOM" => raw.rhoatom.values.len(),
                "PP_DIJ" => raw.nonlocal.dij.values.len(),
                name if name.starts_with("PP_BETA.") => {
                    let index: usize = name[8..]
                        .parse()
                        .map_err(|_| invalid(name, "invalid projector tag"))?;
                    let node = raw
                        .nonlocal
                        .betas
                        .get(index - 1)
                        .ok_or_else(|| invalid(name, "missing projector"))?;
                    match node {
                        PpBetaNode::Beta1(b) | PpBetaNode::Beta2(b) | PpBetaNode::Beta3(b) => {
                            b.values.len()
                        }
                    }
                }
                _ => continue,
            };
            if actual != declared {
                return Err(invalid(&t.name, "size attribute mismatch"));
            }
        }
    }
    let original = raw.clone();
    let header = raw.header;
    let element = header.element.trim().to_string();
    if element.is_empty() || element.len() > 2 || !element.bytes().all(|b| b.is_ascii_alphabetic())
    {
        return Err(invalid("element", "expected one element symbol"));
    }
    if !header.z_valence.is_finite() || header.z_valence <= 0.0 {
        return Err(invalid("z_valence", "must be finite and positive"));
    }
    if raw.mesh.mesh != header.mesh_size {
        return Err(invalid("PP_MESH/mesh", "count mismatch"));
    }
    let mesh = crate::radial::RadialMesh::new(raw.mesh.r.values, raw.mesh.rab.values)
        .map_err(|e| invalid("mesh", &e.to_string()))?;
    if raw.local.values.iter().any(|v| !v.is_finite()) {
        return Err(invalid("PP_LOCAL", "non-finite sample"));
    }
    if raw.rhoatom.values.iter().any(|v| !v.is_finite()) {
        return Err(invalid("PP_RHOATOM", "non-finite sample"));
    }
    if let Some(wfc) = &raw.pswfc {
        use upf::model::PpChiNode;
        for node in &wfc.orbitals {
            let w = match node {
                PpChiNode::Chi1(w) | PpChiNode::Chi2(w) | PpChiNode::Chi3(w) => w,
            };
            if w.values.len() != header.mesh_size
                || w.values.iter().any(|v| !v.is_finite())
                || w.occupation.is_some_and(|x| !x.is_finite() || x < 0.0)
            {
                return Err(invalid("PP_PSWFC", "invalid length, values or occupation"));
            }
        }
    }
    let orbitals = raw.pswfc.as_ref().map_or(0, |w| w.orbitals.len());
    if orbitals != header.number_of_wfc {
        return Err(invalid("number_of_wfc", "count mismatch"));
    }
    let mut projectors = Vec::new();
    for (i, node) in raw.nonlocal.betas.into_iter().enumerate() {
        let (tag_index, b) = match node {
            PpBetaNode::Beta1(b) => (1, b),
            PpBetaNode::Beta2(b) => (2, b),
            PpBetaNode::Beta3(b) => (3, b),
        };
        if tag_index != b.index {
            return Err(invalid("PP_BETA/index", "tag and index disagree"));
        }
        if b.index != i + 1 {
            return Err(invalid(
                "PP_BETA/index",
                "expected consecutive one-based indices",
            ));
        }
        if b.angular_momentum > 3 {
            return Err(UpfError::Unsupported(format!(
                "projector l={}",
                b.angular_momentum
            )));
        }
        if b.angular_momentum > header.l_max {
            return Err(invalid("l_max", "smaller than projector l"));
        }
        if !b.cutoff_radius.is_finite() || b.cutoff_radius < 0.0 {
            return Err(invalid(
                "PP_BETA/cutoff_radius",
                "must be finite and nonnegative",
            ));
        }
        let t = find(&format!("PP_BETA.{}", i + 1))?;
        let cutoff: usize = attr(t, "cutoff_radius_index")?
            .parse()
            .map_err(|_| invalid("cutoff_radius_index", "not an integer"))?;
        if cutoff == 0
            || cutoff > mesh.r().len()
            || (b.values.len() != cutoff && b.values.len() != mesh.r().len())
        {
            return Err(invalid(
                "cutoff_radius_index",
                "inconsistent projector length",
            ));
        }
        if b.values.iter().any(|x| !x.is_finite()) || b.values[cutoff..].iter().any(|x| *x != 0.0) {
            return Err(invalid(
                "PP_BETA",
                "non-finite samples or nonzero samples beyond cutoff",
            ));
        }
        if mesh.r()[0] == 0.0 && b.values[0] != 0.0 {
            return Err(invalid("PP_BETA", "r*beta must vanish at origin"));
        }
        projectors.push(Projector {
            l: b.angular_momentum,
            radial_times_r: b.values[..cutoff].to_vec(),
            cutoff,
        });
    }
    if projectors.len() != np {
        return Err(invalid("number_of_proj", "count mismatch"));
    }
    let dij: Vec<_> = raw
        .nonlocal
        .dij
        .values
        .into_iter()
        .map(|v| Hartree::from(Rydberg(v)).0)
        .collect();
    if dij.len() != np * np || dij.iter().any(|x| !x.is_finite()) {
        return Err(invalid("PP_DIJ", "dimension mismatch or non-finite value"));
    }
    for i in 0..np {
        for j in 0..np {
            if (dij[i * np + j] - dij[j * np + i]).abs() > 1e-12 {
                return Err(UpfError::NonSymmetricD);
            }
            if projectors[i].l != projectors[j].l && dij[i * np + j] != 0.0 {
                return Err(invalid(
                    "PP_DIJ",
                    "coupling between different angular momenta",
                ));
            }
        }
    }
    Ok(Pseudopotential {
        element,
        z_valence: header.z_valence,
        functional,
        mesh,
        local: raw
            .local
            .values
            .into_iter()
            .map(|v| Hartree::from(Rydberg(v)).0)
            .collect(),
        projectors,
        dij,
        raw: original,
    })
}
