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
/// Atomic pseudo-wavefunction, chi(r) = r R(r), on the full radial mesh.
/// Occupation is informational and does not weight an atomic starting basis.
#[derive(Debug, Clone)]
pub struct AtomicWavefunction {
    pub l: usize,
    pub occupation: Option<f64>,
    pub label: Option<String>,
    pub radial_times_r: Vec<f64>,
}
#[derive(Debug)]
pub struct Pseudopotential {
    element: String,
    z_valence: f64,
    functional: String,
    pub(crate) mesh: crate::radial::RadialMesh,
    pub(crate) local: Vec<f64>,
    projectors: Vec<Projector>,
    atomic_wavefunctions: Option<Vec<AtomicWavefunction>>,
    dij: Vec<f64>,
    raw: upf::UpfData,
}
impl Pseudopotential {
    /// Reader compatibility model in UPF units, never a numerical-kernel input.
    /// Original facts come from header_facts. Missing informational attributes
    /// use reader-only neutral defaults; numbered tags above 3 use enum variant
    /// 3 with their original index/order and every sample preserved.
    pub fn raw(&self) -> &upf::UpfData {
        &self.raw
    }
    pub fn atomic_wavefunctions(&self) -> Option<&[AtomicWavefunction]> {
        self.atomic_wavefunctions.as_deref()
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
    /// Immutable spherical model core density, electrons/Bohr^3 (not 4*pi*r^2*rho).
    pub fn core_density(&self) -> Option<&[f64]> {
        self.raw.nlcc.as_ref().map(|n| n.values.as_slice())
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
    start: usize,
    end: usize,
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
    scan_tags(text, false)
}
fn scan_tags(text: &str, inspection: bool) -> Result<Vec<Tag>, UpfError> {
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
        let start_offset = text.len() - rest.len();
        let end_offset = start_offset + end + 1;
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
                "PP_DIJ" => Some("PP_NONLOCAL"),
                n if n.starts_with("PP_BETA.") => Some("PP_NONLOCAL"),
                n if n.starts_with("PP_CHI.") => Some("PP_PSWFC"),
                "PP_INFO" | "PP_HEADER" | "PP_MESH" | "PP_LOCAL" | "PP_NONLOCAL" | "PP_RHOATOM"
                | "PP_PSWFC" | "PP_NLCC" => Some("UPF"),
                _ if inspection => stack.last().map(String::as_str),
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
                start: start_offset,
                end: end_offset,
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
/// Header facts for inspection only, including unsupported files. These facts do
/// not validate numerical arrays and must never authorize an engine calculation.
#[derive(Debug)]
pub struct HeaderFacts {
    pub version: String,
    pub attributes: BTreeMap<String, String>,
    pub projector_counts: BTreeMap<usize, usize>,
    /// File order, one-based index and angular momentum (index, l).
    pub projector_channels: Vec<(usize, usize)>,
    /// Inspection findings only; successful inspection never authorizes SCF.
    pub blockers: Vec<String>,
}
pub fn header_facts(text: &str) -> Result<HeaderFacts, UpfError> {
    let metadata = scan_tags(text, true)?;
    let headers: Vec<_> = metadata.iter().filter(|t| t.name == "PP_HEADER").collect();
    if headers.len() != 1 {
        return Err(invalid("PP_HEADER", "expected exactly one header"));
    }
    let root = metadata
        .first()
        .ok_or_else(|| UpfError::Missing("UPF".into()))?;
    let mut projector_counts = BTreeMap::new();
    let mut projector_channels = Vec::new();
    let mut blockers = Vec::new();
    for tag in metadata.iter().filter(|t| t.name.starts_with("PP_BETA.")) {
        let integer = |key| {
            tag.attrs
                .get(key)
                .and_then(|s| s.trim().parse::<usize>().ok())
        };
        let (Some(index), Some(l)) = (integer("index"), integer("angular_momentum")) else {
            blockers.push(format!(
                "{}.channel: missing or invalid index/angular_momentum; supply valid NC channels",
                tag.name
            ));
            continue;
        };
        if tag.name[8..].parse::<usize>().ok() != Some(index)
            || index != projector_channels.len() + 1
            || l > 3
        {
            blockers.push(format!("{}.channel: malformed or unsupported channel; require consecutive one-based indices and l <= 3", tag.name));
        }
        projector_channels.push((index, l));
        *projector_counts.entry(l).or_insert(0) += 1;
    }
    let a = &headers[0].attrs;
    for (key, capability) in [
        ("is_paw", "PAW augmentation"),
        ("is_ultrasoft", "USPP augmentation"),
        ("has_so", "SOC"),
    ] {
        if a.get(key).and_then(|s| flag(s.trim())) == Some(true) {
            blockers.push(format!("PP_HEADER.{key}: {capability} unsupported; use a supported NC dataset or implement this capability"));
        }
    }
    if a.get("pseudo_type").is_some_and(|s| s.trim() != "NC") {
        blockers
            .push("PP_HEADER.pseudo_type: only NC supported; PAW/USPP require augmentation".into());
    }
    if a.get("relativistic")
        .is_some_and(|s| !matches!(s.trim(), "scalar" | "nonrelativistic" | "no"))
    {
        blockers.push("PP_HEADER.relativistic: full-relativistic/SOC unsupported; use scalar or nonrelativistic NC".into());
    }
    let count = |key| a.get(key).and_then(|s| s.trim().parse::<usize>().ok());
    if count("number_of_proj") != Some(projector_channels.len()) {
        blockers
            .push("PP_HEADER.number_of_proj: channel count mismatch; repair the dataset".into());
    }
    if count("number_of_proj").is_some_and(|n| n > 0)
        && !metadata.iter().any(|t| t.name == "PP_DIJ")
    {
        blockers.push(
            "PP_NONLOCAL.PP_DIJ: missing coupling matrix; supply a complete NC dataset".into(),
        );
    }
    Ok(HeaderFacts {
        version: attr(root, "version")?.into(),
        attributes: headers[0].attrs.clone(),
        projector_counts,
        projector_channels,
        blockers,
    })
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
                | "PP_NLCC"
                | "PP_PSWFC"
                | "PP_CHI.1"
                | "PP_CHI.2"
                | "PP_CHI.3"
        ) && !t.name.starts_with("PP_BETA.")
            && !t.name.starts_with("PP_CHI.")
        {
            return Err(UpfError::Unsupported(t.name.clone()));
        }
        if let Some(unit) = t.attrs.get("units") {
            let expected = match t.name.as_str() {
                "PP_LOCAL" | "PP_DIJ" => "Ry",
                "PP_R" | "PP_RAB" => "Bohr",
                "PP_NLCC" => "Bohr^-3",
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
        .trim()
        .parse()
        .map_err(|_| invalid("number_of_proj", "not an integer"))?;
    if np > metadata.len() {
        return Err(invalid("number_of_proj", "count exceeds document sections"));
    }
    for i in 1..=np {
        find(&format!("PP_BETA.{i}"))?;
    }
    if np > 0 {
        find("PP_DIJ")?;
    }
    // The vendor's serde array decoder separates numbers on spaces only.
    // Replace one ASCII whitespace character with one space: numeric spelling,
    // non-ASCII text and byte offsets are preserved, but line numbers are not.
    let normalized: String = reader_normalization(text, &metadata)?
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
                .trim()
                .parse()
                .map_err(|_| invalid(&t.name, "invalid size attribute"))?;
            let actual = match t.name.as_str() {
                "PP_R" => raw.mesh.r.values.len(),
                "PP_RAB" => raw.mesh.rab.values.len(),
                "PP_LOCAL" => raw.local.values.len(),
                "PP_RHOATOM" => raw.rhoatom.values.len(),
                "PP_NLCC" => raw.nlcc.as_ref().map_or(0, |c| c.values.len()),
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
    let core_enabled = flag(attr(h, "core_correction")?)
        .ok_or_else(|| invalid("core_correction", "invalid boolean"))?;
    if core_enabled != raw.nlcc.is_some() {
        return Err(invalid(
            "PP_NLCC",
            "core_correction and core section disagree",
        ));
    }
    if let Some(core) = &raw.nlcc
        && (core.values.len() != raw.header.mesh_size || core.values.iter().any(|v| !v.is_finite()))
    {
        return Err(invalid(
            "PP_NLCC",
            "mesh size mismatch or non-finite density",
        ));
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
    let mut atomic_wavefunctions = raw.pswfc.as_ref().map(|_| Vec::new());
    if let Some(wfc) = &raw.pswfc {
        use upf::model::PpChiNode;
        for (i, node) in wfc.orbitals.iter().enumerate() {
            let w = match node {
                PpChiNode::Chi1(w) | PpChiNode::Chi2(w) | PpChiNode::Chi3(w) => w,
            };
            let tag = find(&format!("PP_CHI.{}", i + 1))?;
            let l = attr(tag, "l")?
                .trim()
                .parse::<usize>()
                .map_err(|_| invalid("PP_CHI/l", "expected nonnegative integer"))?;
            atomic_wavefunctions
                .as_mut()
                .unwrap()
                .push(AtomicWavefunction {
                    l,
                    occupation: w.occupation,
                    label: tag.attrs.get("label").cloned(),
                    radial_times_r: w.values.clone(),
                });
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
            PpBetaNode::Beta3(b) => (i + 1, b),
        };
        let original_tag = metadata
            .iter()
            .filter(|t| t.name.starts_with("PP_BETA."))
            .nth(i)
            .ok_or_else(|| invalid("PP_BETA/index", "missing original channel"))?;
        if original_tag.name[8..].parse::<usize>().ok() != Some(b.index) || tag_index != b.index {
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
            .trim()
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
        // A finite endpoint sample at r=0 has exactly zero weight in r*u*j_l.
        // Preserve the supplied value (ONCV interpolation can leave a residual);
        // do not introduce a tolerance or alter any nonzero-radius sample.
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
    let matrix_size = np
        .checked_mul(np)
        .ok_or_else(|| invalid("PP_DIJ", "dimension overflow"))?;
    if dij.len() != matrix_size || dij.iter().any(|x| !x.is_finite()) {
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
        atomic_wavefunctions,
        dij,
        raw: original,
    })
}

/// Compatibility normalization for the narrow third-party serde model. This
/// never changes numeric array text, projector indices, D order or occupations.
/// Optional scalar metadata is not used by our radial quadrature: PP_R/PP_RAB
/// are authoritative. Missing cutoff hints remain absent in header_facts.
fn reader_normalization(text: &str, metadata: &[Tag]) -> Result<String, UpfError> {
    fn escaped(s: &str) -> String {
        s.replace('&', "&amp;")
            .replace('"', "&quot;")
            .replace('<', "&lt;")
            .replace('>', "&gt;")
    }
    let mut normalized = text.to_string();
    let mesh_size = metadata
        .iter()
        .find(|t| t.name == "PP_HEADER")
        .and_then(|t| t.attrs.get("mesh_size"))
        .ok_or_else(|| UpfError::Missing("mesh_size".into()))?;
    for tag in metadata.iter().rev() {
        let mut attrs = tag.attrs.clone();
        match tag.name.as_str() {
            "PP_HEADER" => {
                attrs.entry("wfc_cutoff".into()).or_insert("0".into());
                attrs.entry("l_max_rho".into()).or_insert("0".into());
            }
            "PP_MESH" => {
                attrs.entry("mesh".into()).or_insert(mesh_size.clone());
                // Unused informational parameters, not inferred log-grid data.
                for key in ["dx", "xmin", "rmax", "zmesh"] {
                    attrs.entry(key.into()).or_insert("0".into());
                }
            }
            _ => (),
        }
        let mut opening = format!("<{}", tag.name);
        for (key, value) in attrs {
            // UPF numeric attributes commonly have fixed-width leading spaces.
            let value = if matches!(
                key.as_str(),
                "generated"
                    | "author"
                    | "date"
                    | "comment"
                    | "functional"
                    | "element"
                    | "relativistic"
                    | "label"
            ) {
                value.as_str()
            } else {
                value.trim()
            };
            opening.push_str(&format!(" {key}=\"{}\"", escaped(value)));
        }
        opening.push_str(if text[tag.start..tag.end].trim_end().ends_with("/>") {
            "/>"
        } else {
            ">"
        });
        normalized.replace_range(tag.start..tag.end, &opening);
    }
    // The original balanced tree was already checked. Renaming only XML names
    // maps extra numbered nodes to a repeatable enum variant, preserving index.
    // Validate consecutive original suffixes first (including pseudo orbitals).
    for prefix in ["PP_BETA.", "PP_CHI."] {
        for (i, tag) in metadata
            .iter()
            .filter(|t| t.name.starts_with(prefix))
            .enumerate()
        {
            let index = tag.name[prefix.len()..]
                .parse::<usize>()
                .map_err(|_| invalid(prefix, "invalid numbered tag"))?;
            if tag
                .attrs
                .get("index")
                .is_some_and(|v| v.trim().parse::<usize>().ok() != Some(index))
            {
                return Err(invalid("PP_BETA/index", "tag and index disagree"));
            }
            if index != i + 1 {
                return Err(invalid(
                    prefix,
                    "expected consecutive one-based numbered tags",
                ));
            }
            if index > 3 {
                for start in ["<", "</"] {
                    let from = format!("{start}{}", tag.name);
                    let to = format!("{start}{prefix}3");
                    let mut offset = 0;
                    while let Some(pos) = normalized[offset..].find(&from) {
                        let pos = offset + pos;
                        let end = pos + from.len();
                        if normalized[end..]
                            .starts_with(|c: char| c.is_ascii_whitespace() || c == '>' || c == '/')
                        {
                            normalized.replace_range(pos..end, &to);
                            offset = pos + to.len();
                        } else {
                            offset = end;
                        }
                    }
                }
            }
        }
    }
    Ok(normalized)
}
