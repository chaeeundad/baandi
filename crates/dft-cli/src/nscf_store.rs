//! R1 `nscf_bands_v1`: one BAANDICP container per completed (channel, k)
//! unit and a manifest written last (`docs/checkpoint-format.md`, section
//! "NSCF band units"). The engine owns the numerics and re-certifies every
//! resumed unit with a fresh H application; this adapter owns bytes only.
//!
//! Unit file `nscf-c{channel}-k{k:06}.cp`: the common envelope (magic, header
//! length, payload length, JSON header, payload, SHA-256). Header fields are
//! listed in `UnitHeader`; f64 values are IEEE bits. Payload: one section per
//! rank (header `sections`, byte offset/count, array index = rank): u64 LE
//! owned-row count R, R u64 LE global plane-wave indices, then R x bands
//! complex coefficients column-major as (re, im) f64 LE.
use crate::checkpoint_io::{Error, MAX_FILE, MAX_HEADER, envelope, file_bytes, publish_bytes};
use dft_engine::{
    model::ScfError,
    refscf::nscf::{BandStore, Bands},
};
use dft_mpi::Collective;
use dft_numerics::linalg::ComplexMatrix;
use dft_physics::data::Complex64;
use serde::{Deserialize, Serialize};
use sha2::{Digest, Sha256};
use std::{
    collections::BTreeMap,
    path::{Path, PathBuf},
};

pub const KIND: &str = "nscf_bands_v1";
pub const MANIFEST: &str = "nscf-manifest.json";

#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct Section {
    offset: usize,
    count: usize,
}
#[derive(Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
struct UnitHeader {
    schema_version: u32,
    kind: String,
    /// Run identity; every key must match on resume.
    fingerprint: BTreeMap<String, String>,
    channel: usize,
    k_index: usize,
    bands: usize,
    basis_size: usize,
    eigenvalues_ha_bits: Vec<u64>,
    residuals_ha_bits: Vec<u64>,
    sections: Vec<Section>,
}

fn invalid(s: impl Into<String>) -> ScfError {
    ScfError::Numerical(format!("nscf checkpoint: {}", s.into()))
}

/// Every rank learns whether any rank failed; the first local message wins.
fn agreed<T, K: Collective>(local: Result<T, ScfError>, comm: &K) -> Result<T, ScfError> {
    let mut flag = [0.];
    comm.tree_sum_f64(&[f64::from(local.is_err())], &mut flag)
        .map_err(|e| invalid(e.to_string()))?;
    match local {
        Err(e) => Err(e),
        Ok(_) if flag[0] != 0. => Err(invalid("another rank failed")),
        ok => ok,
    }
}

pub struct NscfStore {
    dir: PathBuf,
    fingerprint: BTreeMap<String, String>,
    /// Units loaded from existing files, and units written in this run.
    pub resumed: Vec<(usize, usize)>,
    pub written: Vec<(usize, usize)>,
}

impl NscfStore {
    /// `dir` is created when missing. `fingerprint` is the run identity
    /// (input/UPF/density hashes, solver settings, ranks, backend, ...).
    pub fn new(dir: &Path, fingerprint: BTreeMap<String, String>) -> Result<Self, Error> {
        std::fs::create_dir_all(dir)?;
        Ok(Self {
            dir: dir.to_owned(),
            fingerprint,
            resumed: Vec::new(),
            written: Vec::new(),
        })
    }
    pub fn unit_name(channel: usize, k: usize) -> String {
        format!("nscf-c{channel}-k{k:06}.cp")
    }

    /// Encode one unit from every rank's packed section.
    fn encode(
        &self,
        channel: usize,
        k: usize,
        basis_size: usize,
        bands: &Bands,
        sections: &[Vec<u8>],
    ) -> Result<Vec<u8>, Error> {
        let mut payload = Vec::new();
        let mut refs = Vec::new();
        for s in sections {
            refs.push(Section {
                offset: payload.len(),
                count: s.len(),
            });
            payload.extend_from_slice(s);
        }
        let header = serde_json::to_vec(&UnitHeader {
            schema_version: 1,
            kind: KIND.into(),
            fingerprint: self.fingerprint.clone(),
            channel,
            k_index: k,
            bands: bands.eigenvalues.len(),
            basis_size,
            eigenvalues_ha_bits: bands.eigenvalues.iter().map(|v| v.to_bits()).collect(),
            residuals_ha_bits: bands.residuals.iter().map(|v| v.to_bits()).collect(),
            sections: refs,
        })?;
        if header.len() > MAX_HEADER || header.len() + payload.len() + 56 > MAX_FILE {
            return Err("nscf checkpoint unit exceeds the container quota".into());
        }
        let mut out = crate::checkpoint_io::MAGIC.to_vec();
        out.extend_from_slice(&(header.len() as u64).to_le_bytes());
        out.extend_from_slice(&(payload.len() as u64).to_le_bytes());
        out.extend(header);
        out.extend(payload);
        let digest = Sha256::digest(&out);
        out.extend_from_slice(&digest);
        Ok(out)
    }

    #[allow(clippy::too_many_arguments)]
    fn decode_local(
        &self,
        bytes: &[u8],
        rank: usize,
        ranks: usize,
        channel: usize,
        k: usize,
        rows: &[usize],
        bands: usize,
    ) -> Result<(Vec<f64>, ComplexMatrix), String> {
        let (header, payload) = envelope(bytes).map_err(|e| e.to_string())?;
        let h: UnitHeader = serde_json::from_slice(header).map_err(|e| e.to_string())?;
        if h.schema_version != 1 || h.kind != KIND {
            return Err("unsupported unit kind/version".into());
        }
        if let Some(key) = self
            .fingerprint
            .keys()
            .chain(h.fingerprint.keys())
            .find(|key| self.fingerprint.get(*key) != h.fingerprint.get(*key))
        {
            return Err(format!("fingerprint mismatch on '{key}'"));
        }
        if h.channel != channel || h.k_index != k || h.bands != bands {
            return Err("channel/k/band count mismatch".into());
        }
        if h.sections.len() != ranks
            || h.eigenvalues_ha_bits.len() != bands
            || h.residuals_ha_bits.len() != bands
        {
            return Err("rank count or array length mismatch".into());
        }
        let mut cursor = 0;
        for s in &h.sections {
            if s.offset != cursor || s.count < 8 {
                return Err("section layout".into());
            }
            cursor = s
                .offset
                .checked_add(s.count)
                .filter(|&e| e <= payload.len())
                .ok_or("section range")?;
        }
        if cursor != payload.len() {
            return Err("unreferenced payload bytes".into());
        }
        let s = &h.sections[rank];
        let data = &payload[s.offset..s.offset + s.count];
        let word = |i: usize| u64::from_le_bytes(data[8 * i..8 * i + 8].try_into().unwrap());
        let r = word(0) as usize;
        if r != rows.len() || data.len() != 8 * (1 + r) + 16 * r * bands {
            return Err("owned rows differ from this run's partition".into());
        }
        if (0..r).any(|i| word(1 + i) as usize != rows[i]) {
            return Err("owned row indices differ from this run's partition".into());
        }
        let mut vectors = ComplexMatrix::zeros(r, bands).map_err(|e| e.to_string())?;
        for (n, z) in vectors.as_mut_slice().iter_mut().enumerate() {
            let at = 1 + r + 2 * n;
            *z = Complex64::new(f64::from_bits(word(at)), f64::from_bits(word(at + 1)));
        }
        vectors.validate_finite().map_err(|e| e.to_string())?;
        let values: Vec<f64> = h
            .eigenvalues_ha_bits
            .iter()
            .map(|&b| f64::from_bits(b))
            .collect();
        if values.iter().any(|v| !v.is_finite()) {
            return Err("nonfinite eigenvalue".into());
        }
        Ok((values, vectors))
    }

    /// Rank 0 writes the manifest after every unit (no-clobber; an existing
    /// identical manifest is accepted). Collective.
    pub fn write_manifest<K: Collective>(
        &self,
        comm: &K,
        units: &[(usize, usize)],
    ) -> Result<(), Error> {
        let local = if comm.rank() == 0 {
            (|| -> Result<(), Error> {
                let mut list = Vec::new();
                for &(channel, k) in units {
                    let name = Self::unit_name(channel, k);
                    let bytes = file_bytes(&self.dir.join(&name))?;
                    list.push(
                        serde_json::json!({"file": name, "channel": channel, "k_index": k,
                        "sha256": format!("{:x}", Sha256::digest(&bytes))}),
                    );
                }
                let text = serde_json::to_vec_pretty(&serde_json::json!({
                    "schema_version": 1, "kind": "nscf_bands_v1_manifest",
                    "fingerprint": self.fingerprint, "units": list,
                }))?;
                let path = self.dir.join(MANIFEST);
                if path.exists() {
                    if std::fs::read(&path)? != text {
                        return Err("existing nscf manifest differs from this run".into());
                    }
                    return Ok(());
                }
                publish_bytes(&path, &text)
            })()
            .map_err(|e| invalid(e.to_string()))
        } else {
            Ok(())
        };
        agreed(local, comm).map_err(|e| e.to_string().into())
    }
}

impl<K: Collective> BandStore<K> for NscfStore {
    fn load(
        &mut self,
        comm: &K,
        channel: usize,
        k: usize,
        rows: &[usize],
        bands: usize,
    ) -> Result<Option<(Vec<f64>, ComplexMatrix)>, ScfError> {
        let path = self.dir.join(Self::unit_name(channel, k));
        let local = if path.exists() {
            file_bytes(&path)
                .map_err(|e| e.to_string())
                .and_then(|bytes| {
                    self.decode_local(&bytes, comm.rank(), comm.size(), channel, k, rows, bands)
                })
                .map(Some)
                .map_err(|e| invalid(format!("{}: {e}", path.display())))
        } else {
            Ok(None)
        };
        let loaded = agreed(local, comm)?;
        if loaded.is_some() {
            self.resumed.push((channel, k));
        }
        Ok(loaded)
    }

    fn save(
        &mut self,
        comm: &K,
        channel: usize,
        k: usize,
        global_len: usize,
        bands: &Bands,
    ) -> Result<(), ScfError> {
        let n = bands.eigenvalues.len();
        let r = bands.rows.len();
        // Exact transport: row indices are integers below 2^53 and the
        // coefficients are copied, never combined arithmetically.
        let mut local = Vec::with_capacity(1 + r + 2 * r * n);
        local.push(r as f64);
        local.extend(bands.rows.iter().map(|&g| g as f64));
        local.extend(
            bands
                .eigenvectors
                .as_slice()
                .iter()
                .flat_map(|z| [z.re, z.im]),
        );
        let mut lengths = vec![0.; comm.size()];
        comm.allgatherv_f64(&[local.len() as f64], &vec![1; comm.size()], &mut lengths)
            .map_err(|e| invalid(e.to_string()))?;
        let counts: Vec<usize> = lengths.iter().map(|&l| l as usize).collect();
        let mut all = vec![0.; counts.iter().sum()];
        comm.allgatherv_f64(&local, &counts, &mut all)
            .map_err(|e| invalid(e.to_string()))?;
        drop(local);
        let result = if comm.rank() == 0 {
            let mut sections = Vec::with_capacity(comm.size());
            let mut at = 0;
            for &count in &counts {
                let values = &all[at..at + count];
                let rows = values[0] as usize;
                let mut bytes = Vec::with_capacity(8 * count);
                for (i, v) in values.iter().enumerate() {
                    let word = if i <= rows { *v as u64 } else { v.to_bits() };
                    bytes.extend_from_slice(&word.to_le_bytes());
                }
                sections.push(bytes);
                at += count;
            }
            self.encode(channel, k, global_len, bands, &sections)
                .and_then(|bytes| {
                    publish_bytes(&self.dir.join(Self::unit_name(channel, k)), &bytes)
                })
                .map_err(|e| invalid(e.to_string()))
        } else {
            Ok(())
        };
        agreed(result, comm)?;
        self.written.push((channel, k));
        Ok(())
    }
}
