//! Read-only scalar projection of the shared checkpoint API. Never restores SCF.
use crate::checkpoint_io;
use dft_engine::checkpoint::CheckpointKind;
use serde_json::{Value, json};
use std::{fs::File, io::Read, path::Path};

pub(crate) fn recognizes(path: &Path) -> Result<bool, checkpoint_io::Error> {
    if path.is_dir() {
        return Ok(false);
    }
    // Extensions route damaged/truncated containers through the checkpoint reader;
    // magic also recognizes valid checkpoints with arbitrary names.
    if path.extension().is_some_and(|e| e == "cp") {
        return Ok(true);
    }
    let mut prefix = Vec::with_capacity(8);
    File::open(path)?.take(8).read_to_end(&mut prefix)?;
    Ok(prefix == b"BAANDICP")
}

pub(crate) fn evidence(path: &Path) -> Result<Value, checkpoint_io::Error> {
    let cp = checkpoint_io::read(path)?;
    let (state, capability, history) =
        match cp.kind {
            CheckpointKind::ExactContinuation if cp.result.is_none() => {
                let c = cp
                    .continuation
                    .as_ref()
                    .ok_or("checkpoint.continuation: required for ExactContinuation")?;
                (
                    "nonterminal: process liveness unknown",
                    "exact_continuation",
                    Some(&c.history),
                )
            }
            CheckpointKind::WarmStart
                if cp.continuation.is_none() && cp.runtime.is_none() && cp.result.is_none() =>
            {
                (
                    "warm_start: convergence not recorded",
                    "initial_guess_only",
                    None,
                )
            }
            CheckpointKind::ResultExport if cp.continuation.is_none() && cp.runtime.is_none() => {
                let r = cp
                    .result
                    .as_ref()
                    .ok_or("checkpoint.result: required for ResultExport")?;
                (
                    if r.converged {
                        "converged"
                    } else {
                        "not_converged"
                    },
                    "inspection_only",
                    Some(&r.history),
                )
            }
            _ => return Err(
                "checkpoint.kind: payload must match ExactContinuation, WarmStart, or ResultExport"
                    .into(),
            ),
        };
    let last = history.and_then(|h| h.last());
    let last_iteration = last.map(|h| {
        json!({
            "iteration":h.iteration, "energy_hartree":h.energy.total_hartree(),
            "energy_change_hartree":h.energy_change,
            "density_residual_bohr_minus_3_over_2":h.density_residual,
            "max_eigen_residual_hartree":h.max_eigen_residual,
            "requested_eigen_tolerance_hartree":h.requested_eigen_tolerance,
            "charge_error_relative":h.charge_error,
            "certificate":h.fixed_point_certificate.as_ref().map(|c| json!({
                "density_residual_bohr_minus_3_over_2":c.density_residual,
                "occupied_eigenvalue_change_hartree":c.occupied_eigenvalue_change
            }))
        })
    });
    Ok(json!({
        "source":path, "format":"BAANDICP", "schema_version":cp.schema_version,
        "kind":format!("{:?}",cp.kind), "capability":capability,
        "state":state, "calc_id":cp.calc_id, "sequence":cp.sequence,
        "iteration":cp.continuation.as_ref().map(|c|c.iteration).or_else(||last.map(|h|h.iteration)),
        "state_version":cp.continuation.as_ref().map(|c|c.state_version),
        "execution_mode":cp.continuation.as_ref().map(|c|format!("{:?}",c.config.mode)),
        "saved_runtime_state":cp.runtime.as_ref().map(|r|format!("{:?}",r.state)),
        "policy_revision":cp.runtime.as_ref().map(|r|r.revision),
        "audit_records":cp.runtime.as_ref().map(|r|r.audit.len()),
        "converged":cp.result.as_ref().map(|r|r.converged),
        "result":cp.result.as_ref().map(|r|json!({"energy_hartree":r.energy.total_hartree(),"charge_electrons":r.charge,"gap_hartree":r.gap_hartree})),
        "last_iteration":last_iteration,
        "model":{"atom_count":cp.model.species.len(),"species":cp.model.species,
            "electrons":cp.model.electrons,"cutoff_hartree":cp.model.cutoff_hartree,
            "kpoint_count":cp.model.kpoints.len(),"xc":cp.model.xc,"spin":cp.model.spin},
        "electronic_state":{"grid":cp.electronic_state.grid,"density_values":cp.electronic_state.density.len(),
            "orbital_blocks":cp.electronic_state.orbitals.len(),"bands":cp.electronic_state.occupations.len()},
        "execution":{"build_sha256":cp.execution.build_sha256,"platform":cp.execution.platform,"threads":cp.execution.threads},
        "validation":"Container checksum, schema and decoding checked; restart compatibility and physical accuracy are not established. No SCF executed.",
        "preflight":"Not rerun: original configuration and pseudopotential files are not required for inspection."
    }))
}
