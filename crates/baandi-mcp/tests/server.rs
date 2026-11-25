use baandi_mcp::{Quotas, ReferenceInput, Server, schema};
use dft_engine::runtime::*;
use dft_engine::{
    input::{Atom, Cell, parse_upf},
    *,
};
use serde_json::{Value, json};
use std::time::{Duration, SystemTime};
fn server() -> Server {
    Server::new(env!("CARGO_MANIFEST_DIR"), Quotas::default()).unwrap()
}
fn problem() -> ReferenceInput {
    // Synthetic local potential, no real UPF file and no reference-energy claim.
    let r: Vec<_> = (0..201).map(|i| i as f64 * 0.06).collect();
    let arr = |v: Vec<f64>| {
        v.iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(" ")
    };
    let xml = format!(
        r#"<UPF version="2.0.1"><PP_INFO>Synthetic only</PP_INFO><PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max_rho="0" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_gipaw="F" has_wfc="F" core_correction="F" functional="SLA PZ" z_valence="2" l_max="0" l_local="-1" mesh_size="201" number_of_wfc="1" number_of_proj="0"/><PP_MESH dx="0.06" mesh="201" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH><PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL/><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,
        arr(r.clone()),
        arr(vec![0.06; 201]),
        arr(r.iter().map(|r| -4.0 / (r * r + 1.0).sqrt()).collect()),
        arr(r.iter().map(|r| r * r * (-r * r).exp()).collect())
    );
    let xml = xml.replace("</UPF>", &format!("<PP_PSWFC><PP_CHI.1 l=\"0\" label=\"s\" occupation=\"2\">{}</PP_CHI.1></PP_PSWFC></UPF>", arr(r.iter().map(|r| r * (-0.6 * r * r).exp()).collect())));
    let p = parse_upf(&xml).unwrap();
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    ReferenceInput {
        cell,
        atoms: vec![Atom::new("Xx", [0.; 3])],
        potentials: vec![("Xx".into(), p)],
        options: ScfOptions {
            cutoff_hartree: 2.5,
            bands: Some(2),
            max_iterations: 2,
            ..ScfOptions::default()
        },
        request: dft_engine::refscf::api::Request::new(
            Default::default(),
            vec![dft_engine::input::kpoints::KPoint {
                fractional: [0.; 3],
                weight: 1.,
            }],
        ),
    }
}

fn handshake(s: &mut Server) {
    let r=s.handle(&json!({"jsonrpc":"2.0","id":1,"method":"initialize","params":{"protocolVersion":"2025-06-18","capabilities":{},"clientInfo":{"name":"test","version":"1"}}}).to_string()).unwrap();
    assert_eq!(r["result"]["protocolVersion"], "2025-06-18");
    assert!(
        s.handle(r#"{"jsonrpc":"2.0","method":"notifications/initialized"}"#)
            .is_none()
    );
}
fn call(s: &mut Server, name: &str, args: Value) -> Value {
    s.handle(&json!({"jsonrpc":"2.0","id":2,"method":"tools/call","params":{"name":name,"arguments":args}}).to_string()).unwrap()
}
fn data(r: Value) -> Value {
    assert_eq!(r["result"]["isError"], false, "{r}");
    serde_json::from_str(r["result"]["content"][0]["text"].as_str().unwrap()).unwrap()
}
fn prepared(mode: &str) -> Server {
    let mut s = server();
    s.register_input("validated", problem()).unwrap();
    handshake(&mut s);
    data(call(
        &mut s,
        "create_calculation",
        json!({"input_id":"validated","calc_id":"c","mode":mode}),
    ));
    s
}
fn control(s: &Server, id: &str) -> Value {
    let r = s.runtime("c").unwrap();
    let snap = r.snapshot();
    json!({"calc_id":"c","request_id":id,"observed_state_version":snap.state_version,"observed_policy_revision":snap.policy_revision,"max_iteration":20,"reason":"test"})
}
fn approval() -> Approval {
    Approval {
        id: "host-approval".into(),
        approver: "authenticated-operator".into(),
        expires_at: SystemTime::now() + Duration::from_secs(60),
        allow_revalidation: false,
    }
}
#[test]
fn schema_documents_match_discovery_and_types() {
    assert_eq!(
        schema::tools(),
        serde_json::from_str::<Value>(include_str!("../schemas/tools.json")).unwrap()
    );
    let tools = schema::tools();
    assert_eq!(tools.as_array().unwrap().len(), 11);
    for tool in tools.as_array().unwrap() {
        assert!(!schema::validate(&tool["inputSchema"], &json!({})));
        let mut args = json!({"calc_id":"c"});
        match tool["name"].as_str().unwrap() {
            "create_calculation" => args = json!({"input_id":"v","calc_id":"c","mode":"observe"}),
            "run" => args["steps"] = json!(1),
            "get_events" => args["limit"] = json!(1),
            "pause" | "resume" | "cancel" | "propose_policy_change" => {
                args = json!({"calc_id":"c","request_id":"r","observed_state_version":0,"observed_policy_revision":0,"max_iteration":10,"reason":"fixture"});
                if tool["name"] == "propose_policy_change" {
                    args["mixing_beta"] = json!(100.);
                }
                if tool["name"] == "propose_policy_change" {
                    serde_json::from_value::<schema::Policy>(args.clone()).unwrap();
                } else {
                    serde_json::from_value::<schema::Control>(args.clone()).unwrap();
                }
            }
            _ => {
                serde_json::from_value::<schema::Calculation>(args.clone()).unwrap();
            }
        }
        assert!(schema::validate(&tool["inputSchema"], &args), "{tool}");
        args["approval"] = json!(true);
        assert!(!schema::validate(&tool["inputSchema"], &args));
    }
}
#[test]
fn unsafe_policy_is_rejected_by_reference_runtime_without_mutation() {
    let mut s = prepared("supervised");
    let before = s.runtime("c").unwrap().snapshot();
    let mut args = control(&s, "unsafe");
    args["mixing_beta"] = json!(100.);
    assert_eq!(
        data(call(&mut s, "propose_policy_change", args))["disposition"],
        "pending_host_approval"
    );
    assert!(
        s.approve("c", "unsafe", approval())
            .unwrap_err()
            .contains("beta")
    );
    assert_eq!(before, s.runtime("c").unwrap().snapshot());
    assert!(
        s.runtime("c").unwrap().audit().last().unwrap()["outcome"]
            .as_str()
            .unwrap()
            .contains("beta")
    );
}
#[test]
fn supervised_lifecycle_and_export() {
    let mut s = prepared("supervised");
    for (name, state) in [
        ("pause", RuntimeState::Paused),
        ("resume", RuntimeState::Initialized),
        ("cancel", RuntimeState::Cancelled),
    ] {
        let args = control(&s, name);
        data(call(&mut s, name, args));
        assert_ne!(s.runtime("c").unwrap().state(), state);
        s.approve("c", name, approval()).unwrap();
        assert_eq!(s.runtime("c").unwrap().state(), state);
    }
    let export = data(call(&mut s, "export_result", json!({"calc_id":"c"})));
    assert!(
        export["csv"]
            .as_str()
            .unwrap()
            .starts_with("iteration,energy_hartree\n")
    );
    assert_eq!(export["snapshot"]["energy_hartree"], Value::Null);
}
#[test]
fn observed_denies_control_and_stale_supervised_request_rejects() {
    let mut observed = prepared("observe");
    let args = control(&observed, "pause");
    let r = data(call(&mut observed, "pause", args));
    assert!(
        r["disposition"]
            .as_str()
            .unwrap()
            .contains("ModeDisallowsControl")
    );
    let mut s = prepared("supervised");
    let args = control(&s, "stale");
    data(call(&mut s, "pause", args));
    data(call(&mut s, "step", json!({"calc_id":"c"})));
    s.approve("c", "stale", approval()).unwrap_err();
    assert!(
        s.runtime("c")
            .unwrap()
            .audit()
            .iter()
            .any(|a| a["outcome"] == "StaleState")
    );
}
#[test]
fn bounded_run_events_and_scalar_result() {
    let mut s = prepared("observe");
    let snap = data(call(&mut s, "run", json!({"calc_id":"c","steps":2})));
    assert_eq!(snap["iteration"], 2);
    let events = data(call(&mut s, "get_events", json!({"calc_id":"c","limit":1})));
    assert_eq!(events["events"].as_array().unwrap().len(), 1);
    let export = data(call(&mut s, "export_result", json!({"calc_id":"c"})));
    assert!(
        export["snapshot"]["energy_hartree"]
            .as_f64()
            .unwrap()
            .is_finite()
    );
}
#[test]
fn protocol_and_authority_fail_closed() {
    let mut s = server();
    assert_eq!(s.handle("{").unwrap()["error"]["code"], -32700);
    assert_eq!(
        call(&mut s, "status", json!({"calc_id":"c"}))["error"]["code"],
        -32002
    );
    handshake(&mut s);
    for (tool, args) in [
        ("shell", json!({"command":"echo forbidden"})),
        ("approve", json!({})),
        (
            "create_calculation",
            json!({"input_id":"x","calc_id":"c","mode":"autonomous"}),
        ),
        ("status", json!({"calc_id":"c","path":"/etc/passwd"})),
        ("run", json!({"calc_id":"c","steps":101})),
    ] {
        assert_eq!(call(&mut s, tool, args)["error"]["code"], -32602);
    }
    assert!(s.handle(r#"{"jsonrpc":"2.0","method":"tools/call","params":{"name":"create_calculation","arguments":{}}}"#).is_none());
    assert_eq!(
        s.handle(r#"{"jsonrpc":"2.0","id":{},"method":"ping"}"#)
            .unwrap()["error"]["code"],
        -32600
    );
}
#[test]
fn stdio_frames_notifications_and_quotas() {
    let mut s = server();
    let mut output = vec![];
    s.serve(&b"{\"jsonrpc\":\"2.0\",\"method\":\"ignored\"}\n{\"jsonrpc\":\"2.0\",\"id\":1,\"method\":\"ping\"}\n"[..],&mut output).unwrap();
    assert_eq!(String::from_utf8(output).unwrap().lines().count(), 1);
    let mut s = Server::new(
        env!("CARGO_MANIFEST_DIR"),
        Quotas {
            frame_bytes: 1024,
            requests: 2,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(s.serve(&vec![b'x'; 1025][..], Vec::new()).is_err());
    let mut s = prepared("observe");
    s.register_input("other", problem()).unwrap();
    assert_eq!(
        call(
            &mut s,
            "create_calculation",
            json!({"input_id":"other","calc_id":"c","mode":"observe"})
        )["result"]["isError"],
        true
    );
}

#[test]
fn calculation_and_iteration_quotas_prevent_work() {
    let mut s = Server::new(
        env!("CARGO_MANIFEST_DIR"),
        Quotas {
            calculations: 1,
            total_steps: 0,
            ..Default::default()
        },
    )
    .unwrap();
    s.register_input("only", problem()).unwrap();
    assert!(s.register_input("excess", problem()).is_err());
    handshake(&mut s);
    data(call(
        &mut s,
        "create_calculation",
        json!({"input_id":"only","calc_id":"c","mode":"observe"}),
    ));
    let version = s.runtime("c").unwrap().snapshot().state_version;
    assert_eq!(
        call(&mut s, "step", json!({"calc_id":"c"}))["result"]["isError"],
        true
    );
    assert_eq!(s.runtime("c").unwrap().snapshot().state_version, version);
    let mut s = Server::new(
        env!("CARGO_MANIFEST_DIR"),
        Quotas {
            requests: 1,
            ..Default::default()
        },
    )
    .unwrap();
    assert!(
        s.handle(r#"{"jsonrpc":"2.0","id":1,"method":"ping"}"#)
            .unwrap()
            .get("result")
            .is_some()
    );
    assert_eq!(
        s.handle(r#"{"jsonrpc":"2.0","id":2,"method":"ping"}"#)
            .unwrap()["error"]["code"],
        -32000
    );
}

#[test]
fn valid_policy_and_denial_are_host_bound() {
    let mut s = prepared("supervised");
    let mut args = control(&s, "beta");
    args["mixing_beta"] = json!(0.2);
    data(call(&mut s, "propose_policy_change", args));
    assert_eq!(s.runtime("c").unwrap().policy_revision(), 0);
    s.approve("c", "beta", approval()).unwrap();
    assert_eq!(s.runtime("c").unwrap().policy_revision(), 1);
    assert_eq!(
        s.runtime("c").unwrap().audit().last().unwrap()["outcome"],
        "Applied"
    );
    let args = control(&s, "denied");
    data(call(&mut s, "cancel", args));
    s.deny("c", "denied").unwrap();
    assert_eq!(s.runtime("c").unwrap().state(), RuntimeState::Initialized);
}

#[test]
fn reference_policy_allow_list_and_atomic_safe_point_transaction() {
    let mut s = prepared("supervised");
    for key in [
        "eigensolver",
        "max_iterations",
        "density_tolerance",
        "reset_history",
    ] {
        let mut args = control(&s, key);
        args[key] = json!(1);
        let rejected = call(&mut s, "propose_policy_change", args);
        assert_eq!(rejected["error"]["code"], -32602);
        assert!(
            rejected["error"]["message"]
                .as_str()
                .unwrap()
                .contains("supports only")
        );
    }
    let empty = control(&s, "empty");
    assert_eq!(
        call(&mut s, "propose_policy_change", empty)["error"]["code"],
        -32602
    );
    let mut invalid = control(&s, "invalid_history");
    invalid["mixing_beta"] = json!(0.4);
    invalid["mixing_history"] = json!(0);
    data(call(&mut s, "propose_policy_change", invalid));
    let before = s.runtime("c").unwrap().snapshot();
    assert!(s.approve("c", "invalid_history", approval()).is_err());
    assert_eq!(before, s.runtime("c").unwrap().snapshot());

    let mut args = control(&s, "reference_patch");
    args["mixing_beta"] = json!(0.4);
    args["mixing_history"] = json!(4);
    args["threshold_ry"] = json!(1e-9);
    args["verification"] = json!(false);
    data(call(&mut s, "propose_policy_change", args));
    s.approve("c", "reference_patch", approval()).unwrap();
    assert_eq!(s.runtime("c").unwrap().policy_revision(), 1);
    let mut duplicate_boundary = control(&s, "second_patch");
    duplicate_boundary["mixing_beta"] = json!(0.5);
    data(call(&mut s, "propose_policy_change", duplicate_boundary));
    assert!(
        s.approve("c", "second_patch", approval())
            .unwrap_err()
            .contains("one transaction")
    );
    assert_eq!(s.runtime("c").unwrap().policy_revision(), 1);
    let snap = data(call(&mut s, "step", json!({"calc_id":"c"})));
    assert_eq!(snap["driver"], "reference");
    assert_eq!(snap["completed_maps"], 1);
    assert!(snap["last_step"]["dr2_ry"].as_f64().unwrap().is_finite());
}

#[test]
fn expired_approval_and_duplicate_requests_cannot_apply() {
    let mut s = prepared("supervised");
    let args = control(&s, "pause");
    data(call(&mut s, "pause", args.clone()));
    assert_eq!(call(&mut s, "pause", args)["result"]["isError"], true);
    let mut expired = approval();
    expired.expires_at = SystemTime::UNIX_EPOCH;
    assert!(s.approve("c", "pause", expired).is_err());
    assert_eq!(s.runtime("c").unwrap().state(), RuntimeState::Initialized);
    assert!(s.approve("c", "pause", approval()).is_err());
}

#[test]
fn interleaved_calculations_and_bounded_run_have_identical_reference_maps() {
    let mut s = server();
    for id in ["a", "b"] {
        s.register_input(id, problem()).unwrap();
    }
    handshake(&mut s);
    for id in ["a", "b"] {
        data(call(
            &mut s,
            "create_calculation",
            json!({"input_id":id,"calc_id":id,"mode":"observe"}),
        ));
    }
    let first_a = data(call(&mut s, "step", json!({"calc_id":"a"})));
    let first_b = data(call(&mut s, "step", json!({"calc_id":"b"})));
    assert_eq!(first_a["last_step"], first_b["last_step"]);
    let end_a = data(call(&mut s, "run", json!({"calc_id":"a","steps":100})));
    while !s.runtime("b").unwrap().state().terminal() {
        data(call(&mut s, "step", json!({"calc_id":"b"})));
    }
    let end_b = data(call(&mut s, "status", json!({"calc_id":"b"})));
    for field in [
        "state",
        "completed_maps",
        "iteration",
        "last_step",
        "energy_hartree",
    ] {
        assert_eq!(end_a[field], end_b[field], "{field}");
    }
    // Independent batch entry point, using the same parsed model and profile.
    let input = problem();
    let refs: Vec<_> = input
        .potentials
        .iter()
        .map(|(label, potential)| SpeciesPotential { label, potential })
        .collect();
    let delivered = dft_engine::refscf::api::solve(
        dft_engine::refstart::StartInput {
            cell: &input.cell,
            atoms: &input.atoms,
            potentials: &refs,
            options: &input.options,
        },
        input.request,
    )
    .unwrap();
    let dft_engine::refscf::api::Run::Scalar(run) = delivered.run else {
        panic!("scalar fixture")
    };
    assert_eq!(
        end_a["energy_hartree"].as_f64().unwrap().to_bits(),
        run.outcome.result.total_hartree().to_bits()
    );
}

#[test]
fn policy_after_completed_map_and_paused_step_preserve_boundary() {
    let mut s = server();
    let mut input = problem();
    input.options.max_iterations = 100;
    s.register_input("live", input).unwrap();
    handshake(&mut s);
    data(call(
        &mut s,
        "create_calculation",
        json!({"input_id":"live","calc_id":"c","mode":"supervised"}),
    ));
    loop {
        data(call(&mut s, "step", json!({"calc_id":"c"})));
        let r = s.runtime("c").unwrap();
        assert!(
            !r.state().terminal(),
            "fixture must expose a live safe boundary"
        );
        if r.snapshot().safe_point {
            break;
        }
    }
    let args = control(&s, "pause");
    data(call(&mut s, "pause", args));
    s.approve("c", "pause", approval()).unwrap();
    let before = s.runtime("c").unwrap().snapshot();
    data(call(&mut s, "step", json!({"calc_id":"c"})));
    assert_eq!(before, s.runtime("c").unwrap().snapshot());
    let mut args = control(&s, "policy");
    args["mixing_history"] = json!(3);
    args["threshold_ry"] = json!(1e-10);
    data(call(&mut s, "propose_policy_change", args));
    s.approve("c", "policy", approval()).unwrap();
    assert_eq!(s.runtime("c").unwrap().policy_revision(), 1);
    assert_eq!(
        s.runtime("c").unwrap().snapshot().completed_maps,
        before.completed_maps
    );
    let args = control(&s, "resume");
    data(call(&mut s, "resume", args));
    s.approve("c", "resume", approval()).unwrap();
    data(call(&mut s, "step", json!({"calc_id":"c"})));
    assert_eq!(
        s.runtime("c").unwrap().snapshot().completed_maps,
        before.completed_maps + 1
    );
}
