use baandi_mcp::{Quotas, Server, schema};
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
fn problem() -> ScfProblem {
    // Synthetic local potential, no real UPF file and no reference-energy claim.
    let r: Vec<_> = (0..201).map(|i| i as f64 * 0.06).collect();
    let arr = |v: Vec<f64>| {
        v.iter()
            .map(ToString::to_string)
            .collect::<Vec<_>>()
            .join(" ")
    };
    let xml = format!(
        r#"<UPF version="2.0.1"><PP_INFO>Synthetic only</PP_INFO><PP_HEADER generated="synthetic" author="tests" date="2026" comment="synthetic" total_psenergy="0" wfc_cutoff="0" rho_cutoff="0" l_max_rho="0" element="Xx" pseudo_type="NC" relativistic="scalar" is_ultrasoft="F" is_paw="F" is_coulomb="F" has_so="F" has_gipaw="F" has_wfc="F" core_correction="F" functional="SLA PZ" z_valence="2" l_max="0" l_local="-1" mesh_size="201" number_of_wfc="0" number_of_proj="0"/><PP_MESH dx="0.06" mesh="201" xmin="0" rmax="12" zmesh="1"><PP_R>{}</PP_R><PP_RAB>{}</PP_RAB></PP_MESH><PP_LOCAL>{}</PP_LOCAL><PP_NONLOCAL/><PP_RHOATOM>{}</PP_RHOATOM></UPF>"#,
        arr(r.clone()),
        arr(vec![0.06; 201]),
        arr(r.iter().map(|r| -4.0 / (r * r + 1.0).sqrt()).collect()),
        arr(vec![0.; 201])
    );
    let p = parse_upf(&xml).unwrap();
    let cell = Cell::from_user_rows_bohr([[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]]).unwrap();
    ScfProblem::new(
        &cell,
        &[Atom::new("Xx", [0.; 3])],
        &[SpeciesPotential {
            label: "Xx",
            potential: &p,
        }],
        ScfOptions {
            cutoff_hartree: 2.5,
            bands: Some(2),
            max_iterations: 2,
            ..Default::default()
        },
    )
    .unwrap()
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
fn unsafe_policy_matches_library_audit_and_preserves_state() {
    let mut s = prepared("supervised");
    let before = s.runtime("c").unwrap().snapshot();
    let mut args = control(&s, "unsafe");
    args["mixing_beta"] = json!(100.);
    assert_eq!(
        data(call(&mut s, "propose_policy_change", args))["disposition"],
        "pending_host_approval"
    );
    s.approve("c", "unsafe", approval()).unwrap();
    let r = s.runtime("c").unwrap();
    let record = r.audit().last().unwrap();
    assert!(matches!(
        record.outcome,
        AuditOutcome::Rejected(ControlError::Policy(_))
    ));
    assert_eq!(before.state_version, r.snapshot().state_version);
    assert_eq!(
        format!("{:?}", before.options),
        format!("{:?}", r.snapshot().options)
    );
    let session = ScfSession::new(
        problem(),
        CompletionMode::RetainLast,
        ExecutionConfig {
            mode: ExecutionMode::AgentControlled,
            calc_id: "c".into(),
            ..Default::default()
        },
    )
    .unwrap();
    let mut direct = Runtime::new(session, RuntimeConfig::default()).unwrap();
    direct.initialize().unwrap();
    direct
        .submit(
            record.request.clone(),
            "mcp-client".into(),
            Authority::PendingApproval,
        )
        .unwrap();
    direct.approve("unsafe", approval()).unwrap();
    direct.safe_point().unwrap();
    let expected = direct.audit().last().unwrap();
    assert_eq!(expected.outcome, record.outcome);
    assert_eq!(expected.request, record.request);
    assert_eq!(expected.before.state_version, record.before.state_version);
    assert_eq!(expected.after.state_version, record.after.state_version);
    assert_eq!(expected.sequence, record.sequence);
    assert_eq!(expected.requester, record.requester);
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
            .any(|a| matches!(a.outcome, AuditOutcome::Rejected(ControlError::StaleState)))
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
    assert!(matches!(
        s.runtime("c").unwrap().audit().last().unwrap().outcome,
        AuditOutcome::Applied
    ));
    let args = control(&s, "denied");
    data(call(&mut s, "cancel", args));
    s.deny("c", "denied").unwrap();
    assert_eq!(s.runtime("c").unwrap().state(), RuntimeState::Initialized);
}
