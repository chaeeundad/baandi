use super::*;
use dft_engine::DavidsonPreconditioner as Kind;
#[test]
fn preconditioner_wire_and_config_contract() {
    for kind in [Kind::Default, Kind::Tpa, Kind::RitzShift] {
        let o = ScfOptions {
            preconditioner: kind,
            ..ScfOptions::strict_reference()
        };
        let mut bytes = Vec::new();
        let json = serde_json::to_value(o.encode(&mut bytes).unwrap()).unwrap();
        assert_eq!(json.get("preconditioner").is_none(), kind == Kind::Default);
        let decoded =
            ScfOptions::decode(serde_json::from_value(json).unwrap(), &bytes, &mut 0).unwrap();
        assert_eq!(decoded.preconditioner, kind);
        assert!(
            toml::from_str::<crate::config::Scf>(&format!("preconditioner = {:?}", kind.name()))
                .is_err()
        );
    }
    assert!(Kind::decode("unknown".into(), &[], &mut 0).is_err());
}
