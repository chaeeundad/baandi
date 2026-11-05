//! External licensed data only; never embed real radial samples.
#[test]
#[ignore = "requires external files through DFT_NC_DIR"]
fn reader_diagnostics() {
    let dir = std::path::PathBuf::from(std::env::var_os("DFT_NC_DIR").expect("set DFT_NC_DIR"));
    for element in ["O", "Ni", "Al", "Mg", "Si"] {
        let path = dir.join(format!("{element}.upf"));
        println!("{element}: vendor={:?}", upf::from_file(&path).err());
        let parsed =
            dft_upf::from_file(&path).expect("real NC adapter must load every selected file");
        assert_eq!(parsed.projectors().len(), parsed.raw().nonlocal.betas.len());
        assert_eq!(parsed.core_density().is_some(), element != "Mg");
        println!(
            "{element}: adapter=OK projectors={} wavefunctions={} mesh_points={} nlcc={}",
            parsed.projectors().len(),
            parsed.raw().header.number_of_wfc,
            parsed.mesh().r().len(),
            parsed.core_density().is_some()
        );
    }
}
