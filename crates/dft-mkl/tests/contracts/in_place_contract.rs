// Included as a unit test module so the exact pre-FFI guard runs without MKL.
#[test]
fn caller_buffer_length_contract() {
    for planned in [1, 60, 15625, 46656, 93312, 693] {
        assert!(crate::valid_in_place_len(planned, planned));
        for wrong in [0, planned - 1, planned + 1, usize::MAX] {
            assert!(!crate::valid_in_place_len(planned, wrong));
        }
    }
    assert!(!crate::valid_in_place_len(0, 0));
}
