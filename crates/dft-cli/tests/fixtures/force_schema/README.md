# Format-only fixtures

These documents test schema branches, not scientific forces. The historical
v1 FD envelope comes from the stored CLI occupation fixture; the available
force payload is a deliberately composite format fixture. It must never be
used as a force reference for that envelope's model. Rejected fixtures mutate
one contract at a time, including nonzero numeric arrays under unavailable.

`python3 -I crates/dft-cli/tests/validate_force_schema.py` checks every JSON file.
The Rust CLI tests also validate actual emitted scalar, spin, shared-mu FD and
fixed-moment FD outputs, including cross-array lengths and force sums.
