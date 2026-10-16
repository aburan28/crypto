# Byte-exact nine-arm Linux replay sources

These files copy the five example sources at commit `7b82cb926c`, before
formatting the active examples for the PR's Rust lint gate. The frozen
`LINUX_REPLAY_SPEC_V1.json` and `HELDOUT_NINE_ARM_SPEC_V4.json` hashes refer
to these bytes. The Linux workflow restores them before it verifies source
hashes, compiles, and runs the 108 children. Keep these snapshots immutable;
a changed implementation needs a new spec and replay version.

The matching SHA-256 values, in workflow order, are:

| source | SHA-256 |
| --- | --- |
| `different_j_chi_resultant.rs.txt` | `c3cb3b06fd5d90708334d6b88ad878a780ddbd34979e8da88dfc363e41a74319` |
| `different_j_inverse_velu_oracle.rs.txt` | `185004244350c0a29306e232f9a13d3942cf6630923f2024cdb1bc02a4ddbcbc` |
| `different_j_chi_heldout_check.rs.txt` | `78966d5663e17aa43aeb2b1e3ac44ef605c64241db8fd88f0d79b6e97f8ef314` |
| `different_j_chi_heldout_summary.rs.txt` | `063c599c967fa5eb58a024b2b12c7e5133263868ccb2c862aba8040ced35c73e` |
| `different_j_chi_heldout_meter.rs.txt` | `3e478aa3edea464e97326cd4a4c6e34010a49bc18a3217ce5a809702bf0c2d7a` |
