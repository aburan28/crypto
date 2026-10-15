# Proof map: 21-byte N83 root pairs

Lean toolchain: `leanprover/lean4:v4.34.1`. The checked module is
`RootPair.lean`; run `lean RootPair.lean` from this directory and inspect
`proof_check.log`. Its `#print axioms` lines report only Lean's standard
`propext` and `Quot.sound` for the general arithmetic lemmas, and no axioms
for the closed capacity inequalities.

| Result | Lean declaration | Exact assumptions | Status |
| --- | --- | --- | --- |
| Concatenated roots decode exactly | `root_pair_roundtrip` | `root0 < 2^83`; `root1` arbitrary Nat | Checked |
| Two valid roots fit 21 bytes | `root_pair_fits_21_bytes` | Both roots `< 2^83` | Checked |
| K3240 complete two-key bound fits fixed `2^30` table | `k3240_root_keys_fit_fixed_2pow30_slots` | `stateBound` as defined | Checked |
| K3240 state ID fits 29 bits | `k3240_state_id_fits_29_bits` | `stateBound` as defined | Checked |
| K3240 arrays below 13.6 GiB | `k3240_payload_below_13point6_gibibytes` | 21 bytes per complete state and five bytes per `2^30` slot | Checked |

The source's byte-order code, symmetric S3 root existence for every retained
pair, and exact first-witness equality are checked by the Rust boundary and
six-representative tests and by the full cloud verifier. The bridge from the
Lean `stateBound` arithmetic to that executable traversal and its field
operations remains partial; no theorem here assumes it as a conclusion.
