# Indexed N83 root table proof map

| Result | Lean declaration | Hypotheses | Status |
| --- | --- | --- | --- |
| K2543 source-count maximum fits the 29-bit state identifier | `N83Indexed5.k2543_state_id_fits` | `states ≤ 83*K²`, K=2543 | Built, no axioms |
| K2544 maximum crosses the identifier bound | `N83Indexed5.k2544_state_id_crosses` | Numeric source-count maximum, K=2544 | Built, no axioms |
| The regular two-key request at K2543 fits the `2^31`-slot bucket | `N83Indexed5.k2543_regular_request_fits_bucket` | `keys ≤ 2*states`, existing twofold allocator | Built, no axioms |
| K2544 maximum crosses that bucket | `N83Indexed5.k2544_regular_request_crosses_bucket` | Numeric maximum and allocator request | Built, no axioms |
| Five-byte entries make that bucket 10 GiB | `N83Indexed5.ten_gibibytes_of_slots` | Five bytes per slot | Built, no axioms |

The Rust loop bound `states ≤ 83*K²`, the at-most-two-keys-per-state
insertion rule, the byte encoder, the exact-key callback and first-witness
order are not definitions in this Lean module. They are source and native
test obligations. These conditional capacity lemmas do not constitute an
end-to-end formalization of the Rust index. The pinned toolchain, command,
output and `#print axioms` checks are in `proof_check.json` and
`proof_check.log`. `lean4checker --fresh` was unavailable on this host.
