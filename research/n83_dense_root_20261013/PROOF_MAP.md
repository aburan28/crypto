# Dense packed-root capacity: mathematical proof map

For degree 83 and `K` selected Frobenius orbit representatives, the current
wide producer loops over each of 83 relative rotations for each ordered pair
of representatives with `left <= right`. A retained regular state is therefore
bounded by `83*K^2`. The packed root table requests either `2*states` entries
(regular mode) or `states` entries (dense mode), and its constructor allocates
the least power of two at least twice the request. Packed slots are 16 bytes.

| Result | Lean declaration | Exact hypotheses | Check |
|---|---|---|---|
| K1271 dense request fits `2^28` slots | `k1271_requested_slots_fit_four_gib` | `states : Nat`, `states <= 83*1271^2` | Built, no axioms |
| K1272 maximum state count crosses `2^28` slots | `k1272_maximum_crosses_four_gib` | Numeric maximum, not actual occupancy | Built, no axioms |
| K1798 dense request fits `2^29` slots | `k1798_requested_slots_fit_eight_gib` | `states : Nat`, `states <= 83*1798^2` | Built, no axioms |
| K1799 maximum state count crosses `2^29` slots | `k1799_maximum_crosses_eight_gib` | Numeric maximum, not actual occupancy | Built, no axioms |
| Byte sizes of the two allocated slot buckets | `four_gib_packed_slot_bytes`, `eight_gib_packed_slot_bytes` | `16` bytes per packed slot | Built, no axioms |

Toolchain: `leanprover/lean4:v4.34.1` from `lean-toolchain`; see
`proof_check.log` and `proof_check.json` for the exact build and `#print axioms`
output. `lean4checker --fresh` is unavailable here. The Lean module proves
the arithmetic implications. A mechanized Lean model of the Rust state loop,
its filtering, and the table allocator remains open; source inspection and
measured `regular_states`/`root_table_allocated_bytes` connect the theorem to
the executable for this experiment. The crossover statements concern the
maximum bound and do not assert that a particular measured run crosses it.

## Proposed 75% load rule (next experiment)

`Load75Capacity.lean` formalizes the numeric bound for a proposed request
`ceil(2*states/3)`, represented exactly as `(2*states+2)/3` in `Nat` arithmetic,
before the existing constructor doubles the request and rounds to a power of
two. Under `states <= 83*K^2`, `k2202_request_fits_eight_gib` proves the
pre-rounding request fits `2^29` slots at K2202;
`k2203_maximum_crosses_eight_gib` proves the maximum at K2203 crosses it.
These declarations were built with the same pinned Lean toolchain and audited
with `#print axioms`; the first depends on Lean's standard `propext` foundation
through `Nat.div_le_div_right`, and the second has no axioms. See
`load75_proof_check.json` and `load75_proof_check.log`. The proposed request
rule is not in the measured Rust binary, so this is a conditional capacity
result and no K2202 performance measurement.
