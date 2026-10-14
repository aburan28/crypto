# N83 26-byte state proof map

| Result | Lean declaration | Hypotheses and status |
| --- | --- | --- |
| Label packing is injective and reversible | `N83State26.labels_roundtrip` | `left,right<4096`, `relative<128`; checked. |
| K2940 requests at most `2^30` slots | `N83State26.k2940_slot_request_fits` | Source-loop upper bound and `ceil(4*S/3)` request; checked. |
| K2940 state index fits 29 bits | `N83State26.k2940_state_indices_fit` | State count bounded by the source-loop expression; checked. |
| K2940 two-array payload below 13.7 GiB | `N83State26.k2940_payload_under_13point7_gibibytes` | 26 bytes/state and five bytes/slot; checked. Allocator overhead and actual RSS are empirical. |
| Complete source traversal preserves first witnesses | Open | Encoding/decoding is checked by source tests on real N83 roots; a Lean model of the Rust insertion loop and table probing remains to be constructed. |

`State26.lean` is checked with Lean 4.34.1 pinned by `lean-toolchain`.
`proof_check.log` retains the command output and `#print axioms` audit. The
label theorem uses only Lean's standard `propext`, `Classical.choice`, and
`Quot.sound` foundations; the three numeric theorems use no axioms.
No `sorry`, `admit`, native evaluation, or custom claim axiom appears.
