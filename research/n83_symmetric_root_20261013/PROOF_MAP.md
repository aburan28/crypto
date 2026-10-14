# S3 pair-symmetry proof map

The index's new optional state quotient retains `(left,right,relative)` only
when `left < right`, or when `left == right` and `relative <= floor(83/2)`.
Its mathematical basis is the S3 input exchange plus Frobenius transport.
The exact index claim also requires a bridge to the source's normal-basis
rotation, canonical-key function, root ordering, and first-inserted witness.

| Result | Lean declaration | Exact hypotheses | Status |
| --- | --- | --- | --- |
| S3 input exchange | `N83SymmetricRoot.s3_swap` | Any commutative semiring, all `x,y,z` | Checked |
| Exchanged roots remain roots | `N83SymmetricRoot.s3_swap_root` | Same structure, `s3 x y z = 0` | Checked |
| S3 commutes with an endomorphism | `N83SymmetricRoot.s3_map` | Function preserving addition and multiplication | Checked |
| Roots transport through an endomorphism | `N83SymmetricRoot.s3_map_root` | Above plus zero preservation | Checked |
| K2543 fits the five-GiB root-slot bucket under the quotient state bound | `N83SymmetricRoot.k2543_fits_five_gibibyte_bucket` | `quotientStateBound k = 83 * (k*(k-1)/2) + 42*k`; four slots per state | Checked, conditional on source-count bridge |
| K2543 state indices fit the five-byte slot's 29-bit field | `N83SymmetricRoot.k2543_state_index_fits_twenty_nine_bits` | Same state bound, plus one for the nonzero slot marker | Checked, conditional on source-count bridge |
| K2544 crosses that bucket | `N83SymmetricRoot.k2544_crosses_five_gibibyte_bucket` | Same bound and request rule | Checked, conditional on source-count bridge |
| K2543 state and slot payload fits 12 GiB | `N83SymmetricRoot.k2543_state_and_slot_payload_under_twelve_gibibytes` | 28 bytes per retained state and 5 bytes per `2^30` slots | Checked, conditional on source-count and allocator bridges |
| Reversed N83 index states have the same canonical-root orbit | Open | Identify source `basis.rotate` with the 83-power Frobenius automorphism; prove root-order coverage and source canonicalization | Partial |
| Retaining the earlier state preserves every first exact witness | Open | Formalize source loop order, canonical-root equality, and indexed-table insertion semantics | Partial |

The checked polynomial is
`(xy+xz+yz)*(xy+xz+yz)+xyz`, the characteristic-two S3 equation used by
the root solver. The Lean result is algebraically general and does not
assume the source's field implementation or its canonical encoding. The
source-level unit test checks these missing links for a six-representative
N83 fixture; larger replay and runtime comparisons are separate gates.
The 12-GiB calculation covers the two array payloads only; actual RSS
includes the base, shifted coordinates, hash maps, process code, and allocator
overhead. It is a capacity gate for the next measured run, not an RSS result.

Toolchain: `leanprover/lean4:v4.34.1`, using Lean's bundled
`Init.Grind.Ring.Basic` only. Command: `lean S3Symmetry.lean` from this
directory. `#print axioms` reports only `propext` for the four checked
declarations; the four arithmetic declarations report no axioms. No `sorry`,
`admit`, custom claim axiom, or native proof shortcut
is used. `lean4checker --fresh` was unavailable on the development host.
