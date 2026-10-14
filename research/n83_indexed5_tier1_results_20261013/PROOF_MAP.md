# Five-byte N83 root-table capacity proof map

The regular indexed5 request allocates the next power-of-two slot count at
least `4 * states`: at most two candidate keys per state, then a twofold
table request. The source-loop count is bounded by `83 * K * K` states.
Under these two source obligations, K1798 is the last prefix size whose
regular root table fits the `2^30`-slot, 5 GiB bucket. K1799 moves to the
`2^31`-slot, 10 GiB bucket. This identifies a sharp allocation boundary;
the measured best factor base still depends on target work and full-process
resource use.

| Mathematical statement | Lean 4 declaration | Exact hypothesis and status |
| --- | --- | --- |
| The numeric source bound grows with K | `N83Indexed5Tier1.state_bound_monotone` | Natural `a ≤ b`; built with no axioms |
| K1798 fits `2^30` slots | `N83Indexed5Tier1.k1798_fits_five_gibibyte_bucket` | `4 * (83 * 1798 * 1798) ≤ 2^30`; built with no axioms |
| Every smaller nonnegative K fits that bucket | `N83Indexed5Tier1.all_smaller_k_fit_bucket` | Natural `K ≤ 1798`; built with no axioms |
| K1799 crosses the bucket | `N83Indexed5Tier1.k1799_crosses_five_gibibyte_bucket` | `2^30 < 4 * (83 * 1799 * 1799)`; built with no axioms |
| Every larger K crosses the bucket | `N83Indexed5Tier1.all_larger_k_cross_bucket` | Natural `1799 ≤ K`; built with no axioms |
| Byte size of the two buckets | `N83Indexed5Tier1.five_gibibytes_of_slots`, `N83Indexed5Tier1.ten_gibibytes_of_slots` | Five bytes per slot; built with no axioms |

The written proof applies monotonicity of multiplication on natural numbers
to the quadratic state bound. Direct integer arithmetic establishes the
K1798 and K1799 inequalities, and transitivity extends them to every K on
the corresponding side. Multiplying the slot counts by five gives the two
byte capacities.

The Lean module formalizes the integer implication. It does not define the
Rust source loop, prove its `83 * K * K` state bound, or model the allocator's
power-of-two rounding and actual resident memory. Those remain source and
runtime obligations; no absolute performance optimum follows from the
capacity theorem alone. The pinned toolchain is `leanprover/lean4:v4.34.1`.
`lean Capacity16G.lean` built on 2026-10-13 with Lean commit
`5045d0056413266e57c625dcd7c365b10e377c52`, and every declaration's
`#print axioms` output in `proof_check.log` reported no axioms.
`lean4checker --fresh` was unavailable on this host.
