# Pre-run correction: payload dimensions 3 and 4 require a different factor base

The preregistered extension proposed applying the archived
`CyclotomicFactorSpace(83,s)` to `s=3,4`. Preflight on the archived source
raised `ValueError: s must properly divide p-1`; indeed `83-1=82` has no
divisor 3 or 4. Preserve this failed setup in the evidence; its original
quotient representation remains valid only for `s=2` among the small sizes.
Do not silently substitute a new factor base and call it the same workload.

For a separate, clearly labeled *plain normal-basis subspace* sweep, fix
`x(payload) = sum_{j=0}^{s-1} bit_j(payload) * beta^(2^j)` for `s=2,3,4`,
where `beta` is the first normal generator selected by the archived
`NormalCoordinates` from each seed. These are nested `F_2` subspaces on each
seed and give a genuine controlled size sweep. Each of the three point maps
uses its own Frobenius phase `(0,1,2)`. Retain the same curve, modulus, prime
subgroup, cofactor, generator selection, two frozen seeds and S4 direct
Boolean descent as the original protocol. At each (s,seed), compare the first
proper lexicographic planted target and first random subgroup target using
`Random(seed*1000+830+s)`. Run FES, Z3, native-XOR CryptoMiniSat, the pinned
Boolean Macaulay/F4-style prototype and the pinned F5B prototype for each
input, including setup failure. Use the same 25 s build, 5 s solve, 70 s
process, and 1 GiB address-space limits. Record complete roots and verify
each group relation independently as before; FES checks all 64/512/4096
assignments for sizes 2/3/4, respectively.

The s=2 **quotient** cases from the original run and the native-XOR addition
under `PROTOCOL_EXTENSION.md` remain in their own tables. Every plain
subspace result has `factor_base_policy=normal_subspace` and goes into a
separate immutable results directory; never equate its s=2 curve inputs or
wall-clock figures with the s=2 quotient cases. If this normal-subspace sweep
exceeds its limits, report the bounded failure and stop there. No full-DLP
performance or m=131 transfer is established by this fixed-phase comparison.
