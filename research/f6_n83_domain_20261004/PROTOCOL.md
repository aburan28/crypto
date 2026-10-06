# n83 F6 union-domain width correction

Registered before validation on 2026-10-04. PR #1366 found an exact
332,166-point Frobenius closure in the n83 K0 subgroup, but the existing
three-summand union SAT path currently collects only the low `u64` of
each 83-bit abscissa and its coordinate-domain trie shifts a `u64` by
bit positions above 63. This can exclude valid high-bit coordinates or
mislabel a result as a proof of absence. The hypothesis is that a
`u128` trie and complete two-word coordinate extraction correct that
specific domain predicate without changing any sub-64-bit result.

Frozen validation:

1. Keep the existing exhaustive small-domain trie test for widths 1–5.
2. Add an 83-bit test containing two allowed coordinates with identical
   low 64 bits but different high bits; require both SAT and a neighboring
   high-bit assignment UNSAT.
3. Run the affected unit tests and the n83 F6 geometry tests with the
   same source snapshot. Check rustfmt and diff integrity.

Only widths at most 128 are admitted by this predicate. A wider union
input must report `unsupported` before encoding. No n83 S4 SAT solve,
ordinary relation, CPU speedup, or end-to-end IC claim follows from a
passing width test. A later solver experiment must freeze target inputs,
budgets, and outcomes separately.
