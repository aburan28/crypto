# n=83 K_1 standard-subspace base inventory

Registered before running the inventory on 2026-10-04. This is a structural
factor-base inventory, not a timing comparison or an IC candidate run.

Use the pinned `KoblitzCurve::known_n83_k1` model and the existing
`build_standard_subspace_factor_base` followed by public cofactor
projection. Enumerate exactly dimensions 8, 10 and 12, in that order.
For each, record the parent point count, distinct nonidentity projected
subgroup-usable point count, signed column count, and the exact counting
ceiling `binomial(B+3,4)/r` for uniform subgroup targets. If a construction
fails or runs out of memory, retain the failure. The ceiling is an upper
bound, not an estimate of natural relation yield. No target is solved here.
