# n=83 K_0 confidence-gate base inventory

Registered before the first K_0 inventory run on 2026-10-04. The pinned
curve and subgroup are those of `gate-m83-T001.json`, distinct from the
K_1 stage timing fixture. Both subgroup orders have the checked Sage
primality receipt in this directory. This is a structural inventory, not a
target solve or CPU timing comparison.

Construct the K_0 curve through `known_n83_k0`, verify its registry slug,
then use the standard polynomial-subspace base followed by public cofactor
projection for dimensions 8, 10 and 12, in that order. Record the curve
point count, distinct nonidentity usable projected point count, signed
column count, and `binomial(B+3,4)` for each. Retain failures and partial
rows. The ratio to the frozen 81-bit subgroup order is an exact ceiling on
four-summand coverage for a uniform subgroup target, not a yield estimate.
