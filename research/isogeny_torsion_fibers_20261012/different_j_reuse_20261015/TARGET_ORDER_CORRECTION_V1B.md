# Pre-run correction: public target order must not encode the answers

The committed [V1 protocol](PROTOCOL_V1.md) proposed ordering targets
by `(13+37i) mod 61`. Because both parameters and the rule were public,
the row position itself would disclose each answer scalar to a solver.
No batch target file, batch solver or measured run was produced under
that rule. Keep V1 as the failed preregistration; it is superseded here
before implementation and before any batch measurement.

The corrected public order uses the same 60 nonidentity source subgroup
points, but sorts them by the hexadecimal SHA-256 of the UTF-8 string
`different-j-reuse-target-order-v1:13:<source point key>`, with the
numeric source point key as a tie-breaker. The seed `13` and hash rule
are public, yet they do not disclose the point's discrete logarithm.
The native freeze process computes all 60 points, sorts them, then
writes only source/destination point pairs to the public solver file.
It writes the scalar correspondence to a separate checker file. The
solver must reject an `answer_scalar` or `scalars` field anywhere in its
input and must not read the checker file. The independent checker
reconstructs the hash order using its own point-key formula, verifies
all 60 distinct pairs and the private answer correspondence, and
replays every solver output.

The [V1B machine spec](RUN_SPEC_V1B.json) replaces only V1's
answer-bearing target-order rule. Arms, prefixes, cap, resource
ceilings, cost custody and stop conditions remain as in V1. A
hash-pinned V2 spec will be committed after implementation and before
the first batch run. Sorting by point hashes still permits a deliberate
brute-force lookup in this order-61 group; the solver process is audited
for the absence of such a lookup. This is a public synthetic
known-answer fixture, not a secrecy claim.
