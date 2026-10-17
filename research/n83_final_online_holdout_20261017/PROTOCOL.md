# Twelve sealed N83 points for the final online factor-base comparison

This fixture fixes twelve previously unused public subgroup points before
either pending Tier 1 factor-base panel has produced a selection. It reserves
an independent comparison of the smallest cold base and the strongest online
challengers. No timing, candidate choice, or resource result is inferred from
the point coordinates.

The exact curve is `EC1N83Ckb1a1haacb69ee9f37`: polynomial-basis
`GF(2^83)` with modulus integer `9671406556917033397649557`, equation
`y² + xy = x³ + x² + 1`, prime subgroup order `8569786107849059`, and
cofactor `1128547018`. The generator and subgroup checks are inherited
through the frozen fixture lineage. For index `i = 0..11`, the scalar law is
`1 + (little_endian(SHA256(domain || u32_le(i))) mod (r - 1))`, with domain
`n83-final-online-holdout-twelve-20261017-v1`. Only the public points are
producer inputs; [answers.jsonl](answers.jsonl) is an independent replay
sidecar. `generate.py --verify-sealed` regenerates the exact files, checks
subgroup membership, and proves non-overlap with the precision, finalist,
five-arm, K992 and earlier N83 fixtures represented by its parent chain.

| Frozen item | SHA-256 |
| --- | --- |
| [public targets](targets.jsonl) | `4d167f8b2b85bcf6c599494daf3b3a61210a9dc8c0d65ad3ebcdf07cd6d5c4c7` |
| [answer sidecar](answers.jsonl) | `f94a313ea6e7c414a4fc85c3ec55e384f6c6be14febbe3304d4b9dde1e54e484` |
| parent K992 generator | `514329656510e88d3316b688eef4d15104c64a5be7edb42f56ace269461d1dec` |

The later dispatch must choose and hash exact `IC1` manifests before any
holdout measurement. Its admissible choices come only from the pending
five-arm online panel and K992/K1024/K1107 panel, plus the fixed K992 and
K1107 controls. It must freeze the source and executed binaries, arm order,
caps, Tier 1 resource envelope, verifier, selection rule, and matched rho
policy in a separate dispatch record. Candidate selection must use online
time for one independently verified public point, paired with rho on that
same point; cold process cost is a separate factor-base objective. Run
points 0–5 and 6–11 on independent Tier 1 instances, retain every failed
or noisy pair, and pool strata only after each passes its frozen tier gate.
No point from this fixture is used to choose the candidate list or tune an
implementation before dispatch.
