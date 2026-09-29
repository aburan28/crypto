# Symbolic O-aware width gate: current n131 path refused; n13 μ4 point map passes

The [pre-outcome draft PR #796](https://github.com/aburan28/crypto/pull/796)
froze the final producer, independent verifier, bounded runner, protocol,
inputs and source pins at `7e9f3e3` (FROZEN SHA-256
`44f713e97398a9f3e8bcaad3cf854f207eca302c2798eb8049b0ba0d86b45d4c`).
Its exact-head hash-only and syntax CI passed before the first selected run.
The one cold run passed without a retry. The first raw
[receipt](evidence/receipt.json) is 2,779 bytes, SHA-256
`8dfd7c98397db3fa0be5f7b946ae03e831c1942577c758aabb3d67818320fd30`.
The independent archive replay recomputed the full toy source and μ4 chart
sets from the producer receipt and passed.

| Frozen check | First-run result |
|:--|:--|
| Complete K0 finite point set | 8,011 points; canonical set SHA-256 `3a52f45a2e92c694a70f571815b9fb96c460367bc0794563cc692149b5513652` |
| Complete split-μ4 `X2=1` chart | 8,011 points; canonical set SHA-256 `2f2c4915119723a4889bfb16455a36148036bd7b375c04f96e323cd06c7d0753` |
| Projective total, including unique O chart | 8,012 points; exact sets and inverse roundtrips agree |
| Negative/exception controls | Corrupt `(0,1)` image, fake `X2=0` point, O affine-inverse guard, and `(1,0)` order-4 chain all PASS |
| Unequal n131 m9 width | 131 factor bits + 917 intermediate x bits = 1,048 raw affine-chain bits; `NOT_ADMITTED` |
| Unequal n131 m10 width | 131 + 1,048 = 1,179 bits; `NOT_ADMITTED` |
| Balanced n131 m10 control | 130 + 1,048 = 1,178 bits; `NOT_ADMITTED` |

The normal-form check establishes an exact full-point map on this n13 toy:
`(x,y) ↦ (x²,x²+y,1,x²+x+y)` and `O ↦ (1,1,0,1)`. The independent
verifier used polynomial-product/reduction field arithmetic, Euclid
inversion, complete square-root and Artin-Schreier root tables. The
producer instead used shift/reduce multiplication, Fermat inversion,
trace and half-trace. Both independently enumerated the source and the
normal-form chart and matched every point-set digest. This checks the
conversion and its O exception; no projective addition-law selector or
Boolean circuit was implemented. Kohel's four bidegree-(2,2) law charts
still have exceptional divisors. The separate binary Edwards
complete-formula route likewise needs an exact conversion and circuit
validation for this curve.

The `NOT_ADMITTED` result is a current **interface and evidence** decision.
`Gf2` is limited to degree at most 63; the symbolic field coefficients,
Boolean monomials and SAT wrapper/model assignment use `u64` masks. Both
n131 unequal arms exceed 64 bits in factor inputs alone, before the
full-point y, O, slope and selector bits. No multiword symbolic full-point
exporter, model-to-point map, 32-target toy semantic equivalence for that
exporter, externally certified negative branch panel, or bounded n131
export receipt exists here. The widths do not estimate SAT hardness or
rule out a new encoding. The independent n13 map does not transfer a
solver result to n131.

The next implementation gate is the solver-neutral XOR/AND DAG and
multiword GF(2^n) bit-vector layer described in [DESIGN.md](DESIGN.md).
It must implement the four disjoint full-point addition cases, reject
noncanonical O/off-curve/branch-inconsistent models, and match all 32
fixed #781/#785 Q+T statuses with independent model replay and negative
certificates before a capped n131 **export-only** attempt. Relation
support, rank, collection cost and any rho comparison remain unset.

The measured run's total wrapper wall was 3.138 s, below the 60 s cap.
Producer/independent-verifier child walls were 2.643/0.487 s, CPU
2.113/0.420 s, and cumulative child `ru_maxrss` upper readings
29,343,744/29,605,888 bytes, below 256 MiB. Both exits were zero.
The exact argv, UTC boundaries, stdout/stderr and hashes for every
artifact are in the raw receipt. These are toy map/width-audit costs,
not export, solver or attack costs. The retained
[pre-outcome CI failure](FAILURES.md) concerned an absent empty evidence
directory in checkout; it produced no selected outcome and was corrected
before this first run.
