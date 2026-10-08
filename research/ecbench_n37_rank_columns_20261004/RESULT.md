# Decision: K16 lowers counted cold work, while rho still leads

The preregistered six-base n37 sweep selected **K16** by minimum total cold
charged group-addition equivalents (GAE) across 40 verified measured runs.
It built 1,184 actual usable factor-base points before signed-Frobenius
folding, 16 effective columns, full rank and all 16 verified column logs.
The decision is `COUNTED_ENGINEERING_LEAD` against the fixed K42 base:
K16/K42 is 0.34898, with a workload-clustered 95% interval
[0.34334, 0.35866]. This is a 65.1% reduction in the **counted lower
bound**, exceeding the frozen 20% gate. The comparator used the same eight
public points, algorithm seeds, five measured repetitions and limits.

| One-target arm | Actual base points | Folded columns | Verified measured runs | Mean cold S lower bound | Counted IC/rho diagnostic |
|:--|--:|--:|--:|--:|--:|
| Strong signed-Frobenius rho | — | — | 40/40 | 0.419 | 1.000 |
| IC K4 | 296 | 4 | 40/40 | 4.199 | 10.013 |
| IC K8 | 592 | 8 | 40/40 | 2.133 | 5.087 |
| IC K12 | 888 | 12 | 40/40 | 2.201 | 5.249 |
| **IC K16** | **1,184** | **16** | **40/40** | **2.111** | **5.034** |
| IC K24 | 1,776 | 24 | 40/40 | 2.386 | 5.691 |
| IC K42 | 3,108 | 42 | 40/40 | 6.050 | 14.426 |
| Identical K42 control | 3,108 | 42 | 40/40 | 6.050 | 14.426 |

The K16 `IC1` identity is
`IC1N37Ckb0fb1184PDP3mitmfrobeniuscountedRCsampleLAgaussTDpdpISO0h44f5af6dc772`.
The [candidate claims](candidate_claims) give the exact manifests and base
digests for all six arms; these are descriptive claims because the primary
wall-time isolation gate is not met. The rho and every IC answer passed
scalar replay, and the Linux x86-64 independent auditor reproduced all 320
measured runs exactly. All warm-ups and failure rows are retained; there
were no failures or exhausted searches in this session. The K42 A/A counted
ratio is exactly 1.000.

The selected K16 mean lower-bound cost is 31,312.49 GAE of reusable
preparation plus 748.32 GAE of target-dependent work, versus K42's
91,712.76 plus 156.89. Smaller bases reduce the quadratic folded-table
cost but can make rank collection and target decomposition harder. K4
averaged 48.675 target attempts; K8 averaged 7.375; K16 averaged 1.3;
K24 and K42 each averaged 1. K8 is close to K16: the paired K16/K8
counted ratio is 0.98967, 95% interval [0.92717, 1.05694]. The frozen
minimum-sum rule selects K16, but this sample does not distinguish those
two sizes statistically. A preregistered untouched confirmation should
therefore compare K8 and K16 rather than treat the exact K as universal.

The same-point K16/rho *counted lower-bound ratio* is 5.0344, with a
workload-clustered interval [4.3462, 5.8664]. It is a diagnostic ratio
between two incomplete cost estimates, **not a lower bound on true
IC/rho cost or a verified speedup**. Field arithmetic, hashing, allocation
and modular combination remain unpriced for both arms. Every producer
interval is macOS L0, so online wall-time speedup is `null`; the isolated
CPU receipt required by `AGENTS.md` is absent. The online five-phase clocks,
all cold phase charges, memory and raw statuses remain in the sealed
[session](sessions/mac_arm64_l0_01). The reported online GAE is a counted
stage diagnostic, not a substitute for the primary online wall metric.

This narrows one implementation choice at n37. It does not show that
factor-base descent, degree-263 isogeny transport, F4/F5/SAT PDP, or
ECC2K-130 at n131 benefits from K16. The next decision gate is an untouched
K8/K16 target panel with a complete native-work price and auditable host
isolation, then a larger-field scaling check before transfer to n131.
The [protocol](PROTOCOL.md), [spec](SPEC.json), [plan](PLAN.json),
[raw session](sessions/mac_arm64_l0_01), [local audit](AUDIT.json),
[independent receipt](independent_validation/RECEIPT.json), [saved paired
comparisons](sessions/mac_arm64_l0_01/comparisons), [decision data](RESULT.json)
and [native analyzer](../../examples/n37_rank_columns_analyze.rs) fix the
inputs and checks behind this conclusion.
