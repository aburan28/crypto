# Complete volcano order measurements

All **46 curve cardinalities** and **9200 native point orders** match the independent PARI results. Every P and Q has order **18014494963342201**. Reference order reduction uses **12381 field multiplications including squares per point**. Curve-count CPU time spans **164–218 ms** under L0.

[Inspected run report](RUN_REPORT.pdf) · [Inspected results report](RESULTS_REPORT.pdf) · [Per-curve min/mean/max table](per_curve.csv) · [All per-problem summaries](per_problem.json) · [Definitions and limits](../../ORDER_PROTOCOL.md)

Point-order costs exclude precomputed cardinality/factorization; they are not solver field counters. GP timers have 1 ms resolution, and default internal threads=5. No fastest or hardness claim. Full-group invariants may differ, while the benchmark point orders remain the same.
