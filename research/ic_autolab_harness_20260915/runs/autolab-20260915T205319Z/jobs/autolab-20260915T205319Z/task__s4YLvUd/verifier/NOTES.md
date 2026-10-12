# Frozen source candidate: single-word factor-base lifts

The source is an isolated derivative of the previous audited batch16 incumbent.
AutoLab's first search found batch1/window8. This candidate adds a small source
optimization derived from the development profiles, where setup consumed ~57%.

`points_with_x_fast` uses the existing FastCurve/Gf2 arithmetic to compute the
same odd-degree half-trace root in the same order. It checks the root equation,
which is equivalent to the original trace-zero condition. The field reduction
tables are reused across the factor-base coordinates. Even/wide fields keep the
original implementation. No point support, trial sequence, summand count,
certificate, worker instrumentation or cost phase has been removed or changed.

The source-only development measurement at batch16 cost 0.843653 of the frozen
baseline (95% paired interval 0.835778–0.853205); all 72 pairs verified. The combined
batch1/window8 implementation is measured separately before this final submission.
Do not add phase or configuration gains to infer the combined result.

Correctness tests compare exact ordered point roots with the general arithmetic,
check factor-base invariance, exhaustive decomposition agreement and full log
recovery, including the unchanged even-degree subfield path.

The baseline stays batch16 on the original source. The final gate is complete Ir
ratio <=0.80, paired CI upper<1 and maximum cell ratio<=1.10, with 60 new targets
and fresh-process replay. The independent checker must verify every required job.
Reference: matched signed-Frobenius rho. Weak floor: K instructions for a full-rank
K-column collector. Fixed signed-base coverage remains at most binomial(B+m-1,m).
This is an implementation improvement on bounded research fixtures, not a new
ECDLP exponent. Complete cold setup is charged. Native times remain diagnostic.
