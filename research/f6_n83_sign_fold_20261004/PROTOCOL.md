# n=83 F6 sign-folded pair-sum index experiment

Registered before execution on 2026-10-04. The K0 cofactor-projected factor
base is closed under point negation: its 4,054 usable points form 2,027
sign pairs. The hypothesis is that each nonexceptional pair sum and its
negative can share one stored representative and one target-query
denominator inversion. The candidate will build an exact x-keyed quotient
index with the positive and negative source pair, then test both target
residuals for each representative. It must reject construction if the input
point list is not distinct and sign-closed. Exceptional x=0 and infinity
cases must use the reference group law. Every returned witness must replay
in the full curve group.

Compare against `F6WidePairIndex` at commit
`2a0eb9b2e3fc6c24319754b2a4d77c2e01e1202a`. Freeze K0
`icv1-f2m83-tm6151469093347-debefd74`, its standard polynomial-subspace
cofactor-projected bases at dimensions 8, 10, and 12 (258, 1,048, 4,054
actual usable points), and public T001 (x=`355fb5df7a905f16921eb`,
y=`5900a390f42d290f1bbe`). The full index has 8,219,485 unordered pairs.
Run native release probes for both arms with the same compiler and resource
conditions. First prove exact candidate behavior against brute-force small
groups, including a planted witness and a miss; then run one small-base probe
per arm and full-base runs in baseline, candidate, candidate, baseline order.
Preserve all exits, stdout/stderr, source/binary hashes, index build time,
target-query time, sign-quotient count, peak RSS, and witness replay.

Retain the candidate only if it passes all correctness controls, lowers
full-base query time by at least 20% in the two-run median, has no more than
10% small-base query regression, no more than 2x full index-build cost, and
no more than baseline peak RSS. A failed gate means revert code but keep
the protocol, exact patch and all run receipts. An unisolated host makes
timings exploratory diagnostics only; this path's dimension-12 four-sum
coverage ceiling remains `4.662e-12`, so neither an F6 end-to-end nor a
one-target IC speedup follows from this experiment.
