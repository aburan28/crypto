# Bounded planted-witness continuation

Registered after the five root pairs and before any deeper search on
2026-10-05. The root gate passed: all ten runs returned `exhausted`
with identical exact counters, while delayed children used much less
memory. This continuation asks whether that scheduling change lets the
same planted n83 three-summand block produce an actual verified witness
within a small, fixed search budget. It does not estimate natural
relation yield.

Use the exact K0 curve, dimension-12 source base, planted indices
`[0,2,4]`, two-word table and solver source from the root candidate.
Add a native example that calls `wide_groebner_decompose` with
`m=3`, `node_budget=128`, and `RAYON_NUM_THREADS=1`. This thread count
makes search order and the budget deterministic; its wall time is only
a construction diagnostic, never compared with the multi-threaded
root baseline. Apply `gtimeout -k 5s 120s` and capture one complete JSON
result, exit status, stderr, and observed peak RSS. Treat an exhausted
search, timeout, OOM, or unverified output as failure of this gate and
preserve it.

Acceptance requires a solver-returned witness whose source points sum
to the planted source target, exact subgroup projection `R=4S`, and
the source-to-subgroup four-torsion bridge checked independently in the
existing unit test. This proves only a planted three-summand block, not
an eight-summand ordinary-query relation or a complete IC candidate.
