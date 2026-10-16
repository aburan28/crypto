# IC improvement round 4: results (2026-10-16)

One follow-up to [`../ic_improve3_20261015/README.md`](../ic_improve3_20261015/README.md), run under [`RULES.md`](RULES.md).

- Each candidate is measured against the incumbent component it replaces (round 3's R2).
- Comparisons are paired, with training (curves 5-7) and held-out (curves 0-4) splits, q = 2^12..2^22.
- Every logarithm and every root set is verified.

| # | Candidate vs incumbent | Verdict | Key numbers | Landed |
|---|---|---|---|---|
| [01](01-genus3-relations/README.md) | Genus-3 IC relation phase vs round 3's R2 (`r2:lpsgebw2`). **C1:** field-level splitting test (Lucas test of the Cardano resolvent's cubic character), identical relations. **C2:** 32 lockstep walks sharing one inversion. **C3:** 70% oversampling + clique removal + SGE cx 32. **C4:** C1+C2+C3 with kappa re-tuned. | C1 **ADVANCES** (exponent, constant, primitive); C2 ADVANCES (constant); C3 TRADES; C4 TRADES (constant ADVANCES, exponent regresses) | <ul><li>C1: whole 0.756x at 2^22, exponent Delta -0.015 (12..22), -0.034 (16..22); local exponent 21-22 1.364 vs 1.408; per-trial log slope 9.75 to 2.71 units/bit.</li><li>C2: 0.959x; C1+C2: 0.715x, local 1.356.</li><li>C3: 0.979x (LA 0.78x, relation 1.06x), gain shrinks with q.</li><li>C4: 0.690x at 2^22 (cheapest), but its steep kappa schedule raises the exponent (+0.0045 over 16..22).</li><li>Gate 0 41/41; every log verified.</li></ul> | example only |

## What changed in the picture

- **The per-trial log factor is mostly gone.** C1 replaces the cubic-polynomial exponentiation by a Lucas ladder in
  F_q; the log-dependent part of the cost per trial falls from about 9.8 to 2.7 units per bit, and the fitted
  exponent falls by 0.015-0.034. This was the relation-phase term of round 3's decomposition.
- **Distance remaining to 4/3:** fitted over 16..22, 1.368 (C1+C2) and 1.374 (C1) against 1.409; local 21-22 1.356
  and 1.364 against 1.408.
- **Handing LA a better matrix helps the constant, not the exponent.** Oversampling and clique removal cut LA by
  22-44%, most at small q; the saving shrinks with q, so the LA excess (round 3's largest term) remains.
- **Re-tuning kappa on training again produced a steeper schedule** and a higher exponent, as in rounds 2 and 3.

## Caveats

- **Host:** shared 4-vCPU container, isolation L0; a container restart killed three training jobs (rerun).
- **Measured unit:** operation counts; wall times descriptive.
- **Scale:** 2^23 was not run: 2^22 jobs peaked at about 3.9 GB of the 4.5 GB cap; memory is the limit.
- **Conductor:** the control plane was unreachable.
