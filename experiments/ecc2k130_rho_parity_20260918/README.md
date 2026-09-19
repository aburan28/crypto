# ECC2K-130 toy-suite rho parity receipt

Cited from this directory. Protocol:
[`RESEARCH_ECC2K130_RHO_PARITY.md`](../../RESEARCH_ECC2K130_RHO_PARITY.md).

Host `ip-172-31-19-103`, 2026-09-18 (iterations 0–1) and 2026-09-19
(iterations 2–5). Unit: exclusive group operations. Not n=131.

**Iterations 0–2 under-charged collection and descent**: the `|F| − 1`
point additions a windowed probe makes before its lookups were not in
`G_IC`. Their gates are withdrawn (protocol §10). Their `summary.json`
files carry the counts (`collection_trials`, `window`,
`descent_trials`, `factor_base_size`) from which the corrected `G_IC`
is derived, and that derivation equals the iteration-3 rerun.

| iteration | what | cells | all-cases gate | class |
|---|---|---|---|---|
| 0 | folded-row rule + walked collection, two-pass table | 5 | reported unmet (α 1.05–1.31); corrected α 16.5–23.1 | engineering, superseded |
| 1 | one-pass table | 5 | reported **met** (α 0.76–0.92); corrected α 16.1–22.7 | engineering, superseded |
| 2 | twelve-rung ladder `2^{11}`–`2^{39}`, caps raised | 12 | reported unmet on 3 rungs; corrected α 8.4–612 | superseded |
| 3 | summand additions charged | 12 | unmet on every rung, α 8.4–612 | **accounting** |
| 4 | collection batch 8 | 12 | unmet, α 4.8–616 | engineering |
| 5 | batch 8 + `|F| = ⌈(2r)^{1/3}⌉` | 12 | unmet, α 4.8–218 | engineering |

Every iteration: all pairs verified `[d]G = Q` on both arms
(45/45 for 0–1, 108/108 for 2–5).

`table.txt` is the runner's per-pair print. `summary.json` is the
machine receipt (from iteration 3 it carries
`g_ic_before_summand_correction`, `collection_summand_additions`,
`descent_summand_additions`, and a `fit` block with the least-squares
slopes). `run.log` is the raw session; it is not a cost claim.
