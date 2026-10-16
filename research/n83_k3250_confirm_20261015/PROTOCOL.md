# N83 K3250 vs K3300 confirmation panel (DRAFT, not frozen)

Status: **planned, no data collected.** Tracking PR; the frozen protocol and results land here.

## Question
`research/n83_k_neighborhood_results_20261015/RESULT.md` kept K3300: K3250 was nominally faster
(33.546 vs 34.239 ms, paired -0.692 ms, 95% block interval [-2.109, 0.717]) but the interval
includes zero. Is K3250 faster on a fresh, independently allocated panel?

## Design (to freeze before any data)
- Arms: root11 K3250 and K3300, same code and source hash as the K-neighborhood run.
- Fresh public subgroup targets (not workload `12e14f274f8c`), fresh process per arm, order
  balanced ABBA (see #1894 for the K3300/K3350 version of that balance).
- Statistic: paired K3250-minus-K3300 mean, contiguous-block bootstrap, same rule as before:
  promote K3250 only if the upper 95% endpoint is below 0.
- Every target independently replayed (scalar and four-point relation); failures and caps recorded.
- Report negative or inconclusive outcomes as they are; do not edit the committed
  K-neighborhood record, link to it.

## Out of scope
A K3250 win here is a stage comparison on one cloud instance, not an end-to-end claim.
