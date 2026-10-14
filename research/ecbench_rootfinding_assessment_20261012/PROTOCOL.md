# Root-finding point decomposition in ecbench

Frozen before execution on 2026-10-12 (America/Los_Angeles).

## Question and boundary

For one previously unseen public target on each registered curve, compare the counted cold complete-solve cost of exact two-summand index calculus using direct subtraction or a reusable pair-sum table against matched rho. The prime arm is the algorithmic search in Algorithm 1 of `rootfinding_decomposition.pdf`; the Koblitz table arm is the strategy in Algorithm 2. The ecbench implementations use their own arithmetic and are not a timing reproduction of the paper's optimized code. The paper input has SHA-256 `743d3d3e7e3aab6f980813c0944c1d2cbacc906e075e7283de0abc7b923b7c01`.

The reference is `rho.negation` on the prime curve and `rho.signed_frobenius_strong` on the Koblitz curve. The generic floor is recorded by ecbench for each workload. Both IC arms use the same factor base and full-rank checked linear algebra. Each curve has eight distinct public targets and three measured rounds after one warmup, with identical target and algorithm seeds across arms.

## Predictions and interpretation

The direct subtraction oracle has no pair-table setup; the binary pair table should reduce target PDP work but may not recover its setup cost in a cold one-target solve. A complete run requires ecbench scalar verification and an exact replay audit. A wrong answer, exhaustion, timeout, or missing matched arm leaves that workload incomplete. Counts, phases, factor-base size, setup, unpriced work, and memory are reported separately. Local macOS runs earn L0, so wall time is descriptive only. These two sizes do not support a scaling law or a P-192–P-521 conclusion.

## Execution boundary

The exact specs are `specs/prime20.json` and `specs/koblitz17.json`. Keep the source, specs, target seed, arms, full-rank mode, trial cap, and timeout unchanged after seeing outcomes. Run each spec once into a new session directory; retain all statuses. Audit with `ecbench verify --replay 12`, then compare paired arms and tabulate. No claim promotion follows directly from a run.
