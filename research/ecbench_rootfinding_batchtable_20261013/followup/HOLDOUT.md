# Independent factor-base confirmation

The frozen sweep `ECBS1hdabf729ffad4` selected B48 by mean cold counted
GAE: B48 7,290.083, B40 7,597.708, B64 8,082.625, B80 9,919.125.
All 256 sweep executions verified and the 192 measured executions replayed
identically (audit receipt SHA-256
`1a24f267ecec04dd337c2c64f603273f8a48e4058ab4c2a67d84d991c2e04dd1`).
The B48/B40 sweep interval includes 1, so B40 is also carried into the
confirmation as a prespecified near competitor. This file and the holdout
spec are written before any holdout target is measured.

The holdout fixes the registered 20-bit prime curve, eight new public
targets under target seed 20261017, one warmup and three measured rounds,
seed 20261017, B40/B48/B64 target-pinned batched pair-table arms, matched
`rho.negation` and its A/A repeat. All three IC arms have the same guarded
walk and max_trials 5000000. Preserve failures and audit all measured runs.
Primary metric: paired cold counted GAE, including setup, collection and
verification. Secondary: target trials, hit rate, phase and native counts,
memory, and L0 diagnostic timing. Selection is corroborated only if B48
retains a lower paired counted ratio than B64 on the new targets; B40 is
reported independently. The lower-bound GAE ratio does not establish a
full-cost or wall-time advantage.
