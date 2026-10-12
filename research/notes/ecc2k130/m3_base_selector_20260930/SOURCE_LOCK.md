# Source lock for the charged degree-7 m3 base selector

This lock follows the [preregistered protocol](PROTOCOL.md). Its PR must be
merged before the first candidate, challenge Q or held-out target is generated.
The lock pins the exact producer, independent verifier, score, tests, workflow,
prior-panel inputs and transitive arithmetic source files in `FROZEN.json`.
The producer and verifier each check every pinned byte before doing work.

The selected arm constructs and scores two quota-matched bases per seed and
role before Q exists. Candidate zero is the matched control. The producer
records all 64 policy cells (four seeds, two holdouts, two arms, four mapped
policies), their 512 case trails, per-candidate phase charges, cold field
operations, modular row operations and full-call timings. The independent
verifier regenerates candidates and held-out streams, brute-enumerates all
triples, eliminates every rank trajectory, checks scalar recovery using a
separate bit-polynomial point law, and reconstructs cold-cost ledgers.

The source-lock CI compiles the implementation, checks pinned hashes and runs
an archived-base score regression. It does not execute `produce.py`, so it
cannot generate new fixtures. Once this PR is merged, run one producer outcome
under the protocol caps, retain failures and partial receipts, pin the raw
`result.json` SHA-256 in an evidence manifest, then run `verify.py` with that
manifest. Archive all outputs and the narrow decision in a separate PR.

`cell_timings` include the audit tail; cold comparisons use
`cold_cost_to_rank_or_512` and `cold_mod_r_ops_to_rank_or_512`. No CPU-time
crossover follows from full-call timings. The score's rank is computed from
the first witness per supported source target, with no challenge coefficient.
