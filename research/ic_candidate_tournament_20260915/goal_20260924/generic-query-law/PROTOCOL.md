# Independent generic query-law replay

Follow PR 803 with a bounded, independent check of every exported query
coefficient against the frozen job, including unsuccessful and identity queries.
The current group checker proves witnesses and query retention but cannot detect
a different sampler, seed, or position when a replacement witness is also valid.

Implement the pinned rand 0.8.8 / rand_core 0.6.4 / rand_chacha 0.3.1 integer
stream in standard-library Python. Cross-check against the upstream StdRng vector
and a Rust executable using the real pinned crates: seeds 0, 1, 2026092556 and
u64::MAX, 257 draws per stream, subgroup bounds 2, 31, 127, 65587, 1439393 and
u64::MAX. Include trial indices across 64-query walks and wrapping-u64 keys.
Do not assume modulo reduction is rand's bounded integer sampling rule.

Replay swept and windowed collection and sampled and 64-lane walked descent.
Derive dispatch from the frozen worker job and actual base length, not from a
report's claimed sampler name. Limit this adapter to the existing one-target,
odd-degree 5..31 generic worker; other producers require their own adapters.
Bind the receipt to the job, query ledger and checker identity. Reject wrong
seeds, coefficients, chronology, batch partitions and resource limits, even
when the altered group relation remains true. Retain negative test results.

Replay the already archived PR 803 controls without changing them. Add bounded
fresh n9 and n13 controls, using fixture/algorithm seed 2026092556, including a
pair-table window, ignored windows, two final-LA choices and failed preparation.
Use one supplied public point and one Rayon thread. Preserve all worker output,
inputs, source/dependency/binary hashes and receipts. Extend Linux integration
so native and profiled reports must have identical independently replayed query
laws. Stop on a stream mismatch, correctness failure or lost attempt.

This is an accounting/correctness admission test. It changes no solver algorithm,
does not qualify implementation speed, does not replace source/base/matrix or
exclusive-phase admission, and consumes neither remaining improvement round.
All comparative online/cold/operation costs and speedups remain null.
