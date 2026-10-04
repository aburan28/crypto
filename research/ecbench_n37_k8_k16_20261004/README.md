# Untouched K8/K16 confirmation

The [preregistered protocol](PROTOCOL.md) and [frozen ecbench spec](SPEC.json)
were pushed at `61f3946c` before the new target plan or any solve. The
[native plan](PLAN.json), generated afterward, contains 16 unique public
targets and 384 executions. Its 16 workload IDs and exact target points have
empty intersections with the preceding eight-target panel. The plan file has
SHA-256 `af6f80d8b3a9fdafdd3f78de8c062d7a4ab2268017e93fbdaac3abf3209afbf4`;
the spec file has SHA-256
`96dbd9b4b77ae876b057c5ae5e6bc121095aaa1a303a1d7513699e43c215fdca`.

Measurement, replay, Callgrind profiles, and the decision are pending. This
branch will retain all failures and land the result and canonical scoreboard
update in the same pull request.
