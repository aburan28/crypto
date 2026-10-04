# n37 shared-rank column sweep

The [frozen protocol](PROTOCOL.md) and [specification](SPEC.json) were pushed in commit
`5508af3a` before planning targets or measuring outcomes. The unmodified solver
implementation is commit `694ac8706d6b5d7ee4083b9fe195f9a25efeecd6`.
The subsequent [native decision analyzer](../../examples/n37_rank_columns_analyze.rs)
only reads sealed records, saved comparisons, candidate manifests and both
audit receipts. Its [machine-readable decision](RESULT.json) and [interpretation](RESULT.md)
are separate from the preregistration.

The measured `ecbench` binary has SHA-256
`5f44420ad64b772ec287d52281a838dc4b04f0af27a1e95b1a926a84c21c71e6`.
The [plan](PLAN.json) has SHA-256
`80dd313447299ed193e3baf0a40ccb47323702b8c214fcc3a9c6f0696258ec7d`.
The complete [session](sessions/mac_arm64_l0_01) is
`ECBS1h9a1d373bbe54`, with 384 verified executions and 320 measured
executions across eight public one-target workloads. The raw
`records.jsonl` SHA-256 is
`82c8b31640074b2f829d88e546bf487ee4b598edf4ef26f779984d625ed7de3c`.
All 320 measured records replayed identically in the [local audit](AUDIT.json)
and in the [independent Linux audit](independent_validation/RECEIPT.json).
The latter is from [CI run 37182953423](https://github.com/aburan28/crypto/actions/runs/37182953423)
on a different host class, and its receipt has SHA-256
`efb14a5cb1c60291accc4d28518178f6a7d7ed5305178459936252b346b94029`.

The saved [candidate claims](candidate_claims) hold the exact physical point
inventory and immutable `IC1` manifests for each base. The
[comparison files](sessions/mac_arm64_l0_01/comparisons) retain every paired
operation and exploratory wall result, including K42 A/A and K8 versus K16.
The [native table](TABLE.txt) is a quick index to those rows. `RESULT.json`
is reproducible with:

```sh
cargo build --release --locked --example n37_rank_columns_analyze
target/release/examples/n37_rank_columns_analyze \
  research/ecbench_n37_rank_columns_20261004/sessions/mac_arm64_l0_01 \
  research/ecbench_n37_rank_columns_20261004/candidate_claims \
  research/ecbench_n37_rank_columns_20261004/AUDIT.json \
  research/ecbench_n37_rank_columns_20261004/independent_validation/RECEIPT.json \
  /tmp/n37-rank-columns-result.json
cmp /tmp/n37-rank-columns-result.json research/ecbench_n37_rank_columns_20261004/RESULT.json
```

The selected K16 base is a **counted engineering lead over this K42 base**,
not an admissible IC/rho speedup. Mac CPU timing has L0 isolation, and both
arms leave native costs unpriced. The fully priced cold ratio and the primary
one-target online speedup remain unknown. No n131 inference follows.
