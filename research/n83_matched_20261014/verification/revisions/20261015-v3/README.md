# Verifier-only correction for copied floating-point metadata

This revision preserves the frozen solver/controller sources, executables,
protocol, campaign files and failed v2 audit. It changes no measured workload,
timing, scalar, witness, resource envelope or experimental result.

The completed campaign's v2 audit failed because `equal_json` required literal
numeric token equality between original records and copies that the Rust
controller had parsed and serialized again. Native read-only diagnosis of all
32 processes found127 differing numeric tokens: all127 were finite f64 values,
all127 differed by exactly one decoded ULP, and none were equal decoded doubles.
No integer, string, key or type disagreement was found. See
`diagnose_numeric_bindings.cpp` and `diagnosis.json` for every exact path and
pair of original tokens.

For example, AA04-A's original summary `timing_ms.targets_total` is
`3723.0182520000003`, while its observation copy contains `3723.018252`.
AA01-B's copied process diagnostic `semantic.producer_phase_gap_ms` contains
`0.00019999999949504857`, while the original execution receipt contains
`0.00019999999949504854`. This is one-ULP parsing/reserialization drift, not
merely two spellings of an identical decoded double. Requiring decoded-double
equality alone would still reject all127 cases.

`verify_matched_campaign.cpp` is a new source copy of v2. Its original
`equal_json` and every arithmetic, subgroup, scalar, witness, matrix, resource,
hash and complete-panel gate remain unchanged. Only four observation-copy
bindings call `equal_reserialized_json`: embedded execution receipt, copied
process costs, copied producer summary, and copied producer target records.

The new helper permits exactly equal or adjacent finite decoded f64 values
only at explicit context-dependent diagnostic metadata paths. Both differing
tokens must satisfy complete JSON float grammar with a decimal point or
exponent. Different integers, integer-to-float substitutions, strings,
booleans, seeds, counters, nanosecond phases/residuals, coordinates, scalars,
witnesses, identities and hashes remain exact. A metadata field name under an
unrecognized parent does not receive this allowance. Protocol, configuration,
matrix and A/A comparisons still call the original exact comparator.

The pre-existing eight-epsilon checks for arithmetic-derived preparation
splits and nanosecond-to-millisecond conversion remain unchanged. They are
separate from the new, stricter one-ULP observation-copy rule. Two decoded
ULPs and a one-nanosecond timing change are rejected by the copy comparator.

The native self-test passes4 positive cases,30 negative cases and442 actual
authenticated-copy comparisons. Negatives include two ULPs, a one-nanosecond
change in seconds and milliseconds at the2100-second cap, adjacent integers
above2^53, one-unit u128 phase changes, scalar/hash/seed/boolean changes,
missing/extra keys, reordered witnesses, invalid ancestor contexts, changed
types, malformed numeric tokens and nonfinite values. See
`binding_self_test.cpp`, `self_test_receipt.json` and retained build/run logs.
The first self-test compile attempt failed on a test-only string-to-C-string
conversion; its log is preserved. No campaign run used that attempt.

Commands run from the repository root; no solver is invoked:

```sh
g++ -std=c++17 -O2 -Wall -Wextra -I research/n83_matched_20261014/verification research/n83_matched_20261014/verification/revisions/20261015-v3/diagnose_numeric_bindings.cpp -lgmpxx -lgmp -o /workspace/n83_matched/verification/diagnose_numeric_bindings_v3
/workspace/n83_matched/verification/diagnose_numeric_bindings_v3 /workspace/n83_matched/campaign research/n83_matched_20261014/verification/revisions/20261015-v3/diagnosis.json
g++ -std=c++17 -O2 -Wall -Wextra -I research/n83_matched_20261014/verification research/n83_matched_20261014/verification/revisions/20261015-v3/binding_self_test.cpp -lgmpxx -lgmp -o /workspace/n83_matched/verification/binding_self_test_v3
/workspace/n83_matched/verification/binding_self_test_v3 research/n83_matched_20261014/verification/revisions/20261015-v3/self_test_receipt.json /workspace/n83_matched/campaign
g++ -std=c++17 -O2 -Wall -Wextra -I research/n83_matched_20261014/verification research/n83_matched_20261014/verification/revisions/20261015-v3/verify_matched_campaign.cpp -lgmpxx -lgmp -o /workspace/n83_matched/verification/verify_matched_campaign_v3
```

The new checker retains v2's seven-positional-input CLI and create-new output
rule. Its full campaign replay must be run and metered separately after source
review; passing the binding self-tests does not establish scalar/witness replay
or successful complete-panel admission. Its receipt records the new checker
and numeric-helper hashes and the original v2 checker source hash. Any later
collector correction also needs a separate source/version; the original
controller and raw observation files must remain unchanged.
