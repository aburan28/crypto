# Adaptive F6-IC pair-sum index pilot: negative

The preregistered pilot in [PROTOCOL.md](PROTOCOL.md) completed on October 6,
2026. This is an **exploratory, unisolated Mac measurement**, not a promoted
CPU speedup. All 12 fresh-process runs finished with exit code zero, a
verified recovered scalar, and five exclusive online phases summing exactly
to the reported online wall interval. [measurements.jsonl](measurements.jsonl)
contains one row per candidate, workload, and run ID; `runs/` holds all raw
stdout, stderr, status, timestamps, and SHA-256 hashes.

The frozen n17 curve has a 62-point usable base and 29 folded columns. T1's
public point is `(73407,129763)` with recovered scalar 4785 and workload
`146a1e9ee3c8`; T7's point is `(98625,98119)` with recovered scalar 2391
and workload `ced1677f0976`. The inherited F4 control, original F6-IC,
and pair-index F6-IC used the same worker binary
`96d13172e41242b6efd67bef5206bcc331183f7bb7f64791b94b279b52924001`
and the same prepared mathematical state. The exact candidate IDs and
source/input hashes are in [FREEZE.tsv](FREEZE.tsv) and `candidates/`.

| Target | Repetition | Inherited F4 online ms | Original F6-IC online ms | Pair-index F6-IC online ms | F4 / pair | F6 / pair | Pair geometric additions | Pair index builds |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| T1 | 1 | 3.623 | 3.049 | 2.982 | 1.215× | 1.023× | 633 | 0 |
| T1 | 2 | 3.629 | 2.919 | 2.867 | 1.266× | 1.018× | 633 | 0 |
| T7 | 1 | 242.416 | 152.824 | 154.185 | 1.572× | 0.991× | 62,295 | 10 |
| T7 | 2 | 211.201 | 155.107 | 164.691 | 1.282× | 0.942× | 62,295 | 10 |

The unchanged F6-IC A/A online ratio across the two repetitions is 1.045×
on T1 and 1.015× on T7; inherited F4 varied 1.148× on T7. These two
repetitions do not establish a statistical performance bound, and the host
failed the isolation gate. Peak RSS was not captured and is `null` in the
measurement rows. No rho reference was measured in this pilot, so the
IC-versus-rho online speedup and `S = total_operations / sqrt(r)` are unknown.

Original F6-IC counted 633 geometric additions on T1 and 79,635 on T7.
The new arm counted 633 and 62,295, respectively: 0% reduction on T1 and
21.77% on T7. T7 triggered ten pair-index builds, one on each failed PDP
attempt; each index contains `62 × 63 / 2 = 1,953` unordered pairs and is
discarded when its attempt ends. The last, successful attempt used the
legacy path. The pair arm had 300 indexed lookups across T7. All arms used
one T1 attempt and 11 T7 attempts; the F6 arms performed 9 and 691 Boolean
reductions respectively. The pair index is exact, including repeated
summands and both compatible orderings, as checked by the focused native
test and independent scalar replay in every run.

**Decision:** The pilot fails both preregistered T7 thresholds: the
geometric-addition cut is below 25%, and neither exploratory complete
F4/pair ratio reaches 2×. The T1 ≤5% addition-increase condition passes.
Do not expand this candidate to the eight-target panel or claim a 2× F6
speedup. The observed ten duplicate builds motivate a separately
preregistered shared-index variant; its construction cost and memory must
remain visible, and it needs a new exact candidate ID and paired runs.

Reproduce the derivation with `sh derive.sh`; its consistency check is
[DERIVATION_CHECK.json](DERIVATION_CHECK.json). The first runner invocation
stopped before timing when sandbox `sysctl` access failed; its receipt is in
`preflight-failure-1/`. The runner was committed with a fallback before the
12 measured executions.
