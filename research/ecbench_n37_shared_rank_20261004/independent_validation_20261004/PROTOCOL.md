# Independent replay gate for the n37 shared-rank session

This protocol was committed before the independent replay. It validates the
already frozen, single-public-target session
`ECBS1h33c3e5e4f75d` from PR #1313. It does not measure a new target or
change the candidate, factor base, rank policy, or strong-rho reference.

Hypothesis: a release `ecbench` binary built from this PR on the GitHub-hosted
Ubuntu 24.04 x86-64 runner will reproduce every one of the session's 15
measured records bit for bit. The source session is macOS arm64, host class
`ECBENV2hd82681268e96`; the auditor must report a different non-null host
class. Replay tests the answers, factor-base identity, operation counts,
unpriced-work counters, phase accounting, and derived figures. The original
session's five stage times remain the only recorded online intervals.

The exact frozen inputs are
`research/ecbench_n37_shared_rank_20261004/sessions/mac_arm64_l0_02/`.
Their SHA-256 hashes are:

| File | SHA-256 |
| --- | --- |
| `session.json` | `6887bdf844ca09d74fc26c5f1145d5ebc88b368cdc26749460798012a9e4b006` |
| `spec.json` | `62108306f34b7362c98a3a9c7661f4ec98f8e0db82dca17540bd889ef9497a27` |
| `host.json` | `f1b0404543df2f17326a2acfaaf7538a7261f2369be510a23ecab793dc71abf0` |
| `plan.json` | `b432feebdf3e4f44390054fdcfa733b2fead29f2d5305077394d56f0f2e8928a` |
| `records.jsonl` | `8046061f5078f254627a343998b810c7e570fc9fe575f403466114e70a0b1693` |

The session's own source commit is `2da691940907c56ccd705fbf450336e71fbb0ee4`;
its source and binary hashes are retained in the session records. The
independent runner builds the current PR head and executes
`ecbench verify --replay-all --exit-code --out RECEIPT` against these
unchanged files. The CI job pins the session ID, original host class,
records hash, 15 reproduced runs, and passing audit. It uploads the full
receipt even if later unrelated CI steps fail. The resulting receipt and
GitHub Actions run URL will be added to this PR after execution, preserving
the exact auditor binary and host class.

Success requires all 15 replays to reproduce, every audit check to pass,
and the auditor host class to differ. Any failure, timeout, missing receipt,
or same-class replay stays a failed observation; it does not count as an
independent certificate. One complete successful replay is the stop rule.
No CPU wall-time ratio can be promoted from this hosted VM: PR #1313's
Mac timing remains L0, and this replay measures deterministic counts, not
isolated performance. An isolated, same-point Linux session and native-work
pricing remain separate follow-on gates.
