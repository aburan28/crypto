# cairn registration

The study's two reproducible results are registered as unfunded (`reward 0`)
**replay objectives** on a cairn ledger. Each claim's artifact holds exactly
what the pinned command prints. cairn re-ran the command and accepted the
claim only because every declared field came out identical.

The ledger is the session's local one (`~/.cairn/cairn.jsonl`, the path in
this repository's `.mcp.json`). It lives in an ephemeral cloud container and
is not synced to any network. `ledger-snapshot.jsonl` is that ledger,
byte for byte, after the last claim. Post the objective files to any node to
register them there.

| objective file | objective id | claim | verdict |
|---|---|---|---|
| `objective-pdp-split-replay.json` | `sha256:f4375085…` | `sha256:24dd83be…` | **invalid_spec**: cairn refuses `branches_timed_out` as a reproducible field because it is host-dependent. Kept as recorded. |
| `objective-pdp-split-replay-v2.json` | `sha256:9cedb98e…` | `sha256:c4cff5cb…` | **unavailable**: run with the node root `~/.cairn`, where `target/release/pdp_split_dreg` does not exist (`cwd` resolves against the objective root). |
| `objective-pdp-split-replay-v2.json` | `sha256:9cedb98e…` | `sha256:efa4bc2c…` | **accept**: replay reproduced 5 declared fields |
| `objective-rho-session-replay.json` | `sha256:b6a0a401…` | `sha256:91b16c3b…` | **accept**: replay reproduced 6 declared fields |

- **What the accepted PDP claim reproduces.** Cell `3:0:13:4` of the
  amended grid (`runs-v2/3-0-13-4.jsonl`): 24 runs; per run, the solving
  degree, F4 word XORs, verified count and exactness (`profile`); 5,680,776,626
  XORs in total; 6 verified decompositions; every run complete and exact.
- **What the accepted rho claim reproduces.** `ecbench verify --replay-all`
  of `sessions/rho-01` returns session `ECBS1he23809badf5b`, status
  `complete`, `ok`, 384 records, 384 verified, and the hashes of its five
  files.
- **Wall time is in no field**, by cairn's rule and this study's.

## How the claims were checked

No `bwrap` is installed in this container, so cairn's reference verifier
ran the commands unjailed. The accepted runs used a scratch node root that
held what the pinned commands name: the two release binaries built from
this branch and a copy of `sessions/rho-01`. The commands were:

```sh
CAIRN_EPOCH_SECONDS=5 cairn --log ~/.cairn/cairn.jsonl --root ROOT \
  try cairn/objective-pdp-split-replay-v2.json \
  --submitter pdp-split-dreg-session --artifact cairn/artifact-pdp-split-replay.json
CAIRN_EPOCH_SECONDS=5 cairn --log ~/.cairn/cairn.jsonl --root ROOT \
  try cairn/objective-rho-session-replay.json \
  --submitter pdp-split-dreg-session --artifact cairn/artifact-rho-session-replay.json
```

**These acceptances are this host replaying itself.** They are not an
independent replication. A node on another machine would provide one.
