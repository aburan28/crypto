# Dedicated binary execution

The installed `isogeny` executable runs independently of Cargo and the checkout.
For this campaign, the native supervisor copies the clean executable into its
fresh output root and launches these two commands concurrently:

```sh
./isogeny walk --curve p192 --ell 2341 --count 1050 --timeout 1800 --out ./p192
./isogeny walk --curve p224 --ell 1117 --count 1050 --timeout 1800 --out ./p224
```

Each command must receive a fresh output directory. The first 50 steps per class
reproduce the prior campaign and are excluded from additional-map counts.
Each later step advances only after exact independent arithmetic replay.
The CLI seals `walk.json` after every verified step, even during a running walk.

The campaign utility has separate launch and read-only integrity status modes:

```sh
isogeny-campaign launch /absolute/path/isogeny /fresh/campaign/root /absolute/path/PROTOCOL.md /absolute/path/old-hundred-study
isogeny-campaign status /absolute/path/old-hundred-study /absolute/path/campaign/root
isogeny-campaign report /absolute/path/old-hundred-study /absolute/path/campaign/root /fresh/report/directory
isogeny-campaign manifest /absolute/path/frozen-study
```

`status` checks the frozen earlier kernels and all sealed current replay hashes,
prints actual additional counts, and does not construct maps or repeat the
algebraic replay. `report` produces a snapshot, diagram and PDF in a fresh
directory. Preserve snapshots rather than replacing an earlier observed state.

`launch` records the clean CLI version, binary/protocol hashes, supervisor PID,
stdout/stderr, and child exit receipts. Its process waits for both walks.
It writes `FINAL_SUMMARY.json` after both walks exit, checking the final unique
additional count even if no chat turn is active. The final integrity receipt is
separate from archiving, canonical registration and PR delivery.
Run it detached with `nohup` when starting a long campaign. It invokes macOS
`caffeinate` when available to keep the host awake for its lifetime.

An additional full-chain replay is available as a separate command:

```sh
./isogeny verify-walk --input ./p192
./isogeny verify-walk --input ./p224
```

Replay of a nonempty sealed prefix may pass while `complete` is false. It does
not discharge the 1050-step request or the additional-map target. Do not invoke
it on actively changing output when expecting a stable whole-chain receipt.

The source-controlled supervisor is a research accounting utility. The CLI
remains the dedicated released runtime tool; ordinary users run `isogeny`
directly. ReleaseMe builds remain in the existing PR's release workflow.
