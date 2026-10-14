# Running campaign and frozen initial snapshot

**Superseded live location:** the scratch run below is no longer active. Its later
outputs and root-gate receipt were unavailable after scratch cleanup. See
[the recovery record](RECOVERY_20261011.md) for the active durable restart and
acceptance accounting. The initial snapshot below remains a historical record.

The requested **at least 2000 additional maps remains open**. Two real native
walks are running, each requesting 1050 independently replayed steps. The initial
snapshot contains 6 P192 and 19 P224 sealed steps; all 25 reproduce earlier maps,
so the snapshot's additional count is **0 / 2000**. It is a dated execution
snapshot, not a completed search. Use the live status command for later counts.

The persistent launch's LAUNCH.json was written at
`2026-10-10T22:31:46-0700`. Supervisor PID at launch: `26190`; its direct children
include both native walk processes and macOS sleep prevention. Output root:
`/private/tmp/isogeny-2000-more-20261010/run`. Persistent execution session: `78353`.

```sh
/private/tmp/isogeny-2000-more-20261010/installed/isogeny-campaign status \
  /Volumes/SSD990/crypto/worktrees/codex-p192-p224-large-degree-isogenies-20261008/research/p192_p224_hundred_large_prime_isogenies_20261010 \
  /private/tmp/isogeny-2000-more-20261010/run
```

The supervisor writes each child's exit receipt and automatically produces
`FINAL_SUMMARY.json` when both walks end. A successful final map count does not
itself archive the output, update canonical registry models, or merge the PR.
Those delivery requirements remain open until the results are synchronized.

The runtime CLI's clean source is `ecc2b2250c86037b8ab7ba656ac36b17ca31a5d7`.
The running supervisor's frozen source is
`316150cbabd74a53603f233afd44fd7f9f40250d`. [LAUNCH.json](evidence/LAUNCH.json)
binds both binaries, the protocol and the preceding campaign manifest. The later
artifact-inventory command is an archival helper, outside measured map execution.

The [initial report](snapshots/initial/RESULTS.md),
[exact records and phase sums](snapshots/initial/SUMMARY.json),
[visual](snapshots/initial/SEARCH.svg) and
[reviewed PDF](snapshots/initial/REPORT.pdf) are frozen together. They show the
excluded prefix work rather than treating it as additional yield. The PDF's two
pages were rendered and visually checked. The initial snapshot finds no new
graphable model: the canonical catalogue and route graph stay at the earlier
449-model, hundred-map campaign. Full earlier route identities and raw maps
remain in the preceding committed study. Canonical registration is still
required when archiving new findings from this run.

The first shell-background launch returned before its child initialized and
left no launch manifest, run directory or map output. The corrected launch uses
`nohup` in a persistent execution session; process and manifest checks confirmed
both workers remained active across later tool calls. No maps or timings are
credited to the failed launch. A first mutation-control invocation also returned
127 because an old executable path was absent; its diagnostic is retained.
The saved native control source was then rebuilt and the control run passed.

The first P224 constructor receipt reports 6,914,048 sampled peak resident bytes,
zero unavailable memory samples, and a successful child exit. Constructor RSS
and timeout checks remain enabled. Replay is a separate measured phase. The
shared host and other active work prevent a matched throughput claim.

The research SSD had 1.8 GiB free before setup and later fell to 1.0 GiB during
other host activity. The system disk had 18 GiB free. Raw walks and build/log
scratch therefore use the system disk. A temporary absolute path is live storage,
not final durable delivery: archive raw new maps before claiming the request
fully delivered.

During execution the research SSD fell to 32 MiB free. The inactive CLI cache
and the completed 3.2 GiB all-targets cache belonging to this worktree were moved
to the system disk, preserving their files. This recovered 3.6 GiB on the SSD;
11 GiB remained on the system disk at that check. The active root-gate executable
and both copied walk executables were left at their running paths. The cache
moves do not alter the frozen producer binary or its source identity.

PR 1606 remains ready for review. Its preceding head passed all four dedicated
isogeny release build profiles. Other hosted checks failed, including copied
historical replay lockfiles rejected by `--locked`, a canonical mirror HTTP 404,
and a formatting diagnostic in `prime_fast.rs`. The new head must be checked
after pushing. The released binary has not been published through ReleaseMe.
