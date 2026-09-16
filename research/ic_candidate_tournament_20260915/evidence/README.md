# Frozen tournament evidence

The archives preserve every raw profile, receipt, correctness certificate,
native process output, source snapshot, frozen evaluator and measured executable.
Readable contracts, fixtures, decisions, phase summaries and reports also live
under `../runs/`. Failed preparations, the batch development screen, source
adapter controls and intermediate candidate sources are retained in
`development-and-validation.tar.zst`.

Cargo build caches, Python bytecode and empty operation locks are omitted.
The measured executables remain inside the archives because their hashes are
part of the sealed experiment. Restoring evidence does not execute them.

## Restore and audit

Requirements: Python 3.11+ and the `zstd` command. The original streams were
produced with zlib 1.3. Run from the repository root:

```bash
python3 research/ic_candidate_tournament_20260915/evidence/restore.py

python3 research/ic_candidate_tournament_20260915/runs/round-0005-batch16/evaluator/tournament.py \
  verify --round research/ic_candidate_tournament_20260915/runs/round-0005-batch16
```

Use `--archive round-0006` to restore only the single-target parity round,
`--archive round-0005-batch16` for the 16-target round, or
`--out /absolute/path/to/evidence` to restore elsewhere. Verification uses the
frozen Python evaluator and checker, including all source/artifact hashes,
every recovered scalar, each phase cost, stage summaries and the final decision.
It needs neither Valgrind nor a Rust rebuild. Repeat for rounds `round-0002`,
`round-0003b`, `round-0004` and `round-0006` to audit their complete records.

The restoration program checks each archive SHA-256 from [manifest.json](manifest.json)
before extraction. It refuses unsafe paths and different existing files. Existing
identical files are retained, so restoration is repeatable. About 1.1 GB of
restored files are needed for the complete record, in addition to the archives.

## Packing a new round

`pack.py --round NAME` writes `NAME.tar.zst` in this format from a finished
round directory and records it in [manifest.json](manifest.json). It refuses a
profile whose gzip stream zlib cannot reproduce exactly, excludes build caches
and the operation lock, and normalises member order, modes and timestamps.
Round 0006 was packed with it; the earlier archives predate the script.

## Lossless profile packing

Individually compressed Callgrind profiles share extensive content, but a tar of
those gzip files exceeds GitHub's per-file limit. Each profile is expanded inside
the zstd tar so compression can share that content. Its original 10-byte gzip
header, 8-byte trailer and SHA-256 are retained in tar PAX fields. All other files
are stored unchanged.

`restore.py` reconstructs the original gzip bytes using raw level-9 DEFLATE and
checks their original hash. Every profile was checked for exact reconstruction
while packing. An incompatible zlib encoder fails explicitly; receipts and
their hashes are never rewritten. **Use `restore.py`, not plain tar extraction:**
expanded members retain their original filenames inside the archive.

## Source scope and paths

The measurements bind the archived source, lockfiles, compiler settings and
executables. The repository had unrelated uncommitted source at measurement time;
the complete frozen snapshots preserve that dependency context. They are evidence,
not changes to the library's current defaults. The cumulative candidate change is
[WINNER.patch](../campaign_20260916/WINNER.patch), relative to
`runs/round-0002/source`. Do not apply it to an arbitrary newer library revision.

Historical JSON and logs retain original absolute `/home/ubuntu/crypto/...` paths
as provenance. After relocating, use the same suffix below this research directory.
The 16-target winner source is
`runs/round-0005-batch16/source_candidates/combined_descent/source`; its configuration
is in [winner-config.json](../runs/round-0005-batch16/winner-config.json). The
single-target round-0006 candidates derive from that source through the
committed `campaign_20260916/round6-*.patch` files, which
`campaign_20260916/round6_candidates.py` reapplies onto the restored snapshot;
the frozen copies the measurements bind are inside `round-0006.tar.zst`.
Generate a new registry with `tournament.py propose --from-round ...` to obtain
local paths before preparing a new experiment with a new seed.

The top-level worker in this PR emits field metadata explicitly so it compiles
against current `main` without unrelated serialization derives. Frozen worker
copies and all measured results remain unchanged; a newly built worker requires
its own matched measurements before any performance claim.
