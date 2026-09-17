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

Use `--archive round-0008` to restore only the certified single-target round
(`--archive round-0009` for the same pipeline on the leaner shared job,
`--archive round-0010` for the lean static executable measured under the
evaluator that no longer forks itself, `--archive round-0011` for the
general arithmetic in words, `--archive round-0007` and
`--archive round-0006` for the strict-win and parity rounds it built on),
`--archive round-0005-batch16` or `--archive round-0006-batch16` for the
16-target rounds, or `--out /absolute/path/to/evidence` to restore elsewhere.
`round-0006-batch16` (incumbent retained; index-calculus admission enforced by
per-target descent certificates) verifies with its own frozen evaluator under
`runs/round-0006-batch16/evaluator/`. New rounds are packed with
`pack.py --name ROUND runs/ROUND`, the inverse of `restore.py`, which checks
every profile's gzip reconstruction before writing it and appends the archive
to `manifest.json`. Verification uses the
frozen Python evaluator and checker, including all source/artifact hashes,
every recovered scalar, each phase cost, stage summaries and the final decision.
It needs neither Valgrind nor a Rust rebuild. Repeat for rounds `round-0002`,
`round-0003b`, `round-0004`, `round-0006`, `round-0007`, `round-0008`, `round-0009`, `round-0010` and `round-0011` to audit their complete records.

The restoration program checks each archive SHA-256 from [manifest.json](manifest.json)
before extraction. It refuses unsafe paths and different existing files. Existing
identical files are retained, so restoration is repeatable. About 1.1 GB of
restored files are needed for the complete record, in addition to the archives.

## Packing a new round

`pack.py --name ROUND runs/ROUND` writes `ROUND.tar.zst` in this format and
appends it to [manifest.json](manifest.json). It refuses a profile whose gzip
stream zlib cannot reproduce exactly, excludes build caches and the operation
lock, and refuses to overwrite an existing archive or manifest entry. Rounds
0006-batch16 and 0006 through 0009 were packed with it; the earlier archives
predate the script.

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
the frozen copies the measurements bind are inside `round-0006.tar.zst`. Round-0007's
candidates derive from round-0006's winner source the same way
(`round7-*.patch`, `round7_candidates.py`, `round-0007.tar.zst`), and round-0008's
certified incumbent from round-0007's winner source (`round8-tiny2_cert.patch`,
`round8_candidates.py`, `round-0008.tar.zst`), round-0009's baseline from
round-0008's source (`round9-fastcurve.patch`, `round9_candidates.py`,
`round-0009.tar.zst`), and round-0010's baseline and ablation control from
round-0009's source (`round10-lean.patch`, `round10-lean_stdprobe.patch`,
`round10_candidates.py`, `round-0010.tar.zst`; the snapshots carry the
`.cargo/config.toml` that selects the static link), and round-0011's
baseline and IC-only challenger from round-0010's source
(`round11-wordfield.patch`, `round11-fastio.patch`, `round11_candidates.py`,
`round-0011.tar.zst`).
Generate a new registry with `tournament.py propose --from-round ...` to obtain
local paths before preparing a new experiment with a new seed.

The top-level worker in this PR emits field metadata explicitly so it compiles
against current `main` without unrelated serialization derives. Frozen worker
copies and all measured results remain unchanged; a newly built worker requires
its own matched measurements before any performance claim.
