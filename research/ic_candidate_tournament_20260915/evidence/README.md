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
general arithmetic in words, `--archive round-0012` for the musl static
executable with the arena allocator, `--archive round-0013` for the scalar
products in projective coordinates, `--archive round-0014` for the worker
built as one optimisation unit, `--archive round-0015` for the scan and the
resolution limit of the native gate, `--archive round-0016` for the widened
eight-cell panel on which `beats_rho_strict` fails, `--archive round-0017` for the
orbit-representative certificate that restores it on all eight, `--archive round-0007` and
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
`round-0003b`, `round-0004`, `round-0006`, `round-0007`, `round-0008`, `round-0009`, `round-0010`, `round-0011`, `round-0012`, `round-0013`, `round-0014`, `round-0015`, `round-0016`, `round-0017`, `round-0018`, `round-0018b` and `round-0019` to audit their complete records.

The restoration program checks each archive SHA-256 from [manifest.json](manifest.json)
before extraction. It refuses unsafe paths and different existing files. Existing
identical files are retained, so restoration is repeatable. About 2.0 GB of
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
`round-0011.tar.zst`), and round-0012's baseline, challenger and two ablation
controls from round-0011's source (`round12-musl.patch`, `round12-arena.patch`,
`round12-fastio.patch`, `round12_candidates.py`, `round-0012.tar.zst`; the
snapshots carry the `.cargo/config.toml` that selects the target and link),
and round-0013's baseline and two IC-only challengers from round-0012's source
(`round13-ld.patch`, `round13-canon.patch`, `round13-fastio.patch`,
`round13_candidates.py`, `round-0013.tar.zst`), and round-0014's baseline and
two IC-only challengers from round-0013's source (`round14-lto.patch`,
`round14-arena-tests.patch`, `round14-canon.patch`, `round14-fastio.patch`,
`round14_candidates.py`, `round-0014.tar.zst`), and round-0015's two IC-only
challengers from the unchanged round-0014 source (`round15-scan.patch`,
`round15-fastio.patch`, `round15_candidates.py`, `round-0015.tar.zst`). Round-0016
changed no source at all — both arms are the trees round-0015 sealed — so it ships
its analysis instead: `round16_resolution.py`, `round16_cell_census.py`,
`round16_floor_probe.py` with `round16-floor-probe.patch`, `round16_legacy_subset.py`
and `round-0016.tar.zst`, whose 2,142 receipts carry both the eight-cell verdict and,
through the legacy-subset script, the five-cell one computed from the same trials. Round-0017
ships `round17-orbits.patch`, `round17-rows.patch`, `round17-oracle-orbits.patch` (the checker's
second certificate format), `round17_candidates.py`, `round17_base_sweep.py`, `round17_measure.py`
and `round-0017.tar.zst`, whose 2,340 receipts include both certificate formats side by side on
the development and selection fixtures where `scan_io` and `orbits` both ran.

The first `round-0017.tar.zst` (sha256 `48a120ed19c5f2141cc1f5c7d879dd721b6705a82760b581a9fc56f888dc15ae`,
committed in b9b0fee) was packed before `report.py` had produced its outputs: the
report's support audit did not yet understand the orbit-named certificate and
failed, and the post-run chain masked the failure. That archive holds the same
2,340 receipts and the same decision but no `REPORT.md`, `measurements.json` or
`admission.json`. It is superseded, not rewritten: the archive the manifest now
names was packed from the identical receipts after the report ran, and the
earlier hash is recorded here so the two cannot be confused.
Generate a new registry with `tournament.py propose --from-round ...` to obtain
local paths before preparing a new experiment with a new seed.

The top-level worker in this PR emits field metadata explicitly so it compiles
against current `main` without unrelated serialization derives. Frozen worker
copies and all measured results remain unchanged; a newly built worker requires
its own matched measurements before any performance claim.

## Single-target continuation

Additional complete archives: `round-0006b-single`, `round-0007-single-policy`,
and `round-0008-single-implementation`. Restore each with `--archive NAME` and
run its own frozen evaluator as above. `single-target-development` preserves
the failed round-0006 preparation, all candidate snapshots, plans, patches,
controllers and equivalence-test logs. Build caches are excluded.

[Continuation results](../single_target_20260916/RESULTS.md) include retained
smoke rejections and distinguish promotion from parity and strict beating.
[Archive validation](../single_target_20260916/archive-validation.json) records
successful fresh-directory restoration and full audits of all three rounds.
[The single-target winner](../single_target_20260916/WINNER.json) identifies its
source/configuration; its cumulative patch is relative to round-0002's source.

## Round 0018 comes in two archives, and both are kept

`round-0018.tar.zst` (sha256 `11379eafac1827b652a9aca6e5beb31f7ea4af4533d311b2c21bf4e46eb1624c`) is a
complete, fully verified tournament that answers the wrong question: its incumbent is round 0017's
*incumbent* rather than round 0017's promoted winner, because `prepare` was given
`runs/round-0017/source` instead of `runs/round-0017/source_candidates/orbits/source`. Its incumbent
worker hashes to `7ca9953d…` against the winner's `e9f263b8…`, and its baseline source carries no
`factor_base_orbits`. Every arm-versus-incumbent number in it therefore conflates round 0018's two levers
with round 0017's certificate change, and its `beats_rho_strict=False` is not a statement about round
0017, whose winner was never in it.

`round-0018b.tar.zst` (sha256 `99efc1b7b0a0e65e58a068e960188197be7c1683517c934dfb36712446c0a36d`) is the
same pre-registration, the same seed and the same arms against the incumbent that pre-registration names,
and it is the round of record. The first archive is kept because its 2,340 receipts are real measurements
of those arms against an older baseline, and because a superseded run that is deleted cannot be checked.

## Round 0019 is the promoted record

`round-0019.tar.zst` (sha256 `11cf99af0cf018ab970d883e141c4224d1336b1bfa0da5795b786321a81f8054`,
68,485 files) holds the round that promoted `both` and carries the campaign's
strict-win record: `beats_rho_strict` on all eight cells in both metrics on
both final stages, seed 2026092119, audit VERIFIED over 2,646 trial receipts
and 452 source files.

It is the first round with a per-cell confirmation allocation — forty fixtures
at `n23a1`, twelve at the other seven, 124 cases against the flat panel's 96 —
so its `contract.json` carries `confirmation_cases_per_cell` and
`confirmation_allocation`. Rounds before it have neither key and are read
exactly as before; the flat profile is what `prepare` produces when
`--confirmation-cases` is absent.
