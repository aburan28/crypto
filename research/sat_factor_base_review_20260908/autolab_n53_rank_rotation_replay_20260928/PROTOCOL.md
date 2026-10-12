# Released n53 rotation freeze: replay across later shared-source changes

Status: replay repair. It does not rerun the n53 producer, change the
[recorded result](../autolab_n53_rank_rotation_20260925/RESULT.md), or admit
today's producer as the released one.

## Trigger

The [released protocol](../autolab_n53_rank_rotation_20260925/PROTOCOL.md)
froze fourteen source SHA-256s in its `FROZEN.json` and measured exactly once
(Actions run 36218672105) from release source head
`9759251e79797478386e0dad503da4286574a0b8`. Merging main into the PR branch
(`a5790fe2`) brought in main's opt-in `legacy_rank_fixture_lcg_v1`
factor-base selection (main commit `70800583`) and other shared-file edits.
The unmodified `check_protocol.py` then fails in the `hash-only` job:

| source key | path | frozen / release sha256 | live sha256 at `a5790fe2` |
|---|---|---|---|
| `rust_producer` | `examples/koblitz_s5_sat_instance.rs` | `99f0b311…9e61e071` | `1c7f4fbd…30d83378` |
| `cargo_toml` | `Cargo.toml` | `1435eec1…b0b6b6b8` | `eee9a6a6…4a610238` |

The producer change is real: besides the opt-in schedule, the default base
header now carries a `selection_mode` key. Re-pinning the frozen hashes to
today's bytes would silently substitute a different producer for the one that
produced the archived result, so the frozen hashes stay as they are. The
other twelve sources and all five inputs live in this research thread's own
directories and still match.

## Replay rule

`historical_snapshots/` holds a gzip copy of each drift-exposed source read
from `git show 9759251e:<path>`: `rust_producer`, `cargo_toml`, and the PR
`workflow` (whose `hash-only` job this repair edits, so it is frozen too).
Each decompresses, under a 1 MiB cap, to the exact SHA-256 recorded for its key
in the original `FROZEN.json`; a mismatch is a failure.

`replay.py` imports the **unmodified** original `check_protocol.py` and
`verify_archive.py`, and for the duration of the check replaces only those
three entries of `check_protocol.sources()` with the decompressed snapshot
paths. Everything else, including the checker's own hash, the schedule
regeneration, the certified-base header, the git ancestry checks, the archive
rehash, and the independent audit replay, runs exactly as frozen against the
live tree. `replay.py` also pins the original `FROZEN.json` bytes, so its
status, release heads, and every other frozen hash cannot move under it.

The `hash-only` job runs `--mode freeze`, `--mode preflight`, the fail-closed
controls in `test_replay.py`, and `--mode archive` on the committed evidence
bundle. The replay reports the live digests of the pinned keys as diagnostics;
they are not new frozen hashes. The one-shot `outcome` job is unchanged; its
release preflight still reads live sources and the committed archive already
blocks a second attempt, so it fails closed.

## Scope

No producer build, n53 run, rank claim, cost comparison, or transfer claim is
added. Success means only that the released freeze and the committed first
attempt still verify from the bytes that were frozen.
