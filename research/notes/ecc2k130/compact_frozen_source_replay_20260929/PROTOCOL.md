# Compact-orbit frozen-source replay

The compact S3-batch and swap-quotient evaluations (`compact_s3_batch_20260929`,
`compact_swap_quotient_20260929`) pin, in their `FROZEN.json` `source_sha256`,
the bytes of shared files that later PRs legitimately edit:
`src/cryptanalysis/koblitz_index_calculus.rs`, `Cargo.toml`, `src/lib.rs`, the
workflow files, and others. Their `run_panel.py` and `verify_panel.py` hash
the checkout they run from, so the first PR to touch any of those files broke
the evaluation, even though nothing about the frozen evaluation had changed.

Following the historical-snapshot replay used for the rotated solver-admission
receipt (PR #801) and the n53 rank-rotation replay (PR #823), the workflows
now build and run from a materialized copy of the checkout instead of
re-pinning hashes:

1. `materialize.py --freeze <FROZEN.json> --dest <dir>` copies the sparse
   checkout (without `.git` and `target`) to `<dir>`.
2. For each pinned file whose live bytes differ from the freeze, it writes the
   committed snapshot `historical_snapshots/<sha256>.gz` instead. The snapshot
   is accepted only if `SNAPSHOTS.json` lists it under the same path, its gzip
   and decompressed SHA-256 and length match, and the decompressed hash equals
   the frozen digest. A drifted file without a matching snapshot fails closed.
3. Every pinned file in `<dir>` is re-hashed against the freeze.

The frozen `run_panel.py` and `verify_panel.py` then run unmodified from
`<dir>`, so their own source checks pass against the frozen bytes. `GIT_DIR`
points at the real checkout, so `host.json` still records the commit under
test. `FROZEN.json`, `SOURCE_FROZEN.json`, fixtures, and all recorded
measurements are unchanged.

`SNAPSHOTS.json` snapshots every pinned file outside the two panel directories,
taken from main at `source_commit`, where each matched its frozen digest.

## Amendment 2026-10-12: restore the unpinned source from the base commit

**What broke.** Every replay workflow that builds from the materialized tree
(`compact-s3-batch`, `compact-swap-quotient`, `compact-s3-prefilter`,
`compact-ir-ledger`, `compact-ir-cold-gap`, `ecc2k130-disjoint-cold-q`,
`ecc2k130-disjoint-cold-v2`) failed at `cargo build` on every recent pull
request, with 13 library errors (an unresolved `redis`, methods missing from
the pinned index-calculus crate, `include_str!` files the sparse checkout
does not hold). No measured arm ran, so none of these failures is an
evaluation result.

**Why.** A freeze pins only the files its own pull request touched: 20 to 25
paths, of which 4 are in `src/`. The other ~515 files of `src/`, the data files
`src/` embeds, and the Cargo lock came from the live checkout. Once `main`
gained a module that needs a newer pinned file, or a dependency the pinned
`Cargo.toml` lacks, the mix stopped compiling. A resync of the lock copy on
2026-10-09 (`932304d24`) made the mismatch visible: the lock the run used had
been replaced by one for a different `Cargo.toml`.

**What changed.** The frozen tree is now the base commit's tree with the pinned
files laid over it.

- `BASES.json` names, per freeze, the commit whose tree was current when the
  evaluation was built: `evaluation_source_commit` when that commit is in
  `main`'s history (the prefilter freeze: `2c207b41`), else `source_base_main`
  (`3a4e9859` for the S3 batch, `acc83e53` for the swap quotient). The S3
  batch's `evaluation_source_commit` (`5430011`) is not in history, so its
  unpinned files come from its base, and the 4 pinned `src/` files match that
  base byte for byte.
- `bases/<commit>.tar.gz` holds that commit's `src/`, the data files `src/`
  embeds, and `compact_shared_log_20260925/Cargo.lock`, as `git archive`
  wrote them. `BASES.json` records each archive's SHA-256 and length, the
  paths it holds and the paths absent at the commit, and the git object id of
  each path, so a reviewer re-derives it with
  `git rev-parse <commit>:<path>`.
- `materialize` wipes those paths from the copied checkout (so a live file
  the base never had cannot leak in), unpacks the archive, then applies the
  pinned snapshots exactly as before and re-hashes every pinned file.
  A freeze with no base, an archive whose hash or length differs, a member
  outside the base paths, or two freezes with different bases fail closed.
- The tool is the native crate `frozen-materialize/` (std, `sha2`, `flate2`,
  `tar`; its own `Cargo.lock`). It replaces `materialize.py` in the workflows,
  as the repository's no-Python rule for replay tooling requires.
  `materialize.py` and `test_materialize.py` stay unchanged for
  `paired_fullrank_20260925/verify_committed.py`, which imports them.
- `frozen-materialize snapshot-base` regenerates an archive from a full
  checkout. It is deterministic for a given commit and `flate2` version; the
  committed archive, not a regeneration, is what the workflows trust.

**What did not change.** Every `FROZEN.json`, `SOURCE_FROZEN.json`, fixture,
recorded measurement, snapshot and pin is untouched. The panel and verifier
scripts run unmodified from the materialized tree. The receipt schema is
`compact-frozen-source-materialization-v2` (v1 plus `pinned_matching_base`
and a `base` block).

**Local check, 2026-10-12.** From a sparse root shaped like the workflows'
checkouts, the three distinct freezes materialize and build with their
workflows' `cargo build --release --locked` commands; the S3 batch's example
tests pass in its tree. The hosted runs on the pull request that lands this
amendment are the check that the full panels still reproduce.
