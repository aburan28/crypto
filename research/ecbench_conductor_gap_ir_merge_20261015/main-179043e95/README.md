# Latest-main publication validation — 2026-10-15

Validated head: `0af35c7491e5f0725fa07d6e426d6d7930577f41`, tree
`2dcf20eef50729e3bd5bef7fd0c442517b5add51`. This normally merges main
`179043e95a80084ebc32d083a74bd0c78d99fc45` into the prior validated merge
`efb72d1954b2c687bdf79fa83a34aaffdc265ecb`.

This is correctness and publication validation. No new experimental target,
solve or instruction profile was produced. The original frozen census and
all earlier successful and failed validation records remain unchanged.

| Check | Observed result | Evidence |
|---|---|---|
| Full release library, one test thread | 4,105 passed, 0 failed, 117 ignored | `library-tests.log` |
| Isolated strict-session test | 1 passed, 0 failed, 0 ignored | `strict-session.log` |
| Native driver example | 2 passed, 0 failed | `driver-tests.log` |
| Release library cargo check | exit 0 | `cargo-check.log` |
| Newly compiled native archive analysis | exit 0; both analysis files byte-identical | `native-replay.log`, `native-replay.status`, `native-replay.sha256` |
| Site | 48 passed | `site-tests.log` |
| Full and branch-diff curve-name lint | 0 problems | `curve-names.log` |
| Official generated leaderboard, browser and panel index | all freshness checks passed; 303 panels current | publication checks before commit |

The validation wrapper ran from 2026-10-15T02:18:37Z through
2026-10-15T03:20:41Z. It used release opt-level 3, default features, three
compilation jobs, 256 codegen units and one library test thread. These
builds supplied no new experimental instruction counts. The actual strict
test execution is recorded separately, including its one-test summary.

Validation began while the non-source merge conflicts were being reconciled.
`validation.revisions` records both parents. The full compiled-source object
inventory in `validation.source-index` exactly matches `source-index-final`
from the committed head; there are no unstaged source changes. `validate.sh`,
the `.command` files and `replay.sh` retain the commands. Native analysis used
a separate scratch copy of the archived inputs, leaving the original
analysis files untouched.

Both parents' timeline points are retained exactly. The benchmark index is
their exact union of 22 entries. Generated views retain main's 110 sessions,
198 methods, 373 factor bases and 13,402 yields. The canonical timeline,
embedded copy and table date this PR's four points to 2026-10-15.

All 257 files in the two original study trees retain their original Git
objects; `frozen-final.tree` and `original-final.tree` match the before-merge
inventories in the parent validation directory. Native replay again yields
analysis JSON SHA-256 `c211baf8a91f79d4d03773a7de041aa69c58af7e12ba12ccfe6fab8c4d12b871`
and Markdown SHA-256 `bf6ed84a1ecdb23b9c815322d57e34b4a55a47d4c3d1be672823a418f891e597`.

The earlier 4,088-pass/two-failure run and the later 4,102-pass run remain
historical evidence. This later result does not erase those executions or
establish intrinsic hardness. Current-head GitHub CI and reviews are still
required before the PR may merge.

The required Conductor scope check was attempted again and exited 127 because
the CLI is unavailable. No guard or test was disabled to publish this record.
