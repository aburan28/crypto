# Root-finding relation collection validation receipt

This receipt is for the source based on `51f977323e43a37d2719791354fd32ef0da6cb1d`
plus the five Rust source edits in this follow-up. The pre-archive binary
SHA-256 was
`1ff5b6d5f717f79cdfc557fde0017b0fac7e40502fe26a28bd7cdb2a99354904`
for `target/release/ecbench`; the release library test binary SHA-256 was
`ee8ce469bd00a619aedf48b73bd2a6f55ce18fa742b8293ecea53e9dc414b53a`.
The five-file `git diff --binary` SHA-256 was
`0d07c120d91f77bf609a26fe421590f27bf1db6cf0174cb7c6a308bd91049b54`.
The archive commit binds the final source paths separately.

## Scientific execution gates

| Session | Executions | Measured full replay | Exact paired check |
| --- | ---: | ---: | --- |
| [Target-pinned sweep](sessions/prime20-pinned-sweep-20261013/) | 256/256 verified | [192/192 identical](replay-all-sweep.json) | [24/24 full-rank first pins and 24/24 subtraction/table matches](sweep-pair-check.json) |
| [Unseen-target holdout](sessions/prime20-pinned-holdout-20261013/) | 160/160 verified | [120/120 identical](replay-all-holdout.json) | [rho A/A exact](holdout-pair-check.json) |
| [Native one-lane walk](sessions/prime20-native-walk-20261013/) | 128/128 verified | [96/96 identical](replay-all-native.json) | [24/24 native/generic relation-stream and answer matches](native-pair-check.json) |

The native walk's direct 4,096-step seeded `(point, a, b)`
[test](native-walk-correctness.txt) passed. The independent-lane
[test](native-batch-correctness.txt) passed at widths 1, 2, 4, 8 and 16.
The [pair-table occupancy](pair-table-occupancy.txt) and
[lane kernel profile](native-lane-profile.txt) tests passed; the latter's
short L0 timings are diagnostic only. The release `ecbench` build passed.
`git diff --check` on the hand-edited Rust and Markdown paths passed;
the full staged check reports trailing spaces from Matplotlib SVG path
data and terminal blank lines in preserved raw test logs. Those generated
artifacts remain unmodified.

## Repository release suite

`cargo test --offline --locked --release --lib` ended with **4,065 passed,
13 failed, 116 ignored**, exit 101, after 1,487.59 seconds. Preserve the
[full output](cargo-test-release-lib-parallel.log). Four failures were in
`cryptanalysis::jv_cover::tests`; nine were in
`cryptanalysis::pollard_collab` loopback tests. No new native walk or lane
batch test failed.

The four `jv_cover` cases passed in a focused serial
`cargo test --offline --locked --release --lib cryptanalysis::jv_cover::tests:: -- --test-threads=1`:
[14 passed, one ignored](cargo-test-release-lib-jv-serial.log), exit 0.
The loopback group passed with socket access using
`cargo test --offline --locked --release --lib cryptanalysis::pollard_collab:: -- --test-threads=1`:
[37 passed](cargo-test-release-lib-pollard-loopback.log), exit 0. These
focused runs classify the observed failures at this revision; they do not
turn the full suite into a clean pass.

A preceding serial full-suite attempt was interrupted after prolonged CPU
contention to run the repository's standard parallel command. Its
[partial output](cargo-test-release-lib-serial-interrupted.log) is retained
as an interrupted infrastructure attempt, not a test result.

All ecbench wall records in this follow-up are L0. The counted GAE fields
leave several native operations unpriced, and replay took place on the same
host. See the [open complete-cost problem](../../../problems/rootfinding-complete-cost-calibration.md).
