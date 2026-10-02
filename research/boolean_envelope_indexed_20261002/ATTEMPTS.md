# Native qualification attempt journal

GitHub run 36949246631, attempt 1, source
`0a5b4e8f577d2f7f5badd3f03192f289713553a1`, passed the native Rust
correctness job and built the release worker and retained baseline on Linux
ARM64. Its campaign worker did not start. Immediately after compilation and
tests, the repository isolation controller measured CPU PSI `some avg10` 6.18,
above the frozen 5.0 limit. Other-process CPU was only 0.02 seconds over the
two-second settle interval, below that separate 0.20-second bound.

The original 29-member artifact is preserved in `failed_discovery_01` with
manifest SHA-256 `1761609502cb824531fb929ccdee4f01d1f1cb9ed0ebe6913eba635dcad93f80`.
It contains its binaries, source, build/test receipts, refusal stderr, empty
worker stdout, and explicit failed status. No observation is admitted or
pooled. No holdout was touched.

The next revision adds a bounded native readiness wait before the unchanged
isolation tool's independent locked preflight. Every readiness observation is
retained. No solve runs during this wait, and its duration is outside every
per-arm clock. The 10% CPU limit, 5.0 PSI limit, seeds, construction algorithms,
matrix caps and primary 2x gate are unchanged. The revision requires a fresh
complete discovery; this failed attempt remains failed.

GitHub run 36952544952, attempt 1, source
`8aa793e6673fd00e7e9b1e73324b7ef0133acd9c`, completed the new native
discovery. One two-second readiness observation passed, followed by the
isolation controller's unchanged preflight. The whole campaign used CPU 3 on
Linux ARM64: 12.6452819 wall seconds, 0.22 other-process CPU seconds, and
`contended=false`. The 36-member artifact is preserved unchanged in
`qualified_discovery_01`, manifest SHA-256
`c27dd012263c9627082431e2ea5d16725e56e49d558f80d1d29c3d5b8fa0eb38`.
Its result reports 128 complete cells, 8,960 A/B observations, 2,560 A/A
observations and 244,800 oracle-verified output matrices. Native replay and
all manifest hashes passed. The indexed arm cleared the incremental gate in
16/32 discovery groups; the degree-only arm cleared 0/32. Neither cleared
a 2x group. This is an admitted discovery-stage diagnostic, not a full or
cryptanalytic result. The fresh holdouts remain unused.

The comparator audit in `BOUNDARY_CORRECTION.md` found an already merged
packed-direct constructor on the same finite matrix families that the fixed
five-arm roster omitted. Its old-host times are not pooled with this Linux
ARM64 run. Because the new arms passed no dramatic discovery group even
against that weaker roster, no full run is dispatched for this candidate.
The unexecuted holdout seeds remain unused; the primary full gate and any
matched ratio to packed direct stay null.
