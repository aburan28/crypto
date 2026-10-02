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
