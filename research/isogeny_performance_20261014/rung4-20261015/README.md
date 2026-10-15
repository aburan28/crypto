# Fourth BMSS candidate: exact-output validation

Source `9a55ce6d2` (the complete SHA is in SOURCE_REVISION.txt), native binary
SHA256 `1e375ac4a49d317ca8375b6a7cc22c9c4b99b06102e41afea1847c00c454c8ab`.
The complete standalone source archive is frozen beside the executable. Its
crate files compare byte for byte with the original checkout used for the
passing standalone release suite and strict all-target Clippy.

The change separates direct small-precision Weierstrass recurrence from the
batch-inverted recurrence, removing per-iteration Option dispatch and unused
denominator conversion. This is a timing hypothesis prompted by the retained
AMD small-case regressions. No timing gain is established for this candidate.
The established first-three-rung CPU measurements remain a separate series.

A complete 96-case run covers all eight BMSS methods, four field widths and
prime degrees 7/31/101. The strict native reader compares it with the frozen
original baseline and the previous rung. Both exact full-coefficient output
comparisons pass. All four operation counters match the previous rung in all
96 cells. Setup, sigma and verification remain excluded from constructor counts.
The shared Mac supplied no qualified timing samples.

The first comparison command accidentally selected the retained header-only
`counts-complete.jsonl` from an old aborted launch. The strict reader rejected
it with `incomplete run` before creating a comparison. The accepted cumulative
comparison uses the complete 98-line `counts-valid.jsonl`; neither older file
was changed. The new run and both accepted comparisons are complete.

The next CPU wave independently measures original/candidate and
previous/candidate, five A/A and five interleaved A/B pairs in each series.
Settling lasts five seconds with the same CPU-rate and PSI admission limits.
All censored attempts stay visible; no partial series is promoted.
