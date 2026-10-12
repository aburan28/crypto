# Four equal-size n=41 base windows: no measured selection opportunity

The one-shot [hosted run 36827211814](https://github.com/aburan28/crypto/actions/runs/36827211814)
used merged `main` `4978478120ae83a322c221e59d65a358d5510f29`, the
[preregistered protocol](../base_window_screen_20261001/PROTOCOL.md), the
[frozen source lock](https://github.com/aburan28/crypto/pull/1129), and the
[merged held-out Q](https://github.com/aburan28/crypto/pull/1131) on
`icv1-f2m41-tm2308219-7f48b14a`. All 45
scheduled arms completed. Hosted and separate macOS independent replay agree
byte for byte: 40 full-rank compact traces, 40,960 compact target logarithms,
and 5,120 same-Q rho logarithms pass their group-law, four-point witness,
input, source and command checks. The quiet/PSI preflight, 79 one-second
contention samples, and all four A/A gates pass. **The frozen decision is
`NO_WINDOW_OPPORTUNITY` for these four deterministic 255-column bases on five
n=41 blocks.** It is a scoped negative batch diagnostic, not an ECC2K-130
no-go or a one-target speedup claim.

| Same-Q arm, K=255 and L=1,024 | Median charged CPU seconds | Paired CPU / rho, 95% log-t interval | A/B median, 95% interval |
| :-- | --: | :-- | :-- |
| Matched signed-Frobenius rho, 32 walks | 1.1416 | 1.000 reference | n/a |
| Window 0, accepted orbits 0&ndash;254 | 1.5707 | 1.388 [1.338, 1.423] | 0.988 [0.957, 1.032] |
| Window 1, accepted orbits 255&ndash;509 | 1.5916 | 1.367 [1.339, 1.456] | 1.006 [0.991, 1.050] |
| Window 2, accepted orbits 510&ndash;764 | 1.5962 | 1.398 [1.363, 1.435] | 1.022 [0.988, 1.063] |
| Window 3, accepted orbits 765&ndash;1019 | 1.6116 | 1.420 [1.363, 1.488] | 1.005 [0.946, 1.069] |

The unit is Linux `wait4` child user-plus-system CPU seconds on one AMD EPYC
7763 host. Each compact arm pays for a **fresh generator scan plus a fresh
complete compact process**: materialization, one 255-column S3 index,
full-rank relation collection, linear solve, all 1,024 Q descents, scalar
checks and output. The frozen executable uses W64 blocked prefilter and the
same n=41 curve representation as the strong normal-basis rho child. Each
window ran A and B on the same block of Q; its charged CPU is their geometric
mean, paired to that block's rho CPU. Median CPU seconds and paired median
ratios need not divide exactly. Process-launch wall time was recorded
separately and supplies no speed gain. The source used `rustc 1.98.1` and the
pinned v2 compact SHA-256
`702a0a05709bc14bc10bafbb08edb2a2f5f86794e970d84c317da3bc65cdcf38`;
it does not silently substitute today's different compact example.

The optimistic **post hoc** per-block minimum of the four charged window
ratios was `1.346, 1.385, 1.414, 1.343, 1.366`. Its median is **1.366**, with
a descriptive 95% log-t interval **[1.334, 1.408]**. Every minimum exceeds
one, and the interval's lower endpoint exceeds one, exactly satisfying the
precommitted stop rule. This lower envelope pays one full scan and index for
the winning arm but grants a hypothetical selector *zero* scoring overhead.
For these four same-source implementations on these five Q blocks, a
nonnegative-cost selector that chooses one complete window policy per block
cannot turn the observed CPU totals into a rho crossover. It says nothing
about other bases, an amortized shared scan,
different PDP solvers, degree-263 descendant-native bases, or n=131.

Generator median CPU rose only from 0.0292 s at window 0 to 0.0348 s at
window 3; the corresponding compact medians were 1.5415&ndash;1.5889 s. The
lowest fixed-window paired median still needs about **26.8%** lower complete
CPU merely to reach its same-block rho median ratio of one. Moving among
these four raw-x orbit windows did not produce that improvement. The earlier
[orbit-disjoint v2 result](../disjoint_cold_v2_outcome_20261001/RESULT.md)
remains a separate six-cell fixed-policy observation; its 1.314 n41 ratio
used different Q and a different hosted CPU and is not an A/B arm of this run.

The [sealed manifest](evidence_run_36827211814/MANIFEST.json) lists SHA-256
for all **419** raw members; its deterministic 34,192,237-byte archive has
SHA-256 `8ff5e0b5aed4c3f2d785ceee97db8666267e864af752bfe583cd70e1d0816a76`.
The manifest itself has SHA-256
`7af1742cde1081d311b0dc947a344e9daea7f321f127216f91e319a6248602ac`.
The source/input freeze, build logs and binary hashes, every child receipt,
failure field, scan trace, rank and target file, host record, and isolation
samples are retained. The hosted and macOS replay receipts are byte-identical
SHA-256 `4d2472f926d2847ed82bfdaa5e3182c5d1a349a4b98211dcfb43df880f0b5c5b`.
All four A/A medians lie within the frozen [0.90, 1.10] band and each interval
contains one. The isolation record reports zero contended samples; its largest
sampled other-process CPU use was 0.10 seconds, at the fixed limit. The
two-second preflight saw 0.02 other CPU seconds, 2.3% CPU PSI some and 0.0%
memory PSI some. Eighty-two user threads still had affinity allowing the
reserved CPUs, as on the earlier hosted panel, although the monitor observed
none exceeding the fixed contention gate. It cannot exclude shorter VM noise,
so the A/A test bounds only observed drift under this host and protocol.

The next evidence-ranked route is **natural high-arity descendant-native PDP
yield at equal useful base size**, followed by a non-eager index/query policy
whose setup and every target-dependent phase are charged. Changed descendant
endomorphism orders and discriminants are mathematical inputs to that study,
not evidence that its natural PDP yield is higher. The current result has no
common group-operation-equivalent `S`, n=53/83/131 transfer, Certicom
logarithm, or method-wide speed claim. The repository's primary one-target
IC-versus-rho acceptance gate remains separate from this 1,024-target batch
screen.
