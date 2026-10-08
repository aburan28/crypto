# Full-point pair query: correctness passed, performance unconfirmed

The protocol was committed at `12a746d5a0af1b107a299431a5a5ed9d9e37a889`
and opened as [PR #1068](https://github.com/aburan28/crypto/pull/1068)
before the candidate was implemented or measured. This is its **public-Q
development** result, not a held-out or matched-rho timing result. Both
backends use the same W64 swap/Frobenius index, blocked root prefilter,
constructed base, rank seed 7, and point-only Q files previously published
in the compact-orbit point panel. The extra indexed full-point construction
is inside the candidate's index timer. No fixture scalar enters either
producer.

Independent Python group-law replay passed both arms on all three cells:
694 full-rank relations and 2,052 target logarithms in total, with every
target re-added from its four base points and `[d]G=Q` checked against the
verifier-only fixture. Each arm reached rank K in exactly K attempts. Both
arms have the same base hash, state/root counts, index S3 counts, rank logs,
query-candidate counts and recovered scalars. The candidate materialized a
rational full point for every indexed state and used **zero** query fallbacks
in these runs. The exhaustive GF(2^5) unit test separately compares the
group-law root pair with S3 for every ordered pair of rational points,
including zero/equal-x fallbacks. n37's different root order selected four
different target witnesses and changed 502 probe counts, but all selected
four-point sums and scalars replayed. No rank row changed.

| Public-Q cell | K | Query calls, rank + targets | Added pair-group adds | Extra y payload | Control internal total ms | Point internal total ms | Point/control diagnostic |
|:--|--:|--:|--:|--:|--:|--:|--:|
| n37 / 1,024 | 42 | 3,136 + 70,272 | 65,310 | 0.52 MB | 418.96 | 645.44 | 1.541 |
| n41 / 1 | 85 | 1,067,520 + 5,312 | 296,310 | 2.37 MB | 195.20 | 188.80 | 0.967 |
| n53 / 1 | 220 | 10,423,168 + 10,368 | 2,565,420 | 20.52 MB | 2,225.00 | 2,108.16 | 0.947 |

These numbers are one successive local macOS ARM64 release process per arm.
The producer's internal total timer excludes final summary and base-dump
serialization; no wait4 CPU, A/A drift, isolation monitor, peak RSS, matched
rho, or confidence interval was collected. They **cannot** support an attack
speed or even a stable engineering-speed claim. In particular, the earlier
unoptimized n37 development pair went the opposite way (817.89 ms control,
742.64 ms point), so the n37 release regression must remain visible and not
be selected away. The n41/n53 local differences are small and cannot by
themselves close the retained 24.7%/32.1% batch break-even distances or
the n41/n53 single-target gaps. The candidate nonetheless passes the
preregistered correctness gate and shows a plausible positive local
direction at two rank-heavy cells; the next decision needs a fresh-Q,
eligible cold CPU panel.

The [replay receipt](DEVELOPMENT_RECEIPT.json) gives every input and raw
output hash and each arm's independent rank receipt. The `development/`
directory retains control and point raw base, rank, target and summary
files for all three release cells, plus the earlier n37 debug outputs.
`python3 research/notes/ecc2k130/point_sum_query_20260930/check_development.py`
replays the two implementations independently from a full checkout and
checks the committed receipt and file manifest. The frozen source, binary,
host and command details are in `DEVELOPMENT_MANIFEST.json`. No candidate
has been promoted as faster; the canonical scoreboard remains unchanged
until an eligible full-process comparison exists.
