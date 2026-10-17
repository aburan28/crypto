# K1280/root11 leads the current-source exposed-point N83 controls

K1280/root11 dynamic recovered the previously exposed N83 public point in
**4.73325 ms online**, with **212,480** distinct usable factor-base points,
**1,280** folded columns, and a **1,342,177,280-byte** packed root allocation.
Its cold receipt was **146.960 s**. Full independent arithmetic replay passed
for the base, four-point target relation, scalar, five online phases, and
exact candidate identity:
`IC1N83Ckb1a1fb212480PDP4rootRCguidedLAgaussTDdirectISO0h69261fd42dac`.
The [local IC1 run rows](local_ic1_rows.json) bind both small-base controls
and K2816 to the exact exposed one-target workload `eebf98ee1031`; the
K1280 run ID is
`IC1N83Ckb1a1fb212480PDP4rootRCguidedLAgaussTDdirectISO0h69261fd42dacWeebf98ee1031R1`.

| Same-source, same-binary exposed-point control | K1248/root11 | K1280/root11 | K2543/root11 | K2816/root11 | K3072/root11 |
| --- | ---: | ---: | ---: | ---: | ---: |
| Actual usable base points | 207,168 | 212,480 | 422,138 | 467,456 | 509,952 |
| Folded columns | 1,248 | 1,280 | 2,543 | 2,816 | 3,072 |
| Exact root entries | 129,268,902 | 135,983,147 | 536,737,856 | 658,165,181 | 783,270,873 |
| Packed root allocation, bytes | 1,342,177,280 | 1,342,177,280 | 5,368,709,120 | 5,368,709,120 | 5,368,709,120 |
| Local receipt wall, s | 143.224 | 146.960 | 401.115 | 431.024 | 679.051 |
| One-target online, ms | 18.725 | **4.733** | 18.983 | 25.301 | 31.908 |
| Full arithmetic replay | Pass | Pass | Pass | Pass | Pass |

K1280 used 13,395 target probes and spent 4.609583 ms in point decomposition,
0.000042 ms on target query, 0.017750 ms on relation check, 0.000791 ms on
descent, and 0.105000 ms on scalar replay. Its setup stage was 78.850 s,
including exact-index construction, and rank stage was 67.969 s. K1248 used
10,622 target probes and spent 18.493792 ms in point decomposition; its
setup/rank stages were 64.641/78.385 s. Both completed under a 600 s cap.

K1280's observed online interval was **4.01 times shorter** than K2543's
on this one point; its receipt wall was **2.73 times shorter**. These are
shared-host exposed-point diagnostics, and a one-point interval has no
uncertainty estimate. The primary candidate comparison now needs K1248 and
K1280 alongside K2543, K3072, and the separately selected finalist on
previously unused points in alternating Tier 1 same-point IC–rho pairs. The
ten-point fixture for that comparison was sealed before this result.

The enumerated-set SHA-256 and producer base hash for each small base exactly
match the earlier 2026-10-16 cloud grid's K1248 and K1280 manifests. Their
target probe counts also match that grid (10,622 and 13,395). The new IC1
candidate IDs differ because their manifests bind the current source commit;
the earlier cloud and this shared-host run have different timing conditions.

All five controls used source SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
and local release binary SHA-256
`d71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906`.
K1248's base hash was
`b42195208a24d3d222c27c36a02f6b0b6075df076c7c9b3a2e099b7485bbdf54`;
K1280's was
`1d348fe3486147274bab0a15aae072d3ac2570a35883f356e5c802b29c58bf10`.
The public point, known-answer sidecar and frozen runner are pinned in
[the protocol](PROTOCOL.md).

The two small-base raw runs, including base dumps, target witnesses, exact
commands, flags and local replay scripts, are committed as
[raw.tar.gz](raw.tar.gz),
SHA-256
`f2b9cbc454705add4b8691785fefc8a4498764cc8066e0e6e294ce931f16f685`.
With the archive extracted to a scratch directory, the
[versioned verifier](verify_archive.py) checks its digest and independently
replays both bases, relations, scalars, online phases, and candidate IDs:

```sh
mkdir -p /private/tmp/n83-smallk-replay
tar -xzf research/n83_smallk_current_source_20261017/raw.tar.gz -C /private/tmp/n83-smallk-replay
python3 research/n83_smallk_current_source_20261017/verify_archive.py research/n83_smallk_current_source_20261017/raw.tar.gz /private/tmp/n83-smallk-replay
```

The [visual report](reports/run-20261017-current-source-controls/REPORT.md)
compares all current-source exposed-point controls. Canonical IC/rho graphs
are unchanged because this local control has no new tier-admitted fresh-point
IC/rho pair.
