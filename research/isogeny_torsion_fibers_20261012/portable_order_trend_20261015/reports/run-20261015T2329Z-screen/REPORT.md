# V1 bounded larger-order screen run

Start: 2026-10-15 23:29 UTC. Host: Darwin arm64, shared L0. The
[frozen protocol and hash manifest](../../PROTOCOL_V1.md) fixed the
ascending `223 <= p <= 607`, `p = 1 mod 3` screen and 180-second/512 MiB
caps before execution. Exact commands:

```sh
gtimeout 180 target/order-batch-check/release/examples/different_j_trend_screen screen-trend > SCREEN_V1.jsonl 2> SCREEN_V1.stderr
target/order-batch-check/release/examples/different_j_trend_check SCREEN_V1.jsonl > SCREEN_CHECK_V1.json 2> SCREEN_CHECK_V1.stderr
```

Requested/completed: 30/30 eligible primes, eight accepted rows, 64/64
seeded map witnesses per accepted row, 512/512 independent checker
replays. Both exit codes were zero and both stderr files empty. The complete
[screen rows](../../SCREEN_V1.jsonl), [checker receipt](../../SCREEN_CHECK_V1.json),
and [graph source](../../SCREEN_V1.svg) are retained. The bars show the
eight accepted subgroup orders. Blue bars have `r=p`; orange bars have
`r != p`. These are construction results only; no new target solve or
method timing was measured in this run.

![Accepted subgroup orders](GRAPH.svg)

The [V2 panel protocol](../../PANEL_PROTOCOL_V2.md) preserves this V1 output
and freezes comparable prime-to-characteristic cells for the next run.
