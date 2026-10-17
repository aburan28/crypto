# K992 remains the local N83 cold-frontier candidate under a matched repeat

K992/root11 completed the full cold process in 75.841 seconds on the exposed
control point. A matched K992/root22 repeat completed in 76.330 seconds with
the same enumerated base; the two layouts are effectively tied in this
one-point local screen. K960/root11 was close in cold time but used 680,614
target probes, versus 117,653 at K992. K928/root11 spent 74.417 seconds in
the guided rank stage, making its complete cold process 98.488 seconds.

| Configuration | Usable base points | Root index, s | Rank stage, s | Cold process, s | Target online, ms | Target probes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| K928/root11 | 154,048 | 23.696 | 74.417 | 98.488 | 48.323 | 168,560 |
| K960/root11 | 159,360 | 25.764 | 51.121 | 77.389 | 176.344 | 680,614 |
| K992/root11 | 164,672 | 27.728 | 47.751 | **75.841** | 30.587 | 117,653 |
| K992/root22 | 164,672 | 27.038 | 48.886 | 76.330 | 32.964 | 117,653 |

The [run report](RUN.pdf), [results report](RESULT.pdf), [graph](screen.svg)
and [graph PDF](screen.pdf) show the index/rank balance and
point-dependent target lookup. This is an L0 shared Apple M4 Pro diagnostic
on one previously exposed public point. Its cold-time ranking is not a Tier 1
or fresh-point selection. In an earlier local screen under a different host
load, K992/root11 took 95.223 seconds; the 75.841-second matched repeat shows
why a single local wall interval does not determine the winning configuration.
The sealed ten-point Tier 1 K992/K1024/K1107 comparison addresses that gate.

Every run used `examples/koblitz_orbit_dlp_fast.rs` at SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`
and the executed binary at SHA-256
`d71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906`.
The exposed point file SHA-256 was
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`.
All runs used the same 75%-load dynamic root index, symmetric S3 state,
sparse canonical map, four guided rank threads and one target thread; only K
or the root11/root22 slot layout changed. Each fresh process rebuilt its
base, index, and factor-base logs before solving the point.

The [independent arithmetic replay](replay.json) rebuilt every base point,
summed each four-point target relation to the public point, checked the
recovered scalar, and checked the exclusive online phase sum; all four rows
passed. The [raw archive](raw.tar.gz) has SHA-256
`02993c7b3b7b7bf85af28bf999acfc48b41df2d61afa9408dfde3d9f440d5f55`.
Its 28 members include complete base dumps, process output, stderr, result
receipts, replay, and graph artifacts. `verify_archive.py` checked every
member against the [manifest](manifest.json) and its local counterpart.
