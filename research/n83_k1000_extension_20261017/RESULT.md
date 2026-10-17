# K984, K1000 and K1016 add no faster exposed-point cold control

Three adjacent root11 factor bases reached full rank and solved the same
previously exposed N83 point, but all three took longer than the preceding
K992/root11 local repeat at 75.841 seconds. K1016 was the shortest of this
extension at 94.452 seconds; K1000 took 115.501 seconds, and K984 took
148.974 seconds. The exact target-online intervals were 285.381, 114.163,
and 364.063 milliseconds respectively. All three passed independent base,
four-point relation, scalar, and exclusive-phase replay.

| Configuration | Usable base points | Root index, s | Rank stage, s | Cold process, s | Target online, ms | Target probes |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| K984/root11 | 163,344 | 46.775 | 101.329 | 148.974 | 364.063 | 1,139,099 |
| K1000/root11 | 166,000 | 36.006 | 78.951 | 115.501 | 114.163 | 318,308 |
| K1016/root11 | 168,656 | 36.007 | 57.710 | 94.452 | 285.381 | 992,582 |

The [run report](RUN.pdf), [results report](RESULT.pdf), [graph](screen.svg)
and [graph PDF](screen.pdf) show these measurements
against a dashed K992 local reference from the preceding study. That reference
was recorded in a different local interval; the line is context, not a paired
speedup estimate. These are L0 shared Apple M4 Pro diagnostics on one exposed
point, with no host-isolation or fresh-point selection gate. The sealed ten-point
Tier 1 K992/K1024/K1107 panel remains the selection experiment.

Every process used `examples/koblitz_orbit_dlp_fast.rs` at SHA-256
`8a2fe55dd90dc33892c197d27260182bf783fc8212fb485d8297daf2bbba8e12`,
executed binary SHA-256
`d71a9662c9d4e2dde58443bfc646a05737c6b9ca6b3d5e831b632ce09384c906`,
and exposed point SHA-256
`2c068d2f659588edcdbaa685fc73bb146e2bfee3be123d377d537658f8021c6d`.
All used the dynamic 75%-load root11 index, symmetric S3 state, sparse
canonical map, four rank threads and one target thread. Each fresh process
built the base, index and logs before its one-target solve.

The [independent arithmetic replay](replay.json) records three
`PASS_INDEPENDENT_ARITHMETIC_REPLAY` rows. The [raw archive](raw.tar.gz), SHA-256
`1ebe3cc221b03ce2003078806f4dc2d8fc6f8c26a7976b4368386f62daced5c3`,
contains 23 members with complete base dumps, stdout, stderr, target results,
receipts, replay and graph assets. `verify_archive.py` checked every member
against the [manifest](manifest.json).
