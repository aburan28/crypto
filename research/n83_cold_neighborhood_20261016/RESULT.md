# N83 cold neighborhood: K1152 leads the seven-arm screen

On the frozen exposed public point, K1152/root22 completed the full cold
process in **148.919 s**, compared with 150.531 s for K1024/root22 in the
same cloud allocation. K1152/dynamic root11 followed at 149.586 s. All seven
IC arms and the same-point strong-rho answer passed independent base,
four-point relation and scalar replay. The K1024/dynamic A/A bookends were
150.611 s and 150.397 s, a 0.142% difference against the frozen 10% gate;
each arm had zero new cgroup throttle, OOM and OOM-kill events. The verifier
reported `SCREEN_COMPLETE` and advanced the two K1152 layouts to the four
pre-frozen reserve points. The selected factor base remains null until that
confirmation.

| Arm | Usable base points | Cold process, s | One-target online, ms | Peak RSS, GiB | Replay |
| --- | ---: | ---: | ---: | ---: | --- |
| K1024/dynamic root11 A | 169,984 | 150.611 | 69.302 | 1.187 | pass |
| K896/dynamic root11 | 148,736 | 169.562 | 69.460 | 1.069 | pass |
| K1152/dynamic root11 | 191,232 | 149.586 | 314.148 | 1.944 | pass |
| K1271/dynamic root11 | 210,986 | 153.311 | 223.458 | 2.077 | pass |
| K1024/root22 | 169,984 | 150.531 | 68.910 | 1.630 | pass |
| K1152/root22 | 191,232 | **148.919** | 318.498 | 2.505 | pass |
| K1024/dynamic root11 B | 169,984 | 150.397 | 69.007 | 1.186 | pass |

The repeated K1024/dynamic mean was 150.504 s, within 0.027 s of
K1024/root22. At K1024 the dynamic 11-byte state layout used about 0.443
GiB less peak RSS than root22 in this allocation. At K1152 the dynamic
layout used about 0.561 GiB less peak RSS and took 0.667 s longer in cold
process time. These small timing differences are descriptive until the
reserved-point comparison. The same-point rho online interval was 639.053
ms; this one-point screen ran rho after the IC arms and is not an alternating
Tier 1 online comparison. The cloud runner exposed four AMD EPYC 9V74 CPUs
and 15.61 GiB host memory.

The [canonical IC1 row record](ic1_rows_38883248405.json) has six immutable
candidate manifests, seven runs and workload `eebf98ee1031`. The
K1152/root22 candidate is
`IC1N83Ckb1a1fb191232PDP4rootRCguidedLAgaussTDdirectISO0ha2b0336762f1`;
K1152/dynamic root11 is
`IC1N83Ckb1a1fb191232PDP4rootRCguidedLAgaussTDdirectISO0h0ce483b7d6c3`.
The K1024/dynamic and K1271/dynamic IDs match the earlier exact-index
screen because their factor bases and executed producer source are unchanged;
the two K1024 repetitions have separate run IDs. The candidate manifests
identify the earlier exact producer revision, while each new run row records
this screen's execution commit and binary hash.

The [analysis](analysis_38883248405.json), [source receipt](source_38883248405.json)
and [host receipt](host_38883248405.json) are checked in. The complete raw
artifact is [GitHub Actions run 38883248405](https://github.com/aburan28/crypto/actions/runs/38883248405),
artifact ID `12039579271`. Its exact ZIP is archived at
`/Volumes/SSD990-2/n83-archive-20261016/raw/n83-cold-neighborhood-38883248405.zip`,
SHA-256 `323d03355f16c688e7d16656f6c8be27a00e9f655f22f34aed381c0a4759f985`.
The analysis SHA-256 is
`176271318ee5d3dd78a0bf6d212deb2a663c7ff9670d244c878862d3c04fc4bb`;
independent local replay reproduced it exactly. The executed IC binary
SHA-256 is
`5287c6d1de3b71a8e2044c07db26ae31ac79cd2a42af3607469c2119a7c16d28`.
`materialize.py` checks the raw archive, analysis, source, workload and each
complete base dump before writing candidate and run records; rerunning it is
idempotent. The full serial release library suite on this branch passed
4,148 tests, with zero failures and 117 ignored. Its log is preserved at
`/Volumes/SSD990-2/n83-archive-20261016/logs/n83-cold-neighborhood-libtest.log`,
SHA-256 `45b4661903f9dfbe7dd246c3314792c29169770724d971903b64a7c4ca69cee0`.

Next, compare the two K1152 layouts with the earlier K1024/root22 and
K1536/root22 finalists on all four frozen reserve points. Preserve all
per-point cold intervals, A/A controls, same-point rho receipts and failures.
After a factor base is selected, freeze five new distinct points for the
alternating Tier 1 online comparison.
