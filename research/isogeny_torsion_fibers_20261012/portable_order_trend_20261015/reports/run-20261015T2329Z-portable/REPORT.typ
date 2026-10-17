#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)
= V1 portable-count derivation run
*2026-10-15 23:29 UTC · Darwin arm64, shared L0.* One deterministic native
count derivation read the immutable order-229 run and independent replay.
It checked all 84 arm-round-prefix records and all 1,785 verified scalar
outputs. Exit code was zero; stderr was empty. Source, binary and input
digests and exact command are in `FREEZE_V1.json` and `REPORT.md`.

$S_("group") = ("adds" + "doubles") / (N sqrt(229))$ is a counted-group
lower bound. Map field work, table/hash operations, row arithmetic and
process overhead do not yet have a pinned conversion. The graph shows
N=64 medians of three complete cold children per arm; bar lengths are log
transformed and exact values appear beside them. The source result report
holds every count axis and matched ratio.

#image("GRAPH.svg", width: 100%)
