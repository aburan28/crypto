#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)
= V2 complete three-cell target panel run
*2026-10-15 23:41 UTC · Darwin arm64, shared L0.* The frozen 64-target
seven-arm panel completed at $r=421,433,457$ with three rotated rounds
and four fresh cold prefixes per arm. All 252 children exited zero.
The independent Rust checker replayed all 5,355 scalar outputs and
37,061 relation witnesses. The exact source/binary/input digests, caps,
commands, raw child receipts and replay hashes are linked from
`REPORT.md` and `BATCH_RUN_SPEC_V1.json`.

The graph includes the independently checked order-229 reference.
It plots median $S_("group")$ per target at N=64; exact counted
additions and doublings include common work, setup, queries and final
verification. Three rounds provide min–max values in `TREND_V1.json`.
Field and nonarithmetic work remains separately counted or unpriced;
macOS L0 timing is descriptive.

#image("GRAPH.svg", width: 100%)
