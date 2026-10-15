#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)
= V1 bounded larger-order screen run
*2026-10-15 23:29 UTC · Darwin arm64, shared L0.* The native screen
processed all 30 eligible primes in the frozen 223–607 range. Eight
different-$j$ rows passed the algebraic predicates and 64 seeded map
witnesses each; a separate Rust implementation checked 512/512 map
witnesses, group orders and subgroup factors. Both processes exited zero
with empty stderr. `SCREEN_V1.jsonl`, `SCREEN_CHECK_V1.json` and
`FREEZE_V1.json` hold records, exact commands and hashes.

The bars display every accepted subgroup order. Blue means $r=p$;
orange means $r != p$. These are construction observations; this run
did not include new target solves or method timing. The post-screen
prime-to-characteristic selection is separately versioned in
`PANEL_PROTOCOL_V2.md`.

#image("GRAPH.svg", width: 100%)
