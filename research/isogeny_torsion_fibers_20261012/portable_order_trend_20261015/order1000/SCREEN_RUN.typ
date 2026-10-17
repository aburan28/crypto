#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)

= Extended different-j order screen

*Run and independent replay, 2026-10-15.* The committed protocol fixed
prime range 619–1601, $p = 1 mod 3$, selection $r >= 1009$ and $r != p$,
the first accepted row, a 600-second wall cap and 512 MiB resident cap.
Host: Darwin arm64 shared L0. The exact Rust source/binary hashes,
commands, raw rows, exit codes and checker receipt are in
`PROTOCOL_V1.md` and `SCREEN_RUN.md`.

All 68 candidate rows were emitted. Eight passed the construction
conditions. Independent Rust formulas replayed 512/512 seeded map-law
witnesses and reproduced the first selected cell:
$p=1237, r=1201, h=1275$, with $j_("source")=1066$ and
$j_("destination")=727$. Both screen and checker exited 0.

#image("SCREEN_V1.svg", width: 100%)

The screen establishes a bounded construction and fixed selection.
Exhaustive rational-fiber, subgroup map-law and complete seven-arm
same-target recovery remain separate checks. Group-count pricing,
complete portable S, isolated wall timing and independent-host
performance have their own gates.
