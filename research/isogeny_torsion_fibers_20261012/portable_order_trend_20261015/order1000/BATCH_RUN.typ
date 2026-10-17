#set page(paper: "a4", margin: 17mm)
#set text(size: 10pt)

= Order-1,201 seven-arm cold batch: run record

The committed `RUN_SPEC_V2.json` froze the first selected cell
$(p,r,h)=(1237,1201,1275)$, 64 public targets, seven arms, three
rotated rounds, prefixes 1/4/16/64, caps and exact source/binary/input
hashes before child measurement. The initial target-freeze attempt
failed with an unlisted `j_pair()` prime before writing targets; its
V1 stderr and exit are retained. Corrected V2 target freeze exited 0.

The exact command was `zsh .../order1000/run_order1201_v2.zsh`.
All 84/84 cold child exits were 0. The independent Rust checker
recomputed 1,785/1,785 scalar outputs and 12,486 relation witnesses,
including map identities, rows and recovery. Run SHA-256:
`d314764e87988b2dc12675d199612f50abc538d1782761c562c2c74be5e34d89`.
The replay certificate, per-child stdout/stderr/meter/exit files,
commands, host and source links are in `BATCH_RUN.md`.

#image("TREND_FIVE_V2.svg", width: 100%)

The figure shows the complete *group-counted* cold cost per target,
not complete portable S. Fp counts overlap group calls, and native
hash/table/row/allocation/I/O work lacks a pinned common price.
Shared Darwin arm64 timing is L0 descriptive. An L2-or-better
physical-core wall session and different-host replay are distinct
gates. The V1 graph with an overlapping label is preserved; V2 only
repairs the label position.
