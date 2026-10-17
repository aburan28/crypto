#set page(paper: "us-letter", margin: 0.7in)
#set text(font: "Arial", size: 10pt)
#text(size: 19pt, weight: "bold")[Four-bin direct inverse-Vélu order screen]

*2026-10-16 · bounded structural input selection*

The frozen native screen examined 513 eligible base primes in the range
1238 through 10007. Its independent Rust checker recalculated all rows
and replayed 5,440 sampled map identities. It found first eligible cells
at subgroup orders 2,437, 5,569, and 9,679. The fourth bin had no
eligible cell in this bounded prime range.

#image("SCREEN_GRAPH_V1.svg", width: 100%)

The three selected cells are structural inputs. Exact subgroup targets,
complete relation recovery, portable price, and isolated performance
remain separate gates. The fourth bin needs a separately frozen range
extension. No canonical ECDLP frontier row changes from this screen.

Evidence: `FOUR_GATE_PROTOCOL_V1.md`, `SCREEN_FREEZE_V1.json`,
`SCREEN_V1.jsonl`, `SCREEN_V1.check.json`, and `SCREEN_REPORT_V1.md`.
The raw JSONL SHA-256 is
`ef232a54b22938239819fd6b14182d1e49de5778e9c00542ce257c54097608a5`.
