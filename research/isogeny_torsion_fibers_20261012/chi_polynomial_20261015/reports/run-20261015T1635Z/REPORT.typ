#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)

= Exhaustive χ-polynomial decomposition check

*Run:* EXHAUSTIVE_V1, 2026-10-15. *Host:* shared Darwin, L0. *Status:* 22/22 exact-set cases matched independent pair enumeration.

The source is $E: y^2=x^3+1$ over $F_(61^2)$, subgroup order $r=19$, with 16 signed projected base points. The four-cell fixture SHA-256 is `1a63cc22a22e0fe27d44359325d502038d26d3ac616374fc5558c041939218b7`. Native source commit: `8cde58c74`; source SHA-256: `8591945b8085e300caad38f5bfabaa7db5ba676826558761089b4559a0cf7e9e`.

#image("GRAPH.svg", width: 100%)

The enumerative solver checked 352 base-polynomial evaluations and 2,992 pair equations. It sent 370 candidates to exact group addition: 144 pair appearances were accepted and 226 were filtered. Ten canonical χ-key collisions required retained lift information. The 144 appearances include the order-19 identity check twice, once as `subgroup/0` and again as `special/0`; the 19 subgroup cases account for 136 unordered base pairs. `T` and `2T` returned empty sets in both implementations.

The raw records are `../../EXHAUSTIVE_V1.jsonl` (SHA-256 `b10720711ccb95e62797a464017a55cf345f0c43cef4b44ae9476f4b9e873638`) and `../../EXHAUSTIVE_V1.stderr`. The exact command and identity are in `REPORT.md` beside this PDF. This is a correctness gate for an enumerative polynomial solver. Larger-cell relation solves, complete costs, optimized algebraic root finding, and isolated timing remain pending.
