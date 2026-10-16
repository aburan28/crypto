#set page(paper: "us-letter", margin: 0.7in)
#set text(font: "Arial", size: 9.5pt)
#text(size: 18pt, weight: "bold")[Direct χ chord resultant: exhaustive pair control]

*2026-10-16 · frozen different-j order-229 cell*

The eliminated chord polynomial was intersected with the signed
factor-base χ polynomial. All roots were lifted and checked by full
source addition. A separate Rust checker rebuilt the curve and
enumerated every pair with independent field arithmetic.

#image("RESULTANT_GRAPH_V1.svg", width: 100%)

Across all 229 subgroup points and the two torsion exceptions, the
returned unordered pair set matched direct enumeration in 231 of 231
cases. There were 136 distinct valid pairs and 132 cases with no
decomposition. The oracle made 3640 resultants, 7410 gcds, 94 split
attempts and 6598 exact candidate additions; 6358 candidates were
rejected. The metered repeat used 7,094,272 bytes peak RSS, below the
512 MiB ceiling, and its pair rows matched the first run exactly.

The complete relation solve, portable primitive price and isolated
wall comparison have not been measured by this exactness run. Sources:
`RUN_SPEC_V1.json`, `EXHAUSTIVE_V1.jsonl`, `CHECKER_SPEC_V1.json`,
`EXHAUSTIVE_V1.check.json`, `RESULTS_V1.md`.
