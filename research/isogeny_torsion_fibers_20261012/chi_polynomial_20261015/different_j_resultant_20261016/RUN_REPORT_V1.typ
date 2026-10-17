#set page(paper: "us-letter", margin: 0.7in)
#set text(font: "Arial", size: 9.5pt)
#text(size: 18pt, weight: "bold")[Run receipt: direct χ chord resultant exactness]

*Run ID:* `different_j_resultant_20261016/exhaustive-v1` \
*Host:* Darwin arm64, shared L0 \
*Input SHA-256:* `f6c3126b4d95e3c457c9e552471bc2cdff322de0e0a9de9fbba701d7c05897ba`

#image("RESULTANT_GRAPH_V1.svg", width: 100%)

All 231 requested point cases completed and a separate Rust checker
matched every pair set. The run counted 3,640 resultants, 7,410 gcds,
4,422 roots, 6,598 candidate additions and 6,358 rejected candidates.
The metered repeat used 7,094,272 bytes peak RSS and took 0.39 seconds
real time on a shared host.

The first meter command failed because sandboxed macOS denied
`sysctl kern.clockrate`; its records are preserved. The escalated
repeat passed replay. These observations do not give a cross-method
ratio. Sources: `RUN_REPORT_V1.md`, `RUN_SPEC_V1.json`,
`EXHAUSTIVE_V1.jsonl`, `EXHAUSTIVE_V1.check.json`.
