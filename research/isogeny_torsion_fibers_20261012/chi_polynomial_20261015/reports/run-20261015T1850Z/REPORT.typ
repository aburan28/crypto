#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= Public-only χ root input and exactness gate
*Run:* 2026-10-15 18:50 UTC. *Protocol:* frozen V3. *Host:* Darwin arm64.

Native code removed exactly four `answer_scalar` fields from the private checker fixture, producing one redacted public JSONL for the solver. The executable rejected the private file with exit 101 before timing and relation collection. The separate checker verified that public and private objects differ only by those four fields. The public-only exactness run used the redacted file.

Direct source S3 and mapped quotient S3 roots each matched independently enumerated signed-pair sets on 22/22 smallest-cell targets: 19 subgroup points and three exceptional points. There were 144 decompositions and 2,992 independent pair-addition and quotient-homomorphism checks. Source and quotient polynomial calls were 320 and 288, respectively; neither arm hit its cap.

#image("GRAPH.svg", width: 100%)

This is a custody and bounded exactness result for the named fixtures and the 16-point signed base at p=61. It does not formalize the universal summation-polynomial criterion. Public fixture SHA-256: `bfe512594006202477845574a1778b6b41a7b9acfddae1ca45a0478756fb2a4e`; exactness raw SHA-256: `e899f9ebc015f8f051b4b0a0919c032b8542601cfc06977629b19a98e53b5110`.
