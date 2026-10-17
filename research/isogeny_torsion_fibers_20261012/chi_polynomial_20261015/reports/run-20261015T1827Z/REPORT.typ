#set page(paper: "a4", margin: 19mm)
#set text(size: 10pt)
= Direct Semaev root and quotient exactness gate
*Run:* 2026-10-15 18:27 UTC. *Inputs:* frozen V2 specification and public `p=61,r=19` fixture.

The native source arm solves the third summation polynomial as a quadratic in the second signed-base abscissa. The quotient arm first applies the explicit degree-three map to base points and target. Both check each root by substitution and each candidate decomposition by exact source addition. Infinity, vertical, tangent, linear and all-root cases are handled explicitly.

Both arms returned exactly the independently enumerated signed-pair sets on 22/22 targets: 19 subgroup points and three exceptional points. There were 144 decompositions and 2,992 independent pair-addition and quotient-homomorphism checks. The source arm made 320 polynomial calls and 288 square-root attempts; the quotient arm made 288 and 256. Neither hit the 10⁷-call cap. Polynomial-root candidates rejected by source addition reflect the missing sign in an abscissa root, not returned decompositions.

#image("GRAPH.svg", width: 100%)

The exhaustive claim is bounded to the stated field, 16-point signed base and 22 targets. Raw SHA-256: `571a6774214064b531a1dad7a50ae969b6bb3c2c4103c87d893acf4bbdbb58a5`. The paired relation solve and separate scalar replay have their own records.
