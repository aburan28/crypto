#set page(paper: "a4", margin: 20mm)
#set text(size: 10pt)

= Different-j degree-three candidate screen

*Run:* SCREEN_V1, 2026-10-15. *Status:* first selector pass at $p=67$. *Host:* shared Darwin, L0.

The frozen selector started with $E: y^2=x^3-21x+37$ and rational $T=(3,1)$, searched ascending primes in $[61,271]$ with $p=1 mod 3$, and required a multiplicity-one subgroup prime $r≥19$. At $p=61$, $n_2=3675=3 dot 5^2 dot 7^2$ had no eligible factor. At $p=67$, both curves had $n_2=4575=3 dot 5^2 dot 61$, $r=61$, and cofactor 75. Source coefficients were $(46,37)$ with $j=7$; destination coefficients were $(53,25)$ with $j=11$.

#image("GRAPH.svg", width: 100%)

The implemented Vélu sum map had kernel $\{O,(3,1),(3,66)\}$ and passed 64 sampled on-curve and homomorphism checks. Different j-values exclude an isomorphism for the nonsingular pair; the map's full-point check and DLP transfer are separate pending gates. The raw record is `../../SCREEN_V1.jsonl`, SHA-256 `be2094594f2ca0e7ec49d221a8d9865b5d5ecf640922cd0541520f414ac29707`. The command, source hash, and source link are in `REPORT.md` beside this PDF. See the MIT 18.783 lecture slides on j-invariants for the isomorphism invariant.
