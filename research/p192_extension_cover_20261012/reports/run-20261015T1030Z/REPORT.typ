#set page(paper: "a4", margin: 18mm)
#set text(font: "Libertinus Serif", size: 10pt)
= P-192 Type I fixed-origin Jacobian reduction

All 16 exact quotient-coordinate base classes and 64 accepted guided known-scalar queries have one degree-three divisor representative over the P-192 prime field relative to a common three-point origin. Native PARI checks 160 cubic-section cuts. The guided collector accepts 64 of 259 scalar trials and finds zero three-sum relations and rank 0/16 for the fixed base.

#image("RESULT.svg", width: 100%)

== Construction and proof

For a prior degree-three pair $[B-A]$ and fixed origin $T$, a cubic section cuts $B+T+U$ with $deg U=6$. The unique cubic through $A+U$ cuts $A+U+D$ with $deg D=3$. Both are sections of the same line bundle on the smooth plane quartic, so $[D-T]=[B-A]$. Lean checks this cancellation as `P192ExtensionCover.fixed_origin_cubic_class_identity` under the section-divisor hypotheses; its axiom audit reports only `propext`.

== Evidence and interpretation

The native replay reads `type_i_fixed_origin_native_input.gp` and checks resultants, squarefree residuals, quotient-ring $y$ coordinates, rank-nine interpolation, and all 80 degree-three representatives. Its log ends `NATIVE_FIXED_ORIGIN_PASS base=16 guided=64 cuts=160 prime_bits=192`. Earlier Sage candidate generation and its independent replay are archived with their exact runtime receipts; native replay is the acceptance check.

The character law admits 224 of the fixed base's 816 distinct three-sum targets. Of 259 seeded scalar attempts, 128 fail the first fiber condition, 67 fail the second, and 64 are accepted. All 64 accepted queries retain $q=a(-2G,-2G,-2G)$ and have fixed-origin representatives. None gives a base three-sum relation; the incidence rank is 0/16. The certified ideal-uniform conditioned hit bound is $224/L < 2^(-182)$ for this base and query policy.

Source-linked details, exact field and subgroup parameters, and reproduction commands: `research/p192_extension_cover_20261012/reports/run-20261015T1030Z/REPORT.md`.
