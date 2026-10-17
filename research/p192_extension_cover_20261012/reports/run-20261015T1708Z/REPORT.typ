#set page(paper: "a4", margin: 17mm)
#set text(font: "Libertinus Serif", size: 10pt)
= P-192 Type I signed-base three-sum census

Sign closure of the fixed 16-class Type I quotient base gives 32 signed classes and leaves 16 logarithm columns. Exact enumeration gives 5,504 distinct three-sum targets, of which 1,344 satisfy the existing two-character query law. The positive base gave 816 and 224; accepted support grows sixfold. Native Rust and independent PARI arithmetic agree. The same 64 frozen known-scalar queries give zero matches and rank 0/16.

#image("RESULT.svg", width: 100%)

== Construction

For each existing fixed-origin class $[D_i-T]$, add its negative $[T-D_i]$. Both have explicit degree-three divisor pairs. In the elliptic quotient the set is $plus.minus B_i$ for 16 indices. Each signed triple uses coefficients on those same 16 unknown logarithms. Lean checks the additive-homomorphism equation in `P192ExtensionCover.signed_three_log_equation` and reports only `propext` in its axiom audit.

== Exact census and query law

There are 5,984 unordered signed triples. Of these, 512 have a visible cancellation and map to 32 singleton targets. The remaining 5,472 targets are distinct and disjoint from the singleton targets. The character cells have 1,344, 1,416, 1,392, and 1,352 targets in `++`, `+-`, `-+`, and `--` order. The accepted-scalar lower bound from the existing P-192 density certificate gives the fixed-support bound $1344/L < 2^(-179)$; the positive-base result was $224/L < 2^(-182)$.

The 64 frozen accepted scalar queries satisfy $q=a(-2G)$ and have zero signed-base matches. Extra signed entries increase triple-index work. This run measures exact support and correctness, not online speed. The negative divisor pairs still need the two-cubic reduction to one fixed-origin positive degree-three representative each.

Source-linked evidence, exact inputs and reproduction commands are in `research/p192_extension_cover_20261012/reports/run-20261015T1708Z/REPORT.md`.
