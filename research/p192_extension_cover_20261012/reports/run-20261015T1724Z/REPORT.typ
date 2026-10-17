#set page(paper: "a4", margin: 18mm)
#set text(font: "Libertinus Serif", size: 10pt)
= P-192 Type I negative classes at one fixed origin

All 16 negatives of the previously reduced Type I base classes now have positive degree-three divisor representatives relative to the same origin $T$. Swapping the old input pair changes its class from $[B-A]$ to $[A-B]$. Native PARI constructs two cubic sections per class, and a separate checker replays all 32 cuts, residual coordinates, squarefreeness, disjointness, and rank-nine interpolation.

#image("REDUCTION.svg", width: 100%)

== Class identity

For the swapped pair, cubic $G$ cuts $A+T+U$ with $deg U=6$. Cubic $F$ cuts $B+U+D_-$ with $deg D_-=3$. Their common line-bundle class implies $[D_--T]=[A-B]=-[B-A]$. Lean checks this implication in `negative_fixed_origin_cubic_class_identity` under the two section-divisor equalities. The concrete P-192 quartic and sections are verified by native finite-field arithmetic.

== Exact replay

The frozen input contains the quartic, three anchor points, and 16 original degree-three pairs. Each negative reduction succeeds with the first denominator trial $[1,1,1,1]$. The producer writes all cubic coefficients and triangular $x/y$ residual polynomials. The separate verifier recomputes both resultants and the unique numerator cubic, ending with `NEGATIVE_FIXED_ORIGIN_VERIFY_PASS classes=16 cuts=32 prime_bits=192`.

This construction makes the 32 signed classes available at one fixed origin. It does not change the earlier fixed-base three-sum support or the 64-query zero-rank result. The next collector must measure ordinary-query yield, useful rank, and signed-index cost on a larger frozen workload.

Exact input and source hashes, commands, proof boundary, and links: `research/p192_extension_cover_20261012/reports/run-20261015T1724Z/REPORT.md`.
