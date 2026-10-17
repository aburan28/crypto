#set page(paper: "a4", margin: 17mm)
#set text(font: "Libertinus Serif", size: 10pt)
= P-192 Type I fixed-origin mixed-base collector

All 160 signed quotient classes in the frozen Type I base now have fixed-origin degree-three genus-three Jacobian representatives. Exact enumeration gives 333,184 distinct nonzero-row three-sum targets, 83,874 of which pass the two-character geometric query law. A steered known-scalar law produces 64 accepted three-sum hits, but all are known-only and give rank 0/16. A disjoint fresh accepted-scalar law also gives rank 0/16 on 64 queries.

#image("RESULT.svg", width: 100%)

== Geometry and equations

The base has 32 signed unknown classes and 128 signed guided classes whose scalars are recorded. Every class is represented as $[D-T]$ with rational effective degree-three divisors on the P-192 Type I plane quartic and one fixed origin $T$. For each of the 64 new negative guided classes, native PARI constructs two cubic sections and independently checks their resultants, residual points, squarefreeness, disjointness, and rank-nine interpolation. With the 80 earlier positive and 16 earlier negative classes, the evidence covers 320 cubic cuts.

Lean checks the conditional class identity and the mixed known-log equation. A triple with one unknown and two known entries gives an unknown-log row after subtracting the known contributions. A triple of three known entries gives no unknown row. The concrete P-192 divisor and quotient-map hypotheses are checked by the recorded arithmetic.

== Exact census

There are 695,520 unordered triples: 357,760 known-only, 2,048 with zero unknown row, and 335,712 with a nonzero row. The latter map to 333,184 distinct targets. Their character cells are `++`: 83,874; `+-`: 82,884; `-+`: 83,246; `--`: 83,180. The known-only triples map to 349,696 targets; there is zero overlap between the two target sets in this frozen base. The accepted rank-bearing support is exactly 62.40625 times the earlier 32-point signed base, with the same 16 unknown columns and a larger triple index.

The steered law selects 64 accepted known-only targets and retains each scalar equation. All 64 are hits by construction and produce zero unknown-log rows. The independent SHA-256-counter law obtains 64 fresh accepted scalars in 246 attempts and has zero rank-bearing matches. The certified lower bound $L=1569275433846670190958947353987675545133315652035359680560$ for the ideal uniform accepted-scalar law gives fixed-base hit probability at most $83874/L < 2^(-173)$; Lean checks this integer inequality. The deterministic sample is recorded separately from that distributional bound.

Source links, exact commands, binary and input hashes, per-phase costs, and the proof map are in `research/p192_extension_cover_20261012/reports/run-20261015T1800Z/REPORT.md`.
