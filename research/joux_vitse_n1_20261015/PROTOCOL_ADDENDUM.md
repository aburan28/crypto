# Protocol addendum (2026-10-15, written before its run reported)

**Why.**  After the frozen runs started, the scoreboard showed that the
repository already measures the plain and 2-torsion `k = 5` Joux–Vitse test
(`src/cryptanalysis/jv_quintic.rs`, `jv_quintic_edwards.rs`;
`research/notes/index-calculus/RESEARCH_K5_TORSION_JOUX_VITSE.md`, 2026-09-29/30):
`C″ = 3.0·10¹⁰` per plain test (F4 to degree 24, certificate at degree 19,
`8,826 × 8,444`).  This study's instrument (`joux_vitse.rs`) duplicates the
plain construction and measures `5.6·10⁹` per replayed trial and
`1.07·10¹⁰` per full F4 run (exploratory probe, p = 31).  The survey that
should have found `jv_quintic` before building was incomplete; that is
recorded in the README, and the duplicate is kept because it is the
instrument of record of the frozen runs.

**Added run (attribution control, not a claim).**  The committed
`jv_quintic` example, unchanged, on the current `f4_fp` engine, at its own
first rung:

```
cargo run --release --example jv_quintic -- --variant weierstrass --sizes 271 --seeds 1 \
  --residuals 3 --constructed 2 --max-degree 24 --budget 900 \
  --json research/joux_vitse_n1_20261015/data/jv_quintic_rerun/weierstrass_p271_s1.json
```

**Reading, fixed in advance.**  If its `C″` on today's engine is close to
`3.0·10¹⁰`, the gap to this study's full F4 (`1.07·10¹⁰`) is a difference
between the two drivers (degree bound 36 against 24, system build) and is
investigated before any attribution.  If it is close to `1.07·10¹⁰`, the
factor `≈ 2.8` is the F4 engine's own progress since 2026-09-29 (Gebauer–
Möller criteria, packed monomials, commits `3d4634596`, `cba784a41`) and
only the replay's remaining factor (`≈ 1.9`) is credited to the trace.
Either way the `k = 5` crossover is re-derived on the measured constant and
classed engineering.
